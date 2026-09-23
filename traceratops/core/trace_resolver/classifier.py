"""Modular barcode-multiplicity classification for trace resolution.

The likelihood classifier models every designed barcode, including barcodes
which were not observed in a trace.  Its parameters are consequently fitted
once to a complete trace table rather than separately to individual traces.
"""

from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, gammaln, logsumexp
from scipy.stats import chi2

from .results import TraceMultiplicity


class TraceClassification(str, Enum):
    UNCHANGED = "unchanged"
    CLEAN_ONE = "clean_one"
    RESOLVE_TWO = "resolve_two"


@dataclass(frozen=True)
class ClassificationThresholds:
    """Provisional candidate thresholds intended for simulation calibration."""

    min_repeated_barcodes: int = 3
    min_repeated_fraction: float = 0.3
    candidate_rule: str = "both"

    def __post_init__(self):
        if self.min_repeated_barcodes < 1:
            raise ValueError("min_repeated_barcodes must be positive")
        if not 0 <= self.min_repeated_fraction <= 1:
            raise ValueError("min_repeated_fraction must be between 0 and 1")
        if self.candidate_rule not in {"both", "either"}:
            raise ValueError("candidate_rule must be 'both' or 'either'")


def classify_trace(multiplicity, thresholds=None):
    """Classify a trace without embedding a fixed biological interpretation."""
    if not isinstance(multiplicity, TraceMultiplicity):
        multiplicity = TraceMultiplicity.from_barcodes(multiplicity)
    thresholds = thresholds or ClassificationThresholds()
    if multiplicity.n_repeated_barcodes == 0:
        return TraceClassification.UNCHANGED
    tests = (
        multiplicity.n_repeated_barcodes >= thresholds.min_repeated_barcodes,
        multiplicity.fraction_repeated_barcodes >= thresholds.min_repeated_fraction,
    )
    is_two = all(tests) if thresholds.candidate_rule == "both" else any(tests)
    return TraceClassification.RESOLVE_TWO if is_two else TraceClassification.CLEAN_ONE


def build_multiplicity_matrix(
    table, trace_column="Trace_ID", barcode_column="Barcode #"
):
    """Return ``(trace_ids, barcode_ids, counts)`` for the full barcode universe.

    The universe is inferred from all barcode identities in *table*.  Thus an
    absent barcode gets an explicit zero in each trace's row, rather than being
    silently omitted from its likelihood.
    """
    trace_ids = np.unique(np.asarray(table[trace_column]))
    barcode_ids = np.unique(np.asarray(table[barcode_column]))
    counts = np.zeros((len(trace_ids), len(barcode_ids)), dtype=int)
    trace_index = {value: i for i, value in enumerate(trace_ids)}
    barcode_index = {value: i for i, value in enumerate(barcode_ids)}
    for trace, barcode in zip(table[trace_column], table[barcode_column]):
        counts[trace_index[trace], barcode_index[barcode]] += 1
    return trace_ids, barcode_ids, counts


def barcode_heterogeneity_pvalue(counts):
    """Test whether dataset-wide excess multiplicity is uniform by barcode.

    A G test is used on ``sum(max(K_i - 1, 0))``.  With no excess events there
    is no evidence of heterogeneity.  This aggregate test deliberately avoids
    interpreting stable barcode effects as a biological off-target mechanism.
    """
    counts = np.asarray(counts, dtype=float)
    burden = np.maximum(counts - 1, 0).sum(axis=0)
    total = burden.sum()
    if burden.size < 2 or total == 0:
        return 1.0
    expected = total / burden.size
    positive = burden > 0
    statistic = 2 * np.sum(burden[positive] * np.log(burden[positive] / expected))
    return float(chi2.sf(statistic, burden.size - 1))


def _barcode_log_likelihood(counts, n_polymers, detection_efficiency, lambdas):
    """Per-trace log likelihood under the Binomial-plus-Poisson model."""
    counts = np.asarray(counts, dtype=int)
    lambdas = np.broadcast_to(np.asarray(lambdas, dtype=float), counts.shape)
    terms = []
    p = detection_efficiency
    for detected_copies in range(n_polymers + 1):
        valid = counts >= detected_copies
        off_target = counts - detected_copies
        log_binomial = (
            gammaln(n_polymers + 1)
            - gammaln(detected_copies + 1)
            - gammaln(n_polymers - detected_copies + 1)
            + detected_copies * np.log(p)
            + (n_polymers - detected_copies) * np.log1p(-p)
        )
        log_poisson = (
            off_target * np.log(lambdas)
            - lambdas
            - gammaln(off_target + 1)
        )
        terms.append(np.where(valid, log_binomial + log_poisson, -np.inf))
    return logsumexp(np.stack(terms), axis=0).sum(axis=1)


@dataclass(frozen=True)
class LikelihoodTraceScores:
    """Likelihood diagnostics and decision for one trace."""

    log_likelihood_single: float
    log_likelihood_doublet: float
    log_likelihood_ratio: float
    posterior_doublet_probability: float
    classification: TraceClassification


class LikelihoodMultiplicityClassifier:
    """Dataset-fitted maximum-likelihood singlet/doublet classifier.

    The model currently assumes one genuine detection probability ``p`` shared
    by all barcodes.  In experimental data, stable barcode-specific detection
    efficiency differences can therefore also contribute to the excess-count
    heterogeneity statistic; it should not be interpreted as measuring only
    off-target biology.

    Barcode-specific nuisance rates use an L2 penalty
    ``strength * sum((lambda_i - lambda_global)**2)`` where ``lambda_global``
    is the global-model estimate.  This permits strong hotspots while keeping
    weakly informed rates near the stable shared estimate.
    """

    def __init__(
        self,
        off_target_model="auto",
        heterogeneity_alpha=0.01,
        lambda_regularization=10.0,
        posterior_threshold=0.5,
    ):
        if off_target_model not in {"global", "barcode-specific", "auto"}:
            raise ValueError(
                "off_target_model must be global, barcode-specific, or auto"
            )
        if not 0 < heterogeneity_alpha < 1:
            raise ValueError("heterogeneity_alpha must be between zero and one")
        if lambda_regularization < 0:
            raise ValueError("lambda_regularization must be non-negative")
        if not 0 <= posterior_threshold <= 1:
            raise ValueError("posterior_threshold must be between zero and one")
        self.requested_off_target_model = off_target_model
        self.heterogeneity_alpha = heterogeneity_alpha
        self.lambda_regularization = lambda_regularization
        self.posterior_threshold = posterior_threshold

    @staticmethod
    def _decode(parameters, n_rates):
        p = expit(parameters[0])
        doublet_prior = expit(parameters[1])
        lambdas = np.exp(parameters[2 : 2 + n_rates])
        return p, doublet_prior, lambdas

    def _fit_model(self, counts, n_rates, global_rate=None):
        n_barcodes = counts.shape[1]
        rough_rate = max(float(np.maximum(counts - 2, 0).mean()), 0.01)
        upper_rate = max(5.0, float(counts.max() + 1))
        bounds = [(-9, 9), (-9, 9)] + [
            (np.log(1e-8), np.log(upper_rate))
        ] * n_rates

        def objective(parameters):
            p, prior, rates = self._decode(parameters, n_rates)
            lambdas = np.repeat(rates, n_barcodes) if n_rates == 1 else rates
            single = _barcode_log_likelihood(counts, 1, p, lambdas)
            doublet = _barcode_log_likelihood(counts, 2, p, lambdas)
            mixture = logsumexp(
                np.vstack((np.log1p(-prior) + single, np.log(prior) + doublet)),
                axis=0,
            )
            penalty = 0.0
            if global_rate is not None:
                penalty = self.lambda_regularization * np.sum(
                    (rates - global_rate) ** 2
                )
            return -mixture.sum() + penalty

        starts = ((0.2, 0.2), (0.4, 0.8), (0.7, 0.3), (0.85, 0.8), (0.55, 0.5))
        rate_start = global_rate if global_rate is not None else rough_rate
        best = None
        for p_start, prior_start in starts:
            initial = np.r_[
                np.log(p_start / (1 - p_start)),
                np.log(prior_start / (1 - prior_start)),
                np.repeat(np.log(max(rate_start, 1e-8)), n_rates),
            ]
            fitted = minimize(objective, initial, method="L-BFGS-B", bounds=bounds)
            if best is None or fitted.fun < best.fun:
                best = fitted
        if best is None or not np.isfinite(best.fun):
            raise RuntimeError("Multiplicity likelihood optimization failed")
        return self._decode(best.x, n_rates)

    def fit(self, counts, barcode_ids=None):
        """Fit dataset-level parameters from a trace-by-barcode count matrix."""
        counts = np.asarray(counts, dtype=int)
        if counts.ndim != 2 or not counts.size or np.any(counts < 0):
            raise ValueError("counts must be a non-empty, non-negative 2D matrix")
        self.barcode_ids = np.asarray(
            np.arange(counts.shape[1]) if barcode_ids is None else barcode_ids
        )
        if len(self.barcode_ids) != counts.shape[1]:
            raise ValueError("barcode_ids must match the count matrix columns")
        p, prior, global_rates = self._fit_model(counts, 1)
        self.global_off_target_rate = float(global_rates[0])
        self.heterogeneity_pvalue = barcode_heterogeneity_pvalue(counts)
        use_specific = self.requested_off_target_model == "barcode-specific" or (
            self.requested_off_target_model == "auto"
            and self.heterogeneity_pvalue < self.heterogeneity_alpha
        )
        if use_specific:
            p, prior, rates = self._fit_model(
                counts, counts.shape[1], self.global_off_target_rate
            )
            self.off_target_model = "barcode-specific"
        else:
            rates = global_rates
            self.off_target_model = "global"
        self.detection_efficiency = float(p)
        self.doublet_prior = float(prior)
        self.lambdas = np.repeat(rates, counts.shape[1]) if len(rates) == 1 else rates
        self.used_barcode_specific_rates = use_specific
        return self

    def score(self, counts, has_repeated_barcodes=None):
        """Score and classify one or more count vectors."""
        counts = np.asarray(counts, dtype=int)
        one = counts.ndim == 1
        matrix = counts[None, :] if one else counts
        single = _barcode_log_likelihood(
            matrix, 1, self.detection_efficiency, self.lambdas
        )
        doublet = _barcode_log_likelihood(
            matrix, 2, self.detection_efficiency, self.lambdas
        )
        log_odds = (
            doublet - single
            + np.log(self.doublet_prior)
            - np.log1p(-self.doublet_prior)
        )
        posterior = expit(log_odds)
        repeated = np.any(matrix > 1, axis=1)
        if has_repeated_barcodes is not None:
            repeated = np.broadcast_to(has_repeated_barcodes, len(matrix))
        output = []
        for ll1, ll2, probability, has_repeat in zip(
            single, doublet, posterior, repeated
        ):
            if probability >= self.posterior_threshold:
                classification = TraceClassification.RESOLVE_TWO
            elif has_repeat:
                classification = TraceClassification.CLEAN_ONE
            else:
                classification = TraceClassification.UNCHANGED
            output.append(
                LikelihoodTraceScores(
                    float(ll1),
                    float(ll2),
                    float(ll2 - ll1),
                    float(probability),
                    classification,
                )
            )
        return output[0] if one else output
