import numpy as np
from astropy.table import Table

from traceratops import trace_splitter
from traceratops.core.trace_resolver import (
    LikelihoodMultiplicityClassifier,
    TraceClassification,
    barcode_heterogeneity_pvalue,
    build_multiplicity_matrix,
)
from traceratops.core.trace_resolver.classifier import _barcode_log_likelihood


def _simulate(seed, n_traces=500, n_barcodes=20, p=0.65, pi=0.35, rates=0.03):
    rng = np.random.default_rng(seed)
    polymers = 1 + (rng.random(n_traces) < pi)
    true = rng.binomial(polymers[:, None], p, (n_traces, n_barcodes))
    return true + rng.poisson(np.broadcast_to(rates, (n_traces, n_barcodes)))


def test_full_multiplicity_matrix_contains_unobserved_barcodes():
    table = Table(
        rows=[("a", 1), ("a", 1), ("b", 2)], names=("Trace_ID", "Barcode #")
    )
    trace_ids, barcodes, counts = build_multiplicity_matrix(table)

    assert list(trace_ids) == ["a", "b"]
    assert list(barcodes) == [1, 2]
    assert counts.tolist() == [[2, 0], [0, 1]]


def test_zero_observation_barcodes_change_the_likelihood():
    observed_only = _barcode_log_likelihood([[1]], 1, 0.6, [0.1])[0]
    with_unobserved = _barcode_log_likelihood([[1, 0]], 1, 0.6, [0.1, 0.1])[0]

    assert with_unobserved != observed_only


def test_high_efficiency_singlets_and_doublets_are_distinguished():
    training = _simulate(1, p=0.9, pi=0.5, rates=0.01)
    classifier = LikelihoodMultiplicityClassifier("global").fit(training)

    assert classifier.score(np.ones(20, dtype=int)).classification in {
        TraceClassification.UNCHANGED,
        TraceClassification.CLEAN_ONE,
    }
    assert (
        classifier.score(np.full(20, 2)).classification
        == TraceClassification.RESOLVE_TWO
    )


def test_low_efficiency_and_no_repeat_traces_are_still_scored():
    classifier = LikelihoodMultiplicityClassifier("global").fit(
        _simulate(2, p=0.25, pi=0.75, rates=0.01)
    )
    score = classifier.score(np.r_[np.ones(7, dtype=int), np.zeros(13, dtype=int)])

    assert np.isfinite(score.log_likelihood_single)
    assert np.isfinite(score.log_likelihood_doublet)
    assert 0 <= score.posterior_doublet_probability <= 1


def test_positive_lambda_makes_singlet_repetitions_possible():
    classifier = LikelihoodMultiplicityClassifier("global").fit(
        _simulate(3, rates=0.12)
    )
    score = classifier.score(np.r_[2, np.ones(12, dtype=int), np.zeros(7, dtype=int)])

    assert np.isfinite(score.log_likelihood_single)


def test_zero_off_target_rate_has_finite_fit_scores_and_posteriors():
    counts = _simulate(10, n_traces=800, p=0.7, pi=0.4, rates=0.0)
    classifier = LikelihoodMultiplicityClassifier("global").fit(counts)
    scores = classifier.score(counts)

    assert np.isfinite(classifier.global_off_target_rate)
    assert classifier.global_off_target_rate < 1e-4
    assert all(np.isfinite(score.log_likelihood_single) for score in scores)
    assert all(np.isfinite(score.log_likelihood_doublet) for score in scores)
    assert all(
        np.isfinite(score.posterior_doublet_probability) for score in scores
    )


def test_global_fit_recovers_nuisance_and_mixture_parameters():
    classifier = LikelihoodMultiplicityClassifier("global").fit(
        _simulate(4, n_traces=1200, p=0.6, pi=0.7, rates=0.08)
    )

    assert abs(classifier.global_off_target_rate - 0.08) < 0.04
    assert abs(classifier.doublet_prior - 0.7) < 0.12


def test_fit_supports_either_mixture_component_dominating():
    singlet_dominated = LikelihoodMultiplicityClassifier("global").fit(
        _simulate(5, pi=0.15)
    )
    doublet_dominated = LikelihoodMultiplicityClassifier("global").fit(
        _simulate(6, pi=0.85)
    )

    assert singlet_dominated.doublet_prior < 0.5
    assert doublet_dominated.doublet_prior > 0.5


def test_auto_model_detects_only_strong_heterogeneity():
    homogeneous = _simulate(7, n_traces=800, rates=0.04)
    hotspot_rates = np.full(20, 0.02)
    hotspot_rates[3] = 0.8
    hotspot = _simulate(8, n_traces=800, rates=hotspot_rates)

    homogeneous_fit = LikelihoodMultiplicityClassifier("auto").fit(homogeneous)
    hotspot_fit = LikelihoodMultiplicityClassifier("auto").fit(hotspot)

    assert homogeneous_fit.off_target_model == "global"
    assert hotspot_fit.off_target_model == "barcode-specific"
    assert hotspot_fit.heterogeneity_pvalue < 0.01
    assert hotspot_fit.lambdas[3] > np.median(hotspot_fit.lambdas)
    assert barcode_heterogeneity_pvalue(np.minimum(homogeneous, 1)) == 1.0


def test_specific_rates_are_regularized_toward_global_rate():
    rates = np.full(20, 0.03)
    rates[0] = 0.5
    counts = _simulate(9, n_traces=250, rates=rates)
    unregularized = LikelihoodMultiplicityClassifier(
        "barcode-specific", lambda_regularization=0
    ).fit(counts)
    regularized = LikelihoodMultiplicityClassifier(
        "barcode-specific", lambda_regularization=100
    ).fit(counts)

    weak_barcode = 1
    regularized_distance = abs(
        regularized.lambdas[weak_barcode] - regularized.global_off_target_rate
    )
    unregularized_distance = abs(
        unregularized.lambdas[weak_barcode] - unregularized.global_off_target_rate
    )
    assert regularized_distance <= unregularized_distance
    assert regularized.lambdas[0] > regularized.lambdas[weak_barcode]


def test_cli_exposes_likelihood_selection_and_keeps_threshold_default():
    parser = trace_splitter.parse_arguments()
    assert parser.parse_args([]).multiplicity_classifier == "threshold"
    args = parser.parse_args(
        [
            "--multiplicity-classifier",
            "likelihood",
            "--likelihood-off-target-model",
            "barcode-specific",
        ]
    )
    assert args.multiplicity_classifier == "likelihood"
    assert args.likelihood_off_target_model == "barcode-specific"
