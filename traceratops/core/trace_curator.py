"""Empirical, conservative curation of already resolved chromatin traces.

No barcode ordering or numeric identity is interpreted as genomic position.
See ``docs/source/api/trace_curator.md`` for the fitting and input contract.
"""

import warnings
from array import array
from dataclasses import dataclass, field
from itertools import combinations, islice
from math import ceil, floor, log10, nan
from types import MappingProxyType
from typing import Any, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.table import Table

__all__ = [
    "TraceCuratorModel",
    "TraceCurationIteration",
    "TraceCurationResult",
    "fit_curator_model",
    "curate_trace",
]


@dataclass(frozen=True)
class TraceCuratorModel:
    """Fitted references and calibration metadata; created by fit_curator_model.

    Distributions are sorted, immutable float arrays, keyed by exact separation
    or an unordered frozenset of barcode identities. Positions share one genomic
    coordinate system. Calibration counts exclude unscoreable traces.
    """

    reference_distributions: Mapping
    genomic_positions: Mapping
    global_threshold: float
    reference_mode: str
    minimum_reference_observations: int
    trace_fpr: float
    n_reference_traces: int
    n_reference_pairs: int
    n_calibration_traces: int
    n_calibration_unscoreable: int
    calibration_source: str
    chromosome: Optional[object] = None
    n_crossfit_folds: Optional[int] = None
    n_reference_keys: int = field(init=False)
    n_supported_reference_keys: int = field(init=False)
    median_reference_observations: float = field(init=False)
    min_reference_observations_actual: int = field(init=False)
    fraction_reference_keys_below_minimum: float = field(init=False)

    def __post_init__(self):
        # bytes backing prevents callers from re-enabling array writes.
        references = {}
        for key, values in self.reference_distributions.items():
            values = np.asarray(values, dtype=float)
            if values.ndim != 1 or not np.isfinite(values).all():
                raise ValueError("Reference arrays must be finite and one-dimensional")
            if (values < 0).any() or (np.diff(values) < 0).any():
                raise ValueError("Reference arrays must be nonnegative and sorted")
            references[key] = np.frombuffer(values.tobytes(), dtype=float)
        object.__setattr__(
            self, "reference_distributions", MappingProxyType(references)
        )
        object.__setattr__(
            self,
            "genomic_positions",
            MappingProxyType(_validated_positions(self.genomic_positions)),
        )
        _positive_integer(
            self.minimum_reference_observations, "minimum_reference_observations"
        )
        if self.reference_mode not in {"separation", "barcode_pair"}:
            raise ValueError("reference_mode must be separation or barcode_pair")
        if not np.isfinite(self.global_threshold) or self.global_threshold < 0:
            raise ValueError("global_threshold must be finite and nonnegative")
        if not np.isfinite(self.trace_fpr) or not 0 < self.trace_fpr < 1:
            raise ValueError("trace_fpr must be in (0, 1)")

        sizes = np.asarray([len(sample) for sample in references.values()], dtype=int)
        supported = int(np.count_nonzero(sizes >= self.minimum_reference_observations))
        for name, value in {
            "n_reference_keys": len(sizes),
            "n_supported_reference_keys": supported,
            "median_reference_observations": (
                float(np.median(sizes)) if len(sizes) else nan
            ),
            "min_reference_observations_actual": int(sizes.min()) if len(sizes) else 0,
            "fraction_reference_keys_below_minimum": (
                (len(sizes) - supported) / len(sizes) if len(sizes) else nan
            ),
        }.items():
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class TraceCurationIteration:
    """Compact audit of one state, including the terminal stopping state."""

    iteration: int
    n_localizations_before: int
    C_top3_before: float
    global_threshold: float
    globally_abnormal: bool
    n_scoreable_pairs: int
    candidate_spot_id: Optional[object] = None
    candidate_barcode: Optional[object] = None
    loo_top3: float = nan
    n_tied_top3: int = 0
    top5_tiebreak_used: bool = False
    loo_top5: float = nan
    n_tied_top5: int = 0
    removal_accepted: bool = False
    stop_reason: str = ""


@dataclass(frozen=True)
class TraceCurationResult:
    """Curation outcome with independent pandas copies preserving input columns.

    DataFrames are caller-owned mutable output; diagnostics are immutable.
    Removed rows are in removal order. NaN costs denote unscoreable states.
    """

    curated_trace: pd.DataFrame
    removed_localizations: pd.DataFrame
    iterations: Tuple[TraceCurationIteration, ...]
    initial_cost: float
    final_cost: float
    global_threshold: float
    n_removed: int
    terminal_state: str
    stop_reason: str


def _positive_integer(value, name, allow_zero=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        raise ValueError(
            f"{name} must be {'nonnegative' if allow_zero else 'positive'}"
        )


def _validated_positions(mapping):
    try:
        positions = {barcode: float(position) for barcode, position in mapping.items()}
    except (TypeError, ValueError) as error:
        raise ValueError(
            "Genomic mapping values must be finite numeric positions"
        ) from error
    values = tuple(positions.values())
    if not np.isfinite(values).all() or len(set(values)) != len(values):
        raise ValueError("Genomic mapping positions must be finite and unique")
    return positions


def _frame(trace):
    if isinstance(trace, Table):
        return trace.to_pandas()
    if isinstance(trace, pd.DataFrame):
        return trace.copy()
    raise TypeError("Traces must be pandas DataFrames or astropy Tables")


def _barcode_column(frame):
    if "Barcode" in frame and "Barcode #" in frame:
        if not frame["Barcode"].equals(frame["Barcode #"]):
            raise ValueError("Barcode and Barcode # columns disagree")
    return "Barcode #" if "Barcode #" in frame else "Barcode"


def _validate_trace(trace, genomic_positions=None, chromosome=None):
    frame = _frame(trace)
    barcode_column = _barcode_column(frame)
    required = {"x", "y", "z", "Spot_ID", barcode_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required trace columns: {sorted(missing)}")
    for column in (barcode_column, "Spot_ID"):
        if frame[column].isna().any():
            raise ValueError(f"{column} must not contain missing identities")
        if frame[column].duplicated().any():
            message = (
                "Duplicated Barcode: resolve duplicated localizations before curation"
                if column == barcode_column
                else "Spot_ID must be unique within each trace"
            )
            raise ValueError(message)
    if "Trace_ID" in frame:
        if frame["Trace_ID"].isna().any() or frame["Trace_ID"].nunique() > 1:
            raise ValueError("Expected one trace with a single nonmissing Trace_ID")
    if "Chrom" in frame and len(frame):
        if frame["Chrom"].isna().any() or frame["Chrom"].nunique() != 1:
            raise ValueError("Trace genomic positions must lie on one chromosome")
        if chromosome is not None and frame["Chrom"].iloc[0] != chromosome:
            raise ValueError("Trace chromosome disagrees with the fitted model")
    try:
        xyz = frame[["x", "y", "z"]].to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("x, y, z coordinates must be finite numeric values") from error
    if not np.isfinite(xyz).all():
        raise ValueError("x, y, z coordinates must be finite numeric values")
    barcodes = frame[barcode_column].to_numpy()
    try:
        if "Genomic_Position" in frame:
            positions = frame["Genomic_Position"].to_numpy(dtype=float)
        elif genomic_positions is not None:
            positions = np.asarray(
                [genomic_positions[b] for b in barcodes], dtype=float
            )
        else:
            raise ValueError(
                "Genomic_Position or explicit genomic_positions mapping is required"
            )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "Finite genomic coordinates are required for every barcode"
        ) from error
    if not np.isfinite(positions).all() or len(np.unique(positions)) != len(positions):
        raise ValueError("Genomic positions must be finite and unique within a trace")
    if genomic_positions is not None:
        for barcode, position in zip(barcodes, positions):
            if barcode in genomic_positions and genomic_positions[barcode] != position:
                raise ValueError(
                    "Genomic_Position disagrees with supplied/fitted positions"
                )
    return frame, barcodes, positions, xyz


def _traces(source):
    """Reiterate table groups without retaining copied traces during fitting."""
    if isinstance(source, (pd.DataFrame, Table)):
        frame = source.to_pandas() if isinstance(source, Table) else source
        if "Trace_ID" not in frame:
            raise ValueError(
                "Trace_ID is required for a combined reference/calibration table"
            )
        if frame["Trace_ID"].isna().any():
            raise ValueError("Trace_ID must not be missing")
        yield from (
            trace for _, trace in frame.groupby("Trace_ID", sort=False, observed=True)
        )
    else:
        yield from source


def _pair_key(mode, barcodes, positions, left, right):
    if mode == "separation":
        separation = float(abs(positions[right] - positions[left]))
        if not np.isfinite(separation):
            raise ValueError("Genomic separation exceeds finite numeric range")
        return separation
    return frozenset((barcodes[left], barcodes[right]))


def _empirical_anomaly(sample, observed):
    lower = (np.searchsorted(sample, observed, side="right") + 1) / (len(sample) + 1)
    upper = (len(sample) - np.searchsorted(sample, observed, side="left") + 1) / (
        len(sample) + 1
    )
    return -log10(min(1.0, 2.0 * min(lower, upper)))


def _pair_anomalies(barcodes, positions, xyz, references, mode, minimum):
    edges = []
    for left, right in combinations(range(len(barcodes)), 2):
        sample = references.get(_pair_key(mode, barcodes, positions, left, right))
        if sample is None or len(sample) < minimum:
            continue
        with np.errstate(over="ignore", invalid="ignore"):
            distance = float(np.linalg.norm(xyz[right] - xyz[left]))
        if not np.isfinite(distance):
            raise ValueError("Spatial distance exceeds finite numeric range")
        score = _empirical_anomaly(sample, distance)
        if np.isfinite(score):
            edges.append((left, right, score))
    return tuple(sorted(edges, key=lambda edge: edge[2], reverse=True))


def _top_cost(edges, k, without=None):
    values = list(
        islice(
            (
                score
                for left, right, score in edges
                if without not in (left, right) and np.isfinite(score)
            ),
            k,
        )
    )
    return float(sum(values)) if values else nan


def _loo_scores(edges, n_localizations, k, candidates=None):
    before = _top_cost(edges, k)
    nodes = range(n_localizations) if candidates is None else candidates
    return {node: before - _top_cost(edges, k, without=node) for node in nodes}


def _maximal_candidates(scores):
    finite = {node: score for node, score in scores.items() if np.isfinite(score)}
    if not finite:
        return ()
    best = max(finite.values())
    return tuple(node for node, score in finite.items() if score == best)


def _select_candidate(edges, n_localizations):
    primary = _loo_scores(edges, n_localizations, 3)
    tied = _maximal_candidates(primary)
    secondary = {}
    remaining = ()
    if len(tied) > 1:
        secondary = _loo_scores(edges, n_localizations, 5, candidates=tied)
        remaining = _maximal_candidates(secondary)
    winner = (
        tied[0] if len(tied) == 1 else (remaining[0] if len(remaining) == 1 else None)
    )
    return winner, primary, tied, secondary, remaining


def _empirical_threshold(values, fpr):
    ordered = np.sort(np.asarray(values, dtype=float))
    if not len(ordered):
        raise ValueError(
            "No scoreable clean calibration traces; insufficient reference support"
        )
    allowed = floor(fpr * len(ordered))
    return float(ordered[max(0, len(ordered) - allowed - 1)])


def _fit_reference(traces, genomic_positions, reference_mode):
    """Build one reference using compact double buffers, never float lists."""
    positions_by_barcode = dict(genomic_positions)
    observations: dict[object, array] = {}
    n_reference_traces = 0
    n_reference_pairs = 0
    chromosome = None
    for trace in traces:
        frame, barcodes, positions, xyz = _validate_trace(trace, positions_by_barcode)
        if "Chrom" in frame and len(frame):
            chrom = frame["Chrom"].iloc[0]
            if chromosome is not None and chrom != chromosome:
                raise ValueError(
                    "References must share one chromosome coordinate system"
                )
            chromosome = chrom
        for barcode, position in zip(barcodes, positions):
            positions_by_barcode[barcode] = float(position)
        n_reference_traces += 1
        for left, right in combinations(range(len(frame)), 2):
            key = _pair_key(reference_mode, barcodes, positions, left, right)
            with np.errstate(over="ignore", invalid="ignore"):
                distance = float(np.linalg.norm(xyz[right] - xyz[left]))
            if not np.isfinite(distance):
                raise ValueError("Spatial distance exceeds finite numeric range")
            if key not in observations:
                observations[key] = array("d")
            observations[key].append(distance)
            n_reference_pairs += 1
    if len(set(positions_by_barcode.values())) != len(positions_by_barcode):
        raise ValueError(
            "Each barcode must have a unique genomic position in the model"
        )
    references = {}
    # Pop converted buffers immediately so compact buffers and NumPy copies do
    # not coexist for the entire reference population.
    while observations:
        key, buffer = observations.popitem()
        sample = np.frombuffer(buffer, dtype=np.float64).copy()
        del buffer
        sample.sort()
        references[key] = sample
    return (
        references,
        positions_by_barcode,
        chromosome,
        n_reference_traces,
        n_reference_pairs,
    )


def _calibration_costs(traces, references, positions, chromosome, mode, minimum):
    """Yield only scalar costs; each pairwise scoring state dies before yield."""
    for trace in traces:
        frame, barcodes, genomic, xyz = _validate_trace(trace, positions, chromosome)
        cost = nan
        if len(frame) >= 3:
            cost = _top_cost(
                _pair_anomalies(barcodes, genomic, xyz, references, mode, minimum), 3
            )
        yield cost


def _fold_traces(population, fold, n_folds, *, held_out):
    """Round-robin by deterministic trace order; no hash or random seed."""
    for index, trace in enumerate(_traces(population)):
        if (index % n_folds == fold) == held_out:
            yield trace


def _record_costs(costs, destination):
    unscoreable = 0
    for cost in costs:
        if np.isfinite(cost):
            destination.append(cost)
        else:
            unscoreable += 1
    return unscoreable


def _calibration_threshold(costs, trace_fpr):
    minimum_tail_sample = ceil(1 / trace_fpr)
    if len(costs) < minimum_tail_sample:
        raise ValueError(
            f"Too few scoreable calibration traces ({len(costs)}) to estimate the "
            f"requested empirical tail (trace_fpr={trace_fpr:g}); at least "
            f"{minimum_tail_sample} finite calibration scores are required. "
            "Supply more traces or check reference support."
        )
    if len(costs) * trace_fpr < 5:
        warnings.warn(
            f"Only {len(costs)} scoreable calibration traces are available. "
            f"The requested {100 * trace_fpr:g}% tail is estimable but may be "
            "unstable (fewer than five expected tail observations).",
            UserWarning,
            stacklevel=3,
        )
    return _empirical_threshold(costs, trace_fpr)


def fit_curator_model(
    traces,
    *,
    genomic_positions=None,
    trace_fpr=0.01,
    minimum_reference_observations=20,
    reference_mode="separation",
    calibration_traces=None,
    n_crossfit_folds=5,
) -> TraceCuratorModel:
    """Fit references and a cross-fitted empirical trace-abnormality cutoff.

    The fitting population must contain a majority of valid resolved polymers;
    individual aberrant traces need not be identified. Cross-fitting prevents
    self-inclusion, but contamination robustness has not been quantified.
    By default trace i belongs to fold i % n_crossfit_folds in collection order
    (first-appearance Trace_ID order for combined tables). Each fold is scored
    against references from the other folds, then the final reference is fitted
    once on all traces. Combined DataFrames/Tables and reiterable collections
    work; one-pass iterators require separate calibration_traces.

    Explicit calibration_traces bypass cross-fitting and are intended for an
    independent high-quality control population. In either mode the requested
    tail requires ceil(1 / trace_fpr) finite calibration scores; fewer than five
    expected tail observations emit a statistical warning. Only compact pair
    buffers, the current reference, and scalar calibration costs are retained.
    """
    if not np.isfinite(trace_fpr) or not 0 < trace_fpr < 1:
        raise ValueError("trace_fpr must be in (0, 1)")
    _positive_integer(minimum_reference_observations, "minimum_reference_observations")
    if reference_mode not in {"separation", "barcode_pair"}:
        raise ValueError("reference_mode must be separation or barcode_pair")
    genomic_positions = (
        _validated_positions(genomic_positions) if genomic_positions is not None else {}
    )
    costs = array("d")
    unscoreable = 0
    if calibration_traces is None:
        if iter(traces) is traces:
            raise ValueError(
                "One-pass reference iterators require separate calibration_traces; default cross-fitting requires reiterable input"
            )
        # Convert a combined Astropy table only once, without materializing a
        # list of copied single traces or any calibration scoring histories.
        population = traces.to_pandas() if isinstance(traces, Table) else traces
        n_traces = sum(1 for _ in _traces(population))
        _positive_integer(n_crossfit_folds, "n_crossfit_folds")
        if not 2 <= n_crossfit_folds <= n_traces:
            raise ValueError(
                f"n_crossfit_folds must be between 2 and the number of available traces ({n_traces})"
            )
        for fold in range(n_crossfit_folds):
            references, positions, chromosome, _, _ = _fit_reference(
                _fold_traces(population, fold, n_crossfit_folds, held_out=False),
                genomic_positions,
                reference_mode,
            )
            unscoreable += _record_costs(
                _calibration_costs(
                    _fold_traces(population, fold, n_crossfit_folds, held_out=True),
                    references,
                    positions,
                    chromosome,
                    reference_mode,
                    minimum_reference_observations,
                ),
                costs,
            )
            del references
        threshold = _calibration_threshold(costs, trace_fpr)
        references, positions, chromosome, n_reference_traces, n_reference_pairs = (
            _fit_reference(
                _traces(population),
                genomic_positions,
                reference_mode,
            )
        )
        calibration_source = "crossfit"
    else:
        references, positions, chromosome, n_reference_traces, n_reference_pairs = (
            _fit_reference(
                _traces(traces),
                genomic_positions,
                reference_mode,
            )
        )
        unscoreable = _record_costs(
            _calibration_costs(
                _traces(calibration_traces),
                references,
                positions,
                chromosome,
                reference_mode,
                minimum_reference_observations,
            ),
            costs,
        )
        threshold = _calibration_threshold(costs, trace_fpr)
        calibration_source = "separate_traces"
    return TraceCuratorModel(
        reference_distributions=references,
        genomic_positions=positions,
        global_threshold=threshold,
        reference_mode=reference_mode,
        minimum_reference_observations=minimum_reference_observations,
        trace_fpr=float(trace_fpr),
        n_reference_traces=n_reference_traces,
        n_reference_pairs=n_reference_pairs,
        n_calibration_traces=len(costs),
        n_calibration_unscoreable=unscoreable,
        calibration_source=calibration_source,
        chromosome=chromosome,
        n_crossfit_folds=int(n_crossfit_folds) if calibration_traces is None else None,
    )


def curate_trace(
    trace, model: TraceCuratorModel, *, max_removals=3
) -> TraceCurationResult:
    """Curate one resolved trace, preserving at least three localizations.

    Exact top3 ties invoke top5 only within the tied set. Remaining ties abstain.
    The trace is rescored after every accepted removal. Missing relationships
    are ignored, but deleting a node must leave at least one scoreable pair.
    The fitted genomic mapping supplies positions when the column is absent.
    """
    _positive_integer(max_removals, "max_removals", allow_zero=True)
    current, barcodes, positions, xyz = _validate_trace(
        trace, model.genomic_positions, model.chromosome
    )
    removed: list[pd.DataFrame] = []
    audit: list[TraceCurationIteration] = []
    initial_cost = nan
    while True:
        edges = (
            _pair_anomalies(
                barcodes,
                positions,
                xyz,
                model.reference_distributions,
                model.reference_mode,
                model.minimum_reference_observations,
            )
            if len(current) >= 3
            else ()
        )
        cost = _top_cost(edges, 3) if len(current) >= 3 else nan
        if not audit:
            initial_cost = cost
        abnormal = bool(np.isfinite(cost) and cost > model.global_threshold)
        reason = ""
        if len(current) < 3:
            reason = "too_few_localizations"
        elif not np.isfinite(cost):
            reason = "insufficient_reference"
        elif not abnormal:
            reason = "global_cost_normal"
        elif len(current) <= 3:
            reason = "minimum_trace_size_reached"
        elif len(removed) >= max_removals:
            reason = "max_removals_reached"
        diagnostics: dict[str, Any] = {}
        candidate = None
        if not reason:
            candidate, primary, tied, secondary, remaining = _select_candidate(
                edges, len(current)
            )
            diagnostics.update(
                n_tied_top3=len(tied),
                top5_tiebreak_used=len(tied) > 1,
                n_tied_top5=len(remaining),
                loo_top3=primary[tied[0]] if tied else nan,
                loo_top5=secondary[remaining[0]] if remaining else nan,
            )
            if candidate is None:
                reason = "ambiguous_candidate" if tied else "no_scoreable_candidate"
            else:
                diagnostics.update(
                    candidate_spot_id=current["Spot_ID"].iloc[candidate],
                    candidate_barcode=barcodes[candidate],
                )
        audit.append(
            TraceCurationIteration(
                iteration=len(audit),
                n_localizations_before=len(current),
                C_top3_before=cost,
                global_threshold=model.global_threshold,
                globally_abnormal=abnormal,
                n_scoreable_pairs=len(edges),
                removal_accepted=not bool(reason),
                stop_reason=reason,
                **diagnostics,
            )
        )
        if reason:
            break
        removed.append(current.iloc[[candidate]].copy())
        keep = np.arange(len(current)) != candidate
        current = current.iloc[np.flatnonzero(keep)].copy()
        barcodes, positions, xyz = barcodes[keep], positions[keep], xyz[keep]
    terminal = (
        "unscoreable"
        if not np.isfinite(cost) and not removed
        else (
            "abnormal_unresolved"
            if abnormal or not np.isfinite(cost)
            else "curated" if removed else "normal"
        )
    )
    removed_frame = pd.concat(removed) if removed else current.iloc[:0].copy()
    return TraceCurationResult(
        current,
        removed_frame,
        tuple(audit),
        initial_cost,
        cost,
        model.global_threshold,
        len(removed),
        terminal,
        reason,
    )
