#!/usr/bin/env python3
"""Cycle-4 benchmark for safe iterative trace curation.

This is an exploratory analysis, not the production ``trace_curator`` CLI.
It separates global trace abnormality from candidate actionability and permits
an abnormal trace to stop unresolved rather than forcing another deletion.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import os
import resource
import sys
import time
from collections import Counter
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from types import SimpleNamespace
from typing import Mapping, Sequence

import matplotlib
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIXED_FPRS = (0.001, 0.01, 0.05)
POLICIES = ("A_global", "B_absolute", "B_relative", "C_cap")
OUTPUTS = (
    "cycle4_trace_metrics.ecsv",
    "cycle4_condition_metrics.ecsv",
    "cycle4_iteration_audit.ecsv",
    "cycle4_thresholds.ecsv",
    "cycle4_runtime_profile.ecsv",
)
_PROGRESS_EPOCH = time.monotonic()
_PROGRESS_CONTEXT = ""
_GROUP_PEAK_MB = 0.0


def _rss_mb() -> float:
    """Live resident MiB on Linux; process high-water MiB elsewhere."""
    global _GROUP_PEAK_MB
    try:
        pages = int(Path("/proc/self/statm").read_text().split()[1])
        rss = pages * os.sysconf("SC_PAGE_SIZE") / 1024**2
    except (OSError, IndexError, ValueError):
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss = peak / (1024**2 if sys.platform == "darwin" else 1024)
    _GROUP_PEAK_MB = max(_GROUP_PEAK_MB, rss)
    return rss


def _progress(message: str) -> None:
    """Emit a timestamped progress message immediately to stderr."""
    elapsed = time.monotonic() - _PROGRESS_EPOCH
    print(
        f"[{elapsed:10.1f}s] RSS={_rss_mb():.1f} MiB {_PROGRESS_CONTEXT} {message}",
        file=sys.stderr,
        flush=True,
    )


def _pair_key(left, right) -> tuple[str, str]:
    """Return an unordered pair without imposing a numeric barcode type."""
    return tuple(sorted((str(left), str(right))))


def _cycle3():
    path = Path(__file__).with_name("analyze_trace_curator_loo_validation.py")
    spec = importlib.util.spec_from_file_location("_trace_curator_cycle3", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write(frame: pd.DataFrame, path: Path) -> None:
    safe = frame.copy()
    for column in safe.select_dtypes(include="object"):
        safe[column] = safe[column].fillna("").astype(str)
    Table.from_pandas(safe).write(path, format="ascii.ecsv", overwrite=True)


def fold_seed_sets(
    seeds: Sequence[int], evaluation_seed: int, calibration_seed: int | None = None
) -> tuple[set[int], set[int], set[int]]:
    evaluation = {int(evaluation_seed)}
    calibration = set() if calibration_seed is None else {int(calibration_seed)}
    reference = set(map(int, seeds)) - evaluation - calibration
    if reference & evaluation or reference & calibration or evaluation & calibration:
        raise AssertionError("Reference/calibration/evaluation seed leakage")
    return reference, calibration, evaluation


def empirical_threshold(values: Sequence[float], fpr: float) -> float:
    """Return the largest clean cutoff allowing at most ``fpr`` exceedances."""
    finite = np.sort(np.asarray(values, float)[np.isfinite(values)])
    if not len(finite):
        return math.nan
    allowed = math.floor(float(fpr) * len(finite))
    return float(finite[max(0, len(finite) - allowed - 1)])


def empirical_anomaly(reference: Sequence[float], observed: float) -> float:
    """Score against a reference array sorted when its model was fitted."""
    sample = np.asarray(reference, float)
    lower = (np.searchsorted(sample, observed, side="right") + 1) / (len(sample) + 1)
    upper = (len(sample) - np.searchsorted(sample, observed, side="left") + 1) / (
        len(sample) + 1
    )
    return -math.log10(min(1.0, 2.0 * min(lower, upper)))


@dataclass(frozen=True)
class TraceScores:
    """One shared pairwise scoring structure used by every policy."""

    trace: pd.DataFrame
    edges: tuple[tuple[int, int, float], ...]
    cost: float
    loo: np.ndarray
    cost_without: np.ndarray


@dataclass(frozen=True)
class CompactTraceStats:
    """Only the three scalars required for exact clean calibration."""

    cost: float
    best_loo: float
    best_relative_loo: float

    @classmethod
    def from_scores(cls, state: TraceScores):
        finite = state.loo[np.isfinite(state.loo)]
        best = float(np.max(finite)) if len(finite) else math.nan
        relative = best / max(state.cost, 1e-12) if np.isfinite(best) else math.nan
        return cls(state.cost, best, relative)


_GATE_KEYS = (
    "fraction_global_abnormal",
    "fraction_global_and_absolute_actionable",
    "fraction_global_and_relative_actionable",
)


class CleanGateCounter:
    """Constant-size joint gate counts for one operating point."""

    def __init__(self, thresholds):
        self.thresholds = thresholds
        self.count = 0
        self.passed = np.zeros(3, dtype=np.int64)

    def add(self, stats: CompactTraceStats):
        abnormal = bool(
            np.isfinite(stats.cost) and stats.cost > self.thresholds["global_threshold"]
        )
        self.count += 1
        self.passed += (
            abnormal,
            abnormal and stats.best_loo > self.thresholds["absolute_threshold"],
            abnormal
            and stats.best_relative_loo > self.thresholds["relative_threshold"],
        )

    def rates(self):
        return dict(
            zip(
                _GATE_KEYS,
                self.passed / self.count if self.count else np.full(3, np.nan),
            )
        )


def _positions(trace: pd.DataFrame, allow_uniform_fallback: bool) -> np.ndarray:
    if "Genomic_Position" in trace:
        positions = trace["Genomic_Position"].to_numpy(float)
        if np.isfinite(positions).all() and len(np.unique(positions)) == len(positions):
            return positions
        raise ValueError("Genomic_Position must be finite and unique within a trace")
    if not allow_uniform_fallback:
        raise ValueError(
            "Genomic_Position is required; use --assume-uniform-barcode-spacing "
            "only for the historical simulation"
        )
    # Barcode determines rank only, never genomic distance.
    barcode = trace["Barcode"].to_numpy()
    order = np.argsort(barcode, kind="stable")
    positions = np.empty(len(trace), float)
    positions[order] = np.arange(1, len(trace) + 1, dtype=float)
    return positions


def build_trace_scores(
    trace: pd.DataFrame,
    reference,
    reference_mode: str = "separation",
    allow_uniform_fallback: bool = False,
) -> TraceScores:
    """Build pair anomalies once and calculate exact top-3 LOO costs.

    Top-three values after deleting each node are selected from the cached
    sorted edges. No pair distance or reference lookup is repeated per policy.
    """
    if reference_mode not in {"separation", "barcode_pair"}:
        raise ValueError("reference_mode must be separation or barcode_pair")
    current = trace.copy().reset_index(drop=True)
    positions = _positions(current, allow_uniform_fallback)
    xyz = current[["x", "y", "z"]].to_numpy(float)
    barcodes = current["Barcode"].to_numpy()
    source = getattr(reference, reference_mode)
    edges = []
    for left in range(len(current)):
        for right in range(left + 1, len(current)):
            key = (
                float(abs(positions[right] - positions[left]))
                if reference_mode == "separation"
                else _pair_key(barcodes[left], barcodes[right])
            )
            sample = source.get(key)
            if sample is None or len(sample) < reference.minimum_observations:
                continue
            distance = float(np.linalg.norm(xyz[right] - xyz[left]))
            edges.append((left, right, empirical_anomaly(sample, distance)))
    edges.sort(key=lambda edge: edge[2], reverse=True)
    cost = float(sum(edge[2] for edge in edges[:3])) if edges else math.nan
    without = np.full(len(current), np.nan)
    if edges:
        for node in range(len(current)):
            kept = (edge[2] for edge in edges if node not in edge[:2])
            values = []
            for value in kept:
                values.append(value)
                if len(values) == 3:
                    break
            if values:
                without[node] = sum(values)
    return TraceScores(current, tuple(edges), cost, cost - without, without)


def calibrate_thresholds(clean_stats: Sequence[CompactTraceStats]) -> pd.DataFrame:
    """Calibrate unchanged empirical cutoffs from compact clean statistics."""
    samples = {
        "global": [stats.cost for stats in clean_stats],
        "absolute": [stats.best_loo for stats in clean_stats],
        "relative": [stats.best_relative_loo for stats in clean_stats],
    }
    rows = []
    for fpr in FIXED_FPRS:
        global_threshold = empirical_threshold(samples["global"], fpr)
        row = {
            "trace_fpr": fpr,
            "global_threshold": global_threshold,
            "absolute_threshold": empirical_threshold(samples["absolute"], fpr),
            "relative_threshold": empirical_threshold(samples["relative"], fpr),
            "observed_clean_global_fpr": float(
                np.mean(np.asarray(samples["global"]) > global_threshold)
            ),
            "n_clean_traces": len(clean_stats),
        }
        row.update(clean_gate_rates(clean_stats, row))
        rows.append(row)
    return pd.DataFrame(rows)


def clean_gate_rates(
    clean_stats: Sequence[CompactTraceStats], thresholds: Mapping[str, float]
) -> dict[str, float]:
    """Report joint gate rates without retaining flags or full trace states."""
    counter = CleanGateCounter(thresholds)
    for stats in clean_stats:
        counter.add(stats)
    return counter.rates()


def run_policy(
    trace: pd.DataFrame,
    reference,
    thresholds: Mapping[str, float],
    policy: str,
    *,
    max_removals: int = 3,
    min_remaining: int = 3,
    reference_mode: str = "separation",
    allow_uniform_fallback: bool = False,
    scorer=build_trace_scores,
) -> tuple[pd.DataFrame, dict]:
    """Run one policy and return its explainable audit and terminal metrics."""
    if policy not in POLICIES:
        raise ValueError(f"Unknown policy: {policy}")
    current = trace.copy()
    removed = []
    audit = []
    while True:
        state = scorer(current, reference, reference_mode, allow_uniform_fallback)
        abnormal = bool(
            np.isfinite(state.cost) and state.cost > thresholds["global_threshold"]
        )
        stop = ""
        if not np.isfinite(state.cost):
            stop = "insufficient_reference"
        elif not abnormal:
            stop = "global_cost_normal"
        elif len(current) <= min_remaining:
            stop = "minimum_trace_size_reached"
        elif policy == "C_cap" and len(removed) >= max_removals:
            stop = "max_removals_reached"
        finite = np.flatnonzero(np.isfinite(state.loo))
        best_score = float(np.max(state.loo[finite])) if len(finite) else math.nan
        best = finite[state.loo[finite] == best_score] if len(finite) else np.array([])
        n_tied_best = int(len(best))
        candidate_unique = n_tied_best == 1
        candidate = int(best[0]) if len(best) else None
        if not stop and candidate is None:
            stop = "no_scoreable_candidate"
        elif not stop and not candidate_unique:
            stop = "ambiguous_candidate"
        absolute = float(state.loo[candidate]) if candidate is not None else math.nan
        after = (
            float(state.cost_without[candidate]) if candidate is not None else math.nan
        )
        relative = (
            absolute / max(state.cost, 1e-12) if candidate is not None else math.nan
        )
        actionable = candidate is not None and candidate_unique
        if policy == "B_absolute":
            actionable = actionable and absolute > thresholds["absolute_threshold"]
        elif policy == "B_relative":
            actionable = actionable and relative > thresholds["relative_threshold"]
        if not stop and not actionable:
            stop = "candidate_not_actionable"
        accepted = not bool(stop)
        row = state.trace.iloc[candidate] if candidate is not None else None
        audit.append(
            {
                "Trace_ID": str(trace["Trace_ID"].iloc[0]),
                "policy": policy,
                "iteration": len(removed),
                "n_localizations_before": len(current),
                "C_before": state.cost,
                "global_threshold": thresholds["global_threshold"],
                "global_abnormal_before": abnormal,
                "candidate_Barcode": row["Barcode"] if row is not None else math.nan,
                "candidate_S_loo": absolute,
                "candidate_rank": 1 if candidate is not None else math.nan,
                "n_tied_best": n_tied_best,
                "candidate_unique": candidate_unique,
                "candidate_actionable": actionable,
                "C_after_candidate_removal": after,
                "absolute_C_drop": absolute,
                "relative_C_drop": relative,
                "removal_accepted": accepted,
                "candidate_is_corrupted": (
                    bool(row.get("is_corrupted", False)) if row is not None else False
                ),
                "stop_reason": stop,
            }
        )
        if not accepted:
            break
        spot_id = str(row["Spot_ID"])
        removed.append((spot_id, bool(row.get("is_corrupted", False))))
        current = current[current["Spot_ID"].astype(str) != spot_id].copy()
    n_true = int(trace.get("is_corrupted", False).astype(bool).sum())
    true_removed = sum(value for _, value in removed)
    false_removed = len(removed) - true_removed
    reason = audit[-1]["stop_reason"]
    terminal_cost = float(audit[-1]["C_before"])
    terminal_abnormal = bool(audit[-1]["global_abnormal_before"])
    terminal_global_state = (
        "unscoreable"
        if not np.isfinite(terminal_cost)
        else ("abnormal" if terminal_abnormal else "normal")
    )
    exact_recovery = len(removed) == n_true and true_removed == n_true
    false_clean = terminal_global_state == "normal" and n_true > true_removed
    overcurated = false_removed > 0
    abnormal_unresolved = terminal_abnormal and reason in {
        "candidate_not_actionable",
        "ambiguous_candidate",
        "max_removals_reached",
        "minimum_trace_size_reached",
        "insufficient_reference",
        "no_scoreable_candidate",
    }
    if overcurated:
        terminal_outcome = "overcurated"
    elif exact_recovery:
        terminal_outcome = "exact_recovery"
    elif false_clean:
        terminal_outcome = "false_clean"
    elif abnormal_unresolved:
        terminal_outcome = "abnormal_unresolved"
    else:
        terminal_outcome = "stopped"
    return pd.DataFrame(audit), {
        "Trace_ID": str(trace["Trace_ID"].iloc[0]),
        "policy": policy,
        "n_true_corruptions": n_true,
        "n_removed": len(removed),
        "true_removed": true_removed,
        "false_removed": false_removed,
        "true_corruption_recall": true_removed / n_true if n_true else math.nan,
        "removal_precision": true_removed / len(removed) if removed else math.nan,
        "exact_recovery": exact_recovery,
        "false_clean": false_clean,
        "overcurated": overcurated,
        "true_corruptions_remaining": n_true - true_removed,
        "unresolved_true_corruption": n_true > true_removed,
        "abnormal_unresolved": abnormal_unresolved,
        "terminal_C": terminal_cost,
        "terminal_global_abnormal": terminal_abnormal,
        "terminal_global_state": terminal_global_state,
        "terminal_outcome": terminal_outcome,
        "stop_reason": reason,
        "n_iterations": len(removed),
    }


def summarize(trace_metrics: pd.DataFrame, group: Sequence[str]) -> pd.DataFrame:
    rows = []
    for key, frame in trace_metrics.groupby(list(group), sort=False, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        rows.append(
            {
                **dict(zip(group, key)),
                "n_traces": len(frame),
                "true_corruption_recall": frame["true_removed"].sum()
                / max(frame["n_true_corruptions"].sum(), 1),
                "removal_precision": frame["true_removed"].sum()
                / max(frame["n_removed"].sum(), 1),
                "exact_recovery_rate": frame["exact_recovery"].mean(),
                "mean_genuine_localizations_removed": frame["false_removed"].mean(),
                "fraction_with_1plus_false_removal": (
                    frame["false_removed"] >= 1
                ).mean(),
                "fraction_with_2plus_false_removals": (
                    frame["false_removed"] >= 2
                ).mean(),
                "maximum_false_removals": frame["false_removed"].max(),
                "mean_true_corruptions_remaining": frame[
                    "true_corruptions_remaining"
                ].mean(),
                "fraction_unresolved_true_corruption": frame[
                    "unresolved_true_corruption"
                ].mean(),
                "fraction_abnormal_unresolved": frame["abnormal_unresolved"].mean(),
                "mean_removals": frame["n_removed"].mean(),
                "fraction_with_1plus_removal": (frame["n_removed"] >= 1).mean(),
                "fraction_with_2plus_removals": (frame["n_removed"] >= 2).mean(),
                "maximum_removals": frame["n_removed"].max(),
                "fraction_terminal_global_abnormal": frame[
                    "terminal_global_abnormal"
                ].mean(),
                "mean_iterations": frame["n_iterations"].mean(),
                "stop_reason_distribution": ";".join(
                    f"{key}:{int(value)}"
                    for key, value in frame["stop_reason"].value_counts().items()
                ),
                "removal_count_distribution": ";".join(
                    f"{int(k)}:{int(v)}"
                    for k, v in frame["n_removed"].value_counts().sort_index().items()
                ),
            }
        )
    return pd.DataFrame(rows)


def _reference_with_coordinates(
    cycle3, baseline: pd.DataFrame, minimum: int, mode: str
):
    if mode not in {"separation", "barcode_pair"}:
        raise ValueError("reference_mode must be separation or barcode_pair")
    values = {}
    for _, trace in baseline.groupby(["simulation_id", "Trace_ID"], sort=False):
        positions = _positions(trace, allow_uniform_fallback=False)
        xyz = trace[["x", "y", "z"]].to_numpy(float)
        barcode = trace["Barcode"].to_numpy()
        for left in range(len(trace)):
            for right in range(left + 1, len(trace)):
                key = (
                    float(abs(positions[right] - positions[left]))
                    if mode == "separation"
                    else _pair_key(barcode[left], barcode[right])
                )
                values.setdefault(key, []).append(
                    float(np.linalg.norm(xyz[right] - xyz[left]))
                )
    fitted = {key: np.sort(np.asarray(value, float)) for key, value in values.items()}
    return SimpleNamespace(
        separation=fitted if mode == "separation" else {},
        barcode_pair=fitted if mode == "barcode_pair" else {},
        minimum_observations=minimum,
    )


SUMMARY_GROUP = (
    "dataset_type",
    "policy",
    "n_barcodes",
    "K",
    "detection_efficiency",
    "displacement",
)
_SUM_FIELDS = (
    "true_removed",
    "n_true_corruptions",
    "n_removed",
    "exact_recovery",
    "false_removed",
    "true_corruptions_remaining",
    "unresolved_true_corruption",
    "abnormal_unresolved",
    "terminal_global_abnormal",
    "n_iterations",
)
_MEANS = {
    "exact_recovery_rate": "exact_recovery",
    "mean_genuine_localizations_removed": "false_removed",
    "mean_true_corruptions_remaining": "true_corruptions_remaining",
    "fraction_unresolved_true_corruption": "unresolved_true_corruption",
    "fraction_abnormal_unresolved": "abnormal_unresolved",
    "mean_removals": "n_removed",
    "fraction_terminal_global_abnormal": "terminal_global_abnormal",
    "mean_iterations": "n_iterations",
}


class ConditionAccumulator:
    """Exact integer sums/counts per condition; never retain trace rows."""

    def __init__(self):
        self.groups = {}

    def add(self, row):
        key = tuple(row[column] for column in SUMMARY_GROUP)
        if key not in self.groups:
            self.groups[key] = {
                "count": 0,
                "sums": Counter(),
                "stops": Counter(),
                "removals": Counter(),
                "false_removals": Counter(),
            }
        group = self.groups[key]
        group["count"] += 1
        group["sums"].update({field: int(row[field]) for field in _SUM_FIELDS})
        group["stops"][row["stop_reason"]] += 1
        group["removals"][int(row["n_removed"])] += 1
        group["false_removals"][int(row["false_removed"])] += 1

    def frame(self):
        rows = []
        for key, group in self.groups.items():
            n, sums = group["count"], group["sums"]
            row = {
                **dict(zip(SUMMARY_GROUP, key)),
                "n_traces": n,
                "true_corruption_recall": sums["true_removed"]
                / max(sums["n_true_corruptions"], 1),
                "removal_precision": sums["true_removed"] / max(sums["n_removed"], 1),
                **{name: sums[field] / n for name, field in _MEANS.items()},
            }
            for field, distribution in (
                ("false_removal", group["false_removals"]),
                ("removal", group["removals"]),
            ):
                row[f"fraction_with_1plus_{field}"] = (
                    sum(v for k, v in distribution.items() if k >= 1) / n
                )
                row[f"fraction_with_2plus_{field}s"] = (
                    sum(v for k, v in distribution.items() if k >= 2) / n
                )
                row[f"maximum_{field}s"] = max(distribution)
            row["stop_reason_distribution"] = ";".join(
                f"{k}:{v}" for k, v in group["stops"].most_common()
            )
            row["removal_count_distribution"] = ";".join(
                f"{k}:{v}" for k, v in sorted(group["removals"].items())
            )
            rows.append(row)
        columns = [
            *SUMMARY_GROUP,
            "n_traces",
            "true_corruption_recall",
            "removal_precision",
            "exact_recovery_rate",
            "mean_genuine_localizations_removed",
            "fraction_with_1plus_false_removal",
            "fraction_with_2plus_false_removals",
            "maximum_false_removals",
            "mean_true_corruptions_remaining",
            "fraction_unresolved_true_corruption",
            "fraction_abnormal_unresolved",
            "mean_removals",
            "fraction_with_1plus_removal",
            "fraction_with_2plus_removals",
            "maximum_removals",
            "fraction_terminal_global_abnormal",
            "mean_iterations",
            "stop_reason_distribution",
            "removal_count_distribution",
        ]
        return pd.DataFrame(rows, columns=columns)


class CostTrajectories:
    """Finite cost sums/counts by policy and iteration for the trajectory plot."""

    def __init__(self):
        self.values = {}

    def add(self, audit):
        corrupted = audit[audit["dataset_type"] == "corrupted_evaluation"]
        for row in corrupted.itertuples():
            key = (row.policy, int(row.iteration))
            total, count = self.values.get(key, (0.0, 0))
            if pd.notna(row.C_before):
                total += row.C_before
                count += 1
            self.values[key] = (total, count)

    def frame(self):
        return pd.DataFrame(
            [
                {
                    "policy": policy,
                    "iteration": iteration,
                    "C_before": total / count if count else math.nan,
                }
                for (policy, iteration), (total, count) in self.values.items()
            ],
            columns=["policy", "iteration", "C_before"],
        )


class BufferedEcsvWriter:
    """Bounded row buffer around the earlier benchmark's chunk writer."""

    def __init__(self, writer, max_rows=1024):
        self.writer = writer
        self.max_rows = max_rows
        self.pending = []
        self.barcode_numeric = None

    def write(self, frame):
        for row in frame.to_dict("records"):
            self.pending.append(row)
            if len(self.pending) >= self.max_rows:
                self.flush()

    def flush(self):
        if self.pending:
            frame = pd.DataFrame(self.pending)
            # These audit fields can become NaN after a scoreable first chunk.
            for column in ("candidate_rank", "candidate_Barcode"):
                if column in frame and pd.api.types.is_numeric_dtype(frame[column]):
                    frame[column] = frame[column].astype(float)
            if "candidate_Barcode" in frame and self.barcode_numeric is False:
                frame["candidate_Barcode"] = (
                    frame["candidate_Barcode"].fillna("").astype(str)
                )
            self.writer.write(frame)
            self.pending.clear()

    def finish(self):
        self.flush()
        self.writer.finish()


def benchmark(
    root: Path,
    output: Path,
    minimum_observations: int = 20,
    reference_mode: str = "separation",
    assume_uniform: bool = False,
    max_removals: int = 3,
    min_remaining: int = 3,
    profile_traces: int = 0,
    profile_calibration_traces: int = 0,
) -> None:
    """Run the cycle-3 folds with bounded evaluation/output retention.

    Input tables and exact empirical references remain group/condition sized;
    calibration retains three scalars per clean trace, never scoring structures.
    """
    global _PROGRESS_CONTEXT, _GROUP_PEAK_MB
    started = time.monotonic()
    rss_start = _rss_mb()
    cycle3 = _cycle3()
    cycle1 = cycle3._cycle(1)
    manifest = yaml.safe_load((root / "sweep_manifest.yaml").read_text())
    conditions = manifest.get("conditions", [])
    if manifest.get("dry_run") or not conditions:
        raise ValueError("A completed cycle-3 sweep manifest is required")
    output.mkdir(parents=True, exist_ok=True)
    trace_writer = BufferedEcsvWriter(cycle1._EcsvChunkWriter(output / OUTPUTS[0]))
    audit_writer = BufferedEcsvWriter(cycle1._EcsvChunkWriter(output / OUTPUTS[2]))
    threshold_writer = cycle1._EcsvChunkWriter(output / OUTPUTS[3])
    aggregates = ConditionAccumulator()
    trajectories = CostTrajectories()
    timings = dict.fromkeys(
        (
            "reference_fitting_seconds",
            "calibration_scoring_seconds",
            "clean_evaluation_seconds",
            "corrupted_condition_scoring_seconds",
            "plotting_writing_seconds",
        ),
        0.0,
    )
    build_trace_scores_calls = 0
    reference_cache_hits = reference_cache_misses = 0
    runtime_rows = []

    def score_trace(trace, reference, category):
        nonlocal build_trace_scores_calls
        score_started = time.monotonic()
        result = build_trace_scores(trace, reference, reference_mode, assume_uniform)
        timings[category] += time.monotonic() - score_started
        build_trace_scores_calls += 1
        _rss_mb()
        return result

    def evaluate_trace(
        identity, trace, reference, primary, meta, category, initial_state=None
    ):
        # This function boundary prevents score caches/states escaping a trace.
        cache = {}
        if initial_state is not None:
            cache[tuple(sorted(trace["Spot_ID"].astype(str)))] = initial_state

        def shared_scorer(current, model, mode, fallback):
            key = tuple(sorted(current["Spot_ID"].astype(str)))
            if key not in cache:
                cache[key] = score_trace(current, model, category)
            return cache[key]

        for policy in POLICIES:
            audit, metrics = run_policy(
                trace,
                reference,
                primary,
                policy,
                max_removals=max_removals,
                min_remaining=min_remaining,
                reference_mode=reference_mode,
                allow_uniform_fallback=assume_uniform,
                scorer=shared_scorer,
            )
            write_started = time.monotonic()
            audit = audit.assign(simulation_id=identity[0], **meta)
            record = {"simulation_id": identity[0], **meta, **metrics}
            audit_writer.barcode_numeric = pd.api.types.is_numeric_dtype(
                trace["Barcode"]
            )
            audit_writer.write(audit)
            trace_writer.write(pd.DataFrame([record]))
            aggregates.add(record)
            trajectories.add(audit)
            timings["plotting_writing_seconds"] += time.monotonic() - write_started

    groups = {}
    for item in conditions:
        n, _ = cycle3._condition_metadata(item)
        groups.setdefault((n, float(item["detection_efficiency"])), []).append(item)
    processed = 0
    for (n_barcodes, efficiency), items in groups.items():
        group_started = time.monotonic()
        _PROGRESS_CONTEXT = f"group B={n_barcodes} p={efficiency:g}"
        _GROUP_PEAK_MB = 0.0
        group_rss_start = _rss_mb()
        group_processed_start = processed
        group_hits_start, group_misses_start = (
            reference_cache_hits,
            reference_cache_misses,
        )
        _progress("start group")
        # Load only this efficiency's baselines, preserving the existing
        # n_barcodes-column selection and uniform-coordinate mapping semantics.
        parts = []
        for item in manifest["simulations"]:
            if not np.isclose(float(item["detection_efficiency"]), efficiency):
                continue
            frame = cycle1.load_baseline(root / "simulations" / str(item["id"]), item)
            if (
                int(frame.get("n_barcodes", pd.Series([n_barcodes])).iloc[0])
                != n_barcodes
            ):
                continue
            if "Genomic_Position" not in frame and assume_uniform:
                mapping = {
                    b: i + 1 for i, b in enumerate(sorted(frame["Barcode"].unique()))
                }
                frame["Genomic_Position"] = frame["Barcode"].map(mapping)
            parts.append(frame)
        baseline = pd.concat(parts)
        del parts, frame
        seeds = sorted({int(item["seed"]) for item in items})
        reference_cache = {}

        def fit_reference(seeds, baseline, label):
            nonlocal reference_cache_hits, reference_cache_misses
            key = (reference_mode, tuple(sorted(map(int, seeds))))
            if key in reference_cache:
                reference_cache_hits += 1
                _progress(f"reference cache hit: {label}; seeds={key[-1]}")
                return reference_cache[key]
            reference_cache_misses += 1
            fit_started = time.monotonic()
            _progress(f"{label} reference fit start; seeds={key[-1]}")
            model = _reference_with_coordinates(
                cycle3,
                baseline[baseline["seed"].isin(seeds)],
                minimum_observations,
                reference_mode,
            )
            fit_elapsed = time.monotonic() - fit_started
            timings["reference_fitting_seconds"] += fit_elapsed
            reference_cache[key] = model
            _progress(f"{label} reference fit end; elapsed={fit_elapsed:.3f}s")
            return model

        for evaluation_seed in seeds:
            evaluation_started = time.monotonic()
            fold_context = (
                f"group B={n_barcodes} p={efficiency:g} eval={evaluation_seed}"
            )
            _PROGRESS_CONTEXT = fold_context
            _progress("start evaluation seed")
            ref_seeds, _, _ = fold_seed_sets(seeds, evaluation_seed)
            reference = fit_reference(ref_seeds, baseline, "outer")
            clean_stats = []
            calibration_seeds = []
            nested_reference_sets = []
            for calibration_seed in seeds:
                if calibration_seed == evaluation_seed:
                    continue
                _PROGRESS_CONTEXT = f"{fold_context} calibration={calibration_seed}"
                _progress("start calibration seed")
                nested, _, _ = fold_seed_sets(seeds, evaluation_seed, calibration_seed)
                calibration_seeds.append(calibration_seed)
                nested_reference_sets.append(
                    f"{calibration_seed}:{','.join(map(str, sorted(nested)))}"
                )
                nested_ref = fit_reference(nested, baseline, "nested")
                calibration_groups = baseline[
                    baseline["seed"] == calibration_seed
                ].groupby(["simulation_id", "Trace_ID"], sort=False)
                if profile_calibration_traces:
                    calibration_groups = islice(
                        calibration_groups, profile_calibration_traces
                    )
                calibration_index = 0
                for calibration_index, (_, clean) in enumerate(
                    calibration_groups, start=1
                ):
                    clean_stats.append(
                        CompactTraceStats.from_scores(
                            score_trace(
                                clean, nested_ref, "calibration_scoring_seconds"
                            )
                        )
                    )
                    del clean
                    if calibration_index % 250 == 0:
                        _progress(f"clean calibration traces={calibration_index}")
                _progress(f"end calibration seed; traces={calibration_index}")
                del nested_ref, calibration_groups
            thresholds = calibrate_thresholds(clean_stats)
            del clean_stats
            threshold_records = thresholds.to_dict("records")
            gates = [CleanGateCounter(row) for row in threshold_records]
            primary = (
                thresholds[np.isclose(thresholds["trace_fpr"], 0.01)].iloc[0].to_dict()
            )
            clean_meta = {
                "condition": f"clean-evaluation-{n_barcodes}-{efficiency:g}-{evaluation_seed}",
                "seed": evaluation_seed,
                "n_barcodes": n_barcodes,
                "K": 0,
                "detection_efficiency": efficiency,
                "displacement": 0.0,
                "dataset_type": "clean_evaluation",
            }
            _PROGRESS_CONTEXT = f"{fold_context} clean_evaluation"
            evaluation_groups = baseline[baseline["seed"] == evaluation_seed].groupby(
                ["simulation_id", "Trace_ID"], sort=False
            )
            clean_index = 0
            for clean_index, (identity, clean) in enumerate(evaluation_groups, start=1):
                # Gate evaluation covers all held-out traces, as before. Once
                # the profiling policy cap is reached, only the gates continue.
                state = score_trace(clean, reference, "clean_evaluation_seconds")
                stats = CompactTraceStats.from_scores(state)
                for gate in gates:
                    gate.add(stats)
                if not profile_traces or processed < profile_traces:
                    evaluate_trace(
                        identity,
                        clean,
                        reference,
                        primary,
                        clean_meta,
                        "clean_evaluation_seconds",
                        state,
                    )
                    processed += 1
                del state, clean
                if clean_index % 250 == 0:
                    _progress(
                        f"held-out clean traces={clean_index}; policy traces={processed}"
                    )
            _progress(
                f"end held-out clean evaluation; traces={clean_index}; policy traces={processed}"
            )
            del evaluation_groups
            threshold_writer.write(
                pd.DataFrame(
                    [
                        {
                            "n_barcodes": n_barcodes,
                            "detection_efficiency": efficiency,
                            "evaluation_seed": evaluation_seed,
                            "reference_seeds": ",".join(map(str, sorted(ref_seeds))),
                            "calibration_seeds": ",".join(
                                map(str, sorted(calibration_seeds))
                            ),
                            "nested_reference_seed_sets": ";".join(
                                nested_reference_sets
                            ),
                            "evaluation_seeds": str(evaluation_seed),
                            **{
                                f"evaluation_{key}": value
                                for key, value in gate.rates().items()
                            },
                            **threshold,
                        }
                        for threshold, gate in zip(threshold_records, gates)
                    ]
                )
            )
            if not profile_traces or processed < profile_traces:
                for item in items:
                    if int(item["seed"]) != evaluation_seed:
                        continue
                    condition_id = str(
                        item.get("id", item.get("condition_id", "unknown"))
                    )
                    _PROGRESS_CONTEXT = f"{fold_context} condition={condition_id}"
                    _progress("start corrupted condition")
                    observations = cycle3._load_observations(cycle1, root, item)
                    if "Genomic_Position" not in observations and assume_uniform:
                        mapping = {
                            b: i + 1
                            for i, b in enumerate(
                                sorted(observations["Barcode"].unique())
                            )
                        }
                        observations["Genomic_Position"] = observations["Barcode"].map(
                            mapping
                        )
                    condition_groups = observations.groupby(
                        ["simulation_id", "Trace_ID"], sort=False
                    )
                    meta = {
                        **cycle3._trace_metadata(item),
                        "dataset_type": "corrupted_evaluation",
                    }
                    condition_index = 0
                    for identity, trace in condition_groups:
                        if profile_traces and processed >= profile_traces:
                            del trace
                            break
                        evaluate_trace(
                            identity,
                            trace,
                            reference,
                            primary,
                            meta,
                            "corrupted_condition_scoring_seconds",
                        )
                        del trace
                        processed += 1
                        condition_index += 1
                        if condition_index % 250 == 0:
                            _progress(
                                f"corrupted traces={condition_index}; policy traces={processed}"
                            )
                    _progress(
                        f"end corrupted condition; traces={condition_index}; policy traces={processed}"
                    )
                    del condition_groups, observations
                    if profile_traces and processed >= profile_traces:
                        break
            _PROGRESS_CONTEXT = fold_context
            _progress(
                f"end evaluation seed; elapsed={time.monotonic() - evaluation_started:.3f}s"
            )
            del reference
            if profile_traces and processed >= profile_traces:
                break
        # Drop both cache ownership and the input slice at every boundary,
        # including profiling exits. No per-trace caches exist in this scope.
        cache_size = len(reference_cache)
        reference_cache.clear()
        del baseline
        _PROGRESS_CONTEXT = f"group B={n_barcodes} p={efficiency:g}"
        group_rss_end = _rss_mb()
        runtime_rows.append(
            {
                "scope": "group",
                "n_barcodes": n_barcodes,
                "detection_efficiency": efficiency,
                "profile_traces": processed - group_processed_start,
                "elapsed_seconds": time.monotonic() - group_started,
                "rss_start_mb": group_rss_start,
                "rss_end_mb": group_rss_end,
                "rss_peak_mb": _GROUP_PEAK_MB,
                "reference_cache_hits": reference_cache_hits - group_hits_start,
                "reference_cache_misses": reference_cache_misses - group_misses_start,
                "reference_cache_size_before_clear": cache_size,
                "reference_cache_size_after_clear": len(reference_cache),
            }
        )
        _progress(
            f"end group; elapsed={time.monotonic() - group_started:.3f}s; reference cache cleared"
        )
        if profile_traces and processed >= profile_traces:
            break
    output_started = time.monotonic()
    trace_writer.finish()
    audit_writer.finish()
    threshold_writer.finish()
    _write(aggregates.frame(), output / OUTPUTS[1])
    plot_results(output, trajectories)
    timings["plotting_writing_seconds"] += time.monotonic() - output_started
    elapsed = time.monotonic() - started
    runtime_rows.append(
        {
            "scope": "total",
            "profile_traces": processed,
            "profile_calibration_traces": profile_calibration_traces,
            "elapsed_seconds": elapsed,
            "seconds_per_trace_all_policies": elapsed / max(processed, 1),
            "seconds_per_1000_traces": elapsed / max(processed, 1) * 1000,
            "projected_hours_for_100000_traces": elapsed
            / max(processed, 1)
            * 100000
            / 3600,
            **timings,
            "build_trace_scores_calls": build_trace_scores_calls,
            "reference_cache_hits": reference_cache_hits,
            "reference_cache_misses": reference_cache_misses,
            "rss_start_mb": rss_start,
            "rss_end_mb": _rss_mb(),
            "rss_peak_mb": max(
                [_GROUP_PEAK_MB, *(row["rss_peak_mb"] for row in runtime_rows)]
            ),
        }
    )
    _write(pd.DataFrame(runtime_rows), output / OUTPUTS[4])
    _PROGRESS_CONTEXT = ""
    _progress(
        "runtime summary: "
        + ", ".join(f"{key}={value:.3f}s" for key, value in timings.items())
        + f", build_trace_scores_calls={build_trace_scores_calls}, "
        f"reference_cache_hits={reference_cache_hits}, reference_cache_misses={reference_cache_misses}"
    )


def plot_results(output: Path, trajectories: CostTrajectories | None = None) -> None:
    condition = Table.read(output / OUTPUTS[1], format="ascii.ecsv").to_pandas()
    if trajectories is None:
        trajectories = CostTrajectories()
        # ECSV data uses a space delimiter with quoted strings. Replotting must
        # also remain bounded; never read the completed audit through Astropy.
        for chunk in pd.read_csv(
            output / OUTPUTS[2],
            comment="#",
            sep=" ",
            chunksize=4096,
            keep_default_na=False,
            na_values=["nan"],
            usecols=["dataset_type", "policy", "iteration", "C_before"],
        ):
            trajectories.add(chunk)
    performance = condition[condition["dataset_type"] == "corrupted_evaluation"]
    specs = [
        ("true_corruption_recall", "recall_vs_displacement.png"),
        ("removal_precision", "precision_vs_displacement.png"),
        ("exact_recovery_rate", "exact_recovery_vs_displacement.png"),
        ("mean_genuine_localizations_removed", "false_removals_vs_displacement.png"),
        ("fraction_abnormal_unresolved", "abnormal_unresolved_vs_displacement.png"),
        ("mean_iterations", "iterations_vs_displacement.png"),
    ]
    for metric, name in specs:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for policy, frame in performance.groupby("policy", sort=False):
            line = frame.groupby("displacement")[metric].mean().sort_index()
            ax.plot(line.index, line.values, marker="o", label=policy)
        ax.set(xlabel="Displacement (µm)", ylabel=metric)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(output / name, dpi=160)
        plt.close(fig)
    clean = condition[condition["dataset_type"] == "clean_evaluation"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if len(clean):
        values = clean.groupby("policy")["fraction_with_1plus_false_removal"].mean()
        ax.bar(values.index, values.values)
        ax.tick_params(axis="x", rotation=20)
    ax.set(ylabel="Clean traces with any removal")
    fig.tight_layout()
    fig.savefig(output / "clean_any_removal.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    corrupted_audit = trajectories.frame()
    for policy, frame in corrupted_audit.groupby("policy", sort=False):
        line = frame.groupby("iteration")["C_before"].mean().sort_index()
        ax.plot(line.index, line.values, marker="o", label=policy)
    ax.set(xlabel="Iteration", ylabel="Mean C(T)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output / "cost_trajectories.png", dpi=160)
    plt.close(fig)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-root", type=Path, default=Path("benchmark_curator_loo_validation")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-reference-observations", type=int, default=20)
    parser.add_argument(
        "--reference-mode", choices=("separation", "barcode_pair"), default="separation"
    )
    parser.add_argument("--assume-uniform-barcode-spacing", action="store_true")
    parser.add_argument("--max-removals", type=int, default=3)
    parser.add_argument("--min-remaining-localizations", type=int, default=3)
    parser.add_argument("--profile-traces", type=int, default=0)
    parser.add_argument(
        "--profile-calibration-traces",
        type=int,
        default=0,
        help="profiling-only deterministic per-seed cap on clean calibration traces",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_arguments()
    benchmark(
        args.benchmark_root.resolve(),
        args.output_dir.resolve(),
        args.minimum_reference_observations,
        args.reference_mode,
        args.assume_uniform_barcode_spacing,
        args.max_removals,
        args.min_remaining_localizations,
        args.profile_traces,
        args.profile_calibration_traces,
    )


if __name__ == "__main__":
    main()
