#!/usr/bin/env python3
"""Cycle-5 targeted validation of exact LOO-top3 ties, using cycle-3 data.

This is exploratory analysis. It does not change the production API or cycle 4.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from collections import Counter, OrderedDict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def _cycle4():
    path = Path(__file__).with_name("analyze_trace_curator_iterative_stopping.py")
    spec = importlib.util.spec_from_file_location("_tie_resolution_cycle4", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


C4 = _cycle4()
METHODS = ("stop_unresolved", "loo_top5", "k1_local", "one_step_lookahead")
POLICIES = dict(
    zip(METHODS, ("C_cap_stop", "C_cap_top5", "C_cap_k1", "C_cap_lookahead"))
)
EVENT_KEYS = ("dataset_type", "K", "detection_efficiency", "displacement")
OUTPUTS = (
    "cycle5_tie_events.ecsv",
    "cycle5_tie_method_scores.ecsv",
    "cycle5_tie_resolution_performance.ecsv",
    "cycle5_iterative_performance.ecsv",
    "cycle5_clean_control_performance.ecsv",
    "cycle5_k1_control_performance.ecsv",
    "cycle5_runtime_profile.ecsv",
)


def primary_candidates(state):
    """Exact maximal finite LOO-top3 scores; no tolerance or row-order winner."""
    finite = np.flatnonzero(np.isfinite(state.loo))
    if not len(finite):
        return np.array([], dtype=int)
    return finite[state.loo[finite] == np.max(state.loo[finite])]


@dataclass(frozen=True)
class Resolution:
    candidates: tuple[int, ...]
    scores: tuple[float, ...]
    winner: int | None
    left: tuple[float, ...]
    right: tuple[float, ...]
    left_spots: tuple[str, ...]
    right_spots: tuple[str, ...]

    @property
    def n_scoreable(self):
        return int(np.isfinite(self.scores).sum())

    @property
    def status(self):
        if not self.n_scoreable:
            return "unscoreable"
        return (
            "scoreable"
            if self.n_scoreable == len(self.candidates)
            else "partially_scoreable"
        )


def resolve_tie(
    state, method, scorer, reference, reference_mode="separation", assume_uniform=False
):
    """Rank only the primary tied set. Labels never enter secondary scoring."""
    if method not in METHODS:
        raise ValueError(f"Unknown tie method: {method}")
    tied = primary_candidates(state)
    if len(tied) < 2:
        raise ValueError("Secondary scoring requires an exact primary LOO-top3 tie")
    scores = np.full(len(tied), np.nan)
    left, right = scores.copy(), scores.copy()
    left_spots, right_spots = [""] * len(tied), [""] * len(tied)
    if method == "loo_top5":
        edges = sorted(
            (edge for edge in state.edges if np.isfinite(edge[2])),
            key=lambda edge: edge[2],
            reverse=True,
        )
        before = sum(edge[2] for edge in edges[:5]) if edges else math.nan
        for offset, node in enumerate(tied):
            kept = [edge[2] for edge in edges if node not in edge[:2]][:5]
            if kept:
                scores[offset] = before - sum(kept)
    elif method == "k1_local":
        positions = C4._positions(state.trace, assume_uniform)
        # These edges already contain the exact empirical anomaly from C4's
        # reference mode, minimum-support check, and genomic separation.
        edge_scores = {tuple(sorted((a, b))): value for a, b, value in state.edges}
        for offset, node in enumerate(tied):
            available_left = np.flatnonzero(positions < positions[node])
            available_right = np.flatnonzero(positions > positions[node])
            if len(available_left):
                neighbor = int(available_left[np.argmax(positions[available_left])])
                left[offset] = edge_scores.get(
                    tuple(sorted((int(node), neighbor))), math.nan
                )
                left_spots[offset] = str(state.trace.iloc[neighbor]["Spot_ID"])
            if len(available_right):
                neighbor = int(available_right[np.argmin(positions[available_right])])
                right[offset] = edge_scores.get(
                    tuple(sorted((int(node), neighbor))), math.nan
                )
                right_spots[offset] = str(state.trace.iloc[neighbor]["Spot_ID"])
            flanks = np.array([left[offset], right[offset]])
            finite = flanks[np.isfinite(flanks)]
            if len(finite):
                scores[offset] = float(np.mean(finite))
    elif method == "one_step_lookahead":
        for offset, node in enumerate(tied):
            spot = str(state.trace.iloc[node]["Spot_ID"])
            reduced = state.trace[state.trace["Spot_ID"].astype(str) != spot]
            rescored = scorer(reduced, reference, reference_mode, assume_uniform)
            finite = rescored.loo[np.isfinite(rescored.loo)]
            if len(finite):
                scores[offset] = float(np.max(finite))
    finite = np.flatnonzero(np.isfinite(scores))
    best = finite[scores[finite] == np.max(scores[finite])] if len(finite) else []
    winner = int(tied[best[0]]) if len(best) == 1 else None
    return Resolution(
        tuple(map(int, tied)),
        tuple(map(float, scores)),
        winner,
        tuple(left),
        tuple(right),
        tuple(left_spots),
        tuple(right_spots),
    )


class TraceScorer:
    """Bounded, trace-local LRU; provisional lookahead states never escape it."""

    def __init__(self, builder=None, max_states=64):
        self.builder = builder or C4.build_trace_scores
        self.max_states = max_states
        self.cache = OrderedDict()
        self.last_state = None
        self.calls = 0
        self.hits = 0
        self.peak_states = 0

    def __call__(self, trace, reference, mode, fallback):
        key = tuple(sorted(trace["Spot_ID"].astype(str)))
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
        else:
            self.calls += 1
            self.cache[key] = self.builder(trace, reference, mode, fallback)
            if len(self.cache) > self.max_states:
                self.cache.popitem(last=False)
            self.peak_states = max(self.peak_states, len(self.cache))
        self.last_state = self.cache[key]
        return self.last_state

    def clear(self):
        self.cache.clear()
        self.last_state = None


class ProfileLimit(Exception):
    """Abort an incomplete profile trace rather than biasing one policy variant."""


class EventBudget:
    def __init__(self, limit=0):
        self.limit = limit
        self.used = 0

    def consume(self):
        if self.limit and self.used >= self.limit:
            raise ProfileLimit()
        self.used += 1


class TieSession:
    """Compact secondary results for at most the states of one capped trace."""

    def __init__(self, scorer, reference, mode, fallback, budget):
        self.scorer, self.reference = scorer, reference
        self.mode, self.fallback, self.budget = mode, fallback, budget
        self.results = {}

    def resolve(self, state, method):
        key = tuple(sorted(state.trace["Spot_ID"].astype(str)))
        if key not in self.results:
            self.budget.consume()
            self.results[key] = {}
        if method not in self.results[key]:
            self.results[key][method] = resolve_tie(
                state, method, self.scorer, self.reference, self.mode, self.fallback
            )
        return self.results[key][method]


def terminal_metrics(trace, state, removed, reason, global_threshold):
    """Cycle-4 terminal accounting, with the same reason precedence/outcomes."""
    n_true = int(trace["is_corrupted"].astype(bool).sum())
    true_removed = sum(truth for _, truth in removed)
    false_removed = len(removed) - true_removed
    abnormal = bool(np.isfinite(state.cost) and state.cost > global_threshold)
    global_state = (
        "unscoreable"
        if not np.isfinite(state.cost)
        else ("abnormal" if abnormal else "normal")
    )
    exact = len(removed) == n_true and true_removed == n_true
    false_clean = global_state == "normal" and n_true > true_removed
    unresolved = bool(
        abnormal
        and reason
        in {
            "candidate_not_actionable",
            "ambiguous_candidate",
            "max_removals_reached",
            "minimum_trace_size_reached",
            "insufficient_reference",
            "no_scoreable_candidate",
        }
    )
    outcome = (
        "overcurated"
        if false_removed
        else (
            "exact_recovery"
            if exact
            else (
                "false_clean"
                if false_clean
                else "abnormal_unresolved" if unresolved else "stopped"
            )
        )
    )
    return {
        "Trace_ID": str(trace["Trace_ID"].iloc[0]),
        "policy": "C_cap",
        "n_true_corruptions": n_true,
        "n_removed": len(removed),
        "true_removed": true_removed,
        "false_removed": false_removed,
        "true_corruption_recall": true_removed / n_true if n_true else math.nan,
        "removal_precision": true_removed / len(removed) if removed else math.nan,
        "exact_recovery": exact,
        "false_clean": false_clean,
        "overcurated": false_removed > 0,
        "true_corruptions_remaining": n_true - true_removed,
        "unresolved_true_corruption": n_true > true_removed,
        "abnormal_unresolved": unresolved,
        "terminal_C": float(state.cost),
        "terminal_global_abnormal": bool(abnormal),
        "terminal_global_state": global_state,
        "terminal_outcome": outcome,
        "stop_reason": reason,
        "n_iterations": len(removed),
    }


def run_integrated(
    trace,
    reference,
    thresholds,
    method,
    *,
    scorer=None,
    resolver=None,
    reference_mode="separation",
    assume_uniform=False,
):
    """C_cap with secondary ranking only at an otherwise actionable exact tie."""
    scorer = scorer or TraceScorer()
    if method not in METHODS:
        raise ValueError(f"Unknown tie method: {method}")
    if method == "stop_unresolved":
        audit, metrics = C4.run_policy(
            trace,
            reference,
            thresholds,
            "C_cap",
            max_removals=3,
            min_remaining=3,
            reference_mode=reference_mode,
            allow_uniform_fallback=assume_uniform,
            scorer=scorer,
        )
        invocations = int(metrics["stop_reason"] == "ambiguous_candidate")
        return metrics, {
            "tie_invocations": invocations,
            "tie_resolved": 0,
            "tie_correct": 0,
            "had_primary_tie": bool(invocations),
        }
    resolver = resolver or (
        lambda state, selected: resolve_tie(
            state, selected, scorer, reference, reference_mode, assume_uniform
        )
    )
    current, removed = trace.copy(), []
    counts = {
        "tie_invocations": 0,
        "tie_resolved": 0,
        "tie_correct": 0,
        "had_primary_tie": False,
    }
    while True:
        state = scorer(current, reference, reference_mode, assume_uniform)
        abnormal = bool(
            np.isfinite(state.cost) and state.cost > thresholds["global_threshold"]
        )
        reason = ""
        if not np.isfinite(state.cost):
            reason = "insufficient_reference"
        elif not abnormal:
            reason = "global_cost_normal"
        elif len(current) <= 3:
            reason = "minimum_trace_size_reached"
        elif len(removed) >= 3:
            reason = "max_removals_reached"
        if reason:
            break
        tied = primary_candidates(state)
        if not len(tied):
            reason = "no_scoreable_candidate"
            break
        candidate = int(tied[0]) if len(tied) == 1 else None
        if len(tied) >= 2:
            counts["tie_invocations"] += 1
            counts["had_primary_tie"] = True
            resolution = resolver(state, method)
            candidate = resolution.winner
            if candidate is None:
                reason = "ambiguous_candidate"
                break
            if candidate not in tied:
                raise ValueError("Secondary winner must belong to the primary tied set")
            counts["tie_resolved"] += 1
            counts["tie_correct"] += int(
                bool(state.trace.iloc[candidate]["is_corrupted"])
            )
        row = state.trace.iloc[candidate]
        spot = str(row["Spot_ID"])
        removed.append((spot, bool(row["is_corrupted"])))
        current = current[current["Spot_ID"].astype(str) != spot].copy()
    return (
        terminal_metrics(trace, state, removed, reason, thresholds["global_threshold"]),
        counts,
    )


def tie_rows(state, iteration, thresholds, meta, session):
    """Offline events reached by the conservative cycle-4 path only."""
    tied = primary_candidates(state)
    truth = state.trace["is_corrupted"].astype(bool).to_numpy()
    n_true = int(truth[tied].sum())
    tie_class = (
        "no_true_candidate"
        if not n_true
        else (
            "all_tied_candidates_true"
            if n_true == len(tied)
            else "single_true_candidate_in_tie" if n_true == 1 else "mixed_tie"
        )
    )
    event = {
        **meta,
        "Trace_ID": str(state.trace["Trace_ID"].iloc[0]),
        "iteration": iteration,
        "n_localizations": len(state.trace),
        "C_top3": state.cost,
        "global_threshold": thresholds["global_threshold"],
        "tie_size": len(tied),
        "tied_Spot_IDs": json.dumps(
            state.trace.iloc[tied]["Spot_ID"].astype(str).tolist()
        ),
        "tied_Barcodes": json.dumps(
            state.trace.iloc[tied]["Barcode"].astype(str).tolist()
        ),
        "n_true_corruptions_remaining": int(truth.sum()),
        "n_true_candidates": n_true,
        "tie_class": tie_class,
        "oracle_resolvable": n_true > 0,
        "all_tied_candidates_true": n_true == len(tied),
        "mixed_tie": 0 < n_true < len(tied),
        "single_true_candidate_in_tie": n_true == 1,
    }
    rows, summaries = [], []
    for method in METHODS:
        result = session.resolve(state, method)
        selected = (
            state.trace.iloc[result.winner] if result.winner is not None else None
        )
        decision = {
            "method": method,
            "method_scoreable": bool(result.n_scoreable),
            "method_scoreability": (
                "conservative_stop" if method == "stop_unresolved" else result.status
            ),
            "n_scoreable_candidates": result.n_scoreable,
            "unique_winner": result.winner is not None,
            "selected_Spot_ID": (
                str(selected["Spot_ID"]) if selected is not None else ""
            ),
            "selected_Barcode": (
                str(selected["Barcode"]) if selected is not None else ""
            ),
            "selected_is_true_corruption": (
                bool(selected["is_corrupted"]) if selected is not None else False
            ),
        }
        summaries.append({**event, **decision})
        for offset, node in enumerate(result.candidates):
            candidate = state.trace.iloc[node]
            rows.append(
                {
                    **event,
                    **decision,
                    "candidate_Spot_ID": str(candidate["Spot_ID"]),
                    "candidate_Barcode": str(candidate["Barcode"]),
                    "candidate_is_true_corruption": bool(candidate["is_corrupted"]),
                    "secondary_score": float(result.scores[offset]),
                    "candidate_scoreable": bool(np.isfinite(result.scores[offset])),
                    "left_anomaly": float(result.left[offset]),
                    "right_anomaly": float(result.right[offset]),
                    "left_Spot_ID": result.left_spots[offset],
                    "right_Spot_ID": result.right_spots[offset],
                }
            )
    return event, rows, summaries


def ratio(numerator, denominator):
    return numerator / denominator if denominator else math.nan


class EventAccumulator:
    """Fixed-size event sums per stratum plus pooled results by dataset type."""

    def __init__(self):
        self.groups = {}

    def add(self, decision):
        for aggregation in ("stratum", "pooled"):
            key = (
                aggregation,
                decision["dataset_type"],
                decision["method"],
                decision["K"] if aggregation == "stratum" else -1,
                decision["detection_efficiency"] if aggregation == "stratum" else -1.0,
                decision["displacement"] if aggregation == "stratum" else -1.0,
            )
            values = self.groups.setdefault(key, Counter())
            values.update(
                {
                    "n": 1,
                    "size": decision["tie_size"],
                    "oracle": int(decision["oracle_resolvable"]),
                    "all_true": int(decision["all_tied_candidates_true"]),
                    "mixed": int(decision["mixed_tie"]),
                    "scoreable": int(decision["method_scoreable"]),
                    "resolved": int(decision["unique_winner"]),
                    "correct": int(
                        decision["unique_winner"]
                        and decision["selected_is_true_corruption"]
                    ),
                }
            )

    def frame(self):
        rows = []
        for key, v in self.groups.items():
            rows.append(
                {
                    **dict(
                        zip(
                            (
                                "aggregation",
                                "dataset_type",
                                "method",
                                "K",
                                "detection_efficiency",
                                "displacement",
                            ),
                            key,
                        )
                    ),
                    "n_ambiguous_events": v["n"],
                    "n_scoreable_events": v["scoreable"],
                    "n_uniquely_resolved_events": v["resolved"],
                    "n_correct_resolutions": v["correct"],
                    "n_incorrect_resolutions": v["resolved"] - v["correct"],
                    "mean_tie_size": v["size"] / v["n"],
                    "oracle_resolvable_fraction": v["oracle"] / v["n"],
                    "fraction_tie_set_containing_true": v["oracle"] / v["n"],
                    "fraction_all_tied_candidates_true": v["all_true"] / v["n"],
                    "fraction_mixed_ties": v["mixed"] / v["n"],
                    "method_coverage": v["scoreable"] / v["n"],
                    "fraction_uniquely_resolved": v["resolved"] / v["n"],
                    "precision_among_resolved": ratio(v["correct"], v["resolved"]),
                    "incorrect_resolution_rate": (v["resolved"] - v["correct"])
                    / v["n"],
                    "effective_correct_resolution_rate": v["correct"] / v["n"],
                    "remaining_unresolved_fraction": 1 - v["resolved"] / v["n"],
                }
            )
        return (
            pd.DataFrame(rows, columns=EVENT_COLUMNS)
            if rows
            else empty_frame(EVENT_COLUMNS)
        )


EVENT_COLUMNS = [
    "aggregation",
    "dataset_type",
    "method",
    "K",
    "detection_efficiency",
    "displacement",
    "n_ambiguous_events",
    "n_scoreable_events",
    "n_uniquely_resolved_events",
    "n_correct_resolutions",
    "n_incorrect_resolutions",
    "mean_tie_size",
    "oracle_resolvable_fraction",
    "fraction_tie_set_containing_true",
    "fraction_all_tied_candidates_true",
    "fraction_mixed_ties",
    "method_coverage",
    "fraction_uniquely_resolved",
    "precision_among_resolved",
    "incorrect_resolution_rate",
    "effective_correct_resolution_rate",
    "remaining_unresolved_fraction",
]


class IterativeAccumulator:
    def __init__(self):
        self.base = C4.ConditionAccumulator()
        self.ties = {}

    def add(self, row):
        self.base.add(row)
        key = tuple(row[k] for k in C4.SUMMARY_GROUP)
        values = self.ties.setdefault(key, Counter())
        values.update(
            {
                key: int(row[key])
                for key in (
                    "tie_invocations",
                    "tie_resolved",
                    "tie_correct",
                    "had_primary_tie",
                )
            }
        )
        values["tie_unresolved"] += int(row["stop_reason"] == "ambiguous_candidate")

    def frame(self):
        frame = self.base.frame()
        for column in (
            "tie_break_invocations_per_trace",
            "fraction_tie_invocations_resolved",
            "tie_break_precision",
            "primary_tie_frequency",
            "fraction_tie_unresolved",
            "mean_false_removals",
            "n_tie_invocations",
            "n_resolved_tie_invocations",
            "n_correct_tie_resolutions",
        ):
            frame[column] = pd.Series(dtype=float)
        for index, row in frame.iterrows():
            v = self.ties[tuple(row[k] for k in C4.SUMMARY_GROUP)]
            frame.loc[index, "n_tie_invocations"] = v["tie_invocations"]
            frame.loc[index, "n_resolved_tie_invocations"] = v["tie_resolved"]
            frame.loc[index, "n_correct_tie_resolutions"] = v["tie_correct"]
            frame.loc[index, "tie_break_invocations_per_trace"] = (
                v["tie_invocations"] / row["n_traces"]
            )
            frame.loc[index, "fraction_tie_invocations_resolved"] = ratio(
                v["tie_resolved"], v["tie_invocations"]
            )
            frame.loc[index, "tie_break_precision"] = ratio(
                v["tie_correct"], v["tie_resolved"]
            )
            frame.loc[index, "primary_tie_frequency"] = (
                v["had_primary_tie"] / row["n_traces"]
            )
            frame.loc[index, "fraction_tie_unresolved"] = (
                v["tie_unresolved"] / row["n_traces"]
            )
            frame.loc[index, "mean_false_removals"] = row[
                "mean_genuine_localizations_removed"
            ]
        metrics = (
            "true_corruption_recall",
            "exact_recovery_rate",
            "mean_false_removals",
            "fraction_with_1plus_false_removal",
            "mean_removals",
        )
        baseline_keys = [key for key in C4.SUMMARY_GROUP if key != "policy"]
        baseline = {
            tuple(row[k] for k in baseline_keys): row
            for _, row in frame.iterrows()
            if row["policy"] == "C_cap_stop"
        }
        for metric in metrics:
            frame[f"delta_{metric}_vs_stop"] = [
                row[metric] - baseline[tuple(row[k] for k in baseline_keys)][metric]
                for _, row in frame.iterrows()
            ]
        return frame


META_FIELDS = [
    "simulation_id",
    "condition",
    "seed",
    "n_barcodes",
    "K",
    "detection_efficiency",
    "displacement",
    "dataset_type",
    "evaluation_seed",
    "reference_seeds",
    "calibration_seeds",
    "nested_reference_seed_sets",
]
TIE_FIELDS = META_FIELDS + [
    "Trace_ID",
    "iteration",
    "n_localizations",
    "C_top3",
    "global_threshold",
    "tie_size",
    "tied_Spot_IDs",
    "tied_Barcodes",
    "n_true_corruptions_remaining",
    "n_true_candidates",
    "tie_class",
    "oracle_resolvable",
    "all_tied_candidates_true",
    "mixed_tie",
    "single_true_candidate_in_tie",
]
SCORE_FIELDS = TIE_FIELDS + [
    "method",
    "method_scoreable",
    "method_scoreability",
    "n_scoreable_candidates",
    "unique_winner",
    "selected_Spot_ID",
    "selected_Barcode",
    "selected_is_true_corruption",
    "candidate_Spot_ID",
    "candidate_Barcode",
    "candidate_is_true_corruption",
    "secondary_score",
    "candidate_scoreable",
    "left_anomaly",
    "right_anomaly",
    "left_Spot_ID",
    "right_Spot_ID",
]
TRACE_FIELDS = META_FIELDS + [
    "method",
    "Trace_ID",
    "policy",
    "n_true_corruptions",
    "n_removed",
    "true_removed",
    "false_removed",
    "true_corruption_recall",
    "removal_precision",
    "exact_recovery",
    "false_clean",
    "overcurated",
    "true_corruptions_remaining",
    "unresolved_true_corruption",
    "abnormal_unresolved",
    "terminal_C",
    "terminal_global_abnormal",
    "terminal_global_state",
    "terminal_outcome",
    "stop_reason",
    "n_iterations",
    "tie_invocations",
    "tie_resolved",
    "tie_correct",
    "had_primary_tie",
]
_STRING_FIELDS = set(
    META_FIELDS[:2]
    + [
        "aggregation",
        "dataset_type",
        "reference_seeds",
        "calibration_seeds",
        "nested_reference_seed_sets",
        "Trace_ID",
        "tied_Spot_IDs",
        "tied_Barcodes",
        "tie_class",
        "method",
        "method_scoreability",
        "selected_Spot_ID",
        "selected_Barcode",
        "candidate_Spot_ID",
        "candidate_Barcode",
        "left_Spot_ID",
        "right_Spot_ID",
        "policy",
        "terminal_global_state",
        "terminal_outcome",
        "stop_reason",
    ]
)


def empty_frame(columns):
    return pd.DataFrame(
        {
            column: pd.Series(dtype=str if column in _STRING_FIELDS else float)
            for column in columns
        }
    )


class DetailWriter(C4.BufferedEcsvWriter):
    """Cycle-4 bounded writer with a readable typed schema even for zero ties."""

    def __init__(self, cycle1, path, columns, max_rows=512):
        super().__init__(cycle1._EcsvChunkWriter(path), max_rows=max_rows)
        self.columns = columns
        self.barcode_numeric = False

    def write_rows(self, rows):
        if rows:
            self.write(pd.DataFrame(rows, columns=self.columns))

    def finish(self):
        if not self.pending and not self.writer._has_rows:
            C4._write(empty_frame(self.columns), self.writer.path)
        else:
            super().finish()


def add_coordinates(frame, assume_uniform):
    if "Genomic_Position" not in frame and assume_uniform:
        mapping = {b: i + 1 for i, b in enumerate(sorted(frame["Barcode"].unique()))}
        frame["Genomic_Position"] = frame["Barcode"].map(mapping)
    return frame


def select_conditions(manifest, cycle3=None):
    cycle3 = cycle3 or C4._cycle3()
    selected = []
    classes = set()
    for item in manifest["conditions"]:
        n, k = cycle3._condition_metadata(item)
        if n != 25 or k not in (1, 2, 3):
            continue
        p = next(
            (
                value
                for value in (0.5, 0.8)
                if np.isclose(float(item["detection_efficiency"]), value)
            ),
            None,
        )
        d = next(
            (
                value
                for value in (0.4, 0.8)
                if np.isclose(float(item["displacement_um"]), value)
            ),
            None,
        )
        if p is not None and d is not None:
            selected.append(item)
            classes.add((k, p, d))
    required = {(k, p, d) for k in (1, 2, 3) for p in (0.5, 0.8) for d in (0.4, 0.8)}
    if classes != required:
        raise ValueError(
            f"Missing cycle-5 primary/K1 condition classes: {sorted(required - classes)}"
        )
    return selected


def load_group_baseline(
    cycle3, cycle1, root, manifest, n_barcodes, efficiency, assume_uniform
):
    sizes = {
        str(item["simulation_id"]): cycle3._condition_metadata(item)[0]
        for item in manifest["conditions"]
    }
    parts = []
    for item in manifest["simulations"]:
        if not np.isclose(float(item["detection_efficiency"]), efficiency):
            continue
        known = item.get("n_barcodes", sizes.get(str(item["id"])))
        if known is not None and int(known) != n_barcodes:
            continue
        frame = cycle1.load_baseline(root / "simulations" / str(item["id"]), item)
        n = int(
            item.get("n_barcodes", sizes.get(str(item["id"]), frame["Barcode"].max()))
        )
        if n == n_barcodes:
            frame["n_barcodes"] = n
            parts.append(add_coordinates(frame, assume_uniform))
    if not parts:
        raise ValueError(
            f"No matching baseline population for B={n_barcodes}, p={efficiency}"
        )
    return pd.concat(parts, ignore_index=True)


class ReferenceCache:
    """Owned and cleared by exactly one (B, detection_efficiency) group."""

    def __init__(self, cycle3, baseline, minimum, mode):
        self.cycle3, self.baseline, self.minimum, self.mode = (
            cycle3,
            baseline,
            minimum,
            mode,
        )
        self.models = {}
        self.hits = self.misses = 0

    def fit(self, seeds):
        key = tuple(sorted(map(int, seeds)))
        if key not in self.models:
            self.misses += 1
            self.models[key] = C4._reference_with_coordinates(
                self.cycle3,
                self.baseline[self.baseline["seed"].isin(seeds)],
                self.minimum,
                self.mode,
            )
        else:
            self.hits += 1
        return self.models[key]

    def clear(self):
        self.models.clear()
        self.baseline = None


def calibrate_fold(
    baseline, seeds, evaluation_seed, references, mode, assume_uniform, progress
):
    compact, provenance = [], []
    for calibration_seed in seeds:
        if calibration_seed == evaluation_seed:
            continue
        reference_seeds, _, _ = C4.fold_seed_sets(
            seeds, evaluation_seed, calibration_seed
        )
        provenance.append(
            f"{calibration_seed}:{','.join(map(str, sorted(reference_seeds)))}"
        )
        nested = references.fit(reference_seeds)
        progress(f"calibration={calibration_seed} start")
        for index, (_, trace) in enumerate(
            baseline[baseline["seed"] == calibration_seed].groupby(
                ["simulation_id", "Trace_ID"], sort=False
            ),
            start=1,
        ):
            compact.append(
                C4.CompactTraceStats.from_scores(
                    C4.build_trace_scores(trace, nested, mode, assume_uniform)
                )
            )
            del trace
            if index % 250 == 0:
                progress(f"calibration={calibration_seed} traces={index}")
        del nested
    if not compact:
        raise ValueError(f"No calibration traces for evaluation seed {evaluation_seed}")
    thresholds = C4.calibrate_thresholds(compact)
    primary = thresholds[np.isclose(thresholds["trace_fpr"], 0.01)].iloc[0].to_dict()
    if not np.isfinite(primary["global_threshold"]):
        raise ValueError(
            "Clean calibration produced an unscoreable primary global threshold"
        )
    return primary, ";".join(provenance)


def analyze_trace(
    trace,
    reference,
    thresholds,
    meta,
    budget,
    stage="both",
    mode="separation",
    assume_uniform=False,
    builder=None,
):
    """Atomic per-trace work: a profile cutoff never emits partial policy sets."""
    scorer = TraceScorer(builder)
    session = TieSession(scorer, reference, mode, assume_uniform, budget)
    try:
        audit, baseline = C4.run_policy(
            trace,
            reference,
            thresholds,
            "C_cap",
            max_removals=3,
            min_remaining=3,
            reference_mode=mode,
            allow_uniform_fallback=assume_uniform,
            scorer=scorer,
        )
        terminal = scorer.last_state
        events, scores, decisions = [], [], []
        baseline_tie = baseline["stop_reason"] == "ambiguous_candidate"
        if baseline_tie:
            event, scores, decisions = tie_rows(
                terminal, baseline["n_iterations"], thresholds, meta, session
            )
            events.append(event)
        metrics = []
        if stage == "both":
            counts = {
                "tie_invocations": int(baseline_tie),
                "tie_resolved": 0,
                "tie_correct": 0,
                "had_primary_tie": baseline_tie,
            }
            metrics.append(
                {
                    **meta,
                    **baseline,
                    **counts,
                    "method": METHODS[0],
                    "policy": POLICIES[METHODS[0]],
                }
            )
            for method in METHODS[1:]:
                result, counts = run_integrated(
                    trace,
                    reference,
                    thresholds,
                    method,
                    scorer=scorer,
                    resolver=session.resolve,
                    reference_mode=mode,
                    assume_uniform=assume_uniform,
                )
                metrics.append(
                    {
                        **meta,
                        **result,
                        **counts,
                        "method": method,
                        "policy": POLICIES[method],
                    }
                )
        return (
            events,
            scores,
            decisions,
            metrics,
            {
                "score_builds": scorer.calls,
                "score_cache_hits": scorer.hits,
                "peak_trace_states": scorer.peak_states,
            },
        )
    finally:
        scorer.clear()


class Monitor:
    def __init__(self):
        self.started = time.monotonic()
        self.peak = 0.0
        self.group_peak = 0.0
        self.traces = 0
        self.events = 0
        self.context = ""

    def rss(self):
        value = C4._rss_mb()
        self.peak = max(self.peak, value)
        self.group_peak = max(self.group_peak, value)
        return value

    def progress(self, message):
        print(
            f"[{time.monotonic() - self.started:8.1f}s] RSS={self.rss():.1f} MiB "
            f"{self.context} traces={self.traces} ambiguous_events={self.events} {message}",
            file=sys.stderr,
            flush=True,
        )


def benchmark(
    root,
    output,
    *,
    assume_uniform=False,
    reference_mode="separation",
    minimum_observations=20,
    profile_ambiguous_events=0,
    stage="both",
    cycle4_results=None,
):
    manifest_path = root / "sweep_manifest.yaml"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Existing cycle-3 benchmark is required: {manifest_path}. No benchmark is generated."
        )
    manifest = yaml.safe_load(manifest_path.read_text())
    if (
        manifest.get("dry_run")
        or not manifest.get("simulations")
        or not manifest.get("conditions")
    ):
        raise ValueError(
            "A completed cycle-3 sweep manifest with simulations and conditions is required"
        )
    cycle3 = C4._cycle3()
    selected = select_conditions(manifest, cycle3)
    if cycle4_results is not None and not cycle4_results.is_dir():
        raise FileNotFoundError(
            f"Cycle-4 provenance directory not found: {cycle4_results}"
        )
    cycle1 = cycle3._cycle(1)
    output.mkdir(parents=True, exist_ok=True)
    writers = (
        DetailWriter(cycle1, output / OUTPUTS[0], TIE_FIELDS),
        DetailWriter(cycle1, output / OUTPUTS[1], SCORE_FIELDS),
        DetailWriter(cycle1, output / "cycle5_trace_metrics.ecsv", TRACE_FIELDS),
    )
    thresholds_writer = C4.BufferedEcsvWriter(
        cycle1._EcsvChunkWriter(output / "cycle5_thresholds.ecsv")
    )
    event_aggregate, iterative_aggregate = EventAccumulator(), IterativeAccumulator()
    monitor, budget = Monitor(), EventBudget(profile_ambiguous_events)
    runtime, totals = [], Counter()
    rss_start = monitor.rss()
    limited = False
    for efficiency in (0.5, 0.8):
        group_started = time.monotonic()
        monitor.group_peak = 0.0
        monitor.context = f"group B=25 p={efficiency:g}"
        group_start_rss = monitor.rss()
        monitor.progress("start group")
        baseline = load_group_baseline(
            cycle3, cycle1, root, manifest, 25, efficiency, assume_uniform
        )
        # Seed membership comes from the full C4 group, not the narrow subset.
        group_items = [
            item
            for item in manifest["conditions"]
            if cycle3._condition_metadata(item)[0] == 25
            and np.isclose(float(item["detection_efficiency"]), efficiency)
        ]
        seeds = sorted({int(item["seed"]) for item in group_items})
        if len(seeds) < 3:
            raise ValueError(
                "At least three replicate seeds are required for nested calibration"
            )
        references = ReferenceCache(
            cycle3, baseline, minimum_observations, reference_mode
        )
        group_traces_start, group_events_start = monitor.traces, monitor.events
        try:
            for evaluation_seed in seeds:
                fold_context = f"group B=25 p={efficiency:g} eval={evaluation_seed}"
                monitor.context = fold_context
                monitor.progress("start evaluation seed")
                reference_seeds, _, _ = C4.fold_seed_sets(seeds, evaluation_seed)
                reference = references.fit(reference_seeds)
                thresholds, nested_sets = calibrate_fold(
                    baseline,
                    seeds,
                    evaluation_seed,
                    references,
                    reference_mode,
                    assume_uniform,
                    monitor.progress,
                )
                fold_meta = {
                    "evaluation_seed": evaluation_seed,
                    "reference_seeds": ",".join(map(str, sorted(reference_seeds))),
                    "calibration_seeds": ",".join(
                        map(str, [s for s in seeds if s != evaluation_seed])
                    ),
                    "nested_reference_seed_sets": nested_sets,
                }
                thresholds_writer.write(
                    pd.DataFrame(
                        [
                            {
                                "n_barcodes": 25,
                                "detection_efficiency": efficiency,
                                **fold_meta,
                                **thresholds,
                            }
                        ]
                    )
                )

                def process(frame, meta):
                    for identity, trace in frame.groupby(
                        ["simulation_id", "Trace_ID"], sort=False
                    ):
                        trace_meta = {
                            **meta,
                            **fold_meta,
                            "simulation_id": str(identity[0]),
                        }

                        def builder(*args):
                            state = C4.build_trace_scores(*args)
                            monitor.rss()
                            return state

                        try:
                            events, scores, decisions, metrics, diagnostics = (
                                analyze_trace(
                                    trace,
                                    reference,
                                    thresholds,
                                    trace_meta,
                                    budget,
                                    stage,
                                    reference_mode,
                                    assume_uniform,
                                    builder,
                                )
                            )
                        except ProfileLimit:
                            monitor.progress(
                                "profile budget reached; incomplete trace excluded"
                            )
                            return False
                        writers[0].write_rows(events)
                        writers[1].write_rows(scores)
                        writers[2].write_rows(metrics)
                        for decision in decisions:
                            event_aggregate.add(decision)
                        for row in metrics:
                            iterative_aggregate.add(row)
                        totals.update(
                            {
                                key: value
                                for key, value in diagnostics.items()
                                if key != "peak_trace_states"
                            }
                        )
                        totals["peak_trace_states"] = max(
                            totals["peak_trace_states"],
                            diagnostics["peak_trace_states"],
                        )
                        monitor.traces += 1
                        monitor.events += len(events)
                        del trace
                        if monitor.traces % 250 == 0 or (
                            events and monitor.events % 100 == 0
                        ):
                            monitor.progress("progress")
                    return True

                monitor.context = f"{fold_context} condition=held-out-clean"
                clean_meta = {
                    "condition": f"clean-evaluation-25-{efficiency:g}-{evaluation_seed}",
                    "seed": evaluation_seed,
                    "n_barcodes": 25,
                    "K": 0,
                    "detection_efficiency": efficiency,
                    "displacement": 0.0,
                    "dataset_type": "clean_control",
                }
                monitor.progress("start clean control")
                if not process(
                    baseline[baseline["seed"] == evaluation_seed], clean_meta
                ):
                    limited = True
                    break
                for item in selected:
                    if int(item["seed"]) != evaluation_seed or not np.isclose(
                        float(item["detection_efficiency"]), efficiency
                    ):
                        continue
                    monitor.context = f"{fold_context} condition={item['id']}"
                    monitor.progress("start condition")
                    observations = add_coordinates(
                        cycle3._load_observations(cycle1, root, item), assume_uniform
                    )
                    meta = {
                        **cycle3._trace_metadata(item),
                        "dataset_type": (
                            "k1_control"
                            if cycle3._condition_metadata(item)[1] == 1
                            else "primary"
                        ),
                    }
                    complete = process(observations, meta)
                    del observations
                    if not complete:
                        limited = True
                        break
                    monitor.progress("end condition")
                del reference
                monitor.context = fold_context
                monitor.progress("end evaluation seed")
                if limited:
                    break
        finally:
            before_clear = len(references.models)
            references.clear()
            # The process closure captures this cell; releasing it also releases
            # the outer model after a profiling exception in clean controls.
            reference = None
            del baseline
            runtime.append(
                {
                    "scope": "group",
                    "n_barcodes": 25,
                    "detection_efficiency": efficiency,
                    "elapsed_seconds": time.monotonic() - group_started,
                    "n_traces": monitor.traces - group_traces_start,
                    "n_stage1_ambiguous_events": monitor.events - group_events_start,
                    "rss_start_mb": group_start_rss,
                    "rss_end_mb": monitor.rss(),
                    "rss_peak_mb": monitor.group_peak,
                    "reference_cache_hits": references.hits,
                    "reference_cache_misses": references.misses,
                    "reference_cache_before_clear": before_clear,
                    "reference_cache_after_clear": len(references.models),
                }
            )
            monitor.progress("end group; reference cache cleared")
        if limited:
            break
    for writer in writers:
        writer.finish()
    thresholds_writer.finish()
    event_frame, iterative_frame = event_aggregate.frame(), iterative_aggregate.frame()
    C4._write(event_frame, output / OUTPUTS[2])
    for dataset, name in (
        ("primary", OUTPUTS[3]),
        ("clean_control", OUTPUTS[4]),
        ("k1_control", OUTPUTS[5]),
    ):
        C4._write(
            iterative_frame[iterative_frame["dataset_type"] == dataset], output / name
        )
    plot_results(output, event_frame, iterative_frame)
    runtime.append(
        {
            "scope": "total",
            "elapsed_seconds": time.monotonic() - monitor.started,
            "n_traces": monitor.traces,
            "n_stage1_ambiguous_events": monitor.events,
            "n_unique_tie_states_calculated": budget.used,
            "profile_ambiguous_events": profile_ambiguous_events,
            "profile_limit_reached": limited,
            "stage": stage,
            "reference_mode": reference_mode,
            "assume_uniform_barcode_spacing": assume_uniform,
            "minimum_reference_observations": minimum_observations,
            "benchmark_root": str(root),
            "cycle4_results_provenance": str(cycle4_results or ""),
            "cycle4_results_role": "provenance_only",
            "rss_start_mb": rss_start,
            "rss_end_mb": monitor.rss(),
            "rss_peak_mb": monitor.peak,
            **totals,
            "reference_cache_hits": sum(row["reference_cache_hits"] for row in runtime),
            "reference_cache_misses": sum(
                row["reference_cache_misses"] for row in runtime
            ),
        }
    )
    C4._write(pd.DataFrame(runtime), output / OUTPUTS[6])
    write_report(output, event_frame, iterative_frame, runtime[-1])
    monitor.progress("complete; conservative default remains stop_unresolved")


def plot_results(output, events, iterative):
    import matplotlib.pyplot as plt

    specs = [
        (
            events[
                (events["dataset_type"] == "primary")
                & (events["aggregation"] == "stratum")
            ],
            "method",
            "precision_among_resolved",
            "tie_resolution_precision.png",
        ),
        (
            events[
                (events["dataset_type"] == "primary")
                & (events["aggregation"] == "stratum")
            ],
            "method",
            "fraction_uniquely_resolved",
            "tie_resolution_coverage.png",
        ),
        (
            events[
                (events["dataset_type"] == "primary")
                & (events["aggregation"] == "stratum")
            ],
            "method",
            "effective_correct_resolution_rate",
            "effective_correct_resolution.png",
        ),
    ]
    for metric, name in (
        ("exact_recovery_rate", "iterative_exact_recovery.png"),
        ("true_corruption_recall", "iterative_recall.png"),
        ("mean_false_removals", "iterative_false_removals.png"),
        ("fraction_abnormal_unresolved", "ambiguous_unresolved_fraction.png"),
    ):
        specs.append(
            (iterative[iterative["dataset_type"] == "primary"], "policy", metric, name)
        )
    for frame, label, metric, name in specs:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
        for ax, k in zip(axes, (2, 3)):
            subset = frame[frame["K"] == k]
            for method, group in subset.groupby(label, sort=False):
                ordered = group.sort_values(["detection_efficiency", "displacement"])
                x = [
                    f"p={p:g}, d={d:g}"
                    for p, d in zip(
                        ordered["detection_efficiency"], ordered["displacement"]
                    )
                ]
                ax.plot(x, ordered[metric], marker="o", label=method)
            ax.set(title=f"K={k}", ylabel=metric)
            ax.tick_params(axis="x", rotation=25)
            if len(subset):
                ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(output / name, dpi=160)
        plt.close(fig)


def write_report(output, events, iterative, runtime):
    def number(value):
        return f"{value:.4f}" if np.isfinite(value) else "unavailable"

    pooled = events[
        (events["aggregation"] == "pooled") & (events["dataset_type"] == "primary")
    ]
    primary_events = int(pooled.iloc[0]["n_ambiguous_events"]) if len(pooled) else 0
    oracle = pooled.iloc[0]["oracle_resolvable_fraction"] if len(pooled) else math.nan
    lines = [
        "Cycle-5 targeted tie-resolution validation",
        f"Benchmark: {runtime['benchmark_root']}",
        f"Stage: {runtime['stage']}; profile-limited: {runtime['profile_limit_reached']}",
        f"Stage-1 ambiguous events analyzed (including controls): {runtime['n_stage1_ambiguous_events']}",
        f"Primary ambiguous events analyzed: {primary_events}; oracle-resolvable fraction: {number(oracle)}",
        f"Unique tie states calculated (including integrated paths): {runtime['n_unique_tie_states_calculated']}",
    ]
    for method in METHODS:
        selected = pooled[pooled["method"] == method]
        if not len(selected):
            lines.append(
                f"{method}: precision and coverage unavailable (no primary tie events)"
            )
            continue
        row = selected.iloc[0]
        lines.append(
            f"{method}: events={row['n_ambiguous_events']}, resolved={row['n_uniquely_resolved_events']}, "
            f"precision={number(row['precision_among_resolved'])}, score coverage={number(row['method_coverage'])}, "
            f"unique resolution={number(row['fraction_uniquely_resolved'])}, "
            f"effective correct={number(row['effective_correct_resolution_rate'])}"
        )
    lines.append("Integrated results by stratum (differences relative to C_cap_stop):")
    if not len(iterative):
        lines.append("Integrated validation was not run or has no completed traces.")
    for _, row in iterative.iterrows():
        lines.append(
            f"{row['dataset_type']} {row['policy']} K={row['K']} p={row['detection_efficiency']:g} d={row['displacement']:g}: "
            f"recall={row['true_corruption_recall']:.4f} (delta={row['delta_true_corruption_recall_vs_stop']:+.4f}), "
            f"exact={row['exact_recovery_rate']:.4f} (delta={row['delta_exact_recovery_rate_vs_stop']:+.4f}), "
            f"mean false={row['mean_false_removals']:.4f} (delta={row['delta_mean_false_removals_vs_stop']:+.4f}), "
            f"any removal={row['fraction_with_1plus_removal']:.4f}, >=2={row['fraction_with_2plus_removals']:.4f}, "
            f"max removals={row['maximum_removals']}, terminal abnormal={row['fraction_terminal_global_abnormal']:.4f}, "
            f"tie frequency={number(row['primary_tie_frequency'])}, tie-break precision={number(row['tie_break_precision'])}"
        )
    lines += [
        f"Elapsed seconds: {runtime['elapsed_seconds']:.3f}; sampled peak RSS MiB: {runtime['rss_peak_mb']:.1f}",
        "No automated production recommendation. Stop unresolved remains the conservative default.",
        "Inspect precision, recovery, false removals, clean controls, and every p/d stratum together.",
    ]
    (output / "cycle5_report.txt").write_text("\n".join(lines) + "\n")


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-root", type=Path, default=Path("benchmark_curator_loo_validation")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--assume-uniform-barcode-spacing", action="store_true")
    parser.add_argument(
        "--reference-mode", choices=("separation", "barcode_pair"), default="separation"
    )
    parser.add_argument("--minimum-reference-observations", type=int, default=20)
    parser.add_argument(
        "--cycle4-results",
        type=Path,
        help="provenance only; all scoring/thresholds are reconstructed from original data",
    )
    parser.add_argument(
        "--profile-ambiguous-events",
        type=int,
        default=0,
        help="cap unique tie states across offline/integrated paths; exclude incomplete traces",
    )
    parser.add_argument(
        "--stage",
        choices=("offline", "both"),
        default="both",
        help="offline stage alone, or offline plus integrated validation (default)",
    )
    args = parser.parse_args(argv)
    if args.profile_ambiguous_events < 0 or args.minimum_reference_observations < 1:
        parser.error(
            "profile limit must be nonnegative and minimum observations must be positive"
        )
    return args


def main():
    args = parse_arguments()
    try:
        benchmark(
            args.benchmark_root.resolve(),
            args.output_dir.resolve(),
            assume_uniform=args.assume_uniform_barcode_spacing,
            reference_mode=args.reference_mode,
            minimum_observations=args.minimum_reference_observations,
            profile_ambiguous_events=args.profile_ambiguous_events,
            stage=args.stage,
            cycle4_results=args.cycle4_results,
        )
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
