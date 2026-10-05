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
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table

FIXED_FPRS = (0.001, 0.01, 0.05)
POLICIES = ("A_global", "B_absolute", "B_relative", "C_cap")
OUTPUTS = (
    "cycle4_trace_metrics.ecsv",
    "cycle4_condition_metrics.ecsv",
    "cycle4_iteration_audit.ecsv",
    "cycle4_thresholds.ecsv",
    "cycle4_runtime_profile.ecsv",
)


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
    sample = np.sort(np.asarray(reference, float))
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
                else tuple(sorted((int(barcodes[left]), int(barcodes[right]))))
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


def calibrate_thresholds(clean_states: Sequence[TraceScores]) -> pd.DataFrame:
    """Calibrate global and secondary safeguards from clean traces only."""
    samples = {"global": [], "absolute": [], "relative": []}
    for state in clean_states:
        finite = np.flatnonzero(np.isfinite(state.loo))
        samples["global"].append(state.cost)
        if len(finite):
            best = finite[np.argmax(state.loo[finite])]
            samples["absolute"].append(state.loo[best])
            samples["relative"].append(state.loo[best] / max(state.cost, 1e-12))
    rows = []
    for fpr in FIXED_FPRS:
        global_threshold = empirical_threshold(samples["global"], fpr)
        rows.append(
            {
                "trace_fpr": fpr,
                "global_threshold": global_threshold,
                "absolute_threshold": empirical_threshold(samples["absolute"], fpr),
                "relative_threshold": empirical_threshold(samples["relative"], fpr),
                "observed_clean_global_fpr": float(
                    np.mean(np.asarray(samples["global"]) > global_threshold)
                ),
                "n_clean_traces": len(clean_states),
            }
        )
    return pd.DataFrame(rows)


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
        candidate = int(finite[np.argmax(state.loo[finite])]) if len(finite) else None
        if not stop and candidate is None:
            stop = "no_scoreable_candidate"
        absolute = float(state.loo[candidate]) if candidate is not None else math.nan
        after = (
            float(state.cost_without[candidate]) if candidate is not None else math.nan
        )
        relative = (
            absolute / max(state.cost, 1e-12) if candidate is not None else math.nan
        )
        actionable = candidate is not None
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
                "candidate_actionable": actionable,
                "C_after_candidate_removal": after,
                "absolute_C_drop": absolute,
                "relative_C_drop": relative,
                "removal_accepted": accepted,
                "candidate_is_corrupted": bool(row.get("is_corrupted", False))
                if row is not None
                else False,
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
    return pd.DataFrame(audit), {
        "Trace_ID": str(trace["Trace_ID"].iloc[0]),
        "policy": policy,
        "n_true_corruptions": n_true,
        "n_removed": len(removed),
        "true_removed": true_removed,
        "false_removed": false_removed,
        "true_corruption_recall": true_removed / n_true if n_true else math.nan,
        "removal_precision": true_removed / len(removed) if removed else math.nan,
        "exact_recovery": len(removed) == n_true and true_removed == n_true,
        "true_corruptions_remaining": n_true - true_removed,
        "unresolved_true_corruption": n_true > true_removed,
        "abnormal_unresolved": reason == "candidate_not_actionable",
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
    separation, pairs = {}, {}
    sep_values, pair_values = {}, {}
    for _, trace in baseline.groupby(["simulation_id", "Trace_ID"], sort=False):
        positions = _positions(trace, allow_uniform_fallback=False)
        xyz = trace[["x", "y", "z"]].to_numpy(float)
        barcode = trace["Barcode"].to_numpy()
        for left in range(len(trace)):
            for right in range(left + 1, len(trace)):
                distance = float(np.linalg.norm(xyz[right] - xyz[left]))
                sep_values.setdefault(
                    float(abs(positions[right] - positions[left])), []
                ).append(distance)
                pair_values.setdefault(
                    tuple(sorted((int(barcode[left]), int(barcode[right])))), []
                ).append(distance)
    separation = {key: np.asarray(value) for key, value in sep_values.items()}
    pairs = {key: np.asarray(value) for key, value in pair_values.items()}
    return SimpleNamespace(
        separation=separation,
        barcode_pair=pairs,
        minimum_observations=minimum,
    )


def benchmark(
    root: Path,
    output: Path,
    minimum_observations: int = 20,
    reference_mode: str = "separation",
    assume_uniform: bool = False,
    max_removals: int = 3,
    min_remaining: int = 3,
    profile_traces: int = 0,
) -> None:
    """Run cycle 4 against the existing cycle-3 manifest and datasets."""
    started = time.monotonic()
    cycle3 = _cycle3()
    cycle1 = cycle3._cycle(1)
    manifest = yaml.safe_load((root / "sweep_manifest.yaml").read_text())
    conditions = manifest.get("conditions", [])
    if manifest.get("dry_run") or not conditions:
        raise ValueError("A completed cycle-3 sweep manifest is required")
    baselines = {}
    for item in manifest["simulations"]:
        frame = cycle1.load_baseline(root / "simulations" / str(item["id"]), item)
        if "Genomic_Position" not in frame and assume_uniform:
            mapping = {
                b: i + 1 for i, b in enumerate(sorted(frame["Barcode"].unique()))
            }
            frame["Genomic_Position"] = frame["Barcode"].map(mapping)
        baselines[str(item["id"])] = frame
    output.mkdir(parents=True, exist_ok=True)
    traces_out, audits_out, threshold_out = [], [], []
    groups = {}
    for item in conditions:
        n, _ = cycle3._condition_metadata(item)
        groups.setdefault((n, float(item["detection_efficiency"])), []).append(item)
    processed = 0
    for (n_barcodes, efficiency), items in groups.items():
        seeds = sorted({int(item["seed"]) for item in items})
        baseline = pd.concat(
            [
                f
                for f in baselines.values()
                if int(f.get("n_barcodes", pd.Series([n_barcodes])).iloc[0])
                == n_barcodes
                and np.isclose(f["detection_efficiency"].iloc[0], efficiency)
            ]
        )
        for evaluation_seed in seeds:
            ref_seeds, _, _ = fold_seed_sets(seeds, evaluation_seed)
            reference = _reference_with_coordinates(
                cycle3,
                baseline[baseline["seed"].isin(ref_seeds)],
                minimum_observations,
                reference_mode,
            )
            clean_states = []
            for calibration_seed in seeds:
                if calibration_seed == evaluation_seed:
                    continue
                nested, _, _ = fold_seed_sets(seeds, evaluation_seed, calibration_seed)
                nested_ref = _reference_with_coordinates(
                    cycle3,
                    baseline[baseline["seed"].isin(nested)],
                    minimum_observations,
                    reference_mode,
                )
                for _, clean in baseline[baseline["seed"] == calibration_seed].groupby(
                    ["simulation_id", "Trace_ID"], sort=False
                ):
                    clean_states.append(
                        build_trace_scores(
                            clean, nested_ref, reference_mode, assume_uniform
                        )
                    )
            thresholds = calibrate_thresholds(clean_states)
            for threshold in thresholds.to_dict("records"):
                threshold_out.append(
                    {
                        "n_barcodes": n_barcodes,
                        "detection_efficiency": efficiency,
                        "evaluation_seed": evaluation_seed,
                        "reference_seeds": ",".join(map(str, sorted(ref_seeds))),
                        "evaluation_seeds": str(evaluation_seed),
                        **threshold,
                    }
                )
            primary = (
                thresholds[np.isclose(thresholds["trace_fpr"], 0.01)].iloc[0].to_dict()
            )
            for item in items:
                if int(item["seed"]) != evaluation_seed:
                    continue
                observations = cycle3._load_observations(cycle1, root, item)
                if "Genomic_Position" not in observations and assume_uniform:
                    mapping = {
                        b: i + 1
                        for i, b in enumerate(sorted(observations["Barcode"].unique()))
                    }
                    observations["Genomic_Position"] = observations["Barcode"].map(
                        mapping
                    )
                for identity, trace in observations.groupby(
                    ["simulation_id", "Trace_ID"], sort=False
                ):
                    if profile_traces and processed >= profile_traces:
                        break
                    meta = cycle3._trace_metadata(item)
                    score_cache = {}

                    def shared_scorer(current, model, mode, fallback):
                        key = tuple(sorted(current["Spot_ID"].astype(str)))
                        if key not in score_cache:
                            score_cache[key] = build_trace_scores(
                                current, model, mode, fallback
                            )
                        return score_cache[key]

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
                        audit = audit.assign(simulation_id=identity[0], **meta)
                        audits_out.append(audit)
                        traces_out.append(
                            {"simulation_id": identity[0], **meta, **metrics}
                        )
                    processed += 1
                if profile_traces and processed >= profile_traces:
                    break
            if profile_traces and processed >= profile_traces:
                break
        if profile_traces and processed >= profile_traces:
            break
    trace_metrics = pd.DataFrame(traces_out)
    audit = pd.concat(audits_out, ignore_index=True) if audits_out else pd.DataFrame()
    condition = summarize(
        trace_metrics,
        ["policy", "n_barcodes", "K", "detection_efficiency", "displacement"],
    )
    _write(trace_metrics, output / OUTPUTS[0])
    _write(condition, output / OUTPUTS[1])
    _write(audit, output / OUTPUTS[2])
    _write(pd.DataFrame(threshold_out), output / OUTPUTS[3])
    elapsed = time.monotonic() - started
    _write(
        pd.DataFrame(
            [
                {
                    "profile_traces": processed,
                    "elapsed_seconds": elapsed,
                    "seconds_per_trace_all_policies": elapsed / max(processed, 1),
                    "projected_hours_for_100000_traces": elapsed
                    / max(processed, 1)
                    * 100000
                    / 3600,
                }
            ]
        ),
        output / OUTPUTS[4],
    )
    plot_results(output)


def plot_results(output: Path) -> None:
    condition = Table.read(output / OUTPUTS[1], format="ascii.ecsv").to_pandas()
    audit = Table.read(output / OUTPUTS[2], format="ascii.ecsv").to_pandas()
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
        for policy, frame in condition.groupby("policy", sort=False):
            line = frame.groupby("displacement")[metric].mean().sort_index()
            ax.plot(line.index, line.values, marker="o", label=policy)
        ax.set(xlabel="Displacement (µm)", ylabel=metric)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(output / name, dpi=160)
        plt.close(fig)
    clean = condition[condition["K"] == 0]
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
    for policy, frame in audit.groupby("policy", sort=False):
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
    )


if __name__ == "__main__":
    main()
