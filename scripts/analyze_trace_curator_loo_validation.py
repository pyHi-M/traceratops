#!/usr/bin/env python3
"""Final, focused validation of the empirical ``loo_top3_all`` model.

This program is benchmark code, not the production trace curator.  It performs
replicate-held-out reference fitting and nested clean-trace calibration while
processing one condition at a time.
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table

FIXED_FPRS = (0.001, 0.01, 0.05)
EPSILON = 1e-12
OUTPUTS = (
    "loo_single_anomaly_performance.ecsv",
    "loo_multi_anomaly_single_pass.ecsv",
    "loo_iterative_recovery.ecsv",
    "loo_iterative_thresholded_performance.ecsv",
    "loo_confidence_metrics.ecsv",
    "loo_confidence_precision_coverage.ecsv",
    "loo_clean_trace_thresholds.ecsv",
    "loo_reference_summary.ecsv",
    "loo_clean_trace_evaluation.ecsv",
    "loo_confidence_distributions.ecsv",
)
_START = time.monotonic()


def _cycle(number: int):
    name = (
        "analyze_trace_curator_benchmark.py"
        if number == 1
        else "analyze_trace_curator_attribution_benchmark.py"
    )
    spec = importlib.util.spec_from_file_location(
        f"_trace_curator_cycle{number}", Path(__file__).with_name(name)
    )
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


def _progress(**values) -> None:
    detail = "; ".join(f"{key}={value}" for key, value in values.items())
    print(
        f"[loo validation {time.monotonic() - _START:8.1f}s] {detail}",
        file=sys.stderr,
        flush=True,
    )


def fold_seed_sets(
    seeds: Sequence[int], evaluation_seed: int, calibration_seed: int | None = None
):
    evaluation = {int(evaluation_seed)}
    calibration = set() if calibration_seed is None else {int(calibration_seed)}
    reference = set(map(int, seeds)) - evaluation - calibration
    if reference & calibration or reference & evaluation or calibration & evaluation:
        raise AssertionError("Reference/calibration/evaluation replicate leakage")
    return reference, calibration, evaluation


def confidence_metrics(scores: Sequence[float], epsilon: float = EPSILON) -> dict:
    """Return stable winner metrics.  Ratios with a near-zero runner-up are finite."""
    values = np.asarray(scores, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "S1": math.nan,
            "S2": math.nan,
            "absolute_margin": math.nan,
            "relative_margin": math.nan,
            "score_ratio": math.nan,
            "winner_is_unique": False,
        }
    ordered = np.sort(values)[::-1]
    first = float(ordered[0])
    second = float(ordered[1]) if len(ordered) > 1 else 0.0
    margin = first - second
    return {
        "S1": first,
        "S2": second,
        "absolute_margin": margin,
        "relative_margin": margin / max(abs(first), epsilon),
        "score_ratio": first / max(second, epsilon),
        "winner_is_unique": bool(len(ordered) == 1 or first > second),
    }


def trace_threshold(values: Sequence[float], fpr: float) -> float:
    """Largest-null empirical cutoff permitting at most ``fpr`` exceedances."""
    finite = np.sort(np.asarray(values, float)[np.isfinite(values)])
    if not len(finite):
        return math.nan
    allowed = int(math.floor(float(fpr) * len(finite)))
    return float(finite[max(0, len(finite) - allowed - 1)])


def _top3(values: Sequence[float]) -> float:
    finite = np.asarray(values, float)
    finite = finite[np.isfinite(finite)]
    return (
        float(np.sort(finite)[-min(3, len(finite)) :].sum())
        if len(finite)
        else math.nan
    )


def empirical_anomaly(reference: Sequence[float], observed: float) -> float:
    """Finite-sample corrected, two-sided empirical tail score."""
    sample = np.sort(np.asarray(reference, float))
    lower = (np.searchsorted(sample, observed, side="right") + 1) / (len(sample) + 1)
    upper = (len(sample) - np.searchsorted(sample, observed, side="left") + 1) / (
        len(sample) + 1
    )
    return -math.log10(min(1.0, 2.0 * min(lower, upper)))


def score_trace(
    trace: pd.DataFrame, reference, minimum_localizations: int = 3
) -> tuple[pd.DataFrame, dict]:
    """Calculate exactly C(T)-C(T\\i) for every localization in a trace."""
    trace = trace.sort_values(["Barcode", "Spot_ID"]).reset_index(drop=True)
    xyz = trace[["x", "y", "z"]].to_numpy(float)
    barcodes = trace["Barcode"].to_numpy(float)
    edges: list[tuple[int, int, float]] = []
    for left in range(len(trace)):
        for right in range(left + 1, len(trace)):
            sample = reference.separation.get(
                float(abs(barcodes[right] - barcodes[left]))
            )
            if sample is None or len(sample) < reference.minimum_observations:
                continue
            distance = float(np.linalg.norm(xyz[right] - xyz[left]))
            edges.append((left, right, empirical_anomaly(sample, distance)))
    before = _top3([edge[2] for edge in edges])
    scores = np.full(len(trace), np.nan)
    if len(trace) >= minimum_localizations and np.isfinite(before):
        for node in range(len(trace)):
            after = _top3(
                [value for left, right, value in edges if node not in (left, right)]
            )
            if np.isfinite(after):
                scores[node] = before - after
    result = trace.copy()
    result["loo_top3_all"] = scores
    result["score_status"] = np.where(
        np.isfinite(scores), "ok", "insufficient_reference"
    )
    metrics = confidence_metrics(scores)
    valid = np.flatnonzero(np.isfinite(scores))
    winner = int(valid[np.argmax(scores[valid])]) if len(valid) else None
    summary = {
        "C_top3": before,
        **metrics,
        "top_candidate_barcode": barcodes[winner] if winner is not None else math.nan,
        "top_candidate_spot_id": str(result.iloc[winner]["Spot_ID"])
        if winner is not None
        else "",
        "n_observed_localizations": len(result),
        "n_scorable_localizations": len(valid),
        "n_true_corruptions": int(result.get("is_corrupted", False).astype(bool).sum()),
    }
    return result, summary


def iterative_scores(
    trace: pd.DataFrame, reference, maximum_steps: int = 3
) -> list[dict]:
    """Remove and rescore with a fixed reference, retaining explicit states."""
    current = trace.copy()
    rows = []
    for step in range(maximum_steps):
        if len(current) < 3:
            break
        scored, summary = score_trace(current, reference)
        candidate = scored.loc[
            scored["Spot_ID"].astype(str) == str(summary["top_candidate_spot_id"])
        ]
        if candidate.empty or not np.isfinite(summary["S1"]):
            break
        row = candidate.iloc[0]
        rows.append(
            {
                "iteration": step,
                **summary,
                "predicted_spot_id": str(row["Spot_ID"]),
                "predicted_barcode": row["Barcode"],
                "predicted_is_true_corruption": bool(row.get("is_corrupted", False)),
                "remaining_true_corruptions_before_step": int(
                    current.get("is_corrupted", False).astype(bool).sum()
                ),
                "n_remaining_localizations": len(current),
            }
        )
        current = current.loc[
            current["Spot_ID"].astype(str) != str(row["Spot_ID"])
        ].copy()
    return rows


def step_null_scores(
    clean: pd.DataFrame, reference, maximum_steps: int = 3
) -> list[dict]:
    rows = []
    for key, trace in clean.groupby(["simulation_id", "Trace_ID"], sort=False):
        for row in iterative_scores(trace, reference, maximum_steps):
            rows.append(
                {
                    "simulation_id": key[0],
                    "Trace_ID": key[1],
                    "step": row["iteration"],
                    "S1": row["S1"],
                }
            )
    return rows


def _condition_metadata(item: Mapping) -> tuple[int, int]:
    n_barcodes = item.get("n_barcodes", item.get("barcode_count"))
    k = item.get(
        "corrupted_barcodes_per_trace", item.get("n_corruptions", item.get("k"))
    )
    if n_barcodes is None or k is None:
        raise ValueError(
            f"Condition {item.get('id')} lacks n_barcodes or corrupted_barcodes_per_trace"
        )
    return int(n_barcodes), int(k)


def _load_observations(cycle1, root: Path, item: Mapping) -> pd.DataFrame:
    frame = cycle1.load_observations(root / str(item["directory"]), item)
    truth = cycle1._frame(root / str(item["directory"]) / "simulated.ground_truth.ecsv")
    extras = [
        name
        for name in ("corruption_target_rank", "realized_displacement_um")
        if name in truth
    ]
    if extras:
        frame = frame.merge(
            truth[["Spot_ID", *extras]], on="Spot_ID", how="left", validate="one_to_one"
        )
    n_barcodes, k = _condition_metadata(item)
    frame["n_barcodes"] = n_barcodes
    frame["K"] = k
    return frame


def _trace_metadata(item: Mapping) -> dict:
    n_barcodes, k = _condition_metadata(item)
    return {
        "condition": str(item["id"]),
        "seed": int(item["seed"]),
        "n_barcodes": n_barcodes,
        "K": k,
        "detection_efficiency": float(item["detection_efficiency"]),
        "displacement": float(item["displacement_um"]),
    }


def _evaluate_trace(
    trace: pd.DataFrame, reference, item: Mapping, thresholds: Mapping
) -> tuple[dict, list[dict], list[dict], list[dict], pd.DataFrame]:
    meta = _trace_metadata(item)
    ident = {
        "simulation_id": str(trace["simulation_id"].iloc[0]),
        "Trace_ID": str(trace["Trace_ID"].iloc[0]),
    }
    scored, summary = score_trace(trace, reference)
    valid = scored[np.isfinite(scored["loo_top3_all"])].copy()
    valid = valid.sort_values("loo_top3_all", ascending=False)
    values = valid["loo_top3_all"].to_numpy(float)
    valid["rank"] = [1 + int(np.sum(values > value)) for value in values]
    truths = valid[valid["is_corrupted"].astype(bool)]
    true_ranks = truths["rank"].to_numpy(int)
    top_score = float(values[0]) if len(values) else math.nan
    n_top_tied = int(np.sum(values == top_score)) if len(values) else 0
    selected_trace = bool(
        trace.get("selected_for_corruption", False).astype(bool).any()
    )
    confidence = {
        **meta,
        **ident,
        **summary,
        # A tied maximum is ambiguous, not a correct attribution chosen by row order.
        "top_candidate_is_true": bool(
            len(valid) and summary["winner_is_unique"] and valid.iloc[0]["is_corrupted"]
        ),
        "scoreable": bool(len(valid)),
        "true_target_ranks": ",".join(map(str, true_ranks.tolist())),
        "n_tied_for_top_score": n_top_tied,
        "selected_trace": selected_trace,
    }
    primary = thresholds.get((0, 0.01), math.nan)
    if meta["K"] == 1 and selected_trace:
        target_rank = int(true_ranks[0]) if len(true_ranks) == 1 else math.nan
        target_score = (
            float(truths.iloc[0]["loo_top3_all"]) if len(truths) == 1 else math.nan
        )
        n_tied_at_target = (
            int(np.sum(values == target_score)) if np.isfinite(target_score) else 0
        )
        single = {
            **meta,
            **ident,
            "threshold": primary,
            "target_sensitivity": bool(
                len(truths) and truths.iloc[0]["loo_top3_all"] > primary
            ),
            "unique_top1_attribution": bool(
                len(valid)
                and valid.iloc[0]["is_corrupted"]
                and summary["winner_is_unique"]
            ),
            "top1_unique": bool(
                len(true_ranks) == 1 and target_rank == 1 and n_tied_at_target == 1
            ),
            "top2_conservative_attribution": bool(
                len(true_ranks) == 1 and int(np.sum(values >= target_score)) <= 2
            ),
            "top2_conservative": bool(
                len(true_ranks) == 1 and int(np.sum(values >= target_score)) <= 2
            ),
            "n_tied_at_target_rank": n_tied_at_target,
            "true_target_rank": target_rank,
            "reciprocal_rank": 1 / target_rank if np.isfinite(target_rank) else 0.0,
            "collateral_call_rate": float(
                (
                    valid.loc[~valid["is_corrupted"].astype(bool), "loo_top3_all"]
                    > primary
                ).mean()
            )
            if len(valid)
            else math.nan,
            "attribution_coverage": bool(len(truths) == 1),
        }
        for fpr in FIXED_FPRS:
            single[f"target_sensitivity_fpr_{fpr:g}"] = bool(
                len(truths)
                and truths.iloc[0]["loo_top3_all"] > thresholds.get((0, fpr), math.nan)
            )
    else:
        single = {}
    n_true = int(trace["is_corrupted"].astype(bool).sum())
    k = meta["K"]

    def top_n_fraction(count: int, include_ties: bool) -> float:
        if not n_true or len(valid) < count:
            return math.nan
        cutoff = values[count - 1]
        chosen = valid[valid["loo_top3_all"] >= cutoff]
        if not include_ties and len(chosen) > count:
            return math.nan
        return float(chosen["is_corrupted"].astype(bool).sum() / n_true)

    if len(valid) >= k:
        boundary = float(values[k - 1])
        boundary_group = valid[valid["loo_top3_all"] == boundary]
        including_ties = valid[valid["loo_top3_all"] >= boundary]
        boundary_is_tied = len(including_ties) > k
        conservative = including_ties if not boundary_is_tied else valid.iloc[0:0]
    else:
        boundary_group = including_ties = conservative = valid.iloc[0:0]
        boundary_is_tied = False
    single_pass = {
        **meta,
        **ident,
        "fraction_true_in_top1": top_n_fraction(1, False),
        "fraction_true_in_top2": top_n_fraction(2, False),
        "fraction_true_in_top3": top_n_fraction(3, False),
        "fraction_true_in_top1_including_ties": top_n_fraction(1, True),
        "fraction_true_in_top2_including_ties": top_n_fraction(2, True),
        "fraction_true_in_top3_including_ties": top_n_fraction(3, True),
        "precision_top_k": float(conservative["is_corrupted"].astype(bool).mean())
        if len(conservative)
        else math.nan,
        "recall_top_k": float(conservative["is_corrupted"].astype(bool).sum() / n_true)
        if n_true
        else math.nan,
        "precision_top_k_including_ties": float(
            including_ties["is_corrupted"].astype(bool).mean()
        )
        if len(including_ties)
        else math.nan,
        "recall_top_k_including_ties": float(
            including_ties["is_corrupted"].astype(bool).sum() / n_true
        )
        if n_true
        else math.nan,
        "top_k_boundary_is_tied": boundary_is_tied,
        "n_tied_at_k_boundary": len(boundary_group),
        "top_k_conservative": not boundary_is_tied and len(conservative) == k,
        "exact_set_recovery_top_k": bool(
            not boundary_is_tied
            and len(conservative) == n_true
            and conservative["is_corrupted"].astype(bool).all()
        ),
        "mean_true_rank": float(np.mean(true_ranks)) if len(true_ranks) else math.nan,
    }
    iterative = [
        {**meta, **ident, "selected_trace": selected_trace, **row}
        for row in iterative_scores(trace, reference)
    ]
    localization = scored.assign(**meta, **ident)
    return (
        single,
        ([single_pass] if selected_trace else []),
        iterative,
        [confidence],
        localization,
    )


def _summarize_iterative(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return rows
    keys = [
        "condition",
        "seed",
        "n_barcodes",
        "K",
        "detection_efficiency",
        "displacement",
    ]
    result = []
    for key, group in rows.groupby(keys, sort=False):
        trace_keys = ["simulation_id", "Trace_ID"]
        for trace_key, trace in group.groupby(trace_keys, sort=False):
            if not bool(trace["selected_trace"].iloc[0]):
                continue
            trace = trace.sort_values("iteration")
            truths = trace["predicted_is_true_corruption"].astype(bool).to_numpy()
            k = int(trace["K"].iloc[0])
            recovered = np.cumsum(truths)
            record = {
                **dict(zip(keys, key)),
                "simulation_id": trace_key[0],
                "Trace_ID": trace_key[1],
            }
            for step in range(3):
                record[f"precision_step_{step + 1}"] = (
                    bool(truths[step]) if step < len(truths) else math.nan
                )
                record[f"cumulative_recall_step_{step + 1}"] = (
                    min(int(recovered[step]), k) / k
                    if step < len(recovered)
                    else (int(recovered[-1]) / k if len(recovered) else 0.0)
                )
            record["exact_recovery_after_k"] = bool(
                len(truths) >= k and truths[:k].all()
            )
            record["first_k_all_true"] = record["exact_recovery_after_k"]
            positions = np.flatnonzero(recovered >= k)
            record["false_selections_before_all_recovered"] = (
                int((positions[0] + 1 - k)) if len(positions) else math.nan
            )
            result.append(record)
    return pd.DataFrame(result)


def _thresholded(rows: pd.DataFrame, thresholds: Mapping) -> pd.DataFrame:
    result = []
    keys = [
        "condition",
        "seed",
        "n_barcodes",
        "K",
        "detection_efficiency",
        "displacement",
        "simulation_id",
        "Trace_ID",
    ]
    for key, trace in rows.groupby(keys, sort=False):
        trace = trace.sort_values("iteration")
        selected = []
        stop = 0
        stop_reason = "maximum_steps"
        for row in trace.itertuples():
            stop = int(row.iteration)
            if row.S1 <= thresholds.get((int(row.iteration), 0.01), math.nan):
                stop_reason = "below_threshold"
                break
            if not bool(row.winner_is_unique):
                stop_reason = "ambiguous_tie"
                break
            selected.append(bool(row.predicted_is_true_corruption))
            stop = int(row.iteration) + 1
        k = int(trace["remaining_true_corruptions_before_step"].iloc[0])
        true = sum(selected)
        false = len(selected) - true
        result.append(
            {
                **dict(zip(keys, key)),
                "selected_trace": bool(trace["selected_trace"].iloc[0]),
                "target_recall": true / k if k else math.nan,
                "false_removal_rate": false / max(len(selected), 1),
                "over_removal_after_all_true_gone": max(0, len(selected) - k)
                if true == k
                else 0,
                "under_removal_remaining_missed": k - true,
                "stopped_at_iteration": stop,
                "stop_reason": stop_reason,
                "n_provisional_removals": len(selected),
            }
        )
    return pd.DataFrame(result)


def confidence_precision_coverage(
    frame: pd.DataFrame, group_columns: Sequence[str] = ()
) -> pd.DataFrame:
    """Descriptive curves only; quantiles are reported, never selected as rules."""
    rows = []
    groups = (
        frame.groupby(list(group_columns), sort=False, dropna=False)
        if group_columns
        else [((), frame)]
    )
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        metadata = dict(zip(group_columns, key))
        for metric in (
            "S1",
            "absolute_margin",
            "relative_margin",
            "score_ratio",
            "C_top3",
        ):
            values = group[metric].to_numpy(float)
            finite = np.isfinite(values)
            for quantile in np.linspace(0, 0.95, 20):
                cutoff = (
                    float(np.quantile(values[finite], quantile))
                    if finite.any()
                    else math.nan
                )
                kept = finite & (values >= cutoff)
                rows.append(
                    {
                        **metadata,
                        "metric": metric,
                        "quantile": quantile,
                        "cutoff": cutoff,
                        "coverage": float(kept.mean()),
                        "precision": float(
                            group.loc[kept, "top_candidate_is_true"].astype(bool).mean()
                        )
                        if kept.any()
                        else math.nan,
                        "n": int(kept.sum()),
                    }
                )
    return pd.DataFrame(rows)


def _aggregate_condition(frame: pd.DataFrame, metadata: Mapping) -> dict:
    """Collapse a condition-sized per-trace frame immediately."""
    record = {**metadata, "n_traces": len(frame)}
    if frame.empty:
        return record
    for column in frame.select_dtypes(include=["number", "bool"]).columns:
        if column not in metadata:
            record[column] = frame[column].mean()
    return record


def confidence_distribution_summary(
    frame: pd.DataFrame, metadata: Mapping
) -> list[dict]:
    rows = []
    for correct, group in frame.groupby("top_candidate_is_true", dropna=False):
        for metric in ("S1", "S2", "absolute_margin"):
            values = group[metric].dropna().to_numpy(float)
            for quantile in np.linspace(0, 1, 21):
                rows.append(
                    {
                        **metadata,
                        "top_candidate_is_true": bool(correct),
                        "metric": metric,
                        "quantile": quantile,
                        "value": float(np.quantile(values, quantile))
                        if len(values)
                        else math.nan,
                        "n": len(values),
                    }
                )
    return rows


def run(
    root: Path,
    output: Path,
    minimum_observations: int = 20,
    write_localization_scores: bool = False,
) -> None:
    cycle1, cycle2 = _cycle(1), _cycle(2)
    manifest = yaml.safe_load((root / "sweep_manifest.yaml").read_text())
    if (
        manifest.get("dry_run")
        or not manifest.get("simulations")
        or not manifest.get("conditions")
    ):
        raise ValueError(
            "Manifest must be a completed sweep with simulations and conditions"
        )
    output.mkdir(parents=True, exist_ok=True)
    baselines = {}
    simulation_sizes = {
        str(item["simulation_id"]): _condition_metadata(item)[0]
        for item in manifest["conditions"]
    }
    for item in manifest["simulations"]:
        frame = cycle1.load_baseline(root / "simulations" / str(item["id"]), item)
        n = int(
            item.get(
                "n_barcodes",
                simulation_sizes.get(str(item["id"]), frame["Barcode"].max()),
            )
        )
        frame["n_barcodes"] = n
        baselines[str(item["id"])] = frame
    groups: dict[tuple[int, float], list[Mapping]] = defaultdict(list)
    for item in manifest["conditions"]:
        n, _ = _condition_metadata(item)
        groups[(n, float(item["detection_efficiency"]))].append(item)
    single_rows = []
    multi_rows = []
    iterative_rows = []
    thresholded_rows = []
    precision_rows = []
    distribution_rows = []
    clean_evaluation_rows = []
    reference_rows = []
    threshold_rows = []
    writer = (
        cycle1._EcsvChunkWriter(output / "loo_localization_scores.ecsv")
        if write_localization_scores
        else None
    )
    detail_writer = cycle1._EcsvChunkWriter(output / "loo_iterative_diagnostics.ecsv")
    confidence_writer = cycle1._EcsvChunkWriter(output / "loo_confidence_metrics.ecsv")
    for (n_barcodes, efficiency), items in groups.items():
        seeds = sorted({int(item["seed"]) for item in items})
        if len(seeds) < 3:
            raise ValueError("At least three replicate seeds are required")
        baseline = pd.concat(
            [
                f
                for f in baselines.values()
                if int(f["n_barcodes"].iloc[0]) == n_barcodes
                and np.isclose(f["detection_efficiency"].iloc[0], efficiency)
            ],
            ignore_index=True,
        )
        for evaluation_seed in seeds:
            ref_seeds, _, eval_seeds = fold_seed_sets(seeds, evaluation_seed)
            reference = cycle2.fit_reference(
                baseline[baseline["seed"].isin(ref_seeds)], minimum_observations
            )
            for separation, values in reference.separation.items():
                reference_rows.append(
                    {
                        "n_barcodes": n_barcodes,
                        "detection_efficiency": efficiency,
                        "evaluation_seed": evaluation_seed,
                        "reference_seeds": ",".join(map(str, sorted(ref_seeds))),
                        "calibration_seeds": ",".join(
                            map(str, sorted(set(seeds) - eval_seeds))
                        ),
                        "evaluation_seeds": str(evaluation_seed),
                        "genomic_separation": separation,
                        "n_observations": len(values),
                        "sufficient": len(values) >= minimum_observations,
                    }
                )
            nulls: dict[int, list[float]] = defaultdict(list)
            for calibration_seed in seeds:
                if calibration_seed == evaluation_seed:
                    continue
                nested, cal, ev = fold_seed_sets(
                    seeds, evaluation_seed, calibration_seed
                )
                nested_ref = cycle2.fit_reference(
                    baseline[baseline["seed"].isin(nested)], minimum_observations
                )
                clean = baseline[baseline["seed"] == calibration_seed]
                for row in step_null_scores(clean, nested_ref):
                    nulls[row["step"]].append(row["S1"])
                del nested_ref, clean
            thresholds = {
                (step, fpr): trace_threshold(values, fpr)
                for step, values in nulls.items()
                for fpr in FIXED_FPRS
            }
            for (step, fpr), value in thresholds.items():
                threshold_rows.append(
                    {
                        "n_barcodes": n_barcodes,
                        "detection_efficiency": efficiency,
                        "evaluation_seed": evaluation_seed,
                        "step": step,
                        "trace_fpr": fpr,
                        "threshold": value,
                        "n_clean_traces": len(nulls[step]),
                        "reference_seeds": ",".join(map(str, sorted(ref_seeds))),
                        "calibration_seeds": ",".join(
                            map(str, sorted(set(seeds) - eval_seeds))
                        ),
                        "evaluation_seeds": str(evaluation_seed),
                    }
                )
            for item in items:
                if int(item["seed"]) != evaluation_seed:
                    continue
                _, k = _condition_metadata(item)
                _progress(
                    n_barcodes=n_barcodes,
                    K=k,
                    efficiency=efficiency,
                    evaluation_seed=evaluation_seed,
                    condition=item["id"],
                    iteration="0-2",
                )
                observations = _load_observations(cycle1, root, item)
                condition_single = []
                condition_multi = []
                condition_iterative = []
                condition_confidence = []
                for _, trace in observations.groupby(
                    ["simulation_id", "Trace_ID"], sort=False
                ):
                    single, multi, iteration, confidence, localization = (
                        _evaluate_trace(trace, reference, item, thresholds)
                    )
                    if single:
                        condition_single.append(single)
                    condition_multi.extend(multi)
                    condition_iterative.extend(iteration)
                    condition_confidence.extend(confidence)
                    if writer is not None:
                        writer.write(localization)
                metadata = _trace_metadata(item)
                single_frame = pd.DataFrame(condition_single)
                multi_frame = pd.DataFrame(condition_multi)
                iterative_frame = pd.DataFrame(condition_iterative)
                confidence_frame = pd.DataFrame(condition_confidence)
                if not single_frame.empty:
                    single_rows.append(_aggregate_condition(single_frame, metadata))
                if not multi_frame.empty:
                    multi_rows.append(_aggregate_condition(multi_frame, metadata))
                if not iterative_frame.empty:
                    detail_writer.write(iterative_frame)
                    per_trace_iterative = _summarize_iterative(iterative_frame)
                    if not per_trace_iterative.empty:
                        iterative_rows.append(
                            _aggregate_condition(per_trace_iterative, metadata)
                        )
                    per_trace_thresholded = _thresholded(iterative_frame, thresholds)
                    thresholded_record = _aggregate_condition(
                        per_trace_thresholded, metadata
                    )
                    thresholded_record["trace_fpr"] = 0.01
                    clean_stopping = per_trace_thresholded[
                        ~per_trace_thresholded["selected_trace"].astype(bool)
                    ]
                    thresholded_record["fraction_clean_traces_with_any_removal"] = (
                        float((clean_stopping["n_provisional_removals"] > 0).mean())
                        if len(clean_stopping)
                        else math.nan
                    )
                    thresholded_record["mean_false_removals_per_clean_trace"] = (
                        float(clean_stopping["n_provisional_removals"].mean())
                        if len(clean_stopping)
                        else math.nan
                    )
                    thresholded_rows.append(thresholded_record)
                if not confidence_frame.empty:
                    confidence_writer.write(confidence_frame)
                    precision_rows.extend(
                        confidence_precision_coverage(
                            confidence_frame,
                            ("n_barcodes", "K", "detection_efficiency", "displacement"),
                        ).to_dict("records")
                    )
                    distribution_rows.extend(
                        confidence_distribution_summary(confidence_frame, metadata)
                    )
                    clean_confidence = confidence_frame[
                        ~confidence_frame["selected_trace"].astype(bool)
                    ]
                    for fpr in FIXED_FPRS:
                        clean_record = {
                            **metadata,
                            "trace_fpr": fpr,
                            "threshold": thresholds.get((0, fpr), math.nan),
                            "clean_trace_flag_rate": float(
                                (
                                    clean_confidence["S1"]
                                    > thresholds.get((0, fpr), math.nan)
                                ).mean()
                            )
                            if len(clean_confidence)
                            else math.nan,
                            "n_clean_traces": len(clean_confidence),
                        }
                        if not iterative_frame.empty:
                            clean_record["stopping_trace_fpr"] = 0.01
                            clean_record["fraction_clean_traces_with_any_removal"] = (
                                thresholded_record[
                                    "fraction_clean_traces_with_any_removal"
                                ]
                            )
                            clean_record["mean_false_removals_per_clean_trace"] = (
                                thresholded_record[
                                    "mean_false_removals_per_clean_trace"
                                ]
                            )
                        clean_evaluation_rows.append(clean_record)
                del observations
                del single_frame, multi_frame, iterative_frame, confidence_frame
                gc.collect()
    if writer is not None:
        writer.finish()
    detail_writer.finish()
    confidence_writer.finish()
    single = pd.DataFrame(single_rows)
    multi = pd.DataFrame(multi_rows)
    iterative = pd.DataFrame(iterative_rows)
    thresholds_frame = pd.DataFrame(threshold_rows)
    thresholded = pd.DataFrame(thresholded_rows)
    precision = pd.DataFrame(precision_rows)
    for frame, name in (
        (single, OUTPUTS[0]),
        (multi, OUTPUTS[1]),
        (iterative, OUTPUTS[2]),
        (thresholded, OUTPUTS[3]),
        (precision, OUTPUTS[5]),
        (thresholds_frame, OUTPUTS[6]),
        (pd.DataFrame(reference_rows), OUTPUTS[7]),
        (pd.DataFrame(clean_evaluation_rows), OUTPUTS[8]),
        (pd.DataFrame(distribution_rows), OUTPUTS[9]),
    ):
        _write(frame, output / name)
    plot_results(output)
    _progress(status="complete", output=output)


def _read(path: Path) -> pd.DataFrame:
    return Table.read(path, format="ascii.ecsv").to_pandas()


def plot_results(output: Path) -> None:
    """Create decision plots solely from compact ECSV outputs."""
    missing = [name for name in OUTPUTS if not (output / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing compact output(s): {', '.join(missing)}")
    single = _read(output / OUTPUTS[0])
    iterative = _read(output / OUTPUTS[2])
    stopping = _read(output / OUTPUTS[3])
    curve = _read(output / OUTPUTS[5])
    distribution = _read(output / OUTPUTS[9])

    def line_plot(frame, x, y, hue, name, ylabel):
        fig, ax = plt.subplots(figsize=(8, 5))
        if not frame.empty:
            for key, g in frame.groupby(hue, sort=False):
                means = g.groupby(x)[y].mean().sort_index()
                ax.plot(means.index, means.values, marker="o", label=str(key))
        ax.set(xlabel=x, ylabel=ylabel)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(output / name, dpi=160)
        plt.close(fig)

    line_plot(
        single,
        "displacement",
        "target_sensitivity",
        ["n_barcodes", "detection_efficiency"],
        "k1_sensitivity.png",
        "Sensitivity at 1% trace FPR",
    )
    line_plot(
        single,
        "displacement",
        "unique_top1_attribution",
        ["n_barcodes", "detection_efficiency"],
        "k1_unique_top1.png",
        "Unique top-1 accuracy",
    )
    long = []
    for row in iterative.itertuples():
        for step in range(1, 4):
            long.append(
                {
                    "K": row.K,
                    "step": step,
                    "recall": getattr(row, f"cumulative_recall_step_{step}"),
                    "precision": getattr(row, f"precision_step_{step}"),
                }
            )
    long = pd.DataFrame(long)
    line_plot(
        long, "step", "recall", ["K"], "iterative_recall.png", "Cumulative recall"
    )
    line_plot(
        long, "step", "precision", ["K"], "iterative_precision.png", "Step precision"
    )
    line_plot(
        iterative,
        "K",
        "exact_recovery_after_k",
        ["n_barcodes"],
        "exact_recovery.png",
        "Exact recovery",
    )
    line_plot(
        stopping,
        "K",
        "false_removal_rate",
        ["n_barcodes"],
        "thresholded_false_removals.png",
        "False-removal rate",
    )
    line_plot(
        curve,
        "coverage",
        "precision",
        ["metric"],
        "precision_coverage.png",
        "Precision",
    )
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, metric in zip(axes, ("S1", "S2", "absolute_margin")):
        metric_frame = distribution[distribution["metric"] == metric]
        for correct, g in metric_frame.groupby("top_candidate_is_true"):
            line = g.groupby("quantile")["value"].mean().sort_index()
            ax.plot(line.values, line.index, label=f"correct={correct}")
        ax.set(xlabel=metric, ylabel="Empirical quantile")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output / "winner_confidence_distributions.png", dpi=160)
    plt.close(fig)


def plot_only(output: Path) -> None:
    plot_results(output)
    _progress(status="plots complete", output=output)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-root", type=Path, default=Path("benchmark_curator_loo_validation")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-reference-observations", type=int, default=20)
    parser.add_argument("--write-localization-scores", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_arguments()
    if args.plot_only:
        plot_only(args.output_dir.resolve())
    else:
        run(
            args.benchmark_root.resolve(),
            args.output_dir.resolve(),
            args.minimum_reference_observations,
            args.write_localization_scores,
        )


if __name__ == "__main__":
    main()
