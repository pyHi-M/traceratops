#!/usr/bin/env python3
"""Compare interpretable anomaly scores on a trace-curator benchmark.

This is an exploratory benchmark utility, not a production ``trace_curator``
command.  Reference models and operating-point thresholds are learned without
using the replicate being evaluated.  See
``docs/source/contribute/trace_curator_benchmark_analysis.md`` for details.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table
from sklearn.metrics import average_precision_score, roc_auc_score

CONTEXTS = ("k1", "k2", "k3", "all")
FIXED_FPRS = (0.001, 0.01, 0.05)
REQUIRED_MAIN = {"Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z"}
REQUIRED_TRUTH = {
    "Spot_ID",
    "Input_Trace_ID",
    "Barcode #",
    "selected_for_corruption",
    "is_corrupted",
    "injected_displacement_um",
}


def _frame(path: Path) -> pd.DataFrame:
    """Read an ECSV into pandas without losing string identifiers."""
    return Table.read(path, format="ascii.ecsv").to_pandas()


def _require(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")


def load_observations(condition_dir: Path, metadata: Mapping) -> pd.DataFrame:
    """Load and strictly join one condition's observed and evaluation tables."""
    main_path = condition_dir / "simulated.ecsv"
    truth_path = condition_dir / "simulated.ground_truth.ecsv"
    main, truth = _frame(main_path), _frame(truth_path)
    _require(main, REQUIRED_MAIN, main_path)
    _require(truth, REQUIRED_TRUTH, truth_path)
    if main["Spot_ID"].duplicated().any() or truth["Spot_ID"].duplicated().any():
        raise ValueError(f"Spot_ID must be unique in {condition_dir}")
    truth_columns = list(REQUIRED_TRUTH)
    joined = main.merge(
        truth[truth_columns],
        on="Spot_ID",
        how="left",
        validate="one_to_one",
        suffixes=("", "_truth"),
        indicator=True,
    )
    if len(joined) != len(main) or (joined["_merge"] != "both").any():
        raise ValueError(f"Main/truth Spot_ID mismatch in {condition_dir}")
    if not (
        joined["Trace_ID"].astype(str) == joined["Input_Trace_ID"].astype(str)
    ).all():
        raise ValueError(f"Trace identity mismatch in {condition_dir}")
    if not (joined["Barcode #"] == joined["Barcode #_truth"]).all():
        raise ValueError(f"Barcode identity mismatch in {condition_dir}")
    if joined.duplicated(["Trace_ID", "Barcode #"]).any():
        raise ValueError(f"Duplicate barcode within a trace in {condition_dir}")
    joined = joined.drop(columns=["_merge", "Input_Trace_ID", "Barcode #_truth"])
    joined = joined.rename(columns={"Barcode #": "Barcode"})
    joined["condition"] = str(metadata["id"])
    joined["simulation_id"] = str(metadata["simulation_id"])
    joined["replicate"] = int(metadata["replicate_index"])
    joined["seed"] = int(metadata["seed"])
    joined["detection_efficiency"] = float(metadata["detection_efficiency"])
    joined["displacement"] = float(metadata["displacement_um"])
    return joined


def load_baseline(path: Path, metadata: Mapping) -> pd.DataFrame:
    """Load one clean simulator output for reference fitting."""
    frame = _frame(path / "simulated.ecsv")
    _require(frame, REQUIRED_MAIN, path / "simulated.ecsv")
    if frame.duplicated(["Trace_ID", "Barcode #"]).any():
        raise ValueError(f"Reference contains duplicate barcodes: {path}")
    frame = frame.rename(columns={"Barcode #": "Barcode"})
    frame["simulation_id"] = str(metadata["id"])
    frame["seed"] = int(metadata["seed"])
    frame["detection_efficiency"] = float(metadata["detection_efficiency"])
    # Explicit provenance makes accidental use of corrupted conditions testable.
    frame["is_corrupted"] = False
    return frame


def clean_reference_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Return only rows eligible for reference fitting."""
    if "is_corrupted" not in frame:
        raise ValueError("Reference rows require an is_corrupted provenance column")
    return frame.loc[~frame["is_corrupted"].astype(bool)].copy()


def select_context(barcodes: Sequence[int], target: int, context: str) -> list[int]:
    """Select nearest genomic neighbours on each side of ``target``."""
    if context not in CONTEXTS:
        raise ValueError(f"Unknown context {context!r}; expected one of {CONTEXTS}")
    values = np.asarray(barcodes, dtype=int)
    candidates = np.flatnonzero(values != int(target))
    if context == "all":
        return candidates.tolist()
    k = int(context[1:])
    left = candidates[values[candidates] < target]
    right = candidates[values[candidates] > target]
    left = left[np.argsort(target - values[left])][:k]
    right = right[np.argsort(values[right] - target)][:k]
    return np.concatenate((left, right)).tolist()


def empirical_tail_probability(values: Sequence[float], value: float) -> float:
    """Return a finite-sample, two-sided empirical tail probability.

    Add-one smoothing is applied separately to ``P(X <= value)`` and
    ``P(X >= value)``; consequently the result is positive and at most one.
    """
    sample = np.sort(np.asarray(values, dtype=float))
    if not len(sample) or not np.isfinite(value):
        return math.nan
    lower = (np.searchsorted(sample, value, side="right") + 1) / (len(sample) + 1)
    upper = (len(sample) - np.searchsorted(sample, value, side="left") + 1) / (
        len(sample) + 1
    )
    return float(min(1.0, 2.0 * min(lower, upper)))


@dataclass
class ReferenceModel:
    """Clean empirical pair-distance distributions."""

    separation: dict[int, np.ndarray]
    barcode_pair: dict[tuple[int, int], np.ndarray]
    residual: dict[int, tuple[float, float, int]]
    minimum_observations: int = 20

    def values(self, strategy: str, left: int, right: int) -> np.ndarray | None:
        key = (
            abs(int(right) - int(left))
            if strategy == "separation"
            else tuple(sorted((int(left), int(right))))
        )
        values = getattr(self, strategy).get(key)
        return (
            values
            if values is not None and len(values) >= self.minimum_observations
            else None
        )

    def residual_stats(self, separation: int) -> tuple[float, float, int] | None:
        """Return the requested or nearest well-supported residual bin."""
        populated = {
            key: value
            for key, value in self.residual.items()
            if value[2] >= self.minimum_observations
        }
        if not populated:
            return None
        separation = abs(int(separation))
        key = (
            separation
            if separation in populated
            else min(populated, key=lambda candidate: abs(candidate - separation))
        )
        return populated[key]


def fit_reference(
    frame: pd.DataFrame, minimum_observations: int = 20
) -> ReferenceModel:
    """Fit pair distributions using explicitly clean rows only."""
    clean = clean_reference_rows(frame)
    by_sep: dict[int, list[float]] = defaultdict(list)
    by_pair: dict[tuple[int, int], list[float]] = defaultdict(list)
    for _, trace in clean.groupby("Trace_ID", sort=False):
        barcodes = trace["Barcode"].to_numpy(dtype=int)
        coords = trace[["x", "y", "z"]].to_numpy(dtype=float)
        for i in range(len(trace)):
            for j in range(i + 1, len(trace)):
                distance = float(np.linalg.norm(coords[i] - coords[j]))
                by_sep[abs(int(barcodes[j]) - int(barcodes[i]))].append(distance)
                by_pair[tuple(sorted((int(barcodes[i]), int(barcodes[j]))))].append(
                    distance
                )
    separation = {key: np.asarray(value) for key, value in by_sep.items()}
    barcode_pair = {key: np.asarray(value) for key, value in by_pair.items()}
    residual = {
        key: (float(np.mean(value)), max(float(np.std(value)), 1e-6), len(value))
        for key, value in separation.items()
    }
    return ReferenceModel(separation, barcode_pair, residual, minimum_observations)


def reference_summary(model: ReferenceModel, fold: Mapping) -> pd.DataFrame:
    rows = []
    for strategy in ("separation", "barcode_pair"):
        for key, values in getattr(model, strategy).items():
            rows.append(
                {
                    **fold,
                    "reference_strategy": strategy,
                    "reference_key": str(key),
                    "n_observations": len(values),
                    "sufficient": len(values) >= model.minimum_observations,
                }
            )
    return pd.DataFrame(rows)


def _aggregate(anomalies: np.ndarray) -> dict[str, float]:
    ordered = np.sort(anomalies)
    return {
        "mean": float(np.mean(anomalies)),
        "median": float(np.median(anomalies)),
        "maximum": float(ordered[-1]),
        "second_largest": float(ordered[-2]) if len(ordered) >= 2 else math.nan,
        "fraction_95": float(np.mean(anomalies >= -math.log10(0.05))),
        "count_95": float(np.sum(anomalies >= -math.log10(0.05))),
        "fraction_99": float(np.mean(anomalies >= -math.log10(0.01))),
        "count_99": float(np.sum(anomalies >= -math.log10(0.01))),
    }


def score_observations(frame: pd.DataFrame, model: ReferenceModel) -> pd.DataFrame:
    """Score all localizations using coordinates and barcode identities only."""
    output = []
    feature_columns = ["Trace_ID", "Spot_ID", "Barcode", "x", "y", "z"]
    # Copy only feature columns before scoring: evaluation labels cannot affect scores.
    features = frame[feature_columns]
    diagnostics = frame.drop(columns=["x", "y", "z"]).set_index("Spot_ID")
    for _, trace in features.groupby("Trace_ID", sort=False):
        barcodes = trace["Barcode"].to_numpy(dtype=int)
        coords = trace[["x", "y", "z"]].to_numpy(dtype=float)
        for i, row in enumerate(trace.itertuples(index=False)):
            barcode = int(row.Barcode)
            left_values = barcodes[barcodes < barcode]
            right_values = barcodes[barcodes > barcode]
            common = {
                "nearest_left_barcode": (
                    int(left_values.max()) if len(left_values) else math.nan
                ),
                "nearest_right_barcode": (
                    int(right_values.min()) if len(right_values) else math.nan
                ),
                "left_genomic_gap": (
                    int(barcode - left_values.max()) if len(left_values) else math.nan
                ),
                "right_genomic_gap": (
                    int(right_values.min() - barcode) if len(right_values) else math.nan
                ),
                "has_left_flank": bool(len(left_values)),
                "has_right_flank": bool(len(right_values)),
            }
            metadata = diagnostics.loc[str(row.Spot_ID)].to_dict()
            metadata.update(
                {
                    "Trace_ID": str(row.Trace_ID),
                    "Spot_ID": str(row.Spot_ID),
                    "Barcode": barcode,
                }
            )
            for context in CONTEXTS:
                selected = select_context(barcodes, barcode, context)
                base = {
                    **metadata,
                    **common,
                    "context": context,
                    "n_context_requested": len(selected),
                }
                residuals = []
                residual_supports = []
                for j in selected:
                    stats = model.residual_stats(abs(barcode - int(barcodes[j])))
                    if stats:
                        distance = float(np.linalg.norm(coords[i] - coords[j]))
                        residuals.append(((distance - stats[0]) / stats[1]) ** 2)
                        residual_supports.append(stats[2])
                output.append(
                    {
                        **base,
                        "model": "spatial_residual",
                        "score": float(np.mean(residuals)) if residuals else math.nan,
                        "n_context": len(residuals),
                        "min_reference_n": min(residual_supports, default=0),
                        "score_status": "ok" if residuals else "insufficient_reference",
                    }
                )
                for strategy in ("separation", "barcode_pair"):
                    anomalies, supports = [], []
                    for j in selected:
                        values = model.values(strategy, barcode, int(barcodes[j]))
                        if values is None:
                            continue
                        distance = float(np.linalg.norm(coords[i] - coords[j]))
                        tail = empirical_tail_probability(values, distance)
                        anomalies.append(-math.log10(tail))
                        supports.append(len(values))
                    aggregates = (
                        _aggregate(np.asarray(anomalies))
                        if anomalies
                        else {
                            name: math.nan
                            for name in (
                                "mean",
                                "median",
                                "maximum",
                                "second_largest",
                                "fraction_95",
                                "count_95",
                                "fraction_99",
                                "count_99",
                            )
                        }
                    )
                    for aggregation, score in aggregates.items():
                        output.append(
                            {
                                **base,
                                "model": f"empirical_{strategy}_{aggregation}",
                                "score": score,
                                "n_context": len(anomalies),
                                "min_reference_n": min(supports, default=0),
                                "score_status": (
                                    "ok"
                                    if np.isfinite(score)
                                    else "insufficient_reference"
                                ),
                            }
                        )
    result = pd.DataFrame(output)
    result = result.rename(
        columns={
            "is_corrupted": "ground_truth_corrupted",
            "injected_displacement_um": "ground_truth_displacement",
        }
    )
    return result


def fixed_fpr_threshold(clean_scores: Sequence[float], target_fpr: float) -> float:
    """Choose a conservative empirical upper-tail threshold (calls use ``>``)."""
    values = np.sort(np.asarray(clean_scores, dtype=float))
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan
    allowed = int(math.floor(target_fpr * len(values)))
    return (
        float(values[max(0, len(values) - allowed - 1)])
        if allowed
        else float(values[-1])
    )


def _rates(group: pd.DataFrame, threshold: float) -> dict[str, float]:
    valid = group[np.isfinite(group["score"])]
    labels = valid["ground_truth_corrupted"].astype(bool).to_numpy()
    calls = valid["score"].to_numpy() > threshold
    tp, fp = int(np.sum(calls & labels)), int(np.sum(calls & ~labels))
    n_pos, n_neg = int(labels.sum()), int((~labels).sum())
    selected = valid["selected_for_corruption"].astype(bool).to_numpy()
    return {
        "true_positive_rate": tp / n_pos if n_pos else math.nan,
        "false_positive_rate": fp / n_neg if n_neg else math.nan,
        "precision": tp / (tp + fp) if tp + fp else math.nan,
        "n_corrupted": n_pos,
        "n_genuine": n_neg,
        "n_called": int(calls.sum()),
        "selected_target_call_rate": (
            float(np.mean(calls[selected])) if selected.any() else math.nan
        ),
        "n_selected_targets": int(selected.sum()),
    }


def evaluate(
    scores: pd.DataFrame, calibration: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate rankings and thresholds learned from independent replicates."""
    keys = ["detection_efficiency", "model", "context"]
    thresholds = {}
    for key, group in calibration.groupby(keys, sort=False):
        for fpr in FIXED_FPRS:
            thresholds[(*key, fpr)] = fixed_fpr_threshold(group["score"], fpr)
    fixed, overall = [], []
    condition_keys = [
        "condition",
        "replicate",
        "seed",
        "detection_efficiency",
        "displacement",
        "model",
        "context",
    ]
    for key, group in scores.groupby(condition_keys, sort=False):
        labels = group["ground_truth_corrupted"].astype(bool).to_numpy()
        valid = np.isfinite(group["score"].to_numpy())
        auc = (
            roc_auc_score(labels[valid], group.loc[valid, "score"])
            if valid.sum() and len(np.unique(labels[valid])) == 2
            else math.nan
        )
        ap = (
            average_precision_score(labels[valid], group.loc[valid, "score"])
            if valid.sum() and labels[valid].any()
            else math.nan
        )
        overall.append(
            dict(
                zip(condition_keys, key),
                roc_auc=auc,
                precision_recall_auc=ap,
                n_scored=int(valid.sum()),
            )
        )
        for fpr in FIXED_FPRS:
            threshold = thresholds.get((key[3], key[-2], key[-1], fpr), math.nan)
            fixed.append(
                {
                    **dict(zip(condition_keys, key)),
                    "target_fpr": fpr,
                    "threshold": threshold,
                    **_rates(group, threshold),
                }
            )
    return pd.DataFrame(overall), pd.DataFrame(fixed)


def _write(frame: pd.DataFrame, path: Path) -> None:
    safe = frame.copy()
    for column in safe.select_dtypes(include="object"):
        safe[column] = safe[column].fillna("").astype(str)
    Table.from_pandas(safe).write(path, format="ascii.ecsv", overwrite=True)


def plot_outputs(
    scores: pd.DataFrame, fixed: pd.DataFrame, output: Path
) -> pd.DataFrame:
    primary = fixed[np.isclose(fixed["target_fpr"], 0.01)].copy()
    primary["detection_rate"] = primary["true_positive_rate"]
    zero_displacement = np.isclose(primary["displacement"], 0)
    primary.loc[zero_displacement, "detection_rate"] = primary.loc[
        zero_displacement, "selected_target_call_rate"
    ]
    summary_keys = ["detection_efficiency", "displacement", "model", "context"]
    summary = (
        primary.groupby(summary_keys, dropna=False)["detection_rate"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    summary["sem"] = summary["std"] / np.sqrt(summary["count"])
    # Keep figures readable: the residual baseline plus multi-relationship means.
    shown = summary[
        summary["model"].isin(
            [
                "spatial_residual",
                "empirical_separation_mean",
                "empirical_barcode_pair_mean",
            ]
        )
    ]
    for efficiency, data in shown.groupby("detection_efficiency"):
        fig, ax = plt.subplots(figsize=(8, 5))
        for (model, context), line in data.groupby(["model", "context"]):
            line = line.sort_values("displacement")
            ax.errorbar(
                line["displacement"],
                line["mean"],
                yerr=line["sem"],
                marker="o",
                label=f"{model}/{context}",
            )
        ax.set(
            xlabel="Injected displacement (µm)",
            ylabel="Sensitivity at 1% target FPR",
            ylim=(-0.02, 1.02),
        )
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(output / f"sensitivity_efficiency_{efficiency:g}.png", dpi=160)
        plt.close(fig)
    context_order = {name: index for index, name in enumerate(CONTEXTS)}
    context_plot = shown[shown["displacement"] == shown["displacement"].max()].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    for model, data in context_plot.groupby("model"):
        values = data.groupby("context")["mean"].mean()
        ordered = sorted(values.index, key=context_order.get)
        ax.plot(ordered, values[ordered], marker="o", label=model)
    ax.set(xlabel="Context", ylabel="Mean sensitivity at largest displacement")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output / "performance_by_context.png", dpi=160)
    plt.close(fig)
    representative = scores[
        (scores["model"] == "empirical_separation_mean") & (scores["context"] == "k3")
    ]
    fig, ax = plt.subplots(figsize=(8, 5))
    displacements = sorted(representative["displacement"].unique())
    if displacements:
        positive = [value for value in displacements if value > 0]
        categories = [("clean", representative["displacement"] == displacements[0])]
        if positive:
            categories.extend(
                [
                    (
                        f"corrupted {positive[0]:g} µm",
                        (representative["displacement"] == positive[0])
                        & representative["ground_truth_corrupted"].astype(bool),
                    ),
                    (
                        f"corrupted {positive[-1]:g} µm",
                        (representative["displacement"] == positive[-1])
                        & representative["ground_truth_corrupted"].astype(bool),
                    ),
                ]
            )
        for label, selector in categories:
            values = representative.loc[selector, "score"].dropna()
            ax.hist(
                values,
                bins=50,
                density=True,
                histtype="step",
                label=label,
            )
    ax.set(xlabel="Anomaly score", ylabel="Density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "score_distributions.png", dpi=160)
    plt.close(fig)
    zero = scores[
        (scores["displacement"] == 0)
        & (scores["model"] == "empirical_separation_mean")
        & (scores["context"] == "k3")
    ]
    fig, ax = plt.subplots(figsize=(7, 5))
    for selected, label in (
        (False, "ordinary clean"),
        (True, "selected zero-displacement"),
    ):
        values = zero.loc[
            zero["selected_for_corruption"].astype(bool) == selected, "score"
        ].dropna()
        ax.hist(values, bins=50, density=True, histtype="step", label=label)
    ax.set(xlabel="Anomaly score", ylabel="Density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "zero_displacement_sanity.png", dpi=160)
    plt.close(fig)
    return summary


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-reference-observations", type=int, default=20)
    return parser.parse_args(argv)


def run(root: Path, output: Path, minimum_observations: int = 20) -> None:
    manifest_path = root / "sweep_manifest.yaml"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Benchmark manifest not found: {manifest_path}")
    manifest = yaml.safe_load(manifest_path.read_text())
    if manifest.get("dry_run"):
        raise ValueError("Cannot analyze a dry-run benchmark")
    conditions = manifest.get("conditions")
    simulations = manifest.get("simulations")
    if not conditions or not simulations:
        raise ValueError("Manifest must contain non-empty simulations and conditions")
    output.mkdir(parents=True, exist_ok=True)
    baselines = {
        str(item["id"]): load_baseline(root / "simulations" / str(item["id"]), item)
        for item in simulations
    }
    conditions_by_efficiency: dict[float, list[Mapping]] = defaultdict(list)
    for item in conditions:
        conditions_by_efficiency[float(item["detection_efficiency"])].append(item)
    all_scores, all_calibration, summaries = [], [], []
    for efficiency, condition_items in conditions_by_efficiency.items():
        seeds = sorted({int(item["seed"]) for item in condition_items})
        if len(seeds) < 3:
            raise ValueError(
                f"Efficiency {efficiency} needs at least 3 replicates for nested out-of-sample calibration"
            )
        efficiency_baselines = pd.concat(
            [
                frame
                for frame in baselines.values()
                if np.isclose(frame["detection_efficiency"].iloc[0], efficiency)
            ],
            ignore_index=True,
        )
        for evaluation_seed in seeds:
            reference_rows = efficiency_baselines[
                efficiency_baselines["seed"] != evaluation_seed
            ]
            model = fit_reference(reference_rows, minimum_observations)
            fold = {
                "detection_efficiency": efficiency,
                "evaluation_seed": evaluation_seed,
                "reference_seeds": ",".join(
                    map(str, sorted(set(reference_rows["seed"])))
                ),
            }
            summaries.append(reference_summary(model, fold))
            calibration_parts = []
            for calibration_seed in seeds:
                if calibration_seed == evaluation_seed:
                    continue
                nested_reference = efficiency_baselines[
                    ~efficiency_baselines["seed"].isin(
                        [evaluation_seed, calibration_seed]
                    )
                ]
                nested_model = fit_reference(nested_reference, minimum_observations)
                calibration_frame = efficiency_baselines[
                    efficiency_baselines["seed"] == calibration_seed
                ].copy()
                calibration_frame["condition"] = "calibration"
                calibration_frame["replicate"] = calibration_seed
                calibration_frame["displacement"] = 0.0
                calibration_frame["selected_for_corruption"] = False
                calibration_frame["injected_displacement_um"] = 0.0
                calibration_parts.append(
                    score_observations(calibration_frame, nested_model)
                )
            calibration = pd.concat(calibration_parts, ignore_index=True)
            all_calibration.append(calibration.assign(evaluation_seed=evaluation_seed))
            for item in condition_items:
                if int(item["seed"]) != evaluation_seed:
                    continue
                observations = load_observations(root / str(item["directory"]), item)
                all_scores.append(score_observations(observations, model))
    scores = pd.concat(all_scores, ignore_index=True)
    calibration_all = pd.concat(all_calibration, ignore_index=True)
    # Each evaluation fold needs its own independently calibrated thresholds.
    performances, fixed_tables = [], []
    for seed, fold_scores in scores.groupby("seed", sort=False):
        fold_calibration = calibration_all[calibration_all["evaluation_seed"] == seed]
        overall, fixed = evaluate(fold_scores, fold_calibration)
        performances.append(overall)
        fixed_tables.append(fixed)
    performance, fixed = pd.concat(performances, ignore_index=True), pd.concat(
        fixed_tables, ignore_index=True
    )
    plot_summary = plot_outputs(scores, fixed, output)
    _write(scores, output / "localization_scores.ecsv")
    _write(performance, output / "model_performance.ecsv")
    _write(fixed, output / "performance_at_fixed_fpr.ecsv")
    _write(
        pd.concat(summaries, ignore_index=True), output / "reference_model_summary.ecsv"
    )
    _write(plot_summary, output / "plot_aggregate_values.ecsv")


def main() -> None:
    args = parse_arguments()
    run(
        args.benchmark_root.resolve(),
        args.output_dir.resolve(),
        args.minimum_reference_observations,
    )


if __name__ == "__main__":
    main()
