#!/usr/bin/env python3
"""Evaluate whole-trace source-attribution models on benchmark_curator_v1.

This is an exploratory, bounded-memory benchmark analysis, not a production
curator.  Genomic order is always supplied explicitly; the opt-in uniform
barcode mapping exists only for the historical homogeneous simulation.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

FIXED_FPRS = (0.001, 0.01, 0.05)
STRONG_EDGE = -math.log10(0.05)
ML_FEATURES = (
    "k1_separation_mean",
    "k2_separation_mean",
    "all_separation_mean",
    "loo_separation_top3_mean",
    "edge_separation_maximum",
    "edge_separation_second_largest",
    "edge_separation_mean",
    "edge_separation_fraction_strong",
    "edge_separation_incident_vs_rest",
    "left_anomaly",
    "right_anomaly",
    "both_flanks_available",
    "left_genomic_gap",
    "right_genomic_gap",
    "bridge_separation",
    "n_detected_loci",
)
FORBIDDEN_FEATURES = {
    "ground_truth_corrupted",
    "selected_for_corruption",
    "ground_truth_displacement",
    "injected_displacement_um",
    "displacement",
    "condition",
    "seed",
    "simulation_id",
}
_START = time.monotonic()


@lru_cache(maxsize=1)
def _cycle1():
    """Load cycle 1 by path without turning exploratory scripts into a package."""
    path = Path(__file__).with_name("analyze_trace_curator_benchmark.py")
    spec = importlib.util.spec_from_file_location("trace_curator_cycle1", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _progress(message: str) -> None:
    print(f"[cycle 2 {time.monotonic() - _START:8.1f}s] {message}", file=sys.stderr, flush=True)


def _read(path: Path) -> pd.DataFrame:
    return Table.read(path, format="ascii.ecsv").to_pandas()


def _write(frame: pd.DataFrame, path: Path) -> None:
    safe = frame.copy()
    for column in safe.select_dtypes(include="object"):
        safe[column] = safe[column].fillna("").astype(str)
    Table.from_pandas(safe).write(path, format="ascii.ecsv", overwrite=True)


def load_genomic_coordinates(
    path: Path | None, barcodes: Iterable, assume_uniform_barcode_spacing: bool
) -> pd.DataFrame:
    """Return a validated barcode-to-coordinate mapping.

    Uniform barcode positions are deliberately unavailable unless explicitly
    enabled for the historical simulation benchmark.
    """
    if path is None:
        if not assume_uniform_barcode_spacing:
            raise ValueError(
                "Provide --genomic-coordinates; for benchmark_curator_v1 only, "
                "explicitly use --assume-uniform-barcode-spacing"
            )
        mapping = pd.DataFrame(
            {"Barcode": sorted(set(barcodes), key=lambda value: float(value))}
        )
        mapping["Genomic_Position"] = mapping["Barcode"].astype(float)
        return mapping
    mapping = _read(path) if path.suffix.lower() == ".ecsv" else pd.read_csv(path)
    mapping = mapping.rename(columns={"Barcode #": "Barcode"})
    required = {"Barcode", "Genomic_Position"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"Coordinate mapping requires columns {sorted(required)}")
    if mapping["Barcode"].duplicated().any() or mapping["Genomic_Position"].duplicated().any():
        raise ValueError("Barcode and Genomic_Position values must each be unique")
    if not np.isfinite(mapping["Genomic_Position"].astype(float)).all():
        raise ValueError("Genomic_Position values must be finite")
    missing = set(barcodes).difference(set(mapping["Barcode"]))
    if missing:
        raise ValueError(f"Coordinate mapping is missing barcodes: {sorted(missing)}")
    return mapping[["Barcode", "Genomic_Position"]].copy()


def attach_genomic_coordinates(frame: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    result = frame.drop(columns=["Genomic_Position"], errors="ignore").merge(
        mapping, on="Barcode", how="left", validate="many_to_one"
    )
    if result["Genomic_Position"].isna().any():
        raise ValueError("Every localization must have an explicit genomic coordinate")
    return result


def select_genomic_context(
    genomic_positions: Sequence[float], target_index: int, context: str
) -> list[int]:
    """Select nearest detected loci on each side in genomic-coordinate order."""
    if context not in {"k1", "k2", "k3", "all"}:
        raise ValueError(f"Unknown context: {context}")
    positions = np.asarray(genomic_positions, dtype=float)
    target = positions[target_index]
    candidates = np.flatnonzero(np.arange(len(positions)) != target_index)
    if context == "all":
        return candidates.tolist()
    k = int(context[1:])
    left = candidates[positions[candidates] < target]
    right = candidates[positions[candidates] > target]
    left = left[np.argsort(target - positions[left])][:k]
    right = right[np.argsort(positions[right] - target)][:k]
    return np.concatenate([left, right]).tolist()


def nearest_genomic_flanks(positions: Sequence[float], target_index: int) -> tuple[int | None, int | None]:
    selected = select_genomic_context(positions, target_index, "k1")
    target = float(positions[target_index])
    left = [index for index in selected if positions[index] < target]
    right = [index for index in selected if positions[index] > target]
    return (left[0] if left else None, right[0] if right else None)


def genomic_separation(left: float, right: float) -> float:
    return abs(float(right) - float(left))


def _coordinate_key(value: float) -> float:
    return float(value)


@dataclass
class ReferenceModel:
    separation: dict[float, np.ndarray]
    barcode_pair: dict[tuple[object, object], np.ndarray]
    bridge: dict[tuple[float, float], np.ndarray]
    minimum_observations: int = 20

    def edge_values(self, strategy: str, barcode_a, barcode_b, position_a, position_b):
        key = (
            _coordinate_key(genomic_separation(position_a, position_b))
            if strategy == "separation"
            else tuple(sorted((barcode_a, barcode_b), key=str))
        )
        values = getattr(self, strategy).get(key)
        return values if values is not None and len(values) >= self.minimum_observations else None


def _trace_key(frame: pd.DataFrame):
    return ["simulation_id", "Trace_ID"] if "simulation_id" in frame else "Trace_ID"


def fit_reference(frame: pd.DataFrame, minimum_observations: int = 20) -> ReferenceModel:
    """Fit clean edge and joint bridge distributions using genomic coordinates."""
    clean = frame.loc[~frame["is_corrupted"].astype(bool)]
    separations, pairs, bridges = defaultdict(list), defaultdict(list), defaultdict(list)
    for _, trace in clean.groupby(_trace_key(clean), sort=False):
        trace = trace.sort_values("Genomic_Position")
        coords = trace[["x", "y", "z"]].to_numpy(float)
        positions = trace["Genomic_Position"].to_numpy(float)
        barcodes = trace["Barcode"].to_numpy()
        for i in range(len(trace)):
            for j in range(i + 1, len(trace)):
                distance = float(np.linalg.norm(coords[i] - coords[j]))
                separations[_coordinate_key(positions[j] - positions[i])].append(distance)
                pairs[tuple(sorted((barcodes[i], barcodes[j]), key=str))].append(distance)
        for i in range(1, len(trace) - 1):
            gaps = (_coordinate_key(positions[i] - positions[i - 1]), _coordinate_key(positions[i + 1] - positions[i]))
            bridges[gaps].append(
                [np.linalg.norm(coords[i] - coords[i - 1]), np.linalg.norm(coords[i + 1] - coords[i]), np.linalg.norm(coords[i + 1] - coords[i - 1])]
            )
    return ReferenceModel(
        {key: np.asarray(values) for key, values in separations.items()},
        {key: np.asarray(values) for key, values in pairs.items()},
        {key: np.asarray(values) for key, values in bridges.items()},
        minimum_observations,
    )


def empirical_anomaly(values: Sequence[float], value: float) -> float:
    sample = np.sort(np.asarray(values, float))
    below = np.searchsorted(sample, value, side="left")
    at_or_below = np.searchsorted(sample, value, side="right")
    lower = (at_or_below + 1) / (len(sample) + 1)
    upper = (len(sample) - below + 1) / (len(sample) + 1)
    return -math.log10(min(1.0, 2 * min(lower, upper)))


def trace_cost(edge_scores: Sequence[float], variant: str) -> float:
    values = np.asarray(edge_scores, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return math.nan
    if variant == "mean":
        return float(values.mean())
    if variant == "trimmed_mean":
        ordered = np.sort(values)
        trim = int(math.floor(0.1 * len(ordered)))
        return float(ordered[trim : len(ordered) - trim].mean()) if trim else float(ordered.mean())
    if variant == "top3_mean":
        return float(np.sort(values)[-min(3, len(values)) :].mean())
    if variant == "strong_fraction":
        return float(np.mean(values >= STRONG_EDGE))
    raise ValueError(f"Unknown trace cost: {variant}")


def leave_one_out_improvements(
    n_nodes: int, edges: Sequence[tuple[int, int, float]], variant: str
) -> np.ndarray:
    """Calculate C(T) - C(T without i) exactly from trace edges."""
    before = trace_cost([edge[2] for edge in edges], variant)
    result = []
    for node in range(n_nodes):
        after = trace_cost([score for left, right, score in edges if node not in (left, right)], variant)
        result.append(before - after if np.isfinite(after) else math.nan)
    return np.asarray(result)


def edge_concentration_features(n_nodes: int, edges: Sequence[tuple[int, int, float]]) -> pd.DataFrame:
    rows = []
    all_scores = np.asarray([score for _, _, score in edges], float)
    for node in range(n_nodes):
        incident = np.asarray([score for left, right, score in edges if node in (left, right)], float)
        elsewhere = np.asarray([score for left, right, score in edges if node not in (left, right)], float)
        ordered = np.sort(incident)
        rows.append(
            {
                "mean": float(np.mean(incident)) if len(incident) else math.nan,
                "maximum": float(ordered[-1]) if len(ordered) else math.nan,
                "second_largest": float(ordered[-2]) if len(ordered) >= 2 else math.nan,
                "fraction_strong": float(np.mean(incident >= STRONG_EDGE)) if len(incident) else math.nan,
                "count_strong": int(np.sum(incident >= STRONG_EDGE)),
                "incident_vs_rest": (float(np.mean(incident) - np.mean(elsewhere)) if len(incident) and len(elsewhere) else math.nan),
                "concentration": (float(np.mean(ordered[-min(3, len(ordered)) :]) - np.mean(elsewhere)) if len(ordered) and len(elsewhere) else math.nan),
                "n_incident": len(incident),
                "n_edges": len(all_scores),
            }
        )
    return pd.DataFrame(rows)


def bridge_score(trace: pd.DataFrame, index: int, model: ReferenceModel) -> dict:
    """Score a localization jointly against its true nearest genomic flanks."""
    positions = trace["Genomic_Position"].to_numpy(float)
    left, right = nearest_genomic_flanks(positions, index)
    base = {"both_flanks_available": left is not None and right is not None}
    if left is None or right is None:
        return {**base, "score": math.nan, "status": "insufficient_bridge_context", "left_index": left, "right_index": right}
    coords = trace[["x", "y", "z"]].to_numpy(float)
    gaps = (_coordinate_key(positions[index] - positions[left]), _coordinate_key(positions[right] - positions[index]))
    values = model.bridge.get(gaps)
    result = {**base, "left_index": left, "right_index": right, "left_genomic_gap": gaps[0], "right_genomic_gap": gaps[1]}
    if values is None or len(values) < model.minimum_observations:
        return {**result, "score": math.nan, "status": "insufficient_reference"}
    observed = np.asarray([np.linalg.norm(coords[index] - coords[left]), np.linalg.norm(coords[right] - coords[index]), np.linalg.norm(coords[right] - coords[left])])
    # Sum of squared robust standardized joint residuals retains covariance-free
    # interpretability and evaluates the three bridge distances together.
    center = np.median(values, axis=0)
    scale = np.maximum(np.median(np.abs(values - center), axis=0) * 1.4826, 1e-6)
    return {**result, "score": float(np.mean(((observed - center) / scale) ** 2)), "status": "ok"}


def _edge_tables(trace: pd.DataFrame, model: ReferenceModel, strategy: str):
    coords = trace[["x", "y", "z"]].to_numpy(float)
    positions = trace["Genomic_Position"].to_numpy(float)
    barcodes = trace["Barcode"].to_numpy()
    edges = []
    lookup = {}
    for i in range(len(trace)):
        for j in range(i + 1, len(trace)):
            values = model.edge_values(strategy, barcodes[i], barcodes[j], positions[i], positions[j])
            if values is None:
                continue
            score = empirical_anomaly(values, float(np.linalg.norm(coords[i] - coords[j])))
            edges.append((i, j, score))
            lookup[(i, j)] = score
    return edges, lookup


def score_trace(trace: pd.DataFrame, model: ReferenceModel) -> pd.DataFrame:
    """Build feature-only scores; evaluation labels are never consulted."""
    trace = trace.sort_values("Genomic_Position").reset_index(drop=True)
    positions = trace["Genomic_Position"].to_numpy(float)
    output = trace[["Trace_ID", "Spot_ID", "Barcode", "Genomic_Position"]].copy()
    for strategy in ("separation", "barcode_pair"):
        edges, lookup = _edge_tables(trace, model, strategy)
        concentration = edge_concentration_features(len(trace), edges)
        for context in ("k1", "k2", "all"):
            scores = []
            for i in range(len(trace)):
                selected = select_genomic_context(positions, i, context)
                values = [lookup.get(tuple(sorted((i, j)))) for j in selected]
                finite = [value for value in values if value is not None]
                scores.append(float(np.mean(finite)) if finite else math.nan)
            output[f"{context}_{strategy}_mean"] = scores
        for variant in ("mean", "trimmed_mean", "top3_mean", "strong_fraction"):
            output[f"loo_{strategy}_{variant}"] = leave_one_out_improvements(len(trace), edges, variant)
        for feature in ("mean", "maximum", "second_largest", "fraction_strong", "count_strong", "incident_vs_rest", "concentration"):
            output[f"edge_{strategy}_{feature}"] = concentration[feature].to_numpy()
    bridge_rows = [bridge_score(trace, index, model) for index in range(len(trace))]
    output["bridge_separation"] = [row["score"] for row in bridge_rows]
    output["bridge_status"] = [row["status"] for row in bridge_rows]
    output["both_flanks_available"] = [row["both_flanks_available"] for row in bridge_rows]
    output["left_genomic_gap"] = [row.get("left_genomic_gap", math.nan) for row in bridge_rows]
    output["right_genomic_gap"] = [row.get("right_genomic_gap", math.nan) for row in bridge_rows]
    left_anomaly, right_anomaly = [], []
    _, separation_lookup = _edge_tables(trace, model, "separation")
    for index in range(len(trace)):
        left, right = nearest_genomic_flanks(positions, index)
        left_anomaly.append(separation_lookup.get(tuple(sorted((index, left)))) if left is not None else math.nan)
        right_anomaly.append(separation_lookup.get(tuple(sorted((index, right)))) if right is not None else math.nan)
    output["left_anomaly"], output["right_anomaly"] = left_anomaly, right_anomaly
    output["n_detected_loci"] = len(trace)
    return output


def score_observations(frame: pd.DataFrame, model: ReferenceModel) -> pd.DataFrame:
    parts = []
    feature_frame = frame[[column for column in frame.columns if column not in FORBIDDEN_FEATURES and column not in {"is_corrupted", "selected_for_corruption", "injected_displacement_um"}] + (["simulation_id"] if "simulation_id" in frame else [])]
    # score_trace consumes coordinates/identity only; labels are reattached below.
    for _, trace in feature_frame.groupby(_trace_key(feature_frame), sort=False):
        scored = score_trace(trace, model)
        scored["simulation_id"] = trace["simulation_id"].iloc[0] if "simulation_id" in trace else ""
        parts.append(scored)
    result = pd.concat(parts, ignore_index=True)
    keys = ["simulation_id", "Spot_ID"] if "simulation_id" in frame else ["Spot_ID"]
    evaluation_columns = keys + [column for column in frame.columns if column not in result.columns and column not in {"x", "y", "z", "Genomic_Position"}]
    return result.merge(frame[evaluation_columns], on=keys, how="left", validate="one_to_one")


def long_scores(features: pd.DataFrame, include_ml: bool = True) -> pd.DataFrame:
    models = {
        "baseline_separation_k1": "k1_separation_mean",
        "baseline_barcode_pair_k1": "k1_barcode_pair_mean",
        "loo_separation_mean": "loo_separation_mean",
        "loo_separation_trimmed_mean": "loo_separation_trimmed_mean",
        "loo_separation_top3_mean": "loo_separation_top3_mean",
        "loo_separation_strong_fraction": "loo_separation_strong_fraction",
        "loo_barcode_pair_top3_mean": "loo_barcode_pair_top3_mean",
        "edge_separation_concentration": "edge_separation_concentration",
        "edge_barcode_pair_concentration": "edge_barcode_pair_concentration",
        "bridge_separation": "bridge_separation",
    }
    if include_ml and "ml_probability" in features:
        models["ml_logistic_regression"] = "ml_probability"
    id_columns = [column for column in ("condition", "replicate", "seed", "simulation_id", "Trace_ID", "Spot_ID", "Barcode", "detection_efficiency", "displacement", "is_corrupted", "selected_for_corruption") if column in features]
    parts = []
    for model_name, column in models.items():
        part = features[id_columns].copy()
        part["model"], part["score"] = model_name, features[column]
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def ml_training_matrix(frame: pd.DataFrame) -> np.ndarray:
    if FORBIDDEN_FEATURES.intersection(ML_FEATURES):
        raise AssertionError("Evaluation labels or metadata entered ML_FEATURES")
    return frame.loc[:, ML_FEATURES].astype(float).to_numpy()


def fit_ml(training: pd.DataFrame):
    pipeline = make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler(), LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=0))
    pipeline.fit(ml_training_matrix(training), training["is_corrupted"].astype(int))
    return pipeline


def assert_seed_split(training: pd.DataFrame, evaluation_seed: int) -> None:
    if evaluation_seed in set(training["seed"].astype(int)):
        raise ValueError("Evaluation seed leaked into ML training")


def attribution_metrics(scores: pd.DataFrame, threshold: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return per-trace attribution, collateral, and rank-distribution rows."""
    attribution, collateral, ranks = [], [], []
    group_columns = [column for column in ("condition", "replicate", "seed", "detection_efficiency", "displacement", "model", "simulation_id", "Trace_ID") if column in scores]
    for key, group in scores.groupby(group_columns, sort=False):
        metadata = dict(zip(group_columns, key if isinstance(key, tuple) else (key,)))
        valid = group[np.isfinite(group["score"])].copy()
        bad = valid[valid["is_corrupted"].astype(bool)]
        if len(bad) != 1:
            continue
        # Stable ranking gives deterministic handling; ties receive minimum rank.
        bad_score = float(bad["score"].iloc[0])
        rank = int(1 + np.sum(valid["score"].to_numpy() > bad_score))
        good = valid[~valid["is_corrupted"].astype(bool)]
        margin = bad_score - float(good["score"].max()) if len(good) else math.nan
        calls = good["score"].to_numpy() > threshold
        attribution.append({**metadata, "rank": rank, "top1": rank == 1, "top2": rank <= 2, "score_margin": margin, "n_scored_in_trace": len(valid)})
        collateral.append({**metadata, "threshold": threshold, "collateral_calls": int(calls.sum()), "zero_collateral_calls": not calls.any()})
        ranks.append({**metadata, "corrupted_rank": rank})
    return pd.DataFrame(attribution), pd.DataFrame(collateral), pd.DataFrame(ranks)


def fixed_fpr_threshold(values: Sequence[float], fpr: float) -> float:
    values = np.sort(np.asarray(values, float)); values = values[np.isfinite(values)]
    if not len(values): return math.nan
    allowed = int(math.floor(fpr * len(values)))
    return float(values[max(0, len(values) - allowed - 1)]) if allowed else float(values[-1])


def performance(scores: pd.DataFrame, thresholds: Mapping[tuple[str, float], float]):
    overall, fixed = [], []
    keys = ["condition", "replicate", "seed", "detection_efficiency", "displacement", "model"]
    for key, group in scores.groupby(keys, sort=False):
        metadata = dict(zip(keys, key)); valid = group[np.isfinite(group.score)]
        labels = valid.is_corrupted.astype(bool).to_numpy(); values = valid.score.to_numpy()
        overall.append({**metadata, "roc_auc": roc_auc_score(labels, values) if len(np.unique(labels)) == 2 else math.nan, "precision_recall_auc": average_precision_score(labels, values) if labels.any() else math.nan, "n_scored": len(valid)})
        for fpr in FIXED_FPRS:
            threshold = thresholds.get((str(key[-1]), fpr), math.nan); calls = values > threshold
            tp, fp = np.sum(calls & labels), np.sum(calls & ~labels)
            fixed.append({**metadata, "target_fpr": fpr, "threshold": threshold, "sensitivity": tp / labels.sum() if labels.sum() else math.nan, "observed_fpr": fp / (~labels).sum() if (~labels).sum() else math.nan, "precision": tp / calls.sum() if calls.sum() else math.nan, "n_corrupted": int(labels.sum()), "n_genuine": int((~labels).sum())})
    return pd.DataFrame(overall), pd.DataFrame(fixed)


def _aggregate_replicates(frame: pd.DataFrame, value_columns: Sequence[str]) -> pd.DataFrame:
    keys = [column for column in ("model", "detection_efficiency", "displacement") if column in frame]
    rows = []
    replicate_keys = [*keys, *[column for column in ("seed", "replicate") if column in frame]]
    replicate_values = frame.groupby(replicate_keys, sort=False)[list(value_columns)].mean().reset_index()
    for key, group in replicate_values.groupby(keys, sort=False):
        base = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        for value in value_columns:
            data = group[value].dropna().astype(float)
            rows.append({**base, "metric": value, "mean": data.mean(), "sd": data.std(ddof=1), "sem": data.sem(ddof=1), "n_replicates": len(data)})
    return pd.DataFrame(rows)


def plot_results(fixed: pd.DataFrame, attribution: pd.DataFrame, collateral: pd.DataFrame, output: Path) -> None:
    selected = fixed[np.isclose(fixed.target_fpr, 0.01)]
    plots = [(selected, "sensitivity", "Sensitivity", "cycle2_sensitivity.png"), (attribution, "top1", "Top-1 attribution", "cycle2_top1_attribution.png"), (attribution, "top2", "Top-2 attribution", "cycle2_top2_attribution.png"), (collateral, "collateral_calls", "Mean collateral calls", "cycle2_collateral_calls.png"), (attribution, "score_margin", "Median score margin", "cycle2_score_margin.png")]
    for frame, value, ylabel, filename in plots:
        if frame.empty: continue
        figure, axis = plt.subplots(figsize=(9, 5))
        replicate = frame.groupby(["model", "detection_efficiency", "displacement", "seed"])[value].mean().reset_index()
        summary = replicate.groupby(["model", "detection_efficiency", "displacement"])[value].agg(["mean", "sem"]).reset_index()
        for (model, efficiency), group in summary.groupby(["model", "detection_efficiency"]):
            axis.errorbar(group.displacement, group["mean"], yerr=group["sem"], marker="o", label=f"{model}; eff={efficiency:g}")
        axis.set(xlabel="Displacement (µm)", ylabel=f"{ylabel} (mean ± SEM across seeds)"); axis.legend(fontsize=6, ncol=2); figure.tight_layout(); figure.savefig(output / filename, dpi=150); plt.close(figure)
    representative = attribution[np.isclose(attribution.displacement, 0.2)]
    if not representative.empty:
        figure, axis = plt.subplots(figsize=(8, 5))
        for model, group in representative.groupby("model"):
            counts = group["rank"].value_counts(normalize=True).sort_index()
            axis.step(counts.index, counts.values, where="mid", label=model)
        axis.set(xlabel="Rank of corrupted localization", ylabel="Fraction of traces", title="Rank distribution at 0.2 µm")
        axis.legend(fontsize=6, ncol=2); figure.tight_layout(); figure.savefig(output / "cycle2_rank_distribution.png", dpi=150); plt.close(figure)


def run(root: Path, output: Path, genomic_coordinates: Path | None = None, assume_uniform_barcode_spacing: bool = False, minimum_observations: int = 20, write_localization_scores: bool = False) -> None:
    """Run seed-held-out analysis while retaining at most one condition's scores."""
    cycle1 = _cycle1(); manifest = yaml.safe_load((root / "sweep_manifest.yaml").read_text())
    simulations, conditions = manifest["simulations"], manifest["conditions"]
    # Inspect baseline headers to establish the explicit mapping once.
    raw_baselines = {str(item["id"]): cycle1.load_baseline(root / "simulations" / str(item["id"]), item) for item in simulations}
    all_barcodes = set().union(*(set(frame.Barcode) for frame in raw_baselines.values()))
    mapping = load_genomic_coordinates(genomic_coordinates, all_barcodes, assume_uniform_barcode_spacing)
    baselines = {key: attach_genomic_coordinates(frame, mapping) for key, frame in raw_baselines.items()}
    output.mkdir(parents=True, exist_ok=True)
    overall_parts, fixed_parts, attribution_parts, collateral_parts, rank_parts, references, ml_importance = [], [], [], [], [], [], []
    detail_writer = cycle1._EcsvChunkWriter(output / "cycle2_localization_scores.ecsv") if write_localization_scores else None
    by_efficiency = defaultdict(list)
    for item in conditions: by_efficiency[float(item["detection_efficiency"])].append(item)
    for efficiency, items in by_efficiency.items():
        seeds = sorted({int(item["seed"]) for item in items})
        efficiency_baselines = pd.concat([frame for frame in baselines.values() if np.isclose(frame.detection_efficiency.iloc[0], efficiency)], ignore_index=True)
        for evaluation_seed in seeds:
            _progress(f"efficiency={efficiency:g}, held-out seed={evaluation_seed}")
            reference_rows = efficiency_baselines[efficiency_baselines.seed != evaluation_seed]
            model = fit_reference(reference_rows, minimum_observations)
            for strategy in ("separation", "barcode_pair", "bridge"):
                for key, values in getattr(model, strategy).items(): references.append({"detection_efficiency": efficiency, "evaluation_seed": evaluation_seed, "reference_strategy": strategy, "reference_key": str(key), "n_observations": len(values), "sufficient": len(values) >= minimum_observations})
            def build_ml_training(excluded_seeds, feature_model):
                # Retain only the compact fixed ML matrix, label, and split key;
                # full condition score frames are discarded immediately.
                compact_parts = []
                for train_item in items:
                    if int(train_item["seed"]) in excluded_seeds: continue
                    train_observations = attach_genomic_coordinates(cycle1.load_observations(root / str(train_item["directory"]), train_item), mapping)
                    train_features = score_observations(train_observations, feature_model)
                    compact_parts.append(train_features[[*ML_FEATURES, "is_corrupted", "seed"]].copy())
                    del train_observations, train_features
                return pd.concat(compact_parts, ignore_index=True)

            training = build_ml_training({evaluation_seed}, model); assert_seed_split(training, evaluation_seed)
            learner = fit_ml(training)
            coefficients = learner.named_steps["logisticregression"].coef_[0]
            for index, coefficient in enumerate(coefficients[: len(ML_FEATURES)]): ml_importance.append({"detection_efficiency": efficiency, "evaluation_seed": evaluation_seed, "feature": ML_FEATURES[index], "coefficient": coefficient})
            # Nested cross-fit: each calibration seed is excluded from its
            # statistical reference and its ML learner, as is evaluation_seed.
            calibration_long_parts = []
            for calibration_seed in seeds:
                if calibration_seed == evaluation_seed: continue
                nested_rows = efficiency_baselines[~efficiency_baselines.seed.isin([evaluation_seed, calibration_seed])]
                nested_model = fit_reference(nested_rows, minimum_observations)
                nested_training = build_ml_training({evaluation_seed, calibration_seed}, nested_model)
                nested_learner = fit_ml(nested_training)
                calibration = efficiency_baselines[efficiency_baselines.seed == calibration_seed].copy()
                calibration["condition"] = "clean_calibration"; calibration["replicate"] = calibration_seed; calibration["displacement"] = 0.0; calibration["selected_for_corruption"] = False
                calibration_features = score_observations(calibration, nested_model)
                calibration_features["ml_probability"] = nested_learner.predict_proba(ml_training_matrix(calibration_features))[:, 1]
                calibration_long_parts.append(long_scores(calibration_features))
                del nested_rows, nested_model, nested_training, nested_learner, calibration, calibration_features
            calibration_long = pd.concat(calibration_long_parts, ignore_index=True)
            thresholds = {(name, fpr): fixed_fpr_threshold(group.score, fpr) for name, group in calibration_long.groupby("model") for fpr in FIXED_FPRS}
            del calibration_long, calibration_long_parts, training
            for item in [entry for entry in items if int(entry["seed"]) == evaluation_seed]:
                _progress(f"scoring {item['id']}")
                observations = attach_genomic_coordinates(cycle1.load_observations(root / str(item["directory"]), item), mapping)
                features = score_observations(observations, model); features["ml_probability"] = learner.predict_proba(ml_training_matrix(features))[:, 1]
                scores = long_scores(features); overall, fixed = performance(scores, thresholds); overall_parts.append(overall); fixed_parts.append(fixed)
                for model_name, group in scores.groupby("model"):
                    att, col, ranks = attribution_metrics(group, thresholds[(str(model_name), 0.01)]); attribution_parts.append(att); collateral_parts.append(col); rank_parts.append(ranks)
                if detail_writer is not None:
                    detail_writer.write(features)
                del observations, features, scores, overall, fixed
    if detail_writer is not None:
        detail_writer.finish()
    performance_table = pd.concat(overall_parts, ignore_index=True); fixed_table = pd.concat(fixed_parts, ignore_index=True)
    attribution_table = pd.concat(attribution_parts, ignore_index=True); collateral_table = pd.concat(collateral_parts, ignore_index=True); rank_table = pd.concat(rank_parts, ignore_index=True)
    _write(performance_table, output / "cycle2_model_performance.ecsv"); _write(fixed_table, output / "cycle2_fixed_fpr_performance.ecsv"); _write(attribution_table, output / "cycle2_attribution_performance.ecsv"); _write(collateral_table, output / "cycle2_collateral_calls.ecsv"); _write(rank_table, output / "cycle2_rank_distribution.ecsv"); _write(pd.DataFrame(references), output / "cycle2_reference_summary.ecsv"); _write(pd.DataFrame(ml_importance), output / "cycle2_ml_feature_importance.ecsv")
    aggregate = pd.concat([_aggregate_replicates(fixed_table[np.isclose(fixed_table.target_fpr, .01)], ["sensitivity", "observed_fpr", "precision"]), _aggregate_replicates(attribution_table, ["top1", "top2", "rank", "score_margin"]), _aggregate_replicates(collateral_table, ["collateral_calls", "zero_collateral_calls"])], ignore_index=True)
    _write(aggregate, output / "cycle2_aggregate_performance.ecsv"); plot_results(fixed_table, attribution_table, collateral_table, output)
    _progress(f"complete: {output}")


def parse_arguments(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--benchmark-root", required=True, type=Path); parser.add_argument("--output-dir", required=True, type=Path); parser.add_argument("--genomic-coordinates", type=Path); parser.add_argument("--assume-uniform-barcode-spacing", action="store_true", help="Simulation-only: use numeric barcode values as uniformly spaced polymer positions."); parser.add_argument("--minimum-reference-observations", type=int, default=20); parser.add_argument("--write-localization-scores", action="store_true"); return parser.parse_args(argv)


def main() -> None:
    args = parse_arguments(); run(args.benchmark_root.resolve(), args.output_dir.resolve(), args.genomic_coordinates.resolve() if args.genomic_coordinates else None, args.assume_uniform_barcode_spacing, args.minimum_reference_observations, args.write_localization_scores)


if __name__ == "__main__": main()
