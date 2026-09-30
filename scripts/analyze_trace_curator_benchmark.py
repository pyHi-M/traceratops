#!/usr/bin/env python3
"""Compare interpretable anomaly scores on a trace-curator benchmark.

This is an exploratory benchmark utility, not a production ``trace_curator``
command.  Reference models and operating-point thresholds are learned without
using the replicate being evaluated.  See
``docs/source/contribute/trace_curator_benchmark_analysis.md`` for details.
"""

from __future__ import annotations

import argparse
import io
import math
import os
import resource
import sys
import tempfile
import time
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

_START_TIME = time.monotonic()


def _progress(message: str) -> None:
    """Print a timestamped progress message immediately.

    Progress goes to stderr so redirected result output remains clean.  Explicit
    flushing is important when this script is run non-interactively, where the
    normal stream buffering previously made a long analysis appear to hang.
    """
    elapsed = time.monotonic() - _START_TIME
    try:
        resident_pages = int(Path("/proc/self/statm").read_text().split()[1])
        rss = resident_pages * os.sysconf("SC_PAGE_SIZE") / (1024**2)
        memory_label = "RSS"
    except (FileNotFoundError, IndexError, OSError, ValueError):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        memory_label = "peak RSS"
    print(
        f"[trace benchmark {elapsed:8.1f}s, {memory_label} {rss:,.0f} MiB] {message}",
        file=sys.stderr,
        flush=True,
    )


class _EcsvChunkWriter:
    """Incrementally write identically shaped frames to one ECSV file.

    Pairwise diagnostics are by far the largest output of this analysis.  The
    old implementation retained every diagnostic dictionary until the very end,
    causing memory use to grow for the entire run.  ECSV has one commented
    metadata header and a plain-text data section, so subsequent Astropy-rendered
    chunks can safely contribute only their data rows.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self.temporary_path = path.with_name(f".{path.name}.tmp")
        self.columns: list[str] | None = None
        self._has_rows = False
        self.temporary_path.unlink(missing_ok=True)

    @staticmethod
    def _render(frame: pd.DataFrame) -> list[str]:
        safe = frame.copy()
        for column in safe.select_dtypes(include="object"):
            safe[column] = safe[column].fillna("").astype(str)
        stream = io.StringIO()
        Table.from_pandas(safe).write(stream, format="ascii.ecsv")
        return stream.getvalue().splitlines(keepends=True)

    def write(self, frame: pd.DataFrame, chunk_rows: int = 50_000) -> None:
        if frame.empty:
            return
        columns = list(frame.columns)
        if self.columns is None:
            self.columns = columns
        elif columns != self.columns:
            raise ValueError("Pairwise diagnostic columns changed between chunks")
        for start in range(0, len(frame), chunk_rows):
            lines = self._render(frame.iloc[start : start + chunk_rows])
            mode = "a" if self._has_rows else "w"
            if self._has_rows:
                # The first non-comment line is the repeated column-name row.
                header_index = next(
                    index
                    for index, line in enumerate(lines)
                    if not line.startswith("#")
                )
                lines = lines[header_index + 1 :]
            with self.temporary_path.open(mode) as stream:
                stream.writelines(lines)
            self._has_rows = True

    def finish(self) -> None:
        if not self._has_rows:
            _write(pd.DataFrame(), self.temporary_path)
        os.replace(self.temporary_path, self.path)


class _CalibrationAccumulator:
    """Spool the score-only calibration data needed to derive thresholds."""

    def __init__(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="trace-benchmark-calibration-"
        )
        self._keys: dict[tuple[str, str], Path] = {}

    def add(self, scores: pd.DataFrame) -> None:
        for (model, context), group in scores.groupby(["model", "context"], sort=False):
            key = (str(model), str(context))
            path = self._keys.setdefault(
                key,
                Path(self._temporary_directory.name) / f"scores-{len(self._keys)}.bin",
            )
            values = group["score"].to_numpy(dtype=np.float64)
            with path.open("ab") as stream:
                values[np.isfinite(values)].tofile(stream)

    def thresholds(self) -> dict[tuple[str, str, float], float]:
        result = {}
        for (model, context), path in self._keys.items():
            values = np.fromfile(path, dtype=np.float64)
            for fpr in FIXED_FPRS:
                result[(model, context, fpr)] = fixed_fpr_threshold(values, fpr)
            del values
        return result

    def close(self) -> None:
        self._temporary_directory.cleanup()


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


def trace_group_key(frame: pd.DataFrame) -> str | list[str]:
    """Return a trace key that cannot merge independent simulations."""
    if "simulation_id" in frame.columns:
        return ["simulation_id", "Trace_ID"]
    return "Trace_ID"


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


def empirical_distance_diagnostics(
    values: Sequence[float], value: float
) -> tuple[float, float]:
    """Return directional empirical percentile and two-sided tail probability.

    The percentile uses the midpoint of ties and add-one smoothing. The tail
    probability uses add-one-smoothed inclusive lower and upper tails, so it is
    positive and at most one. Percentiles below/above 0.5 indicate unusually
    short/long distances, respectively.
    """
    sample = np.sort(np.asarray(values, dtype=float))
    if not len(sample) or not np.isfinite(value):
        return math.nan, math.nan
    below = np.searchsorted(sample, value, side="left")
    at_or_below = np.searchsorted(sample, value, side="right")
    percentile = (below + 0.5 * (at_or_below - below) + 0.5) / (len(sample) + 1)
    lower = (at_or_below + 1) / (len(sample) + 1)
    upper = (len(sample) - below + 1) / (len(sample) + 1)
    tail = min(1.0, 2.0 * min(lower, upper))
    return float(percentile), float(tail)


def empirical_tail_probability(values: Sequence[float], value: float) -> float:
    """Return the two-sided component of :func:`empirical_distance_diagnostics`."""
    return empirical_distance_diagnostics(values, value)[1]


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

    def residual_stats(
        self, separation: int, fallback_mode: str = "none"
    ) -> tuple[tuple[float, float, int] | None, bool]:
        """Return an exact residual bin, optionally using an explicit fallback."""
        if fallback_mode not in {"none", "nearest"}:
            raise ValueError("fallback_mode must be 'none' or 'nearest'")
        populated = {
            key: value
            for key, value in self.residual.items()
            if value[2] >= self.minimum_observations
        }
        if not populated:
            return None, False
        separation = abs(int(separation))
        if separation in populated:
            return populated[separation], False
        if fallback_mode == "none":
            return None, False
        key = min(populated, key=lambda candidate: abs(candidate - separation))
        return populated[key], True


def fit_reference(
    frame: pd.DataFrame, minimum_observations: int = 20
) -> ReferenceModel:
    """Fit pair distributions using explicitly clean rows only."""
    clean = clean_reference_rows(frame)
    by_sep: dict[int, list[float]] = defaultdict(list)
    by_pair: dict[tuple[int, int], list[float]] = defaultdict(list)
    for _, trace in clean.groupby(trace_group_key(clean), sort=False):
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


def score_observations(
    frame: pd.DataFrame,
    model: ReferenceModel,
    residual_fallback: str = "none",
    relationship_rows: list[dict] | None = None,
) -> pd.DataFrame:
    """Score all localizations using coordinates and barcode identities only."""
    output = []
    feature_columns = ["Trace_ID", "Spot_ID", "Barcode", "x", "y", "z"]
    if "simulation_id" in frame.columns:
        feature_columns.insert(0, "simulation_id")
    # Copy only feature columns before scoring: evaluation labels cannot affect scores.
    features = frame[feature_columns]
    diagnostic_key = (
        ["simulation_id", "Spot_ID"] if "simulation_id" in frame.columns else "Spot_ID"
    )
    diagnostics = frame.drop(columns=["x", "y", "z"]).set_index(diagnostic_key)
    for _, trace in features.groupby(trace_group_key(features), sort=False):
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
            diagnostic_id = (
                (str(row.simulation_id), str(row.Spot_ID))
                if "simulation_id" in features.columns
                else str(row.Spot_ID)
            )
            metadata = diagnostics.loc[diagnostic_id].to_dict()
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
                    "both_flanks_available": bool(
                        len(left_values) and len(right_values)
                    ),
                }
                residuals = []
                residual_attempted_supports = []
                residual_fallbacks = 0
                for j in selected:
                    separation = abs(barcode - int(barcodes[j]))
                    stats, used_fallback = model.residual_stats(
                        separation, residual_fallback
                    )
                    exact_stats = model.residual.get(separation)
                    residual_attempted_supports.append(
                        stats[2]
                        if used_fallback and stats is not None
                        else (exact_stats[2] if exact_stats is not None else 0)
                    )
                    if stats:
                        distance = float(np.linalg.norm(coords[i] - coords[j]))
                        residuals.append(((distance - stats[0]) / stats[1]) ** 2)
                        residual_fallbacks += int(used_fallback)
                output.append(
                    {
                        **base,
                        "model": "trace_splitter_like_residual",
                        "score": float(np.mean(residuals)) if residuals else math.nan,
                        "n_context_used": len(residuals),
                        "minimum_reference_support": min(
                            residual_attempted_supports, default=0
                        ),
                        "n_insufficient_reference": len(selected) - len(residuals),
                        "n_reference_fallbacks": residual_fallbacks,
                        "reference_fallback_mode": residual_fallback,
                        "score_status": "ok" if residuals else "insufficient_reference",
                    }
                )
                for strategy in ("separation", "barcode_pair"):
                    anomalies, attempted_supports = [], []
                    for j in selected:
                        other_barcode = int(barcodes[j])
                        reference_key = (
                            abs(barcode - other_barcode)
                            if strategy == "separation"
                            else tuple(sorted((barcode, other_barcode)))
                        )
                        raw_values = getattr(model, strategy).get(reference_key)
                        attempted_supports.append(
                            len(raw_values) if raw_values is not None else 0
                        )
                        values = model.values(strategy, barcode, other_barcode)
                        distance = float(np.linalg.norm(coords[i] - coords[j]))
                        if values is None:
                            if relationship_rows is not None:
                                relationship_rows.append(
                                    {
                                        **metadata,
                                        "context": context,
                                        "reference_strategy": strategy,
                                        "other_barcode": other_barcode,
                                        "genomic_separation": abs(
                                            barcode - other_barcode
                                        ),
                                        "distance": distance,
                                        "empirical_percentile": math.nan,
                                        "tail_probability": math.nan,
                                        "pair_anomaly_score": math.nan,
                                        "anomaly_direction": "insufficient_reference",
                                        "reference_support": (
                                            len(raw_values)
                                            if raw_values is not None
                                            else 0
                                        ),
                                        "relationship_status": "insufficient_reference",
                                    }
                                )
                            continue
                        percentile, tail = empirical_distance_diagnostics(
                            values, distance
                        )
                        anomalies.append(-math.log10(tail))
                        if relationship_rows is not None:
                            relationship_rows.append(
                                {
                                    **metadata,
                                    "context": context,
                                    "reference_strategy": strategy,
                                    "other_barcode": other_barcode,
                                    "genomic_separation": abs(barcode - other_barcode),
                                    "distance": distance,
                                    "empirical_percentile": percentile,
                                    "tail_probability": tail,
                                    "pair_anomaly_score": -math.log10(tail),
                                    "anomaly_direction": (
                                        "short" if percentile < 0.5 else "long"
                                    ),
                                    "reference_support": len(values),
                                    "relationship_status": "ok",
                                }
                            )
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
                                "n_context_used": len(anomalies),
                                "minimum_reference_support": min(
                                    attempted_supports, default=0
                                ),
                                "n_insufficient_reference": len(selected)
                                - len(anomalies),
                                "n_reference_fallbacks": 0,
                                "reference_fallback_mode": "none",
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


def evaluate_condition(
    scores: pd.DataFrame,
    thresholds: Mapping[tuple[str, str, float], float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate one condition using its fold's precomputed thresholds."""
    overall, fixed = [], []
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
        metadata = dict(zip(condition_keys, key))
        overall.append(
            {
                **metadata,
                "roc_auc": auc,
                "precision_recall_auc": ap,
                "n_scored": int(valid.sum()),
            }
        )
        for fpr in FIXED_FPRS:
            threshold = thresholds.get((str(key[-2]), str(key[-1]), fpr), math.nan)
            fixed.append(
                {
                    **metadata,
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
    representative: pd.DataFrame, fixed: pd.DataFrame, output: Path
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
                "trace_splitter_like_residual",
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
    zero = representative[representative["displacement"] == 0]
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
    parser.add_argument(
        "--write-localization-scores",
        action="store_true",
        help="Write the very large per-localization score table (default: disabled).",
    )
    parser.add_argument(
        "--write-pairwise-diagnostics",
        action="store_true",
        help="Write the very large pairwise diagnostic table (default: disabled).",
    )
    parser.add_argument(
        "--residual-separation-fallback",
        choices=("none", "nearest"),
        default="none",
        help=(
            "Fallback for missing residual separation bins (default: none; "
            "nearest is explicit and reported in diagnostics)."
        ),
    )
    return parser.parse_args(argv)


def run(
    root: Path,
    output: Path,
    minimum_observations: int = 20,
    residual_fallback: str = "none",
    write_localization_scores: bool = False,
    write_pairwise_diagnostics: bool = False,
) -> None:
    _progress(f"Starting analysis in {root}")
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
    _progress(
        f"Loaded manifest with {len(simulations)} simulations and "
        f"{len(conditions)} conditions"
    )
    output.mkdir(parents=True, exist_ok=True)
    baselines = {}
    for index, item in enumerate(simulations, start=1):
        simulation_id = str(item["id"])
        _progress(f"Loading baseline {index}/{len(simulations)}: {simulation_id}")
        baselines[simulation_id] = load_baseline(
            root / "simulations" / simulation_id, item
        )
    conditions_by_efficiency: dict[float, list[Mapping]] = defaultdict(list)
    for item in conditions:
        conditions_by_efficiency[float(item["detection_efficiency"])].append(item)
    performances, fixed_tables, summaries, representative_parts = [], [], [], []
    score_writer = (
        _EcsvChunkWriter(output / "localization_scores.ecsv")
        if write_localization_scores
        else None
    )
    relationship_writer = (
        _EcsvChunkWriter(output / "pairwise_relationship_scores.ecsv")
        if write_pairwise_diagnostics
        else None
    )
    rows_processed = 0
    for efficiency, condition_items in conditions_by_efficiency.items():
        seeds = sorted({int(item["seed"]) for item in condition_items})
        _progress(
            f"Processing efficiency {efficiency:g}: {len(condition_items)} "
            f"conditions across {len(seeds)} folds"
        )
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
            _progress(
                f"Fitting reference model for efficiency {efficiency:g}, "
                f"evaluation seed {evaluation_seed}"
            )
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
            calibration_accumulator = _CalibrationAccumulator()
            try:
                for calibration_seed in seeds:
                    if calibration_seed == evaluation_seed:
                        continue
                    _progress(
                        f"Scoring calibration seed {calibration_seed} for evaluation "
                        f"seed {evaluation_seed}"
                    )
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
                    calibration_scores = score_observations(
                        calibration_frame,
                        nested_model,
                        residual_fallback=residual_fallback,
                    )
                    rows_processed += len(calibration_scores)
                    calibration_accumulator.add(calibration_scores)
                    _progress(
                        f"Spooled {len(calibration_scores):,} calibration scores; "
                        f"{rows_processed:,} score rows processed overall"
                    )
                    del calibration_scores, calibration_frame, nested_model
                thresholds = calibration_accumulator.thresholds()
            finally:
                calibration_accumulator.close()
            _progress(
                f"Derived {len(thresholds):,} fixed-FPR thresholds for evaluation "
                f"seed {evaluation_seed}; raw calibration scores discarded"
            )
            evaluation_items = [
                item for item in condition_items if int(item["seed"]) == evaluation_seed
            ]
            for condition_index, item in enumerate(evaluation_items, start=1):
                _progress(
                    f"Scoring condition {condition_index}/{len(evaluation_items)} "
                    f"for seed {evaluation_seed}: {item['id']}"
                )
                observations = load_observations(root / str(item["directory"]), item)
                relationships: list[dict] | None = (
                    [] if relationship_writer is not None else None
                )
                condition_scores = score_observations(
                    observations,
                    model,
                    residual_fallback=residual_fallback,
                    relationship_rows=relationships,
                )
                rows_processed += len(condition_scores)
                overall, fixed = evaluate_condition(condition_scores, thresholds)
                performances.append(overall)
                fixed_tables.append(fixed)
                representative_parts.append(
                    condition_scores.loc[
                        (condition_scores["model"] == "empirical_separation_mean")
                        & (condition_scores["context"] == "k3"),
                        [
                            "displacement",
                            "ground_truth_corrupted",
                            "selected_for_corruption",
                            "score",
                        ],
                    ].copy()
                )
                if score_writer is not None:
                    score_writer.write(condition_scores)
                relationship_count = 0
                if relationship_writer is not None and relationships is not None:
                    relationship_table = pd.DataFrame(relationships).rename(
                        columns={
                            "is_corrupted": "ground_truth_corrupted",
                            "injected_displacement_um": "ground_truth_displacement",
                        }
                    )
                    relationship_count = len(relationship_table)
                    relationship_writer.write(relationship_table)
                _progress(
                    f"Finished condition {item['id']}: "
                    f"{len(condition_scores):,} scores and "
                    f"{relationship_count:,} pairwise diagnostics; "
                    f"{rows_processed:,} score rows processed overall"
                )
                del condition_scores, observations
    if score_writer is not None:
        score_writer.finish()
    if relationship_writer is not None:
        relationship_writer.finish()
    performance, fixed = pd.concat(performances, ignore_index=True), pd.concat(
        fixed_tables, ignore_index=True
    )
    representative = pd.concat(representative_parts, ignore_index=True)
    _progress("Generating plots")
    plot_summary = plot_outputs(representative, fixed, output)
    _progress("Writing result tables")
    _write(performance, output / "model_performance.ecsv")
    _write(fixed, output / "performance_at_fixed_fpr.ecsv")
    _write(
        pd.concat(summaries, ignore_index=True), output / "reference_model_summary.ecsv"
    )
    _write(plot_summary, output / "plot_aggregate_values.ecsv")
    _progress(f"Analysis complete; results written to {output}")


def main() -> None:
    args = parse_arguments()
    run(
        args.benchmark_root.resolve(),
        args.output_dir.resolve(),
        args.minimum_reference_observations,
        args.residual_separation_fallback,
        args.write_localization_scores,
        args.write_pairwise_diagnostics,
    )


if __name__ == "__main__":
    main()
