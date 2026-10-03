#!/usr/bin/env python3
"""Explore culprit attribution on the frozen ``benchmark_curator_v1`` data.

This is deliberately not production ``trace_curator``.  It uses only empirical
references indexed by genomic separation and preserves cycle 1's nested,
replicate-held-out reference fitting and threshold calibration.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import sys
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

FIXED_FPRS = (0.001, 0.01, 0.05)
EDGE_EXCESS_THRESHOLD = -math.log10(0.05)
_START = time.monotonic()
MODELS = {
    "empirical_separation_mean_k1": ("baseline", "k1", "mean"),
    "empirical_separation_maximum_k1": ("baseline", "k1", "maximum"),
    "loo_sum_all": ("loo", "all", "sum"),
    "loo_top3_all": ("loo", "all", "top3"),
    "star_top2_mean_all": ("star", "all", "top2"),
    "star_top3_mean_all": ("star", "all", "top3"),
    "star_excess_sum_all": ("star", "all", "excess"),
    "star_anomalous_fraction_all": ("star", "all", "fraction"),
    "empirical_two_flank_bridge": ("bridge", "k1", "empirical"),
}
COMPACT_OUTPUTS = (
    "attribution_model_performance.ecsv",
    "attribution_performance_at_fixed_fpr.ecsv",
    "culprit_ranking_performance.ecsv",
    "collateral_fpr_by_barcode_distance.ecsv",
    "culprit_rank_distribution.ecsv",
    "clean_trace_false_culprit_behavior.ecsv",
    "attribution_reference_summary.ecsv",
)


def fold_seed_sets(
    seeds: Sequence[int], evaluation_seed: int, calibration_seed: int | None = None
) -> tuple[set[int], set[int], set[int]]:
    """Return disjoint reference, calibration, and evaluation provenance sets."""
    evaluation = {int(evaluation_seed)}
    calibration = set() if calibration_seed is None else {int(calibration_seed)}
    reference = {int(seed) for seed in seeds}.difference(evaluation | calibration)
    if (reference & evaluation) or (reference & calibration) or (evaluation & calibration):
        raise AssertionError("Reference/calibration/evaluation replicate leakage")
    return reference, calibration, evaluation


def _cycle1():
    path = Path(__file__).with_name("analyze_trace_curator_benchmark.py")
    spec = importlib.util.spec_from_file_location("_frozen_trace_curator_cycle1", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _progress(message: str) -> None:
    print(f"[attribution {time.monotonic()-_START:8.1f}s] {message}", file=sys.stderr, flush=True)


def _write(frame: pd.DataFrame, path: Path) -> None:
    safe = frame.copy()
    for column in safe.select_dtypes(include="object"):
        safe[column] = safe[column].fillna("").astype(str)
    Table.from_pandas(safe).write(path, format="ascii.ecsv", overwrite=True)


def _read_compact_outputs(output: Path) -> dict[str, pd.DataFrame]:
    """Read and validate the compact results needed to reproduce the plots."""
    missing = [name for name in COMPACT_OUTPUTS if not (output / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing compact benchmark output(s) in {output}: {', '.join(missing)}"
        )
    return {
        name: Table.read(output / name, format="ascii.ecsv").to_pandas()
        for name in COMPACT_OUTPUTS
    }


def _trace_keys(frame: pd.DataFrame) -> list[str]:
    return ["simulation_id", "Trace_ID"] if "simulation_id" in frame else ["Trace_ID"]


@dataclass
class EmpiricalReference:
    separation: dict[float, np.ndarray]
    bridge: dict[tuple[float, float], tuple[np.ndarray, np.ndarray]]
    minimum_observations: int = 20


def fit_reference(frame: pd.DataFrame, minimum_observations: int = 20) -> EmpiricalReference:
    """Fit separation and two-flank references from clean traces only."""
    clean = frame.loc[~frame["is_corrupted"].astype(bool)]
    edge_values: dict[float, list[float]] = defaultdict(list)
    bridge_values: dict[tuple[float, float], list[tuple[float, float]]] = defaultdict(list)
    for _, trace in clean.groupby(_trace_keys(clean), sort=False):
        trace = trace.sort_values("Barcode")
        xyz = trace[["x", "y", "z"]].to_numpy(float)
        barcode = trace["Barcode"].to_numpy(float)
        for i in range(len(trace)):
            for j in range(i + 1, len(trace)):
                edge_values[float(barcode[j] - barcode[i])].append(float(np.linalg.norm(xyz[j]-xyz[i])))
        for i in range(1, len(trace)-1):
            key = (float(barcode[i]-barcode[i-1]), float(barcode[i+1]-barcode[i]))
            geometry = bridge_coordinates(xyz[i - 1], xyz[i], xyz[i + 1])
            if geometry is not None:
                bridge_values[key].append(geometry)
    bridges = {key: tuple(np.asarray(values, float)[:, column] for column in range(2)) for key, values in bridge_values.items()}
    return EmpiricalReference({key: np.asarray(values, float) for key, values in edge_values.items()}, bridges, minimum_observations)


def empirical_anomaly(reference: Sequence[float], observed: float) -> float:
    """Finite-sample corrected, two-sided empirical anomaly score."""
    sample = np.sort(np.asarray(reference, float))
    lower = (np.searchsorted(sample, observed, side="right") + 1) / (len(sample) + 1)
    upper = (len(sample) - np.searchsorted(sample, observed, side="left") + 1) / (len(sample) + 1)
    return -math.log10(min(1.0, 2.0 * min(lower, upper)))


def bridge_coordinates(
    left: np.ndarray, candidate: np.ndarray, right: np.ndarray
) -> tuple[float, float] | None:
    """Represent a candidate by projection and offset relative to its flanks."""
    axis = np.asarray(right, float) - np.asarray(left, float)
    squared_span = float(np.dot(axis, axis))
    if squared_span <= np.finfo(float).eps:
        return None
    candidate_from_left = np.asarray(candidate, float) - np.asarray(left, float)
    fractional_projection = float(np.dot(candidate_from_left, axis) / squared_span)
    projected = np.asarray(left, float) + fractional_projection * axis
    perpendicular_distance = float(np.linalg.norm(np.asarray(candidate, float) - projected))
    return fractional_projection, perpendicular_distance


def _context(barcodes: np.ndarray, index: int, name: str) -> list[int]:
    candidates = np.flatnonzero(np.arange(len(barcodes)) != index)
    if name == "all":
        return candidates.tolist()
    left = candidates[barcodes[candidates] < barcodes[index]]
    right = candidates[barcodes[candidates] > barcodes[index]]
    chosen = []
    if len(left): chosen.append(int(left[np.argmax(barcodes[left])]))
    if len(right): chosen.append(int(right[np.argmin(barcodes[right])]))
    return chosen


def trace_cost(scores: Sequence[float], variant: str) -> float:
    values = np.asarray(scores, float)
    values = values[np.isfinite(values)]
    if not len(values): return math.nan
    if variant == "sum": return float(values.sum())
    if variant == "top3": return float(np.sort(values)[-min(3, len(values)):].sum())
    raise ValueError(variant)


def leave_one_out_improvement(n_nodes: int, edges: Sequence[tuple[int, int, float]], variant: str) -> np.ndarray:
    """Return cost improvement; the sum variant equals incident-edge sum.

    Consequently ``loo_sum_all`` is retained as the most transparent LOO
    baseline, but is not an attribution principle independent of a star sum.
    ``loo_top3_all`` remains distinct because removing a node can change which
    edges enter the whole-trace top-three set.
    """
    before = trace_cost([score for _, _, score in edges], variant)
    result = np.full(n_nodes, np.nan)
    for node in range(n_nodes):
        retained = [score for left, right, score in edges if node not in (left, right)]
        after = trace_cost(retained, variant)
        if np.isfinite(after): result[node] = before - after
    return result


def star_scores(n_nodes: int, edges: Sequence[tuple[int, int,float]], variant: str) -> np.ndarray:
    result = np.full(n_nodes, np.nan)
    for node in range(n_nodes):
        values = np.asarray([score for left, right, score in edges if node in (left, right)], float)
        if not len(values): continue
        if variant in {"top2", "top3"}:
            count = 2 if variant == "top2" else 3
            # Repeated evidence: nodes with fewer than two supported edges are unscorable.
            if len(values) < 2: continue
            result[node] = np.sort(values)[-min(count, len(values)):].mean()
        elif variant == "excess":
            result[node] = np.maximum(values-EDGE_EXCESS_THRESHOLD, 0).sum()
        elif variant == "fraction":
            result[node] = np.mean(values >= EDGE_EXCESS_THRESHOLD)
        else: raise ValueError(variant)
    return result


def bridge_score(trace: pd.DataFrame, index: int, reference: EmpiricalReference) -> tuple[float, str]:
    """Score candidate projection and perpendicular offset relative to flanks."""
    if index == 0 or index == len(trace)-1:
        return math.nan, "missing_flank"
    barcode = trace["Barcode"].to_numpy(float)
    xyz = trace[["x", "y", "z"]].to_numpy(float)
    key = (float(barcode[index]-barcode[index-1]), float(barcode[index+1]-barcode[index]))
    distributions = reference.bridge.get(key)
    if distributions is None or any(len(values) < reference.minimum_observations for values in distributions):
        return math.nan, "insufficient_reference"
    observed = bridge_coordinates(xyz[index - 1], xyz[index], xyz[index + 1])
    if observed is None:
        return math.nan, "degenerate_flank_geometry"
    # Both quantities depend on the candidate.  The score is transparent and
    # empirical without asserting a multivariate parametric distribution.
    return float(np.mean([empirical_anomaly(values, value) for values, value in zip(distributions, observed)])), "ok"


def score_trace(trace: pd.DataFrame, reference: EmpiricalReference) -> pd.DataFrame:
    trace = trace.sort_values("Barcode").reset_index(drop=True)
    xyz = trace[["x", "y", "z"]].to_numpy(float)
    barcode = trace["Barcode"].to_numpy(float)
    edges, lookup = [], {}
    for i in range(len(trace)):
        for j in range(i+1, len(trace)):
            values = reference.separation.get(float(barcode[j]-barcode[i]))
            if values is None or len(values) < reference.minimum_observations: continue
            score = empirical_anomaly(values, float(np.linalg.norm(xyz[j]-xyz[i])))
            edges.append((i, j, score)); lookup[(i,j)] = score
    scores: dict[str, np.ndarray] = {}
    for aggregator in ("mean", "maximum"):
        values = []
        for i in range(len(trace)):
            incident = [lookup[tuple(sorted((i,j)))] for j in _context(barcode, i, "k1") if tuple(sorted((i,j))) in lookup]
            values.append((np.mean(incident) if aggregator == "mean" else np.max(incident)) if incident else math.nan)
        scores[f"empirical_separation_{aggregator}_k1"] = np.asarray(values)
    scores["loo_sum_all"] = leave_one_out_improvement(len(trace), edges, "sum")
    scores["loo_top3_all"] = leave_one_out_improvement(len(trace), edges, "top3")
    for variant in ("top2", "top3", "excess", "fraction"):
        name = {"top2":"star_top2_mean_all", "top3":"star_top3_mean_all", "excess":"star_excess_sum_all", "fraction":"star_anomalous_fraction_all"}[variant]
        scores[name] = star_scores(len(trace), edges, variant)
    bridge = [bridge_score(trace, i, reference) for i in range(len(trace))]
    scores["empirical_two_flank_bridge"] = np.asarray([item[0] for item in bridge])
    identities = [column for column in ("simulation_id","Trace_ID","Spot_ID","Barcode") if column in trace]
    parts = []
    for model, values in scores.items():
        part = trace[identities].copy(); part["model"] = model; part["score"] = values
        part["score_status"] = [item[1] for item in bridge] if model == "empirical_two_flank_bridge" else np.where(np.isfinite(values), "ok", "insufficient_reference")
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def score_observations(frame: pd.DataFrame, reference: EmpiricalReference) -> pd.DataFrame:
    parts = [score_trace(trace, reference) for _, trace in frame.groupby(_trace_keys(frame), sort=False)]
    scored = pd.concat(parts, ignore_index=True)
    keys = [column for column in ("simulation_id", "Spot_ID") if column in frame]
    metadata = [column for column in ("condition","replicate","seed","detection_efficiency","displacement","is_corrupted","selected_for_corruption","injected_displacement_um") if column in frame]
    return scored.merge(frame[keys+metadata], on=keys, how="left", validate="many_to_one")


def fixed_fpr_threshold(values: Sequence[float], fpr: float) -> float:
    values = np.sort(np.asarray(values, float)); values = values[np.isfinite(values)]
    if not len(values): return math.nan
    allowed = int(math.floor(fpr*len(values)))
    return float(values[max(0, len(values)-allowed-1)]) if allowed else float(values[-1])


def collateral_distance(barcode, target_barcode) -> str:
    distance = abs(float(barcode)-float(target_barcode))
    if distance >= 4: return ">=4"
    return str(int(distance))


def rank_trace(group: pd.DataFrame) -> dict:
    valid = group[np.isfinite(group.score)]
    target = valid[valid.is_corrupted.astype(bool)]
    if len(target) != 1:
        return {"scoreable": False, "rank": math.nan, "top1": False,
                "top1_unique": False, "top2": False,
                "top2_including_ties": False, "top2_conservative": False,
                "n_tied_at_target_rank": 0, "reciprocal_rank": 0.0}
    target_score = float(target.score.iloc[0])
    rank = int(1 + np.sum(valid.score.to_numpy() > target_score))
    tied = int(np.sum(valid.score.to_numpy() == target_score))
    top1_unique = rank == 1 and tied == 1
    top2_including_ties = rank <= 2
    # Conservatively require the complete set ranked at least as highly as the
    # target to fit into two slots.  A three-way maximum therefore fails.
    top2_conservative = int(np.sum(valid.score.to_numpy() >= target_score)) <= 2
    return {"scoreable": True, "rank": rank, "top1": top1_unique,
            "top1_unique": top1_unique, "top2": top2_conservative,
            "top2_including_ties": top2_including_ties,
            "top2_conservative": top2_conservative,
            "n_tied_at_target_rank": tied, "reciprocal_rank": 1/rank}


def evaluate(scores: pd.DataFrame, thresholds: Mapping[tuple[str,float],float]):
    keys = ["condition","replicate","seed","detection_efficiency","displacement","model"]
    overall, fixed, ranks, collateral = [], [], [], []
    for key, group in scores.groupby(keys, sort=False):
        meta = dict(zip(keys,key)); valid = group[np.isfinite(group.score)]
        labels = valid.is_corrupted.astype(bool).to_numpy(); values = valid.score.to_numpy()
        trace_keys = _trace_keys(group)
        trace_state = (
            group.groupby(trace_keys, sort=False, dropna=False)["is_corrupted"]
            .any().rename("trace_has_target").reset_index()
        )
        valid_with_state = valid.merge(
            trace_state, on=trace_keys, how="left", validate="many_to_one"
        )
        trace_has_target = valid_with_state.trace_has_target.astype(bool).to_numpy()
        overall.append({**meta,"roc_auc":roc_auc_score(labels,values) if len(np.unique(labels))==2 else math.nan,"precision_recall_auc":average_precision_score(labels,values) if labels.any() else math.nan,"n_scored":len(valid)})
        for fpr in FIXED_FPRS:
            threshold=thresholds.get((meta["model"],fpr),math.nan); calls=valid.score.to_numpy()>threshold
            # A corrupted trace has exactly one positive-displacement target; r=0
            # selected targets remain clean and contribute to clean-trace FPR.
            target = labels; clean = ~trace_has_target; collateral_mask = trace_has_target & ~target
            n_targets = int(group.is_corrupted.astype(bool).sum())
            fixed.append({**meta,"target_fpr":fpr,"threshold":threshold,"target_sensitivity":np.sum(calls & target)/n_targets if n_targets else math.nan,"clean_trace_fpr":np.mean(calls[clean]) if clean.any() else math.nan,"collateral_fpr":np.mean(calls[collateral_mask]) if collateral_mask.any() else math.nan,"selected_zero_call_rate":np.mean(calls[valid.selected_for_corruption.astype(bool).to_numpy()]) if np.isclose(meta["displacement"],0) and valid.selected_for_corruption.astype(bool).any() else math.nan})
        primary=thresholds.get((meta["model"],.01),math.nan)
        for trace_key, trace in group.groupby(trace_keys,sort=False):
            trace_meta={name:value for name,value in zip(trace_keys,trace_key if isinstance(trace_key,tuple) else (trace_key,))}
            corrupted=trace[trace.is_corrupted.astype(bool)]
            if len(corrupted)==1:
                ranks.append({**meta,**trace_meta,**rank_trace(trace)})
                target_barcode=corrupted.Barcode.iloc[0]
                for row in trace.loc[~trace.is_corrupted.astype(bool)].itertuples():
                    if np.isfinite(row.score): collateral.append({**meta,**trace_meta,"barcode_distance":collateral_distance(row.Barcode,target_barcode),"called":bool(row.score>primary)})
    return pd.DataFrame(overall),pd.DataFrame(fixed),pd.DataFrame(ranks),pd.DataFrame(collateral)

def _aggregate(frame: pd.DataFrame, metrics: Sequence[str]) -> pd.DataFrame:
    keys=["model","detection_efficiency","displacement"]
    rows=[]
    for key,group in frame.groupby(keys,sort=False):
        meta=dict(zip(keys,key))
        replicate=group.groupby("seed")[list(metrics)].mean(numeric_only=True)
        for metric in metrics:
            values=replicate[metric].dropna()
            rows.append({**meta,"metric":metric,"mean":values.mean(),"sd":values.std(),"sem":values.sem(),"n_replicates":len(values)})
    return pd.DataFrame(rows)


def summarize_ranking(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize culprit ranks without hiding unscoreable-target failures."""
    keys = ["model", "detection_efficiency", "displacement"]
    rows = []
    for key, group in frame.groupby(keys, sort=False):
        scoreable = group[group.scoreable.astype(bool)]
        rows.append(
            {
                **dict(zip(keys, key)),
                "n_corrupted_traces": len(group),
                "attribution_coverage": group.scoreable.mean(),
                # top-k values are False for an unscoreable target, making these
                # unconditional attribution accuracies.
                "top1_attribution_accuracy": group.top1_unique.mean(),
                "top1_unique_accuracy": group.top1_unique.mean(),
                "top2_attribution_accuracy": group.top2_conservative.mean(),
                "top2_conservative_accuracy": group.top2_conservative.mean(),
                "top2_including_ties_accuracy": group.top2_including_ties.mean(),
                "mean_n_tied_at_target_rank": group.n_tied_at_target_rank.mean(),
                "mean_true_target_rank": scoreable["rank"].mean(),
                "median_true_target_rank": scoreable["rank"].median(),
                "mean_reciprocal_rank": group.reciprocal_rank.mean(),
            }
        )
    return pd.DataFrame(rows)


def _plot_lines(frame: pd.DataFrame, metric: str, output: Path, filename: str, ylabel: str) -> None:
    if frame.empty or metric not in frame: return
    aggregate=_aggregate(frame,[metric])
    fig,axis=plt.subplots(figsize=(10,6))
    for (model,efficiency),line in aggregate.groupby(["model","detection_efficiency"]):
        line=line.sort_values("displacement")
        axis.errorbar(line["displacement"],line["mean"],yerr=line["sem"],marker="o",label=f"{model}; eff={efficiency:g}")
    axis.set(xlabel="Injected displacement (µm)",ylabel=ylabel); axis.legend(fontsize=6,ncol=2)
    fig.tight_layout(); fig.savefig(output/filename,dpi=160); plt.close(fig)


def plot_results(fixed: pd.DataFrame, ranking: pd.DataFrame, collateral: pd.DataFrame, clean_maxima: pd.DataFrame, output: Path) -> None:
    primary=fixed[np.isclose(fixed["target_fpr"],.01)]
    for metric,name,label in (("target_sensitivity","target_sensitivity.png","Target sensitivity at 1% FPR"),("clean_trace_fpr","clean_trace_fpr.png","Clean-trace FPR"),("collateral_fpr","collateral_fpr.png","Collateral FPR")):
        _plot_lines(primary,metric,output,name,label)
    _plot_lines(ranking,"top1",output,"top1_culprit_attribution.png","Unique top-1 attribution accuracy")
    if not collateral.empty:
        keys = ["model", "detection_efficiency", "displacement", "barcode_distance"]
        totals = collateral.groupby(keys)[["n_calls", "n_genuine"]].sum().reset_index()
        totals["collateral_fpr"] = totals["n_calls"] / totals["n_genuine"]
        fig,axis=plt.subplots(figsize=(10,6))
        for (model,distance),line in totals.groupby(["model","barcode_distance"]):
            axis.plot(line["displacement"],line["collateral_fpr"],marker="o",label=f"{model}; Δbarcode={distance}")
        axis.set(xlabel="Injected displacement (µm)",ylabel="Collateral FPR at 1%",ylim=(-.01,1.01)); axis.legend(fontsize=5,ncol=3)
        fig.tight_layout(); fig.savefig(output/"collateral_fpr_by_barcode_distance.png",dpi=160); plt.close(fig)
    if not ranking.empty:
        fig,axis=plt.subplots(figsize=(9,5))
        for model,group in ranking.groupby("model"):
            ranks=group["rank"].dropna(); counts=ranks.value_counts(normalize=True).sort_index()
            axis.step(counts.index,counts.values,where="mid",label=model)
        axis.set(xlabel="Rank of true target",ylabel="Fraction of corrupted traces"); axis.legend(fontsize=6,ncol=2)
        fig.tight_layout(); fig.savefig(output/"true_target_rank_distribution.png",dpi=160); plt.close(fig)
    if not clean_maxima.empty:
        fig,axis=plt.subplots(figsize=(9,5))
        for model,group in clean_maxima.groupby("model"):
            axis.hist(group["maximum_score"].dropna(),bins=40,density=True,histtype="step",label=model)
        axis.set(xlabel="Maximum localization score in a clean trace",ylabel="Density"); axis.legend(fontsize=6,ncol=2)
        fig.tight_layout(); fig.savefig(output/"clean_trace_false_culprit_maxima.png",dpi=160); plt.close(fig)


def run(root: Path, output: Path, minimum_observations: int=20, write_localization_scores: bool=False) -> None:
    """Run bounded-memory nested cross-validation, one condition at a time."""
    cycle1=_cycle1(); manifest_path=root/"sweep_manifest.yaml"
    if not manifest_path.is_file(): raise FileNotFoundError(manifest_path)
    manifest=yaml.safe_load(manifest_path.read_text())
    if manifest.get("dry_run"): raise ValueError("Cannot analyze a dry-run benchmark")
    simulations,conditions=manifest.get("simulations"),manifest.get("conditions")
    if not simulations or not conditions: raise ValueError("Manifest requires simulations and conditions")
    output.mkdir(parents=True,exist_ok=True)
    baselines={str(item["id"]):cycle1.load_baseline(root/"simulations"/str(item["id"]),item) for item in simulations}
    by_efficiency: dict[float,list[Mapping]]=defaultdict(list)
    for item in conditions: by_efficiency[float(item["detection_efficiency"])].append(item)
    performance_parts=[]; fixed_parts=[]; ranking_parts=[]; collateral_parts=[]; clean_max_parts=[]; reference_rows=[]
    writer=cycle1._EcsvChunkWriter(output/"attribution_localization_scores.ecsv") if write_localization_scores else None
    for efficiency,items in by_efficiency.items():
        seeds=sorted({int(item["seed"]) for item in items})
        if len(seeds)<3: raise ValueError("At least three replicates are required for nested calibration")
        baseline=pd.concat([frame for frame in baselines.values() if np.isclose(frame.detection_efficiency.iloc[0],efficiency)],ignore_index=True)
        for evaluation_seed in seeds:
            _progress(f"efficiency={efficiency:g}; evaluation seed={evaluation_seed}")
            reference_seeds, _, _ = fold_seed_sets(seeds, evaluation_seed)
            reference_data=baseline[baseline.seed.isin(reference_seeds)]
            reference=fit_reference(reference_data,minimum_observations)
            for key,values in reference.separation.items(): reference_rows.append({"detection_efficiency":efficiency,"evaluation_seed":evaluation_seed,"reference_type":"separation","reference_key":str(key),"n_observations":len(values),"sufficient":len(values)>=minimum_observations})
            for key,values in reference.bridge.items(): reference_rows.append({"detection_efficiency":efficiency,"evaluation_seed":evaluation_seed,"reference_type":"bridge","reference_key":str(key),"n_observations":len(values[0]),"sufficient":len(values[0])>=minimum_observations})
            calibration={name:[] for name in MODELS}
            # Exactly cycle 1's nested protection: both evaluation and calibration
            # replicates are excluded while scoring a clean calibration replicate.
            for calibration_seed in seeds:
                if calibration_seed==evaluation_seed: continue
                nested_seeds, _, _ = fold_seed_sets(
                    seeds, evaluation_seed, calibration_seed
                )
                nested_data=baseline[baseline.seed.isin(nested_seeds)]
                nested_reference=fit_reference(nested_data,minimum_observations)
                clean=baseline[baseline.seed==calibration_seed].copy()
                clean["condition"]="calibration"; clean["replicate"]=calibration_seed; clean["displacement"]=0.; clean["selected_for_corruption"]=False; clean["injected_displacement_um"]=0.
                nested_scores=score_observations(clean,nested_reference)
                for name,group in nested_scores.groupby("model"): calibration[name].extend(group.score.dropna().tolist())
                del nested_scores,clean,nested_reference
            thresholds={(name,fpr):fixed_fpr_threshold(values,fpr) for name,values in calibration.items() for fpr in FIXED_FPRS}
            del calibration
            for item in items:
                if int(item["seed"])!=evaluation_seed: continue
                _progress(f"scoring {item['id']}")
                observations=cycle1.load_observations(root/str(item["directory"]),item)
                scores=score_observations(observations,reference)
                overall,fixed,ranks,collateral=evaluate(scores,thresholds)
                performance_parts.append(overall); fixed_parts.append(fixed)
                if not ranks.empty: ranking_parts.append(ranks)
                if not collateral.empty: collateral_parts.append(collateral)
                clean=scores.groupby([*_trace_keys(scores),"model"],sort=False).filter(lambda group:not group.is_corrupted.astype(bool).any())
                maxima=clean.groupby([*_trace_keys(clean),"model"],sort=False).score.max().reset_index(name="maximum_score")
                maxima["condition"]=str(item["id"]); maxima["seed"]=evaluation_seed; maxima["detection_efficiency"]=efficiency; maxima["displacement"]=float(item["displacement_um"])
                clean_max_parts.append(maxima)
                if writer is not None: writer.write(scores)
                del observations,scores,overall,fixed,ranks,collateral,clean,maxima
    if writer is not None: writer.finish()
    performance=pd.concat(performance_parts,ignore_index=True); fixed=pd.concat(fixed_parts,ignore_index=True)
    ranking=pd.concat(ranking_parts,ignore_index=True) if ranking_parts else pd.DataFrame()
    collateral=pd.concat(collateral_parts,ignore_index=True) if collateral_parts else pd.DataFrame()
    clean_maxima=pd.concat(clean_max_parts,ignore_index=True)
    ranking_summary=summarize_ranking(ranking) if not ranking.empty else pd.DataFrame()
    collateral_summary=collateral.groupby(["condition","replicate","seed","detection_efficiency","displacement","model","barcode_distance"],sort=False).called.agg(["mean","sum","count"]).reset_index().rename(columns={"mean":"collateral_fpr","sum":"n_calls","count":"n_genuine"}) if not collateral.empty else pd.DataFrame()
    _write(performance,output/"attribution_model_performance.ecsv")
    _write(fixed,output/"attribution_performance_at_fixed_fpr.ecsv")
    _write(ranking_summary,output/"culprit_ranking_performance.ecsv")
    _write(collateral_summary,output/"collateral_fpr_by_barcode_distance.ecsv")
    _write(ranking,output/"culprit_rank_distribution.ecsv")
    _write(clean_maxima,output/"clean_trace_false_culprit_behavior.ecsv")
    _write(pd.DataFrame(reference_rows),output/"attribution_reference_summary.ecsv")
    plot_results(fixed,ranking,collateral_summary,clean_maxima,output)
    _progress(f"complete: {output}")


def plot_only(output: Path) -> None:
    """Regenerate plots from a completed run's compact ECSV outputs."""
    frames = _read_compact_outputs(output)
    plot_results(
        frames["attribution_performance_at_fixed_fpr.ecsv"],
        frames["culprit_rank_distribution.ecsv"],
        frames["collateral_fpr_by_barcode_distance.ecsv"],
        frames["clean_trace_false_culprit_behavior.ecsv"],
        output,
    )
    _progress(f"plots complete: {output}")


def parse_arguments(argv: Sequence[str]|None=None) -> argparse.Namespace:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root",type=Path,default=Path("benchmark_curator_v1"))
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--minimum-reference-observations",type=int,default=20)
    parser.add_argument("--write-localization-scores",action="store_true",help="Write large localization diagnostics (disabled by default).")
    parser.add_argument("--plot-only",action="store_true",help="Regenerate plots from compact ECSV outputs without rerunning the benchmark.")
    return parser.parse_args(argv)


def main() -> None:
    args=parse_arguments()
    output = args.output_dir.resolve()
    if args.plot_only:
        plot_only(output)
    else:
        run(args.benchmark_root.resolve(),output,args.minimum_reference_observations,args.write_localization_scores)


if __name__ == "__main__": main()
