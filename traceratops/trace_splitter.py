#!/usr/bin/env python3
"""Resolve ambiguous chromatin traces with spatial or barcode-aware methods.

By default, traces whose radius of gyration is greater than the population mean
plus one standard deviation are split with K-means.  ``--split-all`` disables
that size filter.  HDBSCAN can be selected as an alternative and its distance
epsilon can either be supplied or estimated from the input traces. Beam search
uses barcode multiplicity, genomic ordering, and empirical spatial continuity
while allowing detections to remain unassigned.
"""

import argparse
import inspect
import os
import select
import sys
import uuid

import numpy as np
from astropy.table import Column, Table, vstack
from scipy.spatial.distance import pdist
from sklearn.cluster import HDBSCAN, KMeans

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.trace_resolver import (
    ClassificationThresholds,
    EmpiricalDistanceModel,
    TraceClassification,
    TraceMultiplicity,
    TraceResolver,
    classify_trace,
)
from traceratops.script_banner import print_script_banner


def parse_arguments():
    """Return the command-line argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Path to the input trace file.")
    parser.add_argument(
        "--output", help="Path to save the modified trace file (default: *_split.ecsv)."
    )
    parser.add_argument(
        "--std_threshold",
        type=float,
        default=1.0,
        help="Standard deviations above mean Rg used to select traces (default: 1.0).",
    )
    parser.add_argument(
        "--split-all",
        action="store_true",
        help="Cluster every trace instead of selecting traces by radius of gyration.",
    )
    parser.add_argument(
        "--method",
        "--clustering-method",
        dest="clustering_method",
        choices=("kmeans", "hdbscan", "beam"),
        default="kmeans",
        help="Resolution method (default: kmeans). --clustering-method is an alias.",
    )
    parser.add_argument(
        "--num_clusters",
        type=int,
        default=2,
        help="Number of K-means clusters (default: 2).",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=3,
        help="Minimum HDBSCAN cluster size (default: 3).",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=3,
        help="HDBSCAN core-point neighborhood size (default: 3).",
    )
    parser.add_argument(
        "--cluster-selection-method",
        choices=("eom", "leaf"),
        default="eom",
        help="HDBSCAN cluster selection method (default: eom).",
    )
    parser.add_argument(
        "--cluster-selection-epsilon",
        type=float,
        help="HDBSCAN distance threshold (default: median per-trace Otsu estimate).",
    )
    parser.add_argument(
        "--allow-single-cluster",
        action="store_true",
        help="Allow HDBSCAN to return a single cluster.",
    )
    beam = parser.add_argument_group("Barcode-aware beam resolution")
    beam.add_argument(
        "--history-mode",
        choices=("nearest", "multi"),
        default="multi",
        help="Polymer history used for transition scoring (default: multi).",
    )
    beam.add_argument(
        "--distance-score",
        choices=("residual", "likelihood"),
        default="residual",
        help="Empirical distance cost (default: residual).",
    )
    beam.add_argument("--history-length", type=int, default=3)
    beam.add_argument("--beam-width", type=int, default=100)
    beam.add_argument(
        "--rejection-cost",
        type=float,
        default=16.0,
        help="Cost for leaving one detection unassigned (default: 4.0).",
    )
    beam.add_argument("--minimum-polymer-size", type=int, default=2)
    beam.add_argument(
        "--minimum-confidence",
        type=float,
        default=0.05,
        help="Minimum non-probabilistic score-gap confidence for a split.",
    )
    beam.add_argument("--split-min-repeated-barcodes", type=int, default=3)
    beam.add_argument("--split-min-repeated-fraction", type=float, default=0.3)
    beam.add_argument("--candidate-rule", choices=("both", "either"), default="both")
    beam.add_argument("--model-min-observations", type=int, default=3)
    beam.add_argument("--variance-floor", type=float, default=1e-6)
    beam.add_argument(
        "--diagnostics-output",
        help="Beam diagnostics ECSV path (default: *_split_diagnostics.ecsv).",
    )
    parser.add_argument(
        "--pipe", help="Input trace-file list from stdin.", action="store_true"
    )
    return parser


def create_dict_args(args):
    p = {"trace_files": [], "pipe": args.pipe}
    if args.pipe:
        if select.select([sys.stdin], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["trace_files"] = [args.input]
    return p


def generate_unique_id():
    """Generate a unique trace identifier."""
    return str(uuid.uuid4())


def compute_radius_of_gyration(coords):
    """Compute the radius of gyration for a trace."""
    center_of_mass = np.mean(coords, axis=0)
    return np.sqrt(np.mean(np.sum((coords - center_of_mass) ** 2, axis=1)))


def _coordinates(trace):
    """Return coordinates in the dtype expected by HDBSCAN's Cython routines.

    ECSV files commonly store coordinates as float32.  Some scikit-learn
    HDBSCAN tree paths fail on float32 input when a non-zero cluster selection
    epsilon is combined with ``allow_single_cluster``.  Converting at this
    boundary also guarantees a C-contiguous matrix for both clusterers.
    """
    return np.ascontiguousarray(
        np.column_stack((trace["x"], trace["y"], trace["z"])), dtype=np.float64
    )


def _otsu_threshold(values, bins=256):
    """Compute an Otsu threshold without requiring an image-processing package."""
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.all(values == values[0]):
        return float(values[0]) if values.size else 0.0
    counts, edges = np.histogram(values, bins=min(bins, max(2, values.size)))
    centers = (edges[:-1] + edges[1:]) / 2
    weights_left = np.cumsum(counts)
    weights_right = counts.sum() - weights_left
    means_left = np.cumsum(counts * centers) / np.maximum(weights_left, 1)
    reverse_sum = np.cumsum((counts * centers)[::-1])[::-1]
    means_right = reverse_sum / np.maximum(weights_right + counts, 1)
    variance = (
        weights_left[:-1]
        * weights_right[:-1]
        * (means_left[:-1] - means_right[1:]) ** 2
    )
    return float(centers[np.argmax(variance)])


def estimate_hdbscan_epsilon(groups):
    """Estimate epsilon as the median per-trace Otsu pair-distance threshold."""
    thresholds = []
    for trace in groups:
        coords = _coordinates(trace)
        if len(coords) > 1:
            thresholds.append(_otsu_threshold(pdist(coords)))
    return float(np.median(thresholds)) if thresholds else 0.0


def _create_hdbscan(**parameters):
    """Create HDBSCAN while remaining compatible with scikit-learn 1.3+.

    Recent scikit-learn releases expose ``copy`` and warn when it is omitted;
    older supported releases do not accept it.
    """
    if "copy" in inspect.signature(HDBSCAN).parameters:
        parameters["copy"] = True
    return HDBSCAN(**parameters)


def _fit_hdbscan(coords, parameters):
    """Fit HDBSCAN, working around an upstream epsilon-tree conversion bug."""
    try:
        return _create_hdbscan(**parameters).fit_predict(coords)
    except TypeError as error:
        scalar_conversion_error = (
            "only 0-dimensional arrays can be converted to Python scalars"
        )
        if (
            scalar_conversion_error not in str(error)
            or parameters["cluster_selection_epsilon"] == 0
        ):
            raise

        # Some scikit-learn/NumPy combinations fail inside epsilon_search for
        # particular condensed trees. HDBSCAN itself remains usable without
        # epsilon-based merging, so retry only this known upstream failure.
        print(
            "! Warning: scikit-learn HDBSCAN failed while applying "
            "cluster_selection_epsilon; retrying this trace without "
            "epsilon-based cluster merging."
        )
        fallback_parameters = parameters.copy()
        fallback_parameters["cluster_selection_epsilon"] = 0.0
        return _create_hdbscan(**fallback_parameters).fit_predict(coords)


def split_large_traces(
    trace_table,
    std_threshold=1.0,
    num_clusters=2,
    split_all=False,
    clustering_method="kmeans",
    min_cluster_size=3,
    min_samples=3,
    cluster_selection_method="eom",
    cluster_selection_epsilon=None,
    allow_single_cluster=False,
):
    """Split selected traces in place with K-means or HDBSCAN.

    HDBSCAN noise points (label ``-1``) retain their original trace ID, so no
    detections are discarded.  A trace is counted as split only when the
    clustering creates at least two resulting trace IDs.
    """
    grouped = trace_table.data.group_by("Trace_ID")
    groups = list(grouped.groups)
    rg_values = [compute_radius_of_gyration(_coordinates(t)) for t in groups]
    mean_rg, std_rg = np.mean(rg_values), np.std(rg_values)
    rg_threshold = mean_rg + std_threshold * std_rg
    print(
        f"$ Mean Rg: {mean_rg:.3f}, Std Rg: {std_rg:.3f}, Threshold: {rg_threshold:.3f}"
    )

    epsilon = cluster_selection_epsilon
    if clustering_method == "hdbscan" and epsilon is None:
        epsilon = estimate_hdbscan_epsilon(groups)
        print(f"$ Estimated HDBSCAN cluster selection epsilon: {epsilon:.3f}")

    new_data = trace_table.data.copy()
    num_splits = 0
    for trace, rg in zip(groups, rg_values):
        original_id = trace["Trace_ID"][0]
        coords = _coordinates(trace)
        if not (split_all or rg > rg_threshold):
            continue
        if clustering_method == "kmeans":
            if len(coords) <= num_clusters:
                continue
            labels = KMeans(
                n_clusters=num_clusters, random_state=42, n_init=10
            ).fit_predict(coords)
        else:
            if len(coords) < min_cluster_size:
                continue
            hdbscan_parameters = {
                "min_cluster_size": min_cluster_size,
                "min_samples": min_samples,
                "metric": "euclidean",
                "cluster_selection_method": cluster_selection_method,
                "cluster_selection_epsilon": epsilon,
                "allow_single_cluster": allow_single_cluster,
            }
            labels = _fit_hdbscan(coords, hdbscan_parameters)

        cluster_labels = np.unique(labels[labels >= 0])
        number_of_parts = len(cluster_labels) + int(np.any(labels == -1))
        if number_of_parts < 2:
            continue

        original_indices = np.flatnonzero(trace_table.data["Trace_ID"] == original_id)
        for label in cluster_labels:
            new_id = generate_unique_id()
            cluster_indices = original_indices[np.flatnonzero(labels == label)]
            new_data["Trace_ID"][cluster_indices] = new_id
        num_splits += 1

    print(f"$ Number of traces split: {num_splits}/{len(groups)}")
    trace_table.data = new_data


def _diagnostic_row(
    trace_id, trace, multiplicity, rg, requested_n_polymers, result=None
):
    """Build one stable diagnostics record for an input trace."""
    if result is None:
        status = "unchanged"
        inferred = 1
        best = alternative = gap = confidence = np.nan
        n_unassigned = n_ambiguous = 0
    else:
        status = result.status
        inferred = result.inferred_n_polymers
        best = result.best_score
        alternative = (
            result.alternative_score if result.alternative_score is not None else np.nan
        )
        gap = result.raw_score_gap if result.raw_score_gap is not None else np.nan
        confidence = (
            result.confidence_score if result.confidence_score is not None else np.nan
        )
        n_unassigned = result.n_unassigned
        n_ambiguous = result.n_ambiguous_barcodes
    return {
        "input_trace_id": str(trace_id),
        "status": status,
        "n_input_localizations": len(trace),
        "n_unique_barcodes": multiplicity.n_unique_barcodes,
        "n_repeated_barcodes": multiplicity.n_repeated_barcodes,
        "fraction_repeated_barcodes": multiplicity.fraction_repeated_barcodes,
        "maximum_barcode_multiplicity": multiplicity.maximum_barcode_multiplicity,
        "n_excess_detections": multiplicity.n_excess_detections,
        "radius_of_gyration": rg,
        "requested_n_polymers": requested_n_polymers,
        "inferred_n_polymers": inferred,
        "best_score": best,
        "alternative_score": alternative,
        "raw_score_gap": gap,
        "confidence_score": confidence,
        "n_unassigned": n_unassigned,
        "fraction_unassigned": n_unassigned / len(trace) if len(trace) else 0.0,
        "n_ambiguous_barcodes": n_ambiguous,
    }


def resolve_traces(
    trace_table,
    thresholds=None,
    history_mode="multi",
    distance_score="residual",
    history_length=3,
    beam_width=100,
    rejection_cost=4.0,
    minimum_polymer_size=2,
    minimum_confidence=0.05,
    model_minimum_observations=3,
    variance_floor=1e-6,
):
    """Resolve barcode duplicates and merged pairs, returning diagnostics.

    Surviving rows retain their original ``Spot_ID`` exactly. A confident split
    receives new ``Trace_ID`` values; rejected detections are omitted.
    """
    required = {"Trace_ID", "Barcode #", "x", "y", "z"}
    missing = required.difference(trace_table.data.colnames)
    if missing:
        raise KeyError(f"Trace table is missing required columns: {sorted(missing)}")
    thresholds = thresholds or ClassificationThresholds()
    model = EmpiricalDistanceModel(
        minimum_observations=model_minimum_observations,
        variance_floor=variance_floor,
    ).fit(trace_table.data)
    if model.fallback_source == "nearest_candidate_pairs":
        print(
            "! Warning: no clean traces were available for the empirical "
            "distance model; using nearest candidate pairs as a weak fallback."
        )
    resolver = TraceResolver(
        model,
        history_mode=history_mode,
        distance_score=distance_score,
        history_length=history_length,
        beam_width=beam_width,
        rejection_cost=rejection_cost,
        minimum_polymer_size=minimum_polymer_size,
        minimum_confidence=minimum_confidence,
    )

    output_groups = []
    diagnostics = []
    for trace in trace_table.data.group_by("Trace_ID").groups:
        trace_id = trace["Trace_ID"][0]
        multiplicity = TraceMultiplicity.from_barcodes(trace["Barcode #"])
        rg = compute_radius_of_gyration(_coordinates(trace))
        classification = classify_trace(multiplicity, thresholds)
        if classification == TraceClassification.UNCHANGED:
            output_groups.append(trace.copy())
            diagnostics.append(_diagnostic_row(trace_id, trace, multiplicity, rg, 0))
            continue

        n_polymers = 2 if classification == TraceClassification.RESOLVE_TWO else 1
        result = resolver.resolve(trace, n_polymers)
        diagnostics.append(
            _diagnostic_row(trace_id, trace, multiplicity, rg, n_polymers, result)
        )
        if result.status == "removed_ambiguous":
            continue
        if n_polymers == 1:
            output_groups.append(trace[result.assignments == 0].copy())
        else:
            for polymer in range(2):
                resolved = trace[result.assignments == polymer].copy()
                if len(resolved):
                    new_id = generate_unique_id()
                    resolved.replace_column(
                        "Trace_ID", Column([new_id] * len(resolved), name="Trace_ID")
                    )
                    output_groups.append(resolved)

    if output_groups:
        output = vstack(output_groups, metadata_conflicts="silent")
        output.meta = trace_table.data.meta.copy()
    else:
        output = trace_table.data[:0].copy()
    # Spot_ID is copied from source rows and is never generated or transformed.
    trace_table.data = output
    diagnostic_table = Table(rows=diagnostics)
    diagnostic_table.meta["comments"] = [
        "confidence_score is a normalized score gap, not a probability",
        f"history_mode={history_mode}",
        f"distance_score={distance_score}",
        f"rejection_cost={rejection_cost}",
    ]
    diagnostic_table.meta["genomic_source"] = model.genomic_source
    diagnostic_table.meta["fallback_source"] = model.fallback_source
    return diagnostic_table


def main():
    print_script_banner(__file__, __doc__)
    args = parse_arguments().parse_args()
    trace_files = create_dict_args(args)["trace_files"]
    if not trace_files or trace_files == [None]:
        print(
            "! Error: did not find any trace file to analyze. Please provide one using --input or --pipe."
        )
        return
    print(
        f"\n{len(trace_files)} trace files to process= {' '.join(map(str, trace_files))}"
    )
    for trace_file in trace_files:
        output = args.output or f"{os.path.splitext(trace_file)[0]}_split.ecsv"
        trace_table = ChromatinTraceTable()
        trace_table.load(trace_file)
        if args.clustering_method == "beam":
            selection = "traces selected by barcode multiplicity"
        else:
            selection = (
                "all traces"
                if args.split_all
                else f"traces with Rg > mean + {args.std_threshold} * std_dev"
            )
        print(f"Applying {args.clustering_method} resolution on {selection}...")
        if args.clustering_method == "beam":
            thresholds = ClassificationThresholds(
                args.split_min_repeated_barcodes,
                args.split_min_repeated_fraction,
                args.candidate_rule,
            )
            diagnostics = resolve_traces(
                trace_table,
                thresholds=thresholds,
                history_mode=args.history_mode,
                distance_score=args.distance_score,
                history_length=args.history_length,
                beam_width=args.beam_width,
                rejection_cost=args.rejection_cost,
                minimum_polymer_size=args.minimum_polymer_size,
                minimum_confidence=args.minimum_confidence,
                model_minimum_observations=args.model_min_observations,
                variance_floor=args.variance_floor,
            )
            diagnostics_output = args.diagnostics_output or (
                f"{os.path.splitext(trace_file)[0]}_split_diagnostics.ecsv"
            )
            diagnostics.write(diagnostics_output, format="ascii.ecsv", overwrite=True)
        else:
            split_large_traces(
                trace_table,
                args.std_threshold,
                args.num_clusters,
                args.split_all,
                args.clustering_method,
                args.min_cluster_size,
                args.min_samples,
                args.cluster_selection_method,
                args.cluster_selection_epsilon,
                args.allow_single_cluster,
            )
        trace_table.save(output)


if __name__ == "__main__":
    main()
