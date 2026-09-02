#!/usr/bin/env python3
"""Split unusually large traces, or every trace, into spatial clusters.

By default, traces whose radius of gyration is greater than the population mean
plus one standard deviation are split with K-means.  ``--split-all`` disables
that size filter.  HDBSCAN can be selected as an alternative and its distance
epsilon can either be supplied or estimated from the input traces.
"""

import argparse
import os
import select
import sys
import uuid

import numpy as np
from scipy.spatial.distance import pdist
from sklearn.cluster import HDBSCAN, KMeans

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
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
        "--clustering-method",
        choices=("kmeans", "hdbscan"),
        default="kmeans",
        help="Clustering algorithm (default: kmeans).",
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
        coords = np.column_stack((trace["x"], trace["y"], trace["z"]))
        if len(coords) > 1:
            thresholds.append(_otsu_threshold(pdist(coords)))
    return float(np.median(thresholds)) if thresholds else 0.0


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
    rg_values = [
        compute_radius_of_gyration(np.column_stack((t["x"], t["y"], t["z"])))
        for t in groups
    ]
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
        coords = np.column_stack((trace["x"], trace["y"], trace["z"]))
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
            labels = HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                metric="euclidean",
                cluster_selection_method=cluster_selection_method,
                cluster_selection_epsilon=epsilon,
                allow_single_cluster=allow_single_cluster,
            ).fit_predict(coords)

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
        selection = (
            "all traces"
            if args.split_all
            else f"traces with Rg > mean + {args.std_threshold} * std_dev"
        )
        print(f"Applying {args.clustering_method} clustering on {selection}...")
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
