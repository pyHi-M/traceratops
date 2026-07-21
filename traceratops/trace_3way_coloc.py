#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute three-way co-localization frequencies with bootstrapping.

A script for performing three-way co-localization analysis in 4M (Multi-way Measurement of
Molecular interactions in space) data. This tool analyzes spatial co-localization between
a specified anchor barcode and all possible pairs of other barcodes in 3D chromatin trace data.

The script calculates three-way co-localization frequencies based on a distance cutoff and
performs bootstrapping to estimate statistical confidence (mean and standard error). It
generates a heatmap showing the frequency of interaction between the anchor and all possible
pairs of other barcodes.

This is particularly useful for analyzing higher-order chromatin organization and complex
spatial relationships in microscopy data.
"""

import argparse
import itertools
import os
import select
import sys

import matplotlib.pyplot as plt
import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.script_banner import print_script_banner


def _matrix_tick_fontsize(n_barcodes):
    """Return a readable tick-label size for dense barcode matrices."""
    if n_barcodes <= 0:
        return 8
    return max(4, min(8, int(450 / n_barcodes)))


def _anchor_boundary_index(sorted_barcodes, anchor_barcode):
    """Return the matrix boundary where an omitted anchor barcode would sort."""
    return np.searchsorted(np.asarray(sorted_barcodes), anchor_barcode) - 0.5


def _trace_barcode_arrays(trace_table, anchor_barcode, distance_cutoff):
    """Precompute per-trace barcode presence and anchor co-localization arrays."""
    all_barcodes = np.sort(np.unique(trace_table["Barcode #"]))
    other_barcodes = all_barcodes[all_barcodes != anchor_barcode]
    n_barcodes = len(other_barcodes)
    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(other_barcodes)}

    trace_groups = trace_table.group_by("Trace_ID").groups
    n_traces = len(trace_groups)
    present = np.zeros((n_traces, n_barcodes), dtype=bool)
    colocated = np.zeros((n_traces, n_barcodes), dtype=bool)

    cutoff_squared = distance_cutoff * distance_cutoff

    for trace_idx, trace in enumerate(trace_groups):
        trace_barcodes = np.asarray(trace["Barcode #"])
        anchor_mask = trace_barcodes == anchor_barcode

        # The original denominator only counts traces where the anchor and both
        # partner barcodes are present.
        if not np.any(anchor_mask):
            continue

        unique_barcodes = np.unique(trace_barcodes)
        partner_barcodes = unique_barcodes[unique_barcodes != anchor_barcode]
        partner_indices = [
            barcode_to_idx[b] for b in partner_barcodes if b in barcode_to_idx
        ]
        present[trace_idx, partner_indices] = True

        anchor_xyz = np.column_stack(
            (
                np.asarray(trace["x"])[anchor_mask],
                np.asarray(trace["y"])[anchor_mask],
                np.asarray(trace["z"])[anchor_mask],
            )
        )

        for barcode in partner_barcodes:
            barcode_idx = barcode_to_idx.get(barcode)
            if barcode_idx is None:
                continue

            barcode_mask = trace_barcodes == barcode
            barcode_xyz = np.column_stack(
                (
                    np.asarray(trace["x"])[barcode_mask],
                    np.asarray(trace["y"])[barcode_mask],
                    np.asarray(trace["z"])[barcode_mask],
                )
            )
            deltas = anchor_xyz[:, None, :] - barcode_xyz[None, :, :]
            distances_squared = np.einsum("ijk,ijk->ij", deltas, deltas)
            colocated[trace_idx, barcode_idx] = np.any(
                distances_squared < cutoff_squared
            )

    return other_barcodes, present, colocated


def _weighted_pair_counts(present, colocated, weights=None):
    """Return pair denominators and co-localized counts for optional trace weights."""
    if weights is None:
        present_values = present.astype(np.int64)
        colocated_values = colocated.astype(np.int64)
    else:
        weights = np.asarray(weights, dtype=np.int64)
        present_values = present.astype(np.int64) * weights[:, None]
        colocated_values = colocated.astype(np.int64) * weights[:, None]

    total_counts = present_values.T @ present.astype(np.int64)
    colocated_counts = colocated_values.T @ colocated.astype(np.int64)
    return total_counts, colocated_counts


def _pair_frequencies_from_counts(barcodes, total_counts, colocated_counts):
    """Convert pair count matrices into the public dictionary representation."""
    frequencies = {}
    for i, barcode1 in enumerate(barcodes):
        for j in range(i + 1, len(barcodes)):
            total = total_counts[i, j]
            barcode2 = barcodes[j]
            frequencies[(barcode1, barcode2)] = (
                colocated_counts[i, j] / total if total > 0 else 0
            )
    return frequencies


def compute_threeway_colocalization(trace_table, anchor_barcode, distance_cutoff):
    """
    Computes the frequency of three-way co-localization between an anchor barcode
    and all possible pairs of other barcodes.

    Parameters:
    ----------
    trace_table : ChromatinTraceTable
        Table containing chromatin trace data
    anchor_barcode : int
        The anchor barcode number
    distance_cutoff : float
        Distance threshold for considering barcodes as co-localized (in µm)

    Returns:
    -------
    dict
        Dictionary with (barcode1, barcode2) tuples as keys and co-localization frequencies as values
    """
    barcodes, present, colocated = _trace_barcode_arrays(
        trace_table, anchor_barcode, distance_cutoff
    )
    total_counts, colocated_counts = _weighted_pair_counts(present, colocated)
    return _pair_frequencies_from_counts(barcodes, total_counts, colocated_counts)


def bootstrap_threeway_colocalization(
    trace_table, anchor_barcode, distance_cutoff, n_bootstrap=100
):
    """
    Performs bootstrapping to estimate mean and SEM of three-way co-localization frequencies.

    Parameters:
    ----------
    trace_table : ChromatinTraceTable
        Table containing chromatin trace data
    anchor_barcode : int
        The anchor barcode number
    distance_cutoff : float
        Distance threshold for considering barcodes as co-localized (in µm)
    n_bootstrap : int
        Number of bootstrap iterations

    Returns:
    -------
    tuple
        (mean_frequencies, sem_frequencies) dictionaries with barcode pairs as keys
    """
    barcodes, present, colocated = _trace_barcode_arrays(
        trace_table, anchor_barcode, distance_cutoff
    )
    n_traces = present.shape[0]
    n_barcodes = len(barcodes)
    n_pairs = n_barcodes * (n_barcodes - 1) // 2

    if n_bootstrap <= 0 or n_pairs == 0 or n_traces == 0:
        empty_frequencies = _pair_frequencies_from_counts(
            barcodes,
            np.zeros((n_barcodes, n_barcodes), dtype=np.int64),
            np.zeros((n_barcodes, n_barcodes), dtype=np.int64),
        )
        return empty_frequencies, empty_frequencies.copy()

    samples = {pair: [] for pair in itertools.combinations(barcodes, 2)}

    # Resampling via trace weights preserves duplicate draws while avoiding the
    # expensive reconstruction of an Astropy table for every bootstrap cycle.
    for _ in range(n_bootstrap):
        sampled_indices = np.random.randint(0, n_traces, size=n_traces)
        weights = np.bincount(sampled_indices, minlength=n_traces)
        total_counts, colocated_counts = _weighted_pair_counts(
            present, colocated, weights=weights
        )
        frequencies = _pair_frequencies_from_counts(
            barcodes, total_counts, colocated_counts
        )
        for pair, frequency in frequencies.items():
            samples[pair].append(frequency)

    pair_means = {
        pair: np.mean(pair_samples) for pair, pair_samples in samples.items()
    }
    pair_sems = {
        pair: np.std(pair_samples) / np.sqrt(n_bootstrap)
        for pair, pair_samples in samples.items()
    }

    return pair_means, pair_sems

def plot_threeway_matrix(
    pair_means,
    pair_sems,
    anchor_barcode,
    output_file,
    distance_cutoff=0.2,
    vmin=None,
    vmax=None,
    output_format="png",
):
    """
    Creates a heatmap of three-way co-localization frequencies using matplotlib.

    Parameters:
    ----------
    pair_means : dict
        Dictionary with (barcode1, barcode2) tuples as keys and mean frequencies as values
    pair_sems : dict
        Dictionary with (barcode1, barcode2) tuples as keys and SEM values as values
    anchor_barcode : int
        The anchor barcode number
    output_file : str
        Output file name for the plot
    """
    # Get all unique barcodes from the pairs
    all_barcodes = set()
    for b1, b2 in pair_means.keys():
        all_barcodes.add(b1)
        all_barcodes.add(b2)

    # Sort barcodes for consistent matrix indexing
    sorted_barcodes = sorted(all_barcodes)
    n_barcodes = len(sorted_barcodes)

    # Create empty matrices for means and SEMs
    mean_matrix = np.zeros((n_barcodes, n_barcodes))
    sem_matrix = np.zeros((n_barcodes, n_barcodes))

    # Create a mapping from barcode to matrix index
    barcode_to_idx = {b: i for i, b in enumerate(sorted_barcodes)}

    # Fill the matrices with the computed values
    for (b1, b2), mean_val in pair_means.items():
        i, j = barcode_to_idx[b1], barcode_to_idx[b2]
        mean_matrix[i, j] = mean_val
        mean_matrix[j, i] = mean_val  # Mirror the matrix (symmetric)

        sem_val = pair_sems[(b1, b2)]
        sem_matrix[i, j] = sem_val
        sem_matrix[j, i] = sem_val  # Mirror the matrix (symmetric)

    # Create a custom colormap from white to dark blue
    cmap = "RdBu"  # LinearSegmentedColormap.from_list('white_to_blue', ['#FFFFFF', '#0343DF'])

    tick_fontsize = _matrix_tick_fontsize(n_barcodes)
    label_fontsize = 11
    title_fontsize = 12
    colorbar_tick_fontsize = 9
    colorbar_label_fontsize = 10

    # Create the figure and subplots for the mean frequencies
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the mean heatmap using matplotlib
    im = ax.imshow(
        mean_matrix,
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin if vmin is not None else 0,
        vmax=(
            vmax
            if vmax is not None
            else (0.9 * mean_matrix.max() if mean_matrix.max() > 0 else 1)
        ),
    )

    # Set up the axes with the correct labels
    ax.set_xticks(np.arange(len(sorted_barcodes)))
    ax.set_yticks(np.arange(len(sorted_barcodes)))
    ax.set_xticklabels(sorted_barcodes, fontsize=tick_fontsize)
    ax.set_yticklabels(sorted_barcodes, fontsize=tick_fontsize)

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)

    # Add perpendicular guide lines where the omitted anchor barcode would sort.
    # The anchor is excluded from the matrix, so draw on the boundary between
    # neighboring barcode cells rather than looking for the anchor in the labels.
    anchor_boundary_idx = _anchor_boundary_index(sorted_barcodes, anchor_barcode)
    ax.axhline(
        y=anchor_boundary_idx, color="black", linestyle="-", linewidth=1.5, alpha=0.7
    )
    ax.axvline(
        x=anchor_boundary_idx, color="black", linestyle="-", linewidth=1.5, alpha=0.7
    )

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Co-localization frequency", fontsize=colorbar_label_fontsize)
    cbar.ax.tick_params(labelsize=colorbar_tick_fontsize)

    # Add title and labels
    ax.set_title(
        f"3-way co-localization with anchor {anchor_barcode}\n(distance cutoff: {distance_cutoff} µm)",
        fontsize=title_fontsize,
    )
    ax.set_xlabel("Barcode #", fontsize=label_fontsize)
    ax.set_ylabel("Barcode #", fontsize=label_fontsize)

    # Adjust layout and saves npy matrix and image
    plt.tight_layout()
    output_filename = f"{output_file.split('.')[0]}_anchor_{anchor_barcode}"

    np.save(f"{output_filename}.npy", mean_matrix)

    plt.savefig(f"{output_filename}.{output_format}", dpi=300)
    print(f"Saved three-way co-localization heatmap to: {output_filename}")
    plt.close()

    # Create the figure and subplots for the standard errors
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the SEM heatmap using matplotlib
    im = ax.imshow(
        np.abs(sem_matrix),
        interpolation="nearest",
        cmap="YlOrRd",
        vmin=0,
        vmax=np.abs(sem_matrix).max() if np.abs(sem_matrix).max() > 0 else 0.3,
    )

    # Set up the axes with the correct labels
    ax.set_xticks(np.arange(len(sorted_barcodes)))
    ax.set_yticks(np.arange(len(sorted_barcodes)))
    ax.set_xticklabels(sorted_barcodes, fontsize=tick_fontsize)
    ax.set_yticklabels(sorted_barcodes, fontsize=tick_fontsize)

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Standard error", fontsize=colorbar_label_fontsize)
    cbar.ax.tick_params(labelsize=colorbar_tick_fontsize)

    # Add title and labels
    ax.set_title(
        f"Errors for 3-way co-localization with anchor {anchor_barcode}\n"
        f"(distance cutoff: {distance_cutoff} µm)",
        fontsize=title_fontsize,
    )
    ax.set_xlabel("Barcode #", fontsize=label_fontsize)
    ax.set_ylabel("Barcode #", fontsize=label_fontsize)

    # Adjust layout and save
    plt.tight_layout()
    sem_output_filename = (
        f"{output_file.rsplit('.', 1)[0]}_anchor_{anchor_barcode}_sem.{output_format}"
    )
    plt.savefig(sem_output_filename, dpi=300)
    print(f"Saved SEM heatmap to: {sem_output_filename}")
    plt.close()


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, help="Path to input trace table (ECSV format)."
    )
    parser.add_argument(
        "--anchors",
        type=int,
        nargs="+",
        required=True,
        help="List of anchor barcode numbers.",
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        required=False,
        default=0.2,
        help="Distance cutoff for co-localization. Default = 0.2 um",
    )

    parser.add_argument(
        "--vmin", type=float, default=None, help="Minimum value for colormap scale."
    )
    parser.add_argument(
        "--vmax", type=float, default=None, help="Maximum value for colormap scale."
    )

    parser.add_argument(
        "--bootstrapping_cycles",
        type=int,
        default=10,
        help="Number of bootstrap iterations.",
    )
    parser.add_argument(
        "--output", default="threeway_coloc_plot.png", help="Output file for the plot."
    )
    parser.add_argument(
        "--output_format",
        choices=["png", "svg", "pdf"],
        default="png",
        help="Output image format. Default = png.",
    )
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )

    return parser


def get_trace_files(args):
    trace_files = []
    if args.pipe:
        if select.select([sys.stdin], [], [], 0.0)[0]:
            trace_files = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        trace_files = [args.input]
    return args, trace_files


def main():
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    args = parser.parse_args()
    _, trace_files = get_trace_files(args)

    if len(trace_files) > 0:
        for trace_file in trace_files:
            print(f">> Processing file: {trace_file}")

            # Initialize and load trace table
            trace = ChromatinTraceTable()
            trace.initialize()
            trace.load(trace_file)

            print(f"Using distance cutoff: {args.cutoff} µm")
            print(f"Performing {args.bootstrapping_cycles} bootstrap iterations")

            for anchor in args.anchors:
                print(f"\nRunning analysis for anchor: {anchor}")

                # Run the bootstrap analysis
                pair_means, pair_sems = bootstrap_threeway_colocalization(
                    trace.data,
                    anchor,
                    args.cutoff,
                    n_bootstrap=args.bootstrapping_cycles,
                )

                # Create the plots
                trace_file_basename = os.path.basename(trace_file)
                output_filename = f"{args.output.rsplit('.', 1)[0]}_{trace_file_basename.rsplit('.', 1)[0]}"
                plot_threeway_matrix(
                    pair_means,
                    pair_sems,
                    anchor,
                    output_filename,
                    distance_cutoff=args.cutoff,
                    vmin=args.vmin,
                    vmax=args.vmax,
                    output_format=args.output_format,
                )

    else:
        print("\nNo trace files were detected")


if __name__ == "__main__":
    main()
