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
import select
import sys

import matplotlib.pyplot as plt
import numpy as np

# Removed seaborn import
from tqdm import tqdm

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


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
    # Dictionary to store co-localization counts: (barcode1, barcode2) -> [co-localized count, total count]
    threeway_interactions = {}

    # Get all unique barcodes
    all_barcodes = np.unique(trace_table["Barcode #"])
    other_barcodes = [b for b in all_barcodes if b != anchor_barcode]

    # Generate all possible pairs of barcodes (excluding the anchor)
    barcode_pairs = list(itertools.combinations(other_barcodes, 2))

    # Initialize the counts dictionary
    for pair in barcode_pairs:
        threeway_interactions[pair] = [0, 0]

    # Group traces by Trace_ID
    trace_groups = trace_table.group_by("Trace_ID").groups

    # Process each trace separately
    for trace in tqdm(trace_groups, desc="Processing traces"):
        # Get positions of the anchor barcode in this trace
        anchor_positions = trace[trace["Barcode #"] == anchor_barcode]

        # Skip if anchor is not present in this trace
        if len(anchor_positions) == 0:
            continue

        # Get the barcodes present in this trace
        barcodes_in_trace = np.unique(trace["Barcode #"])

        # Generate all pairs of barcodes in this trace (excluding the anchor)
        pairs_in_trace = [
            p
            for p in barcode_pairs
            if p[0] in barcodes_in_trace and p[1] in barcodes_in_trace
        ]

        # For each pair, check if they co-localize with the anchor
        for barcode1, barcode2 in pairs_in_trace:
            barcode1_positions = trace[trace["Barcode #"] == barcode1]
            barcode2_positions = trace[trace["Barcode #"] == barcode2]

            # Skip if either barcode is missing in this trace
            if len(barcode1_positions) == 0 or len(barcode2_positions) == 0:
                continue

            # Increment the total count for this pair
            threeway_interactions[(barcode1, barcode2)][1] += 1

            # Check distances between anchor and barcode1
            distances_anchor_barcode1 = np.linalg.norm(
                np.array(
                    [
                        anchor_positions["x"],
                        anchor_positions["y"],
                        anchor_positions["z"],
                    ]
                ).T[:, None]
                - np.array(
                    [
                        barcode1_positions["x"],
                        barcode1_positions["y"],
                        barcode1_positions["z"],
                    ]
                ).T,
                axis=-1,
            )

            # Check distances between anchor and barcode2
            distances_anchor_barcode2 = np.linalg.norm(
                np.array(
                    [
                        anchor_positions["x"],
                        anchor_positions["y"],
                        anchor_positions["z"],
                    ]
                ).T[:, None]
                - np.array(
                    [
                        barcode2_positions["x"],
                        barcode2_positions["y"],
                        barcode2_positions["z"],
                    ]
                ).T,
                axis=-1,
            )

            # Check if both barcodes co-localize with the anchor
            anchor_barcode1_coloc = np.any(distances_anchor_barcode1 < distance_cutoff)
            anchor_barcode2_coloc = np.any(distances_anchor_barcode2 < distance_cutoff)

            # If both barcodes co-localize with the anchor, increment the co-localized count
            if anchor_barcode1_coloc and anchor_barcode2_coloc:
                threeway_interactions[(barcode1, barcode2)][0] += 1

    # Compute frequencies
    threeway_frequencies = {
        pair: (count[0] / count[1] if count[1] > 0 else 0)
        for pair, count in threeway_interactions.items()
    }

    return threeway_frequencies


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
    # Dictionary to store bootstrap samples for each barcode pair
    pair_samples = {}

    # Get all unique trace IDs for bootstrapping
    trace_ids = np.unique(trace_table["Trace_ID"])

    # Run bootstrap iterations
    for _ in tqdm(range(n_bootstrap), desc="Bootstrapping"):
        # Sample traces with replacement
        sampled_traces = np.random.choice(trace_ids, size=len(trace_ids), replace=True)

        # Create a new table with only the sampled traces
        sampled_table = trace_table[np.isin(trace_table["Trace_ID"], sampled_traces)]

        # Compute three-way co-localization for this bootstrap sample
        threeway_frequencies = compute_threeway_colocalization(
            sampled_table, anchor_barcode, distance_cutoff
        )

        # Store the results
        for pair, frequency in threeway_frequencies.items():
            if pair not in pair_samples:
                pair_samples[pair] = []
            pair_samples[pair].append(frequency)

    # Compute mean and standard error of the mean (SEM) for each pair
    pair_means = {pair: np.mean(samples) for pair, samples in pair_samples.items()}

    pair_sems = {
        pair: np.std(samples) / np.sqrt(n_bootstrap)
        for pair, samples in pair_samples.items()
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
    ax.set_xticklabels(sorted_barcodes, fontsize=10)
    ax.set_yticklabels(sorted_barcodes, fontsize=10)

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)

    # Add perpendicular lines for the anchor barcode
    if anchor_barcode in barcode_to_idx:
        anchor_idx = barcode_to_idx[anchor_barcode]

        # Horizontal line across the anchor barcode row
        ax.axhline(y=anchor_idx, color="black", linestyle="-", linewidth=2, alpha=0.7)

        # Vertical line across the anchor barcode column
        ax.axvline(x=anchor_idx, color="black", linestyle="-", linewidth=2, alpha=0.7)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Co-localization frequency", fontsize=12)

    # Add title and labels
    ax.set_title(
        f"3-way co-localization with anchor {anchor_barcode}\n(distance cutoff: {distance_cutoff} µm)",
        fontsize=14,
    )
    ax.set_xlabel("Barcode #", fontsize=14)
    ax.set_ylabel("Barcode #", fontsize=14)

    # Adjust layout and saves npy matrix and image
    plt.tight_layout()
    output_filename = f"{output_file.split('.')[0]}_anchor_{anchor_barcode}"

    np.save(f"{output_filename}.npy", mean_matrix)

    plt.savefig(f"{output_filename}.png", dpi=300)
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
    ax.set_xticklabels(sorted_barcodes, fontsize=14)
    ax.set_yticklabels(sorted_barcodes, fontsize=14)

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(sorted_barcodes), 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Standard error", fontsize=12)

    # Add title and labels
    ax.set_title(
        f"Errors for 3-way co-localization with anchor {anchor_barcode}\n"
        f"(distance cutoff: {distance_cutoff} µm)",
        fontsize=15,
    )
    ax.set_xlabel("Barcode #", fontsize=14)
    ax.set_ylabel("Barcode #", fontsize=14)

    # Adjust layout and save
    plt.tight_layout()
    sem_output_filename = f"{output_file.split('.')[0]}_anchor_{anchor_barcode}_sem.png"
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
    parser = parse_arguments()
    args = parser.parse_args()
    trace_files = get_trace_files(args)

    if len(trace_files) > 0:
        for trace_file in trace_files:
            print(f"Processing file: {trace_file}")

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
                plot_threeway_matrix(
                    pair_means,
                    pair_sems,
                    anchor,
                    args.output,
                    distance_cutoff=args.cutoff,
                    vmin=args.vmin,
                    vmax=args.vmax,
                )

    else:
        print("\nNo trace files were detected")


if __name__ == "__main__":
    main()
