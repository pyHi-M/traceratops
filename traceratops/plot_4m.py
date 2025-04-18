#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute barcode interaction frequencies with bootstrapping.

A script for performing 4M (Multi-way Measurement of Molecular interactions in space) analysis
on chromatin trace data. This tool analyzes spatial colocalization between a specified anchor
barcode and all other barcodes in 3D chromatin trace data.

The script calculates colocalization frequencies based on a distance cutoff and performs
bootstrapping to estimate statistical confidence (mean and standard error). It generates
a plot showing the frequency of interaction between the anchor barcode and all others.

This is particularly useful for analyzing chromatin organization, DNA-DNA interactions,
and spatial proximity relationships in microscopy data.
"""
import argparse
import select
import sys

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


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
        help="Distance cutoff for colocalization. Default = 0.2 um",
    )
    parser.add_argument(
        "--bootstrapping_cycles",
        type=int,
        default=10,
        help="Number of bootstrap iterations.",
    )
    parser.add_argument(
        "--output", default="colocalization_plot.png", help="Output file for the plot."
    )
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument("--x_min", type=int, default=None, help="xscale minimum")
    parser.add_argument("--x_max", type=int, default=None, help="xscale maximum")
    return parser


def compute_colocalization(trace_table, anchor_barcode, distance_cutoff):
    """Computes the frequency of colocalization between the anchor barcode and all other barcodes."""
    # Make sure we're dealing with a single anchor barcode
    if isinstance(anchor_barcode, list):
        raise TypeError(
            "compute_colocalization expects a single anchor_barcode, not a list"
        )
    barcode_interactions = {}
    trace_groups = trace_table.group_by("Trace_ID").groups
    for trace in tqdm(trace_groups, desc="Processing traces"):
        anchor_positions = trace[trace["Barcode #"] == anchor_barcode]
        other_barcodes = np.unique(trace["Barcode #"])
        for barcode in other_barcodes:
            if barcode == anchor_barcode:
                continue
            other_positions = trace[trace["Barcode #"] == barcode]
            if len(anchor_positions) > 0 and len(other_positions) > 0:
                distances = np.linalg.norm(
                    np.array(
                        [
                            anchor_positions["x"],
                            anchor_positions["y"],
                            anchor_positions["z"],
                        ]
                    ).T[:, None]
                    - np.array(
                        [
                            other_positions["x"],
                            other_positions["y"],
                            other_positions["z"],
                        ]
                    ).T,
                    axis=-1,
                )
                colocalized = np.any(distances < distance_cutoff)
                if barcode not in barcode_interactions:
                    barcode_interactions[barcode] = [0, 0]
                barcode_interactions[barcode][1] += 1  # Total occurrences
                if colocalized:
                    barcode_interactions[barcode][0] += 1  # Co-localized occurrences
    # Compute frequencies
    barcode_frequencies = {
        barcode: (count[0] / count[1] if count[1] > 0 else 0)
        for barcode, count in barcode_interactions.items()
    }
    barcode_frequencies[anchor_barcode] = 0.0
    return barcode_frequencies


def bootstrap_colocalization(
    trace_table, anchor_barcodes, distance_cutoff, n_bootstrap=100
):
    """Performs bootstrapping to estimate mean and SEM of colocalization frequencies for multiple anchors."""
    # If a single anchor is provided, convert to a list for consistent processing
    if not isinstance(anchor_barcodes, list):
        anchor_barcodes = [anchor_barcodes]
    barcode_samples = {anchor: {} for anchor in anchor_barcodes}
    trace_ids = np.unique(trace_table["Trace_ID"])
    for _ in tqdm(range(n_bootstrap), desc="Bootstrapping"):
        sampled_traces = np.random.choice(trace_ids, size=len(trace_ids), replace=True)
        sampled_table = trace_table[np.isin(trace_table["Trace_ID"], sampled_traces)]
        # Process each anchor barcode separately
        for anchor in anchor_barcodes:
            colocalization = compute_colocalization(
                sampled_table, anchor, distance_cutoff
            )
            for barcode, frequency in colocalization.items():
                if barcode not in barcode_samples[anchor]:
                    barcode_samples[anchor][barcode] = []
                barcode_samples[anchor][barcode].append(frequency)
    # Calculate means and SEMs for each anchor-barcode pair
    barcode_means = {
        anchor: {
            barcode: np.mean(values)
            for barcode, values in barcode_samples[anchor].items()
        }
        for anchor in anchor_barcodes
    }
    barcode_sems = {
        anchor: {
            barcode: np.std(values) / np.sqrt(n_bootstrap)
            for barcode, values in barcode_samples[anchor].items()
        }
        for anchor in anchor_barcodes
    }
    # Get unique barcodes for each anchor
    unique_barcodes = {}
    for anchor in anchor_barcodes:
        # Make sure all anchors have the same set of barcodes for consistent plotting
        unique_barcodes[anchor] = sorted(barcode_samples[anchor].keys())
    # Create consistent vectors for plotting (with same length for all anchors)
    mean_frequencies = {}
    sem_frequencies = {}
    for anchor in anchor_barcodes:
        mean_frequencies[anchor] = [
            barcode_means[anchor].get(b, 0) for b in unique_barcodes[anchor]
        ]
        sem_frequencies[anchor] = [
            barcode_sems[anchor].get(b, 0) for b in unique_barcodes[anchor]
        ]
    return barcode_means, barcode_sems


def plot_frequencies(
    mean_frequencies, sem_frequencies, anchor_barcodes, output_file, x_min=0, x_max=0
):
    """Plots colocalization frequencies for multiple anchors separately."""
    # Make sure anchor_barcodes is a list for consistent processing
    if not isinstance(anchor_barcodes, list):
        anchor_barcodes = [anchor_barcodes]
    print(f"\n$ Processing anchors: {anchor_barcodes}")
    for anchor in anchor_barcodes:
        plt.figure(figsize=(8, 4))  # Reduce figure size
        # Get the barcodes for this anchor
        barcodes = sorted(mean_frequencies[anchor].keys())
        mean_values = [mean_frequencies[anchor][b] for b in barcodes]
        sem_values = [sem_frequencies[anchor][b] for b in barcodes]
        plt.errorbar(
            barcodes,
            mean_values,
            yerr=sem_values,
            fmt="o-",
            capsize=5,
            label=f"Anchor {anchor}",
        )
        plt.axvline(
            anchor, color="red", linestyle="--", label=f"Anchor {anchor}"
        )  # Highlight anchor
        plt.xlabel("Barcode #", fontsize=13)
        plt.ylabel("Colocalization frequency", fontsize=13)
        plt.title(f"4M plot for anchor: {str(anchor)}", fontsize=15)
        plt.xticks(fontsize=10, rotation=90)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=8)
        # Set x-axis limits
        _x_min = x_min if x_min is not None else min(barcodes) - 1
        _x_max = x_max if x_max is not None else max(barcodes) + 1
        plt.xlim(_x_min, _x_max)
        plt.grid(True)
        output_filename = f"{output_file.split('.')[0]}_anchor_{anchor}.png"
        plt.savefig(output_filename)
        plt.close()  # Close the figure to avoid memory issues with many anchors


def main():
    parser = parse_arguments()
    args = parser.parse_args()

    trace_files = []
    if args.pipe:
        if select.select(
            [
                sys.stdin,
            ],
            [],
            [],
            0.0,
        )[0]:
            trace_files = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        trace_files = [args.input]

    if len(trace_files) > 0:

        for trace_file in trace_files:
            print(f"\n$ Processing trace file: {trace_file}")

            trace = ChromatinTraceTable()
            trace.initialize()

            trace.load(trace_file)
            barcode_means, barcode_sems = bootstrap_colocalization(
                trace.data, args.anchors, args.cutoff, args.bootstrapping_cycles
            )
            plot_frequencies(
                barcode_means,
                barcode_sems,
                args.anchors,
                args.output,
                x_min=args.x_min,
                x_max=args.x_max,
            )

    else:
        print("\nNo trace files were detected")


if __name__ == "__main__":
    main()
