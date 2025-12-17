#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze chromatin trace files.
"""

import argparse
import collections
import select
import sys
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity

from traceratops.core.chromatin_trace_table import ChromatinTraceTable

font = {"weight": "normal", "size": 22}
matplotlib.rc("font", **font)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-F", "--rootFolder", default=".", help="Folder with images")
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument(
        "--plotXYZ",
        default=False,
        help="Plots XYZ traces for all ROIs",
        action="store_true",
    )
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument(
        "--format",
        default="png",
        choices=["png", "svg"],
        help="Output image format (png or svg)",
    )
    parser.add_argument(
        "--rgMin",
        type=float,
        default=0.02,
        help="Minimum Rg value to include in the analysis (default: 0.02).",
    )
    parser.add_argument(
        "--rgBarcodes",
        type=int,
        nargs="+",
        help="List of barcode IDs to include when computing the radius of gyration.",
    )
    return parser


def create_dict_args(args):
    p = {}
    p["input"] = args.input
    p["rootFolder"] = args.rootFolder
    p["plotXYZ"] = args.plotXYZ
    p["format"] = args.format
    p["rgMin"] = args.rgMin
    p["rgBarcodes"] = args.rgBarcodes

    p["trace_files"] = []
    if args.pipe:
        p["pipe"] = True
        if select.select(
            [
                sys.stdin,
            ],
            [],
            [],
            0.0,
        )[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False
        p["trace_files"] = [p["input"]]

    return p


def get_barcode_statistics(trace, output_filename="test_barcodes.png"):
    """
    Calculate and visualize barcode statistics from trace data.

    The function calculates:
        - Number of barcodes per trace
        - Number of unique barcodes per trace
        - Number of repeated barcodes per trace

    Parameters
    ----------
    trace : astropy.table.Table
        Trace table in ASTROPY Table format.
    output_filename : str
        Output figure filename including path and extension.

    Returns
    -------
    None
        The function saves the output figure but does not return any values.
    """
    trace_by_ID = trace.group_by("Trace_ID")

    trace_lengths = list()
    trace_unique_barcodes = list()
    trace_repeated_barcodes = list()
    number_unique_barcodes = list()
    number_repeated_barcodes = list()

    for sub_trace_table in trace_by_ID.groups:
        trace_lengths.append(len(sub_trace_table))

        unique_barcodes = list(set(sub_trace_table["Barcode #"]))
        trace_unique_barcodes.append(unique_barcodes)
        number_unique_barcodes.append(len(unique_barcodes))

        repeated_barcodes = [
            item
            for item, count in collections.Counter(sub_trace_table["Barcode #"]).items()
            if count > 1
        ]
        trace_repeated_barcodes.append(repeated_barcodes)
        number_repeated_barcodes.append(len(repeated_barcodes))

    distributions = [trace_lengths, number_unique_barcodes, number_repeated_barcodes]
    axis_x_labels = [
        "$N_{barcodes}$",
        "$N_{unique-barcodes}$",
        "$N_{repeated-barcodes}$",
    ]
    number_plots = len(distributions)

    fig = plt.figure(constrained_layout=True)
    im_size = 8
    fig.set_size_inches((im_size * number_plots, im_size))
    gs = fig.add_gridspec(1, number_plots)
    axes = [fig.add_subplot(gs[0, i]) for i in range(number_plots)]

    for axis, distribution, xlabel in zip(axes, distributions, axis_x_labels):
        axis.hist(distribution, alpha=0.3)
        axis.set_xlabel(xlabel, fontsize=30)
        axis.set_ylabel("counts", fontsize=30)
        axis.set_title(
            f"n = {str(len(distribution))} | median = {str(np.median(distribution))}",
            fontsize=20,
        )

    fig.suptitle("Trace statistics", fontsize=40)

    plt.savefig(output_filename)


def plot_neighbor_distances(trace, output_filename="neighbor_distances.png"):
    """
    Calculate and visualize distances between consecutive neighboring barcodes.

    This function computes the mean and standard deviation of X, Y, and Z distances
    between strictly consecutive neighboring barcodes and generates histograms for each dimension.

    Parameters
    ----------
    trace : ChromatinTraceTable
        Trace table, instance of the ChromatinTraceTable Class.
    output_filename : str
        The filename for the output figure file.

    Returns
    -------
    tuple
        (mean_dx, mean_dy, mean_dz, std_dx, std_dy, std_dz):
        Mean and standard deviation values for X, Y, and Z distances.
    """
    trace_table = trace.data
    trace_by_ID = trace_table.group_by("Trace_ID")

    dx_all, dy_all, dz_all = [], [], []

    for sub_trace_table in trace_by_ID.groups:
        # Sort by Barcode #
        sorted_trace = sub_trace_table[np.argsort(sub_trace_table["Barcode #"])]

        # Get barcodes and coordinates
        barcodes = sorted_trace["Barcode #"].data
        x_coords = sorted_trace["x"].data
        y_coords = sorted_trace["y"].data
        z_coords = sorted_trace["z"].data

        # Iterate and calculate distances only for strictly consecutive barcodes
        for i in range(len(barcodes) - 1):
            if barcodes[i + 1] == barcodes[i] + 1:  # Ensure strict consecutiveness
                dx_all.append(x_coords[i + 1] - x_coords[i])
                dy_all.append(y_coords[i + 1] - y_coords[i])
                dz_all.append(z_coords[i + 1] - z_coords[i])

    # Compute mean and standard deviation
    mean_dx, std_dx = (np.mean(dx_all), np.std(dx_all)) if dx_all else (0, 0)
    mean_dy, std_dy = (np.mean(dy_all), np.std(dy_all)) if dy_all else (0, 0)
    mean_dz, std_dz = (np.mean(dz_all), np.std(dz_all)) if dz_all else (0, 0)

    # Create figure with three histograms
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Distances between consecutive neighboring barcodes", fontsize=25)

    data = [dx_all, dy_all, dz_all]
    labels = [r"$\Delta x$, um", r"$\Delta y$, um", r"$\Delta z$, um"]
    colors = ["blue", "green", "red"]
    means = [mean_dx, mean_dy, mean_dz]
    stds = [std_dx, std_dy, std_dz]

    for i, (ax, dist, label, mean_val, std_val, color) in enumerate(
        zip(axes, data, labels, means, stds, colors)
    ):
        ax.hist(dist, bins=30, alpha=0.7, color=color, edgecolor="black")
        ax.set_xlabel(label)
        ax.set_ylabel("Counts")
        ax.set_title(
            f"Mean: {mean_val:.3f}\nStd: {std_val:.3f}", fontsize=10
        )  # Smaller title

    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"$ Saved neighbor distances plot: {output_filename}")

    return mean_dx, mean_dy, mean_dz, std_dx, std_dy, std_dz


def compute_radius_of_gyration(coords):
    """Computes the radius of gyration (Rg) for a given trace."""

    center_of_mass = np.mean(coords, axis=0)
    return np.sqrt(np.mean(np.sum((coords - center_of_mass) ** 2, axis=1)))


def plot_radius_of_gyration(
    trace,
    output_filename="radius_of_gyration.png",
    min_rg_threshold=0.02,
    selected_barcodes=None,
):
    """Plot the distribution of radius of gyration (Rg) across traces.

    Three panels are generated:
        - Histogram with single-species (one Gaussian) fit
        - Histogram with two-species (two Gaussian mixture) fit
        - Kernel density estimation with experimental datapoints as rug plot

    Parameters
    ----------
    trace : ChromatinTraceTable
        Trace table containing barcode coordinates.
    output_filename : str
        Path to save the output plot.
    min_rg_threshold : float, optional
        Minimum Rg value to keep when plotting. Defaults to 0.02.
    selected_barcodes : list[int] | None, optional
        If provided, only positions from these barcodes are used to compute Rg.
    """

    trace_table = trace.data
    if selected_barcodes is not None:
        trace_table = trace_table[np.isin(trace_table["Barcode #"], selected_barcodes)]

    trace_by_ID = trace_table.group_by("Trace_ID")

    rg_values = []
    for sub_trace_table in trace_by_ID.groups:
        coords = np.vstack((sub_trace_table["x"], sub_trace_table["y"], sub_trace_table["z"])).T
        if len(coords) > 0:
            rg_values.append(compute_radius_of_gyration(coords))

    if len(rg_values) == 0:
        print("! No traces available to compute radius of gyration.")
        return

    rg_values = np.array(rg_values)
    rg_values = rg_values[rg_values >= min_rg_threshold]

    if len(rg_values) == 0:
        print("! No Rg values passed the minimum threshold.")
        return
    min_rg, max_rg = np.min(rg_values), np.max(rg_values)
    if min_rg == max_rg:
        min_rg -= 0.5
        max_rg += 0.5
    x_range = np.linspace(min_rg, max_rg, 500)

    fig, axes = plt.subplots(1, 3, figsize=(24, 8), constrained_layout=True)
    fig.suptitle("Radius of gyration across traces", fontsize=20)

    # Panel 1: Single-species (one Gaussian) fit
    mu, sigma = np.mean(rg_values), np.std(rg_values)
    sigma = sigma if sigma > 0 else 1e-6
    axes[0].hist(
        rg_values, bins=50, density=True, alpha=0.5, color="tab:blue", label="Rg samples"
    )
    axes[0].plot(x_range, norm.pdf(x_range, mu, sigma), color="black", lw=2, label="1-species fit")
    axes[0].set_xlabel("Rg (um)", fontsize=20)
    axes[0].set_ylabel("Density", fontsize=20)
    axes[0].legend(fontsize=10)
    axes[0].set_title(f"$n$ = {len(rg_values)} | $\mu$ = {mu:.3f}, $\sigma$ = {sigma:.3f}")

    # Panel 2: Two-species Gaussian mixture fit
    axes[1].hist(
        rg_values, bins=50, density=True, alpha=0.5, color="tab:green", label="Rg samples"
    )
    if len(rg_values) > 1:
        try:
            gmm = GaussianMixture(n_components=2, random_state=42)
            gmm.fit(rg_values.reshape(-1, 1))
            gmm_density = np.exp(gmm.score_samples(x_range.reshape(-1, 1)))
            axes[1].plot(x_range, gmm_density, color="black", lw=2, label="2-species fit")

            means = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_).flatten()
            weights = gmm.weights_
            for mean_val, std_val, weight in zip(means, stds, weights):
                component_pdf = weight * norm.pdf(x_range, mean_val, std_val)
                axes[1].plot(x_range, component_pdf, lw=1.5, linestyle="--")

            axes[1].set_title(
                " | ".join(
                    [
                        f"$\mu_{i+1}$={mean_val:.3f}, $\sigma_{i+1}$={std_val:.3f}, $w_{i+1}$={weight:.2f}"
                        for i, (mean_val, std_val, weight) in enumerate(
                            zip(means, stds, weights)
                        )
                    ]
                )
            )
        except ValueError as exc:  # pragma: no cover - defensive programming
            print(f"! Could not fit two-species model: {exc}")
            axes[1].set_title("Two-species fit unavailable")
    else:
        axes[1].set_title("Two-species fit unavailable (insufficient data)")

    axes[1].set_xlabel("Rg (um)", fontsize=20)
    axes[1].set_ylabel("Density", fontsize=20)
    axes[1].legend(fontsize=10)

    # Panel 3: Kernel density estimation with rug plot
    bandwidth = 1.06 * np.std(rg_values) * (len(rg_values) ** (-1 / 5)) if len(rg_values) > 1 else 0.1
    bandwidth = bandwidth if bandwidth > 0 else 0.1
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
    kde.fit(rg_values.reshape(-1, 1))
    kde_density = np.exp(kde.score_samples(x_range.reshape(-1, 1)))

    axes[2].plot(x_range, kde_density, color="black", lw=2, label="KDE")
    axes[2].scatter(rg_values, np.zeros_like(rg_values), color="tab:red", marker="|", s=200, alpha=0.7)
    axes[2].set_xlabel("Rg (um)", fontsize=20)
    axes[2].set_ylabel("Density", fontsize=20)
    axes[2].legend(fontsize=10)
    axes[2].set_title(f"Bandwidth = {bandwidth:.3f}")

    plt.savefig(output_filename)
    print(f"$ Saved radius of gyration plot: {output_filename}")


def barcode_detection_efficiency(
    trace, output_prefix="barcode_detection_efficiency", format="png"
):
    """
    Analyze and visualize barcode detection efficiency across all traces.

    This function computes the frequency at which each barcode is detected across
    all traces and performs bootstrap analysis to estimate confidence intervals.
    The results are visualized as violin plots.

    Parameters
    ----------
    trace : ChromatinTraceTable
        Trace table, instance of the ChromatinTraceTable Class.
    output_prefix : str
        Prefix for the output filename (without extension).

    Returns
    -------
    None
        The function saves the output figure but does not return any values.
    """
    trace_table = trace.data
    trace_groups = trace_table.group_by("Trace_ID").groups
    n_traces = len(trace_groups)
    bootstrap_iterations = 1000

    print(f"$ Calculating overall barcode detection across {n_traces} traces...")

    barcode_presence = defaultdict(list)
    all_barcodes = set(trace_table["Barcode #"])

    # Track barcode presence per trace
    for sub_trace in trace_groups:
        present = set(sub_trace["Barcode #"])
        for barcode in all_barcodes:
            barcode_presence[barcode].append(1 if barcode in present else 0)

    # Bootstrap detection frequencies (generate bootstrapped mean values)
    detection_bootstrap_distributions = {}
    for barcode, detections in barcode_presence.items():
        detections = np.array(detections)
        boot_means = [
            np.mean(np.random.choice(detections, size=n_traces, replace=True))
            for _ in range(bootstrap_iterations)
        ]
        detection_bootstrap_distributions[barcode] = boot_means

    # Plotting
    sorted_barcodes = sorted(detection_bootstrap_distributions.keys())
    num_barcodes = len(sorted_barcodes)
    max_per_row = 50
    n_rows = (num_barcodes - 1) // max_per_row + 1

    fig_height = 6 * n_rows
    fig = plt.figure(figsize=(24, fig_height))
    gs = GridSpec(n_rows, 1, figure=fig)

    for row in range(n_rows):
        ax = fig.add_subplot(gs[row, 0])
        start = row * max_per_row
        end = min(start + max_per_row, num_barcodes)
        barcodes_row = sorted_barcodes[start:end]
        data = [detection_bootstrap_distributions[bc] for bc in barcodes_row]

        positions = np.arange(1, len(barcodes_row) + 1)
        ax.violinplot(data, showmedians=True, positions=positions)
        ax.set_xticks(positions)
        ax.set_xticklabels(barcodes_row)
        ax.set_xlim(0.5, len(barcodes_row) + 0.5)
        ax.set_ylabel("Detection frequency", fontsize=30)
        ax.set_xlabel("barcode IDs", fontsize=30)
        ax.set_ylim(0, 1)
        ax.set_title(
            f"Number of traces = {n_traces} | barcodes {start + 1}-{end}", fontsize=20
        )

    fig.tight_layout()
    fig.savefig(f"{output_prefix}.{format}")
    print(f"$ Exporting barcode detection plot to: {output_prefix}.{format}")


def analyze_trace(
    trace,
    trace_file,
    plotXYZ=False,
    format="png",
    min_rg_threshold=0.02,
    rg_barcodes=None,
):
    """
    Perform comprehensive analysis on a chromatin trace file.

    This function serves as a launcher for various analysis functions:
    - Barcode statistics analysis
    - Barcode detection efficiency analysis
    - Neighbor distance analysis
    - Barcode frequency analysis

    Parameters
    ----------
    trace : ChromatinTraceTable
        Trace table, instance of the ChromatinTraceTable Class.
    trace_file : str
        File name of trace table in ecsv format.
    plotXYZ : bool, optional
        Flag to control whether XYZ traces should be plotted. Default is False.
    format : str, optional
        Output file format for figures ('png' or 'svg'). Default is 'png'.
    min_rg_threshold : float, optional
        Minimum Rg value to include when computing the radius of gyration. Default is 0.02.
    rg_barcodes : list[int] | None, optional
        Barcodes to include when computing the radius of gyration. If None, all barcodes are used.

    Returns
    -------
    None
    """
    trace_table = trace.data

    print(f"$ Number of spots in trace file: {len(trace_table)}")

    # Get base filename without extension
    base_filename = trace_file.split(".")[0]

    # Calculate trace statistics
    output_filename = f"{base_filename}_trace_statistics.{format}"
    get_barcode_statistics(trace_table, output_filename)

    # Plot barcode detection per ROI with bootstrapped errors
    barcode_detection_efficiency(
        trace, output_prefix=base_filename + "_barcode_detection", format=format
    )

    # Compute and plot neighbor distances
    neighbor_distances_output = f"{base_filename}_first_neighbor_distances.{format}"
    mean_dx, mean_dy, mean_dz, std_dx, std_dy, std_dz = plot_neighbor_distances(
        trace, neighbor_distances_output
    )
    print(
        f"$ Mean distances between neighboring barcodes: X={mean_dx:.3f}, Y={mean_dy:.3f}, Z={mean_dz:.3f}"
    )

    # Plot radius of gyration distribution
    rg_output = f"{base_filename}_radius_of_gyration.{format}"
    plot_radius_of_gyration(
        trace,
        output_filename=rg_output,
        min_rg_threshold=min_rg_threshold,
        selected_barcodes=rg_barcodes,
    )

    # Plots how often barcodes are repeated in a single trace
    collective_barcode_stats = trace.barcode_statistics(trace_table)
    trace.plots_barcode_statistics(
        collective_barcode_stats,
        file_name=f"{base_filename}_relative_barcode_frequencies",
        kind="matrix",
        format=format,
    )


def process_traces(p):
    """
    Process a list of trace files and analyze each individually.

    This function iterates through the list of trace files, loads each one,
    and performs analyses on each trace file.

    Parameters
    ----------
    p : dict
        Dictionary containing processing parameters:
        - trace_files: List of trace files to process
        - plotXYZ: Flag to control whether XYZ traces should be plotted
        - format: Output image format (png or svg)

    Returns
    -------
    None
    """
    trace_files = p["trace_files"]

    if len(trace_files) > 0:
        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

        # iterates over traces in folder
        for trace_file in trace_files:
            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            if p["plotXYZ"]:
                print(f"> Plotting traces for {trace_file}")
                trace.plots_traces(
                    [trace_file.split(".")[0], "_traces_XYZ", f".{p['format']}"],
                    pixel_size=[0.1, 0.1, 0.25],
                )

            print(f"> Analyzing traces for {trace_file}")
            analyze_trace(
                trace,
                trace_file,
                plotXYZ=p["plotXYZ"],
                format=p["format"],
                min_rg_threshold=p["rgMin"],
                rg_barcodes=p["rgBarcodes"],
            )

    else:
        print(
            "! Error: did not find any trace file to analyze. Please provide one using --input or --pipe."
        )


def main():
    """
    Main function to execute the trace analyzer script.

    This function parses command-line arguments and initiates the trace analysis process.

    Returns
    -------
    None
    """
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)

    # [loops over lists of datafolders]
    process_traces(p)

    print("Finished execution")


if __name__ == "__main__":
    main()
