#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze chromatin trace files.
"""

from traceratops.script_banner import print_script_banner
import argparse
import select
import sys
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde

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
        "--output_format",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output image format. Default = png.",
    )
    parser.add_argument(
        "--neighbor_distance_range",
        default="1",
        help=(
            "Hexbin axis range for previous-vs-next neighbor distance plots. "
            "Use 'auto' for automatic scaling or provide a positive maximum "
            "distance in µm. Default = 1, which plots 0 to 1 µm."
        ),
    )
    return parser


def create_dict_args(args):
    p = {}
    p["input"] = args.input
    p["rootFolder"] = args.rootFolder
    p["plotXYZ"] = args.plotXYZ
    p["format"] = args.output_format
    p["neighbor_distance_range"] = args.neighbor_distance_range

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
    # trace_repeated_barcodes = list()
    number_unique_barcodes = list()
    # number_repeated_barcodes = list()

    for sub_trace_table in trace_by_ID.groups:
        trace_lengths.append(len(sub_trace_table))

        unique_barcodes = list(set(sub_trace_table["Barcode #"]))
        trace_unique_barcodes.append(unique_barcodes)
        number_unique_barcodes.append(len(unique_barcodes))

    distributions = [trace_lengths, number_unique_barcodes]
    axis_x_labels = [
        "$N_{barcodes}$",
        "$N_{unique-barcodes}$",
    ]

    number_plots = len(distributions)

    fig = plt.figure(constrained_layout=True)
    im_size = 12
    fig.set_size_inches((im_size * number_plots, im_size))
    gs = fig.add_gridspec(1, number_plots)
    axes = [fig.add_subplot(gs[0, i]) for i in range(number_plots)]
    bins = np.arange(1, np.max(number_unique_barcodes))

    for axis, distribution, xlabel in zip(axes, distributions, axis_x_labels):
        axis.hist(distribution, bins=bins, alpha=0.3)
        axis.set_xlabel(xlabel, fontsize=30)
        axis.set_ylabel("counts", fontsize=30)
        axis.set_title(
            f"n = {str(len(distribution))} | median = {str(np.median(distribution))}",
            fontsize=20,
        )

    fig.suptitle("Trace statistics", fontsize=40)

    plt.savefig(output_filename)


def _get_genomic_barcode_order(trace_table):
    """Return barcode order after sorting unique barcodes by genomic coordinates."""
    genomic_columns = ["Chrom", "Chrom_Start", "Chrom_End"]
    missing_columns = [
        col for col in genomic_columns if col not in trace_table.colnames
    ]
    if missing_columns:
        raise ValueError(
            "Trace table is missing genomic coordinate columns required for "
            f"neighbor ordering: {', '.join(missing_columns)}"
        )

    barcode_coordinates = {}
    for row in trace_table:
        barcode = row["Barcode #"]
        genomic_coordinate = (row["Chrom"], row["Chrom_Start"], row["Chrom_End"])
        if (
            barcode in barcode_coordinates
            and barcode_coordinates[barcode] != genomic_coordinate
        ):
            raise ValueError(
                f"Barcode {barcode} has multiple genomic coordinates: "
                f"{barcode_coordinates[barcode]} and {genomic_coordinate}"
            )
        barcode_coordinates[barcode] = genomic_coordinate

    return [
        barcode
        for barcode, _ in sorted(
            barcode_coordinates.items(),
            key=lambda item: (item[1][0], item[1][1], item[1][2]),
        )
    ]


def _resolve_neighbor_distance_limits(distance_range, distance_values):
    """Return shared hexbin axis limits for neighbor-distance density plots."""
    if isinstance(distance_range, str) and distance_range.lower() == "auto":
        if len(distance_values) == 0:
            return (0, 1)
        max_distance = float(np.nanmax(distance_values))
        return (0, max_distance if max_distance > 0 else 1)

    try:
        max_distance = float(distance_range)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "neighbor_distance_range must be 'auto' or a positive number in µm"
        ) from exc

    if max_distance <= 0:
        raise ValueError("neighbor_distance_range must be greater than 0 µm")

    return (0, max_distance)


def _plot_neighbor_hexbin(ax, previous_distances, next_distances, title, axis_limits):
    """Plot density of previous-versus-next neighbor distances as a hexbin map."""
    ax.set_xlabel(r"$d(i+1,i)$, µm", fontsize=18, labelpad=4)
    ax.set_ylabel(r"$d(i,i-1)$, µm", fontsize=18, labelpad=4)
    ax.set_title(title, fontsize=15, pad=6)
    ax.set_xlim(axis_limits)
    ax.set_ylim(axis_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", labelsize=16)

    if len(previous_distances) == 0 or len(next_distances) == 0:
        ax.text(0.5, 0.5, "No consecutive triplets", ha="center", va="center")
        return None

    return ax.hexbin(
        next_distances,
        previous_distances,
        gridsize=45,
        mincnt=1,
        cmap="cubehelix_r",
        extent=(*axis_limits, *axis_limits),
    )


def plot_neighbor_distances(
    trace, output_filename="neighbor_distances.png", neighbor_distance_range="1"
):
    """
    Calculate and visualize distances between consecutive neighboring barcodes.

    This function computes X, Y, and Z distances between barcodes that are
    consecutive both numerically and in genomic-coordinate order. It also plots
    hexbin density maps comparing each barcode's distance to its next and
    previous genomic neighbors.

    Parameters
    ----------
    trace : ChromatinTraceTable
        Trace table, instance of the ChromatinTraceTable Class.
    output_filename : str
        The filename for the output figure file.
    neighbor_distance_range : str or float
        Hexbin axis range. Use "auto" to scale from the data, or provide a
        positive maximum distance in µm for fixed axes from 0 to that value.

    Returns
    -------
    tuple
        (mean_dx, mean_dy, mean_dz, std_dx, std_dy, std_dz):
        Mean and standard deviation values for X, Y, and Z distances.
    """
    trace_table = trace.data
    trace_by_ID = trace_table.group_by("Trace_ID")
    genomic_barcode_order = _get_genomic_barcode_order(trace_table)
    genomic_barcode_rank = {
        barcode: rank for rank, barcode in enumerate(genomic_barcode_order)
    }

    dx_all, dy_all, dz_all = [], [], []
    xy_previous_all, xy_next_all = [], []
    xyz_previous_all, xyz_next_all = [], []

    for sub_trace_table in trace_by_ID.groups:
        # Sort barcodes according to ascending genomic coordinates.
        sorted_trace = sub_trace_table[
            np.argsort(
                [
                    genomic_barcode_rank[barcode]
                    for barcode in sub_trace_table["Barcode #"]
                ]
            )
        ]

        barcodes = sorted_trace["Barcode #"].data
        x_coords = sorted_trace["x"].data
        y_coords = sorted_trace["y"].data
        z_coords = sorted_trace["z"].data
        neighbor_distances = {}

        for i in range(len(barcodes) - 1):
            is_numeric_neighbor = barcodes[i + 1] == barcodes[i] + 1
            is_genomic_neighbor = (
                genomic_barcode_rank[barcodes[i + 1]]
                == genomic_barcode_rank[barcodes[i]] + 1
            )
            if is_numeric_neighbor and is_genomic_neighbor:
                dx = x_coords[i + 1] - x_coords[i]
                dy = y_coords[i + 1] - y_coords[i]
                dz = z_coords[i + 1] - z_coords[i]
                dx_all.append(dx)
                dy_all.append(dy)
                dz_all.append(dz)
                neighbor_distances[barcodes[i]] = (
                    np.hypot(dx, dy),
                    np.sqrt(dx**2 + dy**2 + dz**2),
                )

        for i in range(1, len(barcodes) - 1):
            previous_barcode = barcodes[i - 1]
            current_barcode = barcodes[i]
            if (
                previous_barcode in neighbor_distances
                and current_barcode in neighbor_distances
            ):
                previous_xy, previous_xyz = neighbor_distances[previous_barcode]
                next_xy, next_xyz = neighbor_distances[current_barcode]
                xy_previous_all.append(abs(previous_xy))
                xy_next_all.append(abs(next_xy))
                xyz_previous_all.append(abs(previous_xyz))
                xyz_next_all.append(abs(next_xyz))

    # Compute mean and standard deviation
    mean_dx, std_dx = (np.mean(dx_all), np.std(dx_all)) if dx_all else (0, 0)
    mean_dy, std_dy = (np.mean(dy_all), np.std(dy_all)) if dy_all else (0, 0)
    mean_dz, std_dz = (np.mean(dz_all), np.std(dz_all)) if dz_all else (0, 0)

    all_neighbor_distances = np.array(
        xy_previous_all + xy_next_all + xyz_previous_all + xyz_next_all
    )
    axis_limits = _resolve_neighbor_distance_limits(
        neighbor_distance_range, all_neighbor_distances
    )

    fig = plt.figure(figsize=(18, 13))
    outer = fig.add_gridspec(
        2, 1,
        height_ratios=[2.2, 1],
        hspace=0.55,
        left=0.07,
        right=0.96,
        top=0.88,
        bottom=0.08,
    )

    top_gs = outer[0].subgridspec(
        1, 4,
        width_ratios=[1, 0.035, 1, 0.035],
        wspace=0.35,
    )

    bottom_gs = outer[1].subgridspec(
        1, 3,
        wspace=0.75,
    )

    hexbin_axes = [
        fig.add_subplot(top_gs[0, 0]),
        fig.add_subplot(top_gs[0, 2]),
    ]

    colorbar_axes = [
        fig.add_subplot(top_gs[0, 1]),
        fig.add_subplot(top_gs[0, 3]),
    ]

    hist_axes = [
        fig.add_subplot(bottom_gs[0, 0]),
        fig.add_subplot(bottom_gs[0, 1]),
        fig.add_subplot(bottom_gs[0, 2]),
    ]

    fig.suptitle(
        "Distances between consecutive genomic neighboring barcodes",
        fontsize=22,
        y=0.97,
    )

    data = [dx_all, dy_all, dz_all]
    labels = [r"$\Delta x$, µm", r"$\Delta y$, µm", r"$\Delta z$, µm"]
    colors = ["blue", "green", "red"]
    means = [mean_dx, mean_dy, mean_dz]
    stds = [std_dx, std_dy, std_dz]

    for ax, dist, label, mean_val, std_val, color in zip(
        hist_axes, data, labels, means, stds, colors
    ):
        ax.hist(dist, bins=30, alpha=0.7, color=color, edgecolor="black")
        ax.set_xlabel(label, fontsize=18)
        ax.set_ylabel("Counts", fontsize=18)
        ax.set_title(f"Mean: {mean_val:.3f}\nStd: {std_val:.3f}", fontsize=14, pad=6)
        ax.tick_params(axis="both", labelsize=16)

    xy_hexbin = _plot_neighbor_hexbin(
        hexbin_axes[0],
        xy_previous_all,
        xy_next_all,
        "XY neighbor-distance density",
        axis_limits,
    )
    xyz_hexbin = _plot_neighbor_hexbin(
        hexbin_axes[1],
        xyz_previous_all,
        xyz_next_all,
        "XYZ neighbor-distance density",
        axis_limits,
    )
    hexbin_axes[1].set_ylabel("")

    if xy_hexbin is not None:
        colorbar = fig.colorbar(xy_hexbin, 
                                cax=colorbar_axes[0],
                                )                   
        colorbar.ax.set_title("Counts", fontsize=12, pad=6)
        colorbar.ax.tick_params(labelsize=12)
    else:
        colorbar_axes[0].axis("off")
    if xyz_hexbin is not None:
        colorbar = fig.colorbar(xyz_hexbin, 
                                cax=colorbar_axes[1],
                                )                   
        colorbar.ax.set_title("Counts", fontsize=12, pad=6)
        colorbar.ax.tick_params(labelsize=12)
    else:
        colorbar_axes[1].axis("off")

    plt.savefig(output_filename)
    plt.close(fig)
    print(f"$ Saved neighbor distances plot: {output_filename}")

    return mean_dx, mean_dy, mean_dz, std_dx, std_dy, std_dz


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


def compute_kde(x, y, grid_size=200):
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)

    xi, yi = np.meshgrid(
        np.linspace(x.min(), x.max(), grid_size),
        np.linspace(y.min(), y.max(), grid_size),
    )

    zi = kde(np.vstack([xi.ravel(), yi.ravel()]))
    return zi.reshape(xi.shape)


def plot_kde_projections(trace_table, output_filename, target_ratio=0.5):
    """
    Plot KDE projections (XY, XZ, YZ) from trace data.

    Parameters
    ----------
    trace_table : astropy.table.Table
    output_filename : str
    target_ratio : float
        Controls mild Z stretching
    """
    with matplotlib.rc_context({"font.size": 10}):
        x = np.array(trace_table["x"])
        y = np.array(trace_table["y"])
        z = np.array(trace_table["z"])

        # === Z scaling (mild) ===
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        z_range = z.max() - z.min()

        xy_range = 0.5 * (x_range + y_range)

        z_scale = (xy_range * target_ratio) / z_range
        z_scaled = z * z_scale

        def _plot(ax, x, y, xlabel, ylabel, title, show_z_ticks=False):
            zi = compute_kde(x, y)

            im = ax.imshow(
                zi,
                origin="lower",
                extent=[x.min(), x.max(), y.min(), y.max()],
            )

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.set_aspect("equal")
            ax.grid(alpha=0.3)

            if show_z_ticks:
                z_ticks = np.arange(np.floor(z.min()), np.ceil(z.max()) + 1, 1)
                ax.set_yticks(z_ticks * z_scale)
                ax.set_yticklabels([f"{int(t)}" for t in z_ticks])

            return im

        # === FIGURE ===
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(2, 3, width_ratios=[1, 1, 0.05], hspace=0.05, wspace=0.15)

        ax_xy = fig.add_subplot(gs[:, 0])
        ax_xz = fig.add_subplot(gs[0, 1])
        ax_yz = fig.add_subplot(gs[1, 1])
        cax = fig.add_subplot(gs[:, 2])

        im = _plot(ax_xy, x, y, "X (µm)", "Y (µm)", "XY projection (KDE)")

        _plot(
            ax_xz, x, z_scaled, "X (µm)", "Z (µm)", "XZ projection", show_z_ticks=True
        )

        _plot(
            ax_yz, y, z_scaled, "Y (µm)", "Z (µm)", "YZ projection", show_z_ticks=True
        )

        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("Probability density")

        fig.suptitle("Spatial KDE projections", fontsize=30)

        fig.subplots_adjust(
            left=0.07, right=0.92, top=0.88, bottom=0.08, wspace=0.25, hspace=0.15
        )
        plt.savefig(output_filename)
        plt.close(fig)

        print(f"$ Saved KDE projection plot: {output_filename}")


def analyze_trace(
    trace, trace_file, plotXYZ=False, format="png", neighbor_distance_range="1"
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
        Output file format for figures ('png', 'svg', or 'pdf'). Default is 'png'.
    neighbor_distance_range : str or float, optional
        Hexbin axis range. Use "auto" to scale from the data, or provide a
        positive maximum distance in µm for fixed axes from 0 to that value.

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
        trace,
        neighbor_distances_output,
        neighbor_distance_range=neighbor_distance_range,
    )
    print(
        f"$ Mean distances between neighboring barcodes: X={mean_dx:.3f}, Y={mean_dy:.3f}, Z={mean_dz:.3f}"
    )

    # Plots how often barcodes are repeated in a single trace
    collective_barcode_stats = trace.barcode_statistics(trace_table)
    trace.plots_barcode_statistics(
        collective_barcode_stats,
        file_name=f"{base_filename}_relative_barcode_frequencies",
        kind="matrix",
        format=format,
    )

    # === KDE spatial projections ===
    kde_output = f"{base_filename}_kde_projections.{format}"
    plot_kde_projections(trace_table, kde_output)


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
        - format: Output image format (png, svg, or pdf)

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
                neighbor_distance_range=p["neighbor_distance_range"],
            )

    else:
        print(
            "! Error: did not find any trace file to analyze. Please provide one using --input or --pipe."
        )


def main():
    print_script_banner(__file__, __doc__)
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
