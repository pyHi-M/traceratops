#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot median physical distance as a function of genomic distance.

The script reads one or more chromatin trace tables and, for each trace,
collects all pairwise 3D distances between spots with different barcodes. It
aggregates those measurements by barcode pair, calculates the median physical
separation and bootstrap error bars for each pair, then plots those medians
against the genomic separation between the barcode loci. A power-law fit is
added to the plot with a bootstrap confidence interval.
"""

from traceratops.script_banner import print_script_banner
import argparse
import os
import select
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", help="Path to input trace table (ECSV, CSV, or 4DN)."
    )
    parser.add_argument(
        "-O",
        "--output",
        default="physical_vs_genomic_distance.png",
        help="Output plot filename. Default = physical_vs_genomic_distance.png.",
    )
    parser.add_argument(
        "--output_data",
        help="Optional CSV filename for the aggregated barcode-pair statistics.",
    )
    parser.add_argument(
        "--output_format",
        choices=["png", "svg", "pdf"],
        default=None,
        help="Output image format. Default = inferred from --output, or png.",
    )
    parser.add_argument(
        "--bootstrap_cycles",
        type=int,
        default=1000,
        help="Number of bootstrap cycles used for median error bars. Default = 1000.",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.95,
        help="Confidence level for error bars and fit band. Default = 0.95.",
    )
    parser.add_argument(
        "--fit_bootstrap_cycles",
        type=int,
        default=1000,
        help="Number of parameter samples used for the fit confidence band. Default = 1000.",
    )
    parser.add_argument(
        "--distance_unit",
        default="µm",
        help="Physical distance unit label. Default = µm.",
    )
    parser.add_argument(
        "--genomic_unit",
        choices=["bp", "kb", "Mb"],
        default="kb",
        help="Genomic distance unit for plotting. Default = kb.",
    )
    parser.add_argument(
        "--min_observations",
        type=int,
        default=1,
        help="Minimum observations required to plot a barcode pair. Default = 1.",
    )
    parser.add_argument("--x_min", type=float, default=None, help="Minimum x-axis value.")
    parser.add_argument("--x_max", type=float, default=None, help="Maximum x-axis value.")
    parser.add_argument("--y_min", type=float, default=None, help="Minimum y-axis value.")
    parser.add_argument("--y_max", type=float, default=None, help="Maximum y-axis value.")
    parser.add_argument(
        "--pipe", help="Read trace file paths from stdin (pipe).", action="store_true"
    )
    return parser


def create_dict_args(args):
    trace_files = []
    if args.pipe:
        if select.select([sys.stdin], [], [], 0.0)[0]:
            trace_files = [line.rstrip("\n") for line in sys.stdin if line.strip()]
        else:
            print("Nothing in stdin")
    elif args.input:
        trace_files = [args.input]

    output_format = args.output_format
    if output_format is None:
        ext = os.path.splitext(args.output)[1].lower().lstrip(".")
        output_format = ext if ext in ["png", "svg", "pdf"] else "png"

    return {
        "trace_files": trace_files,
        "output": args.output,
        "output_data": args.output_data,
        "output_format": output_format,
        "bootstrap_cycles": args.bootstrap_cycles,
        "confidence": args.confidence,
        "fit_bootstrap_cycles": args.fit_bootstrap_cycles,
        "distance_unit": args.distance_unit,
        "genomic_unit": args.genomic_unit,
        "min_observations": args.min_observations,
        "x_min": args.x_min,
        "x_max": args.x_max,
        "y_min": args.y_min,
        "y_max": args.y_max,
    }


def genomic_scale(unit):
    return {"bp": 1.0, "kb": 1_000.0, "Mb": 1_000_000.0}[unit]


def locus_midpoint(row):
    return (float(row["Chrom_Start"]) + float(row["Chrom_End"])) / 2.0


def collect_distances(trace_table):
    """Collect physical distances by barcode pair across all traces."""
    distances = {}
    genomic_distances = {}
    for trace in trace_table.group_by("Trace_ID").groups:
        for i in range(len(trace) - 1):
            for j in range(i + 1, len(trace)):
                barcode_i = int(trace[i]["Barcode #"])
                barcode_j = int(trace[j]["Barcode #"])
                if barcode_i == barcode_j:
                    continue
                pair = tuple(sorted((barcode_i, barcode_j)))
                physical_distance = np.linalg.norm(
                    np.array([trace[i]["x"], trace[i]["y"], trace[i]["z"]], dtype=float)
                    - np.array([trace[j]["x"], trace[j]["y"], trace[j]["z"]], dtype=float)
                )
                genomic_distance = abs(
                    locus_midpoint(trace[i]) - locus_midpoint(trace[j])
                )
                distances.setdefault(pair, []).append(float(physical_distance))
                genomic_distances.setdefault(pair, float(genomic_distance))
    return distances, genomic_distances


def bootstrap_median_interval(values, cycles=1000, confidence=0.95, rng=None):
    values = np.asarray(values, dtype=float)
    if len(values) == 1 or cycles <= 0:
        median = float(np.median(values))
        return median, 0.0, 0.0
    if rng is None:
        rng = np.random.default_rng(42)
    samples = rng.choice(values, size=(cycles, len(values)), replace=True)
    medians = np.median(samples, axis=1)
    alpha = (1.0 - confidence) / 2.0
    low, high = np.quantile(medians, [alpha, 1.0 - alpha])
    median = float(np.median(values))
    return median, max(0.0, median - low), max(0.0, high - median)


def summarize_pairs(
    distances, genomic_distances, bootstrap_cycles, confidence, min_observations
):
    rows = []
    rng = np.random.default_rng(42)
    for pair, values in distances.items():
        if len(values) < min_observations:
            continue
        median, err_low, err_high = bootstrap_median_interval(
            values, bootstrap_cycles, confidence, rng
        )
        rows.append(
            {
                "barcode_1": pair[0],
                "barcode_2": pair[1],
                "genomic_distance_bp": genomic_distances[pair],
                "median_physical_distance": median,
                "error_low": err_low,
                "error_high": err_high,
                "n_observations": len(values),
            }
        )
    rows.sort(key=lambda row: row["genomic_distance_bp"])
    return rows


def power_law(x, coefficient, exponent):
    return coefficient * np.power(x, exponent)


def fit_power_law(x, y, yerr, confidence=0.95, cycles=1000):
    positive = (x > 0) & (y > 0)
    x_fit = x[positive]
    y_fit = y[positive]
    yerr_fit = yerr[positive]
    fallback_sigma = (
        np.nanmedian(yerr_fit[yerr_fit > 0]) if np.any(yerr_fit > 0) else 1.0
    )
    sigma = np.where(yerr_fit > 0, yerr_fit, fallback_sigma)
    params, covariance = curve_fit(
        power_law,
        x_fit,
        y_fit,
        p0=(np.nanmedian(y_fit) / np.nanmedian(x_fit) ** 0.3, 0.3),
        sigma=sigma,
        absolute_sigma=False,
        maxfev=10000,
    )
    x_line = np.linspace(np.min(x_fit), np.max(x_fit), 300)
    y_line = power_law(x_line, *params)
    rng = np.random.default_rng(42)
    try:
        sampled_params = rng.multivariate_normal(params, covariance, size=cycles)
        sampled_curves = np.array(
            [power_law(x_line, *sample) for sample in sampled_params]
        )
        alpha = (1.0 - confidence) / 2.0
        low, high = np.quantile(sampled_curves, [alpha, 1.0 - alpha], axis=0)
    except (ValueError, np.linalg.LinAlgError):
        low = y_line
        high = y_line
    return params, x_line, y_line, low, high


def write_summary_csv(rows, path, genomic_unit):
    import csv

    scale = genomic_scale(genomic_unit)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "barcode_1",
                "barcode_2",
                f"genomic_distance_{genomic_unit}",
                "genomic_distance_bp",
                "median_physical_distance",
                "error_low",
                "error_high",
                "n_observations",
            ],
        )
        writer.writeheader()
        for row in rows:
            output_row = dict(row)
            output_row[f"genomic_distance_{genomic_unit}"] = (
                row["genomic_distance_bp"] / scale
            )
            writer.writerow(output_row)


def plot_summary(
    rows, output, output_format, genomic_unit, distance_unit, confidence, fit_cycles, limits
):
    scale = genomic_scale(genomic_unit)
    x = np.array([row["genomic_distance_bp"] / scale for row in rows], dtype=float)
    y = np.array([row["median_physical_distance"] for row in rows], dtype=float)
    yerr = np.array(
        [[row["error_low"] for row in rows], [row["error_high"] for row in rows]],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        markersize=4,
        markerfacecolor="white",
        markeredgecolor="blue",
        ecolor="cornflowerblue",
        elinewidth=0.8,
        capsize=2,
        linestyle="none",
        alpha=0.8,
        label="Barcode-pair medians",
    )

    if np.count_nonzero((x > 0) & (y > 0)) >= 2:
        symmetric_yerr = np.mean(yerr, axis=0)
        params, x_line, y_line, low, high = fit_power_law(
            x, y, symmetric_yerr, confidence, fit_cycles
        )
        ax.plot(
            x_line,
            y_line,
            "r--",
            linewidth=2,
            label=f"Fit: y = {params[0]:.3g} x^{params[1]:.3g}",
        )
        ax.fill_between(
            x_line,
            low,
            high,
            color="red",
            alpha=0.2,
            label=f"{confidence:.0%} fit CI",
        )

    ax.set_xlabel(f"Genomic distance, {genomic_unit}")
    ax.set_ylabel(f"Median distance, {distance_unit}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)
    if limits["x_min"] is not None or limits["x_max"] is not None:
        ax.set_xlim(limits["x_min"], limits["x_max"])
    if limits["y_min"] is not None or limits["y_max"] is not None:
        ax.set_ylim(limits["y_min"], limits["y_max"])
    fig.tight_layout()
    fig.savefig(output, format=output_format, dpi=300)
    plt.close(fig)


def runtime(**kwargs):
    trace_files = kwargs["trace_files"]
    if len(trace_files) < 1:
        print("! Error: no trace file provided. Use --input or --pipe.")
        return 0

    all_distances = {}
    all_genomic_distances = {}
    for trace_file in trace_files:
        trace_table = ChromatinTraceTable()
        table = trace_table.load(trace_file)
        distances, genomic_distances = collect_distances(table)
        for pair, values in distances.items():
            all_distances.setdefault(pair, []).extend(values)
            all_genomic_distances.setdefault(pair, genomic_distances[pair])

    rows = summarize_pairs(
        all_distances,
        all_genomic_distances,
        kwargs["bootstrap_cycles"],
        kwargs["confidence"],
        kwargs["min_observations"],
    )
    if len(rows) == 0:
        print("! Error: no barcode pairs with enough observations were found.")
        return 0

    if kwargs["output_data"]:
        write_summary_csv(rows, kwargs["output_data"], kwargs["genomic_unit"])

    plot_summary(
        rows,
        kwargs["output"],
        kwargs["output_format"],
        kwargs["genomic_unit"],
        kwargs["distance_unit"],
        kwargs["confidence"],
        kwargs["fit_bootstrap_cycles"],
        {
            "x_min": kwargs["x_min"],
            "x_max": kwargs["x_max"],
            "y_min": kwargs["y_min"],
            "y_max": kwargs["y_max"],
        },
    )
    return len(rows)


def main():
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    if len(sys.argv) == 1:
        print("Error: No argument provided")
        print("Redirecting to `--help` option:")
        parser.print_help()
        return
    args = parser.parse_args()
    p = create_dict_args(args)
    runtime(**p)


if __name__ == "__main__":
    main()
