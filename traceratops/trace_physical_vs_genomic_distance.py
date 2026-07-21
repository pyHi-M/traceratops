#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
`trace_physical_vs_genomic_distance.py` summarizes how physical chromatin-trace distances change with genomic separation. It reads a chromatin trace table, computes pairwise distances between all barcode loci in each trace, bins those distances by genomic distance, and writes a CSV table that can also be plotted as a log-log physical-vs-genomic-distance curve.

1. The script loads an input trace table supported by `ChromatinTraceTable` (`.ecsv`, `.dat`, `.4dn`, or `.csv`).
2. Each trace is converted into a NumPy array with one row per trace, one column per barcode, and X/Y/Z coordinates in the final dimension.
3. Genomic separation is computed from barcode## Description midpoint positions, using `(Chrom_Start + Chrom_End) / 2`, and is reported in kilobase pairs (kbp).
4. Physical pairwise distances are computed for:
   - full 3D Euclidean distance, unless `--no_3d` is used;
   - projected X-axis distance;
   - projected Y-axis distance;
   - projected Z-axis distance.
5. Physical distances above `--dist_threshold` are replaced by `NaN` before summarization.
6. Genomic distances are split into `--gen_dist_bins` equally spaced bins, and the median physical distance is reported for each bin and axis.


"""

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import unique
from scipy import stats
from scipy.spatial.distance import pdist
from tqdm import tqdm

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.script_banner import print_script_banner

AXES = {"X": 0, "Y": 1, "Z": 2}


def extract_traces_numpy(table: ChromatinTraceTable) -> np.ndarray:
    """Extract trace coordinates with shape ``(num_traces, num_barcodes, 3)``."""
    unique_traces = unique(table.data, keys="Trace_ID")
    unique_barcodes = unique(table.data, keys="Barcode #")["Barcode #"].data
    barcode_to_index = {barcode: index for index, barcode in enumerate(unique_barcodes)}

    traces = np.full((len(unique_traces), len(unique_barcodes), 3), np.nan)
    for trace_index, trace in enumerate(table.data.group_by("Trace_ID").groups):
        barcodes = np.asarray(trace["Barcode #"].data)
        coordinates = np.column_stack(
            (
                np.asarray(trace["x"].data),
                np.asarray(trace["y"].data),
                np.asarray(trace["z"].data),
            )
        )
        barcode_indices = np.fromiter(
            (barcode_to_index[barcode] for barcode in barcodes),
            dtype=np.intp,
            count=len(barcodes),
        )
        traces[trace_index, barcode_indices, :] = coordinates

    return traces


def compute_genomic_dist_map(table: ChromatinTraceTable) -> np.ndarray:
    """Compute barcode-pair genomic separations in kilobase pairs."""
    barcode_table = unique(table.data, keys="Barcode #")
    mean_pos = (barcode_table["Chrom_Start"].data + barcode_table["Chrom_End"].data) / 2
    return pdist(mean_pos.reshape(-1, 1), metric="euclidean") / 1000


def compute_inter_loci_genomic_dist(table: ChromatinTraceTable) -> pd.DataFrame:
    """Compute genomic separations between consecutive barcode loci in kilobase pairs."""
    barcode_table = unique(table.data, keys="Barcode #")
    mean_pos = (barcode_table["Chrom_Start"].data + barcode_table["Chrom_End"].data) / 2
    distances = np.diff(mean_pos) / 1000
    return pd.DataFrame(
        {
            "locus": [f"locus_{index}_{index + 1}" for index in range(len(distances))],
            "genomic distance (kbp)": distances.astype(int),
        }
    )


def compute_all_pwd(
    traces: np.ndarray, show_progress: bool = True
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Compute 3D and per-axis pairwise physical distances for every trace."""
    n_traces, n_barcodes, n_dimensions = traces.shape
    n_pairs = n_barcodes * (n_barcodes - 1) // 2
    pwd = np.empty((n_traces, n_pairs))
    single_axis_pwd = {axis: np.empty((n_traces, n_pairs)) for axis in AXES}

    iterator = tqdm(
        range(n_traces),
        desc="Calculating physical distances",
        disable=not show_progress,
    )
    for trace_index in iterator:
        pwd[trace_index] = pdist(traces[trace_index])
        for axis, dimension in AXES.items():
            if dimension < n_dimensions:
                single_axis_pwd[axis][trace_index] = pdist(
                    traces[trace_index, :, dimension].reshape(-1, 1)
                )

    return pwd, single_axis_pwd


def compute_physical_vs_genomic_distance(
    pairwise_dist: np.ndarray,
    gen_dist: np.ndarray,
    gen_dist_bins: int,
    dist_threshold: float = np.inf,
    experiment: Optional[str] = None,
    axis: str = "3D",
    show_progress: bool = True,
) -> pd.DataFrame:
    """Bin physical distances by genomic distance and return median distance summaries."""
    merged_dist = np.column_stack(
        (np.tile(gen_dist, pairwise_dist.shape[0]), pairwise_dist.reshape(-1))
    )
    merged_dist[merged_dist[:, 1] > dist_threshold, 1] = np.nan

    gen_bins = np.linspace(np.nanmin(gen_dist), np.nanmax(gen_dist), gen_dist_bins + 1)
    records = []
    iterator = tqdm(
        range(len(gen_bins) - 1),
        desc=f"Binning physical vs genomic distance for {axis}",
        disable=not show_progress,
    )
    for bin_index in iterator:
        lower = gen_bins[bin_index]
        upper = gen_bins[bin_index + 1]
        if bin_index == len(gen_bins) - 2:
            selected = (merged_dist[:, 0] >= lower) & (merged_dist[:, 0] <= upper)
        else:
            selected = (merged_dist[:, 0] >= lower) & (merged_dist[:, 0] < upper)
        values = merged_dist[selected, 1]
        median_distance = np.nanmedian(values) if values.size else np.nan
        midpoint = (lower + upper) / 2
        records.append(
            {
                "genomic distance (kbp)": midpoint,
                "log10 genomic dist (kbp)": (
                    np.log10(midpoint) if midpoint > 0 else np.nan
                ),
                "median euclidean distance (nm)": median_distance,
                "log10 median dist (nm)": (
                    np.log10(median_distance) if median_distance > 0 else np.nan
                ),
                "n_data": np.count_nonzero(~np.isnan(values)),
                "experiment": (
                    f"exp_{experiment}" if experiment is not None else "exp_None"
                ),
                "axis": axis,
            }
        )

    return pd.DataFrame.from_records(records)


def fit_power_law(
    df: pd.DataFrame, confidence: float = 0.95, n_points: int = 200
) -> Optional[Dict[str, np.ndarray]]:
    """Fit ``physical_distance = coefficient * genomic_distance ** exponent``.

    The fit is estimated by ordinary least squares in log10 space. Confidence
    intervals describe the fitted mean response transformed back to linear space.
    """
    fit_df = df.dropna(
        subset=["genomic distance (kbp)", "median euclidean distance (nm)"]
    )
    fit_df = fit_df.loc[
        (fit_df["genomic distance (kbp)"] > 0)
        & (fit_df["median euclidean distance (nm)"] > 0)
    ]
    if len(fit_df) < 3:
        return None

    x_log = np.log10(fit_df["genomic distance (kbp)"].to_numpy(dtype=float))
    y_log = np.log10(fit_df["median euclidean distance (nm)"].to_numpy(dtype=float))
    slope, intercept = np.polyfit(x_log, y_log, 1)

    x_fit_log = np.linspace(np.nanmin(x_log), np.nanmax(x_log), n_points)
    y_fit_log = intercept + slope * x_fit_log

    residuals = y_log - (intercept + slope * x_log)
    degrees_of_freedom = len(x_log) - 2
    residual_std_error = np.sqrt(np.sum(residuals**2) / degrees_of_freedom)
    x_mean = np.mean(x_log)
    sum_squared_x = np.sum((x_log - x_mean) ** 2)
    if sum_squared_x == 0:
        return None
    t_value = stats.t.ppf((1 + confidence) / 2, degrees_of_freedom)
    mean_se = residual_std_error * np.sqrt(
        (1 / len(x_log)) + ((x_fit_log - x_mean) ** 2 / sum_squared_x)
    )

    lower_log = y_fit_log - t_value * mean_se
    upper_log = y_fit_log + t_value * mean_se
    return {
        "coefficient": 10**intercept,
        "exponent": slope,
        "x_fit": 10**x_fit_log,
        "y_fit": 10**y_fit_log,
        "lower": 10**lower_log,
        "upper": 10**upper_log,
        "x_fit_log": x_fit_log,
        "y_fit_log": y_fit_log,
        "lower_log": lower_log,
        "upper_log": upper_log,
    }


def plot_log_distance_graph(
    dist_df: pd.DataFrame, saving_filename: Union[str, Path]
) -> None:
    """Plot log10 median physical distance versus log10 genomic distance."""
    import seaborn as sns

    axes = [axis for axis in ["3D", "X", "Y", "Z"] if axis in set(dist_df["axis"])]
    sns.set_context("paper")
    fig, axs = plt.subplots(
        nrows=len(axes),
        ncols=1,
        squeeze=False,
        figsize=(6.5, 2.2 * len(axes)),
        sharex=True,
        constrained_layout=True,
    )
    palette = dict(
        zip(
            dist_df["experiment"].dropna().unique(),
            sns.color_palette(n_colors=dist_df["experiment"].nunique()),
        )
    )
    for row_index, axis_name in enumerate(axes):
        ax = axs[row_index, 0]
        df = dist_df.loc[dist_df["axis"] == axis_name]
        sns.scatterplot(
            data=df,
            x="log10 genomic dist (kbp)",
            y="log10 median dist (nm)",
            hue="experiment",
            palette=palette,
            s=16,
            linewidth=0,
            legend=False,
            ax=ax,
        )
        for experiment, exp_df in df.groupby("experiment", dropna=False):
            fit = fit_power_law(exp_df)
            if fit is None:
                continue
            label = (
                f"{experiment}: a={fit['coefficient']:.2g}, " f"b={fit['exponent']:.3f}"
            )
            ax.plot(
                fit["x_fit_log"],
                fit["y_fit_log"],
                color="r",
                lw=1.2,
                label=label,
            )
            ax.fill_between(
                fit["x_fit_log"],
                fit["lower_log"],
                fit["upper_log"],
                color="r",
                alpha=0.18,
                linewidth=0,
            )
        ax.set_title(axis_name, fontsize=11, pad=3)
        ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=9)
        ax.legend(fontsize=7, frameon=False, loc="best")
        if row_index < len(axes) - 1:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("log10 genomic dist (kbp)", fontsize=10)
    fig.supylabel("log10 median dist (nm)", fontsize=10)
    print(f"> Exporting figure to: {saving_filename}")
    plt.savefig(saving_filename, dpi=150, bbox_inches="tight")
    plt.close(fig)


def calculate_physical_vs_genomic_distance(
    input_file: Union[str, Path],
    output_csv: Union[str, Path],
    interloci_csv: Optional[Union[str, Path]] = None,
    plot_file: Optional[Union[str, Path]] = None,
    gen_dist_bins: int = 50,
    dist_threshold: float = np.inf,
    experiment: Optional[str] = None,
    include_3d: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """Load a trace table, calculate R(s) summaries, and write output files."""
    table = ChromatinTraceTable()
    table.load(str(input_file))
    traces = extract_traces_numpy(table)
    gen_dist = compute_genomic_dist_map(table)
    pwd, axis_pwds = compute_all_pwd(traces, show_progress=show_progress)

    results = []
    if include_3d:
        results.append(
            compute_physical_vs_genomic_distance(
                pwd,
                gen_dist,
                gen_dist_bins,
                dist_threshold,
                experiment,
                "3D",
                show_progress,
            )
        )
    for axis, pairwise_dist in axis_pwds.items():
        results.append(
            compute_physical_vs_genomic_distance(
                pairwise_dist,
                gen_dist,
                gen_dist_bins,
                dist_threshold,
                experiment,
                axis,
                show_progress,
            )
        )

    result = pd.concat(results, ignore_index=True)
    result.to_csv(output_csv, index=False)
    if interloci_csv is not None:
        compute_inter_loci_genomic_dist(table).to_csv(interloci_csv, index=False)
    if plot_file is not None:
        plot_file = input_file.split(".")[0] + plot_file
        plot_log_distance_graph(result, plot_file)
    return result


def parse_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="Input chromatin trace table (.ecsv, .dat, .4dn, or .csv).",
    )
    parser.add_argument(
        "--output",
        required=False,
        default="binned_physical_vs_genomic_distances.csv",
        help="Output CSV for binned physical-vs-genomic distances.",
    )
    parser.add_argument(
        "--interloci_output",
        help="Optional CSV for consecutive-locus genomic distances.",
    )
    parser.add_argument(
        "--plot",
        default="_physical_vs_genomic_plot.png",
        help="Optional output filename for a log-log distance plot.",
    )
    parser.add_argument(
        "--gen_dist_bins",
        type=int,
        default=50,
        help="Number of genomic-distance bins. Default: 50",
    )
    parser.add_argument(
        "--dist_threshold",
        type=float,
        default=np.inf,
        help="Maximum physical distance retained before replacing with NaN. Default: infinity",
    )
    parser.add_argument(
        "--experiment",
        help="Experiment or replicate label to include in the output CSV.",
    )
    parser.add_argument(
        "--no_3d",
        action="store_true",
        help="Only calculate X, Y, and Z axis distances.",
    )
    parser.add_argument("--quiet", action="store_true", help="Disable progress bars.")
    return parser


def main() -> None:
    print_script_banner(__file__, __doc__)
    args = parse_arguments().parse_args()
    calculate_physical_vs_genomic_distance(
        input_file=args.input,
        output_csv=args.output,
        interloci_csv=args.interloci_output,
        plot_file=args.plot,
        gen_dist_bins=args.gen_dist_bins,
        dist_threshold=args.dist_threshold,
        experiment=args.experiment,
        include_3d=not args.no_3d,
        show_progress=not args.quiet,
    )
    print(f"Saved physical-vs-genomic distance table to {args.output}")


if __name__ == "__main__":
    main()
