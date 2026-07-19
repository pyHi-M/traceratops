#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate physical distances as a function of genomic distances from a chromatin trace table."""

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import unique
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


def plot_log_distance_graph(
    dist_df: pd.DataFrame, saving_filename: Union[str, Path]
) -> None:
    """Plot log10 median physical distance versus log10 genomic distance."""
    import seaborn as sns

    axes = [axis for axis in ["3D", "X", "Y", "Z"] if axis in set(dist_df["axis"])]
    fig, axs = plt.subplots(
        nrows=len(axes), ncols=1, squeeze=False, figsize=(8, 4 * len(axes))
    )
    for row_index, axis_name in enumerate(axes):
        df = dist_df.loc[dist_df["axis"] == axis_name]
        sns.scatterplot(
            data=df,
            x="log10 genomic dist (kbp)",
            y="log10 median dist (nm)",
            hue="experiment",
            ax=axs[row_index, 0],
        )
        axs[row_index, 0].set_title(axis_name)
        if len(dist_df["experiment"].unique()) == 1:
            legend = axs[row_index, 0].get_legend()
            if legend is not None:
                legend.remove()
    print(f"> Exporting figure to: {saving_filename}")
    plt.savefig(saving_filename, dpi=100, bbox_inches="tight")
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
