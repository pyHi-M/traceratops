#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot first/second moments and Gaussianity index from trace pairwise relative positions.

For each pair of barcode, you calculate the distribution of positions of one barcode relatively to the other one (so that gives you a 3D distribution centered on 0).
Then, from each distribution and along each spatial dimension, you just calculate 3 / kurtosis. 

It equals 1 for a Gaussian and < 1 for a mixture of Gaussian. kurtosis = fourth central moment divided by square variance
And that's it => You get a matrix for each dimension. You can optionally average them over x and y (and z if it's not too noisy).


"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import pairwise_distances
from astropy.table import unique
from tqdm.contrib import tzip
from scipy.stats import kurtosis
from scipy.spatial.distance import cdist
from tqdm.auto import tqdm

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to input trace file (.ecsv/.dat/.4dn)")
    parser.add_argument("-O", "--output", default="plots", help="Output folder")
    parser.add_argument("--plot_format", default="png", choices=["png", "pdf", "svg"], help="Output plot format")
    parser.add_argument(
        "--avg_dims",
        default="xyz",
        choices=["xy", "xyz"],
        help="Average GaussianityIndex over dimensions: xy, or xyz",
    )

    return parser

def _barcode_ticks(barcodes):
    n = len(barcodes)
    if n <= 15:
        idx = np.arange(n)
    else:
        idx = np.linspace(0, n - 1, 6, dtype=int)
    return idx, [str(barcodes[i]) for i in idx]


def plot_maps(mean_distance, variance, gaussianity, barcodes, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

    items = [
        (mean_distance, "Distance matrix", "Mean distance"),
        (variance, "Variance matrix", "Variance"),
        (gaussianity, "Gaussianity matrix", "Gaussianity index"),
    ]

    ticks, labels = _barcode_ticks(barcodes)

    for i, (ax, (mat, title, cbar_label)) in enumerate(zip(axes, items)):
        if mat.ndim == 3:
            mat = safe_nanmean_last_axis(mat)

        if i == 2:
            im = ax.imshow(mat, cmap="RdBu", vmin=0, vmax=1)
        else:
            im = ax.imshow(mat, cmap="RdBu")

        ax.set_title(title)
        ax.set_xlabel("Barcode #")
        ax.set_ylabel("Barcode #")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label)
        

    fig.savefig(out_png, dpi=200)
    plt.close(fig)



def build_distance_matrix_fast_min(trace_table, distance_threshold=np.inf, dtype=np.float32, coord_order='xyz'):
    """
    Fast version of build_distance_matrix for mode='min' only.
    """

    data = trace_table.data

    unique_trace_ids = unique(data, keys="Trace_ID")["Trace_ID"].data
    number_matrices = len(unique_trace_ids)

    unique_barcodes = unique(data, keys="Barcode #")["Barcode #"].data
    number_unique_barcodes = len(unique_barcodes)

    print(f"$ Found {number_unique_barcodes} barcodes and {number_matrices} traces.", "INFO")

    # Fast barcode -> matrix index mapping
    barcode_to_index = {bc: i for i, bc in enumerate(unique_barcodes)}

    sc_matrix = np.full(
        (number_unique_barcodes, number_unique_barcodes, number_matrices),
        np.nan,
        dtype=dtype,
    )

    data_traces = data.group_by("Trace_ID")

    print("> Processing traces...", "INFO")

    for itrace, trace in enumerate(tqdm(data_traces.groups, total=number_matrices)):

        barcodes = np.asarray(trace["Barcode #"].data)
        barcode_indices = np.array([barcode_to_index[bc] for bc in barcodes])

        if coord_order=='xyz':
            coords = np.column_stack([
                np.asarray(trace["x"].data, dtype=dtype),
                np.asarray(trace["y"].data, dtype=dtype),
                np.asarray(trace["z"].data, dtype=dtype),
            ])
        elif coord_order=='xy':
            coords = np.column_stack([
                np.asarray(trace["x"].data, dtype=dtype),
                np.asarray(trace["y"].data, dtype=dtype),
                #np.asarray(trace["z"].data, dtype=dtype),
            ])

        pwd_matrix = cdist(coords, coords)

        # Ignore same barcode pairs
        valid = barcode_indices[:, None] != barcode_indices[None, :]

        # Apply distance threshold
        valid &= pwd_matrix < distance_threshold

        if not np.any(valid):
            continue

        ii, jj = np.where(valid)

        row_idx = barcode_indices[ii]
        col_idx = barcode_indices[jj]
        distances = pwd_matrix[ii, jj].astype(dtype)

        # Efficiently assign the minimum distance for duplicate barcode combinations
        current = sc_matrix[:, :, itrace]

        # np.minimum.at does not handle NaNs as desired, so initialize selected NaNs to inf
        missing = np.isnan(current[row_idx, col_idx])
        current[row_idx[missing], col_idx[missing]] = np.inf

        np.minimum.at(current, (row_idx, col_idx), distances)

        # Convert untouched inf values back to NaN
        current[np.isinf(current)] = np.nan

    return sc_matrix, unique_barcodes

def get_mean_distance_map(trace_table, coord_order='xyz'):

    print(f"$ Number of spots in trace file: {len(trace_table.data)}")
 
    sc_matrix, barcodes = build_distance_matrix_fast_min(trace_table,coord_order=coord_order)

    sc_matrix_mean_distance = np.nanmean(sc_matrix, axis=2)

    sc_matrix_var_distance = np.nanvar(sc_matrix, axis=2)

    sc_matrix_kurtosis = kurtosis(
        sc_matrix,
        axis=2,
        nan_policy='omit',   # ignore NaNs
        fisher=False          # 0 = Gaussian baseline (recommended)
    )

    print(f"$ Number of bins in map : {sc_matrix.shape}")

    return sc_matrix_mean_distance, sc_matrix_var_distance, sc_matrix_kurtosis, barcodes

def main():
    args = parse_arguments().parse_args()
    os.makedirs(args.output, exist_ok=True)

    trace_table = ChromatinTraceTable()
    trace_table.load(args.input)

    sc_matrix_mean_distance, sc_matrix_var_distance, sc_matrix_kurtosis, barcodes = get_mean_distance_map(
        trace_table, coord_order=args.avg_dims
    )
    gaussianity = 3/sc_matrix_kurtosis

    stem = os.path.splitext(os.path.basename(args.input))[0]
    np.save(os.path.join(args.output, f"{stem}_mean_distance.npy"), sc_matrix_mean_distance)
    np.save(os.path.join(args.output, f"{stem}_variance.npy"), sc_matrix_var_distance)
    np.save(os.path.join(args.output, f"{stem}_gaussianity.npy"), gaussianity)

    out_plot = os.path.join(args.output, f"{stem}_pwd_moments.{args.plot_format}")
    plot_maps(sc_matrix_mean_distance, sc_matrix_var_distance, gaussianity, barcodes, out_plot)

    print(f"Saved: {out_plot}")



if __name__ == "__main__":
    main()
