#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot first/second moments and Gaussianity index from trace pairwise relative positions."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to input trace file (.ecsv/.dat/.4dn)")
    parser.add_argument("-O", "--output", default="plots", help="Output folder")
    parser.add_argument("--plot_format", default="png", choices=["png", "pdf", "svg"], help="Output plot format")
    parser.add_argument(
        "--avg_dims",
        default="none",
        choices=["none", "xy", "xyz"],
        help="Average GaussianityIndex over dimensions: none, xy, or xyz",
    )
    parser.add_argument(
        "--coords",
        default="xyz",
        choices=["xyz", "zyx"],
        help="Coordinate order used in trace file. xyz means columns x,y,z. zyx swaps first and last axes.",
    )
    return parser


def build_conformations(trace_table_data, coord_order="xyz"):
    barcodes = np.array(sorted({int(b) for b in trace_table_data["Barcode #"]}))
    traces = np.array(sorted({str(t) for t in trace_table_data["Trace_ID"]}))

    b_to_i = {b: i for i, b in enumerate(barcodes)}
    t_to_i = {t: i for i, t in enumerate(traces)}

    conf = np.zeros((len(traces), len(barcodes), 3), dtype=float)
    conf_mask = np.zeros_like(conf, dtype=float)

    x = np.array(trace_table_data["x"], dtype=float)
    y = np.array(trace_table_data["y"], dtype=float)
    z = np.array(trace_table_data["z"], dtype=float)
    if coord_order == "zyx":
        coords = np.stack([z, y, x], axis=1)
    else:
        coords = np.stack([x, y, z], axis=1)

    for row_idx in range(len(trace_table_data)):
        tid = str(trace_table_data["Trace_ID"][row_idx])
        bid = int(trace_table_data["Barcode #"][row_idx])
        ti = t_to_i[tid]
        bi = b_to_i[bid]
        conf[ti, bi, :] = coords[row_idx]
        conf_mask[ti, bi, :] = 1.0

    return conf, conf_mask, barcodes


def compute_moments_and_gaussianity(conf, conf_mask):
    pair_vect = conf[:, :, None, :] - conf[:, None, :, :]
    pair_mask = conf_mask[:, :, None, :] * conf_mask[:, None, :, :]

    valid = np.sum(pair_mask, axis=0)
    valid_safe = np.where(valid > 0, valid, np.nan)

    mean_rel = np.sum(pair_vect * pair_mask, axis=0) / valid_safe
    second_moment = np.sum((pair_vect * pair_mask) ** 2, axis=0) / valid_safe
    variance = np.maximum(second_moment - mean_rel**2, 0.0)

    fourth_moment = np.sum((pair_vect * pair_mask) ** 4, axis=0) / valid_safe
    gaussianity = np.divide(
        3.0 * (second_moment**2),
        fourth_moment,
        out=np.full_like(fourth_moment, np.nan),
        where=fourth_moment > 0,
    )

    return mean_rel, variance, gaussianity


def safe_nanmean_last_axis(arr):
    valid = np.sum(~np.isnan(arr), axis=-1)
    summed = np.nansum(arr, axis=-1)
    return np.divide(summed, valid, out=np.full_like(summed, np.nan, dtype=float), where=valid > 0)


def summarize_to_scalar_map(mean_rel, variance, gaussianity, avg_dims):
    mean_distance = np.linalg.norm(mean_rel, axis=-1)
    variance_total = np.sum(variance, axis=-1)

    if avg_dims == "xy":
        gaussianity_map = safe_nanmean_last_axis(gaussianity[:, :, :2])
    elif avg_dims == "xyz":
        gaussianity_map = safe_nanmean_last_axis(gaussianity)
    else:
        gaussianity_map = gaussianity

    return mean_distance, variance_total, gaussianity_map


def _barcode_ticks(barcodes):
    n = len(barcodes)
    if n <= 15:
        idx = np.arange(n)
    else:
        idx = np.linspace(0, n - 1, 6, dtype=int)
    return idx, [str(int(barcodes[i])) for i in idx]


def plot_maps(mean_distance, variance, gaussianity, barcodes, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

    items = [
        (mean_distance, "Distance matrix", "Mean distance"),
        (variance, "Variance matrix", "Variance"),
        (gaussianity, "Gaussianity matrix", "Gaussianity index"),
    ]

    ticks, labels = _barcode_ticks(barcodes)

    for ax, (mat, title, cbar_label) in zip(axes, items):
        if mat.ndim == 3:
            mat = safe_nanmean_last_axis(mat)
        im = ax.imshow(mat, cmap="coolwarm")
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


def main():
    args = parse_arguments().parse_args()
    os.makedirs(args.output, exist_ok=True)

    trace_table = ChromatinTraceTable()
    trace_table.load(args.input)

    conf, conf_mask, barcodes = build_conformations(trace_table.data, coord_order=args.coords)
    mean_rel, variance, gaussianity = compute_moments_and_gaussianity(conf, conf_mask)
    mean_distance, variance_total, gaussianity_map = summarize_to_scalar_map(
        mean_rel, variance, gaussianity, args.avg_dims
    )

    stem = os.path.splitext(os.path.basename(args.input))[0]
    np.save(os.path.join(args.output, f"{stem}_mean_distance.npy"), mean_distance)
    np.save(os.path.join(args.output, f"{stem}_variance.npy"), variance_total)
    np.save(os.path.join(args.output, f"{stem}_gaussianity.npy"), gaussianity_map)

    out_plot = os.path.join(args.output, f"{stem}_pwd_moments.{args.plot_format}")
    plot_maps(mean_distance, variance_total, gaussianity_map, barcodes, out_plot)
    print(f"Saved: {out_plot}")


if __name__ == "__main__":
    main()
