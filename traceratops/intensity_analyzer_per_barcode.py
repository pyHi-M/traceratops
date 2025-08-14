#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Split intensity plots by 'Barcode #' (from LocalizationTable).
For each Barcode #, generates a 3-panel figure:
1) Histogram counts of flux (Inside trace vs Off trace) + crossing point
2) Scatter zcentroid vs flux, colored by membership
3) Histogram of all localizations + elbow method visualization

Inspired by the original script that plots aggregate distributions:contentReference[oaicite:1]{index=1}.
"""

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from traceratops.core.chromatin_trace_table import (  # :contentReference[oaicite:2]{index=2}
    ChromatinTraceTable,
)
from traceratops.core.localization_table import (  # :contentReference[oaicite:3]{index=3}
    LocalizationTable,
)

font = {"weight": "normal", "size": 12}
matplotlib.rc("font", **font)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-L", "--localization", required=True, help="Localizations file path"
    )
    parser.add_argument("-T", "--trace", required=True, help="Traces file path")
    parser.add_argument(
        "-o", "--outdir", default="plots_by_barcode", help="Output directory"
    )
    parser.add_argument(
        "--nbins", type=int, default=40, help="Number of histogram bins"
    )
    return parser


def get_roi(trace_file):
    """Mimic ROI extraction of the source script:contentReference[oaicite:4]{index=4}."""
    basename = trace_file.split(".")[-2]
    return basename.split("_")[-1]


def compute_crossing_x(centers, counts_in, counts_out):
    """Find first crossing between counts_in and counts_out via sign change + linear interpolation."""
    diff = counts_in - counts_out
    # ignore bins with both zeros to avoid spurious sign changes
    mask_valid = ~((counts_in == 0) & (counts_out == 0))
    diff = diff.copy()
    diff[~mask_valid] = 0

    sign_change = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    if len(sign_change) == 0:
        return None, None

    idx = sign_change[0]
    x0, x1 = centers[idx], centers[idx + 1]
    y0, y1 = diff[idx], diff[idx + 1]
    if (y1 - y0) == 0 or (x1 - x0) == 0:
        return None, None

    cross_x = x0 - y0 * (x1 - x0) / (y1 - y0)

    # interpolate counts_in at cross_x for plotting the marker
    ci0, ci1 = counts_in[idx], counts_in[idx + 1]
    cross_y = ci0 + (ci1 - ci0) * ((cross_x - x0) / (x1 - x0))
    return cross_x, cross_y


def elbow_knee_x(centers, values):
    """
    Elbow via max distance to straight line between endpoints (Kneedle-like geometric heuristic).
    Works best for monotonically decreasing-ish histograms; robust enough for our purpose.
    """
    x = centers.astype(float)
    y = values.astype(float)
    if len(x) < 2:
        return None, None, None

    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])
    denom = np.linalg.norm(p2 - p1)
    if denom == 0:
        return None, None, None

    # point-to-line distances
    vec = p2 - p1
    distances = np.abs(np.cross(vec, (np.vstack([x, y]).T - p1))) / denom
    elbow_idx = int(np.argmax(distances))
    return x[elbow_idx], y[elbow_idx], distances


def plot_one_barcode(df_barcode, trace_ids, nbins, title_suffix, outpath):
    """
    Reproduce 3-panel figure for a single Barcode #:
    0) histogram counts_in/out + crossing
    1) scatter flux vs zcentroid
    2) histogram all + elbow
    """
    # Membership to traces via Buid ~ Spot_ID (as in source code):contentReference[oaicite:5]{index=5}
    buids = np.array(df_barcode["Buid"]).astype(str)
    trace_ids = np.asarray(trace_ids).astype(str)
    is_in_trace = np.isin(buids, trace_ids)

    colors = np.where(is_in_trace, "red", "blue")

    zcentroid = np.array(df_barcode["zcentroid"])
    flux = np.array(df_barcode["flux"])

    # Bins
    fmin = np.nanmin(flux)
    fmax = np.nanmax(flux)
    if not np.isfinite(fmin) or not np.isfinite(fmax):
        return  # skip empty/invalid subset
    if np.isclose(fmin, fmax):
        fmin, fmax = fmin - 0.5, fmax + 0.5
    bins = np.linspace(fmin, fmax, nbins + 1)
    centers = 0.5 * (bins[1:] + bins[:-1])

    counts_in, _ = np.histogram(flux[is_in_trace], bins=bins)
    counts_out, _ = np.histogram(flux[~is_in_trace], bins=bins)
    count_all, _ = np.histogram(flux, bins=bins)

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ax = axes.ravel()

    # --- Plot 1: two distributions + crossing point (inspired by source fig[0]):contentReference[oaicite:6]{index=6}
    ax[0].plot(centers, counts_in, label="Inside trace", linewidth=2)
    ax[0].plot(centers, counts_out, label="Off trace", linewidth=2)
    ax[0].set_xlabel("flux")
    ax[0].set_ylabel("Localization count")
    ax[0].set_yscale("log")
    ax[0].grid(True, alpha=0.3)
    ax[0].set_title(f"Distribution by flux — {title_suffix}")
    ax[0].legend()

    cx, cy = compute_crossing_x(centers, counts_in, counts_out)
    if cx is not None:
        ax[0].axvline(cx, color="black", linestyle="--", alpha=0.6)
        if np.isfinite(cy):
            ax[0].scatter([cx], [max(cy, 1e-9)], color="black", zorder=5)
            annot_y = max(cy * 1.5, cy + 1.0)
        else:
            annot_y = ax[0].get_ylim()[1] * 0.7
        ax[0].annotate(
            f"Cross: {cx:.3f}",
            xy=(cx, annot_y),
            xytext=(cx, annot_y),
            textcoords="data",
            fontsize=10,
        )

    # --- Plot 2: scatter (inspired by source fig[1]):contentReference[oaicite:7]{index=7}
    ax[1].scatter(flux, zcentroid, c=colors, alpha=0.2, s=5)
    ax[1].set_title("Red: Inside trace | Blue: Off trace")
    ax[1].set_xlabel("flux")
    ax[1].set_ylabel("zcentroid")

    # --- Plot 3: all counts + elbow (inspired by source fig[2]):contentReference[oaicite:8]{index=8}
    ax[2].plot(centers, count_all, label="All localizations", linewidth=2)
    ax[2].set_xlabel("flux")
    ax[2].set_ylabel("Localization count")
    ax[2].grid(True, alpha=0.3)
    ax[2].set_title("All localizations + elbow")

    ex, ey, distances = elbow_knee_x(centers, count_all)
    if ex is not None:
        # overlay distances on twin axis
        ax2 = ax[2].twinx()
        ax2.plot(centers, distances, linestyle="--", alpha=0.7, label="Elbow distance")
        ax2.set_ylabel("Distance to baseline")
        # mark elbow
        ax[2].axvline(ex, color="red", linestyle="--", alpha=0.6)
        if np.isfinite(ey):
            ax[2].scatter([ex], [ey], color="red", zorder=5)
            txt_y = ey * 1.1 if ey > 0 else (ax[2].get_ylim()[1] * 0.1)
        else:
            txt_y = ax[2].get_ylim()[1] * 0.1
        ax[2].annotate(
            f"Elbow: {ex:.3f}",
            xy=(ex, txt_y),
            xytext=(ex, txt_y),
            textcoords="data",
            fontsize=10,
        )

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    args = parse_arguments().parse_args()

    # Load tables (same API as source script):contentReference[oaicite:9]{index=9}
    loc_table = LocalizationTable()
    barcode_map, _ = loc_table.load(args.localization)
    print(f"> Analyzing localizations for {args.localization}")

    trace_table = ChromatinTraceTable()
    trace_data = trace_table.load(args.trace)
    print(f"> Analyzing traces for {args.trace}")

    trace_ids = trace_data[
        "Spot_ID"
    ]  # used to flag Inside trace:contentReference[oaicite:10]{index=10}
    roi = get_roi(args.trace)

    # Ensure expected columns exist
    required_cols = {"Barcode #", "Buid", "zcentroid", "flux"}
    missing = required_cols - set(barcode_map.columns)
    if missing:
        raise KeyError(f"Missing required columns in barcode_map: {missing}")

    # Group by Barcode # and plot one figure per group
    # for barcode_value, df_barcode in barcode_map.group_by("Barcode #"):
    grouped = barcode_map.group_by("Barcode #")
    for key_vals, group in zip(grouped.groups.keys, grouped.groups):
        barcode_value = key_vals["Barcode #"]
        df_barcode = group
        safe_barcode = str(barcode_value).replace("/", "_").replace(" ", "_")
        outname = f"localization_stats_{roi}_barcode-{safe_barcode}.png"
        outpath = os.path.join(args.outdir, outname)
        title_suffix = f"Barcode #{barcode_value}"
        print(f"> Plotting {title_suffix} -> {outpath}")
        plot_one_barcode(df_barcode, trace_ids, args.nbins, title_suffix, outpath)

    print("Finished execution")


if __name__ == "__main__":
    main()
