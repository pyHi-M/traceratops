#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will load a localizations table and a trace table to analyze intensity.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.localization_table import LocalizationTable

font = {"weight": "normal", "size": 12}
matplotlib.rc("font", **font)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-L", "--localization", help="Localizations file path")
    parser.add_argument("-T", "--trace", help="Traces file path")
    return parser


def plot_distribution_fluxes(barcode_map, filename_list, trace_ids, nbins=40):

    # Create boolean mask: True if loc inside trace
    buids = barcode_map["Buid"].astype(str)
    trace_ids = np.asarray(trace_ids).astype(str)
    is_in_trace = np.isin(buids, trace_ids)

    colors = np.where(is_in_trace, "red", "blue")

    # initializes figure
    fig, axes = plt.subplots(1, 3)
    ax = axes.ravel()
    fig.set_size_inches((15, 5))

    # initializes variables
    zcentroid = barcode_map["zcentroid"]
    flux = barcode_map["flux"]

    # Bins
    fmin = np.nanmin(flux)
    fmax = np.nanmax(flux)
    if np.isclose(fmin, fmax):
        # degenerate case: a single flux level -> creates 2 bins around it
        fmin, fmax = fmin - 0.5, fmax + 0.5
    bins = np.linspace(fmin, fmax, nbins + 1)
    centers = 0.5 * (bins[1:] + bins[:-1])

    # Comptes
    counts_in, _ = np.histogram(flux[is_in_trace], bins=bins)
    counts_out, _ = np.histogram(flux[~is_in_trace], bins=bins)
    count_all, _ = np.histogram(flux, bins=bins)

    # plots data
    ax[0].plot(centers, counts_in, label="Inside trace", linewidth=2)
    ax[0].plot(centers, counts_out, label="Off trace", linewidth=2)
    ax[0].set_xlabel("flux")
    ax[0].set_ylabel("Localization count")
    ax[0].set_yscale("log")
    ax[0].set_title("Distribution of localizations by flux")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    # ---------- POINT D'INTERSECTION ----------
    diff = counts_in - counts_out
    sign_change = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]

    if len(sign_change) > 0:
        idx = sign_change[0]
        # interpolation linéaire pour une valeur plus précise
        x0, x1 = centers[idx], centers[idx + 1]
        y0, y1 = diff[idx], diff[idx + 1]
        cross_x = x0 - y0 * (x1 - x0) / (y1 - y0)
        cross_y = counts_in[idx] + (counts_in[idx + 1] - counts_in[idx]) * (
            (cross_x - x0) / (x1 - x0)
        )

        ax[0].axvline(cross_x, color="black", linestyle="--", alpha=0.6)
        ax[0].scatter([cross_x], [cross_y], color="black", zorder=5)
        ax[0].annotate(
            f"Cross: {cross_x:.2f}",
            xy=(cross_x, cross_y),
            xytext=(cross_x, cross_y * 1.5),
            arrowprops=dict(facecolor="black", arrowstyle="->"),
            fontsize=10,
        )

    ax[1].scatter(flux, zcentroid, c=colors, alpha=0.2, s=5)
    ax[1].set_title("Red: Inside trace | Blue: off trace")
    ax[1].set_xlabel("flux")
    ax[1].set_ylabel("zcentroid")

    # ---------- PLOT 3 : sum ----------
    # sum_in = np.cumsum(counts_in)
    # sum_out = np.cumsum(counts_out)

    ax[2].plot(centers, count_all, label="Inside trace", linewidth=2)
    ax[2].set_xlabel("flux")
    ax[2].set_ylabel("Localization count")
    # ax[2].set_yscale("log")
    ax[2].set_title("Distribution of all localizations")
    ax[2].grid(True, alpha=0.3)

    # ---------- ELBOW METHOD ----------
    # Normalisation pour éviter l'effet d'échelle
    y = count_all.astype(float)
    x = centers.astype(float)

    # On suppose que le "coude" est l'endroit où la distance à la ligne reliant
    # le premier et le dernier point est maximale
    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])
    distances = []
    for i in range(len(x)):
        p = np.array([x[i], y[i]])
        # distance point-ligne (p1,p2)
        d = np.abs(np.cross(p2 - p1, p1 - p)) / np.linalg.norm(p2 - p1)
        distances.append(d)
    distances = np.array(distances)
    elbow_idx = np.argmax(distances)
    elbow_x = x[elbow_idx]
    elbow_y = y[elbow_idx]

    # Tracé de la courbe de distance (pour visualisation du coude)
    ax2 = ax[2].twinx()
    ax2.plot(
        x, distances, color="orange", linestyle="--", alpha=0.7, label="Elbow distance"
    )
    ax2.set_ylabel("Distance to baseline")
    ax2.legend(loc="upper right")

    # Annotation de la valeur du coude sur le graphe principal
    ax[2].axvline(elbow_x, color="red", linestyle="--", alpha=0.6)
    ax[2].scatter([elbow_x], [elbow_y], color="red", zorder=5)
    ax[2].annotate(
        f"Elbow: {elbow_x:.2f}",
        xy=(elbow_x, elbow_y),
        xytext=(elbow_x, elbow_y * 1.1),
        arrowprops=dict(facecolor="black", arrowstyle="->"),
        fontsize=10,
    )

    # saves figure
    fig.tight_layout()
    filename = "".join(filename_list)
    print(f"save to : {filename}")
    fig.savefig(filename)
    plt.close(fig)


def get_roi(trace_file="Trace_3D_barcode_mask-mask1_ROI-13.ecsv"):
    basename = trace_file.split(".")[-2]
    return basename.split("_")[-1]


def process_intensities(loc_file, trace_file):
    loc_table = LocalizationTable()
    barcode_map, _ = loc_table.load(loc_file)
    print(f"> Analyzing localizations for {loc_file}")

    trace_table = ChromatinTraceTable()
    trace_data = trace_table.load(trace_file)
    print(f"> Analyzing traces for {trace_file}")

    # Récupération des IDs de traces
    trace_ids = trace_data["Spot_ID"]
    roi = get_roi(trace_file)
    localization_stats_file = [
        "localization_table_stats_",
        roi,
        ".png",
    ]
    plot_distribution_fluxes(barcode_map, localization_stats_file, trace_ids)


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    # [loops over lists of datafolders]
    process_intensities(args.localization, args.trace)
    print("Finished execution")


if __name__ == "__main__":
    main()
