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

    # plots data
    ax[0].plot(centers, counts_in, label="Inside trace", linewidth=2)
    ax[0].plot(centers, counts_out, label="Off trace", linewidth=2)
    ax[0].set_xlabel("flux")
    ax[0].set_ylabel("Localization count")
    ax[0].set_yscale("log")
    ax[0].set_title("Distribution of localizations by flux")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()

    ax[1].scatter(flux, zcentroid, c=colors, alpha=0.2, s=5)
    ax[1].set_title("Red: Inside trace | Blue: off trace")
    ax[1].set_xlabel("flux")
    ax[1].set_ylabel("zcentroid")

    # ---------- PLOT 3 : sum ----------
    sum_in = np.cumsum(counts_in)
    sum_out = np.cumsum(counts_out)

    ax[2].plot(centers, sum_in, label="Inside trace", linewidth=2)
    ax[2].plot(centers, sum_out, label="Off trace", linewidth=2)
    ax[2].set_xlabel("flux")
    ax[2].set_ylabel("Cumulative number")
    ax[2].set_title("Cumul localizations")
    ax[2].grid(True, alpha=0.3)
    ax[2].legend()

    # saves figure
    fig.tight_layout()
    fig.savefig("".join(filename_list))
    plt.close(fig)


def process_intensities(loc_file, trace_file):
    loc_table = LocalizationTable()
    barcode_map, _ = loc_table.load(loc_file)
    print(f"> Analyzing localizations for {loc_file}")

    trace_table = ChromatinTraceTable()
    trace_data = trace_table.load(trace_file)
    print(f"> Analyzing traces for {trace_file}")

    # Récupération des IDs de traces
    trace_ids = trace_data["Spot_ID"]

    localization_stats_file = [
        loc_file.split(".")[0],
        "localization_table_stats",
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
