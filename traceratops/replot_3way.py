#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Replot 3-way co-localization matrices from .npy files.

It visualizes the frequency of co-localization between barcodes with reference to an anchor barcode, highlighting the anchor's position with perpendicular lines on the heatmap.
"""
import argparse
import csv
import os
import re
import select
import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", nargs="*", help="List of input .npy matrix files.", default=[]
    )
    parser.add_argument("--pipe", action="store_true", help="Read filenames from stdin")
    parser.add_argument(
        "--vmin", type=float, default=None, help="Minimum color scale value"
    )
    parser.add_argument(
        "--vmax", type=float, default=None, help="Maximum color scale value"
    )
    parser.add_argument(
        "--cmap", default="RdBu", help="Matplotlib colormap (default: RdBu)"
    )
    parser.add_argument(
        "--label_map_file", help="Text file with barcode numbers per row"
    )
    return parser


def plot_threeway_matrix(
    pair_means,
    anchor_barcode,
    output_file="3way",
    distance_cutoff=0.2,
    vmin=None,
    vmax=None,
    cmap="RdBu",
    label_map=None,
):
    all_barcodes = set()
    for b1, b2 in pair_means.keys():
        all_barcodes.add(b1)
        all_barcodes.add(b2)
    sorted_barcodes = sorted(all_barcodes)
    n_barcodes = len(sorted_barcodes)
    mean_matrix = np.zeros((n_barcodes, n_barcodes))
    barcode_to_idx = {b: i for i, b in enumerate(sorted_barcodes)}
    idx_to_barcode = {
        i: label_map[str(i)] if label_map else str(b)
        for i, b in enumerate(sorted_barcodes)
    }
    for (b1, b2), mean_val in pair_means.items():
        i, j = barcode_to_idx[b1], barcode_to_idx[b2]
        mean_matrix[i, j] = mean_val
        mean_matrix[j, i] = mean_val
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(
        mean_matrix,
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin if vmin is not None else 0,
        vmax=(
            vmax
            if vmax is not None
            else (0.9 * mean_matrix.max() if mean_matrix.max() > 0 else 1)
        ),
    )
    tick_labels = [idx_to_barcode[i] for i in range(n_barcodes)]
    ax.set_xticks(np.arange(n_barcodes))
    ax.set_yticks(np.arange(n_barcodes))
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_yticklabels(tick_labels, fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    ax.set_xticks(np.arange(-0.5, n_barcodes, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_barcodes, 1), minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1)
    previous_anchor_barcode_list = [
        int(label_map[x]) - int(anchor_barcode) for x in label_map.keys()
    ]
    # Keep only values that are less than 0 (i.e., label_map[x] < anchor_barcode)
    negative_differences = [
        (i, diff) for i, diff in enumerate(previous_anchor_barcode_list) if diff < 0
    ]
    if negative_differences:
        index_closest_to_zero = min(negative_differences, key=lambda x: abs(x[1]))[0]
        previous_anchor_barcode = label_map[str(index_closest_to_zero)]
    else:
        previous_anchor_barcode = (
            None  # Or handle the "no lower value" case appropriately
        )
    print(f"> previous_anchor_barcode: {previous_anchor_barcode}")
    if str(previous_anchor_barcode) in idx_to_barcode.values():
        anchor_idx = list(idx_to_barcode.values()).index(str(previous_anchor_barcode))
        ax.axhline(
            y=anchor_idx + 0.5, color="black", linestyle="-", linewidth=2, alpha=0.7
        )
        ax.axvline(
            x=anchor_idx + 0.5, color="black", linestyle="-", linewidth=2, alpha=0.7
        )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Co-localization frequency", fontsize=12)
    ax.set_title(
        f"3-way co-localization with anchor {anchor_barcode}\n(distance cutoff: {distance_cutoff} µm)",
        fontsize=14,
    )
    ax.set_xlabel("Barcode #", fontsize=14)
    ax.set_ylabel("Barcode #", fontsize=14)
    plt.tight_layout()
    output_filename = f"{output_file.split('.')[0]}_anchor_{anchor_barcode}_replot"
    plt.savefig(f"{output_filename}.png", dpi=300)
    print(f"Saved three-way co-localization heatmap to: {output_filename}")
    plt.close()


def extract_anchor(filename):
    match = re.search(r"anchor_(\d+)", filename)
    return int(match.group(1)) if match else -1


def load_list(file_name, anchor):
    with open(file_name, newline="", encoding="utf-8") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        idx_to_barcode = {}
        index = 0
        for row in spamreader:
            if anchor != int(row[0]):
                idx_to_barcode[str(index)] = row[0]
                index += 1
    return idx_to_barcode


def load_label_map(filepath, anchor):
    try:
        idx_to_barcode = load_list(filepath, anchor)
        return idx_to_barcode
    except Exception as e:
        print(f"Error loading barcode list file: {e}")
        return None


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    matrix_files = list(args.input)
    if args.pipe and select.select([sys.stdin], [], [], 0.0)[0]:
        matrix_files += [
            line.strip() for line in sys.stdin if line.strip().endswith(".npy")
        ]
    if not matrix_files:
        print("No matrix files provided.")
        return
    for npy_file in matrix_files:
        if not os.path.exists(npy_file):
            print(f"File not found: {npy_file}")
            continue
        print(f"Loading matrix: {npy_file}")
        matrix = np.load(npy_file)
        anchor = extract_anchor(npy_file)
        print(f"$ anchor: {anchor}")
        label_map = (
            load_label_map(args.label_map_file, anchor) if args.label_map_file else None
        )
        if label_map is not None:
            # Remove anchor from label_map
            for k, v in list(label_map.items()):
                if v == str(anchor):
                    label_map.pop(k)
                    break
        if anchor == -1:
            print(f"Could not determine anchor from filename: {npy_file}")
            continue
        n = matrix.shape[0]
        barcodes = list(range(n))
        pair_means = {}
        for i in range(n):
            for j in range(i + 1, n):
                pair_means[(barcodes[i], barcodes[j])] = matrix[i, j]
        output_base = os.path.splitext(npy_file)[0]
        plot_threeway_matrix(
            pair_means,
            anchor,
            output_file=output_base,
            distance_cutoff=0.2,
            vmin=args.vmin,
            vmax=args.vmax,
            cmap=args.cmap,
            label_map=label_map,
        )


if __name__ == "__main__":
    main()
