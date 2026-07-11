#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load a trace file and a number of numpy masks and assign them labels
"""

import argparse
import os
import select
import sys

import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.script_banner import print_script_banner


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Input trace file")
    parser.add_argument(
        "--mask_file", help="Input mask image file. Expected format: NPY"
    )
    parser.add_argument(
        "--pixel_size", help="Lateral pixel size un microns. Default = 0.1"
    )
    parser.add_argument("--label", help="Label to add to trace file. Default=labeled")

    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument(
        "--output_format",
        choices=["png", "svg", "pdf"],
        default="png",
        help="Output image format. Default = png.",
    )

    return parser


def create_dict_args(args):
    p = {}
    p["trace_files"] = []
    if args.input:
        p["trace_files"].append(args.input)
    if args.mask_file:
        p["mask_file"] = args.mask_file
    else:
        print(">> ERROR: you must provide a filename with a mask file")
        sys.exit(-1)
    if args.pixel_size:
        p["pixel_size"] = args.pixel_size
    else:
        p["pixel_size"] = 0.1
    if args.label:
        p["label"] = args.label
    else:
        p["label"] = "labeled"
    p["output_format"] = args.output_format
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
    if len(p["trace_files"]) < 1:
        print(">> ERROR: you must provide a filename with a trace file")
        sys.exit(-1)
    return p


def load_mask(mask_file):
    # accept both NPY and TIFF formats for the mask file
    # but have to be in 2d
    if mask_file.endswith(".npy"):
        data_2d = np.load(mask_file, allow_pickle=False).squeeze()
    elif mask_file.endswith(".tif") or mask_file.endswith(".tiff"):
        data_2d = tf.imread(mask_file)
    else:
        print(
            f"ERROR: mask file format not supported. Please provide a NPY or TIFF file. Provided file: {mask_file}"
        )
        sys.exit(-1)
    if data_2d.ndim != 2:
        print(
            f"ERROR: mask file must be 2D. Provided file has {data_2d.ndim} dimensions."
        )
        sys.exit(-1)
    print(f"$ mask image file read: {mask_file}")
    return data_2d


def assign_masks(trace, mask_file, label="labeled", pixel_size=0.1):
    # [checks if mask file exists for the file to process]
    if os.path.exists(mask_file):
        # load mask
        data_2d = load_mask(mask_file)

        # === ensure label column is large enough ===
        existing_labels = [str(v) for v in trace.data["label"]]
        max_existing_len = (
            max(len(v) for v in existing_labels) if existing_labels else 1
        )
        new_max_len = max_existing_len + 1 + len(label)  # comma + new label
        trace.data["label"] = trace.data["label"].astype(f"<U{new_max_len}")

        # matches traces and masks
        index = 0
        labeled_trace = []

        for trace_row in trace.data:
            x_int = int(trace_row["x"] / pixel_size)
            y_int = int(trace_row["y"] / pixel_size)

            if "x" in trace_row["label"]:
                trace_row["label"] = "_"

            # labels are appended as comma separated lists.
            # Thus a trace can have multiple labels
            if data_2d[x_int, y_int] >= 1:
                trace_row["label"] = str(trace_row["label"]) + "," + label
                index += 1
                labeled_trace.append(trace_row["Trace_ID"])

        unique_traces_labeled = set(labeled_trace)
        print(
            f"\n> {index} trace rows out of {len(trace.data)} were associated to mask {label}. Unique traces: {len(unique_traces_labeled)}"
        )
    else:
        print(f"ERROR: No mask image file found with name: {mask_file}")
        sys.exit(-1)
    return trace


def plot_mask_assignment(
    trace, mask_file, output_file, pixel_size=0.1, label="labeled"
):
    """
    Plot assigned vs non-assigned points over mask image.

    Parameters
    ----------
    trace : ChromatinTraceTable
    mask_file : str
    output_file : str
    pixel_size : float
    label : str
    """

    # load mask
    img = load_mask(mask_file)

    # === extract coordinates ===
    x = np.array(trace.data["y"]) / pixel_size
    y = np.array(trace.data["x"]) / pixel_size

    labels = np.array(trace.data["label"])

    # points with label
    mask_labeled = np.array([label in str(lab) for lab in labels])

    x_labeled = x[mask_labeled]
    y_labeled = y[mask_labeled]

    x_not = x[~mask_labeled]
    y_not = y[~mask_labeled]

    # === plot ===
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.imshow(img, cmap="gray")

    ax.scatter(x_labeled, y_labeled, s=5, c="red", label=label, alpha=0.6)

    ax.scatter(x_not, y_not, s=5, c="cyan", label=f"not:{label}", alpha=0.6)

    ax.set_title("Mask assignment")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"$ Saved mask assignment plot: {output_file}")


def process_traces(
    trace_files=[], mask_file="", label="labeled", pixel_size=0.1, output_format="png"
):
    print(
        "\n{} trace files to process= {}".format(
            len(trace_files), "\n".join(map(str, trace_files))
        )
    )
    if trace_files:
        # iterates over traces in folder
        for trace_file in trace_files:
            trace = ChromatinTraceTable()
            trace.initialize()
            # reads new trace
            trace.load(trace_file)
            trace = assign_masks(trace, mask_file, label=label, pixel_size=pixel_size)
            base_name, _ = os.path.splitext(trace_file)

            outputfile = f"{base_name}_{label}.ecsv"
            trace.save(outputfile, comments=label)
            print(f"$ Saved output trace file at: {outputfile}")

            plot_file = f"{base_name}_{label}_mask_plot.{output_format}"
            plot_mask_assignment(
                trace,
                mask_file=mask_file,
                output_file=plot_file,
                pixel_size=pixel_size,
                label=label,
            )


def main():
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)

    print("=" * 10 + "Started execution" + "=" * 10)
    process_traces(
        trace_files=p["trace_files"],
        mask_file=p["mask_file"],
        label=p["label"],
        pixel_size=p["pixel_size"],
        output_format=p["output_format"],
    )
    print("=" * 9 + "Finished execution" + "=" * 9)


if __name__ == "__main__":
    main()
