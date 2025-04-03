#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
**Chromatin Trace Filtering Utility**

This script processes and filters chromatin trace files based on various criteria:

- Spatial coordinates (x, y, z)
- Minimum number of barcodes per trace
- Duplicate spot removal
- Barcode-specific filtering
- Label-based filtering
- Localization intensity thresholding

The script can process single files or multiple files via pipe input.

**Outputs**

1. Filtered trace file (.ecsv) with naming convention:
    - [original_filename]_[output_tag]_[label_tag].ecsv
2. For intensity filtering:
    - Histogram plots of localization intensities (before and after filtering)
3. For duplicate barcode cleaning:
    - Statistics plots (.png) showing number of spots with same barcode per trace

**Examples**

.. code-block:: bash

    # Basic filtering with spatial constraints and minimum barcode requirement
    $ trace_filter --input Trace.ecsv --z_min 4 --z_max 5 --y_max 175 --output zy_filtered --n_barcodes 3

    # Remove duplicate spots and specific barcodes
    $ trace_filter --input Trace.ecsv --clean_spots --remove_barcode 1,3,5

    # Keep only traces with a specific label
    $ trace_filter --input Trace.ecsv --keep_label region1

    # Remove traces with a specific label
    $ trace_filter --input Trace.ecsv --remove_label region1

    # Filter by localization intensity
    $ trace_filter --input Trace.ecsv --localization_file Localizations.ecsv --intensity_min 1000

    # Process multiple files via pipe
    $ ls *Trace.ecsv | trace_filter --pipe --n_barcodes 3

**Usage**
"""

import argparse
import sys

import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable

# from traceratops.core.localization_table import LocalizationTable


def check_required_arg(args, parser):
    """
    Check that no required arguments are supplied.
    If it's the case, we make a normal code exit with the print of usage help before.
    """
    if not args.input and not args.pipe:
        print("Error: No argument provided. You must use '--input <file>' or '--pipe'.")
        print("Redirecting to `--help` option:\n")
        parser.print_help()
        sys.exit(0)


def parse_arguments():
    parser = argparse.ArgumentParser(add_help=True, description=__doc__)
    psr_basic = parser.add_argument_group(
        "Basic arguments", description="*One of --input or --pipe is required.*"
    )
    psr_basic.add_argument(
        "--input", help="Name of input trace file (ECSV format).", default=None
    )
    psr_basic.add_argument(
        "-O", "--output", help="Tag to add to the output file.", default="filtered"
    )
    psr_basic.add_argument(
        "--pipe",
        help="inputs Trace file list from stdin (for batch processing)",
        action="store_true",
    )

    psr_opt = parser.add_argument_group("Filtering options")
    psr_opt.add_argument(
        "--n_barcodes",
        help="Minimum number of barcodes by trace to keep.",
        default=2,
        type=int,
    )
    psr_opt.add_argument(
        "--clean_spots",
        help="Removes both spots with same UID and barcodes repeated in a single trace.",
        action="store_true",
    )
    psr_opt.add_argument(
        "--remove_barcode",
        help="Comma-separated list of barcode IDs to remove (e.g., ``1,2,3``)",
        default=None,
    )
    psr_opt.add_argument(
        "--remove_label",
        help="Provide a label name to remove traces with this label.",
        default=None,
    )
    psr_opt.add_argument(
        "--keep_label",
        help="Select traces containing this label, removes all other traces.",
        default=None,
    )

    psr_intensity = parser.add_argument_group("Intensity filtering")
    psr_intensity.add_argument(
        "--localization_file", default=None, help="Name of input localizations file."
    )
    psr_intensity.add_argument(
        "--intensity_min",
        type=float,
        default=0.0,
        help="Minimum intensity threshold for localizations.",
    )

    psr_coord = parser.add_argument_group(
        "Coordinate filtering", description="Coordinate limits (default: 0 to infinity)"
    )
    psr_coord.add_argument(
        "--z_min", help="Z minimum for a localization.", default=0, type=float
    )
    psr_coord.add_argument(
        "--z_max", help="Z maximum for a localization.", default=np.inf, type=float
    )
    psr_coord.add_argument(
        "--y_min", help="Y minimum for a localization.", default=0, type=float
    )
    psr_coord.add_argument(
        "--y_max", help="Y maximum for a localization.", default=np.inf, type=float
    )
    psr_coord.add_argument(
        "--x_min", help="X minimum for a localization.", default=0, type=float
    )
    psr_coord.add_argument(
        "--x_max", help="X maximum for a localization.", default=np.inf, type=float
    )

    return parser


def args_coord_to_dict(args):
    return {
        "x_min": args.x_min,
        "y_min": args.y_min,
        "z_min": args.z_min,
        "x_max": args.x_max,
        "y_max": args.y_max,
        "z_max": args.z_max,
    }


def get_files_from_args(args):
    """Read input directly without using select.select()"""
    if args.pipe:
        trace_files = sys.stdin.read().strip().split("\n")
        return trace_files
    else:
        return [args.input]


def check_file_number(trace_files):
    """checks number of trace files"""
    if len(trace_files) < 1:
        print(
            "! Error: no trace file provided. Please either use pipe or the --input option to provide a filename."
        )
        return 0
    elif len(trace_files) == 1:
        print(f"\n$ trace files to process: {trace_files}")
    else:
        f2p = "\n".join(map(str, trace_files))
        print(f"\n$ {len(trace_files)} trace files to process: \n{f2p}")


def filter_duplicat(remove_duplicate_spots, trace, trace_file):
    if remove_duplicate_spots:
        # remove duplicated UID spots
        trace.remove_duplicates()
        # removes barcodes in traces where they are repeated
        trace.filter_repeated_barcodes(trace_file)
    return trace


def filter_barcode_number(n_barcodes, trace, comments):
    # filters trace by minimum number of barcodes
    if n_barcodes > 1:
        trace.filter_traces_by_n(minimum_number_barcodes=n_barcodes)
        comments.append("filt:n_barcodes>" + str(n_barcodes))
    return trace, comments


def filter_label(label_to_keep, label_to_remove, trace):
    if label_to_keep is not None:
        trace.trace_keep_label(label_to_keep)
        file_tag = "_" + label_to_keep
    elif label_to_remove is not None:
        trace.trace_remove_label(label_to_remove)
        file_tag = "_not-" + label_to_remove
    else:
        file_tag = ""
    return trace, file_tag


def runtime(
    trace_files=[],
    n_barcodes=2,
    coord_limits=dict(),
    tag="filtered",
    remove_duplicate_spots=False,
    remove_barcode=None,
    label_to_keep="",
    label_to_remove="",
):
    if len(trace_files) <= 0:
        print("No trace file found to process!")
        return len(trace_files)

    # iterates over traces
    for trace_file in trace_files:
        trace = ChromatinTraceTable()
        trace.initialize()
        comments = list()
        # reads new trace
        trace.load(trace_file)

        trace = filter_duplicat(remove_duplicate_spots, trace, trace_file)
        trace, comments = filter_barcode_number(n_barcodes, trace, comments)

        # filters trace by coordinate
        for coord in ["x", "y", "z"]:
            coor_min = coord_limits[coord + "_min"]
            coor_max = coord_limits[coord + "_max"]

            if coor_min > 0.0 or coor_max != np.inf:
                trace.filter_traces_by_coordinate(
                    coor=coord,
                    coor_min=coor_min,
                    coor_max=coor_max,
                )
                comments.append("filt:{}<{}>{}".format(coor_min, coord, coor_max))

        if remove_barcode is not None:
            bc_list = remove_barcode.split(",")
            print(f"\n$ Removing barcodes: {bc_list}")
            for bc in bc_list:
                trace.remove_barcode(bc)

        trace, file_tag = filter_label(label_to_keep, label_to_remove, trace)

        # saves output trace
        outputfile = trace_file.split(".")[0] + "_" + tag + file_tag + ".ecsv"
        trace.save(outputfile, comments=", ".join(comments))
        print(f"$ Saved output trace file at: {outputfile}")
    return len(trace_files)


def main():
    print("=" * 10 + "Started execution" + "=" * 10)
    # [parsing arguments]
    parser = parse_arguments()
    args = parser.parse_args()
    check_required_arg(args, parser)

    trace_files = get_files_from_args(args)
    check_file_number(trace_files)

    # [loops over lists of datafolders]
    n_traces_processed = runtime(
        trace_files=trace_files,
        n_barcodes=args.n_barcodes,
        coord_limits=args_coord_to_dict(args),
        tag=args.output,
        remove_duplicate_spots=args.clean_spots,
        remove_barcode=args.remove_barcode,
        label_to_keep=args.keep_label,
        label_to_remove=args.remove_label,
        localizations_file=args.localization_file,
        intensity_min=args.intensity_min,
    )

    print(f"Processed <{n_traces_processed}> trace file(s)\n")
    print("=" * 9 + "Finished execution" + "=" * 9)


if __name__ == "__main__":
    main()
