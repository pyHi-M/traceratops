#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads a chromatin trace file and computes basic statistics:
   - Number of unique ROIs
   - Number of unique chromatin traces
"""

from traceratops.script_banner import print_script_banner
import argparse
import os
import select
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", help="Path to the input trace file.")
    input_group.add_argument(
        "--pipe", action="store_true", help="Read trace filenames from stdin (pipe)."
    )
    return parser


def get_trace_files(args):
    if args.pipe:
        if select.select([sys.stdin], [], [], 0.0)[0]:
            trace_files = [line.strip() for line in sys.stdin if line.strip()]
        else:
            print(
                "Error: No filenames received from stdin. Provide input with --pipe or use --input."
            )
            sys.exit(1)
    else:
        trace_files = [args.input]

    return trace_files


def compute_trace_statistics(trace_file):
    trace_table = ChromatinTraceTable()
    trace_table.load(trace_file)

    if trace_table.data is None or len(trace_table.data) == 0:
        print("Error: The trace file is empty or could not be loaded.")
        sys.exit(1)

    # Compute statistics
    num_unique_rois = len(set(trace_table.data["ROI #"]))
    num_unique_traces = len(set(trace_table.data["Trace_ID"]))
    num_unique_barcodes = len(set(trace_table.data["Barcode #"]))

    print(f"Statistics for {trace_file}:")
    print(f"- Number of unique ROIs: {num_unique_rois}")
    print(f"- Number of unique chromatin traces: {num_unique_traces}")
    print(f"- Number of unique barcodes: {num_unique_barcodes}")


def main():
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    args = parser.parse_args()
    trace_files = get_trace_files(args)

    processed_files = 0
    for trace_file in trace_files:
        if not os.path.exists(trace_file):
            print(f"Error: The file {trace_file} does not exist.")
            continue

        compute_trace_statistics(trace_file)
        processed_files += 1

    if processed_files == 0:
        print("Error: No valid trace files to process.")
        sys.exit(1)


if __name__ == "__main__":
    main()
