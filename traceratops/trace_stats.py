#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads a chromatin trace file and computes basic statistics:
- Number of unique ROIs
- Number of unique chromatin traces

Usage:
    python trace_stats.py --input trace_file.ecsv
"""

import argparse
import os
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to the input trace file.")
    return parser


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
    parser = parse_arguments()
    args = parser.parse_args()
    if not os.path.exists(args.input):
        print(f"Error: The file {args.input} does not exist.")
        sys.exit(1)

    compute_trace_statistics(args.input)


if __name__ == "__main__":
    main()
