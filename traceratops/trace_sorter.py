#!/usr/bin/env python3
"""Sort complete chromatin traces by their size or radius of gyration.

Rows belonging to a trace remain together and retain their original order.
"""

import argparse
import os
import select
import sys

import numpy as np

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.script_banner import print_script_banner
from traceratops.trace_splitter import compute_radius_of_gyration


def parse_arguments():
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", help="Path to the input trace file.")
    input_group.add_argument(
        "--pipe", action="store_true", help="Read trace filenames from stdin (pipe)."
    )
    parser.add_argument("-O", "--output", help="Path to the sorted output trace file.")
    parser.add_argument(
        "--sort_by",
        required=True,
        choices=("number_of_spots", "radius_of_gyration"),
        help="Trace-level quantity on which to sort.",
    )
    parser.add_argument(
        "--ascending",
        action="store_true",
        help="Sort from smallest to largest (default: descending).",
    )
    return parser


def get_trace_files(args):
    """Return filenames supplied directly or over standard input."""
    if not args.pipe:
        return [args.input]
    if select.select([sys.stdin], [], [], 0.0)[0]:
        return [line.strip() for line in sys.stdin if line.strip()]
    return []


def sorted_output_filename(input_filename):
    """Insert ``_sorted`` before an input filename's extension."""
    basename, extension = os.path.splitext(input_filename)
    return f"{basename}_sorted{extension}"


def sort_trace_table(trace_table, sort_by, ascending=False):
    """Sort a table by a trace-level property, modifying it in place."""
    data = trace_table.data
    trace_ids = list(dict.fromkeys(data["Trace_ID"]))
    if not trace_ids:
        return trace_table

    indices_by_id = {
        trace_id: np.flatnonzero(data["Trace_ID"] == trace_id) for trace_id in trace_ids
    }

    if sort_by == "number_of_spots":
        values = {trace_id: len(indices_by_id[trace_id]) for trace_id in trace_ids}
    elif sort_by == "radius_of_gyration":
        values = {}
        for trace_id in trace_ids:
            trace = data[indices_by_id[trace_id]]
            coordinates = np.vstack((trace["x"], trace["y"], trace["z"])).T
            values[trace_id] = compute_radius_of_gyration(coordinates)
    else:
        raise ValueError(f"Unknown sorting criterion: {sort_by}")

    ordered_ids = sorted(trace_ids, key=values.__getitem__, reverse=not ascending)
    ordered_indices = np.concatenate(
        [indices_by_id[trace_id] for trace_id in ordered_ids]
    )
    trace_table.data = data[ordered_indices]
    return trace_table


def main():
    """Load, sort, and save each requested chromatin trace table."""
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    args = parser.parse_args()
    trace_files = get_trace_files(args)
    if not trace_files:
        parser.error("No filenames received from stdin.")
    if args.output and len(trace_files) > 1:
        parser.error("--output can only be used when processing one input file.")

    for trace_file in trace_files:
        trace_table = ChromatinTraceTable()
        trace_table.load(trace_file)
        sort_trace_table(trace_table, args.sort_by, args.ascending)
        output = args.output or sorted_output_filename(trace_file)
        trace_table.save(output)


if __name__ == "__main__":
    main()
