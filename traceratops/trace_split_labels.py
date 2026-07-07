#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Split a trace file into two files based on presence of a label.

Example:
    Input:
        Trace_filtered.ecsv

    Output:
        Trace_filtered_Pdx1.ecsv
        Trace_filtered_not:Pdx1.ecsv
"""

from traceratops.script_banner import print_script_banner
import argparse
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


# === ARGUMENTS ===
def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input",
        help="Input trace file (ECSV)",
        required=False,
        default=None,
    )

    parser.add_argument(
        "--label",
        help="Label to split on (e.g. Pdx1)",
        required=True,
    )

    parser.add_argument(
        "--pipe",
        help="Read trace files from stdin",
        action="store_true",
    )

    return parser


def get_files(args):
    if args.pipe:
        return sys.stdin.read().strip().split("\n")
    else:
        return [args.input]


# === SPLIT FUNCTION ===
def split_trace_by_label(trace, label):

    data = trace.data

    mask_keep = [label in str(lab) for lab in data["label"]]

    # === NEW OBJECTS ===
    trace_keep = ChromatinTraceTable()
    trace_keep.initialize()
    trace_keep.data = data[mask_keep]

    trace_not = ChromatinTraceTable()
    trace_not.initialize()
    trace_not.data = data[[not m for m in mask_keep]]

    # === copy metadata ===
    trace_keep.data.meta = trace.data.meta.copy()
    trace_not.data.meta = trace.data.meta.copy()

    return trace_keep, trace_not


# === MAIN PROCESS ===
def process(trace_files, label):

    if len(trace_files) == 0:
        print("No input files provided")
        return

    print(f"\n$ Processing {len(trace_files)} trace file(s)")

    for trace_file in trace_files:

        print(f"\n> Processing: {trace_file}")

        trace = ChromatinTraceTable()
        trace.initialize()
        trace.load(trace_file)

        trace_keep, trace_not = split_trace_by_label(trace, label)

        base = trace_file.split(".")[0]

        out_keep = f"{base}_{label}.ecsv"
        out_not = f"{base}_not:{label}.ecsv"

        trace_keep.save(out_keep, comments=f"split:{label}")
        trace_not.save(out_not, comments=f"split:not:{label}")

        print(f"$ Saved: {out_keep} ({len(trace_keep.data)} rows)")
        print(f"$ Saved: {out_not} ({len(trace_not.data)} rows)")


# === ENTRYPOINT ===
def main():
    print_script_banner(__file__, __doc__)

    print("=" * 10 + " Started execution " + "=" * 10)

    parser = parse_arguments()
    args = parser.parse_args()

    if not args.input and not args.pipe:
        print("Error: provide --input or --pipe")
        parser.print_help()
        sys.exit(0)

    trace_files = get_files(args)

    process(trace_files, args.label)

    print("=" * 10 + " Finished execution " + "=" * 10)


if __name__ == "__main__":
    main()
