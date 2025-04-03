#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This just takes a list of trace files and merges them together

$ ls Trace*.ecsv | trace_merge.py

outputs

ChromatinTraceTable() object and output .ecsv formatted file with assembled trace tables.
"""

import argparse
import os
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Simpler version of trace_combinator.
        This just takes a list of trace files and merges them together.

        Usage example:

        ``ls trace*.ecsv | trace_merge``

        or

        ``trace_merge.py --traces <folder_path_with_trace_files>``
        """
    )
    parser.add_argument(
        "-T", "--traces", help="Input folder with traces to merge", default=None
    )
    parser.add_argument(
        "--pipe",
        help="input trace file list from stdin (pipe)",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "-N", "--name", help="Output file name", default="merged_traces.ecsv"
    )
    parser.add_argument("-F", "--folder", help="Output folder", default=os.getcwd())
    return parser


def check_required_arg(args, parser):
    """
    Check that no required arguments are supplied.
    If it's the case, we make a normal code exit with the print of usage help before.
    """
    if not args.traces and not args.pipe:
        print("Error: No argument provided. You must use '--traces' or '--pipe'")
        print("Redirecting to `--help` option:\n")
        parser.print_help()
        sys.exit(0)


def get_files_from_args(args):
    """Read input directly without using select.select()"""
    if args.pipe:
        trace_files = sys.stdin.read().strip().split("\n")
        return trace_files
    else:
        return os.listdir(args.traces)


def appends_traces(traces, trace_files):
    new_trace = ChromatinTraceTable()
    # iterates over traces in folder
    for trace_file in trace_files:
        # reads new trace
        new_trace.load(trace_file)
        # adds it to existing trace collection

        traces.append(new_trace.data)
        traces.number_traces += 1
        print(f" $ appended trace file with {len(new_trace.data)} traces")
    print(f" $ Merged trace file will contain {len(traces.data)} traces")
    return traces


def load_traces(trace_files=[]):
    traces = ChromatinTraceTable()
    traces.initialize()
    if len(trace_files) > 1:
        # user provided a list of files to concatenate
        traces = appends_traces(traces, trace_files)
    print(f"Read and accumulated {traces.number_traces} trace files")
    return traces


def create_out_folder(folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
        print(f"Folder created: {folder_path}")


def main():
    # [parsing arguments]
    parser = parse_arguments()
    args = parser.parse_args()
    check_required_arg(args, parser)
    trace_files = get_files_from_args(args)

    print(f"Number of trace files to merge: {len(trace_files)}")
    if len(trace_files) < 2:
        raise ValueError("\nNothing to process...\n")

    create_out_folder(args.folder)
    traces = load_traces(trace_files)
    traces.save(
        args.name,
        comments="appended_trace_files=" + str(traces.number_traces),
    )

    print("Finished execution")


if __name__ == "__main__":
    main()
