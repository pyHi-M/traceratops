#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to plot one or multiple traces in 3D

Takes a trace file and either:
    - exports the first N traces to PDB with --number_traces
    - exports one trace to PDB with --selected_trace
    - exports all traces to PDB with --all
    - saves output coordinates for selected traces in PDB format so they can be loaded by other means including https://www.rcsb.org/3d-view, pymol, or nglviewer.
"""

import argparse
import os
import select
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.core.him_matrix_operations import write_xyz_2_pdb
from traceratops.core.io_manager import create_folder, load_barcode_dict
from traceratops.script_banner import print_script_banner


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument(
        "-n",
        "--number_traces",
        help="Number of traces treated. Ignored when --all or --selected_trace is used.",
        type=int,
    )
    parser.add_argument(
        "-N", "--N_barcodes", help="minimum_number_barcodes. Default = 2"
    )
    parser.add_argument("--selected_trace", help="Selected trace for analysis")
    parser.add_argument(
        "--barcode_type_dict",
        help="Json dictionary linking barcodes and atom types (MUST BE 3 characters long!). ",
    )
    parser.add_argument(
        "--all", help="plots all traces in trace file", action="store_true"
    )
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument(
        "-O", "--output", help="Tag to add to the output file. Default = filtered"
    )

    return parser


def create_dict_args(args):
    p = {}
    if args.input:
        p["input"] = args.input
    else:
        p["input"] = None

    if args.output:
        p["output"] = args.output
    else:
        p["output"] = "PDBs"

    if args.N_barcodes:
        p["N_barcodes"] = int(args.N_barcodes)
    else:
        p["N_barcodes"] = 2

    if args.number_traces is not None:
        if args.number_traces < 1:
            raise ValueError("--number_traces must be at least 1")
        p["number_traces"] = args.number_traces
    else:
        p["number_traces"] = 10

    if args.selected_trace:
        p["selected_trace"] = args.selected_trace
    else:
        p["selected_trace"] = "all"

    if args.barcode_type_dict:
        p["barcode_type_dict"] = args.barcode_type_dict
    else:
        p["barcode_type_dict"] = "barcode_type_dict.json"

    if args.all:
        p["select_traces"] = "all"
    elif args.selected_trace:
        p["select_traces"] = "selected"
    else:
        p["select_traces"] = "first"

    p["trace_files"] = []
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
        p["trace_files"] = [p["input"]]

    return p


def runtime(
    N_barcodes=2,
    trace_files=None,
    selected_trace="fa9f0eb5-abcc-4730-bcc7-ba1da682d776",
    barcode_type=None,
    folder_path="./PDBs",
    select_traces="one",
    number_traces=2,
):
    # gets trace files
    trace_files = [] if trace_files is None else trace_files
    barcode_type = {} if barcode_type is None else barcode_type

    print(
        f"$ Will retrieve {number_traces} traces using the criterium: {select_traces}"
    )
    if len(trace_files) > 0:
        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

        # iterates over traces in folder
        for trace_file in trace_files:
            trace = ChromatinTraceTable()
            trace.initialize()

            # reads new trace
            trace.load(trace_file)

            # filters trace
            # trace.filter_traces_by_n(minimum_number_barcodes=N_barcodes)

            # indexes traces by Trace_ID
            trace_table = trace.data
            trace_table_indexed = trace_table.group_by("Trace_ID")
            total_traces = len(trace_table_indexed.groups)
            if select_traces == "first":
                traces_to_process = min(number_traces, total_traces)
            elif select_traces == "all":
                traces_to_process = total_traces
            else:
                traces_to_process = 1

            print(f"$ number of traces available: {total_traces}")
            print(f"$ number of traces to process: {traces_to_process}")

            # iterates over traces
            for idx, single_trace in enumerate(trace_table_indexed.groups):
                trace_id = single_trace["Trace_ID"][0]
                flag = False

                if select_traces == "selected" and trace_id == selected_trace:
                    flag = True
                elif select_traces == "all":
                    flag = True
                elif select_traces == "first" and idx < number_traces:
                    flag = True

                if flag:
                    print("Converting trace ID: {}".format(trace_id))

                    # sorts by barcode
                    new_trace = single_trace.copy()
                    new_trace = new_trace.group_by("Barcode #")
                    # ascii.write(new_trace['Barcode #', 'x','y','z'], selected_trace+'.ecsv', overwrite=True)

                    write_xyz_2_pdb(
                        folder_path + os.sep + trace_id + ".pdb",
                        new_trace,
                        barcode_type,
                    )
    else:
        print("No trace file found to process!")

    return len(trace_files)


def main():
    print_script_banner(__file__, __doc__)

    # [parsing arguments]
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)
    # [loops over lists of datafolders]
    barcode_type = load_barcode_dict(p["barcode_type_dict"])

    # creates output folder
    output_folder = p["output"]
    folder_path = os.path.join(os.getcwd(), output_folder)

    create_folder(folder_path)

    n_traces_processed = runtime(
        N_barcodes=p["N_barcodes"],
        trace_files=p["trace_files"],
        selected_trace=p["selected_trace"],
        barcode_type=barcode_type,
        folder_path=folder_path,
        select_traces=p["select_traces"],
        number_traces=p["number_traces"],
    )

    print(f"Processed <{n_traces_processed}> trace file(s)")
    print("Finished execution")


if __name__ == "__main__":
    main()
