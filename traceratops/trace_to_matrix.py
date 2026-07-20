#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script converts a tracefile to a NUMPY array containing the pair-wise distances between barcodes, for each trace in the tracefile.
The dimensions of the NUMPY matrix is therefore: N x N xN_traces, where N is the number of barcodes and N_traces the number of traces.
"""

import argparse
import select
import sys

import numpy as np

from traceratops.core.build_matrix import BuildMatrix
from traceratops.script_banner import print_script_banner


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-F", "--outputFolder", help="Output folder, Default: PWD")
    parser.add_argument("--input", help="Name of input trace file.")
    parser.add_argument(
        "--distance_threshold",
        help="Threshold for the maximum distance allowed. Default: np.inf",
    )

    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of parallel workers for per-trace matrix calculation. Use -1 for all available CPUs. Default: 1",
    )
    parser.add_argument(
        "--plot_histograms",
        action="store_true",
        help="Calculate and save PWD KDE histograms. Disabled by default because this is often the slowest step.",
    )

    return parser


def create_dict_args(args):
    p = {}
    if args.outputFolder:
        p["rootFolder"] = args.outputFolder
    else:
        p["rootFolder"] = None

    if args.input:
        p["input"] = args.input
    else:
        p["input"] = None

    if args.distance_threshold:
        p["distance_threshold"] = float(args.distance_threshold)
    else:
        p["distance_threshold"] = np.inf

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

    p["n_jobs"] = args.n_jobs
    p["plot_histograms"] = args.plot_histograms

    p["colormaps"] = {
        "Nmatrix": "Blues",
        "PWD_KDE": "terrain",
        "PWD_median": "terrain",
        "contact": "coolwarm",
    }

    return p


def runtime(
    trace_files=[],
    colormaps=dict(),
    distance_threshold=np.inf,
    outputFolder=None,
    n_jobs=1,
    plot_histograms=False,
):
    if len(trace_files) < 1:
        print(
            "! Error: no trace file provided. Please either use pipe or the --input option to provide a filename."
        )
        return 0
    elif len(trace_files) == 1:
        print("\n$ trace file to process= {}".format(trace_files))
    else:
        print(
            "\n{} trace files to process= {}".format(
                len(trace_files), "\n".join(map(str, trace_files))
            )
        )

    if len(trace_files) > 0:
        # iterates over traces in folder
        for trace_file in trace_files:
            # converts trace to matrix
            param = dict()
            acq_params_dict = {
                "zBinning": 2,
                "pixelSizeXY": 0.1,
                "pixelSizeZ": 0.25,
            }
            new_matrix = BuildMatrix(param, acq_params_dict, colormaps=colormaps)
            new_matrix.launch_analysis(
                trace_file,
                distance_threshold=distance_threshold,
                outputFolder=outputFolder,
                n_jobs=n_jobs,
                plot_histograms=plot_histograms,
            )

    return len(trace_files)


def main():
    print_script_banner(__file__, __doc__)
    # [parsing arguments]
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)

    # [loops over lists of datafolders]
    n_traces_processed = runtime(
        trace_files=p["trace_files"],
        colormaps=p["colormaps"],
        distance_threshold=p["distance_threshold"],
        outputFolder=p["rootFolder"],
        n_jobs=p["n_jobs"],
        plot_histograms=p["plot_histograms"],
    )

    print(f"Processed <{n_traces_processed}> trace(s)")
    print("Finished execution")


if __name__ == "__main__":
    main()
