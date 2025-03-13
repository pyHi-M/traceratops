#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates and plots matrices (PWD and proximity) from:
    - a file with single-cell PWD matrices in Numpy format
    - a file with the unique barcodes used


Example:
$ figureHiMmatrix.py -T Trace_3D_barcode_KDtree_ROI:4_Matrix_PWDsc_matrix.npy -U Trace_3D_barcode_KDtree_ROI:4_Matrix_uniqueBarcodes.ecsv

Options:
    - plottingFileExtension: format of figure
    - cScale: value of the max of the cScale used to plot the matrix
    - cmap: name of cmap
    - mode: indicated the plotting mode, either ["proximity"] or ["KDE", "median"] for PWD matrix.
    - outputFolder: name of outputfolder. 'plots' is the default

Left to do:
    - need to implement a way to select a subset of chromatin traces...
"""


import argparse
import os
import sys

import numpy as np

from traceratops.core.him_matrix_operations import plot_matrix
from traceratops.core.plotting_functions import gets_matrix


def parse_arguments():
    # [parsing arguments]
    parser = argparse.ArgumentParser(
        add_help=True,
        description="""calculates and plots matrices (PWD and proximity) from:
    - a file with single-cell PWD matrices in Numpy format
    - a file with the unique barcodes used""",
    )
    parser_required = parser.add_argument_group(
        "Required arguments", description="*These both arguments are required.*"
    )
    parser_required.add_argument(
        "-M", "--matrix", help="Filename of single-cell PWD matrices in NPY format"
    )
    parser_required.add_argument(
        "-B", "--barcodes", help="csv file with a simple list of unique barcodes (int)"
    )

    parser_advanced = parser.add_argument_group(
        "Advanced arguments",
        description="[Optional] Advanced args to personalize outputs",
    )
    parser_advanced.add_argument("-O", "--output", help="Folder for outputs")
    parser_advanced.add_argument(
        "--plot_format", help="By default: png. Other options: svg, pdf, png"
    )
    parser_advanced.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
    )
    parser_advanced.add_argument("-T", "--threshold", help="Proximity threshold in µm")
    parser_advanced.add_argument(
        "--mode",
        help="Mode used to calculate the mean distance. Can be either 'median', 'KDE' or 'proximity'. Default: proximity",
    )
    parser_advanced.add_argument(
        "--norm",
        help="Matrix normalization mode. If activate, remove NaN values before compute statistics on bin.",
        action="store_true",
    )
    parser_visu = parser.add_argument_group(
        "Visualization arguments", description="[Optional] Custom visualization"
    )
    parser_visu.add_argument("--c_min", help="Colormap min scale. Default: 0")
    parser_visu.add_argument("--c_max", help="Colormap max scale. Default: automatic")
    parser_visu.add_argument(
        "--c_map",
        help="Colormap (see: matplotlib > colormaps > diverging). Default: coolwarm",
    )
    parser_visu.add_argument("--fontsize", help="Size of fonts to be used in matrix")
    return parser


def check_required_arg(args, parser):
    """
    Check that no required arguments are supplied.
    If it's the case, we make a normal code exit with the print of usage help before.
    """
    if not args.barcodes or not args.matrix:
        print(
            "Error: No argument provided. You must use '--matrix <file>' and '--barcodes'."
        )
        print("Redirecting to `--help` option:\n")
        parser.print_help()
        sys.exit(0)


def create_dict_args(args):
    run_parameters = {}

    if args.matrix:
        run_parameters["scPWDMatrix_filename"] = args.matrix
    else:
        print(
            ">> ERROR: you must provide a filename with the single cell PWD matrices in Numpy format"
        )
        sys.exit(-1)

    if args.barcodes:
        run_parameters["uniqueBarcodes"] = args.barcodes
    else:
        print(">> ERROR: you must provide a CSV file with the unique barcodes used")
        sys.exit(-1)

    if args.output:
        run_parameters["outputFolder"] = args.output
    else:
        run_parameters["outputFolder"] = "plots"

    if args.fontsize:
        run_parameters["fontsize"] = args.fontsize
    else:
        run_parameters["fontsize"] = 22

    if args.threshold:
        run_parameters["proximity_threshold"] = float(args.threshold)
    else:
        run_parameters["proximity_threshold"] = 0.25

    if args.c_max:
        run_parameters["cMax"] = float(args.c_max)
    else:
        run_parameters["cMax"] = 0.0

    if args.c_min:
        run_parameters["cMin"] = float(args.c_min)
    else:
        run_parameters["cMin"] = 0.0

    if args.plot_format:
        run_parameters["plottingFileExtension"] = "." + args.plot_format
    else:
        run_parameters["plottingFileExtension"] = ".png"

    if args.shuffle:
        run_parameters["shuffle"] = args.shuffle
    else:
        run_parameters["shuffle"] = 0

    run_parameters["pixelSize"] = 1

    if args.c_map:
        run_parameters["cmap"] = args.c_map
    else:
        run_parameters["cmap"] = "coolwarm"

    if args.mode:
        run_parameters["dist_calc_mode"] = args.mode
    else:
        run_parameters["dist_calc_mode"] = "proximity"

    if run_parameters["dist_calc_mode"] == "proximity":
        run_parameters["cmtitle"] = "proximity frequency"
    else:
        run_parameters["cmtitle"] = "distance, µm"

    if args.norm:
        run_parameters["matrix_norm_mode"] = "nonNANs"
    else:
        run_parameters["matrix_norm_mode"] = "n_cells"

    return run_parameters


# %%

# =============================================================================
# MAIN
# =============================================================================


def main():
    print(">>> Producing HiM matrix")

    parser = parse_arguments()
    args = parser.parse_args()
    check_required_arg(args, parser)
    run_parameters = create_dict_args(args)

    if not os.path.exists(run_parameters["outputFolder"]):
        os.mkdir(run_parameters["outputFolder"])
        print("Folder created: {}".format(run_parameters["outputFolder"]))

    (
        sc_matrix,
        uniqueBarcodes,
        n_cells,
        outputFileName,
        nan_matrix,
    ) = gets_matrix(
        run_parameters,
        scPWDMatrix_filename=run_parameters["scPWDMatrix_filename"],
        uniqueBarcodes=run_parameters["uniqueBarcodes"],
    )

    meansc_matrix, fileNameEnding = plot_matrix(
        sc_matrix,
        uniqueBarcodes,
        run_parameters["pixelSize"],
        1,
        outputFileName,
        "log",
        figtitle="Map: " + run_parameters["dist_calc_mode"],
        mode=run_parameters["dist_calc_mode"],  # median or KDE
        clim=run_parameters["cMax"],
        c_min=run_parameters["cMin"],
        n_cells=n_cells,
        c_m=run_parameters["cmap"],
        cmtitle=run_parameters["cmtitle"],
        filename_addon="_"
        + run_parameters["dist_calc_mode"]
        + "_"
        + run_parameters["matrix_norm_mode"],
        filename_extension=run_parameters["plottingFileExtension"],
        font_size=run_parameters["fontsize"],
        proximity_threshold=run_parameters["proximity_threshold"],
        nan_matrix=nan_matrix,
    )

    print("Output figure: {}".format(outputFileName))

    # saves output matrix in NPY format
    outputFileName = outputFileName + fileNameEnding
    np.save(outputFileName, meansc_matrix)
    print("Output data: {}.npy".format(outputFileName))

    print("\nDone\n\n")


if __name__ == "__main__":
    main()
