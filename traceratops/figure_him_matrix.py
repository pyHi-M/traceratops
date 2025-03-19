#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates and plots matrices (PWD and proximity) from:
    - a file with single-cell PWD matrices in Numpy format
    - a file with the unique barcodes used
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
    - a file with the unique barcodes used

Outputs:
    - Matrix (NxN) of mean pairwise distance, stored in a NPY file.
    - Visualization of this matrix in PNG format (or PDF/SVG)
    - Shadow matrix of NaN (Numpy null value) percentage for each pairwise distance
    """,
    )
    parser_required = parser.add_argument_group(
        "Required arguments", description="*These both arguments are required.*"
    )
    parser_required.add_argument(
        "-M",
        "--matrix",
        help="Filename of single-cell PWD matrices in NPY format",
        default=None,
    )
    parser_required.add_argument(
        "-B",
        "--barcodes",
        help="csv file with a simple list of unique barcodes (int)",
        default=None,
    )

    parser_advanced = parser.add_argument_group(
        "Advanced arguments",
        description="[Optional] Advanced args to personalize outputs",
    )
    parser_advanced.add_argument(
        "-O", "--output", help="Folder for outputs", default="plots"
    )
    parser_advanced.add_argument(
        "--plot_format",
        help="Available options: svg, pdf, png",
        default="png",
    )
    parser_advanced.add_argument(
        "--shuffle",
        help="Provide shuffle vector: 0,1,2,3... of the same size or smaller than the original matrix. No spaces! comma-separated!",
        default=None,
    )
    parser_advanced.add_argument(
        "--mode",
        help="Mode used to calculate the mean distance. Can be either 'median', 'KDE' or 'proximity'",
        default="proximity",
    )

    parser_proximity = parser.add_argument_group(
        "Proximity arguments",
        description="[Optional] Only for proximity ```--mode```",
    )
    parser_proximity.add_argument(
        "-T", "--threshold", help="Proximity threshold in µm", default=0.25, type=float
    )
    parser_proximity.add_argument(
        "--norm",
        help="Matrix normalization mode. If activate, remove NaN values before compute statistics on bin.",
        action="store_true",
    )

    parser_visu = parser.add_argument_group(
        "Visualization arguments", description="[Optional] Custom visualization"
    )
    parser_visu.add_argument(
        "--c_min",
        help="Colormap min scale: automatic mode by default (detects the first frequency higher than 0)",
        default=-1.0,
        type=float,
    )
    parser_visu.add_argument(
        "--c_max",
        help="Colormap max scale: automatic mode by default (detects the highest frequency value)",
        default=0.0,
        type=float,
    )
    parser_visu.add_argument(
        "--c_map",
        help="Colormap (see: matplotlib > colormaps > diverging)",
        default="coolwarm",
    )
    parser_visu.add_argument(
        "--fontsize", help="Size of fonts to be used in matrix", default=22
    )
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


# =============================================================================
# MAIN
# =============================================================================


def main():
    print(">>> Producing HiM matrix")

    parser = parse_arguments()
    args = parser.parse_args()
    check_required_arg(args, parser)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
        print(f"Folder created: {args.output}")
    matrix_norm_mode = "nonNANs" if args.norm else "n_cells"
    (
        sc_matrix,
        uniqueBarcodes,
        n_cells,
        outputFileName,
        nan_matrix,
    ) = gets_matrix(
        scPWDMatrix_filename=args.matrix,
        uniqueBarcodes=args.barcodes,
        out_folder=args.output,
        proximity_threshold=args.threshold,
        shuffle=args.shuffle,
        dist_calc_mode=args.mode,
        matrix_norm_mode=matrix_norm_mode,
    )

    base_filename = "_" + args.mode
    base_filename += "_norm" if args.norm else ""
    cmtitle = "proximity frequency" if args.mode == "proximity" else "distance, µm"

    meansc_matrix, fileNameEnding = plot_matrix(
        sc_matrix,
        uniqueBarcodes,
        1,
        1,
        outputFileName,
        "log",
        figtitle="Mode: " + args.mode,
        mode=args.mode,
        clim=args.c_max,
        c_min=args.c_min,
        n_cells=n_cells,
        c_m=args.c_map,
        cmtitle=cmtitle,
        filename_addon=base_filename,
        filename_extension="." + args.plot_format,
        font_size=args.fontsize,
        proximity_threshold=args.threshold,
        nan_matrix=nan_matrix,
        matrix_norm_mode=matrix_norm_mode,
    )

    print("Output figure: {}".format(outputFileName))

    # saves output matrix in NPY format
    outputFileName = outputFileName + fileNameEnding
    np.save(outputFileName, meansc_matrix)
    print("Output data: {}.npy".format(outputFileName))

    print("\nDone\n\n")


if __name__ == "__main__":
    main()
