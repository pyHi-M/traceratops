#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates and plots matrices (PWD and proximity) from:
    - a file with single-cell PWD matrices in Numpy format
    - a file with the unique barcodes used
"""


import argparse
import itertools
import os
import sys

import numpy as np

from traceratops.core.him_matrix_operations import (
    calculate_contact_probability_matrix,
    calculate_ensemble_pwd_matrix,
    plot_him_matrix,
    plot_nan_matrix,
)


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
        "-R",
        "--remove_nan",
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


def create_output_folder(folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    print(f"Output path: {folder_path}")


def load_matrix(matrix_path):
    if not os.path.exists(matrix_path):
        raise ValueError(f"File not found: {matrix_path}")
    print(f"$ Matrix loaded: {matrix_path}")
    return np.load(matrix_path)


def load_barcodes(barcodes_path):
    unique_barcodes = list(np.loadtxt(barcodes_path, delimiter=" "))
    unique_barcodes = [int(x) for x in unique_barcodes]
    print(f"$ Unique barcodes loaded: {unique_barcodes}")
    return unique_barcodes


def shuffle_matrix(matrix, index):
    new_size = len(index)
    new_matrix = np.zeros((new_size, new_size, len(matrix[0][0])))

    if new_size > matrix.shape[0]:
        raise ValueError(
            f"Error: shuffle size {new_size} is larger than matrix dimensions {matrix.shape[0]}\nShuffle: {index}"
        )
    for i, j in itertools.product(range(new_size), range(new_size)):
        if index[i] < matrix.shape[0] and index[j] < matrix.shape[0]:
            new_matrix[i, j] = matrix[index[i], index[j]]
        else:
            raise ValueError(
                f"Out of index; matrix.shape[0]: {matrix.shape[0]} |i: {i} |index[i]: {index[i]} |j: {j} |index[j]: {index[j]}"
            )
    return new_matrix


def new_shuffle_matrix(shuffle_csl, barcode_list, sc_matrix):
    index = [barcode_list.index(int(i)) for i in shuffle_csl.split(",")]
    new_barcode_list = [barcode_list[i] for i in index]
    sc_matrix_shuffled = shuffle_matrix(sc_matrix, index)
    return new_barcode_list, sc_matrix_shuffled


def merge_matrices(mode, matrices, threshold=None, remove_nan=None):
    print(f"$ averaging method: {mode}")
    if mode == "proximity":
        print("$ calculating contact probability matrix")
        sc_matrix, nan_matrix = calculate_contact_probability_matrix(
            matrices,
            1,
            threshold=threshold,
            remove_nan=remove_nan,
        )
    else:
        nan_matrix = None
        cells_to_plot = range(matrices.shape[2])
        sc_matrix, _ = calculate_ensemble_pwd_matrix(
            matrices, 1, cells_to_plot, mode=mode
        )
    return sc_matrix, nan_matrix


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    check_required_arg(args, parser)
    create_output_folder(args.output)
    sc_matrices = load_matrix(args.matrix)
    u_barcodes = load_barcodes(args.barcodes)
    if args.shuffle:
        u_barcodes, sc_matrices = new_shuffle_matrix(
            args.shuffle, u_barcodes, sc_matrices
        )
    matrix_to_plot, nan_matrix = merge_matrices(
        args.mode, sc_matrices, args.threshold, remove_nan=args.remove_nan
    )
    cmtitle = "proximity frequency" if args.mode == "proximity" else "distance, µm"
    n_cells = sc_matrices.shape[2]
    if args.mode == "proximity":
        plot_nan_matrix(
            nan_matrix,
            u_barcodes,
            input_filename=args.matrix,
            output_folder=args.output,
            file_format=args.plot_format,
            n_cells=n_cells,
            font_size=args.fontsize,
            remove_nan=args.remove_nan,
        )

    plot_path = plot_him_matrix(
        matrix_to_plot,
        u_barcodes,
        input_filename=args.matrix,
        output_folder=args.output,
        file_format=args.plot_format,
        mode=args.mode,
        n_cells=n_cells,
        font_size=args.fontsize,
        proximity_threshold=args.threshold,
        remove_nan=args.remove_nan,
        cmtitle=cmtitle,
        c_min=args.c_min,
        c_max=args.c_max,
        c_m=args.c_map,
    )

    # saves output matrix in NPY format
    np.save(plot_path, matrix_to_plot)
    print(f"Output data: {plot_path}.npy")


if __name__ == "__main__":
    main()
