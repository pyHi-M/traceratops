#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load a localization table and plot localization quality-control distributions.

The script produces a figure with:
- the number of localizations per barcode
- the snr distribution per barcode
- scatterplot of the snr versus z
- scatterplot of roundness versus skew
"""

import argparse
import os

from traceratops.core.localization_table import LocalizationTable
from traceratops.script_banner import print_script_banner


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-L",
        "--localization_file",
        required=True,
        help="Localization file path.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help=(
            "Output plot file name. Defaults to the name used by "
            "LocalizationTable.plot_distribution_fluxes."
        ),
    )
    parser.add_argument(
        "--output_format",
        "-f",
        "--format",
        dest="output_format",
        choices=["png", "svg", "pdf"],
        default="png",
        help="Output plot format. Default = png.",
    )
    return parser


def create_dict_args(args):
    p = {}
    p["localization_file"] = args.localization_file
    p["output_file"] = args.output_file
    p["format"] = args.output_format

    print("Input parameters\n" + "-" * 15)
    for item in p.keys():
        print("{}-->{}".format(item, p[item]))

    return p


def get_output_file(output_file, output_format):
    if output_file is None:
        output_file = "localization_distribution_fluxes"

    output_root, _ = os.path.splitext(output_file)
    return f"{output_root}.{output_format}"


def run(p):
    localization_table = LocalizationTable()
    barcode_map, _ = localization_table.load(p["localization_file"])
    output_file = get_output_file(p["output_file"], p["format"])

    localization_table.plot_distribution_fluxes(barcode_map, [output_file])

    print("Finished execution")


def main():
    print_script_banner(__file__, __doc__)
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)
    run(p)


if __name__ == "__main__":
    main()
