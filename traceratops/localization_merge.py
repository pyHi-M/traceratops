#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A utility script for merging multiple localization files into a single consolidated output.
This script is designed to work with localization tables containing spatial coordinates
and related metadata for microscopy or imaging analysis.

The script takes localization files as input via stdin (e.g., through piping or redirection)
and merges them into a single output file. It preserves all data from the original files
while combining them into one comprehensive table.
"""

import argparse
import os
import select
import sys

from traceratops.core.localization_table import LocalizationTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output_file",
        help="Output File name. Default = merged_localizations.ecsv",
    )
    parser.add_argument("-O", "--output_folder", help="Output File name. Default = ./")
    return parser


def create_dict_args(args):
    p = {}
    if args.output_folder:
        p["outputFolder"] = args.output_folder
    else:
        p["outputFolder"] = "."
    if args.output_file:
        p["output_file"] = args.output_file
    else:
        p["output_file"] = "merged_localizations.ecsv"
    p["loc_files"] = []
    if select.select(
        [
            sys.stdin,
        ],
        [],
        [],
        0.0,
    )[0]:
        p["loc_files"] = [line.rstrip("\n") for line in sys.stdin]
    else:
        print("Nothing in stdin!\n")
    print("Input parameters\n" + "-" * 15)
    for item in p.keys():
        print("{}-->{}".format(item, p[item]))
    return p


def appends_traces(loc_files):
    new_loc_table = LocalizationTable()
    number_loc_tables = 0
    # iterates over traces in folder
    for loc_file in loc_files:
        print(f"$ loc file to process: {loc_file}")
        # reads new trace
        new_table, _ = new_loc_table.load(loc_file)
        # adds it to existing trace collection
        if number_loc_tables == 0:
            collected_tables = new_table
        else:
            collected_tables = new_loc_table.append(collected_tables, new_table)
        number_loc_tables += 1
        print(f" $ appended loc file with {len(new_table)} localizations")
    print(f" $ Merged loc file will contain {len(collected_tables)} localizations")
    return collected_tables, number_loc_tables


def load_localizations(loc_files=[]):
    if len(loc_files) > 1:
        # user provided a list of files to concatenate
        collected_tables, number_loc_tables = appends_traces(loc_files)
    print(f"Read and accumulated {number_loc_tables} localization files")
    return collected_tables, number_loc_tables


def run(p):
    print("\n" + "-" * 80)
    localizations = LocalizationTable()
    # [ creates output folder]
    if not os.path.exists(p["outputFolder"]):
        os.mkdir(p["outputFolder"])
        print("Folder created: {}".format(p["outputFolder"]))
    # loads and merges traces
    collected_tables, number_loc_tables = load_localizations(loc_files=p["loc_files"])
    # saves merged trace table
    output_file = p["output_file"]
    localizations.save(
        output_file,
        collected_tables,
        comments="appended_loc_files=" + str(number_loc_tables),
    )
    print(f"$ Saved merged file to: {output_file}")
    print("Finished execution")


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)
    print("loc_files{}".format(len(p["loc_files"])))
    if len(p["loc_files"]) < 1:
        print("\nNothing to process...\n")
    else:
        run(p)


if __name__ == "__main__":
    main()
