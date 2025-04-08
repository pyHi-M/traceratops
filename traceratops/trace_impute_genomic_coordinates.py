#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assign genomic coordinates to a chromatin trace table.
It assigns genomic coordinates (Chrom, Chrom_Start, Chrom_End) from the BED file to each row
in the trace table based on the 'Barcode #' column.
"""

import argparse
import select
import sys

from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Path to the input trace file (ECSV format).")
    parser.add_argument(
        "--bed",
        required=True,
        help="Path to the BED file containing genomic coordinates.",
    )
    parser.add_argument("--output", help="Path to save the updated trace file.")
    parser.add_argument(
        "--pipe", help="inputs Trace file list from stdin (pipe)", action="store_true"
    )
    parser.add_argument(
        "--auto-continue",
        action="store_true",
        help="Automatically continue processing even with unmatched barcodes",
    )

    return parser


def create_dict_args(args):
    p = dict()
    p["trace_files"] = []
    p["auto_continue"] = args.auto_continue

    if args.pipe:
        p["pipe"] = True
        if select.select([sys.stdin], [], [], 0.0)[0]:
            p["trace_files"] = [line.rstrip("\n") for line in sys.stdin]
        else:
            print("Nothing in stdin")
    else:
        p["pipe"] = False
        if args.input:
            p["trace_files"] = [args.input]
        else:
            print("Error: No input file specified. Use --input or --pipe.")
            sys.exit(1)

    return p


def load_bed_file(bed_file):
    """Loads the BED file into a dictionary mapping barcode numbers to genomic coordinates.
    Handles files with inconsistent tab spacing by using regex splitting."""
    bed_dict = {}

    with open(bed_file, "r") as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue

            # Split on any number of whitespace characters
            # This handles inconsistent tabs/spaces more robustly
            fields = line.strip().split()

            # Ensure we have exactly 4 fields
            if len(fields) != 4:
                print(f"Warning: Skipping malformed line: {line.strip()}")
                continue

            try:
                chrom = fields[0]
                chrom_start = int(fields[1])
                chrom_end = int(fields[2])
                barcode = int(fields[3])

                bed_dict[barcode] = {
                    "Chrom": chrom,
                    "Chrom_Start": chrom_start,
                    "Chrom_End": chrom_end,
                }
            except ValueError as e:
                print(f"Warning: Skipping line with invalid data types: {line.strip()}")
                print(f"Error: {e}")

    if not bed_dict:
        raise ValueError("No valid entries found in the BED file.")

    print(f"Successfully loaded {len(bed_dict)} barcode mappings from BED file.")
    return bed_dict


def impute_genomic_coordinates(trace_file, bed_dict, output_file, p):
    """Updates the Chrom, Chrom_Start, and Chrom_End columns in the trace file based on the BED file."""
    trace_table = ChromatinTraceTable()
    trace_table.load(trace_file)

    if trace_table.data is None or len(trace_table.data) == 0:
        print("Error: The trace file is empty or could not be loaded.")
        return

    unmatched_barcodes = set()
    matched_count = 0
    total_count = 0

    for row in trace_table.data:
        barcode = row["Barcode #"]
        total_count += 1

        if barcode in bed_dict:
            row["Chrom"] = bed_dict[barcode]["Chrom"]
            row["Chrom_Start"] = bed_dict[barcode]["Chrom_Start"]
            row["Chrom_End"] = bed_dict[barcode]["Chrom_End"]
            matched_count += 1
        else:
            unmatched_barcodes.add(barcode)

    if unmatched_barcodes:
        missing_percent = (len(unmatched_barcodes) / total_count) * 100
        print(
            f"Warning: {len(unmatched_barcodes)} unique barcodes ({missing_percent:.1f}% of rows) couldn't be matched:"
        )

        # Only show up to 10 unmatched barcodes to avoid cluttering the output
        if len(unmatched_barcodes) <= 10:
            print(
                f"  Unmatched barcodes: {', '.join(map(str, sorted(unmatched_barcodes)))}"
            )
        else:
            print(
                f"  First 10 unmatched barcodes: {', '.join(map(str, sorted(list(unmatched_barcodes))[:10]))}"
            )
            print(f"  ... and {len(unmatched_barcodes) - 10} more")

        # Ask the user if they want to continue if more than 10% of barcodes are unmatched
        if missing_percent > 10 and not p.get("auto_continue", False):
            response = input(
                "More than 10% of barcodes couldn't be matched. Continue anyway? (y/n): "
            )
            if response.lower() != "y":
                print("Operation aborted by user.")
                return

    print(
        f"Successfully matched {matched_count}/{total_count} rows ({(matched_count/total_count)*100:.1f}%)"
    )
    trace_table.save(
        output_file,
        trace_table.data,
        comments=f"Genomic coordinates imputed from BED file. {matched_count}/{total_count} rows matched.",
    )
    print(f"Updated trace file saved to {output_file}")


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    p = create_dict_args(args)
    bed_dict = load_bed_file(args.bed)

    trace_files = p["trace_files"]
    if len(trace_files) > 0:
        print(
            f"\n{len(trace_files)} trace files to process= {', '.join(map(str, trace_files))}"
        )

        # Iterate over all trace files
        for trace_file in trace_files:
            output_file = (
                args.output
                if args.output
                else trace_file.replace(".ecsv", "_imputed.ecsv")
            )
            print(f"\nProcessing file: {trace_file}")
            impute_genomic_coordinates(trace_file, bed_dict, output_file, p)
            print(f"Completed: {output_file}")
    else:
        print(
            "! Error: did not find any trace file to analyze. Please provide one using --input or --pipe."
        )


if __name__ == "__main__":
    main()
