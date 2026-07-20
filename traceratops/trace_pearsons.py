#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare multiple chromatin trace tables by computing pairwise distances
between barcode combinations and quantifying similarity via Pearson correlation.These can be different fields of view of a single acquisition or different biological replicates.

This tool is useful for analyzing the structural similarity between different chromatin
trace datasets, helping to identify patterns and relationships in chromatin organization
across multiple samples or conditions.

The script calculates, for each tracefile, the median pairwise distance map.
Then it calculates the bin-by-bin Pearson correlation between all tracefiles provided.
The result is therefore a matrix of Pearson correlations.


"""

import argparse
import itertools
import os
import select
import sys
from collections import defaultdict
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.script_banner import print_script_banner


def silent_load_trace(trace, path):
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        trace.load(path)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        nargs="+",
        help="Paths to trace files to compare.",
    )
    input_group.add_argument(
        "--pipe", action="store_true", help="Read input filenames from stdin (pipe)."
    )
    parser.add_argument(
        "-O",
        "--output",
        default=".",
        help="Output folder name for the correlation matrix plot (filename: trace_correlation_matrix.png)",
    )
    parser.add_argument(
        "--vmin", type=float, default=-10, help="Minimum value for colormap scaling"
    )
    parser.add_argument(
        "--vmax", type=float, default=10, help="Maximum value for colormap scaling"
    )
    parser.add_argument(
        "--output_format",
        choices=["png", "svg", "pdf"],
        default="png",
        help="Output image format. Default = png.",
    )
    parser.add_argument("--verbose", action="store_true", help="Increase verbosity")
    return parser


def get_trace_files(args):
    if args.pipe:
        if select.select([sys.stdin], [], [], 0.0)[0]:
            trace_files = [line.strip() for line in sys.stdin if line.strip()]
        else:
            print(
                "Error: No filenames received from stdin. Provide input with --pipe or use --input."
            )
            sys.exit(1)
    else:
        trace_files = args.input

    return trace_files


def find_unique_substrings(filenames):
    """
    Find unique substrings in each filename that distinguish them from all others.

    This function identifies the minimal substring that makes each filename unique
    within the given set, without relying on predefined patterns.

    Parameters:
    ----------
    filenames : list
        List of filenames to analyze

    Returns:
    -------
    dict
        Dictionary mapping each filename to its unique identifying substring
    """

    # Create a list of tuples with each filename and its character-by-character comparison
    file_chars = []
    for filename in filenames:
        file_chars.append((filename, list(filename)))

    # Dictionary to store unique identifiers
    unique_identifiers = {}

    for i, (filename, chars) in enumerate(file_chars):
        # Compare with all other filenames
        differing_positions = []

        for j, (other_filename, other_chars) in enumerate(file_chars):
            if i == j:  # Skip self-comparison
                continue

            # Find differing positions
            min_len = min(len(chars), len(other_chars))
            for pos in range(min_len):
                if chars[pos] != other_chars[pos]:
                    differing_positions.append(pos)

            # Handle case where one filename is a prefix of another
            if len(chars) != len(other_chars):
                for pos in range(min_len, max(len(chars), len(other_chars))):
                    differing_positions.append(pos)

        # Get all unique positions where this file differs from others
        unique_positions = sorted(set(differing_positions))

        if not unique_positions:
            unique_identifiers[filename] = (
                ""  # No unique part (should not happen with different filenames)
            )
            continue

        # Find consecutive ranges of differing positions
        ranges = []
        current_range = [unique_positions[0]]

        for pos in unique_positions[1:]:
            if pos == current_range[-1] + 1:
                current_range.append(pos)
            else:
                ranges.append(current_range)
                current_range = [pos]

        ranges.append(current_range)  # Add the last range

        # Find the most significant range (typically the longest or most meaningful)
        # For simplicity, we'll use the longest range
        longest_range = max(ranges, key=len)

        # Extract the unique substring
        start = longest_range[0]
        end = longest_range[-1] + 1
        unique_substring = filename[start:end]

        # Clean up the substring - remove partial words or patterns
        # This attempts to find natural word boundaries around the unique part
        expanded_start = start
        expanded_end = end

        # Expand backward to include a word boundary or common delimiter
        while expanded_start > 0 and filename[expanded_start - 1] not in [
            " ",
            "_",
            "-",
            ".",
        ]:
            expanded_start -= 1

        # Expand forward to include a word boundary or common delimiter
        while expanded_end < len(filename) and filename[expanded_end] not in [
            " ",
            "_",
            "-",
            ".",
        ]:
            expanded_end += 1

        # Use the expanded substring if it's not too much longer
        if (expanded_end - expanded_start) <= 2 * (end - start):
            unique_substring = filename[expanded_start:expanded_end]

        # remove extension
        unique_substring = unique_substring.split(".")[0]

        unique_identifiers[filename] = unique_substring.strip()

    return unique_identifiers


def accumulate_distances(trace_data):
    """
    Calculate pairwise distances between barcodes across all traces.

    Parameters:
    ----------
    trace_data : ChromatinTraceTable
        Table containing trace data with barcode positions

    Returns:
    -------
    dict
        Dictionary mapping barcode pairs to their median distances
    """
    distances = defaultdict(list)
    trace_groups = trace_data.group_by("Trace_ID").groups

    for trace in trace_groups:
        barcode_positions = {
            row["Barcode #"]: np.array([row["x"], row["y"], row["z"]])
            for row in trace
            if not np.isnan(row["x"])
        }
        for bc1, bc2 in combinations(barcode_positions, 2):
            p1 = barcode_positions[bc1]
            p2 = barcode_positions[bc2]
            dist = np.linalg.norm(p1 - p2)
            key = tuple(sorted((bc1, bc2)))
            distances[key].append(dist)

    return {key: np.median(vals) for key, vals in distances.items()}


def compare_distance_maps(distance_maps):
    """
    Calculate Pearson correlation between distance maps from different trace files.

    Parameters:
    ----------
    distance_maps : dict
        Dictionary mapping filenames to their distance maps

    Returns:
    -------
    tuple
        (files, corr_matrix) where files is a list of filenames and
        corr_matrix is the correlation matrix
    """
    files = list(distance_maps.keys())
    all_keys = set(
        itertools.chain.from_iterable([dm.keys() for dm in distance_maps.values()])
    )

    # Create vectors for each file over all keys
    vectors = {}
    for fname in files:
        vec = []
        for key in sorted(all_keys):
            vec.append(distance_maps[fname].get(key, np.nan))
        vectors[fname] = np.array(vec)

    # Build correlation matrix
    n = len(files)
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            x = vectors[files[i]]
            y = vectors[files[j]]
            mask = ~np.isnan(x) & ~np.isnan(y)
            if np.sum(mask) > 1:
                corr, _ = pearsonr(x[mask], y[mask])
            else:
                corr = np.nan
            corr_matrix[i, j] = corr

    return files, corr_matrix


def plot_correlation_matrix(
    files, matrix, output_filename="trace_correlation_matrix.png", vmin=-10, vmax=10
):
    """
    Plot a correlation matrix between files with unique identifiers as labels.

    Parameters:
    ----------
    files : list
        List of filenames to use
    matrix : numpy.ndarray
        Square correlation matrix of the files
    output_filename : str
        Name of the output file (default: trace_correlation_matrix.png)

    Returns:
    -------
    None
        Saves the plot as a PNG file
    """
    # Get unique identifiers for each file
    unique_identifiers = find_unique_substrings(files)

    # Create labels for the plot
    labels = [os.path.basename(unique_identifiers[f]) for f in files]

    title_fontsize = 16
    label_fontsize = 12
    tick_fontsize = 10

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the matrix with dynamic color range
    if vmin == -10:
        vmin = np.min(matrix)
    if vmax == 10:
        vmax = np.max(matrix)
    im = ax.imshow(matrix, cmap="RdBu", interpolation="nearest", vmin=0, vmax=1)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pearson Correlation", fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    # Set tick labels
    ax.set_xticks(range(len(files)))
    ax.set_yticks(range(len(files)))
    ax.set_xticklabels(labels, rotation=90, fontsize=tick_fontsize)
    ax.set_yticklabels(labels, fontsize=tick_fontsize)

    # Add axis labels
    ax.set_xlabel("Files", fontsize=label_fontsize)
    ax.set_ylabel("Files", fontsize=label_fontsize)

    # Add title
    ax.set_title("Trace Table Similarity Matrix", fontsize=title_fontsize)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"$ Saved correlation matrix as {output_filename}")

    np.save(os.path.splitext(output_filename)[0] + ".npy", matrix)
    print(
        f"$ Saved correlation matrix data in NPY format: {os.path.splitext(output_filename)[0] + '.npy'}"
    )

    plt.close()


def main():
    print_script_banner(__file__, __doc__)
    """
    Main function that executes the trace comparison workflow.
    """
    parser = parse_arguments()
    args = parser.parse_args()
    trace_files = get_trace_files(args)

    print(f"Analyzing {len(trace_files)} trace files...")

    # Calculate distance maps for each file
    distance_maps = {}
    for fpath in trace_files:
        print(f"Processing {os.path.basename(fpath)}")
        trace = ChromatinTraceTable()
        if args.verbose:
            trace.load(fpath)
        else:
            silent_load_trace(trace, fpath)
        distance_maps[fpath] = accumulate_distances(trace.data)

    # Compare distance maps and generate correlation matrix
    files, corr_matrix = compare_distance_maps(distance_maps)

    if args.verbose:
        print("\nPearson Correlation Matrix:")
        print(corr_matrix)

    # Plot and save the correlation matrix
    plot_correlation_matrix(
        files,
        corr_matrix,
        output_filename=os.path.join(
            args.output, f"trace_correlation_matrix.{args.output_format}"
        ),
        vmin=args.vmin,
        vmax=args.vmax,
    )


if __name__ == "__main__":
    main()
