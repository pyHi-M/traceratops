"""Helpers for printing a standard startup banner for command-line scripts."""

from __future__ import annotations

from pathlib import Path

_SCRIPT_DESCRIPTIONS = {
    "collect_files.py": "Collect one file per subdirectory using pattern matching.",
    "localization_analyzer.py": "Load a localization table and plot localization quality-control distributions.",
    "localization_merge.py": "Merge multiple localization files into a single consolidated output.",
    "plot_3way_coloc.py": "Replot 3-way co-localization matrices from .npy files.",
    "plot_4m.py": "Compute barcode interaction frequencies with bootstrapping.",
    "plot_bootstrapping.py": "Perform bootstrapping analysis on PWD distance matrices.",
    "plot_compare2matrices.py": "Compare two proximity matrices.",
    "plot_him_matrix.py": "Calculate and plot PWD and proximity matrices.",
    "plot_matrix_comparison.py": "Compare PWD matrices from two experiments.",
    "trace_3way_coloc.py": "Compute three-way co-localization frequencies with bootstrapping.",
    "trace_analyzer.py": "Analyze chromatin trace files.",
    "trace_assign_mask.py": "Assign mask labels to traces using numpy masks.",
    "trace_filter.py": "Filter chromatin trace files.",
    "trace_filter_advanced.py": "Apply advanced filtering to chromatin trace files.",
    "trace_genomic_coordinates.py": "Assign genomic coordinates to a chromatin trace table.",
    "trace_import_from_fofct.py": "Convert a FOF-CT CSV file to a pyHiM trace table in ECSV format.",
    "trace_merge.py": "This script will merge trace tables provided as inputs.",
    "trace_pearsons.py": "Compare chromatin trace tables by computing pairwise distances.",
    "trace_plot.py": "Plot one or multiple traces in 3D.",
    "trace_split_labels.py": "Split a trace file into two files based on the presence of a label.",
    "trace_splitter.py": "Split chromatin traces using K-means clustering when radius of gyration exceeds a threshold.",
    "trace_stats.py": "Compute basic statistics for chromatin trace files.",
    "trace_to_matrix.py": "Convert a trace file to a matrix using pyHiM core routines.",
}


def print_script_banner(script_file: str, description: str | None) -> None:
    """Print a standard banner with the script name and a short description."""
    script_name = Path(script_file).name
    summary = _SCRIPT_DESCRIPTIONS.get(script_name, "") or _first_description_line(
        description
    )

    print(f"------- Running {script_name} --------")
    if summary:
        print(summary)
    print()


def _first_description_line(description: str | None) -> str:
    """Return the first non-empty line of a script docstring."""
    if description is None:
        return ""

    for line in description.splitlines():
        stripped_line = line.strip()
        if stripped_line:
            return stripped_line

    return ""
