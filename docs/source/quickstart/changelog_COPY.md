# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0] - Unreleased

The major version reflects compatibility changes, including Python >=3.11.

### Added
- `trace_physical_vs_genomic_distance` computes binned physical distances against genomic separation, with per-axis/3D summaries, CSV output, and log-log power-law fits with confidence intervals (#146, #148, #168).
- `trace_sorter` sorts whole traces by localization count or radius of gyration, preserving within-trace row order (#175). This tool is documented as development status.
- `trace_splitter` supports HDBSCAN and `--split-all`, in addition to the existing K-means method (#177).
- Opt-in barcode-aware beam resolution in `trace_splitter` can classify and clean candidate singlets, separate candidate doublets, and leave incompatible detections unassigned. Threshold and likelihood multiplicity classifiers, reliability warnings, and diagnostic ECSV outputs are available (#183).
- Opt-in candidate-specific duplicate ambiguity scoring records localization scores and decisions, preserving `Spot_ID` (#184). Score-gap confidence is not a calibrated probability.
- `trace_filter` accepts combined localization-quality thresholds and saves the localizations retained in the final traces (#130–#132).
- PyMOL helper commands and a tutorial support displaying exported PDB traces in a grid and coloring barcodes (#141). PyMOL is a separate application, not an installed runtime dependency.

### Changed
- **Python 3.11 or newer is required. Python 3.9 and 3.10 are no longer supported** (#205).
- **`trace_filter_advanced` is removed**; update workflows to the documented `trace_filter` and `trace_splitter` commands (#156). There is no automatic one-to-one migration for every legacy option.
- Localization outputs use `.ecsv`; legacy `.dat` inputs remain readable with a deprecation warning and `.dat` output requests are rewritten to `.ecsv` (#144).
- Neighbor-distance analysis in `trace_analyzer` requires consistent `Chrom`, `Chrom_Start`, and `Chrom_End` annotations and counts neighbors in genomic order. Previous/next neighbor distances are shown as hexbin density maps (#127, #128).
- `trace_to_matrix --n_jobs` enables per-trace parallel matrix construction; vectorized minimum-distance handling preserves duplicate-barcode semantics (#125). No general speedup factor is claimed.
- Plotting commands accept consistent input/output, barcode, colormap and color-limit flags, retaining documented legacy aliases (#161). Bootstrap matrices use the shared Hi-M plotting utility (#163).
- Matrix comparison difference/ratio, mixed-triangle and rank-sum plots use the shared matrix renderer; colors and appearance change, including the rank-sum colormap (#204).
- `trace_plot` now exports the first ten traces by default; use `--number_traces`, `--selected_trace`, or `--all` explicitly for scripted selection (#142, #175). The former minimum-barcode prefilter is no longer applied in this export path.
- `trace_physical_vs_genomic_distance` separates `--data_output` CSV naming from `--output`/`-O` figure naming and supports PNG/SVG/PDF export (#168).
- uv is the recommended installer. Development installation uses `uv pip install -e ".[dev]"`; `pyproject.toml` is authoritative and `uv.lock` is removed and ignored (#137, #205).

### Fixed
- Restore the installed `trace_export_to_fofct` entry point and obtain author metadata from the installed distribution instead of a missing checkout-only COPYRIGHT.txt file. A console-command round-trip test checks numerical coordinates and identities.
- FOF-CT import respects declared XYZ units, converting nanometers to micrometers and retaining micron coordinates (#129). Re-importing nanometer data can therefore correct distances compared with earlier results.
- `trace_plot --number_traces` produces the requested first-N PDB exports instead of no output (#142).
- Three-way contact plots use the correct partner-barcode labels and anchor guides; contact/bootstrap calculations reuse per-trace presence arrays while retaining the documented denominator (#159).
- Colormap access in trace plotting is compatible with maintained Matplotlib versions (#173).
- Roundness and spot-pixel-percentage quality filters use maximum thresholds; other documented quality metrics use minimum thresholds (#131).

### Security and maintenance
- Declare the existing mask-assignment command's `tifffile` runtime dependency so a clean installation can import and run `trace_assign_mask`.
- Require Pillow >=13.0.0 for reported loader/decoder, Ghostscript, TIFF, GIF, FITS and PDF security fixes; retain Matplotlib and Pillow (#205).
- Keep pytest>=8.4.2 and pre-commit>=4.3.0 in development extras only (#205).
- Test Python 3.11/3.12/3.13 on Linux/macOS/Windows using uv. Add weekly, dependency-related PR, and manual pip-audit workflows with separate runtime/development reports (#205).
- Expand scientific regression, CLI argument and PNG/SVG/PDF export coverage. Progress bars and quiet controls are available for trace splitting (#189).
- Reorganize tutorials and document ECSV schemas, quality thresholds, physical/genomic distances, plotting, and installation (#143, #145, #151, #171, #172).
- Internal import ordering, executable permissions and shared plotting implementation changes are not new scientific tools (#138, #164, #204).

[Unreleased]: https://github.com/pyHi-M/traceratops/compare/v2.0.0...dev
[2.0.0]: https://github.com/pyHi-M/traceratops/compare/v1.0.1...v2.0.0

## [1.0.1]

### Added
- `collect_files` script to gather localization files from nested acquisition folders
- `trace_split_labels` script to split traces by assigned mask labels
- `localization_analyzer` script for localization quality-control analysis
- Support for assigning barcode genomic coordinates from BED files in `trace_genomic_coordinates`
- Pipe/stdin input support across trace-processing scripts, including `trace_stats`
- New tutorials for trace merging, quality control, filtering, label splitting, matrix visualization, multi-way colocalization, dataset comparison, and pipe usage
- Test coverage for localization tables, trace filtering, trace merging, matrix plotting, matrix conversion, and file collection

### Fixed
- NumPy 2 compatibility for NaN handling in matrix and plotting workflows
- Localization table formatting and comment-field handling when reading, merging, and updating localization tables
- Missing comments in chromatin trace tables
- `trace_filter --clean_spots` test discovery and behavior coverage
- FOF-CT import/export compatibility and standardized input option naming
- Multi-worker output collisions by including the trace file name in generated output paths

### Changed
- Refactored localization table handling and homogenized command-line arguments across scripts
- Renamed `trace_impute_genomic_coordinates` to `trace_genomic_coordinates`
- Default spot-table output format is now `ecsv`
- `plot_3way_coloc` can write plots in user-selected output formats
- Updated documentation structure, script pages, pipe examples, release guidance, and project logo

### Removed
- `trace_correct_coordinates`: superseded by `trace_genomic_coordinates`
- `localization_cp_files`: superseded by `collect_files`
- `analyze_localizations`: replaced by `localization_analyzer`
- `plot_single_cell`: unused & deprecated script
- `plot_n_him_matrices`: unused & deprecated script

## [0.5.0] - 2025-04-17

### Added
- This CHANGELOG.md file

### Fixed
- figure_him_matrix: Apply delta margin to avoid massive tick generation for narrow color ranges (Matplotlib <3.7)

### Changed
* rename `figure_him_matrix` to `plot_him_matrix`
* rename `compare_pwd_matrices` to `plot_matrix_comparison`
* rename `plot_3way` to `trace_3way_coloc`
* rename `replot_3way` to `plot_3way_coloc`

### Removed
- process_him_matrix: unused & deprecated script
- trace_combinator: unused & deprecated script
- pwd_matrix_2_pdb: unused & deprecated script


## [0.4.0] - 2025-04-09

### Added
- Doc: Reliability status for each script (development < beta-test < stable)
- trace_filter: [localization intensity filtering](https://traceratops.readthedocs.io/en/latest/tutorials/filter_intensity.html)
- Import from pyHiM **post-processing** scripts:
    * trace_to_matrix
    * trace_stats
    * trace_splitter
    * trace_plot
    * trace_pearsons
    * trace_genomic_coordinates
    * trace_import_from_fofct
    * trace_export_to_fofct
    * trace_correct_coordinates
    * trace_combinator
    * trace_assign_mask
    * trace_analyzer
    * pwd_matrix_2_pdb
    * process_him_matrix
    * localization_merge
    * localization_cp_files
    * compare_pwd_matrices
    * analyze_localizations
- Import from pyHiM **plot** scripts:
    * replot_3way.py
    * plot_bootstrapping.py
    * plot_4m.py
    * plot_3way.py
    * plot_single_cell.py
    * plot_n_him_matrices.py
    * plot_compare2matrices.py

## [0.3.0] - 2025-03-25

### Added
- `trace_merge` script (import from pyHiM)
    - prevent ROI number merge conflict

### Fixed
- figure_him_matrix: change default proximity normalization (By default, NaN values per bin are removed before compute statistics for proximity).

### Changed
- figure_him_matrix: `--remove_nan` become `--keep_nan`

## [0.2.0] - 2025-03-20

### Added
- Installation guide
- `trace_filter` script (import from pyHiM)
- [Documentation available online](https://traceratops.readthedocs.io)
- `figure_him_matrix` script (import from pyHiM)
    - nan_matrix: plot with him_matrix the nan percentage
    - add a threshold based on NaN percentage [#26](https://github.com/pyHi-M/traceratops/issues/26)

### Fixed
- Fix the `--clean_spots` option of `trace_filter` to ensure that a spot_ID is unique to a single trace_ID
- Windows compatibility
- figure_him_matrix: fix `--shuffle` option

### Changed
- `figure_him_matrix`: `--norm` option become `--remove_nan`

## [0.1.0] - 2025-03-07

### Added
- Documentation can be build locally
- Project can be installed locally with `pip install .` or `pip install .[dev]`

## [Template] - YYYY-MM-DD

### Added

### Fixed

### Changed

### Removed
