# ℹ️ Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
    * trace_impute_genomic_coordinates
    * trace_import_from_fofct
    * trace_export_to_fofct
    * trace_filter_advanced
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
