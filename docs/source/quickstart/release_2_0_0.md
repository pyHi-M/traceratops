# Traceratops 2.0.0

Traceratops 2.0.0 extends chromatin-trace quality control, barcode-aware trace
resolution, and physical-versus-genomic distance analysis. It also updates
installation and dependency security. Existing scripts should be reviewed for
the compatibility changes below before reprocessing data.

## Highlights

- Resolve candidate merged traces with an opt-in barcode-aware beam-search
  method, with threshold or likelihood singlet/doublet classification,
  unassigned-localization handling, and diagnostic tables.
- Relate physical distances to genomic separation with the new
  `trace_physical_vs_genomic_distance` command, CSV summaries and power-law plots.
- Order complete traces by localization count or radius of gyration with the new
  `trace_sorter` command, which retains development status.

## New tools and improvements

`trace_splitter` adds HDBSCAN, `--split-all`, barcode-aware beam resolution,
and optional per-candidate ambiguity diagnostics for duplicate barcode cleanup.
The default method remains K-means. Beam resolution and likelihood classification
must be selected explicitly. Score-gap confidence is not a calibrated probability;
reliability warnings and documented limitations should guide interpretation.

`trace_to_matrix --n_jobs` supports parallel per-trace matrix construction.
Plotting commands use more consistent options while retaining documented legacy
aliases. Bootstrap and matrix-comparison plots share the Hi-M rendering utility.
Three-way contact plots have corrected barcode labels and anchor guides. PDB
export has explicit first-N, selected-trace and all-trace modes, with a new PyMOL
helper and tutorial for viewing exported traces.

## Quality control and reproducibility

`trace_filter` combines localization-quality thresholds and writes a localization
table matching the spots retained in the final traces. Roundness and
spot-pixel-percentage filters use maximum thresholds.

`trace_analyzer` uses genomic annotations for neighbor-distance analysis and
shows previous/next-neighbor density maps. FOF-CT import now converts declared
nanometer coordinates to micrometers; this fixes previously incorrect distances
for those inputs. Regression tests cover matrix construction, trace processing,
three-way contact calculations, and representative PNG/SVG/PDF figure exports.

## Installation and compatibility

- **Python 3.11 or newer is required. Python 3.9 and 3.10 are no longer supported.**
- Python 3.11, 3.12 and 3.13 are tested on Linux, macOS and Windows.
- uv is the recommended installer. In a cloned repository, create an environment
  with `uv venv --python 3.11`, then install with `uv pip install -e ".[dev]"`.
- Dependencies are specified in `pyproject.toml`; no universal `uv.lock` is
  committed. Pillow >=13.0.0 addresses the reported security advisories.
- Automated dependency security monitoring audits runtime and development
  environments weekly and on dependency-related pull requests.
- `trace_filter_advanced` has been removed. Migrate to the documented
  `trace_filter` and `trace_splitter` workflows; not every legacy option has a
  direct replacement.
- Localization outputs use `.ecsv`; `.dat` input remains readable with a
  deprecation warning. Update scripts that expect `.dat` output filenames.
- Neighbor-distance analysis needs consistent `Chrom`, `Chrom_Start` and
  `Chrom_End` annotations. Figure styling and some output names have changed.
- `trace_plot` exports the first ten traces by default and no longer applies its
  former minimum-barcode prefilter. Set selection options explicitly in scripts.

## Bug fixes and maintenance

Corrected first-N PDB export, FOF-CT unit conversion, three-way plot annotation,
maintained-Matplotlib colormap access, and the missing `tifffile` dependency
needed by the mask-assignment command. The installed FOF-CT export command now
has a callable entry point and reads author metadata from the installed package,
so export no longer depends on a missing checkout-only file. Tutorials and ECSV format references
have been updated. Development tools are no longer runtime dependencies.

The unmerged empirical trace curator (PR #201), simulator and exploratory
curator benchmarks are **not included** in this release.

## Full changelog

[Detailed changelog](changelog_COPY.md) ·
[Full comparison: v1.0.1...v2.0.0](https://github.com/pyHi-M/traceratops/compare/v1.0.1...v2.0.0)
