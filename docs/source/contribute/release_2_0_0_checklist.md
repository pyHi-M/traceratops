# 2.0.0 release assessment and checklist

## Baseline and version recommendation

Assessed on 2026-10-28. The latest published GitHub release is
[v1.0.1](https://github.com/pyHi-M/traceratops/releases/tag/v1.0.1), published
2026-07-07, at `d6baa2f4c850143f19058e8993917c8d86ef587b`.
The inspected `main` tip is `261a48fe92498abcb154891fe6e5a5ca8219b610`,
40 commits after that tag. The release preparation starts at `dev`
`c9cdc058b7525fc660119521020e117fffa5d043`, which includes merged security
PR #205. The pre-preparation version remains 1.0.1 despite these changes.

Main and dev have diverged histories: 18 main-only and 19 dev-only commits.
Some changes have been merged/squashed separately; commit counts are not feature
counts. The release inventory uses the v1.0.1-to-dev source diff and the 41
post-tag commits reachable from dev, with their merged PR descriptions.
The dev-to-main diff was inspected separately. Do not overwrite main or assume
fast-forward integration; the final integration PR must review the actual merge.

**Recommend 2.0.0**, with publication tag **v2.0.0**. The root changelog and
`traceratops/_version.py` explicitly prescribe Semantic Versioning. Dropping
Python 3.9/3.10, removing a command and changing some filename/selection behavior
justify a major version. Although earlier releases grouped substantial changes
under 1.0.1, this recommendation follows the declared convention rather than
continuing that patch-level practice. The release date remains Unreleased until
publication timing is approved. Only pyproject.toml stores the package version;
Sphinx and runtime version reporting read installed metadata.

## Verified user-facing inventory

| Area/tool | Status relative to v1.0.1 | Evidence |
| --- | --- | --- |
| Physical/genomic distances | New CLI, binned axis/3D summaries, power-law plots, independent CSV/figure output controls | #146, #148, #168; `trace_physical_vs_genomic_distance.py` |
| trace_sorter | New whole-trace sorting by size or radius of gyration; documented development status | #175; `trace_sorter.py` |
| trace_splitter | Existing K-means extended with HDBSCAN, split-all, beam resolution, threshold/likelihood multiplicity classification, reliability diagnostics, candidate ambiguity and progress/quiet controls | #177, #183, #184, #189; `core/trace_resolver/`, `trace_splitter.py` |
| trace_filter | Existing tool extended with combined localization-quality thresholds and filtered localization output; roundness/spot-pixel-percentage use upper bounds | #130–#132; `filter_by_localization_metrics` |
| trace_analyzer | Existing tool; genomic neighbor validation, adjacent-distance hexbin maps and plot layout changes | #127, #128; `plot_neighbor_distances` |
| trace_to_matrix | Existing tool; vectorized minimum mode and optional per-trace parallelism | #125; `core/build_matrix.py` |
| plot_him_matrix | Existing tool; canonical input/colormap aliases and rendering improvements | #128, #140, #161 |
| plot_bootstrapping | Existing tool; canonical flags and shared Hi-M rendering/output handling | #161, #163 |
| plot_compare2matrices | Existing tool; canonical flags and common comparison renderer | #161, #204 |
| trace_3way_coloc | Existing calculation; corrected partner labels/anchor guides and reused presence arrays for contacts/bootstraps | #159; `trace_3way_coloc.py` |
| plot_3way_coloc | Existing replot tool; no new calculation claimed; existing output-format support predates v1.0.1 | Source diff: import ordering only |
| plot_4m | Existing tool; not a new capability in this release | Source diff: import ordering/permissions only |
| localization_analyzer | Already introduced before v1.0.1; QC documentation/format handling improved | v1.0.1 changelog; #143, #144; module diff is import ordering |
| trace_pearsons | Existing statistical comparison; explanatory documentation improved, no new statistical method | Module diff: docstring/import ordering |
| trace_assign_mask / trace_split_labels | Existing tools; no new scientific behavior claimed | Module diffs: import ordering/permissions; missing tifffile dependency corrected during preparation |
| FOF-CT/localization I/O | Declared XYZ-unit conversion, `.dat` input deprecation and `.ecsv` output enforcement | #129, #144; core table readers/writers |
| PDB/PyMOL | Fixed first-N export and selection; new PyMOL helper and tutorial | #141, #142, #175 |
| Documentation | Data schemas, renamed/reorganized notebooks, plotting and QC guidance | #143, #145, #151, #171, #172 |
| Dependencies/security/install | Python >=3.11, Pillow >=13, uv install, dev-only tools, no lockfile, nine test jobs, separate installed-environment security audits | #137, #205 |
| Internal maintenance | Import order, executable flags and common plotting implementation are not separate scientific features | #138, #164, #204 |

No performance factor or blanket reproducibility improvement is inferred from
refactoring. Beam-search confidence is a non-probabilistic score gap, and the
candidate mode is opt-in. K-means and threshold classification remain defaults.

## Open work excluded from this candidate

At assessment, PR #201 (`development/trace_curator` → dev) is open and unmerged.
It proposes an empirical curator API/production CLI; neither exists in the
candidate tree. PRs #192, #194, #195 (draft benchmark work) and #185 (simulator)
target feature/trace_curator and are also unmerged. These are not release-note
features or silently merged prerequisites. Maintainers should explicitly confirm
this release scope. If #201 is required, defer the release, merge it through its
own review, then repeat inventory, tests, audits and packaging validation.

## Compatibility and migration

- Require Python >=3.11; do not install this release on Python 3.9 or 3.10.
- Replace removed `trace_filter_advanced` invocations using documented filtering
  and splitting commands; validate intended behavior rather than mechanically
  translating every legacy flag.
- Expect localization `.ecsv` outputs, including when a `.dat` output was
  requested. `.dat` inputs emit deprecation warnings.
- Supply consistent Chrom/Chrom_Start/Chrom_End annotations for genomic-neighbor
  analysis. Existing datasets with placeholder coordinates need appropriate
  annotations for scientifically meaningful ordering.
- Nanometer FOF-CT coordinates are converted to micrometers. Corrected numerical
  distances can differ from results obtained with the old import behavior.
- Set PDB-export selection explicitly: first ten is now the CLI default, and
  the prior minimum-barcode prefilter is not applied in that export path.
- Canonical plotting flags retain documented aliases, but shared-renderer figure
  appearance, bootstrap font defaults, and some generated filenames differ.
- uv uses pyproject.toml without a universal lockfile; fresh installations may
  resolve different compatible dependency releases. Scikit-learn >=1.3 is needed
  for the integrated HDBSCAN implementation. Pillow >=13.0.0 is the security floor.
- PyMOL is external software. Its interactive visualization is not covered by
  the automated Python regression suite.

## Validation evidence

Validation results are recorded in the release-preparation PR. The starting dev
commit's [nine-job test run](https://github.com/pyHi-M/traceratops/actions/runs/40635248930)
passed on Python 3.11/3.12/3.13 across Ubuntu, macOS and Windows.

Release-preparation validation must separately cover:

- Full pytest regressions, including trace I/O, filtering, splitting/resolution,
  distance/contact matrices, three-way contacts and PNG/SVG/PDF exports.
- Every declared installed console entry point, with CLI help distinguished from
  actual numerical tests; targeted mask TIFF/NPY loading after dependency repair.
- Wheel and sdist build, Twine metadata checks, metadata/entry-point/payload
  inspection, clean wheel installation outside the source checkout and installed
  distribution version agreement.
- Runtime and development pip-audit reports, including the added tifffile
  dependency; errors and vulnerable findings must not be suppressed.
- Candidate CI and changed-file pre-commit checks. Existing static-analysis or
  documentation-build warnings must be disclosed rather than attributed to
  scientific validation.

The first fresh-install check found trace_assign_mask could not import tifffile.
The preparation adds this already-required runtime dependency and its Sphinx
mock entry; it does not modify mask algorithms. Installed wheel checks also found a missing
FOF-CT export main callable and a missing checkout-only COPYRIGHT.txt metadata
file. Preparation adds the wrapper and reads author names from installed package
metadata, without changing coordinate conversion. A real installed-console
export/import test validates coordinates and identities. The initial distribution build
and Twine checks passed; setuptools warned that the existing table-form license
metadata is deprecated and must be modernized before its 2027-02-18 deadline.
This warning does not invalidate the build and is not changed as unrelated
license metadata work.

## Prepared dev-to-main PR

Open this only after the release-preparation PR has been reviewed and merged.

**Title:** Release traceratops 2.0.0

**Description:**

Release 2.0.0 collects the verified changes since v1.0.1: barcode-aware trace
resolution and multiplicity classification, HDBSCAN splitting, physical/genomic
distance analysis, whole-trace sorting, localization-quality filtering, matrix
construction and visualization improvements, FOF-CT unit correction and PDB
export fixes. The unmerged curator/simulator/benchmark work is excluded.

Breaking changes include Python >=3.11 (dropping 3.9/3.10), removal of
trace_filter_advanced, enforced localization .ecsv output, genomic-neighbor
annotation requirements and changed PDB selection. uv installation, Pillow >=13,
explicit tifffile and scheduled runtime/development security audits accompany
the release. See the detailed changelog and prepared GitHub release notes.

Attach final release-preparation CI, scientific test, security audit and package
validation results, then verify this PR's own checks on the actual merge result.
Call out unresolved static-analysis warnings, manual PyMOL/private-data QA and
scope decisions. Do not merge automatically. Merge into main, tagging and
publication require explicit approval. No tags/releases/packages were published
as part of preparation.

## Approval and publication checklist

- [ ] Confirm 2.0.0 and exclusion of unmerged PR #201 and related work.
- [ ] Review validation evidence and any remaining readiness limitations.
- [ ] Approve and merge the preparation PR into dev.
- [ ] Open/finalize dev → main PR with its own passing CI/audits and actual release date.
- [ ] Explicitly approve the main merge; verify the resulting main SHA/version.
- [ ] Explicitly approve creating v2.0.0 on that main SHA; check for duplicates.
- [ ] Explicitly approve publication with the prepared notes and tagged-source artifacts.
- [ ] Verify tag, published release, metadata and artifacts; no PyPI upload.
- [ ] Confirm security-audit.yml is on default main, dispatch an audit and verify scheduling/notifications.

Follow the exact sequence in the [release guide](release_guide.md).

### Local readiness findings

The full starting-candidate regression suite passed 187 tests with six warnings.
A targeted FOF-CT installed-command round-trip is added for the packaging repair.
The final PR checks must include that test and the final dependency declarations.
Both development and wheel-runtime dependency audits found no known vulnerabilities,
including tifffile 2026.3.3 on Python 3.11.

`mypy traceratops` reported 90 errors in 19 files before the export repair. These
are existing scientific-source typing findings; they are not a passing static
check and have not been silently waived. Codespell also flags existing source
identifiers (for example Buid and CONECT) and spellings. New release prose passed
the spelling check. Maintainers must triage these findings against the main-PR
checklist rather than treating numerical test success as static-analysis success.
No broad typing/spelling cleanup is included. The preparation PR remains draft
pending readiness review, even if automated scientific/packaging checks pass.
Manual interactive PyMOL and private acquisition-data QA were not executed.
