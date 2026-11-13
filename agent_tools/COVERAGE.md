# Coverage and limitations

Reviewed snapshot: `pyHi-M/traceratops`, branch `dev`, commit
`d767666dfff6be0eca410cca1cc6acc79970140e`, package 2.0.0.
Verification date: 2026-11-13. No scientific source, CLI interface, packaging
configuration, or existing expected test result was changed.

## Coverage

| Surface | Discovered / cataloged / verified |
|---|---|
| Registered executable names | **29 / 29 / 29 executed examples** |
| Distinct implementations | **25 / 25 / 25 executed examples** |
| Deprecated aliases | 4, each separately cataloged and executed |
| Analytical Python callables | **76 cataloged**, source/signature reviewed; not all individually execution-tested |
| Public-name Python inventory | 310 function/method candidates classified in `python_surface.yaml` |
| CLI example syntax | 29/29 parsed with exact flags, abbreviation disabled during validation |
| Syntax-only CLI examples | 0 |

The Python count includes functions and methods. A source docstring counts as
source documentation, not a stability promise. Public analytical methods on
explicitly exported resolver classes and documented PyMOL commands are also
included. Parser factories, launchers, I/O-only lifecycle methods, filename/report
utilities, private names and unexported resolver internals are not recommended
analysis APIs. Every excluded candidate has a reason in the inclusion ledger.
The package root itself exports no public analysis functions. The tutorials
primarily demonstrate CLI commands rather than a consolidated Python API.

Discovery used `pyproject.toml`, all package Python modules, documented executable
scripts, source docstrings, dedicated command pages, tutorial notebook code and
existing tests. All installed console entry points match the checkout. PyMOL's
script helpers are Python/PyMOL operations, not additional installed CLI names.
Project installation/release/documentation scripts are maintenance infrastructure.

## Validation results

- `python agent_tools/validate.py`: **PASS**. Checks complete command inventory,
  installed distribution entry points, exact argparse actions/defaults/choices/types/
  arity/groups, referenced source paths, source fingerprints, example syntax,
  tool/task/example/schema references, duplicate names/IDs/mapping keys, scientific
  descriptions, Python inclusion inventory/signatures/docstrings, and the two
  critical scientific routing regressions.
- `python -m pytest agent_tools/test_validation.py -q`: **15 passed**. Mutation
  cases prove detection of missing/invented commands, changed parser fields,
  invalid example flags, broken task/example references, changed source/signatures,
  absent capabilities and duplicate mappings.
- `python agent_tools/examples/run_examples.py --work-dir /tmp/traceratops-examples`:
  **29 passed**. Each runs an installed executable in a fresh directory. Checks
  include expected nonempty artifacts, selected semantic assertions and unchanged
  input fixture bytes. Plot-only tools receive artifact checks, not exhaustive
  visual-content validation. Evidence: `examples/verification.json`.
- Existing suite with the virtual environment's executable directory on PATH:
  **188 passed, 3 warnings**. Initial invocation without this PATH had 41 executable-
  lookup failures and 147 passes; the environment correction resolved them. No test
  expectations or package code were edited. Remaining warnings concern all-NaN
  ensemble reduction and Trace_ID string truncation in existing splitter tests.

Use an installed editable checkout: `python` above must be the interpreter for
that installation. Example runner uses its scripts directory, so it does not rely
on a shell activation. Test execution here used Python 3.11.16. Runtime dependency
versions are recorded in `examples/verification.json`. Nothing was installed for
this assignment; validation uses the standard library and the project's existing
runtime dependencies. Tests use the project's existing pytest dependency.

## Data and verification boundaries

No experimental data were required for the final CLI smoke tests. Small synthetic
trace/localization tables, genomic BED, metadata and masks are generated on demand;
the repository's 4-barcode matrix fixture is reused. Suitable synthetic bootstrap
input is generated through `BuildMatrix.build_distance_matrix`, and replot input
through `compute_threeway_colocalization`. The sparse repository fixture initially
failed bootstrapping because a pair had fewer than two finite observations; this
is a data requirement, not an absent feature.

The final suite is not exhaustive across CLI parameters, alternate file formats,
large biological datasets, stochastic estimators or operating systems. Bootstrap
cycles are deliberately small and no exact random results are asserted. Spatial
splitting methods, beam likelihood submodes, intensity/quality combinations, TIFF
masks, every visualization option and external PyMOL execution are not all covered
by these examples. Existing tests provide additional coverage but do not make the
catalog a biological validation study. No CLI example remains blocked for missing
data. PyMOL is not installed as a required traceratops dependency; its Python
functions were reviewed without importing or executing PyMOL.

The generic trace schema lists formats handled by the trace reader. Individual
commands may require ECSV explicitly (notably FOF-CT export), and downstream plots
need specific columns even when a loader can read the file. Source documentation
for many Python functions is sparse or historically inconsistent; null documented
return types are intentional. Parameters without argparse choices can still have
restricted runtime meanings. Runtime fallback source rules are evidence, not an
exhaustive symbolic analysis of arbitrary code.

## Important source/documentation discrepancies

| Operation | Source-verified behavior and implication |
|---|---|
| `trace_stats` | Counts unique barcode identities as well as traces/ROIs; module summary omits barcode count. Counts Trace_ID globally, not (ROI,Trace_ID). |
| `trace_analyzer` | Two trace-statistics histograms distinguish row counts and unique barcode counts. Histogram bins omit the upper end for some datasets; small inputs need care. `--rootFolder` does not relocate outputs. `--output_format` is valid; the documented `--format` is not. |
| Analyzer neighbor QC | Requires consistent Chrom/Chrom_Start/Chrom_End. Neighbors must be consecutive both numerically and in genomic order; no missing-column fallback. “KDE projections” are actually centered hexbin plots. Detection uncertainty is displayed as bootstrap violins. |
| `trace_to_matrix` | Uses minimum repeated-barcode pair distance. Default is `mode=min` despite the core method docstring saying mean. No pixel/nm rescaling in that path. Outputs go beside input unless outputFolder is supplied. |
| Matrix barcode sidecar | `_uniqueBarcodes.ecsv` is plain numeric text written with np.savetxt, not an Astropy ECSV document. Preserve this distinction for np.loadtxt consumers. |
| `trace_physical_vs_genomic_distance` | CSV/plot labels say nm, but the pipeline does not rescale loaded coordinates. Distances and filtering remain in input units. Genomic separations are interval-midpoint kbp; no chromosome stratification. Duplicate observations overwrite previous coordinates. |
| `plot_bootstrapping` | SciPy bootstrap uses np.mean, then reports median of bootstrap means and their SD. Filenames/titles say median. Missing pairs need at least two observations. Actual filenames include underscore-separated mode and color limits. |
| `plot_4m` | Anchor-to-partner pairwise proximity, not four-way contact. Resampled IDs are subset with isin, losing duplicate draws. SEM is std of cycle frequencies divided by sqrt(cycles). This uncertainty procedure should be reviewed before scientific use. |
| `trace_3way_coloc` | Trace bootstrap preserves repeated draws; SEM still divides bootstrap SD by sqrt(cycles). Each partner can be near a different duplicate anchor observation; no partner-partner proximity requirement. Output prefix includes input stem. |
| `plot_3way_coloc` | A label map is effectively needed despite being parser-optional. Title cutoff is hardcoded to 0.2 µm. Frequency NPY alone does not preserve barcode mapping or original cutoff. |
| Matrix comparisons | Require matched barcode ordering. `plot_matrix_comparison` loads NPY despite “trace file” help; violin values are log-transformed. Both rank-sum implementations use scipy.stats.ranksums, not paired signed-rank tests. |
| Shuffle conventions | plot_him_matrix resolves barcode identities, despite index-oriented help. Other matrix plotting paths use zero-based indices. |
| `plot_compare2matrices` | Proximity uses the total-trace denominator even when matrix_norm_mode names another mode. Difference plot maps normalize=none to maximum normalization, while mixed plot does not normalize. Ratio is log2. Wilcoxon-named NPY stores log10(p), with no multiple-testing correction. |
| `trace_merge` | --pipe defaults true with no disabling option, bypassing advertised folder input. Use ≥2 stdin paths. Conflicting ROI values are remapped; append does not regenerate Trace_ID/Spot_ID. |
| `localization_merge` | Reads filenames from stdin without a pipe flag. One-file input reaches uninitialized variables; use ≥2 files. |
| `trace_assign_mask` | Labels rows using mask[x,y], not entire traces. Explicit --pixel_size remains a string and can fail numeric division; default 0.1 is numeric. Input coordinates need valid mask indices. Label strings containing x are reset. |
| Label partition/selection | Substring matching per row; a Trace_ID can appear in both partitions. Negative split filename contains a colon, problematic on Windows. The trace_remove_label docstring reverses its actual selection behavior. |
| FOF-CT export | BED and metadata file are needed at runtime despite parser-optional flags. Default metadata file is cwd/parameters.json. Default CSV is written in cwd; overwrite protection uses bitwise ~ on a boolean and does not prevent overwriting. |
| FOF-CT roundtrip | Export omits Mask_id and label; import supplies placeholders. This does not recover original annotations. |
| Older output stem handling | Several commands use split('.') rather than a final-extension operation; dots in parent paths may truncate output prefixes. |
| Documentation aliases | Some examples use singular --anchor; exact parser flag is --anchors. Default argparse abbreviation may accept singular, but catalog examples deliberately use exact flags. |

These are observations from source review and smoke execution, not changes to
the analytical implementations. Numerical estimators and behavior remain intact.

## Incomplete documentation and scientific uncertainty

`trace_split_labels` has no dedicated preexisting command documentation page or
command-specific test; the catalog describes its actual row-level behavior.
Several other commands have no dedicated test module, shown in `INVENTORY.md`.
Source-documented Python methods sometimes promise tables in Returns sections
while mutating self.data and returning None. These conflicts are flagged and literal
documented return text is kept separate from annotations and implementation notes.

Useful names without public callable documentation include
`calculate_contact_probability_matrix`, `calculate_nan_matrix`, `bootstraps_matrix`,
`Wilcoxon_matrix`, and several plotting/reduction helpers. Their CLI behavior was
reviewed where relevant, but they are not promoted as documented Python APIs.
The inclusion ledger records the remaining gaps rather than claiming stability.
The short documentation of `distances_2_coordinates` does not establish accuracy,
missing-value support or positive eigenvalues; this entry is explicitly unexecuted.
Matrix-derived barcode detection counts require a finite pair and can miss isolated
singletons, unlike table-derived presence. Missing-pair radius-of-gyration behavior
also deserves scientific review before interpretation as a conventional complete-
trace Rg.

Beam resolution confidence is a normalized score gap, not a calibrated probability;
its multiplicity likelihood posterior is separate. Dataset representativeness,
detection heterogeneity and nuisance-rate assumptions affect interpretation.
These assumptions are described and exposed through diagnostics, not independently
validated here. Synthetic smoke tests cannot establish calibration or biological
accuracy. Other scientific uncertainties are attached directly to tool entries.

## Overlap and recommended improvements (not implemented)

- Keep `trace_stats` for global distinct counts and `trace_analyzer` for within-trace distributions. Add a brief parser epilog linking these tasks explicitly.
- Publish a parser/tool registry with structured metadata and consistent input/output path semantics. Require runtime-essential inputs in argparse, support exact flags consistently, and retain explicit deprecated aliases.
- Separate analysis estimators from names: clarify bootstrap statistic/error definitions, R(s) units, 4M dimensionality, rank-sum versus signed-rank tests, three-way detection denominators and duplicated-anchor semantics. Review estimators separately with biological validation; do not silently change them through documentation work.
- Harmonize duplicate-localization policy across distance generation, Pearson comparisons, scaling and co-localization, or expose and record the selected policy.
- `plot_him_matrix`, automatic trace_to_matrix plots and plot_compare2matrices all reduce PWD stacks but differ in defaults. Record normalization/threshold/units in output metadata. A future agent or MCP adapter should use catalog IDs and explicit contracts, not infer equivalence from output labels.
- `plot_matrix_comparison` and `trace_pearsons` both report Pearson correlation, but differ in inputs, barcode matching and duplicate handling. `plot_compare2matrices` adds per-pair tests and different normalization. Keep them distinct in task routing.
- Align `plot_4m`, `trace_3way_coloc` and `plot_3way_coloc` naming and output metadata; add cutoff/barcode sidecars and valid standalone replot defaults.
- Fix documented path/overwrite/typing issues in a separate implementation task. Add noninteractive missing-data policies, structured errors and result manifests; a zero exit code currently does not always establish successful analysis.
- Document a supported Python API boundary with schemas, return types and mutability, plus tiny fixtures for each public operation. Retain direct parser validation in CI, run example smoke tests on dependency updates, and review scientific text whenever source fingerprints change.

The catalog is portable YAML 1.2/JSON and Markdown; it has no PantheonOS,
vector database, embedding, API-service, orchestration or MCP dependency.
A later MCP adapter can reuse IDs, parameter snapshots, task distinctions and
verified argv examples without changing traceratops algorithms.
