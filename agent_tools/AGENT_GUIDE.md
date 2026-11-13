# traceratops: agent tool selection

traceratops post-processes chromatin tracing tables: quality control, curation,
polymer resolution, matrices, statistics and visualization. It does not localize
spots from raw images. This guide targets the `dev` snapshot recorded in
`catalog.yaml`; validate against the installed package before use.

## Lookup

| Scientific request | Select |
|---|---|
| Count dataset-wide distinct traces, ROIs, barcode identities | `trace_stats` |
| **Plot distribution of unique barcodes per trace**; trace completeness QC | **`trace_analyzer`** |
| Plot localizations per trace, barcode multiplicity, detection efficiency, neighbor distances, centered projections | `trace_analyzer` |
| SNR/shape QC, localization detections per barcode | `localization_analyzer` |
| Filter coordinate bounds, minimum distinct barcode count, quality metrics, labels, duplicate spots | `trace_filter` |
| Append trace tables / append localization tables | `trace_merge` / `localization_merge` |
| Collect one file per ROI subdirectory | `collect_files` |
| Resolve duplicated/merged polymers; spatial splitting | `trace_splitter` (choose method explicitly) |
| Sort by spot count or radius of gyration | `trace_sorter` |
| Label trace rows with a 2D mask / partition label-positive and negative rows | `trace_assign_mask` / `trace_split_labels` |
| Add genomic intervals or substitute barcode IDs from BED | `trace_genomic_coordinates` |
| Physical versus genomic distance, R(s), power-law fit | `trace_physical_vs_genomic_distance` |
| Export traces as PDB pseudoatoms for a viewer | `trace_plot` |
| Generate per-trace 3D pairwise-distance stack and QC matrices | `trace_to_matrix` |
| Plot proximity frequency, ensemble median/KDE distance, missing-data map | `plot_him_matrix` |
| Bootstrap pairwise **mean** distances, median of bootstrap estimates and bootstrap SD | `plot_bootstrapping` |
| Pearson comparisons across trace datasets, matching barcode identities | `trace_pearsons` |
| Two matrix ensembles: scatter/Pearson, log violin distributions | `plot_matrix_comparison` |
| Matrix difference/log2 ratio, mixed triangles, rank-sum comparison | `plot_compare2matrices` |
| Pairwise anchor proximity profile (4M style) | `plot_4m` |
| Anchor plus two partner loci, three-way proximity | `trace_3way_coloc` |
| Replot existing three-way NPY frequency matrix | `plot_3way_coloc` |
| ECSV → FOF-CT / FOF-CT → ECSV | `trace_export_to_fofct` / `trace_import_from_fofct` |

Deprecated names: `figure_him_matrix` → `plot_him_matrix`;
`compare_pwd_matrices` → `plot_matrix_comparison`; `plot_3way` →
`trace_3way_coloc`; `replot_3way` → `plot_3way_coloc`.

## Retrieve exact syntax

Read only the needed entry, then its referenced example:

```python
import json
from pathlib import Path
base = Path("agent_tools")  # repository or deployed catalog directory
catalog = json.loads((base / "catalog.yaml").read_text())
tool = next(t for t in catalog["cli_tools"] if t["id"] == "cli.trace_analyzer")
print(json.dumps(tool, indent=2))
```

Use `tasks.yaml` for synonyms, scientific distinctions and ordered workflows.
Python-only analytical capabilities are in `catalog["python_apis"]`, including
PWD-based radius of gyration, coordinate inference, matrix-based barcode detection,
localization comparison and PyMOL helpers. A source-documented name is not a
stability guarantee; inspect its input contract and caveats.

## Execute existing implementations

1. Match the scientific quantity and observational unit before selecting a tool.
2. Read `runtime_requirements`, referenced input schema, parser, outputs and limitations. An argparse-optional input may still be required downstream. Parser defaults can differ from runtime fallbacks.
3. Prefer exact canonical command names and flags; do not invent Python APIs or flags. `--help` can confirm installed syntax; it cannot establish scientific meaning.
4. Use copied inputs and an isolated working directory. Some outputs follow input paths; an output option may be a tag rather than a path. Save stdout/stderr.
5. Run through the Python interpreter tool using subprocess argv, with `shell=False`:

```python
import subprocess
result = subprocess.run(
    ["trace_analyzer", "--input", "trace.ecsv", "--output_format", "png"],
    cwd="/absolute/analysis/directory", capture_output=True, text=True,
    timeout=300, check=False,
)
```

For filename-list input, set `input="a.ecsv\nb.ecsv\n"` and use the documented
pipe flag. `localization_merge` reads stdin without a flag; `trace_merge` also
needs at least two stdin paths. Confirm the executable is on PATH.

## Scientific distinctions and caveats

- Dataset-wide distinct barcodes ≠ distinct barcodes per trace ≠ localizations per trace. Repeated barcode detections increase localizations without adding identities. Completeness additionally requires an expected barcode universe.
- Distance stacks have length units. Proximity maps contain threshold-dependent fractions. `nonNANs` and `n_cells` normalization use different denominators.
- `trace_to_matrix` uses minimum duplicate-pair distances. `trace_pearsons` uses the last repeated barcode observation before ensemble medians.
- Filtering removes rows; spatial splitting reassigns traces; beam resolution may reject rows. Label splitting partitions rows, not necessarily whole traces.
- `plot_4m` is pairwise, not four-way. Three-way analysis requires both partners near an anchor; partners need not be near each other or the same duplicate anchor observation.
- Physical-versus-genomic CSV headings say nm, but values retain loaded coordinate units. Bootstrap distance filenames say median, but the resampled statistic is mean. See catalog before interpreting either output.
- Matrix comparisons require matched barcode ordering. “Wilcoxon” outputs use rank-sum tests; the NPY contains log10 p-values. Difference and mixed plots handle normalization differently. `plot_4m` bootstrap loses repeated draws. Review limitations before statistical interpretation.

## Custom code and reporting

Use custom Python for orchestration, input validation, format adaptation or a
quantity absent from both CLI operations and documented Python APIs. First inspect
the relevant tasks, canonical tools and Python entries; a failed run or missing
synonym does not prove missing functionality. State the precise gap and any
scientific assumptions before implementing a new estimator.

Report executable/version, catalog commit, exact argv and stdin paths, working
directory, inputs and units, parameters/normalization, exit code, output paths,
and checks performed. For errors report stdout/stderr and failed assumptions;
separate syntax, missing-data, dependency and scientific-compatibility failures.
A zero exit code alone is insufficient: verify expected artifacts and their meaning.
