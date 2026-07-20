# trace_physical_vs_genomic_distance

**Reliability status**: `development`

`trace_physical_vs_genomic_distance.py` summarizes how physical chromatin-trace distances change with genomic separation. It reads a chromatin trace table, computes pairwise distances between all barcode loci in each trace, bins those distances by genomic distance, and writes a CSV table that can also be plotted as a log-log physical-vs-genomic-distance curve.

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_physical_vs_genomic_distance.parse_arguments
   :prog: trace_physical_vs_genomic_distance
```

## Description

1. The script loads an input trace table supported by `ChromatinTraceTable` (`.ecsv`, `.dat`, `.4dn`, or `.csv`).
2. Each trace is converted into a NumPy array with one row per trace, one column per barcode, and X/Y/Z coordinates in the final dimension.
3. Genomic separation is computed from barcode## Description
 midpoint positions, using `(Chrom_Start + Chrom_End) / 2`, and is reported in kilobase pairs (kbp).
4. Physical pairwise distances are computed for:
   - full 3D Euclidean distance, unless `--no_3d` is used;
   - projected X-axis distance;
   - projected Y-axis distance;
   - projected Z-axis distance.
5. Physical distances above `--dist_threshold` are replaced by `NaN` before summarization.
6. Genomic distances are split into `--gen_dist_bins` equally spaced bins, and the median physical distance is reported for each bin and axis.

## Output Files

For each trace file analyzed, the script generates:

1. A CSV file containing one row per genomic-distance bin and axis:

- `genomic distance (kbp)`: midpoint of the genomic-distance bin.
- `log10 genomic dist (kbp)`: log10-transformed bin midpoint.
- `median euclidean distance (nm)`: median physical distance for barcode pairs in the bin.
- `log10 median dist (nm)`: log10-transformed median physical distance.
- `n_data`: number of non-NaN physical distances used in the bin.
- `experiment`: experiment label from `--experiment`, or `exp_None` when no label is supplied.
- `axis`: one of `3D`, `X`, `Y`, or `Z`.

If `--interloci_output` is supplied, the script also writes a CSV containing genomic distances between consecutive barcode loci.


2. `[tracefile]_physical_vs_genomic_plot.png`: a figure with 4 panels with the physical versus genomic distance plot, with a power law fit and confidence intervals.

[](../../_static/merged_traces_filtered_physical_vs_genomic_plot.png)


### Plot and power-law fit

When `--plot` is provided, the script saves a stacked log-log plot with one panel per calculated axis. The x-axis is shared and shown only on the bottom panel, while a single shared y-axis label spans the panels to avoid label overlap.

For every axis panel and experiment, the script fits the binned median distances to a power-law model:

```text
physical_distance = a * genomic_distance^b
```

The fit is performed by linear regression in log10 space. The plotted fit line is shown together with a 95% confidence interval for the fitted mean response. Each panel legend reports both fitted coefficients:

- `a`: the scale coefficient.
- `b`: the power-law exponent.

## Example

```bash
trace_physical_vs_genomic_distance.py \
  --input traces_KC_AB_merged.ecsv \
  --output binned_physical_vs_genomic_distances.csv \
  --interloci_output interloci_genomic_distances.csv \
  --plot _physical_vs_genomic_plot.png \
  --gen_dist_bins 50 \
  --dist_threshold 1000 \
  --experiment KC_AB
```

The plot filename is appended to the input stem by default. For the example above, the figure is saved as `traces_KC_AB_merged_physical_vs_genomic_plot.png`.
