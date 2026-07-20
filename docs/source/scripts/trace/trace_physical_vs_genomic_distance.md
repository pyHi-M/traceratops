# trace_physical_vs_genomic_distance

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_physical_vs_genomic_distance.parse_arguments
   :prog: trace_physical_vs_genomic_distance
```


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
