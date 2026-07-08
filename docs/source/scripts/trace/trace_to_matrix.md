# trace_to_matrix

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_to_matrix.parse_arguments
   :prog: trace_to_matrix
```


## Notes on default outputs

By default, `trace_to_matrix` now skips the all-pairs PWD KDE histogram figure because this calculation is often the slowest step for datasets with many barcodes. The default run still produces the single-cell PWD matrix, unique barcode list, N-matrix, Hi-M contact matrix, PWD median plot, PWD KDE plot, and N-matrix plot.

To also calculate and save `*_Matrix_PWDhistograms.png`, opt in explicitly:

```bash
trace_to_matrix --input traces.ecsv --plot_histograms
```

For larger datasets, you can combine histogram generation with per-trace parallel matrix construction:

```bash
trace_to_matrix --input traces.ecsv --n_jobs 6 --plot_histograms
```
