# trace_to_matrix

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_to_matrix.parse_arguments
   :prog: trace_to_matrix
```

## Description

`trace_to_matrix` converts a tracefile to a NUMPY array containing the pair-wise distances between barcodes, for each trace in the tracefile. The dimensions of the NUMPY matrix is therefore: N\*N\*N_traces, where N is the number of barcodes and N_traces the number of traces.


## Usage

```bash
trace_to_matrix --input traces.ecsv
```

## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_Matrix_PWDscMatrix.npy`: a NUMPY matrix file described above,

2. `[tracefile]_Matrix_uniqueBarcodes.ecsv`: the unique barcode list in ecsv format,

3. `[tracefile]_Matrix_Nmatrix`: a matrix of the number of times each combination of barcodes is detected (N-matrix),
[](../../_static/merged_traces_filtered_split_Matrix_Nmatrix.png)

4. `[tracefile]_Matrix_HiMmatrix`: the Hi-M contact matrix calculated using a fixed threshold (default=0.25 microns),
[](../../_static/merged_traces_filtered_split_Matrix_HiMmatrix.png)

5. `[tracefile]_Matrix_PWDmatrixKDE`:the PWD KDE plot.
[](../../_static/merged_traces_filtered_split_Matrix_PWDmatrixKDE.png)

6. `[tracefile]_Matrix_PWDmatrixMedian`: the PWD median plot
[](../../_static/merged_traces_filtered_split_Matrix_PWDmatrixMedian.png)


### additional outputs

`trace_to_matrix` can also produce a figure with the distributions of distances for each combination of barcodes. This figure is not produced by default as it is often the slowest step for the processing of datasets with many barcodes.

To also calculate and save `*_Matrix_PWDhistograms.png`, opt in explicitly:

```bash
trace_to_matrix --input traces.ecsv --plot_histograms
```

For larger datasets, you can combine histogram generation with per-trace parallel matrix construction:

```bash
trace_to_matrix --input traces.ecsv --n_jobs 6 --plot_histograms
```
