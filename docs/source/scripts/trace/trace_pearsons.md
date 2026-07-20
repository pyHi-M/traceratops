# trace_pearsons

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_pearsons.parse_arguments
   :prog: trace_pearsons
```
## Description

This script is intended to compare the results from different experiments. These can be different fields of view of a single acquisition or different biological replicates. The script calculates, for each tracefile, the median pairwise distance map. Then it calculates the bin-by-bin Pearson correlation between all tracefiles provided. The result is therefore a matrix of Pearson correlations.

## Usage

```bash
$ ls *ecsv | trace_pearsons --pipe
```


## Examples
```bash
$ ls *ecsv | trace_pearsons [options]
$ find . -name "*.ecsv" | trace_pearsons [options]
```


## Output Files

1. `trace_correlation_matrix.png`: A Pearson correlation matrix plot comparing all input trace tables.

![trace_correlation_matrix](https://github.com/user-attachments/assets/ae4c19f7-3638-4e56-9901-5a206d0e64d6)

2. `trace_correlation_matrix.npy`: Values of the correlation matrix.


## Notes

- Input files must be in ECSV format compatible with ChromatinTraceTable
- The script identifies unique parts of filenames to create readable labels in the plot
- Correlation is calculated based on the spatial distances between barcode pairs
