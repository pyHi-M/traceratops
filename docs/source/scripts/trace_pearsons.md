# trace_pearsons

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_pearsons.parse_arguments
   :prog: trace_pearsons
```

## Examples
```bash
$ ls *ecsv | trace_pearsons [options]
$ find . -name "*.ecsv" | trace_pearsons [options]
```


## Output
- A Pearson correlation matrix plot comparing all input trace tables
- Terminal output showing the numerical correlation values


![trace_correlation_matrix](https://github.com/user-attachments/assets/ae4c19f7-3638-4e56-9901-5a206d0e64d6)

## Notes
- Input files must be in ECSV format compatible with ChromatinTraceTable
- The script identifies unique parts of filenames to create readable labels in the plot
- Correlation is calculated based on the spatial distances between barcode pairs
