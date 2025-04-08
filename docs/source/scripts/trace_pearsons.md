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

## Notes
- Input files must be in ECSV format compatible with ChromatinTraceTable
- The script identifies unique parts of filenames to create readable labels in the plot
- Correlation is calculated based on the spatial distances between barcode pairs
