# trace_stats

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_stats.parse_arguments
   :prog: trace_stats
```

## Output Files

The output is on the terminal. An example is provided below:

```bash
Statistics for merged_traces.ecsv:
- Number of unique ROIs: 3
- Number of unique chromatin traces: 6115
- Number of unique barcodes: 86
```

## Examples

```trace_stats.py --input traces_KC_AB_merged.ecsv```

Terminal output:
```
$ Importing table from pyHiM format
Successfully loaded trace table: traces_KC_AB_merged.ecsv
Statistics for traces_KC_AB_merged.ecsv:
- Number of unique ROIs: 26
- Number of unique chromatin traces: 7573
```

## Notes

- See the command-line reference above for the complete option list.
