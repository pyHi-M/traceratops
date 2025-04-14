# trace_stats

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_stats.parse_arguments
   :prog: trace_stats
```

## Example

```trace_stats.py --input traces_KC_AB_merged.ecsv```

Terminal output:
```
$ Importing table from pyHiM format
Successfully loaded trace table: traces_KC_AB_merged.ecsv
Statistics for traces_KC_AB_merged.ecsv:
- Number of unique ROIs: 26
- Number of unique chromatin traces: 7573
```
