# trace_filter_advanced

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_filter_advanced.parse_arguments
   :prog: trace_filter_advanced
```

## Output Files

For each trace file analyzed, the script generates:

1. `[output]/[tracefile]_[output].ecsv`: Filtered trace table saved in the folder selected with `--output` (default folder/tag: `filtered`).

2. `[output]/trace_stat_[tag].[format]`: Trace-statistics diagnostic plot written during advanced filtering.

3. `[output]/pairwise_distance_stat.[format]`: Pairwise-distance threshold diagnostic plot written when pairwise-distance statistics are saved.

4. `[output]/duplicated_bc_pwd_stat.[format]`: Diagnostic plot for duplicated-barcode pairwise-distance statistics.

## Examples

```bash
trace_filter_advanced --help
```

## Notes

- See the command-line reference above for the complete option list.
