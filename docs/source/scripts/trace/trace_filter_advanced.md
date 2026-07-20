# trace_filter_advanced

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_filter_advanced.parse_arguments
   :prog: trace_filter_advanced
```

## Output Files

- A filtered trace table saved in the output folder selected with `--output` (default folder/tag: `filtered`).
- Diagnostic plots such as trace statistics, pairwise-distance statistics, and duplicated-barcode statistics, using `--output_format`.

## Examples

```bash
trace_filter_advanced --help
```

## Notes

- See the command-line reference above for the complete option list.
