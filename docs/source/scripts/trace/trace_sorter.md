# trace_sorter

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_sorter.parse_arguments
   :prog: trace_sorter
```

## Example

```
trace_sorter --input original_traces.ecsv --sort_by number_of_spots
```

Given a chromatin trace table, this script:
   - Sorts complete traces by either `number_of_spots` or `radius_of_gyration`.
   - Sorts in descending order by default; `--ascending` reverses that order.
   - Keeps every trace's rows together and preserves their internal order.
   - Saves to the path supplied with `--output`, or appends `_sorted` to the input basename.
