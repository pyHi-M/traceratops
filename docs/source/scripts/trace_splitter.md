# trace_splitter

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_splitter.parse_arguments
   :prog: trace_splitter
```

## Example

```
trace_splitter --input original_traces.ecsv --std_threshold 1.5 --num_clusters 3
```

Given a chromatin trace table, this script:
   - Computes radius of gyration (Rg) for all traces.
   - Identifies traces with Rg larger than `mean + N * std_dev`.
   - Uses K-means to split these traces into `num_clusters`.
   - Saves the modified trace table with updated Trace_IDs.
