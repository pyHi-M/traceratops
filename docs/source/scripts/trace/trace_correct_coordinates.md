# trace_correct_coordinates

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_correct_coordinates.parse_arguments
   :prog: trace_correct_coordinates
```

## Example usage

```
trace_correct_coordinates --input traces.ecsv --output traces_corrected.ecsv --max_iter 10 --tolerance 0.01
```

## Algorithm

1. Compute the initial center of mass (CoM) for each chromatin trace.
2. Iteratively adjust the z-coordinates of each barcode to minimize its distance to the CoM.
3. Update the trace table and save the corrected version.
