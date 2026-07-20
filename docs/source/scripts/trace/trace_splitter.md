# trace_splitter

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_splitter.parse_arguments
   :prog: trace_splitter
```

## Description

This script will calculate the radii of gyration of all traces in a tracefile, estimate the median and standard deviation of the distribution to identify outliers which may arise from the merging of two traces. It will then target these outliers and perform a K-means clustering to break the trace. Each time a trace is split, it will be replaced in the tracefile and new identifiers will be provided.

Traces with R_g higher than the `median` + `std_threshold` will be treated as outliers. `std_threshold` is therefore an input argument.

Clusterization will be performed assuming 2 traces by default. This number can be modified using `num_clusters` as input argument.

## Usage

```
trace_splitter --input path/to/your/trace_file.ecsv
```
## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_split.ecsv`: The main output is a new trace file containing the original un-split traces + the split traces.


## Example

```
trace_splitter --input original_traces.ecsv --std_threshold 1.5 --num_clusters 3
```

Given a chromatin trace table, this script:
   - Computes radius of gyration (Rg) for all traces.
   - Identifies traces with Rg larger than `mean + N * std_dev`.
   - Uses K-means to split these traces into `num_clusters`.
   - Saves the modified trace table with updated Trace_IDs.
