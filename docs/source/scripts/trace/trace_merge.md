# trace_merge

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_merge.parse_arguments
   :prog: trace_merge
```

## Description
This script merges two or more trace files.


## Usage

```
ls *ecsv | trace_merge
```

## Output Files

The script produces as output the merged trace file:

1. `merged_traces.ecsv`: this is the default name for the output tracefile if non is provided as argument.
