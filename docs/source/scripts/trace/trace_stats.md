# trace_stats

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_stats.parse_arguments
   :prog: trace_stats
```
## Description

This script provides general statistics on a trace file.

## Usage

```
trace_stats.py --input path/to/your/trace_file.ecsv
```

## Output Files

The output is on the terminal. An example is provided below:

```bash
Statistics for merged_traces.ecsv:
- Number of unique ROIs: 3
- Number of unique chromatin traces: 6115
- Number of unique barcodes: 86
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
