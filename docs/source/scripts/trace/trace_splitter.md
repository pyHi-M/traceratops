# trace_splitter

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_splitter.parse_arguments
   :prog: trace_splitter
```

## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_split.ecsv`: Trace table containing the original un-split traces plus the split traces, with updated `Trace_ID` values. When `--output` is provided for a single input, that exact path is used instead.

## Examples

```console
trace_splitter --input original_traces.ecsv --std_threshold 1.5 --num_clusters 3
trace_splitter --input original_traces.ecsv --split-all
trace_splitter --input original_traces.ecsv --clustering-method hdbscan --min-cluster-size 3 --min-samples 3
```

Given a chromatin trace table, this script:
   - Computes radius of gyration (Rg) for all traces.
   - By default, identifies traces with Rg larger than `mean + N * std_dev`.
     `--split-all` instead sends every trace to the selected clustering method.
   - Uses K-means by default, preserving the behavior of earlier releases. Select
     HDBSCAN with `--clustering-method hdbscan`.
   - Saves the modified trace table with updated Trace_IDs.

## Clustering options

K-means uses `--num_clusters` (default: 2). Its deterministic random seed and
initialization settings remain unchanged.

HDBSCAN provides these controls, based on the corresponding scikit-learn
parameters:

- `--min-cluster-size` (default: 3): minimum number of detections in a cluster.
- `--min-samples` (default: 3): number of neighboring samples required for a
  point to be considered a core point. Larger values make clustering more
  conservative.
- `--cluster-selection-method` (default: `eom`): use excess-of-mass (`eom`) or
  select leaf clusters (`leaf`).
- `--cluster-selection-epsilon`: distance below which clusters are merged. If
  omitted, trace_splitter computes an Otsu threshold from the pairwise distances
  in each trace and uses their median, following the automatic thresholding
  strategy of the prototype script.
- `--allow-single-cluster`: permit HDBSCAN to select one cluster. It is disabled
  by default.

HDBSCAN labels some detections as noise. Noise detections are retained under
their original `Trace_ID`; trace_splitter never discards input rows. Traces with
fewer detections than `--min-cluster-size` are left unchanged.

Coordinates read from ECSV are converted to a contiguous 64-bit floating-point
array before clustering. This is intentional: scikit-learn's compiled HDBSCAN
tree traversal can fail for 32-bit ECSV coordinates when a non-zero selection
epsilon is used together with `--allow-single-cluster`. If a compatible input
still encounters scikit-learn's known scalar-conversion failure in its epsilon
tree traversal, trace_splitter reports a warning and retries that trace with
epsilon-based merging disabled. Other HDBSCAN errors are not suppressed.

## Notes

- See the command-line reference above for the complete option list.
- Without `--split-all` and `--clustering-method`, selection by radius of
  gyration and K-means clustering behave as in previous releases.
