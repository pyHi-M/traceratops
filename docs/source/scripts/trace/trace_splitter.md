# trace_splitter

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_splitter.parse_arguments
   :prog: trace_splitter
```

## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_split.ecsv`: Resolved trace table. Legacy methods retain all
   input rows and update `Trace_ID` values for split clusters. Beam mode can
   omit rejected detections or globally ambiguous traces. When `--output` is
   provided for a single input, that exact path is used instead.
2. `[tracefile]_split_diagnostics.ecsv`: Per-input-trace diagnostics produced
   by beam mode. Use `--diagnostics-output` to choose another path.

## Examples

```console
trace_splitter --input original_traces.ecsv --std_threshold 1.5 --num_clusters 3
trace_splitter --input original_traces.ecsv --split-all
trace_splitter --input original_traces.ecsv --method hdbscan --min-cluster-size 3 --min-samples 3
trace_splitter --input original_traces.ecsv --method beam --history-mode multi --distance-score residual
```

Given a chromatin trace table, this script:
   - Computes radius of gyration (Rg) for all traces.
   - By default, identifies traces with Rg larger than `mean + N * std_dev`.
     `--split-all` instead sends every trace to the selected clustering method.
   - Uses K-means by default, preserving the behavior of earlier releases. Select
     HDBSCAN or barcode-aware resolution with `--method hdbscan` or
     `--method beam`. The former `--clustering-method` spelling remains an alias.
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

HDBSCAN labels some detections as noise. In this legacy method, noise detections
are retained under their original `Trace_ID`; K-means and HDBSCAN do not discard
input rows. Traces with fewer detections than `--min-cluster-size` are unchanged.

## Barcode-aware beam resolution

Beam resolution uses repeated barcode identities to classify each input trace
as unchanged, a one-polymer duplicate-cleanup candidate, or a possible merged
two-polymer trace. The provisional classification thresholds are controlled by
`--split-min-repeated-barcodes`, `--split-min-repeated-fraction`, and
`--candidate-rule`. These settings are intended to be calibrated with
simulations; they are not authoritative biological thresholds. Radius of
gyration is reported as supporting information but is not required for beam
resolution.

Spatial compatibility is learned from input traces that have at most one
detection per barcode. Genomic midpoint separation is used when `Chrom`,
`Chrom_Start`, and `Chrom_End` are present; otherwise the resolver uses barcode
index separation. A separation with fewer than `--model-min-observations`
borrows the nearest sufficiently populated separation. If none exists, the
pooled clean-trace distance distribution is used. If the file has no clean
trace, a deliberately weak fallback is estimated from the nearest candidate
pair between detected barcode identities. `--variance-floor`
prevents degenerate distance samples from producing unstable scores.

History selection and distance cost are independent:

- `--history-mode nearest` scores against the last assigned localization.
- `--history-mode multi` averages over up to `--history-length` prior assigned
  localizations.
- `--distance-score residual` uses squared standardized distance residuals.
- `--distance-score likelihood` uses Gaussian negative log likelihood.

The defaults are `multi` and `residual`. Missing intermediate barcodes have no
separate penalty: their connection is scored at the corresponding larger
genomic separation. Each inferred polymer can contain at most one localization
for a barcode, and polymers need not have equal completeness.

`--rejection-cost` is the explicit algorithmic cost of leaving a localization
unassigned. Unassigned beam detections are omitted from the output. Surviving
`Spot_ID` values are copied exactly; only `Trace_ID` changes after a successful
split. This differs intentionally from legacy HDBSCAN, whose noise rows remain
under the original trace ID for backward compatibility.

For a two-polymer candidate, a normalized gap below `--minimum-confidence`
causes the complete input trace to be omitted with status
`removed_ambiguous`. For one-polymer cleanup, localizations at an unresolved
barcode are omitted while the rest of the trace is retained with status
`cleaned`.

Beam mode writes one diagnostics row per input trace to
`[tracefile]_split_diagnostics.ecsv`, or to `--diagnostics-output`. Diagnostics
include barcode multiplicity, Rg, inferred polymer count, best and alternative
scores, raw score gap, unassigned detections, and ambiguous barcode count.
`confidence_score` is the raw gap divided by the sum of the absolute best and
alternative scores. It is a non-probabilistic ranking diagnostic, **not** an
estimated probability.

Coordinates read from ECSV are converted to a contiguous 64-bit floating-point
array before clustering. This is intentional: scikit-learn's compiled HDBSCAN
tree traversal can fail for 32-bit ECSV coordinates when a non-zero selection
epsilon is used together with `--allow-single-cluster`. If a compatible input
still encounters scikit-learn's known scalar-conversion failure in its epsilon
tree traversal, trace_splitter reports a warning and retries that trace with
epsilon-based merging disabled. Other HDBSCAN errors are not suppressed.

## Notes

- See the command-line reference above for the complete option list.
- Without `--split-all` and `--method`, selection by radius of
  gyration and K-means clustering behave as in previous releases.
