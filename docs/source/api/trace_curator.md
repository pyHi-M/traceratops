# Trace curation Python API

`trace_curator` operates on a resolved single polymer after reconstruction or
trace splitting. It asks whether removing a localization explains several
spatially inconsistent relationships. Duplicate barcode assignments must be
resolved before calling this API.

```python
from traceratops.core.trace_curator import fit_curator_model, curate_trace

# A population of resolved traces, most expected to be valid.
# Each trace contains explicit Genomic_Position values.
model = fit_curator_model(traces)
result = curate_trace(trace, model)
print(result.terminal_state)
print(result.removed_localizations)
print(result.iterations)
```

## Input contract

Single traces can be pandas DataFrames or Astropy Tables. Required columns are
`x`, `y`, `z`, `Spot_ID`, and `Barcode #` (the ECSV convention) or `Barcode`.
Coordinates must be finite numeric values in consistent spatial units.
Identities must be nonmissing; `Spot_ID` and barcode identities must each be
unique within the trace. Arbitrary string barcodes are supported. If both
barcode columns are present, they must agree. `Trace_ID` is optional for a single
trace, but when present must identify exactly one trace.

Genomic positions must be supplied through `Genomic_Position` or the explicit
`genomic_positions={barcode: position, ...}` fitting argument. The model retains
this mapping so subsequent traces can omit the position column for known
barcodes. Supplied columns must agree with the fitted mapping. Positions must
be finite, unique by barcode, and use one genomic coordinate system and unit.
If `Chrom` is present, it must identify one chromosome and agree across fitting,
calibration, and curation. Chromosome-free inputs are the caller's assertion of
one coordinate system. Numeric barcode identities and row order never supply
positions. BED midpoint conversion is the caller's responsibility.

Fitting accepts a combined DataFrame/Table with `Trace_ID`, or a reiterable
collection of single traces with a stable iteration order. Default cross-fitting
requires reiterable input and rejects one-pass iterators. A one-pass reference
iterator is supported only with explicit separate calibration traces. Empty and
fewer-than-three-localization curation inputs return `unscoreable`, provided their
schema is valid. Fit fails clearly if too few
calibration traces can be scored to estimate the requested empirical tail.

## Public functions

```python
fit_curator_model(
    traces, *, genomic_positions=None, trace_fpr=0.01,
    minimum_reference_observations=20, reference_mode="separation",
    calibration_traces=None, n_crossfit_folds=5,
) -> TraceCuratorModel

curate_trace(trace, model, *, max_removals=3) -> TraceCurationResult
```

All input traces must be resolved single polymers. The fitting population is
expected to contain a majority of valid chromatin traces; individual aberrant
traces do not need to be identified beforehand.
The default reference mode pools spatial distances by **exact absolute genomic
separation**. No separation binning, rounding, interpolation, or sparse-bin
fallback is performed. Optional `barcode_pair` mode uses unordered barcode
identity pairs. Unsupported pairs are ignored, including references below the
minimum observation count. Missing relationships are not assigned zero scores.

## Two fitting modes

The default call:

```python
model = fit_curator_model(traces)
```

performs deterministic K-fold cross-fitting, with `n_crossfit_folds=5` by
default. Assign trace number `i` to fold `i % n_crossfit_folds`. Combined tables
use first-appearance `Trace_ID` order; individual-trace collections use collection
order, even without `Trace_ID`. There is no shuffle, random seed, or hashing.
Fold sizes differ by at most one. The fold count must be an integer between two
and the number of available traces; it is not automatically reduced.

For each held-out fold, fit an empirical reference on the other folds and score
each held-out trace against it. Each trace is held out exactly once, and no
trace contributes pair distances to the reference used for its own calibration
score. Pool only finite held-out `C_top3` values to derive the threshold. Discard
each fold reference before fitting the next. After calibration, fit the final
production reference once on **all** input traces. No in-sample calibration mode
is used silently.

Cross-fitting prevents self-inclusion but does not make the reference estimator
robust to arbitrarily high contamination. The prototype assumes most traces
are valid. **Contamination robustness has not yet been quantified**, and no
supported contamination percentage is claimed. The nominal FPR specifies an
empirical calibration tail, not a guaranteed error rate on unseen traces.

When an independent high-quality reference/control population is available:

```python
model = fit_curator_model(traces, calibration_traces=clean_control_traces)
```

fit the production reference from `traces`, score the supplied independent
controls against it, and derive the threshold from those control scores. This
bypasses cross-fitting; `n_crossfit_folds` is unused and recorded as `None`.
Independence of explicitly supplied controls is the caller's responsibility.

In **both modes**, require `0 < trace_fpr < 1` and at least
`ceil(1 / trace_fpr)` finite calibration scores. At the default `0.01`, fewer
than **100** scoreable traces raises a clear `ValueError`, even if the total
input count is large. Sparse reference support or too-small traces can reduce
the scoreable sample. If `N * trace_fpr < 5`, emit a `UserWarning` via
`warnings.warn`: the tail is estimable but may be unstable because there are
fewer than five expected tail observations. At 1%, this warning applies to
100–499 finite scores. This is a statistical criterion, not a biological
contamination claim. No warning is emitted at or above five expected tail
observations.

Only finite costs enter the cutoff; unscoreable calibration traces are counted
separately. Pair observations use compact `array("d")` double buffers, converted
once per reference key to sorted NumPy arrays; scoring uses binary searches
without re-sorting. Buffers are released as keys are converted. No collection of
fold-specific references, trace-score histories, or pairwise calibration tables
is retained. Scalar calibration costs also use a compact double buffer. Final
reference arrays retain their immutable backing.

## Scientific decision rule

For a distance `d` and sorted reference sample `s` of length `n`:

```python
lower = (searchsorted(s, d, side="right") + 1) / (n + 1)
upper = (n - searchsorted(s, d, side="left") + 1) / (n + 1)
anomaly = -log10(min(1, 2 * min(lower, upper)))
```

`C_top3` is the sum of the three largest finite pair anomalies (or all available
scores if fewer than three). No finite pair gives an unavailable cost, never
zero. Traces smaller than three localizations are unscoreable.

The calibration cutoff uses `allowed = floor(trace_fpr * n_calibration)` and
`sorted_costs[max(0, n_calibration - allowed - 1)]`, with the finite calibration
sample count. Abnormality requires strictly
`C_top3 > global_threshold`; equality is normal.

For an abnormal trace, each finite `LOO_top3` is the current cost minus the cost
without that localization. If removing a localization leaves no scoreable
relationship, its LOO is unavailable. A unique maximum selects the candidate.
An exact top3 tie invokes `LOO_top5` **only among those tied candidates**. A
unique top5 maximum resolves the tie; a remaining tie stops unresolved. There
is no tolerance-based tie or input-order tie breaking, and top5 never replaces
the global top3 statistic.

Every accepted removal triggers rescoring and the global test again. At least
three localizations are preserved. The default removal cap is three; a
configurable nonnegative `max_removals` allows inspection with zero removals.
There is no additional candidate confidence or actionability threshold.

## Returned objects

`TraceCuratorModel` is frozen and contains:

- `reference_distributions`: immutable mapping to sorted, read-only NumPy arrays;
- `genomic_positions`: immutable barcode-to-position mapping;
- `global_threshold`, `reference_mode`, `minimum_reference_observations`, `trace_fpr`;
- `n_reference_traces`, `n_reference_pairs`, `n_calibration_traces` (scoreable),
  `n_calibration_unscoreable`, `calibration_source` (`"crossfit"` or
  `"separate_traces"`), `n_crossfit_folds` (or `None`), and optional `chromosome`;
- `n_reference_keys`, `n_supported_reference_keys`,
  `median_reference_observations`, `min_reference_observations_actual`, and
  `fraction_reference_keys_below_minimum`.

Support diagnostics describe the **final full reference**, counting distinct
keys equally. A supported key has at least `minimum_reference_observations`;
median and minimum include all keys, including unsupported ones. These compact
summaries do not retain an additional per-key table. With an empty reference,
counts and the minimum are zero, and median/fraction are NaN.

`TraceCurationResult` contains independent pandas copies `curated_trace` and
`removed_localizations`, preserving input columns, original indices, and stable
`Spot_ID` values. Retained rows keep input order; removed rows follow removal
order. The frozen result also contains `iterations`, `initial_cost`,
`final_cost`, `global_threshold`, `n_removed`, `terminal_state`, and `stop_reason`.
Its DataFrames are mutable caller-owned outputs; immutable diagnostics are a
tuple of frozen `TraceCurationIteration` objects.

Terminal states are `normal` (unchanged and below/equal cutoff), `curated`
(removals and terminal normal), `abnormal_unresolved`, or `unscoreable`
(unavailable initial cost). A previously abnormal trace that loses reference
support after curation remains unresolved. Costs unavailable for scoring use
NaN. Stop reasons are `global_cost_normal`, `too_few_localizations`,
`insufficient_reference`, `minimum_trace_size_reached`, `max_removals_reached`,
`no_scoreable_candidate`, or `ambiguous_candidate`.

Each iteration, including the final stopping state, records zero-based
`iteration`, `n_localizations_before`, `C_top3_before`, `global_threshold`,
`globally_abnormal`, `n_scoreable_pairs`, candidate `candidate_spot_id` and
`candidate_barcode`, `loo_top3`, `n_tied_top3`, `top5_tiebreak_used`, `loo_top5`,
`n_tied_top5`, `removal_accepted`, and `stop_reason`. When attribution is attempted,
LOO fields record maximal scores even if a tie prevents selection. Unused
scores are NaN and absent candidate identities are None. Full pairwise tables
are not included.

This prototype provides the reusable Python API only. CLI loading, batch ECSV
output, model persistence, and parallel execution are future work.
