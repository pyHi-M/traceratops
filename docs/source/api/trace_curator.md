# Trace curation Python API

`trace_curator` operates on a resolved single polymer after reconstruction or
trace splitting. It asks whether removing a localization explains several
spatially inconsistent relationships. Duplicate barcode assignments must be
resolved before calling this API.

```python
from traceratops.core.trace_curator import fit_curator_model, curate_trace

# A collection of clean, resolved single-trace DataFrames or Astropy Tables.
# Each trace contains explicit Genomic_Position values.
model = fit_curator_model(reference_traces)
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
collection of single traces. A one-pass reference iterator requires separate
calibration traces. Empty and fewer-than-three-localization curation inputs
return `unscoreable`, provided their schema is valid. Fit fails clearly if no
clean calibration trace can be scored.

## Public functions

```python
fit_curator_model(
    traces, *, genomic_positions=None, trace_fpr=0.01,
    minimum_reference_observations=20, reference_mode="separation",
    calibration_traces=None,
) -> TraceCuratorModel

curate_trace(trace, model, *, max_removals=3) -> TraceCurationResult
```

All fitting and calibration traces must be clean, resolved single polymers.
The default reference mode pools spatial distances by **exact absolute genomic
separation**. No separation binning, rounding, interpolation, or sparse-bin
fallback is performed. Optional `barcode_pair` mode uses unordered barcode
identity pairs. Unsupported pairs are ignored, including references below the
minimum observation count. Missing relationships are not assigned zero scores.

For independent calibration pass held-out clean `calibration_traces`. Otherwise
the fitting traces calibrate an **in-sample** cutoff; the nominal FPR describes
that empirical calibration sample and is not a guarantee on unseen traces.
Only finite, scoreable clean costs enter the cutoff. Unscoreable calibration
traces are counted separately. Reference arrays are sorted once during fitting;
scoring uses binary searches without sorting them again. Fitting retains pair
observations and scalar calibration costs, without a history of trace-score
objects or pairwise calibration tables.

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

The clean cutoff uses `allowed = floor(trace_fpr * n_clean)` and
`sorted_costs[max(0, n_clean - allowed - 1)]`. Abnormality requires strictly
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
  `n_calibration_unscoreable`, `calibration_source`, and optional `chromosome`.

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
