# Trace-curator benchmark: cycle 5

`scripts/analyze_trace_curator_tie_resolution.py` tests only whether exact
LOO-top3 ties can be safely resolved. It is standalone exploratory analysis;
it changes neither the production API nor the cycle-4 script. The conservative
default remains `stop_unresolved`. Improved recall alone is insufficient grounds
for recommending a tie breaker.

## Source data and folds

Use the existing completed `benchmark_curator_loo_validation/` sweep. The script
never generates data or invokes the simulator. It requires all eight primary
classes: B=25, K=2/3, p=0.5/0.8, d=0.4/0.8 µm. It also requires the four matching
K1 classes and evaluates held-out clean baselines for both efficiencies.
Missing classes fail explicitly rather than silently producing a partial study.
All replicate seeds come from the complete cycle-4 `(B, p)` condition group,
including conditions outside this narrow evaluation subset.

Baseline membership follows cycle-3 provenance: simulation metadata, then the
condition-derived simulation-size mapping, then the observed maximum barcode
only if both metadata sources are absent. Known nonmatching counts are excluded
before ECSV loading. The current group never supplies a simulation's count.

The script imports cycle-4 empirical scoring, coordinate handling, reference
fitting, exact seed exclusions, compact calibration, and summary utilities.
For evaluation seed e, clean calibration seed c is scored against a reference
excluding both e and c. The primary global threshold is the unchanged 1% cutoff.
Original benchmark observations and the reconstructed held-out reference are
always authoritative. `--cycle4-results` optionally records a directory as
provenance; it does not preselect traces, replace thresholds, or supply scores.

## Stage 1: offline tie events

Each trace first runs the actual cycle-4 `C_cap` policy. Only a terminal
`ambiguous_candidate` stop is a stage-1 event. A tie encountered after the removal
cap, at minimum trace size, or on a normal/unscoreable trace is not an eligible
event. This preserves cycle-4 stopping precedence and does not change its path.
There can be at most one offline event per trace: cycle 4 stops at its first
eligible ambiguity. Iterations 0, 1, and 2 can all contribute events.

Secondary scoring is restricted to candidates with exactly equal maximal finite
LOO-top3 scores. No tolerance is introduced and no winner is chosen by row order.

- `stop_unresolved` leaves the ambiguity unresolved.
- `loo_top5` ranks only the primary tied candidates by the sum of the five largest
  finite pair anomalies before removal minus that sum after removal.
- `k1_local` uses the nearest detected genomic neighbor on each side of each
  candidate, identified by `Genomic_Position`. It takes the mean of the finite
  available flank anomalies already scored by cycle 4. An unscoreable nearest
  flank is not replaced by a more distant detected neighbor. One scoreable flank
  is sufficient; neither scoreable flank makes the candidate unscoreable.
- `one_step_lookahead` provisionally removes each primary tied candidate, runs the
  established cycle-4 scorer on the reduced trace, and ranks by its maximum
  finite LOO-top3 score. It does not rank by immediate after-cost and does not
  recursively search. No finite reduced LOO makes that candidate unscoreable.

Secondary selection uses only finite scores. A unique finite maximum resolves
an event even when another candidate is unscoreable; those events are explicitly
marked `partially_scoreable` and remain subject to measured precision. A secondary
tie or no scoreable candidate remains unresolved. Truth labels enter only
post-selection diagnostics and evaluation, never scoring or winner selection.

Tie classes form the following mutually exclusive labels:
`no_true_candidate`, `all_tied_candidates_true`,
`single_true_candidate_in_tie`, and otherwise `mixed_tie`. The separate
`mixed_tie` boolean and mixed-tie fraction include single-true mixed sets.
The oracle-resolvable fraction is the fraction of events with at least one
true culprit inside the primary tied set, not a claim that a method can find it.

## Stage 2: integrated capped policies

The four variants are `C_cap_stop`, `C_cap_top5`, `C_cap_k1`, and
`C_cap_lookahead`. `C_cap_stop` uses the actual cycle-4 implementation. Each
other variant retains the same global gate, unique primary decisions, minimum
remaining size of three, maximum of three removals, and terminal classifications.
It invokes its secondary method only at an eligible primary tie. Secondary
methods can encounter additional ties on their own subsequent paths; these are
included in integrated tie-invocation metrics, not in the offline event table.
The stage-1 event denominator consequently differs from the integrated denominator.
Conservative stopping at an eligible tie counts as one unresolved tie invocation.

## Commands

First profile the offline analysis:

```bash
python scripts/analyze_trace_curator_tie_resolution.py \
  --benchmark-root benchmark_curator_loo_validation \
  --output-dir curator_cycle5_tie_profile \
  --assume-uniform-barcode-spacing \
  --reference-mode separation \
  --stage offline \
  --profile-ambiguous-events 500
```

Run stage 1 without profiling limits:

```bash
python scripts/analyze_trace_curator_tie_resolution.py \
  --benchmark-root benchmark_curator_loo_validation \
  --output-dir curator_cycle5_offline \
  --assume-uniform-barcode-spacing \
  --reference-mode separation \
  --stage offline
```

Run both stages on the same narrow data selection:

```bash
python scripts/analyze_trace_curator_tie_resolution.py \
  --benchmark-root benchmark_curator_loo_validation \
  --output-dir curator_cycle5_tie_resolution \
  --assume-uniform-barcode-spacing \
  --reference-mode separation
```

Use the spacing fallback only for the historical regular-polymer data. Inputs
with real genomic coordinates should supply `Genomic_Position`. The established
`barcode_pair` reference mode remains supported without changing genomic flank
selection.

`--profile-ambiguous-events` caps unique primary tie states requiring analysis,
counting offline and integrated paths together. Secondary results for the same
trace state are reused. If a trace needs another state beyond the budget, all
its pending event/metric output is discarded so policy comparisons contain
identical completed traces. The runtime profile distinguishes calculated tie
states from emitted offline events and marks the cutoff. Calibration is never
truncated: exact original thresholds require all clean calibration traces.
Profiling output is a deterministic partial sample, not final scientific results.
Do not use profiling limits for a production decision.

## Outputs and denominators

- `cycle5_tie_events.ecsv`: one row per offline ambiguity, identifiers and fold
  provenance, costs/threshold, iteration, tie size and JSON-encoded Spot_ID/
  Barcode sets, remaining truth count, and oracle classification.
- `cycle5_tie_method_scores.ecsv`: one row per event/method/tied candidate,
  scoreability, secondary score, unique selected identity/truth, and separate
  k1 left/right neighbor identities and anomaly values. Selected identity is
  empty for unresolved decisions. Candidate/winner barcodes are identity strings.
- `cycle5_tie_resolution_performance.ecsv`: stratum results by method/K/p/d and
  pooled results separately for primary data and controls. Counts provide explicit
  denominators. Score coverage is events with at least one finite secondary
  score divided by all ambiguous events. Unique-resolution coverage is resolved
  events divided by all events. Precision is correct divided by resolved events.
  Incorrect-resolution rate is incorrect divided by all events; effective correct
  resolution is correct divided by all events. Unresolved fraction is unresolved
  divided by all events. Undefined denominators produce NaN, not zero precision.
- `cycle5_iterative_performance.ecsv`: primary K2/K3 recall, removal precision,
  exact recovery, false-removal severity, abnormal/unresolved and tie-unresolved
  fractions, removals, tie invocation/resolution/precision metrics, and differences
  in recall/recovery/false removals relative to `C_cap_stop` in each stratum.
- `cycle5_clean_control_performance.ecsv`: the same accounting for held-out clean
  controls, including any removal, mean removals, at least two removals, maximum
  removals, and terminal abnormal fraction.
- `cycle5_k1_control_performance.ecsv`: K1 recall, precision, exact recovery,
  primary tie frequency, and integrated tie metrics.
- `cycle5_runtime_profile.ecsv`: group and total timings, sampled RSS, completed
  trace/event counts, reference cache hits/misses and reset sizes, per-trace cache
  diagnostics, selection/stage/profiling provenance. Score-build diagnostics
  describe completed evaluation traces; calibration is included in elapsed time.
- `cycle5_thresholds.ecsv`: reconstructed primary 1% cutoffs and complete seed
  provenance, including each nested reference seed set.
- `cycle5_trace_metrics.ecsv`: streamed metrics for each complete trace/policy.
- `cycle5_report.txt`: event numbers, oracle ceiling, method precision and both
  coverage definitions, integrated stratum differences, clean/K1 behavior, runtime
  and sampled peak RSS. It never automatically recommends a method for production.

The seven figures compare methods with conservative stopping separately for K2
and K3, retaining all four p/d strata. Precision is undefined for the conservative
baseline because it never resolves a tie. Offline-only runs retain readable,
empty stage-2 tables/figures; empty files mean that stage was not evaluated.

## Memory and validation

Detailed outputs reuse the cycle-4 ECSV writer with 512-row buffers and atomic
finalization. Zero-event outputs still have readable schemas. Only compact
condition/event aggregates remain across traces. Per-trace scored states use a
64-entry LRU shared across methods, including reduced lookahead states, and are
cleared even at a profiling cutoff. Secondary-result caches live for one trace.
Reference caches are cleared at each `(25, p)` boundary and profiling exit.
No full scoring structures or audit DataFrames accumulate across traces.

Baseline groups, the current condition's input table, exact empirical reference
arrays, and three compact calibration scalars per clean trace still scale with
their input populations. This is bounded evaluation retention, not out-of-core
reference fitting. Progress every 250 completed traces or 100 offline ambiguities
includes elapsed time, RSS, condition, evaluation seed, trace count and ambiguity
count. Calibration progress occurs every 250 clean traces. Linux RSS is live
resident memory; other platforms use the cycle-4 high-water fallback. Peaks are
sampled at scoring/progress boundaries, not continuously observed OS peaks.

Run focused tests with:

```bash
python -m pytest test/test_trace_curator_tie_resolution.py \
  test/test_trace_curator_iterative_stopping.py -q
```

Unit cases exercise exact tie detection, secondary ties and permutation invariance,
restricted top5 candidates, genomic k1 and missing support, reduced-state lookahead,
truth-blind selection, unchanged non-tied empirical decisions, clean controls,
cap precedence, fold exclusion, deterministic output, readable empty schemas,
streaming/cache bounds, and fair profiling cutoffs. The driver test mocks input
readers with tiny in-memory cases; it does not create a benchmark sweep or use
unit-case results as scientific evidence.

The historical dataset was unavailable in the development workspace. The CLI
was checked to fail explicitly on its missing manifest without creating data.
Historical ambiguous-event counts, oracle-resolvable fractions, tie-breaker
precision/coverage, recovery/false-removal changes, clean-control behavior,
runtime and peak RSS remain unmeasured. No production recommendation can be
made from these unit tests.

Development validation: all 48 cycle-5 cases pass. Across the five benchmark
analysis suites, 120 tests pass and the pre-existing cycle-3
`test_iterative_rescoring_can_find_second_culprit` fails. Cycle-4 and production
code are unchanged by this work.

Before considering production, inspect the full, unprofiled results: precision
around or above 90% in strong-signal K2/K3 conditions, clear recovery improvements,
no meaningful false-removal increase, essentially unchanged clean behavior,
and no collapse in any p/d stratum. Retain unresolved stopping if these joint
requirements are not supported.
