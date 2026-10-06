# Trace-curator benchmark: cycle 4

`scripts/analyze_trace_curator_iterative_stopping.py` is an exploratory
benchmark, not the production `trace_curator` CLI. It tests whether global
trace stopping needs a secondary deletion safeguard after cycle 3 established
iterative LOO-top3 ranking.

## Policies and safety interpretation

Every policy uses the same global cost (`C_top3`) and candidate ranking
(`C(T) - C(T without i)`). They differ only in whether a proposed removal is
accepted:

- **A_global** removes the best candidate while the global cost exceeds its
  clean-derived threshold.
- **B_absolute** additionally requires the best absolute LOO improvement to
  exceed its held-out-clean threshold.
- **B_relative** instead requires the fractional cost drop to exceed its
  held-out-clean threshold.
- **C_cap** follows A but stops after `--max-removals` (three by default).

The B policies deliberately distinguish abnormality from actionability. An
anomalous trace whose best candidate fails the safeguard stops with
`candidate_not_actionable` and is reported as **`abnormal_unresolved`**. It is
not called clean, and good localizations are not forcibly deleted to normalize
the trace.

Tied best LOO scores are never broken by input row order. The audit records
`n_tied_best` and `candidate_unique`; a globally abnormal trace with a tied
winner stops as `ambiguous_candidate` and remains abnormal/unresolved.

All policies share cached pairwise anomaly/LOO structures for identical trace
states. Candidate removal selects the next top-three nonincident cached edges;
it does not reconstruct every candidate's pair list. References and thresholds
are fit once per fold rather than once per policy.

## Leakage and genomic coordinates

Global and actionability thresholds are calibrated at 0.1%, 1%, and 5%
trace-level FPR from held-out clean traces. For evaluation seed `e`, every
calibration seed `c` is scored with a reference excluding both `e` and `c`.
The primary policy comparison uses the 1% thresholds. Threshold provenance is
written explicitly to `cycle4_thresholds.ecsv`.

The absolute and relative safeguards are independently calibrated from the
same held-out clean folds; they are not conditional quantiles fitted only among
globally abnormal traces. The threshold table therefore reports, for both
calibration and clean evaluation traces, `fraction_global_abnormal`,
`fraction_global_and_absolute_actionable`, and
`fraction_global_and_relative_actionable`. It also records calibration seeds
and each calibration seed's actual nested reference-seed set.

`Genomic_Position` controls genomic order and separation. Barcode identity is
retained only as identity and for the barcode-pair reference. The historical
regular-polymer cycle-3 data may use the explicit
`--assume-uniform-barcode-spacing` fallback; the fallback assigns polymer ranks
after sorting barcode IDs and never treats their numeric differences as genomic
distances. Both separation-pooled and barcode-pair reference modes remain
available through `--reference-mode`.

## Profile, then run

First profile a representative subset:

```bash
python scripts/analyze_trace_curator_iterative_stopping.py \
  --benchmark-root benchmark_curator_loo_validation \
  --output-dir curator_cycle4_memory_profile \
  --assume-uniform-barcode-spacing \
  --reference-mode separation \
  --profile-traces 2000 \
  --profile-calibration-traces 200
```

Both limits are profiling-only: the latter deterministically takes the first
200 clean traces from each calibration seed, so a small run does not first pay
the full calibration cost. Do not use either option for final results.
`cycle4_runtime_profile.ecsv` records elapsed time, component timings,
reference-cache hits and misses, per-group RSS samples, the number of score
builds, seconds per trace for all policies, and a transparent 100,000-trace projection. Progress and
fit timings are also flushed to stderr throughout the run. Inspect this
projection before the full run; if the relevant manifest projection exceeds
6–8 hours, profile and optimize rather than launching it.

Each outer fold also sends the held-out seed's clean baseline traces through
all four complete iterative policies with the primary 1% thresholds and the
outer reference model. These rows are identified by `K = 0`,
`displacement = 0`, and `dataset_type = clean_evaluation`. This evaluation
does not refit the reference, and its already-computed first state seeds the
shared per-trace scoring cache.

Exact full command for the existing cycle-3 historical benchmark:

```bash
python scripts/analyze_trace_curator_iterative_stopping.py \
  --benchmark-root benchmark_curator_loo_validation \
  --output-dir curator_cycle4_results \
  --assume-uniform-barcode-spacing \
  --reference-mode separation \
  --max-removals 3 \
  --min-remaining-localizations 3
```

Omit the simulation fallback and provide `Genomic_Position` in experimental
inputs. A confirmatory pair-reference run changes only
`--reference-mode barcode_pair`.

## Tables

- `cycle4_iteration_audit.ecsv` contains one row per trace/policy/iteration:
  trace size, global cost and threshold, abnormal flag, proposed barcode and
  LOO score, tied-winner count and uniqueness, actionability, counterfactual
  after-cost, absolute/relative drop, acceptance, candidate truth (evaluation
  only), and stop reason.
- `cycle4_trace_metrics.ecsv` reports removal counts, true and false removals,
  recall, precision, exact recovery, remaining corruption, unresolved state,
  terminal `C(T)` and abnormal flag, compact terminal outcome, terminal reason,
  and iterations for every trace and policy. Direct flags distinguish exact
  recovery, false-clean under-curation, abnormal-unresolved stopping, and
  over-curation.
- `cycle4_condition_metrics.ecsv` aggregates recall, removal precision, exact
  recovery, false-removal severity, under-curation, abnormal-unresolved rate,
  iterations, and the removal-count distribution by benchmark condition. For
  held-out clean policy rows it directly reports the fraction with one or two
  removals, mean and maximum removals, terminal-abnormal fraction,
  abnormal-unresolved fraction, and stop-reason distribution.
- `cycle4_thresholds.ecsv` records all three operating points, clean sample
  size, and disjoint reference/evaluation seed provenance.
- `cycle4_runtime_profile.ecsv` contains `scope = group` rows with start/end/
  sampled-peak RSS, trace counts, cache hits/misses, and cache sizes before/after
  clearing, plus a `scope = total` row with component timings, seconds per
  1,000 traces, and projected runtime. Select the total row for overall timing;
  group rows intentionally leave total-only fields empty.

The focused figures cover recall, precision, exact recovery, false genuine
removals, clean probability of any removal, abnormal-unresolved frequency,
iterations/removals, and global-cost trajectories. `clean_any_removal.png`
uses only the explicit `clean_evaluation` policy rows, not any optional K=0
conditions in the cycle-3 manifest. Condition tables retain K, detection
efficiency, displacement, and barcode count so comparisons can be stratified
rather than pooled.


## Memory retention and verification

The memory redesign preserves the empirical anomaly formula, sorted reference
arrays, C_top3 and LOO-top3 definitions, thresholds, nested seed exclusions,
1% primary operating point, all four policies, tie/ambiguity behavior, terminal
classifications, genomic-coordinate semantics, and removal limits. No scientific
or statistical behavior was intentionally changed.

Previously, `clean_states` retained a DataFrame, all pair anomalies, and LOO arrays
for every calibration trace. `evaluation_clean` retained every held-out trace and
its full initial state. `audits_out` and `traces_out` grew throughout the run,
reference caches survived completed groups, and plotting reloaded the full audit.
These collections have been removed:

- Calibration retains only `CompactTraceStats(cost, best_loo,
  best_relative_loo)` for the current fold. Full states are discarded immediately,
  and compact samples are released after computing exact empirical thresholds.
- Held-out clean traces are scored sequentially. Three fixed-size gate counters
  cover all operating points, and each trace immediately runs all four policies.
- Per-trace score caches live inside one function call and are released after its
  four policies. Full calibration/evaluation states never escape into outer lists.
- Audit and trace metrics reuse cycle 1's ECSV chunk writer with separate buffers
  of at most 1,024 rows. Data is appended to hidden `.tmp` ECSV files during the
  run and atomically renamed to the requested filenames on finalization.
- Condition summaries keep integer sums, maxima and distributions per condition
  key rather than individual trace rows. Cost trajectories keep sums/counts by
  policy and iteration. Standalone replotting reads at most 4,096 audit rows at
  once instead of loading the full file through Astropy.
- Only the requested reference mode is fitted. Sorted arrays are reused directly
  by empirical scoring. The reference cache belongs to one `(n_barcodes,
  detection_efficiency)` group and is cleared on both normal and profiling exits.
  Baseline tables are selected with cycle-3 barcode-count provenance: simulation
  metadata first, then the condition-derived simulation-size mapping, then the
  observed maximum barcode when neither source exists. Known nonmatching counts
  are filtered before ECSV loading; counts never come from the current group.
  Baseline tables are released at the group boundary.
- Calibration profiling uses lazy `islice` rather than materializing grouped
  DataFrames before slicing. No `gc.collect()` is needed for the retention fix.

The evaluation/output working set is independent of the number of traces already
processed, for fixed trace sizes and condition keys. Input loading still holds the
current baseline group and one corrupted condition; exact empirical references
scale with the reference population, and calibration stores three scalars per
clean calibration trace. This is not an out-of-core rewrite of input loading or
reference fitting. Increasing input/reference/calibration population sizes can
therefore increase the initial working set even though evaluation retention is
bounded. Allocators may retain freed arenas at group boundaries.

Every progress line includes elapsed time, RSS in MiB, and group/fold context;
calibration and corrupted-condition lines also include the seed or condition ID.
Progress is logged every 250 traces and at group, seed, and condition boundaries.
Linux RSS comes from `/proc/self/statm`; other platforms fall back to the process
high-water value from `resource`. `rss_peak_mb` is the maximum sampled RSS,
sampled at scoring and progress boundaries, not a continuously observed OS peak.

For profiling, `--profile-traces` limits traces passed through all policies.
Clean gate evaluation still covers the full held-out seed, preserving the former
profiling semantics. `--profile-calibration-traces` limits each calibration seed.
Neither limit belongs in final scientific results. Compare successive RSS samples
within a group after reference fitting and initial buffer allocation; a small
allocator warm-up is expected, growth proportional to processed traces is not.

### Validation in this development environment

The historical `benchmark_curator_loo_validation` dataset was unavailable. A
synthetic five-barcode dataset with three seeds and 1,000 clean/corrupted traces
per seed was profiled with both limits (`2000` and `200`). RSS samples were:

| Policy traces processed | RSS (MiB) |
| ---: | ---: |
| 500 | 188.1 |
| 1,000 | 188.2 |
| 1,500 | 191.0 |
| 2,000 | 191.0 |

The sampled group peak was 191.0 MiB; the total sampled peak including plotting
was 203.9 MiB. Elapsed time was 44.8 seconds, or 22.4 seconds per 1,000 traces
including calibration, reference fitting and output. The limited run used one
outer fold: 0 reference-cache hits, 3 misses, and cache size 3 before clearing / 0
after clearing. These measurements demonstrate retention behavior on synthetic
inputs, not runtime or memory requirements for the historical benchmark.

A direct tiny-benchmark comparison against the original script matched all
scientific tables: 72 trace-metric rows, 92 audit rows, 8 condition-summary rows,
and 9 threshold rows. Tests also cover compact/full-state calibration equivalence
including NaNs, streamed summary/trajectory equivalence, chunk finalization,
unused reference modes, lazy profiling, bounded live scoring/audit objects,
reference release across groups, and deterministic complete tiny benchmarks in
both reference modes. The two existing assertion issues were corrected in tests:
the cap test now avoids simultaneously reaching the minimum trace size, and the
NumPy/math logarithm comparison permits floating-point rounding.

Run the focused regression suite with:

```bash
python -m pytest test/test_trace_curator_iterative_stopping.py \
  test/test_trace_curator_loo_validation.py \
  test/test_trace_curator_benchmark_analysis.py \
  test/test_trace_curator_attribution_benchmark.py -q
```

Final checks: all 34 cycle-4 tests passed. The four-script regression suite
reported 72 passed and one existing cycle-3 failure,
`test_iterative_rescoring_can_find_second_culprit`. That failure reproduces in
isolation, and its implementation and test are unchanged from the fetched feature
branch. Black formatting, isort import checks, flake8 on the changed Python files,
and `git diff --check` passed.

The mixed-population regression uses five- and seven-barcode simulations at the
same detection efficiency and checks reference-fitting and calibration membership
against cycle-3's provenance logic. It covers explicit simulation metadata,
condition-derived counts, and the observed-barcode fallback, and verifies that
known nonmatching simulations are excluded before ECSV loading. All three cases
fail against the pre-fix implementation and pass after the correction.
