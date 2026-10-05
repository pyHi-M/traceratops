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
  --output-dir curator_cycle4_profile \
  --assume-uniform-barcode-spacing \
  --profile-traces 1000
```

`cycle4_runtime_profile.ecsv` records elapsed time, seconds per trace for all
policies, and a transparent 100,000-trace projection. Inspect this projection
before the full run; if the relevant manifest projection exceeds 6–8 hours,
profile and optimize rather than launching it.

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
  LOO score, actionability, counterfactual after-cost, absolute/relative drop,
  acceptance, candidate truth (evaluation only), and stop reason.
- `cycle4_trace_metrics.ecsv` reports removal counts, true and false removals,
  recall, precision, exact recovery, remaining corruption, unresolved state,
  terminal reason, and iterations for every trace and policy.
- `cycle4_condition_metrics.ecsv` aggregates recall, removal precision, exact
  recovery, false-removal severity, under-curation, abnormal-unresolved rate,
  iterations, and the removal-count distribution by benchmark condition.
- `cycle4_thresholds.ecsv` records all three operating points, clean sample
  size, and disjoint reference/evaluation seed provenance.
- `cycle4_runtime_profile.ecsv` supplies measured and projected runtime.

The focused figures cover recall, precision, exact recovery, false genuine
removals, clean probability of any removal, abnormal-unresolved frequency,
iterations/removals, and global-cost trajectories. Condition tables retain K,
detection efficiency, displacement, and barcode count so comparisons can be
stratified rather than pooled.
