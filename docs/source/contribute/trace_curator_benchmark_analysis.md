# Exploratory trace-curator benchmark analysis

`scripts/analyze_trace_curator_benchmark.py` is a development utility for
comparing candidate scores. It is deliberately **not** installed as a
`traceratops` command, does not curate or delete localizations, and does not
generate simulations.

## Example

```bash
python scripts/analyze_trace_curator_benchmark.py \
  --benchmark-root benchmark_curator_v1 \
  --output-dir curator_model_analysis
```

The input must have the `sweep_manifest.yaml`, `simulations/SIM*`, and
`conditions/C*` layout produced by `trace_curator_benchmark.sh`. At least three
replicate seeds per detection efficiency are required. The reader joins the
observed and ground-truth tables by `Spot_ID` and fails on missing required
columns, duplicate barcode identities, or inconsistent trace/barcode identity.

## Leakage controls

For each evaluation replicate, reference distributions are fitted from clean
`simulations/` baselines belonging to every *other* seed at the same detection
efficiency. Conditions are never used as reference data. Fixed-FPR thresholds
use a nested cross-fit: each non-evaluation calibration replicate is scored by
a model excluding both that replicate and the evaluation replicate. Thus the
evaluated localization contributes neither to its reference distribution nor
to its operating-point calibration. Coordinates and barcode identity are
copied into a feature-only table before scoring; corruption fields are attached
afterward solely for evaluation.

Whenever rows from multiple baselines are concatenated, traces are grouped by
the globally unique `(simulation_id, Trace_ID)` pair. This prevents repeated
trace labels from independent simulation replicates from being combined.

## Scores

All scores increase with anomalousness.

* `trace_splitter_like_residual` is a residual-only baseline inspired by the
  distance cost used by `trace_splitter`. It is **not** the beam-search resolver.
  For every usable relationship it computes the squared standardized residual
  from the clean mean and standard deviation at that exact genomic barcode
  separation, then takes the mean over context. Missing separation bins are not
  substituted by default. `--residual-separation-fallback nearest` explicitly
  enables nearest-bin substitution; the output records the mode and number of
  fallbacks used.
* `empirical_separation_*` uses the clean empirical distance distribution for
  the absolute barcode separation.
* `empirical_barcode_pair_*` uses the clean distribution for that exact
  unordered barcode pair. A distribution with fewer than
  `--minimum-reference-observations` (20 by default) is not extrapolated.

The empirical two-sided probability is
`min(1, 2 * min((# <= d + 1)/(n + 1), (# >= d + 1)/(n + 1)))`; add-one
smoothing prevents zero probabilities. Pair anomaly is `-log10(probability)`.
The exported directional percentile uses the midrank of ties with the same
finite-sample smoothing: `(# < d + 0.5 * # == d + 0.5)/(n + 1)`. Values below
0.5 indicate short-distance anomalies and values above 0.5 indicate
long-distance anomalies.
Suffixes select mean, median, maximum, second-largest, or the count/fraction of
relationships outside the empirical 95% or 99% interval. `k1`, `k2`, and `k3`
select that many nearest detected barcodes independently on each genomic side;
`all` selects all other detections. A missing side is allowed and reported. A
missing score means zero context or insufficient reference support, never
evidence that a localization is normal.

`count_95` and `count_99` are diagnostics whose scale increases with the number
of usable relationships. They must not be compared directly between context
sizes. The corresponding `fraction_95` and `fraction_99` scores are the
context-size-normalized alternatives.

## Outputs

* `localization_scores.ecsv`: one row per localization, model, and context,
  including score status, usable/requested context counts, insufficient-pair
  count, minimum reference support, fallback count/mode, nearest-flank and
  both-flanks diagnostics, condition metadata, and evaluation-only ground truth.
* `pairwise_relationship_scores.ecsv`: one row per evaluated empirical
  relationship and context. It reports distance, reference support, empirical
  percentile, two-sided tail probability, pair anomaly, and short/long
  direction. Insufficient references remain explicit rows with missing scores.
* `reference_model_summary.ecsv`: support for every separation and barcode-pair
  distribution in each evaluation fold.
* `model_performance.ecsv`: per-condition/per-replicate ROC AUC and PR AUC.
* `performance_at_fixed_fpr.ecsv`: thresholds calibrated at 0.1%, 1%, and 5%
  target clean-localization FPR, plus observed sensitivity, FPR, precision, and
  denominators.
* `plot_aggregate_values.ecsv`: replicate-first mean, SD, count, and SEM used by
  plots. Error bars are SEM across independent replicate results.
* `sensitivity_efficiency_*.png`, `performance_by_context.png`,
  `score_distributions.png`, and `zero_displacement_sanity.png`: compact
  diagnostics for the requested primary and sanity-check comparisons.

The supplied generator does not confine singlet polymers to the stated nuclear
diameter and does not clip displaced positions. Results must therefore not be
interpreted as testing nuclear-boundary consistency. The zero-displacement arm
has no positive `is_corrupted` labels; selected targets are retained explicitly
for the false-positive sanity plot.
