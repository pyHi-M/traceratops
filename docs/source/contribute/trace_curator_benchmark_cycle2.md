# Trace-curator benchmark: cycle 2

`scripts/analyze_trace_curator_benchmark_cycle2.py` is an exploratory analysis
of source attribution. It is not an installed command, does not curate traces,
and neither generates nor changes simulations. The cycle-1 script remains the
reproducible baseline.

## Run the frozen benchmark

The historical homogeneous simulation has regular polymer beads whose numeric
barcode values encode polymer order. The fallback sorts those identifiers and
assigns equally spaced, one-based polymer ranks; it does not use the identifier
itself as a genomic coordinate. That assumption is enabled explicitly:

```bash
python scripts/analyze_trace_curator_benchmark_cycle2.py \
  --benchmark-root benchmark_curator_v1 \
  --output-dir curator_model_analysis_cycle2 \
  --assume-uniform-barcode-spacing
```

For any experimental or nonuniform design, provide an ECSV or CSV mapping
instead. It must contain unique `Barcode #` (or `Barcode`) and
`Genomic_Position` columns:

```bash
python scripts/analyze_trace_curator_benchmark_cycle2.py \
  --benchmark-root benchmark_curator_v1 \
  --output-dir curator_model_analysis_cycle2 \
  --genomic-coordinates barcode_coordinates.ecsv
```

Barcode identity is never treated as genomic position unless the simulation-only
flag is present. All contexts select the nearest *detected* loci in genomic
coordinate order; bridge gaps and separation references use the mapped values.

## Models and leakage controls

The analysis retains genomic-separation and exact barcode-pair k1 baselines.
It compares leave-one-localization-out changes in four whole-trace costs, edge
concentration, and a robust nearest-flank bridge residual. The bridge score
standardizes its three distances marginally with median/MAD estimates and does
not model their covariance; it is not a full joint probabilistic model. A deterministic
logistic-regression comparator combines a fixed, documented feature list.

The separation-pooled reference currently requires an exact genomic-separation
match. This is directly appropriate for the regularly spaced frozen benchmark.
Arbitrary experimental probe spacing will require a later, explicitly evaluated
binning, interpolation, or smoothing design. Exact barcode-pair references do
not have this limitation and remain available.

The outer fold holds out a complete simulation seed and every displacement
condition derived from it. References and supervised training use other seeds.
Threshold calibration is nested: each calibration seed and the evaluation seed
are excluded from its reference and ML training. In addition, every ML training
seed is scored with a separately fitted statistical reference that excludes that
training seed's own clean baseline. All displacement conditions for a seed stay
together. Labels and the condition, seed, displacement, and simulation identity
are not classifier features.
The ML result is only an exploratory upper bound: a homogeneous polymer, one
fixed corruption mechanism, one displaced localization, and fixed radial
displacements do not establish generalization to experimental errors.

Conditions are scored and evaluated one at a time. Detailed localization output
is disabled by default; `--write-localization-scores` enables it. Full
observation, condition-feature, and long-score tables are released before
advancing; only the compact ML columns described below may remain cached.

Within one detection-efficiency block, deterministic caches reuse reference
models keyed by their sorted reference-seed tuple and compact ML feature tables
keyed by condition identity plus that exact tuple. Cached feature tables contain
only the fixed ML features, label, seed, and reference-seed provenance. The
cache never retains observations or long-score tables. Compact features beyond
the in-memory limit are stored in a temporary disk-backed cache, and all cache
state is released before the next efficiency. Progress output distinguishes
fits/scores from cache reuse and reports fit, scoring, and hit counts per block.

## Result tables

* `cycle2_model_performance.ecsv`: per-condition ROC and PR AUC.
* `cycle2_fixed_fpr_performance.ecsv`: conditional and effective sensitivity,
  scoreability coverage, observed FPR, and precision at 0.1%, 1%, and 5%
  clean-localization FPR. Effective sensitivity counts unscored corruptions as
  undetected.
* `cycle2_attribution_performance.ecsv`: per-corrupted-trace top-1/top-2 result,
  attribution coverage, minimum rank, unique and tie-inclusive top-1 results,
  number tied at the best score, score margin, and scored trace size. Unique
  top-1 is the primary attribution metric; unconditional metrics count an
  unscoreable corrupted localization as failure.
* `cycle2_collateral_calls.ecsv`: genuine calls per corrupted trace at 1% FPR.
* `cycle2_rank_distribution.ecsv`: the corrupted localization rank for plotting.
* `cycle2_reference_summary.ecsv`: fold-specific support for genomic-separation,
  barcode-pair, and bridge references.
* `cycle2_ml_feature_importance.ecsv`: standardized logistic coefficients.
* `cycle2_aggregate_performance.ecsv`: mean, SD, SEM, and replicate count used
  to summarize fixed-FPR, attribution, and collateral metrics.

Line plots first average within a seed and report mean ± SEM across independent
seeds. The aggregate table likewise reports replicate-first mean, SD, and SEM.
The rank-distribution table supports condition-specific rank plots without
requiring the much larger optional localization table.
