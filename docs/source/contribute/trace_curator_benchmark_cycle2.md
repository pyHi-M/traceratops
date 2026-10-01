# Trace-curator benchmark: cycle 2

`scripts/analyze_trace_curator_benchmark_cycle2.py` is an exploratory analysis
of source attribution. It is not an installed command, does not curate traces,
and neither generates nor changes simulations. The cycle-1 script remains the
reproducible baseline.

## Run the frozen benchmark

The historical homogeneous simulation has regular polymer beads whose numeric
barcode values encode polymer order. That assumption is enabled explicitly:

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
concentration, and a joint nearest-flank bridge residual. A deterministic
logistic-regression comparator combines a fixed, documented feature list.

The outer fold holds out a complete simulation seed and every displacement
condition derived from it. References and supervised training use other seeds.
Threshold calibration is nested: each calibration seed and the evaluation seed
are excluded from its reference and ML training. Labels and the condition, seed,
displacement, and simulation identity are not classifier features.
The ML result is only an exploratory upper bound: a homogeneous polymer, one
fixed corruption mechanism, one displaced localization, and fixed radial
displacements do not establish generalization to experimental errors.

Conditions are scored and evaluated one at a time. Detailed localization output
is disabled by default; `--write-localization-scores` enables it. Calibration,
condition features, and scores are released before advancing to the next fold
or condition.

## Result tables

* `cycle2_model_performance.ecsv`: per-condition ROC and PR AUC.
* `cycle2_fixed_fpr_performance.ecsv`: sensitivity, observed FPR, and precision
  at 0.1%, 1%, and 5% clean-localization FPR.
* `cycle2_attribution_performance.ecsv`: per-corrupted-trace top-1/top-2 result,
  rank, score margin, and scored trace size.
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
