# trace_splitter

**Reliability status**: `stable`

`trace_splitter` curates repeated barcode assignments, identifies traces
consistent with one versus two underlying polymers, resolves candidate doublets
into separate traces, and can reject spatially incompatible localizations.

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_splitter.parse_arguments
   :prog: trace_splitter
```

## Recommended barcode-aware workflow

The recommended workflow uses `--method beam` with
`--multiplicity-classifier likelihood`. Its default nuisance model is the
validated `--likelihood-off-target-model global`. Classification and spatial
resolution are separate stages.

### 1. Multiplicity classification

A **singlet** contains one underlying polymer; a **doublet** contains two
polymers merged into one reconstructed trace. For each barcode identity
represented in the dataset, including zero detections of that barcode in
individual traces, the likelihood classifier models

```text
T_i ~ Binomial(m, p)
A_i ~ Poisson(lambda)
K_i = T_i + A_i
```

Here `m` is one or two polymers, `p` is the dataset-wide genuine detection
efficiency, `lambda` is the off-target detection rate, and `K_i` is the observed
number of detections for barcode `i`. The dataset fit estimates `p`, doublet
fraction `pi`, and the global `lambda`. Per-trace singlet and doublet likelihoods
produce a posterior doublet probability; the default decision threshold is
`0.5`. The resulting states are `UNCHANGED` (no repeated barcode and classified
singlet), `CLEAN_ONE` (one polymer requiring duplicate cleanup), and
`RESOLVE_TWO` (two-polymer spatial resolution).

### Why `global` is the recommended default

In the tested simulation regimes, likelihood classification substantially
outperformed fixed repeated-barcode thresholds at low detection efficiency while
retaining high singlet specificity. The global nuisance model remained stable
across barcode-efficiency SD values from 0 to 0.15. In contrast, `auto`
increasingly selected barcode-specific nuisance rates as genuine barcode
efficiencies varied, sometimes substantially reducing sensitivity and balanced
accuracy at low mean efficiency. These conclusions are limited to the tested
simulation regimes; nevertheless, `global` is the validated and recommended
default.

### 2. Beam-search spatial resolution

Beam search learns an empirical spatial-distance model conditioned on genomic
separation (or barcode-index separation when genomic coordinates are absent).
The current default settings are:

```text
history_mode = multi
distance_score = residual
beam_width = 100
rejection_cost = 16
doublet_posterior_threshold = 0.5
likelihood_off_target_model = global
minimum_confidence = 0.05
```

The benchmark series exercised settings including `history_mode=multi`,
`distance_score=residual`, and `rejection_cost=16`. The
`minimum_confidence=0.05` value shown here is the current software default; it
should not be interpreted as benchmark-validated because the benchmarks
generally used `minimum_confidence=0`.

`history-mode=multi` scores a candidate against up to `--history-length` prior
localizations; residual scoring uses squared standardized distance residuals.
`--beam-width` controls the retained search hypotheses. One-polymer cleanup
chooses among repeated detections, whereas two-polymer resolution assigns
compatible detections to two polymers. Incompatible detections may remain
unassigned: `--rejection-cost` controls that tradeoff, and
`--minimum-polymer-size` prevents undersized inferred polymers. A two-polymer
solution below `--minimum-confidence` is removed as ambiguous. Confidence is a
normalized score gap, not a probability.

The empirical model borrows nearby genomic separations when a separation has
fewer than `--model-min-observations`; `--variance-floor` avoids degenerate
scores. If no clean traces exist, a warning announces the weaker
nearest-candidate-pairs fallback. Missing intermediate barcodes incur no
separate penalty.

## Optional nuisance models

- `global` uses one nuisance rate shared by every barcode. It is the recommended
  default.
- `barcode-specific` fits one rate per barcode, regularized toward the global
  estimate by `--likelihood-lambda-regularization`. Use it when barcode-specific
  off-target contamination is independently expected.
- `auto` tests barcode-specific excess multiplicity with
  `--likelihood-heterogeneity-alpha` and may switch to barcode-specific rates.
  Genuine barcode-specific detection efficiencies can produce the same pattern,
  so use this mode cautiously.

## Reliability warnings

Two distinct dataset-level warnings do not alter posterior probabilities or
classification decisions:

1. **Low-information / identifiability.** With few barcode identities
   represented in the dataset and low fitted `p`, singlets and doublets contain
   fundamentally too few informative observations to separate reliably.
   Simulations found that 5 represented barcodes around `p=0.3` could approach
   chance-level balanced accuracy; 10 represented barcodes at low `p` also lost
   reliability, especially for strongly imbalanced mixtures. In the tested
   regimes, 25 or more represented barcodes were robust across a broad range of
   doublet fractions.
2. **`auto` nuisance-model ambiguity.** A dataset may contain ample information
   but genuine barcode-efficiency differences can mimic barcode-specific
   off-target rates. Benchmark v6 found that `global` stayed stable as efficiency
   SD increased while `auto` selected barcode-specific nuisance more often and
   could lose sensitivity and balanced accuracy, particularly at low mean `p`.
   Therefore, if `auto` selects barcode-specific rates at fitted `p <= 0.35`, one
   caution is printed for the dataset. Prefer `global` unless independent
   evidence supports barcode-specific contamination.

## Threshold classifier

`--multiplicity-classifier threshold` remains fully supported and is the
backward-compatible classifier default. It combines
`--split-min-repeated-barcodes`, `--split-min-repeated-fraction`, and
`--candidate-rule {both,either}`. This classifier is simple and deterministic,
but its repeated-barcode cutoffs depend strongly on detection efficiency and
are generally less sensitive to doublets at low `p`.

## Alternative spatial clustering methods

`--method kmeans` and `--method hdbscan` remain available for generic spatial
partitioning, but they are not the recommended barcode-aware singlet/doublet
workflow. These older methods select extended candidate traces using radius of
gyration (`mean + --std_threshold * SD`); `--split-all` bypasses selection.
K-means uses the predefined `--num_clusters`. HDBSCAN uses density clustering,
supports `--min-cluster-size`, `--min-samples`, selection method and epsilon,
and can identify noise points. For backward compatibility, its noise rows retain
the original `Trace_ID`; these methods do not reject input rows.

## CLI examples

Minimal recommended likelihood workflow (uses `global` implicitly):

```console
trace_splitter \
  --input traces.ecsv \
  --output traces_split.ecsv \
  --method beam \
  --multiplicity-classifier likelihood
```

The same selection with the recommendation explicit:

```console
trace_splitter \
  --input traces.ecsv \
  --method beam \
  --multiplicity-classifier likelihood \
  --likelihood-off-target-model global
```

Threshold classification:

```console
trace_splitter --input traces.ecsv --method beam \
  --multiplicity-classifier threshold --split-min-repeated-barcodes 3 \
  --split-min-repeated-fraction 0.3 --candidate-rule both
```

Alternative spatial clustering:

```console
trace_splitter --input traces.ecsv --method kmeans --num_clusters 2
trace_splitter --input traces.ecsv --method hdbscan \
  --min-cluster-size 3 --min-samples 3
```

## Diagnostics output

Beam mode writes `[tracefile]_split_diagnostics.ecsv` (override with
`--diagnostics-output`). Each input trace reports requested and inferred polymer
number, repeated-barcode counts and fractions, likelihood ratio, posterior
doublet probability, fitted `p`, fitted `pi`, global nuisance rate, selected
nuisance model, heterogeneity p-value, unassigned-localization fraction, and
classifier reliability level. Dataset-level reliability and `auto` warning
levels/messages are ECSV metadata. When barcode-specific rates are used, the
per-barcode rates are stored as the
`classifier_barcode_off_target_rates` metadata mapping rather than table
columns. Diagnostics also contain spatial best/alternative scores, score gap,
confidence, and ambiguity counts.

## Output files

1. `[tracefile]_split.ecsv` contains resolved traces. Beam mode can omit rejected
   detections or ambiguous traces; `Spot_ID` is preserved and successful splits
   receive new `Trace_ID` values. A single-file `--output` selects an exact path.
2. `[tracefile]_split_diagnostics.ecsv` contains the beam diagnostics described
   above.

Coordinates are converted to contiguous 64-bit floats before clustering to
avoid known scikit-learn HDBSCAN failures with some float32 ECSV inputs. The
legacy no-method behavior remains radius-of-gyration selection followed by
K-means for backward compatibility.
