# trace_filter

**Reliability status**: `stable`

## Localization quality filtering

`trace_filter` can remove spots by comparing values in a localization table
against user-provided minimum thresholds. Pass the localization table with
`--localization_file`, then add any combination of quality filters:

- `--intensity_min`: minimum spot intensity. This single argument supports both
  localization-table formats: it filters on `mean_intensity` when that column is
  present in newer tables, or on `peak` for legacy tables.
- `--snr_min`: minimum `snr`.
- `--spot_pixel_percentage_min`: minimum `spot_pixel_percentage`.
- `--skew_min`: minimum `skew`.
- `--patch_size_min`: minimum `patch_size`.
- `--roundness_min`: minimum `roundness`.
- `--object_class_min`: minimum `object_class`; use `--object_class_min 1` to
  keep spots classified as spots (`object_class >= 1`) and remove background
  localizations (`object_class == 0`).

When multiple quality filters are provided, a spot is kept only if it satisfies
all requested minimum thresholds.

Example:

```bash
trace_filter \
  --input Trace.ecsv \
  --localization_file Localizations.ecsv \
  --intensity_min 1000 \
  --snr_min 5 \
  --object_class_min 1
```

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_filter.parse_arguments
   :prog: trace_filter
```
