# trace_filter

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_filter.parse_arguments
   :prog: trace_filter
```

## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_[output_tag]_[label_tag].ecsv`: Filtered trace table. The output tag comes from `--output` (default: `filtered`), and the label tag is appended only when `--keep_label` or `--remove_label` is used.

2. `[localizationfile]_localization_intensities.[format]`: Intensity-distribution plot for the input localization table. This file is written when `--localization_file` and `--intensity_min` are provided.

3. `[tracefile]_filtered_intensities.[format]`: Intensity-distribution plot for localizations kept after intensity filtering. This file is written when `--localization_file` and `--intensity_min` are provided.

4. `[tracefile]_[output_tag]_[label_tag]_localizations.ecsv`: Filtered localization table containing localizations whose `Buid` values remain in the final filtered trace table. This file is written when `--localization_file` is provided.

5. `[tracefile]_[output_tag]_[label_tag]_localization_distribution_fluxes.[format]`: Quality-control plot summarizing the filtered localization table. This file is written when `--localization_file` is provided.

## Examples

```bash
trace_filter \
  --input Trace.ecsv \
  --localization_file Localizations.ecsv \
  --intensity_min 1000 \
  --snr_min 5 \
  --object_class_min 1 \
  --spot_pixel_percentage_max 0.7 \
  --roundness_max 0.8
```

## Notes

### Removal of barcodes by intensity

```bash
trace_filter --input path/to/your/trace_file.ecsv
```

`trace_filter` can remove spots by comparing values in a localization table
against user-provided quality thresholds. Pass the localization table with
`--localization_file`, then add any combination of quality filters.

Minimum-threshold filters remove spots whose value is smaller than the provided
threshold:

- `--intensity_min`: minimum spot intensity. This single argument supports both
  localization-table formats: it filters on `mean_intensity` when that column is
  present in newer tables, or on `peak` for legacy tables.
- `--snr_min`: minimum `snr`.
- `--skew_min`: minimum `skew`.
- `--patch_size_min`: minimum `patch_size`.
- `--object_class_min`: minimum `object_class`; use `--object_class_min 1` to
  keep spots classified as spots (`object_class >= 1`) and remove background
  localizations (`object_class == 0`).

Maximum-threshold filters remove spots whose value is larger than the provided
threshold:

- `--spot_pixel_percentage_max`: maximum `spot_pixel_percentage`.
- `--roundness_max`: maximum `roundness`.

When multiple quality filters are provided, a spot is kept only if it satisfies
all requested thresholds. In other words, it must be greater than or equal to
each requested minimum and less than or equal to each requested maximum.

When a localization table is provided, `trace_filter` also writes a filtered
localization table next to the filtered trace table. This table contains only the
localizations whose `Buid` values are present as `Spot_ID` values in the final
filtered trace table, after all trace filters have run. The same filtered
localization table is used for the `*_localization_distribution_fluxes` plot, so
that the spot statistics describe the spots kept in the trace table.
