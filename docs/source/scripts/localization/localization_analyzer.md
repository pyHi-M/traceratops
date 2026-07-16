# localization_analyzer

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.localization_analyzer.parse_arguments
   :prog: localization_analyzer
```

## What it does

`localization_analyzer` loads one localization table and generates a quality-control
summary plot for the detected barcode localizations. The plot includes:

- SNR distribution per barcode.
- SNR versus z-position scatter plot, colored by object class.
- Number of detections per barcode.
- Roundness versus skew scatter plot, colored by barcode.

The input table is read with `LocalizationTable.load`, so the script supports the
same localization formats as the localization table reader. Use `.ecsv` for
localization tables; legacy `.dat` files still load for now but emit a
deprecation warning because support will be discontinued in a future traceratops
release. `.4dn` files are also supported.

## Usage examples

Create the default PNG output, `localization_distribution_fluxes.png`:

```bash
localization_analyzer \
    --localization_file localizations_3D_barcode.ecsv
```

Choose an explicit output file:

```bash
localization_analyzer \
    --localization_file localizations_3D_barcode.ecsv \
    --output_file qc/localization_qc.png
```

Save the plot as SVG:

```bash
localization_analyzer \
    --localization_file localizations_3D_barcode.ecsv \
    --output_file qc/localization_qc \
    --format svg
```

## Notes

- `--localization_file` is required.
- `--format` accepts `png` or `svg` and defaults to `png`.
- When `--output_file` is omitted, the output is named
  `localization_distribution_fluxes.<format>`.
- If `--output_file` includes an extension, it is replaced by the selected
  `--format` value.
