# plot_3way

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_3way.parse_arguments
   :prog: plot_3way
```

## Usage

```bash
plot_3way --input TRACE_FILE.ecsv --anchor BARCODE_NUMBER [options]
cat file_list.txt | plot_3way --pipe --anchor BARCODE_NUMBER [options]
```

## Examples

1. Analyze a single trace file with default parameters:
   ```bash
   plot_3way --input traces.ecsv --anchor 42
   ```

2. Analyze with custom distance cutoff and more bootstrap cycles:
   ```bash
   plot_3way --input traces.ecsv --anchor 42 --cutoff 0.25 --bootstrapping_cycles 100
   ```

## Output

- A PNG image with a heatmap showing three-way co-localization frequencies between the
  anchor barcode and all possible pairs of other barcodes
- The output filename will include the anchor barcode number
  (e.g., "threeway_coloc_plot_anchor_42.png")

## Notes

- The script requires the ChromatinTraceTable class from the matrixOperations module
- Input files must be in ECSV format compatible with ChromatinTraceTable.load()
- Bootstrapping is used to estimate mean and standard error of co-localization frequencies
- The distance cutoff is in micrometers (µm)
