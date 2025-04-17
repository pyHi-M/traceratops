# trace_3way_coloc

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_3way_coloc.parse_arguments
   :prog: trace_3way_coloc
```

## Usage

```bash
trace_3way_coloc --input TRACE_FILE.ecsv --anchor BARCODE_NUMBER [options]
cat file_list.txt | trace_3way_coloc --pipe --anchor BARCODE_NUMBER [options]
```

## Examples

1. Analyze a single trace file with default parameters:
   ```bash
   trace_3way_coloc --input traces.ecsv --anchor 42
   ```

2. Analyze with custom distance cutoff and more bootstrap cycles:
   ```bash
   trace_3way_coloc --input traces.ecsv --anchor 42 --cutoff 0.25 --bootstrapping_cycles 100
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
