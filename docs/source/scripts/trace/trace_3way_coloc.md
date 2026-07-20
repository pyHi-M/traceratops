# trace_3way_coloc

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_3way_coloc.parse_arguments
   :prog: trace_3way_coloc
```

## Output Files

For each trace file analyzed, the script generates:

1. `[output]_anchor_[anchor].npy`: Mean three-way co-localization matrix for the selected anchor barcode.

2. `[output]_anchor_[anchor].[format]`: Heatmap of three-way co-localization frequencies between the anchor barcode and all possible pairs of other barcodes.

3. `[output]_anchor_[anchor]_sem.[format]`: Heatmap of the standard error of the mean (SEM) for the bootstrapped three-way co-localization frequencies.

## Examples


1. Analyze a single trace file with default parameters:
   ```bash
   trace_3way_coloc --input traces.ecsv --anchor 42
   ```

2. Analyze with custom distance cutoff and more bootstrap cycles:
   ```bash
   trace_3way_coloc --input traces.ecsv --anchor 42 --cutoff 0.25 --bootstrapping_cycles 100
   ```

## Notes

- The script requires the ChromatinTraceTable class from the matrixOperations module
- Input files must be in ECSV format compatible with ChromatinTraceTable.load()
- Bootstrapping is used to estimate mean and standard error of co-localization frequencies
- The distance cutoff is in micrometers (µm)
