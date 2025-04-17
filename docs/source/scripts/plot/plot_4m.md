# plot_4m

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_4m.parse_arguments
   :prog: plot_4m
```


## Usage

```bash
plot_4m --input TRACE_FILE.ecsv --anchor BARCODE_NUMBER [options]
cat file_list.txt | python plot_4m --pipe --anchor BARCODE_NUMBER [options]
find . -name "*.ecsv" | python plot_4m --pipe --anchor BARCODE_NUMBER [options]
```

## Examples

Analyze a single trace file with default parameters:

`plot_4m --input traces.ecsv --anchor 42`

Analyze with custom distance cutoff and more bootstrap cycles:

`plot_4m --input traces.ecsv --anchor 42 --cutoff 0.25 --bootstrapping_cycles 100`

Process multiple trace files in batch mode:

`cat trace_files.txt | plot_4m --pipe --anchor 42 --output batch_results.png`

Process all ECSV files in a directory:

`find ./data -name "*.ecsv" | plot_4m --pipe --anchor 42`

## Output

A PNG image with a plot showing colocalization frequencies between the anchor barcode
and all other barcodes, including error bars derived from bootstrapping
The output filename will be modified to include the anchor barcode number
(e.g., "colocalization_plot_anchor_42.png")

## Notes

The script requires the ChromatinTraceTable class from the matrixOperations module
Input files must be in ECSV format compatible with ChromatinTraceTable.load()
Bootstrapping is used to estimate mean and standard error of colocalization frequencies
The distance cutoff is in micrometers (µm)
Higher numbers of bootstrapping cycles increase statistical confidence but require more computation time
