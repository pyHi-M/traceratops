# Using the `--pipe` argument

The `--pipe` flag lets traceratops commands read lists of trace files from standard input. This is handy when you want to reuse shell tools such as `find`, `ls`, or `cat` to build flexible file lists without typing them all on the command line. Commands accept one filename per line.

## Filter traces discovered with `find`

Use `find` to collect all trace tables in a directory tree, then pipe the list into `trace_filter` for batch processing:

```
find /data/experiment_42 -name "*.ecsv" \
  | trace_filter --pipe --n_barcodes 3 --output merged_traces_filtered.ecsv
```

This scans every `.ecsv` file under `/data/experiment_42`, filters them, and writes the combined output to `merged_traces_filtered.ecsv`.

## Plot interactions from a saved list

If you already have a text file listing the traces you want to analyze, feed it to plotting commands with `cat` and `--pipe`:

```
cat trace_files.txt | plot_4m --pipe --anchor 18 --cutoff 0.20 --output coloc_anchor18.png
```

Here `trace_files.txt` contains one trace path per line. The command produces a 4M interaction plot for anchor barcode `18` with a `0.20` cutoff and saves it as `coloc_anchor18.png`.

## Analyze traces with other tools

Most traceratops CLI tools that accept trace files support `--pipe`, so you can mix and match pipelines. For example, to run quality analysis across all traces in the current directory:

```
ls *.ecsv | trace_analyzer --pipe --output trace_summary.ecsv
```

This reads every listed trace, computes the analysis, and stores the summary in `trace_summary.ecsv`.
