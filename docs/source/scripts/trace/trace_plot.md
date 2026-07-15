# trace_plot

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_plot.parse_arguments
   :prog: trace_plot
```


## Examples

### Export the first N traces

By default, `trace_plot` exports a limited number of traces from the input trace
file to PDB files. Use `--number_traces` (or `-n`) to choose how many traces to
export. For example, this exports only the first 100 traces to the default
`PDBs/` output folder:

```bash
$ trace_plot --input Trace_3D_barcode_KDtree_ROI:1.ecsv --number_traces 100
```

The output files are named after each exported `Trace_ID`.

```bash
$ trace_plot --input Trace_3D_barcode_KDtree_ROI:1.ecsv --number_traces 100 --output first_100_PDBs
```

this exports the first 100 traces to the folder `first_100_PDBs/`.

### Export one selected trace

```bash
$ ls Trace_3D_barcode_KDtree_ROI:1.ecsv | trace_plot --pipe --selected_trace 5b1e6f89-0362-4312-a7ed-fc55ae98a0a5
```

this pipes the file 'Trace_3D_barcode_KDtree_ROI:1.ecsv' into trace_plot and then selects a trace for conversion.

When `--selected_trace` is used, `--number_traces` is ignored.

### Export all traces

```bash
$ trace_plot --input Trace_3D_barcode_KDtree_ROI:1.ecsv --all
```

this exports all traces in the trace file. When `--all` is used,
`--number_traces` is ignored.



## Format for json dict

Please use the following format for the json dictionary to link barcode identities with different ATOM names in the PDB file:

```{"12": "C  ", "18": "C  ", "9": "P  "}```

keys provide barcode names in the trace file, these should be attributed to 3 character codes



## pymol

Some useful pymol commands:


```
set grid_mode,1
color green,  (name C*)
color red, (name P*)
```
