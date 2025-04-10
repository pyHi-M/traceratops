# trace_plot

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_plot.parse_arguments
   :prog: trace_plot
```


## Examples

```bash
$ ls Trace_3D_barcode_KDtree_ROI:1.ecsv | trace_plot --pipe --selected_trace 5b1e6f89-0362-4312-a7ed-fc55ae98a0a5
```

this pipes the file 'Trace_3D_barcode_KDtree_ROI:1.ecsv' into trace_plot and then selects a trace for conversion.

```bash
$ trace_plot --input Trace_3D_barcode_KDtree_ROI:1.ecsv --all
```

this plots all traces in the trace file.



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
