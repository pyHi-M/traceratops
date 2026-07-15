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


## Visualizing traces in pymol

`trace_plot` generates a folder with PDB files containing the 3D coordinates of the traces. Each barcode is assigned a different ATOM name, which is then used in pymol to color code barcodes. Instead, if the user provided a json dictionary (see section below), these are used as ATOM names.

Once you generated the folder with the `PDB` files, it is time to load them in pymol.

For this you need to first install pymol. Then in the top box you will find a terminal where you can write text, starting with `PyMOL> `. In this terminal you will write the following commands:


```bash
run /home/user/Repositories/traceratops/traceratops/pymol_script.py

load_pdb_grid /home/user/data/outputs/PDBs
  
color_all_barcodes()
  
spin_grid
```

This assumes you have installed traceratops at `/home/user/Repositories/`. Otherwise select a different folder name.

This also assumes your PDBs are in folder `/home/user/data/outputs/PDBs`, if they are in a different location correct accordingly.

The first line will load the functions from a python script in traceratops. The second command will load 100 structures from the PDB file and represent them in a grid. The third line will color code each barcode with a different color. The final line will make the structures spin.

If you want to colorcode your barcodes using a user defined code, see next section.

## Format for json dict

Please use the following format for the json dictionary to link barcode identities with different ATOM names in the PDB file:

```{"12": "C  ", "18": "C  ", "9": "P  "}```

keys provide barcode names in the trace file, these should be attributed to 3 character codes. In this example, the barcode 12 is assigned the ATOM name `C`, barcode 18 is assigned `C` and atom 9 is assigned `P`. This provides a useful way of colorcoding your barcodes using a user defined coloring code.


In this case you can run the first two lines of the script above in python (e.g. run ... and load_pdb_grid ...) and then apply separate colors for each barcode as follows as follows:


```
set grid_mode,1
color green,  (name C*)
color red, (name P*)
```


