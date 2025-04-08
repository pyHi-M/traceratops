# trace_import_from_fofct

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_import_from_fofct.parse_arguments
   :prog: trace_import_from_fofct
```

## Usage Example:
To convert a CSV file back to the ECSV format using the specified BED and JSON files, you would run:

```
trace_import_from_fofct --fofct_file output.csv --bed_file barcode.bed --output_file Trace_3D_barcode_KDtree_ROI-5.ecsv
```

If the `--output_file` argument is not provided, the script will save the ECSV file with the same name as the input CSV file but with an `.ecsv` extension.
