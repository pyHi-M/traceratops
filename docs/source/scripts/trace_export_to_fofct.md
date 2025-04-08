# trace_export_to_fofct

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_export_to_fofct.parse_arguments
   :prog: trace_export_to_fofct
```

## Example BED file

The bed file should have no header.
```
chr2L	2343645	2356099	5
chr2L	2356147	2369783	9
chr2L	2369828	2381912	13
chr2L	2381947	2393854	17
chr2L	2393892	2405589	21
```

## Notes

For this first version, we don't consider, from the pyHiM trace table:
- the "mask_id" column: it's can be linked to the "Cell_ID" column but sometimes it's not the case, two "mask_id" can be linked to the same Cell_ID.
- the "label" column: usually it's linked to RNA species but it's like a global mask for many cells, so it's not RNA spots.

The output CSV file will have the following columns:
- Spot_ID
- Trace_ID
- X
- Y
- Z
- Chrom
- Chrom_Start
- Chrom_End
- Extra_Cell_ROI_ID ("ROI #" in the pyHiM trace table)

We need as run arguments:
- the path to the ECSV file
- the path to the BED file
- the path to the JSON file (optional)
- the path to the output CSV file (optional)

## Example

```
trace_export_to_fofct --ecsv_file /path/to/Trace_3D_barcode_KDtree_ROI-5.ecsv --bed_file /path/to/barcode.bed --json_file /path/to/parameters.json --output_file /path/to/output.csv
```
