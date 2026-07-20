# trace_genomic_coordinates

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_genomic_coordinates.parse_arguments
   :prog: trace_genomic_coordinates
```

## Output Files

- An updated ECSV trace table containing genomic coordinate columns populated from the BED file.
- With `--output`, the table is written to that path; otherwise the input stem is reused with a genomic-coordinate suffix.

## Examples

```sh
trace_genomic_coordinates --input trace_file.ecsv --bed bed_file.bed --output output_file.ecsv
```

## Notes

### BED file format

```{warning}
**No header!**
```

### 4-column BED (default behavior)

```text
chrX            14785864        14789298                1
chrX            14789398        14792430                2
chrX            14792433        14795380                3
chrX            14795381        14798629                7
chrX            14799003        14802202                8
chrX            14802255        14805371                9
chrX            14805412        14809056                10
chrX            14809057        14812112                11
chrX            14812113        14814817                12
chrX            14814823        14818656                13
chrX            14824792        14829604                14
```

The last column should contain only numbers and should match the inputs in the Trace table which are themselves taken
from the filenames processed by pyHiM.

### 5-column BED (barcode substitution)

If you provide a fifth column, trace_genomic_coordinates will substitute barcode names in the Trace table: the 4th
column must match the barcode ID in the trace table, and the 5th column is the value that will be used after
substitution.

```text
chrX            14785864        14789298                1              HoxA1
chrX            14789398        14792430                2              HoxA2
chrX            14792433        14795380                3              HoxA3
chrX            14795381        14798629                7              HoxA7
chrX            14799003        14802202                8              HoxA8
chrX            14802255        14805371                9              HoxA9
chrX            14805412        14809056                10             HoxA10
chrX            14809057        14812112                11             HoxA11
chrX            14812113        14814817                12             HoxA12
chrX            14814823        14818656                13             HoxA13
chrX            14824792        14829604                14             HoxA14
```
