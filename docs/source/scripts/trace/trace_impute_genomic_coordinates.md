# trace_impute_genomic_coordinates

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_impute_genomic_coordinates.parse_arguments
   :prog: trace_impute_genomic_coordinates
```


## Usage example

```sh
trace_impute_genomic_coordinates --input trace_file.ecsv --bed bed_file.bed --output output_file.ecsv
```

## BED file format

```{warning}
**No header!**
```

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
from the filenames processed by pyHiM
