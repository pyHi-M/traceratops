# trace_assign_mask

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.trace_assign_mask.parse_arguments
   :prog: trace_assign_mask
```

## Note

Use `trace_assign_mask` to assign specific *labels* to chromatin traces in a trace table.

`trace_assign_mask` will load a trace file and a number of NUMPY-formatted mask files and assign labels. If a trace falls within a mask, then the mask label will be assigned to the corresponding column of the trace table. If a trace falls *at the same time* within multiple masks, multiple labels will be appended to the corresponding column of the trace table. If a trace falls within no mask, then the label column of the trace table will be kept empty.

## Invoke

```bash
$ trace_assign_mask --input trace_file.ecsv --mask_file my_mask.npy --label mymask
```

This will apply the label `mymask` to traces falling within the masks of the file `my_mask.npy`. The output will be a trace file with the extension `labeled`.



Multiple mask files can be provided using piping.

```bash
$ ls my_traces*.ecsv | trace_assign_mask --mask_file my_mask.npy --pipe  --label mymask
```

In this case the `mymask` will be applied to multiple trace files.
