# plot_3way_coloc

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_3way_coloc.parse_arguments
   :prog: plot_3way_coloc
```

## Output Files

For each input matrix, the script generates:

1. `[matrixfile]_anchor_[anchor]_replot.[format]`: Heatmap replot of the three-way co-localization matrix. The anchor value is extracted from the input filename.

## Examples

```bash
plot_3way_coloc --input path/to/matrix1.npy path/to/matrix2.npy
```

## Notes

- See the command-line reference above for the complete option list.
