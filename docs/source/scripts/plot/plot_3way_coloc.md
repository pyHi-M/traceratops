# plot_3way_coloc

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_3way_coloc.parse_arguments
   :prog: plot_3way_coloc
```

## Output Files

- One heatmap image per input three-way co-localization matrix.
- The output is written next to each input `.npy` file using the selected `--output_format`.

## Examples

```bash
plot_3way_coloc --input path/to/matrix1.npy path/to/matrix2.npy
```

## Notes

- See the command-line reference above for the complete option list.
