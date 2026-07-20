# plot_matrix_comparison

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_matrix_comparison.parse_arguments
   :prog: plot_matrix_comparison
```

## Output Files

The script generates:

1. `[output]_scatter_plot.[format]`: Scatter plot comparing the values from the two input matrices and reporting the Pearson correlation.

2. `[output]_violin_plot.[format]`: Violin plot comparing the distributions of values from the two input matrices.

## Examples

```bash
plot_matrix_comparison --help
```

## Notes

- See the command-line reference above for the complete option list.
