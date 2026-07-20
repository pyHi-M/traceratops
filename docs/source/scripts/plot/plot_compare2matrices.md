# plot_compare2matrices

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_compare2matrices.parse_arguments
   :prog: plot_compare2matrices
```

## Output Files

The script generates:

1. `[outputFolder]/Fig_[input1]_[mode]_[matrix_norm_mode]_[cMax].[format]`: Matrix plot for the first input matrix.

2. `[outputFolder]/Fig_[input2]_[mode]_[matrix_norm_mode]_[cMax].[format]`: Matrix plot for the second input matrix.

3. `[outputFolder]/Fig_[input1]_difference.[format]`: Difference or ratio plot comparing the two matrices.

4. `[outputFolder]/Fig_[input1]_mixed_matrix.[format]`: Mixed matrix plot with one input in the upper triangle and the other in the lower triangle.

5. `[outputFolder]/Fig_[input1]_mixed_matrix.npy`: NPY file containing the mixed matrix values.

6. `[outputFolder]/Fig_[input1]_Wilcoxon.[format]`: Wilcoxon statistical comparison plot. This file is written for non-proximity distance modes.

7. `[outputFolder]/Fig_[input1]_Wilcoxon.npy`: NPY file containing the Wilcoxon statistical comparison values. This file is written for non-proximity distance modes.

## Examples

```bash
plot_compare2matrices --help
```

## Notes

- See the command-line reference above for the complete option list.
