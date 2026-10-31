# plot_bootstrapping

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_bootstrapping.parse_arguments
   :prog: plot_bootstrapping
```

## Output Files

The script generates:

1. `[output]/Fig_[matrixfile]bootstrapping_median.npy`: Bootstrapped median matrix values.

2. `[output]/Fig_[matrixfile]bootstrapping_median.[format]`: Heatmap of the bootstrapped median matrix.

3. `[output]/Fig_[matrixfile]bootstrapping_std_median.npy`: Standard-error matrix values from bootstrapping.

4. `[output]/Fig_[matrixfile]bootstrapping_std_median.[format]`: Heatmap of the standard-error matrix.

## Examples

```
plot_bootstrapping --input Trace_Trace_all_ROIs_filtered_beta_mask0_exp_ND_1_PDX1LR_Matrix_PWDscMatrix.npy --barcodes Trace_Trace_all_ROIs_filtered_exocrine_mask0_exp_ND_3_PDX1LR_Matrix_uniqueBarcodes.ecsv --N_bootstrap 100
```

## Notes

- See the command-line reference above for the complete option list.
