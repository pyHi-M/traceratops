# plot_bootstrapping

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.plot_bootstrapping.parse_arguments
   :prog: plot_bootstrapping
```

## Example

```
plot_bootstrapping --input Trace_Trace_all_ROIs_filtered_beta_mask0_exp_ND_1_PDX1LR_Matrix_PWDscMatrix.npy -U Trace_Trace_all_ROIs_filtered_exocrine_mask0_exp_ND_3_PDX1LR_Matrix_uniqueBarcodes.ecsv --N_bootstrap 100
```
