# figure_him_matrix

```{eval-rst}
.. argparse::
   :ref: traceratops.figure_him_matrix.parse_arguments
   :prog: figure_him_matrix
```


## Example

Here is examples usage of figure_him_matrix:

### Default (proximity)
```bash
figure_him_matrix -M PWDscMatrix.npy -B unique_barcodes.ecsv
```

<p align="center">
    <img src="../_static/Fig_PWDscMatrix_proximity_0.01-0.10.png" width="45%">
    <img src="../_static/Fig_PWDscMatrix_proximity_norm_0.82-0.95_nan.png" width="45%">
</p>


### Default normalized VS. KDE (c_map: Spectral)


```bash
figure_him_matrix -M PWDscMatrix.npy -B unique_barcodes.ecsv --remove_nan
```

**VS.**

```bash
figure_him_matrix -M PWDscMatrix.npy -B unique_barcodes.ecsv --mode KDE --c_map Spectral
```

<p align="center">
    <img src="../_static/Fig_PWDscMatrix_proximity_norm_0.20-0.59.png" width="45%">
    <img src="../_static/Fig_PWDscMatrix_KDE_0.21-0.37.png" width="45%">
</p>
