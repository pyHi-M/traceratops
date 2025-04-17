# process_him_matrix

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.process_him_matrix.parse_arguments
   :prog: process_him_matrix
```

## Output

sc_matrix_collated: 3D npy matrix. PWD matrix for single cells. Axes:0-1 barcodes, Axis:2, cellID
unique_barcodes: npy array. list of unique barcodes
SClabeledCollated: npy array. binary label indicating if cell is in pattern or not. Axis:0 cellID

## Note

This script performs the post-processing of one or more datasets previously analysed with *pyHiM*, defined in the `folders2Load.json` file.

It performs the following operations:
- Merges datasets from different experiments.
- Calculates and plots ensemble pairwise distance (PWD) matrix.
- Calculates and plots the inverse of the PWD matrix.
- Calculates and plots contact probability matrix for each dataset.
- Calculates and plots ensemble contact probability matrix.
- Calculates and plots tensemble 3-way contact probability matrix for the set of anchors defined in the `folders2Load.json` file.
- Optional: Reads MATLAB single-cell PWD matrices and performs all previous operations.
