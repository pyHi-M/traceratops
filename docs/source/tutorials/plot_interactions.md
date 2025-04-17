# Plotting interactions from a trace table

## 4M plots

To make 4M plots for a cutoff distance of `0.20` microns, for anchor `18` using `10` bootstrapping cycles:

```
plot_4m --input  merged_traces_filtered.ecsv --anchor 18 --cutoff 0.20 --bootstrapping_cycles 10 --output coloc_plot.png
```

![coloc_plot_anchor_18](https://github.com/user-attachments/assets/cff0e957-e8d6-44ef-ab61-fe87c4bb7f0a)


Result:

## 3-way interactions

To plot the 3-way interactions for a cutoff distance of `0.20` microns, for anchor `25` using `2` bootstrapping cycles:

```
plot_3way --input merged_traces_filtered.ecsv --cutoff 0.20 --bootstrapping_cycles 2 --anchor 25
```

The results look like this:

![threeway_coloc_plot_anchor_26](https://github.com/user-attachments/assets/ad0cbc1a-3d5d-4f7e-920d-eb0fdaac9882)
