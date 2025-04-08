# trace_analyser

**Reliability status**: `development`

*Read a trace table and derive several different analysis, including the distribution of number barcodes per trace: total, duplicated, and non-duplicated. It also exports the Rg distribution.*


```{eval-rst}
.. argparse::
   :ref: traceratops.trace_assign_mask.parse_arguments
   :prog: trace_assign_mask
```

## Description
A Python script for analyzing chromatin trace files. The script loads trace files and generates statistical analysis and visualizations of various properties:
- Number of barcodes detected per trace
- Number of duplicated barcodes
- Distribution of localizations in X, Y, Z coordinates
- Distances between consecutive neighboring barcodes
- Visual representations of traces
- Barcode detection efficiency analysis with bootstrapped error estimates

## Usage

```bash
# Analyze a single trace file
trace_analyzer --input path/to/your/trace_file.ecsv

# Analyze multiple trace files using pipe
find /path/to/traces -name "*.ecsv" | trace_analyzer --pipe

# Specify a different root folder
trace_analyzer --rootFolder /custom/data/folder --input trace_file.ecsv

# Generate SVG output instead of PNG
trace_analyzer --input trace_file.ecsv --format svg

# Generate XYZ plots in addition to statistics
trace_analyzer --input trace_file.ecsv --plotXYZ
```

## Output Files

For each trace file analyzed, the script generates:

1. `[tracefile]_trace_statistics.[format]`: Histograms of barcode statistics
2. `[tracefile]_first_neighbor_distances.[format]`: Histograms of distances between consecutive neighboring barcodes
3. `[tracefile]_barcode_detection.[format]`: Plot of barcode detection efficiency
4. `[tracefile]_relative_barcode_frequencies`: Barcode statistics file
5. `[tracefile]_traces_XYZ.[format]`: Visual representation of the traces (if --plotXYZ is set)

![Trace_3D_barcode_mask:DAPI_ROI:3_filtered_traces_XYZ_ROI3](https://github.com/pyHi-M/pyHiM/assets/341757/2b3f32f2-d9a6-41c8-98b7-372cc60a0439)

![Trace_3D_barcode_mask:DAPI_ROI:3_filtered_trace_statistics](https://github.com/pyHi-M/pyHiM/assets/341757/281cf895-d043-422c-a7c1-5fc7dcbbf857)

![Trace_3D_barcode_mask:DAPI_ROI:3_filtered_xyz_statistics](https://github.com/pyHi-M/pyHiM/assets/341757/c33e6a4b-0678-4f1e-b1cd-8834e0779560)
