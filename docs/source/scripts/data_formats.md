# Data formats

Traceratops reads and writes tabular data as [Astropy ECSV](https://docs.astropy.org/en/stable/io/ascii/ecsv.html) files. ECSV files are plain-text tables with a YAML metadata header, followed by a whitespace-delimited data table.

A typical ECSV file has these parts:

- `# %ECSV 1.0`: identifies the file as ECSV version 1.0.
- `# ---`: starts the YAML header.
- `# datatype`: lists every column name and data type.
- `# meta`: optional table-level metadata. Traceratops commonly stores processing comments here, such as filtering steps, coordinate units, genome assembly, or merge history.
- `# schema: astropy-2.0`: identifies the Astropy ECSV schema.
- The first non-comment line is the table header. Column names containing spaces or special characters are quoted, for example `"ROI #"`.
- Remaining lines are one row per spot, trace localization, or detected object.

## Trace table

A trace table stores spots that have been assigned to chromatin traces. Each row is one spot in one trace.

### Example header

```text
# %ECSV 1.0
# ---
# datatype:
# - {name: Spot_ID, datatype: string}
# - {name: Trace_ID, datatype: string}
# - {name: x, datatype: float32}
# - {name: y, datatype: float32}
# - {name: z, datatype: float32}
# - {name: Chrom, datatype: string}
# - {name: Chrom_Start, datatype: int64}
# - {name: Chrom_End, datatype: int64}
# - {name: 'ROI #', datatype: int64}
# - {name: Mask_id, datatype: int64}
# - {name: 'Barcode #', datatype: int64}
# - {name: label, datatype: string}
# meta: !!omap
# - comments: [xyz_unit=micron, genome_assembly=mm10]
# schema: astropy-2.0
Spot_ID Trace_ID x y z Chrom Chrom_Start Chrom_End "ROI #" Mask_id "Barcode #" label
```

### Columns

| Column | Type | Meaning |
| --- | --- | --- |
| `Spot_ID` | `string` | Unique identifier for the spot/localization. This usually matches the localization table `Buid` value when the trace was built from a localization table. |
| `Trace_ID` | `string` | Unique identifier for the chromatin trace. Rows with the same `Trace_ID` belong to the same reconstructed trace. |
| `x` | `float32` | Spot x coordinate. In traceratops trace tables this is normally in micrometers when the metadata contains `xyz_unit=micron`. |
| `y` | `float32` | Spot y coordinate, normally in the same unit as `x`. |
| `z` | `float32` | Spot z coordinate, normally in the same unit as `x`. |
| `Chrom` | `string` | Chromosome or contig associated with the barcode, for example `chr2L`. This can be imputed from a BED file. |
| `Chrom_Start` | `int64` | Genomic start coordinate for the barcode interval. |
| `Chrom_End` | `int64` | Genomic end coordinate for the barcode interval. |
| `ROI #` | `int64` | Region-of-interest identifier from the acquisition or processing workflow. |
| `Mask_id` | `int64` | Segmentation mask identifier associated with the trace or spot. Depending on the dataset, this often represents the cell or nuclear mask used during trace construction. |
| `Barcode #` | `int64` | Barcode identifier for the genomic locus or probe round. This value is used to order spots along the assayed genomic region and to match BED barcode definitions. |
| `label` | `string` | Optional label assigned to the trace or spot, for example after mask assignment or label splitting. It may be empty when no label has been assigned. |

### Metadata comments

The `meta.comments` list records table-level processing information. Common entries include:

- `xyz_unit=micron`: coordinate units for `x`, `y`, and `z`.
- `genome_assembly=...`: genome assembly used for genomic coordinates, for example `mm10`.
- `filt:...`: filtering operations applied to the table.
- `Genomic coordinates imputed from BED file. ... rows matched.`: genomic coordinate annotation history.
- `appended_trace_files=...`: number of trace files merged into the table.

## Localization table

A localization table stores detected spots before or alongside trace assignment. Each row is one detected localization and its quality-control measurements.

### Example header

```text
# %ECSV 1.0
# ---
# datatype:
# - {name: Buid, datatype: string}
# - {name: 'ROI #', datatype: int64}
# - {name: 'CellID #', datatype: int64}
# - {name: 'Barcode #', datatype: int64}
# - {name: id, datatype: int64}
# - {name: zcentroid, datatype: float32}
# - {name: xcentroid, datatype: float32}
# - {name: ycentroid, datatype: float32}
# - {name: snr, datatype: float32}
# - {name: spot_pixel_percentage, datatype: float32}
# - {name: skew, datatype: float32}
# - {name: patch_size, datatype: int64}
# - {name: object_class, datatype: int64}
# - {name: mean_intensity, datatype: float32}
# - {name: flux, datatype: float32}
# - {name: roundness, datatype: float32}
# meta: !!omap
# - comments: [filtered, filtered]
# schema: astropy-2.0
Buid "ROI #" "CellID #" "Barcode #" id zcentroid xcentroid ycentroid snr spot_pixel_percentage skew patch_size object_class mean_intensity flux roundness
```

### Columns

| Column | Type | Meaning |
| --- | --- | --- |
| `Buid` | `string` | Unique barcode localization identifier. This is the localization-table identifier matched to `Spot_ID` in trace tables. |
| `ROI #` | `int64` | Region-of-interest identifier from the acquisition or processing workflow. |
| `CellID #` | `int64` | Cell identifier assigned during segmentation or localization processing. |
| `Barcode #` | `int64` | Barcode identifier for the genomic locus or probe round. |
| `id` | `int64` | Row-local object identifier for the detected localization within the source table or barcode image. |
| `zcentroid` | `float32` | Detected spot centroid in z, typically in image-coordinate units from the localization pipeline. |
| `xcentroid` | `float32` | Detected spot centroid in x, typically in image-coordinate units from the localization pipeline. |
| `ycentroid` | `float32` | Detected spot centroid in y, typically in image-coordinate units from the localization pipeline. |
| `snr` | `float32` | Signal-to-noise ratio for the detected spot. Higher values indicate stronger signal relative to noise. |
| `spot_pixel_percentage` | `float32` | Percentage of pixels in the detection patch assigned to the spot. This can be used to remove overly broad or poorly segmented detections. |
| `skew` | `float32` | Intensity-distribution skewness measurement for the detected spot. |
| `patch_size` | `int64` | Size of the image patch, in pixels, used to measure or classify the object. |
| `object_class` | `int64` | Classification label for the detected object. In traceratops filtering, `object_class >= 1` is commonly used to keep spot-like objects and remove background class `0`. |
| `mean_intensity` | `float32` | Mean intensity measured for the detected localization. |
| `flux` | `float32` | Integrated or fitted spot signal reported by the localization workflow. In some outputs this can match the SNR-like score used by upstream detection. |
| `roundness` | `float32` | Shape descriptor for the detected spot. Values are used by quality filters to remove objects with undesired shape. |

### Metadata comments

The `meta.comments` list records table-level processing information. For localization tables, entries such as `filtered` indicate that one or more filtering steps were applied before writing the file.

## Relationship between trace and localization tables

Trace tables and localization tables can be used together for quality filtering and provenance:

- `trace.Spot_ID` corresponds to `localization.Buid`.
- `ROI #` and `Barcode #` appear in both formats and help match spots to acquisition regions and barcode rounds.
- Trace coordinates are stored as `x`, `y`, and `z`, while localization centroids are stored as `xcentroid`, `ycentroid`, and `zcentroid`.
- Localization quality columns such as `snr`, `mean_intensity`, `spot_pixel_percentage`, `skew`, `patch_size`, `object_class`, and `roundness` can be used to filter spots, which in turn filters the corresponding trace-table rows.
