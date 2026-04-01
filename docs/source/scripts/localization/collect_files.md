# collect_files

**Reliability status**: `stable`

```{eval-rst}
.. argparse::
   :ref: traceratops.collect_files.parse_arguments
   :prog: collect_files
```

## Why ?

When merging data from multiple ROIs, you need to collect one file per subdirectory
into a single folder. Two problems arise:

1. **Localization files** have the same name in every ROI (`localizations_3D_barcode.dat`),
   so they cannot be copied to the same folder without renaming.
2. **Trace files** have different names (the ROI number varies), but you need to select
   only the right files and ignore others (e.g. Matrix files).

`collect_files` solves both problems with a single tool using fixed-length pattern matching.

## Two Matching Modes

### Exact match (no `--variable-part`)

When all files share the same name, omit `--variable-part`. Files are automatically
renamed by inserting the subdirectory name before the extension:

```bash
collect_files --root data/RUT \
    --example-file "localizations_3D_barcode.dat" \
    --copy-to collected/
```

Result: `localizations_3D_barcode_013_ROI.dat`, `localizations_3D_barcode_014_ROI.dat`, etc.

### Variable match (with `--variable-part`)

When filenames differ by a fixed-length substring (e.g. the ROI number), specify
`--variable-part` so any characters of the same length are accepted at that position:

```bash
collect_files --root data/RUT \
    --example-file "Trace_3D_barcode_mask-mask0_ROI-13.ecsv" \
    --variable-part "13" \
    --copy-to collected/
```

This matches `ROI-14.ecsv`, `ROI-25.ecsv`, etc. but rejects `ROI-021.ecsv`
(different length) and `_Matrix_uniqueBarcodes.ecsv` (different total length).

## Notes

- Each immediate subdirectory of `--root` is scanned recursively for exactly one match.
- If a subdirectory has zero or multiple matches, the script stops with a clear error.
- Use `--force` to skip subdirectories with no match.
- Use `--rename` to insert the subdirectory name in the output filename (automatic in exact mode).
- Original file metadata (timestamps, permissions) is preserved.

## Replaces

This script replaces the former `localization_cp_files` script and the `find -exec cp`
pattern previously used for trace files.
