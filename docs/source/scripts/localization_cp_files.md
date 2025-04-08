# localization_cp_files

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.localization_cp_files.parse_arguments
   :prog: localization_cp_files
```

## Why ?

At the moment the localization files form multiple ROIs have the same name. Thus we cannot copy them
onto the same folder to merge.

This script copies the files to a destination folder while renaming them to include
the unique folder name from their original path.

## Usage Examples

Basic usage with shell wildcard expansion:
```bash
localization_cp_files --files /path/to/*/data/file.dat --destination_folder /output/path
```

For instance:
```bash
localization_cp_files --files /home/marcnol/grey/ProcessedData_2024/Experiment_76_David_DNAFISH_HiM_HRKit_G1E_20240711/*/localize_3d/data/localizations_3D_barcode.dat --destination_folder .
```

With explicit file list:
```bash
localization_cp_files --files /path/to/001/data/file.dat /path/to/002/data/file.dat --destination_folder /output/path
```

## Notes
- The script automatically detects which part of the path varies (usually a folder number)
- Original file metadata (timestamps, permissions) is preserved
- Destination folder is created if it doesn't exist
