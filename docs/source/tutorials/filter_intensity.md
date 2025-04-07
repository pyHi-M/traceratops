# filtering trace tables by spot intensity

## copies localization files

Use a generic path (with `*`) to give a set of files:

```
localization_cp_files.py --files /user/DATA/exp_7/*/localizations.dat --destination_folder .
```

```{note}
You can check your generic path with:
`ls /user/DATA/exp_7/*/localizations.dat`
```

## merges localization files

```
ls *.dat | localization_merge.py
```

## filters repeated spots while keeping only the one with highest intensity

```
trace_filter.py --input merged_traces.ecsv --clean_spots --localization_file merged_localizations.ecsv
```

## Compare results

### Removing all repeated barcodes
![merged_traces_filtered_trace_statistics](https://github.com/user-attachments/assets/c355e091-d8d1-4194-9497-1fe65f6c5f64)
![merged_traces_filtered_barcode_detection](https://github.com/user-attachments/assets/fa6f0ab2-9ccb-45b5-a827-7715b8474921)

### Keeping only spot with highest intensity
![merged_traces_filtered_loc_trace_statistics](https://github.com/user-attachments/assets/4061a3bb-63b9-4f00-a166-b2e8427de8d9)
![merged_traces_filtered_loc_barcode_detection](https://github.com/user-attachments/assets/556e1695-030a-41e9-ab8a-91fb4ca6c876)
