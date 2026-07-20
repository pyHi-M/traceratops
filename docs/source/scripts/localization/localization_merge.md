# localization_merge

**Reliability status**: `development`

```{eval-rst}
.. argparse::
   :ref: traceratops.localization_merge.parse_arguments
   :prog: localization_merge
```

## Output Files

The script generates:

1. `merged_localizations.ecsv`: Merged localization table in ECSV format. This is the default filename when `--output_file` is omitted.

2. `[output_folder]/[output_file]`: Merged localization table written to the destination selected with `--output_folder` and `--output_file`.

## Examples

```
$ cat file_list.txt | localization_merge [options]
$ ls *.ecsv | localization_merge [options]
$ find . -name "*.ecsv" | localization_merge [options]
```

1. Merge all .ecsv files in the current directory and save to the default output:

   ```$ ls *.ecsv | localization_merge```

2. Merge specific files and save with a custom name:

    `$ echo -e "file1.ecsv\nfile2.ecsv\nfile3.ecsv" | localization_merge -o combined.ecsv`

3. Merge files and save to a specific directory:

    `$ find ./data -name "loc_*.ecsv" | localization_merge -O ./results -o merged.ecsv`

4. Process files listed in a text file:

    `$ cat files_to_merge.txt | localization_merge`

## Notes

- The script requires the LocalizationTable class from the imageProcessing module
- Input files must be in a format readable by the LocalizationTable.load() method
- Output is saved in ECSV (Enhanced Character Separated Values) format
- The script reports the number of localizations in each file and the final merged file
