# Source inventory

Snapshot: dev `d767666dfff6be0eca410cca1cc6acc79970140e`.

29 registered command names; 25 implementations. All registered modules expose `parse_arguments()` returning argparse parsers.

| Command | Entry point | Documentation | Tests |
|---|---|---|---|
| trace_filter | `traceratops.trace_filter:main` | docs/source/scripts/trace/trace_filter.md | test/test_trace_filter.py |
| trace_merge | `traceratops.trace_merge:main` | docs/source/scripts/trace/trace_merge.md | test/test_trace_merge.py |
| figure_him_matrix | `traceratops.plot_him_matrix:main` | docs/source/scripts/plot/plot_him_matrix.md | test/test_plot_him_matrix.py |
| plot_him_matrix | `traceratops.plot_him_matrix:main` | docs/source/scripts/plot/plot_him_matrix.md | test/test_plot_him_matrix.py |
| trace_to_matrix | `traceratops.trace_to_matrix:main` | docs/source/scripts/trace/trace_to_matrix.md | test/test_trace_to_matrix.py |
| trace_stats | `traceratops.trace_stats:main` | docs/source/scripts/trace/trace_stats.md | No dedicated command test |
| trace_splitter | `traceratops.trace_splitter:main` | docs/source/scripts/trace/trace_splitter.md | test/test_trace_splitter.py |
| trace_sorter | `traceratops.trace_sorter:main` | docs/source/scripts/trace/trace_sorter.md | test/test_trace_sorter.py |
| trace_plot | `traceratops.trace_plot:main` | docs/source/scripts/trace/trace_plot.md | test/test_trace_plot.py |
| trace_pearsons | `traceratops.trace_pearsons:main` | docs/source/scripts/trace/trace_pearsons.md | No dedicated command test |
| trace_genomic_coordinates | `traceratops.trace_genomic_coordinates:main` | docs/source/scripts/trace/trace_genomic_coordinates.md | No dedicated command test |
| trace_physical_vs_genomic_distance | `traceratops.trace_physical_vs_genomic_distance:main` | docs/source/scripts/trace/trace_physical_vs_genomic_distance.md | test/test_trace_physical_vs_genomic_distance.py |
| trace_import_from_fofct | `traceratops.trace_import_from_fofct:main` | docs/source/scripts/trace/trace_import_from_fofct.md | No dedicated command test |
| trace_export_to_fofct | `traceratops.trace_export_to_fofct:main` | docs/source/scripts/trace/trace_export_to_fofct.md | No dedicated command test |
| trace_assign_mask | `traceratops.trace_assign_mask:main` | docs/source/scripts/trace/trace_assign_mask.md | No dedicated command test |
| trace_analyzer | `traceratops.trace_analyzer:main` | docs/source/scripts/trace/trace_analyzer.md | No dedicated command test |
| localization_merge | `traceratops.localization_merge:main` | docs/source/scripts/localization/localization_merge.md | No dedicated command test |
| localization_analyzer | `traceratops.localization_analyzer:main` | docs/source/scripts/localization/localization_analyzer.md | No dedicated command test |
| compare_pwd_matrices | `traceratops.plot_matrix_comparison:main` | docs/source/scripts/plot/plot_matrix_comparison.md | test/test_plot_matrix_comparison.py |
| plot_matrix_comparison | `traceratops.plot_matrix_comparison:main` | docs/source/scripts/plot/plot_matrix_comparison.md | test/test_plot_matrix_comparison.py |
| replot_3way | `traceratops.plot_3way_coloc:main` | docs/source/scripts/plot/plot_3way_coloc.md | No dedicated command test |
| plot_3way_coloc | `traceratops.plot_3way_coloc:main` | docs/source/scripts/plot/plot_3way_coloc.md | No dedicated command test |
| plot_bootstrapping | `traceratops.plot_bootstrapping:main` | docs/source/scripts/plot/plot_bootstrapping.md | test/test_plot_bootstrapping.py |
| plot_4m | `traceratops.plot_4m:main` | docs/source/scripts/plot/plot_4m.md | No dedicated command test |
| plot_3way | `traceratops.trace_3way_coloc:main` | docs/source/scripts/trace/trace_3way_coloc.md | test/test_trace_3way_coloc.py |
| trace_3way_coloc | `traceratops.trace_3way_coloc:main` | docs/source/scripts/trace/trace_3way_coloc.md | test/test_trace_3way_coloc.py |
| plot_compare2matrices | `traceratops.plot_compare2matrices:main` | docs/source/scripts/plot/plot_compare2matrices.md | test/test_plot_compare2matrices_args.py |
| collect_files | `traceratops.collect_files:main` | docs/source/scripts/localization/collect_files.md | test/test_collect_files.py |
| trace_split_labels | `traceratops.trace_split_labels:main` | No dedicated page | No dedicated command test |

## Python surface and inclusion policy

The package root exports no analysis API. Tutorials primarily invoke CLI commands. Source docstrings document public analytical functions and methods in core modules and CLI modules; trace_resolver additionally defines an explicit __all__. The catalog includes source-documented analytical callables and public analytical methods on exported resolver classes. Source documentation is not a stability guarantee. CLI launchers, parser factories, filename/report utilities, I/O-only methods, private names, and internal resolver submodule helpers are excluded as recommended interfaces. The inventory below records candidates, including exclusions, so additions can be reviewed.

| Source | Callable | Docstring summary |
|---|---|---|
| traceratops/collect_files.py:60 | `FilePattern.total_length` |  |
| traceratops/collect_files.py:63 | `FilePattern.matches` | Return True if *filename* matches this pattern. |
| traceratops/collect_files.py:81 | `build_pattern` | Derive a :class:`FilePattern` from *example_file*. |
| traceratops/collect_files.py:121 | `scan_subdirectories` | Scan immediate subdirectories of *root* for files matching *pattern*. |
| traceratops/collect_files.py:144 | `dest_name` | Compute the destination filename for a collected file. |
| traceratops/collect_files.py:155 | `check_collisions` | Return human-readable collision descriptions, if any. |
| traceratops/collect_files.py:185 | `format_report` | Build a deterministic, sorted, human-readable scan report. |
| traceratops/collect_files.py:210 | `collect_files` | Collect one file per subdirectory and copy to *copy_to*. |
| traceratops/collect_files.py:308 | `parse_arguments` | Build and return the argument parser (used by sphinx-argparse). |
| traceratops/collect_files.py:362 | `main` | Entry point.  Returns 0 on success, 1 on business-logic failure. |
| traceratops/core/build_matrix.py:64 | `BuildMatrix.initialize_parameters` |  |
| traceratops/core/build_matrix.py:81 | `BuildMatrix.calculate_pwd_single_mask` | Calculates PWD between barcodes detected in a given mask. For this: |
| traceratops/core/build_matrix.py:178 | `BuildMatrix.build_distance_matrix` | Builds pairwise distance matrix from a coordinates table |
| traceratops/core/build_matrix.py:255 | `BuildMatrix.calculate_n_matrix` |  |
| traceratops/core/build_matrix.py:266 | `BuildMatrix.get_output_prefix` | Build the output filename prefix for a trace file. |
| traceratops/core/build_matrix.py:289 | `BuildMatrix.plots_all_matrices` | Plots all matrices after analysis |
| traceratops/core/build_matrix.py:418 | `BuildMatrix.save_matrices` |  |
| traceratops/core/build_matrix.py:435 | `BuildMatrix.launch_analysis` | run analysis for a chromatin trace table. |
| traceratops/core/build_matrix.py:486 | `BuildMatrix.run` |  |
| traceratops/core/chromatin_trace_table.py:27 | `read_table_from_ecsv` | Read an astropy Table saved as an ``ecsv`` file. |
| traceratops/core/chromatin_trace_table.py:35 | `save_table_to_ecsv` | Save an astropy table into an ``ecsv`` file. |
| traceratops/core/chromatin_trace_table.py:47 | `generate_pyhim_identifier` | Return a pyHiM-style UUID identifier. |
| traceratops/core/chromatin_trace_table.py:52 | `relabel_pyhim_identifiers` | Relabel 4DN identifiers to pyHiM Astropy UUID nomenclature. |
| traceratops/core/chromatin_trace_table.py:69 | `micrometer_scale_from_xyz_unit` | Return the factor needed to convert an XYZ unit to micrometers. |
| traceratops/core/chromatin_trace_table.py:86 | `random_label_cmap` | Generates a random colormap similar to stardist (so you don't have to import this library just to do it). |
| traceratops/core/chromatin_trace_table.py:119 | `ChromatinTraceTable.initialize` |  |
| traceratops/core/chromatin_trace_table.py:156 | `ChromatinTraceTable.load` | Loads a trace table from a .ecsv or .4dn file. |
| traceratops/core/chromatin_trace_table.py:182 | `ChromatinTraceTable.remove_empty_comments` |  |
| traceratops/core/chromatin_trace_table.py:191 | `ChromatinTraceTable.remove_duplicate_comments` |  |
| traceratops/core/chromatin_trace_table.py:200 | `ChromatinTraceTable.save` | Saves the trace table in the appropriate format (.ecsv or .4dn). |
| traceratops/core/chromatin_trace_table.py:373 | `ChromatinTraceTable.prevent_roi_conflict` |  |
| traceratops/core/chromatin_trace_table.py:387 | `ChromatinTraceTable.load_bed_file` | Loads the BED file into a dictionary mapping barcode numbers to genomic coordinates. |
| traceratops/core/chromatin_trace_table.py:431 | `ChromatinTraceTable.impute_genomic_coordinates` | Updates the Chrom, Chrom_Start, and Chrom_End columns in the trace file based on the BED file. |
| traceratops/core/chromatin_trace_table.py:482 | `ChromatinTraceTable.append` | appends <table> to self.data |
| traceratops/core/chromatin_trace_table.py:499 | `ChromatinTraceTable.filter_traces_by_coordinate` | This function will remove the spots that are outside coordinate limits |
| traceratops/core/chromatin_trace_table.py:554 | `ChromatinTraceTable.filter_by_intensity` | Filters localizations in the trace file based on intensity from the localization table. |
| traceratops/core/chromatin_trace_table.py:562 | `ChromatinTraceTable.filter_by_localization_metrics` | Filter trace rows by one or more localization-table quality thresholds. |
| traceratops/core/chromatin_trace_table.py:687 | `ChromatinTraceTable.barcode_statistics` | calculates the number of times a barcode is repeated in a trace for all traces in a trace table |
| traceratops/core/chromatin_trace_table.py:726 | `ChromatinTraceTable.plots_barcode_statistics` | plots the collecive_bracode stats (see previous function) |
| traceratops/core/chromatin_trace_table.py:826 | `ChromatinTraceTable.trace_remove_label` | This function will remove traces that do not contain the word 'label' in the 'label' column |
| traceratops/core/chromatin_trace_table.py:857 | `ChromatinTraceTable.trace_keep_label` | This function will remove traces that do not contain the word 'label' in the 'label' column |
| traceratops/core/chromatin_trace_table.py:888 | `ChromatinTraceTable.filter_repeated_barcodes` | This function will remove the barcodes that are present more than once in a trace. |
| traceratops/core/chromatin_trace_table.py:983 | `ChromatinTraceTable.remove_duplicates_loc` | Removes duplicated barcodes within each trace. |
| traceratops/core/chromatin_trace_table.py:1068 | `ChromatinTraceTable.remove_duplicates` | removes duplicated (identical) spots |
| traceratops/core/chromatin_trace_table.py:1122 | `ChromatinTraceTable.remove_barcode` | Removes a specific barcode from a trace table |
| traceratops/core/chromatin_trace_table.py:1171 | `ChromatinTraceTable.filter_traces_by_n` | Removes rows in trace table with less than `minimum_number_barcodes` barcodes |
| traceratops/core/chromatin_trace_table.py:1239 | `ChromatinTraceTable.plots_traces` | This function plots 3 subplots (xy, xz, yz) with the localizations. |
| traceratops/core/him_matrix_operations.py:35 | `AnalysisHiMMatrix.load_data` | loads dataset |
| traceratops/core/him_matrix_operations.py:117 | `AnalysisHiMMatrix.plot_2d_matrix_simple` |  |
| traceratops/core/him_matrix_operations.py:180 | `AnalysisHiMMatrix.n_cells_loaded` |  |
| traceratops/core/him_matrix_operations.py:196 | `AnalysisHiMMatrix.retrieve_sc_matrix` | retrieves single cells that have the label requested |
| traceratops/core/him_matrix_operations.py:232 | `load_list` |  |
| traceratops/core/him_matrix_operations.py:244 | `list_sc_to_keep` |  |
| traceratops/core/him_matrix_operations.py:266 | `plot_matrix` |  |
| traceratops/core/him_matrix_operations.py:376 | `shuffle_matrix` |  |
| traceratops/core/him_matrix_operations.py:394 | `decodes_trace` | from a trace entry, provides Numpy array with coordinates, barcode and trace names |
| traceratops/core/him_matrix_operations.py:421 | `write_xyz_2_pdb` | Write a chromatin trace to a PDB file using pseudoatoms. |
| traceratops/core/him_matrix_operations.py:567 | `write_xyz_2_pdb_legacy` |  |
| traceratops/core/him_matrix_operations.py:703 | `distances_2_coordinates` | Infer coordinates from distances |
| traceratops/core/him_matrix_operations.py:749 | `is_notebook` | This function detects if you are running on an ipython console or in the shell. |
| traceratops/core/him_matrix_operations.py:774 | `is_nan_symmetric_matrix` | Return True when a 3D single-cell matrix is symmetric, treating NaNs as equal. |
| traceratops/core/him_matrix_operations.py:779 | `plot_distance_histograms` |  |
| traceratops/core/him_matrix_operations.py:869 | `adjust_colorbar` | Adjusts the colorbar to display correct limits and tick values. |
| traceratops/core/him_matrix_operations.py:922 | `plot_single_matrix` |  |
| traceratops/core/him_matrix_operations.py:956 | `get_matrix_title` |  |
| traceratops/core/him_matrix_operations.py:975 | `adjust_cmin_cmax` |  |
| traceratops/core/him_matrix_operations.py:985 | `get_plot_path` |  |
| traceratops/core/him_matrix_operations.py:1004 | `plot_nan_matrix` |  |
| traceratops/core/him_matrix_operations.py:1040 | `plot_him_matrix` |  |
| traceratops/core/him_matrix_operations.py:1090 | `calculate_nan_matrix` |  |
| traceratops/core/him_matrix_operations.py:1104 | `calculate_contact_probability_matrix` |  |
| traceratops/core/him_matrix_operations.py:1145 | `find_optimal_kernel_width` |  |
| traceratops/core/him_matrix_operations.py:1155 | `retrieve_kernel_density_estimator` | Gets the kernel density function and maximum from a distribution of PWD distances |
| traceratops/core/him_matrix_operations.py:1206 | `distribution_maximum_kernel_density_estimation` | calculates the kernel distribution and its maximum from a set of PWD distances |
| traceratops/core/him_matrix_operations.py:1280 | `get_rg_from_pwd` | Calculates the Rg from a 2D pairwise distance matrix |
| traceratops/core/him_matrix_operations.py:1324 | `get_detection_eff_barcodes` | Return the detection efficiency of all barcodes. |
| traceratops/core/him_matrix_operations.py:1361 | `get_barcodes_per_cell` | Returns the number of barcodes that were detected in each cell of sc_matrix_collated. |
| traceratops/core/him_matrix_operations.py:1377 | `get_coordinates_from_pwd_matrix` |  |
| traceratops/core/him_matrix_operations.py:1398 | `sort_cells_by_number_pwd` |  |
| traceratops/core/him_matrix_operations.py:1421 | `kde_fit` |  |
| traceratops/core/him_matrix_operations.py:1430 | `calculate_ensemble_pwd_matrix` | performs a KDE or median to calculate the max of the PWD distribution |
| traceratops/core/io_manager.py:10 | `create_folder` | Create folder with `makedirs` from os module. |
| traceratops/core/io_manager.py:27 | `load_barcode_dict` | Loads a barcode type dict JSON |
| traceratops/core/io_manager.py:55 | `save_json` | Save a python dict as a JSON file |
| traceratops/core/localization_table.py:24 | `warn_legacy_dat_extension` | Warn when a legacy ``.dat`` localization table is used as input. |
| traceratops/core/localization_table.py:30 | `normalize_localization_output_path` | Return an ECSV output path for non-4DN localization tables. |
| traceratops/core/localization_table.py:192 | `LocalizationTable.load` | Loads barcode_map |
| traceratops/core/localization_table.py:240 | `LocalizationTable.remove_empty_comments` |  |
| traceratops/core/localization_table.py:249 | `LocalizationTable.remove_duplicate_comments` |  |
| traceratops/core/localization_table.py:253 | `LocalizationTable.save` | Saves output table |
| traceratops/core/localization_table.py:294 | `LocalizationTable.append` | appends <table> to self.data |
| traceratops/core/localization_table.py:311 | `LocalizationTable.plot_distribution_fluxes` | This function will plot: |
| traceratops/core/localization_table.py:493 | `LocalizationTable.plot_intensity_distribution` | Plots the distribution of intensity values from the localization table |
| traceratops/core/localization_table.py:512 | `LocalizationTable.plots_localizations` | This function plots 3 subplots (xy, xz, yz) with the localizations. |
| traceratops/core/localization_table.py:573 | `LocalizationTable.compares_localizations` | Compares the localizations of two barcode tables |
| traceratops/core/localization_table.py:638 | `decode_rois` |  |
| traceratops/core/localization_table.py:648 | `build_color_dict` |  |
| traceratops/core/localization_table.py:657 | `plots_localization_projection` | This function will produce the scatter plot and add title |
| traceratops/core/localization_table.py:683 | `read_table_from_ecsv` | Read an astropy Table saved as an ``ecsv`` file. |
| traceratops/core/localization_table.py:691 | `create_output_table` |  |
| traceratops/core/plotting_functions.py:21 | `bootstrapping` |  |
| traceratops/core/plotting_functions.py:36 | `bootstraps_matrix` |  |
| traceratops/core/plotting_functions.py:66 | `get_position_from_unique_barcode_list` |  |
| traceratops/core/plotting_functions.py:70 | `gets_matrix` |  |
| traceratops/core/plotting_functions.py:138 | `Wilcoxon_matrix` |  |
| traceratops/core/plotting_functions.py:158 | `plot_Wilcoxon_matrix` |  |
| traceratops/core/plotting_functions.py:209 | `normalize_matrix` |  |
| traceratops/core/plotting_functions.py:235 | `plot_2d_matrix_simple` |  |
| traceratops/core/plotting_functions.py:300 | `plot_matrix_difference` |  |
| traceratops/core/plotting_functions.py:355 | `plot_mixed_matrix` |  |
| traceratops/core/trace_resolver/beam.py:25 | `beam_search` |  |
| traceratops/core/trace_resolver/classifier.py:43 | `classify_trace` | Classify a trace without embedding a fixed biological interpretation. |
| traceratops/core/trace_resolver/classifier.py:58 | `build_multiplicity_matrix` | Return ``(trace_ids, barcode_ids, counts)`` for the full barcode universe. |
| traceratops/core/trace_resolver/classifier.py:77 | `barcode_heterogeneity_pvalue` | Test whether dataset-wide excess multiplicity is uniform by barcode. |
| traceratops/core/trace_resolver/classifier.py:166 | `assess_auto_nuisance_reliability` | Assess benchmark-informed ambiguity without changing fitted results. |
| traceratops/core/trace_resolver/classifier.py:195 | `assess_likelihood_classifier_reliability` | Assess an empirical warning regime without changing classification. |
| traceratops/core/trace_resolver/classifier.py:324 | `LikelihoodMultiplicityClassifier.fit` | Fit dataset-level parameters from a trace-by-barcode count matrix. |
| traceratops/core/trace_resolver/classifier.py:363 | `LikelihoodMultiplicityClassifier.score` | Score and classify one or more count vectors. |
| traceratops/core/trace_resolver/model.py:53 | `EmpiricalDistanceModel.positions` |  |
| traceratops/core/trace_resolver/model.py:68 | `EmpiricalDistanceModel.fit` |  |
| traceratops/core/trace_resolver/model.py:138 | `EmpiricalDistanceModel.for_separation` |  |
| traceratops/core/trace_resolver/resolver.py:34 | `TraceResolver.resolve` |  |
| traceratops/core/trace_resolver/resolver.py:100 | `TraceResolver.resolve_one_polymer_candidates` | Resolve duplicates using a fresh constrained search per candidate. |
| traceratops/core/trace_resolver/results.py:20 | `TraceMultiplicity.from_barcodes` |  |
| traceratops/core/trace_resolver/results.py:48 | `ResolutionResult.n_unassigned` |  |
| traceratops/core/trace_resolver/results.py:52 | `ResolutionResult.fraction_unassigned` |  |
| traceratops/core/trace_resolver/scoring.py:27 | `PolymerScorer.transition_cost` |  |
| traceratops/localization_analyzer.py:20 | `parse_arguments` |  |
| traceratops/localization_analyzer.py:48 | `create_dict_args` |  |
| traceratops/localization_analyzer.py:61 | `get_output_file` |  |
| traceratops/localization_analyzer.py:69 | `run` |  |
| traceratops/localization_analyzer.py:79 | `main` |  |
| traceratops/localization_merge.py:22 | `parse_arguments` |  |
| traceratops/localization_merge.py:35 | `create_dict_args` |  |
| traceratops/localization_merge.py:63 | `appends_traces` |  |
| traceratops/localization_merge.py:82 | `load_localizations` |  |
| traceratops/localization_merge.py:90 | `run` |  |
| traceratops/localization_merge.py:110 | `main` |  |
| traceratops/plot_3way_coloc.py:23 | `parse_arguments` |  |
| traceratops/plot_3way_coloc.py:50 | `plot_threeway_matrix` |  |
| traceratops/plot_3way_coloc.py:137 | `extract_anchor` |  |
| traceratops/plot_3way_coloc.py:142 | `load_list` |  |
| traceratops/plot_3way_coloc.py:154 | `load_label_map` |  |
| traceratops/plot_3way_coloc.py:163 | `main` |  |
| traceratops/plot_4m.py:30 | `parse_arguments` |  |
| traceratops/plot_4m.py:72 | `compute_colocalization` | Computes the frequency of colocalization between the anchor barcode and all other barcodes. |
| traceratops/plot_4m.py:121 | `bootstrap_colocalization` | Performs bootstrapping to estimate mean and SEM of colocalization frequencies for multiple anchors. |
| traceratops/plot_4m.py:175 | `plot_frequencies` | Plots colocalization frequencies for multiple anchors separately. |
| traceratops/plot_4m.py:224 | `main` |  |
| traceratops/plot_bootstrapping.py:22 | `parse_arguments` |  |
| traceratops/plot_bootstrapping.py:94 | `create_dict_args` |  |
| traceratops/plot_bootstrapping.py:174 | `plot_results` |  |
| traceratops/plot_bootstrapping.py:205 | `main` |  |
| traceratops/plot_compare2matrices.py:24 | `parse_arguments` |  |
| traceratops/plot_compare2matrices.py:101 | `create_dict_args` |  |
| traceratops/plot_compare2matrices.py:191 | `gets_ensemble_matrix` |  |
| traceratops/plot_compare2matrices.py:234 | `main` |  |
| traceratops/plot_him_matrix.py:26 | `parse_arguments` |  |
| traceratops/plot_him_matrix.py:141 | `check_required_arg` | Check that no required arguments are supplied. |
| traceratops/plot_him_matrix.py:155 | `create_output_folder` |  |
| traceratops/plot_him_matrix.py:161 | `load_matrix` |  |
| traceratops/plot_him_matrix.py:168 | `load_barcodes` |  |
| traceratops/plot_him_matrix.py:175 | `shuffle_matrix` |  |
| traceratops/plot_him_matrix.py:193 | `new_shuffle_matrix` |  |
| traceratops/plot_him_matrix.py:200 | `merge_matrices` |  |
| traceratops/plot_him_matrix.py:215 | `apply_nan_threshold` |  |
| traceratops/plot_him_matrix.py:224 | `apply_triangular_mask` | Apply triangular mask to matrix. |
| traceratops/plot_him_matrix.py:258 | `main` |  |
| traceratops/plot_matrix_comparison.py:28 | `parse_arguments` |  |
| traceratops/plot_matrix_comparison.py:57 | `create_dict_args` |  |
| traceratops/plot_matrix_comparison.py:101 | `plot_result` |  |
| traceratops/plot_matrix_comparison.py:139 | `calculates_pearson_correlation` |  |
| traceratops/plot_matrix_comparison.py:144 | `remove_zeros` |  |
| traceratops/plot_matrix_comparison.py:149 | `plots_distributions` | makes violin plots containing each dataset. |
| traceratops/plot_matrix_comparison.py:181 | `parses_matrix_to_vector` |  |
| traceratops/plot_matrix_comparison.py:197 | `load_matrix` |  |
| traceratops/plot_matrix_comparison.py:201 | `calculates_ensemble_matrices` |  |
| traceratops/plot_matrix_comparison.py:224 | `common_member` |  |
| traceratops/plot_matrix_comparison.py:229 | `filters_zero_values` | finds the list of common indices in these arrays that contain nonzero values |
| traceratops/plot_matrix_comparison.py:243 | `main_script` |  |
| traceratops/plot_matrix_comparison.py:260 | `main` |  |
| traceratops/pymol_script.py:37 | `color_all_barcodes` |  |
| traceratops/pymol_script.py:56 | `load_pdb_grid` | Load the first PDB files in a folder and display them in a PyMOL grid |
| traceratops/pymol_script.py:125 | `spin_grid` | Create and play a continuously looping rotation movie. |
| traceratops/script_banner.py:33 | `print_script_banner` | Print a standard banner with the script name and a short description. |
| traceratops/trace_3way_coloc.py:138 | `compute_threeway_colocalization` | Computes the frequency of three-way co-localization between an anchor barcode |
| traceratops/trace_3way_coloc.py:164 | `bootstrap_threeway_colocalization` | Performs bootstrapping to estimate mean and SEM of three-way co-localization frequencies. |
| traceratops/trace_3way_coloc.py:225 | `plot_threeway_matrix` | Creates a heatmap of three-way co-localization frequencies using matplotlib. |
| traceratops/trace_3way_coloc.py:399 | `parse_arguments` |  |
| traceratops/trace_3way_coloc.py:449 | `get_trace_files` |  |
| traceratops/trace_3way_coloc.py:461 | `main` |  |
| traceratops/trace_analyzer.py:30 | `parse_arguments` |  |
| traceratops/trace_analyzer.py:61 | `create_dict_args` |  |
| traceratops/trace_analyzer.py:90 | `get_barcode_statistics` | Calculate and visualize barcode statistics from trace data. |
| traceratops/trace_analyzer.py:235 | `plot_neighbor_distances` | Calculate and visualize distances between consecutive neighboring barcodes. |
| traceratops/trace_analyzer.py:438 | `barcode_detection_efficiency` | Analyze and visualize barcode detection efficiency across all traces. |
| traceratops/trace_analyzer.py:520 | `plot_kde_projections` | Plot hexbin projections (XY, XZ, YZ) from trace-centered data. |
| traceratops/trace_analyzer.py:625 | `analyze_trace` | Perform comprehensive analysis on a chromatin trace file. |
| traceratops/trace_analyzer.py:696 | `process_traces` | Process a list of trace files and analyze each individually. |
| traceratops/trace_analyzer.py:754 | `main` |  |
| traceratops/trace_assign_mask.py:20 | `parse_arguments` |  |
| traceratops/trace_assign_mask.py:44 | `create_dict_args` |  |
| traceratops/trace_assign_mask.py:84 | `load_mask` |  |
| traceratops/trace_assign_mask.py:105 | `assign_masks` |  |
| traceratops/trace_assign_mask.py:147 | `plot_mask_assignment` | Plot assigned vs non-assigned points over mask image. |
| traceratops/trace_assign_mask.py:199 | `process_traces` |  |
| traceratops/trace_assign_mask.py:231 | `main` |  |
| traceratops/trace_export_to_fofct.py:29 | `parse_arguments` |  |
| traceratops/trace_export_to_fofct.py:49 | `get_trace_files` |  |
| traceratops/trace_export_to_fofct.py:64 | `load_trace_ecsv_file` |  |
| traceratops/trace_export_to_fofct.py:71 | `load_barcode_bed_file` |  |
| traceratops/trace_export_to_fofct.py:85 | `load_json_file` |  |
| traceratops/trace_export_to_fofct.py:92 | `check_metadata` |  |
| traceratops/trace_export_to_fofct.py:100 | `get_output_file` |  |
| traceratops/trace_export_to_fofct.py:113 | `find_chrom_info` |  |
| traceratops/trace_export_to_fofct.py:129 | `assign_chrom_info` |  |
| traceratops/trace_export_to_fofct.py:146 | `remove_unused_columns` |  |
| traceratops/trace_export_to_fofct.py:152 | `rename_pyhim_columns` |  |
| traceratops/trace_export_to_fofct.py:161 | `write_fof_ct_csv` |  |
| traceratops/trace_export_to_fofct.py:173 | `get_initials` |  |
| traceratops/trace_export_to_fofct.py:183 | `get_authors_from_copyright` | Read author names from installed metadata, without checkout-only files. |
| traceratops/trace_export_to_fofct.py:199 | `get_header_comments` |  |
| traceratops/trace_export_to_fofct.py:222 | `convert_ecsv_to_csv` |  |
| traceratops/trace_export_to_fofct.py:247 | `main` |  |
| traceratops/trace_filter.py:62 | `check_required_arg` | Check that no required arguments are supplied. |
| traceratops/trace_filter.py:74 | `parse_arguments` |  |
| traceratops/trace_filter.py:175 | `args_coord_to_dict` |  |
| traceratops/trace_filter.py:186 | `args_quality_filters_to_dict` |  |
| traceratops/trace_filter.py:207 | `get_files_from_args` | Read input directly without using select.select() |
| traceratops/trace_filter.py:216 | `check_file_number` | checks number of trace files |
| traceratops/trace_filter.py:230 | `filter_duplicate` |  |
| traceratops/trace_filter.py:244 | `filter_barcode_number` |  |
| traceratops/trace_filter.py:252 | `filter_label` |  |
| traceratops/trace_filter.py:264 | `get_filtered_localizations` | Return localizations whose ``Buid`` is still present in the trace table. |
| traceratops/trace_filter.py:277 | `get_file_root` | Return a file path without its final extension. |
| traceratops/trace_filter.py:282 | `runtime` |  |
| traceratops/trace_filter.py:420 | `main` |  |
| traceratops/trace_genomic_coordinates.py:18 | `parse_arguments` |  |
| traceratops/trace_genomic_coordinates.py:39 | `create_dict_args` |  |
| traceratops/trace_genomic_coordinates.py:61 | `load_bed_file` | Loads the BED file into a dictionary mapping barcodes to genomic coordinates. |
| traceratops/trace_genomic_coordinates.py:123 | `impute_genomic_coordinates` | Updates the Chrom, Chrom_Start, Chrom_End, and Barcode # columns in the trace file based on the BED file. |
| traceratops/trace_genomic_coordinates.py:186 | `main` |  |
| traceratops/trace_import_from_fofct.py:26 | `parse_arguments` |  |
| traceratops/trace_import_from_fofct.py:40 | `get_trace_files` |  |
| traceratops/trace_import_from_fofct.py:55 | `read_column_names_from_csv` |  |
| traceratops/trace_import_from_fofct.py:66 | `load_csv_file` |  |
| traceratops/trace_import_from_fofct.py:76 | `load_barcode_bed_file` |  |
| traceratops/trace_import_from_fofct.py:84 | `validate_columns` |  |
| traceratops/trace_import_from_fofct.py:94 | `find_barcode_id` |  |
| traceratops/trace_import_from_fofct.py:108 | `add_missing_columns` |  |
| traceratops/trace_import_from_fofct.py:131 | `rename_columns_for_ecsv` |  |
| traceratops/trace_import_from_fofct.py:140 | `convert_csv_to_ecsv` |  |
| traceratops/trace_import_from_fofct.py:148 | `main` |  |
| traceratops/trace_merge.py:29 | `silent_load` |  |
| traceratops/trace_merge.py:37 | `parse_arguments` |  |
| traceratops/trace_merge.py:61 | `check_required_arg` | Check that no required arguments are supplied. |
| traceratops/trace_merge.py:73 | `get_files_from_args` | Read input directly without using select.select() |
| traceratops/trace_merge.py:82 | `appends_traces` |  |
| traceratops/trace_merge.py:98 | `load_traces` |  |
| traceratops/trace_merge.py:109 | `create_out_folder` |  |
| traceratops/trace_merge.py:115 | `main` |  |
| traceratops/trace_pearsons.py:35 | `silent_load_trace` |  |
| traceratops/trace_pearsons.py:43 | `parse_arguments` |  |
| traceratops/trace_pearsons.py:76 | `get_trace_files` |  |
| traceratops/trace_pearsons.py:91 | `find_unique_substrings` | Find unique substrings in each filename that distinguish them from all others. |
| traceratops/trace_pearsons.py:202 | `accumulate_distances` | Calculate pairwise distances between barcodes across all traces. |
| traceratops/trace_pearsons.py:235 | `compare_distance_maps` | Calculate Pearson correlation between distance maps from different trace files. |
| traceratops/trace_pearsons.py:280 | `plot_correlation_matrix` | Plot a correlation matrix between files with unique identifiers as labels. |
| traceratops/trace_pearsons.py:351 | `main` |  |
| traceratops/trace_physical_vs_genomic_distance.py:38 | `extract_traces_numpy` | Extract trace coordinates with shape ``(num_traces, num_barcodes, 3)``. |
| traceratops/trace_physical_vs_genomic_distance.py:64 | `compute_genomic_dist_map` | Compute barcode-pair genomic separations in kilobase pairs. |
| traceratops/trace_physical_vs_genomic_distance.py:71 | `compute_inter_loci_genomic_dist` | Compute genomic separations between consecutive barcode loci in kilobase pairs. |
| traceratops/trace_physical_vs_genomic_distance.py:84 | `compute_all_pwd` | Compute 3D and per-axis pairwise physical distances for every trace. |
| traceratops/trace_physical_vs_genomic_distance.py:109 | `compute_physical_vs_genomic_distance` | Bin physical distances by genomic distance and return median distance summaries. |
| traceratops/trace_physical_vs_genomic_distance.py:164 | `fit_power_law` | Fit ``physical_distance = coefficient * genomic_distance ** exponent``. |
| traceratops/trace_physical_vs_genomic_distance.py:217 | `plot_log_distance_graph` | Plot log10 median physical distance versus log10 genomic distance. |
| traceratops/trace_physical_vs_genomic_distance.py:290 | `get_plot_output_file` | Return the plot output path with the requested file extension. |
| traceratops/trace_physical_vs_genomic_distance.py:311 | `calculate_physical_vs_genomic_distance` | Load a trace table, calculate R(s) summaries, and write output files. |
| traceratops/trace_physical_vs_genomic_distance.py:366 | `parse_arguments` |  |
| traceratops/trace_physical_vs_genomic_distance.py:424 | `main` |  |
| traceratops/trace_plot.py:24 | `parse_arguments` |  |
| traceratops/trace_plot.py:54 | `create_dict_args` |  |
| traceratops/trace_plot.py:116 | `runtime` |  |
| traceratops/trace_plot.py:195 | `main` |  |
| traceratops/trace_sorter.py:19 | `parse_arguments` | Build the command-line parser. |
| traceratops/trace_sorter.py:42 | `get_trace_files` | Return filenames supplied directly or over standard input. |
| traceratops/trace_sorter.py:51 | `sorted_output_filename` | Insert ``_sorted`` before an input filename's extension. |
| traceratops/trace_sorter.py:57 | `sort_trace_table` | Sort a table by a trace-level property, modifying it in place. |
| traceratops/trace_sorter.py:87 | `main` | Load, sort, and save each requested chromatin trace table. |
| traceratops/trace_split_labels.py:24 | `parse_arguments` |  |
| traceratops/trace_split_labels.py:49 | `get_files` |  |
| traceratops/trace_split_labels.py:57 | `split_trace_by_label` |  |
| traceratops/trace_split_labels.py:80 | `process` |  |
| traceratops/trace_split_labels.py:111 | `main` |  |
| traceratops/trace_splitter.py:56 | `parse_arguments` | Return the command-line argument parser. |
| traceratops/trace_splitter.py:201 | `create_dict_args` |  |
| traceratops/trace_splitter.py:213 | `generate_unique_id` | Generate a unique trace identifier. |
| traceratops/trace_splitter.py:218 | `compute_radius_of_gyration` | Compute the radius of gyration for a trace. |
| traceratops/trace_splitter.py:257 | `estimate_hdbscan_epsilon` | Estimate epsilon as the median per-trace Otsu pair-distance threshold. |
| traceratops/trace_splitter.py:310 | `split_large_traces` | Split selected traces in place with K-means or HDBSCAN. |
| traceratops/trace_splitter.py:524 | `resolve_traces` | Resolve barcode duplicates and merged pairs, returning diagnostics. |
| traceratops/trace_splitter.py:759 | `main` |  |
| traceratops/trace_stats.py:18 | `parse_arguments` |  |
| traceratops/trace_stats.py:28 | `get_trace_files` |  |
| traceratops/trace_stats.py:43 | `compute_trace_statistics` |  |
| traceratops/trace_stats.py:62 | `main` |  |
| traceratops/trace_to_matrix.py:18 | `parse_arguments` |  |
| traceratops/trace_to_matrix.py:45 | `create_dict_args` |  |
| traceratops/trace_to_matrix.py:93 | `runtime` |  |
| traceratops/trace_to_matrix.py:137 | `main` |  |

## Other executable files

All top-level Python files with main functions are represented by registered entry points. `pymol_script.py` is a documented PyMOL helper script rather than a standalone CLI; its public PyMOL functions are cataloged separately. `_version.py` and `script_banner.py` are metadata/presentation infrastructure. Installation, release, and documentation shell scripts are project maintenance tools, not scientific traceratops commands.

## Proposed schema (used by the catalog)

Versioned root: schema_version, repository (URL, branch, commit), input_schemas, cli_tools, python_apis. CLI records: id, command, entry_point, alias_of, description, scientific_purpose, operations, concepts, select_for, inputs, runtime_requirements, parser (arguments and mutually-exclusive groups), runtime_defaults, outputs, modifies_input_files, examples, limitations, references, verification. Python records additionally contain module, callable, signature, documented_return_type, documentation, and caveats. Tasks contain synonyms, distinctions, tool IDs, and ordered workflow steps. Examples use argv arrays plus stdin, never shell interpolation. Parser defaults are separate from runtime transformations; null is an actual parser default, not an unknown value.

## Final inclusion ledger

`python_surface.yaml` classifies all 310 non-underscore function/method names found by AST traversal: 76 analytical callables are cataloged; every other candidate has an explicit exclusion reason. This includes undocumented public-name functions, which are documentation gaps rather than invented public interfaces. All 29 installed command entry points were executed.
