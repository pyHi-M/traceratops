traceratops documentation
=========================

**Date**: |today| **Version**: |release|

traceratops is an open-source Python toolbox for reproducible post-processing, quality control, analysis, and visualization of reconstructed chromatin tracing datasets. It provides a modular trace-level workflow for harmonizing chromatin trace tables from pyHiM and community-standard formats such as the `4DN FISH Omics Format <https://fish-omics-format.readthedocs.io/en/latest/>`_, assessing trace completeness and barcode detection quality, filtering and annotating curated datasets, generating single-cell and ensemble distance or proximity matrices, analyzing multi-locus spatial interactions, and producing comparative, publication-ready visualizations that support standardized reuse of imaging-based spatial genomics data.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   quickstart/installation
   quickstart/changelog_COPY

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   Getting started with notebooks<tutorials/start_with_ipynb>
   Using the --pipe argument<tutorials/using_pipe>
   Merge Multi-ROI<tutorials/tutorial_01_merge_multi_roi.ipynb>
   Quality Control<tutorials/tutorial_02_quality_control.ipynb>
   Filter Thresholds<tutorials/tutorial_03_filter_thresholds.ipynb>
   Filter Duplicates<tutorials/tutorial_04_filter_duplicates.ipynb>
   Split Traces<tutorials/tutorial_05_split_traces.ipynb>
   Classify Traces<tutorials/tutorial_06_assign_masks_split_labels.ipynb>
   Matrix Visualization<tutorials/tutorial_07_matrix_visualization.ipynb>
   Multiway Co-localization<tutorials/tutorial_08_multiway_coloc.ipynb>
   Compare Datasets<tutorials/tutorial_09_compare_datasets.ipynb>
   Physical vs genomic Distance<tutorials/tutorial_10_physical_vs_genomic_distances.ipynb>

.. toctree::
   :maxdepth: 1
   :caption: Scripts

   Data formats<scripts/data_formats>
   scripts/localization
   scripts/trace
   scripts/matrix
   scripts/plot


.. toctree::
   :maxdepth: 1
   :caption: Contribute

   contribute/CONTRIBUTING
   contribute/dev_installation
   contribute/how_to_document
   contribute/pr_checklists
   contribute/release_guide
