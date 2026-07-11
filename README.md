# traceratops

*A toolbox for trace analysis, visualization, and quality assessment for chromatin-trace data.*

traceratops provides the post-processing layer for [pyHiM](https://github.com/pyHi-M/pyHiM) multiplexed DNA-FISH workflows. It focuses on the trace and localization tables produced after image processing, and helps users clean, merge, inspect, transform, and visualize those data in formats compatible with the [4DN FISH Omics Format](https://fish-omics-format.readthedocs.io/en/latest/).

## Documentation

The full documentation is available on Read the Docs:

- [traceratops documentation](https://traceratops.readthedocs.io/en/latest/)
- [Installation guide](https://traceratops.readthedocs.io/en/latest/quickstart/installation.html)
- [Tutorial notebooks](https://traceratops.readthedocs.io/en/latest/tutorials/start_with_ipynb.html)

## What traceratops does

traceratops is organized as a set of command-line tools and Python modules for working with pyHiM post-processing outputs:

- **Trace processing**: merge traces from multiple regions of interest, filter traces and localizations, split traces by labels or criteria, assign masks, add genomic coordinates, and compute descriptive statistics.
- **Format conversion**: import and export trace tables using the 4DN FISH Omics Format for Chromatin Tracing (FOF-CT).
- **Matrix generation and analysis**: convert traces into pairwise-distance or proximity matrices, compare matrices, and support downstream chromatin-structure analyses.
- **Visualization**: plot Hi-M matrices, matrix comparisons, bootstrapping results, multi-way co-localization, and other trace-derived summaries.
- **Localization quality control**: collect, merge, and analyze localization tables before or alongside trace-level analysis.

The package is intended for users who already have pyHiM-style localization or trace tables and want a documented, scriptable post-processing workflow.

## Installation

traceratops supports Linux, macOS, and Windows. We recommend installing it in a dedicated `uv` environment. For a development installation, use the bundled installer:

```bash
curl -O https://raw.githubusercontent.com/pyHi-M/traceratops/main/install_traceratops_uv.bash
bash install_traceratops_uv.bash
source $HOME/Repositories/traceratops/.venv/bin/activate
```

To install manually with `uv`, clone the repository and install the development environment:

```bash
git clone https://github.com/pyHi-M/traceratops.git
cd traceratops
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install -e ".[dev]"
```

After installation, traceratops command-line tools such as `trace_filter`, `trace_merge`, `trace_to_matrix`, `trace_stats`, `plot_him_matrix`, and `localization_analyzer` should be available in the active environment.

## Quick start

A practical way to start is to run the first tutorial notebook:

```bash
jupyter notebook docs/source/tutorials/tutorial_01_merge_multi_roi.ipynb
```

Additional tutorials cover quality control, threshold filtering, duplicate removal, trace splitting, mask assignment, matrix visualization, multi-way co-localization, dataset comparison, and use of the `--pipe` argument.

## Major dependencies

traceratops is written in Python and relies on the scientific Python ecosystem. Major runtime dependencies include:

- [NumPy](https://numpy.org/) for numerical arrays and matrix operations.
- [pandas](https://pandas.pydata.org/) for tabular trace and localization data.
- [Astropy](https://www.astropy.org/) for table I/O, including ECSV workflows.
- [SciPy](https://scipy.org/) and [scikit-learn](https://scikit-learn.org/) for scientific computing and analysis utilities.
- [Matplotlib](https://matplotlib.org/) and [Seaborn](https://seaborn.pydata.org/) for plotting and visualization.
- [tqdm](https://tqdm.github.io/) for progress reporting.

Development and documentation extras include [pytest](https://docs.pytest.org/), [pre-commit](https://pre-commit.com/), [mypy](https://mypy-lang.org/), [Sphinx](https://www.sphinx-doc.org/), [sphinx-rtd-theme](https://sphinx-rtd-theme.readthedocs.io/), [sphinx-argparse](https://sphinx-argparse.readthedocs.io/), [MyST Parser](https://myst-parser.readthedocs.io/), and [sphinx-panels](https://sphinx-panels.readthedocs.io/).

## License

traceratops is distributed under the [GNU General Public License v3.0](LICENSE). The GPLv3 is a copyleft free-software license: you may use, study, share, and modify the software, but redistributed copies or derivative works must preserve the same license terms and provide the corresponding source code as required by the license.

The third-party packages used by traceratops are distributed under their own licenses. Please consult each dependency's project page or package metadata for the exact license terms that apply to those dependencies.

## Authors

traceratops is authored by:

- Marcelo Nollmann — `marcelo.nollmann@cbs.cnrs.fr`
- Xavier Devos — `xavier.devos@cbs.cnrs.fr`

The project is part of the pyHi-M software ecosystem developed around pyHiM analysis and post-processing workflows.
