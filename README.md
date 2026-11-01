# traceratops

Requires Python **3.11 or newer**. The test matrix covers Python 3.11, 3.12 and
3.13 on Linux, macOS and Windows. Stable users install from **PyPI**, preferably
with **pip**; developers use **uv** with an editable **dev** checkout.

*A toolbox for trace analysis, visualization, and quality assessment for chromatin-trace data.*

Raw chromatin tracing experiments consist of hundreds to thousands of images that need to be analyzed to extract and decode spot localizations for each genomic loci, for example using [pyHiM](https://github.com/pyHi-M/pyHiM). This analysis produces trace tables containing the genomic identities and 3D localizations of each localized spot for each chromatin trace (trace table).

`traceratops` is a software package providing multiple tools for the post-processing of trace tables. These tools help users clean, merge, inspect, transform, and visualize chromatin tracing data in formats compatible with the [4DN FISH Omics Format](https://fish-omics-format.readthedocs.io/en/latest/).

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

For stable use, install from **PyPI**, preferably with pip in a virtual environment.
Python **3.11 or newer** is required; cloning the repository is unnecessary.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install traceratops
```

uv users can instead run `uv venv --python 3.11` followed by
`uv pip install traceratops`. See the [complete installation guide](docs/source/quickstart/installation.md)
for Windows instructions, verification and updates.

Developers use **uv** with an editable checkout of **dev** and
`uv pip install -e ".[dev]"`. The [development guide](docs/source/contribute/dev_installation.md)
covers checkout setup, tests, updates and the optional installer.

## Quick start

A practical way to start is to run the first tutorial notebook:

```bash
jupyter notebook docs/source/tutorials/tutorial_merge_multi_roi.ipynb
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
