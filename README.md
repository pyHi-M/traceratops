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

Stable releases are distributed through **PyPI**. Python **3.11 or newer** is
required; you do not need to clone this repository for a stable installation.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install traceratops
python -c "from importlib.metadata import version; print(version('traceratops'))"
trace_filter --help
```

On Windows PowerShell, create the environment with `py -3.11 -m venv .venv` and
activate with `.\.venv\Scripts\Activate.ps1`; Command Prompt uses
`.venv\Scripts\activate.bat`. See the [installation guide](docs/source/quickstart/installation.md)
for commands that do not require activation.

Stable users who prefer [uv](https://docs.astral.sh/uv/) can use:

```bash
uv venv --python 3.11
uv pip install traceratops
uv run --no-project trace_filter --help
```

The initial PyPI publication of 2.0.0 is pending separate maintainer action.
The [installation guide](docs/source/quickstart/installation.md) includes an
existing-tag source alternative until PyPI publication succeeds.

### Development version

Use uv with an editable checkout of **dev**:

```bash
git clone https://github.com/pyHi-M/traceratops.git
cd traceratops
git switch dev
uv venv --python 3.11
uv pip install -e ".[dev]"
uv run --no-project pytest -q
```

See the [development guide](docs/source/contribute/dev_installation.md) for the
optional Bash installer and updating a checkout. `pyproject.toml` specifies
library dependencies; no universal `uv.lock` is committed. Runtime/development
security audits run in GitHub Actions.

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
