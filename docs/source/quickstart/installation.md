# Installation

Traceratops requires **Python 3.11 or newer**. Python 3.11, 3.12 and 3.13 are
tested on Linux, macOS and Windows.

## Stable releases: install from PyPI

For stable use, install the published package in a virtual environment. You do
**not** need to clone the GitHub repository or install development tools.

On Linux or macOS:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install traceratops
```

Use `python3` instead of `python` if that is your Python >=3.11 executable.
On Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install traceratops
```

For Windows Command Prompt, activate with `.venv\Scripts\activate.bat`.
If PowerShell activation is restricted, use `.\.venv\Scripts\python.exe -m pip
install traceratops` directly; changing system execution policy is unnecessary.

## Stable releases with uv

After [installing uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uv venv --python 3.11
uv pip install traceratops
uv run --no-project python -c "from importlib.metadata import version; print(version('traceratops'))"
uv run --no-project trace_filter --help
```

These commands target `.venv` without shell-specific activation.

## Check and update a stable installation

In the activated environment:

```bash
python -c "from importlib.metadata import version; print(version('traceratops'))"
trace_filter --help
python -m pip check
```

Without activation, use `.venv/bin/python` / `.venv/bin/trace_filter` on Linux or
macOS, or `.venv\Scripts\python.exe` / `.venv\Scripts\trace_filter.exe` on Windows.
To update, run `python -m pip install --upgrade traceratops` (or
`uv pip install --upgrade traceratops`). To select a particular published release,
use `traceratops==2.0.0` instead of the unversioned name.

## Development version

For the latest development version, editable installation, tests and contribution
tools, follow the [development installation guide](../contribute/dev_installation.md).
It explicitly selects the **dev** branch and uses `uv pip install -e ".[dev]"`.
