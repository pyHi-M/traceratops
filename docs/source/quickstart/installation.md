# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|Yes|Yes|

## Install uv (virtual environment manager)

*We use `uv` environments to avoid version problems with other application dependencies and to keep installation commands reproducible.*

If `uv` is not installed yet, install it with:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
```

## Automated installation

The fastest way to install traceratops from source is to run the bundled uv installer:

```bash
curl -O https://raw.githubusercontent.com/pyHi-M/traceratops/main/install_traceratops_uv.bash
bash install_traceratops_uv.bash
source $HOME/Repositories/traceratops/.venv/bin/activate
```

## Manual installation

```bash
mkdir -p $HOME/Repositories
cd $HOME/Repositories
git clone https://github.com/pyHi-M/traceratops.git
cd traceratops
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install -e ".[dev]"
```
