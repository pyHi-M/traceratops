#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

REPO_DIR="${HOME}/Repositories"
PYTHON_VERSION="3.11"
REPO_URL="https://github.com/pyHi-M/traceratops.git"

# -----------------------------
# Install uv if missing
# -----------------------------

if ! command -v uv >/dev/null 2>&1; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="${HOME}/.local/bin:${PATH}"
else
    echo "uv already installed."
fi

# -----------------------------
# Clone repository
# -----------------------------

mkdir -p "${REPO_DIR}"
cd "${REPO_DIR}"

if [ ! -d "traceratops" ]; then
    git clone "${REPO_URL}"
else
    echo "traceratops repository already exists."
fi

# -----------------------------
# Create fresh uv environment
# -----------------------------

cd "${REPO_DIR}/traceratops"

if [ -d ".venv" ]; then
    echo "Removing existing .venv..."
    rm -rf .venv
fi

uv venv .venv --python "${PYTHON_VERSION}"
source .venv/bin/activate

# -----------------------------
# Install traceratops with development dependencies
# -----------------------------

uv pip install -e ".[dev]"

# -----------------------------
# Test installation
# -----------------------------

python -c "import numpy, scipy, astropy; print('numpy', numpy.__version__); print('scipy', scipy.__version__); print('astropy', astropy.__version__)"
python -c "import traceratops; print('traceratops OK')"
trace_filter --help >/dev/null

cat <<EOM

Installation complete.
Activate with:
source ${REPO_DIR}/traceratops/.venv/bin/activate

Test traceratops with:
trace_filter --help
EOM
