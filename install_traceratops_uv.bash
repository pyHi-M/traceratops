#!/usr/bin/env bash
set -euo pipefail

RESET_VENV=false
if [ "$#" -gt 1 ]; then
    echo "Usage: bash install_traceratops_uv.bash [--reset-venv]" >&2
    exit 2
fi
case "${1:-}" in
    "") ;;
    --reset-venv) RESET_VENV=true ;;
    --help)
        echo "Usage: bash install_traceratops_uv.bash [--reset-venv]"
        echo "Existing environments are reused; --reset-venv explicitly deletes and recreates .venv."
        exit 0
        ;;
    *)
        echo "Unknown option: $1" >&2
        exit 2
        ;;
esac

# -----------------------------
# User settings
# -----------------------------

REPO_DIR="${TRACERATOPS_REPO_DIR:-${HOME}/Repositories}"
PYTHON_VERSION="3.11"
REPO_URL="https://github.com/pyHi-M/traceratops.git"
REPO_BRANCH="${TRACERATOPS_BRANCH:-dev}"
git check-ref-format --branch "${REPO_BRANCH}" >/dev/null

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
    git clone --branch "${REPO_BRANCH}" "${REPO_URL}"
else
    echo "traceratops repository already exists."
fi

# -----------------------------
# Create or reuse the uv environment
# -----------------------------

cd "${REPO_DIR}/traceratops"

# Do not switch a user's dirty checkout or silently update its history.
if [ "$(git branch --show-current)" != "${REPO_BRANCH}" ]; then
    if [ -n "$(git status --porcelain)" ]; then
        echo "Commit or stash local changes before switching to ${REPO_BRANCH}." >&2
        exit 1
    fi
    git switch -- "${REPO_BRANCH}"
fi

if [ "${RESET_VENV}" = true ]; then
    echo "Explicit --reset-venv requested: removing .venv..."
    rm -rf -- .venv
fi

if [ -e ".venv" ] || [ -L ".venv" ]; then
    if [ ! -f ".venv/bin/activate" ]; then
        echo "Existing .venv is not an activatable environment; left untouched." >&2
        echo "Inspect it before explicitly choosing --reset-venv." >&2
        exit 1
    fi
    echo "Reusing existing .venv."
else
    uv venv .venv --python "${PYTHON_VERSION}"
fi
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
