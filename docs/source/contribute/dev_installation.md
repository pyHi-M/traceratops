# Installation [DEV]

Use **uv** for development on Python **3.11 or newer**. Stable users should follow
the [PyPI installation guide](../quickstart/installation.md) instead.

## Manual development installation

[Install uv](https://docs.astral.sh/uv/getting-started/installation/), then select
`dev` explicitly: cloning alone selects the repository's default `main` branch.

```bash
git clone https://github.com/pyHi-M/traceratops.git
cd traceratops
git switch dev
uv venv --python 3.11
uv pip install -e ".[dev]"
```

`-e` installs traceratops in editable mode: Python source edits in this checkout
are immediately visible in the environment. `[dev]` installs test, pre-commit,
typing and documentation tools, which stable users do not need.

## Verify and run tests

```bash
uv run --no-project python -c "from importlib.metadata import version; print(version('traceratops'))"
uv run --no-project trace_filter --help
uv pip check
uv run --no-project pytest -vv -s
```

`uv run --no-project` targets `.venv` without activation and does not create a
lockfile. It works on Linux, macOS and Windows. To activate manually, use
`source .venv/bin/activate` in Bash, `.\.venv\Scripts\Activate.ps1` in Windows
PowerShell, or `.venv\Scripts\activate.bat` in Command Prompt.

## Update an existing checkout

Commit or stash local changes before switching branches or pulling. For a clean
tracking checkout of `dev`:

```bash
git fetch origin
git switch dev
git pull --ff-only origin dev
uv pip install -e ".[dev]"
uv run --no-project pytest -q
```

Reinstall after dependency, entry-point or version metadata changes. Work on a
separate feature branch when contributing; do not pull dev over uncommitted work.
Recreate `.venv` when changing Python versions or checking a fresh resolution.

## Optional Bash development installer

On Linux/macOS (or a suitable Bash environment):

```bash
curl -O https://raw.githubusercontent.com/pyHi-M/traceratops/dev/install_traceratops_uv.bash
bash install_traceratops_uv.bash
```

The installer defaults to branch `dev`, Python 3.11 and
`$HOME/Repositories/traceratops`. It recreates `.venv`, installs `.[dev]` editable,
and checks imports and `trace_filter`. It selects the requested branch but does
not pull an existing checkout; update that checkout explicitly as above.

Set `TRACERATOPS_BRANCH` to choose another existing branch, for example
`TRACERATOPS_BRANCH=main bash install_traceratops_uv.bash`. Use `TRACERATOPS_REPO_DIR` to override the parent repository directory.
For an existing checkout, switching requires a clean working tree; local changes are never
stashed or discarded automatically. The script is for source/development setup,
not the recommended stable PyPI installation.

## [Recommended] Setup pre-commit in local

- check if it's well installed

  `pre-commit --version`

- install command of the file ".pre-commit-config.yaml" inside ".git/hooks/pre-commit"

  `pre-commit install`

- test pre-commit without any commit

  `pre-commit run --all-files`

```{note}
- First time, it can take few minutes because it need to install all dependencies.
- Then, you can re-run this command when you want to simulate action that are execute, by hidden way, when you try to make a commit.
- For info, if one step of the pre-commit fail, your commit fail and you need to fix pre-commit error to be allow to commit.
```

### Tips

- [Optional] fix strange issue or warning

  `pre-commit autoupdate --repo https://github.com/pre-commit/pre-commit-hooks`

- [Optional] Update pre-commit file
  - `pre-commit clean`
  - `pre-commit autoupdate`
  - `pre-commit install`

### Codespell

- [Very optional] Activate locally codespell inside .pre-commit-config.yaml
  - Uncomment codespell section
  - Ignore locally modification on the `.pre-commit-config.yaml` file:
    ```git update-index --assume-unchanged .pre-commit-config.yaml```
  - To undo: ```git update-index --no-assume-unchanged .pre-commit-config.yaml```
  - to see 'assume unchanged' files ("-v" option ==> use lowercase letters)
    ```git ls-files -v```

## Strange error

- If a script name seems to refer to another package:

  - Identify the error package name.
  - If it's not in the dependencies of your project, it's a "ghost".
  - To exorcise it:
    - Identify exe path : `which pylint` (for example `/home/user/.local/bin/pylint`)
    - Remove it : `rm ~/.local/bin/pylint`
  - Re-install your package to be sure `uv pip install -e ".[dev]"`

## Dependency maintenance and security

`pyproject.toml` is the authoritative dependency specification. As a library,
traceratops does not commit `uv.lock`; `uv pip install -e ".[dev]"` resolves the
specified requirements without a lockfile. Use a fresh environment when checking
compatibility. To run tests without shell-specific activation (including Windows):

```bash
uv run --no-project pytest -vv -s
uv run --no-project trace_filter --help
```

The dependency security audit runs weekly, manually, and on pull requests to
`dev` or `main` that change dependency or audit configuration. Separate runtime
and development jobs audit the installed dependency snapshots with `pip-audit`
using the PyPI vulnerability service. Findings fail the job; logs and JSON
artifacts contain package versions, vulnerability identifiers and available fixes.
The audit does not change dependencies or merge fixes. The development report
includes runtime dependencies as well as development tools.

Maintainers should enable GitHub Actions notifications in their personal
**Settings → Notifications → System → Actions** (email and/or web, optionally
failed workflows only), and watch the repository's Actions activity. Enable the
repository dependency graph and Dependabot alerts under **Settings → Advanced
Security** (or **Code security**), and enable personal Dependabot alert
notifications. Alert access depends on repository permissions. No external
notification service or credentials are needed.

Scheduled workflows run from the default branch: the audit workflow must reach
the default branch before weekly runs start. Until then, use the pull request
checks; manual dispatch also requires the workflow on the default branch.
GitHub may disable scheduled workflows in inactive public repositories after
60 days; maintainers can re-enable them from the Actions tab.
