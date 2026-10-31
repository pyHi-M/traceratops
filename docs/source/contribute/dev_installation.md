# Installation [DEV]

Requires Python **3.11 or newer**. The test matrix covers Python 3.11, 3.12 and 3.13 on Linux, macOS and Windows. `uv` is the recommended installation tool.

The recommended development setup uses [`uv`](https://docs.astral.sh/uv/) for both the virtual environment and package installation. The repository includes `install_traceratops_uv.bash`, which installs `uv` if needed, clones traceratops, creates a fresh `.venv`, installs traceratops with development dependencies, and runs a small import/CLI check.

## Automated uv installation

From any terminal, run:

```bash
curl -O https://raw.githubusercontent.com/pyHi-M/traceratops/main/install_traceratops_uv.bash
bash install_traceratops_uv.bash
```

By default, the script:

- clones traceratops into `$HOME/Repositories/traceratops`;
- creates a fresh uv environment at `$HOME/Repositories/traceratops/.venv`;
- uses Python 3.11;
- installs traceratops in editable mode with the `dev` optional dependencies;
- checks that core scientific dependencies and the `trace_filter` command are available.

Activate the environment after installation with:

```bash
source $HOME/Repositories/traceratops/.venv/bin/activate
```

## Manual uv installation

If you prefer to run the commands yourself, install `uv` first if it is missing:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
```

Then clone and install traceratops:

```bash
mkdir -p $HOME/Repositories
cd $HOME/Repositories
git clone https://github.com/pyHi-M/traceratops.git
cd traceratops
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install -e ".[dev]"
```


## Check the installation

```bash
python -c "import traceratops; print('traceratops OK')"
trace_filter --help
```

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
