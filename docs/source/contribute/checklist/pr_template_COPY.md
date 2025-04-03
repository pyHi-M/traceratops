# Pull Request template (COPY)

## Type of Change

<label><input type="checkbox"> Bug fix</label>
<label><input type="checkbox"> New feature</label>
<label><input type="checkbox"> Documentation update</label>
<label><input type="checkbox"> Test update</label>
<label><input type="checkbox"> Code refactoring</label>
<label><input type="checkbox"> Other (please specify):</label>

## Tests

Tests can be performed in two ways:

- [Updating the private repository for manual testing](https://traceratops.readthedocs.io/en/latest/contribute/checklist/manual_test_checklist.html)
- A `test/test_*.py` file for automation with [`pytest`](https://docs.pytest.org/en/stable/) [(cheatsheet)](https://cheatography.com/hvid2301/cheat-sheets/pytest-usage/)

Tests are required to:

- Validate the new feature
- Confirm that the bug no longer occurs

## [DEV] Checklist

**Only for Pull Request on `dev` branch.**

<label><input type="checkbox"> [Documentation has been updated (tutorial)](https://traceratops.readthedocs.io/en/latest/contribute/how_to_document.html)</label>
<label><input type="checkbox"> [Feature reliability status is updated](https://traceratops.readthedocs.io/en/latest/contribute/how_to_document.html#reliability-status)</label>
<label><input type="checkbox"> **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit</label>
<label><input type="checkbox"> [Code has been reviewed](https://traceratops.readthedocs.io/en/latest/contribute/checklist/code_review.html)</label>
<label><input type="checkbox"> Checks of GitHub Actions pass successfully</label>

## [MAIN] Checklist

**Only for Pull Request on `main` branch.**

<label><input type="checkbox"> Version number has been updated inside `pyproject.toml` [(Depending on the type of change)](https://semver.org/)</label>
<label><input type="checkbox"> The status of at least one feature are updated to `stable` OR it's a **hotfix** branch</label>
<label><input type="checkbox"> Manual tests are checked</label>
<label><input type="checkbox"> **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit</label>
<label><input type="checkbox"> `codespell` pass successfully</label>
<label><input type="checkbox"> `mypy` pass successfully</label>
<label><input type="checkbox"> Checks of GitHub Actions pass successfully</label>
