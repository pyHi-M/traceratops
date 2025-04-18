## Type of Change

- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Test update
- [ ] Code refactoring
- [ ] Other (please specify):

## Tests

Tests can be performed in two ways:

- [Updating the private repository for manual testing](https://traceratops.readthedocs.io/en/latest/contribute/checklist/manual_test_checklist.html)
- A `test/test_*.py` file for automation with [`pytest`](https://docs.pytest.org/en/stable/) [(cheatsheet)](https://cheatography.com/hvid2301/cheat-sheets/pytest-usage/)

Tests are required to:

- Validate the new feature
- Confirm that the bug no longer occurs

## [DEV] Checklist

**Only for Pull Request on `dev` branch.**

- [ ] [Documentation has been updated (tutorial)](https://traceratops.readthedocs.io/en/latest/contribute/how_to_document.html)
- [ ] [Feature reliability status is updated](https://traceratops.readthedocs.io/en/latest/contribute/how_to_document.html#reliability-status)
- [ ] **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit
- [ ] [Code has been reviewed](https://traceratops.readthedocs.io/en/latest/contribute/checklist/code_review.html)
- [ ] Checks of GitHub Actions pass successfully

## [MAIN] Checklist

**Only for Pull Request on `main` branch.**

*For more details, follow the [release guide](https://traceratops.readthedocs.io/en/latest/contribute/release_guide.html).*

- [ ] Version number has been updated inside `pyproject.toml` [(Depending on the type of change)](https://semver.org/)
- [ ] The status of at least one feature are updated to `stable`
- [ ] Manual tests are checked
- [ ] `codespell` does not identify any major problems
- [ ] `mypy` does not identify any major problems
- [ ] Checks of GitHub Actions pass successfully
