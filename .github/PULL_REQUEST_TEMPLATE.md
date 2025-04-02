## Type of Change

- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Code refactoring
- [ ] Other (please specify):

## Tests

Tests can be performed in two ways:

- A `test/test_*.py` file for automation with [`pytest`](https://docs.pytest.org/en/stable/) [(cheatsheet)](https://cheatography.com/hvid2301/cheat-sheets/pytest-usage/)
- Updating the private repository for manual testing

Tests are required to:

- Validate the new feature
- Confirm that the bug no longer occurs

## [DEV] Checklist

**Only for Pull Request on `dev` branch.**

- [ ] Tests are available
- [ ] [Documentation has been updated (tutorial)](https://traceratops.readthedocs.io/en/main/contribute/how_to_document.html)
- [ ] [Argument status are updated](https://traceratops.readthedocs.io/en/main/contribute/how_to_document.html#inside-your-python-script)
- [ ] **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit
- [ ] [Code has been reviewed](https://traceratops.readthedocs.io/en/main/contribute/code_review.html)
- [ ] GitHub Actions checks pass successfully

## [MAIN] Checklist

**Only for Pull Request on `main` branch.**

- [ ] Version number has been updated inside `pyproject.toml` [(Depending on the type of change)](https://semver.org/)
- [ ] The status of at least one feature are updated to `stable` (for each argument involved) OR it's a **hotfix** branch
- [ ] Manual tests are checked
- [ ] **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit
- [ ] GitHub Actions checks pass successfully
