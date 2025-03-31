## Type of Change

- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Code refactoring
- [ ] Other (please specify):

## Version Update

[Depending on the type of change](https://semver.org/), please update the version number in the `pyproject.toml` file.

- **Bug fix** (*patch*): `X.Y.Z` → `X.Y.Z+1`
- **New feature** (*minor*): `X.Y.Z` → `X.Y+1.0`
- **Breaking change** (*major*): `X.Y.Z` → `X+1.0.0`

## Tests

Tests can be performed in two ways:

- A `test/test_*.py` file for automation with [`pytest`](https://docs.pytest.org/en/stable/) [(cheatsheet)](https://cheatography.com/hvid2301/cheat-sheets/pytest-usage/)
- Updating the private repository for manual testing

Tests are required to:

- Validate the new feature
- Confirm that the bug no longer occurs

## Additional Comments

<!-- Add any extra information that might be helpful to the reviewers. -->

## Checklist

- [ ] Version number has been updated
- [ ] Tests are available
- [ ] [Documentation has been updated (tutorial)]()
- [ ] **Squash commits**: I understand that the detailed commit history of this branch will be squashed into a single commit
- [ ] [Code has been reviewed]()
- [ ] GitHub Actions checks pass successfully
