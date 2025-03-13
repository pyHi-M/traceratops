# Contributing guide

## Branches

- `main` is the official branch for users, must be **always stable**.
- `dev` is the official branch for developpers, must be stable. It's a temporary state for the code, waiting to test in real large dataset before sending code on `main`.
- `<type>/<short_task_name>` is the template of the temporary branch for any development. `<type>` can be:
    - `feat` for development of a new feature
    - `fix` for bug fix
    - `doc` for documentation update
    - `test` for adding missing tests
    - `chore` for any technical task that has no impact on the code's functionality (to keep the project clear and organized)
        Can be dependancy update, bash script to make some build shortcut, files for continuous integration, ...

## Workflow

1. Create an issue to describe what you want to do (can be update, it can be used as a documentation)
2. From the latest version of `dev` branch, create a new branch named like `<type>/<short_task_name>`
3. Write a test (that must fail) to specify the feature behavior or reproduce the bug to fix.
4. Code the minimum to pass the test
5. Enjoy this advance by contemplating your code for a few moments...
    - Isn't there a simpler way of doing this?
    - Are all cases covered by the test?
    - Are you sure there are no hidden side effects?
    - If it's necessary, go back to the step 3.
6. Document what you've just done before you forget the subtleties.
7. Update the software version
8. Squash your commits to have clear and coherent messages.
9. Rebase your branch with `dev`
10. Test the whole software with real data
11. Rebase on `main`

## Update version number MAJOR.MINOR.PATCH

We following the [Semantic Versioning convention](https://semver.org/).

When a bugfix is rebase on dev branch, make "PATCH += 1".

When a new feature is rebase on dev branch, make "MINOR += 1".

When you break compatibility with the previous versions, make "MAJOR += 1".

There is three places where version number needs to be updated:
- `docs/source/conf.py::release`
- `pyproject.toml::project::version`
- `traceratops::_version.py::__version__`
