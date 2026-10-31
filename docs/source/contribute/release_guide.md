# Release guide

Traceratops uses Semantic Versioning and `vMAJOR.MINOR.PATCH` tags. GitHub releases
and PyPI publication are separate, maintainer-approved steps. The publishing
workflow never creates tags or changes branches.

## A. Prepare a release

1. Integrate approved feature-branch changes into `dev`. Review the actual code
   changes since the previous release and identify unfinished work to exclude.
2. Prepare and review the version in `pyproject.toml`, changelog/date and public
   release notes. Keep `[Unreleased]` for future changes and synchronize the
   documentation's changelog copy. Runtime and Sphinx read the installed version.
3. Run relevant scientific tests, the supported OS/Python CI matrix, dependency
   audits, pre-commit and documentation checks. Build the wheel and sdist, run
   `twine check`, and verify a clean wheel installation and console commands.
   Review any remaining warnings and coverage limits before approval.

The [2.0.0 assessment](release_2_0_0_checklist.md) is retained as historical
validation evidence, not a checklist for future releases.

## B. Integrate into main

Open a reviewed PR targeting `main`. If histories diverge, use an integration
branch from main, merge dev into it and resolve conflicts there. Compare the
result with both source branches and validate the actual integration result.

After approval, use the repository's normal **merge-commit** strategy to preserve
history. Verify the resulting version and release contents on main. Synchronize
main back into dev through a reviewed merge as appropriate; do not reset or
force-push the long-lived branches.

## C. Create the GitHub release

Confirm the exact approved commit on `main` and its package version. Create the
matching `vMAJOR.MINOR.PATCH` tag on that commit and publish the corresponding
GitHub release with the public release notes. Check for an existing tag/release
first: **never overwrite or move an existing release tag**.

A GitHub release does not automatically publish a package to PyPI.

## D. Publish to PyPI

Publication currently requires **manual dispatch** of `publish-pypi.yml` from
**main**, providing the existing release tag. In GitHub Actions:

1. Select **Build and publish Python distributions** and choose **Run workflow**.
2. Select branch **main** and enter the approved existing tag in `release_tag`.
3. Inspect the resolved tag/commit SHA, version and validated wheel/sdist
   artifacts, including filenames and hashes.
4. If the **pypi** environment requires approval, approve only when a real upload
   is intended: **environment approval authorizes publication to PyPI**, not a
   dry run. Inspect the publishing job's result.

The workflow validates the tag's format, main ancestry and package version, then
builds from that exact tagged source using infrastructure on current main. It
checks distributions and a clean wheel installation before a separate job
publishes with **Trusted Publishing**. **No PyPI API token is required.** PR runs
validate packages but cannot enter the publishing job; publication is not
triggered automatically by tags or GitHub releases.

## E. Verify publication

- Check the published version on PyPI, its Python requirement, package metadata,
  wheel/sdist files and attestations. Compare file hashes with the validated
  workflow artifacts and confirm the GitHub release points to the approved tag.
- In a fresh environment, install the exact published version and check the
  installed version, `trace_filter --help` and `python -m pip check`. Follow the
  [stable installation guide](../quickstart/installation.md) for pip, uv and
  Windows commands. Local wheel tests alone do not verify a PyPI installation.
- Confirm successful Actions checks and continued dependency security monitoring.

Published PyPI files cannot be overwritten. If an upload fails or only some files
appear, inspect the existing files, hashes and logs before retrying; do not
blindly rerun or enable `skip-existing`. A source correction generally requires
a new version and tag through the reviewed release process, never a moved tag.
