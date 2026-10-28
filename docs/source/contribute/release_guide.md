# Release guide

The project uses Semantic Versioning and `vMAJOR.MINOR.PATCH` tags. Prepare and
validate changes on `release/*`, review them into `dev`, then review `dev` into
`main`. The published tag must identify the approved release commit on `main`.
Do not merge into `main`, create a tag, or publish a release without explicit
maintainer approval. Do not publish packages to PyPI as part of this procedure.

## 1. Assess the release baseline

Fetch `main`, `dev` and tags. Inspect the latest **published GitHub release**,
resolve its tag to a commit, and compare that commit with both branches. Do not
assume the highest tag, published release and current `main` are identical.
Review merged history, relevant source diffs, open PRs and CI results.

For the proposed 2.0.0 release, see the
[assessment and validation checklist](release_2_0_0_checklist.md) and
[prepared release notes](../quickstart/release_2_0_0.md).

## 2. Prepare a release branch and PR into dev

Create `release/2.0.0` from current `origin/dev`. Update the authoritative version
in `pyproject.toml`, the root changelog, its documentation copy, and release notes.
`traceratops._version` and Sphinx obtain the installed package version through
package metadata; there is no second version constant to edit.

Run the existing scientific suite, the supported OS/Python CI matrix, dependency
audits, and package build/installation checks. Distinguish numerical regression
coverage from imports and CLI help checks. Do not change algorithms or defaults
to make checks pass. Open a release-preparation PR targeting `dev`; do not merge
it automatically. If further features are merged meanwhile, reassess the release
notes and rerun affected checks.

## 3. After approval, merge preparation into dev

After review, a maintainer merges the preparation PR into `dev`. Open the prepared
`dev` → `main` PR and confirm its actual commit range and release description.
The preparation PR's checks do not replace checks on this final integration PR.
There is no need to reset or force-push the temporary branch.

## 4. Verify the dev-to-main PR

Require successful test, security and existing hygiene checks. Confirm that the
version and release notes match the reviewed scope, no unapproved features have
entered, and the changelog covers all changes since the published release tag.
Keep the previous tag untouched. Confirm `v2.0.0` and its release do not already
exist. Resolve any merge conflicts on a reviewed branch, not directly on `main`.

## 5. Approval boundary: merge into main

Only after **explicit approval**, merge the reviewed `dev` → `main` PR, preferably
with the repository's existing merge-commit convention. Fetch `origin/main` and
verify that the merged commit contains the approved version and all validated
changes. If main has advanced again, identify the intended release commit
explicitly rather than tagging an unreviewed tip.

## 6. Approval boundary: tag the main release commit

Only after **explicit approval**, use the exact verified main commit SHA:

```bash
git fetch origin main --tags
git show <approved-main-sha>:pyproject.toml
git merge-base --is-ancestor <approved-main-sha> origin/main
git ls-remote --tags origin refs/tags/v2.0.0
# Stop if that tag already exists; do not overwrite it.
git tag -a v2.0.0 <approved-main-sha> -m "Release v2.0.0"
git push origin refs/tags/v2.0.0
```

Check package version 2.0.0 before creating the tag. Do not tag `dev` before the
main merge: a tag on dev need not identify the final main release commit.

## 7. Approval boundary: publish the GitHub release

After **explicit approval**, create the GitHub release using the existing
`v2.0.0` tag and the prepared release notes. Set the actual release date in the
changelog through the reviewed preparation/integration PR before tagging.
Do not duplicate or overwrite a tag/release. No automated release or package
publishing workflow is currently configured.

Build from the exact tagged source with `python -m build`, validate the sdist and
wheel with `twine check`, and install the wheel in a fresh environment. If binary
artifacts are to be attached, attach these tagged-source builds, not an earlier
local preparation build. GitHub source archives alone are not wheel artifacts.
Do not upload anything to a package index.

## 8. Verify publication and monitoring

- Confirm that the release is published (not draft), its tag resolves to the
  approved main SHA, and the package metadata declares 2.0.0 and Python >=3.11.
- Confirm the release title, notes, comparison link and any attached artifacts.
- Check the security audit workflow exists on `main`, the default branch.
  Weekly schedules run from the default branch; until it contains the workflow,
  dependency-related PR checks supply audit coverage. Once present, manually
  dispatch an audit and check the next weekly run.
- Enable native GitHub Actions failure notifications and Dependabot alerts as
  described in the [development guide](dev_installation.md#dependency-maintenance-and-security).
- Keep approval and CI evidence in the release PR. Delete temporary branches only
  after confirming publication; do not delete `dev` or change existing tags.
