# Release guide

Traceratops uses Semantic Versioning and `vMAJOR.MINOR.PATCH` tags. GitHub releases
automatically start the PyPI publishing workflow; a protected environment can
require maintainer approval before upload. The workflow never creates tags or
changes branches.

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

Only publishing a **stable, non-draft, non-prerelease** GitHub release triggers
automatic publication. Creating a tag alone, editing a release or merging workflow
changes does not trigger an upload or republish the existing v2.0.0 package.

## D. Publish to PyPI

The `release.published` event starts **Build and publish Python distributions**
(`publish-pypi.yml`). No second manual launch is needed for ordinary releases.
The workflow resolves the existing stable tag, checks its commit against the
release event, verifies main ancestry and package version, then builds from that
exact commit using validation infrastructure from current main. PR runs only
build and validate; they cannot publish or receive publishing OIDC permissions.

Inspect the run summary's tag, commit SHA and version, plus the validated
wheel/sdist artifacts and hashes. If the **pypi** environment requests approval,
**approval authorizes a real PyPI upload**, not a dry run. Confirm that the tag
still identifies the approved commit before approving. A separate job publishes
only that run's validated artifacts using **Trusted Publishing**; no PyPI API
token is required. The release is available on PyPI only after that job succeeds.

### Manual fallback

For an exceptional case, select this workflow in GitHub Actions, choose **Run
workflow**, select branch **main** and enter the existing stable tag in
`release_tag`. Review the same validation and approval evidence. Manual and
automatic runs use the same tag validator and concurrency group. Do not dispatch
while another attempt is pending or blindly rerun a failed upload.

### Environment deployment policy

Before enabling automatic publication, a repository administrator must open
**Settings → Environments → pypi → Deployment branches and tags**. With
**Selected branches and tags**, retain the **Branch** rule `main` for manual runs
and add a **Tag** rule `v*` for release events. A branch rule alone does not allow
tag refs. Retain required reviewers and all existing protection rules; this PR
does not change settings or the existing Trusted Publisher registration.

The `v*` deployment rule is broader than the stable-tag validator: validation
still rejects prereleases and commits outside main. Keep release-tag creation
and workflow changes restricted to trusted maintainers and review deployments;
the environment pattern alone does not prove main ancestry. Future release tags
must include this workflow version for automatic triggering; manual dispatch
from current main remains available for older tags.

Never move release tags. A tag moved before validation fails against the release
event commit. Builds remain pinned to the validated SHA even if a tag changes
later; reviewers must recheck it before approval. Runs for the same tag are
serialized without cancellation, but duplicate events or manual retries can
still attempt a second upload. Existing files fail visibly (`skip-existing:
false`); there is no automatic upload retry.

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
