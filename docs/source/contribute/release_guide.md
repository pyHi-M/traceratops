# Release and PyPI publishing guide

Traceratops uses Semantic Versioning and `vMAJOR.MINOR.PATCH` release tags.
Stable packages are distributed through PyPI; development installations use an
editable `dev` checkout with uv. GitHub releases and PyPI publication are separate
steps. Every merge, tag/release creation and publication requires maintainer
approval; the workflow never creates or moves Git tags or branches.

## 1. Develop and validate on feature branches

Start from current `origin/dev`, create a dedicated feature branch, and install
with `uv pip install -e ".[dev]"`. Run scientific regressions, the supported
OS/Python matrix, dependency audits, documentation and packaging checks.
Open a PR into `dev`; do not change scientific algorithms to satisfy maintenance
checks. The distribution workflow also builds and smoke-tests PR source without
publishing or requesting an OIDC token.

## 2. Integrate into dev and prepare release metadata

After review, integrate the feature into `dev`. Review changes since the last
published release, including actual source changes and open/unmerged PRs. Prepare
the version, changelog/date and public release notes on a release branch and
review them into dev. Preserve `[Unreleased]` for future work. Only
`pyproject.toml` defines the version; runtime and Sphinx read installed metadata.
The [2.0.0 assessment](release_2_0_0_checklist.md) retains the historical validation
and scope decisions; its uncompleted historical checklist is not current branch
or publication state.

## 3. Integrate the release into main

Review dev → main through a PR. If histories diverge, create an integration
branch from main, normally merge dev into it, and resolve conflicts there while
preserving both histories. Compare the resulting tree with both source branches.
Require scientific tests, security audits, pre-commit, documentation and package
checks on the actual result. Do not assume main equals the last release tag.

After approval, merge into main with a normal merge commit. Merge the resulting
main history back into dev through a reviewed PR to avoid recurring divergence.
Do not reset/force-push either long-lived branch.

## 4. Create the Git tag and GitHub release

After explicit approval, identify the exact reviewed main commit and confirm its
package version. Check for an existing tag/release before creating anything.
Create `vMAJOR.MINOR.PATCH` on that approved main commit, then publish the GitHub
release using the public release notes. Do not replace a tag or rewrite a release
that already exists. GitHub source archives do not themselves publish to PyPI.

**v2.0.0 already exists** at
`904b31118d648089afe1d3a4b73b32d60d1f257d`, with a published GitHub release.
Do not repeat these creation steps for 2.0.0. The first PyPI publication uses the
existing tag as described below, after the publishing infrastructure reaches main.

## 5. One-time Trusted Publisher setup (manual)

The maintainer's personal PyPI account will initially own `traceratops`. In that
account, configure a **pending Trusted Publisher** for the not-yet-created
project. If the project has already been created by that account, configure its
Trusted Publisher instead. Verify name availability/ownership at registration;
stop if another project owns the name. A missing public project does not reserve
the name. No project or settings are created by this PR.

Use these exact values in PyPI's GitHub publisher form:

| Field | Value |
| --- | --- |
| PyPI project name | `traceratops` |
| GitHub owner | `pyHi-M` |
| GitHub repository | `traceratops` |
| Workflow filename | `publish-pypi.yml` (filename only) |
| GitHub environment | `pypi` |

Keep the owning account secure, with PyPI's required account authentication/2FA.
No long-lived PyPI API token or GitHub secret is needed. Pending publishers allow
the first successful authorized upload to create the project under that account.
They do not establish name availability in advance.

A repository maintainer must also create the GitHub environment **pypi** under
**Settings → Environments**, before publication, and configure:

- Deployment branches restricted to **main only**, not arbitrary branches or tags.
- Required reviewers and prevention of self-review where available; protect main
  and require review of workflow changes. Repository plan/settings may affect
  available protections.
- No PyPI username/password/API-token secret. The job uses short-lived OIDC
  credentials tied to this repository, workflow filename and environment.

These are manual maintainer actions requiring account/repository permissions.
The workflow's own event/ref/repository checks supplement environment protections;
they are not a replacement for them. Do not approve deployment from an unexpected
workflow revision. Pin changes to third-party Actions must be reviewed as code.

## 6. Publish an existing tag manually

Only an authorized maintainer should run publication, after checking the project,
publisher, environment settings and release approval. The workflow must already
be on default **main**. In GitHub Actions:

1. Select **Build and publish Python distributions** (`publish-pypi.yml`).
2. Choose **Run workflow**, select branch **main**, and enter the existing tag,
   for example **v2.0.0**.
3. Inspect the build logs, resolved commit SHA, version and distribution artifact.
4. Approve the protected **pypi** environment deployment if configured, after
   checking the requested tag and artifacts. Inspect the publishing result.

Equivalent CLI invocation, for use only after publication approval:

```bash
gh workflow run publish-pypi.yml --repo pyHi-M/traceratops --ref main -f release_tag=v2.0.0
```

The workflow accepts only stable `vMAJOR.MINOR.PATCH` tags without leading zeros.
It rejects missing tags, commits outside main ancestry, project/version mismatches
and malformed distributions. Manual dispatch from a different branch/repository
fails before building. It does not accept a free-form source SHA or create refs.

The **infrastructure checkout runs from current main**, while a **separate source
checkout uses the exact validated tag commit SHA**. Thus v2.0.0 can be published
even though its source commit predates this workflow and validation scripts.
Current main source is not substituted for the requested tag. The first 2.0.0
PyPI description and metadata therefore retain the README/license metadata from
the existing tag; these new documentation changes are not patched into old
release artifacts. Future releases can include updated metadata through normal
versioned source changes.

The build job has read-only repository permissions and no OIDC permission. It
runs validator tests, builds an sdist and a wheel from that sdist, runs Twine,
checks metadata/source payloads, then installs the wheel into a fresh Python
3.11 environment and verifies version and installed commands. It uploads only
the two validated distribution files as a GitHub Actions artifact.

The separate publishing job runs **only** for manual dispatch on canonical main,
uses environment `pypi`, and alone receives `id-token: write`. It downloads only
that run's artifact. It has no source checkout, build or dependency-install step.
The pinned official PyPA action publishes to PyPI using Trusted Publishing and
produces digital attestations. PRs never enter this job or consume prior PR
artifacts for publication. There is no automatic tag/release-triggered upload.

### First publication of v2.0.0

After this infrastructure has been reviewed into dev and integrated into main:

1. Register the pending publisher and configure/protect `pypi` manually.
2. Confirm existing v2.0.0 still resolves to
   `904b31118d648089afe1d3a4b73b32d60d1f257d` and remains on main ancestry.
3. Confirm PyPI has no existing 2.0.0 files and the intended account controls or
   can create the project. Resolve a naming conflict before proceeding.
4. After explicit publication approval, dispatch from **main** with
   `release_tag=v2.0.0`, inspect the build, and approve the environment deployment.
5. Perform the post-publication checks below. Do not recreate the GitHub release,
   move its tag or change its source/version as a publishing workaround.

## 7. Post-publication verification and failure handling

Check `https://pypi.org/project/traceratops/2.0.0/` (substitute the release version),
the owning account, Python requirement, wheel/sdist filenames, metadata and
attestations. Compare PyPI SHA-256 file hashes with the validated run artifacts.
In a **new environment**, test the actual published package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --index-url https://pypi.org/simple traceratops==2.0.0
python -c "from importlib.metadata import version; assert version('traceratops') == '2.0.0'"
trace_filter --help
python -m pip check
```

Use the [Windows/uv instructions](../quickstart/installation.md) as appropriate.
Before the first publication these index commands cannot succeed; local wheel
validation is not evidence that PyPI installation has already been tested.
After successful first-publication verification, remove the temporary pending
publication notices from README/installation docs through a reviewed PR.
Confirm native Actions failure notifications/Dependabot settings and continued
weekly dependency audits on main.

**PyPI filenames/versions cannot be overwritten.** The action uses
`skip-existing: false` so duplicates fail visibly. On failed or partial uploads,
inspect PyPI's actual files, hashes, attestations and workflow logs before any
retry. Do not blindly rerun, enable skip-existing, retag, or overwrite artifacts.
If some files were accepted, decide the recovery with the owning maintainer/PyPI
support; a corrected source generally requires a new version and tag through the
normal release process. Never move v2.0.0 to repair a publication failure.
