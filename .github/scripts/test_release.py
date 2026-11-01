"""Negative release checks use mocks/temporary artifacts, never repository tags."""

import io
import json
import os
import subprocess
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import release


class TagValidationTests(unittest.TestCase):
    def test_invalid_or_injected_tags_are_rejected_before_git(self):
        for tag in (
            "2.0.0",
            "v02.0.0",
            "v2.0",
            "v2.0.0rc1",
            "v2.0.0; echo bad",
            "../main",
        ):
            with self.subTest(tag=tag), patch.object(release, "git") as git:
                with self.assertRaises(ValueError):
                    release.validate_tag(tag)
                git.assert_not_called()

    def test_missing_tag_is_rejected(self):
        with patch.object(
            release, "git", side_effect=subprocess.CalledProcessError(128, "git")
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                release.validate_tag("v2.0.0")

    def test_non_main_ancestry_is_rejected(self):
        with (
            patch.object(release, "git", return_value="a" * 40),
            patch.object(
                release.subprocess,
                "run",
                side_effect=subprocess.CalledProcessError(1, "git"),
            ),
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                release.validate_tag("v2.0.0")

    def test_version_or_project_mismatch_is_rejected(self):
        for name, version in (("traceratops", "2.0.1"), ("different", "2.0.0")):
            project = f'[project]\nname="{name}"\nversion="{version}"\n'
            with (
                patch.object(release, "git", side_effect=["a" * 40, project]),
                patch.object(release.subprocess, "run"),
                self.assertRaises(ValueError),
            ):
                release.validate_tag("v2.0.0")

    def test_exact_existing_tag_commit_and_metadata_are_used(self):
        project = '[project]\nname="traceratops"\nversion="2.0.0"\n'
        with (
            patch.object(release, "git", side_effect=["a" * 40, project]) as git,
            patch.object(release.subprocess, "run") as ancestry,
        ):
            self.assertEqual(
                release.validate_tag("v2.0.0"), {"sha": "a" * 40, "version": "2.0.0"}
            )
            self.assertEqual(
                git.call_args_list[0].args,
                ("rev-parse", "--verify", "refs/tags/v2.0.0^{commit}"),
            )
            self.assertEqual(
                git.call_args_list[1].args, ("show", f"{'a' * 40}:pyproject.toml")
            )
            self.assertIn("refs/remotes/origin/main", ancestry.call_args.args[0])


class EventSelectionTests(unittest.TestCase):
    """Safe event harness: real authorization/tag logic, mocked read-only Git calls."""

    def setUp(self):
        self.sha = "a" * 40
        self.project = '[project]\nname="traceratops"\nversion="2.0.0"\n'
        self.event = {
            "action": "published",
            "release": {"tag_name": "v2.0.0", "draft": False, "prerelease": False},
        }

    def select(
        self,
        name="release",
        ref="refs/tags/v2.0.0",
        event=None,
        repository="pyHi-M/traceratops",
    ):
        return release.select_event(
            name, repository, ref, self.sha, self.event if event is None else event
        )

    def test_release_tag_ref_resolves_real_tag_not_payload_target(self):
        self.event["release"]["target_commitish"] = "dev"
        with (
            patch.object(release, "git", side_effect=[self.sha, self.project]) as git,
            patch.object(release.subprocess, "run") as ancestry,
        ):
            result = self.select()
        self.assertEqual(
            result,
            {
                "publish": "true",
                "release_tag": "v2.0.0",
                "sha": self.sha,
                "version": "2.0.0",
            },
        )
        self.assertEqual(
            git.call_args_list[0].args,
            ("rev-parse", "--verify", "refs/tags/v2.0.0^{commit}"),
        )
        ancestry.assert_called_once_with(
            [
                "git",
                "merge-base",
                "--is-ancestor",
                self.sha,
                "refs/remotes/origin/main",
            ],
            check=True,
        )

    def test_manual_main_resolves_tag_independently_of_infrastructure_sha(self):
        with (
            patch.object(release, "git", side_effect=["b" * 40, self.project]),
            patch.object(release.subprocess, "run"),
        ):
            result = self.select(
                "workflow_dispatch",
                "refs/heads/main",
                {"inputs": {"release_tag": "v2.0.0"}},
            )
        self.assertEqual(result["sha"], "b" * 40)
        self.assertEqual(result["publish"], "true")

    def test_pr_workflow_changes_cannot_republish_existing_v2(self):
        with patch.object(release, "validate_tag") as validate:
            result = self.select("pull_request", "refs/pull/123/merge")
        self.assertEqual(result, {"publish": "false", "sha": self.sha})
        validate.assert_not_called()

    def test_non_publication_events_and_edited_existing_release_rejected(self):
        for name in ("push", "schedule", "release"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.select(name, event={**self.event, "action": "edited"})

    def test_manual_dev_and_unexpected_repository_rejected_before_tag(self):
        with patch.object(release, "validate_tag") as validate:
            with self.assertRaises(ValueError):
                self.select("workflow_dispatch", "refs/heads/dev")
            for name in ("release", "workflow_dispatch"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    self.select(name, "refs/heads/main", repository="other/traceratops")
        validate.assert_not_called()

    def test_draft_prerelease_and_incomplete_payload_rejected(self):
        for changes in ({"draft": True}, {"prerelease": True}, {"draft": None}):
            event = {**self.event, "release": {**self.event["release"], **changes}}
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.select(event=event)

    def test_release_branch_ref_and_malformed_event_sha_rejected(self):
        with self.assertRaises(ValueError):
            self.select(ref="refs/heads/dev")
        self.sha = "untrusted-target"
        with self.assertRaises(ValueError):
            self.select()

    def test_invalid_tag_for_both_publication_paths(self):
        for name, ref, event in (
            (
                "release",
                "refs/tags/v2.0.0rc1",
                {
                    **self.event,
                    "release": {**self.event["release"], "tag_name": "v2.0.0rc1"},
                },
            ),
            (
                "workflow_dispatch",
                "refs/heads/main",
                {"inputs": {"release_tag": "v2.0.0rc1"}},
            ),
        ):
            with self.subTest(name=name), patch.object(release, "git") as git:
                with self.assertRaises(ValueError):
                    self.select(name, ref, event)
                git.assert_not_called()

    def test_missing_and_moved_release_tags_fail_closed(self):
        with (
            patch.object(
                release, "git", side_effect=subprocess.CalledProcessError(128, "git")
            ),
            self.assertRaises(subprocess.CalledProcessError),
        ):
            self.select()
        with (
            patch.object(release, "git", return_value="b" * 40),
            self.assertRaisesRegex(ValueError, "moved"),
        ):
            self.select()

    def test_outside_main_and_version_mismatch_rejected_via_event(self):
        with (
            patch.object(release, "git", return_value=self.sha),
            patch.object(
                release.subprocess,
                "run",
                side_effect=subprocess.CalledProcessError(1, "git"),
            ),
            self.assertRaises(subprocess.CalledProcessError),
        ):
            self.select()
        with (
            patch.object(
                release,
                "git",
                side_effect=[self.sha, self.project.replace("2.0.0", "2.0.1")],
            ),
            patch.object(release.subprocess, "run"),
            self.assertRaises(ValueError),
        ):
            self.select()

    def test_event_cli_reads_json_context_and_writes_nonpublishing_pr_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = Path(directory) / "event.json"
            output = Path(directory) / "output"
            payload.write_text(json.dumps(self.event))
            environment = {
                "GITHUB_EVENT_NAME": "pull_request",
                "GITHUB_REPOSITORY": "pyHi-M/traceratops",
                "GITHUB_REF": "refs/pull/123/merge",
                "GITHUB_SHA": self.sha,
                "GITHUB_EVENT_PATH": str(payload),
            }
            with (
                patch.dict(os.environ, environment),
                patch(
                    "sys.argv", ["release.py", "event", "--github-output", str(output)]
                ),
                patch("sys.stdout", new_callable=io.StringIO),
            ):
                release.main()
            self.assertEqual(output.read_text(), f"publish=false\nsha={self.sha}\n")


class DistributionValidationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.source = Path(self.temporary.name) / "source"
        self.dist = Path(self.temporary.name) / "dist"
        (self.source / "traceratops").mkdir(parents=True)
        self.dist.mkdir()
        (self.source / "pyproject.toml").write_text(
            '[project]\nname="traceratops"\nversion="2.0.0"\nrequires-python=">=3.11"\n'
        )
        for file in ("README.md", "LICENSE", "traceratops/__init__.py"):
            (self.source / file).write_text("example\n")
        self.metadata = (
            b"Name: traceratops\nVersion: 2.0.0\nRequires-Python: >=3.11\n\n"
        )
        self.write_artifacts()

    def write_artifacts(self, wheel_metadata=None, package_bytes=b"example\n"):
        with zipfile.ZipFile(
            self.dist / "traceratops-2.0.0-py3-none-any.whl", "w"
        ) as wheel:
            wheel.writestr(
                "traceratops-2.0.0.dist-info/METADATA", wheel_metadata or self.metadata
            )
            wheel.writestr("traceratops/__init__.py", package_bytes)
        with tarfile.open(self.dist / "traceratops-2.0.0.tar.gz", "w:gz") as sdist:
            for file in self.source.rglob("*"):
                if file.is_file():
                    sdist.add(
                        file,
                        arcname="traceratops-2.0.0/"
                        + file.relative_to(self.source).as_posix(),
                    )
            info = tarfile.TarInfo("traceratops-2.0.0/PKG-INFO")
            info.size = len(self.metadata)
            sdist.addfile(info, io.BytesIO(self.metadata))

    def test_matching_artifacts_pass(self):
        self.assertEqual(
            release.validate_distributions(self.source, self.dist), {"version": "2.0.0"}
        )

    def test_wrong_wheel_version_is_rejected(self):
        self.write_artifacts(wheel_metadata=self.metadata.replace(b"2.0.0", b"2.0.1"))
        with self.assertRaises(ValueError):
            release.validate_distributions(self.source, self.dist)

    def test_altered_package_payload_is_rejected(self):
        self.write_artifacts(package_bytes=b"altered\n")
        with self.assertRaises(ValueError):
            release.validate_distributions(self.source, self.dist)

    def test_unexpected_extra_artifact_is_rejected(self):
        (self.dist / "unexpected.txt").write_text("not a distribution")
        with self.assertRaises(ValueError):
            release.validate_distributions(self.source, self.dist)


if __name__ == "__main__":
    unittest.main()
