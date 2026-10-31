"""Negative release checks use mocks/temporary artifacts, never repository tags."""

import io
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
