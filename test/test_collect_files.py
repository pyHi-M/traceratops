"""Tests for collect_files module."""

from __future__ import annotations

from pathlib import Path

import pytest

from traceratops.collect_files import (
    FilePattern,
    build_pattern,
    collect_files,
    dest_name,
    main,
    scan_subdirectories,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EXAMPLE = "Trace_3D_barcode_mask-mask0_ROI-13.ecsv"
VARPART = "13"
LOC_EXAMPLE = "localizations_3D_barcode.ecsv"


def _make_tree(tmp_path: Path, layout: dict[str, list[str]]) -> Path:
    """Create a directory tree under tmp_path/root.

    *layout* maps subdirectory names to lists of filenames.  Filenames
    containing ``/`` create nested subdirectories automatically.

    Returns the root directory.
    """
    root = tmp_path / "root"
    root.mkdir()
    for subdir_name, filenames in layout.items():
        subdir = root / subdir_name
        subdir.mkdir(parents=True, exist_ok=True)
        for fname in filenames:
            fpath = subdir / fname
            fpath.parent.mkdir(parents=True, exist_ok=True)
            fpath.write_text(f"content of {fname}")
    return root


# ---------------------------------------------------------------------------
# FilePattern unit tests
# ---------------------------------------------------------------------------


class TestFilePattern:
    def test_matches_exact_variable(self):
        p = FilePattern(prefix="ROI-", suffix=".ecsv", variable_length=2)
        assert p.matches("ROI-13.ecsv")
        assert p.matches("ROI-99.ecsv")
        assert p.matches("ROI-ab.ecsv")

    def test_rejects_different_length(self):
        """Variable zone must have exactly the same length."""
        p = FilePattern(prefix="ROI-", suffix=".ecsv", variable_length=2)
        assert not p.matches("ROI-021.ecsv")  # 3 chars instead of 2
        assert not p.matches("ROI-1.ecsv")  # 1 char instead of 2

    def test_rejects_wrong_prefix(self):
        p = FilePattern(prefix="ROI-", suffix=".ecsv", variable_length=2)
        assert not p.matches("FOO-13.ecsv")

    def test_rejects_wrong_suffix(self):
        p = FilePattern(prefix="ROI-", suffix=".ecsv", variable_length=2)
        assert not p.matches("ROI-13.csv")

    def test_exact_match(self):
        """variable_length=0 means exact string matching."""
        p = FilePattern(prefix="data.dat", suffix="", variable_length=0)
        assert p.matches("data.dat")
        assert not p.matches("data.txt")
        assert not p.matches("xdata.dat")


# ---------------------------------------------------------------------------
# build_pattern
# ---------------------------------------------------------------------------


class TestBuildPattern:
    def test_nominal(self):
        p = build_pattern(EXAMPLE, VARPART)
        assert p.prefix == "Trace_3D_barcode_mask-mask0_ROI-"
        assert p.suffix == ".ecsv"
        assert p.variable_length == 2

    def test_exact_match_when_no_variable_part(self):
        p = build_pattern(LOC_EXAMPLE)
        assert p.prefix == LOC_EXAMPLE
        assert p.suffix == ""
        assert p.variable_length == 0
        assert p.matches(LOC_EXAMPLE)
        assert not p.matches("other.dat")

    def test_variable_not_found(self):
        with pytest.raises(ValueError, match="not found"):
            build_pattern(EXAMPLE, "ZZZZZ")

    def test_variable_appears_multiple_times(self):
        with pytest.raises(ValueError, match="appears 2 times"):
            # "3" appears twice in the example filename ("3D" and "13")
            build_pattern(EXAMPLE, "3")

    def test_empty_variable_part(self):
        with pytest.raises(ValueError, match="must not be empty"):
            build_pattern(EXAMPLE, "")


# ---------------------------------------------------------------------------
# scan_subdirectories
# ---------------------------------------------------------------------------


class TestScanSubdirectories:
    def test_unique_matches(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["Trace_3D_barcode_mask-mask0_ROI-15.ecsv"],
            },
        )
        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        assert len(result.unique) == 2
        assert not result.empty
        assert not result.ambiguous

    def test_empty_subdir(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["unrelated.txt"],
            },
        )
        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        assert len(result.unique) == 1
        assert result.empty == ["sub_b"]

    def test_ambiguous_subdir(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "sub_a": [
                    "Trace_3D_barcode_mask-mask0_ROI-14.ecsv",
                    "Trace_3D_barcode_mask-mask0_ROI-15.ecsv",
                ],
            },
        )
        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        assert not result.unique
        assert "sub_a" in result.ambiguous
        assert len(result.ambiguous["sub_a"]) == 2

    def test_recursive_search(self, tmp_path: Path):
        """Files nested deep inside a subdirectory are found."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a/deep/nested": [
                    "Trace_3D_barcode_mask-mask0_ROI-14.ecsv",
                ],
            },
        )
        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        assert "sub_a" in result.unique

    def test_only_immediate_subdirs_as_units(self, tmp_path: Path):
        """Nested directories are NOT treated as separate collection units."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        # Create a nested dir with its own match — should NOT appear as
        # a separate unit, only sub_a is a unit.
        nested = root / "sub_a" / "inner"
        nested.mkdir()
        (nested / "Trace_3D_barcode_mask-mask0_ROI-22.ecsv").write_text("x")

        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        # sub_a has two matches (one at top, one nested) → ambiguous
        assert "sub_a" in result.ambiguous
        # "inner" should NOT appear as a separate collection unit
        assert "inner" not in result.unique
        assert "inner" not in result.ambiguous
        assert "inner" not in result.empty

    def test_length_mismatch_not_matched(self, tmp_path: Path):
        """A file with a different-length variable zone is not matched."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": [
                    "Trace_3D_barcode_mask-mask0_ROI-021.ecsv",  # 3 chars
                ],
            },
        )
        pattern = build_pattern(EXAMPLE, VARPART)
        result = scan_subdirectories(root, pattern)

        assert not result.unique
        assert result.empty == ["sub_a"]

    def test_exact_match_scan(self, tmp_path: Path):
        """Exact-match pattern finds the right file in subdirectories."""
        root = _make_tree(
            tmp_path,
            {
                "013_ROI/localize_3d/data": [LOC_EXAMPLE],
                "014_ROI/localize_3d/data": [LOC_EXAMPLE],
            },
        )
        pattern = build_pattern(LOC_EXAMPLE)
        result = scan_subdirectories(root, pattern)

        assert len(result.unique) == 2
        assert not result.empty


# ---------------------------------------------------------------------------
# dest_name
# ---------------------------------------------------------------------------


class TestDestName:
    def test_no_rename(self):
        p = Path("/some/dir/data.ecsv")
        assert dest_name(p, "013_ROI", rename=False) == "data.ecsv"

    def test_rename(self):
        p = Path("/some/dir/data.ecsv")
        assert dest_name(p, "013_ROI", rename=True) == "data_013_ROI.ecsv"

    def test_rename_localization_ecsv(self):
        p = Path("/x/localizations_3D_barcode.ecsv")
        assert (
            dest_name(p, "013_ROI", rename=True)
            == "localizations_3D_barcode_013_ROI.ecsv"
        )


# ---------------------------------------------------------------------------
# collect_files — variable-part mode (integration tests)
# ---------------------------------------------------------------------------


class TestCollectFilesVariableMode:
    def test_nominal_success(self, tmp_path: Path):
        """1. Nominal success with a unique match per subdirectory."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["Trace_3D_barcode_mask-mask0_ROI-15.ecsv"],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        assert len(copied) == 2
        assert dest.is_dir()
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").exists()
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-15.ecsv").exists()

    def test_fail_empty_subdir_no_force(self, tmp_path: Path):
        """2. Failure when a subdirectory has no match and --force is off."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["unrelated.txt"],
            },
        )
        dest = tmp_path / "output"

        with pytest.raises(FileNotFoundError, match="no matching files"):
            collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        # Nothing was copied
        assert not dest.exists() or not list(dest.iterdir())

    def test_success_with_force(self, tmp_path: Path):
        """3. Success with --force despite empty subdirectories."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["unrelated.txt"],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, EXAMPLE, dest, variable_part=VARPART, force=True)

        assert len(copied) == 1
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").exists()

    def test_fail_ambiguous(self, tmp_path: Path):
        """4. Failure when a subdirectory has multiple matches."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": [
                    "Trace_3D_barcode_mask-mask0_ROI-14.ecsv",
                    "Trace_3D_barcode_mask-mask0_ROI-15.ecsv",
                ],
            },
        )
        dest = tmp_path / "output"

        with pytest.raises(ValueError, match="multiple matches"):
            collect_files(root, EXAMPLE, dest, variable_part=VARPART)

    def test_fail_variable_not_found(self, tmp_path: Path):
        """5. Failure when variable_part is not in example_file."""
        root = _make_tree(tmp_path, {"sub_a": ["anything.txt"]})
        dest = tmp_path / "output"

        with pytest.raises(ValueError, match="not found"):
            collect_files(root, EXAMPLE, dest, variable_part="ZZZZZ")

    def test_fail_variable_multiple_occurrences(self, tmp_path: Path):
        """6. Failure when variable_part appears multiple times."""
        root = _make_tree(tmp_path, {"sub_a": ["anything.txt"]})
        dest = tmp_path / "output"

        with pytest.raises(ValueError, match="appears 2 times"):
            collect_files(root, EXAMPLE, dest, variable_part="3")

    def test_no_match_different_length(self, tmp_path: Path):
        """7. No match when the variable zone has a different length."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-021.ecsv"],
            },
        )
        dest = tmp_path / "output"

        with pytest.raises(FileNotFoundError):
            collect_files(root, EXAMPLE, dest, variable_part=VARPART)

    def test_fail_collision_existing_file(self, tmp_path: Path):
        """8. Failure when a collected file collides with one in copy-to."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        dest = tmp_path / "output"
        dest.mkdir()
        (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").write_text("old")

        with pytest.raises(FileExistsError, match="collisions"):
            collect_files(root, EXAMPLE, dest, variable_part=VARPART)

    def test_creates_copy_to(self, tmp_path: Path):
        """9. Destination directory is created automatically."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        dest = tmp_path / "deep" / "nested" / "output"
        assert not dest.exists()

        collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        assert dest.is_dir()
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").exists()

    def test_recursive_search_in_subdirs(self, tmp_path: Path):
        """10. Files deep inside a subdirectory are found and copied."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a/tracing/data": [
                    "Trace_3D_barcode_mask-mask0_ROI-14.ecsv",
                ],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        assert len(copied) == 1

    def test_only_immediate_subdirs(self, tmp_path: Path):
        """11. Only immediate children of root are collection units."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        # A file directly in root should be ignored (root itself is not a unit)
        (root / "Trace_3D_barcode_mask-mask0_ROI-99.ecsv").write_text("x")

        dest = tmp_path / "output"
        copied = collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        assert len(copied) == 1
        names = [p.name for p in copied]
        assert "Trace_3D_barcode_mask-mask0_ROI-14.ecsv" in names
        assert "Trace_3D_barcode_mask-mask0_ROI-99.ecsv" not in names

    def test_metadata_preserved(self, tmp_path: Path):
        """Copy preserves file content (shutil.copy2)."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        source = root / "sub_a" / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv"
        source.write_text("precious data")

        dest = tmp_path / "output"
        collect_files(root, EXAMPLE, dest, variable_part=VARPART)

        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").read_text() == (
            "precious data"
        )

    def test_collision_between_collected_files(self, tmp_path: Path):
        """Two subdirs yielding the same basename should fail."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b/nested": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
            },
        )
        dest = tmp_path / "output"

        with pytest.raises(FileExistsError, match="collisions"):
            collect_files(root, EXAMPLE, dest, variable_part=VARPART)

    def test_explicit_rename_with_variable_part(self, tmp_path: Path):
        """--rename works with --variable-part too."""
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["Trace_3D_barcode_mask-mask0_ROI-15.ecsv"],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, EXAMPLE, dest, variable_part=VARPART, rename=True)

        assert len(copied) == 2
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14_sub_a.ecsv").exists()
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-15_sub_b.ecsv").exists()


# ---------------------------------------------------------------------------
# collect_files — exact-match mode (no variable_part)
# ---------------------------------------------------------------------------


class TestCollectFilesExactMode:
    def test_nominal_success(self, tmp_path: Path):
        """Exact match collects and auto-renames identical filenames."""
        root = _make_tree(
            tmp_path,
            {
                "013_ROI/localize_3d/data": [LOC_EXAMPLE],
                "014_ROI/localize_3d/data": [LOC_EXAMPLE],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, LOC_EXAMPLE, dest)

        assert len(copied) == 2
        assert (dest / "localizations_3D_barcode_013_ROI.ecsv").exists()
        assert (dest / "localizations_3D_barcode_014_ROI.ecsv").exists()

    def test_preserves_content(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {"sub_a/deep": [LOC_EXAMPLE]},
        )
        source = root / "sub_a" / "deep" / LOC_EXAMPLE
        source.write_text("important data")

        dest = tmp_path / "output"
        collect_files(root, LOC_EXAMPLE, dest)

        assert (dest / "localizations_3D_barcode_sub_a.ecsv").read_text() == (
            "important data"
        )

    def test_collision_with_existing_file(self, tmp_path: Path):
        """Collision with pre-existing file in dest is caught."""
        root = _make_tree(
            tmp_path,
            {"sub_a/deep": [LOC_EXAMPLE]},
        )
        dest = tmp_path / "output"
        dest.mkdir()
        (dest / "localizations_3D_barcode_sub_a.ecsv").write_text("old")

        with pytest.raises(FileExistsError, match="collisions"):
            collect_files(root, LOC_EXAMPLE, dest)

    def test_rejects_wrong_filename(self, tmp_path: Path):
        """Exact match does not match similar but different filenames."""
        root = _make_tree(
            tmp_path,
            {"sub_a": ["localizations_3D_barcode.ecsv.bak"]},
        )
        dest = tmp_path / "output"

        with pytest.raises(FileNotFoundError):
            collect_files(root, LOC_EXAMPLE, dest)

    def test_recursive_search(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {"sub_a/deep/nested": [LOC_EXAMPLE]},
        )
        dest = tmp_path / "output"

        copied = collect_files(root, LOC_EXAMPLE, dest)
        assert len(copied) == 1

    def test_force_skips_empty_subdirs(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "013_ROI/data": [LOC_EXAMPLE],
                "014_ROI": ["unrelated.txt"],
            },
        )
        dest = tmp_path / "output"

        copied = collect_files(root, LOC_EXAMPLE, dest, force=True)
        assert len(copied) == 1

    def test_fail_empty_subdir_no_force(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "013_ROI/data": [LOC_EXAMPLE],
                "014_ROI": ["unrelated.txt"],
            },
        )
        dest = tmp_path / "output"

        with pytest.raises(FileNotFoundError, match="no matching files"):
            collect_files(root, LOC_EXAMPLE, dest)


# ---------------------------------------------------------------------------
# CLI (main) tests
# ---------------------------------------------------------------------------


class TestMain:
    def test_variable_part_mode(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {"sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"]},
        )
        dest = tmp_path / "output"
        code = main(
            [
                "--root",
                str(root),
                "--example-file",
                EXAMPLE,
                "--variable-part",
                VARPART,
                "--copy-to",
                str(dest),
            ]
        )
        assert code == 0
        assert (dest / "Trace_3D_barcode_mask-mask0_ROI-14.ecsv").exists()

    def test_exact_mode_no_variable_part(self, tmp_path: Path):
        """CLI without --variable-part uses exact match + auto rename."""
        root = _make_tree(
            tmp_path,
            {
                "013_ROI/data": [LOC_EXAMPLE],
                "014_ROI/data": [LOC_EXAMPLE],
            },
        )
        dest = tmp_path / "output"
        code = main(
            [
                "--root",
                str(root),
                "--example-file",
                LOC_EXAMPLE,
                "--copy-to",
                str(dest),
            ]
        )
        assert code == 0
        assert (dest / "localizations_3D_barcode_013_ROI.ecsv").exists()
        assert (dest / "localizations_3D_barcode_014_ROI.ecsv").exists()

    def test_failure_exit_code(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {"sub_a": ["unrelated.txt"]},
        )
        dest = tmp_path / "output"
        code = main(
            [
                "--root",
                str(root),
                "--example-file",
                EXAMPLE,
                "--variable-part",
                VARPART,
                "--copy-to",
                str(dest),
            ]
        )
        assert code == 1

    def test_force_flag(self, tmp_path: Path):
        root = _make_tree(
            tmp_path,
            {
                "sub_a": ["Trace_3D_barcode_mask-mask0_ROI-14.ecsv"],
                "sub_b": ["unrelated.txt"],
            },
        )
        dest = tmp_path / "output"
        code = main(
            [
                "--root",
                str(root),
                "--example-file",
                EXAMPLE,
                "--variable-part",
                VARPART,
                "--copy-to",
                str(dest),
                "--force",
            ]
        )
        assert code == 0

    def test_nonexistent_root(self, tmp_path: Path):
        code = main(
            [
                "--root",
                str(tmp_path / "nope"),
                "--example-file",
                EXAMPLE,
                "--copy-to",
                str(tmp_path / "output"),
            ]
        )
        assert code == 1
