#!/usr/bin/env python3
"""Collect one file per subdirectory using pattern matching.

Walks the immediate subdirectories of a root folder, searches recursively
for files whose name matches an example, and copies exactly one match per
subdirectory to a destination folder.

Two matching modes:

**Exact match** (no ``--variable-part``)::

    collect_files --root data/RUT \\
        --example-file "localizations_3D_barcode.dat" \\
        --copy-to collected/

Files are automatically renamed with the subdirectory name
(``localizations_3D_barcode_013_ROI.dat``) since they all share the
same name.

**Variable match** (with ``--variable-part``)::

    collect_files --root data/RUT \\
        --example-file "Trace_3D_barcode_mask-mask0_ROI-13.ecsv" \\
        --variable-part "13" \\
        --copy-to collected/

Matches ``ROI-14.ecsv`` but rejects ``ROI-021.ecsv`` (different length).
"""

from __future__ import annotations

from traceratops.script_banner import print_script_banner
import argparse
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Business logic
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FilePattern:
    """Filename pattern with an optional fixed-length variable zone.

    When *variable_length* is 0, the pattern performs exact matching
    (the candidate must equal *prefix*).  Otherwise a candidate matches
    when its total length equals ``len(prefix) + variable_length +
    len(suffix)`` and the prefix and suffix are identical.
    """

    prefix: str
    suffix: str
    variable_length: int

    @property
    def total_length(self) -> int:
        return len(self.prefix) + self.variable_length + len(self.suffix)

    def matches(self, filename: str) -> bool:
        """Return True if *filename* matches this pattern."""
        return (
            len(filename) == self.total_length
            and filename[: len(self.prefix)] == self.prefix
            and filename[len(self.prefix) + self.variable_length :] == self.suffix
        )


@dataclass
class ScanResult:
    """Outcome of scanning subdirectories for pattern matches."""

    unique: dict[str, Path] = field(default_factory=dict)
    ambiguous: dict[str, list[Path]] = field(default_factory=dict)
    empty: list[str] = field(default_factory=list)


def build_pattern(example_file: str, variable_part: str | None = None) -> FilePattern:
    """Derive a :class:`FilePattern` from *example_file*.

    When *variable_part* is ``None``, the pattern matches the exact
    filename.  Otherwise it matches filenames that are identical to
    *example_file* except at the position of *variable_part*, where any
    characters of the same length are accepted.

    Raises
    ------
    ValueError
        If *variable_part* is the empty string, or appears zero or more
        than one time in *example_file*.
    """
    if variable_part is None:
        return FilePattern(prefix=example_file, suffix="", variable_length=0)

    if not variable_part:
        raise ValueError("variable_part must not be empty")

    count = example_file.count(variable_part)
    if count == 0:
        raise ValueError(
            f"Variable part {variable_part!r} not found in "
            f"example file {example_file!r}"
        )
    if count > 1:
        raise ValueError(
            f"Variable part {variable_part!r} appears {count} times in "
            f"example file {example_file!r} (must appear exactly once)"
        )

    pos = example_file.index(variable_part)
    return FilePattern(
        prefix=example_file[:pos],
        suffix=example_file[pos + len(variable_part) :],
        variable_length=len(variable_part),
    )


def scan_subdirectories(root: Path, pattern: FilePattern) -> ScanResult:
    """Scan immediate subdirectories of *root* for files matching *pattern*.

    Each immediate child directory of *root* is an independent collection
    unit.  Files are searched **recursively** within each subdirectory.
    """
    result = ScanResult()

    subdirs = sorted(d for d in root.iterdir() if d.is_dir())
    for subdir in subdirs:
        matches = sorted(
            f for f in subdir.rglob("*") if f.is_file() and pattern.matches(f.name)
        )
        if len(matches) == 0:
            result.empty.append(subdir.name)
        elif len(matches) == 1:
            result.unique[subdir.name] = matches[0]
        else:
            result.ambiguous[subdir.name] = matches

    return result


def dest_name(source: Path, subdir_name: str, rename: bool) -> str:
    """Compute the destination filename for a collected file.

    When *rename* is True, the subdirectory name is inserted before the
    extension: ``data.ecsv`` from ``013_ROI`` becomes ``data_013_ROI.ecsv``.
    """
    if not rename:
        return source.name
    return f"{source.stem}_{subdir_name}{source.suffix}"


def check_collisions(
    files: dict[str, Path], copy_to: Path, rename: bool = False
) -> list[str]:
    """Return human-readable collision descriptions, if any.

    Checks for:
    * collisions between collected files (same destination name)
    * collisions with files already present in *copy_to*
    """
    problems: list[str] = []
    seen: dict[str, str] = {}  # dest filename -> subdir name

    for subdir_name, path in sorted(files.items()):
        name = dest_name(path, subdir_name, rename)
        if name in seen:
            problems.append(
                f"  {name} from {subdir_name!r} collides with "
                f"file from {seen[name]!r}"
            )
        seen[name] = subdir_name

        if (copy_to / name).exists():
            problems.append(
                f"  {name} from {subdir_name!r} collides with "
                f"existing file in {copy_to}"
            )

    return problems


def format_report(result: ScanResult) -> str:
    """Build a deterministic, sorted, human-readable scan report."""
    lines: list[str] = []

    if result.unique:
        lines.append(f"Matched ({len(result.unique)}):")
        for name in sorted(result.unique):
            lines.append(f"  {name} -> {result.unique[name]}")

    if result.empty:
        lines.append(f"No match ({len(result.empty)}):")
        for name in result.empty:
            lines.append(f"  {name}")

    if result.ambiguous:
        lines.append(f"Ambiguous ({len(result.ambiguous)}):")
        for name in sorted(result.ambiguous):
            paths = result.ambiguous[name]
            lines.append(f"  {name} ({len(paths)} matches):")
            for p in paths:
                lines.append(f"    {p}")

    return "\n".join(lines)


def collect_files(
    root: Path,
    example_file: str,
    copy_to: Path,
    variable_part: str | None = None,
    force: bool = False,
    rename: bool = False,
) -> list[Path]:
    """Collect one file per subdirectory and copy to *copy_to*.

    Parameters
    ----------
    root : Path
        Root directory whose immediate children are the collection units.
    example_file : str
        Example filename.  When *variable_part* is given it must contain
        exactly one occurrence of it; otherwise an exact match is used.
    copy_to : Path
        Destination folder (created if missing).
    variable_part : str or None
        The substring that varies across subdirectories.  When ``None``,
        exact matching is used and files are automatically renamed with
        the subdirectory name to avoid collisions.
    force : bool
        When True, subdirectories with no match are silently skipped.
    rename : bool
        When True, the subdirectory name is inserted before the file
        extension (e.g. ``data.dat`` from ``013_ROI`` becomes
        ``data_013_ROI.dat``).  Automatically enabled when
        *variable_part* is ``None``.

    Returns
    -------
    list[Path]
        Paths of the copied files in *copy_to*.

    Raises
    ------
    ValueError
        Bad arguments or ambiguous matches.
    FileNotFoundError
        Empty subdirectories without *force*, or nothing to copy.
    FileExistsError
        Name collisions in *copy_to*.
    """
    # Exact match → rename is mandatory (all files share the same name).
    if variable_part is None:
        rename = True

    pattern = build_pattern(example_file, variable_part)
    result = scan_subdirectories(root, pattern)

    report = format_report(result)
    if report:
        print(report)

    if result.ambiguous:
        print("\nDecision: nothing copied.")
        raise ValueError(
            f"{len(result.ambiguous)} subdirectory(ies) have multiple matches. "
            "Use a more specific --variable-part to disambiguate."
        )

    if result.empty and not force:
        print("\nDecision: nothing copied.")
        raise FileNotFoundError(
            f"{len(result.empty)} subdirectory(ies) have no matching files. "
            "Use --force to skip them."
        )

    if not result.unique:
        print("\nDecision: nothing copied.")
        raise FileNotFoundError("No files found to copy.")

    copy_to.mkdir(parents=True, exist_ok=True)

    collisions = check_collisions(result.unique, copy_to, rename)
    if collisions:
        print("\nDecision: nothing copied.")
        raise FileExistsError(
            "Name collisions in destination:\n" + "\n".join(collisions)
        )

    copied: list[Path] = []
    for subdir_name, source in sorted(result.unique.items()):
        dest = copy_to / dest_name(source, subdir_name, rename)
        shutil.copy2(source, dest)
        copied.append(dest)

    print(f"\nCopied {len(copied)} file(s) to {copy_to}")
    return copied


# ---------------------------------------------------------------------------
# CLI layer
# ---------------------------------------------------------------------------


def parse_arguments() -> argparse.ArgumentParser:
    """Build and return the argument parser (used by sphinx-argparse)."""
    parser = argparse.ArgumentParser(
        prog="collect_files",
        description=(
            "Collect exactly one file per subdirectory of ROOT by matching "
            "EXAMPLE_FILE.  Without --variable-part the match is exact and "
            "files are automatically renamed with the subdirectory name.  "
            "With --variable-part a fixed-length pattern match is used."
        ),
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
        help="Root directory to scan (immediate subdirectories are collection units)",
    )
    parser.add_argument(
        "--example-file",
        required=True,
        help="Example filename to search for in each subdirectory",
    )
    parser.add_argument(
        "--variable-part",
        default=None,
        help=(
            "Substring in EXAMPLE_FILE that varies across subdirectories. "
            "When omitted, exact filename matching is used and files are "
            "automatically renamed with the subdirectory name."
        ),
    )
    parser.add_argument(
        "--copy-to",
        type=Path,
        required=True,
        help="Destination folder for collected files (created if missing)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Skip subdirectories with no match instead of failing",
    )
    parser.add_argument(
        "--rename",
        action="store_true",
        help=(
            "Rename collected files by inserting the subdirectory name "
            "before the extension (e.g. data.dat from 013_ROI -> "
            "data_013_ROI.dat). Automatic when --variable-part is omitted."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point.  Returns 0 on success, 1 on business-logic failure."""
    print_script_banner(__file__, __doc__)
    args = parse_arguments().parse_args(argv)

    if not args.root.is_dir():
        print(f"Error: {args.root} is not a directory", file=sys.stderr)
        return 1

    try:
        collect_files(
            root=args.root,
            example_file=args.example_file,
            copy_to=args.copy_to,
            variable_part=args.variable_part,
            force=args.force,
            rename=args.rename,
        )
    except (ValueError, FileNotFoundError, FileExistsError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
