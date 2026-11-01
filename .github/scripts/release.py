"""Validate release inputs and artifacts; never create refs or upload packages."""

import argparse
import importlib.metadata
import json
import re
import subprocess
import sys
import sysconfig
import tarfile
import tomllib
import zipfile
from email.parser import BytesParser
from pathlib import Path

TAG_PATTERN = r"v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)"


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def validate_tag(tag):
    if not re.fullmatch(TAG_PATTERN, tag):
        raise ValueError("Release tag must be vMAJOR.MINOR.PATCH without leading zeros")
    sha = git("rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}")
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", sha, "refs/remotes/origin/main"],
        check=True,
    )
    project = tomllib.loads(git("show", f"{sha}:pyproject.toml"))["project"]
    if project["name"] != "traceratops" or project["version"] != tag[1:]:
        raise ValueError("Tag version and tagged traceratops package metadata disagree")
    return {"sha": sha, "version": project["version"]}


def validate_distributions(source, dist):
    source, dist = Path(source), Path(dist)
    project = tomllib.loads((source / "pyproject.toml").read_text())["project"]
    version = project["version"]
    if project["name"] != "traceratops" or not re.fullmatch(TAG_PATTERN, "v" + version):
        raise ValueError("Expected a stable traceratops version")
    wheels, sdists = list(dist.glob("*.whl")), list(dist.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1 or len(list(dist.iterdir())) != 2:
        raise ValueError("Expected exactly one wheel and one source distribution")
    with zipfile.ZipFile(wheels[0]) as wheel, tarfile.open(sdists[0]) as sdist:
        metadata_paths = [
            n for n in wheel.namelist() if n.endswith(".dist-info/METADATA")
        ]
        if len(metadata_paths) != 1:
            raise ValueError("Expected one wheel metadata record")
        wheel_metadata = BytesParser().parsebytes(wheel.read(metadata_paths[0]))
        prefix = f"traceratops-{version}/"
        sdist_metadata = BytesParser().parsebytes(
            sdist.extractfile(prefix + "PKG-INFO").read()
        )
        for metadata in (wheel_metadata, sdist_metadata):
            if (
                metadata["Name"] != "traceratops"
                or metadata["Version"] != version
                or metadata["Requires-Python"] != project["requires-python"]
            ):
                raise ValueError("Distribution metadata differs from source metadata")
        for file in (source / "traceratops").rglob("*.py"):
            name = file.relative_to(source).as_posix()
            expected = file.read_bytes()
            if (
                wheel.read(name) != expected
                or sdist.extractfile(prefix + name).read() != expected
            ):
                raise ValueError(f"Distribution source mismatch: {name}")
        for name in ("pyproject.toml", "README.md", "LICENSE"):
            if sdist.extractfile(prefix + name).read() != (source / name).read_bytes():
                raise ValueError(f"Source distribution mismatch: {name}")
    return {"version": version}


def validate_installed(version):
    import traceratops
    from traceratops._version import __version__

    distribution = importlib.metadata.distribution("traceratops")
    if distribution.version != version or __version__ != version:
        raise ValueError("Installed version does not match built source")
    if (
        Path(traceratops.__file__).resolve().parent
        != Path(distribution.locate_file("traceratops")).resolve()
    ):
        raise ValueError("Import resolved outside the installed distribution")
    direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
    if direct_url.get("dir_info", {}).get("editable"):
        raise ValueError("Expected a wheel installation, not an editable checkout")
    entries = [e for e in distribution.entry_points if e.group == "console_scripts"]
    if not {"trace_filter", "trace_splitter", "plot_him_matrix"} <= {
        entry.name for entry in entries
    }:
        raise ValueError("Representative console entry points are missing")
    for entry in entries:
        if not callable(entry.load()):
            raise ValueError(f"Invalid console entry point: {entry.name}")
        subprocess.run(
            [str(Path(sysconfig.get_path("scripts")) / entry.name), "--help"],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    print(f"Installed traceratops {version}; {len(entries)} console commands passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    tag = commands.add_parser("tag")
    tag.add_argument("tag")
    tag.add_argument("--github-output", type=Path)
    dist = commands.add_parser("distributions")
    dist.add_argument("--source", type=Path, required=True)
    dist.add_argument("--dist", type=Path, required=True)
    dist.add_argument("--github-output", type=Path)
    installed = commands.add_parser("installed")
    installed.add_argument("version")
    args = parser.parse_args()
    if args.command == "tag":
        outputs = validate_tag(args.tag)
    elif args.command == "distributions":
        outputs = validate_distributions(args.source, args.dist)
    else:
        validate_installed(args.version)
        return
    print(json.dumps(outputs))
    if args.github_output:
        with args.github_output.open("a") as output:
            for key, value in outputs.items():
                output.write(f"{key}={value}\n")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
