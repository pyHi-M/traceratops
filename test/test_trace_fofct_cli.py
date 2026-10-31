"""Exercise the installed FOF-CT conversion commands and their metadata."""

import json
import os
import subprocess
import sysconfig
from email.utils import getaddresses
from importlib.metadata import metadata as package_metadata
from pathlib import Path

import numpy as np
from astropy.table import Table


def _command(name):
    suffix = ".exe" if os.name == "nt" else ""
    return str(Path(sysconfig.get_path("scripts")) / f"{name}{suffix}")


def test_installed_fofct_export_import_roundtrip(tmp_path):
    source = tmp_path / "trace.ecsv"
    bed = tmp_path / "barcodes.bed"
    metadata = tmp_path / "parameters.json"
    exported = tmp_path / "trace.csv"
    restored = tmp_path / "restored.ecsv"
    table = Table(
        {
            "Spot_ID": ["spot-a", "spot-b"],
            "Trace_ID": ["trace-a", "trace-a"],
            "x": [1.25, 2.5],
            "y": [3.0, 4.0],
            "z": [0.1, 0.2],
            "Chrom": ["chr1", "chr1"],
            "Chrom_Start": [100, 200],
            "Chrom_End": [150, 250],
            "ROI #": [1, 1],
            "Barcode #": [1, 2],
            "Mask_id": [-1, -1],
            "label": ["None", "None"],
        }
    )
    table.write(source, format="ascii.ecsv")
    bed.write_text("chr1\t100\t150\t1\nchr1\t200\t250\t2\n")
    metadata.write_text(
        json.dumps(
            {
                "genome_assembly": "test",
                "experimenter_name": "Test User",
                "experimenter_contact": "test@example.org",
            }
        )
    )
    result = subprocess.run(
        [
            _command("trace_export_to_fofct"),
            "--input",
            str(source),
            "--bed_file",
            str(bed),
            "--json_file",
            str(metadata),
            "--output_file",
            str(exported),
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert "##XYZ_unit=micron" in exported.read_text()
    installed_authors = getaddresses(
        package_metadata("traceratops").get_all("Author-email", [])
    )
    expected_authors = []
    for name, _ in installed_authors:
        parts = name.split()
        if not parts:
            continue
        initials = "-".join(part[0] for part in parts[0].split("-"))
        expected_authors.append(
            f"{parts[-1]}, {initials};" if len(parts) > 1 else f"{name};"
        )
    author_header = next(
        line
        for line in exported.read_text().splitlines()
        if line.startswith("#Software_Authors:")
    )
    assert author_header == "#Software_Authors: " + " ".join(expected_authors)
    result = subprocess.run(
        [
            _command("trace_import_from_fofct"),
            "--input",
            str(exported),
            "--bed_file",
            str(bed),
            "--output_file",
            str(restored),
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    output = Table.read(restored, format="ascii.ecsv")
    for name in ("x", "y", "z", "Chrom_Start", "Chrom_End", "Barcode #", "ROI #"):
        np.testing.assert_array_equal(output[name], table[name])
    assert list(output["Spot_ID"]) == list(table["Spot_ID"])
    assert list(output["Trace_ID"]) == list(table["Trace_ID"])
