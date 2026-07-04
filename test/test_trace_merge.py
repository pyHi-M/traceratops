import filecmp
import os
import subprocess
import uuid

TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_merge", "IN")
OUTPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_merge", "OUT")


def test_spots_2():
    out_file = "merged_trace_1_2.ecsv"
    # Run script with CLI
    result = subprocess.run(
        f"cd {INPUT_DIR} && ls trace_*_spots_2.ecsv | trace_merge --name {out_file}",
        capture_output=True,
        text=True,
        shell=True,  # Allows shell commands like `|`
    )
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    gen_file = os.path.join(INPUT_DIR, out_file)
    expected_file = os.path.join(OUTPUT_DIR, out_file)
    assert os.path.exists(gen_file)
    assert os.path.exists(expected_file)
    assert filecmp.cmp(
        gen_file, expected_file, shallow=False
    ), f"Difference detected between {gen_file} and {expected_file}"
    os.remove(gen_file)


def test_merge_conflict():
    out_file = "merged_trace_1_2_3.ecsv"
    # Run script with CLI
    result = subprocess.run(
        f"cd {INPUT_DIR} && ls trace_*.ecsv | trace_merge --name {out_file}",
        capture_output=True,
        text=True,
        shell=True,  # Allows shell commands like `|`
    )
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    gen_file = os.path.join(INPUT_DIR, out_file)
    expected_file = os.path.join(OUTPUT_DIR, out_file)
    assert os.path.exists(gen_file)
    assert os.path.exists(expected_file)
    assert filecmp.cmp(
        gen_file, expected_file, shallow=False
    ), f"Difference detected between {gen_file} and {expected_file}"
    os.remove(gen_file)


def test_merge_4dn_numeric_spot_id_with_ecsv_spot_id(tmp_path):
    from astropy.table import Table, vstack
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    fofct_file = tmp_path / "numeric_spot_ids.4dn"
    fofct_file.write_text(
        "##FOF-CT_version=v0.1\n"
        "##Table_namespace=4dn_FOF-CT_core\n"
        "##genome_assembly=GRCm38\n"
        "##XYZ_unit=nm\n"
        "##columns=(Spot_ID, Trace_ID, X, Y, Z, Chrom, Chrom_Start, Chrom_End)\n"
        "1000000,500365,1275.7,1817.9,5362.4,chr13,55945001,55955000\n"
    )

    trace = ChromatinTraceTable()
    fofct_table = trace.load(str(fofct_file))
    ecsv_table = Table(
        rows=[("existing", "trace-a", 1.0, 2.0, 3.0, "chr13", 1, 2, 0, -1, 1, "None")],
        names=fofct_table.colnames,
    )

    merged = vstack([ecsv_table, fofct_table])

    assert uuid.UUID(str(fofct_table["Spot_ID"][0]))
    assert uuid.UUID(str(fofct_table["Trace_ID"][0]))
    assert len(merged) == 2


def test_4dn_conversion_relabels_ids_to_pyhim_nomenclature(tmp_path, monkeypatch):
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    generated_ids = iter(
        [
            uuid.UUID("00000000-0000-0000-0000-000000000001"),
            uuid.UUID("00000000-0000-0000-0000-000000000002"),
            uuid.UUID("00000000-0000-0000-0000-000000000003"),
            uuid.UUID("00000000-0000-0000-0000-000000000004"),
            uuid.UUID("00000000-0000-0000-0000-000000000005"),
        ]
    )
    monkeypatch.setattr(
        "traceratops.core.chromatin_trace_table.uuid.uuid4",
        lambda: next(generated_ids),
    )

    fofct_file = tmp_path / "numeric_trace_ids.4dn"
    fofct_file.write_text(
        "##FOF-CT_version=v0.1\n"
        "##Table_namespace=4dn_FOF-CT_core\n"
        "##genome_assembly=GRCm38\n"
        "##XYZ_unit=nm\n"
        "##columns=(Spot_ID, Trace_ID, X, Y, Z, Chrom, Chrom_Start, Chrom_End)\n"
        "1000000,500365,1275.7,1817.9,5362.4,chr13,55945001,55955000\n"
        "1000001,500365,1276.7,1818.9,5363.4,chr13,55955001,55965000\n"
        "1000002,500366,1277.7,1819.9,5364.4,chr13,55965001,55975000\n"
    )

    trace = ChromatinTraceTable()
    fofct_table = trace.load(str(fofct_file))

    assert list(fofct_table["Spot_ID"]) == [
        "00000000-0000-0000-0000-000000000001",
        "00000000-0000-0000-0000-000000000002",
        "00000000-0000-0000-0000-000000000003",
    ]
    assert list(fofct_table["Trace_ID"]) == [
        "00000000-0000-0000-0000-000000000004",
        "00000000-0000-0000-0000-000000000004",
        "00000000-0000-0000-0000-000000000005",
    ]
