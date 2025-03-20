import filecmp
import os
import subprocess

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
