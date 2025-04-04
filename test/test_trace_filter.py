import filecmp
import os
import subprocess

import pytest

TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_filter", "IN")
OUTPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_filter", "OUT")


def check_script_run_normally(result, gen_out_path, filtered, expected):
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    assert os.path.exists(gen_out_path), f"Output file {filtered} isn't created"
    assert filecmp.cmp(
        gen_out_path, expected, shallow=False
    ), f"Difference detected between {gen_out_path} and {expected}"


input_files = [
    f
    for f in os.listdir(INPUT_DIR)
    if f.endswith(".ecsv") and "_filtered" not in f and "trace" in f
]


@pytest.mark.parametrize("input_file", input_files)
def test_trace_filter_input(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)

    # Determine output file that should be generated by the script
    filtered_filename = input_file.replace(".ecsv", "_filtered.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    # Delete old filtered file if exist to avoid conflict
    if os.path.exists(generated_output_path):
        os.remove(generated_output_path)

    # Run script with CLI
    result = subprocess.run(
        ["trace_filter", "--input", input_path],
        capture_output=True,
        text=True,
    )
    check_script_run_normally(
        result, generated_output_path, filtered_filename, expected_output_path
    )
    os.remove(generated_output_path)


def test_missing_arguments():
    """Test the case when no arguments are provided to trace_filter.py."""
    result = subprocess.run(
        ["trace_filter"],  # Run the script without arguments
        capture_output=True,
        text=True,
    )

    # Capture standard output
    output = result.stdout

    # Assertions
    assert (
        "Error: No argument provided" in output
    ), "Expected error message not found in output"
    assert (
        "Redirecting to `--help` option:" in output
    ), "Expected redirection message not found in output"
    assert "usage:" in output, "Expected argparse usage information not found"

    # Ensure the script exits successfully (exit code 0)
    assert (
        result.returncode == 0
    ), "Script should exit normally (code 0) when missing arguments"


forpipe_files = [
    f for f in os.listdir(INPUT_DIR) if f.endswith(".txt") and "forpipe" in f
]


@pytest.mark.parametrize("input_file", forpipe_files)
def test_trace_filter_pipe(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    # Run script with CLI
    result = subprocess.run(
        f"cd {INPUT_DIR} && cat {input_path} | trace_filter --pipe",
        capture_output=True,
        text=True,
        shell=True,  # Allows shell commands like `|`
    )
    # Print the output explicitly so pytest displays it
    print("\n===== STDOUT =====")
    print(result.stdout)
    print("===== STDERR =====")
    print(result.stderr)
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    file1 = os.path.join(INPUT_DIR, "one_trace_four_spots_filtered.ecsv")
    file2 = os.path.join(INPUT_DIR, "two_traces_seven_spots_filtered.ecsv")
    if os.path.exists(file1):
        os.remove(file1)
    if os.path.exists(file2):
        os.remove(file2)


input_files = [f for f in os.listdir(INPUT_DIR) if "one_trace_four_spots.ecsv" in f]


@pytest.mark.parametrize("input_file", input_files)
def test_trace_filter_output(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    output_arg = "custom_out"
    # Determine output file that should be generated by the script
    filtered_filename = input_file.replace(".ecsv", f"_{output_arg}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    # Delete old filtered file if exist to avoid conflict
    if os.path.exists(generated_output_path):
        os.remove(generated_output_path)

    # Run script with CLI
    result = subprocess.run(
        ["trace_filter", "--input", input_path, "--output", output_arg],
        capture_output=True,
        text=True,
    )

    check_script_run_normally(
        result, generated_output_path, filtered_filename, expected_output_path
    )
    os.remove(generated_output_path)


input_files = [f for f in os.listdir(INPUT_DIR) if "duplicate_spot" in f]


@pytest.mark.parametrize("input_file", input_files)
def test_clean_spots(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    output_arg = "cleaned"
    # Determine output file that should be generated by the script
    filtered_filename = input_file.replace(".ecsv", f"_{output_arg}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)
    # Delete old filtered file if exist to avoid conflict
    if os.path.exists(generated_output_path):
        os.remove(generated_output_path)
    # Run script with CLI
    result = subprocess.run(
        [
            "trace_filter",
            "--input",
            input_path,
            "--output",
            output_arg,
            "--clean_spots",
        ],
        capture_output=True,
        text=True,
    )
    check_script_run_normally(
        result, generated_output_path, filtered_filename, expected_output_path
    )
    os.remove(generated_output_path)
    in_path_base = input_path.split(".")[0]
    os.remove(in_path_base + "_before_filtering.png")
    os.remove(in_path_base + "_filtered.png")


@pytest.mark.parametrize("command", ["remove_label", "keep_label"])
def test_label(command):
    input_file = "two_traces_seven_spots.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    # Determine output file that should be generated by the script
    end = "not-label_1" if command == "remove_label" else "label_1"
    filtered_filename = input_file.replace(".ecsv", f"_filtered_{end}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)
    # Delete old filtered file if exist to avoid conflict
    if os.path.exists(generated_output_path):
        os.remove(generated_output_path)
    # Run script with CLI
    result = subprocess.run(
        [
            "trace_filter",
            "--input",
            input_path,
            f"--{command}",
            "label_1",
        ],
        capture_output=True,
        text=True,
    )
    check_script_run_normally(
        result, generated_output_path, filtered_filename, expected_output_path
    )
    os.remove(generated_output_path)


@pytest.mark.parametrize(
    ("bc_list", "output"),
    [("999", "not_present"), ("1", "one"), ("3,1,5", "list")],
    ids=["not present", "one", "list"],
)
def test_remove_barcode(bc_list, output):
    input_file = "remove_barcode.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    filtered_filename = input_file.replace(".ecsv", f"_{output}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)
    # Delete old filtered file if exist to avoid conflict
    if os.path.exists(generated_output_path):
        os.remove(generated_output_path)
    # Run script with CLI
    result = subprocess.run(
        [
            "trace_filter",
            "--input",
            input_path,
            "--remove_barcode",
            bc_list,
            "--output",
            output,
        ],
        capture_output=True,
        text=True,
    )
    check_script_run_normally(
        result, generated_output_path, filtered_filename, expected_output_path
    )
    os.remove(generated_output_path)
