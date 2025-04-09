import filecmp
import os
import subprocess

import pytest

TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_filter", "IN")
OUTPUT_DIR = os.path.join(TESTS_DIR, "data", "trace_filter", "OUT")


def remove_file_if_exists(path):
    if os.path.exists(path):
        os.remove(path)


def run_trace_filter(args, shell=False):
    """
    If shell == True: Allows shell commands like `|`
    """
    return subprocess.run(args, capture_output=True, text=True, shell=shell)


def check_output(result, gen_out_path, filtered, expected):
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    assert os.path.exists(gen_out_path), f"Output file {filtered} isn't created"
    assert filecmp.cmp(
        gen_out_path, expected, shallow=False
    ), f"Difference detected between {gen_out_path} and {expected}"


def _test_trace_filter_common(
    cmd_args,
    generated_output_path,
    filtered_filename,
    expected_output_path,
    shell=False,
):
    """Helper to run trace_filter, check output, and clean generated files."""
    remove_file_if_exists(generated_output_path)
    result = run_trace_filter(cmd_args, shell=shell)
    check_output(result, generated_output_path, filtered_filename, expected_output_path)
    remove_file_if_exists(generated_output_path)


# ==== INPUT FILES SETUP ====

input_files = [
    f
    for f in os.listdir(INPUT_DIR)
    if f.endswith(".ecsv") and "_filtered" not in f and "trace" in f
]


@pytest.mark.parametrize("input_file", input_files)
def test_trace_filter_input(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    filtered_filename = input_file.replace(".ecsv", "_filtered.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    _test_trace_filter_common(
        ["trace_filter", "--input", input_path],
        generated_output_path,
        filtered_filename,
        expected_output_path,
    )


def test_missing_arguments():
    """Test the case when no arguments are provided to trace_filter.py."""
    result = run_trace_filter(["trace_filter"])  # Run the script without arguments

    output = result.stdout
    assert (
        "Error: No argument provided" in output
    ), "Expected error message not found in output"
    assert (
        "Redirecting to `--help` option:" in output
    ), "Expected redirection message not found in output"
    assert "usage:" in output, "Expected argparse usage information not found"
    assert (
        result.returncode == 0
    ), "Script should exit normally (code 0) when missing arguments"


# ==== PIPE FILES SETUP ====

forpipe_files = [
    f for f in os.listdir(INPUT_DIR) if f.endswith(".txt") and "forpipe" in f
]


@pytest.mark.parametrize("input_file", forpipe_files)
def test_trace_filter_pipe(input_file):
    result = run_trace_filter(
        f"cd {INPUT_DIR} && cat {input_file} | trace_filter --pipe",
        shell=True,
    )

    print("\n===== STDOUT =====\n", result.stdout)
    print("===== STDERR =====\n", result.stderr)

    assert result.returncode == 0, f"Runtime error: {result.stderr}"

    # Clean up generated files
    remove_file_if_exists(os.path.join(INPUT_DIR, "one_trace_four_spots_filtered.ecsv"))
    remove_file_if_exists(
        os.path.join(INPUT_DIR, "two_traces_seven_spots_filtered.ecsv")
    )


input_files = [f for f in os.listdir(INPUT_DIR) if "one_trace_four_spots.ecsv" in f]


@pytest.mark.parametrize("input_file", input_files)
def test_trace_filter_output(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    output_arg = "custom_out"
    filtered_filename = input_file.replace(".ecsv", f"_{output_arg}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    _test_trace_filter_common(
        ["trace_filter", "--input", input_path, "--output", output_arg],
        generated_output_path,
        filtered_filename,
        expected_output_path,
    )


input_files = [f for f in os.listdir(INPUT_DIR) if "duplicate_spot" in f]


@pytest.mark.parametrize("input_file", input_files)
def test_clean_spots(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    output_arg = "cleaned"
    filtered_filename = input_file.replace(".ecsv", f"_{output_arg}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    _test_trace_filter_common(
        [
            "trace_filter",
            "--input",
            input_path,
            "--output",
            output_arg,
            "--clean_spots",
        ],
        generated_output_path,
        filtered_filename,
        expected_output_path,
    )

    base_path = os.path.splitext(input_path)[0]
    remove_file_if_exists(base_path + "_before_filtering.png")
    remove_file_if_exists(base_path + "_filtered.png")


@pytest.mark.parametrize("command", ["remove_label", "keep_label"])
def test_label(command):
    input_file = "two_traces_seven_spots.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    end = "not-label_1" if command == "remove_label" else "label_1"
    filtered_filename = input_file.replace(".ecsv", f"_filtered_{end}.ecsv")
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    _test_trace_filter_common(
        ["trace_filter", "--input", input_path, f"--{command}", "label_1"],
        generated_output_path,
        filtered_filename,
        expected_output_path,
    )


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

    _test_trace_filter_common(
        [
            "trace_filter",
            "--input",
            input_path,
            "--remove_barcode",
            bc_list,
            "--output",
            output,
        ],
        generated_output_path,
        filtered_filename,
        expected_output_path,
    )
