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
    """If shell == True: Allows shell commands like `|`"""
    return subprocess.run(args, capture_output=True, text=True, shell=shell)


def check_output(result, gen_out_path, filtered, expected):
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    assert os.path.exists(gen_out_path), f"Output file {filtered} isn't created"
    assert filecmp.cmp(
        gen_out_path, expected, shallow=False
    ), f"Difference detected between {gen_out_path} and {expected}"


def _test_trace_filter_common(
    input_file, args, suffix="_filtered", clean_png=False, shell=False
):
    """Helper to run trace_filter, check output and clean generated files."""
    base_name = os.path.splitext(input_file)[0]
    filtered_filename = f"{base_name}{suffix}.ecsv"
    generated_output_path = os.path.join(INPUT_DIR, filtered_filename)
    expected_output_path = os.path.join(OUTPUT_DIR, filtered_filename)

    remove_file_if_exists(generated_output_path)
    result = run_trace_filter(args, shell=shell)
    check_output(result, generated_output_path, filtered_filename, expected_output_path)
    remove_file_if_exists(generated_output_path)

    if clean_png:
        all_suffix = ["before_filtering", "filtered", "filtered_intensities"]
        for suf in all_suffix:
            remove_file_if_exists(os.path.join(INPUT_DIR, f"{base_name}_{suf}.png"))
        # just test_intensity():
        remove_file_if_exists(
            os.path.join(INPUT_DIR, "intensity_localization_intensities.png")
        )


# ==== FILE LISTS ====

INPUT_FILES = os.listdir(INPUT_DIR)

trace_input_files = [
    f
    for f in INPUT_FILES
    if f.endswith(".ecsv") and "_filtered" not in f and "trace" in f
]
forpipe_files = [f for f in INPUT_FILES if f.endswith(".txt") and "forpipe" in f]
one_trace_files = [f for f in INPUT_FILES if "one_trace_four_spots.ecsv" in f]
duplicate_spot_files = [f for f in INPUT_FILES if "duplicate_spot" in f]


# ==== TESTS ====


@pytest.mark.parametrize("input_file", trace_input_files)
def test_trace_filter_input(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    args = ["trace_filter", "--input", input_path]
    _test_trace_filter_common(input_file, args)


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


@pytest.mark.parametrize("input_file", forpipe_files)
def test_trace_filter_pipe(input_file):
    result = run_trace_filter(
        f"cd {INPUT_DIR} && cat {input_file} | trace_filter --pipe",
        shell=True,
    )

    print("\n===== STDOUT =====\n", result.stdout)
    print("===== STDERR =====\n", result.stderr)

    assert result.returncode == 0, f"Runtime error: {result.stderr}"

    remove_file_if_exists(os.path.join(INPUT_DIR, "one_trace_four_spots_filtered.ecsv"))
    remove_file_if_exists(
        os.path.join(INPUT_DIR, "two_traces_seven_spots_filtered.ecsv")
    )


@pytest.mark.parametrize("input_file", one_trace_files)
def test_trace_filter_output(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    args = ["trace_filter", "--input", input_path, "--output", "custom_out"]
    _test_trace_filter_common(input_file, args, suffix="_custom_out")


@pytest.mark.parametrize("input_file", duplicate_spot_files)
def test_clean_spots(input_file):
    input_path = os.path.join(INPUT_DIR, input_file)
    args = [
        "trace_filter",
        "--input",
        input_path,
        "--output",
        "cleaned",
        "--clean_spots",
    ]
    _test_trace_filter_common(input_file, args, suffix="_cleaned", clean_png=True)


@pytest.mark.parametrize("command", ["remove_label", "keep_label"])
def test_label(command):
    input_file = "two_traces_seven_spots.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    args = ["trace_filter", "--input", input_path, f"--{command}", "label_1"]
    suffix = (
        "_filtered_not-label_1" if command == "remove_label" else "_filtered_label_1"
    )
    _test_trace_filter_common(input_file, args, suffix=suffix)


@pytest.mark.parametrize(
    ("bc_list", "output"),
    [("999", "not_present"), ("1", "one"), ("3,1,5", "list")],
    ids=["not present", "one", "list"],
)
def test_remove_barcode(bc_list, output):
    input_file = "remove_barcode.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    args = [
        "trace_filter",
        "--input",
        input_path,
        "--remove_barcode",
        bc_list,
        "--output",
        output,
    ]
    _test_trace_filter_common(input_file, args, suffix=f"_{output}")


@pytest.mark.parametrize("min", ["4", "5"])
def test_n_barcodes(min):
    input_file = "two_traces_seven_spots.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    args = ["trace_filter", "--input", input_path, "--n_barcodes", min, "--output", min]
    _test_trace_filter_common(input_file, args, suffix=f"_{min}")


def test_xyz_min_max():
    input_file = "xyz_min_max.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)

    args = [
        "trace_filter",
        "--input",
        input_path,
        "--x_min",
        "80",
        "--x_max",
        "200",
        "--y_min",
        "30",
        "--y_max",
        "70",
        "--z_min",
        "5",
        "--z_max",
        "8",
    ]

    _test_trace_filter_common(input_file, args)


def test_intensity():
    input_file = "duplicate_spot.ecsv"
    input_path = os.path.join(INPUT_DIR, input_file)
    intensity_path = os.path.join(INPUT_DIR, "intensity.ecsv")

    args = [
        "trace_filter",
        "--input",
        input_path,
        "--clean_spots",
        "--localization_file",
        intensity_path,
        "--intensity_min",
        "555",
        "--output",
        "intensity",
    ]

    _test_trace_filter_common(input_file, args, suffix="_intensity", clean_png=True)
