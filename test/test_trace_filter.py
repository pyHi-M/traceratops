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

# Keep parametrized inputs explicit so generated files from interrupted or failed
# test runs do not become test cases on the next invocation. Several tests write
# outputs back into INPUT_DIR before cleaning them up; if a run is interrupted,
# those stale ``.ecsv``/``.png`` files must not pollute collection.
INPUT_FILES = sorted(os.listdir(INPUT_DIR))

trace_input_files = [
    "one_trace_four_spots.ecsv",
    "trace_3D_barcode_KDtree_ROI-5.ecsv",
    "two_traces_seven_spots.ecsv",
]
forpipe_files = ["forpipe1file.txt", "forpipe2files.txt"]
one_trace_files = ["one_trace_four_spots.ecsv"]
duplicate_spot_files = ["duplicate_spot.ecsv", "duplicate_spot_id.ecsv"]


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


def test_localization_intensity_column_prefers_mean_intensity():
    from astropy.table import Table
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    localization_table = Table(
        rows=[(10.0, 20.0)],
        names=("mean_intensity", "peak"),
    )

    assert (
        ChromatinTraceTable._get_localization_intensity_column(localization_table)
        == "mean_intensity"
    )


def test_localization_intensity_column_accepts_legacy_peak():
    from astropy.table import Table
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    localization_table = Table(rows=[(20.0,)], names=("peak",))

    assert (
        ChromatinTraceTable._get_localization_intensity_column(localization_table)
        == "peak"
    )


def test_filter_by_localization_metrics_combines_thresholds():
    from astropy.table import Table
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    trace = ChromatinTraceTable()
    trace.data = Table(
        rows=[
            ("spot-1", "trace-a", 1),
            ("spot-2", "trace-a", 2),
            ("spot-3", "trace-a", 3),
            ("spot-4", "trace-a", 4),
        ],
        names=("Spot_ID", "Trace_ID", "Barcode #"),
    )
    localizations = Table(
        rows=[
            ("spot-1", 10.0, 5.0, 1, 0.9),
            ("spot-2", 9.0, 4.0, 1, 0.8),
            ("spot-3", 11.0, 6.0, 0, 0.7),
            ("spot-4", 12.0, 7.0, 1, 0.6),
        ],
        names=("Buid", "mean_intensity", "snr", "object_class", "roundness"),
    )

    intensities_kept = trace.filter_by_localization_metrics(
        trace,
        localizations,
        {
            "intensity": 10.0,
            "snr": 5.0,
            "object_class": 1,
            "roundness": 0.6,
        },
    )

    assert list(trace.data["Spot_ID"]) == ["spot-1", "spot-4"]
    assert intensities_kept == [10.0, 12.0]


def test_clean_spots_preserves_reused_spot_ids_across_traces():
    from astropy.table import Table
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    trace = ChromatinTraceTable()
    trace.data = Table(
        rows=[
            ("1", "trace-a", 1, 0.0, 0.0, 0.0),
            ("1", "trace-b", 2, 1.0, 1.0, 1.0),
            ("2", "trace-b", 3, 2.0, 2.0, 2.0),
            ("2", "trace-b", 4, 3.0, 3.0, 3.0),
        ],
        names=("Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z"),
    )

    trace.remove_duplicates()

    assert len(trace.data) == 2
    assert list(trace.data["Trace_ID"]) == ["trace-a", "trace-b"]
    assert list(trace.data["Spot_ID"]) == ["1", "1"]


def test_filter_traces_by_n_handles_empty_trace_table():
    from astropy.table import Table
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable

    trace = ChromatinTraceTable()
    trace.data = Table(
        names=("Spot_ID", "Trace_ID", "Barcode #"), dtype=(str, str, int)
    )

    trace.filter_traces_by_n(minimum_number_barcodes=4)

    assert len(trace.data) == 0
