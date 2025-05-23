import os
import subprocess

import numpy as np
import pytest

TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(TESTS_DIR, "data", "plot_him_matrix", "IN")
OUTPUT_DIR = os.path.join(TESTS_DIR, "data", "plot_him_matrix", "OUT")
INPUT_NPY = INPUT_DIR + os.sep + "n_cells_250_pwd_sc_matrix.npy"
INPUT_ECSV = INPUT_DIR + os.sep + "unique_barcodes.ecsv"
FAKE_MATRIX_PATH = ""
FAKE_BARCODES_PATH = ""


# Once before every tests of this file
@pytest.fixture(scope="module", autouse=True)
def setup_module():
    global FAKE_MATRIX_PATH, FAKE_BARCODES_PATH
    FAKE_MATRIX_PATH, FAKE_BARCODES_PATH = save_fake_data(INPUT_DIR)


def create_test_matrix():
    ab = [0.3, 0, 0, 0]  # Majority of NaN
    ac = [0.2, 0.2, 0.2, 0]  # Same for median and KDE
    ad = [0.1, 0.1, 0.9, 0.9]  # Diff for median and KDE
    bc = [0.6, 2.0, 2.0, 0.0]  # Majority of NaN after threshold
    bd = [0.0, 3.0, 0.0, 0.0]  # full of NaN after threshold
    cd = [0.1, 0.4, 0.4, 0.7]
    mtx = np.zeros((4, 4, 4))
    # Set diagonal to NaN
    mtx[0][0] = np.full(4, np.nan)
    mtx[1][1] = np.full(4, np.nan)
    mtx[2][2] = np.full(4, np.nan)
    mtx[3][3] = np.full(4, np.nan)
    # ab
    mtx[0][1] = np.array(ab)
    mtx[1][0] = np.array(ab)
    # ac
    mtx[0][2] = np.array(ac)
    mtx[2][0] = np.array(ac)
    # ad
    mtx[0][3] = np.array(ad)
    mtx[3][0] = np.array(ad)
    # bc
    mtx[1][2] = np.array(bc)
    mtx[2][1] = np.array(bc)
    # bd
    mtx[1][3] = np.array(bd)
    mtx[3][1] = np.array(bd)
    # cd
    mtx[2][3] = np.array(cd)
    mtx[3][2] = np.array(cd)

    mtx[mtx == 0] = np.nan
    return mtx


def save_fake_data(directory):
    mtx_name = "bc_4_cells_4.npy"
    mtx_path = os.path.join(directory, mtx_name)
    bc_name = "bc_4_unique.csv"
    bc_path = os.path.join(directory, bc_name)
    if not os.path.exists(mtx_path):
        mtx = create_test_matrix()
        np.save(mtx_path, mtx)
    if not os.path.exists(bc_path):
        barcodes = ["1", "3", "13", "999"]
        with open(bc_path, "w") as f:
            f.write("\n".join(barcodes))
    return mtx_path, bc_path


def compare_npy(gen_out_path, expected):
    gene = np.load(gen_out_path)
    expe = np.load(expected)
    return np.testing.assert_array_equal(gene, expe)


def check_script_run_normally(result, gen_out_path, expected):
    assert result.returncode == 0, f"Runtime error: {result.stderr}"
    assert os.path.exists(gen_out_path), f"Output file {gen_out_path} isn't created"
    # assert compare_npy(gen_out_path, expected), f"Difference detected between {gen_out_path} and {expected}"
    compare_npy(gen_out_path, expected)


def delete_paths(path_list: list) -> None:
    for p in path_list:
        if os.path.exists(p):
            os.remove(p)


def test_real_data():
    # Determine output file that should be generated by the script
    matrix_filename = "Fig_n_cells_250_pwd_sc_matrix_proximity_0.00-0.10.npy"
    png_filename = "Fig_n_cells_250_pwd_sc_matrix_proximity_0.00-0.10.png"
    nan_filename = "Fig_n_cells_250_pwd_sc_matrix_nan%_0.81-0.96.png"
    generated_png_path = os.path.join(INPUT_DIR, png_filename)
    generated_nan_path = os.path.join(INPUT_DIR, nan_filename)
    generated_matrix_path = os.path.join(INPUT_DIR, matrix_filename)
    expected_matrix_path = os.path.join(OUTPUT_DIR, matrix_filename)

    # Delete old files if exist to avoid conflict
    delete_paths([generated_png_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            INPUT_NPY,
            "-B",
            INPUT_ECSV,
            "-O",
            INPUT_DIR,
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )
    assert os.path.exists(generated_png_path)
    assert os.path.exists(generated_nan_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_png_path, generated_matrix_path, generated_nan_path])


missing_args_cases = [
    ["plot_him_matrix"],
    ["plot_him_matrix", "-M", INPUT_NPY],
    ["plot_him_matrix", "-B", INPUT_ECSV],
]


@pytest.mark.parametrize(
    "bad_args", missing_args_cases, ids=["BOTH", "--barcodes", "--matrix"]
)
def test_missing_arguments(bad_args):
    """Test the case when no arguments are provided to plot_him_matrix.py"""
    result = subprocess.run(
        bad_args,  # Run the script without arguments
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


formats = ["png", "svg", "pdf"]


@pytest.mark.parametrize("format", formats)
def test_plot_format(format):
    # Determine output file that should be generated by the script
    expected_nan_name = "Fig_bc_4_cells_4_nan%_0.25-0.75"
    generated_nan_path = os.path.join(INPUT_DIR, f"{expected_nan_name}.{format}")
    expected_name = "Fig_bc_4_cells_4_proximity_0.25-0.75"
    generated_matrix_path = os.path.join(INPUT_DIR, expected_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, expected_name + "." + format)
    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    _ = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
            "--plot_format",
            format,
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )

    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    assert os.path.exists(generated_nan_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


configs = ["median_0.20-3.00", "KDE_0.20-3.00", "proximity_0.25-0.75"]


@pytest.mark.parametrize("config", configs, ids=["median", "KDE", "proximity"])
def test_mode(config):
    # Determine output file that should be generated by the script
    mode_name = config.split("_")[0]
    expected_nan_name = "Fig_bc_4_cells_4_nan%_" + config.split("_")[1]
    generated_nan_path = os.path.join(INPUT_DIR, f"{expected_nan_name}.png")
    expected_name = "Fig_bc_4_cells_4_" + config
    generated_matrix_path = os.path.join(INPUT_DIR, expected_name + ".npy")
    expected_matrix_path = os.path.join(OUTPUT_DIR, expected_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, f"{expected_name}.png")

    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
            "--mode",
            mode_name,
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )
    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


def test_visualization_args():
    # Determine output file that should be generated by the script
    expected_name = "Fig_bc_4_cells_4_proximity_0.10-0.60"
    generated_matrix_path = os.path.join(INPUT_DIR, expected_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, expected_name + ".png")
    generated_nan_path = os.path.join(INPUT_DIR, "Fig_bc_4_cells_4_nan%_0.25-0.75.png")

    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    _ = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
            "--c_min",
            "0.1",
            "--c_max",
            "0.6",
            "--c_map",
            "jet",
            "--fontsize",
            "12",
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )

    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


thresholds = ["0.5", "0", "1", "-1"]


@pytest.mark.parametrize("threshold", thresholds)
def test_thresholds(threshold):
    # Determine output file that should be generated by the script
    expected_nan_name = "Fig_bc_4_cells_4_nan%"
    generated_nan_path = os.path.join(INPUT_DIR, f"{expected_nan_name}_0.25-0.75.png")
    expected_name = f"Fig_bc_4_cells_4_proximity_T{float(threshold):.1f}"
    if threshold == "1":
        expected_name += "_0.25-1.00"
    elif threshold == "0.5":
        expected_name += "_0.25-0.75"
    elif threshold == "0" or threshold == "-1":
        expected_name += "_nan-0.00"
    else:
        raise ValueError("Unexpected case")
    generated_matrix_path = os.path.join(INPUT_DIR, expected_name + ".npy")
    expected_matrix_path = os.path.join(OUTPUT_DIR, expected_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, f"{expected_name}.png")

    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
            "--threshold",
            threshold,
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )

    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    assert os.path.exists(generated_nan_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


def test_norm():
    # Determine output file that should be generated by the script
    expected_nan_name = "Fig_bc_4_cells_4_nan%_norm"
    generated_nan_path = os.path.join(INPUT_DIR, f"{expected_nan_name}_0.25-0.75.png")
    expected_name = "Fig_bc_4_cells_4_proximity_norm"
    expected_name += "_0.25-1.00"
    generated_matrix_path = os.path.join(INPUT_DIR, expected_name + ".npy")
    expected_matrix_path = os.path.join(OUTPUT_DIR, expected_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, f"{expected_name}.png")

    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
        ],
        capture_output=True,
        text=True,
    )

    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    assert os.path.exists(generated_nan_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


shuffle_cases = [
    ("1,3,13,999", "0.76"),
    ("1,3,13", "0.77"),
    ("3,13,999", "0.78"),
    ("3,999,1,13", "0.79"),
]
case_names = [
    "Change nothing",
    "Remove one barcode",
    "Remove first barcode",
    "Shuffle all",
]


@pytest.mark.parametrize(("shuffle_case", "c_max"), shuffle_cases, ids=case_names)
def test_shuffle(shuffle_case, c_max):
    # Determine output file that should be generated by the script
    expected_nan_name = "Fig_bc_4_cells_4_nan%_0.25-0.75"
    generated_nan_path = os.path.join(INPUT_DIR, f"{expected_nan_name}.png")
    expected_matrix_name = "Fig_bc_4_cells_4_proximity_0.25-" + c_max
    generated_matrix_path = os.path.join(INPUT_DIR, expected_matrix_name + ".npy")
    expected_matrix_path = os.path.join(OUTPUT_DIR, expected_matrix_name + ".npy")
    generated_plot_path = os.path.join(INPUT_DIR, f"{expected_matrix_name}.png")

    # Delete old files if exist to avoid conflict
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            FAKE_MATRIX_PATH,
            "-B",
            FAKE_BARCODES_PATH,
            "-O",
            INPUT_DIR,
            "--shuffle",
            shuffle_case,
            "--c_min",
            "0.25",
            "--c_max",
            c_max,
            "--keep_nan",
        ],
        capture_output=True,
        text=True,
    )

    assert os.path.exists(generated_plot_path)
    assert os.path.exists(generated_matrix_path)
    assert os.path.exists(generated_nan_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_plot_path, generated_matrix_path, generated_nan_path])


def test_matplotlib_v3_5_1_error():
    # Determine output file that should be generated by the script
    matrix_filename = "Fig_n_cells_250_pwd_sc_matrix_proximity_norm_T0.2_0.15-0.57.npy"
    png_filename = "Fig_n_cells_250_pwd_sc_matrix_proximity_norm_T0.2_0.15-0.57.png"
    nan_filename = "Fig_n_cells_250_pwd_sc_matrix_nan%_norm_0.81-0.96.png"
    generated_png_path = os.path.join(INPUT_DIR, png_filename)
    generated_nan_path = os.path.join(INPUT_DIR, nan_filename)
    generated_matrix_path = os.path.join(INPUT_DIR, matrix_filename)
    expected_matrix_path = os.path.join(OUTPUT_DIR, matrix_filename)

    # Delete old files if exist to avoid conflict
    delete_paths([generated_png_path, generated_matrix_path, generated_nan_path])

    # Run script with CLI
    result = subprocess.run(
        [
            "plot_him_matrix",
            "-M",
            INPUT_NPY,
            "-B",
            INPUT_ECSV,
            "-O",
            INPUT_DIR,
            "-T",
            "0.20",
            "--mode",
            "proximity",
            "--c_min",
            "0.15",
        ],
        capture_output=True,
        text=True,
    )
    assert os.path.exists(generated_png_path)
    assert os.path.exists(generated_nan_path)
    check_script_run_normally(result, generated_matrix_path, expected_matrix_path)
    delete_paths([generated_png_path, generated_matrix_path, generated_nan_path])
