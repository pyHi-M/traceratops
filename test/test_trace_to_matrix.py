import numpy as np
from astropy.table import Table

from traceratops.core.build_matrix import BuildMatrix
from traceratops.trace_to_matrix import create_dict_args, parse_arguments


def test_trace_to_matrix_output_prefix_uses_output_folder_and_input_stem(tmp_path):
    input_trace = tmp_path / "input" / "Trace_3D_barcode_mask-mask0_ROI-51.ecsv"
    output_dir = tmp_path / "output"
    input_trace.parent.mkdir()
    input_trace.touch()

    matrix_builder = BuildMatrix({}, {})

    assert matrix_builder.get_output_prefix(input_trace, output_dir) == str(
        output_dir / "Trace_3D_barcode_mask-mask0_ROI-51_Matrix"
    )
    assert output_dir.is_dir()


def test_trace_to_matrix_output_prefix_defaults_to_input_folder(tmp_path):
    input_trace = tmp_path / "input" / "Trace_3D_barcode_mask-mask0_ROI-51.ecsv"
    input_trace.parent.mkdir()
    input_trace.touch()

    matrix_builder = BuildMatrix({}, {})

    assert matrix_builder.get_output_prefix(input_trace) == str(
        input_trace.parent / "Trace_3D_barcode_mask-mask0_ROI-51_Matrix"
    )


class _TraceTable:
    def __init__(self, data):
        self.data = data


def test_build_distance_matrix_parallel_matches_serial_for_duplicate_barcodes():
    data = Table(
        {
            "Trace_ID": [1, 1, 1, 1, 2, 2, 2],
            "Barcode #": [10, 20, 20, 30, 10, 20, 30],
            "x": [0.0, 1.0, 2.0, 4.0, 0.0, 3.0, 6.0],
            "y": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    )
    matrix_builder = BuildMatrix({}, {})
    matrix_builder.trace_table = _TraceTable(data)

    matrix_builder.build_distance_matrix(n_jobs=1)
    serial_matrix = matrix_builder.sc_matrix.copy()

    matrix_builder.build_distance_matrix(n_jobs=2)

    np.testing.assert_equal(matrix_builder.sc_matrix, serial_matrix)


def test_trace_to_matrix_optional_arguments_default_and_parse():
    parser = parse_arguments()

    default_args = create_dict_args(parser.parse_args([]))
    parallel_args = create_dict_args(
        parser.parse_args(["--n_jobs", "-1", "--plot_histograms"])
    )

    assert default_args["n_jobs"] == 1
    assert default_args["plot_histograms"] is False
    assert parallel_args["n_jobs"] == -1
    assert parallel_args["plot_histograms"] is True
