from traceratops.core.build_matrix import BuildMatrix


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
