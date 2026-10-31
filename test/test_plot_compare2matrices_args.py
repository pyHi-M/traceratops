from traceratops.plot_compare2matrices import create_dict_args, parse_arguments


def test_plot_compare2matrices_accepts_homogeneous_argument_names():
    parser = parse_arguments()
    args = parser.parse_args(
        [
            "--input1",
            "matrix_1.npy",
            "--input2",
            "matrix_2.npy",
            "--barcodes",
            "barcodes.ecsv",
            "--output",
            "plots",
            "--c_min",
            "0.1",
            "--c_max",
            "0.8",
            "--cmap",
            "viridis",
        ]
    )

    run_parameters = create_dict_args(args)

    assert run_parameters["input1"] == "matrix_1.npy"
    assert run_parameters["input2"] == "matrix_2.npy"
    assert run_parameters["uniqueBarcodes"] == "barcodes.ecsv"
    assert run_parameters["outputFolder"] == "plots"
    assert run_parameters["cMin"] == 0.1
    assert run_parameters["cMax"] == 0.8
    assert run_parameters["cmap"] == "viridis"


def test_plot_compare2matrices_keeps_legacy_argument_names():
    parser = parse_arguments()
    args = parser.parse_args(
        [
            "--input1",
            "matrix_1.npy",
            "--input2",
            "matrix_2.npy",
            "--uniqueBarcodes",
            "barcodes.ecsv",
            "--outputFolder",
            "plots",
            "--cMin",
            "0.1",
            "--cMax",
            "0.8",
            "--c_map",
            "viridis",
        ]
    )

    run_parameters = create_dict_args(args)

    assert run_parameters["uniqueBarcodes"] == "barcodes.ecsv"
    assert run_parameters["outputFolder"] == "plots"
    assert run_parameters["cMin"] == 0.1
    assert run_parameters["cMax"] == 0.8
    assert run_parameters["cmap"] == "viridis"
