from traceratops.plot_bootstrapping import create_dict_args, parse_arguments


def test_plot_bootstrapping_accepts_homogeneous_argument_names():
    parser = parse_arguments()
    args = parser.parse_args(
        [
            "--input",
            "matrix.npy",
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

    assert run_parameters["scPWDMatrix_filename"] == "matrix.npy"
    assert run_parameters["uniqueBarcodes"] == "barcodes.ecsv"
    assert run_parameters["outputFolder"] == "plots"
    assert run_parameters["cMin"] == 0.1
    assert run_parameters["cMax"] == 0.8
    assert run_parameters["cmap"] == "viridis"


def test_plot_bootstrapping_keeps_legacy_input_argument_names():
    parser = parse_arguments()
    args = parser.parse_args(
        [
            "--input",
            "matrix.npy",
            "--uniqueBarcodes",
            "barcodes.ecsv",
            "--outputFolder",
            "plots",
        ]
    )

    run_parameters = create_dict_args(args)

    assert run_parameters["scPWDMatrix_filename"] == "matrix.npy"
    assert run_parameters["uniqueBarcodes"] == "barcodes.ecsv"
    assert run_parameters["outputFolder"] == "plots"
