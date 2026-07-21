from traceratops.plot_bootstrapping import (
    create_dict_args,
    get_adaptive_fontsize,
    get_adaptive_tick_labels,
)


class Args:
    input = "matrix.npy"
    uniqueBarcodes = "barcodes.csv"
    outputFolder = None
    fontsize = None
    axisLabel = False
    cMax = None
    cMin = None
    cMax_std = None
    output_format = "png"
    cMin_std = None
    shuffle = None
    cmap = None
    cmap_std = None
    N_bootstrap = None


def test_create_dict_args_converts_fontsize_to_numeric():
    args = Args()
    args.fontsize = "7.5"

    run_parameters = create_dict_args(args)

    assert run_parameters["fontsize"] == 7.5


def test_adaptive_fontsize_keeps_requested_size_for_small_matrices():
    assert get_adaptive_fontsize(matrix_size=10, max_fontsize=9) == 9


def test_adaptive_fontsize_scales_down_for_large_matrices():
    assert get_adaptive_fontsize(matrix_size=100, max_fontsize=9) == 5.0


def test_adaptive_fontsize_has_readable_lower_bound():
    assert get_adaptive_fontsize(matrix_size=1000, max_fontsize=9) == 4.0


def test_adaptive_tick_labels_keeps_all_labels_for_small_matrices():
    barcodes = ["1", "2", "3"]

    assert get_adaptive_tick_labels(barcodes) == barcodes


def test_adaptive_tick_labels_thins_labels_for_large_matrices():
    barcodes = [str(i) for i in range(1, 81)]

    labels = get_adaptive_tick_labels(barcodes)

    assert len(labels) == len(barcodes)
    assert labels[:8] == ["1", "", "", "", "5", "", "", ""]
    assert sum(bool(label) for label in labels) == 20
