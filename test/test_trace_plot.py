import os

from traceratops.trace_plot import parse_arguments, create_dict_args, runtime

TESTS_DIR = os.path.dirname(os.path.realpath(__file__))
TRACE_FILE = os.path.join(
    TESTS_DIR, "data", "trace_filter", "IN", "two_traces_seven_spots.ecsv"
)


def test_number_traces_exports_first_n_traces(tmp_path):
    runtime(
        trace_files=[TRACE_FILE],
        folder_path=str(tmp_path),
        select_traces="first",
        number_traces=1,
    )

    generated = sorted(path.name for path in tmp_path.glob("*.pdb"))
    assert generated == ["aaaaaaaa.pdb"]


def test_number_traces_argument_selects_first_traces_mode():
    parser = parse_arguments()
    args = parser.parse_args(["--input", TRACE_FILE, "--number_traces", "1"])

    parsed = create_dict_args(args)

    assert parsed["number_traces"] == 1
    assert parsed["select_traces"] == "first"


def test_all_overrides_number_traces(tmp_path):
    runtime(
        trace_files=[TRACE_FILE],
        folder_path=str(tmp_path),
        select_traces="all",
        number_traces=1,
    )

    generated = sorted(path.name for path in tmp_path.glob("*.pdb"))
    assert generated == ["aaaaaaaa.pdb", "bbbbbbbb.pdb"]
