import numpy as np
from astropy.table import Table

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.trace_physical_vs_genomic_distance import (
    calculate_physical_vs_genomic_distance,
    compute_all_pwd,
    compute_genomic_dist_map,
    extract_traces_numpy,
)


def _trace_table():
    table = ChromatinTraceTable()
    table.data = Table(
        {
            "Trace_ID": [1, 1, 1, 2, 2, 2],
            "Barcode #": [0, 1, 2, 0, 1, 2],
            "x": [0.0, 3.0, 6.0, 0.0, 4.0, 8.0],
            "y": [0.0, 4.0, 8.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "Chrom_Start": [0, 1000, 3000, 0, 1000, 3000],
            "Chrom_End": [1000, 2000, 4000, 1000, 2000, 4000],
        }
    )
    return table


def test_extract_traces_numpy_and_pairwise_distances():
    traces = extract_traces_numpy(_trace_table())
    pwd, axis_pwds = compute_all_pwd(traces, show_progress=False)

    assert traces.shape == (2, 3, 3)
    np.testing.assert_allclose(pwd[0], [5.0, 10.0, 5.0])
    np.testing.assert_allclose(axis_pwds["X"][1], [4.0, 8.0, 4.0])


def test_genomic_distances_are_reported_in_kbp():
    gen_dist = compute_genomic_dist_map(_trace_table())

    np.testing.assert_allclose(gen_dist, [1.0, 3.0, 2.0])


def test_calculate_physical_vs_genomic_distance_writes_csv_outputs(tmp_path):
    input_file = tmp_path / "traces.ecsv"
    output_file = tmp_path / "rs.csv"
    interloci_file = tmp_path / "interloci.csv"
    _trace_table().save(input_file)

    result = calculate_physical_vs_genomic_distance(
        input_file=input_file,
        output_csv=output_file,
        interloci_csv=interloci_file,
        gen_dist_bins=3,
        dist_threshold=1000,
        show_progress=False,
    )

    assert output_file.is_file()
    assert interloci_file.is_file()
    assert set(result["axis"]) == {"3D", "X", "Y", "Z"}
