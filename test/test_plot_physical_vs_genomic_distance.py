import numpy as np
from astropy.table import Table

from traceratops.plot_physical_vs_genomic_distance import (
    collect_distances,
    runtime,
    summarize_pairs,
)


def make_trace_table():
    return Table(
        rows=[
            ("s1", "t1", 0.0, 0.0, 0.0, "chr1", 0, 1000, 1, 1, 1, "label"),
            ("s2", "t1", 3.0, 0.0, 0.0, "chr1", 10000, 11000, 1, 1, 2, "label"),
            ("s3", "t1", 0.0, 4.0, 0.0, "chr1", 20000, 21000, 1, 1, 3, "label"),
            ("s4", "t2", 0.0, 0.0, 0.0, "chr1", 0, 1000, 1, 1, 1, "label"),
            ("s5", "t2", 5.0, 0.0, 0.0, "chr1", 10000, 11000, 1, 1, 2, "label"),
            ("s6", "t2", 0.0, 6.0, 0.0, "chr1", 20000, 21000, 1, 1, 3, "label"),
        ],
        names=(
            "Spot_ID",
            "Trace_ID",
            "x",
            "y",
            "z",
            "Chrom",
            "Chrom_Start",
            "Chrom_End",
            "ROI #",
            "Mask_id",
            "Barcode #",
            "label",
        ),
    )


def test_collect_distances_and_summarize_pairs():
    distances, genomic_distances = collect_distances(make_trace_table())

    assert distances[(1, 2)] == [3.0, 5.0]
    assert distances[(1, 3)] == [4.0, 6.0]
    assert np.isclose(genomic_distances[(1, 2)], 10000.0)
    assert np.isclose(genomic_distances[(1, 3)], 20000.0)

    rows = summarize_pairs(distances, genomic_distances, 100, 0.95, 1)
    pair_to_row = {(row["barcode_1"], row["barcode_2"]): row for row in rows}
    assert pair_to_row[(1, 2)]["median_physical_distance"] == 4.0
    assert pair_to_row[(1, 3)]["median_physical_distance"] == 5.0
    assert pair_to_row[(1, 2)]["n_observations"] == 2


def test_runtime_creates_plot_and_summary(tmp_path):
    input_path = tmp_path / "traces.ecsv"
    output_path = tmp_path / "distance_plot.png"
    summary_path = tmp_path / "distance_summary.csv"
    make_trace_table().write(input_path, format="ascii.ecsv")

    n_rows = runtime(
        trace_files=[str(input_path)],
        output=str(output_path),
        output_data=str(summary_path),
        output_format="png",
        bootstrap_cycles=10,
        confidence=0.95,
        fit_bootstrap_cycles=10,
        distance_unit="µm",
        genomic_unit="kb",
        min_observations=1,
        x_min=None,
        x_max=None,
        y_min=None,
        y_max=None,
    )

    assert n_rows == 3
    assert output_path.exists()
    assert summary_path.exists()
    assert "genomic_distance_kb" in summary_path.read_text()
