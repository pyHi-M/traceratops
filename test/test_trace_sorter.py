import numpy as np
from astropy.table import Table

from traceratops.core.chromatin_trace_table import ChromatinTraceTable
from traceratops.trace_sorter import sort_trace_table


def test_sort_trace_table_by_spots_and_radius_of_gyration():
    table = ChromatinTraceTable()
    table.data = Table(
        {
            "Trace_ID": ["small", "large", "large", "large", "wide", "wide"],
            "Barcode #": [1, 1, 2, 3, 1, 2],
            "x": [0.0, 0.0, 0.1, 0.2, 0.0, 10.0],
            "y": np.zeros(6),
            "z": np.zeros(6),
        }
    )

    sort_trace_table(table, "number_of_spots")
    assert list(table.data["Trace_ID"]) == [
        "large",
        "large",
        "large",
        "wide",
        "wide",
        "small",
    ]

    sort_trace_table(table, "radius_of_gyration", ascending=True)
    assert list(table.data["Trace_ID"]) == [
        "small",
        "large",
        "large",
        "large",
        "wide",
        "wide",
    ]
