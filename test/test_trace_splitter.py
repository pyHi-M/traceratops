import numpy as np
from astropy.table import Table

from traceratops import trace_splitter
from traceratops.core.chromatin_trace_table import ChromatinTraceTable


def _table(traces):
    rows = []
    for trace_id, xs in traces.items():
        rows.extend((trace_id, i, x, 0.0, 0.0) for i, x in enumerate(xs))
    result = ChromatinTraceTable()
    result.data = Table(rows=rows, names=("Trace_ID", "Barcode #", "x", "y", "z"))
    return result


def test_argument_defaults_preserve_existing_behavior():
    args = trace_splitter.parse_arguments().parse_args([])
    assert args.split_all is False
    assert args.clustering_method == "kmeans"
    assert args.num_clusters == 2


def test_split_all_applies_kmeans_to_otherwise_unselected_traces(monkeypatch):
    table = _table({"aaaa": [0, 1, 10], "bbbb": [0, 1, 10]})
    ids = iter(("c001", "c002", "c003", "c004"))
    monkeypatch.setattr(trace_splitter, "generate_unique_id", lambda: next(ids))

    trace_splitter.split_large_traces(table, 1.0, 2, split_all=True)

    assert set(table.data["Trace_ID"]) == {"c001", "c002", "c003", "c004"}


def test_default_radius_filter_only_splits_large_outlier(monkeypatch):
    table = _table(
        {
            "aaaa": [0, 0.1, 0.2],
            "bbbb": [0, 0.1, 0.2],
            "wide": [0, 10, 20],
        }
    )
    ids = iter(("c001", "c002"))
    monkeypatch.setattr(trace_splitter, "generate_unique_id", lambda: next(ids))

    trace_splitter.split_large_traces(table, 1.0, 2)

    assert list(table.data["Trace_ID"]).count("aaaa") == 3
    assert list(table.data["Trace_ID"]).count("bbbb") == 3
    assert set(table.data["Trace_ID"][-3:]) == {"c001", "c002"}


def test_hdbscan_parameters_and_noise_preservation(monkeypatch):
    table = _table({"orig": [0, 1, 9, 10, 100]})
    captured = {}

    class FakeHDBSCAN:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def fit_predict(self, coords):
            return np.array([0, 0, 1, 1, -1])

    monkeypatch.setattr(trace_splitter, "HDBSCAN", FakeHDBSCAN)
    ids = iter(("one1", "two2"))
    monkeypatch.setattr(trace_splitter, "generate_unique_id", lambda: next(ids))

    trace_splitter.split_large_traces(
        table,
        split_all=True,
        clustering_method="hdbscan",
        min_cluster_size=2,
        min_samples=4,
        cluster_selection_method="leaf",
        cluster_selection_epsilon=2.5,
        allow_single_cluster=True,
    )

    assert captured == {
        "min_cluster_size": 2,
        "min_samples": 4,
        "metric": "euclidean",
        "cluster_selection_method": "leaf",
        "cluster_selection_epsilon": 2.5,
        "allow_single_cluster": True,
    }
    assert list(table.data["Trace_ID"]) == ["one1", "one1", "two2", "two2", "orig"]


def test_hdbscan_epsilon_is_estimated_when_omitted(monkeypatch):
    table = _table({"orig": [0, 1, 9, 10]})
    captured = {}

    class FakeHDBSCAN:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def fit_predict(self, coords):
            return np.full(len(coords), -1)

    monkeypatch.setattr(trace_splitter, "HDBSCAN", FakeHDBSCAN)
    monkeypatch.setattr(trace_splitter, "estimate_hdbscan_epsilon", lambda groups: 3.25)
    trace_splitter.split_large_traces(
        table, split_all=True, clustering_method="hdbscan"
    )
    assert captured["cluster_selection_epsilon"] == 3.25
