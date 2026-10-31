import numpy as np
import pytest
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


def _assert_tables_equal(actual, expected):
    """Compare Astropy tables exactly while treating paired NaNs as equal."""
    assert actual.colnames == expected.colnames
    for name in actual.colnames:
        actual_values = np.asarray(actual[name])
        expected_values = np.asarray(expected[name])
        assert actual_values.dtype == expected_values.dtype
        if np.issubdtype(actual_values.dtype, np.inexact):
            np.testing.assert_allclose(
                actual_values,
                expected_values,
                rtol=0,
                atol=0,
                equal_nan=True,
                err_msg=f"Column {name!r} differs",
            )
        else:
            np.testing.assert_array_equal(
                actual_values,
                expected_values,
                err_msg=f"Column {name!r} differs",
            )


def test_argument_defaults_preserve_existing_behavior():
    args = trace_splitter.parse_arguments().parse_args([])
    assert args.split_all is False
    assert args.clustering_method == "kmeans"
    assert args.num_clusters == 2
    assert args.history_mode == "multi"
    assert args.distance_score == "residual"
    assert args.one_polymer_ambiguity_mode == "global"
    assert args.duplicate_minimum_confidence is None
    assert args.candidate_diagnostics_output is None
    assert args.quiet is False


def test_quiet_disables_progress_bars(monkeypatch):
    table = _table({"aaaa": [0, 1, 10], "bbbb": [0, 1, 10]})
    calls = []

    def fake_tqdm(iterable, **kwargs):
        calls.append(kwargs)
        return iterable

    monkeypatch.setattr(trace_splitter, "tqdm", fake_tqdm)

    trace_splitter.split_large_traces(table, split_all=True, show_progress=False)

    assert [call["desc"] for call in calls] == [
        "Computing trace radii",
        "Applying kmeans",
    ]
    assert all(call["disable"] for call in calls)


def test_method_and_legacy_clustering_method_are_aliases():
    parser = trace_splitter.parse_arguments()
    assert parser.parse_args(["--method", "beam"]).clustering_method == "beam"
    assert (
        parser.parse_args(["--clustering-method", "hdbscan"]).clustering_method
        == "hdbscan"
    )


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
    table.data["x"] = table.data["x"].astype(np.float32)
    table.data["y"] = table.data["y"].astype(np.float32)
    table.data["z"] = table.data["z"].astype(np.float32)
    captured = {}

    class FakeHDBSCAN:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def fit_predict(self, coords):
            captured["coords_dtype"] = coords.dtype
            captured["coords_contiguous"] = coords.flags.c_contiguous
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
        "coords_dtype": np.dtype(np.float64),
        "coords_contiguous": True,
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
    monkeypatch.setattr(
        trace_splitter,
        "estimate_hdbscan_epsilon",
        lambda groups, show_progress=True: 3.25,
    )
    trace_splitter.split_large_traces(
        table, split_all=True, clustering_method="hdbscan"
    )
    assert captured["cluster_selection_epsilon"] == 3.25


def test_hdbscan_retries_known_sklearn_epsilon_failure(monkeypatch, capsys):
    calls = []

    class FakeHDBSCAN:
        def __init__(self, **kwargs):
            calls.append(kwargs)

        def fit_predict(self, coords):
            if len(calls) == 1:
                raise TypeError(
                    "only 0-dimensional arrays can be converted to Python scalars"
                )
            return np.array([0, 0, 1, 1])

    monkeypatch.setattr(trace_splitter, "HDBSCAN", FakeHDBSCAN)
    parameters = {
        "min_cluster_size": 2,
        "min_samples": 2,
        "metric": "euclidean",
        "cluster_selection_method": "eom",
        "cluster_selection_epsilon": 0.884,
        "allow_single_cluster": True,
    }

    labels = trace_splitter._fit_hdbscan(np.zeros((4, 3)), parameters)

    assert list(labels) == [0, 0, 1, 1]
    assert calls[0]["cluster_selection_epsilon"] == 0.884
    assert calls[1]["cluster_selection_epsilon"] == 0.0
    assert "retrying this trace" in capsys.readouterr().out


def test_hdbscan_does_not_hide_unrelated_type_error(monkeypatch):
    class FakeHDBSCAN:
        def __init__(self, **kwargs):
            pass

        def fit_predict(self, coords):
            raise TypeError("unrelated failure")

    monkeypatch.setattr(trace_splitter, "HDBSCAN", FakeHDBSCAN)
    parameters = {
        "cluster_selection_epsilon": 1.0,
    }

    with pytest.raises(TypeError, match="unrelated failure"):
        trace_splitter._fit_hdbscan(np.zeros((4, 3)), parameters)


def _duplicate_table():
    result = ChromatinTraceTable()
    result.data = Table(
        rows=[
            ("original-0", "trace-a", 1, 0.0, 0.0, 0.0),
            ("genuine-2", "trace-a", 2, 1.0, 0.0, 0.0),
            ("offtarget-2", "trace-a", 2, 20.0, 0.0, 0.0),
            ("original-3", "trace-a", 3, 2.0, 0.0, 0.0),
        ],
        names=("Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z"),
    )
    return result


def test_default_and_explicit_global_modes_reproduce_the_same_output():
    default = _duplicate_table()
    explicit = _duplicate_table()

    default_diagnostics = trace_splitter.resolve_traces(default)
    explicit_diagnostics = trace_splitter.resolve_traces(
        explicit, one_polymer_ambiguity_mode="global"
    )

    _assert_tables_equal(default.data, explicit.data)
    _assert_tables_equal(default_diagnostics, explicit_diagnostics)


def test_candidate_diagnostics_preserve_spot_ids_and_barcode_identity():
    table = _duplicate_table()
    candidate_rows = []

    trace_splitter.resolve_traces(
        table,
        one_polymer_ambiguity_mode="candidate",
        duplicate_minimum_confidence=0.0,
        candidate_diagnostics=candidate_rows,
        model_minimum_observations=1,
    )

    assert {row["Spot_ID"] for row in candidate_rows} == {
        "genuine-2",
        "offtarget-2",
    }
    assert {row["Barcode #"] for row in candidate_rows} == {2}
    assert {row["Input_Trace_ID"] for row in candidate_rows} == {"trace-a"}
    assert sorted(row["candidate_rank"] for row in candidate_rows) == [1, 2]
    assert sum(row["selected_candidate"] for row in candidate_rows) == 1


def test_candidate_mode_does_not_change_two_polymer_output_or_diagnostics(
    monkeypatch,
):
    rows = []
    for barcode in range(1, 5):
        rows.append((f"left-{barcode}", "doublet", barcode, float(barcode), 0.0, 0.0))
        rows.append((f"right-{barcode}", "doublet", barcode, float(barcode), 10.0, 0.0))

    def make_table():
        result = ChromatinTraceTable()
        result.data = Table(
            rows=rows,
            names=("Spot_ID", "Trace_ID", "Barcode #", "x", "y", "z"),
        )
        return result

    identifiers = iter(("global-one", "global-two", "global-one", "global-two"))
    monkeypatch.setattr(trace_splitter, "generate_unique_id", lambda: next(identifiers))
    global_table = make_table()
    candidate_table = make_table()
    arguments = {
        "thresholds": trace_splitter.ClassificationThresholds(2, 0.2, "both"),
        "beam_width": 200,
        "rejection_cost": 20,
        "minimum_confidence": 0.001,
        "model_minimum_observations": 1,
    }

    global_diagnostics = trace_splitter.resolve_traces(global_table, **arguments)
    candidate_diagnostics = trace_splitter.resolve_traces(
        candidate_table, one_polymer_ambiguity_mode="candidate", **arguments
    )

    _assert_tables_equal(global_table.data, candidate_table.data)
    _assert_tables_equal(global_diagnostics, candidate_diagnostics)
    assert global_diagnostics.meta == candidate_diagnostics.meta
