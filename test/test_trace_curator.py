"""Independent fixtures for empirical top3 curation and exact top5 attribution."""

from dataclasses import FrozenInstanceError, replace
from itertools import combinations
from math import log10

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from traceratops.core import trace_curator as curator
from traceratops.core.trace_curator import (
    curate_trace,
    fit_curator_model,
)


def trace(n=6, positions=None, trace_id="trace"):
    positions = (
        np.arange(n, dtype=float) if positions is None else np.asarray(positions)
    )
    return pd.DataFrame(
        {
            "Spot_ID": [f"spot-{i}" for i in range(n)],
            "Trace_ID": [trace_id] * n,
            "Barcode #": [f"locus-{i}" for i in range(n)],
            "Genomic_Position": positions,
            "x": positions.astype(float),
            "y": np.zeros(n),
            "z": np.zeros(n),
            "extra": ["preserved"] * n,
        }
    )


@pytest.fixture
def model():
    return fit_curator_model(
        [trace() for _ in range(99)], calibration_traces=[trace()] * 500
    )


def edges(values):
    return tuple(sorted(values, key=lambda edge: edge[2], reverse=True))


# Expected numbers encode the inspected empirical/top3/top5 conventions without
# importing any external implementation or retaining its scripts.
TIE_EDGES = edges(
    [
        (0, 1, 9.0),
        (0, 4, 6.0),
        (1, 3, 6.0),
        (1, 5, 6.0),
        (2, 4, 6.0),
        (3, 5, 6.0),
        (0, 3, 4.0),
        (1, 2, 4.0),
        (4, 5, 3.0),
        (2, 5, 2.0),
        (1, 4, 1.0),
        (2, 3, 1.0),
        (3, 4, 1.0),
        (0, 2, 0.0),
        (0, 5, 0.0),
    ]
)


@pytest.mark.parametrize(
    "observed,expected",
    [
        (-1.0, log10(2.5)),
        (1.0, log10(1.25)),
        (2.0, 0.0),
        (3.0, log10(1.25)),
        (5.0, log10(2.5)),
    ],
)
def test_empirical_tail_equivalence(observed, expected):
    assert curator._empirical_anomaly(
        np.array([1.0, 2.0, 2.0, 4.0]), observed
    ) == pytest.approx(expected)


def test_top3_loo_and_top5_equivalence():
    assert curator._top_cost(TIE_EDGES, 3) == 21.0
    assert curator._loo_scores(TIE_EDGES, 6, 3) == {
        0: 3.0,
        1: 3.0,
        2: 0.0,
        3: 0.0,
        4: 0.0,
        5: 0.0,
    }
    winner, primary, tied, secondary, remaining = curator._select_candidate(
        TIE_EDGES, 6
    )
    assert winner == 1
    assert tied == (0, 1)
    assert secondary == {0: 5.0, 1: 8.0}
    assert remaining == (1,)
    assert primary[winner] == 3.0
    assert curator._top_cost(TIE_EDGES, 5) == 33.0


def test_unique_top3_winner_does_not_calculate_top5(monkeypatch):
    graph = edges(
        [(0, 1, 8.0), (0, 2, 7.0), (1, 3, 6.0), (0, 4, 5.0), (2, 4, 4.0), (3, 4, 1.0)]
    )
    original = curator._loo_scores
    calls = []

    def record(*args, **kwargs):
        calls.append(args[2])
        return original(*args, **kwargs)

    monkeypatch.setattr(curator, "_loo_scores", record)
    assert curator._select_candidate(graph, 5)[0] == 0
    assert calls == [3]


def test_top5_only_scores_primary_tied_candidates(monkeypatch):
    original = curator._loo_scores
    calls = []

    def record(*args, **kwargs):
        calls.append((args[2], kwargs.get("candidates")))
        return original(*args, **kwargs)

    monkeypatch.setattr(curator, "_loo_scores", record)
    curator._select_candidate(TIE_EDGES, 6)
    assert calls == [(3, None), (5, (0, 1))]


def test_top5_tie_unresolved():
    graph = edges([(i, j, 1.0) for i, j in combinations(range(4), 2)])
    winner, _, tied, _, remaining = curator._select_candidate(graph, 4)
    assert winner is None
    assert len(tied) == len(remaining) == 4


def test_partial_top_cost_and_unavailable_loo():
    graph = ((0, 1, 2.0), (0, 2, 1.0))
    assert curator._top_cost(graph, 3) == 3.0
    scores = curator._loo_scores(graph, 3, 3)
    assert np.isnan(scores[0])
    assert scores[1] == 2.0
    assert scores[2] == 1.0
    assert np.isnan(curator._top_cost((), 3))
    assert curator._maximal_candidates({0: np.nan}) == ()


def test_fit_sorted_arrays_and_metadata():
    frames = [trace(4, positions=[0, 10, 30, 70]) for _ in range(3)]
    for factor, frame in zip([3.0, 1.0, 2.0], frames):
        frame["x"] *= factor
    with pytest.warns(UserWarning, match="may be unstable"):
        fitted = fit_curator_model(
            frames,
            minimum_reference_observations=3,
            calibration_traces=frames,
            trace_fpr=1 / 3,
        )
    np.testing.assert_array_equal(
        fitted.reference_distributions[10.0], [10.0, 20.0, 30.0]
    )
    assert set(fitted.reference_distributions) == {10.0, 20.0, 30.0, 40.0, 60.0, 70.0}
    assert fitted.n_reference_traces == 3
    assert fitted.n_reference_pairs == 18
    assert fitted.n_calibration_traces == 3
    assert fitted.calibration_source == "separate_traces"
    for sample in fitted.reference_distributions.values():
        assert np.all(np.diff(sample) >= 0)
        with pytest.raises(ValueError):
            sample.setflags(write=True)
    with pytest.raises(TypeError):
        fitted.reference_distributions[1.0] = np.array([1.0])
    with pytest.raises(FrozenInstanceError):
        fitted.trace_fpr = 0.5


@pytest.mark.parametrize(
    "fpr,expected", [(0.0, 99.0), (0.01, 98.0), (0.05, 94.0), (0.5, 49.0)]
)
def test_empirical_threshold_convention(fpr, expected):
    assert curator._empirical_threshold(np.arange(100.0), fpr) == expected


def test_actual_clean_threshold_calibration():
    base = trace(4)
    clean = []
    for offset in [0.0, 0.5, 1.0, 2.0, 3.0]:
        frame = base.copy()
        frame.loc[0, "y"] = offset
        clean.append(frame)
    with pytest.warns(UserWarning, match="may be unstable"):
        fitted = fit_curator_model([base] * 30, calibration_traces=clean, trace_fpr=0.2)
    a = log10((90 + 1) / 2)
    b = log10((60 + 1) / 2)
    c = log10((30 + 1) / 2)
    assert fitted.global_threshold == pytest.approx(a + b + c)
    assert fitted.n_calibration_traces == 5
    assert fitted.calibration_source == "separate_traces"
    assert curate_trace(clean[-1], fitted).terminal_state == "normal"


def test_minimum_support_ignores_sparse_pairs(model):
    sparse = replace(
        model, reference_distributions={1.0: np.ones(19), 2.0: np.full(20, 2.0)}
    )
    frame, b, p, xyz = curator._validate_trace(trace(4), sparse.genomic_positions)
    graph = curator._pair_anomalies(
        b, p, xyz, sparse.reference_distributions, "separation", 20
    )
    assert [(i, j) for i, j, _ in graph] == [(0, 2), (1, 3)]
    assert curate_trace(frame, sparse).terminal_state == "normal"
    unavailable = replace(model, minimum_reference_observations=1000)
    result = curate_trace(frame, unavailable)
    assert result.terminal_state == "unscoreable"
    assert result.stop_reason == "insufficient_reference"
    assert np.isnan(result.final_cost)


def test_scoring_never_sorts_reference_arrays(model, monkeypatch):
    monkeypatch.setattr(
        np, "sort", lambda *a, **k: pytest.fail("sorting during scoring")
    )
    assert curate_trace(trace(), model).terminal_state == "normal"


def test_barcode_pair_mode_and_string_barcodes():
    fitted = fit_curator_model(
        [trace(4)] * 20,
        reference_mode="barcode_pair",
        calibration_traces=[trace(4)] * 500,
    )
    assert frozenset(("locus-0", "locus-1")) in fitted.reference_distributions
    assert curate_trace(trace(4), fitted).terminal_state == "normal"


def test_genomic_mapping_is_reused_and_column_alias(model):
    frame = (
        trace()
        .drop(columns="Genomic_Position")
        .rename(columns={"Barcode #": "Barcode"})
    )
    assert curate_trace(frame, model).terminal_state == "normal"
    fitted = fit_curator_model(
        [frame] * 20,
        genomic_positions=dict(model.genomic_positions),
        calibration_traces=[frame] * 500,
    )
    assert fitted.global_threshold == model.global_threshold


def test_table_input_and_combined_reference_table():
    frames = [trace(trace_id=f"trace-{i}") for i in range(20)]
    fitted = fit_curator_model(
        Table.from_pandas(pd.concat(frames)), calibration_traces=[frames[0]] * 500
    )
    result = curate_trace(Table.from_pandas(frames[0]), fitted)
    pd.testing.assert_frame_equal(result.curated_trace, frames[0])
    assert result.n_removed == 0
    assert result.iterations[-1].stop_reason == "global_cost_normal"


def test_end_to_end_equivalent_deletion_and_rescoring(model):
    frame = trace()
    frame.loc[2, "y"] = 100.0
    frame.index = [7] * len(frame)  # Row indices cannot be used as identities.
    original = frame.copy()
    result = curate_trace(frame, model)
    # The three strongest relationships have 495, 495 and 396 samples.
    expected = 2 * log10(496 / 2) + log10(397 / 2)
    assert result.initial_cost == pytest.approx(expected)
    assert result.removed_localizations["Spot_ID"].tolist() == ["spot-2"]
    assert result.n_removed == 1
    assert result.terminal_state == "curated"
    assert result.final_cost == 0.0
    first, last = result.iterations
    assert first.candidate_spot_id == "spot-2"
    assert first.candidate_barcode == "locus-2"
    assert first.loo_top3 == pytest.approx(expected)
    assert first.n_tied_top3 == 1
    assert not first.top5_tiebreak_used
    assert first.removal_accepted
    assert last.n_localizations_before == 5
    assert not last.globally_abnormal
    pd.testing.assert_frame_equal(frame, original)
    result.curated_trace.loc[:, "extra"] = "edited"
    pd.testing.assert_frame_equal(frame, original)


@pytest.mark.parametrize(
    "order", [list(range(6)), [5, 4, 3, 2, 1, 0], [2, 5, 0, 4, 1, 3]]
)
def test_deletion_is_row_order_independent(model, order):
    frame = trace()
    frame.loc[2, "y"] = 100.0
    result = curate_trace(frame.iloc[order], model)
    assert result.removed_localizations["Spot_ID"].tolist() == ["spot-2"]
    assert result.terminal_state == "curated"


@pytest.mark.parametrize(
    "order", [list(range(6)), [5, 4, 3, 2, 1, 0], [2, 5, 0, 4, 1, 3]]
)
def test_top5_selected_identity_equivalence(model, monkeypatch, order):
    def pair_scores(barcodes, *args):
        identities = [int(b.split("-")[1]) for b in barcodes]
        index = {identity: offset for offset, identity in enumerate(identities)}
        if 1 not in index:
            return edges(
                [
                    (index[i], index[j], 0.0)
                    for i, j, _ in TIE_EDGES
                    if i in index and j in index
                ]
            )
        return edges([(index[i], index[j], s) for i, j, s in TIE_EDGES])

    monkeypatch.setattr(curator, "_pair_anomalies", pair_scores)
    result = curate_trace(trace().iloc[order], model)
    assert result.removed_localizations["Spot_ID"].tolist() == ["spot-1"]
    assert result.n_removed == 1
    assert result.terminal_state == "curated"
    decision = result.iterations[0]
    assert decision.C_top3_before == 21.0
    assert decision.loo_top3 == 3.0
    assert decision.n_tied_top3 == 2
    assert decision.top5_tiebreak_used
    assert decision.loo_top5 == 8.0
    assert decision.n_tied_top5 == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_unresolved_tie_never_chooses_row_order(model, monkeypatch, reverse):
    def symmetric(barcodes, *args):
        return edges([(i, j, 1.0) for i, j in combinations(range(len(barcodes)), 2)])

    monkeypatch.setattr(curator, "_pair_anomalies", symmetric)
    frame = trace(4)
    result = curate_trace(frame.iloc[::-1] if reverse else frame, model)
    assert result.terminal_state == "abnormal_unresolved"
    assert result.stop_reason == "ambiguous_candidate"
    assert result.n_removed == 0
    decision = result.iterations[0]
    assert decision.top5_tiebreak_used
    assert decision.n_tied_top3 == decision.n_tied_top5 == 4
    assert decision.candidate_spot_id is None


@pytest.mark.parametrize(
    "maximum,expected,state",
    [
        (0, 0, "abnormal_unresolved"),
        (1, 1, "abnormal_unresolved"),
        (3, 3, "abnormal_unresolved"),
        (4, 4, "curated"),
    ],
)
def test_multiple_removals_rescored_and_cap_respected(
    model, monkeypatch, maximum, expected, state
):
    calls = []

    def descending(barcodes, *args):
        identities = [int(b.split("-")[1]) for b in barcodes]
        calls.append(tuple(identities))
        weights = [100.0, 10.0, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0]
        return edges(
            [
                (i, j, weights[identities[i]] + weights[identities[j]])
                for i, j in combinations(range(len(barcodes)), 2)
            ]
        )

    monkeypatch.setattr(curator, "_pair_anomalies", descending)
    result = curate_trace(trace(8), model, max_removals=maximum)
    assert result.n_removed == expected
    assert result.terminal_state == state
    assert result.removed_localizations["Spot_ID"].tolist() == [
        f"spot-{i}" for i in range(expected)
    ]
    assert len(calls) == expected + 1
    assert [len(call) for call in calls] == list(range(8, 8 - expected - 1, -1))
    if state == "abnormal_unresolved":
        assert result.stop_reason == "max_removals_reached"


def test_default_three_removal_limit(model, monkeypatch):
    def descending(barcodes, *args):
        w = [10.0 ** (8 - int(b.split("-")[1])) for b in barcodes]
        return edges([(i, j, w[i] + w[j]) for i, j in combinations(range(len(w)), 2)])

    monkeypatch.setattr(curator, "_pair_anomalies", descending)
    result = curate_trace(trace(8), model)
    assert result.n_removed == 3
    assert result.stop_reason == "max_removals_reached"


def test_abnormal_minimum_trace_size(model):
    frame = trace(3)
    frame.loc[0, "y"] = 100.0
    result = curate_trace(frame, model)
    assert result.terminal_state == "abnormal_unresolved"
    assert result.stop_reason == "minimum_trace_size_reached"
    assert result.n_removed == 0


@pytest.mark.parametrize("n", [0, 1, 2])
def test_small_and_empty_traces_unscoreable(model, n):
    result = curate_trace(trace(n), model)
    assert result.terminal_state == "unscoreable"
    assert result.stop_reason == "too_few_localizations"
    assert result.n_removed == 0


def test_no_scoreable_candidate(model, monkeypatch):
    monkeypatch.setattr(curator, "_select_candidate", lambda *a: (None, {}, (), {}, ()))
    frame = trace()
    frame.loc[0, "y"] = 100.0
    result = curate_trace(frame, model)
    assert result.terminal_state == "abnormal_unresolved"
    assert result.stop_reason == "no_scoreable_candidate"


@pytest.mark.parametrize(
    "mutation,match",
    [
        (lambda f: f.assign(**{"Barcode #": ["same"] * len(f)}), "resolve duplicated"),
        (lambda f: f.assign(Spot_ID=["same"] * len(f)), "Spot_ID must be unique"),
        (lambda f: f.assign(x=np.nan), "coordinates must be finite"),
        (lambda f: f.assign(z=np.inf), "coordinates must be finite"),
        (lambda f: f.assign(x="bad"), "coordinates must be finite"),
        (lambda f: f.drop(columns="x"), "Missing required"),
        (lambda f: f.assign(Genomic_Position=0.0), "finite and unique"),
        (lambda f: f.assign(Genomic_Position=np.inf), "finite and unique"),
        (lambda f: f.assign(Genomic_Position=np.arange(len(f)) + 1.0), "disagrees"),
        (lambda f: f.assign(Trace_ID=np.arange(len(f))), "single nonmissing"),
        (lambda f: f.assign(Spot_ID=None), "missing identities"),
        (lambda f: f.assign(Chrom=["chr1", "chr2"] * 3), "one chromosome"),
    ],
)
def test_invalid_trace_rejected(model, mutation, match):
    with pytest.raises(ValueError, match=match):
        curate_trace(mutation(trace()), model)


def test_missing_genomic_coordinates_rejected():
    frame = trace().drop(columns="Genomic_Position")
    with pytest.raises(ValueError, match="genomic coordinates"):
        fit_curator_model([frame] * 20)


def test_inconsistent_positions_between_reference_traces_rejected():
    changed = trace()
    changed["Genomic_Position"] += 1
    with pytest.raises(ValueError, match="disagrees"):
        fit_curator_model([trace(), changed], calibration_traces=[trace()] * 500)


def test_missing_model_position_for_new_barcode_rejected(model):
    frame = trace().drop(columns="Genomic_Position")
    frame.loc[0, "Barcode #"] = "unknown"
    with pytest.raises(ValueError, match="genomic coordinates"):
        curate_trace(frame, model)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"trace_fpr": -1},
        {"trace_fpr": 0},
        {"trace_fpr": 1},
        {"trace_fpr": np.nan},
        {"minimum_reference_observations": 0},
        {"minimum_reference_observations": 1.5},
        {"reference_mode": "invalid"},
    ],
)
def test_invalid_fit_settings(kwargs):
    with pytest.raises(ValueError):
        fit_curator_model([trace()], **kwargs)


@pytest.mark.parametrize("maximum", [-1, 1.5, True])
def test_invalid_removal_limit(model, maximum):
    with pytest.raises(ValueError):
        curate_trace(trace(), model, max_removals=maximum)


def test_no_scoreable_calibration_fails_clearly():
    with pytest.raises(ValueError, match=r"Too few scoreable calibration traces \(0\)"):
        fit_curator_model([trace()] * 100, minimum_reference_observations=1000)


def test_calibration_excludes_unscoreable_traces(model):
    fitted = fit_curator_model(
        [trace()] * 20, calibration_traces=[trace()] * 500 + [trace(2)]
    )
    assert fitted.n_calibration_traces == 500
    assert fitted.n_calibration_unscoreable == 1


def test_streaming_references_require_separate_calibration():
    with pytest.raises(ValueError, match="One-pass"):
        fit_curator_model(iter([trace()] * 20))
    fitted = fit_curator_model(
        iter([trace()] * 20), calibration_traces=iter([trace()] * 500)
    )
    assert fitted.n_reference_traces == 20
    assert fitted.n_calibration_traces == 500


def test_real_coordinates_allow_multiple_removals():
    base = trace(10)
    fitted = fit_curator_model([base] * 99, calibration_traces=[base] * 500)
    frame = base.copy()
    frame.loc[0, "y"] = 100.0
    frame.loc[5, "y"] = 600.0
    result = curate_trace(frame, fitted)
    assert result.removed_localizations["Spot_ID"].tolist() == ["spot-5", "spot-0"]
    assert result.n_removed == 2
    assert result.terminal_state == "curated"
    assert result.final_cost == 0.0
    assert [i.n_localizations_before for i in result.iterations] == [10, 9, 8]


def test_real_coordinates_top3_tie_resolved_by_top5():
    base = trace(10)
    fitted = fit_curator_model([base] * 99, calibration_traces=[base] * 500)
    frame = base.copy()
    frame.loc[1, "y"] = 200.0
    frame.loc[5, "y"] = 600.0
    result = curate_trace(frame, fitted)
    assert result.removed_localizations["Spot_ID"].tolist() == ["spot-5", "spot-1"]
    decision = result.iterations[0]
    assert decision.n_tied_top3 == 2
    assert decision.top5_tiebreak_used
    assert decision.n_tied_top5 == 1
    assert result.terminal_state == "curated"


def test_equal_threshold_is_normal(model):
    frame = trace()
    frame.loc[2, "y"] = 100.0
    cost = curate_trace(frame, model, max_removals=0).initial_cost
    equal = replace(model, global_threshold=cost)
    result = curate_trace(frame, equal)
    assert result.terminal_state == "normal"
    assert result.n_removed == 0


def test_model_constructor_rejects_unsorted_reference(model):
    with pytest.raises(ValueError, match="sorted"):
        replace(model, reference_distributions={1.0: np.array([2.0, 1.0])})


def test_calibration_does_not_calculate_loo(monkeypatch):
    monkeypatch.setattr(
        curator, "_loo_scores", lambda *a, **k: pytest.fail("LOO during fitting")
    )
    with pytest.warns(UserWarning, match="may be unstable"):
        fitted = fit_curator_model([trace()] * 100)
    assert fitted.n_calibration_traces == 100


def test_conflicting_barcode_columns_rejected(model):
    frame = trace().assign(Barcode="different")
    with pytest.raises(ValueError, match="columns disagree"):
        curate_trace(frame, model)


def test_chromosome_consistency():
    base = trace().assign(Chrom="chr1")
    fitted = fit_curator_model([base] * 20, calibration_traces=[base] * 500)
    with pytest.raises(ValueError, match="chromosome disagrees"):
        curate_trace(base.assign(Chrom="chr2"), fitted)
    with pytest.raises(ValueError, match="share one chromosome"):
        fit_curator_model(
            [base, base.assign(Chrom="chr2")], calibration_traces=[base] * 500
        )


def test_combined_reference_table_requires_trace_id():
    with pytest.raises(ValueError, match="Trace_ID is required"):
        fit_curator_model(trace().drop(columns="Trace_ID"))


@pytest.mark.parametrize("position", [float("inf"), float("nan"), "invalid"])
def test_unused_invalid_mapping_position_rejected(position):
    with pytest.raises(ValueError, match="Genomic mapping"):
        fit_curator_model([trace()] * 20, genomic_positions={"unused": position})


def test_top_cost_ignores_nonfinite_edges():
    graph = ((0, 1, float("inf")), (0, 2, 2.0), (1, 2, float("nan")))
    assert curator._top_cost(graph, 3) == 2.0


def test_numeric_spot_identity_is_not_coerced_by_mixed_numeric_row():
    base = trace().drop(columns=["Trace_ID", "extra"])
    base["Barcode #"] = np.arange(len(base))
    base["Spot_ID"] = np.arange(2**53 + 1, 2**53 + 1 + len(base), dtype=np.int64)
    fitted = fit_curator_model([base] * 20, calibration_traces=[base] * 500)
    frame = base.copy()
    frame.loc[2, "y"] = 100.0
    result = curate_trace(frame, fitted)
    identity = base["Spot_ID"].iloc[2]
    assert result.iterations[0].candidate_spot_id == identity
    assert isinstance(result.iterations[0].candidate_spot_id, np.integer)
    assert result.removed_localizations["Spot_ID"].tolist() == [identity]


@pytest.mark.parametrize("n_folds", [2, 3, 5, 11])
def test_crossfit_disjoint_balanced_and_deterministic(monkeypatch, n_folds):
    frames = tuple(trace(3, trace_id=f"trace-{i}") for i in range(11))
    original_fit = curator._fit_reference
    original_costs = curator._calibration_costs
    fits, held_out = [], []
    active_training = set()

    def record_fit(traces, *args):
        nonlocal active_training
        training = list(traces)
        active_training = {f["Trace_ID"].iloc[0] for f in training}
        fits.append(tuple(f["Trace_ID"].iloc[0] for f in training))
        return original_fit(training, *args)

    def record_costs(traces, *args):
        fold = list(traces)
        identities = tuple(f["Trace_ID"].iloc[0] for f in fold)
        assert set(identities).isdisjoint(active_training)
        held_out.append(identities)
        yield from original_costs(fold, *args)

    monkeypatch.setattr(curator, "_fit_reference", record_fit)
    monkeypatch.setattr(curator, "_calibration_costs", record_costs)
    first = fit_curator_model(
        frames,
        n_crossfit_folds=n_folds,
        trace_fpr=0.5,
        minimum_reference_observations=1,
    )
    first_fits, first_held = list(fits), list(held_out)
    fits.clear()
    held_out.clear()
    second = fit_curator_model(
        frames,
        n_crossfit_folds=n_folds,
        trace_fpr=0.5,
        minimum_reference_observations=1,
    )
    assert fits == first_fits
    assert held_out == first_held
    assert len(fits) == n_folds + 1  # Exactly K fits, then one final fit.
    assert fits[-1] == tuple(f"trace-{i}" for i in range(11))
    assert sorted(identity for fold in held_out for identity in fold) == sorted(
        fits[-1]
    )
    assert max(map(len, held_out)) - min(map(len, held_out)) <= 1
    for fold in range(n_folds):
        assert held_out[fold] == tuple(f"trace-{i}" for i in range(fold, 11, n_folds))
    assert first.global_threshold == second.global_threshold == 0.0
    assert first.n_crossfit_folds == n_folds
    assert first.calibration_source == "crossfit"
    assert first.n_calibration_traces == 11
    assert first.n_calibration_unscoreable == 0
    assert first.n_reference_traces == 11
    for key in first.reference_distributions:
        np.testing.assert_array_equal(
            first.reference_distributions[key], second.reference_distributions[key]
        )


def test_default_threshold_is_crossfitted_and_differs_from_self_inclusion():
    frames = []
    for factor in range(1, 101):
        frame = trace(3)
        frame["x"] *= factor
        frames.append(frame)
    with pytest.warns(UserWarning, match="100.*1% tail.*may be unstable"):
        fitted = fit_curator_model(frames, reference_mode="barcode_pair")
    # Each training fold contains 80 samples per barcode pair. Both extrema
    # fall outside their held-out references, giving finite-sample p=2/81.
    expected_crossfit = 3 * log10(81 / 2)
    expected_in_sample = 3 * log10(101 / 4)
    assert fitted.global_threshold == pytest.approx(expected_crossfit)
    values = list(
        curator._calibration_costs(
            frames,
            fitted.reference_distributions,
            fitted.genomic_positions,
            None,
            "barcode_pair",
            20,
        )
    )
    in_sample = curator._empirical_threshold(values, 0.01)
    assert in_sample == pytest.approx(expected_in_sample)
    assert fitted.global_threshold != in_sample
    assert fitted.n_crossfit_folds == 5
    assert fitted.n_reference_pairs == 300
    assert all(len(sample) == 100 for sample in fitted.reference_distributions.values())


@pytest.mark.parametrize("folds", [0, 1, -2, 12, 2.5, True])
def test_invalid_crossfit_fold_count_rejected(folds):
    with pytest.raises(ValueError, match="n_crossfit_folds"):
        fit_curator_model([trace(3)] * 11, n_crossfit_folds=folds, trace_fpr=0.5)


@pytest.mark.parametrize("n_folds", [2, 5])
def test_insufficient_crossfitted_sample_rejected(n_folds):
    with pytest.raises(ValueError, match="Too few scoreable.*99.*at least 100"):
        fit_curator_model([trace(3)] * 99, n_crossfit_folds=n_folds)


def test_calibration_sample_rule_counts_finite_scores_not_input_size():
    frames = [trace(3)] * 99 + [trace(2)] * 10
    with pytest.raises(ValueError, match="Too few scoreable.*99.*at least 100"):
        fit_curator_model(frames)


@pytest.mark.parametrize(
    "fpr,n,minimum", [(0.01, 99, 100), (0.05, 19, 20), (0.3, 3, 4)]
)
def test_tail_sample_requirement_generalizes(fpr, n, minimum):
    with pytest.raises(
        ValueError, match=f"at least {minimum} finite calibration scores"
    ):
        fit_curator_model(
            [trace(3)] * 20, calibration_traces=[trace(3)] * n, trace_fpr=fpr
        )


@pytest.mark.parametrize("fpr,n", [(0.01, 100), (0.01, 499), (0.05, 20)])
def test_low_but_valid_sample_warns(fpr, n):
    with pytest.warns(UserWarning, match=f"Only {n}.*estimable but may be unstable"):
        fitted = fit_curator_model([trace(3)] * n, trace_fpr=fpr)
    assert fitted.n_calibration_traces == n


@pytest.mark.parametrize("fpr,n", [(0.01, 500), (0.05, 100)])
def test_five_expected_tail_observations_do_not_warn(fpr, n):
    import warnings

    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter("always")
        fitted = fit_curator_model([trace(3)] * n, trace_fpr=fpr)
    assert not observed
    assert fitted.n_calibration_traces == n


def test_independent_calibration_bypasses_crossfit(monkeypatch):
    original_fit = curator._fit_reference
    fit_calls = []

    def record_fit(traces, *args):
        fit_calls.append(1)
        return original_fit(traces, *args)

    monkeypatch.setattr(curator, "_fit_reference", record_fit)
    monkeypatch.setattr(
        curator,
        "_fold_traces",
        lambda *a, **k: pytest.fail("cross-fitting independent controls"),
    )
    base = trace(3)
    controls = base.copy()
    controls["x"] *= 2
    fitted = fit_curator_model(
        [base] * 20, calibration_traces=[controls] * 500, n_crossfit_folds=0
    )
    assert fit_calls == [1]
    assert fitted.calibration_source == "separate_traces"
    assert fitted.n_crossfit_folds is None
    assert fitted.global_threshold == pytest.approx(2 * log10(41 / 2) + log10(21 / 2))
    assert fitted.n_reference_traces == 20
    assert fitted.n_calibration_traces == 500


def test_final_reference_support_diagnostics():
    fitted = fit_curator_model(
        [trace(4), trace(3), trace(3)],
        minimum_reference_observations=5,
        calibration_traces=[trace(3)] * 500,
    )
    assert fitted.n_reference_traces == 3
    assert fitted.n_reference_pairs == 12
    assert fitted.n_reference_keys == 3
    assert fitted.n_supported_reference_keys == 1
    assert fitted.median_reference_observations == 4.0
    assert fitted.min_reference_observations_actual == 1
    assert fitted.fraction_reference_keys_below_minimum == pytest.approx(2 / 3)


def test_compact_pair_buffers_match_float_list_reference(monkeypatch):
    from array import array

    buffers = []

    def record_buffer(typecode):
        buffer = array(typecode)
        buffers.append(buffer)
        return buffer

    monkeypatch.setattr(curator, "array", record_buffer)
    frames = [trace(4, positions=[0, 10, 30, 70]) for _ in range(3)]
    expected = {}
    for factor, frame in zip([3.0, 1.0, 2.0], frames):
        frame["x"] *= factor
        xyz = frame[["x", "y", "z"]].to_numpy()
        for i, j in combinations(range(4), 2):
            key = float(
                abs(
                    frame["Genomic_Position"].iloc[i]
                    - frame["Genomic_Position"].iloc[j]
                )
            )
            expected.setdefault(key, []).append(float(np.linalg.norm(xyz[i] - xyz[j])))
    references, _, _, n_traces, n_pairs = curator._fit_reference(
        frames, {}, "separation"
    )
    assert len(buffers) == len(expected)
    assert all(
        isinstance(buffer, array) and buffer.typecode == "d" for buffer in buffers
    )
    assert sum(map(len, buffers)) == n_pairs == 18
    assert n_traces == 3
    for key, values in expected.items():
        np.testing.assert_array_equal(references[key], np.sort(values))


def test_fold_references_are_released_before_next_fit(monkeypatch):
    import weakref

    original_fit = curator._fit_reference
    previous_arrays = []

    def record_fit(*args):
        assert all(reference() is None for reference in previous_arrays)
        result = original_fit(*args)
        previous_arrays[:] = [weakref.ref(sample) for sample in result[0].values()]
        return result

    monkeypatch.setattr(curator, "_fit_reference", record_fit)
    fitted = fit_curator_model([trace(3)] * 100, trace_fpr=0.05)
    assert fitted.n_crossfit_folds == 5
    assert all(reference() is None for reference in previous_arrays)
    for sample in fitted.reference_distributions.values():
        with pytest.raises(ValueError):
            sample.setflags(write=True)


@pytest.mark.parametrize("as_table", [False, True])
def test_combined_crossfit_uses_first_appearance_order(as_table):
    identifiers = ["z", "a", "m", "b", "c", "x", "d", "w", "e", "v"]
    frame = pd.concat([trace(3, trace_id=identity) for identity in identifiers])
    population = Table.from_pandas(frame) if as_table else frame
    # Main API converts Astropy once; helper receives pandas for combined folds.
    prepared = population.to_pandas() if as_table else population
    for fold in range(3):
        held = list(curator._fold_traces(prepared, fold, 3, held_out=True))
        assert [trace["Trace_ID"].iloc[0] for trace in held] == identifiers[fold::3]
    fitted = fit_curator_model(
        population, n_crossfit_folds=3, trace_fpr=0.5, minimum_reference_observations=1
    )
    assert fitted.n_reference_traces == fitted.n_calibration_traces == 10


def test_combined_astropy_population_is_converted_once(monkeypatch):
    population = Table.from_pandas(
        pd.concat([trace(3, trace_id=f"t-{i}") for i in range(10)])
    )
    original = Table.to_pandas
    calls = []

    def record(table, *args, **kwargs):
        calls.append(1)
        return original(table, *args, **kwargs)

    monkeypatch.setattr(Table, "to_pandas", record)
    fit_curator_model(population, trace_fpr=0.5, minimum_reference_observations=1)
    assert calls == [1]


def test_collection_without_trace_ids_supports_crossfit():
    frames = tuple(trace(3).drop(columns="Trace_ID") for _ in range(10))
    fitted = fit_curator_model(frames, trace_fpr=0.5, minimum_reference_observations=1)
    assert fitted.calibration_source == "crossfit"
    assert fitted.n_calibration_traces == 10


def test_crossfit_unscoreable_count_and_reference_support():
    frames = [trace(3)] * 100 + [trace(2)] * 3
    with pytest.warns(UserWarning, match="Only 100"):
        fitted = fit_curator_model(frames)
    assert fitted.n_calibration_traces == 100
    assert fitted.n_calibration_unscoreable == 3
    assert fitted.n_reference_traces == 103
    assert fitted.n_reference_pairs == 303
    assert fitted.n_reference_keys == fitted.n_supported_reference_keys == 2
    assert fitted.min_reference_observations_actual == 100
    assert fitted.median_reference_observations == 151.5
    assert fitted.fraction_reference_keys_below_minimum == 0.0


def test_crossfit_support_uses_training_fold_not_final_population():
    # Final all-data pairs would have 100 observations; fold references have
    # only 80. No held-out cost may use the final reference to pass this gate.
    with pytest.raises(ValueError, match=r"Too few scoreable calibration traces \(0\)"):
        fit_curator_model(
            [trace(3)] * 100,
            reference_mode="barcode_pair",
            minimum_reference_observations=90,
        )
