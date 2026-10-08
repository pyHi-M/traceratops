"""Cycle-5 unit cases; no synthetic benchmark sweep or scientific data generation."""

import importlib.util
import sys
import weakref
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "analyze_trace_curator_tie_resolution.py"
)
SPEC = importlib.util.spec_from_file_location("cycle5", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)
C4 = analysis.C4
THRESHOLDS = {
    "global_threshold": 5.0,
    "absolute_threshold": 2.0,
    "relative_threshold": 0.5,
}


def trace(n=6, corrupted=()):
    return pd.DataFrame(
        {
            "Trace_ID": "t",
            "Spot_ID": [f"p{i}" for i in range(n)],
            "Barcode": [90, 2, 40, 7, 18, 8, 6][:n],
            "Genomic_Position": np.arange(n) * 10.0,
            "x": np.arange(n, dtype=float),
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": [i in corrupted for i in range(n)],
        }
    )


def state(frame=None, loo=None, edges=(), cost=8.0):
    frame = trace() if frame is None else frame.reset_index(drop=True)
    loo = (
        np.array([3.0, 3.0, 0.0, 0.0, 0.0, 0.0])
        if loo is None
        else np.asarray(loo, float)
    )
    return C4.TraceScores(frame, tuple(edges), cost, loo, cost - loo)


def stable_builder(frame, reference, mode, fallback):
    # One unique winner per state, independent of row order; three removals then normal.
    remaining = frame["Spot_ID"].astype(str).tolist()
    first = min(remaining)
    scores = np.array([2.0 if spot == first else 0.0 for spot in remaining])
    cost = 8.0 if len(frame) > 3 else 4.0
    return state(frame, scores, cost=cost)


def meta(k=2, dataset="primary"):
    return {
        "simulation_id": "s",
        "condition": "condition",
        "seed": 11,
        "n_barcodes": 25,
        "K": k,
        "detection_efficiency": 0.5,
        "displacement": 0.4,
        "dataset_type": dataset,
        "evaluation_seed": 11,
        "reference_seeds": "12,13",
        "calibration_seeds": "12,13",
        "nested_reference_seed_sets": "12:13;13:12",
    }


def test_exact_tie_detection_and_nonfinite_scores():
    assert analysis.primary_candidates(state()).tolist() == [0, 1]
    close = state(loo=[3.0, np.nextafter(3.0, 0.0), np.nan, np.inf, -np.inf, 0.0])
    assert analysis.primary_candidates(close).tolist() == [0]
    assert not len(analysis.primary_candidates(state(loo=[np.nan] * 6)))


@pytest.mark.parametrize("method", analysis.METHODS)
def test_secondary_ties_never_use_row_order(method):
    edges = [(a, b, 1.0) for a in range(6) for b in range(a + 1, 6)]
    current = state(edges=edges)

    def score(frame, *args):
        return state(frame, np.ones(len(frame)))

    for frame in (current.trace, current.trace.iloc[::-1]):
        ordered = frame.reset_index(drop=True)
        tied_state = state(ordered, np.full(6, 3.0), edges=edges)
        result = analysis.resolve_tie(tied_state, method, score, object())
        assert result.winner is None


def test_top5_only_ranks_primary_tied_set():
    current = state(
        edges=[(0, 2, 5.0), (0, 3, 4.0), (1, 2, 3.0), (2, 3, 2.0), (2, 4, 1.0)]
    )
    result = analysis.resolve_tie(current, "loo_top5", None, None)
    assert result.candidates == (0, 1)
    assert result.scores == (9.0, 3.0)
    assert (
        result.winner == 0
    )  # Node 2 has larger top5 drop but is outside the primary tie.
    with pytest.raises(ValueError, match="exact primary"):
        analysis.resolve_tie(state(loo=[4, 3, 2, 1, 0, 0]), "loo_top5", None, None)


def test_k1_neighbors_follow_genomic_positions_and_missing_flanks():
    frame = trace()
    frame["Genomic_Position"] = [30.0, 0.0, 50.0, 10.0, 40.0, 20.0]
    current = state(frame, edges=[(0, 5, 2.0), (0, 4, 4.0), (1, 3, 1.0), (0, 1, 99.0)])
    result = analysis.resolve_tie(current, "k1_local", None, None)
    assert result.left_spots == ("p5", "")
    assert result.right_spots == ("p4", "p3")
    assert result.scores == (3.0, 1.0)
    assert np.isnan(result.left[1])
    empty = analysis.resolve_tie(state(frame), "k1_local", None, None)
    assert empty.n_scoreable == 0 and empty.winner is None


def test_k1_does_not_skip_unscoreable_nearest_detected_neighbor():
    current = state(edges=[(0, 2, 50.0), (1, 2, 2.0)])
    result = analysis.resolve_tie(current, "k1_local", None, None)
    assert np.isnan(
        result.scores[0]
    )  # p1 is nearest, so distant p2 is not substituted.
    assert result.scores[1] == 2.0 and result.winner == 1
    assert result.status == "partially_scoreable"


def test_lookahead_rescores_each_reduced_trace_and_uses_next_loo():
    calls = []

    def rescored(frame, *args):
        missing = set(trace()["Spot_ID"]) - set(frame["Spot_ID"])
        calls.append(missing)
        best = 7.0 if missing == {"p1"} else 2.0
        return state(frame, [best] + [0.0] * (len(frame) - 1), cost=5.0)

    result = analysis.resolve_tie(state(), "one_step_lookahead", rescored, None)
    assert calls == [{"p0"}, {"p1"}]
    assert result.scores == (2.0, 7.0) and result.winner == 1

    def unscoreable(frame, *args):
        return state(frame, [np.nan] * len(frame))

    result = analysis.resolve_tie(state(), "one_step_lookahead", unscoreable, None)
    assert result.winner is None and result.status == "unscoreable"


@pytest.mark.parametrize("method", analysis.METHODS)
@pytest.mark.parametrize("corrupted", [(), (0,), (0, 1, 2)])
def test_non_tied_decisions_and_clean_controls_match_cycle4(method, corrupted):
    frame = trace(corrupted=corrupted)
    expected = C4.run_policy(
        frame, object(), THRESHOLDS, "C_cap", scorer=stable_builder
    )[1]

    def forbidden(*args):
        pytest.fail("Secondary method invoked outside a primary tie")

    actual, counts = analysis.run_integrated(
        frame, object(), THRESHOLDS, method, scorer=stable_builder, resolver=forbidden
    )
    assert actual == expected
    assert counts["tie_invocations"] == 0


def test_cap_stopping_precedes_secondary_tie_at_three_removals():
    frame = trace(n=7, corrupted=(0, 1, 2))

    def builder(current, *args):
        if len(current) == 4:
            return state(current, [3.0] * 4, cost=8.0)
        return stable_builder(current, *args)

    def forbidden(*args):
        pytest.fail("Tie breaker invoked after reaching removal cap")

    for method in analysis.METHODS:
        result, counts = analysis.run_integrated(
            frame, object(), THRESHOLDS, method, scorer=builder, resolver=forbidden
        )
        assert result["n_removed"] == 3
        assert result["stop_reason"] == "max_removals_reached"
        assert counts["tie_invocations"] == 0


def test_integrated_tie_resolution_respects_maximum_and_unresolved_outcome():
    def tied(current, *args):
        return state(current, [3.0] * len(current), cost=8.0)

    def resolve(current, method):
        candidates = tuple(range(len(current.trace)))
        return analysis.Resolution(
            candidates,
            tuple(range(len(candidates))),
            0,
            tuple([np.nan] * len(candidates)),
            tuple([np.nan] * len(candidates)),
            tuple([""] * len(candidates)),
            tuple([""] * len(candidates)),
        )

    result, counts = analysis.run_integrated(
        trace(n=7), None, THRESHOLDS, "loo_top5", scorer=tied, resolver=resolve
    )
    assert result["n_removed"] == 3 and result["stop_reason"] == "max_removals_reached"
    assert counts["tie_invocations"] == counts["tie_resolved"] == 3
    result, counts = analysis.run_integrated(
        trace(), None, THRESHOLDS, "k1_local", scorer=tied
    )
    assert (
        result["abnormal_unresolved"] and result["stop_reason"] == "ambiguous_candidate"
    )
    assert result["n_removed"] == 0


@pytest.mark.parametrize("cost", [4.0, np.nan])
def test_global_stop_never_invokes_secondary_method(cost):
    def builder(frame, *args):
        return state(frame, cost=cost)

    def forbidden(*args):
        pytest.fail(
            "Secondary scores calculated for a globally normal/unscoreable trace"
        )

    result, counts = analysis.run_integrated(
        trace(), None, THRESHOLDS, "loo_top5", scorer=builder, resolver=forbidden
    )
    assert result["n_removed"] == 0 and counts["tie_invocations"] == 0


@pytest.mark.parametrize(
    "corrupted,classification,mixed",
    [
        ((), "no_true_candidate", False),
        ((0,), "single_true_candidate_in_tie", True),
        ((0, 1), "all_tied_candidates_true", False),
    ],
)
def test_oracle_classifications_and_method_denominators(
    corrupted, classification, mixed
):
    current = state(trace(corrupted=corrupted))
    scorer = analysis.TraceScorer(lambda frame, *args: state(frame, [1.0] * len(frame)))
    session = analysis.TieSession(
        scorer, None, "separation", False, analysis.EventBudget()
    )
    event, rows, decisions = analysis.tie_rows(current, 0, THRESHOLDS, meta(), session)
    assert event["tie_class"] == classification and event["mixed_tie"] is mixed
    assert len(rows) == 8 and len(decisions) == 4
    aggregate = analysis.EventAccumulator()
    for decision in decisions:
        aggregate.add(decision)
    pooled = aggregate.frame().query("aggregation == 'pooled'")
    assert (pooled["n_ambiguous_events"] == 1).all()
    assert (pooled["oracle_resolvable_fraction"] == bool(corrupted)).all()


def test_event_metric_precision_and_effective_correct_use_different_denominators():
    aggregate = analysis.EventAccumulator()
    for resolved, correct in ((True, True), (True, False), (False, False)):
        aggregate.add(
            {
                **meta(),
                "method": "loo_top5",
                "tie_size": 2,
                "oracle_resolvable": True,
                "all_tied_candidates_true": False,
                "mixed_tie": True,
                "method_scoreable": True,
                "unique_winner": resolved,
                "selected_is_true_corruption": correct,
            }
        )
    row = aggregate.frame().query("aggregation == 'pooled'").iloc[0]
    assert row["precision_among_resolved"] == 0.5
    assert row["effective_correct_resolution_rate"] == 1 / 3
    assert row["incorrect_resolution_rate"] == 1 / 3
    assert row["fraction_uniquely_resolved"] == 2 / 3


def test_secondary_scores_are_truth_blind():
    edges = [(0, 2, 5.0), (0, 3, 4.0), (1, 2, 3.0), (2, 3, 2.0), (2, 4, 1.0)]

    def scorer(frame, *args):
        return state(frame, [2.0] * len(frame))

    for method in analysis.METHODS:
        before = analysis.resolve_tie(state(trace(), edges=edges), method, scorer, None)
        after = analysis.resolve_tie(
            state(trace(corrupted=(0, 1, 2, 3, 4, 5)), edges=edges),
            method,
            scorer,
            None,
        )
        assert before.winner == after.winner
        np.testing.assert_array_equal(before.scores, after.scores)


def test_profile_limit_discards_partial_policy_sets_and_releases_states():
    live = weakref.WeakValueDictionary()

    def builder(frame, *args):
        result = state(
            frame, [3.0] * len(frame), edges=[(0, 2, 5.0), (0, 3, 4.0), (1, 2, 3.0)]
        )
        live[id(result)] = result
        return result

    budget = analysis.EventBudget(1)
    with pytest.raises(analysis.ProfileLimit):
        analysis.analyze_trace(
            trace(n=7), None, THRESHOLDS, meta(), budget, builder=builder
        )
    assert budget.used == 1
    assert not live


def test_trace_cache_is_bounded_and_freed():
    scorer = analysis.TraceScorer(stable_builder, max_states=2)
    live = weakref.WeakValueDictionary()
    for i in range(30):
        frame = trace().assign(Spot_ID=lambda data: data["Spot_ID"] + f"-{i}")
        result = scorer(frame, None, "separation", False)
        live[id(result)] = result
        del result
        assert len(scorer.cache) <= 2 and len(live) <= 2
    scorer.clear()
    assert not live


def test_streamed_details_and_zero_event_outputs_are_readable(tmp_path):
    cycle1 = C4._cycle3()._cycle(1)
    for rows in (0, 50):
        writer = analysis.DetailWriter(
            cycle1, tmp_path / f"events-{rows}.ecsv", analysis.TIE_FIELDS, max_rows=8
        )
        for i in range(rows):
            writer.write_rows(
                [{**meta(), "Trace_ID": f"t {i}", "iteration": 0, "tie_size": 2}]
            )
            assert len(writer.pending) < 8
            if i == 7:
                assert (
                    len(Table.read(writer.writer.temporary_path, format="ascii.ecsv"))
                    == 8
                )
        writer.finish()
        assert len(Table.read(writer.writer.path, format="ascii.ecsv")) == rows


def manifest_conditions():
    return [
        dict(
            id=f"c-{seed}-{p}-{k}-{d}",
            simulation_id=f"s-{seed}-{p}",
            n_barcodes=25,
            corrupted_barcodes_per_trace=k,
            seed=seed,
            detection_efficiency=p,
            displacement_um=d,
            directory="unused",
        )
        for p in (0.5, 0.8)
        for seed in (11, 12, 13)
        for k in (1, 2, 3)
        for d in (0.4, 0.8)
    ]


def test_narrow_condition_selection_and_fail_closed_missing_classes():
    conditions = manifest_conditions()
    extras = [
        dict(conditions[0], n_barcodes=50),
        dict(conditions[0], displacement_um=0.2),
        dict(conditions[0], corrupted_barcodes_per_trace=0),
    ]
    assert analysis.select_conditions({"conditions": conditions + extras}) == conditions
    with pytest.raises(ValueError, match="Missing cycle-5"):
        analysis.select_conditions(
            {
                "conditions": [
                    item
                    for item in conditions
                    if not (
                        item["corrupted_barcodes_per_trace"] == 3
                        and item["detection_efficiency"] == 0.8
                        and item["displacement_um"] == 0.8
                    )
                ]
            }
        )


def test_missing_historical_data_never_generates_a_benchmark(tmp_path):
    with pytest.raises(FileNotFoundError, match="No benchmark is generated"):
        analysis.benchmark(tmp_path / "missing", tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("stage", ["offline", "both"])
def test_end_to_end_in_memory_folds_determinism_and_bounded_reference_scope(
    tmp_path, monkeypatch, stage
):
    """Mock only the input readers; no synthetic sweep/ECSV inputs are generated."""
    import yaml

    conditions = manifest_conditions()
    simulations = [
        dict(id=f"s-{seed}-{p}", seed=seed, detection_efficiency=p, n_barcodes=25)
        for p in (0.5, 0.8)
        for seed in (11, 12, 13)
    ]
    root = tmp_path / "existing-input-interface"
    root.mkdir()
    (root / "sweep_manifest.yaml").write_text(
        yaml.safe_dump({"simulations": simulations, "conditions": conditions})
    )
    cycle3 = C4._cycle3()
    cycle1 = cycle3._cycle(1)

    def load_baseline(path, item):
        return trace().assign(
            simulation_id=item["id"],
            seed=item["seed"],
            detection_efficiency=item["detection_efficiency"],
        )

    def load_observations(module, path, item):
        return trace(
            corrupted=(0, 1) if item["corrupted_barcodes_per_trace"] >= 2 else (0,)
        ).assign(simulation_id=item["simulation_id"])

    fits = []
    live = weakref.WeakValueDictionary()
    prior_efficiency = [None]

    class Reference:
        pass

    def fit(module, baseline, minimum, mode):
        p = baseline["detection_efficiency"].iloc[0]
        if p != prior_efficiency[0]:
            assert not live, "Completed group's reference models escaped"
            prior_efficiency[0] = p
        reference = Reference()
        live[id(reference)] = reference
        fits.append((p, set(baseline["seed"])))
        return reference

    def build(frame, reference, mode, fallback):
        # Clean calibration and held-out controls remain globally normal.
        if not frame["is_corrupted"].any():
            return state(frame, [1.0] * len(frame), cost=4.0)
        edges = [(0, 2, 5.0), (0, 3, 4.0), (1, 2, 3.0), (2, 3, 2.0), (2, 4, 1.0)]
        return state(
            frame, [3.0, 3.0] + [0.0] * (len(frame) - 2), edges=edges, cost=8.0
        )

    monkeypatch.setattr(C4, "_cycle3", lambda: cycle3)
    monkeypatch.setattr(cycle3, "_cycle", lambda number: cycle1)
    monkeypatch.setattr(cycle1, "load_baseline", load_baseline)
    monkeypatch.setattr(cycle3, "_load_observations", load_observations)
    monkeypatch.setattr(C4, "_reference_with_coordinates", fit)
    monkeypatch.setattr(C4, "build_trace_scores", build)
    monkeypatch.setattr(analysis, "plot_results", lambda *args: None)
    for output in (tmp_path / "first", tmp_path / "second"):
        analysis.benchmark(root, output, minimum_observations=1, stage=stage)
        assert not live
    for filename in analysis.OUTPUTS[:6]:
        pd.testing.assert_frame_equal(
            Table.read(tmp_path / "first" / filename, format="ascii.ecsv").to_pandas(),
            Table.read(tmp_path / "second" / filename, format="ascii.ecsv").to_pandas(),
        )
    # Exact exclusion sets: every outer fit omits one seed; nested fits omit two.
    assert all(len(seeds) in (1, 2) for _, seeds in fits)
    for p in (0.5, 0.8):
        assert {frozenset(seeds) for efficiency, seeds in fits if efficiency == p} == {
            frozenset({11}),
            frozenset({12}),
            frozenset({13}),
            frozenset({11, 12}),
            frozenset({11, 13}),
            frozenset({12, 13}),
        }
    thresholds = Table.read(
        tmp_path / "first" / "cycle5_thresholds.ecsv", format="ascii.ecsv"
    ).to_pandas()
    for row in thresholds.itertuples():
        assert str(row.evaluation_seed) not in str(row.reference_seeds).split(",")
        for nested in row.nested_reference_seed_sets.split(";"):
            calibration, seeds = nested.split(":")
            assert calibration not in seeds.split(",")
            assert str(row.evaluation_seed) not in seeds.split(",")
    controls = Table.read(
        tmp_path / "first" / analysis.OUTPUTS[4], format="ascii.ecsv"
    ).to_pandas()
    assert (controls["mean_removals"] == 0).all()
    assert bool(len(controls)) == (stage == "both")
    runtime = Table.read(
        tmp_path / "first" / analysis.OUTPUTS[6], format="ascii.ecsv"
    ).to_pandas()
    assert (runtime.query("scope == 'group'")["reference_cache_after_clear"] == 0).all()
    assert runtime.iloc[-1]["peak_trace_states"] <= 64
    # Profiling must stop without emitting any partially compared trace.
    analysis.benchmark(
        root,
        tmp_path / "profile",
        minimum_observations=1,
        profile_ambiguous_events=1,
        stage=stage,
    )
    rows = Table.read(
        tmp_path / "profile" / "cycle5_trace_metrics.ecsv", format="ascii.ecsv"
    ).to_pandas()
    if len(rows):
        assert (
            rows.groupby(["simulation_id", "condition", "Trace_ID"]).size() == 4
        ).all()
    profile = Table.read(
        tmp_path / "profile" / analysis.OUTPUTS[6], format="ascii.ecsv"
    ).to_pandas()
    assert profile.iloc[-1]["n_unique_tie_states_calculated"] == 1
    assert profile.iloc[-1]["profile_limit_reached"]
    assert not live


def test_mixed_oracle_classification_with_multiple_true_candidates():
    current = state(trace(corrupted=(0, 1)), loo=[3.0, 3.0, 3.0, 0.0, 0.0, 0.0])
    session = analysis.TieSession(
        lambda frame, *args: state(frame, [1.0] * len(frame)),
        None,
        "separation",
        False,
        analysis.EventBudget(),
    )
    event, scores, _ = analysis.tie_rows(current, 1, THRESHOLDS, meta(), session)
    assert event["tie_class"] == "mixed_tie" and event["n_true_candidates"] == 2
    assert event["mixed_tie"] and not event["single_true_candidate_in_tie"]
    assert len(scores) == 12


@pytest.mark.parametrize("method", analysis.METHODS)
@pytest.mark.parametrize("reference_mode", ["separation", "barcode_pair"])
def test_real_empirical_non_tied_path_is_identical(method, reference_mode):
    frame = trace(corrupted=(0,))
    frame.loc[0, "x"] = 30.0
    reference = SimpleNamespace(
        separation={
            float(separation): np.linspace(
                separation / 10 - 0.5, separation / 10 + 0.5, 40
            )
            for separation in (10, 20, 30, 40, 50)
        },
        barcode_pair={
            C4._pair_key(
                frame.iloc[a]["Barcode"], frame.iloc[b]["Barcode"]
            ): np.linspace(b - a - 0.5, b - a + 0.5, 40)
            for a in range(len(frame))
            for b in range(a + 1, len(frame))
        },
        minimum_observations=1,
    )
    thresholds = dict(THRESHOLDS, global_threshold=3.0)
    expected = C4.run_policy(
        frame, reference, thresholds, "C_cap", reference_mode=reference_mode
    )[1]
    actual, counts = analysis.run_integrated(
        frame, reference, thresholds, method, reference_mode=reference_mode
    )
    assert counts["tie_invocations"] == 0
    assert actual == expected


def test_figure_outputs_use_compact_aggregates_and_handle_no_events(tmp_path):
    analysis.plot_results(
        tmp_path,
        analysis.EventAccumulator().frame(),
        analysis.IterativeAccumulator().frame(),
    )
    assert len(list(tmp_path.glob("*.png"))) == 7


def test_finite_secondary_winner_survives_input_permutation():
    current = state(
        edges=[(0, 2, 5.0), (0, 3, 4.0), (1, 2, 3.0), (2, 3, 2.0), (2, 4, 1.0)]
    )
    result = analysis.resolve_tie(current, "loo_top5", None, None)
    winner = current.trace.iloc[result.winner]["Spot_ID"]
    order = [1, 0, 5, 4, 3, 2]
    indices = {node: offset for offset, node in enumerate(order)}
    permuted = state(
        current.trace.iloc[order],
        loo=current.loo[order],
        edges=[(indices[a], indices[b], value) for a, b, value in current.edges],
    )
    result = analysis.resolve_tie(permuted, "loo_top5", None, None)
    assert permuted.trace.iloc[result.winner]["Spot_ID"] == winner
