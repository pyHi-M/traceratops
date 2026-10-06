import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "analyze_trace_curator_iterative_stopping.py"
)
SPEC = importlib.util.spec_from_file_location("cycle4", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


def trace(corrupted=(), barcodes=(90, 2, 40, 7, 18)):
    return pd.DataFrame(
        {
            "Trace_ID": "t",
            "Spot_ID": [f"p{i}" for i in range(len(barcodes))],
            "Barcode": barcodes,
            "Genomic_Position": np.arange(len(barcodes), dtype=float) * 10,
            "x": np.arange(len(barcodes), dtype=float),
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": [i in corrupted for i in range(len(barcodes))],
        }
    )


class ScriptedScorer:
    def __init__(self, costs, winners, drops=None):
        self.costs = costs
        self.winners = winners
        self.drops = drops or [2.0] * len(costs)
        self.calls = 0

    def __call__(self, frame, reference, mode, fallback):
        step = min(self.calls, len(self.costs) - 1)
        self.calls += 1
        cost = self.costs[step]
        loo = np.zeros(len(frame))
        winner_id = self.winners[step]
        matches = np.flatnonzero(frame["Spot_ID"].astype(str) == winner_id)
        if len(matches):
            loo[matches[0]] = self.drops[step]
        without = cost - loo
        return analysis.TraceScores(
            frame.reset_index(drop=True), (), cost, loo, without
        )


def tied_scorer(frame, reference, mode, fallback):
    loo = np.array([3.0, 3.0, *([0.0] * (len(frame) - 2))])
    return analysis.TraceScores(frame.reset_index(drop=True), (), 9.0, loo, 9.0 - loo)


THRESHOLDS = {
    "global_threshold": 5.0,
    "absolute_threshold": 1.0,
    "relative_threshold": 0.15,
}


def test_clean_trace_below_threshold_has_no_removal():
    audit, result = analysis.run_policy(
        trace(), object(), THRESHOLDS, "A_global", scorer=ScriptedScorer([4], ["p0"])
    )
    assert result["n_removed"] == 0
    assert audit.iloc[-1]["stop_reason"] == "global_cost_normal"


def test_one_bad_localization_is_removed_and_exact():
    _, result = analysis.run_policy(
        trace((2,)),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([8, 4], ["p2", "p0"]),
    )
    assert result["true_removed"] == 1 and result["exact_recovery"]


def test_two_bad_localizations_are_recomputed():
    scorer = ScriptedScorer([9, 7, 4], ["p1", "p3", "p0"])
    audit, result = analysis.run_policy(
        trace((1, 3)), object(), THRESHOLDS, "A_global", scorer=scorer
    )
    assert list(audit.loc[audit["removal_accepted"], "candidate_Barcode"]) == [2, 7]
    assert result["true_removed"] == 2 and scorer.calls == 3


def test_global_cost_can_remain_abnormal_after_correct_removal():
    audit, _ = analysis.run_policy(
        trace((1, 3)),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([9, 7, 4], ["p1", "p3", "p0"]),
    )
    assert audit.iloc[1]["global_abnormal_before"]


def test_abnormal_without_actionable_candidate_is_unresolved():
    audit, result = analysis.run_policy(
        trace(),
        object(),
        THRESHOLDS,
        "B_absolute",
        scorer=ScriptedScorer([8], ["p0"], [0.5]),
    )
    assert result["abnormal_unresolved"] and result["n_removed"] == 0
    assert audit.iloc[-1]["stop_reason"] == "candidate_not_actionable"


def test_tied_best_candidate_is_ambiguous_and_not_removed():
    audit, result = analysis.run_policy(
        trace((0,)), object(), THRESHOLDS, "A_global", scorer=tied_scorer
    )
    assert result["n_removed"] == 0
    assert result["abnormal_unresolved"]
    assert result["terminal_global_abnormal"]
    assert result["terminal_C"] == 9.0
    assert result["terminal_outcome"] == "abnormal_unresolved"
    assert audit.iloc[-1]["n_tied_best"] == 2
    assert not audit.iloc[-1]["candidate_unique"]
    assert audit.iloc[-1]["stop_reason"] == "ambiguous_candidate"


def test_maximum_removal_safeguard():
    audit, result = analysis.run_policy(
        trace((0, 1, 2)),
        object(),
        THRESHOLDS,
        "C_cap",
        max_removals=2,
        min_remaining=2,
        scorer=ScriptedScorer([10, 9, 8], ["p0", "p1", "p2"]),
    )
    assert result["n_removed"] == 2
    assert audit.iloc[-1]["stop_reason"] == "max_removals_reached"


def test_minimum_remaining_localization_safeguard():
    audit, result = analysis.run_policy(
        trace((0,)),
        object(),
        THRESHOLDS,
        "A_global",
        min_remaining=5,
        scorer=ScriptedScorer([10], ["p0"]),
    )
    assert result["n_removed"] == 0
    assert audit.iloc[-1]["stop_reason"] == "minimum_trace_size_reached"


def test_false_removal_and_unresolved_metrics():
    _, result = analysis.run_policy(
        trace((1, 3)),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([9, 4], ["p0", "p1"]),
    )
    assert result["false_removed"] == 1
    assert result["true_corruptions_remaining"] == 2
    assert result["unresolved_true_corruption"]
    assert not result["exact_recovery"]
    assert result["overcurated"]
    assert result["terminal_outcome"] == "overcurated"


def test_false_clean_terminal_classification():
    _, result = analysis.run_policy(
        trace((1,)),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([4], ["p1"]),
    )
    assert result["false_clean"]
    assert not result["terminal_global_abnormal"]
    assert result["terminal_outcome"] == "false_clean"


def test_clean_policy_metrics_capture_multiple_iterative_false_removals():
    _, result = analysis.run_policy(
        trace(),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([9, 8, 4], ["p0", "p1", "p2"]),
    )
    frame = pd.DataFrame(
        [
            {
                "dataset_type": "clean_evaluation",
                "policy": "A_global",
                "n_barcodes": 5,
                "K": 0,
                "detection_efficiency": 0.5,
                "displacement": 0.0,
                **result,
            }
        ]
    )
    summary = analysis.summarize(
        frame,
        [
            "dataset_type",
            "policy",
            "n_barcodes",
            "K",
            "detection_efficiency",
            "displacement",
        ],
    ).iloc[0]
    assert result["n_removed"] == 2
    assert summary["fraction_with_1plus_removal"] == 1
    assert summary["mean_removals"] == 2
    assert summary["fraction_with_2plus_removals"] == 1
    assert summary["maximum_removals"] == 2
    assert summary["fraction_terminal_global_abnormal"] == 0
    assert "global_cost_normal:1" in summary["stop_reason_distribution"]


def test_no_scoreable_candidate_and_insufficient_reference_reasons():
    def no_candidate(frame, reference, mode, fallback):
        values = np.full(len(frame), np.nan)
        return analysis.TraceScores(
            frame.reset_index(drop=True), (), 8.0, values, values
        )

    audit, _ = analysis.run_policy(
        trace(), object(), THRESHOLDS, "A_global", scorer=no_candidate
    )
    assert audit.iloc[-1]["stop_reason"] == "no_scoreable_candidate"
    audit, _ = analysis.run_policy(
        trace(),
        object(),
        THRESHOLDS,
        "A_global",
        scorer=ScriptedScorer([np.nan], ["p0"]),
    )
    assert audit.iloc[-1]["stop_reason"] == "insufficient_reference"


def test_calibration_excludes_evaluation_seed():
    reference, calibration, evaluation = analysis.fold_seed_sets([1, 2, 3], 3, 2)
    assert reference == {1} and calibration == {2} and evaluation == {3}


def test_genomic_positions_not_barcode_identity_control_separation():
    frame = trace()
    reference = SimpleNamespace(
        separation={10.0: np.linspace(0.5, 1.5, 20)},
        barcode_pair={},
        minimum_observations=1,
    )
    state = analysis.build_trace_scores(frame, reference)
    adjacent = {(left, right) for left, right, _ in state.edges}
    assert {(0, 1), (1, 2), (2, 3), (3, 4)}.issubset(adjacent)
    assert (0, 2) not in adjacent


def test_barcode_pair_reference_remains_supported():
    frame = trace(barcodes=("locus-z", "locus-a", "locus-q"))
    reference = SimpleNamespace(
        separation={},
        barcode_pair={("locus-a", "locus-z"): np.linspace(0.5, 1.5, 20)},
        minimum_observations=1,
    )
    assert len(analysis.build_trace_scores(frame, reference, "barcode_pair").edges) == 1


def test_presorted_reference_preserves_empirical_anomaly_values():
    reference = np.array([3.0, 0.5, 2.0, 4.5, 1.0])

    def legacy_empirical_anomaly(values, observed):
        sample = np.sort(np.asarray(values, float))
        lower = (np.searchsorted(sample, observed, side="right") + 1) / (
            len(sample) + 1
        )
        upper = (len(sample) - np.searchsorted(sample, observed, side="left") + 1) / (
            len(sample) + 1
        )
        return -np.log10(min(1.0, 2.0 * min(lower, upper)))

    sorted_reference = np.sort(reference)
    for observed in (0.25, 0.5, 1.5, 3.0, 5.0):
        np.testing.assert_allclose(
            analysis.empirical_anomaly(sorted_reference, observed),
            legacy_empirical_anomaly(reference, observed),
            rtol=1e-15,
        )


def test_reference_model_arrays_are_sorted_once():
    baseline = pd.concat(
        [
            trace().assign(simulation_id="s1"),
            trace().assign(
                simulation_id="s2",
                Trace_ID="t2",
                x=lambda frame: frame["x"][::-1].to_numpy(),
            ),
        ],
        ignore_index=True,
    )
    reference = analysis._reference_with_coordinates(None, baseline, 1, "separation")
    for samples in (*reference.separation.values(), *reference.barcode_pair.values()):
        assert np.all(samples[:-1] <= samples[1:])


def test_execution_is_deterministic():
    kwargs = dict(
        trace=trace((2,)),
        reference=object(),
        thresholds=THRESHOLDS,
        policy="A_global",
        scorer=ScriptedScorer([8, 4], ["p2", "p0"]),
    )
    first = analysis.run_policy(**kwargs)[1]
    kwargs["scorer"] = ScriptedScorer([8, 4], ["p2", "p0"])
    assert first == analysis.run_policy(**kwargs)[1]


def test_clean_calibration_produces_all_operating_points():
    frame = trace()
    states = [
        analysis.TraceScores(frame, (), float(i), np.arange(5), float(i) - np.arange(5))
        for i in range(1, 101)
    ]
    result = analysis.calibrate_thresholds(
        [analysis.CompactTraceStats.from_scores(state) for state in states]
    )
    assert set(result["trace_fpr"]) == set(analysis.FIXED_FPRS)
    assert set(result["n_clean_traces"]) == {100}
    assert {
        "fraction_global_abnormal",
        "fraction_global_and_absolute_actionable",
        "fraction_global_and_relative_actionable",
    }.issubset(result.columns)


def test_shared_score_cache_pattern_avoids_duplicate_reference_scoring():
    raw = ScriptedScorer([4], ["p0"])
    cache = {}

    def shared(frame, reference, mode, fallback):
        key = tuple(sorted(frame["Spot_ID"]))
        if key not in cache:
            cache[key] = raw(frame, reference, mode, fallback)
        return cache[key]

    for policy in analysis.POLICIES:
        analysis.run_policy(trace(), object(), THRESHOLDS, policy, scorer=shared)
    assert raw.calls == 1


def legacy_calibration(states):
    """Independent oracle for the former full-state calibration formulas."""
    samples = {"global": [], "absolute": [], "relative": []}
    for state in states:
        finite = np.flatnonzero(np.isfinite(state.loo))
        samples["global"].append(state.cost)
        if len(finite):
            best = finite[np.argmax(state.loo[finite])]
            samples["absolute"].append(state.loo[best])
            samples["relative"].append(state.loo[best] / max(state.cost, 1e-12))
    rows = []
    for fpr in analysis.FIXED_FPRS:
        thresholds = {
            f"{key}_threshold": analysis.empirical_threshold(values, fpr)
            for key, values in samples.items()
        }
        flags = [[], [], []]
        for state in states:
            finite = state.loo[np.isfinite(state.loo)]
            best = np.max(finite) if len(finite) else np.nan
            relative = best / max(state.cost, 1e-12) if np.isfinite(best) else np.nan
            abnormal = bool(
                np.isfinite(state.cost) and state.cost > thresholds["global_threshold"]
            )
            flags[0].append(abnormal)
            flags[1].append(abnormal and best > thresholds["absolute_threshold"])
            flags[2].append(abnormal and relative > thresholds["relative_threshold"])
        rows.append(
            {
                "trace_fpr": fpr,
                **thresholds,
                "observed_clean_global_fpr": np.mean(
                    np.asarray(samples["global"]) > thresholds["global_threshold"]
                ),
                "n_clean_traces": len(states),
                **dict(zip(analysis._GATE_KEYS, map(np.mean, flags))),
            }
        )
    return pd.DataFrame(rows)


def test_compact_calibration_matches_full_state_oracle_and_releases_states():
    import weakref

    rng = np.random.default_rng(13)
    states = [
        analysis.TraceScores(trace(), (), float(cost), loo, cost - loo)
        for cost, loo in zip(rng.uniform(0, 10, 200), rng.normal(size=(200, 5)))
    ]
    states.extend(
        [
            analysis.TraceScores(
                trace(), (), np.nan, np.full(5, np.nan), np.full(5, np.nan)
            ),
            analysis.TraceScores(trace(), (), 0.0, np.zeros(5), np.zeros(5)),
            analysis.TraceScores(
                trace(), (), 2.0, np.full(5, np.nan), np.full(5, np.nan)
            ),
        ]
    )
    expected = legacy_calibration(states)
    refs = [weakref.ref(state) for state in states]
    compact = [analysis.CompactTraceStats.from_scores(state) for state in states]
    del states
    assert all(ref() is None for ref in refs)
    pd.testing.assert_frame_equal(analysis.calibrate_thresholds(compact), expected)


def test_streamed_condition_summaries_match_summarize():
    rng = np.random.default_rng(97)
    rows = []
    accumulator = analysis.ConditionAccumulator()
    for i in range(211):
        removed, true, false = (int(x) for x in rng.integers(0, 4, 3))
        row = {
            "dataset_type": "clean_evaluation" if i % 3 else "corrupted_evaluation",
            "policy": analysis.POLICIES[i % 4],
            "n_barcodes": 5,
            "K": i % 3,
            "detection_efficiency": 0.5,
            "displacement": 0.0,
            "true_removed": true,
            "n_true_corruptions": true + 2,
            "n_removed": removed,
            "false_removed": false,
            "exact_recovery": i % 2 == 0,
            "true_corruptions_remaining": 2,
            "unresolved_true_corruption": True,
            "abnormal_unresolved": i % 5 == 0,
            "terminal_global_abnormal": i % 7 == 0,
            "n_iterations": removed,
            "stop_reason": ("normal", "ambiguous", "cap")[i % 3],
        }
        rows.append(row)
        accumulator.add(row)
    expected = analysis.summarize(pd.DataFrame(rows), analysis.SUMMARY_GROUP)
    actual = accumulator.frame()[expected.columns]
    pd.testing.assert_frame_equal(actual, expected)
    # Repetition changes counts, not the number of retained condition buckets.
    size = len(accumulator.groups)
    for _ in range(50):
        for row in rows:
            accumulator.add(row)
    assert len(accumulator.groups) == size


def test_audit_and_trace_writers_flush_bounded_chunks(tmp_path):
    from astropy.table import Table

    cycle1 = analysis._cycle3()._cycle(1)
    for name in (analysis.OUTPUTS[0], analysis.OUTPUTS[2]):
        underlying = cycle1._EcsvChunkWriter(tmp_path / name)
        writer = analysis.BufferedEcsvWriter(underlying, max_rows=8)
        for i in range(25):
            writer.write(pd.DataFrame([{"Trace_ID": f"trace {i}", "score": float(i)}]))
            assert len(writer.pending) < 8
            if i == 7:
                assert underlying.temporary_path.exists()
                assert (
                    len(Table.read(underlying.temporary_path, format="ascii.ecsv")) == 8
                )
        writer.finish()
        result = Table.read(tmp_path / name, format="ascii.ecsv").to_pandas()
        assert len(result) == 25
        assert result["Trace_ID"].iloc[-1] == "trace 24"
        assert not writer.pending


def test_unused_reference_mode_is_not_built(monkeypatch):
    baseline = trace().assign(simulation_id="s")

    def forbidden_pair(*args):
        raise AssertionError("Unused barcode pair reference was constructed")

    monkeypatch.setattr(analysis, "_pair_key", forbidden_pair)
    reference = analysis._reference_with_coordinates(None, baseline, 1, "separation")
    assert reference.separation and not reference.barcode_pair
    monkeypatch.undo()
    reference = analysis._reference_with_coordinates(None, baseline, 1, "barcode_pair")
    assert reference.barcode_pair and not reference.separation
    assert all(isinstance(key, tuple) for key in reference.barcode_pair)


def make_tiny_benchmark(root, n_traces=3, efficiencies=(0.5,), n_barcodes=5):
    """Actual ECSV simulator inputs; also usable for local profiling."""
    import yaml
    from astropy.table import Table

    simulations, conditions = [], []
    for efficiency in efficiencies:
        for seed in (11, 12, 13):
            simulation_id = f"s-{efficiency}-{seed}"
            simulation = {
                "id": simulation_id,
                "seed": seed,
                "detection_efficiency": efficiency,
                "n_barcodes": n_barcodes,
            }
            simulations.append(simulation)
            frames = []
            for i in range(n_traces):
                frame = trace(barcodes=tuple(range(1, n_barcodes + 1))).rename(
                    columns={"Barcode": "Barcode #"}
                )
                frame["Trace_ID"] = f"t-{i}"
                frame["Spot_ID"] = [f"t-{i}-p-{j}" for j in range(n_barcodes)]
                frame["x"] *= 0.9 + 0.03 * seed + 0.02 * (i % 7)
                frames.append(frame.drop(columns="is_corrupted"))
            clean = pd.concat(frames, ignore_index=True)
            simulation_dir = root / "simulations" / simulation_id
            simulation_dir.mkdir(parents=True)
            Table.from_pandas(clean).write(
                simulation_dir / "simulated.ecsv", format="ascii.ecsv"
            )
            condition_id = f"c-{efficiency}-{seed}"
            condition_dir = root / condition_id
            condition_dir.mkdir()
            observations = clean.copy()
            corrupted = observations["Spot_ID"].str.endswith("p-2")
            observations.loc[corrupted, "y"] += 5.0
            Table.from_pandas(observations).write(
                condition_dir / "simulated.ecsv", format="ascii.ecsv"
            )
            truth = observations[["Spot_ID", "Trace_ID", "Barcode #"]].rename(
                columns={"Trace_ID": "Input_Trace_ID"}
            )
            truth["selected_for_corruption"] = True
            truth["is_corrupted"] = corrupted
            truth["injected_displacement_um"] = np.where(corrupted, 5.0, 0.0)
            Table.from_pandas(truth).write(
                condition_dir / "simulated.ground_truth.ecsv", format="ascii.ecsv"
            )
            conditions.append(
                {
                    "id": condition_id,
                    "directory": condition_id,
                    "simulation_id": simulation_id,
                    "seed": seed,
                    "replicate_index": seed - 11,
                    "n_barcodes": n_barcodes,
                    "corrupted_barcodes_per_trace": 1,
                    "detection_efficiency": efficiency,
                    "displacement_um": 5.0,
                }
            )
    (root / "sweep_manifest.yaml").write_text(
        yaml.safe_dump(
            {
                "simulations": simulations,
                "conditions": conditions,
            }
        )
    )


@pytest.mark.parametrize("reference_mode", ["separation", "barcode_pair"])
def test_tiny_benchmark_deterministic_and_streamed_outputs(tmp_path, reference_mode):
    from astropy.table import Table

    root = tmp_path / "inputs"
    make_tiny_benchmark(root)
    first, second = tmp_path / "first", tmp_path / "second"
    analysis.benchmark(
        root, first, minimum_observations=1, reference_mode=reference_mode
    )
    analysis.benchmark(
        root, second, minimum_observations=1, reference_mode=reference_mode
    )
    for name in analysis.OUTPUTS[:4]:
        pd.testing.assert_frame_equal(
            Table.read(first / name, format="ascii.ecsv").to_pandas(),
            Table.read(second / name, format="ascii.ecsv").to_pandas(),
        )
    metrics = Table.read(first / analysis.OUTPUTS[0], format="ascii.ecsv").to_pandas()
    assert len(metrics) == 3 * 3 * 2 * 4
    assert set(metrics["dataset_type"]) == {"clean_evaluation", "corrupted_evaluation"}
    expected = analysis.summarize(metrics, analysis.SUMMARY_GROUP)
    actual = Table.read(first / analysis.OUTPUTS[1], format="ascii.ecsv").to_pandas()
    pd.testing.assert_frame_equal(actual[expected.columns], expected)
    # Standalone replotting streams the audit too.
    analysis.plot_results(first)


def test_many_traces_keep_states_audits_and_reference_cache_bounded(
    tmp_path, monkeypatch
):
    import weakref

    from astropy.table import Table

    root = tmp_path / "inputs"
    make_tiny_benchmark(root, n_traces=50, efficiencies=(0.5, 0.8))
    live_states, live_audits, live_references = (
        weakref.WeakValueDictionary() for _ in range(3)
    )
    peaks = {"states": 0, "audits": 0}
    raw_build, raw_policy, raw_fit = (
        analysis.build_trace_scores,
        analysis.run_policy,
        analysis._reference_with_coordinates,
    )
    previous_group = None
    reference_groups = []

    def build(*args, **kwargs):
        state = raw_build(*args, **kwargs)
        live_states[id(state)] = state
        peaks["states"] = max(peaks["states"], len(live_states))
        return state

    def policy(*args, **kwargs):
        audit, metrics = raw_policy(*args, **kwargs)
        live_audits[id(audit)] = audit
        peaks["audits"] = max(peaks["audits"], len(live_audits))
        return audit, metrics

    class Reference:
        pass

    def fit(*args, **kwargs):
        nonlocal previous_group
        group = analysis._PROGRESS_CONTEXT.split(" eval=")[0]
        if group != previous_group:
            assert not live_references, "Completed group's references escaped the cache"
            reference_groups.append(group)
            previous_group = group
        reference = Reference()
        reference.__dict__.update(raw_fit(*args, **kwargs).__dict__)
        live_references[id(reference)] = reference
        return reference

    raw_calibrate = analysis.calibrate_thresholds

    def calibrate(stats):
        assert all(isinstance(item, analysis.CompactTraceStats) for item in stats)
        assert not live_states, "Calibration retained full scoring states"
        return raw_calibrate(stats)

    monkeypatch.setattr(analysis, "build_trace_scores", build)
    monkeypatch.setattr(analysis, "run_policy", policy)
    monkeypatch.setattr(analysis, "_reference_with_coordinates", fit)
    monkeypatch.setattr(analysis, "calibrate_thresholds", calibrate)
    monkeypatch.setattr(analysis, "plot_results", lambda *args: None)
    analysis.benchmark(root, tmp_path / "out", minimum_observations=1)
    assert peaks["states"] <= 12  # Four policy paths, at most three states each.
    assert peaks["audits"] <= 2
    assert not live_states and not live_audits and not live_references
    assert len(reference_groups) == 2
    runtime = Table.read(
        tmp_path / "out" / analysis.OUTPUTS[4], format="ascii.ecsv"
    ).to_pandas()
    groups = runtime[runtime["scope"] == "group"]
    assert (groups["reference_cache_size_after_clear"] == 0).all()
    assert (groups["rss_peak_mb"] >= groups["rss_start_mb"]).all()
    assert runtime.iloc[-1]["profile_traces"] == 600


def test_profile_calibration_limit_is_lazy_and_effective(tmp_path, monkeypatch):
    import inspect

    source = inspect.getsource(analysis.benchmark)
    compact_source = "".join(source.split())
    assert "islice(calibration_groups,profile_calibration_traces)" in compact_source
    assert "list(calibration_groups)" not in compact_source
    root = tmp_path / "inputs"
    make_tiny_benchmark(root, n_traces=10)
    lengths = []
    raw = analysis.calibrate_thresholds

    def calibrate(stats):
        lengths.append(len(stats))
        return raw(stats)

    monkeypatch.setattr(analysis, "calibrate_thresholds", calibrate)
    monkeypatch.setattr(analysis, "plot_results", lambda *args: None)
    analysis.benchmark(
        root,
        tmp_path / "out",
        minimum_observations=1,
        profile_traces=1,
        profile_calibration_traces=2,
    )
    assert lengths == [4]  # Two calibration seeds, first two traces from each.


def test_incremental_cost_trajectories_match_full_audit():
    accumulator = analysis.CostTrajectories()
    chunks = []
    for i in range(9):
        chunk = pd.DataFrame(
            {
                "dataset_type": ["corrupted_evaluation"] * 4,
                "policy": ["A_global", "B_relative"] * 2,
                "iteration": [0, 0, 1, 1],
                "C_before": [i + 1.0, np.nan, 2.0, 3.0],
            }
        )
        chunks.append(chunk)
        accumulator.add(chunk)
    expected = (
        pd.concat(chunks)
        .groupby(["policy", "iteration"], sort=False)["C_before"]
        .mean()
    )
    actual = accumulator.frame().set_index(["policy", "iteration"])["C_before"]
    pd.testing.assert_series_equal(actual, expected)


@pytest.mark.parametrize("barcode_numeric", [True, False])
def test_audit_chunks_preserve_missing_candidate_types(tmp_path, barcode_numeric):
    from astropy.table import Table

    writer = analysis.BufferedEcsvWriter(
        analysis._cycle3()._cycle(1)._EcsvChunkWriter(tmp_path / "audit.ecsv"),
        max_rows=1,
    )
    writer.barcode_numeric = barcode_numeric
    for value in (2 if barcode_numeric else "locus q", np.nan):
        writer.write(
            pd.DataFrame(
                [
                    {
                        "candidate_Barcode": value,
                        "candidate_rank": 1.0 if pd.notna(value) else np.nan,
                    }
                ]
            )
        )
    writer.finish()
    result = Table.read(tmp_path / "audit.ecsv", format="ascii.ecsv").to_pandas()
    assert len(result) == 2
    assert result.iloc[0]["candidate_Barcode"] == (2 if barcode_numeric else "locus q")
    assert pd.isna(result.iloc[1]["candidate_rank"])
    assert pd.isna(result.iloc[1]["candidate_Barcode"])


@pytest.mark.parametrize("size_source", ["simulation", "condition", "observed"])
def test_mixed_barcode_populations_match_cycle3_membership(
    tmp_path, monkeypatch, size_source
):
    import re

    import yaml

    root = tmp_path / "inputs"
    simulations, conditions = [], []
    for n in (5, 7):
        population = root / f"population-{n}"
        make_tiny_benchmark(population, n_traces=2, n_barcodes=n)
        manifest = yaml.safe_load((population / "sweep_manifest.yaml").read_text())
        for item in manifest["simulations"]:
            original_id = item["id"]
            item["id"] = f"B{n}-{original_id}"
            (root / "simulations").mkdir(exist_ok=True)
            (population / "simulations" / original_id).rename(
                root / "simulations" / item["id"]
            )
            if size_source != "simulation":
                item.pop("n_barcodes")
            simulations.append(item)
        for item in manifest["conditions"]:
            item["id"] = f"B{n}-{item['id']}"
            item["directory"] = str(Path(f"population-{n}") / item["directory"])
            item["simulation_id"] = f"B{n}-{item['simulation_id']}"
            if size_source == "observed":
                # No baseline simulation ID in the condition-derived map:
                # exercise cycle-3's last-resort observed barcode fallback.
                item["simulation_id"] = f"observations-{item['simulation_id']}"
            conditions.append(item)
    (root / "sweep_manifest.yaml").write_text(
        yaml.safe_dump({"simulations": simulations, "conditions": conditions})
    )

    cycle3 = analysis._cycle3()
    cycle1 = cycle3._cycle(1)
    simulation_sizes = {
        str(item["simulation_id"]): cycle3._condition_metadata(item)[0]
        for item in conditions
    }
    # Independent oracle: the exact population assignment/selection from
    # cycle-3 run(), including the loader's absence of n_barcodes provenance.
    cycle3_baselines = []
    for item in simulations:
        frame = cycle1.load_baseline(root / "simulations" / item["id"], item)
        assert "n_barcodes" not in frame
        frame["n_barcodes"] = int(
            item.get(
                "n_barcodes",
                simulation_sizes.get(str(item["id"]), frame["Barcode"].max()),
            )
        )
        cycle3_baselines.append(frame)
    expected = pd.concat(cycle3_baselines, ignore_index=True)
    loaded, fitted, calibrated = [], [], []
    raw_load = cycle1.load_baseline
    raw_fit = analysis._reference_with_coordinates
    raw_score = analysis.build_trace_scores

    def context():
        n = int(re.search(r"group B=(\d+)", analysis._PROGRESS_CONTEXT)[1])
        evaluation = re.search(r"eval=(\d+)", analysis._PROGRESS_CONTEXT)
        calibration = re.search(r"calibration=(\d+)", analysis._PROGRESS_CONTEXT)
        return (
            n,
            int(evaluation[1]) if evaluation else None,
            int(calibration[1]) if calibration else None,
        )

    def load(path, item):
        n, _, _ = context()
        if size_source != "observed":
            count = int(item.get("n_barcodes", simulation_sizes.get(item["id"])))
            assert count == n, "Nonmatching ECSV was loaded instead of prefiltered"
        loaded.append((n, item["id"]))
        return raw_load(path, item)

    def fit(module, baseline, minimum, mode):
        n, evaluation, calibration = context()
        seeds, _, _ = cycle3.fold_seed_sets((11, 12, 13), evaluation, calibration)
        population = expected[
            (expected["n_barcodes"] == n)
            & np.isclose(expected["detection_efficiency"], 0.5)
            & expected["seed"].isin(seeds)
        ]
        columns = ["simulation_id", "Trace_ID", "Spot_ID", "n_barcodes"]
        pd.testing.assert_frame_equal(
            baseline[columns].sort_values(columns).reset_index(drop=True),
            population[columns].sort_values(columns).reset_index(drop=True),
        )
        fitted.append((n, evaluation, calibration))
        return raw_fit(module, baseline, minimum, mode)

    def score(frame, reference, mode, fallback):
        n, _, calibration = context()
        if calibration is not None:
            population = expected[
                (expected["n_barcodes"] == n) & (expected["seed"] == calibration)
            ]
            assert set(frame["simulation_id"]).issubset(
                set(population["simulation_id"])
            )
            assert set(frame["Spot_ID"]).issubset(set(population["Spot_ID"]))
            assert (frame["n_barcodes"] == n).all()
            calibrated.append((n, calibration))
        return raw_score(frame, reference, mode, fallback)

    monkeypatch.setattr(analysis, "_cycle3", lambda: cycle3)
    monkeypatch.setattr(cycle3, "_cycle", lambda number: cycle1)
    monkeypatch.setattr(cycle1, "load_baseline", load)
    monkeypatch.setattr(analysis, "_reference_with_coordinates", fit)
    monkeypatch.setattr(analysis, "build_trace_scores", score)
    monkeypatch.setattr(analysis, "plot_results", lambda *args: None)
    analysis.benchmark(root, tmp_path / "out", minimum_observations=1)
    assert {n for n, _, _ in fitted} == {5, 7}
    assert len(fitted) == 12  # Six unique outer/nested reference models per group.
    assert (
        len(calibrated) == 24
    )  # Two traces x two calibration seeds x three folds x two sizes.
    assert len(loaded) == (12 if size_source == "observed" else 6)
