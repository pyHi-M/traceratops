import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

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
        upper = (
            len(sample) - np.searchsorted(sample, observed, side="left") + 1
        ) / (len(sample) + 1)
        return -np.log10(min(1.0, 2.0 * min(lower, upper)))

    sorted_reference = np.sort(reference)
    for observed in (0.25, 0.5, 1.5, 3.0, 5.0):
        assert analysis.empirical_anomaly(
            sorted_reference, observed
        ) == legacy_empirical_anomaly(reference, observed)


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
    result = analysis.calibrate_thresholds(states)
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
