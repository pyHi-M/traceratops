import importlib.util
import gc
import sys
import weakref
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).parents[1] / "scripts" / "analyze_trace_curator_benchmark.py"
SPEC = importlib.util.spec_from_file_location("curator_analysis", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


def _reference_frame(outlier=False):
    rows = []
    for trace_number in range(30):
        for barcode in range(1, 5):
            x = float(barcode + trace_number / 100)
            if outlier and trace_number == 0 and barcode == 2:
                x = 100.0
            rows.append(
                (
                    f"t{trace_number}",
                    f"s{trace_number}-{barcode}",
                    barcode,
                    x,
                    0.0,
                    0.0,
                    False,
                )
            )
    return pd.DataFrame(
        rows, columns=["Trace_ID", "Spot_ID", "Barcode", "x", "y", "z", "is_corrupted"]
    )


def _evaluation(xs=(1.0, 2.0, 3.0, 4.0), labels=(False,) * 4):
    frame = pd.DataFrame(
        {
            "Trace_ID": ["eval"] * 4,
            "Spot_ID": [f"e{i}" for i in range(4)],
            "Barcode": [1, 2, 3, 4],
            "x": xs,
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": labels,
            "selected_for_corruption": labels,
            "injected_displacement_um": [1.0 if value else 0.0 for value in labels],
            "condition": "toy",
            "simulation_id": "sim",
            "replicate": 1,
            "seed": 1,
            "detection_efficiency": 1.0,
            "displacement": 1.0,
        }
    )
    return frame


def test_reference_excludes_corrupted_rows():
    frame = _reference_frame()
    frame.loc[
        (frame.Trace_ID == "t0") & (frame.Barcode == 2), ["x", "is_corrupted"]
    ] = [1000.0, True]
    model = analysis.fit_reference(frame, minimum_observations=1)
    assert model.separation[1].max() < 10


def test_reference_does_not_merge_reused_trace_ids_across_simulations():
    frame = pd.DataFrame(
        {
            "simulation_id": ["sim-a", "sim-a", "sim-b", "sim-b"],
            "Trace_ID": ["trace-1"] * 4,
            "Spot_ID": ["a1", "a2", "b1", "b2"],
            "Barcode": [1, 2, 1, 2],
            "x": [0.0, 1.0, 100.0, 102.0],
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": False,
        }
    )

    model = analysis.fit_reference(frame, minimum_observations=1)

    np.testing.assert_array_equal(model.separation[1], [1.0, 2.0])
    np.testing.assert_array_equal(model.barcode_pair[(1, 2)], [1.0, 2.0])


def test_scoring_does_not_merge_reused_trace_ids_across_simulations():
    first = _evaluation()
    second = _evaluation(xs=(10.0, 11.0, 12.0, 13.0))
    first["simulation_id"] = "sim-a"
    second["simulation_id"] = "sim-b"
    second["Spot_ID"] = [f"other-{value}" for value in second["Spot_ID"]]
    combined = pd.concat([first, second], ignore_index=True)

    scores = analysis.score_observations(
        combined, analysis.fit_reference(_reference_frame(), 1)
    )
    all_context = scores[
        (scores.model == "trace_splitter_like_residual") & (scores.context == "all")
    ]

    assert len(all_context) == 8
    assert (all_context.n_context_requested == 3).all()


def test_held_out_observation_does_not_enter_reference():
    reference = _reference_frame()
    evaluated = _evaluation(xs=(1, 200, 3, 4), labels=(False, True, False, False))
    before = analysis.fit_reference(reference, 1).separation[1].copy()
    analysis.score_observations(evaluated, analysis.fit_reference(reference, 1))
    np.testing.assert_array_equal(
        before, analysis.fit_reference(reference, 1).separation[1]
    )
    assert before.max() < 10


def test_empirical_tail_probability_known_distribution():
    sample = [1, 2, 3, 4, 5]
    assert analysis.empirical_tail_probability(sample, 3) == 1.0
    assert analysis.empirical_tail_probability(sample, 100) == 2 / 6
    assert analysis.empirical_tail_probability(
        [], 3
    ) != analysis.empirical_tail_probability([], 3)


def test_empirical_percentile_preserves_anomaly_direction():
    sample = [1, 2, 3, 4, 5]
    short_percentile, short_tail = analysis.empirical_distance_diagnostics(sample, 0)
    long_percentile, long_tail = analysis.empirical_distance_diagnostics(sample, 6)
    assert short_percentile < 0.5 < long_percentile
    assert short_tail == long_tail


def test_residual_reference_requires_exact_bin_unless_fallback_is_explicit():
    model = analysis.fit_reference(_reference_frame(), minimum_observations=1)
    model.residual.pop(2)
    assert model.residual_stats(2) == (None, False)
    stats, used_fallback = model.residual_stats(2, fallback_mode="nearest")
    assert stats is not None
    assert used_fallback


def test_context_selection_k_and_all():
    barcodes = [1, 3, 4, 7, 9, 12, 15, 20]
    assert analysis.select_context(barcodes, 9, "k1") == [3, 5]
    assert analysis.select_context(barcodes, 9, "k2") == [3, 2, 5, 6]
    assert analysis.select_context(barcodes, 9, "k3") == [3, 2, 1, 5, 6, 7]
    assert analysis.select_context(barcodes, 9, "all") == [0, 1, 2, 3, 5, 6, 7]


def test_missing_genomic_flank_is_recorded_and_available_side_is_used():
    result = analysis.score_observations(
        _evaluation(), analysis.fit_reference(_reference_frame(), 1)
    )
    first = result[
        (result.Spot_ID == "e0")
        & (result.model == "trace_splitter_like_residual")
        & (result.context == "k2")
    ].iloc[0]
    assert not first.has_left_flank
    assert first.has_right_flank
    assert first.n_context_used == 2
    assert first.n_context_requested == 2
    assert first.n_insufficient_reference == 0
    assert not first.both_flanks_available


def test_obvious_outlier_has_larger_multi_relationship_score():
    model = analysis.fit_reference(_reference_frame(), 1)
    central = analysis.score_observations(_evaluation(), model)
    outlier = analysis.score_observations(_evaluation(xs=(1, 100, 3, 4)), model)
    selector = lambda frame: frame[
        (frame.Spot_ID == "e1")
        & (frame.model == "empirical_separation_mean")
        & (frame.context == "all")
    ].score.iloc[0]
    assert selector(outlier) > selector(central)


def test_zero_context_and_insufficient_reference_are_explicit():
    single = _evaluation().iloc[:1].copy()
    model = analysis.fit_reference(_reference_frame(), minimum_observations=10_000)
    result = analysis.score_observations(single, model)
    assert result.score.isna().all()
    assert (result.score_status == "insufficient_reference").all()
    assert (result.n_context_used == 0).all()
    assert (result.n_insufficient_reference == 0).all()


def test_evaluation_labels_do_not_change_scores():
    model = analysis.fit_reference(_reference_frame(), 1)
    clean = _evaluation(labels=(False,) * 4)
    relabeled = _evaluation(labels=(True,) * 4)
    one = analysis.score_observations(clean, model).sort_values(
        ["Spot_ID", "model", "context"]
    )
    two = analysis.score_observations(relabeled, model).sort_values(
        ["Spot_ID", "model", "context"]
    )
    np.testing.assert_allclose(one.score, two.score, equal_nan=True)


def test_fixed_fpr_threshold_is_conservative_on_toy_data():
    values = np.arange(100, dtype=float)
    threshold = analysis.fixed_fpr_threshold(values, 0.05)
    assert np.mean(values > threshold) == 0.05
    assert analysis.fixed_fpr_threshold(values, 0.001) == 99


def test_scoring_is_deterministic():
    model = analysis.fit_reference(_reference_frame(), 1)
    first = analysis.score_observations(_evaluation(), model)
    second = analysis.score_observations(_evaluation(), model)
    pd.testing.assert_frame_equal(first, second)


def test_pairwise_output_includes_directional_and_insufficient_diagnostics():
    model = analysis.fit_reference(_reference_frame(), minimum_observations=1)
    model.barcode_pair.pop((1, 2))
    relationships = []
    analysis.score_observations(
        _evaluation(xs=(0, 20, 20, 4)), model, relationship_rows=relationships
    )
    pairwise = pd.DataFrame(relationships)
    assert {
        "empirical_percentile",
        "tail_probability",
        "pair_anomaly_score",
        "anomaly_direction",
        "reference_support",
        "relationship_status",
    }.issubset(pairwise.columns)
    missing = pairwise[
        (pairwise.reference_strategy == "barcode_pair")
        & (pairwise.Barcode == 1)
        & (pairwise.other_barcode == 2)
    ]
    assert (missing.relationship_status == "insufficient_reference").all()
    assert missing.empirical_percentile.isna().all()
    assert set(pairwise.anomaly_direction).issuperset({"short", "long"})


def test_run_releases_condition_score_frames_between_conditions(tmp_path, monkeypatch):
    root = tmp_path / "benchmark"
    output = tmp_path / "output"
    root.mkdir()
    seeds = [11, 22, 33]
    simulations = [
        {
            "id": f"simulation-{seed}",
            "seed": seed,
            "detection_efficiency": 1.0,
        }
        for seed in seeds
    ]
    conditions = [
        {
            "id": f"condition-{seed}-{index}",
            "simulation_id": f"simulation-{seed}",
            "seed": seed,
            "replicate_index": index,
            "detection_efficiency": 1.0,
            "displacement_um": float(index),
            "directory": f"condition-{seed}-{index}",
        }
        for seed in seeds
        for index in range(2)
    ]
    (root / "sweep_manifest.yaml").write_text(
        analysis.yaml.safe_dump(
            {"dry_run": False, "simulations": simulations, "conditions": conditions}
        )
    )

    def fake_baseline(path, metadata):
        return pd.DataFrame(
            {
                "seed": [metadata["seed"]],
                "detection_efficiency": [metadata["detection_efficiency"]],
            }
        )

    def fake_observations(path, metadata):
        return pd.DataFrame(
            {
                "condition": [metadata["id"]],
                "replicate": [metadata["replicate_index"]],
                "seed": [metadata["seed"]],
                "detection_efficiency": [metadata["detection_efficiency"]],
                "displacement": [metadata["displacement_um"]],
            }
        )

    condition_frames = []

    def fake_scores(frame, model, residual_fallback="none", relationship_rows=None):
        is_calibration = frame["condition"].iloc[0] == "calibration"
        result = pd.DataFrame(
            {
                "condition": frame["condition"],
                "replicate": frame["replicate"],
                "seed": frame["seed"],
                "detection_efficiency": frame["detection_efficiency"],
                "displacement": frame["displacement"],
                "model": "empirical_separation_mean",
                "context": "k3",
                "score": [0.0 if is_calibration else 1.0],
                "ground_truth_corrupted": [not is_calibration],
                "selected_for_corruption": [not is_calibration],
            }
        )
        if not is_calibration:
            condition_frames.append(weakref.ref(result))
            gc.collect()
            # At most the immediately preceding condition can still be referenced
            # while Python evaluates the next score_observations() call.
            assert sum(reference() is not None for reference in condition_frames) <= 2
        return result

    monkeypatch.setattr(analysis, "load_baseline", fake_baseline)
    monkeypatch.setattr(analysis, "load_observations", fake_observations)
    monkeypatch.setattr(analysis, "fit_reference", lambda frame, minimum: object())
    monkeypatch.setattr(
        analysis,
        "reference_summary",
        lambda model, fold: pd.DataFrame([{**fold, "n_observations": 1}]),
    )
    monkeypatch.setattr(analysis, "score_observations", fake_scores)
    monkeypatch.setattr(
        analysis,
        "plot_outputs",
        lambda representative, fixed, destination: pd.DataFrame(
            [{"representative_rows": len(representative)}]
        ),
    )

    analysis.run(root, output, minimum_observations=1)
    gc.collect()

    assert len(condition_frames) == len(conditions)
    assert all(reference() is None for reference in condition_frames)
    assert not (output / "localization_scores.ecsv").exists()
    assert not (output / "pairwise_relationship_scores.ecsv").exists()
