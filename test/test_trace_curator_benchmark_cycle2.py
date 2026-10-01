import gc
import importlib.util
import sys
import weakref
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy.table import Table

SCRIPT = Path(__file__).parents[1] / "scripts" / "analyze_trace_curator_benchmark_cycle2.py"
SPEC = importlib.util.spec_from_file_location("curator_cycle2", SCRIPT)
cycle2 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cycle2
SPEC.loader.exec_module(cycle2)


def toy_trace(xs=(0.0, 1.0, 2.0, 3.0)):
    return pd.DataFrame(
        {
            "Trace_ID": "trace",
            "Spot_ID": ["a", "b", "c", "d"],
            "Barcode": [17, 3, 41, 2],
            "Genomic_Position": [100.0, 130.0, 300.0, 900.0],
            "x": xs,
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": False,
        }
    )


def reference(minimum=1):
    parts = []
    for index in range(20):
        frame = toy_trace()
        frame["Trace_ID"] = f"r{index}"
        frame["x"] += index / 100
        parts.append(frame)
    return cycle2.fit_reference(pd.concat(parts), minimum)


def test_out_of_order_barcode_ids_are_preserved_as_identifiers():
    mapping = pd.DataFrame({"Barcode": [17, 3, 41], "Genomic_Position": [100, 130, 300]})
    joined = cycle2.attach_genomic_coordinates(pd.DataFrame({"Barcode": [3, 41, 17]}), mapping)
    assert joined.Genomic_Position.tolist() == [130, 300, 100]


def test_nearest_neighbours_follow_genomic_coordinates_not_barcode_ids():
    positions = toy_trace().Genomic_Position.to_numpy()
    assert cycle2.select_genomic_context(positions, 1, "k1") == [0, 2]


def test_nonuniform_spacing_is_used_for_separation():
    model = reference()
    assert {30.0, 170.0, 600.0}.issubset(model.separation)
    assert 24.0 not in model.separation  # barcode-number difference is irrelevant


def test_leave_one_out_is_cost_before_minus_cost_after():
    edges = [(0, 1, 1.0), (0, 2, 5.0), (1, 2, 2.0)]
    scores = cycle2.leave_one_out_improvements(3, edges, "mean")
    assert scores[0] == pytest.approx(np.mean([1, 5, 2]) - 2.0)


def test_obvious_bad_node_has_largest_leave_one_out_improvement():
    edges = [(0, 1, 0.1), (0, 2, 0.1), (0, 3, 0.1), (1, 2, 8), (1, 3, 9), (2, 3, 0.1)]
    scores = cycle2.leave_one_out_improvements(4, edges, "top3_mean")
    assert np.nanargmax(scores) == 1


def test_edge_concentration_identifies_multiple_incident_anomalies():
    edges = [(0, 1, 0.1), (0, 2, 0.1), (0, 3, 0.1), (1, 2, 6), (1, 3, 7), (2, 3, 0.1)]
    features = cycle2.edge_concentration_features(4, edges)
    assert features.concentration.idxmax() == 1
    assert features.loc[1, "count_strong"] == 2


def test_bridge_uses_true_left_and_right_genomic_flanks():
    trace = toy_trace()
    result = cycle2.bridge_score(trace, 1, reference())
    assert (result["left_index"], result["right_index"]) == (0, 2)
    assert (result["left_genomic_gap"], result["right_genomic_gap"]) == (30.0, 170.0)


def test_bridge_requires_both_flanks():
    result = cycle2.bridge_score(toy_trace(), 0, reference())
    assert result["status"] == "insufficient_bridge_context"
    assert np.isnan(result["score"])


def test_attribution_top1_top2_rank_and_margin():
    scores = pd.DataFrame(
        {
            "condition": "c", "seed": 1, "model": "m", "Trace_ID": "t",
            "Spot_ID": ["good1", "bad", "good2"], "score": [1.0, 2.0, 3.0],
            "is_corrupted": [False, True, False],
        }
    )
    attribution, _, ranks = cycle2.attribution_metrics(scores, 1.5)
    assert attribution.iloc[0]["rank"] == 2
    assert not bool(attribution.iloc[0]["top1_unique"])
    assert bool(attribution.iloc[0]["top2"])
    assert attribution.iloc[0]["score_margin"] == -1
    assert ranks.iloc[0]["corrupted_rank"] == 2


def test_collateral_calls_exclude_corrupted_localization():
    scores = pd.DataFrame(
        {"condition": "c", "seed": 1, "model": "m", "Trace_ID": "t", "score": [5, 4, 0], "is_corrupted": [True, False, False]}
    )
    _, collateral, _ = cycle2.attribution_metrics(scores, 1)
    assert collateral.iloc[0]["collateral_calls"] == 1


def test_ml_seed_split_rejects_all_conditions_from_evaluation_seed():
    training = pd.DataFrame({"seed": [1, 1, 2, 2], "condition": ["a", "b", "c", "d"]})
    cycle2.assert_seed_split(training, 3)
    with pytest.raises(ValueError, match="leaked"):
        cycle2.assert_seed_split(training, 2)


def test_labels_and_evaluation_metadata_are_not_model_features():
    assert not cycle2.FORBIDDEN_FEATURES.intersection(cycle2.ML_FEATURES)


def test_uniform_mapping_is_explicit_simulation_only():
    with pytest.raises(ValueError, match="--genomic-coordinates"):
        cycle2.load_genomic_coordinates(None, [1, 2], False)
    result = cycle2.load_genomic_coordinates(None, [2, 1], True)
    assert result.Genomic_Position.tolist() == [1.0, 2.0]


def test_uniform_mapping_uses_order_ranks_not_numeric_barcode_values():
    result = cycle2.load_genomic_coordinates(None, [41, 3, 17], True)
    assert result.Barcode.tolist() == [3, 17, 41]
    assert result.Genomic_Position.tolist() == [1.0, 2.0, 3.0]


def test_attribution_reports_unscoreable_corruption_as_unconditional_failure():
    scores = pd.DataFrame(
        {"condition": "c", "seed": 1, "model": "m", "Trace_ID": "t", "score": [np.nan, 4.0], "is_corrupted": [True, False]}
    )
    attribution, _, _ = cycle2.attribution_metrics(scores, 1)
    row = attribution.iloc[0]
    assert row.attribution_coverage == 0
    assert not bool(row.top1_unique)
    assert not bool(row.top2)
    assert np.isnan(row.top1_given_scoreable)


def test_attribution_distinguishes_tied_and_unique_top1():
    tied = pd.DataFrame(
        {"condition": "c", "seed": 1, "model": "m", "Trace_ID": "t", "score": [5.0, 5.0, 1.0], "is_corrupted": [True, False, False]}
    )
    attribution, _, _ = cycle2.attribution_metrics(tied, 10)
    row = attribution.iloc[0]
    assert bool(row.top1_including_ties)
    assert not bool(row.top1_unique)
    assert row.n_tied_best == 2


def test_performance_reports_score_coverage_and_effective_sensitivity():
    scores = pd.DataFrame(
        {
            "condition": "c", "replicate": 1, "seed": 1,
            "detection_efficiency": 1.0, "displacement": 1.0, "model": "m",
            "score": [2.0, np.nan, 2.0, np.nan],
            "is_corrupted": [True, True, False, False],
        }
    )
    _, fixed = cycle2.performance(scores, {("m", fpr): 1.0 for fpr in cycle2.FIXED_FPRS})
    row = fixed.iloc[0]
    assert row.n_corrupted_total == 2 and row.n_corrupted_scored == 1
    assert row.corrupted_coverage == 0.5 and row.genuine_coverage == 0.5
    assert row.sensitivity_given_scored == 1.0
    assert row.effective_sensitivity == 0.5


def test_long_scores_exposes_all_k1_baselines_and_barcode_pair_loo_variants():
    features = pd.DataFrame({
        "k1_separation_mean": [1], "k1_separation_maximum": [1],
        "k1_barcode_pair_mean": [1], "k1_barcode_pair_maximum": [1],
        "loo_separation_mean": [1], "loo_separation_trimmed_mean": [1],
        "loo_separation_top3_mean": [1], "loo_separation_strong_fraction": [1],
        "loo_barcode_pair_mean": [1], "loo_barcode_pair_trimmed_mean": [1],
        "loo_barcode_pair_top3_mean": [1], "loo_barcode_pair_strong_fraction": [1],
        "edge_separation_concentration": [1], "edge_barcode_pair_concentration": [1],
        "bridge_separation": [1],
    })
    names = set(cycle2.long_scores(features, include_ml=False).model)
    assert {"baseline_separation_k1_mean", "baseline_separation_k1_maximum", "baseline_barcode_pair_k1_mean", "baseline_barcode_pair_k1_maximum"}.issubset(names)
    assert {"loo_barcode_pair_mean", "loo_barcode_pair_trimmed_mean", "loo_barcode_pair_top3_mean", "loo_barcode_pair_strong_fraction"}.issubset(names)


def test_ml_training_features_crossfit_each_training_seed(monkeypatch):
    baselines = pd.DataFrame({"seed": [1, 2, 3, 4], "is_corrupted": False})
    items = [{"seed": seed} for seed in (1, 2, 3, 4)]
    mapping = pd.DataFrame({"Barcode": [1], "Genomic_Position": [1.0]})

    monkeypatch.setattr(cycle2, "fit_reference", lambda rows, minimum: set(rows.seed))
    monkeypatch.setattr(
        cycle2,
        "score_observations",
        lambda observations, model: pd.DataFrame({**{name: [0.0] for name in cycle2.ML_FEATURES}, "is_corrupted": [False], "seed": observations.seed}),
    )
    training = cycle2.build_crossfit_ml_training(
        items, {4}, baselines, 1,
        lambda item: pd.DataFrame({"Barcode": [1], "seed": [item["seed"]]}),
        mapping,
    )
    for row in training.itertuples():
        assert str(row.seed) not in row.feature_reference_seeds.split(",")
        assert "4" not in row.feature_reference_seeds.split(",")


def test_scoring_is_deterministic():
    first = cycle2.score_trace(toy_trace(), reference())
    second = cycle2.score_trace(toy_trace(), reference())
    pd.testing.assert_frame_equal(first, second)


def test_condition_score_frames_can_be_released_in_streaming_architecture():
    references = []
    for _ in range(3):
        frame = cycle2.score_trace(toy_trace(), reference())
        references.append(weakref.ref(frame))
        del frame
        gc.collect()
    assert all(item() is None for item in references)


def _write_ecsv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    Table.from_pandas(frame).write(path, format="ascii.ecsv", overwrite=True)


def test_run_end_to_end_tiny_four_seed_benchmark(tmp_path):
    root, output = tmp_path / "benchmark", tmp_path / "output"
    simulations, conditions = [], []
    for replicate, seed in enumerate((11, 22, 33, 44), start=1):
        simulation_id = f"sim-{seed}"
        simulations.append({"id": simulation_id, "seed": seed, "detection_efficiency": 1.0})
        baseline_rows = []
        for trace_index in range(2):
            for barcode in range(1, 5):
                baseline_rows.append({"Spot_ID": f"b-{trace_index}-{barcode}", "Trace_ID": f"t{trace_index}", "Barcode #": barcode, "x": float(barcode + trace_index * .1), "y": float(seed % 3) * .01, "z": 0.0})
        _write_ecsv(pd.DataFrame(baseline_rows), root / "simulations" / simulation_id / "simulated.ecsv")
        for displacement in (0.0, 1.0):
            condition_id = f"c-{seed}-{displacement:g}"
            directory = f"conditions/{condition_id}"
            main_rows, truth_rows = [], []
            for trace_index in range(2):
                for barcode in range(1, 5):
                    spot = f"s-{trace_index}-{barcode}"
                    selected = barcode == 2
                    corrupted = selected and displacement > 0
                    main_rows.append({"Spot_ID": spot, "Trace_ID": f"t{trace_index}", "Barcode #": barcode, "x": float(barcode + trace_index * .1 + (displacement * 5 if corrupted else 0)), "y": float(seed % 3) * .01, "z": 0.0})
                    truth_rows.append({"Spot_ID": spot, "Input_Trace_ID": f"t{trace_index}", "Barcode #": barcode, "selected_for_corruption": selected, "is_corrupted": corrupted, "injected_displacement_um": displacement if corrupted else 0.0})
            _write_ecsv(pd.DataFrame(main_rows), root / directory / "simulated.ecsv")
            _write_ecsv(pd.DataFrame(truth_rows), root / directory / "simulated.ground_truth.ecsv")
            conditions.append({"id": condition_id, "simulation_id": simulation_id, "seed": seed, "replicate_index": replicate, "detection_efficiency": 1.0, "displacement_um": displacement, "directory": directory})
    root.mkdir(exist_ok=True)
    (root / "sweep_manifest.yaml").write_text(yaml.safe_dump({"dry_run": False, "simulations": simulations, "conditions": conditions}))

    cycle2.run(root, output, assume_uniform_barcode_spacing=True, minimum_observations=1)

    expected = {
        "cycle2_model_performance.ecsv", "cycle2_fixed_fpr_performance.ecsv",
        "cycle2_attribution_performance.ecsv", "cycle2_collateral_calls.ecsv",
        "cycle2_rank_distribution.ecsv", "cycle2_reference_summary.ecsv",
        "cycle2_ml_feature_importance.ecsv", "cycle2_aggregate_performance.ecsv",
    }
    assert expected.issubset({path.name for path in output.iterdir()})
    fixed = Table.read(output / "cycle2_fixed_fpr_performance.ecsv", format="ascii.ecsv").to_pandas()
    assert {"corrupted_coverage", "sensitivity_given_scored", "effective_sensitivity"}.issubset(fixed.columns)
    assert "ml_logistic_regression" in set(fixed.model)
