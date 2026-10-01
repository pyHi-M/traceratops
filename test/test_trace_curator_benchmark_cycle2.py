import gc
import importlib.util
import sys
import weakref
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
    assert not bool(attribution.iloc[0]["top1"])
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
