import numpy as np

from traceratops.trace_simulator import (
    SimulationConfig,
    random_unit_vectors,
    simulate_traces,
)


def _config(**changes):
    values = dict(
        n_traces=12,
        n_barcodes=10,
        detection_efficiency=1.0,
        mislocalization_fraction=1.0,
        mislocalization_displacement=0.2,
        mislocalization_edge_margin=2,
        seed=41,
    )
    values.update(changes)
    return SimulationConfig(**values)


def test_fraction_zero_has_no_corruption():
    _, truth = simulate_traces(_config(mislocalization_fraction=0))
    assert not np.any(truth["Mislocalization_Selected_Trace"])
    assert not np.any(truth["Is_Mislocalized"])


def test_zero_displacement_records_target_and_preserves_coordinates():
    observed, truth = simulate_traces(_config(mislocalization_displacement=0))
    targets = truth[truth["Is_Mislocalized"]]
    assert len(targets) == 12
    for axis in "xyz":
        np.testing.assert_array_equal(
            targets[f"PreCorruption_Observed_{axis}"],
            targets[f"PostCorruption_Observed_{axis}"],
        )
    assert np.all(targets["Realized_Displacement_um"] == 0)
    assert np.all(targets["selected_for_corruption"])
    assert not np.any(targets["is_corrupted"])
    assert len(observed) == 120


def test_one_spot_moves_without_changing_identity_or_multiplicity():
    observed, truth = simulate_traces(_config())
    detected_truth = truth[truth["Detection_Source"] == "polymer"]
    for trace_id in set(observed["Trace_ID"]):
        rows = detected_truth[detected_truth["Input_Trace_ID"] == trace_id]
        assert np.sum(rows["Is_Mislocalized"]) == 1
        changed = np.linalg.norm(
            np.column_stack([rows[f"PostCorruption_Observed_{a}"] for a in "xyz"])
            - np.column_stack([rows[f"PreCorruption_Observed_{a}"] for a in "xyz"]),
            axis=1,
        )
        assert np.sum(changed > 0) == 1
        np.testing.assert_allclose(changed[rows["Is_Mislocalized"]], 0.2, atol=1e-12)
        output = observed[observed["Trace_ID"] == trace_id]
        assert len(np.unique(output["Barcode #"])) == len(output)
        assert set(output["Spot_ID"]) == set(rows["Spot_ID"])
        assert set(output["Barcode #"]) == set(rows["Barcode"])
    assert np.all(
        (truth[truth["Is_Mislocalized"]]["Barcode"] >= 3)
        & (truth[truth["Is_Mislocalized"]]["Barcode"] <= 8)
    )


def test_reproducible_and_latent_truth_independent_of_corruption():
    observed_a, truth_a = simulate_traces(_config())
    observed_b, truth_b = simulate_traces(_config())
    assert observed_a.as_array().tobytes() == observed_b.as_array().tobytes()
    assert truth_a.as_array().tobytes() == truth_b.as_array().tobytes()
    _, control = simulate_traces(_config(mislocalization_fraction=0))
    for axis in "xyz":
        np.testing.assert_array_equal(
            truth_a[f"Latent_{axis}"], control[f"Latent_{axis}"]
        )
        np.testing.assert_array_equal(
            truth_a[f"PreCorruption_Observed_{axis}"],
            control[f"PreCorruption_Observed_{axis}"],
        )


def test_paired_displacements_reuse_targets_and_directions():
    _, low = simulate_traces(_config(mislocalization_displacement=0.1))
    _, high = simulate_traces(_config(mislocalization_displacement=0.4))
    low = low[low["selected_for_corruption"]]
    high = high[high["selected_for_corruption"]]
    assert list(low["Spot_ID"]) == list(high["Spot_ID"])
    low_delta = np.column_stack(
        [low[f"Observed_{a}"] - low[f"Uncorrupted_Observed_{a}"] for a in "xyz"]
    )
    high_delta = np.column_stack(
        [high[f"Observed_{a}"] - high[f"Uncorrupted_Observed_{a}"] for a in "xyz"]
    )
    np.testing.assert_allclose(high_delta, 4 * low_delta, atol=1e-12)


def test_legacy_benchmark_ground_truth_schema_is_present():
    _, truth = simulate_traces(_config(n_traces=1))
    expected = {
        "Spot_ID",
        "Input_Trace_ID",
        "GroundTruth_Group_ID",
        "GroundTruth_Polymer_ID",
        "Detection_Source",
        "Barcode #",
        "Occurrence_Index",
        "True_x",
        "True_y",
        "True_z",
        "Observed_x",
        "Observed_y",
        "Observed_z",
        "Uncorrupted_Observed_x",
        "Uncorrupted_Observed_y",
        "Uncorrupted_Observed_z",
        "selected_for_corruption",
        "is_corrupted",
        "injected_displacement_um",
        "realized_displacement_um",
        "realized_error_from_true_um",
    }
    assert expected <= set(truth.colnames)


def test_selected_trace_without_eligible_detection_is_reported_not_resampled():
    observed, truth, _, summary = simulate_traces(
        _config(n_traces=3, detection_efficiency=0), include_auxiliary=True
    )
    assert len(observed) == 0
    assert len(truth) == 0
    assert np.all(summary["selected_for_corruption"])
    assert not np.any(summary["corruption_target_available"])
    assert set(summary["status"]) == {"no_eligible_detected_barcode"}


def test_edge_margin_can_make_all_detected_barcodes_ineligible():
    _, truth = simulate_traces(_config(n_barcodes=4, mislocalization_edge_margin=2))
    assert np.all(truth["Mislocalization_Selected_Trace"])
    assert not np.any(truth["Mislocalization_Applied"])


def test_unit_vector_generator_is_isotropic():
    directions = random_unit_vectors(np.random.default_rng(8), 30000)
    np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1, atol=1e-14)
    # Uniform sphere components have mean zero and second moment 1/3.
    np.testing.assert_allclose(directions.mean(axis=0), 0, atol=0.015)
    np.testing.assert_allclose((directions**2).mean(axis=0), 1 / 3, atol=0.015)
