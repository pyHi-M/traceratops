import gc
import importlib.util
import sys
import weakref
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "analyze_trace_curator_loo_validation.py"
)
SPEC = importlib.util.spec_from_file_location("loo_validation", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


class Reference:
    minimum_observations = 1
    separation = {float(i): np.linspace(i - 0.1, i + 0.1, 101) for i in range(1, 8)}


def trace(outliers=(2,)):
    x = np.arange(7, dtype=float)
    for index in outliers:
        x[index] += 20
    return pd.DataFrame(
        {
            "simulation_id": "s",
            "Trace_ID": "t",
            "Spot_ID": [f"p{i}" for i in range(7)],
            "Barcode": np.arange(7),
            "x": x,
            "y": 0.0,
            "z": 0.0,
            "is_corrupted": [i in outliers for i in range(7)],
        }
    )


def test_exact_loo_top3_calculation_and_multiple_ground_truth():
    scored, summary = analysis.score_trace(trace((2, 5)), Reference())
    xyz = scored[["x", "y", "z"]].to_numpy(float)
    barcode = scored["Barcode"].to_numpy(float)
    edges = []
    for i in range(7):
        for j in range(i + 1, 7):
            edges.append(
                (
                    i,
                    j,
                    analysis.empirical_anomaly(
                        Reference.separation[float(barcode[j] - barcode[i])],
                        np.linalg.norm(xyz[j] - xyz[i]),
                    ),
                )
            )
    before = analysis._top3([v for _, _, v in edges])
    expected = [
        before - analysis._top3([v for i, j, v in edges if node not in (i, j)])
        for node in range(7)
    ]
    np.testing.assert_allclose(scored["loo_top3_all"], expected)
    assert summary["C_top3"] == before and summary["n_true_corruptions"] == 2


def test_iterative_rescoring_can_find_second_culprit():
    rows = analysis.iterative_scores(trace((2, 5)), Reference())
    assert rows[0]["predicted_is_true_corruption"]
    assert rows[1]["predicted_is_true_corruption"]
    assert rows[1]["remaining_true_corruptions_before_step"] == 1


def test_seed_leakage_and_step_specific_thresholds():
    ref, cal, evaluation = analysis.fold_seed_sets([1, 2, 3, 4, 5], 5, 4)
    assert ref == {1, 2, 3} and cal == {4} and evaluation == {5}
    assert analysis.trace_threshold(range(1000), 0.01) == 989
    clean = pd.concat([trace(()).assign(Trace_ID=f"t{i}") for i in range(3)])
    rows = analysis.step_null_scores(clean, Reference())
    assert set(row["step"] for row in rows) == {0, 1, 2}


def test_margin_ratio_edges_and_ties_are_controlled():
    zero = analysis.confidence_metrics([1.0, 0.0])
    assert (
        np.isfinite(zero["score_ratio"]) and zero["score_ratio"] == 1 / analysis.EPSILON
    )
    tied = analysis.confidence_metrics([2.0, 2.0, 1.0])
    assert tied["absolute_margin"] == 0 and tied["relative_margin"] == 0
    assert tied["score_ratio"] == 1 and not tied["winner_is_unique"]


def test_k_and_barcode_metadata_and_condition_memory_release():
    for k in (1, 2, 3):
        n, value = analysis._condition_metadata(
            {"id": "c", "n_barcodes": 25, "corrupted_barcodes_per_trace": k}
        )
        assert (n, value) == (25, k)
    frame = trace()
    reference = weakref.ref(frame)
    del frame
    gc.collect()
    assert reference() is None


def test_plot_only_uses_written_compact_tables(tmp_path):
    frames = {
        analysis.OUTPUTS[0]: pd.DataFrame(
            {
                "displacement": [0.2],
                "target_sensitivity": [1.0],
                "unique_top1_attribution": [1.0],
                "n_barcodes": [25],
                "detection_efficiency": [0.5],
            }
        ),
        analysis.OUTPUTS[1]: pd.DataFrame({"K": [2]}),
        analysis.OUTPUTS[2]: pd.DataFrame(
            {
                "K": [2],
                "n_barcodes": [25],
                "cumulative_recall_step_1": [0.5],
                "cumulative_recall_step_2": [1.0],
                "cumulative_recall_step_3": [1.0],
                "precision_step_1": [1.0],
                "precision_step_2": [1.0],
                "precision_step_3": [0.0],
                "exact_recovery_after_k": [1.0],
            }
        ),
        analysis.OUTPUTS[3]: pd.DataFrame(
            {"K": [2], "n_barcodes": [25], "false_removal_rate": [0.0]}
        ),
        analysis.OUTPUTS[4]: pd.DataFrame(
            {
                "S1": [2.0],
                "S2": [1.0],
                "absolute_margin": [1.0],
                "top_candidate_is_true": [True],
            }
        ),
        analysis.OUTPUTS[5]: pd.DataFrame(
            {"metric": ["absolute_margin"], "coverage": [1.0], "precision": [1.0]}
        ),
        analysis.OUTPUTS[6]: pd.DataFrame({"threshold": [1.0]}),
        analysis.OUTPUTS[7]: pd.DataFrame({"n_observations": [10]}),
    }
    for name, frame in frames.items():
        analysis._write(frame, tmp_path / name)
    analysis.plot_only(tmp_path)
    assert (tmp_path / "precision_coverage.png").is_file()


def test_tiny_end_to_end_manifest_writes_required_outputs(tmp_path):
    root = tmp_path / "benchmark"
    simulations = []
    conditions = []
    for seed in (1, 2, 3):
        simulation_id = f"base-{seed}"
        simulations.append(
            {
                "id": simulation_id,
                "seed": seed,
                "n_barcodes": 5,
                "detection_efficiency": 0.5,
            }
        )
        baseline = pd.DataFrame(
            [
                {
                    "Spot_ID": f"b-{seed}-{trace_id}-{barcode}",
                    "Trace_ID": f"t{trace_id}",
                    "Barcode #": barcode,
                    "x": float(barcode) + trace_id / 100,
                    "y": 0.0,
                    "z": 0.0,
                }
                for trace_id in range(3)
                for barcode in range(5)
            ]
        )
        baseline_dir = root / "simulations" / simulation_id
        baseline_dir.mkdir(parents=True)
        analysis._write(baseline, baseline_dir / "simulated.ecsv")

        condition_id = f"condition-{seed}"
        condition_dir = root / "conditions" / condition_id
        condition_dir.mkdir(parents=True)
        observed = baseline.copy()
        observed["Spot_ID"] = [f"e-{seed}-{i}" for i in range(len(observed))]
        target = (observed["Trace_ID"] == "t0") & (observed["Barcode #"] == 2)
        observed.loc[target, "x"] += 20
        truth = pd.DataFrame(
            {
                "Spot_ID": observed["Spot_ID"],
                "Input_Trace_ID": observed["Trace_ID"],
                "Barcode #": observed["Barcode #"],
                "selected_for_corruption": target,
                "is_corrupted": target,
                "corruption_target_rank": np.where(target, 1, 0),
                "injected_displacement_um": np.where(target, 0.8, 0.0),
                "realized_displacement_um": np.where(target, 0.8, 0.0),
            }
        )
        analysis._write(observed, condition_dir / "simulated.ecsv")
        analysis._write(truth, condition_dir / "simulated.ground_truth.ecsv")
        conditions.append(
            {
                "id": condition_id,
                "directory": f"conditions/{condition_id}",
                "simulation_id": simulation_id,
                "replicate_index": seed,
                "seed": seed,
                "n_barcodes": 5,
                "corrupted_barcodes_per_trace": 1,
                "detection_efficiency": 0.5,
                "displacement_um": 0.8,
            }
        )
    root.mkdir(exist_ok=True)
    (root / "sweep_manifest.yaml").write_text(
        yaml.safe_dump(
            {"dry_run": False, "simulations": simulations, "conditions": conditions}
        )
    )
    output = tmp_path / "output"
    analysis.run(root, output, minimum_observations=1)
    assert all((output / name).is_file() for name in analysis.OUTPUTS)
    assert (output / "loo_iterative_diagnostics.ecsv").is_file()
