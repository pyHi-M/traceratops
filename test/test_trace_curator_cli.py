"""Exercise the installed command and its transparent persistence/ECSV path."""

import hashlib
import json
import shutil
import subprocess
import sysconfig
from dataclasses import fields, replace

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table, vstack

from traceratops.core.trace_curator import (
    TraceCurationIteration,
    curate_trace,
    fit_curator_model,
)
from traceratops.core.trace_curator_io import load_curator_model, save_curator_model


def trace(identity="trace", n=4, barcode_column="Barcode #", numeric_ids=False):
    table = Table(
        {
            "Trace_ID": [identity] * n,
            "Spot_ID": (
                np.arange(2**53 + 1, 2**53 + 1 + n, dtype=np.int64)
                if numeric_ids
                else [f"spot-{i}" for i in range(n)]
            ),
            barcode_column: [f"locus-{i}" for i in range(n)],
            "Genomic_Position": [0.0, 10.0, 30.0, 70.0][:n],
            "x": np.arange(n, dtype=float),
            "y": np.zeros(n),
            "z": np.zeros(n),
            "extra": np.arange(n, dtype=np.int16),
            "Chrom": ["chr1"] * n,
        }
    )
    table.meta["experiment"] = "preserve me"
    return table


def population(n=100, **kwargs):
    return vstack([trace(identity=f"trace-{i}", **kwargs) for i in range(n)])


def write(table, path):
    table.write(path, format="ascii.ecsv", overwrite=True)
    return str(path)


def read(path):
    return Table.read(path, format="ascii.ecsv")


def cli(*args, success=True):
    executable = shutil.which("trace_curator", path=sysconfig.get_path("scripts"))
    assert executable is not None, "Install traceratops to test its console entry point"
    result = subprocess.run(
        [executable, *map(str, args)],
        capture_output=True,
        text=True,
    )
    if success:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode != 0
        assert "Traceback" not in result.stderr
    return result


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    directory = tmp_path_factory.mktemp("curator-model")
    with pytest.warns(UserWarning, match="may be unstable"):
        model = fit_curator_model(population())
    prefix = directory / "model"
    save_curator_model(model, prefix)
    return model, prefix


def assert_same_table(left, right):
    assert left.colnames == right.colnames
    for name in left.colnames:
        assert left[name].dtype == right[name].dtype
        mask = np.ma.getmaskarray(left[name])
        np.testing.assert_array_equal(mask, np.ma.getmaskarray(right[name]))
        np.testing.assert_array_equal(left[name][~mask], right[name][~mask])


@pytest.mark.parametrize(
    "arguments",
    [("--help",), ("fit", "--help"), ("curate", "--help"), ("run", "--help")],
)
def test_help(arguments):
    result = cli(*arguments)
    assert "usage:" in result.stdout
    assert result.stderr == ""
    assert "--genomic-positions" not in result.stdout


@pytest.mark.parametrize("command", [[], ["fit"], ["curate"], ["run"]])
def test_required_invocation_arguments(command):
    assert "required" in cli(*command, success=False).stderr


def test_fit_default_crossfit_and_summary(tmp_path):
    source = write(population(), tmp_path / "input.ecsv")
    prefix = tmp_path / "model"
    result = cli("fit", "--input", source, "--output-model", prefix)
    assert prefix.with_suffix(".npz").exists()
    assert prefix.with_suffix(".json").exists()
    model = load_curator_model(prefix)
    assert model.calibration_source == "crossfit"
    assert model.n_crossfit_folds == 5
    assert model.n_calibration_traces == 100
    assert "may be unstable" in result.stderr
    for field in fields(model):
        if field.name not in {"reference_distributions", "genomic_positions"}:
            assert f"{field.name}: {getattr(model, field.name)}" in result.stdout


def test_fit_independent_controls(tmp_path):
    source = write(population(20), tmp_path / "input.ecsv")
    control = write(population(), tmp_path / "controls.ecsv")
    prefix = tmp_path / "model"
    cli(
        "fit",
        "--input",
        source,
        "--calibration-input",
        control,
        "--output-model",
        prefix,
    )
    model = load_curator_model(prefix)
    assert model.calibration_source == "separate_traces"
    assert model.n_crossfit_folds is None
    assert model.n_reference_traces == 20


@pytest.mark.parametrize("mode", ["separation", "barcode_pair"])
@pytest.mark.parametrize(
    "barcode_values",
    [
        ["alpha", "beta", "gamma", "delta"],
        [1, 2, 3, 4],
        ["1", "01", "a/b", "α"],
        [1.25, 2.5, 3.75, 4.5],
    ],
)
def test_model_roundtrip(mode, barcode_values, tmp_path):
    base = trace(numeric_ids=True)
    base["Barcode #"] = barcode_values
    with pytest.warns(UserWarning):
        model = fit_curator_model([base] * 100, reference_mode=mode)
    prefix = tmp_path / "roundtrip"
    save_curator_model(model, prefix)
    loaded = load_curator_model(prefix.with_suffix(".json"))
    assert loaded.global_threshold == model.global_threshold
    assert dict(loaded.genomic_positions) == dict(model.genomic_positions)
    for field in fields(model):
        if field.name not in {"reference_distributions", "genomic_positions"}:
            assert getattr(loaded, field.name) == getattr(model, field.name)
    assert set(loaded.reference_distributions) == set(model.reference_distributions)
    for key, sample in model.reference_distributions.items():
        np.testing.assert_array_equal(loaded.reference_distributions[key], sample)
        with pytest.raises(ValueError):
            loaded.reference_distributions[key].setflags(write=True)
    with np.load(prefix.with_suffix(".npz"), allow_pickle=False) as archive:
        assert all(archive[name].dtype == np.dtype("float64") for name in archive.files)
    base["y"][2] = 100.0
    original, restored = curate_trace(base, model), curate_trace(base, loaded)
    assert original.removed_localizations.equals(restored.removed_localizations)
    assert original.curated_trace.equals(restored.curated_trace)
    assert original.terminal_state == restored.terminal_state == "curated"
    assert original.initial_cost == restored.initial_cost
    assert original.final_cost == restored.final_cost


def test_mixed_typed_barcode_identity_persistence(tmp_path, fitted):
    model = replace(
        fitted[0],
        reference_mode="barcode_pair",
        genomic_positions={"1": 0.0, 1: 10.0, False: 30.0},
        reference_distributions={
            frozenset(("1", 1)): np.array([1.0, 2.0]),
            frozenset(("1", False)): np.array([3.0]),
        },
    )
    save_curator_model(model, tmp_path / "typed")
    loaded = load_curator_model(tmp_path / "typed")
    assert loaded.genomic_positions == model.genomic_positions
    assert set(loaded.reference_distributions) == set(model.reference_distributions)


def test_multi_trace_curation_outputs_and_core_equivalence(tmp_path, fitted):
    model, prefix = fitted
    base = trace(identity="normal", numeric_ids=True)
    bad = trace(identity="curated", numeric_ids=True)
    bad["y"][2] = 100.0
    unresolved = trace(identity="unresolved", n=3, numeric_ids=True)
    unresolved["y"][0] = 100.0
    small = trace(identity="small", n=2, numeric_ids=True)
    source = vstack([base, bad, unresolved, small])
    source["optional"] = MaskedColumn(
        np.arange(len(source), dtype=np.int32),
        mask=[True] + [False] * (len(source) - 1),
    )
    input_path = write(source, tmp_path / "input.ecsv")
    paths = [
        tmp_path / f"{name}.ecsv"
        for name in ("out", "removed", "summary", "diagnostic")
    ]
    cli(
        "curate",
        "--input",
        input_path,
        "--model",
        prefix,
        "--output",
        paths[0],
        "--removed-output",
        paths[1],
        "--summary-output",
        paths[2],
        "--diagnostics-output",
        paths[3],
    )
    curated, discarded, summaries, diagnostics = map(read, paths)
    assert_same_table(curated, source[[i for i in range(len(source)) if i != 6]])
    assert curated.meta == source.meta
    assert discarded["Spot_ID"].tolist() == [bad["Spot_ID"][2]]
    assert discarded["Trace_ID"].tolist() == ["curated"]
    assert discarded["curator_removal_iteration"].tolist() == [0]
    assert discarded["curator_candidate_barcode"].tolist() == ["locus-2"]
    assert summaries["terminal_state"].tolist() == [
        "normal",
        "curated",
        "abnormal_unresolved",
        "unscoreable",
    ]
    assert summaries["n_removed"].tolist() == [0, 1, 0, 0]
    assert summaries["n_input_localizations"].tolist() == [4, 4, 3, 2]
    assert summaries["n_output_localizations"].tolist() == [4, 3, 3, 2]
    assert summaries["n_iterations"].tolist() == [1, 2, 1, 1]
    assert diagnostics.colnames == [
        "Trace_ID",
        *(field.name for field in fields(TraceCurationIteration)),
    ]
    assert len(diagnostics) == 5
    assert diagnostics["removal_accepted"].tolist() == [
        False,
        True,
        False,
        False,
        False,
    ]
    assert np.ma.getmaskarray(diagnostics["candidate_spot_id"]).tolist() == [
        True,
        False,
        True,
        True,
        True,
    ]
    assert diagnostics["candidate_spot_id"].dtype.kind in "iu"
    assert diagnostics["candidate_spot_id"][1] == bad["Spot_ID"][2]
    for offset, (identity, group) in enumerate(
        source.to_pandas().groupby("Trace_ID", sort=False)
    ):
        result = curate_trace(group, model)
        for name in (
            "terminal_state",
            "stop_reason",
            "initial_cost",
            "final_cost",
            "global_threshold",
            "n_removed",
        ):
            actual, expected = summaries[name][offset], getattr(result, name)
            assert actual == expected or (
                isinstance(expected, float) and np.isnan(actual) and np.isnan(expected)
            )
        for row, iteration in zip(
            diagnostics[diagnostics["Trace_ID"] == identity], result.iterations
        ):
            for field in fields(TraceCurationIteration):
                value = getattr(iteration, field.name)
                if value is None:
                    assert np.ma.is_masked(row[field.name])
                elif isinstance(value, float) and np.isnan(value):
                    assert np.isnan(row[field.name])
                else:
                    assert row[field.name] == value


@pytest.mark.parametrize("barcode_column", ["Barcode", "Barcode #"])
def test_row_order_independent_curation_preserves_retained_order(
    tmp_path, fitted, barcode_column
):
    frame = trace(barcode_column=barcode_column)
    frame["y"][2] = 100.0
    frame = frame[[3, 2, 0, 1]]
    source = write(frame, tmp_path / "input.ecsv")
    cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        "--removed-output",
        tmp_path / "removed.ecsv",
    )
    assert read(tmp_path / "out.ecsv")["Spot_ID"].tolist() == [
        "spot-3",
        "spot-0",
        "spot-1",
    ]
    assert read(tmp_path / "removed.ecsv")["Spot_ID"].tolist() == ["spot-2"]


def test_interleaved_trace_rows_use_first_appearance_order(tmp_path, fitted):
    source = vstack([trace(identity="z"), trace(identity="a")])[
        [0, 4, 1, 5, 2, 6, 3, 7]
    ]
    path = write(source, tmp_path / "input.ecsv")
    cli(
        "curate",
        "--input",
        path,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
    )
    assert_same_table(
        read(tmp_path / "out.ecsv"), vstack([trace(identity="z"), trace(identity="a")])
    )


def test_empty_outputs_have_typed_schemas(tmp_path, fitted):
    source = write(trace()[:0], tmp_path / "empty.ecsv")
    paths = [
        tmp_path / f"{name}.ecsv"
        for name in ("out", "removed", "summary", "diagnostic")
    ]
    cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        paths[0],
        "--removed-output",
        paths[1],
        "--summary-output",
        paths[2],
        "--diagnostics-output",
        paths[3],
    )
    assert all(len(read(path)) == 0 for path in paths)
    assert read(paths[0]).colnames == trace().colnames
    assert "curator_removal_iteration" in read(paths[1]).colnames


@pytest.mark.parametrize("independent", [False, True])
def test_run_matches_explicit_fit_then_curate(tmp_path, independent):
    data = population(20 if independent else 100)
    data["y"][2] = 100.0
    source = write(data, tmp_path / "input.ecsv")
    calibration = (
        ["--calibration-input", write(population(), tmp_path / "controls.ecsv")]
        if independent
        else []
    )
    cli("fit", "--input", source, "--output-model", tmp_path / "explicit", *calibration)

    def outputs(prefix):
        return [
            "--output",
            tmp_path / f"{prefix}-out.ecsv",
            "--removed-output",
            tmp_path / f"{prefix}-removed.ecsv",
            "--summary-output",
            tmp_path / f"{prefix}-summary.ecsv",
            "--diagnostics-output",
            tmp_path / f"{prefix}-diagnostic.ecsv",
        ]

    cli(
        "curate",
        "--input",
        source,
        "--model",
        tmp_path / "explicit",
        *outputs("explicit"),
    )
    cli(
        "run",
        "--input",
        source,
        "--output-model",
        tmp_path / "run",
        *outputs("run"),
        *calibration,
    )
    for suffix in ("out", "removed", "summary", "diagnostic"):
        assert_same_table(
            read(tmp_path / f"run-{suffix}.ecsv"),
            read(tmp_path / f"explicit-{suffix}.ecsv"),
        )
    assert (
        load_curator_model(tmp_path / "run").global_threshold
        == load_curator_model(tmp_path / "explicit").global_threshold
    )


@pytest.mark.parametrize(
    "column", ["Trace_ID", "Genomic_Position", "x", "y", "z", "Spot_ID", "Barcode #"]
)
@pytest.mark.parametrize("command", ["fit", "curate"])
def test_missing_required_columns(tmp_path, fitted, column, command):
    data = trace()
    data.remove_column(column)
    source = write(data, tmp_path / "input.ecsv")
    arguments = (
        ["--output-model", tmp_path / "model"]
        if command == "fit"
        else ["--model", fitted[1], "--output", tmp_path / "out.ecsv"]
    )
    result = cli(command, "--input", source, *arguments, success=False)
    assert column in result.stderr
    if column == "Genomic_Position":
        assert "preprocessing tool" in result.stderr
    assert not (tmp_path / "out.ecsv").exists()


@pytest.mark.parametrize(
    "mutation,message",
    [
        ("barcode", "resolve duplicated"),
        ("spot", "Spot_ID must be unique"),
        ("alias", "columns disagree"),
        ("genomic_duplicate", "finite and unique"),
        ("genomic_model", "disagrees"),
        ("nonfinite", "coordinates must be finite"),
        ("chrom_mixed", "one chromosome"),
        ("chrom_model", "chromosome disagrees"),
    ],
)
def test_invalid_curation_inputs(tmp_path, fitted, mutation, message):
    data = trace()
    if mutation == "barcode":
        data["Barcode #"][1] = data["Barcode #"][0]
    if mutation == "spot":
        data["Spot_ID"][1] = data["Spot_ID"][0]
    if mutation == "alias":
        data["Barcode"] = ["other"] * len(data)
    if mutation == "genomic_duplicate":
        data["Genomic_Position"][1] = 0.0
    if mutation == "genomic_model":
        data["Genomic_Position"] += 1
    if mutation == "nonfinite":
        data["x"][0] = np.nan
    if mutation == "chrom_mixed":
        data["Chrom"][1] = "chr2"
    if mutation == "chrom_model":
        data["Chrom"] = ["chr2"] * len(data)
    source = write(data, tmp_path / "input.ecsv")
    result = cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        success=False,
    )
    assert message in result.stderr


def test_invalid_population_mapping_and_sample_count(tmp_path):
    data = population()
    data["Genomic_Position"][4] = 1.0
    source = write(data, tmp_path / "bad.ecsv")
    assert (
        "disagrees"
        in cli(
            "fit",
            "--input",
            source,
            "--output-model",
            tmp_path / "model",
            success=False,
        ).stderr
    )
    source = write(population(99), tmp_path / "few.ecsv")
    assert (
        "Too few scoreable"
        in cli(
            "fit",
            "--input",
            source,
            "--output-model",
            tmp_path / "model",
            success=False,
        ).stderr
    )
    assert not (tmp_path / "model.json").exists()


def test_unsupported_mode_and_negative_removal_limit(tmp_path, fitted):
    source = write(trace(), tmp_path / "input.ecsv")
    assert (
        "invalid choice"
        in cli(
            "fit",
            "--input",
            source,
            "--output-model",
            tmp_path / "model",
            "--reference-mode",
            "unknown",
            success=False,
        ).stderr
    )
    assert (
        "nonnegative"
        in cli(
            "curate",
            "--input",
            source,
            "--model",
            fitted[1],
            "--output",
            tmp_path / "out.ecsv",
            "--max-removals",
            "-1",
            success=False,
        ).stderr
    )


@pytest.mark.parametrize(
    "column",
    [
        "curator_removal_iteration",
        "curator_candidate_barcode",
        "curator_initial_cost",
        "curator_global_threshold",
    ],
)
def test_provenance_collision_rejected(tmp_path, fitted, column):
    data = trace()
    data[column] = np.zeros(len(data))
    source = write(data, tmp_path / "input.ecsv")
    result = cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        "--removed-output",
        tmp_path / "removed.ecsv",
        success=False,
    )
    assert column in result.stderr
    assert not (tmp_path / "out.ecsv").exists()


def test_output_collisions_and_existing_output_policy(tmp_path, fitted):
    source = write(trace(), tmp_path / "input.ecsv")
    assert (
        "conflicts"
        in cli(
            "curate",
            "--input",
            source,
            "--model",
            fitted[1],
            "--output",
            source,
            success=False,
        ).stderr
    )
    assert (
        "conflicts"
        in cli(
            "curate",
            "--input",
            source,
            "--model",
            fitted[1],
            "--output",
            tmp_path / "same.ecsv",
            "--summary-output",
            tmp_path / "same.ecsv",
            success=False,
        ).stderr
    )
    (tmp_path / "out.ecsv").write_text("replace this")
    cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
    )
    assert_same_table(read(tmp_path / "out.ecsv"), trace())


@pytest.mark.parametrize(
    "damage", ["json", "version", "checksum", "array", "object_array", "missing_file"]
)
def test_malformed_model_fails_clearly(tmp_path, fitted, damage):
    save_curator_model(fitted[0], tmp_path / "broken")
    json_path, npz_path = tmp_path / "broken.json", tmp_path / "broken.npz"
    metadata = json.loads(json_path.read_text())
    if damage == "json":
        json_path.write_text("invalid json")
    if damage == "version":
        metadata["format_version"] = 99
        json_path.write_text(json.dumps(metadata))
    if damage == "checksum":
        npz_path.write_bytes(b"broken")
    if damage in {"array", "object_array"}:
        with np.load(npz_path, allow_pickle=False) as archive:
            arrays = {name: archive[name] for name in archive.files}
        name = next(iter(arrays))
        arrays[name] = (
            np.array(["unsafe"], dtype=object)
            if damage == "object_array"
            else np.array([2.0, 1.0])
        )
        np.savez(npz_path, **arrays)
        metadata["npz_sha256"] = hashlib.sha256(npz_path.read_bytes()).hexdigest()
        json_path.write_text(json.dumps(metadata))
    if damage == "missing_file":
        npz_path.unlink()
    source = write(trace(), tmp_path / "input.ecsv")
    result = cli(
        "curate",
        "--input",
        source,
        "--model",
        tmp_path / "broken",
        "--output",
        tmp_path / "out.ecsv",
        success=False,
    )
    assert "Cannot load curator model" in result.stderr


def test_unsupported_identity_type_fails_clearly(tmp_path, fitted):
    model = replace(fitted[0], genomic_positions={tuple(["tuple"]): 0.0})
    with pytest.raises(ValueError, match="Model identities must be"):
        save_curator_model(model, tmp_path / "unsupported")


def test_string_candidate_identity_and_empty_stop_reason(tmp_path, fitted):
    source = trace()
    source["y"][2] = 100.0
    path = write(source, tmp_path / "input.ecsv")
    cli(
        "curate",
        "--input",
        path,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        "--diagnostics-output",
        tmp_path / "diagnostic.ecsv",
    )
    diagnostic = read(tmp_path / "diagnostic.ecsv")
    assert diagnostic["candidate_spot_id"][0] == "spot-2"
    assert not np.ma.is_masked(diagnostic["candidate_spot_id"][0])
    assert np.ma.is_masked(diagnostic["candidate_spot_id"][1])
    assert diagnostic["stop_reason"][0] == ""
    assert not np.ma.is_masked(diagnostic["stop_reason"][0])


def test_zero_removal_requested_outputs(tmp_path, fitted):
    source = write(trace(), tmp_path / "input.ecsv")
    cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        "--removed-output",
        tmp_path / "removed.ecsv",
    )
    removed = read(tmp_path / "removed.ecsv")
    assert len(removed) == 0
    assert "curator_candidate_barcode" in removed.colnames
    assert removed["Spot_ID"].dtype.kind in "US"


def test_calibration_input_requires_genomic_position(tmp_path):
    source = write(population(20), tmp_path / "input.ecsv")
    control = trace()
    control.remove_column("Genomic_Position")
    path = write(control, tmp_path / "control.ecsv")
    result = cli(
        "run",
        "--input",
        source,
        "--calibration-input",
        path,
        "--output-model",
        tmp_path / "model",
        "--output",
        tmp_path / "out.ecsv",
        success=False,
    )
    assert "Genomic_Position" in result.stderr
    assert not (tmp_path / "model.json").exists()


def test_multiple_removals_follow_api_order(tmp_path):
    positions = np.arange(10, dtype=float)
    base = Table(
        {
            "Trace_ID": ["trace"] * 10,
            "Spot_ID": [f"spot-{i}" for i in range(10)],
            "Barcode": [f"locus-{i}" for i in range(10)],
            "Genomic_Position": positions,
            "x": positions.copy(),
            "y": np.zeros(10),
            "z": np.zeros(10),
        }
    )
    with pytest.warns(UserWarning):
        model = fit_curator_model([base] * 100)
    save_curator_model(model, tmp_path / "model")
    base["y"][0], base["y"][5] = 100.0, 600.0
    source = write(base, tmp_path / "input.ecsv")
    cli(
        "curate",
        "--input",
        source,
        "--model",
        tmp_path / "model",
        "--output",
        tmp_path / "out.ecsv",
        "--removed-output",
        tmp_path / "removed.ecsv",
    )
    result = curate_trace(base, model)
    removed = read(tmp_path / "removed.ecsv")
    assert (
        removed["Spot_ID"].tolist()
        == result.removed_localizations["Spot_ID"].tolist()
        == ["spot-5", "spot-0"]
    )
    assert removed["curator_removal_iteration"].tolist() == [0, 1]


def test_parser_defaults_match_core_api():
    import inspect

    from traceratops.trace_curator import parse_arguments

    args = parse_arguments().parse_args(
        ["run", "--input", "in.ecsv", "--output", "out.ecsv", "--output-model", "model"]
    )
    for name in (
        "trace_fpr",
        "minimum_reference_observations",
        "reference_mode",
        "n_crossfit_folds",
    ):
        assert (
            getattr(args, name)
            == inspect.signature(fit_curator_model).parameters[name].default
        )
    assert (
        args.max_removals
        == inspect.signature(curate_trace).parameters["max_removals"].default
    )


def test_shared_ecsv_helpers_retain_old_imports(tmp_path):
    from traceratops.core import chromatin_trace_table, io_manager

    assert chromatin_trace_table.read_table_from_ecsv is io_manager.read_table_from_ecsv
    assert chromatin_trace_table.save_table_to_ecsv is io_manager.save_table_to_ecsv
    path = tmp_path / "shared.ecsv"
    io_manager.save_table_to_ecsv(trace(), path)
    assert_same_table(io_manager.read_table_from_ecsv(path), trace())


def test_missing_trace_identity_is_not_silently_dropped(tmp_path, fitted):
    data = trace()
    data["Trace_ID"] = [1.0, 1.0, np.nan, 1.0]
    source = write(data, tmp_path / "input.ecsv")
    result = cli(
        "curate",
        "--input",
        source,
        "--model",
        fitted[1],
        "--output",
        tmp_path / "out.ecsv",
        success=False,
    )
    assert "Trace_ID must not be missing" in result.stderr
    assert not (tmp_path / "out.ecsv").exists()
