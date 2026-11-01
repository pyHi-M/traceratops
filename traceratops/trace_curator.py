"""Fit empirical references and curate resolved chromatin traces in ECSV files."""

import argparse
import os
import tempfile
from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.table import Column, MaskedColumn, Table

from traceratops.core.io_manager import (
    read_table_from_ecsv,
    save_table_to_ecsv,
)
from traceratops.core.trace_curator import (
    TraceCurationIteration,
    curate_trace,
    fit_curator_model,
)
from traceratops.core.trace_curator_io import (
    load_curator_model,
    model_paths,
    save_curator_model,
)

_PROVENANCE = (
    "curator_removal_iteration",
    "curator_candidate_barcode",
    "curator_initial_cost",
    "curator_global_threshold",
)
_SUMMARY = (
    "terminal_state",
    "stop_reason",
    "initial_cost",
    "final_cost",
    "global_threshold",
    "n_removed",
    "n_input_localizations",
    "n_output_localizations",
    "n_iterations",
)


def parse_arguments():
    """Build the production fit/curate/run parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("fit", "curate", "run"):
        command = commands.add_parser(
            name,
            help={
                "fit": "Fit and save a reference model",
                "curate": "Curate with a saved model",
                "run": "Fit and curate the same input population",
            }[name],
        )
        command.add_argument(
            "--input", required=True, help="Input ECSV with genomic positions"
        )
        if name in {"fit", "run"}:
            command.add_argument(
                "--output-model",
                required=True,
                help="Output model prefix (.npz and .json)",
            )
            command.add_argument("--trace-fpr", type=float, default=0.01)
            command.add_argument(
                "--minimum-reference-observations", type=int, default=20
            )
            command.add_argument(
                "--reference-mode",
                choices=("separation", "barcode_pair"),
                default="separation",
            )
            command.add_argument("--n-crossfit-folds", type=int, default=5)
            command.add_argument(
                "--calibration-input", help="Independent calibration-control ECSV"
            )
        if name in {"curate", "run"}:
            command.add_argument("--output", required=True, help="Curated ECSV output")
            command.add_argument(
                "--removed-output", help="Removed rows with provenance"
            )
            command.add_argument("--summary-output", help="One summary row per trace")
            command.add_argument(
                "--diagnostics-output",
                help="One row per iteration, including final state",
            )
            command.add_argument("--max-removals", type=int, default=3)
        if name == "curate":
            command.add_argument(
                "--model", required=True, help="Model prefix, .json, or .npz path"
            )
    return parser


def _read_input(path):
    table = read_table_from_ecsv(path)
    if "Genomic_Position" not in table.colnames:
        raise ValueError(
            "Input trace table must contain 'Genomic_Position'. Add genomic coordinates using the existing traceratops preprocessing tool before running trace_curator."
        )
    missing = {"Trace_ID", "Spot_ID", "x", "y", "z"} - set(table.colnames)
    if not ({"Barcode", "Barcode #"} & set(table.colnames)):
        missing.add("Barcode or Barcode #")
    if missing:
        raise ValueError(
            f"Input trace table is missing required columns: {', '.join(sorted(missing))}"
        )
    if (
        np.ma.getmaskarray(table["Trace_ID"]).any()
        or pd.isna(np.asarray(table["Trace_ID"])).any()
    ):
        raise ValueError("Trace_ID must not be missing")
    return table


def _check_paths(args):
    inputs = [args.input]
    if getattr(args, "calibration_input", None):
        inputs.append(args.calibration_input)
    if getattr(args, "model", None):
        inputs.extend(model_paths(args.model))
    outputs = [
        getattr(args, name, None)
        for name in ("output", "removed_output", "summary_output", "diagnostics_output")
    ]
    if getattr(args, "output_model", None):
        outputs.extend(model_paths(args.output_model))
    resolved_inputs = {Path(path).resolve() for path in inputs}
    resolved_outputs = set()
    for filename in filter(None, outputs):
        path = Path(filename).resolve()
        if path in resolved_inputs or path in resolved_outputs:
            raise ValueError(
                f"Output path conflicts with an input or another output: {filename}"
            )
        if path.is_dir() or not path.parent.is_dir():
            raise ValueError(
                f"Output must be a file in an existing directory: {filename}"
            )
        resolved_outputs.add(path)


def _identity_column(name, values, source):
    mask = [value is None for value in values]
    fill = "" if source.dtype.kind in "US" else 0
    return MaskedColumn(
        [fill if value is None else value for value in values],
        dtype=source.dtype,
        mask=mask,
        name=name,
    )


def _curate_table(table, model, args):
    if args.removed_output:
        conflicts = set(_PROVENANCE) & set(table.colnames)
        if conflicts:
            raise ValueError(
                f"Removal provenance columns already exist: {', '.join(sorted(conflicts))}"
            )
    barcode_name = "Barcode #" if "Barcode #" in table.colnames else "Barcode"
    frame = table.to_pandas()
    curated_indices: list[int] = []
    removed_indices: list[int] = []
    summaries: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    provenance: list[tuple] = []
    for trace_id, trace in frame.groupby("Trace_ID", sort=False, observed=True):
        result = curate_trace(trace, model, max_removals=args.max_removals)
        curated_indices.extend(result.curated_trace.index.tolist())
        if args.removed_output:
            removed_indices.extend(result.removed_localizations.index.tolist())
            provenance.extend(
                (
                    iteration.iteration,
                    iteration.candidate_barcode,
                    result.initial_cost,
                    result.global_threshold,
                )
                for iteration in result.iterations
                if iteration.removal_accepted
            )
        if args.summary_output:
            summaries.append(
                {
                    "Trace_ID": trace_id,
                    **{name: getattr(result, name) for name in _SUMMARY[:6]},
                    "n_input_localizations": len(trace),
                    "n_output_localizations": len(result.curated_trace),
                    "n_iterations": len(result.iterations),
                }
            )
        if args.diagnostics_output:
            diagnostics.extend(
                {
                    "Trace_ID": trace_id,
                    **{
                        field.name: getattr(iteration, field.name)
                        for field in fields(TraceCurationIteration)
                    },
                }
                for iteration in result.iterations
            )
    outputs = {args.output: table[np.asarray(curated_indices, dtype=int)]}
    if args.removed_output:
        removed = table[np.asarray(removed_indices, dtype=int)]
        for offset, (name, dtype) in enumerate(
            zip(_PROVENANCE, (int, table[barcode_name].dtype, float, float))
        ):
            removed.add_column(
                Column([row[offset] for row in provenance], name=name, dtype=dtype)
            )
        outputs[args.removed_output] = removed
    if args.summary_output:
        summary = Table()
        summary["Trace_ID"] = _identity_column(
            "Trace_ID", [row["Trace_ID"] for row in summaries], table["Trace_ID"]
        )
        for name in _SUMMARY:
            dtype = (
                str
                if name in {"terminal_state", "stop_reason"}
                else int if name.startswith("n_") else float
            )
            summary[name] = np.asarray([row[name] for row in summaries], dtype=dtype)
        outputs[args.summary_output] = summary
    if args.diagnostics_output:
        diagnostic = Table()
        sources = {
            "Trace_ID": table["Trace_ID"],
            "candidate_spot_id": table["Spot_ID"],
            "candidate_barcode": table[barcode_name],
        }
        for name in (
            "Trace_ID",
            *(field.name for field in fields(TraceCurationIteration)),
        ):
            values = [row[name] for row in diagnostics]
            if name in sources:
                diagnostic[name] = _identity_column(name, values, sources[name])
            else:
                dtype = (
                    str
                    if name == "stop_reason"
                    else (
                        bool
                        if name
                        in {
                            "globally_abnormal",
                            "top5_tiebreak_used",
                            "removal_accepted",
                        }
                        else (
                            int
                            if name == "iteration" or name.startswith("n_")
                            else float
                        )
                    )
                )
                diagnostic[name] = (
                    MaskedColumn(np.asarray(values, dtype=dtype))
                    if name == "stop_reason"
                    else np.asarray(values, dtype=dtype)
                )
        outputs[args.diagnostics_output] = diagnostic
    return outputs


def _write_output(table, path):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=Path(path).parent, delete=False) as output:
            temporary = Path(output.name)
        for column in table.itercols():
            if isinstance(column, MaskedColumn):
                column.info.serialize_method["ecsv"] = "data_mask"
        save_table_to_ecsv(table, temporary)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _print_model(model):
    for field in fields(model):
        if field.name not in {"reference_distributions", "genomic_positions"}:
            print(f"{field.name}: {getattr(model, field.name)}")


def main(argv=None):
    """Run the installed trace_curator command without scientific duplication."""
    parser = parse_arguments()
    args = parser.parse_args(argv)
    try:
        _check_paths(args)
        table = _read_input(args.input)
        if args.command == "curate":
            model = load_curator_model(args.model)
        else:
            calibration = (
                _read_input(args.calibration_input) if args.calibration_input else None
            )
            model = fit_curator_model(
                table,
                trace_fpr=args.trace_fpr,
                minimum_reference_observations=args.minimum_reference_observations,
                reference_mode=args.reference_mode,
                n_crossfit_folds=args.n_crossfit_folds,
                calibration_traces=calibration,
            )
        outputs = _curate_table(table, model, args) if args.command != "fit" else {}
        if args.command != "curate":
            save_curator_model(model, args.output_model)
            _print_model(model)
        for path, output in outputs.items():
            _write_output(output, path)
    except (ValueError, TypeError, OSError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
