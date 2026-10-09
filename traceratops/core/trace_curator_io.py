"""Transparent JSON/NPZ persistence for fitted trace-curator models.

Identity support follows scalar ECSV conventions: strings, integers, finite
floats, and booleans. Arrays never use object dtype or pickle serialization.
"""

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any
from zipfile import BadZipFile

import numpy as np

from traceratops.core.trace_curator import TraceCuratorModel

_MODEL_FIELDS = (
    "global_threshold",
    "reference_mode",
    "minimum_reference_observations",
    "trace_fpr",
    "n_reference_traces",
    "n_reference_pairs",
    "n_calibration_traces",
    "n_calibration_unscoreable",
    "calibration_source",
    "n_crossfit_folds",
)


def model_paths(prefix):
    """Return NPZ and JSON paths from a prefix or either model filename."""
    prefix = Path(prefix)
    if prefix.suffix in {".json", ".npz"}:
        prefix = prefix.with_suffix("")
    return Path(str(prefix) + ".npz"), Path(str(prefix) + ".json")


def _identity(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, str):
        kind = "str"
    elif isinstance(value, bool):
        kind = "bool"
    elif isinstance(value, int):
        kind = "int"
    elif isinstance(value, float) and np.isfinite(value):
        kind = "float"
    else:
        raise ValueError(
            "Model identities must be strings, integers, finite floats, or booleans"
        )
    return {"type": kind, "value": value}


def _read_identity(record):
    if not isinstance(record, dict) or set(record) != {"type", "value"}:
        raise ValueError("Invalid typed identity")
    kind, value = record["type"], record["value"]
    types = {"str": str, "bool": bool, "int": int, "float": float}
    if kind not in types or type(value) is not types[kind]:
        raise ValueError("Invalid identity type/value")
    if kind == "float" and isinstance(value, float) and not np.isfinite(value):
        raise ValueError("Nonfinite identity")
    return value


def _checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_curator_model(model, prefix):
    """Save versioned typed metadata and numeric NPZ arrays, replacing outputs.

    Each file is replaced atomically. The checksum detects mismatched file
    pairs (the two-file replacement itself is not a filesystem transaction).
    """
    npz_path, json_path = model_paths(prefix)
    arrays, keys = {}, []
    for offset, (key, sample) in enumerate(model.reference_distributions.items()):
        name = f"reference_{offset}"
        encoded: Any
        if model.reference_mode == "separation":
            encoded = float(key)
        else:
            encoded = sorted(
                (_identity(barcode) for barcode in key),
                key=lambda x: json.dumps(x, sort_keys=True),
            )
            if len(encoded) != 2:
                raise ValueError(
                    "Barcode-pair reference keys must contain two identities"
                )
        keys.append({"array": name, "key": encoded})
        arrays[name] = sample
    metadata = {
        "format_version": 1,
        "model": {
            name: (
                value.item()
                if isinstance(value := getattr(model, name), np.generic)
                else value
            )
            for name in _MODEL_FIELDS
        },
        "chromosome": None if model.chromosome is None else _identity(model.chromosome),
        "genomic_positions": [
            {"barcode": _identity(barcode), "position": float(position)}
            for barcode, position in model.genomic_positions.items()
        ],
        "references": keys,
    }
    temporary = []
    try:
        with tempfile.NamedTemporaryFile(dir=npz_path.parent, delete=False) as output:
            temporary.append(Path(output.name))
            np.savez(output, **arrays)
        metadata["npz_sha256"] = _checksum(temporary[0])
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=json_path.parent, delete=False
        ) as output:
            temporary.append(Path(output.name))
            json.dump(metadata, output, ensure_ascii=False, allow_nan=False, indent=2)
            output.write("\n")
        os.replace(temporary[0], npz_path)
        os.replace(temporary[1], json_path)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)


def load_curator_model(prefix):
    """Load a checked, pickle-free model; reconstruct immutable references."""
    npz_path, json_path = model_paths(prefix)
    try:
        with open(json_path, encoding="utf-8") as source:
            metadata = json.load(source)
        if (
            set(metadata)
            != {
                "format_version",
                "model",
                "chromosome",
                "genomic_positions",
                "references",
                "npz_sha256",
            }
            or type(metadata["format_version"]) is not int
            or metadata["format_version"] != 1
        ):
            raise ValueError("Unsupported or malformed model metadata version")
        parameters = metadata["model"]
        if set(parameters) != set(_MODEL_FIELDS):
            raise ValueError("Missing or unknown model parameters")
        for name in (
            "minimum_reference_observations",
            "n_reference_traces",
            "n_reference_pairs",
            "n_calibration_traces",
            "n_calibration_unscoreable",
        ):
            if type(parameters[name]) is not int or parameters[name] < 0:
                raise ValueError(f"Invalid model count: {name}")
        folds = parameters["n_crossfit_folds"]
        if parameters["calibration_source"] == "crossfit":
            if (
                type(folds) is not int
                or not 2 <= folds <= parameters["n_reference_traces"]
            ):
                raise ValueError("Invalid cross-fitting metadata")
        elif parameters["calibration_source"] != "separate_traces" or folds is not None:
            raise ValueError("Invalid calibration strategy")
        positions = {}
        for record in metadata["genomic_positions"]:
            barcode = _read_identity(record["barcode"])
            if barcode in positions:
                raise ValueError("Duplicated genomic mapping identity")
            positions[barcode] = record["position"]
        if _checksum(npz_path) != metadata["npz_sha256"]:
            raise ValueError(
                "NPZ checksum mismatch: model files are damaged or mismatched"
            )
        references = {}
        with np.load(npz_path, allow_pickle=False) as archive:
            names = [record["array"] for record in metadata["references"]]
            if len(set(names)) != len(names) or set(archive.files) != set(names):
                raise ValueError("Reference arrays do not match metadata")
            for record in metadata["references"]:
                key = record["key"]
                if parameters["reference_mode"] == "separation":
                    if (
                        isinstance(key, bool)
                        or not isinstance(key, (int, float))
                        or not np.isfinite(key)
                        or key <= 0
                    ):
                        raise ValueError("Invalid genomic separation key")
                    key = float(key)
                else:
                    key = frozenset(_read_identity(item) for item in key)
                    if len(key) != 2:
                        raise ValueError("Invalid barcode-pair key")
                if key in references:
                    raise ValueError("Duplicated reference key")
                sample = archive[record["array"]]
                if sample.dtype != np.dtype("float64") or sample.ndim != 1:
                    raise ValueError(
                        "Reference arrays must be one-dimensional float64 arrays"
                    )
                references[key] = sample
        chromosome = metadata["chromosome"]
        return TraceCuratorModel(
            reference_distributions=references,
            genomic_positions=positions,
            chromosome=None if chromosome is None else _read_identity(chromosome),
            **parameters,
        )
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        AttributeError,
        EOFError,
        BadZipFile,
    ) as error:
        raise ValueError(f"Cannot load curator model '{prefix}': {error}") from error
