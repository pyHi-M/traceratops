#!/usr/bin/env python3
"""Simulate single-polymer traces, optionally mislocalizing one detected spot."""

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table


@dataclass(frozen=True)
class SimulationConfig:
    """Parameters for single-polymer trace simulation (distances are microns)."""

    n_traces: int = 100
    n_barcodes: int = 25
    detection_efficiency: float = 0.8
    localization_error_xy: float = 0.02
    localization_error_z: float = 0.02
    bond_length: float = 0.15
    bond_length_sd: float = 0.02
    nucleus_diameter: float = 1.5
    mislocalization_fraction: float = 0.0
    mislocalization_displacement: float = 0.0
    mislocalization_edge_margin: int = 0
    seed: int = 0

    def validate(self):
        if self.n_traces < 1 or self.n_barcodes < 1:
            raise ValueError("n_traces and n_barcodes must be positive")
        if not 0 <= self.detection_efficiency <= 1:
            raise ValueError("detection_efficiency must be between 0 and 1")
        if not 0 <= self.mislocalization_fraction <= 1:
            raise ValueError("mislocalization_fraction must be between 0 and 1")
        if self.mislocalization_displacement < 0:
            raise ValueError("mislocalization_displacement must be nonnegative")
        if self.mislocalization_edge_margin < 0:
            raise ValueError("mislocalization_edge_margin must be nonnegative")
        if min(self.localization_error_xy, self.localization_error_z) < 0:
            raise ValueError("localization errors must be nonnegative")
        if self.bond_length < 0 or self.bond_length_sd < 0:
            raise ValueError("bond lengths must be nonnegative")
        if self.nucleus_diameter <= 0:
            raise ValueError("nucleus_diameter must be positive")


def random_unit_vectors(rng, size=1):
    """Draw directions uniformly from the unit sphere using normalized normals."""
    vectors = rng.normal(size=(size, 3))
    norms = np.linalg.norm(vectors, axis=1)
    while np.any(norms == 0):  # vanishingly rare, but avoids an undefined direction
        vectors[norms == 0] = rng.normal(size=(np.sum(norms == 0), 3))
        norms = np.linalg.norm(vectors, axis=1)
    return vectors / norms[:, None]


def _polymer(rng, config):
    positions = np.zeros((config.n_barcodes, 3), dtype=float)
    radius = config.nucleus_diameter / 2
    for index in range(1, config.n_barcodes):
        length = max(0.0, rng.normal(config.bond_length, config.bond_length_sd))
        # Rejection confines the walk without changing an accepted bond length.
        for _ in range(10000):
            candidate = positions[index - 1] + length * random_unit_vectors(rng)[0]
            if np.linalg.norm(candidate) <= radius:
                positions[index] = candidate
                break
        else:
            raise RuntimeError("could not place a polymer bond inside the nucleus")
    return positions


def simulate_traces(config, include_auxiliary=False):
    """Return ``(observed, ground_truth)`` Astropy tables.

    Corruption selection happens after detection. A selected trace with no detected
    barcode inside the configured margins remains unchanged and has
    a ``no_eligible_detected_barcode`` trace-summary status.
    """
    config.validate()
    observed_rows, truth_rows, polymer_rows, summary_rows = [], [], [], []
    margin = config.mislocalization_edge_margin

    for trace_index in range(config.n_traces):
        # Independent streams keep latent polymers and ordinary observation noise
        # unchanged when only corruption parameters are varied.
        streams = np.random.SeedSequence([config.seed, trace_index]).spawn(5)
        polymer_rng, detection_rng, noise_rng, selection_rng, corruption_rng = (
            np.random.default_rng(stream) for stream in streams
        )
        trace_id = f"trace-{trace_index + 1}"
        latent = _polymer(polymer_rng, config)
        detected = detection_rng.random(config.n_barcodes) < config.detection_efficiency
        pre = latent + noise_rng.normal(
            scale=(config.localization_error_xy,) * 2 + (config.localization_error_z,),
            size=latent.shape,
        )
        post = pre.copy()
        selected = bool(selection_rng.random() < config.mislocalization_fraction)
        eligible = np.flatnonzero(
            detected
            & (np.arange(config.n_barcodes) >= margin)
            & (np.arange(config.n_barcodes) < config.n_barcodes - margin)
        )
        target = (
            int(corruption_rng.choice(eligible)) if selected and len(eligible) else None
        )
        if target is not None:
            post[target] += (
                config.mislocalization_displacement
                * random_unit_vectors(corruption_rng)[0]
            )

        summary_rows.append(
            {
                "Input_Trace_ID": trace_id,
                "selected_for_corruption": selected,
                "corruption_target_available": target is not None,
                "selected_Barcode #": target + 1 if target is not None else -1,
                "selected_Spot_ID": (
                    f"{trace_id}-spot-{target + 1}" if target is not None else ""
                ),
                "requested_displacement_um": config.mislocalization_displacement,
                "status": (
                    "not_selected"
                    if not selected
                    else (
                        "applied"
                        if target is not None
                        else "no_eligible_detected_barcode"
                    )
                ),
            }
        )

        for index in range(config.n_barcodes):
            barcode = index + 1
            spot_id = f"{trace_id}-spot-{barcode}" if detected[index] else ""
            is_target = index == target
            realized = (
                float(np.linalg.norm(post[index] - pre[index])) if is_target else 0.0
            )
            group_id = trace_id
            polymer_id = f"{trace_id}-polymer-1"
            polymer_rows.append(
                {
                    "Input_Trace_ID": trace_id,
                    "GroundTruth_Group_ID": group_id,
                    "GroundTruth_Polymer_ID": polymer_id,
                    "Barcode #": barcode,
                    "True_x": latent[index, 0],
                    "True_y": latent[index, 1],
                    "True_z": latent[index, 2],
                }
            )
            common = {
                "Input_Trace_ID": trace_id,
                "Spot_ID": spot_id,
                "GroundTruth_Group_ID": group_id,
                "GroundTruth_Polymer_ID": polymer_id,
                "Detection_Source": "polymer",
                "Barcode #": barcode,
                "Occurrence_Index": 0,
                "True_x": latent[index, 0],
                "True_y": latent[index, 1],
                "True_z": latent[index, 2],
                "Observed_x": post[index, 0],
                "Observed_y": post[index, 1],
                "Observed_z": post[index, 2],
                "Uncorrupted_Observed_x": pre[index, 0],
                "Uncorrupted_Observed_y": pre[index, 1],
                "Uncorrupted_Observed_z": pre[index, 2],
                "selected_for_corruption": is_target,
                "is_corrupted": is_target and realized > 0,
                "injected_displacement_um": (
                    config.mislocalization_displacement if is_target else 0.0
                ),
                "realized_displacement_um": realized,
                "realized_error_from_true_um": float(
                    np.linalg.norm(post[index] - latent[index])
                ),
                # More explicit aliases retained for the trace-curator specification.
                "Barcode": barcode,
                "Latent_x": latent[index, 0],
                "Latent_y": latent[index, 1],
                "Latent_z": latent[index, 2],
                "Is_Mislocalized": is_target,
                "Mislocalization_Selected_Trace": selected,
                "Mislocalization_Applied": target is not None,
                "Requested_Displacement_um": (
                    config.mislocalization_displacement if is_target else 0.0
                ),
                "Realized_Displacement_um": realized,
                "PreCorruption_Observed_x": (
                    pre[index, 0] if detected[index] else np.nan
                ),
                "PreCorruption_Observed_y": (
                    pre[index, 1] if detected[index] else np.nan
                ),
                "PreCorruption_Observed_z": (
                    pre[index, 2] if detected[index] else np.nan
                ),
                "PostCorruption_Observed_x": (
                    post[index, 0] if detected[index] else np.nan
                ),
                "PostCorruption_Observed_y": (
                    post[index, 1] if detected[index] else np.nan
                ),
                "PostCorruption_Observed_z": (
                    post[index, 2] if detected[index] else np.nan
                ),
            }
            if detected[index]:
                truth_rows.append(common)
                observed_rows.append(
                    {
                        "Spot_ID": spot_id,
                        "Trace_ID": trace_id,
                        "x": post[index, 0],
                        "y": post[index, 1],
                        "z": post[index, 2],
                        "Barcode #": barcode,
                    }
                )

    observed_names = ("Spot_ID", "Trace_ID", "x", "y", "z", "Barcode #")
    truth_names = tuple(common) if config.n_barcodes else ()
    observed = Table(rows=observed_rows, names=observed_names)
    truth = Table(rows=truth_rows, names=truth_names)
    observed.meta["comments"] = ["xyz_unit=micron", "single-polymer simulation"]
    truth.meta["comments"] = ["xyz_unit=micron", "one row per detected barcode"]
    truth.meta["mislocalization_fraction"] = config.mislocalization_fraction
    truth.meta["mislocalization_displacement_um"] = config.mislocalization_displacement
    truth.meta["mislocalization_edge_margin"] = config.mislocalization_edge_margin
    truth.meta["seed"] = config.seed
    polymer_truth = Table(rows=polymer_rows)
    summary = Table(rows=summary_rows)
    if include_auxiliary:
        return observed, truth, polymer_truth, summary
    return observed, truth


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Observed ECSV output path.")
    parser.add_argument(
        "--ground-truth-output",
        help="Ground-truth ECSV path (default: <output stem>.ground_truth.ecsv).",
    )
    parser.add_argument("--n-traces", type=int, default=100)
    parser.add_argument("--n-barcodes", type=int, default=25)
    parser.add_argument("--detection-efficiency", type=float, default=0.8)
    parser.add_argument(
        "--localization-error-xy",
        type=float,
        default=0.02,
        help="XY localization SD in microns.",
    )
    parser.add_argument(
        "--localization-error-z",
        type=float,
        default=0.02,
        help="Z localization SD in microns.",
    )
    parser.add_argument(
        "--bond-length", type=float, default=0.15, help="Mean bond length in microns."
    )
    parser.add_argument(
        "--bond-length-sd", type=float, default=0.02, help="Bond-length SD in microns."
    )
    parser.add_argument(
        "--nucleus-diameter",
        type=float,
        default=1.5,
        help="Spherical nucleus diameter in microns.",
    )
    parser.add_argument(
        "--mislocalization-fraction",
        type=float,
        default=0.0,
        help="Independent probability that each input trace enters the corruption arm.",
    )
    parser.add_argument(
        "--mislocalization-displacement",
        type=float,
        default=0.0,
        help="Fixed radial displacement in microns (not a Gaussian SD).",
    )
    parser.add_argument(
        "--mislocalization-edge-margin",
        type=int,
        default=0,
        help="Number of barcode indices excluded at each genomic endpoint.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed controlling simulation, detection, target choice, and direction.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    output = Path(args.output)
    ground_truth = (
        Path(args.ground_truth_output)
        if args.ground_truth_output
        else output.with_name(output.stem + ".ground_truth.ecsv")
    )
    values = vars(args).copy()
    values.pop("output")
    values.pop("ground_truth_output")
    config = SimulationConfig(**values)
    observed, truth, polymer_truth, summary = simulate_traces(
        config, include_auxiliary=True
    )
    observed.write(output, format="ascii.ecsv", overwrite=True)
    truth.write(ground_truth, format="ascii.ecsv", overwrite=True)
    polymer_truth.write(
        output.with_name(output.stem + ".polymer_ground_truth.ecsv"),
        format="ascii.ecsv",
        overwrite=True,
    )
    summary.write(
        output.with_name("corruption.summary.ecsv"),
        format="ascii.ecsv",
        overwrite=True,
    )
    with output.with_name(output.stem + ".arguments.yaml").open("w") as stream:
        # JSON is a YAML 1.2 subset and avoids adding a simulator-only dependency.
        json.dump(asdict(config), stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
