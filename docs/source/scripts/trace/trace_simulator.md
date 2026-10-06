# trace_simulator

`trace_simulator` generates single-polymer traces and complete benchmark ground
truth. Polymer bonds are sampled from the configured normal bond-length model,
the walk is confined to a spherical nucleus, barcodes are detected independently,
and the existing XY/Z localization-error parameters are applied to detections.

## Controlled mislocalization benchmark

Each trace independently enters the corruption arm with probability
`--mislocalization-fraction`. Target selection occurs **after detection** and is
uniform among detected barcodes after excluding
`--mislocalization-edge-margin` barcode indices at both genomic endpoints. For
25 barcodes and margin 2, the eligible barcode identities are 3 through 23.

Exactly one eligible spot is translated by

```
post = pre + --mislocalization-displacement * isotropic_unit_direction
```

Thus displacement is a fixed radial distance in microns, not a Gaussian scale.
A zero displacement still records the selected spot as mislocalized. If a trace
enters the arm but has no eligible detected barcode, it is left unchanged and all
its trace-level corruption summary records `selected_for_corruption=True`,
`corruption_target_available=False`, and status
`no_eligible_detected_barcode`; another trace is not substituted.

The observed ECSV retains one row at most per barcode. Its detection ground truth
preserves the benchmark-runner schema: final coordinates are `Observed_x/y/z`,
pre-perturbation coordinates are `Uncorrupted_Observed_x/y/z`, and corruption is
described by `selected_for_corruption`, `is_corrupted`,
`injected_displacement_um`, `realized_displacement_um`, and the distinct
`realized_error_from_true_um`. A zero-distance target is selected but is not
corrupted. More explicit trace-curator names are included as aliases.

The separate polymer ground truth contains every latent barcode, including
undetected ones. The seed controls polymers, detections, arm assignment, target
choice, and direction. Separate deterministic random streams ensure that runs
with the same seed and detection efficiency select the same traces, barcodes,
and directions at every displacement magnitude; only the amplitude changes.

## Example

```bash
trace_simulator --output small_benchmark.ecsv --n-traces 20 --n-barcodes 25 \
  --detection-efficiency 0.5 --localization-error-xy 0.02 \
  --localization-error-z 0.02 --bond-length 0.15 --bond-length-sd 0.02 \
  --nucleus-diameter 1.5 --mislocalization-fraction 0.2 \
  --mislocalization-displacement 0.2 --mislocalization-edge-margin 2 --seed 123
```

This writes the runner-compatible `small_benchmark.ecsv`,
`small_benchmark.ground_truth.ecsv`,
`small_benchmark.polymer_ground_truth.ecsv`, and
`small_benchmark.arguments.yaml`, plus
`corruption.summary.ecsv` for trace-level outcomes. Supplying an
output named `simulated.ecsv` therefore produces the standard benchmark names.
