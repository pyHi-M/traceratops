import pytest
from astropy.table import Table

from traceratops.core.localization_table import LocalizationTable


def test_plot_distribution_fluxes_skips_legacy_pyhim_table(tmp_path, capsys):
    legacy_table = Table(
        rows=[
            (
                "1",
                1,
                1,
                1,
                1,
                0.0,
                0.0,
                0.0,
                0.5,
                0.6,
                0.7,
                3,
                1.0,
                10.0,
                20.0,
                0.1,
            )
        ],
        names=(
            "Buid",
            "ROI #",
            "CellID #",
            "Barcode #",
            "id",
            "zcentroid",
            "xcentroid",
            "ycentroid",
            "sharpness",
            "roundness1",
            "roundness2",
            "npix",
            "sky",
            "peak",
            "flux",
            "mag",
        ),
    )
    output_path = tmp_path / "legacy_distribution.png"

    LocalizationTable().plot_distribution_fluxes(legacy_table, [str(output_path)])

    assert "Please update pyHiM" in capsys.readouterr().out
    assert not output_path.exists()


def _localization_table():
    return Table(
        rows=[
            (
                "1",
                1,
                1,
                1,
                "spot-1",
                0.0,
                0.0,
                0.0,
                5.0,
                0.5,
                0.1,
                3,
                1,
                10.0,
                20.0,
                0.9,
            )
        ],
        names=(
            "Buid",
            "ROI #",
            "CellID #",
            "Barcode #",
            "id",
            "zcentroid",
            "xcentroid",
            "ycentroid",
            "snr",
            "spot_pixel_percentage",
            "skew",
            "patch_size",
            "object_class",
            "mean_intensity",
            "flux",
            "roundness",
        ),
    )


def test_load_dat_warns_about_legacy_extension(tmp_path, capsys):
    path = tmp_path / "localizations_3D_barcode.dat"
    _localization_table().write(path, format="ascii.ecsv")

    with pytest.warns(FutureWarning, match=".dat extension are deprecated"):
        table, barcodes = LocalizationTable().load(str(path))

    assert len(table) == 1
    assert list(barcodes) == [1]
    assert (
        "will be discontinued in a future traceratops release"
        in capsys.readouterr().out
    )


def test_save_dat_extension_writes_ecsv_instead(tmp_path):
    requested_path = tmp_path / "merged_localizations.dat"
    actual_path = tmp_path / "merged_localizations.ecsv"

    LocalizationTable().save(str(requested_path), _localization_table())

    assert actual_path.exists()
    assert not requested_path.exists()
