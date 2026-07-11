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
