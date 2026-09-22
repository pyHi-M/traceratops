import numpy as np
import matplotlib.pyplot as plt
from traceratops.core.localization_table import LocalizationTable

lt = LocalizationTable()
barcode_map, _ = lt.load("/home/jb/Desktop/Aileen/experiment_231_analysis/024/localize_3d/data/localizations_3D_barcode.ecsv")
barcode_id = barcode_map["Barcode #"]
snr = barcode_map["snr"]

# unique_barcodes = np.sort(np.unique(barcode_id))
# print(f"unique barcodes in file: {len(unique_barcodes)}")
#
# print(type(barcode_map["snr"]))
# mask = getattr(barcode_map["snr"], "mask", None)
# if mask is not None:
#     print("masked entries:", mask.sum(), "/", len(mask))
#     print("fill_value:", barcode_map["snr"].fill_value)


unique_barcodes = np.sort(np.unique(barcode_id))[:10]
positions = np.arange(len(unique_barcodes))
snr_by_barcode = [
    np.asarray(snr[barcode_id == bc], dtype=float) for bc in unique_barcodes
]
snr_by_barcode = [vals[np.isfinite(vals)] for vals in snr_by_barcode]
# snr_by_barcode = [snr[barcode_id == bc] for bc in unique_barcodes]

for bc, vals in zip(unique_barcodes, snr_by_barcode):
    arr = np.ma.asarray(vals)  # keeps mask info if present, unlike plain np.asarray
    print(bc, "n =", len(arr), "masked =", np.ma.count_masked(arr),
          "min/max =", arr.min(), arr.max())

fig, ax = plt.subplots(figsize=(10, 6))
parts = ax.violinplot(snr_by_barcode, positions=positions, widths=0.8,
                       showmeans=False, showmedians=True, showextrema=True)

print("bodies returned:", len(parts["bodies"]), "  expected:", len(snr_by_barcode))

ax.set_xticks(positions)
ax.set_xticklabels(unique_barcodes)
plt.savefig("/home/jb/Desktop/Aileen/experiment_231_analysis/024/debug_violin_first10.png", dpi=150)