"""Small deterministic synthetic inputs. Scientific outputs are produced by traceratops.

No experimental data or generated figures are committed. Units: µm; genome: synthetic.
"""
import json
from pathlib import Path
import numpy as np
from astropy.table import Table


def make_inputs(folder, make_bootstrap=False):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(17)
    rows = []
    for trace in range(16):
        for barcode in range(1, 6):
            if trace % 4 == barcode % 4 and trace > 3:
                continue
            coord = np.array([0.5 + barcode * 0.04, 0.6, 0.7]) + rng.normal(0, 0.025, 3)
            rows.append((f'spot-{trace}-{barcode}', f'trace-{trace:02}', *coord,
                         'chrSynthetic', barcode * 10000, barcode * 10000 + 100,
                         trace % 2 + 1, 1, barcode, 'positive' if trace % 2 else 'negative'))
    table = Table(rows=rows, names=['Spot_ID','Trace_ID','x','y','z','Chrom','Chrom_Start','Chrom_End','ROI #','Mask_id','Barcode #','label'])
    table.meta['comments'] = ['xyz_unit=micron', 'genome_assembly=synthetic']
    table.write(folder / 'trace.ecsv', format='ascii.ecsv', overwrite=True)
    table2 = table.copy()
    table2['Trace_ID'] = [f'second-{v}' for v in table2['Trace_ID']]
    table2['Spot_ID'] = [f'second-{v}' for v in table2['Spot_ID']]
    table2['x'] = table2['x'] * 1.1
    table2.write(folder / 'trace2.ecsv', format='ascii.ecsv', overwrite=True)
    loc = Table({'Buid':table['Spot_ID'], 'ROI #':table['ROI #'], 'CellID #':table['Mask_id'],
                 'Barcode #':table['Barcode #'], 'id':np.arange(len(table)),
                 'xcentroid':table['x'], 'ycentroid':table['y'], 'zcentroid':table['z'],
                 'snr':rng.uniform(10,20,len(table)), 'spot_pixel_percentage':rng.uniform(.1,.4,len(table)),
                 'skew':rng.uniform(.1,.8,len(table)), 'patch_size':np.full(len(table),9),
                 'object_class':np.ones(len(table),dtype=int), 'mean_intensity':rng.uniform(1000,2000,len(table)),
                 'flux':rng.uniform(500,1000,len(table)), 'roundness':rng.uniform(.1,.6,len(table))})
    loc.meta['comments'] = ['synthetic fixture']
    loc.write(folder / 'localizations.ecsv', format='ascii.ecsv', overwrite=True)
    loc2 = loc.copy(); loc2['Buid'] = [f'second-{v}' for v in loc2['Buid']]
    loc2.write(folder / 'localizations2.ecsv', format='ascii.ecsv', overwrite=True)
    (folder / 'barcodes.bed').write_text(''.join(f'chrSynthetic\t{b*10000}\t{b*10000+100}\t{b}\n' for b in range(1,6)))
    (folder / 'metadata.json').write_text(json.dumps({'genome_assembly':'synthetic','experimenter_name':'Synthetic test','experimenter_contact':'test@example.org'}))
    np.savetxt(folder / 'barcodes.txt', np.arange(1,6), fmt='%d')
    np.save(folder / 'mask.npy', np.ones((32,32),dtype=int))
    # Reuse repository's small distance-stack fixture (4 barcodes x 4 traces).
    root = Path(__file__).resolve().parents[2]
    matrix = np.load(root / 'test/data/plot_him_matrix/IN/bc_4_cells_4.npy')
    np.save(folder / 'matrix.npy', matrix)
    np.save(folder / 'matrix2.npy', matrix * 1.1)
    np.savetxt(folder / 'matrix_barcodes.txt', [1,2,3,4], fmt='%d')
    # Replot input comes from the actual public analysis function, not fabricated results.
    from traceratops.trace_3way_coloc import compute_threeway_colocalization
    frequencies = compute_threeway_colocalization(table, 1, .2)
    threeway = np.zeros((4,4))
    for (a,b), value in frequencies.items():
        threeway[a-2,b-2] = threeway[b-2,a-2] = value
    np.save(folder / 'threeway_anchor_1.npy', threeway)
    for name in ['roi1','roi2']:
        dest = folder / 'collection' / name
        dest.mkdir(parents=True,exist_ok=True)
        (dest / 'sample.ecsv').write_bytes((folder/'trace.ecsv').read_bytes())

    # Dense-enough stack for bootstrap: use the package implementation itself.
    from traceratops.core.build_matrix import BuildMatrix
    from traceratops.core.chromatin_trace_table import ChromatinTraceTable
    if make_bootstrap:
        builder = BuildMatrix({})
        builder.trace_table = ChromatinTraceTable()
        builder.trace_table.data = table
        builder.build_distance_matrix()
        np.save(folder / 'bootstrap.npy', builder.sc_matrix)
