#!/usr/bin/env python3
"""Run actual installed executables in isolated folders; retain logs outside deliverables.

python agent_tools/examples/run_examples.py --work-dir /tmp/traceratops-examples
python agent_tools/examples/run_examples.py --work-dir /tmp/traceratops-examples --command trace_stats
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_inputs import make_inputs


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(name):
    candidate = Path(sysconfig.get_path('scripts')) / (name + ('.exe' if os.name == 'nt' else ''))
    if not candidate.is_file():
        raise FileNotFoundError(f'Installed executable missing: {candidate}')
    return str(candidate)


def semantic_checks(command, folder, stdout):
    import numpy as np
    from astropy.table import Table
    read = lambda name: Table.read(folder/name, format='ascii.ecsv')
    if command == 'trace_stats':
        assert 'Number of unique chromatin traces: 16' in stdout
        assert 'Number of unique barcodes: 5' in stdout
        assert 'Number of unique ROIs: 2' in stdout
    if command == 'trace_filter':
        result = read('trace_filtered.ecsv')
        for trace in result.group_by('Trace_ID').groups:
            assert len(set(trace['Barcode #'])) >= 4
    if command == 'trace_merge':
        assert len(read('merged.ecsv')) == len(read('trace.ecsv')) + len(read('trace2.ecsv'))
    if command == 'localization_merge':
        assert len(read('merged_localizations.ecsv')) == 2 * len(read('localizations.ecsv'))
    if command == 'trace_split_labels':
        positive, negative, original = read('trace_positive.ecsv'), read('trace_not:positive.ecsv'), read('trace.ecsv')
        assert len(positive)+len(negative) == len(original)
        assert set(positive['Spot_ID']).isdisjoint(set(negative['Spot_ID']))
    if command == 'trace_to_matrix':
        matrix = np.load(next((folder/'matrices').glob('*PWDscMatrix.npy')))
        assert matrix.shape == (5,5,16)
        np.testing.assert_allclose(matrix,matrix.transpose(1,0,2),equal_nan=True)
        assert np.isnan(matrix[np.arange(5),np.arange(5),:]).all()
        counts = np.load(next((folder/'matrices').glob('*Nmatrix.npy')))
        np.testing.assert_array_equal(counts,np.isfinite(matrix).sum(axis=2))
    if command == 'trace_sorter':
        original, result = read('trace.ecsv'), read('trace_sorted.ecsv')
        assert set(original['Spot_ID']) == set(result['Spot_ID'])
        groups = list(dict.fromkeys(result['Trace_ID']))
        sizes = [sum(result['Trace_ID']==t) for t in groups]
        assert sizes == sorted(sizes,reverse=True)
    if command == 'trace_import_from_fofct':
        original, result = read('trace.ecsv'), read('restored.ecsv')
        for col in ['x','y','z','Barcode #','ROI #']:
            np.testing.assert_allclose(original[col],result[col])
    if command == 'trace_assign_mask':
        assert all('maskpositive' in str(v) for v in read('trace_maskpositive.ecsv')['label'])
    if command in ['trace_3way_coloc','plot_3way']:
        matrix = np.load(folder/'threeway_trace_anchor_1.npy')
        assert matrix.shape == (4,4)
        assert np.isfinite(matrix).all() and (matrix>=0).all() and (matrix<=1).all()
    if command == 'trace_pearsons':
        matrix = np.load(folder/'trace_correlation_matrix.npy')
        assert matrix.shape == (2,2)
        np.testing.assert_allclose(np.diag(matrix),[1,1])
    if command == 'collect_files':
        assert digest(folder/'collection/roi1/sample.ecsv') == digest(folder/'collected/sample_roi1.ecsv')


def run_one(example, work):
    command = example['argv'][0]
    folder = work / command
    # Never delete user directories. Each run creates a fresh unique directory.
    folder.mkdir()
    make_inputs(folder, make_bootstrap=command == "plot_bootstrapping")
    tracked = {p:p.read_bytes() for p in folder.rglob('*') if p.is_file()}
    env = dict(os.environ, MPLBACKEND='Agg', MPLCONFIGDIR=str(work/'mpl'), XDG_CACHE_HOME=str(work/'cache'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = {'id':example['id'],'command':command,'argv':example['argv'],'stdin':example['stdin'],'status':'failed'}
    started = time.monotonic()
    try:
        if command == 'trace_import_from_fofct':
            prep = subprocess.run([executable('trace_export_to_fofct'),'--input','trace.ecsv','--bed_file','barcodes.bed','--json_file','metadata.json','--output_file','export.csv'],cwd=folder,env=env,capture_output=True,text=True,timeout=90)
            assert prep.returncode == 0, prep.stderr
        proc = subprocess.run([executable(command),*example['argv'][1:]],cwd=folder,env=env,input=example['stdin'] or '',capture_output=True,text=True,timeout=120)
        (folder/'stdout.log').write_text(proc.stdout)
        (folder/'stderr.log').write_text(proc.stderr)
        result['returncode'] = proc.returncode
        assert proc.returncode == 0, proc.stderr[-1500:]
        checked = []
        for pattern in example['expected_outputs']:
            if pattern.startswith('stdout: '):
                assert pattern[8:] in proc.stdout
            else:
                paths = list(folder.glob(pattern))
                assert paths, f'No output matches {pattern}'
                assert all(p.stat().st_size>0 for p in paths), f'Empty output: {pattern}'
                checked += [str(p.relative_to(folder)) for p in paths]
        semantic_checks(command,folder,proc.stdout)
        assert all(p.read_bytes()==value for p,value in tracked.items()), 'Input fixture changed'
        result.update(status='passed',outputs=checked,checks='exit, nonempty expected artifacts, semantic assertions, unchanged input fixtures')
    except Exception as exc:
        result['error'] = str(exc)
    result['elapsed_seconds'] = round(time.monotonic()-started,2)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--command',action='append')
    parser.add_argument('--jobs',type=int,default=2)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True,exist_ok=True)
    # Subfolder protects previous runs and prevents stale output assertions.
    import tempfile
    work = Path(tempfile.mkdtemp(prefix='run-',dir=work))
    examples = json.loads((Path(__file__).parent/'examples.yaml').read_text())['examples']
    if args.command:
        examples = [e for e in examples if e['argv'][0] in args.command]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda e: run_one(e,work),examples))
    report = {'python':sys.version.split()[0], 'results':results}
    (work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    for r in results:
        print(r['command'],r['status'],r.get('error',''),flush=True)
    print('Evidence:',work/'results.json')
    return int(any(r['status']!='passed' for r in results))

if __name__=='__main__':
    raise SystemExit(main())
