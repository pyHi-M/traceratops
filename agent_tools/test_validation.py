"""Regression tests for detecting catalog drift, without changing analysis code."""
import json
from pathlib import Path
import shutil
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parent))
from validate import load, validate

ROOT = Path(__file__).resolve().parents[1]


def test_repository_catalog_is_consistent():
    assert validate(ROOT) == []


@pytest.fixture
def catalog_copy(tmp_path):
    shutil.copytree(ROOT/'agent_tools',tmp_path/'agent_tools')
    for name in ['traceratops','pyproject.toml','docs','test']:
        (tmp_path/name).symlink_to(ROOT/name,target_is_directory=(ROOT/name).is_dir())
    return tmp_path


@pytest.mark.parametrize('mutation,expected',[
    ('missing_command','inventory'),('extra_command','inventory'),
    ('duplicate_id','Duplicate tool'),('default','argparse specification drift'),
    ('choices','argparse specification drift'),('required','argparse specification drift'),
    ('bad_example','invalid syntax'),('unknown_task','unknown tool'),
    ('unknown_preferred','unknown tool'),('missing_example','unknown example'),
    ('api_signature','API signature drift'),('source_hash','scientific source drift'),
    ('empty_capability','missing scientific capability'),
])
def test_detects_drift(catalog_copy,mutation,expected):
    base=catalog_copy/'agent_tools'
    cat=load(base/'catalog.yaml'); tasks=load(base/'tasks.yaml'); examples=load(base/'examples/examples.yaml')
    if mutation=='missing_command':cat['cli_tools'].pop()
    if mutation=='extra_command':
        extra=dict(cat['cli_tools'][0],id='cli.invented',command='invented');cat['cli_tools'].append(extra)
    if mutation=='duplicate_id':cat['cli_tools'][1]['id']=cat['cli_tools'][0]['id']
    if mutation=='default':cat['cli_tools'][0]['parser']['arguments'][1]['default']='invented'
    if mutation=='choices':cat['cli_tools'][0]['parser']['arguments'][0]['choices']=['invented']
    if mutation=='required':cat['cli_tools'][0]['parser']['arguments'][0]['required']=True
    if mutation=='bad_example':examples['examples'][0]['argv']+=['--invented-option']
    if mutation=='unknown_task':tasks['tasks'][0]['tools']=['cli.invented']
    if mutation=='unknown_preferred':tasks['tasks'][0]['preferred_tool']='cli.invented'
    if mutation=='missing_example':cat['cli_tools'][0]['examples']=['example.invented']
    if mutation=='api_signature':cat['python_apis'][0]['signature']='invented()'
    if mutation=='source_hash':
        key=next(iter(cat['cli_tools'][0]['dependency_source_sha256']))
        cat['cli_tools'][0]['dependency_source_sha256'][key]='incorrect'
    if mutation=='empty_capability':cat['cli_tools'][0]['scientific_purpose']=''
    for name,data in [('catalog.yaml',cat),('tasks.yaml',tasks),('examples/examples.yaml',examples)]:
        (base/name).write_text(json.dumps(data))
    assert any(expected in error for error in validate(catalog_copy))


def test_duplicate_mapping_keys_rejected(tmp_path):
    path=tmp_path/'duplicate.yaml';path.write_text('{"same":1,"same":2}')
    with pytest.raises(ValueError,match='Duplicate key'):
        load(path)
