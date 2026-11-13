#!/usr/bin/env python3
"""Validate catalog against checkout and installed traceratops. Python 3.11+ stdlib only.

Run from any directory: python agent_tools/validate.py
The .yaml files use JSON syntax (a YAML 1.2 subset), so no YAML dependency is needed.
"""
import argparse
import ast
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inspect_cli import parser_for, parser_spec


def load(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f'Duplicate key: {key}')
            result[key] = value
        return result
    return json.loads(path.read_text(), object_pairs_hook=unique)


def validate(root, installed=True):
    base = root / 'agent_tools'
    cat, tasks, examples = (load(base / name) for name in
                            ('catalog.yaml', 'tasks.yaml', 'examples/examples.yaml'))
    scripts = tomllib.loads((root / 'pyproject.toml').read_text())['project']['scripts']
    errors = []
    def check(condition, message):
        if not condition:
            errors.append(message)
    cli, apis = cat['cli_tools'], cat['python_apis']
    ids = [t['id'] for t in cli + apis]
    check(len(ids) == len(set(ids)), 'Duplicate tool identifiers')
    check(len(cli) == len({t['command'] for t in cli}), 'Duplicate command definitions')
    check({t['command']: t['entry_point'] for t in cli} == scripts,
          'Catalog command inventory differs from project.scripts')
    installed_scripts = {}
    if installed:
        dist = importlib.metadata.distribution('traceratops')
        installed_scripts = {e.name: e.value for e in dist.entry_points if e.group == 'console_scripts'}
        check(installed_scripts == scripts, 'Installed entry points differ from checkout')
    examples_by_id = {e['id']: e for e in examples['examples']}
    check(len(examples_by_id) == len(examples['examples']), 'Duplicate example IDs')
    all_ids = set(ids)
    for tool in cli:
        tid = tool['id']
        check(bool(tool['scientific_purpose'].strip()) and bool(tool['operations']), f'{tid}: missing scientific capability')
        for ref, expected_hash in tool.get('dependency_source_sha256', {}).items():
            check(hashlib.sha256((root/ref).read_bytes()).hexdigest() == expected_hash,
                  f'{tid}: scientific source drift: {ref}')
        parser = parser_for(tool['entry_point'])
        check(parser_spec(parser) == tool['parser'], f'{tid}: argparse specification drift')
        module_path = Path(importlib.util.find_spec(tool['entry_point'].split(':')[0]).origin).resolve()
        check(module_path == (root / tool['entry_point'].split(':')[0].replace('.', '/')).with_suffix('.py').resolve(), f'{tid}: imported module is outside checkout: {module_path}')
        if tool.get('alias_of'):
            check(tool['alias_of'] in all_ids, f'{tid}: unknown alias target')
        check(all(i['schema_ref'] in cat['input_schemas'] for i in tool['inputs']), f'{tid}: unknown input schema')
        check(bool(tool['examples']), f'{tid}: no example')
        for eid in tool['examples']:
            check(eid in examples_by_id, f'{tid}: unknown example {eid}')
        for ref in tool['references']:
            check((root / ref.split(':')[0]).is_file(), f'{tid}: broken reference {ref}')
    by_command = {t['command']: t for t in cli}
    for example in examples['examples']:
        check(example['tool_id'] in all_ids, f"Unknown example tool: {example['tool_id']}")
        command, *args = example['argv']
        check(command in by_command, f'Unknown example executable: {command}')
        if command not in by_command:
            continue
        parser = parser_for(by_command[command]['entry_point'])
        parser.allow_abbrev = False  # examples must use exact verified option spellings
        try:
            parser.parse_args(args)
        except SystemExit as exc:
            check(False, f"{example['id']}: invalid syntax ({exc.code})")
    tids = [t['id'] for t in tasks['tasks']]
    check(len(tids) == len(set(tids)), 'Duplicate task IDs')
    for task in tasks['tasks']:
        for tid in task['tools'] + task.get('sequence', []) + [task['preferred_tool']]:
            check(tid in all_ids, f"{task['id']}: unknown tool {tid}")
    task_tools = {tid for task in tasks['tasks'] for tid in task['tools']}
    check(all(t['id'] in task_tools or t.get('alias_of') for t in cli), 'Canonical CLI missing from task index')
    surface = load(base/'python_surface.yaml')['callables']
    actual_surface = set()
    for source in (root/'traceratops').rglob('*.py'):
        module = str(source.relative_to(root).with_suffix('')).replace('/', '.')
        for node in ast.parse(source.read_text()).body:
            names = [(node.name+'.'+m.name, m) for m in node.body if isinstance(m, ast.FunctionDef)] if isinstance(node, ast.ClassDef) else [(node.name,node)] if isinstance(node, ast.FunctionDef) else []
            actual_surface.update((module,name) for name,n in names if not n.name.startswith('_'))
    check(actual_surface == {(s['module'],s['callable']) for s in surface}, 'Python public-name inventory drift; classify additions/removals')
    check({s['catalog_id'] for s in surface if s['catalog_id']} == {a['id'] for a in apis}, 'Python inclusion ledger differs from catalog')
    for api in apis:
        path = root / api['source']
        check(hashlib.sha256(path.read_bytes()).hexdigest() == api['source_sha256'], f"Python source drift: {api['id']}")
        check(path.is_file(), f"Missing API source {path}")
        nodes = ast.parse(path.read_text()).body
        node = None
        for part in api['callable'].split('.'):
            node = next((n for n in nodes if getattr(n, 'name', None) == part), None)
            if node is None:
                break
            nodes = getattr(node, 'body', [])
        check(node is not None, f"API no longer exists: {api['id']}")
        if isinstance(node, ast.FunctionDef):
            signature = f'{node.name}({ast.unparse(node.args)})'
            if node.returns:
                signature += ' -> ' + ast.unparse(node.returns)
            check(signature == api['signature'], f"API signature drift: {api['id']}")
            check(hashlib.sha256((ast.get_docstring(node) or '').encode()).hexdigest() == api['source_docstring_sha256'], f"API docstring drift: {api['id']}")
    # Scientific routing regression checks: these must stay distinct.
    routing = {t['id']: t for t in tasks['tasks']}
    check(routing['unique_barcode_distribution']['preferred_tool'] == 'cli.trace_analyzer', 'Barcode-distribution routing regression')
    check(routing['dataset_counts']['preferred_tool'] == 'cli.trace_stats', 'Dataset-count routing regression')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--source-only', action='store_true', help='Skip installed distribution metadata check; imports still require runtime dependencies')
    args = parser.parse_args()
    errors = validate(args.root, installed=not args.source_only)
    for error in errors:
        print('ERROR:', error)
    if errors:
        return 1
    print('PASS: CLI inventory, installed entry points, parsers, examples, references, tasks and Python signatures/docstrings')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
