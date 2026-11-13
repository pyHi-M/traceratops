"""Inspect real argparse objects, never formatted --help output."""
import argparse
import importlib
import math


def portable(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [portable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): portable(v) for k, v in value.items()}
    return str(value)


def parser_for(entry_point):
    module, _ = entry_point.split(':')
    return importlib.import_module(module).parse_arguments()


def parser_spec(parser):
    arguments = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        arguments.append({
            'flags': action.option_strings, 'dest': action.dest,
            'required': action.required, 'default': portable(action.default),
            'choices': portable(action.choices), 'nargs': action.nargs,
            'type': getattr(action.type, '__name__', None),
            'action': type(action).__name__, 'const': portable(action.const),
            'help': action.help,
        })
    return {'arguments': arguments, 'mutually_exclusive_groups': [
        {'required': group.required, 'members': [a.dest for a in group._group_actions]}
        for group in parser._mutually_exclusive_groups]}
