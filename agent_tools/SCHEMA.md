# Catalog contract, version 1.0

Files ending in `.yaml` use JSON syntax, a YAML 1.2 subset. Read with Python's
standard `json` module or a YAML 1.2 reader. No new runtime dependency is required.
UTF-8, no duplicate mapping keys, IDs stable across releases. Changes to scientific
behavior require review, not only parser snapshot refresh.

`catalog.yaml`:

- `repository`: upstream URL, reviewed branch and exact commit. A snapshot is not a promise about a later `dev` revision.
- `input_schemas`: reusable format/schema contracts, referenced by `inputs[].schema_ref`.
- `cli_tools[]`: one entry per installed command name, including deprecated aliases. `alias_of` links to the canonical command; prefer canonical names.
- `id`: `cli.<exact executable>` or `python.<module beneath traceratops>.<qualified callable>`.
- `description`, `scientific_purpose`, `operations`, `concepts`, `select_for`: discovery fields. The task index adds synonyms and distinctions.
- `parser.arguments[]`: direct argparse introspection. `flags` includes all aliases; `dest`, `action`, `type`, `nargs`, `required`, `const`, `choices`, `default`, `help` preserve actual parser semantics. `type: null` means no argparse conversion. `choices: null` means the parser does not constrain values, not that arbitrary values work scientifically.
- `parser.mutually_exclusive_groups[]`: member destinations and whether one is required. A group requirement is not an individual argument requirement.
- `runtime_requirements`: inputs enforced downstream or needed scientifically, even when argparse marks them optional.
- `runtime_default_rules`: literal fallback statements from simple `if args.<dest> ... else ...` branches in `create_dict_args`. These are source evidence, not code for an agent to execute. They are intentionally separate from parser defaults. They are not an exhaustive symbolic interpretation of runtime code; inspect source references for complex conditions. Infinity is encoded as string `inf`.
- `outputs[]`: file patterns, formats and meanings. `modifies_input_files: false` describes ordinary execution with derived output paths, not a guarantee against an explicitly chosen output path overwriting the input. `write_caveat` is mandatory.
- `examples[]`: IDs in `examples/examples.yaml`.
- `references[]`: paths relative to the repository root; optional `:line` suffix.
- `verification`: source review, parser inspection and example execution are separate evidence levels. Executing one example does not validate every option or the scientific estimator.
- `python_apis[]`: separate callable entries with source signature (including `self`), optional class initialization signature, literal documented return type, return annotation, parameter documentation, input contract and caveats. Null return documentation means undocumented, not `None`. Source-docstring hashes detect drift without duplicating full docstrings. Source documentation does not establish a package-wide public stability promise.

`tasks.yaml`: tasks have stable IDs, request text, preferred tool ID, alternate
synonyms, `distinguish_from`, all relevant `tools`, and an ordered `sequence` for
workflows. Sequence entries use the output of the preceding operation: agents must
supply concrete paths and choose scientifically appropriate parameters. A synonym
is a retrieval cue, not a claim of exact scientific equivalence.

`examples/examples.yaml`: argv arrays are passed directly to subprocess; stdin is
newline-separated filenames, not file content. Run in the directory prepared by
`make_inputs.py`. `expected_outputs` includes globs or `stdout: ` assertions.
Every example documents objective, setup, requirements, command, outputs, checks,
limitations and verification status. Outputs stay in an external work directory.

`examples/verification.json` records observed exit codes, checked outputs and
runtime versions, without generated images or datasets. Initial failed trials
remain in the evidence summary where scientifically relevant.

Maintenance:

1. Inspect `pyproject.toml`, public executable scripts, source docstrings, documentation and tests; update `INVENTORY.md`.
2. Update semantic fields by source review. Compare the installed parser with `inspect_cli.parser_spec(inspect_cli.parser_for(entry_point))`. Do not rewrite scientific text from help alone.
3. Update tasks and examples, then run `python agent_tools/validate.py` and `python -m pytest agent_tools/test_validation.py`.
4. Run representative commands with `python agent_tools/examples/run_examples.py --work-dir /tmp/traceratops-examples`. Review logs and output meanings; record evidence and update `COVERAGE.md`.
5. Review new Python analytical names against inclusion policy. Public source docstrings count as documentation; private helpers and parser/I/O/reporting infrastructure are not recommended analysis APIs. Functions without documentation are listed as gaps, not invented interfaces.
