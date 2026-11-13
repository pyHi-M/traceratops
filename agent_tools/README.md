# Scientific agent tools for traceratops

Start with [AGENT_GUIDE.md](AGENT_GUIDE.md). Detailed tool records are in
[catalog.yaml](catalog.yaml), scientific routing in [tasks.yaml](tasks.yaml), and
executed examples in [examples/](examples/README.md).

[INVENTORY.md](INVENTORY.md) records discovery and scope;
[python_surface.yaml](python_surface.yaml) records Python inclusion/exclusion;
[SCHEMA.md](SCHEMA.md) defines the data contract and maintenance workflow;
[COVERAGE.md](COVERAGE.md) reports verification and scientific limitations.

The `.yaml` files use JSON syntax, which is valid YAML 1.2. A standard-library JSON
reader is sufficient. The catalog is tied to an exact reviewed `dev` commit, not
to an assumed installed version. Distribute this directory with the matching
traceratops checkout or installation; the validation and fixture generator use
repository source/tests for evidence.

In a virtual environment with traceratops installed from this checkout:

```bash
python agent_tools/validate.py
python -m pytest agent_tools/test_validation.py -q
python agent_tools/examples/run_examples.py --work-dir /tmp/traceratops-examples
```

For CI, the validator and regression tests are lightweight; the full examples
also run plotting commands. Set `MPLBACKEND=Agg` in headless CI, and use writable
`MPLCONFIGDIR`/`XDG_CACHE_HOME` paths if the account has no writable home. Existing
CLI tests expect the environment's scripts directory on PATH.

Source drift is a review signal: update scientific descriptions and examples
before refreshing fingerprints. The validator never executes catalog-provided
Python code or shell strings. Running a trusted installed parser does import its
module, so validation should run against trusted traceratops installations.
