# Contributing to gbdraw

Thanks for your interest in improving gbdraw. This guide covers the practical
steps for reporting a bug, proposing a change, and getting a pull request
merged.

## Before you start

- Search [existing issues](https://github.com/satoshikawato/gbdraw/issues)
  and [pull requests](https://github.com/satoshikawato/gbdraw/pulls) first;
  someone may already be working on the same thing.
- For a substantial change (a new CLI option, a new track type, a rendering
  behavior change), open an issue describing the problem before writing code.
  Small fixes and documentation corrections can go straight to a pull
  request.
- Read [`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md).

## Reporting a bug

Open a [GitHub issue](https://github.com/satoshikawato/gbdraw/issues/new) with:

- the `gbdraw` version (`gbdraw -h` prints it) and installation method
  (Bioconda, PyPI, or source),
- the exact command or Python snippet you ran,
- what you expected versus what happened, and
- a minimal input file, or a small excerpt from one, when the issue depends
  on specific input data.

## Development setup

```bash
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw
python -m pip install -e ".[dev]"
```

Add the optional export dependency if your change touches PNG/PDF/EPS/PS
output:

```bash
python -m pip install -e ".[dev,export]"
```

Non-SVG export also needs system Cairo/Pango libraries; see
[Installation](./docs/INSTALL.md#4-source-installation-for-development) for
platform-specific package names.

## Running the checks CI runs

```bash
# Lint (blocks CI on failure)
ruff check gbdraw/

# Fast test suite (excludes slow/regression-heavy tests)
pytest tests/ -v -m "not slow"

# Full suite, including slow tests (only run on push to main in CI)
pytest tests/ -v
```

Run a narrower slice while iterating:

```bash
pytest tests/test_regression.py -v
pytest tests/ -v -k "test_circular_basic"
pytest tests/ -v -m "circular"      # or -m "linear"
```

If your change affects the web app, see
[`gbdraw/web/CLAUDE.md`](./gbdraw/web/CLAUDE.md) for the JavaScript test
suite and Playwright checks.

## Changing rendered SVG output

Many tests compare generated SVGs against tracked references in
`tests/reference_outputs/`. If your change intentionally alters diagram
geometry:

```bash
# 1. Confirm the differences are the ones you intended
pytest tests/test_output_comparison.py::TestOutputComparison -v

# 2. Regenerate references only after reviewing that output
pytest tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs -v

# 3. Review the diff and re-verify
git diff -- tests/reference_outputs/
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Do not update reference outputs to make a failing test pass without first
confirming the new geometry is correct; the reference files are the
project's record of intended rendering behavior.

## Changing documentation

`docs/` follows a [Diátaxis](https://diataxis.fr/) structure (Tutorials,
How-to guides, Reference, Explanation); see
[`docs/DOCS.md`](./docs/DOCS.md) for the index and
[`docs/internal/DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md`](./docs/internal/DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md)
for the rationale behind it. Most commands and code blocks in Tutorials and
How-to guides are checked by the test suite against the current CLI, GUI,
and Python API:

```bash
pytest tests/ -v -k documentation
```

Keep a new page inside the existing structure (pick Tutorial, How-to,
Reference, or Explanation based on what the reader is trying to do) rather
than adding a new top-level document. Internal planning and audit documents
belong in `docs/internal/`, not alongside user-facing pages.

## Style

- Python: type hints throughout (`from __future__ import annotations`),
  `snake_case` for functions/modules, `PascalCase` for classes,
  `UPPER_SNAKE_CASE` for constants, and the `logging` module for diagnostics
  (no bare `print`).
- Prefer extending an existing configurator, drawer, or group module over
  adding a new parallel code path; see the Architecture section of
  [`CLAUDE.md`](./CLAUDE.md) for how the CLI, Web UI, and Python API are
  meant to converge on shared internals.
- Don't add speculative options, fallbacks, or abstractions for
  hypothetical future use; keep changes scoped to the problem at hand.

## Submitting a pull request

1. Create a branch from `main`.
2. Make your change with tests that cover it (a bug fix without a
   regression test is unlikely to be merged).
3. Run `ruff check gbdraw/` and `pytest tests/ -v -m "not slow"` locally.
4. Open a pull request against `main` describing what changed and why.
   Link the issue it resolves, if any.
5. CI runs the fast test suite on Python 3.10–3.12, a CairoSVG-specific job,
   and a lint job; all three must pass before merge.

By contributing, you agree that your contribution is licensed under the
project's [MIT License](./LICENSE.txt).
