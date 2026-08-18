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

Public documentation has four routes: [Tutorials](./docs/TUTORIALS/README.md),
[Technical documentation](./docs/REFERENCE/README.md),
[FAQ](./docs/FAQ.md), and [Gallery](./docs/GALLERY.md). See
[`docs/DOCS.md`](./docs/DOCS.md) for the index and the
[current implementation plan](./docs/internal/DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md)
for the ownership and evidence rules.

Use the fewest pages that answer distinct reader questions. Before adding a
public page, identify the reader question, the existing owner, and whether the
change should keep, merge, delete, or create a page. A separate GUI, CLI, or
Python test scenario does not require a separate public page. Technical
documentation owns exact contracts; FAQ answers choices and troubleshooting by
linking to those contracts.

Tutorial commands, public code blocks, and internal evidence scenarios are
checked against the current CLI, GUI, and Python API:

```bash
pytest tests/ -v -k documentation
```

Extend an existing owner before adding a new top-level document. Internal
plans, audits, and executable evidence registries belong in `docs/internal/`,
not alongside user-facing pages.

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

## Branch roles

- Work-branch base: the latest merged `origin/dev`.
- Pull-request target: `dev`.
- Integration branch: `dev`.
- Release branch: `main`.

Agent-authored work branches must start from `origin/dev` without tracking
`dev` or `main`. Maintainers promote `dev` to `main` through a pull request
after the integration checks pass. Feature work must not target `main`
directly.

## Submitting a pull request

1. Fetch `origin` and create a non-tracking work branch from `origin/dev`:

   ```bash
   git fetch origin
   git switch --no-track -c <branch-name> origin/dev
   ```

2. Make your change with tests that cover it (a bug fix without a
   regression test is unlikely to be merged).
3. Run `ruff check gbdraw/` and `pytest tests/ -v -m "not slow"` locally.
4. Open a pull request against `dev` describing what changed and why.
   Link the issue it resolves, if any.
5. Complete exactly one declaration under `Architecture impact` in the pull
   request template. Architecture-bearing changes must include the concrete
   before and after sets required by the
   [architecture fitness-function ratchet](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md).
6. For an architecture-bearing change, before merge, record the permalink to
   the architecture owner's structured decision for the exact final head SHA.
   The author declaration is evidence, not approval; a later commit requires a
   new owner decision.
7. Wait for all required CI checks to pass before merge.

By contributing, you agree that your contribution is licensed under the
project's [MIT License](./LICENSE.txt).
