# PR 319 CI, Cloudflare, and CodeQL remediation plan

## Objective

Restore every failing check observed for PR 319 and make the Cloudflare Pages build reproducible from repository-owned commands. Completion requires focused evidence for each root cause, regenerated documentation artifacts from their authoritative recipes, and a broader regression gate.

## Confirmed failure map

| Surface | Root cause | Required correction |
| --- | --- | --- |
| Cloudflare Pages | `tools/prepare_browser_wheel.py` invokes `python -m build`, but Cloudflare's automatic `pip install .` does not install the `dev` extra that owns `build`. | Build the browser wheel through pip's standard wheel command so build-system requirements come from `pyproject.toml` in an isolated backend environment. Do not add a build frontend to runtime dependencies. |
| GitHub Actions `test` matrix | Python Playwright is installed, but the Chromium executable is not. Three parametrized browser assertions fail before rendering. | Install Python Playwright's Chromium browser in every matrix job before tests. |
| GitHub Actions `test` matrix | The workflow overrides the repository timeout with 120 seconds, shorter than a recipe test's own 180-second subprocess allowance and the repository default of 300 seconds. | Remove the fast-test timeout override and use the `pyproject.toml` default. Keep test-owned timeouts unchanged. |
| GitHub Actions `losat-cache-browser-acceptance` | The runner imports `gbdraw.session_io`, but the dedicated job does not install the project. | Install the repository and core dependencies with `python -m pip install -e .` before preparing the browser wheel. |
| Documentation contracts | The interactive SVGs for H-CLI-13 and T-PY-08 embed an older copy of `standalone-interactivity-assets.js`. | Regenerate both scenarios through their authoritative recipe runners and commit only the changed generated SVGs. |
| CodeQL alert 41 | One-pass removal of multi-character HTML-like tags can expose a new tag after replacement. | Move comparison-record label formatting into the existing linear-comparison owner and strip innermost tags until the value reaches a fixed point. Cover nested-tag emergence. |
| CodeQL alert 42 | A test treats any occurrence of `raw.githubusercontent.com` as a trusted-host match. | Extract HTTP(S) URLs, parse them with `urllib.parse.urlsplit`, and compare the hostname exactly. Cover deceptive hosts and path substrings. |

The GitHub Actions, Cloudflare, and CodeQL failures are independent. A shared symptom is incomplete preparation or validation at a boundary; no single workaround can close all of them.

## Constraints

- Preserve unrelated worktree changes.
- Keep `build` out of runtime dependencies and do not restore a duplicate `requirements.txt` owner.
- Do not weaken or skip browser tests.
- Do not increase test-owned timeout assertions.
- Do not hand-edit generated SVGs or files under `tests/reference_outputs/`.
- Do not publish, push, deploy, or dismiss CodeQL alerts in this implementation session.

## Execution phases

### Phase 1 — Repository-owned build and CI preparation

Status: completed

Owned files:

- `tools/prepare_browser_wheel.py`
- `tests/test_web_packaging.py`
- `.github/workflows/test.yml`
- `recipe/build.sh`
- `recipe/meta.yaml`

Changes:

1. Replace the browser-wheel subprocess's dependency on `python -m build` with `python -m pip wheel --no-deps`, using an isolated build backend and a temporary wheel directory.
2. Lock the new command contract in the existing packaging unit test and remove the stale `build` availability skip from browser-wheel preparation.
3. Install Python Chromium in the regular test matrix.
4. Let fast tests inherit the repository's 300-second pytest timeout.
5. Install gbdraw in the dedicated LOSAT cache acceptance job.
6. Keep conda builds offline by disabling pip build isolation at that caller and remove its now-unused `python-build` host dependency.

Evidence:

- Focused packaging unit test passes.
- `python tools/prepare_browser_wheel.py` succeeds.
- `python tools/prepare_cloudflare_pages.py` succeeds in the ignored deployment output directory and produces the expected bundle.
- Workflow diff contains the three intended setup/timeout changes and no skipped coverage.

### Phase 2 — CodeQL findings

Status: completed

Owned files:

- `gbdraw/web/js/app/linear-comparisons.js`
- `gbdraw/web/js/app/app-setup.js`
- `tests/web/linear-comparisons.test.mjs`
- `tests/test_tutorial_documentation_contracts.py`

Changes:

1. Give `linear-comparisons.js` ownership of comparison-record display labels and use fixed-point tag removal.
2. Reuse that formatter from `app-setup.js` and add a Node regression for nested tags and empty fallback behavior.
3. Replace the raw-GitHub substring condition with exact parsed-host matching and add a Python regression for lookalike URLs.

Evidence:

- `node tests/web/linear-comparisons.test.mjs` passes.
- Focused tutorial documentation contract tests pass.
- Neither flagged source pattern remains.

### Phase 3 — Authoritative artifact regeneration

Status: completed

Owned files:

- `docs/images/h-cli-13/cli_export.interactive.svg`
- `docs/images/t-py-08/interactive_human_mitochondrion.interactive.svg`

Changes:

1. Run the H-CLI-13 runner without `--check`.
2. Run the T-PY-08 runner without `--check`.
3. Review the generated diffs and verify that only the embedded standalone interactivity payload changed.

Evidence:

- Both scenario runners pass with `--check` after regeneration.
- The all-Python recipe check passes, covering the previously hidden T-PY-08 drift.

### Phase 4 — Regression gates and audit

Status: completed

Verification order:

1. Packaging, CodeQL regression, and scenario-focused checks.
2. The three browser-rendered linear-definition cases with Chromium available.
3. The clean external Python recipe regeneration test under the repository timeout.
4. Web JavaScript tests and relevant packaging/documentation test modules.
5. Fast test suite, lint, and targeted diff checks when time and environment allow.

Completion evidence must record the exact command, result, and any environment limitation. A phase remains pending or in progress until its required behavior is observed.

## Evidence log

### Baseline diagnosis

- Phase: planning and root-cause isolation
- Behavior: mapped each red check to the first actionable failure in the supplied logs and reproduced both stale recipe artifacts locally with `--check`.
- Evidence: Cloudflare stops at `No module named build`; LOSAT stops at `No module named 'gbdraw'`; the test matrix reports one stale H-CLI-13 artifact, three missing-Chromium failures, and one 120-second pytest timeout. An independent all-Python recipe check also exposes stale T-PY-08 output.
- Deviations: none.
- Remaining risk: the source-controlled pip-wheel replacement and both CodeQL patterns still require implementation and focused verification.

### Phase 1 evidence

- Phase: repository-owned build and CI preparation
- Behavior: browser-wheel creation no longer imports the optional `build` frontend; Cloudflare uses pip's isolated PEP 517 backend, conda opts into its already-provisioned build environment, and both GitHub jobs install what their tests import or launch.
- Evidence: the isolated browser-wheel command completed and produced `gbdraw-0.14.0b0-py3-none-any.whl`; `python tools/prepare_cloudflare_pages.py` completed and created `dist/cloudflare-pages`; the no-isolation conda path also completed; the focused packaging selection passed 8 tests.
- Deviations: the shared conda caller required a new explicit `--no-build-isolation` option and removal of its unused `python-build` host dependency. This preserves offline conda builds while leaving Cloudflare isolated.
- Remaining risk: the updated workflow itself can only be exercised by the next GitHub Actions run; the local browser and LOSAT acceptance checks passed in Phase 4.

### Phase 2 evidence

- Phase: CodeQL findings
- Behavior: comparison labels remove nested complete tags to a fixed point and then remove unmatched angle brackets; documentation URL policy compares parsed hostnames rather than substrings.
- Evidence: `node tests/web/linear-comparisons.test.mjs` passed; both edited JavaScript modules passed `node --check`; all 13 tests in `tests/test_tutorial_documentation_contracts.py` passed; ruff passed for the Python contract file; the two originally flagged source patterns are absent.
- Deviations: a final single-character angle-bracket pass was added after fixed-point tag removal so an unmatched opener such as `<script` cannot survive.
- Remaining risk: GitHub must rerun CodeQL to close the remote alerts.

### Phase 3 evidence

- Phase: authoritative artifact regeneration
- Behavior: both stale interactive SVGs now embed the current standalone interactivity payload.
- Evidence: H-CLI-13 and T-PY-08 both passed their focused `--check`; the all-Python `--all --check` run verified every Python scenario, including T-PY-08 and all later scenarios.
- Deviations: H-CLI-13 also regenerated EPS, PDF, and PS creation timestamps. Their normalized recipe payloads passed, and those date-only diffs were removed so the final artifact change contains only the two intended interactive SVGs.
- Remaining risk: none for recipe freshness; the focused browser and timeout regressions passed in Phase 4.

### Phase 4 evidence

- Phase: regression gates and final audit
- Behavior: the formerly failing CI paths reach their intended assertions when their required runtime is present, and the full fast suite is green.
- Evidence: the clean external recipe test passed in 79.09 seconds under the repository's 300-second default; the three browser-rendered definition-gap cases passed in 6.71 seconds outside the Codex sandbox; the dedicated LOSAT cache browser runner passed both Playwright cases in 1.2 minutes; all 62 Node tests passed; the fast Python suite finished with 2,904 passed, 17 skipped, and 6 deselected in 15 minutes 23 seconds; production ruff, workflow YAML parsing, shell/Python syntax, reference-output immutability, and the targeted diff check passed.
- Deviations: Chromium cannot launch inside the Codex seccomp sandbox (`sandbox_host_linux.cc ... Operation not permitted`), so the same browser command was rerun outside the sandbox as repository guidance requires.
- Remaining risk: GitHub Actions, hosted CodeQL, and Cloudflare deployment must rerun after the branch is pushed. The full-worktree `git diff --check` still reports pre-existing trailing-whitespace/line-ending changes in unrelated user-owned files; the complete remediation file set passes its targeted check.

## Verification record

Repository-owned build paths:

```bash
# passed; isolated PEP 517 browser wheel
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python tools/prepare_browser_wheel.py

# passed; conda-style offline build environment
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python tools/prepare_browser_wheel.py --no-build-isolation

# passed; created dist/cloudflare-pages
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python tools/prepare_cloudflare_pages.py

# 8 passed, 30 deselected
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest tests/test_web_packaging.py -q \
  -k "prepare_browser_wheel or cloudflare or conda_build"
```

CodeQL and generated-artifact regressions:

```bash
# passed
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  node tests/web/linear-comparisons.test.mjs

# 13 passed
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest tests/test_tutorial_documentation_contracts.py -q

# both passed
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --check
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python docs/recipes/run_python_scenarios.py --scenario T-PY-08 --check

# passed all Python scenarios
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python docs/recipes/run_python_scenarios.py --all --check
```

Original CI failures and broad gates:

```bash
# 1 passed in 79.09 seconds
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest \
  tests/test_python_howto_recipe_contracts.py::test_python_evidence_recipes_regenerate_from_a_clean_external_context \
  -q

# 3 passed, 15 deselected; run outside the Codex sandbox
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest tests/test_linear_definition_alignment.py \
  -k browser_rendered_definition_gap -q

# 2 Playwright tests passed in 1.2 minutes; run outside the Codex sandbox
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python tests/run_losat_cache_browser_acceptance.py

# 62 passed
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  node --test tests/web/*.test.mjs

# 2,904 passed, 17 skipped, 6 deselected; run outside the Codex sandbox
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest tests/ -v --tb=short -m "not slow"

# passed
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

Syntax and scoped diff audit:

```bash
# all passed
bash -n recipe/build.sh
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m py_compile tools/prepare_browser_wheel.py
env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -c "import yaml; data=yaml.safe_load(open('.github/workflows/test.yml', encoding='utf-8')); assert isinstance(data, dict)"
git diff --check -- \
  .github/workflows/test.yml \
  tools/prepare_browser_wheel.py \
  recipe/build.sh recipe/meta.yaml \
  tests/test_web_packaging.py \
  gbdraw/web/js/app/linear-comparisons.js \
  gbdraw/web/js/app/app-setup.js \
  tests/web/linear-comparisons.test.mjs \
  tests/test_tutorial_documentation_contracts.py \
  docs/images/h-cli-13/cli_export.interactive.svg \
  docs/images/t-py-08/interactive_human_mitochondrion.interactive.svg
```
