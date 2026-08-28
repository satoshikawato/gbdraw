# AGENTS.md

Guidance for automated agents working in this repository.

## Read First
- コード1行増やすごとに技術負債が増えると心得よ。
- See `CLAUDE.md` for project-wide guidance.
- If working on the web UI, also read `gbdraw/web/CLAUDE.md`.

## Agent Branching

- Before starting work that may create a branch or commit, fetch `origin` and
  create a fresh work branch from the latest `origin/dev`:
  `git switch --no-track -c <branch-name> origin/dev`.
- All agent-authored commits must be made on a work branch derived from
  `origin/dev`. Do not base agent work on `main` or `master`.
- Do not configure a work branch to track `main` or `dev`. Before committing or
  pushing, verify the current branch and its upstream; publish only to the
  same-named remote work branch.
- Never commit directly to, or push directly to, `main` or `dev` unless the user
  explicitly authorizes that exact direct target in the current request.

## Project Summary

- `gbdraw` is a Python 3.10+ bioinformatics tool for publication-quality genome diagrams.
- Outputs: SVG, PNG, PDF, EPS, PS.
- Main package lives in `gbdraw/`; tests in `tests/`; web UI in `gbdraw/web/index.html` with JS modules under `gbdraw/web/js/` (app entry points in `gbdraw/web/js/app/` with submodules grouped by feature).

## Common Commands

```bash
# Fast tests (skip slow)
pytest tests/ -v -m "not slow"

# Full tests
pytest tests/ -v

# Compare generated SVGs with tracked references (read-only)
pytest tests/test_output_comparison.py::TestOutputComparison -v

# Update tracked references only for an intentional, reviewed geometry change
pytest tests/test_output_comparison.py::TestGenerateReferences --update-reference-outputs -v

# Lint
ruff check gbdraw/

# Prepare the generated browser wheel for offline web packaging/tests
python tools/prepare_browser_wheel.py

# Refresh the cache-bust token when preparing a deployable web bundle
python tools/prepare_browser_wheel.py --refresh-cache-bust

# Build
python -m build
```

- Allow at least 30 minutes for test commands before treating them as timed
  out, and monitor long runs incrementally. Keep shorter test-owned timeout
  assertions unchanged unless the task explicitly requires changing them.

## Browser / Playwright Checks

- Do not conclude that browser testing is unavailable just because `node_modules/`, `package.json`, or `@playwright/test` is missing at the repo root. This workspace may have Playwright installed through Python/conda instead.
- Check both paths when browser verification matters:
  - `command -v playwright && playwright --version`
  - `python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"`
- The JavaScript specs under `tests/web/*.playwright.spec.js` require Node's `@playwright/test`. Verify it with `node -e "console.log(require.resolve('@playwright/test'))"` before trying to run those specs.
- If Node's `@playwright/test` is unavailable, use Python Playwright for targeted browser checks instead of skipping browser verification.
- In Codex/agent sandboxes, Chromium may fail with `sandbox_host_linux.cc ... Operation not permitted`. When that happens, rerun the same local browser check with the required sandbox escalation rather than reporting that Playwright is unavailable.

## Expectations When Editing

- Audit production-code and test diffs separately before handoff. Optimize for
  fewer change points, branches, and duplicated behavior rather than raw line counts.
- Follow the [architecture fitness-function ratchet](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
  for every architecture-bearing change. Ordinary non-increasing changes use
  concise owner/path evidence and remove superseded paths in the same change.
  Complete before/after OE, PE, and CB sets are reserved for the exception
  conditions defined by that policy.
- Keep the web UI as a single-page app with no build step; `gbdraw/web/index.html` hosts HTML/CSS/templates and loads ES modules from `gbdraw/web/js/` (`app.js` entry with `app/`, `services/`, `utils/`).
- Keep larger UI modules split into focused subfolders under `gbdraw/web/js/app/` (for example `legend/`, `legend-layout/`, `feature-editor/`) and keep the `create*` entry points in the top-level `app/*.js` files.
- If adding CDN dependencies, update the CSP in `gbdraw/web/index.html`.
- Treat `tests/reference_outputs/` as read-only during normal tests. If diagram output changes intentionally, regenerate it with `--update-reference-outputs`, review the SVG diff, and rerun `TestOutputComparison`.
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
- Treat `gbdraw/web/gbdraw-<version>-py3-none-any.whl` as a generated, gitignored asset. Prepare it when tests or packaging need it, but do not commit it.

## Product Behavior Semantics Review

For changes affecting Web runtime or normative Web behavior contracts:

- Follow [`docs/internal/PRODUCT_IMPACT_RATCHET.md`](docs/internal/PRODUCT_IMPACT_RATCHET.md).
- Identify every reachable user-visible difference or continuation and separate
  the product outcome from the owner or path used to implement it.
- Search explicit supported behavior, compatibility commitments, authoritative
  domain rules, and accepted base-branch decisions. Current code and tests are
  evidence, not automatic product authority. Cite a `BD-###` only when it exists
  in the base branch.
- Before implementation, use the same policy's developer preflight when a
  proposed change may alter a material user effect but no registered concern is
  triggered. Use
  [`docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`](docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md)
  when evidence or Product judgment is required.
- Compare option realization with AND-of-OR requirements. Matching option IDs do
  not prove that independent behavior contributions remain satisfied.
- Keep the common path automatic for changes with no registered Product Impact
  delta or unresolved outcome.
- Stop the affected convergence when materially different user-visible outcomes
  remain unresolved. Do not select a product option autonomously.
- Present the Product Decision Owner with stable outcome choices, decision
  routes, and the documented `PRODUCT_DECISION` response template.
- After an explicit human choice, serialize only that outcome and show the
  generated machine representation for review. Do not infer missing rationale,
  retirement intent, or accepted residual risk, and do not broaden the choice.
- Candidate authority never authorizes the same candidate runtime.

## Showcase Figure Quality

- Treat `examples/gbdraw_social_preview.png` as an owner-maintained asset. Automated agents must not regenerate, replace, or edit it.
- Treat every public figure that demonstrates a feature as a finished example, even when the underlying change is small.
- Start from a realistic Gallery-quality recipe or session when one exists, and change only the setting being demonstrated where practical.
- Keep the labels, legend, color rules, record metadata, quantitative tracks, and comparison context that make the source figure useful. Remove an element only when the example has a clear reason to omit it.
- Keep minimal smoke diagrams in tests. Do not use them as tutorial images, Gallery entries, release examples, or other public showcases.
- Render and visually inspect the final artifact at a readable scale. Confirm that its documented command, session, or reproduction recipe generates the displayed figure.

## Pull request communication

- Before creating, editing, or reviewing a PR title or maintainer-facing body,
  read `.agents/skills/write-clear-pull-request/SKILL.md`.
- Write for a maintainer who has not read the implementation plan or earlier
  PRs. Internal labels are metadata, not explanations.
- Keep exact identifiers in backticks and explain their concrete role.
- Prepare the complete body file first, then run `gh pr create` or a
  wording-changing `gh pr edit` as a separate command so the hook can inspect it.
- Do not use `gh pr create --fill`, `--fill-first`, or `--fill-verbose`.

## Completion Handoff

- After completing an implementation, treat the session as one commit and provide a proposed commit title and summary in English.
