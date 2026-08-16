# AGENTS.md

Guidance for automated agents working in this repository.

## Read First
- コード1行増やすごとに技術負債が増えると心得よ。
- See `CLAUDE.md` for project-wide guidance.
- If working on the web UI, also read `gbdraw/web/CLAUDE.md`.
- Before moving, centralizing, splitting, or replacing existing Web behavior,
  read and apply
  `.agents/skills/refactor-gbdraw-web-safely/SKILL.md`.
- Behavior-preserving Web refactors establish base characterization tests
  before production edits. The new implementation must not serve as its own
  behavioral oracle.

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

# Web behavior-preserving refactor guardrails
npm run test:web:refactor-guards

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

## Behavior-preserving Web refactors

- Characterize the base before editing production code.
- Keep renderer-owned state and editor-owned overrides as separate provenance
  domains.
- Treat validated Feature catalogs, extracted Features, biological Features,
  orthogroups, sequences, SVG, and Results as large borrowed owners. Do not
  JSON-clone or JSON-serialize them as defensive copying.
- Keep mounted SVG and selected Result ownership synchronous. If ownership must
  cross `await`, use and revalidate an explicit lease.
- Use production-connected counters and independent oracles. Dummy counters and
  comparison of two paths using the same new helper are not evidence.
- Run `npm run test:web:refactor-guards` before handoff.
- Require a separate adversarial review for cross-cutting Web refactors.
- Do not declare completion while a required CI failure is red or unclassified.

## Expectations When Editing

- Audit production-code and test diffs separately before handoff. Optimize for
  fewer change points, branches, and duplicated behavior rather than raw line counts.
- Keep the web UI as a single-page app with no build step; `gbdraw/web/index.html` hosts HTML/CSS/templates and loads ES modules from `gbdraw/web/js/` (`app.js` entry with `app/`, `services/`, `utils/`).
- Keep larger UI modules split into focused subfolders under `gbdraw/web/js/app/` (for example `legend/`, `legend-layout/`, `feature-editor/`) and keep the `create*` entry points in the top-level `app/*.js` files.
- If adding CDN dependencies, update the CSP in `gbdraw/web/index.html`.
- Treat `tests/reference_outputs/` as read-only during normal tests. If diagram output changes intentionally, regenerate it with `--update-reference-outputs`, review the SVG diff, and rerun `TestOutputComparison`.
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
- Treat `gbdraw/web/gbdraw-<version>-py3-none-any.whl` as a generated, gitignored asset. Prepare it when tests or packaging need it, but do not commit it.

## Showcase Figure Quality

- Treat `examples/gbdraw_social_preview.png` as an owner-maintained asset. Automated agents must not regenerate, replace, or edit it.
- Treat every public figure that demonstrates a feature as a finished example, even when the underlying change is small.
- Start from a realistic Gallery-quality recipe or session when one exists, and change only the setting being demonstrated where practical.
- Keep the labels, legend, color rules, record metadata, quantitative tracks, and comparison context that make the source figure useful. Remove an element only when the example has a clear reason to omit it.
- Keep minimal smoke diagrams in tests. Do not use them as tutorial images, Gallery entries, release examples, or other public showcases.
- Render and visually inspect the final artifact at a readable scale. Confirm that its documented command, session, or reproduction recipe generates the displayed figure.

## Completion Handoff

- After completing an implementation, treat the session as one commit and provide a proposed commit title and summary in English.
