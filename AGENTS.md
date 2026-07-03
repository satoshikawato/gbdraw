# AGENTS.md

Guidance for automated agents working in this repository.

## Read First
- コード1行増やすごとに技術負債が増えると心得よ。
- See `CLAUDE.md` for project-wide guidance.
- If working on the web UI, also read `gbdraw/web/CLAUDE.md`.

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

# Lint
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503

# Prepare the generated browser wheel for offline web packaging/tests
python tools/prepare_browser_wheel.py

# Refresh the cache-bust token when preparing a deployable web bundle
python tools/prepare_browser_wheel.py --refresh-cache-bust

# Build
python -m build
```

## Browser / Playwright Checks

- Do not conclude that browser testing is unavailable just because `node_modules/`, `package.json`, or `@playwright/test` is missing at the repo root. This workspace may have Playwright installed through Python/conda instead.
- Check both paths when browser verification matters:
  - `command -v playwright && playwright --version`
  - `python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"`
- The JavaScript specs under `tests/web/*.playwright.spec.js` require Node's `@playwright/test`. Verify it with `node -e "console.log(require.resolve('@playwright/test'))"` before trying to run those specs.
- If Node's `@playwright/test` is unavailable, use Python Playwright for targeted browser checks instead of skipping browser verification.
- In Codex/agent sandboxes, Chromium may fail with `sandbox_host_linux.cc ... Operation not permitted`. When that happens, rerun the same local browser check with the required sandbox escalation rather than reporting that Playwright is unavailable.

## Expectations When Editing

- Keep the web UI as a single-page app with no build step; `gbdraw/web/index.html` hosts HTML/CSS/templates and loads ES modules from `gbdraw/web/js/` (`app.js` entry with `app/`, `services/`, `utils/`).
- Keep larger UI modules split into focused subfolders under `gbdraw/web/js/app/` (for example `legend/`, `legend-layout/`, `feature-editor/`) and keep the `create*` entry points in the top-level `app/*.js` files.
- If adding CDN dependencies, update the CSP in `gbdraw/web/index.html`.
- If diagram output changes, update reference SVGs in `tests/reference_outputs/`.
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
- Treat `gbdraw/web/gbdraw-<version>-py3-none-any.whl` as a generated, gitignored asset. Prepare it when tests or packaging need it, but do not commit it.
