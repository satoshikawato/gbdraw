# AGENTS.md

Guidance for automated agents working in this repository.

## Read First

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

# Build
python -m build
```

## Expectations When Editing

- Keep the web UI as a single-page app with no build step; `gbdraw/web/index.html` hosts HTML/CSS/templates and loads ES modules from `gbdraw/web/js/` (`app.js` entry with `app/`, `services/`, `utils/`).
- Keep larger UI modules split into focused subfolders under `gbdraw/web/js/app/` (for example `legend/`, `legend-layout/`, `feature-editor/`) and keep the `create*` entry points in the top-level `app/*.js` files.
- If adding CDN dependencies, update the CSP in `gbdraw/web/index.html`.
- If diagram output changes, update reference SVGs in `tests/reference_outputs/`.
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
