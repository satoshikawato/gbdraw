# AGENTS.md

Guidance for automated agents working in this repository.

## Read First

- See `CLAUDE.md` for project-wide guidance.
- If working on the web UI, also read `gbdraw/web/CLAUDE.md`.

## Project Summary

- `gbdraw` is a Python 3.10+ bioinformatics tool for publication-quality genome diagrams.
- Outputs: SVG, PNG, PDF, EPS, PS.
- Main package lives in `gbdraw/`; tests in `tests/`; web UI in `gbdraw/web/index.html`.

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

- Keep `gbdraw/web/index.html` as a single-file SPA; no new build step.
- If adding CDN dependencies, update the CSP in `gbdraw/web/index.html`.
- If diagram output changes, update reference SVGs in `tests/reference_outputs/`.
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
