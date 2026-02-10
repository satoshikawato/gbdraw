# Preprint Readiness Checklist

This document tracks preprint-oriented readiness items for the current
codebase and documentation.

Baseline version:
- Display version: 0.9.0-beta
- Package version (PEP 440): 0.9.0b0

## 1) Stability and API polish

Status: mostly complete.

- `sys.exit()` is confined to CLI entry points (`gbdraw/cli.py`); library
  modules raise exceptions.
- `gbdraw/config/toml.py` raises `ConfigError` with actionable messages when
  config loading fails.

## 2) Documentation consistency (Quickstart and Tutorials)

Status: covered.

- Linear selectors are documented across:
  - `docs/QUICKSTART.md`
  - `docs/TUTORIALS/1_Customizing_Plots.md`
  - `docs/TUTORIALS/2_Comparative_Genomics.md`
  - `docs/TUTORIALS/3_Advanced_Customization.md`
  - `docs/CLI_Reference.md`
- Covered options:
  - `--record_id`
  - `--reverse_complement`
  - `--region`

## 3) Roadmap alignment

Status: tracked in `docs/ROADMAP.md`.

- Short-term items include tick label configuration and feature-wise SVG
  grouping.
- Long-term items include label-overlap handling, custom tracks, and PAF/MAF
  support.
