# H-CLI-13 raster reliability implementation plan

Status: corrected locally after remote failure; replacement CI verification pending

Date: 2026-08-11

## 1. Objective

Keep H-CLI-13 useful on GitHub-hosted runners without treating PNG glyph-edge
rasterization as a stable cross-environment artifact. The recipe must prove
that all six exports succeed. Its PNG contract is limited to successful
decode, expected dimensions, RGBA and alpha properties, and the absence of a
large visual change.

Exact SVG and Interactive SVG comparison remains the geometry and content
oracle. PDF, EPS, and PS retain their existing normalized body comparisons.

## 2. Remote failure that corrected the original plan

GitHub Actions run `31399138435`, job `93489423523`, used runner image
`ubuntu-24.04 20260720.247.2`. The bundled-font Fontconfig probe passed, but
H-CLI-13 still failed with:

- 3,541 changed RGB pixels;
- maximum RGB channel delta 129;
- bounding box `(76, 59, 975, 929)`;
- non-identical alpha values.

The run used CairoSVG 2.9.0, cairocffi 1.7.1, and Pillow 12.3.0. Its native
stack included Cairo 1.18.0, Fontconfig 2.15.0, FreeType 2.13.2, Pixman 0.42.2,
and libpng 1.6.43.

The exact failure was reproduced locally with the repository font config and
bundled fonts by using that complete native stack. The earlier 71-pixel result
changed Cairo alone while retaining newer Fontconfig, FreeType, Pixman, and
libpng versions. It was not a reproduction of the hosted runner.

The original decision to preserve the 100-pixel, delta-64, area-256, and exact
alpha limits is therefore superseded. Font isolation controls the selected
font files, but it does not pin the native rasterizer closure.

## 3. Contract decision

Keep the Linux bundled-font environment because it prevents substitution of a
different font file. Replace only H-CLI-13's PNG comparison with these checks:

1. The export command exits successfully and creates exactly the six declared
   files.
2. Pillow verifies and fully loads the PNG.
3. Pillow identifies it as PNG, its dimensions match the exact SVG master, and
   its mode is RGBA.
4. Its alpha channel contains both transparent and opaque pixels.
5. The generated and tracked PNGs are alpha-composited over black and white.
   For both composites, the largest per-channel RMS difference must be at most
   4.0.

The black and white composites expose both color and transparency changes. RMS
measures whole-image drift without requiring native libraries to place every
glyph-edge sample on the same pixel. Exact alpha equality, raw changed-pixel
counts, maximum single-pixel delta, and a full-resolution difference bounding
box are not part of this contract.

The shared `compare_raster_images` helper is unchanged because Gallery capture
tests still rely on its stricter alpha and bounding-box behavior.

## 4. Threshold evidence

The tracked PNG and the image reproduced with the job's native stack are both
1168 x 973 RGBA images. The runner-equivalent difference has maximum composite
RMS 2.240.

| Comparison against the tracked PNG | Maximum composite RMS |
| --- | ---: |
| Runner-equivalent native raster output | 2.240 |
| Recolor 6,411 gray pixels | 9.614 |
| Shift the figure by one pixel | 19.970 |
| Recolor the dominant blue | 34.735 |
| Blank image | 47.964 |

The 4.0 limit leaves measured headroom for the runner output and rejects every
tested material change. Exact SVG comparison independently catches changes to
layout, text, and colors before the PNG check is reached.

## 5. Files and ownership

- `docs/recipes/run_cli_scenarios.py`
  - owns H-CLI-13 PNG decode, structure, alpha, and visual checks;
  - retains the scenario-local Fontconfig environment.
- `tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py`
  - covers runner-like RGB and alpha noise;
  - rejects corruption, wrong format or dimensions, wrong mode, collapsed
    alpha range, blank output, recoloring, and geometry movement.
- `docs/internal/HCLI13_BUNDLED_FONT_RASTER_IMPLEMENTATION_PLAN_2026-08-10.md`
  - records the corrected failure boundary and verification evidence.

Production render code, public APIs, shared screenshot comparisons, published
H-CLI-13 artifacts, and `tests/reference_outputs/` remain unchanged.

## 6. Implementation phases

### Phase 1: replace the H-CLI-13 PNG oracle

Status: completed locally

- Remove the raw pixel-count, maximum-delta, bounding-box, and exact-alpha
  comparison for H-CLI-13 only.
- Verify and decode with Pillow.
- Compare black and white composites with the measured RMS limit.

### Phase 2: replace obsolete threshold tests

Status: completed locally

- Add a runner-like case with 3,541 RGB changes, delta 129, and changed alpha.
- Add structural and material-regression negative cases.
- Keep PDF, EPS, PS, SVG, Fontconfig, and clean-directory tests intact.

### Phase 3: verify the complete recipe partition

Status: completed locally

- Run the focused test module.
- Run H-CLI-13 from a clean checkout with `--check`.
- Run the standard non-heavy recipe partition.
- Confirm no public or reference artifact changed.

### Phase 4: verify on GitHub Actions

Status: pending

- Record a new GitHub Actions result separately. A local pass is not evidence
  of a remote pass.

## 7. Acceptance gates

```bash
/home/kawato/micromamba/bin/python -m pytest \
  tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py -q

PATH=/home/kawato/micromamba/bin:$PATH \
  /home/kawato/micromamba/bin/python \
  docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --check

/home/kawato/micromamba/bin/python -m pytest \
  tests/ -m "recipe and not recipe_heavy and not slow" -q

git diff --exit-code -- \
  docs/images/h-cli-13/ \
  tests/reference_outputs/
```

When the shared worktree contains unrelated fixture changes, run the
regeneration gates from `git archive HEAD` with only this plan's modified
source, test, and plan files overlaid. Do not refresh a tracked artifact to
make the check pass.

## 8. Evidence log

| Date | Evidence | Result |
| --- | --- | --- |
| 2026-08-10 | Original hosted-runner failure | PNG-only failure measured at 3,541 pixels and delta 129 |
| 2026-08-10 | Cairo-only experiment | 71-pixel result; later found not to reproduce the full runner stack |
| 2026-08-11 | Run 31399138435, job 93489423523 | Font isolation passed; the same PNG comparison failed again |
| 2026-08-11 | Full native-stack reproduction | Reproduced the job's RGB, alpha, delta, and bounding-box metrics exactly |
| 2026-08-11 | Composite RMS measurements | Runner 2.240; tested material changes 9.614 to 47.964 |
| 2026-08-11 | New helper against runner-equivalent artifact | Accepted |
| 2026-08-11 | Focused PNG contract tests | 2 passed |
| 2026-08-11 | Clean target module | 18 passed |
| 2026-08-11 | Clean H-CLI-13 six-format check | All six exports verified |
| 2026-08-11 | Clean standard recipe partition | 184 passed, 2,754 deselected |
| 2026-08-11 | Ruff on modified Python files | Passed |
| 2026-08-11 | Diff and whitespace checks | No semantic public or reference artifact change; passed |
| 2026-08-11 | Independent visual regression review | Passed with no blockers |

The shared workspace has pre-existing CRLF-only changes in documentation
artifacts, fixtures, and reference SVGs. Those changes are not owned by this
plan and must not be overwritten.

## 9. Completion rule

Set this document to `Status: completed` only after the local gates pass, the
diff review confirms that no generated artifact changed, and a replacement
GitHub Actions run passes. Until then, remote verification remains pending.
