# H-CLI-13 bundled-font raster implementation plan

Status: implemented and locally verified; remote CI verification pending

Date: 2026-08-10

## 1. Objective

Make H-CLI-13 reproducible on GitHub-hosted Ubuntu runners without weakening
its existing artifact comparisons. The documented command must still generate
SVG, Interactive SVG, PNG, PDF, EPS, and PS in one real invocation, and stale
or materially changed public artifacts must still fail.

This plan corrects the environment remediation recorded in
`CI_TEST_RELIABILITY_AND_RUNTIME_IMPLEMENTATION_PLAN_2026-08-10.md`. It does
not change that plan's test partitioning or runtime work.

## 2. Failure evidence

The 2026-08-10 GitHub Actions log in `job-logs3.txt` has one primary failure:
H-CLI-13 passes exact static SVG and Interactive SVG comparison, then its PNG
comparison reports 3,541 changed RGB pixels, maximum channel delta 129,
bounding box `(76, 59, 975, 929)`, and non-identical alpha. The run ends with
178 passed and 1 failed. It is not a timeout or package-install failure.

Both attempted environment fixes were no-ops:

- `ubuntu-latest` and `ubuntu-24.04` selected the same Ubuntu 24.04 runner
  family at the time of the run. A hosted-runner label does not pin the weekly
  image or its native raster stack.
- `apt-get install fonts-liberation` reported that version `1:2.1.5-3` was
  already installed and installed zero packages.

The repository currently has two font owners:

- text measurement explicitly loads the Liberation Sans files bundled under
  `gbdraw/data`;
- CairoSVG receives SVG family names and resolves the raster font through the
  host Fontconfig configuration.

The exact SVG pass localizes the failure after layout and SVG assembly. The
wide PNG and alpha difference is consistent with an uncontrolled raster font
selection/configuration boundary, not with the already measured Cairo-only
noise.

## 3. Controlled reproduction

A Fontconfig configuration was tested with this order:

1. load the platform's default `conf.d` rendering and alias rules;
2. reset every font directory contributed by those rules;
3. add only the repository's `gbdraw/data` font directory;
4. place Fontconfig caches under an isolated XDG cache directory.

The tested configuration resolved all four requested Liberation Sans faces to
the tracked TTF files:

- `LiberationSans-Regular.ttf`;
- `LiberationSans-Bold.ttf`;
- `LiberationSans-Italic.ttf`;
- `LiberationSans-BoldItalic.ttf`.

Rendering the tracked H-CLI-13 SVG in separate processes produced:

| Renderer | Difference from tracked PNG |
| --- | --- |
| Cairo 1.18.4 with isolated bundled fonts | exact byte match |
| Cairo 1.18.0 with isolated bundled fonts | 71 RGB pixels, maximum delta 62, 13 x 11 bounding box, exact alpha |

The existing PNG allowance is 100 pixels, channel delta 64, bounding-box area
256, and exact alpha. The controlled Cairo 1.18.0 result fits that allowance
without changing any threshold. The Cairo 1.18.4 SHA-256 is the tracked PNG
SHA-256, `6e41714a670f1ddc31c2739fd66b5fcd98d6d3a0952a14e274916c9317af9d2d`.

A minimal bundled-only configuration was rejected during investigation. It
omitted platform rendering rules and changed 7,916 alpha pixels. The default
rules must be loaded before resetting font directories.

## 4. Decision

Keep every current H-CLI-13 comparison and make font resolution owned by the
recipe runner.

- Add `docs/recipes/h-cli-13-fonts.conf`.
- On Linux, the configuration uses the literal
  `<include ignore_missing="no" prefix="default">conf.d</include>`, then resets
  inherited font directories, adds only `gbdraw/data`, and uses an XDG cache
  directory. Other platforms retain their current environment.
- For H-CLI-13 on Linux only, the runner overwrites `FONTCONFIG_FILE`, removes
  inherited `FONTCONFIG_PATH`, `FONTCONFIG_SYSROOT`, `HOME`, `FC_DEBUG`, and
  `FC_DBG_MATCH_FILTER`, and points `XDG_CACHE_HOME` and `XDG_CONFIG_HOME` at a
  separate temporary directory. Removing `HOME` prevents the default
  `50-user.conf` rule from loading deprecated `~/.fonts.conf*` files without
  repurposing the process home directory.
- On Linux, fail before rendering unless `fc-match` is available, resolves
  Regular, Bold, Italic, and Bold Italic to the corresponding repository TTF,
  and reports the expected effective antialias, hinting, autohint, hint-style,
  subpixel, and LCD-filter properties. Accept the Fontconfig RGB value used by
  headless Ubuntu and the no-subpixel value reported by Wayland-backed
  image-surface sessions; both reproduced the tracked PNG. Include unexpected
  paths and property values in the error.
- Preserve exact SVG and Interactive SVG comparison, bounded PNG comparison,
  and normalized PDF/EPS/PS body comparison.
- Remove both no-op `fonts-liberation` installation steps from the workflow.
- Do not regenerate a tracked artifact when the controlled configuration
  reproduces it within the existing contract. An unexpected artifact diff is
  a stop condition, not permission to refresh the golden.

The repository configuration fixes the font-file owner. The existing narrow
PNG allowance continues to cover the separately measured Cairo 1.18.0 versus
1.18.4 raster noise. This plan does not claim to make every native library
bit-identical.

## 5. Rejected alternatives

### Increase the PNG tolerance

Rejected. The failure exceeds every threshold independently. Widening the
threshold could hide real color, geometry, or alpha regressions.

### Install a system font

Rejected. It was already installed, and a package name does not prove which
font file Fontconfig selected. It also leaves measurement and conversion with
different font owners.

### Replace binary comparisons with signatures or provenance only

Rejected after independent plan review. A provenance hash can detect
post-publication mutation but cannot detect that a checked-in derivative became
stale after converter-only behavior changed. The current PNG and normalized
PDF/EPS/PS comparisons retain that coverage.

### Add a canonical renderer container

Deferred. A digest-pinned container is the correct escalation if isolated font
selection still cannot keep the observed native differences inside the
existing bounds. The repository has no artifact-generation container owner,
so it is not the smallest justified change while the tested font boundary is
sufficient.

## 6. Files and ownership

- `docs/recipes/h-cli-13-fonts.conf`
  - isolated font directories and cache declaration.
- `docs/recipes/run_cli_scenarios.py`
  - H-CLI-13 environment construction and fail-fast font resolution.
- `tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py`
  - configuration, hostile environment, and resolution regression tests.
- `.github/workflows/test.yml`
  - removal of two no-op apt steps.
- `docs/internal/CI_TEST_RELIABILITY_AND_RUNTIME_IMPLEMENTATION_PLAN_2026-08-10.md`
  - evidence-backed correction to its earlier remote diagnosis.

Production rendering modules, public APIs, comparison thresholds, published
H-CLI-13 artifacts, and `tests/reference_outputs/` are out of scope unless a
verification gate proves the controlled reproduction evidence wrong.

## 7. Implementation phases

### Phase 1: add the isolated Fontconfig owner

Status: completed

- Add the tested include-reset-dir-bundled-dir configuration.
- Give the cache a deterministic subdirectory name under `XDG_CACHE_HOME`.
- Keep the file limited to font discovery; do not duplicate Fontconfig's
  platform hinting and antialias rules.

### Phase 2: apply it only to H-CLI-13

Status: completed

- Create a separate temporary XDG directory outside the scenario workdir so
  `assert_exact_workdir_files()` cannot see cache files.
- On Linux, assign, rather than default, `FONTCONFIG_FILE`, `XDG_CACHE_HOME`,
  and `XDG_CONFIG_HOME`; remove inherited search-path, sysroot, home, and debug
  overrides. Leave non-Linux environments unchanged.
- Resolve and verify all four Liberation Sans styles before starting gbdraw on
  Linux, including the effective raster properties required by the tracked
  golden.
- Leave every other scenario environment unchanged.

### Phase 3: remove the false CI remediation

Status: completed

- Delete both apt update/install font steps.
- Keep the existing Ubuntu runner selection and dependency/test commands.
- Do not add a second font-install or environment path.

### Phase 4: test and record evidence

Status: completed locally; remote GitHub Actions run pending

- Prove the checked-in configuration exposes only bundled font files.
- Start with hostile inherited `FONTCONFIG_FILE`, `FONTCONFIG_PATH`,
  `FONTCONFIG_SYSROOT`, `XDG_CACHE_HOME`, and `XDG_CONFIG_HOME` values and prove
  the runner replaces or removes them as specified.
- Keep the existing PNG threshold and PDF/EPS/PS negative tests unchanged.
- Run H-CLI-13 with `--check` and confirm all six formats pass.
- Run the standard recipe partition and confirm reference outputs and public
  H-CLI-13 artifacts are unchanged.

## 8. Required tests

- The config parses and contains the literal compiled-default `conf.d` include
  before `reset-dirs`, followed by the relative bundled font directory.
- On Linux, `fc-list` under the isolated environment returns paths only below
  `gbdraw/data`.
- Regular, Bold, Italic, and Bold Italic resolve to the expected filenames.
- Effective values are antialias true, hinting true, autohint false,
  `hintslight`, RGB or no subpixel order, and the default LCD filter; this
  proves the included rules were loaded rather than merely checking XML order.
- A hostile inherited Fontconfig environment, including a process home and
  Fontconfig debug variables, is overwritten or removed on Linux.
- A non-Linux environment is returned unchanged.
- A wrong resolved path raises `RecipeContractError` before gbdraw executes.
- The existing accepted 71-pixel Cairo case still passes.
- The existing 101-pixel, channel-delta-65, scattered-bounding-box, alpha,
  mode, and dimension cases still fail.
- PDF/EPS/PS renderer metadata normalization still rejects changed content.
- The clean-directory H-CLI-13 check still creates and verifies exactly the six
  declared formats.

## 9. Acceptance gates

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

Completion requires the focused tests, clean-directory H-CLI-13 check, and
available broader local gates to pass. Remote GitHub-hosted verification must
be recorded separately after a workflow run; it must not be inferred from
local results.

## 10. Stop conditions

- The isolated configuration resolves any requested Liberation Sans face
  outside `gbdraw/data`.
- Cairo 1.18.0 versus 1.18.4 exceeds the existing PNG bounds after font
  isolation.
- SVG or Interactive SVG changes.
- A public H-CLI-13 artifact must be regenerated merely to make the test pass.
- The fix requires changing a production render API or global process
  environment.

If a stop condition is reached, do not widen the comparison. Amend this plan
to introduce a digest-pinned artifact-generation container with exact native
library and font hashes.

## 11. Evidence log

| Date | Phase | Evidence | Result |
| --- | --- | --- | --- |
| 2026-08-10 | Diagnosis | `job-logs3.txt`, first failure and summary | Confirmed one deterministic PNG-only failure |
| 2026-08-10 | Diagnosis | workflow apt output | Confirmed `fonts-liberation` installed zero packages |
| 2026-08-10 | Controlled reproduction | bundled-font config, Cairo 1.18.4 | Exact tracked PNG SHA-256 reproduced |
| 2026-08-10 | Controlled reproduction | same config, Cairo 1.18.0 | 71 pixels / 62 max / 13 x 11 / exact alpha |
| 2026-08-10 | Plan review | independent read-only critic | Rejected provenance-only oracle; retained native body comparisons |
| 2026-08-10 | Focused tests | target module excluding clean-fixture regeneration | 11 passed |
| 2026-08-10 | Clean checkout test | target module with current runner, tests, and config | 17 passed in 21.65 seconds |
| 2026-08-10 | Clean H-CLI-13 check | real six-format generation and published comparison | All six formats verified |
| 2026-08-10 | Standard recipe gate | `recipe and not recipe_heavy and not slow` | 183 passed in 84.65 seconds |
| 2026-08-10 | Static validation | Ruff, workflow YAML, Fontconfig XML, diff check | Passed |
| 2026-08-10 | Implementation review | independent read-only scoped diff review | No correctness blockers |

The shared workspace contains pre-existing CRLF-only changes in documentation
artifacts, fixtures, and reference SVGs. Direct regeneration tests therefore
report unrelated fixture-manifest failures in that dirty tree. The clean
checkout gates above used `git archive HEAD` plus only the files owned by this
plan. They did not update public H-CLI-13 artifacts. The tracked PNG and PDF
remain byte-identical to `HEAD`; semantic diffs under the public H-CLI-13 and
reference-output paths are empty when end-of-line-only changes are ignored.

## 12. Completion rule

Set this document to `Status: completed` only after the applicable local gates
pass, every phase has recorded evidence, and production, test, workflow,
documentation, and generated-artifact diffs have been reviewed separately.
