# Web audit sweeps

These are the manual sweeps from the 2026-09-30 Web GUI audit
(`docs/internal/web-gui-audit-20260930/`). They are too slow or too broad
for CI, so they live here as tools. CI runs the guards that were derived from
them (G-A to G-J); these sweeps cover the long tail. Run them during the periodic
audit before a `dev` to `main` promotion (see
[`docs/internal/WEB_PERIODIC_AUDIT.md`](../../docs/internal/WEB_PERIODIC_AUDIT.md)),
or when you investigate a finding.

A sweep records evidence; a failing sweep is not a release block by itself.
Each confirmed problem becomes an audit row with a severity, and then a guard
or a recorded decision not to add one.

## Setup

```bash
python tools/prepare_browser_wheel.py          # the app loads the local wheel
export GBDRAW_WEB_TEST_PORT=<free port>         # one port per concurrent run
export GBDRAW_AUDIT_OUT=/path/to/evidence       # default: <os tmpdir>/gbdraw-audit
```

Node `@playwright/test` must resolve from the repository root (`node_modules`).
Every sweep runs with one worker under its own config:

```bash
npx playwright test -c tools/audit/playwright.audit.config.cjs <file-or-filter>
```

The config serves the repository root with `python3 -m http.server`, like
`playwright.config.js`. Evidence never goes into the repository.

## Sweeps

| File | What it checks | Options | Evidence (`$GBDRAW_AUDIT_OUT/...`) |
| --- | --- | --- | --- |
| `parity-sweep.audit.spec.js` | GUI Generate against the Source recipe: a baseline and 31 Circular (HmmtDNA) and 25 Linear (MjeNMV + MelaMJNV + BLAST table) option probes. Each probe starts from a freshly opened app and saves the Result SVG, the recipe command, and the Run info helper files. | `AUDIT_ONLY=probe1,probe2`, `AUDIT_INPUT=<file in tests/test_inputs>` | `parity/circular-<input>/`, `parity/linear-MJNV/` |
| `parity_replay.py` | Replays each saved recipe with this checkout's CLI and compares the CLI SVG with the GUI SVG through `tests/utils/svg_compare.compare_svgs`. It ignores the same binding attributes as the Gallery publication parity check, plus the root `baseProfile`. Probes without a recipe are listed as `NO_RECIPE` with the GUI's reason. | probe directory and input files | `<probe dir>/replay-report.json` |
| `xss-sweep.audit.spec.js` | A synthetic record with script-like DEFINITION, ORGANISM, qualifiers, and file name. Checks the preview, search, the feature popup, Run info, static and interactive SVG export, and the exported interactive SVG opened on its own. Any dialog fails the test. | `XSS_FILE=<GenBank>`, `XSS_FNAME=<upload name>` | `xss/` |
| `viewport-sweep.audit.spec.js` | Horizontal overflow of the app (open, after upload, after Generate) at 390, 768, 1280, and 1920 px, and of the Gallery page and its tabs at 390, 768, and 1280 px. Records the elements past the viewport and saves screenshots. | `AUDIT_WIDTHS=390,768` | `viewport/` |
| `gallery-roundtrip.audit.spec.js` | Every Gallery session: Load, Save, Load the saved file in a fresh context, and Save again. User-owned state (the G-C snapshot from `tests/web/helpers/app-lifecycle.cjs`) and the saved documents must match, apart from save timestamps. | `AUDIT_SESSIONS=name1,name2` | `gallery-roundtrip/` |
| `session-load-timing.audit.spec.js` | Session load time and the order of Worker starts, Worker messages, and session lifecycle events after the load starts. | `AUDIT_SESSIONS=name1,name2` | `session-load-timing/` |

Example: run the full parity sweep, then replay it.

```bash
npx playwright test -c tools/audit/playwright.audit.config.cjs parity-sweep
python tools/audit/parity_replay.py "$GBDRAW_AUDIT_OUT/parity/circular-HmmtDNA_gbk" tests/test_inputs/HmmtDNA.gbk
python tools/audit/parity_replay.py "$GBDRAW_AUDIT_OUT/parity/linear-MJNV" \
  examples/MjeNMV.gb examples/MelaMJNV.gb examples/MjeNMV.MelaMJNV.tblastx.out
```

The parity sweep reopens the app for every probe, so the full run takes longer
than the others. The 2026-09-30 audit estimated about 20 minutes for both modes.
Use `AUDIT_ONLY` to repeat single probes.

## Investigation helpers

`helpers/cdp-debug.cjs` has two Chrome DevTools Protocol helpers for a throwaway
spec:

- `captureExceptions(page)` records every exception thrown while armed,
  including caught ones that never reach `page.on('pageerror')`.
- `captureAtBreakpoint(page, urlRegex, line, expression)` evaluates an
  expression each time a source line runs.

`helpers/audit-common.cjs` holds the shared helpers: Generate with the Source
recipe, helper-file download, Session load and save, and overflow measurement.
It reuses `tests/web/helpers/app-lifecycle.cjs` and does not copy it.

## Keeping the sweeps current

The sweeps drive the app through `window.__GBDRAW_APP__` and visible labels. If
a UI change breaks a sweep, fix the sweep in the same PR as the UI change or in
the next periodic audit. Do not move a sweep into CI. Instead, turn the bug
class it found into a small guard under `tests/` (see W8 in the audit folder).
