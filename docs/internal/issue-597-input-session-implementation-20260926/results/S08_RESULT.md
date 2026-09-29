# Issue #597 / S08 integrated acceptance — 2026-09-29

Status: **implementation and integration committed for review; integrated acceptance remains incomplete**. The two expressly permitted observations remain classified as measured FAIL where they occur: real large Session heartbeat maximum above 500 ms and one Python diagram Worker on Load of a real Linear saved preview with non-default saved configuration. The original gzip is lost and native structured-clone wire/copy bytes remain UNAVAILABLE. Other failed gates below have no waiver. BUG-01 was not started.

## Checkout, authority, and integration

Work used only `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovered-20260928` on `fix/issue-597-input-session-20260926`; raw evidence is in the sibling `issue597-S05-recovery-evidence-20260928/S08/`. All browser, pytest, and temporary files were under `/home/kawato/gbdraw-issue597-s08-scratch`, never `/tmp`. The shared checkout and Issue #619 work were preserved. Intake classification and exact remote/PR state are in `S08/intake/`. Actual `origin/dev` was `57cef3ba47f4b7790a9145f2ce6988422a00710e`; S08 integrated it in separate merge `c127d361e1372b5c8771f8a22d8bc81e0da8bd93` before editing. PD-OI-044/045/046, Product Contract revision 26, and trusted `dev` are the behavior and privacy authority; no decision or guard file was edited.

Six conflict resolutions retain both #597 and trusted `dev`:

| File | Resolution |
| --- | --- |
| `gbdraw/web/index.html` | Kept #597 discovery, availability, and disclosure bindings with `dev`'s preview layout hint and orthogroup panel sizing. |
| `gbdraw/web/js/app/run-analysis.js` | Kept #597 reflow/operation guards and `dev` decoration-continuity capture. |
| `gbdraw/web/js/services/config.js` | Kept #597 import candidate/Worker ownership and `dev` saved active-mode restoration. |
| `gbdraw/web/js/services/error-normalization.js` | Kept #597 error vocabulary and `dev` `DECORATION_CONTINUITY`. |
| `tests/web/feature-selection.test.mjs` | Retained both layout-target and pairwise SVG expectations. |
| `tests/web/right-drawer.playwright.spec.js` | Retained the #597 rollback expectation and `dev` right-drawer assertions. |

S08 removed the duplicate floating zoom toolbar created by that merge, retaining `dev`'s docked controls and #597's disabled Reset Layout during a Session operation. `Source records` now renders the normalized error summary without automatically exposing the filename (PD-OI-046); the exact NO_RECORDS summary and absent filename have a browser assertion. `exportSession` returns the normalized error with its failure status, matching the existing Save lifecycle contract. Known Session size, unsupported gzip, and unavailable import Worker errors now retain actionable structured codes instead of `UNKNOWN`. When a protein run computes reusable evidence but a same-row layout has no displayed comparisons, its empty canonical comparison request is committed without projecting a nonexistent typed comparison resource. Nonempty comparison requests still use the existing projection validation.

The F3 status after an operation error with no Result is reproducible on trusted `dev`: `Invalid settings · Canonical resource record-1-genbank is missing`. `compareGenerationIntent` checks the resource before the no-Result state, and `projectGenerationIntent` exposes that exception. This is a separate status semantics issue candidate; changing `Invalid settings` to `Not generated` would change the user's stated cause, so S08 does not choose that Product outcome. Repro: open without a Result, trigger an operation error, then inspect the generation status (`S07/observe-errors-current/observation.json`).

## Acceptance matrix

The authoritative definitions are in `MASTER_PLAN.md` § Acceptance (D-01–D-04, S-01–S-05, A-01, W-01). `S08/` paths below are raw evidence. Earlier S03–S07 results are historical evidence; reuse as current PASS is limited to unchanged source, input, environment, and assertion conditions.

| ID | Current S08 evidence and scope | Classification |
| --- | --- | --- |
| D-01 | `S08/focused-browser-final.log`: native GenBank/DDBJ and GFF+FASTA, packaged Worker path, and helper settlement among 13 passing cases. | PASS for exercised inputs |
| D-02 | `S08/focused-browser-final.log`: one/two/duplicate/DDBJ records, incomplete pair, invalid/Retry/Replace/Remove, delayed completion, and mode transitions; exact F1 summary and filename privacy assertion. S03/S07 History observations are historical after integration. | PASS for current focused cases; History reuse conditional |
| D-03 | Saved preview deferred/Inspect path in PR browser contracts and `S07/observe-sessions` for prior source. Real Linear Load in `S08/real-session/observation.json`. | Partial; specific Worker 1 accepted only for real non-default saved config |
| D-04 | `S08/d04-browser.log`: 3 current-head PASS for 1280 px and 390 px one-record disclosure, keyboard/focus/manual close, grid, and explicit batch. | PASS for exercised current cases |
| S-01 | `S08/focused-node.log` and `S08/focused-browser-final.log` cover JSON/gzip, current/historical/CLI/settings-only active mode, draft/committed divergence, and failed import rollback. Exact replay fails below. | Partial / replay FAIL |
| S-02 | `S08/focused-browser-final.log` Save single-flight/settlement and failed-import rollback, focused Node unsafe keys, PR browser rollback. The full comparison contract retains rotation and threaded Save failures. | Partial / two comparison FAILs |
| S-03 | `S08/vibrio-performance.log`, `S08/real-session/observation.json`; S06 multi-run stage/heap/heartbeat evidence applies only to its source/fixture/environment. Native clone wire/copy bytes remain UNAVAILABLE. | Incomplete |
| S-04 | Vibrio source/saved semantic hashes and CLI cross-surface PASS in `S08/vibrio-performance.log`; real Load/Save/Generate PASS in `S08/real-session/observation.json`. Current real-data field/strict-SVG equivalence was not rerun, and `tests/test_run_info_exact_replay.py` FAILs. Whole Session equality is not inferred from SVG parity. | Incomplete / exact replay FAIL |
| S-05 | Current Vibrio Save max 358.3 ms PASS against unchanged 500 ms limit. S06 real-data max-heartbeat observations remain FAIL under the user’s narrow allowance; S08 real Load/Save/Generate probe did not measure heartbeat. | Vibrio PASS; real historical FAIL/current unmeasured |
| A-01 | `S08/architecture-base-final.log` and postcommit `S08/architecture-head.log`, plus separate diff review. Trusted-base Gate PASS; Review REQUIRED. | Gate PASS, human review pending |
| W-01 | Branch, upstream, non-force push, and local/remote SHA in postcommit `S08/remote-receipt.json`. | Commit-time pending; inspect remote receipt |

The scenario map covers Circular/Linear, single/grid/batch, native/restored/inactive/settings-only, draft/committed divergence, History/rollback, and full replay through the listed focused/PR tests and historical scoped evidence. A listed historical case is not promoted to current PASS when source or condition differs.

## Verification and artifacts

| Command / artifact | S08 result |
| --- | --- |
| Focused Node Session/discovery/error tests | `S08/focused-node.log`: 16 PASS after the same-row fix. |
| Full Web Node | `S08/full-web-node-final.log`: 1,111 PASS on the final runtime/test source (352.973 s). Earlier `S08/full-web-node.log` also had 1,111 PASS before the same-row fix and is historical only. |
| `pytest tests/ -m 'not slow'` (native `--basetemp`) | `S08/pytest-not-slow.log`: 6,682 PASS, 17 skipped, 11 deselected, 2 FAIL. Comparison browser wrapper had four failing subcases before the same-row fix; its post-fix direct run has 14 PASS / 2 FAIL. Exact SVG replay has a separate sub-pixel mismatch. |
| D-04 narrow viewport browser contracts | `S08/d04-browser.log`: 3 PASS on the final source at 1280 px and 390 px, including grid and explicit batch. |
| PR browser contracts | `S08/pr-browser-contracts-final.log`: 19 PASS on final S08 source with native Playwright output directory. Initial NTFS output run had an unrelated trace ENOENT, retained separately. |
| Complete-record comparison contracts | `S08/comparison-contracts-after-fix.log`: 14 PASS / 2 FAIL after the same-row fix. Both all-record same-row cases now pass. Record rotation falls back to Python BLASTP generation because its committed generated-protein recipe lacks a typed resource; the browser Python runtime has no LOSAT/BLAST+ binary (`S08/rotation-debug.log`). The threaded Save case intermittently lost its execution context; a later unchanged isolated run passed (`S08/threaded-isolated-after-diag.log`), so the full 16-case invocation remains FAIL. |
| Read-only output references | `S08/output-comparison.log`: 16 PASS; references unchanged. |
| Guide/capture Python contracts | `S08/capture-contract-pytest.log`: 104 PASS on final capture recipes. These static contracts do not replace the three failed screenshot `--check` runs. |
| `ruff check gbdraw/` | `S08/ruff.log`: PASS. |
| Vibrio performance | `S08/vibrio-performance.log`: 1 PASS, Save 5.141 s, heap delta 145,727,423 bytes, max heartbeat 358.3 ms, no diagram Worker construction. This does not override the real-data 500 ms FAIL. |
| Reconstructed real full-pairwise Session | SHA `1a89693457e8bbe56a99eb2565e3f7a45e597d8aa6eaa47b802415aa98d05808`, 134,471,286 gzip bytes, 12 records. `S08/real-session/observation.json`: Load 32.954 s / 12 records, Save 17.805 s / 133,519,086 gzip bytes / SHA `67f6ad5efd54f6eb584fe25619c03ccb85284961e1efd55bc28266f71829e39b`, Generate 209.191 s / status `ok` / one Result / no page errors. The initial Load built one diagram Worker for the non-default saved config; the same Worker performed Generate. The original gzip SHA `d3cafef9664bff958aec4932eea6e264ae872cce89762a7cbe617447af975388` is unavailable. |
| Screenshot recipes | `S08/capture/`: T-GUI-05 and H-GUI-16 `--check` PASS; T-GUI-06, T-GUI-10, and T-GUI-12 `--check` FAIL on freshly captured finished-figure pixels. Their generators complete, inspect SVG semantics, and retain source figures; direct raster comparisons are in `S08/capture/` and native candidate paths. No screenshot was manually edited and no screenshot threshold was changed. Fresh-capture differences are local to the figure in sampled pairs: T-GUI-06 first diagram 210 pixels / max channel delta 72, T-GUI-10 result 3,612 / 255, T-GUI-12 result 1,553 / 192 (`S08/capture/raster-diff-summary.json`). T-GUI-10’s obsolete floating-search drag was removed after the search became docked; its remaining mismatch is within the figure. Failed `--check` logs are retained. |
| Architecture | `node tools/check-web-change-budget.mjs --base origin/dev`: Gate PASS, Review REQUIRED; no hard violation. Exact-head report is `S08/architecture-head.log`, produced after commit. |

The strict SVG replay failure reproduces in isolated pytest and after forcing its CLI to import the S08 checkout (`S08/replay-focused*.log`); 44 SVG attributes across 31 elements differ (max numeric delta 6.83e-13; `S08/replay-svg-diff.json`), although the native CLI output retains its pinned SHA. No test comparison or output reference was relaxed. F3 and this replay failure are separate from the two user-accepted observations.

Source/input/environment fingerprints: `S08/source-input-environment.json` (Python 3.13.3, Node 26.8.2, Chromium 149.0.7827.55, Pyodide 0.29.0). The regenerated public screenshots and their owner recipes were reviewed separately from production and tests; finished T-GUI-05/06 figures show labels, legend, quantitative tracks, and comparison context at readable scale. No tracked reference output, social preview, `dist`, `gbdraw.egg-info`, or generated browser wheel changed.

## Remaining boundaries and handoff

S08 does not declare complete integrated acceptance while the exact replay and guide capture checks fail, native structured-clone bytes remain unmeasured, and ordinary architecture review is pending. Reproduce the failed comparison, replay, and capture gates from their exact logs and fix their causes without changing pinned SVG output, acceptance thresholds, or mapped contracts. The next maintainer should inspect the final SHA's policy report, remote receipt, PR required checks, and this matrix before considering merge. No exception packet is generated because no positive OE/PE/CB exception condition was identified; `Review REQUIRED` is an ordinary human review signal, not an exception approval.

Proposed commit title: **Complete input and Session regression acceptance for issue 597**. Summary: Integrate trusted `dev`, repair discovery/Session errors and empty comparison projection, refresh public guide captures, and record the integrated gates and remaining failures.
