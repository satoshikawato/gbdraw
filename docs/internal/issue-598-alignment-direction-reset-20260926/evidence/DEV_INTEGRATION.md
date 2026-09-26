# Issue #598 direction / Apply / Reset integration into dev

Date: 2026-09-27 (Asia/Tokyo). This report owns the dependency integration requested for Issue #602 Session 06. It does not mark the separate Issue #598 S02/S03 writer sessions complete or start S04/S07, deployment, or release preparation.

## Provenance and branch boundary

- Base: `origin/dev@af5d942af60353dda199aa487da9152a3576b3fe`.
- Work branch: `integrate/issue-598-runtime-dev-20260926`, created from that base without tracking dev/main.
- Dedicated checkout: `/tmp/gbdraw-issue598-runtime-dev-20260926`.
- Published S01: `d5edc4bb00b0d793d94361202c31a270f78a4aa0`; retained by ordinary merge `6335fb905a60e8addf10f89413f44f444e2e7238`.
- Implementation inputs: read-only, stable snapshot of 29 uncommitted files in `/tmp/gbdraw-598-s02-20260926@0453380087cbae640cf0c3b9792a1b7d1d023213`. Each tracked baseline blob matched the integration candidate before copying; every source file remained unchanged throughout the snapshot. Source work was not described as a completed session. The source writer's checkout, branch, processes, and remote branch were not modified.
- Snapshot manifest SHA-256: `cc89f48bee4d8b528d7282bbdbb8eecd5c758a8dd3b512d652fe4df6194d0375`. The exact input hashes are below. Subsequent integration fixes are independently owned by this candidate.
- The shared checkout and `fix/issue-602-linear-live-edit-20260926@94faf8a98eddf823f5ebfed859d324e258daa689` were preserved. The dependency is published on its own same-named branch and integrated through a normal PR to dev. Session 06 must synchronize from the resulting dev commit rather than importing this working branch directly.

## Current Product authority

The formal base contract remains revision 21, SHA-256 `62e3a9c08ceb64acc349ace81a97a9c87dc45db4d18187a3b78e6e131b2dc595`. Every receipt field (concern, scenarioRevision, choice, rationale, mustPreserve, mayRetire, acceptedResidualRisk, owner, decisionDate) for the four Issue #598 decisions matched origin/dev: PD-OI-027/029/031/034 revisions 5/3/5/5. PD-OI-035 revision 3 and PD-OI-039 revision 2 are independent additional requirements. The latter selects `EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`; the historical revision-1 Match acceptance was not rewritten or restored. The supplied historical intersection packet records the owner's revision-2 receipt already present in dev; it confers no candidate runtime approval.

Preflight: `IMPLEMENT_EXISTING_AUTHORITY`. Existing Product choices are realized; no new Product decision, exception, compatibility waiver, or self-authorization is requested. Other formal contract bytes are unchanged, including PD-OI-038 and PD-OI-044/045. Issue #602 S01–S05 remain on their original work branch and are not included in this dependency PR.

| Required contribution | Verified behavior |
| --- | --- |
| Explicit exclusive directions | Keep, all-right, all-left, Custom; Custom-only per-record Keep/right/left; exact reference participates. Real Gallery right mode reverses only the minority reference, left mode reverses the four selected targets. |
| Local intent and identity | Select/Skip and Custom changes, including keyboard and width transitions, retain exact source identities, local choices, before/after arrows and unchanged reasons; additional Worker posts are zero. Source strand is not rewritten. |
| Final validation and retry | One final helper batch per Apply; material facts changing the preview require another explicit Apply. Validation/render failure retains draft and underlying error, previous artifact and directions. Cancel/stale cannot commit. |
| Atomic artifact and History | Successful Align/Reset use the existing canonical candidate admission and one History transaction. Undo/Redo restores plan, directions, placements, receipt and SVG. |
| Reset | Default positions-only preserves directions. Combined Reset restores only actual changes from the latest Align, including reference; bound original evidence is consumed by either scope. Missing historical evidence and empty modern changes have distinct reasons. Failed Reset retains evidence for retry. |
| Persistence and geometry | Fresh compressed Save/Load preserves receipt and SVG; corrupt binding is rejected. Reset and manual Reverse use original comparison evidence with zero LOSAT jobs. Five record anchors remain centered with ribbons retained; measured offset <= 0.000123 px. |
| Review presentation | Same nonmodal review and SVG, compact dock with independent list and reachable footer, wide drag/focus/pan/zoom preserved. Existing CSS 40rem query is the only compact boundary. Editor closes through its owner, selected tab stays, reopening is disabled with reason during narrow review, and is explicit afterward. |
| Saved artifact boundaries | Review markers/draft are absent from saved SVG/session. PNG/PDF bytes agree with review open/closed after removal of format-only timestamps/identifiers. |

## Integration fixes and owner/path evidence

- `record-display-options.js` retains the existing committed source identity binding and directly resolves typed canonical targets for alignment intent. Fresh loaded sessions no longer require discovery or Generate before a reference-only direction Apply. Rollback restores pending orientation and existing display drafts; source replacement still rejects stale bindings.
- `normalize_current_session_artifacts` is the existing shared Python current-writer boundary. Both direct API writes and historical CLI rewrites initialize an absent alignment reset receipt to explicit null there. Readers still reject missing current-format evidence and malformed/stale receipts; historical evidence is never invented. The duplicate CLI-only initialization was removed.
- The same record transform owner, canonical session-request projection, lazy diagram Worker, Python typed planner/renderer, run-analysis Result admission and History remain the sole paths. The old Match controller/checkbox and run-analysis `canonicalStateOverride` overlay path are removed. No second renderer, decision owner, Result admission, review controller, compatibility reader, schema version, dependency, or visibility/tab ref is added.
- CSS compact status is read from the existing container query custom property through the existing lifecycle. Separate JS 640px classification was removed. Compact transitions terminate a wide drag; the compact header has no drag affordance. Drawer visibility-only controls do not create an empty History revision.
- Existing plain-text record-label formatting is reused; Custom controls have explicit record-bound accessible names.
- This is architecture-bearing, ordinary non-increasing owner/path work: reviewed changed capabilities retain one owner and one canonical path. OE and PE do not increase; CB is unchanged. No exception trigger or hard-invariant waiver applies. Persisted receipt validation is in the existing session contract, with only actual deltas and source binding, not direction policy or a general transform snapshot.

## Verification and reproducible commands

Environment: Linux/WSL, Python 3.13.3, Node 26.8.2; Node Playwright 1.61.1 and Python Playwright 1.61.0 both available. Real Chromium used an isolated checkout and free local test port; no shared editable install. The browser wheel was generated normally and remained ignored; cache-bust token unchanged. Sandbox host mount failure required the same local checks to run with escalation.

```bash
python tools/prepare_browser_wheel.py --no-build-isolation
node --test tests/web/*.test.mjs
GBDRAW_WEB_TEST_PORT=47200 node node_modules/@playwright/test/cli.js test tests/web/similarity-alignment-ui.playwright.spec.js tests/web/right-drawer.playwright.spec.js --workers=1 --output=/tmp/issue598-integration-final-browser-all
python -m pytest tests/test_api_session.py tests/test_session_compat.py tests/test_session_io.py -q
ruff check gbdraw/
python tools/gallery_artifact_manifest.py
node tools/check-web-change-budget.mjs --base origin/dev
# After commit:
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check
```

- Final fast Web: **711 passed**, including the architecture tests and the added genuine fresh-Load source-binding regression. Focused owner rerun: **82 passed**.
- Final current/legacy session regression set: **287 passed**. Source-bound Python planning/helper tests previously passed **371** and the untouched S01 integration baseline passed **275**.
- Full fast Python run completed naturally in 25 minutes: **6262 passed, 17 skipped, 11 deselected, 6 failed** before the writer fix. All six failures were current alignment fixtures or the missing shared writer initialization. The complete affected session files were then rerun successfully (287 above); no later whole-suite pass is claimed. Untouched render/reference and other test evidence is reused within its unchanged scope.
- Final focused browser result: **26 passed**. These journeys cover released v40/v44 inputs, real Gallery comparisons, Custom/Select/Skip, exact reference-only direction Apply, retries, History, Save/fresh Load, Reset, exports, and right-drawer continuations.
- Strict 390×740/844 geometry: canvas width 370px; visible heights **202px / 254px**, measured after actual header and overlapping Generate-bar intersections. Canvas pointer and both Apply/Cancel centers hit their intended targets. Internal screenshots were visually inspected. They are QA artifacts, not public showcases. The source recipe's labels, metadata, legend, five records and comparison ribbons remain intact.
- Lint, whitespace and Gallery artifact manifest: pass. Web policy: **Gate PASS / Review REQUIRED**; no blocker. Review REQUIRED reflects architecture/session changes and size, not an inferred approval.

Production, test, documentation and generated-data diffs were reviewed separately. Current Gallery session plus identical test input gain only `editorState.alignmentResetReceipt: null`; the manifest changes only that file's byte count and hash (2,890,972 bytes, SHA-256 `a32f880bf0281a8d9e47b475be1d0a2259f5d4b40ddad49ec054b6b58ee51465`). No figure, reference output, owner-maintained social preview, generated wheel, dependency/policy/guard or CI file is committed.

## Scope limits and continuation

This dependency validation does not replace Issue #602 Session 06's final matrix after its own S05 layout is synchronized. Compact Editor regression and compact review acceptance must be checked together there. It does not claim all staging, S04 cross-version/release readiness, real OS keyboard/physical browser zoom, or the separate known Legend-override retry issue is resolved. Public documentation is Session 07 and remains out of scope.

Proposed commit title: **Complete alignment direction and Reset runtime for dev**.
Summary: Integrate explicit local direction choices and source-bound Reset evidence through the existing atomic canonical render path; fix fresh-load Apply and current/legacy session writers, and verify browser, persistence and geometry continuations.

## Input snapshot hashes

| Path | SHA-256 |
| --- | --- |
| `gbdraw/api/record_planning.py` | `f603aa55e25021161996728336a8b4207512812ba3a88c1f7a880635ba91c016` |
| `gbdraw/api/request_render.py` | `b03cec3bc44143323cbdf562b629499afd7541c4ecf2195ae54761e03c51f203` |
| `gbdraw/session_io.py` | `5fb3a38dda507a502d5c201a0c80b886511fc59d2a779cbffa872e34e0957be3` |
| `gbdraw/web/gallery/artifact-manifest.json` | `00a627dfadba536aa2a5b3f582df24184bee9e3704099f5a83885ab33959146e` |
| `gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json` | `a32f880bf0281a8d9e47b475be1d0a2259f5d4b40ddad49ec054b6b58ee51465` |
| `gbdraw/web/index.html` | `3a5b1148ce70b85e66577a5264d48285fa968ee1ea23770158be1cdcf88c75ec` |
| `gbdraw/web/js/app/app-setup.js` | `f56deb4f0f693373f42929be811d373d12eb773c808fefd8044538c7367844d5` |
| `gbdraw/web/js/app/record-display-options.js` | `3794b98c4d7df163ac46a5d8ed6a8f2c1c17f87a2a55b7954ff5b02bdefb311f` |
| `gbdraw/web/js/app/right-drawer.js` | `190169e486a6f5fec7c88690767dea5bc5852395e9b2396c472f78567cbd547c` |
| `gbdraw/web/js/app/run-analysis.js` | `20b7de3202e38f5e2c201fd0468d4102a6f7dc65b03c74841c1fe5b26bcf1ffd` |
| `gbdraw/web/js/app/similarity-alignment.js` | `9cceb78ed792a71b22623fc8d31f0aef115808c8bd52d123535bf3679b33d751` |
| `gbdraw/web/js/services/config.js` | `1c96a5de66a06ebd9c4f821b0d74e4ebf948d1a10c6466550aa9085eb2d7566d` |
| `gbdraw/web/js/services/error-normalization.js` | `55a473698ac967c716f1888ac812436ad3d33442040eee7339a518928d600330` |
| `gbdraw/web/js/services/history-snapshot.js` | `997169a7d9353cdd0d733d274de4c70105099c5921aa10f8681891a72cc80f01` |
| `gbdraw/web/js/services/session-active-config-contract.js` | `7d12b467e91d89f2244cb9fc485569c8fd2282878708753adc56df8355a17ac6` |
| `gbdraw/web/js/services/session-authority.js` | `f9042293b5c4ec7c67491b46a7f1f3c897f7ae9c443ac0808d7ecbf9deb5bfdc` |
| `gbdraw/web/js/services/session-request.js` | `814938b39fb15640e4f20ab1fcd14a24ef2cee76811d7a0a47fa819d041badc6` |
| `gbdraw/web/js/state.js` | `cff284514953c26b7acf3e46de0ed1f5e159e900bffc172db6ab0fecece40018` |
| `tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json` | `c5d51aadd979565a3146dedbfe94e6a6d9c5948d19e3d0ee5d8bc57501168e9c` |
| `tests/test_record_planning.py` | `d7374d343210dc84d84656611165d0a59236020f9d19c94b20897609cbe9d72f` |
| `tests/web/error-normalization.test.mjs` | `05ed1e0a8ac2d92169f4f667aa5cdeeea34fdc59f4b9873ac39c508be2faa099` |
| `tests/web/history.test.mjs` | `09a40f445f387f3aa3ab56dc4d5731f5d8b9790e187607f07afea3eae3cd9391` |
| `tests/web/record-display-options.test.mjs` | `e001bc4ecb0e04163ae4616b00d90a22f328611b3219e63df5a3cfb4cc6ef576` |
| `tests/web/run-analysis-simple-path.test.mjs` | `a43989937c3d654076ddcd38b858624d37e9b3072cd21c54dd2bbf03bb4b8961` |
| `tests/web/session-losat-cache-validation.test.mjs` | `7f21c5c0c05ed131f903a17d453230b52dc94ed46c832b4d69e8b7d1ea889aba` |
| `tests/web/similarity-alignment-actions.test.mjs` | `8e1ddf41b5082d9c8819958f05d81a4ba3858c0833549aafe365f8f4763b18ff` |
| `tests/web/similarity-alignment-ui.playwright.spec.js` | `bd923d8723f674daac431e46a4387373c854c618f3079c9367c79d65175585dc` |
| `docs/internal/issue-598-alignment-direction-reset-20260926/decisions/05_REVIEW_MATCH_INTERSECTION.md` | `4f3635975403b31a21575ad5ebf7feb9595f78c5376fe9c2c68920919e27f0d1` |
| `tests/web/alignment-reset-receipt.test.mjs` | `adaf00855b558051090dfb5cce78fd2bbf46142c7fc5fdbf29f59a7b2f4b6bff` |
