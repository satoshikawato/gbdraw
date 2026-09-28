# Issue #597 S06 follow-up result — 2026-09-28

Status: **local follow-up implementation and verification completed; S-03 remains INCOMPLETE** because native structured-clone wire/copy bytes are still `null` / UNAVAILABLE. Heartbeat max ≤500 ms remains FAIL where shown below. The user allowed only that miss and the one Python diagram Worker used by real Linear saved-preview Load with non-default saved config. S07, S08, and BUG-01 were not started. No commit, push, PR, remote merge, or deployment was performed. This document supplements [S06_RESULT.md](./S06_RESULT.md) and does not replace its recorded measurements.

## Source, refs, and preservation

- Checkout: `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovered-20260928`. Raw evidence: `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovery-evidence-20260928/S06-followup-20260928/` (below: `F/`). No work or evidence was stored in `/tmp`. Native-Linux pytest scratch was used only for case-sensitive temporary files.
- Intake state: HEAD and same-named remote were `ee6b6207845cc861ae428b4890a091509e1736b6`, and `MERGE_HEAD` was `007388567222638b707fbb16fe82dbeba61551c9`. There were 193 staged, 58 unstaged, and 21 untracked paths. The staged patch is byte-identical to the S06 preflight (`F/preflight-*`).
- Remote refs: actual `dev` was `ff5b58d716d1b14dea3e48689ef3045c95b6ce44` at intake and at final reading (`F/remote-refs-intake.txt`, `F/remote-refs-final.txt`). Its only advance after `8f2d8ed` is Issue #599 documentation. The Product Contract is byte-identical, so nothing was integrated. The fetch updated only `FETCH_HEAD`.
- The shared checkout and the Issue #619 changes were untouched. Pre-edit copies of every follow-up file are in `F/pre-followup-sources/`.
- Fixture: the reconstructed `real-full-pairwise.gbdraw-session.json.gz` with SHA-256 `1a89693457e8bbe56a99eb2565e3f7a45e597d8aa6eaa47b802415aa98d05808`. It holds 12 DNA records, 62,089 biological features, and 144 full pairwise entries. It contains 198 MB of LOSAT raw text, a 118 MB resources section, and a 103 MB editorState section. The original gzip is lost, and no PASS is claimed for it. Whole Session equality remains false.
- Environment: Chromium headless 149.0.7827.55, Pyodide 0.29.0, Playwright 1.61.0, Python 3.13.3, WSL2 with 31 GB RAM. The browser wheel `4382ac20…` matches all 215 checkout Python files. Other Playwright jobs ran on the host during some probes and are noted where they matter.

## Wasm limit: measured owners and risk

A diagnostic probe (`F/probe-wasm-stages.py`) served the checkout unchanged except for test-served wrappers. It recorded Pyodide `HEAP8.buffer.byteLength` before and after each Worker message, each Python helper, and each render phase. Linear memory never shrinks, so each reading is the high-water reached so far. The reported heap limit is **4,294,901,760 bytes** (`getHeapMax`). RSS, V8 heap, and transfer bytes were not used as Wasm proxies.

| Scenario (same Worker) | Before follow-up | After K1 |
| --- | --- | --- |
| After saved-preview Load | 90,439,680 | 90,439,680 |
| First Generate high-water | 1,819,934,720 / 1,821,704,192 / 1,919,418,368 (S06: 1,919,352,832) | 1,694,760,960 ×3 |
| `convertLosatpPairsToGenomicPayload` peak | 1,537.5–1,537.7 MB (3 runs) | 1,118.8 MB (3 runs) |
| Same request Generate ×2 more | unchanged (1,819,934,720) | not repeated |
| Collinear `minAnchors` change after Generate | 2,335,506,432 | 2,349,531,136 – 2,349,727,744 |
| Three consecutive changes (G2–G4) | G3 crashed the page (see below) | 2,349,727,744 each; plateau; headroom **1,945,174,016** |

Sub-steps of the baseline convert, from `F/wasm-stages-substeps-1`:
- 187.7 MB → 270.3 MB: `pairs.json` load
- → 389.3 MB: parsing 144 tables
- → 863.2 MB: `build_orthogroup_collinearity_blocks_from_hits` (+474 MB)
- → 1,118.8 MB: `encode_canonical_typed_resource` (142,494,014 bytes)
- → 1,537.7 MB: parsing that canonical result back into Python objects to `json.dumps` a 155,962,190-character cache string

After one Generate, the Worker retains the following:
- The LOSATP converted cache (one entry) and the filtered cache (63 entries, 39,051,330 deep bytes).
- The prepared-input cache. Its size estimates are decoded 142,494,014, parsed 228,761,165, resolved 228,761,165, and interactive 315,158,174 bytes; these are estimates, not object sizes.

Remaining concrete Wasm risk: a derived Collinear/Similarity recomputation after a render runs its working set on top of about 1.2–1.3 GB of retained Worker state. It reaches about 2.35 GB here (55% of the limit). K1 does not change this path. The high-water is bounded for this fixture and did not ratchet over three changes. The working set grows with all-scope pair count, so larger record sets remain unmeasured. The prepared-input cache is a documented one-project cache, and changing its bound would trade measured render time; it was not changed.

## Main-thread crash found during repeated Generate

Repeated Generate crashed the page twice before the fix. The first was a renderer loss in `wasm-stages-baseline-1-interrupted`, cause unrecorded. The second was a page crash in `wasm-stages-param-repeat-k1-1` G3, with page V8 used at 3.97 GB in the last sample while Wasm was 2.35 GB and the Worker V8 heap was 0.08 GB. A Wasm growth failure raises `MemoryError` rather than crashing the page, so this crash was not the Worker Wasm limit.

The owner was located with a V8 sampling heap profile after diagnostic GC (`F/page-heap-owners-1`) and a WeakRef-plus-heap-snapshot retainer path on Vibrio (`F/retainers-vibrio-1`):
- Used heap after GC grew 1.32 → 2.12 → 2.90 GB over two Generates. `jsHeapSizeLimit` is 4,395,630,592.
- Each Generate kept about 0.8 GB alive. This included about 190 MB of base64 canonical comparison resource and about 170 MB of catalog-admission clones.
- Retainer path: `latestCliHelperFiles[0].buildData` → closure `buildCanonicalReplayText` context → `committedArtifactHandle` (the pre-Generate rollback handle) → `runtimeState.canonical.committedCanonicalSession.renderRequest`.
- Each handle's runtime state holds the previous helper files, so the chain retained every earlier artifact.

## Implemented changes (production +72/−60 lines, `F/followup-production.patch`)

1. **K1 — `app/python-helpers.js`, existing LOSATP converted cache.**
   - A converted entry is now the summary JSON (pairs, provenance, cache) plus the canonical `(kind, bytes)` pair.
   - One helper, `_web_losatp_payload_json`, assembles both the fill output and the hit output:
     - With a canonical path, it writes the encoder's bytes and appends `canonicalResource`. The Web text is identical to the previous `json.dumps` result.
     - Without a canonical path, it inlines the same canonical JSON value. This is the Python-test/API path.
   - Removed paths: `json.loads` of the canonical result before cache serialization, and the cache-hit re-parse and re-encode of the canonical result.
2. **F1 — `app/run-analysis.js`, existing deferred canonical replay builder.**
   - The builder is created by a module-level `createCanonicalReplayTextBuilder`. The Generate scope passes it only the published request, the resources, and a holder.
   - The holder is filled with `results` and `featureCatalog` immediately after the candidate is destructured. Helper files are built after that and published only on activation, so lazy semantics, content, and failure isolation are unchanged.
   - Removed: the inline closure that captured the Generate scope.

Architecture evidence (ordinary, non-increasing):
- Both changes are private decomposition inside existing owners.
- They add no module, persisted format, compatibility path, Worker path, validator, renderer, transport, full graph clone, extra base64 conversion, or fallback. OE/PE/CB changes are 0/0/0.
- The trusted checker against `ff5b58d` gives **Gate PASS / Review REQUIRED**, with 0 import cycles and no privileged change (`F/architecture-gate-ff5b58d.log`). The review requirement reflects the inherited S05/S06 scope.
- Rollback: revert the two hunks listed above.

Product Impact classification is `IMPLEMENT_EXISTING_AUTHORITY` with `NO_USER_VISIBLE_DIFFERENCE`:
- Outputs are unchanged, and memory is the only difference.
- PD-OI-022 reuse, PD-OI-016/OIPC-C07 failure isolation, PD-OI-045, and the Circular admission and rollback behavior keep their existing owners.

Tests (+33/−5, `F/followup-tests.patch`):
- `tests/test_protein_colinearity.py` adds Web canonical-path miss/hit parity for the canonical bytes and the payload. The assertion passes against both the old and new helper (`F/canonical-parity-test-against-pre-followup.log`).
- `tests/web/run-info.test.mjs` failed on the S06 native decoder, both with and without this follow-up (`F/run-info-failure-pre-followup.log`), because it observed decodes through `atob`. It now counts the resource owner's `base64DecodeCount` structural metric. The asserted counts (lazy 1 → 2, native 0) are unchanged. The test is not a mapped contract.

## Results after the changes

- **Repeated Generate** (`F/wasm-stages-param-repeat-k1f1-1`): four Generates with three consecutive `minAnchors` changes all returned `ok`. The page did not crash, Wasm plateaued at 2,349,727,744, and the sampled page V8 peak was 3.61 GB (4.05 GB before). With diagnostic GC (`F/page-heap-owners-f1-1`), used heap was 1.32 → 1.55 → 1.77 → 1.77 GB. On Vibrio, the previous request, catalog, and Result are collected (`F/retainers-vibrio-f1-1`).
- **Science and continuation** (`F/final-full-pair-generate`):
  - Generate returned 12 records with 30,993 rendered and 62,089 biological features.
  - A CDS visibility edit, Undo, Redo, Save, fresh Load with its override, and export all completed. The fresh export exactly matches the edited export.
  - Strict XML: regenerated and Undo SVGs match the pinned CLI SVG, and the fresh SVG matches the edited SVG.
  - The generated Save equals S06's generated Save in every section except `createdAt` (`F/generated-save-vs-s06.json`). The audit result is identical to S06.
  - The CLI helper ZIP from the old and new `run-analysis.js` has the same entries, and the replay JSON is equal except `createdAt` (`F/cli-replay-parity-vibrio-3`).
- **Saved-preview equivalence**: all three final saves match the pinned CLI SVG under strict XML comparison, and all added-metadata checks pass (`F/final-saved-strict.json`, `F/final-added-metadata-validation.json`). Whole-content equality stays false.
- **Cancellation**: cancel after `session-candidate-prepared` returns `canceled` with identical canonical, request, Result, and History and the pending flag cleared (`F/final-private-cancel.json`).
- **Browsing while pending** (`F/browsing-observation-final.json`): during Load and during Save, scroll, search input, zoom, and animation frames all advanced, with 0 page errors and 0 external requests. The Save search status stayed 0/0, which proves input and paint, not a changed search result.

## Save process RSS

The final code does not change the Save path.

Owners from `F/save-owners-profile-1` and `F/save-owners-baseline-2`:
- Renderer RSS rises from 3.23 to 3.52 GB. The page ArrayBuffer backing store rises from 0.25 to 0.63 GB (S06 runs: 0.22 → 0.82 GB), because JSON and gzip chunks pass through native streams and are freed only by a major GC. The browser process rises from 0.10 to 0.29 GB during download handoff.
- The main-thread CPU profile shows native compression and stream work at about 13 s, 6.4 s idle, about 1.5 s of JSON generation, and 0.3 s of GC.
- Reusing chunk buffers would corrupt data behind a JS identity TransformStream (the test harness uses one), so it was not tried.

Measured and rejected candidates (one run each, compared with baseline 4.03 and 4.06 GB):

| Yield variant | RSS peak | Compression | p95 |
| --- | --- | --- | --- |
| `MessageChannel` | 4.85 GB | 14.0 s | 166 ms |
| `scheduler.yield()` | 4.83 GB | 14.3 s | 177 ms |

Less idle time delays GC. The earlier Blob-collector, 256 KiB-decode, and 64 KiB-JSON candidates were not repeated.

Final three-run pipeline (`F/final-pipeline-3/metrics.json`, same recipe as S06):

| Fixture / operation | Wall s | Heartbeat p95 ms | Heartbeat max ms | RSS peak GB | ≤500 ms |
| --- | --- | --- | --- | --- | --- |
| Real Load | 31.57 / 31.51 / 31.16 | 107 / 109 / 109 | 4,230 / 4,643 / 4,375 | 4.037 / 3.620 / 4.004 | FAIL (allowed) |
| Real Save | 18.09 / 16.74 / 18.10 | 128 / 128 / 114 | 693 / 738 / 850 | 4.236 / 4.249 / 4.030 | FAIL (allowed) |
| Vibrio Load | 6.60 / 7.60 / 6.72 | 172 / 165 / 150 | 1,144 / 1,174 / 1,155 | 1.394 / 1.347 / 1.393 | FAIL (allowed) |
| Vibrio Save | 4.07 / 3.95 / 3.97 | 142 / 184 / 150 | 178 / 196 / 229 | 1.584 / 1.556 / 1.564 | PASS |

Real Save RSS of 4.03–4.25 GB overlaps S06 (4.22–4.33 GB) and S05 (4.29–4.46 GB). No Save memory reduction is claimed. The shorter walls than S06 come from host conditions on unchanged Load/Save code, not from these changes, so no speed improvement is claimed. Real Load constructed one diagram Worker (allowed), and Vibrio constructed zero.

## Verification, failures, and reuse conditions

- **Python**: helper tests (`test_collinearity`, `test_protein_colinearity`, `test_reverse_complement_feature_identity`) gave 402 passed, 1 skipped. `test_session_io` gave 224 passed. `ruff check gbdraw/` passed.
  - An earlier helper run on `/mnt/c` failed 1 LOSAT-path test because that filesystem is case-insensitive (`F/helpers-pytest.log` history). The native-Linux rerun passed.
  - One combined run was killed at session end with exit 137. It was rerun as the separate passing runs above.
- **Node**: 25 focused files with 139 tests all passed (`F/focused-node-final.log`), including the Session, base64, run-info, run-analysis, and LOSAT cache tests. The full Web Node suite, the full pytest suite, and Node Playwright specs were **not** rerun. S06's Node Playwright launch failures are unchanged and not called PASS.
- **Interrupted diagnostics** are kept as recorded:
  - `wasm-stages-baseline-1-interrupted` was killed by a mistaken `pkill` that matched the probe itself.
  - `cli-replay-parity-vibrio-2-recipe-error` had a route-handler argument bug that was fixed in `-3`.
- **Remaining**:
  - Native structured-clone bytes are UNAVAILABLE.
  - Heartbeat max FAILs, and the lost original gzip is not measured.
  - Wasm on inputs larger than this fixture is not measured.
  - At the steady state, the page still holds two about-190 MB base64 copies of the canonical comparison resource. This is bounded and was not investigated.
  - Human review is still required.
- **Reuse conditions**: S06 evidence may be reused only where its sources are unchanged. `run-analysis.js` and `python-helpers.js` changed, so the Generate, continuation, CLI replay, Wasm, and pipeline evidence above was re-measured on the final sources recorded in `F/source-and-evidence-manifest.json` (`python-helpers.js` `0296050d…`, `run-analysis.js` `311315717…`).

Proposed commit title: **Bound Worker comparison memory and release prior Generate artifacts**

Short summary: **Keep canonical comparison results as encoded bytes in the Worker conversion cache, stop the CLI replay builder from retaining earlier Generate artifacts, and verify Wasm headroom, repeated Generate, Save memory owners, and scientific equivalence on the reconstructed real fixture.**
