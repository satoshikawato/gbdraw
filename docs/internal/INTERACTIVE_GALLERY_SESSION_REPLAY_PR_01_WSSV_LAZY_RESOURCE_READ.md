# PR 1 instruction: WSSV lazy session resource read

- Parent plan: [Interactive Gallery session replay remediation plan](./INTERACTIVE_GALLERY_SESSION_REPLAY_REMEDIATION_PLAN_2026-08-22.md)
- Suggested branch: `fix/wssv-lazy-session-resource`
- Base: latest `origin/dev`
- Dependency: none; PR 1 may proceed in parallel with PR 2
- Scope: one Web owner-boundary correction plus focused and real-browser regression evidence
- Gallery artifacts and schema versions: no change expected

PR 1 must establish the final isolated WSSV full-browser owner. PR 2 may invoke it from CI and add its shared fingerprint assertion in place, but must not move, copy, or delete the WSSV full-render case.

## 1. Required outcome

The bundled WSSV session must Generate from its 20 embedded FASTAs while they remain frozen metadata-only lazy resource views.

`services/file-content-cache.js::readFileText()` remains the only FASTA readability owner. `run-analysis.js` may choose an explicit string fast path; otherwise it delegates to that reader, which already owns lazy views, native File/Blob input, File-like `.text()` input, and invalid-input failure.

Do not:

- add payload methods or eager data to `SessionResourceFileView`;
- convert adopted resources to native `File` objects or add another resource-type dispatcher;
- request external FASTAs, rerun LOSAT, or hide a cache miss with fallback execution;
- change parsing, sequence transforms, session/request schemas, or Gallery artifacts.

## 2. Baseline and production change

Start from the current base and preserve unrelated work:

```bash
git fetch origin
git switch --no-track -c fix/wssv-lazy-session-resource origin/dev
git status --short

node --test \
  tests/web/session-resource-backing.test.mjs \
  tests/web/file-content-cache.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs \
  tests/web/run-analysis-simple-path.test.mjs
python -m pytest tests/test_wssv_gallery_fastas.py -q
```

Using the final regression rather than a disposable harness, record that import exposes 20 FASTAs as views without own `text`, `arrayBuffer`, `data`, or `resourceId` fields, then Generate fails at `Input file is not available for browser FASTA extraction.` Do not log sequence or base64 content. If current `origin/dev` no longer reproduces this boundary, reconcile the plan with the merged implementation instead of adding a duplicate fix.

In `gbdraw/web/js/app/run-analysis.js`, remove only the `file?.text` preconditions from:

- `extractLosatFastaFast()`;
- `extractAllLosatFastaFast()`.

Both paths must reduce to:

```text
sourceText = explicit text when it is a string
             otherwise await readFileText(file)
```

Do not replace the guards with another type/capability branch. Preserve GenBank/FASTA parsing, record selection, region and reverse-complement transforms, record ordering/IDs, and canonical-length calculation. Let `readFileText()` own low-level invalid-input errors; add stage context only once at an existing Generate orchestration boundary and never include payload data.

Expected unchanged files include:

- `gbdraw/web/js/services/file-content-cache.js`;
- `gbdraw/web/js/services/session-resource-backing.js`;
- the WSSV session, SVGs, and thumbnail;
- session/request version constants.

If either shared service appears insufficient, stop and identify the missing existing contract before expanding scope. This fix should remove two duplicated decisions, not add a compatibility path.

## 3. Regression evidence

### 3.1 Focused owner tests

In `tests/web/wssv-gallery-fastas.test.mjs`, remove the synthetic `withText()` wrapper. Use `adoptCurrentSessionResources()` and give its `sessionResourceTable` to `projectCanonicalSessionRequest()` so `buildRestoredMatchSequenceSources()` consumes the real frozen views.

Retain evidence for:

- one circular reference and 20 ordered homology sources;
- metadata-only views and decode only on actual read;
- filename, source index, record ID, canonical length, hash, and provenance.

Do not weaken the exact 20-resource checks in `tests/test_wssv_gallery_fastas.py`.

In `tests/web/run-analysis-simple-path.test.mjs`, or the nearest existing `createRunAnalysis()` owner suite, drive the public orchestration boundary through both private extraction paths:

- Circular conservation uses lazy reference/comparison views;
- Linear sequence preparation uses a lazy view with selector, region, and reverse settings;
- explicit string input performs no file materialization;
- missing/invalid input propagates the shared-reader failure.

Keep corrupt-base64 and checksum-mismatch coverage at the shared resource boundary. Assert observable request/executor inputs and resource-read metrics; do not export private helpers or duplicate their parsing logic.

### 3.2 Final isolated browser owner

Create these retained test surfaces in PR 1:

- `tests/web/contracts/wssv-full-generation.serial.spec.js`, the sole WSSV full-render owner;
- `playwright.wssv.config.js`, modeled on the isolated Vibrio config with `fullyParallel: false`, `workers: 1`, failure-retained trace, and a measured timeout with headroom;
- npm script `test:web:wssv-generate`.

Do not also place WSSV full generation in functional smoke, `gallery-session-regeneration.playwright.spec.js`, or a common Gallery replay spec. PR 2 must call the dedicated script from its WSSV CI job instead of relocating the test.

The spec must use a fresh browser context and the real application lifecycle:

1. import `WSSV_genome_comparison.gbdraw-session.json` through `importSession()`;
2. assert successful saved-preview restore, metadata-only 20-FASTA order, and no diagram Worker construction during import;
3. install a probe that rejects unexpected external/cache-miss LOSAT execution;
4. call `app.runAnalysis()` and await its settlement;
5. assert success, no error log, one committed current Result, and absence of the old missing-file message;
6. assert 20 ordered conservation groups/rings with expected labels and colors, plus sanitized interactive SVG admission.

Do not create a temporary SVG comparator in PR 1. PR 2 adds initial/generated parity through its single shared comparator in this same spec; the test location and full-render body remain stable.

## 4. Verification and acceptance

Run narrow gates first:

```bash
node --test \
  tests/web/session-resource-backing.test.mjs \
  tests/web/file-content-cache.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs \
  tests/web/run-analysis-simple-path.test.mjs
python -m pytest tests/test_wssv_gallery_fastas.py -q

python tools/prepare_browser_wheel.py
npm run test:web:wssv-generate
```

Then run the broader gates justified by changing shared Generate orchestration:

```bash
npm run test:web:functional-smoke
node tools/check-web-change-budget.mjs
python -m pytest tests/ -v -m "not slow"
```

Allow the browser and Python suites at least 30 minutes. If Chromium hits the host sandbox restriction, rerun the same browser command with escalation. Record exact commands/results and WSSV elapsed time; elapsed time is evidence, not a performance claim without a threshold.

Accept PR 1 only when all of the following are observed:

- [ ] the final tests fail against the two pre-fix guards and pass after their removal;
- [ ] import stays lazy and Worker-free, and exactly 20 embedded FASTAs remain ordered and hash-verified;
- [ ] Generate succeeds with all 20 ordered rings/labels/colors and no selection, external request, or unexpected LOSAT execution;
- [ ] `readFileText()` is the sole FASTA readability owner and invalid inputs remain safe and clear;
- [ ] `wssv-full-generation.serial.spec.js` is the only WSSV full-render owner;
- [ ] no Gallery artifact, generated wheel, package-lock dependency, or schema version changes;
- [ ] every listed gate passes, or an unrun/failed gate remains explicitly unverified rather than marked complete.

## 5. Architecture evidence

Capability scope: browser text-resource readability used by FASTA extraction. It is complete because `readFileText()` owns supported representations and the only narrower WSSV-path decisions are the two named extraction guards; tests cover both callers, the real adopted representation, and the browser path.

```text
Owners before = {
  services/file-content-cache.js::readFileText,
  app/run-analysis.js::extractLosatFastaFast capability guard,
  app/run-analysis.js::extractAllLosatFastaFast capability guard
}
Owners after = {
  services/file-content-cache.js::readFileText
}

Canonical path before = FASTA extraction -> narrower .text guard -> shared reader
Canonical path after  = FASTA extraction -> shared reader

Compatibility paths before = none
Compatibility paths after  = none
Superseded paths removed    = both .text guards
New production modules/dependencies = none
```

The Playwright spec/config are verification infrastructure, not semantic owners. Expected reviewed result: owner excess decreases; path excess does not increase; compatibility burden is unchanged. Scientific-output evidence is the unchanged WSSV manifest plus the ordered 20-ring assertion; no reference refresh is authorized.

## 6. Handoff

Review production and test diffs separately and confirm `git status --short` contains no artifacts or unrelated edits. The PR body must record the architecture sets, exact gate results, elapsed browser time, and remaining risk.

Suggested commit title: `Fix lazy WSSV FASTA reads during generation`

Suggested summary:

```text
- route WSSV FASTA extraction through the shared lazy-resource reader
- verify both extraction paths and the real metadata-only resource views
- retain one isolated full-browser WSSV generation owner without artifact changes
```
