<!-- Raw design report of workstream W8 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

I checked the test suite against every bug ID, read-only, in DEV (dev 4c89bab1), skipping tests/web/audit-*. None of the nearest tests would catch any of the listed bugs. Three reasons repeat across them:
- They use single-Result fixtures.
- Their oracles compare copies of the same artifact, never a live edit against a fresh Generate.
- They assert state maps or delta numbers, not the rendered outcome.

**Tier legend** (from playwright*.config.js and .github/workflows/test.yml):
- **smoke**: Playwright tagged @pr-smoke. Runs on PRs to dev (web-pr-smoke). Only about 15 tests carry the tag.
- **full**: Playwright functional-full. Runs only on push to dev or workflow_dispatch (test.yml:470-514), never on PRs. This includes tests/web/contracts/*.playwright.spec.js.
- **node**: tests/web/*.test.mjs. Runs on PRs (web-contracts-pr, when the ci-impact plan requires it) and on dev push (browser job).
- **manual**: tests/web/*.playwright.py. No config's testMatch picks them up and pytest does not collect them, so they are not in CI.
- **perf**: *.performance.playwright.spec.js. Runs only in the perf job on dev push.

Of all the nearest tests below, only right-drawer.playwright.spec.js:702 (desktop case) and mode-transition-result.playwright.spec.js:30 are smoke.

### Batch / grid question
- **No Playwright spec loads a batch, edits, and then switches Result.**
- The only batch + edit + switch coverage is tests/web/decoration-continuity.playwright.py:197-244. It sets up two batch Results, drags the legend/title on each via `selectedResultIndex`, runs Generate twice, then Session. It only does drags, and it is manual.
- Batch specs with no edit:
  - gui-audit-regressions.playwright.spec.js:143-175 switches Result with `selectOption('1')` during PDF/interactive export and checks only the export.
  - composite-session-resources.playwright.spec.js:166-200 ('grid-batch-grid') checks resources only; no edit, no switch.
  - circular-record-presentation.playwright.spec.js:224-258 checks control applicability only.
- Grid specs that edit, but where grid is one Result:
  - composite-placement-regeneration.playwright.spec.js:35-184 (placement and label edits, then Generate).
  - definition-replay-visual-state.playwright.spec.js:10-54 (composite case).
- Unit-level batch coverage:
  - candidate-render.test.mjs:161-168 (only `callerTransforms` per Result).
  - composition-layout.test.mjs:829-847 (batch decoration mapping).
  - feature-catalog.test.mjs:727-776 (batch record_idx rebasing).
  - preview-runtime.test.mjs:181-190 (Result 0 is flushed on `selectResult(1)`).
- `multi_record_canvas` otherwise appears only as a setup value in request/session tests.

### Class 1: multi-Result edits

**FE-01** (label overrides cleared by syncLabelEditor)
- mode-transition-editor-state.playwright.spec.js:16-62 (full): label text survives mode Undo/Redo and Save/Load/Generate, Circular and Linear. Asserts intent maps and that the SVG contains the text.
- Same file :65-89 (label visibility off, Circular only) and :141-169 (feature hide, then mode change and Generate, with no label override present).
- multipart-label-regeneration.playwright.spec.js:17-167 (full): tobacco, single Result. Checks selected, mounted and exported copies agree, plus label text after Generate.
- contracts/session-regenerate-intent.playwright.spec.js:1889-2122 with its edit helper at :1137-1265 (full). It sets label text, label visibility and a hide, then Generates. It only compares override maps between two Generates (:2010-2012, :2099-2101). The pre-Generate check (:1940-1948) is the only assertion that the label edits exist; nothing asserts they survive Generate.
- right-drawer.playwright.spec.js:702-834 (smoke): live label edits, checked in the DOM and in Result content, with no Generate afterwards.
- Missing: the triggers themselves. No test hides the labeled feature and Generates, switches batch Result, or round-trips the Circular record selector. `labelOverrideContextKey` is captured in helpers/mode-transition.cjs:91 but no test asserts on it.

**FE-02** (scoped edits apply only to the displayed Result)
- contracts/active-result-edit-transaction.playwright.spec.js:238-435 (full, retries 0): "Apply to all tRNA" through Undo/Redo, Save, fresh Load, Generate and export. It is single-Result HmmtDNA and only inspects `results[selectedResultIndex]` (:111).
- right-drawer.playwright.spec.js:702 (smoke): single Result.
- Unit tests (node), all single-Result:
  - svg-style-completion.test.mjs:6-10 (applySpecificRulesToSvg writes `results[0]`).
  - feature-visibility-actions.test.mjs:60, 81-110 (fixture has one Result, 'one.svg').
  - candidate-render.test.mjs:18-158 (admission `resultNames:['diagram.svg']`).
  - preview-runtime.test.mjs:181-190 checks that Result 0 is flushed on switch, not that Result 1 receives the edit.
- Missing: any batch Result 2 check after a scoped color, hide or legend edit, whether live or after Generate.

**FE-03** (Features drawer lists record 1 while Result 2 is shown)
- feature-catalog.test.mjs:727-776 (node): batch catalog gives record_idx [0,1] and the right featureRecordIds. It never ties `selectedFeatureRecordIdx` to `selectedResultIndex`.
- gui-audit-regressions.playwright.spec.js:143-175 (full): switches to Result 2 for export only; the drawer is never opened.
- Nothing tests `filteredFeatures` (state.js:711-725); a grep finds no reference. Every drawer spec in right-drawer.playwright.spec.js loads a single-Result session.

### Class 2: live edit differs from Generate

**IN-01** (updateDefinitionText uses the whole file)
- definition-replay-visual-state.playwright.spec.js:10-54 (full): edits the definition font on a CLI session (single source and a 2-source grid). Asserts group font-size = 19 and that selected, mounted and exported copies agree. It never compares against a Generate, and never uses a record selector, region or reverse complement.
- circular-record-presentation.playwright.spec.js:37-164 (full): record selection + reverse complement + region, but only Generate-vs-Generate equality (:148). It uses `labels_mode:'none'`, makes no live definition edit after Generate, and reads a hard-coded `results[0]` (:21-35).
- Same file :173-222: sets `def_font_size += 1` (:213) with a record selected, but no Result exists, so the code path returns early. Only the disclosure UI is asserted.
- definition-layout-completion.test.mjs:52-86 (node): stale-callback guard with a stubbed helper and fake `files.c_gb:{}`; the payload is never inspected. :88-117 covers watcher scheduling.
- tests/test_web_feature_metadata.py:137-231: helper tests with whole-file input only.
- helpers/retained-visual-journeys.py:175 (C01 definition counter) is manual; the archive it needs is not in the repo.

**GE-02** (stroke live edit writes invalid values)
- history-generated-authority.playwright.spec.js:180-205 (full): sets `block_stroke_width = 2` (a valid value) in Circular and asserts the mounted `stroke-width` is "2". It never checks Result content, never tries empty or invalid values, and never compares with Generate.
- decoration-continuity.playwright.py:136-142 (manual): sets block_stroke_color, then Generates, and asserts deltas only.
- session-request.test.mjs:2172-2186 (node) covers request projection, not live edits.
- No unit test calls `applyStylesToSvg`.

**PV-10** (Linear in-place legend move differs from post-Generate layout)
- composition-layout-real.playwright.spec.js:659-791 (full): Linear in-place side switches (:673-678, :733-745). Asserts validity, containment (`viewBoxContainmentError < 1` at :616) and preserved deltas, but never compares with a Generate at the same side.
- composition-layout.playwright.spec.js:316-325 (full): legacy session moved in place; metadata only.
- composition-runtime-parity.test.mjs:12-66 (node): JS replan vs Python oracle on synthetic metadata. It covers the placement plan, not the entry reflow in `reflowSingleLegendLayout`.
- linear-typography.playwright.spec.js:95 asserts the form value only. source-legend-reconciliation.playwright.spec.js:68 asserts entry colors after Generate.

### Class 3: live edits not carried by Generate

**PV-01** (dragged legend/title leaves canvas after Position change + Generate)
- composition-layout.test.mjs:778-800 (node) asserts the current behaviour. With a new automatic baseline of [300,400] it requires the delta to stay [40,20] and the legend x to equal 340 (:787, :796). It never checks that the legend stays inside the canvas.
- decoration-continuity.playwright.py:136-145 (manual): sets legend 'left', Generates, then asserts deltas are preserved and the automatic translation changed. No canvas containment check.
- composition-layout-real.playwright.spec.js:683-745 (full): drag, then an in-place side switch, with a containment check. No Generate after the drag; the Circular branch never changes side.
- run-analysis-simple-path.test.mjs:603, 847-870 (node): stubbed `captureDecorationContinuity` failure rollback.

**PV-02** (renamed featureless legend entries revert on Generate)
- feature-color-actions.test.mjs:739-746 (node): featureless rename in a fake DOM; asserts text and data-legend-key only. No Generate or plan.
- candidate-render.test.mjs:85-90 (node): the legendRename case assumes originalCaption 'CDS' differs from caption 'Genes'. It never models the overwrite of originalCaption at color-actions.js:666.
- source-legend-reconciliation.playwright.spec.js:130 (full): renames feature-backed 'rRNA' (rule path), which does survive.
- composition-layout-real.playwright.spec.js:346-380: renames a manually added entry; Save/Load, no Generate.
- No test renames GC content or GC skew.

**PV-03** (legend sort lost on Generate)
- color-captions.playwright.spec.js:98-111 (full): calls `sortLegendEntries()`, then Generate, but asserts only color via `.find` (:111, :157). Line 96 even sorts both arrays before comparing.
- legend-sync.test.mjs:385-409 (node) covers the originalLegendOrder inventory, not user sort.
- session-operation-consistency.playwright.spec.js:75 checks only the busy rejection.
- No test asserts legend order after Generate.

**PV-07** (canvas padding reset to 0 on Generate)
- No test sets canvas padding and then Generates.
- Nearest:
  - composition-layout.playwright.spec.js:84-88 (padding applied to a synthetic SVG).
  - history.test.mjs:995, 1115 (fake-state restore).
  - session-draft-authority.test.mjs:998 (setup only).
  - right-drawer.playwright.spec.js:198-205 (toggles the padding panel only).
  - session-operation-consistency.playwright.spec.js:76 (busy).

**CO-08** (group names bound to unstable og_* ids)
- gui-audit-regressions.playwright.spec.js:118-141 (full): rename and description, mode round trip, Generate, all with the same inputs, so the ids never change. Asserts the override maps are equal.
- orthogroup-computation-cache.test.mjs:12-58 (node) asserts that lookup fails after an id change, i.e. it asserts id-keyed identity.
- orthogroups-stable-identity.test.mjs uses fixed og ids and tests feature-to-group membership.
- right-drawer.playwright.spec.js:1279 (row rendering) and interactive-svg-v3.playwright.spec.js:1128 (override forwarding to export) are also near.
- Missing: re-inference that renumbers ids (changed inputs, order or selection). `pruneOrthogroupOverrides` is untested.

### Class 7: History

**SE-01** (checkpoint undo/redo then Linear→Circular throws)
- history-generated-authority.playwright.spec.js:69-125 (full): Undo/Redo of the Generate A/B checkpoint, Circular. Asserts selected/mounted/exported agreement and the request. No mode switch afterwards.
- mode-record-identity.playwright.spec.js:70-150 and mode-transition-result.playwright.spec.js:64-106 (full): Undo/Redo of mode switches or intent edits only, never a checkpoint Undo first.
- history.test.mjs:1052-1069 (node): fake-state checkpoint round trip. history.test.mjs:1850-1880 replaces `applyEditorStateData` with a fake and asserts `adoptCatalog === true`, so the real clone path in config.js:1169-1171 is never exercised. No test covers `featureStateFromCatalog` after a restore.

**SE-02** (checkbox/radio via label click not undoable)
- draft-placement-capability.playwright.spec.js:102-114 (full): `uncheck()` on the Separate Strands input itself, then Undo/Redo.
- history-inputs.test.mjs:58-90 (node, selects only) and history-inputs.playwright.spec.js:47-135 (full, select via keyboard and pointer).
- Missing: a click on the wrapping `<label>`/`<span>` (index.html:2878-2880).

**SE-03** (click while a text field is focused is dropped or merged)
- history-inputs.test.mjs:158-212 (node): select/text transitions, but with an explicit focusout dispatched between them.
- history.test.mjs:124-160 (node): pending transaction, then a checkpoint.
- Missing: pointerdown before focusout, where `begin()` returns the already-active transaction (history.js:368).

**SE-04** (shortcuts don't fire when a select is focused)
- No test. Nothing in tests/web imports history-shortcuts.js or presses Ctrl/Meta+Z/Y.

**GE-06** (Ctrl+Z during Generate)
- session-operation-consistency.test.mjs:28-33 (node) asserts that `sessionOperationAvailability()` returns null while `processing` is true, so it asserts undo is allowed mid-Generate.
- session-operation-consistency.playwright.spec.js:183-249 (full): only Save/Load are rejected during a gated Generate.
- history-generated-authority.playwright.spec.js:69 undoes after Generate completes. generation-feedback.playwright.spec.js:247 does no Undo mid-run.
- No keyboard Undo during Generate.

**FE-04** (product/protein-ID hide reappears after Redo/unrelated Undo)
- mode-transition-editor-state.playwright.spec.js:141-169 (full): "Exact product" hide, mode round trip, Generate. Asserts the featureVisibility/visibilityRules maps only; no DOM check, no Undo/Redo.
- source-visibility-reconciliation.playwright.spec.js:36-75 (full): Undo/Redo, but only with a feature-scoped override and asserting the map. :76 has a manual product rule but no Undo.
- feature-visibility-actions.test.mjs:81-110, 117-133 (node): feature-scope command apply/revert; the protein_id scope is only listed in the dialog.
- feature-visibility-reconciliation.test.mjs:5-44 (node) looks at overrides only, not rules.

### Test hook mechanism
- gbdraw/web/js/services/runtime-test-hooks.js reads `globalThis.__GBDRAW_TEST_HOOKS__`. It exports:
  - `runtimeTestHooksEnabled()`.
  - `recordStructuralMetric(name, value=1, detail)`, which calls `onStructuralMetric({name, value, ...detail})`.
  - `recordSessionLifecycleEvent(name, detail)`, which calls `onSessionLifecycleEvent({name, timestamp, ...detail})`.
- Other hooks read the global directly:
  - diagram-generation.js:393 `beforeDiagramGenerationResponse` (an async gate or throw, used by failed-generate-bindings, history-generated-authority:317 and right-drawer:847).
  - services/history.js:98 `onHistoryDiagnostic`.
- Emitters: app-setup, label-override-table:124, svg-actions, feature-search/preview-actions, decoration-continuity, preview-runtime:529 (`previewBinderInvocationCount`), run-analysis:1009 (`canonicalCandidateExecutionCount`), config, diagram-generation, diagram-resource-staging, feature-catalog, history-snapshot:996-998, history, session-import-client, session-request, session-resource-backing, session-resources, svg-result-ingestion, ui.
- Metrics are counts only (parse, serialize, scan, worker), sometimes with phase/resultIndex. No hook exposes content or visual equivalence.
- Tests also use `window.__GBDRAW_APP__` and `window.__GBDRAW_HISTORY__` (app-setup.js:1214).

### session-regeneration-contract helper
- helpers/session-regeneration-contract.cjs:27-229 `evaluateSessionRegenerationContract({profile, metrics, events, resultCount, state})` supports six profiles: saved-preview, first-generate, second-generate, failure-before-activation, failure-after-activation, stale-completion.
- It checks exact or maximum structural-metric counts (worker, renderer, sanitize, parse, candidate, activation, rollback, binder, ready receipts), about 20 lifecycle-event orderings, and a few state booleans (`selectedResultCount`, `currentCatalog`, `interactiveProbePassed`, `lastSuccessfulOwnerUnchanged`, `restoredPreviewReady`, `newerArtifactUnchanged`).
- It performs no SVG, content or visual equivalence check.
- session-regeneration-contract.test.mjs (node) only tests the evaluator with synthetic inputs. Its only real consumer is hepatoplasmataceae-session-regeneration.performance.playwright.spec.js:523, 608-633 (perf).
- Existing equivalence tools you could reuse:
  - helpers/svg-visual-semantics.mjs `compareVisualSemantics`, which visual-state.cjs:158-172 only uses for selected/mounted/exported copies.
  - tests/utils/svg_compare.py, called from session-regenerate-intent.playwright.spec.js:1110-1135, which only compares Generate vs Generate and explicitly expects the live draft to differ from Generate (:1986).

### Patterns across the suite
1. **Single-Result fixtures.** helpers/mode-transition.cjs:5-8 seeds (HmmtDNA, lambda) and single-Result unit fixtures back almost every edit and history test. No CI test does batch + edit + switch.
2. **Self-comparing oracles.** Checks compare selected vs mounted vs exported copies of one artifact (visual-state.cjs, multipart-label :77-95, mode-transition-result :16-27, active-result-edit :164-178). None compares a live edit against a fresh Generate.
3. **Tests that assert the current, buggy behaviour:**
   - composition-layout.test.mjs:796 (PV-01).
   - session-operation-consistency.test.mjs:28-33 (GE-06).
   - orthogroup-computation-cache.test.mjs:43-45 (CO-08).
4. **State, not render.** Override maps, deltas and `toContain(color)` stand in for rendered order, containment or survival after Generate (color-captions :96/:111, mode-transition-editor-state).
5. **Happy paths only.** Valid stroke value 2, direct input clicks, explicit focusout ordering, no select-focused shortcuts, no record selector/region/reverse complement combined with a live edit.
6. **Mocked collaborators in history unit tests.** Faked `applyEditorStateData` and `applyConfigData` hide the real catalog-admission and reconcile paths. No Playwright test chains a checkpoint Undo/Redo into a mode switch.
7. **Tiering.** Nearly all nearest Playwright tests run only post-merge (full). The .playwright.py scripts and retained journeys are not wired into CI.
