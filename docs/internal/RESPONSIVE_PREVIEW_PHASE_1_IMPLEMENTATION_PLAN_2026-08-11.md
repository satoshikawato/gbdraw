# Work package E: Responsive Preview Phase 1 implementation plan

- Date: 2026-08-11
- Status: planned
- Baseline: `docs_renovation` / `6f14e2c4fd2a`
- Source: [`gbdraw v0.14.0 Release Roadmap and Codex Implementation Brief`](gbdraw_v0.14.0_codex_roadmap.md), Work package E
- Prerequisites: Work packages C and D; Web performance remediation Phases 1 and 2

Related plans and contracts:

- [Work package C implementation plan](WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package D implementation plan](WORK_PACKAGE_D_MANUAL_FEATURE_PLACEMENT_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Web application performance remediation plan](WEBAPP_PERFORMANCE_REMEDIATION_IMPLEMENTATION_PLAN_2026-08-10.md)
- [Web application maintenance contract](../../gbdraw/web/CLAUDE.md)
- [Session and request compatibility history](../SESSION_COMPATIBILITY.md)

This document turns Work package E into executable phases. At creation time it changes documentation only. Production code, tests, Gallery assets, browser wheels, and reference SVGs have not been changed by this plan.

The planning audit began at `a2c69de112ae`. The shared worktree advanced through the Work package C and D planning commits while this document was being written; the final plan was rechecked against `6f14e2c4fd2a` and its current C/D roadmap contracts.

## 1. Outcome

Responsive Preview Phase 1 keeps the current direct preview edits that already have stable semantic targets and adds automatic renderer updates for a narrow geometry scope:

- per-record display start from Work package C, for Circular 12 o'clock and Linear left-edge anchoring;
- manual feature placement;
- `feature_overlap_tolerance_bp` if Work package D retains it;
- existing label-layout reconciliation after direct label edits.

The automatic path is a renderer path, not a shortened Generate path. It starts from the last committed canonical request/resources, overlays an allowlisted render-safe draft, calls the existing typed renderer, stages the complete candidate, and commits only the current revision. LOSAT, input discovery, comparison resolution, and analysis-cache mutation are unreachable from this path.

Phase 1 does not promise automatic rendering for every figure setting. Additional fields join the allowlist only after they satisfy the same request, persistence, history, race, and performance gates.

## 2. Current architecture and constraints

The implementation must account for these existing owners.

| Area | Current owner | Constraint for Work package E |
| --- | --- | --- |
| Full Generate | `gbdraw/web/js/app/run-analysis.js` | `runAnalysisInternal()` validates inputs, may run LOSAT, builds the canonical request, renders, stages, and commits. It remains the explicit analysis path. |
| Existing label reflow | `runLabelReflow()` in `run-analysis.js` | It calls `runAnalysisInternal({ runMode: 'reflow' })`. A cache miss can reach LOSAT, and reflow replaces Results without adopting a full catalogue/canonical bundle. It cannot become the general preview path. |
| Diagram worker | `gbdraw/web/js/services/diagram-generation.js` | One render may be active. Concurrent calls are rejected. Cancellation terminates the shared worker and loses the warm Pyodide runtime. |
| Candidate staging | `gbdraw/web/js/app/candidate-render.js` | `prepareCandidateRenderCommit()` already validates catalogue-derived feature state, sanitizes candidate SVG, reapplies stable overrides, and serializes before live mutation. Reuse this boundary. |
| Committed canonical authority | `gbdraw/web/js/services/config.js` | `committedCanonicalSession` contains the last successful request/resources but is module-private. Expose an immutable render basis; do not deep-clone unchanged large resource payloads on every update. |
| Direct preview editing | `gbdraw/web/js/app/preview-runtime.js` and feature/label/legend/layout action modules | Fill, visibility, and stroke are first-class runtime patches. Other mutations are distributed. Do not add another DOM/serialization authority. |
| Session projection | `gbdraw/web/js/services/session-request.js` | This is the only Web-to-canonical request boundary. Render-safe overlays belong here, not in `config.js` or an argv-shaped builder. |
| History | `gbdraw/web/js/services/history.js`, `history-snapshot.js`, and `app/history-inputs.js` | Current snapshots include large generated artifacts. Automatic rendering cannot start until the performance plan separates user-intent history from artifact snapshots. |
| Preview preference | `autoLabelReflowEnabled` and `paletteInstantPreviewEnabled` | Two controls currently describe overlapping behaviour. Phase 1 replaces their visible ownership with one automatic-layout preference. |
| Draft/committed mismatch | scattered generated-mode comparisons and watchers | There is no generic draft revision, render-pending fingerprint, or analysis-stale state. Some draft mode changes clear editor state while the old Result remains visible. |

The measured large-session costs are recorded in [`WEBAPP_PERFORMANCE_REMEDIATION_IMPLEMENTATION_PLAN_2026-08-10.md`](WEBAPP_PERFORMANCE_REMEDIATION_IMPLEMENTATION_PLAN_2026-08-10.md). In particular, history begin/commit and repeated SVG sanitization are already expensive enough that automatic updates would amplify them.

## 3. Fixed product decisions

### D1. Effects are composable

Each admitted edit has an explicit policy with three independent booleans:

```text
patchNow
needsRender
invalidatesAnalysis
```

The registry is exhaustive for admitted fields and fail-closed for unknown fields. An unknown field does not trigger automatic work. It remains pending for explicit Generate.

Track a normalized full drawing-intent fingerprint in addition to the render-safe and analysis fingerprints. Every drawing-control action reports its normalized semantic key before policy lookup. A `manual-only` or unknown key therefore creates visible pending divergence even though it dispatches no automatic work. An accepted allowlisted render advances only the committed render-safe projection; it cannot clear unrelated manual-only or unknown divergence. Explicit Generate is the only operation that adopts the complete draft fingerprint.

Do not infer the policy from a deep watch of all reactive state. Call the coordinator from the action that owns the semantic change, or watch a small normalized projection where a control has no action owner.

### D2. Initial effect register

| Edit | `patchNow` | `needsRender` | `invalidatesAnalysis` | Automatic admission | Phase 1 rule |
| --- | ---: | ---: | ---: | --- | --- |
| Feature/palette fill | yes | no | no | admitted | Patch existing semantic feature targets and update the established colour-rule/editor owner. |
| Existing feature stroke colour | yes | no | no | admitted | Keep the current editor overlay and Result serialization path. Do not invent a new public stroke setting. |
| Existing feature stroke width | yes | no | no | direct-only | Patch immediately. `needsRender` is always false in Phase 1; a later phase may reclassify it only with bounds/clipping evidence and a reproducible owner. |
| Visibility: hide an existing target | yes | yes | no | admitted | Hide immediately, then reconcile labels/legend/layout through render-only work when admitted by the current action. |
| Visibility: show a target absent from the DOM | no | yes | no | admitted | Renderer fallback is required because the omitted element cannot be patched into existence. |
| Label text/visibility on an existing label | yes | yes | no | admitted | Patch optimistically and use the no-analysis renderer for label geometry. |
| Existing legend/composition/canvas movement | yes | no | no | direct-only | Keep only current operations with an established session/export owner. |
| Per-record display start | no | yes | no | admitted | Release-blocking render-only field from Work package C in both Circular and Linear modes. |
| Manual placement | no | yes | no | admitted | Release-blocking render-only field from Work package D. |
| Feature-overlap tolerance | no | yes | no | conditional | Include only if retained by Work package D. |
| Annotation colour/label or legend create/delete/rename | no | no | no | manual-only | Explicit Generate in Phase 1. Add a direct patcher only after semantic metadata/catalogue commit is designed. |
| Font, arrow, track, annotation-lane, scale/tick geometry | no | no | no | manual-only | Explicit Generate in Phase 1. These are later allowlist candidates. |
| Input sequence, record endpoint, LOSAT program/task, search genetic code, raw-evidence arguments | no | no | yes | stale-only | Mark analysis stale and require explicit Generate/LOSAT. |

`direct-only`, `conditional`, `manual-only`, and `stale-only` are explicit admission states, not hidden fallback branches. Tests must name the state or prove that no automatic render occurs.

### D3. Committed Result and draft intent stay separate

The committed bundle contains:

- Result list and selected Result;
- canonical request/resources and their Web resource bindings;
- feature catalogue and extracted/biological feature state;
- render metadata and output topology;
- editor bindings and stable selection identity;
- export and session-save authority.

Draft controls may be newer. Transient intent, render, direct-edit, analysis, and lifecycle revisions; the three fingerprints; scheduler tokens; status; and viewport snapshots are not part of the canonical request.

A successful direct patch advances only its semantic keys in the committed intent baseline. A mixed patch/render edit can therefore have matching text/colour intent while its render-safe geometry remains pending. A render-only commit advances only its allowlisted render-safe keys. This field-level adoption prevents either path from clearing unrelated divergence.

### D4. Analysis stale pauses automatic geometry rendering

When an analysis-owned field changes, set `analysis-stale` and leave the committed Result visible. Proven-safe direct patches may still affect that Result. Do not overlay new geometry settings onto old analysis evidence in Phase 1. Show `Generate and update` when renderer input is stale without a LOSAT run, and `Run LOSAT and update` when the active comparison workflow requires new LOSAT evidence. Only that successful explicit action clears the state.

### D5. One visible automatic-update preference

The UI label is `Auto-update layout`. It controls renderer-based geometry reconciliation. Proven-safe direct patches remain immediate when the preference is off.

Use `autoUpdatePreviewEnabled` as the runtime owner and `Auto-update layout` as the visible label, but keep `ui.autoLabelReflow` as the sole persisted wire key in Phase 1. Its name is a documented compatibility alias whose meaning expands to the consolidated automatic-layout policy. Do not also write `ui.autoUpdatePreview`; rename the wire key only in a future evidence-backed session format. Read the compatibility key from every supported historical session version that stored it. The inventory must cover at least versions 29–33 and 39–40 and verify whether versions 27–28 carried it.

The old `ui.paletteInstantPreviewEnabled` field no longer controls behaviour. Migrate applied/pending palette state first, then omit the flag from runtime/history/current-session output. Every supported historical session with saved applied/pending palette state must load without changing the committed Result during hydration.

Represent the runtime preference with provenance, not as an ambiguous false boolean:

- `unset` for a new/Reset state or a restored session with no supported preference field;
- `set` after the first successful manual Generate defaults it on, the user changes the control, or a supported session restores either on or off.

Before any Result exists, `unset` performs no automatic render. The first successful manual Generate changes only `unset` to `set`/on. It must not overwrite an explicit or restored off choice, and a failed/cancelled Generate must not change it.

Serialize `ui.autoLabelReflow` only when provenance is set; write the Boolean value for either on or off. Omit the key while provenance is `unset`, including a field-absent historical session saved again without a successful Generate. Reset restores `unset`. The automatic-update preference is session UI metadata, not figure-edit history: toggling it creates no Undo item, and history snapshots omit the preference entirely.

This UI cleanup does not justify a request-schema change. At Phase 0, audit the C/D writer status. If session 41 is not main- or release-backed, carry the compatibility key and palette normalization through that one coordinated writer. If it is already main- or release-backed and changing the current-session UI inventory requires a format change, allocate the next session version from compatibility evidence while keeping canonical request schema 6. Never silently redefine a released session 41 contract.

### D6. Scheduler semantics are fixed

- trailing-edge debounce: 400 ms;
- capacity: one active renderer request and one replaceable latest pending draft;
- ordinary supersession: discard the obsolete response; do not terminate the worker;
- manual Generate: invalidate queued preview work and take the next renderer slot;
- render/lifecycle and replay-revision checks: immediately before staging, after staging, and before commit;
- session import, Reset, mode change, committed Result replacement, and disposal: invalidate queued revisions through the same coordinator.

Full analysis may prepare its inputs while an obsolete preview render finishes, but every diagram-worker job enters one lifecycle broker. The broker serializes renderer and feature-extraction jobs because explicit cancellation terminates the shared worker. Manual-workflow dependencies and the manual renderer retain priority, interactive feature extraction remains queued by stable job identity, and automatic rendering uses the one replaceable latest slot. A preview request must never cause the existing concurrent-request error.

Use separate revision domains. A direct-only patch advances the intent/direct-edit revisions but not the render revision. Candidate preparation snapshots and replays the latest direct-edit state. If that revision changes during staging, discard the candidate and automatically requeue the still-pending latest render-safe draft; never commit a stale replay and never lose the geometry request. Render-safe, analysis, and lifecycle changes keep their own invalidation rules.

`Update layout` is a manual invocation of the same no-analysis executor, not Generate. It bypasses debounce, dispatches the latest render-safe draft exactly once, and does not change an off or slow-paused preference. Duplicate clicks for the same fingerprint coalesce. If another worker job is active, keep one latest pending layout request. Analysis-stale state dispatches nothing. Failure retains the prior bundle and pending draft so one later click can retry. `Resume auto-update` clears only the slow pause and schedules pending work through the normal 400 ms rule. Manual layout updates do not contribute slow-guard samples.

### D7. Automatic commit preserves context

All fallible parsing, catalogue validation, sanitization, direct-edit replay, identity reconciliation, and canonical adoption preparation happen before the first live assignment. The commit then performs no fallible work.

Candidate replay is exhaustive for the admitted direct-edit register. Palette/default and per-feature fill, feature stroke, visibility, label text/visibility, and supported legend/composition/canvas transforms are either represented in the candidate request or reprojected from their established semantic editor owner. Replay uses stable identity and a declared geometry-compatibility rule; it never copies a stale DOM subtree. Missing or incompatible required replay data fails candidate preparation and leaves the prior committed bundle unchanged.

Automatic commit preserves the selected Result index, zoom, pan, and stable feature/annotation selection where the candidate catalogue still contains the identity. Feature selection uses its stable biological identity; annotation selection uses the composite record, track, set, and annotation ID. Rebind the mounted SVG and popup after commit, and clear only missing identities. Manual Generate may keep its current reset policy.

### D8. Renderer completion is history-neutral

One continuous user interaction creates one history item. The accepted renderer result is a derived artifact update and creates no history item. After the performance-plan history separation, history may refresh its current artifact reference without pushing an entry.

Undo/redo restores user intent, advances the revision domains required by the restored effects, and requests at most one render. Auto-update off leaves `Changes pending`; it does not render during Undo/Redo.

### D9. Slow-update guard is deterministic

Do not count the first render after worker initialization or session hydration. Pause automatic layout updates for the current committed artifact after two consecutive accepted updates each take more than 2,000 ms across worker rendering plus candidate staging/commit. Failed, stale, and superseded requests do not contribute a sample.

Pause state is transient. Direct patches remain active. The UI offers `Update layout` and `Resume auto-update`.

### D10. Analytics counts accepted outcomes

If analytics is enabled, emit an automatic `render_result` only after the current candidate commits. Debounce replacement, stale/superseded completion, cancellation, and direct patches emit no render outcome. Use bounded trigger values only; do not include field values, biological identifiers, labels, coordinates, or raw errors.

## 4. Surface and authority matrix

| Surface | Change in Work package E | Evidence required |
| --- | --- | --- |
| Typed core/request | No new preview field. C/D semantic fields must already exist in the canonical request. | Existing C/D request decode and non-default consumer tests pass. |
| Python API | No preview scheduler API. Display-start/placement parity remains owned by C/D. | No public signature drift from E. |
| CLI | No preview option. Display-start/placement parity remains owned by C/D. | No help/default drift from E. |
| Web state | Add one preview preference and transient coordinator state. | Defaults, Reset, session hydration, and no-auto-render-on-load tests. |
| Canonical Web projection | Add a strict copy-on-write function that overlays only registered render-safe fields onto an immutable committed basis. | Unknown paths reject; committed inputs/evidence remain unchanged; unchanged large payloads are not deep-cloned before worker dispatch. |
| Persistence | Keep `ui.autoLabelReflow` as the single compatibility wire key for the renamed preference; remove the Palette toggle only after migrating its applied/pending state; normalize Web and Python/CLI current-session writers. | Current-writer round-trip, migration fixtures for every supported historical preference/palette shape, no retired key in current output/history, and draft/committed mismatch tests. |
| Rendering | Reuse the typed diagram worker and candidate sanitizer/commit path. | Circular/Linear and single/grid/batch topology tests. |
| History | Store user intent, not a new artifact snapshot for renderer completion. | One continuous edit is one entry; Undo/Redo schedules once. |
| Export | Continue using the committed Result; direct patches flush through their established owner. | Clean/interactive SVG and raster export reflect accepted state only. |
| Documentation/Gallery | Update public material only when the final control is part of a documented workflow. | Regenerated evidence, not hand-edited Gallery artifacts. |

## 5. Target module boundaries

### `app/responsive-preview.js`

Keep the `createResponsivePreview*` entry point at the top level. It owns:

- effect-policy dispatch;
- intent, render, direct-edit, analysis, and lifecycle revisions plus normalized fingerprints;
- debounce and latest-wins scheduling;
- automatic/manual renderer arbitration;
- status and slow-pause state;
- invalidation on import, Reset, mode/result replacement, and disposal;
- manual `Update layout` and resume actions.

Split pure policy or scheduler helpers under `app/responsive-preview/` only if the top-level module would otherwise mix distinct responsibilities. Do not create a general state framework.

### `services/session-request.js`

Add the strict render-only projection at the existing canonical boundary. It receives an opaque immutable committed basis plus a normalized allowlisted delta. It clones mutable request and resource-index structures copy-on-write, while reusing immutable unchanged resource payloads until the existing worker dispatch performs its required structured clone/transfer. It never reads mutable files, caches, or reactive comparison state.

### `services/config.js`

Expose a read-only accessor for an immutable committed render basis. Keep mutable module-local authority private, and do not return objects a caller can mutate. Adopt a new canonical bundle only after candidate staging succeeds and the revision remains current.

### Shared candidate commit

Factor the manual Generate candidate assignment into a shared owner only where it removes duplicated commit logic between manual and automatic rendering. The shared function accepts an interaction policy:

- manual: existing selection/viewport reset behaviour;
- automatic: stable-identity selection and viewport preservation.

Both use the same catalogue validation, sanitizer, exhaustive admitted direct-edit replay, metadata projection, and rollback contract. Replay covers palette/default and per-feature fill, stroke, visibility, label text/visibility, and supported legend/composition/canvas transforms through their semantic owners. Remove the superseded Results-only reflow commit.

### `services/diagram-generation.js`

Keep one worker client and the existing typed request payloads. Add one lifecycle broker for renderer and feature-extraction jobs so only one worker job is active and queued jobs survive an explicitly cancelled worker generation. Do not add a second worker, a second Pyodide runtime, or an argv request.

## 6. Execution phases

### Phase 0: prerequisite and baseline gate

Status: pending

#### Required state

1. Work package C has frozen final per-record display-start semantics and record-aligned provenance. Before Phase 1 production work, its canonical schema 6 fields and Circular/Linear tests must have landed.
2. Work package D has frozen final manual-placement identity and request representation. Before Phase 1 production work, its schema 6 fields and tests must have landed. Record whether overlap tolerance remains in release scope.
3. Performance-remediation Phase 1 has separated normal setting history from generated artifact snapshots and recorded passing evidence.
4. Performance-remediation Phase 2 has established one shared SVG ingestion/sanitize/commit boundary and recorded passing evidence.
5. Audit whether the coordinated C/D session writer is still pending/branch-local or already main-/release-backed. Record whether UI normalization can remain in session 41 or needs the next evidence-backed session version; do not change request schema 6 for this UI metadata.

Do not duplicate those prerequisite implementations in this plan. Execute their owning plans first when they remain pending.

#### Baseline work

- Record branch, HEAD, dirty files, and pre-existing failures.
- Confirm Node and Python Playwright availability.
- Capture current counts for label reflow: LOSAT executor calls, cache probes, input reads, worker renders, sanitize passes, Result serializations, history begin/commit time, and viewport reset.
- Measure cold and warm renderer-only candidates on a small diagram and the WSSV large-session fixture.
- Confirm the current Palette instant preview and Auto Reflow runtime, history, Web writer, and Python/CLI pass-through fields and visible controls.

#### Commands

```bash
git status --short --untracked-files=all
git diff --stat
git rev-parse --short=12 HEAD
command -v playwright && playwright --version
node -e "console.log(require.resolve('@playwright/test'))"
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
python tools/prepare_browser_wheel.py
node --test tests/web/preview-runtime.test.mjs tests/web/history.test.mjs tests/web/history-inputs.test.mjs
node --test tests/web/run-analysis-simple-path.test.mjs tests/web/session-request.test.mjs tests/web/session-authority.test.mjs
node --test tests/web/svg-sanitization.test.mjs tests/web/diagram-generation-worker.test.mjs
```

#### Completion gate

- All five prerequisites have named evidence, including the current-session writer decision.
- Counters distinguish analysis, worker rendering, candidate staging, serialization, and history costs.
- Any pre-existing failure is recorded before production changes.

#### Stop condition

Stop this plan if C/D semantic fields are not stable or either performance P0 boundary is absent. Continue through the owning prerequisite plan rather than adding an E-only workaround.

### Phase 1: effect policy, revisions, and state authority

Status: pending

#### Owners

- new `gbdraw/web/js/app/responsive-preview.js`
- optional pure helpers under `gbdraw/web/js/app/responsive-preview/`
- `gbdraw/web/js/state.js`
- `gbdraw/web/js/services/session-request.js`
- new `tests/web/responsive-preview-policy.test.mjs`

#### Work

1. Encode the D2 register as a pure policy. Unknown field/action keys return an explicit fail-closed result.
2. Define normalized full drawing-intent, render-safe, and analysis fingerprints. Use existing stable input/resource identities or mutation revisions; do not reread or rehash uploaded file contents on every control event. Exclude transient UI state, cache metadata, Results, and viewport.
3. Add monotonic intent, render, direct-edit, analysis, and lifecycle revisions plus transient status without persisting them. Define exactly which effects advance each domain.
4. Make every drawing-control action report a normalized semantic key before effect lookup. Add explicit actions for render-safe, analysis-invalidating, manual-only, unknown, and mixed patch/render changes. Use a tested control/action inventory instead of a deep watch over the entire state object.
5. Derive `up-to-date`, `changes-pending`, and `analysis-stale` from all three fingerprints and the committed bundle rather than manually clearing strings in many callers. Use priority `analysis-stale` > manual/render pending > up to date.

#### Tests

- every initial allowlist entry returns the expected effect combination;
- label text combines patch and render;
- visibility is target-presence aware;
- analysis-owned fields never request automatic rendering;
- manual-only and unknown fields fail closed, set `Changes pending`, and dispatch no automatic render;
- an admitted render advances only the render-safe committed projection and does not clear unrelated manual-only/unknown divergence;
- explicit Generate adopts the full drawing-intent fingerprint and clears pending state;
- fingerprint order is deterministic and ignores transient state;
- fingerprint evaluation performs no file read, resource serialization, or content hash;
- Undo/Redo restoration advances the required intent/effect revisions without a history entry from the coordinator;
- a direct-only action advances the direct-edit revision without erasing or independently dispatching an existing render-safe divergence.

#### Completion gate

The policy is pure, exhaustive for Phase 1, and has one owner. No worker, file, cache, DOM, or session API is called by the policy tests.

### Phase 2: shared candidate transaction

Status: pending

#### Owners

- `gbdraw/web/js/app/candidate-render.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/legend-layout/canvas-actions.js`
- `gbdraw/web/js/app/legend-layout/composition-actions.js`
- a focused shared commit module under `gbdraw/web/js/app/` only if needed by both manual and automatic paths
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/app/watchers.js`
- `tests/web/run-analysis-simple-path.test.mjs`
- `tests/web/svg-sanitization.test.mjs`

#### Work

1. Separate candidate preparation from live assignment in the manual path.
2. Stage Results, validated feature catalogue, extracted/biological feature state, selector safety data, orthogroup/collinear metadata, track geometry, generation metadata, stable selection reconciliation, canonical artifacts, and export/run metadata.
3. Perform all fallible work before live assignment. Preserve the existing manual rollback on late failure.
4. Add an automatic interaction policy that preserves selected Result, stable selection, zoom, and pan.
5. Make draft mode watchers retain committed catalogue/editor state while the old Result remains authoritative.
6. Keep manual Generate behaviour and default rendering bytes unchanged.
7. Define and test the replay contract for every admitted direct-edit family. Reject candidate preparation if a required semantic replay cannot be applied safely; do not silently lose the edit.
8. Move base-config capture and composition-metadata parsing from the post-`svgContent` watcher into candidate preflight over the staged sanitized SVG. The mounted watcher consumes the prevalidated projection. Any unavoidable post-mount binding failure must roll back the entire logical transaction before success status, preference transition, or analytics.
9. Snapshot the direct-edit revision with replay. If it changes before commit, reject the candidate without live mutation and return the still-pending render-safe draft to the scheduler.

#### Tests

- manual Generate retains existing successful and rollback behaviour;
- candidate sanitizer/catalogue/adoption failure changes no committed owner;
- automatic policy restores only identities present in the candidate;
- annotation selection with the same local ID in different record/track/set scopes restores only the composite match;
- automatic policy preserves viewport; manual policy keeps its declared reset behaviour;
- draft mode/input changes do not clear committed catalogue/editor state;
- one accepted Result crosses the shared sanitizer boundary at most once, as required by the performance plan;
- palette/default and per-feature fill, stroke, visibility, label text/visibility, and supported legend/composition/canvas transforms survive a geometry render;
- missing or incompatible required replay data rolls back the entire candidate;
- injected base-config/composition parsing failure mutates no live owner, and an injected post-mount binding failure restores the complete prior bundle.
- fill/stroke/visibility/label or supported composition edits immediately before and during staging cannot be overwritten; a replay-revision mismatch commits nothing and requests one latest replacement render.

#### Completion gate

Manual and future automatic rendering can call one candidate transaction without a Results-only branch. Existing reference outputs are unchanged.

#### Stop condition

Stop if the factor changes manual output topology, reference SVG bytes, session authority, or cancellation rollback. Fix the common boundary before adding automatic scheduling.

### Phase 3: dedicated no-analysis render-only execution

Status: pending

#### Owners

- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/app/responsive-preview.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/diagram-generation.js`
- `gbdraw/web/js/app/feature-metadata-extraction.js`
- `tests/web/responsive-preview.test.mjs`
- `tests/web/session-request.test.mjs`
- `tests/web/linear-multi-record.playwright.spec.js`

#### Work

1. Add an opaque read-only accessor for the committed canonical render basis. Render-only execution needs request/resources, not a deep clone of inactive `webFiles` bindings.
2. Add a strict copy-on-write render-only overlay function at `session-request.js`. It accepts only the Phase 1 allowlist and rejects unknown paths, mode mismatches, missing stable identities, invalid values, and topology changes.
3. Preserve committed input resources and resolved comparison evidence. Label/resource overlays may update only generated editor tables already owned by the canonical request.
4. Call `runDiagramGeneration()` directly with the complete candidate request/resources. Do not call `runAnalysisInternal()` or any LOSAT/helper preparation function.
5. Require the worker response to include a complete feature catalogue. Stage through Phase 2 and adopt canonical artifacts only while the render/lifecycle revisions and direct-edit replay revision remain current.
6. Move label-layout reconciliation onto this entry point. Remove the superseded `runMode: 'reflow'` path after all callers migrate.
7. After a successful candidate commit, update only the admitted render-safe baseline. Preserve any full-intent divergence caused by manual-only or unknown fields.

#### Hard no-analysis tests

Run with empty caches and spies that throw on:

- Circular and Linear LOSAT executors;
- cache lookup/probe/promotion;
- mutable input serialization or file reads;
- comparison-plan resolution;
- main-thread helper Pyodide initialization;
- analysis-resource mutation.

Per-record display start, placement, retained tolerance, and label reconciliation must still render successfully. The candidate request must reuse committed input/evidence resources and change only allowlisted render fields or generated editor tables.

Replace the existing Linear browser assertion that expects label reflow to wait on `File.arrayBuffer()`. Block the mutable file read without releasing it and prove that render-only dispatch and commit complete from committed resources while the read count stays zero.

Add a mixed-divergence fixture: change a manual-only font/track field, then an admitted display-start field. The automatic candidate uses the committed font/track value, applies display start, and leaves `Changes pending` until explicit Generate commits the complete draft.

Add a large-resource fixture that proves projection does not deep-clone unchanged resource payloads before worker dispatch and cannot mutate the committed basis through a returned reference. Extend the same counter through accepted canonical adoption and current-history artifact refresh: one automatic update must not deep-clone the complete unchanged WSSV resource bundle at any point in the transaction.

#### Completion gate

The module import/call graph and runtime spies prove that render-only execution cannot reach analysis. A cache hit is not used as proof.

#### Stop condition

Do not add an `isRenderOnly` boolean to the full analysis monolith. If committed artifacts cannot reproduce a renderer request without inspecting mutable draft inputs, fix committed authority first.

### Phase 4: single-flight scheduler and manual arbitration

Status: pending

#### Owners

- `gbdraw/web/js/app/responsive-preview.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/services/diagram-generation.js`
- `gbdraw/web/js/app/feature-metadata-extraction.js`
- `tests/web/responsive-preview.test.mjs`
- `tests/web/diagram-generation-worker.test.mjs`

#### Work

1. Implement 400 ms trailing debounce with injectable clock/timers for tests.
2. Keep one active render and replace any pending draft with the latest normalized snapshot.
3. Mark an active obsolete request stale without terminating the worker. After it settles, run the latest pending request.
4. Route the final manual Generate renderer dispatch through the same coordinator. Manual work invalidates queued automatic work and owns the next slot.
5. On successful manual Generate, adopt the complete current draft fingerprint and clear pending automatic/manual work. On manual failure, retain the prior committed Result and all divergence states.
6. Invalidate revisions during session import, Reset, mode change, committed Result replacement, and disposal.
7. Put renderer and feature-extraction jobs behind one worker lifecycle broker. Keep at most one worker job active. Preserve queued jobs by stable identity across worker generations; manual-workflow dependencies/manual render run before automatic preview, and obsolete automatic work keeps latest-wins semantics.
8. Route the existing Generate Cancel action through the coordinator. Cancelling a queued manual workflow removes only its owned queued jobs. Cancelling its active feature-extraction or renderer job settles the UI and releases the slot without waiting for Python, terminates that worker generation, and retains the prior committed bundle. Because no unrelated worker job is active, queued unrelated feature extraction survives and runs after worker recreation. Analysis-specific cancellation still stops the LOSAT/comparison job it owns. Reset/import and ordinary supersession mark work stale without termination; disposal and fatal recovery may terminate.
9. Check job identity/cancellation before staging, after staging, and before commit. Cancellation during candidate preparation or mount rollback cannot publish status, change preference provenance, or emit analytics.
10. Route `Update layout` through the no-analysis coordinator with the D6 immediate/deduplicated semantics. Keep it distinct from both full Generate and preference resume.

#### Tests with fake timers/deferred promises

- one change waits 399 ms and runs at 400 ms;
- repeated changes inside the window produce one render with the latest snapshot;
- edits during an active render leave one latest pending render;
- stale completion never reaches candidate commit;
- automatic updates never cause the concurrent-request error;
- manual Generate drops queued work and gets the next slot;
- the first successful manual Generate turns an unset preference on, while a failed/cancelled Generate and an explicit or restored off choice remain unchanged;
- session import/Reset/disposal invalidates queued and active revisions;
- Cancel before manual dispatch removes only that workflow's queued jobs without terminating an unrelated active worker job;
- Cancel during an owned active feature-extraction or diagram-render job settles the manual UI and releases the broker slot within 250 ms in the injected-clock browser fixture, terminates the worker generation, and commits nothing;
- Cancel after the worker response or during candidate staging/commit preserves the prior complete bundle and a subsequent Generate/render-only interaction succeeds;
- an unrelated feature-extraction job queued behind the cancelled job survives worker recreation, settles once, and populates its cache; no Promise is left pending or double-settled;
- direct calls to the old global cancellation path are absent outside the coordinator/worker lifecycle boundary;
- ordinary supersession never terminates the worker; explicit Cancel terminates only when its owned job is active under the one-job broker invariant.
- auto off edit dispatches zero renders, one `Update layout` click dispatches exactly one immediate render, and a duplicate click or click during an active job leaves at most one latest request;
- slow-paused `Update layout` performs one render without resuming or adding a slow sample; `Resume auto-update` uses normal debounce;
- analysis stale dispatches zero render-only work, and a failed manual layout update remains pending for a one-click retry.

#### Completion gate

Scheduling behaviour is deterministic without real sleeps. All renderer and feature-extraction jobs pass through one shared-worker lifecycle broker, while one coordinator owns renderer policy.

### Phase 5: direct-patch convergence and session preparation

Status: pending

#### Owners

- `gbdraw/web/js/app/preview-runtime.js`
- existing feature colour/visibility/stroke and label action modules
- existing legend/composition/layout action modules in the admitted scope
- `gbdraw/web/js/app/results.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/session.py`
- `gbdraw/session_io.py`
- `gbdraw/api/session_compat.py`
- `gbdraw/cli_utils/session.py`
- `tests/web/preview-runtime.test.mjs`
- relevant session authority/draft tests
- `tests/test_session_io.py`
- `tests/test_session_compat.py`
- `tests/test_api_session.py`
- `tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py`

#### Work

1. Route admitted immediate mutations through `preview-runtime` or the single ingestion/commit boundary established by the performance plan.
2. Apply each user action to every affected Result or semantic owner needed for Circular batch/inactive outputs. Do not update only the mounted SVG when session/export can select another Result.
3. Keep one serialization/Result update per direct action.
4. Rename the runtime owner to `autoUpdatePreviewEnabled`, but serialize only the documented `ui.autoLabelReflow` compatibility key. Migrate it from every supported historical version that wrote it and preserve explicit on/off provenance.
5. Preserve applied/pending palette state from every supported historical shape and derive draft/committed mismatch during hydration without patching or automatically rendering.
6. After palette migration, remove `paletteInstantPreviewEnabled` from state, Results actions, watchers, config, history snapshots, reset, session-authority inventory, and all current Web-session output.
7. Normalize Python/CLI load→save and typed-adjunct paths so a current session keeps the preference compatibility key but does not pass the retired Palette flag through arbitrary `config` or `ui` cloning. Apply the Phase 0 version decision rather than silently changing a release-backed writer.
8. Do not add annotation or general legend-content patching in this phase.

#### Tests

- fill/stroke/visibility/label actions patch only admitted semantic targets;
- hide/show falls back to rendering when the target is absent;
- label edit is immediate and schedules one reconciliation;
- feature-label discovery excludes annotation text from direct label patching and reconciliation triggers;
- direct patches remain immediate with auto layout off;
- the current writer round-trip plus fixtures for the oldest and newest supported historical preference shapes restore explicit on/off choices through the single compatibility key;
- supported historical applied/pending palette fixtures retain their committed-versus-pending meaning;
- current Web-session and Python/CLI load→save output contain `ui.autoLabelReflow` only for set provenance and contain neither `paletteInstantPreviewEnabled` nor `ui.autoUpdatePreview`;
- history snapshots contain none of `autoLabelReflow`, `autoUpdatePreview`, or `paletteInstantPreviewEnabled`;
- version-40-to-current CLI replay and typed-adjunct tests preserve the compatibility preference plus applied/pending palette state while removing both `config.paletteInstantPreviewEnabled` and `ui.paletteInstantPreviewEnabled`;
- session hydration performs zero automatic renders;
- Circular batch and inactive Results do not diverge from export/session state;
- annotation and unclassified legend operations do not enter direct patching.

#### Completion gate

The admitted direct patches and preference/session semantics have one internal owner, supported historical sessions retain preference and committed/pending palette meaning, hydration performs no automatic render, and no duplicated direct-mutation/sanitization pipeline exists. Do not expose the consolidated UI until Phase 6 passes.

### Phase 6: history, viewport, UI exposure, slow guard, and analytics

Status: pending

#### Owners

- `gbdraw/web/js/app/history-inputs.js`
- `gbdraw/web/js/services/history.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/app/responsive-preview.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/index.html`
- preview navigation/selection owners used by the existing Result mount
- analytics wrapper/event owner if Work package I has landed
- `tools/verify_gui_offline.py`
- `tests/test_web_packaging.py`
- history, session, and responsive-preview tests
- new `tests/web/responsive-preview.playwright.spec.js`

#### Work

1. Keep input focus/change coalescing for continuous text/number edits and add explicit coverage for range/spinner interactions used by C/D.
2. Add the performance-plan mechanism that adopts a derived artifact into current history state without pushing a user-visible entry.
3. Make Undo/Redo schedule one render for restored intent when auto layout is on; otherwise mark pending.
4. Keep the auto-update preference, runtime status, tokens, revisions, viewport snapshots, and slow-pause state out of figure-edit history. Persist only the preference through the D5 optional session compatibility key.
5. Snapshot zoom, pan, selected Result, and stable selection before automatic commit; restore them after the candidate mounts.
6. Implement the D9 slow guard with injected timing evidence.
7. Emit analytics only for accepted current commits and genuine current-candidate failures.
8. Replace visible Palette instant preview and Auto Reflow controls with `Auto-update layout` plus help text explaining that direct edits remain immediate.
9. Add the non-blocking statuses: `Up to date`, `Updating…`, `Changes pending`, `Analysis out of date`, failure with prior Result retained, and slow pause.
10. Show `Update layout` only for a render-safe pending draft while automatic updates are off/paused. Show the active-intent `Generate and update` or `Run LOSAT and update` action for analysis stale.
11. Do not display the full-page Generate overlay for automatic work.
12. Replace the offline smoke test's deferred/instant Palette assertions with the consolidated contract: Palette/direct patches are immediate even when automatic layout is off, a render-safe geometry edit becomes pending or updates according to the preference, and session hydration dispatches no render. Rename the packaging test so it no longer promises the removed Palette toggle behaviour.
13. Extend the offline verifier so the same network-blocked smoke flow can serve `gbdraw/web/` extracted from a built wheel, not only the source `WEB_ROOT`. Keep wheel member inspection, but do not treat filenames alone as runtime proof.

#### Tests

- 8 → 9 → 10 → 11 in one interaction is one history item;
- renderer completion adds no item;
- toggling `Auto-update layout` adds no history item and Undo/Redo does not change the preference;
- Undo/Redo schedules exactly once and does not loop;
- auto off performs zero worker renders during edit and Undo/Redo;
- auto off/paused browser flows show `Update layout`; one click runs exactly once through the no-analysis path, duplicate clicks coalesce, failure remains retryable, and analysis stale hides/disables the action and runs zero times;
- failed/stale update leaves saved/exported committed artifacts unchanged;
- zoom/pan and valid stable selection survive automatic commit;
- missing selection is cleared without disturbing unrelated editor state;
- cold first sample is ignored; two accepted updates above 2,000 ms pause; Resume clears the transient pause;
- scheduler churn emits no analytics result;
- automatic render status is accessible and non-blocking;
- analysis stale shows the correct explicit action for the active workflow and performs no geometry render;
- the two superseded visible toggles are absent and one `Auto-update layout` owner remains;
- the initial unset state enables automatic layout only after a successful manual Generate, while explicit/restored off remains off;
- fresh/Reset unset and a field-absent historical session omit `ui.autoLabelReflow` on save; save/load retains `unset`, a failed Generate leaves it unset, and the first successful Generate enables it;
- the offline smoke test exercises the consolidated policy without accessing removed state fields.

#### Completion gate

History, session, export, viewport, UI status, and analytics all describe the same accepted Result and user intent. The consolidated UI is exposed only after these checks pass.

### Phase 7: integration evidence and diff audit

Status: pending

#### Browser matrix

| Mode/topology | Required cases |
| --- | --- |
| Circular single | display start, label reconciliation, placement, direct patches, selection/viewport |
| Circular grid | per-record display start/placement, stable result topology |
| Circular batch | inactive output preservation, result selection, session/export |
| Linear | per-record display start, manual placement, labels, uploaded comparison, LOSAT-backed committed comparison |
| Session replay | auto preference on/off, pending draft, analysis stale, no render during hydration |
| Races | rapid edits, direct patch during staging, manual Generate/Cancel during render, queued feature extraction across worker recreation, Reset/import/disposal |
| Failures | worker error, invalid catalogue, sanitizer failure, canonical adoption failure, stale response |

#### Focused gate

```bash
node --test tests/web/responsive-preview-policy.test.mjs tests/web/responsive-preview.test.mjs
node --test tests/web/preview-runtime.test.mjs tests/web/run-analysis-simple-path.test.mjs
node --test tests/web/history.test.mjs tests/web/history-inputs.test.mjs
node --test tests/web/session-authority.test.mjs tests/web/session-draft-authority.test.mjs tests/web/session-request.test.mjs
node --test tests/web/svg-sanitization.test.mjs tests/web/diagram-generation-worker.test.mjs
python tools/prepare_browser_wheel.py
npx playwright test tests/web/responsive-preview.playwright.spec.js --project=chromium --workers=1
npx playwright test tests/web/linear-multi-record.playwright.spec.js --project=chromium --workers=1
npx playwright test tests/web/webapp-performance.playwright.spec.js --project=chromium --workers=1
pytest tests/test_web_packaging.py -q
pytest tests/test_session_io.py tests/test_session_compat.py tests/test_api_session.py tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py -q
python tools/verify_gui_offline.py smoke-test
```

If Node Playwright is unavailable, run the equivalent targeted flow with Python Playwright. If Chromium is blocked by the agent sandbox, rerun the same check with the required approval.

#### Broader gate

Run after the focused gate passes:

```bash
node --test tests/web/*.test.mjs
pytest tests/ -v -m "not slow"
python -m build
python tools/verify_gui_offline.py inspect-wheel dist/gbdraw-0.14.0b0-py3-none-any.whl
python tools/verify_gui_offline.py smoke-test --wheel dist/gbdraw-0.14.0b0-py3-none-any.whl
```

The `--wheel` flow extracts into a temporary directory, serves only the packaged `gbdraw/web/` tree, blocks runtime network access, and cleans up after the browser closes. Prepare the browser wheel before any browser run that loads the worker. Run the built-wheel smoke after the final Python/package changes; a source-tree smoke does not substitute for it.

Do not refresh the cache-bust token unless preparing a deployable Web bundle.

#### Structural and performance success criteria

- LOSAT executor calls, cache probes, mutable input reads, helper-Pyodide initialization, and analysis-resource mutations are zero for direct/render-only updates.
- One debounce burst dispatches one worker render.
- At most one render is active and one latest draft is pending.
- At most one shared diagram-worker job is active; explicit manual cancellation releases its owned slot promptly and queued unrelated feature extraction survives worker recreation.
- A candidate commits only when its render/lifecycle revisions and replayed direct-edit revision are current; a replay mismatch queues one latest replacement.
- One accepted Result crosses the shared full-string sanitizer at most once.
- One render-only transaction performs no full deep clone of unchanged committed resource payloads during projection, worker preparation, canonical adoption, or current-history artifact refresh.
- One direct action serializes and updates an affected Result at most once.
- One continuous input interaction creates one history item.
- WSSV history and SVG-ingestion measurements remain within the completed performance-plan gates.
- Slow-pause behaviour follows D9 exactly.
- Automatic rendering itself changes no renderer geometry or tracked reference SVG.

#### Diff audit

Review separately:

1. production JavaScript/HTML;
2. Node and browser tests;
3. documentation;
4. generated artifacts.

Confirm that `tests/reference_outputs/`, Gallery sessions/media, `dist/`, `gbdraw.egg-info/`, and the generated browser wheel have no unintended changes. C/D may have their own reviewed geometry/reference changes; Work package E orchestration does not justify another refresh.

#### Completion gate

All focused tests, the applicable broader gate, browser matrix, performance counters, offline verification, and diff audit have recorded evidence. Only then may this plan status become `completed`.

## 7. Dependency and ownership order

```text
Work package C final request fields
  + Work package D final request/identity fields
  + Performance plan Phase 1 history separation
  + Performance plan Phase 2 SVG ingestion/commit
    ↓
E Phase 1 effect policy
    ↓
E Phase 2 shared candidate transaction
    ↓
E Phase 3 no-analysis executor
    ↓
E Phase 4 scheduler/manual arbitration
    ↓
E Phase 5 direct patches and session preparation
    ↓
E Phase 6 history/viewport/UI exposure/slow guard
    ↓
E Phase 7 integration evidence
    ↓
Performance plan Phase 5 reassessment
    ↓
Work package F navigation polish
    ↓
final documentation and Gallery recapture
```

Performance-remediation Phases 3 and 4 are mostly independent. Its Phase 5 must be reassessed after direct-patch convergence and must not create a second colour-update path.

## 8. Stop conditions and risks

Stop the current phase and fix its owner when any of these occurs:

- render-only execution cannot prove zero LOSAT/cache/input-discovery work;
- a request uses mutable draft files or comparison state instead of committed resources/evidence;
- an unknown field enters automatic rendering;
- a candidate mutates committed state before all fallible validation finishes;
- rollback does not restore Result, catalogue, metadata, canonical artifacts, selection, viewport, and export authority together;
- stale render/lifecycle or direct-edit replay revisions can commit;
- more than one shared diagram-worker job becomes active, or an unrelated queued feature-extraction job is lost across explicit cancellation;
- ordinary supersession terminates the worker;
- history again clones Results/catalogue per input event or renderer completion pushes an entry;
- projection, canonical adoption, or current-history artifact refresh deep-clones unchanged large resources for each accepted update;
- session hydration starts automatic rendering;
- the current-session UI inventory would change after session 41 is main-/release-backed without assigning the next evidence-backed session version;
- a second sanitizer, serializer, request builder, renderer, or state authority appears;
- performance targets can be met only by weakening sanitization, disabling large-diagram behaviour, reducing Undo, or dropping metadata;
- E unexpectedly changes Python/CLI defaults, render geometry, output topology, reference SVGs, or Gallery artifacts.

Main residual risks are large batch Results, restoring stable annotation selection, interaction between an in-flight preview and manual LOSAT completion, and supported historical sessions containing pending palette drafts. Cover each with an explicit fixture before marking the owning phase complete.

## 9. Scope cuts and non-goals

Cut in this order when necessary:

1. overlap-tolerance automatic rendering if Work package D cuts the setting;
2. slow-pause tuning, while keeping manual `Update layout`;
3. later allowlist candidates such as fonts, arrows, tracks, annotation lanes, and ticks;
4. direct annotation and general legend-content patching, which are already outside the release-blocking minimum.

Do not cut the C/D display-start and placement path, existing safe direct patches, no-analysis boundary, 400 ms single-flight/latest-wins scheduler, atomic revision/replay commit, history-neutral derived rendering, session/export authority, or minimum viewport preservation.

This plan does not add:

- a frontend framework, build step, or state manager;
- a second worker/Pyodide renderer;
- an argv-shaped Web request;
- automatic LOSAT confirmation;
- general live rendering for every control;
- a new annotation data model or DOM patcher;
- Work package F cursor-centred zoom or zoom-to-selection controls;
- a new session/request schema solely for preview runtime state;
- reference-output or Gallery regeneration without a separate reviewed owner.

## 10. Evidence ledger

For each completed phase, append a record in this document:

```text
Phase:
Behavior implemented:
Production owners changed:
Evidence commands and results:
Browser fixtures/viewports:
Performance counters:
Reference/Gallery status:
Deviations from plan:
Remaining risk:
```

Do not mark a phase complete from code inspection alone when its gate names runtime, browser, session, race, or performance evidence.

## 11. Definition of done

Responsive Preview Phase 1 is complete only when all of the following are true:

1. C/D semantic settings and performance P0 prerequisites have passing evidence.
2. The effect registry is explicit, composable, and fail-closed.
3. Proven-safe direct edits remain immediate and session/export-reproducible.
4. Per-record display start, placement, retained tolerance, and label reconciliation use the dedicated no-analysis renderer path.
5. Empty-cache/spying tests prove zero analysis, cache, input-read, and analysis-resource work.
6. One scheduler owns automatic and manual diagram-render dispatch with 400 ms debounce, one active request, one latest pending draft, separated revision domains, and current-revision/replay commit.
7. Failed, cancelled, stale, and superseded work preserves the complete committed bundle.
8. Renderer completion is history-neutral, and Undo/Redo reconciles restored intent once.
9. Session load migrates every supported historical preference/palette shape, preserves explicit off, removes the retired Palette flag from current outputs, and restores draft/committed mismatch without automatically rendering during hydration.
10. Automatic commit preserves Result topology, zoom/pan, and valid stable selection across Circular single/grid/batch and Linear.
11. Slow-pause and analytics outcome semantics match this plan.
12. Focused tests, applicable broader tests, offline verification, performance evidence, and separate diff audits pass without unintended reference or Gallery changes.
13. Every admitted direct-edit family survives a later geometry update through canonical representation or stable-identity replay.
14. One shared-worker broker serializes renderer and feature-extraction jobs; explicit Cancel releases its owned active slot promptly and preserves queued unrelated work across worker recreation.
