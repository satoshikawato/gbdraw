# Web application maintenance

Read the repository-level `CLAUDE.md` first. This file describes the intended
Web architecture and its ownership boundaries. It deliberately avoids line
numbers and implementation-size snapshots because those become stale quickly.

## Non-negotiable properties

- The Web UI is a single-page app with no JavaScript build step.
- `index.html` owns markup, local styles, and templates. ES modules live under
  `gbdraw/web/js/`.
- All runtime assets are served from the same origin. Do not add a CDN or other
  runtime network dependency.
- Genome inputs and rendered results stay in the browser. Generation must not
  upload user data.
- Generated SVG is sanitized before insertion into the preview.
- The browser wheel is generated and gitignored. Prepare it when packaging or
  wheel-dependent tests need it; never edit it by hand.

## Runtime data flow

The canonical generation path is:

```text
index.html
  -> app.js and app/app-setup.js
  -> reactive state
  -> services/session-request.js builds one canonical schema-9 render request
  -> app/run-analysis.js validates and orchestrates optional LOSAT work
  -> services/diagram-generation.js dispatches the request and resources
  -> workers/diagram-generation-worker.js owns Pyodide rendering
  -> typed request decoding and render_request()
  -> sanitized preview, interactive editing, export, and session save
```

Do not introduce a second argv-shaped generation contract. Fresh generation,
saved sessions, and replay must converge on the typed render-request boundary.
JavaScript owns browser state and resource bytes; Python owns request validation,
planning, loading, diagram assembly, and rendering.

For an owner or canonical-path change, apply the
[Product Impact Ratchet](../../docs/internal/PRODUCT_IMPACT_RATCHET.md) and trace
the affected architecture subjects through their mapped requirements, user
effects, checkpoints, authority, and behavior contracts. Select the supported
behavior before choosing its implementation location. Preserve every jointly
required contribution; matching a high-level option ID alone is not evidence
of preservation.

The application has one Pyodide runtime, owned by the lazy diagram Worker.
Typed helper operations, generation-time feature extraction, and diagram
rendering share that runtime; Python must not return to the main thread.
The Worker may retain a bounded cache of parsed biological sources, resolved
records, and interactive metadata. Key it only by validated resource tokens and
the semantic fields used to prepare those values. Every Generate still decodes
and validates the request, builds the drawing, writes the SVG, and builds the
feature catalog. Publish cache fills only after the full render succeeds, and
release them with the Worker; never cache drawings, catalogs, SVG, or final
results.

App-shell readiness means Vue is mounted and required local UI assets and
palette definitions are available. It is independent of diagram-Worker
readiness. Loading a saved preview that needs no Python must leave the Worker
unconstructed. The first Python-backed helper or render operation constructs
and initializes the Worker, and later operations reuse it. Callers await their
own explicit success, error, or canceled result. Tests and capture tools start
the operation they need and await its settlement; they never wait for a global
Worker-ready flag.

LOSAT workers are separate from the diagram-generation worker. Threaded LOSAT
needs a cross-origin-isolated page, and `losatThreadingPrecondition()` in
`services/losat.js` is the one check. An explicit Threaded choice without it
fails with `LOSAT_THREADING_UNAVAILABLE`; Serial and Auto keep the single-thread
path.

## Module ownership

This table records current ownership. Changes to an owner or canonical path
must follow the repository
[architecture fitness-function ratchet](../../docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md),
using concise owner/path evidence for ordinary non-increasing changes and
complete before/after sets only for defined exceptions. Remove superseded paths
in the same change.

| Area | Owner |
|---|---|
| Vue mount and exports | `js/app.js` |
| Composition and dependency wiring | `js/app/app-setup.js` |
| Reactive state and computed values | `js/state.js` |
| Each diagram mode's drawing (settings, editor edits, derived values) | `js/state.js` (`drawings`, `activeDrawing`; R15) |
| Registry of what a mode's Session slice holds | `gbdraw/web_support/mode_scoped_settings.py`, generated into `js/mode-scoped-settings.generated.js` by `tools/generate_mode_scoped_settings.py` |
| Session 46 mode slices and the split of a Session 27–44 draft into them | `js/services/mode-scoped-migration.js` (Python twin: `split_draft_into_modes` in `gbdraw/session_io.py`) |
| Generate-button orchestration | `js/app/run-analysis.js` |
| Canonical request/session projection and equivalence | `js/services/session-request.js` |
| Current active-config defaults, inventory, and validation | `js/services/session-active-config-contract.js` |
| Historical Gallery session migration | `js/services/gallery-session-migration.js` |
| Per-feature edits and Feature placements by source identity (draft key, row checks, request rows, notices) | `js/services/feature-placement.js` |
| Rendered-ID feature edits of Sessions 44 and older, and Feature placement drafts of Sessions 41–44, onto mode-scoped identity rows | `js/services/feature-edit-migration.js` |
| Gallery publication preparation, finalization, and readiness | `js/services/gallery-session-publication.js` |
| Save/load coordination | `js/services/config.js` |
| History transactions and availability | `js/services/history.js`, `js/app/history-inputs.js`, `js/services/history-snapshot.js` |
| Error producer contract and user-facing wording | `js/utils/error-normalization.js` (`diagnosticError`, `normalizeUserFacingError`) |
| Numeric draft projection | `js/utils/optional-positive-number.js` (`projectOptionalNumber`) |
| GenBank header reading | `js/services/genbank-header.js` |
| Managed Depth slot reconciliation | `js/services/depth-track-state.js` (`reconcileManagedDepthSlots`), called by `changeCircularDepthSources` and `changeLinearDepthSources` |
| Editor-intent projection onto a displayed Result | `js/app/app-setup.js` (`projectMountedEditorIntent`) |
| Render-worker client and lifecycle | `js/services/diagram-generation.js` |
| Pyodide typed rendering | `js/workers/diagram-generation-worker.js` |
| LOSAT dispatch and workers | `js/services/losat.js`, `js/workers/losat-*` |
| SVG/PNG/PDF downloads | `js/services/export.js` |
| SVG Result admission and sanitization | `js/services/svg-result-ingestion.js` |
| Feature-editor entry point | `js/app/feature-editor.js` |
| Feature-editor helpers | `js/app/feature-editor/` |
| Right-side editor drawer state and transitions | `js/services/right-drawer-state.js`; the drawer controller is `js/app/right-drawer.js` |
| Legend entry point | `js/app/legend.js` |
| Legend helpers | `js/app/legend/` |
| Legend/diagram positioning | `js/app/legend-layout.js`, `js/app/legend-layout/` |
| Help tips | `js/components.js` (`HelpTip`), template in `index.html` |
| Public Gallery inventory and asset generation | `tools/prepare_interactive_gallery_assets.py::EXAMPLES` |
| Gallery refresh process, staging, and replacement | `tools/refresh_gallery_sessions.py` |
| Gallery data, media, and tutorials | `gallery/` |

During an artifact History transaction, the captured `before` checkpoint is the
current checkpoint until the `after` capture replaces it. Do not retain an older
equivalent current checkpoint alongside the transaction copy.

Generate is an already-applied immutable artifact replacement. It must not use
full checkpoint cloning or signing, and its pre-state is the sole rollback authority.

Keep top-level `create*` entry points in `js/app/*.js`. Put larger,
single-purpose helpers in the matching subfolder instead of growing another
general utility module.

## Reactive availability and live-edit invariants

Treat a UI selection from a dynamic domain as valid state, not as an unchecked
preference.

- Keep one canonical visibility value and one canonical selected value. Do not
  add per-tab visibility flags or other mirrors of those values.
- Give one feature owner all open, toggle, close, reset, restore, availability,
  and reconciliation transitions. Templates, lifecycle code, session code, and
  global events call those transitions; they do not assign the refs.
- Use the same availability predicate for enabled rendering and action
  resolution. An unavailable or unknown request resolves to a deterministic
  available fallback. An enabled action must not silently return without a
  state change.
- Keep the selected value valid while the control is open or closed. Reconcile
  capability loss synchronously. Capability gain does not replace a currently
  valid selection.
- Close and Escape change visibility only. Successful document replacement
  resets transient UI state. Failed replacement restores source data first,
  then restores and reconciles transient UI state against that data.
- A watcher may invoke the owner reconciler for external mutations, but watcher
  execution is not an invariant mechanism. Bulk operations may suppress or
  coalesce watchers and must call the same transition explicitly when needed.
- Derive mount and visibility conditions from canonical state. Do not introduce
  a second ref solely to trigger editor synchronization.

The right-side editor is a live-edit boundary. A single Feature, Label, Legend,
or visibility edit must commit its canonical override state before returning;
Generate, session save, and drawer close are never live-edit triggers. Route the
mutation through the owning editor action so History records the same operation.
If a mounted target exists, update it and the current Result synchronously before
optional geometry reflow. If canonical output requires creating or removing
geometry and no mounted target exists, queue the owning automatic rerender in
the same action and replace the Result when it completes. That rerender renders
the committed Session plus the current editor-intent tables
(`projectCommittedEditorIntent` in `services/session-request.js`) and never reads
draft settings (R1). A rerender failure must keep any direct edit already
applied and report the failure. Metadata-only edits
with no static SVG target, such as Similarity group names and descriptions,
update canonical state only.

## Design rules R1-R15

Each rule names one owner and the guard that enforces it. Cite the rule id in PR
descriptions and tests. Change the code and its guard together, and lower a
shrink-only baseline in the change that removes the last site it counted.

### R1: Three writers of the Result

The current Result changes only through:

- (a) the Generate compiler (`compilePlanBundle` in `app/candidate-render.js`,
  then `services/svg-result-ingestion.js`) re-applying editor intent with the
  same executor;
- (b) a composition edit that keeps parity with the Python renderer;
- (c) an automatic rerender from the committed Session and declared projection
  fields.

A setting that none of these applies is **Applies on Generate**: say so beside
the control, write only the draft, and leave the Result unchanged. Never patch
the Result from a draft or with a helper's SVG fragment. A new Python helper that
returns SVG text (`DIAGRAM_HELPER_OPERATIONS` in
`services/diagram-worker-protocol.js`) needs architecture review.

Guards:

- the single `app/candidate-render.js` to `services/svg-result-ingestion.js`
  edge (`canonical-path.current-result-admission` in
  `tools/web-architecture-rules.json`);
- the `Result content commit` owner list in `tools/web-change-policy.json`:
  an unlisted Result writer fails `node tools/check-web-change-budget.mjs`, and a
  production change may only remove entries;
- `tests/web/history-generated-authority.playwright.spec.js` and the label-reflow
  test in `tests/web/gui-audit-20260930-editor.playwright.spec.js` (a draft edit
  leaves the Result unchanged until Generate).

### R2: Lifetime of editor intent

Editor overrides are created, pruned, or deleted only by an explicit Reset or
Import, Undo or Redo, a Session replacement, the owner reconcile inside a
successful source-replacing Generate (`pruneUnmatchedFeatureOverrides` in
`services/feature-visibility.js`), or **Remove N unmatched feature edits**
(`removeUnresolvedFeatureEdits`). Result selection, mount, record selection,
mode change, hiding, and reflow never touch them. Per-feature edits
(`featureOverrides`) and Feature placements live in the drawing of the mode
they were made in (R15) and are keyed by
`JSON.stringify([recordKey, biologicalFeatureId])`. Both modes can use the same
record key for the same feature (`record-1` in a Gallery or Python Session);
each drawing holds its own row, so a request, the live projection, and every
reconcile read only the rows of their own drawing and records, and the other
drawing keeps its rows. Only `services/feature-placement.js` builds a draft
key. A selected-feature annotation target (`featureIdentity`) belongs to its
drawing the same way, and `annotationOptionsPayload` in
`services/annotation-state.js` gives a request only the targets of its records.
A row or target of a Session 41–44 that names its mode (`scope`) keeps that
field only until the Session 46 split moves it into that mode's drawing.

Guards: `tests/web/non-edit-state-preservation.playwright.spec.js` with
`tests/web/non-edit-state-diff.test.mjs` (user-owned state is identical before
and after non-edit operations); the disjoint-Result binding test in
`tests/web/feature-label-visual-unit.test.mjs`;
`tests/web/mode-scoped-feature-edits.test.mjs` (no draft row reaches a request,
reconcile, or TSV load of the other mode, and no other module builds a draft
key) and `tests/web/mode-scoped-feature-edits.playwright.spec.js`;
`tests/web/annotations.test.mjs` and the OV-21 test in
`tests/web/annotation-feature-identity.playwright.spec.js` (annotation targets).

### R3: One projection per domain

A live action, a History apply, and a Result display call the same projection:
`projectMountedEditorIntent` (palette, rules, visibility, labels), which calls
the domain projections `projectPaletteAndRules` in `app/app-setup.js` (palette
and specific rules) and `projectFeatureVisibility` in
`app/feature-editor/visibility-actions.js` (feature visibility and the labels of
the features it hides or shows),
`orderLegendEntries` in `services/legend-svg.js` (legend order), and
`resolveFeatureDrawn` in `services/feature-visibility.js` (whether a feature is
drawn, as Python's `should_render_feature` answers; the Features list and
Search features read it through `listFeatureRows`). A displayed batch Result whose shared legend entries already follow
the legend order keeps its order (`orderLegendEntries` with `keepFollowed`),
so the entries only that Result draws keep their places; a Result last shown
with another order also receives the default order. A History legend step made
on another batch Result projects only its shared legend intent onto the
displayed Result (`reconcileLegendEntries` with the step's other side, B19), so
no Result gains or loses an entry that only one Result draws. The displayed
population comes from the mounted Result's committed metadata
(`renderedFeatureIdentities`), not from a second selection ref. Per-feature
edits reach a Result through its catalog features' identities
(`resultRenderedFeatures` in `services/feature-catalog.js`), so every copy of
a feature in every Result shows the same edit. A dialog's
reactive object holds display values only, never a copy of an owner's data.

One projection means one call shape: a domain's projection function is called
from outside its owner in one form, not re-implemented or re-sequenced per
caller. The distinct call shapes per domain may only decrease (R13).

A live edit shows on the Result what the next Generate draws (PD-OI-066). A
change to a live-edit path is done when the case for its edit kind in
`tests/web/live-generate-parity.playwright.spec.js` passes with Auto Reflow on
and off and with the popup opened and not. A difference the Result may show
until Generate is listed with its decision in
`tests/web/contracts/live-generate-parity-allowed.json`. A mismatch the matrix
finds is a bug: the case stays `test.fail` with its finding id until the fix
makes it pass.

Guards: `tests/web/gui-audit-20260930-editor.playwright.spec.js` (a Result shows
the edits made on another Result, also after Undo, Save, and Load),
`tests/web/feature-visibility-actions.test.mjs`,
`tests/web/feature-color-actions.test.mjs`,
`tests/web/live-generate-parity.playwright.spec.js` (the matrix of edit kinds
and editor states) with `tests/web/live-generate-parity-legend.playwright.spec.js`
and `tests/web/live-generate-parity-track-data.playwright.spec.js` (Legend edits,
track-data changes, and mode switches): a live edit equals the next Generate;
and the `projection-shapes` baseline in
`tests/web/owner-graph-baseline.test.mjs`.

### R4: A fast path matches the canonical reader or declines

A browser shortcut for a fact Python owns either equals the canonical loader or
hands the decision to the Worker. `services/genbank-header.js` is the only JavaScript
GenBank header reader, and it declines what it cannot read exactly. A Worker
helper that answers what Generate reads calls the canonical loader. A JavaScript
check that anticipates a Python match (letter case, color domains) uses Python's
equivalence.

Guard: the shared vectors in `tests/fixtures/record_metadata_inference_cases.json`
run through `tests/test_record_metadata.py` (loader and Worker helper) and
`tests/web/record-metadata-inference.test.mjs` (JavaScript). Add a vector, not a
special case, when the paths disagree.

### R5: Comparison frames and reuse

Comparison evidence uses three frames: F, the search frame (the selected and
cropped record on its source strand, which is the LOSAT output); V, the display
frame (F with the effective reverse complement applied); and Src, input-file
coordinates. Raw rows and **Save Raw LOSAT TSV** stay in F. Every comparison
table is read in F, and only the Python planner projects orientation to V.
Human-readable coordinates (popups, FASTA headers) lead with Src. The Linear
orientation owner is the File card's `region_reverse`, written by
`app/record-display-options.js`; `services/record-display-model.js`
(`effectiveRecordReverseComplement`) derives the effective orientation.

Reuse compares every input that is knowable without extraction as data; an
invalidation event is only an optimization. LOSAT job planning is one pure plan
for execution and the estimate: `planLosatSourceJobs` (`services/linear-sources.js`)
and `buildLosatJobSpecs` (`services/linear-comparisons.js`).

Guards: `tests/web/linear-sources.test.mjs` (one-file packaging equals separate
files; job plan) and `tests/web/gui-audit-20260930-comparisons.playwright.spec.js`
(reuse and orientation).

### R6: The producer owns the failure meaning

A failure the user can fix is thrown with a code and a bounded context:
`diagnosticError(code, context, { stage, operation })` in JavaScript,
`GbdrawError(..., diagnostic=)` in Python. Wording lives only in the normalizer
(`normalizeUserFacingError`). A locator in the context (Sequence N, Line N,
Column N, Track row N, Depth series N, Available band, Feature) appears in the
summary, and an offered action must be one the user can press. Never add a
classifier that matches message text.

Guards:

- `tests/web/error-producer-coverage.test.mjs`: unclassified `throw new Error`
  sites per validation owner and the message-table sizes are shrink-only.
- `tests/test_web_error_producer_coverage.py`: the Python adapter tables are
  shrink-only, and the Python vocabulary is a subset of the Web wording owner.
- `tests/web/error-normalization.test.mjs`: table-driven producer cases and
  code, reason, and context-key parity.

### R7: Python's typed layer owns value checks

The Web projection never converts a value: a blank is `null`, a finite number is
that number, and anything else is a typed `INPUT_INVALID`
(`projectOptionalNumber`). A check that must run before Python (comparison
thresholds before LOSAT) evaluates the domains generated from
`gbdraw/mode_profiles.py` into `js/mode-profiles.generated.js`
(`resolveComparisonThresholds` in `js/mode-profiles.js`). Generate never rewrites
the draft.

Guards: `tests/web/option-input-integrity.test.mjs` (every projected numeric
draft field either rejects `NaN`, `1e-50x`, `-5`, `0`, and `12.5` with a typed
diagnostic or sends them literally; the count of draft assignments on the
Generate path in `app/run-analysis.js` may only decrease) and
`tests/fixtures/option_domain_vectors.json` with
`tests/test_option_domain_parity.py` (the CLI, the typed request, and the Web give
the same field and reason).

### R8: The Multi-Record Canvas owns placement only

Legend rows, Depth presence, and slot composition come from the same function as
the single-record path. "No input" (`None`) differs from "input without a cell"
(`[]`).

Guard: `tests/test_circular_multi_record_parity.py` (a one-record canvas equals
the single-record path in slot geometry, legend rows, and definition lines).

### R9: One reader per input format

BLAST outfmt 6 tables are read only by `read_comparison_table` in
`gbdraw/io/comparisons.py`, with one normalizer
(`normalize_comparison_dataframe`). An unreadable row raises `ValidationError`
and is never skipped. GenBank headers follow R4.

Guard: `tests/test_comparison_tables.py` (`test_outfmt_table_columns_have_one_reader`
scans the Python and Web JavaScript sources for a second mapping onto
`COMPARISON_COLUMNS`, and the shared vectors run through every entry point).

### R10: A watcher does not repair state

A transition that changes the inputs of a derived structure calls the owner's
reconcile explicitly, through the port the composition root injected or by
returning to the root, which calls the reconcile (R13); it does not hold the
other owner's object. Depth sources go through `changeCircularDepthSources` and
`changeLinearDepthSources`, and the suppress controls through
`applyCircularSuppressControlsToSlots`. One concept has one builder
(`resetCircularTrackSlotsToPreset` for the simple-controls stack). Display values
are looked up by identity (`slotId`), never by position.

Every edit of a feature-slot input (Track Preset, Track Layout, Separate
Strands, Use custom stack, and the custom stack rows) runs through one
transition in `app/feature-editor/placement-actions.js`: the composition root
passes it to the stack editors as the `changeTrackLayout` port, each editor
routes its feature-slot edits through it (`featureSlotEdits` in
`app/track-slot-edits.js`), and a template passes `$event` so Cancel restores
the control. An edit that would leave a lane Feature placement of the
drawing undrawable asks first (Owner decision Q3); Reset applies the edit
and removes those rows as one History step. The control's History adapter
records an edit that loses no lane, and restores install state as is (R11).

Every edit of a Feature Visibility rule (the Features panel and the popup's
product and protein ID scopes) runs through `editFeatureVisibilityRules` in
`app/feature-editor/visibility-actions.js`, which writes the rules, projects
them, and follows with the labels as one History step.

Guards: `tests/web/depth-slot-lifecycle.test.mjs` (the same cases through the
Circular and Linear editors), `tests/web/circular-track-slots.test.mjs`,
`tests/web/track-slot-display.test.mjs`, and
`tests/web/track-layout-transition.test.mjs` (every writer of a feature-slot
input is that transition, a listed restore, or a listed reconcile) with
`tests/web/track-layout-transition.playwright.spec.js`, and
`tests/web/feature-visibility-rule-writers.test.mjs` (every writer of the
Feature Visibility rules is that transition, a History restore, a Session load,
or Reset Settings).

### R11: History transaction boundaries

- A transaction belongs to one owner, a control or a gesture. Starting another
  owner's transaction first settles the open one. A legend or diagram drag passes
  a per-gesture owner.
- `settlePendingIntent` in `services/history.js` is the one settlement of the
  open intent. `begin`, `beginCheckpoint`, `beginArtifactReplacement`,
  `runUndoableCommand`, `undo`, and `redo` share it.
- A discrete control begins in the capture phase of the event that commits its
  value: `change` for a checkbox, radio, or file input, `click` for a button. A
  text control begins on focus. The pointer, the label, the keyboard, and a file
  picker opened by another button record the same one step
  (`app/history-inputs.js`).
- A file input whose change ends with asynchronous work is
  `data-history-managed` and owns one `runUndoable` step that commits after
  that work: a `file-uploader` `afterChange` import, and **Add Seq**
  (`addCircularConservationComparisonFile`), whose step holds the ring record
  label read by the Python reader.
- Undo, Redo, and the header buttons share one reactive availability
  (`historyAvailability`), busy while an artifact replacement or checkpoint is
  open or an action owns the open intent.
- The History intent holds every file binding by reference
  (`buildIntentFilesData` in `services/history-snapshot.js`), the BLAST rows and
  comparison sequences of a LOSAT-cache replay included: its ring rows and the
  managed track slots' `series_key` name those rows. Only the generated Linear
  comparisons stay in artifact checkpoints.
- The feature catalog of the current Results (Generate's, or the automatic
  rerender's, which replaces it with its Results) is held by reference in checkpoints, history
  entries, and Session rollback, and is never cloned through JSON.
  `state.featureCatalog` is null or an admitted catalog (`admittedFeatureCatalog`
  in `services/config.js`). The two paragraphs under Module ownership govern
  artifact transactions.
- A restore (Undo, Redo, the failed Session Load rollback) installs a copy of
  the captured state as is: `applyConfigData(..., { resolveTrackPlacements: false })`
  and `applyEditorStateData(..., { normalized: true })`. Unset values stay unset;
  resolving them belongs to the request builder, not to the restore.
- The History intent holds neither the selected Result nor the offsets read
  from the mounted root (`buildUiIntentData` in
  `services/history-snapshot.js`). Composition offsets are recorded per Result
  by committed identity (`captureCompositionIntent` in `app/legend-layout.js`),
  and Undo and Redo restore the Result a step was made on through
  `commitResultEdit` in `app/preview-runtime.js`, also while another Result is
  displayed.
- The Result picker is navigation (`data-history-ignore`), so a Result switch
  records no step, although the extracted legend entries follow the displayed
  Result.

Guards:

- `tests/web/history-inputs.test.mjs` and
  `tests/web/history-inputs.playwright.spec.js`: control kind, input means, focus
  state, and last-in-first-out Undo; each changed control adds exactly one step.
- `tests/web/history-generated-authority.playwright.spec.js`: checkpoint Undo and
  Redo, a mode round trip, and Undo while Generate runs.
- `tests/web/multi-result-edit-matrix.playwright.spec.js`: Undo and Redo of
  drags on one batch Result while another Result is displayed (B17), a
  Result switch after a legend sort that records no step (B18), Undo and Redo
  of a legend sort made on another Result (B19), and Sort by default reaching
  another Result (B20).
- `tests/web/history-config-restore.test.mjs`: catalog identity through Undo,
  Redo, and Session rollback, and unset settings that stay unset through intent
  and checkpoint Undo and Redo (also `session-draft-authority.test.mjs` for the
  rollback).
- `tests/web/non-edit-state-preservation.playwright.spec.js` (G-H rows): Undo^n
  Redo^n through mixed edits and a mode round trip, and Reset Settings Undo and
  Redo, leave user-owned state unchanged.
- `tests/web/session-operation-consistency.test.mjs`: Undo availability during
  Generate.

### R12: Source recipe and accessibility

- Run Info builds the Source recipe argv and reads it back with the CLI's split
  rules (`assertCliSlotTokenLossless` in `app/run-info.js`, following
  `gbdraw/tracks/parsing.py`). A recipe that does not read back identically is
  unavailable with `sourceRecipe.unavailableReason`; it is never emitted lossy.
- Give every input an accessible name that matches its visible label. A
  placeholder or a title that changes with state is not a name. A help tip is a
  focusable button outside the `<label>`, and the control it explains references
  the tip text through `aria-describedby`.

Guards: `tests/web/run-info.test.mjs` (including slot labels with `,` or ` #` and
a Linear scale font without a ruler-label font) and
`tests/web/accessibility.playwright.spec.js` (every visible input in both modes
has an author-provided name that does not change with its state, and every help
tip is a reachable disclosure).

### R13: Owner layering and ports

Modules sit in layers, lowest first: leaves (`utils/`, `config.js`,
`web-ux-profile.js`, `mode-profiles.js`, `mode-profiles.generated.js`,
`mode-scoped-settings.generated.js`), then
state-free `services/`, then
`state.js`, then the two state-bound services (`services/config.js`,
`services/reset.js`), the only services that import `state.js`, then owner
modules (`app/*.js` and `app/*/` factories), then the composition roots
(`app/app-setup.js`, `app/feature-editor.js`, `app/legend.js`,
`app/legend-layout.js`), then `app.js`, `components.js`, and the template. A
module imports its own layer and lower ones; a state-free module never imports
`state.js`, and an owner module never imports a composition root.

- An owner never imports, holds, or calls another owner object. It receives
  ports: a port is one function, in one direction, named for the reaction it
  performs (`applyFeatureVisibilityToLabels`), not for its provider
  (`labelActions`). A factory that takes a whole owner object outside a
  composition root is a whole-object port.
- A reaction of owner A to a change in owner B is wired by a composition root:
  as a port injected into B, or as the root calling A's projection after B's
  transition returns. B does not call A's object, and A does not call back into
  B.
- A closure passed at composition that resolves to an owner created later in
  the same root (a forward-reference closure) is a composition error: reorder
  the creation, or register the port after both owners exist. A closure at any
  depth counts, and so does a late-bound function variable that the root
  reassigns after the owner exists.
- `state.<name> = <function>` outside `state.js` is a backdoor: a value an owner
  needs from a higher layer arrives as a port or on an existing data path
  (Result metadata, the admitted catalog).
- A heavy derived value (`rulePreparation`) is triggered from one place per
  flow: the composition root's projection, the Generate compiler, or the owner
  of its input. An owner that only reads the result awaits `pending` rather
  than triggering it. A trigger site includes a call of a producer port (a
  producer method handed to an owner as a function, such as
  `rulePreparation.run` or `rulePreparation.runDrawn`) in the module that
  receives it.
- Each projection domain has one projection function and one call shape
  outside its owner (R3).

The baseline of the current violations is recorded in
`tests/web/owner-graph-baseline.test.mjs` and may only shrink: a runtime change
removes the subjects it fixes from the baseline in the same pull request, and
an addition is an authority-only change (`tools/web-design-rule-guards.json`,
WEB_CHANGE_POLICY.md "Design-rule co-change"). The detectors are
`tools/web-owner-graph-detectors.mjs` (injection edges, forward-reference
closures `owner-graph.forward-closure.v2`, state backdoors, whole-object ports,
projection call shapes, heavy derived trigger sites
`heavy-derived.trigger-site.v2`, and import direction
`layer.import-direction.v1`, which ranks modules with `webLayerOf`; no upward
import is allowed, and none is recorded);
`node tools/report-web-owner-graph.mjs --at worktree` prints the current
subjects. Fix an upward import by moving the imported code down or taking a
port.

Run `node --test tests/web/owner-graph-baseline.test.mjs` before and after a
change to `gbdraw/web/js`. A new subject means the change couples owners:
restructure it instead of recording the subject. Renaming a parameter or a call
so that a subject disappears is not a fix; the pull request states which owner
no longer holds which object, and which port replaced it. A reaction of one
editor domain to another follows the decision record's Reaction owner, Channel,
and Projection (`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`); without
a record, ask the Owner before implementing it.

Guards: `tests/web/owner-graph-baseline.test.mjs` (every observed subject is in
the baseline with a count no higher than recorded, per detector, and no module
imports from a higher layer) and
`tests/web/owner-graph-detectors.test.mjs` (the detectors).

### R14: Typed boundaries

Modules under `js/` except `workers/` are checked by `tsc --noEmit` with JSDoc
types (`tests/web/types/tsconfig.json`); there is no build step and no
TypeScript source. A checked module starts with the line `// @ts-check`. The
modules not yet checked are listed in `UNCHECKED_MODULES` in
`tests/web/typed-boundaries.test.mjs`, which may only shrink: the pull request
that checks a module removes it from the list, and an addition is an
authority-only change. A new module is checked from its first commit.

- A checked module declares the type of every parameter of its exported
  `create*` and `setup*` factories. A port (R13) is a function type declared by
  the owner that receives it; a composition root imports it with `@import`. A
  type naming another owner's factory result (`ReturnType<typeof createX>`) is
  a whole-object port.
- A persisted contract (the Session, the canonical render request, the feature
  catalog) has one typedef, in its writer beside its version constant
  (`SESSION_VERSION`, `CANONICAL_REQUEST_SCHEMA`, `FEATURE_CATALOG_SCHEMA`),
  for the current format only. A reader takes unvalidated data and returns that
  typedef.
- Types are JSDoc in the module that owns them, imported with
  `/** @import { T } from './x.js' */`, never a `.d.ts` or a types module under
  `gbdraw/web/`. A type import follows the layers of R13. The globals of
  vendored scripts and the test hooks are declared once, as `any`, in
  `tests/web/types/web-globals.d.ts`.
- `@ts-ignore`, `@ts-expect-error`, and `@ts-nocheck` are not used: fix the
  type, cast with JSDoc, or fix the code.
- `noImplicitAny` is not yet on in the guard config.
  `tests/web/types/tsconfig.no-implicit-any.json` is the guard config plus
  `noImplicitAny`, and `NO_IMPLICIT_ANY_BASELINE` in
  `tests/web/no-implicit-any-ratchet.test.mjs` caps its diagnostics across
  `js/`. A pull request keeps the count at or below the cap and may lower the
  cap to the new count: a new parameter or variable declares its type, or the
  collection a callback iterates is typed. Raising the cap is an
  authority-only change.

Guards: `tests/web/typed-boundaries.test.mjs` (the checked set, the compiler
run, declared factory parameters, no suppression, type-import direction) and
`tests/web/no-implicit-any-ratchet.test.mjs` (the `noImplicitAny` cap). They
run the `typescript` version pinned in `package.json` and fail until `npm ci`
has installed it.

### R15: Drawing context and mode-scoped settings

Each diagram mode has one drawing, `state.drawings.circular` and
`state.drawings.linear` (`createDrawingState` in `state.js`). A drawing holds
that mode's settings (`form`, `adv`, the LOSAT search fields, the comparison
plan, record display, the layout slot), its editor edits (per-feature edits and
placements, colors and color rules, label tables and filters, Legend edits,
annotation sets, canvas padding), and the values derived from them, under their
former `state` names. Nothing else holds them, and the two drawings share no
member object, so a setting or an edit made in one mode never reaches the other
(PD-OI-086).

- `state` keeps only what both modes share: the input files, the LOSAT caches,
  the LOSAT execution and thread settings, app settings (Auto Reflow, PNG DPI,
  palette Instant Preview), the Session title, the shown Result and its
  artifacts (each mode keeps its own Result, `services/artifact-slot.js`), and
  transient flags. Every `state` key and drawing member is in exactly one class
  in `tests/web/helpers/drawing-state.mjs`. A value both modes should share is
  a Product decision: record it in the Product Contract before moving it out
  of the drawing.
- A state-free service never resolves the active drawing: it takes the
  `drawing` it reads, or the mode whose drawing it reads (for example the
  mode of the Result it handles). An owner resolves `state.activeDrawing()`, the
  drawing of the shown mode, once per action, or the drawing of the mode its
  action is about, and passes it down. Generate resolves its drawing once, at
  the entry of `runAnalysis`; its request, validation, and commit use that
  drawing only. A reflow or rerender of a committed Result uses the drawing of
  that Result's mode.
- A mode switch (`transitionDiagramMode`) writes no drawing (R10). History
  stays one manager, and a switch is one step (R11).
- A Session 46 saves each drawing as `modes.circular` and `modes.linear`.
  `gbdraw/web_support/mode_scoped_settings.py` is the one registry of what a
  slice holds (domain, path, modes, row key, migration token);
  `tools/generate_mode_scoped_settings.py` writes it into
  `js/mode-scoped-settings.generated.js`, which is never edited by hand. Load
  rebuilds each drawing from the projection of that mode's committed request,
  then the slice's stored values, then the mode's defaults.
  `services/mode-scoped-migration.js` splits a Session 27–44 draft into the two
  slices, and `split_draft_into_modes` in `gbdraw/session_io.py` is its Python
  twin. A Session 45 is not read.

Guards: `tests/web/drawing-context.test.mjs` (the key classes, one drawing per
mode, no drawing key read from `state`, no `activeDrawing` identifier in
`services/`, and Generate resolving it once); `tests/web/per-mode-drawings.test.mjs` (no
shared member object, and each mode's request carries its own drawing's
values); `tests/web/mode-scoped-session.test.mjs` (the Session 46 round trip,
and a Session 45 or a field outside the registry rejected);
`tests/web/mode-split-vectors.test.mjs` with `tests/test_mode_split_vectors.py`
(the JavaScript and Python splits agree on shared vectors);
`tests/test_mode_scoped_settings.py` (the registry rows, and the generated file
equal to its Python source); and
`tests/web/mode-scoped-settings.playwright.spec.js` (settings, Legend edits,
Depth, and History per mode in the browser).

## Computation ownership

Follow [Computation ownership](../../docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#computation-ownership)
(CW-01 to CW-06). Status, selection, and help rendering do not call the
canonical request builder, generated-table serialization, or full feature
metadata construction.

Current owners in the initial automated scope:
- The bulk label projection is owned by
  `services/label-override-table.js` (`buildBulkLabelProjection`).
  `services/session-request.js` builds it once per request
  (`requestLabelProjection`) and writes its identity-row label text and its
  `* * label` table rows into the canonical request.
- Observation uses `services/runtime-test-hooks.js` (`recordStructuralMetric`),
  History `getDiagnostics()`, and the Worker tracking in
  `tests/web/helpers/app-lifecycle.cjs`.

## Request and session boundary

`services/session-request.js` is the projection boundary between reactive UI
state and the persisted/rendered model. It owns:

- canonical request schema and resource descriptors;
- canonical publication form, resource identity, and request equivalence;
- Circular `single`, `grid`, and `batch` grouping;
- Linear record/group topology;
- output-prefix and per-batch output projection;
- explicit track-slot projection;
- per-feature edit and Feature placement rows (`diagramOptions.featureOverrides`
  and `featurePlacements`, built by `services/feature-placement.js`);
- the layout-preference slot that a committed request sets in its mode's
  drawing (`modes.<mode>.ui.layoutPreferences` in a Session 46).

`services/config.js` coordinates user-facing save and load actions. It must not
grow a parallel model of render fields. Compatibility migrations are reader
concerns: normalize supported old data once, then use the current model. Do not
write retired fields into new sessions.

`services/session-active-config-contract.js` is the DOM-free current writer
contract used by normal restore and Gallery publication. Historical Gallery
migration stays in `services/gallery-session-migration.js`; it does not clean up
current sessions. `services/gallery-session-publication.js` may align an admitted
Gallery draft with its committed render intent, but ordinary import keeps its
saved Result and active draft separate.

Explicit track slots are authoritative when enabled. Legacy flat controls are
compatibility inputs, not a second source of truth. Preserve empty positions in
per-record depth and comparison inputs when those positions carry alignment
meaning.

See `docs/REFERENCE/session-and-request-compatibility.md` for accepted versions
and migration limits, and `docs/internal/SESSION_FORMAT_HISTORY.md` for the
version-by-version history.

## Modes and output topology

- Circular `single` renders one record.
- Circular `grid` renders several records in one figure; a one-record grid is
  valid.
- Circular `batch` renders one output per record and has one resolved output
  target per item.
- Linear renders one ordered multi-record layout. Record groups and explicit
  comparison pairs affect topology and must survive session round trips.

The preview and download layers consume render results; they must not infer a
different grouping from the number of uploaded files.

## Assets, CSP, and privacy

Vue, Pyodide, styling, fonts, icons, DOMPurify, and export libraries are vendored
under `gbdraw/web/vendor/` or otherwise packaged locally. Keep the Content
Security Policy in `index.html` aligned with actual same-origin needs. Adding an
external host to the CSP is not a substitute for vendoring a dependency.

Do not log genome sequence, full uploaded file contents, or generated
comparison rows. Error messages should identify the failed stage and safe
resource label without exposing private data.

Generated SVG must pass through the shared sanitization profile before preview
or interactive editing. Event-handler attributes, scripts, foreign content,
and unsafe URL schemes remain forbidden.

## Changing a setting

Before adding or changing a control, trace its complete lifecycle:

1. Declare its reactive state and default in the appropriate state/setup owner.
2. Bind the visible control and its help text. Give each input an accessible
   name that matches its visible label; a placeholder or a title that changes
   with state is not a name (R12).
3. Project it once into the canonical request. The projection does not convert
   a value: a blank becomes `null`, a number stays that number, and anything
   else is a typed `INPUT_INVALID`. Generate does not rewrite the draft (R7).
4. Decode and validate it in the typed Python request layer if it is new.
5. Preserve it through session save/load when it is user-owned state.
6. Add focused state/request/session tests and a browser assertion when visual
   behavior changes.
7. Update Gallery tutorial captures when the control is part of a documented
   workflow.

For a dynamically available control, also identify its availability, fallback,
reconciliation, and live-edit owner under the invariants above. For a setting
the Result shows, choose its application path under R1: a live-edit owner, or
**Applies on Generate** with the note beside the control.

Do not add a diagram argv builder. Project each control once at the canonical
request boundary instead of duplicating it in configuration and session modules.

## Gallery ownership

`tools/prepare_interactive_gallery_assets.py::EXAMPLES` is the public example
inventory. `gbdraw/web/gallery/examples.json`, session artifacts, source/example
SVGs, thumbnails, and `artifact-manifest.json` are generated projections.
`tools/refresh_gallery_sessions.py` is the supported unfiltered owner command;
hosted builders verify these checked-in bytes and do not regenerate them.

`gbdraw/web/gallery/index.html` reads `gbdraw/web/gallery/examples.json`.
Tutorial instructions live under
`gbdraw/web/gallery/tutorials/`; screenshots and thumbnails live under
`gbdraw/web/gallery/media/` and `thumbnails/`.

Use `tools/capture_gallery_tutorial_screenshots.py` for declarative tutorial
captures. Data-dependent captures must load the example's own session, declare
the expected app state, and prove that the controls or data identity named by
the instruction are visible in the final crop.

Read `.claude/skills/web-gallery-screenshot-maintenance/SKILL.md` before editing
Gallery tutorials or screenshots.

## Local build and verification

```bash
gbdraw gui
python tools/prepare_browser_wheel.py
python tools/prepare_browser_wheel.py --refresh-cache-bust
python -m build
pytest tests/ -v -m "not slow"
npm ci
node --test tests/web/typed-boundaries.test.mjs
```

Refresh the cache-bust token only when preparing a deployable bundle.

When browser verification matters, check both Playwright installations:

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

The JavaScript specs under `tests/web/` require Node's `@playwright/test`. If it
is unavailable, use Python Playwright for focused browser checks. In an agent
sandbox, Chromium may fail with `sandbox_host_linux.cc ... Operation not
permitted`; rerun the same check with the required sandbox escalation.

For offline packaging, test that the wheel version matches `pyproject.toml`,
all runtime URLs are local, CSP allows the required worker/runtime behavior, and
the app reaches a ready state without external network access.

## Debugging principles

- Treat request validation failures as boundary errors, not reasons to bypass
  the typed request path.
- Throw a failure the user can fix with `diagnosticError` in JavaScript or
  `diagnostic=` in Python, carrying a code and a bounded context (R6).
- Do not add a classifier that matches message text.
- Keep worker errors structured enough for the main thread to distinguish
  initialization, resource transfer, validation, LOSAT, render, and export
  failures.
- Revoke object URLs and terminate superseded workers so repeated generation
  does not leak browser memory. A Cancel with no active Worker request keeps the
  warm Worker (`cancelDiagramGeneration` in `services/diagram-generation.js`).
- Preserve cancellation and stale-result guards when changing async code.
- Verify both Circular and Linear modes after shared state or request changes.

Related project documentation:

- `CLAUDE.md`
- `docs/REFERENCE/typed-requests.md`
- `docs/REFERENCE/session-and-request-compatibility.md`
- `docs/internal/SESSION_FORMAT_HISTORY.md`
- `docs/REFERENCE/interactive-svg-and-semantic-hooks.md`
