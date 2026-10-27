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
the semantic fields used to prepare those values. The same cache may also hold
the render-stage results that the last successful render used: derived data,
group contents with their serialized text, and feature-catalog parts, including
the catalog item of a Result. Key each of these by the complete inputs of the
stage that built it: resource and record identities, the settings it declares,
and every other argument. Every Generate still decodes and validates the
request, runs the whole render path, and writes a complete SVG document; it
writes it as ordered pieces: the root start tag, each top-level element, and
the end tag. A placed group is its start tag, its children text, and that
text's revision (its sha256). Joined in order, the pieces are the document. The
Generate rebuilds every part whose inputs changed. A Generate publishes the
decoded resources and parsed sources it built (they depend only on the
resource bytes and their key) when it ends, also when its render fails, unless
the request changed them; it publishes resolved records and render-stage
results only after the full render succeeds. The blastp protein extraction
helper and the Circular ring search-frame helper
(`convertRingTablesToSearchFrame`) also publish, when they end, the decoded
resources and parsed sources they built for the files they staged, unless they
changed them; any other helper lookup publishes nothing. A build that failed is never kept. Release every cache fill with the Worker; never reuse a final
reply.

The reply carries a piece's children text only when the request did not name
its revision among the segments the page holds (`heldSegments`). The page
sanitizes each children text once per revision, keeps the sanitized text on the
Result it commits, and joins the sanitized pieces in order. For the pieces
Python writes, which are whole SVG elements with escaped text and no comments,
CDATA or processing instructions, the joined text equals the sanitization of
the whole document; tests verify it.

Every Generate leaves the page holding a complete feature catalog equal to the
one a Generate with an empty cache returns. The reply carries that catalog
unless the request names the catalog revision the page holds for the drawing
and this Generate reused that same catalog build. A revision names one catalog
build in one Worker lifetime: a build id with the Worker's epoch, never a
content hash. (A segment revision is the opposite: only a content hash.) For a
matching revision the reply says the catalog is unchanged; the Worker does not
serialize, hash, or send it, and the page keeps the catalog object it admitted
for that revision without admitting it again. A changed catalog is always sent
and admitted. Session save, export, and History read the kept object as they
read a returned one.

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
| The drawings by ID and the shown drawing (settings, editor edits, derived values) | `js/state.js` (`drawingOf`, `drawingList`, `activeDrawingId`, `activeDrawing`, and `drawingOfMode` until B; R15) |
| The slices of a drawing's editor edits and settings, their only writers, and their install, capture, and reset for Session Load, History, and Reset | `js/services/drawing/` (R15) |
| Typed Vue reactivity for drawing members (`ref`, `reactive`, `computed`) | `js/utils/reactivity.js` (R14) |
| Registry of what a drawing holds in a Session, by mode | `gbdraw/web_support/mode_scoped_settings.py`, generated into `js/mode-scoped-settings.generated.js` by `tools/generate_mode_scoped_settings.py` |
| Session 46 mode slices and the split of a Session 27–44 draft into them; once the Session format holds a list of drawings, Python alone converts a Session 27–44, once at Load | `js/services/mode-scoped-migration.js` until then (Python: `split_draft_into_modes` in `gbdraw/session_io.py`) |
| Generate-button orchestration | `js/app/run-analysis.js` |
| Canonical request/session projection and equivalence | `js/services/session-request.js` |
| The displayed Result's text: marking its mounted root changed, the one write of the root into the Result, and the reader of the Results (R1) | `js/app/preview-runtime.js` (`commitActiveResultEdit`, `currentResults`) |
| When each setting reaches the Result: the application class and note of every settings control (R1) | `js/services/setting-application.js` (`SETTING_APPLICATION`; `SettingNote` in `js/components.js` renders the note) |
| Current active-config defaults, inventory, and validation | `js/services/session-active-config-contract.js` |
| Historical Gallery session migration; its Session 27–44 conversion holds until the Session format holds a list of drawings; then Python alone converts a Session 27–44, once at Load | `js/services/gallery-session-migration.js` until then |
| Per-feature edits and Feature placements by source identity (draft key, row checks, request rows, notices) | `js/services/feature-placement.js` |
| Rendered-ID feature edits of Sessions 44 and older, and Feature placement drafts of Sessions 41–44, onto mode-scoped identity rows; `hash=` annotation targets of Sessions 40–44 (R-7); drawn-hash `hash` color rules of Sessions 44 and older, renamed to source hashes; this conversion holds until the Session format holds a list of drawings; then Python alone converts a Session 27–44, once at Load | `js/services/feature-edit-migration.js` until then (Python twins in `gbdraw/session_io.py`: `migrate_session_feature_edits`, `migrate_session_annotation_targets`, `migrate_session_hash_rules`; the CLI replay of a request also uses their mapper `source_hash_selector_value`) |
| Check of scoped Feature placements in the current writer's config (`validateScopedFeaturePlacements`; not a migration, so it stays when the migrations go) | `js/services/feature-edit-migration.js` (imported by `services/session-active-config-contract.js`) |
| Gallery publication preparation, finalization, and readiness | `js/services/gallery-session-publication.js` |
| Save/load coordination | `js/services/config.js` |
| Load's reading of saved color and stroke edits (a Legend row's `colorOverrides` and `strokeOverrides`, a feature's stroke overrides) for Session Load and the Session 46 Legend reading | `js/services/style-overrides.js` (Python: `_legend_color_overrides`, `_stroke_override_map`, `_stroke_width` in `gbdraw/session_io.py`) |
| History transactions and availability: one History per drawing and the project stream behind one router, with a shared byte budget and file retention | `js/services/history.js` (`createHistoryRouter`), `js/app/history-inputs.js`, `js/services/history-snapshot.js` |
| Error producer contract and user-facing wording | `js/utils/error-normalization.js` (`diagnosticError`, `normalizeUserFacingError`) |
| Numeric draft projection | `js/utils/optional-positive-number.js` (`projectOptionalNumber`) |
| GenBank header reading | `js/services/genbank-header.js` |
| Managed Depth slot reconciliation | `js/services/depth-track-state.js` (`reconcileManagedDepthSlots`), called by `changeCircularDepthSources` and `changeLinearDepthSources` |
| Editor intent shown on the displayed Result: the edit transition, the port, and the preparation of a History step or a Result display (R3) | `js/app/app-setup.js` (`editEditorIntent`, `showEditorIntent`, `projectMountedEditorIntent`); it compiles with `compilePlanBundle` in `js/app/candidate-render.js` and reconciles with `reconcileMountedResult` in `js/services/svg-result-ingestion.js` |
| The drawing switch, the shown drawing, and the artifact of each drawing that is not shown | `js/app/app-setup.js` (`selectDrawing`); the artifacts wait by drawing ID in the slot store (`createArtifactSlotStore` in `js/services/artifact-slot.js`, created by `js/services/history-snapshot.js`), which the switch and Session Load change only through its named operations; the one writer of `activeDrawingId` is `installShownDrawing` in `js/services/config.js` |
| Render-worker client and lifecycle | `js/services/diagram-generation.js` |
| Pyodide typed rendering | `js/workers/diagram-generation-worker.js` |
| The Worker's package recipe, and the Pyodide snapshot of that runtime (capture, checks, restore; served only by a test server) | `js/workers/pyodide-snapshot.js` (builder: `tools/build_pyodide_snapshot.mjs`) |
| LOSAT dispatch and workers | `js/services/losat.js`, `js/workers/losat-*` |
| SVG/PNG/PDF downloads | `js/services/export.js` |
| SVG Result admission and sanitization | `js/services/svg-result-ingestion.js` |
| Feature-editor entry point | `js/app/feature-editor.js` |
| Feature-editor helpers | `js/app/feature-editor/` |
| Right-side editor drawer state and transitions | `js/services/right-drawer-state.js`; the drawer controller is `js/app/right-drawer.js` |
| Settings sections: the catalog (key, modes, title, icon), the section index entries, and the section a click on a drawn element opens | `js/services/settings-sections.js` (`SETTINGS_SECTIONS`, `settingsEntriesOf`, `drawnElementSection`); `SettingsSectionIndex` in `js/components.js` renders the index, template in `index.html` |
| Settings pane display (the stacked settings or the icon rail), the entry docked beside the rail, and the section in the popup beside the Result with its placement and availability | `js/services/settings-pane-state.js` (`showSettingsRail`, `showSettingsStacked`, `dockSettingsSection`, `floatSettingsSection`, `moveSettingsPopup`, `returnSettingsSection`); the controller is `js/app/settings-pane.js` (`createSettingsPaneController`; `settingsPopupAvailable` decides when the popup can open and closes it when it no longer can) |
| Browser-local preferences (the settings pane display and the UI language, PD-OI-106; the browser remembers them, not the Session) | `js/services/browser-preferences.js` (`readSettingsPaneDisplay`, `writeSettingsPaneDisplay`, `readUiLanguage`, `writeUiLanguage`); the header Language menu, the stored language at boot, and the translation feedback link beside the menu (the shown language and `GBDRAW_VERSION` of `js/config.js`) are `js/app/ui-language.js` (`chooseUiLanguage`, `showStoredUiLanguage`, `translationFeedbackUrl`) |
| Web UI text translation (`t`, `tn`, `tx`, `<i18n-text>`, locale catalogs, the shown UI locale); user-facing error wording stays with `normalizeUserFacingError` (R6) and is translated from its output, never by matching message text | `js/utils/i18n.js`, `js/utils/i18n-catalogs/` |
| The slot renderer names and ids by which the Web app finds drawn groups of a Result (scale, Linear length bar, dinucleotide and depth tracks, plot title; Python writes them) | `js/services/drawn-groups.js` (`TICK_GROUPS`, `DINUCLEOTIDE_CONTENT_GROUPS`, `DINUCLEOTIDE_SKEW_GROUPS`, `DEPTH_GROUPS`, `LENGTH_BAR_ID`, `PLOT_TITLE_ID`, `findTrackGroups`, `trackGroupSelector`) |
| Fold state of the Linear File cards | `js/services/linear-sources.js` (`setLinearSourceFolded`, `restoreFoldedLinearSources`); the set is `linearFoldedSources` in `js/state.js` |
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
If a mounted target exists, update it synchronously before optional geometry
reflow: a paint edit (R2) through the port `showEditorIntent` (R3), a label
edit through the label owner (`app/feature-editor/label-actions.js`). The edit
is in the current Result when the action returns: its text follows the mounted
root on demand (R1). A stroke action, a Legend row color,
and a Legend row deletion are each one transition, `editEditorIntent` in
`app/app-setup.js`: the action writes the intent of the drawing it resolved,
through that drawing's paint slice (R15), the port shows it when it changed,
and both are one step in that drawing's History (R11). The other paint edits
are History steps of their own that show through the same port (R2). If
canonical output requires creating or removing geometry and no mounted target exists,
queue the owning automatic rerender in the same action and replace the Result
when it completes. That rerender renders
the committed Session plus the current editor-intent tables
(`projectCommittedEditorIntent` in `services/session-request.js`). Of the draft
settings it reads only the values an action changed in the live rows of
`SETTING_APPLICATION` (`services/setting-application.js`), and no other draft
setting (R1). A rerender failure must keep any direct edit already
applied and report the failure. Metadata-only edits
with no static SVG target, such as Similarity group names and descriptions,
update canonical state only.

## Design rules R1-R15

Each rule names one owner and the guard that enforces it. Cite the rule id in PR
descriptions and tests. Change the code and its guard together. A count
baseline is a ceiling: the observed count stays at or below it, lowering it is
optional (its test names what may be lowered), and raising it is an
authority-only change.

### R1: Three writers of the Result

The current Result changes only through:

- (a) the editor intent, compiled by the Generate compiler
  (`compilePlanBundle` in `app/candidate-render.js`) and applied by the
  executor in `services/svg-result-ingestion.js`: at Generate admission onto
  the new Result, and for a live paint edit, a History step, or a Result
  display onto the displayed Result through the port `showEditorIntent` and
  `reconcileMountedResult` (R3). The label elements belong to the preview
  binder, which applies the label operations in its label step, and a live
  label edit updates them through the label owner
  (`app/feature-editor/label-actions.js`);
- (b) a composition edit that keeps parity with the Python renderer;
- (c) an automatic rerender from the committed Session and declared projection
  fields, plus the draft values an action changed in the live rows of
  `SETTING_APPLICATION` (`services/setting-application.js`). It commits those
  values to the committed Session (`services/config.js`) only when it succeeds
  and the committed Session has not changed since it started. If the committed
  Session changed, it is stale and runs again from the current committed
  Session; if it fails, it commits nothing and does not run again.

The displayed Result's text (`content`) follows its mounted root on demand. An
edit by (a) or (b) changes the mounted root and marks it changed
(`commitActiveResultEdit` in `app/preview-runtime.js`), without serializing it.
`currentResults` in `app/preview-runtime.js` returns the Results; when the
mounted root changed since it was last written, it first serializes the root
once and writes the text into its Result, if the Result is still the one the
root shows. Code that reads a Result's text, or holds the Results to restore
them later (Save, a History checkpoint, the capture of a generated artifact or
of a drawing's artifact slot), reads them through `currentResults`; the release
of the mounted Result (a Result switch, a drawing switch, a Session Load, a
Generate) writes it the same way first. Export reads the mounted root. No other
code writes a mounted root into a Result.

A setting that none of these applies is **Applies on Generate**: say so beside
the control, write only the draft, and leave the Result unchanged. Never patch
the Result from a draft value outside the declared live rows or with a helper's
SVG fragment. The one exception is the drawn-track hide and show that
`PD-OI-104` keeps (`web.layout.track-toggle-application`): `applyTrackVisibility`
in `app/svg-styles.js` hides or shows a track the Result draws at once (rows
`track-visibility` and `depth-sources`); a track the Result lacks and the
layout apply on Generate. A Generate reply's
pieces are not a helper fragment: they are the Generate compiler's input and
are admitted whole into one new Result. A new Python helper that
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
  of a setting that applies on Generate leaves the Result unchanged until
  Generate);
- `tests/web/preview-runtime.test.mjs`: an edit serializes nothing; the first
  `currentResults` after edits serializes once and writes once, and a second
  serializes nothing; a release writes the departing Result's edits only while
  it is still current; an unchanged root keeps the Result's bytes.

### R2: Lifetime of editor intent

Editor overrides are created, pruned, or deleted only by an explicit Reset or
Import, Undo or Redo, a Session replacement, the owner reconcile inside a
successful Generate (`pruneUnmatchedFeatureOverrides` in
`services/feature-visibility.js`), **Remove N unmatched feature edits**
(`removeUnresolvedFeatureEdits`), or **Delete settings** of a record whose
drawing is turned off, which deletes that record's per-feature edits, Feature
placements, and the annotations that target that record. Result selection,
mount, record selection (including turning a record's drawing on or off), a
drawing switch, hiding, and reflow never touch them. A source-replacing Generate keeps
the rows of records that are turned off. The owner reconcile prunes only in a
source-replacing Generate; in every successful Generate it moves the rows of a
record whose request key changed with its view (a Circular record drawn alone or
on the Multi-Record Canvas, a Linear File card that draws all of its records or
one of them) to its new key, and keeps the rows of a record that the request does
not draw under that record's identity (`recordKeyMoves` in
`services/record-draw-selection.js`): its per-feature edits, Feature placements,
feature strokes, and annotation targets. It also brings the rows a drawn
record has under its other keys (its key drawn alone or on the Multi-Record
Canvas, its identity, or the record key of the Linear card it was split from)
into its key, field by field: the row moved from its previous key, or else the
row already under its key, wins. The `hash` color rules of those records move
with the specific rules the Generate commits (`moveInstanceRuleRecords` in
`services/feature-selector.js`); of two rules that then name one feature
instance, the record's own stays and the other is dropped. Per-feature edits
(`featureOverrides`) and Feature placements live in the drawing they were made
in (R15) and are keyed by
`JSON.stringify([recordKey, biologicalFeatureId])`. Two drawings can use the
same record key for the same feature (`record-1` in a Gallery or Python
Session); each drawing holds its own row, so a request, the live projection,
and every reconcile read only the rows of their own drawing and records, and
every other drawing keeps its rows. Only `services/feature-placement.js` builds a draft
key. A selected-feature annotation target (`featureIdentity`) belongs to its
drawing the same way, and `annotationOptionsPayload` in
`services/annotation-state.js` gives a request only the targets of its records.
A row or target of a Session 41–44 that names its mode (`scope`) keeps that
field only until the Session 46 split moves it into that mode's drawing.

The paint intent of a drawing is its paint slice (`services/drawing/paint.js`,
R15): its Legend rows (`legendEntries`, with the `deletedLegendEntries`,
`dormantLegendEntries`, and `addedLegendCaptions` lists) and its color and
stroke edits (`legendColorOverrides`, `legendStrokeOverrides`,
`featureColorOverrides`, `featureStrokeOverrides`). Only that module writes
them. The edits that call its writers are History steps that show what they
change through `showEditorIntent` (R3):

- a stroke action, a Legend row color, and a Legend row deletion, each
  wrapped by `editEditorIntent`;
- the rule commit (`commitOnce` in `app/feature-editor/rule-actions.js`): a
  feature fill, a Legend row fill that rules draw, and the Legend rows, names,
  and colors the rules draw;
- a Change Scope choice;
- until each drawing binds its own auxiliary file cards, the Specific rules file card Remove for a drawing that is not shown, which
  calls `moveLegendRecords` (`takeBackRuleRows` in
  `app/feature-editor/rule-actions.js`, which Delete settings also uses for the
  `hash` color rows of a turned-off record's feature instances);
- a Legend rename, the Restore of deleted Legend rows, and the Sort and Move
  of Legend rows (a rename of a row that no rule or feature names writes the
  row's name in the step the Legend name input or the popup's rename action
  holds open, its own `rename` route);
- an edit of another owner that retires the paint edits of what it removed,
  in that edit's step: track data that no longer names a Legend row, a drawn
  Depth, GC content, or GC skew track that the draft hides or shows again (its
  Legend rows leave the Legend list, a renamed one waiting in
  `dormantLegendEntries`, and return to it; PD-OI-104), and **Delete
  settings** of a record whose drawing is turned off.

Undo and Redo, a Session load or replacement, Reset Settings, Generate, and the
binding of a Result install or rebuild it outside these steps (Python's Legend
rows and colors, the edits a source-replacing Generate prunes, and the feature
strokes a Generate moves to a record's new key). The Legend
name input also writes an unchanged copy of the list outside a step, so that it
shows the row's name again after a refused or cancelled rename. No other code
calls the paint slice's writers; a new caller, or a caller that moves to
another route, is an authority-only change.

Guards: `tests/web/feature-visibility-rule-writers.test.mjs`
(`PAINT_INTENT_WRITERS`: each function that writes the paint intent, directly,
by calling a slice function with the paint or visibility-rules slice handle,
or by calling a writer it imports from that slice's module whatever handle it
passes, by its name, by module, with its route; a slice reader is one only by
the `SLICE_READERS` list in `tests/web/helpers/drawing-state.mjs`, which the
isolation table in `tests/web/per-mode-drawings.test.mjs` proves changes no
drawing, and a new reader is an authority-only change; a write through a row
read from the list, a slice writer handed on as a value, one reached through another
module's re-export, and `p.SLICE.method(` on a namespace import are not seen); `tests/web/non-edit-state-preservation.playwright.spec.js` with
`tests/web/non-edit-state-diff.test.mjs` (user-owned state is identical before
and after non-edit operations); the disjoint-Result binding test in
`tests/web/feature-label-visual-unit.test.mjs`;
`tests/web/mode-scoped-feature-edits.test.mjs` (no draft row reaches a request,
reconcile, or TSV load of the other mode, and no other module builds a draft
key) and `tests/web/mode-scoped-feature-edits.playwright.spec.js`;
`tests/web/annotations.test.mjs` and the OV-21 test in
`tests/web/annotation-feature-identity.playwright.spec.js` (annotation targets).

### R3: One projection per domain

A live paint edit, a History apply, and a Result display show the paint and
visibility domains of the editor intent through one port,
`showEditorIntent({ domains })` in `app/app-setup.js`. It
compiles the operations of the shown domains from the active drawing with the
Generate compiler (`compilePlanBundle` in `app/candidate-render.js`, running
only the stages those domains need) and reconciles the displayed Result with
`reconcileMountedResult` in `services/svg-result-ingestion.js`, the executor of
Generate admission: an attribute of a shown domain that no operation sets
returns to the value Python drew, so Undo, Redo, and a Result display apply no
inverse. A live edit reaches the port through its action (R2). A History apply and a Result display first prepare the
domains whose intent changed with `projectMountedEditorIntent`, which calls the
domain preparations `projectPaletteAndRules` in `app/app-setup.js` (the
palette's track colors and the specific rules' matches) and
`projectFeatureVisibility` in `app/feature-editor/visibility-actions.js`
(feature visibility and the labels of the features it hides or shows), and then
shows them in one call; the labels follow in the preview binder's label step.
The other domain functions are
`orderLegendEntries` in `services/legend-svg.js` (legend order) and
`resolveFeatureDrawn` in `services/feature-visibility.js` (whether a feature is
drawn, as Python's `should_render_feature` answers; the Features list and
Search features read it through `listFeatureRows`). A displayed batch Result whose shared legend entries already follow
the legend order keeps its order (`orderLegendEntries` with `keepFollowed`),
so the entries only that Result draws keep their places; a Result last shown
with another order also receives the default order. A History legend step made
on another batch Result writes the displayed Result's Legend list from the rows
it shows and the step's shared Legend intent (`adoptRestoredLegend` in
`app/legend/entry-actions.js`, with the step's other side, B19), so no Result
gains or loses an entry that only one Result draws. The displayed
population comes from the mounted Result's committed metadata
(`renderedFeatureIdentities`), not from a second selection ref. Per-feature
edits reach a Result through its catalog features' identities
(`resultRenderedFeatures` in `services/feature-catalog.js`), so every copy of
a feature in every Result shows the same edit. A dialog's
reactive object holds display values only, never a copy of an owner's data.

One projection means one call shape: a domain's projection function is called
from outside its owner in one form, not re-implemented or re-sequenced per
caller. The count of distinct call shapes per domain stays at or below its
baseline (R13): lowering an entry is optional, and adding or raising one is an
authority-only change.

A live edit shows on the Result what the next Generate draws (PD-OI-066). A
parity case or a step of the random walk compares a live path with Python's
Generate from the same draft; a case that compares one Web path with another
through the same executor is not kept. The pull-request random walk
(`tests/web/live-generate-random-walk.playwright.spec.js`, over
`tests/web/helpers/live-random-walk.mjs`) draws every edit kind of its list on
every run, with a fixed seed, and after each step checks that the Result equals
a fresh Generate and that the drawing the step did not edit keeps its slices.
A change to a live-edit path is done when the walk passes, and, for a path or
an editor state the walk does not draw, when the case for its edit kind in
`tests/web/live-generate-parity.playwright.spec.js` passes with the popup
opened and not, and with Auto Reflow on and off when the path reads it
(`queueLabelReflow` without `force`). A matrix case leaves that spec only
when the walk draws its kind in the states the case pinned, which the walk
spec asserts. Each path that shows the editor intent (a live edit, Undo, Redo,
a Result display, a drawing switch, Reset Settings, Session Load, and the
automatic rerender) has rows in the work allowlist (`WORK_ALLOWLIST` in
`tests/web/helpers/live-work-allowlist.mjs`), each pinning the compile stages
the path runs, the Worker requests it sends, and its compile count; the work
guard in the parity spec runs the rows, and the walk names the step kind that
draws each row's path, or none. A difference the Result may show until
Generate is listed with its decision in
`tests/web/contracts/live-generate-parity-allowed.json`, which the cases and
the walk read; adding or changing an entry is an authority-only change. A
mismatch a case or the walk finds is a bug: it stays an expected failure with
its finding id until the fix makes it pass.

Guards: `tests/web/gui-audit-20260930-editor.playwright.spec.js` (a Result shows
the edits made on another Result, also after Undo, Save, and Load),
`tests/web/feature-visibility-actions.test.mjs`,
`tests/web/feature-color-actions.test.mjs`,
`tests/web/live-generate-random-walk.playwright.spec.js` with
`tests/web/helpers/live-random-walk.mjs` (every edit kind of the walk on every
run, a step kind for each work allowlist row, and each matrix case the walk
replaced drawn in that case's states),
`tests/web/live-generate-parity.playwright.spec.js` (the work guard over
`tests/web/helpers/live-work-allowlist.mjs`, the matrix cases the walk does not
replace, and named regression cases) with
`tests/web/live-generate-parity-legend.playwright.spec.js` and
`tests/web/live-generate-parity-track-data.playwright.spec.js` (named
regression cases of Legend edits, track-data changes, and drawing switches): a
live edit equals the next Generate; and the `projection-shapes` baseline in
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
coordinates. Raw rows and **Save Raw LOSAT TSV** stay in F, except for Circular
similarity rings: their LOSAT database is the whole reference file (uncropped),
so their raw rows and **Save Raw LOSAT TSV** are that search's output in Src,
and the Python owner `ring_rows_in_search_frame` converts them to F before they
become ring tables. Every comparison
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

- `tests/web/error-producer-coverage.test.mjs`: the total of unclassified
  `throw new Error` sites in all modules under `js/`, workers included, and
  each message-table size stay at or below their baselines.
- `tests/test_web_error_producer_coverage.py`: the total of `ValidationError`
  sites of the package that neither carry a diagnostic nor classify and each
  Python adapter table size stay at or below their baselines, and the Python
  vocabulary is a subset of the Web wording owner.
- `tests/web/error-normalization.test.mjs`: table-driven producer cases and
  code, reason, and context-key parity.

### R7: Python's typed layer owns value checks

The Web projection never converts a value: a blank is `null`, a finite number is
that number, and anything else is a typed `INPUT_INVALID`
(`projectOptionalNumber`). A check that must run before Python (comparison
thresholds before LOSAT) evaluates the domains generated from
`gbdraw/mode_profiles.py` into `js/mode-profiles.generated.js`
(`resolveComparisonThresholds` in `js/mode-profiles.js`). Generate never rewrites
the draft settings; the editor edits it moves to a record's new key are R2's.

Guards: `tests/web/option-input-integrity.test.mjs` (every projected numeric
draft field either rejects `NaN`, `1e-50x`, `-5`, `0`, and `12.5` with a typed
diagnostic or sends them literally; the count of draft assignments on the
Generate path in `app/run-analysis.js` stays at or below its baseline, where
lowering the baseline is optional and raising it is an authority-only change)
and
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

Logic that Python and the Web each implement (a twin: a reader, migration,
normalizer, layout, or constant with one copy on each side) has a row in
`tools/twins.json`: its concern, both owners, the shared vector, the reason,
and the condition that retires it. The registry only shrinks: a row leaves
when its retirement condition is met or one side is deleted. Adding a row,
or adding an owner to a row (a further copy of a registered twin), is an
authority-only change. Removing or renaming an owner, or editing a row's
reason, retirement condition, or a vector that is not part of its id, is
not; a `vector:` row's vector is part of its id. A new `vector:` or
`function:` row goes in before its twin, in that authority-only pull request,
with `"pending": true`: the guard does not require a pending row's twin to be
detected or its owners and vector to exist. The pull request that adds the
twin removes `pending`, which is not an authority-only change; a pending row
whose twin is detected fails the guard, and a pending row whose twin is
abandoned is removed.

Guards: `tests/test_comparison_tables.py` (`test_outfmt_table_columns_have_one_reader`
scans the Python and Web JavaScript sources for a second mapping onto
`COMPARISON_COLUMNS`, and the shared vectors run through every entry point),
`tests/web/genbank-header-reader-owner.test.mjs` (no Web module other than
`services/genbank-header.js` matches a GenBank header keyword), and
`tests/web/twin-registry.test.mjs` (the detector `twin.unregistered.v1` in
`tools/web-twin-detectors.mjs` finds a JavaScript test and a Python test that
read the same vector file, and a watch-list function name defined both in
`gbdraw/` Python and in `js/services/`; each such twin has a row that is not
pending, and the owners and vector of every row that is not pending exist).

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
`app/feature-editor/visibility-actions.js`, which writes the drawing's rules
through its visibility-rules slice (R15), projects them, and follows with the
labels as one History step.

No watcher writes a drawing. A transition that changes a drawing calls the
reaction of that change explicitly: the palette draft sync after a color edit,
the canvas padding after a padding edit, the label override rows after a bulk
label edit, the Circular multi-record layout defaults after the multi-record
toggle, the queued palette when Instant Preview turns on
(`setPaletteInstantPreview` in `app/results.js`), and the annotation record
bindings at the end of a record discovery and a drawing switch
(`rebindAnnotationRecords` in `app/app-setup.js`). An auxiliary table import
(default colors, color table, qualifier priority, whitelist, blacklist) is a
transition of its file input (`importAuxiliaryTable` in `app/watchers.js`)
that writes through the slices of each drawing it changes. Until each drawing
binds its own auxiliary file cards, a card's file is shared by every drawing:
choosing, replacing, or removing it takes the earlier file's rows out of every
drawing, and the new file's rows go into the shown drawing. Once each drawing binds its own
auxiliary file cards, a card shows the shown drawing's file, a file chosen in it goes into the
shown drawing, and its Remove changes only that drawing's rows (PD-OI-090).

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

- Each drawing has its own History, and the project has one project stream,
  behind one router (`createHistoryRouter` in `services/history.js`) that the
  composition root creates. A transaction, `undo`, `redo`, and availability go
  to the shown drawing's History; a step closes in the History that opened it,
  also after a drawing switch. Undo and Redo walk only the shown drawing's
  steps and the project steps, one step per operation, and never undo a step of
  another drawing; the Undo or Redo of a drawing's step never switches the shown drawing
  (PD-OI-090, OIC-030): Undo takes the newer top step of the two, Redo the
  older. A new project step clears every Redo; a new drawing step clears its
  own and the project's.
- A step that changes a value of a drawing belongs whole to that drawing; the
  Session-level values it also changes are part of its entry and are restored
  with it where they are still as it left them. A project step is one of:
  an operation on the drawing list (adding, duplicating, renaming, moving, or
  deleting a drawing); Reset all drawings (until the Reset dialog offers Reset
  this drawing, Reset Settings resets all drawings); until each drawing binds
  its own auxiliary file cards, a card's import, Remove, and Retry (R10); once
  a genome input file can be replaced for every drawing that uses it, that
  replacement (PD-OI-090); or a step that changes only Session-level values
  (the Session title, Auto Reflow, PNG DPI, Instant Preview, the LOSAT
  execution and threads, the Feature panel tab, and the drafts; R15-23
  Q-DR-8). Such a step records only
  the Session-level changes made since it opened. A Session-level change made
  outside any step marks nothing unsaved and rides with the next step; the
  drawing switch's own transient reset is no change. Each operation is one step. Reset
  this drawing is one step of the shown drawing, Reset all drawings is the only step that
  resets more than the shown drawing, and the Reset dialog's Cancel records nothing
  (OIC-030). Deleting a drawing that has a Result or edits asks first, and the last
  drawing cannot be deleted. A new drawing, a Duplicate included, starts with an empty
  History (PD-OI-090).
- Turning Instant Preview on is one project step with a part for each
  drawing: it applies every drawing's queued palette, and its Undo brings each
  queue back (R15-23 Q-DR-5, OV-1160). Turning it off changes only the flag.
- A step captures and restores the artifact slot values of a drawing that is
  not shown (the applied palette, the alignment Reset receipt) in that
  drawing's own slot, through the one port of `services/history-snapshot.js`
  keyed by drawing ID (`readSlotValues`, `writeArtifactValues`). It never
  writes them to the shown Result (OV-1168, OV-1195).
- The Histories share one byte budget, which evicts the steps of the drawings
  that are not shown first, then project steps, then the shown drawing's,
  oldest first within each; and one file retention, which keeps a file while
  any History or binding names it.
- A transaction belongs to one owner, a control or a gesture. Starting another
  owner's transaction first settles the open one. A legend or diagram drag passes
  a per-gesture owner.
- `settlePendingIntent` in `services/history.js` is the one settlement of the
  open intent. `begin`, `beginCheckpoint`, `beginArtifactReplacement`,
  `runUndoableCommand`, `undo`, and `redo` share it.
- A discrete control begins in the capture phase of the event that commits its
  value: `change` for a checkbox, radio, file input, or select, `click` for a
  button. A text control begins on focus. The pointer, the label, the keyboard,
  and a file picker opened by another button record the same one step
  (`app/history-inputs.js`).
- A file input whose change ends with asynchronous work is
  `data-history-managed` and owns one step that commits after that work: a
  `file-uploader` `afterChange` import, and **Add Seq**
  (`addCircularConservationComparisonFile`), whose step holds the ring record
  label read by the Python reader. Add Seq is a step of its drawing
  (`runUndoable`); an auxiliary file card's import and its Retry are one
  project step (`runProjectStep`) while the card's file is shared by every
  drawing (R10); once each drawing binds its own cards, they are steps of the shown drawing.
- An undoable action started while another action still runs on the open step
  applies at once and names that step (`runUndoable`); nothing waits. What an
  action writes after another caller committed its step is a step of its own
  under the action's name (`commitIntent`). That step is not recorded when the
  action ends while a checkpoint is open (it captures the write), while an
  artifact replacement is open, or while an Undo or Redo is being applied, nor
  when, since its step was committed, a baseline was set, an artifact
  replacement committed or was declined (`shouldCommit`), or a checkpoint ended
  without committing.
- A call the running action makes itself keeps the step: one made
  synchronously inside its `fn`, or one that passes `joinsOpenStep` (the Label
  Not Shown choice, and the rule table commit of a popup, Features panel, or
  dialog color action). A `joinsOpenStep` call commits under the step name of
  the action that started on the step with its label, and under the step's
  current name when no action has that label.
- Undo, Redo, and the header buttons share one reactive availability
  (`historyAvailability`) of the shown drawing's History, busy while an artifact replacement or checkpoint is
  open or an action owns the open intent. Undo and Redo, and their enabled state, also read one check of the
  step they would take, the top of their stack (`topStepAvailability`).
- A drawing's History intent holds, by reference (`buildIntentFilesData` in
  `services/history-snapshot.js`), the file bindings of the input slots its
  drawing uses and the auxiliary table files its own steps imported, the BLAST
  rows and comparison sequences of a LOSAT-cache replay included: its ring rows and the
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
  Result. A drawing switch (`selectDrawing`, R15) is navigation as well: it
  settles the departing drawing's open step and records no step.

Guards:

- `tests/web/history-inputs.test.mjs` and
  `tests/web/history-inputs.playwright.spec.js`: control kind, input means, focus
  state, and last-in-first-out Undo; each changed control adds exactly one step.
- `tests/web/history-generated-authority.playwright.spec.js`: checkpoint Undo and
  Redo, a drawing round trip that records no step, and Undo while Generate runs.
- `tests/web/history-per-drawing.test.mjs`: Undo and Redo walk the shown
  drawing's steps and the project steps and never change the other drawing (a
  hash of each drawing per step); a switch settles the departing drawing's open
  step in that drawing and records none; Reset Settings, an Instant Preview
  step are one project step each, and while a card's file is shared by every
  drawing (R10) so is a card import's Retry; a Session-level
  edit is a project step; the project step limit and the slot port of a
  drawing that is not shown.
- `tests/web/multi-result-edit-matrix.playwright.spec.js`: Undo and Redo of
  drags on one batch Result while another Result is displayed (B17), a
  Result switch after a legend sort that records no step (B18), Undo and Redo
  of a legend sort made on another Result (B19), and Sort by default reaching
  another Result (B20).
- `tests/web/history-config-restore.test.mjs`: catalog identity through Undo,
  Redo, and Session rollback, and unset settings that stay unset through intent
  and checkpoint Undo and Redo (also
  `tests/web/session-draft-authority.test.mjs` for the rollback).
- `tests/web/non-edit-state-preservation.playwright.spec.js` (G-H rows): Undo^n
  Redo^n through mixed edits and a drawing round trip, and Reset Settings Undo and
  Redo, leave user-owned state unchanged.
- `tests/web/session-operation-consistency.test.mjs`: Undo availability during
  Generate.
- `tests/web/history.test.mjs` (the action cases) and the Label TSV case in
  `tests/web/gui-audit-regressions.playwright.spec.js`: an edit made while an
  action runs is its own step under its own name, and Undo then Redo keeps both
  writes; the cases where the later write is in no step.
- `tests/web/specific-color-actions.test.mjs` and
  `tests/web/feature-label-visual-unit.test.mjs`: the step name of each
  `joinsOpenStep` call site through the real History manager.

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

Modules sit in layers, lowest first: leaves (`utils/`, `contracts/`,
`config.js`, `web-ux-profile.js`, `mode-profiles.js`,
`mode-profiles.generated.js`, `mode-scoped-settings.generated.js`), then
state-free `services/` (among them the drawing slices in `services/drawing/`:
a slice module imports only leaves and other slice modules), then
`state.js`, then the two state-bound services (`services/config.js`,
`services/reset.js`), the only services that import `state.js`, then owner
modules (`app/*.js` and `app/*/` factories), then the composition roots
(`app/app-setup.js`, `app/feature-editor.js`, `app/legend.js`,
`app/legend-layout.js`), then `app.js`, `components.js`, and the template. A
module imports its own layer and lower ones; a state-free module never imports
`state.js`, and an owner module never imports a composition root. `contracts/`
is a leaf: like every leaf its modules import only leaves.

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
- `state` is frozen (`Object.freeze` in `state.js`): no module adds or
  replaces a `state` member, so a value an owner needs from a higher layer
  arrives as a port or on an existing data path (Result metadata, the admitted
  catalog), never as a function assigned into `state`. A drawing's slice
  members are written only by their slice module (R15).
- A heavy derived value (`rulePreparation`) is triggered from one place per
  flow: the composition root's projection, the Generate compiler, or the owner
  of its input. An owner that only reads the result awaits `pending` rather
  than triggering it. A trigger site includes a call of a producer port (a
  producer method handed to an owner as a function, such as
  `rulePreparation.run` or `rulePreparation.runDrawn`) in the module that
  receives it.
- Each projection domain has one projection function and one call shape
  outside its owner (R3).
- The identity and membership of a Legend row or a drawn feature come from its
  identity (Python's Legend row key, the feature ids, `recordKey`), never from
  a fill or stroke that another owner writes or from a value read off the
  Result (its fill, stroke, Legend key attribute, transform, or text). The
  detector counts these joins (`paint-join`, `shown-join`); a join on a
  caption is a report-only inventory.

The baseline of the current violations is recorded in
`tests/web/owner-graph-baseline.test.mjs`: the count of a subject stays at or
below its entry, and a subject without an entry has none. Lowering or removing
an entry is optional (the test names what may be lowered). Adding or raising an
entry is an authority-only change (`tools/web-design-rule-guards.json`,
WEB_CHANGE_POLICY.md "Design-rule co-change"). The detectors are
`tools/web-owner-graph-detectors.mjs` (injection edges, forward-reference
closures `owner-graph.forward-closure.v2`, whole-object ports,
projection call shapes, heavy derived trigger sites
`heavy-derived.trigger-site.v2`, identity read from a display value
`owner-graph.identity-from-display.v1` (a join on a fill or stroke, and a join
on the values a Result shows), and import direction
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
imports from a higher layer),
`tests/web/owner-graph-detectors.test.mjs` (the detectors),
`tests/web/owner-graph-replay.test.mjs` (the detectors report the couplings
that past merges introduced; it reads the revision history and runs in the
`web-change-budget` job), `tests/web/drawing-context.test.mjs` (`state` is
frozen; tsc refuses a write to a member, R14), and
`tests/web/drawing-slices.test.mjs` (a slice module imports only leaves and
other slice modules).

### R14: Typed boundaries

Modules under `js/` except `workers/`, and `workers/pyodide-snapshot.js`, are
checked by `tsc --noEmit` with JSDoc types (`tests/web/types/tsconfig.json`); there is no build step and no
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
  catalog) has one typedef of what its writer writes, for the current format
  only. The Session and the canonical render request are declared once, in
  Python. `tools/generate_web_contracts.py` writes their typedefs, and the
  constants of the persisted contracts (among them `SESSION_VERSION`,
  `CANONICAL_REQUEST_SCHEMA`, and `FEATURE_CATALOG_SCHEMA`), to the modules
  `js/contracts/*.generated.js`. `js/contracts/` holds only the modules that
  tool writes. They are never edited by hand, and no other module declares a
  typedef or a constant that one of them declares. A hand-written Web copy of
  a Python constant or schema literal that the Web does not import from the
  generated modules is a `constant:` row of `tools/twins.json` (R9) until it
  does. The feature catalog's
  typedef is in its writer, `services/feature-catalog.js`.
- A reader of an unvalidated document takes `unknown`, or the generated input
  type whose values are `unknown` (`SessionInput`, and `SessionEnvelope` once
  `format` and `version` are read). The validator of each part returns that
  part's type. `GbdrawSession` is the writer's type; no reader returns it. A
  reader that still takes `Record<string, any>` is counted in
  `EXPLICIT_ANY_BASELINE`.
- Types are JSDoc in the module that owns them, imported with
  `/** @import { T } from './x.js' */`, never a `.d.ts` or a types module under
  `gbdraw/web/` other than the generated modules of `js/contracts/`. A type
  import follows the layers of R13. The globals of vendored scripts and the
  test hooks are declared once, as `any`, in
  `tests/web/types/web-globals.d.ts`.
- `@ts-ignore`, `@ts-expect-error`, and `@ts-nocheck` are not used: fix the
  type, cast with JSDoc, or fix the code.
- `noImplicitAny` is on in the guard config: a parameter or variable
  declares its type, or the collection a callback iterates is typed.
- Explicit `any` in JSDoc (`any`, `*`, or a type that resolves to `any`, such
  as an alias of it) is type debt. `EXPLICIT_ANY_BASELINE` in
  `tests/web/explicit-any-ratchet.test.mjs` records it per module: a module
  stays at or below its entry, and new code declares a type or uses
  `unknown`. A module without an entry, including a new module, has none.
  Lowering or removing an entry is optional. Adding or raising an entry is an
  authority-only change.
- Distinct identity values (a Python key, a display caption, a swatch color, a
  feature id, a drawing ID) take distinct named types, not one shared `string`
  alias. A drawing's identity is its `DrawingId`; its mode is its `mode`
  member and is never parsed from the ID.
- `state.js` and `services/drawing/` create a drawing's members with
  `utils/reactivity.js`, which types Vue's `ref`, `reactive`, and `computed`.
  Other modules see a slice's members as a readonly view (`DeepReadonly`), so a
  write outside the slice module does not type-check. The view is a type only,
  not a Vue `readonly` proxy: the vendored production build drops such a write
  silently.

Guards: `tests/web/typed-boundaries.test.mjs` (the checked set, the compiler
run, declared factory parameters, no suppression, type-import direction, no
second declaration of a generated typedef or constant) and
`tests/web/explicit-any-ratchet.test.mjs` (explicit `any` per module). They
run the `typescript` version pinned in `package.json` and fail until `npm ci`
has installed it. `tests/test_web_contracts_generated.py` checks that the
generated modules equal the tool's output and that `js/contracts/` holds
nothing else.

### R15: Drawing context and mode-scoped settings

A project has one or more drawings, held by ID in the drawing store of
`state.js`, a `ReadonlyMap` that `createDrawingState` fills and only
`state.drawingList()`, `state.drawingOf(id)`, and, until the Session format holds
a list of drawings, `state.drawingOfMode(mode)` read. `state.activeDrawingId` names the shown
drawing, and `state.activeDrawing()` returns it. A drawing's ID (`DrawingId`)
comes from the kind it is created as (`circular`, `linear`, `linear-2`, …) and never
changes, and its `mode` is fixed at creation; once the Session format holds a list of
drawings, several drawings can have one mode (PD-OI-090). A drawing holds its settings (`form`, `adv`, the LOSAT
search fields, the comparison plan, record display, the layout slot), its
editor edits (per-feature edits and placements, colors and color rules, label
tables and filters, Legend edits, annotation sets, canvas padding), and the
values derived from them, under their former `state` names. Nothing else holds
them, and two drawings share no member object, so a setting or an edit made in
one drawing never reaches another (PD-OI-090).

- The members are grouped into slices, one module each in `services/drawing/`:
  paint, feature edits, color rules, visibility rules, label tables, palette,
  annotation sets, group names, canvas, and settings. Each member is in exactly
  one slice. A slice's writers take its handle, `drawing.slices.<name>`
  (`services/drawing/slice-handle.js`).
- A slice module is the only writer of its editor edits. Its writers are named
  for the edit they make, and other modules read the slice through the
  drawing's readonly view (R14).
- Each slice module also installs, captures, and resets its slice. Session
  Load, the History restore and capture, Save, and Reset Settings call those
  three functions and no other copy; Reset Settings resets every slice through
  its own `reset`.
- Settings stay written by their controls and owner transitions (R7, R10).
  Only their install, capture, and reset go through the settings slice.

- `state` keeps only what every drawing shares: the input files (once the Session format
  holds a list of drawings, only their content: each drawing binds the files, records, and
  regions it draws), the LOSAT
  caches, the LOSAT execution and thread settings, app settings (Auto Reflow,
  PNG DPI, palette Instant Preview), the Session title, the shown drawing's
  Result and artifacts (the artifacts of each drawing that is not shown wait
  in its slot in the slot store, `services/artifact-slot.js`), and transient
  flags. Every `state` key is in
  exactly one class in `tests/web/helpers/drawing-state.mjs`. A value that
  drawings should share is a Product decision: record it in the Product
  Contract before moving it out of the drawing.
- A service never resolves the active drawing: it takes the `drawing` it
  reads (Session Save and the legacy configuration Load take it as
  `shownDrawing` from their caller). An owner resolves `state.activeDrawing()`
  once per action, or the drawing its action is about by its `DrawingId`
  (`state.drawingOf`), and passes it down. No helper resolves it again. Until
  the Session format holds a list of drawings, a site about the drawing of one
  mode resolves it through `state.drawingOfMode(mode)`, the one lookup by mode,
  whose callers may only decrease. Generate resolves its
  drawing once, at the entry of `runAnalysis`; its request, validation, and
  commit use that drawing only. A reflow or rerender of a committed Result uses
  the drawing that Result belongs to, and installs through that drawing's
  artifact owner.
- The drawing switch is `selectDrawing` in `app/app-setup.js`. It and a
  Session replacement (Load, a new Session) show a drawing through the one
  writer of `activeDrawingId`, `installShownDrawing` in `services/config.js`.
  The switch is navigation: it writes no drawing value, copies or resets
  nothing, and records no History step (PD-OI-090, R10, R11).
- The switch waits while Generate, a label rerender, Session Save or Load, a
  History restore, an open History transaction, or the matching of a color
  rule table edit runs (OV-1100, PD-OI-090 scenario revision 2).
- Before it switches, it settles the departing drawing: its open step, its
  pending rule-pattern drafts, the feature selection, and the shown Result's
  projection. It then keeps that drawing's artifact in its slot, installs the
  arriving drawing's slot, and resets transient UI. Reset Settings and the
  Instant Preview step change the slot of a drawing that is not shown in that
  slot, and History restores it through its slot port (R11); none of them
  writes it to the shown Result.
- Each drawing has its own History (R11).
- Until the Session format holds a list of drawings, the store holds exactly
  the drawings `circular` and `linear` (read through `state.drawingList()`,
  `state.drawingOf(id)`, and `state.drawingOfMode(mode)`), and Save writes them as
  `modes.circular` and `modes.linear`.
- Once the Session format holds a list of drawings (Session 49), a new Session starts with
  one Circular drawing. Save Session writes, as one document, every drawing that has a
  file, a changed setting, an edit, or a Result, the drawing that is shown, and only the
  files a saved drawing uses. Load installs every drawing or changes nothing, shows each
  saved Result without Generate, and shows the drawing that was shown at Save (PD-OI-090).
- A Session 46 saves each drawing as `modes.circular` and `modes.linear`.
  `gbdraw/web_support/mode_scoped_settings.py` is the one registry of what a
  slice holds (domain, path, modes, row key, migration token);
  `tools/generate_mode_scoped_settings.py` writes it into
  `js/mode-scoped-settings.generated.js`, which is never edited by hand. Load
  rebuilds each drawing from the projection of its committed request,
  then its stored values, then its mode's defaults.
  Until the Session format holds a list of drawings, `services/mode-scoped-migration.js`
  splits a Session 27–44 draft into the two slices, and `split_draft_into_modes` in
  `gbdraw/session_io.py` is its Python twin; then Python alone converts a Session 27–44,
  once at Load. Once the Session format holds a list of drawings, a Session 46, 47, or 48
  moves to the list without Python, and Load makes the drawing of a mode that has no
  Result only when the Session shows that mode was used
  (its input files, a Save while it was shown, or a setting that differs from the default
  and from the other mode's). A Session 45 is not read.

Guards: `tests/web/drawing-context.test.mjs` (the key classes, the drawings that
hold the drawing keys, no drawing key read from `state`, no `activeDrawing` identifier in
`services/`, and Generate resolving it once);
`tests/web/drawing-slices.test.mjs` (each member in exactly one slice, no
module or template binding writing a slice member outside its slice,
`install(capture(x))` giving back `x` also after `reset`, and
`installShownDrawing` the one writer of `activeDrawingId`);
`tests/web/drawing-of-mode-cap.test.mjs` (the `drawingOfMode` calls within
`DRAWING_OF_MODE_CAP`, which may only fall, and only `ID_BY_MODE_OWNERS`
naming a drawing ID by a mode; raising the cap or adding an owner is an
authority-only change);
`tests/web/per-mode-drawings.test.mjs` (no
shared member object, and each mode's request carries its own drawing's
values); `tests/web/mode-scoped-session.test.mjs` (Save and Load of each drawing in
the current format, a Session 46 read into its drawings, and a Session 45 or a field
outside the registry rejected);
`tests/web/mode-split-vectors.test.mjs` with `tests/test_mode_split_vectors.py`
(the Web Load of a Session 27–44 and Python's split agree on shared vectors);
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

Guards: `tests/web/computation-ownership.test.mjs` and
`tests/web/computation-ownership.playwright.spec.js`.

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
- the layout-preference slot that a committed request sets in its
  drawing (`modes.<mode>.ui.layoutPreferences` in a Session 46).

`services/config.js` coordinates user-facing save and load actions. It must not
grow a parallel model of render fields. Compatibility migrations are reader
concerns: normalize supported old data once, then use the current model. Do not
write retired fields into new sessions.

`services/session-active-config-contract.js` is the DOM-free current writer
contract used by normal restore and Gallery publication. Historical Gallery
migration stays in `services/gallery-session-migration.js` until the Session format holds a
list of drawings, then Python alone converts a Session 27–44, once at Load; it does not clean up
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
- A region drawing, made with + New drawing from a feature, a whole record, or coordinates,
  draws the region of the same record of the same file in source coordinates, 1-based and
  inclusive; regions of several records make one Linear drawing with one region per record.
  A region that crosses the origin of a circular record is refused with the reason. The
  source drawing does not change (PD-OI-090).

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
the Result shows, choose its application path under R1 and record it in its
`SETTING_APPLICATION` row: a live-edit owner, a live row that the automatic
rerender (c) applies, or **Applies on Generate** with the note beside the
control.

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
python -m build
pytest tests/ -v -m "not slow"
npm ci
node --test tests/web/typed-boundaries.test.mjs
```

`tools/prepare_cloudflare_pages.py` writes the wheel's content-hash token into the hosted bundle.

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
