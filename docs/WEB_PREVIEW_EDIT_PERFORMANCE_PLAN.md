# Web Preview Edit Performance Plan

## Goal

Make preview editing feel proportional to the edit size.

Changing 1 feature should update 1 feature. Changing 71 selected features should
perform one bulk edit over 71 features. Undo and redo should restore the same
small edit intent, not reload the whole session.

The target is to remove the current hot path where feature edits repeatedly
clone, serialize, store, and reinitialize the full SVG preview.

## Problem Summary

The current implementation treats the serialized SVG result as both the preview
view and part of the editable application state.

This makes small edits expensive:

- Bulk feature actions call `history.runUndoable()` in
  `gbdraw/web/js/app/app-setup.js`.
- `history.runUndoable()` captures a full before and after snapshot in
  `gbdraw/web/js/services/history.js`.
- `buildHistorySnapshot()` includes config, UI, files, results, features,
  editor state, orthogroup state, and run state in
  `gbdraw/web/js/services/history-snapshot.js`.
- `serializeResults()` serializes the current SVG from the live DOM in
  `gbdraw/web/js/services/config.js`.
- `serializeCleanSvg()` clones the full SVG, strips transient classes, and
  serializes it with `XMLSerializer` in
  `gbdraw/web/js/services/svg-serialization.js`.
- Updating `results[].content` retriggers `watch(svgContent)` in
  `gbdraw/web/js/app/watchers.js`, which redoes legend extraction, drag setup,
  diagram setup, and feature handler indexing.
- Feature visibility bulk edits currently apply one feature at a time; each
  per-feature preview update serializes the full SVG and writes `results`.

The result is an accidental O(selected features * full SVG size) path for some
bulk edits, plus full-session snapshot work before and after every undoable
action.

## Design Principles

- Keep the source of truth small. Store edit intent as structured data, not as a
  serialized SVG snapshot.
- Treat SVG DOM as a view/cache. It may be updated immediately for preview
  responsiveness, but it should not be the primary edit model.
- Make bulk actions bulk-native. Do not implement bulk by looping through
  single-feature actions with full persistence side effects.
- Keep history separate from session save/load. Undo/redo records should be
  compact edit commands; session export may still build a complete project file.
- Prefer explicit update calls over broad watchers. If a feature fill changes,
  update feature fills; do not run the full SVG replacement pipeline.
- Add code only when it removes more complexity elsewhere. New helpers must
  replace duplicated paths or isolate a real responsibility.

## Non-Goals

- Do not redesign diagram generation or Pyodide execution.
- Do not change CLI output behavior.
- Do not remove session save/load.
- Do not make a second editing system beside the existing feature, legend,
  label, and layout editors. The plan should consolidate edit behavior, not
  fork it.

## Target Architecture

### 0. Active Preview Runtime Boundary

Split the current result into two explicit representations:

- Persisted result record: `results[].name` and `results[].content`, used for
  generation output, session save/load, export fallbacks, and remounting a
  stored result.
- Active preview runtime: the live SVG root currently mounted in the preview,
  per-root DOM indexes, a dirty flag, and serialization helpers.

`results[].content` must not be used as the edit event bus. A feature edit
updates the structured edit model and the active SVG DOM, then marks the active
runtime dirty. It does not rewrite `results[].content` just to trigger Vue
watchers.

Only the active result can be dirty. Hidden results are represented by their
stored `results[].content`; they are not edited while hidden. Result switching
therefore has a hard boundary: flush the active dirty runtime before changing
`selectedResultIndex`, then mount the next result from its stored SVG text.

The runtime can be small:

```js
{
  resultIndex,
  svg,
  dirty: false,
  dirtyReasons: new Set(),
  indexes: {
    features: null,
    legend: null,
    pairwiseMatches: null,
    orthogroupComparisons: null
  }
}
```

### 1. Structured Edit Model

The preview edit state should be the source of truth:

- `featureColorOverrides`
- `featureVisibilityRules`
- `featureStrokeOverrides`
- label text and visibility overrides
- legend entry edits
- canvas padding and layout offsets

Direct SVG DOM changes are a projection of that state for preview
responsiveness. They are not the canonical state for undo/redo or session save.

Rule-based edits must be represented by their canonical rule state, not only by
the list of currently affected SVG element IDs. For example, a visibility
command should restore `featureVisibilityRules` and the visibility override
cache, then update the affected DOM elements from that restored state.

### 2. Single-Stack Command History

History should remain one logical undo/redo stack, but entries may have different
storage strategies:

```js
{
  type: 'snapshot',
  label: 'Generate diagram',
  before,
  after
}
```

```js
{
  type: 'command',
  label: 'Change selected feature visibility',
  apply() {},
  revert() {},
  estimateBytes() {}
}
```

`undo()` and `redo()` dispatch by entry type. Snapshot entries continue to call
`applyHistorySnapshot()`. Converted preview edits call command `apply()` and
`revert()` and must not call `buildHistorySnapshot()` or
`applyHistorySnapshot()`.

Do not add a second history stack. A separate stack would make mixed operations
ambiguous, especially when a user generates a diagram, edits features, changes a
file input, and then undoes those actions in order.

Converted actions should enter through a dedicated API such as
`history.runUndoableCommand(label, buildCommand)`. They should not call
`history.runUndoable()`, because that API is snapshot-based by design.

### 3. Explicit Preview Mutation Layer

Introduce a small preview runtime/mutation layer that owns the active SVG DOM,
dirty state, and current DOM indexes:

- `mountResultSvg(resultIndex, svgText)`
- `markActiveResultDirty(reason)`
- `flushActiveResult({ stripTransient: true })`
- `applyFeatureFillChanges(changes)`
- `applyFeatureVisibilityChanges(changes)`
- `applyFeatureStrokeChanges(changes)`
- `applyLegendChanges(changes)`
- `invalidatePreviewIndexes(reason)`
- `getFeatureElements(featureId)`

This layer should not know about Pyodide, file uploads, or the session schema.
It may accept callbacks for reading/writing the active result content, but its
core responsibility is to update the current preview efficiently and record
whether serialized output is stale.

### 4. DOM Index Lifecycle

Build feature, legend, pairwise-match, and orthogroup comparison indexes when an
SVG root is mounted. Reuse them only while the operation is known not to change
the indexed structure.

Index invalidation must be explicit:

- SVG root remount or regeneration: invalidate all indexes and delegated handler
  state.
- Feature fill, visibility, and stroke edits: keep indexes.
- Selection class updates: keep indexes and touch only changed selected IDs.
- Label reflow or feature structure changes: invalidate feature and label-related
  indexes.
- Legend entry add/remove/reorder/reflow: invalidate legend indexes and any
  layout index that depends on legend geometry.
- SVG normalization that changes IDs or definitions: invalidate affected indexes
  before reuse.

Delegated event handlers should be installed per SVG root. Small edits should
not rebuild handlers. Structural edits should invalidate indexes, not silently
reuse stale caches.

### 5. Deferred Serialization

Serialize the full SVG only when a consumer requires SVG text:

- download/export, if the export path cannot use a live SVG clone directly;
- save session;
- active result switch;
- explicit regeneration result capture;
- optional idle-time cache refresh.

Before session export, call `flushActiveResult()` once, then let
`serializeResults()` read stored result content. Before result switching, flush
the active result, clear the active runtime, and mount the newly selected result
from stored content.

Normal Apply, Undo, and Redo should update structured state and current DOM only.

## Implementation Plan

### Phase 0: Establish The Active Preview Runtime

This phase creates the boundary that the rest of the work depends on.

- Add a small preview runtime module for the active SVG root, dirty flag, and DOM
  indexes.
- Add `mountResultSvg(resultIndex, svgText)` and use it when generated, imported,
  or selected result SVG text is mounted.
- Add `flushActiveResult()` that serializes the active SVG once and writes it
  back to `results[activeIndex].content` only when the active runtime is dirty.
- Replace direct result switching paths with an explicit `selectResult(index)`
  flow that flushes the old active result before changing the selected result.
  Keep a watcher only as a defensive fallback, not as the primary switch
  mechanism.
- Call `flushActiveResult()` before session export. Export paths that clone the
  live SVG may keep doing so, but any path that reads stored SVG text must flush
  first.
- Keep loading old sessions unchanged: restored `results[].content` remains the
  source used to mount a result.

Acceptance criteria:

- Feature edits can mark the active runtime dirty without rewriting
  `results[].content`.
- Saving a session after an unflushed edit includes the edit.
- Switching results preserves edits to the previous active result and mounts the
  next result from stored content.
- No hidden result can be dirty; dirty state is active-runtime-only.

### Phase 1: Stop Per-Feature Serialization In Bulk Visibility

Use the new runtime boundary for the smallest high-impact preview edit.

- Add a bulk visibility preview function that accepts all affected feature IDs
  and updates their SVG elements in one DOM pass.
- Update `setSelectedFeaturesVisibility()` so it computes all rule changes first,
  updates `featureVisibilityRules` once, rebuilds the visibility override cache
  once, updates the SVG once, marks the active runtime dirty once, and queues
  label reflow once.
- Remove the per-feature call path where `applyVisibilityPreviewBySvgId()`
  serializes and writes `results` for each selected feature.
- Keep existing single-feature visibility behavior by routing it through the
  same bulk helper with one ID.

Acceptance criteria:

- Hiding 71 selected features performs no SVG serialization during Apply.
- `watch(svgContent)` is not triggered once per selected feature.
- The visual result and label reflow behavior remain unchanged.

### Phase 2: Add Single-Stack Command History

- Extend `createHistoryManager()` so the undo and redo stacks can store both
  snapshot entries and command entries.
- Add a command API such as `runUndoableCommand(label, buildCommand)`.
- Keep `runUndoable()` as the snapshot path for unconverted operations.
- Start with selected feature visibility because it is the clearest performance
  problem.
- Store canonical before/after state in commands. For visibility, store the
  before and after `featureVisibilityRules` data plus the affected feature IDs
  used for efficient DOM updates.
- On apply/revert, restore canonical state, rebuild the relevant derived cache,
  update only affected DOM elements, and mark the active runtime dirty.
- Do not allow a converted action to store both a snapshot entry and a command
  entry.

Acceptance criteria:

- Converted feature visibility Apply, Undo, and Redo do not call
  `buildHistorySnapshot()`.
- Converted feature visibility Undo/Redo does not call `applyHistorySnapshot()`.
- Snapshot and command entries share one undo/redo order.
- Undo/redo restores `featureVisibilityRules`, the override cache, and visible
  SVG DOM.

### Phase 3: Convert Feature Color And Stroke Commands

- Convert selected feature color and stroke bulk actions to command entries.
- Store canonical before/after override state, not only DOM attributes.
- Apply/revert through the preview mutation layer so DOM updates and dirty
  marking are centralized.
- Keep single-feature color/stroke edits routed through the same command
  machinery where possible.

Acceptance criteria:

- Converted color and stroke Apply, Undo, and Redo do not serialize the full SVG.
- Undo/redo latency depends on changed feature count, not full SVG size or total
  session size.
- Legend state, feature state, and exported/saved SVG remain correct.

### Phase 4: Make SVG Initialization Explicit

- Rename or split the current broad `watch(svgContent)` responsibilities into
  explicit flows:
  - `mountGeneratedSvg()` or `mountResultSvg()` for new SVG roots;
  - `refreshEditedSvgParts(reason)` for small DOM edits that require a specific
    follow-up.
- Move feature handler indexing and delegated event setup to SVG mount only.
- Run legend extraction only when legend structure changes or a new SVG is
  mounted.
- Run diagram drag setup only when the SVG root or draggable layout elements
  change.
- Keep `watch(svgContent)` focused on root remounts. Small edits should not
  update `svgContent`.

Acceptance criteria:

- Feature fill, visibility, and stroke edits do not rebuild delegated feature
  handlers.
- Console timing for `attachSvgFeatureHandlers()` appears on generation or SVG
  mount, not on every converted preview edit.
- Preview pan, zoom, feature popup, hover, legend drag, layout edit, and search
  still work after edits.

### Phase 5: Centralize Index Invalidation And Selection Diffs

- Move feature, legend, pairwise-match, and orthogroup comparison index ownership
  into the preview runtime layer.
- Add explicit invalidation calls for structural operations, including label
  reflow, legend structure changes, SVG root replacement, and SVG ID
  normalization.
- Replace `syncFeatureSelectionClasses()` full-SVG walks with a diff-based
  update:
  - remove classes from IDs that left the selection;
  - add classes to IDs that entered the selection;
  - update the anchor class only for the old and new anchor IDs.
- Use the preview runtime feature index instead of
  `svg.querySelectorAll('[data-gbdraw-feature-id]')` for selection class
  updates.
- Keep full geometry scans only for rectangle and range selection because those
  operations genuinely need rendered positions.

Acceptance criteria:

- Toggling one selected feature touches only that feature's rendered elements
  and the old/new anchor elements.
- Fill, visibility, and stroke edits keep DOM indexes valid.
- Structural edits invalidate the affected indexes before later reuse.
- Rectangle/range selection behavior remains rendered-geometry based.
- Transient selection classes are still stripped from exported/saved SVG.

### Phase 6: Migrate Remaining Preview Edit Paths And Remove Replaced Paths

After the converted feature operations are stable, migrate the remaining preview
edit paths in risk order:

- label text and visibility overrides;
- legend cosmetic edits;
- legend structural edits;
- canvas padding and layout offsets;
- diagram and legend drag persistence.

Delete or simplify old paths as they are replaced:

- immediate `serializeCleanSvg()` calls from preview edit actions;
- bulk loops that call single-feature persistence side effects;
- snapshot-history usage for converted preview edits;
- watcher work that exists only because `results[].content` was used as the edit
  event bus.

This phase is required. Keeping both old and new paths would violate DRY and
make future feature edits slower to reason about.

## File-Level Touch Points

Expected primary files:

- `gbdraw/web/js/app/preview-runtime.js` (new, or equivalent small module)
- `gbdraw/web/js/services/history.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/export.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/svg-serialization.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/results.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/app/feature-editor/svg-actions.js`
- `gbdraw/web/js/app/feature-editor/color-actions.js`
- `gbdraw/web/js/app/feature-editor/visibility-actions.js`
- `gbdraw/web/js/app/feature-selection.js`
- `gbdraw/web/js/app/legend/*`
- `gbdraw/web/js/state.js`

Avoid adding a large framework. A small command-history extension and a small
preview-runtime/mutation module should be enough if existing responsibilities
are collapsed as they migrate.

## Test Plan

Add or update focused web tests:

- Active preview runtime edits do not rewrite `results[].content` until
  `flushActiveResult()`.
- Session save flushes the active dirty result exactly once and includes the
  latest live DOM edits.
- Result switching flushes the previous active result before mounting the next
  result.
- Bulk visibility over many selected features updates all targets and records
  one command history entry.
- Visibility commands restore canonical `featureVisibilityRules` state for
  feature, product/protein_id, and orthogroup scopes.
- Bulk color and stroke edits preserve legend/feature state and undo/redo
  correctly.
- Mixed snapshot and command history entries undo/redo in one chronological
  order.
- Undo/redo restores structured state and visible SVG DOM without rebuilding a
  full session snapshot for converted actions.
- Exported SVG strips transient selection/search classes.
- Index invalidation tests cover root remount, fill/visibility/stroke edits,
  selection diffs, label reflow, and legend structure changes.

Manual checks:

- Use a large linear comparison with thousands of features and many pairwise
  matches.
- Select tens or hundreds of features.
- Apply Hide, Show, color, and stroke changes.
- Undo and redo repeatedly.
- Confirm feature popup, hover summaries, pan/zoom, search, legend drag, layout
  edit, and export still work.

## Performance Checks

Use the existing console timing groups during migration.

Expected direction:

- Feature bulk Apply should not log repeated `watch(svgContent)` timings.
- Converted Undo/Redo should not run `buildHistorySnapshot()` or
  `applyHistorySnapshot()`.
- `attachSvgFeatureHandlers()` should run on generation/SVG mount, not on
  every converted edit.
- Full SVG serialization should appear only at export, session save, active
  result switch, regeneration capture, or optional idle-cache boundaries.

## Risks

- Some existing features may rely on `results[].content` being up to date
  immediately. Phase 0 must route stored-SVG consumers through
  `flushActiveResult()` and allow live-SVG consumers to clone the active DOM.
- Direct `selectedResultIndex` mutations can bypass the flush-before-switch
  boundary. Replace UI entry points with `selectResult(index)` and keep the
  watcher as a fallback only.
- Active-only dirty state is intentionally narrower than per-result dirty state.
  If future UI adds hidden-result editing, it will need a separate structured
  edit model per result instead of DOM-based dirty flushing.
- Command entries must store canonical before/after state. Storing only DOM
  element IDs is not enough for rule-based edits such as product/protein_id or
  orthogroup visibility.
- Legend edits can change SVG structure more than feature fill/visibility
  edits. Convert them after the simpler feature edit commands are stable.
- Label reflow may legitimately rebuild parts of the SVG. Keep it separate from
  simple visibility/fill/stroke updates.
- Mixed snapshot and command history during migration must stay in one stack.
  Keep the fallback explicit and remove snapshot usage as soon as converted
  operations cover the common preview edit paths.

## Definition Of Done

- Common preview edits are proportional to the number of changed items.
- `results[].content` is no longer used as the preview edit event bus.
- Active dirty result flushing is centralized and happens before session save,
  stored-SVG export paths, and result switching.
- Apply, Undo, and Redo for feature color, visibility, and stroke do not require
  full SVG serialization.
- Converted Undo/Redo does not restore full session snapshots.
- Converted command history and unconverted snapshot history share one undo/redo
  stack and one chronological order.
- DOM indexes are reused only when valid and invalidated explicitly after
  structural edits.
- Session save/load and exports remain correct.
- The code has fewer implicit side effects: SVG generation, preview mutation,
  history, and session serialization have separate responsibilities.
- Replaced snapshot and watcher-based paths are removed instead of left as
  duplicate maintenance burden.
