# Web Feature Multi-Selection Plan

## Goal

Add multi-feature selection to the Web preview so users can select individual
features, contiguous rendered feature ranges, and rectangular groups before
applying feature-editing actions.

The first complete implementation must let the selected features receive one
color and legend caption, one visibility setting, and one stroke setting through
feature-editing actions. Selection-only UI is not a complete milestone.

The feature must preserve the current default behavior:

- plain feature click opens the existing feature popup;
- hover summaries continue to work;
- preview pan, zoom, layout edit, legend drag, and search highlights remain
usable;
- exported SVG/PNG/PDF output must not contain transient selection styling.

## Confirmed Interaction Decisions

- `Ctrl` click on Windows/Linux, and `Cmd` click on macOS, toggles one feature
  in the selection.
- `Shift` click selects the range between the anchor feature and clicked
  feature.
- `Shift` click range selection is allowed only within the same record scope.
  If the anchor and clicked feature belong to different records, do not modify
  the selection.
- `Shift` drag creates a rectangular selection marquee and selects features
  intersecting the rectangle.
- `Shift` click range ordering follows the rendered left-to-right order on the
  screen, not raw genomic coordinate order.
- Selection is cleared on regeneration, result switch, mode change, SVG
  replacement, `Escape`, and preview background click.
- Transient selection classes must be stripped anywhere live SVG is serialized
  from the first phase that introduces those classes.

## Implementation Principles

- Keep the selection engine small and preview-scoped. Do not create a parallel
  editing system.
- Reuse existing feature color, legend synchronization, visibility, stroke, and
  history paths.
- Implement color/legend-caption, visibility, and stroke bulk editing before
  calling the feature complete.
- Prefer rendered SVG geometry already available in the browser over duplicating
  layout math.
- Add helpers only where they remove repeated logic or prevent transient preview
  state from leaking into persisted output.

## Current State

- Preview feature interactions are delegated from the root SVG in
  `gbdraw/web/js/app/feature-editor/svg-actions.js`.
- Feature SVG elements are discoverable through `FEATURE_SELECTOR` and stable
  identity normalization:
  - `data-gbdraw-feature-id`
  - fallback `id^="f"` paths/polygons/rects
- `buildFeatureElementIndex()` already indexes SVG elements by `svg_id`.
- `extractedFeatures` and `featuresBySvgId` already provide feature metadata
  needed for selection scope and editing:
  - `svg_id`
  - `id`
  - `fileIdx`
  - `record_idx`
  - `record_id`
  - `displayRecordId`
  - `start`
  - `end`
- Preview search already applies transient SVG classes in
  `gbdraw/web/js/app/feature-search/preview-svg.js`; feature selection should
  follow the same pattern.
- Preview panning starts from `gbdraw/web/js/app/ui.js`, but it already skips
  feature targets. `Shift` drag selection still needs explicit coordination so
  the gesture cannot be interpreted as pan or layout drag.

## UX Model

### Plain Click

Plain click keeps the current behavior:

- open the clicked feature popup;
- clear any transient multi-selection;
- set the clicked feature as the range anchor.

This keeps existing users from seeing persistent selections after ordinary
feature inspection.

### Background Click

Clicking empty preview space clears the selection and closes open popups. Clicks
on form controls, floating preview controls, feature elements, pairwise matches,
labels, legends, and layout-edit targets must not be treated as background
clicks.

### Ctrl/Cmd Click

Modifier click toggles one feature:

- if the feature is unselected, add it;
- if the feature is selected, remove it;
- close any open feature or pairwise match popup;
- update the range anchor to the clicked feature when it remains selected;
- do not open the feature popup.

Use `event.ctrlKey || event.metaKey` so macOS users get the expected `Cmd`
behavior.

### Shift Click

Shift click performs range selection:

- if no anchor exists, select only the clicked feature and set it as anchor;
- if anchor exists and both features have the same record scope, select all
  features between them in rendered left-to-right order;
- if the record scopes differ, leave the selection unchanged and show a compact
  status message such as `Range selection is limited to one record.`;
- do not open the feature popup.

Range selection should add to the existing selection when `Ctrl`/`Cmd` is also
held. With `Shift` alone, replace the current selection with the range.

### Shift Drag Rectangle

Shift drag starts a rectangular marquee selection when the pointer moves beyond
a small threshold, for example 4 px.

Behavior:

- `Shift` mouse down records the drag origin.
- A small movement under the threshold still behaves as `Shift` click.
- Once the threshold is crossed, draw a fixed-position marquee overlay inside
  the preview area.
- On pointer up, select all rendered feature elements whose screen bounding box
  intersects the marquee.
- `Ctrl`/`Cmd` plus `Shift` drag adds the intersecting features to the existing
  selection. `Shift` drag without `Ctrl`/`Cmd` replaces the selection.
- Do not include hidden features or elements with empty feature identity.
- Close popups and suppress hover summaries while the marquee is active.
- After a committed marquee drag, suppress the browser-generated follow-up
  `click` so it cannot also open a popup or run `Shift` click range selection.

The rectangle selection is based on screen coordinates using
`getBoundingClientRect()`. This avoids having to invert the SVG transform and
works with current pan and zoom.

## Record Scope

Create a stable helper for feature record scope:

```js
const getFeatureSelectionScope = (feature) => [
  feature.fileIdx ?? '',
  feature.record_idx ?? '',
  feature.displayRecordId || feature.record_id || ''
].join('::');
```

Range selection requires equal scope between anchor and target. This prevents
ambiguous cross-record "between" behavior.

Individual `Ctrl`/`Cmd` selection may include features from multiple records.
Rectangle selection may include features from multiple records because it is a
spatial selection, not an ordered range. Bulk edit actions must therefore work
on an arbitrary selected feature list.

## Feature Order

For `Shift` click range selection, order only rendered features in the same
record scope by their current screen position.

Suggested ordering:

```js
[
  featureScreenCenterX,
  featureScreenCenterY,
  String(feature.type || ''),
  String(feature.svg_id || '')
]
```

Use `extractedFeatures` as the source list, but filter to features whose
`svg_id` exists in the current SVG element index. This avoids selecting
features that are extracted but not rendered under the current settings.
Compute screen centers from the union of each feature's rendered SVG element
bounding boxes. This keeps reversed linear records, region extraction, pan, and
zoom aligned with what the user sees. Do this only when a range selection is
requested; do not maintain a continuously updated geometry cache.

## State Model

Add focused selection state in `gbdraw/web/js/state.js`:

```js
const selectedFeatureIds = ref(new Set());
const selectedFeatureAnchorId = ref('');
const featureSelectionStatus = ref('');
const featureSelectionSuppressNextClick = ref(false);
const featureSelectionDrag = reactive({
  active: false,
  committed: false,
  startX: 0,
  startY: 0,
  currentX: 0,
  currentY: 0,
  additive: false
});
```

Use a `Set` for identity operations, but always replace it with a new `Set`
when updating so Vue watchers and computed values observe the change.

Derived values:

- `selectedFeatureCount`
- `selectedFeatures`
- `hasFeatureSelection`

The selected feature list should resolve IDs through `featuresBySvgId`, falling
back to `extractedFeatures` only when needed.

Selection lifecycle:

- clear before starting a new generation and after a successful generation;
- clear when `selectedResultIndex` changes;
- clear when `mode` changes;
- clear when `svgContent` is replaced by a non-incremental SVG update;
- clear on preview background click;
- clear on `Escape` when no marquee is active and no popup owns the key event.

Pan, zoom, search, hover, and opening/closing the feature drawer must not clear
selection by themselves.

## Module Structure

Keep `svg-actions.js` from growing into a general selection engine by adding a
focused module:

- `gbdraw/web/js/app/feature-selection.js`

Responsibilities:

- normalize feature selection scope;
- sort rendered features for range selection;
- update selected feature ID sets;
- apply and strip transient SVG selection classes;
- manage rectangular selection pointer events and overlay state;
- clear selection on preview lifecycle events;
- suppress the first click after a committed marquee drag;
- expose small action functions to `app-setup.js` and `svg-actions.js`.

If this module grows too large, split SVG class helpers into:

- `gbdraw/web/js/app/feature-selection/preview-svg.js`

This mirrors the existing feature search structure.

## SVG Styling

Add transient CSS classes in `gbdraw/web/index.html`:

- `.gbdraw-feature-selected`
- `.gbdraw-feature-selection-anchor`
- `.gbdraw-feature-selection-candidate`
- `.feature-selection-marquee`
- `.feature-selection-status`

Use visible but restrained styling:

- selected feature: blue stroke, slightly stronger stroke width, no fill
  change;
- anchor feature: blue/purple accent or stronger shadow;
- marquee: translucent blue fill with dashed border;
- status: compact overlay near the search control or lower preview controls.

Selection styling must use classes only. Do not write permanent stroke/fill
attributes to selected feature elements.

Implement `stripFeatureSelectionClasses(svg)` at the same time as these classes.
Use it anywhere live SVG DOM is serialized, including instant preview updates,
SVG/PNG/PDF downloads, session data, and standalone interactive SVG generation.

## Event Integration

Update `createFeatureSvgActions()` in
`gbdraw/web/js/app/feature-editor/svg-actions.js`:

1. Detect modifier click before opening a feature popup.
2. Route `Ctrl`/`Cmd` click to `toggleFeatureSelection()`.
3. Route `Shift` click to `selectFeatureRange()`.
4. Keep plain click on existing `openFeatureEditorForFeature()`.
5. Suppress hover summary while `featureSelectionDrag.active` is true.
6. At the start of the delegated click handler, consume
   `featureSelectionSuppressNextClick` with `preventDefault()` and
   `stopPropagation()` when it is set.
7. Treat non-feature, non-pairwise preview clicks as background clicks that clear
   selection.

Add pointer handling for rectangle selection at the SVG or preview container
level:

1. `pointerdown` with `Shift` records the origin and captures the pointer.
2. `pointermove` beyond threshold marks the gesture as marquee selection.
3. `pointerup` commits the rectangle.
4. `pointercancel`, `Escape`, and `blur` cancel the rectangle.

Use capture-phase listeners only where needed to prevent pan from receiving the
same `Shift` drag gesture.

## Pan And Layout Edit Coordination

Update `gbdraw/web/js/app/ui.js`:

- if `event.shiftKey` and SVG content exists, do not start pan;
- keep normal pan behavior unchanged for non-Shift drags.

Update layout drag guards if needed:

- layout edit should ignore `Shift` drag gestures so selection can win;
- legend drag remains plain drag only.

The intended priority is:

1. form controls and floating UI controls;
2. feature selection modifier gestures;
3. feature/pairwise popup click handling;
4. layout edit drag;
5. pan.

## Selection Toolbar

Add a compact preview overlay shown only when `selectedFeatureCount > 0`:

- selected count;
- Clear button;
- Open first selected feature button;
- Apply color/caption entry point using the existing color and legend editing
  path;
- Apply visibility entry point using the existing visibility rule path;
- Apply stroke entry point using the existing stroke editing path.

The first complete implementation must include bulk color/legend-caption,
visibility, and stroke application. Keep the toolbar compact; do not add a
second full feature editor.

## Bulk Edit Integration

Existing editor actions already contain paths for feature color, legend,
visibility, and stroke updates. Reuse those paths instead of creating a second
bulk editing mechanism.

Target behavior:

- selected features can receive one color and legend caption together;
- selected features can have visibility set together;
- selected features can have stroke set together.

Implementation approach:

1. Add a bulk action that receives `selectedFeatures`.
2. For color and legend caption, route through existing multi-feature
   synchronization helpers where possible.
3. For visibility and stroke, route through the existing rule/action paths
   rather than adding a second editing model.
4. Keep generated rules hash-based per selected `svg_id`.
5. Preserve undo/history semantics by wrapping one bulk edit in one history
   transaction.
6. Clear selection only after a successful destructive or broad operation if
   keeping the selection would be misleading.

Do not add broad regex rules from selected features automatically. Selection is
an explicit per-feature operation, so the generated rule should remain
feature-hash based.

Avoid duplicating the color editing flow. If current helpers are private inside
`color-actions.js`, expose the smallest action needed for selected feature
color/caption updates rather than copying rule-generation logic.
Apply the same rule to visibility and stroke: expose the smallest reusable action
needed for selected features.

## Search Highlight Interaction

Search classes and selection classes can coexist.

Rules:

- search dimming should not make selected features unreadable;
- active search match and selected feature should both be visible;
- clearing search must not clear feature selection;
- clearing feature selection must not clear search.

If CSS specificity conflicts appear, give selected features enough specificity
to keep their stroke visible while respecting search dimming for non-selected
features.

## Export And Persistence

Selection state is preview-only.

Requirements:

- session save should not persist `selectedFeatureIds` unless a later workflow
  explicitly needs resumable selection;
- SVG serialization for `results.value` must not capture selection classes;
- download/export helpers must strip selection classes from live SVG clones;
- generated standalone interactive SVG should not include current Web preview
  selection state.

Implement a helper similar to search cleanup:

```js
stripFeatureSelectionClasses(svg)
```

Call this helper anywhere live SVG DOM is serialized for persistent output.
This is not a late cleanup task; it ships with the first selection styling
phase.

## Accessibility And Keyboard

Add lightweight keyboard support:

- `Escape` clears active marquee if dragging;
- `Escape` clears selection when no marquee is active and no popup owns the
  key event;
- Clear button has an accessible label;
- status text should be visible in the DOM, not only console output.

Do not add complex keyboard list navigation in the first implementation.

## Implementation Phases

### Phase 1: Selection State And Classes

Files:

- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/feature-selection.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/services/export.js`
- live-SVG serialization paths such as
  `gbdraw/web/js/app/feature-editor/svg-actions.js`
- lifecycle wiring paths such as `gbdraw/web/js/app/run-analysis.js`,
  `gbdraw/web/js/app/watchers.js`, and `gbdraw/web/js/app/app-setup.js`

Tasks:

1. Add selection state and derived computed values.
2. Add helper functions for scope, ordering, Set replacement, and class
   toggling.
3. Add `stripFeatureSelectionClasses(svg)` and call it before live SVG
   serialization.
4. Add selected/anchor/marquee/status CSS.
5. Add lifecycle clear actions for regeneration, result switch, mode change, SVG
   replacement, preview background click, and `Escape`.
6. Add tests for pure helper functions if a suitable JS test harness exists;
   otherwise keep helpers simple and verify through browser checks.

Acceptance:

- applying selection classes to a current SVG visually highlights selected
  features;
- clearing selection removes all selection classes;
- search classes and selection classes can coexist;
- instant preview serialization and download serialization strip selection
  classes.

### Phase 2: Ctrl/Cmd Click And Shift Click

Files:

- `gbdraw/web/js/app/feature-editor/svg-actions.js`
- `gbdraw/web/js/app/feature-selection.js`
- `gbdraw/web/js/app/app-setup.js`

Tasks:

1. Route modifier clicks before the existing popup open path.
2. Implement individual toggle selection.
3. Implement same-record range selection using rendered left-to-right order.
4. Add non-blocking status for rejected cross-record range selections.
5. Keep plain click popup behavior unchanged.

Acceptance:

- `Ctrl`/`Cmd` click toggles one feature without opening a popup;
- `Shift` click selects a same-record range in screen left-to-right order;
- `Shift` click across records does not change the selection;
- plain click still opens the feature popup and clears prior selection.

### Phase 3: Shift Drag Rectangle

Files:

- `gbdraw/web/js/app/feature-selection.js`
- `gbdraw/web/js/app/ui.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`

Tasks:

1. Add pointer event handling for `Shift` drag.
2. Add marquee overlay rendering.
3. Select intersecting feature bounding boxes on commit.
4. Suppress the browser-generated follow-up click after a committed drag.
5. Add cancellation on `Escape`, pointer cancel, and window blur.
6. Coordinate with pan and layout drag guards.

Acceptance:

- `Shift` drag shows a rectangle and selects intersecting features;
- `Shift` drag does not pan the preview;
- `Shift` plus `Ctrl`/`Cmd` drag adds to the current selection;
- releasing a marquee drag does not open a feature popup or run Shift-click
  range selection;
- no stale marquee remains after pointer cancel or blur.

### Phase 4: Selection Toolbar And Bulk Color/Caption

Files:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/feature-selection.js`
- `gbdraw/web/js/app/feature-editor/color-actions.js`

Tasks:

1. Add a compact selection toolbar overlay.
2. Show selected count.
3. Add Clear action.
4. Add Open first selected feature action using existing
   `openFeatureEditorForFeature()`.
5. Add one bulk color/legend-caption action for `selectedFeatures`.
6. Wrap the bulk edit in one undoable history transaction.

Acceptance:

- users can see how many features are selected;
- users can clear selection without hunting for a keyboard shortcut;
- opening one selected feature reuses the existing popup;
- selected features can be assigned one color/caption in one operation;
- generated specific rules remain per-feature hash rules;
- existing single-feature color popup behavior is unchanged.

### Phase 5: Bulk Visibility/Stroke And Regression Checks

Files:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/feature-editor/visibility-actions.js`
- `gbdraw/web/js/app/legend/stroke-actions.js`
- `gbdraw/web/js/app/feature-selection.js`

Tasks:

1. Expose bulk visibility through the current visibility rule path.
2. Expose bulk stroke through the current stroke action path.
3. Keep each bulk operation as one undoable history transaction.
4. Add targeted packaging or browser checks where feasible.

Acceptance:

- selected features can receive one visibility setting in one operation;
- selected features can receive one stroke setting in one operation;
- bulk visibility and stroke actions reuse existing rule paths and do not duplicate
  logic;
- downloaded SVG/PNG/PDF does not contain selection styling;
- session save/load does not restore transient selection;
- existing web packaging tests still pass.

## Manual Verification Checklist

- Generate a linear diagram with multiple records.
- Plain click a feature and confirm the popup opens.
- `Ctrl` click two features and confirm both are highlighted and no popup opens.
- `Cmd` click behavior should be checked on macOS or by simulating `metaKey`.
- `Shift` click within one record selects the rendered left-to-right range.
- `Shift` click from one record to another leaves selection unchanged and shows
  status.
- `Shift` drag selects features intersecting the rectangle.
- `Shift` drag does not pan the preview.
- Releasing a `Shift` drag does not open a popup or trigger range selection.
- Click empty preview space and confirm selection clears.
- Switch result, regenerate, change mode, and confirm selection clears.
- Search for a feature while a selection is active; both visual states remain
  understandable.
- Apply one color/caption to selected features and confirm one undo restores the
  previous state.
- Apply one visibility setting to selected features and confirm one undo restores
  the previous state.
- Apply one stroke setting to selected features and confirm one undo restores the
  previous state.
- Export SVG and confirm selection classes are absent.
- Save and load session and confirm selection is not restored.

## Risks

- `getBoundingClientRect()` can be expensive on very large SVGs. Mitigate by
  computing intersections only on pointer up, not every pointer move.
- Search dimming and selected stroke styling can conflict. Keep CSS scoped and
  test both active states together.
- Serializing live SVG after selection could leak transient classes. Centralize
  cleanup helpers and call them from export and instant-preview persistence
  paths.
- `Shift` is used by browser text selection and some OS/browser shortcuts. The
  gesture is limited to the SVG preview area to reduce interference.

## Non-Goals

- Persisting selected feature sets in sessions.
- Adding polygon/lasso selection.
- Adding keyboard range navigation.
- Selecting labels, legends, pairwise ribbons, or non-feature tracks.
- Creating broad regex rules automatically from a selected set.
