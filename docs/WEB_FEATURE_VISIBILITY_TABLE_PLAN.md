# Web Feature Visibility Table Plan

## Goal

Add first-class Web UI support for `--feature_visibility_table` in the left
Features panel, and make the Editor tab write to the same rule set.

The Web app already supports per-feature visibility from the feature popup by
generating hash-based rows for `/web_feature_visibility_table.tsv`. This plan
turns that hidden mechanism into an editable rule table without creating a
second visibility system.

## Principles

- **Debt control:** every added line must either expose required behavior,
  remove duplicate state, or prevent a real inconsistency.
- **KISS:** keep one rule format: the existing TSV columns
  `record_id`, `feature_type`, `qualifier`, `value`, `action`.
- **YAGNI:** do not add a new policy engine, grouped presets, import wizards, or
  bulk rule templates until users need them.
- **DRY:** serialize, import, export, session-save, and Editor updates through
  the same rule helpers.
- **Single priority model:** the visible `featureVisibilityRules` array order is
  the serialization order and the Python matching priority. Do not apply hidden
  source-based reordering.
- **SOLID:** separate responsibilities:
  - UI components edit rule rows.
  - feature-editor actions convert one clicked feature to one hash rule.
  - run-analysis serializes the current rules into CLI TSV input.
  - Python remains the source of truth for final matching semantics.

## Current State

- CLI/Python accept `--feature_visibility_table`.
- The table columns are fixed:
  `record_id`, `feature_type`, `qualifier`, `value`, `action`.
- Actions normalize to `show`, `hide`, or `suppress`.
- First matching rule wins.
- The feature popup currently stores per-feature visibility in
  `featureVisibilityOverrides` and serializes those overrides as hash rows:

```text
*    *    hash    ^<feature_svg_id>$    show|hide|suppress
```

## Proposed Data Model

Add one canonical Web state array. This array is the source of truth for Web
visibility state:

```js
featureVisibilityRules = [
  {
    id,
    source,       // "manual" | "editor" | "file"
    featureId,    // svg_id for editor/hash rules, otherwise ""
    label,        // optional display-only feature label
    recordId,
    featureType,
    qualifier,
    value,
    action        // "show" | "hide" | "suppress"
  }
]
```

The array order is significant. It is the UI order, the generated TSV order, and
therefore the Python "first matching rule wins" order.

Compatibility layer:

- Keep accepting old saved sessions with `featureVisibilityOverrides`.
- During session load, convert those overrides into `source: "editor"` hash
  rules.
- Convert `featureVisibilityOverrides` to a derived compatibility cache
  immediately. It must never be directly edited, serialized, saved to new
  sessions, or used by run-analysis.
- If existing call sites still need the object shape during the same change,
  rebuild it from `featureVisibilityRules` after rule changes. The cache may
  represent editor/hash feature rows only; broad manual rules that cannot be
  represented as `{svg_id: mode}` stay only in `featureVisibilityRules`.
- When both `featureVisibilityRules` and old `featureVisibilityOverrides` are
  present in a session, load `featureVisibilityRules` and ignore the old map.

## Left Panel UI

Place a **Feature Visibility** subsection inside the existing **Features** card.
This is feature rendering behavior, not color styling, so it should not live in
the Colors card.

Controls:

- compact rule list with Up, Down, Delete buttons
- row editor fields:
  - `record_id`: text input, default `*`
  - `feature_type`: text input with datalist suggestions from `featureKeys`,
    extracted feature types, and `*`; arbitrary values must be allowed for GFF3
    and imported TSV compatibility
  - `qualifier`: text input with datalist suggestions for common values
    `product`, `gene`, `locus_tag`, `hash`, `location`, `record_location`
  - `value`: regex text input
  - `action`: select `Show`, `Hide`, `Suppress`
- Add rule button
- Download TSV button
- optional Load TSV control if it can reuse the existing file-import pattern
  with little code

Display behavior:

- The visible row order is the effective matching order.
- Editor-created hash rows should show the feature label/location when known.
- Hash rows remain editable as normal TSV rules, but the UI should make their
  feature-specific nature visible.
- Do not implement live Python-regex validation. Basic empty-field checks are
  enough; Python will report authoritative regex errors on Generate.
- Left-panel manual and imported rule edits apply on the next Generate. Do not
  attempt broad live preview in JavaScript.

## Editor Tab Synchronization

The Editor tab's Feature visibility dropdown should update
`featureVisibilityRules` directly.

Mapping:

- `Default` removes the editor hash rule for the clicked feature.
- `On` creates or updates an editor hash rule with `action: "show"`.
- `Off` creates or updates an editor hash rule with `action: "hide"`.
- `Suppress` creates or updates an editor hash rule with `action: "suppress"`.
- Existing editor hash rules keep their current array position when updated.
- New editor hash rules are inserted at the top of `featureVisibilityRules` so
  exact clicked-feature overrides win by default, but the user can move them
  with the same Up/Down controls as any other row.

Rule shape for Editor writes:

```js
{
  source: "editor",
  featureId: feat.svg_id,
  label: bestAvailableFeatureLabel,
  recordId: "*",
  featureType: "*",
  qualifier: "hash",
  value: `^${escapeRegexLiteral(feat.svg_id)}$`,
  action
}
```

Ordering:

- `featureVisibilityRules` order is the only priority model.
- Serialization must not regroup rows by `source`.
- Manual row order remains user-controlled.
- File-imported rows preserve file order relative to their insertion point.

This keeps "this exact feature" overrides stronger than broad rules such as
"show all CDS" by default, while still making the override order explicit and
user-controlled.

## Serialization

Create small helpers, likely near `gbdraw/web/js/app/file-imports.js` or a new
focused module under `gbdraw/web/js/app/feature-visibility/`:

- `normalizeFeatureVisibilityAction(value)`
- `normalizeFeatureVisibilityRule(raw)`
- `parseFeatureVisibilityRules(text)`
- `serializeFeatureVisibilityRules(rules)`
- `upsertEditorFeatureVisibilityRule(rules, feat, action)`
- `removeEditorFeatureVisibilityRule(rules, featureId)`
- `getEditorFeatureVisibilityMode(rules, featureId)`
- `buildFeatureVisibilityOverrideCache(rules)`

Keep these helpers free of Vue-specific dependencies so they are easy to test.

Parser behavior should intentionally match Python semantics:

- skip blank lines and `#` comments
- skip a header row
- reject rows with missing required fields
- reject rows with extra columns
- accept action aliases supported by Python:
  `show/on/display/include/true/1`, `hide/off/false/0`, and
  `suppress/exclude`
- normalize saved/serialized actions to canonical `show`, `hide`, or `suppress`

Generated TSV should include only valid non-default rows:

```text
record_id<TAB>feature_type<TAB>qualifier<TAB>value<TAB>action
```

Do not write a header unless existing Python behavior requires it. The parser
already tolerates a header row, but the current generated Web tables are
headerless.

## Run Analysis Integration

Replace the current `featureVisibilityOverrides` serialization with
serialization of `featureVisibilityRules`.

Generate `/web_feature_visibility_table.tsv` only when at least one valid rule
exists, then pass:

```text
--feature_visibility_table /web_feature_visibility_table.tsv
```

Keep the generated file path unchanged so existing run-info/session file mapping
continues to work.

The same table path must still be passed to browser protein extraction for
LOSAT/BLASTp caches, preserving the current `suppress` behavior.
The LOSAT/BLASTp cache key should be based on the serialized TSV content, not on
the derived compatibility cache.

## Session, History, And Reset

Session data:

- Save `featureVisibilityRules`.
- Load `featureVisibilityRules` when present.
- Migrate old `featureVisibilityOverrides` to editor hash rules when rules are
  absent.
- Do not save `featureVisibilityOverrides` in new sessions.

Undo/redo:

- Visibility edits from the left panel and Editor tab should be undoable.
- Prefer wrapping the existing action boundary instead of adding a parallel
  history mechanism.

Reset:

- Reset Settings clears `featureVisibilityRules`.
- Reset Settings also rebuilds/clears the derived `featureVisibilityOverrides`
  cache.
- Feature panel reset, if added later, should only clear these rules when the
  user explicitly requests visibility reset.

## Implementation Sequence

1. Add parser/serializer/normalizer/cache helpers and focused unit tests.
2. Add `featureVisibilityRules` state.
3. Convert runtime `featureVisibilityOverrides` use to a derived compatibility
   cache rebuilt from `featureVisibilityRules`.
4. Wire session/history/reset to save and restore `featureVisibilityRules` only,
   while migrating old `featureVisibilityOverrides` on load.
5. Change Editor visibility actions to upsert/remove editor hash rules and keep
   current instant SVG preview for those exact-feature changes.
6. Update run-analysis TSV generation and LOSAT/BLASTp cache keys to use the
   serialized rule array.
7. Add the left-panel Feature Visibility UI inside the Features card.
8. Add import/export controls only if they reuse the helper code cleanly.
9. Run focused Web tests and one packaging/config test.

## Tests

Preferred tests:

- parser accepts headerless TSV rows.
- parser skips blank/comment lines and header rows.
- parser accepts Python-compatible action aliases and serializes canonical
  actions.
- parser rejects missing required fields and extra columns consistently with
  Python behavior.
- serializer preserves rule order and emits `show`, `hide`, `suppress`.
- Editor `On/Off/Suppress/Default` upserts or removes the expected hash rule.
- new Editor hash rules insert at the top, and existing Editor hash rules keep
  their position when updated.
- session load migrates old `featureVisibilityOverrides`.
- new session save includes `featureVisibilityRules` and omits
  `featureVisibilityOverrides`.
- derived compatibility cache rebuilds from editor/hash rules.
- run-analysis emits `--feature_visibility_table` and the generated TSV content.

Use existing Web test patterns. Do not add broad screenshot tests for this
change unless a regression cannot be covered otherwise.

## Non-Goals

- Do not change Python visibility semantics.
- Do not invent a new visibility table format.
- Do not make feature type selection (`-k`) depend on this UI.
- Do not make JavaScript duplicate the full Python regex/matching engine.
- Do not refactor unrelated feature editor color or label code.
- Do not add presets or automatic rule suggestions.

## Open Decisions

- Whether TSV import belongs in the first implementation pass. It is useful, but
  manual add/edit plus export covers the primary workflow with less code.
