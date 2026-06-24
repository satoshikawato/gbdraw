# Session JSON Unification Plan

## Goal

Make `*.gbdraw-session.json` the single save/load format for the web app.

The session file should behave as a project file: loading it should restore every
user-controlled setting and enough state to continue work without remembering
which controls were changed before saving.

## Current State

The web UI currently exposes two JSON workflows:

- `Save Config` / `Load Config`
- `Save Session` / `Load Session`

Session export already embeds `config: buildConfigData()`, so the session format
is a superset of the config format. This makes the separate config workflow
mostly redundant and confusing because it is unclear which JSON file fully
restores the work.

Global feature stroke settings such as `adv.block_stroke_width` and
`adv.block_stroke_color` are already present in the session's `config.adv`
payload and are converted by CLI session loading. However, post-generation
stroke edits made from the feature or legend editors are tracked separately in
web state and are not currently stored as structured session state.

Both the web app and CLI currently use session version `27`. The CLI validates
against its current version and only supports explicitly listed versions, so the
version bump must update both sides together and keep old session files readable.

## Definition Of "All Settings"

Save every user-controlled intent, not internal runtime objects.

Persist these categories:

- Generation settings: `form`, `adv`, `losat`, palette, color table edits,
  specific rules, qualifier priority rules, filters, and circular conservation
  settings.
- Inputs: uploaded GenBank, GFF3, FASTA, depth, BLAST, color, filter, and
  qualifier-priority files, including per-record linear input settings.
- UI preferences that affect resuming work: active mode, input type, zoom,
  selected result, canvas pan/padding, layout preferences, title and legend
  position caches, DPI, palette preview state, and active editor tab.
- Generated results: current SVG results, including the currently displayed SVG
  after direct DOM edits.
- Post-generation feature edits: feature fill colors, feature visibility,
  feature stroke overrides, label text overrides, label visibility, and label
  override context.
- Post-generation legend edits: legend order, added and deleted entries, renamed
  captions, legend colors, and legend stroke overrides.
- Stroke edits: global `adv.*stroke*` settings and structured post-generation
  stroke override state for legend-level and feature-level edits. These
  structured overrides are required so edits can survive regeneration, not only
  preview restoration from serialized SVG.
- Orthogroup metadata edits and selection state.
- Optional cached computation state, such as LOSAT raw cache entries, with the
  existing large-file confirmation behavior.

Do not persist these runtime-only values:

- DOM nodes, Vue refs, worker handles, Pyodide state, timers, promises, hover
  state, drag-in-progress state, loading flags, and transient errors.
- Derived values that are fully recomputed from persisted inputs and settings.

## Target UX

- Header JSON actions should be simplified to `Save Session`, `Load Session`,
  `Reset Settings`, and the session title control.
- `Save Config` and `Load Config` should be removed from the visible UI.
- `Load Session` should remain backward-compatible with old `gbdraw-session`
  files and legacy config JSON. If the imported JSON is not a
  `gbdraw-session` envelope but validates as a legacy config payload, load it as
  a settings-only session.
- New downloads should use `*.gbdraw-session.json` only.

## Schema Plan

Increase the web `SESSION_VERSION` and Python `CURRENT_SESSION_VERSION` from
`27` to `28`.

Keep version `27` readable:

```python
CURRENT_SESSION_VERSION = 28
SUPPORTED_SESSION_VERSIONS = frozenset({27, 28})
```

The web importer should likewise accept versions `27` and `28`, reject future
versions, and normalize older sessions through a small migration path before
applying them. The initial `27 -> 28` migration should add a default
`editorState` object when it is absent and leave existing `config`, `ui`,
`files`, `results`, `features`, `orthogroupState`, and `losatCache` data in
place.

Keep the existing top-level shape:

```json
{
  "format": "gbdraw-session",
  "version": 28,
  "createdAt": "...",
  "title": "...",
  "config": {},
  "ui": {},
  "files": {},
  "results": [],
  "features": {},
  "editorState": {},
  "orthogroupState": {},
  "losatCache": {}
}
```

Add structured editor state. A compact top-level `editorState` object is
preferable to adding more unrelated keys to `features`.

```json
{
  "editorState": {
    "legend": {
      "entries": [],
      "deletedEntries": [],
      "originalOrder": [],
      "originalColors": {},
      "colorOverrides": {},
      "strokeOverrides": {},
      "addedCaptions": []
    },
    "featureStrokes": {
      "overrides": {
        "stableFeatureKey": {
          "strokeColor": "#222222",
          "strokeWidth": 1.2,
          "originalStrokeColor": "#222222",
          "originalStrokeWidth": 1
        }
      }
    },
    "originalSvgStroke": {
      "color": null,
      "width": null
    }
  }
}
```

`featureStrokes.overrides` is part of the first implementation. Use the same
stable feature identity strategy as existing feature color and visibility
overrides wherever possible. If the current UI only applies a feature stroke as
a direct SVG edit, add the small bridge needed to record that edit in structured
state at the time it is made.

## Restore Precedence

Use a fixed precedence rule to avoid drift between serialized SVG and
structured editor state:

- `results[].content` is the source of truth for the immediate preview restored
  after loading a session.
- `editorState` is the source of truth for post-generation edit intent that
  must be re-applied after regeneration, legend reflow, layout refresh, or any
  operation that rebuilds SVG content.
- If serialized SVG and `editorState` disagree, prefer `editorState` the next
  time edits are re-applied. Do not mutate the restored preview immediately
  unless the normal post-load layout refresh requires it.

## Implementation Steps

### 1. Centralize session serialization

In `gbdraw/web/js/services/config.js`:

- Keep `buildConfigData()` as the generation-settings serializer.
- Add helpers such as `buildEditorStateData()` and `applyEditorStateData()`.
- Use these helpers from `exportSession()` and `importSession()`.
- Normalize imported editor state defensively, just as config import already
  normalizes `adv`, track slots, depth tracks, and LOSAT settings.
- Add `normalizeSessionData()` or equivalent migration logic so version `27`
  sessions are upgraded to the current in-memory shape before application.
- Reject future session versions with a clear error. Accept only explicitly
  supported versions.

In `gbdraw/session_io.py`:

- Set `CURRENT_SESSION_VERSION = 28`.
- Set `SUPPORTED_SESSION_VERSIONS = frozenset({27, 28})`.
- Keep CLI loading focused on `config` and `files`; ignore web-only
  `editorState` safely.

### 2. Persist missing editor state

Include these state fields in session export:

- `legendStrokeOverrides`
- `legendColorOverrides`
- `deletedLegendEntries`
- `originalLegendOrder`
- `originalLegendColors`
- `addedLegendCaptions`
- `originalSvgStroke`
- `editorState.featureStrokes.overrides`

On import:

- Clear the corresponding reactive objects before applying imported values.
- Validate colors and numeric stroke widths.
- Rehydrate `addedLegendCaptions` as a `Set`.
- Preserve compatibility with old sessions where `editorState` is absent.
- Do not overwrite restored serialized SVG content after import unless a layout
  refresh is required.
- Apply feature stroke overrides to structured state using the same stable keys
  used by feature color and visibility overrides.
- When a feature stroke edit is made through the UI, record it in structured
  state immediately instead of relying on the SVG snapshot alone.

### 3. Remove visible config JSON workflow

In `gbdraw/web/index.html`:

- Remove the `Save Config` and `Load Config` buttons.
- Remove the config file input from the visible toolbar.
- Keep `Reset Settings`, `Save Session`, `Load Session`, and session title.

In `gbdraw/web/js/app/app-setup.js`:

- Remove `exportConfig` and `importConfig` from the template return surface if
  they become unused.

In `gbdraw/web/js/services/config.js`:

- Stop exporting `exportConfig()` for UI use.
- Keep legacy config import logic as an internal compatibility path used by
  session import.

### 4. Add legacy config compatibility through Load Session

Update `importSession()`:

- If `data.format === 'gbdraw-session'`, validate the version, migrate if
  needed, and use the normal session path.
- Otherwise, validate it as a legacy config payload with a dedicated
  `isLegacyConfigPayload(data)` check.
- Apply it with `applyConfigData(data)`.
- Restore palette state with `restorePaletteStateAfterConfigImport()`.
- Show a message such as "Legacy configuration loaded. Save as a session to use
  the current format."

This allows existing users to open old `gbdraw_config.json` files without
keeping two save buttons.

`isLegacyConfigPayload(data)` should be conservative:

- Accept only plain objects without a `format` field.
- Require at least one known config key such as `form`, `adv`, `losat`,
  `colors`, `palette`, `rules`, `qualifierPriorityRules`, `filterMode`,
  `whitelist`, `blacklistText`, `blastSource`, `losatProgram`, or
  `circularConservation`.
- Reject objects that primarily look like malformed session envelopes, such as
  objects with `files`, `results`, `features`, `ui`, or `editorState` but no
  valid `format`.
- Run the same circular and linear track-slot validation used by current config
  import before applying the payload.

### 5. Guard against future save omissions

Add tests that fail when a new default setting is not covered by session
serialization.

Recommended tests:

- A static or lightweight JS-source test that asserts every key from
  `createDefaultForm()` and `createDefaultAdv()` is serialized under
  `config.form` and `config.adv`.
- A source-level test that confirms the toolbar no longer contains `Save Config`
  or `Load Config`.
- A source-level test that confirms `legendStrokeOverrides` and
  `originalSvgStroke` are present in session export and import.
- A source-level or unit test that confirms feature stroke overrides are present
  in session export and import.
- A Python CLI session test confirming `CURRENT_SESSION_VERSION` is `28`,
  `SUPPORTED_SESSION_VERSIONS` includes both `27` and `28`, and
  `block_stroke_width` and `block_stroke_color` from `config.adv` still become
  CLI arguments.
- A Web importer test, or a source-level fallback test if browser coverage is
  not available, confirming version `27` sessions without `editorState` load
  through the migration path.
- A legacy config test confirming `Load Session` accepts a valid legacy config
  payload and rejects malformed session-like JSON without `format`.

If Playwright coverage is practical, add a browser test:

1. Open the web app.
2. Change `Block Stroke Color` and `Block Stroke Width`.
3. Save a session.
4. Load that session.
5. Assert the controls still show the changed values.
6. Generate or inspect the SVG and assert the stroke values are present.

### 6. Keep CLI compatibility intact

`gbdraw/session_io.py` already reads GUI session `config.adv` and maps
`block_stroke_width`, `block_stroke_color`, `line_stroke_*`,
`axis_stroke_*`, and other settings to CLI options.

During implementation:

- Update Python `CURRENT_SESSION_VERSION` to match the web `SESSION_VERSION`
  value of `28`.
- Keep `SUPPORTED_SESSION_VERSIONS` broad enough to accept version `27` and
  `28` sessions.
- Do not rename existing `config.form` or `config.adv` keys.
- Do not move generation settings out of `config`.
- Only add web editor state around the existing schema.
- Update the existing web/Python version synchronization test so it checks the
  current version and explicit support for the previous version.

## Acceptance Criteria

- The web UI has only one user-facing JSON save/load workflow: Session.
- Saving and loading a session restores all generation settings, including
  `Block Stroke Color` and `Block Stroke Width`.
- Saving and loading a session restores uploaded inputs, including circular
  GenBank/GFF3/FASTA/depth files, linear per-record file settings, BLAST files,
  color/filter/qualifier-priority files, and circular conservation files when
  present.
- Saving and loading a session restores UI resume state, including active mode,
  input type, zoom, selected result, canvas pan/padding, title and legend
  position caches, DPI, palette preview state, and active editor tab.
- Saving and loading a session restores current SVG results for immediate
  preview, including direct SVG edits captured in `results`.
- Saving and loading a session restores post-generation feature edits, including
  fill colors, visibility, feature stroke overrides, label text overrides,
  label visibility, and label override context.
- Saving and loading a session restores post-generation legend stroke edits.
- Saving and loading a session restores post-generation legend order, added and
  deleted entries, renamed captions, legend color overrides, and legend stroke
  overrides.
- Saving and loading a session restores orthogroup metadata edits and selection
  state.
- Saving and loading a session preserves optional LOSAT raw cache entries when
  the user confirms large embedded data.
- Legacy config JSON files can still be loaded through `Load Session`.
- Old version `27` session files without `editorState` still load successfully
  in both the web app and CLI.
- New session files remain usable by CLI session loading.
- Tests cover serialization of core settings, editor state, legacy config import,
  and version `27`/`28` compatibility.

## Risks And Mitigations

- Session files may grow larger as more state is stored.
  Keep the existing large embedded-data confirmation and only store structured
  editor state, not duplicated DOM snapshots beyond `results`.
- Direct SVG edits and structured editor state can drift.
  Treat structured state as the source for regeneration and serialized SVG as
  the source for preview restoration.
- Removing config buttons may surprise existing users.
  Keep legacy config import through `Load Session` and use a clear success
  message.

## Non-Goals

- Do not change the CLI TOML configuration model.
- Do not remove the `config` object from session JSON.
- Do not persist runtime-only browser or worker state.
- Do not make sessions deterministic byte-for-byte because embedded files,
  timestamps, and current SVG serialization make that impractical.
