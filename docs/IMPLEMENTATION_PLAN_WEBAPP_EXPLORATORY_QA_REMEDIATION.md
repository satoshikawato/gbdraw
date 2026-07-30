# Implementation plan: Web app exploratory QA remediation

Date: 2026-07-30

Status: completed (2026-07-31)

Intended reader: Codex in the next implementation session. This is an internal
execution plan, not user-facing documentation.

Do not run `avoid-ai-writing` over this file. Preserve operational detail and
machine-oriented redundancy that helps a later Codex session execute the work.

## Objective

Fix the confirmed P1, P2, and P3 findings from the July 30 exploratory Web app
audit. The work covers Circular and Linear Custom Track Slots, render result
topology, session round trips, editor state, resource I/O, standalone
interactive SVG metadata, and controls that are currently ineffective or
hidden.

The next session should treat this file as its task checklist. Do not replace it
with a shorter plan before implementation. Update `Status` only after all exit
criteria pass.

## Repository constraints

- Read `AGENTS.md`, `CLAUDE.md`, and `gbdraw/web/CLAUDE.md` before editing.
- Keep the Web app as a no-build single-page app.
- Use the typed schema-5 request path. Do not add an argv-shaped Web render
  path.
- Do not add external runtime dependencies.
- Preserve unrelated worktree changes. Inspect production and test diffs
  separately.
- Do not edit `CLAUDE.md` or a `SKILL.md` as part of the bug fixes described
  here.
- Do not update `tests/reference_outputs/` unless a reviewed default geometry
  change requires it. The Circular geometry fixes below must leave defaults
  unchanged.
- Use synthetic inputs or these safe fixtures for browser checks:
  - `tests/test_inputs/NC_001416.gb`
  - `tests/test_inputs/HmmtDNA.gbk`
  - `tests/test_inputs/BGC0000708.gbk`
  - `tests/test_inputs/BGC0000709.gbk`
- Do not introduce a dependency on WSSV or another pathogen fixture. Generate
  small GenBank, depth TSV, annotation, and comparison fixtures in the test
  temporary directory when an existing safe input is insufficient.
- Prepare the gitignored browser wheel after Python changes and before the
  affected Playwright checks:

  ```bash
  python tools/prepare_browser_wheel.py
  ```

- Finish the implementation as one session/commit, as required by
  `AGENTS.md`. Work packages below are review boundaries, not separate commits.

## Audit baseline

The audit ran against the current branch and browser wheel
`gbdraw-0.14.0b0-py3-none-any.whl`.

- Focused Python track/request tests: 125 passed.
- Node Web tests: 43 passed, 1 failed.
- Chromium Playwright: 41 passed, 8 failed.
- Confirmed automated product failures:
  - one Linear diagram appears as two Results;
  - the Vibrio Gallery interactive SVG exceeds the 40 MiB test limit;
  - a standalone feature FASTA action can expose an internal hash ID.
- Largest Gallery artifact baseline:
  - `vibrio-harveyi-group-collinear.svg`: 89,250,706 bytes;
  - decoded interactive metadata: approximately 277 MB;
  - `vibrio-harveyi-group-collinear.gbdraw-session.json.gz`: 96,434,798
    bytes;
  - expanded Vibrio session JSON: 450,464,802 bytes.
- Some other failures are stale expectations or missing fixtures. One
  version-33 consecutive-import failure was not fully classified.

Rerun the baseline at the start of the next session. Do not assume the counts
are unchanged after other worktree changes.

## Finding ledger

| ID | Priority | Finding | Work package |
| --- | --- | --- | --- |
| P1-01 | P1 | One logical diagram produces a plain SVG Result and an interactive SVG Result. | A |
| P1-02 | P1 | Editing a Linear or Circular slot ID replaces the keyed Vue row after the first character. | B |
| P1-03 | P1 | Closing and reopening Circular Custom Track Slots resets the stack. | B |
| P1-04 | P1 | UI mutations can create slot stacks that fail only in Python. | C |
| P1-05 | P1 | Inactive stacks and disabled rows are written to raw config but lost on session load. | D |
| P1-06 | P1 | Six visible Circular geometry fields never reach the render request. | F |
| P1-07 | P1 | Generate clears feature stroke overrides and a failed Generate clears the last good preview. | E |
| P1-08 | P1 | A retained depth file can re-enable Depth after `Show Depth` is off. | C |
| P1-09 | P1 | A manually added Circular Pairwise comparison slot can succeed while drawing no ring. | C |
| P2-01 | P2 | Generate serializes inactive-mode files and rereads active files several times. | H |
| P2-02 | P2 | The largest Gallery interactive SVG is 89.25 MB, its decoded metadata is about 277 MB, and its session is 96.43 MB gzip/450.46 MB expanded. | I |
| P2-03 | P2 | Per-mode comparison/axis profiles live only in a private Map and are not saved. | D |
| P2-04 | P2 | The 50 MB session warning counts bytes that are not saved and omits bytes that are saved. | H |
| P2-05 | P2 | Standalone single-feature FASTA can retain an internal `h_...` header. | G |
| P2-06 | P2 | Structured errors can display `[object Object]` and raw Python tracebacks. | G |
| P2-07 | P2 | `none`, null, and named colors appear as black in native color inputs. | G |
| P3-01 | P3 | The supported Spacer renderer is filtered out of both Custom Track UIs. | J |
| P3-02 | P3 | Several supported Annotation track options have no Web controls. | J |
| QA-01 | Follow-up | Existing test-contract drift and the version-33 import case need classification. | K |

### P1-04 cases that must all be covered

- Linear permits more than one enabled Features slot.
- Circular permits a Features duplicate in a state where a visible feature
  underlay requires exactly one enabled Features slot.
- Depth can be added without a corresponding logical depth source.
- Annotation can be added with no annotation set and an empty `set_id`.
- An overlay Annotation can point at a missing, disabled, annotation, ticks, or
  Spacer anchor.
- Moving an Annotation to overlay hard-codes `anchor_slot='features'`.
- Changing a renderer can invalidate a previously valid row without a local
  warning.
- Filtering disabled Circular rows does not rebase the Axis index.
- A source-less `sequence_conservation` row is silently skipped.

## Contracts to establish before editing

These decisions are part of the plan. Change them only if source inspection
finds a hard incompatibility, and record the reason in this file.

### Result and export contract

- One logical diagram equals one preview Result.
- Circular single and grid produce one Result.
- Circular batch produces one Result per rendered record.
- Linear produces one Result for the ordered layout.
- The Web render request asks Python for plain SVG.
- `Download SVG` serializes the committed preview.
- `Download Interactive SVG` adds standalone interactivity in the browser from
  the committed preview and committed feature metadata.
- CLI and Python API `interactive_svg` output behavior does not change.
- A supported old Web session containing `.svg` plus `.interactive.svg` for the
  same prefix loads as one logical preview. Prefer the plain SVG as the editable
  preview.

### Track editor state contract

- Panel disclosure, custom-stack activation, and explicit Reset are separate
  actions.
- A slot has a transient immutable editor key that is not its public `id`.
- The public slot `id` remains editable and is the canonical anchor/reference
  identity.
- Disabled rows remain in Web draft state and session state but are absent from
  the render request.
- Axis position is stored against the draft array and rebased against the
  enabled array when a request is emitted.
- The shared JavaScript validator owns the UI/draft constraints enumerated in
  P1-04, P1-08, and P1-09. The request builder must call it; disabled Add
  buttons and row warnings are not the only enforcement point.
- The Python typed decoder and renderer remain the final authority for
  canonical schema and render invariants. Do not duplicate every Python
  validation rule in JavaScript.

### Session authority contract

- `renderRequest` and its referenced resources own the last committed render.
- Web config owns editable UI drafts, including inactive Custom Track stacks
  and disabled rows.
- `editorState` owns committed editor overrides.
- Mode-profile snapshots are explicit session state, not an unexported Map.
- Web file bindings preserve user-owned files in both modes. Generate transfers
  only active render resources; Save Session preserves inactive-mode inputs
  once without duplicating their bytes.
- Current-writer sessions store each file payload once in `resources`.
  `webFiles` binds UI file slots to resource IDs. Do not write a second base64
  copy under legacy `files`.
- Older sessions that store raw `files` remain readable.

Current Web `SESSION_VERSION` and Python `CURRENT_SESSION_VERSION` are both 39.
This work changes the writer contract. Version selection is based on the last
release-backed writer, not merely the largest number in the worktree:

1. Recheck both constants, the release history, and
   `docs/IMPLEMENTATION_PLAN_316_LOSAT_OPT_OUT.md`.
2. If 39 is still the last release-backed writer, make the combined next
   contract version 40. If issue 316 already introduced an unreleased version
   40 on the branch, extend that same version-40 contract and its one migration
   from 39; do not create version 41.
3. If version 40 has been released or otherwise emitted as a supported public
   writer, use version 41 and retain version 40 as a supported input.
4. Change `gbdraw/session_io.py` and
   `gbdraw/web/js/services/config.js` together. Their current-version and
   supported-version sets must remain identical under
   `tests/test_session_io.py`.
5. Add only release-backed predecessor versions to the supported set. Do not
   support a branch-only intermediate contract.

### Render transaction contract

- `results`, selected Result, feature catalogs, editor selections, and editor
  overrides describe the last committed successful render.
- Generate builds candidate state without clearing committed state.
- Validation, worker, cancellation, stale-result, and post-processing failures
  discard candidate state and leave committed state unchanged.
- A successful render commits the candidate atomically.
- Feature fill and stroke overrides reapply by stable feature identity. A
  positional `file0_f2` key alone is not sufficient across regeneration.
- Unmatched overrides remain unapplied and may be pruned only after a successful
  commit confirms that their source feature no longer exists.

### P3 publication decision

- Expose Spacer in both Custom Track UIs.
- Expose Annotation `marks`, `lane_gap_px`, `padding_px`, and `cover_anchor`.
- Keep `style_override` API/session-only in this change. The full
  `RegionAnnotationStyle` contains nested hatch and label style objects and
  does not have a compact, understandable row editor. Preserve imported
  `style_override` values without data loss. Do not silently discard them when
  another Annotation field is edited.
- Record the API-only `style_override` decision in the final implementation
  summary. A later dedicated style-editor task can publish it.

## Dependency order

Use this implementation order:

1. Work package 0: capture the baseline and add focused failing tests.
2. Work packages B and C: stabilize row identity, panel state, and the common
   draft validator.
3. Work package D: establish session authority, synchronized Python/Web writer
   versions, resource bindings, history, and mode-profile persistence.
4. Work package I, phase 1: define the single normalized feature-catalog schema
   and builder used by Web preview state and standalone interactive SVG.
5. Work package A: switch Web preview output to one plain SVG per logical
   diagram and return the phase-1 catalog beside it.
6. Work package H: consume that bridge catalog, remove post-render parsing,
   deduplicate file payloads, and calculate exact session size. H depends on D,
   I phase 1, and A.
7. Work package E: make Generate transactional using the complete Result and
   catalog response from H. Do not implement E against the old reparse path.
8. Work packages F, G, and J: finish geometry, error/color/FASTA, Spacer, and
   Annotation work. F and J depend on C and D.
9. Work package I, phase 2: use the same catalog for metadata v3, regenerate
   Gallery artifacts, and enforce SVG/session size limits.
10. Work package K: classify stale tests and run the complete verification
    matrix.

Do not start P3 UI expansion before the common validator and session draft
round trip work. Do not switch the Web request to plain SVG until the phase-1
catalog builder is available; otherwise `RequestRenderResult.interactive_context`
is null and the editor/export metadata path regresses.

## Work package 0: baseline and regression scaffolding

### Inspect

- `git status --short`
- `git diff --stat`
- existing changes in every file before editing it
- current browser wheel version and source version
- available Playwright installations

### Run

```bash
node --test tests/web/*.test.mjs
pytest tests/test_circular_track_slots.py tests/test_linear_track_slots.py \
  tests/test_web_request_render.py -q
npx playwright test --project=chromium --workers=3
```

If Chromium fails with the known sandbox error, rerun the same command with the
required sandbox escalation.

### Add failing tests before production changes

Prefer focused owner tests. Add new files only where an existing test file has
no clear owner.

- `tests/web/session-request.test.mjs`
  - plain SVG output for new Web requests;
  - Circular Axis rebase after disabled rows;
  - depth resources omitted when neither simple nor custom Depth is active;
  - all six Circular geometry values change the request.
- Add `tests/web/track-slot-validation.test.mjs`
  - every P1-04, P1-08, and P1-09 UI-representable invalid state;
  - row/global issue ownership;
  - enabled-slot projection and Axis rebase;
  - unmanaged comparison/conservation rejection;
  - confirmation that Python-only semantic validation is not duplicated.
- `tests/web/mode-profiles.test.mjs`
  - export/import of both mode snapshots and managed flags.
- `tests/web/session-authority.test.mjs`
  - active canonical state plus inactive/disabled Web draft state.
- `tests/web/feature-color-actions.test.mjs`
  - stable feature override keys and reapply behavior.
- Add `tests/web/custom-track-slots.playwright.spec.js`
  - multi-character slot ID typing in both modes;
  - Circular Close/Open preservation;
  - explicit Reset;
  - early row validation;
  - inactive and disabled session round trips.
- Add `tests/web/render-state.playwright.spec.js`
  - one logical Result;
  - successful regeneration preserves overrides;
  - failed regeneration preserves the last preview and editor state;
  - interactive SVG download remains functional.
- Add or extend a resource-I/O unit test with a counting `File` implementation.
- Add a Python bridge test in `tests/test_web_request_render.py` for exactly one
  preview result per logical output.

Keep the initial regression tests focused on observed behavior. Do not encode a
large implementation object shape unless that shape is part of a persisted
contract.

## Work package A: one Result per logical diagram

Findings: P1-01.

### Production files

- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web_support/request_render.py`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/export.js`
- `gbdraw/api/render.py` only if a Web-only fix cannot preserve the public
  export contract

### Changes

1. Change new Web render output from `formats: ['interactive_svg']` to
   `formats: ['svg']`.
2. In the Web bridge, select the base SVG for each `RequestRenderResult`.
   Do not append every `.svg` path as a separate preview.
3. Leave `save_figure_to(..., formats=['interactive_svg'])` behavior unchanged
   for CLI/Python callers.
4. Add a pure session-result normalizer:
   - group `name.svg` and `name.interactive.svg`;
   - keep one preview entry;
   - keep unrelated batch prefixes separate;
   - do not collapse a legitimately distinct prefix just because it contains
     the word `interactive`.
5. Apply the normalizer when loading supported old sessions and when writing a
   current session.
6. Confirm `downloadInteractiveSVG()` builds from the selected committed plain
   SVG and feature metadata.
7. Do not enrich or write an interactive SVG for ordinary Web preview.
   Retain the in-memory phase-1 feature catalog required by H through the
   explicit Web bridge flag; metadata production is not a second Result.

### Acceptance

- Circular single: 1 Result.
- Circular one-record grid: 1 Result.
- Circular two-record grid: 1 Result.
- Circular two-record batch: 2 Results.
- Linear one or multiple records: 1 Result.
- Current sessions contain only those logical Results.
- An old session with paired plain/interactive Results displays one preview per
  logical prefix.
- Each Result has exactly one matching phase-1 catalog item in the bridge
  response.
- Plain and interactive downloads both work.
- Feature popup/search in the downloaded interactive SVG still works.

### Focused verification

```bash
node --test tests/web/session-request.test.mjs
pytest tests/test_web_request_render.py -v
npx playwright test tests/web/render-state.playwright.spec.js \
  --project=chromium
```

## Work package B: stable slot rows and non-destructive panel state

Findings: P1-02 and P1-03.

### Production files

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/linear-track-slots.js`

### Changes

1. Add a transient editor-key allocator owned by each track-slot editor.
   A `WeakMap<object, string>` plus a monotonic counter is preferred because it
   cannot leak into canonical requests or sessions.
2. Expose `trackSlotEditorKey(slot)` to the template and use it as the Vue key.
   The array index may be included for the Axis pseudo-row only, not for a slot
   row.
3. Add `circularTrackSlotsPanelOpen` beside
   `linearTrackSlotsPanelOpen`.
4. Make the Circular caret button toggle only panel disclosure.
5. Add an explicit `Use custom stack` control for Circular, matching Linear.
6. Change `setCircularTrackSlotsEnabled(true)`:
   - initialize from simple controls only when the draft stack is empty;
   - never reset a nonempty draft.
7. Keep `Reset` as the only action that replaces the current stack from simple
   controls.
8. Keep disclosure state out of the canonical request. It may be saved in
   `ui` if existing panel disclosure state is already persisted; otherwise it
   may remain transient.

### Acceptance

- Typing `custom_features` into either slot ID field produces the complete
  string without losing focus.
- Reordering, renderer changes, and duplicate/remove operations do not change a
  row's editor key.
- Circular Close/Open repeated ten times does not mutate the stack.
- Disable/enable custom stack does not mutate a nonempty stack.
- Reset replaces the stack once and is undoable through existing history if
  track-stack edits currently participate in history.

## Work package C: Custom Track validation, binding, Axis, and depth intent

Findings: P1-04, P1-08, and P1-09.

### Production files

- `gbdraw/web/js/app/track-slot-validation.js`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/linear-track-slots.js`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/index.html`

### Validator design

Add one pure validation entry point that returns both a normalized projection
and structured row/global issues. Suggested shape:

```js
validateCustomTrackPlan({
  mode,
  slots,
  axisIndex,
  trackType,
  depthTrackCount,
  annotationSetIds,
  visibleFeatureUnderlays,
  conservationSeries
})
```

Return:

```js
{
  draftSlots,
  enabledSlots,
  emittedAxisIndex,
  rowIssues: new Map(),
  globalIssues: []
}
```

Throw a typed boundary error only in the request builder. UI callers use the
same issue list to disable actions and display row messages.

### Structural checks

- unique nonempty public IDs;
- supported renderer and side;
- geometry syntax;
- Linear maximum of one enabled Features slot;
- Circular exact-one Features requirement only when current visible feature
  underlays require it;
- valid overlay anchor:
  - exists;
  - enabled;
  - drawable;
  - not Annotation;
  - not Ticks or Spacer;
- valid Annotation `set_id` from the current annotation-set inventory;
- valid Depth logical `track_index`;
- source-bound `sequence_conservation`;
- renderer-specific params after renderer change.

### UI behavior

1. Disable Add Annotation when there are no annotation sets.
2. Disable Add Depth when no logical depth source exists.
3. Disable Add/Duplicate Features when it would violate the current mode's
   Features constraint.
4. Replace free-text overlay anchor entry with a select of eligible enabled
   slots. Preserve an imported unknown value visibly as invalid until the user
   selects a valid anchor.
5. When an Annotation is moved to overlay:
   - choose the only eligible anchor automatically;
   - otherwise leave it unbound and show a row error;
   - never hard-code `features`.
6. Remove unmanaged `sequence_conservation` from the general Circular Add
   menu. Comparison Series remains its owner.
7. Keep imported unmanaged conservation rows visible but invalid; do not
   silently discard them.
8. Show concise row errors before Generate.
9. Block worker startup if any enabled row or Axis boundary is invalid.

### Axis rebase

Generalize the existing Linear enabled-slot Axis helper or add an equivalent
Circular helper. Define Axis as the boundary after `axisIndex` draft rows.
When disabled rows are filtered, the emitted index is the number of enabled
draft rows before that boundary.

Test:

- no disabled rows;
- disabled row before Axis;
- disabled row after Axis;
- all rows before Axis disabled;
- all rows disabled;
- Axis at 0;
- Axis at `slots.length`;
- an on-Axis Features/Annotation row if represented specially.

### Depth intent

Compute `depthRequested` before `buildDepthResources()`:

- simple mode: `form.show_depth`;
- custom mode: at least one enabled Depth slot;
- inactive depth files alone do not request a depth resource.

Session file preservation is separate from render-resource inclusion. Turning
off `Show Depth` must not delete the uploaded file.

### Acceptance

- All P1-04 cases fail in JavaScript before worker initialization.
- Errors identify the row and field without a Python traceback.
- Disabling a row can make the draft valid without deleting the row.
- Circular Axis index is in the enabled-slot range and preserves the intended
  inside/outside boundary.
- Turning `Show Depth` off while retaining a file produces no Depth track and
  no depth resource in the render request.
- Managed Circular comparison Series still produces its rings.

## Work package D: session authority, inactive state, and mode profiles

Findings: P1-05 and P2-03.

### Production files

- `gbdraw/session_io.py`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/mode-profiles.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/gallery-session-migration.js`
- `tests/test_session_io.py`
- `tests/web/session-authority.test.mjs`
- `tests/web/history.test.mjs`
- `tests/web/history-inputs.test.mjs`
- `docs/SESSION_COMPATIBILITY.md`

### Session writer/version

Advance the current Web session writer as described under the authority
contract. Change Python `CURRENT_SESSION_VERSION`, Python
`SUPPORTED_SESSION_VERSIONS`, Web `SESSION_VERSION`, and the Web supported set
as one operation. Add the immediately preceding release-backed writer to the
migration set. Keep versions 27-33 and request schemas 1, 2, and 5 only where
the current compatibility document requires them.

Add tests that:

- read both constants from source and require equality;
- require identical supported-version sets in Python and Web;
- promote the last release-backed writer exactly once;
- validate the new current-writer inventory, including `webFiles`, draft
  stacks, profiles, editor state, and the absence of legacy duplicated
  `files`;
- reject future and branch-only intermediate versions;
- load a version-39 fixture through Python and Web and preserve its active
  behavior.

### Draft-state merge

On load:

1. Validate and project the canonical request.
2. Restore the last committed render state.
3. Overlay Web-owned draft state from current-writer config:
   - `circular_track_slots_enabled`;
   - complete `circular_track_slots`, including disabled rows;
   - Circular draft Axis index;
   - `linear_track_slots_enabled`;
   - complete `linear_track_slots`, including disabled rows;
   - Linear draft Axis index;
   - simple Circular geometry values.
4. Validate that an active draft's enabled projection matches the canonical
   request written by the same current writer.
5. Fail closed on a mismatch in a current-writer session. For an older session
   without draft state, synthesize a draft from the canonical request.

Do not use canonical null tracks to overwrite a stored inactive stack.

### Mode profiles

Add explicit export/import APIs to `createModeProfileStateManager()`:

```js
manager.exportState()
manager.importState(payload, activeMode, adv)
```

Persist both mode snapshots and managed flags. Include:

- e-value;
- minimum bitscore;
- identity;
- alignment length;
- mode-specific Axis color fields already owned by the manager.

On old sessions without profiles, derive the active-mode snapshot from restored
`adv` and initialize the other mode from current defaults.

### History and Reset

Include both mode profiles, both complete Custom Track drafts, draft Axis
indexes, activation flags, and P3 Annotation/Spacer params in the history
snapshot. Undo/Redo must restore them without passing through the canonical
request projection.

- The Custom Track `Reset` button resets only that mode's draft stack and Axis
  to the values derived from its simple controls. It does not reset the other
  mode, mode profiles, or uploaded files.
- Settings Reset resets both mode profiles and both mode drafts to defaults,
  including new P3 params. Preserve the existing file-retention scope; add an
  assertion rather than changing it incidentally.
- Undo after a Custom Track Reset restores inactive rows, disabled rows, order,
  Axis, and nested `style_override`. Redo reapplies the reset.
- A Save/Load operation clears or establishes history exactly as the current
  documented behavior requires; it must not leave snapshots that refer to the
  pre-load session.

### Web file bindings

Move current-writer Web file persistence to resource references:

- Circular primary GenBank or GFF3+FASTA;
- Circular depth matrix;
- Circular comparison BLAST/FASTA sources;
- shared color/filter/qualifier files;
- Linear record files and depth rows;
- Linear uploaded comparison files;
- current canonical comparison artifacts where they are user-restorable state.

Use stable record/comparison UIDs in `webFiles`. Active canonical request
resources and Web file bindings must refer to the same resource entry rather
than duplicate bytes.

Define file-payload identity as byte length plus SHA-256 of the cached bytes:

- the same `File` object used in several roles has one resource payload;
- separate `File` objects with byte-identical content also share one payload;
- equal filename/size/date without equal content never deduplicates;
- compute the digest from `readFileBytes()` without a second `arrayBuffer()`
  call;
- each `webFiles` binding retains its own display name, MIME type, and
  `lastModified`, so sharing payload bytes does not erase UI identity;
- canonical render references and inactive Web bindings may point to the same
  resource ID;
- resource IDs remain opaque and collision-safe; do not expose the digest as a
  biological or user-facing identifier.

Test every file category listed above, including two separately constructed
byte-identical `File` objects and two same-metadata/different-byte objects.

### Acceptance

- Inactive Linear and Circular stacks survive Save, fresh page, Load, and
  reactivation.
- Disabled rows survive with ID, renderer, geometry, params, and order intact.
- The draft Axis index survives and emits the correct enabled Axis index.
- Circular identity 88 survives saving in Linear mode and switching back after
  Load.
- Inactive-mode input files survive a lossless Web session round trip.
- Each byte-distinct file payload appears once in current-writer session
  resources; binding-specific names and timestamps survive.
- Version-39 sessions retain their previous active behavior.
- Supported versions 27-33 still follow the documented migration limits.
- Unknown current-writer fields still fail the authority inventory.
- Python and Web current/supported version sets are identical.
- Undo/Redo and both Reset scopes follow the rules above for inactive,
  disabled, mode-profile, Spacer, and Annotation state.

## Work package E: transactional Generate and stable editor overrides

Finding: P1-07.

### Production files

- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/app/feature-editor/color-actions.js`
- `gbdraw/web/js/app/legend/stroke-actions.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/state.js`

### Changes

1. Stop clearing these at manual-run start:
   - `results`;
   - selected Result;
   - feature catalogs;
   - feature fill/stroke overrides;
   - editor selection;
   - committed legend/editor state.
2. Create a candidate render object scoped to the generation token.
3. Populate candidate Results, geometry metadata, feature metadata, record IDs,
   and run info.
4. Reapply editor overrides to candidate SVGs by stable biological/rendered
   feature identity.
5. Commit the candidate only after:
   - worker success;
   - response validation;
   - feature metadata completion;
   - SVG sanitization;
   - override reapplication.
6. On validation error, worker error, cancellation, stale result, or
   post-processing error, discard the candidate.
7. Keep error status separate from committed preview state.
8. Preserve current cancellation snapshots, but make normal failure follow the
   same state-preservation guarantee.

### Stable override identity

Use a compound identity:

- canonical persisted `recordKey` namespaces the record; do not use array index
  or accession because accessions may repeat and records may be reordered;
- `biologicalFeatureId` is built from the source feature before region clipping
  or reverse-complement display transforms, including source coordinate parts,
  strand, type, and a deterministic disambiguator for otherwise identical
  features;
- rendered SVG IDs and split segments map back to that biological key;
- fill/stroke overrides are stored against
  `recordKey + biologicalFeatureId`, not a rendered array position.

Reuse the existing stable hash only if it satisfies those rules. Otherwise
change the phase-1 metadata builder once and use the same identity in the
editor, session, and standalone metadata.

Store a positional key only as a compatibility alias during migration. For
supported version-39 editor overrides:

- one unique catalog match migrates to the compound key;
- no match remains preserved as an unresolved legacy override but is not
  applied to an unrelated feature;
- multiple matches or a key collision remain unresolved and produce a bounded
  migration warning; never guess.

### Acceptance

- Changing a feature stroke to `#ff00ff`, width 5, then Generate preserves the
  visible stroke and stored override.
- Fill overrides follow the same rule.
- A Generate failure caused by an invalid Annotation leaves:
  - SVG content;
  - Result list;
  - selected Result;
  - editor selection;
  - fill/stroke overrides;
  - legend edits
  unchanged.
- A successful render with a different input does not apply an old override to
  an unrelated feature with the same positional index.
- Overrides follow their biological feature when records are reordered, a
  region selector clips the display, or reverse-complement display is toggled.
- Duplicate accessions and two records containing the same feature index remain
  isolated by `recordKey`.
- Version-39 migration tests cover one match, no match, multiple matches, and a
  compound-key collision.

## Work package F: make Circular geometry controls effective

Finding: P1-06.

### Fields

- `feature_width_circular`
- `depth_width_circular`
- `gc_content_width_circular`
- `gc_content_radius_circular`
- `gc_skew_width_circular`
- `gc_skew_radius_circular`

### Production files

- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/index.html`
- Python files only if the existing structured Circular slot contract cannot
  represent the required scalar

### Changes

1. Validate every non-null value as finite and positive.
2. In simple Circular mode, only materialize explicit slots when at least one
   geometry shortcut is set.
3. Build the same implicit slot order the renderer would use for current
   Features/Ticks/Depth/GC controls.
4. Patch:
   - Features width in px;
   - every active Depth slot width in px;
   - GC content radius factor and width in px;
   - GC skew radius factor and width in px.
5. Emit the resulting structured slots and correct Axis index.
6. In active Custom Track mode, explicit row geometry wins. Disable the six
   simple fields with a direct explanation rather than applying two competing
   geometries.
7. `Reset` copies the simple geometry values into the new custom draft.
8. Keep all-null behavior canonical-request equivalent to the current default.

Prefer structured slot objects at the Web typed-request boundary if string
specification would lose nested or typed values. Keep string parsing for
compatibility input.

### Tests

- Each field independently changes the canonical request.
- Width/radius values survive current session round trip.
- Default/all-null request remains unchanged.
- Invalid zero, negative, `NaN`, and infinite values fail before worker launch.
- Python geometry metadata or final SVG geometry proves that each value affects
  the intended track.
- A combined simple-mode fixture with Features, Ticks, two Depth tracks, GC
  content, and GC skew retains every implicit slot in the same order and emits
  the correct Axis after geometry materialization.
- Combination tests cover underlay/overlay Annotation, a managed conservation
  Series, one-record Circular, multi-record grid, and batch. Geometry shortcuts
  must not drop, duplicate, or rebind those tracks.
- Enabling Custom Track mode disables the six shortcuts without mutating the
  explicit draft; Reset copies the shortcut values exactly once.
- Default reference SVGs remain unchanged.

## Work package G: errors, color adapters, and public FASTA identity

Findings: P2-05, P2-06, and P2-07.

### Error normalization

Production files:

- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/workers/diagram-generation-worker.js`
- `gbdraw/web/js/app/app-setup.js`

Add one bounded safe-text normalizer that accepts strings, Errors, Python error
objects, nested `{summary, details}`, and arbitrary objects. Never rely on
implicit object string conversion.

User-facing error rules:

- summary contains stage/type plus a safe message;
- details contain bounded safe diagnostic sections;
- raw traceback/stack is not displayed in the normal UI;
- cleanup notes remain attached without replacing the primary error;
- no uploaded sequence or full file content is included;
- `[object Object]` is impossible.

Keep raw traceback only in a development diagnostic channel if one already
exists and does not expose user data. Do not add a new public traceback panel.

### Color input adapters

Production files:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/color-utils.js`
- relevant feature-editor and advanced-control modules

Rules:

- null means Auto/Inherited, not black;
- `none` means No stroke/fill, not black;
- named colors resolve to a hex picker display without mutating state until the
  user edits;
- native `input[type=color]` receives only `#rrggbb`;
- an explicit Auto/None control restores null or `none`.

Cover global block/line/scale controls and the feature editor, including
`repeat_region`.

### FASTA identity

Production file:

- `gbdraw/web/js/services/standalone-interactivity-assets.js`

Always rebuild single-feature FASTA from:

- a public display protein/record ID;
- current feature description;
- normalized sequence.

Do not return an existing FASTA blob unchanged based only on a denylist scan.
Apply the same public identity to download filenames. Keep compatibility with
literal `\\n`, real newlines, and pre-v2 metadata.

### Acceptance

- Annotation and duplicate-Features errors have readable summaries.
- No normal error details contain `Traceback` or `[object Object]`.
- `strokeColor='none'` displays No stroke and remains `none` until edited.
- null global colors display Auto.
- Single-feature nucleotide and amino-acid Copy/Download never expose an
  internal `h_...` ID.
- Similarity-group bulk FASTA remains correct.

## Work package H: resource reads and exact session sizing

Findings: P2-01 and P2-04.

### Production files

- add `gbdraw/web/js/services/file-content-cache.js` if no existing owner fits
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/workers/diagram-generation-worker.js`
- `gbdraw/web_support/request_render.py`
- `gbdraw/api/request_render.py`
- `gbdraw/web_support/feature_metadata.py`

### File content cache

Use a `WeakMap<File, Promise<Uint8Array>>`:

- `readFileBytes(file)` reads once;
- `readFileText(file)` decodes cached bytes;
- base64 encoding consumes cached bytes;
- Pyodide staging consumes cached bytes;
- worker transfer clones only at the detaching transfer boundary;
- failed reads evict their cache entry;
- reset does not need to clear a WeakMap.

Do not cache by filename/size/date alone.

### Split render and session serialization

Replace overloaded `serializeFiles()` with separate operations:

- `serializeActiveRenderFiles(mode, state)`:
  - only files reachable from the active canonical request;
  - called by Generate;
  - no inactive-mode reads.
- `buildSessionResources(state, committedRequest)`:
  - all user-owned files for both modes;
  - called only by Save Session;
  - one resource payload per byte-length/SHA-256 identity defined in D;
  - emits `webFiles` resource bindings.

Remove the legacy current-writer raw `files` duplicate after the session version
advance. Keep its decoder for supported older sessions.

### Remove post-render rereads

Implement the Work package I phase-1 catalog builder first. Add a narrowly
scoped `include_feature_catalog=True` path to `render_request()` or its internal
prepared-render helper. The Web bridge uses it; existing CLI/Python callers
retain the false/default behavior.

When the flag is true:

1. Build the interactive/feature context from the same
   `PreparedDiagramRequest`, prepared inputs, records, and comparison metadata
   used for drawing.
2. Retain that context on each `RequestRenderResult` even though the requested
   output format is plain SVG.
3. Do not ask `save_figure_to()` to emit an interactive SVG.
4. Normalize the retained context through the one v3 catalog builder.

The current Web bridge response becomes:

```text
{
  results: [{ name, content }],
  metadata: {
    trackSlotGeometry,
    featureCatalog: {
      schema: 3,
      items: [{
        resultIndex,
        resultName,
        recordKeys,
        features,
        biologicalFeatures,
        orthogroups,
        annotations,
        comparisonMatches
      }]
    }
  }
}
```

Return exactly one catalog item for each logical Result and match by
`resultIndex` plus `resultName`. Empty arrays are allowed; a missing or
malformed catalog is not.

Then remove the ordinary successful path that calls
`extractGeneratedDiagramFeatures()` and reparses the input after rendering.
Do not retain a current-generation reparse fallback. Old sessions regenerate
through the current worker. A cache-skewed worker response without
`metadata.featureCatalog` fails the candidate render with a reload message and
leaves the committed preview intact.

For a simple non-LOSAT render, do not stage all uploaded inputs into the
main-thread Pyodide filesystem. Initialize/stage the main runtime only for a
helper or LOSAT operation that actually requires it.

### Exact warning

1. Assemble the actual current session payload.
2. Serialize/compress it once using the same bytes that will be downloaded.
3. Use the final gzip Blob size for the 50 MB download warning.
4. The confirmation message must say compressed session size.
5. Avoid a second `JSON.stringify`, gzip pass, or Blob copy after confirmation.
6. A cheap preflight warning for obvious raw input size may remain, but label it
   as input size and never substitute it for the final artifact size.

### Instrumented acceptance

- Circular Generate with an unused Linear 8 MiB file calls
  `arrayBuffer()` zero times on that Linear file.
- Linear Generate with unused Circular files calls zero reads on them.
- Each active primary File is read once on the main thread.
- The successful render path does not parse the same GenBank/GFF3 input again
  for feature extraction.
- Plain-SVG Web rendering returns one complete phase-1 feature catalog item per
  Result without writing an interactive SVG.
- A missing/malformed catalog fails transactionally and never invokes the old
  extraction helper.
- Save Session preserves inactive files while reading each uncached file once.
- Active resources are not duplicated in the saved session.
- The warning size equals the Blob that would be downloaded.
- Generate remains cancellable and transferred buffers are not reused after
  detachment.

Do not set wall-clock assertions in general CI. Track deterministic read counts,
payload bytes, and runtime initialization counts. Record representative local
timings in the final handoff.

## Work package I: compact interactive SVG metadata and Gallery artifacts

Finding: P2-02.

### Production files

- `gbdraw/web_support/feature_metadata.py`
- `gbdraw/web_support/request_render.py`
- `gbdraw/api/request_render.py`
- `gbdraw/session_io.py`
- `gbdraw/render/interactive_svg.py`
- `gbdraw/web/js/services/standalone-interactivity.js`
- `gbdraw/web/js/services/standalone-interactivity-assets.js`
- relevant feature/session catalog helpers
- `tools/refresh_gallery_sessions.py`
- `tests/test_refresh_gallery_sessions.py`

### Phase 1: shared catalog contract

This phase must land before A and H. Define one v3 catalog normalizer that
accepts the prepared Python interactive context and returns the exact
`metadata.featureCatalog.items[]` shape specified in H. It must not read an
input file or rebuild a `SeqRecord`.

Normalize v3 metadata:

- one biological feature payload per feature;
- rendered `features` contain DOM/render binding plus a reference to the
  biological payload, not a second complete payload;
- exclude non-interactive `source` features from the biological catalog;
- store nucleotide and amino-acid sequence once;
- do not repeat `translation` in qualifiers when the amino-acid sequence field
  already contains it;
- store whole-record sequence at most once if a supported action requires it;
- orthogroups and matches reference feature IDs instead of copying feature
  payloads;
- omit null/empty/default fields;
- keep public display IDs separate from internal matching IDs.

Use the same normalizer in Python-generated and browser-generated interactive
SVGs. The browser receives the normalized catalog from the render bridge and
embeds that catalog when downloading an interactive SVG; it does not construct
a second biological catalog from DOM state. Do not create two v3 wire formats.

### Phase 2: standalone encoding and functional requirements

Advance the interactive feature-popup metadata writer from v2 to v3. The
embedded runtime must read v2 and v3. Existing standalone v1/v2 SVGs remain
self-contained and do not need conversion.

Preserve:

- feature search;
- rich and simple popups;
- qualifier display;
- nucleotide/amino-acid Copy and Download;
- pairwise-match popups;
- similarity-group member tables and bulk FASTA;
- annotation interaction;
- old v2 decoding.

### Size tests

Add a metadata composition test with many repeated rendered/biological features
that proves payload references are used.

Keep the existing Gallery 40 MiB maximum. Add component measurements for the
largest example:

- total interactive SVG bytes;
- compressed metadata bytes;
- decoded metadata bytes;
- rendered feature count;
- biological feature count;
- gzip and expanded session bytes;
- serialized bytes by top-level session component: `results`, `resources`,
  `webFiles`, feature catalog/editor state, `losatCache`, and
  `losatDerivedCache`.

Required regenerated artifact targets:

- interactive SVG: below 40 MiB (41,943,040 bytes);
- decoded v3 metadata: at most 200,000,000 bytes;
- Vibrio session gzip: at most 90,000,000 bytes;
- Vibrio expanded session: at most 400,000,000 bytes.

After regeneration, set deterministic regression limits with headroom, but do
not loosen them above 220,000,000 decoded metadata bytes, 95,000,000 session
gzip bytes, or 420,000,000 expanded session bytes. If the measured artifact
cannot meet the required targets, P2-02 remains unresolved; do not simply raise
the old limits.

Add structural assertions:

- current-session `results` contains one plain SVG per logical diagram and no
  paired interactive Result;
- each byte-distinct resource payload appears once;
- a biological feature's qualifiers and sequences appear only in its catalog
  entry, not copied into rendered features, orthogroups, or matches;
- the required 59-pair current protein cache remains complete;
- component byte totals are reported when a limit fails.

### Gallery regeneration

Regenerate the affected interactive SVG and session through their owning
generator. Do not hand-edit generated SVG/session bytes.

Check:

- `gbdraw/web/gallery/examples/vibrio-harveyi-group-collinear.svg`
- its source SVG/session inputs;
- generated session, thumbnail, and `examples.json` references only if the
  generator changes them.

No screenshot recapture is required if the visible diagram and tutorial steps
do not change. Run the Gallery screenshot maintenance checks only if a visible
capture or tutorial state changes.

## Work package J: P3 Spacer and Annotation controls

Findings: P3-01 and P3-02.

### Spacer

Production files:

- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/linear-track-slots.js`
- `gbdraw/web/index.html`

Remove the UI-only Spacer filter. Keep runtime capability gating. Expose:

- Linear height and spacing;
- Circular width/radius/gaps through the existing generic geometry controls;
- ID, enabled state, order, side, and legend behavior consistent with other
  non-anchorable rows.

Spacer cannot be an Annotation anchor.

### Annotation track options

Expose:

- `marks`: checkboxes for line, bracket, band, highlight; empty selection means
  all marks and serializes as null/absent;
- `lane_gap_px`: nonnegative number, default 3;
- `padding_px`: nonnegative number, default 2;
- `cover_anchor`: overlay-only Boolean, default false.

Extend both slot serializers and normalizers. Preserve imported
`style_override` objects unchanged even though no style editor is added.
If current string slot specs cannot round-trip nested `style_override`, emit
structured slot objects in current schema-5 Web requests. Keep string specs as
compatibility input.

### Automated test ownership

- `tests/web/session-request.test.mjs`
  - serialize Spacer and all four published Annotation options in both modes;
  - omit default/null values canonically;
  - preserve an imported nested `style_override` unchanged after editing
    another field.
- `tests/web/session-authority.test.mjs`
  - Save, fresh-page Load, and re-save enabled, disabled, and inactive P3 rows;
  - deep-compare `style_override` before and after the round trip.
- `tests/web/custom-track-slots.playwright.spec.js`
  - add, edit, reorder, disable, reset, and reload Spacer/Annotation rows;
  - verify the controls remain connected to the selected row.
- `tests/test_circular_track_slots.py`,
  `tests/test_linear_track_slots.py`, and `tests/test_annotations.py`
  - prove Spacer changes resolved layout;
  - prove marks, lane gap, padding, and cover-anchor values reach the renderer
    and alter the intended geometry without changing defaults.
- `tests/web/history.test.mjs`
  - Undo/Redo P3 option edits and Custom Track Reset with nested imported style
    data intact.

### Acceptance

- Spacer can be added, rendered, reordered, disabled, saved, and loaded in both
  modes.
- Spacer changes visible/resolved geometry and does not become an anchor.
- Annotation mark filtering changes rendered marks.
- Lane gap and padding affect resolved annotation layout.
- `cover_anchor` affects overlay coverage only.
- All options survive a fresh-page session round trip.
- Imported `style_override` remains deep-equal after an unrelated Annotation
  edit, history operation, save, fresh-page load, and re-save.

## Work package K: test-contract cleanup, compatibility triage, and final pass

Finding: QA-01.

### Classify the previously failing tests

Do not update expectations merely to get green tests.

- Old CLI-argument expectations:
  - verify whether `lastRunInfo` intentionally reports session replay only;
  - update tests only after confirming the current product contract.
- Gallery session regeneration:
  - do not mutate a current schema-5 fixture to schema 2 without required
    legacy fields;
  - build a valid representative legacy fixture.
- Missing fixture:
  - replace the dependency with a synthetic safe fixture or remove the test if
    it no longer represents supported behavior.
- Version-33 consecutive import:
  - reproduce with a focused timeout and inspect pending promises/dialog state;
  - if a valid fifth import is ignored or hangs, fix the product;
  - if the fixture is not a valid supported session, correct the fixture;
  - retain the rollback checks from the four preceding invalid imports.
- Node Gallery migration failure:
  - align the fixture with the current writer's optional color-table contract;
  - do not invent a resource reference that the current session does not own.

Keep these changes logically separated in the diff from product code and new
regression tests, even though the session ends with one commit.

### Documentation generated by behavior changes

Update only documents whose contract changes:

- `docs/SESSION_COMPATIBILITY.md` for the new Web session writer and migration;
- `docs/TUTORIALS/8_Interactive_SVG_Sessions.md` for interactive metadata v3 if
  it names the current schema;
- user-facing Web help text where control behavior changes.

Do not add this audit's implementation lessons to CLAUDE or skill files.

## Full test matrix

### Node

```bash
node --test tests/web/*.test.mjs
```

Required focused coverage:

- result format and topology;
- slot structural/resource validation;
- Circular Axis rebase;
- session draft/canonical merge;
- mode-profile export/import;
- resource read counts;
- exact session Blob size;
- public FASTA identity;
- error formatting;
- color adapters;
- metadata v2/v3 compatibility.

### Python

```bash
pytest tests/test_web_request_render.py -v
pytest tests/test_circular_track_slots.py -v
pytest tests/test_linear_track_slots.py -v
pytest tests/test_annotations.py -v
pytest tests/test_api_request_render.py -v
pytest tests/test_session_io.py -v
pytest tests/test_gallery_session_semantics.py -v
pytest tests/test_refresh_gallery_sessions.py -v
```

Add focused tests rather than expanding unrelated reference-output tests.

### Browser

```bash
npx playwright test tests/web/custom-track-slots.playwright.spec.js \
  --project=chromium
npx playwright test tests/web/render-state.playwright.spec.js \
  --project=chromium
npx playwright test tests/web/depth-track-session.playwright.spec.js \
  --project=chromium
npx playwright test tests/web/linear-multi-record.playwright.spec.js \
  --project=chromium
npx playwright test tests/web/gallery-session-regeneration.playwright.spec.js \
  --project=chromium
npx playwright test tests/web/gallery-tutorial.playwright.spec.js \
  --project=chromium
npx playwright test --project=chromium --workers=3
```

Browser assertions must use synthetic/safe inputs and must not depend only on
console logs. Assert state, request payload, Result count, visible SVG
attributes, downloaded content, and session JSON as appropriate.

### Static/build

```bash
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/services/session-request.js
node --check gbdraw/web/js/services/config.js
node --check gbdraw/web/js/services/standalone-interactivity.js
node --check gbdraw/web/js/services/standalone-interactivity-assets.js
ruff check gbdraw/
python -m build
```

Run:

```bash
pytest tests/ -v -m "not slow"
```

after the focused suites pass.

## Manual browser checklist

Use a fresh browser context for session tests.

### Circular

- [ ] Generate `HmmtDNA.gbk`; confirm one Result.
- [ ] Type a multi-character slot ID.
- [ ] Close/Open the panel; confirm no state change.
- [ ] Disable/enable the custom stack; confirm no state change.
- [ ] Explicit Reset changes the stack.
- [ ] Disable rows before Axis; Generate succeeds with intended placement.
- [ ] Attempt Depth without a source; row error appears before worker startup.
- [ ] Attempt Annotation without a set; row error appears.
- [ ] Turn off Show Depth while keeping a file; no Depth ring appears.
- [ ] Change each of the six geometry fields and inspect the intended ring.
- [ ] Save inactive/disabled draft state; load in a fresh page; reactivate it.
- [ ] Add Spacer and Annotation advanced options; round-trip them.
- [ ] Confirm unmanaged Pairwise comparison cannot be added.

### Linear

- [ ] Generate `NC_001416.gb`; confirm one Result.
- [ ] Type a multi-character slot ID.
- [ ] Attempt to duplicate Features; action is disabled or locally rejected.
- [ ] Create a valid custom stack and round-trip enabled and disabled rows.
- [ ] Edit feature stroke; regenerate; confirm it remains.
- [ ] Trigger a validation failure; confirm the old preview/editor remain.
- [ ] Load `BGC0000708.gbk` and `BGC0000709.gbk`; confirm record rows,
  comparison binding, Save/Load, and regeneration.
- [ ] Add Spacer and Annotation advanced options; round-trip them.

### Export and performance

- [ ] Download plain SVG.
- [ ] Download interactive SVG and test popup/search/FASTA offline.
- [ ] Confirm public FASTA headers.
- [ ] Confirm Generate reads no inactive-mode file.
- [ ] Confirm Save Session preserves inactive inputs without byte duplication.
- [ ] Confirm the displayed session-size warning matches the downloaded gzip
  Blob.
- [ ] Inspect the largest Gallery artifact size without displaying sensitive
  biological results.

## Exit criteria

All items must be true:

- [x] Every P1, P2, and P3 ledger row has production code and an automated
  regression test or an explicit API-only product decision.
- [x] One logical diagram produces one preview Result.
- [x] Interactive SVG export still provides all supported actions.
- [x] Custom Track row editing, panel disclosure, activation, and Reset are
  independent and stable.
- [x] Every UI-representable P1-04, P1-08, and P1-09 case listed in this plan
  is rejected before worker startup; Python remains final authority for other
  semantic validation.
- [x] Inactive/disabled drafts and mode profiles survive a fresh-page session
  round trip.
- [x] Supported older sessions retain their documented behavior.
- [x] Python and Web current/supported Session version sets match, and the last
  release-backed writer has an explicit migration.
- [x] The six Circular geometry controls alter request and rendered geometry.
- [x] Successful Generate preserves editor overrides.
- [x] Failed Generate preserves the last successful preview and editor state.
- [x] Stable overrides remain isolated across duplicate accessions, reordered
  records, region display, and reverse-complement display.
- [x] Inactive-mode read count during Generate is zero.
- [x] The bridge returns the complete phase-1 catalog from prepared render data
  and the successful path does not reparse inputs.
- [x] Current session resources contain one copy per file payload.
- [x] Session warning uses actual downloadable gzip size.
- [x] No FASTA action exposes internal IDs.
- [x] No normal error panel shows `[object Object]` or a Python traceback.
- [x] Auto/None colors are not displayed as explicit black.
- [x] Largest checked-in interactive Gallery SVG is below 40 MiB.
- [x] Vibrio decoded metadata is at most 200,000,000 bytes; its session is at
  most 90,000,000 gzip bytes and 400,000,000 expanded bytes.
- [x] Spacer and selected Annotation options are usable in both modes and
  session-safe.
- [x] Undo/Redo, Custom Track Reset, and Settings Reset preserve their declared
  scopes for inactive/disabled/P3/profile state.
- [x] Node suite passes.
- [x] Focused Python suites pass.
- [x] Full Chromium Playwright suite passes, or any remaining external
  environment failure is documented with evidence.
- [x] `pytest tests/ -v -m "not slow"` passes.
- [x] `ruff check gbdraw/` passes.
- [x] Build succeeds.
- [x] No unintended `tests/reference_outputs/` diff exists.
- [x] Production and test diffs were reviewed separately.

## Final handoff

The next session's final response must include:

- implemented work-package list;
- unresolved item, if any, with exact blocker;
- before/after Result counts;
- before/after resource read counts;
- before/after largest Gallery SVG and session sizes;
- session writer version and migration coverage;
- focused and full test results;
- generated artifact changes;
- proposed commit title and English summary.

Suggested commit title:

```text
fix(web): repair track state, sessions, and render performance
```
