# Multi-drawing project: design (synthesis)

Revision 4.1, 2026-10-07 (§8.1 and §9 status refreshed). Orchestrator: session gbdraw-09. Code base read: `origin/dev` `b355b8ae`
(Session 46 shape: Phase E `ov80/plan.md` rev 3 §4, Owner answers QA-QC of 2026-10-07).
**Target release: v0.15.0** for drawings (Owner R5-1: 「Multiple drawing projectは次のバージョン(v0.15.0)に持ち越せばいいから、このセッションはまずデザインに徹してください。」).
**Two per-mode steps land in 0.14.0 first** and this design builds on them (§0.3):
1. **E1 + Q0** (this campaign, Owner R6-1): each mode keeps its own Result; dev's Session 45 gains one
   optional field `otherModeResult` so a Session saves both Results. 45 never ships: PR-1 rewrites it as 46.
2. **Per-mode settings store** (Phase E, Owner decision to gbdraw-9b; QC(A): 0.14.0 waits for it): every
   drawing-level setting and all Legend edits per mode, levels from `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` `level_P`;
   **Session 46** = top-level `modes.{circular,linear}` slices + the committed set + `otherModeResult`
   (dev-only 45 dropped, no 45 reader anywhere).
Drawings in v0.15.0 are therefore **Session 47**, read from 46 (each mode's slice and Result become one
drawing, §3.1) and from 27-44 (through Phase E's split, §3.2).
This document decides. The evidence (file:line) lives in the part files it cites. The level table and its notes are next to this file; the other part files, the audits, the findings list and the Owner replies stay in the campaign folder `gbdraw-baselines/multi-drawing-project-20261007/`, outside the repository:

| Part | File | Scope |
| --- | --- | --- |
| py | `design-parts/py.md` | Python document, API, CLI, validation, migration, performance |
| inputs | `design-parts/inputs.md` | resource table, per-drawing input drafts, Depth series, caches |
| contract | `design-parts/contract.md` | Product Contract re-check, receipt draft, format history, fixtures, docs |
| runtime | `design-parts/runtime.md` | Web store, reader conversion, History, switch, caches, PR sequence, §7 early per-drawing Results |
| ux | `design-parts/ux.md` (+ `ux/` screenshots and mockups) | layout, creation, delete, Reset, Export, words |
| levels | `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv`, `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.md` | level of every setting (270 rows), migration token per row |
| audits | `owner-coupling-phase-e-20261005/ov80/audit/{form-adv,editor,mechanisms}.md` | today's model |

Owner replies are verbatim in `owner/owner-replies.txt` (R1-x, R2-x, R3-x); questions in
`owner/questions-batch-*.txt`.

## 0. Decisions

### 0.3 Version path (agreed with Phase E, 2026-10-07)

| Step | Release | Owner | Session | Shape change |
| --- | --- | --- | --- | --- |
| E1 + Q0 | 0.14.0 | gbdraw-09 | 45 (dev-only, extended; never released) | optional top-level `otherModeResult` = the other mode's committed set (request, Results, catalog, run metadata, CLI invocation, generated Legend refs, applied palette, alignment receipt); one shared resource table; the explicit `transitionDiagramMode` |
| per-mode store PR-0a/0b/0c, PR-1 | 0.14.0 | gbdraw-9b | 46 | top-level `modes.{circular,linear}`, each slice = today's `config` (minus `modeProfiles` and LOSAT execution), `features`, the Legend part of `editorState` (+ `featureStrokes`), four `ui` keys (`layoutPreferences` slot, `canvasPadding`, `pendingPalette*`, `linearTypographyLinked`); top-level `config`/`features` retired; R2 rows lose `scope`; committed set top-level; `otherModeResult` unchanged; `ui.losatExecution`. Old Sessions: JS split on Web Load + Python twin (tokens copy/own/profile/layout/show-if-source/depth/result-mode/by-scope/by-side/by-leaf/by-binding) |
| drawings (this design) | v0.15.0 | this campaign | 47 | `drawings[]`; each 46 mode slice + its committed set → one drawing; the other mode's drawing only when it shows use (R1-2, rule in §3.1) |

Phase E sequence: layering A → **E1** → layering B → PR-0a → 0b → 0c → PR-1 → layering C → D.
Owner answers for 0.14.0 (Phase E QA-QC): one History in which a mode switch is a step, and Reset Settings
resets both modes (QA); no copy-from-other-mode action (QB); 0.14.0 is held for PR-1 (QC); PD-OI-063 is read
per mode. Per-drawing History and "Reset this drawing (+ others if checked)" (R3-4) arrive with drawings.


### 0.1 Owner decisions

| # | Decision | Source |
| --- | --- | --- |
| D1 | One Session (project) holds shared inputs and N drawings. Each drawing has its own mode, all settings, Legend edits, History and Result. Several drawings of one mode later. | 2026-10-07 |
| D2 | Python owns the project document (read, validate, migrate, write; CLI and API render one or all drawings). The Web owns drawing editing and UI. One typed contract. | 2026-10-07 |
| D3 | The per-mode store of `ov80/plan.md` rev 2 is not built; its "every reader names its context" step is reused for drawings. | 2026-10-07 |
| P | **Drawings share only input bytes and computed caches.** Everything the user edits belongs to one drawing: settings, palette, color rules, qualifier priority, label tables and filters, annotation sets, comparison plans, Legend edits, feature edits, decorations. Copies between drawings are explicit actions. | R1-1 |
| C-new | A new drawing starts blank: no files, its mode's defaults, no Legend edits. | R1-1, R3-2 |
| C-menu | "+ New drawing": (1) New Circular / New Linear (blank); (2) Duplicate (files, settings, Legend and feature edits, Result; History empty); (3) New with the style of a drawing (choose mode; settings and style; no files, no Legend edits, no per-feature edits; to the other mode only settings both modes have); (4) later: New from a selected region (same file, a region; style inherited; length-dependent values adapted: feature height, label font, stroke widths, GC/skew widths, windows). | R1-1, R3-2 |
| M-B | An old Session (27-44) opens with drawing A (the saved Result's mode) and drawing B (the other mode) only when the file shows that mode was used (inputs of that mode, a changed setting of that mode, saved while that mode was shown). No saved value is lost. | R1-2 |
| M-L | Old Legend edits go only to the drawing with the saved Result (drawing A). This replaces the per-mode answer "copy Legend edits to both". Other shared values are still copied to every created drawing. | R1-3 |
| M-Py | Opening an old-format Session may start Python once (about 8 s observed) to convert it; the saved picture shows at once; editing waits. New-format Sessions and the Gallery open without Python. PD-OI-044 gets scenario revision 2 for old files only. | R1-4 |
| I-1 | Input file bytes are stored once per project (content-addressed); which files, records and regions a drawing draws is drawing state. Another drawing can pick a file from "Session files" without uploading it again. | R2-1 |
| I-2 | Replacing a file that other drawings use asks: "this drawing only / all N drawings that use it". | R2-2 |
| CLI | New `gbdraw render --session FILE [--drawing ID ...] [--list_drawings]`; `gbdraw circular|linear --session` keep working when the file has one drawing of that mode; drawing IDs come from the kind (`circular`, `linear`, `linear-2`, fixed at creation); several drawings rendered together get `_<id>` in output names. | R2-3 |
| K | PD-OI-086 (and PD-OI-044 rev 2) enter the Contract through a Contract-only PR first (Review REQUIRED, no auto-merge), before the implementation. | R2-4 |
| U-1 | A drawing bar under the header spans the settings pane and the Preview. A fresh Session is one Circular drawing plus a dashed "Linear" offer tab that creates a blank Linear drawing on first click. The Circular/Linear switch goes. Empty drawings are not saved. | R3-1 |
| U-2 | Delete asks for confirmation when the drawing has a Result or edits; the header Undo covers the active drawing's steps and drawing-list steps (add, duplicate, rename, move, delete), newest first; activating a drawing is navigation, not a step. | R3-3 |
| U-3 | Reset Settings resets the active drawing, with an "Also reset the other N drawings" checkbox (default off); files, record selection, Depth files and the Result stay. | R3-4 |
| PRI | **First runtime deliverable of v0.15.0: each mode/drawing keeps its own Result.** A Generate never replaces another drawing's Result; switching back shows that drawing's own Result with its drags and edits. | Owner via Phase E, 2026-10-07 (「すごくイヤだねこのbehavior. なくしたいね」) |
| SAVE | **A Session file saves every Result.** Per-mode Results never ship with a Save that drops a visible Result; they ship together with a format that stores every drawing's Result. | R4-1 (「それは完全にregressionじゃん。ダメだよこれは。Session fileはすべてのResultを保存して。」) |
| OV-104 | The OV-104 stopgap was not to be merged while E1 was a day away (R4-2); with the move to v0.15.0 Phase E re-asks the Owner for 0.14.0. | R4-2, R5-1 |
| REGION | "New drawing from selected region" (with length-adapted style) is in the first drawings UI release, not later. | R4-3 |
| REGION-1 | A region drawing carries the source's look for that record: settings, palette and rules (sizes adapted), Legend edits, per-feature edits and annotations inside the region, that record's Depth; uploaded comparison tables are not copied. | R7-1 (`design-parts/region.md` 2.6) |
| REGION-2 | A hand-set size is kept only where its Auto is unchanged at the new length class, tier and mode; otherwise it returns to Auto and the dialog lists it ("Label font size 10 → 24"); "Adapt sizes to the region length" (on by default) can be turned off. Auto values follow the drawing's own length. | R7-2 (`region.md` 2.4 B) |
| REGION-3 | No origin-spanning regions in the first release: the dialog explains and disables Create; margins stop at 1 / record end. | R7-3 (`region.md` 2.3 R-a) |
| REL | The whole project ships in **v0.15.0**; this session designs. | R5-1 |

Earlier per-mode answers re-read for drawings: Q1 (new mode starts with defaults) = C-new. Q2 (a
Legend edit changes the Result on screen) holds trivially: a drawing shows only its own Result. Q3, Q6
(all Depth, GC and Legend settings separate) are subsumed by P. Q4 (Depth files per mode) = I-1 (bytes
shared, assignment per drawing). Q5 (Show Depth on only where that mode has a Depth file) is the
migration rule in 3.3. "Old Sessions copy shared values to both" holds for settings; M-L overrides it for
Legend edits. Q7 (Contract in the implementation PR) is replaced by K.

### 0.2 Decisions settled from code or precedent (not asked)

| Decision | Why |
| --- | --- |
| Drawings are Session **47** | 46 is taken by the 0.14.0 per-mode store (Phase E); a new number for a new shape (`contract.md` K3) |
| The 47 reader reads **46** (released in 0.14.0) and 27-44 through Phase E's split | one migrator chain in Python: 27-44 → 46 (Phase E's Python twin) → 47 (drawings); 45 is never read (dev-only, rejected since PR-1) |
| The Web reads **46 and 47 without Python** (JS fast path; Python twin + shared vectors); only 27-44 go through the Python migrator on Load (M-Py, R1-4) | every 0.14.0 Session would otherwise wait ~8 s for the Python Worker on its first v0.15.0 Load; 46 → 47 is a structural move whose only table (MODE-ONLY fields) is already generated for the Web (`mode-scoped-settings.generated.js`); D2 holds because Python owns the vectors (Phase E's twin pattern, #908) |
| Format name and extension stay (`gbdraw-session`, `.gbdraw-session.json[.gz]`); the word "Session" stays for the file and the project, "drawing" names its parts | `py.md` 2.1, `ux.md` U17 |
| A drawing's mode is fixed at creation | `ux.md` U3-a; removes every mode swap (D3) |
| Inside a drawing, today's field names stay (`suppress_gc`, `circular_track_slots`, ...); a drawing holds only its mode's fields | `runtime.md` 2.2: renaming churns the R10 guard and ~30 specs for no behaviour gain, and later would need another bump |
| Record display (crop, start, reverse) is per drawing (C4) | P; PD-OI-027's "record" reads as "the record entry of a drawing" (`contract.md` 1.3) |
| Depth series, including names, are per drawing; a series cannot exist without a source | Q3, P; `inputs.md` 3.4 |
| Annotation sets per drawing, with an explicit "Copy to drawing…" later | P; `inputs.md` 3.5 |
| Genetic code per record is drawing state (the record entry), not a project fact | P: no invisible cross-drawing effect; `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` rows marked project under P are overridden |
| Files no drawing references are not written on Save | `inputs.md` OQ-IN4; today's behaviour |
| LOSAT execution and thread settings are an app preference (`ui`); search settings and hit limits are drawing state | `contract.md` K6; PD-OI-017, OIC-015 |
| Drawing names unique per Session (case-insensitive) | `ux.md` U6: export names and CLI selectors use them |
| Export stays "the displayed Result of the active drawing"; "Export all drawings (ZIP)" later | `ux.md` U10, `runtime.md` RQ4 |
| Drawing switch disabled while Generate or a Session operation runs (as the mode switch today) | `runtime.md` RQ5 |
| LOSAT raw caches and the protein manifest are persisted **per drawing** in 47; a content-keyed project cache is a follow-up | their keys are request-local (`queryRecordInstanceKey`, `recordInstances`, `edgeKey`), `py.md` 1.10; `inputs.md` 3.6 (union manifest) is the follow-up once keys are content-derived (`recordAnalysisId`) |

## 1. Document (version 47, owned by Python)

```jsonc
{
  "format": "gbdraw-session",
  "version": 47,
  "createdAt": "2026-10-07T00:00:00+00:00",
  "title": "Vibrio",                                  // optional; the project (Session title)
  "resources": {                                       // project: every input byte once
    "genbank-1": { "kind": "genbank", "name": "genbank-1-Vnig.gbk", "originalName": "Vnig.gbk",
                   "type": "text/plain", "size": 23345678, "lastModified": 0,
                   "sha256": "9f2c…", "encoding": "base64", "data": "…" }
  },
  "ui": { "activeDrawingId": "linear", "losatExecution": { "...": "46 ui.losatExecution, unchanged" },
          "...": "46 app-level ui keys (downloadDpi, autoLabelReflow, paletteInstantPreviewEnabled)" },
  "drawings": [                                        // array order = tab order = render order
    {
      "id": "linear", "name": "Linear", "mode": "linear",
      "renderRequest": { "schema": 9, "mode": "linear", "...": "..." },   // or null: no Generate yet
      "results": [ { "name": "out", "content": "<svg …>" } ],
      "webFiles": { "bindings": { "...": "this drawing's input draft (resource IDs)" } },
      "config": { "form": {}, "adv": {}, "...": "this drawing's settings and editor intent" },
      "editorState": { "legend": {}, "featureCatalog": { "schema": 5, "items": [] }, "...": "" },
      "features": { "featureOverrides": {} },
      "ui": { "...": "this drawing's view state: 46 slice keys (layoutPreferences slot, canvasPadding, pendingPalette*, linearTypographyLinked) + zoom, pan, selectedResultIndex, generated*, appliedPalette*, featurePanelTab, input type" },
      "orthogroupState": {}, "runMetadata": {}, "cliInvocation": {},
      "losatCache": { "entries": [] }, "losatDerivedCache": { "entries": [] },
      "proteinIdentityManifest": { "schema": 2, "...": "" }
    }
  ]
}
```

Rules (Python validator; `py.md` 2.8, `inputs.md` 3.2-3.3):

- **Project fields**: `format`, `version`, `createdAt`, `title`, `resources`, `ui`, `drawings`. Nothing
  else. `drawings` is non-empty; IDs match `^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$` (≤ 40 chars) and are
  unique; names are non-empty and unique case-insensitively; `ui.activeDrawingId` names a drawing.
- **Resource descriptor**: `{kind, name, type, size, encoding, data}` required, plus `sha256` (required
  in 47), `originalName`, `lastModified`. `originalName` replaces `webFiles.resourceOriginalNames`. One
  validator owns the descriptor (today two disagree on `checksum`: OV-103). Every resource is referenced by
  some drawing's request, input draft or Result metadata; unreferenced resources are not written.
- **Drawing fields**: today's per-diagram top-level fields, unchanged in content, moved into the drawing
  (`py.md` 2.1 option A). `renderRequest.mode == drawing.mode`; every request resource reference
  (`resourceId`, `gffResourceId`, `fastaResourceId`) resolves. No request ⇒ `results == []`, catalog
  `null`, no `runMetadata` / `cliInvocation` (this replaces the settings-only variant).
- **Draft rows have no `scope`** (already so in 46 slices): the drawing is the scope. Placements, feature
  overrides and record-display drafts are keyed `[recordKey, biologicalFeatureId]` and validated against
  `drawing.mode`. The `scope` rows of 41-44 are split by Phase E's `by-scope`/`by-side` tokens (3.2).
- **Input draft** (`webFiles.bindings`) is per drawing and keeps today's slot shape for the drawing's
  mode, with resource IDs (`inputs.md` 3.3). Unifying Circular and Linear input shapes is not needed.
- **Depth**: one `depthSeries[]` per drawing, `{id, label, color, height, ticks…, sources: shared |
  perRecord}`; each series has a source; `byRecord` keys are this drawing's record keys (`inputs.md`
  3.4). This replaces the index-aligned `adv.depth_tracks` + source matrices.
- **Comparisons**: the Linear comparison plan and its file bindings live in the Linear drawing; Circular
  conservation series refer to resource IDs (not file signatures).
- **Unmanaged config overrides** are validated against the drawing's own mode (PD-OI-009; removes
  OV-106).

**Resource allocation** (one rule, Python-owned, mirrored in the Web by shared vectors; `inputs.md` 3.2,
`py.md` 2.4): the `SessionResourceTable` (`gbdraw/session_resources.py`, PY-C core, PR #910) is used while
encoding every drawing's request. Equal bytes (sha256 + size) held by the project or another request reuse
the ID; inside one request each input keeps a resource of its own (the planner reads each source file as
one genome, `losat_source_ids`); otherwise the preferred ID if free, else `<id>-<n>`; input files use `<kind>-<n>`;
request-generated tables use `<drawingId>-<today's id>`; a clashing materialization name becomes `X.2`
with ring labels pinned at encode time. IDs stay stable for the project's life, so no writer rewrites
references; only the 27-44 reader merges and rewrites once. PR #910 already replaced the Python codec
allocator and the attach dedupe; X replaces the Web save-time alias rewrite (`services/session-resources.js`)
with a JS mirror bound by vectors.

**One-drawing project vs today's canonical Session**: a lossless wrapper. For one drawing, the request
and resource bytes are identical to today's (PR #910 pins this for five fixtures). Guard: every Session fixture, upgraded, renders drawing A identically to replaying the old file.

## 2. Python API and CLI

From `py.md` 3.2-3.3 (adopted):

```python
class SessionDrawingSelectionError(SessionError): ...
@dataclass(frozen=True)
class SessionDrawing:           id: str; name: str; mode: Literal["circular", "linear"]; has_canonical_request: bool
@dataclass(frozen=True)
class SessionDrawingSpec:       request: DiagramRequest | None = None; mode=None; id=None; name=None; state=None
SessionDocument.drawings -> tuple[SessionDrawing, ...]; .drawing(selector=None, *, mode=None); .active_drawing_id
SessionDocument.mode / .has_canonical_request   # only drawing; raises SessionDrawingSelectionError when several
session_to_request(materialized, *, drawing=None) -> DiagramRequest
render_session(materialized, *, drawing=None) -> RequestRenderResult | CircularBatchRenderResult
render_session_drawings(materialized, *, drawings=None, output_prefix=None, formats=None, overwrite=None)
    -> dict[str, RequestRenderResult | CircularBatchRenderResult]
build_session_document(request=None, *, drawings=None, title=None, created_at=None, active_drawing=None, web_file_inventory=None)
save_session_document(path, request=None, *, drawings=None, ..., overwrite=False)
upgrade_session_document(document | mapping | path, *, temporary_directory=None) -> SessionDocument
```

- The beginner API (`gbdraw` root: `read_genbank`, `draw_circular`, `draw_linear`, `Diagram`) does not
  change; projects are an integration feature of `gbdraw.api`.
- `SessionDocument` is the only validator call; compat functions take a validated drawing view (drops
  three re-validations and a per-render deep copy).
- Rendering several drawings materializes each resource once and parses it once inside one
  `PreparedBiologicalInputCache` transaction whose retention is project-scoped with a byte budget; a
  counting test asserts one parse per resource.
- CLI (R2-3): `gbdraw render --session P [--drawing ID ...] [--list_drawings] [-o O] [-f F] [--overwrite]
  [--save_session | --session_output PATH]`; default: every drawing with a request; drawings without one
  are skipped with a notice (naming one is an error). `gbdraw circular|linear --session P [--drawing ID]`:
  without `--drawing`, the only drawing of that mode, else an error listing candidates. One drawing keeps
  today's output names; several get `<base>_<id>`; a Circular batch inside still appends `_<n>`. All
  target paths (and the sidecar) are preflighted together before any write. `--save_session` replaces
  only the rendered drawings' request, Results, catalog, run metadata and LOSAT artifacts; the document
  is validated **before** rendering (no partial output, the OV-102 symptom). Run Info replay becomes
  `gbdraw render --session <file> --drawing <id>`.

## 3. Migration to drawings (Python owns it)

**Chain.** One chain, owned by Python:
- **27-44 → 46** is Phase E's split (0.14.0). Its registry is `gbdraw/web_support/mode_scoped_settings.py`, with the Python twin of `services/mode-scoped-migration.js` and the vectors `mode-split-vectors.json`. It applies R1-1, R1-3 and `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` `level_P`.
- **46 → 47** is this design's step (3.1).
- **45** is never read: it was dev-only and is rejected since PR-1.

Readers:
- `SessionDocument` (lazy view) and `upgrade_session_document()` (CLI save path) run the whole chain.
- The Web reads 46 and 47 without Python (0.2): the 46 → 47 step has a JS twin with shared vectors (`tests/fixtures/sessions/drawings-vectors.json`, read by pytest and a node test).
- 27-44 go through the Worker helper `migrateSessionDocument` (M-Py, R1-4). PD-OI-044 revision 2 covers only these.
- The Web's JS migrations for 27-44 are deleted in X once PY-F has proven parity: Phase E's split, `services/feature-edit-migration.js`, the legacy paths of `services/config.js`, and `services/gallery-session-migration.js`.

Evidence of the old contracts:
- first-parent `main` `fe6861f0` writes 44;
- tag `0.13.0` writes 30;
- the `0.14.0` tag, once cut, writes 46.

`contract.md` 2.1 has the version table.

### 3.1 Steps 46 → 47

1. **Committed sets.** Mode m's set is the top level when `renderRequest.mode == m`, otherwise `otherModeResult` when its `mode == m`. A settings-only 46 (`renderRequest: null`) has no set, and `otherModeResult` is forbidden there.
2. **Drawing A** is the mode of the top-level `renderRequest`; for a settings-only Session, `ui.mode`.
3. **Drawing B** (the other mode, R1-2) is created when any of these holds:
   - a. B has a committed set (`otherModeResult`);
   - b. `ui.mode == B` (the Session was saved while B was shown);
   - c. a non-null binding in B's input slots of `webFiles.bindings`. Depth alone counts; the placeholder empty `linearSeqs` row does not;
   - d. a registry row of `modes[B]` whose value differs both from B's default and from `modes[A]`'s value for the same row. A MODE-ONLY row is compared with the default only. The comparison is canonical JSON, per row of Phase E's registry. A keyed-row domain (registry `key`: placements, record display, per-feature edits, annotation sets) is compared entry by entry: an entry of B counts only when A has no equal entry under the same key, so the Main placement rows that `by-side` writes to both slices never count. The wording is the same as Phase E plan §4.6, so the shared vectors agree.

   Why d is not just "non-default":
   - Phase E's `copy` tokens write a 27-44 file's shared values into both slices.
   - A plain "non-default" test would therefore give every old Session that was re-saved in 0.14.0 an unused drawing B.
   - A value equal to A's carries nothing that A does not keep.

   Residual: a 0.14.0 user who styled B exactly like A, with no Result, input or display in B, gets no drawing B. "New with the style of A" recreates it. On a 46 produced by the split, rule d reproduces each 27-44 signal of 3.2:
   - `own` puts class-M values only in their mode;
   - `profile` puts inactive profile values in B;
   - `layout`, `by-scope` and `by-side` (lane rows) write only B;
   - Main placement rows, `copy`, shared `by-leaf` leaves and `by-binding` rows write both slices equally.

   One divergence: a `managed: false` profile value equal to A's flat value no longer counts as use.
4. **Drawing m** is built from two parts, and the other mode's MODE-ONLY fields are trimmed (`runtime.md` K, registry `modes`):
   - `modes[m]`: `config`, `features`, `editorState.legend` + `featureStrokes`, and the slice's `ui` keys;
   - m's committed set: request, Results, catalog, `runMetadata`, `cliInvocation`, `orthogroupState`, `legacyArtifacts`, the Result-bound Legend refs (`originalOrder`, `originalColors`, `originalSvgStroke`), the alignment receipt, and the view keys (`generated*`, `appliedPalette*`, `selectedResultIndex`; zoom and pan for A).
5. **Input drafts.** `webFiles.bindings` is split by slot mode: the Circular slots (`c_*`) go to the Circular drawing, the Linear slots (`linearSeqs`, comparison and Depth bindings) to the Linear drawing. Resource IDs are kept (`inputs.md` 3.3). Resources that no drawing references are dropped.
6. **Resources** move to the project table. IDs are kept, `sha256` is computed, equal bytes are merged (the only reference rewrite in the design), and `originalName` comes from `webFiles.resourceOriginalNames`.
7. **LOSAT raw caches, the derived cache and the protein manifest** (top level in 46) are copied to each created drawing. This is lossless (content-keyed entries), and each drawing's next Generate prunes them.
8. **Project.** It keeps `title`. `ui.activeDrawingId` is the drawing of `ui.mode`, which always exists by 3b. The app-level `ui` keys (`losatExecution`, `downloadDpi`, `autoLabelReflow`, `paletteInstantPreviewEnabled`) go to the project `ui`.
9. **IDs and names.** IDs are `circular` and `linear`; names are "Circular" and "Linear". The array order is A, then B.

### 3.2 27-44: what the chain must produce (composition oracle)

Phase E's split produces 46, and 3.1 then makes drawings. The composed result must equal the rules below, which were revision 3's direct 27-44 rules (`MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` `migration_from_v44`). A composition test asserts it over F1-F3 and every single-mode fixture. The negative control: a single-mode fixture gives one drawing.

Two rules come from Phase E, as stricter forms of revision 3:
- Unmanaged overrides go `by-leaf`. A no longer keeps the other mode's leaves (PD-OI-009, OV-106).
- Depth is copied, then reconciled per mode. This equals the zip with the mode's source column; the vectors check it.

1. Normalizations come first:
   - today's normalizations: catalog expansion, `migrate_persisted_web_state_field_names` (with #908's placement-draft step), the legacy Linear comparison draft, and the payload adapters;
   - every request is promoted to schema 9 (the ≤42 frame and <44 legacy alignment steps, PD-OI-030, PD-OI-073).
2. **A** is the committed mode: `renderRequest.mode`; for settings-only 42/44, `ui.mode`; for 27-30, `session_mode()`.
3. **B** exists when any of these exists for the other mode (`contract.md` K2, `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.md` §1, `inputs.md` 3.9.3):
   - a non-null binding in that mode's slots (Depth alone counts; the placeholder `linearSeqs` row does not);
   - a `modeProfiles.profiles[other]` field with `managed: false`;
   - a non-default value of a class-M field;
   - a non-default `ui.layoutPreferences` slot;
   - a non-default track-slot stack, comparison plan, conservation series or LOSATP hit limit;
   - a `recordDisplayDrafts` row with that `scope`;
   - a lane placement on that mode's side;
   - `ui.mode` = that mode.
4. **Values** follow the tokens:
   - own-mode values go to their mode;
   - both: shared values are copied to every created drawing;
   - active (A only): the Result set and the Legend editor state (R1-3), decoration drags, and rendered-ID per-feature edits;
   - the flat plot title and fonts of ≤42 and early 44 go to A only;
   - placement drafts: a lane row goes to its side's mode, a Main row to both;
   - Show Depth is on only in a drawing with a Depth file (Q5);
   - annotation targets go to the drawing of their scope.
5. **Resources** are merged and hashed as in 3.1.6. In 27-30, embedded files become resources.

### 3.3 Fixtures (evidence: Web "Save Session" bytes from first-parent main or a release tag)

- **Single-mode shapes:** 29 existing fixtures (`contract.md` 2.3).
- **E0** (pushed as `test/session-two-mode-fixtures` `2703664f`; PR-1 cherry-picks it). Phase E plan §4.7 gives the per-row oracle (expected 46 and 46 → 47 outcome) for F1-F3 and the two #908 lane fixtures; the 46 → 47 vectors use it: F1 → B Circular (rules b, c, d), F2 → B Linear (d), F3 → B Circular (d, E-value `1e-3`), neither lane fixture nor any single-mode fixture → no B.
  - F1 `two-mode-project.v44`: both modes used; the Linear Result was saved while Circular was shown; Legend edits; Depth and GC values; the same file in both modes; record drafts of both scopes. Its placement is Main, because a lane placement cannot be saved while Circular is shown on main (OV-118).
  - F2 `inactive-class-m.v44`: drawing B from settings alone.
  - F3 `two-mode-thresholds.v42`: 5-field profile and a flat title.
- **Other v44 fixtures:**
  - the #908 placement fixtures `feature-placements-{circular,linear}.v44`;
  - Phase E's new main-written v44 fixture (Legend edits, Depth, GC Percent, a "This feature only" fill, a record-bound annotation).
- **v46 fixtures** (written by 0.14.0 after promotion; added in X):
  - F4 `resaved-single-mode.v46`: a Circular-only v44 with a non-default palette, opened and re-saved in 0.14.0. It is the negative control for 3.1 rule 3d, which must give one drawing.
  - F5 `two-results.v46`: both Results (`otherModeResult`), saved while B was shown.

## 4. Web runtime

From `runtime.md` (adopted; it rejects `mechanisms.md` 6.6's field accessor views with evidence).

- **Store (S-C)**: `project` (title, resource table, caches) + ordered `drawings` + `activeDrawingId`;
  each `DrawingState` holds today's drawing keys under today's names (105 keys: 25 draft settings, 37
  editor intent, 43 artifacts) plus the per-drawing state that lives outside `state.js` today
  (`committedCanonicalSession`, CLI helper files, match-sequence registry). Every `state` key is classified
  (drawing / project / ui / transient) by a guard.
- **Readers name a drawing**: services and asynchronous flows (request build, Generate, rerender, reflow,
  History, Session, Reset, export, Gallery) take an explicit `DrawingState`; UI owners resolve
  `state.activeDrawing()` once per action; the template binds setup-level aliases
  (`form: computed(() => activeDrawing.value.form)`); Generate pins its drawing at entry. An unconverted
  reader fails loudly (the flat keys leave `state`).
- **Switch** `selectDrawing(id)` is the only writer of `activeDrawingId`: refuses while busy; settles the
  departing drawing's open History intent; flushes the active Result's live DOM edits into its `results`;
  clears transient UI; sets the ID; the mount path binds the drawing's Result and projects its editor
  intent. No watcher writes state on a switch (`runtime.md` 3.3 disposes of every watcher).
- **History**: one manager per drawing behind a router with today's API, plus the drawing-list steps
  (U-2: add, duplicate, rename, move, delete) on a small project stack walked together with the active
  drawing's stack, newest first. One shared 200 MiB budget evicts inactive drawings' oldest steps first;
  file retention is the union over all stacks and bindings (fixes the latent per-manager release); each
  artifact handle is counted once (OV-112).
- **Results**: per drawing; the Preview shows only the active drawing's selected Result; the batch
  picker stays inside a drawing. `generatedMode`, the "retained Result" logic and every
  `generatedMode !== mode` guard disappear, and with them OV-83/OV-100, OV-104, OV-105, OV-107, OV-108,
  OV-113.
- **Mechanisms removed**: already in 0.14.0 by Phase E PR-1: the mode-profile manager and
  `config.modeProfiles`, R2 `scope` keys, `UNUSED_MODE_FRESH_FIELDS`, the Show Depth watcher and the
  profile transition. In v0.15.0: the `modes`/committed-set/`otherModeResult` top level (47), the
  `generatedMode` guards (W2), E1's per-mode artifact slots (folded into `DrawingState`), the mode toggle
  buttons (drawing bar), the per-mode discovery purge. Kept: `mode-profiles.generated.js` and
  `mode-scoped-settings.generated.js` (defaults and MODE-ONLY rows of a new drawing), `hitLimitsByMode`
  (LOSATP modes inside a Linear drawing).
- **Session Save/Load**: Save writes every non-empty drawing (one consistent document, PD-OI-045);
  Load builds the candidate project off-line and installs it by reference (the 40-domain
  snapshot/rollback goes); 46 and 47 files use a JS fast path bound to Python by a generated
  constants module and shared valid/invalid vectors (R4); 27-44 files go through the Python migrator in
  the Worker (M-Py).
- **Caches**: one Worker for all drawings; Generate is one at a time project-wide. Content tokens
  (sha256) and a retain set covering every drawing's bound inputs replace the "last request only" caches
  in the transport, the Worker and `PreparedBiologicalInputCache` (W5 / IN-A, measured first).
- **Memory**: N drawings each hold a live Result and catalog (Vibrio catalog ≈ 34 MB); X is gated on a
  measured Vibrio Circular + Linear project.

## 5. UX (from `ux.md` and the Owner's answers)

- **Drawing bar** (U-1) under the header, above both panes; tabs show the mode icon, the name and a
  hollow ring when there is no Result yet (no derived "needs Generate" status: PD-OI-037 rev 2). Fresh
  Session: `[O Circular] [: = Linear (+) :] [+ v]`. Phone: one picker row. Accessible names stay
  "Circular" and "Linear" for the first drawings, so tutorials ("Select Linear") and most of the 125 test
  selectors survive.
- **"+ New drawing"** (C-menu): New Circular, New Linear (blank); Duplicate; New with the style of
  "<active>" (choose mode); "From a record or region…" (first release, R4-3).
- **New drawing from a region** (`design-parts/region.md`; R7-1..R7-3):
  - Entry points:
    - "New drawing from this feature…" in the feature popup header;
    - "New drawing from selection…" in the selection toolbar (Ctrl/Shift+click, Shift+drag). A selection across several records gives a Linear drawing with one region per record;
    - "+ ▸ From a record or region…" with three choices: this feature, a whole record (an entry), or coordinates.
    - No base-pair drag on the Preview; annotation items and search multi-select come later.
  - Region rule: source coordinates, 1-based and inclusive. The region is the hull of the selected feature parts; on a circular record, the shortest arc. A 1,000 bp margin is added on each side and clamped to the record ends. A selection that crosses the origin is refused with an explanation (R7-3).
  - Mode: Linear by default; the source's mode for a whole record. Circular only for one record, because grid and batch keep no crop.
  - Carry-over (R7-1):
    - the source's look for that record, as listed in REGION-1;
    - LOSAT-sourced comparison edges between included records (rerun at Generate);
    - uploaded tables are dropped, and the dialog lists them.
  - Sizes (R7-2):
    - Python owns the rule in `gbdraw/auto_sizes.py`, published as `autoSizes` in `mode-profiles.generated.js` with shared vectors. The same function serves "New with the style of" across modes (N6).
    - The dialog lists each value that returns to Auto.
    - "Create and generate" is the primary button.
  - Known limits:
    - the 50 kb class cliff;
    - a 1 kb GC window on short regions;
    - features cut at the region edge are drawn with an arrow at the cut.

    All three are follow-up candidates (`region.md` §5).
- **Drawing menu** on a tab: Rename…, Duplicate, New with this style…, Move left/right, Delete
  drawing…; drag reorders. Names default to the mode ("Linear 2"), "<name> copy", "<feature> region".
  The last drawing cannot be deleted.
- **Empty state of a new drawing**: names the drawing and mode, says it starts from defaults, and offers
  "Choose its input files": Session files (N)…, Upload…; picking a file is an explicit copy (C-new).
- **Inputs** (I-1, I-2): each drawing keeps its Input Genomes panel; every uploader offers "Session
  files"; replacing a file used elsewhere asks "this drawing / all N drawings".
- **Reset** (U-3), **Delete and Undo** (U-2), **Export** (active drawing; ZIP later), **Save/Load**
  (all drawings; opens the drawing active at save), **Gallery** (opens the drawing the example names),
  **Generate** overlay names the drawing.
- Text changes: `ux.md` 3.3.

## 6. Level of every setting

`MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` column `level_P` is normative, with one override: the genetic code per record
is `drawing` (0.2). Counts under P: drawing 203, project 21 (resources and caches, Session title, drawing
order), ui 28 (18 view, 8 app, 2 session), retired 7, constant 11. Retire candidates and fields without a
clear home: `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.md` §6 (legacy Depth tick fields, `losatProgram` alias, flat `depth_color`,
`circular_definition_interval`, the Circular `pairwise_match_style` snapshot, `config.cliOptions.rawArgs`,
the four naming fields — Session title = project; drawing name = tab and default download name;
`prefix` = export override; plot title = diagram text).

## 7. Contract

- Phase E's per-mode store record lands first (0.14.0) and takes **PD-OI-086** (Contract rev 33; P0a
  sent to the Owner 2026-10-07). The drawings record of this design takes the **next free id**
  (PD-OI-087 or later, rev 34 or later), concern
  `session.project-drawings`, choice `A / PROJECT-WITH-DRAWINGS`; it extends the per-mode record from
  "per mode" to "per drawing" and states the per-drawing reading of every record that still holds
  (`contract.md` 1.3: about 60 records V/R).
- **PD-OI-044 scenario revision 2** (M-Py, R1-4), with the drawings record in v0.15.0 (Phase E: 0.14.0 needs
  none, its split is JavaScript): "preview-only Load の Python Worker 0" holds for Session 46 and later;
  27-44 Loads may start Python. Draft: `owner/receipt-draft-pd-oi-044-rev2.txt` (its "version 46 より前"
  wording is now exact).
- **OIC-028** (acceptance): edit isolation between drawings, a switch writes nothing, Save → fresh Load
  restores every drawing and its Result, each older fixture migrates with no saved value lost, CLI/API render
  one or all drawings with distinct paths.
- S00 decision 2's scope is superseded (already by the 0.14.0 per-mode record); its migration rule stays.
- R13 reaction introduced: replacing or removing a project input used by other drawings (I-2); owner = the
  drawing-list / input owner named in X; channel `port`.
- Route (K): P0a Contract-only PR (Review REQUIRED, no auto-merge) at the start of the v0.15.0 work, after
  the Owner approves the receipt; P0b `gbdraw/web/CLAUDE.md` wording (R2, R11, R14, Module ownership,
  app-shell, drawing-switch rule) separately; P0c `tools/web-change-policy.json` owner expansion if X adds an
  owner. All before X.
- **Relation to PD-OI-086** (Contract rev 33, merged into dev as #909 at `d474ac01`; receipt SHA-256 `ee7ef6a5…`). The rationale of
  PD-OI-086 says the drawings inherit its shape. The drawings record replaces four of its terms:
  - one History in which a mode switch is a step;
  - Reset Settings resets both modes;
  - 27-44 Load without Python;
  - the limit of two diagrams (Circular and Linear).

  Open for P0a, and the contract agent checks the precedent:
  - (a) **recommended**: a new concern `session.project-drawings`, with a `Supersedes` line that names those terms of PD-OI-086 scenario revision 1;
  - (b) PD-OI-086 scenario revision 2 under concern `web.mode.scoped-settings`.

  (a) keeps the project, input and CLI terms in their own concern.
- Receipt draft: `owner/receipt-draft.txt` (revision 2 of the draft: 46 → 47 chain, rule 3.1-3d, what
  PD-OI-086 already retired; not sent to the Owner yet; re-checked against the approved PD-OI-086 text
  before P0a).

## 8. PR sequence

### 8.1 v0.14.0 (bug-fix track; this campaign, beside Phase E's per-mode store)

Status on 2026-10-07 (dev `46067d77`).

| # | PR | Content | Files (main) | Status |
| --- | --- | --- | --- | --- |
| 1 | #908 OV-102 | Python placement-draft migration for 41-44 CLI replay; pre-render draft validation; shared vector; main-saved fixtures | `session_io.py`, `cli_utils/session.py`, tests | merged |
| 2 | #910 PY-C core | seeded Python resource table + reference walker; CLI re-save keeps IDs | `session_request_codec.py`, `cli_utils/session.py`, `session_io.py`, tests | merged |
| 3 | #913 R0 OV-119 | Depth TSV positions through the record coordinate map (crop window, reverse complement) | `gbdraw/analysis/depth.py`, tests, input-formats doc | merged |
| 4 | E1 + Q0 | per-mode artifact slots, `transitionDiagramMode`, switch refused during Generate, reflow or a History restore (OV-116), per-Result cache live identities, History restore port; `otherModeResult` writer/reader (Web) and validator/CLI (Python); one Load pipeline for both sets. Fixes OV-104, OV-107, OV-83/100, OV-108, OV-113, OV-115, OV-116, OV-117 by construction | ~13 Web files + Python session files, tests | review fixes (REVIEW-1) in progress; lands after Phase E's layering B (Owner); Review REQUIRED. Written as 45; PR-1 rewrites it as 46 |
| 5 | small fixes, each STANDARD | #918 OV-114 (CLI keeps v40-44 rendered-ID edits; Python twin of the JS migration, shared vectors); #919 OV-112 (History counts each artifact once); #920 OV-131 (Label TSV that applies to no label keeps the edits); OV-130 (actionable errors for an other-mode setting and a stale Reset alignment receipt); OV-132 (Python twin of the `hash=` → `featureIdentity` annotation-target move, shared vectors); OV-134 (own-key lookups in `feature-edit-migration.js`); OV-103 (resource `checksum` accepted and verified once in both Python loaders) | by finding | #918-#920 in CI; the rest open in the order #918 → OV-132 → OV-134, OV-103, ≤ 2-3 of ours in CI |
| — | E0 fixtures | F1-F3 (v44/v42 two-mode Sessions from main), branch `test/session-two-mode-fixtures` | fixtures only | pushed; for Phase E's PR-1 migration tests |

OV-106 and OV-135 are fixed in Phase E's PR-1. OV-133 moves to v0.15.0 with PY-F.

Phase E's track (for reference): OV-83/108 guard → layering A → (E1) → layering B → PR-0a/0b/0c (readers
name a mode; `state.drawings.{circular,linear}` + `DrawingState`, guard `tests/web/drawing-context.test.mjs`)
→ PR-1 (per-mode storage, Session 46, split + Python twin, Gallery regenerated) → layering C → D. Its
Contract-only P0a (PD-OI-086, rev 33) is with the Owner.

### 8.2 v0.15.0 (drawings)

| # | PR | Content | Notes |
| --- | --- | --- | --- |
| 1 | P0a | Contract: drawings record (PD-OI-087+), OIC-028, PD-OI-044 rev 2 if still needed | Owner approval, no auto-merge |
| 2 | P0b | `gbdraw/web/CLAUDE.md` wording | authority-only |
| 3 | PY-B | API: drawing views generalized beyond `circular`/`linear` (`SessionDrawing`, `SessionDrawingSpec`, `render_session_drawings`) | format-neutral |
| 4 | PY-D | `gbdraw render --session [--drawing] [--list_drawings]`, joint preflight, project-scoped parse cache with measurements | format-neutral |
| 5 | PY-F | port the remaining JS-only 27-44 migrations to Python with parity vectors (`feature-edit-migration.js`, legacy `config.js` paths, `gallery-session-migration.js`; Phase E's split already has its Python twin); the composition test of 3.2 | before X |
| 5b | R1 | `gbdraw/auto_sizes.py` (size class, window and tick tiers, Auto values per setting; replaces three inline window-tier copies in `api/diagram.py`), `autoSizes` in the generated profiles, Web Auto hints and `AUTO_SETTING_FIELDS` (replaces the hand copies in `auto-value-display.js`, `circular-track-slots.js`); reference SVGs unchanged | format-neutral; touches `session-request.js` → serialize with W1 |
| 6 | IN-A / W5 | content tokens (sha256) and a retain set for transport, Worker and `PreparedBiologicalInputCache`; LOSAT caches out of the artifact owner set | measured first |
| 7 | W1 | generalize `state.drawings` (PR-0a, keyed by mode) to drawing ids; fold E1's artifact slots into `DrawingState`; guard `drawing-context.test.mjs` | behaviour-neutral |
| 8 | W2 | History router per drawing + drawing-list stack; `selectDrawing` generalizes `transitionDiagramMode`; switch becomes navigation | behaviour-neutral for N = 2 |
| 9 | X | Session 47: `drawings[]`, per-drawing input drafts with resource IDs, `depthSeries`, the 46 → 47 step (Python + JS twin, vectors), Python for 27-44 on Load, removals (top-level `modes`/committed set/`otherModeResult`, `generatedMode`, the 27-44 JS migrations), Gallery and fixtures regenerated (+ F4, F5), memory gate | one format PR |
| 9b | R2 | `derive_region_drawing` / `RegionSelection` (Python API) + vectors | after PY-B / X |
| 10 | W3 + R3 | drawing bar, creation menu (New / Duplicate / New with style / From a record or region — R4-3), rename, move, delete + Undo, Reset dialog, Session files picker, Replace dialog (I-2); R3 = region dialog, popup and selection-toolbar entries, `services/region-drawing.js`, carry-over (R7-1), size rule (R7-2) | first drawings UI release |
| 11 | W4 | Export all (ZIP), Generate all, Gallery multi-drawing examples, tutorial disposition | later |

## 9. Findings

The campaign's findings list (OV-100 to OV-135) stays in the campaign folder. Status per finding:

- Fixed or fixing in 0.14.0 (this campaign): OV-102 (#908), OV-103, OV-112 (#919), OV-114 (#918), OV-119
  (#913), OV-130, OV-131 (#920), OV-132, OV-134; by E1: OV-104, OV-107, OV-113, OV-115, OV-116, OV-117.
- Phase E: OV-100 (= OV-83), OV-105, OV-108, OV-109; OV-106 and OV-135 in PR-1.
- Drawings (v0.15.0): OV-101 (Show Depth cleared on switch), OV-110 (unpruned composition deltas), OV-111
  (reflow outside the owner set; W1), OV-133 (CLI re-save of Sessions 31-39 drops rendered-ID edits; PY-F).
- OV-118: main only; already fixed on dev.

## 10. Risks

- **X size** (format + store + removals + Gallery): mitigated by landing W0, W1a-c, W2 and the Python
  PRs first (format-neutral), the field-name decision, and per-mode input shapes kept in X.
- **Migration completeness**: parity vectors before JS migrations are deleted; F1-F3; "no saved
  non-default leaf lost" test over all 29+ fixtures; SVG identity of drawing A vs old replay.
- **Old-file Load latency** (M-Py): measured in X; the migrator receives the document skeleton only, not
  resource bytes, Results or caches (PD-OI-045 performance).
- **Memory** with N drawings: shared History budget, live-artifact metric, gate on Vibrio C+L.
- **Promotion**: no dev → main promotion between the first PR that writes 47 and X's completion; the 47
  writer lands only in X (`py.md` 6). PY-B/PY-D/PY-F stay format-neutral (they read 46 and older).
- **Rule 3.1-3d**: comparison per registry row needs Phase E's registry to stay the single row list;
  a row added outside it would be ignored by the B signal. Guard: the 46 → 47 vectors enumerate every
  registry row once.
