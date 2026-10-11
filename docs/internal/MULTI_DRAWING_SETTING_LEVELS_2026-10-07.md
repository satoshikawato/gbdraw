# Part: levels — one level table for every setting and piece of state

Author: levels part agent, 2026-10-07. Table: `MULTI_DRAWING_SETTING_LEVELS_2026-10-07.tsv` (270 rows, 11 columns).
Generator and per-row decisions: `design-parts/levels/build_levels.py` in the campaign folder (outside the repository) (re-run to rebuild the TSV).

Inputs: `ov80/audit/form-adv.tsv` (140 rows) and `editor.tsv` (114 rows) at audit base `d61b2b0f`;
`design-parts/inputs.md`; `mechanisms.md` section 6. I did not re-read the audit evidence line by line
at `b355b8ae`; I re-checked only the points the table depends on (PD-OI-027 text, `recordDisplayDrafts`
scope in Session 44, which names exist in `fe6861f0`, the `createDefault*` defaults). Everything else is
"by code reading" through the audits. Nothing was run in a browser.

## 1. What the table says

### Columns and vocabulary

- `source`: `form-adv` (140 rows, all of them), `editor` (103 rows: the 114 audit rows minus 11 pointer
  rows already covered by a form-adv row), `inputs` (27 rows: input and resource state that neither audit
  lists as its own row, plus the new drawing fields).
- `today`: the audit classification, unchanged. `audit_level`: the audit's `level_proposal` text, unchanged
  (inputs rows: the editor audit level they refine, or `n/a`).
- `level_P`, `level_L`: `project`, `drawing`, `ui`, plus two non-state values: `retired` (a mechanism that
  disappears) and `constant` (generated code constant that seeds a drawing default).
- `ui` rows carry a tag in `notes`: `[ui:app]` app or machine preference (not saved per drawing),
  `[ui:view]` view or transient state of one drawing (never a setting), `[ui:session]` project navigation.
- `default_for_new_drawing`: a new empty drawing of mode M takes the mode-M default of
  `createDefaultForm/Adv` and `mode-profiles.generated.js`. Under Q1 it has no Legend or feature edits.
  Scale-dependent fields are stored as Auto (null).
- `on_duplicate`: same-mode Duplicate copies. The scale-dependent fields (`window_size`, `step_size`,
  `depth_window_size`, `depth_step_size`, `depth_min`, `depth_max`, `scale_interval`) reset to Auto only
  when the new drawing changes region or extent (Zoom to region), not on a plain Duplicate. Record-keyed
  state is remapped to the new drawing's record entries (fresh Linear uids). Result-derived state is
  `copied with the Result`.
- `migration_from_v44` tokens:
  - `active`: the drawing of the committed mode of the old Session (`renderRequest.mode`, else `ui.mode`),
    called drawing A.
  - `own-mode`: the drawing of the field's own mode. It is A if that mode is the committed mode, else
    drawing B. A split pair goes half to each drawing.
  - `both`: copied to every created drawing (Owner Q6 "copy shared values").
  - `project`, `ui`, `dropped`, each with a reason.
- Drawing B is created only on explicit evidence that the other mode was used. Proposed trigger (extends
  H3): a bound input in that mode's slots (incl. Depth only); a `modeProfiles.profiles[other]` field with
  `managed=false`; a row scoped to that mode (placement draft side, `recordDisplayDrafts.scope`);
  a non-default own-mode field of that mode, its track-slot stack, comparison plan or conservation series.
  If B is not created, own-mode values of the other mode were default or unreachable and are dropped.

### Duplicate as the other mode (not a table column)

A drawing's mode is fixed, so "Duplicate as Linear" creates a new drawing. Proposed rule from the `today`
column: `SHARED` fields copy, `PER-MODE` split pairs take the target mode default, `MODE-ONLY` fields of the
source mode do not exist in the target. Exceptions that take the target default although SHARED today
(different Auto algorithm or scale per mode): `label_font_size`, `axis_stroke_width`, `scale_interval`,
`show_scale`, `legend_font_size`, `legend_box_size`, `block_stroke_width`, `line_stroke_width`,
`gc_content_tick_font_size`, and all scale-dependent fields. Legend edits are not copied (Q1).
This is a UX choice for the Owner (question N6).

## 2. Counts per level (270 rows)

| Source | Rows | P: drawing | P: project | P: ui | P: retired | P: constant |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| form-adv | 140 | 121 | 0 | 2 | 6 | 11 |
| editor | 103 | 67 | 11 | 25 | 0 | 0 |
| inputs | 27 | 15 | 10 | 1 | 1 | 0 |
| **P total** | **270** | **203** | **21** | **28** | **7** | **11** |
| **L total** | **270** | **193** | **31** | **28** | **7** | **11** |

L moves exactly 10 rows from `drawing` to `project` (section 3); every other row has the same level.
Of the 28 `ui` rows: 18 `[ui:view]`, 8 `[ui:app]`, 2 `[ui:session]`.
`project` under P (21): resource table, sha256, originalName, the five kinds of input bytes (Circular,
Linear, Depth, conservation, comparison), colour-table bytes, label-filter bytes, LOSAT caches (raw, derived,
protein raw), protein manifest, record-list and parsed-source caches, GC series cache, record genetic codes
(two rows), Session title, drawing order.

Migration tokens in the table: `own-mode` 105, `both` 58, `dropped` 44, `active` 38, `project` 16, `ui` 7,
other 2 (`absent in v44`, `own-mode by leaf path`).

## 3. (1) Rows where P and L differ

Ten rows; they form three Owner questions. In L the project item is used by reference with a per-drawing
detach or override. In P it is a drawing item, and consistency comes from Duplicate and explicit
"Copy ... to other drawings" / "Apply to all drawings" actions.

| Row | P | L | Question |
| --- | --- | --- | --- |
| Palette selection | drawing | project | N1 style library |
| Custom palette colors (incl. `gc_content`, `skew_high`, `skew_low`: N3) | drawing | project | N1 |
| Specific color rules and captions (pattern rules only; hash rules stay drawing in both) | drawing | project | N1 |
| Qualifier priority rules | drawing | project | N1 (L has a "?" here: it changes label text) |
| Annotation sets and items | drawing | project | N-annot |
| Annotation styles | drawing | project | N-annot |
| Annotation target record | drawing | project | N-annot |
| Selected-feature annotation target | drawing | project | N-annot |
| Legend title (annotation set) | drawing | project | N-annot |
| History stack | drawing (+ project-operation stack) | project (one stack, drawing-tagged steps) | N-history |

Evidence behind each side:

- **Style library (N1).**
  - For L: a whole genome and its zoom should show one gene in one colour; rules classify the data and a rule
    matching nothing draws no row, so sharing cannot fail Generate (`editor.md` section 6, 8).
  - Against L: today per-figure edits are stored as rules: "This feature only" fill is a `hash` rule and a
    feature-type Legend rename is a set of hash rules (22 tRNA rules in the probe). A library would make
    them global. L therefore must split `config.rules` into pattern rules (library) and hash rules
    (drawing). It also needs a store, a reference and detach per drawing, and a visible "this drawing
    differs" marker.
  - P cost: colours drift between figures of one paper unless the user presses "Copy palette and rules to
    other drawings". P is additive-friendly: a library can be added later as data at the project boundary
    (`mechanisms.md` 6.10) without moving drawing items; removing a library later is harder.
  - Owner's rule "explicit overrides win; ask when invisible": under L, editing the library from drawing A
    changes drawings that are not on screen.
- **Annotation sets (N-annot).**
  - For L: a region is a fact about a sequence (`editor.md` section 8).
  - Against L: the target binding is by record key today, which fails in the other mode
    (`ANNOTATION_TARGET`, probe in `editor.md` section 5); marks, lanes, styles and the Legend caption are
    presentation (Q6). Feature targets survive copying because `biologicalFeatureId` is content-based
    (`inputs.md` section 1.4), so "Copy set to drawing" is cheap under P. `inputs.md` recommends P (OQ-IN5 A).
  - Either principle needs: a target names its record by source identity (resource + record id) and a
    drawing skips an annotation whose record it does not draw, with a notice.
- **History (N-history).**
  - P: matches D1 literally ("each drawing has its own History"); project operations (Duplicate, Delete,
    reorder, replace-in-all-drawings) need a second, small stack and a second Undo.
  - L: one stack. R11 already tags a step with its Result and Undo restores that Result while another is
    displayed (`commitResultEdit`). Undo first shows the drawing it changes. A drawing switch records no step.
  - History is not saved in the Session (today or later), so this has no format cost.

My lean: P for palette, rules, qualifier priority and annotation sets (add the explicit copy actions), and
decide History separately with the Owner (D1 favours P, usability favours L).

## 4. (2) Rows whose `audit_level` I changed, and why

`audit_level` in the TSV keeps the audit text; `level_P` and `level_L` hold mine. Rows where `level_P`
differs from the audit:

**A. The orchestrator principle P replaces the audit's project placement (L keeps the audit value):**
Palette selection, Custom palette colors, Specific color rules and captions, Qualifier priority rules,
Annotation sets and items, Annotation styles, Annotation target record (audit drawing; L project),
Selected-feature annotation target (same), Legend title (annotation set), History stack.

**B. Changed in both principles:**

1. **Similarity group names, descriptions, dormant names: project -> drawing.** Group ids come from one
   analysis run of one record set with one set of hit limits, so a name only means something for the
   drawing whose analysis produced the id. By code reading only (ids are run output; I did not prove they
   are not content hashes): verify before building. P rule: every edit belongs to one drawing.
2. **`depth_tracks` (and editor "Legend title (Depth series)"): drawing for the whole list.** The audit TSV
   says `drawing`, its section 7 says split (identity and label project, style drawing). `inputs.md`
   section 3.4 shows a project series cannot hold per-record sources and a shared label contradicts Q3/Q6.
   One list per drawing, `depthSeries[]` with its sources inside (inputs row).
3. **Audit "unclear" rows decided as drawing:**
   - `species`, `strain`, `circular_record_label`, `circular_record_subtitle`: drawing is the typed
     override; what the file says (organism, strain, inferred definition, `file_definition`) is a project
     cache of the parsed record (inputs row).
   - `circular_reverse`: drawing. PD-OI-027 says direction is "independent record state"; the record entry
     is drawing state, so there is no conflict (contract text re-read).
   - `window_size`, `step_size`, `depth_*`, `scale_interval`: drawing, stored as Auto (form-adv section 7).
   - `block_stroke_color`, `line_stroke_color`, `depth_color`, `features`, `feature_shapes`, `nt`,
     `min_bitscore`, `evalue`, `identity`, `alignment_length`: drawing. L may extend the library to the
     three colour fields; I kept them per drawing in both (they are settings, not palette).
4. **Relabelled, not moved:** 17 rows. `MODE_PROFILE_DATA.*` (11 rows) become `constant`: they seed
   defaults and are not state. `modeProfiles.profiles.*.values` and `.managed`, `losatProgram` alias,
   `depth_large_tick_interval`, `depth_small_tick_interval`, `depth_tick_font_size` become `retired`
   (6 rows; the 7th retired row is `webFiles.bindings`, inputs). `modeProfiles.schema/activeMode` is `ui`
   `[ui:session]` (audit already said ui).
5. **Audit rows split in two (bytes vs binding):** "Circular input files and type", "Linear File cards
   (files)", "Depth files", "Circular conservation files": the bytes stay `project`; which resource a
   drawing draws is a new inputs row at `drawing`. "Circular conservation settings" keeps `drawing`; the
   genetic code of a subject sequence is split out as `project` (fact about a sequence). Linear record
   labels: the typed override is `drawing`; the derived file facts are a `project` cache.
   "LOSAT caches": display rows (`losatCacheInfo`) are split to `drawing` (inputs.md section 3.6).
6. **Migration deviation (not a level change):** old Legend edits go to drawing A only, not to both
   (Q6 says "copy shared values to both"). They are Result-bound captions: copying them to drawing B,
   which has no Result, reintroduces OV-80 (a style on a row the drawing does not draw fails Generate).
   Owner question N7.

**C. Unchanged but contested:** LOSAT genetic code per record is `project` (a fact about a sequence, P
principle) while `inputs.md` section 3.3 keeps it in the drawing's record entry. Reconcile; an edit is
visible in every drawing that uses the record, so it should ask first.

## 5. (3) Per-mode mechanisms that disappear

Under both principles (a drawing has one fixed mode), with row counts from the table:

| Mechanism | Rows | What goes | Evidence (audit) |
| --- | ---: | --- | --- |
| Profile swap | 17 `profile-swap` (9 fields + 8 data) | `mode-profiles.js` state manager and transition; `config.modeProfiles` and `managed`; the mode watcher and the `semanticFileWatchersSuppressed` branch; `PROFILE_MANAGED_ADV_FIELDS`. Stays: `mode-profiles.generated.js` and `mode_profiles.py` as default source | `mode-profiles.js:7-26,268-294`; `app-setup.js:2106-2112` |
| Split pairs | 22 form/adv fields (`split-fields`) + 2 editor rows | one field per drawing: `show_gc`/`suppress_gc` (one polarity), `show_skew`/`suppress_skew`, label spacing, track layout, scale fonts; pairs with different domains (label placement, feature height vs width, GC/skew/Depth width vs height, `labels_mode`/`show_labels_linear`) become mode-only fields of their drawing | form-adv section 1, 4 |
| Split slot collections | 8 form/adv rows (+1 editor) | one `track_slots`, `_enabled`, `_axis_index`, `_schema_version` per drawing; a set rename rewrites one stack | form-adv rows `*_track_slots*` |
| R2 scope keys | 6 `scope-key` rows | `[scope, recordKey, featureId]` keys and `scope` in `recordDisplayDrafts`; the drawing is the scope; the "either mode" lane check; the Session 41-44 placement side rule only in the reader | `feature-placement.js:80-83`; `record-display-options.js:7` |
| `layoutPreferences` per mode | 2 form/adv rows (3 entries) | accessors `form.legend`, `adv.plot_title_position`, `ui.layoutPreferences`; Circular single vs multi becomes two fields or a grouping default (UX) | `state.js:277-300` |
| `generatedMode` guards | 7 `reset-on-switch` rows + the guards below | `generatedMode`, `resultCatalogFeatures` mode check, label-editor dormancy, rule-edit refusal, placement `targetsFor`, History track-visibility skip (OV-66), `UNUSED_MODE_FRESH_FIELDS`, the Linear track-toggle bug, mount key `svg-${mode}-...`, `matchSequenceRegistry.reset()` on switch | `state.js:748-765`; `watchers.js:328-368`; `label-actions.js:719-735,802-820`; `rule-actions.js:300-305`; `placement-actions.js:82-102`; `app-setup.js:2886-2890` |
| Cross-mode couplings | rows named in form-adv 5a-5j | Show Depth repair watcher, `label_rendering` rewrite, index-aligned `depth_tracks`, decoration-continuity mode clause, annotation target mode key | `app-setup.js:2099-2104`; `config.js:2045`; `session-request.js:1480-1490`; `decoration-continuity.js:98-106` |
| Duplicated controls | 8 depth display fields, `show_scale`, `scale_interval`, `axis_stroke_width`, `separate_strands` | the same field no longer appears in two mode cards | form-adv 5d |

Stays under both: `hitLimitsByMode` (keyed by LOSATP mode, PD-OI-002), `resolveComparisonThresholds`,
Python `mode_profiles.py`.

What differs by principle:

- **P adds** explicit actions: Duplicate (copies everything), "Copy palette and rules to other drawings",
  "Copy annotation set to drawing", "Apply to all drawings" for values that must be comparable (`depth_max`,
  GC percent range). One History stack per drawing plus a project-operation stack.
- **L adds** a project style library (store, reference, detach, "this drawing differs"), a pattern/hash split
  of `config.rules`, project annotation sets with per-drawing display/style override in the annotation track
  slot, and one project History stack with drawing-tagged steps. L removes the explicit copy actions for
  palette, rules and annotations.

## 6. (4) Fields with no clear home, or that look unused (retire candidates)

| Item | Why |
| --- | --- |
| `depth_large_tick_interval`, `depth_small_tick_interval`, `depth_tick_font_size` | No control; legacy fallback for series 1 (`config.js:918-924`). Marked `retired`; fold into `depthSeries[0]` on read. |
| `losatProgram` adv alias | Compatibility alias; the real Linear field is the editor row "LOSAT program". Marked `retired`. |
| `depth_color` | UI edits per-series colour; the flat value is series 1's fallback. Candidate to merge into `depthSeries[0].color`. |
| `circular_definition_interval` | No control and no Linear field though Python has `objects.definition.linear.interval`. Reachable only by import. Keep as mode-only, or add a control, or retire. |
| `pairwise_match_style` Circular snapshot | Profile-swapped but only Linear reads it; the request fallback (`ribbon`) differs from the Linear Web default (`curve`). Linear-only field. |
| `*_track_slots_schema_version` | Format marker of a stack; becomes the drawing schema version if Python owns the drawing schema. |
| `circular_grouping_intent` | Written as a side effect of `circular_record_selector` (single or auto). Candidate to derive. |
| `suppressCircularMultiRecordDefaults`, `depthTrackUiCounts.circular` | Runtime-only flags (not saved); derivable once defaults are applied at drawing creation and the series count is the list length. |
| `config.cliOptions.rawArgs` | Preserved CLI arguments the Web never reads; overlaps Python `cliInvocation`. |
| `rich_feature_popup` | Live-only viewer preference; `ui`, not diagram state. |
| `prefix`, drawing name, plot title, Session title | Four naming fields with unclear separation. Proposal: Session title = project; drawing name = tab and default download name; `prefix` = drawing export override; plot title = diagram text. |
| Genetic code per record | Home disputed between project (P principle) and the drawing record entry (`inputs.md`); see section 4C. |
| `show_scale` | One field, two meanings (Circular tick ring, Linear bar/ruler, and Linear managed axis colour). Fine per drawing, but the Linear axis colour coupling should be explicit. |

## 7. Owner question seeds produced by the table

- **N1 style library** (palette, custom colours incl. GC/skew colours N3, pattern colour rules, qualifier
  priority): P (per drawing + copy actions) or L (project library with detach). Recommendation: P.
- **N-annot annotation sets:** P (per drawing + "Copy set to drawing") or L (project sets, per-drawing
  display/style). Recommendation: P (agrees with `inputs.md` OQ-IN5 A).
- **N-history:** per-drawing stacks plus a project-operation stack (P), or one project stack with
  drawing-tagged steps (L). No strong recommendation; D1 is literal about P.
- **N6 Duplicate as the other mode:** what it copies (section 1 rule); recommendation: SHARED fields except
  the listed exceptions, no Legend edits.
- **N7 old Legend edits:** drawing A only (recommended) or both drawings (literal Q6).
- **N8 genetic code per record:** project fact (P principle) or part of the drawing's record entry
  (`inputs.md`). Recommendation: project, with a choice dialog on edit when other drawings use the record.

## 8. Risks and test implications

- The audits read `d61b2b0f`; main-line behaviour of Session 44 comes from `docs/SESSION_COMPATIBILITY.md`
  and `git show fe6861f0`. Per-feature edits in v44 are rendered-ID maps, not scoped rows (R2 scope is
  dev-only 45), so migration rows say `active`, not `own-mode`, for them.
- Rows marked `both` multiply: a v44 Session with an explicit shared value now stores it in two drawings.
  Tests: for every `both` field, a v44 fixture with a non-default value must show the value in both
  drawings; for every `own-mode` field, only in its own drawing. The table is the oracle for these tests.
- Rows marked `active` for Legend and drag edits depend on the A/B split rule. A v44 Session that opened
  with a Linear draft beside a Circular Result opens on B (`ui.mode`).
- Scale-dependent fields stored as Auto: the Web Auto hints must follow Python's per-drawing resolution
  (form-adv section 7), otherwise a zoomed drawing shows a stale hint.
- The 140 form-adv rows are all drawing level, so Python's drawing schema carries them; the Session
  format must keep `form`/`adv` objects per drawing (about 240 template bindings stay through the
  accessor-view plan, `mechanisms.md` 6.6).

## 9. PR implications

- The table is the input of the migration-vector fixtures (P1) and of the drawing draft schema:
  `drawing` rows define the draft shape, `project` rows the project block, `ui` rows the `ui` blocks,
  `retired` rows the code to delete in X, `constant` rows what stays generated.
- Can land before the format bump: nothing in this part; all level changes are format or state-model
  changes (X, W1). Authority wording (R2, S00 decision 2) follows the chosen principle.
- If L is chosen, add a project `styles` block and `annotationSets` project block to P1 and a library
  store to X; if P, neither is needed and the copy actions land in W3.
