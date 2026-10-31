# Session format history

Maintainer record of how gbdraw session files, canonical render requests, and
saved LOSAT results changed, version by version. Until 0.14.0 this file was
`docs/SESSION_COMPATIBILITY.md`; older plans and pull requests cite it by that
path and its line numbers. The public
[session and request compatibility reference](../REFERENCE/session-and-request-compatibility.md)
documents current support and what a user does with an older file, including
the [retired inputs](../REFERENCE/session-and-request-compatibility.md#retired-inputs).
Release notes record when a format changed. A Session version changes when a
breaking change enters `main`; a later breaking change made before that version
reaches `main` goes into its entry
([Persisted contract co-change](./WEB_CHANGE_POLICY.md#persisted-contract-co-change)).

## Unreleased: Session 50 stores a project and its list of drawings, the Linear wrap length, and italic gene names

Session version 50 replaces the two diagram-mode slices of Sessions 46 to 49
with a list of drawings, so a project can hold several drawings of one mode
(PD-OI-090), it stores the Linear wrap length in each Linear drawing, and it
stores italic gene names in each drawing. Session 49 reached `main` with 0.14.3
as the mode-scoped Session below, so all three changes take 50 (CLAUDE.md, "Persisted-format compatibility").

The project keeps `format`, `version`, `createdAt`, `title`, `resources`, `ui`
(the project's view keys and `ui.activeDrawingId`), and `legacyArtifacts`. Each
`drawings[]` entry holds one drawing: `id`, `name`, `mode`, `renderRequest` (or
null, with no `results`), `results`, `editorState`, `features`, `config`, `ui`,
`webFiles.bindings` (only the slots of its mode, and the five file cards),
`orthogroupState`, `runMetadata`, `cliInvocation`, `cliOptions`, `losatCache`,
`losatDerivedCache`, and `proteinIdentityManifest`. Array order is tab order. The
top-level `modes`, `otherModeResult`, `renderRequest`, `results`, `config`,
`features`, `files`, and `ui.mode` are retired. There is no settings-only form: a
drawing without a request has no Results.

The owner of the layout is `gbdraw/session_drawings.py`
(`read_session_drawings`, `write_session_drawings`); the admission band is
`SessionChecked49` in `gbdraw/session_contract.py`. Python and the Web app read
Session 50, and Save writes it. Reading a Session 46 to 49 moves each mode's
part into its own drawing (`circular`, `linear`), the rules of
`services/session-drawings-upgrade.js` in the Web app and of the 46-49 step of
`upgrade_session_document` in Python; Sessions 27 to 44 are converted by Python
once at Load. A Session 49's Circular `adv.show_definition` becomes the
Circular drawing's `config.adv.show_definition`. Fixtures:
`tests/fixtures/sessions/`.

The wrap length is the config leaf `canvas.linear.wrap_length` (`"none"`,
`"auto"`, or a positive length in bp) and a Linear drawing's
`config.form.linear_wrap_length`, an `own` row of the mode-scoped registry
(`gbdraw/web_support/mode_scoped_settings.py`) whose replay flag is
`--wrap_length`. Wrapped rows are drawn now: Python, the CLI and a loaded Session
draw a record longer than the length over several rows, so a writer can produce
a wrapped file.

Measured against `main`'s build (Session 49): it refuses any file of version 50
with "Session version 50 is newer than this gbdraw supports (49)", so it reads
neither the list nor a wrap length. The same build, given version 49, refuses
`modes.linear.config.form.linear_wrap_length` (the slice's `form` is closed,
`_validate_mode_slice_fields`) and a `configOverrides` entry
`canvas.linear.wrap_length` (`UNKNOWN_CONFIG_PATH`), and it accepts a CLI Session
with a stored `canvas.linear.wrap_length` leaf and silently draws the record on
one row (`from_dict` ignores the unknown leaf). That is why a stored length needs
the version.

A Session up to 49 has no length and loads with no wrap (`"none"`); no reader
converts it. A Web **Save Session** keeps a length loaded from a Session in the
request's `configOverrides`; the Web has no control for it yet. The Gallery and
test-input Sessions were refreshed to 50 (`tools/refresh_gallery_sessions.py`).
Positive fixtures of Session 49: `tests/fixtures/sessions/*.v49.*`.

Session 50 also adds the config leaf `labels.<mode>.italic_gene_labels` (true or
false, default false) in both modes and each drawing's key
`config.adv.italic_gene_labels` (either mode), a `copy` row of the mode-scoped registry next to
`label_protection`, whose replay flag is `--italic_gene_labels`. A Web writer
puts `labels.<mode>.italic_gene_labels: true` in the request's
`configOverrides` only when the setting is on, so a Session with it off has the
same request as before. The key reached no `main` Session before 50 and shares
its version with the drawings and the wrap length; a Session up to 49 has no key,
and reading it gives every drawing the packaged default (off); no reader converts
it.

## Unreleased: Session 49 stores whether a Circular slice shows the centre definition

Session version 49 adds `adv.show_definition` (boolean, Circular **Show
Definition**, PD-OI-109) to the Circular slice (`modes.circular.config.adv`),
an `own` row of the mode-scoped registry
(`gbdraw/web_support/mode_scoped_settings.py`). The mode-slice `adv` shape is
closed and Session 48 reached `main` with 0.14.0, whose reader refuses the key,
so the key needs a new version. The canonical request carries the setting as
the config override `objects.definition.circular.show` only when it is `false`
(the Web writer, `--hide_definition`, or a Python `config_overrides`); a
request that shows the definition is the request an earlier writer wrote, and
Gallery publication compares it unchanged. A Session up to 48 has no key and
loads with the definition shown (`trueUnlessFalse` in
`services/config-field-readers.js`); no reader converts it. Load projects a
request's `false` into the slice (`projectCanonicalSessionRequest`).

The checks that listed each Web-draft version by hand (the settings-only form
and Web binding schema 2, in `gbdraw/session_io.py` and
`services/session-authority.js`) now accept every version from 46 to the
current one, so Session 48 keeps its settings-only form. The Gallery and
test-input Sessions were refreshed to 49 (`tools/refresh_gallery_sessions.py`).
Positive fixtures of Session 48: `tests/fixtures/sessions/*.v48.*`.

## Unreleased: Collinear defaults and value sets

No Session or request version changes; the generated request description and
`LEGACY_REQUEST_TABLES` change (`tools/generate_web_contracts.py`).

- Merge orientation: the typed values are `strand` and `order`
  (`CollinearityMergeOrientation` in `gbdraw/analysis/collinearity.py`). The
  stored value set of `collinearityParams.parameters.mergeOrientation` keeps
  `either`, which released writers stored (it was the Web and Python default)
  and which drew as `strand`; its spellings read it as `strand` with no row
  (`LosslessCollinearityParameters.__post_init__`,
  `readCollinearMergeOrientation` in `services/losat-normalization.js`).
  Writers write `strand`, in every LOSATP mode: outside Collinear the Web
  writer writes the stored value as read. The schema 1-2 `standard` reading
  uses `strand`.
- Anchor mode: `one_to_one` and its spellings (`top1`, `mutual_best`, ...) are
  removed from `collinearityAnchorMode`. A stored one reads as `rbh`
  (`RETIRED_COLLINEARITY_ANCHOR_MODES`): Python logs a warning
  (`read_retired_anchor_mode`, applied before the schema 9 description check
  by `_as_decoded_settings`); the Web request repair takes the description
  default `rbh` with a `setting-replaced` row, and the Web draft reader
  (`readCollinearAnchorMode`) does the same. That description default is the
  repair target, not a writer's empty form: the Web repair replaces every
  unreadable stored anchor mode (an unknown word as well) by `rbh` with the
  row, where Python refuses an unknown word. For a retired spelling,
  `upgrade_session_document` writes the replacement in the requests and drafts of an
  older Session, so a Python upgrade and a CLI re-save hold `rbh` there. `cliInvocation.args`,
  derived entry identities, and other draft values stay as stored.
- Inference: the CLI and Python default of `collinear_infer_orthogroups`
  becomes off, as on the Web. A request without `collinearInferOrthogroups`
  still reads it as on (`_ABSENT_PROTEIN_SETTINGS`). A Session below 31
  replays its recorded command line with `--collinear_infer_orthogroups on`
  added when it names Collinear without the flag
  (`_keep_recorded_collinear_inference`), and a pre-canonical Web Session
  replays with the flag written from its draft.
- Search between: Collinear with inference off searches the displayed pairs
  whatever `collinearitySearchScope` holds (`effective_collinear_search_scope`
  in `gbdraw/comparisons/losat_jobs.py`, `effectiveCollinearSearchScope` in
  `services/linear-comparisons.js`), because without similarity groups the
  blocks read only those pairs. The stored value is unchanged and written as
  before. Raw cache keys hold no scope, so the raw entries of those pairs are
  reused; the derived identity keeps the stored scope and names the searches
  of the run, so a derived entry of an inference-off **All records** run is
  recomputed once from the raw cache. The drawn result is the same.
- Hit limits: the Web Collinear defaults of `hitLimitsByMode.collinear` become
  unbounded (`null`), as on the CLI; a draft that stores a limit keeps it.
- Derived entries: the default merge orientation text is now `strand`, so a
  derived entry written with `either` in its identity is recomputed once from
  the raw cache (token v4, not yet on `main`).

Positive fixtures: the Gallery Session requests now store `strand`, so the
request-level `either` comes from the `v0.14.3` file
(`tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json` at that tag) and from
the fixture below (the Gallery drafts still hold `either` in
`collinearMergeOrientation`); `tests/test_collinear_option_surfaces.py` and
`tests/web/stored-request-read.test.mjs` read the stored fixture
`tests/retired_sessions/BGC0000708-BGC0000713.v49-one-to-one.gbdraw-session.json.gz`.
Release tag `v0.14.3` offered `one_to_one` (`gbdraw/web/index.html`,
`<option value="one_to_one">`; the request literal listed it), so a released
writer could store it. No file in this repository was written by a 0.14.3 run
with it selected (users' saved files may hold it): the fixture is the Session 49 file of
that tag (`tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json`, request schema 9)
with both copies of the setting a user could choose written in, as the select stores
them: `renderRequest.comparisons[0].settings.collinearityAnchorMode` and the draft key
`modes.linear.config.losat.blastp.collinearAnchorMode` (the request derives from the
draft) (provenance in `BGC0000708-BGC0000713.v49-one-to-one.provenance.json`). It sits in
`tests/retired_sessions/`, outside the folders the scans of stored Sessions read, because
they refuse a value outside the description.

## Unreleased: derived LOSATP token v5 for a CDS across the origin in the Web

No Session, request, or derived envelope version changes. The Web converter
now receives each CDS's `origin_arc` and `record_length` and the record lengths
(OV-2368), so a Web result with a CDS across a circular record's origin has new
rows for the same inputs. The derived identity names the bindings, settings, and
raw keys, not the protein map, so a derived entry that 0.14.4 wrote would be a
hit with the old rows. The token becomes `derived-option-conformance-v5`:
earlier derived entries are cache misses and are recomputed once from the raw
cache at the next Generate; no search runs again and no reader reads them. The
protein map itself is never saved: it lives in the page's per-file memory
(`readOncePerFile`), and no Session field holds it, so `proteinMapSchema` and the
Session stay as they are. The comparison tables of a saved request keep their
rows until that Generate.

## Unreleased: derived LOSATP token v4 for chained Collinear blocks

No Session, request, or derived envelope version changes. Collinear blocks are
now the longest compatible chains of anchors (OV-2294..2296), and the block
rows no longer carry `collinearity_block_evalue` (OV-2300), so a derived entry
written before would replay the old blocks. The derived LOSATP identity token
becomes `derived-option-conformance-v4`: earlier derived entries are cache
misses and are recomputed once from the raw cache; no search runs again and no
reader reads them. Saved `collinearity-result` resources keep their shape:
`CollinearityBlock` still has `blockEvalue`, `label`, and the `syntenic` kind,
which the block caller no longer writes.

## Unreleased: LOSAT runtime records without `cli` or `bundled`

No Session or request version changes. The source tree no longer holds a LOSAT
build, so new CLI and Python writers no longer write the `cli` dialect field in
a `losatCache` entry's `runtime` record and cannot write `source: "bundled"`.
A new record is `{kind, version, source, path, program}` with `source` one of
`explicit`, `conda`, `managed`, or `path` (`LosatRuntimeRecord` in
`gbdraw/comparisons/losat_runtime.py`). Older records still load unchanged:
Python copies `runtime` as stored (`LosatRawCache._add_loaded_entry` and
`_add_search_entry`), and the Web run information reads only `source` and
`path` (`gbdraw/web/js/app/run-info.js`). The positive fixture is
`gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz`,
whose 12 records are `{source: "bundled", path: "gbdraw/bin/linux-x86_64/losat",
cli: "v1"}`; `tests/web/run-info.test.mjs` and `tests/web/losat-cache.test.mjs`
pin how the Web reads them.

## Unreleased: the palette `hearbal_garden` reads as `herbal_garden`

No Session or request version changes. Releases through v0.14.2 shipped the
packaged palette as `hearbal_garden` (`gbdraw/data/color_palettes.toml`). A
Session names a palette in `config.palette` (Web), `ui.appliedPaletteName` and
`ui.pendingPaletteName` (Web), and
`renderRequest.diagramOptions.colors.defaultColorsPalette`. The table is now
`herbal_garden`; `PALETTE_NAME_ALIASES` (`gbdraw/io/colors.py`, written to
`js/contracts/schema.generated.js` by `tools/generate_web_contracts.py`) maps
the former name to it, through `current_palette_name` in Python and
`currentPaletteName` in the Web.

- Python: `packaged_color_palettes()` looks a former name up and lists only
  current names (render, `require_palette_name`, the Session mode split).
  `ColorOptions` stores the current name, so every request (CLI, Python API, a
  decoded Session or request) draws it and a Session written from it names it.
  A CLI Session's `cliInvocation.args` keeps the command line as typed.
- Web: the palette slice (`services/drawing/palette.js`) stores the current
  name of the drawing's palette and of a queued palette, from Load, History,
  and every switch. `normalizePaletteDefinitions` adds the former name as a
  non-enumerable lookup key, so a raw-name lookup (the mode split of a
  Gallery publication, `services/gallery-session-publication.js`) finds it
  and the palette list shows only `herbal_garden`; a Session 27 to 44 is
  converted by Python once at Load (`drawings_document_from_legacy`). Load copies
  `ui.appliedPaletteName` and passes `renderRequest` through as saved
  (`services/config.js`), so a Web Save keeps the former name there until the
  next Generate; `app/results.js` compares the applied palette by its current
  name, and Python reads the request through `ColorOptions`.

Fixture: `tests/fixtures/sessions/HmmtDNA-hearbal-garden-cli.v48`, written by
v0.14.2.

## Unreleased: per-feature rows follow their record to a new key

No Session or request version changes, and no reader changes. A per-feature
row (`featureOverrides`, `featurePlacementOverrides`, `featureStrokeOverrides`,
an annotation's `featureIdentity`, and a `hash` color rule
`<feature id>@<record key>`) is still keyed by the request record key of the
view it was last drawn in: Circular `record-<k>` on the Multi-Record Canvas and
`circular-<record id>-#<k>` for a record drawn alone, Linear `<uid>:<k>` for a
File card that draws several records and `<uid>` for one record. Each
successful Web Generate moves the rows of a record whose key changed to its new
key (`recordKeyMoves` in `js/services/record-draw-selection.js`; the rows are
written by `pruneUnmatchedFeatureOverrides`, the `hash` rules by the
Generate's rule commit), and the rows of a record that the request does not
draw wait under its identity (`record-<k>`, or `<uid>:<k>` in a multi-record
File). Each Generate also brings the rows a drawn record has under its other
keys (Circular: its key drawn alone or on the canvas; Linear: its
identity) into its key,
field by field: the row moved from its previous key, or else the row already
under its key, wins; of two `hash` rules that then name one feature instance,
that record's own rule stays and the other is dropped. A Session of 48 or earlier is read as before; its rows move
at its first Generate. Limit: edits made on a Linear card that draws every
record of its file (written by gbdraw 0.13.0) do not follow its records after
the first Generate splits the card into one card per record. The only request-key change: a Circular record drawn
alone no longer keeps a stored `record-<n>` key that names another record.
Fixture: `tests/fixtures/sessions/record-identity-circular-single.v48`, written
by `v0.14.2`.

## Unreleased: one comparison-match attribute in saved Results

No Session or request version changes. On the first-parent line of `main`,
from cc707a5be (#254, 2026-06-18, Session 27) until 7ec9ce7cc (#304,
2026-07-19 13:01, Session 32), the writers named a comparison match only
`data-gbdraw-pairwise-match-id`. That range includes gbdraw 0.12.0, 0.12.1 and
0.13.0 (Sessions 29 and 30), Session 31, and Session 32 from be13a5c03
(2026-07-17). From 7ec9ce7cc the writers added `data-gbdraw-match-id` with the same value, and now
they write only `data-gbdraw-match-id`. The legacy name is a row of
`LEGACY_SVG_ATTRIBUTE_RENAMES` (`gbdraw/render/interactive_svg.py`, written to
`js/contracts/schema.generated.js` by `tools/generate_web_contracts.py`). The
Web legacy Result admission (`admitLegacyImportedResults`, through
`normalizeSvgResultIds`) and `enrich_svg` rename it once; an element that
already has `data-gbdraw-match-id` keeps it. Every other reader reads only
`data-gbdraw-match-id`. A Session of 40 or later never holds a Result with only
the legacy name: such a Result can be saved only after a Generate, and the
Python upgrade drops Results without a feature catalog. Fixtures:
`tests/fixtures/sessions/BGC0000708-BGC0000713.v29`/`.v30` and
`cli-linear-protein.v30`; no fixture holds a Session 27, 28, 31 or 32 comparison, which takes
the same legacy admission (below Session 40).

## Unreleased: Session 48 stores protein results in the search frame, what each file card loaded, and the feature label style

Session version 48 stores the protein results of a reverse-complemented Linear
record in the record's search frame (the selected and cropped record on its
source strand), as it already stored comparison tables (R5, OV-1145): saved
LOSATP tables (`precomputedProteinComparison`), `orthogroupResult` members, and
`collinearityResult` anchors and block orientations. Sessions 46 and 47 stored
them after the reverse complement, in the drawn frame. Ranks, block IDs, and
orthogroup IDs are unchanged.

One owner, `gbdraw.api.record_planning.project_source_bound_results`, moves the
results between frames, item by item, by their tags: a table row by its
`*_view_feature_svg_id`, an orthogroup member by its span against its
feature's spans on both strands, and an anchor side by its view feature ID.
An anchor side without tags (Session 39 anchors have no `*_view_feature_svg_id`
and no `*_feature_index`) follows its record: on a reversed record
`x -> L + 1 - x` (`gbdraw.linear_comparison.reverse_unbound_endpoint_span`).
A block's orientation flips when exactly one side moves. The planner projects
the stored results into the drawn frame; every Python writer (CLI
`--session_output` and `--save_session`, `save_session_document`) and the Web
converter at Generate store the search frame.

Readers convert a Session up to 47 once (`DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION`,
`drawn_frame_protein_results_to_search_frame` in `gbdraw/api/session_compat.py`).
A replay draws as before. An item that matches neither strand of its record
is kept and drawn as stored; a writer of Session 48 refuses it with
`SessionConversionError`, naming the first item. `upgrade_session_document`
rewrites a Session 46 or 47 to 48 with each converted resource's bytes
replaced under its stored ID; a CLI re-save writes the same bytes.
`gbdraw/web_support/stored_protein_frames.py` is the Worker boundary of the
same conversion. Gallery publication refuses a Session 46 or 47 that holds
such results; `tools/refresh_gallery_sessions.py` rewrites the Gallery.
`BGC0000708-BGC0000713.v47.gbdraw-session.json.gz` (main's Gallery Session)
and `majanivirus_collinear-rc.v47.gbdraw-session.json.gz` (main's converter,
Collinear) are the positive fixtures, each with its replay SVG.

Main's Web writers of Session 39, and of Session 40 before `3bca0e8d9`
(2026-08-11), stored a reversed record's drawn-frame feature ID in a saved
LOSATP table's `*_feature_svg_id` and wrote no `*_view_feature_svg_id`
(OV-1183). A reader of a Session up to 40
(`DRAWN_FEATURE_ID_SESSION_MAX_VERSION`, `decode_session_request`) marks a
table that has rows and source feature indexes and IDs but no view ID columns
(`LinearComparison.drawn_feature_ids`); the
projection reads its stored IDs as views and takes each source ID from the
feature index, so the table converts like a Session 47 one; the projection
clears the mark. A Session writer moves a still-marked table to the search
frame first (`resolve_session_write_request`), and the request codec refuses
one. A table with view ID columns, and every table of a later Session, is read
strictly.
`BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz` and
`BGC0000708-BGC0000709-rc.v40-early.gbdraw-session.json.gz` (main's pairwise
writer at `cdd310134` and `7aad9e3e0`) are the positive fixtures.

The derived LOSATP identity token becomes `derived-option-conformance-v3`:
derived entries written before are cache misses and are recomputed once from
the raw cache, while the stored results are still reused, so no search runs
again.

Session 48 also stores what each Web file card loaded into a mode slice, so
that **Remove** on the card gives the settings as if the file had never been
loaded (OV-1280, D-19). One owner, `gbdraw/web/js/services/file-imports.js`,
writes and reads it; the card key (`t_color`, `whitelist`,
`qualifier_priority`, `d_color`, `blacklist`) is the identity:

- A row that a file added carries its card's mark (`markedRows`, `isCardRow`): a `rules`
  row keeps `fromFile: true`, the mark Sessions have stored since main's
  2026-01-14 writer (D-ORIGIN-1), set only by the rules-file import and the CLI
  `-t` reading: a built-in preset's rows carry no mark and are the user's
  (D-21); a `whitelist` or `qualifierPriorityRules`
  row carries `importedFrom` with its card key. A priority row whose hand
  order the file replaced keeps that order in `handOrder`. A row edited by
  hand after the load loses its mark and is the user's.
- `config.fileImports` (registered in `gbdraw/web_support/mode_scoped_settings.py`)
  keeps, for `d_color`, each feature type's file color and the hand color it
  replaced (`handColor`, `null` for the palette's), and for `blacklist`, the
  file's words and the text before the load. Remove takes a value back only
  while it still holds the file's value.

A Session up to 47 has no `importedFrom` and no `fileImports`: its whitelist
and priority rows read as hand rows, and Remove keeps them; its `fromFile`
rules rows are the Specific rules card's, as they were, including the rows of a
preset, which the earlier writer stored with `fromFile` too. A Session 48 without a
draft (CLI `--session_output`, the Python API) reads the rows its `-t`,
whitelist, and priority files, and its `-d` and blacklist files, put into the
request as that card's rows (`FILE_IMPORT_ORIGIN_SESSION_MIN_VERSION` in
`gbdraw/session_io.py`); a CLI re-save keeps the draft's marks.
`file-card-imports-circular.v48.gbdraw-session.json.gz` (the Web writer after
all five cards and hand edits, D-ORIGIN-2) is the positive fixture.

Session 48 also stores the feature label style of each diagram mode as
`adv.label_line` (`none`, `under`) and `adv.label_protection` (`none`, `halo`,
`box`) in each mode's slice (`modes.<m>.config.adv`), the `copy` row of the
mode-scoped registry. The mode-slice `adv` shape is closed, so these keys need
a new Session version; Session 48 has not reached `main` (whose reader, Session
47, refuses a Session 48 file), so they go into its entry. The canonical
request carries them as `labels.<mode>.line` and `labels.<mode>.protection`
only when they differ from the mode's packaged style, so a request of the
default style is the request an earlier writer wrote (Gallery publication
compares it unchanged). A Session without them loads with the packaged
defaults of each mode (Linear `under` and `box`, Circular `none` and `halo`,
`gbdraw.labels.label_style`); a stored value outside its list loads as the
mode's default. The config value `labels.stroke_width` is removed: a stored
config that holds it loads and draws the new leader width (0.07 em, at least
1 px); an override of it fails with the unknown-path error. No other field
changed.

## Unreleased: Session 47 keys the Legend edits by the row key

Session version 47 stores each mode slice's Legend color and stroke edits
(`editorState.legend.colorOverrides` and `strokeOverrides`) under the row key:
the caption gbdraw drew for the row, kept in the entry's `originalCaption`
(the entry's `caption` when it has none). Session 46 stored them under the
name the row shows (`caption`). Keep both in the Legend rename dialog keeps the
typed name (PD-OI-061 scenario revision 5), so two rows can show one name; the
shown name no longer identifies a row. A rename changes only the entry's
`caption`; its edits stay under the key. Within one Legend list the entries are
unique by row key, not by caption.

Session 47 also adds the diagram font of each mode slice,
`modes.<mode>.config.adv.font_family` (the **Font** section, PD-OI-098): the
CSS `font-family` list a drawing's text uses, `null` for the default
(Liberation Sans). A Session stores only the name, never the font file
(D-02). A Session older than 47 has no such field; when its draft holds the
font as the preserved override `objects.text.font_family` (a Python API
request that set it, saved again by the Web app), the Web Load moves that value into
the field. `font-family-preserved.v44.gbdraw-session.json.gz` and
`font-family-preserved.v46.gbdraw-session.json.gz` are the positive fixtures.
No other field changed.

The bytes of a Session 46 cannot be read as Session 47: after a rename of
GC skew to `Z` and of GC content to `GC skew`, a color stored under `GC skew`
belongs to the GC content row in Session 46 and would name the GC skew row in
Session 47. Every reader of a Session older than 47 therefore gives each row
key the value Session 46 applied to that row (the Web compile of first-parent
`main` 9839de8e: `draftLegendRows` and `compileDirectEditorMutationPlan` in
`app/candidate-render.js`). It reads the stores as Session 46 Load read them
(keys trimmed, a value without a usable color or width left out) and the rows
as Session 46 Load kept them (an entry with a caption and a usable color, the
first of each caption in its list). An edit reaches the row whose caption it
is stored under, a listed row before a dormant one and a deleted row as its
Restore shows it, and moves to that row's key (`moved`). An edit under a name
that no row shows reached Python's row of that name and stays under its key:
`unshown` when it is the key of a renamed row, else `unmatched`. Session 46
applied the color edits in their order, then the strokes in the order of the
color edits' keys and then their own, so of several edits that reach one row
the later wins and each earlier one is a `collision` (dropped); a stroke field
(`strokeColor`, `strokeWidth`) that no later edit sets keeps the earlier value,
whose `collision` then names the row's key as `toKey`.
Each reader returns these outcomes as data (`{store, fromKey, toKey, outcome}`,
`toKey` null when dropped); the Web Load computes them before it commits.
The Web Load (`sessionLegendEditsByRowKey`
in `services/mode-scoped-migration.js`, before any other reader), the CLI
re-save (`build_session_json`), the Python migration of an older Session
(`_split_session_adjunct`), and `upgrade_session_document` (a Session 46 keeps
its drawings) read them this way; `legend_edits_by_row_key` in
`gbdraw/session_io.py` is the Python twin of `legendEditsByRowKey`.
`tests/fixtures/sessions/legend-row-key-vectors.json` pins both, and
`legend-renamed-rows.v46.gbdraw-session.json.gz` is the positive Session 46
fixture written by first-parent `main` 9839de8e (#1003), which holds Session
46. The Legend inventory and generated colors (`originalOrder`,
`originalColors`) were already keyed by Python's key and are unchanged, as is
`addedCaptions` (rule row keys).

The multi-drawing format planned as "Session 47" in PD-OI-090 becomes Session
48 (Owner decision of 2026-10-13, quoted in the Product Contract revision
that records PD-OI-061 scenario revision 5).

## Unreleased: Session 46 keeps the settings and the Result of each diagram mode

In the Web app, Circular and Linear each keep their own settings, edits, and
Result. A setting or an edit made in one mode, and a Generate, stay in that
mode; **Reset Settings** resets both modes, and one Undo restores both.
Switching modes changes no setting and shows that mode's Result with its moves
and edits, or the empty Preview.

**Drafts.** Session version 46 keeps each mode's draft in `modes.circular` and
`modes.linear`. A slice holds the mode's `config` (the former `config` without
`modeProfiles`, `cliOptions`, `paletteInstantPreviewEnabled`,
`adv.rich_feature_popup`, and the LOSAT execution settings), its `features`
edits, its Legend and stroke edits (`editorState.legend` and
`editorState.featureStrokes`), and its `ui` values: the mode's layout slot
(`layoutPreferences`), Canvas padding, the pending palette, the Linear
typography link, and the Features-list record. A slice holds only these
fields, and may omit any of them: Load reads an omitted value as that mode's
default, after the committed request of that mode's Result. App-level settings
moved to the top level: `ui.losatExecution` (LOSAT run mode, thread budget,
threads per job, parallel workers), `ui.richFeaturePopup`, and `cliOptions`
(the options of the CLI command that wrote the Session, unchanged). Per-feature
edit rows, Feature placement rows, record display rows, and `featureIdentity`
annotation targets no longer name a mode (`scope`): the slice is the mode, and a
row is keyed by `recordKey` and `biologicalFeatureId`. A slice's
`editorState.legend.entries` lists the Legend rows of the mode's Result, then
the renamed rows its last Generate did not draw (for example with GC content off
or Show Depth off), each marked `"dormant": true`; the rename and the styles
under the new name apply again when a Generate draws the row. A Session 46 that holds
the top-level `config` or `features`, a former `editorState` or `ui` draft
field, or a field outside these slices is rejected
(`INPUT_INVALID {field: schema, reason: FIELDS}`). A settings-only Session keeps
the slice of its shown mode. The CLI and the Python API write no `modes`; a CLI
re-save keeps the Session's `modes` and `cliOptions`.

**Results.** The optional top-level `otherModeResult` holds the other mode's
Result set. In a Session the Web app writes, the top-level `renderRequest`,
`results`, `editorState.featureCatalog`, `runMetadata`, and `cliInvocation` hold
the shown mode's set when it has a Result, otherwise the other mode's set;
`otherModeResult` holds the remaining set only when it has a Result. Its fields
mirror the top-level ones: `renderRequest` (the other mode, the same request
schema), `results` (at least one), `editorState` (`featureCatalog`,
`alignmentResetReceipt`, and the Result's generated Legend order, colors, and
stroke defaults), `ui` (the selected Result and the generated legend, title,
and palette), `runMetadata`, and `cliInvocation`. Its request names resources
in the one top-level `resources` table. A settings-only Session
(`renderRequest: null`) cannot hold the field. **Load Session** shows the saved
mode (`ui.mode`) when that mode has a Result, otherwise the mode that has one;
the other mode's Result waits for the mode button. A Gallery Session shows one
mode, holds that mode's slice only, and has no `otherModeResult`.

**Older Sessions.** Loading a Session 27–44 splits its one draft into the two
slices: a shared setting goes to both modes, a setting of one mode to that
mode, and the Legend edits and per-feature edits to the mode of the saved
Result. The draft of the mode that was not shown when the Session was saved
(`config.modeProfiles`) gives that mode its title, font sizes, and comparison
thresholds; a Depth or Show Depth value goes to a mode only as far as that
mode's Depth files reach, and an annotation bound to one mode's record goes to
that mode. Session 45 was a development-only format and is rejected; load its
source files instead.

On the command line, `gbdraw circular --session` and `gbdraw linear --session`
render the set of their own mode, at the top level or in `otherModeResult`. A
re-save with `--save_session` or `--session_output` writes Session 46. It
replaces the subcommand's set with the new render where that set is, keeps the
other set, and splits an older Session's draft into `modes` as the Web app
splits it. It keeps `ui.mode`, so the Web app still opens on the saved mode. It
also changes these parts that both sets share:

- `losatCache` holds the entries that the render returns, and
  `losatDerivedCache` keeps the Session's own entries unchanged. The CLI builds
  no derived entries, so a Session saved by a fresh CLI render has an empty
  derived cache.
- `proteinIdentityManifest` and `legacyArtifacts` are replaced by the render's.
  A Circular render, one diagram or a batch, returns the Session's own manifest
  and `legacyArtifacts`, so the protein entries it keeps still resolve. With
  both drawings, `legacyArtifacts` come from the Linear one.
- A Linear LOSATP render tries only the pending legacy protein candidates of
  the record pairs its protein search plans. A promoted candidate moves to
  `losatCache`; an unmatched one stays pending, with a warning, so a later Web
  Generate can still promote it; a rejected one stays. Legacy derived evidence
  is dropped once the Session has current derived entries, or when the render
  promoted the last pending candidate (OV-1132, OV-1143). The Web Save follows
  the same rule, pinned by `tests/fixtures/legacy_protein_artifact_vectors.json`.
- When the render migrates legacy protein IDs (only Sessions 44 and older have
  them), the protein references are rewritten throughout the Session.
- Resources that neither set's request nor the Web files name are dropped.

## Unreleased: Web Load of 0.13.0 Sessions

Session version 46 is unchanged. Sessions 27–30, written by gbdraw 0.13.0 and
earlier, have no canonical `renderRequest`. The Web app loads such a Session
from its saved settings, with its saved preview; the CLI replays it, and the
typed-session bridge does not convert it. The Web writers of Sessions 27–33
saved every Circular Custom Track Slots row with `spacing: null`, a field the
current slots do not read. When the Session's Custom Track Slots are off, that
null carries no setting and Load drops it, so the 0.13.0 Gallery Sessions load.
When the slots of a Session 27–30 are on, or a row holds a `spacing` value or
another retired field (`strict`, `compress`, `reserve`, `placement`,
`inner_radius`, `outer_radius`, or `gap_after`), Load fails with a message that
names the field and the track row, and the previous Session stays loaded. The
current Custom Track Slots use `radius`, `width`, `inner_gap_px`,
`outer_gap_px`, `side`, and `z`.

## Unreleased: Web Load of Session 31–39 table rows

Session version 46 is unchanged. The Web writers of Sessions 31–39 stored the
Default colors, Label whitelist, and Qualifier priority tables without
normalizing their cell values, so a value typed with a tab became extra cells
and a value with a line break became a short row. Since the Web table readers
reject a row with the wrong number of cells, Load rejected such a Session. Load
now reads these three tables of a Session 31–39 as the current writer writes
them: the extra cells of a row join its last column with one space, and a row
without its required columns is dropped (`feature_type` and `qualifier` for
Label whitelist, `feature_type` and `priorities` for Qualifier priority, and
`feature_type` and a color for Default colors). The Load notice names each
table and line read this way. CLI replay and re-save (`--session`,
`gbdraw render --session`) and `upgrade_session_document` read them the same
way and log a warning that names the lines; a re-saved Session holds the table
as read. Sessions 40 and later are read as before, and table file imports
still reject such a row. When another table of a Session fails to load, the
message names the table, for example `Session table: Specific colors.`

## Unreleased: Session 46 and Web feature edits by source identity

Session version 46 stores the Web app's Feature visibility, Label visibility,
and label text edits in each mode's `features.featureOverrides`: one row per
original-source feature, keyed and named by `recordKey` and
`biologicalFeatureId`, with the same fields as request
`diagramOptions.featureOverrides` and the label's original text
(`labelSourceText`). Feature placement drafts (`config.featurePlacementOverrides`
of the slice) are keyed the same way. Both modes can use the same record key;
each mode's slice holds its own rows, so a request carries only the rows of its
own mode. The four rendered-ID maps (`featureVisibilityOverrides`,
`labelTextFeatureOverrides`, `labelTextFeatureOverrideSources`,
`labelVisibilityOverrides`) are rejected in Session 46. The Web app sends these
rows in the request, loads a request whose `featureOverrides` array is not
empty, and draws them in the live preview by identity. An annotation made from
selected features is saved with a `featureIdentity` target in the slice of the
mode it was made in. Loading a Session 40–44 moves an annotation
target with one `hash=` feature selector, as a selection made it, to a
`featureIdentity` target in the mode of the Session's diagram only when the
figure cannot change: the record the target binds is drawn without a crop,
reverse complement, or rotation, and the hash names exactly one feature of the
Session's saved feature catalog, in that record. Load reports how many targets
moved. Every other annotation target of an older Session, including one of a
Session without a saved catalog, loads unchanged and keeps its meaning: the
hash of the drawn feature. The CLI moves the same targets in
`config.annotationSets` when it replays a Session 40–44 with `--session_output`
or `--save_session`, and logs how many moved. As after a Web Load, the
written request keeps the `hash=` targets that drew the figure.

The feature catalog is schema 5. Each drawn feature records the hash,
location, and record location it was drawn with (`drawnSelector`), which live
rule matching uses. A schema 3 or 4 catalog reads as schema 5 with no selector
values. Until the next Generate, a feature whose rendered ID carries its source
hash was drawn with its source coordinates, so live matching and **Load Label
TSV** use its source values; on other features (cropped, reverse-complemented,
or rotated records) live matching leaves `location` and `record_location` rules
to Generate, and **Load Label TSV** declines a table with such rows and says
why.

Loading a Session 44, or a Web Session 31–33, moves each rendered-ID edit onto
its feature. In a Session 40–44, a rendered ID in the saved feature catalog
names its feature; otherwise its `_record_<n>` and `__instance_` suffixes are
removed and the edit moves only when exactly one feature remains. A Session
before 40 has no catalog that Load reads, so Load reads its GenBank sources
again with its crops and orientations and matches each rendered ID's drawn
hash and record position. Without readable sources it uses the saved feature
metadata, which names only features of records drawn without a crop or
reverse complement. Other edits are dropped and Load reports how many. If the
diagram runtime cannot start or fails while reading the sources, Load fails
and keeps the current Session, so no edit is dropped. An older
Feature visibility edit hid every feature with the same hash, such as each copy
of a duplicated record; it now applies only to the feature that was edited, and
Load reports how many edits the next Generate draws differently for this.
Each moved edit goes to the slice of the Session's diagram's mode. The CLI moves the
edits in the same way when it replays the Session with `--session_output` or
`--save_session`, as does `upgrade_session_document()`, and logs these counts:
a Session 40–44 through its saved catalog, a Session 31–39 through its GenBank
sources read again with the crops and orientations Load uses. When a source
cannot be read, the CLI uses the saved feature metadata; it does not read the
features of the saved SVG, which Load tries before the saved metadata, so an
edit that only those features name is dropped and counted. The replay of a
Session 27–30 still drops these edits. A Feature placement
draft of a Session 41–44 reached every request with its record key: a lane
placement goes to the slice of its side's mode, and a Main placement to both
slices. The CLI applies the same mapping when it replays such a Session with
`--session_output` or `--save_session`.

## Unreleased: request schema 9 and feature identity overrides

Session version 44 is unchanged. Canonical request schema 9 adds the required
`diagramOptions.featureOverrides` array and the `featureIdentity` annotation
target. Each override row names one original-source feature by `recordKey` and
`biologicalFeatureId` and sets its Feature visibility, Label visibility, or
label text (see
[Feature identity overrides](../REFERENCE/typed-requests.md#feature-identity-overrides)).
Schema-8 requests are read as requests without overrides, and current saves
write schema 9. The Web app writes an empty array and cannot yet load a request
whose array is not empty or that has a `featureIdentity` target; it reports
this instead of dropping the edits. An exact Feature placement whose identity
the source does not have no longer fails the render; it becomes a feature
identity notice.

## Unreleased: value checks, derived label maps, and definition spacing

Session version 44 and request schema 8 are unchanged. Replaying a Session, or
generating from it, applies the value checks shared with the CLI and Python API:
window, step, and depth window or step must be positive integers, the
dinucleotide must be two letters from `A`, `C`, `G`, `T`, and `U` (`U` counts
as `T`), font sizes must be greater than zero, and stroke widths must be zero
or greater. A Session holding another value fails with the field or setting
named instead of drawing an empty or flat track. The Web app sends these values
as typed and no longer replaces a rejected value with Auto or a default.

A Web Session older than Session 40 loads with its saved preview, but its Result
has no current feature metadata. **Save Session** then asks for one **Generate
Diagram** and offers it; after that Generate, Save writes Session 44. Loading
such a Session and saving it directly, as 0.13.0 allowed, is retired. Session
40 and later Sessions save unchanged.

CLI Sessions and Sessions saved from `main` can contain the label maps
`whitelist_map`, `priority_map`, and `label_override_rules` under
`labels.filtering`. They are compiled from the label tables during a render and
are not settings. Load no longer keeps them as a preserved
`labels.filtering.raw` setting, and attaching a label table recompiles them, so
Qualifier Priority and whitelist edits take effect.

Web Sessions no longer write `features.labelOverrideContextKey`. Readers ignore
it in Session 44 files that contain it; label edits no longer depend on which
Result was displayed when the Session was saved.

A Circular request that sets `objects.definition.circular.font_size` without
`objects.definition.circular.interval` uses the font size plus 2, truncated to
an integer, as the definition line interval. This restores the 0.13.0 Web and
CLI spacing. The rule is applied when the overrides are applied, so the stored
request is unchanged and replay writes the same overrides back.

## Unreleased: records turned off for drawing

Session version 46 and request schema 9 are unchanged. A slice of `modes` may
hold `config.recordsOff`, the list of record keys the Web app does not draw: a
Linear card uid, or a Circular source selector `#N`. A Session without it, and
`[]`, read as all records on. The Web app writes it only when a record is off; the render request
holds only the records that are on, so the request has no new field. A CLI
`--session` replay draws only the records that are on, and a CLI re-save keeps
`recordsOff` and the settings of the records that are off.

## Unreleased: default label of a precomputed similarity ring

Session version 46 and request schema 9 are unchanged. A Circular ring from a
BLAST table (`--conservation_blast`, Python `ComparisonRingTrackOptions(source=...)`)
without a label is now labelled with the table's file name without the last
extension, as in the web app. CLI and Python Sessions now store
`conservationLabels` for every precomputed ring, as Web Sessions already did.
A Session without them was written before this change and drew the full file
name, which its replay still draws. A Session 27–30 replay now takes that name
from the Session's file bindings instead of drawing the temporary copy's name
(`arg3-<file>`). Web Load of such a Session also shows the full file name; it
showed the name without the extension before.

## Unreleased: a scale interval of 0 or less

Session version 46 and request schema 9 are unchanged. `--scale_interval` and
`objects.scale.interval` in the Python option objects (`config` and
`config_overrides`) reject a scale interval of 0 or less, and the web app's
Scale Interval field starts at 1. A render request and a Session read a value
of 0 or less as the automatic interval, as every earlier writer drew it:
in `renderRequest.diagramOptions.config`, in `configOverrides` (including the
flat `scale_interval` that earlier Web Sessions wrote), and in a Session 27–30
`cliInvocation`. Web Load shows the stored value in the field, **Generate
Diagram** draws the automatic interval, and the Source recipe omits
`--scale_interval`. A value of 0 or less typed into the field behaves the same.

## Unreleased: CLI and Python LOSATN / TLOSATX results

Session version 44 and request schema 8 are unchanged. A Linear run with
`--losat losatn` or `--losat tlosatx` (Python `losat_search` with those
programs) saves what the web app saves: one `nucleotideBlast` resource per
compared record pair with the raw search-frame rows, and one schema 2
`losatCache` entry per pair with the web raw key and the non-key `runtime`
record. Entries that the web app searches carry `runtime` too, as
`{kind: "losat", source: "wasm", version: null, program}`. A CLI or Python
runtime records only its executable's name as `path` (for example `losat`,
`LOSAT`, or `blastn`), so the runtime record carries no local directory. The
runtime once bundled in a source checkout recorded
`gbdraw/bin/<platform>/losat` (see "LOSAT runtime records without `cli` or
`bundled`" above). A CLI Session's
`cliInvocation.args` still hold `--losat_bin` and the input and output paths
as typed. The saved request carries the resolved comparisons, not the search
intent, so replay needs no LOSAT runtime. A request that still carries the
search intent cannot be encoded; resolve or render it first.

## Unreleased: CLI LOSATP results in Sessions

Session version 47 and the request schema are unchanged. A Linear LOSATP run
that writes a Session (`--session_output`, `--save_session`, a Gallery
rebuild) stores the results it computed: the request's
`generatedProteinComparison` marker has mode `none`. Similarity groups stores
the `precomputedProteinComparison` tables and the `orthogroupResult`; Collinear
stores only the `collinearityResult` (no tables); Pairwise stores nothing, so
its first Generate still searches (OV-1068). The first Web Generate after Load
reuses what is stored instead of searching again.

## Unreleased: one LOSAT cache entry per raw key

Session version 44 is unchanged. The Web app saved a raw key once per row, so a
Session with two Circular rings of one sequence repeated the key and could not
be loaded. It now saves one `losatCache` entry per raw key, named by the first
row, as the CLI does. After a load, the cache list shows one row per raw key
until the next Generate. A version 39 or later Session that repeats a raw key
is rejected with an `INPUT_INVALID` Session diagnostic.

## Unreleased: comparison rows in the search frame

Session version 44 and request schema 8 are unchanged. Comparison rows stored
in a Session (`nucleotideBlast` resources, uploaded tables, generated LOSATN
text, and `linear_comparisons`) use the search frame of each record as the
Session persists it: the selected and cropped record, 1-based, on its source
strand. A record saved with `presentation.reverseComplement` or a reversed
region keeps its rows in the search frame, and the planner projects the
orientation when it draws. A stored row outside 1..L of its record stops the
replay with a `COMPARISON_INPUT` error whose reason is `SEARCH_FRAME`.

Sessions saved from `main` (version 42 and older) stored the rows of a
reverse-complemented Linear record after the reverse complement. Both readers
convert them once at Load with `x -> L + 1 - x`, where L is the length of the
selected and cropped record:

- CLI replay (`--session`) converts the rows of the adapted request.
- The Web app rewrites the stored table bytes, both uploaded tables and
  generated LOSATN text, so a Session saved after Load is current. A table
  downloaded from that Session therefore differs from the file first uploaded
  to `main`. Saved LOSAT raw cache entries were already in the search frame
  and are unchanged.

A CLI sidecar of `-b` with `--reverse_complement` or a reversed region embeds
the reversed record as a sequence without `reverseComplement`, so the writer
stores the `-b` table rewritten into that sequence's coordinates; replay draws
the ribbons of the original run. CLI sidecars written by `main` already stored
the rows that way and replay unchanged.

Feature-bound rows (rows with `*_feature_index` and `*_feature_svg_id`, such as
saved LOSATP rows) follow the same search-frame rule. Sessions written by the
`main` CLI or Python API (`--session_output`, `save_session_document`) can
instead hold a reversed record's feature-bound rows after the reverse
complement, with `*_view_feature_svg_id` columns naming the drawn features
(OV-399). The planner reads each row's frame from its view IDs, so these rows
draw once and unchanged. Until Session 48 a re-save kept them as stored and
Gallery `finalize` reported them as a difference from the Web Session; a
reader of a Session up to 47 converts them to the search frame (OV-1145).

Linear Sessions saved by the `main` Web app write
`orthogroupState.selectedOrthogroupAlignmentFeature: ""` when no alignment
target is selected. Readers now treat the empty string as no target instead of
rejecting the Session.

## Session 44: typed Similarity alignment display state

The current Session 44 writer uses canonical request schema 8 to move
Linear Similarity Group alignment out of generated-protein pipeline settings. The Linear request layout now owns
finite `x` and `y` base translations keyed by `recordKey`, plus an optional
resolved `SimilarityAlignmentPlan`. The plan records the exact reference feature,
one validated decision per displayed record, and its rationale. Record
presentation or region state owns orientation, and anchor centers are projected
from the resolved record display. The nested plan is schema 2. Its unreleased
schema-1 predecessor has no current reader; the Session remains version 44 and
the canonical request remains schema 8.

Current writers never emit `align_orthogroup_feature`,
`alignOrthogroupFeature`, or the former Session-only
`selectedOrthogroupAlignmentFeature` copy. A Python
`SimilarityAlignmentReference` is resolved before saving; the Session stores
only its resolved plan. Supported request schemas 1, 2, 5,
6, and 7 retain a version-bounded reader for the old string. The reader keeps
that value private
until the historical selection can be materialized from saved stable feature and
orthogroup metadata; a successful save writes only schema-8 plan and translation
state. Loading a saved preview remains Worker-lazy.

## Session 44: feature anchors and independent record transforms

Session 44 stores source feature anchor capability metadata in feature catalog
schema 4 and preserves per-record display start, absolute orientation, and
feature-anchor provenance in one editable record-display draft. Released
Session 42/catalog 3 documents remain readable. Exact single-part locations
are migrated conservatively; ambiguous compound locations require Generate
again before feature-based rotation. Development-only Session 43 is rejected.

Session 44 also stores the selected visibility mode for Linear **Accession** and
**Length / Coordinates** independently. Each field is `auto`, `show`, or `hide`.
Auto shows the field while every effective rendered row contains one record and
hides that field diagram-wide when any row contains two or more records.
Disabled Record Layout ignores dormant shared-row values, and changing row
placement does not rewrite the selected mode.

Session 44 mode profiles also hold the per-mode Plot Title, Plot Title font
size, and Definition font size. A flat value without a per-mode entry, as in
Session 42 and earlier Session 44 files, belongs to the active mode only; the
inactive mode starts from fresh defaults. In every accepted version a missing
`config.linearRecordLayout` means Arrange in rows is on, and a missing
`ui.linearTypographyLinked` means linked while both font sizes are equal.

Canonical request schema 8 retains the effective booleans introduced in schema 7.
Version-42 editable `true` values migrate to Show, `false` values migrate to
Hide, and missing booleans migrate to historical Show. A selected-mode field,
when present, takes precedence over the retired boolean. Current writers emit
only the selected fields. Saved Results and committed render requests are not
regenerated during Load; the migrated selection takes effect on the next
Generate.

## Session 42: settings before the first source

Web **Save Session** can preserve settings before any biological input is loaded.
Session 42 adds the explicit `renderRequest: null` document variant with empty
`results` and a null `editorState.featureCatalog`. The existing config, mode
profiles, editor preferences and file-binding owners retain the settings and
auxiliary resources. No empty sequence or synthetic request is persisted.

Both active and inactive biological inputs exclude this variant, as do committed
render artifacts. Missing requests and missing or dangling resources remain
errors. Loading settings-only replaces the old Session, including its source,
Result and committed request, through the existing atomic import transaction.
After loading a real source, Generate and History use the usual render path.

Request schema 7, catalog schema 3 and bindings schema 2 are unchanged. Bindings
schema 1 and supported sessions 27–33/39–41 retain their existing reader support;
session 41 still requires a canonical request. All new full Sessions also use
version 42. Readers whose maximum version is 41 reject these new files.

Python document loading and materialization accept settings-only documents.
CLI replay, typed request conversion and rendering require a biological request
and reject this variant explicitly; they never fall back to legacy argv.

## Supported versions

Current writers emit one session and request format:

| Format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 48 | 27–33, 39–42, 44, and 46–48 |
| Canonical `renderRequest` | 9 | 1, 2, 5, 6, 7, 8, and 9 |
| Web file bindings | 2 | 1; 2 in sessions 41–42, 44, and 46–48 |

Previously written Session 44 documents with request schema 7 or 8 and feature
catalog schema 4 remain readable. Current saves write request schema 9;
the saved preview is retained until the next Generate.

Session versions 34–38 and 45 and canonical request schemas 3–4 were
development-only formats. They were never released on the supported history and
are rejected.

The public typed-session bridge can convert full session versions 31–33, 39–42, and 44 to
a typed request. Versions 27–30 do not contain a canonical `renderRequest`, so
the bridge does not convert them: the CLI replays them with the same
`circular` or `linear` subcommand that created the session, and the Web app
loads them from their saved settings.

`render_session()` is the compatibility boundary for canonical session replay.
It migrates supported persisted artifacts into `CurrentRequestArtifacts`, then
calls the same current typed renderer used by fresh requests. Fresh
`render_request()` calls never parse old session envelopes or rewrite retired
protein identifiers.

Both `.gbdraw-session.json` and lossless gzip-compressed
`.gbdraw-session.json.gz` files are accepted. The Web app writes the compressed
form by default; the CLI writes uncompressed JSON unless the output name ends in
`.gz`.

## Current request ownership

Canonical request schemas 5 and 6 record Circular grouping explicitly as
`single`, `grid`, or `batch`. Schema 6 also records each input's cardinality;
schemas 1, 2, and 5 decode records as `exactly_one`. A single diagram or grid
has one output object. A Circular batch has one resolved output object per
record. `renderRequest.output.prefix` is the output-prefix owner.

The Web projects a selectorless Linear schema-5 card to explicit `all` when it
is saved with schemas 6–8. Legacy multi-record Web inputs already have explicit
selectors, so this preserves the embedded source records shown by the card.

Current sessions keep mode-specific layout values under
`ui.layoutPreferences`. Supported older parallel title, legend, and grouping
fields are migrated when read and are not written again.

The Web writer stores file bytes once under `resources`. `webFiles` binds those
resources to active and inactive input controls, so files shared by the
committed request and an editable draft are not copied into a second payload.
For Linear comparisons, `webFiles.bindings.linearComparisons` contains only a
stable comparison-edge ID and its file binding. Endpoint, source, inclusion,
and filename metadata are not duplicated there. Version 39 sessions with
legacy embedded `files` remain readable.

`renderRequest` owns the last committed render. Web config retains inactive
Custom Track stacks, disabled rows, draft Axis positions, and per-mode
comparison profiles. `editorState.featureCatalog` holds the schema-3 feature
catalog used by the saved preview and editor. Current sessions store one base
SVG Result per logical diagram; readers collapse paired base and interactive
Results from supported older sessions.

Linear comparison intent is an independent editable draft under
`config.linearComparisonPlan`. Its mode is `none`, `adjacent`, or `selected`;
it also stores the default adjacent source and stable edge metadata. Placement
alone remains under `config.linearRecordLayout`. Current writers do not store a
global `blastSource`, nested layout comparisons, per-record BLAST files, or
per-record LOSAT filenames. The editable plan can therefore opt out without
changing the last committed comparison artifacts in `renderRequest`.

Supported pre-40 Web sessions migrate directly to this plan. A disabled or
absent legacy layout retains the old adjacent LOSAT/upload behavior, an enabled
explicit list becomes `selected`, and an authoritative empty explicit list
becomes `none`. Legacy per-record uploads and custom filenames are attached to
their original positional gap by stable record UID. CLI-only replay sessions
do not gain a synthetic Web comparison draft: they load with **No comparison**.
A CLI-only session written with `--losat` keeps the adjacent
LOSATP comparison its CLI drew. The accepted session versions remain 27–33,
39–42, and 44.

## Retired inputs

The current table of retired CLI, Python, and table names is in the public
[reference](../REFERENCE/session-and-request-compatibility.md#retired-inputs).

The private `__gbdraw_legacy_spacing` key is read only from canonical request
schemas 1–2 and is never written by schemas 5–6. Pixel spacing migrates to
`inner_gap_px` and `outer_gap_px`. Factor-based spacing can be replayed but
cannot be saved losslessly in the current format; replace it with explicit pixel
gaps before saving a migrated session.

## Saved protein-comparison results

After a source replacement, the Web writer saves only protein raw cache entries
that resolve through the saved identity manifest. Valid entries for inactive
search settings remain eligible; saving does not clear the live cache or History.

A protein-search cache hit requires the same amino-acid sequences, selected
protein set, record-instance and feature bindings, query/subject direction,
program, and meaningful search arguments. Upload filenames, modification times,
resource names, and display-only aliases are not part of that identity.

Current sessions use these independent payload schemas:

| Saved payload | Schema |
|---|---:|
| Protein raw search result | 4 |
| Derived protein comparison | 3 |
| Protein identity manifest | 2 |
| Nucleotide raw search result | 2 |
| Canonical typed analysis resource | 3 |

Typed resource readers 1 and 2 preserve their saved ortholog path tuples as
explicit collections in the current model. Schema 3 writers store lossless DAGs
for newly inferred paths and explicit collections for supplied legacy corpora.
Counts are exact decimal strings. Request 8 and derived envelope 3 remain
unchanged; derived identity includes `pathRepresentation` to prevent
reusing an old analysis payload as a current helper result. Older releases that
only support typed schemas 1 and 2 cannot read the new typed resources.

Derived envelope 3 stores each displayed pair's comparison table as
`canonicalTsv`, the text of the Python canonical table encoder (D-66). Entries
written from #312 (17e2c9dee) until then carry `rows` instead: they still load.
In the Web, any derived result read from a Session (old `rows` entries
included) is recomputed once from the raw cache.
A zero-hit pair is accepted in either form, a header-only `canonicalTsv` or
`rows: []`, but not both. Session version and the envelope schema are unchanged.
A zero-hit entry (no hits) may also hold its mode's result member,
`orthogroupResult` or `collinearityResult`, as the Web writes it, when that
member decodes as the mode's typed result; it and every other member must hold
no values. Python Load and the Web Load admit the same entries (OV-1127), except one
known difference (OV-1149): an `orthogroupResult` whose typed `value.fields` is
empty is admitted by the Web Load, which checks only the envelope, and refused
by the Python Load, which decodes the member.

OV-1063 (D-90): one Python owner (`gbdraw.analysis.protein_artifacts`) builds a
derived entry's identity, key, provenance, and pair numbering for the Web
lookup, the Web converter, and the Python API. The identity's `semantics` token
was `derived-option-conformance-v2`, so entries keyed by the earlier
identities are cache misses and are recomputed; no reader or migrator reads
them. Session version 47 and derived envelope schema 3 were unchanged. Session
48 raised the token to `-v3` (search-frame protein results).

Generated protein FASTA, raw QUERY/SUBJECT fields, protein maps, and derived
references use deterministic session-internal handles of the form
`h_[a-z2-7]{26}`. Those handles bind a record instance to a complete CDS
identity. They are not user-facing protein names.

Sessions 27–33 may contain schema-2 protein candidates and schema-1 derived
evidence. Current readers keep those candidates separate from current cache
entries. Generation verifies the full FASTA content, program and arguments,
direction, and one-to-one feature mapping. A verified candidate is promoted
without rerunning LOSAT; an unverifiable candidate is a cache miss only for its
record pair.

**Save Raw LOSAT TSV** resolves generated protein handles through the identity
manifest immediately before download. It writes readable, percent-encoded
`protein_id`, `locus_tag`, GFF `ID`, or location-based aliases, while preserving
comments, row order, columns 3–12, numeric spelling, and line endings. The whole
download fails if any handle cannot be resolved. User-uploaded comparison TSV is
never rewritten.

## Typed-session resource lifetime

`materialize_session()` decodes embedded resources to temporary paths. Those
paths are valid only inside the active materialization context:

```python
from gbdraw.api import (
    load_session_document,
    materialize_session,
    render_session,
)

document = load_session_document("figure.gbdraw-session.json")
with materialize_session(document, output_directory="out") as materialized:
    result = render_session(materialized)
```

Do not retain a decoded resource path after the `with` block ends. Use
`with_request_output()` inside the same context when replay needs a different
prefix, output directory, format, or overwrite policy.

## Session 41 and request schema 7

The joint format adds requested record rotation, exact non-Auto feature placement,
and canvas overlap tolerance. Web config retains editable drafts separately from
the successful request and Result. Session 40 keeps its authority rules; schema 6
keeps cardinality and row inheritance. See the [current field and migration
contract](../REFERENCE/session-and-request-compatibility.md#record-rotation-and-feature-placement).

## Web binding schema 2

Session 41 writers assemble `webFiles.bindings` schema 2 while retaining request
schema 7. The binding inventory includes `c_gb`, which may be null, an ordinary
File binding, or one explicit composite. A composite stores `kind: "composite"`,
an ordered `components` list of at least two ordinary bindings, and the logical
File's `name`, `type` and `lastModified`. Ordinary bindings store `resourceId`
and those three metadata fields. Components refer to the existing `resources`
table; repeated references retain their positions. Nested composites and
composites outside `c_gb` are rejected.

Schema 1 remains readable with its existing metadata defaults and independent
File-array semantics. A validating serializer preserves a supplied schema-1
document; a newly assembled Web inventory uses schema 2. Sessions without an
explicit binding retain request-derived initialization. Schema 2 requires
session 41; older readers that support schema 1 reject schema 2 explicitly.

Canonical CLI replay with `--session_output` transports schema-1 and schema-2
draft bindings through the existing Web file inventory. The sidecar contains
regenerated committed inputs and Results alongside the original draft component
bytes, metadata and ordered occurrences, with destination-safe resource IDs.
Explicit binding slots override historical direct lists even when null or empty.
The frozen main-writer witness and its provenance are in
`tests/fixtures/sessions/single.v41-bindings1.json` and the adjacent provenance file.
