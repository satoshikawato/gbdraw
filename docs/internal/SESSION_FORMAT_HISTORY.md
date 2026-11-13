# Session format history

Maintainer record of how gbdraw session files, canonical render requests, and
saved LOSAT results changed, version by version. Until 0.14.0 this file was
`docs/SESSION_COMPATIBILITY.md`; older plans and pull requests cite it by that
path and its line numbers. The public
[session and request compatibility reference](../REFERENCE/session-and-request-compatibility.md)
documents current support and what a user does with an older file, including
the [retired inputs](../REFERENCE/session-and-request-compatibility.md#retired-inputs).
Release notes record when a format changed. Each heading below starts with the
release that first shipped its entry, or with `Unreleased:` while the entry is
only on `dev`. A Session version changes when a breaking change enters `main`;
a later breaking change made before that version reaches `main` goes into its
entry ([Persisted contract co-change](./WEB_CHANGE_POLICY.md#persisted-contract-co-change)).

## Unreleased: a Session is a `.gbdraw` container file

A Session is one ZIP file, `<title>.gbdraw`, instead of one JSON file (PD-OI-127).
Its first entry `gbdraw-session.json` is the Session version 51 document without
input bytes, and each input is the entry `resources/<id>` (a gzip input stored as
it is, any other input compressed), so no input is embedded as base64 and the
Web app can write and read the inputs one at a time. The Session version stays
51: the document keeps its fields, and only the file around it changes. The
checksum field is not written for a container (the ZIP's CRC-32 and sizes catch
damage); the older JSON Sessions keep reading it.

The JSON form of version 51 (`.gbdraw-session.json[.gz]`) never reached `main`,
so it is retired without a reader: Python refuses it with "A Session version 51
file must be a .gbdraw file, not JSON: {path}", and the Web app with "Load a
Session file that this version of gbdraw supports.". Every writer (Web Save
Session, `save_session_document`, `--session_output`, `--save_session`) writes
only the container and refuses a name that ends with `.json` or `.json.gz`
before writing anything. Every reader decides by content, not by name: it opens
a container of any name and the JSON Sessions 27 to 33, 39 to 42, 44, and 46 to
50, plain or gzip. A damaged container and a ZIP without `gbdraw-session.json`
are refused with their own messages. gbdraw 0.14.5 and older read the new file
as invalid JSON. An older Session converts with
`gbdraw render --session old.json.gz --session_output new.gbdraw`. The inputs of
a container have no size limit on the command line and in Python; the Web app
limits only `gbdraw-session.json` to 512 MiB, and the limits of the JSON
Sessions (200 MB, 512 MB expanded) stay. The public description is
[The Session file](../REFERENCE/session-and-request-compatibility.md#the-session-file).

## Unreleased: Session 51 stores drawings; request schema 10 writes record windows

Session version 51 replaces the two diagram-mode slices of Sessions 46 to 50
with a list of drawings, so a project can hold several drawings of one mode
(PD-OI-090). Session 50 reached `main` with 0.14.5 as the mode-scoped Session
below, so the list takes 51 (CLAUDE.md, "Persisted-format compatibility"). A
Session 50 is read like 46 to 49: each drawing keeps the Linear wrap length,
italic gene names and the identity colour scale of its mode. A file with
version 50 and `drawings` was never released and is refused (Python: "Session
version 50 requires a canonical renderRequest object."; Web Load:
`INPUT_INVALID`, reason `FIELDS`).

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
`SessionChecked51` in `gbdraw/session_contract.py`. Python and the Web app read
Session 51, and Save writes it. Reading a Session 46 to 50 moves each mode's
part into its own drawing (`circular`, `linear`), the rules of
`services/session-drawings-upgrade.js` in the Web app and of the 46-50 step of
`upgrade_session_document` in Python; Sessions 27 to 44 are converted by Python
once at Load. A Session 49's Circular `adv.show_definition` becomes the
Circular drawing's `config.adv.show_definition`. Fixtures:
`tests/fixtures/sessions/`.

Session 51 also adds the record source kind `fasta` (PD-OI-068 rev 2): a FASTA
given without GFF3, named by `fastaResourceId` and stored as a resource of kind
`fasta`. Its records have sequence and no features, and the file makes no
topology claim. The Web app writes it for the input type `fasta` (the third
source choice, `config.inputType` / `webFiles.bindings.c_fasta` or a Linear
card's `fasta` slot), and the CLI writes it for `--fasta` without `--gff` and for
a records table row with only a `fasta` path. This adds no request schema:
the kind is a new value of the record source's `kind`. Fixture:
`tests/fixtures/sessions/cli-circular-fasta.v51.gbdraw`.

Session 51 also stores a Circular drawing's input as File cards (circular
multi-file input). Its `webFiles.bindings` holds `circularSources`, a list of
`{uid, gb, gff, fasta}` cards; each of `gb`, `gff` and `fasta` is a binding or
null, and the uids are distinct. The CLI writes one card per `--gbk` value, or
per `--gff`/`--fasta` pair, with the uid of the CLI record key, `record-<n>`.
The Circular slice's `recordsOff`, `form.circular_record_selector` and
`adv.multi_record_positions[].selector` name records by identity: `<uid>` for a
file of one record, `<uid>:<k>` for the `k`-th record of a file of several; a
`#N` that the move of an older Session could not place stays. A
Session 51 does not write `c_gb`, `c_gff`, `c_fasta` or the composite, and
refuses a Circular drawing that holds them. Reading a Session 41 to 50 makes one
card per file (the components of a composite, in order; the first card keeps
`c_gff` and `c_fasta`) with uid `record-<n>` (`source-<n>` when the request
drew the files as one source), and replaces a `#N` selector by the identity of
that record where the record counts are known; otherwise the selector stays
`#N`, and the Web app places it over the records of the cards as they were
loaded, taken together in card order
(`circular_sources_of_mode_scoped` in `gbdraw/session_drawings.py`; the Web
twin is `circularSourcesOfModeScoped`). Session 51 and binding schema 2 keep
their numbers. Fixtures: the `v51` files that hold a Circular drawing in
`tests/fixtures/sessions/`.

Session 51 also stores each similarity ring's comparison genome as a list of files. A Circular drawing's
`webFiles.bindings.c_conservation_fastas` and `c_conservation_sequence_sources` hold one list of bindings per ring, in
the ring's file order (an empty list for a ring without files), and `webFiles.conservationLosatFastaSources` and
`conservationSequenceSources` one list of resource IDs per ring (`conservation-losat-fasta-files-<i>-<n>`). A Session 50
or older wrote one binding or resource ID, or null, per ring; the Web app reads it as a list of that one file or an
empty list (`ringGenomeFileLists` in `gbdraw/web/js/services/drawing/inputs.js`), and Python keeps the stored
shape (`drawing_resource_ids` in `gbdraw/session_resources.py` reads both). The ring of one file keeps its series key,
so a Session 50 ring cache entry stays valid. Fixtures: `tests/fixtures/sessions/synthetic_conservation.*` (a Session
40, read by `tests/web/load-silent-rows.test.mjs`) and the schema-5 request in
`tests/web/session-losat-cache-validation.test.mjs`. Session 51 has not reached `main`, so the version does not change.

Session 51 also adds two keys to each record of a drawing's `runMetadata.trackSlotGeometry` (PD-OI-115):
`recordKey`, the request's record key, and `axis`, where the record's axis lies in the frame of its drawn group
(Circular: `radiusPx`, the feature band `bandInnerPx` to `bandOuterPx`, `startAngleDeg`, `degreesPerBp`; Linear: one
row each with `xStartPx`, `xEndPx`, `yPx`, `bandTopPx`, `bandBottomPx`, `displayStart`, `displayEnd`) and which
source base each drawn base shows (`source`: `anchorBase`, `recordBase`, `step`, `sourceLength`, `drawnLength`, and
`isCircular`, whether a selected span may cross the record's origin: the record is circular and has no region
applied, the annotation resolver's rule).
Python writes them (`gbdraw/render/track_slot_metadata.py`); the Web app's Select region reads them
(`services/region-axis.js`). The Session shape keeps each geometry record opaque (`Record<string, unknown>`), so
the version does not change. A saved Result without them loads as before, and Select region stays off until the next
Generate; no reader adds them. Fixture: `tests/fixtures/sessions/HmmtDNA-region-axis.v51-web.*`.

A drawing's `runMetadata` is declared once, in `SessionRunMetadata` (`gbdraw/session_contract.py`): the geometry's
header and the rows of `annotationWarnings`, `featureIdentityNotices` and `comparisonWarnings`. One function builds it
(`gbdraw/api/run_metadata.py`), and every writer stores it as built: a Web Save, and the gbdraw command's run and
re-save of a Session. Each writes the three notice lists, empty when there is no notice, and `trackSlotGeometry` when
a Result has one. This adds no field and keeps the version: a Session that omits a list (an earlier Web Save writes
only `annotationWarnings` when the others are empty) loads it as an empty list, as before.

Measured against `main`'s build (0.14.5, 033a3aa89, Session 50): `load_session_document` on each of the four
`tests/fixtures/sessions/drawings-*.v51.*` fixtures raises `SessionVersionError`: "Session version 51 is newer
than this gbdraw supports (50).". The wrap length, italic gene names and identity colour scale reached `main` in
Session 50 (entry below); Session 51 only moves them from the slices into each drawing.

The full `config` that the CLI and Python store in a Session 51 request no longer has `canvas.linear.canvas_padding`,
a Linear canvas leaf that no code read (the canvas padding is `canvas.padding`). Readers still load a stored config
that has it (`GbdrawConfig.from_dict` reads named keys only) and draw the same figure; an override of that path is
checked as a number and dropped (`RETIRED_CONFIG_LEAVES` in `gbdraw/config/modify.py`). The generated contracts, the
Session version and the request schema do not change for it. 0.14.5 requires the key in a stored config, but it refuses a
Session 51 by its version before it reads the config. Fixtures that store the key stay as they are, for example
`tests/fixtures/sessions/cli-linear-wrapped.v50.gbdraw-session.json.gz`.

Session 51 also adds the config leaves `canvas.figure_width` and `canvas.figure_height` (`"auto"`, the default, or
a length with a unit, issue #1427) to a drawing's stored full config; a CLI or Python writer sets them as config
overrides. The request schema keeps config opaque, and Session 51 has not reached `main`, so the version does not
change. A Session without the leaves reads as auto, and no reader converts it. `main`'s build (0.14.5) ignores the
leaves in a stored full config and refuses them as a `configOverrides` entry ("Unknown config override path").
The Web's Figure size field (Layout, both modes) writes each drawing's `config.form.figure_width` and
`config.form.figure_height`, `copy` rows of the mode-scoped registry: `null` for Auto, else the amount and unit as
typed (`"180mm"`). A set side is the request leaf; Auto writes none, so an Auto drawing's request is unchanged, and
its saved form holds the two keys as `null`. A Session without the form keys reads as Auto. Load reads a stored leaf
the form cannot hold as stored (`"180 MM"`) with Python's configuration check and drops one Python refuses with a
Load repair row. A sized figure's composition metadata carries `requestedCanvas` (`widthPx`, `heightPx`), so the Web
re-plan after a Legend edit keeps the requested canvas; an Auto figure's metadata has no such key.

Session 51 also adds `windowLosatFilenames` to an edge of a Linear drawing's `config.linearComparisonPlan` (OV-2742):
the raw LOSAT file names typed on the window pairs of that pair of records, keyed by the unit edge key
(`<query unit key>-><subject unit key>`, a unit key `<record key>~w<n>` for a window of a record drawn in several).
The Web writes the key only when a window pair has a typed name (an edge without one holds no key), Remove window
moves each name with its window or drops it, and the name edit drops the entries of window pairs that the shown
Result no longer lists for that record pair. Python admits it from Session 51
(`WINDOW_LOSAT_FILENAMES_SESSION_VERSION`) as a map of non-empty keys to non-empty names
(`_validate_comparison_plan_edge` in `gbdraw/session_io.py`) and refuses it in a Session 46 to 50; the Web loads such
a file without it, with a Load row. No renderer reads it. A Session without the key reads as no typed names, so every window pair takes its default name.
Session 51 has not reached `main`, so the version does not change. Positive fixture:
`tests/fixtures/sessions/web-linear-window-pair-names.v51.gbdraw`.

A Linear request that states no `diagramOptions.output`, or a null `plotTitlePosition`, draws with Python's
`LinearOutputOptions()` sides: legend `right`, title `bottom` (OV-2491, Owner decision R17-38). The generated
`ABSENT_LINEAR_OUTPUT` (`schema.generated.js`, from `describe_absent_linear_output()` in
`gbdraw/session_request_codec.py`) carries them to the Web's Load, which reads every drawing that holds a request
and no layout of its own with its request's sides, opened or not (OV-2493). A new Web Linear drawing keeps its legend
and title at the bottom and states them in what it saves. The format does not change: Python and the CLI drew such a
request this way before; only the Web's reading changes. Positive fixture:
`tests/fixtures/sessions/python-linear-absent-output.v51.gbdraw` (Python `save_session_document`).

Session 51 also adds an optional `gradient` to a region annotation's style in the canonical request
(`diagramOptions.annotations.sets[].defaultStyle` and `.annotations[].style`; PR-B, the Web's Sequence search, receipt
`web.sequence-search` rev 1): `startColor` and `endColor`, the colours at the start and at the end of what the region
stands for (a search hit: its query); `queryStart` and `queryEnd`, where the region's two ends lie on that scale as
fractions from 0 to 1 in query order (a hit over a part of its query), written only when they are not 0 and 1; and
`direction`, `forward` or `reverse`, which way that order runs on the record's source strand. A `band` or `highlight` mark draws it along the region in place of its fill (`gbdraw/render/drawers/annotation_marks.py`,
`gradient_steps`), and the Legend draws one bar per caption and pair of colours, from "start" to "end". The writer
writes the key only when a style has a gradient, so every request without one encodes byte for byte as before; the
reader reads an absent or null key as no gradient. This adds no request schema: the key is a new optional field
of the request description (`CanonicalAnnotationGradient`). Session 51 has not reached `main`, so the version does
not change; `main`'s build (0.14.5) refuses every Session 51 file. Positive fixture:
`tests/fixtures/sessions/python-search-hits-gradient.v51.gbdraw` (Python `save_session_document`).

Session 51 also adds an optional `keepHighlightOverlay` to an annotation set (`diagramOptions.annotations.sets[]`,
and the Web draft `config.annotationSets[]`; the same receipt, its custom-stack line): `true` on a set that
Sequence search wrote (`gbdraw/annotations/search_hits.py`, the only writer). Its highlights are then drawn in front of
the feature track, covering its features, also in an explicit track stack (Use custom stack) that has no enabled slot naming the
set: `prepare_annotation_track_slots` adds the set's automatic highlight overlay after the stack's slots, in front of
the stack's feature track. A slot that names the set takes its place; other sets keep following the stack. The writer
writes the key only when it is true, so every other set encodes as before; the reader reads an absent key as false.
This adds no request schema (a new optional field of `CanonicalAnnotationSet`). The same positive fixture holds it.

Session 51 also adds an optional `strand` (`"+"` or `"-"`) to an annotation row (`diagramOptions.annotations.sets[].annotations[]`,
and the Web draft `config.annotationSets[].annotations[]`; the same receipt, its strand line): the strand on the
record's source that a Sequence search hit lies on (`gbdraw/annotations/search_hits.py`, the only writer). Where the
row's highlight covers a feature track that separates the strands (`cover_anchor`), it covers that strand's lane as
the record is drawn (a reverse-complemented record draws the source's plus strand as its minus strand); a row without
it covers the whole feature track, as before. The writer writes the key only when it is set, so every other row
encodes as before; the reader reads an absent key as no strand. This adds no request schema (a new optional field of
`CanonicalRegionAnnotation`). The same positive fixture holds it.

Request schema 10 writes a record's spans as `windows`, a list of `{start, end, reverseComplement, column}` in the
record's source coordinates, in place of `region` (`{selector, start, end, reverseComplement}` or null). A record
drawn as one span has one window with a null `column` (the presentation's `gridColumn` places it), and a record
with no crop has none. A crop names its record in the record's own `selector`, not in the region. Schemas 1 to 9
write `region`; the decoder reads it as one window and a selector inside it as the record's selector, and a stored
`region` is read the same way wherever a payload is read undecoded (`stored_record_windows` in
`gbdraw/session_request_codec.py`). A stored `region` has exactly the keys of `STORED_REGION_FIELDS` in that module, which the
Web reads from `schema.generated.js`; Python refuses any other key, and so does Web Load (OV-3057). Session 51 has not reached `main`, so the version does not change and its
requests are written at schema 10: the fixtures `tests/fixtures/sessions/*.v51*` and the documentation and Gallery
Sessions of Session 51 carry it, and Sessions 50 and older keep their schema 9 requests. `main`'s build (0.14.5)
refuses every Session 51 file by its version.

Request schema 10 also writes `diagramOptions.conservationFastaFiles` as one list of file references per similarity
ring: the files of the ring's comparison genome, in the listed order (a genome split over several files, such as a
chromosome and its plasmids), and an empty list for a ring without one. Schemas 1 to 9 wrote one file reference or
null per ring; the decoder (`_decode_option_value` in `gbdraw/session_request_codec.py`) and the Web promotion
(`promoteRingFilesToLists` in `gbdraw/web/js/services/session-request.js`) read it as a list of that one file or an
empty list (`tests/fixtures/ring-genome-files-vectors.json`, read by both). The file references of ring `i` are
allocated as `conservation-fasta-files-<i>-<n>`. This change entered schema 10 before schema 10 reached `main`, so
the schema does not change again.

Request schema 10 also lets a Circular request draw several records one after another on one circle (issue #200):
`renderRequest.layout` may carry `recordLayout` (`"grid"` or `"one_circle"`), `oneCircleGapDegrees`,
`oneCircleMinLength` and `oneCircleRecordBand`. The writer writes all four only for `"one_circle"`, so a grid layout
encodes as before; the reader reads an absent key as its default (`grid`, 1.0 degree, 1000 bp, band on) and refuses the
keys in a request of schema 9 or older (`ONE_CIRCLE_LAYOUT_OPTIONS` and `RECORD_LAYOUT_SCHEMA` in
`gbdraw/session_request_codec.py`). The Web app does not read them yet. This adds no request schema (schema 10 has not
reached `main`). `tests/test_session_request_codec.py` holds the positive and refusal cases.

A Circular drawing of Session 51 keeps the window of its record that its Start/End fields cannot hold, one window
across the origin of a circular record (start after end, PD-OI-122), in `config.form.circular_windows`: a list of
`{start, end, reverseComplement, column}` as the request writes it, empty for a record drawn whole or in one ordered
window (`circular_region_start`/`circular_region_end` hold that one). It is a mode-scoped row of the Circular slice
(`gbdraw/web_support/mode_scoped_settings.py`). Load fills it from the request (`projectCanonicalSessionRequest`), so
a Session that the Python API saved with `RecordWindow(16000, 100)` loads with the window, and Save and Generate write
it back as the request's window. A stored config without the key reads as an empty list; a value that is not a list
reads as empty, with its Load repair row. The version does not change (Session 51 has not reached `main`).

A Circular grid or batch drawing of Session 51 keeps the crop that a command line or Python API Session gave each of
its records (PD-OI-044 scenario revision 5) in `config.form.circular_record_ranges`: a list of
`{source, record, windows}`, where `source` is the uid of the record's File card, `record` names the record in that
card's file as the request did (its record ID, `#k`, or empty for the record of a file of one record), and `windows`
are the request's windows of that record. Load fills it from the request's grid or batch records; the next Generate
and Save write each record's windows from it while the record is drawn in a grid or batch. A row that names no card,
no record or no window with both ends is dropped with its Load repair row; a stored config without the key reads as
an empty list. `config.form.circular_range_fields_record` names the record (its identity) whose crop the one-record
crop fields were filled from, or is empty, so a change made there stays through All records and draws when that
record is chosen alone again; a stored config without the key reads as empty. The version does not change
(Session 51 has not reached `main`).

A Circular drawing of Session 51 stores the Ring colors setting of its comparison rings (PD-OI-110 scenario
revision 3) in `config.circularConservation`: `colorScheme` (`"same"`, the default, `"ten"`, or `"per_five"`) and
`reassignOnReorder` (`false` by default). Each `series[]` entry's `color` is the ring's own (hand-set) colour, a
lowercase `#rrggbb`, or `null` when the ring takes the scheme's colour at its `schemeIndex` (an integer from 0; under
`same` the ring draws the configured colour and the request sends `''`). Every colour that a Session up to 50 stores
reads as hand-set, since its Web writer filled every ring, so it opens with the colours it was saved with; a missing
`colorScheme` reads as `same`, a missing `schemeIndex` as the ring's place in the list, and an unset or `''` colour
as `null`. A Session that the command line or the Python API saved, which has no Web draft, opens with its
`objects.conservation.color_scheme` as `colorScheme` and its rings without a colour as `null`, each with its drawn
position (the request's slot stack order) as `schemeIndex`; the config path is a Web-managed path, so it is no
preserved setting and the Web never writes it back. A Session 50 that a gbdraw 0.14.4 or 0.14.5 Web Save wrote after it
opened such a Session kept the scheme as a preserved override of its Circular draft
(`unmanagedConfigOverrides["objects.conservation.color_scheme"]`); Load moves it into `colorScheme` with no dialog row.
That Web had filled every ring with `#4e79a7` (OV-3054): while the draft's stored request is a Circular request that
names no ring colour (`diagramOptions.conservationColors`) and whose slot stack the projection can read, a ring colour
`#4e79a7` reads as `null` with the drawn position of its request ring (the request's slot stack order, as for a
Session without a Web draft) as `schemeIndex`. A series entry finds its request ring by the key of the ring's file
among the request's ring files (`ringKeysOfFiles` over `conservationBlastFiles` and their resources), else by the
file's name (the entries a CLI Session's projection wrote hold an empty `sourceKey`), never by its place in the series,
which the user may have changed before the Save; a ring whose colour stays its own takes its request ring's drawn
position as `schemeIndex` too, and an entry no request ring names keeps its colour and its index unless a ring holds
that index (then the first free one), so no two rings share an index. A request that names
ring colours, a request whose slot stack cannot be read (the projection refuses it), and a draft with no Circular
request (a Session 46-50 whose one request is Linear, a flat draft whose request is Linear, or none), keep the draft's
ring colours as its rings' own. Sessions up to 50 (gbdraw 0.14.5) kept the
copy numbers of a ring's `sourceKey` (the order of copies of one file, same name, size and date) in the file order
before the Generate, so after a copy moved the keys of two copies are swapped against their places (OV-3238). Load
needs no repair: a key two copies share is ambiguous, so each series entry finds its ring by its stored place
(`sourceIndex`, `conservationSourcesOfSeries`). The version does not change (Session 51 has not reached `main`).

The mode-scoped registry names the first Session version of each row (`since`, 46 unless a later Session added the
row). `config.form.figure_width`, `config.form.figure_height`, `config.form.circular_windows`,
`config.form.circular_record_ranges` and `config.form.circular_range_fields_record` are `since` 51: no
Session 46-50 writer wrote them, so a Session 46-50 slice that holds one is refused by Python
(`slice_containers(version)` in `_validate_mode_slices`, and the band-46 admission shape) and dropped by Web Load with
a `setting-replaced` row. Python's step of a Session 27-44 through the Session 50 layout writes none of them
(`mode_slice_of_version`), so the drawing takes their defaults. The rows Sessions 48-50 added carry no `since` yet.

## Unreleased: LOSATN and TLOSATX raw-cache args carry the genome pair's search space

No Session or request version changes. `LOSAT_SEARCH_FLAGS` (generated from `LOSAT_OPTION_FLAGS` in
`gbdraw/comparisons/losat_runtime.py`) gains `--searchsp`, spelled `-searchsp` for LOSAT. The CLI
(`gbdraw.comparisons.linear_losat`) and the Web (`buildLosatArgs` in `gbdraw/web/js/app/run-analysis.js`, from
the Worker helper `planGenomeScales`, which runs `linear_genome_scales` in `gbdraw/api/request_render.py` on the
drawing's draft request and answers each scale as a decimal string, with the warnings of genome pairs searched
per file because their search space is past what LOSAT reads) write it into a LOSATN or TLOSATX raw entry's
`args` when a drawn record of either genome is not the whole of its genome; the value is `genome_search_spaces` in
`gbdraw/comparisons/losat_jobs.py`. The args
are part of the raw cache key, so a raw entry an earlier writer stored for such a pair no longer matches and that
pair is searched again. A whole one-record pair writes no
`--searchsp` (when every drawn record of both genomes is the whole of a one-record file), so its key and its
stored entries are unchanged. Readers accept the new word as any other search
argument; no migration is needed. Vectors: `searchSpaceCases` in `tests/fixtures/losat_job_plan_cases.json`, read
by `tests/test_losat_job_plan.py` and `tests/web/linear-sources.test.mjs`.

LOSATP raw entries follow the same rule with `--dbsize` (spelled `-dbsize` for LOSAT; PD-OI-018 revision 6): the
CLI and the Web (the same `planGenomeScales` reply) write it into a LOSATP raw entry's `args` when a drawn record of
either genome is not the whole of its genome; the value is
`genome_database_sizes` in `gbdraw/comparisons/losat_jobs.py`, the subject genome's protein letters
(`source_protein_letters` in `gbdraw/analysis/protein_colinearity.py`), which the request planner computes
(`_losatp_database_sizes` in `gbdraw/api/request_render.py`). A LOSATP raw entry stored for such a pair before no
longer matches and is searched again; a whole one-record pair's key is unchanged.

Circular similarity ring raw entries (flow `circular-conservation`) follow the same rule with `--searchsp`
(PD-OI-018 revision 6): the CLI (`resolve_circular_conservation_losat` in `gbdraw/comparisons/circular_losat.py`) and
the Web (the ring flow of `gbdraw/web/js/app/run-analysis.js`, from the Worker helper `planRingScales`, which runs
`plan_ring_scales_json` in `gbdraw/web_support/comparison_sequences.py` on the reference files and each ring's
files) write the search space of the reference genome and the ring's genome (`ring_search_spaces`) into the entry's
`args`, unless both genomes are one file of one record. A ring entry stored before for such a ring no longer matches
and is searched again; a one-record pair's key is unchanged. Vectors: `tests/web/circular-ring-subject-oracle.py`,
read by `tests/web/run-analysis-simple-path.test.mjs`, and `tests/test_ring_genome_evalue.py`.

## Unreleased: request schema 10 stores feature density tracks

Request schema 10 (not on `main` yet, so its number stays) adds the optional
`diagramOptions.featureDensityTracks`: a list of tracks, each
`{id, query: {featureType, qualifier, value, minLength, maxLength}, metric,
style, window, step, color, caption}`, every key written (`null` where unset).
`metric` is `count` or `coverage`, `style` `heatmap` or `histogram`. The writer
is `_encode_feature_density_tracks` in `gbdraw/session_request_codec.py`, the
typed owner `FeatureDensityTrack` (`gbdraw/analysis/feature_density.py`) and
`FeatureQuery` (`gbdraw/features/query.py`). An absent key or `null` reads as no
tracks, as every earlier request does. A track slot names a track by id
(`feature_density@track=<id>`) in the existing slot lists. The version of
Session 51 does not change (it has not reached `main`); each drawing's draft
gains `config.adv.feature_density_tracks`, a mode-scoped row of Session 51
(`gbdraw/web_support/mode_scoped_settings.py`, `since` 51): the tracks as the
request writes them, empty for a drawing without tracks. Load fills it from the
request (`projectCanonicalSessionRequest`), Generate and Save write it back as
the request's `featureDensityTracks` (the key only when the list is not empty),
so a CLI Session's tracks survive a Web Generate and Save. A stored config
without the key reads as no tracks; an entry that leaves the request's track
shape is dropped with its Load repair row (`services/feature-density-tracks.js`).
Fixture: `tests/fixtures/sessions/HmmtDNA-feature-density.v51-cli.gbdraw`.

## Unreleased: depth and GC Percent axis label styles

No Session or request version changes. The config leaves `objects.depth.axis_labels` and
`objects.gc_content.axis_labels` (`"range"` or `"numbered"`) say how a depth or GC Percent axis is labelled: one
`[min–max]` range label, the default since D-16 (Owner 2026-11-08), or a number on each large tick, as every writer
before drew it. Each Session 51 drawing's draft holds them as `config.adv.depth_axis_labels` and
`config.adv.gc_content_axis_labels` (Web: **Numbered Ticks**), and every writer of `configOverrides` (the Web, and
Python for a request with overrides) now states both leaves (D-24).

A reader takes an absent style for an older writer's: a stored full config (`migrate_legacy_full_config`), stored
`configOverrides` (`_decode_config_overrides`), and a stored draft (`readScalarAxisLabels` in
`gbdraw/web/js/services/config-field-readers.js`) without it read as `"numbered"`, the value
`LEGACY_SCALAR_AXIS_LABELS` of `gbdraw/session_legacy_values.py`, which `tools/generate_web_contracts.py` carries to
the Web. A request without config and without overrides draws the default. Evidence: tag v0.14.5 (no writer states
the style); positive fixtures `tests/fixtures/sessions/REGEX-scalar-axes.v50-cli0145.gbdraw-session.json.gz` (a
full config the v0.14.5 CLI wrote) and `HmmtDNA-wrap-5kb.v50-web0145.gbdraw-session.json.gz` (overrides and a
draft a v0.14.5 Web Save wrote), read by `tests/test_scalar_axis_tick_labels.py`; the Web Load records of the
legacy fixtures (`legacy-request-load.expected.json`) hold the numbered drafts.

## Unreleased: unset collinearity parameters are written as their defaults

No Session or request version changes. A `generatedProteinComparison` holds `settings.collinearityParams`, which may be
`null` or a `lossless` object; every Python reader reads `null` as `LosslessCollinearityParameters()`. main's codec
(`_encode_pipeline_value` in `gbdraw/session_request_codec.py`) wrote `null` for options that leave the parameters
unset: Sessions the Python API saves, and the drawing Find in other files derives from a drawing without a LOSATP
comparison. The Web and the CLI always wrote the object. The Web comparison classifier
(`services/imported-comparison-intent.js`) could not hold `null` in the Comparison controls, so such a comparison
opened as "Saved comparison preserved as read-only" (OV-3382).

The codec now writes the defaults as the object, and the request description gives that object as the field's
default (`describe_canonical_request()`, `request-shape.generated.js`, the JSON Schema), so a repair of an unreadable
value fills it too. The Web reads a stored `null` as that described default (`readStoredCollinearityParams` in
`services/losat-normalization.js`). A request without saved Web controls (a CLI or Python API Session, a drawing Find
in other files creates) that lists the pairs of a LOSATP Pairwise pipeline now loads them as the Selected pairs of the
Comparison controls; before, it loaded the adjacent pairs (`projectCanonicalSessionRequest` in
`services/session-request.js`). The `null` reader stays: main writes it at request schema 9. Fixture:
`tests/fixtures/sessions/python-api-linear-pairs.v50.gbdraw-session.json.gz` (main's Python API, `null`), run by
`tests/web/session-cli-compatibility.test.mjs`.

## Unreleased: the Web reads request text patterns as Python does

No Session or request version changes, and no writer changes. Six text values of a canonical request have a
pattern in `describe_canonical_request()` (`gbdraw/session_request_codec.py`): a record's `recordKey` (a character
that is not a space), the two forms of an annotation color (`#rgb` or `#rrggbb`, and a color name, each with outer
spaces), a region annotation `fill` of spaces only, and the Circular and Linear `dinucleotide`. Python reads a space
as `re` does: U+0009–U+000D, U+001C–U+0020, U+0085, U+00A0 and the other Unicode spaces, but not U+FEFF. The Web
Load (`readStoredRenderRequest`, through `shapeViolations` in `gbdraw/web/js/utils/contract-shape.js`) read the same
pattern text as ECMA-262 does, which counts U+FEFF as a space and U+001C–U+001F and U+0085 as none. Before this
change, a stored request whose `recordKey` was U+FEFF alone was dropped at Load, and one whose `recordKey` was
U+001C alone loaded and was refused at Generate; a `dinucleotide` with U+001C around it was replaced by its
default, and one with U+FEFF around it was kept.

`tools/generate_web_contracts.py` now writes each pattern of `request-shape.generated.js` in the ECMA-262 spelling
that `docs/REFERENCE/schemas/gbdraw-render-request.schema.json` already used (`ecma_types`), so the Web Load keeps
and repairs exactly what Python accepts and refuses. No gbdraw writer creates such a value; a file holds one only
when it was edited by hand. Vectors: `requestTypes` in `tests/fixtures/contract_shape_vectors.json` (U+001C and
U+FEFF as a `recordKey` and around an annotation color), run by `tests/test_contract_shape_vectors.py` and
`tests/web/contract-shape.test.mjs`.

## Unreleased: a Circular track gap below zero, or a number that is not finite, reads as its default at Web Load

No Session or request version changes, and no writer changes. A Circular track slot's `innerGapPx` and `outerGapPx`
now have the bound `minimum: 0` in `describe_canonical_request()` (`gbdraw/session_request_codec.py`): the domain that
`CircularRequestTrackOptions` already checks when Python reads a request (`parse_optional_pixel` in
`gbdraw/tracks/circular.py`, zero allowed), and that the CLI's `--track_slot ...@inner_gap_px=-1` refuses.
`tools/generate_web_contracts.py` carries it to `request-shape.generated.js` and
`docs/REFERENCE/schemas/gbdraw-render-request.schema.json`. Before this change, the Web Load
(`readStoredRenderRequest`) checked only that a gap was a number, and a gap below zero reached the request projection,
which refused the whole Session with an unrecognized error (OV-3237). The Load now replaces such a gap by its default,
`null` (Auto), with a `setting-replaced` row, as it does for a gap that is no number.

The two structural checkers (`shape_violations` in `gbdraw/contract_shapes.py`, `shapeViolations` in
`gbdraw/web/js/utils/contract-shape.js`) also no longer take a number that is not finite at a `number` field. JSON has
none, but Python's json reads `1e400` as inf (and `NaN`, `Infinity`), and `JSON.parse` reads `1e400` as Infinity. An
`integer` field refused it already. For every `number` field of the request description, Python's reader already
refused such a value (its typed check reads a `float` as a finite number, and each option's `__post_init__` checks its
own); the Web Load kept it, so `1e400` in a gap also ended in the unrecognized error. The Load now repairs it like any
other value outside the description, with a row. No writer writes such a number: the Python request encoder runs the
same typed check before it writes, the typed-resource and derived-result writers use `allow_nan=False`, and the Web
writer's `JSON.stringify` writes none; no stored Session, fixture, Gallery or documentation JSON in the repository
holds one. The published JSON Schemas take `1e400` as a number, as JSON Schema reads any number text
(`tests/test_contract_json_schema.py` lists the vector).

Python, the CLI and `gbdraw render --session` still refuse such a request; a region annotation coordinate or lane of
`1e400` is now their decoding error, not an `OverflowError` (OV-3390, `gbdraw/annotations/models.py`). No gbdraw
writer writes a gap below zero; a file holds one only when it was edited by hand. Vectors: `circular-slot-gap-below-zero`
and `circular-slot-gap-of-zero` in `requestTypes`, and the three `1e400` cases in `jsonDifferences`, of
`tests/fixtures/contract_shape_vectors.json`, run by `tests/test_contract_shape_vectors.py`,
`tests/web/contract-shape.test.mjs` and `tests/test_contract_json_schema.py`; Load test:
`tests/web/stored-request-read.test.mjs`.

## Unreleased: the Web reads a Record Order row as Python does

No Session or request version changes, and no writer changes. `RECORD_ROW_PATTERN`
(`gbdraw/layout/record_placement.py`, carried to the Web by `tools/generate_web_contracts.py` as the
`schema.generated.js` export of the same name) is the row of a `layout.multiRecordPositions` token, the text after its
last `@`, as `parse_record_row_position` reads it: spaces around it, an optional `+`, at most 4300 ASCII digits (the
digits `int()` converts by default, leading zeros counted; the pattern carries the limit), and a value above zero. The
Web read that row with `Number()`, so `#1@1.0`, `#1@0x1` and `#1@1e0` drew row 1 on the Web while
`gbdraw render --session` refused them (OV-3313). The Web reads every row with this pattern (`recordOrderRows`,
`services/record-draw-selection.js`), from the token as stored (JavaScript's `trim()` and Python's `strip()` take other
characters). The renderer reads only the order of a list's rows, and reads rows above 2^53 exactly; the Web reads the
rows of one list as numbers in that order: each row's value when every row is a safe integer, else each row's rank.

Load refuses a request for its Record Order where the CLI refuses it, and nowhere else. A Circular grid request, or a
Linear request of schema 6 or later, whose token has a row Python refuses (these three, and `0`, `-1`, text, or more
than 4300 digits) leaves no request to draw, as a layout its owner refuses does, and the Load repair dialog names the
token's place (`renderRequest.layout.multiRecordPositions[<i>]`). A one-circle Circular request keeps its tokens: the
renderer draws one circle without reading them. A schema 2 to 5 Linear request turns each card's `#<card>@<row>` into
its `gridRow`; Python and the Web read the row text (the whole token when it holds no `@`) with the same pattern
(`record_row_number`, `recordRowNumbers`). Python read it as ASCII digits only, so a hand-edited `+2` or ` 2 ` gave no
`gridRow` and now gives row 2; the Web read it with `Number()`, so `2.0`, `0x2` and `2e0` gave row 2 and now give none.
When some cards take a row and some none, the CLI refuses the request ("If one Linear record has grid_row, every
record must have it."), and Load now leaves it undrawn too, naming the first card's token without a row (or the token
list when a card has no token); when no card takes a row, both draw the cards in their order. Schemas 6 and 7 keep
their Linear tokens, which the renderer reads by their selectors: the Web promoted them to each card's `gridRow` by the
token's place, so `#2@1`, `#1@2` drew record 1 first on the Web and record 2 first in the CLI, and a row the CLI
refuses gave no row instead of refusing the request. No stored file in the repository holds a schema 6 or 7 Linear
token.

When Load leaves out a Session 51 drawing's request, the drawing stays: its saved diagram, Results, feature list and
run metadata are left out, and it keeps its settings and the input files they name (`withStoredResultSet`). A drawing
saved by the Web or the CLI names its input files in its bindings and keeps them; a Gallery publication names them only
in its request and keeps none. The dialog said "The file loads as settings only; the saved Results, the GenBank, GFF
and FASTA inputs, the feature list and the caches are left out", which is true of a Session 27 to 50 only; a Session 51
drawing's row now says what it keeps (`drawingStays` on the `request-dropped` row, OV-3472). On `dev`, a row of `0` or
`-1` read as row 1 with a row naming `config.adv.multi_record_positions[<i>].row` when the drawing held no Web draft,
and was not read when it held one. `#1@+1`, `#1@01` and `#1@ 1 ` read as row 1, as before. No gbdraw writer writes
another row; a file holds one only when it was edited by hand, and every stored token in the repository's fixtures,
Gallery and documentation reads. Vector: `tests/fixtures/record_order_row_vectors.json`, run by
`tests/test_web_contracts_generated.py`, `tests/test_canonical_request_shape.py`,
`tests/web/record-draw-selection.test.mjs` and `tests/web/session-request.test.mjs`; Load test:
`tests/web/stored-request-read.test.mjs`.

## Unreleased: a `#+N` record selector reads as record N

No Session or request version changes, and no writer changes. `RECORD_INDEX_SELECTOR_PATTERN`
(`gbdraw/io/record_select.py`, carried to the Web by `tools/generate_web_contracts.py` as the `schema.generated.js`
export of the same name) now allows one `+` after `#` and its spaces, so `#+2` and `# +02` name record 2 wherever
gbdraw reads a `#N` record selector, as Linear Record Order (`--multi_record_position`) always read them (OV-2988).
Before, the move of a Session 41 to 50 (`gbdraw/session_drawings.py`, Web `services/session-drawings-upgrade.js`)
kept such a selector as it was, a key no record has; it now places it at the record's identity, as it places `#2`.
`#+ 2`, `#-1`, `#+0` and `#++2` stay as they were. No gbdraw writer writes a `+`; a file holds one only when it was
edited by hand. Vector: the OV-2988 case of `tests/fixtures/sessions/drawings-vectors.json`, run by
`tests/test_session_drawings_vectors.py` and `tests/web/session-drawings-upgrade.test.mjs`.

A hand-edited Circular Record select `#+2` of a Session 41 to 50 matched, on `dev`, a record whose ID is `#+2` (the
Record select compares a value with the record's place first, then with its ID); it now names record 2, as `#2` did.

The Web Load also reads a stored record selector of kind `recordId` whose value reads as `#N` (`#2`, `#+2`) by its
kind (OV-3311). A Linear card, the Circular Record select and a Circular record range hold their record as text,
which reads `#N` as a record number, so on `dev` a stored record ID `#2` drew record 2 in the Web while Python drew
the record with that ID (with the pattern above, `#+2` would too). Load now lists the record's file once through the resources it adopted (`readStoredRecordIdPlaces`,
`services/record-discovery.js`) and gives the text the record's place `#k`; a file it cannot list, or one that holds
no single record with that ID, keeps the text and adds a `record-id-unresolved` row to the Load repair dialog. A
Session with no such selector reads no file for this. No gbdraw writer stores such a selector (the Web and Python
writers name such a record by its place). Test: `tests/web/stored-record-id-places.test.mjs`.

## 0.14.5: Session 50 stores the Linear wrap length, italic gene names, and the identity colour scale

Session version 50 adds the config leaf `canvas.linear.wrap_length` (`"none"`,
`"auto"`, or a positive length in bp) and the Linear slice's
`form.linear_wrap_length`, an `own` row of the mode-scoped registry
(`gbdraw/web_support/mode_scoped_settings.py`) whose replay flag is
`--wrap_length`. Wrapped rows are drawn now: Python, the CLI and a loaded Session
draw a record longer than the length over several rows, so a writer can produce
a wrapped file. Session 49 reached `main` with 0.14.3, so the format needs a new
version. Measured against `main`'s build: it refuses a Session with
`modes.linear.config.form.linear_wrap_length` (the slice's `form` is closed,
`_validate_mode_slice_fields`) and a `configOverrides` entry
`canvas.linear.wrap_length` (`UNKNOWN_CONFIG_PATH`), and it accepts a CLI Session
with a stored `canvas.linear.wrap_length` leaf and silently draws the record on
one row (`from_dict` ignores the unknown leaf). An unwrapped CLI Session holds
the leaf as `"none"`, which an older build ignores and draws the same.

A Session up to 49 has no key and loads with no wrap (`"none"`); no reader
converts it. A Web **Save Session** keeps a length loaded from a Session in the
request's `configOverrides`; the Web has no control for it yet. The Gallery and
test-input Sessions were refreshed to 50 (`tools/refresh_gallery_sessions.py`);
they differ from Session 49 in the version and the `"none"` leaf. Positive
fixtures of Session 49: `tests/fixtures/sessions/*.v49.*`.

Session 50 also adds the config leaf `labels.<mode>.italic_gene_labels` (true or
false, default false) in both modes and the slice key
`adv.italic_gene_labels`, a `copy` row of the mode-scoped registry next to
`label_protection`, whose replay flag is `--italic_gene_labels`. A Web writer
puts `labels.<mode>.italic_gene_labels: true` in the request's
`configOverrides` only when the setting is on, so a Session with it off has the
same request as before. The key reached no `main` Session before 50 and shares
its version with the wrap length; a Session up to 49 has no key and loads with
it off, and no reader converts it.

Session 50 also adds the config leaves `objects.blast_match.identity_scale`
(Linear) and `objects.conservation.identity_scale` (Circular), `"auto"` or
`[MIN, MAX]` percent identities, and the slice key `adv.identity_scale`
(`null` or `[MIN, MAX]`), a `copy` row of the mode-scoped registry next to
`identity`, whose replay flag is `--identity_scale MIN MAX`. A Web writer puts
the mode's leaf in the request's `configOverrides` only when a scale is set, so
a Session without one has the same request as before. `main`'s build refuses
that override (`UNKNOWN_CONFIG_PATH`) and the slice key (the slice's `adv` is
closed); it ignores the leaf in a stored full config and draws the threshold
scale, as a Session without a scale draws. The keys reached no `main` Session
before 50 and share its version; a Session up to 49 has none and loads without
a fixed scale, and no reader converts it.

## 0.14.3: Session 49 stores whether a Circular slice shows the centre definition

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

## 0.14.4: Collinear defaults and value sets

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

## 0.14.5: derived LOSATP token v6 for a CDS across the origin in the Web and for score ties ordered by CDS

No Session, request, or derived envelope version changes. 0.14.4 writes token
`derived-option-conformance-v4`; the token becomes `derived-option-conformance-v6`
for two reasons, and no released writer wrote v5. Earlier derived entries are
cache misses and are recomputed once from the raw cache at the next Generate;
no search runs again and no reader reads them.

- A CDS across the origin in the Web (OV-2368). The Web converter now receives
  each CDS's `origin_arc` and `record_length` and the record lengths, so a Web
  result with a CDS across a circular record's origin has new rows for the same
  inputs. The derived identity names the bindings, settings, and raw keys, not
  the protein map, so an entry that 0.14.4 wrote would be a hit with the old
  rows. The protein map itself is never saved: it lives in the page's per-file
  memory (`readOncePerFile`), and no Session field holds it, so
  `proteinMapSchema` and the Session stay as they are. The comparison tables of
  a saved request keep their rows until that Generate.
- Score ties ordered by CDS (OV-2370, OV-2380). A score tie between protein hits
  (a Pairwise or member hit limit, a best hit, a group's representative, the
  bitscore fit, the order of edges) was broken on the runtime handle, which
  hashes the record instance key: the CLI's record key or the Web page's
  sequence uid. Equal inputs gave other links on the CLI and the Web, and on two
  Web pages. Ties now follow the CDS order (record, then location), so a derived
  entry written before has rows of the old order for the same identity.

## 0.14.4: derived LOSATP token v4 for chained Collinear blocks

No Session, request, or derived envelope version changes. Collinear blocks are
now the longest compatible chains of anchors (OV-2294..2296), and the block
rows no longer carry `collinearity_block_evalue` (OV-2300), so a derived entry
written before would replay the old blocks. The derived LOSATP identity token
becomes `derived-option-conformance-v4`: earlier derived entries are cache
misses and are recomputed once from the raw cache; no search runs again and no
reader reads them. Saved `collinearity-result` resources keep their shape:
`CollinearityBlock` still has `blockEvalue`, `label`, and the `syntenic` kind,
which the block caller no longer writes.

## 0.14.3: LOSAT runtime records without `cli` or `bundled`

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

## 0.14.3: the palette `hearbal_garden` reads as `herbal_garden`

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

## 0.14.3: per-feature rows follow their record to a new key

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

## 0.14.2: one comparison-match attribute in saved Results

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

## 0.14.0: Session 48 stores protein results in the search frame, what each file card loaded, and the feature label style

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

## 0.14.0: Session 47 keys the Legend edits by the row key

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

## 0.14.0: Session 46 keeps the settings and the Result of each diagram mode

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

## 0.14.0: Web Load of 0.13.0 Sessions

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

## 0.14.0: Web Load of Session 31–39 table rows

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

## 0.14.0: Session 46 and Web feature edits by source identity

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

## 0.14.0: request schema 9 and feature identity overrides

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

## 0.14.0: value checks, derived label maps, and definition spacing

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

## 0.14.0: records turned off for drawing

Session version 46 and request schema 9 are unchanged. A slice of `modes` may
hold `config.recordsOff`, the list of record keys the Web app does not draw: a
Linear card uid, or a Circular source selector `#N`. A Session without it, and
`[]`, read as all records on. The Web app writes it only when a record is off; the render request
holds only the records that are on, so the request has no new field. A CLI
`--session` replay draws only the records that are on, and a CLI re-save keeps
`recordsOff` and the settings of the records that are off.

## 0.14.0: default label of a precomputed similarity ring

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

## 0.14.0: a scale interval of 0 or less

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

## 0.14.0: CLI and Python LOSATN / TLOSATX results

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

## 0.14.0: CLI LOSATP results in Sessions

Session version 47 and the request schema are unchanged. A Linear LOSATP run
that writes a Session (`--session_output`, `--save_session`, a Gallery
rebuild) stores the results it computed: the request's
`generatedProteinComparison` marker has mode `none`. Similarity groups stores
the `precomputedProteinComparison` tables and the `orthogroupResult`; Collinear
stores only the `collinearityResult` (no tables); Pairwise stores nothing, so
its first Generate still searches (OV-1068). The first Web Generate after Load
reuses what is stored instead of searching again.

## 0.14.0: one LOSAT cache entry per raw key

Session version 44 is unchanged. The Web app saved a raw key once per row, so a
Session with two Circular rings of one sequence repeated the key and could not
be loaded. It now saves one `losatCache` entry per raw key, named by the first
row, as the CLI does. After a load, the cache list shows one row per raw key
until the next Generate. A version 39 or later Session that repeats a raw key
is rejected with an `INPUT_INVALID` Session diagnostic.

## 0.14.0: comparison rows in the search frame

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

The current writer and the versions current readers accept are in the table of
the public
[session and request compatibility reference](../REFERENCE/session-and-request-compatibility.md).

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
A Session 51 has no composite: the move of a Session 41 to 50 reads it and
makes one File card (`circularSources`) per component.

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
