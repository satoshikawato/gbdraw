# Session compatibility fixtures

`BGC0000708-BGC0000713.v39.gbdraw-session.json.gz` is the unchanged session
JSON from first-parent `main` commit
`17e2c9dee32724219aa7c96a02d183280ffbe438`, compressed with `gzip -n -9`.
Its decompressed SHA-256 is
`9407365a3d5684490f1e72d81826101e661b9c0c362ec0b4e0b164af8e1b50b5`.

The schema-v2 fixture is the older supported compatibility case. Its expected
projection is recorded in the adjacent `.expected.json` file.

`lambda_basic_linear.v31-schema1.gbdraw-session.json.gz` and
`HmmtDNA_basic_circular.v32-schema1.gbdraw-session.json.gz` are the positive
fixtures of request schema 1, copied unchanged and compressed with
`gzip -n -9`. The first is the version 31 Linear Gallery Session from
first-parent `main` commit `10d3a3d28b3c9faa42db01c8bd0bd36b9be8433c`
(`gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json`); its
decompressed SHA-256 is
`15c8f32338453ede2c11e9976a7d972483483b1e2e85820d13398d3304326596`. The second
is the version 32 Circular Gallery Session from first-parent `main` commit
`be13a5c030bbfba9e874cbaf22a74d08015532a9`
(`gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json`); its
decompressed SHA-256 is
`e565f29c358a0c372d967a29a2c58654c87417fa120800a90a05f6e1528b3429`.

`../legacy-render-request-vectors.json` is the shared vector of the Web lift
of request schemas 1 and 2 and Python's decoder
(`tests/web/legacy-render-request.test.mjs`,
`tests/test_legacy_render_request_lift.py`). Its `requests` are every distinct
request schema 1 and 2 `renderRequest` stored in a Gallery Session of
first-parent `main` history (23 requests in 47 Session blobs, from `10d3a3d2`
to the last schema-2 writer before `17e2c9de`), as stored, each with the
commit, path, blob and Session version of every place that holds it; the four
committed (33, 2) fixtures are listed by reference (`fixtureReferences`). Its
`cases` are forms the corpus does not hold, written by hand: a corpus request
with edits, what the lift writes and lists, and whether Python decodes both
alike, the finding of a deliberate difference, or that Python refuses both;
a case may instead expect the lift to refuse a form no writer stored (a
field of a later schema, or an output without `outputPrefix`). Its SHA-256 is
`e0aa2c6e0fa29580418b1e94e8f74ed96aa03e77fb94664e1c5956143cc38c18`.
`legacy-request-load.expected.json` records what the Web Load and Save Session
do with the six full legacy fixtures (`tests/web/legacy-request-sessions.test.mjs`).

`HmmtDNA_basic_circular.issue-469.json.gz` preserves the unmodified JSON from
`0f00436da728402d06d4a8fdce80cff63b552488:gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json`.
Its decompressed SHA-256 is
`69786dd18f7a431441085fad3e7d20ed9b54ed2f0b5aeac9685dbb3a704e6da1`.
`test_run_info_exact_replay.py` applies only the title and unused-resource changes
specified in [Issue #469](https://github.com/satoshikawato/gbdraw/issues/469).

`HmmtDNA_basic_circular.v44-schema7.json.gz` preserves the released version 44,
request schema 7, feature catalog schema 4 Gallery session from
`4d1cf93514d0f75fa7a0ee32c1c4b2f176e03682:gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json`
before its owner refresh to request schema 8. It is compressed with `gzip -n -9`.
Its decompressed SHA-256 is
`d5758ff5fbd22c3a9ecb277f716b88c7766d7aeda81a7d2e45ca68ff9c88985e`.

`BGC0000708-BGC0000713.v30.gbdraw-session.json.gz` preserves the unmodified
version 30 Gallery session from release tag `0.13.0`
(`git show 0.13.0:gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json`),
compressed with `gzip -n -9`. Its decompressed SHA-256 is
`a4793e572600d62767d9bb6be645952e8a515d3dc80ec1ce096e082217b91276`. Its
`cliInvocation.args` use the retired LOSATP flags, so
`test_losatp_option_names.py` replays them through the legacy argv rewrite.

`cli-linear-protein.v30.replay.svg` is the SVG written by replaying
`cli-linear-protein.v30.gbdraw-session.json.gz` (see its `.provenance.json`) with
`gbdraw linear --session ... -o replay -f svg` before the LOSATP option rename. The
rename must not change it. It was redrawn when the Legend took the drawn text width
(Owner V-D03, 2026-10-24; only the Legend moved) and again for the OV-165 labels, ticks, definition and canvas change. SHA-256:
`d10b2c96e88477a1b1ae6b53988591fb143a7c86eced9129dcab132517370e55`.

`conservation-fasta.v39.gbdraw-session.json.gz` is a version 39 Circular CLI
session written by first-parent `main` commit
`cdd31013` (`git archive cdd31013 gbdraw`) for
`gbdraw circular --gbk HmmtDNA.gbk --conservation_blast danio-human.tlosatx.tsv
--conservation_fasta NC_002333.2.fna --conservation_reference subject
--conservation_labels Danio --identity 40 -o argv-v39 -f svg --save_session`,
compressed with `gzip -n -9`. Its decompressed SHA-256 is
`a885bdad7bf8778c0836c686710111f85871ffecbdd7da3e460b971f8701884a`. It is the
positive fixture for the retired `--conservation_fasta` flag (design D18).

`composite-circular-three-files.v44-schema8.gbdraw-session.json.gz` is the
Circular CLI Session for three single-record GenBank files (see its
`.provenance.json`). Circular has one GenBank File, so the Session binds the
three files as one composite `c_gb` whose components keep each file's name and
bytes. `composite-session-resources.playwright.spec.js` loads, saves, regenerates
and replays it. Its decompressed SHA-256 is
`9c371d16a5b74bebffd28d9efdff16e8678d0fdecd1bfe1c07a326f33078571b`.

`feature-edits-crop-rc.v44.gbdraw-session.json.gz`,
`feature-edits-circular-copies.v44.gbdraw-session.json.gz`,
`feature-edits-circular.v33.gbdraw-session.json.gz`, and
`feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz` are Web **Save
Session** downloads, kept unchanged, from first-parent `main` commits
`fe6861f0` (Session 44) and `b05a6bb8` (Session 33). They hold Feature
visibility, Label visibility, and label text edits keyed by rendered feature
ID: a cropped and a reverse-complemented Linear record (Sessions 44 and 33),
two Circular records with the same ID, and a Circular Session 33. The
Sessions 33 have no feature catalog. They are the positive fixtures for the
readers that move those edits onto source identities (Session 46). The steps,
inputs, and hashes are in `feature-edits.provenance.json`.

`forced-label-underlay-strokes.v44.gbdraw-session.json.gz` is a Web **Save
Session** download, kept unchanged, from first-parent `main` commit `fe6861f0`
(Session 44). Its Result shows a feature stroke, a Legend row stroke, and a
Legend row color without `data-gbdraw-base-*` records. It is the positive
fixture for the Load reader that records them from the strokes of the parts no
edit reached, `originalColors`, and catalog fills (`createSavedResultPlan`),
which runs for Sessions older than 46. The steps, inputs, and hashes are in
`forced-label-underlay-strokes.provenance.json`.

`forced-label-underlay-legend-rows.v44.gbdraw-session.json.gz` is a Web **Save
Session** download, kept unchanged, from first-parent `main` commit `fe6861f0`
(Session 44). Its Result shows a renamed, an added and no deleted Legend row,
without Legend structure records: the deleted row is not in its bytes. It is
the positive fixture for loading such a Result with the same Legend, and for
the Restore that asks Python to draw a row the bytes lack (U3a O-2). The steps,
inputs, and hashes are in `forced-label-underlay-legend-rows.provenance.json`.

`whitelist-tab-keyword.v39.gbdraw-session.json.gz` is a Web **Save Session**
download, kept unchanged, from first-parent `main` commit `17e2c9de`
(Session 39). Its Label whitelist rule was typed with a tab in the keyword, and
the Session 39 writer stored the row as four cells. It is the positive fixture
for the reader that reads Session 31–39 Default colors, Label whitelist, and
Qualifier priority rows as the current writer writes them. The steps, inputs,
and hashes are in `whitelist-tab-keyword.provenance.json`.

`selected-feature-annotations.v44.gbdraw-session.json.gz` is a Web **Save
Session** download, kept unchanged, from first-parent `main` commit `fe6861f0`
(Session 44). Its three Region Annotations were made from selected features,
so each names its feature by `hash=`: a feature whose hash names one source
feature, one of two CDS at the same coordinates, and a feature of a
reverse-complemented record. It is the positive fixture for the reader that
moves such targets onto source identities (Session 46). The steps, inputs,
and hashes are in `selected-feature-annotations.provenance.json`.

`hash-color-rule-linear-rc.v30.gbdraw-session.json.gz` is a Web **Save
Session** download from release tag `0.13.0` (Session 30), compressed with
`gzip -n -9`. Its one "This feature only" color rule names a CDS of a
reverse-complemented Linear record by its rendered ID, the hash of the drawn
feature with `_record_2`. It is the positive fixture for the reader that names
such `hash` rules by source hash. The steps, inputs, and hashes are in
`hash-color-rule-linear-rc.provenance.json`.

`mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz` and
`mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz` are Web
**Save Session** downloads from release tag `0.13.0` (Session 30), compressed
with `gzip -n -9`. Each was saved after Generate and a switch to the other
mode, so its Result is of one mode while `ui.mode` names the other; the
Result's mode is in `cliInvocation.mode` (0.13.0 writes no
`ui.generatedMode`). Release tag `0.12.0` (Session 29) writes the same
field (`gbdraw/web/js/services/config.js:1935` and
`gbdraw/web/js/app/run-info.js:151` at that tag). Both releases write it only
for a Generate whose file bindings the CLI can replay, so a Session saved
after, for example, a Linear LOSAT run has none. They are the positive
fixtures for the reader of a Session 27–30 Result's mode (OV-677). The steps,
inputs, and hashes are in `mode-switch.provenance.json`.

`lambda_basic_linear.v44-schema8.gbdraw-session.json.gz` and
`HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz` preserve the released
version 44, request schema 8, feature catalog schema 4 Gallery Sessions from
first-parent `main` commit `fe6861f0`
(`git show fe6861f0:gbdraw/web/gallery/sessions/<name>.gbdraw-session.json`),
compressed with `gzip -n -9`, before the Gallery refreshes after it. Their
decompressed SHA-256 values are
`46d72e44f6f54c0fbf6c9c93c806a2f11570e1d024fa3f7552312042637abbad` (lambda) and
`e532e83bf78d6ad63cfb228336b7ff3dcb3b82076a399da787e771b69abce649` (HmmtDNA).
`gallery-session-publication.test.mjs` promotes them to the current writer.

`two-mode-project.v44.gbdraw-session.json.gz` is a Web **Save Session**
download, kept unchanged, from first-parent `main` commit `fe6861f0`
(Session 44). Both modes were used: Circular holds `TESTA.gb` with a Depth TSV
and Circular title, font, legend, GC, and Depth settings; Linear holds `TESTA.gb`
(same Depth TSV) and `TESTB.gb` (reverse complemented) with Legend edits, a color
rule, a feature placement, and a label edit. The Linear Result is the committed
Result while `ui.mode` is `circular`, both modes have a staged record-display
row, and the inactive Linear profile holds an edited plot title. It is the
positive fixture for the reader that splits one Session 27–44 draft into
drawings, and, with its Linear Legend row stroke and colors, for the Sessions
40–45 path of the Load reader above (`savedResultEdits`); `TESTA.gb` and the
Depth TSV are stored once for both modes.
`inactive-class-m.v44.gbdraw-session.json.gz` is a Web **Save Session**
download from the same commit: Circular has `TESTA.gb` and a Result, and Linear
has no inputs but an edited plot title, Accession, Length, and legend position.
It is the control that shows whether inactive-mode values survive. The steps,
inputs, deviations, and hashes are in `two-mode-project.provenance.json`.

`two-mode-thresholds.v42.gbdraw-session.json.gz` is a Web **Save Session**
download, kept unchanged, from first-parent `main` commit `3fd50841`
(Session 42). Linear holds `TESTA.gb` and `TESTB.gb` with the flat title
`FLAT_TITLE`; the Circular E-value and Identity thresholds, which Session 42
stores as a five-field `modeProfiles` entry, were edited while Circular held no
inputs. It is the positive fixture for the Session 40–42 reader. The steps,
inputs, and hashes are in `two-mode-thresholds.provenance.json`.

`mode-split-vectors.json` holds the shared cases of the split of a Session 27–44
draft into Session 46 mode slices: a fixture above or an inline flat draft, the
split context, and the JSON pointers the split writes (`expect`) or leaves to
Load (`expectAbsent`). `expectedModes` are the slices of the JavaScript split as
Web Load makes it (`tests/web/mode-split-vectors.test.mjs`); the Python split
(`tests/test_mode_split_vectors.py`) must equal them. `loadDefaults` are the
values Session 46 Load gives the pointers the split leaves absent. `loadExpect`
(node only) replaces an `expect` value where Web Load sends Python the value a
drawing read rather than the stored one (a carried field, such as `"20"` read
as `20`); the Python split of the stored draft keeps `expect`. A Session
27–39 case has no `expectedModes`: Web Load reads such a Session's settings from
its request, which the Python split does not.

`scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz`,
`scale-interval-negative-linear-api.v44.gbdraw-session.json.gz`,
`scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz`,
`scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz`, and
`scale-interval-zero-circular-web.v33.gbdraw-session.json.gz` hold a scale
interval of 0 or less, written by the first-parent `main` commit `fe6861f0`
CLI (in `diagramOptions.config`) and typed API (in `configOverrides`, where
`main`'s Web Save also writes it), by the release tag `0.13.0` CLI (in
`cliInvocation`), and by the first-parent `main` commit `b05a6bb8` Web Save
(the flat `configOverrides.scale_interval` of Session 31–3x). Those writers drew such a value as the automatic interval.
They are the positive fixtures for the request decoder, which still reads it
as automatic while the CLI and the Python option objects reject it as input
(D-04, `test_scale_interval_domain.py`). The commands
and hashes are in `scale-interval.provenance.json`.

`precomputed-ring-unlabelled-cli.v44.gbdraw-session.json.gz` and
`precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz` are Circular CLI
Sessions with one `--conservation_blast` ring and no `--conservation_labels`,
written by first-parent `main` commit `fe6861f0` and by release tag `0.13.0`.
Those writers stored no ring label and drew the table's file name. They are the
positive fixtures for the readers that keep that label after the default
became the file name without the extension (D-03,
`test_losat_circular_rings.py`). The commands and hashes are in
`precomputed-ring.provenance.json`.

`legacy-font-size-zero-linear-cli.v30.gbdraw-session.json.gz` and
`legacy-font-size-zero-circular-cli.v30.gbdraw-session.json.gz` are Linear and
Circular CLI Sessions written by release tag `0.13.0` with `--label_font_size 0`
(and, Circular, `--legend_font_size 0`) in `cliInvocation`. That writer drew
them. They are the positive fixtures for the legacy replay, which reads a
stored font size of 0 or less as Auto with the Web's notice (OV-1282, D-PYB-7,
`test_stored_font_size.py`). The commands and hashes are in
`legacy-font-size.provenance.json`.

`cli-linear-tables.v40-early.gbdraw-session.json.gz`,
`cli-linear-tables.v40.gbdraw-session.json.gz`,
`cli-linear-tables.v41.gbdraw-session.json.gz`,
`cli-linear-tables.v42.gbdraw-session.json.gz`, and
`cli-linear-tables.v44.gbdraw-session.json.gz` are CLI Sessions of the
`linear_tables` case of `tests/fixtures/cli_session_cross_surface/cases.json`,
written by first-parent `main` commits `8228ffab` (Session 40 with a full
option-derived `config`), `4ea96685` (Session 40, `config: {"adv": {}}`), and
`4e8c9380` (Session 41, the same), `3fd50841` (Session 42, request schema 7, no top-level
`config`), and `fe6861f0` (Session 44, request schema 8, feature catalog schema 4)
from the repository root. The first one's
schema 3 catalog stores no source feature index. They are the positive fixtures
for reading a CLI-written Session 40 or 41 `config` as no draft, for loading the
request-only Session 42 and 44 files, and for inferring the index (OV-269). The v40-early file is also the only positive fixture for a
Session 40+ Result without composition metadata, which Load gives the legacy
composition (OV-273); as a CLI Session it holds no manual positions (no
`ui` offsets). Commands, inputs, and hashes are in each `.provenance.json`.

`legend-renamed-rows.v46.gbdraw-session.json.gz` is the Session that the Web
app of first-parent `main` commit `9839de8e` wrote (bytes unchanged; see its
`.provenance.json`) after renaming the Legend rows GC skew (+) to `Z` and GC
content to `GC skew (+)`, coloring the latter red, adding a stroke to `Z`, and
renaming GC skew (-) to `W`, coloring it blue, and deleting it. Session 46
stored these edits under the shown captions (`GC skew (+)`, `Z`, `W`). It is
the positive fixture of the Session 47 reading that moves them to the row keys
(`GC content`, `GC skew (+)`, `GC skew (-)`): `tests/test_session_io.py`
upgrades and re-saves it, and `tests/web/mode-split-vectors.test.mjs` loads
it. `legend-row-key-vectors.json` holds the cases that the Web
(`legendEditsByRowKey`) and Python (`legend_edits_by_row_key`) readers share.
The Session's decompressed SHA-256 is
`8efd7a20948b7dcc9e5cb815d2391178d2b7a48c22bf37ece64160385340237a`.

`rule-captions-with-color-id.v46.gbdraw-session.json.gz` and
`rule-captions-with-color-id.specific-colors.tsv` are the Session and the
downloaded color-table TSV that the Web app of first-parent `main` commit
`9839de8e` wrote for one caption in two colors (see their `.provenance.json`),
bytes unchanged. That build saved the captions as `Complex IV [#ffe700]` and
`Complex IV [#1f77b4]`, and the Session holds a Legend stroke stored under the
first text. They are the positive fixtures of the reading that drops a
caption's own color ID (PD-OI-042 scenario revision 2):
`tests/test_specific_color_captions.py` reads both, and
`feature-identity-overrides.playwright.spec.js` loads the Session. The
Session's decompressed SHA-256 is
`0e986cee60a78fab2dffa49b9ea0eebcff34039be512ec26aa5d2397a658ffab`.

`font-family-preserved.v44.gbdraw-session.json.gz` and
`font-family-preserved.v46.gbdraw-session.json.gz` are Sessions written by
first-parent `main` commits `fe6861f06d0eaaecac0ba4ee6b1c3095019789df`
(Session 44) and `9839de8ee` (Session 46): the Python API wrote a request whose
config overrides name `objects.text.font_family`, and the Web app of the same
commit loaded it and saved it again, keeping that path as a preserved override
in its draft. `font-family-preserved.provenance.json` records the writers and
the decompressed SHA-256 values.

`BGC0000708-BGC0000713.v29.gbdraw-session.json.gz` is the unchanged Session
`tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json` of first-parent
`main` commit `10d3a3d28`, the only tracked version 29 Session of a size for a
fixture, compressed with `gzip -n -9`. It is the version 29 base of the
admission description's soundness test (`tests/test_session_contract.py`). Its
decompressed SHA-256 is
`71b923ce1a2b8fb587ed394e51571e5b15f289df58574fbc84734a1b2042b220`.

`BGC0000708-BGC0000713.v47.gbdraw-session.json.gz` is the Gallery Session of
first-parent `main` commit `150d26706` (bytes unchanged, `gzip -n -9`; see its
`.provenance.json`). Its record 5 is reverse complemented, and Session 47
stored that record's LOSATP table (`comparison-canonical-protein-4`) and its
`orthogroupResult` members in the drawn frame.
`majanivirus_collinear-rc.v47.gbdraw-session.json.gz` holds two records of the
majanivirus Gallery Session, the second reverse complemented, with the
`collinearityResult` that `main`'s Web converter (the Python of
`python-helpers.js` at `150d26706`, run natively) computes from that pair's raw
LOSATP entries, stored as `main`'s Web Session writer lays it out;
`tools/build_main_collinear_v47_fixture.py` rebuilds it from a worktree at
`150d26706` (command in its `.provenance.json`). Both are the positive fixtures of the
Session 48 conversion to the search frame (OV-1145): `tests/test_session_compat.py`
replays each and compares it with its `.v47.replay.svg`, the SVG the dev base
`e8b460772` drew (equal to `main`'s replay) with the Linear scale bounds since added to its composition
metadata; the BGC replay was redrawn after the
Linear label packing changed its record rows and again after the label boxes
added white boxes to its turned labels, whose text is unchanged, once more for float last bits when one
owner placed the Linear records, when `above_feature` labels kept their row gap, and when the
gradient Legend block took its title's width and the plot title its painted width, and both
replays when the composition metadata recorded the Linear scale bounds and when the
Legend took the drawn text width and once more for the OV-165 labels, ticks, definition and canvas change (`.v47.provenance.json`), and `tests/test_api_session.py`
upgrades and re-saves the BGC Session.

`BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz` and
`BGC0000708-BGC0000709-rc.v40-early.gbdraw-session.json.gz` hold the first two
records of the BGC Gallery Session, the second reverse complemented, with the
pairwise LOSATP table that `main`'s Web writer computes from that pair's raw
LOSATP entry (the Python of `python-helpers.js` and the JavaScript of
`session-request.js` `addCanonicalTable` at `cdd310134` (Session 39) and
`7aad9e3e0` (Session 40 before `3bca0e8d9`), run natively), stored as
`comparison-resolved-1`, with the Result of `main`'s CLI replay.
`tools/build_main_pairwise_rc_fixture.py` rebuilds each byte for byte (command
in its `.provenance.json`). The table names the reversed record's drawn-frame
features in `*_feature_svg_id` and has no `*_view_feature_svg_id` column; both
are the positive fixtures of that reading (OV-1183) and of the Session 39 Web
Load conversion (OV-1181): `tests/test_session_compat.py` converts and replays
them.

`file-card-imports-circular.v48.gbdraw-session.json.gz` is a Web **Save
Session** download of the current writer (Session 48, OV-1280, D-ORIGIN-2):
Circular holds no genome, the five file cards, and hand edits made before and
after each load. It pins the rows each card marks, the Default colors and
Blacklist words in `config.fileImports`, and the rows a hand edit made the
user's. The steps and inputs are in `file-card-imports.provenance.json`.

`HmmtDNA-hearbal-garden-cli.v48.gbdraw-session.json.gz` is the CLI sidecar that
release v0.14.2 wrote for `-p hearbal_garden`, the former name of the palette
`herbal_garden`, compressed with `gzip -n -9`. Its command, input hash, and
decompressed SHA-256 are in `HmmtDNA-hearbal-garden-cli.v48.provenance.json`;
`tests/test_palette_name_alias.py` replays it.

`record-identity-circular-single.v48.gbdraw-session.json.gz` is a Web **Save
Session** download, kept unchanged, from release tag `v0.14.2` (Session 48).
Circular draws record 2 of a two-record GenBank File alone (Multi-Record Canvas
off), and four features of that record carry a "This feature only" fill, a
label edit, Feature visibility Off and a stroke, all keyed by the record key
0.14.2 gave a record drawn alone (`circular-BGC0000709-#2`). It is the positive
fixture of the Generate that moves a record's rows to the key of its new view
(OV-1949): `record-identity-moves.playwright.spec.js` loads it, turns the canvas
on, and expects the four edits on record 2. The steps, edited features
and hashes are in `record-identity-circular-single.provenance.json`.

`two-mode-linear-shown.v48.gbdraw-session.json.gz` is a Web **Save Session**
download of the released v0.14.0 app (Session 48, tag `v0.14.0`) with both
modes used and Linear shown: Circular holds a genome, a Default colors file
and a changed setting, and has a Result that Session 48 does not keep;
Linear holds two genomes, a LOSAT comparison (one `losatCache` entry), a
Blacklist file and its Result. It is the positive fixture of the move of a
Session 48 to drawings (`tests/test_session_compat.py`,
`drawings-vectors.json` `upgrades`). The route, inputs and checks are in
`two-mode-linear-shown.v48.provenance.json`.

The four `drawings-*.v50.gbdraw-session.json.gz` files are Session 49
(`drawings`, DESIGN §1).
`gbdraw.session_drawings.write_session_drawings(version=49)` wrote each from a
fixture above, upgraded to Session 48 first, leaving out the resources no
drawing names: `circular-only` (one Circular drawing,
`font-family-preserved.v46`), `circular-linear` (both Gallery drawings of
`tests/utils/two_mode_session.py`, Linear active), `linear-pair-losat`
(`q-frame-main-web-losatn.v42` drawn twice: `linear` keeps its LOSATN entry,
`linear-2` has none) and `per-drawing-cards` (`file-card-imports-circular.v48`
drawn twice: `circular` keeps the five file cards, `circular-2` has none).
Each `.provenance.json` names the source, the edits, the script and the
decompressed SHA-256; `tests/test_session_admission_50.py` admits them and
writes each again unchanged.

`hidden-definition-circular.v49.gbdraw-session.json.gz` is a Web **Save
Session** download of the current writer (Session 49, PD-OI-109), written by
`tests/web/circular-hide-definition.playwright.spec.js` and recompressed with
`gzip -n -9`: Circular draws `tests/test_inputs/HmmtDNA.gbk` with the plot title
"Mitochondrial genome" at the top and **Show Definition** off, so its slice holds
`adv.show_definition: false` and its request
`objects.definition.circular.show: false`. Its decompressed SHA-256 is
`5d25d11dea37f278ff2aee041d9d7c92c63c7756bc1e4291ef9a118c91e33a94`.
