# Session compatibility fixtures

`BGC0000708-BGC0000713.v39.gbdraw-session.json.gz` is the unchanged session
JSON from first-parent `main` commit
`17e2c9dee32724219aa7c96a02d183280ffbe438`, compressed with `gzip -n -9`.
Its decompressed SHA-256 is
`9407365a3d5684490f1e72d81826101e661b9c0c362ec0b4e0b164af8e1b50b5`.

The schema-v2 fixture is the older supported compatibility case. Its expected
projection is recorded in the adjacent `.expected.json` file.

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
rename must not change it. SHA-256:
`0ae2dedeb8dbb426dc1dc5fa8d66d5c9c0a3b1945deca39eb7a71db5a4560b1b`.

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
drawings; `TESTA.gb` and the Depth TSV are stored once for both modes.
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
values Session 46 Load gives the pointers the split leaves absent. A Session
27–39 case has no `expectedModes`: Web Load reads such a Session's settings from
its request, which the Python split does not.

`scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz`,
`scale-interval-negative-linear-api.v44.gbdraw-session.json.gz`, and
`scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz` hold a scale
interval of 0 or less, written by the first-parent `main` commit `fe6861f0`
CLI (in `diagramOptions.config`) and typed API (in `configOverrides`, where
`main`'s Web Save also writes it), and by the release tag `0.13.0` CLI (in
`cliInvocation`). Those writers drew such a value as the automatic interval.
They are the positive fixtures for the reader that keeps drawing it so after
fresh input rejects it (D-04, `test_scale_interval_domain.py`). The commands
and hashes are in `scale-interval.provenance.json`.
