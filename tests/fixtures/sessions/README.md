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
`feature-edits-circular-copies.v44.gbdraw-session.json.gz`, and
`feature-edits-circular.v33.gbdraw-session.json.gz` are Web **Save Session**
downloads, kept unchanged, from first-parent `main` commits `fe6861f0`
(Session 44) and `b05a6bb8` (Session 33). They hold Feature visibility, Label
visibility, and label text edits keyed by rendered feature ID: a cropped and
a reverse-complemented Linear record, two Circular records with the same ID,
and a Session 33 without a feature catalog. They are the positive fixtures for
the readers that move those edits onto source identities in Session 45. The
steps, inputs, and hashes are in `feature-edits.provenance.json`.
