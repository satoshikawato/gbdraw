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
