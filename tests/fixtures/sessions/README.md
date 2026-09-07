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
