# Session compatibility fixtures

`BGC0000708-BGC0000713.v39.gbdraw-session.json.gz` is the unchanged session
JSON from first-parent `main` commit
`17e2c9dee32724219aa7c96a02d183280ffbe438`, compressed with `gzip -n -9`.
Its decompressed SHA-256 is
`9407365a3d5684490f1e72d81826101e661b9c0c362ec0b4e0b164af8e1b50b5`.

The schema-v2 fixture is the older supported compatibility case. Its expected
projection is recorded in the adjacent `.expected.json` file.

`feature-style-v40-gallery-minimal.gbdraw-session.json` is a hand-minimized
subset of the version-40 Gallery session at first-parent commit
`c71d998b5362e556ce003ab5c615127dc7fca577`. The source file's SHA-256 is
`5e566c305535bb4caccf50723158462471ecf82ce5f71ecc2fda57f91e80f937`.
The fixture retains the stale canonical default-color table, the matching
editor table, two broad rules, their derived overrides and captions, the
schema-3 catalogue identities, and the corresponding saved SVG fills. The
adjacent expected file is the shared normalized Web/Python projection.

The `.ambiguous-mutation.json` fixture describes the single negative mutation:
duplicating the saved default-color resource. A reader must reject that input
because neither resource is uniquely authoritative.
