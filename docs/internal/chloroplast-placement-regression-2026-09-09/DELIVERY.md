# Session 00 delivery verification

Starting `origin/dev`: `991f8eef6734060d4e8d0d2e6b4294309f57e413`.
Version: `0.14.0`, unreleased.

## Retained fix and ancestry

The existing fix was located on `fix/chloroplast-feature-placement` in
`/tmp/gbdraw-chloroplast-fix`. Its HEAD and merge base with the freshly fetched
`origin/dev` both equal the starting SHA above. The fix was uncommitted; the
original checkout and the unrelated dirty primary checkout were preserved.

A clean worktree at `/tmp/gbdraw-session00-delivery` was created from that exact
`origin/dev`, with branch `fix/session00-chloroplast-delivery` and no upstream.
All 26 retained changed/untracked files were copied after checking SHA-256
identity. No rebase conflict or production adaptation was needed.

The retained regression note, before/after images, final test logs, mutation
sensitivity results, and generated artifact audit correspond to this same base
and fix. The delivery inventory and copied final logs are retained locally in
`/tmp/gbdraw-session00-evidence/`. Production files and generated artifacts are
byte-identical to the existing fix. The only test extension loads the saved
chloroplast session in a fresh browser context; its completion predicate and
180-second bound match the existing initial cold-import check.

## Verification

| Check | Result |
| --- | --- |
| Chloroplast real-data combinations | 8 passed, including reverse-direction placements |
| Placement/radial/composition and documentation contracts | 379 passed, including the 8 chloroplast cases |
| Tracked reference SVG comparisons | 16 passed; references unchanged |
| Web session/history/placement/composition | 18 test files passed |
| Architecture contracts | 137 passed |
| Chromium placement/history/replay | 13 existing scenarios passed |
| Chloroplast cold context and fresh-context restore | 1 passed after adapting the import wait to the fresh context |
| Retained complete fast Python suite | 3969 passed, 17 skipped on the identical production fix and base |
| Ruff and whitespace | Passed |
| Web change policy | Gate PASS; Review REQUIRED for session/output effects |
| Browser wheel | Built as `gbdraw-0.14.0-py3-none-any.whl`; generated and untracked |

The chloroplast browser regression blocks external requests and records zero
attempts. It uses Chromium through Playwright 1.61.1, a 1600 × 1000 desktop
viewport, and a 390 × 844 mobile viewport. It checks `Separate Strands`
on → off → on, labels both → outer → both, repeated generation, placement
availability before Generate after both imports, preserved placement overrides,
and the mounted legend/record bounds. The fresh-load desktop screenshot and
mobile control screenshot were inspected. Mobile coverage establishes the
restored controls and enabled Generate action, not full diagram visibility in
the initial mobile viewport.

Local logs and screenshots: `python-focused.log`, `reference-svg.log`,
`web-session-history.log`, `architecture-contracts.log`, `browser-regression.log`,
`browser-fresh.log`, `browser-fresh/`, and `architecture-policy.log` beneath the
delivery evidence directory. The first browser run passed 13 scenarios and
exceeded the old five-second restore assertion in the new cold context; the
final single-scenario rerun passed in 24 seconds. No production change was
needed for that test-harness adaptation.

The refreshed Gallery session changes only `results` and `runMetadata`.
Its canonical request, resource bytes, and editor state are unchanged. The
three chloroplast artifact-manifest size/hash entries match the retained files.
Existing generated documentation images and their earlier pixel review are
reused; no new Gallery or release-wide audit was performed.

## Ownership and product scope

The five outcomes are explicitly selected by the maintainer's Session 00
instruction. This delivery retains the existing placement, radial-layout,
radial-label, circular assembly, and session/config owners. No second owner,
canonical path, geometry cache, schema, migrator, or dependency is introduced.
The base session authority already projects `runMetadata` as artifact state;
the session coordinator now retains its existing `trackSlotGeometry` field.
Registered architecture subjects and Product Impact authority are unchanged.
The ordinary non-increasing architecture evidence and rollback are recorded in
the original regression note. The policy's Review REQUIRED result is retained
for maintainer review; it is separate from the passing executable Gate.

Session 01, Session 02, LOSAT integration, S11/S12, version changes, release
baseline assignment, and publication are outside this delivery.
