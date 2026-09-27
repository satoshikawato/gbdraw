# S05 — guard integration follow-up

The original S05 acceptance mapping is retained in
[`SESSION_RESULTS/S05.md`](./SESSION_RESULTS/S05.md). This follow-up
records the separate guard prerequisite and the subsequent trusted-base
incorporation. It does not change the four approved outcomes or add runtime.

## History and actual authority

- Original S05 runtime: `9c6a15288cc40d2ddd8496b2b32e2b858e11c329`.
- Guard implementation: `6c98cf713a2dd6791565e6d64dab9a77bc551543`.
- Guard branch authority incorporation: `77105b13009830d5fb951ece097afba532f67969`.
- Separate guard PR: https://github.com/satoshikawato/gbdraw/pull/616
- Guard dev merge: `c922fc38aac78da9be83342c09ac0164ecef6ff6`.
- Accepted #614 incorporation into #600: `b7a1aa19aa5fad7ee253591c2a6f5d3425fd521d`.
- Latest actual authority-integrated base: `c922fc38aac78da9be83342c09ac0164ecef6ff6`.
- Normal #600 incorporation of the integrated guard: `cb2ad027d3f60fde2331a23ad8f9f6274fe863ba`.

All 18 scoped runtime acceptances and required local gates are satisfied.
The original two architecture failures are resolved; this does not assert
runtime PR review, runtime hosted CI, release or dev admission readiness.

All S00–S05 history is preserved. Shared checkout branch and unrelated files
are unchanged. Authority revision 22 / 47 active IDs / 25 receipts are current
observations, not pinned acceptance counts. PD-OI-040–043 retain all nine
receipt fields, five outcome clauses and complete human approval records.
No inferred rationale, retirement or residual risk is added.

## Results and valid reuse

- Guard suite on original latest dev: 192 passed, exit 0.
- Guard suite on fixed original S05 plus test-only overlay: 192 passed, exit 0.
- Guard suite after the #615 authority-only merge: 192 passed, exit 0.
- Integrated #600 required architecture/ratchet suite: 192 passed, exit 0 (`E/runtime-integrated-architecture.log`).
- Formal fast Node contracts after #614: 843 passed, exit 0.
- Four combined annotation/style journeys plus eleven alignment UI journeys:
  15 passed, exit 0 (6.7m), port 42221, workers 1 / retries 0.
- Formal PR smoke: 13 passed, exit 0 (2.9m), port 42222.
- Required non-slow Python gate on the incorporated shared UI: 6642 passed, 17 skipped, 11 deselected, 23 warnings in 801.16s; independently recorded pytest return code 0 (`E/runtime-python-supervised-exit.json`).
- Guard PR required hosted checks: all seven full-plan jobs and both required protected statuses SUCCESS for guard head `77105b13009830d5fb951ece097afba532f67969`; normal merge confirmed.
- Final exact-commit trusted-base gate: recorded externally in
  `/tmp/gbdraw-issue600-guard-integration-evidence/runtime-final-gate.log` and
  the handoff, after the report commit. No self SHA appears in this report.

Executed selections (all from the isolated runtime clone):

```bash
PYTHONPATH=. GBDRAW_WEB_TEST_PORT=42224 python -m pytest tests/ -v -m 'not slow'
node --test tests/web/architecture-contracts.test.mjs \
  tests/web/architecture-ratchet-fixtures.test.mjs \
  tests/web/product-impact-ratchet-fixtures.test.mjs
find tests/web -maxdepth 1 -type f -name '*.test.mjs' \
  ! -name 'architecture-contracts.test.mjs' \
  ! -name 'gallery-session-publication.test.mjs' -print0 | xargs -0 node --test
PYTHONPATH=. GBDRAW_WEB_TEST_PORT=42221 npx playwright test \
  --config=playwright.functional.config.js \
  tests/web/annotation-style-integration.playwright.spec.js \
  tests/web/similarity-alignment-ui.playwright.spec.js \
  --workers=1 --retries=0 \
  --output=/tmp/gbdraw-issue600-guard-integration-evidence/runtime-browser
PYTHONPATH=. GBDRAW_WEB_TEST_PORT=42222 npm run test:web:pr-smoke
```

Evidence root: `/tmp/gbdraw-issue600-guard-integration-evidence` (`E`).
The supervisor itself also exits 0 (`python E/run-pytest-supervised.py`).
The first Python run used `PYTHONPATH=.` and port 42223 on
`b7a1aa19aa5fad7ee253591c2a6f5d3425fd521d`, with source held fixed. It printed
6642 passed but its exec session returned 143; it is not labeled an exit-0 gate.
The sender/cause of that signal is unproven. The supervised rerun on
`cb2ad027d3f60fde2331a23ad8f9f6274fe863ba` uses port 42224 and records the
actual pytest child return code durably. No assertion, selection or test-owned
timeout is changed. The post-first-run dev incorporation is guard tests, authority
and internal documentation only; it changes no Python runtime, pytest owner,
Gallery input, dependency or browser runtime. Authority-sensitive guard
checks are executed again after incorporation.

Python production, native pytest owners, Gallery inputs, references and wheel
source bytes are unchanged from original S05. Its Ruff, read-only reference
comparison (16 passed), recipe (189/189), Gallery (103/103), and Gallery
publication browser parity (9 passed) are valid unchanged-path evidence.
The supervised full run also passes 189/189 recipe, 103/103 Gallery and 16/16 read-only reference cases, bound to exact node IDs in `E/recipe-coverage-supervised.json` and `E/gallery-coverage-supervised.json`; no extra subset invocation is claimed. The accepted #614 shared UI change is
covered by fresh Node, PR smoke, combined journeys, alignment UI and Python
browser-wrapper runs instead of being labeled unchanged-source evidence.

Each latest combined case verifies warning ownership, canonical captions and
rules, native geometry, import/live-edit Undo/Redo, disabled-draft lazy Load,
repeated Generate, original file bytes and actual downloads. All four actual
SVG downloads are byte-identical to the original readable, visually reviewed
S05 figures (`E/figure-equivalence.json`). No public figure is regenerated.
The original 18 acceptance IDs remain supported by their mapped tests plus
these updated integration runs; no acceptance is weakened.

The guard retains the two metadata admission callers and permits only the
pure named validator from the two request/Session callers. Legend replacement
ownership remains required with a five-operator ceiling; positive contractions
and negative growth/missing-owner fixtures pass. Detector, allowlist, workflow,
Product decision source and runtime owner/path remain unchanged.

## Remaining limits and rollback

GitHub auto-merge is disabled. The rejected auto-enable request made no merge;
actual remote state was inspected and the PR was merged normally after its
required checks passed. Repository settings and branch protections are intact.
No direct dev/main push, runtime PR, publication, deployment or tag was invoked.
The #600 runtime itself has not been merged into dev by this guard integration.

Physical assistive-technology/device tests, full functional staging,
performance and release acceptance are not claimed. The inherited Circular
Width/Radius `[object Object]` presentation issue remains outside the four
optional-pixel fields. Original S05 retained geometry, scalar types and bytes.

Tests, production merge diff, reports and generated paths are reviewed
separately; `git diff --check` passes. The test-induced LOSAT executable mode
is restored without byte changes before committing. Only the S05 status and
this follow-up report are staged for the completion record.

Rollback is a normal work-branch revert, preserving Product authority and
S00–S05 history. Reverting the guard PR would reinstate the original two
characterization failures without changing runtime behavior.

English commit title: `docs: record issue 600 guard integration`.
Summary: Record the separately integrated metadata/legend guard prerequisite
and successful S05 revalidation against trusted dev.
