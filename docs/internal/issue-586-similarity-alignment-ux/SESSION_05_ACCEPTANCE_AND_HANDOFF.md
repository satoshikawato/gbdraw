# Session 05 instruction prompt: final acceptance and implementation handoff

## Mission

Audit the complete Issue #586 implementation against the signed Product
outcomes, architecture ratchets, acceptance criteria, and repository gates.
Fix in-scope defects, remove superseded paths, collect reproducible evidence,
and leave the implementation branch ready for human review. Do not create or
edit a pull request unless separately authorized.

## Branch and prerequisites

Use only `issue-586-similarity-alignment-ux-20260924`. Fetch `origin`, verify the
branch/upstream, and confirm:

- Session 00 authority is merged into the branch's `origin/dev` ancestry;
- Sessions 01–04 are present and recorded in the ledger;
- no unresolved session blocker remains;
- unrelated working-tree changes are identified and preserved.

When this prompt is supplied as the session request, it authorizes in-scope
fixes, one final acceptance commit if changes are necessary, and push to the
same-named branch. It does not authorize a PR, merge, deployment, tag, release,
or direct push to `dev`/`main`.

## Read and audit scope

Read the master plan, all session prompts and ledger evidence, the accepted
Product Contract outcomes, architecture/product ratchets, repository guidance,
and the entire branch diff from the accepted authority base. Do not accept a
session's completion claim without inspecting the code and evidence it names.

Review separately:

1. production Python and Web code;
2. tests and fixtures;
3. user/internal documentation;
4. generated Gallery/session artifacts;
5. dependency, CSP, offline, and Worker impact;
6. reference-output impact.

## Product acceptance matrix

Demonstrate every master-plan acceptance criterion:

- immediate busy feedback and single-flight invocation;
- exact reference identity from popup and drawer;
- one Web `Align…` action and one always-open review surface;
- every target row, deterministic default, visible finite reason, replace, and
  Skip;
- local no-Worker candidate/orientation/canvas edits;
- independent per-record orientation and the full known/unknown strand matrix;
- one Python-validated atomic Apply;
- correction/retry without draft loss;
- Cancel/failure/stale preservation of diagram, active plan, last Result, and
  History;
- schema-2 Session/regeneration/repair/Reset/Undo/Redo/reorder behavior;
- strict CLI ambiguity errors and typed Python fully resolved plans;
- keyboard, focus, screen-reader state, canvas interaction, and 390-pixel
  usability;
- explicit non-support for Collinear controls, anchor TSV, scores, support
  counts, and multi-hop inference.

Test recommendation stability under biological-record reorder and presentation
changes. Verify the default candidate is never recalculated from SVG geometry or
JavaScript sorting. Verify a recommendation remains transient until Apply.

## Architecture acceptance

Produce concise owner/path evidence for the non-increasing architecture change:

- one Python semantic owner for eligibility, ordering, recommendation,
  orientation effect, and final validation;
- one JSON adapter without policy;
- one Web controller and local draft;
- one Worker/resolver route and one canonical Result/History admission path;
- one current schema writer;
- only evidence-backed released compatibility readers;
- removed separate Align/Align & orient, global Web mode, auto-apply,
  zero-selection requirement, and discarded-on-retry behavior.

If the ratchet's exception conditions are triggered, supply the complete OE,
PE, and CB evidence required by policy. Never waive a failing Gate.

## Required verification

Run focused tests first, fix failures, then run the full gates. At minimum:

```bash
python -m pytest \
  tests/test_similarity_alignment.py \
  tests/test_similarity_alignment_rendering.py \
  tests/test_similarity_alignment_web_adapter.py \
  tests/test_session_request_codec.py \
  tests/test_session_compat.py \
  tests/test_api_session.py \
  tests/test_session_io.py -v

node --test \
  tests/web/similarity-alignment-actions.test.mjs \
  tests/web/right-drawer.test.mjs \
  tests/web/session-request.test.mjs \
  tests/web/session-authority.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/history-config-restore.test.mjs \
  tests/web/diagram-generation-worker.test.mjs

npx playwright test \
  tests/web/similarity-alignment-ui.playwright.spec.js \
  --project=chromium --workers=1

ruff check gbdraw/
node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
python tools/update_cli_reference_help.py --check
python -m pytest tests/ -v -m "not slow"
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python -m build
git diff --check
```

Monitor long tests incrementally and allow the repository-prescribed timeout.
Do not update reference outputs merely to make a comparison pass. Prepare the
gitignored browser wheel when the browser/offline harness requires it.

Run Gallery reproduction/parity checks for changed generated artifacts. Run an
offline-browser audit only if the implementation changed runtime dependencies,
privacy, bundle composition, or Worker lifecycle; ordinary UI changes alone do
not require that skill.

Perform a final real-browser visual review at desktop and 390-pixel widths.
Confirm candidate details, reasons, effective orientation, busy/error states,
focus, pan/zoom, Apply/Cancel, and retry are legible and operable.

## Final cleanup and ledger

- Search for stale user-visible `Align & orient`, old automatic-apply logic,
  global orientation-mode branches, schema-1 current readers/writers, duplicate
  recommendation logic, and obsolete tests/docs.
- Confirm no generated wheel, cache, temporary screenshot, build directory,
  or unrelated file is staged.
- Update the master ledger with exact commits, commands, results, visual review,
  Product Impact classification, architecture evidence, and any accepted
  residual risk.
- Do not mark the work complete if a required gate or acceptance item is
  missing.

If fixes are required, use an English commit title such as:

```text
Complete Issue 586 similarity alignment acceptance
```

If no tracked change is required, do not create an empty commit. Push only the
implementation branch when it has a new authorized commit.

## Handoff

Report:

- final branch and head commit;
- authority merge commit and implementation session commits;
- user-visible outcome;
- semantic owners and canonical paths;
- schema/compatibility result;
- focused, browser, full-suite, architecture, Product, build, Gallery, and
  output-comparison evidence;
- generated artifacts and visual review;
- remaining risks or blockers;
- an English proposed PR title and concise maintainer-facing summary.

Do not open or edit the PR without explicit permission. If later asked to do so,
load the repository's `write-clear-pull-request` skill and run the required
language check before the remote mutation.
