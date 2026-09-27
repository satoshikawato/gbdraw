# S05 shared error-boundary handoff to the export owner

This is a non-normative acquisition and work-order record. It changes no Product
receipt, authority, runtime or branch belonging to the export session. Read it
with [S05 result](SESSION_05_RESULT.md), the original export owner transfer at
`f0c8128ac9d252bbc1b64074cf995f536f789d03`, and the current accepted contract.

BUG-15/19 has one implementation owner: `fix/issue-601-bug15-bug19`, S02–S05.
Integrated coverage is E01–E07, R01–R06, U01, P01, SCI01 and G01 in the result.
The export session's transferred S03 is not resumed as a duplicate error/regex
implementation. PDF-specific producers, fonts, manifests/assets, snapshot and
filename handling, licensing/distribution and BUG-07 acceptance stay with the
export owner.

## Acquire the published result before shared-file editing

After the S05 result's final verification and ordinary same-named push, fetch
`origin/fix/issue-601-bug15-bug19` from the export session's own isolated clone.
The S05 implementation SHA is the enclosing commit, resolved without relying on
this document's mutable branch tip:

```sh
git log -1 --format=%H origin/fix/issue-601-bug15-bug19 -- docs/internal/issue-601-bug15-bug19-implementation-20260926/SESSION_05_RESULT.md
```

Its exact pushed SHA and matching local/remote/tracking verification are supplied
in the final S05 report and
`/tmp/gbdraw-issue601-s05-lPG7Gv/evidence/SESSION_05_PUBLICATION_HANDOFF.json`.
Read the complete result and any newer handoff before synchronizing by normal
merge. The export branch as a whole was not imported by S05; S05 did not edit or
push that branch. No reset/rebase/amend/force push is needed.

## Shared owners and consumption

| Shared boundary | Current owner / use |
| --- | --- |
| Python producer and native exception identity | Existing domain producer supplies finite cause facts. PDF/font-specific producers stay with export; native CLI/API capture must remain compatible. |
| `web_support/error_adapter.py`, `app/python-helpers.js` | One bounded adapter for real render/helper causes. Use existing helper wrappers before Python exceptions can become strings; retain code/operation/stage/allowed context and bounded cleanup facts. Do not add another raw serializer. |
| `workers/diagram-generation-worker.js`, `services/diagram-generation.js` | Existing one lazy diagram Worker/client transports the same model. Preserve typed operation names, initialization/staging stages, resource ownership, cancellation/currentness and prepared reuse. No second runtime or raw protocol. |
| `services/error-normalization.js` | Sole public wording/model owner. Feed finite producer facts; retain unknown actual stage and stable code. Existing PDF/export codes/actions remain here. Any needed common code/action extension belongs at this owner with focused coverage, never a PDF-local classifier or raw substring fallback. |
| `app/app-setup.js`, `index.html` | Existing composition and OperationError UI. Existing SVG/PNG/PDF caller callbacks own Retry/Use SVG and snapshot choice; reuse initially collapsed Details, explicit displayed-only Copy, manual selection and concrete recovery. Preserve compact Editor/Status and Session/Result boundaries. |
| Generate/Align and Color rule recovery | Existing run-analysis/transaction and `feature-editor/rule-actions.js` plus `pattern-drafts.js`. Export consumes accepted state/current Result and must not move draft, canonical request, Result admission or History ownership. |

The shared-file order remains S02 → S03 → S04 → S05 verified push → export's
shared PDF integration. Acquire the exact pushed S05 commit first, inspect newer
dev and both branches' actual differences, then coordinate editing of the same
files in sequence. This document is not evidence that another process stopped,
a distributed lease, or permission for simultaneous shared-file edits. Same-named
remote branch writers remain exclusive. Independent PDF-only work may proceed
under its own original authorization; S05 does not start it.

## Independent requirements and limits

PD-OI-046 (`web.errors.diagnostic-disclosure`, scenario 1, Choice A) and PD-OI-047
(`web.rules.rejected-pattern-edit-recovery`, scenario 1, Choice A) retain every
original receipt field. The original export receipt is not retired, remapped or
replaced by a second active diagnostic concern. Preserve each independent
requirement, including known corrections, actual stages, private/local-only
processing, native capture, initial/previous Result distinction, primary cause
through cleanup, draft/request/direction/History/retry/cancel isolation,
optional initially collapsed Details with always-visible summary,
keyboard/390px, displayed-only manual Copy and Clipboard fallback. Preserve
Python Color/Label versus JavaScript Search, field Not applied/Retry/Revert,
accepted state separation, one preparation/Worker, nonpersistence and lifecycle.

PDF receipt requirements remain separate: original Unicode/text, selective
font/style, positions/geometry, full features/legend/comparisons, captured
snapshot/filename, PNG DPI, lazy font/retry/glyph validation, searchable/selectable
embedded-font output and distribution/license evidence. S05 acceptance of common
errors does not certify these PDF outcomes or relax their acceptance.

Inherited limits: original BUG-19 audit pattern/build unidentified; standalone
module-fetch failure still requires reload; metadata-free Session and
Legend-override constraints unresolved. S05 local success is not supported-version
matrix, remote CI, human review, exact-dev staging, release or deployment proof.
No PDF runtime PR/integration or Issue #601 closure follows from this handoff.
