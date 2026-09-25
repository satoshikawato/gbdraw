# Similarity Alignment auto-apply and Definition column regression — master plan

Status: implementation present on the fixed branch; Session 04 acceptance record below.
Implementation branch: fix/similarity-alignment-auto-apply-definition-column-20260925
Branch base at creation: origin/dev@ac067f9b5482b480a1edca983709e3321f75157f, fetched 2026-09-25.
Product outcome: Choice A in 00_DECISION_PACK.md was selected by the requester on 2026-09-25. The Product Decision Owner subsequently approved the complete four receipts, serialized in `b772359b` and merged into `origin/dev@69134e28` before dependent runtime changes.

## 1. Purpose and scope

gbdraw is a Python genome-diagram generator with a single-page Web app. In a Linear diagram, a user can select an exact feature in a Similarity Group and run Align. The Python resolver classifies the response as resolved or ambiguous. Current Web code opens Select alignment anchors in both cases. This follow-up changes the normal Align action so a resolved plan is applied directly, while an ambiguous response still opens the full review. A separate, visible Review alignment options… action at both existing reference entry points always opens the full review, including for a resolved plan. The explicit review retains Select, Skip, per-target Match reference direction, candidate facts, canvas selection, and batch Apply.

The same report identifies a separate rendering regression: with Lock Definition Column enabled, Definition text lines have a left anchor within each block, but Align moves the blocks by different record-specific x translations. The lock must keep the row/side Definition blocks on one common left edge while preserving the configured gap from the leftmost drawn sequence.

Scope includes the authority update, Python Linear layout correction, Web controller/UI change, documentation, focused tests, and final gates. No request or Session schema version, alignment recommendation algorithm, orientation rule, Worker owner, or CLI behavior is intended to change.

## 2. Fixed branches and authority gate

- All implementation, test, and documentation work in Sessions 01–04 uses exactly fix/similarity-alignment-auto-apply-definition-column-20260925. It was created with no upstream from freshly fetched origin/dev. Before each commit or push, verify branch, upstream, and target. Push only origin/fix/similarity-alignment-auto-apply-definition-column-20260925.
- Session 00 is the sole exception: durable Product authority must be committed on a separate authority-only branch made from the then-latest origin/dev. Its branch must contain no runtime, tests, or implementation plan. Obtain explicit approval of the complete four Choice A receipts in 00_DECISION_PACK.md before serializing; the requester has chosen A but has not supplied a complete PRODUCT_DECISION response. Do not infer or silently reuse the proposed owner/date/risk terms as a signed receipt.
- Merge of the authority-only change into origin/dev is required before Sessions 02 and 03 change dependent Web runtime. The implementation branch may be prepared and Session 01 may correct the independent Definition regression while this merge is pending. After authority merges, fetch origin/dev and merge it into the fixed implementation branch; do not replace this branch.
- No authority candidate authorizes its own dependent runtime. Normal PR/merge/deploy/tag permissions remain separate from the authorization to push this plan branch.
- Preserve unrelated working-tree content. In particular, docs/internal/gbdraw_meet_gbdraw_implementation_v2/ was untracked before this plan and is outside scope.

## 3. Existing authority and observable behavior

| Concern | Current accepted authority | Intended disposition |
| --- | --- | --- |
| Python anchor resolution | PD-OI-026 revision 2 | Revise review timing only; keep exact reference, stable identity, only-usable and unique-direct-RBH resolution, and ambiguous recommendation semantics. |
| Per-target orientation | PD-OI-027 revision 2 | Keep preserve as default and reverse only for an explicit match request with known opposite strands; make the explicit review action available before commit. |
| Web alignment surface | PD-OI-031 revision 2 | Supersede mandatory review for resolved default Align and allow a named explicit review action. |
| Web choice/retry | PD-OI-034 revision 2 | Supersede mandatory Apply for resolved default Align; retain editable draft and retry on review or automatic failure. |
| Linear Definition | PD-OI-024 revision 1 | Already requires one left edge with Lock on. Correct this regression under existing authority. |

The current public Web reference also says the palette always opens. Update it when the authority-gated runtime is implemented. Existing code/tests are evidence, not authority.

Acceptance semantics:
- Branch solely on the validated Python response status. A resolved response includes an explicit validated plan. An ambiguous response includes transient recommendations but no final plan. The visible 0 need a choice count is not an ambiguity signal.
- Normal Align plus resolved: immediately enter applying state, generate with that plan, then publish one Result/History replacement and a concise summary. Do not briefly mount the palette or run a redundant second resolver call.
- Normal Align plus ambiguous, or explicit Review action plus either status: present the full existing draft. Edits remain local until Apply, which validates the complete batch in Python before rendering.
- A failed automatic render opens the retained resolved draft with an actionable error and retry. A failed explicit Apply retains the draft. Failure, Cancel, stale/superseded work, and resolver errors preserve the previous Result and History.
- An exact reference never moves or reverses; skipped, missing, and unusable targets keep position and orientation; every target keeps vertical position. Idempotent center alignment and source-relative rev indication remain intact.
- Transient guide/badges and draft remain absent from saved Sessions and downloads. Existing plan regeneration, Reset, Undo/Redo, and released Session readers retain their meanings.

## 4. Architecture and smallest-change design

Canonical flow:

~~~text
popup or drawer exact reference
  -> one Web alignment controller in gbdraw/web/js/app/similarity-alignment.js
  -> existing Worker helper / Python resolver and validator
  -> validated resolved plan -> existing applyPlan / runAnalysis / Result admission
     or full local review draft -> existing batch Apply -> same applyPlan
  -> typed render request -> Linear assembly -> sanitized Result
~~~

- Python in gbdraw/layout/similarity_alignment.py owns candidate eligibility, ordering, recommendation, and typed plan validation. The Web controller owns operation state, branch selection, draft, retry, and stale-result admission. The template in gbdraw/web/index.html and wiring in gbdraw/web/js/app/app-setup.js only expose those actions; they do not decide biology or duplicate the state machine.
- Linear layout uses one shared definition-placement rule in gbdraw/layout/linear.py. gbdraw/diagrams/linear/assemble.py calculates final record positions and collision bands; gbdraw/diagrams/linear/builders.py paints groups using the same placement. For Lock on, derive a common row/side Definition origin after final translations, with enough space before the leftmost record. Do not apply each record translation to that locked origin. Local record labels above sequences still move with their records. Lock off still follows row translation. Check both ordinary and multi-record paths.
- SOLID: retain one semantic owner for each policy and keep UI entry points as thin callers. KISS: use response.status and the existing controller transitions; avoid a second planner or a generic workflow framework. DRY: reuse reviewRows, applyPlan, the canonical Result/History path, and shared definition placement for paint and collision. YAGNI: add no schema, global orientation mode, feature flag, compatibility reader, or parallel renderer for this follow-up.
- Under ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md, record concise owner/path evidence for a non-increasing change. If a defined exception is actually triggered, provide the full OE, PE, and CB sets. Under PRODUCT_IMPACT_RATCHET.md, compare the jointly required effects of the four superseded concerns, not just their option IDs.

## 5. Root cause and test gap for Definition

place_linear_definition returns a record-independent x when keep_left is true. After similarity alignment calculates final_translations in Linear assembly, both add_record_definition_group and _record_collision_bands add each record's translation_x to this locked result. The text inside each group is left-anchored, but the groups differ in x. Existing Lock browser tests check Generate without unequal record translations, so they do not catch this case.

Fix both paint and collision using the same final common column and the final leftmost record position. Protect the invariant with an SVG/Chromium assertion of each row/side Definition's actual left edge and configured gap. Include unequal positive and negative translations, one-record rows, mixed multi-record rows, and Lock off. Compare the generated figure after Align, manual record translation, save/load, and regeneration. Retain the current tests for name/subtitle line alignment, center alignment when Lock is off, and text-width measurement.

## 6. Sessions and completion gates

| Session | Instruction file | Scope | Start condition |
| --- | --- | --- | --- |
| 00 | SESSION_00_AUTHORITY.md | Exact receipt and authority-only supersession of PD-OI-026/027/031/034 | Can start now; dependent runtime waits for merged authority. |
| 01 | SESSION_01_DEFINITION_COLUMN.md | Shared locked-column placement and regression coverage | Can start now on the fixed implementation branch. |
| 02 | SESSION_02_AUTO_APPLY_CONTROLLER.md | Resolved automatic path, draft/retry/stale invariants, controller tests | Authority merged into origin/dev and then into the implementation branch. |
| 03 | SESSION_03_REVIEW_ENTRY_AND_DOCS.md | Named explicit review entry in popup/drawer, accessibility, public docs and browser tests | Session 02 complete and authority available. |
| 04 | SESSION_04_ACCEPTANCE_AND_HANDOFF.md | Integrated geometry, workflow, persistence, visual review, and required gates | Sessions 01–03 complete. |

Every session instruction is self-contained for a contributor who has not seen this conversation. At the end of Sessions 00–03, display the complete next-session instruction text in the response for the user to paste into a new Codex session. Session 00 may hand off to Session 01 before authority merges. Session 01 must identify the authority gate before handing off to Session 02.

Focused evidence:
- Python: tests/test_linear_definition_alignment.py, tests/test_linear_multi_record_layout.py, tests/test_similarity_alignment.py, tests/test_similarity_alignment_rendering.py, relevant typed request/Session tests.
- Web unit: tests/web/similarity-alignment-actions.test.mjs plus focused wiring/History tests.
- Browser: tests/web/similarity-alignment-ui.playwright.spec.js and tests/web/linear-multi-record.playwright.spec.js. Verify Node Playwright availability; if absent, use Python Playwright for equivalent targeted checks. Rerun a Chromium sandbox failure with the documented escalation.
- Cases: one usable candidate, unique direct RBH among multiple candidates, missing/unusable targets, genuine ambiguity with preselection, explicit review of a resolved plan, preserve/match direction, automatic failure and retry, cancellation, stale completion, repeated Align, History, Session, and actual SVG Definition alignment.

Final gates after focused failures are fixed: ruff check gbdraw/; node tests/web/architecture-contracts.test.mjs; node tools/check-web-change-budget.mjs; python tools/update_cli_reference_help.py --check; python -m pytest tests/ -v -m "not slow"; python -m pytest tests/test_output_comparison.py::TestOutputComparison -v; python -m build; git diff --check. Prepare the gitignored browser wheel when needed. Allow at least 30 minutes for long pytest runs and monitor incrementally. Review final desktop and narrow-viewport browser results. Run an offline bundle audit only if runtime dependencies, privacy, bundle composition, or Worker lifecycle change. Never modify reference outputs merely to silence a failed comparison.

## 7. Exit criteria and record

The fix is complete only when every accepted behavior above is demonstrated, Definition blocks share a measured left edge after unequal alignment translations, the explicit review action remains accessible, no required gate is failing, and the authority is on the runtime base. Review production, tests, docs, and generated diffs separately. Record exact commits, commands, results, visual observations, architecture owner/path evidence, Product authority, and any remaining risk here or in a session handoff. Preserve unrelated files and do not stage the pre-existing untracked directory.

Planning commit only: this master plan, the Decision Pack, and Sessions 00–04. No runtime source, tests, generated wheel, reference SVG, or Gallery asset belongs in the planning commit.

## 8. Session 03 review entry record

The named work branch contains merged Choice A authority
`origin/dev@69134e28` and Session 02 controller commit `86aa03d9`.
This session adds a visible **Review alignment options…** button beside
**Align…** in the exact-feature popup and Similarity Groups drawer. Both
pass their existing exact reference and group through `app-setup.js` to
the one `similarity-alignment.js` controller with `mode: 'review'`.
The controller, Python resolver, Worker, draft, Result/History admission,
Session schema, and transient canvas overlay owners are unchanged.
The normal resolved path still applies without mounting the palette;
ambiguity and explicit review open the same complete palette. This
implements PD-OI-026/027/031/034 revision 3, including the jointly
required resolved, ambiguous, explicit-review, and retry continuations.

The Web reference now distinguishes those paths. The BGC Gallery
tutorial names both actions, and its popup image was recaptured from the
restored Session using the supported capture tool. Node controller/drawer:
32 passed; Python resolver/Web adapter: 85 passed; Chromium alignment:
7 passed, with the final resolved-path automatic failure/retry check
rerun and passed; Gallery Chromium: 23 passed; documentation contracts:
12 passed; strict Gallery capture check, Ruff, CLI help, architecture
contracts, Web change-budget gate, and diff whitespace check passed.
Desktop and 390-pixel palette screenshots were inspected. The generated
browser wheel remains ignored, and no reference SVG or unrelated worktree
file was changed. Session 04 owns integrated geometry, persistence,
full-suite, and build acceptance.

## 9. Session 04 integrated acceptance record (2026-09-25)

Sections 1–7 retain the original implementation plan and its then-pending
authority gate. This record describes the verified implementation branch.

- Branch: `fix/similarity-alignment-auto-apply-definition-column-20260925`,
  upstream `origin/fix/similarity-alignment-auto-apply-definition-column-20260925`;
  fetched `origin/dev@69134e28` is an ancestor. Authority commit
  `b772359b` (PD-OI-026/027/031/034 revision 3, complete Choice A) was merged
  into that base by `69134e28`, then into this branch by `6e17ac35`.
- Implementation commits: `015d1257` (plan), `2092dee4` (Definition column),
  `86aa03d9` (automatic resolved application), and `7d42f3e7` (explicit
  popup/drawer review and public wording). No unmerged authority candidate
  was used as runtime authority.
- Normal Align applies Python's validated resolved plan with one helper and one
  successful Result/History action, without mounting the palette. Only-usable,
  unique-direct-RBH, missing, and unusable targets retain their Python
  decisions. Genuine ambiguity and either explicit review action open the
  same complete, preselected draft; Select, Skip, Match reference direction,
  and canvas selection remain available. Failed automatic render opens that
  retained draft for correction and Apply retry. Failed Apply, Cancel, stale
  completion, duplicate starts, and supersession preserve the last Result.
  Session save/load, regeneration, source-relative orientation, Reset,
  Undo/Redo, and export cleanliness are covered by the focused Chromium and
  unit journeys.
- Chromium SVG measurement on the regenerated five-record Streptomyces Gallery
  Session with manual x inputs `[-85, 40, -35, 75, 15]` and Lock on: all five
  transformed Definition left edges were exactly `16 px`; final sequence
  translations spanned `-303.35` to `+35.34 px`. The nearest sequence began
  at `593 px`, and the widest Definition ended at `572.09 px`, leaving
  `20.91 px` against the configured `20 px` gap. The separate four-record
  Chromium contract measured two-row and mixed three-row cases, Lock off,
  record-local labels, save/load, and regeneration; locked edge spread was at
  most `1 px`, and sequence clearance was at least `19 px` including browser
  measurement tolerance.
- The five-record Gallery Session regenerated with five species labels, gene
  colors, quantitative comparison ribbons, scale, legend, and title. The
  3200×1200 desktop capture was inspected at readable scale; the 390×740
  popup and review captures showed legible labels, exact reference, both
  actions, candidate reason, Apply, and Cancel within the viewport. The
  palette occupies `x=12..378 px` and `y=12..728 px` at 390×740 px and
  covers most of the figure while open. The list controls remain reachable;
  closing restores the preview, and desktop canvas selection, manual pan/zoom,
  and palette movement remain available. The Product Decision Owner selected
  `A / RETAIN_MOBILE_PALETTE_COVERAGE` for PD-OI-035 scenario revision 2,
  explicitly retiring the literal unobscured-preview requirement and
  simultaneous canvas selection and manual pan/zoom at 390 px while open.
  The revised Session 04 criterion requires legible and reachable palette
  controls and candidate list, with the preview restored on close. The exact
  receipt is serialized in local authority-only commit `5f3f1445` from
  `origin/dev@69134e28`; it has not been merged into the runtime base.
  Existing Chromium 390 px screenshots and browser checks support this
  documentation correction; no runtime, test, or figure bytes changed.
  Final authority-backed acceptance remains pending that merge.

| Gate | Result |
| --- | --- |
| Focused four-file Python command in Session 04 | 184 passed |
| `node --test tests/web/similarity-alignment-actions.test.mjs` | 28 passed |
| Requested two-spec Chromium command, one worker | 43 passed |
| Gallery tutorial and Session regeneration Chromium specs | 25 passed |
| Gallery JSON, strict operation capture check, and Gallery JS syntax | passed; 19 media/operations checked |
| Gallery manifest and published Session verification | passed after refresh |
| Clean-directory H-CLI-07 and tutorial-9 arrow reproduction tests | 1 passed; 2 passed |
| Embedded Chromium comparison contracts with Node Playwright on PATH | 2 passed |
| `ruff check gbdraw/` | passed |
| `node tests/web/architecture-contracts.test.mjs` | 137 passed |
| `node tools/check-web-change-budget.mjs --base origin/dev` (full working tree) | Gate PASS, Review REQUIRED for the reviewed BGC Session artifact |
| `python tools/update_cli_reference_help.py --check` | passed |
| `python -m pytest tests/ -v -m "not slow"` | 6,287 passed, 17 skipped, 11 deselected in 13:18 |
| Read-only `TestOutputComparison` | 16 passed |
| `python -m build` and `git diff --check` | passed |

Owner/path evidence: before this follow-up the Web controller always built a
review draft and the Definition painter/collision pass each reapplied record x
to the locked origin. After it, the existing `layout/similarity_alignment.py`
resolver/validator still decides biology; one `app/similarity-alignment.js`
controller branches on its validated status and reuses `runAnalysis` and the
existing Result/History admission. `layout/linear.py` owns one placement rule,
called by both `diagrams/linear/assemble.py` collision planning and
`diagrams/linear/builders.py` painting after final translations. Current
request/Session projection remains in `services/session-request.js`; no new
schema, compatibility reader, Worker path, dependency, or CSP rule was added.
The superseded unconditional-review path and duplicate locked-column offsets
were removed. OE, PE, and CB do not increase; no exception packet applies.
Product preflight for the implemented automatic alignment and Definition
placement remains `IMPLEMENT_EXISTING_AUTHORITY` under the merged revision-3
Choice A and unchanged PD-OI-024. The 390 px acceptance wording follows the
owner's PD-OI-035 revision-2 receipt, with runtime-base authority pending its
separate merge. No runtime is changed under that unmerged authority. Reverting
the Web behavior would require matching Product authority to be restored first;
the independent Definition correction can be managed separately.

Production, tests, public/internal documentation, and generated media were
reviewed separately. The new popup WebP was compared to the old crop at equal
displayed size and passed strict capture verification. The Web reference and
BGC tutorial describe automatic Align, ambiguity, and optional review. The
first full-suite run exposed two generated SVGs from before the Definition
correction and an external browser-test PATH issue. The clean-directory
H-CLI-07 recipe regenerated `docs/images/h-cli-07/cli_losatp_groups.svg`;
`tools/reproduce_examples.py` regenerated the public tutorial-9 arrow figure;
and `tools/refresh_gallery_sessions.py` regenerated the BGC Gallery Session,
source and interactive SVGs, thumbnail, and artifact manifest. For the H-CLI
and tutorial figures and BGC source, structural comparison found only four
Definition transforms changed in each; the BGC Session also updated the
corresponding collision-band coordinates. The owner-wide Gallery refresh
rewrote nine other interactive SVGs solely to place the existing `outline:
none` style on the base pairwise-match selector rather than its hover rule;
no geometry or content changed in those nine. The Gallery manifest and
published Session verification pass. The H-CLI figure, tutorial-9 figure,
and BGC thumbnail were rendered and visually inspected, and the Gallery
Chromium checks passed after refresh. The browser-contract subprocess passed
with the existing Node Playwright binary on PATH. No SVG under
`tests/reference_outputs/` or `examples/gbdraw_social_preview.png` was
changed. The generated browser wheel (SHA-256
`552abef4f9d89319e4aa31f9f4e576acd77311104e16b5e40ba613b58567c306`)
and build outputs remain ignored; unrelated untracked video-plan files in the
original checkout were preserved. No offline bundle audit was triggered:
dependencies, privacy, bundle composition, and Worker lifecycle are unchanged.
