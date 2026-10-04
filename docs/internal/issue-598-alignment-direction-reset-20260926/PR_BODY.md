The Web guide still describes the retired Match checkbox and says Reset always keeps directions. This PR documents the existing exclusive Keep/right/left/Custom choices and both Reset scopes, reproduces the five-BGC tutorial, and strengthens browser coverage of reference-only direction changes. After merge, readers can choose and undo alignment effects using instructions that match the controls already on `dev`.

Refs #598 (BUG-03 and BUG-04). BUG-17 remains outside this change; this PR does not close the whole issue.

## Change class

- [x] STANDARD

## Purpose

Base: `dev@f5f86634459e0dcd46c1a452e9219fbba635d429`. PR #617 advanced dev while this PR was checking CI, so the branch normally merges that admitted base and preserves every earlier writer. The runtime is byte-identical to the new base. Normative PR admission is `Web base policy (trusted base)` and `PR / gate`. Both succeeded on the initial head; checks for the synchronized head are reported separately.

## User-visible impact

The documentation explains selected-anchor scope, reference participation, unchanged source strands and unknown/Skip exclusions. It distinguishes the reference's logical pre-Align center from screen pixels and explains both Reset scopes, later manual edits, one-use restoration evidence, fresh Load, Undo/Redo and editable retry.

## Changes

- Update existing Web, Session and Python references, release notes, Tutorial and Gallery instructions.
- Regenerate six existing Tutorial PNGs from five checksum-verified MIBiG sources using the existing `T-GUI-04` capture owner. Preserve labels, legend, metadata and comparison context.
- Restore the reference-only Custom browser journey, including pending form preservation and complete Reset/Undo/Redo; strengthen compact canvas hit checks and successful retry error clearing.
- Preserve S02/S04 evidence and add this review's identity, authority and verification record. Runtime, authority, guards and tracked reference SVGs are unchanged relative to the base.

## Verification

Initial review: five accepted receipts match all nine fields; four required ancestors, five source pins and S04's browser-test diff digest are retained. Documentation contracts: **41 passed**; Ruff, original six-image visual review, language and final-commit Web policy passed on the original head.

After the admitted dev runtime advance, newly executed checks cover the shared History/Session/render boundary: **118 focused Node tests passed**; native alignment **83**, documentation **41** and read-only SVG **16** tests passed (**140 total**); **16 browser tests passed** for alignment, fresh Load, both Reset scopes, zero Reset LOSAT jobs, rollback/pending form, History, palette and region. Ruff passes. Reference-center delta remains below 0.000018 logical px, with 77 ribbons and maximum anchor offset 0.0001221 px. The real `T-GUI-04` recipe validates SVG (249702 bytes), TSV (232 rows) and both Reset scopes. Current image/check details are recorded in the integration evidence.

Historical **S04 results**, not new runs or blanket evidence for the changed runtime: Python non-slow **6268 passed, 17 skipped, 11 deselected**; Node **712 passed** plus final guards **192 passed**; separate browser **41 passed**; documentation **41 passed**; read-only SVG **16 passed** and original GUI replay. The admitted #600 writer also records its own **6642-test** native verification; this is prior base evidence, not this PR's new execution. See [dev synchronization evidence](docs/internal/issue-598-alignment-direction-reset-20260926/evidence/PR_DEV_SYNC_20260927.md), [original review](docs/internal/issue-598-alignment-direction-reset-20260926/evidence/REVIEW_20260927.md) and [S04 evidence](docs/internal/issue-598-alignment-direction-reset-20260926/evidence/S04.md).

Not newly executed after synchronization: whole Python/Node/browser suites, slow/performance tests, supported-version matrix, physical zoom or assistive-technology speech. Browser verification uses a clone-local rebuilt wheel and an isolated server port; no earlier session's server is reused. Post-sync remote CI is reported for the new head, independently from the initial green run.

## Risk and review notes

- Gate: PASS, no blocker. Review: CLEAR from the checker, independent of normal maintainer approval.
- **This is not architecture-bearing** relative to the actual base. No owner, canonical path, compatibility reader, schema version, dependency or privileged permission changes; `OE`/`PE`/`CB` scope is unchanged. No `architecture-change` label or exception packet is needed.
- Review the instructions, capture/source evidence and retained browser assertions, especially logical center measurement, pending form, both fresh-Load Reset scopes, zero Reset LOSAT jobs, rollback and ribbon geometry. S04's historical Review REQUIRED is not represented as human approval.
- Recorded environment fields match. A complete historic dependency freeze is unavailable; old broad results remain historical, and fresh scoped tests cover the changed shared boundary. Future integrated `dev` needs exact-SHA staging.

## Rollback

Remove this PR's documentation, media and test changes through a normal follow-up PR, retaining the existing `dev` runtime and authority. The inherited S04 commit has two parents; do not blindly revert its merge or earlier writers' work.

## Conditional Product Impact

- Product-impact role: EVIDENCE_ONLY.
- Product preflight: IMPLEMENT_EXISTING_AUTHORITY; no registered subject delta or unresolved material outcome.
- Authority: PD-OI-027 r5, PD-OI-029 r3, PD-OI-031 r5, PD-OI-034 r5 and PD-OI-039 r2 match Decision Packs 01–05. `EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH` is accepted; PD-OI-035 r3 and independent canvas/keyboard/focus/Editor requirements remain jointly necessary.
- Behavior contracts and A01–A16 evidence: the linked review and S04 reports. No new decision or retirement is requested.

<!-- gbdraw-product-impact-decision:start -->
{"schemaVersion":1,"headSha":"","decisions":[]}
<!-- gbdraw-product-impact-decision:end -->
