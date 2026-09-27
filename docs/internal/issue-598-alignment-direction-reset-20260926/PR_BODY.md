The Web guide still describes the retired Match checkbox and says Reset always keeps directions. This PR documents the existing exclusive Keep/right/left/Custom choices and both Reset scopes, reproduces the five-BGC tutorial, and strengthens browser coverage of reference-only direction changes. After merge, readers can choose and undo alignment effects using instructions that match the controls already on `dev`.

Refs #598 (BUG-03 and BUG-04). BUG-17 remains outside this change; this PR does not close the whole issue.

## Change class

- [x] STANDARD

## Purpose

Base: `dev@c922fc38aac78da9be83342c09ac0164ecef6ff6`. The branch already contains that base; no additional merge is needed. The runtime is byte-identical to it. Normative PR admission is `Web base policy (trusted base)` and `PR / gate`; neither remote check has been run for this unpublished PR.

## User-visible impact

The documentation explains selected-anchor scope, reference participation, unchanged source strands and unknown/Skip exclusions. It distinguishes the reference's logical pre-Align center from screen pixels and explains both Reset scopes, later manual edits, one-use restoration evidence, fresh Load, Undo/Redo and editable retry.

## Changes

- Update existing Web, Session and Python references, release notes, Tutorial and Gallery instructions.
- Regenerate six existing Tutorial PNGs from five checksum-verified MIBiG sources using the existing `T-GUI-04` capture owner. Preserve labels, legend, metadata and comparison context.
- Restore the reference-only Custom browser journey, including pending form preservation and complete Reset/Undo/Redo; strengthen compact canvas hit checks and successful retry error clearing.
- Preserve S02/S04 evidence and add this review's identity, authority and verification record. Runtime, authority, guards and tracked reference SVGs are unchanged relative to the base.

## Verification

New local review: all five accepted receipts match all nine fields on the base; four required ancestors are retained; runtime and browser-test diff hashes match S04; five source pins and six PNGs pass byte/format checks and visual review. Documentation contracts: **41 passed**. Ruff and whitespace: **PASS**. Web policy: **Gate PASS / Review CLEAR**; the mandatory post-commit result is reported with the candidate SHA in the delivery response. PR wording is checked with `tools/check-pr-language.mjs` before publication.

Reused **S04 results**, not new executions: Python non-slow **6268 passed, 17 skipped, 11 deselected**; Node **712 passed** plus final architecture/Product guards **192 passed**; separate browser **41 passed**; read-only SVG comparison **16 passed**; documentation contracts **41 passed**; real `T-GUI-04` replay with six visually reviewed PNGs and verified SVG/TSV exports. Source, tests, fixtures, capture owner, baseline and recorded environment fields are unchanged. See [review evidence](docs/internal/issue-598-alignment-direction-reset-20260926/evidence/REVIEW_20260927.md) and [S04 evidence](docs/internal/issue-598-alignment-direction-reset-20260926/evidence/S04.md).

Not newly executed: broad Python/Node/browser suites, GUI recipe replay, slow/performance tests, supported-version matrix, remote CI, physical zoom or assistive-technology speech. No wheel is generated and no browser server is reused in this review.

## Risk and review notes

- Gate: PASS, no blocker. Review: CLEAR from the checker, independent of normal maintainer approval.
- **This is not architecture-bearing** relative to the actual base. No owner, canonical path, compatibility reader, schema version, dependency or privileged permission changes; `OE`/`PE`/`CB` scope is unchanged. No `architecture-change` label or exception packet is needed.
- Review the instructions, capture/source evidence and retained browser assertions, especially logical center measurement, pending form, both fresh-Load Reset scopes, zero Reset LOSAT jobs, rollback and ribbon geometry. S04's historical Review REQUIRED is not represented as human approval.
- Prior recorded environment fields match; a complete historic dependency freeze and regenerated wheel-byte comparison are unavailable. Reuse is limited to the identified S04 environment and acceptance, not a new staging or release certification. Future integrated `dev` needs exact-SHA staging.

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
