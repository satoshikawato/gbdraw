# Session 00 instruction prompt — authorize Choice A

Paste this entire file into a new Codex session. It is self-contained.

## Mission and exact context

gbdraw Linear Similarity Group alignment currently opens Select alignment anchors for every valid alignment. The requester chose Choice A on 2026-09-25: normal Align should directly apply a fully resolved Python plan; ambiguity still opens the full review; a visible Review alignment options… action in both feature popup and Similarity Groups drawer must open full review on demand. Lock Definition Column has a separate regression to be fixed on the implementation branch.

The fixed implementation branch is fix/similarity-alignment-auto-apply-definition-column-20260925. This session handles Product authority only and must use a separate authority-only branch from the then-latest origin/dev. Do not change dependent runtime, tests, or this plan in the authority branch. The requester selected the outcome, but the four full candidate PRODUCT_DECISION responses in docs/internal/issue-586-alignment-followup-20260925/00_DECISION_PACK.md have not been explicitly approved as exact text. A choice code alone is not the complete receipt required by docs/internal/PRODUCT_IMPACT_RATCHET.md.

## Preparation

1. Read AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md, docs/internal/PRODUCT_IMPACT_RATCHET.md, docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md, docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md, docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md, the Decision Pack, and 01_MASTER_PLAN.md.
2. Inspect git status, current branch/upstream, origin/dev SHA, and the accepted base records PD-OI-026, PD-OI-027, PD-OI-031, and PD-OI-034. Fetch origin first. Preserve unrelated changes and check for newer Product authority or competing decisions.
3. Compare all four proposed Choice A receipts with the accepted records and their jointly required effects, including exact reference, candidate eligibility, orientation, typed validation, review, retry, Result/History, and Session contracts. Confirm the Decision Pack still names the complete product outcome. Correct a draft in the plan branch if evidence makes it inaccurate; do not call a correction a signed decision.

## Decision boundary and authority work

Present the four complete Choice A PRODUCT_DECISION texts to the Product Decision Owner for explicit approval of their exact wording or explicit edits. The owner's 2026-09-25 selection of A authorizes planning the route but is not permission to infer rationale, retirement intent, accepted risk, signer, or decision date. Until complete text is approved, do not serialize it as accepted authority or implement the dependent Web runtime. Independent Session 01 may proceed.

After an explicit complete response arrives, fetch origin/dev and create an authority-only work branch from it with git switch --no-track -c. Serialize only the approved four revised concern records in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md, updating its current revision and supersession notes. Keep unrelated PD-OI records unchanged. Do not add a second decision store, invent a BD record, or include runtime/test changes. Check the current policy for any required mechanical registry edit.

Run node tests/web/architecture-contracts.test.mjs and git diff --check, plus any current Product Contract checker. Review the authority-only diff and the signed wording independently. Commit only on that branch. Push, open a PR, or merge only with authorization for that specific publication step. Dependent Sessions 02–03 start only when this authority is merged into origin/dev and then brought into the fixed implementation branch.

## Handoff

Report base SHA, branch, decision status, exact changed files, verification, publication status, and the remaining authority gate. Do not mark the gate complete merely because an authority candidate exists. At the end of this session, print the entire contents of SESSION_01_DEFINITION_COLUMN.md in a copyable code block as the exact prompt for the user's next Codex session; Session 01 can run while the authority merge is pending.
