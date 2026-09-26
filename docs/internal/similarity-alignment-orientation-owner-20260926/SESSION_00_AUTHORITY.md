# Session 00 instruction prompt — publish Product authority for record-owned orientation

Paste this entire file into a new agent session. It is self-contained.

## Mission

In gbdraw's Web app, Linear Similarity Group alignment offers a per-target
**Match reference direction** checkbox. The checkbox fails whenever a
reversal is needed. The cause is that record orientation is stored both in each
record and in the alignment plan. The Product Decision Owner (`satoshikawato`)
selected the following outcome on 2026-09-26:

- orientation is owned only by records;
- the review offers one Match reference direction option for the whole
  alignment;
- a manual Reverse keeps the active alignment;
- Reset Align does not restore orientation.

On the same date the owner approved the five complete `PRODUCT_DECISION`
receipt texts in
`docs/internal/similarity-alignment-orientation-owner-20260926/00_DECISION_PACK.md`
exactly as written:

- `PD-OI-027` revision 4;
- `PD-OI-028` revision 2;
- `PD-OI-029` revision 2;
- `PD-OI-031` revision 4;
- `PD-OI-034` revision 4.

Commit `61eec6c5`, on the authority-only branch
`authority/similarity-alignment-orientation-owner-20260926` (created from
`origin/dev@22cbcca96f2ef5e20bb45fe397ba4232fb158574`), serializes them into
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` as contract revision 16.

This session verifies that serialization and carries it to a merge into
`origin/dev`. It changes no runtime code. Do not ask the owner to approve the
wording again. Do not change the approved wording; a needed correction requires
a new explicit owner response.

## Preparation

1. Read these files:
   - `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`;
   - `docs/internal/PRODUCT_IMPACT_RATCHET.md`;
   - `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` (on the authority
     branch);
   - the Decision Pack;
   - `docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`
     (sections 1–3).
2. Run `git fetch origin`. Locate the authority branch locally or on `origin`.
   - If the branch or commit `61eec6c5` is missing, recreate it from the
     then-latest `origin/dev` with
     `git switch --no-track -c authority/similarity-alignment-orientation-owner-20260926 origin/dev`.
   - Then serialize the five approved texts following the file's existing
     pattern:
     - a `Supersedes` line;
     - `Normative outcome: exactly the approved PRODUCT_DECISION receipt below.`;
     - the decision source;
     - a JSON receipt reproducing the approved fields exactly.
   - Advance the contract revision from `15` to `16` with one metadata line.
3. Check whether `origin/dev` gained newer authority for the same concerns.
   - If `origin/dev` advanced without touching the contract, merge it into the
     authority branch.
   - If a competing decision appeared, stop and report it.

## Verification

- Confirm that the authority-only diff touches only
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`.
- Confirm the five changed records are exactly the ones listed above.
- Confirm every JSON receipt equals the Decision Pack text field for field.
- Confirm all other records are byte-identical to `origin/dev`.
- Confirm `PD-OI-026`, `PD-OI-030`, and `PD-OI-035` stay valid alongside the
  new records.
- Run:
  - `node tests/web/architecture-contracts.test.mjs`;
  - `node --test tests/web/product-impact-ratchet-fixtures.test.mjs`;
  - `node --test tests/ci/ci-impact-policy.test.mjs`;
  - `node tools/check-web-change-budget.mjs --base origin/dev` (expect Gate
    PASS with Review REQUIRED for a governance change);
  - `git diff --check`.

## Publication

Each of these steps needs explicit authorization for that step:

- pushing the authority branch;
- opening an authority-only pull request into `dev`;
- merging it.

In the pull request body, state that the change is authority-only, list the
five records and their revisions, and cite the owner's approval date. Sessions
02–04 start only after the merge into `origin/dev`.

## Handoff

Report:

- the `origin/dev` SHA, the authority branch, and its commit;
- verification results;
- publication status;
- the remaining gate: merge into `origin/dev`, then merge `origin/dev` into
  `fix/similarity-alignment-orientation-owner-20260926`.

Do not describe the gate as complete before the merge.

At the end, print the entire contents of
`docs/internal/similarity-alignment-orientation-owner-20260926/SESSION_01_ERROR_REPORTING.md`
in a copyable code block as the prompt for the next session. Session 01 can
run while the authority merge is pending.
