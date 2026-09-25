# Session 00 instruction prompt — Product authority for record-owned orientation

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

The owner authorized revising the affected Product records. This session turns
that decision into merged durable authority. It changes no runtime code.

Read first:

- `docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`
  (sections 1–3);
- `docs/internal/similarity-alignment-orientation-owner-20260926/00_DECISION_PACK.md`.

The Decision Pack contains five **unsigned** proposed `PRODUCT_DECISION`
receipts:

- `PD-OI-027` revision 4;
- `PD-OI-028` revision 2;
- `PD-OI-029` revision 2;
- `PD-OI-031` revision 4;
- `PD-OI-034` revision 4.

The owner's selection of the outcome is not approval of the exact wording.

## Preparation

1. Read these files:
   - `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`;
   - `docs/internal/PRODUCT_IMPACT_RATCHET.md`,
     `docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`,
     `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`,
     `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`;
   - the Decision Pack and the master plan.
2. Run `git fetch origin`. Record the `origin/dev` SHA and check for newer
   authority or competing decisions on the same concerns. Preserve unrelated
   working-tree content.
3. Compare each proposed receipt with its current record and with the jointly
   required effects of the unchanged records. Check:
   - exact reference, orientation, and translation;
   - Reset, Undo, and Session behavior;
   - retry and error behavior.
4. Confirm that `PD-OI-026`, `PD-OI-030`, and `PD-OI-035` stay valid without
   revision. If evidence shows a proposal is inaccurate or incomplete, correct
   the Decision Pack on `fix/similarity-alignment-orientation-owner-20260926`
   in a separate commit. Do not call a correction a signed decision.

## Decision boundary

Present the five complete texts to the Product Decision Owner. Ask for explicit
approval of the exact wording or explicit edits, including owner and decision
date. Do not infer or fill any of the following:

- rationale;
- preservation or retirement scope;
- accepted residual risk;
- signer or date.

Until a complete explicit response arrives, serialize nothing and do not start
Sessions 02–04. Session 01 may proceed.

## Authority work, after explicit approval

1. Fetch `origin/dev` and run
   `git switch --no-track -c authority/similarity-alignment-orientation-owner-20260926 origin/dev`.
2. In `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`, replace the five
   records with the approved revisions. Follow the file's existing pattern:
   - a `Supersedes` line naming the prior revision and choice;
   - the normative outcome;
   - the decision source;
   - the JSON receipt reproducing exactly the approved fields.
3. Advance the contract revision from `15` to `16` and add one metadata line
   describing the change. Do not accumulate superseded records. Leave all other
   `PD-OI` records unchanged.
4. Do not add a `BD-###` record or a second decision store, and do not include
   runtime or test changes. Apply any mechanical registry edit the current
   policy requires.
5. Run `node tests/web/architecture-contracts.test.mjs`, `git diff --check`,
   and any current Product Contract checker. Review the diff against the
   approved wording line by line.
6. Commit only on the authority branch.
7. Push the authority branch, open a pull request, or merge only with explicit
   authorization for each of those steps. Sessions 02–04 require the merge into
   `origin/dev`.

## Handoff

Report:

- the `origin/dev` SHA and the authority branch;
- the decision status and exact changed files;
- verification results and publication status;
- the remaining gate: merge into `origin/dev`, then merge `origin/dev` into the
  fix branch.

Do not describe the gate as complete because an authority commit exists.

At the end, print the entire contents of
`docs/internal/similarity-alignment-orientation-owner-20260926/SESSION_01_ERROR_REPORTING.md`
in a copyable code block as the prompt for the next session. Session 01 can
run while the authority merge is pending.
