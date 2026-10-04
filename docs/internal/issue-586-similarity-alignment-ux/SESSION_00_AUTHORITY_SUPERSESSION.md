# Session 00 instruction prompt: serialize the approved Product outcomes

## Mission

Serialize the four signed Issue #586 Product decisions as revision-2
replacements in gbdraw's durable Product Contract. This is an authority-only
session. Do not change runtime code, UI, tests, schemas, or generated assets.

The runtime change is not authorized by a plan document or by an unmerged
candidate authority branch. Dependent implementation begins only after this
authority change has been reviewed and merged into `origin/dev`.

## Branch and publication boundary

Create and use this authority-only branch from the latest `origin/dev`:

```text
issue-586-similarity-alignment-authority-20260924
```

Do not perform this session on the runtime branch. Before creating the branch,
fetch `origin`, inspect local changes, and preserve unrelated work. The expected
Product Contract revision at plan creation is 12; increment the actual current
revision by one without overwriting newer unrelated authority.

When this prompt is supplied as the session request, it authorizes one commit
and a push to the same-named remote authority branch. It does not authorize a
PR, merge, deployment, tag, release, or direct push to `dev`/`main`.

## Read before editing

Read all of these completely:

1. `AGENTS.md` and `CLAUDE.md`;
2. `docs/internal/PRODUCT_IMPACT_RATCHET.md`;
3. `docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`;
4. `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`;
5. `tools/web-product-decisions.json`;
6. section 3 of
   `docs/internal/issue-586-similarity-alignment-ux/IMPLEMENTATION_PLAN.md` on
   branch `issue-586-similarity-alignment-ux-20260924`.

The master plan is intentionally stored on the implementation branch and will
not exist in an authority branch created directly from `origin/dev`. Read it
before switching branches or inspect it without copying it into the authority
branch:

```bash
git show issue-586-similarity-alignment-ux-20260924:docs/internal/issue-586-similarity-alignment-ux/IMPLEMENTATION_PLAN.md
```

The master plan contains the full, signed `PRODUCT_DECISION` receipts. Copy
their content exactly, normalizing only Markdown/HTML transport escaping. Do
not translate, shorten, reinterpret, or add inferred obligations.

## Required authority changes

Replace these existing accepted records with the same `PD-OI` identifiers and
scenario revision 2:

| Record | Concern | New selected outcome |
| --- | --- | --- |
| `PD-OI-026` | `diagram-generation.similarity-alignment.anchor-resolution` | `A / AUTO_PRESELECTED_SUGGESTIONS` |
| `PD-OI-027` | `diagram-generation.similarity-alignment.transform-semantics` | `A / SINGLE_ALIGN_PER_RECORD_ORIENTATION` |
| `PD-OI-031` | `diagram-generation.similarity-alignment.surface-scope` | `A / SINGLE_REVIEW_ALIGNMENT_SURFACE` |
| `PD-OI-034` | `web.similarity-alignment.choice-and-retry` | `A / PRESELECTED_LOCAL_BATCH_REVIEW` |

For each record:

- set scenario revision to 2 and status to `ACCEPTED`;
- identify the prior scenario-revision-1 outcome in a `Supersedes` statement;
- write a concise normative outcome that is no broader than the receipt;
- reproduce the complete signed receipt in the existing JSON evidence format;
- identify `satoshikawato`, `2026-09-24`, and Issue #586 as the decision source;
- state that the serialization adds no terms and becomes runtime authority only
  after merge;
- keep unaffected decisions, including plan lifecycle, reset/history, Session
  compatibility, and canvas interaction, unchanged.

Update the contract revision/changelog using the document's existing format.
Do not create a second decision registry or a `BD-###` entry. Update
`tools/web-product-decisions.json` only if its documented schema requires a
mechanical revision/reference change; do not invent a candidate decision there.

## Required consistency review

Check the resulting active contract as an AND-of-OR set, not merely by matching
IDs. Confirm that the four replacements coexist with:

- the active alignment plan lifecycle;
- immediate-pre-align Reset semantics;
- reader-only released Session compatibility;
- floating palette and canvas-pick behavior;
- generic request/artifact, accessibility, failure-isolation, and architecture
  acceptance contracts.

If another already-merged decision conflicts with any signed receipt, stop and
report the exact conflict. Do not resolve a new Product choice autonomously.

## Verification

Run the Product Contract/decision checks named by the current policy and at
minimum:

```bash
node tests/web/architecture-contracts.test.mjs
git diff --check
```

Review the Product Contract diff independently and verify that no runtime or
test file changed.

## Commit and handoff

Use an English commit title such as:

```text
Supersede similarity alignment UX decisions
```

Push only `issue-586-similarity-alignment-authority-20260924`. Report the branch,
commit, exact changed files, verification results, and the fact that runtime
work remains blocked until the authority commit is merged into `origin/dev`.
Do not mark Session 00 complete in the implementation-plan ledger until the
merge commit is known.
