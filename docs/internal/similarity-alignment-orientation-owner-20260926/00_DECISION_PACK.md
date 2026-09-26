# Product Decision Pack — record-owned orientation for Similarity Group alignment

Status: approved. The Product Decision Owner (`satoshikawato`) explicitly
approved the five receipt texts below exactly as written on 2026-09-26.
Commit `61eec6c5` on the authority-only branch
`authority/similarity-alignment-orientation-owner-20260926` serializes them
as contract revision 16; it becomes runtime authority only after merge into
`origin/dev`. This pack is a
developer-preflight input under the
[Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md). It is not runtime
authority. Active authority is the
[Option Integrity Product Contract](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)
on the runtime base.

## Identity

- Primary concern: `diagram-generation.similarity-alignment.transform-semantics`
  (current `PD-OI-027`, scenario revision `3`), proposed revision `4`.
- Jointly affected concerns:
  - `diagram-generation.similarity-alignment.plan-lifecycle` (`PD-OI-028`,
    revision `1` → `2`);
  - `diagram-generation.similarity-alignment.reset-and-history` (`PD-OI-029`,
    revision `1` → `2`);
  - `diagram-generation.similarity-alignment.surface-scope` (`PD-OI-031`,
    revision `3` → `4`);
  - `web.similarity-alignment.choice-and-retry` (`PD-OI-034`, revision `3` → `4`).
- Concerns checked and left unchanged: `PD-OI-026` (anchor resolution),
  `PD-OI-030` (legacy Session reader), `PD-OI-035` (canvas interaction).
- Discovery lane: developer preflight. The Web Product Impact map has no
  registered Similarity Group alignment concern.
- Prepared from `origin/dev@22cbcca96f2ef5e20bb45fe397ba4232fb158574`
  (fetched 2026-09-26). Implementation branch:
  `fix/similarity-alignment-orientation-owner-20260926`.
- Related closed issues: `#561` (alignment edge cases) and `#586` (alignment UX).

## Trigger

In a Linear diagram, a user selects an exact feature in a Similarity Group,
opens **Review alignment options…**, enables **Match reference direction** for
a target whose anchor lies on the opposite strand, and clicks **Apply**. Apply
always fails with `Alignment generation failed. Review the draft and retry
Apply.` When the target anchor already lies on the same strand, the checkbox
changes nothing. The option therefore never produces a visible effect.

Reproduction on the Gallery example `BGC0000708-BGC0000713`: click
`CAG38712.1` (livA) in the first record, open the review, enable Match
reference direction for *Streptomyces fradiae* ATCC 10745 (anchor neoA on the
opposite strand), and Apply.

The underlying error is
`ValueError: Comparison source feature index conflicts with its view feature ID.`
(`gbdraw/render/groups/linear/pairwise_match.py`). The Web shows only a generic
message because the alignment controller replaces the error published by
`runAnalysis`.

## Root cause

Record orientation currently has two owners:

1. the record's Reverse setting (`linearSeqs[].region_reverse`, projected to
   `presentation.reverseComplement` or `region.reverseComplement` in the render
   request); and
2. the alignment plan (`records[].effectiveReverseComplement`), which Python
   applies late, inside `materialize_similarity_alignment_display()` in
   `gbdraw/api/record_planning.py`.

The browser projects protein-comparison rows into displayed coordinates before
the render request is sent. The projection reads owner 1 through
`buildRegionSpec()` in `gbdraw/web/js/app/run-analysis.js`; each row carries a
`view_feature_svg_id` derived from displayed coordinates. Apply sends the
pre-alignment orientation as owner 1 and puts the reversal only in owner 2.
Python reverses the record after the comparison rows were projected, so the
rows describe the old orientation. The renderer detects the mismatch and
raises. Removing that check would only draw ribbons at wrong positions.

A manual Reverse works because the whole pipeline, including the comparison
projection, reads owner 1.

## Authority search

| Source | Finding | Effect |
| --- | --- | --- |
| Product Impact map / `BD-###` | No registered alignment concern in `tools/web-product-impact-map.json`; no matching record in `tools/web-product-decisions.json`. | Developer preflight. |
| `PD-OI-027` r3 | Requires per-record Match reference direction in review. | Must be superseded. |
| `PD-OI-028` r1 | Requires a manual orientation change to clear the active plan. | Must be superseded. |
| `PD-OI-029` r1 | Requires Reset Align to restore the pre-align orientation. | Must be superseded. |
| `PD-OI-031` r3, `PD-OI-034` r3 | Require per-record orientation choices in the review and draft. | Must be superseded. |
| `PD-OI-026` r3 | Anchor resolution and recommendations are orientation-independent. | Unchanged. |
| `PD-OI-030` r1 | New Sessions preserve the exact resolved plan and effective transform intent. With the new model, anchors persist in the plan and orientation persists in each record. | Unchanged; verify during serialization. |
| `PD-OI-035` r2 | Canvas and 390 px review interaction. | Unchanged. |
| Persisted formats | The typed alignment plan (`layout.similarityAlignment`) and `gbdraw/layout/similarity_alignment.py` are absent from `origin/main` and tag `0.13.0`. Only the legacy `alignOrthogroupFeature` string exists there. | No compatibility reader is required for the plan field change. |

Classification: **PRODUCT_DECISION_REQUIRED**, resolved by the owner's
selection and exact-text approval below.

## Selected outcome

On 2026-09-26 the Product Decision Owner (`satoshikawato`) selected the
following outcome and authorized revising the affected records:

- Orientation is owned only by each record. An alignment plan stores anchors
  and positions only.
- The review offers one **Match reference direction** option for the whole
  alignment, initially off. On Apply it reverses each aligned target whose
  selected anchor has a known displayed strand opposite to the reference
  anchor, by changing that record's Reverse setting in the same way as a manual
  Reverse.
- A manual Reverse after alignment keeps the active plan; the same anchors are
  re-aligned in the new orientation.
- **Reset Align does not restore record orientation.** Undo restores the
  complete previous artifact, including orientation.

On the same date the owner approved the exact receipt wording below without
edits. Serialize only these texts; do not ask for approval again.

## Product Decision Owner response — approved wording

The following five texts were approved exactly as written by `satoshikawato`
on 2026-09-26.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.transform-semantics
Scenario revision: 4
Choice: A / RECORD_OWNED_ORIENTATION_WITH_REVIEW_MATCH
Rationale: Record orientation must have one owner so that features, labels, annotations, and comparison ribbons always follow the same direction. An alignment chooses anchors and positions; matching the reference direction is a one-time change to record orientation made through the same path as a manual Reverse.
Must preserve: Each record's current orientation by default; one explicit Match reference direction option in the opened review, initially off; reversal on Apply only of aligned targets whose selected anchor has a known displayed strand opposite to the reference anchor; unchanged orientation for targets with an unknown strand and for skipped, missing, and unusable targets, with each target's resulting direction shown before Apply; the reference record's position and orientation; every target's vertical position; exact idempotent anchor-center alignment in the resulting orientation; comparison ribbons, features, labels, and annotations drawn in their record's orientation; readable text; an accurate source-relative rev indication; and one atomic validated commitment of orientation and alignment.
May retire: Per-record Match reference direction controls in the review; storage of an orientation policy or effective orientation inside the alignment plan; reversal of records during rendering from alignment-plan data.
Accepted residual risk: Match reference direction applies to every eligible target in the review. A user who wants a different direction for one record changes that record's Reverse setting after Apply; the alignment is kept and recalculated.
Owner: satoshikawato
Decision date: 2026-09-26
```

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.plan-lifecycle
Scenario revision: 2
Choice: A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN
Rationale: An active plan identifies anchors, not directions. A manual orientation change leaves every anchor valid, so keeping the plan spares the user a second alignment after reversing a record.
Must preserve: Survival of the active plan across regeneration after style, label, and canvas-size changes, and across record reorder when stable record identities remain; survival across manual record orientation changes, with the same anchors aligned in the new orientation when the diagram is next generated; explicit clearing and notification after manual record movement, source replacement, crop changes, or record-selector changes; validation before regeneration; the last successful Result while a stale plan is repaired; and explicit reselect, Skip, or Clear actions without automatic anchor substitution.
May retire: Clearing the active plan when the user manually changes a record's orientation.
Accepted residual risk: After a manual orientation change, the reversed record moves horizontally on the next generation so that its anchor stays aligned. A user who wanted the previous position uses Undo or Reset Align.
Owner: satoshikawato
Decision date: 2026-09-26
```

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.reset-and-history
Scenario revision: 2
Choice: A / IMMEDIATE_PREALIGN_POSITION_BASELINE
Rationale: Reset Align removes the alignment's positioning. Record orientation is ordinary record state that Reset Align does not change; normal Undo restores the complete previous artifact, including orientation.
Must preserve: Replacement of the preceding active plan by a new Align; restoration by Reset Align of the record positions immediately before that Align, together with clearing of the active plan; one atomic history transaction for Apply, Reset, and manual clearing; normal Undo of the complete prior artifact, including record orientation and any prior active plan; and no committed history entry after failed, canceled, superseded, or stale work.
May retire: Restoration of record orientation by Reset Align, including orientation changed by Match reference direction.
Accepted residual risk: After an Apply with Match reference direction, Reset Align leaves the reversed records reversed. Returning to the previous orientation requires Undo or the record's Reverse setting.
Owner: satoshikawato
Decision date: 2026-09-26
```

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.surface-scope
Scenario revision: 4
Choice: A / AUTO_APPLY_RESOLVED_WITH_EXPLICIT_REVIEW_SINGLE_MATCH
Rationale: A fully resolved alignment should complete without an unnecessary review, and an opened review should offer direction matching as one clear option instead of repeating it for every target.
Must preserve: The exact selected reference from the feature popup and Similarity Groups drawer; automatic application of a Python-resolved default alignment; a full, accessible review with Select, Skip, candidate details, a resolution summary, one Match reference direction option, and each target's resulting direction when Python reports ambiguity or the user explicitly requests review; typed fully resolved plans, shared Python validation, strict CLI ambiguity rejection, actionable errors that state the underlying failure, and accurate unsupported-feature disclosure.
May retire: Per-record orientation controls in the review.
Accepted residual risk: The review cannot match the direction of only some targets. Per-record exceptions use the record's Reverse setting after Apply.
Owner: satoshikawato
Decision date: 2026-09-26
```

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.choice-and-retry
Scenario revision: 4
Choice: A / AUTO_APPLY_RESOLVED_WITH_REVIEW_RETRY_SINGLE_MATCH
Rationale: The Web should keep automatic application of resolved plans and local review editing, with one draft-level direction option replacing per-target orientation choices.
Must preserve: Independent Select or Skip choices per target and one draft-level Match reference direction option; Python ownership of candidate eligibility, strand facts, and final validation; local no-Worker draft editing in review; one batch validation for an edited draft; correction and retry without losing the draft, with the underlying failure message shown; visible candidate facts and each target's resulting direction before Apply; the last Result, record orientations, and History after failure, Cancel, stale, or superseded work; canvas interaction during review; exposure of the retained draft and retry after an automatic render failure without a partial commit; and plan regeneration, Session, Reset, and Undo/Redo meanings as defined by the current plan-lifecycle and reset-and-history decisions.
May retire: Independent per-target orientation choices in the draft.
Accepted residual risk: The draft cannot match the direction of only some targets. Per-record exceptions use the record's Reverse setting after Apply.
Owner: satoshikawato
Decision date: 2026-09-26
```

## Engineering note

This outcome removes a duplicate owner instead of adding a synchronization
rule. The alternatives are rejected:

- Re-projecting comparison rows in Python after a plan-driven reversal would
  add a second projection path.
- Removing the renderer's ID check would hide the defect and draw misplaced
  ribbons.
