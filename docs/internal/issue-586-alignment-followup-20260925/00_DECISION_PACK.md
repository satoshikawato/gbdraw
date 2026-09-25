# Product Decision Pack — unambiguous Similarity Group alignment

Status: Choice A selected by the requester on 2026-09-25; the complete proposed response texts remain unsigned pending explicit Product Decision Owner approval of their wording.

This pack is a developer-preflight input under the [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md). It is not runtime authority. Choice A is selected for implementation planning; its four complete candidate response texts remain proposals, not signed decisions.

## Identity

- Primary concern: `diagram-generation.similarity-alignment.surface-scope`, proposed scenario revision `3` (current `PD-OI-031` revision `2`).
- Jointly affected concerns for Choice A: `diagram-generation.similarity-alignment.anchor-resolution` (current `PD-OI-026` revision `2`), `diagram-generation.similarity-alignment.transform-semantics` (current `PD-OI-027` revision `2`), and `web.similarity-alignment.choice-and-retry` (current `PD-OI-034` revision `2`). Each would advance to revision `3`; the Python anchor-selection algorithm remains unchanged. Choice B requires supersession of `PD-OI-031` and `PD-OI-034`; it preserves the existing anchor and orientation contracts.
- Discovery lane: developer preflight. The Web Product Impact map has no registered Similarity Group alignment concern.
- Prepared from local `origin/dev@ac067f9b5482b480a1edca983709e3321f75157f`; planning and implementation branch `fix/similarity-alignment-auto-apply-definition-column-20260925` (created at that base).
- Prepared by: Codex, 2026-09-25.
- Related work: issue `#586`; user report of a five-record Linear alignment and a separate Lock Definition Column regression.

## Trigger and classification

The user asks that `Select alignment anchors` stop displaying candidate choices on every Align when no target is ambiguous. The existing workflow opens the full palette even when the Python resolver returns a fully resolved plan. Skipping or shortening that review changes when alignment is committed and whether the user can choose `Match reference direction` before commitment.

Authority search:

| Source | Finding | Effect on this proposal |
| --- | --- | --- |
| Product Impact concern/map | No registered Similarity Group alignment concern in `tools/web-product-impact-map.json`. | Use developer preflight; do not infer an automatic Gate result. |
| Active durable `BD-###` | No matching record in `tools/web-product-decisions.json`. | No separate mapped decision selects the proposed workflow. |
| Static Product Contract | `PD-OI-031` requires a review surface for every Align. `PD-OI-034` requires an editable draft and Apply. `PD-OI-027` provides per-record orientation in review and defaults to preservation. `PD-OI-026` distinguishes automatic resolution from recommendations for ambiguity and promises replacement or Skip before final validation. | Current revision-2 authority selects the current full review. Choice A must supersede all four concerned records; Choice B must supersede the two review-surface records. |
| Domain and integrity rules | The exact reference, independent target resolution, Python eligibility and validation, unchanged missing/unusable/skipped records, and idempotent transforms remain required. | Both choices below must preserve these effects. |
| Released compatibility evidence | No request or Session schema change is proposed. Existing plans, saved Results, and reader compatibility remain subject to their accepted contracts. | Do not rewrite or reinterpret stored plans. |
| Eligible exact-head current decision | None supplied for this proposed head. A local PR decision cannot retire the mandatory review affordance. | Durable authority is required. |
| Current code/tests | `resolveRequest()` always makes a draft and enters `reviewing`. Python's `resolved` response already includes a validated plan. The one-usable-member browser test currently expects a palette and manual Apply. | Evidence of the present path, not authority for continuing it. |

Result: **UNRESOLVED AUTHORITY** for the selected new workflow; there is no conflict among current active authorities. Classification: **PRODUCT_DECISION_REQUIRED**. `IMPLEMENT_EXISTING_AUTHORITY` cannot select a shorter workflow against `PD-OI-031` and `PD-OI-034`; measurements alone cannot choose between the two valid UX outcomes; neither choice below violates a non-waivable rule if its stated safeguards hold. Dependent runtime work waits for explicit Product Decision Owner disposition and merged authority.

## User journey and checkpoints

- Actor/context/goal: a user viewing a Linear Similarity Group diagram selects one exact reference feature and wants the other displayed records aligned with minimal repetitive interaction.
- Entry points: feature popup and Similarity Groups drawer `Align…` action. The clicked or explicitly selected record/feature remains the reference.
- Preconditions: a current diagram and a Python resolver response tied to its exact artifact, group, records, and source features.
- Current major steps: invoke Align, wait for resolution, inspect every target in a full palette, optionally change anchor/orientation or Skip, then Apply.
- Affected checkpoints: immediate busy feedback; resolved-versus-ambiguous branching; pre-commit orientation access; validated generation; Result and History replacement; failure/retry; Session load, regenerate, Reset, Undo/Redo, and export.
- Current observed behavior: the full palette opens for `resolved` and `ambiguous` responses alike. `0 need a choice` is not an ambiguity indicator because recommended ambiguous choices are also preselected.
- Next action after Choice A: inspect the committed summary, Undo or Reset, or use the explicit review action to choose anchors/orientation before a new commit. On failure, correct the retained draft and retry.
- Next action after Choice B: inspect the compact resolved summary, change orientation or Skip, then Apply or Cancel. On failure, correct the retained draft and retry.

## Non-waivable constraints

- Architecture: retain one JS alignment controller, one Python resolver/validator, and the typed render-request path; no second alignment planner or render contract.
- Security/privacy: genome inputs and results stay local; preview SVG still passes the shared sanitizer; no new runtime network dependency.
- Scientific correctness: branch on Python `response.status`, never on candidate count, recommendation count, visible badges, or `unresolvedCount`; `resolved` must carry a validated explicit plan. Unique representative and deterministic candidate 1 remain disclosed suggestions in `ambiguous`, not automatic scientific resolutions.
- Orientation: preserve each target's current direction by default; reverse only on explicit per-target `Match reference direction` with known opposite displayed strands. Preserve exact reference and vertical positions.
- Persistence and recovery: one successful user action updates canonical Result/History; error, stale result, cancellation, or supersession preserves the previous Result/History. Session, regeneration, Reset, Undo/Redo, and export retain their accepted meanings.
- Accessibility: both paths remain reachable by keyboard and narrow viewport, with status and errors announced. An optional review path cannot depend on a hidden mouse gesture.
- Evidence: deterministic resolver fixtures, focused controller/browser tests, and an actual SVG/Chromium check for the Definition regression described below.

## Choice A — AUTO_APPLY_RESOLVED_WITH_EXPLICIT_REVIEW

- Complete outcome: the normal `Align…` action applies a Python `resolved` plan without opening `Select alignment anchors`. An `ambiguous` response opens the existing full review with all target rows, recommendations, Select, Skip, and orientation choices. A clearly named, keyboard-accessible `Review alignment options…` action from the same feature popup and drawer always opens the full review before commitment, including for a resolved plan. Both actions use the same exact reference and resolver. `Match reference direction` defaults off. If automatic generation fails, the retained resolved draft opens for correction and retry; it does not commit a partial plan.
- Preserved: exact reference; automatic only-usable and unique-direct-RBH decisions; independent record outcomes; explicit review access to all candidate and orientation choices; Python validation; summary; existing persistence and recovery contracts.
- Added: one-action application on the common unambiguous path; an explicit way to request full review before committing that path.
- Lost/retired: mandatory review and Apply for `resolved` responses; the single-action-only entry surface. No review choice is retired from the explicit review path.
- Discoverability/accessibility: show distinct accessible action names and busy states together in the popup/drawer; neither route may be available only through hover, modifier keys, or the canvas.
- Canonical state and Undo/Redo: success commits one validated plan and one Result/History action; Undo and Reset return to the pre-align state.
- Session/regeneration/export: existing plan and Result serialization and regenerated/downloaded diagram semantics remain unchanged; transient draft and canvas markers are not saved or exported.
- Validation/error/recovery: the initial `resolved` response supplies the validated plan. Render still uses typed validation. A failed automatic render exposes the preserved draft and actionable error; stale work cannot replace the current Result. The explicit review path validates the edited batch at Apply.
- Scientific output/cache/provenance: default resolved anchor centers and preserved orientation match today's immediately applied draft; no new inference or cache key. Summary records the resolved outcome without presenting a heuristic as biological certainty.
- Performance: the resolved default path avoids a second helper call and human Apply wait; render cost is unchanged. Explicit review retains current costs.
- Compatibility/architecture: no schema change; one controller branches on resolver status and reuses its existing apply/review paths. Rendering remains Python-owned.
- Evidence available/missing: the response schema and resolved-plan validation already exist. Browser evidence for both entry actions, automatic failure/retry, and committed-state parity is missing.
- Residual risk: a user may commit the default alignment before noticing that they wanted to change direction or Skip; explicit review must be conspicuous, and Undo/Reset must remain reliable.
- Route: **DURABLE_AUTHORITY_REQUIRED**. Next action if selected: obtain complete owner responses for all four affected concerns, merge an authority-only change into the runtime base, then implement and test.

## Choice B — COMPACT_CONFIRMATION_FOR_RESOLVED

- Complete outcome: `Align…` opens a compact confirmation when Python returns `resolved`. It shows the exact reference, per-target resolved outcome and effective direction, `Match reference direction` and Skip, plus Apply/Cancel, without expanding candidate lists. An `ambiguous` response opens the existing full candidate review. A resolved target can expand its candidate detail on request before Apply; this explicit expansion uses the same resolver facts and local draft.
- Preserved: one `Align…` entry, a review and Apply before every transform, exact reference, all pre-commit orientation and Skip choices, Python validation, summary, persistence, and failure recovery.
- Added: a shorter default display for already resolved targets, with optional candidate detail.
- Lost/retired: default display of all candidate cards in a resolved review. No confirmation step is retired.
- Discoverability/accessibility: the compact summary exposes outcome, direction, and optional detail by keyboard and on narrow viewports.
- Canonical state and Undo/Redo: nothing changes until Apply; one successful Apply commits one Result/History action. Cancel, Undo, and Reset keep their meanings.
- Session/regeneration/export: saved plan, regenerated diagram, and downloads are unchanged; compact/expanded display is transient.
- Validation/error/recovery: edits stay local until one batch validation at Apply. A retryable failure keeps the draft and compact/expanded state; stale work preserves the last Result.
- Scientific output/cache/provenance: accepted defaults have the same anchors and orientation as the current palette; optional details disclose evidence without promoting suggestions to facts. No new cache policy.
- Performance: the Python helper/render sequence stays as today; fewer visible cards may reduce layout work, but Apply still requires a second helper call and user confirmation.
- Compatibility/architecture: no schema change; the same controller and typed resolver own choices and validation.
- Evidence available/missing: the current draft and validation path exist. Browser evidence for compact display, keyboard expansion, and candidate/effective-direction parity is missing.
- Residual risk: a user may overlook a resolved candidate's details until expanding it, and the extra confirmation remains even when no choice is needed.
- Route: **DURABLE_AUTHORITY_REQUIRED**. Next action if selected: obtain complete owner responses for the surface and Web review concerns, merge authority, then implement and test.

## Comparison

| Dimension | A: automatic resolved path | B: compact resolved confirmation |
| --- | --- | --- |
| Normal resolved action | Commits after Python resolve and render | Waits for Apply |
| Ambiguous action | Full editable palette | Full editable palette |
| Orientation or Skip before commit | Explicit `Review alignment options…` action | In compact confirmation |
| Feedback and accessibility | Immediate busy state, then committed summary; named second action | Immediate busy state, then keyboard-accessible compact summary |
| Canonical Result/History | One successful replacement | One successful replacement after Apply |
| Error and recovery | Failed automatic render opens retained review draft | Failed Apply retains review draft |
| Session/regeneration/export | Existing plan and diagram contract | Existing plan and diagram contract |
| Validation | Validated initial resolved plan, typed render | Batch Python validation at Apply, typed render |
| Performance | One fewer helper call and confirmation on resolved default | Current helper and confirmation cost |
| Main residual risk | Intent to orient or Skip may be discovered after default commit | Repetitive confirmation remains |
| Authority route | Durable supersession | Durable supersession |

## Evidence-first option

Deterministic evidence can verify resolver status, output parity, timing, and the recovery path, but cannot decide whether immediate commitment or an orientation confirmation is preferable. No evidence-only PR or additional runtime/default change is required before the Product Decision Owner chooses between A and B. The implementation PR must provide the missing contracts listed above.

## Independent Lock Definition Column regression

Classification: **IMPLEMENT_EXISTING_AUTHORITY** under `PD-OI-024` (`linear.definition-display`). The lock contract already requires a common left edge. `place_linear_definition(..., keep_left=True)` computes a shared origin, but Linear assembly and drawing add each record's final horizontal translation to that origin. Alignment makes those translations differ, so each Definition block is internally left-aligned but the blocks drift apart. Collision bands use the same shifted origin, hiding the error from bounds consistency checks. Existing browser coverage generates a diagram without different per-record translations.

The repair should place locked row/side Definition groups from one final shared column using post-alignment record positions, preserve the configured gap from the leftmost sequence, and use that same placement in drawing and collision planning. Record-local labels above their sequences continue to follow those sequences; Lock-off definitions continue to follow their rows. Add a browser regression with nonzero, unequal positive and negative translations, checking actual SVG left edges and gap after Align, Generate, save/load, and regeneration. Include single and mixed multi-record rows. This correction can proceed independently after the normal architecture and test gates; it requires no new Product choice.

## Engineering recommendation

The requester selected **Choice A**. It was recommended because it directly removes the repeated candidate/Apply step reported by the user while retaining a visible, explicit route to orientation and Skip before commitment. This is engineering advice, **not Product authority**. Choice B is a valid, smaller change in commitment timing if the owner values mandatory pre-commit orientation review more than the one-action path.

## Product Decision Owner response — proposed wording for Choice A

The following four complete texts are **unsigned proposals**. Approval of the idea alone is not a signed receipt. The Product Decision Owner may approve these exact texts or supply corrections; Codex must serialize only the explicit response. `Owner` and `Decision date` below identify the proposed signer and preparation date, not an assertion that a decision occurred.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.anchor-resolution
Scenario revision: 3
Choice: A / DETERMINISTIC_RESOLUTION_WITH_OPTIONAL_REVIEW
Rationale: Python should keep its deterministic automatic resolutions and disclosed ambiguity recommendations, while a fully resolved default alignment can commit without opening candidate review unless the user requests it.
Must preserve: The exact selected reference; stable record and biological-feature identity; only-usable-candidate and unique-direct-RBH automatic resolution; independent displayed-record treatment; unchanged missing, unusable, and skipped records; visible recommendation reasons and local replacement or Skip in every opened review; final shared Python validation; and independence from viewport, scroll, ribbon geometry, confidence score, supporting-edge count, and multi-hop evidence. Unique representative and deterministic candidate 1 remain transient recommendations, not automatic resolutions.
May retire: The requirement that replacement or Skip be presented before every Python-resolved default alignment. The explicit review action must still present those choices before commitment.
Accepted residual risk: A resolved default plan may commit without per-target inspection. The named explicit review action and reliable Undo/Reset must remain available; ambiguous suggestions still require the full disclosed review.
Owner: satoshikawato
Decision date: 2026-09-25
```

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.surface-scope
Scenario revision: 3
Choice: A / AUTO_APPLY_RESOLVED_WITH_EXPLICIT_REVIEW
Rationale: A fully resolved alignment should complete without making the user inspect and Apply candidate choices that require no decision, while an explicit review action retains control when the user wants to inspect or change the result before commitment.
Must preserve: The exact selected reference from the feature popup and Similarity Groups drawer; a full, accessible review with Select, Skip, per-record orientation, candidate details, and resolution summary when Python reports ambiguity or the user explicitly requests review; typed fully resolved plans, shared Python validation, strict CLI ambiguity rejection, actionable errors, and accurate unsupported-feature disclosure.
May retire: Mandatory opening of the full review and manual Apply for a Python-resolved default alignment; the single-action-only Web entry surface, to permit a named explicit review action.
Accepted residual risk: The default resolved alignment may commit before a user notices that they wanted to change direction or Skip. The explicit review action must be discoverable and accessible, and Undo/Reset must remain reliable.
Owner: satoshikawato
Decision date: 2026-09-25
```

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.choice-and-retry
Scenario revision: 3
Choice: A / AUTO_APPLY_RESOLVED_WITH_REVIEW_RETRY
Rationale: The Web should spend interaction time on genuine ambiguity or requested review, using Python's resolved plan immediately when every target has a deterministic outcome.
Must preserve: Independent target choices, Python ownership of eligibility and final validation, local no-Worker draft editing in review, one batch validation for an edited draft, correction and retry without losing the draft, visible candidate facts and effective orientation in review, last Result and History after failure, Cancel, stale, or superseded work, canvas interaction during review, and existing plan regeneration, Session, Reset, and Undo/Redo meanings. An automatic render failure must expose the retained draft and actionable retry path without a partial commit.
May retire: Presenting and applying an editable draft on every Python-resolved default alignment; a second resolver call solely to confirm an unchanged resolved draft.
Accepted residual risk: The default resolved path commits without an individual pre-commit inspection. A visible explicit review action and reliable Undo/Reset must cover users who need a different anchor, Skip, or direction.
Owner: satoshikawato
Decision date: 2026-09-25
```

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.transform-semantics
Scenario revision: 3
Choice: A / PRESERVE_DEFAULT_WITH_EXPLICIT_REVIEW_ORIENTATION
Rationale: Immediate alignment should preserve current orientation by default, while users who need direction matching can explicitly open a review before committing that choice.
Must preserve: The reference record's position and orientation; every target's vertical position; preservation of current orientation unless Match reference direction is explicitly enabled for that record in review; whole-record reversal only with known opposite displayed anchor strands; preservation for unknown strands; readable text; accurate source-relative rev indication; exact idempotent anchor-center alignment; independent records; and atomic validated commitment.
May retire: Showing per-record orientation controls before every Python-resolved default alignment. No orientation policy or reversal rule is retired.
Accepted residual risk: A user who chooses the default Align may need to Undo or Reset and use explicit review to change orientation. The review action must be visible at both entry points and state its effective per-record outcome.
Owner: satoshikawato
Decision date: 2026-09-25
```

## Product Decision Owner response — proposed wording for Choice B

These are also **unsigned proposals**. Select either the complete Choice A set above or the complete Choice B set below, with any explicit edits. Do not combine their outcome IDs. `PD-OI-026` and `PD-OI-027` remain active without revision if B is selected.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.surface-scope
Scenario revision: 3
Choice: B / COMPACT_RESOLVED_CONFIRMATION
Rationale: A fully resolved alignment should avoid repeatedly displaying candidate cards while keeping one pre-commit confirmation for orientation and Skip.
Must preserve: The exact selected reference from the feature popup and Similarity Groups drawer; one keyboard-accessible review for every Align; visible resolved target outcomes and effective directions; Select and candidate facts through explicit expansion; Skip and per-record orientation before Apply; a full disclosed review for ambiguity; typed fully resolved plans, shared Python validation, strict CLI ambiguity rejection, actionable errors, and accurate unsupported-feature disclosure.
May retire: Showing all candidate cards by default when Python has fully resolved every target. The user must be able to expand candidate details before Apply.
Accepted residual risk: A user may overlook a resolved candidate's details until expanding them, and an extra Apply remains necessary even when no target needs a choice. The compact summary must clearly show each effective result and make expansion accessible.
Owner: satoshikawato
Decision date: 2026-09-25
```

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.choice-and-retry
Scenario revision: 3
Choice: B / COMPACT_RESOLVED_LOCAL_REVIEW
Rationale: The Web should show a shorter review for already resolved targets while retaining local choice editing, validation at Apply, and reliable recovery.
Must preserve: Exact reference identity; record-independent Select or Skip and orientation choices; Python ownership of candidate eligibility and the final plan; local no-Worker edits; one batch Python validation at Apply; correction and retry without losing the draft; candidate facts and recommendation reasons on accessible expansion; visible effective orientation before Apply; last Result and History after failure, Cancel, stale, or superseded work; canvas interaction in full review; and existing plan regeneration, Session, Reset, and Undo/Redo meanings.
May retire: Default expansion of candidate cards for Python-resolved targets. No pre-commit Apply or target choice is retired.
Accepted residual risk: The compact view may hide useful evidence until expansion and still requires confirmation for a resolved plan. Outcome and expansion controls must remain clear on narrow screens and by keyboard.
Owner: satoshikawato
Decision date: 2026-09-25
```

## Authority and implementation sequence after a decision

1. Verify the current remote `origin/dev` authority and ancestry; re-evaluate this pack if newer decisions changed the relevant concerns.
2. Serialize only the owner's explicit complete responses for every affected concern in an authority-only change. Review and merge that authority into `origin/dev` before any dependent runtime candidate is admitted.
3. Implement the selected Web outcome in the existing alignment controller and add focused browser checks for both `resolved` and `ambiguous` responses, explicit orientation, failure/retry, stale work, and History/Session invariants.
4. Fix the independent Definition regression in the shared Linear placement path; verify drawing, collision bands, gap, Lock-off behavior, and the translated browser fixture.
5. Review production, test, documentation, and generated diffs separately; run focused checks and the required project gates. Do not alter tracked reference diagrams unless an intentional, reviewed geometry change requires their regeneration.
