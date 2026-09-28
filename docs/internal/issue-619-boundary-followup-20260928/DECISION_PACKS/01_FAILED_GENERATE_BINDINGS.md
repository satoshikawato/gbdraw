# Conditional Product Decision Pack 01 — bindings after failed Generate

Status: **RESOLVED — IMPLEMENT_EXISTING_AUTHORITY**. S01 reproduced a persisted Session difference after failed Generate. The Product Decision Owner subsequently selected Choice A in the complete `PD-OI-055` receipt, merged into `dev` by `e5f5f1fb`. Choice B below is the unselected engineering proposal. See the [Product Impact Ratchet](../../PRODUCT_IMPACT_RATCHET.md) and [master plan](../MASTER_PLAN.md).

## Identity and trigger

- Concern key: `web.generate.failed-source-binding-continuation` (developer preflight key; not yet a registered concern or `BD-###`).
- Discovery lane: Product Impact Ratchet developer preflight (unmapped material Save continuation).
- Scenario revision: `1` for the approved Choice A outcome and the historical comparison.
- S01 evidence base: `origin/dev` `c2818ce72168e3a35124468e41bb869623ac3148`. S04 authority base: `origin/dev` `775473f4ee7c825c888f705146bcbec71b11cc4b`.
- Proposed implementation branch: `fix/issue-619-boundary-followup-20260928`.
- Related scope: Issue #619 follow-up, failed Generate source preparation (finding 5).
- Prepared by: implementation-plan author. Product Decision Owner receipt: `PD-OI-055`, approved by `satoshikawato` on `2026-09-28`; authority-only commit `82ad5330dd0dc06ddcca13254e0dd615a845f5bd`, merged by `e5f5f1fb`.

A user has a committed diagram and an editable draft. A new Generate begins to resolve record or annotation bindings but fails before admitting a new Result. Existing evidence shows that scalar draft, History, committed request, and Result survive, while some internal bindings may be supplemented. S01 confirmed that subsequent Save Session and fresh Load retain those new bindings. That saved-document difference posed the Product question now resolved by `PD-OI-055` in favor of retaining only validated same-source enrichment.

## Authority search and current evidence

| Source | Result and limit |
| --- | --- |
| Product Impact map / active `BD-###` | No registered `BD-###` is claimed for this exact choice. The static Product Contract now contains accepted `PD-OI-055`, selecting Choice A for scenario revision 1. |
| `PD-OI-037` | Failed Generate keeps previous Result/History and pending draft; it does not explicitly decide generation-only binding enrichment in a subsequent Save. |
| `PD-OI-045` | Save/Load must represent one coherent document, preserve source bytes/cache/evidence/provenance, and recover atomically; it does not authorize a mixed, stale, or dangling binding. |
| `PD-OI-046` | Error guidance and retry must preserve request, draft, Result, History, cancel/stale/superseded continuation. |
| Public Session compatibility | Save preserves committed Result and editable state separately; Load restores saved values. Source bindings must remain valid and source bytes exact. |
| Scientific/integrity rules | Record identity, annotation target, resource identity, and actual source bytes cannot change merely to satisfy a Product choice. |
| Current code/tests | S01 browser evidence on source `807f579082db9e0f60d71f2454f2547e4c179d1f` confirms an actual saved-document change; code and tests remain evidence, not authority. |

Authority search result on the S04 base: **RESOLVED** by `PD-OI-055` (`A / RETAIN_VALIDATED_BINDING_ENRICHMENT`). Procedural classification: `IMPLEMENT_EXISTING_AUTHORITY`. S01 and S04 actual Save/download → fresh Load show one position entry and four annotation binding keys retained after a post-render/pre-admission fault. The keys resolve to the same `NC_001879.2` source record; source hashes, annotation targets, scalar draft, committed request, Result, and History remain unchanged. S04 also verifies export of the old SVG, canceled/stale completion, and a newly enriched key for an explicit RecB target in a two-record source at the same boundary. A concurrent second Generate is rejected as `GENERATION_BUSY`, so it is not a supported supersession route. Choice A permits this complete validated metadata enrichment; it does not permit incorrect, partial, stale, dangling, or wrong-source bindings. Choice B remains below as the historical comparison and must not be implemented under the selected authority. Exact commands and limits are in [S04](../SESSION_RESULTS/S04.md).

## User journey and checkpoints

- Actor/goal: Web user editing a previously generated genome diagram and trying another Generate without losing current work.
- Entry: committed Result visible; draft and biological sources present; Generate starts.
- Checkpoints: source preparation; failure/cancel/stale/superseded settlement; existing Result and controls; Save Session; fresh Load; retry Generate; export.
- Required next actions: inspect a safe error, correct draft, retry, Save the coherent current document, export the still-current Result, or cancel.
- Performance/privacy: no extra full-document clone, genome-byte duplication, private source logging, or extra Worker.

## Choice A — RETAIN_VALIDATED_BINDING_ENRICHMENT

- **Complete outcome:** A failed Generate may retain newly resolved bindings only when they are validated, source-identical, and coherent with the unchanged draft and prior committed artifact. A subsequent Save/Load may include these bindings even though no new Result was admitted.
- Preserved effects: prior Result/request, scalar and other drafts, History, source bytes, retry/cancel/stale guards, old Export, valid source/annotation identity.
- Added effect: later Save may contain more explicit bindings and may avoid rediscovery on retry.
- Lost/retired effect: a strict guarantee that failed Generate leaves all saved binding metadata unchanged.
- Discoverability/accessibility: no new control; error and retry remain available.
- Canonical state/Undo/Redo: validated bindings join the editable document; no Generate History entry or committed request change.
- Session/regeneration/export: Save and fresh Load retain those bindings; the old Result remains the exported diagram; later Generate uses the same source identity.
- Validation/failure: malformed, stale, partial, or wrong-source bindings fail closed; no private data in diagnostics.
- Scientific/cache/performance: scientific identity unchanged; possibly less retry discovery, with bounded cache reuse and no extra Worker.
- Compatibility/architecture: current schema and single source-binding owner; no second reader or generation path.
- Evidence remaining: multi-record target identity, same-boundary cancel/stale/superseded binding persistence, and resource cost. Single-record Save/fresh Load and retry were observed in S01.
- Residual risk: binding metadata may change after an unsuccessful action and surprise users comparing saved documents.
- Route if material: `DURABLE_AUTHORITY_REQUIRED` before runtime that intentionally establishes this as a public continuation.
- Next action if selected: use the existing preparation owner; prove all coherence and identity constraints, then document the visible Save behavior.

## Choice B — COMMIT_BINDINGS_ONLY_WITH_SUCCESSFUL_GENERATE (recommended)

- **Complete outcome:** Binding changes made solely for a Generate attempt stay provisional until the new Result is successfully admitted. Failure, cancel, stale completion, and supersession leave the saved document's semantic bindings as they were immediately before that attempt. Source discovery independently completed before Generate remains valid.
- Preserved effects: prior Result/request, all editable drafts and History, source bytes, existing cache/evidence/provenance allowed by contract, valid retry, old Export, and user-initiated discovery.
- Added effect: Save/fresh Load after failure reconstructs the same pre-attempt binding choices; retry can resolve provisional bindings again or use a validated existing cache.
- Lost/retired effect: implicit promotion of generation-only binding metadata by an unsuccessful attempt.
- Discoverability/accessibility: no new control; existing error/correction/retry path stays visible and keyboard reachable.
- Canonical state/Undo/Redo: successful admission commits the candidate through the existing transaction; failure creates no Generate History entry or request/Result replacement.
- Session/regeneration/export: Save after failure retains coherent pre-attempt bindings, draft, and old Result; a later successful Generate may publish newly resolved bindings and Result together.
- Validation/failure: stale/partial/wrong-source candidates never commit; privacy-safe bounded diagnostics remain.
- Scientific/cache/performance: source/record/annotation identity exact; small provisional binding data, no full-app snapshot or extra Worker, retry cost measured.
- Compatibility/architecture: current schema and existing preparation/commit owners; no parallel binding store, compatibility reader, or global rollback.
- Evidence remaining: multi-record target identity, whether a binding from independent pre-Generate discovery is already public, same-boundary cancel/stale/superseded persistence, retry cost, and cache behavior. S01 located the post-render/pre-admission fault and saved difference.
- Residual risk: a retry may redo bounded binding preparation; measure it and retain valid cache reuse without promoting candidate document state.
- Route if material: `DURABLE_AUTHORITY_REQUIRED` before dependent runtime.
- Next action if selected: merge the exact human-approved durable authority into `dev`, then implement candidate staging at the existing owner on this work branch.

## Comparison

| Dimension | A: retain validated enrichment | B: commit on success |
| --- | --- | --- |
| Immediate visible diagram and controls | Prior Result and draft remain | Prior Result and draft remain |
| Save/fresh Load after failed Generate | May include new validated bindings | Keeps pre-attempt semantic bindings |
| Retry | May reuse committed bindings | Re-resolves or uses validated cache |
| Undo/Redo | No failed Generate entry | No failed Generate entry |
| Export/scientific content | Prior Result; exact source identity | Prior Result; exact source identity |
| Failure/cancel/stale | Candidate may enrich document if valid | Candidate document changes discarded |
| Compatibility | No new schema | No new schema |
| Main risk | Unsuccessful action changes saved metadata | Bounded repeat preparation cost |
| Route if both remain valid | Durable authority | Durable authority |

## Evidence-first work and engineering recommendation

S01 used an actual browser Save/download → fresh Load before and after one fault injected after binding preparation and before Result admission. The S01 result compares source byte hashes, binding referents, annotation targets, draft, request, Result, History, and retry. S04 completed the actual Export, two-record target, and same-boundary cancel/stale Save/fresh Load checks; the existing Linear guard was also rerun. The exact source SHA, fixture, injection point, measurements, and limitations are in S04. No runtime selection, authority, schema, expected-output baseline, deadline, or retry was changed by the S04 commit.

**Historical engineering recommendation: Choice B.** The Product Decision Owner selected Choice A in `PD-OI-055`; the recommendation is superseded and is not implementation authority.

## Historical proposed Choice B response (unselected)

This pre-decision template is retained for the option comparison. The actual complete Choice A receipt and nine-field representation are in `PD-OI-055`; this Choice B text is not a request for another decision.

```text
PRODUCT_DECISION
Concern: web.generate.failed-source-binding-continuation
Scenario revision: 1
Choice: B / COMMIT_BINDINGS_ONLY_WITH_SUCCESSFUL_GENERATE
Rationale: A failed Generate should leave the saved document's binding choices coherent with the prior committed Result and the current editable draft, so Save and retry do not silently adopt changes from an unsuccessful attempt.
Must preserve: The previous Result, canonical request, all editable drafts, History, source bytes, exact record and annotation identity, valid source discovery completed before Generate, retry, Save/Load, Export, cancel/stale/superseded recovery, bounded privacy-safe diagnostics, and one existing Worker and generation path.
May retire: Promotion into saved document state of binding changes made solely by a failed, canceled, stale, or superseded Generate attempt; no supported source discovery, cache reuse, or user edit is retired.
Accepted residual risk: Retry may repeat bounded binding preparation. Measure the cost and retain validated cache reuse without committing an unsuccessful candidate or duplicating full document state.
Owner: <Product Decision Owner identity>
Decision date: <YYYY-MM-DD>
```

The approved Choice A representation was merged separately into `dev` before S04 verification. This Pack does not change or expand that authority.
