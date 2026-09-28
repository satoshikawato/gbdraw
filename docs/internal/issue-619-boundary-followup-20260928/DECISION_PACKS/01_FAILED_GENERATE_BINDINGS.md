# Conditional Product Decision Pack 01 — bindings after failed Generate

Status: **PRODUCT_DECISION_REQUIRED; no Product Decision Owner receipt or approved authority**. S01 reproduced a persisted Session difference after failed Generate. The two outcomes below remain product-valid under the merged authority search. This Pack remains a proposal until an explicit complete human receipt and the required authority sequence. See the [Product Impact Ratchet](../../PRODUCT_IMPACT_RATCHET.md) and [master plan](../MASTER_PLAN.md).

## Identity and trigger

- Concern key: `web.generate.failed-source-binding-continuation` (developer preflight key; not yet a registered concern or `BD-###`).
- Discovery lane: Product Impact Ratchet developer preflight (unmapped material Save continuation).
- Scenario revision: `1` for this proposed outcome comparison.
- Base: `origin/dev` `c2818ce72168e3a35124468e41bb869623ac3148`; refresh the search if dev moves.
- Proposed implementation branch: `fix/issue-619-boundary-followup-20260928`.
- Related scope: Issue #619 follow-up, failed Generate source preparation (finding 5).
- Prepared by: implementation-plan author. Product Decision Owner receipt: **pending**.

A user has a committed diagram and an editable draft. A new Generate begins to resolve record or annotation bindings but fails before admitting a new Result. Existing evidence shows that scalar draft, History, committed request, and Result survive, while some internal bindings may be supplemented. S01 confirmed that subsequent Save Session and fresh Load retain those new bindings. The unresolved Product question is whether this persisted continuation should remain supported or whether binding changes made solely by the failed attempt must be absent.

## Authority search and current evidence

| Source | Result and limit |
| --- | --- |
| Product Impact map / active `BD-###` | No registered decision is claimed for this exact failure-binding choice. Recheck current base before activation. |
| `PD-OI-037` | Failed Generate keeps previous Result/History and pending draft; it does not explicitly decide generation-only binding enrichment in a subsequent Save. |
| `PD-OI-045` | Save/Load must represent one coherent document, preserve source bytes/cache/evidence/provenance, and recover atomically; it does not authorize a mixed, stale, or dangling binding. |
| `PD-OI-046` | Error guidance and retry must preserve request, draft, Result, History, cancel/stale/superseded continuation. |
| Public Session compatibility | Save preserves committed Result and editable state separately; Load restores saved values. Source bindings must remain valid and source bytes exact. |
| Scientific/integrity rules | Record identity, annotation target, resource identity, and actual source bytes cannot change merely to satisfy a Product choice. |
| Current code/tests | S01 browser evidence on source `807f579082db9e0f60d71f2454f2547e4c179d1f` confirms an actual saved-document change; code and tests remain evidence, not authority. |

Authority search result: **UNRESOLVED**; no conflicting active authority found. Procedural classification **now: `PRODUCT_DECISION_REQUIRED`**. On a current writer Session 44 / request schema 8 from the tracked tobacco Gallery fixture, a post-render/pre-admission fault changed `config.adv.multi_record_positions` from `[]` to `[{"selector":"#1","row":1}]` and added `_gbdraw_web_target_record_key` to each of four annotation metadata rows. These five fields were present in the actual failed-attempt download and its fresh Load/re-save; the only other before/after saved-document path difference was `createdAt`. The canonical request, prior Result, embedded resource bytes, scalar draft, History, and annotation record IDs were unchanged. Retry succeeded. The one-record binding key resolved to the same `NC_001879.2` source record; this evidence does not show a different rendered diagram. The changed user-owned Session file and future continuation still distinguish A from B. Merged authority protects coherent Save/Load and failure recovery but does not select whether valid generation-only binding enrichment persists, so neither choice can be selected by the implementation branch. `NOT_ALLOWED` applies to any candidate with stale/dangling/wrong source identity, lost draft/request/Result/History, leaked private data, or weakened required tests. Exact commands, input hashes, and limits are in [S01](../SESSION_RESULTS/S01.md).

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

S01 used an actual browser Save/download → fresh Load before and after one fault injected after binding preparation and before Result admission. The S01 result compares source byte hashes, binding referents, annotation targets, draft, request, Result, History, and retry. Export and same-boundary cancel/stale binding persistence remain S04 acceptance work; the existing Linear cancel/stale guard passed separately in S01. Record exact source SHA, fixture, injection point, measurements, and limitations. An evidence-only commit must not change runtime selection, authority, schema, expected-output baseline, deadlines, or retries.

**Engineering recommendation: Choice B; not an approval.** It aligns generation-only binding changes with the already atomic Result admission using a small provisional candidate. This recommendation is not Product authority and must not be implemented while `PRODUCT_DECISION_REQUIRED` remains unresolved.

## Proposed response for Product Decision Owner, only if activated

If the Product Decision Owner chooses the recommendation after reviewing the updated Pack and S01 evidence, they may send the following completed response in chat. They must supply their own identity and decision date and may edit any proposed term; a choice letter alone is insufficient.

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

After an explicit complete response, serialize **only** the chosen outcome; show the generated machine representation for review. A material public-contract choice requires authority-only integration into `dev` before dependent runtime. Do not infer missing rationale, retirement, risk, owner, or date from this proposal or from plan-branch approval.
