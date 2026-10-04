# Issue #619 follow-up: Session, History, Generate, and replay boundaries

Status: implementation plan, not runtime acceptance or Product authority.
Base: `origin/dev` at `c2818ce72168e3a35124468e41bb869623ac3148` (fetched 2026-09-28).
Implementation branch: `fix/issue-619-boundary-followup-20260928`.
Dedicated worktree: `.worktrees/issue-619-boundary-followup-20260928` at the repository root.

## Purpose and starting state

Issue #619's Circular Width/Radius controls and the related diagnostic and CI repairs are already integrated into the base. Runtime/test [PR #635](https://github.com/satoshikawato/gbdraw/pull/635) and CI/test [PR #633](https://github.com/satoshikawato/gbdraw/pull/633) were merged separately. The exact dev merge SHA above passed [hosted Tests run 36354835871](https://github.com/satoshikawato/gbdraw/actions/runs/36354835871), including Dev staging, and [Gallery run 36354835863](https://github.com/satoshikawato/gbdraw/actions/runs/36354835863), including readiness. Its 373 functional cases had 372 first-attempt passes and one retry pass. The base commit and current repository contracts are the implementation authority; the historical publication record is optional context, not a prerequisite.

This follow-up addresses five bounded findings and one browser-test reliability finding. The former are numbered as in the original [Issue #619 final result](../issue-619-implementation-plan-20260927/SESSION_RESULTS/S03.md):

1. A biological Session saved while the Linear editor is active can reload into the committed Circular mode. Its entire editable config was not proved equal.
2. A nonfinite JavaScript **number** injected into a Circular Width/Radius draft can become `null` through JSON-based History capture. Invalid text such as `"Infinity"` is a different, already retained draft.
4. Restoring a History config domain can invoke the existing `rulePreparation` helper; the necessity and side effects of that call need characterization.
5. Source preparation during a failed Generate can add record/annotation bindings. The proven rollback covers scalar, History, request, and Result; whole-config byte identity was not proved.
6. Native replay has not proved SVG byte identity or full sidecar JSON equality. The published contract expressly permits SVG/text-metric differences across versions and does not promise byte-identical Exact replay across font environments.

The hosted retry pass was `tests/web/preview-navigation.playwright.spec.js`, background pan after text selection and zoom. At zoom 0.3 the pan state moved exactly `(60,80)` and scroll stayed unchanged, while screen displacement was about `(60.5802,80)`. The SVG screen scale also changed between the two samples. Treat this as a focused stability investigation, not evidence that pan semantics changed.

OS Japanese IME testing and repair (original item 3) are explicitly **out of scope** for this plan.

## Supported behavior and acceptance boundaries

| Concern | Existing authority and required outcome | What this plan may change |
| --- | --- | --- |
| Load | [Session and request compatibility](../../REFERENCE/session-and-request-compatibility.md) says Save stores editable state separately from committed Result and Load restores saved values. `PD-OI-045` requires one coherent, atomic Load and failed-load recovery. | Restore the saved active editor mode and draft alongside the saved committed Circular Result when the current writer contains both. Keep settings-only, CLI-origin, historical reader, and failed Load behavior. |
| Circular scalar/History | `PD-OI-048`–`050` and the existing `parseOptionalCircularScalar()` contract admit finite positive px/factor values or Auto. Invalid/incomplete text stays visible and cannot be silently converted to Auto. | Refuse nonfinite **numeric** values before the owning editor action mutates a slot or records History. Keep invalid text in its existing draft path. Do not change generic JSON cloning or Session schema. |
| History restore | `gbdraw/web/CLAUDE.md` requires the existing History and rule-preparation owners, canonical override state, and correct Result after Undo/Redo. | Measure whether `rulePreparation` is necessary on a config restore. Keep one preparation when it changes visible rules; remove only a demonstrated redundant invocation. |
| Failed Generate | `PD-OI-037`, `PD-OI-045`, and `PD-OI-046` preserve the previous Result, request, edit draft, History, source bytes, and recovery actions. | Characterize binding enrichment before selecting a fix. If an observable guarantee is violated, alter the existing preparation/transaction owner only. Do not introduce a whole-app clone or generic rollback. |
| Native replay | The public Session contract promises semantic reconstruction, not cross-environment byte identity. Scientific identity, scalar value/unit, source bytes, and geometry remain exact or bounded by a documented numeric tolerance. | Strengthen the existing semantic comparison evidence and explicitly classify allowed serialization differences. No new byte-equality promise or writer migration. |
| Pan retry | Existing pan behavior is measured through state, scroll, DOM transform, and SVG screen coordinates. | Stabilize the pre-move geometry observation if it is a sampling race; retain the existing displacement tolerance and test-owned deadline. Fix runtime only if a stable reproduction proves a real pan defect. |

The source of Product authority is the latest **merged base**, not this plan, current code alone, or tests alone. `PD-OI-###` denotes the existing static Product Contract, not a `BD-###` decision. Re-run Product Impact preflight against the actual implementation base before each user-visible change. A material outcome not selected by merged authority stops only its dependent runtime work.

## Architecture and design rules

The canonical Web route stays `session-request.js` → `run-analysis.js` → `diagram-generation.js` → one diagram Worker → typed request decoder → sanitized Result. `services/config.js` owns Save/Load coordination; `session-active-config-contract.js` owns current writer validation; `services/history.js` and `history-snapshot.js` own History; `app/track-slot-validation.js` owns Circular scalar meaning; `app/circular-track-slots.js` owns its editor action. `rulePreparation` remains the rule evaluation boundary. Native replay uses the existing typed Session and renderer path.

- **SRP / ISP:** Each fix belongs to the smallest existing owner above. UI code consumes scalar interpretation and Session operations through their existing narrow interfaces. Avoid making the Session coordinator a second scalar parser or making History a source-binding validator.
- **OCP / DIP:** Extend existing contracts where new input evidence requires it. Call the established owner from the editor, instead of adding a mode-specific alternate request path or a new render service.
- **LSP:** A supported current/historical Session and an accepted scalar must retain its existing meaning, failure behavior, and replay path. New guards must distinguish invalid numeric injection from valid or unfinished text drafts.
- **DRY:** One active editor mode, one source of canonical request truth, one scalar parser, one History owner, and one rule-preparation path. Remove a superseded decision point in the same change if one is discovered.
- **KISS / YAGNI:** No new Session version, compatibility reader, generalized snapshot framework, global rollback, parallel Worker, new reactive mirror, or byte-normalizing SVG writer. Do not build a feature for an unproved user effect.

For architecture-bearing changes, record concise before/after owner/path evidence, user effects, checks, and rollback under the [architecture ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md). Exact OE/PE/CB tables are needed only if its exception conditions arise. Keep production, tests, documentation, and generated diffs separate in review.

## Product Decision treatment

Items 1, 2, 4, and 6 are planned as `IMPLEMENT_EXISTING_AUTHORITY`. Item 5 begins as `EVIDENCE_REQUIRED`: binding enrichment may be an internal representation improvement or may alter a saved user's next action. [Conditional Decision Pack 01](DECISION_PACKS/01_FAILED_GENERATE_BINDINGS.md) contains a recommended outcome and a complete response template. It is **not** an approved choice. Activate its `PRODUCT_DECISION_REQUIRED` route only if reproducible evidence leaves two materially different product-valid outcomes after searching merged authority. If the selected outcome changes the public contract, merge authority separately into `dev` before dependent runtime; candidate authority cannot authorize its own runtime. Do not infer a signature from approval of this implementation-plan branch.

No Product Decision is requested for SVG byte identity: the existing published contract does not promise it. A future proposal to add that guarantee is outside this plan and would need its own evidence and preflight.

## Sequential sessions and branch discipline

Each prompt under [`SESSION_PROMPTS/`](SESSION_PROMPTS/) is standalone and names its own outputs. Run them in order; each session ends with one scoped commit and a push to the **same-named remote work branch**. A session that finds no justified runtime change commits its evidence and disposition instead. Never skip a failing acceptance check by changing its deadline, retry count, or oracle.

**S01 must start by checking out `fix/issue-619-boundary-followup-20260928` inside the dedicated worktree above.** The root checkout and other worktrees may belong to other sessions. S01 must inspect status/HEAD/upstream and confirm its branch is based on `c2818ce`; it must not switch the root checkout, reset, clean, prune worktrees, or stop another session's process. S02 and later **reuse the branch first checked out in S01**, in that same worktree. They fetch and fast-forward the same-named remote branch before work; they do not create another branch or worktree for implementation. If another actor is writing that branch, coordinate through Git state rather than overwriting it. Before every push, verify branch and upstream; push only `HEAD:refs/heads/fix/issue-619-boundary-followup-20260928`, with no force.

| Session | Owned result | Completion evidence |
| --- | --- | --- |
| [S01](SESSION_PROMPTS/S01_BASELINE_AND_PREFLIGHT.md) | Reproduce/classify all five findings and pan retry on frozen base; record Product preflight and binding-effect matrix. | `SESSION_RESULTS/S01.md` with commands, source SHA, fixtures, before/after state, privacy-safe evidence, unresolved decision status. |
| [S02](SESSION_PROMPTS/S02_CROSS_MODE_LOAD.md) | Fix active-mode Load at existing Session owner. | Focused current/historical/settings-only/failed Load, Save/fresh Load/re-save, request/Result/draft/History/resource checks. |
| [S03](SESSION_PROMPTS/S03_SCALAR_HISTORY_AND_RULES.md) | Guard nonfinite numeric editor writes and characterize/fix rule preparation. | Node and real-browser Undo/Redo, valid/invalid text, no extra Worker, rule-derived color parity and call count. |
| [S04](SESSION_PROMPTS/S04_FAILED_GENERATE_BINDINGS.md) | Resolve failed Generate binding effect according to S01 evidence and any required approved authority. | Fault-injected failure/retry/cancel/stale matrix and actual Save/fresh Load where relevant. |
| [S05](SESSION_PROMPTS/S05_REPLAY_AND_PAN.md) | Validate semantic native replay; make pan test deterministic if evidence supports a sampling race. | Pinned source/record/scalar/geometry checks; no-retry focused pan repetitions with unchanged tolerance. |
| [S06](SESSION_PROMPTS/S06_INTEGRATION_AND_HANDOFF.md) | Review combined source and complete branch acceptance. | Applicable required gates, exact four-shard inventory, first-attempt/retry distinction, Product/architecture report, final status and handoff. |

`SESSION_RESULTS/` is created by implementation sessions. Each result must be self-contained enough for the next session: exact commit, environment, commands, outcomes, files changed, owner/path effects, known limits, and next action. Do not rely on temporary `/tmp` files as the only evidence. Avoid committing private genome data, downloads, traces, generated wheels, or large logs; retain reproducible commands and bounded summaries in the work branch.

## Verification and stopping rules

Use focused Node/Python/browser checks in S02–S05; check both Node and Python Playwright availability when browser work matters. Run the full functional set and required policy/architecture gates once on the frozen combined candidate in S06, broadening only for a changed source or concrete failure. Preserve all expanded cases, four disjoint shards, two CI retries, the 45-minute CI job budget, and test-owned deadlines. The first-attempt and retry-pass counts must be reported separately. `tests/reference_outputs/`, the social preview, Gallery assets, signed Product authority, generated wheel, and `dist/` are outside this plan unless a separately justified scope change requires them.

S06 completes a reviewable **branch candidate**. Opening a PR, merging into `dev`, deployment, tagging, closing Issue #619, or changing Product authority needs its own applicable authorization. Hosted acceptance after any future merge must be tied to that exact merge SHA; branch or older SHA success is not a substitute.

Proposed implementation commit themes: restore saved active mode on Load; reject nonfinite numeric Circular drafts before History; preserve Generate failure continuation at existing owners; verify semantic replay and stable pan measurement. Each session supplies its own English commit title and short summary.
