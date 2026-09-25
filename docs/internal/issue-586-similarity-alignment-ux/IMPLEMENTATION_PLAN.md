# Issue #586: Similarity Group alignment review UX implementation plan

Status: approved product outcome; implementation not started
Issue: <https://github.com/satoshikawato/gbdraw/issues/586>
Implementation branch: `issue-586-similarity-alignment-ux-20260924`
Branch base at plan creation: `origin/dev@da692d635f609f49bbb31357c5fd529b661c8ee6`
Product Decision Owner: `satoshikawato`
Decision date: `2026-09-24`

## 1. Purpose

Linear diagrams can align records by a feature in a Similarity Group. The
current Web workflow exposes separate `Align` and `Align & orient` actions,
performs work before the initiating controls visibly enter a busy state, and
opens the ambiguity palette with no selection. A user must therefore click one
choice for every ambiguous record even when the application can make a stable,
fully disclosed recommendation.

This change replaces that workflow with one `Align…` action and one review
palette. The palette always opens before a transform is committed, contains a
complete row for every displayed target record, preselects deterministic anchor
recommendations, and offers a per-record `Match reference direction` control.
The draft is edited locally with no Worker round trip. One final `Apply` sends
the complete draft to the shared Python resolver for validation and then
regenerates the diagram atomically.

This document is the master plan for a new contributor. It does not assume
knowledge of the issue discussion or of any earlier implementation session.

## 2. Required branches and authority order

All runtime implementation sessions use exactly this branch:

```text
issue-586-similarity-alignment-ux-20260924
```

Do not create a replacement runtime branch and do not implement this feature on
`dev`, `main`, or an unrelated topic branch. Before each session, fetch the
remote, verify the current branch and upstream, and preserve unrelated changes.
Push only to the same-named remote work branch.

The signed decisions below change four active Product Contract outcomes. The
repository policy requires them to be merged as authority before dependent
runtime code can be accepted. Session 00 therefore uses a separate,
authority-only branch created from the then-current `origin/dev`. It updates the
Product Contract and nothing else. After that authority is reviewed and merged
into `origin/dev`, merge the updated `origin/dev` into the fixed implementation
branch above. Sessions 01–05 must not begin while the superseding authority is
absent from `origin/dev`.

The plan files themselves live on the implementation branch. Their presence is
not authority for runtime behavior.

## 3. Signed product decisions

The following receipts are the complete human decisions. Session 00 must
serialize them without inventing rationale, retirement intent, or accepted
risk.

### 3.1 Anchor resolution

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.anchor-resolution
Scenario revision: 2
Choice: A / AUTO_PRESELECTED_SUGGESTIONS
Rationale: Similarity Group alignment should remain reviewable without requiring one repetitive click for every ambiguous record. A deterministic recommendation gives users an immediately applicable draft while exact identity, visible recommendation reasons, Select, Skip, and final Python validation preserve transparency and control.
Must preserve: The exact selected reference; stable record and biological-feature identity; only-usable-candidate and unique-direct-RBH automatic resolution; independent treatment of every displayed record; unchanged position and orientation for missing, unusable, or skipped records; visible recommendation reasons; the ability to replace every recommendation or select Skip; final validation by the shared Python resolver; and independence from viewport, scroll, ribbon geometry, confidence score, supporting-edge count, and multi-hop evidence.
May retire: The initial unselected state for every ambiguous record; the requirement to click every ambiguous record before Apply; and the prohibition on using a unique representative or deterministic candidate 1 as a disclosed transient recommendation.
Accepted residual risk: A recommended candidate is a convenience heuristic rather than proof of biological superiority, and a user may Apply without inspecting every preselection. The palette must disclose the recommendation basis, permit replacement or Skip, and commit nothing before validated Apply.
Owner: satoshikawato
Decision date: 2026-09-24
```

### 3.2 Transform semantics

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.transform-semantics
Scenario revision: 2
Choice: A / SINGLE_ALIGN_PER_RECORD_ORIENTATION
Rationale: Anchor selection and orientation should be reviewed together so users can preserve or match orientation independently for each record instead of committing to one global orientation mode before seeing the candidates.
Must preserve: The reference record's position and orientation; every target's vertical position; current orientation unless Match reference direction is explicitly enabled for that record; whole-record reversal only when both displayed anchor strands are known and opposite; orientation preservation for unknown strands; readable text; effective source-relative rev indication; exact idempotent anchor-center alignment; independent handling of every displayed record; and atomic validated Apply.
May retire: Separate Align and Align & orient actions; one global orientation mode for all targets; and orientation intent that becomes fixed before the review palette opens.
Accepted residual risk: Always opening the review palette adds one confirmation step, and record-level controls add visual density. Orientation defaults to preservation, unknown strands cannot trigger reversal, and the palette must state each record's effective outcome before Apply.
Owner: satoshikawato
Decision date: 2026-09-24
```

### 3.3 Public surface

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.surface-scope
Scenario revision: 2
Choice: A / SINGLE_REVIEW_ALIGNMENT_SURFACE
Rationale: The Web workflow should expose one predictable Align entry point and one review surface containing anchor and orientation choices, while programmatic surfaces continue to consume strict typed plans.
Must preserve: Exact reference selection in the feature popup and Similarity Groups drawer; Select and Skip; a reviewable resolution summary; typed fully resolved Python plans; strict CLI rejection of unresolved ambiguity; shared Python validation; actionable errors; and accurate disclosure that Collinear alignment controls, anchor TSV, scored inference, and multi-hop automatic selection remain unsupported.
May retire: Separate Web Align and Align & orient buttons; automatic application without opening the review palette; and the requirement to choose global orientation intent at the initiating surface.
Accepted residual risk: A review step is required even when every anchor resolves uniquely. The palette must open promptly, remain keyboard- and narrow-viewport accessible, and provide an immediately applicable default draft.
Owner: satoshikawato
Decision date: 2026-09-24
```

### 3.4 Web choice and retry

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.choice-and-retry
Scenario revision: 2
Choice: A / PRESELECTED_LOCAL_BATCH_REVIEW
Rationale: Users should be able to review all target records, accept deterministic defaults, and adjust anchors or orientation without per-choice Worker delays. One final Apply should validate and generate the complete result.
Must preserve: Exact reference identity; record-independent Select and Skip; Python ownership of candidate eligibility and the final plan; local no-Worker editing of anchor and orientation choices; one batch validation at Apply; correction and retry without losing the draft; visible biological names, coordinates, strand, representative status, and direct evidence; last-Result and History preservation after failure, Cancel, stale, or superseded work; canvas interaction; and existing plan regeneration, Session, Reset, and Undo/Redo meanings.
May retire: The zero-selected initial draft; mandatory per-record selection before Apply; automatic application when no ambiguity exists; separate global orientation actions; and discarding review state after a retryable Apply failure.
Accepted residual risk: Apply still waits for one resolver validation and any required regeneration, and automatically selected ambiguous anchors may be accepted without individual inspection. Recommendation reasons and effective orientation outcomes must remain visible and editable until Apply.
Owner: satoshikawato
Decision date: 2026-09-24
```

## 4. Product behavior specification

### 4.1 Entry points and feedback

- The feature popup and Similarity Groups drawer each expose one `Align…`
  action for the exact selected member. There is no Web `Align & orient` action.
- Invocation keeps the exact selected record and biological feature as the
  reference. A representative member must never silently replace it.
- The initiating action immediately becomes disabled and shows an in-progress
  label such as `Resolving…`. The busy state is also available to assistive
  technology. Repeated activation while the operation is active is ignored.
- Resolver errors are actionable and leave the current diagram, last successful
  Result, History, and active alignment unchanged.

### 4.2 Review draft

- The review palette opens after the initial Python resolution for every valid
  request, including requests with no ambiguity.
- It lists every displayed non-reference record in stable displayed-record
  order. Missing or unusable members appear as unchanged and cannot acquire an
  anchor by UI inference.
- Each usable row has one selected anchor or `Skip`. Deterministic defaults make
  the initial draft immediately applicable.
- A row explains why its current default was selected. The finite reasons are:
  `only usable candidate`, `unique direct RBH`, `unique representative`, and
  `deterministic candidate 1`. The first two remain resolver-owned automatic
  resolutions; the latter two are disclosed transient recommendations for
  otherwise ambiguous records.
- Candidate details retain biological name or feature ID, source coordinates,
  displayed strand, representative status, and direct-evidence information.
- A user can replace any selected candidate or choose `Skip`. Selection changes
  are local UI state and start no Worker.
- Each usable target row has `Match reference direction`, default off. Its
  effective result is stated before Apply: preserve, reverse whole record, or
  preserve because a required strand is unknown. The control cannot imply that
  an unknown strand will reverse a record.
- The existing canvas candidate interaction remains available. Clicking a
  candidate ribbon and selecting its row control invoke the same local draft
  mutation; neither path owns selection policy.

### 4.3 Apply, cancel, retry, and stale work

- `Apply` submits the entire explicit draft once. The Python resolver rechecks
  reference identity, record coverage, candidate eligibility, selection,
  orientation policy, crop/display facts, and the final typed plan.
- No record transform is committed until validation and regeneration both
  succeed. A successful operation enters the existing canonical Result and
  History path as one user action.
- `Cancel` closes the palette without changing the diagram, active plan, last
  Result, or History.
- A retryable validation or generation failure leaves the palette and draft
  intact. The user can correct it and Apply again.
- A stale or superseded completion cannot mutate current state. It also cannot
  erase the last successful Result or History.
- Exact anchor-center alignment is idempotent. Reapplying the same validated
  plan does not drift horizontally or toggle orientation.

### 4.4 Orientation rules

For each target record independently:

| Requested policy | Reference/target displayed strands | Effective outcome |
| --- | --- | --- |
| preserve | any | keep current orientation |
| match reference | both known and equal | keep current orientation |
| match reference | both known and opposite | reverse the whole target once |
| match reference | either unknown | keep current orientation |
| skipped/missing/unusable | any | keep position and orientation |

The reference record never moves or reverses. Target vertical positions never
change. Horizontal translation is calculated from final displayed anchor
centers after any effective reversal. Text remains readable through the
existing rendering path. The inspector's `rev` value remains derived from
source-relative effective orientation, not from the checkbox value alone.

## 5. Architecture and ownership

### 5.1 End-to-end flow

```text
exact selected reference
        |
        v
Web controller starts one resolver job and exposes busy state
        |
        v
Python shared resolver validates facts and returns per-record outcomes,
candidates, stable recommendation, and reason
        |
        v
Web builds one local review draft; user edits selection/orientation locally
        |
        v
Apply sends every explicit row once
        |
        v
Python shared resolver validates and emits one typed schema-2 plan
        |
        v
existing canonical regeneration -> Result admission -> History/Session
```

### 5.2 Single owners

| Responsibility | Owner | Consumers |
| --- | --- | --- |
| Candidate eligibility, canonical ordering, RBH interpretation, recommendation and reason | `gbdraw/layout/similarity_alignment.py` | CLI, typed Python, Web adapter |
| Effective per-record reversal and validated typed plan | same shared Python domain module | render planning, Session codec, Web |
| JSON projection across the Worker boundary | `gbdraw/web_support/similarity_alignment.py` | Web controller |
| Operation state, local review draft, stale-result admission | `gbdraw/web/js/app/similarity-alignment.js` | popup, drawer, palette |
| Markup and presentation | existing `gbdraw/web/index.html` templates/CSS and focused Web view helpers | browser |
| Request/Session serialization | Python and Web request/session codecs at their existing canonical boundaries | save/load/regeneration |
| Record transform/materialization | existing record-planning and diagram assembly owners | all render surfaces |

The JSON adapter may rename or serialize typed fields, but it must not rank
candidates, derive recommendation reasons, or decide reversal. The HTML
template may render state, but it must not duplicate controller transitions or
biological rules. Popup and drawer controls call one controller entry point.

### 5.3 Typed domain changes

Use one per-record orientation enum with two values, conceptually:

```text
preserve
match_reference
```

The submitted record choice carries stable record identity, `Select` or `Skip`,
the exact anchor identity when selected, and the requested orientation policy.
The final plan records both the requested orientation policy and the effective
reverse-complement result. This distinction is necessary because
`match_reference` with equal or unknown strands and `preserve` can all produce
the same effective orientation but are different user intents.

The plan format advances from schema 1 to schema 2 and removes the global
orientation mode from the current representation. Schema 2 contains per-record
requested policy and effective outcome. Reference and skipped rows are
constrained to preservation. Do not add a second plan type or an optional
parallel `mode` path.

An initial resolution may still be incomplete as a final plan, but it returns
enough typed data to construct every review row. A recommendation is not a
committed decision. Only the final explicit batch can produce the plan applied
to the diagram.

### 5.4 Recommendation algorithm

The resolver keeps the existing selection precedence:

1. validate and use an explicit `Select` or `Skip`;
2. select the only usable candidate;
3. select the one distinct candidate with a unique direct reciprocal-best-hit
   relationship to the exact reference;
4. otherwise report ambiguity.

For step 4 only, return one transient recommendation:

1. when exactly one ambiguous candidate is marked representative, recommend it
   with reason `unique representative`;
2. otherwise recommend candidate 1 from the resolver's existing canonical,
   identity-based ordering with reason `deterministic candidate 1`.

Recommendation order must not depend on viewport, scroll, SVG/ribbon geometry,
confidence score, supporting-edge count, or multi-hop evidence. JavaScript must
not reproduce this algorithm. The final Apply converts the preselection into an
explicit choice and the resolver validates it like any user-edited choice.

### 5.5 Persistence compatibility decision

At plan creation, `SIMILARITY_ALIGNMENT_PLAN_SCHEMA = 1` exists on `dev` but its
introducing work is absent from `origin/main` and release tags. The latest
release tag is older than the typed plan. Under the repository's persisted
format policy, this is a branch-owned intermediate schema, not a released
compatibility obligation.

Session 01 must repeat that ancestry/tag audit before editing. If schema 1 is
still unreleased, replace it with schema 2, remove the superseded schema-1
reader/tests/docs, and regenerate branch-owned fixtures through their owning
tooling. Do not create a migration solely for an unreleased intermediate
format. If the audit instead finds that schema 1 has shipped or entered an
accepted compatibility commitment, stop the schema work and prepare a bounded
reader-only migration plan before continuing.

Known branch-owned artifacts to inventory include:

- `gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json`;
- similarly named Session fixtures under `tests/test_inputs/`;
- typed examples and exact codec fixtures that embed plan schema 1.

Never hand-edit generator-owned Gallery Session JSON. Use the existing Gallery
refresh/publication tools and review the generated diff.

### 5.6 State machine

Keep one controller state machine. Equivalent names are acceptable if the
transitions and invariants remain explicit:

```text
idle -> resolving -> reviewing -> applying -> idle
                     ^             |
                     |--- failure--|
```

- `resolving` and `applying` are busy; `reviewing` is locally editable.
- Failure during initial resolution returns to a recoverable idle/error state.
- Failure during Apply returns to `reviewing` with the same draft.
- Cancel from `reviewing` returns to `idle` without a write.
- A monotonically increasing operation token or the existing equivalent rejects
  stale completions.
- Repair/regeneration of an already active plan uses the same plan owner and
  admission rules; it must not create a second hidden alignment workflow.

## 6. SOLID, KISS, DRY, and YAGNI constraints

### SOLID

- **Single responsibility:** Python owns biological resolution; the adapter
  serializes; the controller coordinates; the view renders; persistence codecs
  encode/decode.
- **Open/closed:** add finite orientation/recommendation enums and typed fields
  rather than conditional flags spread across callers.
- **Liskov substitution:** CLI, typed Python, Session replay, and Web all consume
  plans obeying the same validation invariants.
- **Interface segregation:** the Web receives the review projection it needs;
  rendering receives only a fully resolved plan.
- **Dependency inversion:** UI actions depend on the shared resolver boundary,
  never directly on orthogroup internals or SVG geometry.

### KISS

- one Web action, one palette, one draft, one final batch validation, one typed
  plan, and one history path;
- reuse the current Worker operation and palette/canvas infrastructure;
- prefer an explicit small state machine over implicit button-label logic.

### DRY

- do not implement recommendation or strand logic in JavaScript;
- do not maintain separate popup/drawer alignment handlers;
- do not keep schema 1 and schema 2 current writers in parallel when schema 1
  remains unreleased;
- use the same draft mutation for radio controls and canvas candidate clicks.

### YAGNI

Do not add scored ranking, confidence thresholds, support-count ranking,
multi-hop inference, anchor TSV input, Collinear alignment controls, synteny
propagation, a new Worker, a second history stack, a framework, a build step, or
a generalized recommendation engine. These are outside Issue #586.

## 7. Files and likely change points

The exact diff should remain as small as repository ownership permits. Inspect
before editing; do not touch a listed file merely because it appears here.

### Product authority

- `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`
- `docs/internal/PRODUCT_IMPACT_RATCHET.md`
- `tools/web-product-decisions.json`

### Python domain, rendering, and adapter

- `gbdraw/layout/similarity_alignment.py`
- `gbdraw/api/record_planning.py`
- `gbdraw/diagrams/linear/assemble.py`
- `gbdraw/web_support/similarity_alignment.py`
- `gbdraw/api/__init__.py`
- `gbdraw/session_request_codec.py`
- `gbdraw/api/session_compat.py` only where released compatibility actually
  requires it
- `gbdraw/linear.py` and CLI documentation for strict ambiguity behavior

### Web controller, surfaces, and persistence

- `gbdraw/web/js/app/similarity-alignment.js`
- `gbdraw/web/js/app/right-drawer.js`
- the existing feature-popup action owner under
  `gbdraw/web/js/app/feature-editor/`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/index.html`

### Tests and generated owners

- `tests/test_similarity_alignment.py`
- `tests/test_similarity_alignment_rendering.py`
- `tests/test_similarity_alignment_web_adapter.py`
- `tests/test_session_request_codec.py`
- focused Session/API compatibility tests
- `tests/web/similarity-alignment-actions.test.mjs`
- `tests/web/similarity-alignment-ui.playwright.spec.js`
- `tests/web/right-drawer.test.mjs`
- relevant Web Session/history/worker/architecture contract tests
- `tools/refresh_gallery_sessions.py`
- `tools/prepare_interactive_gallery_assets.py`
- generator-owned Gallery Session artifacts only when regenerated by their owner

## 8. Session sequence

| Session | File | Outcome | Runtime gate |
| --- | --- | --- | --- |
| 00 | `SESSION_00_AUTHORITY_SUPERSESSION.md` | Product Contract revision with four exact rev-2 outcomes | Must merge to `origin/dev` before S01 |
| 01 | `SESSION_01_PYTHON_DOMAIN_AND_PLAN_SCHEMA.md` | Python recommendation/orientation owner and schema-2 plan | S00 merged |
| 02 | `SESSION_02_WEB_CONTROLLER_AND_DRAFT.md` | single state machine, local preselected draft, atomic Apply | S01 complete |
| 03 | `SESSION_03_SINGLE_ALIGN_REVIEW_UI.md` | one Align action, accessible full review palette, canvas parity | S02 complete |
| 04 | `SESSION_04_PERSISTENCE_AND_WORKFLOW_INTEGRATION.md` | Session/history/reset/retry, artifact regeneration, docs | S03 complete |
| 05 | `SESSION_05_ACCEPTANCE_AND_HANDOFF.md` | full acceptance, architecture/product evidence, handoff | S04 complete |

Each session prompt is an executable instruction for a fresh contributor. Do not
skip a predecessor or combine authority and dependent runtime in one commit.

## 9. Acceptance criteria

### AC-01: immediate feedback and single-flight operation

From either Web entry point, the exact clicked reference starts one resolver
operation. The action becomes visibly and accessibly busy before awaiting the
Worker. Repeated activation cannot start a duplicate job.

### AC-02: one review surface

The Web exposes only `Align…`. Every valid initial resolution opens one palette,
even when all anchors are uniquely resolved. The reference identity shown in the
palette is the exact initiating feature.

### AC-03: deterministic, disclosed defaults

Every usable target starts with a selected candidate and visible reason. The
same biological input and explicit identities produce the same recommendation
under record reorder, viewport/scroll changes, ribbon layout changes, and score
or support-count changes that do not change eligibility/direct-RBH facts.

### AC-04: independent local editing

Every displayed target can independently select a candidate or Skip and can
independently preserve or request matching direction. Draft edits, including
canvas clicks, create no Worker job. The default draft can be applied without a
mandatory click per record.

### AC-05: orientation correctness

Preserve is the default. Only known opposite displayed anchor strands plus an
explicit per-record match request reverse a whole target. Same or unknown
strands preserve orientation. Reference position/orientation and all target Y
positions remain unchanged. Text, `rev`, crop handling, and anchor centers are
correct after reversal.

### AC-06: Python ownership and atomic Apply

JavaScript submits one explicit batch. Python rejects unknown/stale identities,
ineligible anchors, incomplete record coverage, invalid policy combinations,
and inconsistent crop/display facts. No partial transform is admitted. A valid
plan is fully typed and schema 2.

### AC-07: failure, cancel, stale, and retry isolation

Cancel, validation failure, generation failure, and stale/superseded completion
do not overwrite the current diagram, active plan, last Result, or History. A
retryable Apply failure keeps the draft so it can be corrected and resubmitted.

### AC-08: lifecycle preservation

Successful Apply is one history operation. Regeneration, save/load, Reset Align,
Undo/Redo, record reorder, and active-plan repair retain their existing meanings
with requested and effective orientation preserved exactly.

### AC-09: programmatic surfaces and non-goals

Typed Python consumes only fully resolved typed plans. The CLI remains strict
and noninteractive, rejecting unresolved ambiguity with candidate guidance.
Documentation explicitly states that Collinear controls, anchor TSV, scored
inference, and multi-hop selection are unsupported.

### AC-10: accessible responsive review

The palette is keyboard operable, restores or moves focus deliberately, exposes
busy/status/error semantics, remains usable at a 390-pixel viewport, and keeps
canvas pan/zoom and candidate interaction coherent without becoming a modal
full-screen blocker.

## 10. Verification strategy

### Focused domain and persistence

```bash
python -m pytest \
  tests/test_similarity_alignment.py \
  tests/test_similarity_alignment_rendering.py \
  tests/test_similarity_alignment_web_adapter.py -v

python -m pytest \
  tests/test_session_request_codec.py \
  tests/test_session_compat.py \
  tests/test_api_session.py \
  tests/test_session_io.py -v
```

Required cases include the four recommendation reasons, explicit override and
Skip, record-order invariance, exact identity rejection, the orientation matrix,
idempotent centers, schema-2 round trip, released legacy Session fixtures, and
strict CLI ambiguity rejection.

### Focused Web

```bash
node --test \
  tests/web/similarity-alignment-actions.test.mjs \
  tests/web/right-drawer.test.mjs \
  tests/web/session-request.test.mjs \
  tests/web/session-authority.test.mjs

npx playwright test \
  tests/web/similarity-alignment-ui.playwright.spec.js \
  --project=chromium --workers=1
```

Add focused history, Worker, and Gallery specs when their owners change. If
Node's Playwright package is unavailable, use the repository's Python
Playwright route for equivalent targeted checks. If Chromium fails because of
the agent sandbox, rerun the same check with the required sandbox approval.

### Repository gates

```bash
ruff check gbdraw/
node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
python tools/update_cli_reference_help.py --check
python -m pytest tests/ -v -m "not slow"
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python -m build
git diff --check
```

Do not update reference SVGs unless a reviewed, intentional geometry change
requires it. Prepare the browser wheel when browser tests or offline packaging
need it; the wheel is generated and must not be committed.

## 11. Architecture fitness and Product Impact evidence

This is an architecture-bearing and product-visible change. Every runtime
session must follow:

- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`;
- `docs/internal/PRODUCT_IMPACT_RATCHET.md`;
- `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` after Session 00 merges.

Expected owner/path evidence:

- semantic owner remains the shared Python resolver;
- popup and drawer converge on one Web controller path;
- local draft editing adds no Worker path;
- final Apply uses the existing resolver and canonical generation/admission
  path;
- schema 1 is removed rather than retained as a parallel current path if it is
  still unreleased;
- no superseded Web actions, global mode UI, or auto-apply continuation remains.

For each session, inspect production, tests, documentation, and generated diffs
separately. Update the execution ledger below with evidence rather than merely
marking a checkbox.

## 12. Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| A recommendation looks like scientific proof | Label it as a recommendation and display the finite basis; keep replace/Skip available |
| UI and Python choose differently | Python returns candidate order, recommendation, reason, and final validation; JS never ranks |
| Same effective orientation loses user intent | Persist both requested per-record policy and effective reversal |
| Repeated action races or stale result overwrites state | one state machine, busy guard, operation token, existing Result admission |
| Review palette becomes dense | concise row hierarchy, progressive candidate details, keyboard support, 390 px browser test |
| Schema work grows compatibility debt | repeat release ancestry audit; replace unreleased schema in place; add a reader only for a proven released obligation |
| Generated Gallery artifact is edited by hand | regenerate with owning tools and review manifest/session diff |
| Scope drifts into smart/scored alignment | enforce the explicit YAGNI list and signed exclusions |

## 13. Definition of done

The work is complete only when all of the following are true:

1. the four revision-2 Product decisions are merged into `origin/dev` before
   dependent runtime code;
2. the fixed implementation branch contains that authority and all sessions;
3. Python is the sole owner of eligibility, recommendation, orientation effect,
   and final plan validation;
4. Web has one `Align…` action and one always-reviewed, preselected local draft;
5. per-record orientation is persisted as requested intent plus effective
   outcome in schema 2;
6. Apply is atomic and failure/cancel/stale/retry preserve the existing result
   and draft semantics;
7. current Session, Reset Align, Undo/Redo, regeneration, record reorder, CLI,
   and typed Python contracts pass;
8. branch-owned generated artifacts and user documentation are refreshed by
   their owners;
9. focused tests, architecture/product gates, non-slow suite, read-only output
   comparison, build, and visual browser review pass;
10. no out-of-scope inference, alternate workflow, compatibility branch, or
    framework has been added.

## 14. Execution ledger

| Session | Status | Commit | Evidence / notes |
| --- | --- | --- | --- |
| Plan | complete | initial commit for this directory | Plan and six self-contained prompts created on the fixed implementation branch |
| 00 Authority supersession | pending | — | Must merge to `origin/dev`; record merge commit |
| 01 Python domain and schema | complete | this Session 01 commit | Authority: revision-2 PD-OI-026/027/031/034 merged into `origin/dev` at `762198fc`; repaired dev staging CI passed at `ddbe68f7` (run `36033558444`). Schema audit: introducing commit `0624eb82` absent from `origin/main` first-parent history and all release tags; one schema-1 Gallery Session inventoried for Session 04, none under `tests/test_inputs`; replaced the unreleased current schema-1 writer/reader in place with schema 2. Owner/path: `gbdraw/layout/similarity_alignment.py` owns candidate eligibility, canonical recommendation, requested policy and effective orientation; `gbdraw/web_support/similarity_alignment.py` only projects JSON; existing `record_planning.py` materializes once and existing request/Session codecs serialize. Removed global mode and its current consumers, with no second plan type or JS recommendation algorithm. Verification: focused Python groups 106 + 200 passed, Session I/O 224 passed; Ruff and `git diff --check` passed; Web request, authority, Gallery migration and architecture contracts (137) passed; Web change-budget gate PASS with review required for schema. Production, test, docs and generated diffs reviewed separately; no generated Gallery edit. |
| 02 Web controller and draft | complete | this Session 02 commit | Authority: revision-2 PD-OI-026/027/031/034 is on `origin/dev@762198fc`, and Session 01 schema-2 commit `594084ba` is an ancestor. `similarity-alignment.js` owns one `idle → resolving → reviewing → applying → idle` operation, one row-per-target local draft, operation-token and committed-artifact admission, and retry/repair. Popup and drawer are exact-reference adapters; canvas and palette dispatch the same local Select/Skip action. Python remains the only recommendation/eligibility/orientation owner; the Web validates its schema-2 projection without reranking. Removed the automatic resolved Apply, separate `ready`/`ambiguous` controller paths, and Web global mode. Deterministic Worker fakes show one initial resolver job, zero jobs for local candidate/Skip/orientation edits, and one batch resolver plus one canonical generation call on Apply; duplicate starts/Apply produce no extra jobs. `runAnalysis` retains canonical Result/History admission and `svg-result-ingestion.js` retains SVG admission; stale/failure/Cancel leave prior Result and History intact. Focused Web tests: 27 passed; Chromium alignment spec: 5 passed after rebuilding the ignored browser wheel; architecture contracts: 137 passed; Web change-budget Gate PASS, Review REQUIRED for one computed busy declaration; `git diff --check` passed. Production, tests, and generated diffs reviewed separately; no tracked generated assets changed. Session 03 still owns visual layout, and Session 04 owns Gallery Session regeneration. |
| 03 single review UI | pending | — | — |
| 04 persistence/workflow/docs | pending | — | — |
| 05 final acceptance | pending | — | — |

Update only the row for the session actually executed, except Session 05 may
add final cross-session evidence. A status claim without commands/results and a
reviewed diff is not completion evidence.
