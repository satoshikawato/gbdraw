# Issue #561 Similarity Group alignment — master implementation plan

Status: implementation complete; Product authority merged; S01–S07 evidence complete

This document is the self-contained implementation plan for
[GitHub Issue #561](https://github.com/satoshikawato/gbdraw/issues/561),
"Define expected Align behavior for Similarity Group edge cases (missing,
inparalogs, strand mismatches)." Its reader is assumed to know the gbdraw
repository but not any earlier discussion, review, or planning session.

## 1. Fixed branch and delivery sequence

All implementation sessions in this plan use the branch:

```text
issue-561-similarity-alignment
```

The branch was created without an upstream from
`origin/dev@11aae136694a4433cabc68c0dae31edf77222740` on 2026-09-22.
Authority-only PR
[#570](https://github.com/satoshikawato/gbdraw/pull/570) then merged the six
Product Decisions in section 4 into `dev` at
`a9eaeadd105e0c26e46086626feaa695bdd33c94`. This fixed branch has been
rebased onto that merge commit. Candidate authority on an implementation branch
does not authorize runtime in the same candidate.

Do not create a replacement runtime branch merely because a later session
starts. At the beginning of every session:

1. fetch `origin`;
2. confirm the current branch, HEAD, upstream, and worktree status;
3. confirm that `origin/dev` contains the six accepted authority records;
4. rebase this branch onto the latest `origin/dev` before the first runtime
   commit if its base does not contain them; and
5. preserve unrelated changes and do not commit directly to `dev` or `main`.

The planning base writes Session version `43` and canonical render-request
schema `7`. Those numbers are observations, not reserved allocations. The
implementation is expected to need one new current Session version and one new
request schema because the writer changes from a legacy string to a typed plan.
The session that owns that migration must re-audit the latest base and allocate
the next available versions instead of assuming `44` and `8` remain available.

## 2. Problem and required outcome

The current Linear Similarity Group Align action stores only a group ID. The
clicked feature is resolved in the Web UI, but its identity is discarded before
generation. The Python renderer then selects a group representative or a
score-ranked member for each record. This causes five user-visible problems:

1. clicking a non-representative feature does not guarantee that feature is the
   reference;
2. a record with another valid member can be confused with a record that has no
   member;
3. inparalogs are resolved by an implementation fallback rather than explicit
   evidence or user choice;
4. position alignment and whole-record orientation are not modeled as separate
   intents; and
5. a group string does not record enough information to reproduce, inspect,
   reset, or repair an alignment.

The completed feature must treat alignment as a display operation over one
selected anchor per displayed record. It must not change Similarity Group
membership, representative status, comparison edges, or LOSATP results.

## 3. Terms

- **Reference anchor:** the exact feature selected to initiate alignment. Its
  record remains fixed.
- **Target anchor:** the selected feature in another displayed record.
- **Usable candidate:** a member of the selected group whose canonical feature
  identity resolves uniquely and whose center maps into the current cropped
  display coordinates.
- **Base record translation:** the persistent X/Y translation for a displayed
  record before an active similarity alignment is overlaid.
- **Alignment Plan:** immutable semantic intent containing the reference,
  per-record decisions, rationale, mode, and effective orientation overrides.
- **Effective record transform:** the base translation and presentation plus
  the active plan's derived X translation and optional orientation override.
- **Legacy alignment:** request schema 7 or an older supported schema containing
  `align_orthogroup_feature` / `alignOrthogroupFeature` as a string.

The atomic unit is a displayed biological record identified by stable
`recordKey`, not an uploaded file, row number, list position, or SVG element ID.

## 4. Accepted Product authority

The Product Decision Owner `satoshikawato` accepted all six decisions on
2026-09-22. They must be recorded separately so that one outcome can later be
superseded without changing the others. They are recorded in revision 10 of
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`.

| Authority ID | Concern | Selected outcome |
| --- | --- | --- |
| `PD-OI-026` | `diagram-generation.similarity-alignment.anchor-resolution` | `A / EXPLICIT_DETERMINISTIC_RESOLUTION` |
| `PD-OI-027` | `diagram-generation.similarity-alignment.transform-semantics` | `A / SEPARATE_POSITION_AND_ORIENTATION` |
| `PD-OI-028` | `diagram-generation.similarity-alignment.plan-lifecycle` | `A / PERSISTED_ACTIVE_ALIGNMENT_PLAN` |
| `PD-OI-029` | `diagram-generation.similarity-alignment.reset-and-history` | `A / IMMEDIATE_PREALIGN_BASELINE` |
| `PD-OI-030` | `diagram-generation.similarity-alignment.session-compatibility` | `A / LEGACY_READER_ONLY` |
| `PD-OI-031` | `diagram-generation.similarity-alignment.surface-scope` | `A / WEB_TYPED_CORE_STRICT_CLI` |

### 4.1 Anchor resolution

The exact clicked feature is the reference. Resolve each other displayed record
in this strict order:

1. valid explicit user selection;
2. exactly one usable candidate;
3. exactly one distinct candidate connected directly to the reference by an
   `edgeKind=rbh` edge, treating query/subject direction symmetrically;
4. explicit Select or Skip when multiple candidates remain; or
5. unchanged position and orientation when no usable candidate exists.

Representative status, score, confidence, supporting-edge count, leftmost or
central coordinates, viewport position, visible ribbons, and multi-hop paths
must not select an anchor. Multiple direct RBH candidates remain ambiguous.
Hidden features remain usable when their centers map into the crop. A feature
partly overlapping the crop is unusable when its center is outside it.

### 4.2 Position and orientation

`Align` changes only target X positions. It preserves reference transform,
target Y positions, and all orientations. `Align & orient` may reverse a target
record only when both displayed anchor strands are known and opposite. Unknown
or mixed strand applies position only. Reverse applies canonically to the whole
record while text remains readable. A persistent `rev` indicator is derived
from effective orientation relative to the source.

Alignment is absolute and idempotent. Multiple records in one Linear row are
supported as independent record units.

### 4.3 Active lifecycle

The plan survives Session round trips and ordinary Generate calls after style,
label, or canvas-size changes. Stable record reorder preserves it. Manual
record movement, manual orientation change, source replacement, crop change,
or selector change clears it with a visible reason. Stale reference blocks
regeneration until reselected or cleared. A stale target requires reselection
or explicit Skip. No candidate is silently substituted.

### 4.4 Reset and history

A new Align replaces the prior active plan after materializing the currently
effective record transforms as the new base. Therefore, clearing the new plan
restores the immediately preceding geometry without reviving the prior plan.
Apply, Reset, and manual clear are each one artifact transaction. Normal Undo
can restore the preceding artifact and plan. There is no alignment-specific
deep-history stack. Failed, canceled, superseded, or stale operations commit no
history entry and do not replace the last successful Result.

### 4.5 Session compatibility

Current writers emit only the resolved typed plan. Supported legacy Sessions
remain readable through an isolated reader-only adapter that reproduces their
historical group-ID selection. The normal resolver and current writer must not
call that compatibility path. A successfully materialized legacy alignment is
written only as the current resolved representation. Malformed or unmappable
legacy values produce an actionable error rather than silent substitution.

If a released legacy fixture can be loaded and saved before sufficient metadata
exists to materialize a resolved plan, stop the affected migration session and
prepare a narrow Product Decision Pack for that exact save continuation. Do not
silently drop the value, keep writing the legacy field, or invent a blocking
save rule without authority.

### 4.6 Initial public surfaces

- Web Similarity Groups supplies exact reference selection, `Align`,
  `Align & orient`, ambiguous Select/Skip, an operation summary, and plan
  inspection.
- The typed Python API accepts the resolved plan.
- CLI accepts an exact reference and succeeds only when all records resolve
  uniquely. Ambiguity is an actionable error; CLI does not prompt.
- New CLI and Python calls no longer accept a group ID as permission to select
  a representative implicitly.
- A new anchor TSV, Collinear-mode UI, smart/synteny propagation, and multi-hop
  inference are not part of this implementation.

## 5. Existing authority that must remain intact

These requirements are not new choices and must not be duplicated as new
Product Decisions:

- `PD-OI-016` and `OIPC-C07`: failed, canceled, superseded, and stale work does
  not replace the last successful Result or committed request.
- `PD-OI-018`: drawing-start and reverse display changes reuse compatible raw
  LOSATP evidence; display transforms are not raw-search identity.
- `product.canonical-render-request-boundary`: fresh generation and Session
  replay use one typed request path.
- `OIPC-C04`: request, execution, cache, artifact, and Session semantics agree.
- `OIPC-C05` and `OIPC-C06`: valid intent is not silently lost, and replacement
  or clearing is explicit.

## 6. Current implementation audit

The planning base has the following relevant behavior.

### 6.1 Web loses the exact click

`gbdraw/web/js/app/app-setup.js` resolves
`clickedOrthogroupDetail.currentMember`, but `alignByClickedOrthogroup()` stores
only `orthogroupId` in `selectedOrthogroupAlignmentFeature`. The group drawer in
`gbdraw/web/js/app/orthogroups.js` does the same. Both call full
`runAnalysis()` immediately, so the state is mutated before candidate
validation completes.

`gbdraw/web/js/state.js` owns the alignment as one string.
`gbdraw/web/js/services/session-request.js` serializes it under protein
comparison settings as `alignOrthogroupFeature`, coupling display layout to
analysis intent.

### 6.2 Python chooses implicit members

`gbdraw/diagrams/linear/orthogroup_alignment.py` treats a group ID as a valid
target, chooses its representative or first member, then uses representative
and score ordering per target record. `gbdraw/diagrams/linear/assemble.py`
rejects alignment when multiple records share a row.

The current CLI option `--align_orthogroup_feature` accepts a feature hash or
protein ID but is transported as the same untyped string. The typed API stores
that string in `LinearDiagramOptions.align_orthogroup_feature`.

### 6.3 Reusable foundations already exist

- `gbdraw/web/js/services/feature-identity.js` validates canonical record and
  feature aliases and uniquely resolves a clicked group member.
- `gbdraw/web_support/orthogroup_metadata.py` exposes group members and direct
  `orthologEdges`, including `edgeKind`, query/subject record indexes, and
  protein IDs.
- `gbdraw/layout/record_coordinates.py::RecordDisplayTransform` is the existing
  source/local/display coordinate authority.
- `gbdraw/api/record_planning.py` owns effective per-record reverse
  presentation and record-aligned provenance.
- `gbdraw/web/js/services/history.js::runUndoableArtifactReplacement()` owns
  atomic Result replacement and rollback.
- the existing diagram Worker supports typed helper operations and must remain
  the only Pyodide runtime.
- composition metadata already records per-primary-target X/Y deltas, but its
  array order is not a stable persisted alignment identity.

## 7. Target architecture

```text
Web popup / group drawer          CLI adapter          typed Python caller
             |                        |                        |
             +------------------------+------------------------+
                                      |
                                      v
                     shared pure candidate resolver
                  (stable identities + direct edges only)
                                      |
                      resolved / ambiguous / skipped
                                      |
                                      v
                    immutable SimilarityAlignmentPlan
                       canonical request + Session
                                      |
                     request planning materializes
                     effective orientation per record
                                      |
                                      v
             existing RecordDisplayTransform + Linear placement
                                      |
                         final record translations
                                      |
                                      v
               existing renderer / artifact transaction / Result
```

The boundaries are semantic. They do not require one class or module per box.
Create a new abstraction only when it removes at least two existing execution
paths or supplies a shared Python/Web/CLI contract.

### 7.1 Canonical model

Add a small immutable typed model, expected under
`gbdraw/layout/similarity_alignment.py`, and re-export only the public request
types required by `gbdraw.api`. Exact names may follow current repository
conventions, but the responsibilities must remain:

```text
SimilarityAlignmentPlan
  schema: 1
  mode: position | position_and_orientation
  group_id
  reference: AlignmentAnchorIdentity
  records: tuple[AlignmentRecordDecision, ...]

AlignmentAnchorIdentity
  record_key
  biological_feature_id
  source_feature_index, when needed to disambiguate a duplicate biological ID
  stable_feature_svg_id, as consistency evidence rather than a rendered ID

AlignmentRecordDecision
  record_key
  status: reference | aligned | skipped
  anchor, only for reference/aligned
  rationale: reference | user_selected | only_usable_candidate |
             unique_direct_rbh | skipped_by_user | skipped_no_candidate |
             skipped_unmappable
  effective_reverse_complement, only when the plan overrides base presentation
```

Do not persist list indexes, row numbers, rendered SVG IDs, viewport positions,
or display aliases as identity. Decoder validation rejects duplicate record
keys, duplicate decisions, identities that disagree across supplied aliases,
unknown enum values, and a reference not represented by exactly one reference
decision.

### 7.2 Minimal base translation owner

Add a Linear layout value keyed by stable `recordKey` containing finite X/Y
translations. It is the baseline owner for manual record movement and alignment
Reset. It is not a general affine-transform framework: no scale, rotation,
matrix, animation, constraints, or cross-mode abstraction is added.

Existing composition deltas are converted through one adapter:

- before Align, read the current per-record delta and materialize it into the
  keyed base translation;
- generated SVG composition starts from the request-owned translation and does
  not reapply the same legacy delta;
- after a manual record drag, materialize any active plan, clear it, and write
  the new keyed base translation; and
- Session save/load uses stable record keys, not primary-target array order.

Legend, title, length-bar, and whole-diagram composition deltas retain their
existing owners and formats.

### 7.3 Pure resolver

The resolver consumes only validated values:

- exact reference identity;
- group members indexed by stable record identity;
- current crop/display mappability and displayed strand;
- normalized direct edges; and
- optional explicit choices keyed by record.

It returns an immutable resolution result with decisions and unresolved
candidates. It does not read Vue state, show dialogs, load files, mutate record
presentation, compute SVG pixels, call LOSAT, or choose a fallback based on
render state. Web calls the same Python resolver through an operation on the
existing diagram Worker. JavaScript presents unresolved candidates but does not
reimplement ranking.

### 7.4 Transform calculation

Resolve orientation before anchor centers. Let:

```text
R = automatic reference record X
    + base reference translation X
    + reference anchor center after the effective reference transform

T = automatic target record X
    + target anchor center after the effective target transform

final target translation X = R - T
```

The aligned target's previous X translation is replaced by this absolute value;
its Y translation remains the base Y value. Skipped records retain both base
translations. The reference retains both base translations. Reapplying the
same plan produces the same values.

Canvas extents, ruler position, definitions, tracks, annotations, comparisons,
and composition metadata must use the same final record translations. Do not
translate only feature groups after layout.

### 7.5 Effective orientation

Base orientation remains owned by `RecordPresentation.reverse_complement` and
the existing record-planning transform. An active plan may supply an effective
boolean override for aligned targets. It must not create a second source
orientation flag. Request planning constructs the existing effective display
transform from the base presentation plus that override before calculating
centers or rendering records.

When a plan is materialized because of Reset replacement, new Align, or manual
editing, its effective orientation becomes the base presentation and the plan
is removed. This gives one effective orientation at every render boundary.

### 7.6 Web orchestration

Add one focused Web owner, expected as
`gbdraw/web/js/app/similarity-alignment.js`, with a top-level `create*` entry
point. It coordinates:

- exact reference construction from the popup or drawer;
- the Worker resolver operation;
- the ambiguity draft and candidate highlight;
- Apply/Cancel/Skip;
- active-plan inspection and summary;
- materialize, clear, and Reset operations; and
- lifecycle invalidation callbacks.

`app-setup.js` wires dependencies only. `orthogroups.js` continues to own group
browsing and highlighting, not alignment semantics. `session-request.js`
projects the canonical plan and translations but does not resolve candidates.
`run-analysis.js` performs one existing artifact replacement; it does not own
the resolver policy.

### 7.7 Result and inspection

After Apply, show a concise summary containing aligned, unchanged, skipped, and
reversed record counts. The active-plan inspector exposes the reference and one
row per displayed record with anchor identity and rationale. Candidate details
include feature ID, coordinates, displayed strand, representative/role status,
and direct evidence. Only ambiguous records require interaction.

The `rev` indicator is derived from the effective source-relative orientation;
it is not a separately editable or persisted boolean.

## 8. User workflows

### 8.1 Feature-popup Align

1. User selects an exact feature and chooses `Align` or `Align & orient`.
2. The controller validates the exact stable reference.
3. The Worker resolver returns automatic decisions and unresolved records.
4. If none are unresolved, the controller previews the summary and applies one
   artifact transaction.
5. Otherwise, the dialog visits only ambiguous records. Cancel changes no
   state or Result. Apply is disabled until each ambiguous record has Select or
   Skip.
6. A successful render admits the plan, request, Result, and history entry
   together.

### 8.2 Group-drawer Align

The drawer requires an exact reference selector before enabling either action.
It never submits the group ID as an implicit reference. The selected reference
can be changed before Apply.

### 8.3 Ordinary Generate

Generate validates an active plan before rendering. Style, label, and canvas
changes retain it. Resolver input comes from committed comparison evidence;
Generate does not launch LOSATP solely to satisfy the plan. Missing required
evidence produces a stale or actionable validation result while the current
Result remains.

### 8.4 Reset, Undo, and manual movement

- Reset clears the active overlay and renders the saved base translations and
  base orientations.
- Undo restores the complete previous artifact, including a previous active
  plan when present.
- a new Align materializes the old effective state as the new base, then
  installs the replacement plan;
- manual record drag or orientation edit materializes the effective state,
  clears the plan with a reason, then applies the edit; and
- stable reorder changes ordering only and keeps the key-addressed plan.

### 8.5 Session Load

Load displays the saved preview without starting the diagram Worker. The first
Generate or explicit alignment inspection validates the plan. A stale
reference blocks new generation; a stale target requires Select or Skip. The
saved preview remains the last successful Result until repair succeeds.

## 9. Request and Session migration

The current request field is an analysis-pipeline string. Replace its current
writer projection with typed display state, expected conceptually as:

```json
{
  "linearLayout": {
    "recordTranslations": [
      {"recordKey": "record-1", "x": 0.0, "y": 0.0}
    ],
    "similarityAlignment": {
      "schema": 1,
      "mode": "position",
      "groupId": "og_1",
      "reference": {},
      "records": []
    }
  }
}
```

The exact JSON nesting must follow the latest canonical request conventions.
There is one current writer and one typed decoder. Do not add a parallel
top-level Session-only copy, argv-shaped Worker message, or Web-only plan.

Retain old supported request schemas as readers. Their
`align_orthogroup_feature` value enters a named legacy adapter that invokes the
old representative/score behavior only to reproduce that historical request.
The current typed API, CLI, Web action, and current writer cannot call it.

The migration session must add one representative positive fixture proving a
released legacy contract, negative malformed fixtures, and current round-trip
fixtures. Branch-only intermediate schema/session versions are not supported.

## 10. Public API and CLI

The typed Python entry accepts `SimilarityAlignmentPlan | None`. Remove the
normal current-path `str | None` owner in the same change. Keep a private legacy
decode representation only where an old schema requires it.

Retain `--align_orthogroup_feature` as the initial CLI spelling unless current
CLI conventions give a clearer non-breaking alias, but change its documented
meaning to an exact stable feature/protein reference. A group ID that maps to
multiple candidates is not accepted. CLI uses the shared resolver and:

- proceeds when every participating record is automatically resolved or has no
  usable candidate;
- reports records with ambiguous candidates and suggests exact identifiers;
- performs no interactive prompt; and
- performs no new LOSATP search merely to resolve a saved plan.

Do not add an anchor-selection TSV or a second plan-building API.

## 11. SOLID, KISS, DRY, and YAGNI constraints

### SOLID

- Resolver owns candidate decisions, not UI or geometry.
- Alignment Plan owns semantic state, not file loading or history.
- Existing record transform owns coordinate/orientation projection.
- Existing artifact transaction owns commit, rollback, Undo, and Redo.
- Session codec owns compatibility and current serialization.
- Surface adapters depend on the typed resolver contract and do not implement
  their own policies.

### KISS

- one plan schema;
- one strict priority list;
- one Worker;
- one current request writer;
- one compatibility adapter for released legacy schemas; and
- one active alignment at a time.

### DRY

- Web, CLI, and Python call the same resolver;
- all geometry uses the same final per-record translations;
- `rev` is derived from effective orientation;
- existing feature identity, request, record transform, history, and failure
  admission paths are reused; and
- superseded string writers and implicit current resolvers are removed in the
  same implementation.

### YAGNI

Do not add smart alignment, synteny propagation, multi-hop inference, an anchor
TSV, Collinear UI, multiple simultaneous plans, a constraint solver, a generic
affine transform framework, another Worker, or a new decision store.

## 12. Implementation sessions

The following sessions are sequential unless the dependency column explicitly
permits otherwise. Each has a separate, self-contained instruction file.

| Session | Scope | Depends on | Instruction file |
| --- | --- | --- | --- |
| S01 | typed domain model and pure resolver | merged authority | `ISSUE_561_SIMILARITY_ALIGNMENT_S01_DOMAIN_RESOLVER_INSTRUCTION_PROMPT_2026-09-22.md` |
| S02 | canonical request, base translations, Session migration | S01 | `ISSUE_561_SIMILARITY_ALIGNMENT_S02_REQUEST_SESSION_INSTRUCTION_PROMPT_2026-09-22.md` |
| S03 | Python render transforms, typed API, strict CLI | S01–S02 | `ISSUE_561_SIMILARITY_ALIGNMENT_S03_PYTHON_RENDER_CLI_INSTRUCTION_PROMPT_2026-09-22.md` |
| S04 | Web Worker resolver operation and alignment controller | S01–S03 | `ISSUE_561_SIMILARITY_ALIGNMENT_S04_WEB_CONTROLLER_INSTRUCTION_PROMPT_2026-09-22.md` |
| S05 | popup/drawer/dialog/inspector UI and accessibility | S04 | `ISSUE_561_SIMILARITY_ALIGNMENT_S05_WEB_UI_INSTRUCTION_PROMPT_2026-09-22.md` |
| S06 | lifecycle invalidation, Reset, history, legacy composition bridge | S02, S04–S05 | `ISSUE_561_SIMILARITY_ALIGNMENT_S06_LIFECYCLE_HISTORY_INSTRUCTION_PROMPT_2026-09-22.md` |
| S07 | end-to-end acceptance, documentation, architecture gates | S01–S06 | `ISSUE_561_SIMILARITY_ALIGNMENT_S07_ACCEPTANCE_INSTRUCTION_PROMPT_2026-09-22.md` |

Every session updates section 16 with its base/head SHA, files changed, tests,
manual evidence, remaining risks, and the exact next permitted session.

## 13. Verification matrix

### 13.1 Pure resolver

- exact non-representative reference;
- zero, one, and multiple usable members;
- explicit choice overriding one/member/RBH rules;
- one direct RBH versus multiple RBHs;
- direction-normalized RBH;
- non-RBH and multi-hop edges ignored;
- hidden but mappable member;
- partially cropped member with center outside;
- duplicate or conflicting stable identity rejected;
- deterministic output independent of input ordering.

### 13.2 Geometry and orientation

- Align changes X only and preserves every Y;
- reference transform is byte-for-byte equivalent in the plan result;
- positive/negative, negative/positive, same-strand, and unknown-strand matrix;
- reverse occurs before target-center calculation;
- repeated Generate and repeated application do not accumulate offsets;
- skipped and missing records retain base X/Y/orientation;
- multiple records in one row;
- normalize-length and center-alignment combinations;
- definitions, tracks, ruler, labels, annotations, and comparison edges follow
  the final record transform;
- canvas expands in either X direction without clipping.

### 13.3 State and compatibility

- current request rejects legacy group-only input;
- supported legacy request reproduces historical output through its adapter;
- current writer emits no legacy alignment string;
- current Session round trip retains plan, rationale, base translations, and
  effective orientation;
- Load does not construct the Worker until validation/generation is requested;
- stale reference, stale target, malformed plan, and unknown enums;
- source/crop/selector/manual orientation/manual drag clear the plan;
- stable reorder preserves it;
- new Align, Reset, Undo, Redo, Cancel, failed render, canceled render, and
  stale completion have the accepted transitions;
- no LOSATP job count increase for Align or ordinary plan regeneration.

### 13.4 Surfaces and accessibility

- popup exact reference;
- drawer reference selector;
- automatic no-dialog path;
- ambiguous keyboard and pointer selection;
- Skip and Cancel;
- focus trap, initial focus, Escape behavior, focus return, accessible names,
  and live summary;
- active-plan inspector and persistent `rev` indicator;
- CLI unique success and ambiguous actionable error;
- typed Python request success and validation errors.

### 13.5 Required gates

Run focused tests after each session and the following combined gates in S07:

```bash
ruff check gbdraw/
python -m pytest tests/ -v -m "not slow"
node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
python -m build
```

Use Node Playwright when `@playwright/test` is available. Otherwise run an
equivalent targeted Python Playwright check. If Chromium is blocked by the
agent sandbox, rerun the same check with required escalation. Do not rewrite
tracked reference SVGs unless reviewed geometry changes intentionally require
it; use the documented update command and inspect each diff.

## 14. Architecture and Product Impact review

This is architecture-bearing because it replaces a string semantic owner,
adds typed persisted layout state, and adds a Worker helper operation. The
expected result is non-increasing ownership and path count:

- remove the current string writer and normal implicit resolver;
- retain one isolated legacy reader;
- use the existing request, Worker, record transform, and history paths;
- avoid a JS resolver parallel to Python; and
- avoid a second persisted representation.

Each implementation session records concise owner/path evidence under the
Architecture Fitness Function Ratchet. If the proposed work requires a second
normal runtime path, a new semantic owner, or a compatibility writer, stop and
perform the exception process instead of rationalizing the expansion.

The mapped `product.canonical-render-request-boundary` concern must retain both
effects: current state at Generate submission and equivalent regeneration after
a Session round trip. Matching an option name without both requirements is not
sufficient.

## 15. Completion criteria

Implementation is complete only when all of the following are true:

1. the six authority records are merged into the runtime base;
2. all S01–S07 scopes and evidence are recorded;
3. current Web, CLI, and Python execution use the same typed plan and resolver;
4. the legacy string is reachable only from supported old-schema readers;
5. exact reference, missing member, inparalog, direct RBH, strand, same-row,
   Reset, Undo, Session, and stale cases pass;
6. Align performs no LOSATP or group inference work;
7. failed/canceled/stale work retains the last successful Result;
8. accessibility and real-browser workflows pass;
9. production, tests, docs, and generated diffs receive separate review;
10. Architecture and Product Impact gates pass without an unapproved exception;
    and
11. public documentation states the exact Web, CLI, Python, deferred Collinear,
    and legacy Session behavior.

## 16. Execution ledger

Append one row and a short evidence note at the end of each session. Do not use
private conversation as evidence.

| Session | Status | Base / head | Evidence | Remaining work |
| --- | --- | --- | --- | --- |
| Planning | complete | base `11aae136694a4433cabc68c0dae31edf77222740` | Issue updated `2026-09-21T05:25:23Z`; code and policy audit on 2026-09-22 | S01 |
| Authority | complete | `origin/dev@a9eaeadd105e0c26e46086626feaa695bdd33c94` | PR #570 merged six separate Product Decision receipts as `PD-OI-026`–`PD-OI-031`; branch rebased onto the merge | S01 may start |
| S01 | complete | base `11e49d32accd0b00753ec7df9bcb6fb0c71d62f5`; head is the commit containing this ledger entry | Typed immutable plan/candidate/choice/edge/outcome models and one pure deterministic resolver; focused and existing identity/alignment tests pass | S02 may start; request/Session/runtime wiring remains intentionally absent |
| S02 | complete | base `0624eb821326f0022238b390a1b56f454c51151c`; head is the working-tree candidate containing this ledger entry | Canonical request schema 8 and Session 44 own the typed plan and finite record-keyed X/Y base translations; released v39/v40 legacy evidence materializes through one reader-only adapter; Python, Web, browser, recipe, documentation, and architecture checks pass | S03 may start; current renderer geometry, strict CLI/API resolution, and Web interaction remain intentionally absent |
| S03 | complete | base `9e319a8551e4c9df58cd79277e26fbd399ea2dcf`; head is the working-tree candidate containing this ledger entry | Typed plans now resolve effective orientation and transformed anchor centers during request planning, then drive one final per-record translation through the existing Linear geometry path; strict CLI resolution uses the shared resolver and reuses one completed analysis | S04 may start; Web Worker/controller interaction remains intentionally absent |
| S04 | complete | base `be9d418b363bb43ce062bf9ddbabe081ae90324d`; head is the working-tree candidate containing this ledger entry | One Web controller sends validated current facts through one typed operation on the existing lazy diagram Worker to the S01 resolver, keeps immutable drafts ephemeral, and applies a completed plan through one existing undoable artifact replacement | S05 may start; complete UI/accessibility and S06 lifecycle/Reset behavior remain intentionally absent |
| S05 | complete | base `e3e214fa2f9fbae6804b0c4a44bd41150434d9f0`; head is the working-tree candidate containing this ledger entry | Popup and drawer expose exact-reference `Align` / `Align & orient` journeys through the S04 controller; resolver-only ambiguity UI, preview, atomic Apply, live summary, read-only inspector, and accessible desktop/narrow behavior pass focused and real-browser checks | S06 may start; lifecycle invalidation, Reset/materialization, complete Undo/Redo lifecycle, and the legacy composition bridge remain intentionally deferred |
| S06 | complete | base `fe0e608dccba0676ee5fdf4c730e94a394b79ea8`; head is the working-tree candidate containing this ledger entry | The active plan is one overlay on record-keyed base translations and `RecordPresentation`; sequential Align, Reset, semantic invalidation, stable reorder, manual orientation/drag, complete History restoration, stale Generate repair, current/legacy composition materialization, and lazy Session round trips pass focused and real-browser checks | S07 may start; final acceptance, public documentation, and release work remain intentionally deferred |
| S07 | complete | base `568b072902cb7b0d18aeeabc0dd193f8a97418da`; head is the commit containing this ledger entry | All 11 Issue criteria mapped to accepted authority, production owner, and verification; public docs/examples, browser journeys, required gates, architecture review, and full branch diff pass | None in Issue #561 scope; PR, merge, release, tag, and deployment remain separate |

### S01 evidence — 2026-09-22

- Authority/base: fetched `origin`; both `origin/dev` and this branch contain
  authority merge `a9eaeadd105e0c26e46086626feaa695bdd33c94` with
  `PD-OI-026`–`PD-OI-031`. The S01 commit containing this ledger entry is based
  on `11e49d32accd0b00753ec7df9bcb6fb0c71d62f5`.
- Files: added `gbdraw/layout/similarity_alignment.py` and
  `tests/test_similarity_alignment.py`; updated this ledger only. No renderer,
  Web, request/Session schema, CLI, public documentation, or generated artifact
  changed.
- Verification:
  - `ruff check gbdraw/layout/similarity_alignment.py tests/test_similarity_alignment.py`
    — passed.
  - `python -m pytest tests/test_similarity_alignment.py -v` — 30 passed.
  - `python -m pytest tests/test_protein_colinearity.py -v -k 'orthogroup_alignment' tests/test_record_display_comparisons.py tests/test_linear_duplicate_record_ids.py::test_orthogroup_label_sets_keep_duplicate_ids_separate_by_record_index`
    — 8 passed, 304 deselected.
  - Explicit legacy projected-center and duplicate-record label-identity nodes
    from `tests/test_record_display_comparisons.py` and
    `tests/test_linear_duplicate_record_ids.py` — 2 passed.
  - Six focused `tests/test_web_feature_catalog.py` stable-identity tests — 10
    passed, including their parameterized cases.
- Manual/diff review: production and test additions were reviewed separately.
  The resolver preserves the supplied canonical record order; permutations of
  candidates and edges produce the same result. Its import boundary is checked
  to exclude analysis/LOSAT, diagrams, renderer, Web/UI, BioPython, and pandas.
- Architecture evidence: before S01, new deterministic Similarity Alignment
  resolution had no typed semantic owner or runtime path; after S01,
  `gbdraw/layout/similarity_alignment.py::resolve_similarity_alignment` is its
  single owner and future shared entry. Existing
  `gbdraw/diagrams/linear/orthogroup_alignment.py` remains unchanged solely as
  the current runtime and future bounded legacy behavior required by S01; the
  new types neither import nor call it. No current production entry is added,
  duplicated, or redirected, and no compatibility path is added in S01, so the
  ordinary review is non-increasing for `OE`, `PE`, and `CB`; no Architecture
  exception applies.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; the model and resolver encode
  `PD-OI-026`–`PD-OI-031` without selecting a new Product outcome. S01 has no
  reachable user-visible runtime delta.
- Remaining risks: adapters must still derive unique canonical identities,
  crop-mappable centers, displayed strands, and canonical edge endpoints from
  committed evidence. S02 must add the canonical request/base-translation and
  released-Session migration boundary without creating a second writer or
  allowing current requests to enter the legacy resolver.
- Next permitted session: S02 only.

### S02 evidence — 2026-09-22

- Authority/base: fetched `origin/issue-561-similarity-alignment` and
  `origin/dev`; the worktree started clean with local and remote branch heads at
  S01 commit `0624eb821326f0022238b390a1b56f454c51151c`. S01 was complete and
  S02 was the first pending ledger row.
- Current ownership: canonical request schema 8 stores one
  `SimilarityAlignmentPlan` and one finite X/Y base translation for each stable
  Linear `recordKey`; Session 44 projects that same request. Current Python and
  Web writers reject or omit `align_orthogroup_feature`,
  `alignOrthogroupFeature`, and the former Session-only selected-alignment
  string. `RecordPresentation` remains the base-orientation owner. Renderer
  application, strict CLI/API resolution, and Web interaction are deferred to
  their planned later sessions.
- Legacy boundary: schemas 1, 2, 5, 6, and 7 decode the old string only into a
  private representation. The named Python and Web reader-only adapters use
  saved stable feature/orthogroup metadata to produce schema-1 typed plans and
  zero base translations before a current save; malformed, conflicting, or
  ambiguous metadata fails closed. The historical renderer sees the string
  only at its compatibility materialization boundary.
- Released evidence: the schema-5/session-40 fixture has SHA-256
  `4eb045e070d41243f29d5486ab01be08cf18e84c4ea9a0c28b852f5862e6321d`
  and first-parent-main witness
  `10d3a3d28b3c9faa42db01c8bd0bd36b9be8433c`; its catalog schema 3 has
  stable biological IDs and source feature indexes for all five records. The
  save-before-materialization audit also covers released session 39 fixture
  SHA-256
  `e2e8296807649b4da8af38fd1f13f0af211543b36cb03c8347691ff9583e0093`
  at first-parent-main witness
  `8228ffab272d6ea2a0728ae0e1d925424431b21d`; its five representatives have
  record indexes, source feature indexes, and stable SVG identities. No valid
  released legacy edge lacks the metadata needed for a current write, so the
  S02 Product Decision Pack condition was not triggered.
- Verification:
  - the fast Python suite, excluding separately verified environment-specific
    browser/exact-replay and baseline reference-output cases, completed with
    6,143 passed, 1 skipped, and 11 deselected;
  - focused request/Session/compatibility and recipe/reproduction checks passed,
    including 39 Session compatibility tests, 34 recipe/reproduction tests, and
    an isolated exact-replay test;
  - all changed Node unit suites passed, and the canonical Session CLI browser
    integration passed all four cases;
  - a real Chromium/Pyodide flow loaded the released legacy Gallery Session
    without constructing the Worker, saved schema 44/request 8 with the typed
    five-record plan and no legacy fields, loaded it fresh without constructing
    the Worker, and generated successfully through one Worker;
  - `node tests/web/architecture-contracts.test.mjs` — 137 passed;
  - `node tools/check-web-change-budget.mjs` — Gate PASS, Review REQUIRED,
    zero blocking violations, zero import cycles, and no privileged, dependency,
    vendor, binary, or guard violations; and
  - focused Ruff checks and `git diff --check` passed. The unchanged S01 base
    reproduces the same two label-binding reference-output mismatches, and its
    existing `tests/test_session_io.py` E701 findings are unchanged.
- Diff review: production, compatibility fixtures/generated Session examples,
  tests, and documentation were reviewed separately. The generated examples
  have only the Session 43/request 7 to Session 44/request 8 semantic change
  (plus normal regenerated timestamps); public SVG output is unchanged. The one
  added Web module is the isolated legacy reader;
  the registered current render-request entry and semantic owner remain single,
  and no second Worker or runtime render path was added. This is ordinary
  architecture review with no exception.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; S02 serializes
  `PD-OI-026`–`PD-OI-031` and the documented released compatibility result
  without selecting a new Product outcome.
- Remaining risks: S02 deliberately does not apply the plan in current renderer
  geometry, replace the normal CLI/API string path, or add Web Align lifecycle
  behavior. Those are S03, S04, and S06 responsibilities respectively.
- Next permitted session: S03 only.

### S03 evidence — 2026-09-23

- Authority/base: fetched `origin/issue-561-similarity-alignment` and
  `origin/dev`; the existing worktree started clean with local and remote heads
  at S02 commit `9e319a8551e4c9df58cd79277e26fbd399ea2dcf`. S01 and S02 were complete and
  S03 was the first pending ledger row.
- Orientation and coordinate ownership: `RecordPresentation` remains the one
  base-orientation owner. `materialize_similarity_alignment_display` derives
  the effective boolean once, applies any whole-record reverse once, and then
  projects every selected source-feature center through the existing
  `RecordDisplayTransform`. The request keeps its base records, so repeated
  planning and rendering do not accumulate orientation or translation state.
- Geometry path: `_final_record_translations` implements the documented
  reference-world-X/target-post-transform formula. It preserves reference and
  skipped base X/Y plus every record's Y, replaces only aligned-target X, and
  feeds the resulting placements to records, definitions, tracks, per-record
  rulers, labels, annotations, comparisons, collision/content bounds, length
  bar extents, and composition/track metadata. Same-row records are handled
  independently, and final negative/positive, centered, normalized, and
  same-row bounds determine the canvas without an SVG-only correction path.
- Resolution and public surfaces: the normal renderer no longer collects,
  ranks, or selects alignment members. The S01 shared resolver is the sole
  candidate-decision owner; the CLI adapter only constructs canonical
  candidates/evidence and requires an exact unique feature/protein reference.
  Group IDs are rejected, missing records remain skipped, ambiguity names the
  record and exact candidate IDs, and no prompt is used. The public typed
  request and beginner Python adapter accept `SimilarityAlignmentPlan | None`;
  the former public string field was removed and invalid strings fail during
  request construction.
- Analysis and compatibility: the strict CLI performs the requested
  orthogroup analysis once, resolves the plan from that result, and reuses the
  completed raw/derived artifacts for the final typed render. A supplied or
  active plan starts no LOSATP work solely for resolution. Released old-schema
  strings still enter only the S02 reader adapter, are promoted to a typed
  plan, and then use the same current render path; rendering the released v40
  fixture, saving its current typed sidecar, and rerendering that sidecar
  produces byte-identical SVG.
- Verification:
  - alignment resolver/render, released compatibility, public Python adapter,
    and public-contract tests — 119 passed;
  - Linear track, track-slot, multi-record layout/comparison, and display
    consumer tests — 255 passed;
  - protein/orthogroup and record-display comparison tests — 302 passed,
    1 skipped;
  - typed API/request-render, request codec, API Session, and library forwarding
    tests — 387 passed;
  - Session I/O and focused CLI suites — 326 passed; the only failure is the
    unchanged H-CLI-07 published-SVG freshness mismatch caused exclusively by
    S02's already-recorded additive label-binding metadata;
  - `tests/test_output_comparison.py::TestOutputComparison` — 14 passed and the
    same two S02-base circular label-binding metadata comparisons failed;
    no reference output was regenerated or modified;
  - `ruff check gbdraw/` and focused changed-test Ruff checks — passed;
  - `node tests/web/architecture-contracts.test.mjs` — 137 passed; and
  - `node tools/check-web-change-budget.mjs` — Gate PASS, Review CLEAR, with no
    blocking violation, import cycle, privileged change, dependency change, or
    Web production delta.
- Diff review: production, tests/contract fixture, generated artifacts, and
  documentation were reviewed separately. The public-contract fixture changes
  only for the typed `LinearComparisonOptions` field and exact-ID CLI help.
  The H-CLI-07 recipe was regenerated from a clean directory and its record,
  definition, comparison, ruler, bounds, and composition geometry matched the
  published artifact; the unrelated label-binding-only output was discarded.
  No tracked reference output or other generated artifact changed.
- Architecture evidence: before S03, the current Python runtime still had the
  renderer-owned implicit representative/score selection and alignment-only
  offset/canvas path. After S03, the canonical path is typed request planning
  -> effective `RecordDisplayTransform` -> shared final Linear placement ->
  existing geometry consumers. The superseded renderer selection, ranking,
  offset, canvas-extents, and same-row rejection paths were removed. The S02
  released-schema adapter remains the one bounded compatibility path and now
  converges into the current typed path. No parallel resolver, plan builder,
  orientation owner, render path, or compatibility namespace was added, so
  this is an ordinary non-increasing architecture change with no exception.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; S03 realizes
  `PD-OI-026`–`PD-OI-031`, `PD-OI-016`, `PD-OI-018`, and the canonical-request
  contracts without selecting a materially different Product outcome. No
  Product Decision Pack was required.
- Remaining scope: Web Worker resolver operation and the alignment controller
  are intentionally deferred. S04 may start; S05–S07 must remain pending.
- Next permitted session: S04 only.

### S04 evidence — 2026-09-23

- Authority/base: fetched `origin/issue-561-similarity-alignment` and
  `origin/dev`; the existing worktree started clean with local and remote heads
  at S03 commit `be9d418b363bb43ce062bf9ddbabe081ae90324d`. The accepted
  `PD-OI-026`–`PD-OI-031`, `PD-OI-016`, and `PD-OI-018` outcomes were already
  present, S01–S03 were complete, and S04 was the first pending ledger row.
- Worker and resolver path: `createSimilarityAlignmentActions` is the one new
  alignment orchestration owner. The already-authorized `app-setup.js` Worker
  client injects one typed `resolveSimilarityAlignment` operation into it; the
  operation reuses the existing lazy diagram Worker and calls the S01 Python
  resolver through a strict JSON adapter. No Worker, JavaScript resolver,
  argv-shaped bridge, LOSATP dispatch, or group-inference path was added.
- Identity and draft behavior: popup alignment uses the exact current
  record/biological/source/stable feature identity and rejects non-unique
  matches. Drawer calls without an exact reference are rejected before the
  helper. Only current group members, direct group edges, record/crop/display
  facts, and explicit Select/Skip choices cross the helper boundary. Strict
  response validation rejects unknown fields, enums, combinations, coverage,
  references, and group-member identities. Ambiguous state is deep-frozen and
  ephemeral; Cancel, helper failure, stale completion, and supersession cannot
  update canonical state or the current Result.
- Apply and history: a completed plan plus complete record-keyed base
  translations is injected into the existing Generate candidate without first
  mutating state. The generated artifact owner set now carries the plan and
  translations, so successful Result admission commits plan, canonical request,
  Result, and history together. Render cancellation, error, stale work, and
  late canonical admission failure restore the prior artifact. Compatible
  committed protein evidence is reused; the focused integration test proves no
  additional LOSATP execution or group inference occurs for alignment Apply.
- Verification:
  - focused controller, Worker lifecycle/protocol, orthogroup identity, Session
    request, and run-analysis suites — 6 files passed, including 10 controller
    behavior cases;
  - shared resolver, Web adapter, and embedded Python helper tests — 58 passed;
  - `ruff check gbdraw/` and `git diff --check` — passed;
  - `node tests/web/architecture-contracts.test.mjs` — 137 passed; and
  - `node tools/check-web-change-budget.mjs` — Gate PASS, Review REQUIRED for
    the intended new module/reactive draft and net additions, with zero blocking
    violations, privileged expansions, import cycles, dependency/vendor/binary
    changes, or guard changes.
- Architecture evidence: before S04, popup/drawer actions wrote the legacy
  group string and called Generate directly. After S04, their only active Align
  route is exact identity -> one controller -> injected existing Worker helper
  -> shared Python resolver -> one existing generated-artifact replacement.
  Group browsing, resolver ranking, canonical request encoding, rendering, and
  generic history retain their existing owners. The new controller is the one
  required semantic owner for the new Web interaction; no superseded Align,
  resolver, Worker, render, or commit path remains. This is ordinary review
  with no architecture exception.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; S04 realizes the accepted
  exact-reference, ambiguity, atomic Apply, cancellation, and evidence-reuse
  outcomes without selecting a new Product outcome. No Product Decision Pack
  was required.
- Remaining scope: complete popup/drawer controls, the ambiguity dialog,
  candidate preview, plan inspector, summary, and accessibility are deferred to
  S05. Lifecycle invalidation, Reset/materialization, and the legacy composition
  bridge remain deferred to S06. S05 may start; S06–S07 must remain pending.
- Next permitted session: S05 only.

### S05 evidence — 2026-09-23

- Authority/base: fetched `origin/issue-561-similarity-alignment` and
  `origin/dev`; the existing worktree started clean with local and remote heads
  at S04 commit `e3e214fa2f9fbae6804b0c4a44bd41150434d9f0`. Both heads
  contained S04, S01–S04 were complete, and S05 was the first pending ledger
  row.
- UI journey: the feature popup now exposes separate exact-clicked-feature
  `Align` and `Align & orient` actions with explicit orientation descriptions.
  The drawer has no group-only alignment route: one exact record/feature must
  be selected before either action is enabled. Both entry points call the S04
  controller; `index.html` owns no resolver, decision, plan, persistence,
  render, or history behavior.
- Ambiguity and preview: only resolver-reported ambiguous records enter one
  modal. Native radio choices require Select or Skip for every such record and
  show exact feature ID, coordinates, displayed strand, representative/role,
  and direct evidence without score ranking. Hover/focus delegates exact
  identity preview to the existing SVG highlight owner and restores the prior
  feature or match highlight on cleanup. Auto-resolved operations open no
  modal; Cancel and Escape leave canonical state, Result, and history intact.
- Apply and inspection: one successful Apply remains one existing generated
  artifact/history transaction. A polite live region reports aligned,
  unchanged, explicitly skipped, no-candidate, and reversed counts. The
  read-only active-plan inspector shows the exact reference, each record's
  anchor or Skip, and saved rationale; `rev` is derived from effective
  source-relative orientation without a persisted indicator field. S06 Reset
  and lifecycle behavior was not implemented.
- Accessibility/browser verification: Python Playwright 1.61.0 drove the real
  local Worker UI because Node `@playwright/test` was not installed. The
  representative inparalog fixture passed popup and drawer entry, exact
  selector gating, modal label/description, disabled reason association,
  focus entry/trap/return, native Select/Skip, Escape no-op, candidate preview
  cleanup, one-undo Apply, live summary, and inspector checks at 1600x1000 and
  720x740. Disposable screenshots were reviewed: desktop and narrow dialogs
  kept all candidate, Cancel, and Apply controls reachable; the narrow drawer
  actions were raised above preview controls; the final summary and plan rows
  remained legible. No public screenshot or reference output changed.
- Verification:
  - `node tests/web/similarity-alignment-actions.test.mjs` — 16 passed;
  - focused Worker startup/protocol, run-analysis, History input/core,
    orthogroup identity, and right-drawer suites — passed;
  - shared resolver, rendering, and Web adapter — 59 passed;
  - targeted Python Playwright real-browser acceptance — passed;
  - `ruff check gbdraw/`, JavaScript syntax checks, and `git diff --check` —
    passed;
  - `node tests/web/architecture-contracts.test.mjs` — 137 passed; and
  - `node tools/check-web-change-budget.mjs` — Gate PASS, Review REQUIRED for
    five intended controller presentation refs/computed values and production
    net additions, with zero blocking violations, privileged expansions,
    import cycles, dependency/vendor/binary changes, or guard changes.
- Architecture evidence: the S04 controller remains the sole alignment
  orchestration owner. Popup/drawer/modal markup projects its state and actions;
  the existing diagram Worker client, shared Python resolver, canonical request
  codec, renderer, Result admission, generic history, and SVG highlighting
  retain their existing responsibilities. The group-only drawer route was
  removed in the same change. No alternate Worker, resolver, plan builder,
  persistence path, history stack, framework, or build step was added.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; S05 realizes the already
  accepted exact-reference, distinct-mode, explicit ambiguity, cancellation,
  summary, inspector, and accessibility outcomes. No materially different
  user-visible choice or Product Decision Pack was required.
- Remaining scope: source/crop/selector/orientation/drag invalidation,
  sequential materialization, complete Reset baseline behavior, Undo/Redo
  lifecycle completion, and the legacy composition delta bridge remain S06.
  S06 may start; S07 remains pending.
- Next permitted session: S06 only.

### S06 evidence — 2026-09-23

- Authority/base: fetched `origin/issue-561-similarity-alignment` and
  `origin/dev`; the existing worktree started clean with local and remote heads
  at S05 commit `fe0e608dccba0676ee5fdf4c730e94a394b79ea8`. Both heads
  contained S05, S01–S05 were complete, and S06 was the first pending ledger
  row.
- Lifecycle and history: the S04 controller remains the sole alignment
  orchestration owner. An active plan is one overlay on record-keyed base X/Y
  translations and base `RecordPresentation`. A second Align materializes the
  first plan's effective translations and orientations before installing the
  replacement in one generated-artifact transaction. Reset renders the
  immediate pre-align base without reviving an older plan. Manual orientation
  and record drag materialize once, clear with an accessible reason, and stay
  within the existing input/diagram History transaction. Generated artifact
  snapshots restore plan, base translations, base orientations, and Result;
  failed, canceled, superseded, stale, and no-op paths add no entry.
- Invalidation and regeneration: source replacement, crop, selector, and source
  type changes call explicit controller callbacks from their semantic owners.
  Stable reorder remaps plan decisions and base translations by `recordKey`.
  Ordinary Generate revalidates the plan through the shared Python resolver
  without LOSATP; stale reference offers Reselect/Clear, stale targets require
  Select/Skip, and repair retains the last successful Result and preview.
- Composition and compatibility: current Linear record groups expose final
  renderer-owned X/Y translations with stable `recordKey`. One composition
  adapter combines those values with per-record user deltas and leaves legend,
  title, length-bar, and whole-diagram placement under their existing owners.
  Released pre-recordKey v40 output binds once through the resolved plan's
  stable feature identity and renderer-authored transform; DOM order, row,
  rendered record ID, and array position are not persisted identity. The
  released five-record fixture materialized exact X values after DOM reorder,
  rejected an unmappable anchor explicitly, and constructed no Worker.
- Session and browser verification: current Session Save/fresh Load preserved
  base translations, active plan, rationale, effective orientation, Result,
  and Reset behavior; Load-only preview constructed no Worker. Python
  Playwright 1.61 drove the real UI because Node `@playwright/test` was absent.
  The browser flow covered Align A -> Align B -> Reset B -> Undo -> Redo,
  manual orientation, source/crop/selector reasons, stale reference/target
  repair, record drag, the following non-doubling Generate, and no page errors.
  The disposable 1280x900 screenshot kept the lifecycle notice, both records,
  preview controls, and editor content legible without clipping or overlap; no
  tracked public visual was added.
- Verification:
  - focused resolver, rendering, Web adapter, request codec, Session
    compatibility, and API Session tests — 251 passed;
  - controller — 23 passed; composition, History core/input/config/canonical
    owner, run-analysis, Worker protocol, orthogroup identity, record layout,
    source/selector/drawer, Session request/authority/file/regeneration,
    active-state, metadata, cache/resource, settings-only, and runtime-parity
    suites — passed;
  - canonical Session CLI compatibility — 4 passed;
  - targeted current lifecycle and released-v40 Chromium checks — passed;
  - focused Ruff, JavaScript syntax checks, and `git diff --check` — passed;
  - `node tests/web/architecture-contracts.test.mjs` — 137 passed; and
  - `node tools/check-web-change-budget.mjs` — Gate PASS, Review REQUIRED for
    the intended adapter export, two controller presentation refs, nine touched
    production files, and production net additions, with zero blocking
    violations, privileged expansions, import cycles, dependency/vendor/binary
    changes, or guard changes.
- Baseline/reference evidence: `TestOutputComparison` retained the same two
  previously documented circular label-binding metadata mismatches and passed
  its other 14 cases. `tests/reference_outputs/`,
  `examples/gbdraw_social_preview.png`, Session fixtures, `dist/`, and
  `gbdraw.egg-info/` were not changed; the prepared browser wheel remained a
  gitignored local test artifact.
- Architecture evidence: the existing controller, canonical request/Session
  codec, request planner/renderer, composition adapter, Result admission, and
  generic generated-artifact History transaction retain their single owners.
  Semantic mutation owners only notify the controller. No alternate resolver,
  Worker, history stack, persistence path, affine transform, framework, build
  step, or compatibility writer was added, and no superseded lifecycle path
  remains.
- Product Impact: `IMPLEMENT_EXISTING_AUTHORITY`; S06 realizes
  `PD-OI-026`–`PD-OI-031` lifecycle, Reset/history, stale recovery, and
  reader-only compatibility outcomes. The released fixture had sufficient
  stable feature and transform metadata for the authorized reader-only bridge,
  so no save-continuation choice or Product Decision Pack was required.
- Remaining scope: S06 has no known implementation blocker. S07 final
  end-to-end acceptance, public documentation, and release-readiness work
  remain pending and were not started.
- Next permitted session: S07 only.

### S07 evidence — 2026-09-24

- Authority/base: fetched `origin/issue-561-similarity-alignment` and `origin/dev` before work. The clean work branch and remote started at S06 commit `568b072902cb7b0d18aeeabc0dd193f8a97418da`; that commit was an ancestor of both. Accepted authority merge `a9eaeadd105e0c26e46086626feaa695bdd33c94` was an ancestor of both `origin/dev` and the work branch. S01–S06 were complete; S07 was the first pending row. Reviewed all eight branch commits and the full diff from the accepted authority base; newer unrelated `origin/dev` commits were not used as S07 authority.
- Integration closure: removed the superseded group-string Reset buttons/handler and their `runAnalysis()` route from `index.html`, `app-setup.js`, and `orthogroups.js`. The single current Reset remains in the active-plan inspector and uses the alignment controller's generated-artifact transaction. The right-drawer owner now keeps Similarity groups available while an active plan exists, even if the latest Result has no current group rows; this makes inspection and Reset reachable. Focused drawer and real-browser regressions cover that continuation.

| Issue #561 criterion | Authority | Production owner | Executable / manual evidence |
| --- | --- | --- | --- |
| 1. Exact clicked non-representative reference | `PD-OI-026` | `app/feature-editor/svg-actions.js`, `app/similarity-alignment.js`, `layout/similarity_alignment.py` | `test_exact_non_representative_reference_and_zero_or_one_candidate`; Chromium popup exact `b0` and drawer exact-reference selector |
| 2. Another valid member differs from a missing member | `PD-OI-026` | shared resolver, `web_support/similarity_alignment.py` | pure resolver zero/one/multiple tests; Web adapter crop facts; Chromium 1 aligned / 1 unchanged summary |
| 3. No target anchor leaves offsets/orientation unchanged | `PD-OI-026`, `PD-OI-027` | `api/record_planning.py`, `diagrams/linear/assemble.py` | `test_absolute_translation_formula_preserves_reference_skipped_and_every_y`; `test_final_placements_align_anchors_and_bounds_do_not_clip` |
| 4. Explicit Select overrides automatic choice | `PD-OI-026` | shared resolver and Web controller | pure resolver explicit-choice cases; `similarity-alignment-actions.test.mjs` Select/Skip; Chromium ambiguity Select |
| 5. Inparalog ambiguity requires Select/Skip without heuristic ranking | `PD-OI-026` | shared resolver, Web controller/dialog | direct/multiple RBH, edge-direction, non-RBH/multi-hop, deterministic-order resolver tests; Chromium desktop/narrow dialog, disabled Apply, Escape, focus, candidate preview |
| 6. Unselected copies keep membership, edges, and record-relative positions | `PD-OI-026`, `PD-OI-027` | resolver decision only; `api/record_planning.py`, `diagrams/linear/assemble.py` | pure resolver and rendering skipped/missing assertions; Chromium only selected target aligns; no analysis-group mutation owner was added |
| 7. Align changes target X only | `PD-OI-027` | `api/record_planning.py`, `diagrams/linear/assemble.py` | absolute-translation, same-row/bounds, final geometry-consumer, and strand-matrix tests; Chromium Align summary |
| 8. Align & orient uses displayed strands and whole-record reverse | `PD-OI-027` | `layout/similarity_alignment.py`, `api/record_planning.py`, existing record transform | strand matrix, effective orientation, post-reverse center, text/readability and composition tests; Chromium successful orient action with one usable member |
| 9. `rev` remains derived from source-relative effective orientation | `PD-OI-027`, `PD-OI-028` | `app/similarity-alignment.js`, existing record presentation | `similarity-alignment-actions.test.mjs` inspector/rev cases; S05/S06 real-browser evidence in this ledger |
| 10. Cancel/Reset/Undo/Redo/Session are atomic and offsets do not accumulate | `PD-OI-028`–`PD-OI-030`, `PD-OI-016`, `OIPC-C07` | Web controller, generic `history-snapshot.js` / `run-analysis.js`, request/Session codecs, Result admission | controller cancellation/stale/failed/Reset/reorder tests; rendering idempotence; schema-8/current Session and released-v40 tests; Chromium Cancel snapshot, Save/fresh Load without Worker, UI Reset, Undo/Redo |
| 11. Alignment changes display only; no new LOSATP/group inference | `PD-OI-018`, `PD-OI-031` | existing Worker operation, `web_support/similarity_alignment.py`, CLI adapter | Worker protocol and controller job-count assertions; CLI exact/noninteractive and typed Python tests; reviewed full production diff for deferred paths |

- Extended acceptance: pure resolver tests cover one versus multiple distinct direct RBHs, both edge directions, hidden-but-mappable versus cropped-out centers, duplicate identities, explicit Skip, and input-order independence. Renderer tests cover positive/negative X, same-row records, normalize/center combinations, canvas growth, ruler, definitions, tracks, labels, annotations, comparison edges, reference invariance, and repeated Generate. Node controller and S06 real-browser evidence cover Align A/B/Reset/Undo/Redo, source/crop/selector/manual orientation/drag notices, stable reorder, stale reference/target repair, Result retention, and no added LOSATP work. Current Chromium evidence independently covers desktop/narrow UI, no-dialog one-candidate and ambiguous paths, both actions, inspector reachability, cancellation, Session lazy preview, UI Reset, Undo/Redo, and page/console error and unhandled Promise rejection absence.
- Request/Session boundary: `session_request_codec.py` and Web `session-request.js` remain the current schema-8 typed writer; `session_compat.py` and Web `legacy-similarity-alignment.js` are bounded released-schema readers. Current request rejects group-only input. The legacy string persists only in private compatibility/mutable-intent staging; current request and Session writers emit no legacy alignment string. Released v39/v40/schema-5 fixtures materialize to typed state, malformed/unmappable data fails explicitly, and fresh Load-only preview starts no Worker. Fresh Generate and Session replay both enter the canonical render-request boundary and the existing Result admission path.
- Public surfaces/deferred scope: Web, CLI, and typed Python use `resolve_similarity_alignment` through the existing Worker/Python route. JavaScript presents resolver candidates and stores user choices; it has no separate ranking policy. Strict CLI ambiguous input returns actionable candidate IDs without prompting. No anchor TSV, Collinear alignment UI, smart alignment, synteny propagation, multi-hop inference, new Worker, alternate history stack, affine transform, second persisted plan, framework, or build step was introduced.
- S07 files: four Web production owners (`index.html`, `app-setup.js`, `orthogroups.js`, `right-drawer.js`); five test owners (documentation contract, output comparison, Gallery reproduction, drawer, Chromium journey); eight existing public documentation/navigation/release pages; eight regenerated recipe SVGs and one existing Gallery arrow SVG; this ledger. No Session fixture, tracked reference SVG, social preview, dependency, vendored runtime, or guard changed.
- Documentation: updated existing Web, CLI, typed Python, typed request, Session compatibility, navigation, and release-note owners. The literal `S07-PY-01` typed example ran from a clean temporary directory and is enforced by `test_typed_similarity_alignment_example_runs_from_clean_directory`; generated CLI help passed `python tools/update_cli_reference_help.py --check`. No new FAQ or minimal public smoke figure was needed.
- Public figure review: the first combined gate exposed ten stale public-artifact/arrow-variant checks; none was grandfathered. All eight stale recipe SVGs were regenerated by their documented clean-directory runners. Seven differed only by renderer-owned record metadata; H-CLI-07 also moved record definitions with final record placement. The existing Gallery arrow SVG was reproduced from its Session, retained labels, legend, metadata, quantitative/comparison context, and was rendered and inspected at readable scale; its current-code baseline differs in exactly 152 arrow feature paths, with all other attributes equal. Its circular peer differs in exactly 15 paths. Both refreshed figures remained readable with no clipping. The output comparison strips only the three new nonvisual record-key/translation attributes; focused rendering and composition tests assert those attributes. Tracked reference outputs and `examples/gbdraw_social_preview.png` were unchanged.
- Verification (required gates): `ruff check gbdraw/` — PASS; `python -m pytest tests/ -v -m "not slow"` — 6,184 passed, 17 skipped, 11 deselected, 23 warnings in 1000.00s; `python -m pytest tests/test_output_comparison.py::TestOutputComparison -v` — 16 passed, read-only; `node tests/web/architecture-contracts.test.mjs` — 137 passed; `node tools/check-web-change-budget.mjs` — Gate PASS / Review REQUIRED; `node tools/check-web-change-budget.mjs --base a9eaeadd105e0c26e46086626feaa695bdd33c94` — full branch Gate PASS / Review REQUIRED; `python -m build` — source distribution and wheel built; `npx playwright test tests/web/similarity-alignment-ui.playwright.spec.js --project=chromium --workers=1` — 3 passed.
- Focused verification: `python -m pytest tests/test_similarity_alignment.py tests/test_similarity_alignment_rendering.py tests/test_similarity_alignment_web_adapter.py -q` — 59 passed; `python -m pytest tests/test_session_request_codec.py tests/test_session_compat.py tests/test_api_session.py tests/test_session_io.py tests/test_documentation_contracts.py tests/test_documentation_reference_contracts.py -q` — 432 passed, 10 Biopython warnings; `node --test tests/web/right-drawer.test.mjs tests/web/similarity-alignment-actions.test.mjs tests/web/composition-layout.test.mjs tests/web/history-config-restore.test.mjs tests/web/run-analysis-simple-path.test.mjs tests/web/session-request.test.mjs tests/web/diagram-generation-worker.test.mjs tests/web/session-authority.test.mjs tests/web/session-draft-authority.test.mjs` — 38 passed; the additional focused Web history/record-layout/cache/runtime-parity set — 85 passed; the two Gallery reproduction checks — 2 passed. JavaScript syntax checks, `python tools/update_cli_reference_help.py --check`, and `git diff --check` passed. The final browser wheel was prepared and remains gitignored.
- Browser visual review: disposable `desktop-ambiguity.png`, `narrow-ambiguity.png`, and `narrow-final.png` from the Chromium spec were inspected at readable scale. The dialog and Apply control fit the 720x740 viewport, focus returns on Escape, and the active-plan inspector and Reset remain reachable with zero group rows. No screenshot was added to tracked public documentation.
- Architecture/Product: `IMPLEMENT_EXISTING_AUTHORITY` under `PD-OI-026`–`PD-OI-031`; the active-plan drawer fix realizes accepted inspection/Reset without choosing a new outcome. The shared Python resolver owns candidate decisions, one Web controller owns orchestration, canonical codecs own persistence/compatibility, request planning/renderer own transforms, one composition adapter owns record deltas, Result admission owns failure isolation, and generic History owns artifact transactions. Superseded group-string Reset was removed. S07 working-tree budget Review REQUIRED is solely the removal of one reactive computed declaration; no new owner, export, import cycle, permission, dependency, guard, or runtime path. The complete authority-base differential is Gate PASS / Review REQUIRED for the planned typed-plan/legacy-reader inventories, Session/compatibility paths, 26 production files, 2,871 gross changed lines, and 2,561 net additions; each maps to the single owners listed above. No Architecture or Product exception was triggered.
- Residual: no known accepted-behavior gap after required gates pass. Slow tests were outside the explicitly required `not slow` gate. PR, merge, release, tag, and deployment remain outside S07.

### PR #580 integration evidence — 2026-09-24

- Refreshed `origin/dev@f13cc6709e668b764a13be2ae2f8441f581246d5` and opened PR #580 from `issue-561-similarity-alignment` at S07 head `af9341f88edce035052bb1cbce80a0bb3693b3f2`. Integrated current `dev` into the work branch and resolved 32 conflicts without reverting either branch's separate feature work.
- Current writers use Session 44, request schema 8, and feature catalog schema 4. Previously written Session 44/schema-7/catalog-4 documents remain readable through the existing Session boundary; Gallery publication promotes their request to schema 8 and retains the saved Result. Development-only Session 43 remains rejected, as on `dev`. The CLI legacy alignment projection and catalog promotion both remain active. Python, Node, and real-browser regressions cover this overlap.
- Regenerated conflicting H-CLI-04/07/08/12, H-PY-05, T-PY-08, and Gallery arrow artifacts from their documented clean-directory recipes. Their Session samples report version 44/schema 8/catalog 4 where a catalog exists. The combined CLI and Gallery figures were rendered at 1,800 px and visually inspected; labels, comparison context, legend, and record metadata remained readable.
- Integrated verification: `python -m pytest tests/ -v -m "not slow"` — 6,225 passed, 17 skipped, 11 deselected, 23 warnings in 1305.44s; read-only `TestOutputComparison` — 16 passed; documentation and example reproduction — 35 passed; `node --test tests/web/*.test.mjs` — 665 passed; architecture contracts — 137 passed; Chromium Issue #561 and released Session 44 journeys — 4 passed; Python released Session 44 typed-reader regression — 1 passed; public reference contracts — 12 passed; `ruff check gbdraw/`, CLI help check, and `python -m build` passed. The working-tree Web budget against `origin/dev` with `architecture-change` profile is Gate PASS / Review REQUIRED, with no blockers; final committed-head and remote CI evidence follow in PR #580.
- The first integrated-head PR smoke run passed 12 of 13 browser journeys and exposed one stale test expectation for the current request schema in the AT-skew Session contract. Its expected schema was updated from 7 to 8; the exact previously failing `@pr-smoke` journey passed with the PR smoke configuration. No runtime behavior changed in that follow-up.
- The second integrated PR smoke run passed all 13 browser journeys, then Gallery first-Generate parity exposed a schema-7 generated comparison setting retained during current Generate. The current request writer now omits the legacy-only field. Session import preflight uses the existing schema promotion reader to materialize a non-null legacy selection before Generate; its exact stable feature ID matches the existing feature-catalog migration. A frozen released schema-7 Session import regression checks the typed reference. Gallery parity ignores only the three new nonvisual record-key/translation SVG attributes, matching read-only output comparison.
- The existing BGC Gallery session with legacy `og_1` alignment was replayed through `refresh_gallery_sessions.py` into Session 44/schema 8/catalog 4. Its saved Result, source SVG, interactive SVG, thumbnail, examples metadata, and artifact hashes were regenerated from that same session. The SVG has exactly four intentional definition-placement differences from the older published source after nonvisual metadata is excluded. Its 1,800 px rendering was visually inspected: all five record definitions, comparative ribbons, labels, legend, and title remain readable without clipping. The Gallery command was run from a clean directory and produced schema 8 with the same exact `CAG38695.1` reference; the saved Session reproduces the reviewed Result with its stored comparison and editor state.
- Focused final evidence: all 9 Gallery first-Generate Chromium cases passed; `pytest -m "gallery and not slow"` passed 103 tests; Gallery semantics and refresh tests passed 56; Web unit tests passed 666; architecture contracts passed 137; the frozen released Session import test passed; `ruff`, artifact-manifest verification, and the working-tree Web budget (`PASS / Review REQUIRED`, no blockers) passed. The final package build passed; exact-head remote CI remains the merge gate.
- A further browser sweep of five BGC Session workflows found that switching the saved Similarity Group comparison to pairwise left the old alignment plan active; the second Generate then requested stale-reference repair. Existing comparison-plan, program, and LOSATP-mode mutation routes now notify the single alignment controller on an actual change. The exact failing lazy-Worker and repeated-Generate journey passed after the fix, as did 36 focused comparison/alignment unit tests and the BGC first-Generate parity case. The `3d7169d8` PR head had passed all required GitHub checks before this additional fix; the new head requires its own checks before merge.
- Product/architecture: the integration realizes accepted Issue #561 authority alongside the `dev` feature-anchor/catalog changes, with no new Product outcome. The shared resolver, canonical request/Session codecs, existing Gallery publication owner, current Result admission, and Web alignment controller retain their stated responsibilities. No extra Worker, dependency, guard, or parallel alignment plan was added.

## 17. Handoff rules

- Treat each implementation session as one reviewable commit candidate.
- Do not push, create or merge a PR, publish, tag, or deploy unless the active
  task explicitly authorizes that external action.
- A later session may revise engineering details when current code proves the
  proposal inaccurate, but it may not change the accepted Product outcomes.
- If an accepted outcome cannot be implemented without changing another
  Product effect, stop that convergence and prepare a new Decision Pack.
- After the final implementation, provide an English proposed commit title and
  concise summary as required by repository guidance.
