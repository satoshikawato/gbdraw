# Issue #561 Similarity Group alignment — master implementation plan

Status: implementation plan; Product authority merged; runtime work not started

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
| S01 | pending | | | |
| S02 | pending | | | |
| S03 | pending | | | |
| S04 | pending | | | |
| S05 | pending | | | |
| S06 | pending | | | |
| S07 | pending | | | |

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
