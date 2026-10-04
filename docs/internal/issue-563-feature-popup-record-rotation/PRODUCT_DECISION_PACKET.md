# Product Decision Pack — Feature-popup record rotation

Status: unresolved; non-authoritative working packet

## Identity

- Concern key: `diagram-generation.feature-popup-record-rotation`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `a9eaeadd105e0c26e46086626feaa695bdd33c94`
- Prepared for branch/head:
  `issue-563-feature-popup-record-rotation-20260922` @ `3672a518e09339bffc2cc70cdee42add9103ccb2`
- Related issue: `#563`

## Trigger

Issue #563 and its saved implementation plan propose a new feature-popup
workflow that changes discoverability, editable record-wide transforms,
Session provenance, atomic History, regeneration isolation, and recovery.
Existing authority preserves shared request/Result/History constraints but does
not select whether this new public workflow is supported.

## Authority search

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact concern/map | Canonical request and current Result continuity are authority-covered. | No feature-popup rotation concern/outcome. |
| Active durable `BD-###` | None. | No selection. |
| Option Integrity Product Contract revision 10 | Cross clauses preserve explicit clearing and failure isolation; Similarity alignment decisions are narrowly scoped. | No Issue #563 record. |
| Domain/scientific/integrity rule | Source-coordinate and strand semantics constrain any implementation. | Does not select the affordance or lifecycle. |
| Released compatibility | Sidebar rotation and Session 42/catalog 3 are supported. | Does not authorize the popup workflow. |
| Eligible exact-head current decision | None. | A material new affordance cannot use the narrow PR-local route. |
| Current code/tests | Existing transform, request, Result, History, and LOSAT cache make the proposal feasible. | Evidence only. |
| Issue #563 | Complete proposal and ten acceptance scenarios by `satoshikawato`. | Issue text alone is not durable base authority. |

Result: `UNRESOLVED`

Procedural classification: `PRODUCT_DECISION_REQUIRED`

## User journey

- Actor: Web user
- Context: A fresh or loaded interactive Result contains a source-bound feature
  on a complete effectively circular record.
- Goal: Rotate only that record from the open feature popup and optionally
  orient the clicked feature forward without committing unrelated drafts.
- Entry point: feature search or direct feature popup
- Major steps: open popup; open Record actions; select anchor/offset/orientation
  or feature-end; review source-coordinate preview; Apply or Cancel; inspect the
  fresh Result; optionally Undo/Redo or save/load the Session.
- Checkpoints: popup target binding, preview, Apply submission, Result
  admission, failure/cancel/stale recovery, Undo/Redo, Session round trip.
- Failure/recovery: the old transform and Result remain current; no History
  entry is added.
- Persistence/export: absolute start and orientation plus non-authoritative
  provenance persist; source sequences/annotations and source export do not
  change.
- Next action: continue search/edit/export/save after success; correct an
  explicit reason or regenerate metadata after a disabled/stale state.

## Current observed behavior

Users can edit per-record circular topology and display start in the sidebar
and can use global-selection 5′/midpoint shortcuts. There is no popup-local
3′/offset/orient-forward/feature-end transaction. Normal Generate already
reuses compatible LOSAT evidence after display start and reverse-complement
changes.

## Non-waivable constraints

- Architecture: one canonical request owner, candidate admission path, SVG
  admission owner, History owner, and coordinate-transform owner.
- Security/privacy: local same-origin assets; genome data remains in-browser;
  SVG uses the shared sanitizer.
- Scientific correctness: original source coordinates, exact part traversal,
  non-negative circular modulo, no guessed ambiguous ordering.
- Persisted compatibility: request schema remains 7; released catalog 3 data is
  read conservatively; branch-only v43 is not granted a reader without history
  evidence.
- Performance/resource safety: no forced LOSAT skip; matching raw evidence is
  reused and executor additions remain zero.
- Evidence: all ten Issue scenarios, AC-01–AC-20, failure paths, actual popup
  DOM journey, and target-only request isolation.

## Choice A — POPUP-RECORD-ROTATION

- Complete normative outcome: Adopt the complete Issue #563 outcome and
  AC-01–AC-20 exactly as recorded in the saved implementation plan. The open
  popup's explicit `(recordKey, biologicalFeatureId)` targets one complete
  effectively circular record in either diagram mode. It offers 5′, covered
  midpoint, 3′, signed strand-relative offset, optional absolute
  orient-forward, and distinct feature-end placement with preview, Apply, and
  Cancel. It uses the existing absolute display transform and target-only
  committed request candidate. Apply is one atomic transform/Result History
  transaction; failure/cancel/stale/superseded work is a no-op.
- Preserved effects: source sequence/annotation/qualifiers/biological identity;
  target-external records and layout; pending form edits; shared request,
  Worker, sanitizer, Result, History, and transform paths; LOSAT evidence;
  search and post-generation continuation.
- Added effects: popup-local record action, 3′ and feature-end choices, signed
  offset, optional orientation, source-coordinate preview, precise disabled
  reasons, Session provenance, atomic Undo/Redo.
- Lost effects: none required.
- Retired effects: none required. The sidebar workflow remains available and
  shares the resolver.
- Discoverability/accessibility: rich/simple popup, explicit Record actions,
  keyboard operation, visible reason text, and 390 px support.
- Canonical state update: success writes absolute per-record display start and
  reverse-complement state; provenance is non-authoritative.
- Undo/Redo consequence: one entry restores/reapplies transform, provenance,
  and Result.
- Session/regeneration consequence: Session 44/catalog 4 persists the outcome;
  released catalog 3 is conservative; later Generate uses absolute values.
- Export/artifact consequence: rendered artifacts change geometry only; no
  re-originated source-file export is added.
- Validation/error consequence: invalid offsets and unsafe operations are
  explicit; no truncation, guessing, or global-selection fallback.
- Failure/recovery consequence: previous committed Result and transform remain.
- Scientific-output consequence: display geometry changes through the existing
  transform; biological coordinates and identity do not.
- Cache/provenance consequence: compatible LOSAT evidence is reused; provenance
  never becomes rendering authority.
- Performance consequence: normal comparison planning/cache lookup remains;
  zero additional executor jobs for transform-only changes.
- Compatibility consequence: one bounded released catalog-3 Session reader;
  request schema and Worker protocol remain unchanged.
- Architecture consequence: no owner/path excess; one compatibility-burden
  increase requiring exact-head architecture approval.
- Evidence available/missing: feasibility and zero-job baseline pass; complete
  runtime/browser acceptance is missing until implementation is authorized.
- Residual risk: Must be supplied by the Product Decision Owner; do not infer.
- Route: `DURABLE_AUTHORITY_REQUIRED`
- Next action if selected: supply the complete receipt below, merge an
  authority-only record, explicitly update the runtime branch base, then resume
  Session 02.

## Choice B — RETAIN-SIDEBAR-ONLY

- Complete normative outcome: Keep the existing sidebar display-start controls
  and global-selection 5′/midpoint shortcuts; do not add popup-local record
  rotation, signed offset, 3′, feature-end, or orient-forward actions.
- Preserved effects: all current record rotation, Generate, Result, Session,
  History, comparison, and source behavior.
- Added effects: none.
- Lost effects: none relative to current behavior; the Issue #563 proposed
  workflow remains unavailable.
- Retired effects: the proposed Issue #563 delivery, not an existing supported
  effect.
- Discoverability/accessibility: unchanged.
- Canonical state update, Undo/Redo, Session, export, validation, recovery,
  scientific output, cache, compatibility, and architecture: unchanged.
- Evidence available/missing: current characterization passes.
- Residual risk: Users continue the multi-step popup-to-sidebar workflow and do
  not receive the requested endpoint/orientation semantics.
- Route: no runtime change; close or defer the implementation request.
- Next action if selected: stop after Session 01 and record the Issue
  disposition outside this runtime branch.

## Comparison matrix

| Dimension | Choice A | Choice B |
| --- | --- | --- |
| Entry/discoverability | Direct popup action | Existing sidebar only |
| Immediate feedback | Resolved coordinate/orientation preview | Existing sidebar feedback |
| Canonical update | Target-only atomic candidate | Normal Generate from current form |
| Undo/Redo | Transform + Result in one entry | Existing separate draft/generation behavior |
| Session round trip | Absolute transform + provenance | Existing start/topology draft |
| Validation/error | Operation-specific reasons | Existing shortcut/control reasons |
| Failure recovery | Old transform/Result retained | Existing Result preservation |
| Accessibility | New keyboard/mobile popup surface | No new surface |
| Performance | Existing cache; zero transform-only jobs | Existing cache |
| Compatibility | Catalog 3 reader to catalog 4 | No new path |
| Architecture | Same owners/paths; CB +1 | No delta |
| Route | Durable authority required | No runtime change |

## Engineering recommendation

- Recommended option: Choice A, because it matches the maintainer-authored
  Issue and reuses existing absolute transform/request/admission owners.
- Engineering reason: it can add the requested workflow without a second
  renderer, Worker, request schema, SVG admission path, or History engine.
- This recommendation is not Product authority.

## Product Decision Owner response

The current implementation request indicates an intention to implement Choice
A, but the following complete receipt is still required. Every field must be
provided by an eligible owner; Codex will not infer omitted wording.

```text
PRODUCT_DECISION
Concern: diagram-generation.feature-popup-record-rotation
Scenario revision: 1
Choice: A / POPUP-RECORD-ROTATION
Rationale: <product-level reason>
Must preserve: <effects and affordances>
May retire: <none or explicit scope>
Accepted residual risk: <bounded risk or none>
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

A complete receipt must be serialized in an authority-only change. It cannot
authorize dependent runtime in the same candidate.
