# Option Integrity Product Contract

Status: active Product authority

## Authority metadata

- Contract ID: `OIPC`
- Contract revision: `2`
- Product Decision Owner: `satoshikawato`
- Decision date: `2026-08-28`
- Decision source: explicit Product Decision Owner selection of one (`1`) after
  review of the merged deterministic evidence
- Reviewed candidate:
  `03_INITIAL_OPTION_INTEGRITY_PRODUCT_CONTRACT_CANDIDATE.md`
- Candidate SHA-256:
  `26b41219ca04ff26b56a29e11aa4be74c6030b0a22e386332afc28ac7a80623f`
- Reviewed evidence:
  [`COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md`](./COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md),
  merged by PR `#422` at `878c62ba17c61c45cd0adbd05cbd9fb36306db9d`
- Approved decision IDs: `PD-OI-001`, `PD-OI-002`, `PD-OI-003`,
  `PD-OI-004`, `PD-OI-005`, `PD-OI-006`, `PD-OI-007`, `PD-OI-008`,
  `PD-OI-009`, `PD-OI-010`, `PD-OI-011`, `PD-OI-012`, `PD-OI-013`,
  `PD-OI-014`, `PD-OI-015`, `PD-OI-016`, and `PD-OI-017`
- Initial candidate modification: `PD-OI-014`, as recorded below
- Current revision change: `PD-OI-007`, as recorded below
- Records remaining `EVIDENCE_REQUIRED`: none
- Excluded records: none

This contract owns the user-observable outcomes recorded below. It does not
select source files, classes, module paths, canonical call edges, delivery
order, cache implementation, or test design. Code, tests, fixtures,
screenshots, and historical behavior are evidence, not Product authority.

The mapped concern `product.canonical-render-request-boundary`, scenario
revision `1`, remains governed by its existing selected option
`canonical-typed-request-boundary`. This contract references that authority
where applicable and does not replace or duplicate it. No active `BD-###`
decision governs the concerns recorded here at this revision.

## Interpretation and lifecycle

Active records use only these statuses:

- `ACCEPTED`: the complete outcome is normative.
- `EVIDENCE_REQUIRED`: evidence-only work may proceed, but implementation must
  not select the pending outcome.
- `DEFERRED`: the named capability is outside the current delivery scope.
- `UNSUPPORTED`: the named input or journey must be rejected explicitly.

Authority precedes dependent runtime implementation. A runtime change cites
authority already present on its base and does not change the decision it
implements. Evidence precedes a decision when a record is
`EVIDENCE_REQUIRED`; evidence does not select its own outcome.

To correct an active outcome, use an authority-only replacement. Increment the
scenario revision, identify the prior decision and revision in `Supersedes`,
record the complete replacement outcome, and merge that authority before
changing dependent runtime. Git history retains the former text; the active
contract does not accumulate superseded records.

## Cross-surface clauses

### OIPC-C01: Omitted and explicit values

- Omission uses the public typed default.
- An explicitly supplied default is execution-equivalent to omission.
- An explicit valid non-default changes the documented execution or
  presentation behavior.
- Invalid explicit values are rejected, not silently coerced.

### OIPC-C02: Requested and effective values

When automatic or context-dependent resolution exists, requested intent and
effective execution are both retained. Requested intent is not overwritten by
the effective subtype.

### OIPC-C03: Consume or reject

Every accepted public value reaches its real consumer or is rejected before
execution. A surface must not accept and silently ignore a value.

### OIPC-C04: Request, execution, cache, and artifact agreement

For each applicable field, canonical request intent, resolved execution
values, actual helper invocation, correct stage-specific cache identity,
Session data, and artifact metadata agree. Requested/effective differences are
explicit. The mapped `product.canonical-render-request-boundary` concern
continues to own canonical Web request continuity.

### OIPC-C05: Preservation of valid intent

A valid value is not deleted because a surface lacks an editor. The generic
surface disposition is `EDITABLE`, `READ_ONLY`, `PASS_THROUGH`, or
`UNSUPPORTED`. Imported comparison reconstruction uses the more specific
states `EDITABLE`, `PRESERVED_READ_ONLY`, and `DECISION_REQUIRED` defined in
`PD-OI-008`.

### OIPC-C06: Explicit replacement and clearing

Replacement and removal are explicit actions. Empty controls, missing
properties, inactive modes, failed reconstruction, and failed generation do
not imply deletion.

### OIPC-C07: Failure isolation

Failed, canceled, superseded, or stale generation does not replace the last
successful Result or committed request.

### OIPC-C08: Evidence is not authority

A test or historical implementation that contradicts an active decision is
corrected. Passing evidence does not make incorrect behavior normative.

## Product Decision records

### PD-OI-001: LOSATP Candidate limit fresh default

- Concern key: `diagram-generation.losatp-candidate-limit-default`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: The fresh Candidate limit is `None`. An unbounded request
  executes without an undocumented finite cap.
- Rationale: Raw candidate evidence must not vary because one surface
  substituted a hidden limit.
- Must preserve: Explicit finite values; truthful requested/effective
  metadata; deterministic cancellation and errors.
- May retire: Hidden browser-only finite caps and documentation that calls
  them the Product default.
- Accepted residual risk: Unbounded browser work can be expensive. A
  deterministic warning may request explicit continuation but must not rewrite
  the requested value.
- Acceptance contracts: `OIC-001`, `OIC-005`, `OIC-013`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-002: Candidate-limit scope and GUI placement

- Concern key: `diagram-generation.losatp-candidate-limit-scope`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Candidate limit is one common LOSATP Advanced setting.
  Pairwise, Similarity-group, and Collinear presentation modes do not own,
  duplicate, reset, or rewrite it.
- Rationale: Candidate limit controls raw search evidence rather than one
  presentation.
- Must preserve: One value across presentation changes and Session round
  trips.
- May retire: Presentation-specific Candidate-limit state and duplicate
  controls.
- Accepted residual risk: Advanced placement is less prominent than a
  mode-specific field.
- Acceptance contracts: `OIC-001`, `OIC-002`, `OIC-006`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-003: Pairwise display max hits

- Concern key: `diagram-generation.pairwise-display-max-hits`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Pairwise display max hits has a fresh default of `5` and
  is applied after threshold filtering to Pairwise result selection only.
- Rationale: Presentation density must not change raw search evidence.
- Must preserve: Raw-search cache reuse when only this value changes; explicit
  Session values.
- May retire: Aliasing Pairwise display max hits to Candidate limit or member
  hits.
- Accepted residual risk: Additional qualified hits remain absent from the
  Pairwise view while retained in raw evidence.
- Acceptance contracts: `OIC-002`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-004: Similarity/Collinear member hits per protein

- Concern key: `diagram-generation.member-hits-per-protein`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Member hits per protein has a fresh default of `5` and
  controls directional hits used to construct Similarity-group and Collinear
  evidence.
- Rationale: Group/block construction is derived analysis, distinct from raw
  search and Pairwise display.
- Must preserve: Raw-search cache reuse when only this value changes;
  derived-cache invalidation; explicit Session values.
- May retire: Aliasing this field to Candidate limit or Pairwise display max
  hits.
- Accepted residual risk: Changing it can change grouping/block output and
  must be visible in provenance.
- Acceptance contracts: `OIC-002`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-005: Supported Collinear value domains

- Concern key: `diagram-generation.collinear-value-domains`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Web execution consumes every currently public Collinear
  enum value, including the current equivalents of unit mode `auto`, `cds`,
  `locus`; anchor mode `all`, `one_to_one`, `rbh`; and merge orientation
  `strand`, `order`, `either`. Exact names are verified against the current
  typed API before implementation.
- Rationale: A valid public value must reach execution or be rejected
  explicitly.
- Must preserve: Typed Python validation; requested `auto`; separately
  reported effective resolution.
- May retire: Browser or Worker branches that coerce one valid enum to another.
- Accepted residual risk: Combinations are covered primarily by unit/contract
  tests rather than browser cases.
- Acceptance contracts: `OIC-003`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-006: Collinear search-scope fresh default

- Concern key: `diagram-generation.collinear-search-scope-default`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: The fresh search-scope default is `adjacent`; `all`
  remains an explicit supported value.
- Rationale: Adjacent comparison is the primary Web journey while explicit
  all-pairs analysis remains available.
- Must preserve: Explicit `all` across CLI, Python API, Web, and imported
  Session.
- May retire: Conflicting fresh defaults across surfaces.
- Accepted residual risk: Fresh output can differ from a historical surface
  that used implicit `all`; released Sessions must be preserved or migrated
  explicitly, not reinterpreted.
- Acceptance contracts: `OIC-003`, `OIC-004`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-007: Collinear merge-conflict fresh default

- Concern key: `diagram-generation.collinear-max-conflicts-default`
- Scenario revision: `2`
- Status: `ACCEPTED`
- Normative outcome: The fresh `max_conflicts` default is one (`1`). Explicit
  zero (`0`) remains supported. Omission and explicit one are
  execution-equivalent.
- Rationale: Merged deterministic evidence established the consumer behavior
  of both values, and the Product Decision Owner selected one after reviewing
  that evidence.
- Must preserve: Explicit zero and one; retained singleton anchors; agreement
  among request intent, execution, round trip, and provenance; reproducible
  evidence for the merge-threshold effect.
- May retire: Conflicting fresh omission defaults across surfaces.
- Accepted residual risk: At one, compatible clusters may merge across one
  retained interior singleton where zero keeps them separate. The singleton
  remains in the result, and the selected behavior remains visible in
  provenance.
- Acceptance contracts: `OIC-004`, `OIC-005`.
- Evidence: [`COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md`](./COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md),
  merged by PR `#422` at `878c62ba17c61c45cd0adbd05cbd9fb36306db9d`.
- Supersedes: `PD-OI-007`, scenario revision `1`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-008: Imported comparison reconstruction and resolution

- Concern key: `diagram-generation.imported-comparison-state`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: An imported committed comparison is `EDITABLE`,
  `PRESERVED_READ_ONLY`, or `DECISION_REQUIRED`. Exact reconstruction becomes
  editable without marking a user edit. A valid executable but non-projectable
  comparison remains usable through explicit inheritance. An ambiguous,
  incomplete, or non-executable comparison blocks Generate until explicit
  replacement or clearing. No comparison is silently converted to no
  comparison.
- Rationale: GUI limitations must not destroy valid work or make regeneration
  untrustworthy.
- Must preserve: Last successful Result; committed request; Save/export; exact
  executable data; explicit continuation; failure isolation.
- May retire: Silent empty-comparison fallback, mode-change resets, and
  reconstruction that invents defaults.
- Accepted residual risk: Some imported Sessions require read-only disclosure
  or explicit resolution before regeneration.
- Acceptance contracts: `OIC-007`, `OIC-013`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-009: GUI-unmanaged configuration overrides

- Concern key: `diagram-generation.gui-unmanaged-config-overrides`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Schema-known, valid, safe, active-mode-compatible
  configuration leaves without a GUI editor are preserved losslessly and
  disclosed. Unknown paths, non-leaf paths, unsafe keys, invalid literals, and
  active-mode-incompatible values are rejected explicitly.
- Rationale: Lack of an editor is not authorization to discard a valid typed
  value; generic unknown pass-through is unsafe.
- Must preserve: Managed siblings; valid imported leaves; typed validation;
  explicit reset.
- May retire: Silent dropping of valid GUI-unmanaged values and permissive
  unknown-path pass-through.
- Accepted residual risk: Users can see read-only settings they cannot edit in
  the current GUI.
- Acceptance contracts: `OIC-008`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-010: Circular record selection, crop, and reverse display

- Concern key: `diagram-generation.circular-record-transforms`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Circular Web workflows support selection of one record,
  valid region crop, and reverse-complement display. Arbitrary multi-record
  subset editing is outside this program.
- Rationale: These operations form one coherent single-record preparation
  journey.
- Must preserve: Source record order; stable record identity; deterministic
  selection; coordinate validation; Session round trip; no double reverse.
- May retire: Hidden accepted fields that never reach rendering and duplicate
  display-direction state.
- Accepted residual risk: Arbitrary multi-record subset editing remains
  unavailable.
- Acceptance contracts: `OIC-009`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-011: Linear scale and ruler-label font interaction

- Concern key: `diagram-generation.linear-scale-ruler-fonts`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Scale font size and ruler-label font size are linked by
  default and can be unlinked in Advanced settings. Imported Sessions with
  unequal values open unlinked and retain both. Explicit relink states the
  effect and copies the current scale font size to ruler-label font size once.
- Rationale: The public values are independent while linked fresh state keeps
  the common UI simple.
- Must preserve: Both explicit values; Session round trip; deterministic
  relink behavior.
- May retire: One-field aliasing that overwrites the other value.
- Accepted residual risk: One small linked/unlinked UI state transition is
  added.
- Acceptance contracts: `OIC-011`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-012: Circular record label and subtitle

- Concern key: `diagram-generation.circular-label-subtitle`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Explicit record label and subtitle values affect Circular
  rendering and are editable in the applicable Web journey. Explicit values
  override corresponding inferred title lines; empty values preserve inferred
  output.
- Rationale: Accepted public presentation fields must work or be rejected.
- Must preserve: Existing output when both values are empty; explicit values
  through Session round trip.
- May retire: Parsing and serialization paths that accept fields without
  consuming them.
- Accepted residual risk: Explicit text can change layout and requires targeted
  geometry review.
- Acceptance contracts: `OIC-010`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-013: Web compatibility meaning of `grid_column`

- Concern key: `diagram-generation.grid-column-compatibility`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: The Web normalizes `grid_column` to row-internal render
  ordering. Render-equivalent placement is guaranteed; exact numeric identity
  is not a Web round-trip promise.
- Rationale: Preserving a number with no distinct supported visual effect
  would add compatibility complexity without Product value.
- Must preserve: Relative order and rendered placement.
- May retire: Exact-number assertions without a supported visual difference.
- Accepted residual risk: A Web-exported Session can use different numeric
  columns while rendering equivalently.
- Acceptance contracts: `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-014: Direct non-adjacent Similarity and Collinear display links

- Concern key: `diagram-generation.non-adjacent-grouping-topology`
- Scenario revision: `1`
- Status: `DEFERRED`
- Normative outcome: Similarity groups continue to use all-vs-all
  protein-search evidence across every loaded record. Collinear continues to
  support both **Adjacent pairs** and **All records** evidence scope. An
  uploaded BLAST TSV remains assignable to every supported comparison edge.
  The only deferred capability is drawing a direct link or ribbon between
  arbitrary user-selected non-adjacent display rows for a Similarity-group or
  Collinear result. This Web display deferral does not make the typed API as a
  whole unsupported.
- Rationale: Direct links between arbitrary non-adjacent result rows require
  separate interaction, persistence, and recovery design. Existing evidence
  scopes and supported comparison-edge inputs are independent capabilities.
- Must preserve: Similarity-group all-vs-all evidence; both Collinear evidence
  scopes; uploaded BLAST TSV assignment to supported comparison edges;
  existing adjacent/all display workflows; accurate typed API documentation.
- May retire: Claims that the typed API lacks a capability solely because the
  Web does not expose it.
- Accepted residual risk: The Web cannot directly draw a Similarity-group or
  Collinear link/ribbon between arbitrary user-selected non-adjacent display
  rows until that interaction is designed.
- Acceptance contracts: `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-015: Advanced GUI exposure policy

- Concern key: `diagram-generation.gui-exposure-policy`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: API field existence alone does not require a GUI control.
  A control is added for a defined journey or credible usage evidence. Valid
  GUI-unmanaged values use `READ_ONLY` or `PASS_THROUGH` behavior.
- Rationale: GUI completeness is measured by supported journeys, not field
  count.
- Must preserve: Explicit disposition; lossless valid values; accurate
  unsupported/deferred scope.
- May retire: Field-count parity as an acceptance criterion and a
  repository-wide registry used only to enforce it.
- Accepted residual risk: Some advanced values remain non-editable in the Web.
- Acceptance contracts: `OIC-006`, `OIC-008`, `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-016: Generation failure, cancellation, and stale-result isolation

- Concern key: `diagram-generation.failure-isolation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: A failed, canceled, superseded, or stale Generate attempt
  does not replace the last successful Result or committed canonical request.
  The draft remains available, and the UI identifies that the displayed Result
  is the last successful one.
- Rationale: Failed work must not destroy valid output or present candidate
  state as committed.
- Must preserve: User inputs; last successful artifact; actionable error;
  ability to correct and retry.
- May retire: Pre-validation committed mutation, clearing output on failure,
  and partial-state admission.
- Accepted residual risk: Displayed Result can differ from the current draft;
  the UI discloses this.
- Acceptance contracts: `OIC-013`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-017: Active comparison appearance controls remain reachable

- Concern key: `diagram-generation.comparison-appearance-reachability`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: When a Linear comparison presentation is active,
  applicable Comparison Appearance controls remain reachable. Supported match
  styles include `ribbon` and `curve`, and match-height adjustment remains
  available where visibly effective. Switching analysis modes does not reset
  independent appearance values. Appearance settings affect rendering, not
  raw protein search or grouping semantics.
- Rationale: A supported rendering choice must remain discoverable and
  operable.
- Must preserve: Control reachability; accessibility; Session round trip;
  independence from analysis defaults and cache stages.
- May retire: Hidden active controls and unauthorized mode-change resets.
- Accepted residual risk: Advanced placement can require one additional
  disclosure action if the control remains discoverable.
- Acceptance contracts: `OIC-006`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

## Acceptance contract catalog

| Contract | Required meaning |
| --- | --- |
| `OIC-001` | Candidate limit is truthful; `None` remains unbounded and no hidden cap is applied. |
| `OIC-002` | Candidate, Pairwise display, and member-hit limits are independent and invalidate only the correct stages. |
| `OIC-003` | Every supported Collinear enum reaches the real typed Python analysis path. |
| `OIC-004` | Fresh defaults are consistent across surfaces and explicit imported values win. |
| `OIC-005` | Canonical request, resolved values, actual helper invocation, stage cache identities, and artifact provenance agree. |
| `OIC-006` | Required controls are discoverable, operable, mode-safe, persistent, and accessible. |
| `OIC-007` | Imported comparison intent is never silently cleared; unresolved state has explicit actions. |
| `OIC-008` | Valid GUI-unmanaged config survives and is disclosed; invalid or unknown input is rejected. |
| `OIC-009` | Circular selection, crop, and reverse use one effective path and round-trip without double reversal. |
| `OIC-010` | Circular label/subtitle are effective; empty values preserve inferred output. |
| `OIC-011` | Scale and ruler-label font sizes remain independently public with explicit linked-default behavior. |
| `OIC-012` | `grid_column`, deferred direct-link topology, and GUI surface scope are represented accurately. |
| `OIC-013` | Failed, canceled, or stale Generate preserves the committed request and last successful Result. |
| `OIC-014` | Product authority, Product Impact mapping, Architecture Ratchet, runtime owners, and evidence remain separate. |

## Residual-risk boundary

Accepted residual risk never authorizes a security vulnerability, silent
scientific-output corruption, loss of a must-preserve effect, deterministic
Architecture Ratchet failure, undocumented unbounded performance regression,
cache reuse across different execution semantics, artifact provenance that
disagrees with actual execution, or failure of a required acceptance contract.
