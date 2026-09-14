# Option Integrity Product Contract

Status: active Product authority

## Authority metadata

- Contract ID: `OIPC`
- Contract revision: `3`
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
- Revision 2 change: `PD-OI-007`, as recorded below
- Additional approved decision IDs: `PD-OI-018`, `PD-OI-019`
- Revision 3 addition: `PD-OI-018`, accepted by `satoshikawato` on
  `2026-09-13` after confirming the complete record/search outcome, no feature
  retirement, and the runtime/memory cost of complete comparisons. The initial
  approval and its date above continue to describe `PD-OI-001`–`PD-OI-017`.
  The maintainer subsequently specified the fresh defaults in `PD-OI-019`.
  On the same date, the maintainer replaced the default-five member cap in
  `PD-OI-004` with an unbounded default and clarified `PD-OI-018` to require
  source-file execution, reusable raw evidence after display transforms,
  compact record disclosures, and execution controls in LOSAT Settings.
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
- Scenario revision: `2`
- Supersedes: `PD-OI-004`, scenario revision `1`.
- Status: `ACCEPTED`
- Normative outcome: Member hits per protein defaults to `None` (unbounded)
  in fresh and reset state. A blank Web control means unbounded. Similarity
  groups and Collinear retain every threshold-qualified directional member
  hit without substituting a finite cap. An explicit positive integer remains
  a supported limit. Request, helper execution, provenance, and Session replay
  must distinguish an unbounded value from an explicit finite value.
- Rationale: The maintainer requested removal of the default five-hit cap
  because E-value and the other result thresholds already filter evidence.
  Member selection remains independent of raw search and Pairwise display.
- Must preserve: Raw-search cache reuse when only this value changes;
  derived-cache invalidation; explicit Session values.
- May retire: The default five-member cap and blank-to-five coercion; aliasing
  this field to Candidate limit or Pairwise display max hits.
- Accepted residual risk: Changing it can change grouping/block output and
  must be visible in provenance. Unbounded evidence retains the computation
  and memory cost accepted under `PD-OI-018`.
- Acceptance contracts: `OIC-002`, `OIC-004`, `OIC-005`, `OIC-006`.
- Decision source: Explicit maintainer follow-up requesting no default member
  count cap, while retaining threshold filtering.
- Owner and decision date: `satoshikawato`, `2026-09-13`.

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

### PD-OI-018: Complete Linear records, placement, and comparison scope

- Concern key: `diagram-generation.linear-record-universe-and-search-scope`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome:
  1. Without an explicit record selector or crop, Linear includes every record
     from each GenBank or paired GFF3/FASTA source. Enabling comparisons does
     not shrink that set or select only the first record.
  2. Record Layout exposes one independently editable placement per biological
     record, with an identifiable record label. File count is not record count.
     Exposing the controls preserves the existing row placement. Each uploaded
     source has one file-input card regardless of its record count; GFF3 and its
     paired FASTA remain one source. Per-record controls belong under that
     source or in Record Layout and must not appear as repeated file uploads.
     Each multi-record source's record list starts collapsed behind a compact, single-line
     `Number of records: N` disclosure. Expanding it exposes every record's
     controls; collapsing it changes no record selection or placement.
     Removing or replacing a source updates all of its records without leaving
     hidden records from the former source. Save and fresh Load preserve this
     distinction between source files and biological records.
  3. Adjacent comparison selects the Cartesian product between neighboring
     occupied display rows. Two records above three records means six pairs.
     This applies to nucleotide, translated-nucleotide, and Pairwise protein
     comparisons. Explicit selected subsets and uploaded TSV associations
     retain their meanings.
  4. Similarity groups and Collinear with scope `all` search every loaded
     record against every loaded record, including same-row and non-adjacent
     pairs. Required reverse and within-record evidence is retained. Display
     placement, including a single occupied row, does not restrict this search.
  5. Explicit comparison endpoints stay explicit through decoding and
     rendering, including endpoints whose numeric indices are consecutive.
  6. Save, fresh Load, regeneration, reordering, and cache reuse preserve the
     selected record set, endpoints, and independent placements. Explicit
     selection/cropping, comparison omission, and imported read-only intent
     remain supported. `OIPC-C07` governs failed, canceled, and stale work.
  7. Shared source bytes remain shared. Complete comparisons may require more
     jobs, but no hidden record or pair cap is permitted.
     The LOSAT execution unit is an input source file, not an individual
     record. A source job searches multi-sequence FASTA inputs containing the
     selected records, then routes hits to their record endpoints. Two sources
     in all-vs-all require four directed source jobs, including self and reverse
     searches; eight total records must not become 64 LOSAT invocations.
     TLOSATX records with different explicit genetic codes use compatible
     subsets within each source because each invocation accepts one query and
     one subject genetic code; records sharing those settings remain batched.
     Adjacent and explicit selections restrict retained record-pair evidence,
     not the source-level execution unit. Search arguments, source contents,
     and the actual searched database scope are part of raw-cache identity.
     Progress reports actual source jobs separately from biological record
     pairs. File batching must preserve cancellation and exact hit routing.
     Drawing-start changes and reverse-complement display reuse raw LOSATP
     results whenever source contents, selected biological regions, and search
     settings are unchanged. Only display coordinates and derived presentation
     are updated; display transforms do not become raw-search identity.
- Rationale: Prevent recurrence of incomplete-record search and per-record
  placement regressions reported by the maintainer.
- Must preserve: All seven outcomes above; `PD-OI-006` keeps the fresh Collinear
  scope `adjacent`, and `PD-OI-014` continues to govern evidence versus displayed
  links. A fresh no-comparison document does not gain comparison intent.
  Comparison places Run LOSAT across the top, with No comparison on the left
  and Upload BLAST TSV on the right of the row below. DOM and keyboard Tab
  order follow that visual order on desktop and narrow screens.
  Existing LOSAT Execution, Total threads, Parallel runs, and Threads per run
  controls remain editable in Comparison Settings when LOSAT is active.
  Selecting Run LOSAT opens Settings immediately, exposing the mode and its
  settings without another disclosure click. Restoring active LOSAT intent
  also starts with Settings open. Users can collapse it manually; selecting
  Run LOSAT again reopens it without changing the chosen mode or thread values.
  Fresh and reset Web state defaults Execution to `threaded`. Explicit saved
  `auto`, `serial`, or `threaded` choices remain authoritative on Load.
  The maintainer specified these command order, default, and disclosure
  requirements on `2026-09-14`; the existing execution modes and their support
  checks remain.
- May retire: none.
- Accepted residual risk: Increased computation time and memory from complete
  all-record comparisons. The maintainer explicitly accepted this cost as the
  original behavior. This does not permit silent truncation or hidden caps.
- Acceptance contracts: `OIC-005`, `OIC-006`, `OIC-007`, `OIC-013`, `OIC-015`.
- Decision source: The maintainer explicitly specified Cartesian Adjacent and
  complete all-record scopes, requested durable regression protection, accepted
  the computation/memory cost, and confirmed `May retire: none`, the Owner, and
  the decision date.
- Owner and decision date: `satoshikawato`, `2026-09-13`.

### PD-OI-019: Fresh Web multi-record layout defaults

- Concern key: `diagram-generation.web-multi-record-layout-defaults`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: In a fresh Web document and after Reset Settings, Linear
  defaults to Arrange in rows, grouping records from each source file into the
  same row while retaining independent per-record placement controls. Circular
  defaults to Multiple records in a single canvas (Multi-Record Canvas).
  Repeated record accessions retain distinct feature identities on that canvas.
  Explicitly saved layout choices take precedence on Session import; disabling
  either setting remains supported.
- Rationale: The maintainer explicitly requested these two defaults as part of
  restoring complete multi-record workflows.
- Must preserve: All loaded records; independent placement; explicit saved
  settings; the ability to choose separate rows or disable Circular shared
  canvas; comparison scope as specified in `PD-OI-018`.
- May retire: none. Both existing layout choices remain available.
- Accepted residual risk: Complete all-record comparisons retain the computation
  and memory cost accepted under `PD-OI-018`.
- Acceptance contracts: `OIC-004`, `OIC-006`, `OIC-015`.
- Decision source: Explicit maintainer follow-up requesting same-file Linear
  rows and Circular shared canvas as defaults in the same decision session.
- Owner and decision date: `satoshikawato`, `2026-09-13`.

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
| `OIC-015` | Linear discovery, per-record placement, actual comparison jobs, explicit request endpoints, Session replay, and SVG endpoints retain the complete selected record universe. Adjacent uses neighboring-row Cartesian products; Similarity and Collinear `all` retain complete directed evidence regardless of display rows. |

### OIC-015 required regression coverage

The normal automated PR gate must observe all of the following:

- Multi-record GenBank and GFF3/FASTA discovery retains every record with
  comparisons enabled and disabled; explicit selectors remain exact.
- A two-record source and a three-record source expose two file-input cards
  and five independent placements. A single two-record upload has one file
  card and two record controls. Assert both counts after Save/fresh Load,
  source replacement, source removal, and per-record row changes; filenames
  alone must not collapse distinct input sources. The neighboring 2-by-3 rows
  retain all six Pairwise record combinations through one source-file search
  job. The number of record pairs must not be reported as LOSAT job count.
  The record list initially shows only its single-line count; pointer and
  keyboard expansion expose every record without changing its state.
  Fresh and reset documents use same-file Linear rows and Circular shared
  canvas; explicit saved opt-outs survive Load.
- Similarity and Collinear `all` retain the complete directed record-pair
  matrix through source-file jobs, including self, same-row, and non-adjacent
  evidence. Two sources containing five records require four source jobs and
  cover 25 directed record pairs. A single-row
  layout retains analysis even though it has no between-row links.
- Explicit endpoints survive the typed decoder and reach the intended SVG
  records, including a numerically consecutive pair beside a multi-record row.
- Save, fresh Load, regeneration, and reordering preserve record identity,
  placement, pair mapping, and shared source resources. Repeated execution
  reuses only semantically equivalent cached searches.
- Every LOSAT mode exposes the existing execution and thread controls in
  Comparison Settings; mode switches and Session replay preserve their values.
  Verify Run LOSAT occupies the full top row, with No comparison and Upload
  BLAST TSV side by side below it. At desktop and narrow viewport widths,
  DOM order and actual keyboard Tab traversal match that visual order.
  Pointer and keyboard activation of Run LOSAT open Settings without moving
  focus away from the command; reopening it requires no mode change. Fresh
  and reset Execution is `threaded`, and saved explicit execution modes survive
  Load with active LOSAT settings immediately exposed.
- Changing only a record's drawing start or reverse-complement display issues
  no additional LOSATP search jobs, updates the rendered coordinates/orientation,
  and preserves that reuse through Save and fresh Load.
- Unbounded member selection retains more than five threshold-qualified hits.
  Explicit finite member limits still apply, and changing only that limit
  invalidates derived output while preserving raw search cache reuse.

These are observations of jobs, controls, requests, Sessions, and SVG results;
checking an `All` label or counting uploaded files is insufficient. Restoring
first-record truncation, zipped Adjacent pairing, or positional endpoint
coercion must fail the corresponding regression test. Source changes to these
boundaries require this coverage in the normal PR gate.

## Residual-risk boundary

Accepted residual risk never authorizes a security vulnerability, silent
scientific-output corruption, loss of a must-preserve effect, deterministic
Architecture Ratchet failure, undocumented unbounded performance regression,
cache reuse across different execution semantics, artifact provenance that
disagrees with actual execution, or failure of a required acceptance contract.
