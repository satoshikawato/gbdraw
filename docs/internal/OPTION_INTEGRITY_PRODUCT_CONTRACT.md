# Option Integrity Product Contract

Status: active Product authority

## Authority metadata

- Contract ID: `OIPC`
- Contract revision: `18`
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
- Additional approved decision IDs: `PD-OI-018`, `PD-OI-019`, `PD-OI-020`,
  `PD-OI-021`, `PD-OI-022`, `PD-OI-023`, `PD-OI-024`, `PD-OI-025`,
  `PD-OI-026`, `PD-OI-027`, `PD-OI-028`, `PD-OI-029`, `PD-OI-030`,
  `PD-OI-031`, `PD-OI-032`, `PD-OI-033`, `PD-OI-034`, `PD-OI-035`,
  `PD-OI-036`, `PD-OI-037`, `PD-OI-038`, and `PD-OI-039`
- Revision 3 addition: `PD-OI-018`, accepted by `satoshikawato` on
  `2026-09-13` after confirming the complete record/search outcome, no feature
  retirement, and the runtime/memory cost of complete comparisons. The initial
  approval and its date above continue to describe `PD-OI-001`–`PD-OI-017`.
  The maintainer subsequently specified the fresh defaults in `PD-OI-019`.
  On the same date, the maintainer replaced the default-five member cap in
  `PD-OI-004` with an unbounded default and clarified `PD-OI-018` to require
  source-file execution, reusable raw evidence after display transforms,
  compact record disclosures, and execution controls in LOSAT Settings.
- Revision 4 addition: `PD-OI-020`, requested by `satoshikawato` on
  `2026-09-14` to restore automatic thread allocation when Total threads changes
  and protect that behavior against regression.
- Revision 5 changes: `PD-OI-001`, `PD-OI-002`, `PD-OI-004`, and `PD-OI-018`
  are replaced for the maintainer's `2026-09-14` LOSATP follow-up. New
  `PD-OI-021` records optional Collinear self-search and inference;
  `PD-OI-022` records reuse after cancellation. The final clarification requires
  each mode's edited limits to return when that mode is selected again, not a
  reset on every switch. The maintainer explicitly requested these decisions
  and contracts be recorded to prevent regression. Earlier approvals above
  retain their original scope; this amendment records the new instructions.
- Revision 6 addition: `PD-OI-023`, selected by `satoshikawato` on
  `2026-09-15` through the complete `PRODUCT_DECISION` response for
  `protein-comparison.path-representation`, scenario revision `1`.
  Only PATH-B and the supplied preservation, retirement and risk terms are
  recorded. The earlier decisions retain their scope. Dependent runtime still
  requires this authority on its base; this amendment contains no runtime.
- Revision 7 addition: `PD-OI-024`, selected by `satoshikawato` on
  `2026-09-19` through the complete `PRODUCT_DECISION` response for
  `linear.definition-display`, scenario revision `1`. Only D1-A, D2-P,
  D3-A and the supplied preservation, retirement and risk terms are recorded.
  Earlier decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 8 change: `PD-OI-018` is replaced for scenario revision `3`,
  selected by `satoshikawato` on `2026-09-20` through the complete
  `PRODUCT_DECISION` response for
  `diagram-generation.linear-record-universe-and-search-scope`. The selected
  `LINEAR-FILE-ROW-BLOCK` outcome makes normal-layout File-card order and
  visual row order one operation, blocks File-card moves for custom layouts,
  and records the supplied preservation, retirement, and risk terms. Earlier
  decisions retain their scope. Dependent runtime requires this authority
  merged into its base; this amendment contains no runtime.
- Revision 9 addition: `PD-OI-025`, selected by `satoshikawato` on
  `2026-09-21` through the complete `PRODUCT_DECISION` response for
  `diagram-generation.linear-depth-source-scope-and-discoverability`, scenario
  revision `1`. The selected `FILE-BULK-WITH-RECORD-OVERRIDES` outcome exposes
  common Depth TSV assignment on each Linear File card while preserving sparse
  per-record bindings. Earlier decisions retain their scope. Dependent runtime
  requires this authority merged into its base; this amendment contains no
  runtime.
- Revision 10 additions: `PD-OI-026` through `PD-OI-031`, selected by
  `satoshikawato` on `2026-09-22` through six complete `PRODUCT_DECISION`
  responses for issue `#561`. These additions record only the supplied choices,
  preservation requirements, retirement permissions, and accepted residual
  risks for deterministic Similarity Group alignment. Earlier decisions retain
  their scope. Dependent runtime requires this authority merged into its base;
  this amendment contains no runtime.
- Revision 11 addition: `PD-OI-032`, selected by `satoshikawato` on
  `2026-09-22` through the complete `PRODUCT_DECISION` response for issue
  `#563`. This addition records the selected feature-popup record-rotation
  outcome, preservation requirements, lack of retirement permission, and
  accepted residual risk. Earlier decisions retain their scope. Dependent
  runtime requires this authority merged into its base; this amendment contains
  no runtime.
- Revision 12 additions: `PD-OI-033` through `PD-OI-035`, approved by
  `satoshikawato` on `2026-09-24` through explicit approval of the exact three
  `PRODUCT_DECISION` texts presented for issue `#581`. The receipts below are
  the full scope of these additions. Dependent runtime requires this authority
  merged into its base; this amendment contains no runtime.
- Revision 13 changes: `PD-OI-026`, `PD-OI-027`, `PD-OI-031`, and `PD-OI-034`
  are replaced for scenario revision `2`, selected by `satoshikawato` on
  `2026-09-24` through four complete `PRODUCT_DECISION` responses for issue
  `#586`. The receipts below are the full scope of these replacements. Earlier
  unaffected decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 14 changes: `PD-OI-026`, `PD-OI-027`, `PD-OI-031`, and `PD-OI-034`
  are replaced for scenario revision `3`. `satoshikawato` explicitly
  approved the four complete Choice A `PRODUCT_DECISION` texts on
  `2026-09-25` for the issue `#586` follow-up. The receipts below are the
  full scope of these replacements. Other decisions retain their scope.
  Dependent runtime requires this authority merged into its base; this
  amendment contains no runtime.
- Revision 15 change: `PD-OI-035` is replaced for scenario revision `2`.
  `satoshikawato` supplied the complete `A / RETAIN_MOBILE_PALETTE_COVERAGE`
  `PRODUCT_DECISION` response on `2026-09-25`. The receipt below is the full
  scope of this replacement. Other decisions retain their scope. This
  amendment contains no runtime.
- Revision 16 changes: `PD-OI-027`, `PD-OI-028`, `PD-OI-029`, `PD-OI-031`,
  and `PD-OI-034` are replaced for scenario revisions `4`, `2`, `2`, `4`,
  and `4`. `satoshikawato` explicitly approved the five complete
  `PRODUCT_DECISION` texts as written on `2026-09-26` for the record-owned
  orientation follow-up to issue `#586`. The receipts below are the full
  scope of these replacements. `PD-OI-026`, `PD-OI-030`, `PD-OI-035`, and
  other decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 17 changes: `satoshikawato` signed five complete
  `PRODUCT_DECISION` responses for issue `#602` on `2026-09-26`:
  01=B and 02–05=A. `PD-OI-024` scenario revision `2` replaces only the
  D1 Web fresh/reset default and retains D2-P/D3-A. New `PD-OI-036` through
  `PD-OI-039` serialize the other four receipts. `PD-OI-035` scenario
  revision `3` retains its independent canvas-interaction guarantees and
  supersedes only the scenario-2 mobile coverage exception and close-review
  continuation, using the signed review-presentation receipt in `PD-OI-039`.
  The five receipts below preserve all nine supplied fields exactly. Earlier
  unaffected outcomes retain their scope. This amendment contains no runtime;
  dependent runtime requires all five outcomes and the limited supersession
  merged into its base.
- Revision 18 changes: `PD-OI-027`, `PD-OI-029`, `PD-OI-031`, and
  `PD-OI-034` are replaced for scenario revisions `5`, `3`, `5`, and `5`.
  `satoshikawato` explicitly approved the four complete issue `#598` receipts
  on `2026-09-26`. They define exclusive displayed-direction modes and
  selectable Reset scope. `PD-OI-026`, `PD-OI-028`, `PD-OI-030`,
  `PD-OI-032`, `PD-OI-033`, `PD-OI-035`, and other decisions retain their
  scope. This amendment contains no runtime; dependent implementation must
  have this authority merged into its base. Issue `#598` BUG-17 is excluded
  from implementation at the owner's instruction, with no circular-rotation
  authority change.
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

### PD-OI-001: LOSATP raw-search limit fresh defaults

- Concern key: `diagram-generation.losatp-candidate-limit-default`
- Scenario revision: `2`
- Supersedes: `PD-OI-001`, scenario revision `1`.
- Status: `ACCEPTED`
- Normative outcome: Web exposes the raw LOSATP `max_target_seqs` limit as
  **Max target seqs**. Fresh and reset Collinear starts at `5`; Similarity
  groups starts unbounded (`None`). Pairwise, CLI, and Python omission defaults
  retain their existing meanings. A blank Web value explicitly means unbounded,
  including in Collinear. No hidden cap substitutes for an unbounded request.
- Rationale: The maintainer requested exposure of the actual raw-search limit
  and distinct Collinear/Similarity defaults, with regression protection.
- Must preserve: Explicit finite and unbounded values; truthful requested and
  effective metadata; raw-cache identity; cancellation and errors; saved values.
- May retire: The unbounded fresh Web Collinear default. Unbounded search itself
  remains available.
- Accepted residual risk: The existing unbounded-work cost remains. This
  amendment adds no performance guarantee or additional risk waiver.
- Acceptance contracts: `OIC-001`, `OIC-005`, `OIC-013`, `OIC-017`.
- Decision source: Explicit maintainer instruction to expose Max target seqs,
  default Collinear to `5`, and leave Similarity groups unbounded.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-002: LOSATP mode-specific limit retention and GUI placement

- Concern key: `diagram-generation.losatp-candidate-limit-scope`
- Scenario revision: `2`
- Supersedes: `PD-OI-002`, scenario revision `1`.
- Status: `ACCEPTED`
- Normative outcome: **Max target seqs** is directly editable in LOSATP
  Settings, alongside **Member hits per protein** when member selection applies.
  Each mode remembers its own raw and member limits. The first visit to
  Collinear uses `5`/`5`; the first visit to Similarity groups uses unbounded/
  unbounded. After editing either mode, switching away and back restores its
  edited values. This works repeatedly in both directions and preserves blanks.
  Switching modes neither resets the returning mode to defaults nor copies the
  departing mode's limits into it. Save and fresh Load preserve active and
  inactive limits; Reset Settings clears the drafts to their declared defaults.
- Rationale: The maintainer clarified the required sequence as edit Similarity,
  edit Collinear, return to Similarity's values, then return to Collinear's values.
- Must preserve: Independently edited limits; existing Session values;
  Pairwise display-limit independence; discoverability and keyboard operation.
- May retire: A single raw-limit value shared across presentations and
  Advanced-only placement. There must not be duplicate controls for one active
  raw-search setting.
- Accepted residual risk: No additional risk waiver was supplied. Search
  changes caused by selecting a different saved raw limit remain visible;
  matching raw evidence remains reusable.
- Acceptance contracts: `OIC-001`, `OIC-002`, `OIC-006`, `OIC-017`.
- Decision source: Explicit maintainer correction to remember each mode's
  settings, superseding the earlier same-session request to reset on every switch.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

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
- Scenario revision: `3`
- Supersedes: `PD-OI-004`, scenario revision `2`.
- Status: `ACCEPTED`
- Normative outcome: Fresh and reset Web **Member hits per protein** starts at
  `5` in Collinear and unbounded (`None`) in Similarity groups. CLI/Python
  omission semantics remain unchanged. A blank Web control always means
  unbounded, including in Collinear; an explicit positive integer limits the
  distinct directional subject candidates retained after result filtering.
  Member selection remains independent of raw Max target seqs and Pairwise
  display matches. Collinear consumes it with inference both OFF and ON.
  Per-mode retention follows `PD-OI-002`. Request, helper execution, provenance,
  and Session replay distinguish unbounded and finite choices.
- Rationale: The maintainer requested the same mode defaults and retention for
  both exposed limits, while preserving their distinct scientific roles.
- Must preserve: Raw-search reuse when only member hits changes; correct
  derived invalidation; explicit saved values; threshold filtering.
- May retire: The unbounded fresh Web Collinear member default. Blank-to-five
  coercion and aliasing this field to the raw or Pairwise limits remain prohibited.
- Accepted residual risk: Member selection can change blocks or groups and
  remains visible in provenance. The existing unbounded computation/memory
  allowance is unchanged; no additional waiver is recorded.
- Acceptance contracts: `OIC-002`, `OIC-004`, `OIC-005`, `OIC-006`, `OIC-017`.
- Decision source: Explicit maintainer instruction that Member hits per protein
  receive the same mode defaults and remembered-value behavior as Max target seqs.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

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
- Scenario revision: `3`
- Supersedes: `PD-OI-018`, scenario revision `2`.
- Status: `ACCEPTED`
- Selected option: `LINEAR-FILE-ROW-BLOCK`
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
  4. Similarity groups and Collinear with scope `all` retain every ordered
     between-record comparison, including same-row and non-adjacent pairs.
     Similarity groups always includes within-record evidence. Collinear
     includes it only when **Infer orthogroups with self-comparisons** is ON,
     as specified by `PD-OI-021`; OFF excludes it. Display placement, including
     a single occupied row, does not restrict the selected search scope.
  5. Explicit comparison endpoints stay explicit through decoding and
     rendering, including endpoints whose numeric indices are consecutive.
  6. A File card number and its up/down controls represent both input-source
     order and visual row order in a normal Linear layout. Moving a File moves
     every record belonging to that source as one block. All records in the
     moved File remain in the same File-owned row, and File-owned rows follow
     File-card order. The move is one atomic, undoable draft operation.
     Record identity, source association, selector, crop, reverse complement,
     definition, subtitle, depth binding, feature state, and within-File
     record order remain attached to the same record.

     A normal layout has exactly one occupied row per File and no row shared
     by records from another File. When one File spans multiple rows or one row
     is shared by multiple Files, the File move is unavailable. The application
     explains that custom Record Layout controls visual placement and directs
     the user to Record Layout. A blocked move makes no partial change to File
     order, row assignments, comparisons, cache metadata, or the current
     Result.

     Explicit comparison endpoints remain attached by stable record UID and
     are reindexed without changing their biological endpoints. Adjacent
     comparison is resolved from the new occupied-row adjacency after the File
     move. Compatible raw LOSAT evidence remains reusable; incompatible
     derived comparison artifacts are recalculated.

     Save, fresh Load, regeneration, keyboard operation, and Session replay
     preserve the resulting File and row order. The last successful Result
     remains unchanged until Generate succeeds. Failed, canceled, superseded,
     or stale Generate does not replace it. Explicit selection/cropping,
     comparison omission, and imported read-only intent remain supported.
     `OIPC-C07` governs failed, canceled, and stale work. No new File-order
     state, Session field, compatibility migration, request schema, Worker
     protocol, or rendering path is introduced.
  7. Shared source bytes remain shared. Complete comparisons may require more
     jobs, but no hidden record or pair cap is permitted.
     LOSAT batches compatible records by input source file. A source job searches
     multi-sequence FASTA inputs containing the selected records, then routes
     hits to their record endpoints. Two sources
     in all-vs-all with within-record evidence enabled require four directed
     source jobs, including self and reverse searches; eight total records must
     not become 64 LOSAT invocations merely because records are expanded.
     TLOSATX records with different explicit genetic codes use compatible
     subsets within each source because each invocation accepts one query and
     one subject genetic code; records sharing those settings remain batched.
     Adjacent and explicit selections restrict retained record-pair evidence.
     When Collinear inference is OFF, batching must not search any record
     against itself, even inside a multi-record source. Comparisons within one
     source may therefore require separate jobs; compatible between-source
     comparisons remain batched. Search arguments, source contents,
     and the actual searched database scope are part of raw-cache identity.
     Progress reports actual source jobs separately from biological record
     pairs. File batching must preserve cancellation and exact hit routing.
     Drawing-start changes and reverse-complement display reuse raw LOSATP
     results whenever source contents, selected biological regions, and search
     settings are unchanged. Only display coordinates and derived presentation
     are updated; display transforms do not become raw-search identity.
- Rationale: Prevent recurrence of incomplete-record search and per-record
  placement regressions reported by the maintainer. With Arrange in rows
  enabled by default, changing File order while retaining the File's previous
  absolute row produces a visible result that contradicts the explicit reorder
  action. Source order and visual row order therefore change as one operation.
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
- May retire: Mandatory Collinear within-record evidence when inference is OFF;
  preservation of a moved File's absolute numeric row in a normal
  one-File-per-row layout; File-card movement that changes input source order
  while leaving visible row order unchanged; and any acceptance test that
  disables Arrange in rows before proving the primary default-layout
  File-reorder outcome. Other record coverage, placement, and comparison
  capabilities remain supported.
- Accepted residual risk: Increased computation time and memory from complete
  all-record comparisons. The maintainer explicitly accepted this cost as the
  original behavior. This does not permit silent truncation or hidden caps.
  Moving a File can change the derived Adjacent comparison pair set and can
  require comparison recomputation on the next Generate. Custom layouts require
  the user to edit or normalize Record Layout before File-card reordering is
  available. No silent loss of custom placement or comparison intent is
  accepted.
- Acceptance contracts: `OIC-005`, `OIC-006`, `OIC-007`, `OIC-013`, `OIC-015`.
- Original decision source: The maintainer explicitly specified Cartesian Adjacent and
  complete all-record scopes, requested durable regression protection, accepted
  the computation/memory cost, and confirmed `May retire: none`, the Owner, and
  the decision date.
- Amendment source: The maintainer requested Collinear self-comparison and
  orthogroup inference be optional and default OFF on `2026-09-14`.
- File-row amendment source: The maintainer supplied the complete
  `PRODUCT_DECISION` response reproduced below on `2026-09-20`.
- Owner and decision date: `satoshikawato`, `2026-09-20`.

```json
{
  "concern": "diagram-generation.linear-record-universe-and-search-scope",
  "scenarioRevision": 3,
  "choice": "LINEAR-FILE-ROW-BLOCK",
  "rationale": "A File card number and its up/down controls communicate the source's visual order in a Linear diagram. With Arrange in rows enabled by default, changing File order while retaining the File's previous absolute row number produces a visible result that contradicts the user's explicit reorder action. File reordering must therefore update source order and visual row order as one operation.",
  "mustPreserve": [
    "Each uploaded GenBank file or paired GFF3/FASTA source remains one File card.",
    "Moving a File moves every record belonging to that source as one block.",
    "In the normal layout, every record from the moved File remains in the same File-owned row, and File rows follow File-card order.",
    "Record identity, source association, selector, crop, reverse complement, definition, subtitle, depth binding, feature state, and within-File record order remain attached to the same record.",
    "Explicit comparison endpoints remain attached by stable record UID and are reindexed without changing their biological endpoints.",
    "Compatible raw LOSAT evidence remains reusable; incompatible derived comparison artifacts are recalculated.",
    "Adjacent comparison is resolved from the new occupied-row adjacency after the File move.",
    "The move is one atomic, undoable draft operation.",
    "The last successful Result remains unchanged until Generate succeeds. Failed, canceled, superseded, or stale Generate does not replace it.",
    "Save, fresh Load, regeneration, keyboard operation, and Session replay preserve the resulting File and row order.",
    "No new File-order state, Session field, compatibility migration, request schema, Worker protocol, or rendering path is introduced."
  ],
  "customLayoutRule": "A normal layout means that every File occupies exactly one row and no row is shared by records from another File. If one File spans multiple rows or one row is shared by multiple Files, the File move is unavailable. The application must explain that the custom Record Layout controls visual placement and direct the user to Record Layout. A blocked move makes no partial change to File order, row assignments, comparisons, cache metadata, or the current Result.",
  "mayRetire": [
    "Preservation of a moved File's absolute numeric row in a normal one-File-per-row layout.",
    "File-card movement that changes input source order while leaving the visible row order unchanged.",
    "Any acceptance test that disables Arrange in rows before proving the primary default-layout File-reorder outcome."
  ],
  "acceptedResidualRisk": "Moving a File can change the derived Adjacent comparison pair set and can therefore require comparison recomputation on the next Generate. Custom layouts require the user to edit or normalize Record Layout before File-card reordering becomes available. No silent loss of custom placement or comparison intent is accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-20"
}
```

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


### PD-OI-020: LOSAT total-thread allocation and control agreement

- Concern key: `diagram-generation.losat-total-thread-allocation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome:
  1. Total threads sets the effective total budget. Safe uses half the available
     hardware threads, with a minimum of one; Available uses the hardware
     count. Explicit numeric budgets cannot exceed the available hardware.
  2. Auto Parallel runs and Auto Threads per run recalculate when the budget
     or pending source-job count changes. Auto Threads per run distributes the
     budget across the selected simultaneous runs; multiple jobs do not impose
     a hidden fixed two-thread limit. With four pending LOSATP source jobs and
     both controls on Auto, numeric totals 32, 16, and 2 produce respectively
     4 runs × 8 threads, 4 × 4, and 2 × 1, when hardware permits those totals.
  3. The displayed plan and actual execution use the same allocation rules.
     Before Generate, the controls may use the estimated source-job count;
     cached jobs are omitted from actual execution. Actual run information
     reports the effective execution allocation. Simultaneous runs multiplied
     by threads per run never exceeds the effective total budget.
  4. Explicit Parallel runs and Threads per run choices remain independent
     editable intent. Auto adjusts around the manual choice. A temporary
     effective clamp does not replace the saved manual value or Auto with its
     computed value; any requested/effective difference is visible. Increasing
     the budget makes a preserved manual choice effective again when feasible.
     Save and fresh Load preserve these choices.
- Rationale: Restore the maintainer-requested linkage between Total threads,
  simultaneous runs, and threads per run, and prevent the controls from
  promising an allocation that execution does not use.
- Must preserve: Existing execution modes, browser support checks, fixed
  single-thread-per-run constraints for LOSATN and TLOSATX, source-file job
  batching, cancellation, and explicit saved choices. Scheduling changes do
  not alter search parameters, member-hit limits, or raw-search identity.
- May retire: none.
- Accepted residual risk: The existing computation and memory allowance in
  `PD-OI-018` remains; this decision grants no exception to the selected total
  thread budget and does not promise a fixed speedup for every workload.
- Acceptance contracts: `OIC-005`, `OIC-013`, `OIC-016`.
- Decision source: Explicit maintainer request to restore or implement linked
  automatic allocation and record this regression as a Contract.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-021: Optional Collinear self-search and orthogroup inference

- Concern key: `diagram-generation.collinear-optional-orthogroup-inference`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Web Collinear exposes **Infer orthogroups with
  self-comparisons** as a checkbox, OFF in fresh and reset state. OFF submits
  only between-record searches and builds blocks without orthogroup or paralog
  inference. This exclusion applies even when query and subject records share
  one source file. The chosen anchor, unit, member-hit, threshold, scope, and
  block settings still apply. ON enables within-record evidence and the existing
  orthogroup/paralog inference. Similarity groups retains its inference behavior.
  Session save/load preserves explicit ON/OFF. Older Collinear Sessions that
  omit the choice preserve their historical ON behavior; their saved preview
  and explicit limits are not reinterpreted as fresh defaults. Typed request,
  actual helper call, derived identity, and provenance agree on the choice.
  Matching raw evidence can be reused; inference changes cannot reuse an
  incompatible derived result.
- Rationale: The maintainer identified mandatory self-search and paralog-aware
  grouping during Collinear work and requested a default-OFF checkbox.
- Must preserve: Direct block construction when OFF; existing inference when
  ON; complete cross-record scope; explicit saved choices; raw/derived cache
  correctness. An inactive paralog-link setting does not claim an effect in OFF.
- May retire: Mandatory self-search and orthogroup inference in Web Collinear.
  Both remain available through the checkbox.
- Accepted residual risk: No additional risk waiver was supplied. ON retains
  the existing complete-evidence computation and memory allowance in `PD-OI-018`.
- Acceptance contracts: `OIC-003`, `OIC-005`, `OIC-006`, `OIC-015`, `OIC-018`.
- Decision source: Explicit maintainer request to make Collinear self-comparison
  and orthogroup inference optional with default `false`/OFF.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-022: Reuse completed LOSAT searches after cancellation

- Concern key: `diagram-generation.completed-losat-search-retry`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: After LOSAT raw search completes, canceling downstream
  Collinear work does not require the same raw search to run again on Generate.
  A member-only change recomputes the dependent analysis and reuses matching
  completed raw evidence. Cache reuse still validates the actual inputs and
  raw search settings; changing Max target seqs can require a new search.
  Clear Cache invalidates retained evidence. Session/History replacement must
  not revive unrelated retry data. Failure or cancellation preserves the last
  successful Result and committed request under `OIPC-C07`.
- Rationale: The maintainer reported a completed LOSATP search restarting after
  canceling Collinear work and changing only member hits, and requested a fix
  with durable regression protection.
- Must preserve: Correct raw and derived identities; deterministic cancellation;
  explicit clearing; last successful output. Partial raw batches are not falsely
  admitted as completed searches.
- May retire: Unnecessary raw reruns caused solely by downstream cancellation
  or a member-limit edit.
- Accepted residual risk: No new persistence or cross-reload guarantee is
  introduced; no additional risk waiver was supplied.
- Acceptance contracts: `OIC-002`, `OIC-005`, `OIC-013`, `OIC-019`.
- Decision source: Explicit maintainer request to fix the reported cancellation
  cache behavior and add Product Decisions/Contracts against regression.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-023: Lossless protein paths during normal generation and saving

- Concern key: `protein-comparison.path-representation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `PATH-B / lossless-graph`
- Normative outcome: 通常の生成・保存では、全経路を再現できるgraphを保持する。
  図と解析情報、全経路の内容・順序・ID・shared情報を維持し、対応済み旧ファイルを
  読み込めることと、明示的な旧tuple形式での全経路取得を維持する。
  通常のAPI戻り値と保存形式が常に全経路配列を含む仕様は廃止してよい。
- Decision source: The complete maintainer response to S02 Decision Pack
  revision `1`, prepared at `c46e55d14b6b8fae9e6778a7f04183a4d321185f`.
  The receipt below preserves every supplied field without translating or
  extending its rationale, retirement, risk, owner or date. It is a reviewable
  serialization within this existing authority document, not a new decision
  store or evaluator. Implementation names and schema allocations in the S02
  design remain engineering proposals; this receipt does not freeze them.

```json
{
  "concern": "protein-comparison.path-representation",
  "scenarioRevision": 1,
  "choice": "PATH-B / lossless-graph",
  "rationale": "解析情報を維持しながら、通常の生成・保存での経路展開コストを減らしたい。",
  "mustPreserve": "図と解析情報、全経路の内容・順序・ID・shared情報、対応済み旧ファイルの読み込み、明示的な旧tuple形式での全経路取得。",
  "mayRetire": "通常のAPI戻り値と保存形式が、常に全経路配列を含む仕様。",
  "acceptedResidualRisk": "旧API依存コードの修正、新形式を旧バージョンで読めないこと、明示的な全量取得には大きな時間・メモリが必要になり得ること。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-15"
}
```

### PD-OI-024: Linear Definition alignment and automatic Replicon visibility

- Concern key: `linear.definition-display`
- Scenario revision: `2`
- Supersedes: `PD-OI-024`, scenario revision `1`, only D1's Web fresh/reset
  Lock=false default; D2-P and D3-A remain required in full.
- Status: `ACCEPTED`
- Selected outcomes: `D1-B-WEB-LOCKED-FRESH`, `D2-P`, `D3-A`
- Normative outcome:
  1. **D1-B-WEB-LOCKED-FRESH:** Web fresh/reset selects Lock Definition
     Column=true. Lock=true retains the common left edge and configured gap;
     explicit Lock=false retains common-width centering and row following.
     Single, shared, and mixed rows, the existing accepted `text_anchor`
     domain, explicit saved values, supported old omission meanings, the
     saved Result on Load, and CLI/Python omitted defaults remain supported.
     Linear Layout always explains ON/OFF and application on Generate.
  2. **D2-P — 保存値を保持:** 保存済みSubtitleは、自動／手入力を推測したり、
     Replicon名と文字列が一致したりすることを理由に削除しない。
     読み込みだけでは保存Resultを変えず、Generateで新しい表示契約を適用する。
     不要なSubtitleは利用者が明示的にクリアし、ファイル既定値へ戻る既存の
     継承規則を維持する。手入力Subtitleと行共通・レコード固有ラベルの区別を保つ。
  3. **D3-A — OrganelleもReplicon行で制御:** chromosome、plasmid、organelle
     由来の自動名をShow Repliconで制御し、対象名はオンで一つ、オフでゼロとする。
     候補が競合するときはchromosome→plasmid→organelleの順で一つを選ぶ。
     organelleの表記は既存の自動Subtitle表記を引き継ぐ。
     Web・CLI・Pythonの共通描画に適用し、Show Repliconの既定値falseを維持する。
     自動名のオン／オフは手入力Subtitleの表示を変更しない。
- Decision source: The complete issue `#602` response for
  `linear.definition-display`, scenario revision `2`, signed in full by
  `satoshikawato` on `2026-09-26`, reproduced below without changing any
  supplied field. The retained D2-P/D3-A clauses above are unchanged from
  scenario revision `1`. This is serialization in the existing static
  authority document; dependent runtime requires it merged into its base.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-023`.

```json
{
  "concern": "linear.definition-display",
  "scenarioRevision": 2,
  "choice": "A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A",
  "rationale": "Webで新しく作るLinear比較図ではDefinitionを共通左列にそろえ、行のalignやoffset後も名前を比較しやすくする。",
  "mustPreserve": "Lock=trueの共通左端とconfigured gap、明示Lock=falseの共通幅中央とrow追従、単一/共有/混在行、既存text_anchorの受入範囲、保存Sessionの明示値と対応済み旧省略意味、読込時の保存Result、CLI/Python省略default。PD-OI-024のD2-Pの保存/手入力Subtitleと継承・ラベル区別、およびD3-AのReplicon/Organelle選択順・独立制御・既定falseをすべて維持する。Linear LayoutでON/OFFの違いとGenerate適用を常時説明する。",
  "mayRetire": "D1-Aのうち、Web fresh/resetがLock=falseを初期値として選ぶ部分だけ。OFFの明示操作とCLI/Python既定は退役しない。",
  "acceptedResidualRisk": "新しいWeb図のDefinition外観が従来のfresh図と変わり、Definitionがrowに追従しなくなる。利用者はOFFを選べ、既存Sessionの値と保存Resultは勝手に変更しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-025: Linear Depth source scope and discoverability

- Concern key: `diagram-generation.linear-depth-source-scope-and-discoverability`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `FILE-BULK-WITH-RECORD-OVERRIDES`
- Normative outcome:
  1. Each Linear File card exposes the common Depth TSV assignment without
     requiring the user to expand its record list or any record options.
     Applying or clearing a File-level value updates every record binding for
     that File and logical series as one undoable operation.
  2. Sparse per-record Depth bindings remain editable. A File card distinguishes
     empty, common, and mixed record bindings. Applying a File-level value to a
     mixed series replaces every record binding in that File and series; the UI
     discloses this effect before the action.
  3. Logical-series settings shared across records are presented once rather
     than repeated inside every record card. Per-record controls expose only the
     record-specific source assignment. A single-record File does not receive a
     duplicate record-level uploader for the same binding.
  4. The canonical state remains the record-major Depth matrix. Null cells,
     logical series indexes, per-record overrides, source isolation, one-step
     Undo/Redo, Session round trips, regeneration, and canonical render-request
     semantics remain supported. The outcome introduces no new render path or
     requirement for a new persisted Depth-default field.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-21` for issue `#554`, reproduced below. The
  receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-020`.

```json
{
  "concern": "diagram-generation.linear-depth-source-scope-and-discoverability",
  "scenarioRevision": 1,
  "choice": "FILE-BULK-WITH-RECORD-OVERRIDES",
  "rationale": "Common Depth TSV input should be available once on the File card, while supported sparse per-record bindings remain editable.",
  "mustPreserve": [
    "Record-major Depth matrices",
    "Null cells and logical series indexes",
    "Per-record overrides",
    "Source isolation",
    "One-step Undo/Redo",
    "Session round trips",
    "Regeneration",
    "Existing canonical request semantics"
  ],
  "mayRetire": [
    "Duplicated global Depth settings inside every record card",
    "The need to expand every record before finding Depth input"
  ],
  "acceptedResidualRisk": "Applying or clearing a File-level value replaces or clears every record binding for that File and series; the UI must disclose this and Undo must restore the previous matrix.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-21"
}
```

### PD-OI-026: Deterministic Similarity Group anchor resolution with optional review

- Concern key: `diagram-generation.similarity-alignment.anchor-resolution`
- Scenario revision: `3`
- Supersedes: `PD-OI-026`, scenario revision `2` (`A / AUTO_PRESELECTED_SUGGESTIONS`).
- Status: `ACCEPTED`
- Selected outcome: `A / DETERMINISTIC_RESOLUTION_WITH_OPTIONAL_REVIEW`
- Normative outcome: The exact selected reference and stable biological
  identities govern independent target resolution. Python retains
  only-usable-candidate and unique-direct-RBH automatic resolution and
  disclosed transient recommendations for ambiguity. A fully resolved default
  alignment may commit without opening review. An explicit review still
  exposes candidate reasons, replacement, and Skip before commitment;
  ambiguous suggestions still require review. Missing, unusable, and skipped
  targets remain unchanged, and resolution is independent of viewport,
  geometry, confidence, supporting-edge count, and multi-hop evidence.
- Decision source: The complete four Choice A `PRODUCT_DECISION` texts
  explicitly approved by `satoshikawato` on `2026-09-25` for the issue
  `#586` follow-up, with this concern's receipt reproduced below. This
  serialization adds no terms to the approved receipt and becomes runtime
  authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.anchor-resolution",
  "scenarioRevision": 3,
  "choice": "A / DETERMINISTIC_RESOLUTION_WITH_OPTIONAL_REVIEW",
  "rationale": "Python should keep its deterministic automatic resolutions and disclosed ambiguity recommendations, while a fully resolved default alignment can commit without opening candidate review unless the user requests it.",
  "mustPreserve": "The exact selected reference; stable record and biological-feature identity; only-usable-candidate and unique-direct-RBH automatic resolution; independent displayed-record treatment; unchanged missing, unusable, and skipped records; visible recommendation reasons and local replacement or Skip in every opened review; final shared Python validation; and independence from viewport, scroll, ribbon geometry, confidence score, supporting-edge count, and multi-hop evidence. Unique representative and deterministic candidate 1 remain transient recommendations, not automatic resolutions.",
  "mayRetire": "The requirement that replacement or Skip be presented before every Python-resolved default alignment. The explicit review action must still present those choices before commitment.",
  "acceptedResidualRisk": "A resolved default plan may commit without per-target inspection. The named explicit review action and reliable Undo/Reset must remain available; ambiguous suggestions still require the full disclosed review.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-25"
}
```

### PD-OI-027: Explicit display directions with record-owned orientation

- Concern key: `diagram-generation.similarity-alignment.transform-semantics`
- Scenario revision: `5`
- Supersedes: `PD-OI-027`, scenario revision `4` (`A / RECORD_OWNED_ORIENTATION_WITH_REVIEW_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / EXPLICIT_DISPLAY_DIRECTION_MODES`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.transform-semantics",
  "scenarioRevision": 5,
  "choice": "A / EXPLICIT_DISPLAY_DIRECTION_MODES",
  "rationale": "Users should choose the final displayed direction of selected alignment features instead of making all targets follow a potentially minority-direction reference. Reference identity defines positioning, while record direction remains independent record state.",
  "mustPreserve": "Default Keep current directions; one exclusive Keep, All selected features right-facing, All selected features left-facing or Custom direction mode; Custom per-record Keep/right/left choices; exact reference identity and selected target anchors; bulk direction scope including the reference and only known-strand selected anchors; explicit unchanged reasons for unknown directions and unchanged missing, unusable and skipped records; record-wide absolute orientation updates without editing biological source strands; a fixed pre-Align canvas x of the reference feature center, adjusted record placement as necessary, unchanged vertical placement and exact idempotent selected-anchor center alignment; per-record before/after arrow previews and truthful scope coverage; readable text and feature/label/annotation/ribbon geometry in the same record transform; one atomic validated orientations/placement/plan/Result commit; orientation-independent plans, ordinary Reverse, accepted Reset scope and complete artifact Undo/Redo. Changed final validated directions or reference-center placement require a refreshed preview and another Apply before commitment. Capture actual direction deltas including the reference if a separately authorized reset receipt is enabled.",
  "mayRetire": "The rule that alignment always preserves the reference record's direction and left-edge position; the single reference-relative Match checkbox; the restriction that target exceptions require post-Apply sidebar Reverse. Do not retire exact reference selection, its fixed anchor-center position, the default Keep mode or ordinary Reverse.",
  "acceptedResidualRisk": "Users can confuse display arrows with biological strand annotations or interpret all as every feature in a source. Show selected-anchor scope, reference participation, record names, before/after arrows and unknown exclusions; preserve sources and preview reference left-edge movement while its feature center stays fixed. Custom increases review density and requires keyboard/390 px acceptance. Majority inference, guessed unknown directions, source annotation changes, hidden flips and separate orientation/render owners are not accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-028: Orientation-independent active alignment plan

- Concern key: `diagram-generation.similarity-alignment.plan-lifecycle`
- Scenario revision: `2`
- Supersedes: `PD-OI-028`, scenario revision `1` (`A / PERSISTED_ACTIVE_ALIGNMENT_PLAN`).
- Status: `ACCEPTED`
- Selected outcome: `A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: The complete five `PRODUCT_DECISION` texts for the
  record-owned orientation follow-up to issue `#586`, explicitly approved
  as written by `satoshikawato` on `2026-09-26`, with this concern's
  receipt reproduced below. This serialization adds no terms to the
  approved receipt and becomes runtime authority only after merge into
  the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.plan-lifecycle",
  "scenarioRevision": 2,
  "choice": "A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN",
  "rationale": "An active plan identifies anchors, not directions. A manual orientation change leaves every anchor valid, so keeping the plan spares the user a second alignment after reversing a record.",
  "mustPreserve": "Survival of the active plan across regeneration after style, label, and canvas-size changes, and across record reorder when stable record identities remain; survival across manual record orientation changes, with the same anchors aligned in the new orientation when the diagram is next generated; explicit clearing and notification after manual record movement, source replacement, crop changes, or record-selector changes; validation before regeneration; the last successful Result while a stale plan is repaired; and explicit reselect, Skip, or Clear actions without automatic anchor substitution.",
  "mayRetire": "Clearing the active plan when the user manually changes a record's orientation.",
  "acceptedResidualRisk": "After a manual orientation change, the reversed record moves horizontally on the next generation so that its anchor stays aligned. A user who wanted the previous position uses Undo or Reset Align.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-029: Selectable position and alignment-direction Reset with artifact history

- Concern key: `diagram-generation.similarity-alignment.reset-and-history`
- Scenario revision: `3`
- Supersedes: `PD-OI-029`, scenario revision `2` (`A / IMMEDIATE_PREALIGN_POSITION_BASELINE`).
- Status: `ACCEPTED`
- Selected outcome: `A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.reset-and-history",
  "scenarioRevision": 3,
  "choice": "A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE",
  "rationale": "Alignment can change both placement and record direction, so users must be able to choose whether Reset removes only positioning or also restores the directions actually changed by the latest Align. The restoration scope must be explicit without turning alignment plans into direction owners or replaying unrelated edits.",
  "mustPreserve": "An explicit default Reset positions command that restores positions immediately before the latest successful Align, clears its active plan and keeps all current directions; an additional Reset positions and alignment direction changes command that uses the same position baseline and restores absolute before-Align direction only for records whose direction actually changed in that Align; unchanged direction for every record not reversed by that Align, including the reference when unchanged and all later manual edits on unaffected records; include a reference record in restoration only if an independently authorized alignment direction choice actually reversed it; visible target names, count, current and restored directions, and disclosure that later manual direction edits on restoration targets are replaced by the combined command; replacement of the reset receipt by each new successful Align, with no first-Align or source-orientation fallback; an orientation-independent plan and ordinary record-owned current orientation; source-bound validated restoration information captured from actual successful before/after states, retained across style regeneration and stable reorder, saved and freshly loaded with new Sessions, and cleared atomically with plan invalidation or successful Reset. Missing old restoration information leaves positions Reset available and direction restoration unavailable with an explicit reason, never guessed; an empty modern delta means no Align direction changes. Both Reset scopes consume the active plan and receipt, preserve unrelated settings and pending form edits, and use one canonical artifact transaction. Undo/Redo restores or reapplies complete artifacts including directions, plan and receipt; failed, canceled, stale, superseded, preview-readiness or History-finalization work creates no committed history and preserves the prior artifact and receipt. Original sources, biological identities, target-external settings, readable text, record-consistent feature/label/ribbon geometry, existing canonical request, rendering, Worker, sanitizer/admission and History owners, compatible comparison reuse and zero additional Reset LOSAT jobs remain intact. Restoration receipts never select rendering orientation.",
  "mayRetire": "The rule that Reset always retains directions changed by alignment and that their restoration is available only through ordinary Undo or manual Reverse. Do not retire the position-only choice or either existing recovery workflow.",
  "acceptedResidualRisk": "The combined scope intentionally replaces later manual direction edits on records actually reversed by the latest Align. Bound this to a visible target and before-direction preview, an explicit position-only alternative and one-operation Undo. A positions-only Reset consumes the same active plan and receipt; switching to combined afterward requires Undo of that Reset first. Sessions without historical direction evidence cannot restore it and must show that limit while retaining positions Reset. A compact Session/History receipt adds bounded state-maintenance cost; require source/plan binding, atomic lifecycle, old-information/no-op/re-Align/session/failure/geometry/job regression coverage, and keyboard/390 px acceptance. Guessing missing history, reversing unrelated records, partial restoration and new direction/render/History owners are not accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-030: Reader-only legacy Similarity Group alignment compatibility

- Concern key: `diagram-generation.similarity-alignment.session-compatibility`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / LEGACY_READER_ONLY`
- Normative outcome:
  1. New Sessions persist the exact resolved alignment plan and do not write the
     legacy `alignOrthogroupFeature` group-ID string. Saving and loading a new
     Session preserves the exact reference, per-record resolved anchors, skips,
     effective transform intent, and reset baseline required by the other
     accepted decisions.
  2. Existing Sessions remain loadable through a bounded reader-only adapter
     that reproduces their historical implicit group-resolution behavior. The
     legacy resolver is unavailable to new alignment requests and to the normal
     writer path.
  3. A new alignment or supported edit converts the loaded state to the new
     resolved representation. The writer never downgrades a resolved plan to the
     legacy representation.
  4. Malformed or unsupported legacy values produce an explicit error without
     automatic substitution. The last successful Result remains visible when
     migration or validation fails.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-22` for issue `#561`, reproduced below.
  The receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.

```json
{
  "concern": "diagram-generation.similarity-alignment.session-compatibility",
  "scenarioRevision": 1,
  "choice": "A / LEGACY_READER_ONLY",
  "rationale": "Existing Sessions must remain loadable, but new Sessions must not perpetuate the ambiguous group-ID representation. Compatibility therefore belongs in a bounded reader-only adapter rather than the normal writer and runtime path.",
  "mustPreserve": "Reader-only reproduction of existing Sessions; the exact resolved plan in new Session round trips; explicit errors for malformed or unsupported legacy values; the last successful Result on migration or validation failure; and conversion to the new representation after a new Align or supported edit.",
  "mayRetire": "Writing the legacy alignOrthogroupFeature string in new Sessions; normal-runtime use of the legacy group resolver; and downgrade writing from a resolved plan to the ambiguous legacy representation.",
  "acceptedResidualRisk": "Replaying an old Session remains dependent on an isolated legacy resolver and can retain its historical implicit selection until the user creates a new alignment.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-22"
}
```

### PD-OI-031: Automatic resolved alignment with explicit display-direction review

- Concern key: `diagram-generation.similarity-alignment.surface-scope`
- Scenario revision: `5`
- Supersedes: `PD-OI-031`, scenario revision `4` (`A / AUTO_APPLY_RESOLVED_WITH_EXPLICIT_REVIEW_SINGLE_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.surface-scope",
  "scenarioRevision": 5,
  "choice": "A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW",
  "rationale": "Keep automatic resolved alignment uncomplicated while opened reviews expose mutually exclusive final display direction outcomes, including reversing only a minority reference through an all-right or all-left choice.",
  "mustPreserve": "Exact popup/drawer reference selection; automatic Python-resolved default alignment in Keep mode without forced review; an accessible ambiguity-required or explicit review with candidate facts, Select/Skip, resolution summary, one exclusive Keep/right/left/Custom direction selection and per-record resulting arrows; clear known-strand selected-anchor scope including the reference and unknown/skipped exclusion; shared typed validation, actionable underlying errors, strict CLI ambiguity rejection, existing CLI/API defaults, ordinary record controls, keyboard operation and the existing accepted narrow-screen palette coverage limitation. Apply persists absolute record transforms and placements, never review policies in the alignment plan.",
  "mayRetire": "An opened review exposing only a single reference-relative Match checkbox, permanent preservation of reference direction during an explicit direction operation, and mandatory post-Apply correction for per-record exceptions. Keep automatic/default and explicit review entry points.",
  "acceptedResidualRisk": "Explicit direction choices have more outcomes than the default path and can move the reference record's left edge while its selected feature center stays fixed. Use arrow-based labels, a single radio group, truthful target previews and Custom disclosure. This decision does not redesign narrow-screen palette coverage or introduce new CLI direction flags.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-032: Feature-popup rotation for one circular record

- Concern key: `diagram-generation.feature-popup-record-rotation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / POPUP-RECORD-ROTATION`
- Normative outcome:
  1. The open feature popup targets exactly one source-bound feature through its
     explicit stable record and biological-feature identity. It never falls
     back to a global selection. The operation is available for a complete,
     effectively circular record in either Circular or Linear diagram mode.
  2. Record actions expose the selected feature's 5-prime base, covered
     midpoint, and 3-prime base; a signed strand-relative offset; optional
     absolute forward orientation; and a distinct feature-end placement.
     Preview and resolution use original source coordinates, exact multipart
     traversal, and non-negative circular wrapping. Feature-end placement does
     not collapse into the 3-prime-base anchor.
  3. Apply updates only the target record's absolute display start and, when
     requested, absolute reverse-complement state. Leaving orientation off
     preserves its current value. Repeating the same operation is idempotent;
     every target-external record and layout value remains unchanged.
  4. Apply derives a target-only candidate from the last committed request, so
     unrelated pending form edits remain pending and are neither applied nor
     discarded. The existing sidebar workflow remains available and resolves
     the same display-transform meaning.
  5. Successful Apply admits the fresh Result, absolute transform, and
     non-authoritative provenance as one artifact-history transaction. One Undo
     or Redo restores or reapplies them together. Cancel and failed, stale, or
     superseded work leave the prior Result, transform, provenance, and History
     unchanged.
  6. New Sessions persist the absolute transform and provenance in Session 44,
     catalog 4. A bounded reader conservatively accepts released catalog 3;
     provenance never becomes rendering authority. Manual display-start or
     orientation changes clear stale anchor provenance without changing the
     effective transform.
  7. Invalid offsets and unsafe, ambiguous, fuzzy, cropped, linear, stale, or
     otherwise unsupported operations expose operation-specific reasons instead
     of truncating, guessing, or substituting another target. Duplicate record
     identifiers and split feature fragments retain stable source-bound
     identity.
  8. Feature, label, tick, depth, statistics, and comparison geometry follow the
     same record display transform. Source sequence, annotation, qualifiers,
     biological identity, and source-file export remain unchanged. Compatible
     LOSAT evidence is reused with zero additional executor jobs for a
     transform-only operation.
  9. Feature search and post-generation continuation remain available. Record
     actions are keyboard-operable in both rich and simple popup surfaces, show
     visible reason text, preserve the search query and stable target across
     Result replacement, and remain usable at a 390 px viewport.
  10. The request remains schema 7, and the existing Worker protocol and
      rendering path remain unchanged. The implementation adds no second
      request owner, Worker path, SVG admission path, History engine, or record
      rotation engine.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-22` for issue `#563`, reproduced below.
  The receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.

```json
{
  "concern": "diagram-generation.feature-popup-record-rotation",
  "scenarioRevision": 1,
  "choice": "A / POPUP-RECORD-ROTATION",
  "rationale": "Feature popupから対象featureを基準にrecordを直接回転できるようにし、sidebarとの往復や手動座標計算を減らす。source-coordinate preview、target-only適用、atomic Undo/Redoによって、操作結果を予測可能かつ安全にする。",
  "mustPreserve": "Source sequence、annotation、qualifiers、biological identity、対象外recordのtransformとlayout、未適用のform edits、既存sidebar workflow、canonical request owner、Worker経路、SVG sanitizer/admission経路、ResultとHistoryのowner、RecordDisplayTransform、LOSAT evidence reuse、searchおよびpost-generation workflow、failure/cancel/stale/superseded時の直前Resultとrecord transform。",
  "mayRetire": "none",
  "acceptedResidualRisk": "Popup UIおよびSession catalog compatibility pathの追加に伴う限定的なUI・保守負担を受容する。この負担は既存ownerの再利用、catalog 3からcatalog 4への単一のbounded reader、390 px・keyboard acceptance、AC-01～AC-20、およびfull regression gatesで制限する。科学的意味の変更、source dataの変更、global-selection fallback、追加LOSAT executor jobは受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-22"
}
```

### PD-OI-033: Feature-popup record-actions presentation

- Concern key: `web.feature-popup.record-actions-presentation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / EDIT_DISCLOSURE`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the exact text of all
  three `PRODUCT_DECISION` receipts presented for issue `#581` on
  `2026-09-24`. This serialization adds no terms to that approval and cannot
  authorize dependent runtime until merged into its base.

```json
{
  "concern": "web.feature-popup.record-actions-presentation",
  "scenarioRevision": 1,
  "choice": "A / EDIT_DISCLOSURE",
  "rationale": "通常のfeature確認・編集をすぐ始められる高さに保ちつつ、record回転をpopup内から見つけて使えるようにする。richとsimpleの両popupで同じ操作を提供する。",
  "mustPreserve": "開いたfeatureだけを対象とする回転、既存のanchor・offset・orientation・feature-end操作、適用前preview、操作できない理由の表示、keyboardと390 pxでの到達性、既存sidebar操作、成功時の一体的なResultとUndo/Redo、Cancel・失敗時の直前Resultとrecord transform。",
  "mayRetire": "featureを開くたびに回転フォームがタブより上へ自動展開する動作、およびフォーム内Cancelがfeature popup全体を閉じる動作。",
  "acceptedResidualRisk": "専用Recordタブより操作の分類は目立ちにくい。Editの上部に明確な見出しと開閉ボタンを置き、狭い画面とkeyboardで到達できることを確認する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-24"
}
```

### PD-OI-034: Local exclusive display-direction review with retry

- Concern key: `web.similarity-alignment.choice-and-retry`
- Scenario revision: `5`
- Supersedes: `PD-OI-034`, scenario revision `4` (`A / AUTO_APPLY_RESOLVED_WITH_REVIEW_RETRY_SINGLE_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "web.similarity-alignment.choice-and-retry",
  "scenarioRevision": 5,
  "choice": "A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY",
  "rationale": "One exclusive direction selection should govern local candidate and direction previews, and retry must retain user intent without combining a global policy with per-row overrides or committing an unseen reference reversal.",
  "mustPreserve": "Local no-Worker mode/custom/candidate/Skip editing; one direction tagged-union state with custom row values only in Custom; Python ownership of candidate eligibility and final source/display facts; one final batch validation per Apply attempt; one resolver for visible preview and final absolute orientations/reference-center placement; unchanged unknown/skipped/missing/unusable targets with reasons; an updated review and another Apply if final validated output differs; editable mode/custom choices after validation or render failure with the underlying message; prior Result, directions, placement and History after failed, canceled, stale or superseded work; canvas interaction, retained initially-Keep review after automatic render failure, atomic artifact Undo/Redo, Session/regeneration and the independently accepted Reset contract. Policies remain transient, and successful committed direction changes are actual record-state deltas, including the reference when changed.",
  "mayRetire": "The single draft-level Match flag, reference-relative direction as the only bulk operation, and post-Apply-only per-record exceptions. Do not add a second validation, rendering or History path or persisted direction policy.",
  "acceptedResidualRisk": "A final changed direction/placement receipt may require another Apply. Keep prior artifacts and local intent, show the new preview and limit re-review to actual output differences. Unknown directions are never guessed and stale reference/source bindings reject explicitly.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-035: Similarity-alignment canvas interaction

- Concern key: `web.similarity-alignment.canvas-interaction`
- Scenario revision: `3`
- Supersedes: `PD-OI-035`, scenario revision `2`
  (`A / RETAIN_MOBILE_PALETTE_COVERAGE`), only its permission for review
  coverage of the 390 px Preview, retirement of simultaneous canvas picking
  and pan/zoom, and continuation requiring review closure to inspect the
  diagram. Those permissions are no longer active.
- Status: `ACCEPTED`
- Normative outcome: the independent scenario-2 preservation requirements
  remain unchanged:

  正確な feature identity に基づく選択、候補一覧からの keyboard radio 選択と Skip、390 px での操作と適切な focus、描画されない候補の一覧からの選択、図上位置だけによる自動選択の禁止、デスクトップでの canvas 選択と手動 pan／zoom、レビューを閉じた後のプレビュー、ガイド・番号・draft が Result・download・Session に混入しないこと。

  The compact review presentation, visible canvas and simultaneous canvas
  interaction are governed by `PD-OI-039`, scenario revision `1`, in full.
  This record and `PD-OI-039` are jointly required; presentation does not
  replace identity, keyboard/Skip, non-rendered-candidate, desktop canvas,
  focus, or overlay-exclusion guarantees. No candidate may be selected from
  canvas position alone. The previous mobile coverage rationale and risk
  are retained in Git history, not as an active coverage exception.
- Decision source: the signed issue `#602` response for
  `web.similarity-alignment.review-presentation`, scenario revision `1`,
  reproduced in `PD-OI-039`, supplies exactly the rationale, preservation,
  limited retirement, accepted residual risk, owner, and date for this
  limited supersession. No separate human choice or rationale is inferred
  for the canvas-interaction concern. Other transform, plan, Reset, History,
  and alignment outcomes, including `PD-OI-031`/`PD-OI-034`, retain their scope.
  Dependent runtime requires this supersession merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-026`.

### PD-OI-036: Linear record-label Auto visibility and disclosure

- Concern key: `linear.record-label-auto-visibility`
- Scenario revision: `2`
- Status: `ACCEPTED`
- Selected outcome: `B / AUTO-FRESH-RESET-WITH-DISCLOSURE`
- Scenario revision `2` is the signed scenario; no prior record for this
  concern exists in the base Contract. Fresh/reset remains independently
  Auto for Accession and Length; default Show (01-A) is not selected.
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-022`.

```json
{
  "concern": "linear.record-label-auto-visibility",
  "scenarioRevision": 2,
  "choice": "B / AUTO-FRESH-RESET-WITH-DISCLOSURE",
  "rationale": "共有行の図では簡潔な既定表示を維持し、情報が非表示になる理由とShowへの変更先を配置操作の場所で明示する。",
  "mustPreserve": "fresh/resetの独立Auto、Show/Hideの明示値、diagram-wide Auto解決、休眠行除外、既存Sessionと保存Result、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。",
  "mayRetire": "なし。",
  "acceptedResidualRisk": "共有行でAccession/Lengthが非表示になる結果自体は残る。説明を見落とす可能性があるため、layout操作場所とLabelsの両方で実効値と変更先を示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-037: Derived edit-application feedback

- Concern key: `web.edit-application-feedback`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / DERIVED-APPLICATION-STATUS`
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-005`, `OIC-013`, `OIC-014`, `OIC-024`.

```json
{
  "concern": "web.edit-application-feedback",
  "scenarioRevision": 1,
  "choice": "A / DERIVED-APPLICATION-STATUS",
  "rationale": "利用者が編集後の図と次回Generateの変更を区別できるよう、操作の適用タイミングと生成設定の未適用状態を事実から表示する。",
  "mustPreserve": "操作単位のLive edit、Applies on Generate、Apply requiredを区別し、左Palette Instant Previewと右Alignment reviewを例外なく正しく分類する。canonical即時commitと必要時自動rerender、reviewのlocal draft、target-only操作とpending設定の分離、対応済みoverride継承、atomic Generate、失敗/Cancel/stale時の旧ResultとHistory、Undo/Redo、SessionのResult/draft分離、Exportの現在Result出力を維持する。Pendingとlive applying/errorを独立に示し、invalid/unknownをAppliedとしない。Generateによる配置再計算とzoom reset、Save/Exportの意味を事前に説明する。",
  "mayRetire": "製品の適用タイミング・保存・復旧・編集機能は退役しない。全DOM座標や手動位置を再生成後にも無条件に保持する保証は新設しない。",
  "acceptedResidualRisk": "生成intent比較の漏れや比較基準の誤更新は誤表示を生みうる。生成・live commit・履歴・Sessionの一致テストを必須とし、根拠不足はunknownとして表示する。Status目的のWorker呼出、genome byte読取り/hash、SVG/checkpoint cloneは受け入れない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-038: Compact Editor presentation

- Concern key: `web.editor.compact-presentation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / DOCKED-COMPACT-EDITOR`
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-025`.

```json
{
  "concern": "web.editor.compact-presentation",
  "scenarioRevision": 1,
  "choice": "A / DOCKED-COMPACT-EDITOR",
  "rationale": "狭いPreviewでも即時編集の変化を図で確認できるよう、図とEditorを上下の領域へ配置する。",
  "mustPreserve": "同じSVGとEditor、全tabと同期可用性、canonical live commitと必要時rerender、既存History/Session/Export、camera操作、keyboard、Close/Escapeのvisibility-only意味、選択tab、Result置換/失敗復旧。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、Editor内容を独立scrollさせ、Close/headerとtoolbarを操作可能にする。短いviewport/soft keyboardでは全操作へscrollで到達できる。wideのside drawerを維持する。",
  "mayRetire": "狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。編集機能や保存意味は退役しない。",
  "acceptedResidualRisk": "上下分割で図とEditor listの縦領域が短くなり、list scrollが増える。実操作のpointer/keyboard/browser検証を必須とし、複製Preview・SVG clone・第二editorによる回避は受け入れない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-039: Compact Similarity-alignment review presentation

- Concern key: `web.similarity-alignment.review-presentation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / DOCKED-COMPACT-ALIGNMENT-REVIEW`
- Supersession scope: replaces only the scenario-2 mobile coverage
  exception and close-review continuation in `PD-OI-035`; its independent
  guarantees remain jointly required in scenario revision `3`. Narrow free
  drag and concurrent Editor opening during review may retire only as stated
  in this receipt. Wide drag and non-modal canvas interaction remain required.
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-026`.

```json
{
  "concern": "web.similarity-alignment.review-presentation",
  "scenarioRevision": 1,
  "choice": "A / DOCKED-COMPACT-ALIGNMENT-REVIEW",
  "rationale": "狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。",
  "mustPreserve": "PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、1つのMatch reference directionと各targetの結果方向、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。",
  "mayRetire": "狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続だけ。候補や方向の選択、Apply前draft、failure/retryは退役しない。",
  "acceptedResidualRisk": "狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## Acceptance contract catalog

| Contract | Required meaning |
| --- | --- |
| `OIC-001` | Candidate limit is truthful; `None` remains unbounded and no hidden cap is applied. |
| `OIC-002` | Candidate, Pairwise display, and member-hit limits are independent and invalidate only the correct stages. |
| `OIC-003` | Every supported Collinear enum reaches the real typed Python analysis path. |
| `OIC-004` | Fresh defaults follow the declared surface and mode rules; explicit imported values win. |
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
| `OIC-015` | Linear discovery, per-record placement, actual comparison jobs, explicit request endpoints, Session replay, and SVG endpoints retain the complete selected record universe. Adjacent uses neighboring-row Cartesian products; Similarity and Collinear `all` retain complete directed between-record evidence regardless of display rows; Collinear within-record evidence follows the inference checkbox. |
| `OIC-016` | LOSAT Auto allocation follows Total threads; displayed and executed budgets agree, and manual intent survives temporary clamps and Session replay. |
| `OIC-017` | Web raw/member defaults are 5/5 in Collinear and unbounded/unbounded in Similarity. Each mode restores its own edits repeatedly, including blanks; Session round trips retain both modes; Reset Settings restores defaults. |
| `OIC-018` | Collinear inference defaults OFF; actual raw jobs exclude every self-comparison, including within multi-record source batches, and the real Python path skips orthogroup inference. ON retains the existing inference; request, cache, provenance, and legacy Session interpretation agree. |
| `OIC-019` | Completed raw searches survive downstream cancellation for matching retries; member-only edits do not rerun LOSAT. Raw-setting/input changes, Clear Cache, and Session/History replacement prevent incompatible reuse; the committed Result remains intact. |
| `OIC-020` | Linear File cards expose common Depth TSV assignment without expanding records. File-level apply and clear update only that File and logical series as one undoable operation; empty, common, and mixed states remain truthful. Per-record sparse overrides, logical indexes, canonical requests, Session replay, and regeneration remain unchanged. |
| `OIC-021` | Feature-popup record rotation uses the explicit popup target and source coordinates, changes only one effectively circular record through a target-only atomic transaction, preserves pending edits and the prior artifact on every no-op path, round-trips the absolute transform, and reuses compatible LOSAT evidence without additional executor jobs. |
| `OIC-022` | `PD-OI-036`: independent fresh/reset Auto, explicit Show/Hide, diagram-wide resolution excluding dormant rows, saved Session/Result and Undo/Redo remain supported; Layout and Labels disclose the effective Auto result, affected fields, reason, scope, next Generate effect, and route to Record Labels. |
| `OIC-023` | `PD-OI-024`: only Web fresh/reset selects Lock ON; explicit OFF, saved values and supported omission meanings, saved Result on Load, CLI/Python defaults, accepted anchors, D2-P and D3-A remain supported. ON/OFF and Generate application are always explained in Linear Layout. |
| `OIC-024` | `PD-OI-037`: operation-level Live edit, Applies on Generate, and Apply required are truthful, including Palette Instant Preview and Alignment review. Pending and live applying/error are independent; invalid/unknown is never Applied. Generate/live commit/History/Session observations agree, with no status-only Worker, genome-byte read/hash, or SVG/checkpoint clone. |
| `OIC-025` | `PD-OI-038`: one SVG and Editor retain live commit/rerender, all tabs, availability, History/Session/Export, camera, keyboard, visibility-only Close/Escape, selected tab and Result recovery. At 390×844/740 the canvas uses the available full width and at least 200 px height; content scroll, reachable header/Close/toolbar, short-viewport/soft-keyboard access, and wide side drawer remain required. Pointer/keyboard/browser verification is required; duplicated Preview/SVG/editor is not accepted. |
| `OIC-026` | `PD-OI-035` and `PD-OI-039`: identity, keyboard Select/Skip, non-rendered candidates, no position-only selection, desktop canvas, focus and transient overlay exclusion remain required with all PD-OI-031/034 outcomes. Compact review retains visible, operable canvas at full available width and at least 200 px height at 390×844/740, scrollable candidates and reachable Apply/Cancel, local no-Worker draft edits, atomic batch validation, failure/error/retry and artifact/orientation/History recovery. Narrow review closes Editor through its owner while retaining tab, disables reopening with a reason until review ends, then permits explicit reopen; wide drag remains. Browser verification must show presentation changes leave draft and Result unchanged. |

These new acceptance entries are obligations for dependent runtime work, not
claims of completed runtime or browser verification by this authority amendment.
The signed receipts remain the complete outcome; existing acceptance contracts
and independent preserved guarantees remain jointly required.

### OIC-020 required regression coverage

The normal automated PR gate must observe all of the following:

- A multi-record Linear GenBank File exposes its Depth TSV assignment while its
  record list and record options remain closed. Applying one file binds the same
  logical series to every record in that File and does not affect another File.
- One per-record replacement produces a truthful mixed File state. Applying or
  clearing the File-level value then replaces or clears every record cell in
  that File and series, and one Undo restores the complete prior matrix.
- Empty cells and later logical columns do not shift when a source is cleared.
  Same-named independent files remain distinct, including after Session
  restoration.
- Save, fresh Load, canonical request construction, generation, and subsequent
  regeneration preserve common and mixed bindings without a new Session schema,
  request schema, Worker protocol, or rendering path.

### OIC-021 required regression coverage

| ID | Required observation |
| --- | --- |
| `AC-01` | The popup alone rotates the target feature's record and never consults another or global selection. |
| `AC-02` | An effectively circular record resolves the same source request in Circular and Linear diagram modes. |
| `AC-03` | Every non-target record and layout value remains unchanged, including same-file multi-record inputs. |
| `AC-04` | Positive and negative offsets resolve relative to feature direction and wrap correctly. |
| `AC-05` | Orientation intent is absolute and idempotent; leaving it off preserves the current value. |
| `AC-06` | The 3-prime base anchor and feature-end placement remain distinct operations. |
| `AC-07` | Multipart, origin-spanning, and odd/even covered midpoints follow exact covered traversal. |
| `AC-08` | Unstranded, ambiguous, fuzzy, cropped, linear-topology, and stale cases expose operation-specific reasons. |
| `AC-09` | Duplicate record IDs and split fragments retain stable source-bound identity. |
| `AC-10` | One Undo or Redo restores or reapplies origin, orientation, provenance, and Result together. |
| `AC-11` | Save and fresh Load restore the absolute transform and provenance. |
| `AC-12` | Failed, canceled, stale, and superseded rendering preserves the previous Result and transform. |
| `AC-13` | Unrelated pending edits remain pending and are neither applied nor discarded. |
| `AC-14` | A transform-only operation adds zero LOSAT executor jobs and reuses compatible raw evidence. |
| `AC-15` | Feature, label, tick, depth, statistics, and comparison geometry use the same transform. |
| `AC-16` | Manual display-start or orientation changes clear stale anchor provenance. |
| `AC-17` | Cancel changes no draft, Result, transform, provenance, or History state. |
| `AC-18` | Search query and stable target re-identification survive replacement; popup actions remain keyboard- and 390 px-accessible. |
| `AC-19` | Request schema 7, the Worker protocol, and the renderer path do not expand. |
| `AC-20` | Product Impact and Architecture Ratchet evidence remain reviewable and all required gates pass. |

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
- Similarity and Collinear `all` with inference ON retain the complete directed
  record-pair matrix through source-file jobs, including self, same-row, and
  non-adjacent evidence. Two sources containing five records require four
  source jobs and cover 25 directed record pairs. Collinear `all` with inference
  OFF retains all 20 between-record directions and searches no record against
  itself, including within multi-record sources. A single-row layout retains
  the selected analysis even though it has no between-row links.
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


## OIC-016 acceptance evidence

- With four source jobs, change Total threads between 32, 16, and 2 through the
  visible control and observe both Auto labels recalculate as specified above.
- Exercise one explicit concurrency value with automatic per-run threads and
  one explicit per-run thread value with automatic concurrency. Verify the
  effective budget bound after decreasing and increasing Total threads.
- Preserve explicit choices and Auto through Save and fresh Load, including a
  temporarily clamped value. Show the effective clamp without rewriting intent.
- Observe real threaded LOSATP dispatch and completion on a small workload;
  compare the worker allocation and runtime report with the displayed plan.
- Keep fixed one-thread-per-run programs fixed. Changing only scheduling
  settings preserves the existing raw-search cache identity.

## OIC-017–OIC-019 regression evidence owners

These tests protect the outcomes; they do not supply Product authority.
The mode-restoration, inference-toggle, and cancellation/retry browser
observations below must run in the normal automated PR gate. Full dev staging
provides additional coverage and does not replace this pre-merge requirement.

| Contract | Executable regression owners |
| --- | --- |
| `OIC-017` | [`comparison-ui.playwright.spec.js`](../../tests/web/comparison-ui.playwright.spec.js), `comparison controls drive appearance and current Session round trips`: initial values, repeated mode restoration, saved inactive values, blanks, and Reset Settings. |
| `OIC-018` | [`linear-multi-record.playwright.spec.js`](../../tests/web/linear-multi-record.playwright.spec.js), `Collinear inference checkbox skips self searches and reuses matching evidence`; [`linear-sources.test.mjs`](../../tests/web/linear-sources.test.mjs); [`losat-settings.test.mjs`](../../tests/web/losat-settings.test.mjs); [`test_collinearity.py`](../../tests/test_collinearity.py); [`test_session_request_codec.py`](../../tests/test_session_request_codec.py); [`test_api_request_render.py`](../../tests/test_api_request_render.py). Protect actual job endpoints, no inference call when OFF, member selection, both scopes, all anchor modes, and request/derived identities. |
| `OIC-019` | [`linear-multi-record.playwright.spec.js`](../../tests/web/linear-multi-record.playwright.spec.js), `protein raw cache survives cancellation and derived options preserve search identity`; [`run-analysis-simple-path.test.mjs`](../../tests/web/run-analysis-simple-path.test.mjs). Protect cancellation, repeated retry, member edits, raw-setting/input changes, clearing, and late artifact failure. |
