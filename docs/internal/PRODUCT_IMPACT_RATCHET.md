# gbdraw Product Impact Ratchet

Status: normative repository policy

This policy defines how a registered architecture change is traced to supported
product behavior and how a developer handles material user-effect ambiguity
that is not yet mapped to an architecture subject. It governs product-impact
classification, procedural preflight, authority resolution, Decision Pack
content, and product-decision eligibility. The repository now implements the
version 1 map, decision authority, evaluator, bounded decision source, trusted
checker integration, and two-rule hard pilot described here. Unmapped
developer preflight remains a procedural intake into this same model; it is not
a second evaluator or machine-readable decision store.

## Scope and delegated authority

The Product Impact Ratchet joins existing policies without replacing them:

- [`ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`](./ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
  remains authoritative for semantic-owner and canonical-path definitions,
  `OE`, `PE`, `CB`, architecture convergence, exceptions, and architecture Gate
  behavior.
- [`WEB_CHANGE_POLICY.md`](./WEB_CHANGE_POLICY.md) remains authoritative for
  the delivery lifecycle, trusted admission, independent Gate and Review
  results, workflow ownership, and self-authorization separation.
- [`gbdraw/web/CLAUDE.md`](../../gbdraw/web/CLAUDE.md) records the current Web
  runtime boundaries and explicit supported behavior.
- This policy owns architecture-to-behavior traceability, product-impact
  classification, developer preflight for unmapped material user-effect
  ambiguity, Product Impact Packet and Decision Pack semantics, and
  product-decision eligibility.

A Product Impact decision cannot waive an architecture, security, performance,
scientific-output, or behavior-test failure.

## Objective

Preserve or improve owner/path de-duplication while making the affected user
journey, supported behavior, authority, decision, evidence, and residual risk
explicit. Keep the common path automatic and stop only disputed product
outcomes.

Maintainers select canonical behavior. The canonical semantic owner and
production path are implementation consequences. Files are not presented as
product options unless they embody materially different user-visible outcomes.

Equivalent owner or path substitution is an architecture decision when it
preserves every mapped requirement, effect, checkpoint, and contract. A change
between materially different supported outcomes is a product decision.

## Product model

### Architecture subject

An architecture subject is the canonical identity emitted by a versioned
detector for one source fact. Examples include a definition path and a
canonical entry edge. Subject strings come from the detector's stable encoder.

### Product concern

A product concern is the unit that relates architecture subjects to one
user-facing capability or continuity requirement. It has a stable key and a
scenario revision.

### User journey and checkpoint

A user journey records an actor, context, goal, and major steps. A journey
checkpoint is a specific point where the evaluator compares behavior, such as
request submission, a fresh Result, Session round trip, regeneration, export,
or recovery after failure.

### User effect

A user effect is a stable, user-observable result, available continuation,
preserved state, or contract at a checkpoint. A function or file name alone is
not a user effect.

### Behavior option and realization requirement

A behavior option is one product-level outcome. A realization requirement is
an independent contribution needed for that outcome. Requirements may cover
semantic behavior, an affordance, or compatibility.

Realization uses AND-of-OR semantics:

```text
option realized = all requirements satisfied
requirement satisfied = any approved subject active
```

Subjects within one requirement are declared substitutes. Separate
requirements are jointly necessary. A shared option ID does not make separate
requirements interchangeable.

### Product authority

Product authority is a recognized source that selects the supported behavior
option. Static supported behavior, a domain or scientific rule, an active
durable decision, or an eligible current maintainer decision may supply it.

### Product Decision Owner

A Product Decision Owner is a repository maintainer named in the base-branch
maintainer allowlist. This role is distinct from an ordinary gbdraw application
user. The Product Decision Owner selects the product outcome, product-level
rationale, intentional retirement scope, and accepted residual risk.

### Durable behavior decision

A durable behavior decision is an active base-branch `BD-###` record that may
govern later pull requests for the same concern and scenario revision.

### Current maintainer decision

A current maintainer decision is an exact-head declaration in one pull request.
It can authorize only a mapped, behavior-preserving transition for that pull
request and is not authority for later changes.

### Human product-decision response

A human product-decision response is the Product Decision Owner's short choice,
rationale, preservation and retirement intent, and risk statement. It is not CI
input.

### Behavior contract

A behavior contract is a referenced test or acceptance check that protects a
mapped requirement or checkpoint. Contracts are evidence, not product
authority.

### Product Impact Packet and Decision Pack

A Product Impact Packet is the deterministic report for one changed concern. A
Decision Pack is the additional, outcome-oriented report produced when human
product judgment is required.

The non-normative working template is
[`PRODUCT_DECISION_PACKET_TEMPLATE.md`](./PRODUCT_DECISION_PACKET_TEMPLATE.md).
It applies to both registered Product Impact changes and developer preflight
for unmapped ambiguity. This policy, not the template, defines the lifecycle.

### Procedural preflight classification

A procedural preflight classification selects the developer's next action when
the proposed change exposes a product question:

```text
IMPLEMENT_EXISTING_AUTHORITY
EVIDENCE_REQUIRED
PRODUCT_DECISION_REQUIRED
NOT_ALLOWED
```

These classifications do not replace or combine the impact class, authority
resolution, evaluation observation, or per-choice Decision Pack route.

### Complete coverage

Complete coverage means every allowed subject for the registered architecture
rule appears in at least one requirement candidate, every observed changed
subject is mapped, and all affected checkpoints and blocking requirements have
the required contracts. Partial coverage may report known mappings but cannot
claim automatic safety.

## Actors and the human decision interface

The roles are separate:

- A gbdraw end user supplies product needs and usability feedback but does not
  edit repository authority through ordinary application use.
- The Product Decision Owner selects a product outcome and states the rationale,
  preservation or retirement scope, and accepted residual risk.
- A developer or Codex analyzes alternatives, implements and tests the chosen
  behavior, and serializes an explicit human decision. Codex must not select,
  infer, or broaden the outcome.
- The trusted checker validates the serialized representation. It never ranks
  options or infers product preference.
- A reviewer confirms that the machine representation matches the human choice
  and that the cited evidence covers the proposed head.

One maintainer may act as both developer and Product Decision Owner, but the
roles and sequence remain separate: prepare the Decision Pack as developer,
issue the explicit human receipt as Product Decision Owner, serialize only
that receipt, and implement only after required durable authority has merged.

Version 1 uses a GitHub Actions Summary and a Codex-mediated response. It adds
no dedicated decision UI, external service, general end-user vote, or comment
bot. Humans do not enter `headSha`, architecture subjects, requirement
references, evidence references, or JSON.

The Decision Pack supplies a ready-to-copy response with resolved concern and
revision fields:

```text
PRODUCT_DECISION
Concern: <concern key/title>
Scenario revision: <revision from Decision Pack>
Choice: <A-F and/or stable option ID>
Rationale: <product-level reason>
Must preserve: <effects/affordances>
May retire: <none or explicit scope>
Accepted residual risk: <bounded risk or none>
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

CI does not parse this response. After the human chooses an outcome, Codex may
convert it to bounded machine JSON in the pull request body or prepare an
authority-only pull request. The trusted checker parses only that machine
representation. Missing rationale, retirement intent, or accepted risk remains
unresolved; Codex does not fill it by inference. Owner and decision date are
also explicit for a durable authority receipt.

The initial bounded authority set may use one receipt only when it identifies
the exact candidate content and every included, excluded, or modified record:

```text
PRODUCT_DECISION_SET
Candidate file: <candidate file>
Candidate SHA-256: <exact digest>
Approved decision IDs: <complete explicit list>
Records remaining EVIDENCE_REQUIRED: <complete explicit list>
Excluded or modified records: <none or complete wording>
Rationale / Must preserve / May retire / Accepted residual risk:
  Approved exactly as written for the listed records, except for the explicit modifications above.
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

This bounded receipt is a human review aid for the initial authority-only
bootstrap, not a new machine-readable store. A missing digest, complete ID
accounting, owner, or date leaves the authority unresolved.

## Product Decision intake

The lifecycle has two trigger lanes:

```text
Lane A: a registered Product Impact architecture subject changes
  -> deterministic Product Impact evaluation

Lane B: developer preflight finds unmapped material user-effect ambiguity
  -> procedural intake into the same Product Decision model
```

Both lanes use one authority search, procedural classification, Decision Pack,
Product Decision Owner response, storage rule, supersession rule, and
trusted-base sequence. Lane B does not add an evaluator state machine, JSON
registry, or automatic Gate result.

### Mandatory developer preflight

Perform preflight before implementation when a proposed change may alter any of
the following and merged authority does not select one complete outcome:

- a fresh public or persisted default;
- accepted names, types, enum domains, ranges, aliases, or omission semantics;
- GUI discoverability, editability, read-only preservation, or unsupported
  disposition;
- Session reconstruction, replacement, clearing, migration, or replay;
- Result replacement, failed-generation recovery, cancellation, supersession,
  or stale completion;
- compatibility retention or retirement;
- scientific-output selection or interpretation;
- cache reuse that can change output correctness;
- a user actor, goal, entry point, checkpoint, failure path, or next available
  action;
- which product behavior survives semantic-owner or canonical-path
  convergence; or
- a current implementation accident that might otherwise become a durable
  promise.

The developer decides whether Product judgment is required but does not
silently select the outcome.

### Cases that normally need no Product Decision

Use `IMPLEMENT_EXISTING_AUTHORITY` without a new Product Decision when:

- runtime or a test plainly violates unambiguous merged authority;
- a private refactor preserves every product effect, checkpoint, public
  contract, compatibility commitment, and next action;
- a duplicate owner or path is removed while the same authorized outcome
  remains;
- a named deterministic integrity, safety, specification, security, or
  scientific rule permits only one correct result; or
- documentation or evidence is corrected to match authority without changing
  the product outcome.

An incorrect test does not turn a bug into a product choice. Correct the test
and implementation against merged authority.

### Required procedural classification

Select exactly one classification and explain why the others do not apply:

- `IMPLEMENT_EXISTING_AUTHORITY`: one merged authority or non-waivable rule
  selects the complete outcome.
- `EVIDENCE_REQUIRED`: deterministic evidence is needed before choosing among
  product-valid outcomes.
- `PRODUCT_DECISION_REQUIRED`: two or more product-valid outcomes remain after
  available evidence.
- `NOT_ALLOWED`: the candidate violates an architecture, security, integrity,
  scientific-output, persisted-compatibility, or required-evidence constraint
  that a Product Decision cannot waive.

`PRODUCT_DECISION_REQUIRED` is mandatory when the authority search and required
evidence leave two or more materially different product-valid outcomes. It is
also mandatory before intentionally selecting a new outcome or retiring an
existing default, affordance, continuation, compatibility promise, or
scientific-output behavior that merged authority does not already decide.

When the classification is `EVIDENCE_REQUIRED` or
`PRODUCT_DECISION_REQUIRED`, dependent runtime work stops.

## Independent state axes

Product Impact results keep three questions separate.

Impact class describes the user-visible change:

```text
NO_USER_VISIBLE_DIFFERENCE
AFFORDANCE_PRESERVED
PRODUCT_CHANGE
RETIREMENT
```

Authority resolution identifies the recognized decision source:

```text
STATIC_AUTHORITY
DOMAIN_AUTHORITY
DURABLE_DECISION
CURRENT_MAINTAINER_DECISION
UNRESOLVED
CONFLICT
```

Evaluation observation describes the result:

```text
NOT_APPLICABLE
CONFORMING
ORDINARY_REGRESSION
UNRESOLVED_DECISION
AUTHORITY_CONFLICT
INSUFFICIENT_EVIDENCE
```

Older disposition names must not be kept as a second model that combines these
axes.

## Strict no-difference rule

`NO_USER_VISIBLE_DIFFERENCE` is valid only when all of these conditions hold:

- the before and after realized option sets are identical and nonempty;
- every requirement has the same satisfaction value before and after;
- no affordance or compatibility requirement is lost;
- the stable user-effect sets are identical;
- no checkpoint loses effect or contract coverage;
- mapped contracts are adequate for the affected checkpoints and requirements;
- complete coverage contains no unmapped changed subject.

Option-ID equality alone is insufficient. A provider may be replaced by
another approved subject in the same requirement, but the full requirement,
effect, checkpoint, and contract comparison still applies.

## Authority resolution

Recognized authority has four forms:

- Static authority is an explicit supported behavior or compatibility promise
  in a base-branch normative or runtime-contract document.
- Domain authority is a named scientific, specification, or correctness rule
  that permits only one valid outcome.
- A durable decision is an active base-branch `BD-###` decision for the exact
  concern and scenario revision.
- A current maintainer decision is an eligible exact-head declaration for one
  behavior-preserving pull request.

Current code establishes active source facts. Tests establish what is checked.
Neither is automatic product authority. If recognized authorities select
incompatible options, the result is `CONFLICT` and `AUTHORITY_CONFLICT`; the
checker does not choose one by precedence.

Search all relevant base-branch sources rather than treating this list as a
precedence order: an applicable Product Impact concern and its static
authority, an active `BD-###`, the static Product Contract when it exists, a
named domain or non-waivable rule, released compatibility evidence, and an
eligible exact-head decision for its existing narrow route. Current code,
tests, fixtures, screenshots, reports, and historical behavior remain evidence
only. At most one active authority selects one concern and scenario revision;
incompatible active authorities are a conflict, not a tie to resolve by file
order, timestamp, current implementation, or convenience.

## Durable decisions

The durable authority path is:

```text
tools/web-product-decisions.json
```

It contains active decisions and the base maintainer allowlist. At most one
active decision may exist for a concern and scenario revision. Replacement
receives a new `BD-###`; the former record is removed from the active file and
remains in Git history.

A material product change, public-contract change, affordance or compatibility
retirement requires durable base authority before runtime implementation.
Codex cannot author an accepted product choice without an explicit Product
Decision Owner disposition.

### Future static Product Contract authority

The exact preauthorized future static Product authority path for broad or
currently unmapped durable Option Integrity outcomes is:

```text
docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
```

This path is bounded to durable Product outcomes that are not already
faithfully owned by the Product Impact map and `BD-###` store. It is not a
parallel evaluator or machine-readable decision store. Mapped concerns continue
to use the existing authority model, and documentation or tests may restate or
protect an outcome but may not become a competing authority.

The path declaration selects no Option Integrity outcome. The file is absent
at this policy revision and therefore supplies no authority. The declaration
also provides no executable protection until a later checker-only pull request
adds exact-path enforcement and that checker change has merged.

## Current decisions

The optional pull-request-local mechanism is narrower than durable authority.
A nonempty current decision is eligible only when all of these conditions hold:

- the pull request author is in the base maintainer allowlist;
- the declaration binds to the exact event head SHA;
- the base map marks the concern as `decision-required`;
- `acceptedImpactClass` is `AFFORDANCE_PRESERVED`;
- the selected option is fully realized after the change;
- before and after stable effects match;
- no affordance or compatibility requirement is retired;
- changed requirement references exactly match evaluator output;
- `rationale` is nonempty and product-level;
- evidence references and residual risk are nonempty, but neither replaces a
  mapped contract.

The Product Decision Owner supplies the outcome, rationale, preservation intent,
and risk. Codex derives the head SHA and exact architecture, requirement, and
evidence references. The generated machine block must be shown for review.

A current decision supplies no future authority and cannot authorize
`PRODUCT_CHANGE` or `RETIREMENT`. A synchronize event changes the head SHA and
makes the declaration stale.

## Contracts and evidence integrity

Contract sensitivity is one of:

```text
BASE_FAIL
MUTATION_FAIL
DIRECT_ASSERTION
MANUAL_ACCEPTANCE
```

Contract execution is one of:

```text
PR_GATE
DEV_STAGING
MANUAL
```

Hard coverage requires automated `PR_GATE` evidence for every blocking
requirement and affected checkpoint. Manual evidence may supplement that proof
but cannot by itself establish hard automatic safety.

An affected runtime pull request cannot modify, delete, rename, or replace its
mapped contract and use that change as the sole proof of hard safety. When a
mapped reference must change, first merge an evidence-only pull request, then
update the authority in a separate authority-only pull request. Unrelated test
changes are unaffected.

## Evidence before decision

Use `EVIDENCE_REQUIRED` when deterministic evidence can resolve or materially
bound the choice before Product judgment. An evidence-only pull request:

- states the exact question and compared explicit outcomes;
- uses deterministic representative inputs;
- records method, commands, raw observations, interpretation limits,
  performance where material, and reproducibility data;
- does not change the disputed default, Product authority, runtime selection,
  compatibility promise, or expected-output baseline;
- may state a clearly non-authoritative engineering recommendation; and
- does not exclude a product-valid alternative merely because it is harder to
  implement.

Evidence does not select a Product outcome. When evidence still leaves two or
more product-valid outcomes, prepare the Decision Pack and obtain an explicit
Product Decision Owner receipt.

## Authority correction and supersession

Correct implementation or evidence against existing authority without creating
a new Product Decision merely because current code or a test is wrong. To
replace a wrong durable Product outcome, prepare evidence or a Decision Pack,
obtain a new explicit receipt, and merge an authority-only supersession before
dependent runtime changes.

Exactly one authority remains active for the concern. A Product Contract
replacement increments its scenario revision, identifies the prior decision
and revision in `Supersedes`, and replaces the complete active outcome; Git
history preserves the former text. A mapped durable replacement creates a new
active `BD-###`, removes the prior active record, and relies on Git history for
the former record. Resolve authority conflicts in an authority-only pull
request rather than establishing precedence through ordering or timestamps.

## Trusted-base admission

Runtime admission uses only:

```text
base checker
base detectors
base architecture rules
base product-impact map
base durable decisions and maintainer allowlist
base source facts
head source facts
bounded inert PR-body decision
```

Candidate maps, rules, durable decisions, checker code, workflows, or
replacement evidence may be validated as inert data. They cannot authorize the
same candidate runtime. Trusted workflows must not import or execute candidate
checker or detector code.

## Rule-level rollout

Each mapped architecture rule declares:

```text
coverage: partial | complete
enforcement: report-only | hard
```

`hard` requires `complete`. The initial pilot covers only:

```text
canonical-path.render-request
semantic-owner.render-request
```

Both pilot rules have complete `hard` enforcement after their report-only audit
and evidence sequence. Every later concern starts `report-only`; its evidence
is reviewed before the mapped rules move to hard enforcement.

## Gate and Review

Malformed authority, schema violations, unsafe parsing, and broken
self-authorization separation are policy failures. Report-only runtime
observations do not change the executable exit status.

Under hard enforcement, `ORDINARY_REGRESSION`, `UNRESOLVED_DECISION`,
`AUTHORITY_CONFLICT`, and `INSUFFICIENT_EVIDENCE` fail the Gate. `CONFORMING`
passes. Review remains independent and is required for architecture-bearing or
product-decision work.

## Low-friction applicability

Product Impact evaluation runs only for:

- a registered architecture subject delta;
- a nonempty current-decision block that requires validation;
- an authority-only Product Impact proposal that requires schema validation.

A production deletion, new module, label, churn threshold, or
compatibility-like name does not independently create a Product Impact decision
requirement. Ordinary pull requests with no registered delta leave the decision
array empty and perform no Product Impact form work.

Developer preflight for an unmapped material user effect is the procedural Lane
B defined above. It does not expand executable Product Impact applicability,
populate the current-decision array, or create a new automatic Gate result.

## Promotion of unmapped concerns into automation

Promote a procedurally handled concern into deterministic Product Impact
mapping only when it recurs or is predictably high risk, a stable architecture
subject exists, and a deterministic detector can identify that subject
faithfully. Actor, goal, checkpoints, effects, requirements, and behavior
contracts must be complete; false-positive and false-negative behavior must be
bounded; and report-only evaluation must first demonstrate useful signal.

The map references the active Product authority rather than duplicating it. A
manual durable contract remains preferable to an inaccurate detector or a
one-off automatic rule.

## Product Impact Packet and Decision Pack

The Product Impact Packet reports the concern, changed rules and subjects,
affected journeys and checkpoints, before and after option realization,
requirement changes, stable effects, authority, contracts, residual risk,
impact class, authority resolution, and evaluation observation. A concern is
reported once even when multiple rules changed.

When the observation is `UNRESOLVED_DECISION`, the Decision Pack also includes:

- the actor, context, goal, action or completion, and affected checkpoints;
- the before and after options and a requirement realization matrix;
- preserved, added, and lost stable effects;
- the authority search and its result;
- contract evidence, gaps, and residual risk;
- removal consequences and the allowed next action for every choice;
- a prefilled `PRODUCT_DECISION` response template.

Each selectable outcome has a stable choice code, a stable option ID when one
exists, its preserved, added, and lost effects, and one route:

```text
PR_LOCAL_ALLOWED
DURABLE_AUTHORITY_REQUIRED
EVIDENCE_REQUIRED
NOT_ALLOWED
```

`PR_LOCAL_ALLOWED` is limited to an existing mapped outcome with equivalent
effects. Material change, a new combined or distinct workflow, or retirement
uses `DURABLE_AUTHORITY_REQUIRED`. Missing decision evidence or contracts uses
`EVIDENCE_REQUIRED`. An outcome that violates architecture, security, domain,
or other non-waivable requirements is `NOT_ALLOWED`.

Choices describe product outcomes rather than filenames. Comparison dimensions
include entry point and discoverability, immediate feedback, canonical state
update, Undo/Redo, persistence, regeneration, export, validation and error,
recovery, performance, accessibility, compatibility, and the next available
action when relevant.
