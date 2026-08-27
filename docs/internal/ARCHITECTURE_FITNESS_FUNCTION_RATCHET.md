# gbdraw architecture fitness-function ratchet

Status: normative repository policy

This document defines how gbdraw evaluates architecture-bearing changes. It is
the authority for ratchet definitions, review evidence, exceptions, and merge
decisions. `CLAUDE.md` owns the repository architecture, and
`gbdraw/web/CLAUDE.md` records the current Web runtime path and module owners.
Those documents link here when a change can alter an owner or path.

The normative CI lifecycle, independent Gate and Review results, workflow
ownership, and self-authorization boundaries are defined in
[`WEB_CHANGE_POLICY.md`](./WEB_CHANGE_POLICY.md). A Review result can require
human attention without changing the executable Gate result.

Architecture-to-behavior traceability, product-impact classification, Decision
Pack semantics, and product-decision eligibility are defined in
[`PRODUCT_IMPACT_RATCHET.md`](./PRODUCT_IMPACT_RATCHET.md). That policy does not
redefine semantic owner, canonical path, `OE`, `PE`, or `CB`, and it cannot
waive an architecture Gate failure.

The initial executable scope is the Web application. Human review applies to
architecture-bearing changes across the repository. The ratchet protects
convergence without treating fewer files or fewer lines as an architecture
goal.

## Architecture models and reflexion results

The ratchet compares declared intent with facts derived from source:

```text
Intended architecture model
  = this normative policy
  + inert architecture-rule registry
  + privileged path policy
  + exact accepted violations, when present

Source architecture model
  = trusted detector results for imports and registered evidence
  + executable ownership and canonical-entry markers
  + exact compatibility evidence

Reflexion result
  = pure comparison of intended authority with source facts
```

The intended model is split by concern. No file independently owns the entire
model. The source model contains only facts produced by registered detectors.
Heuristic names may be reported, but they are not evidence of semantic
ownership.

A reflexion result uses these source observations:

- `CONFORMING`: the required owner, path, or relation is present.
- `DIVERGENT`: source contains an unapproved owner, relation, path, or other
  registered violation.
- `ABSENT_REQUIRED`: a required owner, path, or relation is missing.

Frozen rules add baseline relations for accepted, new, and fixed violations.
The report fields and decisions are defined under [Enforcement modes](#enforcement-modes-and-results).

## Definitions

### Architecture fitness function

An architecture fitness function is a deterministic or explicitly reviewed
test of whether a selected architectural characteristic remains within its
intended constraints. Each fitness function must state its architectural
meaning, owner, scope, and decision rule.

A fitness function may be a hard invariant, an exact frozen-rule ratchet, a
differential metric, or a non-blocking trend diagnostic. A count alone is not
an architecture fitness function.

### Module

A module is a physical source file. Module count is not an architecture
maturity measure by itself.

A module may be added as the single owner of a new capability, as private
decomposition under an existing owner, or as a consolidation that removes the
superseded decision points in the same change.

### Semantic owner

A semantic owner is a location that independently determines at least one of
these properties:

- canonical representation or default value;
- validity, normalization, or migration;
- a state mutation invariant or resource lifecycle;
- equivalence or fallback selection;
- an authoritative transformation between representations.

A caller that invokes an owner without deciding these matters is not another
semantic owner.

### Privileged operator

A privileged operator is a location permitted to execute a sensitive
operation, such as constructing a Worker, replacing the current Result,
building the canonical render request, staging resource payloads, or admitting
SVG.

Several privileged operator locations can be valid while one semantic owner
defines the operation. The existing `allowedPrivilegedOwners` key in
`tools/web-change-policy.json` names an operator-location allowlist despite its
historical key name.

### Importer

An importer depends on another module's interface. High fan-in to one canonical
owner can be desirable. Uncontrolled privileged fan-out from feature modules is
not.

### Canonical path

A canonical path is the one supported production path by which behavior reaches
its authoritative boundary. Adapters may translate surface input, but they must
not create an independent schema, validation rule, normalization rule,
lifecycle, or fallback path.

### Compatibility path

A compatibility path is an active reader, migrator, branch, fallback, dual
schema, retired-field promotion, or protocol switch that exists only to accept
or preserve a superseded contract.

Persisted-format compatibility requires evidence that the old contract existed
in the first-parent history of `main` or in a release tag, a representative
positive fixture, a named schema namespace, and a removal condition when
removal is possible.

### Private decomposition

Private decomposition is a helper under an existing owner that introduces no
new mutable authority, validation, normalization, migration, lifecycle,
canonical path, privileged dependency, or unnecessary public export.

### Architecture-bearing change

A change is architecture-bearing when it adds, removes, moves, duplicates, or
redefines any of the following:

- a semantic owner or canonical path;
- a compatibility path;
- a privileged operator or importer;
- reactive source-of-truth state;
- a Worker or resource lifecycle;
- a parser or interpreter for an externally owned semantic system;
- a public contract or persisted representation;
- an abstraction that replaces or parallels an execution path.

Line count and file count do not decide whether a change is
architecture-bearing.

### Accepted deterministic violation

An accepted deterministic violation is one stable, exact violation identity
that is temporarily tolerated so a valuable rule can be enabled on an existing
base. It is not a privileged permission, metric snapshot, broad ignore,
compatibility contract, or approval of the current design.

### Active detector evidence contract

An active detector evidence contract is a versioned detector ID referenced by
an active rule. The ID fixes the subject category, source scope, positive and
negative classification semantics, and canonical subject grammar.

Detector code is reviewed trusted code. Fixtures and unchanged-base output
comparisons characterize it, but cannot prove semantic equivalence for
arbitrary source. During the initial rollout, an executable change to an active
detector or any transitive shared detection helper requires a new detector ID.

### Pure evaluator

The pure evaluator receives normalized source facts and already-loaded inert
authority as data. It returns deterministic schema, set, graph, and decision
results. It does not discover source facts or read, write, or own repository
authority.

## Reviewed architecture-debt model

For capability \(c\), semantic owner excess is:

\[
OE_c = \max(0, O_c - T_c)
\]

\(O_c\) is the number of independent repository locations that determine the
capability's semantics. \(T_c\) is the required number of owners, normally one.
When semantics are delegated to an external platform, \(T_c\) can be zero
inside the repository.

Total owner excess is:

\[
OE = \sum_c OE_c
\]

For behavior \(b\), canonical path excess is:

\[
PE_b = \max(0, P_b - 1)
\]

\(P_b\) is the number of independently meaningful production paths that
implement the behavior. Total path excess is:

\[
PE = \sum_b PE_b
\]

Compatibility burden is:

\[
CB = \text{number of active compatibility paths}
\]

The reviewed vector is:

\[
D_{reviewed} = (OE, PE, CB)
\]

Do not reduce this vector to a weighted score. A decrease in one component must
not hide an increase in another.

`OE`, `PE`, and `CB` are gbdraw-specific architecture fitness functions used in
review. Static analysis must not claim to calculate their complete
repository-wide values. A file-level detector cannot identify every independent
decision within a file, and an unregistered decision may have no deterministic
marker. Machine reports must name the observation they actually count, such as
registered definition locations or observed canonical entry edges.

### Exact changed-scope arithmetic

Full changed-scope sets and arithmetic are exception evidence. They are required
only when at least one of these conditions applies:

- `delta(OE) > 0`;
- `delta(PE) > 0`;
- `delta(CB) > 0` for a new persisted compatibility path;
- an accepted deterministic violation is added;
- multiple canonical owners or paths are intentionally retained; or
- a hard invariant receives a time-bounded waiver.

The last case records governance evidence; it does not turn a failing
deterministic check into a passing one. Authority and checker changes still
follow the separation and bootstrap order in `WEB_CHANGE_POLICY.md`.

For an exception, humans select the changed capabilities, behaviors,
compatibility namespaces, and complete before and after sets. Static analysis
may provide observations for that review, but it cannot infer complete
repository-wide `OE`, `PE`, or `CB` values.

Declare one owner-excess row for every changed capability. Each row has a stable
row ID, the semantic owner sets before and after, `O before`, `T before`, `O
after`, and `T after`. The declared sets and inputs must reproduce:

\[
OE_{before} = \max(0, O_{before} - T_{before})
\]

\[
OE_{after} = \max(0, O_{after} - T_{after})
\]

\[
\operatorname{delta}(OE) = OE_{after} - OE_{before}
\]

Declare one path-excess row for every changed behavior. Each row has a stable
row ID, the canonical production path sets before and after, `P before`, and `P
after`. The declared sets and inputs must reproduce `PE before`, `PE after`, and
`delta(PE)` using the path-excess formula above.

Declare one compatibility-burden row for every changed compatibility namespace.
Each compatibility path has a stable ID. The before and after ID sets must
reproduce `CB before`, `CB after`, and `delta(CB) = CB after - CB before`.

Sum every declared row for each component to obtain the exact changed-scope
totals:

```text
OE before -> OE after; delta(OE) = OE after - OE before
PE before -> PE after; delta(PE) = PE after - PE before
CB before -> CB after; delta(CB) = CB after - CB before
```

Every displayed value and delta is an integer. Multiple changed capabilities,
behaviors, or compatibility namespaces require multiple rows; sum those rows
explicitly. If an exception packet has no changed rows for a component, declare
`none` and record `0 -> 0; delta = 0` for that component.

## Required author evidence

An ordinary non-increasing change records concise decision evidence:

- what owner, path, or responsibility changes;
- its canonical location before and after;
- the superseded location removed, or the reason it remains;
- user-visible behavior verification;
- applicable deterministic checks; and
- rollback.

This concise form covers private decomposition, a one-for-one owner move, path
consolidation, deletion, and unchanged architecture when none of the exception
conditions above applies. It does not require an empty multi-row table or a
dedicated maintainer approval permalink.

For an `ARCHITECTURE_EXCEPTION`, the author records the exact changed-scope rows
and totals. The declaration includes:

- the capability or behavior and why the selected scope is complete;
- semantic owner sets, `O`, `T`, `OE`, and `delta(OE)` before and after for each
  changed capability;
- canonical production path sets, `P`, `PE`, and `delta(PE)` before and after for
  each changed behavior;
- compatibility path stable-ID sets, `CB`, and `delta(CB)` before and after for
  each changed compatibility namespace;
- the exact changed-scope totals obtained by summing all declared rows;
- superseded semantic owners, canonical production paths, and compatibility
  paths in three separate lists;
- the classification of each new module;
- the best non-exception alternative and why it was not selected;
- an expiry or concrete removal condition;
- persisted-compatibility, performance, scientific-output, and deterministic
  checker evidence when applicable; and
- an explicit maintainer decision on the exact proposed head.

Use `none` for an empty exception set. Do not replace complete sets with
unsupported repository-wide totals. `<= 0`, `non-positive`, or `yes` alone is
insufficient exception evidence. A favorable delta over an incomplete scope can
hide a parallel owner or path and is invalid.

## Merge acceptance

Architecture merge acceptance uses the independent Gate and Review axes:

```text
Gate = PASS only when every applicable deterministic condition passes
Review = REQUIRED for architecture-bearing or exception work
```

The executable Gate requires all applicable conditions below:

```text
intended behavior is verified by its test or output gate
AND newly introduced registered deterministic violations = 0
AND an implementation PR does not expand the accepted-violation set
AND first-party dependency cycles, including self-imports, = 0
AND unauthorized privileged surface = 0
AND new production dependencies = 0
AND gated performance or scientific-output regressions = 0
AND guard and authority separation remains intact
```

For an ordinary non-increasing change, human review confirms the concise
evidence, complete owner/path scope, and these bounds:

```text
delta(OE) <= 0
AND delta(PE) <= 0
AND delta(CB) <= 0
AND every superseded owner or path is removed in the same change
AND the declared changed-capability scope is complete
```

These inequalities define ordinary acceptance, not a mandatory table format.
Reviewers may ask for more detail when concise evidence does not establish the
bounds or completeness.

An exception may proceed only with the complete packet and explicit maintainer
decision described above. Positive debt, a new persisted compatibility path,
retained multiple owners or paths, or an accepted deterministic violation is
never authorized by a label, a passing size review, incomplete arithmetic, or
silence. A hard invariant remains Gate-failing until a separately authorized
and correctly sequenced executable authority change represents the disposition.

A new capability may begin with one new semantic owner without increasing owner
excess. It must also begin with one canonical path.

A change is architecture-convergent when supported behavior is preserved or
improved, the vector does not worsen, and at least one component decreases.

## Evidence examples

These examples use changed-scope sets. They do not assert repository-wide
totals.

### Ownership convergence

```text
Capability: render request construction
Owners before:
  {services/session-request.js, app/legacy-request.js}
O before: 2
T before: 1
OE before: 1
Owners after:
  {services/session-request.js}
O after: 1
T after: 1
OE after: 0
Superseded owner removed:
  app/legacy-request.js
Review result:
  OE 1 -> 0; delta(OE) = -1
```

Moving an owner can preauthorize the destination in inert authority, but the
runtime change must remove the old owner when it activates the new one.

### Canonical path switch

```text
Behavior: Generate reaches render-request construction
Paths before:
  {app/run-analysis.js -> services/session-request.js}
P before: 1
PE before: 0
Paths after:
  {app/generate.js -> services/session-request.js}
P after: 1
PE after: 0
Superseded path removed:
  app/run-analysis.js -> services/session-request.js
Review result:
  PE 0 -> 0; delta(PE) = 0
```

The authority may list both edges during a staged move. Source must still have
exactly one active edge.

### Compatibility contraction

```text
Schema namespace: saved-session
Compatibility paths before:
  {saved-session:v4-to-v5}
CB before: 1
Compatibility paths after:
  none
CB after: 0
Removed path:
  saved-session:v4-to-v5
Review result:
  CB 1 -> 0; delta(CB) = -1
```

A new compatibility ID requires the persisted-format evidence named in the
definition. A branch-only intermediate format is rewritten to the current
format instead of gaining a migration.

### Exact accepted-violation contraction

```text
Accepted subjects before:
  {canonical-path.example:app/old.js->services/example.js}
Observed violations after:
  none
Accepted subjects after:
  none
Authority change:
  remove exactly the fixed subject in the same PR
```

Keeping the fixed subject is stale authority and fails. Removing it while its
violation remains is an invalid authority change and fails.

## Authority and trust boundaries

| Concern | Owner |
| --- | --- |
| Definitions, formulas, exceptions, and acceptance | `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` |
| Repository architecture and cross-surface convergence | `CLAUDE.md` |
| Current Web runtime path and module ownership | `gbdraw/web/CLAUDE.md` |
| Executable registered source-fact detection | `tools/web-architecture-detectors.mjs` |
| Pure schema, set, graph, and decision mechanics | `tools/web-architecture-evaluation.mjs` |
| Intended owner and canonical-path rules | `tools/web-architecture-rules.json` |
| Exact accepted violations, when enabled and nonempty | `tools/web-architecture-violations.json` |
| Privileged operator and importer permissions | `tools/web-change-policy.json` |
| Git and file I/O, trusted orchestration, and reporting | `tools/check-web-change-budget.mjs` |
| Cross-cutting executable contracts | `tests/web/architecture-contracts.test.mjs` |
| Fixture-heavy ratchet mechanics | `tests/web/architecture-ratchet-fixtures.test.mjs` |
| Human change declaration | `.github/pull_request_template.md` |
| Architecture-exception approval | Structured manual maintainer review |

`tools/check-web-change-budget.mjs` is the sole Web CLI and CI entry point. A
single entry point does not require one monolithic implementation module. The
detector owns executable source evidence; the evaluator owns only pure
mechanics; inert JSON owns intended paths, counts, modes, and exceptions.

Detector code must not contain intended paths, allowed edges, counts,
enforcement modes, baseline eligibility, accepted violations, compatibility
authorization, or privileged permission arrays. The evaluator must not contain
source scanners, filesystem access, Git calls, environment reads, reporting, or
intended authority.

## Trusted-base evaluation

Runtime conformance is evaluated with trusted target-branch artifacts:

```text
base detector implementation
+ base architecture rules
+ base privileged policy
+ proposed head source facts derived by the base detector
```

Candidate head rules are read only as inert JSON for schema and authority-delta
validation. The trusted-base workflow never imports or executes a proposed head
detector. A head detector, rule, policy, workflow, or accepted-violation change
cannot authorize runtime changed by the same pull request.

Trusted evaluation proceeds in this order:

1. Load the base checker, detector catalog, rules, and privileged policy.
2. Derive base and head source facts with the base detector implementation.
3. Evaluate runtime conformance against base authority.
4. Parse candidate head authority as inert data.
5. Classify the authority direction and validate it against source evidence.
6. Derive the decision and report it in stable rule and subject order.

An authority expansion, relaxation, tightening, contraction, or
reinterpretation is reviewed separately from runtime observation.

## Initial rule schema

The initial registry supports at most three rules and exactly two discriminated
rule kinds:

- `single-semantic-owner`
- `single-canonical-entry-edge`

The cap and kinds belong to schema version 1. Raising the cap or adding a kind
requires an evidence-backed schema plan. The registry is inert JSON. It must not
contain regular expressions, JavaScript, executable expressions, wildcard
paths, detector-narrowing symbols, or a generic catch-all rule kind.

For `single-semantic-owner`, `exactDefinitionCount` is `1`.
`allowedDefinitionPaths` may include an inactive preauthorized destination, but
source must contain exactly one active definition.

For `single-canonical-entry-edge`, `exactActiveEdgeCount` is `1`. Each allowed
edge has a normalized `from` and `to` path. The detector emits the complete set
of observed edges in its versioned rule-specific scope. Conformance is:

```text
observedEdges is a subset of allowedEdges
AND size(observedEdges) = 1
```

Zero observed edges is `ABSENT_REQUIRED`. An unapproved edge is `DIVERGENT`.
Two active edges are `DIVERGENT`, even when both were preauthorized. An allowed
but inactive edge does not fail. This staging supports an authority-only
preauthorization, one atomic runtime switch, and a later removal of the inactive
old edge without permitting parallel canonical paths.

First-party dependency-cycle enforcement is separate and does not consume a
registry slot.

## Enforcement modes and results

Rules use one of three modes:

- `HARD`: a clean invariant; any divergent or absent-required observation
  fails.
- `REPORT_ONLY`: divergent or absent-required observations are reported but do
  not fail by themselves.
- `FROZEN`: exact accepted base violations pass, new violations fail, and fixed
  authority must contract.

For each rule and relevant subject, the checker emits five independent fields:

```text
observation:
  CONFORMING
  DIVERGENT
  ABSENT_REQUIRED

mode:
  HARD
  FROZEN
  REPORT_ONLY

baselineRelation:
  NOT_APPLICABLE
  ACCEPTED
  NEW
  FIXED

authorityResolution:
  NOT_APPLICABLE
  RETAINED
  EXACT_CONTRACTION
  INVALID_CHANGE

decision:
  PASS
  REPORT
  FAIL
```

`decision` is derived by one pure function from the first four fields. It is not
assigned or stored independently.

| Mode and observation | Baseline relation | Authority resolution | Decision |
| --- | --- | --- | --- |
| `HARD` + `CONFORMING` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `HARD` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `FAIL` |
| `REPORT_ONLY` + `CONFORMING` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `REPORT_ONLY` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `REPORT` |
| `FROZEN` + `CONFORMING`, with no fixed accepted subject | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `ACCEPTED` | `RETAINED` | `PASS` |
| `FROZEN` + active accepted violation after its entry was removed | `ACCEPTED` | `INVALID_CHANGE` | `FAIL` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `NEW` | `NOT_APPLICABLE` | `FAIL` |
| `FROZEN` + fixed accepted subject | `FIXED` | `RETAINED` | `FAIL` |
| `FROZEN` + fixed accepted subject | `FIXED` | `EXACT_CONTRACTION` | `PASS` |
| `FROZEN` + invalid store expansion, removal, or mismatch | as derived | `INVALID_CHANGE` | `FAIL` |

Any other combination is an internal checker error and fails closed. Hard and
report-only rules never consult accepted violations and always use
`NOT_APPLICABLE` for baseline relation and authority resolution.

The checker reports measurable inventories as:

```text
before
added
removed
after
delta
```

New modules, exports, watchers, reactive declarations, and heuristic names stay
report-only unless a separately registered deterministic rule gives them an
architectural meaning and decision rule.

## Frozen-rule protocol

Frozen mechanics are conditional. Phase 0 must select at least one materially
valuable deterministic rule with exact existing debt. Otherwise frozen
activation fails closed, accepted-store parsing and evaluation remain absent,
and `tools/web-architecture-violations.json` must not exist.

Phase 0 on 2026-08-18 selected no existing-debt rule. The initial rollout
therefore keeps frozen mechanics unavailable and the accepted-violation store
absent. Adding the first later existing-debt rule requires a new evidence-backed
checker plan.

When frozen mechanics are available, for rule \(r\):

- \(V_r(x)\) is the exact set of violations in revision \(x\).
- \(A_r(x)\) is the exact accepted set recorded for revision \(x\).

```text
New violations:
  N_r = V_r(head) \ A_r(base)

Fixed accepted violations:
  F_r = A_r(base) \ V_r(head)

Contracted fixed authority:
  C_r = F_r \ A_r(head)

Retained stale authority:
  R_r = F_r intersection A_r(head)

Invalid removal of a still-active accepted subject:
  U_r = (A_r(base) intersection V_r(head)) \ A_r(head)
```

A normal implementation PR passes a store change only when:

```text
N_r = empty
A_r(head) = A_r(base) \ F_r
C_r = F_r
R_r = empty
U_r = empty
A_r(head) \ A_r(base) = empty
```

CI never creates, expands, updates, or refreezes the store. Each stored entry
must have a registered rule key, a canonical stable subject without a line
number, a rationale, tracking evidence, and a removal condition when known.
Duplicates, wildcards, regular expressions, broad file exceptions, metric
totals, and privileged permissions are invalid. Delete the file when its final
entry contracts. When frozen mechanics are available, absence of the file means
an empty accepted set.

Accepted-set expansion, when supported, is an authority-only recording decision
for exact untouched-base observations from an already merged eligible
report-only rule. A later authority-only change tightens that rule to frozen.

## Detector and rule onboarding

Every new rule follows a detector-first sequence:

1. Add a versioned rule-specific detector in a checker-only PR. Include a fixed
   fixture corpus and untouched-base characterization. Do not change authority.
2. Add the rule to the inert registry in a later authority-only PR.
3. Activate only the mode supported by untouched-base evidence.

For a clean-base rule, activate `hard` and do not create a violation-store
entry.

For an existing-debt rule after frozen mechanics are available:

1. Merge its detector.
2. Add the rule as `report-only`.
3. Record exact untouched-base observations in a separate authority-only PR.
4. Tighten the rule to `frozen` in a later authority-only PR.
5. Remove each accepted entry in the same runtime PR that fixes its subject.

An active detector ID is immutable as an executable evidence contract. To
change its semantics:

1. Add `*.v2` beside `*.v1` in a checker-only PR and characterize both.
2. Migrate inert authority to `*.v2` in a separate authority-only PR after
   untouched-base validation.
3. Remove the now-unreferenced `*.v1` in a later checker-only PR.

Do not modify a shared helper transitively used by active v1 while adding v2.
Use a version-specific helper when needed. If the canonical subject grammar
changes, migrate rule references and affected accepted subjects together with
an exact old-to-new mapping that neither adds nor drops debt.

## Authority directions

| Change | Direction and condition |
| --- | --- |
| Add a rule after its detector exists | expansion or activation |
| Remove a rule | relaxation, not contraction |
| Change a detector ID, rule kind, or subject encoding | reinterpretation through an already merged versioned detector |
| Add an allowed definition path or edge | expansion or preauthorization |
| Change an exact count | unsupported in schema version 1 |
| `hard` to `frozen` or `report-only` | relaxation |
| `frozen` to `report-only` | relaxation |
| `baselineEligible: false` to `true` | relaxation |
| Remove an allowed definition path | contraction only when head source does not require it |
| Remove an allowed edge | contraction only when the edge is inactive in head source |
| `report-only` to `frozen` or `hard` | tightening only when the destination gate passes |
| `frozen` to `hard` | tightening only when no accepted violation remains |
| `baselineEligible: true` to `false` | tightening only when no accepted violation remains |

An implementation PR cannot combine detector mechanics with authority
activation, expansion, relaxation, tightening, or reinterpretation. There is no
bootstrap exception.

## Manual architecture review

Ordinary architecture-bearing changes use the concise evidence and normal PR
review described above. They do not require a separate reviewer-comment
template or a dedicated approval permalink.

For an `ARCHITECTURE_EXCEPTION`, the author or automated agent prepares a
maintainer review packet before requesting a decision. It contains the exact
head SHA, complete before and after evidence, arithmetic, alternative,
expiry/removal condition, completed CI results, known limitations, and a fully
populated, ready-to-post version of the structured comment below. The maintainer
reviews the packet and then posts, edits, or rejects the comment manually as
their own decision. Preparing the packet is not approval, and an automated agent
must not post the decision.

The architecture owner reviews the final exception head and posts this
structured comment manually:

```markdown
Architecture decision: APPROVED
Reviewed head SHA: <full commit SHA>
Reviewed changed-scope rows:
- OE row <ID>: owners before <set> -> owners after <set>; O before <integer> -> O after <integer>; T before <integer> -> T after <integer>; OE before <integer> -> OE after <integer>; delta(OE) = <integer>
- PE row <ID>: paths before <set> -> paths after <set>; P before <integer> -> P after <integer>; PE before <integer> -> PE after <integer>; delta(PE) = <integer>
- CB row <ID>: stable IDs before <set> -> stable IDs after <set>; CB before <integer> -> CB after <integer>; delta(CB) = <integer>
Changed-scope totals:
- OE before <integer> -> OE after <integer>; delta(OE) = <integer>
- PE before <integer> -> PE after <integer>; delta(PE) = <integer>
- CB before <integer> -> CB after <integer>; delta(CB) = <integer>
Superseded semantic owners:
- <none, or exact set>
Superseded canonical production paths:
- <none, or exact set>
Superseded compatibility paths, with stable IDs:
- <none, or exact stable-ID set>
Scope completeness decision:
- <why the declared changed-scope rows and sets are complete>
Persisted-compatibility exception:
- <none, or exact exception and evidence>
Limitations:
- <none, or explicit non-authorizing limitation>
Expiry or removal condition:
- <date, release, issue, or measurable condition>
```

The exception section of the PR template stores a permalink to the decision.
`Reviewed head SHA` must equal the pull request's current head SHA. Any later
commit invalidates the exception decision and requires a new comment.

An automated agent may prepare evidence but must not post the approval. The
author declaration cannot substitute for the exception decision. If the sole
architecture maintainer is also the author, that maintainer may post the
separate decision after reviewing the final head.

This is an exception review requirement, not an ordinary CI status. Repository
tests may protect template anchors, but do not inspect GitHub comments, comment
authors, or live head SHAs. Workflow presence and branch protection do not prove
that an exception decision exists or matches the merge candidate.

## External precedents

These sources explain the adopted terminology and mechanisms. They are neither
repository authorities nor dependencies.

- Thoughtworks' [architectural fitness function](https://www.thoughtworks.com/en-us/radar/techniques/architectural-fitness-function)
  provides the continuous integrity-check concept.
- Murphy, Notkin, and Sullivan's [software reflexion model](https://www.cs.ubc.ca/~murphy/papers/rm/fse95.html)
  provides the intended-model, source-model, and comparison separation.
- ArchUnit's [freezing rules](https://www.archunit.org/userguide/html/000_Index.html#_freezing_arch_rules)
  provide the accepted-existing, new, and fixed violation pattern. gbdraw keeps
  exact identities and forbids CI refreezing.
- SonarQube's [Clean as You Code](https://docs.sonarsource.com/sonarqube-server/2025.2/user-guide/clean-as-you-code/about-new-code/)
  and [NDepend quality gates](https://www.ndepend.com/docs/quality-gates) provide
  the differential new-debt principle. gbdraw does not import their scoring
  models.
- [ArchSteer](https://www.archsteer.com/) provides a product-level analogue for
  declared, distributed, and actual architecture. gbdraw keeps governance
  repository-local.
- [dependency-cruiser](https://github.com/sverweij/dependency-cruiser) and
  [CodeScene change coupling](https://codescene.io/docs/guides/technical/change-coupling.html)
  inform dependency and historical diagnostics. The initial ratchet adds no
  external governance dependency and keeps historical coupling non-blocking.
