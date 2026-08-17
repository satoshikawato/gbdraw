# gbdraw architecture fitness-function ratchet implementation plan

Status: revised proposal, revision 5
Created: 2026-08-17
Revised: 2026-08-18
Supersedes: revision 4 of this plan and `GBDRAW_ARCHITECTURE_MATURITY_RATCHET_IMPLEMENTATION_PLAN_2026-08-17.md` revision 1
Work-branch base: latest merged `origin/dev`
Pull-request target and integration branch: `dev`
Release branch: `main`, promoted from `dev` under the repository's promotion checks
Target repository path for this plan: `docs/internal/GBDRAW_ARCHITECTURE_FITNESS_FUNCTION_RATCHET_IMPLEMENTATION_PLAN_2026-08-17.md`
Intended permanent authority after implementation: `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`

## 1. Objective

Implement a repository-local **architecture fitness-function ratchet** that makes architectural convergence permanent, reviewable, and partially machine-enforced.

“Architecture maturity” is the desired outcome. The mechanism is a set of project-specific architectural fitness functions evaluated differentially against the target branch. The ratchet must not equate maturity with fewer modules or fewer lines of code. It must prevent growth in the number of independent places that define the same meaning, realize the same canonical behavior, preserve obsolete representations, or exercise privileged operations without an explicit architecture decision.

The implementation must establish six layers:

1. one normative architecture fitness-function document;
2. trusted executable mechanics split into a detector catalog that derives registered source facts and a pure evaluator that classifies those facts; neither contains owner, threshold, enforcement, or exception authority;
3. one inert, machine-readable intended-architecture rule registry for the critical Web capabilities already enforced by repository policy;
4. an optional, version-controlled store of exact accepted deterministic violations, used only when a valuable rule cannot initially pass on the existing base;
5. deterministic CI checks for the portions that can be measured reliably;
6. a pull-request declaration and explicit maintainer review that record the architectural judgment static analysis cannot make safely.

The system must use a differential “no new architecture debt” rule rather than requiring an immediate cleanup of all historical debt. It must also make contractions visible and prevent a fixed violation from being silently reintroduced.

This initiative is a guard-and-process change. It must not re-architect Web runtime behavior or modify `gbdraw/web/js/` production modules.

## 2. External precedents and explicit adoption decisions

This design uses established architecture-governance ideas but keeps the implementation repository-local.

### 2.1 Architectural fitness functions

Thoughtworks describes an architectural fitness function as an objective integrity assessment that continuously validates and preserves selected architectural characteristics. gbdraw adopts that term for deterministic checks such as owner-location constraints, canonical call-path contracts, dependency-cycle checks, privileged-surface allowlists, and performance/parity gates.

The normative document must call `OE`, `PE`, and `CB` **gbdraw-specific architectural fitness functions**. They are project-specific measures, not proposed industry standards.

### 2.2 Software reflexion models

Software reflexion models compare an intended high-level model with a source-derived model and report where they conform or diverge. gbdraw adopts that separation:

```text
Intended architecture model
  = normative rule
  + inert architecture-rule registry
  + privileged path policy
  + explicit accepted-violation store, when present

Source architecture model
  = trusted detector implementation applied to imports
  + executable ownership markers
  + canonical entry-point markers
  + reactive and lifecycle evidence
  + exact compatibility evidence

Reflexion result
  = pure evaluation of intended authority and source facts
  + conforming
  + divergent
  + absent
  + accepted existing violation
  + newly introduced violation
  + fixed violation
```

The checker must not claim to derive the complete architecture automatically. It derives only the source facts registered with deterministic evidence.

### 2.3 Frozen-rule ratchets

ArchUnit's `FreezingArchRule` records existing violations, reports only new violations on later runs, and removes fixed violations from the known set so they cannot regress. gbdraw adopts this behavior conceptually for exact deterministic architecture violations.

The gbdraw implementation must differ in two important ways:

- CI never writes, creates, expands, or “refreezes” the accepted-violation store;
- the store contains stable exact violation identities, not line numbers, metric totals, broad regular expressions, or heuristic name matches.

If the current base has no accepted deterministic violations, no store file is created. Absence of the file means an empty accepted set.

### 2.4 Differential quality gates

SonarQube's Clean as You Code model and NDepend's baseline-aware quality gates focus enforcement on new or worsened debt rather than making every existing defect block all work. gbdraw adopts the same differential principle:

```text
existing explicitly accepted debt may remain temporarily
new architecture debt fails
removed debt cannot be reintroduced
```

The implementation must not import SonarQube, NDepend, or their scoring models.

### 2.5 Agent-oriented architecture governance

ArchSteer is a close product-level analogue: it distinguishes declared, distributed, and actual architecture, applies net-new conformance gates, and distributes rules to coding agents. gbdraw adopts the conceptual separation:

```text
Declared   = normative document and inert rule registry
Distributed = AGENTS.md, CLAUDE.md, Web guidance, PR template
Actual      = checker-derived source model and CI report
```

gbdraw must not add ArchSteer or another external architecture-governance dependency. The existing repository checker remains the single enforcement entry point; inert rule and policy files own intended authority.

### 2.6 Dependency and historical diagnostics

Dependency-cruiser demonstrates declarative JavaScript dependency rules, while CodeScene demonstrates change-coupling analysis from version history. gbdraw adopts only the relevant techniques:

- first-party import cycles and selected dependency boundaries are deterministic merge gates;
- temporal coupling, change amplification, and rework rate are deferred trend diagnostics;
- no new npm, Python, SaaS, or runtime dependency is added for this initiative.

### 2.7 Adoption summary

| Existing idea | Adopt | Do not adopt |
| --- | --- | --- |
| Architectural fitness functions | Continuous project-specific integrity checks | A universal architecture score |
| Reflexion models | Intended/source separation and orthogonal result fields | Claims of complete automatic architecture recovery |
| ArchUnit freezing | Exact accepted violations, no-new regression, contraction | CI auto-write, broad refreeze, line-number identities |
| Clean as You Code / baseline gates | Differential enforcement against the target branch | Generic maintainability scores as architecture truth |
| Agent governance | Short rules distributed to agent entry points | A second external policy owner |
| Dependency analysis | Cycles and explicit forbidden/required boundaries | A new dependency-cruiser installation |
| Change coupling | Historical diagnostic after enough PR history | Immediate per-PR hard gate |

## 3. Existing authority to preserve

This plan operationalizes existing repository principles rather than introducing a competing architecture.

The current repository already requires the following:

- treat every added line as additional maintenance burden;
- optimize for fewer change points, branches, and duplicated behavior rather than raw line counts;
- add an abstraction only when it unifies at least two real execution paths and removes the superseded paths in the same change;
- converge CLI, Python API, Web UI, saved sessions, and replay on canonical typed boundaries;
- normalize compatibility input once and do not write retired fields into current formats;
- maintain explicit module ownership and privileged-operation allowlists;
- keep Web runtime changes separate from guard and policy changes.

The new ratchet must refer to those authorities and clarify their common rule. It must not duplicate all current Web ownership tables or runtime data-flow documentation.

## 4. Problem statement

The current controls are useful but incomplete.

### 4.1 Terminology is overloaded

`module`, semantic owner, privileged operator, and importer are sometimes discussed as if they were the same thing.

They are not:

- a **module** is a physical source file;
- a **semantic owner** decides what a value, state, contract, normalization, or invariant means;
- a **privileged operator** is a call site permitted to execute a sensitive operation;
- an **importer** depends on another module's interface.

The existing `allowedPrivilegedOwners` JSON key is operationally an allowlist of privileged operator locations. It must remain unchanged in this initiative to avoid a policy-schema migration, but the permanent documentation must state the distinction explicitly.

### 4.2 The checker reports growth more readily than contraction

The current Web change checker reports new modules, exports, `create*` declarations, reactive state, watchers, resource-like names, compatibility-like names, session object keys, dependencies, privileged operators, and privileged importers.

For several categories it reports only additions. It does not consistently display:

- total before;
- additions;
- removals;
- total after;
- delta.

A contraction therefore receives less visible evidence than an expansion.

### 4.3 Critical capability definitions are duplicated

The checker and `tests/web/architecture-contracts.test.mjs` contain overlapping descriptions of privileged capabilities, import targets, and detection patterns.

Those definitions can drift. Executable detection needs one implementation owner, intended semantic-owner and canonical-path constraints need a separate inert authority owner, and the allowlist of permitted privileged paths must remain in `tools/web-change-policy.json`. A single module must not combine regexes or scan functions with kind-specific intended paths, edges, counts, enforcement modes, or baseline eligibility.

### 4.4 Semantic maturity cannot be inferred safely from names alone

A new module named `normalizer`, `manager`, or `parser` may be correct or incorrect. A module count cannot determine ownership quality.

Likewise, a regex cannot reliably decide whether a helper merely observes browser-resolved SVG state or independently reimplements SVG, CSS, geometry, or compatibility semantics.

The permanent rule therefore needs both machine-enforced contracts and explicit human evidence.

## 5. Definitions

The normative document created by this plan must define the following terms exactly once.

### 5.1 Architecture fitness function

A deterministic or explicitly reviewed test of whether a selected architectural characteristic remains within its intended constraints.

A fitness function may be:

- a hard invariant that must be zero;
- a frozen-rule ratchet that permits only explicitly accepted existing violations;
- a differential metric that must not worsen;
- a trend diagnostic that informs later architecture work but does not block an individual PR.

A metric is not an architecture fitness function merely because it can be counted. It must have a stated architectural meaning, owner, scope, and decision rule.

### 5.2 Module

A physical source file. Module count is not an architecture maturity metric by itself.

A module may be added when it is one of the following:

1. the single owner of a real new capability;
2. private decomposition under an existing owner;
3. consolidation that replaces multiple existing decision points and removes the superseded implementations in the same change.

### 5.3 Semantic owner

A location that independently determines at least one of the following:

- canonical representation;
- default value;
- validity;
- normalization or migration;
- state mutation invariant;
- resource lifecycle;
- equivalence;
- fallback selection;
- authoritative transformation between representations.

A caller that invokes an owner without independently deciding these matters is not another semantic owner.

### 5.4 Privileged operator

A location that is authorized to perform a sensitive operation, such as constructing a Worker, replacing the current Result, building the canonical render request, staging resource payloads, or admitting SVG.

Several operator locations may be legitimate even when there is only one semantic owner.

### 5.5 Importer

A module that imports an owner or privileged service. High fan-in to a canonical owner is often desirable. Uncontrolled privileged fan-out from feature modules is not.

### 5.6 Canonical path

The one supported production path by which a behavior reaches its authoritative boundary.

Adapters may translate surface-specific input, but they must not create an independent schema, validation rule, normalization rule, lifecycle, or fallback path.

### 5.7 Compatibility path

Any active branch, reader, migrator, fallback, dual schema, retired-field promotion, or protocol switch that exists only to accept or preserve a superseded contract.

Persisted-format compatibility is allowed only under the repository's existing evidence rules.

### 5.8 Private decomposition

A helper module under an existing owner that does not introduce any of the following:

- new mutable authority;
- independent validation;
- independent normalization or migration;
- independent lifecycle;
- a second canonical path;
- a new privileged dependency;
- an unnecessary public export.

### 5.9 Architecture-bearing change

A change is architecture-bearing when it adds, removes, moves, duplicates, or redefines any of the following:

- a semantic owner;
- a canonical path;
- a compatibility path;
- a privileged operator or importer;
- reactive source-of-truth state;
- a Worker or resource lifecycle;
- a parser or interpreter for an externally owned semantic system;
- a public contract or persisted representation;
- an abstraction that replaces or parallels an execution path.

Line count and file count alone do not make a change architecture-bearing.

### 5.10 Intended architecture model

The declared architecture against which source evidence is compared. For the initial Web scope it consists of:

- the normative ratchet;
- the inert architecture-rule registry;
- current Web ownership documentation;
- privileged operator/importer policy;
- exact accepted deterministic violations, if any.

No single artifact independently owns the entire model. Each artifact owns one distinct concern listed in the authority table.

### 5.11 Source architecture model

The deterministic source facts extracted from the selected repository revision by the trusted detector implementation, including registered imports, executable owner markers, canonical call markers, dependency graph edges, and exact compatibility evidence.

Heuristic names may be reported alongside the source model but are not authoritative facts about semantic ownership.

### 5.12 Reflexion result

The comparison result between intended and source architecture:

- **conforming**: intended relation or owner is present exactly as required;
- **divergent**: source contains an unapproved relation, owner, path, or violation;
- **absent**: an intended required owner, path, or relation is missing;
- **accepted existing violation**: an exact source violation is present in the approved violation store;
- **new violation**: a source violation is absent from the approved store;
- **fixed violation**: an approved violation is no longer present in source and should be removed from the store.

### 5.13 Accepted deterministic violation

A stable, exact violation identity that is temporarily tolerated so a valuable architecture rule can be enabled on an existing codebase.

It is not:

- a permission to perform a privileged operation;
- a metric snapshot;
- a regex ignore pattern;
- a general exception for a file;
- a compatibility contract;
- evidence that the architecture is desirable.

Every accepted violation requires a rule key, stable subject identity, rationale, tracking issue or equivalent evidence, and removal condition when known.

### 5.14 Active detector evidence contract

A versioned detector ID referenced by any active rule. The ID fixes the detector's subject category, match scope, positive and negative classification semantics, and canonical subject grammar. A semantically different detector is a new contract and therefore requires a new ID.

Detector implementation is trusted checker code. The ratchet can prevent detector code from authorizing runtime or authority in the same PR, but it cannot mechanically prove semantic equivalence for arbitrary executable changes. Fixtures and untouched-base output comparisons are characterization evidence, not proof. During the initial rollout, any executable change to an active detector or its transitive shared detection helpers requires a new versioned ID. Same-ID edits are limited to comments, formatting, and other changes that leave executable tokens unchanged.

### 5.15 Pure evaluator

Authority-free executable mechanics that receive normalized source facts and already-loaded inert authority as data, then return deterministic validation, set, graph, and decision results. The evaluator neither discovers source facts nor reads, writes, or owns repository authority.

## 6. Normative maturity model

### 6.1 Semantic owner excess

For capability \(c\):

\[
OE_c = \max(0, O_c - T_c)
\]

where:

- \(O_c\) is the number of independent repository locations that determine the semantics of capability \(c\);
- \(T_c\) is the required number of semantic owners.

For most gbdraw capabilities, \(T_c = 1\).

When semantics are deliberately delegated to an external platform, \(T_c = 0\) inside the repository. For example, browser-resolved SVG and CSS semantics must not be reimplemented by a Node test harness merely to compare output.

Total semantic owner excess is:

\[
OE = \sum_c OE_c
\]

### 6.2 Canonical path excess

For behavior \(b\):

\[
PE_b = \max(0, P_b - 1)
\]

where \(P_b\) is the number of independently meaningful production paths that implement behavior \(b\).

Total canonical path excess is:

\[
PE = \sum_b PE_b
\]

### 6.3 Compatibility burden

\[
CB = \text{number of active compatibility paths}
\]

This includes supported migrations as well as accidental fallbacks. A justified persisted-format migration may increase `CB` temporarily, but it requires repository-history evidence, a positive fixture, a named schema namespace, and an explicit removal condition when removal is possible.

### 6.4 Reviewed architecture-debt vector

\[
D_{reviewed} = (OE, PE, CB)
\]

Do not collapse this vector into a weighted scalar score. A single score would permit a reduction in one dimension to hide a harmful increase in another.

`OE`, `PE`, and `CB` remain gbdraw-specific architecture fitness functions, but they are review concepts derived from explicit before/after sets for the changed capability scope. CI must not claim to calculate repository-wide `OE`, `PE`, or `CB`. File-level matches cannot distinguish two independent decisions hidden in one file, and an unregistered semantic decision may not have a deterministic marker.

Machine output must use names that state the actual observation, such as:

```text
Registered Authority Location Count
Registered Canonical Contract Violations
Registered Compatibility Path Count
```

For an architecture-bearing change, the author must enumerate the concrete before/after owner, canonical-path, and compatibility-path sets and explain why the declared capability scope is complete. Bare totals such as `OE before = 3` are not sufficient evidence.

### 6.5 Architecture acceptance rule

Merge acceptance has two independent layers:

```text
MERGEABLE =
  DETERMINISTIC_GATE_PASSED
  AND ARCHITECTURE_REVIEW_PASSED
```

`DETERMINISTIC_GATE_PASSED` means:

```text
intended behavior is verified by the applicable test or output gate
AND registered deterministic architecture violations newly introduced = 0
AND the accepted-violation set does not expand in an implementation PR
AND first-party dependency cycles remain zero, including self-import cycles
AND unauthorized privileged surface remains zero
AND no production dependency is added
AND no unapproved performance or scientific-output regression is introduced where a gate exists
AND guard and authority separation remains intact
```

`ARCHITECTURE_REVIEW_PASSED` means the pull request records, for every changed capability, the concrete before/after sets and a maintainer explicitly confirms:

```text
delta(OE) <= 0
AND delta(PE) <= 0
AND (
      delta(CB) <= 0
      OR the persisted-compatibility exception is satisfied
    )
AND every superseded owner or path is removed in the same change
AND the declared changed-capability scope is complete
```

The automated agent or pull-request author may prepare the evidence but cannot treat the author declaration as architecture approval. Approval means the structured PR comment defined in PR A, posted manually by the maintainer who owns the repository architecture and bound to the exact reviewed head SHA. An automated agent must not post it. If the sole architecture maintainer is also the PR author, the separate explicit architecture-owner comment remains required, but an independent second maintainer is not required.

For a real new capability, one new semantic owner may be added without increasing owner excess. The new capability must begin with one owner and one canonical path.

A change is **architecture-convergent** when behavior and supported capability are preserved or improved, the vector does not worsen, and at least one debt component decreases.

### 6.6 Differential architecture-debt rule

The ratchet compares the proposed head to the target-branch base. It does not require an unrelated historical cleanup in every PR.

This section defines the frozen-rule contract, but the initial rollout implements it only when Phase 0 selects at least one materially valuable existing-debt rule. Otherwise PR F2, PR G1, and PR G2 are skipped, `frozen` activation fails closed as unavailable, the violation store remains absent, and the design remains documented for a later evidence-backed checker-only plan.

The accepted-violation set algebra applies only to a rule whose enforcement mode is `frozen`. For a frozen deterministic rule \(r\):

- \(V_r(x)\) is the exact set of violations detected in revision \(x\);
- \(A_r(x)\) is the exact set of accepted violations recorded for revision \(x\).

New violations are:

\[
N_r = V_r(head) \setminus A_r(base)
\]

Fixed accepted violations are:

\[
F_r = A_r(base) \setminus V_r(head)
\]

For the proposed head authority, classify fixed entries as:

\[
C_r = F_r \setminus A_r(head)
\]

\[
R_r = F_r \cap A_r(head)
\]

where \(C_r\) is contracted fixed authority and \(R_r\) is retained stale authority. Removing an accepted entry while its violation still exists is invalid:

\[
U_r = (A_r(base) \cap V_r(head)) \setminus A_r(head)
\]

The merge gate requires \(N_r = \varnothing\). A normal implementation PR that changes the store passes only when all of the following are true:

```text
A_r(head) = A_r(base) \ F_r
C_r = F_r
R_r = empty
U_r = empty
A_r(head) \ A_r(base) = empty
```

This is the exact same-PR safe-contraction path: every accepted subject fixed by the head is removed, no still-active subject is removed, and no authority is added. A fixed entry that remains in the head store is stale authority and fails. When the final entry is contracted, the store file is deleted rather than committed as an empty object.

For a `report-only` rule, the checker reports `V_r(head)` as signals and does not compare it with the accepted-violation store:

```text
signals = V_r(head)
decision = PASS when signals is empty, otherwise REPORT
```

A `hard` rule is never baseline-eligible and never consults the accepted-violation store.

An implementation PR may not add entries to \(A_r\). A later authority-only recording PR may expand \(A_r\) only for exact untouched-base signals from an already-merged eligible report-only rule, as specified in PR G1. The evaluator reports the authority-delta classification separately from source observation so a source fix and its exact store contraction can pass atomically without allowing an unrelated authority edit.

### 6.7 Frozen accepted-violation protocol

The optional store path is:

```text
tools/web-architecture-violations.json
```

Absence means an empty accepted set.

If present, the schema must contain stable exact entries, for example:

```json
{
  "schemaVersion": 1,
  "acceptedViolations": [
    {
      "rule": "canonical-path.example-rule",
      "subject": "app/example.js -> services/example.js",
      "tracking": "issue-or-ADR-reference",
      "rationale": "Existing violation exposed when the rule was introduced.",
      "removalCondition": "Delete after callers converge on services/example.js."
    }
  ]
}
```

Requirements:

- `rule` must reference a registered deterministic rule key;
- `subject` must be a stable rule-specific identity and must not contain a source line number;
- the already-merged detector for the rule must provide one canonical subject encoder and grammar, and the stored subject must round-trip through it exactly;
- duplicate entries are invalid;
- wildcard, catch-all, and regex subjects are invalid;
- CI must never create or update the file automatically;
- a normal implementation PR cannot expand the accepted set;
- an accepted-set expansion is a separate authority-only recording decision under PR G1;
- fixed entries are reported as stale and must be removed by a safe contraction PR;
- privileged path permissions remain in `tools/web-change-policy.json` and are never copied into this store;
- persisted-format compatibility contracts remain governed by their own evidence rules and are never converted into generic ignored violations.

The store is a violation ledger, not a metric baseline. Repository-wide counts, LOC, module totals, and heuristic-name inventories must not be committed to it.

### 6.8 Future deterministic-rule onboarding

A new deterministic rule must keep executable detection and intended architecture in separate pull requests. The trusted-base checker evaluates a pull request's runtime using:

```text
base detector implementation
+ base architecture rules
+ base privileged policy
```

The proposed head rule registry is read only as inert JSON for schema and authority-delta validation. It never authorizes runtime in the same pull request. The proposed head detector module is never executed by the trusted-base workflow.

Authority changes have these directions:

| Change | Direction |
| --- | --- |
| add a rule after its detector exists | authority expansion/activation |
| remove a rule | relaxation, not contraction |
| change a rule's detector ID | reinterpretation; migrate to an already-merged versioned detector ID in a separate authority-only PR |
| change a rule's kind or subject encoding | reinterpretation; introduce a versioned detector, then reconcile the rule and any affected store identities in a dedicated authority-only migration |
| add an allowed definition path | expansion or preauthorization |
| add an allowed canonical edge | expansion or preauthorization |
| change `exactDefinitionCount` or `exactActiveEdgeCount` | unsupported in schema version 1; requires a later schema plan |
| `hard` to `frozen` or `report-only` | relaxation |
| `frozen` to `report-only` | relaxation |
| `baselineEligible: false` to `true` | relaxation |
| remove an allowed definition path | contraction, only when head source no longer requires it |
| remove an allowed canonical edge | contraction, only when the edge is inactive on head source |
| `report-only` to `frozen` or `hard` | tightening, only when the destination-mode gate passes |
| `frozen` to `hard` | tightening, only when no accepted violation remains |
| `baselineEligible: true` to `false` | tightening, only when no accepted violation remains |

Every new rule follows a detector-first path:

1. add a new versioned rule-specific detector in a checker-only PR, with characterization tests and no rule-authority change; do not extend executable semantics under an active ID;
2. add the rule to the inert authority registry in a later authority-only PR;
3. activate only the mode justified by untouched-base evidence.

The rules JSON references a stable detector ID. It must not contain regular expressions, JavaScript source, executable expressions, symbol patterns, or other data that can narrow detection. Rule-specific matching code remains in the trusted detector module.

An active detector ID is an immutable evidence contract. Its subject category, matching semantics, source scope, and canonical subject encoding must not change while any base rule references it. Detector IDs are versioned from their first introduction, for example:

```text
semantic-owner.render-request.v1
canonical-path.render-request.v1
```

The fixed detector corpus and current untouched `origin/dev` output characterize the contract but cannot prove equivalence on future syntax. During the initial rollout, any executable edit inside an active detector or a transitive shared detection helper uses a new ID even when the edit is believed to preserve behavior. Same-ID changes are limited to non-executable text or formatting that leaves executable tokens unchanged. The standard sequence is:

1. add `*.v2` beside `*.v1` in a checker-only PR and characterize both outputs;
2. migrate the inert rule from `*.v1` to `*.v2` in a later authority-only PR after untouched-base validation;
3. remove unreferenced `*.v1` in a later checker-only cleanup PR.

Adding v2 must not mutate a shared helper transitively used by active v1. Keep the v1 executable path unchanged and introduce a version-specific helper when necessary; temporary detector-internal duplication is removed with v1 in step 3. The new ID must not silently reuse the old subject grammar. If the subject grammar changes, step 2 is a dedicated authority-only migration that changes the rule reference and every affected accepted-violation identity atomically, with an exact old-to-new subject mapping validated against untouched base source. It must not add or drop debt. A detector-only PR is therefore a trusted checker-code review boundary, not an automatic proof of semantic immutability.

#### Clean-base rule

When the untouched target branch has no violations:

1. merge the detector-only PR;
2. add the rule in an authority-only PR as `hard`;
3. let the trusted-base checker validate the candidate inert JSON against the already-merged detector and untouched base;
4. do not create a violation-store entry.

#### Existing-debt rule

When the untouched target branch has exact existing violations and Phase 0 selects the rule as materially valuable:

1. **Detect:** merge the rule-specific detector in a checker-only PR.
2. **Observe:** add the rule in an authority-only PR with `report-only` mode and stable violation identities.
3. **Enable mechanics:** merge conditional PR F2 without changing rule or exception authority.
4. **Record:** after the report-only rule has produced untouched-base evidence, add only the observed exact base violations to `tools/web-architecture-violations.json` in a separate authority-only PR G1.
5. **Enforce:** in later authority-only PR G2, change only the rule's mode from `report-only` to `frozen`.
6. **Contract:** when runtime work fixes an accepted violation, remove the corresponding store entry in the same PR through the narrow safe-contraction path.

The registry schema recognizes:

```text
hard         cannot be baselined
frozen       exact accepted base violations allowed; new violations fail
report-only  produces evidence but does not authorize merge failure by itself
```

Schema recognition does not activate an unavailable mode. Until conditional PR F2 has merged, candidate authority that selects `frozen` fails closed and no accepted-violation store is permitted.

A checker implementation and authority data must never be introduced, activated, expanded, relaxed, or reinterpreted in the same PR. Initial rollout follows the staged sequence in Section 10; it has no bootstrap exception that combines detector and authority.

## 7. Enforcement tiers and machine-readable reflexion output

The repository must not claim that static analysis computes complete `OE`, `PE`, or `CB`. Enforcement is divided into four tiers, followed by a common reflexion-report contract.

### 7.1 Hard invariants

These fail immediately and are not baseline-able:

- unauthorized privileged operator or importer;
- a first-party Web import cycle;
- an unapproved production dependency;
- runtime/guard self-authorization;
- a missing required canonical owner or entry point;
- a new registered semantic-authority location when the rule is declared clean;
- a registered single canonical entry edge that is absent, unapproved, or active in parallel with another candidate edge;
- an unapproved performance or scientific-output regression where a gate exists.

### 7.2 Frozen deterministic violations

These are exact violations of valuable deterministic rules that existed when the rule was introduced and could not be fixed immediately.

This tier is active only when conditional PR F2 has merged. A clean initial rollout does not implement or silently accept frozen violations.

The checker must classify each as:

- accepted existing;
- newly introduced;
- fixed/stale.

Only exact entries in `tools/web-architecture-violations.json`, when the file exists, are accepted. New violations fail. A fixed entry is a failing stale-authority condition unless the same PR removes exactly that entry through the safe-contraction path. When the final accepted violation is removed, the store file must be deleted rather than retained empty.

### 7.3 Differential metrics and reports

For measurable inventories, CI reports:

```text
before
added
removed
after
delta
```

This applies to:

- registered semantic-authority locations;
- registered canonical-path evidence;
- privileged operator permission entries;
- privileged importer permission entries;
- production modules;
- exports;
- `create*` declarations;
- reactive declarations;
- watcher calls;
- compatibility-like names;
- resource-like names;
- session object keys;
- first-party static import cycles;
- production files and churn.

Only metrics with a separately defined decision rule may fail the build.

### 7.4 Trend and diagnostic metrics

The following remain non-blocking in the initial rollout:

- module count and LOC;
- consumer privileged fan-out;
- conventional coupling metrics such as `Ca`, `Ce`, CCD, or NCCD;
- temporal change coupling;
- change amplification;
- rework and revert rate;
- performance-smoke trends.

These metrics can identify investigation targets but do not independently prove semantic ownership quality.

### 7.5 Reflexion report format

Do not encode source state, enforcement mode, baseline relation, authority-delta handling, and decision in one expanding enum. For each registered deterministic rule and relevant subject, the checker emits five report fields:

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

The first four fields are independent inputs. `decision` is derived from them by one pure function; it is not independently assigned or stored. `baselineRelation: FIXED` means the accepted subject is no longer a source violation. For a frozen subject, `authorityResolution` states what the head accepted-violation store did with the relevant base entry; it is not a source observation or the general rule-authority delta classifier. This allows an exact same-PR contraction to pass while retained stale authority still fails.

PR E establishes this five-field output envelope and implements the `HARD` and `REPORT_ONLY` rows; their baseline and authority fields are always `NOT_APPLICABLE`. Conditional PR F2 adds the `FROZEN` rows, accepted-store classification, and non-`NOT_APPLICABLE` authority resolutions. If F2 is skipped, passing `FROZEN` to evaluation is an unavailable-mode error that fails closed.

The required decision table is:

| Mode and observation | Baseline relation | Authority resolution | Decision |
| --- | --- | --- | --- |
| `HARD` + `CONFORMING` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `HARD` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `FAIL` |
| `REPORT_ONLY` + `CONFORMING` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `REPORT_ONLY` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `REPORT` |
| `FROZEN` + `CONFORMING`, with no fixed accepted subject | `NOT_APPLICABLE` | `NOT_APPLICABLE` | `PASS` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `ACCEPTED` | `RETAINED` | `PASS` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` after its accepted entry was removed | `ACCEPTED` | `INVALID_CHANGE` | `FAIL` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `NEW` | `NOT_APPLICABLE` | `FAIL` |
| `FROZEN` + `CONFORMING` for a formerly accepted subject | `FIXED` | `RETAINED` | `FAIL` |
| `FROZEN` + `CONFORMING` for a formerly accepted subject | `FIXED` | `EXACT_CONTRACTION` | `PASS` |
| `FROZEN`, with an invalid accepted-store expansion, removal, or mismatch | As derived | `INVALID_CHANGE` | `FAIL` |

Any combination outside this table is an internal checker error and fails closed. `EXACT_CONTRACTION` is valid only when the whole proposed store delta satisfies the set equality in Section 6.6. Hard and report-only rules must not consult the accepted-violation store and therefore always use `authorityResolution: NOT_APPLICABLE`. Rule-registry expansions, relaxations, tightenings, and reinterpretations remain separate candidate-authority decisions under Section 6.8.

The output must identify:

- rule key and stable subject when applicable;
- intended owner/path/relation from the base rule authority;
- observed source evidence derived by the base detector;
- the five fields above;
- accepted violation identity, if applicable.

CI will continue to report, but not automatically reject solely on the basis of:

- a new module;
- a new export;
- a new watcher;
- a new reactive declaration;
- names containing `legacy`, `fallback`, `manager`, `cache`, `token`, `handle`, `journal`, `protocol`, or `normalize`;
- a compatibility-like branch that requires human interpretation.

The pull-request declaration supplies the evidence for unregistered semantic ownership and compatibility burden. The recorded maintainer architecture review owns the judgment; a filled template alone is not approval.

## 8. Target authority layout

| Concern | Permanent owner |
| --- | --- |
| Ratchet definitions, formulas, exceptions, and acceptance rule | `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` |
| Repository-wide architecture and cross-surface convergence | `CLAUDE.md` |
| Current Web runtime path and module ownership | `gbdraw/web/CLAUDE.md` |
| Executable source-evidence detection only | `tools/web-architecture-detectors.mjs` |
| Pure schema, differential-set, graph, result, and conditional accepted-set evaluation | `tools/web-architecture-evaluation.mjs` |
| Intended semantic-owner/canonical-path rules and enforcement modes | `tools/web-architecture-rules.json` |
| Exact accepted deterministic violations, only after conditional F2 and G1 | `tools/web-architecture-violations.json` |
| Permitted privileged operator and importer paths | `tools/web-change-policy.json` |
| Sole Web change CLI/CI entry point, Git/file I/O, orchestration, and reporting | `tools/check-web-change-budget.mjs` |
| Cross-cutting architecture and workflow contracts | `tests/web/architecture-contracts.test.mjs` |
| Fixture-heavy ratchet schema, algebra, graph, and decision tests | `tests/web/architecture-ratchet-fixtures.test.mjs` |
| Policy workflow and guard-separation explanation | `docs/internal/WEB_CHANGE_POLICY.md` |
| Human change declaration and evidence | `.github/pull_request_template.md` |
| Human architecture approval | Explicit recorded manual maintainer review |
| Contributor entry point | `CONTRIBUTING.md` |
| Agent entry point | `AGENTS.md` |

The same formulas and definitions must not be copied into every entry-point document. Entry points contain a short rule and a link to the normative owner.

## 9. Non-goals

This initiative must not:

- modify production files under `gbdraw/web/js/`, `gbdraw/web/vendor/`, or `gbdraw/web/index.html`;
- repair current runtime ownership problems;
- rename `allowedPrivilegedOwners` in the policy schema;
- create a universal catalog of every feature in gbdraw;
- infer semantic ownership from identifier names alone;
- reject every module, export, watcher, or reactive declaration;
- create a second Web checker;
- equate one CLI/CI entry point with one monolithic implementation module;
- combine executable detector logic with intended owner paths, edges, counts, enforcement modes, baseline eligibility, or accepted violations;
- put filesystem, Git, environment, reporting, or intended-authority ownership in the pure evaluator;
- put regular expressions, JavaScript source, executable expressions, or detector-narrowing symbol patterns in the inert rule registry;
- turn the initial registry into a generic architecture-rule DSL or register more than three initial rules;
- change an active detector's evidence semantics without introducing a new versioned detector ID;
- implement accepted-store parsing/evaluation or frozen-rule mechanics in the initial rollout without a materially valuable existing-debt rule selected by Phase 0;
- execute a proposed head detector module from the trusted-base workflow;
- add an npm or Python dependency;
- introduce a build step;
- create a committed repository-wide metric snapshot that will immediately become stale;
- use the exact accepted-violation store for totals, heuristic signals, or broad ignores;
- add a weighted maturity or conformance score;
- install ArchSteer, ArchUnit, dependency-cruiser, CodeScene, SonarQube, NDepend, or another external governance dependency;
- mix lagging historical metrics into the initial PR merge gate;
- claim that the source model is a complete automatic reconstruction of the architecture;
- present registered file-level observations as repository-wide `OE`, `PE`, or `CB`;
- broaden the task into unrelated branching, documentation, or runtime cleanup.

If implementation appears to require a production runtime edit, stop and create a separate implementation plan for that runtime change.

## 10. Required implementation sequence

Use separate pull requests. Every work branch starts from the latest merged `origin/dev`, every implementation pull request targets `dev`, and `main` remains the release branch promoted from `dev`. Do not infer the pull-request target merely from the work-branch base.

The sequence deliberately protects future authority paths before creating them, introduces executable detection before inert authority, and activates authority only after the trusted-base checker can validate it. There is no bootstrap exception that combines checker implementation and authority.

Whenever a phase edits prose intended for contributors or maintainers, including `CONTRIBUTING.md`, the PR template, repository guidance, the normative document, or `WEB_CHANGE_POLICY.md`, apply `.agents/skills/avoid-ai-writing/SKILL.md` and record the final verification. This internal implementation plan itself remains exempt under repository guidance.

### Prerequisite PR 0A — Align the integration branch and CI coverage

#### Purpose

Make `dev` a real protected integration target before any ratchet PR depends on pull-request checks.

#### Modify

- `.github/workflows/test.yml`:
  - run `pull_request` checks for both `dev` and `main`;
  - preserve the existing `dev` push promotion-range behavior;
  - preserve existing job names used by branch protection.
- `.github/workflows/web-base-policy.yml`:
  - run `pull_request_target` for both `dev` and `main`;
  - continue to check out only the trusted base and fetch the head as Git data;
  - never execute head code.
- `CONTRIBUTING.md`:
  - distinguish work-branch base, PR target, integration branch, and release branch;
  - state that agent work branches start from `origin/dev` and target `dev`;
  - state how maintainers promote `dev` to `main`.
- `tests/web/architecture-contracts.test.mjs`:
  - assert both workflows cover `dev` and `main` PRs;
  - retain the trusted-base safety assertions.

#### Bootstrap condition

This PR cannot be protected by the `dev` pull-request trigger it is introducing. Record that limitation explicitly. Run the exact checks locally, require maintainer review of the workflow diff, and after merge confirm that the next trivial guard-only PR receives both named checks on `dev`. Do not weaken `pull_request_target` permissions or check out head code to avoid this one-time condition.

#### Verification

```bash
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

#### Acceptance criteria

- `dev` and `main` PRs trigger both the normal and trusted-base workflows;
- work branch, PR target, integration branch, and release branch are documented separately;
- existing `dev` promotion-range checks remain intact;
- no checker, runtime, policy allowlist, or generated artifact changes;
- the exact required branch-protection action for `dev` is completed or recorded as an unresolved external gate.

#### Proposed commit title

```text
ci: cover dev pull requests with trusted Web checks
```

### Prerequisite PR 0B — Pre-register future guard and authority paths

#### Purpose

Teach the already-trusted base checker how to classify every file the ratchet will later create. The files remain absent in this PR.

#### Modify

##### `tools/check-web-change-budget.mjs`

Add the following future paths to `guardPaths`:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
.github/pull_request_template.md
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
tools/web-architecture-rules.json
tools/web-architecture-violations.json
tests/web/architecture-ratchet-fixtures.test.mjs
```

Add only executable implementation paths to `checkerImplementationPaths`:

```text
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
```

Add only inert authority paths to `authorityPaths`:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
tools/web-architecture-rules.json
tools/web-architecture-violations.json
```

Do not classify `tools/web-architecture-rules.json` as checker implementation. Do not classify either executable mechanics module as authority. The future PR template and fixture test are guard paths but are neither checker implementation nor intended machine authority.

##### `tests/web/architecture-contracts.test.mjs`

Add fixture coverage proving:

- runtime cannot co-change with any future guard path;
- checker implementation cannot co-change with any future authority path;
- the future files need not exist for the checker to pass;
- adding an unrelated unregistered path does not acquire authority accidentally;
- the evaluator is classified with checker implementation, while its fixture test is guard-only and the classifications remain disjoint;
- the future PR template is guard-only, so it cannot be weakened in a runtime PR.

#### Verification

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

Run the checker against explicit base/head refs before handoff.

#### Acceptance criteria

- every future detector, evaluator, fixture-test, PR-template, and authority path is protected before it exists;
- no future file is added yet;
- checker implementation and authority classifications are disjoint;
- no workflow, runtime, policy allowlist, or dependency changes.

#### Proposed commit title

```text
ci: reserve architecture ratchet guard paths
```

### Phase 0 — Baseline and reconciliation

#### Purpose

Establish an evidence-backed baseline after PR 0A and PR 0B have merged and before adding authority or detector files.

#### Actions

1. Fetch `origin` and record the exact `origin/dev` commit SHA.
2. Create a temporary work branch from `origin/dev` targeting `dev`.
3. Read:
   - `AGENTS.md`;
   - `CLAUDE.md`;
   - `gbdraw/web/CLAUDE.md`;
   - `CONTRIBUTING.md`;
   - `docs/internal/WEB_CHANGE_POLICY.md`;
   - `tools/web-change-policy.json`;
   - `tools/check-web-change-budget.mjs`;
   - `tools/web-change-source.mjs`;
   - `tests/web/architecture-contracts.test.mjs`;
   - `.github/workflows/test.yml`;
   - `.github/workflows/web-base-policy.yml`.
4. Run without editing:
   ```bash
   node tools/check-web-change-budget.mjs
   node --test tests/web/architecture-contracts.test.mjs
   node --test tests/web/*.test.mjs
   ```
5. Create untracked evidence under `/tmp/gbdraw-architecture-ratchet-baseline/`:
   - base commit;
   - current capability and policy keys;
   - current privileged operator permission entries by capability;
   - current privileged importer permission entries by target;
   - current detector patterns and their exact source matches;
   - current first-party import-cycle count, including self-edge inspection;
   - current report-only inventory totals;
   - proposed rule keys, detector IDs, source evidence, and intended constraints;
   - the proposed discriminated rule kind and detector output subject category for each rule;
   - whether each proposed rule is clean-base, existing-debt, or not deterministically enforceable;
   - the fixed detector fixture corpus and untouched-base output used as characterization evidence, without treating them as proof of semantic equivalence;
   - whether at least one materially valuable existing-debt rule is selected and therefore conditional PR F2 is required;
   - exact test commands and results.
6. Reconcile definitions duplicated between the checker and architecture tests before extraction.
7. For each human-reviewed changed capability example, record concrete before/after owner, canonical-path, and compatibility-path sets rather than totals alone.
8. Verify both named Web checks run for a `dev` pull request and verify or record the external branch-protection requirement.
9. Do not create a rule registry or accepted-violation store in this phase.

#### Gate

Phase 0 passes only when:

- current checker and architecture contracts pass on the recorded base;
- untouched `origin/dev` has zero first-party Web static import cycles, including self-imports;
- the duplicated detection definitions are fully inventoried;
- every proposed detector has a versioned stable ID and produces deterministic source evidence with a fixed subject category and encoder;
- every proposed rule is classified as clean-base, report-only existing-debt, or not deterministically enforceable;
- the initial candidate inventory contains no more than three rules, normally one semantic-owner rule and one canonical-path rule, with any third rule justified by Phase 0 evidence;
- proposed detector code contains no intended paths, counts, modes, eligibility, or exceptions;
- proposed rules use only the required kind-specific schema, contain no executable expressions or detector-narrowing patterns, and do not require a generic DSL;
- the F2 trigger is recorded as `required` with the selected existing-debt rule keys or `skipped` with evidence that no such initial rule exists;
- `dev` PR CI coverage is observed, and missing branch-protection authority remains explicitly unresolved rather than claimed complete.

If any first-party static import cycle exists, do not baseline or freeze it. Stop the cycle portion of this initiative, create a separate runtime-cleanup plan, and revise this plan before continuing to PR F1. The other baseline findings do not make a nonzero-cycle base eligible for the hard cycle gate.

### PR A — Establish the human authority and PR declaration

Create a fresh branch from updated `origin/dev` and target `dev`.

#### Add `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`

The normative document must contain:

- scope, authority, definitions, precedents, and the intended/source/reflexion distinction;
- the `OE`, `PE`, and `CB` review model and non-scalar vector;
- the two-layer `MERGEABLE` rule;
- the prohibition on claiming complete machine-calculated `OE`, `PE`, or `CB`;
- concrete before/after set evidence requirements and scope-completeness rationale;
- semantic owner versus privileged operator;
- detector implementation versus inert rule authority versus accepted exceptions;
- sole CLI orchestration versus pure detector and evaluation mechanics;
- trusted-base evaluation order and prohibition on executing head detectors;
- detector code as a reviewed trust root, characterization evidence as non-proof, and the detector-add/rule-migrate/detector-remove sequence for executable changes;
- the two initial discriminated rule kinds and the maximum-three-rule rollout scope;
- exact single-canonical-entry-edge semantics, including preauthorization without parallel active paths;
- hard, frozen, and report-only modes with the five report fields and pure decision table;
- the Phase 0 trigger that makes frozen and accepted-store mechanics conditional;
- accepted-violation protocol and prohibition on CI refreezing;
- ownership, canonical-path, compatibility, and safe-contraction examples using explicit sets rather than unsupported totals;
- the detector-first rule-onboarding and authority-direction procedures.

Do not duplicate the current module ownership table from `gbdraw/web/CLAUDE.md`.

#### Add `.github/pull_request_template.md`

Require exactly one author declaration:

```markdown
## Purpose

## Verification

## Architecture impact

Select exactly one:

- [ ] This is not architecture-bearing. Rationale:
- [ ] This is architecture-bearing; the evidence below is complete.

For an architecture-bearing change:

- Capability or behavior:
- Why this changed-capability scope is complete:
- Semantic owners before:
- Semantic owners after:
- Canonical production paths before:
- Canonical production paths after:
- Compatibility paths before, with stable IDs:
- Compatibility paths after, with stable IDs:
- Superseded owner/path removed:
- New module classification:
- Persisted-compatibility evidence and removal condition:
- Performance/scientific-output evidence, when material:
- Deterministic checker evidence:
- Maintainer architecture decision comment permalink, required before merge:

- [ ] No second parser, normalizer, state mirror, lifecycle owner, or canonical pipeline was added.
- [ ] Every superseded implementation was removed in this change.
- [ ] No accepted deterministic violation was added in this implementation PR.
```

Do not ask for unsupported repository-wide `OE`, `PE`, or `CB` numbers. A completed author declaration is evidence, not architecture approval.

The architecture owner must review the exact proposed head and post a structured pull-request comment in this form:

```markdown
Architecture decision: APPROVED
Reviewed head SHA: <full commit SHA>
Reviewed capability scope:
- <capability and stable before/after set references>
Decision:
- delta(OE) <= 0
- delta(PE) <= 0
- delta(CB) <= 0, or: <persisted-compatibility exception and evidence>
- superseded owners and paths removed: yes
Limitations:
- <none, or explicit non-authorizing limitation>
```

The pull-request template stores the permalink to that comment. The decision is valid only while `Reviewed head SHA` exactly equals the PR head SHA; any later commit invalidates it and requires a new owner comment. An automated agent must not post the approval, and the author declaration cannot substitute for it. A sole architecture maintainer who is also the PR author may post the separate owner comment after reviewing the final head.

This remains a manual architecture-review gate. Workflow presence and branch protection do not by themselves prove that the structured comment exists, was posted by the architecture owner, or matches the final head. The normative document must record this limitation rather than implying automatic enforcement.

#### Modify `tests/web/architecture-contracts.test.mjs`

Require the created PR template to retain these exact process anchors without freezing all prose:

```text
Architecture impact
This is not architecture-bearing
This is architecture-bearing
Semantic owners before
Semantic owners after
Canonical production paths before
Canonical production paths after
Maintainer architecture decision comment permalink
```

The test must fail if the template is deleted or any anchor disappears. It does not inspect GitHub comments, comment authors, or live head SHAs.

#### Modify with short references only

- `AGENTS.md`: non-negotiable ratchet rule and normative link.
- `CLAUDE.md`: short repository architecture-ratchet reference.
- `gbdraw/web/CLAUDE.md`: distinguish current ownership documentation from the rule for changing ownership.
- `CONTRIBUTING.md`: require the architecture declaration and recorded maintainer decision; preserve the branch model established in PR 0A.

Because these files are human-facing, apply `.agents/skills/avoid-ai-writing/SKILL.md` before editing and verify the final prose.

#### Verification and acceptance

```bash
git diff --check
node --test tests/web/architecture-contracts.test.mjs
python -m pytest tests/ -k documentation -m "not slow"
```

- one normative definition owner exists;
- all entry points link rather than copy the model;
- the template collects concrete sets and separates author evidence from an exact-head, structured maintainer approval;
- the template is a protected guard path and its minimum process anchors are executable contracts;
- the manual nature of the maintainer-comment gate is explicit;
- the base checker treats the new normative document as authority;
- no runtime, checker implementation, policy allowlist, or generated artifact changes.

#### Proposed commit title

```text
docs: define the architecture fitness-function ratchet
```

### PR B — Extract pure source-evidence detectors

#### Add `tools/web-architecture-detectors.mjs`

Create one dependency-free ES module containing executable detection only. It may contain private generic scan helpers and stable rule-specific detector functions, for example:

```javascript
export const WEB_ARCHITECTURE_DETECTORS = Object.freeze({
  'semantic-owner.render-request.v1': detectRenderRequestAuthority,
  'canonical-path.render-request.v1': detectCanonicalRenderRequestPath
});
```

It may contain masked-source regexes and normalized import-edge extraction. Cycle classification belongs to the pure evaluator when PR F1 adds it. The detector module must not contain:

- intended definition paths or allowed canonical edges;
- exact or maximum counts;
- enforcement modes;
- baseline eligibility;
- accepted violations;
- compatibility authorization;
- privileged permission arrays.

Avoid stateful global `.test()` behavior. Use Node built-ins only, preferably none beyond the shared source parser.

#### Modify checker and tests

- Make the checker and architecture tests consume the same detector implementation.
- Keep existing inline intended constraints active temporarily; do not add the new rules JSON yet.
- Characterize exact current capability keys, import targets, operator matches, and source masking before removing duplicated detector code.
- Treat each exported detector ID as the immutable evidence contract from Section 6.8. Add a fixed fixture corpus covering normalized facts, subject category, and canonical subject encoding.
- Record fixed-corpus and untouched-source output equivalence as characterization evidence for the initial extraction, not as proof for arbitrary future code.
- State that any later executable edit to an active detector or its transitive shared helpers requires a new versioned ID. Same-ID changes must leave executable tokens unchanged.
- Require a new detector version to preserve the old version's transitive executable path, using a temporary version-specific helper when shared logic must change.
- Retain independent expected-result assertions in tests; sharing detector implementation must not turn every test into a tautology.
- Remove only superseded detector implementations in this PR.

The transitional inline intended constraints remain the sole active authority until PR E. The future rules file is absent, so there are not two active rule authorities.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

- detection has one executable owner;
- detector IDs are versioned, stable, and rule-specific; executable evolution uses a new ID, while output-equivalence fixtures remain characterization evidence;
- the handoff identifies detector code as a trusted checker-code review boundary rather than claiming mechanical semantic proof;
- current source matches and policy coverage are byte-for-byte or set-for-set equivalent to the base characterization;
- no intended paths, thresholds, modes, or exceptions moved into detector code;
- no runtime, inert authority, policy allowlist, workflow, or dependency changes.

#### Proposed commit title

```text
refactor: centralize Web architecture detectors
```

### PR C — Add inert rule-schema and trusted candidate validation

Extend the existing checker through a pure evaluation module; do not add a second CLI/CI entry point and do not add the rules file yet.

#### Add `tools/web-architecture-evaluation.mjs`

Create a dependency-free pure module for:

- strict rule-schema validation;
- detector-output/rule-kind compatibility checks;
- authority-delta classification;
- the five-field output envelope and decision mechanics added by later phases;
- sorted set operations and, when added in PR F1, strongly connected component evaluation;
- accepted-store validation and set algebra only when conditional PR F2 is required.

The module receives already-loaded plain data and normalized source facts and returns plain deterministic results. It must not read files, run Git, inspect the environment, write output, own intended paths or thresholds, or import the rule and violation JSON files directly. `tools/check-web-change-budget.mjs` remains the sole CLI/CI entry point and owns Git/file I/O, trusted-base orchestration, error presentation, and summary reporting. `tools/web-architecture-detectors.mjs` remains the only owner of registered source-fact extraction.

PR C implements only the rule schema, detector-kind compatibility, stable collection, and authority-delta primitives needed for candidate validation. PR E adds `HARD` and `REPORT_ONLY` decision rows, PR F1 adds differential set and SCC mechanics, and conditional PR F2 adds accepted-store and `FROZEN` mechanics. This assigns one final owner without pulling later enforcement into PR C.

#### Rule schema support

Support an optional strict JSON file at:

```text
tools/web-architecture-rules.json
```

Absence means no new registry rules during the staged rollout. Schema version 1 must define a `rules` array stored in ascending rule-key order. Kind-specific path and edge arrays must also use canonical sorted order. Version 1 supports only a discriminated union of the kinds needed for the initial rollout:

```json
{
  "schemaVersion": 1,
  "rules": [
    {
      "key": "semantic-owner.render-request",
      "kind": "single-semantic-owner",
      "detector": "semantic-owner.render-request.v1",
      "allowedDefinitionPaths": ["services/session-request.js"],
      "exactDefinitionCount": 1,
      "enforcement": "hard",
      "baselineEligible": false
    },
    {
      "key": "canonical-path.render-request",
      "kind": "single-canonical-entry-edge",
      "detector": "canonical-path.render-request.v1",
      "allowedEdges": [
        {
          "from": "app/run-analysis.js",
          "to": "services/session-request.js"
        }
      ],
      "exactActiveEdgeCount": 1,
      "enforcement": "hard",
      "baselineEligible": false
    }
  ]
}
```

`single-semantic-owner` accepts only its common fields plus `allowedDefinitionPaths` and `exactDefinitionCount`, whose only valid version 1 value is `1`. This permits preauthorizing a location move without permitting two simultaneous owners.

`single-canonical-entry-edge` accepts only its common fields plus `allowedEdges` and `exactActiveEdgeCount`, whose only valid version 1 value is `1`. The detector emits the complete set `observedEdges` within its immutable rule-specific scope. Conformance is:

```text
observedEdges is a subset of allowedEdges
AND size(observedEdges) = exactActiveEdgeCount
```

An empty observed set is `ABSENT_REQUIRED`. Any unapproved edge is `DIVERGENT`. Two active candidate edges are `DIVERGENT` even when both were preauthorized. An allowed but inactive edge is not missing and does not fail. This permits an authority-only PR to add a future edge, a later runtime PR to switch atomically from the old edge to the new edge, and a later authority-only PR to remove the inactive old edge without ever permitting parallel canonical paths.

A detector declares its normalized output subject category, and the evaluator rejects a detector/rule pairing whose category does not match the rule kind. Version 1 must not provide a generic `allowedEvidencePaths`, `operation`, pattern, predicate, or catch-all kind.

The strict schema must:

- reject duplicate rule keys, duplicate paths, and duplicate edges;
- reject unknown detector IDs and unsupported kinds or modes;
- reject hard rules marked baseline-eligible and frozen rules marked ineligible;
- reject paths outside the normalized `gbdraw/web/js/` source scope;
- reject unknown fields;
- reject regex source, JavaScript source, executable expressions, wildcard paths, and detector-narrowing symbol patterns;
- require every referenced detector to provide the matching subject category and a stable subject encoder for its violation-producing rule kind;
- require nonempty `allowedEdges`, canonical `(from, to)` order, and `exactActiveEdgeCount: 1` for every `single-canonical-entry-edge` rule;
- reject a registry containing more than three rules during the initial rollout.

The three-rule cap is part of schema version 1. Raising the cap or adding a kind requires a later evidence-backed plan and schema-version migration; it is not an ordinary rule-authority edit.

#### Trusted-base evaluation order

After PR E activates the registry, runtime conformance always uses base detector + base rules + base privileged policy. Parse proposed head rules only as inert JSON and classify their authority delta. Never import or execute a proposed head detector module.

During PR C and PR D, the current inline intended constraints substitute for base rules and remain the sole active runtime authority. Candidate rules are evaluated against untouched base source for admission evidence but do not participate in runtime authorization until PR E removes the inline constraints and activates the unchanged base JSON. This staged inactive file is not a second active rule owner.

Candidate validation must reject:

- rules referring to detectors absent from the trusted base;
- an authority expansion combined with runtime, checker, detector, test-rule, policy, or workflow changes;
- a relaxation combined with any implementation change;
- a rule whose claimed clean base is actually divergent;
- a candidate `frozen` activation while conditional PR F2 is absent from the trusted base;
- a transition outside the authority-direction table in Section 6.8.

An authority contraction or tightening must still conform on head source under the trusted base detector. Do not let a head rule change authorize the same head runtime.

#### Tests and documentation

- Add `tests/web/architecture-ratchet-fixtures.test.mjs` for fixture-heavy schema, authority-delta, pure-evaluation, and malicious-input cases. Keep cross-cutting workflow and current-architecture contracts in `tests/web/architecture-contracts.test.mjs`.
- Add fixture tests for malformed and executable-looking JSON, wrong kind-specific fields, duplicate or unsorted edges, an active-edge count other than one, detector subject-category mismatch, unknown detectors, unavailable frozen activation, more than three rules, direction classification, and self-authorization attempts.
- Prove the trusted checker ignores a malicious head detector.
- Prove the pure evaluator produces no file, Git, environment, or output side effects and that the checker delegates evaluation without creating a second entry point.
- Update `docs/internal/WEB_CHANGE_POLICY.md` to explain detector/rule/policy separation and candidate inert-data handling.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

- candidate rules can be validated without executing head code;
- rule kinds are closed and kind-specific rather than a generic evidence-path schema;
- the canonical-edge schema distinguishes allowed preauthorization from the exactly one active edge required for conformance;
- pure evaluation is isolated from I/O, authority ownership, and source-fact extraction;
- missing rules file remains valid during rollout;
- no rule authority file is added;
- no runtime, policy allowlist, workflow, or dependency changes.

#### Proposed commit title

```text
ci: validate inert Web architecture rules
```

### PR D — Add the initial intended-architecture rules

Add `tools/web-architecture-rules.json` in an authority-only PR.

Use Phase 0 evidence to register only critical rules with deterministic detectors already merged in PR B:

- register no more than three rules; normally register one `single-semantic-owner` rule and one `single-canonical-entry-edge` rule;
- use a third rule only when Phase 0 demonstrates material value and deterministic evidence; it must be another instance of one of those two kinds rather than a reason to widen the schema;
- keep the dependency-cycle hard gate separate; it does not consume a registry slot;
- use `hard` for a clean-base invariant;
- use `report-only` for a valuable rule with exact existing debt;
- include an existing-debt rule only when Phase 0 marks conditional PR F2 as required;
- omit a rule that is not deterministically enforceable;
- do not use `frozen` yet when accepted violations would be required;
- do not add `tools/web-architecture-violations.json`.

The trusted-base checker from PR C validates the candidate JSON against the untouched base detector. This PR must not change checker code, detectors, tests, runtime, policy, workflows, or documentation.

After merge and before PR E, the file is protected inert intended authority awaiting activation; the current inline constraints remain the sole active runtime gate.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

- every rule references an already-merged detector;
- the registry has at most three rules, uses only the two version 1 kinds, and explains any third rule in Phase 0 evidence;
- every kind-specific definition path, allowed edge, and exact active count matches recorded intended architecture;
- every canonical rule has exactly one observed edge on untouched base, all observed edges are allowed, and inactive preauthorized edges are not treated as missing;
- every hard rule passes on untouched `origin/dev`;
- every report-only rule exposes exact stable subjects without failing;
- the PR contains inert authority data only.

#### Proposed commit title

```text
chore: declare initial Web architecture rules
```

### PR E — Activate the base rule authority and normalized result model

After PR D merges, update the checker, pure evaluator, and tests without changing the rules JSON.

#### Checker behavior

- Consume the unchanged base `tools/web-architecture-rules.json` as intended architecture.
- Remove the superseded inline intended constraints in the same PR.
- Keep `tools/web-change-policy.json` as the sole privileged path-permission owner.
- Implement the five-field output envelope and the `HARD` and `REPORT_ONLY` rows of the pure decision function from Section 7.5 in `tools/web-architecture-evaluation.mjs`.
- Emit `baselineRelation: NOT_APPLICABLE` and `authorityResolution: NOT_APPLICABLE` for both supported modes.
- Fail closed when a rule requests `FROZEN` before conditional PR F2 has supplied its mechanics.
- Fail closed on impossible field combinations.
- Report rule and subject rows in stable rule-key/subject order.
- Count file-level registered authority locations without calling them complete `OE`.
- Continue selected canonical-path enforcement through explicit registered evidence, not whole-program path inference.

#### Tests

Put the table-driven decision and malformed-combination cases in `tests/web/architecture-ratchet-fixtures.test.mjs`; keep only end-to-end activation and cross-owner contracts in `tests/web/architecture-contracts.test.mjs`. Cover every allowed result combination and representative forbidden combinations. Prove:

- hard conforming passes and hard divergent/absent fails;
- report-only conforming passes and report-only divergent/absent reports without failing;
- no accepted store is loaded or consulted by hard or report-only modes;
- source observation, mode, baseline relation, authority resolution, and decision cannot contradict one another;
- a frozen input fails closed as unavailable rather than being interpreted as hard or report-only;
- a single canonical entry edge passes, zero edges fail as absent, an unapproved edge fails, two simultaneously active allowed edges fail, and an inactive preauthorized edge does not fail;
- removing the old inline authority does not change current policy coverage or rule results.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

- the rules JSON is the only active owner of registered intended paths, counts, and modes;
- the detector module is the only active owner of executable source matching;
- the evaluator contains only pure mechanics and derives every supported decision from the other four fields;
- frozen and accepted-store mechanics are absent unless the recorded Phase 0 trigger requires later PR F2;
- no temporary dual active authority remains;
- no runtime, rules authority, policy allowlist, workflow, or dependency changes.

#### Proposed commit title

```text
ci: enforce declared Web architecture rules
```

### PR F1 — Add differential reports and cycle enforcement

This PR is required. Extend the existing checker and pure evaluator without changing runtime, rules authority, privileged policy, or workflows. Keep Git/file I/O and reporting in the checker; keep stable set operations and SCC evaluation in the evaluator. Start only after Phase 0 has proved that untouched `origin/dev` has zero first-party static import cycles, including self-imports.

#### Differential reports

For each measurable inventory, report `before`, `added`, `removed`, `after`, and `delta`. At minimum include:

- Registered Authority Location Count;
- registered canonical-contract results;
- privileged operator permission entries by capability;
- privileged importer permission entries by target;
- production modules, exports, `create*` declarations, reactive declarations, watcher calls;
- compatibility-like and resource-like names and session object keys;
- first-party static import cycles;
- production files, additions, deletions, gross churn, and net additions.

Label heuristic categories `report-only`. Never label these inventories complete `OE`, `PE`, or `CB`. Add a compact summary table to terminal output and `GITHUB_STEP_SUMMARY`.

#### Dependency-cycle gate

Build the first-party static import graph for `gbdraw/web/js/**/*.js` with the existing parser. Ignore bare imports for cycle analysis while preserving their separate dependency rejection. Resolve the same relative forms as architecture tests.

A cycle exists when:

```text
SCC node count > 1
OR
SCC node count = 1 AND the node has an edge to itself
```

Fail on both cases. Report sorted SCC nodes and sorted internal edges. Dependency cycles are hard and cannot enter the accepted-violation store.

#### Other reports

- Report distinct privileged import-target fan-out per production module; keep it report-only initially.
- Report compatibility-like identifiers as before/after/delta; require the PR declaration to enumerate actual compatibility paths with stable IDs.
- Preserve persisted-format history and positive-fixture rules.

#### Required black-box fixtures

Prove at least:

1. current base architecture passes;
2. an extra registered authority location fails;
3. an unpreauthorized authority move fails;
4. self-import `a -> a` fails;
5. two-node and three-node cycles fail;
6. an acyclic graph passes;
7. multiple cycles report stable nodes and edges;
8. a compatibility-like name reports without failing;
9. a private helper with no authority passes;
10. ordinary and architecture budgets retain exact-limit behavior;
11. guard/authority separation rejects self-authorization;
12. zero, one, unapproved, and parallel active canonical entry edges behave as specified;
13. an inactive preauthorized canonical edge does not count as missing;
14. a proposed head detector or rule cannot authorize head runtime.

Execute the trusted-base checker binary from the base fixture revision when testing self-authorization. Do not accidentally import detector modules from the outer repository instead of the fixture revision.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
node tools/check-web-change-budget.mjs \
  --base "$(git merge-base origin/dev HEAD)" \
  --head HEAD
```

- deterministic inventories display both directions;
- every registered hard or report-only rule emits the five normalized report fields with non-applicable baseline and authority resolution;
- self-loops and multi-node cycles fail;
- report-only signals remain non-blocking;
- no accepted-store parsing/evaluation or frozen-rule implementation is added in this required PR;
- no runtime, rule, policy, workflow, or dependency changes.

#### Proposed commit title

```text
ci: add architecture deltas and cycle enforcement
```

### Conditional PR F2 — Add frozen-rule mechanics for selected existing debt

Run this PR only when Phase 0 selected at least one materially valuable existing-debt rule and PR D registered it as `report-only` with `baselineEligible: true`. Skip it when all initial rules are clean or intentionally remain report-only. Do not implement these mechanics for hypothetical future use.

Extend only `tools/check-web-change-budget.mjs`, `tools/web-architecture-evaluation.mjs`, and the two architecture test files. Do not change runtime, rule authority, the privileged policy, workflows, documentation, or create the violation store.

#### Frozen accepted-violation mechanics

- enable `FROZEN` as an available enforcement mode in trusted-base evaluation;
- support optional `tools/web-architecture-violations.json` with the strict schema from Section 6.7;
- treat absence as an empty accepted set and never write the store from CI;
- require rule keys and subjects to match base rules and stable rule-specific encoders;
- reject wildcard, duplicate, invented, line-number, malformed, and unknown-rule entries;
- implement the Section 6.6 accepted-set algebra in the pure evaluator;
- emit every `FROZEN` row from the Section 7.5 decision table;
- fail new frozen violations;
- make `FIXED + RETAINED` fail and `FIXED + EXACT_CONTRACTION` pass;
- require deletion of the store file when the final entry is contracted;
- permit store expansion only in later authority-only PR G1 for exact untouched-base signals from an already-merged eligible report-only rule;
- permit same-PR runtime plus store contraction only when removed entries exactly equal all fixed head subjects, the Section 6.6 equality holds, and no authority expands or other guard changes;
- never authorize privileged paths, compatibility contracts, cycles, dependencies, or trusted-base integrity violations through the store.

#### Required fixtures

Prove at least:

1. absent store behaves as an empty set;
2. frozen accepted, new, and fixed subjects produce every specified five-field combination and decision;
3. hard and report-only rules never load or consult the store;
4. invalid store entries fail;
5. exact store contraction passes;
6. retained fixed authority fails;
7. removal of a still-active accepted subject fails;
8. mixed contraction and expansion fails;
9. authority-only store recording accepts only exact untouched-base report-only signals;
10. a proposed head detector, rule, or store cannot authorize head runtime;
11. fixture execution uses the checker, detectors, and evaluator from the fixture's trusted base rather than the outer repository.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
node tools/check-web-change-budget.mjs \
  --base "$(git merge-base origin/dev HEAD)" \
  --head HEAD
```

- frozen activation is unavailable before this PR and available afterward;
- the full accepted/new/fixed and authority-resolution model is table-tested;
- no rule, violation-store authority, runtime, policy, workflow, documentation, or dependency change occurs;
- `tools/web-architecture-violations.json` remains absent until PR G1.

#### Proposed commit title

```text
ci: add frozen architecture violation support
```

### Conditional PR G1 and PR G2 — Freeze a selected rule with existing debt

Run these only after conditional PR F2. Skip both when F2 is skipped or every selected existing-debt rule remains intentionally report-only.

#### PR G1 — Record exact accepted violations

Add `tools/web-architecture-violations.json` in an authority-only PR. For every entry verify:

- the rule exists in base rules, is `report-only`, and is baseline-eligible;
- the detector and report-only rule were already merged;
- the exact subject exists on untouched `origin/dev`;
- the subject is stable across line movement and follows the rule-specific encoder;
- rationale, tracking evidence, and removal condition are present;
- it is not a privileged permission, compatibility contract, cycle, dependency, or heuristic signal.

Do not change checker, detector, rules, tests, runtime, policy, workflows, or documentation. The trusted-base checker validates the proposed store as inert data.

Proposed commit title:

```text
chore: record accepted Web architecture violations
```

#### PR G2 — Tighten the rule from report-only to frozen

After G1 merges, change only `tools/web-architecture-rules.json`. Change the selected rule's mode from `report-only` to `frozen`; do not edit checker code, tests, detectors, or the store. The generic frozen behavior was implemented in conditional PR F2.

The PR passes only when every current exact violation is accepted, every accepted entry is current, and every unlisted violation fails under the candidate-mode validation.

Proposed commit title:

```text
chore: activate frozen Web architecture rules
```

### Phase 7 — Repository-setting verification and rollout

Confirm that both protected branches use the expected checks:

- `dev`: `Web change budget`, `Web base policy (trusted base)`, and normal test requirements;
- `main`: the same trusted checks plus the repository's release/promotion requirements.

Do not infer branch protection from workflow presence. If the implementation agent cannot inspect or change settings, record the exact manual action and leave this phase incomplete.

For the first architecture-bearing runtime PR after required PR F1 and, when used, PR G2:

- inspect the terminal and GitHub summaries manually;
- confirm no false hard failure and no silent report-only failure;
- confirm the author enumerates concrete before/after sets and scope completeness;
- record the structured maintainer architecture decision comment and its permalink;
- verify the comment's full reviewed head SHA exactly matches the merge candidate; invalidate and repeat the decision after any new commit;
- confirm neither head rules nor head detectors authorize the same runtime.

### Phase 8 — Deferred trend metrics after sufficient history

Do not implement this phase in the initial sequence unless specifically authorized. After at least 20 relevant merged runtime PRs, consider median/P95 files and churn, follow-up fix or revert rate, temporal coupling, performance-smoke regression rate, and retry-dependent pass rate.

These remain lagging diagnostics. Use Git and existing GitHub metadata, add no production dependency, commit no manually maintained snapshot, create no scalar score, and keep historical analytics separate from the canonical checker unless a later plan establishes a deterministic merge requirement.

## 11. File-level change ledger

### New files

| File | PR | Purpose |
| --- | --- | --- |
| `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` | A | Normative rule |
| `.github/pull_request_template.md` | A | Human architecture declaration |
| `tools/web-architecture-detectors.mjs` | B | Executable source-evidence detection only |
| `tools/web-architecture-evaluation.mjs` | C | Pure schema, authority-delta, set, graph, and staged decision mechanics |
| `tests/web/architecture-ratchet-fixtures.test.mjs` | C | Fixture-heavy ratchet behavior and adversarial cases |
| `tools/web-architecture-rules.json` | D | Inert intended-architecture and enforcement authority |
| `tools/web-architecture-violations.json` | G1, conditional | Exact accepted deterministic violations; omit when empty |

### Modified files

| File | PR | Allowed change |
| --- | --- | --- |
| `AGENTS.md` | A | Short non-negotiable rule and link |
| `CLAUDE.md` | A | Short repository architecture reference |
| `gbdraw/web/CLAUDE.md` | A | Link current ownership table to ratchet |
| `CONTRIBUTING.md` | 0A, A | Branch model, PR declaration, and maintainer-decision requirement |
| `.github/workflows/test.yml` | 0A | Run normal PR checks for `dev` and `main` |
| `.github/workflows/web-base-policy.yml` | 0A | Run trusted-base PR checks for `dev` and `main` |
| `tools/check-web-change-budget.mjs` | 0B, B, C, E, F1, F2 conditional | Path protection, detector use, trusted orchestration, rule enforcement, I/O, and reports |
| `tools/web-architecture-evaluation.mjs` | E, F1, F2 conditional after creation in C | Pure supported-mode decisions, deltas, SCCs, and conditional accepted-set algebra |
| `tests/web/architecture-contracts.test.mjs` | 0A, 0B, A, B, C, E, F1, F2 conditional | Branch coverage, guard separation, PR-template anchors, detector characterization, and cross-cutting integration contracts |
| `tests/web/architecture-ratchet-fixtures.test.mjs` | E, F1, F2 conditional after creation in C | Decision-table, authority-delta, schema, graph, and conditional violation-store fixtures |
| `docs/internal/WEB_CHANGE_POLICY.md` | C | Detector/rule/policy terminology and trusted candidate handling |

### Files that must remain unchanged unless a separate plan is approved

```text
gbdraw/web/index.html
gbdraw/web/js/**
gbdraw/web/vendor/**
tools/web-change-policy.json
package.json
package-lock.json
pyproject.toml
```

The workflows may change only in PR 0A. A policy allowlist contraction discovered during implementation must be submitted separately under the existing contraction process. It is not an excuse to expand this plan.

## 12. PR checklist behavior

The PR template is evidence, not the source of truth.

Reviewers must reject an architecture-bearing PR when:

- the author marks it as non-architecture-bearing despite changing an owner, path, lifecycle, parser, persisted schema, or privileged surface;
- the declared changed-capability scope has no completeness rationale;
- concrete before/after owner, canonical-path, or compatibility-path sets are missing;
- unsupported totals are supplied in place of concrete sets;
- a new abstraction leaves the old path in place;
- a new parser duplicates browser, Python, library, or protocol semantics already owned elsewhere;
- a compatibility path lacks evidence;
- an accepted deterministic violation is added in an implementation PR or lacks a stable exact identity;
- a fixed accepted violation is retained instead of being contracted;
- the CI report shows a registered authority expansion;
- a deterministic check is bypassed or weakened in the same PR;
- executable code in an active detector or its transitive shared detection helpers changes under the same ID instead of using a versioned detector migration;
- a canonical entry rule permits zero, unapproved, or multiple simultaneous active edges;
- the protected PR template loses a required architecture-evidence anchor;
- the initial registry exceeds three rules or introduces a generic/catch-all rule kind;
- frozen or accepted-store mechanics are implemented without the recorded Phase 0 existing-debt trigger;
- a head detector, rule, or accepted-violation change attempts to authorize the same head runtime;
- no explicit maintainer architecture decision is recorded before merge;
- the decision comment is unstructured, lacks a permalink, or names a reviewed head SHA different from the merge candidate.

Reviewers must not reject a PR merely because:

- a focused helper module was added;
- line count increased;
- a canonical owner gained more importers;
- a large module was decomposed without changing semantics;
- a test module became larger to cover a real invariant;
- repository-wide `OE`, `PE`, or `CB` totals are absent when concrete changed-scope sets are present;
- an optional violation-store file is absent because the accepted set is empty.

## 13. Evidence and handoff requirements for every phase

Each phase handoff must include:

- base commit SHA;
- head commit SHA;
- changed-file list separated into production, guard, test, and documentation files;
- checker summary;
- exact test commands;
- pass/fail results;
- known limitations;
- proposed commit title and summary;
- confirmation that no unrelated files were modified;
- for every architecture-bearing PR, the structured approval-comment permalink and proof that its reviewed head SHA equals the final merge candidate.

Additional evidence by phase:

- PR 0A: observed `dev` and `main` workflow triggers, unchanged trusted-base checkout model, and branch-protection state;
- PR 0B: future path classification matrix, including PR-template, evaluator, and fixture-test paths, and self-authorization fixture results;
- Phase 0: exact `origin/dev` SHA, zero-cycle evidence including self-imports, selected rule classifications, and the explicit required/skipped decision for conditional PR F2;
- PR A: normative trust-boundary text, PR-template anchor results, and the documented manual limitation of maintainer-comment enforcement;
- PR B: detector before/after characterization, fixed-corpus/current-base output equivalence labeled as evidence rather than proof, versioned-ID inventory, policy-key coverage, and confirmation that detector code contains no authority;
- PR C: discriminated-schema fixtures, exact single-canonical-entry semantics, unavailable-frozen rejection, detector-kind compatibility, evaluator-purity evidence, authority-direction classification, and proof that head detector code is not executed;
- PR D: exact rule inventory, proof that the initial count is at most three, justification for any third rule, untouched-base result for every hard or report-only rule, and consistency with the F2 trigger;
- PR E: five-field envelope evidence for hard/report-only modes, unavailable-frozen failure, exact canonical-edge cases, and removal of superseded inline authority;
- PR F1: before/add/remove/after/delta report, self-loop and multi-node cycle results, and fixture-local trusted execution evidence;
- PR F2 when required: accepted/new/fixed and authority-resolution subjects, full frozen decision-table results, exact contraction, invalid-store cases, and fixture-local trusted execution evidence;
- PR G1/G2 when used: exact report-only observation, accepted subjects, and the later authority-only mode transition;
- every guard phase: confirmation that no privileged policy allowlist entry changed.

For an initial existing-debt rule when frozen mechanics are absent, include this five-step evidence chain:

```text
detector-only PR
  -> report-only rule authority PR
  -> conditional frozen-mechanics PR F2
  -> accepted-violation authority PR G1
  -> frozen rule-authority tightening PR G2
```

After F2 has merged, later existing-debt rules reuse its generic mechanics and require the detector-only, report-only authority, exact accepted-violation authority, and frozen-authority tightening PRs from Section 6.8. If F2 was skipped, adding the first later existing-debt rule requires a new evidence-backed plan before accepted-store code is introduced.

## 14. Stop conditions

Stop implementation and report the conflict when any of the following occurs:

1. the current `origin/dev` authority differs materially from the files audited for this plan;
2. the checker and architecture test cannot share detector implementation without changing characterized source evidence;
3. a proposed hard rule would fail on the untouched base or a proposed report-only rule cannot emit stable exact subjects;
4. a required change touches production runtime;
5. a new dependency appears necessary;
6. the trusted-base workflow would execute untrusted head code;
7. the implementation requires simultaneous checker and authority-policy modification;
8. a proposed automatic metric cannot distinguish a semantic owner from a caller with acceptable confidence;
9. compatibility evidence cannot be tied to first-parent history or a release tag;
10. a deterministic violation cannot be assigned a stable identity without line numbers or a broad regex;
11. the proposed accepted-violation store would duplicate privileged permissions, compatibility contracts, or metric snapshots;
12. a rule would require automatic CI refreezing or self-updating authority data;
13. a rules JSON design requires regexes, JavaScript, symbol patterns, or other detector-narrowing executable semantics in authority data;
14. the `dev` pull-request workflows do not run the normal and trusted-base checks after PR 0A;
15. the trusted-base checker cannot evaluate head runtime exclusively with base detectors, base rules, and base policy;
16. an executable change to an active detector or transitive shared helper cannot be staged under a new versioned ID;
17. the pure evaluator would need filesystem, Git, environment, reporting, source-detection, or intended-authority ownership;
18. the initial useful registry would require more than three rules, a third rule kind, or a generic rule DSL;
19. untouched `origin/dev` contains a first-party static import cycle, including a self-import;
20. the canonical-edge detector cannot emit the complete in-scope candidate edge set needed for subset and exact-count evaluation;
21. frozen or accepted-store parsing/evaluation mechanics would be added without a Phase 0-selected materially valuable existing-debt rule.

Do not solve a stop condition by adding a fallback, compatibility branch, duplicate parser, or temporary alternate checker.

## 15. Completion criteria

This implementation is complete only when all of the following are true:

- the normative fitness-function ratchet document is merged and referenced by agent, repository, Web, contributor, and PR entry points;
- the document identifies the intended/source/reflexion model and the adopted external precedents without making an external tool authoritative;
- the PR template collects concrete changed-scope sets, a scope-completeness rationale, and a structured maintainer-decision permalink without requiring repository-wide manual metrics;
- the PR template is a protected guard path, retains every required process anchor, and does not claim that branch protection automatically enforces the manual decision comment;
- every final architecture approval names the exact merge-candidate head SHA, and a later commit invalidates the approval;
- critical capability detection has one executable owner containing no intended paths, counts, modes, eligibility, or exceptions;
- every active detector ID is a reviewed checker-code trust contract; any executable change to it or its shared helpers is staged through detector-add, authority-migrate, and detector-remove PRs;
- detector fixtures and untouched-source comparisons are presented as characterization evidence, not mechanical proof of equivalence;
- pure evaluation has one I/O-free mechanics owner containing no source detection or intended authority;
- intended registered architecture has one inert JSON authority with at most three initial rules, only the two approved kind-specific schemas, and no executable or detector-narrowing semantics;
- every `single-canonical-entry-edge` rule allows preauthorized inactive edges but requires exactly one active allowed edge and rejects zero, unapproved, or parallel active edges;
- the checker, architecture contracts, and ratchet fixture tests consume the shared mechanics while retaining independent expected-result assertions;
- the normative document, PR template, detector, evaluator, fixture test, rules, and optional violation store were protected in the base checker before creation;
- the trusted-base checker never executes a proposed head detector and never uses head rules to authorize the same head runtime;
- the checker reports before/add/remove/after/delta for deterministic inventories;
- registered semantic-authority expansion fails CI;
- registered rules emit observation, mode, baseline relation, authority resolution, and derived decision fields; required rollout modes use non-applicable baseline and authority resolution;
- report-only divergent/absent results report without failing, and hard/report-only rules never consult accepted violations;
- when Phase 0 selects existing debt, conditional PR F2 implements frozen mechanics before G1/G2, new violations fail, retained fixed authority fails, and same-PR exact contraction passes only under the Section 6.6 equality;
- when Phase 0 selects no materially valuable existing-debt rule, PR F2/G1/G2 are skipped, frozen activation fails closed, no accepted-store parsing/evaluation is added, and no violation store exists;
- existing-debt rules require separate detector, report-only authority, frozen-mechanics, accepted-violation authority, and frozen-authority PRs;
- first-party Web dependency cycles, including self-imports, fail CI;
- heuristic signals remain clearly report-only;
- CI does not present registered file-level observations as complete `OE`, `PE`, or `CB`;
- architecture-bearing PRs require explicit maintainer review of concrete before/after sets;
- existing privileged allowlists and dependency restrictions still pass unchanged;
- no Web runtime production file was changed by this initiative;
- no new dependency was added;
- no accepted-violation store exists when the accepted set is empty, or every stored entry is an exact justified base violation;
- `dev` and `main` workflow coverage is verified and required branch checks are confirmed; a merely documented external action remains incomplete;
- all superseded duplicated capability definitions have been removed.

## 16. Expected final architecture

```text
Declared intent
  docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
        |
        +--> short references from AGENTS.md, CLAUDE.md,
        |    gbdraw/web/CLAUDE.md, CONTRIBUTING.md
        |
        +--> .github/pull_request_template.md records human evidence
        |      protected guard path with tested process anchors
        +--> explicit manual maintainer review records architecture approval

Intended machine model
  tools/web-architecture-rules.json   (inert authority)
        |
        +--> tools/web-change-policy.json
        |      owns permitted privileged paths
        |
        +--> tools/web-architecture-violations.json
               optional exact accepted violations only after F2 and G1

Source-derived model
  tools/web-change-source.mjs
        |
        +--> tools/web-architecture-detectors.mjs
               trusted executable detection only
               versioned reviewed evidence contracts
               imports, markers, canonical evidence

Pure reflexion mechanics
  tools/web-architecture-evaluation.mjs
        |
        +--> receives normalized source facts and inert authority
        +--> schema / authority delta / differential sets / SCC
        +--> accepted-set algebra only when conditional F2 is required
        +--> observation / mode / baselineRelation
        |    / authorityResolution / decision
        +--> decision is derived by one pure table from four inputs
        |
        +--> tests/web/architecture-ratchet-fixtures.test.mjs

Reflexion orchestration and enforcement
  tools/check-web-change-budget.mjs
        |
        +--> sole CLI/CI entry; Git and file I/O; stable reporting
        +--> trusted detectors + pure evaluator + base authority
        |
        +--> tests/web/architecture-contracts.test.mjs

CI
  Web change budget and trusted-base workflows on dev and main
        |
        +--> reject deterministic expansion
        +--> permit exact accepted historical violations only when F2/G1/G2 run
        +--> expose contractions
        +--> report heuristic and trend signals separately
```

The intended result is not a smaller module count or a higher synthetic score. It is a codebase in which one semantic change normally requires editing one explicit owner and one canonical path, while compatibility, privileged authority, and deterministic architecture violations cannot expand silently.

## 17. External references

These references justify terminology and mechanism choices. They do not become repository authorities or dependencies.

1. Thoughtworks, “Architectural fitness function,” Technology Radar.
   https://www.thoughtworks.com/en-us/radar/techniques/architectural-fitness-function

2. Gail C. Murphy, David Notkin, and Kevin Sullivan, “Software Reflexion Models: Bridging the Gap Between Source and High-Level Models,” ACM SIGSOFT FSE, 1995.
   https://www.cs.ubc.ca/~murphy/papers/rm/fse95.html

3. ArchUnit User Guide, “Freezing Arch Rules.”
   https://www.archunit.org/userguide/html/000_Index.html#_freezing_arch_rules

4. SonarQube Server documentation, “Clean as You Code” and new-code quality gates.
   https://docs.sonarsource.com/sonarqube-server/2025.2/user-guide/clean-as-you-code/about-new-code/

5. NDepend documentation, “Quality Gates and Build Failure.”
   https://www.ndepend.com/docs/quality-gates

6. ArchSteer, architecture governance and net-new conformance ratchet.
   https://www.archsteer.com/

7. dependency-cruiser, declarative JavaScript dependency validation.
   https://github.com/sverweij/dependency-cruiser

8. CodeScene documentation, change coupling as a historical architecture diagnostic.
   https://codescene.io/docs/guides/technical/change-coupling.html
