# gbdraw architecture fitness-function ratchet implementation plan

Status: revised proposal, revision 3
Created: 2026-08-17
Revised: 2026-08-17
Supersedes: revision 2 of this plan and `GBDRAW_ARCHITECTURE_MATURITY_RATCHET_IMPLEMENTATION_PLAN_2026-08-17.md` revision 1
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
2. one executable detector catalog that derives registered source facts but contains no owner, threshold, enforcement, or exception authority;
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
  = conforming
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

Those definitions can drift. Executable detection needs one implementation owner, intended semantic-owner and canonical-path constraints need a separate inert authority owner, and the allowlist of permitted privileged paths must remain in `tools/web-change-policy.json`. A single module must not combine regexes or scan functions with allowed evidence paths, target counts, enforcement modes, or baseline eligibility.

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

The automated agent or pull-request author may prepare the evidence but cannot treat its own declaration as architecture approval. In this repository, approval means an explicit recorded decision by the maintainer who owns the repository architecture; it does not imply that an independent second maintainer must exist.

For a real new capability, one new semantic owner may be added without increasing owner excess. The new capability must begin with one owner and one canonical path.

A change is **architecture-convergent** when behavior and supported capability are preserved or improved, the vector does not worsen, and at least one debt component decreases.

### 6.6 Differential architecture-debt rule

The ratchet compares the proposed head to the target-branch base. It does not require an unrelated historical cleanup in every PR.

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

The merge gate requires \(N_r = \varnothing\), and every entry in \(F_r\) is stale authority that must be removed through the safe-contraction path.

For a `report-only` rule, the checker reports `V_r(head)` as signals and does not compare it with the accepted-violation store:

```text
signals = V_r(head)
decision = PASS when signals is empty, otherwise REPORT
```

A `hard` rule is never baseline-eligible and never consults the accepted-violation store.

An implementation PR may not add entries to \(A_r\). A later guard-only preauthorization PR may expand \(A_r\) only with explicit evidence. A safe contraction may remove entries in \(F_r\), subject to the same trusted-base and guard-separation protections used by the privileged-path policy.

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
- an accepted-set expansion is a separate guard-only authority decision;
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
| change a rule's detector ID, kind, or subject encoding | reinterpretation; replace through a separately staged detector/rule sequence |
| add an allowed evidence path | expansion or preauthorization |
| increase a target count | relaxation |
| `hard` to `frozen` or `report-only` | relaxation |
| `frozen` to `report-only` | relaxation |
| `baselineEligible: false` to `true` | relaxation |
| remove an allowed evidence path | contraction, only when head source no longer requires it |
| reduce a target count | tightening, only when head source conforms |
| `report-only` to `frozen` or `hard` | tightening, only when the destination-mode gate passes |
| `frozen` to `hard` | tightening, only when no accepted violation remains |
| `baselineEligible: true` to `false` | tightening, only when no accepted violation remains |

Every new rule follows a detector-first path:

1. add or extend a stable rule-specific detector in a checker-only PR, with characterization tests and no rule-authority change;
2. add the rule to the inert authority registry in a later authority-only PR;
3. activate only the mode justified by untouched-base evidence.

The rules JSON references a stable detector ID. It must not contain regular expressions, JavaScript source, executable expressions, symbol patterns, or other data that can narrow detection. Rule-specific matching code remains in the trusted detector module.

#### Clean-base rule

When the untouched target branch has no violations:

1. merge the detector-only PR;
2. add the rule in an authority-only PR as `hard`;
3. let the trusted-base checker validate the candidate inert JSON against the already-merged detector and untouched base;
4. do not create a violation-store entry.

#### Existing-debt rule

When the untouched target branch has exact existing violations:

1. **Detect:** merge the rule-specific detector in a checker-only PR.
2. **Observe:** add the rule in an authority-only PR with `report-only` mode and stable violation identities.
3. **Record:** after the report-only rule has produced untouched-base evidence, add only the observed exact base violations to `tools/web-architecture-violations.json` in a separate authority-only PR.
4. **Enforce:** in a later authority-only PR, change only the rule's mode from `report-only` to `frozen`. Generic mode behavior is already implemented; no checker edit accompanies this transition.
5. **Contract:** when runtime work fixes an accepted violation, remove the corresponding store entry in the same PR through the narrow safe-contraction path.

The registry schema must distinguish at least:

```text
hard         cannot be baselined
frozen       exact accepted base violations allowed; new violations fail
report-only  produces evidence but does not authorize merge failure by itself
```

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
- a registered canonical path with more paths than its target;
- an unapproved performance or scientific-output regression where a gate exists.

### 7.2 Frozen deterministic violations

These are exact violations of valuable deterministic rules that existed when the rule was introduced and could not be fixed immediately.

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

Do not encode observation, enforcement mode, baseline relation, and decision in one expanding enum. For each registered deterministic rule and relevant subject, the checker emits four orthogonal fields:

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

decision:
  PASS
  REPORT
  FAIL
```

`decision` is derived from the other fields by one pure function; it is not independently assigned or stored. A fixed accepted subject is represented by `baselineRelation: FIXED` and `decision: FAIL`; `FIXED` is not a source observation.

The required decision table is:

| Mode and observation | Baseline relation | Decision |
| --- | --- | --- |
| `HARD` + `CONFORMING` | `NOT_APPLICABLE` | `PASS` |
| `HARD` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `FAIL` |
| `REPORT_ONLY` + `CONFORMING` | `NOT_APPLICABLE` | `PASS` |
| `REPORT_ONLY` + `DIVERGENT` or `ABSENT_REQUIRED` | `NOT_APPLICABLE` | `REPORT` |
| `FROZEN` + `CONFORMING`, with no fixed accepted subject | `NOT_APPLICABLE` | `PASS` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `ACCEPTED` | `PASS` |
| `FROZEN` + `DIVERGENT` or `ABSENT_REQUIRED` | `NEW` | `FAIL` |
| `FROZEN`, accepted subject no longer present | `FIXED` | `FAIL` until contracted |

Any combination outside this table is an internal checker error and fails closed. In particular, hard and report-only rules must not consult the accepted-violation store.

The output must identify:

- rule key and stable subject when applicable;
- intended owner/path/relation from the base rule authority;
- observed source evidence derived by the base detector;
- the four fields above;
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
| Intended semantic-owner/canonical-path rules and enforcement modes | `tools/web-architecture-rules.json` |
| Exact accepted deterministic violations, only when needed | `tools/web-architecture-violations.json` |
| Permitted privileged operator and importer paths | `tools/web-change-policy.json` |
| Web change inventory, delta report, and deterministic gates | `tools/check-web-change-budget.mjs` |
| Executable architecture invariants and checker fixtures | `tests/web/architecture-contracts.test.mjs` |
| Policy workflow and guard-separation explanation | `docs/internal/WEB_CHANGE_POLICY.md` |
| Human change declaration and evidence | `.github/pull_request_template.md` |
| Human architecture approval | Explicit recorded maintainer review |
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
- combine executable detector logic with intended owner paths, target counts, enforcement modes, baseline eligibility, or accepted violations;
- put regular expressions, JavaScript source, executable expressions, or detector-narrowing symbol patterns in the inert rule registry;
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
tools/web-architecture-detectors.mjs
tools/web-architecture-rules.json
tools/web-architecture-violations.json
```

Add only executable implementation paths to `checkerImplementationPaths`:

```text
tools/web-architecture-detectors.mjs
```

Add only inert authority paths to `authorityPaths`:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
tools/web-architecture-rules.json
tools/web-architecture-violations.json
```

Do not classify `tools/web-architecture-rules.json` as checker implementation and do not classify `tools/web-architecture-detectors.mjs` as authority.

##### `tests/web/architecture-contracts.test.mjs`

Add fixture coverage proving:

- runtime cannot co-change with any future guard path;
- checker implementation cannot co-change with any future authority path;
- the future files need not exist for the checker to pass;
- adding an unrelated unregistered path does not acquire authority accidentally.

#### Verification

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

Run the checker against explicit base/head refs before handoff.

#### Acceptance criteria

- every future authority path is protected before it exists;
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
   - whether each proposed rule is clean-base, existing-debt, or not deterministically enforceable;
   - exact test commands and results.
6. Reconcile definitions duplicated between the checker and architecture tests before extraction.
7. For each human-reviewed changed capability example, record concrete before/after owner, canonical-path, and compatibility-path sets rather than totals alone.
8. Verify both named Web checks run for a `dev` pull request and verify or record the external branch-protection requirement.
9. Do not create a rule registry or accepted-violation store in this phase.

#### Gate

Phase 0 passes only when:

- current checker and architecture contracts pass on the recorded base;
- the duplicated detection definitions are fully inventoried;
- every proposed detector has a stable ID and produces deterministic source evidence;
- every proposed rule is classified as clean-base, report-only existing-debt, or not deterministically enforceable;
- proposed detector code contains no intended paths, counts, modes, eligibility, or exceptions;
- proposed rules contain no executable expressions or detector-narrowing patterns;
- `dev` PR CI coverage is observed, and missing branch-protection authority remains explicitly unresolved rather than claimed complete.

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
- trusted-base evaluation order and prohibition on executing head detectors;
- hard, frozen, and report-only modes with the orthogonal result fields;
- accepted-violation protocol and prohibition on CI refreezing;
- examples listed in revision 2, updated to use explicit sets rather than unsupported totals;
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
- Maintainer architecture decision link, required before merge:

- [ ] No second parser, normalizer, state mirror, lifecycle owner, or canonical pipeline was added.
- [ ] Every superseded implementation was removed in this change.
- [ ] No accepted deterministic violation was added in this implementation PR.
```

Do not ask for unsupported repository-wide `OE`, `PE`, or `CB` numbers. A completed author declaration is evidence, not architecture approval.

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
- the template collects concrete sets and separates author evidence from maintainer approval;
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
  'semantic-owner.render-request': detectRenderRequestAuthority,
  'canonical-path.session-request': detectCanonicalSessionRequestPath
});
```

It may contain masked-source regexes and import-graph logic. It must not contain:

- allowed evidence paths;
- target counts;
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
- detector IDs are stable and rule-specific;
- current source matches and policy coverage are byte-for-byte or set-for-set equivalent to the base characterization;
- no intended paths, thresholds, modes, or exceptions moved into detector code;
- no runtime, inert authority, policy allowlist, workflow, or dependency changes.

#### Proposed commit title

```text
refactor: centralize Web architecture detectors
```

### PR C — Add inert rule-schema and trusted candidate validation

Extend the existing checker; do not add a second checker and do not add the rules file yet.

#### Rule schema support

Support an optional strict JSON file at:

```text
tools/web-architecture-rules.json
```

Absence means no new registry rules during the staged rollout. Schema version 1 must define a `rules` array stored in ascending rule-key order. Kind-specific path arrays must also use canonical sorted order. Each rule includes only inert intended-architecture data such as:

```json
{
  "key": "semantic-owner.render-request",
  "kind": "semantic-authority",
  "detector": "semantic-owner.render-request",
  "allowedEvidencePaths": ["services/session-request.js"],
  "targetLocationCount": 1,
  "enforcement": "hard",
  "baselineEligible": false
}
```

The exact kind-specific schema may differ, but it must:

- reject duplicate rule keys and duplicate paths;
- reject unknown detector IDs and unsupported kinds or modes;
- reject hard rules marked baseline-eligible and frozen rules marked ineligible;
- reject paths outside the normalized `gbdraw/web/js/` source scope;
- reject unknown fields;
- reject regex source, JavaScript source, executable expressions, wildcard paths, and detector-narrowing symbol patterns;
- require every referenced detector to provide a stable subject encoder for its violation-producing rule kind.

#### Trusted-base evaluation order

After PR E activates the registry, runtime conformance always uses base detector + base rules + base privileged policy. Parse proposed head rules only as inert JSON and classify their authority delta. Never import or execute a proposed head detector module.

During PR C and PR D, the current inline intended constraints substitute for base rules and remain the sole active runtime authority. Candidate rules are evaluated against untouched base source for admission evidence but do not participate in runtime authorization until PR E removes the inline constraints and activates the unchanged base JSON. This staged inactive file is not a second active rule owner.

Candidate validation must reject:

- rules referring to detectors absent from the trusted base;
- an authority expansion combined with runtime, checker, detector, test-rule, policy, or workflow changes;
- a relaxation combined with any implementation change;
- a rule whose claimed clean base is actually divergent;
- a transition outside the authority-direction table in Section 6.8.

An authority contraction or tightening must still conform on head source under the trusted base detector. Do not let a head rule change authorize the same head runtime.

#### Tests and documentation

- Add fixture tests for malformed and executable-looking JSON, unknown detectors, direction classification, and self-authorization attempts.
- Prove the trusted checker ignores a malicious head detector.
- Update `docs/internal/WEB_CHANGE_POLICY.md` to explain detector/rule/policy separation and candidate inert-data handling.

#### Verification and acceptance

```bash
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
git diff --check
```

- candidate rules can be validated without executing head code;
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

- use `hard` for a clean-base invariant;
- use `report-only` for a valuable rule with exact existing debt;
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
- every allowed evidence path and target count matches recorded intended architecture;
- every hard rule passes on untouched `origin/dev`;
- every report-only rule exposes exact stable subjects without failing;
- the PR contains inert authority data only.

#### Proposed commit title

```text
chore: declare initial Web architecture rules
```

### PR E — Activate the base rule authority and normalized result model

After PR D merges, update the checker and tests without changing the rules JSON.

#### Checker behavior

- Consume the unchanged base `tools/web-architecture-rules.json` as intended architecture.
- Remove the superseded inline intended constraints in the same PR.
- Keep `tools/web-change-policy.json` as the sole privileged path-permission owner.
- Implement the four orthogonal result fields and one pure decision function from Section 7.5.
- Fail closed on impossible field combinations.
- Report rule and subject rows in stable rule-key/subject order.
- Count file-level registered authority locations without calling them complete `OE`.
- Continue selected canonical-path enforcement through explicit registered evidence, not whole-program path inference.

#### Tests

Use table-driven tests for every allowed result combination and representative forbidden combinations. Prove:

- hard conforming passes and hard divergent/absent fails;
- report-only conforming passes and report-only divergent/absent reports without failing;
- accepted-store state is not consulted by hard or report-only modes;
- source observation, mode, baseline relation, and decision cannot contradict one another;
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
- no temporary dual active authority remains;
- no runtime, rules authority, policy allowlist, workflow, or dependency changes.

#### Proposed commit title

```text
ci: enforce declared Web architecture rules
```

### PR F — Add deltas, frozen violations, and cycle enforcement

Extend the existing checker without changing runtime, rules authority, privileged policy, or workflows.

#### Differential reports

For each measurable inventory, report `before`, `added`, `removed`, `after`, and `delta`. At minimum include:

- Registered Authority Location Count;
- registered canonical-contract results;
- accepted, new, and fixed deterministic violation subjects;
- privileged operator permission entries by capability;
- privileged importer permission entries by target;
- production modules, exports, `create*` declarations, reactive declarations, watcher calls;
- compatibility-like and resource-like names and session object keys;
- first-party static import cycles;
- production files, additions, deletions, gross churn, and net additions.

Label heuristic categories `report-only`. Never label these inventories complete `OE`, `PE`, or `CB`. Add a compact summary table to terminal output and `GITHUB_STEP_SUMMARY`.

#### Frozen accepted-violation protocol

Support optional `tools/web-architecture-violations.json` with the strict schema from Section 6.7:

- absence means an empty set;
- head data is inert JSON and is never written by CI;
- rule keys and subjects must match the base rules and stable rule-specific encoders;
- wildcard, duplicate, invented, line-number, malformed, or unknown-rule entries fail;
- new frozen violations fail;
- fixed accepted entries use `baselineRelation: FIXED` and fail until contracted;
- deleting the final entry requires deleting the file;
- store expansion is permitted only in a separate authority-only PR for exact untouched-base signals from an already-merged report-only baseline-eligible rule;
- same-PR runtime plus store contraction is permitted only when removed entries exactly equal fixed head subjects and no authority expands or other guard changes;
- the store never authorizes privileged paths, compatibility contracts, cycles, dependencies, or trusted-base integrity violations.

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
12. absent store behaves as empty;
13. frozen accepted/new/fixed subjects produce the specified field combinations and decisions;
14. hard and report-only rules never use the store;
15. invalid store entries fail;
16. pure store contraction requires exact fixed subjects;
17. guard-only expansion accepts only exact untouched-base report-only signals;
18. a proposed head detector or rule cannot authorize head runtime.

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
- every registered rule emits the normalized result fields;
- new/fixed frozen violations fail correctly;
- self-loops and multi-node cycles fail;
- report-only signals remain non-blocking;
- no accepted store is created unless the later conditional authority phase needs it;
- no runtime, rule, policy, workflow, or dependency changes.

#### Proposed commit title

```text
ci: complete the architecture fitness-function ratchet
```

### Conditional PR G1 and PR G2 — Freeze a rule with existing debt

Skip both PRs when every registered rule is clean or remains intentionally report-only.

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

After G1 merges, change only `tools/web-architecture-rules.json`. Change the selected rule's mode from `report-only` to `frozen`; do not edit checker code, tests, detectors, or the store. The generic frozen behavior was already implemented in PR F.

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

For the first architecture-bearing runtime PR after PR F or G2:

- inspect the terminal and GitHub summaries manually;
- confirm no false hard failure and no silent report-only failure;
- confirm the author enumerates concrete before/after sets and scope completeness;
- record the maintainer architecture decision;
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
| `tools/check-web-change-budget.mjs` | 0B, B, C, E, F | Path protection, detector use, inert authority validation, rule enforcement, results, deltas, violations, and cycles |
| `tests/web/architecture-contracts.test.mjs` | 0A, 0B, B, C, E, F | Branch coverage, separation, characterization, schema, result, violation, and cycle fixtures |
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
- a head detector, rule, or accepted-violation change attempts to authorize the same head runtime;
- no explicit maintainer architecture decision is recorded before merge.

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
- confirmation that no unrelated files were modified.

Additional evidence by phase:

- PR 0A: observed `dev` and `main` workflow triggers, unchanged trusted-base checkout model, and branch-protection state;
- PR 0B: future path classification matrix and self-authorization fixture results;
- PR B: detector before/after characterization, policy-key coverage, and confirmation that detector code contains no authority;
- PR C: inert schema fixtures, authority-direction classification, and proof that head detector code is not executed;
- PR D: exact rule inventory and untouched-base result for every hard or report-only rule;
- PR E: normalized result decision-table evidence and removal of superseded inline authority;
- PR F: before/add/remove/after/delta report, accepted/new/fixed subjects, self-loop and multi-node cycle results, and fixture-local trusted execution evidence;
- PR G1/G2 when used: exact report-only observation, accepted subjects, and the later authority-only mode transition;
- every guard phase: confirmation that no privileged policy allowlist entry changed.

For any future rule with existing debt, also include the four-step evidence chain:

```text
detector-only PR
  -> report-only rule authority PR
  -> accepted-violation authority PR
  -> frozen rule-authority tightening PR
```

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
15. the trusted-base checker cannot evaluate head runtime exclusively with base detectors, base rules, and base policy.

Do not solve a stop condition by adding a fallback, compatibility branch, duplicate parser, or temporary alternate checker.

## 15. Completion criteria

This implementation is complete only when all of the following are true:

- the normative fitness-function ratchet document is merged and referenced by agent, repository, Web, contributor, and PR entry points;
- the document identifies the intended/source/reflexion model and the adopted external precedents without making an external tool authoritative;
- the PR template collects concrete changed-scope sets, a scope-completeness rationale, and a maintainer-decision link without requiring repository-wide manual metrics;
- critical capability detection has one executable owner containing no intended paths, counts, modes, eligibility, or exceptions;
- intended registered architecture has one inert JSON authority containing no executable or detector-narrowing semantics;
- the checker and architecture contracts consume the same detector implementation while retaining independent expected-result assertions;
- the normative document, detector, rules, and optional violation store were protected in the base checker before creation;
- the trusted-base checker never executes a proposed head detector and never uses head rules to authorize the same head runtime;
- the checker reports before/add/remove/after/delta for deterministic inventories;
- registered semantic-authority expansion fails CI;
- registered rules emit orthogonal observation, mode, baseline relation, and derived decision fields;
- report-only divergent/absent results report without failing, and hard/report-only rules never consult accepted violations;
- new deterministic violations fail and fixed accepted violations are exposed as failing stale authority until contracted;
- existing-debt rules require separate detector, report-only authority, accepted-violation authority, and frozen-authority PRs;
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
        +--> explicit maintainer review records architecture approval

Intended machine model
  tools/web-architecture-rules.json   (inert authority)
        |
        +--> tools/web-change-policy.json
        |      owns permitted privileged paths
        |
        +--> tools/web-architecture-violations.json
               optional exact accepted violations only

Source-derived model
  tools/web-change-source.mjs
        |
        +--> tools/web-architecture-detectors.mjs
               trusted executable detection only
               imports, markers, canonical evidence

Reflexion and enforcement
  tools/check-web-change-budget.mjs
        |
        +--> observation / mode / baselineRelation / decision
        +--> decision is derived by one pure table
        |
        +--> tests/web/architecture-contracts.test.mjs

CI
  Web change budget and trusted-base workflows on dev and main
        |
        +--> reject deterministic expansion
        +--> permit only exact accepted historical violations
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
