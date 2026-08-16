# Web change containment implementation plan

Status: planned; implementation not started

Date: 2026-08-16

## 1. Objective and completion boundary

Turn the Web change-budget and ownership gate into a durable containment
mechanism. An `architecture-change` label must select a larger but finite
profile; it must never remove the absolute budget or integrity checks.

The completed change must:

- enforce production file, gross churn, net addition, module, state-owner, and
  existing-module export limits;
- allow a small stateless leaf module without encouraging existing files to
  grow indefinitely;
- keep uncertain names and state-domain count as review signals rather than
  treating them as architecture proof;
- prevent Web runtime and governance authority from changing together;
- require a falsifiable scope statement and independent oracle in Web runtime
  pull requests;
- let privileged allowlists ratchet downward safely; and
- make the checks unavoidable on protected `main`.

No new Skill, checker, workflow, dependency, generalized policy engine, or
cross-revision test framework is in scope. This plan does not authorize a
commit, push, pull request, merge, or GitHub settings change. Request explicit
authorization when execution reaches an external action.

## 2. Baseline and binding corrections

Reference `main` at plan creation:

```text
43c14ddec49835458f6bac88aa134a7d755c89b1
Merge pull request #338: Add Web change-budget and ownership gate
```

The current branch was `dev` at
`6298b009610a14f1f5484833bd71eae12d344d2b`, with a clean worktree. Its only
work beyond the reference `main` is the focused contraction task document;
the checker behavior has not yet been changed.

At this baseline:

- `architecture-change` drops all ordinary violations;
- line growth uses net additions but has no gross-churn limit;
- every new JavaScript module, export, `create*`, reactive declaration, and
  watcher increase is an ordinary violation;
- governance separation covers a short exact-path set;
- all allowlist path contractions fail;
- the trusted-base workflow reads head files as Git data without executing
  head code;
- no pull-request template exists; and
- GitHub reported `main` as unprotected with required checks disabled.

Phase 0 must recheck these facts against the then-current base. Do not repeat a
phase that has already landed; record equivalent evidence instead.

The following corrections are binding:

1. State-domain count is not a hard budget. One correct live-edit action may
   legitimately cross canonical editor state, History, and mounted SVG/Result.
   Require an atomicity explanation and report cross-domain risk.
2. `Manager`, `Journal`, `Protocol`, `Transaction`, `Lease`, and similar names
   remain report-only. Vocabulary is not architecture.
3. An ordinary pull request may add one stateless leaf module.
4. Removing at least two existing production paths remains the default for a
   consolidation abstraction, not an absolute rule for new capability or a
   security boundary.
5. CI validates the presence and structure of PR evidence; review evaluates
   whether the evidence is true and independent.
6. Commit count is not a stop condition. Cumulative diff, failed oracle, or a
   newly required independent state/rollback/persistence/lifecycle authority
   is a stop condition.
7. Implementation plans and evidence logs are not governance authorities and
   are not classified as such merely because their filenames contain `PLAN`.

## 3. Target enforcement contract

### 3.1 Initial finite profiles

Add `schemaVersion: 2` and `changeBudgets` to
`tools/web-change-policy.json` while retaining the existing owner/importer
allowlists.

| Limit | Ordinary | Architecture |
| --- | ---: | ---: |
| Production files | 8 | 12 |
| Gross churn (`additions + deletions`) | 800 | 1,500 |
| Net additions | 100 | 400 |
| New JavaScript modules | 1 | 2 |
| New state-owner paths | 0 | 1 |
| New exports on existing modules | 0 | 1 |
| New production dependencies | 0 | 0 |
| Changed vendored runtime files | 0 | 0 |
| Added binary runtime files | 0 | 0 |

The profile keys are:

```json
{
  "schemaVersion": 2,
  "changeBudgets": {
    "ordinary": {
      "maxProductionFiles": 8,
      "maxGrossChurn": 800,
      "maxNetAdditions": 100,
      "maxNewModules": 1,
      "maxNewStateOwnerPaths": 0,
      "maxNewExportsOnExistingModules": 0,
      "maxNewProductionDependencies": 0,
      "maxChangedVendorPaths": 0,
      "maxAddedBinaryRuntimePaths": 0
    },
    "architecture": {
      "maxProductionFiles": 12,
      "maxGrossChurn": 1500,
      "maxNetAdditions": 400,
      "maxNewModules": 2,
      "maxNewStateOwnerPaths": 1,
      "maxNewExportsOnExistingModules": 1,
      "maxNewProductionDependencies": 0,
      "maxChangedVendorPaths": 0,
      "maxAddedBinaryRuntimePaths": 0
    }
  }
}
```

The checker must understand the version-1 base during the governance pull
request that introduces version 2. Once version 2 is present on the base,
missing profile sections or keys are non-waivable policy-shape violations.
When Pull Request B runs against a version-1 base, the head checker uses
immutable compatibility defaults equal to the table above and validates the
proposed version-2 policy only as data. The proposed policy must not set the
limits that govern its own pull request.
Review the values after 10-15 representative Web pull requests; adjust them
only through another governance-only pull request. Do not add telemetry or an
adaptive budget.

### 3.2 Counting semantics

Keep the production scope at `gbdraw/web/index.html`, `gbdraw/web/js/**`, and
`gbdraw/web/vendor/**`.

- Gross churn is textual additions plus deletions. Binary paths remain a
  separate count.
- Keep the current no-rename diff behavior. Pure moves and line-ending rewrites
  consume churn because they materially reduce diff reviewability.
- Added `.js`, `.mjs`, and `.cjs` paths count as modules.
- Exports in an added module are covered by the module count; only exports
  newly exposed from an existing path count against the export limit.
- A unique path is a new state-owner path when it gains a top-level `create*`
  in `gbdraw/web/js/app/*.js`, a `ref`/`shallowRef`/`reactive` declaration, or
  additional `watch`/`watchEffect` calls.
- `computed` is derived state. Report new computed values but do not count them
  as independent state owners.
- Privileged owner/importer expansion remains governed by the base allowlist;
  do not count it a second time as state ownership.

The selected profile supplies every limit. Replace the current
`architectureChange ? [] : ordinaryViolations` behavior; the label must never
empty the budget violations.

### 3.3 Hard integrity and governance separation

The following remain non-waivable:

- runtime plus governance in one diff;
- unapproved privileged expansion against the base policy;
- a proposed contraction that excludes an active head owner/importer;
- removal or renaming of a base capability/import-target key; and
- both profiles' absolute limits.

Implement a local `isGovernancePath(path)` predicate, without a glob package,
covering:

```text
AGENTS.md
CLAUDE.md
gbdraw/web/CLAUDE.md
.agents/skills/**
.claude/**
.githooks/**
.github/workflows/**
.github/pull_request_template.md
tools/check-web-change-budget.mjs
tools/web-change-source.mjs
tools/web-change-policy.json
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
tools/check-*policy*
tools/*change-budget*
```

User documentation, ordinary behavior tests, implementation plans, and
evidence logs are negative classification cases.

### 3.4 PR metadata and instruction-level rules

Add `.github/pull_request_template.md` with a delimited Web section:

```text
Change kind:
Primary behavior/invariant:
State domains touched:
Why these domains must change atomically:
Domains explicitly out of scope:
Oracle command or reproducer:
Baseline revision and observed result:
Head revision and observed result:
Normative source or N/A with reason:
```

Offer checkboxes for mounted Result, canonical editor state, History, session,
Worker/Generate/reflow lifecycle, Legend/layout geometry, prepared biological
input/cache, presentation-only, and an explicitly named other domain. The list
is review vocabulary, not a closed ontology.

Extend the existing checker to read `GITHUB_EVENT_PATH` only when a PR event
has a Web runtime diff. Parse the body as data. Reject missing values,
placeholders, and bare `N/A`; report selected domains and warn, but do not fail
solely because three or more domains are selected. A local or non-PR run skips
metadata validation while continuing to check the diff.

Add `edited` to the two existing workflow activity lists. Keep the workflow,
job, and check names unchanged.

Add a short section to `gbdraw/web/CLAUDE.md` that requires:

- a falsifiable shared-cause hypothesis;
- one primary user-visible invariant;
- regression evidence that fails on known-bad and passes on head;
- base/head external equivalence for a behavior-preserving refactor;
- a normative fixture/specification for genuinely new behavior;
- supporting tests not to serve as their own primary oracle;
- immediate stop on hard-budget excess or a newly required independent state,
  rollback, persistence, or lifecycle authority; and
- evidence-backed `IMPLEMENTED`, `NO PRODUCTION CHANGE REQUIRED`,
  `ARCHITECTURE HYPOTHESIS REJECTED`, or `SPLIT REQUIRED` dispositions.

A `SPLIT REQUIRED` result must identify the verified call graph, falsified
hypothesis, independently valuable slices, allowed files, and oracle for each
slice. Arbitrarily serializing one oversized design is not compliance.

## 4. Delivery plan

### External Gate 0: protect `main`

Status: pending; explicit external authorization required

Owner: GitHub Ruleset or branch-protection settings.

Behavior:

- require pull requests;
- require `Web change budget`, `Web base policy (trusted base)`, and the
  repository's other required test checks;
- require an up-to-date head or equivalent merge-queue behavior;
- block force pushes and branch deletion; and
- configure no bypass actor, including administrator bypass.

Evidence:

- authenticated Ruleset/branch-protection JSON or settings capture with the
  branch target, checks, and bypass state;
- public branch metadata reports `protected: true`; and
- a normal PR receives both named policy checks.

Do not test protection by attempting a direct or force push.

Dependency: explicit user authorization.

### Phase 0: revalidate the execution base

Status: pending

Owner: read-only repository and GitHub inspection.

Behavior: record the exact base, dirty files, current policy semantics, and
remote protection state. Replace stale assumptions before editing.

Evidence commands:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git rev-parse origin/main
git diff --stat
node tools/check-web-change-budget.mjs
```

Dependency: none locally; refreshing remote state may require credentials or
network approval.

### Phase 1 / Pull Request A: bootstrap safe allowlist contraction

Status: pending

Owned behavior files:

```text
tools/check-web-change-budget.mjs
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
```

This plan file may record status/evidence but does not expand the behavior
scope. Follow the focused contraction task already stored in `docs/internal/`.

Behavior:

- keep base policy authoritative for expansion;
- parse a changed head policy only as data;
- allow later removal of inactive path entries;
- reject active owner/importer contraction;
- reject missing base capability/import-target keys; and
- preserve trusted-base execution and runtime/guard separation.

Do not modify the real allowlist, policy schema, workflow, runtime, dependency,
or source-parser helper in this phase.

Acceptance evidence:

- revision-based safe contraction PASS;
- active owner and importer contraction FAIL and non-waivable;
- same-change runtime plus policy expansion still fails against the base
  allowlist, not only against runtime/guard separation;
- deletion of an empty base key fails;
- focused tests and working-tree checker pass; and
- policy contents, workflows, dependencies, and Web runtime have no diff.

Commands:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
git diff --exit-code -- tools/web-change-policy.json gbdraw/web
```

Dependency: Phase 0. External Gate 0 must be active before this pull request is
published or merged. Merge this pull request before Pull Request B.

Proposed commit title: `Allow Web ownership policy to ratchet downward`

### Phases 2-4 / Pull Request B: activate bounded containment

Status: pending

Owned behavior files:

```text
tools/check-web-change-budget.mjs
tools/web-change-policy.json
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
.github/pull_request_template.md
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
gbdraw/web/CLAUDE.md
```

This plan file may record status/evidence. Web runtime, root guidance,
dependencies, generated artifacts, and reference outputs are out of scope.

Behavior:

1. Add the version-2 finite profiles, profile selection, gross churn, leaf
   module/export semantics, and state-owner-path count.
2. Complete governance-path classification and make every mixed change
   non-waivable.
3. Add PR metadata validation and the `edited` workflow event without adding a
   job or changing check names.
4. Add only the short stop/oracle section to Web guidance and update the one
   existing policy document.

Acceptance matrix:

| Case | Ordinary | Architecture |
| --- | --- | --- |
| Exact numeric limits | PASS | PASS |
| Any active limit plus one | FAIL | FAIL |
| Ordinary excess within architecture limits | FAIL | PASS |
| Large net-zero rewrite over gross limit | FAIL | FAIL |
| One stateless leaf module | PASS | PASS |
| One new state-owner path | FAIL | PASS |
| Two new state-owner paths | FAIL | FAIL |
| New export in added module | counted only as module | counted only as module |
| New export in existing module | FAIL | one path PASS |
| New `computed` or `Manager`-like name only | PASS + report | PASS + report |
| Runtime plus any governance owner | FAIL | FAIL |
| Unapproved expansion or unsafe contraction | FAIL | FAIL |
| Missing runtime PR metadata | FAIL | FAIL |
| Three domains with a complete reason | PASS + warning | PASS + warning |

Required metadata fixtures include missing block, placeholders, bare `N/A`, a
complete regression, a complete new behavior, a governance-only event, a local
non-PR run, and an event body containing shell-like text that remains inert.

Workflow contract tests must continue to prove one trusted base checkout,
read-only permissions, base-checker execution, and no head checkout or code
execution.

Dependencies: Pull Request A merged; External Gate 0 active before merge.

Proposed commit title: `Enforce bounded Web change containment`

### Phase 5: local and historical acceptance

Status: pending

Owner: test execution and separated diff review; no new behavior edits.

Required commands:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
git diff --stat
```

The historical #337 revisions predate `tools/web-change-policy.json`, so the
new checker cannot be invoked directly against them without inventing a policy
injection mode. Do not add one. Instead:

1. prove the checker behavior with committed temporary-repository fixtures;
2. measure the historical runtime diffs read-only with `git diff --numstat`;
3. record that the first #337 commit is 21 production files,
   `+1,398/-1,048`, gross 2,446; and
4. record that the final #337 diff is 35 production files,
   `+5,501/-2,138`, gross 7,639.

Recompute those measurements from:

```text
base: 3a098d8f25bb308833b7c0d82c6edd27e8e9a485
first commit: ac797e59be47afcc20c370225394da21f03b7b0a
final head: 476fb10abfa1f588c20415cf5c1fd58ef57ae9f0
```

Both historical diffs exceed the architecture absolute limits. Historical
measurements are acceptance evidence, not permanent fixtures containing #337
code.

Final diff evidence must show:

- empty Web runtime, dependency, generated, and reference-artifact diffs;
- tests changed only for the architecture-policy contracts;
- governance changes limited to the declared owners; and
- no duplicate checker, policy, workflow, or instruction owner.

Dependency: Phases 2-4 implemented.

### External Gate 1: merge and verify remote enforcement

Status: pending; explicit external authorization required

Owner: GitHub pull requests, required checks, and protected branch state.

Behavior: merge Pull Request A before Pull Request B, then verify the merged
base enforces schema version 2, the finite architecture profile, metadata, and
all existing integrity checks.

Evidence:

- both pull requests merged in order with required checks green;
- authenticated protection evidence still shows no bypass actor;
- public branch metadata remains protected; and
- a post-merge read-only checker and focused test run pass.

Do not use a direct push, force push, or branch deletion as a test.

Dependencies: explicit publish/merge authorization, External Gate 0, and all
local Phase 5 evidence.

### Optional Pull Request C: remove a stale allowlist entry

Status: optional; not required for plan completion

Only after runtime has independently removed an owner/importer may a separate
guard-only pull request remove its inactive path entry. It must not change
runtime, budget values, policy keys, or checker semantics. Do not invent a
target merely to exercise contraction.

## 5. Evidence ledger

Update a phase only after its named gate passes.

| Phase | Status | Evidence | Deviations | Remaining risk |
| --- | --- | --- | --- | --- |
| External Gate 0 | pending | None | None | `main` remains bypassable until verified |
| Phase 0 | pending | None | None | Baseline may move |
| Phase 1 / PR A | pending | None | None | All contractions still fail |
| Phases 2-4 / PR B | pending | None | None | Architecture waiver and incomplete governance remain |
| Phase 5 | pending | None | None | Local and historical gates not run |
| External Gate 1 | pending | None | None | Local evidence cannot prove remote enforcement |

For each completed phase, record:

```text
Phase:
Behavior implemented:
Evidence command/artifact and result:
Deviations from plan:
Remaining risk:
```

## 6. Completion rule

Mark this plan complete only when:

- all required local phases have exact evidence;
- Pull Requests A and B have merged in order;
- protected `main`, required checks, and the empty bypass set are verified;
- the synthetic profile tests pass and both historical #337 measurements
  exceed the enforced architecture limits;
- Web runtime, dependencies, generated artifacts, and reference outputs have
  no plan-owned diff; and
- no required ledger row remains pending.

Missing credentials or authorization leave the corresponding external gate
pending. They do not convert local results into remote evidence. Pull Request C
does not block completion unless a specific stale entry becomes an explicit
additional objective.
