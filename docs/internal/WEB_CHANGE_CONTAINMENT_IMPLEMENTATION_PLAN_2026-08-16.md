# Web change containment implementation plan

Status: planned; implementation not started

Date: 2026-08-16

## 1. Objective and completion boundary

Turn the Web change-budget and ownership gate into a durable containment
mechanism. An `architecture-change` label must select a larger but finite
profile; it must never remove the absolute budget or integrity checks.

The completed change must:

- enforce production file, gross churn, net addition, module, state-owner, and
  exported-name limits;
- allow one small, single-consumer ordinary module only when it has no detected
  state-owner or privileged-authority signal;
- cap added-module size and prevent existing JavaScript modules from growing
  past a fixed line ceiling or beyond an already-higher baseline;
- keep uncertain names and state-domain count as review signals rather than
  treating them as architecture proof;
- prevent Web runtime and governance authority from changing together;
- prevent a checker or source-parser change from approving
  `tools/web-change-policy.json` or `.github/workflows/**` in the same pull
  request;
- require a falsifiable scope statement and independent oracle in Web runtime
  pull requests, with independent approval required for governance changes;
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
3. An ordinary pull request may add one small, single-consumer module only when
   the checker finds no state-owner, privileged-owner, privileged-importer, or
   re-export signal. This lexical contract does not prove semantic
   statelessness, and neither the checker nor documentation may claim that it
   does.
4. Removing at least two existing production paths remains the default for a
   consolidation abstraction, not an absolute rule for new capability or a
   security boundary.
5. CI validates the presence and structure of PR evidence. A required reviewer
   who did not author the latest reviewable push evaluates whether the evidence
   is true and independent. Checker self-attestation is not sufficient for a
   governance change.
6. Commit count is not a stop condition. Cumulative diff, failed oracle, or a
   newly required independent state/rollback/persistence/lifecycle authority
   is a stop condition.
7. Implementation plans and evidence logs are not governance authorities and
   are not classified as such merely because their filenames contain `PLAN`.
8. Every scope and churn claim uses a recorded base SHA and head SHA. Plain
   working-tree `git diff` output is supplemental evidence and cannot prove a
   pull-request range.
9. The checker/source parser and `tools/web-change-policy.json` or
   `.github/workflows/**` land in separate pull requests. The merged
   trusted-base checker must validate the later policy and workflow activation
   as data.

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
| Lines in each added JavaScript module | 200 | 300 |
| Existing JavaScript module line ceiling | 500 | 500 |
| Exported names in each added module | 3 | 5 |
| New state-owner paths | 0 | 1 |
| Newly exported names on existing modules | 0 | 1 |
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
      "maxAddedModuleLines": 200,
      "existingModuleLineCeiling": 500,
      "maxExportsPerAddedModule": 3,
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
      "maxAddedModuleLines": 300,
      "existingModuleLineCeiling": 500,
      "maxExportsPerAddedModule": 5,
      "maxNewStateOwnerPaths": 1,
      "maxNewExportsOnExistingModules": 1,
      "maxNewProductionDependencies": 0,
      "maxChangedVendorPaths": 0,
      "maxAddedBinaryRuntimePaths": 0
    }
  }
}
```

Pull Request B1 installs immutable compatibility defaults and hard ceilings
equal to the table above without changing the version-1 policy. When Pull
Request B2 later proposes version 2 against that version-1 base, the merged B1
checker uses the compatibility defaults for the pull request and validates the
proposed version-2 policy only as data. The proposed policy must not set the
limits that govern its own pull request.

Before any budget comparison, version-2 policy parsing must fail closed unless:

- `schemaVersion` is exactly the number `2`;
- `changeBudgets`, `ordinary`, and `architecture` are non-array objects;
- the top level contains exactly `schemaVersion`, `changeBudgets`,
  `allowedPrivilegedImporters`, and `allowedPrivilegedOwners`;
- both profiles contain exactly the required keys shown above;
- every limit is a non-negative JavaScript safe integer;
- every architecture limit is greater than or equal to its ordinary limit;
- new-production-dependency, changed-vendor, and added-binary limits remain
  zero in both profiles; and
- no value exceeds the immutable hard ceiling compiled into the merged checker.

Reject missing, unknown, mistyped, negative, fractional, non-safe, or
out-of-order values as non-waivable policy-shape violations. Once version 2 is
present on the base, version 1 and missing profile data are also non-waivable.
Changing a hard ceiling requires a checker-only pull request with independent
approval, followed by a later policy-only pull request; one pull request cannot
change both layers.

Review the values after 10-15 representative Web pull requests; adjust them
within the immutable ceilings only through another policy-only pull request.
That later calibration is not part of initial completion. Open a separate
governance task with its own base/head evidence rather than silently editing
this completed plan. Do not add telemetry or an adaptive budget.

### 3.2 Counting semantics

Keep the production scope at `gbdraw/web/index.html`, `gbdraw/web/js/**`, and
`gbdraw/web/vendor/**`.

- Gross churn is textual additions plus deletions. Binary paths remain a
  separate count.
- Keep the current no-rename diff behavior. Pure moves and line-ending rewrites
  consume churn because they materially reduce diff reviewability.
- Added `.js` and `.mjs` paths count as modules. An added `.cjs` production path
  or detected `module.exports`/`exports.*` assignment is a non-waivable
  integrity violation because the browser application is an ES-module SPA and
  the export scanner does not certify CommonJS surfaces.
- Each added module must satisfy its selected per-module line and exported-name
  limits. Count final public names, including aliases, `default`, and namespace
  re-exports. Any newly introduced export statement that the checker cannot
  reduce to a finite set of public names is a non-waivable violation.
- An ordinary added module must have exactly one literal relative production
  importer in the head graph, must not be re-exported, and must have no detected
  state-owner or privileged-authority signal. Zero importers, multiple
  importers, non-literal loading, or a re-export fail the ordinary profile.
- On an existing module, count each newly exposed public name, not each changed
  path. Two names added to one file count as two. A newly added bare
  `export * from` in any production module is a non-waivable violation; do not
  pretend its public cardinality is known. Named and namespace re-exports use
  their final exported names.
- For an existing JavaScript module, require
  `afterLines <= max(beforeLines, existingModuleLineCeiling)`. A module already
  above 500 lines may change or shrink, but it may not gain lines. This is a
  ratchet, not a demand to split every existing large module in this work.
- A unique path is a new state-owner path when it gains a top-level `create*`
  in `gbdraw/web/js/app/*.js`, a `ref`/`shallowRef`/`reactive` declaration, or
  additional `watch`/`watchEffect` calls.
- `computed` is derived state. Report new computed values but do not count them
  as independent state owners.
- Privileged owner/importer expansion remains governed by the base allowlist;
  do not count it a second time as state ownership.
- State-owner and naming scans are lexical evidence only. Reports and policy
  documentation must say `no detected state-owner signal`, not `stateless`.

The selected profile supplies every limit. Replace the current
`architectureChange ? [] : ordinaryViolations` behavior; the label must never
empty the budget violations.

### 3.3 Hard integrity and governance separation

The following remain non-waivable:

- runtime plus governance in one diff;
- co-change of a checker/source-parser path with either
  `tools/web-change-policy.json` or `.github/workflows/**` in one diff;
- unapproved privileged expansion against the base policy;
- a proposed contraction that excludes an active head owner/importer;
- removal or renaming of a base capability/import-target key;
- invalid policy shape or a proposed value above a hard ceiling;
- deletion or renaming of a required checker, policy, workflow, or architecture
  contract path;
- removal of required workflow events or named checks, permissions broader than
  read-only, or execution/check-out of pull-request head code in the trusted
  base workflow;
- an added or newly introduced CommonJS production surface, a new bare wildcard
  re-export, or newly introduced export syntax that the checker cannot
  enumerate; and
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
.github/CODEOWNERS
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

Add `.github/pull_request_template.md` with exactly one delimited Web section:

```text
<!-- gbdraw-web-evidence:start -->
Change kind:
Primary behavior/invariant:
State domains touched:
Why these domains must change atomically:
Domains explicitly out of scope:
Oracle command or reproducer:
Baseline revision and observed result:
Head revision and observed result:
Normative source or N/A with reason:
<!-- gbdraw-web-evidence:end -->
```

Offer checkboxes for mounted Result, canonical editor state, History, session,
Worker/Generate/reflow lifecycle, Legend/layout geometry, prepared biological
input/cache, presentation-only, and an explicitly named other domain. The list
is review vocabulary, not a closed ontology.

Extend the existing checker to read `GITHUB_EVENT_PATH` only when a PR event
has a Web runtime diff. Parse the body as data. Reject missing values,
duplicate blocks or fields, placeholders, bare `N/A`, no selected domain, and
baseline/head revision fields that do not contain the event's current base and
head SHAs. Treat the text after the first colon on each field line as inert
data, so additional colons and shell-like content remain inert. Report selected
domains and warn, but do not fail solely because three or more domains are
selected. A local or non-PR run skips metadata validation while continuing to
check the diff.

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

### External Gate 0A: protect `main` and require independent review

Status: pending; explicit external authorization required

Owner: GitHub Ruleset or branch-protection settings.

Behavior:

- require pull requests;
- require at least one approving review;
- dismiss stale approvals after a new reviewable push;
- require approval of the most recent reviewable push by someone other than its
  author;
- require `Web change budget`, `Web base policy (trusted base)`, and the
  repository's other required test checks;
- require an up-to-date head or equivalent merge-queue behavior;
- block force pushes and branch deletion; and
- configure no bypass actor, including administrator bypass.

Identify an existing account or team, distinct from the pull-request author,
that can own and review the governance paths. Record that principal before
Pull Request B1 creates `.github/CODEOWNERS`. If no such reviewer exists, the
durable-governance objective is blocked; do not substitute checker
self-approval.

Evidence:

- authenticated Ruleset/branch-protection JSON or settings capture with the
  branch target, required review count, stale-review behavior, latest-push
  approval rule, checks, and bypass state;
- the recorded governance reviewer exists and can review the repository;
- public branch metadata reports `protected: true`; and
- a normal PR receives both named policy checks.

Do not test protection by attempting a direct or force push.

Dependency: explicit user authorization.

### Phase 0: revalidate the execution base

Status: pending

Owner: read-only repository and GitHub inspection.

Behavior: record the exact fetched base, current head, their merge base, dirty
and untracked files, cumulative branch diff, current policy semantics, module
line distribution, and remote protection state. Replace stale assumptions
before editing. The recorded SHA pair, not the mutable name `origin/main`, owns
later scope evidence.

Evidence commands:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git rev-parse origin/main
git merge-base origin/main HEAD
git diff --name-status origin/main...HEAD
git diff --numstat origin/main...HEAD
git diff --stat origin/main...HEAD
node tools/check-web-change-budget.mjs
find gbdraw/web/js -type f -name '*.js' -print0 | xargs -0 wc -l | sort -n
```

At plan review, the tree contained 138 JavaScript modules: 71 were at or below
200 lines, 85 were at or below 300, 100 were at or below 500, and 38 were above
500. Phase 0 must recompute the distribution. The 500-line existing-module rule
is non-retroactive: files already above it become non-growing baselines.

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
scope. Follow
`docs/internal/gbdraw_next_codex_prompt_allow_ownership_policy_contraction.md`
for contraction behavior. This plan supersedes that task's bare `git diff`
commands: pull-request scope evidence must use the recorded base and head SHAs.

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
PR_A_BASE_SHA="$(git merge-base origin/main HEAD)"
PR_A_HEAD_SHA="$(git rev-parse HEAD)"
node tools/check-web-change-budget.mjs --base "$PR_A_BASE_SHA" --head "$PR_A_HEAD_SHA"
git diff --check "$PR_A_BASE_SHA"..."$PR_A_HEAD_SHA"
git diff --name-status "$PR_A_BASE_SHA"..."$PR_A_HEAD_SHA"
git diff --stat "$PR_A_BASE_SHA"..."$PR_A_HEAD_SHA"
git diff --exit-code "$PR_A_BASE_SHA"..."$PR_A_HEAD_SHA" -- \
  tools/web-change-policy.json tools/web-change-source.mjs .github/workflows \
  gbdraw/web package.json package-lock.json pyproject.toml
git status --short --untracked-files=all
```

Record `PR_A_BASE_SHA` and `PR_A_HEAD_SHA` in the ledger. The first three
commands are pre-commit working-tree checks; they do not replace the explicit
revision-range run.

Dependency: Phase 0. External Gate 0A must be active before this pull request
is published or merged. Merge this pull request before Pull Request B1.

Proposed commit title: `Allow Web ownership policy to ratchet downward`

### Phase 2 / Pull Request B1: bootstrap the trusted containment checker

Status: pending

Owned behavior files:

```text
tools/check-web-change-budget.mjs
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
.github/CODEOWNERS
```

This plan file may record status/evidence. Web runtime, root guidance,
policy data, workflows, dependencies, generated artifacts, reference outputs,
and the source-parser helper are out of scope.

Behavior:

1. Add immutable version-1 compatibility defaults and hard ceilings equal to
   section 3.1, but do not modify the real version-1 policy.
2. Add fail-closed version-2 schema validation, selected-profile enforcement,
   gross churn, module-size ratchets, ordinary single-consumer checks,
   exported-name counting, and state-owner-path counting.
3. Reject added CommonJS modules, newly introduced CommonJS export assignments,
   and new bare wildcard re-exports in any production module.
4. Complete governance-path classification, reject co-change of the
   checker/source parser with the JSON policy or workflows, and validate
   required guard/workflow structure from trusted base code rather than relying
   only on head tests.
5. Add PR metadata parsing and fixtures, but activate runtime metadata failure
   only when the base policy is version 2. This avoids an undocumented
   intermediate requirement between B1 and B2.
6. Add `.github/CODEOWNERS` for every governance path using the reviewer
   principal recorded by External Gate 0A.
7. Update `docs/internal/WEB_CHANGE_POLICY.md` with both the active version-1
   compatibility behavior and the version-2 activation condition, so the
   checker never lands with undocumented enforcement.

With a version-1 base, required workflow activity types remain `opened`,
`synchronize`, `reopened`, `labeled`, and `unlabeled`. When the head proposes a
valid version-2 policy, the trusted B1 checker also requires `edited` in both
head workflows. Once version 2 is on the base, all six activity types are
required.

Acceptance matrix:

| Case | Ordinary | Architecture |
| --- | --- | --- |
| Exact numeric limits | PASS | PASS |
| Any active limit plus one | FAIL | FAIL |
| Ordinary excess within architecture limits | FAIL | PASS |
| Large net-zero rewrite over gross limit | FAIL | FAIL |
| One added module at 200/300 lines | PASS | PASS |
| One added module at 201/301 lines | FAIL | FAIL |
| Module with exactly one importer and no owner signal | PASS | PASS |
| Module with zero or two importers | FAIL | PASS |
| Existing module below ceiling ending at 500 lines | PASS | PASS |
| Existing module below ceiling ending at 501 lines | FAIL | FAIL |
| Existing module above ceiling unchanged/shrunk by line count | PASS | PASS |
| Existing module above its baseline by one line | FAIL | FAIL |
| Added module with 3/5 exported names | PASS | PASS |
| Added module export limit plus one | FAIL | FAIL |
| One new state-owner path | FAIL | PASS |
| Two new state-owner paths | FAIL | FAIL |
| One new exported name in an existing module | FAIL | PASS |
| Two new exported names in one existing module | FAIL | FAIL |
| New bare wildcard re-export in any production module | FAIL | FAIL |
| Added `.cjs` or detected CommonJS export assignment | FAIL | FAIL |
| Export syntax whose public names cannot be enumerated | FAIL | FAIL |
| New `computed` or `Manager`-like name only | PASS + report | PASS + report |
| Runtime plus any governance owner | FAIL | FAIL |
| Checker/source parser plus JSON policy or workflow | FAIL | FAIL |
| Required guard deletion or weakened trusted workflow | FAIL | FAIL |
| Unapproved expansion or unsafe contraction | FAIL | FAIL |
| Missing runtime PR metadata | FAIL | FAIL |
| Three domains with a complete reason | PASS + warning | PASS + warning |

Policy-shape fixtures must reject wrong schema versions; missing and unknown
keys; arrays where objects are required; strings, `null`, negative values,
fractions, non-safe values, and exponent overflow; architecture values below
ordinary values; non-zero dependency/vendor/binary limits; and any value above
the compiled hard ceiling. A version-1 base must use the immutable defaults,
and a proposed version-2 policy must never govern its own pull request.

Required metadata fixtures include a missing or duplicate block, duplicate
fields, placeholders, bare `N/A`, no selected domain, stale base/head SHAs, a
complete regression, a complete new behavior, a governance-only event, a local
non-PR run, and an event body containing colons and shell-like text that remains
inert.

Workflow contract tests must continue to prove one trusted base checkout,
read-only permissions, base-checker execution, and no head checkout or code
execution. Trusted-base checker fixtures must also prove that deleting a
required guard, dropping an activity type or named check, broadening workflow
permissions, or combining the checker/source parser with the JSON policy or a
workflow fails even when `architecture-change` is active.

Commands:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
PR_B1_BASE_SHA="$(git merge-base origin/main HEAD)"
PR_B1_HEAD_SHA="$(git rev-parse HEAD)"
node tools/check-web-change-budget.mjs --base "$PR_B1_BASE_SHA" --head "$PR_B1_HEAD_SHA"
git diff --check "$PR_B1_BASE_SHA"..."$PR_B1_HEAD_SHA"
git diff --name-status "$PR_B1_BASE_SHA"..."$PR_B1_HEAD_SHA"
git diff --stat "$PR_B1_BASE_SHA"..."$PR_B1_HEAD_SHA"
git diff --exit-code "$PR_B1_BASE_SHA"..."$PR_B1_HEAD_SHA" -- \
  tools/web-change-policy.json tools/web-change-source.mjs \
  .github/workflows \
  .github/pull_request_template.md gbdraw/web AGENTS.md CLAUDE.md \
  package.json package-lock.json pyproject.toml tests/reference_outputs
git status --short --untracked-files=all
```

Record `PR_B1_BASE_SHA` and `PR_B1_HEAD_SHA`. The forbidden-range check must
show no policy-data, source-parser, workflow, Web, root-guidance, dependency,
generated-artifact, or reference-output diff.

Dependencies: Pull Request A merged; External Gate 0A active before publication
or merge.

Proposed commit title: `Bootstrap trusted Web containment checks`

### External Gate 0B: require governance Code Owner review

Status: pending; explicit external authorization required

Owner: GitHub Ruleset or branch-protection settings and the merged
`.github/CODEOWNERS` file.

Behavior: after Pull Request B1 merges, require Code Owner review in the
protected-branch ruleset. Verify that `.github/CODEOWNERS` assigns every
governance path listed in section 3.3 to the recorded independent reviewer and
that a governance pull request cannot merge without that approval. Keep the
Gate 0A review, check, latest-push, and no-bypass rules.

Evidence:

- authenticated protection settings show required Code Owner review;
- the merged CODEOWNERS content covers every section 3.3 governance path,
  including CODEOWNERS itself; and
- an existing or explicitly authorized harmless governance-only pull request
  requires approval from the recorded owner after its latest push.

Do not test this by weakening a production rule. Use settings evidence and a
harmless governance-only pull request or an existing qualifying review.

Dependency: explicit user authorization and Pull Request B1 merged.

### Phases 3-4 / Pull Request B2: activate version-2 containment

Status: pending

Owned behavior files:

```text
tools/web-change-policy.json
tests/web/architecture-contracts.test.mjs
.github/pull_request_template.md
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
gbdraw/web/CLAUDE.md
```

This plan file may record status/evidence. The checker, source-parser helper,
`.github/CODEOWNERS`, Web runtime, root guidance, dependencies, generated
artifacts, and reference outputs are out of scope.

Behavior:

1. Write the exact version-2 profiles from section 3.1. The merged B1 checker,
   not the proposed policy, validates the policy and governs this pull request.
2. Add the delimited PR evidence template and add `edited` to both existing
   workflow activity lists without adding a workflow, job, or check name.
3. Update only the existing workflow contract assertions needed for `edited`
   and the canonical version-2 policy shape. Do not change checker fixtures or
   enforcement semantics in this pull request.
4. Add the short stop/oracle section to `gbdraw/web/CLAUDE.md`. The existing
   policy document already describes the conditional activation from B1 and
   must not change in this pull request.

Acceptance evidence:

- the trusted B1 checker accepts the exact proposed policy; the committed B1
  fixtures already prove rejection of malformed values and hard-ceiling
  excess;
- checker/source-parser and Web runtime diffs are empty;
- both workflows retain their names, jobs, required events, read-only
  permissions, and trusted-base execution contract while adding `edited`;
- the template and guidance match the metadata fixtures already merged in B1;
- the required test commands pass; and
- an independent Code Owner approves the exact head after the latest
  reviewable push.

Record `PR_B2_BASE_SHA` and `PR_B2_HEAD_SHA`. Use them for every diff, checker,
and scope claim; plain working-tree output remains supplemental.

Commands:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
node tools/check-web-change-budget.mjs
PR_B2_BASE_SHA="$(git merge-base origin/main HEAD)"
PR_B2_HEAD_SHA="$(git rev-parse HEAD)"
node tools/check-web-change-budget.mjs --base "$PR_B2_BASE_SHA" --head "$PR_B2_HEAD_SHA"
git diff --check "$PR_B2_BASE_SHA"..."$PR_B2_HEAD_SHA"
git diff --name-status "$PR_B2_BASE_SHA"..."$PR_B2_HEAD_SHA"
git diff --stat "$PR_B2_BASE_SHA"..."$PR_B2_HEAD_SHA"
git diff --exit-code "$PR_B2_BASE_SHA"..."$PR_B2_HEAD_SHA" -- \
  tools/check-web-change-budget.mjs tools/web-change-source.mjs \
  docs/internal/WEB_CHANGE_POLICY.md .github/CODEOWNERS \
  gbdraw/web/index.html gbdraw/web/js gbdraw/web/vendor \
  AGENTS.md CLAUDE.md package.json package-lock.json pyproject.toml \
  tests/reference_outputs
git status --short --untracked-files=all
```

Dependencies: Pull Request B1 merged and External Gate 0B active before this
pull request is published or merged.

Proposed commit title: `Activate bounded Web change containment`

### Phase 5: local and historical acceptance

Status: pending

Owner: test execution and separated diff review; no new behavior edits.

Required commands:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
node tools/check-web-change-budget.mjs
ACCEPTANCE_BASE_SHA="<recorded PR_B2 base SHA>"
ACCEPTANCE_HEAD_SHA="<recorded PR_B2 head SHA>"
node tools/check-web-change-budget.mjs --base "$ACCEPTANCE_BASE_SHA" --head "$ACCEPTANCE_HEAD_SHA"
git diff --check "$ACCEPTANCE_BASE_SHA"..."$ACCEPTANCE_HEAD_SHA"
git diff --name-status "$ACCEPTANCE_BASE_SHA"..."$ACCEPTANCE_HEAD_SHA"
git diff --numstat "$ACCEPTANCE_BASE_SHA"..."$ACCEPTANCE_HEAD_SHA"
git diff --stat "$ACCEPTANCE_BASE_SHA"..."$ACCEPTANCE_HEAD_SHA"
git status --short --untracked-files=all
```

Replace the two angle-bracket values with recorded immutable SHAs before
running the block. Do not recompute them from a branch name after review. Audit
the A, B1, and B2 ranges separately; a clean worktree after commit or merge is
not scope evidence.

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

- for each of the A, B1, and B2 SHA ranges, empty out-of-scope runtime,
  dependency, generated, and reference-artifact diffs;
- tests changed only in the phase that owns the corresponding contract;
- governance changes limited to each phase's declared owners;
- no A, B1, or B2 range co-changed a checker/source-parser path with the JSON
  policy or a workflow;
- no duplicate checker, policy, workflow, or instruction owner.

Dependency: Phases 2-4 implemented.

### External Gate 1: merge and verify remote enforcement

Status: pending; explicit external authorization required

Owner: GitHub pull requests, required checks, and protected branch state.

Behavior: merge Pull Request A, Pull Request B1, and Pull Request B2 in that
order, then verify the merged base enforces schema version 2, immutable hard
ceilings, the finite architecture profile, metadata, required governance
review, and all existing integrity checks.

Evidence:

- all three pull requests merged in order with required checks and required
  independent reviews;
- authenticated protection evidence still shows required Code Owner review,
  latest-push approval, and no bypass actor;
- public branch metadata remains protected;
- a post-merge read-only checker and focused test run pass; and
- a normal Web runtime pull request receives both named checks and enforces the
  version-2 metadata contract.

Do not use a direct push, force push, or branch deletion as a test.

Dependencies: explicit publish/merge authorization, External Gates 0A and 0B,
and all local Phase 5 evidence.

### Optional Pull Request C: remove a stale allowlist entry

Status: optional; not required for plan completion

Only after runtime has independently removed an owner/importer may a separate
guard-only pull request remove its inactive path entry. It must not change
runtime, budget values, policy keys, or checker semantics. Do not invent a
target merely to exercise contraction. It depends on Pull Request B2 and both
external protection gates and requires the same Code Owner approval.

## 5. Evidence ledger

Update a phase only after its named gate passes.

| Phase | Status | Evidence | Deviations | Remaining risk |
| --- | --- | --- | --- | --- |
| External Gate 0A | pending | None | None | `main` lacks verified independent-review enforcement |
| Phase 0 | pending | None | None | Baseline may move |
| Phase 1 / PR A | pending | None | None | All contractions still fail |
| Phase 2 / PR B1 | pending | None | None | The trusted checker does not enforce bounded profiles |
| External Gate 0B | pending | None | None | Governance paths lack verified Code Owner enforcement |
| Phases 3-4 / PR B2 | pending | None | None | Version 2 and metadata remain inactive |
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
- Pull Requests A, B1, and B2 have merged in order;
- protected `main`, required checks, independent review, Code Owner review,
  latest-push approval, and the empty bypass set are verified;
- the fail-closed schema, module-size ratchet, single-consumer, exported-name,
  profile, and governance-integrity fixtures pass;
- both historical #337 measurements exceed the enforced architecture limits;
- the ledger records immutable base/head SHAs and range evidence for A, B1, and
  B2; a clean worktree is not accepted as a substitute;
- Web runtime, dependencies, generated artifacts, and reference outputs have
  no plan-owned diff; and
- no required ledger row remains pending.

Missing credentials or authorization leave the corresponding external gate
pending. They do not convert local results into remote evidence. Pull Request C
does not block completion unless a specific stale entry becomes an explicit
additional objective.
