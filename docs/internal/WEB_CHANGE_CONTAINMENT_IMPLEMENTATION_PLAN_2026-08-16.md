# Web change containment implementation plan

Status: implementation complete and verified at `915c2226`; pull request #339
is open and merge remains pending

Date: 2026-08-16

## 1. Objective and completion boundary

Prevent an `architecture-change` label from removing Web change containment.
Both ordinary and architecture changes must remain inside finite, reviewable
production file and churn limits.

The first durable implementation is intentionally narrow. It must:

- give ordinary and architecture changes fixed file, gross-churn, and
  net-addition limits;
- keep new production dependencies, changed vendored runtime files, and added
  binary runtime files at zero;
- prevent Web runtime and its governance controls from changing together;
- prevent the checker or its source parser from changing together with the
  policy data or policy workflows;
- preserve base-authoritative privileged owner/importer expansion checks;
- allow only safe, later allowlist contraction.

This plan does not attempt to infer architecture quality from module size,
importer count, export count, naming, or reactive syntax. Those signals may be
reported, but they are not merge gates in the initial implementation.

No commit, push, pull-request update, merge, or GitHub settings change is
authorized by this document. Obtain explicit authorization for each external
action.

## 2. Binding simplifications

The previous version of this plan proposed a versioned policy language,
duplicated active limits and compiled ceilings, module-graph and export
cardinality gates, a 500-line ratchet, machine-parsed PR evidence containing
current SHAs, mandatory Code Owner review, and an empty bypass set. Those
mechanisms are removed from the initial scope.

The following decisions are binding:

1. Budget values live in `tools/check-web-change-budget.mjs` for the initial
   release. `tools/web-change-policy.json` remains an owner/importer allowlist;
   it does not gain `schemaVersion`, profiles, or budget data.
2. The architecture label selects a larger finite profile. It never clears
   budget or integrity violations.
3. Existing new-module, export, `create*`, `ref`/`reactive`, watcher, and
   naming scans are report-only. Module-line, importer-cardinality, and
   state-domain gates or scanners are not added.
4. No PR-body parser is added. A workflow records the actual base and head SHAs;
   contributors do not copy changing GitHub event SHAs into the PR body.
5. Branch protection is outside this implementation's completion boundary.
6. Reconsider additional hard gates only after reviewing 10-15 representative
   Web pull requests. Add one only when the observations show a concrete failure
   that the rule detects with acceptably low false positives.

## 3. Enforcement contract

### 3.1 Fixed finite profiles

The checker owns these constants:

| Limit | Ordinary | Architecture |
| --- | ---: | ---: |
| Production files | 8 | 12 |
| Gross churn (additions + deletions) | 800 | 1,500 |
| Net additions | 100 | 400 |

The selected profile governs every Web runtime pull request. Exact-limit cases
pass; any limit plus one fails. The `architecture-change` label selects the
architecture column and has no other budget-waiver behavior.

The following limits are zero and non-waivable in both profiles:

- new production dependencies, including a newly introduced bare production
  import;
- changed files under `gbdraw/web/vendor/`; and
- added binary files in the Web runtime scope.

### 3.2 Counting semantics

Production scope remains:

```text
gbdraw/web/index.html
gbdraw/web/js/**
gbdraw/web/vendor/**
```

- Gross churn is textual additions plus deletions in production scope.
- Net additions are production additions minus production deletions.
- Binary paths are counted separately and do not receive invented line counts.
- Keep `--no-renames` behavior. A move or line-ending rewrite consumes churn
  because it enlarges the review surface.
- Pull-request enforcement uses the exact base and head revisions supplied by
  the workflow. A working-tree run is a local supplement, not PR-range evidence.

### 3.3 Non-waivable integrity rules

Preserve the current runtime/guard separation for the existing guard set:

```text
tools/check-web-change-budget.mjs
tools/web-change-source.mjs
tools/web-change-policy.json
docs/internal/WEB_CHANGE_POLICY.md
tests/web/architecture-contracts.test.mjs
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

A diff containing both a production path and any path in this set fails.

Also reject a diff that changes either checker implementation path:

```text
tools/check-web-change-budget.mjs
tools/web-change-source.mjs
```

together with either authority-data/workflow group:

```text
tools/web-change-policy.json
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

This is a small explicit separation rule, not a generalized governance-path
classifier.

The following existing allowlist rules remain non-waivable:

- the base policy is authoritative for privileged owner/importer expansion;
- a proposed policy cannot authorize a runtime expansion in the same pull
  request;
- a proposed contraction must retain every base capability/import-target key;
  and
- a proposed contraction must still cover every active owner/importer in the
  head runtime.

Safe contraction remains a two-pull-request sequence:

```text
runtime pull request removes an owner/importer use
-> later guard-only pull request removes the inactive allowlist path
```

### 3.4 Report-only review signals

Keep the existing report sections for new modules, exports, `create*` owners,
reactive declarations, watchers, compatibility-like names, session fields, and
cache/token/handle/journal/protocol/manager-like names. They must not contribute
to `enforcedViolations`.

Do not add these proposed scanners or limits in the initial implementation:

- per-module line ceilings or a 500-line ratchet;
- single-consumer or importer-cardinality enforcement;
- per-module or existing-module export cardinality;
- state-owner path budgets;
- CommonJS or re-export cardinality parsing introduced only to support those
  budgets; or
- semantic claims such as `stateless` derived from lexical evidence.

### 3.5 PR evidence

Do not add an exact-SHA PR-body contract or an `edited` workflow trigger for
this plan. The checker summary already records the revisions it evaluated.

A short, human-reviewed PR description should state:

```text
Primary invariant:
Known-bad revision or normative source:
Oracle command or reproducer:
Observed before/after:
Primary state domain:
```

For genuinely new behavior, a normative source may replace a known-bad
revision. CI checks field values only if a later, separately justified task
demonstrates that missing evidence is a recurring problem. Reviewers, not a
field-presence checker, decide whether the oracle is relevant and independent.

## 4. Delivery plan

### Phase 0: baseline and plan reduction

Status: complete

Observed baseline:

```text
main/base: 43c14ddec49835458f6bac88aa134a7d755c89b1
current implementation commit: f91770cb694e72c830c46ca3cfe3b0945c73a80f
```

The old plan introduced no Web runtime change. Its proposed permanent policy
surface was larger than needed for the containment objective, so this revision
removes those unimplemented mechanisms. The completed one-off contraction task
prompt is removed rather than retained as a second execution authority.

### Phase 1: safe allowlist contraction

Status: implementation complete and verified

Owned behavior files:

```text
tools/check-web-change-budget.mjs
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
```

Implemented behavior:

- base policy remains authoritative for expansion;
- an inactive allowlist path may be removed later;
- an active owner/importer exclusion fails;
- a base capability/import-target key cannot disappear; and
- runtime/guard separation remains unchanged.

Recorded local evidence at `f91770cb694e72c830c46ca3cfe3b0945c73a80f`:

- `node --test tests/web/architecture-contracts.test.mjs`: 24/24 passed;
- the checker passed for base
  `43c14ddec49835458f6bac88aa134a7d755c89b1` and head
  `f91770cb694e72c830c46ca3cfe3b0945c73a80f`; and
- the range contained no Web runtime, workflow, dependency, or real policy-data
  change.

The Phase 1 head `665777bda24ab8a9f8634f079f7c01225227c199`
passes the checker against base
`43c14ddec49835458f6bac88aa134a7d755c89b1`. The current combined pull-request
range through `915c2226e645033487b5eb90c890d425e5fb6278` also passes and contains
only this plan and the three behavior-owner files. It contains no Web runtime,
workflow, dependency, policy JSON, generated-artifact, or reference-output
change. Rerun the range if the pull-request base or head changes before merge.

### External Gate A: protect `main`

Status: out of scope by repository-owner direction on 2026-08-16

Branch protection is not required for this implementation or its completion.
Any future settings change remains a separate external action requiring explicit
authorization.

### Phase 2: enforce the minimal finite profiles

Status: implementation complete and verified at
`915c2226e645033487b5eb90c890d425e5fb6278`; merge remains pending in pull
request #339

Dependencies: Phase 1 behavior is present.

The implementation is committed above the Phase 1 changes and is included in
pull request #339. Per the repository owner's direction on 2026-08-16, separate
pull requests and an A-then-B merge order are not completion requirements.

Owned behavior files:

```text
tools/check-web-change-budget.mjs
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
```

Behavior:

1. Add the fixed ordinary and architecture constants from section 3.1.
2. Add gross-churn calculation and reporting.
3. Replace the current `architectureChange ? [] : ordinaryViolations`
   behavior with selected-profile enforcement.
4. Move production dependency, vendor, and binary violations into the
   non-waivable integrity set.
5. Stop enforcing new-module, export, `create*`, reactive, and watcher signals;
   continue reporting them.
6. Add only the explicit checker/source-parser co-change rule from section 3.3.
7. Do not modify `tools/web-change-policy.json`, workflows, Web runtime,
   dependencies, or generated/reference artifacts.

Required focused cases:

| Case | Ordinary | Architecture |
| --- | --- | --- |
| Exact file/gross/net limits | PASS | PASS |
| Any selected profile limit plus one | FAIL | FAIL |
| Ordinary excess within architecture limits | FAIL | PASS |
| Net-zero rewrite above gross limit | FAIL | FAIL |
| New module/export/reactive signal within direct limits | PASS + report | PASS + report |
| New production dependency | FAIL | FAIL |
| Changed vendored runtime file | FAIL | FAIL |
| Added binary runtime file | FAIL | FAIL |
| Runtime plus guard path | FAIL | FAIL |
| Checker/source parser plus policy/workflow | FAIL | FAIL |
| Unapproved expansion or unsafe contraction | FAIL | FAIL |

Verification:

```bash
node --check tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
node tools/check-web-change-budget.mjs --base "<PR base SHA>" --head "<PR head SHA>"
git diff --check "<PR base SHA>"..."<PR head SHA>"
git diff --name-status "<PR base SHA>"..."<PR head SHA>"
git diff --numstat "<PR base SHA>"..."<PR head SHA>"
git status --short --untracked-files=all
```

The final range must contain no policy JSON, workflow, Web runtime, dependency,
generated-artifact, or reference-output change.

Recorded local evidence:

- `node --check tools/check-web-change-budget.mjs`: passed;
- `node --test tests/web/architecture-contracts.test.mjs`: 30/30 passed;
- `node --test tests/web/*.test.mjs`: 229/229 passed;
- ordinary and architecture exact-range checker runs from `43c14dde` through
  `915c2226` both passed with zero production files and zero production churn;
- exact boundary fixtures cover both profiles for file count, gross churn, and
  net additions, including exact-limit and limit-plus-one cases; and
- `git diff --check 43c14dde...915c2226` passed. The range contains only the
  three behavior-owner files and this plan.

Exact pull-request evidence now exists for pull request #339 at base
`43c14ddec49835458f6bac88aa134a7d755c89b1` and head
`915c2226e645033487b5eb90c890d425e5fb6278`.

### Phase 3: historical acceptance

Status: complete

Use read-only `git diff --numstat` measurements rather than adding a historical
policy-injection framework. Reconfirm:

```text
#337 first commit: 21 production files, gross churn 2,446
#337 final diff: 35 production files, gross churn 7,639
```

Both must exceed the architecture profile. Temporary-repository fixtures, not
historical source copies, prove exact checker boundary behavior.

Recomputed read-only evidence:

```text
base: 3a098d8f25bb308833b7c0d82c6edd27e8e9a485
first commit: ac797e59be47afcc20c370225394da21f03b7b0a
  21 production files; 1,398 additions; 1,048 deletions; gross 2,446
final head: 476fb10abfa1f588c20415cf5c1fd58ef57ae9f0
  35 production files; 5,501 additions; 2,138 deletions; gross 7,639
```

Both ranges exceed the architecture limits of 12 production files and 1,500
lines of gross churn. The historical tip is not an ancestor of the current
`HEAD`; no #337 commit or file content was applied to the working tree. The
current pull-request range is also verified.

### Follow-up observation

Status: optional; does not block initial completion

After 10-15 representative Web pull requests, review checker summaries for:

- false positives or repeated near misses in files/gross/net;
- module growth and importer/export/state-owner warning patterns; and
- attempted dependency, vendor, binary, authority, or governance expansion.

Open a separate task for any policy change. Do not add adaptive telemetry. A
warning becomes a hard gate only with recorded examples showing that it detects
a real containment failure and does not predictably encourage wrapper files or
misplaced ownership.

## 5. Evidence ledger

Update a row only after its named evidence exists.

| Phase | Status | Evidence | Remaining risk |
| --- | --- | --- | --- |
| Plan reduction | complete | Replaced the 850-line plan with the reduced contract, removed the completed 458-line task prompt, passed `git diff --check`, and confirmed the delivery diff contains only the plan plus the three behavior-owner files | Delivery base or head changes require another range check |
| Phase 0 | complete | Base `43c14dde`; current implementation commit `f91770cb`; no Web runtime change | Base may advance before merge |
| Phase 1 | implementation complete; merge pending | Exact checker PASS from base `43c14dde` to Phase 1 head `665777bd`; the combined PR #339 range also passes and contains no Web runtime, workflow, dependency, or policy JSON change | Rerun if the PR base or head changes; merge remains external |
| External Gate A | out of scope | Repository owner removed branch protection from this plan's completion boundary on 2026-08-16 | Any future settings change requires separate authorization |
| Phase 2 | implementation complete; merge pending | Checker syntax PASS; focused 30/30; Web Node 229/229; ordinary and architecture exact-range checks PASS for `43c14dde...915c2226`; exact and plus-one fixtures cover all three finite limits | Pull request #339 remains external and unmerged |
| Phase 3 | complete | Recomputed #337: first commit 21 files / gross 2,446; final diff 35 files / gross 7,639 | Historical evidence is complete |

## 6. Completion rule

Mark this plan complete only when:

- the implementation has merged;
- exact base/head evidence exists for the delivered range;
- the focused acceptance matrix passes;
- the historical #337 measurements exceed the enforced architecture profile;
- no plan-owned Web runtime, dependency, generated-artifact, or reference-output
  change exists; and
- no required ledger row remains pending.

The optional 10-15 pull-request observation period does not block initial
completion. Branch protection is outside this plan's completion boundary.
