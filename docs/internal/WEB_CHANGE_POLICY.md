# Web change, review, and delivery policy

Status: normative repository policy

This document defines CI admission, Web architecture review, integrated `dev`
staging, release promotion, and deployment semantics. The
[architecture fitness-function ratchet](./ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
owns architecture definitions, arithmetic, and exception evidence. Workflow and
checker files implement parts of this policy; their current behavior is not a
substitute for this authority.

The [Product Impact Ratchet](./PRODUCT_IMPACT_RATCHET.md) owns planned
architecture-to-behavior traceability, product-impact classification, Decision
Pack semantics, and product-decision eligibility. Its rollout uses the trusted
admission and self-authorization boundaries in this document.

## Delivery lifecycle

gbdraw uses four distinct boundaries:

```text
work branch -> pull request to dev
  fast implementation admission

push to dev exact SHA
  complete integrated staging

dev -> main pull request
  trusted release admission

push to main exact SHA
  release verification, build, and deployment
```

A pull request to `dev` answers whether one candidate is safe to integrate. Its
required path is intentionally fast. Heavy jobs may still run on the pull
request during the transition, but complete supported-version, slow, browser,
and Gallery evidence belongs to the exact integrated `dev` SHA.

A `dev` to `main` pull request is a `PROMOTION`, not an ordinary combined
architecture implementation pull request. It contains no implementation
changes. Trusted release admission evaluates source coverage, exact-SHA staging
and readiness evidence, result-tree identity, and security. It does not rerun or
re-review every implementation diff already admitted and staged on `dev`.
Legacy jobs may still execute on promotion pull requests during the staged
migration, but those reruns are not the normative promotion decision model.

A push to `main` verifies and builds the exact release SHA before deployment.
Direct `main` hotfixes are unsupported; urgent changes still enter through
`dev`, complete staging, and use the promotion boundary.

### Stable responsibility owners

| Responsibility | Canonical owner |
| --- | --- |
| Candidate fast and full test execution | `.github/workflows/test.yml` |
| Gallery staging evidence | `.github/workflows/gallery-publication.yml` |
| Trusted PR/promotion admission | `.github/workflows/web-base-policy.yml` |
| Ordinary Web policy CLI | `tools/check-web-change-budget.mjs` |
| Promotion evidence verifier | `tools/check-promotion-readiness.mjs` |
| Release verification/build/deploy | `.github/workflows/deploy_web.yml` |

`tools/check-promotion-readiness.mjs` is the planned checker implementation
owner. It does not exist yet and must be introduced in a checker-only pull
request before a later pull request wires it to a workflow or evidence producer.
The table declares ownership; it does not claim that unwritten or unwired
behavior is already executable.

After the operational protection cutover, the complete required-status set for
pull requests to `dev` is:

```text
Web base policy (trusted base)
PR / gate
```

`PR / gate` aggregates all mandatory fast candidate jobs and fails unless each
succeeds. `Web base policy (trusted base)` evaluates untrusted candidate data
with trusted base code. Legacy leaf jobs may continue to run, but they are not
part of the final required set.

## Gate and Review are independent

Every policy result has two axes:

```text
Gate: PASS | FAIL
Review: CLEAR | REQUIRED
```

`Gate` determines the executable exit status. `Review` identifies required
human attention and does not fail CI by itself. In particular, `Review:
REQUIRED` with `Gate: PASS` exits zero.

Until the checker adopts the final field names, its output maps as follows:

```text
current Result      -> final Gate
current Size review -> one source of final Review
```

Other Review sources remain documented or human-classified until their final
executable aggregation is implemented. Do not claim that the current checker
already emits the final fields.

### Gate FAIL

The following are hard failures when applicable:

- new production dependencies or bare production imports;
- binary runtime additions or vendored runtime changes;
- incomplete or cyclic first-party imports;
- unauthorized privileged owners or importers;
- malformed authority or rules;
- active deterministic architecture violations;
- prohibited runtime and guard changes in one pull request;
- checker implementation combined with authority or an evidence producer in
  one ordinary pull request;
- candidate execution in a trusted workflow;
- promotion topology, source-coverage, exact-SHA evidence, result-tree, or
  security failures; and
- gated performance or scientific-output regressions.

Documentation cannot waive or convert these failures to advisory results. A
time-bounded hard-invariant waiver requires an explicit governance decision and
the exception evidence defined by the ratchet; it does not make a failing
executable check pass. Any executable change needed to represent an approved
waiver follows the checker/authority separation below.

### Review REQUIRED

Human review is required for at least these cases:

- a selected size threshold is exceeded;
- a governance- or authority-only change;
- an architecture-bearing owner, path, or responsibility move;
- a new module, public export, reactive owner, watcher, or resource signal;
- a material performance-baseline change;
- a scientific or reference-output change; or
- a compatibility-path change.

`Review: REQUIRED` alone exits zero and is not permission to ignore a hard
failure.

## Current Web size review

Production scope covers `gbdraw/web/index.html`, `gbdraw/web/js/`, and
`gbdraw/web/vendor/`. Every pull request uses one of two review profiles:

| Profile | Production files | Gross churn | Net additions |
| --- | ---: | ---: | ---: |
| Ordinary | 8 | 800 | 100 |
| Architecture | 12 | 1,500 | 400 |

The ordinary profile is the default. Apply the `architecture-change` label when
a Web implementation adds, removes, moves, duplicates, or redefines an
architecture owner, canonical path, compatibility path, privileged boundary,
or lifecycle responsibility. The label selects the larger architecture review
profile and routes reviewer attention; it does not waive a hard failure.

Changes at or below every selected threshold report `Size review: CLEAR`. An
excess reports `Size review: REQUIRED`, but does not affect `Result` or the exit
status. A size-only excess reports `Result: PASS` and exits with status `0`.

Gross churn is the sum of textual additions and deletions in production scope.
Net additions are additions minus deletions. Binary files are counted
separately. Contributors must not game the measurements by compressing source
onto fewer lines, omitting tests or documentation, moving logic outside
measured paths, retaining dead paths, or splitting an incomplete change.

New modules, exports, `create*` declarations, reactive declarations, watchers,
compatibility-like names, session fields, and resource-like names remain in the
checker summary for review. These signals do not fail the check unless an
independent deterministic rule gives them hard semantics.

## CI diff scope and exact-SHA evidence

Pull request checks compare the pull request base SHA with its head SHA. A push
to `dev` checks only the change integrated by that push, from
`github.event.before` to `github.sha`. A manual `dev` staging run has no push
payload, so it checks the current merge commit from its first parent to `HEAD`.
These ranges determine only the change presented to the Web checker. The full
staging matrix runs from the exact current `dev` checkout.

This scope prevents independently reviewed changes already on `dev` from being
reclassified as one combined runtime-and-guard implementation. Promotion does
not weaken that separation: it accepts aggregation only from trusted exact-SHA
evidence and verified tree identity.

## Self-authorization separation

Checker implementation consists of:

```text
tools/check-web-change-budget.mjs
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
tools/web-change-source.mjs
tools/web-promotion-context.mjs
tools/check-promotion-readiness.mjs
```

Authority and evidence producers consist of:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
docs/internal/WEB_CHANGE_POLICY.md
tools/web-change-policy.json
tools/web-architecture-rules.json
tools/web-architecture-violations.json
.github/workflows/test.yml
.github/workflows/gallery-publication.yml
.github/workflows/web-base-policy.yml
.github/workflows/deploy_web.yml
```

The accepted-violation file is absent when no frozen rule has accepted debt.
Its absence does not remove it from the authority class.

The following rules apply:

- checker implementation does not change with authority or an evidence
  producer in the same ordinary pull request;
- production runtime does not change with guard or CI authority, except the
  already-defined narrow safe contraction below;
- normative authority is declared before executable wiring;
- promotion aggregation is admitted only through trusted exact-SHA evidence,
  not by same-pull-request self-approval; and
- a checker implementation may change alone, but not with the workflow or
  evidence producer that invokes or feeds it.

`tools/check-promotion-readiness.mjs` is checker implementation, not authority
data. Its first implementation therefore belongs in a checker-only pull request;
workflow wiring follows after that implementation is merged and characterized.

### Current guard scope

The ordinary Web guard additionally treats these files as a separated group:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
.github/pull_request_template.md
tools/check-web-change-budget.mjs
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
tools/web-architecture-rules.json
tools/web-architecture-violations.json
tools/web-change-source.mjs
tools/web-change-policy.json
docs/internal/WEB_CHANGE_POLICY.md
tests/web/architecture-contracts.test.mjs
tests/web/architecture-ratchet-fixtures.test.mjs
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

The current executable guard does not yet enumerate every lifecycle evidence
producer listed above. The broader normative separation governs future wiring;
the list here documents the current CLI transition without overstating it.

### Planned Product Impact separation

Product Impact extends the same trust direction without adding a workflow or
required status. The existing `tools/check-web-change-budget.mjs` remains the
ordinary Web policy CLI and trusted orchestration owner. The planned pure
mechanics and bounded decision-source modules are checker implementation:

```text
tools/web-product-impact-evaluation.mjs
tools/web-product-impact-decision-source.mjs
```

The planned Product Impact authority and evidence-producing surfaces are:

```text
docs/internal/PRODUCT_IMPACT_RATCHET.md
tools/web-product-impact-map.json
tools/web-product-decisions.json
.github/pull_request_template.md
mapped behavior contract files
```

These planned files are not executable integration. Each layer is added in the
staged order required by the Product Impact policy.

Runtime admission uses the trusted base checker, detectors, architecture rules,
Product Impact map, durable decisions, maintainer allowlist, and base/head
source facts. Candidate authority is validation-only inert data and cannot
authorize the same candidate runtime. A runtime pull request therefore cannot
change Product Impact authority or durable decisions to approve itself.

Future preauthorization may use only this narrow inert authority bundle:

```text
tools/web-architecture-rules.json
tools/web-product-impact-map.json
tools/web-product-decisions.json
```

The bundle is valid only when changed paths are a subset of those files and no
runtime, checker, workflow, template, or evidence-producing contract changes.
The trusted checker validates candidate authority but continues to evaluate
runtime with base authority.

Mapped evidence follows the same separation. An affected runtime pull request
cannot replace its mapped behavior contract and use that candidate change as
the sole proof of hard safety. Update or add the contract in a prior
evidence-only pull request, then update its authority reference separately.
Unrelated test changes are unaffected.

The runtime/guard exception applies only when `tools/web-change-policy.json` is
the sole changed guard file and the change is a pure policy contraction. The
checker implementation and source parser cannot change with that policy or its
policy workflows.

## Changing a privileged owner or importer

Expansion uses two pull requests in this order:

1. Submit a guard-only preauthorization pull request. Add the proposed owner or
   importer path to `tools/web-change-policy.json`. Do not change Web runtime.
2. Merge the preauthorization pull request. Then submit the architecture-change
   implementation against the updated base. Change runtime without changing
   policy, checker, architecture tests, or policy workflows.

Contraction may use two pull requests in the opposite order:

1. Submit the runtime pull request that removes the owner or importer use.
2. After it merges, submit a guard-only pull request that removes the stale path
   from `tools/web-change-policy.json`.

Alternatively, runtime removal and its corresponding policy contraction may be
submitted together. This narrow path is permitted only when it removes at least
one owner or importer path, adds none, preserves the exact capability and import
target key sets, and leaves every policy array a subset of its base-policy
array. `tools/web-change-policy.json` must be the only changed guard file. A
formatting-only or reordered policy change is not a contraction.

Every capability and import-target key in the base policy must remain, including
keys whose arrays become empty. The candidate runtime must remain covered by
base authority, the proposed policy must cover head runtime, and all unrelated
dependency, vendor, binary, integrity, and architecture gates remain blocking.

## Architecture evidence and authority

The architecture guard separates source detection, pure evaluation, and
intended authority:

| Concern | Owner |
| --- | --- |
| Versioned executable source-fact detection | `tools/web-architecture-detectors.mjs` |
| I/O-free schema, observation, and authority-delta mechanics | `tools/web-architecture-evaluation.mjs` |
| Intended semantic-owner and canonical-path rules | `tools/web-architecture-rules.json` |
| Permitted privileged operator and importer paths | `tools/web-change-policy.json` |
| Git and file I/O, trusted-base orchestration, reporting, and ordinary CLI entry point | `tools/check-web-change-budget.mjs` |

Detectors emit normalized facts and stable subjects. They do not contain
allowed paths, counts, enforcement modes, baseline eligibility, or exceptions.
The evaluator receives loaded plain data and normalized facts. It does not read
files, run Git, inspect the environment, emit reports, detect source facts, or
define intended paths and thresholds.

The rule registry is inert authority data. Candidate rules are parsed strictly
and checked with detectors from the trusted base. Candidate authority is
compared with untouched base source before admission and cannot authorize the
same head runtime.

## Trusted admission

The pull request workflow runs normal candidate code. The separate
`pull_request_target` workflow provides `Web base policy (trusted base)`. It
checks out the base SHA, fetches the candidate head only as Git data, and runs
the checker from the base checkout. It never checks out or executes candidate
code. A proposed head policy is read only as inert data; neither a head policy
nor a head checker can authorize itself. Proposed head detector code is never
imported or executed.

Promotion readiness must use the same trust direction: trusted code verifies
the `dev` source coverage, exact successful staging evidence, result-tree
identity, and security conditions for the exact proposed promotion head. A
promotion PR cannot manufacture its own successful evidence and cannot contain
implementation changes.

The checker reports uncertain naming signals without failing the gate. Hard
ownership checks use declarations, imports, executable owner expressions,
manifest dependencies, and repository file metadata.
