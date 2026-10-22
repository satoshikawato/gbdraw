# Web change, review, and delivery policy

Status: normative repository policy

This document defines CI admission, Web architecture review, integrated `dev`
staging, release promotion, and deployment semantics. The
[architecture fitness-function ratchet](./ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
owns architecture definitions, arithmetic, and exception evidence. Workflow and
checker files implement parts of this policy; their current behavior is not a
substitute for this authority.

The [Product Impact Ratchet](./PRODUCT_IMPACT_RATCHET.md) owns
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
required path runs only the jobs the change needs, as defined in
[selective CI](./SELECTIVE_CI.md). A change to Web runtime, session, Gallery,
LOSAT, shared test, or unclassified paths also runs the full functional
Playwright suite before merge, because existing functional cases are the
checks that runtime changes break most often. Complete supported-version,
slow, browser, and Gallery evidence still belongs to the exact integrated
`dev` SHA.

A `dev` to `main` pull request is a `PROMOTION`, not an ordinary combined
architecture implementation pull request. It contains no implementation
changes. Trusted release admission evaluates source coverage, exact-SHA staging
and readiness evidence, result-tree identity, and security. It does not rerun or
re-review every implementation diff already admitted and staged on `dev`.
Legacy jobs may still execute on promotion pull requests during the staged
migration, but those reruns are not the normative promotion decision model.

A push to `main` runs browser verification for the exact release SHA. Cloudflare
Workers Builds owns production builds and deployment through the existing Git
integration; GitHub Actions does not publish the site.
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
| Release browser verification | `.github/workflows/deploy_web.yml` |
| Hosted bundle construction | `tools/prepare_cloudflare_pages.py` |
| Production deployment | Cloudflare Workers Builds and `wrangler.toml` |

`tools/check-promotion-readiness.mjs` verifies exact-SHA workflow evidence and
promotion result-tree identity. `Promotion / gate` in `web-base-policy.yml` runs
it from trusted `main`, and `Release / gate` uses it to require exact-candidate
Gallery readiness. The table declares ownership; it does not claim that
unwritten or unwired behavior is already executable.

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

### Merging pull requests into `dev`

Branch protection on `dev` requires the two checks above. It does not require
a pull request branch to be up to date with `dev`, so a pull request that
passed on an older base can merge without another CI cycle. Two rules replace
the up-to-date requirement:

- After a pull request that changes how the required checks judge a pull
  request merges into `dev`, update every other open pull request to `dev`
  with `gh pr update-branch <number>`, and merge it only after its required
  checks pass again. Both required checks run trusted code from the pull
  request base, so a result computed on the older base applied the old rules.
  These files are:
  - the checker: `tools/check-web-change-budget.mjs`,
    `tools/web-architecture-detectors.mjs`,
    `tools/web-architecture-evaluation.mjs`,
    `tools/web-owner-graph-detectors.mjs`,
    `tools/web-product-impact-evaluation.mjs`,
    `tools/web-product-impact-decision-source.mjs`,
    `tools/web-change-source.mjs`, `tools/web-promotion-context.mjs`, and
    `tools/check-promotion-readiness.mjs`;
  - the rule data that the checker reads from the base:
    `tools/web-change-policy.json`, `tools/web-architecture-rules.json`,
    `tools/web-architecture-violations.json`,
    `tools/web-product-impact-map.json`, `tools/web-product-decisions.json`,
    and `tools/web-owner-graph.json`;
  - the CI plan: `tools/ci-impact.mjs` and `tools/ci-impact-policy.mjs`;
  - the workflows of the two required checks: `.github/workflows/test.yml` and
    `.github/workflows/web-base-policy.yml`.
- A change to `tools/web-design-rule-guards.json` updates only the open pull
  requests that change a file named by an added, removed, or changed entry.
  The checker reads the registry only to classify changes to those files.
- Any other change needs no update of open pull requests. This covers a guard
  test such as `tests/web/architecture-contracts.test.mjs`, the rule text in
  `gbdraw/web/CLAUDE.md`, and the policy, contract, and template documents of
  the [current guard scope](#current-guard-scope) and
  `docs/internal/SELECTIVE_CI.md`. The checks classify these files by path and
  do not read their content from the base, and integrated `dev` staging runs
  a changed guard test on the merged SHA. If the required checks of such an
  open pull request run again later, the trusted base check compares the
  newer `dev` tree with the pull request head, so it can report those files as
  part of the pull request. Update the pull request then.
- Integrated `dev` staging is the safety net after merge. It runs on each
  exact merged SHA and finds conflicts between pull requests that passed on
  different bases.

A maintainer may queue a merge with
`gh pr merge <number> --auto --merge --match-head-commit <sha>`, where
`<sha>` is the reviewed head commit. GitHub refuses the request unless the
pull request head is exactly `<sha>`, and merges once the required checks
pass. This needs repository auto-merge to be enabled. The two rules above
still apply: update the branch first when a governance, policy, checker, or
authority change merged after the pull request's last run.

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
- binary runtime additions that the trusted base does not list, and vendored
  runtime changes that it does not list (see
  [Updating a vendored runtime file](#updating-a-vendored-runtime-file));
- incomplete or cyclic first-party imports;
- unauthorized privileged owners or importers;
- malformed authority or rules;
- active deterministic architecture violations;
- prohibited runtime and guard changes in one pull request, outside the
  static Product Contract co-change route;
- checker implementation combined with authority or an evidence producer in
  one ordinary pull request;
- a persisted contract change without its history, reference, or fixtures,
  or with a new Session version other than `main`'s or `main`'s plus one, or
  one that cannot be read (see
  [Persisted contract co-change](#persisted-contract-co-change));
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
- a static Product Contract change, including a co-change with its runtime;
- an architecture-bearing owner, path, or responsibility move;
- a new module, public export, reactive owner, watcher, or resource signal;
- a material performance-baseline change;
- a scientific or reference-output change;
- a compatibility-path change; or
- a persisted contract change, or a new generated contract export (see
  [Persisted contract co-change](#persisted-contract-co-change)).

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

Exact-SHA evidence is a successful `Dev staging / gate` or
`Gallery readiness / gate` aggregate from the push run for that SHA. The run may
earn part of it through light-change inheritance: when the push leaves the tree
unchanged, or changes only metadata, documentation, policy documentation, and
leaf tests, and the direct parent has its own successful aggregate, the run
executes only the jobs the change needs and inherits the rest
([selective CI](./SELECTIVE_CI.md#light-change-inheritance)). Promotion
checklist evidence outside these aggregates, such as the release-tier run,
recipe and sweep runs, and generated-artifact checks, may come from an ancestor
of the promotion head under the
[carry-forward rule](./WEB_PERIODIC_AUDIT.md#carrying-evidence-forward). Neither
rule changes `Promotion / gate`, which still requires both aggregates on the
exact promotion head.

## Self-authorization separation

Checker implementation consists of:

```text
tools/check-web-change-budget.mjs
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
tools/web-owner-graph-detectors.mjs
tools/web-twin-detectors.mjs
tools/web-change-source.mjs
tools/web-promotion-context.mjs
tools/check-promotion-readiness.mjs
```

Authority and evidence producers consist of:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
docs/internal/WEB_CHANGE_POLICY.md
gbdraw/web/CLAUDE.md
tools/web-change-policy.json
tools/web-architecture-rules.json
tools/web-architecture-violations.json
tools/web-design-rule-guards.json
tools/web-owner-graph.json
.github/workflows/test.yml
.github/workflows/gallery-publication.yml
.github/workflows/web-base-policy.yml
.github/workflows/deploy_web.yml
```

`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` is also authority; see
[Static Product Contract co-change](#static-product-contract-co-change).
The accepted-violation file is absent when no frozen rule has accepted debt.
Its absence does not remove it from the authority class.

The following rules apply:

- checker implementation does not change with authority or an evidence
  producer in the same ordinary pull request;
- production runtime does not change with guard or CI authority, except the
  already-defined narrow safe contraction below and the static Product
  Contract co-change route;
- normative authority is declared before executable wiring;
- promotion aggregation is admitted only through trusted exact-SHA evidence,
  not by same-pull-request self-approval; and
- a checker implementation may change alone, but not with the workflow or
  evidence producer that invokes or feeds it.

`tools/check-promotion-readiness.mjs` is checker implementation, not authority
data. It changes in checker-only pull requests, separately from the workflows
that invoke it.

### Current guard scope

The ordinary Web guard additionally treats these files as a separated group:

```text
docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
docs/internal/PRODUCT_IMPACT_RATCHET.md
docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
.github/pull_request_template.md
tools/check-web-change-budget.mjs
tools/web-architecture-detectors.mjs
tools/web-architecture-evaluation.mjs
tools/web-product-impact-evaluation.mjs
tools/web-product-impact-decision-source.mjs
tools/web-architecture-rules.json
tools/web-architecture-violations.json
tools/web-product-impact-map.json
tools/web-product-decisions.json
tools/web-change-source.mjs
tools/web-promotion-context.mjs
tools/check-promotion-readiness.mjs
tools/web-change-policy.json
tools/web-design-rule-guards.json
tools/web-owner-graph.json
tools/web-owner-graph-detectors.mjs
tools/web-twin-detectors.mjs
gbdraw/web/CLAUDE.md
docs/internal/WEB_CHANGE_POLICY.md
tests/web/architecture-contracts.test.mjs
tests/web/architecture-ratchet-fixtures.test.mjs
tests/web/product-impact-ratchet-fixtures.test.mjs
tests/web/promotion-readiness.test.mjs
tests/web/web-promotion-context.test.mjs
.github/workflows/gallery-publication.yml
.github/workflows/deploy_web.yml
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

Mapped Product Impact contracts are classified from the trusted base map when
an affected concern is evaluated. They do not need a duplicated static path
list in this policy.

### Design-rule co-change

`gbdraw/web/CLAUDE.md` holds the numbered Web design rules (R1, R2, …). Each rule names the
guard tests that enforce it, and some guards carry an allowlist or a baseline
literal that the rule's text authorizes. `tools/web-design-rule-guards.json`
(schema 2) registers, per rule, the heading, the guards, and those literals.
The checker reads the registry from the trusted base.

- A rule has an `id` and the `heading` of its text. A rule whose text is not
  one of the numbered rules names its file in its own `rulesPath` and takes an
  upper-case id: `TYPE-DEBT` is "Type checking" in
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`, and
  `CW` is "Computation ownership" in `gbdraw/web/CLAUDE.md`.
- A guard is `{ path, rootCause, job, retireWhen }`: the root-cause pattern it
  exists to stop (below, or `new`), the job of `.github/workflows/test.yml`
  that runs it, and the condition under which it is deleted (`never` for a
  behavior test with no planned replacement).
- A literal is `{ path, symbol, kind }`. The kinds are `count` (one number),
  `count-map` (a number per key), `set` (exempt names), `writer-map` (exempt
  writers per module, each with a category), `set-map` (names per path),
  `json-entries` (the objects of a JSON array, by `id`), `json-ids` (the `id`
  set of a JSON array and the length of each entry's lists: a new id or a
  further list item is an expansion, a removed id or list item a contraction,
  and a renamed list item or an edit of another field neither), and
  `py-count`, `py-count-map`, and `py-set` for a module-level literal of a
  Python test.
- A count is a ceiling. Its test asserts observed ≤ recorded, per key for a
  map, and a key without an entry counts as 0. Lowering or removing an entry
  is a contraction and is optional; the test names what may be lowered.
  Raising an entry or adding a key is an expansion. The two R6 site counts
  (`UNCLASSIFIED_THROW_BASELINE`, `UNPROVEN_BASELINE`) are one total each
  (`count`, `py-count`), so a site that moves to another module with its code
  changes nothing. The three maps keyed by module path
  (`EXPLICIT_ANY_BASELINE`, `OUTSIDE_OWNER_BASELINE`, `TYPE_DEBT_BASELINE`)
  stay per module, so a new or moved module starts at 0.

Root-cause patterns (boundary audit of 2026-10-10):

| Pattern | What recurs |
|---|---|
| `P1` | A live edit and Generate compute the same derived value twice. |
| `P2` | A setting or an edit is not confined to its diagram mode. |
| `P3` | A migration or a reader exists on one side only (JavaScript or Python). |
| `P4` | JavaScript derives again what Python decided. |
| `P5` | An identity is read from a displayed value (a caption, a swatch color). |
| `P6` | The meaning of a failure is lost between Python and the Web. |
| `P7` | State is written outside a transition, or depends on timing. |
| `P8` | One surface writes a file that another surface reads differently. |
| `P9` | An asynchronous or cached copy is read as the truth. |
| `P10` | A guard has a blind spot of its own. |

The Gate applies these rules:

- A diff that changes `gbdraw/web/CLAUDE.md` or the registry and any
  production runtime path fails the Gate ("production runtime files and Web
  guard/CI files changed together"). Rule text and registration change in an
  authority-only pull request, before or after the runtime change they
  describe.
- A diff that changes production runtime paths may contract a registered
  allowlist in the same guard test (remove an entry, lower a count) and passes;
  the report lists the contraction.
- A test that holds a registered count literal (`count`, `count-map`,
  `py-count`, `py-count-map`) counts with its own code, so beside production
  runtime paths (for a `py-*` literal, also beside Python package files) it
  changes only inside its registered literals, which stay readable, except
  that a quoted path may follow a rename in the same diff; the counting code
  those tests share (`COUNTING_CODE_PATHS` in
  `tools/check-web-change-budget.mjs`) does not change beside production
  runtime paths at all. Any other change, deletion included, fails the Gate
  with `design-rule.count-scope` and goes to a pull request without
  production code.
- A diff that expands a registered allowlist, or leaves it unreadable, alongside
  production runtime paths fails the Gate with `design-rule.co-change`. For a
  `py-*` literal, a Python file of the package (`gbdraw/**/*.py` outside
  `gbdraw/web/`) counts as production as well. The expansion goes to its own
  authority-only pull request, which is Review REQUIRED ("registered
  design-rule allowlists expanded") and names the rule and the reason for the
  new entry. A new twin row in `tools/twins.json` (R9) cannot follow its twin,
  because the twin guard fails a twin without a row, so it goes in first with
  `"pending": true`, and the runtime pull request that adds the twin removes
  `pending`.
- A malformed registry fails closed and is Review REQUIRED.

Adding a test case to a guard test is neither an expansion nor a contraction
and is not restricted by this rule, except in a count-literal test or the
counting code: there a new case, a failing-first case included, goes in a
pull request without production code. A new allowlist or baseline literal in a
guard test is registered in the same authority-only pull request that adds it.

`tests/web/architecture-contracts.test.mjs` checks the registry in both
directions: every `tests/` path that `gbdraw/web/CLAUDE.md` names is
registered, every registered path is named by its rule, every guard is
selected by the job it declares, and every node test under `tests/web` is
selected by a pull-request job. It reads the selection from the workflow's
commands, `tests/ci/functional-shards.json`, and the pytest markers.

### Product Impact separation

Product Impact extends the same trust direction without adding a workflow or
required status. The existing `tools/check-web-change-budget.mjs` remains the
ordinary Web policy CLI and trusted orchestration owner. The pure
mechanics and bounded decision-source modules are checker implementation:

```text
tools/web-product-impact-evaluation.mjs
tools/web-product-impact-decision-source.mjs
```

The Product Impact authority and evidence-producing surfaces are:

```text
docs/internal/PRODUCT_IMPACT_RATCHET.md
tools/web-product-impact-map.json
tools/web-product-decisions.json
.github/pull_request_template.md
mapped behavior contract files
```

The static Product authority path for broad or currently unmapped durable
Option Integrity outcomes is:

```text
docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
```

The file exists. It was created through the bootstrap order below, and the
checker recognizes its exact path as both guard and authority. It is an inert
static Product authority surface. It adds no workflow, required status,
evaluator, or machine-readable decision store, and no parser reads its prose.

The bootstrap order that created it was:

```text
normative process/path preauthorization
  -> checker-only exact-path enforcement
      -> explicit human Product Decision receipt
          -> authority-only Product Contract creation
              -> dependent runtime implementation
```

#### Static Product Contract co-change

A pull request may change the static Product Contract together with the
runtime, tests, and documentation that implement the changed outcome. All of
these conditions must hold:

- the Contract is modified in place at its exact path. It is not added,
  deleted, renamed, or copied, and a similarly named path is not the Contract;
- the Contract is the only changed guard or authority file. Checker
  implementation, detectors, evaluators, workflows, the Product Impact map,
  durable decisions, architecture rules, accepted violations,
  `tools/web-change-policy.json`, policy documents, the pull request template,
  and guard tests remain unchanged;
- no mapped behavior contract of an affected concern is replaced as the sole
  hard proof (see mapped evidence below); and
- each changed Contract record serializes the Product Decision Owner's
  explicit receipt: choice, rationale, preservation, retirement scope,
  accepted residual risk, owner, and decision date. Missing receipt fields
  leave the change unresolved, and nobody fills them by inference.

The co-change is `Review: REQUIRED`. The machine Gate classifies paths and
evaluates mapped concerns with trusted base authority only. Candidate Contract
text never becomes machine authority. The human reviewer confirms that the
serialized receipt matches the explicit human choice and that the runtime and
tests realize exactly that outcome. The Gate does not claim to verify Product
meaning that no parser reads. Headings that trusted authority references,
such as `OIPC-C07` and `PD-OI-016`, keep their anchors.

Contract-only, runtime-only, and mapped-decision pull requests keep their
existing routes. The route does not extend to any other guard or authority
file.

The route takes effect for a pull request only when that pull request's base
checker implements it. Until then, a Contract change must be isolated from
every other changed path. Once the route is effective, it replaces only one
procedural rule: that Contract authority must merge before dependent runtime.
That rule appears in the Contract lifecycle and in individual records. It does
not change any Product outcome, preservation condition, or compatibility
commitment. The first permitted co-change aligns those procedural Contract
sentences.

The JSON files are inert authority. The checker loads trusted base authority,
validates candidate authority as data, and evaluates affected concerns without
executing candidate modules.

Runtime admission uses the trusted base checker, detectors, architecture rules,
Product Impact map, durable decisions, maintainer allowlist, and base/head
source facts. Candidate machine authority is validation-only inert data and
cannot authorize the same candidate runtime. A runtime pull request therefore
cannot change the Product Impact map or durable decisions to approve itself.
A static Product Contract co-change does not approve itself either. Its
authorization is the reviewed explicit human receipt, not the candidate text.

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

Apart from the static Product Contract co-change, the runtime/guard exception
applies only when `tools/web-change-policy.json` is the sole changed guard file
and the change is a pure policy contraction. The
checker implementation and source parser cannot change with that policy or its
policy workflows.

## Updating a vendored runtime file

A change under `gbdraw/web/vendor/` passes the Gate only when the
`allowedVendorChanges` list in `tools/web-change-policy.json` on the pull
request base names it exactly. The checker reads the list from the base. A
list in the pull request head is checked as data and authorizes nothing.

Each entry names one file:

| Field | Value |
| --- | --- |
| `path` | the exact path under `gbdraw/web/vendor/` |
| `change` | `add`, `replace`, or `remove` (Git status `A`, `M`, or `D`) |
| `sha256` | SHA-256 of the bytes Git stores for the file after the change (`add`, `replace`) or before it (`remove`) |
| `package` | the upstream package name |
| `version` | the upstream version |
| `source` | the `https` URL the bytes come from |

A path appears at most once. A change whose path, kind, or bytes differ from
its entry fails as an unlisted change. A malformed list fails closed, in the
base and in the candidate, and so does a policy file that repeats a key. A new
entry that cannot apply to the base tree fails in the pull request that adds
it: `add` for a path the base has as a file or directory, `replace` or `remove`
for a path it lacks,
`remove` whose hash differs from the base bytes, and `replace` with the bytes
the base already has. The report lists every added and removed entry. A listed
change passes with `Review: REQUIRED`. An added binary file under
`gbdraw/web/vendor/` passes the binary rule only when the trusted base lists
it as `add` and the head bytes hash to the entry's `sha256`; any other added
binary runtime file still fails ("added binary runtime files are not
allowed"). The route takes effect for a pull request only when that pull
request's base checker implements it; a checker-only pull request follows the
one that lists the files.

An update takes two pull requests in this order:

1. Authority only: add the entry to `tools/web-change-policy.json`. Take the
   bytes from the upstream release, check them against the registry's
   published integrity value, and hash them as Git stores them
   (`git cat-file blob <commit>:<path> | sha256sum`; with `eol=lf` this equals
   the upstream file when it has LF line endings). The pull request names the
   reason, such as a security fix. It expands authority, so it waits for the
   Owner's approval.
2. After it merges: the runtime pull request changes the file and the code,
   tests, notices, and packaging that use it, and leaves the list unchanged.
   Changing the list beside a runtime path fails the Gate.

The check compares the pull request base with its head, so an open pull
request that is behind `dev` after step 2 merged shows the listed change in
reverse. When the head does not contain its base and leaves a listed path as
their merge base had it, the checker reports that the branch is behind its
base; update it from `dev`. Any other mismatch reports the kind or the hash.

An entry whose change has merged can only authorize the same bytes again.
Remove it in the next authority-only pull request. Because a path appears once,
the next update of that path replaces its entry.

## Persisted contract co-change

The generated modules in `gbdraw/web/js/contracts/` describe the persisted
Session and request contracts; `python tools/generate_web_contracts.py --check`
keeps them equal to their Python owners. The Gate compares them between the
pull request base and head, without running Python:

- `schema.generated.js`: each `export const` statement, ignoring owner
  comments, JSDoc, and whitespace. A changed or removed export counts unless
  `exemptGeneratedExports` names it. A removed export whose exact value
  reappears under a new name is a rename and needs review only; renaming a
  version export (`SESSION_VERSION`, `SUPPORTED_SESSION_VERSIONS`,
  `CANONICAL_REQUEST_SCHEMA`, `SUPPORTED_CANONICAL_REQUEST_SCHEMAS`) fails,
  because the checker reads them by name and needs a change first.
- `request-shape.generated.js` and `session-shape.generated.js`: the
  description as data, with type references inlined, so renaming or
  reordering types does not count.
- `request.generated.js` and `session.generated.js` (JSDoc only) never count.
- A module or version value that cannot be read counts as changed and fails
  closed. A version export that the base has and the head lacks or cannot
  read fails as "a renamed or removed version export (<name>) needs a checker
  change".
- With `--head`, a branch behind its base is compared from the merge base of
  base and head, so only its own changes count, and only its own diff supplies
  the co-changes. When Git cannot find a merge base, the comparison falls back
  to the base, and a failure may come from a newer `dev` change; the message
  then says to update the branch from `dev`. Without `--head` (the working
  tree), the comparison is always two-dot against the base.

| Kind | Counts when | Also needs |
| --- | --- | --- |
| `session-version` | `SESSION_VERSION` or `SUPPORTED_SESSION_VERSIONS` changes | for a new version, a `## ` heading in the history naming `Session <n>`; while the old version stays supported, a `tests/fixtures/sessions/*.v<old>.*` fixture |
| `request-schema` | `CANONICAL_REQUEST_SCHEMA` or `SUPPORTED_CANONICAL_REQUEST_SCHEMAS` changes | for a new schema, a `## ` heading naming `request schema <n>` |
| `request-description` | the request description changes and the request schema does not | nothing more |
| `session-admission` | the Session admission description changes | nothing more |
| `schema-constant` | any other export changes | nothing more |

Every kind needs `docs/internal/SESSION_FORMAT_HISTORY.md`,
`docs/REFERENCE/session-and-request-compatibility.md`, and a changed path under
`tests/fixtures/` in the diff. A missing item fails the Gate
(`persisted-contract.co-change`), naming the kind, the exports or modules, and
each missing item. When all are present, the Gate passes with
`Review: REQUIRED`. A new export or a new description module is
`Review: REQUIRED` only: confirm that it projects an existing value.

The rule data is `persistedContractCoChange.exemptGeneratedExports` in
`tools/web-change-policy.json`, read from the base. Without the key the rule is
off. Malformed data fails closed in the base and in the candidate. Neither may
exempt a version export, and a candidate may not exempt a name the head does
not export. Adding an exemption, or turning the rule on or off, is
`Review: REQUIRED`; adding an exemption or removing the key expands authority
and waits for the Owner's approval. The two documents and `tests/fixtures/`
are evidence, not guard or authority paths.

A breaking Session change raises `SESSION_VERSION` when it enters `main`; the
unit is `main`, not a release tag. A change is breaking when the `main` build
would refuse a Session that the new writer writes, or would read it
differently. While the head writes a version that `main` does not, a further
breaking change keeps that version: rewrite its history entry and fixtures, and
remove any reader or migrator of a version that never reached `main`. With
`persistedContractCoChange.sessionVersionUnit` set to `"main"` in the base
policy (another value fails closed in the base and the candidate), the Gate
reads `main`'s `SESSION_VERSION` from `refs/remotes/origin/main` for every pull
request with one of the kinds above:

- a pull request that changes `SESSION_VERSION` to anything other than
  `main`'s or `main`'s plus one fails the Gate (`persisted-contract.fold`);
- when the head writes `main`'s version, `Review: REQUIRED` asks to confirm
  that `main` reads what the head writes the same way, or to raise
  `SESSION_VERSION` by one;
- when either version cannot be read, the Gate fails closed
  (`persisted-contract.fold`); fetch `main` and rerun.

Whether a change breaks `main` is the reviewer's judgment.
`CANONICAL_REQUEST_SCHEMA` has no rule tied to `main`; it keeps the rules in
the table. Setting or removing `sessionVersionUnit` is
`Review: REQUIRED`; removing it expands authority and waits for the Owner's
approval. `Promotion / gate` does not repeat the check: every change reaches
`main` through a `dev` pull request.

In the pull request body, under "Risk and review notes", write
"Persisted contract:" and either `N/A` or the kinds the Gate reported, the
history heading, the compatibility reference section, the fixture paths, and,
for an older Session version that stays readable, the first-parent `main`
commit or release tag that wrote its fixture. The Gate does not check reader
behaviour outside the descriptions or that `main`/tag evidence; the reviewer
does. Before opening the pull request, commit, then run
`python tools/generate_web_contracts.py --check` and
`node tools/check-web-change-budget.mjs --base origin/dev --head HEAD`.

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
The only exception is a retired capability key, defined below.

### Splitting a capability

A capability key is split, or renamed, in three pull requests in this order:

1. Authority only: add each new key to `allowedPrivilegedOwners` with every
   owner of the key it replaces. Do not change runtime, detectors, or the
   checker.
2. Checker only: after step 1 merges, replace the old detector entry in
   `tools/web-architecture-detectors.mjs` with one entry per new key and
   update the characterization in `tests/web/architecture-contracts.test.mjs`.
   The old policy key is then retired: no trusted detector reports it, so it
   authorizes nothing.
3. Authority only: after step 2 merges, remove the retired key and reduce each
   new key to the owners that its detector reports.

A retired key may be removed only when the trusted base detector no longer
defines it. Removing a key that the base detector still defines remains a
blocking failure.

The former `Mounted SVG/Result replacement` capability was split this way into:

- `Result content commit`: code that replaces the current Result, that is an
  assignment to `results.value` or `state.results.value`, or a call to
  `flushActiveResult(`; and
- `SVG serialization`: code that serializes a mounted SVG with
  `serializeCleanSvg(`.

Web design rule R1 in `gbdraw/web/CLAUDE.md` limits which code may write the
Result, so the `Result content commit` list only contracts. Do not preauthorize
a new owner for it; route the write through an existing owner instead.

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
