# Web change policy

The Web checker measures production change size separately from blocking
integrity and architecture violations. Production scope covers
`gbdraw/web/index.html`, `gbdraw/web/js/`, and `gbdraw/web/vendor/`. Every pull
request uses one of two review profiles:

| Profile | Production files | Gross churn | Net additions |
| --- | ---: | ---: | ---: |
| Ordinary | 8 | 800 | 100 |
| Architecture | 12 | 1,500 | 400 |

The ordinary profile is the default. The `architecture-change` label selects
the larger architecture review profile; it does not waive a blocking violation.
Changes at or below every selected threshold report `Size review: CLEAR`. An
excess reports `Size review: REQUIRED`, but does not affect `Result` or the exit
status. A size-only excess reports `Result: PASS` and exits with status `0`.

Gross churn is the sum of textual additions and deletions in production scope.
Net additions are additions minus deletions. Binary files are counted separately.
Contributors must not game the measurements by compressing source onto fewer
lines, omitting tests or documentation, moving logic outside measured paths,
retaining dead paths, or splitting an incomplete change.

Blocking violations determine `Result` and the exit status. The following checks
remain blocking:

- new production dependencies and bare production imports;
- additions of binary runtime files and changes under `gbdraw/web/vendor/`;
- prohibited combinations of Web runtime and guard or CI changes;
- checker or source-parser changes combined with authority policy or workflow
  changes;
- unauthorized privileged owners or importers, invalid allowlist changes, and
  incomplete import graphs;
- first-party static import cycles; and
- candidate, trusted-base, or active architecture-rule errors and failures.

The required status names remain `Web change budget` and
`Web base policy (trusted base)` for branch-protection compatibility. GitHub
Actions emits one bounded warning annotation when size review is required; the
step summary contains the complete measurements and reasons.

## CI diff scope

Pull request checks compare the pull request base SHA with its head SHA. A push
to `dev` checks only the change integrated by that push, from
`github.event.before` to `github.sha`. A manual `dev` staging run has no push
payload, so it checks the current merge commit from its first parent to `HEAD`.
These ranges determine only the change presented to the Web checker. The full
staging matrix continues to run from the current `dev` checkout.

This scope prevents independently reviewed changes already on `dev` from being
reclassified as one combined runtime-and-guard change. It does not add a
promotion exception or alter any blocking rule, threshold, authority source, or
checker behavior.

## Guard separation

Except for the pure policy contraction described below, Web runtime files cannot
change in the same pull request as these guard files:

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

The exception applies only when `tools/web-change-policy.json` is the sole changed
guard file. The checker implementation and source parser cannot change together
with that policy or either policy workflow. This prevents one pull request from
changing both the enforcement code and its authority data.

New modules, exports, `create*` declarations, reactive declarations, watchers,
compatibility-like names, session fields, and resource-like names remain in the
checker summary for review. These signals do not fail the check.

## Changing a privileged owner or importer

Expansion uses two pull requests in this order:

1. Submit a guard-only preauthorization pull request. Add the proposed owner or
   importer path to `tools/web-change-policy.json`. Do not change Web runtime
   files in this pull request.
2. Merge the preauthorization pull request. Then submit the architecture-change
   implementation pull request against the updated base branch. Change the Web
   runtime without changing the policy, checker, architecture test, or policy
   workflows. Apply the `architecture-change` label when the implementation needs
   the larger review profile.

Contraction may use two pull requests, in the opposite order:

1. Submit the runtime pull request that removes the owner or importer use. Do not
   edit the policy in that pull request.
2. After the runtime change merges, submit a guard-only pull request that removes
   the stale path from `tools/web-change-policy.json`.

Alternatively, the runtime removal and its corresponding policy contraction may
be submitted in one pull request. A same-pull-request contraction is permitted
only when it removes at least one owner or importer path, adds none, preserves
the exact capability and import-target key sets without additions, removals, or
renames, and leaves every policy array as a subset of its base-policy array.
`tools/web-change-policy.json` must be the only changed guard file. A
formatting-only or reordered policy change is not a contraction.

A contraction removes path entries only. Every capability and import-target key
present in the base policy must remain in the proposed policy, including keys
whose arrays are empty.

This narrow same-pull-request path is safe because:

- no authorization is added;
- the base policy remains authoritative for every runtime expansion;
- the proposed policy must still cover every privileged owner and importer active
  in the head runtime;
- the checker, source parser, architecture test, documentation, and workflows
  remain unchanged; and
- trusted-base execution remains authoritative.

An allowed contraction remains subject to the selected size-review thresholds.
A threshold excess requires review but does not fail the check. Dependency,
vendor, binary, and other integrity checks remain blocking. Any other
runtime/guard combination remains prohibited.

## Architecture evidence and authority

The architecture guard separates source detection, pure evaluation, and intended
authority:

| Concern | Owner |
| --- | --- |
| Versioned executable source-fact detection | `tools/web-architecture-detectors.mjs` |
| I/O-free schema, observation, and authority-delta mechanics | `tools/web-architecture-evaluation.mjs` |
| Intended semantic-owner and canonical-path rules | `tools/web-architecture-rules.json` |
| Permitted privileged operator and importer paths | `tools/web-change-policy.json` |
| Git and file I/O, trusted-base orchestration, reporting, and the CLI entry point | `tools/check-web-change-budget.mjs` |

Detectors emit normalized facts and stable subjects. They do not contain allowed
paths, counts, enforcement modes, baseline eligibility, or exceptions. The
evaluator receives already-loaded plain data and normalized facts. It does not
read files, run Git, inspect the environment, emit reports, detect source facts,
or define intended paths and thresholds.

The optional rule registry is inert authority data. Its absence is valid during
the staged rollout. A proposed registry is parsed strictly and checked with the
detectors already present on the trusted base. Candidate rules are compared with
untouched base source before admission and cannot authorize the same head runtime.
Until rule activation in a later phase, the checker's existing inline constraints
remain the only active runtime authority.

## Trusted base check

The pull request workflow runs the normal tests from the pull request checkout.
The separate `pull_request_target` workflow provides the required
`Web base policy (trusted base)` check. It checks out the base SHA, fetches the
pull request head only as Git data, and runs the checker from the base checkout.
It does not check out or execute pull request code. A proposed head policy is read
only as inert JSON data; neither a head policy nor a head checker can authorize
itself. Proposed head detector code is never imported or executed by this
workflow.

The checker reports uncertain naming signals without failing the check. These
signals include cache, token, handle, journal, protocol, manager, compatibility,
and generic session object key names. Hard ownership checks use declarations,
imports, executable owner expressions, manifest dependencies, and repository
file metadata.
