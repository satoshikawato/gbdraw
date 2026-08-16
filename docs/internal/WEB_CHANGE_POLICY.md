# Web change policy

The Web change checks limit growth in `gbdraw/web/index.html`, `gbdraw/web/js/`,
and `gbdraw/web/vendor/`. Every pull request uses one of two finite profiles:

| Profile | Production files | Gross churn | Net additions |
| --- | ---: | ---: | ---: |
| Ordinary | 8 | 800 | 100 |
| Architecture | 12 | 1,500 | 400 |

The ordinary profile is the default. The `architecture-change` label selects
the architecture profile; it does not waive a failed limit. Exact-limit changes
pass. Gross churn is the sum of textual additions and deletions in production
scope. Net additions are additions minus deletions. Binary files are counted
separately.

Both profiles reject new production dependencies, including new bare production
imports. They also reject every change under `gbdraw/web/vendor/` and every added
binary file in production scope.

## Guard separation

Web runtime files cannot change in the same pull request as these guard files:

```text
tools/check-web-change-budget.mjs
tools/web-change-source.mjs
tools/web-change-policy.json
docs/internal/WEB_CHANGE_POLICY.md
tests/web/architecture-contracts.test.mjs
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

The checker implementation and source parser also cannot change together with
`tools/web-change-policy.json` or either policy workflow. This prevents one pull
request from changing both the enforcement code and its authority data.

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
   the larger finite profile.

Contraction also uses two pull requests, in the opposite order:

1. Submit the runtime pull request that removes the owner or importer use. Do not
   edit the policy in that pull request.
2. After the runtime change merges, submit a guard-only pull request that removes
   the stale path from `tools/web-change-policy.json`.

The base policy remains authoritative for expansion. A proposed policy change is
also checked against the head runtime, so it cannot remove a path that still owns
or imports a privileged capability. Runtime and guard changes remain prohibited
in the same pull request.

A contraction removes path entries only. Every capability and import-target key
present in the base policy must remain in the proposed policy, including keys
whose arrays are empty.

## Trusted base check

The pull request workflow runs the normal tests from the pull request checkout.
The separate `pull_request_target` workflow provides the required
`Web base policy (trusted base)` check. It checks out the base SHA, fetches the
pull request head only as Git data, and runs the checker from the base checkout.
It does not check out or execute pull request code.

The checker reports uncertain naming signals without failing the check. These
signals include cache, token, handle, journal, protocol, manager, compatibility,
and generic session object key names. Hard ownership checks use declarations,
imports, executable owner expressions, manifest dependencies, and repository
file metadata.
