# Web change policy

The Web change checks limit growth in `gbdraw/web/index.html`, `gbdraw/web/js/`,
and `gbdraw/web/vendor/`. The `architecture-change` label can waive the ordinary
file, line, module, dependency, binary file, and vendor budgets. It cannot waive
the privileged capability allowlists or the rule against changing runtime and
guard files in one pull request.

## Changing a privileged owner or importer

Use two pull requests in this order:

1. Submit a guard-only preauthorization pull request. Add the proposed owner or
   importer path to `tools/web-change-policy.json`. Do not change Web runtime
   files in this pull request.
2. Merge the preauthorization pull request. Then submit the architecture-change
   implementation pull request against the updated base branch. Change the Web
   runtime without changing the policy, checker, architecture test, or policy
   workflows. Apply the `architecture-change` label when the implementation also
   needs an ordinary budget exception.

The allowlists are upper bounds. Removing an owner or importer from runtime does
not require a policy edit. Adding an owner or importer fails until its path is in
the base branch allowlist.

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
