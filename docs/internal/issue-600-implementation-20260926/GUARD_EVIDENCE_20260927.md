# Issue 600 — separated architecture guard evidence

This change contains architecture tests and this report only. It is prepared
before the annotation/style runtime is admitted to dev. Runtime, checker,
source parser, workflows, Product authority, allowlists, dependencies and
public/reference artifacts are unchanged.

## Endpoints and authority

- Actual starting origin/dev: `9d967f1c72f730b205420d62d383b73127c2a9b1`.
- Existing runtime branch: `fix/issue-600-annotations-styles-20260926`.
- Fixed S05 runtime endpoint: `9c6a15288cc40d2ddd8496b2b32e2b858e11c329`.
- Evidence-only branch: `fix/issue-600-guard-evidence-20260927`, derived from dev
  with no runtime merge and initially no upstream.
- Local evidence: `/tmp/gbdraw-issue600-guard-integration-evidence` (`E`).

The actual base audit matches PD-OI-040–043 Choice A: all nine receipt fields,
five outcome clauses, human receipts, base references and uniqueness checks.
Revision 21 / 45 active IDs / 23 receipts are observations, not fixed counts.
No new Product decision or approval rationale is introduced.

## Meaning of the two corrections

The old metadata assertion treated every dependency on the pure metadata
module as rendered-identity access. S02 reuses `validateAnnotationWarnings`
from request preparation and Session authority. The corrected test keeps
exactly the two existing identity/admission callers and permits those two
specific additional callers only a single static named validator import.
Identity helpers, mixed bindings, namespace/default imports, dynamic imports,
duplicate imports, unapproved callers and loss of an admission caller fail.
There is no new production validator or mutable owner.

S03 removed an unnecessary rollback Result assignment after making canonical
legend admission atomic. Five exact replacement occurrences became four.
The existing detector and its operator grammar are unchanged. The test keeps
the required legend operator location and the existing maximum of five;
contraction is permitted, growth and missing ownership fail. All other current
source facts remain exact. Fixtures pin detection of four and five actual
operators, reject six and zero, and verify that comments/string text are not
counted. Allowlist and canonical Result admission tests remain unchanged.

This follows the architecture ratchet: a count alone is not an architecture
fitness function. The intended invariant is retained ownership without
privileged operator growth. It is not a permanent choice between two runtime
versions or a new compatibility path. Runtime semantic owners, canonical
paths and compatibility readers are unchanged; OE/PE/CB deltas are zero.

## Execution and separation

Executed in the owned isolated clone with Node 26.8.2:

```bash
node --test tests/web/architecture-contracts.test.mjs \
  tests/web/architecture-ratchet-fixtures.test.mjs \
  tests/web/product-impact-ratchet-fixtures.test.mjs
node --test tests/ci/*.test.mjs
```

- Latest dev plus this evidence: exit 0, 192 passed (`E/guard-base.log`).
- Fixed S05 plus this evidence: exit 0, 192 passed (`E/guard-runtime-probe.log`).
- Actual-base working-tree checker: exit 0, Gate PASS / Review CLEAR,
  no blocking violations (`E/guard-working-gate.log`); diagnostic only.
- CI routing requires full PR coverage: `web-change-budget`, `core-pr`,
  `recipes-standard`, `gallery`, `lint`, `web-contracts-pr`, `web-pr-smoke`.
  Remote required checks are awaited before merge; no unexecuted check is
  labeled PASS in this committed report.

An exploratory selective-plan construction was rejected by trusted policy
with `INVALID_SELECTIVE_PLAN`; full coverage is retained. No routing policy
was changed. CI contract result is recorded in `E/ci-contracts.log`.

WSL bubblewrap cannot start on the mounted app-server socket. The same local
operations ran with required sandbox escalation; no global settings changed.

The same test file is overlaid onto a separate detached worktree at the fixed
S05 SHA. That probe changes no runtime and is supplemental evidence, not
integrated authority or an S05 completion declaration. Latest dev and the
fixed runtime are independently exercised.

The policy guard requires runtime and its evidence producer to be separated.
This PR changes neither checker implementation nor Product authority. The
trusted checker and dependencies are extracted from the actual base into
`E/tools`. The final explicit full guard commit SHA is checked with those tools
and recorded externally, avoiding a self-reference in this report.

Only the two reviewed paths are staged. Tests, report and unchanged generated
paths are audited separately, with `git diff --check`. No reference outputs,
Gallery assets, social preview, wheel, binary or dependency changes belong to
this change. Shared checkout changes are preserved.

## Integration and rollback

The user requested integration of this guard/evidence prerequisite. After
passing required PR checks, merge this separate PR normally into dev. Then
fetch the resulting trusted base, incorporate it normally into the existing
#600 work branch, and rerun the affected gates. Preserve all S00–S05 history.
This guard change alone does not integrate the #600 runtime into dev and does
not assert S05 all-gates completion.

S05 runtime behavior evidence may be reused only where source, inputs,
environment and acceptance remain unchanged. Latest dev also contains #614
alignment accessibility changes; the final runtime continuation must test the
affected UI paths without claiming unchanged-source evidence for them.

Rollback is a normal revert of this evidence-only commit via a work-branch PR.
No reset, force push, direct dev/main write, publication, deployment or tag is
part of this task. Required check and branch-protection settings remain intact.

English commit title: `test: preserve metadata and legend ownership checks`.
Summary: Distinguish pure annotation warning callers from metadata admission
and permit legend operator contraction while rejecting growth and lost owners.
