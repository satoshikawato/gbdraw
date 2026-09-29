# S05 protected characterization dependency

Runtime permission is present in trusted dev `d313b70b9f97c2c1d70f9ae885edbead80b62021` (PR #620).
Two required tests in `tests/web/architecture-contracts.test.mjs` still enumerate the old active Worker set:

- `Worker construction and the diagram-generation client have explicit owners`: expects three construction locations; the preauthorized import client is the fourth, with one constructor.
- `shared privileged detectors preserve the characterized current-source facts`: omits the same exact operator from the source characterization.

The file is an explicit protected guard in `tools/check-web-change-budget.mjs` and `WEB_CHANGE_POLICY.md`.
S05 does not change it, the checker, policy, authority, or mapped contracts. The unchanged mapped Node rendering and
History/admission test bodies still pass; this does not certify the mapped browser suite. Local policy Gate PASS / Review REQUIRED does not override the two test failures.
This is a guard delivery dependency, not missing Worker permission or unresolved Product choice.

[S05-guard-characterization-candidate.patch](S05-guard-characterization-candidate.patch) is an **inert** separate-delivery candidate.
It adds only the already-preauthorized client's exact count to the two inventories when that production module exists.
It preserves the old inventories when absent, rejects extra constructors, leaves allowed owners/importers unchanged,
and leaves every named mapped behavior contract unchanged. It is not applied to the implementation checkout.

Disposable read-only overlays at `/tmp/issue597-S05-resume-evidence/guard-review` exercised the two updated
inventories and the two unchanged mapped contracts against both trusted dev source and S05 candidate source:
4 PASS each. Commands, hashes and limits are recorded in the S05 resume validation manifest.
These local candidate checks are not trusted CI, external review, a waiver, or runtime admission.

Required external route: review this exact patch in a guard-only delivery based on current trusted dev, run the
required checks there, and merge it into dev before admitting the S05 runtime. No policy permission is added.
Since mapped references and mapped test bodies remain unchanged, this candidate requests no authority-reference edit;
if review changes a mapped reference, use the prescribed evidence-only then authority-ref-only sequence.
After the trusted merge, incorporate it into the existing implementation branch and rerun the unchanged guard tests
and applicable policy checks. Do not combine the active guard edit with runtime to self-authorize.

The current request authorizes only the same-named implementation branch's commit/push, not a guard branch,
PR creation/merge, dev push, or deployment. No external mutation for this candidate was performed.

## Existing mapped browser failure

The unchanged mapped browser command reported 1 PASS / 2 FAIL / 2 not run. Feature fill scope never
commits; direct preview edits wait for the same missing color application and reach their existing
600000 ms test timeout. The exact Feature fill contract also fails against a disposable static view
of trusted dev `d313b70b9f97c2c1d70f9ae885edbead80b62021` at the same undo-count assertion.
A separate observational reproduction against both views captures the same alert:
`Cannot apply feature style: Legend entry "tRNA" already exists with a different color.`
The conflicting entry is rejected by `app/legend/entry-actions.js::syncFileLegendEntries`; no rules
or History entry are committed. This behavior is present before the S05 transport change.

No color/legend outcome was selected or implemented in S05, and no mapped browser test was changed.
Resolving this existing normative behavior requires its own scoped repair and unchanged-contract
verification. The inert Worker-characterization patch addresses only the two Node inventories;
it does not repair or waive this browser failure. Raw baseline/candidate diagnostics and hashes
are in the resume validation manifest. The runner printed its failing report, while the enclosing
exec session returned 143 without a shell exit marker; both facts are preserved.

## Divergent draft contract remains unresolved

The independently run unchanged divergent-draft contract fails because Save returns `busy` with
`Updating diagram. Retry after the update finishes.` The S04 reflow coordination test passes and
requires this busy outcome; S05 does not weaken it to satisfy the older immediate-Save assumption.
Three inert test-only experiments (wait for reflow, wait for the availability owner, and assert busy
then retry once) still failed. The last experiment's mapped no-draft case passed, while divergent draft
failed and legacy was skipped; unchanged legacy separately passed. The rejected patch stays outside
the repository at `/tmp/issue597-S05-resume-evidence/rejected-draft-contract-candidate.patch`.
These are diagnostic failures, not a ready contract update or runtime authority.

Further convergence must identify queued reflow/readiness across the Save entry and preserve the
accepted busy/retry behavior and all draft/override comparisons. Any required mapped helper/body
update follows a prior evidence-only delivery, then its authority-reference delivery where required;
it cannot be applied with runtime as self-proof. This dependency is separate from the ready inert
Worker inventory patch and the existing tRNA legend conflict. No protected test, timeout, budget,
Product choice, or registry was changed here.
