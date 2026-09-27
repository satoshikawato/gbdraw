# Issue #597: Session import Worker characterization delivery

The guard implementation and local validation are complete. This delivery updates two exact Worker inventories so the preauthorized `services/session-import-client.js` constructor is counted when its production module exists. When the module is absent, the previous inventory remains required. S05 runtime acceptance remains incomplete. This delivery does not merge the PR or modify the S05 runtime branch.

Dedicated checkout: `/tmp/gbdraw-issue597-session-import-guard-20260927`.
Branch: `fix/issue-597-session-import-guard-20260927`, created without an upstream from actual remote dev/start HEAD `27939faebd4728c5aa52c9524d14ad0c3244190b`.
The actual S05 remote/checkpoint remains `c3ae518b914024d1299c1f3f6487640c8a2fbaaf`; its checkout is `/tmp/gbdraw-issue597-S05.ujOWyl`.
The proposed guard branch did not exist locally or remotely. Only the advanced dev ref was fetched. The shared checkout and all unrelated work were preserved.

The saved dev `252986d096011fcf1a0f5564e940480d3b92844d` was not treated as latest. Current dev includes PR #622's independent Web changes. The existing constructor and codec-import permissions from PR #620 remain present. The saved inert patch applies to current dev without a content change. [Fingerprints](fingerprints.json) record base/S05 source files, candidate SHA-256, checker/policy/authority bytes, mapped bodies, and existing S05 evidence.

## Changed scope

Only these two inventories in `tests/web/architecture-contracts.test.mjs` change:

- `Worker construction and the diagram-generation client have explicit owners`: conditionally includes the exact import client path, with constructor count **1**.
- `shared privileged detectors preserve the characterized current-source facts`: conditionally includes the same path and count in the detector inventory.

The conditions check production module presence. They do not derive expected counts from observations. An extra constructor, another owner path, or a missing constructor in the present module still fails both inventories.

The complete updated test file equals the saved inert patch applied to current base. The named mapped bodies `production rendering crosses the canonical request and Worker boundary` and `History intent and SVG admission have one production ownership path` are byte-identical to base. Their authority references are unchanged. No mapped-reference delivery is required by this patch; any future required reference change must first complete evidence-only then authority-ref-only delivery before dependent runtime admission.

Runtime, checker, permission, policy, Product receipt/authority, timeout, performance budget, and Worker-free assertions are unchanged. Semantic owners, canonical paths, and compatibility paths are unchanged; OE/PE/CB have no delta and no superseded runtime path exists. This is a guard characterization correction against existing authority, with no reachable user-visible effect or new Product choice. Production, tests, added evidence, and generated diffs were reviewed separately.

## Local verification

[Validation](validation.json) records exact commands, working directories, exits, summaries, original raw log paths and SHA-256 digests. Repository copies of negative assertion logs strip trailing whitespace; both original and stored hashes are recorded.

| Check | Result |
| --- | --- |
| Before implementation: base policy checker | exit 0; Gate PASS / Review CLEAR |
| Unmodified trusted dev: two inventories and two mapped bodies | 4 PASS |
| Updated guard checkout: complete architecture suite | 139 PASS, including existing privileged-expansion and self-authorization negatives |
| Architecture/Product/promotion/context fixtures | 78 PASS |
| Updated dev: two inventories and two mapped bodies | 4 PASS; old inventory preserved |
| Disposable exact S05 source overlay: same four tests | 4 PASS; exact import client count 1 |
| S05 overlay: extra constructor, extra owner, missing constructor | Each mutation: 2 expected FAIL, exit 1; both inventories reject it |
| Dev overlay: similarly named extra owner while exact module absent | 2 expected FAIL, exit 1 |
| CI contracts, corrected dependency installation | 65 PASS |
| PR Web Node contracts, excluding architecture suite | 851 PASS |
| Vibrio sequence source coverage | 3 PASS |
| Hepatoplasmataceae and tobacco Session projections | Both exit 0 |

The first CI-contract invocation had 63 PASS / 1 FAIL because the new checkout lacked `@playwright/test`. Its log remains recorded. `npm ci --ignore-scripts` installed the locked development dependencies without changing manifests; the corrected invocation passed. No assertion was relaxed.

The S05 overlay was constructed independently from committed S05 Web JS/index source, current trusted-base tools/tests/workflows, and the patched architecture test. It did not write to the implementation checkout or prior evidence. Prior guard PASS evidence was not reused because dev source advanced. Prior S05 browser/performance failures are retained as historical observations, not current PASS claims.

The implementation-after and postcommit checker invocations use the unchanged base implementation: `node tools/check-web-change-budget.mjs --base origin/dev`, then `node tools/check-web-change-budget.mjs --base origin/dev --head HEAD`. Their logs/results, final commit, remote equality, clean tree, PR URL and observed CI state are captured in `/tmp/issue597-session-import-guard-20260927-evidence/handoff-final.json` after the single commit. This separate handoff avoids a second commit or a self-referential commit hash. The implementation-after checker actually reports Gate PASS / Review CLEAR. Its unchanged implementation adds no executable review reason for this test-only guard edit; normative governance human review remains required. Neither result is trusted CI PASS or a review waiver. The first recording wrapper incorrectly expected Review REQUIRED and failed its own assertion; the checker exited 0, and its original output is retained.

The trusted base classifies the guard test as `ci-only`, so the PR requires the full candidate job set: Web budget/architecture, Core PR, standard recipes, Gallery, lint, Web contracts and PR smoke; `PR / gate` aggregates those results. The independent `Web base policy (trusted base)` executes base code. Local Node results do not claim the Python/browser/Gallery jobs passed; exact-head remote results are recorded at handoff.

## Remaining S05 acceptance

Preserve the original S05 result and historical failures:

- mapped Feature fill/direct edit tRNA legend conflict;
- divergent draft Save reflow busy/retry;
- real gzip transport maximum heartbeat 512–526 ms, above the unchanged 500 ms bound;
- real full Load/Save performance;
- real preview Python Worker count 1, where acceptance requires 0;
- real Generate: `Missing GenBank file`;
- strict CLI SVG root `baseProfile` difference;
- native clone bytes and receiver task duration remain null.

This guard removes only the inventory mismatch after its trusted merge. It does not complete S05, waive any failure, or start S06. Separate scoped deliveries must resolve the remaining acceptance boundaries while retaining original measurements, assertions, Product authority and mapped-contract separation.

## After a maintainer merges the guard PR

These are future integration instructions, not actions executed by this delivery. In the existing S05 checkout, first verify the branch, clean tree and actual runtime remote SHA. Stop if another session has changed them. Fetch only an advanced/missing dev ref. Confirm the delivered guard commit is an ancestor of the merged dev, then merge verified dev into the clean runtime branch; do not cherry-pick a candidate guard to self-authorize runtime.

```bash
cd /tmp/gbdraw-issue597-S05.ujOWyl
git status --short --branch
git ls-remote origin refs/heads/dev refs/heads/fix/issue-597-input-session-20260926
# If dev advanced, or its ref/object is missing:
git fetch --no-tags origin refs/heads/dev:refs/remotes/origin/dev
# Replace the placeholder with the commit in the delivery handoff:
git merge-base --is-ancestor <delivered-guard-commit> origin/dev
git merge --no-edit origin/dev
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
```

Resolve integration conflicts explicitly, preserve dev and S05 contributions, and record the integration SHA separately. Rerun applicable unchanged S05 contracts and acceptance measurements before claiming runtime admission. Any runtime commit/push or further PR requires the scope authorized for that later delivery. Guard merge alone does not authorize S06.

English commit/PR title: **Characterize the preauthorized Session import Worker owner**

Summary: Update the two exact Worker inventories for the already permitted import client, preserve absent-module behavior and mapped contracts, and record dev/S05 source checks and rejected constructor mutations separately from unresolved S05 acceptance.
