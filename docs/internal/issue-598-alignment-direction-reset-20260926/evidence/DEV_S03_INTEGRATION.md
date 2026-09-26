# Published S03 integration into dev

Base: `dev@aab5ad4323d51fc8ce94577b900b1c35e5ad7ffe`.
Source: published S03 `3902ef487c65827b8bd45306abfc4fd02794cd73`.
The source branch/checkout were read only. This dev-derived integration applies
the source commit without replacing PR #612/#613 fixes or the Issue #602 work.
The accompanying `S03.md` is the original source-session evidence, preserved
unchanged; this document records this integration's separate verification.

## Behavior and preservation

- Keep/right/left/Custom remain exclusive and include the exact reference.
  Direction scope, whole-record reversal, unchanged biological source strands,
  exclusions and before/after directions have readable accessible descriptions.
- Existing metadata tokens `undefined`, `unstranded` and `mixed` represent unknown
  strand through the existing parser. No direction is guessed. No-candidate
  records retain automatic missing-candidate provenance instead of creating a
  user Skip selection that Python cannot resolve.
- Reset has an independently scrolling body and reachable footer, keyboard radio
  grouping, target labels/current/restored directions, heading/error focus and
  invoker return through the existing focus owner. Missing historical evidence
  remains distinct from an empty current delta. Both Reset scopes consume the
  same latest receipt and retain the existing atomic artifact path.
- Retains fresh-load canonical source binding, pending intent rollback, plain-text
  labels, canonical label priority, Unicode receipt binding, successful retry
  error clearing and the existing compact CSS-marker/ResizeObserver wiring.
- Index conflicts retained one dialog description and target count. Browser
  conflicts retained both source S03 cases and the independent fresh-load Gallery
  case. The existing History invalidation revision remains in snapshots; Undo
  changes it, while both stacks and complete artifacts are checked.

Formal authority is unchanged: contract revision 21, SHA-256
`62e3a9c08ceb64acc349ace81a97a9c87dc45db4d18187a3b78e6e131b2dc595`.
PD-OI-027/029/031/034 revisions 5/3/5/5 and PD-OI-035 revision 3 / PD-OI-039
revision 2 contribute independently. This implements existing authority;
Match, new schema/readers, guessed orientation and new Product outcomes are not
introduced. S06 presentation remains on its own branch.

## Verification

```sh
node --test tests/web/similarity-alignment-actions.test.mjs tests/web/alignment-reset-receipt.test.mjs tests/web/run-analysis-simple-path.test.mjs tests/web/history.test.mjs
python tools/prepare_browser_wheel.py
NODE_PATH=/mnt/c/users/genom/github/gbdraw/node_modules GBDRAW_WEB_TEST_PORT=47208 node /mnt/c/users/genom/github/gbdraw/node_modules/@playwright/test/cli.js test tests/web/similarity-alignment-ui.playwright.spec.js --workers=1 --output=/tmp/issue598-s03-integration-browser
node tools/check-web-change-budget.mjs --base origin/dev
```

- Focused Web owners: **58 passed**, `/tmp/issue598-s03-integration-focused-node.log`.
- Full alignment browser first run: **10 passed / 1 failed**, 3.2m. The failure
  expected the superseded raw `skipped_by_user` UI text. The corrected case checks
  the readable `Skipped by user.` text and also retains the exact machine code.
- Corrected fresh-load reference case: **1 passed**, 28.6s,
  `/tmp/issue598-s03-integration-reference.log`. All 11 cases are covered across
  these runs; a later single-run whole-suite success is not claimed here.
- New source-backed cases cover unknown/no-candidate/skipped records, source
  immutability, reference-only/target-only receipts, both Reset scopes and Undo,
  keyboard grouping and focus, historical/empty evidence, and separate Apply
  after materially changed final facts plus validation failure/retry.
- Working-tree policy: **Gate PASS / Review CLEAR**. Syntax and whitespace pass.
- Python production bytes are identical to this dev baseline; prior PR #613
  Python session and Core PR evidence applies to unchanged code. No new whole
  Python, physical assistive-technology or release result is claimed.

## Ownership and rollback

The existing projection/parser, review/controller, Reset/artifact and focus
owners remain authoritative. Existing Apply and Reset bindings share the focus
bridge; there is no second controller, visibility store, renderer or Worker.
Owner/path/compatibility sets are non-increasing; no architecture exception or
privileged-rule edit. Production, tests and copied source evidence were reviewed
separately. Generated wheel, screenshots and logs are ignored/uncommitted.

Rollback: revert this S03 integration commit, retaining PR #612/#613 and authority.
