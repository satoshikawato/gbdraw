# Issue #601 — separate PR smoke inventory correction

This CI-only candidate raises the existing bounded PR smoke inventory ceiling
from 13 to 19. It preserves the lower bound of 8, every smoke case's membership
in full functional acceptance, all moved-regression assertions, the 16 mandatory
comparison contracts, and the prohibition on running those contracts twice.
No runtime, workflow, timeout, test tag, authority, checker or allowlist changes
are included. The separate S04 runtime candidate adds no PR smoke tags.

The owner authorized preparation on a separate branch and commit on 2026-09-27.
That authorization does not cover push, PR creation, dev integration or mixing
this correction into S04. This local candidate is ready for independent review.

## Cause and scope

Fresh branch `test/issue-601-pr-smoke-inventory-20260927` is derived without an
upstream from latest fetched `origin/dev`,
`27939faebd4728c5aa52c9524d14ad0c3244190b`. Worktree:
`/tmp/gbdraw-issue601-s04-YKjUuc/ci-inventory-fix`.
Latest dev collects 13 cases and passes the existing inventory check. Issue #601
S03 HEAD `3854fdeee5d4f06e3e34b5a8b9ffdfc8fad554bc` collects 19 and fails
that check before S04 edits. The six existing Issue #601 additions are four
native Generate error cases (Circular/Linear at desktop/390 px), standalone
JavaScript search/word-target parity, and native Python Align failure/retry.
The same 19 remain after latest-dev synchronization. Removing these accepted
functional scenarios to satisfy an obsolete cap would lose existing coverage.

Only `tests/ci/playwright-inventory.test.mjs` owns the correction: its title and
upper bound agree with the existing cumulative inventory. No new execution path
or production architecture is introduced.

## Verification

Evidence directory: `/tmp/gbdraw-issue601-s04-YKjUuc/evidence`.
Node 26.8.2; Node Playwright 1.61.1, using the dedicated dependency installation.

- Unchanged S03 baseline: 1 passed / 1 failed; `baseline-inventory.log`.
- Unchanged latest dev: 2 passed; `ci-inventory-latest-base.log`.
- Separate candidate, `node --test tests/ci/*.test.mjs`: **65 passed**;
  `ci-fix-all.log`.
- Exact corrected inventory test executed against the synchronized S04 working
  tree, `node --test ../ci-inventory-fix/tests/ci/playwright-inventory.test.mjs`:
  **2 passed**; `ci-fix-cumulative-projection.log`. This is companion-candidate
  compatibility evidence, not a claim that S04's unchanged CI test passes.
- Latest-dev trusted checker, archived from the exact base above:
  **Gate PASS / Review CLEAR**; `ci-fix-gate.log`.
- Production, tests, documentation and generated scope reviewed separately;
  whitespace check passes; there are no production or generated changes.

The final commit can be identified by this file and title
`Align PR smoke inventory with accepted Issue 601 coverage`; its SHA is reported
in the handoff instead of being embedded in its own content. Required external
CI and integration are not claimed. S04 publication remains blocked by its own
unchanged 13-case inventory ceiling until the independent CI boundary is resolved.

Proposed commit title: Align PR smoke inventory with accepted Issue 601 coverage

Summary: Account for six existing smoke scenarios while preserving full-suite
membership and comparison-contract separation.
