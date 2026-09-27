# S04 — Existing Color rule pattern drafts

The implementation snapshot below records the pre-publication handoff. The
[publication addendum](#publication-addendum--2026-09-27) records the subsequently
authorized independent CI integration and supersedes the earlier hold.

S04 implements only **KEEP_REJECTED_PATTERN_DRAFT**. A rejected existing Color
rule pattern remains in its field with Not applied, a bounded cause, Retry and
Revert. Python evaluation and the accepted rule, Result, History, Session and
export owners retain their responsibilities. Local runtime verification is
complete. Publication is held at the inherited CI inventory boundary described
below; S05 has not started.

## Checkout, synchronization and authority

Dedicated session root: `/tmp/gbdraw-issue601-s04-YKjUuc`.
The independent clone, index, venv, Node dependencies, browser installation,
wheel, HTTP server and evidence are under this root. The shared workspace and
S02/S03 clones were not modified. The sole Issue #601 remote writer in this
session is this clone.

| Boundary | SHA / state |
| --- | --- |
| Start HEAD and same-named remote | `3854fdeee5d4f06e3e34b5a8b9ffdfc8fad554bc` |
| Initial fetched dev | `d313b70b9f97c2c1d70f9ae885edbead80b62021` |
| First later dev | `27939faebd4728c5aa52c9524d14ad0c3244190b` |
| Normal synchronization merge | `3e2e138ccf1bc2d75d718bed64f30d6e6ab73265` |
| Latest fetched dev | `a1450fdfd0cfb776645da74c142b89f6fc080276` |
| Authority-only normal synchronization merge / S04 parent | `6290d67aeed7aaa0b41d02d46083af8e98bdc68e` |
| Branch / upstream | `fix/issue-601-bug15-bug19` / `origin/fix/issue-601-bug15-bug19` |

Start S03 is an ancestor of the resulting candidate. Formal authority integration
`744be7a5943d4a247d027369629058898aa3f33e` remains an ancestor. The received #602
S07 runtime `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7` was already a work-branch
ancestor at start, but was not an ancestor of initial dev. PR #622 independently
integrated it into dev at `27939fae`; the first synchronization also receives its
later Session/Status, progress and bounded caption-History corrections. Merge
conflicts preserve S03 structured failures and dev progress/transact behavior,
both test sets, and the accepted History owner. The final dev increment is only
PD-OI-051 authority from PR #623; no in-flight comparison runtime is implemented.
No reset, rebase, amend, force push, PR or dev/main publication was performed.

Developer preflight: **IMPLEMENT_EXISTING_AUTHORITY**. PD-OI-046 and PD-OI-047
retain all nine fields, concern keys, revision 1 and their original Choice A.
Their JSON receipts are equal at start, latest dev and candidate; JSON-block
SHA-256 values are respectively
`7d657c657f404f44df458dd0d3b811946fa2b7e73d8e7caa802f2160d626d801` and
`259372a60be2477dc5ae02826e8c9a2b67cd3bc550360258e0cbf1df3befdc26`.
These hashes describe the JSON blocks, not the original signed text hashes.
The product-impact map/decision registry and architecture rules are unchanged
from the synchronized parent. No BD number, new authority, export retirement,
checker, allowlist, workflow or CI change is authored in S04. Received #598 and
#602 independent requirements remain; export handoff `f0c8128a` was read.

The one implementation commit is:

```sh
git show -s --format='%H %s' cf3ffe46f4422edbca0f4d3b0e11a51ff6bfb076
```

Its exact SHA, committed-head checker results and publication state are reported
in the post-commit handoff and external evidence. The synchronization merges are
separate from this implementation commit.

## Owner, canonical path and removed behavior

`feature-editor/rule-actions.js` remains the one rule mutation entry owner.
Its focused `pattern-drafts.js` helper owns only transient display text, pending
revision, normalized error and field identity. A Map keys drafts by the actual
canonical row object; WeakMap field IDs are never persisted. Bulk restoration
recognizes unchanged rows by their index and complete accepted field signature,
with a document epoch preventing restoration after a successful replacement.
It stores neither another canonical rule set nor an artifact snapshot.

Input changes update display state without evaluation. Ordinary field change
and Retry use the same row-bound retry action and existing
`prepareCandidate -> Python captions/normalization -> Python matching -> current
checks -> legend transact -> atomic rules/live SVG/Result commit` path. Existing
prepared reuse, priority, valid targets and one Worker are retained. The helper
adds no evaluator, validator, translator, dispatcher or compatibility reader.
Existing new-rule/preset/TSV recovery stays in its existing owner; it is not a
second owner for field drafts.

The old unconditional pattern `finally` restoration and index-based template
change call are removed. An accepted local mutation retains actual surviving row
objects, so focus and drafts follow row reorder rather than their previous DOM
index. A change from a removed/restored DOM row cannot edit its replacement.
Successful completion clears only its own current token; finishing History
cannot discard a newer keystroke. Transient field/status controls are excluded
from the generic History input adapter.

S03 `error-normalization.js` remains the only public wording/model owner, with
its finite enums, bounds and Python-character position semantics unchanged.
Syntax uses REGEX_SYNTAX and aria-invalid; initialization/resource preparation
failures retain WORKER_INIT/RESOURCE_INVALID and their actual stages. Pending
shows Checking and Not applied, never Applied or non-match. Retry executes once;
Revert clears the display draft, restores the current accepted value and focuses
the field. The live field description states that Save/Generate use the accepted
rule, Export uses current Result, and this draft is not saved.

A one-line correction in existing `diagram-generation.js` preserves the existing
cancel error through cold helper initialization. Without it, the received #602
cold-progress cancellation test was misclassified as WORKER_INIT by S03's adapter.
No transport or diagnostic owner is moved.

## Transitions

| Event | Draft / pending response | Accepted owners |
| --- | --- | --- |
| Keystroke | Update this row's text/revision, clear previous cause; no Python call | Unchanged |
| Normal change / Retry | One existing preparation attempt, pending token, current-row check | Commit only if valid/current |
| Syntax/runtime failure or failed correction/Retry | Keep text and classified bounded field cause | No Result/History mutation |
| Cancel/stale/superseded | Keep text; settle pending without publishing a new failure | No old commit |
| Successful edit | Clear only this attempt's current draft | Existing atomic commit, one History entry |
| Revert | Clear this draft; current accepted value and field focus | No evaluation or History entry |
| Drawer close/reopen or temporary mode change | Keep text/cause, invalidate pending revision explicitly | No invisible-mode old commit |
| Add/reorder/remove another local row | Preserve surviving actual row identities/drafts; invalidate old attempt | Existing transaction |
| Remove target or replace rule set | Reconcile and release deleted row's draft | Existing rule owner |
| Unrelated History restore | Capture/suspend; rebind unchanged accepted rows after owner restore | Existing History meaning |
| Target rule Undo/Redo replacement | Accepted signature changes; discard its draft; Redo does not revive it | Existing History meaning |
| Session replacement starts | Capture/suspend; field actions disabled while import pending | Existing import owner |
| Failed Session preflight or post-reset commit | Await source rollback, then restore matching old-document drafts | Exact prior Result/canonical/History |
| Successful Session/document Generate | Clear drafts and advance document epoch | Last accepted rule only |
| Failed/canceled/stale Generate | Restore surviving same-document drafts | Existing artifact rollback |
| Reset | Explicit clear before reset owner | Existing reset/History behavior |

History intent/checkpoint/artifact restores, Session import, Generate, reset,
mode controls and actual drawer close call explicit transitions from composition.
The mode watcher also covers external mode restoration; lifecycle correctness
does not depend on that watcher alone. Rule preparation's existing snapshot
checks catalog, biology, Result identity, mode, canonical rules, source files,
selected Result and editor evidence; the draft owner additionally checks actual
row, token, document epoch and import availability.

Save/fresh Load is tested in an independent context with zero diagram Worker
construction. Saved JSON contains only the accepted rule, not the sentinel draft.
Current Result survives Save/Load and failed imports. Export's entire SVG XML
matches current accepted Result; Generate uses the accepted rule and releases the
old document's draft. No Session schema or History snapshot ownership changes.

## Verification and evidence

All evidence is under `/tmp/gbdraw-issue601-s04-YKjUuc/evidence`. Environment:
Python 3.13.3, pytest 9.1.1, Node 26.8.2, Node Playwright 1.61.1, dedicated Python
Playwright and Chromium installation. Both Playwright CLI/import paths were
checked. The dedicated HTTP server at port 47804 supplies COOP same-origin and
COEP require-corp on every response, as required by nested Worker tests.
The gitignored browser wheel was prepared before tests, never regenerated during
checks, and no cache-bust token was changed. Bubblewrap's WSL host mount prevented
sandbox process creation; the same local checks ran with required escalation.

Commands use the dedicated venv on PATH, PYTHON pointing to its interpreter,
PYTHONPATH at this clone, PLAYWRIGHT_BROWSERS_PATH under the session root and
GBDRAW_WEB_TEST_PORT=47804. Browser invocations use Node's absolute installed CLI,
not Python's same-named `playwright` command. Long runs were monitored without
shortening test-owned timeouts.

| Check | Result / evidence |
| --- | --- |
| Required top-level fast Web Node selection | **871 passed**, 42.01 s; `node-required-final.log` |
| New native Python draft/race/duplicate/25k/cancel Node regressions | **20 cases**, included in 871; focused `node-final-cancel.log` also has 3 progress cases, all 23 pass |
| Architecture plus actual cumulative CI test selection | **203 passed / 1 inherited inventory failure**, 94.34 s; `node-architecture-ci-after-sync.log` (139 architecture + 65 CI) |
| Core/recipes/Gallery non-slow, non-browser Python, `-n 4` | **6632 passed / 17 skipped / 1 executable-permission failure**, 198 s; `python-core-recipes-gallery.log` |
| Same failed native fixture after dedicated-clone chmod preparation | **1 passed**; `python-fixture-recheck.log`; tracked executable bytes/mode restored unchanged |
| Python Playwright required browser marker after runtime synchronization | **38 passed**, 257.83 s; `python-browser-after-sync.log` |
| Existing real Python parity/label TSV/preset/Unicode/25k browser cases | **8 passed** in `browser-drafts-final.log`; subsequent S04 runs supersede that run's four UI failures |
| Final S04 real Worker/UI/History/Session/races/cold cancellation | **15 passed**, 1.9 min; `browser-s04-final.log` |
| Final four Circular/Linear desktop/390px UI cases | **4 passed**, 43.3 s; `browser-ui-final.log`, screenshots; included in final 15, not counted twice |
| Accepted compact Editor and Align Reset/direction/retry selection | **15 passed / 1 outdated raw-error expectation**; corrected same compact rollback case **1 passed**, 7 s; `browser-received-contracts.log`, `browser-compact-session-final.log` |
| Actual bounded PR smoke suite | **19 passed**, 3.8 min; `browser-pr-smoke.log` |
| Ruff | PASS; `lint-final.log` |
| Whitespace, receipt equality, generated/read-only artifact diffs | PASS; no tracked changes to references, social preview, dist, egg-info or wheel |

Primary reproduction commands:

```sh
node --test $(find tests/web -maxdepth 1 -name '*.test.mjs' ! -name architecture-contracts.test.mjs ! -name gallery-session-publication.test.mjs -print)
node --test tests/web/architecture-contracts.test.mjs tests/ci/*.test.mjs
../venv/bin/python -m pytest tests/ -m 'not slow and not browser' -n 4
../venv/bin/python -m pytest tests/ -m 'browser and not slow' --durations=10
node node_modules/@playwright/test/cli.js test tests/web/specific-color-pattern-drafts.playwright.spec.js --config=playwright.config.js --workers=1
node node_modules/@playwright/test/cli.js test tests/web/python-rule-parity.playwright.spec.js --config=playwright.config.js --workers=1
node node_modules/@playwright/test/cli.js test --config=playwright.pr-smoke.config.js --workers=1
../venv/bin/ruff check gbdraw/
```

The unchanged Python scientific code, native fixtures and recipe/Gallery inputs
permit reuse of the core run and isolated permission recheck; no all-green exit
is invented for the original failed command. Final dev's last increment changes
only authority documentation, so unchanged runtime evidence is retained.

The real 25,000-feature Worker test records 12,500 matches, **one matching helper
call**, synchronous prepared reuse and one Worker, with 564 event-loop ticks,
9064.8 ms preparation and 20.41 ms reuse. Native Node verifies the draft action's
one captions call plus one matching call, with no keystroke/Revert evaluation.
Each UI rejected edit plus failed Retry has exactly four existing helpers
(captions and matching per attempt), no extra Worker, and successful correction
adds exactly one History entry. Actual cold init failure, transport staging
failure, warm cancellation and cold cancellation/retry are separately tested.

Keyboard editing, Details, Copy/manual selection, Retry, Revert and retained
successful-edit focus pass at 1600 and 390 px in both modes. Field IDs associate
status/cause with the actual input. Final Linear 390px screenshot was visually
inspected: cause and recovery controls remain readable/reachable above Generate.
Sentinels appear only in the user-entered field; summary, Details, Copy and
console/page errors contain no private pattern/exception.

Failures corrected without relaxing accepted-owner assertions:

- Latest-dev progress exposed cold helper cancel misclassification; existing
  cancel classification is retained and real cold retry verified.
- Bulk Undo exposed an old DOM field change editing the replacement row at the
  same index. Ordinary change now uses actual-row recovery; all four UI cases
  and their one-entry Undo/Redo checks pass.
- XML empty-tag spelling differed across Session materialization. XML round
  trips retain every node/attribute, with no ignored geometry or tolerance.
  Existing mode-keyed stroke binding is settled before the held-edit baseline;
  the second mode cycle and late response preserve the exact entire Result.
- The first Export fixture used a valid non-matching `protein` rule on mtDNA;
  the final real journey uses matching NADH and compares the whole exported SVG.
- Two existing browser expectations required the old global field alert/raw
  Session exception. They now assert the field model and bounded unknown probe
  respectively. Known REGEX_SYNTAX/UNKNOWN_EXTENSION assertions, every Result,
  rollback, History and geometry check remain; the console assertion strengthens
  to no errors. A deliberately arbitrary injected exception is legitimately
  UNKNOWN; no known validation is reclassified.

## Gate, review and independent CI boundary

Trusted checker files were archived from latest dev
`a1450fdfd0cfb776645da74c142b89f6fc080276`, not modified candidate checker code.
Final working-tree cumulative comparison against that base is **Gate PASS /
Review REQUIRED**, no blocking violations. Cumulative production scope is 30
files, gross churn 1977, net additions 679. S04-only comparison against synchronized
parent `6290d67a` is also **PASS / REQUIRED**, seven production files, 232 additions /
32 deletions, gross 264 / net 200. Reports: `gate-final-working.log` and
`gate-s04-only-working.log`. The committed head is checked separately after this
commit using the exact latest-dev base and reported in the handoff. Review is
human attention for architecture/size/Session effects, not a Gate failure and
not fulfilled by this agent's audit.

Actual initial CI classification was profile=pr / impact=full, requiring
web-change-budget, core-pr, recipes-standard, gallery, lint, web-contracts-pr and
web-pr-smoke. The exact committed candidate is reclassified by latest-dev trusted
`ci-impact.mjs` in the post-commit handoff. Local results do not claim remote CI,
Python 3.11 or Node 20 success.

**Open CI prerequisite:** unchanged S03 HEAD already collects 19 PR smoke cases
against the inventory test's maximum 13 (1 passed / 1 failed at baseline).
Latest dev alone collects 13. S04 adds no PR smoke tags. The full real smoke
suite passes 19, but the actual unchanged cumulative CI inventory test still
fails, so CI admission/publication is not represented as complete.

The owner authorized preparing a separate CI branch/commit, without PR or dev
integration. Candidate `cf3486de68bd34484d75f5ad86c8e730b0ca40a9` on
`test/issue-601-pr-smoke-inventory-20260927` raises only that ceiling to 19 and
preserves the lower bound, full-suite membership and comparison separation.
Its **65 CI tests pass**; the corrected inventory against this cumulative
runtime tree **passes both cases**; trusted exact-base/head **Gate PASS /
Review CLEAR**. It remains local and is not merged into S04. Separate record:
`git show cf3486de:docs/internal/issue-601-pr-smoke-inventory-20260927.md`.
This compatibility result cannot certify S04's unchanged test as passing.
The remaining boundary is independent CI integration followed by final candidate
checks and the authorized same-named normal push. No implicit authorization is
inferred for a PR or dev merge.

## Separate audits, ratchet evidence and handoff

Production audit: one existing canonical mutation path; transient display owner
is focused under rule actions. Bulk callers request transitions instead of
creating their own draft stores. Python validity, Result admission, resource
staging, History snapshots and Session/export representations remain in their
existing owners. The accepted rule survives local mutation by identity; bulk
History remains canonical-data driven. Ordinary non-increasing owner/path
ratchet evidence applies: no owner/path excess, new persisted compatibility,
accepted violation, multiple canonical path or hard waiver is introduced.
The single helper is private decomposition of the rule field's new display
responsibility. Superseded unconditional pattern restoration and index-based
field change are removed; no fallback retains them. Rollback is the S04 commit
alone, retaining received authority and independent synchronization merges.

Test audit: native evaluator and real browser checks complement isolated race
fixtures; owner changes have distinct display/canonical/Result/History assertions.
No production failure is stubbed into successful validation. Timeouts and known
classification, matching, rollback and performance criteria are unchanged.
Docs audit: this result only; earlier records stay historical. Generated audit:
no committed wheel, browser output, temporary evidence, binary mode, reference
or public figure changes. The social preview remains owner-maintained.

S05 remains unstarted. Before its convergence, independently resolve the CI
inventory prerequisite, sync current dev normally, verify exact base/head and
required candidate checks, complete same-named publication, and read this
handoff plus the export owner's record. Export shared-file transfer remains S05.
Inherited limits remain: BUG-19's original audit pattern/build is unavailable;
there is no claim of reproducing that artifact. Export standalone module failure
promise caching still requires reload after module-fetch failure. Metadata-free
Session/Legend override limitations remain unresolved. These local results are
not supported-version, remote-CI, full dev staging or release evidence.

Commit title: **Keep rejected color-rule edits available for correction**

Summary: **Add focused transient draft recovery while preserving Python validation and canonical state.**

## Publication addendum — 2026-09-27

This addendum supersedes the publication hold in the implementation snapshot
above. The owner explicitly authorized an independent CI PR and its integration
into dev. The S04 implementation remains the single commit
`cf3ffe46f4422edbca0f4d3b0e11a51ff6bfb076`; its title and behavior are unchanged.

CI-only PR [#625](https://github.com/satoshikawato/gbdraw/pull/625), head
`3b7d86ae29021bab57df3c470a376287e14e9c69`, passed both required remote statuses,
`Web base policy (trusted base)` and `PR / gate`, before normal merge into dev.
Its dev integration commit is `98c21116f6439d5721e7ea62ae1a21a7cf2d4319`. The ceiling correction is now
accepted base content, rather than runtime self-authorization.

The work branch first received latest dev `494091aa` by normal synchronization
merge `1fd226232557ce2191e504e03bf8cbe9da436a27`. That increment contains the
independent Issue #597 Worker-owner characterization and evidence, not a Session
import runtime transfer. The updated architecture contract passes all **139**
cases against actual S04. A further normal synchronization receives the
independent CI integration and includes only this publication-record update as
local documentation. No S04 runtime or behavior-test assertions were changed by these
synchronizations; the CI inventory correction arrived through accepted dev. No reset, rebase, amend or force push was used.

Final local CI contracts: **65 passed**, including the actual cumulative
19-case inventory against its accepted ceiling. Latest-dev trusted checking of
the exact committed final candidate remains **Gate PASS / Review REQUIRED**;
the trusted CI plan remains **pr / full** with the seven required jobs listed
above. Full browser/Python/Node evidence is reused only for unchanged code,
inputs, environment and acceptance conditions. Runtime, generated artifacts,
scientific references and PD-OI-046/047 receipts remain unchanged from the
verified implementation. These statements do not claim S04 remote CI or a
completed human review.

Exact final synchronization HEAD, trusted base/head reports and the subsequent
normal push verification are recorded outside their own commit in
`/tmp/gbdraw-issue601-s04-YKjUuc/evidence/SESSION_04_PUBLICATION_HANDOFF.json`.
The authorized publication target is only
`origin/fix/issue-601-bug15-bug19`. S04 PR creation or integration is not
included. S05 remains unstarted; its handoff retains the implementation evidence
above and the export owner's `OWNER_HANDOFF_20260927.md`, with inherited limits
unchanged.
