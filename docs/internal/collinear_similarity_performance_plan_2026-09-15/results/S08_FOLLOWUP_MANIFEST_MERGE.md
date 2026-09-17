# S08 follow-up: manifest merge

2026-09-17 JST. **RW-03 complete**: local implementation, regression gates,
real source/package browser checks and final preservation audit passed.

## Pre-implementation classification

**IMPLEMENT_EXISTING_AUTHORITY**: remove duplicate validation while preserving
all user-visible outcomes. The base Product contract OIPC-C04/C06/C07 and
PD-OI-016/022 require correct stage identity, failure isolation and completed
raw retry. `docs/SESSION_COMPATIBILITY.md` requires strict manifest/raw/derived
identity and accepted Session formats. This request explicitly requires the
existing invalid-input message, reload continuation and failure precedence.
Current source/tests establish the callers and mechanism, not new authority.
There is no unresolved product choice or proposed format/compatibility change;
EVIDENCE_REQUIRED, PRODUCT_DECISION_REQUIRED and NOT_ALLOWED do not apply.

Owner/path: Generate remains in `run-analysis.js`; manifest validity and merge
remain in `losat-cache.js`. Move the Generate precheck into the existing merge
owner, with an optional invalid-input message. That option validates all inputs
before merging and preserves invalid-input precedence. The default call retains
interleaved validation/merge, including its early-conflict precedence. The
private boolean proving prevalidation is computed by the real validator inside
the synchronous call; callers cannot supply a skip-validation flag. The final
merged manifest is always validated. No await, input mutation, new owner,
public export, index lifetime, shared cache, persistence or Worker path is added.
The superseded caller check/import is removed in the same change. This is an
ordinary non-increasing architecture change; no full OE/PE/CB exception applies.

Caller inventory: one production invocation, in `run-analysis.js`; direct
helper calls occur only in `tests/web/losat-cache.test.mjs`. Preserve default
helper behavior for external users as well: non-array/empty input, merge order,
conflict messages, invalid-input errors, final-validation error and deep copies.
The optional message applies only to invalid inputs, never to conflicts or an
invalid merged manifest. Raw getter/index, Session reader/writer, extraction and
derived validation stay separate required boundaries.

Verification plan: direct helper fixtures and an instrumented actual Generate
merge block; real validator count R inputs plus one merged check; existing Node
identity/Session/run-analysis/RW-01/RW-02 regressions. Reuse unchanged science,
lifecycle and source-replacement evidence after hash/reachability checks. Run a
small real browser/Worker saved-raw Generate, repeat and invalid-manifest failure
with preserved Result/History and zero search/publication; source and newly
built package, local-only network, desktop/narrow viewport, lazy startup and
Worker disposal. No performance timing samples or new scientific benchmarks.

S08 remains complete; S07 approved; S06 and S08 writer require pre-merge Review;
S05 rejected; S07.5–S07.8 extra optimization remains closed.

## Source and inheritance

Fetched `origin/dev`: `5e0cb0fa1d9920592da431b764fa2c9133b5bfdf`.
Created `.worktrees/s08-manifest-merge-20260917`, branch
`perf/s08-manifest-merge-20260917`, without an upstream from that base.
`git cherry` found nine patch-unique dependencies. `cherry-pick --ff` retained
`3a3b29c7`, `e3f93172`, `603e0787`, `f35f2e63`, `dc79915a`, `ec16110c`,
`f52ed68e`, `02ca8f95`, and `8bf4b98a` in order, without reapplying integrated
Product authority. HEAD remains `8bf4b98ae02ef14f72caaebe3114f7a6fc3e0c75`.
The candidate is this HEAD **plus inherited and new uncommitted changes**.

The source worktree's five modified and 32 untracked files were inventoried and
SHA-256 checked before transfer; all 37 matched in the destination before any
new edits. The [inheritance ledger](data/s08-manifest-merge/inheritance.json)
and [inherited patch](data/s08-manifest-merge/inherited.patch) preserve that
boundary. No generated wheel, environment, private Session or screenshot was
copied as source. The browser wheel and installed package were built afresh.

| Final production owner | SHA-256 |
| --- | --- |
| `losat-cache.js` | `524ec4108c83461efc1cc7e3aa4d718687ea7998126892b125a2fd64041d0ceb` |
| `run-analysis.js` | `62bbe9e886820059c7be9a26e0a6bb265ed1fb4f5bbc3dc4254d2865b72e25e9` |

## Implementation and failure contracts

`mergeProteinIdentityManifests` accepts optional `{ invalidInputMessage }`.
For Generate, it runs the existing validator over the inputs before merging,
stopping at the first invalid input exactly as the previous caller precheck
did. Its local `inputsValidated` value is true only after that real check
succeeds. The synchronous merge therefore does not validate those inputs again.
The option never accepts an externally supplied validation result or identity
index. The default helper still validates and merges each input in order.

The caller supplies exactly:
`Protein comparison metadata could not be validated. Reload the page and try again.`
No broad catch translates conflicts or merged-manifest failures into that error.
The merged manifest's independent validation is unchanged. Schema, program,
outfmt, args, search context, directional hashes, protein-set references,
IDs/order, limits, self/reverse evidence, arithmetic, provenance and Session
formats are unchanged. RW-01 getter checks, RW-02 private index lifetime and
success/error/Cancel release, and the S08 source-replacement writer fix remain.

| Input condition | Generate | Default helper |
| --- | --- | --- |
| Invalid manifest first/middle/last | Original reload message | Original invalid-input error |
| Early conflict then later invalid input | Reload message wins | Earlier conflict wins |
| Valid inputs with conflicting identity | Original protein-set / analysis / instance conflict | Same |
| Valid inputs whose union has duplicate runtime handles | Original invalid-merged error | Same |
| Valid multiple/duplicate inputs | Same ordered values, deep copies | Same |

## Counts and measurement limits

The test instruments the **actual** manifest validator and executes the actual
Generate merge block. The inherited source failed the new count assertion;
[red log](data/s08-manifest-merge/focused-red.log). The candidate passes;
[green log](data/s08-manifest-merge/focused-green.log).

| Normal merge scope | Before RW-03 | After RW-03 |
| --- | ---: | ---: |
| Per-record input validation, R inputs | 2R | R |
| Independent merged validation | 1 | 1 |
| Observed fixture, R=3: input + merged | 6 + 1 | 3 + 1 |

This count excludes extraction, prepared-loop index construction, raw/derived
validation and Session boundaries. RW-01/RW-02 reductions are not added again.
There are **zero timing samples and zero warmups**. Browser/command durations
are verification evidence, not a speedup measurement. S08 Vibrio's 199.641 /
118.767 seconds remain historical whole-operation medians; preparation,
rendering and DOM/History time are not assigned to manifest validation.

## Verification and evidence reuse

[Commands](data/s08-manifest-merge/commands.jsonl) record exact executable gates.
The initial focused red/green commands were `node tests/web/losat-cache.test.mjs`
with output redirected to their named logs, before the ledger wrapper existed.

- The focused Node suite covers input immutability/nested deep copy, all invalid
  positions, three conflict kinds combined with invalid input, invalid merged
  output, default-helper behavior, exact Generate errors, and validator calls.
  Its inherited strict raw/getter and scoped-index lifetime assertions pass.
- [Node gates](data/s08-manifest-merge/node-gates.log): seven related scripts
  passed for raw/derived identity, run-analysis, Session export/reader, settings
  and stable identity. The eighth initially hit `spawnSync python EPERM`;
  [same schema gate with escalation](data/s08-manifest-merge/schema-escalated.log)
  passed all three checks. No product assertion was relaxed.
- [Offline asset check](data/s08-manifest-merge/offline-assets.log) passed.
  No asset URL, CSP, Worker protocol, analytics or runtime dependency changed.
- Architecture against inherited S08: [PASS / Review CLEAR](data/s08-manifest-merge/architecture.log).
  Against latest dev: [PASS / Review REQUIRED](data/s08-manifest-merge/architecture-cumulative.log)
  for the inherited compatibility-bearing work. This does not waive S06/S08 writer Review.
- [Evidence reuse](data/s08-manifest-merge/reuse.json) compares all prior 375
  Python/JS package sources: 373 unchanged, only the two target owners differ.
  Unchanged native science, typed/reference gates, source-replacement Save/fresh
  Load, real legacy migration, Cancel/retry and Circular/lifecycle results reuse
  [RW-01/RW-02's audit](data/s08-followup/audit.json) and its S08 evidence chain.
  The reader/writer, runtime/Worker, algorithms and these helper continuations
  remain unchanged. Their original skips, timing limits and Review obligations
  remain. No historical audit script or evidence was refreshed.

The new browser mode reuses S08's server, real transport observer, Gallery
fixture, preparation, Generate and disposal functions. It adds only a focused
fault-injection scenario, with success responses unmodified. The first attempt
hit the sandbox's local socket EPERM; the same check was escalated. A subsequent
harness attempt toggled display reversal, which correctly reused extraction and
never injected a bad response. That assertion failed, and its log/report are
preserved. The corrected fixture gives identical source bytes a fresh File
owner, then corrupts only the schema of a real Worker extraction response.
No production hook, synthetic success result or new search is substituted.

[Source browser report](data/s08-manifest-merge/browser-source.json.gz) and
[installed-package report](data/s08-manifest-merge/browser-package.json.gz)
pass at 1280×720 and 390×844, each in a fresh browser context. External requests
are blocked from navigation onward; attempted external requests and page errors
are both zero. Saved-preview load creates zero Workers. Saved-raw Generate has
13 raw hits and zero searches, and traverses real extraction/analysis/rendering.
Ordinary Generate repeat consumes the committed typed result; it is not labelled
derived-cache reuse. One Worker remains live across success, failure and retry,
and explicit disposal leaves zero live Workers in every context.

A real extraction response with an invalid manifest produces the exact reload
error, with zero search or render dispatch. Existing Result bytes, selection,
History undo/redo counts, raw/derived cache contents and published identity
manifest remain unchanged. Restoring the original File owner and forcing
saved-raw reanalysis succeeds with zero searches. Full SVG SHA-256, path geometry
and provenance match both the initial run and the archived RW-01/RW-02 saved-raw
run in all four source/package/viewport combinations.

The [final audit](data/s08-manifest-merge/audit.json) verifies all 208 Python
files in the new browser wheel against source and prior evidence, all 466
tracked files selected by the existing package-data manifest against the fresh
installation, and the embedded browser wheel. An early ad hoc check mistakenly
expected every Git-tracked file in the package and stopped at Web `CLAUDE.md`.
The final check derives required assets from `_build_support.py`, verifies
exact installed membership and byte equality, and retains the existing
exclusion of repository instructions, hosted Gallery and hosting-only files.
This is an installed GUI check with a local uploaded Gallery fixture; a hosted
analytics-enabled deployment was not built or tested. Unchanged Circular and
full Cancel/source-replacement lifecycle checks use the audited prior evidence.

All 32 inherited untracked evidence/result files remain byte-identical. The
shared dev, inheritance source, S08 and S07.8 worktrees retain their original
HEAD, status and tracked diff hashes. No unrelated worktree was edited.

## Diff audit, rollback and next work

Inherited commit work, the uncommitted RW-01/RW-02 patch, and RW-03 are reviewed
separately. The final audit writes separate `new-production.patch`,
`new-tests.patch`, and `new-documentation.patch` against the inherited dirty
source, so the earlier loop indentation changes do not obscure this change.

- Production: RW-03 changes only the existing merge helper and its one caller;
  it removes the caller's direct validator import/check. No RW-02 loop change.
- Tests: extend `losat-cache.test.mjs` and add a bounded integration mode to the
  existing browser runner. Existing lifecycle/timing modes are unchanged.
- Documentation: this result, MASTER_PLAN and the RW-03 inventory status.
- Generated evidence: only `data/s08-manifest-merge/`. Fresh package/wheel stay
  ignored; no tracked reference SVG, Gallery asset or social preview changed.

Rollback RW-03 alone by reversing `new-production.patch` (or the corresponding
merge/caller hunks), with its test/docs updates. **Do not restore the two entire
files from HEAD**, which would also lose RW-01/RW-02. The inherited patch and
hash ledger provide the exact rollback boundary. Rebuild packages afterward.
The S08 writer repair is outside these hunks and remains intact.

No additional nearby duplicate with proven same input/lifetime/result was found.
The later index build and Session/derived checks remain separate boundaries;
no new candidate is fabricated from similar-looking validators. Shared 64-entry
LRU, coverage union, numeric array conversion and runtime redesign remain out
of scope. S07.5–S07.8 additional experiments remain closed.

Proposed commit title: **Validate protein manifests once before merging**

Summary: Consolidate Generate input validation in the existing merge owner
while preserving failure precedence, reload guidance and merged validation.
Retain RW-01/RW-02 and Session behavior, and verify raw reuse and failure isolation
with real offline Workers in source and packaged applications.

No new implementation commit, push, PR, integration merge, tag or deployment.

## Commit handoff

The subsequent user request authorizes one local commit containing RW-01,
RW-02 and RW-03 together with their tests, documentation and evidence. The
commit adding this section is the complete handoff source; its parent is
`8bf4b98ae02ef14f72caaebe3114f7a6fc3e0c75`. The earlier HEAD/dirty wording and
historical audit-script HEAD assertions describe the pre-commit evidence
snapshot. Do not rewrite those assertions or rerun historical audits to make
past evidence appear newly generated. Production and test hashes were checked
against the final audit before staging. No push, PR, integration merge, tag or
deployment is authorized by this commit request.
