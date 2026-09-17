# S08 follow-up: raw validation

2026-09-17 JST. **RW-01 and RW-02 complete**, with focused regression, real
browser/Worker lifecycle, offline package verification and a final source audit.
S08 remains complete. S07 approval,
S06 and S08 writer pre-merge Review, and S05 rejection are unchanged.

## Source and inheritance

- Latest fetched `origin/dev`: `5e0cb0fa1d9920592da431b764fa2c9133b5bfdf`.
- Worktree: `.worktrees/s08-raw-validation-20260917`.
- Branch: `perf/s08-raw-validation-20260917`, no upstream, created from that base.
- Inherited HEAD: `8bf4b98ae02ef14f72caaebe3114f7a6fc3e0c75`.
  `git cherry` found nine patch-unique dependencies, and `git cherry-pick --ff`
  retained their original objects in order: `3a3b29c7`, `e3f93172`, `603e0787`,
  `f35f2e63`, `dc79915a`, `ec16110c`, `f52ed68e`, `02ca8f95`, `8bf4b98a`.
  `origin/dev` was the merge base and had no additional commits to reconcile.
- The S08 report's older HEAD/dirty wording is historical: its complete source,
  tests and evidence are now committed in `8bf4b98a`.
- [Inheritance ledger](data/s08-followup/inheritance.json) records the checks.
  The two target owner hashes matched S08 before edits. A full file comparison
  overlapped the first new test edit and reported only that test; the inherited
  Git tree itself is exactly S08. No wheel or environment was copied as source.
- [Final audit](data/s08-followup/audit.json) records final source hashes, the
  aggregate source digest, installed files, wheel, and preservation checks.
  The final candidate is inherited HEAD plus the recorded local diff.
  Aggregate source digest (sorted path→SHA-256 JSON):
  `15550c25f924488ec4feb113993872445a4b222e7b2010f7e5b4eb13a3d2e706`.

| Final production file | SHA-256 |
| --- | --- |
| `gbdraw/web/js/app/losat-cache.js` | `4794a5e5473b8270248ee4a9460694aeada1418b698ad6e17bbad9701d04da24` |
| `gbdraw/web/js/app/run-analysis.js` | `af9f9a346aba3c0e8d814d4b27486e3f744f4fb291ce8c40f55f0384147f271a` |

## Preflight and ownership

Classification: **IMPLEMENT_EXISTING_AUTHORITY**. The base's
`OPTION_INTEGRITY_PRODUCT_CONTRACT.md` OIPC-C04/C06/C07 and PD-OI-022 require
correct identity, explicit replacement/clearing, completed-search retry and
failure isolation. `SESSION_COMPATIBILITY.md` specifies directional protein
identity, strict current schemas and supported legacy promotion. The S08 writer
repair preserves the source-replacement save/load outcome and is inherited intact.
There is no alternative product outcome to choose, new compatibility retirement,
new authority, or unresolved behavior. EVIDENCE_REQUIRED,
PRODUCT_DECISION_REQUIRED and NOT_ALLOWED therefore do not apply to this change.

The existing `run-analysis.js` owner lends the existing `losat-cache.js` index
only during prepared-pair lookup. Raw validation remains owned by
`validateProteinRawEntryReferences`; strict row parsing remains owned by
`rawProteinTextMatchesBindings`. No owner, canonical path, Worker message,
shared cache, persistent state or compatibility path is added. Ordinary
non-increasing owner/path evidence applies; none of the architecture ratchet's
full OE/PE/CB exception conditions is introduced by this follow-up.

The architecture gate against `origin/dev` passes with Review REQUIRED for the
inherited compatibility-bearing change. Registered owner/path, module, export,
privilege and compatibility counts do not increase. This does not waive S06 or
S08 writer Review. See [cumulative gate](data/s08-followup/architecture-escalated.log).
The follow-up-only differential against `8bf4b98a` is **PASS / Review CLEAR**;
[follow-up gate](data/s08-followup/architecture-followup.log).

## RW-01: completed

Removed the second `proteinRuntimeIdSets` call and second TSV scan from
`getCurrentRawLosatCacheEntry`. The preceding reference validator already
validates exactly the same manifest, bindings, runtime IDs and text without an
intervening await or mutation. The getter still checks search context, schema,
program, outfmt, args and directional metadata before reference validation.

The public `proteinRuntimeIdSets` export is retained. A repository search found
no remaining production caller after removing this getter's duplicate, but its
public API is not removed as part of a local runtime optimization.
[Caller inventory](data/s08-followup/getter-callers.log).

The counter regression failed on inherited source with manifest=2 and TSV=2,
then passed with 1 and 1. The getter tests exercise matching forward/reverse,
wrong direction, unknown query/subject IDs, invalid manifest/set/binding,
program/outfmt/args/search-context mismatch, wrong column counts, the complete
shared strict-numeric fixture and empty-text success.

## RW-02: completed

The getter accepts an optional fifth `{ identityIndex }` argument and forwards
it to its existing validator. Four-argument direct callers retain full
validation. The orchestration wrappers also forward the index to the committed
raw cache and completed-search retry cache. Verified legacy promotion uses the
same index. Each actual entry retains its own binding and TSV validation.

The lifetime proof is content/data-flow based, not object identity alone:

1. Before this loop, all selected protein records are extracted and
   `mergeProteinIdentityManifests` creates a new manifest. `mergeManifestMap`
   JSON-clones each incoming value, including nested runtime IDs and metadata;
   there are no mutable aliases to source entries or reactive state.
2. The only assignments to `workingProteinIdentityManifest` are the initial
   owner read and this fresh merge. Protein prepared jobs always pass the merge.
   Between merge and loop completion there are no writes to the merged contents.
3. Worker operations, including legacy reference resolution, cache-key creation
   and legacy raw promotion, receive `cloneJsonData` copies. `getSeqEntry` and
   `getSeqHash` read the separately prepared sequence maps. Promotion changes raw
   entries/candidate transactions, not the merged manifest.
4. Build the index after cache keys are prepared. A `try/finally` contains the
   pair loop, including every await and legacy promotion. Success, preparation
   exception and Cancel release the WeakMap registration. An invalid manifest
   produces no index and an explicit error.
5. Release happens before source search, conversion, rendering and candidate
   publication. The index is never stored in state, retry evidence, Session or
   History. A later run extracts/merges anew. External replacement during an
   await cannot mutate this run's private manifest; existing cancellation/stale
   guards still own publication. No new index can authorize another manifest.

The existing index API is a scoped loan, not a mutation-tracking cache. Callers
must keep its manifest unmodified until release. This change establishes that
condition for the async orchestration owner; it does not add a global guarantee
for arbitrary external callers mutating borrowed objects. Tests execute the
production loop with the real validators, mutate the external source while
awaiting, exercise promotion, and interrupt after one successful lookup. They
verify one build, release on success/error/Cancel, retry acquisition, rejection
with a different manifest, and rejection after mutating the same object once its
old index is released. Real Worker legacy migration is checked separately.

## Retained validation and operation counts

All entry-specific TSV checks remain strict: exactly 12 columns, decimal numeric
syntax, integer-only columns, finite values, both runtime IDs. Empty raw remains
a successful hit. Protein sets, directional binding hashes, record references,
search context, raw and derived identities, limits, self/reverse evidence,
scientific numeric order and provenance are unchanged. Session reader/writer and
derived validation are separate boundaries and remain intact.

Counts below cover successful prepared-loop raw lookups only. They exclude
extraction/merge, Worker helpers, Session import/export and derived validation.
Runtime Sets means the query/subject ID membership Sets, not every Set allocated
inside manifest validation.

| Scope | S08 | RW-01 only | RW-01 + RW-02 |
| --- | ---: | ---: | ---: |
| H valid hits, R manifest records: full manifest checks | 2H | H | 1 |
| Runtime ID Sets | 4H | 2H | R |
| TSV scans | 2H | H | H |
| S08 Vibrio's H=47, R=11: manifest / Sets / TSV | 94 / 188 / 94 | 47 / 94 / 47 | 1 / 11 / 47 |
| Hep ON H=13, R=5: manifest / Sets / TSV | 26 / 52 / 26 | 13 / 26 / 13 | 1 / 5 / 13 |

[Operation-count record](data/s08-followup/operation-counts.json).
The getter and production-loop probes observe these counts on a small fixture
(3 hits, 2 records: 1 manifest check, 2 runtime Sets, 3 TSV scans). Gallery counts
are source-derived applications of that verified rule to the actual directional
table inventory, not a new browser profiler. Invalid entries may fail before a
TSV scan or cause separate committed/retry-cache lookups, so the table is not a
formula for every error path. Legacy promotion still validates its returned raw
entry. On an all-miss run RW-02 adds one manifest validation and R transient
Sets; it retains those Sets until loop exit instead of allocating two at a time.
No peak-memory reduction or all-workload speedup is claimed.

No new performance timing samples, warmups or profiles were necessary. All
command durations in the ledger are verification/build durations. New real raw
searches occur only where lifecycle checks require changed scope/settings,
Clear Cache or biological source replacement; this is not a raw-search timing
matrix. Tests and builds were allowed to overlap correctness checks, so those
elapsed values must not become performance samples.

S08's Vibrio medians (saved raw 199.641 s, derived reuse 118.767 s) remain prior
observations, not removable seconds. There is no matched new whole-Generate
baseline or seconds-improvement claim. Preparation, rendering and DOM/History
are not attributed to these validators. Ordinary Generate consumes committed
typed artifacts; explicit Run LOSAT followed by Generate exercises derived reuse.
Both continuations are checked separately.

## Verification and reused evidence

Commands, failures and exact arguments are in
[commands.jsonl](data/s08-followup/commands.jsonl). New outputs live only in
`data/s08-followup/`; S08 outputs and its historical audit scripts are unchanged.

| Gate | Evidence |
| --- | --- |
| Getter regression and index lifetime | [red](data/s08-followup/rw01-red.log), [RW-01 green](data/s08-followup/rw01-green.log), [partial-loop error/Cancel](data/s08-followup/index-partial-lifetime.log) |
| Raw/derived identity, run-analysis, Session export/reader, settings, stable identity | Seven Node files passed: [log](data/s08-followup/node-final.log) |
| Current Session/request/raw/derived schemas | Three assertions/tests passed: [log](data/s08-followup/schema-gate.log) |
| Real old-schema promotion, save/load/export, failure/Cancel rollback | 22,567 assertions passed: [log](data/s08-followup/legacy-browser.log) |
| Source SPA lifecycle | [report](data/s08-followup/lifecycle.json.gz), [log](data/s08-followup/lifecycle-escalated.log) |
| Installed package offline desktop/mobile | [report](data/s08-followup/offline-package.json.gz), [log](data/s08-followup/offline-package.log) |
| Offline asset inventory | [log](data/s08-followup/offline-assets.log) |
| Final source, S08 parity, wheel/package and preserved worktrees | [audit](data/s08-followup/audit.json) |

The unchanged S08 scientific evidence is reused after checking all 208 Python
source files and the current browser wheel. All 466 installed tracked package
files and the embedded wheel match the new worktree; the shared dev, S08 and
S07.8 HEAD/status/tracked-diff hashes remain unchanged. The S08 native/typed/reference/full
suite and frozen whole-result matrices are not rerun for this JavaScript-only
change. Their original skips, measurement limits and Review requirements remain.
The new browser checks traverse the changed lookup path and production Python
helpers/rendering. They supplement, rather than replace, those scientific oracles.

Offline checks apply `browser-offline-qa`: fresh contexts, local-only requests
blocked from navigation onward, cold saved-preview load, successful generation,
repeat and disposal. The installed package checks cover 1280×720 and 390×844,
Linear protein and Circular, preview readiness and retained Worker reuse.
No hosted deployment bundle or analytics-enabled deployment was built or tested.
No asset, URL, CSP, network dependency or Worker protocol changed. Initial local
server and Node child-process attempts hit sandbox EPERM; the same checks passed
with the required sandbox escalation, without weakening the assertions.

The source lifecycle command passed all 13 steps and a fresh-load context.
Saved-raw Generate has 13 hits and zero searches. Ordinary repeat uses committed
typed resources; explicit Run LOSAT then Generate has one derived hit and zero
searches. Eleven steps matching the S08 archived conditions have identical full
SVG SHA-256, path geometry and provenance. The existing recipe also checks
color/block/member/filter/reversal, additional all-record evidence and reorder.

The real render-response hook delays an actual completed Worker response; Cancel
preserves Result and History and terminates the Worker. Member change and retry
use the 13 completed raw entries with zero searches and a newly created Worker.
Clear Cache empties raw/derived state and the explicit Run LOSAT action performs
13 new searches and reproduces the same SVG, geometry and provenance as the
current retry at unchanged settings. Active raw keys match the provenance; the
whole cache correctly shrinks from 38 entries to 13 by dropping inactive evidence. The first audit attempted to compare
this step with an old S08 prefix labelled clear-cache; that prefix actually used
committed typed resources and had no searches or derived provenance. The audit
corrects that non-equivalent comparison, records the initial failure, and uses
the current before/after Clear Cache invariant without rerunning a workload.
A second audit assertion incorrectly expected the whole raw-cache key set to
survive clearing; it was corrected to distinguish the 13 active keys from the
38 retained entries. Both audit failures are recorded; no runtime change or
extra browser execution was needed.
Biological source replacement performs only changed-source jobs,
History undo/redo restores the generated result, and Save followed by fresh Load
succeeds. The fresh context starts with zero Workers; forced saved-raw reanalysis
has zero searches and reproduces changed-source geometry and provenance. This
specifically guards the S08 writer defect. Every source/installed context has
zero external requests, zero page errors and zero live Workers after disposal.
The installed-package Linear result/provenance also matches the source flow.

## Diff audit, rollback and next candidate

- **Production:** only `losat-cache.js` and `run-analysis.js` differ from S08.
  Most orchestration churn is indentation under `try/finally`; `git diff -w`
  exposes the small substantive change. No Python, writer, Worker or schema edits.
- **Tests:** only `tests/web/losat-cache.test.mjs` changes. It extends direct getter
  coverage and instruments actual source for deterministic operation counts and
  loan-lifetime checks. Async I/O in the loop probe is controlled; real migration
  and scientific execution use the unchanged existing browser runners.
- **Documentation:** this result, MASTER_PLAN and S08_REDUNDANT_WORK. Historical
  S08 reports are preserved; their commit-time context is clarified by additions.
- **Generated evidence:** the new `data/s08-followup/` logs/reports and its
  read-only audit script. Generated wheel, new environment, private Sessions and
  screenshots stay ignored. No references or social preview were changed.

Rollback RW-02 by restoring `run-analysis.js` to `8bf4b98a` and removing only the
optional index argument from the getter, retaining RW-01's duplicate removal.
Rollback both by restoring these two production files to `8bf4b98a`, with matching
test/documentation changes. Neither operation rolls back the inherited S08
writer fix. Rebuild generated packages after any rollback.

RW-03 is documented in the redundant-work inventory: the per-record manifest
precheck immediately before merge repeats the merge owner's own input validation.
It is not changed here. Error wording and ownership must be considered before
removing it. Shared 64-entry LRU, coverage union, numeric array conversion,
render/runtime redesign and closed S07.5–S07.8 optimization remain out of scope.

Proposed commit title: **Avoid repeated validation of protein raw cache entries**

Summary: Remove duplicate raw-entry validation and reuse the existing protein
identity index within each prepared-pair loop. Preserve strict entry checks and
Session compatibility, release the index on every exit, and verify real Worker
retry, source replacement and offline workflows.

No new implementation commit, push, PR, integration merge, tag or deployment was
performed. The nine inherited commits retain their original identities.
