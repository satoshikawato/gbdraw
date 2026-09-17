# S08 follow-up: alias validation

2026-09-17 JST. **RW-04 implemented and regression-verified.** The initial
classification below records the plan before production edits.

## Classification and bounded plan

**IMPLEMENT_EXISTING_AUTHORITY**. Preserve the current accepted manifest domain,
normalized alias grouping, export ordinals, failure messages and continuations.
Base `5e0cb0fa` supplies OIPC-C04/C06/C07 and PD-OI-016/022 in
`OPTION_INTEGRITY_PRODUCT_CONTRACT.md`, and the saved-protein contracts in
`SESSION_COMPATIBILITY.md`. The inherited S08 writer repair is also explicitly
required by this request. Source/tests identify mechanisms, not new authority.
No new outcome, compatibility retirement or unresolved Product choice is
proposed; EVIDENCE_REQUIRED, PRODUCT_DECISION_REQUIRED and NOT_ALLOWED do not
apply to this equivalent local change. No BD decision is invented.

Owner and path remain `run-analysis.js` → `losat-cache.js` manifest validation.
Session reader/writer, raw/derived validation and index construction continue to
call this same validator. This changes no semantic owner, canonical path,
privilege, public signature, cache, persisted data or lifetime. Ordinary
non-increasing architecture evidence applies; no OE/PE/CB exception is needed.

**RW-04:** in `validateProteinIdentityManifest`, each valid feature alias executes
`metadata.displayAlias.normalize('NFC').trim()` in the metadata guard and again
in the next local `const alias`. The input is the same primitive string from
JSON/structured-cloned metadata, the second call follows synchronously, and
neither call nor intervening code writes it. Keep the metadata shape/string
guard, compute the existing `alias` once, then reject an empty alias before
inserting it into `featuresByAlias`. No validation is skipped or loan extended.

Preserve the early returns for malformed metadata and non-string aliases;
empty/whitespace-only aliases still return false before alias grouping, runtime
handle registration and export-ordinal checks. NFC (not NFKC), trimming,
case-sensitive grouping, sorted feature IDs and ordinal rules remain identical.
The manifest and original alias text remain unmodified. Side-effectful accessors
and monkey-patched built-ins are not Session/Worker data; this change introduces
no guarantee or reuse based on object identity for those objects.

Plan: characterize Unicode aliases, collisions, blank/wrong-type metadata,
ordinals and input immutability; count actual normalization operations and
confirm the counter fails on inherited source. Preserve RW-01–03 and all existing
identity/Session/run-analysis regressions. Apply `browser-offline-qa`, using the
existing manifest-merge runner with a bounded blank-alias extraction fault;
check source and freshly installed package at desktop/mobile sizes, saved-raw
hits/searches, result/provenance, failure isolation/retry, offline requests and
lazy Worker disposal. Reuse unchanged scientific and full lifecycle evidence
after source/reachability checks. No timing samples or warmups are planned.

S08 remains complete; S07 approved; S06 and S08 writer require future pre-merge
Review; S05 rejected; S07.5–S07.8 additional experiments remain closed.

## Source and inheritance

The inheritance source was clean at the requested commit
`6d3612534c1cd236d7fa3f130d4940c4536ae00c`, on
`perf/s08-manifest-merge-20260917`, with no upstream. Fetch left `origin/dev` at
`5e0cb0fa1d9920592da431b764fa2c9133b5bfdf`. Created
`.worktrees/s08-local-duplication-20260917` and branch
`perf/s08-local-duplication-20260917` from that base, without an upstream.

All ten commits between dev and the inheritance source were patch-unique.
`cherry-pick --ff` retained their original identities, in order: `3a3b29c7`,
`e3f93172`, `603e0787`, `f35f2e63`, `dc79915a`, `ec16110c`, `f52ed68e`,
`02ca8f95`, `8bf4b98a`, `6d361253`. No integrated Product authority was reapplied.
HEAD remains `6d361253`; RW-04 is the recorded uncommitted difference.
[Inheritance ledger](data/s08-alias-validation/inheritance.json) records this
boundary and the initial Git-area permission failure before successful
inheritance. All 375 prior source hashes and both prior test hashes matched
before new edits. No generated application wheel, environment, private Session or screenshot was
copied as source. Historical dirty/HEAD wording and audit scripts remain intact.

| Final owner | SHA-256 |
| --- | --- |
| `gbdraw/web/js/app/losat-cache.js` | `c36c46e615be45a27e112d28500bd5c2965a03acc1d70b03153e24c93507e4a3` |
| `gbdraw/web/js/app/run-analysis.js` (unchanged) | `62bbe9e886820059c7be9a26e0a6bb265ed1fb4f5bbc3dc4254d2865b72e25e9` |

Final source aggregate SHA-256 (375 sorted path→SHA-256 entries, encoded with
`json.dumps(hashes, sort_keys=True)`):
`ae6e34df07cd3acd5fc38cda9522a8dbeae1d4e9dcd89b24954580bdd6264fcc`.

## Investigation scope and retained boundaries

The bounded scan covered `run-analysis.js` saved-raw extraction, manifest merge,
key preparation, prepared pairs and derived admission; `losat-cache.js` current
getters, manifest/index/raw/derived validators; and their direct Session callers
in `services/config.js`. It stopped after proving RW-04. This is not a claim that
the whole application has no other duplicate work.

| Inspected processing | Finding / disposition |
| --- | --- |
| Alias NFC/trim within one validator iteration | Same string, content, synchronous lifetime and result; remove only the second evaluation (RW-04). |
| Input manifests → merged manifest → loop index | Inputs differ from their union; merge may create globally colliding runtime IDs. Later index construction follows key-helper awaits and owns a separate loan. Keep R input checks, one merged check, and index validation. |
| Repeated `getSeqEntry` / `getSeqHash` | Per-run maps already reuse extraction/hash results; call count is not parse count. Persistent extraction is keyed by File and semantic inputs. No change. |
| Committed raw lookup → completed-search retry → legacy promotion | Different possible entries and conditional continuations. Promotion awaits a Worker response and validates rewritten evidence. Keep entry-specific checks and error/miss behavior. |
| Raw classification → public reference validator | Shape checks recur, but the exported validator must also reject arbitrary direct inputs. Removing them needs a separate contract/mechanism change; not bundled with RW-04. |
| Derived getter/setter and Session reader/writer | Different payloads, publication/admission points and lifetimes; retain all checks. The S08 writer must still filter old source bindings before Save. |

The manifest validator is reached directly by merge, index construction,
`proteinRuntimeIdSets`, raw reference validation, derived validation, and
config's manifest adoption/export. `proteinRuntimeIdSets` has no production
caller but remains a public helper. Generate calls merge once with its reload
message; default direct helper users keep interleaved merge/validation and
early-conflict priority. Session import builds an index before admitting raw and
derived entries; export builds/releases its own index and validates the saved
manifest separately. No caller is bypassed. The change is inside each invocation,
so no content proof is borrowed across an await or a public boundary.

## Implementation, failure contracts and operation counts

Production changes only the alias guard in the existing validator. All previous
checks run in the same order: metadata shape, string type, NFC/trim nonempty,
alias grouping, runtime registration, extra metadata IDs, and sorted ordinal
verification. Invalid manifests still return false. Generate still maps an
invalid input to the exact reload error, before any merge conflict; the default
helper still reports an earlier conflict first. No new exception translation.

The new tests preserve NFC-composed/decomposed equivalence, whitespace trimming,
case and compatibility-character distinctions, zero-width non-whitespace,
wrong/missing metadata types, null/omitted ordinals for unique aliases, ordered
numeric ordinals for colliding aliases, rejection of string/reversed ordinals,
input immutability, raw getter rejection and merge failure priority.

[Red](data/s08-alias-validation/focused-red.log) ran these semantic assertions
against inherited production: they passed; the operation-count assertion failed
with 4 instead of 2. [Green](data/s08-alias-validation/focused-green.log) passes
on final production, including RW-01–03 and index success/error/Cancel/retry.
The probe wraps actual NFC/trim expressions in the imported production module;
it adds no runtime counter API or replacement validator.

| Scope | Before RW-04 | After RW-04 |
| --- | ---: | ---: |
| One successful manifest validation visiting P feature aliases | 2P NFC/trim evaluations | P |
| Observed valid fixture, P=2 | 4 | 2 |
| Manifest validator invocations | unchanged | unchanged |
| Empty/whitespace alias encountered | 1 evaluation then false | 1 evaluation then false |
| Invalid metadata shape/string type | 0 | 0 |

P counts visited features, not distinct alias strings. Invalid manifests can
short-circuit; the successful-input formula does not describe every failure.
The operation is NFC normalization followed by trim, not TSV parsing or a full
manifest traversal. RW-01–03 savings are not added again. There are zero timing
samples, warmups or new profiles; no seconds improvement, memory improvement or
whole-Generate speedup is claimed. S08 whole-operation waiting times do not
represent removable validation time.

RW-01's single getter TSV scan, RW-02's private deep-copy/index scope and every
release, RW-03's R+1 merge validation, empty-raw hits, strict 12-column/numeric/
finite checks and directional bindings all remain. No schema, scientific result,
ID/order, limit, self/reverse evidence, arithmetic order, provenance or Session
compatibility changes. No index enters Session, History, retry evidence or state.

## Verification and evidence reuse

[Command ledger](data/s08-alias-validation/commands.jsonl) records commands,
exit codes and diagnostic timestamps. [Environment](data/s08-alias-validation/environment.json)
records available tools; both Python Playwright and Node `@playwright/test` were
available. The existing Python S08 runner owns this acceptance path.

- [Seven Node scripts](data/s08-alias-validation/node-gates.log) passed: cache,
  derived identity, simple/derived Generate, Session export/reader, LOSAT settings
  and stable identity. [Schema gate](data/s08-alias-validation/schema-escalated.log)
  passed its three checks. New alias checks complement the inherited strict raw,
  merge precedence and index lifetime tests; no existing assertion was relaxed.
- [Architecture, RW-04 only](data/s08-alias-validation/architecture-escalated.log):
  PASS / Review CLEAR. [Cumulative from dev](data/s08-alias-validation/architecture-cumulative.log):
  PASS / Review REQUIRED for inherited Session/compatibility work. No authority,
  detector or registered owner/path changes were made by RW-04.
- The browser wheel was built from this checkout, and the application was
  installed into `.venv/alias-validation-package` using `pip --no-deps
  --no-build-isolation --target`. [Offline asset checks](data/s08-alias-validation/offline-assets.log)
  and focused Python lint passed. Runtime assets, URLs, CSP and packaging code
  are unchanged. No hosted analytics-enabled deployment was built or tested.
- [Source browser](data/s08-alias-validation/browser-source.json.gz) and
  [installed-package browser](data/s08-alias-validation/browser-package.json.gz)
  passed at 1280×720 and 390×844 in fresh contexts, with external requests blocked
  from navigation. All four have zero external requests and zero page errors.
  Preview-only load creates zero Workers. Saved-raw Generate uses the real
  extraction, key, analysis and rendering Worker: 13 hits, zero searches.
  Full SVG, geometry and provenance equal RW-03's corresponding saved-raw and
  retry outputs. Normal repeat uses committed typed artifacts; it is not called
  derived-cache reuse.
- The runner's existing manifest-merge mode accepts an optional blank-alias
  fault; its default schema fault and other modes are preserved. It gives
  identical biological bytes a fresh File owner and changes one alias to
  whitespace in a real extraction response. Successful responses are unchanged.
  Generate returns the exact reload error, dispatches no search/render, and
  preserves Result bytes, selection, History undo/redo, raw/derived caches and
  published manifest. Restoring the original File and forcing saved-raw retry
  succeeds. One Worker survives success/failure/retry; disposal leaves zero.

The first architecture and schema attempts hit sandbox `spawnSync git/python
EPERM`; the first browser attempt hit local socket EPERM before browser startup.
The same commands passed with escalation. Their failure logs remain beside the
successful logs. No test timeout or failure expectation was changed.

[Reuse record](data/s08-alias-validation/reuse.json) matches 374 of RW-03's 375
Python/JS source hashes, including `run-analysis.js`, Session writer/reader,
typed request, Worker and science owners. Only the alias validator differs.
Unchanged native/typed/reference checks, Circular, full Cancel/clear/History and
source-replacement Save→fresh Load owner transitions reuse the RW-03→RW-01/02→S08
evidence chain with its original skips and Review limits. New alias/Session Node
tests and real browser runs cover the changed validator; old results do not
prove that changed code by themselves. Explicit Run LOSAT followed by derived
reuse remains the separately verified RW-01/02 lifecycle continuation. No native
scientific matrix, full S08 timing matrix, long search or whole lifecycle rerun.

The [final audit](data/s08-alias-validation/audit.json) binds source and test
hashes, inherited/current diff categories, historical-file preservation,
worktree isolation, browser parity, wheel Python files and installed package
membership. Its script belongs only to the new evidence directory; no historical
audit script was run or rewritten.

## Diff review, rollback and handoff

The inheritance ledger and final audit distinguish all ten inherited commits
from RW-04. Review the current changes separately:

- Production: [patch](data/s08-alias-validation/new-production.patch), one alias
  guard in `losat-cache.js`; no new module, export, state, cache or helper.
- Tests: [patch](data/s08-alias-validation/new-tests-and-tools.patch), alias
  fixtures/counter in the existing cache test and one fault option in the existing
  browser runner. Existing success paths and historical modes are preserved.
- Documentation: [patch](data/s08-alias-validation/new-documentation.patch),
  this result plus additions to MASTER_PLAN and S08_REDUNDANT_WORK.
- Evidence: only `data/s08-alias-validation/`, including authored command/audit
  scripts and generated logs/reports. Fresh build/install outputs remain ignored.
  No reference output, Gallery image or social preview changed.

Rollback RW-04 with the reverse production patch and matching test/docs removal.
At this handoff, restoring `losat-cache.js` from `6d361253` also removes only this
production change and preserves RW-01–03. Do not restore from the older S08
commit. Rebuild generated packages after rollback; the S08 writer repair stays.

No next RW candidate is registered. Public raw shape checks were observed but
not removed: a future request would need to preserve direct-validator rejection
without new validation-bypass flags or duplicate owners. That investigation is
not an outstanding task. Shared 64-entry LRU, coverage union, numeric array
conversion, additional S07.5–S07.8 optimization and runtime redesign stay closed.

Proposed commit title: **Normalize protein manifest aliases once**

Summary: Reuse each normalized alias within manifest validation while preserving
strict input checks and failure precedence. Verify raw reuse and failure recovery
with offline source and packaged Workers, retaining RW-01–03 and Session behavior.

No new implementation commit, push, PR, integration merge, tag or deployment.
