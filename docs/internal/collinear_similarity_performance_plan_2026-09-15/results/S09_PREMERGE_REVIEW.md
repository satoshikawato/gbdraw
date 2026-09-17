# S09 — S06 compatibility / S08 Session writer pre-merge review

2026-09-17 JST. S09 is a new review stage, not another optimization session.
One writer integrity defect was reproduced and corrected in its existing owner.
Review, focused verification and final preservation audit are complete. **Human S06 compatibility-exception Review and
S08 writer Review remain required; this report grants neither approval.**
S08 remains complete, S07 approved, S05 rejected, and S07.5–S07.8 additional
optimization/performance-adoption work closed.

## Candidate and inheritance

- Fetched base: `5e0cb0fa1d9920592da431b764fa2c9133b5bfdf`; no new origin/dev
  commits since the previous base. The shared dirty dev checkout was not edited.
- Dedicated worktree: `.worktrees/s09-premerge-review-20260917`.
- Branch: `review/s09-premerge-20260917`, no upstream, created from origin/dev.
- Inherited HEAD: `6d3612534c1cd236d7fa3f130d4940c4536ae00c`.
  Ten patch-unique commits were inherited once using `cherry-pick --ff`, without
  reapplying merged Product authority or creating replacement commit identities.

| Commit | Dependency retained |
| --- | --- |
| `3a3b29c7` | S01 baselines and reproducible evidence |
| `e3f93172` | S02 path contracts and oracle |
| `603e0787` | S02 PATH-B handoff |
| `f35f2e63` | S03 HSP aggregation |
| `dc79915a` | S04 support and metadata indexing |
| `ec16110c` | S06 lossless paths and compatibility |
| `f52ed68e` | Approved S07 cluster merge |
| `02ca8f95` | S07.6/S07.7 local reductions |
| `8bf4b98a` | S07.8 candidate, S08 writer repair and integrated verification |
| `6d361253` | RW-01–03 and their tests/evidence |

[Inheritance ledger](data/s09-premerge-review/inheritance.json) contains full
SHAs, actual starting status, protected worktree snapshots, and every inherited
uncommitted file's size/hash. RW-04 comprises **37 uncommitted files**: one
production, two tests/tools, three documents, and 31 evidence files. All 375
source hashes, both RW-04 test hashes, and the 30 entries in its self-excluding
evidence inventory matched before edits. The inventory file itself was also
hashed and transferred. [RW-04 patch](data/s09-premerge-review/inherited-rw04.patch)
is separate from S09. No inheritance conflict occurred.

Inherited aggregate: `ae6e34df07cd3acd5fc38cda9522a8dbeae1d4e9dcd89b24954580bdd6264fcc`.
Final source aggregate: `2786f88d13149749167f4f3fc5fa12c8d6bdba5e611ab06aa3c99c5d232daa6d`.
Both use the same 375 sorted path→SHA-256 entries encoded by
`json.dumps(hashes, sort_keys=True)`. Final `services/config.js`:
`658577d078818af2f13feb2f1e726a4813fcc9adbc7e3379a86878d31a03b5cc`.
The candidate is **HEAD + uncommitted RW-04 + uncommitted S09**, not an approved
exact-head commit. Generated wheel/environment/private saved Sessions were not
inherited as source; packaging was rebuilt locally.

## Review authority and scope

Read MASTER_PLAN, S06, S08, S08_REDUNDANT_WORK and all three RW follow-up reports,
the root/Web instructions, architecture/Product ratchets, publication inventory
and relevant raw evidence. Authority was checked against **base**, not candidate
documentation:

- OIPC revision 6 / **PD-OI-023**, concern
  `protein-comparison.path-representation`, scenario 1, selects PATH-B. It keeps
  scientific output, every path's content/order/ID/shared information, supported
  old-file input and explicit exhaustive tuple access. It permits retiring the
  always-exhaustive normal API/storage representation only.
- Base `docs/SESSION_COMPATIBILITY.md` separates Session/request/bindings and
  protein raw/derived/manifest namespaces and requires exact directional identity.
  Its draft/committed distinction and OIPC-C01/C04/C06/C07 apply to Save and retry.
- PD-OI-016/018/021/022 retain failure isolation, complete scope/search identity,
  inference settings, and completed-raw retry. No limits or scientific rule changed.
- Base architecture registry/Product map retain the canonical request and current
  Result entry/owner subjects. Both independent requirements of each mapped
  option remain satisfied; matching the option ID alone was not used as proof.
  `tools/web-product-decisions.json` has no active decisions; no BD ID is invented.
- The candidate Session compatibility additions describe S06 and the S08 repair;
  they do not authorize their own runtime. OIPC, ratchets, registrations, checker
  and instruction files remain byte-identical to base.

## Finding S09-F01 — medium: invalid metadata silently saved as empty

**Trigger:** a null/malformed live protein manifest with an empty raw cache, or
one containing only nucleotide raw. The current-protein-entry case already failed.
The new regression calls the real exported `exportSession`, with valid remaining
Session state. No normal UI sequence causing the corrupted manifest was established;
this is fault-injected writer-boundary coverage, not a claim of routine data loss.

**Location:** inherited `config.js:2420–2467` (`serializeLosatCache`) and
`:4007–4014` (the final manifest fallback). The raw filter checks missing/invalid
manifest state only for `protein-current` entries. With no such entry, the final
writer replaced invalid metadata with `emptyProteinIdentityManifest()` and returned
`saved`, producing a download. Expected: an explicit valid-manifest error, no
download, and preserved state. A valid empty schema-2 manifest remains saveable.

**Basis and classification:** [pre-edit packet](data/s09-premerge-review/preflight.md),
OIPC-C01/C04/C06 and the requested strict writer contract select one outcome:
`IMPLEMENT_EXISTING_AUTHORITY`. This is an ordinary regression against static
authority; it restores rejection of invalid state and preserves every valid
continuation. It is not NO_USER_VISIBLE_DIFFERENCE for the corrupted-state case.
No Product choice, new authority or compatibility retirement is needed.

**Fix:** `config.js:3945` now throws the existing manifest error at the final
validation point, and `:4010` retains the existing valid-manifest adoption/copy.
Only the silent fallback is removed. The earlier raw filter and its `try/finally`,
reader, derived validation and Generate checks are unchanged. No owner, canonical
path, schema, cache, loan lifetime or protocol is added; ordinary non-increasing
architecture evidence applies. [Production patch](data/s09-premerge-review/new-production.patch).

The [red regression](data/s09-premerge-review/writer-red.log) failed with
“Missing expected rejection”; [green](data/s09-premerge-review/writer-green.log)
passed. The final focused gate adds six invalid empty/nucleotide-only cases,
asserts no compression/download and preserved state, and verifies two successful
valid-empty saves admitted by the current reader. Real browser fault cases use
null and blank-alias manifests in both cache states, then restore valid state
and verify successful saving/replay. No additional runtime defect was established.

## S06 compatibility Review materials

| Boundary reviewed | Current result and evidence |
| --- | --- |
| Default producers / explicit legacy access | `protein_colinearity.py:395–446`, public producers `:6503`/`:6621`: graph default; `materialize_ortholog_paths` and explicit `path_representation="exhaustive"` return real legacy tuple fields. No lazy tuple substitution. |
| DAG / explicit corpus | `ortholog_paths.py`: DAG validation, deterministic edge choice, tied-key ordering, rank/unrank and first-edge IDs; explicit collections preserve supplied paths, custom/duplicate IDs, cycles, order and shared fields without edge closure. |
| Exact counts / normal consumers | Count uses Python integers and canonical decimal strings on wire. Metadata uses count summaries, not implicit iteration. Large-count and no-expansion oracles remain in `test_lossless_ortholog_paths.py`. |
| Typed resources 1/2/3 | `session_request_codec.py:3516–3710`: writer 3; readers 1/2 decode old dataclasses once to explicit collections (schema-1 defaults retained); reader 3 admits tagged DAG/explicit data and rebuilds private indexes. |
| Legacy identity promotion | `api/session_compat.py` rewrites declared constructor fields and rebuilds collection indexes; private DP state is not persisted or migrated as public fields. |
| Web / derived boundary | `run-analysis.js:379/417` includes `lossless-graph-v1` in derived identity and requires current typed resources for helper hits; `session-request.js` admits typed 1/2/3 for supported saved artifacts. Old derived admission is distinct from current-helper reuse. |
| Other namespaces | Session 42, request 7, raw 4, derived envelope 3, manifest 2, bindings 2 and catalog 3 retain their versions. Schema 3 belongs to typed analysis resources, not a new envelope migration. |

The path owner, codec, Session promotion, native derived owner, metadata owner and
Web descriptor owner remain byte-identical to S06. Five path/public-result
definitions in the otherwise changed protein analysis file are AST-identical.
Later scientific-owner edits are covered by S08 final evidence, not claimed
whole-file identical to S06. All 208 current Python sources match RW-01/02's
verified S08 snapshot; no new native scientific run is necessary for this writer fix.

Public compatibility is backed by Git objects, not only current fixtures:
typed-1 BGC artifact at main first-parent `6b89c781` (also `17e2c9de`), typed-2
at `3bca0e8d`, and path/edge exports, the result class accepted by public diagram
options, and SVG path IDs at release `0.13.0`.
The audit checks original artifact/resource hashes and first-parent membership.
Use the positive publication witnesses; no claim that `0.13.0` was the earliest
tag containing path classes is necessary. S02's known v39 fixture-provenance
mismatch is retained: use the actual Git object as its publication witness,
not the fixture README's exact-byte claim. Historical reports are not rewritten.

The [S06 exception declaration](S06.md#architecture-exception-evidence-author-declaration-review-pending)
still describes the current changed scope and applies unchanged:

- OE-PATH/API/TYPED/METADATA/NATIVE-DERIVED/WEB-DERIVED/LEGACY-IDS/WEB-REQUEST:
  the eight complete before/after sets remain one required owner each;
  total **OE 0→0, delta 0**.
- PE-INFER/EXHAUSTIVE/TYPED/METADATA/NATIVE-DERIVED/WEB-DERIVED/LEGACY-IDS/WEB-REQUEST:
  the eight canonical paths remain the declared one-for-one transitions;
  total **PE 0→0, delta 0**. The normal eager DFS/global-sort path is removed.
- CB-TYPED: `{TYPED-1-READ}` → `{TYPED-1-READ,TYPED-2-READ}`, **1→2, +1**.
  CB-PUBLIC: `{}` → `{PUBLIC-TUPLE-INPUT,PUBLIC-TUPLE-OUTPUT}`, **0→2, +2**.
  Total **CB 1→4, +3**. One decoder/explicit adapter owns old typed input; the
  exhaustive adapter uses the same representation owner, not another inference.
- `ortholog_paths.py` replaces the eager owner; no superseded compatibility path
  is removed. Rejecting old typed/tuple values would violate PATH-B. Retention
  has no automatic expiry; removal needs separately approved compatibility retirement.

These are S06 changed-scope sets, not repository-wide totals. S08/RW/S09 do not
add another compatibility path. Human review must accept the complete S06
exception packet on the eventual exact candidate commit. An agent cannot supply
the maintainer decision, accepted risk, approval permalink or reviewed head now.
S06's independent performance acceptance remains unestablished, as recorded.

## S08 writer and RW connection Review materials

`serializeLosatCache` projects displayed and dormant current entries, builds its
own validated index, and checks each protein raw's instance, directional binding,
protein-set and strict TSV references. It excludes obsolete bindings, keeps valid
inactive settings and nucleotide entries, and releases its index in `finally`.
The source-replacement regression keeps old/current/dormant/removed-record cases
and empty raw as a normal hit. It does not clear live state or History.
`validateSessionLosatArtifacts` independently validates the manifest and every
raw/derived entry during fresh Load; S09 does not remove that boundary.

| Retained contract | Review evidence |
| --- | --- |
| RW-01 | Getter has one reference-validation/TSV scan, retaining context, program/outfmt/args, binding, runtime-ID and strict 12-column numeric/finite rejection. Direct calls still validate; empty text is a hit. |
| RW-02 | Generate deep-clones merged manifest contents, lends the index only during the prepared-pair loop (`run-analysis.js:3714–3829`), clones Worker payloads, and releases on success/error/Cancel before search/render/publication. Index never enters Session, History, retry evidence or shared state. |
| RW-03 | Merge validates R inputs and the merged value once; Generate's invalid-input reload message precedes conflicts. Default helper retains interleaved conflict precedence. Its internal proof is established by validation, not a caller skip flag. |
| RW-04 | Each visited alias uses NFC/trim once per validation, with blank/type rejection, Unicode grouping, stable ordinal checks and input immutability. Separate admission boundaries still validate. |

The Node gates execute the inherited production-loop/count/lifetime probes and
Session reader/writer regressions. RW-04 source/package browser evidence remains
valid for unchanged alias/Generate fault isolation and retry. New S09 browser
checks traverse the same validators through Save, fresh Load and real Worker
saved-raw reanalysis. Normal repeated Generate uses committed typed artifacts;
the archived explicit Run LOSAT continuation separately proves derived reuse.

## Verification, reuse and limits

[Commands](data/s09-premerge-review/commands.jsonl) retain argument vectors,
working directory, failures and exit codes. Durations are operational logs,
not performance samples. [Reuse map](data/s09-premerge-review/reuse.json) records
source, fixture, boundary and artifact-hash correspondence.

- Final focused Node gate: eight existing scripts passed (export, cache,
  Generate simple/derived, settings, stable identities, settings-only, request).
  One nonexistent extra filename in that command was ignored by Node; it is not
  counted as a test. The real Session raw/derived reader script passed separately.
  Schema gate's three checks passed after the affected schema check was rerun
  with child-process permissions; the already-passing reader was not repeated.
- Architecture against latest dev: **Gate PASS / Review REQUIRED**, with all four
  registered owner/path rules conforming. S09 does not change Product authority,
  guard code or mapped contract evidence. Review is not cleared by Gate success.
- Python runner/evidence lint and local asset inventory pass. No runtime URL,
  CSP, dependency, Worker protocol or bundling contract changed. Wheel is rebuilt
  from current source; package installed freshly into `.venv/s09-package`.
- Source browser: desktop 1280×720 and mobile 390×844, fresh contexts, external
  requests blocked from navigation. Four invalid-Save cases per viewport preserve
  Result/selection/History/cache/manifest and create no download. Recovery saves
  13 raw entries. Fresh Load constructs no Worker; real regeneration has 13 hits,
  zero searches and matching full SVG, geometry and provenance. Repeat retains
  one Worker; disposal leaves zero. No external requests or page errors.
- [Source browser](data/s09-premerge-review/browser-source.json.gz) and
  [installed-package browser](data/s09-premerge-review/browser-package.json.gz)
  both pass the same desktop/mobile checks: eight contexts total, 16 invalid-Save
  cases and four fresh-load saved-raw regenerations. Chromium 149.0.7827.55.
  Final audit verifies all 208 wheel Python files, 466 installed package files,
  inherited evidence, and unchanged HEAD/status/tracked diffs in 11 protected
  worktrees. Current source/test hashes and separate diff categories are recorded.
- Reuse S08's 5,972-pass native/full gate with 17 skips/11 deselections, typed/API/CLI
  and 16 read-only SVG comparisons; no unchanged science, full S08 lifecycle or
  native timing rerun. Reuse raw source-replacement/History/Cancel/clear and legacy
  promotion evidence with the RW-03/04 validation coverage and their original limits.
- No new raw searches, timing samples, warmups, profiles, reference regeneration,
  hosted deployment bundle or screenshots. Narrow desktop Chromium viewports are
  not physical mobile devices. RW-01–04 seconds improvements remain unmeasured.

Initial Git sequencer creation, Node child-process checks and local browser
socket creation hit sandbox restrictions. The affected operations passed with
required escalation. No assertion, timeout or expected scientific result was
weakened. No external process was stopped.
The first S09 audit incorrectly expected `OrthogroupResult` to be re-exported
from the release's `api/__init__.py`. The actual release exports path/edge there
and accepts the result through `api/options.py`; the result is defined in the
analysis module. The corrected audit verifies these exact separate witnesses.
The failed audit log is retained; this evidence correction changes no runtime.

## Separate diff audit, rollback and next operation

[Final audit](data/s09-premerge-review/audit.json) and the evidence inventory
separate the ten committed dependencies, 37-file RW-04 inheritance, and S09:

- Production: one S09 owner, `config.js`, the explicit invalid-manifest error.
  Inherited `losat-cache.js` alias change is RW-04, not attributed to S09.
- Tests/tools: S09 extends actual export assertions and adds a bounded mode to
  the existing S08 runner/server/observer/import/Generate/save/dispose machinery.
  Previous modes remain intact. Inherited alias tests/fault option remain RW-04.
  [S09 test patch](data/s09-premerge-review/new-tests-and-tools.patch).
- Documentation: this result plus an additive S09 section in MASTER_PLAN only.
  Past results, scripts, audits and measurement outputs are preserved.
  [S09 documentation patch](data/s09-premerge-review/new-documentation.patch).
- Evidence: only `results/data/s09-premerge-review/`, including authored audit/
  command/preflight files and generated logs/reports. Generated packages and
  private Session files remain ignored, not source or public examples.

Reverse `new-production.patch` to roll back S09 alone, with its tests/docs.
That reintroduces S09-F01 but preserves S08 raw filtering and RW-01–04. Reverse
the separate RW-04 production hunk if explicitly withdrawing that change;
do not restore whole owners from an older S08 commit. S06 rollback is a coupled
representation/consumer rollback: an older reader cannot read schema-3 resources
already saved; preserve old inputs and retain/export readable data before such a
rollback. Rebuild generated packages after any source rollback.

Next: review this source-plus-diff candidate and, if accepted, fix its complete
content in a local work-branch commit. Bind the existing verification to those
unchanged bytes and run the required actual-head integration gate/CI. The human
architecture owner must review S06's full exception sets, publication evidence,
limits and removal condition and supply a separate exact-head decision under
the ratchet. The human reviewer must also review S08 filtering plus S09 invalid
Save failure behavior. Any later commit invalidates an exact-head exception
approval. No unresolved new Product outcome requires a choice packet here.

Proposed commit title: **Reject invalid session manifests and record pre-merge review**

Summary: Fail Session saves explicitly when protein identity metadata is invalid,
while preserving valid empty manifests, raw filtering and live state. Record S06
compatibility and S08/RW review evidence, with offline source/package round trips
and the remaining human approval requirements.

No new implementation commit, push, PR, integration merge, tag, deployment or
external message was performed. Inherited commits retain their original identities.
