# Worker resource-cache recovery: 05A4-11

Implementation evidence, based on `dev` at
`c3573a62ae50af367cb2b94810d426745fcd81af`. This change addresses only
05A4-11. It does not establish new Product authority.

## Original reproduction and diagnosis

The retained SESSION 05A4 M23 harness was run before production edits in a
clean worktree. In Linear mode, upload the matching lambda GFF3 and FASTA,
Generate A, replace FASTA with the retained mismatch, Generate B, restore the
original FASTA, and Generate once. The first recovery failed; a second
Generate succeeded after constructing another Worker.

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| `gbdraw/web/tutorial-data/lambda-gff3/NC_001416.gff3` | 36794 | `d53e05de87933104cd26111bca42006cce9b5e903fb5b187740f963b3a2098cb` |
| `gbdraw/web/tutorial-data/lambda-gff3/NC_001416.fna` | 49253 | `80897a7ee6b8aaffbab5442e0daad292592ac74701dbdf35af4b400ae0770ef3` |
| Retained `mismatch.fasta` | 424 | `ac7af7ce9b3d2ad76812b042f0befb3902e050ff0c8788f147a9247a16537d16` |

B is rejection class **C**, a structured Python render error. The Worker
returns `type: run, ok: true, results.error.type: ParseError`, reporting no
matching FASTA record for GFF record `NC_001416.1`. Resource staging and Python
workspace cleanup both complete. The Worker remains alive.

Before B, both caches know `record-1-gff3` as `render-resource-1` and
`record-1-fasta` as `render-resource-2`. B reuses GFF3, transfers 424 FASTA bytes
under `render-resource-3`, and deletes the old FASTA cache file. The main cache
still knows token 2 because artifact finalization never runs for B.

The first restored A therefore omits both byte payloads and requests FASTA
token 2 from a Worker holding token 3. Its exact error is
`Render resource cache miss for 'record-1-fasta'.` This is class **D**, a Worker
staging/protocol error, and terminates the Worker. The second recovery
**does restage the same source bytes**: it transfers all 86,047 bytes to a new
Worker. This is a cache transaction boundary defect, independent of FASTA.

## Owners and transaction boundary

`diagram-resource-staging.js` remains the main transport owner. `prepare()`
compares payload identity and size against `cachedResources`, reuses a token
only on an exact hit, materializes misses into `stagedResources`, and constructs
`nextCachedResources`. Preparation alone publishes nothing. `commit()` replaces
the cache with exactly that manifest; a reset or intervening commit invalidates
an older preparation. `reset()` discards all Worker-specific knowledge.

`diagram-generation-worker.js` remains the raw-byte owner. Its existing
`stageRenderResources()` validates the manifest, creates/replaces cache files,
links them into the request workspace, and prunes resources absent from the
manifest. It then renders and removes the request workspace. Raw cached files
remain until replaced, pruned, or the Worker is terminated.

`diagram-generation.js` now commits main cache knowledge when the existing
successful `run` envelope acknowledges completed Worker staging. This includes
structured Python input errors. The client keeps one active render until that
response arrives, so another render cannot use cache knowledge while staging
is incomplete. Fatal Worker responses, runtime/message errors during rendering,
and cancellation terminate the Worker and reset transport knowledge.

`run-analysis.js` continues to own generated-artifact candidates, Result
admission, readiness, activation, History, finalization, and rollback. The
deferred `finalizeResourcePromotion` callback is removed from both the response
and artifact transaction. A rejected Result can preserve artifact A while both
raw resource caches coherently describe B. Restoring A transfers only bytes
that the Worker no longer owns.

There is **no new Worker message or field**. `run.ok` acknowledges transport
completion; `results.error` independently describes the Python outcome.
No second cache, token namespace, retry, or recovery watcher is introduced.

## Settled failure states

| Class | Main raw-cache knowledge | Worker and cleanup |
| --- | --- | --- |
| A: validation before dispatch | Prior acknowledged manifest | Unchanged; no render workspace |
| B: structured helper/extraction error | Prior acknowledged render manifest | Reused; auxiliary operation does not stage render resources |
| C: structured Python render error | Completed request manifest | Reused; request workspace removed, raw files retained/pruned by manifest |
| D: staging/protocol/runtime failure | Empty after reset | Terminated, including partially staged bytes |
| E: post-render Result admission rejection | Completed request manifest | Reused; prior accepted artifact restored independently |
| F: cancellation or stale response | Empty on termination; otherwise last acknowledged manifest | Cancellation releases Worker; a late response cannot republish disposed ownership |

The invariant is that every omitted payload has the same resource ID, token,
size, and bytes in the live Worker. Knowledge is published only after Worker
acknowledgement. Failed preparation cannot publish; incomplete Worker staging
cannot produce that acknowledgement.

Parsed/prepared Python caches retain their existing transaction rules. Worker
manifest tokens are validated against resource paths and byte sizes in
`web_support/request_render.py`. Changed tokens invalidate dependent parsed,
resolved-record, and interactive-context state. Failed prepared fills remain
unpublished. No Python prepared-cache behavior changes here.

## Recovery metrics

| Operation | Before: transferred bytes / hits | After: transferred bytes / hits |
| --- | --- | --- |
| Original A | 86047 / 0 | 86047 / 0 |
| Rejected B | 424 / 1 | 424 / 1 |
| First restored A | 0 / 2 (false FASTA hit; fails) | 49253 / 1 (succeeds) |
| Second restored A | 86047 / 0 (new Worker) | Not needed |
| Unchanged A after recovery | Not part of original M23 | 0 / 2 |

The original recovery journey changes from two Worker constructions and
initializations to one. Recovery restages FASTA under token 4 and reuses GFF3
token 1. The browser regression also records materialization, base64-decode,
decoded-byte, and Worker lifecycle metrics through the existing test hooks.

Across the original journey through eventual success, transferred bytes are
172,518 before and 135,724 after; materializations are 5 and 4. Both runs record
2 base64 decodes totaling 86,047 bytes. The old run reports 3 cache hits,
including the false FASTA hit; the fixed three-Generate journey reports 2 real
GFF3 hits. The packaged application at 390 × 844 produces the same fixed counts.

## Verification and review

- `tests/web/diagram-resource-staging.test.mjs`: uncommitted preparation,
  acknowledged replacement, unchanged hits, changed payloads, and reset.
- `tests/web/diagram-resource-recovery.test.mjs`: real private Worker staging
  with a deterministic filesystem and the real client/transport; W2–W9,
  partial staging failure, validation/helper/extraction rejection, and
  cancellation during response handling.
- `tests/web/diagram-resource-recovery.playwright.spec.js`: exact retained M23
  input journey, first-attempt recovery, warm reuse, changed-source admission
  rejection, cancellation after staging, Result/mounted/export agreement,
  recovered canonical request, and no external requests, page errors, or
  unhandled rejections. The full functional configuration discovers this case;
  it has no PR-smoke marker and does not increase the ten-case budget.
- Existing Python prepared-cache controls verify changed tokens invalidate
  every dependent layer and failed fills preserve prior entries.

Local gates passed: the complete JavaScript unit suite (553 tests), the final
resource-owner controls (9 recovery and 2 transport tests), 79 Python/package
tests, all 6 existing offline GUI browser contracts, and the dedicated recovery
browser test. All 14 focused browser neighbors pass: generated request History,
composite resource identity and CLI sidecars, source visibility, source legend
reconciliation, and active Result/mounted/export authority. The retained M23
harness passes on both the source tree at
1600 × 1000 and an independently rebuilt wheel at 390 × 844. The package's
three changed production files were byte-compared with the worktree.

Negative controls retain the original failing first retry. Running the new
owner tests against archived base source, with base-style successful-artifact
finalization, rejects the failed-request cache state and reproduces a cache
miss for a pruned overlapping resource. The fixed source passes these controls.
Sandbox Git subprocess and dependency-download restrictions were rerun with
the required access; they are not reported as application regressions.

Product preflight: **IMPLEMENT_EXISTING_AUTHORITY**. Base-branch
`OPTION_INTEGRITY_PRODUCT_CONTRACT.md` OIPC-C04 requires correct cache identity;
OIPC-C07 and PD-OI-016 require preservation of the last successful artifact and
the ability to correct and retry. No alternative Product outcome or retirement
is selected, so neither an evidence-dependent Product choice nor a new Product
decision is required. The implementation respects the non-waivable cache
identity and failure-isolation rules.

This is **architecture-bearing**, with normal human review required. The same
main transport, Worker client, raw Worker cache, and artifact owner remain.
The artifact-to-resource promotion path is removed; resource publication has
one client acknowledgement path. Owner/path excess does not increase; no
compatibility path is added. There is no architecture exception, dependency,
privileged operator, or persisted schema change. Session 41, canonical request
7, feature catalog 3, and Web file binding 2 are unchanged. Reverting this
commit restores the previous transaction timing without migration.
