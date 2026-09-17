# S01 — authority, consumers, and measurement boundaries

Base: `9a4f7e29ab1b99676bb783321f8a9f40f969d04c`. This inventory records
observations on the fetched base; it does not select PATH-A or PATH-B.

## Authority and existing implementation

| Concern | Authority on base | Observed owners / verification |
|---|---|---|
| Candidate and member limits | OIPC revision 5, PD-OI-001/002/004; OIPC-C01/C04/C05 | `losat-settings.js`, `losat-normalization.js`, session projection and typed decoder. `losat-settings.test.mjs`, `test_session_request_codec.py`, `test_api_request_render.py`. Finite values and `None` are distinct from omitted fresh defaults. |
| Complete records, actual searched database | PD-OI-018, OIC-015 | `linear-sources.js::prepareLosatSourceBatches`, `splitLosatSourceResult`; `run-analysis.js` builds directional metadata with the batch `searchContext`. `linear-sources.test.mjs` verifies 8 records in sources of 6+2: 64 directions / 4 jobs ON, 56 directions / 34 jobs OFF. These are planned source jobs, not measured LOSAT process launches. |
| Optional Collinear inference | PD-OI-021, OIC-018; current compatibility reference | `65f231af` is an ancestor of base. Fresh/reset Web OFF, explicit ON/OFF preserved, missing historical value ON. Real Python `build_orthogroup_collinearity_blocks_from_hits(..., infer_orthogroups=False)` uses direct block construction and member filtering; no replacement inference stub in this runner. Existing Python tests cover this path. |
| Completed raw retry | PD-OI-022, OIC-019; OIPC-C07 | `7db0539a` is an ancestor. `run-analysis.js::completedLosatSearch`, `getReusableLosatCacheEntry`, `retainCompletedLosatSearch`. Owner identity invalidates retry on Clear Cache or Session/History replacement. Raw retention occurs after completed search; last committed Result remains transactional. Existing simple-path tests exercise rollback and retry. |
| Legacy promotion retry | Same failure/identity contracts | New since plan: `9a4f7e29` / PR #534 keeps `workingLegacyProteinRawCandidates` consumption local until successful commit. No new Product choice or parallel migration is needed. |
| Batch key validation | OIPC-C04; PD-OI-018/022 | `47f5cebd` is an ancestor. `build_protein_losat_cache_keys_json` validates once, then `build_protein_losat_pair_identity` accepts the typed manifest. Native probe counts this helper's full validation, not all Generate validation. Browser helper returns `{requestId, result}`. |
| Scientific meanings | `docs/REFERENCE/comparison-programs-thresholds-and-results.md`; MASTER §3 | Similarity groups are retained protein relationship groups, not phylogenetic orthogroups. Search scope and displayed links differ. Preserve self/reverse evidence according to inference, HSP union, member roles, diagnostic evidence, IDs and ordering. |
| Persisted formats | `docs/REFERENCE/session-and-request-compatibility.md` | Current writers: Session **42**, canonical request **7**, protein raw **4**, derived **3**, manifest **2**, Web file bindings **2**, typed JSON resource **2**, feature catalog **3**. Session readers accept 27–33 and 39–42; request readers accept 1/2/5/6/7. Verify each namespace separately in S02. |

Authority links: [OIPC](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md),
[Product ratchet](../../PRODUCT_IMPACT_RATCHET.md),
[architecture ratchet](../../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md),
[current compatibility](../../../REFERENCE/session-and-request-compatibility.md).
The detailed `docs/SESSION_COMPATIBILITY.md` and `gbdraw/web/CLAUDE.md` still
contain older writer numbers; the concise current reference and executable
schema owners agree on 42/7. Do not create a reader for an old number merely
because an overview lists it. `origin/main` at `8675f33b` includes the fetched
dev through PR #535. S02 must still trace exact fields to first-parent/release
history and a positive fixture before proposing compatibility changes.
No `BD-###` is invented or required for these OIPC records.

### Developer preflight

S01 is `IMPLEMENT_EXISTING_AUTHORITY`: evidence and test tooling preserve the
runtime, scientific outputs, defaults, cache lifetime, and every saved format.
There is no registered runtime Product Impact delta, no unresolved outcome
needed to produce the baseline, and no architecture exception. `NOT_ALLOWED`
does not apply. Further path-contract investigation is S02 evidence work;
PATH-B remains unselected, so S01 does not implement its runtime. This is not
an approval of the existing exponential representation.

## Input reconstruction

The runner reads tracked compressed Gallery Sessions and embedded GenBank/raw
bytes. It rejects unexpected record cardinality/selectors/regions instead of
silently taking a first record. It rebuilds CDS extraction with saved record
instance keys and verifies runtime IDs, runtime binding hashes, protein-set
hashes, and every saved raw reference with the current Python identity owner.
File names, resource IDs, and biological records are counted separately.

| Saved input | GenBank resources / records / proteins | Raw directional tables | Notes |
|---|---|---:|---|
| `hepatoplasmataceae_collinear` | 5 / 5 / 2,829 | 13 | 98,947 raw rows; 13,720 pass benchmark thresholds. |
| `hepatoplasmataceae_orthogroup` | 5 / 5 / 2,829 | 25 | 183,661 raw rows; 21,935 pass benchmark thresholds. |
| `vibrio-harveyi-group-collinear` | 11 / 11 / 24,027 | 59 | 657,743 raw rows; saved candidate argument is `--max-target-seqs 5`. Used for manifest timing, not a new search. |

These saved resources each contain one resolved record. Their resource count
does not prove how the biological records were grouped into original uploaded
files when the raw search ran. Historical process/job counts are explicitly
`null`; the runner does not relabel old per-record evidence as source-job timing.
The old report's 31-direction Vibrio scenario is a different selection from
the 59 stored tables used here. Source batching is separately exercised through
the real JS planner with two shared files, ON and OFF.

The two native post-search Gallery cases deliberately fix bitscore=50,
evalue=0.01, identity=0, alignment length=0, member max hits=5, inference=ON.
Collinear uses adjacent scope, default lossless parameters (recorded in JSON),
unit `auto`, anchor `rbh`, two related links; Similarity uses all saved evidence
and adjacent display projection. These settings recreate the report workload;
they are not fresh Web defaults. `render-gallery` instead replays the **complete
saved committed recipes**, including saved comparison resources, appearance and
annotations. A saved editable `config.losat.blastp.mode` can differ from committed
render intent, so the two measurements must not be added as one workflow.

All input hashes and saved search settings are in machine-readable reports.
The old untracked `docs/internal/report.md` at the shared checkout was read as
historical context. Its missing `/tmp/gbdraw-mode-audit/audit.py`, times, test
counts and Node import failure are not evidence for this base.

## S02 consumer inventory

| Boundary / owner | Actual consumption and expansion | S02 question |
|---|---|---|
| `protein_colinearity.py::_build_ortholog_paths` | Builds outgoing/incoming maps for `rbh`/`coortholog`; walks all start-to-end paths. Cycle guard appends the current partial path. Deduplicates by protein tuple, chooses lexicographically smallest edge tuple, sorts paths, allocates `group.path_N`, counts shared proteins, assigns each edge the first containing path ID. | Prove reachable DAG structure before changing representation; preserve current cycle behavior until authorized. Exact count alone does not define legacy path-ID rank. |
| `OrthologPath`, `OrthologEdge`, `OrthogroupResult` | Frozen dataclasses; path protein/edge/shared tuples, edge `path_id`, maps `ortholog_paths_by_orthogroup_id` and `ortholog_edges_by_orthogroup_id`. `OrthogroupEdgeSelectionResult` also contains all edge tables and adjacent display tables. | Tuple types, ordering, empty values and IDs are observable. Do not silently substitute generators. |
| Public Python API / assembly | `gbdraw.api` exports `OrthologPath` and `OrthologEdge`; typed options accept `OrthogroupResult`; `api/diagram.py` and public protein selectors pass it to drawing. CLI and typed requests converge here. | Identify actual consumers requiring exhaustive access, and the supported API migration outcome. |
| Protein-to-genomic projection | `_edge_metadata_for_protein_pair`, `_orthogroup_member_count`, `convert_pair_protein_hits_to_genomic_links` search group/edge members and emit `ortholog_path_id`, RBH IDs, roles, representatives, member counts. | Any reverse index must preserve first-match and many-to-many ordering. |
| Collinear analysis | `collinearity.py::_orthogroup_edge_metadata_for_anchor`, `build_orthogroup_collinearity_blocks_from_hits`, comparison converters retain edge path IDs in anchors and block metadata. ON currently asks the selector for adjacent display tables it does not consume. `comparison_pairs=()` exists at selector boundary but is not used here. | Separate inference from display generation without changing anchors, IDs or block merging. |
| Typed resource | `session_request_codec.py::_TYPED_TREE_CLASSES`, `_encode_typed_tree`, `_decode_typed_tree`, `encode_canonical_typed_resource` recursively encode **every dataclass field and sequence item**; readers enforce expected types and fields. | A compact catalog downstream cannot remove upstream typed JSON expansion. Resource schema has its own compatibility history. |
| Web helper | `python-helpers.js::convert_losatp_blastp_pairs_to_genomic_payload` creates genomic rows, typed `orthogroupResult` / `collinearityResult`, serialized groups and provenance; stores JSON in Python converted cache. | Trace all consumers and JS derived retention before eliminating duplicated completed-payload ownership. |
| SVG | `render/groups/linear/pairwise_match.py` maps link `ortholog_path_id` to `data-ortholog-path-id`. `web_support/orthogroup_metadata.py::serialize_orthogroups_payload` materializes `orthologPaths` with `proteinIds`, `edgeIds`, `sharedProteinIds`; render metadata carries this to interactivity. | Preserve required semantic attributes and biological identity; visual geometry alone is insufficient. |
| Feature catalog | `web_support/feature_catalog.py::_ORTHOGROUP_COLLECTION_COUNTS` reduces arrays to `orthologPathCount` / edge counts in catalog schema 3. | Count reduction is currently after full path construction and serialization. JS exact integers above 2^53−1 need an explicit representation decision. |
| Popup and editor | `app/pairwise-match-popup.js` reads `orthologPathCount`, falling back to array length; group details display counts and member metadata. `app/orthogroups.js` uses group identity and naming/edit state. | No evidence found here requiring an exhaustive path browsing UI; absence is not permission to retire public arrays elsewhere. |
| Session/raw/derived/History | `session_io.py`, `session_request_codec.py`, JS session/cache owners validate resource and runtime bindings; committed typed resources, editor catalog and derived payloads can contain path metadata. Raw rows have no path array. | Follow each namespace and load/replay path independently. No S01 writer, migration or compatibility path is added. |

The native oracle hashes the complete ordered dataclass result and DataFrames,
including dtypes, index, tuple/list distinction and exceptional values. It also
hashes complete metadata and typed resources. R=8/12/16 expands and checks all
paths; R=24/32/56 is **formula only**, stored as decimal strings. The formula
applies to the complete forward ordered-record DAG, not arbitrary graphs.

## S03–S07 baseline boundaries and counters

| Stage | Reproducible input / measured operations | Scope and limitations |
|---|---|---|
| HSP | Seeded 1,000 pairs × 3 HSPs, randomized row order. Separate overlap, disjoint, reverse, clamp, duplicate/tie, missing/unknown ID, NaN/Infinity and empty cases. Aggregate calls, input rows and representative-rank calls. | `_aggregate_hsps_by_protein_pair` only; filtering/member selection are separate. Infinity alignment length currently raises `OverflowError`; observation is not independent authority. |
| Sparse support | Core/unassigned 200/200, 400/400, 800/800 with only reciprocal core evidence. Actual `_build_core_support_candidate` calls. A separate incoming-only/local/domain/tie fixture and dense 24-core/24-unassigned fixture have complete result oracles. | Sparse candidates currently grow 40k/160k/640k; dense tests do not promise sublinear dense inference. Existing detailed tests cover assignment reasons, snapshots and record-local competition. |
| Cache | Actual helper shared LRU, actual native parse/filter, 49/64/81 tables plus a converted payload; cold/warm parse counts and retained DataFrame deep bytes. | Retention bytes exclude Python containers/keys, strings, old snapshots and Wasm heap. This is an LRU boundary probe, not full converted-helper timing or a new byte-budget design. Cache telemetry is excluded from the scientific digest so an improvement can compare equal. |
| Paths | Public selector on complete record graph R=8/12/16. All group/member/edge/path metadata hashed; summary counts checked against formula. | No R≥24 enumeration. R=56 theoretical count exceeds JS safe integer range. No compact algorithm yet. |
| Merge | Chain clusters interrupted by isolated anchors: 300/600/1200; conflict function calls and supplied-anchor counts. Separate strict-conflict and reverse examples. | Input counts are measured supplied work, not an automatic proof of future scanned work after early-exit/indexing. S07 should replace counters when the owner changes. |
| Gallery | Parse, filter, aggregate, complete post-search, group metadata, typed resource serialization; full saved typed replay/SVG is independent. | Stages are nested/overlapping: post-search includes aggregation. Do not sum stage times. Inputs to each stage are prepared outside its timed interval. |
| Browser | Production diagram Worker batch-key helper, cold call plus seven warm calls; actual `Worker.postMessage` observation and key equality with native. Real source planner ON/OFF. | Round trip includes transfer, queue, JSON and Pyodide. Synchronous `postMessage` measures enqueue/clone only. No raw search timing or separate pure-Pyodide/transfer estimate. |

## Cache and lifecycle read-set for S05

Raw identity owns protein content/bindings, record direction, program/outfmt,
normalized search args and actual source database scope. Completed retry is
in `run-analysis.js`, separate from the diagram Worker's intermediate cache.
Filtered Python tables depend on raw key and thresholds. Converted Python JSON
also depends on mode, member settings, Collinear settings, record protein cache
keys, view length/reverse, pair order, display flags and raw keys. JS derived
cache has a completed payload owner too; inspect its consumers before adding
retention. `_protein_sort_key` reads record/coordinates, so an inferred cache
cannot casually omit display/order/annotation inputs.

The diagram Worker lifecycle owns the helper namespace; terminating/recreating
it discards those two Python caches. Parsed biological input reuse is separately
owned by `PreparedBiologicalInputCache`, with success-only publication. Do not
merge raw retry, intermediate analysis, prepared inputs, and final Results into
one undocumented cache contract. The existing browser acceptance exercises
actual helper/render/migration/save/reload; Node tests exercise cancel/rollback
orchestration with controlled boundaries. Neither proves the future S05
byte budget, large-input peak memory, or every adversarial browser sequence.

The current `clearLosatCache` action explicitly clears completed retry, JS
raw/derived maps, manifest and legacy candidates. It contains no call to clear
the Python helper stores or terminate the diagram Worker. Do not infer a
Python-cache clearing guarantee from that JS action. S05 must explicitly trace
invalidation for Clear Cache, input and Session/History replacement, stale
completion and Worker recreation when designing its new retention owner.

## Ownership and rollback

Before/after production owners and paths are identical. No production hooks,
fallbacks, dependencies, caches, migrations, public figures, or reference
outputs are changed. Instrumentation monkeypatches exist only during a native
probe call and are restored by `ExitStack`; timing and memory runs use actual
unpatched functions. Browser instrumentation observes transport only. The
runner is one CLI, with one private browser adapter. `benchmark_diagram_layout.py`
remains focused on layout; its measurement conventions informed this runner.
Rollback consists of removing this session's tooling, test and evidence files.
