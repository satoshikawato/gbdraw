# S02 — path semantics, consumers, and publication evidence

Evidence base: `9a4f7e29ab1b99676bb783321f8a9f40f969d04c`; fetched main:
`8675f33bb2f2ac8f174f65e9a36529df942c5333`. This is an observation inventory,
not Product authority. [Design](S02_PATH_CONTRACT.md),
[decision](PATH_DECISION_PACK.md), [handoff](S02.md).

Follow-up: the complete PATH-B receipt was received on 2026-09-15; see
[decision handoff](S02_PRODUCT_DECISION.md). The authority search and preflight
below remain the observations on the named base before that receipt, not a
claim that Product selection is still missing.

## 1. Exact current enumeration

Owner: `gbdraw/analysis/protein_colinearity.py::_build_ortholog_paths`.
Its only production call is `_build_anchor_core_orthogroups`.

| Step | Observed rule |
|---|---|
| Group order | Iterate the input mapping in insertion order; return both maps with those keys, including groups with no paths. |
| Participating edges | Only `rbh` / `coortholog`, with **both IDs present in protein_map**. `render_role` is irrelevant: an assigned `display_edge` coortholog participates. Related, same-record inparalog and record-local paralog edges do not. |
| Edge identity | `_edge_id(e)` is `group:query_record:query_id->subject_record:subject_id:kind`. It excludes scores, role and existing path ID. |
| Node order | `K(p) = (int(record_index), int(start), int(end), str(protein_id))`. Not feature index, strand, edge direction or string-only protein order. |
| Starts | Nodes appearing as endpoints, with no incoming participating edge. If there are **no such nodes anywhere in the group**, use every endpoint node. Sort starts by K. There is no isolated-protein path. |
| Traversal | Sort each outgoing list by `(edge.subject_record_index, K(subject), edge_id)`. Follow stored direction; never reverse an edge here. |
| Terminal | No outgoing participating edge: append the current walk only if it contains at least one edge. |
| Cycle | For each edge whose subject is already in the current protein tuple, append the **current prefix**, if nonempty in edges, without the closing edge. Continue other outgoing branches. A repeated self edge from an empty prefix emits nothing. |
| Dedup | Key is the complete **protein ID tuple**, not edge tuple, path length, score or unordered member set. For equal protein tuples, keep the lexicographically least **whole edge-ID tuple**. |
| Final order | Sort `(protein_tuple, selected_edge_tuple)` by `(tuple(K(p) for p in protein_tuple), selected_edge_tuple)`. Intermediate traversal order is not the final ordering contract. |
| Path ID | One-based final rank: `f"{group_id}.path_{rank}"`; no padding. IDs are result/order dependent, not a biological persistent identifier. Numeric ranks 2 and 10 must not be compared by their ID strings. |
| Shared | Count each protein once per **deduplicated path**; membership in more than one path makes it shared. Each path stores its shared members sorted by K, **not necessarily in traversal order**. Shared does not mean present in every path, high degree, paralogy or shared across groups. |
| Edge first path | Visit selected edge tuples in final path order, recording the first path ID for each edge-ID string. Replace matching original edges, preserving their original input order and all other fields. Unselected/filtered edges retain their previous `path_id`, including non-null values. Duplicate edge objects with the same edge ID get the same first ID. |
| Empty result | With no participating edges, return `tuple(original_edges)` and `()`. With paths, return actual frozen `OrthologPath` objects and tuple-valued fields; no generator is returned. |

The upstream caller sorts original membership edges by query record, subject
record, query ID, subject ID and kind before this call. `_orthogroup_result_from_member_ids`
copies sequences into tuples again. `OrthogroupResult` is frozen but contains
ordinary mutable dictionaries/lists; “frozen” does not make its entire object
graph immutable.

### Cyclic observations that a generic DAG count misses

The executable test `test_cycle_observations_do_not_become_a_dag_fallback` covers:

- `a→b→a`: two prefixes `(a,b)` and `(b,a)` because all nodes become starts.
- That cycle plus a disconnected `c→d`: **only `(c,d)`**. A source elsewhere
  suppresses the all-nodes fallback for the cyclic component.
- `s→a→b`, `b→a`, `b→t`: both `(s,a,b)` and `(s,a,b,t)` are paths.
  A cycle terminal can be a proper prefix of another returned path.
- `a→a` alone: zero paths.

The prototype reports a cycle explicitly. It does not claim these graphs have
the DAG path language and does not silently drop/reorient edges.

## 2. DAG proof and its boundary

This is a proof about **current inference construction**, not all publicly
constructible `OrthologEdge`/`OrthogroupResult` objects.

1. `_dedupe_anchor_core_directional_rows` removes self hits, unknown endpoints,
   nonpositive scores, and duplicate query/subject directions using the existing
   HSP/member ranking. Repeated HSPs do not become separate path transitions.
2. `_select_anchor_core_edges` skips equal-record pairs, requires reciprocal
   membership support, and uses `_canonical_edge_endpoint_ids`. Both strict RBH
   and near-reciprocal coortholog core edges therefore increase record index.
   `selected_by_pair` keeps one edge for each canonical endpoint pair. The
   core subgraph is a DAG regardless of evidence table order or reverse hits.
3. `_build_anchor_core_orthogroups` takes `core_member_snapshot` **before**
   assigning extra proteins. Each initially unassigned protein can be assigned
   once and receives exactly one membership edge to one snapshot member.
   `_best_evidence_between_protein_and_members` examines both directions and
   preserves the selected raw direction. Hence an added coortholog may decrease
   record index. The new protein is nevertheless a leaf in the undirected
   membership subgraph: no other added member can use it as a snapshot member.
   Adding one leaf edge cannot create a directed cycle.
4. Same-record support has kind `same_record_inparalog`. Record-local clusters
   append `record_local_paralog` edges. Both are excluded from path traversal.
   Related/diagnostic edges live in a separate map, also excluded.
5. Multiple roots/sinks and disconnected groups do not change this proof.
   `comparison_pairs`, display reversal and adjacent/all projection occur around
   the inference and do not reverse its stored path edges. Collinear inference
   OFF does not invoke this path construction.

The named real-selector fixture has core `p0→p1`, added **`p4→p0`**, and
incoming-only membership evidence `p1→p5`. It proves that sorting by record
number is **not** a valid topological algorithm for the whole inferred graph.
Other fixtures produce real same-record membership, a separate record-local
group, near-reciprocal core coorthologs, and duplicate HSP evidence. Sixty seeded
eight-protein directional inputs also pass DAG and complete-path checks.

Canonical extraction gives unique protein IDs, matching map keys, and thus
injective K. The prototype explicitly requires that canonical identity boundary.
Public dataclasses/typed payloads are less constrained: the typed decoder checks
fields and types, not graph topology, ID ranking, completeness of a path corpus,
or consistency of shared flags. A positive test round-trips a cyclic edge set,
one deliberately selected path, a custom ID and supplied shared tuple.
Those legacy objects must be preserved as explicit data, not reconstructed
from the edge closure. Forged/noncanonical protein-map aliases are outside the
prototype's ordering proof; S06 must not silently classify them as canonical
input or broaden input rejection at an existing public boundary.

## 3. Consumer inventory

Legend: **ALL** = exhaustive data required by current output type/format;
**COUNT** = exact total; **ONE** = a particular path; **EDGE** = first path ID
on each edge; **SHARED** = per-path shared tuple; **TRANSPORT** = forwarding or
validating data without using path meaning. These describe actual reads; an
ALL transport consumer need not semantically need every path for rendering.

| Boundary and owners | Required information / observed behavior | Hidden work and S06 consequence |
|---|---|---|
| Public exports `gbdraw/api/__init__.py`; dataclasses in `analysis/protein_colinearity.py` | ALL, ONE, EDGE, SHARED: `OrthologPath`, `OrthologEdge`, `OrthogroupResult`; `OrthogroupEdgeSelectionResult.orthogroups`, `ProteinBlastpResult.orthogroups` expose the complete result. External callers can index, iterate, call `len`, inspect exact tuple types and compare dataclasses. | Tuple→lazy substitution is an observable API change. Repository search cannot inventory third-party callers; publication evidence requires keeping an explicit exhaustive route. |
| Public producers `select_rbh_orthogroup_edges_from_directional_hits`, `build_rbh_orthogroup_protein_blastp_comparisons`; Collinear producers | ALL and EDGE currently constructed before their return. Collinear result embeds the orthogroup result. | One inference owner serves both modes. A new compact mode must reach the producer, not begin after return. |
| Input/assembly `api/options.py::LinearDiagramOptions`, `interface.py::LinearOptions`, `api/diagram.py`, `api/request_render.py`, `diagrams/linear/assemble.py` | TRANSPORT; accept precomputed orthogroups and embed/pass metadata in the typed render flow. | Preserve old typed input objects. Normal render must receive the compact result before encoding or interactive metadata construction. |
| `diagrams/linear/orthogroup_alignment.py` | Members, representative flags, feature identities and group IDs for alignment/label eligibility; **no path reads**. | No reason to enumerate for alignment, `orthogroup_top` labels, names or descriptions. |
| `protein_colinearity.py::_edge_metadata_for_protein_pair`, `_build_adjacent_display_edges_by_pair` and genomic converters | EDGE, member roles/counts/RBH identities. Scan membership edges then related edges; forward/reverse endpoint matching and first matching edge are observable. | Do not use selected graph transitions to replace the complete evidence-edge list or reorder it. These consumers do not request a path object. |
| `collinearity.py::_orthogroup_edge_metadata_for_anchor`, `CollinearityAnchor`, `_joined_anchor_values`, comparison converters | EDGE in anchors and block rows; multiple anchor values are joined using existing first-occurrence dedup/order. | Count alone cannot replace the first path ID. Preserve anchors, unblocked anchors and block metadata. No all-path access. |
| `session_request_codec.py::_TYPED_TREE_CLASSES`, `_encode_typed_tree`, `_decode_typed_tree`, `_read_typed_json_resource` | ALL + TRANSPORT. Dataclass fields recursively encode; sequence comprehensions create JSON lists. Decode creates a list and then a tuple for tuple hints. Schemas 1/2 are validated independently of request schema. | First major exponential materialization boundary; `dataclasses.fields` will traverse a field even if drawing never reads it. New class/field schema must explicitly encode graph data. Private serializer does not encode selection wrapper dataclasses. |
| Web `python-helpers.js::convert_losatp_blastp_pairs_to_genomic_payload` | ALL + TRANSPORT: Similarity emits typed `orthogroupResult`; Collinear emits typed `collinearityResult` (kind `result`) containing orthogroups. | `encode…` → bytes → UTF-8 text → `json.loads` → enclosing `json.dumps` → Python converted-cache JSON string → Worker JSON parsing/transport. **Current helper does not additionally emit a rich `orthogroups` array for these nonempty modes.** This corrects the abbreviated S01 inventory. |
| `run-analysis.js::hasRequiredCanonicalAnalysisResource`, `getLosatDerivedCacheEntry`, `setLosatDerivedCacheEntry` | TRANSPORT/shape guard; expects schema 1/2 and tagged `OrthogroupResult` or `CollinearityResult`. Stores completed payload; stages new cache map until commit. | Returning a lazy object here cannot survive JSON. Native rich derived entries lacking the canonical typed resource are not reusable helper hits merely because their reference validation passes. |
| `services/session-request.js::buildComparisons`; canonical/imported comparison adapters | TRANSPORT of typed artifact or File; `resources.addJson` serializes the entire typed resource again. Descriptor kinds `orthogroupResult` / `collinearityResult` and `canonicalJson` route to Python. | Requires graph-capable shape admission and one canonical normalization; no alternate render pipeline. |
| `web_support/orthogroup_metadata.py::serialize_orthogroups_payload` | ALL, EDGE, SHARED in current public/internal rich group payload. `_serialize_ortholog_path` makes three lists. | Used by interactive context **and** native derived artifacts. Removing a downstream array does not remove this allocation. Must offer a summary projection at this owner. |
| `render/interactive_context.py` → `web_support/feature_catalog.py::_normalized_orthogroups` | COUNT for paths; edge counts, members and group names. Starts from rich arrays, copies members and `_sequence` copies collections before `len`. Nonempty arrays become counts, empty arrays are removed without inserting zero. | No normal catalog path browsing or shared-array consumer. Preserve absent/zero display behavior. Compute count upstream and preserve exact decimal text. |
| `render/groups/linear/pairwise_match.py`; `services/svg-sanitization.js` | EDGE as `data-ortholog-path-id`; sanitization permits that attribute. | Static SVG and popup alignment details retain the same string, including joined block values. No geometry change is needed. |
| `render/interactive_svg.py::enrich_svg` | TRANSPORT of **schema-3 catalog**, not full path arrays; serializes catalog into metadata (`gbdraw-interactive-feature-metadata`, interactive schema 3). | Final SVG has reduced counts already. New graph need not be embedded in each SVG: SVG is not currently the exhaustive export. Preserve semantic path attributes and count. |
| `services/feature-catalog.js`, `orthogroup-feature-metadata.js`, `app/session-feature-metadata.js` | TRANSPORT/project members into feature indexes; validate IDs, references and Result alignment. | Current schema-3 group extra count fields accept strings. `admitFeatureCatalog` cloning/adoption must not expand graph collections. |
| `app/pairwise-match-popup.js` | COUNT in Similarity popup and Collinear group detail rows; EDGE in alignment row. Prefers `orthologPathCount`; falls back to `orthologPaths.length`. | `firstText` handles decimal strings exactly. No `Number` conversion is necessary. No UI enumerates paths or consumes `sharedProteinIds`. |
| `app/orthogroups.js` | Members, annotations, representative identity; names/descriptions, member FASTA download and selection/alignment. | No path iterator/count computation. “Download members” is not “download all paths”. |
| `api/request_render.py::_build_current_derived_entries` | ALL + TRANSPORT via a **second** `serialize_orthogroups_payload` call. Native payload contains rich `orthogroups` and pair/provenance data. | Normal native rendering can allocate arrays independently of the helper. Must switch to summary + compact canonical analysis resource using the S05 result owner. |
| `analysis/protein_artifacts.py::validate_current_derived_protein_artifacts`; JS `app/losat-cache.js::validateDerivedProteinReferences` | TRANSPORT/integrity: recursive full object walk; recognize scalar/array protein keys, shared arrays, compound edge IDs and unit IDs. | Graph nodes/transitions must use explicit recognized reference fields, with complete manifest validation. A graph count does not establish identity validity. |
| `session_io.py`, `session.py`, `api/session_compat.py`, `services/config.js`, session resources/file modules | ALL + TRANSPORT where typed resources or derived payloads contain arrays; base64 embedding, JSON parse/stringify, resource hash checks. Preserve current/draft separation and existing old-artifact admission. | Current catalogs omit old duplicate rich group state. Old Sessions can contain both `orthogroupState.groups` and rich derived groups. Do not introduce a new migration for an absent namespace/version. |
| `services/history-snapshot.js`, History owners | TRANSPORT: generated owner snapshots retain artifact references; general fallback group snapshots use JSON cloning. Mutable intent separately stores group edits. | Do not claim every History operation deep-copies the full payload. Both retained and copied objects must remain compact. Failed/stale/canceled replacement keeps the last committed Result. |
| Raw cache / manifest / bindings / raw TSV export | No path language. Raw has search evidence, manifest binds identities, bindings embed source occurrences. | Leave raw identity, job scope, source grouping, candidate/member limits and raw retry ownership unchanged. |

Search scope: all tracked `gbdraw/` Python and JS; symbols `OrthologPath`,
`OrthologEdge`, `OrthogroupResult`, both snake/camel path fields, shared fields,
edge path columns, typed result kinds, and full-result serializers/transports.
Tests, Gallery resources and reference documents were searched separately.
External callers are unknown; this is not evidence that public arrays are unused.

## 4. Published namespaces, not inferred version chains

Reproduction: `python tests/prototypes/path_publication_history.py --output …`.
[s02-publication.json.gz](data/s02-publication.json.gz) contains **90 first-parent
source snapshots**, every local release tag checked, exact Git SHAs, source
hashes, decoded typed-resource hashes, artifact counts and path-field locations.
Historical artifacts are read directly from Git objects; no historical code is
executed and no Gallery artifact is rewritten.

| Namespace | Positive public evidence | Current handling / decision relevance |
|---|---|---|
| Python path/edge/result API | First-parent main `174fab8c` (#255) has class + `gbdraw.api` export. Release tag `0.13.0` has them too; earlier local tags lack the path class. | Actual frozen tuple contract is public. Release source and the S01 positive full-result oracle establish concrete shape; cannot silently swap field values for generators. |
| SVG edge attribute | `0.13.0` Hepatoplasmataceae Similarity Session has 1,846 `data-ortholog-path-id` occurrences in saved SVG. Current main Gallery also has 1,846 (Collinear: 500). | Attribute continuity is separate from presence of all-path metadata. This release Session itself contains **no** rich path-array field; do not claim it does. |
| Rich group paths / old derived 1 | Main first-parent `6b89c781` contains BGC Session 33/request 2 with `orthogroupState.groups[*].orthologPaths` and derived schema 1. Existing schema-v2 fixture is **byte-identical after decompression**: SHA-256 `3f10b998e7f257fa0458fb4310f9c5ae246d3e431b8f54755b16fe506f53dbdd`. | Arrays contain path IDs, protein/edge lists and shared lists. Legacy ID promotion stays in existing manifest/raw-candidate owners. |
| Typed JSON resource 1 | Writer present in main `10d3a3d2` (#287). Main `6b89c781` BGC artifact has `comparison-orthogroups`, kind `orthogroupResult`, tagged `OrthogroupResult`; resource SHA-256 `40868b747828f2f1f9b68f2599f51214727897a2a0e0d39c00602dbd9e858b37`. | This **is** a published path-bearing typed format, even though the enclosing request schema is 2. |
| Typed JSON resource 2 | Writer main `3bca0e8d` (#319); its BGC positive resource hash `295ee47f12dbeaf5d763436b9e279480b51143989ca892863a36ff3451d3b1dc`. Current main Gallery contains typed resource 2 too, with exact file-byte equality to the S02 base. | Readers 1/2 and class/field tags need explicit compatibility if the graph format advances. No existing typed-resource 3 is asserted. |
| Protein raw 4 / derived 3 / manifest 2 | Main `17e2c9de` (#312) BGC Session 39/request 5 has all three and rich path arrays. Git-object positive artifact is included by hash and location. | Raw 4/manifest 2 do not encode paths. Derived 3 is a real compatibility input for a new representation. Legacy raw 2 and derived 1 also have positives; raw 3/derived 2 are not introduced as migration steps. |
| Feature catalog 3 | Main `8228ffab` writer; main `3bca0e8d` BGC positive catalog and current Gallery positive catalogs. Arrays are reduced to numeric counts when nonempty. | Existing Python and JS validators accept a decimal-string count; S02 tests actual catalog admission and both popup routes with `18014398509481985`. A separate new catalog reader is unnecessary for this accepted field representation. |
| Interactive SVG envelope | First-parent main `136d8def` declares the **string** `gbdraw-interactive-feature-popup-v2`; `8228ffab` declares numeric schema 3, which embeds catalog 3. | Separate namespace from typed resource schema 2 and feature catalog 2. No extra old-SVG reader is proposed for path representation. |
| Session | Source snapshots show writers 28–33, then 39 (`17e2c9de`), 40 (`8228ffab`), 41 (`4e8c9380`), 42 (`3fd50841`). Positive full Sessions above and `settings-only.v42.json.gz`. Version 27 remains an existing accepted legacy input; no new version-27 reader is proposed. | Current writer 42; readers 27–33/39–42. Existing current-source compatibility is retained; next envelope version can mark a graph-bearing Session for old-reader rejection. |
| Canonical request | Main writers 1 (`10d3a3d2`), 2 (`0266acb7`), 5 (`17e2c9de`), 6 (`4ea96685`), 7 (`4e8c9380`). Positive full fixtures cover 2/5/6/7; source records establish 1. | Current writer 7/readers 1/2/5/6/7. Descriptor routing can remain schema 7 while its **independently versioned resource** advances. |
| Web file bindings | Current writer 2; positive `settings-only.v42` contains bindings 2, `single.v41-bindings1` contains 1. | Not path storage. Do not advance bindings or convert source occurrences for this work. |

No local release tag later than `0.13.0` exists after fetch. A main merge
establishes first-parent inclusion, not proof of a hosted deployment or an
untagged release. Session 34–38/request 3–4 are explicitly rejected branch-only
formats under the base compatibility reference; the source snapshot history
does not make them public. No reader for these or raw 3/derived 2/catalog 2 is
proposed. Existing feature-catalog “biological authority” schema 1 is a separate
legacy projection boundary, not justification for a 1→2→3 migration chain.

### Provenance discrepancy found

`tests/fixtures/sessions/README.md` says the v39 fixture is unchanged JSON from
`17e2c9de`. Its decoded hash is indeed the documented `9407365a…`, but the
actual Git-object artifact at that commit hashes to `1a99019c…`. The JSONs also
differ (resources, typed resource presence, request, editor/config/derived data),
so this is not just compression/whitespace. The fixture remains a positive
schema-39 example, **not exact bytes from that cited commit**. S02's publication
claim uses the real main artifact and stores both hashes. The unrelated
historical README/fixture is not rewritten in this session.

## 5. Authority and developer preflight

| Inspected authority | Finding |
|---|---|
| `tools/web-product-impact-map.json` | Mapped render-request/current-Result boundaries; no path-representation concern. S02 has no registered subject delta. General render and Result requirements still apply jointly. |
| `tools/web-product-decisions.json` | `decisions: []`; maintainer allowlist contains `satoshikawato`. No `BD-###` is cited or invented. |
| `OPTION_INTEGRITY_PRODUCT_CONTRACT.md`, revision 5 | PD-OI-001/002/004 limits; 018 complete evidence/search scope; 021 inference ON/OFF; 022 completed raw retry. These constrain both choices but do not select tuple versus graph storage. |
| `REFERENCE/comparison-programs-thresholds-and-results.md` | Similarity groups are search-derived relationships, not phylogenetic orthogroups; group IDs are local; scientific result meanings remain. It does not authorize loss of paths or path metadata. |
| `REFERENCE/session-and-request-compatibility.md`, `SESSION_COMPATIBILITY.md`, repository `CLAUDE.md` | Preserve supported historical input and draft/committed replay distinction; require published namespace + positive fixture for compatibility. Neither mandates one future path format forever. |
| `PRODUCT_IMPACT_RATCHET.md`, architecture ratchet, decision template | Lane B applies to a new public/default/persisted outcome. Material change needs durable base authority; compatibility growth needs a separate architecture exception review. No current exact-head PR-local decision is eligible for this unmapped material change. |
| S01 / MASTER recommendation | Evidence and engineering direction only. Neither is an accepted Product outcome. |

**Preflight sequence:** initially `EVIDENCE_REQUIRED` for the proposed
representation change; S02 investigation/prototype is authorized evidence-only
work. Evidence resolves DAG/rank/transport feasibility. Final path-selection
classification is **`PRODUCT_DECISION_REQUIRED`** because default exhaustive
return/storage and default graph return/storage remain materially distinct,
product-valid outcomes. It is not `IMPLEMENT_EXISTING_AUTHORITY`: no authority
selects the complete future contract. It is not `NOT_ALLOWED`: neither complete
choice requires a hidden cap, path loss, security weakening or branch-only
compatibility. Full runtime/browser/architecture acceptance is still an S06
prerequisite, not waived by choosing an option.

The compact algorithm is a new test-only feasibility implementation under
`tests/prototypes/`. Production owner/path/compatibility sets are unchanged in
S02, so its ordinary ratchet evidence is non-increasing; no S02 architecture
exception is being requested. Candidate documents here authorize no runtime.
