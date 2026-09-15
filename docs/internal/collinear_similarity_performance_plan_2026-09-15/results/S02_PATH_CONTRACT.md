# S02 — concrete alternatives and S06 implementation contract

**Engineering design for the selected PATH-B outcome; not Product authority.**
Concern `protein-comparison.path-representation`, scenario revision **1**.
The complete PATH-B receipt was received on 2026-09-15. Authority review and
base integration remain pending; see [decision handoff](S02_PRODUCT_DECISION.md).
Implementation names/schema allocations remain proposals, not new human
commitments. PATH-A below is the unselected comparison alternative.
Base and observations: [inventory](S02_INVENTORY.md).

## 1. Common requirements

Both choices preserve complete retained evidence, membership and annotations,
group order/IDs, path protein sequences, selected edge sequences, final rank,
path IDs, edge first-path IDs, shared tuples, block/anchor membership and
ordering, genomic links, SVG geometry and semantic attributes. No path cap,
sampling, approximate count, edge reversal or score-based path replacement.

All existing biological modes, inference ON/OFF, member/candidate limits, actual
search scope, raw keys, manifest validation, display transforms and raw retry
authority remain. OIPC revision 5 decides those matters; no new choice is needed.
Save/Load/History retain draft versus committed Result, with failure/cancel/stale
completion leaving the last committed Result intact. No new path browsing UI.

“All paths” is the current deduplicated **ordered language**, not all possible
edge walks. The graph must retain non-path evidence edges and original edge
order separately from its canonical transition selection.

## 2. PATH-A / exhaustive-current

- Keep existing Python dataclasses, return types, dictionary/tuple/list shapes,
  defaults, empty mappings/tuples, exception behavior and all-path fields.
- Keep normal writers and schema admission unchanged. Typed resources and rich
  group/derived metadata contain complete arrays; catalog still reduces them.
- Explicit retrieval is existing tuple indexing/iteration. `len(tuple)` gives
  the exact materialized count. No safe count can be computed by parsing an
  already rounded JavaScript Number; larger-than-safe counts require an exact
  source if exposed independently in the future.
- Equivalent implementation improvements may remove duplicate temporary
  copies, but producing/holding/serializing the promised arrays costs at least
  **Ω(L)**. A generator upstream of a JSON array does not remove that cost.
- Preserve the existing cycle-prefix behavior wherever that enumerator is
  invoked; preserve externally supplied tuples exactly during typed replay.
- No new compatibility reader or graph type. Rollback of an equivalent S06
  optimization restores the previous implementation under the same formats.

This is feasible and does not waive a performance regression gate. It explicitly
leaves MASTER §1's elimination of exponential normal output **unachieved**.

## 3. PATH-B / lossless-graph: public access

The following names/signatures define the candidate behavior; no functions in
this section have been installed into `gbdraw.api` by S02.

### Result types and explicit exhaustive continuation

1. Keep `OrthologPath`, `OrthologEdge` and the current **`OrthogroupResult` class**
   as actual frozen dataclasses with their existing fields and constructor/tuple
   contract. They remain accepted precomputed inputs and exhaustive outputs.
2. Add `OrthogroupGraphResult`: the same non-path metadata fields, in their
   existing declaration order, plus
   `path_indexes_by_orthogroup_id: dict[str, OrthologPathCollection]` replacing
   the all-path map. Do not attach a property named
   `ortholog_paths_by_orthogroup_id` which secretly enumerates or returns a lazy
   stand-in for a tuple.
3. Existing inference-producing public functions gain a keyword-only
   `path_representation: Literal['graph', 'exhaustive']`.
   **PATH-B selects `graph` as their new default.** The public selection/result
   wrappers' `.orthogroups` types become
   `OrthogroupResult | OrthogroupGraphResult | None` as appropriate; the
   Collinearity result's nested orthogroup field follows the same rule. Inference
   OFF still returns no inferred graph. Bad mode values raise `ValidationError`.
4. `path_representation='exhaustive'` is an explicit compatibility continuation:
   return the previous dataclass/tuple shapes and values. It is allowed to cost
   Ω(L). Share the inference, filtering and projection owners; this is an output
   materialization request, not another inference pipeline or cache-miss fallback.
5. `materialize_ortholog_paths(result: OrthogroupGraphResult) -> OrthogroupResult`
   constructs the old tuple-valued result. If given an `OrthogroupResult`, it
   returns that already explicit value. Native drawing/encoding never calls this
   function implicitly. Third-party code depending on the old default must add
   the explicit keyword or call this conversion. That migration is a **material
   Product effect**, not merely choosing another implementation class.

The ordinary Python `draw_linear`/typed render, CLI and Web use graph production.
Existing precomputed `OrthogroupResult` input is normalized at the common
analysis/typed-preparation boundary without changing its supplied corpus.
S06 must check public signature overloads, options validation, expected resource
types, and all wrappers together. An old dataclass must never suddenly encode
as an empty-path result just to pass old type checks.

### `OrthologPathCollection` contract

Two immutable, explicitly tagged data forms have one access owner:

- **`dag`**: graph-derived source-to-sink language with current selection/order
  rules. Normal inference produces this form.
- **`explicit`**: an already supplied legacy `tuple[OrthologPath, ...]`, preserving
  order, arbitrary supplied IDs and shared tuples. It is used by the one legacy
  input adapter. It can contain a selected subset, noncanonical IDs or cyclic
  evidence, because current typed readers permit those values. It is not a
  production fallback for failed DAG validation.

| API | Type, order and failure contract |
|---|---|
| `count` | Python arbitrary-precision `int`, exactly the number of emitted entries after protein dedup for `dag`; stored tuple length for `explicit`. Zero is valid. No `__len__` that can overflow `Py_ssize_t`; use `.count`. |
| `path_at(rank: int)` | One-based rank in final order; returns one current `OrthologPath`. Non-int, including bool, raises `TypeError`; rank <1 or >count raises `IndexError`. Explicit form uses stored order and stored ID, not a synthetic renumbering. |
| `path_by_id(path_id: str)` | For DAG, exact `group.path_[1-9][0-9]*`; no leading zeros. Wrong type: `TypeError`; malformed/wrong-group ID: `ValidationError`; valid out-of-range rank: `IndexError`. For explicit form, exact stored ID lookup, absent ID: `KeyError`, duplicate stored ID: `ValidationError` for this new lookup API. Preserve duplicate data in enumeration. |
| `rank_of(protein_ids: tuple[str, ...])` | Exact full source-to-sink sequence; returns one-based `int`. Wrong container/item type: `TypeError`; absent/incomplete/invalid transition: `ValidationError`. Explicit duplicate protein sequences are preserved; ambiguous lookup raises `ValidationError`. |
| `iter_paths()` | Explicit `Iterator[OrthologPath]` in final order, no cap. Each returned value has exactly the old protein/selected-edge/shared tuples. Not an implicit `__iter__` used by generic serializers. Stopping iteration need not materialize the remainder. |
| `containing_count(protein_id: str)` | Exact number of emitted paths containing the node, counting it once per path; unknown node: `KeyError`. DAG counts derive from prefix/suffix DP. Explicit input can compute from its stored corpus. It does not overwrite legacy supplied shared tuples. |
| `first_path_id(edge)` | DAG: earliest rank among paths using the selected edge-ID string; if absent, preserve `edge.path_id` exactly. Explicit legacy input preserves the existing edge field, including a supplied ID absent from its corpus. This is not a new consistency validator for old data. |

Specific retrieval need not allocate other path objects. S06's iterator should
use a lexicographically ordered DFS with a reusable traversal stack; the test
prototype deliberately uses repeated `path_at` for simplicity and is **not**
an optimized full-output streamer. All-output time/output stays Ω(L), even if
count and random access are efficient. No browser paging or download UI is added.

### DAG construction and validation

- Obtain path participation from the same kinds and endpoint membership rules
  as the legacy owner. Keep the complete original evidence list.
- Collapse parallel `(query_id, subject_id)` transitions by minimum **edge-ID
  string**, not score or biological edge-kind rank. The independent choices
  minimize the full edge tuple for each protein sequence.
- Persist each path node's current K order key, including its original protein-ID
  tie-break value. The graph's order cannot depend on later display coordinates,
  annotation edits or dictionary order. Order changes require a new inferred
  result; raw search identity remains unaffected.
- Use Kahn/topological traversal on the actual directed transitions, **not a
  record-index scan**. In canonical inference IDs/order keys are unique.
- A new `dag` payload validates exact tags/fields, canonical decimal counts,
  unique nodes, valid endpoint references, selected transition identities,
  unique/injective order keys, and acyclicity. Recompute summaries and compare
  persisted count/first IDs; reject inconsistent new data with
  `CanonicalRequestDecodingError` at typed decoding (analysis constructor:
  `ValidationError`). New graph input must not invent another protein-identity
  normalization owner or trust a claimed count.
- A cycle in generated inference violates the proved boundary: fail the result
  before success publication and retain the prior Result/raw retry. Do not
  reverse, truncate, return count zero, or run the old exhaustive routine as an
  automatic fallback.
- Old dataclass/typed corpus is different: adapt directly to `explicit`. Do not
  reject its cyclic edges or infer missing paths. S06 must preserve the existing
  accepted input surface for `path_representation='exhaustive'`; noncanonical
  caller-built protein maps must not be silently admitted as canonical DAGs.
  If reinspection finds valid normal producer inputs outside the unique-identity
  proof, resolve that boundary with evidence before enabling the new default.

## 4. Exact count, rank, first-edge ID and shared proof

After transition dedup, let `S(v)` be the number of suffix paths from v to a
sink, `F(v)` the number of source prefixes ending at v, and children be K-sorted.
All arithmetic below is integer arithmetic.

```text
S(sink) = 1                     S(v) = sum(S(child)) otherwise
F(source) = 1                   F(v) = sum(F(parent)) otherwise
P = sum(S(source))              containing(v) = F(v) * S(v)
offset(v, child) = sum(S(earlier sibling))
sourceOffset(source) = sum(S(earlier source))
rank(path) = 1 + sourceOffset(first) + sum(offset(each transition))
```

Unrank subtracts whole sibling subtree counts until it finds the branch
containing the desired rank. It constructs only that path. The K sequence is
injective on canonical protein sequences; the final edge tuple tie-break is
already determined by endpoint transition selection.

First-containing edge rank can be computed **without constructing even one
path per edge**. Let `B(v)` be the minimum zero-based rank offset of a prefix
ending at v:

```text
B(source) = sourceOffset(source)
B(v) = min(B(parent) + offset(parent, v))
firstRank(u→v) = 1 + B(u) + offset(u, v)
```

For a fixed prefix/transition, the least suffix contributes zero additional
offset. Topological relaxation therefore yields the earliest complete path
containing the selected transition. `F(v)*S(v)` counts the Cartesian product
of prefix and suffix choices; their only shared node is v in a DAG. Mark nodes
with product >1 and K-sort a retrieved path's marked nodes to reproduce shared
tuples. Duplicate original edge objects retain equal IDs and identical first
path assignment; nonselected kinds/IDs preserve prior values.

Storage is O(V+E) **graph entries**, plus integer bit storage and order-key/ID
bytes. Sorting costs O(V log V + Σ d(v) log d(v)) comparisons; building edge
indexes is expected linear dictionary work. Counts and offsets use O(V+E)
integer additions/comparisons; shared counts use O(V) big-integer products.
With b-bit counts these are not constant-time operations. Their bit/decimal
conversion cost must be measured, not hidden inside a “linear time” claim.
Path retrieval scans source/child ranges in the prototype; prefix sums plus
binary search are an optional S06 implementation choice, not a new Product
outcome. Path storage itself never contains P objects in ordinary computation.

### Wire integers

Use canonical unsigned decimal **strings** (`0|[1-9][0-9]*`) for counts and ranks
crossing JSON/Worker/Session boundaries. No signs, decimals, exponent syntax,
whitespace, leading zeros, Infinity or floating-point intermediates. Python
computes `int`; JS keeps text for display/serialization and uses `BigInt` only
for arithmetic/comparison. Never JSON-stringify a raw JS BigInt or coerce a
count through `Number`. Existing path IDs already contain decimal rank text.

The complete R=56 DAG has `18014398509481984` paths, a power of two which alone
is a weak rounding test. A disconnected extra two-node component gives
**`18014398509481985`**, which is not exactly representable as a JS Number.
The prototype tests both and 2^53±1. Actual unchanged Python/JS catalog validators
and both popup routes retain the latter count string. Catalog schema 3 can
therefore remain 3: its count field already accepts text. Normal new summaries
omit zero-count display fields as the existing reduced catalog does; collection
`.count` and explicit rich export still distinguish zero from unknown.

## 5. Normal render/save, explicit old data, and schemas

### Proposed normal path

```text
same validated evidence → one inference owner → members + edges + DAG index
    → link/block projection (same first-path strings)
    → summary metadata (count; no path iterator)
    → catalog / interactive SVG / popup
    → compact typed resource / derived artifact / Session
```

`_encode_typed_tree` needs explicit support for graph dataclasses; it must not
reflectively visit a legacy all-path property. `serialize_orthogroups_payload`
needs a summary projection used by interactive context and native derived
generation. Its explicit exhaustive projection remains opt-in and materializes
the old rich shape. `_build_current_derived_entries` must not construct rich
arrays before choosing a compact payload. The helper and native derived writer
should carry the same canonical compact analysis resource and pair/provenance
data using S05's established completed-result owner. Do not add another retained
completed payload, a metadata cache or a second inference owner in S06.

### Candidate typed shape (names are concrete; version allocation is conditional)

Keep descriptor `kind: orthogroupResult`, `encoding: canonicalJson`, and the
Collinearity descriptor/value kinds. The next typed-resource schema admits
`OrthogroupGraphResult` plus its path collection tags. Non-path fields keep
their existing camel-case names and meanings. A `dag` collection contains:

```text
kind: "dag"
ordering: "protein-key-edge-tuple-v1"
orthogroupId: string
nodes: [{proteinId: string, orderKey: [recordIndex, start, end, proteinIdTieBreak]}]
transitions: [{queryProteinId: string, subjectProteinId: string, edgeId: string}]
count: canonical decimal string
```

Nodes and transitions are serialized in deterministic K / endpoint order.
Original full edges, roles, diagnostic evidence, existing path IDs and member
metadata remain in the enclosing result fields. Derived DP tables are rebuilt
on read; do not serialize a second authority for shared counts or ranks. First
path IDs on full edges are verified for new DAG data. `explicit` contains
`kind`, `orthogroupId`, and typed `paths` preserving the old exact corpus;
its count derives from its stored tuple. No edge closure is inferred for it.

| Namespace | PATH-A | PATH-B proposal on the inspected base |
|---|---|---|
| Typed JSON resources | Writer 2/readers 1,2 | Next writer **3**; readers 1/2 adapt explicit input, 3 reads tagged graph/explicit. Preserve schema-1 optional anchor fields. |
| Derived protein payload | Writer 3; legacy 1 separately retained | Next writer **4**, identity includes path representation/order version and upstream raw keys. Published 3 adapts to the current representation once; legacy 1 uses the existing identity-verification lane, never assumed a current cache hit. |
| Session envelope | Writer 42 | Next writer **43**, so older readers reject early. Preserve all currently supported Session input versions and settings-only behavior; normalizing 42 does not require Generate. |
| Canonical request | 7 | Remains 7: descriptor grammar and requested biological settings unchanged; resource version handles the analysis result representation. |
| Catalog / interactive SVG | 3 / 3 | Remain 3 / 3: count strings already survive validation and both popup consumers. New graph data lives in typed artifacts, not duplicated in every SVG/catalog. |
| Raw / manifest / Web bindings | 4 / 2 / 2 | Unchanged; no path array belongs to these identities. |

These numbers are a **conditional allocation**, not frozen implementation
constants: S06 must fetch and inspect its actual base, retain its active writer,
and allocate the next unoccupied version once. A branch-only intermediate must
be rewritten directly to that version before merge, removing its reader/test/
artifact. Do not create 42→43→44 chains if 43 never reaches main/release.

### Compatibility and rollback

Three new persisted source-form obligations are planned for PATH-B: published
typed resource 2, derived 3, and Session 42 become superseded inputs. Typed 1 and
older supported Session/derived lanes already exist and are extended directly
to the new current form, not chained through branch-only formats. The explicit
legacy tuple adapter is shared by Python input and typed reading; count it once
as an implementation owner, not once per caller. Its persisted dispatch branches
still count separately in the architecture compatibility ledger.

S02 adds **zero** production compatibility paths; PATH-A adds zero. PATH-B has
**three new persisted-format handling obligations** in this design, plus the
explicit Python return/input continuation. This is not an invented repository
CB total. S06 must enumerate the actual before/after stable reader/branch IDs,
owner sets and canonical paths for the affected namespaces. A positive CB delta
triggers a complete architecture exception packet and a separate maintainer
review on the final implementation head. Product selection cannot waive it.
Removal condition: retain each published reader while its format remains
supported; remove only under separately accepted retirement authority, or
remove an unmerged candidate version before publication. No expiry date is
invented here.

Before publishing a graph writer, runtime rollback is removal of that candidate
implementation and its branch-owned generated artifacts. After graph-bearing
Sessions are actually published, an old binary rejects them. Keep original old
Sessions/backups; a new reader remains necessary for new files. A reverse export
to old exhaustive data is explicit and costs Ω(L); do not promise that every
large graph can be downgraded within available memory or that a byte-level
version edit is a conversion. Failed explicit materialization reports its error
and does not replace the committed Session/Result. Runtime recovery must retain
new-format readability until separate retirement authority exists.

## 6. Evidence, limitations and S06 start checklist

Executable assets:

- `tests/prototypes/ortholog_path_graph.py`: exact DAG count, rank/unrank,
  first-edge rank, shared counts; no production import or writer.
- `tests/test_ortholog_path_contract.py`: named full oracles, 1,024 five-node DAG
  subsets, 160 seeded eight-node DAGs with duplicate edges and coordinate ties,
  3 named + 60 seeded real inference cases, S01 R=8/12/16 complete result equality,
  cyclic observations, permissive legacy typed round trip, errors, JS big integers
  and actual catalog/popup decimal-string admission.
- `tests/prototypes/path_publication_history.py`: read-only Git/tag/fixture
  evidence regeneration; [machine reports](data/).

The prototype's large R=24/32/56 cases run only the compact index. Normal count,
all first-edge ranks and shared counts create **zero `OrthologPath` objects**.
Separate explicit single-path tests verify selected ranks and the path named by
each edge. This is not a full-selector speedup or a browser heap measurement.
S01's 24 pass / 3 regression / 15 inconclusive same-code timing comparison
remains unsuitable as evidence of a speedup; no threshold is relaxed.

S06 may start dependent runtime only when all of these are true:

1. S05 is complete and its common analysis/cache result owner is identified.
2. A complete explicit Product receipt selects PATH-A or PATH-B for revision 1;
   only that wording has been serialized, reviewed and merged as authority.
   For this unmapped concern the existing OIPC authority route is appropriate;
   do not create a new BD store or PR-local self-authorization block.
3. The latest `origin/dev` contains the accepted durable authority, not just
   these candidate design files. Reclassify changed semantics/owners/schemas and
   refresh positive publication evidence if the base moved.
4. For PATH-B, verify each AND-of-OR requirement in the Decision Pack, including
   the explicit exhaustive API and arbitrary legacy corpus, not only the graph
   count. Prove canonical producer identity/order assumptions on that base.
5. Implement all expansion boundaries above in one convergent path; retire the
   superseded eager normal path. No inference duplicate or silent fallback.
6. Run the same complete small oracles, typed read/write negative/positive cases,
   edge/block/SVG semantic and geometry checks, Session/draft/replay/History and
   raw/derived lifecycle tests. Add real offline browser save/reload/Generate
   tests for the **new** format and actual helper/Worker route, plus Circular
   smoke if shared render/cache changes. S02's scalar Node tests do not replace
   that browser gate.
7. Measure normal render/encode/save at R≥24 without invoking legacy enumeration;
   verify integer text and ordering beyond 2^53, small-input overhead, peak memory,
   cancellation and rollback. Use the S01 runner for comparable timing, with
   unchanged 10% regression / 5% noise thresholds and uncontended repeated runs.
8. Complete any required architecture exception review separately from Product
   approval; publish no writer until required acceptance is complete.

No S03–S08 optimization, new runtime/writer/reader, public artifact, authority
record or schema constant is changed by S02. Pending Product judgment is a
dependency for S06, not a reason to redo the completed S02 evidence.
