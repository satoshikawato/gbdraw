# S09-F01: invalid manifest accepted when no protein raw entry is saved

Recorded before production edits, 2026-09-17 JST.

Severity: medium (persisted metadata integrity). In the inherited final source,
`exportSession` succeeds with a null/malformed live protein manifest if the raw
cache is empty or contains only nucleotide entries. `serializeLosatCache` checks
manifest availability only while filtering protein-current entries; the final
`sessionData` property silently substitutes `emptyProteinIdentityManifest()`.
The real exported `exportSession` Node regression fails with “Missing expected
rejection” (`writer-red.log`). Existing current-protein-raw rejection passes.
No ordinary UI action producing this corruption has been established; the test
injects invalid internal state to exercise the requested writer failure boundary.

Expected: report the existing valid-manifest error and publish no download;
preserve the supplied state. A valid empty manifest remains accepted, including
settings-only and nucleotide workflows. Do not turn bad metadata into successful
empty metadata, weaken the reader, or discard valid inactive raw evidence.

Classification: IMPLEMENT_EXISTING_AUTHORITY. Latest base 5e0cb0fa supplies
OIPC-C01 (reject invalid explicit values), C04 (Session/artifact agreement), C06
(failed reconstruction does not imply deletion), plus strict current Session
identity contracts in docs/SESSION_COMPATIBILITY.md. The S09 request explicitly
requires invalid manifests not to be silently accepted. The newly inherited
writer clarification is explanatory evidence, not authority for itself.
One outcome is determined; no new Product decision or retirement is selected.
EVIDENCE_REQUIRED/PRODUCT_DECISION_REQUIRED/NOT_ALLOWED do not apply to this fix.

Architecture: ordinary non-increasing local correction inside config.js's
existing exportSession owner. Replace the final manifest fallback with an
explicit failure at the existing validation point. Preserve the earlier raw
filter/index try/finally and all reader/derived/Generate boundaries. No new
owner, path, schema, compatibility path, index loan, cache, or protocol.
Registered render-request and Result-admission edges are unchanged. Human S06
exception and S08 writer Review remain pending.

Verification: red/green real-export Node cases for null, wrong schema and missing
analysis references, with empty/nucleotide-only raw caches, no compression or
download on error and state preservation; existing writer/filter and RW gates.
Use the existing browser integration server/observer/Generate/import owners for
source and freshly built installed assets at desktop/mobile, blocked external
requests, failed Save/recovery, valid Save/fresh Load/saved-raw Generate, repeat
and disposal. Reuse unchanged native science and full cancellation lifecycle.
No performance timing, long new search, or historical audit rerun.
