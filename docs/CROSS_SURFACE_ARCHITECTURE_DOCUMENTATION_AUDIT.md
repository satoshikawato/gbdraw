# Cross-Surface Architecture and Documentation Audit

Date: 2026-07-30

## How to read this document

This document records both the state found during the audit and the implementation
that followed it. Text introduced by "At audit time" is historical. Text introduced
by "Current resolution" describes the current working tree.

All findings in this audit are resolved and have focused regression coverage. There
is no remaining implementation phase in this document. Full-suite totals belong in
the final change handoff or CI result so that this audit does not preserve a stale
number.

## Executive summary

At audit time, the package-root API, typed API, CLI, and Web app reached the typed
request planner only after each surface had performed some of its own translation,
loading, or metadata work. Web generation and Web session saving used separate
representations. CLI adapters read domain inputs before constructing a request.
Computed Linear analysis metadata moved through private `Drawing` attributes. The
browser inferred capabilities from Python source text and installed temporary
monkey patches to inject comparison data.

The current implementation gives each boundary one owner:

- Web generation, session saving, and Gallery promotion use the same canonical
  schema-5 request and resource document.
- The browser worker passes that document to a Python resource bridge, which
  materializes bounded inputs, decodes the typed request, calls `render_request()`,
  and reports cleanup failure for its request-scoped workspace.
- CLI adapters construct unresolved typed inputs. The planner owns source loading,
  record expansion, transformations, pairing, naming, and prepared feature inputs.
- Linear assembly returns `LinearDiagramBuildResult`, whose
  `LinearDiagramMetadata` carries computed comparisons, similarity groups, and
  collinearity data to static output and interactive metadata.
- A versioned capability manifest is checked once when the diagram worker starts.
- Released session versions 27–33 are handled by an isolated compatibility adapter.
  Current rendering receives only current request artifacts, and promoted sidecars
  are written from the request that was actually rendered.

The seven release-blocking correctness defects from the original audit are also
fixed: missing Circular LOSAT FASTAs, stale ring labels, invalid promoted Gallery
state, implicit CLI overwrite, cross-surface GC/skew policy drift, stale CLI
reference help, and incorrect sparse-depth documentation.

User-facing Web diagnostics now describe the problem and a user action. They no
longer expose wheel deployment instructions or claim that an option is transported
through the CLI when rendering uses the typed request bridge.

## Scope and status

The audit covered:

- Circular and Linear rendering;
- the package-root Python API and `gbdraw.api`;
- CLI parsing, unresolved request construction, rendering, help, and overwrite
  policy;
- Web generation, session serialization, resources, worker capability negotiation,
  UI defaults, and Gallery tutorials;
- README, CLI Reference, Python guides, tutorials, FAQ, release notes, and developer
  guidance;
- released compatibility inputs and apparently obsolete production paths.

The target was the current working tree, including pre-existing uncommitted changes.
Compatibility was retained only when first-parent `main` or a released artifact
provided evidence for it.

Priority meanings:

| Priority | Meaning |
| --- | --- |
| P0 | A release or deployment blocker involving state loss, unsafe replacement, or a materially false public statement. |
| P1 | A high-impact ownership, convergence, or public-contract defect. |
| P2 | A maintainability, terminology, documentation, or product-parity defect. |
| P3 | Internal cleanup or a compatibility-sensitive deprecation decision. |

Status meanings:

| Status | Meaning |
| --- | --- |
| Resolved / covered | The implementation is complete and focused regression coverage fixes the contract. |
| Resolved / covered (compatibility retained) | The current contract is fixed while a released spelling or API remains as an explicit compatibility alias. |

## Prioritized finding register

| ID | Priority | Status | Finding and current contract |
| --- | --- | --- | --- |
| W-01 | P0 | Resolved / covered | Circular LOSAT sessions retain comparison FASTAs and restore popup sequence sources. |
| W-02 | P0 | Resolved / covered | Circular ring labels, colors, BLAST resources, and FASTAs share displayed-series order. |
| W-07 | P0 | Resolved / covered | Gallery promotion writes migrated current Web fields and validates every bundled session. |
| C-01 | P0 | Resolved / covered | CLI replacement requires `--overwrite`; saved sessions never persist that permission. |
| C-02 | P0 | Resolved / covered | The shared Circular builder owns the horizontal inner-label GC/skew policy. |
| D-01 | P0 | Resolved / covered | Generated help blocks and option-parity tests keep the CLI Reference current. |
| D-02 | P0 | Resolved / covered | Sparse Linear depth documentation matches zero-reserve behavior. |
| A-01 | P1 | Resolved / covered | Web generation and session saving consume one canonical request and resource representation. |
| A-02 | P1 | Resolved / covered | The planner exclusively loads, expands, transforms, pairs, and names unresolved CLI inputs. |
| A-03 | P1 | Resolved / covered | Computed Linear metadata uses `LinearDiagramBuildResult` and `LinearDiagramMetadata`, not private `Drawing` attributes. |
| A-04 | P1 | Resolved / covered | Root and typed interactive metadata use one strict error boundary. |
| A-05 | P1 | Resolved / covered | Typed requests reject configuration paths owned by the other diagram mode. |
| W-03 | P1 | Resolved / covered | A versioned capability manifest replaces Python source inspection and fallback maps. |
| API-01 | P1 | Resolved / covered (compatibility retained) | Typed track-option names are distinct; released short names remain identity aliases. |
| A-08 | P1 | Resolved / covered | Readers accept release-backed versions and current versions, not branch-only formats. |
| D-03 | P1 | Resolved / covered | Basic Gallery tutorials use example-owned, state-checked screenshots. |
| D-04 | P1 | Resolved / covered | Broken local links and CairoSVG requirements are corrected. |
| API-02 | P2 | Resolved / covered (compatibility retained) | Typed depth input has one canonical writer; former constructor forms remain decoder inputs. |
| P-01 | P2 | Resolved / covered | Cross-surface grouping and visual differences are declared in a versioned Web UX profile. |
| P-02 | P2 | Resolved / covered (compatibility retained) | Public comparison and similarity-group names state what is computed; persisted aliases remain. |
| A-06 | P2 | Resolved / covered | Released session promotion is isolated and fails closed on unresolved protein references. |
| A-07 | P2 | Resolved / covered | Prepared feature inputs and interactive inputs are loaded or compiled once per render plan. |
| W-04 | P2 | Resolved / covered | The canonical resource bridge replaces the assembler monkey patch and owns workspace cleanup. |
| W-05 | P2 | Resolved / covered | Web diagnostics and help text use neutral user-facing language without CLI transport details. |
| D-05 | P2 | Resolved / covered | One compatibility reference owns current persisted-format details. |
| D-06 | P2 | Resolved / covered | Web developer guidance describes the current worker and request architecture. |
| C-03 | P3 | Resolved / covered (compatibility retained) | `save_figure` has a release-backed deprecation path through 0.16. |
| W-06 | P3 | Resolved / covered | Unused JS style code and transient dialog aliases are removed. |
| C-04 | P3 | Resolved / covered (compatibility retained) | Released Circular geometry shortcuts remain documented and tested. |

## P0 correctness findings

### W-01: Circular LOSAT FASTA resources

At audit time, a Circular session saved uploaded BLAST data but omitted the LOSAT
comparison FASTAs used to create it. Session import also cleared transient sequence
lookup state. Pairwise-match and comparison-ring popups could therefore lose sequence
actions after a round trip.

Current resolution: the session writer stores the comparison FASTAs as resources and
references them through `webFiles.conservationLosatFastaSources`. Rendering, saving,
projection, and popup-registry restoration use the same displayed comparison order.
The bundled WSSV session embeds all 20 verified comparison inputs, comprising 12
accession-pinned NCBI records and eight prepared assemblies. Record ID, length, raw
hash, LOSAT cache identity, BLAST query ID, and maximum query coordinate checks guard
the repair.

The FASTAs remain Web restore metadata for the uploaded or analysis source. Canonical
Python and CLI replay use the saved comparison tables and do not silently rerun
LOSAT.

### W-02: Circular series labels

At audit time, rendering read labels from `circularConservation.series`, while the
session writer could read an older string field. Editing only the visible series label
could be undone after save and load.

Current resolution: the displayed series entry is the authority for its BLAST
resource, optional FASTA, label, color, and source index. The legacy label string is
used only when no series entry exists. A round-trip regression deliberately leaves
the old string stale and changes display order. Missing optional FASTAs remain
explicit null entries, so history restore, popup source reconstruction, and a later
save cannot shift a FASTA onto the wrong displayed series.

### W-07: Gallery promotion

At audit time, Gallery promotion calculated current replacements for
`depth_tick_interval`, per-track `tick_interval`, and
`collinearMaxGeneGap`, then retained the unmigrated source `config`. Current Web
validation rejected the promoted sessions.

Current resolution: promotion writes the migrated configuration. Python and Web
boundaries reject obsolete fields in a current session. Main-backed legacy sessions
remain valid migration inputs. Browser coverage loads every bundled Gallery session
before exercising the WSSV round trip.

### C-01: explicit overwrite permission

At audit time, Circular and Linear CLI adapters set overwrite permission
unconditionally, while typed and package-root APIs defaulted to refusal.

Current resolution: CLI output uses the shared safe default and exposes
`--overwrite`. Writers persist `overwrite: false`, and decoders do not treat a saved
true value as permission for a later invocation. Output and sidecar collisions are
preflighted for every requested format before rendering. Directories, special files,
dangling symlinks, non-directory parents, and occupied targets without current
overwrite permission fail closed. Missing output directories may still be created by
the normal writer. Package-root, typed, and sidecar writers do not follow a final
symlink that appears after preflight; the sidecar writer uses an exclusively created
same-directory temporary file and a no-replace commit when overwrite is denied.
Batch preflight validates every target before replacing anything, so one invalid
target cannot leave earlier batch items behind. Preflight is set-wide, but
multi-format writes are sequential rather than transactional and completed formats
are not rolled back after a later conversion failure. Record IDs used for implicit
Circular output names must be one filename component; path-like IDs, ASCII control
characters, and Windows-reserved device, stream, or wildcard names require an
explicit output prefix. The browser passes overwrite only for its disposable
workspace output.

### C-02: horizontal inner labels and GC/skew tracks

At audit time, CLI preprocessing hid implicit GC content and skew for horizontal
inner labels. Root and typed API calls did not apply that rule.

Current resolution: the shared Circular builder applies the policy for implicit
tracks. Explicit track slots remain authoritative, so an advanced caller can retain
either track. CLI, Web, package-root API, and typed API regressions cover implicit and
explicit cases.

### D-01 and D-02: executable help and depth prose

At audit time, the CLI Reference omitted live options and copied stale help. Three
depth guides also said that a missing record/series cell retained vertical reserve,
opposite to the implemented geometry.

Current resolution: `tools/update_cli_reference_help.py` refreshes marked help blocks,
and an option-set test compares both mode parsers with the document. The depth guides
now state that a missing cell preserves its logical series index but draws nothing
and reserves no vertical geometry.

## Runtime convergence findings

### A-01: one Web render representation

At audit time, `run-analysis.js` assembled diagram argv for a worker and separately
assembled a camelCase request for session saving. It also imported private CLI parser
helpers and maintained a second set of transport defaults.

Current resolution:
[`buildCanonicalRenderRequest()`](../gbdraw/web/js/services/session-request.js) is
the Web projection used by fresh rendering, session saving, replay downloads, and
Gallery promotion. The
[diagram worker](../gbdraw/web/js/workers/diagram-generation-worker.js) receives
`renderRequest` and `resources`, and the
[Python bridge](../gbdraw/web_support/request_render.py) decodes the request at the
public typed boundary. `run-analysis.js` retains UI validation and analysis staging,
but it no longer constructs diagram argv or imports private parser state. Run Info
exposes one `--session` replay command from the same document.

The canonical Web round trip preserves precomputed protein comparison tables,
similarity-group metadata, source mode, explicit comparison pairs, generated protein
settings, alignment settings, and the selected feature-alignment target. It does not
rerun protein analysis merely to reconstruct a saved render request. After a fresh
render succeeds, the exact resolved comparison artifacts used for that render become
the live session artifacts, so a later ordinary save does not rebuild them from stale
input state. Stored Web-only LOSAT scheduling settings are merged with the projected
request while canonical protein-analysis settings remain authoritative.

### A-02: planner ownership of unresolved inputs

At audit time, CLI adapters read GenBank/GFF records, tables, selectors, regions, and
comparison inputs before constructing typed requests. The planner then repeated part
of that work.

Current resolution: CLI adapters construct `GenBankInputSource`,
`GffFastaInputSource`, `RecordInput`, typed depth inputs, and deferred path values.
`RecordCardinality.EXACTLY_ONE`, `FIRST`, and `ALL` make expansion policy explicit.
The planner loads each unique source once, preserves provenance and source order, then
applies selection, reverse complement, regions, collection presentation, batch
naming, and Linear pair resolution.

Runtime requests may remain unresolved until planning. Canonical schema 5 remains a
materialized persistence format: direct encoding rejects deferred cardinality,
deferred tables, collection transforms, or unresolved output names.
`resolve_request()` provides the materialized projection for writers.

### A-03 and A-04: typed Linear build metadata

At audit time, Linear assembly attached computed similarity groups and related data to
private attributes on `svgwrite.Drawing`. Typed rendering inspected those attributes,
while the package-root API inspected only caller-supplied options. Computed groups
could be absent from an interactive SVG. The two paths also treated metadata failures
differently.

Current resolution: the typed [Linear build path](../gbdraw/api/diagram.py) produces
a `LinearDiagramBuildResult` containing the `Drawing` and
`LinearDiagramMetadata`. The direct compatibility builder retains its `Drawing`
return but is not the preferred public workflow. The metadata includes computed
protein comparison tables, normalized Linear comparisons, similarity groups, and
collinearity results. `PreparedDiagramRequest` and `RequestRenderResult` in the
[request renderer](../gbdraw/api/request_render.py) carry that typed value to
interactive context generation and session artifacts. Production code no longer
attaches or reads computed analysis metadata through private attributes on a drawing.

The package-root API and typed API both call
`build_prepared_interactive_context()`. Interactive metadata is built lazily for
interactive export, and preparation failures raise `ExportError`. Static drawing does
not pay the interactive export cost. Regressions cover groups computed during the
builder call, including a package-root `Diagram.to_svg(interactive=True)` export.

### A-05: mode-owned configuration

At audit time, generic override paths accepted recognized fields from the other
diagram mode and persisted them even though the renderer ignored them.

Current resolution: typed option construction validates mode ownership centrally.
Circular requests reject Linear-only paths, Linear requests reject Circular-only
paths, and shared paths remain valid.

### W-03: explicit runtime capabilities

At audit time, Web startup called `inspect.getsource()` on private CLI parsers and
searched for option-name strings. Large JavaScript fallback maps duplicated expected
support. A help-text edit could therefore alter a capability result.

Current resolution:
[`get_web_runtime_capabilities()`](../gbdraw/web_support/capabilities.py) returns a
versioned, JSON-compatible manifest owned by the Python runtime. The manifest declares
its schema, worker render protocol, supported request schemas, unknown-field policy,
resource encodings, rendering option schema, feature and track renderers, and analysis
artifact schemas. The diagram worker sends it during initialization.
[`validateWebRuntimeCapabilities()`](../gbdraw/web/js/services/runtime-capabilities.js)
compares the complete shape and values once before accepting render work. A mismatch
returns a stable user message and a diagnostic for developers. Source inspection and
option fallback maps are absent.

### A-07: prepare once, consume many times

At audit time, the planner, builder, feature configurator, and interactive metadata
helper could each read or normalize GFF candidates, color tables, visibility rules,
and feature types. Comparison FASTAs could also be parsed more than once for batch
interactive output.

Current resolution:
`PreparedDiagramInputs` owns resolved GFF candidate types and the reusable feature
bundle for a request plan.
[`ResolvedFeatureInputs`](../gbdraw/api/prepared.py) owns loaded default and custom
colors, compiled visibility rules, and their normalized lookup maps.
`PreparedDiagramRequest` carries those values through assembly and interactive
metadata generation. Circular batch items share the same prepared feature inputs.
Comparison FASTA records are memoized on the prepared request and shared across
batch items.

Call-count regressions require each reader, candidate resolver, and compiler to run
once. Downstream tests fail if the diagram or interactive layer tries to read or
compile the source again.

### W-04: canonical resources and workspace lifecycle

At audit time, browser Python helpers temporarily replaced the internal Linear
comparison loader so virtual BLAST inputs could be injected. This global monkey patch
was restored after the call but still bypassed the public request boundary.

Current resolution: the canonical request contains typed comparison resources,
including serialized in-memory DataFrames and their metadata. The Web Python bridge
validates resource identifiers and safe unique filenames, materializes each resource
under a new request-scoped workspace, decodes the request with that resource map, and
calls `render_request()`. No assembler loader is replaced.

The bridge refuses a pre-existing workspace. It removes a newly created workspace
after success and after render or validation failure. Cleanup failure after success
raises a validation error. Cleanup failure during another exception is attached to
the original exception instead of hiding it. The worker likewise preserves a primary
Python or JavaScript error, then attaches proxy-destruction and filesystem-invariant
diagnostics in that order. Focused tests cover these paths and verify that in-memory
comparison metadata reaches the SVG.

### W-05: user-facing Web text

At audit time, some Web errors told users to update a wheel or rebuild a deployment.
Several help tips described CLI suppression, CLI-only options, or generated transport
tables instead of the visible effect.

Current resolution: startup and incompatibility messages ask the user to reload and,
if the problem persists, contact the site administrator. Option help describes what
the control changes in the diagram. Packaging contracts reject the retired wheel,
redeployment, CLI transport, and CLI suppression phrases in user-facing sources.
Developer diagnostics retain the technical mismatch details outside the displayed
message.

## Compatibility and public contract findings

### A-06 and A-08: isolated, evidence-backed session migration

At audit time, current request rendering also carried legacy protein candidate state.
Readers accepted session versions 34–38 and request schemas 3–4 even though those
formats had not reached first-parent `main` or a release.

Current resolution: the current
[request renderer](../gbdraw/api/request_render.py) accepts only
`CurrentRequestArtifacts`. The
[session compatibility adapter](../gbdraw/api/session_compat.py) owns released
persisted-artifact parsing, protein candidate promotion, retired ID rewriting,
derived evidence recovery, and released Linear slot migration. Session versions
27–33 and request schemas 1–2 remain supported alongside current session version 39
and request schema 5. Branch-only versions and schemas are rejected.

The session reader and current renderer share the session-neutral validator in
[`protein_artifacts.py`](../gbdraw/analysis/protein_artifacts.py). It requires the
current derived envelope, rejects retired and feature-analysis IDs, and resolves
every runtime handle through the current protein identity manifest. A strict
zero-hit artifact remains valid without fabricated runtime handles.

Both released replay forms are covered:

- version 27–30 CLI-side sessions migrate old Web configuration fields before
  promotion, render, and sidecar serialization;
- version 31–33 canonical sessions apply the same Web migration, save from
  `rendered.request`, and preserve the migrated adjunct configuration. The released
  version 33 fixture is the positive protein-promotion control.

The resulting `--session_output` is current, omits replay-only output arguments, can
be loaded by the current reader, and can render again. Sidecar targets are preflighted
before diagram output.

Protein references are rewritten only from a verified artifact identity map.
Validation scans DataFrames, comparison objects, similarity-group metadata,
collinearity data, and `align_orthogroup_feature`. It also detects an analysis ID
inside a compound feature string such as
`record@instance|alias~f_<sha256>`. Any retired or feature-analysis reference left
after promotion is rejected. This is a fail-closed postcondition, not a best-effort
cleanup.

### API-01, API-02, P-01, and P-02

At audit time, the package root and `gbdraw.api` exported different track-option
classes under the same names. Depth input had several equivalent current forms.
Surface defaults were implicit, and "conservation" or "orthogroup" text could
overstate the computed evidence.

Current resolution:

- typed integrations use `CircularRequestTrackOptions` and
  `LinearRequestTrackOptions`; former short names remain identity aliases;
- `depth_tracks` writes one `DepthTrackInput` per logical series, while released
  constructor and decoder forms remain accepted compatibility inputs;
- Web UX profile v1 declares grouping, strandedness, legend placement, and other
  intentional surface defaults;
- package-root comparison-ring names and Web "similarity group" text describe the
  displayed BLAST/LOSAT result. Persisted and CLI transport aliases remain stable.

### C-03, C-04, and W-06

Release evidence keeps `save_figure` callable with `DeprecationWarning` until its
scheduled 0.16 removal. Released Circular geometry shortcuts remain supported and
expand to track-slot values; they cannot be combined with explicit slot geometry.
The unimported JS style module and transient dialog aliases had no persistence or
release contract and were removed.

## Documentation findings

### D-03 and D-04: factual examples and links

At audit time, two beginner Gallery tutorials used screenshots from another example.
Two local document links pointed to deleted design files, and README gave conflicting
CairoSVG requirements.

Current resolution: the HmmtDNA and lambda tutorials own task-specific WebP captures.
Each data-dependent operation declares its Gallery session, expected app state, and
the controls or text that must appear in the crop. Shared mode-selector media is
marked as generic. The broken references were removed, and README states that PNG,
PDF, EPS, and PS require CairoSVG.

### D-05: one compatibility reference

At audit time, session versions, cache schemas, migration details, and internal handle
formats were repeated across task tutorials, FAQ, Python API prose, CLI Reference,
and release notes.

Current resolution: [Session and request compatibility](./SESSION_COMPATIBILITY.md)
owns the current support and migration contract. Task documents retain the action,
expected result, and recoverable failures, then link to that reference. Typed
integration material lives in [TYPED_API.md](./TYPED_API.md). Release notes retain
historical changes. Documentation contracts prevent current task guides from
reintroducing internal compatibility terms.

### D-06: Web developer guidance

At audit time, [the Web maintenance guide](../gbdraw/web/CLAUDE.md) described
single-genome Circular behavior, CDN assets, permissive CSP, main-thread Pyodide
rendering, and source locations by approximate line number.

Current resolution: the guide documents module ownership, the canonical schema-5
request, the typed diagram worker, separate LOSAT workers, current Circular
single/grid/batch behavior, Linear topology, local assets, self-only CSP, session
ownership, and browser verification. It avoids source line references.

### Writing quality

The audit found few vocabulary-level machine-writing patterns. The documentation did
not need a broad tone rewrite. The implementation removed copied captions, generic
filler, obsolete transport explanations, and repeated compatibility internals where
they obscured the task. Exact interface labels, option names, file names, and
scientific limitations remain unchanged where accuracy requires them.

## Completed implementation phases

| Phase | Result |
| --- | --- |
| 0: protect user state and correct public facts | Complete. W-01, W-02, W-07, C-01, C-02, D-01, D-02, D-03, and D-04 are resolved and covered. |
| 1: converge runtime boundaries | Complete. A-01 through A-05 and W-03 are resolved and covered. |
| 2: simplify the public contract | Complete. API-01, API-02, P-01, and P-02 have canonical names or writers with explicit compatibility inputs. |
| 3: isolate compatibility and reduce documentation load | Complete. A-06 through A-08, W-04 through W-06, D-05, D-06, C-03, and C-04 are resolved and covered. |

No audit finding is deferred to a later phase. The scheduled 0.16 removal of
`save_figure` is an existing deprecation lifecycle, not unfinished work for this
change.

## Regression coverage

The following contracts are covered:

1. Circular LOSAT FASTAs, including missing-file placeholders, survive Web
   save/load and history restore without shifting popup sequence sources.
2. A label changed only in `series[index].label` survives a session round trip with
   its resource, color, and display position.
3. Circular and Linear live-help option sets match the CLI Reference.
4. Equivalent Circular label options use the same implicit GC/skew policy on every
   surface.
5. Package-root, typed, CLI, and session replay paths preflight every requested
   output and sidecar before rendering. Occupied or non-replaceable targets, invalid
   parents, late final-name collisions, and unsafe implicit record-ID names fail
   closed; multi-format and batch target-set validation precedes construction.
6. Web generation and session saving use one canonical request builder. Tests reject
   diagram argv and private-parser fallback.
7. Canonical Web projection preserves precomputed protein comparisons,
   similarity-group and collinearity results, generated settings, comparison pairs,
   alignment settings, and Web-only runtime settings without rerunning analysis.
   A successful fresh render adopts the same resolved artifacts for later saves.
8. Interactive Linear SVG includes similarity groups computed during assembly, for
   typed and package-root API calls.
9. Prepared GFF candidates, colors, visibility rules, feature types, and comparison
   FASTAs are loaded or compiled once.
10. Worker capability negotiation validates the versioned manifest once and rejects
    a mismatched runtime before rendering.
11. Canonical resources carry in-memory comparison metadata without a loader monkey
    patch, and temporary workspaces follow the defined cleanup behavior.
12. Released version 27–33 sessions produce current, reloadable sidecars. Unresolved
    exact or compound protein references fail closed.
13. Wrong-mode configuration paths are rejected.
14. Web UX profile tests cover intentional grouping, strandedness, legend, and
    multi-record differences.
15. Gallery capture contracts verify data identity, app state, and visible crop
    content rather than file presence alone.
16. User-facing source contracts reject obsolete deployment and CLI transport
    phrases.

## Verification evidence

Focused implementation checks cover these areas:

- typed request planning, record cardinality, load-once preparation, computed Linear
  metadata, interactive export, and batch output preflight;
- released session compatibility, current sidecar replay, protein artifact identity
  rewriting, exact and compound fail-closed references, and alignment targets;
- canonical Web request/resource decoding, comparison metadata, workspace cleanup,
  capability negotiation, session projection, and static transport guards;
- WSSV companion FASTA identity and Gallery session round trips;
- CLI help/reference parity, documentation ownership, Markdown targets, and Gallery
  capture contracts;
- Python lint, package build, Web asset packaging, browser startup/session replay
  checks, reference output comparison, and diff hygiene.

The documentation contract suite and an explicit stale-language scan are rerun after
this file changes. Final full-suite pass, skip, and deselection totals are reported in
the implementation handoff or CI result, not frozen into this audit.

Tracked diagram reference outputs are not regenerated unless an intentional geometry
change is reviewed through the repository's reference-update workflow.

## Appendix: CLI options missing at audit time

Generated help now includes every option in this historical baseline.

Circular:

`--circular_label_spacing`, `--depth_color`, `--depth_large_tick_interval`,
`--depth_log_scale`, `--depth_max`, `--depth_min`,
`--depth_small_tick_interval`, `--depth_step`, `--depth_tick_font_size`,
`--depth_window`, `--hide_depth_axis`, `--hide_depth_ticks`,
`--no_depth_log_scale`, `--session`, `--share_depth_axis`,
`--show_depth_axis`, `--show_depth_ticks`, `--tick_label_font_size`.

Linear:

`--align_orthogroup_feature`, `--collinear_max_diagonal_drift`,
`--definition_line_style`, `--depth_color`, `--depth_large_tick_interval`,
`--depth_log_scale`, `--depth_max`, `--depth_min`,
`--depth_small_tick_interval`, `--depth_step`, `--depth_tick_font_size`,
`--depth_window`, `--hide_accession`, `--hide_depth_axis`,
`--hide_depth_ticks`, `--hide_length`, `--linear_label_spacing`,
`--no_depth_log_scale`, `--protein_blastp_candidate_limit`,
`--protein_blastp_max_hits`, `--record_subtitle`, `--save_session`,
`--session`, `--session_output`, `--share_depth_axis`, `--show_depth_axis`,
`--show_depth_ticks`, `--show_replicon`.
