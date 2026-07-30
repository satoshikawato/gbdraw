# Cross-Surface Architecture and Documentation Audit

Date: 2026-07-30

## Executive summary

The rendering core is more convergent than the entry-point code suggests: the
package-root API, typed API, CLI, and Web app all eventually reach the typed request
planner. The convergence happens too late, however. CLI code still performs domain
preprocessing before it creates a request, and Web generation first translates state
into CLI arguments while session saving translates the same state into a separate
JavaScript request. This creates multiple authorities for defaults, validation,
resources, and names.

The audit found seven release-blocking defects. The follow-up implementation resolved
all seven:

1. Circular Web sessions now preserve the comparison FASTA files used by **Run
   LOSAT**. The bundled WSSV session now contains all 20 exact comparison FASTAs:
   12 accession-pinned NCBI Nucleotide records and eight prepared assemblies supplied
   under `tests/test_inputs/MAGs`. Hash, record-ID, length, LOSAT-cache, and BLAST
   coordinate checks prevent substitution with a merely similar sequence.
2. Circular ring labels edited in the Web UI could revert to an older value after a
   session round trip. The current session writer now uses displayed series order for
   BLAST resources, associated FASTAs, labels, and colors.
3. Gallery promotion computed current Web field names but retained the obsolete
   source `config`, so every regenerated Gallery session failed Web import. Promotion,
   current-session validation, and the generated sessions are fixed in this change.
4. CLI rendering now refuses to replace output unless the current invocation includes
   `--overwrite`; saved sessions do not retain that permission.
5. Horizontal inner labels use one shared GC/skew policy across CLI, Web, and Python.
   Explicit track slots remain authoritative.
6. Generated CLI help blocks and an option-parity test keep the CLI Reference current.
7. The three sparse-depth guides now match the tested zero-reserve behavior.

The selected P1–P3 recommendations were also implemented where they represented a
bounded contract decision: strict interactive metadata, wrong-mode override
rejection, distinct typed track names, one canonical typed depth-track input,
an explicit Web UX profile, accurate comparison/similarity terminology, and
evidence-gated compatibility cleanup. A-01, A-02, A-03, A-06, A-07, W-03, W-04,
W-05, D-03, D-05, and D-06 remain architectural or documentation work outside this
follow-up.

The documentation has little vocabulary-level AI writing smell. Its larger problem is
structural: internal session schemas are repeated in tutorials and FAQ answers,
implementation transport terms leak into user instructions, and several captions or
help messages repeat nearby text without adding information.

## Scope and priority

The audit covered:

- Circular and Linear rendering;
- package-root Python API and `gbdraw.api`;
- CLI parsing, preprocessing, rendering, help, and overwrite policy;
- Web generation, session serialization, capability checks, UI defaults, and Gallery
  tutorials;
- README, CLI Reference, Python API guide, tutorials, FAQ, release notes, and
  developer guidance;
- compatibility and apparently obsolete production paths.

The target was the current working tree, not a clean `HEAD` checkout. It already
contained extensive uncommitted production, test, Gallery, and documentation changes;
the findings describe that state. This work also removed branch-only persisted-format
compatibility, migrated bundled sessions to the current writer, updated the reusable
Python API maintenance Skill and compatibility policy, and fixed the two dangling
references that caused CI to fail.

Priority meanings:

| Priority | Meaning |
| --- | --- |
| P0 | Correct before the next release or public deployment; current behavior can lose user state, overwrite files, or materially misdescribe output. |
| P1 | High-impact convergence or contract defect; schedule in the next implementation cycle. |
| P2 | Maintainability, terminology, documentation load, or product-parity work. |
| P3 | Low-priority internal cleanup or a compatibility-sensitive removal or deprecation decision. |

“Confirmed” means the behavior follows from current code, tests, or a reproduced
comparison. “Decision” means the difference may be intentional but is not expressed
as a coherent cross-surface policy.

## Prioritized finding register

| ID | Priority | Status | Finding |
| --- | --- | --- | --- |
| W-01 | P0 | Resolved | A Circular **Run LOSAT** session omitted its comparison FASTA resources and could not rebuild popup sequence state. |
| W-02 | P0 | Resolved | An edited Circular ring label could revert after save/load. |
| W-07 | P0 | Resolved | Gallery promotion retained obsolete Web field names in otherwise current sessions, making all bundled sessions fail Web import. |
| C-01 | P0 | Resolved | CLI output overwrites existing files without an opt-in; Python defaults are safe. |
| C-02 | P0 | Resolved | Circular inner-label policy changes GC/skew only in CLI/Web, not in Python. |
| D-01 | P0 | Resolved | CLI Reference option inventories and copied help are stale. |
| D-02 | P0 | Resolved | Sparse Linear depth geometry is described incorrectly in three guides. |
| A-01 | P1 | Confirmed | Web maintains three render/session representations and imports private CLI parsers. |
| A-02 | P1 | Confirmed | CLI performs domain preprocessing already owned by the typed planner. |
| A-03 | P1 | Confirmed | Root API can lose computed orthogroup metadata from interactive Linear SVG. |
| A-04 | P1 | Resolved | Root and typed interactive-metadata paths have different failure semantics. |
| A-05 | P1 | Resolved | Typed requests accept wrong-mode configuration paths that persist as no-ops. |
| W-03 | P1 | Confirmed | Web capability detection searches Python source text and duplicates fallback maps. |
| API-01 | P1 | Resolved with aliases | Root and typed APIs export incompatible classes under the same track-option names. |
| A-08 | P1 | Resolved | Readers still accepted session/request schemas that existed only on the development branch. |
| D-03 | P1 | Confirmed | Two Gallery basic tutorials use screenshots that do not show the described task. |
| D-04 | P1 | Resolved | Two local documentation links were broken, and README contradicted itself about CairoSVG. |
| API-02 | P2 | Resolved with canonical writer and compatibility inputs | Typed depth input has singular, plural, and matrix forms for the same concept. |
| P-01 | P2 | Resolved by explicit profile | Multi-record grouping and visual defaults differ by surface. |
| P-02 | P2 | Resolved with aliases | “Conservation” and “orthogroup” public names overstate the implemented inference. |
| A-06 | P2 | Confirmed | Current request rendering also owns a large legacy-session migration pipeline. |
| A-07 | P2 | Confirmed | Table/config normalization and metadata construction repeat work within one render. |
| W-04 | P2 | Confirmed | Web injects virtual BLAST inputs by monkey-patching an internal assembler loader. |
| W-05 | P2 | Confirmed | User-facing Web errors repeat deployment instructions and CLI transport details. |
| D-05 | P2 | Confirmed | Session/cache internals are duplicated across tutorials, FAQ, API docs, and release notes. |
| D-06 | P2 | Confirmed | `gbdraw/web/CLAUDE.md` describes an obsolete UI and runtime architecture. |
| C-03 | P3 | Resolved | `save_figure` remains as a deprecated direct-import compatibility path until 0.16. |
| W-06 | P3 | Resolved | The unimported JS style module and transient dialog aliases were removed. |
| C-04 | P3 | Resolved | Release-backed Circular geometry options are documented as supported shortcuts. |

## P0 findings

### W-01: Circular Run LOSAT FASTA resources were not saved

The Web rendering path reads `files.c_conservation_fastas` when the Circular comparison
source is LOSAT. Before this change, the current session request builder materialized
uploaded BLAST files and optional companion sequence sources but did not add the LOSAT
FASTA array. Session load also cleared the transient match-sequence registry, so
Pairwise match and Conservation ring popups could lose sequence lookup even when their
source inputs were present.

The current writer now stores those FASTAs as resources referenced by
`webFiles.conservationLosatFastaSources`, and the current projector restores
`files.c_conservation_fastas`
([session-request.js](../gbdraw/web/js/services/session-request.js), lines 975–994 and
1779–1788). Session import rebuilds the registry from restored GenBank/FASTA inputs,
including Linear selectors, regions, reverse complements, and the displayed
Conservation ring order
([match-sequences.js](../gbdraw/web/js/app/match-sequences.js), lines 339–415;
[config.js](../gbdraw/web/js/services/config.js), lines 3123–3134).
Restoration now uses the same identity-aware Conservation ordering helper as rendering
and session writing. It does not trust a stale `sourceIndex` when a BLAST file identity
still determines the intended companion FASTA.

Result and limitation:

- newly saved sessions preserve available popup sequence sources and rebuild the
  Pairwise/comparison-ring lookup state after import;
- the bundled WSSV session now embeds all 20 exact comparison FASTAs in displayed
  ring order. Twelve are raw NCBI EFetch records for `NC_003225.3`, `AF440570.1`,
  `NC_075105.1`, `AF369029.2`, `AP027278.1`, `AP027279.1`, `AP027284.1`,
  `AP027286.1`, `AP027288.1`, `KT995471.1`, `KY827813.1`, and `MF768985.1`.
  Eight are the prepared files supplied under `tests/test_inputs/MAGs`;
- `tools/restore_wssv_gallery_fastas.py` records source URLs and pinned raw SHA-256
  values, then checks each FASTA against its record ID, sequence length, LOSAT cache
  identity, BLAST query ID, and maximum BLAST query coordinate before embedding it;
- save/load and the Web request projector retain the 20 FASTAs, 20 BLAST tables,
  displayed labels, colors, and source indexes, so popup sequence actions can be
  rebuilt without the network or the untracked MAG directory;
- an arbitrary historical session that stored neither its sequence nor a verifiable
  source still cannot be reconstructed from BLAST coordinates alone. The WSSV repair
  was possible because all 20 exact inputs were independently identified and
  validated, not inferred from the hits;
- the FASTA references remain Web restore metadata. Python/CLI replay models
  precomputed `conservationBlastFiles`; it does not rerun LOSAT from these Web-only
  sources. Converging that contract remains part of A-01.

Remaining recommendation: model a Circular comparison input as one typed union
covering uploaded BLAST plus optional companion FASTA, or analysis FASTA plus program
settings and derived BLAST. Then Python, CLI, Web generation, and saved-session replay
can consume the same request rather than treating LOSAT FASTAs as Web-only metadata.

### W-02: edited Circular ring labels reverted after a session round trip

The active UI edits `circularConservation.series[index].label`, and rendering reads
those series labels
([index.html](../gbdraw/web/index.html), lines 1875–1882;
[run-analysis.js](../gbdraw/web/js/app/run-analysis.js), lines 3096–3105).
Session saving previously read the older comma/newline string
`circularConservation.labels`, even though colors and rendering already used
`series`. Synchronization was one-way from the legacy string to the series.

The writer now treats `series` as the current-state authority and uses the legacy
string only when a source has no series entry. The same ordered entry determines the
BLAST resource, optional comparison FASTA, label, and color
([conservation-series.js](../gbdraw/web/js/app/conservation-series.js), lines 140–193;
[session-request.js](../gbdraw/web/js/services/session-request.js), lines 926–994).
The regression test deliberately leaves the legacy label string stale, reverses two
displayed entries, and checks the complete save/project round trip
([session-request.test.mjs](../tests/web/session-request.test.mjs), lines 838–958).

### W-07: current Gallery sessions contained obsolete Web state fields

The Gallery promoter migrated `depth_tick_interval`, per-track `tick_interval`, and
`collinearMaxGeneGap` while constructing a schema-5 request, but then spread the
unmigrated source session back into the result. All 11 regenerated Gallery sessions
therefore carried at least one obsolete field under `config`. The Web loader rejects
obsolete names in a current session, so Gallery import failed before restoring its
saved SVG.

The promoter now writes the migrated `config` into its result
([gallery-session-migration.js](../gbdraw/web/js/services/gallery-session-migration.js),
lines 348–407). Python validation now rejects the same three obsolete fields when
reading or writing the current session version, while main-backed v33 sessions remain
valid migration inputs
([session_io.py](../gbdraw/session_io.py), lines 640–746). Gallery and current test
sessions were rewritten with the current names. The browser regression reports the
failing session and dialog text directly, then loads all 11 sessions before checking
the WSSV round trip
([test_web_packaging.py](../tests/test_web_packaging.py), lines 4903–5028).

### C-01: CLI output silently overwrites existing files

At audit time, Circular and Linear CLI adapters constructed
`RenderOutputRequest(overwrite=True)`
([circular.py](../gbdraw/circular.py), lines 122–139;
[linear.py](../gbdraw/linear.py), lines 1945–1952). The typed request default is
`False` ([requests.py](../gbdraw/api/requests.py), lines 176–184), and
`Diagram.save()` also refuses replacement unless `overwrite=True`
([interface.py](../gbdraw/interface.py), lines 474–488).

This is a safety inconsistency, not merely a visual default: the same output name can
destroy a file in CLI use and raise a validation error in Python use.

Resolution: Circular and Linear now use the shared safe default and expose
`--overwrite`. Request and session writers always persist `overwrite: false`, and
decoders treat any stored true value as untrusted compatibility input rather than
permission. Session replay and session-sidecar replacement therefore accept overwrite
permission only from the current invocation. Direct session-document writes also
refuse an existing path unless the caller opts in. Known sidecar collisions are
rejected before rendering, so a failed CLI invocation does not leave a new diagram
behind. The browser's ephemeral worker output passes `--overwrite` explicitly because
it writes disposable in-browser files.

### C-02: the same Circular label configuration changes track visibility by surface

At audit time, Circular CLI code suppressed both GC content and GC skew when
horizontal inner labels were enabled
([circular.py](../gbdraw/circular.py), lines 867–876). Web generation reached this
same policy through the CLI adapter. The typed Circular mode profile kept both tracks
enabled by default
([mode_profiles.py](../gbdraw/mode_profiles.py), lines 107–117), and the root/typed
API path does not apply the CLI-only suppression.

The result is a different SVG track set for equivalent options. Existing tests lock
the CLI behavior but do not compare surfaces
([test_circular_feature_width.py](../tests/test_circular_feature_width.py), lines
497–568).

Resolution: the shared Circular builder suppresses implicit GC content and GC skew
when labels use scope `both` with horizontal placement. CLI-only preprocessing was
removed. Explicit track slots remain authoritative, allowing an advanced caller to
retain either track deliberately. Single-record and multi-record regressions cover
the implicit and explicit cases.

### D-01: CLI Reference does not mirror current help

[CLI_Reference.md](./CLI_Reference.md), line 5, states that it mirrors current command
help. A live option-set comparison found:

- 120 unique Circular long options in current help; 18 are not represented in the
  Circular reference section;
- 133 unique Linear long options in current help; 28 are not represented in the
  Linear reference section;
- the copied Circular help block also omits `--circular_track_axis_index`, although
  that option is mentioned later in prose;
- root help currently shows `--protein_blastp_mode` and `--losatp_threads`, but the
  copied main-command block does not.

The missing lists are in the appendix. This drift is predictable because executable
help and a hand-copied reference are separate authorities.

Resolution: `tools/update_cli_reference_help.py` now replaces marked Circular and
Linear help blocks from the live parsers. `tests/test_cli_reference.py` compares the
documented and executable long-option sets, including the new `--overwrite` switch.
Descriptive recipes remain hand-written.

### D-02: sparse Linear depth documentation states the opposite of tested behavior

The implementation and tests say that a missing record/series cell has no paint and a
zero-height reserve band, allowing the following numeric track to compact
([test_linear_vertical_layout.py](../tests/test_linear_vertical_layout.py), lines
436–491; [test_depth_track.py](../tests/test_depth_track.py), lines 1186–1208).
[FAQ.md](./FAQ.md), line 97, and the release notes state this correctly.

The following text says the missing cell retains a reserve band or does not compact:

- [CLI_Reference.md](./CLI_Reference.md), lines 1128–1132;
- [Tutorial 6](./TUTORIALS/6_Depth_Quantitative_Tracks.md), line 110;
- [Tutorial 7](./TUTORIALS/7_Linear_Layout.md), line 142.

Resolution: the CLI Reference and Tutorials 6 and 7 now use the tested rule: the
logical series index is preserved, while a missing cell draws nothing and reserves no
vertical geometry. Later numeric tracks therefore compact upward without being
renumbered.

## P1 architecture and API findings

### A-01: Web has three representations of one render request

The intended architecture routes every surface through one typed planner and limits
adapters to translation
([CLAUDE.md](../CLAUDE.md), lines 118–127). Current Web generation instead:

1. assembles CLI arguments in the 5,304-line `run-analysis.js` module;
2. sends `mode` and `args` to a worker;
3. imports private `_get_args` functions and reparses those arguments in Python
   ([python-helpers.js](../gbdraw/web/js/app/python-helpers.js), lines 16–17 and
   192–203);
4. separately builds a camelCase JavaScript request for session saving
   ([session-request.js](../gbdraw/web/js/services/session-request.js), lines
   855–998).

Thus Web state is translated to a JS session request and independently to CLI argv,
which is then translated again to a Python typed request. W-01 and W-02 are concrete
drift defects produced by this structure.

Recommendation: send one versioned typed request plus resources to the worker and call
the public request renderer. Serialize that same request for sessions. Derive a
display-only CLI command from the request if the UI needs one.

### A-02: CLI preprocessing duplicates planner ownership

Both CLI adapters eventually call `render_request`, but only after loading records and
tables, resolving selectors and regions, applying reverse-complement policy, and
performing comparison cardinality work. Examples are
[circular.py](../gbdraw/circular.py), lines 812–853 and 1154–1157, and
[linear.py](../gbdraw/linear.py), lines 1552–1644 and 1862–1869.

The typed planner already owns record loading, candidate resolution, selectors,
reverse complement, and regions
([request_render.py](../gbdraw/api/request_render.py), lines 413–557). Linear CLI also
imports the private `_select_linear_comparison_records` helper and applies it before
the planner applies the same selection policy.

Recommendation: CLI should create `GenBankInputSource` or `GffFastaInputSource` and
`RecordInput` values without materializing them. Make record cardinality explicit in
the request rather than silently selecting a record only when comparisons happen to be
present.

### A-03 and A-04: root interactive metadata diverges from typed rendering

The Linear assembler can compute orthogroups during rendering and attach them to the
drawing
([diagram.py](../gbdraw/api/diagram.py), lines 2293–2305). Typed rendering reads that
computed attribute when it builds interactive metadata
([request_render.py](../gbdraw/api/request_render.py), lines 1622–1637). The root API
metadata helper only examines the original input options
([interface.py](../gbdraw/interface.py), lines 766–815), so computed metadata can be
absent from `Diagram.to_svg(interactive=True)`.

The two paths also handled metadata failures differently: typed rendering degraded to
an empty context, while the root helper let the exception fail `draw_*` after the
drawing had been built.

A-04 is resolved. Root and typed APIs now use the same strict interactive-metadata
boundary and raise `ExportError` when requested metadata or comparison FASTA content
cannot be prepared. The root API resolves that context lazily only when exporting
interactive SVG, so static drawing remains unaffected.

A-03 remains open: every surface should consume the prepared diagram and interactive
context produced by one owner. Add a root-API test for similarity groups computed
in-process, not only groups supplied by the caller.

### A-05: wrong-mode configuration paths are accepted and ignored

Mode-specific typed options reject a small set of wrong nested option objects, but
generic config overrides accept any recognized global path. Current requests therefore
accept examples such as:

- Circular: `canvas.linear.track_layout` or
  `objects.axis.linear.stroke_color`;
- Linear: `canvas.circular.track_type`.

The fields persist in the request but the selected renderer ignores them. This gives a
false signal that a valid, accepted option affected output.

Resolution: typed option construction now validates mode ownership centrally and
rejects wrong-mode paths such as `canvas.linear.*` on Circular requests and
`canvas.circular.*` on Linear requests. Shared paths remain accepted, and regressions
cover accepted and rejected overrides.

### W-03: capability detection inspects implementation source text

Web startup calls `inspect.getsource(_get_args)` and searches for option-name
substrings, with large JavaScript fallback maps
([run-analysis.js](../gbdraw/web/js/app/run-analysis.js), lines 1355–1644). One shared
option is detected only because its name happens to appear in unrelated help text.
A help rewrite can therefore change a capability result without changing parser
behavior.

Recommendation: export one versioned, machine-readable capability/request schema from
the installed package. Validate it once during worker initialization; typed request
validation should be the final authority.

### API-01: identical exported class names have incompatible constructors

At audit time, the package root exported `CircularTrackOptions` and
`LinearTrackOptions` with
`slots` and `axis_index`
([interface.py](../gbdraw/interface.py), lines 183–227). `gbdraw.api` exported
different classes with the same names and mode-prefixed fields
([options.py](../gbdraw/api/options.py), lines 199–249). Tests already alias the typed
classes to keep them distinguishable.

Resolution: typed integrations now use `CircularRequestTrackOptions` and
`LinearRequestTrackOptions`. Their former short names remain identity aliases for
compatibility. Package-root `CircularTrackOptions` and `LinearTrackOptions` retain
their beginner-facing constructors, and the Python API guide states the boundary
explicitly.

### A-08: branch-only persisted-format compatibility was treated as public

First-parent history on `origin/main` shows public session versions 27–33 and canonical
request schemas 1–2. The current writer is session version 39 with request schema 5.
Session versions 34–38 and request schemas 3–4 existed only between those states on
the development branch, but Python and Web readers, migration helpers, fixtures, and
tests still treated them as supported inputs.

That code had no users to protect and made the current contract harder to reason about.
This change:

- limits readers to main-backed versions 27–33 plus the current version 39;
- limits request codecs to main-backed schemas 1–2 plus current schema 5;
- rejects branch-only versions and schemas at both Python and Web boundaries;
- removes the Gallery migration option that retained schema 3;
- regenerates bundled Gallery and test sessions with the current writer instead of
  preserving an unpublished intermediate representation.

The compatibility policy is now explicit in [CLAUDE.md](../CLAUDE.md) and the
[Python API maintenance Skill](../docs/skills/maintain-python-api/SKILL.md):
compatibility requires first-parent `main` or release evidence and a representative
fixture. Namespace-local schema numbers are evaluated separately, so current protein
cache schemas 3 and 4 are not confused with retired request schemas 3 and 4.

### D-03 and D-04: basic documentation contains direct factual defects

Gallery screenshots:

- [HmmtDNA_basic_circular.json](../gbdraw/web/gallery/tutorials/HmmtDNA_basic_circular.json),
  lines 45–53, asks for four settings, but the image shows only the preset and shows
  Separate Strands unchecked.
- [lambda_basic_linear.json](../gbdraw/web/gallery/tutorials/lambda_basic_linear.json),
  lines 32–40, asks for `NC_001416.gb`, but the screenshot displays
  `BGC0000708.gbk` and *Streptomyces* metadata.

The screenshot checker verifies presence and metadata, not semantic agreement, so
these images can pass automated checking.

Other defects:

- at the start of the audit, [PYTHON_API.md](./PYTHON_API.md) linked to the deleted
  `PYTHON_SESSION_COMPATIBILITY_MATRIX.md`;
- at the start of the audit,
  [RELEASE_NOTES_0.14.0b0.md](./RELEASE_NOTES_0.14.0b0.md) linked to the deleted
  `ADR_PYTHON_SESSION_API.md`;
- [README.md](../README.md), line 29, implies only EPS/PS need CairoSVG, while line 72
  and executable format validation require it for PNG/PDF/EPS/PS.

The two dangling references were removed rather than restoring 355 lines of deleted
internal design documents. The release notes now state the materialization lifetime
directly. README now states once that PNG, PDF, EPS, and PS require CairoSVG, resolving
D-04. D-03 remains: recapture the two task-specific Gallery views.

## P2 design, parity, and documentation findings

### API and runtime structure

- API-02 (resolved with a canonical writer and compatibility inputs): typed requests now use one
  `DepthTrackInput` per logical series under `depth_tracks`. Its `source` is either
  one path/DataFrame shared by displayed records or a per-record sequence of
  path/DataFrame/`None` values. Style and tick fields live on the same object; Linear
  alone accepts `height`. The former singular, plural, and `depth_track_*` matrices
  remain public constructor and decoder inputs and cannot be mixed with the
  canonical form. Python and Web writers always emit `depthTracks`; both decoders
  round-trip the logical-track representation.

- A-06 (legacy replay in current rendering):
  [request_render.py](../gbdraw/api/request_render.py), roughly lines 689–1453,
  combines current planning with legacy protein candidate promotion, legacy ID
  rewriting, and derived-evidence recovery. `session_io.py` contains related
  normalization. Move old-session migration into a session compatibility adapter so a
  fresh request does not traverse legacy concepts.

- A-07 (repeated work): GFF candidate features, color/visibility tables, and
  interactive metadata can be normalized or read in the planner, builder, and metadata
  helper during one render. Carry normalized inputs in `PreparedDiagramRequest` and
  add call-count tests, following the existing load-once test for conservation FASTA.

- W-04 (internal monkey patch):
  [python-helpers.js](../gbdraw/web/js/app/python-helpers.js), lines 126–190, globally
  replaces `gbdraw.diagrams.linear.assemble.load_comparisons` and restores it later.
  The public request already accepts in-memory comparisons. Pass a resource/DataFrame
  through that boundary instead.

### Cross-surface product decisions

| Concept | Current behavior | Assessment |
| --- | --- | --- |
| Circular grouping default | Root API auto-grids multiple records; CLI defaults to separate batch output; Web profile v1 uses `single` for one record, `batch` for several, and exposes `grid` as an opt-in. | The Web writer reads the profile values directly, and behavior tests cover `single`, `batch`, and `grid`. |
| Separate strands | Web profile v1 defaults to true; CLI and config/API default to false. | Intentional UX-profile difference, covered by a profile test. |
| Legend position | Web profile v1 uses Circular left and Linear bottom; CLI/typed default is right. | Intentional UX-profile difference, documented with the profile. |
| Circular legend choices | CLI and Web expose upper-left, upper-right, lower-left, and lower-right. | Resolved; a Web regression covers all four corners. |
| Circular Web input | Web accepts one GenBank/GBFF file or one GFF/FASTA pair per run, while CLI/API accept sequences. | Documented as a deliberate Web limitation. A multi-record GenBank file remains supported. |
| Annotation styling | API/state support many advanced fields; Web editor exposes a small subset. | The Web guide directs advanced styling to annotation TSV rather than implying form parity. |

### Scientific terminology

The public Circular names `ConservationOptions`, `ConservationTrackOptions`, and
`conservation_*` draw filtered BLAST/LOSAT HSP spans; they do not infer evolutionary
conservation. The comparative-genomics tutorial now says this explicitly
([Tutorial 2](./TUTORIALS/2_Comparative_Genomics.md), lines 122–126), but the Python
guide still calls them “Conservation rings”
([PYTHON_API.md](./PYTHON_API.md), lines 157–183).

Likewise, the protein “orthogroup” mode creates gbdraw similarity groups from its
reciprocal-hit/grouping procedure. User-facing text should not imply that this alone
establishes orthology.

Resolution: the package-root API now exports `ComparisonRingOptions` and
`ComparisonRingTrackOptions`, with `CircularOptions.comparison_rings` as the
canonical field. `CircularOptions.conservation` resolves to the same value for
runtime compatibility without storing a second option; `ConservationOptions` and
`ConservationTrackOptions` are identity
aliases. Lower request/session `conservation_*` fields remain transport compatibility
names. Supplying non-null values for both root option names is rejected, and dataclass replacement uses
the canonical field. Web controls, search, popups, standalone interactive SVGs, Gallery tutorials,
and Gallery runtime assets use “similarity group”. Internal `orthogroup` identifiers
remain stable for session and CLI compatibility.

### User-facing text and AI-writing assessment

Vocabulary-level AI smell is low. The usual generic intensifiers, canned conclusions,
and promotional phrases were not prevalent. Do not perform a broad “humanizing”
rewrite.

The specific writing problems are:

- session versions, cache schemas, internal `h_...` handles, and migration details are
  repeated in [Tutorial 4](./TUTORIALS/4_Protein_Comparisons.md), lines 132–145,
  [Tutorial 8](./TUTORIALS/8_Interactive_SVG_Sessions.md), lines 132–158,
  [FAQ.md](./FAQ.md), lines 62–81, Python API integration prose, and release notes;
- the first 40 lines of the CLI Reference lead with retired names, private migration
  transport, and request schemas before current syntax;
- the beginner Python API guide changes audience at lines 318–363 and begins listing
  planners, private migration transport, and session grouping internals;
- Web errors repeat “Current gbdraw wheel” / “Rebuild and redeploy” instructions in
  user-visible option failures; hosted-app users cannot perform that action;
- Web tips expose transport details such as “suppressed by CLI” and “written as the
  CLI feature visibility table” instead of explaining the visible effect;
- `gbdraw cli` says “Oops! It seems like...” and tells users to choose only Circular
  or Linear even though `gui` is also valid
  ([cli.py](../gbdraw/cli.py), lines 164–170);
- the Python API guide says “three output methods” and visually presents four bullets
  by listing two uses of `to_svg()` separately;
- a few Gallery operation bodies repeat their captions verbatim.

Keep release history in release notes and one compatibility page. Tutorials should
retain only the action, expected result, and user-recoverable failure. The FAQ should
answer why a rerun happened without publishing cache implementation schemas.

### Obsolete developer documentation

[gbdraw/web/CLAUDE.md](../gbdraw/web/CLAUDE.md) currently:

- calls Circular mode single-genome only (lines 120–122);
- describes CDN assets and a permissive CSP while the app uses local vendor assets and
  a self-only policy;
- shows main-thread `pyodide.runPython` instead of the generation worker;
- directs contributors to obsolete line numbers such as `setup() ~1339` and
  `args ~6840`.

Rewrite it around module ownership, the worker boundary, local assets, current
multi-record behavior, and the single request/session contract. Avoid volatile line
numbers.

## P3 compatibility decisions

First-parent and release evidence determined the cleanup boundary:

- `gbdraw.render.export.save_figure` exists in released versions through 0.13.0.
  It remains callable, emits `DeprecationWarning` in 0.14, and is scheduled for
  removal in 0.16. Routine tests and internal references now use
  `save_figure_to` or `render_to_bytes`.
- `gbdraw/web/js/app/annotations/style-actions.js` had no production, template,
  test-fixture, or plugin import in this repository and was deleted.
- `colorScopeDialog.individualLabel*` was transient UI state, not persisted
  session data. Current templates and fixtures use `annotationLabel*`, so the
  aliases were removed.
- `--feature_width`, `--depth_width`, `--gc_content_width`,
  `--gc_content_radius`, `--gc_skew_width`, and `--gc_skew_radius` exist in
  releases 0.11.0 through 0.13.0 and remain active Web inputs. They are supported
  geometry shortcuts that expand to track-slot `w` or `r` values. They cannot be
  combined with explicit `--circular_track_slot` or `--circular_track_table`
  geometry.
- The Web Python wrapper now rejects modes other than `circular` and `linear`
  before deleting output files or dispatching a renderer.

The remaining item is internal validation and ownership cleanup:

- request models, public builders, and assemblers resolve or validate some values more
  than once. Separate unresolved input from a resolved plan before removing validation
  needed by direct advanced APIs.

## Recommended implementation order

### Phase 0: protect user state and correct public facts

1. Keep regression coverage for LOSAT FASTA persistence, popup-registry restoration,
   and `series`-ordered save/load behavior.
2. Keep current-field validation at both Python and Web session boundaries, plus the
   browser check that loads every bundled Gallery session.
3. Add safe CLI overwrite policy.
4. Move or remove the CLI-only label/GC/skew policy.
5. Fix sparse-depth prose, broken links, CairoSVG text, and the two basic Gallery
   tutorials.
6. Generate CLI help blocks or add option-set parity tests.

Follow-up status: WSSV now has all exact inputs; overwrite safety, shared Circular
track policy, sparse-depth prose, local links, CairoSVG text, generated help, and help
parity are complete. The two D-03 basic-tutorial screenshots remain.

### Phase 1: converge runtime boundaries

1. Make Web generation and session save consume one typed request.
2. Make CLI adapters construct unresolved typed inputs without domain preprocessing.
3. Centralize prepared interactive context and computed orthogroup metadata.
4. Reject wrong-mode overrides and replace source-text capability probing.

Follow-up status: A-04 strict interactive metadata and A-05 wrong-mode override
rejection are complete. A-01, A-02, A-03, and W-03 remain.

### Phase 2: simplify the public contract

1. Resolve the duplicate track-option names.
2. Consolidate depth inputs.
3. Define one explicit grouping/default profile across surfaces.
4. Introduce accurate comparison/similarity terminology with aliases.

Follow-up status: all four bounded contract decisions are implemented. Compatibility
decoders and aliases retain the previously released or persisted spellings.

### Phase 3: isolate compatibility and reduce documentation load

1. Move legacy session promotion out of current request rendering.
2. Apply the new main/release evidence gate before retaining any compatibility path.
3. Remove the deprecated export path in 0.16; keep the release-backed Circular
   geometry shortcuts documented and tested.
4. Create one compatibility reference and shorten task tutorials and FAQ answers.
5. Update Web developer guidance and remove confirmed dead JS.

## Required regression coverage

Keep or add these tests before the corresponding changes:

1. **Covered:** Circular LOSAT source FASTA survives Web save/load and rebuilds the popup sequence
   registry. Resource projection and display-ordered registry tests were added here;
   retain a browser-level regeneration check.
2. **Covered:** Editing only `series[index].label` survives a Web session round trip.
3. **Covered:** Circular and Linear live-help long-option sets match the CLI Reference inventory.
4. **Covered:** Equivalent Circular label options produce the same GC/skew group policy in CLI,
   root API, and typed API.
5. **Covered:** CLI refuses existing output without explicit overwrite permission.
6. **Remaining (A-03):** Root interactive Linear SVG includes similarity groups computed during rendering.
7. **Covered:** Wrong-mode config paths are rejected.
8. **Covered:** Cross-surface default-profile snapshots cover strandedness, legend, grouping, and
   other intentionally divergent defaults.
9. **Remaining (D-03):** Gallery screenshot validation records which controls and data identity must be
   visible, not only that an image exists.

## Verification performed

- Compared live `circular --help` and `linear --help` long-option sets with the
  corresponding CLI Reference sections.
- Checked local Markdown targets. The two missing files found at the start of the
  audit were the only broken relative document links and are fixed in this change.
- Verified the WSSV repair with `tests/test_wssv_gallery_fastas.py`: all 20 displayed
  rings have exact companion FASTAs, comprising 12 accession-pinned raw NCBI EFetch
  records and eight prepared MAG assemblies. Record IDs, lengths, raw hashes, LOSAT
  cache hashes, BLAST query IDs, and maximum query coordinates agree.
- Loaded, projected, and resaved the WSSV session through the Web session code. Its 20
  BLAST tables, 20 FASTAs, labels, colors, and source indexes remain in displayed
  order.
- Ran all `tests/web/*.test.mjs` files. Session projection, one-record explicit batch
  grouping, default-profile behavior, popup sequence lookup, and branch-only schema
  rejection passed.
- Ran the focused Gallery/session integration checks: 22 Gallery semantic tests and
  30 selected Web packaging tests passed. All 11 Gallery sessions persist
  `overwrite: false`.
- Verified overwrite refusal before rendering for ordinary Circular and Linear CLI
  runs, canonical session replay in both modes, and Circular's implicit output name.
  Existing sidecars were preserved and no SVG was created.
- Recaptured and visually checked 15 Gallery WebP assets. Strict capture checks passed
  for BGC (`16/16/16`), Hepatoplasmataceae collinear (`10/10/10`),
  Hepatoplasmataceae similarity groups (`10/10/10`), and Majanivirus similarity
  groups (`11/11/11`).
- Ran the complete non-slow test suite after the final preflight and test changes:
  2,226 passed, 19 skipped, and 10 deselected.
- Applied the repository `avoid-ai-writing` Skill in detect mode to modified
  user-facing prose and UI strings, then reviewed the matches manually. No publishing
  placeholders, chatbot citation tokens, promotional vocabulary clusters, or generic
  conclusion patterns remained.
- `ruff check gbdraw/` and `git diff --check` passed.

This audit also changed the current Web session writer/loader and popup sequence
restoration, removed branch-only session/request compatibility, migrated bundled
sessions, fixed the CI-breaking documentation targets, and regenerated Gallery-owned
outputs. Tracked reference outputs were not regenerated as part of this audit.

## Appendix: options missing from the CLI Reference at audit time

The generated help blocks now include every option below. The list is retained as the
baseline that motivated the parity test.

Circular (18):

`--circular_label_spacing`, `--depth_color`, `--depth_large_tick_interval`,
`--depth_log_scale`, `--depth_max`, `--depth_min`,
`--depth_small_tick_interval`, `--depth_step`, `--depth_tick_font_size`,
`--depth_window`, `--hide_depth_axis`, `--hide_depth_ticks`,
`--no_depth_log_scale`, `--session`, `--share_depth_axis`,
`--show_depth_axis`, `--show_depth_ticks`, `--tick_label_font_size`.

Linear (28):

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
