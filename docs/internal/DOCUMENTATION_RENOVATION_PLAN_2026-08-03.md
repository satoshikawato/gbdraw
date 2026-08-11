# gbdraw documentation renovation and executable tutorial plan

Date: 2026-08-03

Status: Implemented and verified

Primary audience for the new tutorial: biologists who want a publication-ready genome diagram without Python or command-line knowledge

Secondary audiences for the full renovation: biologists using reproducible CLI workflows and developers integrating the public Python APIs

Tutorial language and format: English, plain Markdown, reproducible PNG screenshots

## 1. Decision

The first end-user learning path will remain one English browser-GUI tutorial:
`bundled HmmtDNA.gbk -> first Circular preview -> focused customization -> regenerated final figure -> verified SVG download`.

That first path is the entry point, not the full scope. The renovated documentation set will
also provide explicit chapters for Linear diagrams, browser LOSATN/TLOSATX/LOSATP workflows,
precomputed comparisons, quantitative and custom tracks, editing and sessions, publication
exports, the Circular and Linear CLI, and the public Python interfaces. Every major capability
listed in Section 5 must have one canonical chapter, a stable fixture, and an executable proof.

The implementation will follow these rules.

1. The first diagram must appear in Step 2. This satisfies the stricter reading of the
   `document-generate` rule that the first visible result occurs in fewer than three steps.
2. The tutorial will follow one fixed Circular workflow. Linear layouts, comparisons,
   custom colors, quantitative tracks, sessions, and alternative export formats will appear
   only as links after `What you built`.
3. A committed Python Playwright flow will perform every documented browser action from a
   clean context. It will not restore a finished session or mutate application state through
   `window.__GBDRAW_APP__`.
4. Every screenshot will be produced by that flow. A screenshot will be included only when
   it helps the reader choose a setting, locate an action, or verify an outcome.
5. The tutorial fixture will be a version-pinned GenBank file shipped with both the hosted
   and local web app. Public prose will not depend on `tests/test_inputs`.
6. SVG export is part of the standard path. Playwright must capture the download event and
   validate the saved file, not merely click the button.
7. Existing numbered “tutorials” that are task-oriented will be rewritten and indexed as
   How-to guides. Their current public paths will be retained during this renovation to avoid
   unnecessary link breakage.
8. The correctness work identified as P0 in the audit remains release-blocking and will be
   completed before the documentation set is declared renovated.
9. All GUI actions in Tutorials and How-to guides will be executed through Playwright. CLI
   commands will be executed from a clean temporary directory, and Python examples will be
   imported and executed as tests. A surface is never considered documented by testing only a
   different surface.
10. Coverage will be broad at the chapter level and narrow within each chapter. Each chapter
    solves one reader goal; exhaustive controls, option types, defaults, and schemas stay in
    Reference.

This plan applies the local
[`love-me-love-my-docs`](../../.agents/skills/love-me-love-my-docs/SKILL.md)
workflow to application inspection, capture, and verification, and the local
[`gstack/document-generate`](../../../gstack/document-generate/SKILL.md) source to
Diátaxis partitioning and document quality. The repository-local
[`document-generate`](../../.agents/skills/document-generate/SKILL.md) copy was verified as
byte-identical during planning. The tutorial structure also follows the
[official Diátaxis tutorial guidance](https://diataxis.fr/tutorials/).

## 2. Source material and inspected evidence

The plan is based on the following sources.

- [Documentation audit](./DOCUMENTATION_AUDIT_2026-08-03.md), including all P0-P2 findings.
- Repository and web-maintenance guidance in [`CLAUDE.md`](../../CLAUDE.md) and
  [`gbdraw/web/CLAUDE.md`](../../gbdraw/web/CLAUDE.md).
- Browser navigation and controls in [`gbdraw/web/index.html`](../../gbdraw/web/index.html).
- Existing Gallery beginner tutorial data in
  [`HmmtDNA_basic_circular.json`](../../gbdraw/web/gallery/tutorials/HmmtDNA_basic_circular.json).
- Existing screenshot automation in
  [`capture_gallery_tutorial_screenshots.py`](../../tools/capture_gallery_tutorial_screenshots.py)
  and offline generation checks in
  [`verify_gui_offline.py`](../../tools/verify_gui_offline.py).
- Documentation and Gallery contract tests under `tests/` and `tests/web/`.

The planning inspection established the following facts.

| Finding | Evidence | Consequence |
|---|---|---|
| The repository has no committed English Markdown tutorial for the browser GUI. | `docs/DOCS.md` and `docs/TUTORIALS/README.md` route GUI learners to the hosted Gallery; `docs/QUICKSTART.md` is CLI-only. | Add one canonical Markdown GUI tutorial and link it directly from README and the docs landing page. |
| The current Circular Gallery tutorial first generates in manual step 4. | `HmmtDNA_basic_circular.json` lines 275-319. | Reorder the learning path so upload and Generate produce a preview in Step 2. |
| Existing data-dependent Gallery captures load completed sessions and then assign Vue state through JavaScript. | `capture_gallery_tutorial_screenshots.py` lines 1134-1148 and the `evaluate` actions in `HmmtDNA_basic_circular.json`. | Gallery captures are useful visual references but cannot prove that documented user actions work. |
| The existing final image asserts a right legend, but the tutorial never tells the user to select Right; the fresh Circular default is Left. | `HmmtDNA_basic_circular.json` lines 304-308 and `gbdraw/web/js/web-ux-profile.js` lines 3-10. | The new flow must set and assert every value that affects the promised result. |
| The app has a file uploader but no bundled example loader. The Gallery’s “Input downloads” route to NCBI. | `index.html` lines 1082-1093 and `gallery.js` lines 225-257. | Ship a stable raw tutorial fixture and give the tutorial a direct local/same-origin download target. Do not depend on a live NCBI response. |
| The current end-user fixture exists only in `tests/test_inputs/HmmtDNA.gbk`. | Repository fixture inventory and package-data configuration. | Promote it to web tutorial data and update test references atomically. |
| Core form labels are not consistently connected to their controls, and no `data-testid` attributes were found in the web app. | `index.html` file-uploader template lines 5963-5980 and selector census. | Add accessible `for`/`id` or `aria-label` relationships first; add test IDs only where accessible roles are ambiguous. |
| Python Playwright, Node Playwright, and Chromium are available in the inspected environment. | `playwright --version`, Python import, and `require.resolve('@playwright/test')` checks. | A real-browser capture flow is feasible; degraded mode is not currently needed. |
| Existing Gallery media validation and the offline GUI generation smoke test pass. | `capture_gallery_tutorial_screenshots.py --example HmmtDNA_basic_circular --check --strict` and `verify_gui_offline.py smoke-test`. | The runtime and existing assets are a sound base, but the new user-action path still needs its own smoke capture. |
| The current public product has three entry subcommands: `circular`, `linear`, and `gui`; the Python package also exposes simple and typed programmatic interfaces. | `gbdraw/cli.py`, `gbdraw/__init__.py`, `docs/PYTHON_API.md`, and `docs/TYPED_API.md`. | Separate first-result Tutorials and execution harnesses are required for GUI, CLI, and Python rather than treating one surface as proof for another. |
| The GUI exposes browser LOSATN, TLOSATX, and LOSATP; LOSATP includes Pairwise, Similarity groups, and Collinear blocks. | `gbdraw/web/index.html`, `gbdraw/web/js/state.js`, LOSAT worker/services, and Gallery tutorials. | Nucleotide, translated-nucleotide, protein-pairwise, group, and Collinear goals need distinct chapter ownership and real browser runs. |
| The CLI's built-in comparison search path is LOSATP/BLASTP through `--protein_blastp_mode`; nucleotide and translated-nucleotide comparisons are consumed as prepared BLAST/LOSAT tables. | `gbdraw/linear.py`, `gbdraw/circular.py`, `docs/WORKFLOW_GUIDE.md`, and current comparison guides. | A tested surface capability table must prevent GUI LOSATN/TLOSATX instructions from being copied into CLI prose. |
| The 11 Gallery tutorials collectively show beginner Circular/Linear, AT skew, region annotations, multi-record canvases, similarity groups, Collinear blocks, and Circular rings, but most are not canonical raw-input executable documentation. | JSON files under `gbdraw/web/gallery/tutorials/` and the screenshot register. | Reuse their biological recipes and visual targets, not their session-preload/state-injection capture pattern. |

## 3. Concept map

| Item | Definition for this renovation |
|---|---|
| Purpose | Lead a non-programmer through the first browser diagram, then provide executable task chapters for every major GUI, CLI, and Python capability. |
| Key concepts | Annotated inputs, Circular versus Linear layout, generate/regenerate feedback, comparison evidence, tracks and presentation, interaction/session state, and publication/reproducible export. |
| Public surface | Hosted/local GUI, `gbdraw circular`, `gbdraw linear`, `gbdraw gui`, simple Python API, typed requests, supported inputs/comparisons/tracks/editors/sessions/outputs. |
| Dependencies | Packaged fixture suite, offline web assets, Pyodide diagram worker, LOSAT Wasm/native runtime, sanitized SVG preview, export services, Python runtime, and vendored fonts. |
| Downstream documents | README, docs landing page, tutorial index, Gallery beginner entry, export guide, workflow guide, relevant how-to guides. |
| Common failure modes | Worker/runtime not ready, mismatched inputs/IDs, unsupported surface claim, empty comparison after filtering, stale session/settings, blocked download, missing optional exporter, malformed or misnamed output. |
| Design decisions | Put the small mitochondrial GUI path first; use additional small fixtures for Linear/LOSAT/GFF3/depth; show early results; use SVG as the standard publication path; give alternatives separate goal-focused chapters. |

## 4. Target information architecture

Diátaxis roles will be explicit in headings, directories, and navigation. Because the current
Markdown site has no redirect mechanism, the target directory structure will be introduced in
stages and old public paths will temporarily act as compatibility landing pages.

| Reader need | Canonical content after renovation | Required change |
|---|---|---|
| Tutorial: learn a complete workflow | New surface-specific chapters under `docs/TUTORIALS/GUI/`, `CLI/`, and `PYTHON/` | Put the first Circular GUI workflow first, then separate Linear, LOSAT, CLI, and Python learning paths. Every Tutorial has an early result and one uninterrupted successful route. |
| How-to: accomplish a known task | New task chapters under `docs/HOW_TO/GUI/`, `CLI/`, and `PYTHON/` | Give each major task one canonical goal-oriented procedure. Migrate the useful parts of the existing numbered guides without retaining their course-like classification. |
| Reference: look up exact behavior | Target chapters under `docs/REFERENCE/`, migrated from `CLI_Reference.md`, `PYTHON_API.md`, `TYPED_API.md`, `SESSION_COMPATIBILITY.md`, `SVG_SEMANTIC_HOOKS.md`, plus Web/Input/Comparison/Output references | Keep types/defaults/constraints here. Move browser control and format tables out of Tutorials and How-to guides. |
| Explanation: understand choices | Target chapters under `docs/EXPLANATION/`, migrated from `WORKFLOW_GUIDE.md` and focused publication/reproducibility material | Explain interface, layout, comparison, runtime, and publication trade-offs. Remove the user-inaccessible term `Web UX profile v1`. |
| Auxiliary | `GALLERY.md`, `FAQ.md`, `ABOUT.md`, release notes | Call Gallery a curated showcase unless every entry has a canonical recipe. Keep FAQ as support and release notes as history. |
| Internal | Audits, implementation/capture plans, maintenance registers and skills | Move or retain them under `docs/internal/` and `docs/skills/`, with audience and owner stated. Do not expose them through end-user navigation. |

The first navigation block in `README.md`, `docs/DOCS.md`, and
`docs/TUTORIALS/README.md` will route readers by intent:

1. **New to gbdraw?** Open the browser GUI tutorial.
2. **Prefer reproducible commands?** Open the CLI first-diagram tutorial.
3. **Already know the basics?** Choose a How-to guide by task.
4. **Need an exact option or default?** Open Reference.
5. **Need to choose an interface or understand a trade-off?** Open Explanation.

The browser tutorial must be reachable from `README.md` in at most two clicks.

The path migration will be staged. Existing public files such as `docs/QUICKSTART.md` and
`docs/TUTORIALS/1_...` through `9_...` remain at their current URLs for one renovation release
as short compatibility landing pages that point to the new canonical chapter. They must not
duplicate runnable instructions. This preserves inbound links while allowing the public index
to expose a coherent Tutorial/How-to/Reference/Explanation structure.

## 5. Chapter Plan Gate

No capture or chapter-writing implementation may begin until this complete census is reviewed
once. Approval means that every major capability has a canonical owner and that the proposed
fixtures are small enough for repeatable execution. It does not mean that every chapter must be
implemented in one pull request.

### 5.1 Major-feature coverage contract

The documentation scenario manifest proposed in Section 8 will use the following capability
groups as its coverage vocabulary. A release may not claim complete documentation renovation
while any row lacks a canonical chapter and an executable scenario.

| Capability group | Major capabilities that require an owner |
|---|---|
| Entry points | Hosted GUI, local GUI through `gbdraw gui`, Circular CLI, Linear CLI, public Python interface, typed request interface. |
| Inputs | GenBank/EMBL/DDBJ, GFF3 + FASTA pairing, multi-record files, record selection, regions, reverse complement, and TSV manifests. |
| Layout | Single- and multi-record Circular canvases; single- and multi-row Linear diagrams; track placement, rulers, axes, titles, definitions, and legends. |
| Comparisons | Uploaded BLAST/LOSAT tables, browser LOSATN, browser TLOSATX, LOSATP Pairwise, LOSATP Similarity groups, LOSATP Collinear blocks, selected/mixed edges, and Circular similarity rings. |
| Quantitative and annotation tracks | GC content, arbitrary dinucleotide skew including AT skew, depth/numeric tracks, region annotations, and custom Circular/Linear track slots. |
| Feature presentation | Palettes, feature-specific colors, qualifier priority, label filters and overrides, visibility, arrow/rectangle/underlay shapes, strokes, and overlap handling. |
| Interactive work | Preview search, feature/match/group popups, feature and legend editing, layout editing, undo/redo/reset, session save/load, and cached comparison reuse. |
| Outputs | Static SVG, interactive SVG, PNG, PDF, EPS, PS, exported comparison/FASTA data, and reproducible command/session handoff where supported. |
| Programmatic use | Simple Python drawing API, in-memory records and bytes, comparisons/tracks, typed requests, and session round trips. |

### 5.2 Tutorial chapters

These are learning-oriented journeys. They are ordered deliberately: the first chapter is the
standard browser workflow, and optional layouts, comparisons, and programming surfaces follow
only after it. Every Tutorial produces its first diagram by Step 2 or, when runtime preparation
cannot be combined safely, no later than Step 3.

| ID | Proposed English chapter | Audience and outcome | Canonical fixture | Executable proof | Planned images |
|---|---|---|---|---|---:|
| T-GUI-01 | **Create and export your first circular genome diagram** | A non-programmer creates, labels, customizes, and exports a publication SVG. | `HmmtDNA.gbk` | Clean-context Playwright flow and parsed SVG download. | 6 |
| T-GUI-02 | **Create and export your first linear genome diagram** | A non-programmer creates a labeled Lambda map with a ruler and exports SVG. | `NC_001416.gb` | Clean-context Playwright flow, real Linear generation, and parsed SVG download. | 4 |
| T-GUI-03 | **Compare two genomes with LOSATN in the browser** | A biologist creates a nucleotide comparison without installing a search tool. | Whole Lambda `NC_001416.1` and whole Enterobacteria phage DE3 `NC_042057.1`. | Playwright runs browser LOSATN with serial/one-thread settings, verifies the six frozen links, and downloads raw LOSAT TSV and SVG. | 5 |
| T-GUI-04 | **Compare annotated proteins with LOSATP in the browser** | A biologist starts with five annotated records and finishes with protein Similarity groups. | All five whole records in `aminoglycoside-bgc-five`. | Playwright runs LOSATP Similarity groups, opens one group, verifies members, and downloads TSV/SVG. | 5 |
| T-CLI-01 | **Create a reproducible circular diagram from the command line** | A CLI user creates the same small Circular diagram in a clean directory. | `HmmtDNA.gbk` | Subprocess contract checks exit status, filename, SVG structure, accession, features, GC/skew, and ticks. | 0 screenshots; 1 generated figure |
| T-CLI-02 | **Create a reproducible linear diagram from the command line** | A CLI user creates a Lambda map with labels and ruler. | `NC_001416.gb` | Subprocess contract checks exit status, filename, record, axes, labels, and ruler. | 0 screenshots; 1 generated figure |
| T-PY-01 | **Draw and save your first genome diagram from Python** | A Python user loads a record, draws one Circular diagram, and saves it. | `HmmtDNA.gbk` | The exact documented code is imported and executed by pytest; its output is the published figure. | 0 screenshots; 1 generated figure |

T-GUI-03 and T-GUI-04 each teach one comparison model. Pairwise links, Collinear blocks,
TLOSATX, Circular rings, and threshold variations are linked after `What you built` rather than
added as branches inside those Tutorials.

### 5.3 GUI How-to chapters

Every documented browser action in this table is owned by a committed Playwright scenario.
The image estimate is a ceiling, not a quota; screenshots are retained only for a setting choice
or visible verification.

| ID | Proposed English chapter | Reader goal and scope | Primary fixture/scenario | Proof | Images |
|---|---|---|---|---|---:|
| H-GUI-01 | **How to use GenBank and GFF3 + FASTA inputs** | Load GenBank/DDBJ or a matched GFF3 + FASTA pair and diagnose record-ID errors. | Lambda GenBank and the whole-record Lambda GFF3/FNA derivative. | Upload both paths; assert record ID, complete length, CDS count, strands, and validation messages. | 3 |
| H-GUI-02 | **How to arrange multiple circular records** | Put several records on one Circular canvas and choose placement/size behavior. | Four complete, naturally circular mitochondrial RefSeq records: *Homo sapiens*, *Danio rerio*, *Drosophila melanogaster*, and *Caenorhabditis elegans*. Large Vibrio GBFF remains a Gallery case study. | Assert the four accessions, complete-record and circular-topology metadata, grid rows/columns, labels on every record, CDS `gene` labels with no CDS `product` labels, and non-overlap semantics. | 3 |
| H-GUI-03 | **How to arrange linear records, regions, and orientation** | Set rows, order, crop ranges, reverse complement, record gaps, definitions, and ruler. | Lambda/multi-record fixture. | Assert record order, displayed coordinate range, orientation, row topology, and ruler. | 4 |
| H-GUI-04 | **How to use uploaded BLAST results** | Draw selected Linear comparison edges from prepared outfmt 6/7 files, including omitted edges. | Whole Lambda and DE3 records plus their frozen LOSATN/TLOSATX TSV. | Upload the real TSV; assert only selected endpoints and expected links are drawn. | 3 |
| H-GUI-05 | **How to use TLOSATX for translated nucleotide comparisons** | Choose translated-nucleotide evidence when nucleotide identity is insufficient. | Whole Lambda `NC_001416.1` and whole DE3 `NC_042057.1`. | Run browser TLOSATX serially with one thread and a 1,000 bp minimum alignment length; validate 397 raw rows, seven displayed links, endpoints, SVG links, and TSV download. | 3 |
| H-GUI-06 | **How to add circular sequence-similarity rings** | Run or upload comparison evidence and place one or more rings around a Circular reference. | Complete, naturally circular mitochondrial genomes: human `NC_012920.1` as the reference and *Danio rerio* `NC_002333.2`, *Drosophila melanogaster* `NC_024511.2`, and *Caenorhabditis elegans* `NC_001328.1` as three rings; WSSV remains advanced Gallery material. | Assert all four whole-record identities/topologies, three-ring count/order and endpoints, human CDS `gene` labels with no CDS `product` labels, selected reference side, HSP popup, and span-FASTA download. | 4 |
| H-GUI-07 | **How to create protein similarity groups with LOSATP** | Build all-record search-derived groups, inspect a group, and optionally align by a member. | Five BGC records. | Run actual LOSATP; assert 155 features, expected non-zero matches/groups, stable group IDs, popup, and raw TSV. | 4 |
| H-GUI-08 | **How to draw collinear protein-match blocks with LOSATP** | Build ordered blocks and choose evidence scope and orientation/identity color mode. | Five complete Hepatoplasmataceae genomes restored with the qualified All-vs-all cache. | Reuse all 25 cached results; assert all-record evidence scope, adjacent display ribbons, orientation, popup, and two-span FASTA export. | 3 |
| H-GUI-09 | **How to add depth, GC content, and skew tracks** | Add one/multiple depth series, axes, GC percentage, GC skew, or AT skew. | Reduced depth fixture plus HmmtDNA. | Assert logical series count, kind/order, axis scale/ticks, legend, and values. | 4 |
| H-GUI-10 | **How to add region annotations and custom track slots** | Import bracket/highlight/band annotations and place annotation/numeric tracks. | Tobacco chloroplast record + region TSV; HmmtDNA AT-skew. | Assert annotation set/lane, renderer, slot order/side, and visible region labels. | 4 |
| H-GUI-11 | **How to style features, labels, titles, and legends** | Apply palettes/rules, qualifier priority, label filtering/overrides, title, and legend choices. | HmmtDNA or BGC plus repository-managed TSV rules. | Assert control values and SVG fill colors, label text, title, and legend position/order. | 4 |
| H-GUI-12 | **How to control feature visibility, shapes, strokes, and overlaps** | Hide/show selected features and choose arrow, rectangle, underlay, stroke, or overlap behavior. | HmmtDNA plus a small visibility rule table. | Generate from raw input; assert visible counts, element shapes, strokes, and overlap policy. | 3 |
| H-GUI-13 | **How to inspect and edit a finished diagram** | Search features; open feature, pairwise, and group popups; bulk-edit color/visibility/stroke; move layout/legend. | Completed HmmtDNA and BGC scenarios. | Playwright performs preview interactions and asserts edits survive regeneration/export. | 5 |
| H-GUI-14 | **How to save, restore, undo, and reproduce your work** | Use undo/redo/reset, download a session, reload it in a fresh context, and obtain a reproducible handoff. | State created by T-GUI-01 and T-GUI-04. | Capture a real session download, validate current schema/resources, reload fresh, and compare SVG semantics. | 3 |
| H-GUI-15 | **How to export publication and interactive figures** | Export SVG, interactive SVG, PNG, and PDF and choose the right artifact. | Completed T-GUI-01 state. | Click every real export action; check filenames, magic bytes/XML, dimensions, embedded behavior, size, and sanitization. | 2 |

### 5.4 CLI How-to chapters

CLI chapters are executed as commands from an empty temporary working directory. Terminal
screenshots are intentionally omitted; copyable commands, captured stdout/stderr, and the
generated figure are better evidence.

| ID | Proposed English chapter | Scope | Primary executable proof |
|---|---|---|---|
| H-CLI-01 | **How to use GenBank, DDBJ, GFF3, and FASTA inputs** | Input selection, GFF3/FASTA ID matching, multi-contig records, and validation. | Run both GenBank and Lambda GFF3/FNA recipes; assert records, CDS, strands, and failure messages. |
| H-CLI-02 | **How to use record, comparison, conservation, annotation, and track tables** | TSV manifests and their path-resolution behavior. | Execute one minimal recipe per manifest; validate schemas, dependencies, and output semantics. The record-table recipe uses the four complete metazoan mitochondrial genomes; the conservation-table recipe uses three complete-mtDNA TLOSATX comparisons against human mtDNA. Both assert every CDS label comes from `gene`, not `product`. |
| H-CLI-03 | **How to arrange multiple circular records and tracks** | Multi-record canvas positions/size modes and Circular track order/slots, using the four complete metazoan mitochondrial genomes. | Assert accessions, complete-record and circular-topology metadata, grid topology, axis index, track order, and per-record CDS `gene` labels with no CDS `product` labels in generated SVG/run metadata. |
| H-CLI-04 | **How to arrange linear records, regions, orientation, labels, and rulers** | Rows, selectors, crop, reverse complement, track layout, definitions, and ruler. | Execute every recipe in order from one empty directory; assert row/order/range/orientation and scale. |
| H-CLI-05 | **How to draw comparisons from precomputed BLAST results** | Linear `--blast` for complete Lambda/DE3 and Circular `--conservation_blast` for complete metazoan mtDNA. | Assert HSP filtering, query/subject mapping, six Linear Lambda/DE3 links, three human-reference mtDNA rings, and CDS `gene` labels without product labels. |
| H-CLI-06 | **How to run Pairwise protein searches with LOSATP** | `--protein_blastp_mode pairwise`, runtime selection, threads, hit limits, and saved evidence. | Run pinned LOSAT with one thread; assert raw hits and Pairwise links. |
| H-CLI-07 | **How to create protein similarity groups with LOSATP** | `orthogroup` compatibility token, group membership, labels, and alignment target. | Assert stable group membership/IDs, labels, and aligned record placement. |
| H-CLI-08 | **How to draw collinear protein-match blocks with LOSATP** | Collinear scope, anchors, gaps, and color modes. | Assert block count, anchor/order metadata, orientation, and colors. |
| H-CLI-09 | **How to add depth, GC content, and skew tracks** | Single/multiple depth series, axis scaling, GC percentage, skew, and custom order. | Assert series/kind/order, min/max/log scale, axes, ticks, and legend. |
| H-CLI-10 | **How to add region annotations and custom track slots** | Annotation tables and explicit Circular/Linear slot geometry/placement. | Assert annotation IDs/lanes, renderer ownership, slot order, side, and axis placement. |
| H-CLI-11 | **How to set colors, labels, visibility, shapes, and strokes** | Palettes, tables, qualifier priority, filters, overrides, rendering shapes, and overlap behavior. | Run one canonical recipe per goal and compare semantic features rather than entire SVG bytes. |
| H-CLI-12 | **How to save and regenerate sessions safely** | Plain/gzip sessions, `--save_session`, session input, overwrite protection, and version errors. | Round-trip current sessions, reject incompatible inputs, and prove overwrite behavior in a temporary directory. |
| H-CLI-13 | **How to export static and interactive outputs** | SVG/PNG/PDF/EPS/PS and standalone interactive SVG, optional export dependencies, names, and overwrite. | Validate every produced file by signature/parser, expected filename, non-zero size, and interactive controls. |

CLI comparison guidance must state the surface boundary explicitly: built-in CLI search covers
LOSATP/BLASTP through `--protein_blastp_mode`; nucleotide LOSATN and translated-nucleotide
TLOSATX evidence is prepared outside the CLI workflow and read through `--blast`,
`--comparisons_table`, or `--conservation_blast`.

### 5.5 Python How-to chapters

| ID | Proposed English chapter | Scope | Proof |
|---|---|---|---|
| H-PY-01 | **How to draw Circular and multi-record diagrams from Python** | `read_genbank`, `draw_circular`, layout options, and save/bytes output; the four-record example uses complete *H. sapiens*, *D. rerio*, *D. melanogaster*, and *C. elegans* mitochondrial genomes. | Execute the exact snippets; assert all four accessions, complete-record and circular-topology metadata, tracks, per-record CDS `gene` labels with no CDS `product` labels, and output bytes. |
| H-PY-02 | **How to draw Linear diagrams and comparisons from Python** | `draw_linear`, rows, prepared comparisons, and comparison options. | Execute prepared-comparison and multi-row examples; inspect topology and links. |
| H-PY-03 | **How to add tracks, annotations, colors, and labels from Python** | Depth, similarity rings, custom slots, feature/label options. | Type validation plus SVG semantic assertions and the published figure recipe. |
| H-PY-04 | **How to use GFF3, in-memory records, and byte output** | `read_gff`, `SeqRecord`, file-like/in-memory inputs, `to_svg`, and `to_bytes`. | Execute each boundary and validate biological records and output signatures. |
| H-PY-05 | **How to build typed requests and round-trip sessions** | Mode-specific typed requests, resource planning, render results, and session codecs. | Request validation, render, session round-trip, and wrong-mode error tests. |

### 5.6 Reference, Explanation, and Auxiliary chapters

| Role | Proposed English chapter | Canonical ownership |
|---|---|---|
| Reference | **Web app reference** | Panels, accepted inputs, mode-specific availability, defaults, limits, validation, output actions, and keyboard/accessibility behavior. |
| Reference | **Command-line reference** | Generated `circular`, `linear`, and `gui` help plus curated semantics that cannot be inferred from argparse output. |
| Reference | **Input formats and TSV schemas** | GenBank/EMBL/DDBJ, GFF3/FASTA identity rules, BLAST outfmt, depth, records/comparisons/conservation/annotation/track tables. |
| Reference | **Comparison programs, thresholds, and result semantics** | LOSATN/TLOSATX/LOSATP capability matrix, Pairwise/groups/Collinear meanings, filtering, reference sides, and caches. |
| Reference | **Palettes, feature rules, labels, shapes, and track renderers** | Exact tokens, precedence, defaults, and constraints; this is where control/option inventories belong. |
| Reference | **Python API reference** | Stable top-level types, signatures, accepted sources, returns, errors, and version support. |
| Reference | **Typed request reference** | Mode-specific request types, ownership, validation, resources, and results. |
| Reference | **Session and request compatibility** | Current version/support matrix and migration behavior; release notes retain only historical deltas. |
| Reference | **Output format and export reference** | Format support by GUI/CLI/API, optional dependencies, naming, DPI, fonts, overwrite, and limitations. |
| Reference | **Interactive SVG and semantic-hook reference** | Embedded controls/data, search/popups, stable IDs, sanitization, and integration hooks. |
| Reference | **Tutorial fixture provenance** | Accessions, licenses/sources, checksums, derivations, record/feature counts, and intended scenarios. |
| Explanation | **Choose Circular or Linear** | Biological/layout trade-offs and when comparison topology changes the choice. |
| Explanation | **Choose the GUI, CLI, or Python interface** | Interactivity, reproducibility, scale, automation, privacy, and installation trade-offs. |
| Explanation | **Choose a genome-comparison method** | Precomputed evidence versus LOSATN, TLOSATX, LOSATP Pairwise, Similarity groups, Collinear blocks, or Circular rings. |
| Explanation | **Understand tracks, axes, and layout ownership** | Simple presets versus custom slots, axes, overlay/underlay behavior, and why some controls become inactive. |
| Explanation | **Prepare publication and reproducible figures** | Vector/raster choice, fonts, accessibility, provenance, sessions, commands, and visual QA. |
| Explanation | **Understand browser privacy, offline execution, and performance** | Local/hosted execution, external requests, Pyodide/Wasm workers, cache behavior, and practical data-size limits. |
| Auxiliary | **Gallery** | A curated showcase. Only entries with raw stable inputs and an executable recipe may claim reproducibility. |
| Auxiliary | **Palette Explorer** | Visual palette lookup linked to exact palette tokens in Reference; it is not a Tutorial chapter. |

Skip from Tutorials: full control inventories, schema history, every threshold variant, debug or
internal routes, large Gallery case studies, maintainer scripts, and architectural contracts.
Those remain Reference, Explanation, Auxiliary, or Internal content.

### 5.7 Chapter manifest and approval record

Create `docs/scenarios/manifest.json` before implementation. Each entry records:

- stable chapter/scenario ID, title, Diátaxis role, surface, and priority;
- capability tags from Section 5.1;
- canonical fixture IDs and deterministic settings;
- capture or execution script path;
- expected output filenames and semantic assertions;
- screenshot paths and a short reason for each image;
- source document(s) being replaced and canonical destination;
- implementation and review status.

Contract tests will fail on an uncovered capability, duplicate canonical owner, unknown fixture,
unowned screenshot, or public chapter with no executable proof. Chapter Plan Gate approval is
recorded by committing the reviewed manifest; the Markdown plan remains the human-readable
rationale.

## 6. Core browser tutorial specification

### 6.1 Proposed English document outline

```text
# Create and export your first circular genome diagram

## What you'll need
## Step 1: Load the bundled mitochondrial genome
## Step 2: Generate the first diagram
## Step 3: Add a publication label
## Step 4: Make the feature map easier to read
## Step 5: Export the SVG
## What you built
```

The opening paragraph will show the finished figure and state the concrete outcome: a labeled
Circular SVG made from the bundled human mitochondrial GenBank record. It will not begin with
a feature list or an explanation of all available modes.

### 6.2 Fixed data and settings

| Field | Fixed value |
|---|---|
| Fixture | `HmmtDNA.gbk` |
| Accession/version | `NC_012920.1` |
| Sequence length | 16,569 bp |
| Existing file size | 64,640 bytes |
| SHA-256 | `7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f` |
| Mode | Circular |
| Input type | GenBank |
| Output Prefix | `human_mitochondrion` |
| Species | `<i>Homo sapiens</i>` |
| Track Preset | Middle |
| Separate Strands | On |
| Hide GC Content | Off |
| Hide GC Skew | Off |
| Label Mode | Out |
| Legend Position | Right |
| Standard export | SVG |
| Expected filename | `human_mitochondrion.svg` |

### 6.3 Step-by-step acceptance storyboard

| Step | Documented reader action | Visible result | Screenshot and reason | Playwright assertion |
|---:|---|---|---|---|
| 1 | Open the web app, select Circular if necessary, obtain the bundled fixture, and upload `HmmtDNA.gbk`. | The uploader shows a green `HmmtDNA.gbk` selection. | `01-input-ready.png`: confirms the correct file and uploader. | Use the accessible file input with `set_input_files`; assert the visible filename and no error. |
| 2 | Click **Generate Diagram** without changing advanced settings. | **Result Preview** displays the first Circular map. | `02-first-diagram.png`: the required early success result. | Click the real button; wait for the worker/result; assert one SVG, `NC_012920.1`, features, GC content, GC skew, and coordinate ticks. |
| 3 | Enter the output prefix and species, then click **Generate Diagram** again. | The preview center contains *Homo sapiens* and the output name is updated. | `03-publication-label.png`: verifies the literal species markup was interpreted correctly. | Fill labeled inputs; click Generate; assert control values and visible SVG text. |
| 4 | Choose Middle, keep Separate Strands on and both GC tracks visible, choose Labels Out and Legend Right, then regenerate. | External feature labels, two GC tracks, ticks, center definition, and right legend form the final figure. | `04-layout-settings.png` helps reproduce the choice; `04-finished-diagram.png` verifies the result and doubles as the opening hero image. | Select/check through roles or labels; click Generate; assert each DOM control and final SVG semantics. |
| 5 | Click **SVG**. | The browser saves `human_mitochondrion.svg`. | `05-export-svg.png` locates the export action; the actual file check is mechanical rather than photographic. | Wrap the real click in `expect_download`, save to a temporary directory, and validate filename, size, XML root, accession, species, feature groups, GC tracks, and absence of scripts/event handlers. |

Each screenshot will have concise alt text and a caption that tells the reader what to notice.
No screenshot will be added for a control merely because that control exists.

### 6.4 Optional material boundary

Only after `What you built`, add links to the following task guides:

- create a Linear diagram;
- choose a color palette or define feature colors;
- compare two or more genomes;
- add depth or other quantitative tracks;
- save and restore a session;
- export PNG, PDF, EPS, PS, or interactive SVG;
- inspect exact web defaults and input limits.

Do not include abbreviated versions of those procedures in the core tutorial.

### 6.5 Acceptance sketches for the remaining Tutorials

The later Tutorials use the same Diátaxis rules as T-GUI-01 but do not repeat its detailed
storyboard in prose. Their manifest must encode the following successful path before writing.

| Tutorial | Standard path and early result | Finished result and verification |
|---|---|---|
| T-GUI-02 Linear | Step 1 upload `NC_001416.gb` and select Linear/No comparison; Step 2 click Generate and see the first Lambda map; Step 3 add labels and a ruler. | Regenerate, inspect the labeled/ruler result, download the promised SVG, and assert record, coordinates, features, labels, axis, ruler, and legend. |
| T-GUI-03 LOSATN | Step 1 upload whole Lambda `NC_001416.1` and whole DE3 `NC_042057.1`; Step 2 generate a plain two-record Linear map; Step 3 choose browser LOSATN `megablast` with pinned serial/one-thread settings and generate comparison links. | Inspect one of the six links, download raw LOSAT TSV and SVG, and assert program/settings, both complete endpoints, six rows/links, popup spans, and file contents. |
| T-GUI-04 LOSATP | Step 1 upload all five annotated BGC records; Step 2 generate the plain five-record Linear map; Step 3 choose LOSATP Similarity groups with pinned settings and generate protein groups. | Inspect one group, export raw evidence/member sequences where exposed, export SVG, and assert all record IDs, CDS translations, non-zero groups, popup members, and files. |
| T-CLI-01 Circular | Step 1 enter a clean tutorial directory containing the bundled fixture; Step 2 run the fixed `gbdraw circular` command and obtain SVG. | Inspect the named file and figure; automated proof parses accession, record length, features, GC/skew, ticks, label, and legend. |
| T-CLI-02 Linear | Step 1 enter a clean tutorial directory containing the bundled Lambda fixture; Step 2 run the fixed `gbdraw linear` command and obtain SVG. | Inspect the named file and figure; automated proof parses record, features, labels, axis, ruler, and legend. |
| T-PY-01 Python | Step 1 create the minimal script around the bundled HmmtDNA fixture; Step 2 execute it and obtain SVG. | Inspect the figure and returned object/file; pytest executes the literal snippet and validates the same semantics as the published image. |

Each ends with `What you built`, then links to the relevant How-to and Reference chapters.
Comparison Tutorials explain the biological meaning of only the selected evidence type; the
comparison-method Explanation chapter owns cross-method choice and limitations.

## 7. Fixture and package-data work

Use one canonical root, `gbdraw/web/tutorial-data/`, for every small public tutorial fixture.
The CLI and Python chapters read the same repository files; the hosted and local GUIs serve
them at the same relative URL. Do not retain independent copies under `tests/test_inputs`,
Gallery files, and Tutorial folders after migration.

### 7.1 Bundled core fixture set

The first public bundle is intentionally small but covers Circular, Linear, all three browser
LOSAT programs, CLI protein modes, and GFF3 + FASTA. The frozen source files total 981,228
bytes before packaging, within the 1 MiB core budget.

| Fixture ID | Canonical files | Primary chapter ownership | Qualification before freeze |
|---|---|---|---|
| `human-mitochondrion` | `HmmtDNA.gbk`, `HmmtDNA_qualifier_priority.tsv`, small style/visibility rules | Circular GUI/CLI/Python, labels, shapes, GC/skew, sessions, export | Preserve `NC_012920.1`, 16,569 bp, file checksum, expected feature types/counts, and current visual quality. |
| `lambda` | `NC_001416.gb` | Linear GUI/CLI, ruler, regions, orientation | Preserve `NC_001416.1`, 48,502 bp, record/features, and deterministic label/ruler result. |
| `lambda-gff3` | `NC_001416.gff3`, `NC_001416.fna` | GFF3 + FASTA GUI/CLI/API | Preserve the complete natural `NC_001416.1` record without cropping or splitting; freeze 48,502 bp, gene/CDS counts, strands, translations, and checksums. |
| `de3` | `NC_042057.1.gb` | LOSATN/TLOSATX comparison input | Preserve the complete natural `NC_042057.1` record: 42,925 bp, 57 CDS, +37/-20 CDS strands, accession/version, and checksums. Do not substitute another record or a cropped derivative. |
| `lambda-de3-comparison` | `lambda-de3.losatn.tsv`, `lambda-de3.tlosatx.tsv` | Browser LOSATN/TLOSATX, uploaded Linear comparisons, and CLI/Python comparisons | Derive only from the two whole records with pinned browser WASM, serial scheduling, and one thread. Freeze LOSATN's six raw/default rows and TLOSATX's 397 raw rows, 266 default-displayed rows, seven tutorial-displayed rows at 1,000 bp minimum alignment length, arguments, endpoints, orientations, and checksums. |
| `aminoglycoside-bgc-five` | `BGC0000708.gbk`, `BGC0000709.gbk`, `BGC0000711.gbk`, `BGC0000712.gbk`, `BGC0000713.gbk`; color/qualifier TSVs | LOSATP Similarity groups, Collinear blocks, feature styling | Run pinned serial one-thread LOSATP across all five whole records. Record runtime version, thresholds, raw-result checksums, 155 displayed features, and qualified group/block counts. |

The BGC fixture qualification already indicates that a fixed one-thread run is practical for
PR checks: Similarity groups and Collinear modes complete on the order of seconds to
tens of seconds in the inspected environment. Those observations are planning evidence, not a
portable timeout contract; CI timeouts will be derived from repeated pinned runs.

### 7.2 Extended deterministic fixtures

| Fixture ID | Canonical or derived files | Use and policy |
|---|---|---|
| `tobacco-plastome-regions` | Complete `NC_001879.2` as `NC_001879.gbk`, `nicotiana-tabacum-regions.tsv`, `modified_default_colors.tsv`, and `qualifier_priority.tsv` | Region annotations and custom Circular slots. The four-file, 332,572-byte fixture preserves the complete 155,943 bp circular plastome and four named source-coordinate regions. |
| `depth-1kb` | Complete circular `AP027133.1` plus `AP027133.DRR394922.depth-1kb.tsv` | Depth/axis chapters. `tools/build_depth_1kb_fixture.py` reproducibly reduces the 606,194-row, 12.6 MB per-base source into 607 consecutive non-overlapping 1 kbp arithmetic means. Freeze source and derivative checksums, bin starts, min/max, and the plotting contract; do not ship the raw table in the tutorial bundle. The two packaged files total 1,361,007 bytes. |
| `metazoan-mitochondria-four` | Existing `HmmtDNA.gbk` (`NC_012920.1`) plus `NC_002333.2.gb`, `NC_024511.2.gb`, and `NC_001328.1.gb` | Circular grid examples for *Homo sapiens*, *Danio rerio*, *Drosophila melanogaster*, and *Caenorhabditis elegans*. Every input is the complete natural RefSeq mitochondrial genome and is annotated as circular. Never render a cropped, partial, or naturally linear record as Circular merely to demonstrate layout. The three additional records form a 167,108-byte extended-tier increment. |
| `metazoan-mitochondria-comparison` | `danio-human.tlosatx.tsv`, `drosophila-human.tlosatx.tsv`, and `caenorhabditis-human.tlosatx.tsv` | Circular conservation examples. Reproduce the three tables from the complete mtDNA records with pinned native LOSAT 0.1.0, record-specific mitochondrial genetic codes, serial one-thread execution, and two byte-identical runs. Freeze 435 raw rows, 106 displayed rows, 9,813 bp union coverage on the human subject, endpoints, orientations, and checksums. The three derived files add 30,552 bytes. |

Large AP027 depth tables, full Vibrio assemblies, nine-record Majanivirus data, the WSSV
20-series comparison, and large saved sessions remain Gallery/nightly fixtures. They are not
downloaded by beginner chapters and do not run in every pull-request capture. A Gallery chapter
may use them only with a committed raw-input or derivation recipe and a separate size/runtime
budget.

The committed extended tier is 1,891,239 bytes across 12 files: the three additional
metazoan mitochondrial records, three mtDNA TLOSATX tables, the tobacco plastome fixture,
and the AP027133 depth fixture.
It remains below the 5 MiB extended ceiling.

### 7.3 Manifest contract

Add one `gbdraw/web/tutorial-data/manifest.json`. Every raw or derived file records:

- fixture ID, relative path, media/input type, file size, and SHA-256;
- accession and sequence version where applicable;
- organism, record description, record IDs/order, sequence lengths, feature/CDS counts, and strands;
- authoritative source URL, retrieval date, license/provenance note, and any filename normalization;
  legacy files whose acquisition date cannot be recovered use `retrievedOn: null` plus
  `retrievalDateStatus: unknown-legacy`, while retaining the repository-added date and verified
  checksum; files added after 2026-08-04 require an ISO retrieval date;
- derivation script, source fixture/checksum, tool version, arguments, and output statistics;
- intended scenario IDs, package tier (`core`, `extended`, or `gallery-nightly`), and size budget;
- expected search/result semantics when the file is comparison evidence.

Fixture tests parse biological data rather than checking only bytes. They also reject duplicate
canonical copies and public links to `tests/test_inputs`. Update `_WEB_APP_PACKAGE_DATA` in
`gbdraw/_build_support.py`, source/wheel packaging tests, offline asset checks, and deployment
packaging so core fixtures are downloadable from both the hosted app and `gbdraw gui`.

The capture production gate accepts only `localhost` or `127.0.0.1` by default. Fixtures are
public reference data; no account, mutable database seed, or real user data is involved.
External requests are blocked during capture, and no tutorial step depends on a live NCBI,
MIBiG, DDBJ, or ENA response.

## 8. Executable documentation and Playwright capture implementation

### 8.1 Committed layout

```text
docs/
  scenarios/
    manifest.json
  capture/
    README.md
    config.py
    run_all.py
    assertions/
      downloads.py
      svg_semantics.py
      session_semantics.py
    flows/
      tutorials/
        gui_first_circular.py
        gui_first_linear.py
        gui_losatn.py
        gui_losatp_pairwise.py
      how_to/
        inputs.py
        layouts.py
        comparisons.py
        tracks.py
        styling.py
        interactive_sessions.py
        exports.py
  recipes/
    run_cli_scenarios.py
    run_python_scenarios.py
  TUTORIALS/
    GUI/
    CLI/
    PYTHON/
  HOW_TO/
    GUI/
    CLI/
    PYTHON/
  images/
    <scenario-id>/
      <step-id>-<purpose>.png
gbdraw/web/tutorial-data/
  manifest.json
  <fixture-id>/
    <raw-or-derived-files>
```

`python docs/capture/run_all.py` is the one-command core regeneration path. Its default/all mode
runs the selected cumulative GUI tier and every implemented CLI and Python recipe by directly
reusing their `run_scenario` APIs. `--scenario <id>` selects exactly one GUI, CLI, or Python
scenario, `--tier core|extended|nightly` filters GUI captures, and `--check` propagates to the
selected surface. Large Gallery case studies are nightly-only. `README.md` documents dependency
installation, browser/runtime installation, the local server, pinned environment, artifact
review, and expected runtimes.

### 8.2 Environment contract

Pin these values in `docs/capture/config.py`, not per screenshot:

- local base URL and isolation headers;
- Chromium/Playwright version used in CI;
- viewport and device scale factor;
- `en-US` locale, UTC timezone, light color scheme, reduced motion;
- vendored fonts and a no-external-network route policy;
- output image dimensions and file-size budget;
- maximum worker and generation timeouts.

LOSAT scenarios additionally pin program/mode, task, e-value/bitscore/identity/alignment-length
thresholds, serial versus parallel scheduling, runs in flight, threads per run, and the input
order. Screenshot byte identity is not used as the biological correctness oracle; raw search
and SVG semantics are asserted independently.

Add Python Playwright to the development dependency set and align its pinned browser version
with the Playwright version used by CI. A screenshot run from an unpinned browser is diagnostic,
not authoritative.

### 8.3 Selector and action policy

Before the smoke capture, connect visible labels to their controls with stable IDs or
`aria-label`. Use selectors in this order:

1. `get_by_role` with an accessible name;
2. `get_by_label`;
3. `get_by_test_id` only where a unique accessible selector is impractical.

Committed flows may not use bare CSS or text selectors without an explicitly recorded app
finding. `page.evaluate` may inspect SVG semantics or apply capture-only highlighting, but it
may not change form values, inject files, invoke `runAnalysis`, or call export functions.

The first smoke-capture accessibility changes cover:

- Circular mode;
- Circular GenBank file input;
- Output Prefix and Species;
- Track Preset and Separate Strands;
- Hide GC Content and Hide GC Skew;
- Label Mode and Legend Position;
- Generate Diagram;
- Result Preview;
- SVG download.

Before each later scenario cohort, run a selector census for the controls it owns. The extended
census must cover Linear record rows/regions/orientation, GFF3/FASTA inputs, comparison-plan
edges, LOSAT program and protein mode, thresholds, depth/annotation/custom-slot controls,
feature/legend editors, session/history actions, and every supported export button. Adding a
selector does not authorize adding that control to a Tutorial; chapter scope still follows the
reader goal in Section 5.

### 8.4 Surface-specific execution policy

| Surface/content | Executable authority | What counts as success |
|---|---|---|
| GUI Tutorial or How-to | Python Playwright flow from a fresh browser context | Every documented click/fill/select/upload/download occurs through the UI; promised preview, popup, state, and files are asserted. |
| CLI Tutorial or How-to | `docs/recipes/run_cli_scenarios.py` subprocess in a new temporary directory | Exact documented command exits as expected; stdout/stderr and every artifact satisfy the scenario contract. |
| Python Tutorial or How-to | `docs/recipes/run_python_scenarios.py` plus pytest import/execution | The literal documented code executes, returns the stated type, and creates the published artifact. |
| Reference examples | Generated help/schema extraction or focused contract test | Signatures, option names/defaults, tables, and version/capability matrices match implementation authorities. |
| GUI-generated CLI command | Playwright obtains/downloads the command and the CLI runner executes it | The handoff is syntactically valid and produces a semantically equivalent diagram or a documented surface-specific difference. |

No GUI scenario may preload a completed session, assign Vue/reactive state with JavaScript, call
an app method directly, or copy an existing screenshot. Sessions are permitted only in the
chapter that explicitly teaches save/load, after the state has first been built through visible
actions. CLI/Python chapters do not use terminal screenshots; their committed diagram assets
come from their executable recipes.

### 8.5 Workflow gates from `love-me-love-my-docs`

1. **Frame**: audience, English, browser, plain Markdown. Completed in this plan.
2. **Flow census**: the complete chapter list in Section 5. Approval is required before capture code.
3. **Demo data**: the public fixture suite, manifest, checksums, and local-only target. No authentication or
   mutable seed is needed.
4. **Smoke capture**: implement only Step 1 through Step 2, produce
   `02-first-diagram.png`, render the Markdown page, and confirm the image displays. Do not build
   the remaining flow until this passes.
5. **Capture harness**: complete T-GUI-01, then add T-GUI-02, comparison Tutorials, and the
   GUI How-to cohorts in dependency order. Each flow is cumulative within a chapter.
6. **Capture run**: generate the complete screenshot set, run CLI/Python recipes, and triage
   every failed action or unexpected artifact.
7. **Manual written**: write English chapters only around actions their harness has passed.
8. **Verify and report**: run all links, screenshots, recipes, semantic checks, fixture/package
   checks, and exports; document regeneration, runtimes, and CI behavior.

## 9. Existing documentation remediation

### 9.1 P0 correctness workstream

| Audit item | Required implementation | Exit condition |
|---|---|---|
| P0-1 session v40 | Correct release notes, use a stable release-specific statement, update the browser acceptance fixture, and add a version contract test. | Release notes, compatibility reference, implementation constant, and acceptance expectation agree. |
| P0-2 Tutorial 7 | Put every input in one preparation section and execute all blocks from an empty working directory. | Every documented output is generated in order without repository-root assumptions. |
| P0-3 Tutorial 9 | Make raw GenBank plus complete GUI/CLI actions the primary path; demote session import to comparison/shortcut status. | Gallery-quality figures reproduce without loading a session. |
| P0-4 Python API | Split minimal quick start from the finished example; make the documented Python code generate the published image; explain label selection. | Independent code blocks run and the final visual passes label/track QA. |
| P0-5 Tutorial 5 GFF3+FASTA | Replace the 87 bp toy as the main example with the existing accession-pinned Lambda fixture. | Records, CDS count, strands, checksum, output, and image are verified. |
| P0-6 Tutorial 5 similarity ring | Select an evidence-backed HSP threshold after overlap/coverage analysis and regenerate the figure. | Retained HSP count and union coverage are stated and multiple spans are visible at document width. |

These fixes should land before broad copyediting. Do not regenerate public figures by hand;
first establish the executable recipe and expected semantics.

### 9.2 Migrate current public content

Keep their paths in the first renovation release, but change the index and document contracts.

| Current file or surface | New canonical destination and reader goal |
|---|---|
| `QUICKSTART.md` | T-CLI-01. Replace live download drift with the fixed HmmtDNA fixture and one verified command path. |
| Gallery `HmmtDNA_basic_circular` | T-GUI-01. Replace session/state injection with the raw-input Playwright journey; keep Gallery as a linked showcase view. |
| Gallery `lambda_basic_linear` | T-GUI-02. Make “No comparison” and every final layout value explicit and executable. |
| `1_Customizing_Plots.md` | H-GUI-11/H-CLI-11 plus palette/shape Reference. |
| `2_Comparative_Genomics.md` | H-GUI-04/H-CLI-05 and comparison-method Explanation. Do not mix uploaded evidence, LOSAT search, and Circular rings in one procedure. |
| `3_Advanced_Customization.md` | H-GUI-11/H-GUI-12/H-CLI-11; exact rule precedence moves to Reference. |
| `4_Protein_Comparisons.md` | T-GUI-04, H-GUI-07/08, H-CLI-06/07/08, and comparison-program Reference. Pairwise, groups, and Collinear remain distinct goals. |
| `5_Table_Driven_Inputs.md` | H-CLI-02, H-GUI-06/10, and input/TSV Reference. Column contracts leave procedural prose. |
| `6_Depth_Quantitative_Tracks.md` | H-GUI-09/H-CLI-09/H-PY-03 plus quantitative-axis Reference. |
| `7_Linear_Layout.md` | H-GUI-03/H-CLI-04 after the P0 clean-directory repair. |
| `8_Interactive_SVG_Sessions.md` | H-GUI-13/14/15 and H-CLI-12/13; schema history remains only in Session Reference. |
| `9_Feature_Visibility_Shapes.md` | H-GUI-12/H-CLI-11, with raw input as the primary path and session load only as an optional shortcut. |
| `GFF3_FASTA.md` | H-GUI-01/H-CLI-01/H-PY-04 using the fixed whole-record Lambda pair. |
| `EXPORT.md` | H-GUI-15/H-CLI-13, output-format Reference, and publication Explanation. |
| `PYTHON_API.md` | T-PY-01 and H-PY-01 through H-PY-04; exact signatures/defaults stay in Python Reference. |
| `TYPED_API.md` | H-PY-05 plus typed request Reference. |
| `WORKFLOW_GUIDE.md` | The interface/layout/comparison Explanation chapters; remove `Web UX profile v1`. |
| `GALLERY.md` and Gallery JSON | Auxiliary showcase. An entry is called reproducible only when raw fixture, executable scenario, and artifact checks all exist. |

Each How-to guide must:

- start with `How to ...`;
- state prerequisites and exact input source/checksum;
- start from an explicit working state;
- use actionable steps, not a control inventory;
- state expected output and mechanical/visual verification;
- include troubleshooting for known failure modes;
- link to the relevant Reference page instead of duplicating tables or schema history.

Replace the current test that requires exactly numbered tutorials 1-9 with role-aware contracts:

- T-GUI-01 is listed first and optional capabilities follow it;
- GUI Circular, GUI Linear, GUI LOSATN, GUI LOSATP, CLI Circular, CLI Linear, and Python each
  have distinct Tutorial entries;
- each task guide is listed once under How-to;
- H1 labels match index labels;
- all new documents are reachable within two clicks;
- tutorial and how-to links point to Reference where exact contracts live;
- every Section 5.1 capability maps to a scenario-manifest owner;
- public navigation does not include `docs/internal/` or maintenance skills.

### 9.3 P1/P2 cleanup

After the P0 paths are executable:

- make Tutorial 2 and all remaining guides self-contained;
- bind documented commands, figure recipes, and images to one executable source or a tested
  semantic-equivalence contract;
- decide per Gallery entry whether it is a reproducible recipe or a curated showcase;
- replace live/current input downloads in Quickstart with stable fixtures or state expected drift;
- remove `Web UX profile v1` and normalize `hosted web app` to `web app` where deployment is not
  the subject;
- make Installation the authority for privacy and install paths;
- make Session Compatibility the authority for current compatibility;
- separate curated CLI semantics from generated CLI help;
- resolve PyPI/export-install inconsistencies;
- add navigation to long FAQ/Recipes pages and correct filename/typo mismatches;
- move completed plans, audits, and capture registers under `docs/internal/`.

## 10. Verification and CI

### 10.1 GUI capture acceptance

The GUI documentation set passes only when all of the following are true.

- Every scenario starts in a fresh browser context, reaches worker-ready state, and makes no
  unapproved external request.
- Playwright performs every documented GUI action through visible/accessibility-owned controls.
- Each GUI Tutorial displays its first diagram by Step 2 or Step 3 as declared in its manifest.
- Each screenshot is generated by its manifest-owned committed flow, appears at the documented
  path, has meaningful alt text/caption, and has a recorded decision/verification purpose.
- Every final figure is visually inspected at normal document width; comparison and track
  figures are also inspected at a scale where their evidence remains legible.
- T-GUI-01 and T-GUI-02 download real SVGs with the promised names and record/layout semantics.
- LOSATN, TLOSATX, and each LOSATP mode execute against pinned raw inputs. Raw TSV downloads,
  endpoint IDs, row/group/block counts, and diagram links agree with the frozen scenario.
- Uploaded-comparison scenarios prove that only the declared edges and result files are used.
- Session guidance saves a real current session, loads it in a fresh context, and reproduces
  semantically equivalent output without an undocumented search rerun.
- Export guidance captures every supported GUI download and checks its name, non-zero size,
  signature/parser, promised dimensions or embedded controls, and content.
- Exported SVG sanitization checks reject scripts, event handlers, and unsafe URLs where they
  are not part of the explicitly supported interactive-SVG contract.
- The selector audit reports zero unapproved CSS/bare-text selectors.
- The regeneration command succeeds from a clean checkout.

### 10.2 CLI and Python acceptance

- Every CLI scenario runs the exact documented command from an empty temporary directory with
  no undeclared repository-root files, network access, shell aliases, or prior outputs.
- Success, expected failure, overwrite refusal, filenames, and stdout/stderr are asserted.
- SVG outputs are parsed for records, features, tracks, matches, groups/blocks, axes, and labels;
  raster/vector outputs are validated by format-specific parsers or signatures.
- Pinned CLI LOSATP scenarios record the selected runtime/version and run actual Pairwise,
  Similarity-group, and Collinear searches with deterministic inputs and one-thread settings.
- Documented nucleotide/TLOSATX CLI workflows consume the frozen prepared tables and never
  imply that the CLI runs those searches itself.
- Every Python snippet is extracted/imported from its canonical executable source, runs as
  documented, and produces the same artifact used by the page.
- Simple and typed APIs are tested independently; compatibility adapters do not stand in for a
  documented current-interface example.

### 10.3 Documentation contracts

Add or update tests for:

- local Markdown targets and heading fragments;
- two-click discoverability from README;
- Diátaxis section membership and H1/index agreement;
- tutorial Step 2 first-result marker;
- scenario IDs, fixed values, screenshot paths, alt text, and captions;
- complete Section 5.1 capability ownership with no duplicate canonical chapters;
- fixture manifest/checksum, derivation, size tier, and source/wheel package membership;
- P0 session-version authority;
- clean-directory execution of documented CLI/Python blocks;
- semantic figure checks for record, feature, HSP, and track counts;
- surface capability claims, especially GUI versus CLI LOSAT behavior and export availability;
- compatibility landing pages that contain no duplicate runnable recipe;
- absence of public links to `tests/test_inputs` and internal plans.

### 10.4 CI jobs

Use three levels.

1. **Fast documentation contracts** run on every pull request: links/navigation, scenario and
   fixture manifests, CLI/Python snippet extraction, schema/help drift, and cheap semantic checks.
2. **Core executable documentation** runs T-GUI-01/T-GUI-02, CLI/Python Tutorials, and export
   checks when core docs, UI, renderer, fixture, font, packaging, or export files change.
3. **Extended comparison/track capture** runs LOSATN/TLOSATX/LOSATP, prepared comparisons,
   depth, annotations, editing, and sessions on relevant changes and on a scheduled job. Large
   Gallery rebuilds remain nightly/manual with their own runtime and artifact budget.

The capture job must fail on stale screenshots or failed actions. Intentional visual changes
require running `python docs/capture/run_all.py`, reviewing the images at readable scale, and
committing the regenerated assets.

Recommended local verification sequence after implementation:

```bash
python docs/capture/run_all.py --tier core
python docs/capture/run_all.py --tier extended
python docs/capture/run_all.py --check --tier core
python docs/recipes/run_cli_scenarios.py --all
python docs/recipes/run_python_scenarios.py --all
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
python tools/verify_gui_offline.py smoke-test
python -m pytest \
  tests/test_documentation_contracts.py \
  tests/test_documentation_scenario_contracts.py \
  tests/test_documented_recipes.py \
  tests/test_tutorial_documentation_contracts.py \
  tests/test_gallery_capture_contracts.py \
  tests/test_web_packaging.py -q
npm run test:playwright -- tests/web/gallery-tutorial.playwright.spec.js
```

Run the repository’s normal fast test and lint commands before handoff as required by the size
of the implementation diff.

## 11. Implementation order and gates

| Order | Work package | Gate to advance |
|---:|---|---|
| 0 | Confirm Section 5 chapter/feature census, staged-path architecture, priorities, and package-size budget. | Reviewed `docs/scenarios/manifest.json` is committed; no capture code before it. |
| 1 | Close P0-1/P0-2; establish fixture/derivation manifests, package contracts, accessible selectors, and executable-doc test skeletons. | Session authority and clean-directory contracts pass; core fixtures are present in hosted/local packages with checksums. |
| 2 | Build only the T-GUI-01 Step 1-2 smoke path and render `02-first-diagram.png` in Markdown. | Real upload/Generate succeeds, the image displays, and no state injection or external request occurs. |
| 3 | P1-A onboarding: finish T-GUI-01, T-GUI-02, T-CLI-01, T-CLI-02, T-PY-01, and standard SVG export. | All five onboarding journeys pass, show an early result, publish recipe-owned figures, and validate downloads. |
| 4 | P1-B comparisons: qualify the complete-genome nucleotide pair and the five-record BGC LOSATP fixture; add GUI LOSATN/LOSATP Tutorials and TLOSATX, Similarity-group, Collinear, uploaded-BLAST, Circular-ring, and CLI comparison How-tos. | All search/input modes execute on their actual supported surface; raw evidence and diagrams pass semantic checks. |
| 5 | P2-A data/layout: GFF3 + FASTA, multi-record Circular, Linear rows/regions/orientation, and all TSV manifest chapters. | Fixture identity/topology and clean-context GUI/CLI scenarios pass. |
| 6 | P2-B tracks/presentation: depth/GC/skew, annotations, custom slots, colors, labels, visibility, shapes, strokes, titles, legends, axes, and overlap behavior. | Each chapter has one focused outcome; mechanical checks and readable final figures agree. |
| 7 | P2-C interaction/reproducibility: preview editors/popups, sessions/history, all exports, Python/typed API, and generated-command handoff. | Fresh-context session round-trip, all file downloads, API recipes, and cross-surface handoff pass. |
| 8 | Migrate legacy pages, complete Reference/Explanation authority, classify Gallery entries, and clean public navigation/internal records. | No broken links, duplicate runnable authorities, uncovered capabilities, or public links to internal/test-only fixtures. |
| 9 | Run final visual/prose QA, packaging, browser, CLI/API, core/extended capture, and full relevant tests. | All Section 10 and Definition-of-done criteria pass from a clean checkout. |

For end-user prose changes, run the repository-required `avoid-ai-writing` detect pass, make
only targeted edits, and run its final verification after technical correctness is established.

## 12. Requirement traceability

| Requirement | Planned proof |
|---|---|
| Teach one standard end-to-end workflow first | T-GUI-01 is the first public learning path; optional layouts, comparisons, and surfaces are linked only after `What you built`. |
| Visible success within the first three steps | T-GUI-01 and each later Tutorial declare and test an early-result step; the first Circular diagram appears in Step 2. |
| Do not inventory every control | Each Tutorial/How-to fixes only values needed for one outcome; Web/CLI/Input/Comparison/Style References own exhaustive contracts. |
| Execute every documented browser action with Playwright | Every GUI chapter has a cumulative clean-context flow with role/label/test-id selectors and no state injection. |
| Execute CLI and Python actions on their real surfaces | Clean-directory subprocess scenarios and literal Python-code execution own those chapters; Playwright does not pretend to validate them. |
| Generate every screenshot from a committed script | `docs/capture/run_all.py` and manifest-owned flows own all GUI screenshots; CLI/Python diagrams come from committed recipes. |
| Use stable fixture files and deterministic settings | Multi-fixture packaged suite, provenance/derivation/checksums, pinned search/runtime/browser/display settings, and no live data dependency. |
| Use screenshots for decisions or verification | Every scenario manifest entry records a reason for each image; terminal and control-inventory screenshots are omitted. |
| Separate tutorial, how-to, and reference content | Role-based directories/navigation and canonical ownership contracts in Sections 4, 5, and 9. |
| Verify exported files | Playwright/subprocess download or output events plus filename, size, format parser/signature, semantic, and safety assertions. |
| Cover Linear, LOSAT, CLI, GUI, and all other major features | Section 5.1 coverage vocabulary maps to the Tutorial, GUI How-to, CLI How-to, Python, Reference, and Explanation chapter census; CI rejects gaps. |

## 13. Risks and mitigations

| Risk | Mitigation |
|---|---|
| The current uploader and form labels are not robust automation targets. | Fix accessible names before capture; use test IDs only for remaining ambiguity. |
| Pyodide startup makes captures slow or flaky. | Serve locally with required isolation headers, wait on worker readiness, use bounded stage-specific timeouts, and keep the core fixture small. |
| Running every LOSAT mode makes pull-request checks too slow. | Use the small BGC suite, serial one-thread deterministic settings, affected-scenario selection, and core/extended/nightly tiers; never replace real scheduled execution with session preload. |
| A proposed LOSATN/TLOSATX fixture produces empty or unstable evidence. | Make non-zero result, endpoint, row-count, checksum, and visual-legibility qualification a fixture-freeze gate before prose or screenshots are written. |
| Pixel output changes across Chromium/font versions. | Pin Playwright/Chromium and fonts; rely on semantic assertions plus reviewed image diffs rather than screenshot byte hashes alone. |
| The final mitochondrial labels are too dense at document width. | Inspect the generated image at normal width; adjust only the documented label/layout settings or fixture recipe, then recapture. |
| Markdown, Gallery JSON, and screenshots drift independently. | Give every action a stable scenario/step ID; make capture the executable authority; validate Markdown values and Gallery mappings against it. |
| Moving existing guide paths breaks external links. | Keep current public paths as tested one-release compatibility landing pages with one canonical destination and no duplicated runnable blocks. |
| Broad coverage turns individual pages into feature inventories. | Enforce one reader goal per Tutorial/How-to, split comparison modes where decisions differ, and move tokens/defaults/schemas to Reference. |
| GUI and CLI capability claims drift, especially for LOSAT and export. | Maintain a tested surface capability matrix; require an actual owner and execution harness per surface. |
| Tutorial fixtures inflate the browser package. | Track file and compressed package size by fixture tier; keep the 981,228-byte core set below its 1 MiB ceiling and the 1,891,239-byte deterministic extended set below 5 MiB. Keep raw 12–14 MiB depth tables and large Gallery data outside the tutorial bundle. |
| The tutorial accidentally depends on live NCBI or user data. | Ship the fixed public fixture, block external requests, and reject non-local capture targets by default. |
| A click succeeds but the export is missing or malformed. | Treat the download event and parsed file contents as release-blocking assertions. |

## 14. Definition of done

The renovation is complete only when:

- the P0 audit items are closed with their stated tests and visual checks;
- T-GUI-01 is the first English end-user path and follows its exact successful Playwright flow;
- Circular GUI, Linear GUI, LOSATN GUI, LOSATP GUI, Circular CLI, Linear CLI, and Python each
  have their planned Tutorial, and every Section 5.1 capability has a canonical chapter;
- every Tutorial meets its declared Step 2/Step 3 early-result gate;
- every GUI action in public procedural docs runs through Playwright, every CLI command runs in
  a clean temporary directory, and every Python snippet executes literally;
- every GUI screenshot and every CLI/Python figure regenerates from committed automation;
- actual SVG, interactive SVG, PNG, PDF, EPS/PS where supported, session, raw comparison TSV,
  and FASTA exports are captured and validated on their supported surfaces;
- the stable multi-fixture suite is packaged by tier for hosted/local GUI and CLI/API reuse with
  provenance, derivation, checksum, biological semantics, and size contracts;
- Tutorial, How-to, Reference, Explanation, Auxiliary, and Internal content are clearly routed;
- existing task guides no longer claim to be a sequential tutorial course;
- all public Markdown links and headings resolve and new content is reachable within two clicks;
- the tested surface matrix does not imply CLI LOSATN/TLOSATX execution or unsupported exports;
- CI detects uncovered features, broken actions/commands/code, stale screenshots/figures,
  missing fixture data, comparison drift, session drift, and export regressions;
- the final figures have been inspected at readable scale and the documented recipe produces
  the displayed result;
- the proposed commit title and English summary are included in the implementation handoff.

## 15. Current workflow status

```text
Docs Progress:
- [x] Step 1: Frame - audience, English, browser, plain Markdown
- [x] Step 2: Flow census - the 59-chapter GUI, CLI, Python, Reference,
      Explanation, and auxiliary census is approved in the scenario manifest
- [x] Step 3: Demo data - core and extended fixtures are accession-pinned,
      checksum-bound, packaged by tier, and qualified with whole biological records
- [x] Step 4: Smoke capture - the first Circular result is generated in Step 2
      from a fresh browser context without external requests or state injection
- [x] Step 5: Capture harness - committed GUI, CLI, and Python runners own the
      documented actions, figures, downloads, and semantic assertions
- [x] Step 6: Capture run - core, comparison/presentation, interaction, session,
      and export captures are generated and match fresh browser runs
- [x] Step 7: Manual written - canonical Tutorial, How-to, Reference, and
      Explanation chapters are published; retired guide URLs are compatibility routes
- [x] Step 8: Verify + report - browser captures, final visual/prose review,
      relevant test sets, and the implementation handoff are complete
```
