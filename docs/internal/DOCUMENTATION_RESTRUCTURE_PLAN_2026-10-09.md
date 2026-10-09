# Documentation restructure plan (v0.14.0)

Status: active
Created: 2026-10-09
Release target: 0.14.0 (does not hold the release without the Owner's word)
Supersedes: the public page structure, Tutorial surface policy, and per-page
writing in
[DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md](DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md).
Its four public routes, evidence architecture, and scenario manifest schema 2
stay in force.

## Objective

Make the public documentation and the Gallery easy for users to follow. A
reader should find one page per question, see the result first, and meet
contract details and data provenance only where they look for them.

The 2026-08-09 plan fixed the routes (Tutorials, Technical documentation, FAQ,
Gallery). Three problems remain:

1. Each Tutorial project has three pages (web app, command line, Python) that
   repeat the same introduction, inputs, result, and checks: 30 pages and three
   interface indexes for 10 projects. CLAUDE.md already says separate interface
   evidence does not justify separate pages.
2. Pages are written for an auditor. They open with revision pins and
   compatibility caveats ("A tracked image is not, by itself, a reproducible
   recipe."; satkey lists before Step 1) and name operations instead of results
   ("Multi-record LOSATP collinear blocks").
3. Seven pages exist only for old links ("retained temporarily"), and the
   604-line Session compatibility history sits beside the current reference.

The revised Vibrio Gallery entry is the first page in the new style.

## Writing rules for public pages

These rules apply to Tutorials, Gallery entries (`docs/GALLERY.md` and
`gbdraw/web/gallery/`), FAQ, the documentation home, and release notes.
Technical documentation follows rules 3 to 7.

1. **Lead with the question and the result.** The H1 or Gallery title names
   the outcome or the biological question, not the program or mode. The first
   paragraph says what the figure shows and what to look for. The result image
   comes before the steps.
   - Good: "Collinearity analysis of multi-replicon bacterial genomes
     (<i>Vibrio</i> spp.)", "Show where two phage genomes match".
   - Bad: "Multi-record LOSATP collinear blocks", "Compare Lambda and DE3 from
     the command line".
2. **One page per project.** A Tutorial page holds the shared introduction,
   inputs, and result once, then one section per interface in a fixed order:
   `## In the web app`, `## On the command line`, `## In Python`. A project may
   omit an interface it does not support.
3. **Provenance goes to a data section.** Input tables keep the local filename,
   record ID, length, and download link. Revision pins (`sat`, `satkey`),
   checksums, generator details, and source licences go to a closing
   `## About the data` section, or to
   [Get the tutorial inputs](../GETTING_TUTORIAL_DATA.md) when they apply to
   every Tutorial.
4. **Contract caveats go to Technical documentation.** A Tutorial, Gallery
   entry, or FAQ answer states the user-visible consequence in one sentence and
   links to the technical owner. Keep scientific qualifications that change how
   a reader interprets the figure (for example, that similarity links are not
   orthology calls), stated once, next to the figure.
5. **Write for the person doing the task.** Second person, present tense, one
   action per numbered step. Do not use maintainer terms in user prose:
   evidence, contract, owner, canonical, provenance, restoration, consume,
   authoritative, retained temporarily, reader-only, adapter.
6. **Exact names stay exact.** UI labels in bold as displayed, CLI options and
   file names in code, accessions and taxon names as in the source (genus and
   species in italics).
7. **Every visible result stays bound to its generator.** Screenshots come from
   capture flows, figures from the declared command or recipe, code blocks from
   the executable scenarios. Moving prose between pages does not regenerate
   artifacts; a changed figure or step does.

A Gallery entry has a title (rule 1), a description of at most three sentences
(what is shown, the setting that makes it work, what to look for), tags, and a
tutorial whose steps follow the web-gallery-screenshot-maintenance skill.

## Page inventory and dispositions

Reader question, owner, and disposition for every public page. `merge` names
the page that receives the content. Paths are under `docs/` unless shown.

| Page | Reader question | Disposition |
| --- | --- | --- |
| `README.md` (repository root) | What is gbdraw and where do I start? | keep; update links |
| `DOCS.md` | Where is the page for my question? | keep; rewrite route table question-first; drop removed pages |
| `QUICKSTART.md` | What is the fastest first figure? | keep (linked by the 0.13.0 README); router to the two first Tutorials |
| `INSTALL.md` | How do I install gbdraw? | keep |
| `GETTING_TUTORIAL_DATA.md` | How do I download and check Tutorial inputs? | keep; receives revision-pin rules shared by all Tutorials |
| `ABOUT.md` | How do I cite gbdraw? | keep |
| `TUTORIALS/README.md` | Which Tutorial should I follow? | keep; one row and one link per project |
| `TUTORIALS/{GUI,CLI,PYTHON}/README.md` (3) | Which Tutorials exist for my interface? | delete; the index states which interfaces each project covers |
| `TUTORIALS/{GUI,CLI,PYTHON}/<project>.md` (30) | How do I make this figure? | merge into `TUTORIALS/<project>.md` (10 pages, below) |
| `REFERENCE/README.md` and nine topic pages | What exactly does this control, option, schema, or API do? | keep; prose pass in PR-P |
| `REFERENCE/interactive-svg-and-semantic-hooks.md` | What can I query in an interactive SVG? | keep; receives `SVG_SEMANTIC_HOOKS.md` |
| `SVG_SEMANTIC_HOOKS.md` | Same question | merge into `REFERENCE/interactive-svg-and-semantic-hooks.md` |
| `CLI_Reference.md` (generated) | What are all CLI options? | keep (0.13.0 path) |
| `RECIPES.md` | Which command makes this common figure? | keep (0.13.0 path) |
| `SESSION_COMPATIBILITY.md` | Which session versions changed what? | move to `internal/SESSION_FORMAT_HISTORY.md`; the current reference keeps what a user does with an old file |
| `FAQ.md` | Which approach should I choose; why did this fail? | keep; prose pass in PR-P |
| `GALLERY.md` | What can a gbdraw figure look like, and how do I make one like it? | keep; list the Web Gallery entries with the same titles and order, then other figures |
| `PALETTE_EXPLORER.md` | How do I compare palettes? | merge: entry in `GALLERY.md`, colour-accessibility note in `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md` |
| `EXPORT.md`, `GFF3_FASTA.md`, `PYTHON_API.md`, `TYPED_API.md`, `WORKFLOW_GUIDE.md` | none (link routers) | delete; never in a release tag, targets already exist |
| `RELEASE_NOTES_0.14.0.md` | What changed for me in 0.14.0? | keep; rewrite in user language (PR-P) |
| `RELEASE_NOTES_0.14.0b0.md` | What changed in the beta? | keep as a record; not rewritten |
| `examples/color_palette_examples.md` (generated) | Which colours does each palette use? | keep |
| Web Gallery entries (`gbdraw/web/gallery/`, 10) | Show me a finished figure I can open and reproduce | keep; titles and descriptions per the writing rules (PR-V) |

### One-page Tutorials

Each project keeps its current slug, so links change only by the removed
interface folder. The H1 follows rule 1; the implementer may refine wording.

| Page `TUTORIALS/<slug>.md` | Merges | Interfaces |
| --- | --- | --- |
| `first-circular-genome-diagram.md` | T-GUI-01, T-CLI-01, T-PY-01 (`PYTHON/first-genome-diagram.md`) | web, CLI, Python |
| `first-linear-genome-diagram.md` | T-GUI-02, T-CLI-02, T-PY-03 | web, CLI, Python |
| `compare-genomes-losatn.md` | T-GUI-03, T-CLI-07, T-PY-04 | web, CLI, Python |
| `compare-proteins-losatp.md` | T-GUI-04, T-CLI-08, T-PY-05 | web, CLI, Python |
| `build-an-annotated-chloroplast-map.md` | T-GUI-05, T-CLI-06, T-PY-02 | web, CLI, Python |
| `add-precomputed-circular-comparison-rings.md` | T-GUI-06, T-CLI-09, T-PY-06 | web, CLI, Python |
| `compare-proteins-losatp-collinear.md` | T-GUI-08, T-CLI-10, T-PY-07 | web, CLI, Python |
| `create-and-resume-an-interactive-figure.md` | T-GUI-09, T-CLI-11, T-PY-08 | web, CLI, Python |
| `highlight-mitochondrial-features.md` | T-GUI-10, T-CLI-03, T-PY-09 | web, CLI, Python |
| `build-a-quantitative-genome-map.md` | T-GUI-12, T-CLI-05, T-PY-11 | web, CLI, Python |

Page layout:

```text
# <outcome title>
<question and result, 1-3 sentences>
<result image>
## Before you start        inputs table once; files each interface creates
## In the web app          numbered steps, screenshots
## On the command line     working directory, command, check
## In Python               program, run, check
## Next steps
## About the data          revision pins, sources, licences
```

Scenario IDs (`T-GUI-01` and so on), capture flows, recipe runners, and
artifact paths under `docs/images/` stay. Manifest `destination` values point
to the merged page plus its interface anchor (`#in-the-web-app`,
`#on-the-command-line`, `#in-python`); the executable markers
(`<!-- executable:<ID>:start -->`) are keyed by scenario ID, so three of them
can share one page. The manifest changes with the pages:

- `tutorial_projects.<project>` gains the page path and its H1; the page H1
  equals that project title, not the per-scenario `title`.
- `tutorial_project_policy.navigation` describes the interface sections instead
  of the "Choose how to build this figure" table.
- `destination` stays unique per scenario through its anchor; tests compare
  the file part with the page and check the anchor exists.
- `sources` entries that name a deleted page name the page that received its
  content.

### Gallery and Tutorials overlap

Five Web Gallery tutorials rebuild the figure of a docs Tutorial
(`HmmtDNA_basic_circular`, `lambda_basic_linear`, `tobacco-chloroplast`,
`BGC0000708-BGC0000713`, `hepatoplasmataceae_collinear`). For 0.14.0 both stay
and link to each other: the Gallery tutorial starts from the Gallery session in
the browser, the docs Tutorial starts from downloaded inputs on every
interface. Generating one from the other is deferred until after 0.14.0.

## Target sitemap

```text
README.md
docs/DOCS.md                       home: four routes, supporting pages
docs/QUICKSTART.md                 router to the two first Tutorials
docs/TUTORIALS/README.md           10 projects
docs/TUTORIALS/<slug>.md           x10, web app / command line / Python sections
docs/REFERENCE/README.md           + 10 topic pages
docs/CLI_Reference.md              generated option inventory
docs/RECIPES.md                    command templates
docs/FAQ.md
docs/GALLERY.md                    mirrors https://gbdraw.app/gallery/
docs/INSTALL.md  docs/GETTING_TUTORIAL_DATA.md  docs/ABOUT.md
docs/RELEASE_NOTES_0.14.0.md  docs/RELEASE_NOTES_0.14.0b0.md
```

Public Markdown pages under `docs/` fall from 64 to 33.

## Pull requests

Fewer, larger PRs grouped by files. Each PR lists the plan sections it
completes.

| PR | Contents | Depends on |
| --- | --- | --- |
| PR-T | This plan and the CLAUDE.md pointer; one-page Tutorials; `TUTORIALS/README.md`, `DOCS.md`, `QUICKSTART.md`, `README.md` links; delete the interface indexes and the five link routers; merge `PALETTE_EXPLORER.md`; manifest `destination`s, runners, and tests that name Tutorial paths; Gallery tutorial links to Tutorial pages | none |
| PR-V | Web Gallery in the new style: Vibrio entry rotated and regenerated with its new tutorial; titles and descriptions of the other nine entries; `GALLERY.md` rewritten to mirror the Web Gallery | Vibrio search captures after #954 (OV-211 names the record in search results) |
| PR-P | Release notes, FAQ, and Technical documentation prose in user language; `SVG_SEMANTIC_HOOKS.md` merged; `SESSION_COMPATIBILITY.md` moved to `internal/` | PR-T (links) |

PR-T and PR-V both touch `GALLERY.md` links; the second to merge rebases.

## Vibrio Gallery entry (`vibrio-harveyi-group-collinear`)

- Title: "Collinearity analysis of multi-replicon bacterial genomes
  (<i>Vibrio</i> spp.)".
- Description: "Two <i>Vibrio</i> genomes, each with two chromosomes. Each
  chromosome is rotated in gbdraw to start at its replication initiator gene —
  <i>dnaA</i> on chromosome I and <i>rctB</i> on chromosome II — so records
  that NCBI starts at unrelated positions line up. Inversions occur within each
  chromosome, but collinear blocks rarely connect chromosome I to chromosome
  II."
- Recipe change: the records table gains `display_start` (NC_004603.1 7,680;
  NC_004605.1 1,150; NC_022349.1 2,367,125; NC_022359.1 430,232, the start of
  each initiator CDS). Two specific colour rules (`CDS gene ^dnaA$` and
  `CDS product RctB`) with legend captions, and labels shown only on the four
  initiator CDSs. Everything else in the declared command stays.
- Tutorial story: generate a first preview without comparison; search a
  feature by qualifier value (`gene` = `dnaA`, then `product` = `RctB`); open
  each match; under **Layout** choose **Rotate record using this feature** and
  **Start of the record**; show its label; add the two colour rules; then run
  LOSATP **Collinear blocks** and generate once.
- Regenerate the session, interactive SVG, source figure, and thumbnail with
  the Gallery tools; recapture changed media with
  `tools/capture_gallery_tutorial_screenshots.py`.

## Verification

- Every PR: `tests/test_tutorial_documentation_contracts.py`,
  `tests/test_documentation_*`, `tests/test_documented_recipes.py`,
  `tests/test_onboarding_recipe_contracts.py`,
  `tests/test_gui_protein_comparison_capture_contracts.py`,
  `tests/test_python_tutorial_recipe_contracts.py`,
  `tests/test_web_packaging.py::test_gallery_tutorial_links_resolve`, and
  `tests/test_reproduce_examples.py::test_public_markdown_local_targets_exist`
  plus `test_public_figures_have_reproduction_inventory_coverage`. CI runs only
  the recipes job on a docs-only PR, so run the last three locally.
  `GALLERY.md` links written as HTML `<a href>` are not link-checked; check
  them by hand.
- PR-T: the CLI and Python recipe runners for the moved blocks; no
  screenshot or figure regeneration (rule 7).
- PR-V: Gallery JSON and capture checks from the
  web-gallery-screenshot-maintenance skill; the Vibrio session validation in
  `tools/refresh_gallery_sessions.py`; focused Gallery Playwright specs only.
- A detect-only `avoid-ai-writing` pass on changed public prose.
- Local Playwright runs stay targeted; CI runs the full suite.

## Stop conditions

- Do not drop a step, check, or runnable block when merging pages; move it.
- Do not change a figure or screenshot unless its step changed.
- If a project needs a second page to stay readable, record why here before
  adding it.
- Any change to what the product does (not how it is described) leaves this
  plan and follows the Product Impact ratchet.
