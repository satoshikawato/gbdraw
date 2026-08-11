# Tutorial Coverage Expansion Implementation Plan

| Field | Value |
| --- | --- |
| Date | 2026-08-04 |
| Status | Approved for implementation |
| Branch | `docs_renovation` |
| Chapter census | Supersedes [DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md](DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md) |

## Decision

Retired tutorials contain useful end-to-end workflows that the renovated documentation currently exposes mainly as isolated How-to or Reference material. Add seven project-oriented tutorials that teach those workflows without duplicating the focused How-to pages.

Implementation amendment (2026-08-04): omit the proposed Pairwise-only
`T-GUI-07` tutorial. Keep BGC protein comparison as Similarity groups
(all-vs-all), and teach Hepatoplasmataceae Collinear blocks from all-record
LOSATP evidence in `T-GUI-08`. Pairwise remains documented by its focused
How-to and Reference pages.

At the same time:

- remove the public `Tutorial fixture provenance` reference page and all public links to it;
- retain the machine-readable tutorial fixture manifest and its tests as internal documentation infrastructure;
- rename the existing LOSATP tutorial so that its title states that it teaches Similarity groups;
- regroup the Tutorials index by learning goal rather than only by interface;
- keep the nine Retired pages as short compatibility routes during the existing one-release compatibility window.

The implementation is documentation, deterministic fixture data, capture automation, recipe verification, and documentation-contract tests. It does not add product features or change rendering behavior.

## Goals

1. Restore the workflows that were useful in the Retired tutorials but are not yet taught as complete learning journeys.
2. Give every tutorial one concrete outcome and a visible first result no later than Step 2.
3. Make every GUI tutorial reproducible through the committed Playwright capture harness.
4. Make every CLI tutorial executable and semantically verified by the committed recipe runner.
5. Keep How-to and Reference pages as the canonical owners of individual options and schemas.
6. Remove fixture provenance from the reader-facing information architecture without weakening fixture integrity checks.
7. Preserve the Retired URLs as concise redirects to the most relevant current tutorial or How-to page.

## Non-goals

- Do not recreate the Retired pages verbatim.
- Do not turn Tutorials into exhaustive option catalogs.
- Do not duplicate complete How-to procedures inside Tutorials.
- Do not add terminal screenshots; CLI evidence is generated output plus automated assertions.
- Do not regenerate or modify `examples/gbdraw_social_preview.png`.
- Do not change diagram geometry or tracked reference SVGs.
- Do not expose fixture checksums, packaging details, or capture internals in public tutorial prose.

## Documentation model

The expanded documentation continues to follow the Diátaxis split:

| Content type | Responsibility after this change |
| --- | --- |
| Tutorial | A cumulative, guided project with a visible result and a defined finished figure. |
| How-to | The canonical, focused procedure for one task or control. |
| Reference | Exact CLI options, schemas, formats, and behavior contracts. |
| Explanation | Concepts, tradeoffs, and interpretation. |
| Internal infrastructure | Fixture provenance, checksums, capture recipes, and scenario verification. |

New tutorials must link to the relevant How-to and Reference pages for optional variations. They must not acquire `canonical_for` ownership in the scenario manifest; current How-to and Reference chapters retain that ownership.

## Target chapter census

Update [docs/scenarios/manifest.json](../scenarios/manifest.json) from 59 to 65 scenarios:

| Role | Current | Add | Remove | Target |
| --- | ---: | ---: | ---: | ---: |
| Tutorial | 7 | 7 | 0 | 14 |
| How-to | 33 | 0 | 0 | 33 |
| Reference | 11 | 0 | 1 | 10 |
| Explanation | 6 | 0 | 0 | 6 |
| Auxiliary | 2 | 0 | 0 | 2 |
| Total | 59 | 7 | 1 | 65 |

Remove `R-FIXTURE-01`. Keep the current IDs for all other chapters. Assign these IDs to the new tutorials:

| ID | Surface | Path |
| --- | --- | --- |
| `T-GUI-05` | GUI | `docs/TUTORIALS/GUI/build-an-annotated-chloroplast-map.md` |
| `T-GUI-06` | GUI | `docs/TUTORIALS/GUI/add-precomputed-circular-comparison-rings.md` |
| `T-GUI-08` | GUI | `docs/TUTORIALS/GUI/compare-proteins-losatp-collinear.md` |
| `T-GUI-09` | GUI | `docs/TUTORIALS/GUI/create-and-resume-an-interactive-figure.md` |
| `T-CLI-03` | CLI | `docs/TUTORIALS/CLI/highlight-mitochondrial-features.md` |
| `T-CLI-04` | CLI | `docs/TUTORIALS/CLI/build-a-table-driven-genome-comparison.md` |
| `T-CLI-05` | CLI | `docs/TUTORIALS/CLI/build-a-quantitative-genome-map.md` |

Set `plan_source` in the manifest to this plan. Record the earlier renovation plan in the approval note so the audit trail remains explicit.

## Retired coverage map

Use this table when updating the compatibility routes. The primary destination should receive the first link and explain the main project; secondary destinations cover remaining specialized controls.

| Retired route | Primary destination | Secondary coverage |
| --- | --- | --- |
| `1_Customizing_Plots.md` | T-GUI-05 annotated chloroplast map | T-CLI-03 feature presentation and the focused styling How-to pages |
| `2_Comparative_Genomics.md` | T-CLI-04 table-driven genome comparison | T-GUI-06 circular comparison rings and focused comparison How-to pages |
| `3_Advanced_Customization.md` | T-CLI-03 mitochondrial feature highlighting | Feature-label, color-table, shape, and visibility Reference pages |
| `4_Protein_Comparisons.md` | T-GUI-04 LOSATP Similarity groups | T-GUI-08 all-record Collinear blocks |
| `5_Table_Driven_Inputs.md` | T-CLI-04 table-driven genome comparison | TSV schema Reference |
| `6_Depth_Quantitative_Tracks.md` | T-CLI-05 quantitative genome map | T-GUI-05 and focused quantitative-track How-to pages |
| `7_Linear_Layout.md` | T-CLI-04 table-driven genome comparison | Focused linear-layout How-to |
| `8_Interactive_SVG_Sessions.md` | T-GUI-09 interactive figure and saved session | Interactive SVG and session How-to pages |
| `9_Feature_Visibility_Shapes.md` | T-CLI-03 mitochondrial feature highlighting | Focused feature-visibility and shape How-to pages |

## Tutorial specifications

### T-CLI-03: Highlight mitochondrial features without editing the GenBank file

Retired coverage: feature colors and labels; feature visibility and shapes; relevant styling controls from the circular styling tutorial.

Outcome: produce a mitochondrial map whose selected features are recolored, relabeled, filtered, and reshaped entirely through presentation tables and CLI options.

Inputs:

- the existing mitochondrial GenBank fixture used by the current CLI documentation;
- a tracked feature color table;
- a tracked feature presentation table containing label and shape overrides.

Learning sequence:

1. Run a minimal circular command against the raw GenBank record.
2. Open `mitochondrial_features_baseline.svg` as the first visible result.
3. Add qualifier-priority color rules and explain first-match precedence.
4. Add label overrides and a feature whitelist.
5. Demonstrate `show`, `off`, and `exclude_matching` visibility behavior.
6. Apply arrow, rectangle, and underlay shapes, including arrow geometry and strokes.
7. Render `mitochondrial_features_highlighted.svg` and compare it with the baseline.

Automated proof:

- both commands complete through `docs/recipes/run_cli_scenarios.py`;
- expected labels, colors, shapes, hidden features, and output files are asserted semantically;
- the final SVG contains no broken references or unexpected empty layers.

Extension links: focused feature styling, visibility, labels, color-table schema, and CLI option Reference pages.

### T-GUI-05: Build an annotated chloroplast map with custom tracks

Retired coverage: circular styling, custom tracks, quantitative tracks, and layout controls.

Outcome: build a publication-ready tobacco plastome diagram with LSC, SSC, IRa, and IRb annotations, GC content, AT skew, labels, title, and legend.

Inputs:

- the existing tobacco plastome GenBank fixture;
- a small tracked annotation table for LSC, SSC, IRa, and IRb.

Learning sequence:

1. Load the GenBank file from the input panel.
2. Render and inspect the default circular map.
3. Import or enter the four structural-region annotations.
4. Place the regions in deliberate custom slots.
5. Add GC content and AT skew, then tune the quantitative-track presentation.
6. Refine labels, title, and legend.
7. Export the finished SVG and verify the four regions and both numeric tracks.

Screenshots:

1. `01-input-ready.png`
2. `02-first-diagram.png`
3. `03-annotation-table.png`
4. `04-track-settings.png`
5. `05-finished-diagram.png`

Automated proof: the capture validates input identity, record length, region names and coordinates, track types, slot placement, visible title and legend, and successful SVG export.

Extension links: annotation tracks, quantitative tracks, feature styling, labels and legends, and circular layout How-to pages.

### T-CLI-04: Build a multi-row genome comparison from TSV manifests

Retired coverage: precomputed comparison links, TSV-driven inputs, and multi-record linear layout.

Outcome: arrange four majanivirus records in a two-by-two layout and connect the intended pairs with precomputed TBLASTX results.

Inputs:

- four complete majanivirus GenBank records;
- two precomputed TBLASTX tables;
- tracked `records.tsv` and `comparisons.tsv` files with explicit endpoint IDs.

Learning sequence:

1. Render the four records from `records.tsv` without comparison links.
2. Open `table_driven_comparison_baseline.svg` as the first visible result.
3. Read the record IDs, rows, columns, and orientation fields in the manifest.
4. Add `comparisons.tsv` with explicit source and target endpoints.
5. Apply identity and e-value thresholds.
6. Tune link curves, ruler placement, and the shared coordinate scale.
7. Render and inspect `table_driven_comparison.svg`.

Automated proof:

- both TSV files pass schema and referential-integrity checks;
- the four records occupy the expected cells and orientations;
- the two TBLASTX files connect only the declared endpoint pairs;
- the final SVG contains the expected comparison layers, rulers, and shared scale.

Extension links: table schema Reference, precomputed comparison How-to, and linear layout How-to.

### T-GUI-06: Add circular comparison rings from precomputed results

Retired coverage: precomputed comparison files, circular ring layout, and interactive match inspection.

Outcome: draw three circular similarity rings around the human mitochondrial reference and inspect one HSP interactively.

Inputs:

- the existing human mitochondrial GenBank record;
- the existing zebrafish, fruit fly, and nematode TLOSATX tables;
- deterministic FASTA exports for the three comparison records.

Learning sequence:

1. Load the human mitochondrial reference.
2. Render and inspect the default circular map.
3. Add the three precomputed tables and select the reference side.
4. Attach comparison FASTA so both HSP spans can be exported.
5. Arrange and style the three comparison rings.
6. Open an HSP popup and export the reference and comparison spans.
7. Export the final SVG.

Screenshots:

1. `01-input-ready.png`
2. `02-first-diagram.png`
3. `03-ring-settings.png`
4. `04-ring-result.png`
5. `05-hsp-popup.png`

Automated proof: validate the reference record, all three comparison tables, ring count and ordering, selected reference side, HSP metadata, both FASTA spans, and SVG export.

Extension links: circular similarity-ring How-to, comparison input Reference, and interactive SVG How-to.

### T-GUI-08: Find conserved gene order from all-vs-all LOSATP evidence

Retired coverage: collinear block discovery, anchors, orientation, colors, block inspection, and sequence export.

Outcome: search all record pairs for protein evidence, reduce the results to conserved-order blocks between adjacent display rows, and inspect one collinear block.

Inputs: five complete Hepatoplasmataceae GenBank records in Gallery order.

Learning sequence:

1. Load the existing Hepatoplasmataceae all-vs-all Similarity-group session.
2. Inspect the restored all-vs-all result and its five annotated records.
3. Select LOSATP Collinear, set evidence scope to All records, and set the minimum anchor count.
4. Reuse the 25 cached directional/self results from the fixed 32-thread all-vs-all run, reduce all ten record pairs, and interpret orientation and identity colors.
5. Adjust display thresholds and block presentation.
6. Open a block popup and export its envelope sequences.

Screenshots:

1. `01-input-ready.png`
2. `02-first-diagram.png`
3. `03-collinear-settings.png`
4. `04-collinear-result.png`
5. `05-block-popup.png`

Automated proof: assert the 25-entry all-vs-all cache, mode, all-record evidence scope, ten searched record pairs, no raw-search rerun, minimum anchors, block count, orientation and identity metadata, popup coordinates, envelope FASTA export, and final SVG export.

Extension links: Collinear How-to, interpretation Explanation, and export How-to.

### T-CLI-05: Build a quantitative genome map with depth, GC content, and skew

Retired coverage: read depth and numeric tracks, GC content/skew, axes, ticks, windows, and track slots.

Outcome: create a circular bacterial genome map combining depth, GC content, and GC or AT skew on readable, deliberately placed tracks.

Inputs:

- the existing `AP027133` GenBank fixture;
- the existing 1 kb depth table.

Learning sequence:

1. Render the annotated genome without numeric tracks.
2. Open `quantitative_genome_baseline.svg` as the first visible result.
3. Add the depth table and configure its axis and ticks.
4. Add GC content and choose window and step sizes.
5. Add GC or AT skew and explain the zero baseline.
6. Place the three tracks in explicit slots and style them for comparison.
7. Render and inspect `quantitative_genome_map.svg`.

Automated proof: assert record identity, depth row count and coordinate coverage, selected window and step, track types, slot assignments, axis and tick elements, and both output files.

Extension links: multiple and sparse tracks, log scaling, quantitative TSV schema, and focused quantitative-track How-to pages.

### T-GUI-09: Create an interactive figure and reproduce it from a saved session

Retired coverage: interactive SVGs, feature search, popup inspection, session download, session restore, and CLI handoff.

Outcome: create an interactive diagram, inspect a feature, save the complete working state, and reproduce the figure in a fresh browser context.

Inputs: the existing first-circular GUI tutorial fixture.

Learning sequence:

1. Load the GenBank record.
2. Render and inspect the default circular diagram.
3. Find a feature and open its popup.
4. Export an interactive SVG and verify that its controls work offline.
5. Download the compressed session file.
6. Start a fresh browser context, load the session, and regenerate the diagram.
7. Compare the restored result with the original and link to optional CLI `--session` use.

Screenshots:

1. `01-input-ready.png`
2. `02-first-diagram.png`
3. `03-interactive-export.png`
4. `04-feature-search.png`
5. `05-session-download.png`
6. `06-reloaded-result.png`

Automated proof: assert feature search and popup content, interactive SVG download, session gzip structure, successful load in a fresh context, equivalent regenerated state, and final SVG export.

Extension links: interactive SVG How-to, session How-to, session schema Reference, and CLI session handoff.

## Existing tutorial rename

Keep the path and scenario ID of `docs/TUTORIALS/GUI/compare-proteins-losatp.md` (`T-GUI-04`) but change its title from “Compare annotated proteins with LOSATP in the browser” to:

> Create protein Similarity groups with LOSATP in the browser

Update the H1, Tutorials indexes, scenario manifest, capture labels, fixture semantic description, and title-contract tests. Do not change its workflow or capture identity.

## Tutorials navigation

Rework `docs/TUTORIALS/README.md` into this reader path:

1. Start here
2. Compare genomes
3. Build a publication figure
4. Preserve and share work
5. Compatibility routes

The GUI, CLI, and Python indexes remain available for readers who choose by interface. Within each surface index, distinguish short beginner tutorials from multi-feature project tutorials. Each tutorial must appear exactly once in the scenario-defined order and use its exact H1 title.

Recommended grouping:

- Start here: first circular GUI, first linear GUI, first CLI diagrams, first Python figure.
- Compare genomes: LOSATN, Similarity groups, all-record Collinear blocks, precomputed circular rings, table-driven linear comparison.
- Build a publication figure: annotated chloroplast, mitochondrial feature highlighting, quantitative genome map, CLI styling.
- Preserve and share work: interactive SVG and saved-session tutorial.
- Compatibility routes: the nine Retired pages, placed last and labeled as compatibility links.

## Fixture strategy

### Keep fixture provenance internal

Delete `docs/REFERENCE/tutorial-fixture-provenance.md` and remove its links from reader-facing pages. Keep `gbdraw/web/tutorial-data/manifest.json` as the machine-readable authority for fixture names, checksums, accession metadata, and semantic expectations.

Retain and expand `tests/test_tutorial_fixture_manifest.py`. The fixture manifest remains an implementation contract, not a documentation chapter.

Each public tutorial should state only what a learner needs:

- recognizable file name;
- accession or record name where useful;
- expected sequence length or comparison relationship when it helps verify the step;
- a link to download the fixture bundle or named file.

### Add only missing reusable fixtures

Prefer existing tutorial data. Add data only where a browser workflow cannot be completed deterministically with the current bundle.

1. Promote the four majanivirus GenBank records and two TBLASTX tables used by the existing example into the tutorial-data fixture set, or copy them there with stable public names. Add `records.tsv` and `comparisons.tsv` with explicit endpoint IDs. Do not make the tutorial depend on mutable files under `examples/`.
2. Add deterministic FASTA files for zebrafish, fruit fly, and nematode comparison records so T-GUI-06 can export both HSP spans. Generate them once from the existing complete GenBank records and record their checksums in the manifest.
3. Reuse the current mitochondrial, tobacco plastome, `AP027133`, depth, protein-comparison, and first-circular fixtures.

For every added file, update:

- the tutorial-data manifest entry and semantic checks;
- package-data inclusion and offline bundle budget where applicable;
- fixture checksum and record-identity tests;
- the capture or CLI runner fixture resolver.

## Implementation architecture

### Scenario manifest first

Update the chapter plan before writing tutorial prose:

- add all seven IDs with `status: planned` and `canonical_for: []`;
- remove `R-FIXTURE-01`;
- update counts and ordering contracts;
- define `first_result_step: 2` for every new tutorial;
- set GUI execution evidence to Playwright capture and CLI evidence to recipe-runner verification;
- promote each scenario to `implemented` only when its prose, inputs, outputs, automation, and assets all exist;
- promote to `verified` only after the applicable tests and visual inspection pass.

### GUI capture flows

Add focused flows under `docs/capture/flows/tutorials/`:

- `gui_annotated_chloroplast.py`
- `gui_precomputed_circular_rings.py`
- `gui_losatp_collinear.py`
- `gui_interactive_handoff.py`

Register them in `docs/capture/config.py` and `docs/capture/run_all.py`, including screenshot budgets and fixture declarations.

Reuse behavior already exercised by the How-to flows, but do not make one scenario impersonate another. Extract small shared helpers from the existing annotation, quantitative-track, comparison, and session flows when doing so removes duplicated UI operations. Keep raw user actions visible: upload files, change controls, render, open popups, and download outputs through the UI rather than injecting final application state.

Every flow must:

- produce the first plain diagram at Step 2;
- validate fixture identity before interaction;
- wait for stable render completion rather than arbitrary delays;
- assert semantic state before taking each proof screenshot;
- save screenshots under `docs/images/t-gui-05` through `docs/images/t-gui-09`;
- leave no personal data, absolute local paths, or transient browser notifications in published images.

### CLI recipe runner

Extend `docs/recipes/run_cli_scenarios.py` for `T-CLI-03` through `T-CLI-05`.

- Register expected output files and required fixtures.
- Add generated presentation tables for T-CLI-03.
- Add generated or copied `records.tsv` and `comparisons.tsv` for T-CLI-04.
- Reuse the existing quantitative fixture plumbing for T-CLI-05.
- Parameterize existing semantic assertion helpers for feature presentation, quantitative tracks, precomputed comparisons, and linear layout instead of copying scenario-specific branches.
- Support a baseline command and a finished command in each new tutorial while preserving the “one executable bash block per scenario” contract.

The runner must validate output semantics, not merely file existence.

### Public tutorial prose

Every new page must contain:

1. a one-sentence outcome;
2. “What you will build” with the finished result;
3. prerequisites and exact input files;
4. numbered steps with copyable commands or exact GUI actions;
5. a visible first result by Step 2;
6. short “What changed” explanations after meaningful edits;
7. verification cues the learner can see in the output;
8. “What you built” summarizing the completed figure;
9. links to the relevant focused How-to and Reference pages;
10. the deterministic output file names expected by automation.

Use concise English prose consistent with the current documentation. Avoid fixture provenance language, internal scenario IDs in reader-facing paragraphs, and exhaustive option listings.

### Retired compatibility routes

Keep the current nine files under `docs/TUTORIALS/` as short compatibility pages. Update their destination links so that each Retired workflow points first to the closest new tutorial, then to focused How-to or Reference material for remaining details.

Compatibility pages must remain below the existing length limit, contain no executable code blocks, and must not be included among the primary learning paths.

### Remove the public provenance page safely

Remove links before deleting the file. Audit at least:

- `README.md`
- `docs/DOCS.md`
- `docs/GALLERY.md`
- `docs/GFF3_FASTA.md`
- GUI How-to pages that currently cite fixture provenance
- `docs/REFERENCE/README.md`
- `docs/REFERENCE/input-formats-and-tsv-schemas.md`
- `docs/scenarios/manifest.json`
- documentation reference and link-contract tests

Replace links only when the surrounding sentence still needs a learner-facing destination. Otherwise remove the sentence rather than substituting an internal manifest link.

## File-level change inventory

### Add

- Seven tutorial Markdown files listed in the chapter census.
- Four GUI tutorial capture flows.
- GUI screenshot assets under four new scenario directories.
- New tutorial fixture files required for the majanivirus comparison and comparison-span FASTA export.
- Focused contract tests only where existing domain tests cannot express the new invariant cleanly.

### Modify

- `docs/TUTORIALS/README.md`
- `docs/TUTORIALS/GUI/README.md`
- `docs/TUTORIALS/CLI/README.md`
- `docs/TUTORIALS/GUI/compare-proteins-losatp.md`
- the nine Retired compatibility pages
- `docs/scenarios/manifest.json`
- `gbdraw/web/tutorial-data/manifest.json`
- `docs/capture/config.py`
- `docs/capture/run_all.py`
- shared capture helpers as required
- `docs/recipes/run_cli_scenarios.py`
- public pages that link to fixture provenance
- scenario, tutorial, capture, recipe, fixture, and documentation-reference contract tests

### Delete

- `docs/REFERENCE/tutorial-fixture-provenance.md`

Do not delete the fixture manifest or its integrity tests.

## Test plan

### Static documentation contracts

Update or add tests that assert:

- 65 total scenarios and the target role counts;
- all 14 tutorial IDs in manifest order;
- all new tutorial paths exist and have one matching H1;
- every new tutorial has `first_result_step == 2`;
- GUI scenarios name registered capture flows and screenshot budgets;
- CLI scenarios name registered recipe-runner evidence;
- `R-FIXTURE-01` and `tutorial-fixture-provenance.md` no longer occur in public documentation;
- Tutorials indexes list each tutorial exactly once;
- Retired pages remain short and contain no runnable code;
- every relative Markdown link resolves.

### Fixture contracts

Verify:

- checksum, size, accession, record count, and sequence length for every new biological file;
- comparison tables reference the intended query and subject records;
- `records.tsv` IDs are unique;
- `comparisons.tsv` endpoints resolve to declared record IDs;
- comparison FASTA sequences match their source GenBank records;
- the packaged offline tutorial-data bundle remains within its configured size budget.

### CLI execution

Run the three new scenarios individually, then the entire CLI documentation recipe suite. Assert semantic SVG structure and deterministic output names. Do not update tracked rendering reference outputs.

### GUI capture

Run the four new flows individually before running the complete capture suite. Inspect all final screenshots at readable scale and confirm that:

- text is legible;
- no controls obscure the figure;
- the diagram demonstrates the promised feature;
- popup screenshots show the intended record, match, feature, or block;
- download and restore assertions verify actual files rather than toast messages alone.

### Regression suite

At minimum, run:

```bash
python -m pytest \
  tests/test_documentation_scenario_contracts.py \
  tests/test_tutorial_documentation_contracts.py \
  tests/test_documentation_reference_contracts.py \
  tests/test_documentation_capture_contracts.py \
  tests/test_tutorial_fixture_manifest.py \
  tests/test_gui_nucleotide_comparison_capture_contracts.py \
  tests/test_gui_protein_comparison_capture_contracts.py \
  tests/test_gui_tracks_capture_contracts.py \
  tests/test_gui_interactive_capture_contracts.py \
  tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py \
  tests/test_cli_comparison_how_to_recipe_contracts.py \
  tests/test_onboarding_recipe_contracts.py -q
```

Then run the new CLI scenarios, new GUI captures, the full documentation-link check, and the broader fast test suite if the fixture or runner refactor touches shared behavior.

## Implementation phases and gates

### Phase 0: Approve the chapter contract

1. Add the seven planned scenarios and remove `R-FIXTURE-01`.
2. Update count, ordering, title, and evidence contracts.
3. Record the T-GUI-04 title change.

Gate: the scenario manifest validates, and no tutorial prose is started before paths, outcomes, inputs, outputs, and evidence are agreed.

### Phase 1: Prepare deterministic inputs

1. Add only the missing majanivirus and comparison FASTA fixtures.
2. Update checksums, semantic metadata, package inclusion, and budgets.
3. Add fixture-integrity and endpoint-integrity assertions.

Gate: all fixture tests pass offline and each new input has a verified biological identity.

### Phase 2: Refactor shared verification helpers

1. Extract the smallest reusable capture helpers needed by multiple GUI scenarios.
2. Parameterize CLI semantic assertions used by old and new recipes.
3. Keep existing scenario behavior unchanged.

Gate: current capture and recipe contract tests pass before new scenarios are added to the runners.

### Phase 3: Implement CLI tutorials

Implement T-CLI-03, T-CLI-04, and T-CLI-05 with baseline and finished outputs. Add prose only after commands and semantic checks execute successfully.

Gate: each tutorial can be copied into a clean temporary directory and completed using only declared inputs.

### Phase 4: Implement GUI tutorials

Implement T-GUI-05, T-GUI-06, T-GUI-08, and T-GUI-09 one at a time. For each scenario: build the flow, make assertions pass, capture images, visually inspect them, then write prose against the recorded action sequence.

Gate: every numbered screenshot corresponds to a deterministic captured state and the finished image clearly demonstrates the tutorial outcome.

### Phase 5: Integrate navigation and compatibility routes

1. Apply the T-GUI-04 title change everywhere.
2. Rebuild Tutorials indexes by learning goal and surface.
3. Update Retired routes to the new canonical destinations.
4. Add cross-links without duplicating focused How-to procedures.

Gate: index-order, exact-title, uniqueness, compatibility-page, and relative-link tests pass.

### Phase 6: Remove fixture provenance from public docs

1. Remove or rewrite every public citation to the provenance page.
2. Remove `R-FIXTURE-01` from remaining tests and inventories.
3. Delete `docs/REFERENCE/tutorial-fixture-provenance.md`.
4. Confirm that the internal manifest and fixture-integrity tests still cover provenance.

Gate: repository search finds no public reference to the deleted page, and the fixture suite passes unchanged in purpose.

### Phase 7: Full verification and status promotion

1. Run static contracts, fixtures, CLI scenarios, GUI flows, and link checks.
2. Inspect public diffs separately from automation and test diffs.
3. Inspect all generated tutorial figures at readable scale.
4. Change scenario statuses from `planned` to `implemented` and finally `verified` only when their evidence exists.

Gate: all Definition of Done items are satisfied in one internally consistent change set.

## Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Tutorials duplicate How-to material and become costly to maintain. | Keep the tutorial narrative project-oriented, link out for option matrices, and leave `canonical_for` empty. |
| Seven additions overwhelm the Tutorials index. | Group by learning goal, distinguish beginner and project tutorials, and place compatibility routes last. |
| GUI captures become flaky or depend on hidden state. | Use explicit UI actions, semantic waits, fresh browser contexts, and existing tested helpers. |
| Fixture growth makes offline packaging too large. | Reuse current files, add only missing records, record sizes, and enforce the package budget. |
| Biological inputs become ambiguous. | Assert accessions, lengths, checksums, record IDs, and comparison endpoints in the internal manifest tests. |
| Removing the provenance page creates broken links. | Remove links first, run repository-wide search and Markdown link checks, then delete the page. |
| Shared-runner refactoring changes existing documentation behavior. | Refactor before adding scenarios and require all current contract tests to remain green at the phase gate. |
| Baseline and finished commands drift apart from prose. | Execute the exact fenced block through the scenario runner and give both outputs deterministic names. |

## Definition of Done

The work is complete only when:

- all seven new tutorials exist at the approved paths and achieve their stated outcomes;
- T-GUI-04 uses the Similarity groups title everywhere while retaining its path and ID;
- the scenario manifest contains 65 chapters with the target role counts;
- each new tutorial produces a visible result by Step 2;
- all three new CLI recipes execute and pass semantic output checks;
- three new GUI flows run from raw inputs, the Hepatoplasmataceae flow starts from the validated All-vs-all Gallery session and reuses its cached evidence, and all four produce their declared screenshots and pass visual inspection;
- the Tutorials indexes present a coherent learning path and list every tutorial once;
- the nine Retired pages remain concise compatibility routes with updated destinations;
- `docs/REFERENCE/tutorial-fixture-provenance.md` and its public links are gone;
- `gbdraw/web/tutorial-data/manifest.json` and fixture-integrity tests remain the internal provenance authority;
- all fixture, documentation, link, recipe, and capture contract tests pass;
- no unrelated dirty files, generated distribution artifacts, tracked reference SVGs, or owner-maintained showcase images are changed.

## Proposed implementation commit

Title:

`docs: restore advanced tutorial workflows`

Summary:

- add seven reproducible project tutorials covering styling, comparisons, quantitative tracks, layout, and saved sessions;
- expand deterministic fixture, capture, and recipe verification for the new workflows;
- reorganize tutorial navigation and retire the public fixture-provenance chapter while retaining internal integrity checks.
