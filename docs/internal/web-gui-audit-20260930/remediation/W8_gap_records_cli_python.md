<!-- Raw design report of workstream W8 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

All paths are under DEV=/tmp/claude-1000/-mnt-c-Users-genom-GitHub-gbdraw/c9538259-b3ce-4057-baa8-e783da1943ee/scratchpad/dev, at commit 4c89bab1. I ignored the tests/web/audit-*/ directories. I did not modify anything. The only Python I ran was read-only JSON inspection and one in-memory pandas check.

# Class 4: JS fast record discovery vs Python loading

## IN-02 (Prokka empty ACCESSION/VERSION gives record ID 'KEYWORDS')

Nearest tests:
- **`tests/web/record-selector.test.mjs:84-87`**: asserts `parseSequenceRecordText` gives `A1.2` (VERSION wins) and `B1` (ACCESSION only).
- **`tests/web/record-selector.test.mjs:93-99`**: asserts the fast path returns `FAST.1` and makes 0 Worker calls.
- **`tests/web/record-display-options.test.mjs:23-27`**: LOCUS-only chunks; checks topology only.
- **`tests/web/record-metadata-inference.test.mjs:28-97`**: organism, strain and inferred definition. All its chunks have filled `ACCESSION`/`VERSION` (for example `NC_000001.1`).
- **`tests/web/record-display-discovery.playwright.spec.js:19-58`**: the only fast-path vs Worker parity check. At line 51 it compares `discoverSequenceRecords` with the Worker's `list_sequence_records` (Biopython, `gbdraw/web/js/app/python-helpers.js:1604-1643`) using `expect(fast).toEqual(worker)`.
  - Its fixture (lines 4-17) is synthetic, with `ACCESSION   same` / `VERSION     same` filled in.
- **`record-display-discovery.playwright.spec.js:96-111`**: uploads `MjeNMV.gbk` as `.ddbj`. It checks only the fast path against a hard-coded `LC738868.1` and asserts the Worker stayed idle, so there is no parity check here.
- **`tests/test_record_metadata.py:18-36`** with **`tests/fixtures/record_metadata_inference_cases.json`**:
  - The shared vectors have keys `comment` and `definition` only (8 cases, organism/strain formatting).
  - They contain no record-ID vectors.
- **`tests/test_genome_loading.py:20-77`**: its GenBank files are written by `SeqIO.write`, which always fills ACCESSION and VERSION.

What they miss:
- There is no Prokka-style fixture anywhere. A grep of `tests/`, `examples/`, `gbdraw/web/gallery/` and `docs/` for `^(ACCESSION|VERSION)\s*$` finds nothing.
  - The MAG files have placeholder `ACCESSION   ########`, which is not empty (`tests/test_inputs/MAGs/F1.ddbj:3-4`).
- No record-ID shared vector exists between JS and Python.
- The one parity test uses a single synthetic, fully populated header.
- It is Playwright, not tagged `@pr-smoke`, so it runs only in `playwright-functional` on push to dev.

## IN-03 (GFF3+FASTA fast path lists featureless FASTA records)

Nearest tests:
- **`record-display-discovery.playwright.spec.js:113-122`**: `NC_013668.gff3`/`.fasta` gives `[['#1','NC_013668.3']]` and the Worker stays idle.
- **`tests/web/record-selector.test.mjs:121-141`**: the GFF pair goes through a mocked helper (`GffOwned`). It checks only the payload roles.
- **`tests/web/record-display-options.test.mjs:299-320`**: discovery status only.
- **`tests/test_gff_multipart.py:10-39`**: two records, both with GFF features. It asserts FASTA order is kept.
- **`tests/web/linear-multi-record.playwright.spec.js:2800-2830` and `:3456-3470`**: multi-record GFF/FASTA, but every FASTA ID has a GFF feature line.
- **`tests/test_web_error_adapter.py:192`**: the inverse case (GFF record without FASTA gives an error).

What they miss:
- There is no fixture where the FASTA holds a sequence with no GFF features.
  - `tests/test_inputs/NC_013668.fasta` and `gbdraw/web/tutorial-data/lambda-gff3/NC_001416.fna` each have 1 header, and it has features.
  - Every inline GFF test fixture has features for every FASTA ID.
- No test compares against `merge_gff_fasta_records` (`gbdraw/io/genome.py:105-135`).
- A JS-vs-Worker parity test would not catch this either: the Worker helper `list_gff_fasta_records` (`python-helpers.js:1645-1660`) also lists every FASTA record ("List every FASTA record available…"). Only a comparison against `load_gff_fasta` would.

**Parity inventory:** the only JS-vs-Python record-listing comparison is `record-display-discovery.playwright.spec.js:40-51` (GenBank). Nothing compares with `load_gbks`/`load_gff_fasta`.

# Class 8: Web vs CLI compatibility

## GE-03 (Source recipe does not reproduce the diagram)

Nearest tests:
- **`tests/test_run_info_exact_replay.py:36-201`** (marker `browser`):
  - One fixture: `HmmtDNA_basic_circular.issue-469.json.gz` (line 25). It is circular, `grouping single`, 1 record.
  - Runs both `sourceRecipe` and `exactReplay` through the CLI.
  - Pins the SVG SHA (lines 26-33, 178) and does a `compare_svgs` against the browser SVG (lines 196-201).
  - The fixture's overrides are `objects.definition.circular.font_size=18` and `plot_title_font_size=32`. There is no interval, scale or ruler override.
  - 18 and interval 20 are the `config.toml:202-203` defaults, so the CLI's derived `int(18+2)=20` (`gbdraw/circular.py:951-954`) coincides with the Web default.
- **`tests/web/run-info.test.mjs:175-196`**: builds recipes for two sessions and runs them through the CLI.
  - `HmmtDNA_ATskew`: font 28 with an explicit interval 30, which is consistent with font+2.
  - `lambda_basic_linear`: no `objects.scale.font_size` or `ruler_label_font_size`.
  - The harness at `:116-172` asserts only exit code 0 and SVG size > 0. **It does not compare SVGs.**
- **`run-info.test.mjs:388-433`**: slot recipes (`ticks:ticks`, `tick_label_layout`). It checks the parser accepts them and they run.
- **`tests/web/joint-display-placement.playwright.spec.js:363-400`** with **`tests/web/helpers/joint-replay.py:20-73`**:
  - Uses Web defaults and a single record.
  - Compares Source-recipe SVG to Exact-replay SVG, both produced by the CLI, not against the Web preview.
- **`tests/test_circular_feature_width.py:620-671`** (CLI `--definition_font_size 12.5` gives interval 14) and **`tests/test_linear_track_layout.py:1491-1511`** (CLI ruler falls back to `--scale_font_size`):
  - These pin the CLI-only coupling. No test sets this against the Web request, where `session-request.js:1163-1165` sends font and interval separately.
- **`tests/web/linear-typography.test.mjs:18-46`**: only tests Web-state linking of scale and ruler font sizes.

What they miss:
- No case has a non-default definition font with the interval omitted.
- No case has `scale_font_size` set with `ruler_label_font_size` unset.
- The SVG-equality oracle covers exactly one circular single-record default-font fixture. There is no linear or grid SVG comparison.

## SE-06 (CLI Linear+BLAST session opens read-only in Web)

Nearest tests:
- **`tests/web/session-cli-compatibility.test.mjs:47-52` and `:59-115`**: four cases (single circular, composite circular, `linear --gbk mito lambda`, gff).
  - Line 79 asserts `importedComparisonIntent.disposition === 'EDITABLE'`, but no case has comparisons, so the check is trivially true.
- **`session-cli-compatibility.playwright.spec.js:17-22` and `:107-172`**: the same four cases, plus Web Generate and bidirectional replay. No `--blast` or `--comparisons_table`.
- **`tests/web/imported-comparison-intent.test.mjs:15-24` and `:78-88`**: synthetic records `recordKey 'record-a'/'record-b'`. Committed and candidate keys always match.
- **`tests/web/session-draft-authority.test.mjs:371-516`**: BGC gallery session with Web-style `record-1..5` keys.
- **`tests/fixtures/sessions/test_linear_cli_sidecar_reuses0.v40-schema6.json.gz`**: the only CLI-written linear sidecar.
  - It has uid `cli-seq-1`, 1 record and a `generatedProteinComparison`.
  - It is used only for promotion checks (`tests/web/joint-display-placement.test.mjs:167-179`, `tests/test_joint_display_placement_surfaces.py:121`).
- **`tests/test_session_compat.py`** (for example `:1064`): Python only.

What they miss:
- The string `cli-seq` appears in 0 tests.
- There is no CLI Linear session with a nucleotide BLAST import followed by Web Generate.
- The Playwright spec is not `@pr-smoke`, so it runs only on dev push.

## SE-07 (CLI session legend position lost)

Nearest tests:
- The session-cli-compatibility pair above never passes `--legend`.
  - The CLI default `right` (`gbdraw/circular.py:309-312`, `gbdraw/linear.py:646-651`) equals the projection fallback at `session-request.js:4410`.
  - The Playwright `svgSemantics` (lines 36-47) compares record IDs, feature d/fill and text contents only, so legend transforms are invisible to it.
  - After Generate it checks only `after.request.records` length (line 152), not `output.legend` or `form.legend`.
- **`session-draft-authority.test.mjs:438-516`**: a Web session with `ui.layoutPreferences`. It asserts the editor preference `'bottom'` overrides the request's `'right'` (lines 502-506).
- **`depth-track-session.playwright.spec.js:1656-1698`**: `config.form.legend='right'` wins over `output.legend='bottom'`.
- **`session-losat-cache-validation.test.mjs:181, 249, 411, 655`**: Web sessions carrying `config.form`.

What they miss: no test imports a CLI (config-less, ui-less) session with a non-default `--legend`, and none asserts `state.form.legend` or the regenerated legend position.

## TR-07 (slot legend label with ',' or ' #' breaks the Run Info CLI command)

Nearest tests:
- Every `legend_label` in the tests is letters and spaces only: "Genes", "AT skew", "Depth 2", "Sample A/B", "Selected Sample B", "Reviewed region".
  - Locations: `session-request.test.mjs:122-188, 3606-3634`; `test_circular_track_slots.py:384-392, 1854`; `test_linear_track_slots.py:656`; `test_cli_tables.py:259-268`; `test_depth_track.py:1549-1585`; `test_api_session.py:275-358`.
- **`session-request.test.mjs:143-157`**: round-trips through the JS builder and **JS** parser only.
- **`track-slot-validation.test.mjs:723-760`**: JS pixel-field round-trip corpus; it has no text fields.
- **`run-info.test.mjs:388-433`**: no `legend_label` in recipe slots.

What they miss:
- There are no tests for `split_kv_list` or `strip_inline_comment` (`gbdraw/tracks/parsing.py:118-146`).
- No test sends a Web-built slot spec (`circular-track-slots.js:1149`, `linear-track-slots.js:438-439`) through the Python parser.

## CO-07 (Save Raw LOSAT TSV for a reverse-complemented record)

Nearest tests:
- **`tests/web/run-analysis-simple-path.test.mjs:1634-1720`**: nucleotide LOSAT with `region_reverse=true`.
  - Asserts the conversion `queryViewTransform {length:8, reverse:true}`.
  - Asserts the helper-zip BLAST entry contains the substring `MIDDLE\tTHIRD`, which comes from a mocked Worker TSV.
  - It never calls `downloadLosatPair`.
- **`tests/web/losat-cache-migration.playwright.spec.js:81-140`** and **`tests/run_losat_cache_browser_acceptance.py:426-485`**:
  - Assert 12 columns, alias IDs and `firstDataIds`.
  - The fixture is BGC protein comparisons, all `reverseComplement:false`.
  - They run on dev push only (ignored by the PR-smoke and functional configs).
- **`gallery-session-regeneration.playwright.spec.js:96-100`**: 12 fields and no handles; BGC fixture, no reverse complement.
- **`losat-cache.test.mjs:436-444`** and **`run-analysis-derived-cache.test.mjs:89`**: validate only the shape of the `viewTransform.reverse` metadata.

What they miss: no test downloads a raw TSV for a reversed record and checks its coordinates.

## FE-12 (Specific table color 'none': CLI accepts, Web rejects)

Nearest tests:
- **`tests/web/file-imports.test.mjs:95-99`**: rejects `not-a-color`; every other row is hex.
- **`tests/web/color-utils.test.mjs:16-25`**: tests the `'none'` mode helpers only.
- **`tests/test_color_table_parsing.py:11-79`**: 4/5-column hex rows and missing values; no `none`.
- **`tests/web/auxiliary-file-history.playwright.spec.js:6`** and **`python-rule-parity.playwright.spec.js:117-119`**: browser `t_color` imports, hex only.

What they miss:
- The Node harness structurally cannot see this bug. `file-imports.js:54` has `domFreeNamedColor = !globalThis.document?.createElement && /^[a-z]+$/i`, so under Node `'none'` passes.
- In the browser, `resolveBrowserNamedColor` (`color-utils.js:109-124`) returns null for `none`, so the hex regex rejects it.
- No browser test imports a `none` row, and no test compares CLI and Web acceptance of the same TSV.

# Class 9: Python core defects

## CO-05 (outfmt 6 with 13-14 columns misread)

I confirmed in memory that pandas `read_csv(names=<12 columns>)` on a 13-field row makes column 1 the index and shifts every field.

Readers, with no shared owner (only the `COMPARISON_COLUMNS` constant is shared):
- `gbdraw/io/comparisons.py:73-78`
- `gbdraw/session_request_codec.py:3330-3335`
- `gbdraw/api/record_planning.py:1426-1431`
- `gbdraw/analysis/conservation.py:175-180`: the file path uses the same 12-name pattern. The `_coerce` handling at `:147-149` only helps DataFrame inputs; by the time it runs on a file, the shift has already happened.
- `gbdraw/analysis/protein_colinearity.py:3084-3103`: this one is strict and rejects anything other than 12 columns.

Tested column counts are 12 only:
- `tests/test_comparisons.py:41-46` (helper) and `:130-168`
- `tests/test_linear_multi_record_comparisons.py:455-460`
- `tests/test_circular_conservation.py:47-48` (`DataFrame(columns=COMPARISON_COLUMNS)`)
- `tests/fixtures/losat-outfmt6-numeric-contract.json`, used by `tests/test_session_io.py:328-340`
- Every fixture file is 12 columns: `examples/*.tblastx.out`, `tests/test_inputs/*.tblastx.out`, `gbdraw/web/tutorial-data/**/*.l*tsv`.

What they miss: there is no test with more than 12 columns anywhere.

## TR-01 (multi_record_canvas legend ignores custom track slots)

Nearest tests:
- **`tests/test_circular_track_slots.py:1798-1829` and `:1832-1870`**: the only slot legend-content assertions ("AT skew", "AT skew (+)/(-)"). Both use single-record `assemble_circular_diagram_from_record` with `legend="right"`.
  - The file has 0 uses of `multi_record_canvas` or `from_records`.
- **`tests/test_circular_multi_canvas.py:685-720`**: shared legend with no slots. It asserts only that there is one legend group.
- **`test_circular_multi_canvas.py:2290-2350`**: track table goes to a mocked single builder; forwarding only.
- Multi-record plus slot renders all use `legend="none"` and check geometry or IDs:
  - `test_circular_conservation.py:470-502`
  - `test_circular_svg_id_integrity.py:132-150`
  - `test_depth_track.py:1620-1700`
  - `test_svg_id_contract.py:1085-1098`
- `test_api_library_usage.py:611-655` is a mocked forwarding test.
- The linear analogue is tested: `test_depth_track.py:1546-1585` checks a linear multi-record slot `legend_label` shows up as `data-legend-key`. Circular multi-record has no equivalent.

What they miss: nothing asserts legend content for circular grid plus slots, which is the `diagram.py:3114-3128` path (the single-record path is `assemble.py:2766`).

## PV-08 (Web default settings fail with long /organism definitions)

How Web defaults are captured: static checks only, never rendered.
- `gbdraw/web/js/web-ux-profile.js:3-17` and `session-active-config-contract.js:16` (`track_type 'tuckin'`).
- Tested by `tests/web/mode-profiles.test.mjs:67-81, 353-370` (deepEqual on constants).
- `tests/test_documentation_reference_contracts.py:179-199` (string-contains).
- `tests/test_web_mode_profiles.py:15-60` (generated-file check plus source-text asserts).
- `tests/test_mode_profiles.py:95-128` with `tests/fixtures/mode_semantic_parity.json`, which covers thresholds, gc/skew and axis colors. It has no grid, separate-strands or track-type fields.

CLI tests using the Web default set:
- Every `--multi_record_canvas` test (`test_circular_multi_canvas.py:2104-2839`, 7 calls; `test_depth_track.py:2335`) uses `dummy.gb` with mocked loaders or builders. None adds `--separate_strands`.
- The session-cli-compatibility composite case runs `--multi_record_canvas` alone on `Homo sapiens` + `Escherichia phage Lambda`.

Organism names:
- Multi-canvas test organisms are at most 27 characters.
- The only gallery session with the full Web default set is `Vnig_TUMSAT-TG-2018` (grid, 6 records, tuckin, strandedness, gc/skew), with organism `Vibrio nigripulchritudo` (24 characters).
- The long-organism fixtures (about 43-52 characters: SARS-CoV-2, MERS-CoV, NC_007795, LvMJNV, BGC0000709) are all single-record.
- The 2-record GCF gbff files (`GCF_000354175.2`, 44 characters) appear only in Web gallery text checks.

## FE-07 (CDS without /translation does not use cds=True, so GTG/TTG is not translated as M)

Nearest tests:
- **`tests/test_web_feature_metadata.py:348-365`**: `ATGAAA` with table 11 gives `MK`.
- **`:368-385`**: length not divisible by 3 gives a warning.
- **`:303-345`**: the explicit `/translation` path.

What they miss:
- No GTG or TTG start codon anywhere in the tests. An ATG start gives the same answer with or without `cds=True`.
- The other translator, `gbdraw/analysis/protein_colinearity.py:2767`, is not compared against this one.

# Reference outputs, markers and CI

**Reference cases:** `tests/test_output_comparison.py:89-184` has 16 cases.
- Circular (10): all single-record (MjeNMV, and AP027078 for repeat_underlay). They cover separate strands, labels, radial labels, tuckin/middle/spreadout, and no-gc/no-skew combinations.
- Linear (6): basic, gc/skew, separate strands, a 2-genome case, one with 12-column tblastx and `--align_center`, and repeat_underlay.
- Every case uses `--legend none`.
- None uses `--multi_record_canvas`, track slots, comparisons with more than 12 columns, or custom definition or ruler fonts.
- `TestOutputComparison` has no marker, so it runs in core-pr. Generation is `reference_generation` and skipped by default (`conftest.py:68-82`).

**Markers** (`pyproject.toml` markers: slow, regression, circular, linear, reference_generation, recipe, gallery, browser):
- `test_run_info_exact_replay.py` is `browser`.
- `test_color_table_parsing.py` has one test marked regression+circular.
- The comparisons, track-slot, multi-canvas and linear-layout files use only circular/linear markers.
- `test_session_compat`, `test_genome_loading`, `test_record_metadata`, `test_web_feature_metadata`, `test_mode_profiles`, `test_web_mode_profiles` and `test_circular_feature_width` have no markers.

**CI impact** (`tools/ci-impact-policy.mjs`):
- `python-core` (`gbdraw/{api,config,configurators,core,io}/`, `cli|circular|linear.py`, `gbdraw/data/`) requires on PR: web-change-budget, core-pr, gallery, lint, web-contracts-pr, web-pr-smoke (lines 38, 183-186). `renderer` (`diagrams|tracks|legend|…`) requires the same set.
- `gbdraw/session*.py` maps to session-persistence, which adds recipes-standard.
- `gbdraw/analysis/*`, `gbdraw/web_support/*` and `gbdraw/cli_utils/session.py` match no rule, so they get full PR_JOBS.
- JS under `gbdraw/web/js/**` (except session/config/history owners) is `web-runtime`, which does not include core-pr.

What each job runs (`.github/workflows/test.yml`):
- **core-pr** (179-211): `pytest -m "not slow and not (recipe or gallery or browser)" -n auto`, then `git diff --exit-code tests/reference_outputs/`.
- **web-contracts-pr** (374-421): `node --test tests/web/*.test.mjs` (maxdepth 1, so `tests/web/contracts/` is excluded), then `pytest -m "browser and not slow"`.
- **web-pr-smoke** (423-468): only `@pr-smoke`-tagged Playwright tests (15 tags in 11 files). record-display-discovery, session-cli-compatibility and joint-display-placement are not among them.
- **playwright-functional** and **losat-cache-browser-acceptance**: push to dev or workflow_dispatch only.

# Systemic patterns

1. **Single-record, well-formed fixtures only.** No Prokka headers, no featureless FASTA records, no BLAST files wider than 12 columns, no reverse-complemented records in LOSAT downloads, no long organism names on multi-record canvases. Every slot legend label is alphanumeric.
2. **Parity checks compare the wrong pair.**
   - JS fast path vs Worker helper, where the helper shares the GFF assumption.
   - Source recipe vs Exact replay, both CLI.
   - JS slot spec vs JS parser.
   - None compares against the Python loader, Python parser or Web preview SVG, except one default-font circular fixture.
3. **Tests pass because test values equal the defaults.** Definition font 18 gives interval 20; the CLI legend default `right` equals the projection fallback; no scale-font override. So a dropped value is invisible.
4. **Web defaults are asserted as constants, never rendered.** The `--multi_record_canvas --separate_strands` tuckin gc/skew combination only reaches the CLI with mocked builders, `legend="none"`, or a short organism name.
5. **Legend content is rarely asserted on multi-record paths.** Grid and slot tests use `legend="none"`, and every reference output uses `--legend none`.
6. **The Node test environment hides browser-only behaviour** (the DOM-free named-color bypass in FE-12).
7. **Tier gaps.** The Playwright CLI-Web compatibility and discovery-parity specs run only on dev push. JS-only changes don't trigger core-pr, and Python-core PRs don't run the functional Playwright suite.
8. **Readers without a shared owner.** BLAST TSV is read in 4 places with `names=` plus one strict reader. Slot text has two separate builders (JS) and one parser (Python), with no shared vectors.
