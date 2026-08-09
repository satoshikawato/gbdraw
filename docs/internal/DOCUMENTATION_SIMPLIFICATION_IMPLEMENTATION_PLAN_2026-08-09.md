# Documentation simplification implementation plan

Status: complete
Created: 2026-08-09
Release target: 0.14.0
Current authority: this plan supersedes the public-information-architecture,
page-count, mandatory-surface, and branch-only URL-compatibility requirements
in the 2026-08-03 and 2026-08-04 documentation plans. Those files remain as
historical records.

## Objective

Ship 0.14.0 with a small public documentation structure whose pages have clear,
non-overlapping ownership. Keep the current reproducible Tutorials and their
execution evidence. Remove the public Task guides and Concepts & decisions
categories, merge their useful material into Technical documentation or FAQ,
and delete duplicated reader-facing procedures.

The target public routes are:

1. Tutorials
2. Technical documentation
3. FAQ
4. Gallery

Installation, tutorial-data acquisition, release notes, citation, and project
pages remain supporting documents. They do not form additional teaching
taxonomies.

## Verified baseline

- `docs/TUTORIALS/README.md` routes to 10 projects with GUI, CLI, and Python
  variants. These 30 pages are the strongest current learning material and stay
  public.
- `docs/HOW_TO/` contains 33 content pages. Several repeat the same records,
  settings, steps, and results as Tutorials. The useful remainder is mostly
  exact interface behavior, command syntax, input contracts, and
  troubleshooting.
- `docs/EXPLANATION/` contains six short pages. Choice questions fit FAQ;
  renderer defaults, schemas, compatibility details, privacy boundaries, and
  output contracts fit Technical documentation.
- `docs/scenarios/manifest.json` currently treats 81 evidence scenarios as 81
  public chapters, hard-codes role counts, and requires GUI, CLI, and Python
  variants for every Tutorial project.
- The capture and recipe runners provide valuable execution evidence. Their
  current one-to-one link to public pages is the defect; the evidence itself is
  not.
- `.agents/skills/document-generate/SKILL.md` prescribes Diataxis and explicitly
  favors completeness and multiple document types over minimal ownership.
- `.agents/skills/love-me-love-my-docs/` has strong reproducibility rules, but
  its chapter census and information-role rules can recreate the same
  proliferation.

## Editorial decisions

### Tutorials

- Rename the public `Start here` category to `Tutorials`.
- Keep the current 30 Tutorial pages and the outcome-by-surface index.
- Keep `Start here` only as a call to action for the first Circular and Linear
  projects, not as the category name.
- Do not require future projects to have one page for every surface. A separate
  surface page requires a materially different reader journey. Test or capture
  parity does not require public-page parity.
- Remove links whose only purpose is to send a Tutorial to a duplicate Task
  guide. Link to the exact Technical documentation page or FAQ answer instead.

### Technical documentation

- Keep `docs/REFERENCE/` as the storage path and rename its public label and H1
  to `Technical documentation`.
- Exact behavior, inputs, controls, options, schemas, API contracts, sessions,
  outputs, and short operational examples belong here.
- Treat `docs/CLI_Reference.md`, `docs/RECIPES.md`,
  `docs/SESSION_COMPATIBILITY.md`, and `docs/SVG_SEMANTIC_HOOKS.md` as technical
  subdocuments and expose them from the Technical documentation index.
- Extend an existing owner before adding a new public technical page.

### FAQ

- Keep answers short and question-shaped.
- Put layout choice, interface choice, method choice, privacy boundaries, and
  troubleshooting here.
- Link exact values and full contracts to Technical documentation instead of
  duplicating them in FAQ.

### Gallery

- Keep Gallery as outcome discovery and showcase material.
- Gallery entries may link to reproducible Tutorials. They do not become the
  canonical owner of procedures or technical contracts.

## Content disposition

The old category name does not survive as a public route. Each old page has one
of three outcomes: delete duplicated prose, merge unique material into the
owner below, or retain only its executable block as internal scenario evidence.

| Retired material | Public owner |
| --- | --- |
| GUI and CLI sequence inputs, GFF3/FASTA pairing, DDBJ/GenBank input boundaries, and TSV manifests | `REFERENCE/input-formats-and-tsv-schemas.md`, `REFERENCE/web-app.md`, `REFERENCE/command-line.md`, and compact examples in `RECIPES.md` |
| Circular multi-record layout; Linear rows, regions, orientation, labels, and rulers | `REFERENCE/web-app.md`, `REFERENCE/command-line.md`, layout answers in `FAQ.md`, and `RECIPES.md` |
| Uploaded BLAST, LOSATN, TLOSATX, LOSATP Pairwise/Similarity groups/Collinear, and Circular rings | `REFERENCE/comparison-programs-thresholds-and-results.md`, surface references, matching Tutorials, and `RECIPES.md` |
| GC, skew, depth, annotations, track slots, axes, and placement | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md`, surface references, matching Tutorials, and `RECIPES.md` |
| Colors, labels, visibility, shapes, strokes, overlap, titles, and legends | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md`, surface references, matching Tutorials, and `RECIPES.md` |
| Web inspection, editing, undo/reset, sessions, cache reuse, and export | `REFERENCE/web-app.md`, `REFERENCE/session-and-request-compatibility.md`, `REFERENCE/output-formats-and-export.md`, and the interactive-session Tutorial |
| Python drawing, comparisons, tracks, GFF3, in-memory records, bytes, typed requests, and session round trips | `REFERENCE/python-api.md`, `REFERENCE/typed-requests.md`, `REFERENCE/session-and-request-compatibility.md`, and matching Tutorials |
| Circular versus Linear and GUI versus CLI versus Python | concise decision answers in `FAQ.md`, with links to Tutorials and surface references |
| Comparison-method choice | concise answer in `FAQ.md`; the full matrix and limitations in `REFERENCE/comparison-programs-thresholds-and-results.md` |
| Track and outer-composition model | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md` |
| Publication and reproducibility checklist | `REFERENCE/output-formats-and-export.md` |
| Browser privacy, offline execution, cache behavior, and performance boundaries | `REFERENCE/web-app.md` plus a short FAQ answer |

### Capability owners

This map separates public prose ownership from execution evidence. The listed
path and anchor are the primary owner of the exact contract. Tutorials and FAQ
may teach or summarize the same capability, but they link here for exhaustive
facts.

| Public owner | Capabilities |
| --- | --- |
| `REFERENCE/web-app.md#execution-privacy-and-offline-use` | `entry.hosted-gui`, `entry.local-gui` |
| `REFERENCE/command-line.md` | `entry.circular-cli`, `entry.linear-cli` |
| `REFERENCE/python-api.md` | `entry.simple-python`, `program.simple-drawing`, `program.in-memory-records`, `program.bytes`, `program.comparisons-tracks` |
| `REFERENCE/typed-requests.md` | `entry.typed-python`, `program.typed-requests` |
| `REFERENCE/session-and-request-compatibility.md` | `program.session-round-trip`, `interactive.session-save-load`, `interactive.comparison-cache` |
| `REFERENCE/input-formats-and-tsv-schemas.md` | `input.genbank-embl-ddbj`, `input.gff3-fasta`, `input.multi-record`, `input.record-selection`, `input.regions`, `input.reverse-complement`, `input.tsv-manifests` |
| `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#diagram-layout` | `layout.circular-single`, `layout.circular-multi`, `layout.linear-single-row`, `layout.linear-multi-row`, `layout.track-placement`, `layout.rulers`, `layout.axes`, `layout.titles`, `layout.definitions`, `layout.legends` |
| `REFERENCE/comparison-programs-thresholds-and-results.md` | `comparison.uploaded-tables`, `comparison.browser-losatn`, `comparison.browser-tlosatx`, `comparison.losatp-pairwise`, `comparison.losatp-groups`, `comparison.losatp-collinear`, `comparison.selected-mixed-edges`, `comparison.circular-rings` |
| `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations` | `track.gc-content`, `track.dinucleotide-skew`, `track.at-skew`, `track.depth-numeric`, `track.region-annotations`, `track.custom-circular-slots`, `track.custom-linear-slots` |
| `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation` | `presentation.palettes`, `presentation.feature-colors`, `presentation.qualifier-priority`, `presentation.label-filters`, `presentation.label-overrides`, `presentation.visibility`, `presentation.shapes`, `presentation.strokes`, `presentation.overlaps` |
| `REFERENCE/web-app.md#preview-search-and-editor` | `interactive.preview-search`, `interactive.feature-popup`, `interactive.match-popup`, `interactive.group-popup`, `interactive.feature-editing`, `interactive.legend-editing`, `interactive.layout-editing`, `interactive.undo-redo-reset` |
| `REFERENCE/output-formats-and-export.md` | `output.static-svg`, `output.interactive-svg`, `output.png`, `output.pdf`, `output.eps`, `output.ps`, `output.comparison-data`, `output.fasta-data`, `output.reproducible-handoff` |

### Retired-page ledger

The primary owner receives the reusable procedure or decision. Exact option and
schema fragments follow the capability-owner map above. CLI and Python evidence
means the marked executable block and runner-checked TSV literals move to the
single internal evidence registry. GUI evidence remains in its manifest entry,
capture flow, assertions, and screenshots; it needs no replacement page.

| Retired public page | Primary public owner | Evidence after deletion |
| --- | --- | --- |
| `HOW_TO/README.md` and the three surface indexes | `DOCS.md` | none |
| `GUI/use-genbank-and-gff3-fasta-inputs.md` | `REFERENCE/input-formats-and-tsv-schemas.md` | GUI capture flow |
| `GUI/arrange-multiple-circular-records.md` | `REFERENCE/web-app.md#circular-multi-record-canvas` | GUI capture flow |
| `GUI/arrange-linear-records-regions-and-orientation.md` | `REFERENCE/web-app.md#record-selection-and-layout` | GUI capture flow |
| `GUI/use-uploaded-blast-results.md` | `REFERENCE/web-app.md#comparison-surfaces` | GUI capture flow |
| `GUI/use-tlosatx.md` | `REFERENCE/comparison-programs-thresholds-and-results.md` | GUI capture flow |
| `GUI/add-circular-similarity-rings.md` | `TUTORIALS/GUI/add-precomputed-circular-comparison-rings.md` | GUI capture flow |
| `GUI/create-protein-similarity-groups.md` | `TUTORIALS/GUI/compare-proteins-losatp.md` | GUI capture flow |
| `GUI/draw-collinear-protein-blocks.md` | `TUTORIALS/GUI/compare-proteins-losatp-collinear.md` | GUI capture flow |
| `GUI/add-depth-gc-and-skew-tracks.md` | `TUTORIALS/GUI/build-a-quantitative-genome-map.md` | GUI capture flow |
| `GUI/add-region-annotations-and-track-slots.md` | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations` | GUI capture flow |
| `GUI/style-features-labels-titles-and-legends.md` | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation` | GUI capture flow |
| `GUI/control-feature-visibility-shapes-strokes-and-overlaps.md` | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation` | GUI capture flow |
| `GUI/inspect-and-edit-a-diagram.md` | `REFERENCE/web-app.md#preview-search-and-editor` | GUI capture flow |
| `GUI/save-restore-undo-and-reproduce-work.md` | `REFERENCE/session-and-request-compatibility.md` | GUI capture flow |
| `GUI/export-publication-and-interactive-figures.md` | `REFERENCE/output-formats-and-export.md` | GUI capture flow |
| `CLI/use-genbank-ddbj-gff3-and-fasta.md` | `REFERENCE/input-formats-and-tsv-schemas.md` | internal executable registry |
| `CLI/use-input-tables.md` | `REFERENCE/input-formats-and-tsv-schemas.md#manifest-tables` | internal executable registry |
| `CLI/arrange-multiple-circular-records-and-tracks.md` | `REFERENCE/command-line.md#record-selection-and-layout` | internal executable registry |
| `CLI/arrange-linear-records-regions-orientation-labels-and-rulers.md` | `REFERENCE/command-line.md#record-selection-and-layout` | internal executable registry |
| `CLI/draw-precomputed-comparisons.md` | `REFERENCE/comparison-programs-thresholds-and-results.md` | internal executable registry |
| `CLI/run-losatp-pairwise.md` | `REFERENCE/comparison-programs-thresholds-and-results.md` | internal executable registry |
| `CLI/create-protein-similarity-groups.md` | `TUTORIALS/CLI/compare-proteins-losatp.md` | internal executable registry |
| `CLI/draw-collinear-protein-blocks.md` | `TUTORIALS/CLI/compare-proteins-losatp-collinear.md` | internal executable registry |
| `CLI/add-depth-gc-and-skew-tracks.md` | `TUTORIALS/CLI/build-a-quantitative-genome-map.md` | internal executable registry |
| `CLI/add-region-annotations-and-track-slots.md` | `REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields` | internal executable registry |
| `CLI/set-colors-labels-visibility-shapes-and-strokes.md` | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation` | internal executable registry |
| `CLI/save-and-regenerate-sessions.md` | `REFERENCE/session-and-request-compatibility.md` | internal executable registry |
| `CLI/export-static-and-interactive-outputs.md` | `REFERENCE/output-formats-and-export.md` | internal executable registry |
| `PYTHON/draw-circular-and-multi-record-diagrams.md` | `REFERENCE/python-api.md` | internal executable registry |
| `PYTHON/draw-linear-diagrams-and-comparisons.md` | `TUTORIALS/PYTHON/compare-genomes-losatn.md` | internal executable registry |
| `PYTHON/add-tracks-annotations-colors-and-labels.md` | `REFERENCE/python-api.md#layout-and-track-options` | internal executable registry |
| `PYTHON/use-gff3-in-memory-records-and-bytes.md` | `REFERENCE/python-api.md` | internal executable registry |
| `PYTHON/build-typed-requests-and-round-trip-sessions.md` | `REFERENCE/typed-requests.md` | internal executable registry |
| `EXPLANATION/README.md` | `DOCS.md` | none |
| `EXPLANATION/choose-circular-or-linear.md` | `FAQ.md#should-i-use-circular-or-linear` | none |
| `EXPLANATION/choose-gui-cli-or-python.md` | `FAQ.md#should-i-use-the-hosted-app-local-gui-cli-or-python` | none |
| `EXPLANATION/choose-a-genome-comparison-method.md` | `FAQ.md#which-comparison-method-should-i-use` | none |
| `EXPLANATION/understand-tracks-axes-and-layout.md` | `REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#tracks-axes-and-annotations` | none |
| `EXPLANATION/prepare-publication-and-reproducible-figures.md` | `FAQ.md#how-should-i-prepare-a-figure-for-publication` | none |
| `EXPLANATION/browser-privacy-offline-execution-and-performance.md` | `FAQ.md#is-my-data-uploaded-and-can-i-work-offline` | none |

## Evidence architecture

The public Task-guide pages currently double as executable scenario sources.
Preserve the smallest required evidence without presenting those files as user
documentation or moving a second manual under `docs/internal/`.

- Delete all 37 Markdown files under `docs/HOW_TO/` after their public material
  is merged.
- Create one `docs/internal/SCENARIO_EVIDENCE.md` containing the 13 CLI and five
  Python marked executable blocks plus the runner-checked TSV literals. It is a
  machine-owned evidence registry, not reader documentation.
- GUI H-scenarios keep their capture flow, manifest settings and assertions,
  and existing screenshot artifacts. They do not require a source page.
- Keep the `H-GUI-*`, `H-CLI-*`, and `H-PY-*` identifiers, runners, capture
  flows, and verified artifacts where they still test distinct product
  behavior. Stable internal IDs do not define the public information
  architecture.
- Change manifest schema version 2 from a `chapters` inventory to a `scenarios`
  registry.
- Remove per-scenario `canonical_for`. Scenario `capabilities`, assertions, and
  execution results describe evidence coverage; they do not own public prose.
- Add the capability-owner map above to `public_owners` in the manifest. Any
  number of scenarios may exercise a capability, but its exact public contract
  has one owner.
- Change retired How-to entries to the `evidence` role and remove their
  `destination`. CLI and Python evidence entries use `execution.source` to
  point to the internal registry.
- Remove the six Explanation entries. They provide no independent executable
  evidence after their content is merged.
- Remove `required_surfaces`, the no-partial-project rule, and any test that
  treats cross-surface parity as a prerequisite for a public Tutorial. Existing
  verified variants and their semantic parity checks remain.
- Remove exact scenario counts and mandatory surface lists from tests. Retain
  unique IDs, valid paths, fixture declarations, execution contracts,
  screenshot ownership, and current cross-surface output checks.

The deleted HOW_TO and EXPLANATION paths have never existed on `main` or any
release tag. They are branch-only, so this plan explicitly supersedes the old
one-release compatibility-stub requirement for those paths. Released flat
paths such as `QUICKSTART.md`, `GFF3_FASTA.md`, `EXPORT.md`, `PYTHON_API.md`,
`TYPED_API.md`, and `WORKFLOW_GUIDE.md` remain short compatibility routers.

## Governance changes

- Update `CLAUDE.md` and `CONTRIBUTING.md` to make the four-route structure and
  minimum-owner rule current policy.
- Require `reader question`, `existing owner`, and `keep / merge / delete / new`
  before a public page is added.
- State that coverage does not require one page per capability or interface.
- State that `avoid-ai-writing` checks prose style only. It is not an
  information-architecture, pedagogy, accuracy, or duplication acceptance
  gate.
- Remove the project-local generated `document-generate` skill. Its upstream
  template is not present in this repository, and retaining a skill that
  mandates the retired taxonomy would leave a conflicting active rule.
- Revise `love-me-love-my-docs` and its structure reference so evidence
  scenarios may support one public page or no public page, and so public-page
  planning defaults to merging with an existing owner.

## Implementation phases

### Phase 0: Plan and baseline

- [x] Inspect the worktree and preserve unrelated changes.
- [x] Audit current navigation, page inventory, history, rules, manifest, and
  tests.
- [x] Review this plan independently and apply the blocking corrections to
  evidence sources, public ownership, surface policy, URL compatibility, and
  verification scope.

Exit condition: the target routes, content dispositions, evidence boundary,
and verification commands are explicit before another file is changed.

### Phase 1: Change governance first

- [x] Update `CLAUDE.md` and `CONTRIBUTING.md`.
- [x] Remove `.agents/skills/document-generate/SKILL.md`.
- [x] Revise `love-me-love-my-docs/SKILL.md` and
  `references/manual-structure.md`.
- [x] Validate the revised skill with `quick_validate.py`.

Exit condition: no active repository rule requires Diataxis, one page per
capability, or one page per interface.

### Phase 2: Consolidate public content

- [x] Rewrite `README.md`, `docs/DOCS.md`, the Tutorials indexes, and the
  Technical documentation index for the four public routes.
- [x] Merge useful HOW_TO and EXPLANATION material according to the disposition
  table.
- [x] Extract the 18 CLI/Python executable blocks and checked TSV literals into
  `docs/internal/SCENARIO_EVIDENCE.md`, then delete `docs/HOW_TO/`.
- [x] Delete `docs/EXPLANATION/` after all useful content has an owner.
- [x] Replace all public links and breadcrumbs to the retired categories.
- [x] Confirm that public prose contains no source-checkout-only command unless
  the page is explicitly for contributors.

Exit condition: `docs/HOW_TO/` and `docs/EXPLANATION/` no longer exist, and the
top-level teaching navigation uses Tutorials, Technical documentation, FAQ, and
Gallery. Supporting Installation, data, release, citation, and project links
remain directly reachable.

### Phase 3: Separate evidence from public ownership

- [x] Migrate `docs/scenarios/manifest.json` to schema 2.
- [x] Update capture and recipe helpers to read `scenarios`; CLI/Python recipe
  extraction reads `execution.source`, while GUI evidence has no destination.
- [x] Replace page-count, Diataxis-label, mandatory-surface, and mandatory
  Task-guide-link tests with ownership and routing contracts.
- [x] Keep focused tests for literal recipes, capture flows, artifacts, and
  semantic output parity.

Exit condition: adding or retaining executable evidence cannot force a new
public page, and adding a public page requires a distinct reader question and
one canonical owner.

### Phase 4: Verify the release documentation

- [x] Run targeted scenario, Tutorial, reference, link, recipe, and capture
  contract tests.
- [x] Run `pytest tests/ -v -k documentation`.
- [x] Run the affected CLI recipe, Python recipe, GUI capture, fixture,
  reproduction, and onboarding test modules.
- [x] Run the full fast suite, `pytest tests/ -v -m "not slow"`, because recipe,
  capture, fixture, reproduction, and onboarding tests consume the manifest
  outside the documentation test selection.
- [x] Run `test_public_markdown_local_targets_exist` and a stale-token/path scan
  for `HOW_TO/`, `EXPLANATION/`, `Task guides`, `Concepts & decisions`,
  `required_surfaces`, `manifest["chapters"]`, and `canonical_for` in active
  rules, code, tests, and public docs.
- [x] Run a detect-only `avoid-ai-writing` audit on changed public prose, make
  targeted corrections, and audit again.
- [x] Review production, tests, public docs, internal evidence, and generated
  artifacts as separate diffs.
- [x] Record commands, results, and any justified residual risk below.

Exit condition: all required commands pass, public navigation has no retired
routes, and the evidence table contains observed results rather than intended
results.

## Acceptance contracts

1. The README and documentation home expose Tutorials, Technical
   documentation, FAQ, and Gallery. `Start here` is a Tutorial call to action,
   not a category.
2. The current 30 Tutorial pages remain indexed and keep their verified inputs,
   commands or UI actions, expected results, and artifacts.
3. No public Markdown page links to `docs/HOW_TO/` or `docs/EXPLANATION/`.
4. No active repository rule tells contributors to choose a Diataxis quadrant.
5. Every capability in the scenario manifest has exactly one existing public
   owner path and anchor, separately from any scenarios that exercise it.
6. Internal evidence scenarios may share a public owner and do not count as
   public pages.
7. Future Tutorial projects may implement only the surfaces that serve real
   reader journeys.
8. Exact technical facts are stated once; FAQ and Tutorials link to that owner.
9. Changed public prose passes the repository style audit without altering UI
   labels, CLI options, code, paths, accessions, or scientific qualifications.
10. Documentation and focused scenario tests pass from the current worktree.

## Evidence log

| Phase | Command or review | Result |
| --- | --- | --- |
| Baseline | `git status --short` plus scoped path review | The worktree already contained unrelated changes. Documentation-scope paths were inspected before editing, and unrelated changes were preserved. |
| Baseline documentation contracts | `env TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_documentation_scenario_contracts.py tests/test_tutorial_documentation_contracts.py tests/test_documentation_reference_contracts.py tests/test_documentation_contracts.py -q` | 42 passed |
| Plan review | Independent plan critic plus evidence-consumer audit | Blocking issues corrected before implementation |
| Skill validation | `/home/kawato/micromamba/bin/python /home/kawato/.codex/skills/.system/skill-creator/scripts/quick_validate.py .agents/skills/love-me-love-my-docs` | `Skill is valid!` |
| Evidence extraction | Compared all 18 retired executable blocks with `docs/internal/SCENARIO_EVIDENCE.md`; checked H-CLI-02 and H-CLI-11 TSV literals through the runner assertion | One former owner per block, no content mismatch; TSV evidence matched |
| Targeted documentation tests | `env TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_documentation_scenario_contracts.py tests/test_tutorial_documentation_contracts.py tests/test_documentation_reference_contracts.py tests/test_documentation_contracts.py tests/test_documented_recipes.py tests/test_python_howto_recipe_contracts.py tests/test_cli_how_to_recipe_contracts.py tests/test_cli_comparison_how_to_recipe_contracts.py tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py -q -k 'not regenerate'` | 68 passed, 14 deselected |
| Static GUI capture contracts | `env TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_documentation_capture_contracts.py tests/test_gui_interactive_capture_contracts.py tests/test_gui_nucleotide_comparison_capture_contracts.py tests/test_gui_presentation_capture_contracts.py tests/test_gui_tracks_capture_contracts.py tests/test_gui_protein_comparison_capture_contracts.py -q` | 67 passed |
| Documentation test selection | `env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin /home/kawato/micromamba/bin/python -m pytest tests/ -v -k documentation` | 72 passed, 2,854 deselected |
| Affected consumers | CLI and Python recipe, GUI capture, fixture, reproduction, and onboarding modules | 142 of 143 grouped items passed; onboarding alone then passed in 10.84 seconds. The grouped failure was its internal 180-second subprocess timeout, so the timeout was not weakened. |
| Public link audit | `env TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_reproduce_examples.py::test_public_markdown_local_targets_exist -q` | 1 passed |
| Retired-route audit | Contract tests plus exact searches of active rules, code, tests, and public Markdown | No public `HOW_TO/`, `EXPLANATION/`, retired category labels, schema-1 `chapters`, `required_surfaces`, or `canonical_for` dependency remained. Internal `docs/capture/flows/how_to/` names remain implementation-only stable paths. |
| Manifest and runner audit | Parsed schema 2 and compared manifest-derived IDs with the cross-surface and standalone runners | 75 scenarios: 30 Tutorial, 33 evidence, 10 reference, 2 auxiliary; 72 capabilities have one of 12 public owners; 23 CLI and 15 Python verified runner IDs matched. |
| Final prose audit | Detect-only `avoid-ai-writing` scan, targeted correction pass, and added-line re-audit | No added AI markers, stock transitions, public governance narration, or em-dash prose remained. |
| Retired-content audit | Independent comparison of every retired HOW_TO and EXPLANATION page against its public owner and current implementation | Eight initially omitted contracts were restored to five existing owners; final re-audit found no remaining unique current contract loss or duplicate owner. |
| H-CLI-13 rerun | Exact failing node, its 11-test module, 10 hash-seed checks, and six fresh PDFs including four parallel renders | All reruns passed. Raw PDFs differed only in volatile metadata; all normalized payloads matched the tracked PDF at SHA-256 `5a53ed36085a379a2696e49cbb4e1a219b156eb5123718bc5140bce36bcf8c69`. No artifact regeneration was justified. |
| Final focused confirmation | Scenario, Tutorial, Technical documentation, general documentation, manifest-runner, and public-link contracts | 44 passed in 34.93 seconds from the final worktree. |
| Python lint | `ruff check` on the changed capture, recipe, helper, manifest-tool, and documentation-test Python files | All checks passed. |
| Diff review | Separate scoped reviews of governance, public docs, internal evidence, tests, and runner changes | Task-scoped LF diffs passed `git diff --check`. The whole dirty worktree is not a valid whitespace gate because unrelated pre-existing CRLF files are reported wholesale; `--ignore-space-at-eol` review isolated the intended edits in the two touched mixed-EOL files. The internal evidence registry retains one trailing-tab TSV row because those delimiters encode six empty fields. No Tutorial screenshot, reference SVG, or generated public artifact was changed. |
| Full fast suite | `env TMPDIR=/tmp PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin /home/kawato/micromamba/bin/python -m pytest tests/ -v -m "not slow"` | 2,903 passed, 17 skipped, 6 deselected, 4 existing Biopython warnings in 1,353.58 seconds. The clean rerun followed termination of a runaway read-only audit process that had consumed one CPU core during two earlier attempts. |

## Stop conditions

- Do not delete a unique behavior or runnable example until its public owner or
  internal evidence destination is identified.
- Do not regenerate Tutorial screenshots or reference SVGs merely because their
  links move.
- If an extracted executable block cannot be read and run from
  `docs/internal/SCENARIO_EVIDENCE.md`, fix the evidence boundary before
  continuing.
- If a public page needs a new category to remain discoverable, revise the
  four-route landing page instead of recreating a taxonomy layer.
