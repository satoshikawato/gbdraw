# Web Gallery first-Generate parity remediation implementation plan

- Date: 2026-08-17
- Status: proposed; implementation not started
- Planning baseline: `dev` at `6f084dca72671868b74c16d8a58a4af3e0388027`
- Session writer / canonical request baseline: session version 40 / request schema 5
- Scope owner: Interactive Gallery session authoring, browser session replay, Gallery artifacts, and publication gates

> Next-session instruction: execute this document with
> `$execute-plan-with-evidence`. Update phase status only after recording the
> commands and observed evidence required by that phase. Re-audit the baseline
> before editing because Gallery artifacts and Web session contracts may have
> changed after the planning commit.

Revision note: this version makes the baseline audit executable, separates
characterization from the behavior-neutral validator extraction, defines the
scientific-baseline and semantic-SVG contracts, and splits static refresh
validation from authoritative real-browser replay. All implementation phases
remain pending.

## 1. Objective and completion boundary

Make every bundled Interactive Gallery session an executable reproduction of
the figure initially shown to the user.

For every entry in `gbdraw/web/gallery/examples.json`, the following chain must
hold:

```text
Gallery source SVG
        ==
saved Session Result
        ~=
browser preview immediately after import
        ~=
first Generate Diagram result with no UI edits
```

`==` is the existing source/result artifact equality. `~=` is semantic SVG and
scientific-content equivalence, not raw byte equality.

Completion requires all of the following:

1. All 11 bundled sessions import successfully through the real browser path.
2. The first unchanged `Generate Diagram` succeeds for all 11 sessions.
3. The loaded and generated diagrams have equivalent record topology, tracks,
   annotations, features, labels, legends, comparison evidence, and geometry.
4. Every resource needed by the committed render and the active Web draft is
   embedded and resolvable.
5. Gallery refresh rejects a session that violates the shared current-writer
   structural contract, canonical request projection, complete active-resource
   closure, or publication fixed-point contract. The post-refresh real-browser
   replay separately proves compatibility with the complete Web importer; do
   not claim that the static refresh validator is equivalent to the importer.
6. The same real-browser replay workflow runs and succeeds for pull requests to
   both integration branches and blocks deployment. Repository rules must also
   mark its jobs required before claiming that it blocks pull-request merges.
7. General user sessions retain the supported ability to save a committed
   Result together with a different, not-yet-generated draft.
8. Regenerated sessions, SVGs, examples, thumbnails, tutorial media, and
   `examples.json` remain mutually consistent and are visually reviewed.

This plan does not authorize publishing, pushing, deploying, tagging, or
changing repository settings. Obtain separate authorization for those actions.

## 2. Problem statement

The current session has two valid authorities:

- `renderRequest` and `results` own the last committed render and initial
  preview.
- `config` owns the active Web controls and the next Generate.

That split is intentional for ordinary user sessions. A user can save a figure,
change controls without generating, and preserve both the committed Result and
the pending draft. `restoreCurrentWriterActiveConfig()` in
`gbdraw/web/js/services/config.js` therefore correctly gives current-writer
`config` authority over the next Generate.

The Gallery publication pipeline currently violates a different requirement:
published examples must not contain a pending draft. In
`tools/refresh_gallery_sessions.py`, `_refresh_one_session()` replays the saved
canonical request through the CLI and `_merge_refreshed_gallery_artifacts()`
updates the committed request/result/artifacts while retaining the old active
`config`. The refresh tool only contains partial repairs for the legend and
Circular track slots. It can consequently publish a correct-looking preview
whose next Generate uses different settings.

The current validation stack does not see this failure:

- Python session validation checks the session envelope and typed request, but
  not the Web current-writer `config.form` / `config.adv` allowlist.
- `_validate_staged_gallery_session()` checks resource references from
  `renderRequest`, but not all active `webFiles` bindings.
- existing Gallery semantic tests compare saved request/result/source assets;
  they do not traverse browser import, lazy resource restoration, active-state
  projection, LOSAT preflight, Worker staging, and first Generate;
- deployment refreshes Gallery outputs and publishes them without replaying the
  refreshed sessions in Chromium.

## 3. Reproduced baseline evidence

These observations were made against the planning baseline. Re-run them in
Phase 0 and retain the new evidence rather than assuming the bytes are
unchanged.

| Example | Observed failure | Root cause or boundary |
| --- | --- | --- |
| `HmmtDNA_ATskew` | Import and Generate succeed, but the generated SVG differs. The loaded viewBox was approximately `1280.93 x 953.70`; the generated viewBox was approximately `1315.41 x 958.85`. The AT-skew track remains present, so this is not simply a missing slot. | The committed full CLI-derived configuration and the Web-authored next request do not resolve to identical presentation geometry. The request diff includes representation-only changes and at least one effective Circular definition/layout difference. |
| `majanivirus_orthogroup` | Import and Generate succeed, but pairwise-match paths fell from 627 to 262 in the browser audit. The loaded and generated SVGs had 2,346 versus 1,975 paths. | The active draft and the saved comparison/color artifacts are not one reproducible authority. The committed color table has three rules, while `config.rules` is empty. The committed table also contains Web-parser-incompatible named colors `yellow` and `red`. |
| `vibrio-harveyi-group-collinear` | The loaded preview contains comparison groups, ribbons, and a comparison legend. The next request contains zero comparisons and the generated diagram removes them all. | `config.linearComparisonPlan.mode` is `none`, while the committed request represents adjacent LOSATP collinear comparison. The dedicated test currently expects this divergence and is not connected to normal CI. |
| `WSSV_genome_comparison` | Import succeeds; Generate fails with `Input file is not available for browser FASTA extraction.` | All 20 exact FASTAs are already embedded and bound through `webFiles.conservationLosatFastaSources`. `extractLosatFastaFast()` and `extractAllLosatFastaFast()` reject frozen lazy session-resource views because they do not expose `.text()`, even though shared `readFileText()` materializes them correctly. |
| `tobacco-chloroplast` | Import fails with `Current session active configuration contains unknown config.adv field(s): cli_circular_track_order, cli_circular_track_slots.` | The v40 generator-owned artifact retains two retired fields. Removing them in memory makes import and Generate succeed and restores the four region annotations. |

Additional artifact facts at the baseline:

- `HmmtDNA_ATskew` stores the intended five Circular slots in both the request
  and draft: `features`, `gc_content`, `gc_skew`, `a_skew_2`, `ticks`.
- `majanivirus_orthogroup` stores eight precomputed comparison descriptors, an
  `orthogroupResult`, and a generated-pipeline descriptor; the active plan is
  `adjacent + losat`. Its committed specific-color TSV contains
  `#89d1fa`, `yellow`, and `red`, while the tutorial already documents the last
  two as `#ffff00` and `#ff0000`.
- `vibrio-harveyi-group-collinear` retains 59 raw LOSAT cache entries. Its
  committed comparison settings are collinear, adjacent scope,
  orientation-and-identity color, minimum anchors 3, maximum unit gap 2,
  maximum diagonal drift 2, candidate limit 5, and 16 threads.
- `WSSV_genome_comparison` contains 45 resources, including 20 FASTAs and 20
  corresponding comparison results.
- `tobacco-chloroplast` otherwise retains `plastome_regions` with LSC, IRb,
  SSC, and IRa plus three intended slots: features, annotations, and GC content.

## 4. Binding design decisions

### D1. The Gallery adds a publication invariant; it does not redefine sessions

Do not change the general current-writer rule in
`restoreCurrentWriterActiveConfig()`. Do not overwrite every imported current
session's draft with its committed request. Preserve the divergent-draft cases
in `tests/web/contracts/session-regenerate-intent.playwright.spec.js`.

Only generator-owned Gallery sessions must be published with no render-affecting
pending draft.

### D2. Gallery active state must be produced through the Web projection boundary

Use `projectCanonicalSessionRequest()`, `buildCanonicalRenderRequest()`,
`migrateLegacyLinearComparisonDraft()`, and the existing comparison resolver.
Do not create a Python copy of Web `config` rules, a second argv-shaped Web
contract, or a hand-maintained list of every active field.

Add one explicitly named Gallery-publication normalization entry point in
`gbdraw/web/js/services/gallery-session-migration.js`, invoked by a small Node
tool from the refresh pipeline. Its contract is:

1. promote a supported Gallery session to the current format;
2. project the committed request and resources without allowing stale stored
   render-affecting config to overwrite the projection;
3. create an explicit Web comparison draft for a Gallery example, including
   CLI-authored examples, through the existing comparison-plan machinery;
4. rebuild the current Web request/resources from that state;
5. retain non-rendering document metadata and any applied editor/semantic
   override already represented in the committed Result, feature state, or
   interactive metadata, while dropping genuinely pending render changes;
6. retain required resource payloads;
7. validate the resulting current-writer config; and
8. be idempotent.

In particular, audit `config.webEdits` and `features` separately from ordinary
form state. BGC currently stores the applied orthogroup-name override
`og_18 -> livA`. Publication normalization must preserve and reapply that
committed interactive meaning; it must not erase it as if it were an ungenerated
draft.

Do not add example-specific mutations to the browser importer. If an example
needs a publication recipe that cannot be inferred from the canonical request,
place declarative metadata beside the existing `GallerySessionExample`
definition and make it required/tested. Do not scatter example-ID branches
through runtime code.

### D3. Replace partial synchronization with one owner

After the Gallery publication normalizer proves the existing behaviors, remove
or reduce the superseded partial owners in `tools/refresh_gallery_sessions.py`:

- `_sync_legacy_legend_control_with_render_request()`;
- `_sync_circular_track_draft_with_render_request()`; and
- any palette repair whose only purpose is covered by canonical Web projection.

Retain a helper only if a focused test proves it owns a distinct input that the
canonical projection cannot represent. Do not keep both a general normalizer
and several silent per-field fixups.

### D4. One reader owns File-like compatibility

`readFileText()` and `readFileBytes()` in
`gbdraw/web/js/services/file-content-cache.js` own native File and lazy session
resource compatibility. Callers must not pre-reject a resource based on
`.text()` or `.arrayBuffer()` capabilities.

For WSSV, remove the two `.text` capability guards in
`extractLosatFastaFast()` and `extractAllLosatFastaFast()`. Do not add fake
methods to frozen session-resource views and do not add a WSSV-only bypass.

### D5. Current-writer validation stays strict

Do not add `cli_circular_track_order` or `cli_circular_track_slots` to the Web
allowlist and do not add a v40 reader migration for generator-owned branch-only
fields. Rewrite the chloroplast artifact through the Gallery generator and
make the generator call the actual JS current-writer validator.

The validator is not currently a browser-neutral boundary:
`validateCurrentWriterActiveConfig()` lives in `config.js`, derives field names
from live Vue state, and shares an import graph with Gallery migration. Do not
solve this with a second allowlist or permanent global browser shims in the
Node tool. Extract the structural contract into a pure module, for example
`gbdraw/web/js/services/current-writer-active-config.js`:

- move the current-writer domain/shape validation and Circular/Linear slot
  structural validation into that module;
- obtain the allowed `form` and `adv` field inventories from browser-neutral
  default factories or constants;
- if necessary, move `createDefaultForm()` and `createDefaultAdv()` to a pure
  `state-defaults.js` module and re-export them from `state.js` to preserve
  existing imports;
- have `config.js` and the Node publication tool import the same pure validator;
- assert that runtime `Object.getOwnPropertyNames(state.form/state.adv)`,
  including accessor-backed `legend` and `plot_title_position`, exactly match
  the pure contract inventory.

This extraction must be behavior-neutral for ordinary session import/save.
Do not let the new module import `config.js`, Vue state, DOM code, or Gallery
migration.

No session-version or request-schema bump is expected for this remediation.

### D6. Gallery comparison evidence is scientific output

Do not make Majanivirus or Vibrio tests pass by merely accepting fewer or zero
matches. First establish the intended current settings and review the resulting
match identities.

The default decision is to preserve the displayed pre-fix scientific content:

- Majanivirus: preserve the approved nine-record, eight-comparison-group
  similarity-group figure and its match-ID set.
- Vibrio: preserve non-empty adjacent collinear groups, ribbons, and comparison
  legend across all 11 replicons.

If the current algorithm cannot reproduce the existing Majanivirus 627-match
figure, stop at the Phase 4 decision gate. Compare the 627 and 262 results,
their settings, orthogroup identities, and the declared CLI recipe. Obtain
repository-owner approval before changing the public scientific baseline.

Approval must produce a durable, reviewable test input rather than an
unrecorded decision. Add
`tests/web/fixtures/gallery-scientific-baselines.json` as a manually reviewed,
non-generated contract. For every comparison-bearing Gallery example it must
record:

- a schema version and example ID;
- comparison mode, scope, and the settings that determine scientific output;
- sorted comparison, pairwise-match, orthogroup, and collinear-block IDs where
  applicable, plus their counts and set digests;
- the source session/SVG hashes and evidence commit used for the decision; and
- an approval state and durable repository review reference.

Keep the sorted ID sets in reviewable text, even when digests are also stored,
so a failure can report additions and removals. Do not derive this file from the
newly regenerated Gallery artifacts in the same acceptance run. A change to an
approved scientific set must be an explicit source/test-input diff reviewed
before generated sessions or SVGs are refreshed. Until the Majanivirus and
Vibrio entries have an approved state and review reference, their Phase 4 gate
is blocked and the plan cannot be marked complete.

### D7. Semantic SVG equality is the required oracle

Raw SVG bytes can differ because of metadata, serialization order, resource
representation, or harmless numeric spelling. A low-resolution raster hash can
also hide missing semantics. Require both:

1. semantic SVG equivalence; and
2. example-specific scientific identity/invariant checks.

The semantic comparison must include at least:

- viewBox and output topology;
- record IDs, order, grouping, and transforms;
- track slot IDs, renderer, side, order, and geometry;
- feature IDs and visible geometry/style;
- annotation IDs, labels, and styles;
- visible label and legend text plus positions;
- comparison group, pairwise-match, orthogroup, and collinear-block IDs/counts;
- applied feature/orthogroup labels and interactive metadata overrides;
- fill, stroke, visibility, font size, and relevant path/transform attributes.

Normalize color spelling and insignificant numeric serialization only. Do not
allow `exact || canonical || semantic || raster` as the pass condition; semantic
equivalence and scientific invariants are mandatory. Raster output remains a
diagnostic attachment and manual-review aid.

Implement the oracle as a versioned `SemanticSvgSnapshot` contract owned by
`tests/web/helpers/gallery-reproduction.cjs`, not as an ad hoc list of hashes.
Version 1 must define all of the following before it is used as a required gate:

- match elements by stable `id` / `data-gbdraw-*` identity; reject missing or
  duplicate required identities, and use an explicit structural path only for
  elements that are documented as identity-free;
- preserve child order where it affects paint order, grouping, record order, or
  track order; do not globally sort SVG elements;
- parse path and point data into commands and numeric operands, and normalize
  transforms to a documented matrix representation;
- compare resolved visible style, including inline attributes, `style`, CSS,
  and inheritance, for the enumerated drawing properties;
- normalize numeric spelling, units, exponent form, and negative zero without
  silently accepting a changed value. Start with numeric equality after
  parsing; add an attribute-specific tolerance only after cross-platform
  evidence demonstrates the need and records a bound that cannot hide visible
  geometry drift;
- normalize SVG whitespace only where SVG semantics make it insignificant;
- enumerate excluded metadata/editor-only fields rather than ignoring an
  entire subtree; and
- emit the first differing semantic owner/path plus compact added/removed ID
  sets instead of only whole-document digests.

Add oracle contract tests that make one mutation at a time to a representative
SVG: remove/reorder a track, annotation, feature, label, legend, comparison, or
record; change geometry, transform, text, color, visibility, or font size; and
introduce a duplicate semantic ID. Every mutation must fail. Separate positive
fixtures must prove that color spelling and numeric serialization differences
allowed by the contract pass. The required pass condition is:

```text
SemanticSvgSnapshot equality
AND approved example-specific scientific invariants
```

Exact, canonical, and raster comparisons remain diagnostics only.

### D8. The Gallery manifest owns complete test inventory

All 11 current examples and every future example must participate
automatically. Add a required `regeneration_tier` to
`GallerySessionExample` and emit `regenerationTier` in `examples.json`.
Reject a missing or unknown tier. Do not infer test coverage from a mutable file
size threshold and do not maintain an unrelated test allowlist.

Initial tiers:

- `standard`: HmmtDNA basic, Lambda basic, HmmtDNA AT skew, tobacco
  chloroplast, and BGC;
- `heavy`: V. nigripulchritudo, both Hepatoplasmataceae examples,
  Majanivirus, and WSSV;
- `vibrio-special`: Vibrio Harveyi group.

### D9. Use real assets and real browser boundaries

Browser tests must upload the actual `.json` or `.json.gz` path with
`setInputFiles()`. Do not expand large gzip payloads in the Node test process,
mock away the diagram Worker, replace the biological fixture, skip LOSAT cache
hydration, or inject a synthetic `.text()` method.

Run one example per fresh browser context and `workers: 1` within a job. Release
the page/context/Worker after every example. Use CI matrix parallelism between
examples, not concurrent heavy sessions in one browser process.

### D10. Generated Gallery artifacts remain generator-owned

Never hand-edit session JSON, source SVG, displayed SVG, thumbnail WebP, or
`examples.json`. Regenerate them from declared inputs and the publication
normalizer. Review production code, tests, source inputs, sessions, SVGs,
thumbnails, tutorials, and workflow diffs separately.

Do not modify `examples/gbdraw_social_preview.png`.

### D11. Deployment must test the bytes it publishes

The minimal change keeps deploy-time Gallery refresh, but runs the same static
and real-browser corpus gates after refresh and before copying to `public/`.
The long-term cleaner model is to refresh and review Gallery artifacts in a pull
request and deploy those unchanged tracked bytes. That deployment redesign is
not required for this remediation as long as post-refresh browser replay blocks
upload.

### D12. The post-merge session must be a request/config fixed point

Normalizing before CLI replay is insufficient because the CLI-produced session
can replace `renderRequest` during `_merge_refreshed_gallery_artifacts()`. After
the merge, run publication normalization again in memory and require a fixed
point:

1. the stored committed request and the request rebuilt from the stored active
   draft have equal semantic request signatures;
2. a second publication-normalization pass does not change config, request,
   resource bindings, or referenced resource payload hashes; and
3. the browser first-Generate semantic SVG remains equivalent.

The request signature must cover mode/grouping, record source identity and
presentation, effective diagram options, tracks, annotations, comparison
descriptors/settings, layout, and referenced resource content hashes. Exclude
only enumerated, evidence-backed non-drawing fields such as an output filename
or equivalent resource ID/representation. Do not ignore a whole
`diagramOptions`, `comparisons`, `config`, or interactive-metadata subtree.

### D13. Static refresh validation and browser replay are separate gates

Do not reproduce the complete browser importer in Node or Python. The static
refresh gate owns only contracts that already have browser-neutral owners:

1. session and request schema validation;
2. the shared current-writer active-config structure;
3. canonical request projection;
4. complete committed and active resource closure;
5. post-merge request/config/resource fixed-point equality; and
6. publication-normalizer idempotence.

The real-browser gate remains authoritative for session upload/import, gzip and
lazy-resource restoration, Result admission and SVG sanitization, feature
catalog restoration, LOSAT preflight, Worker staging, first Generate, and
loaded/generated semantic equivalence. A static pass must never be described as
proof that the complete Web importer accepts the session.

Local refresh may replace generator-owned working-tree artifacts after its
transactional static gates pass, but no deployment copy or upload may occur
until the real-browser gate has tested those exact post-refresh bytes. If a
future requirement demands that the refresh command itself guarantee full
browser acceptance, add an explicit `--verify-browser` staging workflow that
serves and tests staged artifacts before replacement; do not approximate the
importer with more static allowlists.

## 5. Scope and non-goals

### In scope

- WSSV lazy session-resource FASTA extraction;
- Gallery-only active-draft normalization;
- current-writer and active-resource validation in refresh;
- a browser-neutral shared current-writer validation contract;
- AT skew, Majanivirus, Vibrio, WSSV, and chloroplast artifact repair;
- verification of the other six Gallery entries;
- generated sessions, source/display SVGs, thumbnails, `examples.json`, and
  affected tutorial result media;
- Node, Python, and Playwright regressions;
- a repeatable known-bad browser audit with retained per-example evidence;
- a versioned semantic-SVG oracle and independently reviewed scientific
  baseline fixture;
- required pull-request and deployment gates;
- explicit synchronization of declared Gallery download inputs used by the
  repaired Majanivirus recipe.

### Out of scope unless Phase 4 proves it necessary

- changing general user-session draft semantics;
- a new session or canonical request version;
- changing LOSAT algorithms or match definitions;
- changing Python public APIs or CLI option behavior;
- core renderer geometry changes;
- broad Gallery UI redesign;
- unrelated tutorial screenshot renovation;
- manual edits to reference outputs.

The Circular CLI/codec currently derives
`objects.definition.circular.interval` from the Circular definition font size.
If AT-skew parity reveals a genuinely missing Web rendering option rather than
a projection bug, pause before expanding scope. A new public Web option must be
implemented across state, control, request projection, persistence, typed
decode, CLI/Python parity where applicable, documentation, and tests under the
`change-gbdraw-rendering-surface` contract. Do not hide it in Gallery-only JSON.

## 6. Surface and ownership matrix

| Surface | Current owner | Planned responsibility |
| --- | --- | --- |
| Lazy resource backing | `gbdraw/web/js/services/session-resource-backing.js` | Remain unchanged unless a failing test proves a backing bug. Frozen lazy views remain method-free. |
| File content compatibility | `gbdraw/web/js/services/file-content-cache.js` | Remain the single native/lazy File reader. |
| Generate/LOSAT orchestration | `gbdraw/web/js/app/run-analysis.js` | Remove invalid capability guards; preserve shared error/fallback semantics. |
| Pure current-writer contract | new `gbdraw/web/js/services/current-writer-active-config.js` and, if needed, new `gbdraw/web/js/state-defaults.js` | Own browser-neutral field inventories and structural validation used by browser config and Node publication tooling. |
| Active config import/save | `gbdraw/web/js/services/config.js` | Keep v40 draft authority, import the pure validator, and preserve its public re-export if callers rely on it. |
| Canonical Web projection | `gbdraw/web/js/services/session-request.js` | Remain the only UI-state/request projection boundary. |
| Linear comparison intent | `gbdraw/web/js/app/linear-comparisons.js` | Derive Gallery publication plans through existing normalization/resolution. |
| Gallery promotion | `gbdraw/web/js/services/gallery-session-migration.js` | Add one current publication normalizer; do not change ordinary session import semantics. |
| Node promotion bridge | `tools/promote_gallery_session.mjs` or a focused companion | Normalize/validate current Gallery sessions, handle JSON/gzip deterministically, and report the example name on failure. |
| Gallery refresh | `tools/refresh_gallery_sessions.py` | Invoke publication normalization for every Gallery entry, validate complete resource closure, retain transaction rollback, then regenerate artifacts. |
| Gallery manifest/assets | `tools/prepare_interactive_gallery_assets.py` | Own required replay tier and all generated siblings. |
| Declared Gallery downloads | `tools/prepare_interactive_gallery_assets.py` plus source-to-target metadata | Copy the repaired Majanivirus specific-color TSV from one declared source and validate byte equality. |
| WSSV provenance | `tools/restore_wssv_gallery_fastas.py` | Keep the existing 20-source identity/hash contract. |
| Known-bad browser audit | new `playwright.gallery-audit.config.js` and `tests/web/contracts/gallery-first-generate-audit.serial.spec.js` | Reproduce the five baseline symptoms with real uploads and retain structured evidence without becoming a normal passing regression contract. |
| Python session validation | `gbdraw/session_io.py` | Keep typed/session validation. Do not duplicate Web field allowlists unless a shared generated contract replaces both owners. |
| Browser equivalence helper | new `tests/web/helpers/gallery-reproduction.cjs` | Own the versioned `SemanticSvgSnapshot`, reusable capture, structured diff, scientific-ID extraction, and attachments. |
| Scientific baselines | new `tests/web/fixtures/gallery-scientific-baselines.json` | Own manually reviewed settings and ID sets independently of generator-owned Gallery outputs; never regenerate it as an artifact side effect. |
| Standard browser cases | `tests/web/gallery-session-regeneration.playwright.spec.js` | Cover all `standard` manifest entries. |
| Heavy browser cases | new `tests/web/contracts/gallery-heavy-regeneration.serial.spec.js` | Cover each `heavy` manifest entry using the same helper. |
| Vibrio browser case | `tests/web/contracts/vibrio-full-generation.serial.spec.js` | Preserve existing memory/Worker/determinism coverage while reversing the divergent-draft expectation. |
| Static publication contract | new `tests/web/gallery-publication-contract.test.mjs` | Validate all manifest sessions, post-merge request/draft fixed points, and resource closure without rendering or retaining all expanded sessions at once. |
| CI/deploy | `.github/workflows/test.yml`, `.github/workflows/deploy_web.yml`, `package.json`, Playwright configs | Run the same publication contract before merge and before upload. |

Likely unchanged boundaries:

- `gbdraw/circular.py::_render_output_request()`;
- `gbdraw/session_request_codec.py` canonical encode/decode;
- public Python API and CLI parser;
- `tests/reference_outputs/`;
- Gallery HTML/JS/CSS, unless the manifest field requires a small reader change;
- WSSV FASTA payloads and pinned hashes.

## 7. Execution phases

### Phase 0: re-establish the baseline

Status: pending

#### 0.1 Record the untouched repository baseline

1. Read `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`,
   `.agents/skills/change-gbdraw-rendering-surface/SKILL.md`, and
   `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md`.
2. Fetch and create the implementation branch from the latest `origin/dev` as
   required by repository guidance. Do not commit on `main` or `dev`.
3. Record branch, HEAD, working-tree state, session file hashes, and unrelated
   pre-existing changes.
4. Confirm Node and Python Playwright availability and install the locked Node
   test dependencies when needed.
5. Prepare the browser wheel without refreshing the cache-bust token.

#### 0.2 Add and run an audit-only known-bad browser harness

Before any production/runtime edit, add
`playwright.gallery-audit.config.js` and
`tests/web/contracts/gallery-first-generate-audit.serial.spec.js`. The spec is
an audit/characterization driver, not the desired-behavior regression test. It
must use the same capture primitives that will move into
`tests/web/helpers/gallery-reproduction.cjs` in Phase 1 and must:

- enumerate exactly the five reported examples from a constant beside the
  audit, not from an unrelated file-size rule;
- upload each actual `.json` or `.json.gz` session with `setInputFiles()` in a
  fresh browser context;
- capture import outcome, page/console errors, loaded state, the active-next
  request, first-Generate outcome, and generated state without UI edits;
- write `loaded.svg`, `generated.svg`, `summary.json`, and `errors.json` under
  `/tmp/gbdraw-gallery-first-generate-baseline/<example-id>/`, plus PNG/trace
  diagnostics when available;
- write `/tmp/gbdraw-gallery-first-generate-baseline/index.json` only after all
  five cases settle, listing the commit, example/session hashes, outcome, and
  relative evidence paths;
- include HEAD, session SHA-256, request/schema versions, resource counts and
  hashes, viewBox/path/comparison counts, and the observed symptom in each
  summary; and
- pass in `known-bad` mode only when every example reaches its expected
  baseline symptom. Setup failures, unexpected import/Generate errors, missing
  evidence, or a changed symptom must fail the audit.

Do not keep an `expect known-bad` test in normal CI after the repair. Phase 1
must convert the reusable capture into desired-behavior failing tests; delete
the audit-only expectations when their retained baseline evidence has been
recorded.

#### 0.3 Record static context

Run the existing all-11 static checks and record session sizes and hashes. These
checks are context, not substitutes for the browser audit.

Commands:

```bash
git status --short --untracked-files=all
git rev-parse --abbrev-ref HEAD
git rev-parse HEAD
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
npm ci
node -e "console.log(require.resolve('@playwright/test'))"
python tools/prepare_browser_wheel.py
GBDRAW_GALLERY_AUDIT_MODE=known-bad \
GBDRAW_GALLERY_AUDIT_DIR=/tmp/gbdraw-gallery-first-generate-baseline \
  npx playwright test --config=playwright.gallery-audit.config.js
git rev-parse HEAD > /tmp/gbdraw-gallery-first-generate-baseline/head.txt
sha256sum gbdraw/web/gallery/sessions/*.gbdraw-session.json* \
  > /tmp/gbdraw-gallery-first-generate-baseline/session-sha256.txt
wc -c gbdraw/web/gallery/sessions/*.gbdraw-session.json* \
  > /tmp/gbdraw-gallery-first-generate-baseline/session-sizes.txt
node --test tests/web/*.test.mjs
python -m pytest \
  tests/test_gallery_session_semantics.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_wssv_gallery_fastas.py -q
```

Expected known-bad evidence before repair:

- AT skew semantic/geometry mismatch;
- Majanivirus comparison mismatch;
- Vibrio comparison disappearance;
- WSSV FASTA extraction error;
- chloroplast current-writer import rejection.

If any expected failure has disappeared, inspect the intervening diff and
update this plan's evidence section before continuing. Do not apply a stale
fix.

Phase gate:

- one machine-readable summary names the exact commit and all five session file
  hashes, and every example directory contains the required evidence;
- the five symptoms are reproduced through real upload/import/first Generate,
  or a changed symptom is investigated and this plan is updated before repair;
- the audit-only test/config/helper diff is reviewed separately from production
  code, and no production/runtime or generator-owned artifact has changed;
- existing all-11 static checks and sizes are recorded without being presented
  as browser-import evidence.

### Phase 1: freeze desired behavior and establish shared test contracts

Status: pending

#### 1.1 Convert characterization into desired-behavior browser tests

Move the reusable parts of the Phase 0 audit, `capturePageEvidence()`,
`svgEquivalence()`, and related diagnostics from
`tests/web/contracts/session-regenerate-intent.playwright.spec.js` into
`tests/web/helpers/gallery-reproduction.cjs`.

The helper must:

- serialize the clean committed SVG;
- implement and version the D7 `SemanticSvgSnapshot` contract;
- capture normalized semantic structure and explicit scientific ID sets;
- capture active comparison intent and a canonical-request summary;
- compare loaded versus generated output;
- produce a compact structural diff, not only two hashes;
- attach loaded/generated SVG, PNG, diff JSON, errors, and trace on failure.

Add focused helper tests for every required D7 positive normalization and
negative mutation in `tests/web/gallery-reproduction.test.mjs` before relying
on the oracle for Gallery acceptance. A mutation that removes or changes one
required semantic category must make the helper test fail for that category.

Keep the existing divergent-user-draft tests using the helper or equivalent
coverage. Strengthen the existing no-draft fixture test so it asserts the
source-fixture-to-first-Generate semantic result rather than merely attaching
it, but do not change tests that intentionally prove a user draft can diverge.

Replace the audit-only known-bad expectations with assertions for the desired
import/first-Generate parity and run them against the untouched production
baseline. Record the expected failures for AT skew, Majanivirus, Vibrio, WSSV,
and chloroplast before the validator extraction or any other production edit.
After the retained audit evidence is referenced from the execution record,
delete the audit-only spec and config or convert the config directly into the
final Gallery reproduction config; do not leave a second dormant browser
harness beside the Phase 6 owner.

#### 1.2 Extract the current-writer validator behavior-neutrally

Only after the real-browser desired-behavior failures are retained, extract the
pure current-writer contract described in D5. `config.js` must continue to
export the same validator name if existing tests/imports use it, but its
implementation must delegate to the pure module.

Before and after the extraction, run the existing current-session and
session-draft-authority tests and compare their results. Add focused parity
tests proving the pure validator and the browser `config.js` wrapper
accept/reject the same fixtures, plus an import-cycle/module-import smoke test
for both `config.js` and the future Node publication entry point. This step is a
green, behavior-neutral refactor; do not add Gallery normalization or weaken an
accepted/rejected fixture here.

#### 1.3 Add a manifest-driven static contract

Create `tests/web/gallery-publication-contract.test.mjs` and make it enumerate
the generated manifest. Process one session at a time and release references
before the next, especially for Vibrio.

The pure current-writer contract from 1.2 must already exist before this test
imports Gallery migration. Assert that the runtime form/advanced field
inventory equals the pure inventory.

For each session assert:

- manifest, physical session, generator definition, and tier agree;
- every comparison-bearing example has a schema-valid independent scientific
  baseline entry whose final pass requires `approved` state, sorted IDs,
  settings/source hashes, and a non-empty repository review reference, and
  generated Gallery tooling cannot rewrite that fixture;
- current session version and request schema are correct;
- `validateCurrentWriterActiveConfig()` accepts the stored config;
- `projectCanonicalSessionRequest()` succeeds;
- projected active inputs have backing resources;
- all `renderRequest` and `webFiles` resource references resolve to non-empty
  resources;
- retired current-writer fields are absent;
- a Linear Gallery figure with committed comparisons does not publish an
  inactive `none` plan unless the committed figure genuinely has no comparison;
- the stored post-merge committed request and active-next request have equal
  semantic request signatures under a small enumerated non-drawing exclusion
  list; and
- publication normalization is a request/config/resource fixed point.

At this phase, missing or `pending` scientific approvals and the known stale
publication state are expected contract failures, not harness/setup failures.
Phase 4 owns review and approval of those entries; Phase 1 must not bless them
to make the static suite green.

Use the actual JS validator. Do not transcribe its allowlists into the test.
Make request-signature failures report the first changed owner/path and resource
hash rather than only a whole-object digest.

#### 1.4 Add standard and heavy browser corpus definitions

Expand `tests/web/gallery-session-regeneration.playwright.spec.js` over every
`standard` manifest item. Add
`tests/web/contracts/gallery-heavy-regeneration.serial.spec.js` for `heavy`
items. Keep Vibrio in its purpose-built serial spec.

Every example must execute:

1. upload actual session path;
2. assert import `{status: 'ok'}` and no page/console error;
3. assert a valid loaded preview and expected Result count;
4. capture loaded evidence;
5. invoke the actual `window.__GBDRAW_APP__.runAnalysis()`;
6. assert `{status: 'ok'}`, `errorLog === null`, and a committed Result;
7. capture generated evidence;
8. require semantic SVG and example-invariant equality.

The five desired-behavior cases must already have failed in 1.1 before the
behavior-neutral extraction. Expanding to the complete manifest may reveal
additional failures; record them without weakening the oracle or changing the
scientific baseline.

Phase gate:

- tests fail for the expected reasons, not setup or timeout errors;
- the D7 oracle contract tests detect every required semantic mutation and
  accept only the enumerated representation-only differences;
- the validator extraction is behavior-neutral under before/after evidence and
  has one owner (it must not be repeated in Phase 3);
- a deliberately removed track, annotation, or comparison is detected;
- an ordinary intentionally divergent user session remains supported.

### Phase 2: fix WSSV lazy-resource extraction

Status: pending

Production change:

- remove the `.text()` capability checks from
  `extractLosatFastaFast()` and `extractAllLosatFastaFast()` in
  `gbdraw/web/js/app/run-analysis.js`;
- allow `readFileText()` to produce the existing TypeError when neither native
  nor lazy resource reading is available;
- preserve the current fast parser, Worker fallback, cancellation, cache, and
  error-reporting behavior.

Focused tests:

1. Construct a real frozen session-resource view with no `.text()` method.
2. Exercise both single-record and all-record FASTA extraction through the
   Generate orchestration path.
3. Update `tests/web/wssv-gallery-fastas.test.mjs` so it no longer repairs the
   fixture by adding a fake `text()` method.
4. Run the real WSSV import-to-Generate browser case.
5. Preserve all 20 exact FASTA hashes, labels, colors, order, and BLAST result
   bindings.

Likely files:

```text
gbdraw/web/js/app/run-analysis.js
tests/web/run-analysis-simple-path.test.mjs
tests/web/wssv-gallery-fastas.test.mjs
tests/web/session-resource-backing.test.mjs          # only if coverage belongs here
tests/web/contracts/gallery-heavy-regeneration.serial.spec.js
```

Phase gate:

- the previous WSSV error is absent;
- Generate succeeds without network access or a LOSAT rerun caused by missing
  data;
- 20 conservation rings and their IDs/labels/colors remain present;
- native File and invalid-object error tests still pass.

### Phase 3: make Gallery refresh produce a valid no-draft session

Status: pending

#### 3.1 Add Gallery publication normalization

Add a distinct exported publication function to
`gbdraw/web/js/services/gallery-session-migration.js`. Do not overload ordinary
import or `restoreCurrentWriterActiveConfig()`. The browser-neutral validator
from Phase 1.2 is a dependency of this phase and must not be extracted or
redefined again here.

Implementation requirements:

- accept current and supported historical Gallery sessions;
- ignore stale render-affecting stored config when projecting the committed
  render into a new publication draft;
- preserve non-rendering provenance such as the supported CLI invocation;
- preserve committed feature/editor overrides and `webEdits` only when their
  application is proven by the saved Result/catalog/interactive metadata;
- use `forceWebDraft: true` only in the Gallery-publication path so a
  CLI-authored comparison figure gains explicit reproducible Web intent;
- derive comparison edges/settings with existing comparison functions;
- separate comparison-resource retention from comparison-descriptor authority;
  do not let `preserveComparisonResources()` overwrite the newly resolved
  descriptor list merely because old descriptors exist;
- retain only resources referenced by the committed/publication request,
  active `webFiles`, caches required for first Generate, and documented
  downloads;
- preserve exact WSSV FASTA provenance and deterministic gzip (`mtime=0`);
- validate and return a current, idempotent session.

Update `tools/promote_gallery_session.mjs` or add one focused companion CLI so
Python can invoke this function for a staged file. The tool must support JSON
and gzip, process one session at a time, and include the session path in errors.

#### 3.2 Integrate it transactionally

Change `tools/refresh_gallery_sessions.py` so every Gallery session passes
through publication normalization before CLI rendering. Test-input sessions
outside the Gallery inventory retain their existing semantics unless they are
explicitly declared publication artifacts.

Keep `_gallery_file_transaction()` and the stage-all-before-replace behavior.
If any normalization, render, validation, or asset step fails, restore every
target and publish no partial Gallery.

#### 3.3 Validate the Web contract and full resource closure

This subsection implements only the D13 static refresh gate. After
normalization and again after artifact merge:

- call the shared JS current-writer validator;
- project active files through the real session-request code;
- check resources referenced from `renderRequest`, `webFiles.bindings`,
  `conservationLosatFastaSources`, `conservationSequenceSources`, and other
  current active bindings;
- require present, non-empty payloads with valid encoding/declared size;
- preserve existing Python geometry, catalog, cache-schema, protein identity,
  and Vibrio size-limit validations.

Own persisted active-resource traversal in one Gallery publication helper.
Classify every property emitted by the current writer's representative
`webFiles` output as either a resource reference/container or an explicitly
non-resource field, and fail a contract test when a newly emitted field is
unclassified. Do not scatter separate lists of WSSV, comparison, annotation,
and input bindings across Python and Node. This helper may remain
Gallery-specific; it does not make the static validator equivalent to browser
import.

Then compute the D12 fixed point for every staged Gallery session. Rebuild the
active-next request from the merged session, compare its semantic request
signature with the stored committed request, and run publication normalization
a second time in memory. Fail before replacement if either signature or the
second normalized config/request/resource binding changes. This check applies
to all 11 examples, not only AT skew.

Do not broaden `collectCanonicalResourceIds()` to Web-only fields unless that
shared abstraction has at least two real consumers. A Gallery validation helper
is preferable to changing canonical-request ownership.

The real-browser corpus in Phases 4, 6, and 7 must still upload and import the
post-refresh files. Static validator success is not evidence for Result
admission, SVG sanitization, lazy materialization, Worker staging, or complete
Web importer acceptance.

#### 3.4 Remove superseded partial repairs

Once the full normalizer passes equivalent tests, remove the per-field legend,
track-slot, and palette synchronization paths it replaces. Add an idempotence
test proving a second refresh does not change semantic config/request/resources
apart from allowed timestamps or compression representation.

Likely files:

```text
gbdraw/web/js/services/gallery-session-migration.js
tools/promote_gallery_session.mjs                    # or one focused companion
tools/refresh_gallery_sessions.py
tests/web/gallery-session-migration.test.mjs
tests/web/gallery-publication-contract.test.mjs
tests/test_refresh_gallery_sessions.py
```

The validator/default/state files belong to Phase 1.2 and are dependencies,
not expected Phase 3 change points. Reopen them here only if a failing import or
contract test proves the supposedly behavior-neutral boundary is incomplete.

Phase gate:

- chloroplast's two retired fields are rejected in a negative fixture and absent
  from generated output;
- every staged Gallery config passes the actual JS validator;
- removing a WSSV active FASTA resource fails refresh before files are replaced;
- a second normalizer/refresh pass is semantically idempotent;
- ordinary imported current sessions still preserve intentional drafts.

### Phase 4: converge each reported example

Status: pending

Do not run a blind all-artifact refresh until each decision below has passing
focused evidence.

#### 4.0 Produce and approve scientific baseline evidence

For every comparison-bearing example, extract the old committed and candidate
first-Generate ID sets with the D7 helper. Write a review bundle under `/tmp`
containing sorted IDs, added/removed IDs, settings, counts, digests, source
session/SVG hashes, and representative loaded/generated SVG and PNG evidence.

Populate `tests/web/fixtures/gallery-scientific-baselines.json` only from the
reviewed decision. Majanivirus and Vibrio require explicit repository-owner
approval because their current candidates differ from the public figure. Other
comparison examples may retain their unchanged current sets after review, but
they still require explicit fixture entries. Do not generate or rewrite this
fixture from `refresh_gallery_sessions.py`.

#### 4.1 Tobacco chloroplast

- regenerate through the publication normalizer so
  `cli_circular_track_order` and `cli_circular_track_slots` disappear;
- preserve three track slots and the `plastome_regions` annotation set;
- require IDs `lsc`, `irb`, `ssc`, and `ira`, their labels, coordinate targets,
  and bracket styling;
- prove real import and first Generate success;
- do not weaken the current-writer reader.

#### 4.2 HmmtDNA AT skew

- compare the saved request, active projected request, and generated typed
  request after normalization;
- classify every diff as representation-only, existing Web projection bug, or
  genuinely unsupported CLI presentation state;
- preserve the five-slot order, `nt=AT`, colors `#deaf6e` / `#7294e3`, legend
  label `AT skew`, split feature lanes, tick layout
  `label_in_tick_out`, ajisai palette, qualifier priority, and representative
  labels including `ND1`;
- require semantic geometry equality with the initially displayed figure.

If an unsupported CLI-only value is the only obstacle, stop for the scope
decision described in Section 5. Do not silently refresh the public figure to a
different geometry merely to make the gate green.

#### 4.3 Majanivirus

- replace Web-incompatible named colors in the declared source TSV with exact
  hex equivalents: `yellow` -> `#ffff00`, `red` -> `#ff0000`;
- designate `tests/test_inputs/majani_custom_color_table.tsv` as the source of
  truth for this Gallery download;
- add declarative source-to-target metadata to
  `GallerySessionExample` (or one adjacent Gallery-download manifest) and make
  `prepare_gallery_assets()` copy it to
  `gbdraw/web/gallery/files/majani_custom_color_table.tsv`;
- include that generated download target in `_gallery_mutation_targets()` so a
  later refresh failure restores it transactionally;
- add a byte-equality test between the declared source and public download;
- regenerate the downloadable Gallery copy through
  `python tools/refresh_gallery_sessions.py`, rather than editing both copies;
- project the three specific rules and neutral default CDS palette into the
  active draft;
- preserve nine record labels, eight comparison groups, 152 orthogroups,
  legend captions, and the IDs selected through the D6 baseline decision;
- compare the 627-match saved figure and the current 262-match regenerated
  result using raw cache entries, filters, membership mode, and orthogroup IDs;
- emit the 627-versus-262 review bundle from 4.0 and obtain the decision required
  by D6 if current code cannot reproduce the existing figure;
- record the approved sorted match/orthogroup ID sets, settings signature,
  source hashes, and repository review reference in the independent scientific
  baseline fixture before refreshing its Gallery artifacts.

Likely source/generated files include:

```text
tests/test_inputs/majani_custom_color_table.tsv
gbdraw/web/gallery/files/majani_custom_color_table.tsv
gbdraw/web/gallery/sessions/majanivirus_orthogroup.gbdraw-session.json.gz
tools/prepare_interactive_gallery_assets.py
tools/refresh_gallery_sessions.py
tests/test_refresh_gallery_sessions.py
```

Do not change unrelated color tables unless the Gallery source inventory proves
they feed this artifact.

#### 4.4 Vibrio Harveyi group

- replace the Gallery-only `none` draft with explicit `adjacent + losat`
  comparison intent derived through the existing planner;
- preserve collinear mode, adjacent scope, orientation-and-identity colors,
  candidate/anchor/gap/drift/thread settings, 11 replicons, five record rows,
  record gap 48, record labels/subtitles, and cached raw evidence;
- rebuild derived comparison state from the retained raw cache as normal
  Generate does;
- change `vibrio-full-generation.serial.spec.js` from expecting zero
  comparisons and divergent SVGs to requiring loaded/first/repeated Generate
  equivalence;
- retain its current Worker reuse, memory, timeout, lifecycle, repeated-run,
  and recoverable-error assertions;
- compare the retained public and candidate collinear-block/pairwise-match ID
  sets in the 4.0 review bundle, then record the approved sets, settings
  signature, source hashes, and repository review reference in the independent
  scientific baseline fixture before refreshing its Gallery artifacts.

Update the screenshot register statement that the editable Vibrio session is
CLI-only with `linearComparisonPlan.mode = none`. Tutorial capture actions must
no longer manufacture comparison intent that the downloadable session lacks.

#### 4.5 Verify the remaining seven cases

The publication normalizer affects every example. Run the same contract for
HmmtDNA basic, Lambda, BGC, V. nigripulchritudo, both Hepatoplasmataceae
figures, and WSSV even when their initial symptom was not reported.

Phase gate:

- every focused browser case passes;
- AT retains approved geometry;
- every comparison-bearing example has a complete, reviewable scientific
  baseline entry independent of generated Gallery outputs;
- Majanivirus and Vibrio use explicitly approved scientific baselines with a
  durable repository review reference;
- no unreviewed comparison-count reduction is accepted;
- WSSV and chloroplast retain their exact biological annotations/resources.

### Phase 5: regenerate and visually review the complete Gallery artifact set

Status: pending

After Phase 4 decisions are fixed:

1. Run the transactional all-session refresh.
2. Regenerate every generator-owned sibling from declared inputs.
3. Review session, source SVG, displayed SVG, thumbnail, and manifest diffs as
   separate groups.
4. Confirm the refresh did not modify
   `tests/web/fixtures/gallery-scientific-baselines.json`; an intended
   scientific-baseline change must already have been reviewed in Phase 4.
5. Render old and new affected figures at the same readable scale.
6. Check clipping, label collisions, legend placement, quantitative tracks,
   annotations, comparison ribbons, and interactive metadata.
7. Recapture only tutorial result/control media made stale by the corrected
   session state.
8. Update `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md` before
   replacing screenshot files.
9. Run strict capture checks for every changed tutorial.

Primary artifact families:

```text
gbdraw/web/gallery/sessions/<id>.gbdraw-session.json[.gz]
gbdraw/web/gallery/sources/<id>.svg
gbdraw/web/gallery/examples/<id>.svg
gbdraw/web/gallery/thumbnails/<id>.webp
gbdraw/web/gallery/examples.json
```

Potentially affected tutorial/result media:

```text
gbdraw/web/gallery/media/HmmtDNA_ATskew/manual-09-01-atskew-preview.webp
gbdraw/web/gallery/media/majanivirus_orthogroup/manual-07-01-orthogroup-preview.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-04-00-run-adjacent-losat.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-04-01-search-method-collinear.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-04-02-color-mode-orientation-identity.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-04-03-adjacent-pairs.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-04-04-advanced-collinear.webp
gbdraw/web/gallery/media/vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp
gbdraw/web/gallery/media/WSSV_genome_comparison/manual-09-01-conservation-rings.webp
gbdraw/web/gallery/media/tobacco-chloroplast/manual-08-01-chloroplast-preview.webp
```

This is an audit list, not an instruction to replace every file. Keep a current,
truthful capture unchanged when session repair does not change its documented
state.

Commands:

```bash
python tools/refresh_gallery_sessions.py
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

Phase gate:

- generated inventory has no stale or orphan sibling;
- changed images pass same-size old/new visual review on desktop and mobile;
- tutorial controls agree with the corrected session, especially Vibrio's
  comparison plan;
- `examples/gbdraw_social_preview.png` is untouched.

### Phase 6: make replay a required CI and deployment gate

Status: pending

#### 6.1 Playwright configuration and scripts

Add `playwright.gallery-regeneration.config.js` with:

- Chromium;
- `fullyParallel: false`;
- `workers: 1`;
- up to 1,200,000 ms for the largest case;
- at most one CI retry for non-Vibrio cases;
- `trace: 'retain-on-failure'`;
- the existing local HTTP server/base URL;
- standard, heavy, and Vibrio reproduction specs.

Add one exact single-example selector, for example
`GBDRAW_GALLERY_EXAMPLE_ID`. Resolve it once in
`playwright.gallery-regeneration.config.js`, before Playwright discovers test
files:

- no selector means run the full set allowed by the selected config/spec;
- a selector must name exactly one manifest entry, otherwise configuration
  fails;
- use that entry's `regenerationTier` to set `testMatch` to only the owning
  standard, heavy, or Vibrio spec;
- the selected spec must then define exactly one test for that exact ID;
- a valid selected ID is not an error merely because non-owning specs exist;
  those specs are not loaded by the selected config;
- the static publication test proves that the union of tier IDs equals the
  manifest inventory.

This selector, rather than a broad `--grep @gallery-heavy`, is the CI matrix
boundary that guarantees one heavy example per job.

Add scripts such as:

```json
{
  "test:web:gallery-regeneration": "playwright test --config=playwright.gallery-regeneration.config.js",
  "test:web:gallery-regeneration:heavy": "playwright test --config=playwright.gallery-regeneration.config.js --grep @gallery-heavy",
  "test:web:gallery-regeneration:one": "playwright test --config=playwright.gallery-regeneration.config.js"
}
```

Keep `test:web:vibrio-generate` as a supported focused command, even if it
delegates to the new config/tag.

#### 6.2 Pull-request CI

In `.github/workflows/test.yml`:

- add `dev` to `on.pull_request.branches`; the workflow currently runs for
  pull requests to `main` only even though agent work integrates through
  `dev`;
- let the existing Node job run the static publication contract;
- keep standard five-entry replay in the ordinary functional Playwright job;
- add a required heavy matrix for V. nigripulchritudo, both
  Hepatoplasmataceae examples, Majanivirus, WSSV, and Vibrio;
- use `max-parallel: 2`, one example per job, and a per-job timeout of at least
  45 minutes;
- prepare the browser wheel in each browser job;
- upload Playwright `test-results/`, traces, SVGs, PNGs, and diff JSON on
  failure.

Set `fail-fast: false` for the diagnostic matrix. Add a stable, non-matrix
`gallery-regeneration-required` job with `if: ${{ always() }}` and `needs` on
the static publication job, standard replay job, and the complete
heavy/special matrix job. Its only responsibility is to fail unless every
dependency result is `success`. This stable aggregate status is the branch-rule
contract; individual matrix children remain diagnostic and may change as the
manifest grows.

Using the current workflow job IDs, the dependency shape should be equivalent
to the following (rename only if the underlying job IDs are deliberately
renamed in the same change):

```yaml
gallery-regeneration-required:
  name: Gallery regeneration required
  if: ${{ always() }}
  needs:
    - gallery
    - browser
    - playwright-functional
    - gallery-regeneration-heavy
  runs-on: ubuntu-latest
  steps:
    - name: Require every Gallery dependency
      run: |
        test "${{ needs.gallery.result }}" = "success"
        test "${{ needs.browser.result }}" = "success"
        test "${{ needs.playwright-functional.result }}" = "success"
        test "${{ needs.gallery-regeneration-heavy.result }}" = "success"
```

Here `gallery-regeneration-heavy` is the new matrix job ID and includes all
five `heavy` entries plus the `vibrio-special` entry. GitHub exposes the matrix
job's aggregate result to `needs`, so any failed or cancelled child makes the
stable aggregator fail. Retain the existing Python `gallery` and Node/browser
static owners in the dependency list instead of silently narrowing the gate to
Playwright alone.

The matrix must be derived from or statically checked against the required
manifest tiers so a new example cannot be silently untested.

Each heavy matrix job must set exactly one ID:

```yaml
env:
  GBDRAW_GALLERY_EXAMPLE_ID: ${{ matrix.example }}
```

and run `npm run test:web:gallery-regeneration:one`. Do not invoke the whole
heavy corpus once per matrix row.

#### 6.3 Deployment

In `.github/workflows/deploy_web.yml`, order the build as follows:

```text
install Python and Node dependencies
install Chromium
refresh Gallery sessions/assets transactionally
run static Node and Python Gallery contracts
prepare the cache-busted browser wheel
run all Gallery browser regeneration cases
build the release wheel
assemble public/
upload Pages artifact
deploy
```

No public copy/upload step may run after a failed replay. Test the refreshed
working-tree artifacts, not the checkout that existed before refresh. Record a
post-refresh session/manifest hash inventory before browser replay and confirm
the same bytes are copied into `public/` afterward.

Keep the D13 claims distinct in workflow names and logs: static contracts prove
schema/config/projection/resource/fixed-point validity, while the Playwright
jobs prove complete browser import and first Generate. Neither job may report
the other's scope as its own evidence.

#### External Gate A: required-check repository rules

Status: pending explicit repository-owner authorization

Workflow files can make jobs run and fail, but cannot by themselves make those
jobs required for merge. After the stable `gallery-regeneration-required` job
name is in place, inspect the branch protection/rulesets for both `dev` and
`main`. With explicit authorization, require that one aggregate status on the
branches through which changes merge. Do not require a mutable list of matrix
child names.

If this external change is not authorized, report the exact remaining gap and
do not claim that pull-request merges are blocked. Deployment remains blocked
inside `deploy_web.yml` because replay is part of the build job before upload.

Phase gate:

- standard and every heavy matrix entry run in PR CI;
- pull requests to both `dev` and `main` trigger the workflow;
- the exact selector loads only its owning spec and defines exactly one test;
- the former standalone Vibrio case participates in the aggregate check;
- `gallery-regeneration-required` fails for any failed, cancelled, or skipped
  dependency and becomes merge-required when External Gate A is authorized;
- deployment tests the post-refresh bytes before upload;
- failure artifacts identify the example and the first differing semantic
  category;
- repository rules mark the jobs required, or the handoff explicitly records
  External Gate A as pending authorization rather than completed.

### Phase 7: final verification and cleanup

Status: pending

Run focused checks first, then the complete justified gates.

Focused static/unit checks:

```bash
node --test \
  tests/web/gallery-reproduction.test.mjs \
  tests/web/file-content-cache.test.mjs \
  tests/web/session-resource-backing.test.mjs \
  tests/web/run-analysis-simple-path.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs \
  tests/web/gallery-session-migration.test.mjs \
  tests/web/gallery-publication-contract.test.mjs

python -m pytest \
  tests/test_wssv_gallery_fastas.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_gallery_session_semantics.py -q
```

Browser checks:

```bash
python tools/prepare_browser_wheel.py
npm run test:web:functional-smoke
npm run test:web:gallery-regeneration:heavy
npm run test:web:vibrio-generate
npx playwright test \
  --config=playwright.functional.config.js \
  tests/web/contracts/session-regenerate-intent.playwright.spec.js
```

Broader gates:

```bash
node --test tests/web/*.test.mjs
python -m pytest tests/ -m "gallery and not slow" --durations=30
npm run test:web:gallery-regeneration
python -m pytest tests/ -m "not slow" --durations=30
ruff check gbdraw/
pytest tests/test_output_comparison.py::TestOutputComparison -v
git diff --check
git status --short --untracked-files=all
```

Do not regenerate `tests/reference_outputs/` unless a separately reviewed core
geometry change is intentional. Gallery-only regenerated SVGs are not a reason
to update core reference outputs.

Final diff audit must separate:

1. production/runtime code;
2. test/helper/config code;
3. independently reviewed scientific-baseline fixtures;
4. source input/provenance changes;
5. session JSON/gzip;
6. source/display SVG;
7. thumbnails/tutorial media;
8. manifest/tutorial/register documentation;
9. CI/deploy workflows.

## 8. Per-example acceptance matrix

Counts below describe the planning baseline. For every comparison-bearing
figure, Phase 4 must replace an informal count-only expectation with the
independent approved settings and sorted ID sets in
`tests/web/fixtures/gallery-scientific-baselines.json`. Counts remain readable
summaries, not the scientific acceptance oracle.

| ID | Tier | Required invariants after unchanged first Generate |
| --- | --- | --- |
| `HmmtDNA_basic_circular` | standard | Circular, one record, 37 features, labels/legend, ticks, GC content, and GC skew retained. |
| `lambda_basic_linear` | standard | Linear, one record, `linearComparisonPlan.mode = none`, zero comparison groups, ruler, labels, and left legend retained. |
| `HmmtDNA_ATskew` | standard | Five ordered slots; AT-skew `nt`, colors, legend; split feature lane; ticks; `ND1`; viewBox and path geometry semantically equivalent. |
| `tobacco-chloroplast` | standard | Import succeeds; retired fields absent; one record, 197 features, four region annotations, three intended slots, and bracket labels retained. |
| `BGC0000708-BGC0000713` | standard | Five records, record labels/subtitles, four color rules, qualifier priority, four comparison groups, the approved scientific-baseline match/orthogroup IDs (77 planning-baseline matches), legend, and applied `og_18 -> livA` interactive override retained. |
| `Vnig_TUMSAT-TG-2018` | heavy | Six records, 6,316 features, grid grouping, and positions `#1@1,#2@1,#3@2,#4@2,#5@2,#6@2` retained. |
| `hepatoplasmataceae_collinear` | heavy | Five records, adjacent collinear intent, four visible comparison groups, the approved scientific-baseline match/block IDs (500 planning-baseline matches), GC content/skew, and record-local geometry retained. |
| `hepatoplasmataceae_orthogroup` | heavy | Five records, adjacent similarity-group intent, four visible groups, the approved scientific-baseline orthogroup/match IDs (577 orthogroups and 1,899 planning-baseline matches), and GC tracks retained. |
| `majanivirus_orthogroup` | heavy | Nine record labels, eight comparison groups, 152 orthogroups, three specific colors/captions, neutral CDS palette, and the repository-owner-approved scientific-baseline IDs retained. |
| `WSSV_genome_comparison` | heavy | Reference record, 20 exact FASTAs, 20 BLAST results, 20 ordered rings, approved series/result identities, labels/colors, feature metadata, and successful Generate retained. |
| `vibrio-harveyi-group-collinear` | vibrio-special | Eleven replicons in five rows, gap 48, adjacent collinear plan, comparison groups/ribbons/legend, repository-owner-approved match/block IDs, and repeated-Generate determinism retained. |

## 9. Required negative tests

The implementation is incomplete unless the gates reject all of these cases:

1. a current Gallery config containing an unknown `config.adv` field;
2. a chloroplast session containing either retired `cli_circular_track_*`
   field;
3. a WSSV `webFiles` FASTA reference with missing or empty resource data;
4. a frozen lazy resource view rejected only because it lacks `.text()`;
5. a compared Gallery figure whose active Linear plan is `none`;
6. a missing AT-skew slot or changed slot order/parameters;
7. a removed chloroplast annotation;
8. a Majanivirus/Vibrio comparison set reduced to zero or an unapproved ID set;
9. a committed BGC editor/orthogroup override dropped during publication
   normalization;
10. a new `GallerySessionExample` without a valid regeneration tier;
11. a browser test that sees import/Generate errors but still records a visual
    attachment and passes;
12. a refresh failure that leaves a subset of Gallery artifacts replaced;
13. a deployment that reaches `public/` assembly before post-refresh replay
    succeeds;
14. a session that is normalized before CLI replay but diverges again after
    `_merge_refreshed_gallery_artifacts()`;
15. a mismatch between the pure current-writer field inventory and the live
    runtime `state.form` / `state.adv` properties;
16. an unknown, out-of-tier, or broad-match value passed through
    `GBDRAW_GALLERY_EXAMPLE_ID`;
17. a CI matrix whose exact selected-ID union differs from the manifest
    inventory;
18. a Majanivirus public download whose bytes differ from its declared source
    TSV, including after a forced transactional refresh failure;
19. a selected valid heavy ID that also loads a non-owning standard or Vibrio
    spec, or defines zero/multiple tests; and
20. a stable aggregate CI job that passes when any static, standard, heavy, or
    special dependency failed, was cancelled, or was skipped;
21. a semantic-oracle mutation to identity, order, geometry, transform, text,
    resolved style, or visibility that passes because a canonical or raster
    fallback happened to match;
22. a Gallery refresh or asset-generation run that rewrites the independent
    scientific-baseline fixture; and
23. a post-refresh session that passes the static contract but fails real
    browser upload/import or first Generate, followed by `public/` assembly or
    upload anyway.

## 10. Risks and controls

| Risk | Control |
| --- | --- |
| Normalizing Gallery config accidentally breaks legitimate user drafts | Keep normalization in an explicitly Gallery-only function/tool and retain divergent-session browser tests. |
| CLI-only values are silently lost | Compare old/new semantic output per example; expose a real Web setting or stop for approval instead of storing a hidden field. |
| Comparison algorithms produce a different scientific figure | Review settings and ID sets; require explicit baseline approval for Majanivirus/Vibrio before artifact refresh. |
| A regenerated artifact silently becomes its own scientific expectation | Keep reviewed settings and sorted ID sets in a non-generated fixture, record source hashes/review reference, and fail if refresh modifies it. |
| The semantic SVG helper flakes or hides a meaningful difference | Version the snapshot contract, match stable identities, preserve meaningful order, start with parsed numeric equality, and mutation-test every required category. |
| Heavy sessions exhaust CI memory | One context and one worker per example; matrix jobs; retain gzip upload path; release Worker/context after each test. |
| Tests pass on a mocked resource shape but production fails | Use actual session files, frozen lazy views, `setInputFiles()`, Pyodide, and Worker staging. |
| Static active-config validation drifts from the browser's same structural contract | Extract one browser-neutral JS contract, make browser import and Node publication use it, and assert its field inventory matches live runtime state. |
| Static validation is mistaken for complete importer acceptance | Limit its claim to schema/config/projection/resource/fixed-point contracts and require real upload/import/Generate against the same post-refresh bytes before publication. |
| Refresh publishes partial output | Preserve the existing all-staged transaction and add failure-path tests around normalization and browser gates. |
| CLI merge recreates a stale draft after an earlier normalization pass | Normalize and validate the final merged document, then require request/config/resource fixed-point equality for all 11 sessions. |
| A matrix row accidentally runs the full heavy corpus or silently skips its target | Select one exact manifest ID through `GBDRAW_GALLERY_EXAMPLE_ID`; reject unknown/out-of-tier IDs and statically prove matrix-union coverage. |
| A valid selected ID is rejected by another tier's loaded spec | Resolve the manifest tier in the Playwright config, set `testMatch` to the single owning spec, and require exactly one defined test. |
| The downloadable Majanivirus TSV drifts from the session recipe | Declare one source-to-public target mapping, copy it inside the Gallery transaction, and test byte equality and rollback. |
| Screenshot/tutorial evidence becomes stale | Audit exact session-backed operations, update the operation register, and recapture only changed truthful states. |
| Byte-level noise causes false failures | Compare normalized semantic structure and explicit scientific IDs, with raw/raster evidence only for diagnosis. |
| A new Gallery entry escapes coverage | Require manifest tier metadata and assert the union of standard/heavy/special browser definitions equals the manifest inventory. |
| Dynamic matrix child names drift outside branch protection | Collapse all static/standard/matrix outcomes into one stable aggregate status and require only that status. |
| CI runs but is not actually required for merge | Trigger both `dev` and `main`; separately obtain authorization to require the stable aggregate status, or report that external gate as pending without claiming merge protection. |

## 11. Definition of done

Mark this plan complete only when all boxes are supported by recorded evidence:

- [ ] The Phase 0 real-browser audit records the exact baseline commit, all five
      session hashes, and complete loaded/generated evidence before production
      edits.
- [ ] All 11 real session files import in Chromium.
- [ ] All 11 unchanged first Generates return `{status: 'ok'}` with no
      `errorLog` or page error.
- [ ] All 11 loaded/generated semantic SVG comparisons pass.
- [ ] The versioned semantic-SVG oracle passes its positive normalization cases
      and rejects every required identity/order/geometry/style/text mutation
      without exact/canonical/raster fallback.
- [ ] All per-example scientific invariants pass.
- [ ] Every comparison-bearing example has a complete independently reviewed
      scientific-baseline entry with settings, sorted IDs, hashes, and review
      reference; refresh and asset generation leave the fixture unchanged.
- [ ] WSSV preserves all 20 exact FASTAs, BLAST results, and rings.
- [ ] Chloroplast contains no retired current-writer fields and retains all four
      annotations.
- [ ] AT skew retains its approved five-slot geometry.
- [ ] Majanivirus and Vibrio match explicitly approved comparison baselines.
- [ ] General user-session divergent drafts remain supported.
- [ ] Browser import and Node publication use the same pure current-writer
      validator, and its field inventory equals live runtime state.
- [ ] Static refresh evidence is reported only for schema/config/projection/
      resource/fixed-point contracts; complete importer acceptance is supported
      by real upload/import/Generate evidence for the exact post-refresh bytes.
- [ ] Every final, post-merge Gallery document has equal committed-request and
      active-next-request semantic signatures, including resource hashes.
- [ ] A second publication-normalization pass changes neither config, request,
      resource bindings, nor referenced resource payloads.
- [ ] Publication normalization, declared download synchronization, and
      complete resource validation are transactional.
- [ ] The declared Majanivirus source and public download TSV are byte-equal.
- [ ] Generated sessions, SVGs, thumbnails, manifest, tutorials, and media are
      mutually consistent.
- [ ] Changed public figures are visually inspected at readable desktop and
      mobile sizes.
- [ ] Pull requests to both `dev` and `main` trigger standard and every exact
      heavy/special replay case.
- [ ] The exact selector rejects invalid IDs, and the standard/heavy/special
      union equals all manifest entries.
- [ ] A selected ID loads only its owning spec and defines exactly one test.
- [ ] The stable `gallery-regeneration-required` job fails unless every static,
      standard, heavy, and special dependency succeeds.
- [ ] Repository rules make that aggregate Gallery status merge-required, or the
      handoff truthfully records External Gate A as pending authorization and
      does not claim merge protection.
- [ ] Deployment records post-refresh hashes, replays those artifacts before
      `public/` assembly, and copies the same bytes before upload.
- [ ] `tests/reference_outputs/` and `examples/gbdraw_social_preview.png` are
      unchanged unless separately authorized.
- [ ] Production, tests, source data, generated artifacts, documentation, and
      workflow diffs were reviewed separately.
- [ ] Focused and broader command evidence is recorded in this document or the
      execution handoff.

## 12. Implementation handoff template

At completion, record:

```text
Implementation commit:
Branch and base:
Session/request versions:

Runtime fix evidence:
- WSSV lazy resource:
- chloroplast import:

Gallery publication evidence:
- known-bad baseline audit commit/session hashes/evidence directory:
- SemanticSvgSnapshot version and oracle mutation tests:
- static all-11 contract:
- static validation scope (no complete-importer claim):
- pure validator/runtime inventory parity:
- post-merge request/config/resource fixed point:
- browser all-11 contract:
- scientific-baseline fixture review reference and refresh-immutability proof:
- Majanivirus approved baseline:
- Majanivirus source/public TSV byte equality and rollback:
- Vibrio approved baseline:
- AT geometry decision:

Generated artifact review:
- sessions:
- SVGs:
- thumbnails:
- tutorial media/register:
- reference outputs unchanged:
- social preview unchanged:

CI/deploy evidence:
- dev/main PR triggers:
- exact-ID matrix inventory coverage:
- selected-ID owning-spec/exactly-one-test evidence:
- `gallery-regeneration-required` dependency/result evidence:
- External Gate A repository-rule status/authorization:
- post-refresh hash inventory/browser replay/public copy equality:

Known residual risks or follow-up work:
```

Proposed implementation commit title:

```text
Fix Gallery session regeneration parity
```

Proposed summary:

```text
- align published Gallery drafts with their committed render state
- support lazy session resources during WSSV FASTA extraction
- validate and replay every Gallery session before publication
```
