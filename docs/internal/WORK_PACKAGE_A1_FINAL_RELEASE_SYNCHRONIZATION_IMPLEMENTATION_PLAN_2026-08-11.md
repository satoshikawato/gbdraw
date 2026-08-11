# Work package A1: final release documentation and evidence synchronization implementation plan

- Date: 2026-08-11
- Status: planned; entry gates are not yet satisfied
- Initial audit baseline: `docs_renovation` / `6f14e2c4fd2a`
- Drafting note: `HEAD` moved in the shared worktree while this plan was being
  written; Phase 0 must record the accepted feature-freeze revision and must
  not reuse the initial audit SHA as candidate evidence
- Source: [gbdraw v0.14.0 release roadmap](gbdraw_v0.14.0_codex_roadmap.md), Work package A1
- Prerequisites: Work package A0; frozen Work packages B–I; Gate 0 feature
  freeze; Work package H1 packaging readiness
- Later gates: Work package H-final, J-RC, J-Final, and K-Publication

Related authorities:

- [Documentation simplification plan](DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md)
- [Work package J implementation plan](WORK_PACKAGE_J_QA_COMPATIBILITY_RELEASE_ENGINEERING_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Documentation evidence runner](../capture/README.md)
- [Session compatibility history](../SESSION_COMPATIBILITY.md)
- [Repository working rules](../../AGENTS.md)
- [Web application working rules](../../gbdraw/web/CLAUDE.md)

This plan turns the roadmap A1 contract into executable phases. At creation it
changes internal planning documentation only. It does not implement a release
feature, regenerate a public figure, change a Gallery artifact, build a
distribution, tag a revision, publish a package, deploy the hosted application,
or create an archive.

## 1. Outcome

A1 binds release documentation and reproducibility evidence to observed
software state instead of to an intended future release.

The candidate-synchronization milestone has one source revision and contains:

- public prose that describes the frozen behavior and actual route status;
- generated CLI reference and compatibility facts that agree with code;
- manifest-declared CLI, Python, GUI, and contract evidence;
- reproducible screenshots, figures, Gallery sessions, assets, and operation
  media;
- stable checksum-bound inputs with cross-platform LF checkout behavior;
- a provisional built-wheel clean-install run that exposes documentation or
  packaging mistakes before H-final;
- a ledger naming every command, result, artifact, deviation, and remaining
  external dependency.

A1 does not certify the final wheel or sdist by itself. After A1 stops changing
shipped files, H-final builds the exact artifacts and J-RC/J-Final repeat the
declared installation and reproduction gates against those hashes.

Before H-final, A1 classifies every release identifier owner as shipped or
publication-only. A shipped owner is finalised with either an observed/reserved
identifier or truthful non-live wording before H-final builds. After K has
performed explicitly authorised release, deployment, and archive actions, A1
has one small verification re-entry phase. It may update only publication-only
owners, verifies the observed live version, URL, hash, and DOI, and hands the
result to K-Publication. Any required change to README, package metadata,
release notes, CHANGELOG, `CITATION.cff`, ABOUT, or another shipped file creates
a patch candidate and returns through H-final and J.

## 2. Scope

### 2.1 Goals

1. Record the exact candidate source revision and invalidate evidence when that
   revision or a shipped artifact changes.
2. Make README, installation, tutorial, reference, FAQ, privacy, release, and
   citation statements agree with the frozen implementation.
3. Run literal public CLI/Python recipes from clean directories against an
   installed wheel, not only an editable checkout.
4. Regenerate final GUI and Gallery evidence through its committed owner.
5. Verify package, session, request, CLI, Python, Gallery, and documentation
   version facts from one inventory.
6. Leave a complete handoff for H-final and J, including the gates that must be
   repeated against exact artifacts.
7. Verify observed live identifiers after K release/archive actions and update
   only explicitly inventoried publication-only owners without changing the
   accepted source distribution.

### 2.2 In scope

- `README.md`, `docs/INSTALL.md`, `docs/QUICKSTART.md`, `docs/DOCS.md`;
- Tutorials and technical documentation under `docs/TUTORIALS/` and
  `docs/REFERENCE/`;
- `docs/CLI_Reference.md` and its generator;
- `docs/FAQ.md`, `docs/ABOUT.md`, privacy/analytics documentation;
- `CHANGELOG.md`, release notes, `CITATION.cff`, and citation instructions;
- `docs/scenarios/manifest.json` and manifest-owned capture/recipe evidence;
- documentation screenshots and generated CLI/Python figures;
- Gallery sessions, source SVGs, rendered examples, thumbnails,
  `examples.json`, tutorial JSON, and operation media;
- generator-owned session fixtures under `tests/test_inputs/` when the final
  writer requires their refresh;
- package/session/request/table/resource version inventories;
- provisional and exact-artifact installation/reproduction evidence;
- documentation links, route availability, support statements, final live
  identifier verification, and publication-only closeout.

### 2.3 Out of scope

- implementing or redesigning Work packages B–I;
- creating a fifth public documentation route;
- executing illustrative placeholder fences as if they were complete recipes;
- reacquiring every legacy input when no public reproduction claim requires it;
- refreshing every GUI image on every ordinary pull request;
- manually editing generator-owned Gallery or capture outputs;
- changing `tests/reference_outputs/` without an independently reviewed
  rendering-contract change;
- modifying `examples/gbdraw_social_preview.png`;
- refreshing the Web cache-bust token before a deployable-bundle action;
- tagging, pushing, publishing to TestPyPI/PyPI, updating Bioconda, deploying,
  creating a Zenodo archive, revising the preprint, or submitting a manuscript.

## 3. Planning audit

These facts describe the planning baseline only. Phase 0 must measure them
again; they are not release evidence.

| Area | Observed state at planning | Consequence for A1 |
| --- | --- | --- |
| Branch relation | Local `main` was an ancestor of `docs_renovation`; the development branch was 36 commits ahead at `6f14e2c4` | There was no two-sided divergence, but remote refs and the final candidate were not frozen |
| Package/version | Package `0.14.0b0`, session writer 40, request schema 5 | Work packages C/D plan a coordinated later writer, so current values must not be copied into final prose |
| Public structure | Four owners: Tutorials, Technical documentation, FAQ, Gallery | Preserve the structure unless a distinct reader question has no owner |
| Tutorials | Ten projects across GUI/CLI/Python, 30 public tutorial pages | Final behavior must remain cross-surface where the manifest declares parity |
| Evidence inventory | 75 manifest scenarios: 25 Playwright, 23 CLI recipes, 15 Python recipes, 12 contracts | `run_all.py` covers the 63 executable scenarios; contract tests cover the other 12 |
| Gallery | Eleven ready examples with session/request metadata at the audit point | Final schema/UI changes may require generator-owned session, asset, and media refresh |
| Stable inputs | 36 files in 12 fixture groups with size/SHA/semantic metadata | Checkout bytes, not only logical text, must remain stable |
| Line endings | `* text=auto eol=lf` plus an `autocrlf=true` fixture regression | Retain the policy and rerun it on the candidate |
| Documentation gaps | README omitted PyPI and Python API positioning; INSTALL said three routes while listing four; PyPI was not live | Route state must be factual at candidate and post-release closeout |
| Full final run | Targeted evidence existed, but no one green sweep against final frozen behavior | A1 must run the final nightly and Gallery gates once, not infer completion from older runs |

Known authority defect to resolve before Gallery regeneration:
`gbdraw/web/CLAUDE.md` names `tools/build_web_gallery.py`, which does not exist.
The current generator owners are `tools/refresh_gallery_sessions.py` and
`tools/prepare_interactive_gallery_assets.py`. Phase 0 must verify the current
owners and correct the guidance if the mismatch still exists.

## 4. Fixed decisions

### D1. A0 and A1 are separate gates

A0 establishes the integration and documentation contract before feature work.
A1 synchronizes final content and evidence after feature contracts freeze. A1
does not repeat a completed A0 structural merge. It verifies A0 ancestry and
drift; if structural divergence has appeared, A0 reopens before A1 continues.

### D2. Heading order does not imply execution order

A1 is documented beside A0 because they form one documentation workstream. The
roadmap dependency order controls execution: B–I, H1, and Gate 0 precede A1
candidate synchronization; H-final and J follow it.

### D3. One evidence pass belongs to one revision

Every pass records:

- source commit and dirty-tree boundary;
- package, session, request, support, and scope versions;
- tool/runtime versions;
- exact command or CI run;
- produced artifact names, sizes, and SHA-256 where applicable;
- inspected visual artifacts;
- failures, approved deviations, and remaining risk.

A source change does not automatically invalidate every gate. It invalidates
the gates whose input, output, or claim can change. A version-only change still
invalidates version concordance, artifact build, installed-artifact checks, and
release-facing version text.

### D4. “All documentation commands” has a testable meaning

The release gate covers:

- all 23 manifest-declared CLI recipes;
- all 15 manifest-declared Python recipes;
- all 25 manifest-declared GUI scenarios at their final cumulative tier;
- all 12 manifest-declared contract scenarios;
- clean base-wheel and export-extra installation smokes;
- generated CLI help, Python signatures, session/request examples, local links,
  and applicable static option/syntax contracts.

An illustrative block containing placeholders such as `genome.gb` is not a
literal recipe unless the manifest marks it as executable. It still receives
the static validation owned by its documentation contract.

### D5. Source synchronization and exact-artifact certification are different

A1 may build a provisional wheel after source documentation and packaged Web
assets are synchronized. This catches package-boundary errors early.

H-final later builds the wheel, sdist, and deployable Web bundle once from the
candidate revision. J downloads and verifies those exact files. A1 does not
replace that evidence with a locally rebuilt wheel.

### D6. Generated artifacts have one owner

| Artifact | Owner | A1 rule |
| --- | --- | --- |
| CLI help blocks | CLI parser plus `tools/update_cli_reference_help.py` | Never hand-edit generated help |
| CLI/Python figures | Literal marked block plus recipe runner | Regenerate in clean temporary directories |
| GUI screenshots | Manifest plus `docs/capture/flows/` and `run_all.py` | Capture through the real loopback-only Web flow |
| Gallery sessions/assets | `refresh_gallery_sessions.py` and `prepare_interactive_gallery_assets.py` | Regenerate only after schema/render behavior freezes |
| Gallery operation media | Tutorial JSON metadata plus `capture_gallery_tutorial_screenshots.py` | Use declarative browser capture; strict validation is mandatory |
| Reference SVGs | Reference-output test owner | Read-only unless a reviewed rendering change explicitly authorizes refresh |
| Browser wheel | `tools/prepare_browser_wheel.py` | Generated and ignored; do not commit or refresh cache-bust here |
| Social preview | Maintainer | Never modify |

### D7. Public prose has one semantic owner

| Reader question | Owner |
| --- | --- |
| Where do I start? | `README.md` and `docs/DOCS.md` as short routers |
| Which installation route is available? | `docs/INSTALL.md`; README carries only the concise route summary |
| How do I complete a task? | Tutorials |
| What is the exact API/CLI/session behavior? | Technical documentation |
| What should I choose or how do I recover? | FAQ |
| What finished results are possible? | Gallery |
| What changed in this release? | Release notes and CHANGELOG |
| How should this software/version be cited? | `CITATION.cff` and `docs/ABOUT.md` |

Do not duplicate complete instructions across these owners. Link to the
authority.

### D8. Installation routes use observed states

Use one of these states explicitly:

- `live`: the documented URL/package can be observed and smoke-tested;
- `release candidate`: an exact candidate artifact exists but the stable route
  is not yet public;
- `planned for v0.14.0`: no candidate/live artifact exists;
- `unavailable`: the route was cut or failed qualification.

Before publication, `python -m pip install gbdraw` must not be presented as a
verified stable PyPI route. After K publication, A1 marks it live only after the
observed project URL and clean-install smoke pass.

### D9. Persisted-format facts come from implementation constants

The current package version, session writer/readers, request schema/readers, and
Gallery inventories are extracted from code and contract tests. Do not retain a
branch-only intermediate writer in user documentation or compatibility tables.
The coordinated C/D writer is documented only after its migration and replay
gates pass.

### D10. Candidate citations and live archive identifiers are separate

A1 candidate synchronization may keep the observed current preprint citation
and clearly mark the release archive identifier as pending. K owns tag/archive
creation. If K can reserve a version DOI before the final source transaction,
A1 records that observed reservation before H-final. Otherwise shipped citation
owners use truthful non-live wording that remains correct after publication.
After K records a live DOI, A1 verifies shipped owners and may update only an
inventoried publication-only citation owner before K-Publication. Updating a
shipped citation owner requires a patch candidate and H-final/J
requalification.

### D11. Stable input scope is proportional

Every claimed reproducible input is accession-pinned or repository-tracked and
checksum-bound. New or replaced authoritative inputs require retrieval metadata
and mirror verification. Missing byte-for-byte upstream provenance for an
unchanged legacy fixture is not by itself a `v0.14.0` blocker unless the public
claim depends on retrieving that exact upstream revision.

### D12. J fixes reopen evidence

Every P0/P1 and approved P2 fix identifies its affected A1/H-final/J gates.
Rendering, UI, schema, package-content, public-name, output, or documentation
changes rerun their owner gates. The prior pass remains historical evidence and
is not relabelled as evidence for the new revision.

### D13. One final release-facing source transaction has named owners

The last source transaction before H-final is ordered and single-owner:

1. H changes package/version metadata.
2. A1 changes release notes, citation metadata, and public prose.
3. The declared generators refresh derived files.
4. The candidate revision is recorded and A1 Phases 0–5 close.
5. H-final builds the exact artifacts once.
6. J certifies those artifacts and does not silently take ownership of A1 prose.

Any J-requested source fix reopens the affected owner and repeats the downstream
transaction. J-RC therefore requires the A1 **candidate-synchronization
milestone** (Phases 0–5), not overall A1 completion after publication.

## 5. Ownership and change matrix

| Surface | Primary files or systems | Intended A1 result | Required evidence |
| --- | --- | --- | --- |
| Branch/source | Git refs and candidate commit | A0 relationship verified; one candidate revision recorded | ref commands, status inventory, reviewed diff |
| Package versions | `pyproject.toml`, `gbdraw/version.py`, Web wheel/version owners | One supported version statement | build metadata and version contracts |
| Session/request | Python/Web constants, compatibility docs, Gallery sessions | Writer/readers and artifacts agree | codec/session tests and Gallery inventory |
| README/install | `README.md`, `docs/INSTALL.md` | Factual live/candidate/planned routes | clean-install evidence and route review |
| Navigation | `docs/DOCS.md`, Tutorial index | Four owner model retained | local-link and owner contracts |
| CLI | parser, help generator, CLI reference, recipes | Generated help and literal commands agree | generator check plus all CLI recipes |
| Python | public exports/signatures, reference, recipes | Public signatures/examples agree | contract tests plus all Python recipes |
| Web/tutorials | public GUI tutorials, capture manifest/flows/images | Final control names/actions/results agree | nightly capture check and visual review |
| Gallery | sessions, sources, examples, thumbnails, tutorials/media | Final schema/output/tutorial steps agree | generators, strict media check, browser tests |
| Privacy | Web consent behavior plus public explanation | Frozen retained/cut decision described accurately | offline/network tests and doc contracts |
| Release/citation | CHANGELOG, release notes, ABOUT, CITATION | Shipped wording final before H-final; observed live identifiers verified after K | version/link/citation checks |
| Distribution | H provisional/H-final artifacts | Provisional mistakes caught; exact hashes handed to J | artifact manifest and installed smokes |

## 6. Entry gates

A1 Phase 0 remains blocked until all of the following are observed:

- A0 ledger is accepted against a recorded integration baseline.
- Work packages B–G have completed or recorded an approved release cut.
- C/D have one complete persisted-format decision with no partial writer.
- E/F final public control names and behavior are frozen.
- H1 has a named, available command owner for unique-output wheel/sdist build,
  inspection, and release automation. If the H1 implementation plan or path is
  absent, A1 blocks instead of inventing a replacement release procedure.
- I has a frozen retained implementation or an explicit release cut; its
  privacy text matches that decision.
- Gate 0 has recorded the feature-freeze revision used as the A1 comparison
  boundary.
- The `R-CLI-01` GUI-help assertion has a truthful owner. Either `gbdraw gui
  --help` is non-blocking and covered by a current contract, or the manifest and
  reference claim are narrowed to what is actually tested. The current
  circular/linear-only generator is not evidence for GUI help.
- No unreviewed reference-output update is pending.
- The worktree contains no unidentified generated artifacts or unrelated edits
  that would be swept into an A1 candidate.

If one condition fails, record `blocked` with the owning work package. Do not
weaken A1 acceptance to proceed.

## 7. Execution phases

### Phase 0: lock the candidate boundary

Status: not started

#### Owned state

- A1 header/status and evidence ledger;
- candidate/support/scope manifest references;
- branch and worktree inventory;
- current generator/tool owner inventory.

#### Work

1. Re-read repository instructions and all completed upstream ledgers.
2. Read the accepted A0 ledger to derive `a1_target_ref`; do not assume that its
   integration target is named `main`. Fetch that ref and compare it with the
   recorded A0 baseline.
3. Record `a1_feature_freeze_commit` from Gate 0, then record candidate branch,
   source SHA, intended integration target, merge base, ahead/behind counts,
   dirty files, and untracked files. Candidate evidence diffs start at the
   feature-freeze commit, not at an unrelated whole-branch base.
4. Classify dirty files as A1-owned, upstream-owned, generated, or unrelated.
5. Re-extract package/session/request/table/resource/Gallery versions.
6. Confirm the current CLI/Python/GUI/contract scenario counts from the
   manifest.
7. Verify the actual Gallery generator owners and correct stale developer
   guidance before any Gallery write.
8. Confirm that the named H1 plan and unique-output build/inspection commands
   exist. Stop if ownership is unavailable.
9. Resolve the false-positive risk in `R-CLI-01`: verify a non-blocking,
   directly tested GUI-help contract or narrow the manifest/reference claim.
   Do not invoke a server-starting `gbdraw gui --help` as an unattended probe.
10. Record Python, Node, Playwright, browser, font, locale, timezone, and
    relevant system-library versions.

#### Baseline commands

```bash
git status --short --branch
git fetch --prune origin
a1_target_ref="<accepted A0 integration target ref>"
a1_feature_freeze_commit="<Gate 0 feature-freeze commit>"
git rev-parse --short=12 HEAD
git rev-parse --short=12 "$a1_target_ref"
git merge-base --is-ancestor "$a1_target_ref" HEAD
git rev-list --left-right --count "$a1_target_ref"...HEAD
git diff --name-status "$a1_feature_freeze_commit"..HEAD
git diff --check "$a1_feature_freeze_commit"..HEAD

rg -n 'gui_help_current|gui --help' \
  docs/scenarios/manifest.json \
  tests/test_documentation_reference_contracts.py \
  tools/update_cli_reference_help.py \
  gbdraw/cli.py

command -v playwright
playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

`git fetch` is read-only with respect to repository content but still requires
network access. It does not authorize merge, push, deployment, or release.

#### Exit gate

- Exact source/ref/tool inventory is recorded.
- A0 relation is still valid or A0 has been reopened and completed.
- Gate 0 feature-freeze revision and the H1 command owner are available.
- The GUI-help claim has direct, non-blocking evidence or has been narrowed
  truthfully.
- No unidentified dirty file can enter the candidate.
- All entry-gate owners are complete or have an explicit approved cut.

### Phase 1: synchronize public truth and ownership

Status: not started

#### Owned files

- README, INSTALL, DOCS, QUICKSTART;
- Tutorial/reference indexes and affected pages;
- FAQ, ABOUT, privacy documentation;
- compatibility documentation;
- CHANGELOG, release notes, CITATION;
- public-owner and scenario metadata when ownership actually changed.

#### Work

1. Build a truth matrix with one row per release claim: owner, shipped versus
   publication-only classification, source constant or behavior, candidate
   state, evidence, identifier strategy, and post-release verification need.
2. Apply `keep`, `merge`, `delete`, or `new` to every proposed public-page
   change. A new page requires a reader question that no current owner can
   answer clearly.
3. Correct route counts and route availability. Include Hosted Web, Bioconda,
   PyPI, and source, but distinguish live, candidate, planned, and unavailable.
4. Add Python API positioning to the concise product-surface description where
   missing.
5. Synchronize final feature names/defaults only after their upstream contracts
   pass.
6. Synchronize privacy text with the retained/cut Work package I result.
7. Keep release notes historical and release-specific; keep exact current
   compatibility in its technical owner.
8. Record the Work package D manual-placement Gallery teaching decision as
   `yes` or `no`. Default to no new Gallery example unless it answers a reader
   question that the existing Gallery and tutorials do not. If yes, name the
   Gallery-quality source, session/asset generator, tutorial owner, and visual
   review gate before Phase 3.
9. For every unknown post-release version/URL/DOI, either use truthful non-live
   wording that can ship unchanged or have K reserve the identifier before the
   final source transaction. Do not schedule a shipped-file edit after the
   accepted tag/artifact build.
10. Give `docs/RELEASE_NOTES_0.14.0b0.md` one explicit disposition: retain it as
    immutable beta history, rename it with all inbound links updated, or
    supersede it with a final-version owner. Inventory `docs/DOCS.md`,
    `CHANGELOG.md`, `docs/scenarios/manifest.json`,
    `tests/test_documentation_contracts.py`, and package/link contracts before
    choosing.
11. Search for roadmap-only feature names leaked into public current-capability
    prose.

#### Narrow verification

```bash
python tools/update_cli_reference_help.py --check
rg -n 'RELEASE_NOTES_0\.14\.0b0|0\.14\.0b0' \
  docs/DOCS.md CHANGELOG.md docs/scenarios/manifest.json \
  tests/test_documentation_contracts.py
python -m pytest \
  tests/test_documentation_scenario_contracts.py \
  tests/test_tutorial_documentation_contracts.py \
  tests/test_documentation_reference_contracts.py \
  tests/test_documentation_contracts.py \
  tests/test_documented_recipes.py \
  tests/test_reproduce_examples.py::test_public_markdown_local_targets_exist \
  -q
```

#### Exit gate

- Every public release claim has one owner and evidence source.
- Every identifier owner is classified as shipped or publication-only, and the
  beta release-note file has an explicit reviewed disposition.
- The manual-placement Gallery teaching decision and receiving owner are
  recorded.
- No unavailable route or unimplemented feature is presented as current.
- No unnecessary public page or duplicated full procedure was added.
- Static documentation, owner, compatibility, and local-link contracts pass.

### Phase 2: synchronize CLI/Python reference and executable evidence

Status: not started

#### Owned artifacts

- generated CLI help;
- affected CLI/Python tutorial and technical pages;
- literal marked blocks;
- generated CLI/Python SVGs and declared non-SVG outputs;
- scenario expected-output metadata;
- `tools/reproduce_examples.py` installed-artifact execution mode and its
  focused tests.

#### Work

1. Regenerate CLI help through its generator when the frozen parser changed.
2. Review the help diff against parser/default decisions, then run `--check`.
3. Run all CLI and Python recipe checks from clean temporary directories.
4. When an intentional final behavior change causes a mismatch, regenerate only
   through the runner, review semantics and visual output, then rerun `--check`.
5. Confirm each recipe copies only declared checksum-bound inputs and rejects
   undeclared outputs.
6. Confirm heavy LOSATP scenarios complete; do not treat ordinary CI exclusion
   as final release evidence.
7. Before the Phase 4 candidate freeze, add and test
   `reproduce_examples.py --installed-artifact-prefix <venv-prefix>`. In that
   mode every child `python -m gbdraw.cli` process runs from a runner-created
   directory outside `PROJECT_ROOT`, probes its own module origin, and fails
   unless `gbdraw.__file__` is below the supplied prefix. Preserve absolute
   clean-archive input/manifest resolution and keep ordinary source mode
   unchanged.

#### Commands

```bash
python tools/update_cli_reference_help.py
python tools/update_cli_reference_help.py --check

python docs/recipes/run_cli_scenarios.py --all --check
python docs/recipes/run_python_scenarios.py --all --check
python -m pytest tests/test_reproduce_examples.py \
  -k installed_artifact -q
```

The first help command is a write command. Run it only when parser/help changes
are expected; otherwise use `--check` alone.

#### Exit gate

- All 23 CLI and 15 Python recipes pass from clean directories.
- Every changed generated output has a reviewed source recipe and declared
  inputs.
- CLI help, Python signatures, package defaults, and prose agree.
- The installed-artifact figure-runner mode has a child-context origin test and
  is committed before Phase 4 records or builds the candidate.

### Phase 3: synchronize GUI and Gallery evidence

Status: not started

#### Owned artifacts

- `docs/capture/flows/`, scenario metadata, and `docs/images/`;
- Gallery session JSON, source SVG, rendered example, thumbnail, and
  `examples.json`;
- Gallery tutorial JSON and operation WebP media.
- generator-owned test-input sessions refreshed by
  `refresh_gallery_sessions.py`.

#### Work

1. Prepare the ignored browser wheel without refreshing the cache-bust token.
2. Run the final GUI `--check` first to identify stale evidence.
3. Regenerate only through the scenario owner after confirming the behavior
   change is intentional.
4. Finish with one cumulative nightly check covering all 25 GUI scenarios and
   all CLI/Python recipes invoked by the orchestrator.
5. After the final session/request writer freezes, regenerate Gallery sessions
   and derived assets through `refresh_gallery_sessions.py`.
6. Apply the recorded manual-placement Gallery disposition. When it is `no`,
   record why the existing tutorial/reference owner is sufficient. When it is
   `yes`, generate the example only from the named Gallery-quality source and
   declared owners; do not promote a test smoke diagram.
7. Regenerate affected operation media from tutorial JSON; run strict validation
   across all ready examples.
8. Run Gallery session-replay/tutorial browser tests.
9. Inspect every changed public image at readable scale. Check input identity,
   labels, legends, comparison context, controls, expected state, crop, and
   visual finish.
10. Confirm `examples/gbdraw_social_preview.png` and unrelated reference outputs
   did not change.

#### Commands

```bash
python tools/prepare_browser_wheel.py
python docs/capture/run_all.py --tier nightly --check

# Run these writes only for reviewed, intentionally changed evidence.
python docs/capture/run_all.py --scenario <affected-id> --tier <affected-tier>
python docs/capture/run_all.py --tier nightly --check

# Run only after the final session/request writer and Gallery content freeze.
python tools/refresh_gallery_sessions.py
python tools/capture_gallery_tutorial_screenshots.py \
  --example <affected-gallery-id>
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict

python -m pytest \
  tests/test_gallery_capture_contracts.py \
  tests/test_gallery_session_semantics.py \
  tests/test_refresh_gallery_sessions.py \
  -q

npx playwright test \
  tests/web/gallery-tutorial.playwright.spec.js \
  tests/web/gallery-session-regeneration.playwright.spec.js \
  --project=chromium --workers=1
```

When a sandboxed Chromium launch fails with an operating-system sandbox error,
rerun the same command with the required approval. Do not substitute a static
contract for a required browser flow.

`refresh_gallery_sessions.py` is a generator write. Run it only after recording
the pre-generation diff and confirming that the final writer is authoritative.

#### Exit gate

- The cumulative nightly check passes, including long-running scenarios.
- Gallery sessions/assets/media pass their generator and strict validation.
- All changed visual artifacts have a recorded visual review.
- No protected or unrelated generated artifact changed.

### Phase 4: build and exercise a provisional candidate

Status: not started

This phase uses the H1 build path after A1 source and shipped Web/Gallery files
have stopped changing. It is early fault detection, not H-final/J evidence.

#### Work

1. Record a source revision for the provisional candidate. If a local commit is
   not authorized or the tree is not clean, do not claim a revision-bound pass.
2. Build wheel and sdist once through the H1 path into a unique temporary
   output directory. Resolve exactly one wheel and one sdist by explicit path;
   never glob the repository's historical `dist/` directory.
3. Inspect both artifacts' contents, metadata, size, hashes, entry points,
   notices, native/Wasm/browser
   assets, caches, local paths, credentials, tests, and accidental data.
4. Create a base-only clean virtual environment and prove import, help,
   Circular/Linear SVG, and session replay without development or export
   extras. Run outside the checkout with `PYTHONPATH` unset, user site disabled,
   `pip check` green, and `gbdraw.__file__` below that environment's prefix.
5. Create an export-extra environment and run the complete CLI/Python recipe
   sets against the installed wheel.
6. Inspect the wheel's offline GUI assets. Separately run the repository Web
   bundle's offline startup smoke after preparing its ignored browser wheel;
   H-final/J owns installed-artifact GUI entry-point qualification.
7. Run `reproduce_examples.py --group all --list`. Resolve every missing input
   or remove/narrow the public claim before continuing. Then reproduce the
   declared public figures from a `git archive`-equivalent clean candidate tree
   while importing the installed wheel and using no checkout-only cache/path.
   Use the installed-artifact mode completed in Phase 2. If the candidate's
   `--help` does not expose that mode or its child-context test is not green,
   stop and return to Phase 2; do not modify source after recording/building the
   Phase 4 candidate. A direct parent-process assertion or `PYTHONSAFEPATH`
   alone is not a substitute for this child-mode contract.
8. Record artifact filename, size, SHA-256, environment, command, reproduction
   report, semantic comparison, and visual-review result.

#### Representative Bash commands

Create unique temporary directories; do not reuse `HOME` or a broad workspace
path.

```bash
a1_checkout_root="$(pwd -P)"
test -z "$(git status --porcelain --untracked-files=all)"
a1_candidate_commit="$(git rev-parse HEAD)"
a1_build_root="$(mktemp -d /tmp/gbdraw-a1-build.XXXXXX)"
python -m build --outdir "$a1_build_root/dist"

mapfile -t a1_wheel_paths < <(
  find "$a1_build_root/dist" -maxdepth 1 -type f -name 'gbdraw-*.whl' -print
)
mapfile -t a1_sdist_paths < <(
  find "$a1_build_root/dist" -maxdepth 1 -type f -name 'gbdraw-*.tar.gz' -print
)
test "${#a1_wheel_paths[@]}" -eq 1
test "${#a1_sdist_paths[@]}" -eq 1
wheel_path="${a1_wheel_paths[0]}"
sdist_path="${a1_sdist_paths[0]}"
sha256sum "$wheel_path" "$sdist_path"
python -m zipfile -l "$wheel_path"
tar -tzf "$sdist_path"

a1_base_root="$(mktemp -d /tmp/gbdraw-a1-base.XXXXXX)"
python -m venv "$a1_base_root/venv"
a1_base_python="$a1_base_root/venv/bin/python"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$a1_base_python" -m pip install "$wheel_path"
(
  cd "$a1_base_root"
  env -u PYTHONPATH PYTHONNOUSERSITE=1 "$a1_base_python" -m pip check
  env -u PYTHONPATH PYTHONNOUSERSITE=1 "$a1_base_python" -c \
    'import pathlib, sys, gbdraw; p = pathlib.Path(gbdraw.__file__).resolve(); q = pathlib.Path(sys.prefix).resolve(); assert p.is_relative_to(q), (p, q); print(gbdraw.__version__, p)'
  env -u PYTHONPATH PYTHONNOUSERSITE=1 \
    PATH="$a1_base_root/venv/bin:$PATH" gbdraw --help
  for a1_scenario in T-CLI-01 T-CLI-02 T-CLI-11; do
    env -u PYTHONPATH PYTHONNOUSERSITE=1 \
      PATH="$a1_base_root/venv/bin:$PATH" \
      "$a1_base_python" \
      "$a1_checkout_root/docs/recipes/run_cli_scenarios.py" \
      --scenario "$a1_scenario" --check
  done
)

a1_export_root="$(mktemp -d /tmp/gbdraw-a1-export.XXXXXX)"
python -m venv "$a1_export_root/venv"
a1_export_python="$a1_export_root/venv/bin/python"
env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  "$a1_export_python" -m pip install "${wheel_path}[export]"
(
  cd "$a1_export_root"
  env -u PYTHONPATH PYTHONNOUSERSITE=1 "$a1_export_python" -m pip check
  env -u PYTHONPATH PYTHONNOUSERSITE=1 "$a1_export_python" -c \
    'import pathlib, sys, gbdraw; p = pathlib.Path(gbdraw.__file__).resolve(); q = pathlib.Path(sys.prefix).resolve(); assert p.is_relative_to(q), (p, q); print(gbdraw.__version__, p)'
  env -u PYTHONPATH PYTHONNOUSERSITE=1 \
    PATH="$a1_export_root/venv/bin:$PATH" \
    "$a1_export_python" \
    "$a1_checkout_root/docs/recipes/run_cli_scenarios.py" --all --check
  env -u PYTHONPATH PYTHONNOUSERSITE=1 \
    PATH="$a1_export_root/venv/bin:$PATH" \
    "$a1_export_python" \
    "$a1_checkout_root/docs/recipes/run_python_scenarios.py" --all --check
)

python "$a1_checkout_root/tools/verify_gui_offline.py" inspect-wheel "$wheel_path"
python "$a1_checkout_root/tools/verify_gui_offline.py" check-assets
python "$a1_checkout_root/tools/verify_gui_offline.py" smoke-test

a1_archive_root="$(mktemp -d /tmp/gbdraw-a1-archive.XXXXXX)"
a1_repro_root="$(mktemp -d /tmp/gbdraw-a1-repro.XXXXXX)"
git archive --format=tar \
  --output="$a1_build_root/candidate.tar" "$a1_candidate_commit"
tar -xf "$a1_build_root/candidate.tar" -C "$a1_archive_root"
(
  cd "$a1_repro_root"
  env -u PYTHONPATH PYTHONNOUSERSITE=1 \
    PATH="$a1_export_root/venv/bin:$PATH" \
    "$a1_export_python" \
    "$a1_archive_root/tools/reproduce_examples.py" \
    --group all --list \
    --installed-artifact-prefix "$a1_export_root/venv" \
    --missing-report "$a1_repro_root/missing-inputs.json"
  env -u PYTHONPATH PYTHONNOUSERSITE=1 \
    PATH="$a1_export_root/venv/bin:$PATH" \
    "$a1_export_python" \
    "$a1_archive_root/tools/reproduce_examples.py" \
    --group all \
    --installed-artifact-prefix "$a1_export_root/venv" \
    --output-root "$a1_repro_root/output" \
    --missing-report "$a1_repro_root/reproduction-report.json"
)

python -m pytest "$a1_checkout_root/tests/test_reproduce_examples.py" -q
```

The H1 implementation may wrap these commands, but it must retain the unique
output directory, exact wheel/sdist paths, clean working directory, sanitized
import environment, installed-artifact child-runner mode, and module-origin
assertion from every child working context. Resolve local-wheel extra syntax
explicitly if the platform requires it. Do not replace the installed wheel with
`pip install -e .`.

Windows/macOS exact support evidence belongs to H-final/J. Phase 4 records the
local provisional environment and catches documentation/package integration
errors before that matrix runs.

#### Exit gate

- Provisional wheel and sdist correspond to the recorded source revision.
- Base installation works without accidentally receiving export/dev extras.
- Complete CLI/Python recipes pass against the installed wheel with export
  support where required.
- The public-figure inventory has no unexplained missing input, and clean-archive
  reproduction has a reviewed report and visual/semantic evidence.
- Wheel offline-asset inspection and the separate repository Web-bundle offline
  smoke pass; H-final/J retains installed GUI entry-point qualification.
- No packaging failure is hidden by an editable install.

### Phase 5: integrated A1 candidate gate

Status: not started

#### Work

1. Run documentation, recipe, Gallery, session/request, packaging, offline, and
   local-link tests in their CI-aligned partitions.
2. Run Node and required Playwright tests.
3. Run reference-output comparison read-only.
4. Audit production, test, public documentation, internal documentation,
   generated-artifact, and binary/image diffs separately.
5. Confirm the release scope/support/compatibility manifests and A1 ledger name
   the same candidate revision.
6. Record every deferred P2 and every external K dependency explicitly.

#### Commands

```bash
python -m pytest tests/ \
  -m "not slow and not (recipe or gallery or browser)" \
  --durations=30

python -m pytest tests/ \
  -m "recipe and not slow" \
  --durations=30

python -m pytest tests/ \
  -m "gallery and not slow" \
  --durations=30

python tools/prepare_browser_wheel.py
node --test tests/web/*.test.mjs
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz
python -m pytest tests/ \
  -m "browser and not slow" \
  --durations=30
npx --no-install playwright test --project=chromium --workers=1

python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
git diff --name-status "$a1_feature_freeze_commit"..HEAD
git diff --check "$a1_feature_freeze_commit"..HEAD
git status --short
```

Long-running final documentation scenarios remain mandatory even if the
ordinary `not slow` partitions exclude them. Their Phase 2/3 evidence must refer
to the same revision.

#### Exit gate

- All A1 acceptance rows have observed evidence or an explicitly external later
  owner.
- No unexplained source, generated, image, binary, or package-content diff
  remains.
- The ledger records one candidate revision and its provisional artifact hashes.
- Status may advance to `candidate synchronization complete; H-final/J-RC
  pending`. Do not mark all of A1 complete.

### Phase 6: H-final and J handoff

Status: not started

#### Handoff contract

1. Execute D13's final source transaction: H owns package/version metadata; A1
   owns release notes, citation metadata, and public prose; declared generators
   own derived files.
2. Rerun every Phase 0–5 gate invalidated by that transaction, record the
   resulting revision, and only then close the A1 candidate-synchronization
   milestone. This milestone, not overall post-publication A1 completion,
   satisfies the A1 prerequisite for J-RC.
3. Freeze that recorded revision. H-final builds wheel, sdist, and deployable
   Web bundle once and records exact
   hashes.
4. J-RC downloads those exact files and repeats clean-install, support-matrix,
   offline, documentation/reproduction, and full release gates. J certifies A1
   prose; it does not perform an untracked second final-doc rewrite.
5. If J finds a defect, map the fix to affected A1/H/J gates and create a new
   candidate.
6. A final version/documentation change is a new candidate. Repeat version
   concordance, H-final build/inspection, clean install, full tests, and
   reproduction.
7. A1 ledger links the accepted J-RC/J-Final evidence but does not claim that a
   locally rebuilt artifact is identical.

#### Exit gate

- J evidence names the exact A1 source revision and H-final artifact hashes.
- Every post-A1 fix has a regression test and rerun map.
- No P0/P1 or unapproved P2 remains.
- This gate establishes release readiness only; it does not authorize an RC tag,
  push, upload, deployment, or archive.

### Phase 7: observed live-identifier closeout

Status: not started; blocked until authorized K release/archive actions succeed

#### Owned state

- K's observed tag/package/deployment/archive evidence references;
- the shipped-versus-publication-only identifier matrix;
- any explicitly inventoried publication-only metadata owner;
- the final A1 verification ledger and K-Publication handoff.

README, INSTALL, release notes, CHANGELOG, `CITATION.cff`, ABOUT, package
metadata, and packaged Web assets are verification-only in this phase because
they were frozen before H-final.

#### Work

1. Obtain K's observed tag, commit, package URLs, artifact hashes, hosted bundle
   identity, Bioconda state, and Zenodo version DOI.
2. Compare those values with every shipped and publication-only row in the
   Phase 1 truth matrix. A truthful shipped `pending`/non-live statement does
   not need a cosmetic post-release rewrite.
3. Smoke the live PyPI install and hosted application through K's accepted
   post-publication commands; reference their evidence rather than republishing.
4. Update an identifier only when its truth-matrix row is explicitly
   publication-only, then rerun its deterministic link/citation/version
   contracts. A verification-only closeout with no repository diff is normal.
5. If a live identifier must be inserted into a shipped owner, stop. Create a
   patch candidate and return through the affected A1 source phase, H-final,
   and J instead of changing the accepted release after certification.
6. Confirm no shipped source, behavior, screenshot, recipe output, or generated
   Gallery asset changed during closeout.
7. Hand the exact citation/version state to K-Publication.

#### Exit gate

- Every route marked `live` has observed external evidence.
- Shipped citation/archive wording remains truthful; publication-only owners,
  if any changed, agree on the exact released version and DOI.
- Any publication-only diff and its affected contracts are reviewed; otherwise
  the no-diff verification is recorded.
- A1 may be marked `complete` only now.

## 8. Acceptance matrix

| Gate | Candidate synchronization | Exact H-final/J qualification | Post-release closeout |
| --- | --- | --- | --- |
| Branch/source | A0 drift verified; source revision recorded | J names same accepted revision | Final tag points to accepted revision |
| Version/schema | Source/docs contracts pass | Exact artifacts report same values | Live route/citation reports released values |
| Base install | Provisional local smoke | Supported OS/Python matrix on exact wheel | Live PyPI smoke |
| Export extra | Provisional recipe/export smoke | Exact-artifact matrix | No new build |
| CLI/Python recipes | All 38 pass against provisional wheel | Repeated against exact wheel | No shipped command edit; patch flow otherwise |
| GUI docs | Nightly check and visual review | Reproduction gate at candidate | No shipped documentation edit |
| Gallery | Generator/strict/browser gates | Packaged/offline exact-artifact gate | Hosted smoke only |
| Privacy | Frozen prose agrees with retained/cut behavior | J offline/network/privacy gates | Hosted consent smoke |
| Release notes | Candidate-complete with final truthful identifier strategy | Final-version candidate checked | Verify shipped wording; publication-only update or patch flow |
| Citation | Observed/reserved value or truthful non-live wording | Version concordance | Verify shipped wording; publication-only DOI update or patch flow |

## 9. Evidence ledger

Update this table only after the named evidence exists.

| Phase | Status | Behavior or decision completed | Evidence | Generated-artifact status | Deviations | Remaining risk |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | not started |  |  |  |  |  |
| 1 | not started |  |  |  |  |  |
| 2 | not started |  |  |  |  |  |
| 3 | not started |  |  |  |  |  |
| 4 | not started |  |  |  |  |  |
| 5 | not started |  |  |  |  |  |
| 6 | not started |  |  |  |  |  |
| 7 | not started |  |  |  |  |  |

Each evidence entry must contain:

```text
Phase:
Status:
Candidate revision:
Behavior or decision completed:
Evidence: exact command or CI run, observed result, and inspected artifact
Artifact filenames and SHA-256:
Generated-artifact status:
Deviations:
Remaining risk:
Invalidated by:
```

## 10. Stop conditions

Stop and reopen the owning gate when:

- branch ancestry differs structurally from A0;
- package/session/request/support decisions are still moving;
- H1 cannot build or inspect a provisional artifact;
- a built-wheel failure is reproducible only with an editable-install
  workaround;
- a manifest-declared recipe, long GUI scenario, strict Gallery check, or
  required browser flow fails;
- a screenshot or figure cannot be regenerated through its owner;
- a new input lacks required redistribution/provenance evidence;
- a generated Gallery/session diff cannot be explained by the final writer;
- reference outputs or the protected social preview change unexpectedly;
- an unavailable external route or placeholder identifier would need to be
  described as live;
- a post-release live identifier would require changing a shipped file without
  creating and requalifying a patch candidate;
- integration into the accepted target, push, deployment, publication, archive,
  or submission would occur without explicit authority.

Do not mark the phase complete, weaken the gate, substitute a different
environment, or hide the failure in prose.

## 11. Completion and handoff

A1 has three truthful milestone states:

1. `candidate synchronization complete; H-final/J-RC pending` after Phases 0–5;
2. `exact candidate qualified; release actions pending` after Phase 6 and the
   accepted J gate;
3. `complete` only after Phase 7 live-identifier verification, any
   publication-only closeout, and its affected contracts pass. A repository
   diff is not required for this milestone.

The final handoff reports:

- source and final tag revisions;
- H-final wheel, sdist, and hosted-bundle hashes;
- exact commands/CI runs and results;
- changed public pages and generated artifacts;
- visual-review coverage;
- supported install/browser matrix;
- deferred P2 issues;
- live PyPI, hosted Web, Bioconda, and archive state;
- exact citation/version DOI;
- confirmation that no external action was taken without authorization.

Do not use an earlier candidate-synchronization pass as evidence for a later
source revision or rebuilt artifact.
