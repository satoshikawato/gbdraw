# PR 2 instruction: Gallery publication-ready vertical slice

- Parent plan: [Interactive Gallery session replay remediation plan](./INTERACTIVE_GALLERY_SESSION_REPLAY_REMEDIATION_PLAN_2026-08-22.md)
- Suggested branch: `fix/gallery-publication-ready`
- Base: fresh `origin/dev` after PR 1
- Dependency: PR 1's WSSV lazy-resource fix and focused evidence
- Merge property: this PR must be deployable by itself; do not merge a publication-path change without its refreshed artifacts and blocking gates
- Scope class: persisted-session contract, Gallery publication ownership, complete refresh, browser parity, CI, and hosted-build verification
- Schema decision: keep session version 40 and canonical request schema 5
- Follow-up: tutorial/media audit remains a separate documentation PR after this code-and-artifact slice

## 1. Outcome and boundaries

Deliver one complete vertical slice in which:

1. canonical Gallery versions 31-33/39 migrate to the current format while ordinary v27-30 CLI replay compatibility remains unchanged;
2. a valid current v40 session follows an identity path during admission/version migration, before explicit Gallery-only publication;
3. an invalid current v40 session is rejected instead of silently repaired;
4. every refresh-owned Gallery session passes through one JavaScript publication owner;
5. the first Generate request is equivalent to the committed render intent;
6. all 13 generator-owned sessions and all public Gallery artifacts are refreshed and committed in this PR;
7. all 11 public sessions pass an actual browser import-to-first-Generate parity gate; and
8. GitHub Pages and Cloudflare package only the reviewed checkout bytes.

The work fixes persisted active-state and resource-binding drift. It must not:

- bump the session or request schema;
- change renderer geometry or update `tests/reference_outputs/`;
- change normal current-session draft authority;
- add a current-v40 cleanup reader, fallback, or allowlist;
- add an argv-shaped generation contract or parse raw `cliInvocation.args` as publication semantics;
- add an external dependency or runtime network requirement;
- modify `examples/gbdraw_social_preview.png`; or
- defer generated artifacts and replay gates to another merge.

## 2. Required starting evidence

Create the branch according to repository policy and inspect the dirty tree before changing it:

```bash
git fetch origin
git switch --no-track -c fix/gallery-publication-ready origin/dev
git status --short
```

The following commands and paths exist at the time this instruction was written and are valid baseline checks:

```bash
node --test \
  tests/web/gallery-session-migration.test.mjs \
  tests/web/session-request.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs

python -m pytest \
  tests/test_session_io.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_gallery_session_semantics.py \
  tests/test_wssv_gallery_fastas.py \
  tests/test_web_packaging.py \
  -q

node tools/check-web-change-budget.mjs
```

Before implementation, capture a machine-readable inventory of:

- every current and historical session admission path;
- every function that changes Gallery `config`, `renderRequest`, `resources`, `webFiles`, results, cache, catalog, or geometry;
- the 11 public sessions and two refresh-only test fixtures;
- all generator-owned session/source/example/thumbnail/catalog targets and their hashes;
- the current GitHub Pages and Cloudflare generation/copy paths; and
- baseline wall time and peak RSS for the unfiltered refresh and the largest-session publication projection.

Do not treat the plan's expected defects as observed evidence. Record the actual failures and measurements from the selected base SHA.

## 3. Persisted-session and inventory contracts

### 3.1 Strict version matrix

| Input | Required behavior | Written output |
| --- | --- | --- |
| Gallery canonical v31-33 or v39 | Validate its evidenced historical contract, migrate once through the compatibility owner, then validate as current | v40/schema 5 only |
| Ordinary CLI-only v27-30 | Keep the existing CLI replay compatibility; these versions have no canonical request and do not enter Gallery publication | Existing CLI compatibility output |
| Valid v40/schema 5 | Validate and take an identity path during admission/version migration; cloning or serialization is allowed, cleanup is not | Admission emits a semantically identical v40/schema 5 document |
| Invalid v40/schema 5 | Reject at the owning Web/Python boundary with safe field paths | None |
| Unsupported or branch-only version | Reject | None |

For the two retired paths below, raw v40 input must fail in both Web and Python:

```text
config.adv.cli_circular_track_order
config.adv.cli_circular_track_slots
```

The current field `config.adv.circular_track_slots` remains valid.

Gallery historical publication is limited to versions 31-33 and 39 because they contain a canonical `renderRequest`. Versions 27-30 remain supported by the ordinary CLI-only replay boundary and are neither rejected globally nor forced through a publication flow that requires a committed request. Record the first-parent commit or release tag and a representative positive fixture for every retained historical path. Do not support development-only versions 34-38.

Current v40 bypasses `promoteGallerySessionToCurrent()` entirely. The publication coordinator calls the shared current-active-config validator before any current publication work; after historical promotion it calls the same validator on the produced current config. `promoteGallerySessionToCurrent()` therefore owns historical migration only and must not import the current validator or retain a current cleanup branch.

Identity applies only to admission and version migration. After a Gallery source has passed that boundary, the explicitly invoked Gallery publication coordinator may replace its active draft with writer-native state derived from the committed render intent. That is a generator-only publication transform, not silent current-version migration: it must preserve the committed request semantics, pass canonical equivalence before writing, and remain unreachable from ordinary import.

Move the complete current Web active-config allowlist with its validator into the shared DOM-free contract defined below. Python may reject the two known retired envelope paths without duplicating the complete JavaScript reactive-state schema. A current writer must fail if handed invalid v40 state; it must not sanitize it before validation.

### 3.2 One-time tobacco chloroplast bootstrap

The checked-in tobacco chloroplast v40 is branch-owned invalid input, not a compatibility contract. Repair it once inside the same outer transaction as the strict unfiltered refresh:

1. start an outer Gallery transaction that snapshots every in-scope tracked target, including the invalid chloroplast source;
2. make a temporary staged chloroplast copy outside the tracked artifact path and remove exactly the two retired keys there;
3. produce a structured before/after report proving that no request, decoded resource, Result SVG, feature catalog, annotation, slot, legend, or geometry semantics changed;
4. validate the staged copy as current v40;
5. supply that staged copy as a one-invocation source override without replacing the tracked source;
6. prepare, replay, finalize, and validate all 13 sessions and all public assets in staging;
7. replace the repaired chloroplast and every other target only in the transaction's final replacement phase, restoring the original snapshot on any caught exception; and
8. discard the bootstrap adapter or temporary script before review.

Do not hand-edit the tracked one-line JSON. Do not commit a chloroplast-specific repair helper, a generic current-v40 cleaner, a retained source override, or the invalid branch-owned document as a positive compatibility fixture. A synthetic invalid-v40 negative fixture is appropriate.

The exact bootstrap command is deliberately not prescribed because no such repository command currently exists. Record the actual one-time command and semantic report in the PR evidence; do not present an invented command as a supported tool. Add fault injection after staged bootstrap and before final replacement to prove that the invalid tracked original remains untouched on failure.

### 3.3 Ordinary sessions versus Gallery publication

Normal import preserves the current contract that a saved Result can coexist with an intentional, not-yet-generated active draft. Gallery publication is a generator-only operation that deliberately aligns active state with the committed figure before publishing it.

These are separate paths:

```text
ordinary import
  -> evidence-backed historical migration when required
  -> current validation/restore

Gallery generator
  -> the same admission contract
  -> explicit Gallery publication preparation and finalization
```

Do not add a `publicationMode` boolean to the normal importer and do not route arbitrary user sessions through the Gallery publication owner.

### 3.4 One authoritative public inventory

Keep `tools/prepare_interactive_gallery_assets.py::EXAMPLES` as the authoritative owner of the 11 public examples. Extend its existing `GallerySessionExample` data only when a narrowly justified structured publication recipe is required. Do not introduce another inventory file.

Remove `GALLERY_SESSION_FILES` or derive it mechanically from `EXAMPLES`; it must not remain a second maintained public-session list. Browser selection, generated target enumeration, artifact-manifest generation, and hosted-build verification consume `EXAMPLES` or a checked generated projection whose exact parity is tested.

Keep the refresh-only scope distinct in `TEST_INPUT_SESSION_FILES`, resolving under `tests/test_inputs/` and containing exactly:

```text
AP027280_comparison.gbdraw-session.json
BGC0000708-BGC0000713.gbdraw-session.json
```

The unfiltered refresh set is computed as `EXAMPLES` plus `TEST_INPUT_SESSION_FILES`; do not maintain a third combined tuple. `examples.json` is the generated public projection of `EXAMPLES` and must not become an independently edited membership owner. JavaScript/browser consumers may read that generated projection only when a focused test proves exact membership parity with `EXAMPLES`.

## 4. Semantic ownership after this PR

### 4.1 Current admission and dedicated JavaScript publication owner

Extract `validateCurrentWriterActiveConfig()` and the complete field/default inventory it uses from `config.js` into one DOM-free current-active-config contract module under `gbdraw/web/js/services/`. Normal import and the Gallery publication coordinator both call that module. Move, rather than copy, the pure field/default inventory needed to validate `form` and `adv`; add a parity test against the reactive defaults. Update the internal test/import consumers in this PR and do not retain a `config.js` re-export alias solely as a compatibility shortcut.

The shared contract must not import reactive `state`, `config.js`, `gallery-session-migration.js`, or the publication module. This gives current admission one acyclic, DOM-free owner. The publication coordinator validates current v40 before calling publication logic and validates the current output of any historical migration.

Add a focused module under `gbdraw/web/js/services/`, provisionally `gallery-session-publication.js`. It owns only Gallery publication preparation, finalization, and readiness validation. Keep it separate from `gallery-session-migration.js`, whose responsibility contracts to evidenced historical migration.

The publication module may expose preparation and finalization functions; exact export names are implementation-selected. It must not become a general session service, a normal importer, or a second canonical request builder.

### 4.2 Request construction and equivalence owner

`session-request.js` remains the only canonical Web request projection/build owner. Put the canonical-form/equivalence API in that owner or in an existing narrowly scoped request/resource codec that it owns. Do not put a field-by-field request comparator in the Gallery publication module.

Canonical comparison must:

- recursively include all request fields by default;
- preserve array order and normalize object key order;
- replace a resource reference with its semantic binding path plus `{kind, decodedPayloadSha256}`;
- produce a stable digest and a structured safe-path diff; and
- exclude only fields proven not to affect SVG, feature catalog, editor semantics, or cache identity.

The initial permitted exclusion categories are:

- resource ID spelling after kind and decoded payload hash substitution;
- `createdAt`, result index, and runtime telemetry;
- output prefix, overwrite, and artifact filename; and
- the output-format representation forced by the publication tool.

Every additional exclusion requires evidence and a focused negative test. Do not maintain a manual list of every render field; new fields must be compared automatically unless explicitly classified non-semantic.

### 4.3 Python orchestration boundary

`tools/refresh_gallery_sessions.py` owns process invocation, temporary paths, staging, exception rollback, filesystem replacement, and sequencing. It may invoke existing Python schema/asset validators, but it must not choose publication semantics.

Move semantic decisions out of Python, including:

- legend, circular-slot, and palette synchronization;
- interactive metadata policy;
- provenance/resource merge precedence;
- resource deduplication and binding selection;
- comparison-plan repair; and
- post-replay selection among prepared and CLI-produced config/request/resources.

Remove the superseded per-field helpers rather than adding another repair. Also resolve these existing paths explicitly:

- delete `_load_gallery_refresh_source()` and use fail-closed `load_session()`; invalid current Results/catalogs must not be erased and retried;
- remove the direct refresh-source call to `normalize_current_session_artifacts()`; the ordinary Python writer may retain its existing contract when producing the replay sidecar, but refresh must not normalize source input before strict admission;
- delete `_omit_regenerable_gallery_derived_cache()` and its filename-specific Vibrio branch; retain or omit derived cache only through one general, identity-validated JavaScript finalizer policy, and only if evidence shows omission is still required; and
- delete `_preserve_gallery_cli_invocation()` and `_with_interactive_svg_format()` when it becomes unused; the JavaScript finalizer preserves provenance without rewriting argv or render formats.

The JavaScript publication finalizer owns the one semantic merge after CLI replay. Python passes prepared and replayed documents to that finalizer, receives a validated publication document, stages it, and handles I/O only.

Keep the large-session bridge path-oriented and staged. Before semantic request projection or canonicalization, detach non-projection bulk such as `results`, render caches, feature catalogs, and geometry blobs from the object graph being cloned or hashed; retain validated staged references and reattach them only through finalization. Digest resource payloads incrementally and preserve lazy resource views where projection needs content. Do not add a second whole-session JSON clone/buffer in Python and Node or keep source, prepared, replayed, and final 285 MB-scale documents resident at once.

### 4.4 Existing SVG comparator

Use `tests/utils/svg_compare.py` as the single SVG semantic comparator. Extend or migrate that implementation to provide the structured digest/diff needed by Gallery tests, and route the existing Vibrio and new browser replay checks through it.

Do not add `gallery-svg-fingerprint.cjs` or another independent JavaScript SVG normalizer. Playwright may write initial/generated SVGs to its output directory and invoke the shared Python comparator, as the current Vibrio test already does. Curated per-example anchors remain independent evidence against both SVGs drifting together.

Converge `tests/web/contracts/session-regenerate-intent.playwright.spec.js::svgEquivalence()` as well. Its semantic-equivalence decision must call the shared Python comparator or be removed. Exact-byte, raster, and direct-edit lifecycle probes may remain when their tests need those distinct signals, but they must be named and asserted as such and must not form another semantic-equivalence OR/fallback policy.

### 4.5 Authoritative architecture guidance

Update the module-ownership table in `gbdraw/web/CLAUDE.md` in the same PR. It must name the DOM-free current-active-config contract, historical migration owner, Gallery publication coordinator, canonical request/equivalence owner, and Python orchestration boundary without presenting overlapping owners.

Run the deterministic Web architecture checker against the proposed import graph. Update `tools/web-change-policy.json` only for privileged modules/importers the checker actually classifies as new, using complete and minimal importer sets; do not add blanket or speculative exceptions. Remove superseded registry paths in the same change and include the exact policy delta in the architecture review packet.

## 5. Publication and refresh flow

### 5.1 Admission and publication preparation

For each entry in the computed refresh set, the JavaScript bridge must:

1. apply the strict version matrix;
2. project the committed `renderRequest`, canonical resources, and current bindings through `projectCanonicalSessionRequest()`;
3. construct writer-native config, files, annotations, tracks, rules, layout, and resource views without treating stale saved config as render authority;
4. derive a deterministic current Web comparison plan from committed topology and canonical analysis/resource metadata;
5. rebuild the next request through `buildCanonicalRenderRequest()`;
6. compare its canonical form with the committed request; and
7. emit a current document only if equivalence succeeds.

For current v40 and Gallery publication, `cliInvocation` is opaque provenance. Preserve it without parsing `cliInvocation.args`, using it to select a plan, or making it a fallback authority. An evidence-backed pre-v40 migrator may retain existing argv interpretation only where that historical version's persisted contract and positive fixture prove it is required. Do not expand that legacy interpretation, use it for current v40, or treat it as publication authority.

Derive comparison plans deterministically:

- no committed comparisons becomes `none`;
- a complete adjacent edge set in canonical record order becomes an adjacent plan when canonical analysis metadata identifies the current generated pipeline;
- otherwise use an exact selected plan when current Web controls can represent every committed endpoint and resource binding; and
- fail closed when the committed topology cannot be represented losslessly.

If committed topology and canonical metadata genuinely cannot disambiguate one public Gallery recipe, add one structured recipe to its `EXAMPLES` entry. Keep its schema limited to the missing Web-plan facts, validate it against the committed request, and add a negative conflict test. A recipe is an exceptional data declaration, not a CLI parser or a place to duplicate render options.

### 5.2 Canonical replay and JavaScript finalization

The Python orchestrator writes the prepared session to staging and invokes the existing canonical CLI/session replay. It then invokes the JavaScript publication finalizer through staged paths or bounded semantic envelopes, not unconditional duplicate in-memory documents, with:

- the prepared publication document; and
- the replayed sidecar containing refreshed Results, cache, catalog, geometry, and request/resources.

The finalizer must:

- adopt replay-produced artifacts whose identity is valid;
- apply the one publication metadata policy;
- preserve `cliInvocation` as opaque provenance;
- merge and deduplicate resources by semantic binding and payload identity;
- reject collisions, missing bindings, stale config, or incompatible cache/catalog/geometry;
- rebuild the next request through the canonical request owner; and
- repeat committed/next-request equivalence before returning output.

Do not copy a saved committed comparison list into a rebuilt request merely to make parity pass. Do not copy `results` as proof of preparation equivalence; replay and final validation establish the published Result.

### 5.3 Complete generator transaction

The final supported unfiltered owner command remains:

```bash
python tools/refresh_gallery_sessions.py
```

For the one PR 2 migration run, a temporary outer bootstrap harness may invoke the same unfiltered pipeline with the staged chloroplast override. Its transaction must start before the override is created and must cover refresh, asset generation, and final replacement; it may not mutate the tracked source before entering the command transaction. Remove the harness before final review. Once the repaired artifact is checked in, the supported command above must reproduce the complete output without any bootstrap path.

The successful run covers:

- all 11 public sessions;
- both refresh-only test fixtures;
- all matching source SVGs, Gallery example SVGs, and thumbnails;
- `gbdraw/web/gallery/examples.json`; and
- the generated artifact hash manifest used by hosted-build verification.

All 13 sessions must be prepared, replayed, finalized, and validated before the session replacement phase. A later asset-generation exception must restore every in-scope pre-run file through `_gallery_file_transaction()`.

Describe this guarantee accurately: it is a staged transaction with exception rollback. Replacement consists of multiple filesystem operations, and a process/host crash can interrupt them. Do not call it crash-atomic or claim filesystem atomicity across all targets. Adding a crash-recovery journal is out of scope unless implementation evidence shows it is required for the stated remediation.

Partial `--session`, `--no-assets`, or `--skip-session-refresh` runs are useful diagnostics but are not completion evidence. Do not manually edit generator-owned outputs or preserve dirty in-scope copies over generator results.

## 6. Required tests

### 6.1 Version and ordinary-import contract

Update focused Web and Python suites to prove:

- an evidenced canonical v31-33 or v39 fixture migrates and writes neither retired field;
- ordinary v27-30 CLI replay compatibility remains unchanged and is not admitted to Gallery publication;
- valid v40 follows an identity path during admission/version migration and does not call current cleanup;
- explicit Gallery publication may align an admitted v40 active draft, while ordinary import remains unchanged and canonical committed semantics remain equivalent;
- v40 containing either retired field is rejected, and the safe error names the path;
- current `adv.circular_track_slots` remains valid;
- unsupported/branch-only versions fail;
- the repaired chloroplast imports through the real current path; and
- an ordinary current session still preserves an intentional divergent active draft.

Do not create a second Web retired-field list solely for test parity. Cross-language fixture/case parity is sufficient because the Web complete allowlist and Python narrow envelope guard have different responsibilities.

### 6.2 Publication owner and all-13 refresh

Add a DOM-free publication suite that reads the computed refresh set one entry at a time, or uses one child process per large entry. It must cover:

- all 13 inputs traverse the same preparation and finalization bridge;
- current and canonical historical v31-33/v39 inputs converge on current publication shape without introducing a current compatibility path;
- the input object is not mutated;
- rebuilt and committed request canonical forms match;
- every resource binding resolves to the expected kind and decoded SHA-256;
- checked-in public artifacts need no render-affecting correction;
- raw argv provenance has no effect on prepared semantics;
- a structured recipe is used only for an evidenced ambiguous case and conflicts fail;
- finalization cannot reintroduce stale config or select Python-owned semantics; and
- staged-bootstrap/bridge/replay/finalizer/asset failure triggers whole-scope exception rollback without first replacing the tracked chloroplast source.

Add explicit negative coverage that invalid current Result/catalog input fails without `_load_gallery_refresh_source()` sanitation, refresh never pre-normalizes current source artifacts, provenance argv is preserved rather than rewritten, and derived-cache retention has no filename-specific Vibrio branch.

Negative cases must include comparison topology and filters, record order/presentation/layout, arrow geometry, track side/order/params, legend/title, annotation state, color/priority payload, missing WSSV FASTA binding, unknown current `config.adv`, and resource collisions. These are inputs to one canonical comparator, not separate field-specific repair code.

### 6.3 Checked-in static publication gate

Against all 11 public checked-in sessions, before any publication transform:

1. validate current config;
2. import/adopt resource views;
3. build the first next request from raw active state;
4. compare it with committed intent through the canonical request equivalence owner;
5. verify every resource kind/hash binding; and
6. fail if publication preparation would change any render-affecting field.

This is a Gallery-only gate. It must not impose Gallery publication semantics on arbitrary user sessions.

### 6.4 Actual browser parity for all 11

For every public `EXAMPLES` entry:

```text
new browser context
  -> import exact checked-in session through the user-facing path
  -> assert saved preview visible and diagram Worker still lazy
  -> save initial SVG
  -> Generate once
  -> assert success, one Result replacement, and no unsafe fallback/network
  -> save generated SVG
  -> compare with tests/utils/svg_compare.py
  -> assert curated anchors
  -> Generate a second time where cache/lifecycle behavior matters
  -> close the context and release expanded data
```

Use one common test body selected by `EXAMPLES` metadata for the non-specialized examples. Circular and Linear may be separate serial projects. Preserve PR 1's final WSSV surfaces in place: `tests/web/contracts/wssv-full-generation.serial.spec.js`, `playwright.wssv.config.js`, and `test:web:wssv-generate`. Add the shared-comparator initial/generated parity assertion inside that spec and call the stable npm script from the isolated WSSV CI job; do not move, copy, rename, or delete its full-render body. Keep the existing specialized Vibrio serial spec as Vibrio's sole owner and invert its defect-preserving expectations.

The shared SVG comparator must attach initial/generated SVG, canonical digest data, structured diff, trace, screenshot, and safe stage summaries on failure. Do not include biological payloads in logs or artifacts beyond the already generated SVG under test.

### 6.5 Curated anchors

Retain focused independent assertions at least for:

| Example | Required evidence |
| --- | --- |
| AT skew | slot order `features,gc_content,gc_skew,a_skew_2,ticks`, AT colors/caption, labels, palette |
| Tobacco chloroplast | annotation set `plastome_regions`; `lsc`, `irb`, `ssc`, `ira`; three slots; resources; `upper_left` legend; saved Result |
| Majanivirus | 9 records, committed filters, arrow shaft ratio 0.5, 627 matches, three captions/colors |
| Vibrio | 11 records in order, five rows, 48 px gap, active adjacent LOSAT/collinear intent, observed nonzero comparison topology |
| WSSV | 20 ordered metadata-only FASTA views and conservation rings, cache-only replay, no executor/network fallback |

Observe and record Vibrio's intended comparison count after refresh; do not copy Majanivirus's count into that contract.

## 7. CI and performance gates

### 7.1 Fixed performance gate

On the same CI-equivalent runner, run baseline and final head three times each for at least:

- the complete unfiltered refresh; and
- publication projection of the largest session.

Run baseline and head from separate clean worktrees pinned to their recorded SHAs. Direct every measurement to disposable staging/output paths and discard its outputs; a timing run must not replace the review worktree's Gallery targets. Keep all 12 runs—six per operation—separate from the single unfiltered owner-command run that produces the reviewed generated-artifact diff.

For each operation, record the base/head SHA, environment, exact command or checked-in harness, three wall times, and three peak-RSS values. Compute wall time as the median of three runs and memory as the maximum peak RSS of three runs. Both final-head values must satisfy:

```text
head wall-time median <= 1.25 * baseline wall-time median
head peak-RSS maximum <= 1.25 * baseline peak-RSS maximum
```

The checked-in test or workflow must enforce both comparisons. This `1.25x` decision rule is fixed; implementation may select only the smallest reliable measurement harness. A PR-description measurement alone is not a gate. Do not shorten test-owned timeout assertions to make the gate pass.

The largest-session measurement must also prove the staged/detached design: no new whole-session duplicate clone or whole-payload hash buffer may be introduced merely to pass data between Python, Node, projection, and finalization.

### 7.2 Tiered, change-scoped CI

Implement a conservative, tested change classifier or equivalent workflow conditions:

| Tier | When | Gate |
| --- | --- | --- |
| Fast | Gallery/session/request/resource/hosted-build relevant changes | checked-in static all-11 publication gate |
| Refresh path | refresh generator, publication finalizer, inventory, or generated-artifact ownership changes | computed all-13 preparation/replay/finalization and rollback tests |
| Common browser | Web/session/rendering changes capable of affecting public diagrams | serial Circular and Linear first-Generate parity |
| Large isolated | the same relevant changes | stable WSSV script and specialized Vibrio script, each worker 1 and separately attributable |
| Pre-deploy | every hosted release | read-only static semantics, exact artifact manifest/hash verification, and a machine-enforced successful public all-11 browser result for the exact deployed SHA |

This remediation PR must run all tiers regardless of its own classifier result. Unrelated documentation or Python-only changes may skip the expensive browser replay only when the classifier proves they cannot affect the browser/session/render/packaging path. When classification is uncertain, run the gate.

Pre-deploy may reuse an already successful all-11 run only when the deploy workflow mechanically verifies the required job names and exact commit SHA. A manual PR note, branch-level success, or earlier head is insufficient. If exact-SHA evidence is absent, run all 11 before packaging or block deployment. Apply the same prerequisite to the Cloudflare production trigger rather than allowing its hosted builder to race ahead of the gate.

Do not duplicate WSSV or Vibrio full rendering in functional smoke. Upload large diagnostic artifacts only on failure.

## 8. Hosted builds are read-only consumers

### 8.1 Shared manifest/hash verifier

Generate a deterministic artifact hash manifest as part of the unfiltered generator transaction. Its deployed target set is derived from `EXAMPLES` and includes all public sessions, sources, examples, thumbnails, and `examples.json`.

Add one read-only verifier, with an implementation-selected name, that:

- rejects inventory/physical/catalog/manifest membership drift;
- checks checkout file hashes against the committed manifest;
- checks locally packaged file hashes against the same manifest;
- confirms local package bytes equal checkout bytes;
- for Cloudflare assets intentionally represented by the existing remote-asset manifest, verifies the immutable full-commit URL, target path, and checkout hash without regenerating content; and
- reports safe paths/hashes without rewriting files.

Both hosted builders must invoke this verifier. Do not add separate GitHub Pages and Cloudflare hash logic.

### 8.2 GitHub Pages

Remove `python tools/refresh_gallery_sessions.py` from `.github/workflows/deploy_web.yml`. Run read-only session/publication tests and the shared verifier before copy, then verify the upload tree against the same manifest after copy. Fail if any Gallery target changes during the build.

### 8.3 Cloudflare

Remove Gallery mutation from `tools/prepare_cloudflare_pages.py`, including the hosted-builder `--refresh-gallery` path, its loader, and the `refresh_gallery_sessions` parameter. The explicit developer generator already exists as `python tools/refresh_gallery_sessions.py`; the Cloudflare packager does not need a second generation mode.

Update packaging tests to require read-only verification and to reject any hosted-build Gallery mutation. Cloudflare's existing remote-asset packaging may transform bundle topology where already required, but it must source verified checkout bytes and the verifier must account for the documented remote-manifest outcome without regenerating Gallery content.

## 9. Generated-artifact and visual review

Review production, tests, workflows, and generated outputs separately. For each of the 13 sessions, produce a semantic report containing:

- session/request version, mode, grouping, record order, and layout;
- active, next-request, and committed-request canonical digests;
- decoded resource kind/SHA-256 map and `webFiles` binding summary;
- Result SVG digest;
- cache/catalog/geometry schema and counts; and
- compressed/uncompressed size.

For public examples, compare session Result, source SVG, Gallery SVG, and thumbnail. At readable scale inspect at least AT skew, Majanivirus, Vibrio, WSSV, and tobacco chloroplast, plus a contact sheet of all 11 thumbnails. Record each as `unchanged`, `intended change`, or `unexpected/blocking`.

Most final SVG bytes may remain unchanged when active state is corrected to match an already correct committed Result. Do not force visual churn to demonstrate that generation ran. `tests/reference_outputs/` and `examples/gbdraw_social_preview.png` remain untouched.

## 10. Verification commands

Run the existing focused commands from section 2 first. After implementation, run the complete unfiltered generator exactly once for reviewed artifact production:

```bash
python tools/refresh_gallery_sessions.py
```

The following command shape uses planned names that do not all exist on the base branch. Treat them as provisional, update this section or the PR evidence to the actual checked-in names, and do not claim them as verification until they run successfully:

```bash
node --test \
  tests/web/gallery-session-publication.test.mjs \
  tests/web/gallery-session-migration.test.mjs \
  tests/web/session-request.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs

python -m pytest \
  tests/test_session_io.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_gallery_session_semantics.py \
  tests/test_wssv_gallery_fastas.py \
  tests/test_web_packaging.py \
  -q

python tools/prepare_browser_wheel.py
npm run test:web:gallery-replay -- --project=gallery-circular
npm run test:web:gallery-replay -- --project=gallery-linear
npm run test:web:wssv-generate
npm run test:web:vibrio-generate
```

Run broader existing gates after the focused and browser evidence is green:

```bash
node --test tests/web/*.test.mjs
node tools/check-web-change-budget.mjs
node --test tests/web/architecture-contracts.test.mjs
npm run test:web:functional-smoke
pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
```

The output comparison must pass without regenerating reference outputs. Allow at least 30 minutes and monitor long browser/full-suite runs incrementally.

## 11. Acceptance criteria

### Session and architecture

- [ ] Canonical v31-33/v39 Gallery input migrates from an evidenced contract; ordinary CLI-only v27-30 compatibility remains unchanged and outside publication.
- [ ] Valid v40 validates and follows an identity path at admission/version migration; only the explicit Gallery publication stage may rewrite its active draft.
- [ ] Invalid v40 is rejected by the owning Web/Python boundaries; no current cleanup path remains.
- [ ] Tobacco chloroplast was repaired by a discarded one-time bootstrap and the tracked artifact is valid v40.
- [ ] `gallery-session-migration.js` owns historical migration only.
- [ ] `gbdraw/web/CLAUDE.md` records the final non-overlapping owners, and any `tools/web-change-policy.json` delta is minimal, complete, and checker-backed.
- [ ] One DOM-free current-active-config contract validates both normal current import and publication without a dependency cycle or duplicated field inventory.
- [ ] The dedicated publication module owns Gallery preparation/finalization only.
- [ ] `session-request.js` or its owned codec is the sole request construction/canonical-equivalence owner.
- [ ] Python owns orchestration/I/O/exception rollback and makes no publication merge decision.
- [ ] `_load_gallery_refresh_source`, the refresh-source normalization call, `_omit_regenerable_gallery_derived_cache`, and `_preserve_gallery_cli_invocation` are removed from the refresh path with negative evidence.
- [ ] Current-v40 and Gallery publication treat raw `cliInvocation.args` as opaque provenance; any retained pre-v40 interpretation is unchanged and fixture-proven.
- [ ] `EXAMPLES` is the sole public-11 owner, `TEST_INPUT_SESSION_FILES` is the separate two-fixture owner, and no maintained combined list remains.
- [ ] Normal current-session draft divergence remains intact.

### Artifacts and replay

- [ ] The unfiltered generator refreshes all 13 sessions and all public assets in this PR.
- [ ] The one-time chloroplast bootstrap uses a staged override inside the outer transaction; pre-replacement fault injection leaves the invalid tracked original and all other targets unchanged.
- [ ] Exception rollback is tested; documentation does not claim crash atomicity.
- [ ] All 11 checked-in public sessions are already publication-ready without render-affecting correction.
- [ ] All 11 pass real browser import and first-Generate SVG parity through `tests/utils/svg_compare.py`.
- [ ] WSSV remains in PR 1's stable spec/config/script, Vibrio remains in its specialized spec, and neither has a duplicate functional-smoke render.
- [ ] Curated anchors for AT skew, chloroplast, Majanivirus, Vibrio, and WSSV pass.
- [ ] No unexpected LOSAT/cache-miss executor or network fallback occurs.
- [ ] Production/test/workflow/generated diffs and gzip semantic reports are reviewed separately.

### CI, deploy, and non-goals

- [ ] Baseline/head are measured three times in SHA-pinned clean worktrees with disposable outputs; wall-time median and peak-RSS maximum are each enforced at `<= 1.25x` for unfiltered refresh and largest-session publication.
- [ ] CI is conservative and change-scoped, with WSSV/Vibrio isolated.
- [ ] GitHub Pages performs no Gallery generation.
- [ ] Cloudflare has no hosted-builder Gallery refresh mode.
- [ ] Both production deploy paths mechanically require the exact deployed SHA's successful public all-11 browser gate or run it before packaging.
- [ ] Both builders use the same read-only inventory/manifest/hash verifier; local package bytes are checkout-identical and Cloudflare remote entries identify the same reviewed commit/path/hash.
- [ ] Session version remains 40 and request schema remains 5.
- [ ] No external dependency, runtime network path, renderer geometry, reference SVG, or social-preview change is introduced.
- [ ] The final PR head is deployable; no required artifact or blocking gate is deferred.

## 12. Architecture fitness-function evidence

### Capability: Gallery publication state

Complete scope:

```text
{
  Gallery source admission,
  committed-state projection,
  deterministic Web-plan reconstruction,
  publication metadata/resource policy,
  post-replay semantic finalization,
  readiness validation
}
```

Owners before:

```text
{
  config.js current active-config validation,
  gallery-session-migration.js legacy/current normalization and argv interpretation,
  session-request.js canonical request projection/build,
  session_io.py current Python envelope validation,
  refresh_gallery_sessions.py per-field config repairs,
  refresh_gallery_sessions.py metadata/provenance/resource merge decisions,
  refresh_gallery_sessions.py post-replay artifact merge
}
```

Owners after:

```text
{
  current-active-config contract current Web admission validation,
  gallery-session-migration.js evidenced historical migration,
  gallery-session-publication.js preparation/finalization/readiness policy,
  session-request.js or its owned codec canonical request construction/equivalence,
  session_io.py current Python envelope validation
}
```

Non-semantic orchestration after:

```text
{
  refresh_gallery_sessions.py process/staging/rollback/filesystem orchestration
}
```

These owners have disjoint contracts: admission validates current state, migration transforms only evidenced history (including only fixture-proven legacy argv interpretation), publication decides Gallery preparation/finalization, request code constructs and compares canonical intent, and Python validates its session envelope. None may independently repair another owner's current input.

Canonical paths before:

```text
{
  current -> Python per-field repair -> CLI replay -> Python semantic merge,
  historical -> Node promotion/argv interpretation -> Python repair -> CLI replay
}
```

Canonical path after:

```text
any inventory source -> strict admission -> JS publication preparation
  -> canonical request equivalence -> CLI replay -> JS publication finalization
  -> staged validation -> exception-rollback publication transaction
```

Compatibility paths after:

```text
{
  evidenced canonical Gallery v31-33/v39 migration,
  ordinary CLI-only v27-30 replay outside Gallery publication (unchanged)
}
```

Current-v40 cleanup compatibility paths after: `none`.

Expected ratchet result:

- Gallery publication owner excess decreases;
- separate current/historical refresh paths converge;
- compatibility burden does not increase;
- request and SVG comparison each have one owner;
- no new dependency cycle, current reader fallback, or accepted-violation expansion appears.

Prepare the complete before/after owner, path, compatibility, dependency, and privileged-importer evidence against the final head SHA. Run the deterministic architecture checker and prepare the maintainer review packet; the implementing agent must not self-approve the architecture change.

## 13. Evidence ledger

Keep this table truthful during implementation. `Pending` changes only after the named evidence exists.

| Work item | Status | Required evidence |
| --- | --- | --- |
| Strict version matrix and chloroplast bootstrap | Pending | historical positive/current negative tests, identity test, bootstrap semantic report |
| Inventory and publication owner convergence | Pending | `EXAMPLES`/fixture-set derivation proof, removed-owner list, all-13 DOM-free parity, architecture packet |
| Finalizer and unfiltered generator transaction | Pending | staged replay/finalization evidence, rollback fault tests, full generated diff |
| Public browser replay | Pending | all-11 initial/Generate comparison, curated anchors, failure-artifact contract |
| Performance and CI | Pending | three-run baseline/head table, enforced `1.25x` wall/RSS gate, and change-scope tests |
| Hosted-build read-only contract | Pending | GitHub/Cloudflare mutation removal and shared checkout/package hash evidence |

## 14. Diff review and handoff

Review these classes separately using the actual implementation-selected files:

```bash
git diff -- gbdraw/web/js/services/ gbdraw/session_io.py tools/
git diff -- gbdraw/web/CLAUDE.md tools/web-change-policy.json
git diff -- tests/ package.json playwright*.config.js
git diff -- .github/workflows/
git diff -- gbdraw/web/gallery/sessions/ tests/test_inputs/
git diff -- gbdraw/web/gallery/sources/ gbdraw/web/gallery/examples/
git diff -- gbdraw/web/gallery/thumbnails/ gbdraw/web/gallery/examples.json
git status --short
```

For gzip and binary outputs, supplement the raw diff with the semantic/hash reports from section 9.

Suggested commit title:

```text
Make Gallery sessions publication-ready
```

Suggested summary:

```text
- enforce strict current-session admission and repair the chloroplast artifact once
- converge Gallery preparation and replay finalization on canonical JavaScript owners
- commit the full refresh with all-session browser and read-only deploy gates
```
