# gbdraw v0.14.0 Release Roadmap and Codex Implementation Brief

**Status:** Proposed implementation and release plan  
**Audience:** Codex, the maintainer, and future implementation/review sessions  
**Target release:** `v0.14.0`  
**Current baseline:** repository development version `0.14.0b0`  
**Primary objective:** finish the current figure-construction product, release it through stable installation routes, archive it, update the preprint, and submit the journal manuscript  

> This document defines the final intended scope of `v0.14.0`. It supersedes earlier informal plans where they conflict with this roadmap. Do not expand the release with plasmid restriction analysis, common-feature detection, ORF prediction, chromosome-scale Genome Overview, or sequence editing.

---

## 1. Executive decision

`v0.14.0` is the publication baseline for the first gbdraw software paper.

The release should establish gbdraw as:

> A local-first, browser-native, reproducible environment for constructing publication-quality circular and linear genome diagrams, with CLI and Python interfaces for automation and integration.

The release is not intended to complete every plausible future use case. It should complete the existing small-genome figure workflow and then freeze scope.

The final release additions for `v0.14.0` are:

1. Per-record display-start anchoring: place an arbitrary source coordinate at 12 o'clock in Circular mode or at the left edge in Linear mode.
2. Manual feature placement that remains active independently of automatic overlap resolution, with one explicit lane in each supported direction.
3. A substantially cleaner LOSAT/comparison user interface based on progressive disclosure.
4. Responsive preview Phase 1: retain the proven-safe direct preview edits and add debounced, single-flight renderer updates for record display start, manual placement, and existing label-layout reconciliation through a path that cannot run analysis.
5. Preview navigation improvements, especially zoom-to-selection and cursor-centred zoom.
6. Annotation TSV download from the current Web editor draft, with a bounded Web/Python round-trip contract.
7. PyPI packaging and publishing readiness with exact built-artifact clean-install validation; external publication remains a separately authorised release action.
8. A mandatory hosted/local privacy boundary and, only if all consent, payload, and operations gates pass, minimal opt-in feature analytics on the hosted Web app.
9. Documentation consolidation, evidence-backed release qualification, archival, preprint revision, and journal submission through distinct gates.

After the feature, packaging-readiness, privacy, and candidate-documentation work
is complete, freeze scope and run J-RC. Only an accepted J-RC baseline may be
handed to an explicitly authorised release-candidate action. After the first RC,
allow only release-blocking and specifically approved P2 fixes, and requalify
every changed candidate.

---

## 2. Repository instructions and working rules

Before making changes:

1. Read the repository-root `AGENTS.md` and `CLAUDE.md`.
2. Read `gbdraw/web/CLAUDE.md` before changing the Web application.
3. Before A0, inspect the current `main` and `docs_renovation` branches and preserve unrelated worktree changes. After A0, compare later drift with its recorded integration baseline instead of assuming the relationship is unchanged.
4. Treat the committed canonical typed request and resources, current session format, and last atomically committed Result as the authoritative product contracts.
5. Keep proven-safe direct preview edits, but do not expand SVG-only state. Semantic intent remains history-backed and session/request-reproducible; geometry changes use the renderer.
6. Do not update reference outputs unless a reviewed rendering-contract change requires it.
7. Do not silently add compatibility aliases to fresh APIs. Persisted-session migration and fresh-input acceptance remain distinct concerns.
8. Do not replace this roadmap with a shorter plan before implementation. Track A0 and A1 independently in their status and completion ledgers.
9. A gate pass is bound to one commit, version, frozen scope/support manifest, exact commands or CI runs, and artifact hashes. A later source or artifact change invalidates affected evidence.
10. A gate pass records readiness only. It does not authorise tagging, pushing, publishing, deploying, archiving, submitting, or another external action.

---

## 3. Current product baseline

The following is a planning inventory for the `0.14.0b0` development line, not
release evidence. Some capabilities already exist and others are completed only
when their upstream work-package ledger and J traceability entry contain
observed evidence at the qualification commit:

- A small top-level Python API with `read_genbank`, `read_gff`, `draw_circular`, and `draw_linear`.
- First-party `Diagram` output with SVG/bytes/save methods.
- Mode-specific typed request and render planning.
- Canonical session/request serialization and migration boundaries.
- Circular and Linear genome diagrams.
- Browser GUI, CLI, and Python entry points using the same rendering semantics.
- Direct feature colour, label, visibility, and stroke editing.
- Undo/redo.
- Interactive SVG output.
- Region annotations: highlight, band, line, and bracket.
- Custom track stacks.
- LOSATN, TLOSATX, and LOSATP-based comparison workflows.
- Pairwise, similarity-group, and collinearity presentations.
- Session save/load and reproducibility helper files.
- Deterministic feature identities and SVG metadata.
- Executable documentation work on the `docs_renovation` branch.

`v0.14.0` should stabilise and complete this product rather than introducing a new biological-analysis domain.

---

## 4. Release-level product contracts

### 4.1 Last successfully committed Result remains authoritative

- The current preview, canonical request/resources, feature catalogue, editor state, and export actions describe one committed Result.
- Stage a renderer candidate against an immutable draft revision without replacing the live Result.
- Commit the candidate Result, canonical request/resources, feature catalogue, editor/selection metadata, and export authority together, only if that draft revision is still current.
- A failed, cancelled, stale, superseded, or rejected update leaves the committed bundle and viewport intact.
- Proven-safe direct edits remain history-backed edits to the current Result; they are not partial renderer commits.
- Exports use the committed result, not an unrendered draft.
- The UI must clearly indicate when the draft is newer than the committed preview.

### 4.2 Source coordinates remain authoritative

- Per-record display-start transforms in Circular and Linear mode must not rewrite biological coordinates.
- Feature popups, annotation tables, comparison metadata, and exported source-coordinate values remain in source coordinates.
- Reverse complement, crop, rotation, and rendering transforms must not destroy source feature identity.

### 4.3 Cross-interface parity

A public figure-setting change that affects rendering must be representable through:

- Web state and saved session;
- canonical typed request;
- Python API or typed API;
- CLI option or generated helper table where appropriate.

The Web app may offer convenience interactions, but saved semantic intent must not be trapped in transient DOM state.

Preview effect flags, debounce/scheduler state, revision tokens, and
transaction-only viewport snapshots are Web runtime coordination data. They are
not persisted figure intent. Existing user-owned saved navigation state is a
separate Web-session policy and never enters the canonical render request.
Persist semantic display-start and placement settings through the canonical
request/session and supported programmatic surfaces.

### 4.4 No hidden analysis reruns

- A user edit may apply a proven-safe direct patch, schedule a renderer rebuild, mark analysis evidence stale, or combine those effects. These effects are not mutually exclusive classes.
- Automatic renderer work starts from an immutable snapshot of the committed canonical request/resources and overlays only allowlisted render-safe intent from the triggering render revision. Clone mutable request/resource-index structures copy-on-write; do not deep-clone unchanged large resource payloads for each update.
- The render-only entry point must not invoke LOSAT, rediscover inputs, resolve comparisons, inspect analysis caches, or mutate analysis evidence.
- At most one shared diagram-worker job is active; renderer policy keeps one replaceable latest draft pending. A candidate may commit only when its render/lifecycle revisions and replayed direct-edit revision are current.
- Analysis-invalidating changes show stale state and require an explicit analysis action.

### 4.5 Local-first privacy

- Uploaded biological data remain local to the browser or local installation.
- Feature analytics, if retained by the frozen scope, may be enabled only by an
  explicit analytics-capable hosted build and a current explicit opt-in. It must
  not become a runtime dependency.
- `unknown`, `declined`, expired/invalid consent, and storage failure result in no
  analytics script, request, cookie, queue, or later backfill.
- Revoking consent blocks new dispatch, invalidates pending loading, clears the
  configured GA cookies, and returns the document to a no-tag state.
- Consent is a versioned browser preference separate from figure sessions,
  canonical requests, history, Undo/Redo, and Reset.
- Blocking external network access must not prevent the packaged application from reaching ready state or completing a core figure workflow.
- Source Web, local GUI, packaged browser-wheel Web, CLI, and Python are
  telemetry-incapable for this feature even if browser storage is forged.

---

# 5. Work packages

## Work package A0 — Lock the integration baseline and documentation contract

### Objective

Establish one auditable integration baseline before the remaining feature work proceeds. Record the branch relationship, current package and persisted-format contracts, documentation owners, executable-evidence inventory, and stable-input policy without pretending that the final `0.14` screenshots or release text can already be frozen.

### Current status at the split

Status: **in progress**.

The 2026-08-11 audit at `docs_renovation` commit `6f14e2c4` found that local `main` was an exact ancestor and there was no two-sided Git divergence. The development branch had not been integrated back into `main`, the shared worktree was still changing, and the last observed `origin/main` state required a fresh fetch before any final reconciliation.

The current `0.14.0b0`, session writer 40, and request schema 5 documentation agreed with the implementation at that audit point. The documentation architecture, 30 public tutorials, Gallery metadata, manifest-declared inputs, and executable capture/recipe harness were already strong. LF checkout policy and an `autocrlf` regression now protect tutorial fixture hashes. Known A0 gaps remain: README and installation-route wording must agree, unavailable PyPI publication must not be described as current, and the eventual A1 candidate baseline must be re-recorded after feature work stops.

### Required work

- Fetch and inspect the relevant refs before reconciliation. Record the exact baseline commit and whether the integration branch contains, is contained by, or has diverged from the intended base.
- Reconcile branch history deliberately. Do not overwrite newer implementation changes with stale documentation assumptions, and do not use a destructive rewrite of a shared branch as a shortcut.
- Record the current package version, session writer and supported readers, request schema and accepted versions, table/resource schemas, and Gallery session/request inventory.
- Keep the four public documentation owners—Tutorials, Technical documentation, FAQ, and Gallery—stable unless a distinct unanswered reader question requires an information-architecture change.
- Maintain the manifest of executable scenarios and stable, accession-pinned or repository-tracked inputs. Enforce LF checkout for checksum-bound text fixtures on Windows, WSL, Linux, and macOS Git worktrees.
- Correct current-baseline claims. In particular, distinguish a live installation route from a route planned for the final release and remove statements that present unreleased Work packages B–I as current behaviour.
- Hand every feature-dependent reference, recipe, screenshot, Gallery artifact, release-note item, and citation closeout to A1 rather than regenerating it prematurely.

### Acceptance criteria

- The baseline commit, ref relationship, dirty-tree boundary, package version, session version, request schema, and documentation inventory are recorded with commands or contract-test evidence.
- Current public documentation does not advertise an unavailable installation route or an unimplemented feature as currently usable.
- Manifest-declared inputs pass size and SHA-256 checks after a checkout simulated with `core.autocrlf=true`.
- Current schema and compatibility documentation agree with the code at the recorded baseline.
- A1 has an executable implementation plan with entry gates, owned files, regeneration commands, evidence requirements, invalidation rules, and a completion ledger.

### Scope stop

Do not regenerate the complete screenshot, figure, Gallery, or release-note set while rendering, persisted formats, packaging, or public controls are still changing. Do not use A0 to redesign the established documentation architecture.

---

## Work package A1 — Synchronize the final release documentation and evidence

Implementation plan:
[`WORK_PACKAGE_A1_FINAL_RELEASE_SYNCHRONIZATION_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_A1_FINAL_RELEASE_SYNCHRONIZATION_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

After Work packages B–I are feature-complete, synchronize public documentation,
executable recipes, screenshots, Gallery assets, compatibility statements,
release notes, and reproducibility inputs. J-RC performs the final exact-artifact
verification after H-final. After the final tag/archive actions, A1 performs a
small observed-identifier verification and may update only inventoried
publication-only owners before the K-Publication gate.

### Entry conditions

- A0 is accepted against a recorded integration baseline.
- Work packages B–I have frozen their release-facing behaviour, defaults, and public names.
- Gate 0 has recorded the feature-freeze revision used as the A1 comparison boundary.
- Work packages C and D have one coordinated final request/session writer; no branch-only partial persisted format is being documented.
- Work package H1 can produce provisional wheel, sdist, and canonical hosted
  bundle candidates. Work package I has frozen whether analytics is enabled and
  the matching privacy copy, consent, CSP/network, event, and operator contract;
  an analytics-disabled decision requires the corresponding absence contract.
- The candidate tree is clean apart from explicitly inventoried release-documentation work.

### Required work

- Re-fetch and verify the A0 ancestry/drift record. Integrate approved final documentation changes without rewriting the established implementation history; if structural branch divergence has appeared, reopen the A0 reconciliation gate. Record the exact candidate commit used for every A1 gate.
- Build and install a provisional candidate wheel in clean environments while synchronizing documentation. J-RC repeats the declared recipe and installation evidence against the exact H-final artifacts rather than treating the provisional build as final evidence.
- Regenerate CLI help, reference material, executable outputs, affected screenshots, Gallery sessions/assets/tutorial media, and release examples from the final code and declared inputs.
- Synchronize README, installation guide, quickstart, tutorials, Python API guide, CLI reference, session/request compatibility guide, FAQ, privacy/analytics documentation, release notes, Gallery metadata, and citation instructions.
- Describe Hosted Web, Bioconda, PyPI, and source installation according to their actual release state. Work package K owns publication, deployment, and final archive identifiers; A1 must not claim that those external actions have already happened.
- Record the final commands, results, reviewed artifacts, deviations, and remaining external dependencies in the A1 evidence ledger.
- Invalidate and rerun the affected A1 gates after any Work package J fix changes public behaviour, rendering, schemas, packaging, screenshots, or documented output.
- Before H-final, classify identifier owners as shipped or publication-only and finalize every shipped owner with an observed/reserved identifier or truthful non-live wording. After K records the live tag, package URLs, hosted build, and archive DOI, verify shipped owners, update only inventoried publication-only owners, rerun affected gates, and hand the result to K-Publication. If a shipped owner must change, create a patch candidate and return through H-final and J.

The post-release closeout is limited to publication documentation that is not
part of the accepted source distribution or hosted bundle. A change to shipped
source, package metadata, or a packaged release claim creates a patch candidate
and returns through H-final and Work package J.

### Acceptance criteria

- The A1 candidate-synchronization milestone (implementation-plan Phases 0–5) is complete before J-RC. J-RC does not wait for overall A1 completion, which occurs only after the post-publication verification phase.
- One recorded source commit identifies each provisional synchronization pass, and the J-RC evidence identifies the later exact H-final artifacts used for qualification.
- The base wheel installs and imports in a clean environment; the export extra passes its documented smoke path; the complete manifest-declared CLI and Python recipe sets pass against the built wheel.
- `docs/capture/run_all.py --tier nightly --check` passes against the candidate, including the long-running scenarios, and all changed public figures receive visual review.
- Gallery sessions, generated SVGs, thumbnails, tutorial metadata, and operation media pass their generator-owned regeneration and strict validation gates.
- The package version, session writer/readers, request schema/readers, CLI help, Python signatures, Gallery inventory, and compatibility documentation agree.
- “All documentation commands” means all manifest-declared executable recipes, installation smokes, generated help/signature contracts, and applicable static checks. It does not require executing illustrative placeholder fences as if they were complete recipes.
- README and installation documentation distinguish live, release-candidate, and planned Hosted Web, Bioconda, PyPI, and source routes without presenting an unavailable route as current.
- The A1 evidence ledger is complete, the documentation diff contains no unexplained generated artifacts, and no A1-invalidating Work package J fix remains unverified.
- Post-release verification records observed version/URL/hash/DOI evidence; any repository diff is limited to inventoried publication-only owners and passes the affected link, citation, and documentation contract checks before K-Publication.

### Scope stop

Do not reopen product design, create a fifth public documentation route, execute every illustrative code fence literally, mass-reacquire every legacy input without a provenance need, or run the full GUI capture on every ordinary pull request. Do not tag, publish, deploy, push, create a Zenodo record, or revise the external manuscript under A1; those actions remain in Work package K and require explicit authority.

---

## Work package B — Clean up the LOSAT and comparison UI

### Objective

Make sequence input and ordinary figure construction visually primary. Keep advanced comparison controls available without allowing them to dominate the initial viewport.

### Target information architecture

1. Input sequences.
2. Compact comparison summary.
3. Basic figure settings.
4. Generate or Update preview.
5. Advanced comparison and layout controls.

### Required UI behaviour

Provide a compact top-level comparison choice:

```text
Comparison

[ No comparison | Run LOSAT | Upload BLAST TSV ]

LOSATN · All adjacent records · Default filters    [Settings ▸]
```

Detailed settings should be collapsed initially.

Show only settings relevant to the active comparison mode:

- `No comparison`: no LOSAT-specific settings.
- `LOSATN`: nucleotide search settings only.
- `TLOSATX`: translated nucleotide settings and genetic-code fields only.
- `LOSATP Pairwise`: pairwise protein settings.
- `Similarity groups`: grouping-related options.
- `Collinear blocks`: collinearity-related options.

Move the following to an Advanced/Reproducibility subsection:

- thread count;
- cache details;
- raw result filenames;
- cache clear/download controls;
- job estimates;
- rarely changed search implementation details.

Pair-specific comparison editing should be a separate collapsed section such as:

```text
Selected pairs (3) ▸
```

### Responsive-layout requirement

At a typical desktop viewport such as 1280 × 720:

- the first sequence upload control must be visible without scrolling;
- adding a second sequence should remain discoverable;
- comparison settings must not push the primary input below the initial viewport.

### State contract

- Closing an accordion must not reset values.
- Mode switching must preserve inactive drafts without activating them silently.
- `No comparison` must not initialise LOSAT, inspect comparison files, or revive stale comparison artifacts.
- Session save/load must preserve comparison intent independently from disclosure state.

### Acceptance criteria

- Playwright tests verify the initial viewport.
- `No comparison` performs no LOSAT work.
- Only program-relevant controls are visible.
- Keyboard users can open and close disclosures.
- Session round-trip preserves values.

---

## Work package C — Per-record circular coordinate anchoring

Implementation plan:
[`WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

Allow each complete circular record to use an arbitrary source coordinate as its display start without changing the biological coordinate system. Circular diagrams place it at 12 o'clock; Linear diagrams place it at the left edge and wrap from the source end to source base 1.

### Required product contract

- Diagram mode, biological record topology, and display start are separate settings.
- Topology and display start are owned by each resolved record, not by its source file or Linear source card.
- One shared mode-neutral transform serves Circular and Linear renderers.
- Source sequence, source coordinates, feature identity, annotation data, depth data, and comparison data remain unchanged.
- The no-shift default preserves current output, including reverse-complemented records.
- Cropped records do not accept a modular display start in `v0.14.0`.
- Every coordinate-dependent track, label, annotation, comparison, and interactive binding uses the same transform.
- Web, CLI, typed Python, package-root Python, and saved sessions represent equivalent intent.

### Scope stop

Do not add sequence rewriting, biological re-origining, submission-file rewriting, circular crop semantics, or base-level reconstruction of gapped alignments in `v0.14.0`.

---

## Work package D — Manual feature placement for overlap resolution

Implementation plan:
[`WORK_PACKAGE_D_MANUAL_FEATURE_PLACEMENT_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_D_MANUAL_FEATURE_PLACEMENT_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

Allow users to place an arbitrary source feature on the main lane or one supported displaced lane. Manual placement remains effective whether automatic overlap resolution is on or off. Features without an override remain Auto.

This work has two independent layers:

1. manual placement constraints owned by a source feature identity;
2. automatic placement of only the remaining Auto features.

`resolve_overlaps` controls the second layer. It must not disable or discard the first.

### Minimum release scope

The release-blocking minimum is:

```text
Overlap placement

Auto
Keep on main lane
1 lane outward / above
1 lane inward / below
```

The directional choices are required when the resolved feature slot is bidirectional:

| Resolved feature-slot geometry | Required manual choices |
|---|---|
| Circular `lane_direction=split`, combined strands | Auto, Main, Outward lane 1, Inward lane 1 |
| Linear feature slot `side=overlay`, combined strands | Auto, Main, Above lane 1, Below lane 1 |
| Circular `lane_direction=inside` / `outside` | Auto, Main |
| Linear feature slot `side=above` / `below` | Auto, Main |
| Separate-strand layouts | Auto, Main, where Main means the nominal lane for that strand pool |

Built-in `middle` presets normally resolve to the bidirectional cases above, but custom track slots may override their geometry. Do not decide support from the preset name alone. Validate token and level structure in the decoder, validate direction support after track-slot resolution, and retain the last successful preview on error.

### Resolver and precedence contract

| Requested placement | `resolve_overlaps = false` | `resolve_overlaps = true` |
|---|---|---|
| Auto | Use the normal unresolved lane; overlap is allowed | Use the deterministic automatic allocator |
| Main | Pin the feature to its nominal main lane | Reserve its nominal main lane before placing Auto features |
| Outward/Inward or Above/Below lane 1 | Apply the fixed lane | Reserve the fixed lane before placing Auto features |

- Manual placement always wins over automatic ordering.
- Conflict and Main-lane semantics use the resolved physical strand pool, not a raw strand string. Undefined-strand features remain in the nominal pool chosen by the existing layout policy.
- `Auto` is represented by the absence of an override. Selecting Auto deletes the stored override and recomputes all derived lane assignments from fresh feature objects.
- Do not persist a lane selected by the automatic allocator.
- An explicit feature may overlap an Auto feature when `resolve_overlaps` is off. Two explicit features that request the same physical lane and exceed the active genomic-overlap tolerance are a conflicting request and must fail validation.
- A hidden, underlay, or source-known crop-excluded feature keeps its override as dormant intent. It does not reserve a foreground lane until it becomes a visible foreground feature.
- Requested placement and the resolved/current lane are different concepts. The UI must not present an automatically chosen lane as a manual override.

### Optional extension after the minimum

```text
2 lanes outward
2 lanes inward
```

Lane 2 and beyond require stable sparse-lane geometry and are not part of the release-blocking implementation. Arbitrary pixel dragging is not part of `v0.14.0`.

### Do not expose raw internals

Do not present signed `feature_track_id` values directly. Their interpretation depends on mode, strand separation, and layout policy.

Persist semantic placement intent, for example:

```text
Auto
Main
Outward lane 1
Inward lane 1
```

The persisted value must carry mode-appropriate semantic intent, not an integer lane ID. The renderer may derive a transient `feature_track_id` only after it has resolved the semantic lane and strand pool.

### Genomic overlap tolerance

Treat near-boundary overlap tolerance as a separate setting:

```text
feature_overlap_tolerance_bp = 0
```

- The value is a non-negative integer in base pairs. Boolean values are invalid.
- Two feature envelopes conflict only when their genomic intersection is greater than the configured tolerance.
- `0` preserves the current half-open interval behaviour: touching endpoints do not conflict, and every positive intersection does.
- Use the same predicate for automatic allocation and explicit-placement conflict validation, including Circular origin-spanning intervals.
- The value applies only to genomic feature-lane collision. It does not relax label, annotation, or track-stack pixel clearance.
- Keep the current multipart outer-envelope collision rule. A feature wholly inside a long intron may still overlap the full envelope; manual displacement is the release solution for that case.
- Do not add a percentage threshold or change automatic main-lane priority in this release.

### Placement algorithm

1. Resolve record keys and stable source feature identities before visibility or rendering filters remove candidates.
2. Convert semantic placement intent into a physical lane after resolving the effective feature slot (`lane_direction` for Circular and `side` for Linear), mode, and strand pool.
3. Place explicitly constrained foreground features first.
4. Detect impossible conflicts among explicitly constrained features with the shared tolerance predicate.
5. If `resolve_overlaps` is on, place Auto features around the fixed features using the existing deterministic order and nearest-lane policy.
6. If `resolve_overlaps` is off, leave Auto features on their normal unresolved lanes.
7. Recompute required lane reservation and canvas geometry from the final assignments.
8. Re-run label layout using the final physical lane side rather than biological strand or `resolve_overlaps` as a proxy.

A useful interaction is:

```text
Feature A: Main lane
Feature B: Auto
```

The automatic resolver must then displace B if necessary.

### Persistence and reproducibility

- Store overrides by `(recordKey, biologicalFeatureId)` and retain source-index disambiguation for duplicate biological hashes. Resolve the ID against the original source catalog, carry the resolved source index through the record-aligned render plan, and do not re-hash cropped or reverse-complemented coordinates in the renderer.
- Preserve overrides across regeneration, colour edits, label edits, Circular rotation, reverse complement, record reordering, and session reopen. A feature retained by a crop keeps its source identity; a source-known feature currently outside the crop keeps a dormant override that can reactivate when the crop changes.
- Include placement overrides in the canonical request/session.
- Store only non-Auto overrides in a deterministic sorted array. Do not persist Web map keys, SVG DOM IDs, or resolved lane IDs.
- Provide a typed Python representation and a CLI/Python helper TSV. Human-authored TSV rows may use record and feature selectors, but each row must resolve to exactly one source feature before canonical serialization.
- Work package C Phase 3 owns the next schema 6 / session 41 writer. Fold D's persisted fields into the same atomic writer update with C's record-display fields and migration predicates; never create a D-only or C-incomplete 6/41 format. If C's writer is already main- or release-backed before D is integrated, allocate the next version from compatibility evidence rather than redefining 6/41.

### Error handling

- If two explicitly fixed features collide in an impossible way, show a direct validation message.
- Do not silently override one user's placement with another.
- Reject duplicate overrides, identities absent from the selected original source catalog, ambiguous selectors, unsupported directional/layout combinations, invalid lane levels, and invalid tolerance values. Do not misclassify a source-known feature outside the current crop as stale.
- Do not clear the last successful preview if the draft placement is invalid.

### Acceptance tests

- Resolver off with one manually displaced feature and all Auto features left on nominal lanes.
- Resolver on with an automatically displaced feature pinned back to Main; the previous Auto main-lane feature moves instead.
- Multiple overlaps with explicit features placed before Auto features.
- Positive, negative, and undefined strands in separate- and combined-strand layouts.
- Circular combined-strand resolved `lane_direction=split`, including built-in and custom-slot cases, with manual outward and inward lane 1.
- Linear combined-strand resolved feature slot `side=overlay`, including built-in and custom-slot cases, with manual above and below lane 1.
- Direct errors for unsupported directional choices in one-sided or separate-strand layouts.
- Exact default parity at tolerance 0, endpoint touching, and 1 bp / 2 bp tolerance boundaries.
- Multipart/intron-envelope and Circular origin-spanning cases.
- Final feature geometry, labels, leader lines, adjacent track reservation, comparison corridors, and canvas bounds use the same resolved lanes.
- Rotation, reverse complement, record reordering, crop-retained identity replay, and dormant reactivation for a source-known crop-excluded feature.
- Canonical request, CLI, Python, Web, save/load/re-render, and multi-record parity.
- Undo/redo and Auto-as-override-deletion.
- Invalid candidate placement preserves the last successful preview and does not rerun LOSAT.
- With no placement overrides and tolerance 0, existing reference outputs remain unchanged.

### Scope stop

Do not cut manual placement while `resolve_overlaps` is off or the two lane-1 directions in supported bidirectional combined-strand feature slots. If scope must shrink, defer overlap tolerance and lane 2 or greater before reducing the manual placement contract. Do not implement arbitrary pixel dragging, part-aware collision, or a new automatic priority heuristic in `v0.14.0`.

---

## Work package E — Responsive preview Phase 1

### Objective

Show common figure edits without requiring repeated Generate actions. Automatic preview work must never rerun analysis or rebuild analysis resources from mutable draft inputs.

The executable design, phase ownership, and evidence gates are defined in [Responsive Preview Phase 1 implementation plan](RESPONSIVE_PREVIEW_PHASE_1_IMPLEMENTATION_PLAN_2026-08-11.md).

### Release-blocking Phase 1 scope

Phase 1 must:

- preserve the existing proven-safe direct edit paths for feature/palette fill, stroke, visibility, label text, and existing composition/canvas movement;
- connect per-record display start and manual feature placement to debounced render-only updates, including feature-overlap tolerance if Work package D retains it;
- move existing label-layout reconciliation to the same no-analysis render-only path;
- consolidate the current Palette instant preview and Auto Reflow behaviour under one preview policy and status model;
- preserve committed selection and the minimum usable viewport across an accepted automatic render.

Do not add a general annotation or legend-content DOM patcher in this phase. Font sizes, arrow geometry, track order/geometry, annotation lane placement, and scale/tick geometry may join the render-only allowlist only after the initial scope passes the same persistence, history, race, and performance gates.

### Composable update effects

Each supported edit declares any applicable combination of these effects:

```text
patchNow
needsRender
invalidatesAnalysis
```

The effects are independent. Label text, for example, can patch the mounted SVG immediately and also request a later geometry reconciliation. A field that has not been classified must fail closed: it does not enter automatic preview work and remains pending for an explicit Generate action.

Track full normalized drawing intent separately from the render-safe and analysis projections. A manual-only or unknown edit must show `Changes pending`. Accepting a later allowlisted render must not clear that divergence or overlay the unadmitted field; only explicit Generate adopts the complete draft.

#### Direct patch policy

- Fill and stroke colour may patch existing semantic SVG targets immediately.
- Stroke width may patch an existing target, but any admitted geometry reconciliation remains a separate `needsRender` effect.
- Visibility is presence-aware. Hiding an existing target may patch immediately; showing a target omitted by the renderer requires render-only fallback.
- Label text may patch an existing label immediately and then request label-layout reconciliation.
- Existing legend/composition/canvas movement may keep its current direct path only where its session and export owner already exists.
- Annotation colour, annotation label, legend creation/deletion/rename, and any missing SVG target remain outside the direct-patch allowlist until their semantic metadata and catalogue updates can be committed together.

Direct patches must update the established backing Result/editor owner used by export and session save. The live DOM must not become a second authority.

#### Render-only policy

The render-only path uses this fixed boundary:

```text
latest render-safe draft revision
    +
snapshot of committed canonical request + immutable resources
    ↓
overlay allowlisted render fields at the canonical request boundary
    ↓
wait 400 ms after the last change
    ↓
run the existing typed renderer through the diagram worker
    ↓
stage, validate, sanitize, and commit only if the revision is current
```

This path must call the renderer directly. It must not call the normal analysis orchestrator in a special mode, depend on an analysis cache hit, or introduce another renderer/request contract.

#### Analysis-stale policy

Input sequence changes, record endpoints, LOSAT program/task, genetic code affecting search, and raw-evidence search arguments invalidate analysis. Phase 1 pauses automatic geometry rendering while this state exists and displays:

```text
Inputs or analysis settings changed. Preview data are out of date.
[Generate and update]
```

Use `Run LOSAT and update` instead when the active workflow requires new LOSAT evidence. Proven-safe direct edits may still update the committed preview. Only the explicit Generate/LOSAT action can clear analysis-stale state.

### Scheduler and worker contract

- Use trailing-edge 400 ms debounce for render-only work.
- Keep at most one render active and one replaceable latest draft pending.
- Do not terminate and restart the diagram worker for ordinary superseded keystrokes.
- Discard an obsolete response and render the latest pending draft next.
- Manual Generate removes queued preview work and has priority for the next renderer slot.
- `Update layout` bypasses debounce but uses the same no-analysis coordinator. One click dispatches the latest render-safe draft once, duplicate clicks coalesce, and analysis-stale state dispatches nothing.
- Session import, Reset, mode change, result replacement, and disposal invalidate queued revisions through the same coordinator.
- A direct-only edit does not erase an in-flight geometry request. Replay snapshots carry a direct-edit revision; if it changes during candidate staging, discard that candidate and queue the latest render-safe draft once.
- Serialize renderer and feature-extraction work through one shared-worker lifecycle broker. Manual Cancel terminates only its owned active job, releases the slot promptly, and preserves unrelated queued work for worker recreation. Ordinary supersession does not terminate the worker.

### UI and status

Provide one control for automatic geometry updates:

```text
☑ Auto-update layout
```

Direct colour and other proven-safe patches remain immediate when this control is off. Treat the preference as unset until the user changes it or a session restores it. After the first successful manual Generate, turn it on only when it is still unset; preserve an explicit or restored off choice. Omit the persisted preference key while it remains unset, and keep the preference out of figure-edit Undo history. Session hydration itself must not trigger an automatic render.

When automatic layout updates are off and a render-safe draft is pending, show:

```text
[Update layout]
```

Use a non-blocking preview status near the Result:

```text
Up to date
Updating…
Changes pending
Analysis out of date
Update failed — previous preview retained
Live updates paused — Update layout
```

Do not use the full-page Generate processing overlay for automatic updates.

### Transaction and viewport contract

- Keep the committed preview visible while an automatic candidate is pending.
- Stage all fallible work before mutating live state.
- Precompute base-config capture and composition metadata from the staged sanitized SVG; do not defer those parsers to a post-assignment watcher. If a remaining mount/binding step fails, roll back the complete logical transaction.
- Commit Result list, selected Result, feature catalogue, extracted/biological feature state, render metadata, editor/selection bindings, and committed canonical request/resources as one logical transaction.
- Before commit, reproject every admitted direct-edit family onto the candidate through its established semantic owner. This includes palette/default and per-feature fill, feature stroke, visibility, label text/visibility, and supported legend/composition/canvas transforms. A render-only update must not erase an earlier direct edit.
- Restore selected objects by stable semantic identity and clear only identities absent from the candidate catalogue.
- Preserve zoom and pan on automatic commit. Rebind and reposition the selected object when it still exists; do not run the unconditional manual-Generate viewport reset.
- Failure, cancellation, stale completion, or supersession leaves the previous committed bundle, export authority, selection, and viewport unchanged.
- Draft mode or input changes must not clear editor/catalogue state that still describes the committed Result.

### History behaviour

The user edit owns the history item; automatic rendering is a derived effect and creates no second item. Coalesce a continuous interaction such as font size 8 → 9 → 10 → 11 into one history item. Undo/redo restores the user intent and requests at most one derived render for the restored state.

Before automatic updates are enabled, complete or verify the history/artifact separation and single SVG ingestion/commit boundary in [Web application performance remediation plan](WEBAPP_PERFORMANCE_REMEDIATION_IMPLEMENTATION_PLAN_2026-08-10.md). Auto-update must not clone or serialize a large SVG merely to record each input event.

### Performance fallback

Ignore the first render after worker initialisation when evaluating the live-update guardrail. Pause automatic layout updates for the current diagram after two consecutive accepted render-only updates each exceed 2 seconds, measured across worker rendering plus candidate staging/commit. Keep direct patches active and offer `Update layout` and `Resume auto-update`. The paused state is transient and is not written into the canonical request.

### Acceptance criteria

- Existing proven-safe direct edits remain visually immediate and survive export/session round-trip through their established semantic owners.
- Per-record display start, manual placement, retained feature-overlap tolerance, and admitted label-layout reconciliation update after one 400 ms debounce burst.
- With empty analysis caches and analysis executors configured to fail if called, a render-only update succeeds with zero LOSAT calls, cache probes, mutable input reads, or analysis-resource changes.
- One render is active at a time; the latest pending draft wins and stale responses cannot commit.
- A direct edit made during an in-flight geometry update survives in active and inactive Results; an outdated replay candidate cannot erase the edit or lose the pending geometry change.
- A failed update preserves the previous Result, catalogue, metadata, selection, viewport, export, and committed canonical request/resources.
- Automatic updates off perform no worker render; `Update layout` performs exactly one.
- Analysis-stale state performs no automatic geometry render and requires explicit Generate/LOSAT action.
- One continuous edit creates one history item; renderer completion creates none; undo/redo reconciles the restored state once.
- Circular single/grid/batch and Linear preserve Result topology and inactive batch outputs.
- Automatic commit preserves zoom/pan and stable selection.
- The accepted candidate uses the shared SVG sanitization/commit boundary once and introduces no reference-output change by itself.
- Render-only projection, canonical adoption, and current-history artifact refresh do not deep-clone unchanged large resource payloads for every accepted update.

---

## Work package F — Preview navigation improvements

### Objective

Make large Circular and Linear diagrams easier to inspect without attempting full nucleotide-level semantic zoom in this release.

### Required controls

- Zoom in.
- Zoom out.
- 100%.
- Fit diagram.
- Fit width.
- Reset pan.
- Current zoom percentage.
- Cursor-centred wheel zoom.
- Zoom to selected feature.
- Zoom to selected annotation.

### Behaviour requirements

- Work package E owns the minimum automatic-commit rule: preserve zoom and pan, avoid an unconditional viewport reset, and restore selection by stable identity when the object still exists. This package adds navigation controls and refined zoom behaviour on top of that baseline.
- Avoid the current experience where zooming around a fixed origin causes the inspected region to move away from the cursor.
- Ensure keyboard and button-based navigation works without wheel input.

### Scope stop

Do not implement base-letter rendering, complement tracks, translation tracks, or sequence editing in `v0.14.0`. Those belong to a future Sequence Inspector.

---

## Work package G — Annotation UX finish

Implementation plan:
[`WORK_PACKAGE_G_ANNOTATION_UX_FINISH_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_G_ANNOTATION_UX_FINISH_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

Finish the existing annotation authoring workflow by letting users download the
current editable annotation rows. Do not expand the annotation data model or the
row-oriented TSV schema.

### Required product contract

- Add a visible `Download TSV` action to the Region Annotations panel beside
  `Add set` and `Import TSV`.
- Export the normalized current editor draft, not the last committed Result or
  its canonical render request. The action must work before Generate and while
  runtime network access is blocked.
- Download one file named `annotations.tsv` with the established TSV MIME type.
- Reuse the existing annotation table encoder and text-download service. Do not
  add a second serializer or route the action through Generate, a worker, LOSAT,
  analysis history, or Result commit.
- Disable the action when the total annotation-row count is zero, including both
  zero sets and one or more empty sets, and expose a clear reason in the visible
  or accessible UI. Empty sets do not make a header-only file downloadable.
- Downloading is non-mutating: annotation state, history, Result, catalogue,
  committed request/resources, selection, and viewport remain unchanged.

Recommended top-row layout:

```text
[Add set] [Import TSV] [Download TSV]
```

### Bounded round-trip contract

“Without semantic loss” applies only to valid, non-empty annotation rows emitted
from normalized Web editor state and to the fields representable by the existing
table. For those encoded rows, Web re-import and Python
`read_annotation_table()` must preserve:

- set order, row order, set IDs, and annotation IDs;
- coordinate targets, including record selector, start/end, coordinate space,
  origin wrapping, and out-of-bounds policy;
- feature targets, including record selector, selectors, envelope, and circular
  path;
- record IDs and one-based record-index selectors such as `#2`;
- all four marks, labels, `null`/zero/non-zero lanes, and effective legend text;
- effective stroke, stroke width, dash array, line cap, fill, fill opacity,
  hatch fields, label colour, label font size, label orientation, label
  position, and label offset.

For this contract:

```text
effective style  = annotation.style ?? set.defaultStyle
effective legend = blankToNull(normalizeCell(annotation.legendLabel ?? set.legendLabel))
```

The encoder intentionally materializes those effective values into each row.
Re-import does not need to reconstruct whether a value was owned by the set or
the annotation, but it must produce the same effective value. In particular,
an explicit no-fill value must remain `null`; a blank exported fill cell must
not silently become the Web default `#94a3b8` when another style cell is set.

Text fields are compared after the table's existing cell normalization: replace
tabs/CR/LF with spaces, trim surrounding whitespace on import, and normalize a
blank optional legend to `null`. The existing row schema cannot represent an
empty set, set-level inheritance ownership, annotation/set metadata, or
byte-for-byte preservation of those normalized characters. Those properties
are outside the TSV round-trip contract rather than being described as lossless.

### Deferred selection shortcut

Do not add a toolbar or popup-level `Annotate selection` shortcut in `v0.14.0`.
Keep the existing per-set selected-feature action. The current selection helper
prefers `locus_tag`, then `gene`, then a bare `#N`; duplicate qualifiers can
expand into a grouped span, the fallback is not a stable feature selector, and
duplicate record IDs need catalogue-owned binding. A future shortcut first
needs stable feature-hash targeting, record-catalogue binding, a destination-set
interaction, and duplicate-identity browser tests. It does not widen Work
package E's annotation direct-preview allowlist.

### Acceptance criteria

- `Download TSV` is always visible in the Region Annotations top row and is
  keyboard accessible.
- Zero sets and empty sets show the reasoned disabled state; the first valid row
  enables the action.
- One activation creates exactly one `annotations.tsv` download from the current
  draft with the expected MIME type, header, column order, and terminal newline.
- Downloaded bytes re-import through the Web parser and the Python reader with
  the bounded row semantics above, including `#2`, all four marks, lane zero and
  non-zero lanes, all effective style/legend fields, and explicit no-fill.
- A real browser download can be selected by the existing Import TSV control and
  restores the same bounded semantic projection.
- Downloading before Generate and with runtime network blocked succeeds without
  rendering, analysis, history, canonical request, or committed Result changes.
- Current session/request schema versions, Python/CLI/API contracts, rendered
  SVGs, and reference outputs remain unchanged.

### Non-goals

Do not add:

- a new annotation model, TSV column, or table version;
- preservation of empty sets, set/item style-ownership structure, metadata, or
  byte-exact embedded tab/newline/surrounding-whitespace content through TSV;
- a toolbar/popup `Annotate selection` shortcut or new feature-identity path;
- a Python API, CLI, canonical request, session, worker, renderer, or SVG change;
- new annotation marks or target kinds;
- grouped selected-feature span as a new domain behaviour;
- drag-to-resize annotation endpoints;
- point callouts;
- full nested style editor;
- conversion into GenBank/GFF3 biological features.

---

## Work package H — PyPI release and packaging audit

Implementation plan:
[WORK_PACKAGE_H_PYPI_RELEASE_PACKAGING_AUDIT_IMPLEMENTATION_PLAN_2026-08-11.md](./WORK_PACKAGE_H_PYPI_RELEASE_PACKAGING_AUDIT_IMPLEMENTATION_PLAN_2026-08-11.md).
The roadmap owns the release contract; the plan owns execution details and
evidence status.

### Objective

Make `pip install gbdraw` a tested installation route and ensure that every
published file is one of the exact artifacts accepted by the release gates.

### Ownership and execution points

Work package H has two execution points:

- **H1 — packaging contract and automation:** fix package contents, metadata,
  version owners, artifact checks, clean-install tests, the canonical hosted
  bundle recipe, and protected Trusted Publishing workflows before feature
  freeze;
- **H-final — exact candidate artifacts:** after all shipped-file owners have
  stopped changing the candidate, build one wheel, one sdist, and one deployable
  Web bundle, audit them, record their hashes, and hand those immutable inputs to
  Work package J.

H-final runs for every candidate entering J-RC and again for the final-version
candidate entering J-Final. It produces evidence inputs; it does not certify the
candidate. Work package J owns certification. Work package K alone may create a
tag, invoke TestPyPI/PyPI publication, deploy, archive, or perform another
external action, and each action requires explicit authorisation.

### Distribution contract for `v0.14.0`

Keep one distribution and one import surface:

```text
Distribution name: gbdraw
Import name:       gbdraw
CLI command:       gbdraw
```

Do not introduce `gbdraw-core`, `gbdraw-gui`, or platform-specific wheels in
this release. The PyPI release consists of one `py3-none-any` wheel and one
sdist. The wheel retains the local/offline Web application and contains exactly
one prepared browser wheel.

Remove the Linux x86-64 native LOSAT executable from the universal wheel and
sdist. Native CLI analysis discovers an explicitly supplied LOSAT executable,
then `losat` or NCBI BLAST+ `blastp` on `PATH`, with the existing explicit error
when none is usable. LOSAT Wasm remains a browser asset. Package tests must
reject native executable members in a universal wheel and verify the declared
external fallback on every supported OS family.

Gallery files belong to the hosted Web bundle, not the PyPI wheel or sdist. The
hosted bundle audit covers its Gallery files and immutable remote-assets
manifest separately.

### Base and optional capabilities

The base wheel must support:

- package-root and typed Python APIs;
- Circular and Linear CLI;
- GenBank/DDBJ and GFF3 + FASTA input;
- static and interactive SVG;
- typed requests and session save/replay;
- annotation rendering;
- precomputed comparison rendering;
- local GUI startup from installed package resources with runtime network
  access blocked.

Keep CairoSVG-based non-SVG export in the existing extra:

```bash
python -m pip install "gbdraw[export]"
```

The base-install lane must not receive `cairosvg`, test tools, or other
development dependencies accidentally.

### Artifact inventories and budgets

Gate 0 freezes a required/forbidden inventory and compressed/unpacked size
budget for each artifact. H1 records fresh measurements rather than reusing
ignored files under `dist/`.

- The outer wheel requires Python code and data, entry-point metadata,
  license/notices, local Web assets, fonts, LOSAT Wasm, vendored Pyodide assets,
  and exactly one non-recursive browser wheel. It forbids native executables,
  Gallery files, tests, build caches, credentials, local absolute paths, debug
  output, stale nested wheels, and unrelated generated files.
- The sdist requires the source, build configuration, license/readme material,
  the one prepared browser wheel needed to build the local GUI, and only the
  tests/tools explicitly retained by its inventory. It forbids native LOSAT,
  Gallery files, `dist/`, `build/`, `*.egg-info`, caches, credentials, local
  absolute paths, and stale generated artifacts.
- The hosted bundle requires the prepared browser wheel whose bytes match the
  single nested member in the accepted outer wheel, plus Web runtime assets,
  deployment headers/build stamp, Gallery, and its remote-assets manifest. It
  forbids an outer CPython wheel substituted under the browser-wheel filename,
  an obsolete wheel, undeclared external runtime dependency, credential, cache,
  or maintainer-local path.

Record total and category sizes for Web vendor assets, Pyodide, nested wheels,
LOSAT Wasm, fonts, tutorial data, Gallery/remote assets, and all other material
categories large enough to affect the reviewed budget. An inventory or budget
change requires an explicit Gate 0 update and invalidates prior artifact
evidence.

### H-final build and audit

Run H-final in a newly created empty staging directory. The release tooling must
make the following order fail closed:

1. Record the candidate commit, intended tag, clean or declared worktree state,
   and version/support manifest.
2. Run `python tools/prepare_browser_wheel.py`, then
   `python tools/verify_gui_offline.py check-assets` before the outer build.
3. Build into the fresh staging directory and require exactly one wheel and one
   sdist. Do not use `dist/*.whl` or another wildcard that can select an older
   artifact.
4. Extract the nested browser wheel from the accepted outer wheel without
   changing its bytes, then pass that staged file to the repository-owned
   canonical bundle command. Do not use a separately rebuilt or mutable
   source-tree copy as the hosted-bundle input.
5. Run `twine check --strict`, `verify_gui_offline.py inspect-wheel`, the
   per-artifact inventory/budget audit, metadata/version checks, license/notice
   checks, and local-path/credential/cache/debug scans against exact paths.
6. Write an evidence manifest containing filenames, versions, intended index,
   compressed/unpacked and category sizes, SHA-256 hashes, source commit,
   intended tag, build commands, tool versions, and source-tree before/after
   state.

The candidate fails if the build creates an undeclared tracked-source change.
Only ignored generated browser-wheel/build outputs and explicitly declared
stamping inside the staged Web bundle are permitted. A later shipped-source
change creates a new candidate and invalidates the hashes and affected evidence.

An isolated sdist install may build a derived wheel as a qualification output.
Record it, but never substitute it for the certified release wheel or publish
it.

### Version and publication state machine

H1 must replace hard-coded beta filenames in tests with the shared version
owner where possible and enforce agreement among `pyproject.toml`, including
the final `Development Status` classifier, `gbdraw.version`, the embedded
browser-wheel name and Web config, `recipe/meta.yaml`, `CITATION.cff`, the
CHANGELOG heading and release-note filename/link, generated notices,
distribution metadata, CLI output, and the intended tag in the release
manifest. H validates source concordance; K retains ownership of any Bioconda
publication action.

Use these publication states:

| Candidate | Allowed H action | Allowed K action after the matching J gate |
|---|---|---|
| Development or beta branch | Build, audit, and dry-run only | None |
| `vX.Y.ZrcN` intended tag and matching prerelease metadata | Produce and audit an RC artifact set for J certification | Create the authorised tag, publish the accepted files to TestPyPI, and run staged smokes |
| `vX.Y.Z` intended tag and matching final metadata | Produce and audit a final artifact set for J certification; reject any prerelease version | Create the authorised final tag and publish the accepted files to PyPI |

Before K creates a tag, J checks the intended tag recorded in the manifest.
After tag creation, K verifies the actual tag, commit, version, and artifact
hashes. Published PyPI files are immutable: do not use `skip-existing` or rebuild
a filename. A bad published version follows the documented yank or patch-release
path.

### Trusted Publishing boundary

Separate build/test jobs from publisher jobs:

- build, audit, and matrix jobs have read-only repository permissions and no
  OIDC token;
- TestPyPI and PyPI jobs run only through explicit K invocation, use protected
  `testpypi` and `pypi` environments, and hold job-scoped `id-token: write`;
- publisher jobs download the accepted artifact set, verify its manifest and
  hashes, and never checkout or rebuild the package;
- branch pushes and pull requests cannot reach an upload step;
- RC metadata routes only to TestPyPI, while final non-prerelease metadata routes
  only to PyPI.

Static workflow tests must prove these rules without publishing. H1 may not use
a token-based fallback or perform a test upload as implementation evidence.

### Clean-install and support matrix

Use Python 3.10, 3.11, and 3.12 on Linux, macOS, and Windows unless Gate 0 freezes
a different tested matrix before `rc1`. Probe newer stable CPython versions
before that decision; passing versions may be added, while unverified versions
must not be described as release-tested. Do not add an artificial upper bound
solely to mirror CI. Package metadata, support documentation, and CI must state
the same tested boundary.

Every matrix job downloads the exact H-final wheel, starts outside the checkout
with `PYTHONPATH` unset, installs only the requested variant, runs `pip check`,
and records the installed version plus source artifact SHA-256. The full Python
suite remains a J gate on Linux. H runs a lightweight base-wheel matrix on every
declared OS/Python cell, the full H-AC4 capability smoke on a declared reference
lane, sdist qualification on minimum and maximum supported Linux Python, and
H-AC5 export/installed-GUI smokes on declared representative lanes.

### Acceptance criteria

- **H-AC1 — distribution contract:** the wheel and sdist match their frozen
  inventories and budgets; each contains exactly the expected browser wheel,
  neither contains Gallery or a native executable, and the external LOSAT/
  `blastp` fallback behaves as documented.
- **H-AC2 — protected workflow:** static tests prove candidate/tag/index
  routing, protected environments, job-scoped OIDC, accepted-artifact download
  and hash verification, no publisher rebuild, and no branch/PR upload path.
- **H-AC3 — exact artifact set:** each H-final run yields exactly one wheel, one
  sdist, one deployable Web bundle, and one complete evidence manifest. Strict
  metadata/content/license/budget checks pass and downstream jobs receive the
  same hashes.
- **H-AC4 — base artifact capability:** an artifact-only clean environment
  passes `pip check`, import/version/help, the documented package-root API,
  GenBank/DDBJ and GFF3 + FASTA, both CLI modes, static and interactive SVG,
  typed request/session replay, annotation rendering, and precomputed
  comparison rendering, with explicit output assertions.
- **H-AC5 — optional export and installed GUI:** a separate `[export]`
  environment produces representative PNG and PDF output, and the local GUI is
  served from installed package resources, reaches ready state, generates
  Circular and Linear results, replays a session, and exports while external
  network access is blocked.

After an authorised RC upload, K must download the exact version from TestPyPI
with dependencies resolved separately from production PyPI, verify its hash and
index metadata, and run H-AC4 outside the checkout. K also records the staged
hosted version/hash/privacy smoke. An unresolved failure blocks entry to
J-Final.

### Deferred hardening

Formal SBOM formats, byte-for-byte double builds, custom SLSA or signing
infrastructure, runtime lockfiles, additional interpreter/architecture matrices,
and a separate GUI distribution are not `v0.14.0` blockers unless Gate 0 records
an external requirement or the measured artifact cannot meet the frozen
inventory and size budget.

---

## Work package I — Feature analytics and privacy consent

Implementation plan:
[`WORK_PACKAGE_I_FEATURE_ANALYTICS_PRIVACY_CONSENT_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_I_FEATURE_ANALYTICS_PRIVACY_CONSENT_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

Make every shipped Web artifact privacy-correct, then measure a minimal set of
feature outcomes only if the hosted deployment, consent lifecycle, observed GA4
payload, user-facing notice, and operator configuration all pass their gates.
Analytics is conditional for `v0.14.0`; the local/offline and zero-contact
privacy contracts are not.

### Release decision and remediation boundary

Split the work into three gates:

- **I1 — mandatory privacy and artifact remediation:** remove the current eager
  GA script/config injection, make analytics disabled by default, establish one
  authoritative hosted artifact, and prove that source, local GUI, wheel, and
  analytics-disabled hosted artifacts cannot contact Google;
- **I2 — conditional minimal instrumentation:** add consent-controlled loading
  and the closed event registry only after I1 passes and the retained feature
  contracts have stable accepted-Result owners;
- **I3 — mandatory evidence and operations if I2 is retained:** approve the
  notice, configure the GA property, inspect real transport payloads, run the
  browser state matrix, and record a kill-switch procedure.

This is remediation, not greenfield instrumentation. At the planning baseline,
`tools/prepare_cloudflare_pages.py` injects `gtag.js` and `gtag('config')`
eagerly, while `.github/workflows/deploy_web.yml` copies `gbdraw/web` through a
different recipe. Neither behaviour is an acceptable release baseline.

If any enabled-analytics gate cannot be proven before Gate 0, ship the
analytics-disabled profile. Do not delay the release by weakening a privacy
gate or leave an untested dormant Google path.

### Artifact and deployment contract

| Artifact profile | Hosted capability | Required network behaviour |
| --- | --- | --- |
| Source Web tree, local GUI, and packaged browser wheel | No hosted analytics configuration, Google CSP origin, consent prompt, or effective loader path | Zero analytics script loads, requests, and cookies, including with forged `allowed` storage |
| Hosted, analytics disabled | No measurement ID or Google CSP origin; no analytics notice, event dispatch, or effective loader path | Zero analytics script loads, requests, and cookies |
| Hosted, analytics enabled | Passive same-origin build configuration only; GA may load dynamically after a valid `allowed` decision | Zero Google contact while `unknown` or `declined`; consented requests only while `allowed` |

H1 owns one canonical hosted-bundle recipe and deployment topology. The current
implementation plan selects `tools/prepare_cloudflare_pages.py`,
`dist/cloudflare-pages`, and `gbdraw/web/cloudflare-worker.js` as that path,
subject to a Phase 0 check of the live DNS/deployment owner. The GitHub Pages
workflow must then be disabled as a production publisher or assigned an
explicitly non-production purpose; it must not remain a second route to
`gbdraw.app`. If the live ownership check contradicts this selection, stop and
revise H1 once rather than maintaining two production recipes.

The enabled profile requires both an explicit build flag and a valid measurement
ID. The repository and ordinary build commands default to disabled. Hosted
packaging may widen only the built artifact's CSP to the smallest captured GA
host set; the source `gbdraw/web/index.html` keeps the same-origin/offline
contract. An external analytics host is an optional hosted capability, never a
runtime dependency.

Do not weaken the app's COOP/COEP cross-origin-isolation or worker security
headers to make GA load. If the consented tag cannot coexist with the required
hosted isolation contract, ship analytics disabled.

Before implementing the enabled loader, update `gbdraw/web/CLAUDE.md` with
maintainer approval so its authoritative same-origin rule explicitly preserves
source/local/disabled artifacts while documenting this one consent-gated,
hosted-build-only exception and its mandatory evidence. If that exception is not
approved, retain the disabled profile.

### Hosted-app consent model

Show the notice only in the analytics-enabled hosted profile. Use two equally
clear choices, not `Accept all / Only necessary / Deny all`:

```text
Privacy settings

With your permission, gbdraw uses Google Analytics to understand
which features are used and to improve the app.

[Allow usage analytics]
[Continue without analytics]

Privacy details
```

The application may initialize and process files locally while the decision is
`unknown`; it must not initialize analytics. All features remain available for
both choices. Use a non-modal notice that does not trap or block the application;
ignoring it leaves the state `unknown` and analytics off. Do not add a dismiss
action that silently means allow. A persistent `Privacy settings` action must
provide the same information and make withdrawal no harder than permission.

### Consent state and storage contract

Use one first-party, versioned browser record with these runtime states:

```text
unknown | allowed | declined
```

The stored shape is owned outside diagram state and contains only schema
version, notice/policy version, status, decision time, and expiry. Consent
expires after at most 180 days; a notice-version change, expired or malformed
record, invalid timestamp, unavailable storage, or storage exception resolves
to `unknown` and fails closed. A declined choice remains stored; analytics
cleanup must not erase the evidence needed to keep the next visit declined.

Consent and analytics deduplication state must never enter the canonical
request, `.gbdraw-session` data, history snapshots, Undo/Redo, saved navigation,
or `Reset Settings`. Session import cannot change consent. Per-visit event
deduplication stays in memory and is never transmitted as an identifier.

The state owner must handle reload and cross-tab `storage` events. It must also
invalidate a pending script-load generation when permission is withdrawn. Do
not create `dataLayer`, queue an event, buffer an event, or backfill an event
that occurred while consent was not `allowed`.

Because the planning baseline already loads GA eagerly, the first consent-aware
hosted visit must also clear known legacy GA cookies when the state is `unknown`
or `declined`, without loading the tag or contacting Google. The property owner
must record the disposition of data collected by the prior eager deployment.

### Required runtime behaviour

While `unknown` or `declined`:

- do not load the GA4 script;
- do not send analytics requests;
- do not set GA cookies;
- clear known legacy GA cookies from the prior eager deployment;
- do not initialize `dataLayer` or retain an event for later replay.

If allowed:

- persist the valid decision before any external request;
- initialize the bounded wrapper and load `gtag.js` exactly once, dynamically;
- enable analytics storage while keeping ad storage, ad user data, and ad
  personalization denied;
- disable automatic page views and send only the approved custom events;
- treat script blocking, timeout, and collector failure as non-fatal no-ops for
  every application feature.

If declined:

- persist the declined decision;
- keep the analytics wrapper as a no-op;
- do not contact Google Analytics;
- keep every application feature available.

On withdrawal from `allowed`, synchronously block new dispatch, set the GA
disable flag, invalidate any loader in flight, clear the configured GA cookies,
then overwrite or remove the stored `allowed` record and read it back. Reload
into the no-tag state only after verifying that storage no longer resolves to
`allowed`; a failed `declined` write is safe if verified removal leaves
`unknown`. If that verification fails, do not reload: keep the current document
hard-disabled, show a privacy-state error, and fail the enabled-profile gate. No
request may be initiated between the withdrawal action and unload, and the new
document must contain no tag, `dataLayer`, or GA cookie. If the real tag cannot
satisfy this boundary, analytics is disabled for the release.

Local GUI, CLI, Python, source-checkout Web, and packaged browser-wheel Web are
telemetry-incapable for this feature, not merely opted out by a default value.

### GA transport and property gate

Wrapper validation does not prove what GA4 transmits. Before enabling the
profile:

- call `config` with `send_page_view: false`,
  `allow_google_signals: false`, and
  `allow_ad_personalization_signals: false`;
- disable Enhanced Measurement, Google Signals, advertising features,
  cross-domain linking, user-provided data, and unnecessary data sharing in the
  GA property;
- never set `user_id`, user properties, or a page URL/title/referrer derived
  from application state; use fixed page context and a strict referrer policy;
- capture the actual script and collection requests, including pagehide/beacon
  traffic, and inspect URL, body, headers, cookies, and referrer;
- enumerate and disclose the provider-generated browser, device, cookie,
  network, and approximate-location metadata that remains unavoidable;
- register the retained custom parameters as event-scoped reporting dimensions,
  freeze one schema-version-filtered report/query owner, and reconcile its three
  KPI outputs against a synthetic truth table before using live results;
- record property ID/account owner, access control, retention, deletion,
  incident/kill-switch, prior eager-deployment data disposition, and
  configuration-review evidence.

Automated tests must stub/intercept the collector and must not send test data to
the live property. A separately authorised staging capture supplies the real-tag
transport evidence. If an unapproved automatic event, field, endpoint, or
cookie cannot be suppressed or explicitly accepted and disclosed, retain the
analytics-disabled build.

### User-facing copy

The following remains the engineering draft for the main notice:

```text
Privacy settings

With your permission, gbdraw uses Google Analytics to understand which
features are used and to improve the app.

Your uploaded genome sequences are processed locally in your browser.
We do not send your uploaded files or file names, genome sequences,
annotations, gene names, genomic coordinates, LOSAT input sequences,
or generated figures to Google Analytics.

All features remain available whether or not you allow analytics.

[Allow usage analytics]
[Continue without analytics]

Privacy details
```

`Privacy details` is release-blocking content for an enabled build. It must state
the controller/contact, provider and purpose, exact bounded custom events,
provider-generated metadata, consent and GA storage, retention, withdrawal and
cookie cleanup, recipients/transfers, and deletion/contact route. It must also
list the prohibited biological data below. Maintainer, legal/privacy, and
analytics-property owner sign-off is required; engineering evidence alone does
not make a legal conclusion.

### Closed event contract

The first release uses one fail-closed registry. Unknown event names, parameter
keys, missing required keys, extra keys, and values outside an enum are rejected
locally and not sent. The analytics API accepts only the registry's primitive
enum values. It never accepts an `Error`, reactive/canonical state, request,
resource, session, file, SVG, or free-text value.

Retained events and their complete custom-parameter keys are:

| Event | Allowed custom parameters | Emission contract |
| --- | --- | --- |
| `workflow_start` | `event_schema` | Once per in-memory analytics visit, when the first genuine render attempt dispatches |
| `feature_progress` | `event_schema`, `feature_id`, `stage` | Once per feature/stage/visit; safe eligibility is captured at genuine dispatch, while rendered/exported use the committed Result; `stage = eligible \| rendered \| exported` |
| `render_result` | `event_schema`, `status`, `trigger`, `error_code` | Once per genuine current render attempt; all values are closed enums and success uses `error_code = none` |
| `export_figure` | `event_schema`, `format` | Only after the application successfully constructs and hands off an SVG, interactive SVG, PNG, or PDF export of the committed Result; browser/OS disk completion is not claimed |

`event_schema` is `1`. Initial trigger values are
`analysis_generate | manual_render | automatic_render | session_replay` and
formats are `svg | interactive_svg | png | pdf`. Freeze the exact error-code enum
from control-flow stages; never derive it from user-visible error text.

Emit `render_result` success only after the current candidate commits atomically. Emit failure only for a genuine failure of the current candidate. Debounce coalescing, pending replacement, stale/superseded completion, and cancellation are not render attempts or failures and emit no result event. Direct patches are not renders. A bounded `trigger` value may distinguish explicit analysis generation, manual no-analysis rendering, and automatic no-analysis rendering.

At genuine renderer dispatch, derive a safe eligibility snapshot from validated,
allowlisted capability booleans so an eligible attempt that later fails still
enters the denominator. After atomic acceptance, derive the active/rendered
snapshot. Both store only feature IDs and booleans. Export reads the committed
active snapshot rather than the live draft. Count each feature/stage at most once
per analytics visit. An analytics visit is one top-level main-app document
during one uninterrupted `allowed` epoch; reload or re-allow starts a new visit,
and no visit ID is sent.

Defer `panel_open`, `configured`, session save/load, tutorial, label-editing,
legend-editing, and Gallery events. The Gallery has a separate entry document
and CSP and no truthful tutorial-completion transition. It remains outside the
v0.14 analytics surface.

### Stable feature identifiers

Initial retained candidates are:

```text
record_display_start
manual_feature_placement
region_annotations
custom_tracks
depth_tracks
```

Before I2, give every retained identifier one exact eligibility predicate,
active predicate, genuine-dispatch owner, accepted-Result owner, safe classifier
input, deduplication scope, and metric consumer. Omit an identifier rather than
approximate any of those contracts. `record_display_start` is mode-neutral
because Work package C applies in Circular and Linear mode.

Defer LOSAT program/presentation identifiers from the `v0.14.0` analytics
registry. Their truthful attempt boundary is the validated LOSAT-job dispatch,
which occurs before renderer dispatch, and LOSAT failure therefore cannot share
the renderer-owned denominator above without a separate outcome owner. Do not
measure only the post-analysis survivors or compare that biased ratio with the
retained renderer-owned features. This is an analytics scope cut only; it does
not remove the LOSAT product work in Work package B.

### Data that must never be sent

- biological sequences;
- uploaded file contents;
- filenames or filesystem paths;
- record IDs or accessions;
- species/organism names;
- gene, product, or locus-tag values;
- annotation labels;
- exact genomic coordinates;
- SVG/PDF/PNG content;
- raw LOSAT input or output;
- raw error text or tracebacks;
- session titles.

Do not send record-count or sequence-size buckets in `v0.14.0`. They are derived
biological metadata and are unnecessary for the retained decisions. The closed
event table above is normative; examples do not extend it.

### Product metrics

For each retained feature:

```text
Successful adoption
= feature_progress(rendered) visits
  / feature_progress(eligible) visits
```

```text
Export conversion
= feature_progress(exported) visits
  / feature_progress(rendered) visits
```

Overall render reliability is:

```text
Failure rate
= render_result(status=failure) events
  / all render_result events
```

Exclude debounce replacements, stale/superseded completions, and cancellation from both the failure numerator and render-attempt denominator.

Do not report per-feature failure rates until a safe active-feature snapshot is
explicitly attached to attempt outcomes. Set no adoption target until a stable
baseline exists. Treat analytics as a consenting-user sample, not a census of
all users. The adoption denominator covers consenting visits that reached a
genuine eligible render dispatch; preflight-only failures and visits with no
render dispatch remain outside it. Privacy/network violations have a target of
exactly zero, and collector availability is never a product-health dependency.

### Acceptance criteria

- One canonical hosted bundle and deployment owner is recorded; a second
  production recipe cannot bypass the selected profile.
- Source, local GUI, browser-wheel, and analytics-disabled hosted artifacts make
  zero Google attempts and set zero GA cookies even with forged consent storage.
- The enabled hosted bundle makes zero Google requests and sets zero GA cookies
  for unknown, decline, decline/reload, corrupt/old storage, storage exception,
  revoke, and revoke/reload flows; pre-existing GA cookies are removed and it
  does not backfill pre-consent events.
- Allow, allow/reload, decline-to-allow, double-action, pending-load withdrawal,
  cross-tab transitions, and partial withdrawal-storage failures follow the
  state contract; reload never reactivates an unverified old `allowed` record.
- The actual request and cookie payload after allow matches the closed registry
  and contains none of the sentinel filename, accession, organism, gene, title,
  coordinate, raw-error, sequence, or figure data.
- Consent remains absent from request/session/history/Undo/Reset, and session
  import cannot alter it.
- Eligibility follows the validated current dispatch, rendered progress follows
  the accepted committed Result, automatic-preview scheduler churn emits
  nothing, and completed export events use the committed Result.
- Collector blocking or timeout does not block ready state, Circular/Linear
  generation, save, reload, or export on desktop or narrow profiles.
- Privacy details and property settings have recorded owner approval and list
  both sent and unsent data.
- If any enabled branch fails, the release artifact uses the analytics-disabled
  profile and passes its complete-absence tests.

---

## Work package J — QA, compatibility, and release engineering

Implementation plan:
[`WORK_PACKAGE_J_QA_COMPATIBILITY_RELEASE_ENGINEERING_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_J_QA_COMPATIBILITY_RELEASE_ENGINEERING_IMPLEMENTATION_PLAN_2026-08-11.md)

### Objective

Freeze the agreed release scope and schemas, validate one exact candidate commit
and its built artifacts against declared compatibility and support matrices,
record auditable evidence, and promote only a passing baseline toward
`v0.14.0`.

Work package J prepares and certifies release candidates. It does not itself
authorise a tag, package upload, hosted deployment, archive, or external
publication. Those actions remain in Work package K.

### Responsibility boundary

Feature-owned tests must land with Work packages B–I. Do not defer their unit,
rendering, request/session, or focused browser coverage until J. J owns:

- the requirements-to-test traceability audit;
- cross-package and cross-surface acceptance;
- compatibility and public-contract certification;
- the declared Python, OS, and browser support matrices;
- exact built-artifact, offline, security, and reproducibility gates;
- release-candidate evidence and final handoff to K.

If a mandatory Work package B–I acceptance criterion has no named automated or
reviewed evidence, J is blocked; it does not replace the missing test with a
checkmark.

### Gate 0 — scope, schema, and support freeze

Before candidate testing, record a release-support manifest containing:

- candidate commit SHA and clean/declared worktree state;
- target package version and candidate stage;
- retained and cut conditional features, especially overlap tolerance, preview
  fallback tuning, and analytics, plus the annotation shortcut's deferred status;
- current and accepted session, canonical request, cache, editor metadata,
  table, and interactive-SVG schema versions;
- supported native Python versions and operating systems;
- supported browser engines, versions or policy, and tested desktop/mobile
  profiles;
- wheel, sdist, hosted bundle, documentation, session, and figure artifacts in
  scope;
- native LOSAT and external fallback policy;
- analytics/offline network decision;
- P2 waiver authority and release rollback/yank owner.

For every retained conditional feature, all implementation, persistence,
documentation, and replay tests are mandatory. For every cut feature, remove
the shipped surface and claim and add absence tests. Do not leave a dormant
unreleased schema field or hidden UI path.

### Auditable evidence

Maintain a gate table with, at minimum:

| Field | Required content |
|---|---|
| Requirement | Stable roadmap/work-package acceptance ID |
| Owner | Work package and production owner |
| Execution | CI job or exact command |
| Fixture | Input/session/table plus provenance and hash where applicable |
| Environment | Commit, artifact, Python, OS, browser, tool versions |
| Oracle | Exact semantic, structural, visual, privacy, or performance assertion |
| Evidence | Job URL, log, report, screenshot, trace, artifact name, and hash |
| Status | Pass, fail, or approved P2 waiver with issue and workaround |

Unexpected skips, unexpected xfails, missing matrix jobs, hidden Playwright
retries, or a manual rerun that alone turns red to green do not count as passing
release evidence. Keep diagnostics outside tracked reference-output directories.

### Test layers

#### Unit and typed-contract tests

- coordinate transform;
- source/display coordinate round-trip;
- manual placement precedence and validation with the resolver on and off;
- feature-overlap tolerance at 0, 1, and 2 bp, including origin-spanning intervals;
- preview effect composition and fail-closed gating;
- render-safe draft projection onto cloned committed canonical request/resources;
- debounce, one-active/one-latest-pending scheduling, and stale-revision rejection;
- proof that render-only execution cannot reach LOSAT, analysis caches, mutable input reads, or analysis-resource mutation;
- request decoding;
- session migration;
- public `gbdraw` and `gbdraw.api` exports, signatures, dataclass fields/defaults,
  aliases, and exception semantics;
- CLI help/default/error contract;
- packaging metadata;
- version concordance across package, Web, session/request documentation,
  citation metadata, recipe, wheel/sdist metadata, and the intended release tag
  recorded in the frozen manifest;
- annotation TSV effective-semantic parity across the Web encoder, Web parser,
  and Python reader, including record-index selectors and explicit no-fill;
- versioned consent-state/storage reduction, including invalid/expired data,
  storage exceptions, reload, and cross-tab transitions;
- fail-closed analytics event registry, safe feature classification,
  deduplication, loader-generation invalidation, fixed error-code mapping, and
  synthetic KPI truth-table reconciliation.

#### Rendering tests

- per-record display-start transforms across every Circular and Linear
  coordinate consumer, including Circular wrap/origin and the Linear left edge;
- Circular and Linear manual placement, multipart/intron envelopes, labels, and canvas reservation;
- unchanged default output with no overrides and tolerance 0;
- deterministic output under a named normalization policy;
- structural/semantic SVG comparisons for default-output contracts;
- representative Circular and Linear publication figures rendered at readable
  scale and visually reviewed from declared recipes;
- no reference-output update without a reviewed rendering-contract change.

#### Web tests

- LOSAT progressive disclosure;
- initial viewport;
- annotation TSV download visibility, zero-row disabled reason, current-draft
  filename/MIME/content, browser download-to-import round-trip, multi-record
  selector retention, and zero render/analysis/history/Result side effects;
- immediate behaviour of the direct-patch allowlist;
- automatic rendering for display start, placement, overlap tolerance, and admitted label reconciliation;
- auto-update on/off, manual update, single-flight/latest-wins, and manual Generate priority;
- stale-analysis and stale-preview indicators;
- failed/stale candidate transaction across Result, catalogue, editor/selection metadata, canonical artifacts, export, and viewport;
- undo/redo reconciliation without a renderer-owned history item;
- minimum zoom/pan and stable-selection preservation during automatic commit;
- manual placement save/load and undo/redo;
- invalid fixed-placement conflict with the last successful preview retained;
- zoom controls;
- consent separation from canonical request, session import/export, history,
  Undo/Redo, saved navigation, and Reset;
- unknown, decline, allow, reload, revoke, re-allow, invalid/old storage,
  storage-error, double-action, pending-loader, and cross-tab consent paths;
- zero request/cookie and no-queue/backfill assertions before consent, after
  decline, and after revoke, including pagehide/beacon observation;
- closed payload inspection with biological-data sentinels in URL, body,
  headers, referrer, cookies, and `dataLayer`;
- shared SVG sanitisation, CSP, malicious import/session rejection, and atomic
  failure without replacing live state;
- packaged-app cold startup and core workflow with external network blocked;
- repeated worker use, cancellation, object-URL cleanup, and no biological-data
  canary in network requests.

The general Node Playwright suite is a release gate, not an optional local
check. Run the complete supported Chromium project and the focused smoke matrix
for every additional browser declared supported. A core workflow failure in a
supported browser is not an "uncommon browser" P2.

#### End-to-end acceptance tests

At minimum:

1. Circular figure with labels, GC/skew, annotations, a non-default per-record
   display start, resolver-off manual inward/outward placement, and session replay.
2. Linear multi-record figure with uploaded comparison, resolver-on Main/Auto precedence, and manual above/below placement.
3. Linear LOSATP pairwise workflow.
4. Linear LOSATP Similarity-group and Collinear-block workflows; they may reuse
   one deterministic input/session set but each presentation needs an explicit
   oracle.
5. Session save, close, reopen, render, and export through Web, CLI, and Python
   where supported.
6. H-AC4 against the exact base wheel and exact sdist in clean environments
   outside the source tree, including the declared matrix/basic-smoke and
   reference-lane capability oracles.
7. H-AC5 `[export]` clean install followed by representative PNG/PDF output.
8. H-AC5 source-tree and wheel-extracted local Web startup, Circular/Linear
   generation, session replay, and export with all external network blocked and
   forged `allowed` storage; assert zero external attempts.
9. The exact hosted bundle's unknown, decline/reload, allow/reload,
   allow-to-revoke/reload, decline-to-allow, corrupt/old storage, storage-error,
   loader-race, and cross-tab flows when analytics is retained. Inspect every
   intercepted request and cookie with forbidden-data canaries. If analytics is
   cut, assert the complete absence of hosted configuration, Google CSP origins,
   consent UI, and effective loader/event paths.
10. The enabled hosted bundle with Google endpoints blocked may attempt only the
    consented analytics request and must still reach ready state, generate,
    save, reopen, and export on desktop and narrow profiles.

Each flow must name its fixture, exact assertions, output inventory, timeout,
and evidence. "Representative", "where supported", and "works" are not
sufficient release oracles without the frozen support manifest.

### Compatibility rules

- Finalise the current writer only after the Work package C and D persisted
  fields are complete. If schema 6/session 41 remain the agreed boundary, write
  them once with both contracts; never land or retain a partial C-only or D-only
  format.
- Do not bump the session version for branch-only intermediate formats.
- Do not create multiple unreleased schema versions during implementation.
- Keep old persisted-data migration separate from fresh API acceptance.
- Track session, canonical request, raw/derived cache, identity manifest, editor
  metadata, Web file bindings, annotation/placement tables, Gallery sessions,
  and interactive SVG metadata as separate compatibility namespaces.
- For every distinct supported migration branch, keep an authentic positive
  fixture from first-parent `main` or a release tag, record its source commit,
  original path, and decompressed hash, and assert its expected semantic
  projection. Constructing a current document and changing only its version
  number is not sufficient positive compatibility evidence.
- Test current-writer load/save/reload, supported migration, future and retired
  version rejection, malformed/truncated/resource-hash failure, and rollback
  without partial live-state mutation.
- Effect flags, debounce state, scheduler tokens, pending jobs, slow-pause state,
  and transaction-only viewport snapshots are ephemeral runtime state. They
  must not create canonical request fields or cause a session-version bump.
- Do not conflate transaction viewport snapshots with existing user-owned Web
  navigation state such as saved zoom or pan. Freeze and document the session
  policy for that UI state; it never becomes canonical figure intent.
- Persist only the established preview preference and semantic figure intent.
  Keep one documented compatibility wire key for the renamed automatic-layout
  preference and do not write superseded keys.

### Support and artifact matrix

- Run the complete Python/core suite on Linux for Python 3.10, 3.11, and 3.12.
- Install-test the exact built artifacts on Linux, macOS, and Windows for the
  frozen supported Python matrix.
- Test base installation separately from `[export]`; the base lane must not
  receive CairoSVG accidentally through development dependencies.
- Build wheel, sdist, and deployable Web bundle once per candidate. Downstream
  jobs download and test those exact files instead of rebuilding them.
- A wheel derived only to qualify the exact sdist is a recorded test output. It
  cannot replace the certified release wheel and must never be published.
- Inspect package/archive contents, entry points, license/notices, native
  executables, Wasm/browser assets, generated data, local paths, caches, debug
  files, credentials, and unpublished datasets.
- Record artifact filename, target, version, size, SHA-256, tested invocation,
  fixture, expected result, and support status.
- The production hosted deployment must identify the accepted release commit
  and bundle hash. A later unrelated `main` build is not evidence for the tagged
  release.
- Reproduction claims require checksum-pinned, redistributable inputs and an
  executable recipe or session from a clean archive of the candidate. Public
  figures with unavailable inputs must be removed from the release
  reproducibility claim or explicitly labelled non-reproducible.

### Release gates

#### J-RC — candidate gate

The candidate passes only when:

- Gate 0 scope, schema, compatibility, and support manifests are frozen;
- every retained Work package B–I criterion has named evidence;
- focused, full, rendering, Node, Playwright, offline, clean-install, and
  documentation/reproduction gates pass against the candidate SHA/artifacts;
- no P0/P1 issue remains open;
- every deferred P2 has an owner-approved issue, impact, workaround, and release
  note where user-visible;
- the evidence manifest and artifact hashes are complete.

An authorised `v0.14.0rcN` tag and TestPyPI upload may occur through Work package
K only after J-RC passes.

K then downloads the accepted RC wheel and sdist from TestPyPI outside the
checkout, resolves dependencies separately from production PyPI, verifies both
filenames and hashes plus the version and index metadata, and runs H-AC4 against
the wheel. K also records the staged hosted version/hash/privacy smoke. An
unresolved failure blocks J-Final entry.

#### J-Final — final software-release gate

Any change after an accepted release candidate requires a regression test and a
new candidate or a complete rerun of the affected and full release gates. The
final version-only/documentation commit is still a different candidate. Rerun
the affected A1 candidate-synchronization gates as needed, rerun H-final to
create the exact final artifacts, then pass
version concordance, artifact inspection, clean install, full tests, and
reproduction gates.

J-Final certifies the exact final commit and stored artifacts. Work package K
then tags that same commit, publishes/deploys only the accepted files or bundle,
performs post-publication smoke checks, and records final URLs and hashes. A
failure before publication returns to J; a failure after publication follows
the documented rollback, yank, or patch-release path.

### Release-blocking severity

#### P0 — stop candidate testing and release

- data loss;
- incorrect biological coordinates;
- corrupted session;
- security/privacy failure;
- any Google analytics script load, request, or cookie while consent is unknown,
  declined, expired/invalid, or revoked;
- withdrawal reloads or another tab reactivates an unverified stale `allowed`
  record;
- any analytics attempt from source Web, local GUI, packaged browser-wheel Web,
  CLI, or Python;
- biological content, a forbidden-data sentinel, or unapproved metadata in an
  analytics URL, body, header, referrer, cookie, or automatic event;
- session import, history, Undo/Redo, saved navigation, or Reset changing
  analytics consent;
- package cannot install/import;
- exported figure materially incorrect.

After publication, a P0 requires an immediate rollback/yank assessment.

#### P1 — block promotion to final release

- core workflow crash;
- automatic preview reaches LOSAT, analysis-cache mutation, or mutable input/resource reconstruction;
- concurrent renderer jobs run or an obsolete draft revision commits;
- renderer-owned Result, catalogue, metadata, selection, canonical artifacts, and export authority commit only partially;
- undo/redo produces a Result that does not match the restored user intent;
- automatic rendering unconditionally resets zoom, pan, or a still-valid selection;
- display-start transform inconsistency across coordinate consumers;
- manual placement or overlap tolerance cannot replay;
- last good preview lost after failure;
- a core workflow fails in a declared supported browser or operating system;
- a required upload, Generate, save, export, or privacy control is inaccessible;
- retained analytics emits an incorrect bounded outcome or produces an
  unapproved duplicate without violating the P0 privacy boundary;
- documented command does not work.

#### P2 — defer only by recorded approval

- minor layout defect;
- issue outside the declared browser/support matrix;
- genuinely non-blocking accessibility defect;
- cosmetic inconsistency.

A P2 deferral requires maintainer approval, a tracked issue, user impact,
workaround, and release-note decision. After `rc1`, merge only P0/P1 fixes and
specifically approved P2 fixes. Every such change receives a regression test and
reruns the declared candidate gates; no fix is promoted only because a local
rerun passed.

---

## Work package K — Release, archive, preprint, and journal submission

### Objective

Use one exact software baseline for the package release, archived artifact, preprint revision, and journal submission.

Work package K performs external actions only after the matching Work package J
gate has passed and the maintainer has explicitly authorised the tag, upload,
deployment, archive, or submission. Passing tests is not standing authority to
publish.

Work package J pre-populates the expected commit, intended tag, accepted
artifact/bundle hashes, and required post-action smokes in the K action ledger
defined by the Work package J implementation plan. Work package K owns each
action authority, observed URL/hash, workflow run, and smoke result. The
K-Publication gate cannot pass while an applicable observed field remains
pending.

### Release sequence

```text
complete retained Work packages B–I and H1 readiness
    ↓
scope / schema / support freeze
    ↓
A1 candidate synchronization (implementation-plan Phases 0–5)
    ↓
H-final exact candidate artifacts
    ↓
J-RC candidate gate
    ↓
authorised v0.14.0rcN tag + TestPyPI + staged smoke
    ↓
TestPyPI exact-download H-AC4 + staged version/hash/privacy gate
    ↓
release-blocking fixes → new candidate → rerun J-RC
    ↓
final-version and affected A1 candidate gates rerun as needed
    ↓
H-final exact final artifacts
    ↓
J-Final gate
    ↓
v0.14.0 tag on the accepted commit
    ↓
PyPI publication of the accepted artifacts
    ↓
hosted Web deployment of the accepted bundle
    ↓
post-publication version / hash / workflow smoke
    ↓
Bioconda update
    ↓
Zenodo archive / version DOI
    ↓
A1 observed-identifier verification / publication-only closeout
    ↓
K-Publication gate
    ↓
bioRxiv preprint v2
    ↓
journal submission
```

If observed publication state requires a change to README, package metadata,
release notes, CHANGELOG, `CITATION.cff`, ABOUT, or another shipped file, stop
the sequence and create a patch candidate. Repeat the affected A1 source gates,
H-final, and J before publishing replacement artifacts. Phase 7 does not edit a
certified source distribution in place.

The release manifest binds the package tag, commit, session/request versions,
wheel, sdist, hosted bundle, documentation recipes, and manuscript figures.
Do not rebuild a different artifact between certification and publication. If a
publishing service requires a rebuild, treat the output as a new artifact and
rerun its inspection and smoke gates before promotion.

Deploy the exact J-certified hosted bundle through the single H1-selected
production path. Before declaring deployment complete, compare its build stamp
and hash to the manifest and repeat the enabled-or-disabled privacy smoke against
the live origin. A second workflow must not rebuild or overwrite the accepted
bundle.

### Publication gate

After the final software release is live and verified:

- complete the K action ledger with the authorised actions, observed URLs and
  hashes, workflow runs, and post-action smoke results;
- archive the accepted tag and release artifacts and record the version DOI;
- rerun manuscript-figure recipes from the archived baseline or verify their
  recorded candidate evidence remains byte/semantically identical;
- update preprint, manuscript, supplement, citation, and installation/support
  tables to the archived version;
- verify every capability claim has release evidence and every limitation or
  P2 waiver is disclosed;
- submit only after archive identifiers and exact reproduction instructions are
  final.

### Preprint v2 positioning

Replace a narrow “genome diagram generator” framing with:

> gbdraw is a local-first, browser-native, reproducible environment for constructing publication-quality genome diagrams, with CLI and Python interfaces for automation and integration.

### Claims to emphasise

- ordinary users can complete a figure in the GUI;
- CLI and Python are for automation, batch use, and embedding rather than mandatory completion steps;
- direct semantic editing remains reproducible;
- region annotations, per-record display start, manual overlap placement, comparison results, and every retained conditional feature are saved as project state according to the frozen release-scope manifest;
- uploaded genomic data remain local;
- static and interactive outputs are available;
- the same rendering semantics are shared across Web, CLI, and Python;
- PyPI and Bioconda provide complementary installation routes.

### Manuscript figures

At minimum, include a compact multi-panel figure showing:

1. Web/CLI/Python/session architecture.
2. Publication-quality Circular figure.
3. Linear comparative figure with LOSAT-derived relationships.
4. Direct editing and reproducibility workflow.

Supplementary material should contain:

- task-based competitor comparison;
- installation/support matrix;
- exact reproduction commands;
- browser compatibility;
- representative sessions;
- limitations;
- archive identifiers.

### Journal order

1. `Bioinformatics` — Application Note.
2. `Bioinformatics Advances` — Application Note.
3. JOSS as a fallback focused on software quality and maintainability.

Do not wait for `v0.15.0` if the first journal rejects or desk-rejects the manuscript. Submit the same `v0.14.x` publication baseline to the next suitable journal.

### Release-line branch policy

```text
release/0.14
  before final: release-blocking fixes and approved P2 fixes only
  after final: bug fixes
  reviewer-requested compatibility fixes
  patch releases only

main
  v0.15.0 development
```

Create or designate the `release/0.14` line at Gate 0 or the first accepted RC,
not only after submission. Every pre-final change invalidates the affected J
evidence. The manuscript should cite the exact release and archive DOI. Update
the manuscript to a patch release only if a relevant material defect is fixed.

### AI-assisted development disclosure

Prepare a transparent statement for the cover letter and Methods/Acknowledgements, for example:

> AI-assisted coding agents were used during software implementation and documentation. All architectural decisions, code review, testing, validation, and scientific interpretation were performed by the author.

---

# 6. Dependency order

Use the following order to reduce rework:

1. **A0 — Integration-baseline and documentation-contract lock**
   - record the integration branch relation, current schemas, documentation owners, executable-evidence inventory, and stable-input policy;
   - correct current-baseline claims, including the distinction between live and planned installation routes;
   - do not regenerate all assets yet if rendering contracts will still change.
2. **H1 hosted-artifact topology and I1 mandatory privacy remediation**
   - designate one canonical hosted-bundle recipe and one production deployment
     owner, remove eager GA loading/configuration, and make analytics disabled by
     default;
   - start after A0 and run alongside B–G. I1 must not wait for feature-event
     semantics because the planning baseline already violates the target consent
     boundary;
   - H1 may also resolve universal-wheel/native LOSAT policy, version owners,
     clean-install topology, and protected publishing automation, but it does not
     claim final artifact evidence.
3. **B — LOSAT/comparison UI cleanup**
   - largely independent of rendering geometry;
   - improves the working interface for subsequent testing.
4. **C model/planner work and D identity/allocator work**
   - establish source/display transform and record-aligned provenance before placement overrides depend on them;
   - allow independent implementation where the owners do not overlap, but do not describe C as fully complete while its writer depends on D.
5. **Joint C/D format-integration boundary**
   - bind to source identity through a record-aligned request-plan-to-renderer boundary, fold persisted fields into C's complete schema/session writer, and test under rotation, crop, reverse complement, duplicate record IDs, and custom feature slots;
   - treat C Phase 3 and D Phase 6 as one format-integration boundary rather than landing either partial writer;
   - keep manual placement independent from the Auto resolver and treat overlap tolerance as a separate render-only setting.
6. **Web performance prerequisites, then E — Responsive preview Phase 1**
   - first complete or verify the Web performance plan's history/artifact separation and single SVG ingestion/commit boundary;
   - retain proven-safe direct patches, add composable effects, and route display start, placement, retained overlap tolerance, and existing label reconciliation through a debounced single-flight no-analysis path based on committed canonical request/resources;
   - commit the current candidate atomically and preserve the minimum viewport.
7. **F — Preview navigation**
   - add navigation controls and refined zoom behaviour on top of E's commit-time viewport-preservation baseline.
8. **G — Annotation TSV download**
   - lock the bounded Web/Python table contract and repair any codec drift before
     exposing the download action;
   - use only the current annotation editor draft and existing codec after state
     contracts are stable. The deferred selection shortcut introduces no
     feature-identity dependency.
9. **I2/I3 — Conditional instrumentation, transport proof, or explicit cut**
   - after accepted-Result and export owners are stable, freeze the minimal closed
     event/feature registry and implement it only if analytics remains retained;
   - approve the privacy copy/property settings and pass the real-payload and
     complete browser state matrix, or select the analytics-disabled profile.
10. **Gate 0 — Scope, schema, compatibility, and support freeze**
    - record retained/cut conditional features, the complete C/D writer, supported runtimes, browsers, artifact set, figure inventory, and waiver authority;
    - no further feature or schema expansion is allowed after this boundary.
11. **A1 candidate synchronization — Implementation-plan Phases 0–5**
    - lock one exact candidate commit and run the built-wheel, executable-recipe, final-capture, Gallery, compatibility, and release-document gates in the A1 implementation plan;
    - regenerate examples and screenshots against the final behaviour, and invalidate affected evidence after any later release-blocking fix.
12. **H-final — Exact candidate artifact build and audit**
    - build the wheel, sdist, and deployable Web bundle once; inspect, hash, and pass them unchanged to the release matrix.
13. **J-RC — Full QA and release-candidate certification**
    - run the requirements, compatibility, support, Playwright, offline, artifact-install, and reproducibility gates;
    - require the A1 candidate-synchronization milestone, not overall A1 completion;
    - any accepted fix creates a new candidate and invalidates the affected evidence.
14. **Final-version/A1 source transaction, H-final rerun, then J-Final**
    - apply the final version and synchronize its shipped documentation owners;
    - rerun H-final, then have J-Final certify the resulting exact artifacts
      through version concordance, full tests, clean-install, offline, and
      reproduction gates.
15. **K release actions — Authorised tag, publication, deployment, Bioconda, and archive**
    - publish only J-Final-certified artifacts and record observed URLs, hashes,
      hosted state, and the version DOI.
16. **A1 post-publication verification — Implementation-plan Phase 7**
    - verify shipped identifier wording and update only inventoried
      publication-only owners; a required shipped-file change follows the patch,
      H-final, and J path.
17. **K-Publication — Preprint and journal programme**
    - close the publication checklist, then perform only separately authorised
      preprint and submission actions.

H1/I1 begin early because they remove a current privacy and deployment-path
defect. I2 waits for stable feature and accepted-Result contracts. Neither H1 nor
I may claim final support or artifact evidence before feature freeze. A1 and
H-final must run after every shipped-file owner has stopped changing the
candidate.

---

# 7. Explicit non-goals for v0.14.0

The following must not be added to this release:

- restriction enzyme site analysis;
- plasmid-specific preset or digest simulation;
- common plasmid-feature database;
- FASTA-only ORF prediction;
- biological feature creation/editing in GenBank/GFF3;
- reference-guided annotation;
- TLOSATN/TBLASTN implementation;
- nucleotide sequence editing;
- base-letter semantic zoom;
- cloning or assembly design;
- GFA/assembly-graph rendering;
- chromosome-scale Genome Overview;
- minimap2/SyRI/MCScanX orchestration;
- MCP server as a release blocker;
- desktop packaging as a release blocker;
- cloud accounts or collaborative editing.

If implementation reveals a prerequisite that belongs to one of these domains, expose the minimum internal abstraction needed for the current feature and defer the user-facing domain expansion.

---

# 8. Scope-cut rules

If the schedule or complexity expands, cut in this order:

1. Defer `feature_overlap_tolerance_bp` and lane 2 or greater; retain Auto/Main plus lane-1 outward/inward or above/below in supported bidirectional combined-strand feature slots, including when automatic resolution is off.
2. Defer non-essential preview-navigation polish; retain cursor-centred zoom and zoom-to-feature if possible.
3. Defer automatic performance fallback tuning; retain the preview effect policy, no-analysis boundary, scheduler, stale-response rejection, atomic commit, history semantics, and minimum viewport preservation.
4. Reduce analytics to fewer feature identifiers and the minimum
   `workflow_start`, `feature_progress`, `render_result`, and `export_figure`
   registry. Do not add tutorial/session events, biological size/count buckets,
   or per-feature failure correlation to rescue a weak metric.
5. If the consent lifecycle, approved copy/property settings, real transport
   payload, or canonical hosted-artifact proof cannot pass, cut analytics
   completely and ship the disabled profile. Do not cut I1 privacy remediation
   or its absence/offline tests.

Before Gate 0 closes, convert every conditional item into a binary release-scope
decision. A retained item receives all implementation, compatibility,
documentation, and replay gates. A cut item is removed from shipped UI/API/
session output and release claims and receives an absence test. Do not use
"optional" to leave an untested dormant field or CSP/network path in the final
candidate.

Do **not** cut:

- source-coordinate correctness;
- annotation TSV download with its bounded Web/Python row-semantics contract and
  reasoned zero-row state;
- manual placement while automatic overlap resolution is off;
- lane-1 placement in both supported directions for bidirectional combined-strand feature slots;
- session replay;
- last-successful-render transaction behaviour;
- proven-safe direct preview edits and the deliberately narrow display-start/placement auto-render scope;
- the dedicated no-analysis boundary, single-flight/latest-wins scheduling, and stale-response rejection;
- history-derived rendering and minimum viewport preservation;
- clean PyPI installation;
- exact-artifact inspection, hashes, and supported-platform smoke tests;
- compatibility fixtures and current-writer replay;
- packaged Web offline startup, sanitisation, and forbidden-data network checks;
- zero Google contact in non-allowed states and telemetry incapability in local
  artifacts; either satisfy the full consent/transport contract or ship no
  analytics capability;
- release archive and reproducibility documentation.

---

# 9. Release gates and checklists

Every checked item must point to the Work package J evidence table. A checkbox
without a commit, command/job, environment, oracle, and result is not complete.

## Gate 0 — Scope, schema, and support freeze

- [ ] Candidate commit and worktree state are recorded.
- [ ] Retained/cut status is recorded for every conditional feature and matches code, tests, docs, and event identifiers.
- [ ] Current and accepted session/request/cache/metadata/table versions are recorded by namespace.
- [ ] The complete joint C/D writer version is evidence-allocated once; no branch-only partial version is supported or documented.
- [ ] Native Python/OS and browser support matrices are declared and match package metadata and documentation.
- [ ] The universal wheel/sdist exclude native LOSAT, and the declared external LOSAT/`blastp` fallback is reflected in package contents, tests, and documentation.
- [ ] One hosted-bundle recipe and production deployment owner are frozen;
      analytics is either retained under the complete consent/transport/operator
      contract or disabled with absence tests.
- [ ] Public/release/manuscript figure inventories name reproducible recipes, sessions, inputs, hashes, and redistribution status.
- [ ] P2 waiver approver and rollback/yank/patch owner are recorded.
- [ ] Every mandatory Work package B–I criterion maps to a named test or reviewed artifact.

## J-RC — Release-candidate qualification

### Evidence health

- [ ] No unexpected skip, xfail, missing matrix job, hidden retry-only pass, or unresolved flaky result is counted as evidence.
- [ ] Production, test, documentation, fixture, reference-output, generated-artifact, and release-automation diffs are reviewed separately.
- [ ] No P0 or P1 issue remains open; every deferred P2 has an approved disposition.

### Product and cross-surface behavior

- [ ] Sequence upload is visible in the initial desktop viewport.
- [ ] No-comparison mode performs no comparison work.
- [ ] Per-record display start works for all supported Circular and Linear coordinate consumers.
- [ ] Source coordinates remain correct in popups, annotations, tables, metadata, and outputs.
- [ ] Main and lane-1 directional placement work with automatic overlap resolution on and off and survive regeneration and session replay.
- [ ] If retained, overlap tolerance preserves default behavior at 0 and follows the documented 0/1/2-bp boundary, including origin-spanning intervals.
- [ ] Proven-safe direct preview edits remain immediate and use their established session/export owners.
- [ ] Admitted geometry edits use one debounced render-only update that cannot reach an analysis entry point, mutable input read, or analysis-cache mutation.
- [ ] At most one automatic render is active; latest pending wins and an obsolete response cannot commit.
- [ ] Failed, cancelled, stale, and superseded candidates preserve the complete previous committed bundle and viewport.
- [ ] Auto-update off performs no automatic worker render; manual update performs exactly one.
- [ ] Undo/redo restores user intent and renders at most once without a renderer-owned history item.
- [ ] Automatic commit preserves zoom/pan and every still-valid stable selection.
- [ ] Every retained navigation control has an explicit browser assertion.
- [ ] Annotation TSV export is available from the current editor draft; zero-row
      state is clearly disabled, and valid non-empty rows re-import through Web
      and Python with the bounded effective target/style/legend semantics defined
      in Work package G.

### Session, API, and compatibility

- [ ] The final writer contains the complete display-start and placement contract; no partial intermediate writer remains.
- [ ] One documented automatic-layout preference wire key is written; retired preview keys are absent from current output.
- [ ] Authentic provenance-backed fixtures cover every distinct supported migration branch.
- [ ] Supported sessions load and semantically replay through CLI, Python, and Web according to the documented boundary.
- [ ] Future, branch-only, malformed, truncated, and resource-inconsistent inputs fail closed without replacing live state.
- [ ] Web, CLI, package-root Python, and typed API describe equivalent display-start and placement intent.
- [ ] `gbdraw` and `gbdraw.api` public exports, signatures, dataclass defaults, aliases, and errors match the approved contract snapshot.
- [ ] CLI help, defaults, validation, and documented commands match the approved contract.
- [ ] Deterministic feature identity survives rotation, reverse complement, crop-retained transforms, dormant/reactivated crop state, duplicate record IDs, and duplicate-hash disambiguation.
- [ ] Runtime effect/scheduler/transaction state is absent from canonical figure intent; the separate saved-navigation policy is documented and tested.

### Web, privacy, security, and offline

- [ ] The complete Node unit and supported Playwright suites pass against the candidate.
- [ ] Generated, imported, and restored SVG content passes through the shared sanitization contract.
- [ ] CSP and packaged asset URLs satisfy the frozen per-profile runtime/network
      policy; source/local/disabled artifacts have no Google allowance.
- [ ] Source and wheel-extracted Web cold-start and complete generation, session
      replay, and export with external network blocked and forged `allowed`
      storage, with zero external analytics attempts.
- [ ] Repeated generation/cancellation terminates or reuses workers as designed and releases object URLs without cross-run state leakage.
- [ ] Consent state is versioned and fail-closed for missing, expired, corrupt,
      old-policy, and unavailable storage; it remains separate from
      request/session/history/Undo/Reset.
- [ ] If analytics is retained, no `dataLayer`, queue/backfill, script, cookie, or
      request occurs while unknown or declined; known legacy GA cookies are
      removed without Google contact.
- [ ] Withdrawal initiates no new request before unload, clears GA cookies,
      reloads only after verified non-allowed storage, and produces a new
      document with no tag, `dataLayer`, GA cookie, or request.
- [ ] If analytics is retained, allow/reload, decline/reload, re-allow,
      pending-load withdrawal, cross-tab withdrawal, and blocked-collector flows
      pass against the exact hosted artifact.
- [ ] If analytics is retained, actual automatic and custom payloads contain only
      approved/disclosed metadata and never contain biological-data canaries;
      GA property settings, event-scoped dimensions, schema-filtered KPI
      report/query, synthetic reconciliation, and owner approval are recorded.
- [ ] If analytics is disabled, the hosted bundle has no measurement ID, Google
      CSP origins, consent notice, or effective loader/event path.
- [ ] Analytics unavailability never blocks an application feature.
- [ ] Local GUI, CLI, and Python remain telemetry-off.

### Exact artifacts and support matrix

- [ ] H-AC2 proves the protected workflow, index routing, OIDC, accepted-artifact reuse, and absence of a branch/PR upload path.
- [ ] H-AC3 proves that the wheel, sdist, and canonical deployable Web bundle are built once from
      the candidate and passed unchanged to all downstream jobs and the selected
      production deployment.
- [ ] Artifact filenames, target, version, size, SHA-256, contents, invocation, fixture, and results are recorded.
- [ ] H-AC1 proves metadata, entry points, licenses/notices, vendored assets, native-executable absence, Wasm, fonts, scientific fixtures, Gallery ownership, and external runtime fallback.
- [ ] No cache, credential, local path, debug file, unpublished dataset, recursive browser wheel, or unintended Gallery artifact is shipped.
- [ ] H-AC4 passes against the exact base wheel on the supported OS/Python matrix and the exact sdist on its declared qualification lanes.
- [ ] H-AC5 proves the separate `[export]` outputs and the installed-resource local GUI/offline workflow.
- [ ] Every claimed browser passes the declared full or smoke matrix; untested engines are not advertised as supported.

### Documentation and reproducibility

- [ ] A0 evidence and the A1 candidate commit are recorded.
- [ ] README, installation/support matrix, quickstarts, tutorials, release notes, compatibility reference, privacy text, and CLI/Python references match the candidate.
- [ ] Literal CLI/Python recipes execute from clean directories against the built candidate artifact.
- [ ] Documentation capture and strict Gallery validation pass against the candidate.
- [ ] Candidate manuscript/public figures reproduce from a clean archive using declared inputs and recipes and have readable-scale visual review evidence.
- [ ] Unavailable or non-redistributable figure inputs are resolved or the corresponding reproducibility claim is removed.
- [ ] Citation and archive placeholders are ready without claiming an identifier or publication route that is not yet live.

J-RC pass authorizes only a request for an explicitly approved RC tag,
TestPyPI upload, or staged deployment through Work package K.

## J-Final — Final software-release qualification

- [ ] Every post-RC change has a regression test, invalidated-evidence list, and complete affected/full-gate rerun.
- [ ] The authorised RC's TestPyPI wheel/sdist exact-download hashes, wheel H-AC4, and staged hosted version/hash/privacy gates passed; no unresolved staging failure remains.
- [ ] The final-version commit passes package/intended-tag/version/schema/citation/release-note/Web-wheel concordance.
- [ ] The exact final wheel, sdist, and Web bundle pass the complete artifact, clean-install, offline, and reproduction gates with new hashes.
- [ ] The final support and compatibility manifests are unchanged or explicitly requalified.
- [ ] No P0/P1 remains; every P2 disposition is present in release-facing records where relevant.
- [ ] The final release manifest binds commit, intended tag, artifact hashes, Web build stamp, documentation recipes, and figure evidence.
- [ ] Rollback, yank, patch-release, and post-publication smoke procedures are ready.

J-Final pass authorizes only a request for Work package K to perform the
separately approved final tag, PyPI publication, hosted deployment, Bioconda,
and archive actions.

## K-Publication — Archive and scholarly-publication checklist

- [ ] The final tag points to the J-Final commit and published artifact hashes match the accepted manifest.
- [ ] Clean PyPI install and hosted-app version/hash/privacy smokes pass after publication.
- [ ] The production hosted app was deployed from the accepted bundle, not an unrelated later `main` build.
- [ ] Bioconda metadata and smoke evidence refer to the same release baseline.
- [ ] Zenodo archive contains or references the accepted baseline and the version DOI is recorded.
- [ ] Paper figures reproduce from the final tag/archive.
- [ ] Preprint, manuscript, supplement, citation, support matrix, limitations, archive identifiers, and P2 disclosures agree with the exact release.
- [ ] Cover letter and AI-use disclosure are ready.
- [ ] The maintainer separately authorises the preprint revision and journal submission.

---

# 10. Definitions of done

## `v0.14.0` software release

The software release is complete when:

1. A user can upload ordinary Circular or Linear inputs without advanced LOSAT controls obstructing the primary workflow.
2. A user can anchor each complete circular record by source coordinate in Circular or Linear mode and preserve source-coordinate semantics.
3. A user can keep an overlapping feature on the main lane or move it one supported lane in either direction, even when automatic overlap resolution is off.
4. Existing proven-safe direct edits remain immediate, while admitted geometry edits use the no-analysis single-flight renderer contract.
5. Successful automatic candidates commit the complete renderer-owned bundle and preserve valid context; failed, cancelled, stale, or superseded candidates leave the committed Result intact.
6. Undo/redo produces the Result for restored user intent without adding a renderer-owned history item.
7. The same figure intent is reproducible through a saved session and supported Web, CLI, and Python surfaces.
8. The exact final wheel/sdist install and run on the declared native support matrix, and the exact packaged Web bundle passes the declared browser/offline/security matrix.
9. Hosted analytics either passes the complete opt-in, withdrawal,
   closed-payload, operator, and failure-independence contract, or is disabled
   with no Google-capable artifact path; every local surface is
   telemetry-incapable.
10. Documentation and public/manuscript candidate figures reproduce from the final release baseline.
11. Gate 0, J-RC, and J-Final have complete evidence with no unresolved P0/P1.
12. The separately authorised final tag, PyPI publication, hosted deployment, and required post-publication smokes match the accepted release manifest.

## Publication programme

The publication programme is complete when:

1. the exact final release is archived and its version DOI is recorded;
2. Bioconda and all installation/support claims point to the same baseline;
3. preprint, manuscript, supplement, figures, citations, limitations, and disclosure text match the archive;
4. the K-Publication checklist passes; and
5. the maintainer has authorised and completed the preprint revision and journal submission.

Close feature scope and move the release line to blocker-only work at Gate 0.
After the software release, keep `release/0.14` for patch and
reviewer-requested compatibility fixes and move feature development to
`v0.15.0`. Scholarly publication may continue without holding the software
release definition open.
