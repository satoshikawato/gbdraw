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

The final feature additions for `v0.14.0` are:

1. Per-record circular coordinate anchoring: place an arbitrary source coordinate at 12 o'clock in Circular mode or at the left edge in Linear mode.
2. Manual feature placement that remains active independently of automatic overlap resolution, with one explicit lane in each supported direction.
3. A substantially cleaner LOSAT/comparison user interface based on progressive disclosure.
4. Responsive preview Phase 1: retain the proven-safe direct preview edits and add debounced, single-flight renderer updates for record display start, manual placement, and existing label-layout reconciliation through a path that cannot run analysis.
5. Preview navigation improvements, especially zoom-to-selection and cursor-centred zoom.
6. A small annotation UX finish, principally annotation TSV download.
7. Formal PyPI publication and clean-install validation.
8. Privacy-preserving feature analytics with prior opt-in consent on the hosted Web app.
9. Documentation consolidation, release QA, archival, preprint revision, and journal submission.

After these items are complete, create a release candidate and allow only release-blocking fixes.

---

## 2. Repository instructions and working rules

Before making changes:

1. Read the repository-root `AGENTS.md` and `CLAUDE.md`.
2. Read `gbdraw/web/CLAUDE.md` before changing the Web application.
3. Inspect the current `main` and `docs_renovation` branches and preserve unrelated worktree changes.
4. Treat the committed canonical typed request and resources, current session format, and last atomically committed Result as the authoritative product contracts.
5. Keep proven-safe direct preview edits, but do not expand SVG-only state. Semantic intent remains history-backed and session/request-reproducible; geometry changes use the renderer.
6. Do not update reference outputs unless a reviewed rendering-contract change requires it.
7. Do not silently add compatibility aliases to fresh APIs. Persisted-session migration and fresh-input acceptance remain distinct concerns.
8. Do not replace this roadmap with a shorter plan before implementation. Update the status and completion ledger instead.

---

## 3. Current product baseline

The `0.14.0b0` development line already contains or is expected to contain the following major capabilities:

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

- Circular display rotation must not rewrite biological coordinates.
- Feature popups, annotation tables, comparison metadata, and exported source-coordinate values remain in source coordinates.
- Reverse complement, crop, rotation, and rendering transforms must not destroy source feature identity.

### 4.3 Cross-interface parity

A public figure-setting change that affects rendering must be representable through:

- Web state and saved session;
- canonical typed request;
- Python API or typed API;
- CLI option or generated helper table where appropriate.

The Web app may offer convenience interactions, but saved semantic intent must not be trapped in transient DOM state.

Preview effect flags, debounce/scheduler state, revision tokens, and viewport snapshots are Web runtime coordination data. They are not persisted figure intent. Persist semantic display-start and placement settings through the canonical request/session and supported programmatic surfaces.

### 4.4 No hidden analysis reruns

- A user edit may apply a proven-safe direct patch, schedule a renderer rebuild, mark analysis evidence stale, or combine those effects. These effects are not mutually exclusive classes.
- Automatic renderer work starts from an immutable snapshot of the committed canonical request/resources and overlays only allowlisted render-safe intent from the triggering draft revision. Clone mutable request/resource-index structures copy-on-write; do not deep-clone unchanged large resource payloads for each update.
- The render-only entry point must not invoke LOSAT, rediscover inputs, resolve comparisons, inspect analysis caches, or mutate analysis evidence.
- At most one renderer request is active and one replaceable latest draft is pending. Only the current draft revision may commit.
- Analysis-invalidating changes show stale state and require an explicit analysis action.

### 4.5 Local-first privacy

- Uploaded biological data remain local to the browser or local installation.
- Feature analytics may be enabled only after explicit opt-in on the hosted app.
- Rejecting analytics must result in no analytics request to Google.
- Local GUI, CLI, and Python use remain telemetry-free by default.

---

# 5. Work packages

## Work package A — Consolidate the `0.14` baseline and documentation branch

### Objective

Merge or otherwise reconcile `docs_renovation` with the final implementation baseline before release-candidate testing.

### Required work

- Resolve branch divergence deliberately; do not overwrite newer implementation changes with stale documentation assumptions.
- Regenerate executable documentation outputs against the final `0.14` API and defaults.
- Ensure all examples use stable, accession-pinned or repository-tracked inputs where appropriate.
- Ensure the following documentation agrees:
  - README;
  - installation guide;
  - quickstart;
  - tutorials;
  - Python API guide;
  - CLI reference;
  - session/request compatibility guide;
  - release notes;
  - Gallery tutorial metadata;
  - citation instructions.
- Remove statements that describe unreleased future features as current capabilities.

### Acceptance criteria

- All documentation commands execute successfully in a clean environment.
- All documented screenshots and generated figures are reproducible from the final code.
- The package version, session version, request schema, and documented compatibility table agree.
- The README accurately distinguishes hosted Web, Bioconda, PyPI, and source installs.

### Scope stop

Do not use documentation renovation as a reason to redesign unrelated interfaces after feature freeze.

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
Analysis settings changed. Search results are out of date.
[Run LOSAT and update]
```

Proven-safe direct edits may still update the committed preview. Only the explicit Generate/LOSAT action can clear analysis-stale state.

### Scheduler and worker contract

- Use trailing-edge 400 ms debounce for render-only work.
- Keep at most one render active and one replaceable latest draft pending.
- Do not terminate and restart the diagram worker for ordinary superseded keystrokes.
- Discard an obsolete response and render the latest pending draft next.
- Manual Generate removes queued preview work and has priority for the next renderer slot.
- Session import, Reset, mode change, result replacement, and disposal invalidate queued revisions through the same coordinator.

### UI and status

Provide one control for automatic geometry updates:

```text
☑ Auto-update layout
```

Direct colour and other proven-safe patches remain immediate when this control is off. After the first successful manual Generate, the default is on unless a restored session records the user's preference. Session hydration itself must not trigger an automatic render.

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
- Commit Result list, selected Result, feature catalogue, extracted/biological feature state, render metadata, editor/selection bindings, and committed canonical request/resources as one logical transaction.
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
- A failed update preserves the previous Result, catalogue, metadata, selection, viewport, export, and committed canonical request/resources.
- Automatic updates off perform no worker render; `Update layout` performs exactly one.
- Analysis-stale state performs no automatic geometry render and requires explicit Generate/LOSAT action.
- One continuous edit creates one history item; renderer completion creates none; undo/redo reconciles the restored state once.
- Circular single/grid/batch and Linear preserve Result topology and inactive batch outputs.
- Automatic commit preserves zoom/pan and stable selection.
- The accepted candidate uses the shared SVG sanitization/commit boundary once and introduces no reference-output change by itself.

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

### Objective

Finish the existing annotation workflow without expanding the annotation data model.

### Required addition

Add a visible `Download TSV` action to the Region Annotations panel using the existing annotation table encoder.

Recommended top-row layout:

```text
[Add set] [Import TSV] [Download TSV]
```

### Optional low-risk addition

Expose an `Annotate selection` shortcut near the feature-selection toolbar or feature popup. It should invoke the existing selected-feature annotation flow and ask for a destination set or create a new one.

### Acceptance criteria

- Exported TSV re-imports without semantic loss.
- Multi-record targets retain record selectors.
- Marks, labels, lanes, colours, fill, and hatch settings survive round-trip.
- Disabled/no-annotation state is handled clearly.

### Non-goals

Do not add:

- grouped selected-feature span as a new domain behaviour;
- drag-to-resize annotation endpoints;
- point callouts;
- full nested style editor;
- conversion into GenBank/GFF3 biological features.

---

## Work package H — PyPI release and packaging audit

### Objective

Make `pip install gbdraw` a first-class, tested installation route for Python integration, notebooks, CI, containers, and coding agents.

### Distribution decision for `v0.14.0`

Do not introduce a `gbdraw-core` distribution.

Use:

```text
Distribution name: gbdraw
Import name:       gbdraw
CLI command:       gbdraw
```

The core user experience should be:

```bash
python -m pip install gbdraw
```

```python
from gbdraw import read_genbank, draw_circular

record = read_genbank("input.gb")
diagram = draw_circular(record)
svg = diagram.to_svg()
```

### Base-install capabilities

The base installation should support:

- Python API;
- Circular and Linear CLI;
- GenBank/DDBJ parsing;
- GFF3 + FASTA parsing;
- static SVG;
- interactive SVG where supported;
- sessions and typed requests;
- annotation rendering;
- precomputed comparison rendering.

Keep non-SVG export optional:

```bash
python -m pip install "gbdraw[export]"
```

### Mandatory artifact audit

Build the actual wheel and source distribution, then inspect contents and size:

```bash
python -m build
unzip -l dist/*.whl
tar tf dist/*.tar.gz
```

Measure at least:

- Web vendor assets;
- Pyodide runtime;
- browser-side wheels;
- LOSAT WASM;
- native LOSAT binaries;
- fonts;
- generated Gallery material;
- accidental debug/test files.

### Packaging decision rules

- If the all-in-one wheel is practical and cross-platform, retain one distribution for `v0.14.0`.
- If a platform-specific LOSAT binary makes the wheel contract incorrect, remove that binary from the universal wheel and discover LOSAT externally.
- If GUI assets are a demonstrated installation obstacle, plan a later `gbdraw-gui` companion package; do not block `v0.14.0` unless PyPI publication is impossible.
- Do not rename the main distribution to `gbdraw-core`.

### Release automation

Add a GitHub Actions publishing workflow using PyPI Trusted Publishing/OIDC.

Expected release flow:

1. Checkout the release tag.
2. Run tests.
3. Build wheel and sdist.
4. Inspect or test artifact contents.
5. Install built wheel into clean environments.
6. Publish to TestPyPI for release-candidate validation.
7. Publish the final tag to PyPI.

### Clean-install matrix

Validate:

- Python 3.10;
- Python 3.11;
- Python 3.12;
- Linux;
- macOS;
- Windows.

Smoke tests:

- `import gbdraw`;
- `gbdraw --help`;
- Circular SVG;
- Linear SVG;
- session save/replay;
- local GUI launch if included;
- optional PDF/PNG export with `[export]`.

---

## Work package I — Feature analytics and privacy consent

### Objective

Measure feature adoption and failure without collecting biological content or sending analytics before consent.

### Hosted-app consent model

Use two meaningful choices, not `Accept all / Only necessary / Deny all`:

```text
Privacy settings

With your permission, gbdraw uses Google Analytics to understand
which features are used and to improve the app.

[Allow usage analytics]
[Continue without analytics]

Privacy details
```

### Required behaviour

Before consent:

- do not load the GA4 script;
- do not send analytics requests;
- do not set GA cookies.

If allowed:

- store the choice locally;
- load analytics only after consent;
- enable analytics storage;
- keep advertising-related consent denied;
- send only approved, enumerated events and parameters.

If declined:

- store the choice locally;
- keep the analytics wrapper as a no-op;
- do not contact Google Analytics;
- keep every application feature available.

Provide a persistent `Privacy settings` entry for later changes.

Local GUI, CLI, and Python remain telemetry-off by default.

### Exact user-facing copy

#### Main notice

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

### Event design

Use a small number of generic events with bounded parameters.

Core events:

```text
workflow_start
feature_progress
render_result
export_figure
session_saved
session_loaded
tutorial_begin
tutorial_complete
```

Emit `render_result` success only after the current candidate commits atomically. Emit failure only for a genuine failure of the current candidate. Debounce coalescing, pending replacement, stale/superseded completion, and cancellation are not render attempts or failures and emit no result event. Direct patches are not renders. A bounded `trigger` value may distinguish explicit analysis generation, manual no-analysis rendering, and automatic no-analysis rendering.

`feature_progress.stage` values:

```text
panel_open
configured
rendered
exported
```

Count a feature/stage at most once per session unless a separate repeat-use analysis explicitly requires otherwise.

### Stable feature identifiers

Initial enumeration:

```text
circular_rotation
manual_feature_placement
region_annotations
custom_tracks
depth_tracks
multi_record_layout
label_editing
legend_editing
interactive_svg
session_replay
losatn
tlosatx
losatp_pairwise
losatp_similarity_groups
losatp_collinear
```

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

Allowed properties must be bounded enums or coarse buckets, for example:

```text
mode = circular | linear
input_type = genbank | gff_fasta | session
record_count_bucket = 1 | 2 | 3_5 | 6_10 | 11_plus
sequence_size_bucket = lt_10kb | 10_100kb | 100kb_1mb | 1_10mb | gt_10mb
status = success | failure
error_code = predefined enum
```

### Product metrics

Primary metrics per feature:

```text
Successful adoption
= feature-active successful-render sessions
  / eligible sessions
```

```text
Export conversion
= feature-active export sessions
  / feature-active successful-render sessions
```

```text
Failure rate
= feature-active failed render attempts
  / feature-active render attempts
```

Exclude debounce replacements, stale/superseded completions, and cancellation from both the failure numerator and render-attempt denominator.

Treat analytics as a consenting-user sample, not a census of all users.

### Acceptance criteria

- No Google analytics request before consent.
- No analytics request after rejection.
- Consent can be changed later.
- Feature events are emitted after successful state transitions, not scattered button clicks.
- Automatic preview success is emitted only after the current candidate commits; scheduler churn emits no result event.
- Export events reflect the committed Result.
- Error codes are sanitised enums.
- Privacy documentation lists sent and unsent data.

---

## Work package J — QA, compatibility, and release engineering

### Objective

Turn the feature-complete beta into a release candidate and final publication artifact.

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
- packaging metadata;
- analytics event sanitisation.

#### Rendering tests

- Circular rotation across all track types;
- Circular and Linear manual placement, multipart/intron envelopes, labels, and canvas reservation;
- unchanged default output with no overrides and tolerance 0;
- annotation round-trip;
- deterministic output where expected;
- representative Circular and Linear figures.

#### Web tests

- LOSAT progressive disclosure;
- initial viewport;
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
- privacy consent paths;
- analytics no-op before/after rejection.

#### End-to-end acceptance tests

At minimum:

1. Circular figure with labels, GC/skew, annotations, rotation, resolver-off manual inward/outward placement, and session replay.
2. Linear multi-record figure with uploaded comparison, resolver-on Main/Auto precedence, and manual above/below placement.
3. Linear LOSATP pairwise or collinearity workflow.
4. Session save, close, reopen, render, and export.
5. PyPI wheel install and SVG generation.
6. Hosted app analytics allow and reject flows.

### Compatibility rules

- Recheck supported session versions and current writer version after the final schema decision.
- Do not bump the session version for branch-only intermediate formats.
- Do not create multiple unreleased schema versions during implementation.
- Keep old persisted-data migration separate from fresh API acceptance.
- Effect flags, debounce state, scheduler tokens, pending jobs, slow-pause state, and viewport snapshots are ephemeral runtime state. They must not create canonical request fields or cause a session-version bump. Persist only user-owned preview preference and semantic figure intent.

### Release-blocking severity

#### P0 — must fix

- data loss;
- incorrect biological coordinates;
- corrupted session;
- security/privacy failure;
- package cannot install/import;
- exported figure materially incorrect.

#### P1 — must fix before final release

- core workflow crash;
- automatic preview reaches LOSAT, analysis-cache mutation, or mutable input/resource reconstruction;
- concurrent renderer jobs run or an obsolete draft revision commits;
- renderer-owned Result, catalogue, metadata, selection, canonical artifacts, and export authority commit only partially;
- undo/redo produces a Result that does not match the restored user intent;
- automatic rendering unconditionally resets zoom, pan, or a still-valid selection;
- rotation inconsistency across tracks;
- manual placement or overlap tolerance cannot replay;
- last good preview lost after failure;
- consent rejection still contacts analytics;
- documented command does not work.

#### P2 — fix if practical; otherwise document

- minor layout defect;
- uncommon browser issue;
- non-blocking accessibility defect;
- cosmetic inconsistency.

After `rc1`, only P0/P1 and specifically approved P2 fixes should be merged.

---

## Work package K — Release, archive, preprint, and journal submission

### Objective

Use one exact software baseline for the package release, archived artifact, preprint revision, and journal submission.

### Release sequence

```text
complete work packages
    ↓
feature freeze
    ↓
v0.14.0rc1
    ↓
release-blocking fixes only
    ↓
v0.14.0 tag
    ↓
PyPI publication
    ↓
hosted Web deployment
    ↓
Bioconda update
    ↓
Zenodo archive / version DOI
    ↓
bioRxiv preprint v2
    ↓
journal submission
```

### Preprint v2 positioning

Replace a narrow “genome diagram generator” framing with:

> gbdraw is a local-first, browser-native, reproducible environment for constructing publication-quality genome diagrams, with CLI and Python interfaces for automation and integration.

### Claims to emphasise

- ordinary users can complete a figure in the GUI;
- CLI and Python are for automation, batch use, and embedding rather than mandatory completion steps;
- direct semantic editing remains reproducible;
- region annotations, per-record display start, manual overlap placement, feature-overlap tolerance, and comparison results are saved as project state;
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

### Post-submission branch policy

```text
release/0.14
  bug fixes
  reviewer-requested compatibility fixes
  patch releases only

main
  v0.15.0 development
```

The manuscript should cite the exact release and archive DOI. Update the manuscript to a patch release only if a relevant material defect is fixed.

### AI-assisted development disclosure

Prepare a transparent statement for the cover letter and Methods/Acknowledgements, for example:

> AI-assisted coding agents were used during software implementation and documentation. All architectural decisions, code review, testing, validation, and scientific interpretation were performed by the author.

---

# 6. Dependency order

Use the following order to reduce rework:

1. **A — Baseline/documentation reconciliation**
   - identify the final branch state and current schemas;
   - do not regenerate all assets yet if rendering contracts will still change.
2. **B — LOSAT/comparison UI cleanup**
   - largely independent of rendering geometry;
   - improves the working interface for subsequent testing.
3. **C — Circular coordinate transform**
   - establish source/display transform and record-aligned provenance before placement overrides depend on them;
   - complete C's model/planner work first, but hold its Phase 3 schema 6/session 41 writer until D's persisted fields are ready for the same atomic update.
4. **D — Manual feature placement**
   - bind to source identity through a record-aligned request-plan-to-renderer boundary, fold persisted fields into C's complete schema/session writer, and test under rotation, crop, reverse complement, duplicate record IDs, and custom feature slots;
   - treat C Phase 3 and D Phase 6 as one format-integration boundary rather than landing either partial writer;
   - keep manual placement independent from the Auto resolver and treat overlap tolerance as a separate render-only setting.
5. **E — Responsive preview Phase 1**
   - first complete or verify the Web performance plan's history/artifact separation and single SVG ingestion/commit boundary;
   - retain proven-safe direct patches, add composable effects, and route display start, placement, retained overlap tolerance, and existing label reconciliation through a debounced single-flight no-analysis path based on committed canonical request/resources;
   - commit the current candidate atomically and preserve the minimum viewport.
6. **F — Preview navigation**
   - add navigation controls and refined zoom behaviour on top of E's commit-time viewport-preservation baseline.
7. **G — Annotation TSV download and optional shortcut**
   - low-risk finish after state contracts are stable.
8. **H — Packaging/PyPI audit**
   - build actual artifacts from the feature-complete code.
9. **I — Analytics/privacy**
   - instrument stable feature identifiers and successful state transitions.
10. **A — Final documentation regeneration**
    - regenerate examples and screenshots against the final behaviour.
11. **J — Full QA and release candidate**
12. **K — Release, archive, preprint, and submission**

Analytics/privacy may be developed earlier in parallel, but final event semantics should be bound only after feature contracts are stable.

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
2. Defer `Annotate selection` shortcut; keep only annotation TSV download.
3. Defer non-essential preview-navigation polish; retain cursor-centred zoom and zoom-to-feature if possible.
4. Defer automatic performance fallback tuning; retain the preview effect policy, no-analysis boundary, scheduler, stale-response rejection, atomic commit, history semantics, and minimum viewport preservation.
5. Keep analytics events minimal, but do not drop consent/privacy correctness if analytics remains enabled.

Do **not** cut:

- source-coordinate correctness;
- manual placement while automatic overlap resolution is off;
- lane-1 placement in both supported directions for bidirectional combined-strand feature slots;
- session replay;
- last-successful-render transaction behaviour;
- proven-safe direct preview edits and the deliberately narrow display-start/placement auto-render scope;
- the dedicated no-analysis boundary, single-flight/latest-wins scheduling, and stale-response rejection;
- history-derived rendering and minimum viewport preservation;
- clean PyPI installation;
- privacy consent before analytics;
- release archive and reproducibility documentation.

---

# 9. Release-candidate checklist

## Product

- [ ] Sequence upload is visible in the initial desktop viewport.
- [ ] No-comparison mode performs no comparison work.
- [ ] Per-record display start works for all supported Circular and Linear coordinate consumers.
- [ ] Source coordinates remain correct in popups and outputs.
- [ ] Main and lane-1 directional placement work with automatic overlap resolution on and off and survive regeneration and session replay.
- [ ] Overlap tolerance preserves default behaviour at 0 and follows the documented bp boundary when enabled, unless it was explicitly cut before release freeze.
- [ ] Proven-safe direct preview edits remain immediate and use their established session/export owners.
- [ ] Per-record display start, placement, retained overlap tolerance, and admitted label reconciliation use one debounced render-only update without calling an analysis entry point.
- [ ] At most one automatic render is active; the latest pending draft wins and an obsolete response cannot commit.
- [ ] Failed, cancelled, stale, and superseded updates preserve the previous Result, catalogue, metadata, selection, viewport, export authority, and committed canonical request/resources.
- [ ] Auto-update off performs no automatic worker render, and manual `Update layout` performs one.
- [ ] Undo/redo renders the restored user intent once without adding a renderer-owned history item.
- [ ] Automatic commit preserves zoom/pan and any still-valid stable selection.
- [ ] Zoom-to-feature and core navigation work.
- [ ] Annotation TSV export re-imports losslessly.

## Session and API

- [ ] The current writer version is finalised once and contains the complete Work package C record-display and Work package D placement/tolerance contract; no partial schema 6/session 41 writer exists.
- [ ] Supported older sessions still load according to documented policy.
- [ ] Web, CLI, and Python describe equivalent display-start and placement intent.
- [ ] Deterministic feature identity survives rotation, reverse complement, crop-retained transforms, dormant crop exclusion/reactivation, duplicate record IDs, and duplicate-hash disambiguation.
- [ ] Committed canonical request/resources advance with the accepted Result; effect flags, scheduler state, revision tokens, and viewport snapshots are not persisted figure intent.

## Packaging

- [ ] Wheel and sdist contents are audited.
- [ ] No accidental debug files or generated artifacts are included.
- [ ] Clean PyPI-style install passes on supported Python versions.
- [ ] SVG output works from the base install.
- [ ] Optional export dependency works where supported.
- [ ] CLI and local GUI entry points behave correctly.

## Privacy and analytics

- [ ] GA4 is not loaded before consent.
- [ ] Rejecting analytics sends no analytics request.
- [ ] Allowing analytics sends only enumerated metadata.
- [ ] No biological content is present in analytics payloads.
- [ ] Privacy choice can be changed later.
- [ ] Local GUI/CLI/Python remain telemetry-off by default.

## Documentation

- [ ] README is current.
- [ ] Installation matrix is current.
- [ ] Quickstarts execute.
- [ ] Tutorials execute.
- [ ] Release notes are complete.
- [ ] Session compatibility page is correct.
- [ ] Paper figures reproduce from the release tag.
- [ ] Citation and archive instructions are ready.

## Publication

- [ ] Zenodo release archive is configured.
- [ ] Preprint v2 text matches the exact release.
- [ ] Manuscript and supplement use the archived version.
- [ ] Cover letter and AI-use disclosure are ready.

---

# 10. Definition of done

`v0.14.0` is done when all of the following are true:

1. A user can upload ordinary Circular or Linear inputs without advanced LOSAT controls obstructing the primary workflow.
2. A user can anchor each complete circular record by source coordinate in Circular or Linear mode and preserve source-coordinate semantics.
3. A user can keep an overlapping feature on the main lane or move it one supported lane in either direction, even when automatic overlap resolution is off.
4. Existing proven-safe direct edits remain immediate. Phase 1 automatically renders only per-record display start, manual placement, retained overlap tolerance, and admitted label reconciliation after debounce.
5. Automatic candidates come from the current draft revision plus committed canonical request/resources, run through a single-flight path that cannot invoke LOSAT, and commit only while current.
6. Successful automatic commits preserve the minimum viewport and commit renderer-owned state together. Failed, cancelled, stale, or superseded candidates leave the committed Result intact.
7. Undo/redo produces the Result for the restored user intent without adding a renderer-owned history item.
8. The same figure intent is reproducible through a saved session and the supported programmatic surfaces.
9. `pip install gbdraw` works in a clean supported Python environment.
10. Hosted feature analytics are strictly opt-in and never transmit biological content.
11. Documentation and manuscript figures reproduce from the tagged release.
12. The tagged code is archived, the preprint is revised, and the journal manuscript is submitted without waiting for plasmid, ORF, or Genome Overview development.

Once these conditions are met, stop feature development on the release branch and move future work to `v0.15.0`.
