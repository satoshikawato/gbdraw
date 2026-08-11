# Implementation plan: content-aware diagram layout overhaul

| Field | Value |
| --- | --- |
| Date | 2026-08-04 |
| Repository baseline | HEAD `afff6131`; plan inspected the current dirty working tree on 2026-08-04 |
| Status | In progress; WP0–WP4 and WP6 complete; WP5/WP7/WP8 remain open pending explicit-approval Chromium checks and captures |
| Primary defect | Circular records, external labels, legends, and titles are composed from nominal canvas dimensions instead of measured visible-content bounds |
| Scope | Circular single/grid/batch, Linear, Python/CLI/Web/session replay, interactive editing, SVG/export geometry, reference figures, documentation, and Gallery assets |

## 1. Objective

Replace the current collection of canvas-growth, legend-offset, title-placement,
and browser-side repositioning rules with one content-aware composition contract.
The finished system must place legends next to the content users perceive as the
diagram, center titles on that diagram rather than on the diagram-plus-legend
canvas, preserve manual Web edits, and produce unclipped output in every supported
format.

This is an intentional default-geometry change. It is not a request to shrink
feature tracks, change record scale, remove labels, or add a user-facing spacing
control.

## 2. Decision summary

The implementation will:

1. separate plot-space geometry from final figure composition;
2. measure the authoritative visible bounds before placing docked decorations;
3. add one pure, mode-neutral composition planner under `gbdraw/layout/`;
4. make Circular single-record, Circular grid, and Linear assembly consume that
   planner;
5. treat `left`, `right`, `top`, and `bottom` as docked placements, and
   `upper_left`, `upper_right`, `lower_left`, and `lower_right` as overlays;
6. keep Python as the owner of generated layout policy;
7. replace Web layout heuristics with an interpreter of versioned,
   Python-emitted composition metadata plus explicit user drag deltas;
8. keep current public API signatures, CLI flags, typed request schema 5, and
   session schema 40;
9. intentionally refresh geometry-dependent references and public figures only
   after the new relational acceptance tests and visual review pass.

No compatibility branch may continue to calculate new-output layout using the
old fixed-width rules. Old saved SVG results may be normalized once at the Web
load boundary, then must use the current composition contract.

## 3. Current implementation checkpoint

### 3.1 Measured symptom

The motivating figure is
`docs/images/t-py-06/python_precomputed_circular_rings.svg`.
A local Chromium measurement of the committed SVG gives:

| Metric | Current value |
| --- | ---: |
| Canvas | 2133.1 × 1031.4 px |
| Primary visible content, excluding title and legend | x=190.4…1206.8; width=1016.4 px |
| Legend | x=1433.1…2095.1; width=662.0 px |
| Visible content-to-legend gap | 226.3 px |
| Gap as a share of primary content width | 22.3% |
| Circular axis center | x=700.0 px |
| Total-canvas center used by the title | x=1066.5 px |

The configured side-legend gap itself is only about 33 px. Most of the visible
226 px gap is unused space between the painted content and the fixed 1400 px
Circular base-canvas edge.

A browser audit of the current non-interactive documentation SVGs found:

| Family | Count | Median visible side gap | Range | More than 150 px |
| --- | ---: | ---: | ---: | ---: |
| Circular, side legend | 30 | 195.4 px | 40.2…317.8 px | 21/30 |
| Linear, side legend | 12 | 100.0 px | 78.3…108.1 px | 0/12 |

Circular is the urgent defect. Linear is included because it has a separate
fixed 100 px policy and would otherwise remain a competing layout owner.

### 3.2 Root causes

#### Circular single-record

`gbdraw/canvas/circular.py` currently mixes plot geometry and final composition:

- `CircularCanvasConfigurator.calculate_dimensions` selects the configured
  nominal width (`1400` with labels, `1000` without labels);
- `resolve_circular_side_legend_geometry` derives gaps from legend width and
  canvas height;
- `recalculate_canvas_dimensions` places a side legend after the nominal width,
  not after visible content.

`gbdraw/diagrams/circular/assemble.py` then performs additional, later layout:

- external-label and radial-content canvas growth;
- label-versus-legend collision resolution;
- top/bottom legend placement;
- whole-canvas translation and viewBox synchronization.

`gbdraw/api/diagram.py` adds or moves plot titles after lower-level assembly,
which can expand and translate the canvas again.

#### Circular multi-record

`assemble_circular_diagram_from_records` renders complete subcanvases and then
serializes and reparses their SVG to estimate visible bounds. The helper chain
`_estimate_element_local_x_bounds`,
`_estimate_element_local_y_bounds`,
`_estimate_subcanvas_content_bounds`,
`_estimate_subcanvas_horizontal_insets`, and
`_estimate_subcanvas_vertical_insets`:

- understands only a subset of transforms;
- does not use authoritative text metrics;
- treats numbers found in path data as approximate extents;
- recovers information that was available while the diagram was structured.

The grid, shared legend, and shared title then use another set of spacing
constants and placement branches.

#### Linear

`LinearCanvasConfigurator.recalculate_canvas_dimensions` in
`gbdraw/canvas/linear.py` places a right legend at:

`horizontal_offset + alignment_width + 2 * canvas_padding`.

With the default 50 px padding, the visible gap is designed to be about 100 px.
`assemble_linear_diagram` separately handles record rows, title reserves,
top/bottom legend stacks, final height, and alternate legend-orientation
viewBoxes.

#### Web editor

The generated SVG is laid out in Python, but post-generation Web changes are
laid out again in JavaScript:

- `gbdraw/web/js/app/legend-layout/canvas-actions.js` reconstructs base geometry
  using DOM `getBBox()`;
- `legend-layout/reposition-actions.js` recomputes canvas, title, diagram, and
  legend placement;
- `legend/layout-actions.js` owns browser-side legend reflow and repeats Python
  spacing ratios;
- `drag-actions.js` and `diagram-drag.js` maintain manual transforms;
- imported sessions can silently enter a metadata-free inference path.

The Python and browser text engines need not agree, so fresh generation,
orientation switching, legend editing, session restore, and export can drift.

### 3.3 Duplicated policies to retire

The migration must remove or reduce these to thin applications of one plan:

- `assemble._sync_canvas_viewbox` versus
  `api.diagram._sync_drawing_canvas_size`;
- two whole-canvas translation helpers;
- multiple definition/title local-bound estimators;
- single-, multi-, Linear-, and browser-specific side gaps;
- single- and multi-record top/bottom title/legend stacks;
- Python SVG reparse bounds estimation;
- browser reconstruction of current-output base geometry;
- unversioned `data-horizontal-viewbox` and `data-vertical-viewbox` as the
  primary layout contract.

### 3.4 Dirty-worktree baseline

At planning time the repository has broad pre-existing modifications across
production layout code, Web code, tests, reference SVGs, documentation, and
Gallery assets. This plan was derived from the current working-tree contents;
`afff6131` identifies HEAD, not a claim that those contents match HEAD.

Do not reset or overwrite those changes merely to create a clean baseline. WP0
must classify in-scope diffs, record the actual starting bytes and test results,
and separate them from changes made while executing this plan. Gates that say
“unchanged” mean unchanged from that recorded WP0 baseline, not necessarily
unchanged from `afff6131`.

## 4. Scope and surface contract

### 4.1 Mode and topology matrix

| Mode | Topology | Required result |
| --- | --- | --- |
| Circular | single | One measured primary-content box; one optional title; one optional legend |
| Circular | one-record grid | The same record assembly and bounds contract as single, wrapped by the grid composer |
| Circular | multi-record grid | Direct per-record content bounds, visible-content row/column gaps, one shared title, one shared legend |
| Circular | batch | Each output follows the single-record contract independently |
| Linear | one record | Measured record, labels, tracks, definition, and scale as one primary-content box |
| Linear | several rows/columns | Existing row/scale solver feeds one final primary-content box into the shared composer |
| Web | fresh generation | Must display Python's generated composition unchanged |
| Web | legend/title edit | Must recompute from Python-emitted constraints and current measured decoration size |
| Web | manual drag | Must apply a user delta after automatic composition and preserve it through history/session/export |
| Session replay | current and supported legacy | Must converge on the same current render or normalized editor contract |

### 4.2 Surface decision matrix

| Surface | Current public contract | Planned change |
| --- | --- | --- |
| Package-root Python API | Circular/Linear option objects and `Diagram` methods | No new field; no signature/default change |
| `gbdraw.api` | typed options, requests, render/session helpers, low-level assemblers | Preserve public signatures and return shapes; add internal assembly result only |
| CLI | `--legend` values and defaults | No flag/value/default change |
| Web controls | existing legend and plot-title position controls | No new control; defaults remain mode/profile-specific |
| Canonical request | schema 5 | No schema bump |
| Session | schema 40 and supported readers | No schema bump; generated SVG content gains internal metadata |
| Static SVG | width, height, viewBox, group transforms | Intentional geometry change; semantic IDs/hooks remain stable |
| Interactive SVG | viewport plus standalone interaction metadata | Preserve behavior; composition metadata must survive sanitization/export |
| PNG/PDF/EPS/PS | converter output dimensions | Follow the new SVG viewBox without clipping or scale change |

The Web Linear default (`bottom`) may remain different from the Python/CLI
default (`right`). The shared planner owns what each requested value means; it
does not force all surfaces to choose the same default.

## 5. Non-goals

- Do not alter genome-coordinate scaling, Circular radius, feature width,
  record order, track order, comparison endpoints, or biological filtering.
- Do not redesign the legend's entries, typography, colors, gradients, or
  ordering.
- Do not change the external-label placement algorithm except where overlay
  legends require collision resolution.
- Do not add `legend_gap`, `canvas_padding`, or automatic/manual layout modes to
  Python, CLI, Web, or session schemas in this change.
- Do not crop a finished SVG by reparsing arbitrary path data.
- Do not add renderer-specific composition planners.
- Do not change top-level semantic group IDs or wrap the diagram in a new group
  if that breaks the Web editor's current topology contract.
- Do not hand-edit generated SVGs, sessions, thumbnails, or screenshots.
- Never regenerate, replace, or edit `examples/gbdraw_social_preview.png`.

## 6. Target layout contract

### 6.1 Terminology

| Term | Definition |
| --- | --- |
| Plot space | Local coordinates used while records, tracks, labels, and definitions are built |
| Primary content | Records, tracks, axes, scale, comparisons, definitions, external labels, leaders, and their stroke extents; excludes plot title and legend |
| Decoration | Plot title or legend |
| Dock | A decoration placed outside one edge of primary content |
| Overlay | A decoration placed inside a requested primary-content corner |
| Composition | Placement of primary content and decorations into the final non-negative canvas/viewBox |
| Automatic placement | Planner-produced base transform |
| User delta | Manual Web drag offset applied after automatic placement |

### 6.2 Shared spacing tokens

The first implementation will use one internal, immutable spacing policy:

| Token | Value | Meaning |
| --- | ---: | --- |
| `edge_padding_px` | 16 | Minimum final viewBox clearance around painted bounds |
| `dock_gap_px` | 24 | Visible primary-content-to-docked-legend gap |
| `title_gap_px` | 20 | Primary-content-to-title gap |
| `stack_gap_px` | 20 | Gap when title and legend occupy the same top/bottom side |
| `overlay_clearance_px` | 8 | Minimum overlay-legend clearance from labels/obstacles |

Stroke safety is derived from each item's actual maximum stroke width and is
added to its bounds. It is not another global padding constant.

These are internal rendering defaults, not new user options. Put them in one
`CompositionSpacing` value owned by the composition module. Do not copy the
numbers into Circular, Linear, or JavaScript code. A later product change can
expose spacing only if there is demonstrated user need and the full
cross-surface contract is designed separately.

### 6.3 Geometry invariants

For every generated figure:

1. every painted item, including text, strokes, leaders, and markers, lies
   inside the viewBox with at least its required safety clearance;
2. `width`, `height`, and `viewBox` describe the same final rectangle;
3. docked legend placement is relative to authoritative primary-content bounds:
   - right: `legend.left = primary.right + dock_gap`;
   - left: `legend.right = primary.left - dock_gap`;
   - bottom: `legend.top = primary.bottom + dock_gap`;
   - top: `legend.bottom = primary.top - dock_gap`;
4. a side legend is centered on the primary-content box along the orthogonal
   axis; title lanes do not pull it away from the plot;
5. top/bottom legends are centered on primary content;
6. top/bottom plot titles are centered on primary content, not total canvas;
7. a center-positioned Linear title remains an intentional overlay and does not
   reserve an outer lane;
8. when title and legend request the same top/bottom side, the legend stays
   nearest the primary content and the title occupies the outer lane;
9. if a side legend and a long top/bottom title would intersect, measure or wrap
   the title against primary width first, then move the title outward on its
   declared side; never move the plot away from the legend;
10. changing among docked legend positions does not change record, track, label,
    or comparison geometry; it changes only outer composition transforms;
11. corner placements remain overlays and are the only legend positions allowed
    to invoke label-versus-legend collision resolution;
12. final translation may change absolute SVG coordinates, but relative plot
    geometry and scale remain unchanged;
13. left/right placement is mirror-symmetric for symmetric inputs;
14. empty or `none` legends reserve no space.

### 6.4 Overlay behavior

Map public values to an internal closed type:

- dock: `left`, `right`, `top`, `bottom`;
- overlay: `upper_left`, `upper_right`, `lower_left`, `lower_right`;
- absent: `none`.

For an overlay, keep the requested corner and use this deterministic resolution
order:

1. place at the corner with edge and overlay clearances;
2. for Circular external-label conflicts, return the conflict set to the
   existing mode-owned label solver, which shifts the smallest valid label
   block while preserving order and leader validity;
3. shift the legend within the requested quadrant;
4. expand the corresponding outer edge only when no in-place solution exists.

Port the current label-shift behavior only to this overlay path. Docked legends
must be placed after primary bounds are final, so they cannot collide with
external labels and must not perturb them.

The shared composer may detect obstacle intersections and propose decoration
placements, but it must not learn Circular label ordering or mutate labels.
After a mode-specific label shift, rebuild the affected bounds and submit the
updated request to the pure composer.

### 6.5 Composition sequence

Every mode follows the same sequence:

1. resolve and render plot-space geometry;
2. collect exact primary-content bounds from structured geometry;
3. measure title and legend in their required orientation;
4. resolve dock/overlay placement and title lanes;
5. take the union of all placed bounds;
6. add edge padding;
7. translate all targets into a non-negative final viewBox;
8. apply the plan once to SVG groups and root size;
9. emit internal composition metadata;
10. serialize/export without a later independent layout pass.

## 7. Target architecture and ownership

### 7.1 Reuse the existing AABB type

Extend `gbdraw/layout/spatial.py:Aabb` rather than introducing a parallel bounds
tuple abstraction. Add pure operations as needed:

- `width` and `height`;
- `translated(dx, dy)`;
- `expanded(x, y=None)`;
- `intersects(other, clearance=0)`;
- `union_aabbs(...)`;
- finite-value validation and explicit empty handling.

Keep `AabbIndex` behavior compatible.

### 7.2 Add one mode-neutral composer

Add `gbdraw/layout/composition.py` with frozen internal contracts similar to:

- `CompositionSpacing`;
- `LegendPlacement` with dock/overlay/none semantics;
- `CompositionItem` containing role, authoritative local bounds, and placement
  constraints;
- `CompositionRequest` containing primary, optional legend/title, declared
  positions, and overlay obstacles;
- `CompositionPlacement` containing a role transform and final bounds;
- `CompositionPlan` containing canvas bounds, viewBox, primary bounds,
  placements, and resolved spacing.

The module must be pure: no `svgwrite` objects, DOM parsing, configuration-file
reads, or mode-specific feature logic.

### 7.3 Keep mode-specific measurement with the mode

The shared composer does not guess feature or path bounds.

Circular measurement uses:

- `CircularRadialLayout.outer_content_radius_px` and slot bands;
- exact external-label oriented/AABB geometry already calculated during label
  placement;
- leader endpoints;
- definition metrics;
- annotation and stroke extents.

Linear measurement uses:

- `LinearLayoutPlan` record placements and widths;
- `VerticalBand`/`CollisionBand` results;
- measured definitions and labels;
- track-slot paint bands;
- comparison corridors;
- ruler/scale extents.

Legend measurement remains owned by
`LegendDrawingConfigurator.measure_legend` and the mode-specific legend layout
modules. Change its input from a mutable final canvas to an explicit measurement
constraint:

- side/overlay legends use their intrinsic vertical layout;
- top/bottom legends wrap against primary-content width;
- an intrinsically wider legend may widen the final union, but nominal canvas
  width may not create blank space before it.

Add authoritative local bounds to `LegendMeasurement` so callers stop deriving
`-0.5 * color_rect_size` independently.

### 7.4 Separate plot configurators from final composition

After migration:

- `CircularCanvasConfigurator` owns radius, base plot coordinates, and local
  drawing creation, not docked legend reservation;
- `LinearCanvasConfigurator` owns record-axis/track coordinate inputs, not final
  legend/title canvas expansion;
- `CompositionPlan` is the only owner of final width, height, viewBox, legend
  transform, title transform, and outer translation.

Keep `canvas.circular.width.with_labels` and `without_labels` only as local
plot-space workspace inputs while existing drawers still need those coordinate
origins. They must no longer reserve final output width. Remove or rename those
configuration fields only in a separately approved public-config cleanup; do
not leave them ambiguously described as final canvas dimensions.

If `legend_offset_x`, `legend_offset_y`, or mutable total-size fields must remain
during an intermediate phase, they may only receive values copied from a
`CompositionPlan`. Delete their independent calculations before completion.

### 7.5 Preserve public low-level assembler contracts

`assemble_circular_diagram` and other release-backed low-level entry points keep
their current public return shapes. Add an internal result such as
`CircularAssemblyResult(drawing, legend_measurement, content_bounds, ...)` and
use a private/internal entry point where multi-record assembly needs the extra
geometry. The public wrapper returns the same values it returns today.

Do not attach composition planning to `gbdraw.api` merely to share it. Both CLI
and public APIs already converge on the same internal builder; keep the new
owner below those adapters.

### 7.6 Apply plans without breaking SVG topology

Add one renderer-side application helper that:

- composes translations with existing transforms rather than replacing them;
- keeps current top-level group IDs and semantic attributes;
- updates root size/viewBox once;
- does not insert a wrapper that breaks Web group discovery;
- preserves copied definitions and every local `href`/`url(#...)` reference.

Remove the duplicate whole-canvas translation and size-sync helpers after all
callers use this path.

### 7.7 Versioned internal SVG composition metadata

Generated SVGs used by the Web editor must expose enough data to recompute
decoration placement after legend content or orientation changes without
reconstructing the original plot:

- composition schema version;
- authoritative primary bounds in final SVG coordinates;
- spacing tokens;
- title role, side, and measured local bounds;
- automatic legend/title transforms;
- stable role/group selectors needed to apply a new plan.

Schema 1 also carries two Python-owned policy blocks:

- `legendReflow`, with renderer-derived legend metrics when a legend is
  present and an explicit null value when it is absent;
- `overlayPolicy`, with `quadrantBoundaryRatio`, `candidateScoreOrder`,
  `canvasGrowthCandidateOrder`, and `canvasGrowthScoreOrder`.

The current interpreter validates the complete field set and rejects missing,
null, boolean, string, empty, or non-finite numeric values where numbers are
required. Only the schema-free legacy adapter may synthesize the frozen v0
policy. Python-to-JavaScript oracle tests compare dock, overlay, title-stack,
no-legend, and canvas-growth plans in both runtimes.

Use internal `data-gbdraw-*` attributes or a sanitizer-safe equivalent. These
attributes are not added to the public semantic selector table in
`docs/SVG_SEMANTIC_HOOKS.md`; that document already distinguishes internal
layout bookkeeping from supported semantic hooks.

The Web interpreter may measure an edited legend with `getBBox()`, but it may
not choose new gaps, corner percentages, or canvas-growth policy. Those values
come from metadata emitted by Python.

### 7.8 User offsets are separate from automatic layout

Model the effective browser transform as:

`effective_transform = automatic_composition_transform + user_drag_delta`.

Legend reflow, orientation switching, title movement, history, Reset, session
save/load, and export must preserve that distinction. Reset writes a zero delta
and serializes even when no legend exists. Transform handling must preserve
scale/rotate/matrix/chained transforms instead of replacing them with a plain
`translate(...)`.

## 8. Compatibility decisions

### 8.1 Public API and request compatibility

- Keep `CircularOutputOptions.legend` and `LinearOutputOptions.legend` names,
  defaults, and accepted values.
- Keep low-level assembler signatures and return shapes.
- Keep `gbdraw.api.__all__` and
  `tests/fixtures/public_contract.json` unchanged.
- Keep canonical request schema 5 and request encoding unchanged.
- Keep CLI help and `docs/CLI_Reference.md` option sets unchanged.
- Treat changed SVG dimensions/transforms as the intended behavior change.

If implementation discovers that a public signature or persisted field is
actually required, stop that work package and amend this plan before adding it.

### 8.2 Session compatibility

No session version bump is planned because user-owned state does not change.
`ui.layoutPreferences` continues to own requested semantic positions, and the
saved SVG continues to own manual post-generation edits.

When both sources are present, explicit `ui.layoutPreferences` wins over the
last canonical request because it represents the current editor draft. If it
is absent, import falls back to that request. Legacy fields are migrated once,
and partial preferences are filled from the current mode defaults.

At import:

1. a current SVG with composition schema 1 is validated and used directly;
2. a supported old SVG without the schema enters one explicit legacy
   normalization adapter;
3. the adapter measures known semantic groups once, reads any saved generated
   baseline and `ui.layoutPreferences`, creates an in-memory schema 1 baseline,
   and records that the source was legacy;
4. subsequent edits use only the current interpreter;
5. malformed schema 1 metadata produces an actionable error or requires
   regeneration; it must not silently fall back to the retired heuristic.

Import is transactional. The draft is committed only after `afterLoad` and
asynchronous record discovery/reparse complete. Any error restores owner,
transient, and app-local state, including feature-list scroll position, depth
track controls, and the selected pairwise orthogroup. Suppression flags and
generation counters prevent stale asynchronous completions from mutating the
restored state, and rollback does not reparse restored files.

When a legacy session contains enough generated-baseline evidence, derive the
manual legend/diagram/title deltas by comparing its actual transforms with the
new automatic anchors. When that evidence is absent, exact pre-drag geometry is
not recoverable: preserve the imported appearance as the normalized automatic
baseline, initialize deltas to zero, and make Reset return to that imported
appearance. Test and document this compatibility behavior instead of guessing a
historical automatic position.

Re-saving the session may store the normalized current SVG without changing the
session envelope version.

### 8.3 SVG and interactive compatibility

Preserve:

- deterministic and unique document-local IDs;
- all documented `data-gbdraw-*` semantic hooks;
- record index/ID semantics, including duplicate record IDs;
- `_gbdraw_track_slot_geometry` schema and run-metadata aggregation;
- interactive feature/match/annotation payloads and popup behavior;
- root original-viewBox handling in standalone interactive SVG;
- sanitization of scripts, event handlers, unsafe URLs, and foreign content.

## 9. Implementation work packages

All work packages start as not started. A package may be marked complete only
with the evidence named in its gate.

### WP0 — Freeze baseline and relational requirements

Status: [x] Complete

Execution evidence (2026-08-04):

- Baseline: `HEAD afff613156cf2eb0c7b00bcf639e41c391ca97f6` on
  `docs_renovation`. The broad tracked working-tree status consists of byte/EOL
  changes: the complete in-scope diff is empty with
  `git diff --ignore-space-at-eol -- <in-scope paths>`. Those starting bytes are
  retained as the execution baseline; the untracked plan is retained input and
  `tools/audit_svg_composition.py` is the first plan-created file.
- Raw-byte aggregate SHA-256 baselines: reference SVGs (16)
  `20ee27c770ba63c2a5359210b58541a869801d03030fbb3ce037ed2c296e2592`;
  docs SVGs (53)
  `c88d4b5072848c0eb4bb7a529a6558123e2bceb43bb13394860a659088a03cb9`;
  public contract
  `3dccabde18edddc1cb5caddf69cf7c264792043eacf7e18ffb05d6beaabed09b`;
  Gallery sessions/sources/examples/thumbnails
  `9a2e62ee14bf0b522e33536603ac2224b2a1476071a9aa66603f6b7777990d3b` /
  `aed7371bd7366c7c8e983ea47649ad242cf9eb32846fa66f0fea90e0ec714492` /
  `0f0de9feacfa0083f81631914b2270a12d2d5c8edc0daa51d42a2b1bc513bfdf` /
  `f3a762ced0fcb991b4fb8f678b16163d4b7b2367900da8c485757ae31d1e9c42`;
  Gallery catalog/tutorial media
  `ad5f34b23b8e485f19aff597ba03a7d80f8c119bfa3874f98f4b2c861ad06481` /
  `296c2e2d5193dd64dca0a4d6802e327d4d9ca66d180a209b6b79a2c4e76e5cf6`;
  protected social preview
  `b5e5341fbd8b564626221e79e2c329ba559c1ea58fbb0e3ab65a2f8144d67ee2`.
- Browser audit:
  `python tools/audit_svg_composition.py docs/images --json-out
  /tmp/gbdraw-docs-composition-baseline.json` measured 53 SVGs with no clipped
  primary/title/legend bounds. T-PY-06 measured canvas 2133.1 x 1031.4 px,
  primary x=190.4...1206.8, right gap 226.34 px, and axis center x=700.0.
  The side-legend census reproduced section 3: Circular 30 files, median
  195.36 px, range 40.23...317.76 px, 21 over 150 px; Linear 12 files,
  median 100.0 px, range 78.35...108.11 px, none over 150 px.
- Baseline tests: Circular focused group 256 passed; Linear focused group 188
  passed; API/request/session contracts 441 passed; output comparison 16
  passed. The six section 11.3 Node commands also passed. No pre-existing
  failure was observed in those selected groups. A later direct run of
  `tests/test_public_contract.py::test_public_api_and_cli_contract_snapshot`
  exposed a pre-existing Linear CLI hash mismatch; the same mismatch reproduces
  from a clean `git archive` of the recorded HEAD, while
  `tests/fixtures/public_contract.json` remains byte-identical to the WP0 hash.

Acceptance scenario sources:

| Scenario | Declared source and expected semantic content |
| --- | --- |
| Circular, labels absent | `circular_basic` reference case; one record, feature/GC/skew/tick groups, no external-label group content |
| Circular, horizontal labels | `circular_with_labels` reference and `tests/test_circular_label_placement.py`; label text and leaders retain record identity/order |
| Circular, radial labels | `circular_radial_labels` reference and `tests/test_circular_radial_labels.py`; oriented text AABBs and leaders are present |
| Dense mitochondrial labels | `tests/test_inputs/HmmtDNA.gbk` through T-PY-01 and focused label tests; dense feature labels and legend entries remain present |
| Wide gradient legend | T-PY-06; three comparison rings, 106 hits, feature legend, identity gradient, and plot title |
| Circular one-record/multi-record grid | `tests/test_circular_multi_canvas.py` plus H-PY-01; one shared legend/title and indexed record wrappers |
| Linear one record/multi-row | `linear_basic` / `linear_multi_genome` references and Linear multi-record tests; axes, definitions, record order, and comparisons remain present |
| Titles on supported sides | Circular multi-canvas and Linear track-layout title tests; one title with the requested semantic position |
| No/empty legend | `tests/test_legend_measurement.py` and no-legend Circular/Linear cases; zero measurement and no reserved legend space |

Edit/add:

- add `tools/audit_svg_composition.py`;
- add focused composition scenario fixtures only if needed;
- add baseline assertions to existing layout tests without changing production
  output.

Changes:

1. Capture `git status`, in-scope diffs (including an
   `--ignore-space-at-eol` review), hashes of layout-sensitive outputs, and the
   current focused-test results. Classify each pre-existing in-scope change as
   retained input, superseded implementation, generator-owned output, or
   unrelated/out of scope. Do not clean the tree.
2. Implement a read-only local-browser audit that measures root viewBox,
   primary bounds, legend bounds, title bounds, clipping, edge margins, and
   dock gaps from semantic groups.
3. Report JSON and a compact table for a file or directory.
   Freeze this command contract:

   ```bash
   python tools/audit_svg_composition.py docs/images \
     --json-out /tmp/gbdraw-docs-composition.json
   ```

4. Record the current T-PY-06 values and the full docs Circular/Linear census
   from section 3.
5. Define selectors by documented roles/IDs; do not treat title or legend as
   primary content.
6. Add test fixtures covering:
   - Circular with no external labels;
   - horizontal and radial labels;
   - dense mitochondrial labels;
   - wide gradient legends;
   - Circular one-record and multi-record grids;
   - Linear one record and multi-row records;
   - titles on each supported side;
   - no legend and empty legend.

Gate:

- the audit reproduces section 3 within 1 px;
- WP0 itself does not change production SVG bytes from the recorded starting
  hashes;
- pre-existing focused-test failures and in-scope diffs are recorded separately
  from plan-caused changes;
- each acceptance scenario has a declared source and expected semantic content.

### WP1 — Add pure composition geometry

Status: [x] Complete

Execution evidence (2026-08-04):

- Added immutable `Aabb` operations and the pure composition request/plan
  module with the single shared 16/24/20/20/8 px spacing policy.
- Added authoritative stroke-safe `LegendMeasurement.local_bounds` and changed
  every production/test call site to pass explicit placement and wrap width.
- Pure/measurement verification passed: 82 tests across
  `test_composition_layout.py`, `test_spatial_index.py`, and
  `test_legend_measurement.py`; 266 affected collinearity, SVG-ID, and Circular
  multi-canvas tests also passed.
- `TestOutputComparison` passed 16/16, Ruff and `git diff --check` passed, and
  the reference SVG, docs SVG, and public-contract raw-byte hashes remain the
  recorded WP0 values. The separately recorded pre-existing Linear CLI fixture
  mismatch is unchanged and is not caused by WP1.

Edit/add:

- `gbdraw/layout/spatial.py`;
- new `gbdraw/layout/composition.py`;
- `gbdraw/layout/__init__.py` only if internal exports are useful;
- new `tests/test_composition_layout.py`;
- `tests/test_legend_measurement.py`.

Changes:

1. Extend `Aabb` with the shared operations in section 7.1.
2. Implement the immutable request/plan contracts.
3. Implement dock, overlay, title-lane, union, padding, and final-translation
   logic.
4. Add explicit validation for non-finite sizes, negative spacing, missing
   primary bounds, and unknown placement values.
5. Add local legend bounds to `LegendMeasurement`.
6. Test the planner without `svgwrite` or mode fixtures.

Required tests:

- exact 24 px dock gap on all four sides;
- exact 16 px outer padding for the limiting item;
- title centered on primary content;
- same-side title/legend stacking order;
- tall/wide decoration union;
- left/right and top/bottom symmetry;
- negative local coordinates;
- empty decoration behavior;
- overlay obstacle detection, legend movement, and deterministic canvas-growth
  fallback (label movement remains a Circular test);
- input immutability and repeatability.

Gate:

- pure tests pass;
- Circular/Linear output-comparison status matches the recorded WP0 baseline;
- no existing generated output changes relative to the WP0 hashes;
- no assembler has a second copy of the new spacing constants.

### WP2 — Migrate Circular single-record and batch composition

Status: [x] Complete

Execution evidence (2026-08-04):

- Circular assembly now carries authoritative primary, overlay-obstacle,
  legend, title, and track-slot geometry into one shared composition plan;
  docked legend placement no longer changes record-local label geometry.
- The plan-listed Circular regression group plus the new composition and
  definition-bound cases passed 278/278. Ruff passed for all WP1–WP3
  production and focused-test files. A supplemental feature-width, annotation,
  track-slot, feature-path, and legend-measurement group passed 206/206 after
  its private-helper mocks and transform readers were migrated to the direct
  layout/outer-composition contracts.
- The regenerated T-PY-06 acceptance recipe measured a 24.05 px visible right
  gap, 1734.51 px final width, 0.33 px title-center offset, 390 px radius, and
  no clipped paint. Its three rings, 106 homology hits, record identity, and
  title preset were retained.
- The same SVG converted successfully to PNG (1735 x 995), PDF, EPS, and PS;
  the raster result was visually inspected without edge clipping. Tracked
  reference outputs have no substantive diff from the recorded WP0 state.

Edit:

- `gbdraw/canvas/circular.py`;
- `gbdraw/diagrams/circular/assemble.py`;
- `gbdraw/diagrams/circular/positioning.py`;
- `gbdraw/diagrams/circular/builders.py`;
- `gbdraw/configurators/legend.py`;
- relevant Circular tests.

Changes:

1. Build record/track/definition geometry in plot space.
2. Convert `CircularRadialLayout` plus exact external labels/leaders into one
   primary `Aabb` including stroke safety.
3. Make external-label fitting return bounds instead of expanding the final
   canvas.
4. Measure and place title/legend only after primary bounds are final.
5. Use the shared composer for docked legends and title lanes.
6. Restrict the existing label–legend solver to overlay placements.
7. Center the title on primary content.
8. Apply one final composition plan and attach internal metadata.
9. Keep batch behavior as repeated single-record composition.

Retire in this package:

- `resolve_circular_side_legend_geometry`;
- side branches in `CircularCanvasConfigurator.recalculate_canvas_dimensions`;
- docked use of `_try_shift_labels_away_from_legend`,
  `_try_move_legend_away_from_labels`, and
  `_expand_canvas_for_legend`;
- `_expand_canvas_to_fit_external_labels` as a final-canvas mutation;
- duplicate viewBox synchronization/whole-canvas translation;
- separate single-record top/bottom placement policy.

Preserve and prove:

- radius, feature/track bands, ticks, labels, leaders, definitions, and
  comparison rings keep their relative geometry;
- corner overlays still resolve collisions;
- title and legend each occur at most once;
- interactive and static paths use the same composition.

Gate:

- focused Circular geometry tests pass;
- T-PY-06 meets the quantitative acceptance budget in section 13;
- side changes do not change label positions or record-local geometry;
- rendered SVG, PNG, PDF, EPS, and PS samples are unclipped;
- WP2 has not modified tracked references beyond their recorded WP0 state.

### WP3 — Migrate Circular multi-record composition and delete SVG estimation

Status: [x] Complete

Execution evidence (2026-08-04):

- Multi-record packing consumes each record assembly result's direct visible
  bounds, applies row/column ratios between those bounds, and composes the grid,
  shared legend, and shared title through one global plan. A one-record grid
  follows the same record assembly contract.
- Tests cover unequal and zero gaps, explicit grids, all nine legend positions,
  one shared title/legend, grid-title centering, record-local parity, duplicate
  IDs, and URL/href reference integrity. The combined Circular gate passed
  278/278; an additional planner/renderer/multi-record selection passed
  203/203 during implementation.
- The five SVG estimators and their serialize/reparse, centering, translation,
  multi-only title, and side-legend branches are deleted. Searches confirm no
  SVG serialization or XML reparsing remains in the grid layout path.

Edit:

- `gbdraw/api/diagram.py`;
- `gbdraw/diagrams/circular/assemble.py`;
- `gbdraw/layout/circular.py` as needed for the internal result;
- `tests/test_circular_multi_canvas.py`;
- `tests/test_circular_svg_id_integrity.py`;
- `tests/test_circular_render_context.py`.

Changes:

1. Add an internal Circular record-assembly result carrying authoritative
   content bounds and existing track-slot geometry.
2. Make the grid consume those bounds directly.
3. Pack record gaps between visible bounds while preserving
   `multi_record_column_gap_ratio` and `multi_record_row_gap_ratio`.
4. Compose the grid as the global primary item, then add one shared legend and
   title through the common planner.
5. Make a one-record grid use the same record assembly and composition
   contracts as a single record.
6. Keep record definitions centered on each record axis and preserve mixed-size
   scaling/style harmonization.
7. Preserve ID uniquification when copying subtrees.

Delete after parity is proved:

- `_estimate_element_local_x_bounds`;
- `_estimate_element_local_y_bounds`;
- `_estimate_subcanvas_content_bounds`;
- `_estimate_subcanvas_horizontal_insets`;
- `_estimate_subcanvas_vertical_insets`;
- multi-only side legend spacing branches/constants;
- multi-only title/legend outer-canvas calculations superseded by the composer.

Tests must cover:

- visible row/column gap formulas;
- unequal record sizes and labels;
- explicit record rows/columns;
- zero gap ratios;
- shared legend/title once;
- every legend position;
- duplicate record IDs and reference integrity;
- one-record grid/single record-local parity;
- title centered on grid primary bounds.

Gate:

- no SVG serialization/reparse occurs to lay out a grid;
- current scale-mode and row-placement tests pass;
- direct-bounds tests replace estimator-specific assumptions;
- no copied reference is dangling or duplicated.

### WP4 — Migrate Linear final composition

Status: [x] Complete

Execution evidence (2026-08-04):

- Linear record-row, shared-scale, track-slot, collision-band, ruler, label,
  definition, and comparison geometry remains mode-owned and is collected into
  one structured primary `Aabb`; one shared plan now composes primary targets,
  measured legend, and title and emits schema-1 metadata.
- The full plan-listed Linear group plus label, annotation, and new exact
  composition cases passed 232/232 after the final removal of the fixed
  `2 * canvas_padding` staging offset. The 12 composition cases prove all dock
  gaps, title modes, same-side stacking, metadata, record-local invariance, and
  final track-slot translations. Ruff passed.
- Browser geometry auditing of a dense rotated-label/ruler example and a
  multi-row top-stack example found 0/2 clipped files. The side-legend visible
  gap was 28.45 px (within the 32 px browser tolerance); planner tests remain
  exactly 24 px.
- A representative composed Linear SVG converted to PNG (2366 x 386), PDF,
  EPS, and PS with non-empty recognized artifacts. Dense side-legend and
  multi-row top-stack rasters were visually inspected without clipping or
  record-order drift. Tracked references remain untouched from the WP0 state.

Edit:

- `gbdraw/canvas/linear.py`;
- `gbdraw/diagrams/linear/assemble.py`;
- `gbdraw/diagrams/linear/builders.py`;
- `gbdraw/layout/linear.py`;
- `gbdraw/layout/linear_multi_record.py` only if the result needs an explicit
  primary bounds field;
- Linear layout/legend tests.

Changes:

1. Keep existing record-row, shared-scale, track-slot, collision-band, and
   comparison-corridor solvers.
2. Collect their resolved horizontal and vertical painted extents into one
   primary `Aabb`.
3. Include definitions, labels, ruler/scale, annotations, and stroke extents.
4. Pass primary, title, and measured legend to the common composer.
5. Remove the fixed `2 * canvas_padding` side gap.
6. Center side legends orthogonally on primary content and titles on primary
   width.
7. Preserve `center` title as an overlay.
8. Preserve existing relative order when bottom/top title and legend share a
   side: legend nearest plot, title outward.
9. Emit the same composition metadata contract as Circular.

Retire:

- legend-dependent final sizing in
  `LinearCanvasConfigurator.recalculate_canvas_dimensions`;
- separate title/legend final-height branches;
- duplicate final width extension and viewBox updates;
- new-output dependence on unversioned alternate-viewBox calculations.

Tests must cover:

- one record and multi-row records;
- above/middle/below track layouts;
- rotated/external labels;
- definitions on left and record-local headers;
- scale bar, ruler, and ruler-on-axis;
- quantitative and annotation tracks;
- pairwise, orthogroup, and collinear comparisons;
- top/bottom/left/right/none legend;
- top/bottom/center title;
- same-side legend/title stacks;
- normalized-length behavior where supported.

Gate:

- existing Linear record/track solvers remain semantically unchanged;
- all dock gaps come from `CompositionSpacing`;
- no clipping appears at either axis end or in definitions/labels;
- WP4 has not modified tracked references beyond their recorded WP0 state.

### WP5 — Converge Web editor layout on Python metadata

Status: [~] In progress

Edit/add:

- `gbdraw/render/composition_metadata.py` or an equivalently focused owner;
- `gbdraw/web/js/app/legend-layout/canvas-actions.js`;
- `gbdraw/web/js/app/legend-layout/reposition-actions.js`;
- `gbdraw/web/js/app/legend-layout/transform-utils.js`;
- `gbdraw/web/js/app/legend/layout-actions.js`;
- `gbdraw/web/js/app/legend/drag-actions.js`;
- `gbdraw/web/js/app/legend-layout/diagram-drag.js`;
- `gbdraw/web/js/app/watchers.js`;
- `gbdraw/web/js/app/app-setup.js`;
- `gbdraw/web/js/services/svg-sanitization.js`;
- `gbdraw/web/js/services/svg-serialization.js`;
- new focused Node/browser tests.

Changes:

1. Parse and validate composition schema 1 at preview/session load.
2. Replace Circular percentages, Linear fixed padding, and inferred base
   viewBoxes with the metadata-driven interpreter.
3. Keep browser legend reflow for post-generation entry edits, but source all
   spacing and placement constraints from metadata.
4. Keep palette/comparison-legend updates separate from geometry so a layout
   change cannot regress colors.
5. Preserve user legend and diagram deltas across reflow, position changes,
   history, session serialization, and export.
6. Fix transform composition so scale/rotate/matrix/chained transforms are not
   discarded.
7. Make Reset serialize in both legend and no-legend diagrams.
8. Add the explicit one-time adapter for old SVGs and mark normalized state.
9. Retire the current silent missing-metadata reconstruction for current SVGs.
10. Retire `data-horizontal-viewbox`/`data-vertical-viewbox` from new output
    once every current consumer uses schema 1; retain reader support only in the
    legacy adapter.

Required functional test:

1. render a Linear and a Circular figure;
2. switch legend orientation/side;
3. add or rename a legend entry so dimensions change;
4. move legend, diagram, and title;
5. save, reload, undo/redo, reset, and export;
6. assert positions, user deltas, colors, title alignment, viewBox, and
   sanitization after each boundary.

Gate:

- equivalent fresh Python and Web requests have the same automatic
  composition;
- current SVGs never enter a heuristic fallback;
- old supported sessions remain editable through the explicit adapter;
- browser work remains offline and no new runtime dependency is added.

Execution evidence (updated 2026-08-05):

- added a strict schema-1 browser interpreter and explicit legacy adapter;
  removed current-output percentage/fixed-padding/viewBox heuristics while
  preserving independent primary, legend, and title deltas and chained SVG
  transforms;
- made Python the owner of `legendReflow` and the four-field `overlayPolicy`;
  current metadata rejects incomplete or invalid numeric policy values, while
  only schema-free legacy input receives the frozen v0 policy;
- all 13 focused Node commands pass. The Python-to-JavaScript oracle also
  matches dock, overlay, title-stack, no-legend, and canvas-growth plans in both
  runtimes;
- session import now rolls back atomically across 18 owner snapshots, 54
  transient fields, and four app-local fields. Delayed Circular and Linear
  reparse failures restore scroll, depth-track, selection, and file state
  without reparsing the restored files;
- the required real-app Circular/Linear Playwright journey is implemented and
  passes static syntax/selector checks. Its local loopback/headless Chromium
  execution has not run because explicit approval is still pending;
- targeted Web packaging checks passed 11/11, and the sanitizer allowlists
  exactly the schema-1 composition attributes;
- corrected the pre-existing offline verifier's stale eager-main-Pyodide
  readiness condition to the current lazy-main/ready-worker contract; asset
  validation and the previously approved offline Circular/Linear smoke test
  passed;
- no new dependency or runtime URL was introduced. The
  session-definition-rehydration assertion still fails identically in the
  clean WP0 archive and is not caused by this layout work.

### WP6 — Verify API, CLI, request, session, and export contracts

Status: [x] Complete

Edit only if tests reveal a required correction:

- `gbdraw/api/options.py`;
- `gbdraw/api/requests.py`;
- `gbdraw/api/request_render.py`;
- `gbdraw/session_request_codec.py`;
- `gbdraw/circular.py`;
- `gbdraw/linear.py`;
- Web request/session projection modules;
- contract tests.

Changes:

1. Prove every existing legend/title value reaches the shared planner.
2. Prove Python API, CLI, Web generation, and session replay converge for an
   equivalent request.
3. Prove defaults stay surface-correct.
4. Prove no request/session field or public symbol was added.
5. Verify all export formats use the final viewBox and return real artifacts.

Gate:

- public contract fixture is byte-unchanged from the recorded WP0 state;
- live CLI option sets match `docs/CLI_Reference.md`;
- request schema remains 5 and session writer remains 40;
- current and supported legacy replay tests pass.

Execution evidence (updated 2026-08-05):

- added 38 focused composition-surface contract cases covering every documented
  Circular and Linear legend/title value through typed options, canonical
  schema-5 encode/decode, emitted schema-1 metadata, and fresh-versus-replay
  geometry equivalence;
- validated the exact nine-value Circular and five-value Linear legend sets at
  their existing public option owners without changing signatures or type
  annotations, including Linear CLI rejection of corner values before output;
- produced real SVG, PNG, PDF, EPS, and PS conversions from both diagram modes,
  and verified interactive SVG viewBox/schema preservation;
- the combined API/request/session/interface contract gate passed 625/625
  selected tests after excluding the clean-HEAD orthogroup-sidecar baseline
  failure, and Ruff passed;
- request schema 5, session writer 40, and the public-contract fixture remain
  substantively unchanged from WP0 (the fixture has no diff when EOL-only
  working-tree changes are ignored).

### WP7 — Refresh references, documentation figures, and Gallery assets

Status: [~] In progress

Execution evidence (updated 2026-08-05):

- regenerated all 16 reference SVGs through the owner test and re-ran the 16
  comparisons; an independent Chromium review found no clipping and confirmed
  that record-local paths and topology were preserved;
- regenerated the complete CLI/Python documentation inventory through its
  recipe owners and passed the exact `--check` runs. T-PY-06 now has a 24.05 px
  visible side gap, a 1734.51 px canvas width, a 0.33 px title-center offset,
  the unchanged 390 px radius, and no clipped content;
- refreshed all 11 Gallery sessions, source SVGs, examples, thumbnails, and
  catalog entries through `tools/refresh_gallery_sessions.py`; 86 focused
  Gallery tests and the six strict pre-capture checks passed;
- reproduced and visually approved all 57 manifest-owned public examples. Ten
  recipes remain correctly skipped because their declared source inputs are
  absent, and no available recipe failed. CairoSVG review of
  `NC_001879_color.svg` and `NC_001879_regions.svg` approved all four overlay
  corners without clipping or content overlap;
- all 142 final standalone generated SVGs contain the Python-owned
  `overlayPolicy`. Sixteen embedded SVG results across 17 session artifacts
  contain it as well; T-PY-08 intentionally has `results: []` and therefore no
  embedded composition metadata;
- the protected social preview remains byte-identical to HEAD at SHA-256
  `b5e5341fbd8b564626221e79e2c329ba559c1ea58fbb0e3ab65a2f8144d67ee2`;
- the first GUI owner pass replaced 18 result PNGs and left seven input/settings
  PNGs unchanged. Six outputs for T-GUI-01 and H-GUI-03 passed review. Corrected
  capture flows are ready for the other 12 PNGs, and six Gallery result WebPs
  still require the approved local Chromium run.

Do this only after WP1–WP6 gates pass and the geometry change is approved
visually.

#### Reference SVGs

Owner: `tests/test_output_comparison.py::TEST_CASES`. Workflow:

```bash
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python -m pytest tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs -v
git diff -- tests/reference_outputs/
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Review every changed SVG. Do not accept a bulk rewrite without classifying
viewBox/translation changes separately from record-local geometry.

#### Executable CLI/Python documentation figures

Use `docs/scenarios/manifest.json` as the inventory and the existing recipe
runners as owners:

```bash
python docs/recipes/run_cli_scenarios.py \
  --all --output-root /tmp/gbdraw-layout-cli-review
python docs/recipes/run_python_scenarios.py \
  --all --output-root /tmp/gbdraw-layout-python-review
```

After reviewing the temporary figures, regenerate through the same owners and
then prove exact reproducibility:

```bash
python docs/recipes/run_cli_scenarios.py --all
python docs/recipes/run_python_scenarios.py --all
python docs/recipes/run_cli_scenarios.py --all --check
python docs/recipes/run_python_scenarios.py --all --check
```

Compare at identical displayed size before replacement. T-PY-06 is mandatory
evidence.

#### GUI documentation captures

Audit the manifest before selecting all layout-sensitive scenarios. Current
known candidates include `T-GUI-01`, `T-GUI-02`, `T-GUI-05`, `H-GUI-02`,
`H-GUI-03`, `H-GUI-09`, `H-GUI-10`, and `H-GUI-11`.

Use the real capture flows:

```bash
python docs/capture/run_all.py --scenario <ID> --tier <tier>
python docs/capture/run_all.py --scenario <ID> --tier <tier> --check
```

#### Gallery sessions and rendered assets

The current owner is `tools/refresh_gallery_sessions.py`, which calls the
actual asset generator `tools/prepare_interactive_gallery_assets.py`:

```bash
python tools/refresh_gallery_sessions.py
python -m pytest tests/test_refresh_gallery_sessions.py -v
python -m pytest tests/test_gallery_session_semantics.py -v
python -m pytest tests/test_web_packaging.py -v
```

Do not use `--skip-session-refresh` after renderer geometry changes; it can
reuse stale SVG content from a session. Gallery sessions, source SVGs,
interactive examples, thumbnails, and `examples.json` are generator-owned and
must not be edited by hand.

#### Gallery tutorial screenshots

Create or update
`docs/WEB_GALLERY_DIAGRAM_LAYOUT_RECAPTURE_PLAN.md` before replacements and
record decisions in `docs/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`.
Capture from each exact example session:

```bash
python tools/capture_gallery_tutorial_screenshots.py \
  --example <example-id> --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example <example-id>
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
npx playwright test tests/web/gallery-session-regeneration.playwright.spec.js \
  --project=chromium
```

Compare old and new images side by side at the same displayed size. Preserve
labels, legend entries, color rules, metadata, quantitative tracks, and
comparison context.

#### Older public examples

Respect `tools/reproduce_examples_manifest.py` and its
`MANUALLY_MANAGED_FIGURES` list:

```bash
# Use the manifest-owned reproduction CLI to list and render <figure-id>
# into /tmp/gbdraw-layout-review.
python -m pytest tests/test_reproduce_examples.py -v
```

Do not route specialized browser screenshots, logo assets, the large
owner-managed Vibrio figure, or `examples/gbdraw_social_preview.png` through a
generic replacement.

Gate:

- generated assets reproduce from their declared source;
- every changed public figure is visually approved;
- no owner-managed asset changed;
- production, tests, docs, references, sessions, SVGs, thumbnails, and
  screenshots are reviewed as separate diff groups.

### WP8 — Remove superseded paths and run the full gate

Status: [~] In progress

Execution evidence (updated 2026-08-05):

- removed the retired Circular pre-positioning, Linear canvas-state, SVG
  estimator/reparse, current-output Web reconstruction, fixed side-inset, and
  retired viewBox-attribute paths; exact production searches for the named
  owners return no matches outside the explicit legacy adapter;
- updated the layout explanation and 0.14 release notes. The final focused plan
  gate passed 534/534 tests: 72 composition/legend, 274 Circular, and 188
  Linear. All 13 focused Node commands pass;
- added a deterministic three-case assembly benchmark. With three warmups, 21
  timing samples, and three memory samples, final HEAD-to-current medians are
  216.677074 to 206.170523 ms for single Circular (-4.848944%, 1.145225%
  noise, -36.167821% memory), 277.635626 to 238.481117 ms for multi Circular
  (-14.102840%, 0.820811% noise, +0.459556% memory), and 15.862921 to
  15.417283 ms for multi-record Linear (-2.809306%, 1.714044% noise,
  -1.475103% memory). All three pass the 5%/measured-noise gate;
- reference comparisons passed 16/16; Ruff, browser-wheel preparation, offline
  asset validation, and `python -m build` all pass;
- the final browser-free broad rerun selected 2,984 tests: 2,960 passed, 19
  skipped, and the same five clean-HEAD baseline tests failed; 14 browser or
  slow selections were deselected. No plan-attributed failure remains;
- a separate full packaging run reached 107 passing tests and the known
  clean-HEAD session-definition failure before entering a Playwright case. It
  was interrupted immediately, so no browser result is counted from that run.

The five clean-HEAD baseline failures are:

- `tests/test_public_contract.py::test_public_api_and_cli_contract_snapshot`;
- `tests/test_session_io.py::test_session_sidecar_saves_complete_orthogroup_state`;
- `tests/test_reproduce_examples.py::test_public_figures_have_reproduction_inventory_coverage`;
- `tests/test_web_request_render.py::test_embedded_web_request_preserves_in_memory_comparison_metadata`;
- `tests/test_web_packaging.py::test_web_session_definition_resource_rehydration`.

Changes:

1. Search for every retired helper and spacing constant named in this plan.
2. Remove unused imports, compatibility branches for current output, and
   duplicated tests that assert the retired nominal-width policy.
3. Keep only the explicit supported-old-SVG adapter.
4. Update user-facing explanation and release notes:
   - `docs/EXPLANATION/understand-tracks-axes-and-layout.md`;
   - `docs/RELEASE_NOTES_0.14.0b0.md`.
5. State that default outer composition is tighter and content-aware; do not
   imply a new user setting exists.
6. Run focused, broad, browser, offline, export, and visual checks.

Gate:

- the Definition of Done in section 16 is fully evidenced;
- no superseded layout owner remains;
- the final working tree contains only intentional files.

## 10. Test plan

### 10.1 Pure geometry

Add relation-based tests rather than asserting only final canvas sizes:

- dock gap;
- edge padding;
- union/translation;
- mirror symmetry;
- title and legend alignment;
- same-side stacking;
- overlay collision order;
- invalid/empty inputs;
- deterministic repeated planning.

### 10.2 Circular

Primary files:

- `tests/test_circular_label_placement.py`;
- `tests/test_circular_radial_labels.py`;
- `tests/test_circular_radial_layout.py`;
- `tests/test_circular_render_context.py`;
- `tests/test_circular_feature_width.py`;
- `tests/test_circular_conservation.py`;
- `tests/test_circular_annotation_tracks.py`;
- `tests/test_circular_multi_canvas.py`;
- `tests/test_circular_svg_id_integrity.py`;
- `tests/test_legend_measurement.py`.

Replace dock-specific label-shift expectations with the stronger invariant that
docked legend changes leave label geometry unchanged. Retain the label-shift
tests for overlay legends.

### 10.3 Linear

Primary files:

- `tests/test_linear_vertical_layout.py`;
- `tests/test_linear_multi_record_layout.py`;
- `tests/test_linear_multi_record_comparisons.py`;
- `tests/test_linear_track_layout.py`;
- `tests/test_linear_track_slots.py`;
- `tests/test_linear_definition_alignment.py`;
- `tests/test_linear_label_placement.py`;
- `tests/test_linear_annotation_tracks.py`.

Keep the existing record/track solver assertions. Update only outer-composition
expectations such as the fixed side gap or total viewBox.

### 10.4 API, CLI, request, and session

- `tests/test_api_requests.py`;
- `tests/test_api_request_render.py`;
- `tests/test_api_library_usage.py`;
- `tests/test_dead_api_cleanup.py`;
- `tests/test_session_request_codec.py`;
- `tests/test_api_session.py`;
- `tests/test_session_compat.py`;
- `tests/test_session_io.py`;
- CLI forwarding and reference-help tests.

### 10.5 Web and browser

Node:

- `tests/web/layout-preferences.test.mjs`;
- `tests/web/session-request.test.mjs`;
- `tests/web/session-authority.test.mjs`;
- `tests/web/legend-sync.test.mjs`;
- `tests/web/history.test.mjs`;
- `tests/web/svg-sanitization.test.mjs`;
- `tests/web/composition-layout.test.mjs`;
- `tests/web/legend-layout-actions.test.mjs`;
- `tests/web/composition-runtime-parity.test.mjs`;
- `tests/web/session-draft-authority.test.mjs`;
- `tests/web/run-analysis-simple-path.test.mjs`;
- `tests/web/record-selector.test.mjs`;
- `tests/web/depth-track-state.test.mjs`.

Browser:

- new Circular/Linear legend-edit, drag, reset, round-trip test;
- `tests/web/linear-multi-record.playwright.spec.js`;
- `tests/web/depth-track-session.playwright.spec.js`;
- `tests/web/gallery-session-regeneration.playwright.spec.js`;
- `tests/web/gallery-tutorial.playwright.spec.js`.

### 10.6 Output/reference tests

- `tests/test_output_comparison.py::TestOutputComparison`;
- SVG semantic hook and ID integrity tests;
- static and interactive SVG comparison;
- CairoSVG conversion tests where available.

## 11. Verification sequence

Run in this order so failures identify the responsible boundary.

### 11.1 Pure and focused Python

```bash
python -m pytest tests/test_composition_layout.py tests/test_legend_measurement.py -q
python -m pytest \
  tests/test_circular_label_placement.py \
  tests/test_circular_radial_labels.py \
  tests/test_circular_radial_layout.py \
  tests/test_circular_render_context.py \
  tests/test_circular_multi_canvas.py \
  tests/test_circular_conservation.py \
  tests/test_circular_svg_id_integrity.py -q
python -m pytest \
  tests/test_linear_vertical_layout.py \
  tests/test_linear_multi_record_layout.py \
  tests/test_linear_multi_record_comparisons.py \
  tests/test_linear_track_layout.py \
  tests/test_linear_track_slots.py \
  tests/test_linear_definition_alignment.py -q
```

### 11.2 Contracts

```bash
python -m pytest \
  tests/test_api_requests.py \
  tests/test_api_request_render.py \
  tests/test_api_library_usage.py \
  tests/test_dead_api_cleanup.py \
  tests/test_session_request_codec.py \
  tests/test_api_session.py \
  tests/test_session_compat.py \
  tests/test_session_io.py \
  tests/test_composition_surface_contracts.py \
  -q -k "not test_session_sidecar_saves_complete_orthogroup_state"
```

The excluded orthogroup-sidecar assertion fails unchanged in the recorded clean
HEAD archive and remains listed as a baseline exception rather than a layout
regression.

### 11.3 Node

```bash
node tests/web/layout-preferences.test.mjs
node tests/web/session-request.test.mjs
node tests/web/session-authority.test.mjs
node tests/web/legend-sync.test.mjs
node tests/web/history.test.mjs
node tests/web/svg-sanitization.test.mjs
node tests/web/composition-layout.test.mjs
node tests/web/legend-layout-actions.test.mjs
node tests/web/composition-runtime-parity.test.mjs
node tests/web/session-draft-authority.test.mjs
node tests/web/run-analysis-simple-path.test.mjs
node tests/web/record-selector.test.mjs
node tests/web/depth-track-state.test.mjs
```

### 11.4 Browser

Check both installations first:

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Then run the focused specs. If Node Playwright is absent, run an equivalent
Python Playwright check. If Chromium fails with a sandbox permission error,
rerun the same local check with the required sandbox escalation.

The required real-app composition journey is:

```bash
npx playwright test tests/web/composition-layout-real.playwright.spec.js \
  --project=chromium --workers=1
```

The spec passes static syntax and selector checks. Its local
loopback/headless-Chromium run remains pending explicit approval and is not
counted as passing evidence.

Run the corpus geometry audit against the regenerated documentation figures:

```bash
python tools/audit_svg_composition.py docs/images \
  --json-out /tmp/gbdraw-docs-composition.json
```

### 11.5 References, lint, broad suite, and packaging

```bash
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
python -m pytest tests/ -v -m "not slow"
python tools/prepare_browser_wheel.py
python -m pytest tests/test_web_packaging.py -v
python -m build
```

The targeted packaging and asset checks are complete. The full packaging run
entered a Playwright case after 107 passes and the known clean-HEAD
session-definition failure, so it was interrupted and is not a browser gate.

Run the assembly performance gate in separate processes for a clean HEAD tree
and the working tree, then compare the reports:

```bash
layout_head_dir=$(mktemp -d /tmp/gbdraw-layout-head.XXXXXX)
git archive HEAD | tar -x -C "$layout_head_dir"
python tools/benchmark_diagram_layout.py run \
  --source-root "$layout_head_dir" --label HEAD \
  --output /tmp/gbdraw-layout-head.json
python tools/benchmark_diagram_layout.py run \
  --source-root . --label working-tree \
  --output /tmp/gbdraw-layout-current.json
python tools/benchmark_diagram_layout.py compare \
  --baseline /tmp/gbdraw-layout-head.json \
  --current /tmp/gbdraw-layout-current.json \
  --output /tmp/gbdraw-layout-comparison.json
```

Prepare the generated browser wheel when needed, but never edit or commit it.
Refresh the cache-bust token only when preparing a deployable bundle.

## 12. Manual visual-QA matrix

Inspect SVG at a readable size and at least one raster/vector conversion for
each representative row:

| Case | Required checks |
| --- | --- |
| Circular, labels off, right/left legend | Exact compact gap, symmetry, unchanged radius |
| Circular, dense horizontal labels | No clipping; labels unchanged between dock sides |
| Circular, radial labels | Rotated text and leaders inside viewBox |
| Circular, wide gradient legend | Legend adjacent; title stays plot-centered |
| Circular, top/bottom legend and title | Correct stack order and even spacing |
| Circular corner overlays | Requested corner, no unresolved label collision |
| Circular one-record grid | Same record-local geometry as single |
| Circular mixed-size grid | Visible-content row/column gaps, centered shared title |
| Linear one record | Axis length unchanged; side gap compact |
| Linear rotated/external labels | No label or definition clipping |
| Linear multi-row comparisons | Record order, links, and rows unchanged |
| Linear title+legend on bottom/top | Legend nearest plot; title outward |
| Web legend edit | Reflow uses current measured legend without policy drift |
| Web drag/session/reset | User offsets round-trip and reset exactly |
| Interactive SVG | Search/popups and original viewBox remain correct |
| PNG/PDF/EPS/PS | Same content scale and no crop at every edge |
| Gallery desktop/mobile | Preview controls do not cover labels/title/legend |

For every public figure replacement, compare old and new at the same displayed
size. A tighter canvas is not acceptable if it reduces label readability,
removes useful context, or changes biological content.

## 13. Quantitative acceptance budgets

### 13.1 Motivating T-PY-06 figure

The regenerated
`docs/images/t-py-06/python_precomputed_circular_rings.svg` must satisfy:

- primary-content-to-right-legend visible gap: 24 px ± 1 px;
- final width: no more than 1800 px;
- plot-title center within 1 px of the primary-content horizontal center;
- Circular radius remains 390 px;
- three comparison rings, 106 hits, record identity, labels, legend entries,
  colors, and title text remain unchanged;
- no painted bounds cross the final viewBox.

### 13.2 Documentation corpus

For generated docked side legends, excluding explicitly manually dragged
artifacts:

| Metric | Target |
| --- | ---: |
| Circular median visible gap | ≤ 26 px |
| Circular maximum visible gap | ≤ 32 px |
| Linear median visible gap | ≤ 26 px |
| Linear maximum visible gap | ≤ 32 px |
| Negative gap/overlap | 0 files |
| Clipped primary/title/legend bounds | 0 files |

The small tolerance above the 24 px policy covers browser font/stroke bounding
differences; planner-level tests remain exact.

### 13.3 Geometry preservation

For equivalent pre/post inputs, compare record-local coordinates after removing
the final composition translation:

- Circular axis radius and slot bands: exact within floating-point tolerance;
- feature paths, labels, and leader endpoints: unchanged for dock positions;
- Linear record axis lengths and bp/px scale: unchanged;
- record rows/order/orientation and comparison endpoints: unchanged;
- semantic element counts and documented hooks: unchanged.

## 14. Generated-artifact impact ledger

| Artifact class | Owner | Expected impact | Update policy |
| --- | --- | --- | --- |
| `tests/reference_outputs/*.svg` | output comparison generator | Yes | Regenerate only in WP7 after review |
| `docs/images/**.svg` | scenario capture/recipe runners | Yes for layout-sensitive figures | Generate from declared scenario; never hand-edit |
| Gallery sessions | `tools/refresh_gallery_sessions.py` | Yes, saved result SVG changes | Refresh transactionally |
| Gallery `sources/`, `examples/`, `thumbnails/`, `examples.json` | Gallery asset generator | Yes | Generator-owned overwrite |
| Gallery tutorial WebP | capture metadata/tool | Only previews showing changed result | Exact-session recapture and side-by-side review |
| Older example figures | reproduce manifest | Audit | Respect specialized/manual owners |
| Browser wheel | `tools/prepare_browser_wheel.py` | Local verification only | Generated, gitignored, never commit |
| `examples/gbdraw_social_preview.png` | human owner | None | Never touch |

## 15. Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Font metrics differ between Python, browser, and CairoSVG | Use Python metrics for generation, browser metrics only for edited legend size, metadata-owned spacing, cross-renderer clipping tests |
| Stroke, marker, arc, or rotated-text extents are underestimated | Include half-stroke safety and exact oriented text/leader geometry; add edge-focused fixtures |
| A long title intersects a tall side legend | Measure/wrap title against primary width, then resolve the 2-D intersection by moving only the title outward |
| Overlay behavior regresses while dock logic improves | Keep overlay as a separate explicit solver with its current label-first semantics and dedicated tests |
| Multi-record layout changes record scale or gaps | Carry direct bounds beside the existing scale plan; compare record-local geometry and visible gaps |
| New metadata becomes another public schema | Mark it internal, version it independently, and do not add it to documented semantic hooks |
| Old sessions lack metadata | Normalize once at the compatibility boundary; forbid silent fallback for current metadata |
| Browser transform rewrite loses rotation/scale | Compose transforms or use stable transform components; test chained transforms |
| Layout and recoloring remain coupled in Web code | Split geometry refresh from legend/pairwise color updates before changing behavior |
| Reference regeneration hides record-local regressions | Compare transforms/viewBox separately from normalized record-local SVG geometry |
| Public screenshots become cramped | Review at rendered Gallery size on desktop/mobile and retain breathing room for floating controls |
| Performance regresses from extra measuring | Use already-computed structured geometry, avoid serialize/reparse, and record assembly time/memory before and after |
| Dirty generated assets obscure the intended diff | Inventory first; review production, tests, docs, references, Gallery assets, and screenshots separately |

Performance acceptance:

- no SVG serialization/parsing during Python layout;
- no material regression above 5% in representative single Circular, multi
  Circular, and multi-record Linear assembly median time, unless measurement
  noise is larger and documented;
- no repeated DOM full-tree scan on every drag frame.

## 16. Definition of Done

The overhaul is complete only when:

- one shared pure composer owns final outer layout for Circular and Linear;
- primary-content bounds are authoritative structured data, not recovered from
  finished SVG;
- all docked legend gaps and edge padding come from one spacing policy;
- Circular single, grid, and batch use the same record composition contract;
- Linear uses the same outer composer without changing its record/track solver;
- docked legend position changes do not alter record-local labels or geometry;
- titles are centered on primary content;
- Web current outputs consume versioned Python metadata, while old supported
  results use one explicit adapter;
- manual offsets survive edit/history/session/export and Reset serializes them
  correctly;
- public API, CLI, typed request, and session schemas are unchanged;
- SVG semantic hooks, IDs, interactive behavior, and every export format pass;
- T-PY-06 and the full docs audit meet section 13;
- superseded fixed-width, estimator, and browser-heuristic paths are removed;
- reference and public generated artifacts are regenerated only by their owners
  and visually approved;
- `examples/gbdraw_social_preview.png` and all other manually managed assets are
  unchanged;
- focused tests, full non-slow suite, Ruff, Web packaging, browser checks, and
  build pass;
- the plan-attributed diff, evaluated against the recorded WP0 dirty baseline,
  contains no unrelated files or generated distribution artifacts; unrelated
  pre-existing changes remain untouched.

The five enumerated clean-HEAD failures are baseline exceptions, not passing
results. The explicit-approval browser checks and captures also remain
outstanding, so the overall Definition of Done is not yet marked complete.

## 17. Execution evidence ledger

Update this table while implementing. “Complete” requires the named gate, not
only code presence.

| Work package | Status | Evidence |
| --- | --- | --- |
| WP0 baseline and requirements | Complete | Browser audit reproduced T-PY-06 and the 30 Circular / 12 Linear side-legend census; 256 Circular, 188 Linear, 441 contract, 16 reference, and six Node baseline checks passed; raw-byte artifact hashes recorded above |
| WP1 pure composition geometry | Complete | 82 pure/measurement and 266 affected-call-site tests passed; 16 reference comparisons and all WP0 artifact hashes remained unchanged; Ruff/diff checks passed |
| WP2 Circular single/batch | Complete | Circular gate 278/278 plus supplemental geometry group 206/206 and Ruff passed; T-PY-06 gap 24.05 px, width 1734.51 px, title offset 0.33 px, radius 390 px, no clipping; SVG/PNG/PDF/EPS/PS exports inspected |
| WP3 Circular multi-record | Complete | Direct visible-bounds packing and one global composition plan; combined Circular gate 278/278 and implementation selection 203/203 passed; estimator/reparse paths removed and ID/reference tests passed |
| WP4 Linear | Complete | Linear geometry/composition group 232/232 and Ruff passed; exact 24 px planner gaps and 20 px stacks; two browser-audited figures unclipped; SVG/PNG/PDF/EPS/PS artifacts and dense/multi-row rasters inspected |
| WP5 Web editor convergence | In progress | Strict schema-1 interpreter and legacy adapter; Python-owned reflow/overlay policy with Python-to-JavaScript runtime parity; atomic import rollback; 13 focused Node commands, packaging 11/11, offline assets, and the approved offline Circular/Linear smoke passed. The real-app Playwright journey is statically checked, but its Chromium run awaits approval |
| WP6 cross-surface contracts | Complete | 38 new surface-contract cases plus a 625/625 selected contract gate passed after excluding the clean-HEAD sidecar baseline; exact documented option sets, schema 5/session 40, unchanged public fixture, fresh/replay equivalence, interactive SVG, and SVG/PNG/PDF/EPS/PS exports verified |
| WP7 generated artifacts | In progress | All 16 references, the exact CLI/Python inventory, 11 Gallery asset sets, and 57 available public examples were owner-regenerated and reviewed; all 142 standalone and 16 embedded generated SVGs carry `overlayPolicy`. Six GUI outputs are accepted, while 12 corrected GUI PNGs and six Gallery WebPs await approved Chromium; the protected social preview is hash-identical |
| WP8 cleanup/full gate | In progress | Retired owners are removed; focused Python 534/534 and all 13 Node commands pass; final benchmark changes are -4.848944%, -14.102840%, and -2.809306%; references, Ruff, wheel/assets, and build pass. The browser-free broad rerun has 2,960 passes, 19 skips, and only the five named clean-HEAD failures; Chromium remains |

## 18. Delivery and coordination

- Implement in work-package order; WP2–WP5 may not update public generated
  artifacts independently.
- Keep each work package reviewable and record its focused evidence in section
  17.
- Inspect the working tree before every package and preserve unrelated changes.
- When geometry changes intentionally, do not call the implementation complete
  while tracked references or public generated figures still describe the old
  renderer.
- If the public API/session decision changes, amend sections 4 and 8 before
  implementation continues.
- Finish the implementation as one internally consistent repository change and
  provide an English commit title and short summary.

Proposed implementation commit title:

`layout: compose diagrams from visible content bounds`

Proposed summary:

- introduce one content-aware composition planner for Circular and Linear
  figures;
- replace nominal-canvas and SVG-reparse placement with authoritative bounds;
- converge Web legend editing on versioned render metadata and refresh
  geometry-dependent references and public figures.
