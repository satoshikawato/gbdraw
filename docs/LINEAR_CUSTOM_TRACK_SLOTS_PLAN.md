# Linear Custom Track Slots Plan

## Status

Planned.

## Background

Circular mode already has Custom Track Slots. The public model is implemented in
`gbdraw/tracks/circular.py`, routed through CLI/API options such as
`--circular_track_slot`, and exposed in the web UI through
`gbdraw/web/js/app/circular-track-slots.js`.

Linear mode does not currently have an equivalent model. Its plot-track stack is
implicit and fixed:

- `LinearCanvasConfigurator.set_gc_height_and_gc_padding()` in
  `gbdraw/canvas/linear.py` builds the stack in this order:
  depth tracks, GC content, GC skew.
- `gbdraw/diagrams/linear/positioning.py` positions each plot group by reading
  precomputed offsets from `LinearCanvasConfigurator`.
- `gbdraw/diagrams/linear/assemble.py` renders depth, GC content, and GC skew
  in fixed conditional blocks after the record feature/definition groups.
- `SeqRecordGroup` in `gbdraw/render/groups/linear/seq_record.py` draws the
  linear axis, optional ruler ticks, features, leader lines, and labels as one
  record group.

That design is simple and stable, but it does not let users reorder, duplicate,
hide, or place linear tracks above/below the axis in the same explicit way that
Circular Custom Track Slots allow.

## Goal

Add Linear Custom Track Slots to the Python CLI/API and the no-build-step web
GUI, while preserving existing linear output when the feature is not used.

The target user should be able to:

- enable a Custom Track Slots editor in Linear mode in the web GUI
- reorder linear tracks relative to the axis
- move supported tracks above or below the axis
- add spacer rows to control vertical separation
- duplicate supported numeric tracks where the renderer supports it
- target multiple depth tracks by `track_index`
- save/load web settings with linear slot state
- generate the same CLI args from the web GUI that a CLI user could type

## Design Principles

- **SOLID**: Separate slot parsing, slot normalization, linear stack layout,
  SVG group drawing, CLI/API plumbing, and web UI state. Do not make
  `LinearCanvasConfigurator` parse public slot strings or make drawer classes
  know about CLI syntax.
- **KISS**: Linear slots should solve the concrete linear layout problem:
  ordering and vertical placement around a fixed axis. Do not introduce a
  generic plugin framework, arbitrary x/y coordinates, or pairwise-ribbon slot
  placement in the first implementation.
- **DRY**: Reuse existing helpers where they fit, especially `ScalarSpec`,
  `parse_bool`, `split_kv_list`, depth-track normalization, GC/skew data
  precomputation, and existing SVG group builders. Do not share Circular
  radius-specific logic with Linear just to avoid a small amount of duplicated
  validation code.
- **Backward compatibility**: No SVG reference output should change when
  `linear_track_slots` are not supplied and the web Custom Track Slots toggle is
  off.

## Non-Goals

- Do not make the linear axis removable. In the UI it can appear as an Axis row,
  but it remains a fixed coordinate baseline.
- Do not put BLAST/protein/collinear pairwise comparison ribbons into the first
  slot model. Ribbons depend on spaces between records, not a per-record track
  stack.
- Do not add arbitrary absolute y-coordinate placement.
- Do not make Circular and Linear slots one shared dataclass in the first pass.
  Circular slots are radial and use `radius`/inside/outside semantics; Linear
  slots are vertical and use `height`/above/below semantics.
- Do not migrate old circular saved-session schema under the linear work. Linear
  saved state should get its own schema version.

## User-Facing Model

### Slot Concept

A linear slot is one ordered row in a per-record track stack. The axis is fixed.
Rows before an explicit axis boundary render above the axis. Rows at or after
the axis boundary render below the axis, unless an individual slot explicitly
sets `side`.

The default linear slot list should preserve the current implicit stack:

```text
features
axis
depth track(s)
gc_content
gc_skew
```

In this model, `features` is a supported slot because users expect it to move
relative to numeric tracks. The axis itself is not stored as a slot; it is a
boundary index, mirroring the circular web UI's Axis row approach.

### Supported Renderers

Initial supported renderers:

- `features`
- `dinucleotide_content`
- `dinucleotide_skew`
- `depth`
- `spacer`

Aliases:

- `gc_content` -> `dinucleotide_content`
- `content` -> `dinucleotide_content`
- `gc_skew` -> `dinucleotide_skew`
- `skew` -> `dinucleotide_skew`

Unsupported in the first pass:

- `ticks`: linear ruler ticks are part of the axis/ruler model.
- `pairwise_match`: comparison ribbons are inter-record geometry.
- `definition`: definitions belong in the left column and must stay aligned to
  record axes.
- `legend`, `title`, `length_bar`: global layout elements, not per-record
  tracks.

### CLI Options

Add to `gbdraw linear`:

```text
--linear_track_slot <slot_id>:<renderer>@key=value,key=value
--linear_track_axis_index INT
--linear_track_order features,depth,gc_content,gc_skew
```

`--linear_track_order` is a convenience shortcut for built-in slot ids. It
cannot be combined with `--linear_track_slot`.

Slot examples:

```text
--linear_track_slot features:features@side=above
--linear_track_slot depth_1:depth@track_index=0,h=35px,spacing=8px
--linear_track_slot gc_content:dinucleotide_content@nt=GC,side=below
--linear_track_slot gc_skew:dinucleotide_skew@nt=GC,h=28px
--linear_track_slot gap:spacer@h=12px
--linear_track_axis_index 1
```

Slot-level layout keys:

- `side=above|below|overlay`
- `h` or `height`: positive `ScalarSpec`; for Linear v1 only `px` and unitless
  px should be accepted.
- `spacing`: nonnegative `ScalarSpec`; for Linear v1 only `px` and unitless px
  should be accepted.
- `z`: integer draw order override.
- `enabled`, `show`, or `visible`: boolean.

Renderer params:

- `features`: no v1 renderer params beyond common layout fields.
- `dinucleotide_content`: `nt`, `dinucleotide`, `legend_label`.
- `dinucleotide_skew`: `nt`, `dinucleotide`, `legend_label`.
- `depth`: `track_index`, `legend_label`.
- `spacer`: no renderer params.

Validation:

- Duplicate slot ids are rejected.
- Unknown renderers are rejected.
- Generic layout fields inside `params` are rejected, matching the Circular
  Custom Track Slots schema style.
- `linear_track_axis_index` must be between 0 and number of supplied slots.
- `side=overlay` is allowed only for `features` in v1, and means the feature
  track is centered on the axis using the existing middle-layout semantics.
- `depth` slots require a depth input at generation time.
- Repeated `depth` slots must resolve to valid available depth track indexes.
- Empty slot lists are rejected when explicitly supplied.

## Public Python API

Add a new public track module:

```python
gbdraw/tracks/linear.py
```

Exports:

- `LinearTrackSlot`
- `NormalizedLinearTrackSlot`
- `parse_linear_track_slot()`
- `parse_linear_track_slots()`
- `normalize_linear_track_slots()`
- `normalize_linear_track_slots_with_axis()`
- `linear_track_slots_from_order()`
- `default_linear_track_slots()`

Update:

- `gbdraw/tracks/__init__.py`
- `gbdraw/api/tracks.py`
- `gbdraw/api/__init__.py`
- `gbdraw/api/options.py::TrackOptions`
- `gbdraw/api/diagram.py::assemble_linear_diagram_from_records()`

`TrackOptions` should become mode-neutral:

```python
@dataclass(frozen=True)
class TrackOptions:
    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None
    circular_track_axis_index: int | None = None
    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None
    center_reserved_radius: float | None = None
```

Keep circular fields unchanged.

## Internal Architecture

### 1. Parsing And Normalization

Implement `gbdraw/tracks/linear.py` with a deliberately small API similar to
`gbdraw/tracks/circular.py`, but linear-specific names and field semantics.

Suggested dataclasses:

```python
@dataclass(frozen=True)
class LinearTrackSlot:
    id: str
    renderer: str
    enabled: bool = True
    side: str | None = None
    height: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z: int = 0
    params: Mapping[str, Any] = field(default_factory=dict)

@dataclass(frozen=True)
class NormalizedLinearTrackSlot:
    slot_index: int
    id: str
    renderer: str
    enabled: bool
    side: str
    height: ScalarSpec | None
    spacing: ScalarSpec | None
    z: int
    reserve: bool
    params: Mapping[str, Any]
```

Naming note: use `height`, not `width`, because the linear track's variable
dimension is vertical. The CLI alias `h` keeps slot specs concise.

### 2. Linear Slot Layout Resolver

Add a focused resolver, likely:

```text
gbdraw/diagrams/linear/track_slots.py
```

Responsibilities:

- convert normalized slots into concrete per-renderer track plans
- resolve inherited default heights and spacing from `GbdrawConfig`
- calculate y offsets relative to the record axis
- calculate top and bottom extents needed by each record
- produce compatibility fields currently stored on `LinearCanvasConfigurator`
  such as `depth_track_offsets`, `gc_content_track_offset`,
  `gc_skew_track_offset`, `plot_tracks_height`, and
  `plot_tracks_visual_bottom`

Suggested structures:

```python
@dataclass(frozen=True)
class LinearResolvedTrack:
    slot_index: int
    id: str
    renderer: str
    side: str
    y_offset: float
    height: float
    top_extent: float
    bottom_extent: float
    z: int
    params: Mapping[str, Any]

@dataclass(frozen=True)
class LinearTrackLayout:
    slots: tuple[LinearResolvedTrack, ...]
    top_extent: float
    bottom_extent: float
    plot_tracks_height: float
    plot_tracks_visual_bottom: float
    depth_track_offsets: tuple[float, ...]
    gc_content_track_offset: float
    gc_skew_track_offset: float
```

Resolver rules:

- `below` tracks stack downward from the axis, preserving current defaults.
- `above` tracks stack upward from the axis.
- `features` uses existing feature positioning factors and should not be
  treated like a scalar plot track.
- `overlay` is for feature middle/on-axis behavior.
- `spacer` consumes height and spacing but produces no SVG group.
- `z` affects draw order only, not physical order.

KISS constraint: in v1, the physical stack order is list order plus side. If a
slot specifies a side conflicting with `linear_track_axis_index`, explicit side
wins after validation. The web UI should normally derive sides from the Axis
boundary to avoid conflicts.

### 3. Canvas Configurator Integration

Update `LinearCanvasConfigurator` so it can be configured by either:

- legacy implicit plot tracks, or
- a resolved `LinearTrackLayout`

Avoid making the canvas configurator parse slot specs. It should receive a
resolved layout or normalized simple data from the API/assembly layer.

Target behavior:

- When no custom slots are supplied, `set_gc_height_and_gc_padding()` remains
  the compatibility path.
- When custom slots are supplied, call a new method such as
  `set_linear_track_layout(layout)` that populates the same offset fields used
  by existing positioning helpers.
- Existing `calculate_dimensions()` and later assembly spacing should use
  `linear_track_layout.top_extent` and `.bottom_extent`, not only the legacy
  one-sided `plot_tracks_height`.

### 4. Feature Track Handling

`SeqRecordGroup` currently draws axis, features, labels, and leader lines in one
group. Linear Custom Track Slots need more control over feature placement and
draw order.

Recommended staged refactor:

1. Add optional parameters to `SeqRecordGroup` for a resolved feature side:
   `above`, `below`, or `overlay`.
2. Route that to existing `FeatureDrawer` through `track_layout`:
   - `above` -> existing above layout
   - `below` -> existing below layout
   - `overlay` -> existing middle layout
3. Keep axis drawing in `SeqRecordGroup` initially to minimize blast radius.
4. If z-order tests require it, extract axis/ruler drawing into a separate
   `AxisGroup` only after the first renderer-order tests fail for a concrete
   reason.

This keeps the first implementation smaller while still letting the feature
track move relative to scalar tracks.

### 5. Assembly Changes

Update `gbdraw/diagrams/linear/assemble.py`.

New high-level flow:

1. Normalize `linear_track_slots` if provided.
2. Determine which renderers are enabled from slots.
3. Derive `show_depth`, `show_gc`, and `show_skew` from slot renderers when
   custom slots are active.
4. Build GC/skew/depth data frames only for renderers that need them.
5. Resolve the linear slot layout before record offsets are calculated.
6. Use layout top/bottom extents in:
   - first record top margin calculation
   - inter-record spacing
   - final canvas height
   - comparison ribbon anchors
7. Draw slot-controlled groups by iterating resolved slots sorted by
   `(z, slot_index)`.

Existing conditional draw blocks:

```python
if depth_enabled:
    ...
if canvas_config.show_gc:
    ...
if canvas_config.show_skew:
    ...
```

should be replaced in custom-slot mode by a slot renderer dispatch:

```python
for slot in sorted(layout.slots, key=lambda s: (s.z, s.slot_index)):
    if slot.renderer == "spacer":
        continue
    if slot.renderer == "features":
        add_record_group(...)
    elif slot.renderer == "depth":
        add_depth_group(...)
    elif slot.renderer == "dinucleotide_content":
        add_gc_content_group(...)
    elif slot.renderer == "dinucleotide_skew":
        add_gc_skew_group(...)
```

To preserve compatibility, the legacy non-custom path can remain until the
custom path is stable. If both paths become too duplicated, collapse them by
generating default slots internally and using one renderer dispatch for all
linear diagrams.

### 6. Positioning Helpers

Update `gbdraw/diagrams/linear/positioning.py` carefully.

Preferred approach:

- keep existing `position_depth_group`, `position_gc_content_group`, and
  `position_gc_skew_group` for the legacy path
- add a generic helper:

```python
def position_linear_track_group(
    group: Group,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    track_offset_y: float,
) -> Group:
    ...
```

Then custom-slot rendering can use the generic helper with the resolved
`slot.y_offset`.

This avoids overloading the existing helpers with custom-only parameters.

### 7. Legend Behavior

Keep current legend semantics initially:

- Feature legend is driven by actual features present.
- GC/skew legend entries appear when corresponding renderers are enabled.
- Depth legend entries are synchronized with `record_depth_tracks`.
- `legend_label` slot params should be supported for depth and GC/skew where
  existing legend table code can accept it without large rewrites.

If legend label override requires invasive changes, implement the slot renderer
placement first and track custom legend captions as a follow-up. The parser can
accept `legend_label` from the start without using it immediately.

## Web GUI Plan

### State

Add to `gbdraw/web/js/state.js`:

```javascript
linear_track_slots_enabled: false,
linear_track_slots_schema_version: 1,
linear_track_slots_axis_index: null,
linear_track_slots: createDefaultLinearTrackSlots()
```

Use a separate schema version from circular slots.

### Module

Add:

```text
gbdraw/web/js/app/linear-track-slots.js
```

Responsibilities:

- default slot creation
- normalization
- Axis boundary handling
- move up/down within side
- move above/below across Axis
- add/duplicate/remove slots
- renderer label mapping
- build CLI slot specs
- depth slot synchronization with uploaded depth files

Do not copy the entire circular module. Reuse only small generic concepts:

- `cloneParams`
- optional text normalization
- renderer label helpers
- axis-index clamping
- slot spec builder pattern

If sharing those helpers reduces duplication without coupling circular and
linear semantics, extract them to a small module such as:

```text
gbdraw/web/js/app/track-slot-utils.js
```

### UI

Update the Linear mode branch in `gbdraw/web/index.html`, near the existing
Linear layout controls.

Controls:

- checkbox: `Custom Track Slots`
- reset button: rebuild from current simple Linear controls
- compact Axis stack editor:
  - rows above Axis = above-axis tracks
  - rows below Axis = below-axis tracks
  - icon buttons for move up/down, move above, move below, duplicate, delete
- per-row fields:
  - enabled checkbox
  - slot id
  - renderer select
  - height
  - spacing
  - z
  - `track_index` for depth
  - `nt` for GC content/skew
  - legend label field where accepted

The UI should stay dense and operational, matching the existing web app style.
Do not add a landing-page or explanatory marketing section.

### Run Analysis

Update `gbdraw/web/js/app/run-analysis.js`.

When mode is `linear` and `adv.linear_track_slots_enabled === true`:

- validate wheel support for:
  - `--linear_track_slot`
  - `--linear_track_axis_index`
- suppress conflicting legacy simple layout args if needed
- append `--linear_track_axis_index`
- append one `--linear_track_slot` per enabled slot
- ensure depth files are present when any enabled depth slot exists
- assign missing depth `track_index` values in slot order

Keep existing simple Linear controls working when custom slots are disabled.

### Config Save/Load

Update `gbdraw/web/js/services/config.js`:

- include linear slot state in saved sessions
- validate `linear_track_slots_schema_version`
- reject obsolete linear slot shapes if the schema changes later
- normalize imported slots before writing them into state

Do not reuse circular schema validation messages for linear slots; the error
should name `Linear Custom Track Slots`.

### Watchers

Update `gbdraw/web/js/app/app-setup.js` or `watchers.js`:

- when depth file count changes and linear slots are enabled, ensure enough
  depth slots exist or keep explicit disabled slots untouched
- when `show_gc`, `show_skew`, or `show_depth` changes with slots disabled,
  update the default template used when slots are enabled
- when slots are enabled, do not silently rewrite user-authored ordering except
  for normalization and depth-track synchronization

## CLI/API Plumbing Details

### `gbdraw/linear.py`

Add parser options near existing track layout/depth/GC controls:

```python
parser.add_argument("--linear_track_order", ...)
parser.add_argument("--linear_track_slot", action="append", ...)
parser.add_argument("--linear_track_axis_index", type=int, ...)
```

Validation:

- `--linear_track_order` cannot be combined with `--linear_track_slot`.
- `--linear_track_axis_index` requires order or slots.
- normalize once in CLI validation to catch errors early.
- pass raw specs to API, mirroring circular CLI behavior, unless early
  normalization is needed for validation.

Forward to `assemble_linear_diagram_from_records()`:

```python
linear_track_slots=linear_track_slot_specs,
linear_track_axis_index=args.linear_track_axis_index,
```

### `gbdraw/api/diagram.py`

Add helpers parallel to circular:

- `_parse_linear_track_slot_inputs()`
- `_validate_linear_track_axis_index()`
- `_linear_slots_have_renderer()`
- `_linear_slots_define_renderer()`
- `_dinucleotides_from_linear_slots()`

Adjust `assemble_linear_diagram_from_records()`:

- new parameters `linear_track_slots` and `linear_track_axis_index`
- parse and validate slots
- derive `show_depth`, `show_gc`, and `show_skew` from slots when supplied
- pass parsed slots and axis index into `assemble_linear_diagram()`

## Testing Plan

### Unit Tests

Add `tests/test_linear_track_slots.py`.

Coverage:

- parse slot spec with layout fields
- parse renderer aliases
- reject duplicate ids
- reject unknown renderer
- reject generic layout keys inside params
- normalize sides from `linear_track_axis_index`
- build default slots for combinations of `show_depth`, depth track count,
  `show_gc`, and `show_skew`
- parse `linear_track_order`
- reject invalid axis indexes
- preserve disabled slots by omitting them from normalized output

### API Rendering Tests

Use small in-memory `SeqRecord` objects where possible.

Cases:

- default linear slots produce groups equivalent to existing implicit ordering
  for group ids and transform y positions
- `gc_skew` above axis and `gc_content` below axis
- spacer increases inter-track or inter-record distance
- features below axis with GC above axis
- repeated depth slots select the correct `track_index`
- missing depth input with depth slot raises `ValidationError`
- custom slots disabled path still matches existing reference outputs

### CLI Tests

Add to existing linear tests or `tests/test_linear_track_slots.py`:

- `--linear_track_slot` forwards raw specs to API
- `--linear_track_axis_index` forwards value
- `--linear_track_order` expands/forwards slots
- conflicts are rejected:
  - order plus explicit slots
  - axis index without slots/order
  - invalid renderer

### Web Packaging Tests

Update `tests/test_web_packaging.py`:

- state contains `linear_track_slots_enabled`
- `run-analysis.js` wires `--linear_track_slot`
- config save/load validates linear slot schema
- UI includes Linear Custom Track Slots controls
- JS module can:
  - create default linear slots
  - move slots across Axis
  - build valid CLI specs
  - materialize multiple depth slots from uploaded depth files

### Regression / Reference Tests

Run targeted tests first:

```bash
python -m pytest tests/test_linear_track_layout.py tests/test_depth_track.py tests/test_gc_content_percent.py -q
python -m pytest tests/test_web_packaging.py -q
```

Then run the broader non-slow suite:

```bash
python -m pytest tests/ -m "not slow" -q
```

Reference SVGs should not change unless tests are intentionally added for a new
custom-slot output fixture.

## Implementation Phases

### Phase 1: Python Slot Model

Deliverables:

- `gbdraw/tracks/linear.py`
- exports through `gbdraw/tracks`, `gbdraw/api/tracks`, and `gbdraw/api`
- unit tests for parsing and normalization

Exit criteria:

- linear slot specs parse independently of diagram rendering
- invalid inputs fail with clear errors
- circular tests still pass

### Phase 2: Layout Resolver

Deliverables:

- `gbdraw/diagrams/linear/track_slots.py`
- resolved layout data structures
- compatibility adapter to populate current canvas offset fields

Exit criteria:

- resolver can reproduce current depth/GC/skew offsets from default slots
- resolver reports correct top/bottom extents for above and below stacks

### Phase 3: Assembly Integration

Deliverables:

- API parameters and CLI plumbing
- `assemble_linear_diagram()` accepts parsed slots and axis index
- custom-slot rendering dispatch for supported renderers

Exit criteria:

- API rendering tests pass
- existing reference output tests remain unchanged for non-custom mode

### Phase 4: Web State And CLI Serialization

Deliverables:

- `linear-track-slots.js`
- state additions
- run-analysis argument generation
- config save/load support

Exit criteria:

- web packaging tests prove arguments are wired
- no circular custom slot behavior regresses

### Phase 5: Web UI

Deliverables:

- Linear Custom Track Slots editor in `index.html`
- depth slot sync watchers
- reset/default behavior

Exit criteria:

- editor can create the same default stack as current Linear simple controls
- moving tracks above/below Axis changes generated CLI specs
- multiple depth files materialize multiple depth slots

### Phase 6: Documentation And Final Verification

Deliverables:

- CLI reference update
- tutorial/advanced customization update if needed
- final targeted and broad test runs

Exit criteria:

- users can discover the CLI flags
- web session save/load includes linear slots
- tests pass or known limitations are documented

## Risk Areas

### Record Spacing

Linear record spacing currently assumes the scalar plot stack mostly lives below
the axis. Above-axis custom tracks require record top guards and inter-record
spacing to account for both sides. This is the highest-risk part of the change.

Mitigation:

- make resolver return explicit `top_extent` and `bottom_extent`
- add tests for two-record diagrams with above-axis GC/skew and BLAST ribbons

### Comparison Ribbon Anchors

BLAST/protein ribbons connect spaces between records. If custom tracks alter
record extents, ribbon start/end y positions can overlap tracks.

Mitigation:

- keep ribbons outside the slot model in v1
- compute ribbon anchors from resolved per-record bottom/top extents
- add tests with `--blast` and above/below custom tracks

### Feature Labels

External labels may protrude above or below feature tracks. Moving features
above/below the axis can expose label overlap bugs.

Mitigation:

- continue using `_precalculate_label_dimensions()`
- update label height calculations to use the feature slot's resolved side
- test with `show_labels=all` and custom feature placement

### Draw Order

If `SeqRecordGroup` keeps axis/features/labels together, custom `z` behavior is
limited for feature-vs-axis ordering.

Mitigation:

- define v1 `z` as affecting scalar/spacer/feature group ordering at the parent
  level
- only extract a separate axis group if a real rendering requirement demands it

### Web Duplication

Circular track-slot UI code is large. Copying it would create maintenance
burden.

Mitigation:

- create a lean Linear module
- extract small shared helpers only when needed
- keep Circular semantics out of Linear code

## Open Questions

- Should `features` default to `overlay` or `below` when `track_layout=middle`
  and no custom slots are provided? The compatibility answer is: preserve
  current `track_layout`; default custom slots should inherit it.
- Should `linear_track_order` include `axis` as a token? Preferred answer: no.
  The axis is represented by `linear_track_axis_index`, keeping the stored slot
  list homogeneous.
- Should Linear slots support per-slot GC mode, such as percent vs deviation?
  Preferred answer: not in v1. Use the existing global GC content mode.
- Should `ticks` be a slot later? Possibly, but only after axis/ruler drawing is
  extracted cleanly from `SeqRecordGroup`.

## Acceptance Criteria

- `gbdraw linear` accepts `--linear_track_slot` and
  `--linear_track_axis_index`.
- `assemble_linear_diagram_from_records()` accepts `linear_track_slots`.
- Custom slots can reorder and place features, depth, GC content, GC skew, and
  spacer rows above/below the axis.
- Existing linear output remains unchanged when custom slots are not supplied.
- Web GUI exposes Linear Custom Track Slots and serializes them to CLI options.
- Web saved sessions preserve Linear Custom Track Slots with a schema version.
- Tests cover parser, CLI forwarding, API rendering, depth slot indexing, and
  web argument generation.
