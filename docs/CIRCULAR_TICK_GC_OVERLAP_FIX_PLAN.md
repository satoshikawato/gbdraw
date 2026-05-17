# Circular Tick / GC/Skew Overlap Fix Plan

## Background

In circular `tuckin` layouts, tick marks and tick labels can be hidden under
the GC content track. This is reproducible with EDL933 when separate strands
and GC/skew tracks are enabled.

The issue has two causes:

1. The radial layout resolver computes valid tick geometry, including the
   correct anchor radius and resolved side/length/preset options, but the
   normal preset drawing path only passes the anchor radius into `TickGroup`.
   The renderer then falls back to legacy `label_side`, `tick_side`,
   `tick_length_px`, and `track_preset` behavior, so measured and rendered tick
   footprints can diverge.
2. Preset-generated numeric tracks such as GC content are currently represented
   as ordinary explicit-radius slots. The resolver treats those radii as hard
   pins, so overlap with the tick reserved band can be accepted with a warning
   instead of being resolved.

The fix should make the radial layout resolver the single source of truth for
circular track geometry and reserved radial space. Preset numeric/depth tracks
should keep their legacy radii when those radii fit, but their preset radius
should be treated internally as a preferred anchor that may move and may be
uniformly compressed when the resolver detects a collision with ticks, tick
labels, the center definition band, depth axis adornments, or another reserved
band.

This preferred-anchor behavior is an internal preset layout intent. Do not add
a public `CircularTrackSlot.anchor_policy`, a CLI flag, a config key, or a UI
control for this fix.

## Development Principles

### Scope

- Do not add new public circular slot fields.
- Do not add new CLI flags, top-level config keys, UI controls, or layout
  modes.
- Remove resolver behavior that depends on `params["_preset_generated"]`.
- Fix the interaction between resolved tick geometry and preset numeric/depth
  placement.
- Include depth axis radial footprint in the resolver, but keep it data-frame
  independent.

### Responsibilities

- `radial_layout.py` owns collision checks, same-side ordering, final radial
  placement, and depth reserved-band measurement.
- `presets.py` owns preset expansion and the internal information that preset
  numeric/depth slots have preferred anchors.
- `assemble.py` owns passing resolved geometry to render groups.
- Render groups such as `TickGroup` and `DepthDrawer` draw with the geometry
  they receive.
- `draw_band_px` remains the data/body footprint of a track.
- `reserved_band_px` remains the full radial footprint that must not collide
  with other reserved geometry.

### Simplicity

- Keep default output unchanged when all tracks already fit.
- Preserve preset numeric/depth legacy radii whenever they fit.
- Use deterministic fallback instead of optimal-width search.
- When preferred numeric/depth anchors collide:
  1. try the preferred full-width radii,
  2. try moving the full-width group with inside auto placement,
  3. try uniformly shrinking the group by deterministic scale steps,
  4. raise `ValidationError` if readable minimum widths still cannot fit.
- Do not use per-track width-combination search or priority queues for this
  fix.

## Proposed Implementation

### 1. Always Pass Resolved Tick Options to the Renderer

Update `_draw_resolved_circular_slot()` in
`gbdraw/diagrams/circular/assemble.py`.

The normal preset path currently calls `TickGroup` without the resolved slot
tick options because `use_slot_tick_options` is tied to `user_slot_mode`. The
anchor radius is passed, but `TickGroup` can fall back to legacy side, length,
and preset logic. That can make drawing diverge from resolver measurement.

Change the drawing path so resolved tick options are used in both custom slot
mode and normal preset mode.

Extend `CircularTickLayout` so the typed payload contains every tick drawing
option the renderer needs:

```python
@dataclass(frozen=True)
class CircularTickLayout:
    anchor_radius_px: float
    tick_band_px: RadialBand
    label_band_px: RadialBand | None
    reserved_band_px: RadialBand
    track_preset: str
    label_side: str
    tick_side: str
    tick_length_px: float | None
```

`_tick_layout_from_params()` should normalize both supported preset aliases:

```python
track_preset = normalize_circular_track_preset(
    str(params.get("preset", params.get("track_preset", "tuckin")))
)
```

`_draw_resolved_circular_slot()` should read tick drawing values from
`CircularTickLayout`, not from `resolved_slot.params`:

```python
if isinstance(resolved_slot.payload, CircularTickLayout):
    tick_group_kwargs["label_side"] = resolved_slot.payload.label_side
    tick_group_kwargs["tick_side"] = resolved_slot.payload.tick_side
    tick_group_kwargs["tick_length_px"] = resolved_slot.payload.tick_length_px
    tick_group_kwargs["track_preset"] = resolved_slot.payload.track_preset
```

Delete the temporary mutation block around `add_tick_group_on_canvas()`:

```python
previous_preset = getattr(canvas_config, "circular_track_preset", None)
setattr(canvas_config, "circular_track_preset", ...)
...
setattr(canvas_config, "circular_track_preset", previous_preset)
```

It is acceptable for `canvas_config.circular_track_preset` to remain the
diagram-level preset used by other legacy drawing paths. Tick drawing should no
longer depend on a temporary mutation of that value.

Expected effect:

- The tick renderer receives `label_side`, `tick_side`, `tick_length_px`, and
  `track_preset` from the resolved tick payload.
- `preset` and `track_preset` aliases are normalized once by the resolver.
- Measured and rendered tick footprints use the same values.
- No new public tick option is required.

### 2. Replace Preset Origin Metadata With Internal Placement Intent

Update `gbdraw/diagrams/circular/presets.py`,
`gbdraw/diagrams/circular/assemble.py`, and
`gbdraw/diagrams/circular/radial_layout.py`.

Remove `_preset_generated` from all preset slot params. Do not inspect
`params["_preset_generated"]` anywhere in the resolver. A user-provided
`params={"_preset_generated": True}` must have no layout effect.

Do not add public `CircularTrackSlot.anchor_policy`. Instead, pass preset
preferred-anchor intent through an internal structure used only by the normal
preset assembly path. For example:

```python
@dataclass(frozen=True)
class CircularPresetRadialPlan:
    slots: tuple[CircularTrackSlot, ...]
    preferred_anchor_slot_ids: frozenset[str]
```

Add a preset helper such as `circular_radial_plan_for_preset()` that returns:

- clean `CircularTrackSlot` objects with no `_preset_generated` params,
- `preferred_anchor_slot_ids` containing preset numeric/depth slot ids such as
  `depth`, `gc_content`, and `gc_skew` when those slots have preset radii.

Keep `circular_track_slots_for_preset()` available for existing callers, but
make it return clean slots. The normal assembly path should use the internal
plan and pass its preferred-anchor ids to the resolver. Custom slot mode should
pass no preferred-anchor ids.

Add an internal resolver argument for this intent:

```python
def resolve_circular_radial_layout(
    ...,
    preferred_anchor_slot_ids: Collection[str] = (),
) -> CircularRadialLayout:
    ...
```

Add an internal placement policy to `_RadialSlotIntent`:

```python
PlacementPolicy = Literal["hard", "preferred", "auto", "overlay"]
```

Derive it while building intents:

- `overlay`: `side == "overlay"`.
- `preferred`: slot id is in the internal preset preferred-anchor id set, the
  slot has an explicit radius, `side == "inside"`, and the renderer is numeric
  or depth.
- `hard`: the slot has an explicit radius and is not preferred.
- `auto`: all other non-overlay slots.

If a slot id is marked preferred but the slot has no radius, treat it as `auto`
and optionally log a debug message. This keeps malformed preset intent from
becoming a public validation mode.

Expected effect:

- Preset origin is not encoded through public-looking `params`.
- User custom slots with explicit radius remain hard-pinned by default.
- Only the normal preset path can mark numeric/depth slots as internally
  preferred.
- No public slot API changes are required.

### 3. Enforce Same-Side Order During Placement

Update `resolve_circular_radial_layout()` and its placement helpers in
`gbdraw/diagrams/circular/radial_layout.py`.

The current resolver places explicit anchors first, then uses
`_validate_same_side_order()` as a post-hoc check with a preset-origin
exception. Replace that with candidate-time order enforcement.

Keep `_validate_same_side_order()` as a final safety net, but remove all
`_preset_generated` logic from it. Any final auto/preferred placement that
violates order is a resolver bug and should raise `ValidationError`.

Add a placement window:

```python
@dataclass(frozen=True)
class PlacementWindow:
    inner_px: float
    outer_px: float
```

The window constrains a candidate's packing band or group packing span. The
`occupied` list still remains the source of collision truth for reserved bands.

#### Fixed-Anchor Prepass

In the fixed-anchor prepass, resolve only:

- `placement_policy == "hard"`,
- reserving overlays.

Skip preferred numeric/depth anchors in this pass. They are movable candidates
and must be handled in slot order.

Hard anchors are never moved. If a hard anchor conflicts:

- `strict=True` raises `ValidationError`,
- `strict=False` may warn and preserve the existing pinned-geometry behavior.

#### Inside Order Window

Inside slots are ordered outside-to-inside.

For a current inside slot or preferred group:

- The previous resolved inside slot defines the outer limit:

```text
current/group packing outer <= previous.packing inner - previous.spacing
```

- The next future hard inside slot defines the inner limit:

```text
current/group packing inner >= future_hard.packing outer + last_member.spacing
```

For a preferred group, use the first group member for the previous boundary and
the last group member for the future-hard boundary.

#### Outside Order Window

Outside slots are ordered inside-to-outside.

For a current outside slot:

- The previous resolved outside slot defines the inner limit:

```text
current packing inner >= previous.packing outer + previous.spacing
```

- The next future hard outside slot defines the outer limit:

```text
current packing outer <= future_hard.packing inner - current.spacing
```

#### Candidate Acceptance

Extend inside/outside auto placers so they accept `PlacementWindow`. A candidate
is accepted only if:

- its reserved band does not overlap any occupied reserved band,
- its packing band or group packing span fits the placement window,
- adjacent slots inside a group preserve order and configured spacing.

When a candidate is accepted:

- add each reserving slot's reserved band to `occupied` exactly once,
- update the running same-side boundary from the accepted packing band,
- store the candidate in `resolved_by_index`.

When a candidate is rejected, keep a compact diagnostic reason such as the first
conflicting owner or the violated window. Failure messages should name the slot
or group and the radial window, for example:

```text
Circular track slot 'gc_content' cannot fit inside between 184.2px and 260.0px without overlapping ticks/definition.
Preferred numeric group 'depth,gc_content,gc_skew' cannot fit inside between 90.0px and 260.0px.
```

Expected effect:

- Same-side order is enforced before candidates become final.
- `_validate_same_side_order()` no longer needs preset-origin exceptions.
- Future hard anchors constrain earlier movable preferred/auto slots.
- Hard-pinned user slots remain pinned.

### 4. Place Preferred Numeric/Depth Tracks as Simple Groups

Inside numeric/depth tracks with internal `placement_policy == "preferred"`
should use their radius as a preferred anchor, not a hard pin.

Group only contiguous preferred inside numeric/depth intents in slot order.
Any hard, auto, spacer, non-numeric, outside, or overlay slot breaks the group.
A later non-contiguous preferred slot becomes a separate preferred group and is
placed with the same rules.

Use this deterministic group flow:

1. Measure the whole group at its preferred legacy radii and full widths.
2. Accept those placements if every reserved band fits `occupied`, the group
   packing span fits the `PlacementWindow`, and adjacent group members satisfy
   inside order and spacing.
3. If preferred anchors fail, place the whole group with full widths using
   inside auto placement inside the same window.
4. If full-width auto placement fails, retry auto placement with deterministic
   uniform width scales.
5. Never shrink any numeric/depth width below
   `_min_readable_numeric_width_px()`.
6. If readable minimum widths still cannot fit, raise `ValidationError`.

Do not implement per-track width-combination search. Do not add
`MAX_WIDTH_COMBINATIONS`, priority queues, or "minimum-change" optimality
criteria for this fix.

A simple scale generator is enough. For example:

```python
def _preferred_group_width_scales(intents: Sequence[_RadialSlotIntent]) -> list[float]:
    group_min_scale = max(
        _min_readable_numeric_width_px(intent.renderer, intent.width_px) / intent.width_px
        for intent in intents
        if intent.width_px > LAYOUT_EPSILON
    )
    return _linear_scales_from_1_to_min(group_min_scale, steps=8)
```

The full-width auto attempt covers scale `1.0`; subsequent fallback attempts
can use the remaining descending scales. Widths should be clamped to each
track's readable minimum.

Add small helpers such as:

- `_preferred_numeric_group_from()`
- `_try_place_preferred_numeric_group_at_anchors()`
- `_place_inside_auto_group_with_width_scale()`
- `_group_fits_window_and_order()`

These helpers should reuse existing measurement concepts rather than
duplicating tick or numeric annulus calculations in `assemble.py`.

Expected effect:

- Existing non-overlapping diagrams keep legacy numeric/depth radii and full
  widths.
- EDL933-like diagrams move GC/skew/depth inward only when needed.
- Tight layouts shrink numeric/depth tracks uniformly and deterministically.
- If no readable layout exists, the resolver fails explicitly instead of
  accepting hidden ticks or labels.

### 5. Include Depth Axis Radial Footprint Without `depth_df`

Include depth axis radial footprint in this fix, but keep resolver measurement
independent of `depth_df`.

Depth tick values and exact formatted labels are data-dependent drawing details.
The radial resolver only needs a conservative radial footprint large enough to
avoid hiding other tracks. Use `DepthConfigurator` plus the candidate track
width.

Add an optional resolver argument for this measurement and pass it from
`assemble.py` when depth rendering is enabled:

```python
def resolve_circular_radial_layout(
    ...,
    depth_config: DepthConfigurator | None = None,
) -> CircularRadialLayout:
    ...
```

Add a small shared helper, for example
`gbdraw/layout/circular_depth_axis.py`, that owns the constants and sizing rules
shared by the resolver and the depth drawer:

```python
@dataclass(frozen=True)
class DepthAxisFootprint:
    radial_inner_extra_px: float
    radial_outer_extra_px: float


def resolve_depth_axis_footprint(
    depth_config: DepthConfigurator,
    track_width_px: float,
) -> DepthAxisFootprint:
    ...


def depth_axis_tick_font_size_px(
    depth_config: DepthConfigurator,
    track_width_px: float,
) -> float:
    ...
```

Avoid making the helper import `RadialBand` from `radial_layout.py`. The
resolver can expand its own `RadialBand` using the returned
`DepthAxisFootprint`.

Resolver measurement for depth should look like:

```python
draw_band_px = _band_from_center_width(center_radius_px, depth_track_width_px)
footprint = (
    resolve_depth_axis_footprint(depth_config, depth_track_width_px)
    if intent.renderer == "depth" and depth_config is not None
    else None
)
reserved_band_px = _depth_reserved_band_for_draw_band(draw_band_px, footprint)
```

Conservative footprint rules:

- For GC/skew, `reserved_band_px` remains the same as `draw_band_px`.
- For depth with axis disabled, `reserved_band_px` remains the same as
  `draw_band_px`.
- For depth with axis enabled, reserve at least half the axis stroke width
  radially.
- For depth with axis and ticks enabled, also reserve a conservative radial
  tick-label footprint based on the resolved tick font size.
- Do not add tick line length, label text width, or label gap to the radial
  band in this fix. Those are primarily tangential.
- Recalculate the footprint for each candidate width because default depth tick
  font size depends on track width.
- If tangential label collisions later matter, handle them as a separate
  angle-aware collision model.

Update `DepthDrawer` to use the shared constants and tick-font sizing helper so
the resolver and drawer stay aligned without sharing `depth_df`.

Expected effect:

- Depth plot `draw_band_px` remains the data body width.
- Depth axis/tick/label radial adornments participate in collision prevention
  through `reserved_band_px`.
- The resolver does not need coverage data to measure depth radial footprint.

### 6. Keep Drawing Order Stable

Do not fix the issue by drawing ticks last.

Drawing ticks last would hide the symptom while the layout still contains
overlapping reserved bands. The resolver should produce non-overlapping radial
bands first. Drawing order can remain based on `(z, slot_index)` unless a
separate layering bug is found.

## Test Plan

Add or update focused tests in `tests/test_circular_track_slots.py`,
`tests/test_circular_radial_layout.py`, and
`tests/test_circular_preset_geometry.py`.

Required coverage:

1. Normal preset tick rendering receives resolved tick options.
   - Verify `label_side`, `tick_side`, `tick_length_px`, and `track_preset`
     are passed to `add_tick_group_on_canvas()` outside custom slot mode.
   - Verify the tick render path no longer temporarily mutates
     `canvas_config.circular_track_preset`.
   - Verify both `params={"preset": "..."}` and
     `params={"track_preset": "..."}` resolve to the same
     `CircularTickLayout.track_preset` value and use the same value for
     measuring and drawing.

2. Preset origin metadata is removed.
   - Assert `circular_track_slots_for_preset()` no longer emits
     `_preset_generated` in slot params.
   - Assert a user-provided `params={"_preset_generated": True}` has no layout
     effect and does not turn a user slot into a preferred slot.
   - Assert no `generated_by_preset` or equivalent origin field is emitted on
     `CircularResolvedSlot`.

3. Internal preferred-anchor intent is applied only in normal preset mode.
   - Assert normal preset assembly passes preferred-anchor intent for preset
     numeric/depth slots.
   - Assert custom slot mode does not mark user numeric/depth slots as
     preferred.
   - Assert preset ticks/features with radii remain hard/default anchors.
   - Assert otherwise identical user numeric/depth slots stay pinned.

4. Same-side order is enforced during placement.
   - Add mixed hard/preferred inside tests where a future hard slot creates the
     inner boundary for an earlier preferred group.
   - Add non-contiguous preferred-slot tests and verify preferred groups are
     placed independently while preserving same-side order.
   - Add a fallback case that would fit geometrically but violate same-side
     order, and assert the placer rejects it before final validation.
   - Update `_validate_same_side_order()` tests so it no longer inspects params
     or relies on preset-origin exceptions.

5. EDL933 `tuckin` layout avoids tick/numeric overlap.
   - Use EDL933 with `track_type="tuckin"`, `strandedness=True`,
     `show_gc=True`, and `show_skew=True`.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick
     path band.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick
     label band when labels exist.
   - Include configured slot spacing in overlap assertions.

6. Legacy radii are preserved when no overlap exists.
   - Use an existing smaller/non-conflicting circular record.
   - Assert depth/GC/skew centers match the previous `track_dict`-derived
     radii.
   - Assert widths remain full default widths.

7. Preferred numeric/depth fallback moves before shrinking.
   - Construct a case where preferred full-width radii collide but full-width
     auto placement fits.
   - Assert the group moves and `compressed` remains false.
   - Assert no overlap with ticks, tick labels, definition, feature bands, or
     other numeric/depth tracks.

8. Preferred numeric/depth fallback can uniformly compress.
   - Construct a tight inside layout where full-width preferred and full-width
     auto placement fail, but a uniformly compressed readable stack fits.
   - Assert all compressed preferred group members use the same deterministic
     scale, clamped by readable minimums.
   - Assert no numeric/depth width drops below
     `_min_readable_numeric_width_px()`.
   - Assert the selected scale is the first fitting scale in the deterministic
     scale list.

9. Preferred numeric/depth fallback fails when no readable layout exists.
   - Construct a layout where even readable minimum widths cannot fit.
   - Assert `ValidationError` is raised.
   - Assert the failure is not converted into a warning-only overlap.

10. User-pinned numeric/depth slots remain pinned.
    - Provide an explicit radius in a custom user slot.
    - Verify the resolver does not silently move it.
    - If the slot is strict and conflicts, verify it still raises.
    - If the slot is non-strict and conflicts, verify it preserves the pinned
      placement and follows existing warning behavior.
    - Verify user-provided `CircularTrackSlot(radius=..., compress=True)`
      remains invalid.

11. Depth reserved geometry includes radial axis/tick/label footprint without
    requiring `depth_df`.
    - Add resolver-level coverage for depth with `show_axis=True` and
      `show_ticks=True`.
    - Assert depth `draw_band_px` remains the data plot width.
    - Assert depth `reserved_band_px` expands beyond `draw_band_px` by the
      conservative radial footprint derived from `depth_config` and track
      width.
    - Assert disabling the depth axis returns the reserved band to the data
      plot footprint.
    - Include a compressed-candidate case that recalculates the footprint for
      the compressed width.
    - Assert tangential label width/gap is not added to the radial reserved
      band in this fix.

12. Depth/GC/skew stacks avoid all relevant reserved bands.
    - Include a resolver-level test with depth, GC content, and GC skew
      enabled.
    - Assert depth reserved footprint, GC content, GC skew, ticks, tick labels,
      feature bands, and the center definition band do not overlap.
    - Assert numeric/depth tracks preserve ordering and spacing.

13. SVG-level regression coverage stays focused and stable.
    - Generate the EDL933 `tuckin`, `strandedness=True`, GC/skew SVG through
      the normal assembly path.
    - Assert the final SVG contains the tick path group and tick label
      text/path elements.
    - Do not add pixel-image comparison for this fix; font and CairoSVG
      differences would make the test less stable than resolver geometry
      assertions.

14. Circular preset oracle omits origin metadata.
    - Update `tests/test_circular_preset_geometry.py` and
      `tests/fixtures/circular_preset_oracle/*.json`.
    - Remove `_preset_generated` from captured `params`.
    - Do not include `generated_by_preset` or internal placement policy in
      captured resolved slot data.
    - If the oracle schema changes, bump the fixture schema and document why.

## Acceptance Criteria

- The EDL933 screenshot scenario no longer hides tick marks or tick labels
  under GC content.
- Normal circular presets and custom slot mode both use resolved tick geometry
  for drawing.
- Preset numeric/depth tracks preserve legacy radii when they fit.
- Preset numeric/depth tracks move and uniformly compress only when required to
  avoid overlap.
- Depth `draw_band_px` remains the depth plot body footprint.
- Depth `reserved_band_px` accounts for conservative radial axis/tick/label
  adornments without requiring `depth_df`.
- User-pinned numeric/depth slots remain pinned.
- User-pinned numeric/depth slots still reject `compress=True` when `radius` is
  supplied.
- No public `CircularTrackSlot.anchor_policy`, CLI flag, top-level config key,
  UI control, or new layout mode is added.
- Preset origin is not tracked through `_preset_generated` in public slot params
  or through resolved slot output.
- Preferred/auto placement enforces same-side order during placement, including
  mixed hard/preferred and non-contiguous preferred groups.
- `_validate_same_side_order()` remains a safety net and does not inspect
  params or origin fields.
- Existing non-overlapping diagrams remain visually stable.
- Resolver-level geometry tests cover tick, numeric/depth, depth axis, and
  definition reserved-band non-overlap.
- Focused SVG-level regression tests confirm final tick elements are emitted and
  normal presets pass resolved tick options.
- Relevant circular track tests pass.

## Files Expected to Change

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/presets.py`
- `gbdraw/layout/circular_depth_axis.py` or a similarly neutral shared
  depth-axis geometry helper
- `gbdraw/render/drawers/circular/depth.py`
- `tests/test_circular_track_slots.py`
- `tests/test_circular_radial_layout.py`
- `tests/test_circular_preset_geometry.py`
- `tests/fixtures/circular_preset_oracle/*.json`

Reference SVG files should only be updated if intentional output changes are
observed in snapshot/reference tests.
