# Circular Tick / GC/Skew Overlap Fix Plan

## Background

In circular `tuckin` layouts, tick marks and tick labels can be hidden under the GC content track. This is reproducible with EDL933 when separate strands and GC/skew tracks are enabled.

The issue has two causes:

1. The radial layout resolver computes valid tick geometry, including the correct anchor radius and resolved side/length/preset options, but the normal preset path only passes the anchor radius into `TickGroup`. The renderer then falls back to legacy `label_side`, `tick_side`, `tick_length_px`, and `track_preset` behavior, so the measured and rendered tick footprints can diverge.
2. Preset-generated numeric tracks such as GC content are treated as hard-pinned at their legacy `track_dict` radius, so overlap with the tick reserved band is only warned about instead of being resolved.

The fix should make the radial layout resolver the single source of truth for circular track geometry and reserved radial space. Preset numeric tracks should keep their legacy radii when those radii fit, but their radius should be represented as a preferred anchor that may move and internally compress when the resolver detects a collision with ticks, tick labels, the center definition band, depth axis adornments, or any other already-reserved band. This preferred-anchor behavior should be expressed in the circular track slot API instead of encoded through origin metadata in public-looking `params`.

## Development Principles

### YAGNI

- Do not add new UI controls, top-level config keys, or CLI flags for this fix. The only intentional public surface change is the `CircularTrackSlot.anchor_policy` slot API described below.
- Do not introduce a new layout mode.
- Fix only the incorrect interaction between resolved tick geometry and preset numeric track placement.
- Existing reserved bands such as the center definition band and depth axis labels may be considered by the same resolver logic. That is part of making the resolver authoritative for overlap prevention, not a separate user-facing feature.

### KISS

- Keep the default behavior unchanged when tracks already fit.
- Use existing resolver concepts: `CircularResolvedSlot`, `RadialBand`, `draw_band_px`, `reserved_band_px`, `packing_band_px`, `_place_inside_auto()`, `_place_inside_auto_group()`, `_candidate_widths()`, and `_min_readable_numeric_width_px()`.
- Avoid a separate special-case algorithm for GC content.
- Limit preferred-anchor fallback to inside numeric/depth slots whose slot API explicitly marks the radius as preferred, so the general slot solver stays simple.
- Do not encode preset origin in public slot params or resolved slot output. The API should model the actual layout semantics: hard radius versus preferred radius.

### SOLID

- Keep `radial_layout.py` responsible for resolving collisions and final radial placement.
- Keep `assemble.py` responsible for passing resolved geometry to render groups.
- Keep render groups such as `TickGroup` responsible only for drawing with the geometry they receive.
- Keep `draw_band_px` as the data/body footprint of a track and `reserved_band_px` as the full footprint that must not collide with other reserved geometry.

### DRY

- Reuse existing measurement and placement helpers.
- Avoid duplicating tick or numeric annulus calculations in `assemble.py`.
- Add small predicates/helpers only when they clarify anchor semantics or preferred-anchor placement.

## Proposed Implementation

### 1. Always Pass Resolved Tick Options to the Renderer

Update `_draw_resolved_circular_slot()` in `gbdraw/diagrams/circular/assemble.py`.

Current normal preset behavior calls `TickGroup` without the resolved slot options because `use_slot_tick_options` is tied to `user_slot_mode`. The anchor radius is still passed, but `TickGroup` falls back to legacy side/length/preset behavior, so the rendered tick footprint can differ from the resolver's measured footprint.

Change the call site in `add_record_on_circular_canvas()` so resolved tick options are used for both user slot mode and normal preset mode:

```python
use_slot_tick_options=True
```

Also stop relying on temporarily mutating `canvas_config.circular_track_preset` for tick rendering. The resolver should normalize tick-specific inputs once and store the exact drawing values in `CircularTickLayout`, including the resolved preset alias:

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

`_tick_layout_from_params()` should resolve both supported preset keys:

```python
track_preset = normalize_circular_track_preset(
    str(params.get("preset", params.get("track_preset", "tuckin")))
)
```

`_draw_resolved_circular_slot()` should read tick drawing values from the typed payload, not from `resolved_slot.params`:

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

Expected effect:

- The tick renderer receives `label_side`, `tick_side`, `tick_length_px`, and `track_preset` from `CircularTickLayout`.
- `preset` and `track_preset` aliases are resolved once by the resolver, so measured and rendered tick geometry cannot diverge because of alias handling.
- The rendered tick footprint matches the resolver's measured footprint.
- Tick rendering no longer depends on a mutable `canvas_config.circular_track_preset` side effect.
- No new tick-specific public option is required.

### 2. Represent Preferred Anchors in the Slot API

Update `gbdraw/tracks/circular.py`, `gbdraw/diagrams/circular/presets.py`, and the private radial intent model in `gbdraw/diagrams/circular/radial_layout.py`.

Preset-generated numeric tracks currently have explicit radius values derived from `track_dict`. The first resolver pass treats these as hard-pinned anchors, so a collision with an earlier tick reserved band only emits a warning.

The current `explicit_anchor` flag mixes two meanings:

- The slot input had a radius.
- The radius must not move.

Make the second meaning explicit in the slot API by adding an anchor policy:

```python
CircularTrackAnchorPolicy = Literal["hard", "preferred"]


@dataclass(frozen=True)
class CircularTrackSlot:
    ...
    anchor_policy: CircularTrackAnchorPolicy | None = None


@dataclass(frozen=True)
class NormalizedCircularTrackSlot:
    ...
    anchor_policy: CircularTrackAnchorPolicy | None
```

Normalization rules:

- `anchor_policy=None` preserves existing behavior.
- A slot with `radius is not None` and no policy normalizes to `anchor_policy="hard"`.
- A slot without a radius keeps `anchor_policy=None`.
- `anchor_policy="hard"` is the explicit spelling of the default pinned-radius behavior.
- `anchor_policy="preferred"` requires a radius and is valid only for inside numeric/depth slots for this fix.
- `anchor_policy="preferred"` means the radius is the first placement candidate, not a constraint. If that candidate collides, the resolver may move and internally compress the slot/group.
- `strict=True` does not make a preferred anchor hard. To require an immovable radius, use `anchor_policy="hard"` or omit `anchor_policy`.
- Public `CircularTrackSlot(radius=..., compress=True)` remains invalid. Preferred anchors can use internal fallback compression without setting public `compress=True`.
- Update `parse_circular_track_slot()` to accept `anchor_policy=hard|preferred` as a slot-level field in custom slot strings; do not add any new top-level CLI flag.
- Update the CLI/help documentation for the existing `--circular_track_slot` surface:
  - `gbdraw/circular.py` help text should mention `anchor_policy=hard|preferred`.
  - `docs/CLI_Reference.md` should document that `radius` is hard by default and that `anchor_policy=preferred` is valid only for inside numeric/depth slots.
  - `docs/TUTORIALS/3_Advanced_Customization.md` may add a short Custom Track Slots paragraph or example for preferred numeric/depth anchors.

Preset expansion should remove `_preset_generated` from all slot params. Instead, preset numeric/depth slots should be marked with `anchor_policy="preferred"`:

```python
CircularTrackSlot(
    id=slot_id,
    renderer=renderer,
    side="inside",
    radius=radius,
    width=_scalar_px(_default_numeric_width_px(renderer, context)),
    spacing=_scalar_px(max(1.0, 0.01 * float(context.canvas_config.radius))),
    anchor_policy="preferred",
    params=dict(params or {}),
)
```

Ticks and feature slots should keep the default hard-anchor behavior when they have a radius. This keeps the scope narrow: only numeric/depth tracks with an explicit preferred policy become movable.

Expected effect:

- User-pinned slots remain pinned by default and retain current validation behavior.
- Preset numeric/depth slots no longer overlap ticks or labels just because their legacy radius was explicit.
- The public model documents the semantic difference between "this slot had a radius" and "this radius is immovable".
- `strict` remains a conflict policy for hard anchors only; preferred anchors fail only when no readable non-overlapping layout can be found after preferred-radius and bounded deterministic fallback attempts.
- No resolver-private origin wrapper is needed.

### 3. Use Anchor Policy During Radial Placement

Update `resolve_circular_radial_layout()` and helper validation in `gbdraw/diagrams/circular/radial_layout.py`.

Add `anchor_policy` to `_RadialSlotIntent` from the normalized slot. The resolver should inspect anchor policy, not preset origin:

```python
def _is_preferred_numeric_anchor(intent: _RadialSlotIntent) -> bool:
    return (
        intent.explicit_anchor
        and intent.anchor_policy == "preferred"
        and intent.side == "inside"
        and intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
    )


def _is_hard_anchor(intent: _RadialSlotIntent) -> bool:
    return intent.explicit_anchor and not _is_preferred_numeric_anchor(intent)
```

In the fixed-anchor prepass, place only hard anchors and reserving overlays:

```python
if not _is_hard_anchor(intent) and not (intent.side == "overlay" and intent.reserve):
    continue
```

Preferred numeric/depth anchors should be skipped by the prepass and handled by the ordered placement pass. First try their preferred radius at full width. If it collides, move and optionally compress through the preferred numeric group fallback described below.

The ordered placement pass should enforce same-side order while placing slots, rather than relying on `_validate_same_side_order()` to excuse preferred anchors after the fact. Keep a private placement classification derived from intents:

```python
def _placement_class(intent: _RadialSlotIntent) -> Literal["hard", "preferred", "auto", "overlay"]:
    if intent.side == "overlay":
        return "overlay"
    if _is_preferred_numeric_anchor(intent):
        return "preferred"
    if intent.explicit_anchor:
        return "hard"
    return "auto"
```

Use that classification for placement only; do not add it to `CircularResolvedSlot` output.

Preferred numeric/depth slots should be grouped only when they are contiguous in slot order. Any hard, auto, spacer, non-numeric, outside, or overlay slot breaks the preferred group. A non-contiguous preferred slot becomes a separate group and is placed with the same rules.

For preferred and auto placement, derive a same-side placement window before placing:

```python
@dataclass(frozen=True)
class PlacementWindow:
    inner_px: float
    outer_px: float
```

Inside slots are ordered outside-to-inside:

- A previous resolved inside slot defines the window outer limit: the current slot/group's packing outer radius must be no larger than `previous.packing_band_px.inner_px - previous.spacing_px`.
- A future hard inside slot defines the window inner limit: the current slot/group's innermost packing band must remain outside `future_hard.packing_band_px.outer_px + current_spacing_px`.
- For preferred groups, use the first group member when applying the previous resolved boundary and the last group member when applying the future hard boundary.

Outside slots are ordered inside-to-outside:

- A previous resolved outside slot defines the window inner limit: the current slot/group's packing inner radius must be no smaller than `previous.packing_band_px.outer_px + previous.spacing_px`.
- A future hard outside slot defines the window outer limit: the current slot/group's outermost packing band must remain inside `future_hard.packing_band_px.inner_px - current_spacing_px`.

Extend the inside/outside auto placers so they accept `PlacementWindow(inner_px, outer_px)`. The existing single-boundary `inside_max_outer` / `outside_min_inner` values are not enough for mixed hard/preferred cases. Preferred-radius placement and fallback placement must reject candidates that violate the window, internal group order, or internal group spacing before they are added to `resolved_by_index` or `occupied`.

Keep `occupied` as the source of collision truth for already-reserved radial bands. When a candidate is accepted, update `occupied` exactly once for each reserving slot, then update the running same-side boundary from that candidate's packing band. When a candidate is rejected, keep a compact diagnostic reason such as the first conflicting owner or the violated window. Failure messages should name the slot/group and the radial window, for example:

```text
Circular track slot 'gc_content' cannot fit inside between 184.2px and 260.0px without overlapping ticks/definition.
Preferred numeric group 'depth,gc_content,gc_skew' cannot fit inside between 90.0px and 260.0px.
```

After all placements are resolved, keep `_validate_same_side_order()` as a safety net. It should receive the private placement-class map if needed for diagnostics, but it should not inspect `params["_preset_generated"]`, should not add origin fields to resolved output, and should not use a broad preferred-anchor exception. Any returned auto/preferred placement that violates order is a resolver failure and should raise `ValidationError`. Existing hard-anchor compatibility remains limited to the current pinned-geometry behavior: hard anchors are never silently moved, and `strict=False` may warn for accepted pinned overlaps where the current public API already allows that.

Remove `_preset_generated` from preset slot params entirely, and do not inspect `params["_preset_generated"]` anywhere in the resolver. A user-provided `params={"_preset_generated": True}` must have no layout effect.

### 4. Place Preferred Numeric Tracks as Compressible Groups

Inside numeric/depth tracks with `anchor_policy="preferred"` should use their radius as a preferred anchor, not a hard pin. Preset-generated numeric/depth tracks get this policy by default, and user-provided slots may opt into the same behavior explicitly.

Handle contiguous preferred numeric/depth slots as one stack, for example:

```text
depth -> gc_content -> gc_skew
```

The group-level flow should be:

1. Measure the whole group at its preferred legacy radii and full widths.
2. Accept those placements if every reserved band fits the current occupied bands, every slot fits within the group's same-side inner/outer placement limits, and the numeric bands satisfy group ordering and slot spacing constraints.
3. If the preferred group collides or violates its ordering boundaries, fall back to `_place_inside_auto_group()` with the same inner/outer limits.
4. During fallback, allow internal compression without setting public `CircularTrackSlot.compress=True`.
5. Never compress below `_min_readable_numeric_width_px()`.
6. If the readable minimum still cannot fit, raise `ValidationError` instead of overlapping the tick label, definition, or other reserved bands.

For depth tracks, keep `draw_band_px` as the depth plot body only. Do not inflate `draw_band_px` to include depth axis ticks or labels, and do not shrink depth axis ticks or labels to force them into `draw_band_px`. Instead, expand only `reserved_band_px` for a conservative radial footprint of depth axis adornments.

The resolver should not need `depth_df`. Tick values and exact labels are data-dependent drawing details, while the radial layout only needs enough reserved radial space to avoid hiding other tracks. Use `depth_config` plus the candidate track width to compute a conservative footprint:

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


def depth_reserved_band_for_draw_band(
    draw_band_px: RadialBand,
    footprint: DepthAxisFootprint | None,
) -> RadialBand:
    ...
```

Do not build this as an ad hoc closure in `assemble.py`. Instead, add a small pure shared helper, for example `gbdraw/layout/circular_depth_axis.py`, and use it from the radial resolver and the circular depth drawer where their responsibilities overlap. The shared helper should own axis stroke/tick constants and dynamic tick-font sizing; tick value generation and exact label formatting can remain in the drawer because the resolver does not need them.

Resolver measurement for depth should look like:

```python
draw_band_px = _band_from_center_width(center_radius_px, depth_track_width_px)
footprint = (
    resolve_depth_axis_footprint(depth_config, depth_track_width_px)
    if intent.renderer == "depth" and depth_config is not None
    else None
)
reserved_band_px = depth_reserved_band_for_draw_band(draw_band_px, footprint)
```

The conservative footprint rules are:

- For GC/skew, `reserved_band_px` remains the same as `draw_band_px`.
- For depth with axis disabled, `reserved_band_px` remains the same as `draw_band_px`.
- For depth with axis enabled, reserve at least half the axis stroke width radially.
- For depth with axis and ticks enabled, also reserve a conservative radial tick-label footprint based on the resolved tick font size. This may reserve space even when a particular data set would emit no labels, but it avoids a `depth_df` dependency in the resolver.
- Tick line length, label text width, and label gap are primarily angular/tangential and should not be added to the radial band in this fix.
- Recalculate the reserved band for every candidate width because the default depth tick font size depends on track width.
- If tangential label collisions later matter, handle them as a separate angle-aware collision model.

Implementation notes:

- Do not set preset numeric slots to `compress=True` in the public slot object just to unlock this behavior. Use `anchor_policy="preferred"`.
- Keep public `CircularTrackSlot(radius=..., compress=True)` invalid. Preferred-anchor fallback compression is implicit placement behavior, not the public auto-slot `compress=True` option.
- Add an internal helper or placement flag so only preferred numeric fallback can use compression candidates independent of public `intent.compress`.
- Add a small helper such as `_try_place_preferred_numeric_group_at_anchors()` that:
  - measures each intent at `intent.anchor_offset_px` with full width,
  - rejects the group if any preferred slot lacks an anchor,
  - rejects the group if any reserved band overlaps the current occupied bands,
  - rejects the group if any packing/reserved band exceeds the current inside placement boundaries,
  - rejects the group if adjacent group members violate inside ordering or configured spacing,
  - uses the same `working_occupied` and placement-limit update rules as auto group placement before returning the measured group unchanged.
- Keep `_place_inside_auto_group()` as the fallback placer, but extend it so the caller can provide candidate width lists, a `PlacementWindow`, and a compression policy independent of the public `intent.compress` flag. It does not need to know about preferred anchors; the preferred-anchor helper decides whether fallback is required.
- For preferred-anchor fallback, explore width candidate combinations with a deterministic priority queue ordered by:
  1. lower total shrinkage,
  2. fewer compressed tracks,
  3. lower maximum shrink ratio,
  4. the candidate index tuple, as a final deterministic tie-breaker.
- Add explicit search bounds, for example:

```python
MAX_PREFERRED_GROUP_SIZE = 6
MAX_WIDTH_COMBINATIONS = 512
```

- For the normal preset case of at most three numeric/depth tracks, keep the existing nine candidate widths per track.
- For larger user-provided preferred groups, reduce each track's candidate list before search, for example to full, midpoint, and readable minimum widths, then stop after `MAX_WIDTH_COMBINATIONS` tested combinations.
- If the bounded search cannot find a layout, raise `ValidationError` with a message that names the group and the placement window. Do not fall back to warning-only overlap.
- Add a fixed-width group placer, for example `_try_place_inside_auto_group_with_widths()`, so each candidate combination can be tested without forcing every track to use the same candidate index.
- Keep full-width, legacy-radius output unchanged whenever the preferred group fits.
- Place the group together so one numeric track does not consume space needed by the next numeric track.

Expected effect:

- Existing non-overlapping diagrams keep their legacy numeric track radii.
- EDL933-like diagrams move GC/skew/depth inward only when needed.
- Tight layouts can reduce only the numeric/depth widths that actually need reduction, down to readable minimums, rather than failing too early or shrinking every track equally by default.
- If no readable layout exists, the resolver fails explicitly instead of producing hidden ticks or labels.
- Depth plot width remains stable as a data-track width, while the radial depth axis/tick/label footprint participates in overlap prevention through `reserved_band_px`.

### 5. Preserve User Slot Semantics

Do not relax explicit anchors provided by users unless they explicitly opt into `anchor_policy="preferred"`.

User-provided `CircularTrackSlot(radius=...)` with no `anchor_policy` should continue to behave as pinned geometry. If it conflicts and `strict=True`, it should still raise. If it conflicts and `strict=False`, it may warn as today, but it must not silently move.

User-provided `CircularTrackSlot(radius=..., anchor_policy="preferred")` should opt into movable preferred-anchor behavior. In that mode, `strict=True` does not require the preferred radius to be preserved. The resolver should:

1. try the preferred radius at full width,
2. move and/or internally compress the slot/group if the preferred placement collides,
3. raise `ValidationError` only when no readable non-overlapping placement is found by the bounded deterministic fallback search.

If a user wants `strict=True` to mean "do not move this radius", they should use `anchor_policy="hard"` or omit `anchor_policy`.

User-provided `CircularTrackSlot(radius=..., compress=True)` should remain invalid through public slot validation. Internal compression for preferred numeric/depth tracks must not relax that public API rule.

This keeps the public slot API predictable.

### 6. Keep Drawing Order Stable Unless Tests Prove Otherwise

Do not fix the issue by drawing ticks last.

Drawing ticks last would hide the symptom while the layout still contains overlapping reserved bands. The resolver should produce non-overlapping radial bands first. Drawing order can remain based on `(z, slot_index)` unless a separate layering bug is found.

## Test Plan

Add or update focused tests in `tests/test_circular_track_slots.py` and/or `tests/test_circular_radial_layout.py`.

Required coverage:

1. Normal preset tick rendering receives resolved tick options.
   - Verify `label_side`, `tick_side`, `tick_length_px`, and `track_preset` are passed to `add_tick_group_on_canvas()` outside user slot mode.
   - Verify `canvas_config.circular_track_preset` mutation is no longer used for the tick render path.
   - Verify the renderer receives `track_preset` through the explicit `add_tick_group_on_canvas()` argument.
   - Verify both `params={"preset": "..."}` and `params={"track_preset": "..."}` resolve to the same `CircularTickLayout.track_preset` value and use the same value for measuring and drawing.

2. EDL933 `tuckin` layout avoids tick/numeric overlap.
   - Use EDL933 with `track_type="tuckin"`, `strandedness=True`, `show_gc=True`, and `show_skew=True`.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick path band.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick label band when labels exist.
   - Include configured slot spacing in overlap assertions.

3. Depth reserved geometry includes radial axis/tick/label footprint without changing plot body width.
   - Add resolver-level coverage for depth with `show_axis=True` and `show_ticks=True`.
   - Assert depth `draw_band_px` remains the data plot width.
   - Assert depth `reserved_band_px` expands beyond `draw_band_px` by the conservative radial depth axis/tick/label footprint derived from `depth_config` and track width.
   - Assert disabling the depth axis returns the reserved band to the data plot footprint.
   - Assert the depth reserved band is calculated from the same shared `gbdraw/layout/circular_depth_axis.py` helper/constants used by the circular depth drawer for axis constants and tick-font sizing.
   - Include a case that recalculates the conservative footprint for a compressed candidate width.
   - Assert the resolver does not require `depth_df` to measure the radial depth footprint.
   - Assert tangential label width/gap is not added to the radial reserved band in this fix.

4. Depth/GC/skew stacks avoid all relevant reserved bands.
   - Include a resolver-level test with depth, GC content, and GC skew enabled.
   - Assert depth, radial depth axis/tick/label reserved footprint, GC content, GC skew, ticks, tick labels, feature bands, and the center definition band do not overlap.
   - Assert numeric tracks do not overlap each other and preserve ordering/spacing.

5. Legacy radii are preserved when no overlap exists.
   - Use an existing smaller/non-conflicting circular record.
   - Assert depth/GC/skew centers match the previous `track_dict`-derived radii.
   - Assert widths remain full default widths.

6. User-pinned numeric slots remain pinned.
   - Provide an explicit radius in a user slot with no `anchor_policy`.
   - Verify the resolver does not silently move it.
   - If the slot is strict and conflicts, verify it still raises.
   - Add a regression case where the user slot has `params={"_preset_generated": True}` and verify it is ignored and still treated as user-pinned.
   - Verify user-provided `CircularTrackSlot(radius=..., compress=True)` remains invalid.
   - Verify `anchor_policy="hard"` is equivalent to the default pinned behavior.
   - Verify `anchor_policy="preferred"` opts a numeric/depth slot into movable preferred-anchor behavior.
   - Verify `anchor_policy="preferred", strict=True` still moves away from the preferred radius when required, and raises only when the bounded deterministic fallback cannot find a readable non-overlapping placement.

7. Preferred numeric fallback can compress.
   - Construct a tight inside layout where preferred full-width numeric slots collide but a compressed readable stack fits.
   - Assert at least one numeric width shrinks.
   - Assert tracks that do not need shrinking remain at full width when a smaller-change combination fits.
   - Assert no numeric width drops below `_min_readable_numeric_width_px()`.
   - Assert the selected compressed combination is the minimum-change fitting combination according to the fallback ordering: total shrinkage, compressed track count, maximum shrink ratio, and candidate index tuple.
   - Add a larger user-preferred numeric group test that exercises the reduced candidate list and `MAX_WIDTH_COMBINATIONS` bound without accepting overlap.

8. Preferred numeric fallback fails when no readable layout exists.
   - Construct a layout where even readable minimum widths cannot fit.
   - Assert `ValidationError` is raised.
   - Assert the failure is not converted into a warning-only overlap.

9. Anchor policy replaces preset origin metadata.
   - Assert `circular_track_slots_for_preset()` no longer emits `_preset_generated` in slot params.
   - Assert preset numeric/depth slots emit `anchor_policy="preferred"`.
   - Assert preset ticks/features with radii remain hard/default anchors.
   - Assert otherwise identical numeric/depth slots with no policy or `anchor_policy="hard"` stay pinned.
   - Add a mixed-policy resolver-level test where only the preferred numeric slot can move and the hard pinned slot remains pinned.
   - Assert no `generated_by_preset` field is emitted on `CircularResolvedSlot` or captured in the preset geometry oracle.
   - Update `_validate_same_side_order()` coverage to use private resolver-derived placement classification or diagnostics, not params or resolved-slot fields.
   - Add mixed hard/preferred ordering tests where a future hard inside slot creates an inner boundary for an earlier preferred group.
   - Add a non-contiguous preferred-slot test and verify the preferred groups are placed independently while preserving same-side order.
   - Add a fallback case that would fit geometrically but violate same-side order, and assert the placer rejects it instead of relying on a post-hoc validator exception.

10. Preferred legacy radii are tried before auto fallback.
   - Construct a preset numeric/depth group whose full-width preferred radii fit.
   - Assert the resolved centers match the legacy `track_dict` radii and `compressed` is false.
   - Construct a collision case and assert fallback moves/compresses the group instead of accepting overlap.

11. SVG-level regression coverage stays focused and stable.
   - Generate the EDL933 `tuckin`, `strandedness=True`, GC/skew SVG through the normal assembly path.
   - Assert the final SVG contains the tick path group and tick label text/path elements.
   - Assert normal preset rendering passes resolved tick options to `add_tick_group_on_canvas()` without relying on `canvas_config.circular_track_preset` mutation.
   - Do not add pixel-image comparison for this fix; font and CairoSVG differences would make the test less stable than resolver geometry assertions.

12. Circular preset oracle records anchor policy but omits origin metadata.
   - Update `tests/test_circular_preset_geometry.py` and `tests/fixtures/circular_preset_oracle/*.json`.
   - Remove `_preset_generated` from captured `params`.
   - Do not include `generated_by_preset` in captured resolved slot data.
   - Include public `anchor_policy` only where the oracle captures public slot input or normalized public slot data.
   - If the oracle schema changes, bump the fixture schema and document why.

13. CLI and docs mention the new slot-level policy.
   - Update `gbdraw/circular.py` help text for `--circular_track_slot`.
   - Update `docs/CLI_Reference.md`.
   - Add or update a short Custom Track Slots explanation in `docs/TUTORIALS/3_Advanced_Customization.md` if it helps users understand preferred anchors.

## Acceptance Criteria

- The EDL933 screenshot scenario no longer hides tick marks or tick labels under GC content.
- Normal circular presets and custom slot mode both use the same resolved tick geometry.
- Preset numeric/depth tracks preserve legacy radii when they fit.
- Preset numeric/depth tracks move and compress only when required to avoid overlap.
- Depth `draw_band_px` remains the depth plot body footprint, while depth `reserved_band_px` accounts for conservative radial axis/tick/label adornments through the shared depth-axis footprint helper without requiring `depth_df`.
- User-pinned numeric/depth slots remain pinned.
- User-pinned numeric/depth slots still reject `compress=True` when `radius` is supplied.
- Preferred numeric/depth slots treat `strict=True` as "must find a non-overlapping layout within the bounded deterministic fallback search", not "must keep the preferred radius".
- Circular slot inputs support `anchor_policy="hard" | "preferred"` with backward-compatible defaults.
- Existing `--circular_track_slot` help and user-facing docs describe `anchor_policy=hard|preferred`.
- Preset origin is not tracked through `_preset_generated` in public slot params or through resolved slot output.
- Preferred/auto placement enforces same-side order during placement, including mixed hard/preferred and non-contiguous preferred groups; `_validate_same_side_order()` remains a safety net and does not inspect params or add origin fields to `CircularResolvedSlot`.
- No new UI controls, top-level config keys, CLI flags, or layout modes are added.
- Existing non-overlapping diagrams remain visually stable.
- Resolver-level geometry tests cover tick, numeric, depth, and definition reserved-band non-overlap.
- Focused SVG-level regression tests confirm final tick elements are emitted and normal presets pass resolved tick options.
- Relevant circular track tests pass.

## Files Expected to Change

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/presets.py`
- `gbdraw/tracks/circular.py`
- `gbdraw/circular.py`
- `gbdraw/layout/circular_depth_axis.py` or a similarly neutral shared depth-axis geometry helper
- `gbdraw/render/drawers/circular/depth.py`
- `docs/CLI_Reference.md`
- `docs/TUTORIALS/3_Advanced_Customization.md`
- `tests/test_circular_track_slots.py`
- `tests/test_circular_radial_layout.py`
- `tests/test_circular_preset_geometry.py`
- `tests/fixtures/circular_preset_oracle/*.json`

Reference SVG files should only be updated if intentional output changes are observed in snapshot/reference tests.
