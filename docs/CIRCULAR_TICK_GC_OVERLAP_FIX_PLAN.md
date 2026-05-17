# Circular Tick / GC/Skew Overlap Fix Plan

## Background

In circular `tuckin` layouts, tick marks and tick labels can be hidden under the GC content track. This is reproducible with EDL933 when separate strands and GC/skew tracks are enabled.

The issue has two causes:

1. The radial layout resolver computes a valid tick geometry, but the normal preset path does not pass the resolved tick options into `TickGroup`.
2. Preset-generated numeric tracks such as GC content are treated as pinned at their legacy `track_dict` radius, so overlap with the tick reserved band is only warned about instead of being resolved.

The fix should make the radial layout resolver the single source of truth for circular track geometry and reserved radial space. Preset numeric tracks should keep their legacy radii when those radii fit, but they must become movable and internally compressible when the resolver detects a collision with ticks, tick labels, the center definition band, depth axis adornments, or any other already-reserved band. Preset origin must be resolver-private metadata, not a public-looking `params` key.

## Development Principles

### YAGNI

- Do not add new UI controls, config keys, or CLI flags for this fix.
- Do not introduce a new layout mode.
- Fix only the incorrect interaction between resolved tick geometry and preset numeric track placement.
- Existing reserved bands such as the center definition band and depth axis labels may be considered by the same resolver logic. That is part of making the resolver authoritative for overlap prevention, not a separate user-facing feature.

### KISS

- Keep the default behavior unchanged when tracks already fit.
- Use existing resolver concepts: `CircularResolvedSlot`, `RadialBand`, `draw_band_px`, `reserved_band_px`, `packing_band_px`, `_place_inside_auto()`, `_place_inside_auto_group()`, `_candidate_widths()`, and `_min_readable_numeric_width_px()`.
- Avoid a separate special-case algorithm for GC content.
- Limit preferred-anchor fallback to resolver-known preset numeric/depth groups so the general slot solver stays simple.
- Keep preset-origin handling private to the circular radial layout plumbing. Do not encode origin in public slot params or resolved slot output.

### SOLID

- Keep `radial_layout.py` responsible for resolving collisions and final radial placement.
- Keep `assemble.py` responsible for passing resolved geometry to render groups.
- Keep render groups such as `TickGroup` responsible only for drawing with the geometry they receive.
- Keep `draw_band_px` as the data/body footprint of a track and `reserved_band_px` as the full footprint that must not collide with other reserved geometry.

### DRY

- Reuse existing measurement and placement helpers.
- Avoid duplicating tick or numeric annulus calculations in `assemble.py`.
- Add small predicates/helpers only when they clarify anchor semantics or preset-origin handling.

## Proposed Implementation

### 1. Always Pass Resolved Tick Options to the Renderer

Update `_draw_resolved_circular_slot()` in `gbdraw/diagrams/circular/assemble.py`.

Current normal preset behavior calls `TickGroup` without the resolved slot options because `use_slot_tick_options` is tied to `user_slot_mode`. That causes `TickGroup` to fall back to `legacy` tick geometry and draw the tick band at a different radius from the one measured by the resolver.

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
- No new public option is required.

### 2. Distinguish Hard Anchors From Soft Preset Anchors

Update `resolve_circular_radial_layout()` and the private radial intent model in `gbdraw/diagrams/circular/radial_layout.py`.

Preset-generated numeric tracks currently have explicit radius values derived from `track_dict`. The first resolver pass treats these as pinned anchors, so a collision with an earlier tick reserved band only emits a warning.

The current `explicit_anchor` flag mixes two meanings:

- The slot input had a radius.
- The radius must not move.

Split those semantics internally inside the resolver intent model:

- `explicit_anchor`: the slot input had a radius.
- `generated_by_preset`: the slot came from preset expansion, not from the user slot API.
- `hard_anchor`: a radius that must not move.
- `soft_anchor`: a preset numeric/depth radius that should be tried first, but may move if it collides.

Add `generated_by_preset` only to `_RadialSlotIntent`. Do not add it to `CircularResolvedSlot`, `CircularRadialLayout`, or oracle fixtures. The origin flag is placement metadata, not resolved geometry. Tests should verify behavior rather than assert that origin metadata appears in public-looking resolved output.

Implementation shape:

```python
def _is_soft_preset_numeric(intent: _RadialSlotIntent) -> bool:
    return (
        intent.generated_by_preset
        and intent.explicit_anchor
        and intent.side == "inside"
        and intent.renderer in NUMERIC_CIRCULAR_TRACK_RENDERERS
    )


def _is_hard_anchor(intent: _RadialSlotIntent) -> bool:
    return intent.explicit_anchor and not _is_soft_preset_numeric(intent)
```

In the fixed-anchor prepass, place only hard anchors and reserving overlays:

```python
if not _is_hard_anchor(intent) and not (intent.side == "overlay" and intent.reserve):
    continue
```

Soft preset numeric slots should be skipped by the prepass and handled by the ordered placement pass.

For same-side order validation, keep origin metadata off `CircularResolvedSlot`. Instead, compute a private set from resolver intents and pass it into `_validate_same_side_order()`:

```python
preset_generated_slot_indices = {
    intent.slot_index
    for intent in intents
    if intent.generated_by_preset
}

_validate_same_side_order(
    resolved_slots,
    spacing_by_index,
    preset_generated_slot_indices=preset_generated_slot_indices,
)
```

Then `_validate_same_side_order()` can continue to operate on resolved slots while checking preset-generated status by `slot.slot_index`. This avoids duplicating origin fields in resolved output and avoids relying on `params["_preset_generated"]`.

Expected effect:

- User-pinned slots remain pinned and retain current validation behavior.
- Preset numeric/depth slots are no longer allowed to overlap ticks or labels just because their legacy radius was explicit.
- The code documents the difference between "this slot had a radius" and "this slot is immovable".

### 3. Track Preset Origin as Resolver-Internal Intent Metadata

Do not rely on `_preset_generated` inside `CircularTrackSlot.params` for behavior that affects user slot semantics. `params` is public-facing enough that a user can accidentally or deliberately pass the same key.

Remove `_preset_generated` from preset slot params entirely, and do not inspect `params["_preset_generated"]` anywhere in the resolver.

Do not pass caller-computed `preset_slot_indices`. Having the caller predict normalized `slot_index` values is brittle because `slot_index` is assigned during normalization and disabled slots are skipped after enumeration.

Use a private resolver input wrapper so slot identity and slot origin travel together:

```python
@dataclass(frozen=True)
class _RadialSlotInput:
    slot: CircularTrackSlot
    generated_by_preset: bool = False
```

This wrapper is an implementation detail for circular assembly and radial layout. Do not expose it as a CLI/config option, do not add it to public API exports, and do not serialize it in oracle fixtures.

Plain `CircularTrackSlot` inputs should be treated as user-origin slots:

```python
def _coerce_radial_slot_inputs(
    slots: Sequence[CircularTrackSlot | _RadialSlotInput],
) -> list[_RadialSlotInput]:
    out: list[_RadialSlotInput] = []
    for item in slots:
        if isinstance(item, _RadialSlotInput):
            out.append(item)
        else:
            out.append(_RadialSlotInput(slot=item, generated_by_preset=False))
    return out
```

When normal preset mode expands slots, wrap those slots as preset-origin before calling the resolver:

```python
layout_slots = circular_track_slots_for_preset(...)
resolver_slots = [
    _RadialSlotInput(slot=slot, generated_by_preset=True)
    for slot in layout_slots
]
```

When user slot mode is active:

```python
resolver_slots = [
    _RadialSlotInput(slot=slot, generated_by_preset=False)
    for slot in layout_slots
]
```

Then call:

```python
radial_layout = resolve_circular_radial_layout(
    ...,
    slots=resolver_slots,
)
```

Inside `resolve_circular_radial_layout()`, pass the coerced private inputs into `_slot_intents()`. Normalize the underlying `CircularTrackSlot` list once, then copy origin from the matching input by the normalized slot's original `slot_index`:

```python
slot_inputs = _coerce_radial_slot_inputs(slots)
origin_by_original_index = {
    index: slot_input.generated_by_preset
    for index, slot_input in enumerate(slot_inputs)
}
normalized_slots = normalize_circular_track_slots([item.slot for item in slot_inputs])

generated_by_preset = bool(origin_by_original_index.get(normalized_slot.slot_index, False))
```

Keep that field on `_RadialSlotIntent` only.

Any existing logic that needs to distinguish preset slots from user slots must use `intent.generated_by_preset`, including:

- fixed-anchor conflict handling,
- `_validate_same_side_order()`.

This keeps user-provided `CircularTrackSlot(radius=...)` pinned even if its `params` include `_preset_generated`.

This also keeps mixed-origin slot lists possible without a fragile public index set, while preserving the existing plain-slot behavior for tests and callers that do not need origin metadata.

### 4. Place Preset Numeric Tracks as Preferred, Compressible Groups

Preset-generated inside numeric/depth tracks should use their `track_dict` radii as preferred anchors, not hard pins.

Handle contiguous soft preset numeric/depth slots as one stack, for example:

```text
depth -> gc_content -> gc_skew
```

The group-level flow should be:

1. Measure the whole group at its preferred legacy radii and full widths.
2. Accept those placements if every reserved band fits the current occupied bands, every slot fits inside the current `max_packing_outer_px`, and the numeric bands satisfy group ordering and slot spacing constraints.
3. If the preferred group collides, fall back to `_place_inside_auto_group()`.
4. During fallback, allow internal compression without setting public `CircularTrackSlot.compress=True`.
5. Never compress below `_min_readable_numeric_width_px()`.
6. If the readable minimum still cannot fit, raise `ValidationError` instead of overlapping the tick label, definition, or other reserved bands.

For depth tracks, keep `draw_band_px` as the depth plot body only. Do not inflate `draw_band_px` to include depth axis ticks or labels, and do not shrink depth axis ticks or labels to force them into `draw_band_px`. Instead, expand only `reserved_band_px` for the extra radial footprint of depth axis adornments:

```python
draw_band_px = _band_from_center_width(center_radius_px, depth_track_width_px)
reserved_band_px = _reserved_band_for_numeric_track(
    renderer=intent.renderer,
    draw_band_px=draw_band_px,
    width_px=depth_track_width_px,
    depth_config=depth_config,
)
```

The helper should be conservative and DRY:

- For GC/skew, `reserved_band_px` remains the same as `draw_band_px`.
- For depth with axis/tick labels disabled, `reserved_band_px` remains the same as `draw_band_px`.
- For depth with axis/tick labels enabled, `reserved_band_px` expands enough to cover the axis stroke, major/minor tick lines, label gap, and tick label font footprint used by the depth drawer.
- Extract shared depth axis footprint constants or a pure helper from the depth drawer instead of duplicating magic numbers in the resolver.
- Pass `depth_config` or a small depth-axis metrics object into `resolve_circular_radial_layout()` so the resolver measures the same axis/tick/label visibility settings that the renderer will draw.
- Treat the depth axis footprint as a radial reserved band first. If this later proves too conservative, an angle-aware collision model can be considered separately.

Implementation notes:

- Do not set preset numeric slots to `compress=True` in the public slot object just to unlock this behavior.
- Keep public pinned slot `compress=True` invalid. A user-supplied `radius` means the track is pinned; allowing compression there would make the API ambiguous.
- Add an internal helper or placement flag so only soft preset numeric fallback can use compression candidates without changing user-visible slot semantics.
- Add a small helper such as `_try_place_soft_preset_numeric_group_at_preferred_anchors()` that:
  - measures each intent at `intent.anchor_offset_px` with full width,
  - rejects the group if any preferred slot lacks an anchor,
  - rejects the group if any reserved band overlaps the current occupied bands,
  - rejects the group if any packing/reserved band exceeds the current inside placement boundary,
  - rejects the group if adjacent group members violate inside ordering or configured spacing,
  - uses the same `working_occupied` and `working_outer` update rules as auto group placement before returning the measured group unchanged.
- Keep `_place_inside_auto_group()` as the fallback placer, but extend it so the caller can provide candidate width lists or a compression policy independent of the public `intent.compress` flag. It does not need to know about preferred anchors; the preferred-anchor helper decides whether fallback is required.
- For soft preset fallback, explore width candidate combinations in an order that minimizes visual change:
  1. lower total shrinkage,
  2. fewer compressed tracks,
  3. lower maximum shrink ratio.
- With three preset numeric/depth tracks and the existing nine candidate widths, a bounded exhaustive search is small enough to keep the code straightforward.
- Add a fixed-width group placer, for example `_try_place_inside_auto_group_with_widths()`, so each candidate combination can be tested without forcing every track to use the same candidate index.
- Keep full-width, legacy-radius output unchanged whenever the preferred group fits.
- Place the group together so one numeric track does not consume space needed by the next numeric track.

Expected effect:

- Existing non-overlapping diagrams keep their legacy numeric track radii.
- EDL933-like diagrams move GC/skew/depth inward only when needed.
- Tight layouts can reduce only the numeric/depth widths that actually need reduction, down to readable minimums, rather than failing too early or shrinking every track equally by default.
- If no readable layout exists, the resolver fails explicitly instead of producing hidden ticks or labels.
- Depth plot width remains stable as a data-track width, while depth axis/tick/label footprint participates in overlap prevention through `reserved_band_px`.

### 5. Preserve User Slot Semantics

Do not relax explicit anchors provided by users.

Only resolver-known preset slots should become soft. User-provided `CircularTrackSlot(radius=...)` should continue to behave as pinned geometry. If it conflicts and `strict=True`, it should still raise. If it conflicts and `strict=False`, it may warn as today, but it must not silently move.

User-provided `CircularTrackSlot(radius=..., compress=True)` should remain invalid through the existing public slot validation. Internal compression for soft preset numeric/depth tracks must not relax that public API rule.

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

3. Depth reserved geometry includes axis/tick/label footprint without changing plot body width.
   - Add resolver-level coverage for depth with `show_axis=True` and `show_ticks=True`.
   - Assert depth `draw_band_px` remains the data plot width.
   - Assert depth `reserved_band_px` expands beyond `draw_band_px` by the depth axis/tick/label footprint.
   - Assert disabling depth axis/tick labels returns the reserved band to the data plot footprint.
   - Assert the depth reserved band is calculated from the same shared helper/constants used by the depth drawer.

4. Depth/GC/skew stacks avoid all relevant reserved bands.
   - Include a resolver-level test with depth, GC content, and GC skew enabled.
   - Assert depth, depth axis/tick/label reserved footprint, GC content, GC skew, ticks, tick labels, feature bands, and the center definition band do not overlap.
   - Assert numeric tracks do not overlap each other and preserve ordering/spacing.

5. Legacy radii are preserved when no overlap exists.
   - Use an existing smaller/non-conflicting circular record.
   - Assert depth/GC/skew centers match the previous `track_dict`-derived radii.
   - Assert widths remain full default widths.

6. User-pinned numeric slots remain pinned.
   - Provide an explicit radius in a user slot.
   - Verify the resolver does not silently move it.
   - If the slot is strict and conflicts, verify it still raises.
   - Add a regression case where the user slot has `params={"_preset_generated": True}` and verify it is ignored and still treated as user-pinned.
   - Verify user-provided `CircularTrackSlot(radius=..., compress=True)` remains invalid.

7. Soft preset numeric fallback can compress.
   - Construct a tight inside layout where preferred full-width numeric slots collide but a compressed readable stack fits.
   - Assert at least one numeric width shrinks.
   - Assert tracks that do not need shrinking remain at full width when a smaller-change combination fits.
   - Assert no numeric width drops below `_min_readable_numeric_width_px()`.
   - Assert the selected compressed combination is the minimum-change fitting combination according to the fallback ordering: total shrinkage, compressed track count, maximum shrink ratio.

8. Soft preset numeric fallback fails when no readable layout exists.
   - Construct a layout where even readable minimum widths cannot fit.
   - Assert `ValidationError` is raised.
   - Assert the failure is not converted into a warning-only overlap.

9. Preset origin does not leak through public params or resolved output.
   - Assert `circular_track_slots_for_preset()` no longer emits `_preset_generated` in slot params.
   - Assert resolver inputs wrapped with `_RadialSlotInput(..., generated_by_preset=True)` treat preset numeric/depth anchors as soft.
   - Assert plain `CircularTrackSlot` inputs, and wrapper inputs with `generated_by_preset=False`, leave otherwise identical slots as user-origin hard anchors.
   - Add a mixed-origin resolver-level test where only the wrapped preset-origin numeric slot can move and the user-origin pinned slot remains pinned.
   - Assert no `generated_by_preset` field is emitted on `CircularResolvedSlot` or captured in the preset geometry oracle.
   - Update `_validate_same_side_order()` coverage to use the private `preset_generated_slot_indices` set derived from resolver intents, not params or resolved-slot fields.

10. Preferred legacy radii are tried before auto fallback.
   - Construct a preset numeric/depth group whose full-width preferred radii fit.
   - Assert the resolved centers match the legacy `track_dict` radii and `compressed` is false.
   - Construct a collision case and assert fallback moves/compresses the group instead of accepting overlap.

11. SVG-level regression coverage stays focused and stable.
   - Generate the EDL933 `tuckin`, `strandedness=True`, GC/skew SVG through the normal assembly path.
   - Assert the final SVG contains the tick path group and tick label text/path elements.
   - Assert normal preset rendering passes resolved tick options to `add_tick_group_on_canvas()` without relying on `canvas_config.circular_track_preset` mutation.
   - Do not add pixel-image comparison for this fix; font and CairoSVG differences would make the test less stable than resolver geometry assertions.

12. Circular preset oracle omits private origin metadata.
   - Update `tests/test_circular_preset_geometry.py` and `tests/fixtures/circular_preset_oracle/*.json`.
   - Remove `_preset_generated` from captured `params`.
   - Do not include `generated_by_preset` in captured resolved slot data.
   - If the oracle schema changes, bump the fixture schema and document why.

## Acceptance Criteria

- The EDL933 screenshot scenario no longer hides tick marks or tick labels under GC content.
- Normal circular presets and custom slot mode both use the same resolved tick geometry.
- Preset numeric/depth tracks preserve legacy radii when they fit.
- Preset numeric/depth tracks move and compress only when required to avoid overlap.
- Depth `draw_band_px` remains the depth plot body footprint, while depth `reserved_band_px` accounts for axis/tick/label adornments.
- User-pinned numeric/depth slots remain pinned.
- User-pinned numeric/depth slots still reject `compress=True` when `radius` is supplied.
- Preset origin is tracked only through private resolver input/intent metadata and never through `_preset_generated` in public slot params or through resolved slot output.
- `_validate_same_side_order()` uses private resolver-derived origin information instead of inspecting params or adding origin fields to `CircularResolvedSlot`.
- No new public configuration is added.
- Existing non-overlapping diagrams remain visually stable.
- Resolver-level geometry tests cover tick, numeric, depth, and definition reserved-band non-overlap.
- Focused SVG-level regression tests confirm final tick elements are emitted and normal presets pass resolved tick options.
- Relevant circular track tests pass.

## Files Expected to Change

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/presets.py`
- `gbdraw/render/drawers/circular/depth.py` and/or a nearby shared depth-axis geometry helper if depth axis footprint constants are extracted
- `tests/test_circular_track_slots.py`
- `tests/test_circular_radial_layout.py`
- `tests/test_circular_preset_geometry.py`
- `tests/fixtures/circular_preset_oracle/*.json`

Reference SVG files should only be updated if intentional output changes are observed in snapshot/reference tests.
