# Circular Tick / GC/Skew Overlap Fix Plan

## Background

In circular `tuckin` layouts, tick marks and tick labels can be hidden under the GC content track. This is reproducible with EDL933 when separate strands and GC/skew tracks are enabled.

The issue has two causes:

1. The radial layout resolver computes a valid tick geometry, but the normal preset path does not pass the resolved tick options into `TickGroup`.
2. Preset-generated numeric tracks such as GC content are treated as pinned at their legacy `track_dict` radius, so overlap with the tick reserved band is only warned about instead of being resolved.

The fix should make the radial layout resolver the single source of truth for circular track geometry. Preset numeric tracks should keep their legacy radii when those radii fit, but they must become movable and compressible when the resolver detects a collision.

## Development Principles

### YAGNI

- Do not add new UI controls, config keys, or CLI flags for this fix.
- Do not introduce a new layout mode.
- Fix only the incorrect interaction between resolved tick geometry and preset numeric track placement.

### KISS

- Keep the default behavior unchanged when tracks already fit.
- Use existing resolver concepts: `CircularResolvedSlot`, `RadialBand`, `reserved_band_px`, `packing_band_px`, `_place_inside_auto()`, `_place_inside_auto_group()`, `_candidate_widths()`, and `_min_readable_numeric_width_px()`.
- Avoid a separate special-case algorithm for GC content.

### SOLID

- Keep `radial_layout.py` responsible for resolving collisions and final radial placement.
- Keep `assemble.py` responsible for passing resolved geometry to render groups.
- Keep render groups such as `TickGroup` responsible only for drawing with the geometry they receive.

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

Also stop relying on temporarily mutating `canvas_config.circular_track_preset` for tick rendering. `add_tick_group_on_canvas()` already accepts `track_preset`, so pass the resolved preset directly:

```python
if use_slot_tick_options and "preset" in resolved_slot.params:
    tick_group_kwargs["track_preset"] = str(resolved_slot.params["preset"])
```

Expected effect:

- The tick renderer receives `label_side`, `tick_side`, `tick_length_px`, and `track_preset` from the resolved slot.
- The rendered tick footprint matches the resolver's measured footprint.
- No new public option is required.

### 2. Distinguish Hard Anchors From Soft Preset Anchors

Update `resolve_circular_radial_layout()` and the private radial intent model in `gbdraw/diagrams/circular/radial_layout.py`.

Preset-generated numeric tracks currently have explicit radius values derived from `track_dict`. The first resolver pass treats these as pinned anchors, so a collision with an earlier tick reserved band only emits a warning.

The current `explicit_anchor` flag mixes two meanings:

- The slot input had a radius.
- The radius must not move.

Split those semantics internally:

- `explicit_anchor`: the slot input had a radius.
- `generated_by_preset`: the slot came from normal preset expansion, not from the user slot API.
- `hard_anchor`: a radius that must not move.
- `soft_anchor`: a preset numeric/depth radius that should be tried first, but may move if it collides.

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

Expected effect:

- User-pinned slots remain pinned and retain current validation behavior.
- Preset numeric/depth slots are no longer allowed to overlap ticks or labels just because their legacy radius was explicit.
- The code documents the difference between "this slot had a radius" and "this slot is immovable".

### 3. Track Preset Origin Outside Public Params

Do not rely on `_preset_generated` inside `CircularTrackSlot.params` for behavior that affects user slot semantics. `params` is public-facing enough that a user can accidentally or deliberately pass the same key.

Instead, pass preset origin as resolver-only metadata from `add_record_on_circular_canvas()`.

When normal preset mode expands slots:

```python
layout_slots = circular_track_slots_for_preset(...)
preset_slot_ids = {slot.id for slot in layout_slots}
```

When user slot mode is active:

```python
preset_slot_ids = set()
```

Then call:

```python
radial_layout = resolve_circular_radial_layout(
    ...,
    preset_slot_ids=preset_slot_ids,
)
```

Inside `resolve_circular_radial_layout()`, pass `preset_slot_ids` into `_slot_intents()` and set a private `_RadialSlotIntent.generated_by_preset` field.

This keeps user-provided `CircularTrackSlot(radius=...)` pinned even if its `params` include `_preset_generated`.

### 4. Place Preset Numeric Tracks as Preferred, Compressible Groups

Preset-generated inside numeric/depth tracks should use their `track_dict` radii as preferred anchors, not hard pins.

Handle contiguous soft preset numeric/depth slots as one stack, for example:

```text
depth -> gc_content -> gc_skew
```

The group-level flow should be:

1. Measure the whole group at its preferred legacy radii and full widths.
2. Accept those placements if every reserved band fits the current occupied bands and the numeric bands satisfy slot spacing/order constraints.
3. If the preferred group collides, fall back to `_place_inside_auto_group()`.
4. During fallback, allow compression using the existing `_candidate_widths()` path.
5. Never compress below `_min_readable_numeric_width_px()`.
6. If the readable minimum still cannot fit, raise `ValidationError` instead of overlapping the tick label, definition, or other reserved bands.

Implementation notes:

- Do not set preset numeric slots to `compress=True` in the public slot object just to unlock this behavior.
- Add an internal helper or placement flag so soft preset numeric fallback can use compression candidates without changing user-visible slot semantics.
- Keep full-width, legacy-radius output unchanged whenever the preferred group fits.
- Place the group together so one numeric track does not consume space needed by the next numeric track.

Expected effect:

- Existing non-overlapping diagrams keep their legacy numeric track radii.
- EDL933-like diagrams move GC/skew/depth inward only when needed.
- Tight layouts can reduce numeric/depth widths down to readable minimums rather than failing too early.
- If no readable layout exists, the resolver fails explicitly instead of producing hidden ticks or labels.

### 5. Preserve User Slot Semantics

Do not relax explicit anchors provided by users.

Only resolver-known preset slots should become soft. User-provided `CircularTrackSlot(radius=...)` should continue to behave as pinned geometry. If it conflicts and `strict=True`, it should still raise. If it conflicts and `strict=False`, it may warn as today, but it must not silently move.

This keeps the public slot API predictable.

### 6. Keep Drawing Order Stable Unless Tests Prove Otherwise

Do not fix the issue by drawing ticks last.

Drawing ticks last would hide the symptom while the layout still contains overlapping reserved bands. The resolver should produce non-overlapping radial bands first. Drawing order can remain based on `(z, slot_index)` unless a separate layering bug is found.

## Test Plan

Add or update focused tests in `tests/test_circular_track_slots.py` and/or `tests/test_circular_radial_layout.py`.

Required coverage:

1. Normal preset tick rendering receives resolved tick options.
   - Verify `label_side`, `tick_side`, `tick_length_px`, and `track_preset` are passed to `add_tick_group_on_canvas()` outside user slot mode.
   - Verify `canvas_config.circular_track_preset` mutation is no longer needed for the tick render path.

2. EDL933 `tuckin` layout avoids tick/numeric overlap.
   - Use EDL933 with `track_type="tuckin"`, `strandedness=True`, `show_gc=True`, and `show_skew=True`.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick path band.
   - Assert GC content and GC skew draw/reserved bands do not overlap the tick label band when labels exist.
   - Include configured slot spacing in overlap assertions.

3. Depth/GC/skew stacks avoid all relevant reserved bands.
   - Include a resolver-level test with depth, GC content, and GC skew enabled.
   - Assert depth, GC content, GC skew, ticks, tick labels, feature bands, and the center definition band do not overlap.
   - Assert numeric tracks do not overlap each other and preserve ordering/spacing.

4. Legacy radii are preserved when no overlap exists.
   - Use an existing smaller/non-conflicting circular record.
   - Assert depth/GC/skew centers match the previous `track_dict`-derived radii.
   - Assert widths remain full default widths.

5. User-pinned numeric slots remain pinned.
   - Provide an explicit radius in a user slot.
   - Verify the resolver does not silently move it.
   - If the slot is strict and conflicts, verify it still raises.
   - Add a regression case where the user slot has `params={"_preset_generated": True}` and verify it is still treated as user-pinned.

6. Soft preset numeric fallback can compress.
   - Construct a tight inside layout where preferred full-width numeric slots collide but a compressed readable stack fits.
   - Assert at least one numeric width shrinks.
   - Assert no numeric width drops below `_min_readable_numeric_width_px()`.

7. Soft preset numeric fallback fails when no readable layout exists.
   - Construct a layout where even readable minimum widths cannot fit.
   - Assert `ValidationError` is raised.
   - Assert the failure is not converted into a warning-only overlap.

## Acceptance Criteria

- The EDL933 screenshot scenario no longer hides tick marks or tick labels under GC content.
- Normal circular presets and custom slot mode both use the same resolved tick geometry.
- Preset numeric/depth tracks preserve legacy radii when they fit.
- Preset numeric/depth tracks move and compress only when required to avoid overlap.
- User-pinned numeric/depth slots remain pinned.
- No new public configuration is added.
- Existing non-overlapping diagrams remain visually stable.
- Relevant circular track tests pass.

## Files Expected to Change

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- Possibly `gbdraw/diagrams/circular/presets.py` if `_preset_generated` is removed from preset slot params.
- `tests/test_circular_track_slots.py`
- `tests/test_circular_radial_layout.py`

Reference SVG files should only be updated if intentional output changes are observed in snapshot/reference tests.
