# Circular Track Slots Root Fix Plan

## Purpose

This plan describes the root fix for the current Custom Track Slots bugs:

- Enabling Custom Track Slots can make a `tuckin` diagram look closer to
  `spreadout`.
- GC skew, especially with depth/custom slots, can be packed into the center
  and overlap the circular definition text.

The fix should not be a narrow patch for one renderer. The underlying problem
is that the custom slot path and the legacy circular layout path do not share
the same layout constraints. Custom slots currently enter an independent
auto-packing path, while the legacy path has separate logic for `track_type`,
depth compression, tick label bounds, feature-band avoidance, and definition
space.

## Root Cause

### 1. Custom slots bypass legacy circular semantics

The web UI builds default custom slots without explicit radius or annulus
placement. When these slots are passed to the Python API, `assemble.py` enters
slot mode and calls `resolve_circular_track_slots(..., compatibility_mode=False)`.

In non-compatibility mode the resolver uses legacy centers only as preferred
anchors, then repacks tracks to avoid measured slot collisions. That means a
default custom slot list is no longer equivalent to the normal legacy tuckin
diagram.

The user-visible result is surprising: turning on Custom Track Slots can change
the visual meaning of `track_type=tuckin`, even when the slot list is just the
default feature/tick/GC/skew stack.

### 2. Definition space is not part of slot packing

The slot resolver prevents slot-to-slot radial overlap, but it does not reserve
the circular definition area around the center. It only rejects tracks whose
inner radius goes below zero.

This allows auto-packed inner tracks, especially `gc_skew`, to move too far
inward. With depth enabled, the default slot order can become:

```text
features, ticks, depth, gc_content, gc_skew
```

The final `gc_skew` slot can be packed near the center because the resolver has
no definition guard.

### 3. Depth/GC/skew compression is legacy-only

The legacy path has special logic to compress GC content/skew when depth is
enabled so the center definition area remains usable. Slot mode bypasses that
logic and instead treats depth, GC content, and GC skew as ordinary fixed-width
annuli.

The duplicated logic is the main architecture smell. Any future custom slot
renderer would have to rediscover the same constraints unless the constraints
move into a shared layout contract.

## Design Goal

Custom Track Slots should extend the circular layout system, not replace it.

The target behavior is:

1. `track_type` always keeps its meaning.
2. Enabling Custom Track Slots with the default built-in slots should not change
   the diagram.
3. Additional custom slots should pack predictably while respecting:
   - feature band footprints,
   - tick labels and axis stroke,
   - depth/GC/skew compression constraints,
   - center definition space,
   - explicit user radius/annulus pins.
4. Draw order should remain independent from radial packing order.

## Implementation Strategy

Implement the root fix in three layers:

1. Make slot layout track-type aware and legacy-compatible.
2. Add first-class reserved bands, including definition space.
3. Move depth/GC/skew compression rules into slot layout so legacy and custom
   paths use the same constraints.

Do not solve this only in the web UI. The web UI can improve defaults, but the
Python API and CLI must remain correct when users pass `circular_track_slots`
directly.

## Phase 1: Classify Built-In Legacy-Compatible Slot Sets

### Goal

Default Custom Track Slots should be behaviorally equivalent to the normal
legacy diagram.

### Add helper functions

Add these helpers near circular slot orchestration, likely in
`gbdraw/diagrams/circular/assemble.py` or a new small module under
`gbdraw/diagrams/circular/`:

```python
def _is_default_legacy_slot_stack(
    slots: Sequence[CircularTrackSlot] | None,
    *,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
) -> bool: ...

def _slot_has_custom_geometry(slot: CircularTrackSlot) -> bool: ...

def _slot_has_custom_render_params(slot: CircularTrackSlot) -> bool: ...
```

The stack is legacy-compatible only when all of these are true:

- Slot IDs and renderers match the expected built-in sequence:
  - `features`
  - `ticks`
  - optional `depth`
  - optional `gc_content`
  - optional `gc_skew`
- No slot has explicit radius, inner radius, outer radius, width, gap, or
  non-default `z`.
- No slot has custom legend label or non-default tick side/label side beyond
  the current legacy tick params.
- Dinucleotide params match the high-level `dinucleotide` option.
- No duplicate renderer slots exist.

### Behavior

When the slot stack is legacy-compatible:

- Treat the diagram as legacy layout, not custom auto-pack.
- Still honor the slot-enabled show/hide intent.
- Preserve the current public result that passing slots controls which built-in
  tracks are shown.

There are two possible implementations:

1. Convert `parsed_circular_track_slots` to `None` after deriving
   `show_depth/show_gc/show_skew`.
2. Keep `slot_mode=True`, but pass `compatibility_mode=True` and ensure drawing
   uses exactly the same geometry as the legacy path.

Prefer option 1 for default stacks because it exercises the existing mature
legacy path, including depth compression and definition handling. Use option 2
only if preserving slot IDs for built-in default stacks is required.

### Tests

Add tests in `tests/test_circular_track_slots.py`:

- `test_default_custom_slots_match_legacy_tuckin_svg`
- `test_default_custom_slots_match_legacy_middle_svg`
- `test_default_custom_slots_match_legacy_spreadout_svg`
- `test_default_custom_slots_with_depth_use_legacy_compressed_gc_skew_layout`

The tests should compare relevant geometry rather than full SVG text where
possible:

- GC content group transform/radius.
- GC skew group transform/radius.
- Tick group radius.
- Definition group position.

## Phase 2: Add Definition Reserved Band To Slot Layout

### Goal

When true custom packing is needed, no auto-placed slot may enter the center
definition area.

### Extend layout context

Extend `CircularTrackLayoutContext` in
`gbdraw/diagrams/circular/slot_layout.py`:

```python
@dataclass(frozen=True)
class CircularTrackLayoutContext:
    ...
    reserved_bands_px: tuple[tuple[float, float], ...] = ()
    min_auto_inner_radius_px: float | None = None
```

Use `reserved_bands_px` for hard forbidden annuli. Use
`min_auto_inner_radius_px` as a simple guard for auto-placed slots when a full
band is not necessary.

### Measure definition bounds

Add a helper in `assemble.py`:

```python
def _definition_reserved_radius_px(
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    species: str | None,
    strain: str | None,
    config_dict: dict,
    *,
    cfg: GbdrawConfig,
    plot_title: str | None,
    definition_profile: str,
) -> float:
    ...
```

Implementation options:

- Best option: instantiate `DefinitionGroup`, measure its local bounds, and
  convert half width/height into a conservative radius plus padding.
- Simpler initial option: derive from font size, number of definition lines,
  and maximum text width using existing font metrics helpers.

The reserved area should be conservative:

```text
definition_reserved_band = (0, definition_radius_px + padding_px)
```

Reasonable padding:

```text
max(8 px, 0.02 * base_radius)
```

### Apply during packing

In `resolve_circular_track_slots()`:

- Add reserved bands to `occupied` before auto slots are packed.
- Do not move pinned slots. If a pinned slot overlaps a reserved band, emit a
  warning.
- Auto slots must avoid reserved bands the same way they avoid other occupied
  slots.
- If an auto slot cannot fit above the definition guard, raise a clear error:

```text
circular track slot 'gc_skew' cannot fit without overlapping the center definition; reduce width, remove inner slots, or set an explicit radius.
```

### Tests

Add resolver-level tests:

- Auto-packed GC skew does not enter `(0, reserved_radius)`.
- Pinned slot overlapping definition emits a warning but is not moved.
- Packing raises a helpful error when there is not enough radial space.

Add integration tests:

- Depth + GC + skew + one duplicate skew avoids definition.
- Same test with `track_type=tuckin`.
- Same test with long definition text.

## Phase 3: Share Depth/GC/Skew Compression With Slot Mode

### Goal

The custom slot path should not regress the legacy behavior that compresses
GC/skew tracks to preserve center space when depth is present.

### Refactor current compression logic

Current legacy helpers include:

- `_default_gc_skew_layout_without_depth`
- `_default_gc_skew_gap_px`
- `_resolve_depth_compressed_gc_skew_layout`

Move or adapt this logic so it can produce preferred footprints for slot mode.

Add a helper:

```python
def _resolve_builtin_numeric_slot_preferred_layouts(
    *,
    slots: Sequence[CircularTrackSlot],
    canvas_config: CircularCanvasConfigurator,
    cfg: GbdrawConfig,
    depth_inner_radius_px: float | None,
    definition_reserved_radius_px: float | None,
    track_specs: Sequence[TrackSpec] | None,
) -> dict[str, tuple[float, float]]:
    """Return {slot_id: (center_px, width_px)} for built-in depth/content/skew slots."""
```

This should:

- Preserve explicit user radius/width.
- Apply depth compression only to built-in depth/content/skew slots.
- Keep duplicate custom GC/skew slots in normal auto-packing unless they have
  explicit geometry.
- Respect `definition_reserved_radius_px`.

### Integrate with resolver

There are two viable approaches:

1. Pre-annotate slot placements before calling the resolver.
2. Add `preferred_layouts_by_slot_id` to `CircularTrackLayoutContext`.

Prefer option 2:

```python
preferred_layouts_px: Mapping[str, tuple[float, float]] = field(default_factory=dict)
```

In `resolve_circular_track_slots()`:

- If a slot has no explicit radius/annulus and appears in
  `preferred_layouts_px`, use that center/width as the preferred layout.
- Treat it as auto-movable unless the user explicitly pinned radius/annulus.
- If it conflicts with reserved bands, allow normal packing to move it.

This keeps compression as a preference, not an unsafe pin.

### Tests

Add integration tests that cover:

- Legacy depth compression and custom default slots produce the same GC/skew
  positions.
- True custom slots still avoid compressed built-in slots.
- Explicit `gc_skew@r=...,w=...` is preserved even if it overlaps, with a
  warning rather than silent movement.

## Phase 4: Make Feature And Tick Footprints Authoritative

### Goal

The resolver should pack against the real feature/tick footprint, not a guessed
center plus width.

### Features

The existing slot layout already has `feature_band_offsets_px`, but it should
be treated as authoritative:

- Use `precomputed_feature_dict` whenever slot mode is active.
- Compute full feature-band bounds after overlap resolution and optional
  feature width override.
- Store these bounds in `CircularSlotFootprint`.
- Draw features using `resolved_slot.anchor_radius_px`.

Add tests where:

- `resolve_overlaps=True` expands the feature band and custom GC/skew avoid the
  expanded band.
- `track_type=tuckin` keeps feature bands inward from the anchor.
- `track_type=spreadout` keeps feature bands outward from the anchor.

### Ticks

Tick footprint must include:

- tick marks,
- tick labels,
- optional axis stroke padding.

The current resolver partially measures this. Tighten it so the drawing path
and measurement path use the same tick side, label side, font size, interval,
and track channel override.

Add tests where:

- GC content avoids tick labels.
- `label_side=none` reserves less space than legacy labels.
- `tick_side=inside/outside/both` changes the measured draw annulus.

## Phase 5: Web UI Adjustments

### Goal

The UI should make the new behavior understandable, but correctness must live
in Python.

### Changes

Update `gbdraw/web/js/app/circular-track-slots.js`:

- Keep default generated slots geometry-free.
- Add a small visual state or tooltip explaining:
  - default built-in slots preserve legacy layout,
  - adding duplicate/custom slots enables auto packing.

Update `gbdraw/web/js/app/run-analysis.js`:

- Continue passing `--track_type`.
- Continue passing `--circular_track_slot` when Custom Track Slots is enabled.
- Do not try to encode legacy radii in JS. That would duplicate Python config
  and would drift.

Update config normalization:

- Ensure loaded configs preserve old custom slot arrays.
- No migration should be required for saved configs.

### Tests

Extend `tests/test_web_packaging.py`:

- Confirm web still wires `--track_type` with `--circular_track_slot`.
- Confirm default slot construction remains geometry-free.
- Confirm no CDN/build-step changes are introduced.

## Phase 6: Backward Compatibility Rules

### Public API

Keep these APIs stable:

- `CircularTrackSlot`
- `parse_circular_track_slot`
- `parse_circular_track_slots`
- `default_circular_track_slots`
- `circular_track_slots_from_order`

### Behavior compatibility

Expected changes:

- Default built-in custom slot stacks become closer to legacy output.
- Auto-packed true custom slots may be placed farther outward than before if
  they would overlap definition space.

Preserve:

- Explicit radius/annulus pins.
- Explicit width.
- Duplicate custom slots.
- Slot IDs in generated group IDs for true custom slots.

For pinned overlaps, prefer warnings over automatic movement. A user-specified
radius is an instruction, even if it produces a bad diagram.

## Phase 7: Acceptance Criteria

The fix is complete when all of these hold:

1. `track_type=tuckin` with default Custom Track Slots matches non-custom
   tuckin geometry.
2. `track_type=middle` with default Custom Track Slots matches non-custom
   middle geometry.
3. `track_type=spreadout` with default Custom Track Slots matches non-custom
   spreadout geometry.
4. Depth + GC content + GC skew does not overlap the definition in slot mode.
5. A duplicated custom skew slot is packed without overlapping:
   - features,
   - ticks/tick labels,
   - GC content,
   - depth,
   - definition.
6. Explicitly pinned custom slots are not silently moved.
7. Existing CLI/API tests for track specs and circular output continue passing.
8. Web packaging tests continue passing without introducing a build step.

## Suggested Implementation Order

1. Add tests that reproduce the two reported bugs.
2. Add legacy-compatible slot stack detection.
3. Route default built-in custom stacks through legacy layout.
4. Add definition reserved band to `CircularTrackLayoutContext`.
5. Add resolver tests for reserved bands and pinned overlaps.
6. Refactor depth/GC/skew compression into slot-compatible preferred layouts.
7. Tighten feature and tick footprint measurement tests.
8. Update web tooltip/copy only after Python behavior is correct.
9. Run focused tests:

```bash
python -m pytest tests/test_circular_track_slots.py -q
python -m pytest tests/test_depth_track.py -q
python -m pytest tests/test_circular_feature_width.py -q
python -m pytest tests/test_web_packaging.py -k "web_run_analysis or circular_track" -q
```

10. Run broader fast circular tests before updating any reference SVGs.

## Notes On Reference Outputs

Reference SVG updates should be deferred until the geometry changes are
reviewed. The main expected reference changes are limited to custom-slot cases
and any newly added fixture outputs. Legacy non-custom circular outputs should
not change. If they do, treat it as a regression unless the change is
intentional and explained.

