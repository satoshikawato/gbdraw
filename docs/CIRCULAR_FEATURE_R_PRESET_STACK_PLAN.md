# Circular Feature-R Preset Stack Plan

## Status

Planned.

## Purpose

Simplify circular track layout by redefining `tuckin`, `middle`, and
`spreadout` as presets for the feature stack radius, then stacking every other
track around the resolved feature band.

The intended result is that the web UI no longer needs to group Custom Track
Slots into "Outer tracks", "On-axis tracks", and "Inner tracks". Users should
see one ordered track stack, while the layout engine computes the final radial
placement from a smaller set of rules.

## Current Problem

Circular layout is currently governed by overlapping concepts:

- `track_type` selects `tuckin`, `middle`, or `spreadout`.
- Preset expansion turns that `track_type` into concrete slot geometry.
- Feature layout also has `lane_direction=inside|split|outside`.
- Track slots expose `side=inside|overlay|outside`.
- The radial resolver branches by `side` and then packs tracks separately.
- The web UI mirrors `side` as three sections:
  "Outer tracks", "On-axis tracks", and "Inner tracks".

This is functional, but it makes the model hard to explain. The same visual
intent can be expressed through preset, lane direction, side, and radius.

## Target Model

Use one primary geometric primitive for circular presets:

```text
feature.r = center radius of the resolved feature stack
```

The feature stack means all feature lanes after strand separation and overlap
resolution have been considered. It is not just the center radius of track 0.

Let:

```text
R = circular genome axis radius in px
w = feature lane width in px
g = lane/track spacing in px
n = resolved number of feature lanes
B = feature stack band width in px

B = (n * w) + ((n - 1) * g)
```

Then the presets are:

```text
tuckin:
  feature stack outer edge is inside the axis
  feature_center_px = R - g - (B / 2)
  feature.r = feature_center_px / R

middle:
  feature stack center is on the axis
  feature_center_px = R
  feature.r = 1.0

spreadout:
  feature stack inner edge is outside the axis
  feature_center_px = R + g + (B / 2)
  feature.r = feature_center_px / R
```

The feature band is then:

```text
feature_inner_px = feature_center_px - (B / 2)
feature_outer_px = feature_center_px + (B / 2)
```

Other tracks are placed relative to `feature_inner_px`,
`feature_outer_px`, the fixed axis radius, and the reserved center definition
band.

## Preset Semantics By Feature Mode

### `separate_strands=false`, `resolve_overlaps=false`

There is one feature lane.

```text
n = 1
B = w
```

Preset behavior:

- `tuckin`: one feature lane inside the axis.
- `middle`: one feature lane centered on the axis.
- `spreadout`: one feature lane outside the axis.

### `separate_strands=true`, `resolve_overlaps=false`

There are two feature lanes, one for each strand.

```text
n = 2
B = (2 * w) + g
```

Preset behavior:

- `tuckin`: both strand lanes are inside the axis.
- `middle`: positive strand lane is outside, negative strand lane is inside.
- `spreadout`: both strand lanes are outside the axis.

The lane order inside the stack must remain deterministic:

```text
outer side first: positive, then negative when both are outside
inner side first: negative, then positive when both are inside
middle split: negative inside, positive outside
```

Exact lane ordering can be implemented in one helper and should not be spread
across drawers, label code, and the resolver.

### `separate_strands=false`, `resolve_overlaps=true`

The number of lanes depends on overlap resolution.

```text
n = number of overlap-resolution lanes
B = (n * w) + ((n - 1) * g)
```

Preset behavior:

- `tuckin`: all overlap lanes stack inward from the axis-side boundary.
- `middle`: lane 0 is centered on the axis; overflow lanes expand around it.
- `spreadout`: all overlap lanes stack outward from the axis-side boundary.

For `middle`, preserve the current useful behavior where track 0 remains the
visual anchor and overflow lanes may split across the axis. The important
change is that this is now derived from `feature.r=1.0` plus a feature stack
policy, not from a public preset/lane/side triad.

### `separate_strands=true`, `resolve_overlaps=true`

Positive and negative strand pools may each create multiple lanes.

```text
n = positive_lane_count + negative_lane_count
B = (n * w) + ((n - 1) * g)
```

Preset behavior:

- `tuckin`: both strand lane stacks are inside the axis.
- `middle`: positive strand lanes stack outside, negative strand lanes stack
  inside.
- `spreadout`: both strand lane stacks are outside the axis.

The resolver should measure the resulting full feature band before placing
numeric, depth, tick, spacer, or custom tracks.

## Track Stacking Rules

After resolving the feature band, treat it as an occupied band. Additional
tracks are stacked around it.

### Outer Stack

Outer tracks are placed outward from:

```text
outer_cursor_px = max(R, feature_outer_px) + spacing
```

Each outer track consumes its resolved reserved band and advances:

```text
outer_cursor_px = resolved_reserved_outer_px + spacing
```

If the outer stack exceeds the current canvas radius budget, the canvas may
expand as it does today for outer content.

### Inner Stack

Inner tracks are placed inward from:

```text
inner_cursor_px = min(R, feature_inner_px) - spacing
```

Each inner track consumes its resolved reserved band and advances:

```text
inner_cursor_px = resolved_reserved_inner_px - spacing
```

The inner stack must also avoid the center definition reserved band. If a track
cannot fit, the resolver should either apply explicit compression rules or
raise a clear validation error.

### On-Axis Tracks

The axis itself is fixed and should not be represented as a mutable track slot.

Ticks may still anchor to the axis or to a measured tick band, but this should
be a renderer-specific property rather than a public "On-axis tracks" UI
section. In the final UI, tick placement can be shown as an advanced field
inside the tick row.

## Custom Track Behavior

The preset should not move when users add custom tracks. The preset defines the
feature stack first; custom tracks are then stacked around that resolved band.

### Adding An Outer Custom Track

Given:

```text
feature_inner_px
feature_outer_px
```

An outer custom track is placed after the feature band:

```text
custom_center_px >= feature_outer_px + spacing + (custom_width / 2)
```

Effects by preset:

- `tuckin`: the feature stack stays inside; the custom outer track usually
  lives outside the axis.
- `middle`: the feature stack stays around the axis; the custom outer track
  starts outside the feature band.
- `spreadout`: the feature stack stays outside; the custom outer track starts
  even farther outside, possibly causing canvas expansion.

### Adding An Inner Custom Track

An inner custom track is placed before the feature band:

```text
custom_center_px <= feature_inner_px - spacing - (custom_width / 2)
```

Effects by preset:

- `tuckin`: the feature stack is already inside, so inner custom tracks compete
  with center definition space and may fail sooner.
- `middle`: inner custom tracks use the inner side of the feature band.
- `spreadout`: inner custom tracks can often use the large space between the
  axis/center and the outer feature stack.

### Ordering

The web UI should expose one ordered list. Each row may have an advanced stack
side, but the visual grouping into three sections should be removed.

Recommended public/internal names:

```text
stack_side = outer | inner | auto
```

`side=outside|inside|overlay` can remain as a compatibility field during the
transition, but new UI copy should avoid presenting it as three separate
track-slot categories.

Suggested compatibility mapping:

```text
side=outside -> stack_side=outer
side=inside  -> stack_side=inner
side=overlay -> renderer-specific anchored placement
```

## Web UI Plan

### Remove Three Section UI

Change `gbdraw/web/js/app/circular-track-slots.js` and
`gbdraw/web/index.html` so Custom Track Slots render as a single ordered list.

Remove or retire:

- `SECTION_DEFINITIONS`
- `SECTION_RANKS`
- `circularTrackSlotSections()`
- section-specific add buttons
- section-specific move rules

Add:

- one `circularTrackSlots()` list preserving user order
- one "Add track" control
- one move up/down control that moves across the full list
- optional per-row advanced placement/stack-side control

### Preset Buttons

`Apply Tuckin`, `Apply Middle`, and `Apply Spreadout` should set
`form.track_type` and reset default custom rows without storing preset-derived
fixed px geometry.

The preset summary should describe feature stack geometry, for example:

```text
Tuckin preset: feature stack inside axis
Middle preset: feature stack centered on axis
Spreadout preset: feature stack outside axis
```

When enough preview context is available, it may also show approximate
`feature.r`.

### Saved State

Saved sessions should store:

- `form.track_type`
- ordered slots
- explicit user overrides

Saved sessions should not store record-derived preset values such as measured
tick radius, computed feature stack radius, or px widths unless the user typed
those values.

## Python Layout Plan

### 1. Introduce Feature Stack Measurement

Add a small model in `gbdraw/diagrams/circular/radial_layout.py` or a focused
helper module:

```python
@dataclass(frozen=True)
class CircularFeatureStackMetrics:
    lane_count: int
    lane_width_px: float
    lane_spacing_px: float
    band_width_px: float
    center_radius_px: float
    inner_radius_px: float
    outer_radius_px: float
    lane_centers_by_track_id: Mapping[int, float]
```

The helper should derive lane centers from:

- preset name
- `separate_strands`
- `resolve_overlaps`
- actual `feature_track_id` values
- feature lane width
- spacing
- axis radius

### 2. Redefine Preset Expansion

Update `gbdraw/diagrams/circular/presets.py` so
`circular_feature_slot_defaults_for_preset()` computes feature slot geometry
from the target feature stack rule rather than from `lane_direction` as the
main preset meaning.

During transition, feature slots may still carry a renderer param that tells
legacy feature drawing how lanes should be ordered. That param should be
derived from the new feature stack policy.

### 3. Update Feature Layout

Update `build_circular_feature_layout()` so `anchor_radius_px` represents the
feature stack center. The function should:

- compute full stack `B`
- derive `feature_inner_px` and `feature_outer_px`
- assign lane centers within that band
- stop using preset names as the source of lane placement when explicit stack
  geometry is available

The old `track_type` fallback can remain for compatibility until all callers
pass explicit feature stack policy.

### 4. Update Radial Resolver

The resolver should:

- resolve and reserve the feature band first
- reserve tick bands when they are renderer-anchored
- place remaining outer-stack tracks outward from the feature/axis outer cursor
- place remaining inner-stack tracks inward from the feature/axis inner cursor
- avoid center definition space
- preserve draw order as `(z, slot_index)`

This should reduce the amount of branching on public `side`.

### 5. Keep Compatibility During Migration

Do not remove `CircularTrackSlot.side` in the first implementation pass. It is
public API and many tests/docs currently use it.

Instead:

- keep parsing `side`
- map `side` to stack-side intent
- keep `overlay` only for renderer-specific anchored behavior
- update web UI to avoid exposing the old three-section mental model

Once tests and docs use the new language, consider a later cleanup that
renames the public field or deprecates `side`.

## Files Expected To Change

- `gbdraw/diagrams/circular/presets.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/features/tracks.py` if feature track ids need clearer lane grouping
- `gbdraw/render/drawers/circular/features.py`
- `gbdraw/render/groups/circular/seq_record.py`
- `gbdraw/render/groups/circular/ticks.py`
- `gbdraw/render/groups/circular/gc_content.py`
- `gbdraw/render/groups/circular/gc_skew.py`
- `gbdraw/render/groups/circular/depth.py`
- `gbdraw/labels/circular.py`
- `gbdraw/tracks/circular.py`
- `gbdraw/web/index.html`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/state.js`
- `tests/test_circular_track_slots.py`
- `tests/test_circular_radial_layout.py`
- `tests/test_circular_feature_width.py`
- `tests/test_depth_track.py`
- `tests/test_circular_label_placement.py`
- `tests/test_web_packaging.py`
- circular reference SVGs if output changes intentionally

## Test Plan

### Unit Tests

Add focused tests for feature stack metrics:

- one lane, no overlap resolution
- two lanes from `separate_strands=true`
- multiple lanes from `resolve_overlaps=true`
- split middle behavior
- computed `feature.r` for all presets
- custom inner/outer stacking around the feature band

### Resolver Tests

Cover:

- default `tuckin`, `middle`, `spreadout`
- each combination of `separate_strands` and `resolve_overlaps`
- custom outer track after each preset
- custom inner track after each preset
- center-definition collision with inner custom tracks
- canvas expansion with spreadout plus outer custom tracks
- depth track reserved axis footprint
- tick labels anchored near the axis

### Web Tests

Update packaging/static tests to assert:

- no `Outer tracks`, `On-axis tracks`, or `Inner tracks` section labels
- one ordered Custom Track Slots list exists
- preset buttons still exist
- slot order is preserved
- generated CLI args omit preset-derived geometry unless explicitly edited
- saved config/session schema accepts the new representation

### Reference Output Tests

This change will likely alter circular geometry. For every reference SVG
change, document whether it is intentional:

- feature stack shifted by new preset definition
- custom slots stacked relative to the feature band
- tick or label changes caused by new feature anchor radius

## Implementation Phases

### Phase 1: Internal Feature Stack Model

- Add feature stack metric helper.
- Add tests for metric calculation.
- Keep current rendering behavior unchanged where possible.

### Phase 2: Preset Reinterpretation

- Redefine preset feature geometry as `feature.r`.
- Update feature layout to use feature stack center and band width.
- Keep compatibility params for old feature drawing paths.

### Phase 3: Resolver Stacking

- Make feature band the primary blocker.
- Stack other tracks inward/outward relative to feature band.
- Preserve current validation behavior for impossible layouts.

### Phase 4: Web UI Simplification

- Replace three sections with one ordered list.
- Add per-row advanced stack-side control only where needed.
- Update config/session normalization and CLI arg generation.

### Phase 5: Tests And Documentation

- Update unit and web tests.
- Regenerate reference outputs if geometry changes intentionally.
- Update user-facing docs and CLI reference.

## Open Questions

- Should `middle + resolve_overlaps=true + separate_strands=false` keep the
  current strand-aware split behavior, or should it become a pure symmetric
  stack independent of strand?
- Should default custom tracks use `stack_side=auto`, or should renderer
  defaults choose inner/outer?
- Should `ticks` be modeled as a normal stack item, or remain an axis-anchored
  renderer with optional label/tick side params?
- Should `side` be deprecated publicly after the web UI stops using the
  three-section model?

## Non-Goals

- Do not remove `CircularTrackSlot.side` in the first pass.
- Do not store record-specific computed preset geometry in web sessions.
- Do not silently move user-authored tracks to the opposite stack when they do
  not fit.
- Do not change linear diagram track layout.
- Do not remove existing CLI/API geometry overrides until an explicit
  compatibility decision is made.
