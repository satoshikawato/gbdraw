# Track Spacing Unification Plan

## Status

Planned.

## Goal

Unify user-facing track and label spacing behavior so explicit pixel values mean
the same thing across the web UI, CLI-generated slot strings, and Python layout
code.

The target behavior is:

- Web UI spacing controls use pixels only.
- Explicit user spacing values are honored exactly.
- Automatic/default spacing may still be adjusted by layout code when needed.
- Circular track gaps use physical radial sides: `inner_gap_px` and
  `outer_gap_px`.
- Circular label spacing no longer weakens explicit values through legacy
  dampening.

## User-Facing Principles

### 1. Web UI Uses Pixels

The web UI should display and serialize spacing values in px. For circular
slots, this avoids exposing the lower-level `ScalarSpec` behavior where
unitless values and percentages are resolved against the base circular radius.

Examples emitted by the web UI should use explicit px values:

```text
spacing=4px
inner_gap_px=4
outer_gap_px=6
```

The CLI/API may continue accepting the existing circular `ScalarSpec` grammar
for compatibility, but the UI should not generate ambiguous unitless circular
spacing values.

### 2. Explicit Values Are Strict

When a user explicitly sets a track or label spacing value, layout code must
either honor that value or raise a validation error if the requested layout
cannot fit.

Explicit spacing must not be silently compressed, scaled down, or partially
applied.

Automatic or preset-derived spacing may still be adjusted by layout code. That
keeps dense automatic layouts usable while making manual settings predictable.

## Circular Track Gaps

Circular track slots should expose physical radial gaps:

- `inner_gap_px`: gap on the side closer to the circle center.
- `outer_gap_px`: gap on the side farther from the circle center.

This definition is independent of `side=inside` or `side=outside`.

For an outside track:

```text
center -> inner side -> track -> outer side
```

For an inside track:

```text
center -> inner side -> track -> outer side -> axis
```

Therefore, for `side=inside`, the axis-adjacent side is the physical outer side.
This is expected and should be reflected in UI wording if an axis-relative
control is ever added.

### Relationship To Existing `spacing`

The existing circular slot `spacing` field is ambiguous because it acts as a
single packing gap rather than a side-specific physical gap.

Migration target:

- New UI controls should use `inner_gap_px` and `outer_gap_px`.
- Side-specific gap fields with the `_px` suffix should accept numeric pixel
  values only.
- New internal resolved geometry should carry both values explicitly.
- Existing `spacing` may remain as a compatibility alias during parsing.
- If `spacing` is supplied alone, it should map to both `inner_gap_px` and
  `outer_gap_px` unless a stricter compatibility rule is required by existing
  tests.
- If both `spacing` and a side-specific gap are supplied, the side-specific
  value should win for that side, or validation should reject the ambiguous
  combination. Prefer rejection for new grammar.

### Adjacent Track Clearance

When two circular tracks are adjacent, the required gap between them should be
computed from the physical sides that face each other.

For a track closer to the center followed by a track farther from the center:

```text
required_gap = max(inner_track.outer_gap_px, outer_track.inner_gap_px)
```

This rule should be applied consistently to pinned and auto-placed tracks.

## Linear Track Gaps

Linear mode may use axis-relative or vertical naming in the UI, but it should
still serialize px values.

Preferred public naming for linear-specific controls:

- `upper_gap_px`: gap above the track.
- `lower_gap_px`: gap below the track.

If the implementation needs mode-neutral internal fields, use `gap_before_px`
and `gap_after_px` only inside the linear resolver, where ordering is explicit.
Do not expose `before`/`after` as circular radial terminology.

## Label Spacing

Circular label spacing is separate from track spacing. Track gaps control radial
layout of rings; label spacing controls label-to-label collision avoidance.

Explicit circular label spacing must be honored exactly. The current conservative
legacy dampening behavior should be removed.

Current behavior to retire:

```python
legacy_spacing_px + ((requested_spacing_px - legacy_spacing_px) * 0.25)
```

That behavior makes a user-specified value only partially effective. For
example, if the legacy spacing is `2px` and the user requests `10px`, the layout
uses only `4px`. That is not acceptable for explicit user settings.

Replacement behavior:

- If a user explicitly sets label spacing, use that exact px value.
- If label spacing is omitted, use the default configuration or an automatic
  layout-derived value.
- If the requested explicit label spacing cannot be satisfied, fail with a clear
  validation error instead of silently reducing it.

## Implementation Notes

Likely affected areas:

- `gbdraw/tracks/circular.py`
- `gbdraw/tracks/linear.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/linear/track_slots.py`
- `gbdraw/labels/circular.py`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/linear-track-slots.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/state.js`
- CLI reference and customization tutorials

Testing should cover:

- Web-generated circular slot strings always include px units for spacing/gaps.
- Circular `inner_gap_px` and `outer_gap_px` affect the correct physical radial
  sides for both inside and outside tracks.
- Explicit circular gaps are not compressed in dense inside auto layouts.
- Existing auto/default circular layouts may still compress only non-explicit
  gaps.
- Explicit circular label spacing is applied exactly.
- Reference SVG changes caused by removing label spacing dampening are reviewed
  and updated intentionally.
