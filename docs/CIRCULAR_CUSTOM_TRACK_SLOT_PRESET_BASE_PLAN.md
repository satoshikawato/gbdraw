# Circular Custom Track Slots Preset Base Plan

## Status

Planned.

## Background

Enabling Custom Track Slots in the web UI can turn a valid circular preset into
an invalid custom layout. A typical failing case is EDL933 with the default
`tuckin` layout, separate strands, GC content, and GC skew enabled. The simple
preset renders, but Custom Track Slots can raise:

```text
ValidationError: Circular track slot 'gc_skew' cannot fit inside between 0.0px and 280.5px.
```

The root cause is not the radial resolver itself. The resolver is correctly
rejecting a custom inside stack that cannot fit. The problem is that the web UI
initializes Custom Track Slots with a simplified slot list that is not the same
geometry as the selected circular preset.

Simple circular presets are expanded in Python by
`gbdraw/diagrams/circular/presets.py::circular_track_slots_for_preset()`. That
expansion is record-aware: it uses the record length category, configured
radius, track ratio factors, track ids, tick footprint, strandedness, and
enabled GC/depth/skew tracks. It emits concrete preset slots such as:

```text
features: r=1.00R, w=<preset feature width>
ticks:    r=<measured tick anchor>, w=<measured tick length>
GC:       r=<preset GC radius>, w=<preset GC width>
skew:     r=<preset skew radius>, w=<preset skew width>
```

The web custom editor currently creates default slots in
`gbdraw/web/js/app/circular-track-slots.js` with missing `radius` and `width`
for those same tracks, and `compress=true` on numeric tracks. When sent as
`--circular_track_slot`, these slots enter the custom slot path in
`gbdraw/diagrams/circular/assemble.py`, bypassing Python preset expansion and
the preset preferred-anchor behavior. The layout therefore changes from
"preset geometry that may move/compress if needed" to "auto-pack every numeric
inside track into remaining space".

That difference is enough to fail in center-definition-heavy records, depth
layouts, and multi-record canvases with small records.

## Goals

- Turning on Custom Track Slots must not change a valid simple preset into a
  different auto-packed layout.
- `Reset`, `Apply Tuckin`, `Apply Middle`, and `Apply Spreadout` in the web UI
  must use preset-equivalent geometry semantics.
- Do not bake preset-derived fixed px values into web state, saved sessions, or
  generated `--circular_track_slot` specs.
- Multi-record circular canvas must expand preset geometry per record, not
  reuse one fixed px slot list for every record.
- Explicit user edits must remain authoritative. If a user sets `r`, `w`,
  `spacing`, `side`, order, renderer, or params, that edit must win over the
  inherited preset default.
- Arbitrary invalid custom layouts should still raise a clear
  `ValidationError`; the fix should not silently move user-authored tracks to
  a different side or disable tracks.
- Pure auto custom layout should remain available, but it should be explicit;
  blank web-generated geometry means inherited preset geometry.

## Non-Goals

- Do not make the web editor compute record-specific tick geometry in
  JavaScript.
- Do not store one record's generated px geometry as the saved Custom Track
  Slots state for multi-record diagrams.
- Do not add automatic fallback rules such as "move `gc_skew` outside when it
  fails"; that would make edited layouts unpredictable.
- Do not redesign the Custom Track Slots UI beyond the fields needed to expose
  inherited defaults and explicit overrides.
- Do not change the radial resolver's responsibility: it remains the source of
  truth for final geometry, collision checks, compression, and validation.

## Proposed Design

### 1. Treat Custom Track Slots as Order and Overrides

The corrected model is:

- `track_type` continues to define the default circular preset dimensions.
- Custom Track Slots define track membership, order, and explicit overrides.
- Missing slot geometry does not mean "auto-pack this from scratch"; it means
  "inherit the current preset default for this ordered slot".
- The default custom slot order is the same as the default preset order, so
  turning Custom Track Slots on without edits should render the same dimensions
  as leaving it off.

This is the important distinction:

```text
Bad persisted/generated custom state:
  gc_content:dinucleotide_content@r=0.70,w=74.1px

Good persisted/generated custom state:
  track_type = tuckin
  custom slot order = features, ticks, gc_content, gc_skew
  no explicit r/w unless the user typed an override
```

Custom mode should not force users to specify every radius and width. A slot
row can be only an ordering item until the user changes a field.

### 2. Overlay Slot Order on the Current Preset

Add a helper in `gbdraw/diagrams/circular/presets.py`, for example:

```python
def circular_track_slots_from_preset_order(
    slots: Sequence[CircularTrackSlot],
    preset: str,
    context: CircularPresetContext,
) -> CircularPresetRadialPlan:
    ...
```

The helper should:

- Build the record-specific preset lane templates with
  `circular_track_slots_for_preset(preset, context)`.
- Use the user slot list as the output order.
- For a slot that has no explicit `radius`, use the corresponding preset lane
  radius for its position in the ordered stack.
- For a slot that has no explicit `width`, use the default width for that
  renderer on the current record/canvas.
- For a slot that has no explicit `spacing`, use the preset spacing.
- For renderer params that are not explicitly supplied, inherit preset params:
  - `features`: `lane_direction`
  - `ticks`: `label_side`, `tick_side`, `preset`
  - `dinucleotide_content` / `dinucleotide_skew`: `nt`
  - `depth`: depth renderer defaults as needed
- Preserve explicit user values for `radius`, `width`, `spacing`, `side`,
  `z`, `reserve`, `strict`, `compress`, renderer, and params.
- Return transient layout slots and `preferred_anchor_slot_ids` for inherited
  numeric/depth preset lanes.

The helper must not mutate the original `CircularTrackSlot` objects. Preset
geometry exists only as render-time, record-local defaults.

For default order, the transient layout should be equivalent to the simple
`track_type` preset. For changed order, the slots should occupy the preset lane
sequence in the requested order while keeping renderer-appropriate default
widths. This lets the user change "which track goes outside/inside first"
without manually calculating radii.

### 3. Preserve Pure Explicit Custom Layouts

Some CLI/API users may rely on the current behavior where a geometry-free
custom numeric track is auto-packed. Preserve that by making auto layout an
explicit choice instead of the default interpretation of web-generated rows.

Implementation options:

- Treat web-generated default rows as preset-inheriting rows, and only use
  pure auto placement when the slot explicitly requests it, for example
  `side=inside,compress=true` with no inherited preset lane.
- Or add a small slot-level/internal marker during web generation to distinguish
  "inherited preset geometry" from "pure custom auto".

The preferred behavior for the web UI is preset inheritance. Pure auto remains
available for advanced users, but it is not what Custom Track Slots ON should
mean by default.

### 4. Integrate in Circular Assembly

Update `gbdraw/diagrams/circular/assemble.py`.

Current behavior:

- `circular_track_slots is None`: expand selected `track_type` as a simple
  preset.
- `circular_track_slots is not None`: use supplied slots as a full replacement
  layout.

New behavior:

- `circular_track_slots is None`: unchanged.
- `circular_track_slots is not None`: interpret supplied slots as order and
  explicit overrides on top of `cfg.canvas.circular.track_type`, unless the
  slot explicitly requests pure custom auto geometry.
- Set `layout_slots` to the transient overlay slots.
- Set `preferred_anchor_slot_ids` from inherited preset numeric/depth lanes.
- Keep `canvas_config.circular_track_preset` equal to the active `track_type`.

No new CLI flag is needed for the web fix if Custom Track Slots are redefined
as preset order/overrides. If backward compatibility concerns require a
separate API switch, prefer adding it later after the web behavior is correct.

### 5. Update the Web Custom Track Slots State

Web state should represent order and user overrides:

- Keep `form.track_type` active even when Custom Track Slots are enabled.
- `Reset` rebuilds the slot order from the current simple controls and leaves
  geometry fields blank.
- `Apply Tuckin`, `Apply Middle`, and `Apply Spreadout` update
  `form.track_type` and rebuild the default slot order with blank geometry.
- Do not store preset-derived `radius`, `width`, or `spacing` values.
- Do not send `side`, `strict`, `compress`, or `reserve` unless the user
  explicitly changed them. Default rows should not accidentally become explicit
  geometry overrides.
- If the user enters explicit `radius`, `width`, or `spacing`, send those
  fields and let them override the inherited preset geometry.
- If the user changes order, send the ordered rows with blank geometry.
- If the user changes placement or renderer params, send only the changed
  fields.

This likely requires the web editor to distinguish display defaults from
explicit overrides. For example, the UI may show "Inside (preset)" while the
stored slot `side` remains unset until the user changes placement.

### 6. Update Web Command Generation

Update `gbdraw/web/js/app/run-analysis.js`:

- Continue sending `--track_type form.track_type` in Custom Track Slots mode.
- Continue sending `--circular_track_slot` for the ordered custom rows.
- Build slot specs sparsely: omit inherited geometry and inherited placement.
- Keep suppressing legacy geometry controls such as `--gc_content_width`,
  `--gc_skew_width`, and `--depth_width` while Custom Track Slots are enabled;
  explicit slot fields replace those controls.
- Ensure the default Custom ON command line is equivalent to simple preset mode
  plus an order list, not a full custom auto-pack layout.

Update config/session handling:

- Saved state should round-trip blank inherited geometry fields.
- Old saved configs that contain the current generated defaults with
  `compress=true` and no radius/width should be migrated to blank inherited
  geometry if they match the known default web shape.
- User-authored explicit `compress=true` should remain explicit when it is not
  part of the old default shape.

### 7. Update Documentation

Update user-facing docs after implementation:

- `docs/CLI_Reference.md`: clarify that omitted custom slot geometry inherits
  the selected circular preset unless pure auto geometry is explicitly
  requested.
- `docs/TUTORIALS/1_Customizing_Plots.md`: clarify that Custom Track Slots can
  be used just to reorder tracks.
- `docs/TUTORIALS/3_Advanced_Customization.md`: explain the difference between
  simple presets, preset-inheriting custom order, explicit slot overrides, and
  pure custom auto slots.
- `docs/RECIPES.md`: include a short example of reordering tracks without
  specifying radii or widths.

## Test Plan

### Python Tests

Add or update focused tests in `tests/test_circular_track_slots.py` and
`tests/test_circular_radial_layout.py`.

Required cases:

- EDL933 simple `tuckin` remains valid.
- EDL933 Custom Track Slots ON with default order and blank geometry renders
  successfully and matches the simple `tuckin` preset dimensions.
- EDL933 Custom Track Slots with only slot order changed uses preset-derived
  lane dimensions and does not fall back to full auto packing.
- `middle` and `spreadout` custom order inherit feature lane direction, tick
  side/label side, numeric radii, and widths in the transient Python layout
  slots, without writing those values back to web/custom state.
- Explicit user `r`, `w`, `spacing`, and `side` override inherited preset
  defaults.
- Duplicate extra GC/skew/depth slots that have no preset lane match keep
  advanced custom semantics.
- Depth tracks with inherited preset geometry include the depth reserved radial
  footprint and still avoid center definition overlap.
- Multi-record circular canvas expands inherited preset geometry per record and
  does not reuse one fixed px width for small records.
- Explicit pure-auto slots remain capable of raising the current
  `cannot fit inside` validation error when the requested layout is genuinely
  impossible.

Suggested regression conversion:

- Keep the existing failing EDL933 custom-slots test, but split it:
  - default web-style custom order with blank inherited geometry: succeeds;
  - explicit pure-auto inside GC/skew layout: still raises when unfit.

### Web Packaging Tests

Update `tests/test_web_packaging.py`.

Required assertions:

- Generated `--circular_track_slot` specs do not include preset-derived fixed
  px `r`/`w` values unless the user explicitly typed them.
- Custom Track Slots mode continues to send `--track_type`.
- Default generated slot specs omit inherited `side`, `strict`, `compress`,
  and `reserve` values.
- `Reset` and preset apply controls rebuild order-only custom rows with blank
  geometry.
- Saved config/session round-trips blank inherited geometry.
- Old configs matching the old web default `compress=true` shape are migrated
  to inherited geometry rows.

### Manual Browser Checks

After preparing the browser wheel:

1. Load EDL933.
2. Generate with simple `tuckin`; verify success.
3. Enable Custom Track Slots without editing; generate; verify success and
   geometry matching the simple preset.
4. Apply `middle` and `spreadout`; generate each; verify no unexpected
   auto-packed collapse.
5. Load the multi-record sample with depth data; enable Custom Track Slots and
   Show Depth; verify each record uses appropriate per-record geometry.
6. Explicitly request pure auto inside GC/skew tracks; verify the old
   ValidationError still appears for genuinely unfit layouts.

## Verification Commands

Run focused tests first:

```bash
python -m pytest tests/test_circular_track_slots.py tests/test_circular_radial_layout.py -q
python -m pytest tests/test_web_packaging.py -k "circular_track_slot or web_run_analysis or web_config" -q
```

Then run related circular/depth coverage:

```bash
python -m pytest tests/test_circular_feature_width.py tests/test_depth_track.py tests/test_circular_multi_canvas.py -q
```

If output SVG geometry changes intentionally, review and update reference SVGs:

```bash
python -m pytest tests/test_output_comparison.py -v -m circular
```

For the web wheel:

```bash
python tools/prepare_browser_wheel.py
python -m pytest tests/test_web_packaging.py -q
```

## Acceptance Criteria

- Custom Track Slots ON with no edits renders the same class of layout as the
  selected simple preset.
- Single-record EDL933 no longer fails only because Custom Track Slots were
  enabled.
- Multi-record inherited-preset custom layouts do not reuse fixed px geometry across
  records.
- Web state and saved sessions keep inherited preset geometry blank;
  preset-derived px geometry is never persisted as slot fields.
- Pure custom auto placement remains available and retains current validation
  behavior.
- Explicit user slot overrides are preserved and are not silently replaced by
  inherited preset defaults.

## Risks and Mitigations

- **Risk:** Inheriting preset geometry for blank custom slots could change
  existing CLI custom layouts.
  **Mitigation:** Preserve explicit pure-auto semantics and migrate only the
  known web default shape automatically. Document the new blank-geometry
  inheritance behavior clearly.

- **Risk:** Assigning preset lanes by order could surprise users with duplicate
  tracks.
  **Mitigation:** Use preset lanes only for the built-in/default lane sequence;
  unmatched duplicate tracks keep explicit or advanced custom behavior.

- **Risk:** Web saved sessions from the old simplified defaults may continue to
  fail.
  **Mitigation:** Add a migration or inference path for the exact old default
  slot shapes, and otherwise keep old custom sessions as explicit user layouts.

- **Risk:** UI users may not understand whether blank geometry means auto or
  preset-inherited.
  **Mitigation:** Show inherited placeholders such as "Inside (preset)" and
  make pure auto an explicit advanced choice.

## Implementation Order

1. Update web state/spec generation so default custom rows omit inherited
   geometry and placement fields.
2. Implement the preset-order overlay helper and unit tests for inherited
   geometry, order changes, and override preservation.
3. Integrate the helper into circular assembly for single-record and
   multi-record paths.
4. Convert the EDL933 regression into inherited-default success plus explicit
   pure-auto failure.
5. Update web state, reset/apply behavior, command generation, and config
   round-trip.
6. Update docs and CLI reference.
7. Prepare the browser wheel and run focused browser/manual checks.
