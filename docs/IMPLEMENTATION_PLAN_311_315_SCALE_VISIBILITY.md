# Implementation plan: coordinate-scale visibility

Date: 2026-07-30
Last revalidated: 2026-07-31 against commit `7aab08f8`

Issues:

- [#311: Allow "no ruler/bar"](https://github.com/satoshikawato/gbdraw/issues/311)
- [#315: Allow not showing size bar](https://github.com/satoshikawato/gbdraw/issues/315)

Status: planned; no implementation from this plan is present at the revalidated
baseline

## Current implementation checkpoint

The current code has scale style, dimensions, colors, intervals, and Linear
record-axis ruler placement, but it has no shared scale-visibility field, CLI
flag, or Web control.

Several existing boundaries already provide the required plumbing:

- Canonical dotted override paths and their types are derived from config
  dataclass annotations. Adding `ObjectsScaleConfig.show: bool` automatically
  makes `objects.scale.show` valid for both modes.
- The public Python APIs already accept shared `config_overrides`.
- Canonical request schema 5 generically encodes full configs and dotted
  overrides. The current Web session version is 40. Neither requires a version
  bump or a new legacy flat alias for this setting.
- Current session replay uses the authoritative canonical request. The
  `config.form` copy exists to restore Web controls and to support GUI/CLI
  conversion; it is not the render authority.

The renderer ownership is also narrower than the earlier plan assumed:

- Linear record-axis ruler drawing and its vertical reservation share
  `LinearRecordRenderContext.axis_ruler_enabled`. The bottom `LengthBarGroup` is
  constructed at one assembly site, and its measured height feeds final canvas
  sizing.
- Circular radial layout, legend bounds, and retry logic consume the resolved
  slot list. They do not need a second visibility Boolean. The work is to omit
  an implicit `ticks` slot before planning.
- Implicit Circular slots are materialized in several places: annotation
  preparation, CLI geometry shortcuts, multi-Depth API setup, and Conservation
  API setup. All of those paths must preserve the shared visibility intent.
- Web simple-mode geometry shortcuts also materialize a slot list before the
  Python renderer sees it.

This checkpoint determines the edit scope below. Files that already handle new
dataclass fields or dotted overrides generically are verification targets, not
expected production edit targets.

## Scope and product contract

Implement #311 and #315 as one change. Add one shared setting that controls the
primary genome-coordinate scale in both diagram modes:

```toml
[objects.scale]
show = true
```

The default remains `true`.

When `objects.scale.show` is `false`:

| Mode | Hidden | Retained |
| --- | --- | --- |
| Linear | Bottom scale bar or ruler; record-axis ruler ticks and coordinate labels | Each record's main axis line |
| Circular, implicit tracks | Primary coordinate tick marks, coordinate labels, and their radial reservation | The circular genome axis |
| Circular, explicit track slots | Nothing is inferred from the shared setting; enabled `ticks` slots remain authoritative | The slot-defined layout and circular axis |

GC-content and Depth axes and ticks are unrelated to this setting. Linear
definition text such as `Length / Coord.` also remains independent and continues
to use `objects.definition.linear.show_length`.

Keep `objects.scale.style` as the representation choice `bar | ruler`. Do not
add `none` to that field. Visibility and representation are separate concerns,
and Circular diagrams do not have bar and ruler representations. Validate the
style as an explicit two-value type while editing the scale model so that an
unknown value no longer falls through to a bar.

## Ownership and precedence

`ObjectsScaleConfig.show` owns simple, cross-surface visibility intent.
Renderers consume the resolved configuration before reserving layout space.
Hiding SVG elements after assembly is out of scope because it would retain empty
layout bands.

Circular track provenance matters because a non-`None` slot list is not always
user-authored:

1. With no user-explicit track list, `objects.scale.show` controls whether every
   internally generated or materialized list contains a primary `ticks` slot.
2. User-explicit Python `circular_track_slots`, CLI
   `--circular_track_order`/`--circular_track_slot`/
   `--circular_track_table`, and Web Custom Track Slots remain authoritative.
   Enabled `ticks` slots are drawn; absent or disabled `ticks` slots are not.
3. Annotation, multi-Depth, Conservation, and geometry-shortcut materialization
   remain part of simple/implicit behavior and must not turn the scale back on.
4. The Circular axis is independent of the `ticks` slot and is always retained.

The Web simple control must not compete with Custom Track Slots:

- Disable the simple checkbox while custom slots are active and point the user
  to the `Ticks` slot.
- Base Circular scale-styling control availability on effective ticks: the
  simple setting in simple mode, or enabled custom `ticks` slots in custom mode.
- Preserve the dormant simple value while custom mode is active. It becomes
  effective again when the app returns to simple mode.
- “Reset from simple controls” copies the simple visibility state and therefore
  omits ticks when the simple scale is hidden.
- “Reset to Tuckin/Middle/Spreadout” creates an explicit preset template and
  retains its `ticks` slot.

Expose the simple intent as `--hide_scale` on both CLI subcommands. In Circular
mode, user-explicit track options retain authority. If `--hide_scale` is
combined with a user-explicit list containing an enabled `ticks` slot, emit a
clear warning that the explicit slot wins. Do not warn for absent or disabled
ticks, or for an internally materialized list.

## Implementation work

### 1. Add and validate the configuration contract

Edit:

- `gbdraw/data/config.toml`
- `gbdraw/config/models/objects.py`
- `gbdraw/mode_profiles.py`

Changes:

- Append `show: bool = True` to `ObjectsScaleConfig` after its required fields,
  preserving direct-construction compatibility and valid dataclass field order.
- Read raw config values through the existing compatibility-oriented
  `_bool_from_config` normalizer. Canonical dotted overrides remain strict
  because the dataclass annotation is `bool`.
- Add `ScaleStyle = Literal["bar", "ruler"]` and a matching normalizer following
  the existing pairwise-style pattern.
- Reject unsupported scale styles with `ValidationError`.
- Keep `style = "bar"` and `show = true` as defaults.
- Apply the Linear ruler-axis managed color only when the scale is visible. A
  hidden, dormant `ruler_on_axis` value must not change an otherwise implicit
  axis color.

Verify without an expected production edit:

- `gbdraw/config/modify.py` derives the new path and both types from the model.
- `gbdraw/api/options.py` accepts the path as a shared override in both modes.

Do not make raw config Boolean parsing stricter as part of this feature. That
would reject inputs accepted by existing config readers and is a separate
compatibility decision.

### 2. Project the setting through CLI, API, and session boundaries

Edit:

- `gbdraw/circular.py`
- `gbdraw/linear.py`
- `gbdraw/api/diagram.py`
- `gbdraw/session_io.py`

Changes:

- Add `--hide_scale` to both CLI parsers and project it once as
  `objects.scale.show = false`.
- Keep Linear `--scale_style` limited to `bar | ruler`.
- Treat `--ruler_on_axis` as ineffective when the Linear scale is hidden. Reuse
  the existing ignored-option warning path, set the effective value to false,
  and choose the managed axis color from that effective state.
- In Circular CLI processing, retain provenance between user-explicit
  order/slot/table inputs and the list created by geometry shortcuts. Pass
  `show_ticks=not hide_scale` only to the internally generated list. Warn only
  when a user-explicit list has an enabled `ticks` slot.
- In `gbdraw/api/diagram.py`, pass the resolved visibility into both
  single-record and multi-record materializers:
  - `_default_circular_depth_slots_if_needed()` for multiple Depth tracks;
  - the default slot order created before inserting Conservation tracks.
- Keep the public Python surface small. Callers use
  `config_overrides={"objects.scale.show": False}`; no mode-specific option
  field is required.
- Add common `form.show_scale` projection to CLI-created session Web state.
- When converting GUI session state back to CLI arguments, add `--hide_scale`
  only when `form.show_scale is False`.
- For a CLI-created session containing both `--hide_scale` and
  `--ruler_on_axis`, project `form.linear_ruler_on_axis = false` because the CLI
  has already warned and ignored that option. Do not erase a dormant ruler
  choice from ordinary Web state merely because the user hides the scale.

Verify without an expected production edit:

- `gbdraw/session_request_codec.py` already encodes/decodes full typed configs
  and canonical dotted overrides generically.
- Current session render authority remains
  `renderRequest.diagramOptions.config.objects.scale.show`.
- Do not add `show_scale` to legacy flat request aliases; no released schema
  used that field.

### 3. Remove the Linear scale at the existing layout owners

Edit:

- `gbdraw/layout/linear.py`
- `gbdraw/diagrams/linear/assemble.py`

Changes:

- Include `cfg.objects.scale.show` in
  `LinearRecordRenderContext.axis_ruler_enabled`. The existing predicate then
  suppresses both record-axis ruler drawing and its vertical reservation.
- Add visibility to the sole production `LengthBarGroup` construction
  condition. If hidden, do not create the group, add it to the SVG, or include
  its height in final canvas and legend stacking.
- Keep each record's main axis path unchanged.
- Keep multi-record behavior consistent: there is no shared bottom scale today,
  and any record-local axis ruler is suppressed by the same resolved predicate.

No direct guard is expected in
`gbdraw/render/groups/linear/seq_record.py`; it already consumes
`axis_ruler_enabled`. No second hidden-state branch is expected in
`gbdraw/render/groups/linear/length_bar.py`; production construction is
centralized in the assembler, and scale-style validation makes its unknown-style
fallback unreachable.

### 4. Remove implicit Circular ticks before slot-based planning

Edit:

- `gbdraw/diagrams/circular/assemble.py`

The API and CLI materialization sites are covered in section 2.

Changes:

- Resolve simple-mode tick visibility from `cfg.objects.scale.show`.
- Pass it into annotation default-slot preparation instead of the current
  hard-coded `show_ticks=True`. Automatic annotations must not reinsert ticks.
- Use it for the implicit `show_ticks_track` decision. For a materialized
  user-explicit list, continue to derive visibility solely from enabled
  `ticks` renderers.
- Preserve one resolved slot list through radial planning, radial-label retry,
  and legend-bound calculation. Those consumers already derive the tick
  annulus from the slot list, so no new Boolean plumbing is needed.
- Do not condition Circular axis creation on scale visibility.

No production edit is expected in:

- `gbdraw/diagrams/circular/presets.py`, whose context already supports
  `show_ticks`;
- `gbdraw/tracks/circular.py`, whose default/order helpers already accept
  `show_ticks`;
- `gbdraw/diagrams/circular/radial_layout.py`, which intentionally resolves
  geometry from slots rather than the legacy Boolean argument.

### 5. Add the Web control and preserve simple/custom ownership

Edit:

- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/app/circular-track-slots.js`

Changes:

- Add `form.show_scale`, defaulting to `true`.
- Add a `Show Coordinate Scale` checkbox to both mode-specific
  `Axis & Scale` cards.
- In Linear mode, disable scale style, scale line/color/font controls,
  ruler-on-axis, interval, and ruler label color while hidden. Keep axis color
  and width enabled.
- In Circular simple mode, disable scale interval and tick-label font size while
  hidden. Keep axis color and width enabled.
- In Circular custom mode, disable only the simple checkbox and show a short
  explanation that `Ticks` slots own visibility. Keep scale-styling controls
  enabled exactly when the active custom stack has an effective enabled
  `ticks` slot.
- Include `form.show_scale` in `canUseLinearRulerOnAxis`.
- Gate the ruler control without mutating or zeroing the canonical
  `form.linear_ruler_on_axis` choice in ordinary Web state. The renderer
  predicate makes it ineffective while hidden, and showing the scale again
  restores the dormant choice.
- Add the semantic config key
  `showScale: "objects.scale.show"` to `CONFIG_OVERRIDE_PATHS`. This spelling
  deliberately projects back to `show_scale` through the existing semantic-key
  conversion.
- Emit `[CONFIG_OVERRIDE_PATHS.showScale]: form.show_scale !== false` for both
  modes. Do not use `Boolean(form.show_scale)`, because old and sparse state
  must default to visible.
- Project both full configs and canonical overrides back with
  `show_scale: overrides.show_scale !== false`.
- Include scale visibility in the Web managed Linear axis-color calculation so
  a hidden dormant ruler flag does not emit the ruler-axis color.
- Extend `createDefaultCircularTrackSlots()` with `showTicks = true`.
  Use the simple value for:
  - simple-mode geometry-shortcut materialization;
  - “Reset from simple controls”;
  - the default-template comparison used while reconciling Depth slots.
- Keep default dormant custom state and explicit preset-reset templates ticked.
  Enabling a previously stored custom stack with ticks may therefore show ticks
  even when the dormant simple value is false, as required by explicit
  precedence.
- Confirm ordinary config save/load, canonical session save/load, history
  restore, and global reset preserve or reset the value correctly.

No field-specific production edit is expected in
`gbdraw/web/js/services/config.js`,
`gbdraw/web/js/services/history-snapshot.js`, or
`gbdraw/web/js/services/reset.js`: they clone, merge, snapshot, or recreate the
form generically. Add focused lifecycle tests rather than parallel ownership.

Session version 40 and canonical request schema 5 remain unchanged. Missing
values in older state, full configs, and schema-5 requests resolve to `true`.

### 6. Update user documentation after behavior is implemented

Update the user-facing references that own the affected controls and concepts:

- `docs/CLI_Reference.md`
- `docs/FAQ.md`
- `docs/RECIPES.md`
- `docs/PYTHON_API.md`
- `docs/TUTORIALS/3_Advanced_Customization.md`
- `docs/TUTORIALS/5_Table_Driven_Inputs.md`
- `docs/TUTORIALS/6_Depth_Quantitative_Tracks.md`
- `docs/TUTORIALS/7_Linear_Layout.md`

Audit `docs/TUTORIALS/1_Customizing_Plots.md` for links or summaries affected by
the revised Linear scale workflow. No navigation change is currently expected
in `README.md`, `docs/DOCS.md`, `docs/INSTALL.md`, or `docs/GALLERY.md`; record
them as reviewed rather than editing them without a user-facing need.

Document:

- `--hide_scale` for both modes;
- the Python dotted override;
- the distinction between a coordinate scale and an axis;
- Circular simple-mode behavior and the explicit `ticks` slot override;
- Custom Track Slots control availability and reset semantics;
- the independence of GC-content, Depth, and Linear definition labels.

The following Gallery tutorial captures actually include the `Axis & Scale`
card and require visual review after the checkbox is added:

- `gbdraw/web/gallery/tutorials/lambda_basic_linear.json`
- `gbdraw/web/gallery/tutorials/hepatoplasmataceae_orthogroup.json`
- `gbdraw/web/gallery/tutorials/hepatoplasmataceae_collinear.json`
- `gbdraw/web/gallery/tutorials/BGC0000708-BGC0000713.json`

Regenerate only captures whose visible state or crop changes, following
`docs/skills/web-gallery-screenshot-maintenance/SKILL.md`. Create or update the
task-specific operation register and active capture plan before recapture. If
recapturing the three older tutorial definitions, bring their capture metadata
up to the current deterministic contract (`dataDependent`, exact session,
`assertAppState`, and `visibleControls`) in the same change.

Also audit without presumptive recapture:

- `HmmtDNA_ATskew.json`, whose custom `Ticks` screenshot demonstrates explicit
  authority but does not show the Axis card;
- `vibrio-harveyi-group-collinear.json` and
  `WSSV_genome_comparison.json`, whose settings tables mention scale controls
  although their captures do not show the card.

## Compatibility

- Default rendered SVGs and tracked reference SVGs remain byte-identical:
  `show = true`, `style = "bar"`.
- Existing config files, full-config requests, and sessions that omit
  `objects.scale.show` resolve it to `true`.
- Newly serialized full configs or sessions may contain the new `show` key;
  persisted JSON bytes are not expected to remain identical.
- Canonical request schema 5 and Web session version 40 remain sufficient.
- Existing user-explicit Circular track slots keep their current meaning.
- Scale colors, intervals, widths, and font sizes remain stored while the scale
  is hidden, so showing it again restores prior styling.
- Invalid custom scale styles that previously rendered as bars fail validation.
  This intentionally removes an undocumented fallback.
- Raw config Boolean normalization retains its existing compatibility behavior;
  typed dotted overrides reject non-Booleans.

## Test plan

### Configuration and Python APIs

- `tests/test_config_override_paths.py`
  - Default and missing `show` values resolve to `true`.
  - Boolean overrides are accepted for both modes; non-Booleans are rejected.
  - Raw config and canonical override paths reject invalid scale styles.
- `tests/test_api_requests.py`
  - Both diagram option types accept the shared dotted override.
- `tests/test_python_interface.py`
  - Root `draw_circular` and `draw_linear` calls carry the override through to
    rendering.
- `tests/test_mode_profiles.py`
  - A hidden Linear scale does not trigger ruler-axis managed color defaults.

### CLI and sessions

- `tests/test_linear_track_layout.py`
  - `--hide_scale --ruler_on_axis` warns, resolves the ruler option to false,
    and leaves the main axes.
- `tests/test_circular_track_slots.py`
  - Hidden plus user-explicit enabled ticks warns and retains ticks.
  - Hidden plus absent or disabled explicit ticks does not warn.
- `tests/test_circular_feature_width.py`
  - Geometry-shortcut materialization omits implicit ticks when hidden, while a
    user-specified order retains its ticks.
- `tests/test_session_io.py`
  - Both CLI subcommands store `form.show_scale = false` and the authoritative
    full config value.
  - GUI-to-CLI conversion restores `--hide_scale` in both modes.
  - A CLI session created from both hidden-scale and ruler-on-axis flags stores
    the ignored ruler option as false.
- `tests/test_session_request_codec.py`
  - Both modes round-trip the full config and canonical dotted override through
    generic schema-5 encoding.
  - No new legacy flat alias is introduced.

### Rendering and layout

- `tests/test_linear_track_layout.py`
  - Hidden bottom bar and ruler variants create no scale group and reserve no
    scale height.
  - Record-axis ticks and coordinate labels are hidden.
  - Main record axes remain.
  - Single- and multi-record cases are covered.
  - Default bar and ruler output remains unchanged.
- `tests/test_circular_track_slots.py`
  - An implicit hidden scale has no primary `ticks` group or annular
    reservation.
  - The Circular axis remains.
  - Explicit enabled/absent/disabled tick behavior follows slot authority.
  - GC-content and Depth ticks are unaffected.
- `tests/test_circular_annotation_tracks.py`
  - Automatic annotation slots do not restore implicit ticks.
- `tests/test_depth_track.py`
  - Single- and multi-record multiple-Depth auto slots do not restore ticks.
- `tests/test_circular_conservation.py`
  - Single- and multi-record Conservation insertion does not restore ticks.

### Web state, requests, and browser behavior

- `tests/web/session-request.test.mjs`
  - Both modes emit and project `objects.scale.show`.
  - Missing sparse state and missing schema-5 config values project as `true`.
  - Full Python configs with `show=false` project correctly.
  - Simple geometry shortcuts omit ticks; custom slots remain authoritative.
  - Hidden plus dormant Linear ruler keeps the ordinary managed axis color.
- `tests/web/mode-profiles.test.mjs`
  - New defaults, global reset, and cross-mode parity include `show_scale`.
- `tests/web/history.test.mjs`
  - Undo/redo restores scale visibility.
- `tests/test_web_packaging.py`
  - Static lifecycle wiring and Circular simple-reset/preset-reset distinctions
    are covered without duplicating behavior ownership.
- Add `tests/web/scale-visibility.playwright.spec.js`:
  - Both mode-specific checkboxes affect generated SVGs.
  - Linear `#length_bar` and Circular implicit `#tick` disappear when hidden.
  - Linear record axes and Circular `#Axis` remain.
  - Custom mode disables and explains the simple checkbox.
  - An enabled explicit ticks slot remains visible, and scale styling controls
    remain usable.

## Verification sequence

Run focused checks first:

```bash
pytest tests/test_config_override_paths.py tests/test_api_requests.py \
  tests/test_python_interface.py tests/test_mode_profiles.py -v
pytest tests/test_linear_track_layout.py -v
pytest tests/test_circular_track_slots.py \
  tests/test_circular_annotation_tracks.py \
  tests/test_circular_feature_width.py -v
pytest tests/test_depth_track.py tests/test_circular_conservation.py -v
pytest tests/test_session_io.py tests/test_session_request_codec.py -v
node tests/web/session-request.test.mjs
node tests/web/mode-profiles.test.mjs
node tests/web/history.test.mjs
pytest tests/test_web_packaging.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

Run `node --check` for each changed JavaScript module. Confirm both Node and
Python Playwright availability as described in `AGENTS.md`, then run the focused
browser spec. Prepare the generated browser wheel first if the browser or
packaging path requires it. Finish with the fast suite:

```bash
pytest tests/ -v -m "not slow"
```

Do not regenerate tracked SVG references unless a separately reviewed,
intentional default-geometry change is discovered.

## Acceptance criteria

- Linear diagrams can hide the bottom bar or ruler without leaving a blank
  vertical band.
- Linear record-axis ruler ticks and coordinate labels are hidden, while main
  record axes remain visible.
- Circular implicit layouts, including annotation, multiple-Depth,
  Conservation, and geometry-shortcut paths, omit primary coordinate ticks and
  their radial reservation.
- The Circular main axis remains visible.
- User-explicit Circular `ticks` slots remain authoritative across Python, CLI,
  and Web custom surfaces.
- Circular scale-styling controls follow effective custom ticks rather than a
  dormant simple value.
- GC-content and Depth axes and ticks do not change.
- CLI, Web, Python config overrides, and saved sessions converge on
  `objects.scale.show`.
- Old configs/sessions missing the field retain visible scales.
- Existing default output comparisons pass without reference regeneration.

## Delivery and coordination

Implement this plan as one PR for both issues. Keep it separate from #316.
Both changes touch `gbdraw/web/index.html`, `gbdraw/web/js/state.js`,
`gbdraw/web/js/services/session-request.js`, and Web/session tests. The branch
that lands second must rebase onto the first, resolve ownership-sensitive
conflicts manually, and rerun session-version-40/schema-5 verification. This
scale change does not consume a new session or canonical request version.
