# Implementation plan: coordinate-scale visibility

Date: 2026-07-30

Issues:

- [#311: Allow "no ruler/bar"](https://github.com/satoshikawato/gbdraw/issues/311)
- [#315: Allow not showing size bar](https://github.com/satoshikawato/gbdraw/issues/315)

Status: planned

## Scope and product contract

Implement #311 and #315 as one change. Add one shared setting that controls the
primary genome-coordinate scale in both diagram modes.

The new canonical setting is:

```toml
[objects.scale]
show = true
```

The default remains `true`.

When `objects.scale.show` is `false`:

| Mode | Hidden | Retained |
| --- | --- | --- |
| Linear | Bottom scale bar or ruler; record-axis ruler ticks and coordinate labels | Each record's main axis line |
| Circular, implicit tracks | Primary coordinate tick marks and coordinate labels | The circular genome axis |
| Circular, explicit track slots | Nothing is inferred from the shared setting; enabled `ticks` slots remain authoritative | The slot-defined layout and circular axis |

GC-content and Depth axes and ticks are unrelated to this setting. Linear
definition text such as `Length / Coord.` also remains independent and continues
to use `objects.definition.linear.show_length`.

Keep `objects.scale.style` as the representation choice `bar | ruler`. Do not add
`none` to that field. Visibility and representation are separate concerns, and
Circular diagrams do not have bar and ruler representations. Validate the style
as an explicit two-value type while editing the scale model so that an unknown
value no longer falls through to a bar.

## Ownership and precedence

`ObjectsScaleConfig.show` owns the simple, cross-surface visibility intent.
Renderers consume the resolved configuration before reserving layout space.
Hiding SVG elements after assembly is out of scope because it would retain empty
layout bands.

Circular explicit track slots remain the final authority:

1. If no explicit Circular track list is supplied, `objects.scale.show` controls
   whether the preset planner creates the primary `ticks` slot.
2. If an explicit track list is supplied, enabled `ticks` slots are drawn and
   absent or disabled `ticks` slots are not drawn.
3. The Circular axis is independent of the `ticks` slot and is always retained.

The Web simple control must not compete with Custom Track Slots. Disable the
simple checkbox while custom slots are active and point the user to the `Ticks`
slot. When the app returns to simple mode, the stored simple setting becomes
active again.

The CLI should expose the same simple intent as `--hide_scale` on both
subcommands. In Circular mode, explicit track options retain authority. If
`--hide_scale` is combined with an explicit list that contains an enabled
`ticks` slot, emit a clear warning that the explicit slot wins.

## Implementation work

### 1. Add and validate the configuration contract

Edit:

- `gbdraw/data/config.toml`
- `gbdraw/config/models/objects.py`
- `gbdraw/config/modify.py` only if its path validation requires an update
- `gbdraw/mode_profiles.py`

Changes:

- Add `show: bool = true` to `ObjectsScaleConfig`.
- Parse it with the existing strict Boolean configuration helper.
- Add a `ScaleStyle = Literal["bar", "ruler"]` normalizer.
- Reject unsupported scale styles with `ValidationError`.
- Keep `style = "bar"` and `show = true` as defaults.
- Update the Linear ruler-axis color override so it applies only when the scale
  is shown. A hidden ruler must not change an otherwise implicit axis color.

### 2. Project the setting through CLI and Python boundaries

Edit:

- `gbdraw/circular.py`
- `gbdraw/linear.py`
- `gbdraw/session_request_codec.py`
- `gbdraw/session_io.py`
- `gbdraw/api/options.py` only if shared override validation needs an explicit
  allowlist update

Changes:

- Add `--hide_scale` to both CLI parsers.
- Project the flag once as `objects.scale.show = false`.
- Keep Linear `--scale_style` limited to `bar | ruler`.
- Treat `--ruler_on_axis` as ineffective when `--hide_scale` is present. Reuse
  the existing ignored-option warning path and update its message.
- Preserve `show_scale` through CLI-created sessions and GUI/CLI replay.
- Keep the public Python surface small. Callers can use
  `config_overrides={"objects.scale.show": False}`; no new mode-specific option
  bundle is required.

### 3. Remove the Linear scale before layout measurement

Edit:

- `gbdraw/layout/linear.py`
- `gbdraw/diagrams/linear/assemble.py`
- `gbdraw/render/groups/linear/length_bar.py`
- `gbdraw/render/groups/linear/seq_record.py` only if a direct guard is needed

Changes:

- Include `cfg.objects.scale.show` in the resolved
  `axis_ruler_enabled` property.
- Do not construct `LengthBarGroup` when the scale is hidden.
- Do not reserve the group's height in the canvas or legend layout.
- Keep the record axis path unchanged.
- Make `LengthBarGroup` safe when constructed directly with a hidden scale. It
  must add no bar or ruler elements and report no drawable scale geometry.

### 4. Remove implicit Circular ticks before radial planning

Edit:

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/presets.py` if the preset context needs a named
  resolver
- `gbdraw/circular.py`

Changes:

- Resolve the implicit `show_ticks` value from `cfg.objects.scale.show`.
- Thread that value through preset creation, annotation-created default slots,
  radial layout, legend bounds, and any layout retry. Annotation preparation
  must not reinsert a `ticks` slot after the scale has been hidden.
- Keep explicit slot lists unchanged and derive their tick visibility from
  enabled `ticks` renderers.
- Ensure shortcut-generated implicit lists honor the setting. A
  user-specified `--circular_track_order`, `--circular_track_slot`, or
  `--circular_track_table` remains explicit.
- Do not condition Circular axis creation on scale visibility.

### 5. Add the Web control and preserve it through sessions

Edit:

- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/config.js` only where explicit field validation or
  restoration is required
- `gbdraw/web/js/app/circular-track-slots.js` if simple-mode reset currently
  recreates an unconditional `ticks` slot

Changes:

- Add `form.show_scale`, defaulting to `true`.
- Add a `Show Coordinate Scale` checkbox to both mode-specific
  `Axis & Scale` cards.
- Disable scale-only controls when the scale is hidden. Axis color and width
  controls remain enabled.
- Include `form.show_scale` in `canUseLinearRulerOnAxis`.
- In Circular Custom Track Slots mode, disable the simple checkbox and show a
  short explanation that `Ticks` slots control coordinate ticks.
- Add `objects.scale.show` to the canonical config-override paths.
- Emit it for both modes and project it back to `form.show_scale`.
- Confirm ordinary config save/load, canonical session save/load, history
  restore, and reset preserve the value.

No canonical request schema or Web session version change is expected. Missing
values in existing configuration and schema-5 requests resolve to `true`.

### 6. Update user documentation after behavior is implemented

Edit:

- `docs/CLI_Reference.md`
- `docs/FAQ.md`
- `docs/RECIPES.md`
- `docs/TUTORIALS/1_Customizing_Plots.md`
- `docs/TUTORIALS/5_Table_Driven_Inputs.md`
- `docs/TUTORIALS/7_Linear_Layout.md`

Document:

- `--hide_scale` for both modes.
- The distinction between a coordinate scale and an axis.
- Circular simple-mode behavior and the explicit `ticks` slot override.
- The independence of GC-content, Depth, and Linear definition labels.

Audit Gallery tutorials that capture the `Axis & Scale` cards. Regenerate only
the captures whose visible state or instructions change, following
`docs/skills/web-gallery-screenshot-maintenance/SKILL.md`.

## Compatibility

- Defaults remain byte-compatible: `show = true`, `style = "bar"`.
- Existing config files and sessions that omit `objects.scale.show` keep their
  current output.
- Canonical request schema 5 remains sufficient because config overrides
  already carry shared object settings.
- Explicit Circular track slots keep their current meaning.
- Scale colors, intervals, widths, and font sizes remain stored while the scale
  is hidden, so showing it again restores the prior styling.
- Invalid custom scale styles that previously rendered as bars will fail
  validation. This is an intentional correction of an undocumented fallback.
- Tracked reference SVGs should not change under default settings.

## Test plan

Python tests:

- `tests/test_config_override_paths.py`
  - Default `show` is `true`.
  - Boolean overrides are accepted for both modes.
  - Invalid style values fail.
- `tests/test_linear_track_layout.py`
  - `--hide_scale` omits the bottom scale group.
  - Hidden bar and ruler variants reserve no scale height.
  - Record axes remain.
  - `--hide_scale --ruler_on_axis` produces no axis ticks or coordinate labels.
  - Default bar and ruler output remains unchanged.
- `tests/test_circular_track_slots.py`
  - An implicit hidden scale has no primary `ticks` group.
  - The Circular axis remains.
  - Automatic annotation slots do not restore ticks.
  - An explicit enabled `ticks` slot overrides the shared hidden setting.
  - Explicit absent or disabled ticks still leave the axis.
  - GC-content and Depth ticks are unaffected.
- `tests/test_mode_profiles.py`
  - A hidden Linear scale does not trigger ruler-axis color defaults.
- `tests/test_session_io.py`
  - Both CLI subcommands preserve `--hide_scale` through session conversion.
- `tests/test_session_request_codec.py`
  - The shared Boolean path survives typed request encode/decode.

Web tests:

- `tests/web/session-request.test.mjs`
  - Both modes emit and project `objects.scale.show`.
  - Missing values project as `true`.
  - Custom Circular slots remain authoritative.
- `tests/web/linear-multi-record.playwright.spec.js`, or a smaller focused
  browser spec
  - Both checkboxes affect the generated SVG.
  - Axis elements remain present.
  - The Circular simple control is not presented as authoritative in custom
    slot mode.
- `tests/test_web_packaging.py`
  - Update static source assertions only where the setting lifecycle is
  intentionally covered.

Focused verification:

```bash
pytest tests/test_config_override_paths.py -v
pytest tests/test_linear_track_layout.py -v
pytest tests/test_circular_track_slots.py -v
pytest tests/test_mode_profiles.py -v
pytest tests/test_session_io.py -v
pytest tests/test_session_request_codec.py -v
node tests/web/session-request.test.mjs
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

Run a focused browser check after confirming the available Node or Python
Playwright installation.

## Acceptance criteria

- Linear diagrams can hide the bottom bar or ruler without leaving a blank
  vertical band.
- Linear record-axis ruler ticks and coordinate labels are also hidden.
- Circular implicit layouts can hide primary coordinate ticks and their radial
  reservation.
- Linear and Circular main axes remain visible.
- Explicit Circular `ticks` slots remain authoritative.
- GC-content and Depth axes and ticks do not change.
- CLI, Web, Python config overrides, and saved sessions converge on
  `objects.scale.show`.
- Existing default output comparisons pass without reference regeneration.

## Delivery and coordination

Implement this plan as one PR for both issues. It can be developed in parallel
with #316, but both changes touch `gbdraw/web/index.html`,
`gbdraw/web/js/state.js`, `gbdraw/web/js/services/session-request.js`, and Web
session tests. Merge this smaller scale change first, then rebase the #316 branch
before its final Web and session verification.
