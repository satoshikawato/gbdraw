# GC Content Percent Track Plan

Issue: https://github.com/satoshikawato/gbdraw/issues/211

## Goal

Add an optional GC content rendering mode that draws GC content as a depth-like absolute-value area track. The new mode should support a default 0 to 100 percent axis and user-defined intervals/ranges, without changing the current GC content output by default.

## Design Principles

- **SOLID**: Keep sequence analysis, value scaling, SVG path generation, and UI/CLI option handling separated. Do not make `DepthDrawer` responsible for GC semantics.
- **KISS**: Add one explicit GC content mode instead of introducing a generic plugin framework or broad track abstraction.
- **YAGNI**: Implement percent/range/tick support only. Defer multi-color thresholds, interval bands, and advanced legends until there is a concrete requirement.
- **DRY**: Reuse the existing depth-like area geometry by extracting common scalar-track path and axis helpers where practical, while preserving depth-specific labels and GC-specific labels.

## User-Facing Behavior

The existing GC content track remains the default:

- `deviation`: current behavior, centered around the record/window mean.

Add a new mode:

- `percent`: depth-like filled area from a baseline, scaled by absolute GC percentage.

Default percent mode axis:

- minimum: `0%`
- maximum: `100%`
- large tick interval: `20%`
- small tick interval: none

User-defined range examples:

- `0-100%`: full biological range.
- `30-70%`: emphasize a narrower expected GC range.
- `40-60%` with tick interval `5`: publication-focused view for closely related genomes.

Values outside the configured range should be clipped for drawing, but the original computed GC value should remain available in the data frame.

## Public Options

### Config

Extend `[objects.gc_content]` in `gbdraw/data/config.toml`:

```toml
[objects.gc_content]
mode = "deviation" # deviation | percent
min_percent = 0
max_percent = 100
show_axis = true
show_ticks = true
large_tick_interval = 20
small_tick_interval = "none"
tick_font_size = "auto"
```

The defaults must preserve existing SVG output unless `mode = "percent"` is selected.

### CLI

Add the same options to circular and linear modes:

```text
--gc_content_mode deviation|percent
--gc_content_min_percent FLOAT
--gc_content_max_percent FLOAT
--gc_content_tick_interval FLOAT
--gc_content_large_tick_interval FLOAT
--gc_content_small_tick_interval FLOAT
--gc_content_tick_font_size FLOAT
--show_gc_content_axis
--hide_gc_content_axis
--show_gc_content_ticks
--hide_gc_content_ticks
```

Validation:

- mode must be `deviation` or `percent`.
- min/max percent must be finite numbers.
- min percent must be less than or equal to max percent.
- tick intervals and font size must be positive when supplied.

### Web UI

Add a GC content display mode control in the existing GC/skew advanced area:

- segmented control or select: `Deviation`, `Percent`
- show min/max/tick controls only when `Percent` is active
- pass options through `gbdraw/web/js/app/run-analysis.js`

No new dependency or build step is required.

## Internal Design

### Configuration Model

Extend `ObjectsGcContentConfig` with:

- `mode: Literal["deviation", "percent"]`
- `min_percent: float | None`
- `max_percent: float | None`
- `show_axis: bool`
- `show_ticks: bool`
- `large_tick_interval: float | None`
- `small_tick_interval: float | None`
- `tick_font_size: float | None`

Use parsing helpers similar to the depth config helpers, but name validation errors with `gc_content_*` so users can identify the failing option.

### Configurator

Extend `GcContentConfigurator` rather than creating a separate GC-percent configurator. It already represents GC content drawing settings and is the right owner for mode-specific GC styling.

Add properties:

- `mode`
- `min_percent`
- `max_percent`
- `show_axis`
- `show_ticks`
- `large_tick_interval`
- `small_tick_interval`
- `tick_font_size`

Keep the existing color and opacity behavior.

### Data Preparation

Do not recalculate GC content. `analysis.skew.skew_df()` already provides the `"{nt} content"` column as a fraction from `0.0` to `1.0`.

Add a small helper, likely in `gbdraw/analysis/gc.py` or a new focused module, that converts a GC data frame to scalar plot data:

- `position`
- `value` or `percent`
- `value_normalized`

The helper should:

- read `f"{dinucleotide} content"`
- multiply by `100.0`
- clip only for normalized drawing
- handle empty data frames safely

### SVG Geometry

Extract common depth-like area path generation from:

- `gbdraw/svg/circular_tracks.py::generate_circular_depth_path_desc`
- `gbdraw/svg/linear_tracks.py::calculate_depth_path_desc`

Suggested common functions:

- `generate_circular_scalar_area_path_desc(...)`
- `calculate_linear_scalar_area_path_desc(...)`

Keep existing depth functions as wrappers for backward compatibility and readability.

### Axis Drawing

Do not reuse `depth_axis` IDs or `x` labels for GC. Add GC-specific axis labels:

- circular group id: `gc_content_axis`
- linear group id: `gc_content_axis`
- labels: `0%`, `20%`, ..., `100%`

To avoid duplication, extract a small tick-value helper shared by depth and GC content. Keep label formatting separate:

- depth formatter: `1.5x`
- GC formatter: `50%`

If extraction becomes too invasive, duplicate the minimal axis logic first and refactor after tests pass. The preferred final state is shared tick calculation, separate label formatting.

### Drawers and Groups

Update:

- `gbdraw/render/drawers/circular/gc_content.py`
- `gbdraw/render/drawers/linear/gc_content.py`

Behavior:

- `mode == "deviation"`: call the existing path code unchanged.
- `mode == "percent"`: call the scalar area path code and add the GC content axis when enabled.

Group IDs stay as `gc_content` so styling, legend behavior, and circular track slots continue to work.

### Circular Track Slots

Do not add a new renderer initially. Use the existing `dinucleotide_content` renderer and configure its behavior through `gc_content.mode`.

This keeps existing slot strings valid:

```text
gc_content:dinucleotide_content@nt=GC
```

If per-slot display modes become necessary later, add an optional slot parameter:

```text
gc_content:dinucleotide_content@nt=GC,mode=percent
```

That is out of scope for the first implementation.

## Implementation Steps

1. Extend config model parsing for `objects.gc_content`.
2. Add CLI/config override plumbing for GC content mode and percent-axis options.
3. Extend `GcContentConfigurator`.
4. Add GC percent scalar data helper.
5. Extract shared scalar area SVG path functions.
6. Add GC percent drawing path and axis rendering for circular mode.
7. Add GC percent drawing path and axis rendering for linear mode.
8. Add web UI state and argument generation.
9. Add tests.
10. Update user documentation and CLI reference.

## Test Plan

Unit tests:

- config defaults preserve `mode == "deviation"`.
- invalid GC content ranges and tick intervals raise validation errors.
- GC percent data helper maps `0.5` content to `50.0`.
- clipping affects normalized drawing values, not source percent values.

SVG tests:

- circular default GC output remains unchanged.
- linear default GC output remains unchanged.
- circular percent mode contains `id="gc_content_axis"` and percent labels.
- linear percent mode contains `id="gc_content_axis"` and percent labels.
- percent mode does not emit `id="depth_axis"` for GC content.

CLI tests:

- `--gc_content_mode percent` works in circular and linear commands.
- `--gc_content_min_percent` and `--gc_content_max_percent` validate ordering.

Web tests:

- saved/loaded config includes GC content mode fields.
- generated arguments include GC percent options only when selected.
- no offline packaging regression.

Reference SVGs:

- Do not refresh existing references unless default output changes unexpectedly.
- Add focused new reference or string-based SVG tests for percent mode to avoid unnecessary broad reference churn.

## Out Of Scope

- multi-color GC intervals
- threshold bands
- per-slot GC mode in circular track slots
- logarithmic GC scaling
- new legend semantics for percent ranges

These can be separate issues after the percent mode is stable.

## Risks And Mitigations

- **Risk**: Default output changes accidentally.
  **Mitigation**: Keep `deviation` as the default and add regression tests around existing output.

- **Risk**: Depth terminology leaks into GC SVG or UI.
  **Mitigation**: Use GC-specific IDs, labels, config names, and error messages.

- **Risk**: Circular layout changes because depth-like GC uses different baseline semantics.
  **Mitigation**: Reuse the existing GC content group and track width/radius. Only change the path inside the same reserved band.

- **Risk**: Axis labels collide in narrow circular tracks.
  **Mitigation**: Use automatic font sizing modeled after depth axis sizing and allow the axis to be hidden.

## Recommended First PR Scope

The first PR should implement only:

- config/API/CLI support
- circular and linear percent rendering
- percent axis labels
- focused tests

Web UI can be included in the same PR if the Python surface is stable, but it can also be a follow-up PR if review size becomes too large.
