# Circular Track Slots Unified Resolver Plan

## Goal

Make Circular Track Slots the single source of truth for circular radial layout.
The resolver must place `features`, `ticks`, `depth`,
`dinucleotide_content`, `dinucleotide_skew`, and `spacer` through one model,
without transitional ordered-layout helpers.

The target behavior is:

- Custom Track Slots are authoritative when supplied.
- Presets are only slot-list generators; they do not create a separate layout mode.
- Slot list order is honored consistently as axis-adjacent-to-outward order for
  outside slots and axis-adjacent-to-inward order for inside slots.
- The resolver uses a signed radial offset model internally where the circular
  axis is `0`.
- Public explicit circular geometry lives only on Track Slots. It keeps `r` as
  the existing center/anchor radius scalar and uses `w` for width; `ri` and
  `ro` are removed from all circular radial geometry grammars and data models.
- `side=inside|outside|overlay`, `strict`, `compress`, `reserve`, and
  `spacing` are normalized as slot-level fields with renderer capability
  checks instead of being hidden in renderer params.
- Every slot is normalized before resolver intent creation, so the resolver
  never sees a slot with an unspecified `side`.
- Conflicting renderer-specific side parameters are rejected instead of being
  resolved by precedence.
- TrackSpec is removed from the public and internal diagram configuration
  model. Circular geometry, visibility, side, draw order, and renderer params
  come from Track Slots, preset/default slot generation, or existing dedicated
  non-slot CLI/API/config options. Linear TrackSpec is also removed because it
  is currently parsed but not implemented.
- TrackSpec is not replaced by a visibility-only bridge or a renderer-param
  bridge. Any non-geometry control that is still needed must live in an
  existing dedicated option or become a first-class Track Slot field/renderer
  param.
- Existing circular CLI geometry options such as `--feature_width`,
  `--depth_width`, `--gc_content_radius`, `--gc_content_width`,
  `--gc_skew_radius`, and `--gc_skew_width` are translated into Track Slot
  overrides instead of using the removed TrackSpec bridge.
- Resolved layout entries use one `CircularResolvedSlot` structure for
  `features`, `ticks`, numeric tracks, depth tracks, spacers, and internal
  feature-label reservation slots.
- `CircularResolvedSlot` is a common geometry shell with an optional typed
  renderer payload. Resolved geometry must not be hidden in `params`.
- Resolved layout entries retain `slot_index`, so draw order is stable as
  `(z, slot_index)`.
- `_should_honor_custom_core_slot_order()`, `_resolve_ordered_custom_radial_layout()`, and `slot_layout.py` are removed.

This branch is allowed to make a one-shot breaking change to the experimental
circular Track Slots API and to retire the TrackSpec API even though TrackSpec
is currently wired through API options, CLI geometry overrides, legend and
visibility handling, and the circular resolver. Do not add staged migration
paths or compatibility layers. The one-shot change is acceptable only if every
current TrackSpec path is replaced or removed in the same implementation pass.
Backwards compatibility for the old public helper, `ri`/`ro`,
`CircularTrackPlacement`, `LinearTrackPlacement`, `TrackSpec`,
`parse_track_spec()`, `parse_track_specs()`, legacy kind fallback, or
`gap_after` names is not a goal. Backwards compatibility for old Web-saved
Custom Track Slots state is also not a goal: old saved state may be rejected on
import and discarded/reset on local startup instead of migrated. The
compatibility target is narrower: preset-generated `tuckin`, `middle`, and
`spreadout` layouts should preserve the fixed oracle commit's drawing geometry
unless a specific fix is documented.

## Current Problem

Custom Track Slots currently split across two layout paths:

- The default `resolve_circular_radial_layout()` path treats `features` and
  `ticks` as pre-resolved reference bands, then packs numeric tracks around
  them.
- The `honor_core_slot_order` path calls
  `_resolve_ordered_custom_radial_layout()`, which delegates to
  `slot_layout.resolve_circular_track_slots()` and then reconstructs
  feature/tick layouts.

That creates inconsistent behavior:

- Slot order is only honored for core slots under
  `_should_honor_custom_core_slot_order()`'s heuristic.
- Numeric slots and core slots do not share one resolver model.
- The ordered path does not share the full `side` semantics of the main radial
  resolver.
- The web UI can pass feature and tick slots with `r=1`, pinning both to the
  same radius and making collisions likely.
- Circular Track Slots and TrackSpecs both expose radial geometry, including
  `ri`/`ro` annulus syntax. That gives users two competing layout APIs and
  forces the resolver to reconcile incompatible geometry models.
- Resolved draw order cannot be specified cleanly because the current
  `CircularRadialLayout` separates `features`, `ticks`, and numeric tracks and
  does not keep the original slot index.

The issue is not just UI wording. The Python radial resolver must make Track
Slots authoritative and remove the old helper resolver.

## Files To Change

- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/builders.py`
- `gbdraw/diagrams/circular/presets.py`
- `gbdraw/diagrams/circular/slot_layout.py` (delete)
- `gbdraw/render/groups/circular/labels.py`
- `gbdraw/tracks/circular.py`
- `gbdraw/tracks/__init__.py`
- `gbdraw/tracks/parser.py` (delete after moving shared parsing helpers)
- `gbdraw/tracks/spec.py` (delete after moving `ScalarSpec`)
- `gbdraw/tracks/parsing.py` (new small shared parser helpers, if needed)
- `gbdraw/tracks/scalars.py` (new `ScalarSpec` home, if not kept elsewhere)
- `gbdraw/api/diagram.py`
- `gbdraw/api/options.py`
- `gbdraw/api/tracks.py`
- `gbdraw/api/__init__.py`
- `gbdraw/circular.py`
- `gbdraw/web/index.html`
- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/state.js`
- `tools/capture_circular_preset_oracle.py`
- `tests/test_circular_radial_layout.py`
- `tests/test_circular_track_slots.py`
- `tests/test_circular_feature_width.py`
- `tests/test_depth_track.py`
- `tests/test_circular_label_placement.py` where label reservation behavior is
  affected
- `tests/test_circular_preset_geometry.py` or equivalent preset-oracle tests
- `tests/fixtures/circular_preset_oracle/0228e6fe896768c6dc58513945d935e96578f5d0.json`
- `tests/test_output_comparison.py` or reference-output tests if preset SVG
  geometry changes unexpectedly
- `tests/test_web_packaging.py`
- `docs/CLI_Reference.md`
- `docs/RECIPES.md`
- `docs/TUTORIALS/1_Customizing_Plots.md`
- `docs/TUTORIALS/3_Advanced_Customization.md`
- `docs/CIRCULAR_TRACK_SLOTS_ROOT_FIX_PLAN.md` (retire, archive, or update so
  it no longer describes the deleted public helper)
- Any tests or docs that import `resolve_circular_track_slots()` from the old
  public helper
- Any tests or docs that import `TrackSpec`, `CircularTrackPlacement`,
  `LinearTrackPlacement`, `parse_track_spec()`, or `parse_track_specs()`
- Any tests or docs that use old shortcut slot strings such as `features@...`,
  `ticks@...`, `gc_content@...`, `gc_skew@...`, `depth@...`, or old annulus /
  spacing keys such as `ri`, `ro`, `gap`, and `gap_after`

## Current Tree Inventory

TrackSpec is not unused in the current tree. Removing it in one pass is still
the intended migration, but the implementation must cover these current
integration points:

- Public exports in `gbdraw/tracks/__init__.py`, `gbdraw/api/__init__.py`, and
  `gbdraw/api/tracks.py`.
- Definition and parsing in `gbdraw/tracks/spec.py` and
  `gbdraw/tracks/parser.py`.
- API entry points and bundled options in `gbdraw/api/diagram.py` and
  `gbdraw/api/options.py`, including `DiagramOptions` / `TrackOptions` and the
  linear ignored-TrackSpec warning path.
- Circular CLI compatibility code in `gbdraw/circular.py`, where
  `--feature_width`, `--depth_width`, `--gc_content_radius`,
  `--gc_content_width`, `--gc_skew_radius`, and `--gc_skew_width` are currently
  converted into internal TrackSpec strings.
- Circular assembly and radial layout in `gbdraw/diagrams/circular/assemble.py`,
  `gbdraw/diagrams/circular/radial_layout.py`, and the transitional
  `gbdraw/diagrams/circular/slot_layout.py`.
- Preset generation in `gbdraw/diagrams/circular/presets.py`, which currently
  uses `CircularTrackPlacement` and `ScalarSpec`.
- Web custom-slot state and serialization in `gbdraw/web/js/state.js`,
  `gbdraw/web/js/services/config.js`,
  `gbdraw/web/js/app/circular-track-slots.js`, and
  `gbdraw/web/js/app/run-analysis.js`.
- Focused tests in `tests/test_circular_track_slots.py`,
  `tests/test_circular_radial_layout.py`, `tests/test_circular_feature_width.py`,
  `tests/test_depth_track.py`, `tests/test_circular_label_placement.py`, and
  `tests/test_web_packaging.py`.
- Documentation and older plan documents that mention TrackSpec or old slot
  grammar.

## Public Semantics

### 1. Internal Signed Radial Coordinates

The unified resolver should use a signed radial offset coordinate internally:

- `0` is the circular axis/genome centerline.
- Positive offsets are outside the axis.
- Negative offsets are inside the axis.
- Absolute SVG radius is `axis_radius_px + offset_px`.

This is an internal layout model only. Public custom slot `r` keeps the current
meaning: a center/anchor radius in the existing `ScalarSpec` syntax. Therefore
`r=1` means one times the base circular radius, i.e. the circular axis radius.
The resolver converts public `r` to an internal offset with:

```python
anchor_offset_px = anchor_radius_px - axis_radius_px
```

This keeps the user-facing grammar compatible while still letting the resolver
compare all renderers in one signed coordinate space.

### 2. Slot Order

The slot list is interpreted as placement order within each side group.

For auto placement:

- `side=outside`: earlier slots are closer to the axis, and later slots move
  outward with larger positive offsets.
- `side=inside`: earlier slots are closer to the axis, and later slots move
  farther inward with more negative offsets.
- `side=overlay`: slots are measured at the axis radius unless an explicit `r`
  is supplied. They do not reserve space unless `reserve=true`.

Mixed-side slot lists still preserve each slot's relative intent. The packer
keeps separate inside/outside packing state, but processes unpinned auto slots in
original `slot_index` order so cross-side generated reservations, such as
external feature-label footprints, affect later slots deterministically. Final
resolved entries must retain the original `slot_index`.

This deliberately avoids reverse-packing outside tracks from an unknown outer
boundary. Outside auto slots grow outward from the occupied axis-adjacent band.

Pinned slots and auto slots may be mixed. A pinned slot is a slot with explicit
`r`; an auto slot has no explicit `r` and is placed by the resolver. The
resolver may measure pinned slots first so auto slots can avoid them, but final
same-side radial order must still match slot list order. If a pinned slot forces
the final resolved order to contradict the slot list, raise `ValidationError`.
This order validation applies to pinned-pinned and pinned-auto combinations and
uses each slot's resolved `packing_band`, not the full `reserved_band`.

### 3. Slot Normalization And Side Defaults

All circular slots must pass through a normalization step before resolver intent
creation. Normalization is the only place where missing `side` values are
filled in. Resolver code may assume every enabled slot has a normalized
`side`.

- `features`: `side` and `lane_direction` are normalized together.
  - If neither is supplied, use `lane_direction=inside` and `side=inside`.
  - If only `lane_direction=inside` is supplied, derive `side=inside`.
  - If only `lane_direction=outside` is supplied, derive `side=outside`.
  - If only `lane_direction=split` is supplied, derive `side=overlay` and
    `reserve=true`.
  - If only `side=inside` or `side=outside` is supplied, derive the matching
    `lane_direction`.
  - If only `side=overlay` is supplied, derive `lane_direction=split` and
    `reserve=true`.
  - If both are supplied and conflict, validation raises instead of choosing
    one. `lane_direction=split` is compatible only with `side=overlay`.
- `ticks`: `side` controls the packing channel. `label_side` and `tick_side`
  control the measured tick footprint around the chosen anchor radius. These are
  not aliases for `side`, so they are validated for known values and then
  measured directly. If `side` is not supplied, use `side=inside`.
  Preset-generated tick slots should emit explicit `side`, `r`, `w`,
  `label_side`, `tick_side`, and `preset` values needed to reproduce
  fixed-oracle preset geometry.
  Tick measurement produces both a `packing_band` and a `reserved_band`.
  `packing_band` is used for same-side order and spacing, while
  `reserved_band` includes tick labels and is registered as a global occupied
  band. Therefore `side=inside,label_side=outside` still participates in inside
  order, but outside slots must move outward to avoid the outside label
  footprint.
- Numeric and depth slots use `side` as their packing channel. If `side` is not
  supplied, use `side=inside`.
- `spacer` slots use `side` as their packing channel. If `side` is not
  supplied, use `side=inside`.
- `overlay` slots do not reserve space unless `reserve=true`.
- `side=overlay,reserve=false` means the slot is drawn at the axis anchor but
  does not affect packing. This is allowed but should not be emphasized in the
  web UI.
- `side=overlay,reserve=true` registers the measured reserved band as a global
  occupied band. Inside and outside auto slots must avoid that band. Other
  reserving overlay slots do not auto-shift away from it; overlap follows the
  same warning/strict error policy as pinned reserved-band overlap.
  This is the normalized form for `features` with `lane_direction=split` and is
  required to preserve `middle` preset geometry.

External feature labels are not public custom Track Slots in this migration.
They are generated internally after feature-slot measurement as
`feature_labels_outer` and/or `feature_labels_inner` reservation slots when
external labels are enabled.

- The generated label slots are not emitted by the CLI or web UI and are not
  accepted in user-supplied slot specs.
- They participate in radial packing only through their measured
  `reserved_band`; leader-line and text drawing remain fixed label draw phases.
- Their `reserved_band` is registered as global occupied space when the owning
  feature slot is measured, so auto analysis tracks move around labels instead
  of being silently moved to another side.
- Generated label reservations are blocker intervals owned by the measured
  feature slot. If a generated label reservation collides with an already
  measured auto slot, the monotone packer shifts that auto slot and later
  same-side auto slots away from the axis, then remeasures the shifted slots.
  It should raise `ValidationError` only when no placement exists under pinned
  geometry, side channels, same-side order, spacing, radial bounds, and
  compression rules.
- Pinned user slots that overlap a generated label reservation warn by default
  and raise with `strict=true`, following the same pinned-overlap policy as other
  reserved bands.
- If no valid auto placement remains, raise `ValidationError` with a message
  that points users to move the analysis slot, disable external labels, or supply
  explicit geometry.

Renderer-specific side parameters must not silently contradict the generic
`side` parameter. Invalid values should raise during parsing or normalization,
not later during drawing.

Renderer capability validation must also happen during parsing or
normalization:

- `compress` is valid only for `side=inside` auto `dinucleotide_content`,
  `dinucleotide_skew`, and `depth` slots.
- `compress` supplied on `features`, `ticks`, or `spacer` is rejected.
- `compress` supplied on a pinned numeric/depth slot is rejected.
- `compress` supplied on outside or overlay numeric/depth slots is rejected.

Public input slots should keep layout-affecting options as direct fields, not
inside `params`:

```python
@dataclass(frozen=True)
class CircularTrackSlot:
    id: str
    renderer: str
    enabled: bool = True
    side: str | None = None
    radius: ScalarSpec | None = None
    width: ScalarSpec | None = None
    spacing: ScalarSpec | None = None
    z: int = 0
    strict: bool | None = None
    compress: bool | None = None
    reserve: bool | None = None
    params: Mapping[str, Any] = field(default_factory=dict)
```

After normalization, resolver code consumes a separate internal slot shape with
all generic layout defaults resolved:

```python
@dataclass(frozen=True)
class NormalizedCircularTrackSlot:
    slot_index: int
    id: str
    renderer: str
    enabled: bool
    side: str
    radius: ScalarSpec | None
    width: ScalarSpec | None
    spacing: ScalarSpec | None
    z: int
    strict: bool
    compress: bool
    reserve: bool
    params: Mapping[str, Any]
```

`params` is renderer-specific only. Examples include `lane_direction` for
`features`, `label_side`/`tick_side`/`preset` for `ticks`, and `nt` or
`legend_label` for numeric tracks. Generic layout keys such as `side`,
`strict`, `compress`, `reserve`, `spacing`, `radius`, `width`, and `z` should
not be stored in `params` after parsing or normalization.

### 4. Explicit Geometry

Circular Track Slots are the only public radial geometry API for circular
tracks. They should use only:

- `r`: center/anchor radius using the existing radius scalar syntax.
- `w`: requested draw width.

The CLI string grammar keeps only the explicit full form:

```text
<slot_id>:<renderer>@key=value,key=value
```

Shortcut slot strings such as `features@...`, `ticks@...`, and
`gc_content@...` are removed and must raise parse errors. This keeps slot
identity and renderer selection explicit.

`ri` and `ro` should be removed from circular custom slot grammar, tests, and
public circular placement data models. Annulus syntax is not a good fit for
`features` and `ticks`, because feature lanes and tick labels can reserve a
footprint larger than a simple draw annulus.

Implementation direction:

- Remove `inner_radius` and `outer_radius` from circular radial data models.
- Replace `CircularTrackPlacement` and `LinearTrackPlacement` with direct
  `CircularTrackSlot` fields for circular inputs: `radius`, `width`,
  `spacing`, `side`, `z`, `strict`, `compress`, and `reserve`.
- Remove `TrackSpec`, `CircularTrackPlacement`, `LinearTrackPlacement`,
  `TrackKind`, `LayoutMode`, `parse_track_spec()`, and `parse_track_specs()`
  from the public API and internal diagram configuration path.
- Keep `ScalarSpec` as a generic scalar helper for Track Slots, either by
  moving it out of `tracks/spec.py` or by otherwise decoupling it from the
  removed TrackSpec API.
- Remove `track_specs=` from circular and linear Python API entry points,
  `DiagramOptions`, and any CLI plumbing. Dedicated CLI/API/config options
  remain the way to control non-slot layers.

Explicit behavior:

- `r` pins the slot center/anchor radius.
- `w` pins the requested draw width.
- `r+w` pins the slot anchor and requested draw width without moving them during
  auto packing. The reserved footprint may still be larger for renderers such
  as `features` and `ticks`.
- Explicit pinned slots are measured before auto slots and added to occupied
  bands.
- Pinned `reserved_band` overlap logs a warning by default.
- Pinned `reserved_band` overlap with `strict=true` raises.
- Final same-side resolved radial order must match slot list order. If explicit
  pins make that impossible, raise `ValidationError`.

Spacing behavior:

- Replace the old `gap_after` field and `gap`/`gap_after` grammar with
  `spacing`.
- Public custom slot grammar uses `spacing=<ScalarSpec>`.
- The public and normalized in-memory slot field is `spacing`.
- Resolver intents carry `spacing_px`.
- `spacing` means the minimum gap between `packing_band`s to the next slot in
  final same-side radial order. For a same-side adjacent pair `(A, B)` sorted by
  `slot_index`, the required gap is `A.spacing_px`.
- `spacing` is an auto-packing constraint, not resolved draw geometry. It does
  not need to be stored on `CircularResolvedSlot`, and it must not be stored in
  resolved-slot `params`.
- The default auto spacing is `max(1.0, 0.01 * axis_radius_px)` unless a
  normalized slot supplies an explicit `spacing`.
- Explicit `spacing=0` is allowed.
- Pinned-pinned slots whose `packing_band`s violate the same-side gap log a
  warning by default and raise `ValidationError` if either slot has
  `strict=true`.
- Pinned-auto pairs move the auto slot to satisfy the same-side gap.
- Auto-pinned pairs that cannot satisfy both order and gap raise
  `ValidationError`.
- `side=overlay,reserve=true` does not use same-side spacing. Its reserved band
  is expanded by its own `spacing_px` and registered as a global occupied band.

Remove the web defaults that emit feature and tick slots with `r=1`. That
value currently pins both slots to the axis radius. Preset-generated slots
should emit explicit slot geometry in the existing `r+w` form when they need
stable fixed-oracle placement.

### 5. Presets And Custom Slots

Preset mode is only an input expansion step.

- `circular_track_slots_for_preset()` builds the slot list for `tuckin`,
  `middle`, and `spreadout`.
- The resolver does not branch on preset names after slots are created.
- Preset-generated slots are normal slots once they reach the resolver.
- Preset-generated slots should include the explicit `r`, `w`, `side`,
  `spacing`, and renderer params needed to preserve the fixed-oracle
  drawing geometry for `tuckin`, `middle`, and `spreadout`.
- If depth, GC content, or GC skew is enabled by configuration, the preset
  generator adds those slots in the same fixed-oracle geometry order and spacing
  unless the change is intentional and documented.
- If the user supplies Custom Track Slots, those slots are authoritative.

This means depth or duplicate numeric tracks added later are simply additional
auto slots. The resolver handles them by side, order, and available space.

Preset compatibility should be locked by tests derived from the current `main`
commit, not by guessing values during implementation. The oracle base commit is
fixed as:

```text
0228e6fe896768c6dc58513945d935e96578f5d0
```

Add a small oracle capture workflow that records only preset resolved geometry
from that commit and stores those expected values in versioned test fixtures.
Test execution must compare against the stored fixtures and must not check out
or import a moving `main` branch.

Keep the oracle intentionally narrow:

- Capture `tests/test_inputs/MG1655.gbk` as the baseline representative input.
- Capture one additional stress input with materially different length/label
  behavior, selected from `GCF_015097735.1_ASM1509773v1_genomic.gbff`,
  `AP027280.gb`, or `MjeNMV.gb`. Add more inputs only when an existing
  regression path proves they cover distinct geometry behavior.
- Cover `tuckin`, `middle`, and `spreadout`, `strandedness=true|false`, and a
  small visibility matrix: features+ticks only, GC+skew, and depth+GC+skew.
- Store center/anchor radii, widths, packing bands, reserved bands, tick sides,
  label sides, slot ids, and slot order.
- Do not store full SVG output, text glyph positions, or font-measured label
  bounding boxes in the preset oracle. Those belong in focused rendering or
  label-reservation tests.

The new resolver must match those geometry values unless an intentional visual
fix is documented. Store the oracle as JSON, for example under
`tests/fixtures/circular_preset_oracle/0228e6fe896768c6dc58513945d935e96578f5d0.json`,
and keep the capture script in `tools/` for regeneration. Numeric tolerances
must be explicit in tests: use a tight tolerance such as `abs=1e-6` for pure
layout arithmetic, and avoid font-dependent values in this oracle.

This oracle is an implementation guard, not a complete visual-regression gate.
After the resolver rewrite, generated SVGs and rendered diagrams must still be
reviewed. Reference SVG updates are made only after visual inspection confirms
that any differences are intentional or acceptable.

### 6. TrackSpec Removal And CLI Compatibility

TrackSpec is removed from the API surface and internal diagram configuration
path. It should not be treated as a visibility bridge, renderer-param bridge, or
geometry fallback. Linear TrackSpec is not preserved because it is currently
parsed but ignored.

Do not add a replacement compatibility bridge for visibility-only TrackSpecs or
renderer-param TrackSpecs. TrackSpec's former responsibilities are split into
the APIs that actually own them:

- radial geometry, side, spacing, compression, reservation, and draw order live
  on circular Track Slots;
- renderer-specific circular slot semantics live in Track Slot `params`;
- non-slot visibility and non-slot layers remain controlled by existing config,
  CLI, and API options;
- if removing TrackSpec exposes a genuinely missing non-geometry control, add a
  focused first-class option in the owning subsystem instead of routing it
  through TrackSpec.

- CLI `--track_spec` input is removed or raises a migration-oriented error that
  points users to circular Track Slots or existing dedicated CLI/config options.
- Circular and linear Python API entry points no longer accept `track_specs=`.
  Passing that obsolete keyword to the new function signature may raise Python's
  normal `TypeError`; this is acceptable for removed constructor/function
  arguments.
- `DiagramOptions` and `TrackOptions` no longer expose `track_specs`.
- `resolve_circular_radial_layout()` does not accept or merge TrackSpecs.
- TrackSpec exact-id matching and legacy kind fallback matching are both removed
  from the circular resolver.
- Non-slot layers such as labels, definition text, and legend remain controlled
  by their existing config and dedicated CLI/API options, not by TrackSpec.
- Keep `ScalarSpec` as a generic scalar helper for circular Track Slots; do not
  keep the TrackSpec parser or placement models just to host `ScalarSpec`.

Existing top-level circular geometry CLI options remain supported by translating
them directly into Track Slot overrides before resolver entry:

- `--feature_width` -> `features.width`
- `--depth_width` -> `depth.width`
- `--gc_content_radius` -> `gc_content.radius`
- `--gc_content_width` -> `gc_content.width`
- `--gc_skew_radius` -> `gc_skew.radius`
- `--gc_skew_width` -> `gc_skew.width`

The conversion happens after preset/default slot generation and after
`--circular_track_order` expansion. If the user supplies explicit
`--circular_track_slot` entries, these legacy geometry CLI options are rejected
instead of being merged. The error should tell users to put `r=` or `w=` on the
matching circular track slot. This keeps Custom Track Slots authoritative while
preserving the existing public CLI geometry options in non-custom mode.

### 7. Public API Boundary

The public circular Track Slot API should expose inputs only. Do not replace the
old lower-level helper with another public resolver helper.

Keep public:

- `ScalarSpec`
- `CircularTrackSlot`
- `CircularTrackSlotParseError`
- `parse_circular_track_slot()`
- `parse_circular_track_slots()`
- `default_circular_track_slots()`
- preset slot-list helpers only if they are useful as input generators

Remove from public exports:

- `TrackSpec`
- `TrackKind`
- `LayoutMode`
- `LinearTrackPlacement`
- `CircularTrackPlacement`
- `CircularTrackLayoutContext`
- `CircularSlotFootprint`
- `ResolvedCircularTrackSlot`
- `resolve_circular_track_slots()`
- `TrackSpecParseError`
- `parse_track_spec()`
- `parse_track_specs()`

`CircularResolvedSlot`, `CircularSlotPayload`, `CircularFeatureLabelLayout`,
`_RadialSlotIntent`, and `CircularRadialLayout` are internal resolver
structures. They may live in `radial_layout.py`, but should not be exported
through `gbdraw.tracks` or `gbdraw.api.tracks`.

### 8. Validation And Error Policy

This migration intentionally does not preserve old circular Track Slot
compatibility, but failures must be consistent.

- Circular Track Slot string parsing errors raise `CircularTrackSlotParseError`.
- Invalid `CircularTrackSlot` objects supplied through Python APIs raise
  `ValidationError` after the object has been constructed with the new field
  names.
- Obsolete constructor or function keyword arguments, such as
  `CircularTrackSlot(..., placement=...)`,
  `CircularTrackSlot(..., gap_after=...)`, or `track_specs=...`, are not
  compatibility-handled. Python's normal `TypeError` is acceptable because those
  inputs do not construct a new API object.
- Resolver code assumes normalized valid slots and should not silently repair
  invalid layout values.
- Unknown `side`, `lane_direction`, `tick_side`, or `label_side` values raise
  during parsing or normalization.
- Shortcut slot strings such as `features@...`, `ticks@...`, and
  `gc_content@...` raise during parsing. Use the full
  `<slot_id>:<renderer>@...` form.
- Conflicts such as `features:features@side=outside,lane_direction=inside`
  raise during normalization.
- Unsupported `compress` use, such as
  `features:features@compress=true`, `ticks:ticks@compress=true`, or
  `gc_content:dinucleotide_content@r=0.6,compress=true`, or
  `gc_content:dinucleotide_content@side=outside,compress=true`, raises during
  parsing or normalization.
- Unsupported legacy geometry or spacing keys (`ri`, `ro`, `inner_radius`,
  `outer_radius`, `gap`, `gap_after`) raise with messages that point users to
  `r+w` or `spacing`.
- Any remaining CLI TrackSpec input raises with a message that points users to
  circular Track Slots or dedicated circular CLI/API/config options. This
  includes visibility-only TrackSpecs and non-geometry renderer params, not only
  radial geometry.

## Resolver Design

### 1. Remove The Transitional Resolver Paths

Delete from `assemble.py`:

- `_should_honor_custom_core_slot_order()`
- the `honor_core_slot_order=(...)` callsite

Delete from `radial_layout.py`:

- `_ordered_slot_layout_context()`
- `_resolve_ordered_custom_radial_layout()`
- the `honor_core_slot_order` parameter on `resolve_circular_radial_layout()`
- the branch that calls `_resolve_ordered_custom_radial_layout()`

Delete `gbdraw/diagrams/circular/slot_layout.py` and remove public exports and
tests that target its old API. Do not add a replacement public resolver API;
resolved circular layout is an internal assembly detail.

### 2. Normalize Slots Before Intent Creation

Add a resolver-local or tracks-level normalization step that returns slots with
all layout-affecting defaults resolved.

Normalization must:

- preserve original order by assigning `slot_index` before filtering disabled or
  otherwise non-rendered slots;
- normalize renderer aliases;
- reject shortcut slot strings and require `<slot_id>:<renderer>@...`;
- reject `ri`/`ro` and all annulus aliases in circular Track Slots;
- parse only `spacing`, not `gap` or `gap_after`;
- derive or validate `side` for every enabled slot using the rules in Public
  Semantics;
- keep generic layout data on direct fields and renderer-specific data in
  `params`;
- reject unsupported renderer capability combinations, including `compress` on
  `features`, `ticks`, `spacer`, pinned numeric/depth slots, or outside/overlay
  numeric/depth slots;
- reject any obsolete TrackSpec path before intent creation, if any callsite
  still attempts to pass it internally during the migration;
- reject conflicting side/lane values before measuring;
- keep preset-generated slots as ordinary normalized slots.

### 3. Use Unified Slot Intents

Replace numeric-only `_SlotIntent` with `_RadialSlotIntent`:

- `slot`
- `slot_index`
- `slot_id`
- `renderer`
- `side`
- `anchor_offset_px`
- `width_px`
- `explicit_anchor`
- `explicit_width`
- `spacing_px`
- `params`

The intent list must preserve supplied slot order. Disabled slots are skipped,
but remaining slots keep their original `slot_index`.

### 4. Measure Any Slot At A Candidate Offset

Add one resolver-local measurement function:

```python
def _measure_radial_slot(
    intent: _RadialSlotIntent,
    *,
    anchor_offset_px: float,
    axis_radius_px: float,
    feature_dict,
    canvas_config,
    cfg,
    total_length: int,
    tick_track_channel_override: str | None,
) -> CircularResolvedSlot:
    ...
```

It should compute signed draw/reserved bands for packing comparisons, convert
the final geometry back to absolute SVG radii, and return a unified
`CircularResolvedSlot`. The resolved slot is a common radial geometry shell plus
an optional typed payload for renderers that need richer resolved state:

```python
CircularSlotPayload = (
    CircularFeatureLayout
    | CircularTickLayout
    | CircularFeatureLabelLayout
)

@dataclass(frozen=True)
class CircularResolvedSlot:
    slot_index: int
    id: str
    renderer: str
    side: str
    z: int

    anchor_radius_px: float | None
    anchor_offset_px: float | None
    requested_width_px: float | None
    resolved_width_px: float | None

    packing_band_px: RadialBand | None
    draw_band_px: RadialBand | None
    reserved_band_px: RadialBand | None

    params: Mapping[str, Any]
    payload: CircularSlotPayload | None = None

    explicit_anchor: bool = False
    explicit_width: bool = False
    compressed: bool = False
```

The resolved slot fields represent:

- `slot_index`
- `id`
- `renderer`
- `side`
- `z`
- `anchor_offset_px`
- `anchor_radius_px`
- requested width, if supplied or defaulted before compression
- resolved width after any allowed compression
- packing band in absolute radii, used for same-side order and spacing
- draw band in absolute radii
- reserved band in absolute radii, used as global occupied space
- `params`, containing renderer-specific semantic options only
- optional typed payload for renderers that need richer resolved state
- `explicit_anchor`
- `explicit_width`
- `compressed`

Resolved geometry must not be stored in `params`. Numeric, depth, and spacer
slots use the common geometry shell directly and should have `payload=None`
unless a future renderer needs additional typed resolved state. `features` uses
`CircularFeatureLayout`, `ticks` uses `CircularTickLayout`, and internal
feature-label reservation slots use `CircularFeatureLabelLayout`.

Renderer behavior:

- `features`: build `CircularFeatureLayout` at
  `axis_radius_px + anchor_offset_px`; `packing_band` and `reserved_band` are
  `all_band_px`; store the layout in `payload`.
- `ticks`: build `CircularTickLayout` at
  `axis_radius_px + anchor_offset_px`; `packing_band` is the tick mark band and
  `reserved_band` is the union of tick marks and tick labels; store the layout
  in `payload`.
- numeric: `packing_band`, draw band, and reserved band are the center-width
  band; `payload` is `None`.
- `spacer`: `packing_band` and `reserved_band` are the reserved band only; no
  drawing. `draw_band_px` is `None` and `payload` is `None`.
- internal feature-label reservation: `draw_band_px` is `None`,
  `packing_band_px` is `None`, and `reserved_band_px` is the measured label
  footprint. Store the label layout in `payload` so label draw phases do not
  recompute placement.

This keeps `features` and `ticks` out of special top-level layout fields while
still letting their renderers consume typed resolved state.

### 5. Pack Auto Slots With A Deterministic Monotone Interval Packer

The packer is a one-dimensional monotone interval packer over signed radial
offsets. It is not a general-purpose constraint solver. Auto slots may move only
away from the axis on their side, while pinned slots never move. The packer
places slots against blocker intervals and, when a later-discovered blocker such
as a feature-label reservation collides with an already placed auto slot, it
rewinds to the first affected auto slot on that side and replays placement from
there.

This keeps the algorithm deterministic and avoids bidirectional or cyclic
constraint solving. A layout that can be satisfied by pushing auto slots away
from the axis should succeed. A layout that would require moving an auto slot
toward the axis, moving a pinned slot, violating same-side order, or exceeding
inside bounds fails with `ValidationError`.

The packer keeps two measured bands for each slot:

- signed `packing_band` for same-side order and `spacing_px`;
- signed `reserved_band` for global collision checks.

Hard boundaries:

- Outside slots have no hard outer bound. They may grow outward and update
  `outer_content_radius_px`.
- Inside slots are finite. They must keep absolute radii above `0` and must not
  overlap any global occupied band, including the center definition band when
  present.
- Overlay auto slots start at the axis and do not auto-shift because overlay has
  no inward/outward packing direction. They do not reserve space unless
  `reserve=true`.

Constraint model:

- Pinned slots are fixed variables. Explicit `r` pins the anchor, and explicit
  `w` pins the requested width. Pinned slots may overlap by warning policy, or
  raise when `strict=true`.
- Auto slots are movable variables. Their anchor may be shifted only away from
  the axis in their side's direction.
- Same-side order constraints are hard. For outside slots, later slot
  `packing_band`s must be farther outward than earlier same-side slots. For
  inside slots, later slot `packing_band`s must be farther inward than earlier
  same-side slots.
- Same-side spacing constraints use the earlier slot's `spacing_px`.
- Global exclusion constraints use `reserved_band`s. Auto slots must avoid
  pinned slots, placed/reserving overlays, tick label footprints, internal
  feature-label reservations, the center definition band, and other auto slots'
  reserved bands.
- Feature-label reservations are dependent blocker intervals owned by the
  resolved feature slot. If the owning feature slot moves or is remeasured, its
  old internal label reservation blockers are replaced with reservations
  measured from the same label layout data.
- Tick label footprints are part of the tick slot's `reserved_band` and are
  included every time a tick slot is measured.

Packer flow:

1. Normalize slots and create intents in original slot order.
2. Measure pinned slots and reserving overlay slots at their fixed anchors.
3. Initialize every auto slot at the nearest axis-adjacent anchor allowed by
   same-side order, spacing, and pinned constraints.
4. Place auto slots in `slot_index` order. For each slot, choose the nearest
   valid anchor in its side's outward direction that satisfies same-side order,
   spacing, radial bounds, and all currently known blocker intervals.
5. Measure the slot at that anchor. Tick slots add their full tick-label
   footprint to their own `reserved_band`. Feature slots create or replace their
   owned internal feature-label reservation blockers immediately after
   measurement.
6. If a newly added or replaced blocker collides with already placed auto slots,
   find the earliest colliding auto slot on each affected side, discard measured
   auto slots from that point outward on that side, and replay those slots with
   the expanded blocker set.
7. Replaying one slot may force later same-side auto slots to move farther away
   from the axis to preserve same-side order and spacing.
8. Repeat replay only when a newly measured slot adds or changes blocker
   intervals that affect already placed auto slots.
9. If the inside channel cannot satisfy all constraints at the requested width,
   try allowed compression for compressible auto numeric/depth slots and replay
   with the deterministic compressed-width sequence.
10. If no placement exists within finite inside bounds, compression limits,
    pinned constraints, and same-side order, raise `ValidationError` with the
    blocking slot ids and the suggested remedy.
11. Merge measured slots and sort draw dispatch by `(z, slot_index)`.

The packer must use a deterministic tie-break:

- Prefer the smallest movement away from the axis.
- When two movements are equivalent, prefer preserving the previous anchor.
- When compression is allowed, try requested/default width first, then shrink in
  a deterministic sequence down to the renderer's documented minimum width.
- Use a convergence guard. Reaching the guard is a bug or an under-specified
  replay loop and should raise `ValidationError` with enough diagnostic state to
  debug the involved slots.

Same-side ordering rules:

- Outside slots must remain axis-adjacent-to-outward by `slot_index`.
- Inside slots must remain axis-adjacent-to-inward by `slot_index`.
- Same-side order validation uses final `packing_band`s, not full
  `reserved_band`s.
- Same-side radial order contradictions always raise `ValidationError`,
  regardless of `strict`, because the requested slot order cannot be honored.
- `spacing_px` is the minimum gap from each slot's `packing_band` to the next
  same-side slot's `packing_band` in final radial order. This applies to
  pinned-pinned, pinned-auto, and auto-pinned adjacent pairs.

Global collision rules:

- Global collision checks use `reserved_band`.
- Global occupied space includes pinned slots, placed auto slots, tick label
  footprints, `side=overlay,reserve=true` overlays, internal feature-label
  reservations, and the center definition band.
- A reserved overlay slot with explicit `spacing` expands its global occupied
  band by `spacing_px`.
- Pinned user slots that collide with global occupied space warn by default and
  raise with `strict=true`.
- Inside and outside auto slots must not overlap global occupied space. If no
  valid placement exists, raise `ValidationError`. Reserving overlay auto slots
  do not shift; overlap follows the overlay collision policy above.

Compression:

- Only numeric/depth auto slots can compress.
- `compress=true` allows auto numeric/depth tracks to shrink only after
  normal-width placement fails in the finite inside channel.
- `features`, `ticks`, `spacer`, and pinned numeric/depth slots do not compress;
  specifying `compress=true` for them raises.
- Auto slots must not overlap occupied/reserved bands.
- Compression must be deterministic: try the requested/default width first, then
  shrink down to the renderer's documented minimum width. If the minimum width
  still does not fit, raise `ValidationError`.
- Outside auto slots do not need compression to satisfy radial bounds because the
  outside channel can grow outward; `compress=true` is rejected there to keep the
  semantics narrow.

### 6. Use A Unified Radial Layout Output

Replace the split `features` / `ticks` / `tracks` output with one layout tuple:

```python
@dataclass(frozen=True)
class CircularRadialLayout:
    axis: CircularAxisLayout
    slots: tuple[CircularResolvedSlot, ...]
    definition_reserved_band_px: RadialBand | None
    outer_content_radius_px: float
```

Helper functions may retrieve feature and tick layouts by scanning
`layout.slots`, but those are convenience lookups only. The resolver source of
truth is the unified slot tuple.

The final draw dispatcher should sort all resolved drawable slots by:

```python
key = (resolved_slot.z, resolved_slot.slot_index)
```

Placement order and draw order must remain separate concepts. Radial placement
uses side/order semantics; rendering overlap uses `(z, slot_index)`.

### 7. Draw Phases

Only resolved drawable slots participate in `(z, slot_index)` sorting. Internal
feature-label reservation slots are resolved slots, but they are non-drawable and
exist to carry label occupancy and typed label payloads. Label leaders/text,
definition text, and legend keep fixed diagram phases so their visual layering
remains predictable.

Draw in this order:

1. Axis.
2. External label leader-line underlay, when external labels are shown.
3. Resolved drawable slots sorted by `(z, slot_index)`.
4. External label text overlay, when external labels are shown.
5. Center definition, when enabled.
6. Legend, when enabled.

Label calculation and label drawing should be separated. `prepare_label_list()`
and `assign_leader_start_points()` should run once in assembly, and the measured
external label footprint should be stored in internal label reservation slots.
Rendering should draw external label leaders and external label text in separate
phases from the resolved label payload instead of recomputing placement.
Embedded labels remain part of the `features` slot drawing, so they are drawn
with features and never behind the feature shapes. If the current external label
renderer cannot draw leader lines and text separately, split it as part of this
change. The legend should be drawn last so no slot can overwrite it.

## Web UI Changes

Update `gbdraw/web/index.html`:

- Remove `ri` and `ro` inputs from the Custom Track Slots editor.
- Replace the `gap` input with `spacing`.
- Bind `side`, `strict`, `compress`, and `reserve` as slot-level fields, not
  `slot.params.*`.

Update `gbdraw/web/js/app/circular-track-slots.js`:

- Stop creating default custom slots with `features.radius = '1'`.
- Stop creating default custom slots with `ticks.radius = '1'`.
- Stop setting `ticks.width = '0.02'` unless the user explicitly enters a
  width.
- Rename UI and serialization support from `gapAfter` / `gap_after` to
  `spacing`.
- Emit only the full `<slot_id>:<renderer>@...` string form. Do not emit or
  accept shortcut slot strings.
- Continue preserving semantic params:
  - `features.params.lane_direction`
  - `ticks.params.label_side`
  - `ticks.params.tick_side`
  - `ticks.params.preset`
- Move generic layout settings out of `params` in web state and serialization:
  `side`, `strict`, `compress`, `reserve`, `spacing`, `radius`, `width`, and
  `z` are slot-level fields.
- Update validation and serialization to remove `innerRadius` and
  `outerRadius` for custom circular slots.
- If the UI exposes `r`, label it as center/anchor radius using the existing
  factor/px syntax. Do not describe it as a signed offset.
- Ensure every serialized slot has slot-level placement fields or enough
  renderer-specific semantic params for the Python normalizer to derive `side`.
- Emit the current Custom Track Slots schema version with saved/exported web
  configuration. Do not serialize compatibility aliases for older schemas.

Update `gbdraw/web/js/app/run-analysis.js`:

- When Custom Track Slots are enabled, do not send `--feature_width`,
  `--depth_width`, `--gc_content_radius`, `--gc_content_width`,
  `--gc_skew_radius`, or `--gc_skew_width`.
- The web UI should either disable those simple geometry controls in custom
  mode or raise a client-side validation error before running analysis. It
  should not rely on the CLI to reject the invalid combination.

Update `gbdraw/web/js/services/config.js`:

- Accept only the current Custom Track Slots schema version.
- Reject imported configuration that contains old Custom Track Slot shapes, such
  as `gapAfter`, `gap_after`, `innerRadius`, `inner_radius`, `outerRadius`,
  `outer_radius`, `params.side`, `params.radius`, `params.width`, or
  `placement`.
- When stale Custom Track Slots state is found in local browser storage during
  app startup, discard that saved custom-slot state and reset to the new
  defaults instead of migrating it.
- Surface a clear validation error for explicit config imports so users know the
  file uses a retired experimental Custom Track Slots schema.

Update `gbdraw/web/js/state.js`:

- Store default custom slots with the new slot-level shape, the current schema
  version, and no default `r` or `w` for features/ticks.

Update `tests/test_web_packaging.py`:

- Assert that web serialization emits `spacing`, does not emit `gap`, `ri`, or
  `ro`, and suppresses legacy geometry CLI options when Custom Track Slots are
  enabled.
- Assert that old saved Custom Track Slots shapes are rejected on explicit
  import and reset/discarded from local storage startup state instead of being
  migrated.

## Implementation Steps

1. Add failing tests that define the new behavior.
   - Custom slot order `ticks,features` places ticks outside features without
     `honor_core_slot_order`.
   - Custom slot order `features,ticks` places ticks inside features for
     tuckin-like defaults.
   - Custom order `features,gc_content,ticks,gc_skew` produces visual radial
     order matching the slot list.
   - Reordering `gc_skew,gc_content` swaps numeric radii.
   - Repeated `side=outside` slots respect axis-adjacent-to-outward order.
   - Conflicting feature side settings, such as
     `features:features@side=outside,lane_direction=inside`, are rejected.
   - `lane_direction=split` derives `side=overlay,reserve=true`.
   - `side=overlay,reserve=true` reserves a global occupied band avoided by
     inside and outside auto slots, while other reserving overlays use the
     pinned-overlap warning/strict policy.
   - Slots with omitted `side` are normalized before resolver intent creation.
   - Explicit `r+w` pins a slot and does not get moved.
   - Pinned overlap warns by default and raises with `strict=true`.
   - Pinned-pinned and pinned-auto same-side radial order contradictions raise
     `ValidationError`.
   - Auto-pinned combinations that cannot satisfy both same-side order and gap
     raise `ValidationError`.
   - `spacing` creates the minimum gap to the next same-side slot and replaces
     `gap_after`.
   - Shortcut slot strings such as `features@...`, `ticks@...`, and
     `gc_content@...` are rejected.
   - Unsupported `compress`, including `features:features@compress=true`,
     `ticks:ticks@compress=true`, `spacer:spacer@compress=true`, and pinned
     numeric/depth compression or outside numeric/depth compression, is
     rejected.
   - Resolved slots draw in `(z, slot_index)` order.
   - TrackSpec API inputs and parser entry points are removed or rejected.
   - Existing circular CLI geometry options are translated into slot-level
     overrides in non-custom mode.
   - Existing circular CLI geometry options are rejected when explicit
     `--circular_track_slot` entries are supplied.
   - `ri`/`ro` and `gap`/`gap_after` are rejected with migration-oriented error
     messages.
   - Invalid new-shape `CircularTrackSlot` objects raise `ValidationError`,
     invalid string specs raise `CircularTrackSlotParseError`, and obsolete
     constructor/function kwargs are allowed to fail with Python `TypeError`.
   - Tick label reserved bands displace slots on the opposite side when labels
     cross the tick slot's packing side.
   - External feature labels generate internal non-drawable reservation slots
     whose `reserved_band`s displace auto analysis tracks.
   - External feature-label reservations can displace auto analysis tracks even
     when the owning feature slot appears after those auto tracks in the input
     slot list.
   - Pinned slots that overlap internal label reservations warn by default and
     raise with `strict=true`.
   - External label leaders draw below slots, external label text draws above
     slots, and embedded labels remain inside the features slot rendering.

2. Add compact fixed-main preset geometry oracle tests.
   - Capture geometry from oracle base commit
     `0228e6fe896768c6dc58513945d935e96578f5d0`.
   - Store the captured fixture as JSON, for example under
     `tests/fixtures/circular_preset_oracle/0228e6fe896768c6dc58513945d935e96578f5d0.json`.
   - Keep a regeneration script such as
     `tools/capture_circular_preset_oracle.py`.
   - Do not check out or import `main` during normal test execution.
   - Capture `MG1655.gbk` plus one additional stress input selected from
     `GCF_015097735.1_ASM1509773v1_genomic.gbff`, `AP027280.gb`, or
     `MjeNMV.gb`.
   - Add more oracle inputs only when they cover a distinct preset geometry
     behavior not represented by the first two inputs.
   - Cover `tuckin`, `middle`, and `spreadout`, both strandedness modes, and the
     small visibility matrix: features+ticks only, GC+skew, and depth+GC+skew.
   - Assert feature/tick/numeric centers, widths, packing bands, reserved bands,
     sides, and order in the new resolver.
   - Do not include full SVG output, text glyph positions, or font-measured label
     bounding boxes in this oracle. Cover label reservations with focused tests.
   - Use explicit tight arithmetic tolerances for pure layout values.
   - Treat the oracle as an implementation guard only. Review generated SVGs or
     rendered diagrams after implementation before accepting reference-output
     changes.

3. Remove old circular geometry grammar.
   - Keep `r` and `w`.
   - Keep only the full `<slot_id>:<renderer>@...` slot string grammar.
   - Remove shortcut slot strings such as `features@...`, `ticks@...`, and
     `gc_content@...`.
   - Remove `ri` and `ro` from circular slot parsing, data models, and tests.
   - Remove `TrackSpec`, `CircularTrackPlacement`, `LinearTrackPlacement`,
     `TrackKind`, `LayoutMode`, `parse_track_spec()`, and
     `parse_track_specs()` from public exports and internal diagram paths.
   - Remove `track_specs=` from circular and linear API entry points and
     `DiagramOptions` / `TrackOptions`.
   - Rename `gap_after`/`gap` to `spacing` in parsing, tests, web state, and
     web serialization.
   - Update web serialization to stop emitting `innerRadius` and `outerRadius`.
   - Update docs or API errors to state that annulus geometry is unsupported
     for circular slots and that TrackSpec has been removed.

4. Simplify public circular Track Slot API exports.
   - Delete `gbdraw/diagrams/circular/slot_layout.py`.
   - Remove `resolve_circular_track_slots()`, `CircularTrackLayoutContext`,
     `CircularSlotFootprint`, and `ResolvedCircularTrackSlot` exports from API
     modules.
   - Remove TrackSpec exports and replace circular slot inputs with direct
     `CircularTrackSlot.side`, `.radius`, `.width`, `.spacing`, `.z`,
     `.strict`, `.compress`, and `.reserve` fields.
   - Keep `ScalarSpec` public as a generic scalar helper outside the removed
     TrackSpec API.
   - Keep only input-side Track Slot parsing/default helpers public.
   - Delete or rewrite tests that target `slot_layout.py` internals.

5. Add slot normalization.
   - Assign `slot_index` before filtering.
   - Normalize all enabled slots to a concrete `side`.
   - Validate conflicting `side` / renderer-specific placement params.
   - Keep generic layout data on direct fields and renderer-specific data in
     `params`.
   - Convert public `spacing` to `spacing_px` for intent creation.
   - Reject unsupported `compress` combinations before measuring.
   - Reject any obsolete internal TrackSpec path before intent creation.
   - Raise `CircularTrackSlotParseError` for invalid strings and
     `ValidationError` for invalid new-shape Python API slot objects.

6. Introduce signed-offset intent and measurement helpers in
   `radial_layout.py`.
   - Resolve public slot `r` to `anchor_radius_px`, then convert to
     `anchor_offset_px = anchor_radius_px - axis_radius_px`.
   - Convert measured signed bands back to absolute radius bands for existing
     renderers.
   - Keep existing feature and tick measurement helpers where reusable.
   - Return `CircularResolvedSlot` for every renderer.
   - Keep resolved geometry on common `CircularResolvedSlot` fields and use typed
     payloads only for richer renderer state such as features, ticks, and
     feature labels.

7. Rewrite `resolve_circular_radial_layout()` to build all enabled slots
   through the unified resolver.
   - Measure pinned slots.
   - Pack auto slots with the deterministic signed monotone interval packer.
   - Auto slots move only away from the axis on their side; pinned slots never
     move.
   - Use `packing_band` for same-side order/spacing and `reserved_band` for
     global occupied-space collision checks.
   - Validate final same-side radial order after pinned and auto slots are
     measured.
   - Generate internal feature-label reservation blockers from feature
     measurement. If those blockers collide with already placed auto slots,
     rewind to the earliest affected auto slot on each affected side and replay
     that side outward.
   - Build final `CircularRadialLayout(slots=...)` from measured slots.

8. Update preset expansion and circular CLI geometry compatibility.
   - Make `circular_track_slots_for_preset()` emit any explicit geometry in
     existing `r+w` terms plus `side` and `spacing`.
   - Confirm `tuckin`, `middle`, and `spreadout` retain fixed-oracle
     feature/tick/numeric radii and widths using the oracle tests.
   - Do not add resolver preset branches.
   - Translate `--feature_width`, `--depth_width`, `--gc_content_radius`,
     `--gc_content_width`, `--gc_skew_radius`, and `--gc_skew_width` into
     slot-level overrides after default/preset or `--circular_track_order`
     expansion.
   - Reject those legacy geometry CLI options when explicit
     `--circular_track_slot` entries are supplied.
   - Remove or reject any remaining CLI `--track_spec` input instead of
     converting it.

9. Update label rendering and `assemble.py`.
   - Remove `_should_honor_custom_core_slot_order()`.
   - Remove the `honor_core_slot_order` callsite.
   - Dispatch resolved drawable slots through a single `(z, slot_index)` sort.
   - Exclude internal feature-label reservation slots from drawable-slot sorting.
   - Draw non-slot layers in fixed phases: axis, external label leader-line
     underlay, slots, external label text overlay, definition, legend.
   - Split external label leader-line drawing from external label text drawing.
   - Draw external labels from the resolved internal label payload without
     recomputing placement.
   - Keep embedded labels in the `features` slot rendering path.
   - Draw legend last.
   - Keep label and legend calculations using resolved reserved bands.

10. Update web defaults.
   - Default custom slots should leave `r` and `w` blank unless user-entered.
   - New depth/GC/skew slots should be auto slots.
   - Update `index.html` to remove `ri`/`ro`, replace `gap` with `spacing`,
     and bind generic layout controls to slot-level fields.
   - Update `circular-track-slots.js` to serialize full slot specs only and
     emit `spacing`, not `gapAfter`, `gap_after`, `ri`, or `ro`.
   - Update `run-analysis.js` so Custom Track Slots mode does not send legacy
     circular geometry CLI options.
   - Update `services/config.js` to accept only the current Custom Track Slots
     schema, reject old imported custom-slot state, and reset/discard stale
     local-storage custom-slot state on startup.
   - Update `state.js` to use the new default slot shape and schema version.

11. Run focused tests.
   - `python -m pytest tests/test_circular_radial_layout.py tests/test_circular_track_slots.py -q`
   - `python -m pytest tests/test_circular_preset_geometry.py -q`
   - `python -m pytest tests/test_web_packaging.py -k "circular_track_slot or web_run_analysis" -q`

12. Run broader circular tests if focused tests pass.
    - `python -m pytest tests/ -m circular -q`

## Acceptance Criteria

- No `honor_core_slot_order` parameter remains.
- No `_should_honor_custom_core_slot_order()` remains.
- No `_resolve_ordered_custom_radial_layout()` remains.
- `gbdraw/diagrams/circular/slot_layout.py` is removed.
- Circular slot geometry uses `r` and `w`; `ri` and `ro` are not part of the
  custom slot grammar or circular radial data models.
- Circular slot strings require the full `<slot_id>:<renderer>@...` form;
  shortcut strings such as `features@...`, `ticks@...`, and `gc_content@...`
  are rejected.
- TrackSpec public API and parser entry points are removed:
  `TrackSpec`, `CircularTrackPlacement`, `LinearTrackPlacement`, `TrackKind`,
  `LayoutMode`, `parse_track_spec()`, and `parse_track_specs()` are no longer
  exported.
- Circular and linear public API entry points no longer accept `track_specs=`;
  obsolete kwargs may fail with Python `TypeError`.
- Any remaining CLI `--track_spec` input is removed or rejected with a
  migration-oriented error.
- Existing circular CLI geometry options are translated into slot-level
  overrides in non-custom mode and rejected when explicit custom Track Slots are
  supplied.
- Circular slot spacing uses `spacing`; `gap` and `gap_after` are not part of
  the custom slot grammar.
- Resolver intents use `spacing_px`; no `spacing_after_px` field remains.
- `spacing` is the minimum gap to the next slot in final same-side radial
  order.
- Public `r` keeps the current center/anchor radius meaning; the resolver
  internally places slots in signed radial offsets with axis `0`.
- Every slot reaching resolver intent creation has a normalized `side`.
- Generic layout settings live on direct slot fields; `params` contains only
  renderer-specific settings.
- Unsupported `compress` combinations are rejected. Only `side=inside` auto
  `dinucleotide_content`, `dinucleotide_skew`, and `depth` slots may use
  `compress=true`.
- `side=overlay,reserve=true` registers a global occupied band avoided by inside
  and outside auto slots; other reserving overlays do not auto-shift and use the
  pinned-overlap warning/strict policy.
- Resolved slots expose `packing_band` for same-side order/spacing and
  `reserved_band` for global occupied-space collision checks.
- Auto slots are packed by a deterministic signed monotone interval packer:
  auto slots move only away from the axis on their side, newly discovered
  blockers may cause replay from the earliest affected auto slot on each
  affected side, pinned slots never move, outside can grow outward, and inside
  is bounded by radius `0` and occupied bands.
- Tick label footprints are included in `reserved_band`; slots on the opposite
  side are displaced when tick labels cross the tick slot's packing side.
- External feature labels generate internal non-drawable reservation slots whose
  reserved bands displace auto tracks and collide with pinned tracks according to
  the normal strictness policy.
- External feature-label reservations can displace earlier auto tracks when the
  slot list puts the owning feature slot after those auto tracks.
- Custom Track Slots order affects `features`, `ticks`, and numeric tracks
  consistently.
- Pinned-pinned and pinned-auto final same-side radial order contradictions
  raise `ValidationError`.
- Auto-pinned combinations that cannot satisfy final same-side order and gap
  raise `ValidationError`.
- Presets are implemented as slot-list generation only.
- Preset-generated `tuckin`, `middle`, and `spreadout` layouts match
  the fixed oracle commit's feature/tick/numeric geometry unless a change is
  intentional and documented.
- Preset compatibility is enforced by compact fixed-main oracle tests for
  `MG1655.gbk` plus one additional stress input selected from
  `GCF_015097735.1_ASM1509773v1_genomic.gbff`, `AP027280.gb`, or `MjeNMV.gb`;
  additional inputs are added only when they cover distinct geometry behavior.
- The oracle fixture records base commit
  `0228e6fe896768c6dc58513945d935e96578f5d0`; normal tests do not check out a
  moving `main` branch.
- Circular slot geometry is the only circular radial geometry source; TrackSpec
  is removed instead of merged or ignored.
- No TrackSpec merging, `show=false` bridging, visibility-only TrackSpec bridge,
  renderer-param TrackSpec bridge, exact-id matching, or legacy kind fallback
  remains in the circular resolver.
- `CircularRadialLayout` stores one `slots: tuple[CircularResolvedSlot, ...]`
  collection instead of separate `features`, `ticks`, and `tracks` fields.
- `CircularResolvedSlot` is a common geometry shell with optional typed payloads
  for feature, tick, and feature-label resolved state.
- Resolved geometry is stored on `CircularResolvedSlot` fields, not in `params`;
  `params` remains renderer-specific semantic input.
- Resolved slots keep `slot_index`.
- Slot draw dispatch uses `(z, slot_index)`, while axis, label leader lines,
  label text, definition, and legend use fixed draw phases.
- External label leader lines draw under slots, external label text draws over
  slots from the resolved label payload, embedded labels remain with the features
  slot, and legend draws last.
- Web-generated default custom slots do not pin `features` and `ticks` to the
  same radius.
- Web Custom Track Slots UI no longer exposes `ri`, `ro`, or `gap`, serializes
  `spacing`, stores generic placement fields at slot level, emits only full
  slot specs, and suppresses legacy circular geometry CLI options in custom
  mode.
- Web Custom Track Slots saved/exported state uses the current schema version
  only. Old imported custom-slot state is rejected, and stale local-storage
  custom-slot state is reset/discarded on startup instead of migrated.
- Public circular Track Slot API exports input helpers plus `ScalarSpec` and
  `CircularTrackSlotParseError` only; no public low-level circular slot resolver
  replacement is exported.
- Existing preset geometry remains stable against the fixed oracle commit except
  for intentional, documented fixes.
- Preset geometry oracle tests are implementation guards. Final SVG/reference
  output changes still require visual inspection before acceptance.

## Risks

- This is a one-shot migration of experimental circular Track Slots behavior.
  External compatibility is intentionally not preserved on this private branch.
- Removing TrackSpec, TrackSpec parser helpers, and the public lower-level
  resolver will break experimental callers that used those APIs directly.
- Removing shortcut slot strings will break callers that used `features@...`,
  `ticks@...`, or kind-like slot specs instead of the full
  `<slot_id>:<renderer>@...` form.
- Rejecting old Web-saved Custom Track Slots state means users may need to
  recreate experimental custom-slot layouts after this branch lands. Local
  storage reset must be scoped to custom-slot state so unrelated UI preferences
  are not discarded unnecessarily.
- Translating legacy circular geometry CLI options into slot overrides must
  happen after preset/default or `--circular_track_order` slot generation, or
  the options may silently apply to the wrong layout source.
- Rejecting legacy circular geometry CLI options with explicit
  `--circular_track_slot` is intentionally stricter than previous behavior and
  may require users to move `r`/`w` values into their slot specs.
- Rejecting unsupported `compress` combinations is intentionally stricter than
  ignoring them, but avoids a false impression that features, ticks, spacers, or
  pinned/outside numeric/depth slots can shrink.
- `side=overlay,reserve=true` creates a global occupied band, so middle-preset
  compatibility depends on including the split feature footprint and its
  dependent reservations in the blocker set before finalizing auto anchors.
- Changing outside auto-slot order from outer-to-inner to
  axis-adjacent-to-outward may change diagrams that depended on the old wording.
- The deterministic monotone interval packer may need replay passes because
  later feature-label reservations can force earlier auto slots to move away
  from the axis. It must have a replay guard and clear diagnostics for any
  unexpected non-converging replay loop.
- Pinned-auto combinations are allowed, but order validation must compare final
  same-side `packing_band` order while collision validation compares global
  `reserved_band` occupancy. Mixing these bands up can silently invert the
  requested slot order or leave labels overlapping.
- Tick labels can reserve space on the opposite side of their slot's packing
  side, so outside or inside tracks may move even when no same-side tick track
  exists.
- External feature labels now create internal reservation slots. This avoids a
  hard outside-channel lock, but label measurement must be part of the same
  monotone packing pass as auto-track placement and must use the same label
  layout data that rendering consumes.
- `CircularResolvedSlot.payload` keeps renderer-specific resolved state out of
  `params`, but callers must not pattern-match on payload types outside the
  internal circular renderer path.
- Existing reference SVGs may change if preset expansion does not exactly
  reproduce fixed-oracle explicit feature/tick/numeric geometry. The preset
  oracle intentionally covers geometry only; final SVG/reference-output changes
  still depend on visual inspection and focused rendering tests.
- Feature layout depends on `feature_dict`, so measuring feature slots before
  final rendering must use the same precomputed feature dictionary.
- Tick labels have reserved footprints larger than tick marks; order tests
  should compare `packing_band`, and collision tests should compare
  `reserved_band`, not only draw bands.
- Signed inside/outside packing is easier to reason about, but migration code
  must carefully convert back to absolute radii for existing draw functions.
- Renaming `gap_after` to `spacing` touches Python parsing, dataclasses, web
  state, web serialization, and tests; missing one path can silently drop
  spacing.
- Splitting external label leaders from external label text requires care so
  embedded labels remain in the features slot rendering path.

## Non-Goals

- Do not redesign the web Custom Track Slots UI beyond geometry/default changes.
- Do not add `ri`/`ro` replacement syntax in this change.
- Do not add a second public radial geometry API alongside Track Slots.
- Do not preserve TrackSpec as a compatibility bridge.
- Do not migrate old Web-saved Custom Track Slots state.
- Do not preserve backwards compatibility for old circular Track Slot helper
  APIs, shortcut slot strings, or legacy circular custom slot grammar.
