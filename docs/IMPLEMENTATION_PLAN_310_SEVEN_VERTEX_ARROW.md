# Implementation plan: seven-vertex arrow shape

Date: 2026-07-30

Issue:

- [#310: Feature shape: arrowhead](https://github.com/satoshikawato/gbdraw/issues/310)

Status: planned

## Scope and product contract

Add `arrowhead` as a new feature rendering without changing the existing
`arrow`.

| Rendering | Geometry | Directional | Default use |
| --- | --- | --- | --- |
| `arrow` | Existing five-vertex arrow | Yes | Retain all current feature-type defaults |
| `arrowhead` | New seven-vertex arrow with a narrower shaft | Yes | Opt-in by feature type |
| `rectangle` | Existing rectangular block | No | Unchanged |
| `underlay` | Existing underlay layer policy | No | Unchanged |

The public values remain distinct even though both arrows share head-length and
strand calculations. `arrowhead` is not an alias for `arrow`, and `arrow` is not
reimplemented by sending a seven-vertex polygon with coincident points through
the new serializer. This separation keeps the existing five-vertex SVG path
data stable.

[PyGenomeViz's generalized Arrow patch](https://github.com/moshi4/pyGenomeViz/blob/fb137ded78d56b1d55229309d500aa54c80e3abe/src/pygenomeviz/patches.py#L10-L76)
is a useful reference for the shaft-ratio concept. Pixel-for-pixel PyGenomeViz
compatibility is not a goal.

Keep the default rendering assignments unchanged. CDS, rRNA, tRNA, tmRNA,
ncRNA, and misc_RNA continue to use `arrow`; `repeat_region` continues to use
`underlay`. Users select the new shape explicitly, for example:

```text
--feature_shape CDS=arrowhead
```

## Geometry and tuning contract

Add one shared configuration block:

```toml
[objects.features.arrow_geometry]
head_length_ratio = "auto"
shaft_width_ratio = 0.5
```

The fields have these meanings:

- `head_length_ratio` is the longitudinal head length divided by the full
  rendered feature thickness. It accepts `auto` or a positive finite number.
  It applies to both `arrow` and `arrowhead`.
- `shaft_width_ratio` is the shaft thickness divided by the full head
  thickness. It accepts a finite number in the interval `(0, 1]` and applies
  only to `arrowhead`.
- `arrow` always has an effective shaft ratio of `1.0`. Changing
  `shaft_width_ratio` must not alter its five-vertex outline.
- `auto` uses the current Linear arrow-length parameter and current Circular
  genome-length function. This preserves existing output and gives both arrow
  forms the same head length when they appear in one diagram.

For a numeric `head_length_ratio`, calculate the head length in display space:

```text
head_length_px = full_feature_thickness_px * head_length_ratio
```

Linear paths can use that value directly after coordinate normalization.
Circular paths convert it to an angular span at the feature center radius, then
to a floating-point genomic position. Do not round the converted value to an
integer number of base pairs.

The seven distinct vertices of a positive-strand Linear arrow are:

1. shaft start, upper edge;
2. shaft end, upper edge;
3. head base, upper full-width edge;
4. tip on the feature center line;
5. head base, lower full-width edge;
6. shaft end, lower edge;
7. shaft start, lower edge.

The negative-strand path mirrors this order along the x-axis. On a Circular
track, the shaft uses reduced inner and outer radii. Calculate those radii by
interpolating from the center radius toward the full inner and outer radii by
`shaft_width_ratio`. The head base spans the full track thickness, and the tip
lies on the center radius.

Apply these boundary rules:

- Cap the resolved head length at the length of the terminal feature block.
- Render a new `arrowhead` as a triangle when no positive-length shaft remains.
- Preserve the existing `arrow` short-feature condition and path formatting
  exactly, including its equality boundary.
- A configured `shaft_width_ratio` of `1.0` keeps the semantic value
  `arrowhead` but, given the same head length, produces the same visible
  outline as a five-vertex arrow. The new serializer may remove coincident
  points in this case.
- For multipart `arrowhead` features, draw nonterminal blocks at shaft width
  and draw the head only on the strand-terminal block. Keep connector lines on
  the feature center line.
- Retain the current rectangle fallback for features with an undefined strand.

## Implementation work

### 1. Preserve the selected glyph in the feature model

Edit:

- `gbdraw/features/shapes.py`
- `gbdraw/features/objects.py`
- `gbdraw/features/factory.py`
- `gbdraw/features/__init__.py`
- `gbdraw/configurators/features.py`

Changes:

- Extend `FeatureRendering` and `FEATURE_RENDERING_VALUES` with `arrowhead`.
- Treat both arrow values as directional in
  `resolve_directional_feature_types()`.
- Separate the foreground glyph (`arrow | arrowhead | rectangle`) from the
  `underlay` layer decision.
- Store the resolved glyph on each `FeatureObject`. The current factory loses
  this information by reducing `rendering == "arrow"` to `is_directional`.
- Make `is_directional` a derived compatibility property for
  `arrow | arrowhead`. Do not store a second Boolean source of truth.
- Keep the existing constructor form that accepts `is_directional` as a
  compatibility adapter for current internal and downstream callers. Map
  `True` to `arrow` and `False` to `rectangle`, and reject contradictory glyph
  and Boolean inputs.
- Pass the resolved glyph through all feature-object factory paths, including
  compound features and underlay partitioning.

### 2. Add and validate arrow geometry settings

Edit:

- `gbdraw/data/config.toml`
- `gbdraw/config/models/objects.py`
- `gbdraw/config/modify.py` only if schema-derived override discovery needs an
  update
- `gbdraw/configurators/features.py`

Changes:

- Add a frozen typed `ObjectsFeaturesArrowGeometryConfig` under
  `ObjectsFeaturesConfig`.
- Parse `head_length_ratio` as `auto` or a positive finite float.
- Parse `shaft_width_ratio` as a finite float in `(0, 1]`.
- Reject Booleans, zero, negative values, NaN, infinity, and a shaft ratio
  above `1`.
- Expose both leaves through the existing dotted configuration override
  mechanism.
- Pass the resolved typed settings to both feature drawers rather than reading
  configuration inside SVG path helpers.

### 3. Share metrics while retaining separate path topology

Edit:

- `gbdraw/svg/arrows.py`
- `gbdraw/svg/linear_features.py`
- `gbdraw/svg/circular_features.py`

Changes:

- Add pure helpers for resolved head length, shoulder position, short-feature
  classification, and shaft bounds.
- Keep mode conversion at the geometry boundary: Linear resolves pixels;
  Circular resolves an angular or base-pair span at the center radius.
- Route Circular label-fit checks through the same resolved head metric used by
  the drawer. `gbdraw/labels/circular.py` currently calls the automatic
  Circular arrow-length function independently.
- Rename the existing internal `*arrowhead*` helpers whose output is the legacy
  five-vertex `arrow`. Use unambiguous internal names such as
  `create_arrow_path_linear()` and `generate_circular_arrow_path()`.
- Retain thin compatibility wrappers for currently imported helper names if
  they are part of a supported Python import surface.
- Add dedicated seven-vertex builders for `arrowhead`.
- Keep a dedicated legacy serializer for `arrow`. Do not rebuild its SVG string
  from a generalized point list because numeric spelling, command spacing, or
  duplicate-point removal would change tracked output.

### 4. Render both glyphs in Linear diagrams

Edit:

- `gbdraw/render/drawers/linear/features.py`
- `gbdraw/svg/linear_features.py`
- `gbdraw/canvas/linear.py` only where the existing automatic head length is
  exposed to the shared resolver

Changes:

- Dispatch terminal blocks by `glyph_kind` rather than by
  `is_directional`.
- Preserve the existing five-vertex path for `arrow`.
- Generate the seven-vertex outline for `arrowhead` on positive and negative
  strands.
- Derive shaft y-coordinates from the resolved lane top, center, and bottom so
  separate strands, overlap lanes, and all Linear track layouts use the same
  ratio.
- Add a shaft-width rectangle path for nonterminal multipart blocks of
  `arrowhead`.
- Apply numeric head length after record normalization so diagrams with
  different genome lengths retain the requested visible ratio.

### 5. Render both glyphs in Circular diagrams

Edit:

- `gbdraw/render/drawers/circular/features.py`
- `gbdraw/svg/circular_features.py`
- `gbdraw/svg/arrows.py`
- `gbdraw/labels/circular.py`

Changes:

- Dispatch `arrow`, `arrowhead`, and `rectangle` explicitly for ordinary
  blocks, compound blocks, and coalesced origin-spanning blocks.
- Support both the legacy factor-based path and the explicit
  `inner_radius_px / center_radius_px / outer_radius_px` path.
- Use shaft radii for the body arcs and full feature radii for the two head-base
  corners.
- Retain the existing arc sweep directions, origin wrapping, and the split
  used for arcs over 20 degrees.
- Draw nonterminal multipart `arrowhead` blocks at shaft width.
- Keep undefined-strand features rectangular.
- Use the same resolved head length for drawing and label placement so a
  numeric ratio cannot make embedded-label decisions disagree with the visible
  head.

### 6. Expose the shape and ratios through CLI, Python, and sessions

Edit:

- `gbdraw/cli_utils/common.py`
- `gbdraw/circular.py`
- `gbdraw/linear.py`
- `gbdraw/interface.py`
- `gbdraw/api/options.py`
- `gbdraw/session_io.py`
- `gbdraw/session_request_codec.py`

Changes:

- Accept `TYPE=arrowhead` in every `--feature_shape` parser and help string.
- Add `--arrow_head_length_ratio` to both CLI modes. Accept `auto` or a
  positive finite number.
- Add `--arrowhead_shaft_width_ratio` to both CLI modes. Accept a finite number
  in `(0, 1]`.
- Project both flags to the canonical `objects.features.arrow_geometry`
  override leaves.
- Keep the beginner and typed Python APIs small. `FeatureOptions.shapes`
  selects `arrowhead`; callers tune geometry through the existing
  `config_overrides` mapping.
- Update `_feature_shapes_from_cli_args()` so CLI-created sessions do not drop
  the new rendering.
- Preserve the shape mapping and geometry overrides through canonical request
  encode/decode and CLI/GUI session replay.
- Validate the shape mapping at typed request construction instead of waiting
  for feature extraction.

Canonical request schema 5 and session version 39 already carry the shape map
and arbitrary validated configuration overrides. No format-version change is
expected. Missing geometry fields in an older config, request, or session
resolve to `auto` and `0.5`.

### 7. Add Web controls and runtime negotiation

Edit:

- `gbdraw/web/index.html`
- `gbdraw/web/js/utils/feature-rendering.js`
- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/feature-editor/rule-actions.js` only if its generic
  normalizer needs an update
- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web_support/capabilities.py`
- `gbdraw/web/js/services/runtime-capabilities.js`

Changes:

- Add `Arrowhead (7 vertices)` to the feature-rendering selector. Keep
  `Arrow (5 vertices)` as the selected default for current directional feature
  types.
- Add one `Arrow geometry` control group:
  - head length with `Auto` and a positive numeric ratio;
  - arrowhead shaft width as a ratio or percentage in `(0, 1]`.
- Keep the values available when no feature type currently uses `arrowhead`;
  changing feature assignments later should reuse the stored style.
- Emit the two canonical config override leaves and project them back to Web
  state.
- Preserve them through save/load, reset, history, and canonical artifact
  adoption.
- Add `arrowhead` to the expected `featureRenderings` capability.
- Advance `WEB_RENDER_OPTIONS_SCHEMA` from 1 to 2 because the Web app and wheel
  must agree on the expanded rendering options. Keep
  `WEB_RENDER_PROTOCOL` unchanged.

### 8. Update documentation after the implementation is verified

Edit:

- `docs/CLI_Reference.md`
- `docs/PYTHON_API.md`
- `docs/RECIPES.md`
- `docs/TUTORIALS/1_Customizing_Plots.md`
- `docs/TUTORIALS/9_Feature_Visibility_Shapes.md`

Document:

- The visual and semantic difference between `arrow` and `arrowhead`.
- Shape selection from CLI, Python, and Web.
- The definitions and valid ranges of both ratios.
- The `auto` behavior and short-feature triangle fallback.
- A side-by-side five-vertex and seven-vertex example using the same head
  length.

Update release notes for the target release. Do not regenerate unrelated
Gallery media or tracked default diagrams.

## Compatibility

- Default feature mappings remain `arrow`, so existing commands and API calls
  select the five-vertex shape.
- `head_length_ratio = "auto"` retains the current Linear and Circular
  head-length algorithms.
- The dedicated `arrow` serializers retain existing path text. Tracked
  reference SVGs must pass without regeneration.
- Existing config files and sessions omit the new geometry block and receive
  the defaults.
- Existing `is_directional` consumers continue to treat both arrow forms as
  directional.
- Old Web wheels do not understand `arrowhead`. The rendering option capability
  version and exact feature-rendering list must reject an app/wheel mismatch
  before rendering.
- `shaft_width_ratio = 1.0` changes only the new shape's visible topology. It
  does not rewrite the stored selection from `arrowhead` to `arrow`.

## Test plan

Python unit and integration tests:

- `tests/test_feature_shapes.py`
  - Parse and normalize all four public rendering values.
  - Keep current defaults.
  - Treat `arrow` and `arrowhead` as directional.
  - Preserve `glyph_kind` through factory creation and underlay partitioning.
  - Reject unknown values at CLI and typed API boundaries.
- `tests/test_config_override_paths.py`
  - Confirm both defaults and dotted override paths.
  - Accept `auto` and representative valid ratios.
  - Reject Booleans, nonfinite values, nonpositive values, and shaft ratios
    above `1`.
- Add `tests/test_linear_feature_paths.py`, or equivalent focused coverage.
  - Check exact existing `arrow` path strings.
  - Check seven distinct logical vertices for positive and negative strands.
  - Check numeric head ratios across normalization factors.
  - Check shaft ratios `0.25`, `0.5`, and `1.0`.
  - Check short and exactly-at-boundary features.
  - Check multipart shafts, separate strands, and overlap lanes.
- `tests/test_circular_feature_paths.py`
  - Cover positive and negative strands in both radius APIs.
  - Cover numeric and automatic head lengths.
  - Cover shaft ratios `0.25`, `0.5`, and `1.0`.
  - Cover short features, long arcs, origin-spanning features, multipart
    features, and undefined strands.
  - Confirm label-fit checks use the same resolved head span.
- `tests/test_session_io.py`
  - Preserve `arrowhead` and both CLI ratio flags through session conversion.
- `tests/test_session_request_codec.py`
  - Round-trip `featureShapes: {CDS: "arrowhead"}` and both config overrides
    without changing schema 5.
- `tests/test_web_runtime_capabilities.py`
  - Publish the expanded rendering list and option schema 2.
- `tests/test_output_comparison.py::TestOutputComparison`
  - Confirm every default reference SVG remains unchanged.

Web tests:

- `tests/web/feature-shapes.test.mjs`
  - Normalize `arrowhead`, retain defaults, and reject unknown values.
- `tests/web/session-request.test.mjs`
  - Emit and project the shape map and both geometry overrides.
  - Supply defaults when older sessions omit the values.
- `tests/web/runtime-capabilities.test.mjs`
  - Accept only the expanded option-schema-2 manifest.
- Add a focused browser check that selects `arrowhead` in both modes, changes
  both ratios, saves and restores the session, and verifies the resulting SVG
  path geometry.

Visual and performance checks:

- Add one small opt-in visual matrix with five-vertex and seven-vertex shapes,
  both strands, short and long features, and Linear and Circular modes.
- Inspect SVG and PNG output at ratios `0.25`, `0.5`, and `1.0`.
- Keep existing tracked references read-only. Add a new shape-specific fixture
  only if coordinate assertions do not provide a readable regression signal.
- Measure a large annotated record before and after the change. Keep rendering
  complexity linear in features plus location segments, with no extra full
  feature pass. Use a 5 percent median slowdown as a review signal, not a
  timing assertion in CI.

Focused verification:

```bash
pytest tests/test_feature_shapes.py -v
pytest tests/test_config_override_paths.py -v
pytest tests/test_linear_feature_paths.py -v
pytest tests/test_circular_feature_paths.py -v
pytest tests/test_session_io.py -v
pytest tests/test_session_request_codec.py -v
pytest tests/test_web_runtime_capabilities.py -v
node tests/web/feature-shapes.test.mjs
node tests/web/session-request.test.mjs
node tests/web/runtime-capabilities.test.mjs
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

Run the focused browser check after confirming the available Node or Python
Playwright installation.

## Acceptance criteria

- A long positive- or negative-strand `arrow` has the same five-vertex SVG path
  as before this change.
- A long `arrowhead` with the default shaft ratio has seven distinct logical
  vertices and mirrors correctly by strand.
- Numeric head ratios produce the requested visible proportion in Linear and
  Circular diagrams. `auto` reproduces the existing head-length behavior.
- Shaft ratios `0.25`, `0.5`, and `1.0` produce the specified proportion
  without moving the feature center line.
- Short features do not create a reversed or self-intersecting shaft.
- Multipart, origin-spanning, stranded, overlap-resolved, and undefined-strand
  features follow the geometry rules above.
- CLI, Python config overrides, Web controls, canonical requests, and saved
  sessions resolve to one typed geometry contract.
- Existing default output comparisons pass without updating their reference
  files.
- The Web app rejects a wheel that advertises the old rendering-option
  capability.

## Out of scope

- Changing any default feature type from `arrow` to `arrowhead`.
- Removing or aliasing the existing `arrow`.
- Per-feature-instance, qualifier-based, or table-driven ratio overrides.
- Arbitrary user-defined SVG polygons.
- Applying this shape to comparison ribbons, annotation marks, or nonfeature
  arrows.
- Pixel-for-pixel compatibility with PyGenomeViz.
- Changes to lane assignment, collision resolution, or label placement beyond
  sharing the resolved head length.

## Delivery and coordination

Implement #310 as one PR with separate review checkpoints for legacy output,
new geometry, and user-facing surfaces. The core geometry work can proceed
independently, but the Web and session portion overlaps the plans for #311,
#315, and #316 in `gbdraw/web/index.html`,
`gbdraw/web/js/services/session-request.js`, `gbdraw/web/js/services/config.js`,
runtime capability tests, and documentation. Rebase on whichever of those
changes lands first, then rerun the complete focused verification set.
