# Implementation record: unified arrow geometry

Date: 2026-08-02

Issue:

- [#310: Feature shape: arrowhead](https://github.com/satoshikawato/gbdraw/issues/310)

Status: implemented

## Product decision

Model one directional `arrow` rendering. Do not expose separate shape names for
the full-width and narrow-shaft outlines. The visible topology follows from the
arrow geometry:

| Rendering | Directional | Geometry |
| --- | --- | --- |
| `arrow` | Yes | Global head-length and shaft-width ratios |
| `rectangle` | No | Foreground block |
| `underlay` | No | Full feature-band highlight behind foreground features |

CDS, rRNA, tRNA, tmRNA, ncRNA, and misc_RNA continue to use `arrow` by
default. `repeat_region` continues to use `underlay`; other types default to
`rectangle`.

The earlier branch-only `arrowhead` rendering value and
`--arrowhead_shaft_width_ratio` flag are superseded. They were not merged to
`main` or included in a release tag, so current branch-owned requests, sessions,
tests, examples, and documentation must be rewritten to the unified contract
rather than supported through a permanent migration.

## Geometry contract

Use one shared configuration block:

```toml
[objects.features.arrow_geometry]
head_length_ratio = "auto"
shaft_width_ratio = 1.0
```

- `head_length_ratio` is the longitudinal head length divided by the full
  rendered feature thickness. It accepts `auto` or a positive finite number.
- `shaft_width_ratio` is shaft thickness divided by full head thickness. It
  accepts a finite number in `(0, 1]` and applies to every `arrow` feature.
- `shaft_width_ratio = 1.0` preserves the legacy full-width arrow outline and
  exact default SVG path text.
- A ratio below `1.0` produces a centered narrow shaft and visible head
  shoulders. A long Linear arrow therefore has seven logical outline vertices;
  the vertex count is derived geometry, not a public shape choice.
- `head_length_ratio = "auto"` uses the existing Linear arrow-length parameter
  or Circular genome-length function as its full-width base, then adds the
  thickness removed from a narrowed shaft:

  ```text
  head_length_px = mode_default_px + full_feature_thickness_px * (1 - shaft_width_ratio)
  ```

  At `shaft_width_ratio = 1.0`, this is the exact legacy value.
- A numeric head length resolves in display space:

  ```text
  head_length_px = full_feature_thickness_px * head_length_ratio
  ```

  Circular rendering converts the display length to a floating-point genomic
  span at the feature center radius without rounding to whole base pairs.
- Cap head length at the terminal feature block. When no positive-length shaft
  remains, draw a triangle.
- Multipart arrows draw nonterminal blocks at configured shaft width, place the
  head on the strand-terminal block, and keep connectors on the center line.
- Features with undefined strand retain the rectangle fallback.

## Surface contract

### CLI

Both modes accept:

```text
--feature_shape TYPE=arrow
--arrow_head_length_ratio auto|POSITIVE_FLOAT
--arrow_shaft_width_ratio FLOAT_IN_(0,1]
```

The shaft option defaults to `1.0`. Remove the branch-only second-arrow shape
from feature-shape parsers and help text. Project both geometry flags to the
canonical `objects.features.arrow_geometry` leaves.

### Python API

`FeatureOptions.shapes` accepts `arrow`, `rectangle`, and `underlay`. Keep the
beginner and typed APIs small; callers tune geometry through `config_overrides`:

```python
config_overrides={
    "objects.features.arrow_geometry.head_length_ratio": 1.0,
    "objects.features.arrow_geometry.shaft_width_ratio": 0.5,
}
```

### Web UI

- Offer `Arrow`, `Rectangle`, and `Underlay` in each feature-type rendering
  selector.
- Keep one **Arrow Geometry** group with **Head Length Ratio** and **Shaft Width
  Ratio**.
- Explain that both settings apply globally to every Arrow feature type.
- Use `0.05` spinner increments and a GUI minimum of `0.05`; blank head length
  continues to mean Auto. The typed core continues to accept any valid positive
  finite value.
- Preserve the settings even when no feature type currently uses Arrow.

### Requests and sessions

Canonical request schema 5 and current session format already carry
`featureShapes` and arbitrary validated configuration overrides, so no format
version change is required. Current writers emit only the unified rendering
values and geometry state. Missing geometry in older released artifacts resolves
to `auto` and `1.0`.

The Web rendering-options capability advances to schema 3 because unified shaft
geometry is new. Its exact `featureRenderings` list is `arrow`, `rectangle`, and
`underlay`; app/wheel mismatches are rejected before rendering.

## Renderer routing

Keep the legacy full-width serializer as the `shaft_width_ratio == 1.0` path so
default reference SVG text remains unchanged. Route smaller ratios through the
shouldered-arrow builder. This removes a public shape branch without forcing a
new generalized serializer to rewrite legacy SVG command spelling.

Apply the same routing to:

- positive and negative strands;
- Linear and Circular diagrams;
- explicit lane/radius layouts and legacy factor-based layouts;
- terminal, multipart, and origin-spanning blocks;
- Circular embedded-label fit calculations.

## Documentation examples

Tutorial 9 keeps two Gallery-quality session variants. The Circular figure uses
a global shaft ratio of `0.75`; the Linear figure uses `0.5`:

```text
tutorial_9_arrow_geometry_circular
  -> examples/tutorial-9-arrow-geometry-circular.svg
tutorial_9_arrow_geometry_linear
  -> examples/tutorial-9-arrow-geometry-linear.svg
```

The source Gallery sessions, labels, legends, quantitative tracks, comparison
context, and other layout settings remain unchanged. Regenerate and visually
inspect both outputs before removing the superseded branch-only SVGs.

## Compatibility requirements

- Existing released `arrow` selections retain their meaning.
- Configs and sessions that omit arrow geometry keep the legacy full-width
  appearance through the new `1.0` default.
- Default output-comparison fixtures pass without regeneration.
- Current branch-owned artifacts use the newest unified writer representation;
  do not add a reader migration for the superseded development-only value.
- Renaming the CLI/Web state field must not rename the canonical dotted config
  leaf `objects.features.arrow_geometry.shaft_width_ratio`.

## Verification

Cover at least:

- accepted rendering values and unchanged feature-type defaults;
- config defaults, dotted overrides, and invalid values;
- exact legacy paths at shaft ratio `1.0`;
- centered shafts at `0.25`, `0.5`, and `0.9` in both diagram modes;
- short/equality boundaries and multipart features;
- CLI, Python, Web request, history, reset, and session round trips;
- runtime capability negotiation;
- GUI spinner increments and Auto restoration;
- both regenerated Tutorial 9 figures at a readable scale;
- the read-only default SVG comparison suite.

## Out of scope

- Per-feature-instance or qualifier-specific geometry.
- Per-feature-type shaft ratios.
- Applying feature-arrow geometry to comparison ribbons or annotation marks.
- Arbitrary user-defined polygons.
