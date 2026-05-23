# Label Rendering Policy Plan

## Status

Draft implementation plan for GitHub issue #202, "Optionally disallow
(non)embedded labels".

## Goal

Add an explicit label rendering policy that lets users keep the current
automatic behavior, show only labels that can be embedded in features, or force
labels out of feature bodies.

The default behavior must remain unchanged.

## Terminology

- Embedded label: a feature label drawn on or inside the feature body in the
  existing automatic placement path.
- External label: a non-embedded feature label drawn in an outside label track
  or circular leader-line layer.
- Above-feature label: the linear-only `--label_placement above_feature` mode.
  It currently uses `is_embedded=True` internally, but semantically it is a
  separate placement mode, not an embedded-label policy.

## User-Facing Option

Add a shared option:

```text
--label_rendering {auto,embedded_only,external_only}
```

Behavior:

- `auto`: current behavior. Labels that fit are embedded; labels that do not fit
  are placed externally.
- `embedded_only`: render only labels that the existing placement logic can
  embed. Drop labels that would require external placement.
- `external_only`: render labels externally, including labels that would
  normally be embedded.

Do not add primary negative flags such as `--no_embedded_labels` and
`--no_external_labels`. Two booleans can express an invalid "show no label
placement class" state and make CLI, Web UI, and API behavior harder to explain.
Aliases can be added later if there is a real usability need.

## Linear Above-Feature Relationship

`--label_placement above_feature` remains a separate linear label placement
mode. It supports `--label_rotation`, for example 45 degree slanted labels above
features, and should not be treated as either normal embedded or normal external
placement.

Initial rule:

```text
--label_rendering external_only|embedded_only cannot be used with
--label_placement above_feature
```

`--label_rendering auto` may be accepted with `above_feature`, but it should be a
no-op. The CLI and Web UI should avoid passing it.

Rationale:

- Existing `above_feature` labels are marked `is_embedded=True` for rendering
  and layout reuse.
- User-facing semantics are different: "above feature" means anchored near the
  feature with optional rotation, not "inside the feature body".
- Rejecting mixed modes keeps behavior predictable for the first
  implementation.

## Configuration Model

Add a top-level labels field in `gbdraw/data/config.toml`:

```toml
[labels]
rendering = "auto"
```

Add a typed field to `LabelsConfig`:

```python
rendering: Literal["auto", "embedded_only", "external_only"]
```

Parsing should normalize unknown or empty values to `"auto"` unless the caller is
CLI validation, where invalid CLI values should be rejected by `argparse`.

Relevant files:

- `gbdraw/data/config.toml`
- `gbdraw/config/models/labels.py`
- `gbdraw/config/modify.py`

## CLI Changes

Add `--label_rendering` to both circular and linear parsers.

Circular:

- Keep `--labels {none,out,both}` unchanged.
- `--label_rendering` only applies when `--labels` is not `none`.
- If labels are hidden, do not special-case beyond normal no-label behavior.

Linear:

- Keep `--show_labels {none,all,first}` unchanged.
- Add parser validation:
  - If `--label_placement above_feature` and `--label_rendering` is
    `embedded_only` or `external_only`, return a parser error.
  - `--label_rendering auto` is permitted but does not need to be passed from
    Web UI.

Relevant files:

- `gbdraw/circular.py`
- `gbdraw/linear.py`

## Placement Implementation

Centralize policy interpretation with a small helper rather than duplicating
string checks:

```python
def normalize_label_rendering(value: object) -> Literal["auto", "embedded_only", "external_only"]:
    ...
```

Potential location:

- `gbdraw/labels/policy.py`

### Circular

Modify `prepare_label_list()` after the existing `is_embedded` calculation.

Current decision point:

- `gbdraw/labels/circular.py`, around `_is_label_embedded_in_segment(...)`
- Embedded labels append to `embedded_labels`.
- Non-embedded labels append to `outer_labels` or `inner_labels`.

Policy behavior:

- `auto`: current branch unchanged.
- `embedded_only`: append only if `is_embedded` is true; skip otherwise.
- `external_only`: force `is_embedded = False` before the outer/inner branch.
  Preserve existing outer/inner routing, leader line geometry, and circular
  label collision handling.

### Linear

Modify `prepare_label_list_linear()` at the existing split:

- `force_above_feature`
- `bbox_width_px < longest_segment_length_in_pixels`
- external fallback

Policy behavior:

- `auto`: current behavior unchanged.
- `embedded_only`: keep only the existing embedded-fit branch; skip external
  fallback labels.
- `external_only`: route labels that fit inside features through the existing
  external-label branch instead of adding them to `track_0`.
- `above_feature`: reject non-auto policy before placement, so this function can
  keep the current `force_above_feature` code path unchanged.

Relevant files:

- `gbdraw/labels/circular.py`
- `gbdraw/labels/linear.py`

## Public API

Expose the option through existing configuration override plumbing first:

```python
config_overrides={"label_rendering": "external_only"}
```

This requires adding `label_rendering` to `modify_config_dict()`.

Do not add new top-level keyword arguments to
`assemble_circular_diagram_from_record()` or
`assemble_linear_diagram_from_records()` in the first pass. Those signatures are
already long, and `config_overrides` is the existing extension point.

Option bundles can add a labels-specific dataclass later if label options keep
growing.

## Web UI

Add state:

```javascript
adv.label_rendering = 'auto'
```

Add a select in the Labels card, visible when labels are enabled:

- Auto
- Embedded labels only
- External labels only

Linear UI rule:

- If `adv.label_placement === 'above_feature'`, hide or disable the rendering
  select and force `adv.label_rendering = 'auto'`.

Argument generation:

- Pass `--label_rendering <value>` only when value is not `auto`.
- Include feature detection in the existing option-support checks before
  passing the argument to the Pyodide wheel.

Config/session:

- Ensure `safeDeepMerge` can restore `adv.label_rendering`.
- Normalize unknown imported values to `auto`.
- Increment session version only if the session compatibility policy requires it
  for newly saved state.

Relevant files:

- `gbdraw/web/index.html`
- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`

## Documentation

Update CLI docs after implementation:

- `docs/CLI_Reference.md`
- `docs/TUTORIALS/3_Advanced_Customization.md`, if a short example is useful.

Suggested CLI examples:

```bash
gbdraw circular --gbk genome.gbk --labels --label_rendering embedded_only -o embedded_labels
gbdraw linear --gbk a.gbk b.gbk --show_labels all --label_rendering external_only -o external_labels
gbdraw linear --gbk a.gbk --show_labels --label_placement above_feature --label_rotation 45 -o slanted_labels
```

The third example should not include `--label_rendering`; it demonstrates the
separate above-feature mode.

## Tests

Add focused tests before broad reference-output updates.

Unit-style coverage:

- Circular `auto` returns both embedded and external labels when fixtures allow
  both.
- Circular `embedded_only` drops non-embedded labels.
- Circular `external_only` produces no embedded labels and still returns labels
  through the external pipeline.
- Linear `auto` preserves current embedded/external split.
- Linear `embedded_only` drops external fallback labels.
- Linear `external_only` routes fitting labels externally.
- Linear parser rejects `--label_placement above_feature --label_rendering external_only`.
- Linear parser rejects `--label_placement above_feature --label_rendering embedded_only`.

Reference SVGs:

- Existing reference outputs should not change for default `auto`.
- Add or update reference outputs only for new explicit modes if the existing
  test style requires full SVG comparison.

Likely files:

- `tests/test_circular_label_placement.py`
- A linear label placement test module, or a new
  `tests/test_label_rendering_policy.py` if shared coverage is cleaner.

## Rollout Order

1. Add config model support and `modify_config_dict()` mapping.
2. Add CLI parser options and validation.
3. Implement circular placement policy.
4. Implement linear placement policy and above-feature conflict validation.
5. Add tests for config parsing, CLI validation, and placement behavior.
6. Add Web UI state, controls, argument generation, and import normalization.
7. Update CLI docs and tutorial notes.
8. Run focused tests first, then broader label/web packaging checks if the Web UI
   was changed.

## Open Questions

- Should `external_only` circular labels allow inner labels when `--labels both`
  is active? Recommended: yes. It should force labels out of feature bodies, not
  override the existing outer/inner side policy.
- Should `embedded_only` in circular `--labels both` keep embedded labels on both
  strands? Recommended: yes. Embedded labels do not consume the inner/external
  label arena, so `both` should only affect external side routing.
- Should Web UI call this "Label Rendering" or "Label Fit Policy"? Recommended:
  "Label Rendering" in UI to match CLI, with concise help text explaining
  embedded vs external.
