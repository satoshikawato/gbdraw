[Decision worksheet](./LINEAR_CIRCULAR_ARCHITECTURE_API_DECISION_WORKSHEET.md) | [Python API](./PYTHON_API.md) | [DiagramOptions audit](./DIAGRAM_OPTIONS_AUDIT.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md)

# Linear/Circular architecture and API audit

- Date: 2026-07-28
- Baseline HEAD: `544c10b`
- Branch: `2026-07-20`
- Audited snapshot: current working tree as observed on 2026-07-28
- Scope: Python public API, lower-level API, CLI, web app, request/session layer, configuration, layout, rendering, and contract tests

## Implementation status

This audit describes the 2026-07-28 baseline; the decision worksheet is the
current implementation record. Phase 0, Phase 1, the compatibility cleanup
approved under `O3.api=B`, and the A1/O4 delivery are complete. Phase 2 remains
a separately scoped follow-up. The findings and evidence below retain their
baseline wording so that the reasons for the completed changes remain
reviewable.

Typed requests use `CircularDiagramOptions` and
`LinearDiagramOptions`, with mode-specific `CircularTrackOptions`,
`LinearTrackOptions`, `CircularOutputOptions`, and `LinearOutputOptions`.
`CircularDiagramRequest` explicitly represents `single` and `grid`, while
`CircularBatchRequest` represents `batch` with one output per record.
`CircularRequestPlan`, `CircularBatchRequestPlan`, and `LinearRequestPlan` own
normalized builder selection; root API, fresh CLI/Web generation, current
canonical replay, and legacy internal replay reach those planners. A one-record
grid is valid.

Session version 38 and canonical request schema 5 persist the explicit Circular
grouping and batch output array. Record loading is mode-neutral, with topology
warnings and mode/comparison cardinality policy applied by planners. Active and
public runtime collinearity configuration uses
`LosslessCollinearityParameters`; supported canonical request schemas 1–5
privately migrate legacy `standard` parameter payloads while preserving their
effective fields.

## Executive summary

Linear and Circular should not share one geometry implementation. Polar bands,
inner/outer label placement, Cartesian rows, and pairwise ribbons are genuinely
different domains. The current problems are instead at their shared boundaries:
defaults, validation, request conversion, configuration ownership, SVG identity,
and prepared rendering data.

The audit found one output-affecting consistency issue that should be treated as a
release blocker pending confirmation of the intended defaults, four high-priority
contract or reproducibility issues, and several medium-priority architectural
debts:

1. The same mode does not have one default profile across Python, CLI, and the web
   app. This already changes comparison filtering and visible tracks.
2. There is no single canonical route from user input to a rendered diagram. The
   root API, lower-level API, CLI, web app, and session replay translate overlapping
   option models independently.
3. Mode safety is incomplete. Shared option bundles accept values that one mode
   ignores or rejects later, while comparison thresholds are validated differently
   by surface and mode.
4. Linear SVG generation uses process-randomized IDs and emits duplicate IDs in its
   dual legend. Circular already uses deterministic IDs.
5. Configuration and layout code has multiple sources of truth: raw dictionaries,
   typed models, flat override arguments, mutable canvas attributes, and web state
   caches.

The recommended direction is not a full Linear/Circular merge. It is one canonical,
mode-specific request contract feeding two mode-specific planners, with shared
validated value objects and rendering primitives where their semantics are actually
the same.

## Priority definitions

| Priority | Meaning in this audit |
|---|---|
| **P0 — Critical** | Can silently change scientific or visual output under a normal default workflow. Treat as a release blocker once the intended behavior is confirmed. |
| **P1 — High** | Breaks a public contract, reproducibility, or the ability to keep entry points equivalent. Resolve in the next API/rendering iteration. |
| **P2 — Medium** | Creates material drift risk or blocks safe extension, but the supported main path currently works. Schedule before adding more related options. |
| **P3 — Low** | Naming, compatibility, or documentation cleanup with a bounded current impact. |

## Prioritized findings

| ID | Priority | Finding | Main action |
|---|---|---|---|
| D1 | **P0** | Mode defaults differ across Python, CLI, and web app | Define one versioned default profile per mode and consume it from every surface |
| A1 | **P1** | Public and orchestration APIs have parallel request/option/render paths | Make typed requests the canonical boundary; reduce the other paths to adapters |
| A2 | **P1** | Mode-specific values and comparison thresholds are not validated consistently | Split or discriminate mode bundles and share value-object validation |
| R1 | **P1** | Linear SVG IDs are non-deterministic and duplicated | Use one stable ID factory and namespace dual-orientation elements |
| T1 | **P1** | Tests snapshot each surface but do not assert semantic parity | Add a reviewed cross-surface default and validation matrix |
| C1 | **P2** | Configuration has overlapping schemas and mode-sensitive global fields | Resolve one immutable mode profile from one typed configuration |
| L1 | **P2** | Placement, track-slot, annotation, and legend planning are duplicated inconsistently | Share parsers/planners; keep only geometry-specific solvers separate |
| L2 | **P2** | Circular layout types sit above their consumers, unlike Linear | Move reusable radial geometry contracts into `gbdraw.layout` |
| R2 | **P2** | Modern and legacy feature preparation can render the same feature differently | Require prepared feature layers inside render groups |
| R3 | **P2** | Shared style fields do not have shared rendering semantics | Apply and test the same style contract in both modes |
| W1 | **P2** | Web state has overlapping active and per-mode authorities | Keep one canonical per-mode preference model with computed active values |
| I1 | **P2** | I/O and output naming policy are coupled to mode inconsistently | Separate neutral loading, mode validation, and output policy |
| I2 | **P2** | Circular record collections have different output cardinality by surface | Model grid and batch output explicitly in the canonical request layer |
| U1 | **P3** | Legacy CLI flags and unused shared helpers obscure the intended contract | Deprecate aliases and generate mode-aware flags from declarative specs |

## D1 — P0: defaults are not owned by the mode

The inconsistency is not that Circular and Linear have different defaults. Different
defaults can be appropriate. The defect is that one mode changes behavior depending
on the entry point.

### Current comparison-filter and track defaults

The BLAST-filter columns apply when a comparison or conservation source is present.

| Surface | Circular BLAST filters (`evalue` / `bitscore` / `identity`) | Linear BLAST filters | Circular GC/skew | Linear GC/skew |
|---|---:|---:|---|---|
| Root/lower Python API | `1e-5` / `50` / `70` | `1e-5` / `50` / `70` | on | **on** |
| CLI | `1e-5` / `50` / `70` | `1e-2` / `50` / `0` | on | off |
| Fresh web state | **`1e-2` / `50` / `0`** | `1e-2` / `50` / `0` | on | off |

Evidence:

- The beginner API shares one `Thresholds` default between both modes:
  [`gbdraw/interface.py`](../gbdraw/interface.py#L89-L96) and
  [`gbdraw/interface.py`](../gbdraw/interface.py#L211-L226).
- The lower-level `DiagramOptions` repeats the same shared values:
  [`gbdraw/api/options.py`](../gbdraw/api/options.py#L99-L175).
- The common TOML enables GC and skew globally:
  [`gbdraw/data/config.toml`](../gbdraw/data/config.toml#L5-L12).
- Linear API assembly leaves those global values enabled unless explicit slots
  override them:
  [`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L1983-L2063).
- Linear CLI defines `1e-2`/`50`/`0` and positive `--show_gc`/`--show_skew` flags:
  [`gbdraw/linear.py`](../gbdraw/linear.py#L778-L850).
- Circular CLI uses `1e-5`/`50`/`70`:
  [`gbdraw/cli_utils/common.py`](../gbdraw/cli_utils/common.py#L220-L240).
- Web state starts with the Linear threshold values:
  [`gbdraw/web/js/state.js`](../gbdraw/web/js/state.js#L392-L462).
- The Circular web conservation execution path normalizes against
  `DEFAULT_LINEAR_BLAST_FILTERS` and then compares against Circular defaults. A fresh
  state with a conservation source therefore emits explicit
  `--evalue 1e-2 --identity 0` overrides:
  [`gbdraw/web/js/app/run-analysis.js`](../gbdraw/web/js/app/run-analysis.js#L3123-L3145).

The Python/CLI Linear difference changes both the set of accepted matches and the
presence of two tracks. The web Circular difference changes the comparison filter
relative to Circular CLI and Python. Both are silent and output-affecting.

The web app has a smaller instance of the same ownership problem: its default
feature list omits `misc_RNA`, although the CLI default list and the JavaScript
default-comparison helper include it
([`state.js`](../gbdraw/web/js/state.js#L392-L395),
[`cli-args.js`](../gbdraw/web/js/app/cli-args.js#L1-L9)).
Linear axis color has another visible drift: the Python path retains the TOML
`gray`, while CLI/web resolve an unspecified color to `lightgray`, or `dimgray`
when the ruler is on the axis
([`gbdraw/data/config.toml`](../gbdraw/data/config.toml#L165-L176),
[`gbdraw/linear.py`](../gbdraw/linear.py#L1501-L1521)).

### Recommended resolution

Create a single declarative schema with explicit `CircularDefaults` and
`LinearDefaults`. It should own comparison thresholds, visible default tracks,
feature types, label policy, title position, and output-name policy. Python should
load it directly; the web app should consume a generated/versioned representation.

The least surprising compatibility choice is probably:

- Circular: `1e-5` / `50` / `70`, GC/skew on.
- Linear: `1e-2` / `50` / `0`, GC/skew off.

Those values match the established CLI and current CLI documentation. Confirm this
as a product decision before changing the Python defaults because it is an
output-affecting API change.

## A1 — P1: there is no canonical rendering path

At least four overlapping contracts currently reach the renderers:

- The root `gbdraw` facade exposes mode-specific options, manually translates them
  into lower-level options, and calls builders:
  [`gbdraw/interface.py`](../gbdraw/interface.py#L427-L672).
- `gbdraw.api` exports builders, very large assemblers, canvases, configurators,
  typed requests, and session helpers together:
  [`gbdraw/api/__init__.py`](../gbdraw/api/__init__.py#L1-L184).
- Circular and Linear CLI handlers assemble and save directly, then reconstruct a
  typed request in request-compatible/session cases:
  [`gbdraw/circular.py`](../gbdraw/circular.py#L1113-L1374) and
  [`gbdraw/linear.py`](../gbdraw/linear.py#L1891-L2083).
- Session replay uses the separate `render_request()` route:
  [`gbdraw/api/request_render.py`](../gbdraw/api/request_render.py#L1092-L1297).

The lower-level public assemblers expose the cost of this structure. The Linear
assembler has about 70 parameters, while the Circular assembler exposes internal
parameters such as `_definition_profile`, `_precomputed_depth_df`, and
`_shared_depth_max`
([`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L1777-L1848),
[`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L2356-L2428)).

Output ownership is duplicated as well. `OutputOptions.output_prefix` influences
the drawing filename, while `RenderOutputRequest.output_prefix` controls saved
files. CLI request construction writes the same value into both locations
([`gbdraw/api/options.py`](../gbdraw/api/options.py#L55-L62),
[`gbdraw/api/requests.py`](../gbdraw/api/requests.py#L174-L233)).

This architecture requires every new option to be added to several dataclasses,
translation functions, parsers, serializers, and assembly signatures. D1 is a
concrete example of the resulting drift.

### Recommended resolution

Use this dependency direction:

```text
root API / CLI / web / session import
                  |
                  v
       canonical mode-specific request
                  |
                  v
     CircularPlanner or LinearPlanner
                  |
                  v
       mode-specific diagram assembler
                  |
                  v
         Drawing -> output writer
```

- Make `CircularDiagramRequest` and `LinearDiagramRequest` the canonical normalized
  boundary.
- Make `draw_circular`, `draw_linear`, CLI handlers, and web/session importers thin
  adapters to that boundary.
- Keep `assemble_*` functions internal or explicitly advanced and stop growing their
  public signatures.
- Let `RenderOutputRequest` alone own destination, prefix, format, and overwrite
  policy.
- Represent Circular's default batch-of-files behavior with a small batch request,
  rather than bypassing the canonical renderer.

This can be migrated incrementally. It does not require changing the two geometry
engines at the same time.

## A2 — P1: mode safety and threshold validation are incomplete

The root API is safer than the legacy options at the top level, but shared nested
types still advertise invalid combinations:

- `TrackOptions` stores Circular and Linear axes/slots together:
  [`gbdraw/api/options.py`](../gbdraw/api/options.py#L44-L52).
- `DiagramOptions` contains shared, Circular-only, and Linear-only fields in one
  dataclass. Its mode validator checks selected top-level fields, not nested track
  and output bundles:
  [`gbdraw/api/options.py`](../gbdraw/api/options.py#L99-L236).
- The Circular and Linear builders read only their own nested track fields, so a
  wrong-mode field can be silently ignored:
  [`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L3798-L3800) and
  [`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L3885-L3886).
- Shared `TitleOptions.position` accepts `none`, `center`, `top`, and `bottom`, while
  Circular rejects `center` and Linear rejects `none` later during rendering:
  [`gbdraw/interface.py`](../gbdraw/interface.py#L80-L86) and
  [`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L1358-L1377).
- Shared `DepthTrackOptions.height` is meaningful only for Linear and is rejected
  late for Circular:
  [`gbdraw/interface.py`](../gbdraw/interface.py#L99-L109).

Comparison validation is also asymmetric. Circular CLI validates identity, e-value,
and bitscore; Linear CLI accepts values such as `identity=101`, `identity=nan`,
negative e-value, and `bitscore=nan`
([`gbdraw/circular.py`](../gbdraw/circular.py#L574-L581),
[`gbdraw/linear.py`](../gbdraw/linear.py#L1138-L1213)).
The Python layers validate `alignment_length` but not the other threshold values.
An identity threshold of `100` can reach Linear gradient normalization with a zero
denominator
([`gbdraw/render/groups/linear/pairwise_match.py`](../gbdraw/render/groups/linear/pairwise_match.py#L276-L279)).

Circular layout validation has a similar timing mismatch.
`LinearMultiRecordOptions` validates finite, non-negative values in
`__post_init__`, while `CircularMultiRecordOptions` accepts invalid values until a
deep builder check
([`gbdraw/api/options.py`](../gbdraw/api/options.py#L64-L96),
[`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L1557-L1578)).

### Recommended resolution

- Define a shared validated `ComparisonThresholds` value object and use it at every
  CLI/API/request boundary.
- Split nested option bundles by mode, or use a discriminated mode type that
  validates every nested field at construction.
- Use `CircularTitleOptions` and `LinearTitleOptions`; keep their genuinely
  different value domains.
- Move all layout validation into immutable layout value objects and retain builder
  checks only as defensive assertions.
- Standardize Python/domain-boundary invalid-input failures on `ValidationError`;
  CLI adapters should translate them to `parser.error` or an argparse type error.

This finding extends the earlier
[`DiagramOptions` field audit](./DIAGRAM_OPTIONS_AUDIT.md): that audit established
top-level field ownership, while this audit covers nested bundles and parity across
surfaces.

## R1 — P1: Linear SVG identity is not deterministic or unique

Linear uses Python's process-randomized `hash()` for GC-skew clip IDs and pairwise
legend gradient IDs
([`gbdraw/render/drawers/linear/gc_skew.py`](../gbdraw/render/drawers/linear/gc_skew.py#L42),
[`gbdraw/render/groups/linear/legend.py`](../gbdraw/render/groups/linear/legend.py#L238-L240)).
The same input produced different hashes with `PYTHONHASHSEED=1` and
`PYTHONHASHSEED=2` during this audit.

Circular already uses a deterministic digest for skew clips and deterministic
sanitization for legend gradients
([`gbdraw/render/drawers/circular/gc_skew.py`](../gbdraw/render/drawers/circular/gc_skew.py#L41-L44),
[`gbdraw/render/groups/circular/legend.py`](../gbdraw/render/groups/circular/legend.py#L63-L66)).

Linear also keeps horizontal and vertical legends in the same SVG so the web app
can switch orientation. Both copies use `id="pairwise_legend"` and the same gradient
IDs
([`gbdraw/render/groups/linear/legend.py`](../gbdraw/render/groups/linear/legend.py#L242-L526)).
The dual layout is intentional; duplicate DOM IDs are not.

Consequences include byte-level non-reproducibility, invalid duplicate IDs, and
ambiguous CSS/DOM references. Existing SVG comparison utilities normalize dynamic
IDs, which makes this difficult to detect
([`tests/utils/svg_compare.py`](../tests/utils/svg_compare.py#L179-L190)).
The web app already post-processes gradient and clip references during SVG
repositioning
([`gbdraw/web/js/app/svg-styles.js`](../gbdraw/web/js/app/svg-styles.js#L59-L159)).
That reduces collisions in the web path, but it does not fix CLI/Python SVGs or the
duplicate `pairwise_legend` group IDs at the source.

Introduce one stable SVG ID factory based on semantic inputs. Namespace the two
Linear legend trees with orientation suffixes, including their gradient IDs. Add an
XML-level assertion that every emitted `id` is unique.

## T1 — P1: contract tests do not compare semantics across surfaces

`tests/test_public_contract.py` snapshots the root API and hashes Circular and
Linear parser metadata independently. That detects accidental signature changes but
cannot say whether equivalent defaults or validations agree
([`tests/test_public_contract.py`](../tests/test_public_contract.py#L47-L164)).

The web filter helper test demonstrates Linear values relative to Circular defaults,
but does not run from the real fresh web state. The session test manually supplies
the correct Circular values, bypassing the state that causes D1
([`tests/web/cli-args.test.mjs`](../tests/web/cli-args.test.mjs#L117-L135),
[`tests/web/session-request.test.mjs`](../tests/web/session-request.test.mjs#L283-L373)).
The beginner Python API test suite has a Circular end-to-end drawing test but no
paired default Linear test.

Add a reviewed semantic matrix with one row per mode and columns for:

- root Python options after compilation;
- lower-level typed request;
- CLI parser namespace;
- fresh web state to CLI arguments;
- fresh web state to canonical session request;
- expected visible default tracks and comparison thresholds;
- invalid threshold and wrong-mode option behavior.

The matrix should allow intentional mode differences, but require every surface
within one mode to agree.

## C1 — P2: configuration has multiple overlapping schemas

`CanvasConfig` places `show_gc`, `show_skew`, `show_depth`, and `show_labels` at the
global level. `show_labels` combines Circular's boolean semantics and Linear's
string policy in one union
([`gbdraw/config/models/canvas.py`](../gbdraw/config/models/canvas.py#L168-L200)).
A typo string is effectively on in Circular and off in Linear.

The same configuration then exists in several representations:

- nested TOML/raw dictionaries;
- `GbdrawConfig`;
- the roughly 90 flat parameters of `modify_config_dict`;
- public `apply_config_overrides`;
- constructors that accept both `config_dict` and optional `cfg`;
- resolved layout values dynamically attached to canvas configurators.

Representative dual-input constructors include
[`gbdraw/canvas/circular.py`](../gbdraw/canvas/circular.py#L45-L73),
[`gbdraw/canvas/linear.py`](../gbdraw/canvas/linear.py#L114-L143), and
[`gbdraw/configurators/features.py`](../gbdraw/configurators/features.py#L30-L70).
The flat override schema is maintained separately in
[`gbdraw/config/modify.py`](../gbdraw/config/modify.py#L51-L340).

There are also dead or suspect mappings:

- TOML `norm_factor`/`height` keys without corresponding typed object fields;
- typo `cicular_width_with_labels`;
- BLAST color overrides written to fields absent from the typed
  `ObjectsBlastMatchConfig`.

Resolve raw input once into `GbdrawConfig`, then derive an immutable
`CircularRenderProfile` or `LinearRenderProfile`. Internal constructors should
accept that typed profile only. Generate override paths from the typed schema and
reject unknown keys after a compatibility-warning period.

## L1 — P2: shared planning concepts are duplicated

Several pieces have the same syntax or semantic job but separate implementations:

- Circular and Linear track-slot parsers duplicate normalization and slot-head
  parsing. Circular's parse error lives in the shared parsing module, while Linear
  defines a separate error:
  [`gbdraw/tracks/circular.py`](../gbdraw/tracks/circular.py),
  [`gbdraw/tracks/linear.py`](../gbdraw/tracks/linear.py),
  [`gbdraw/tracks/parsing.py`](../gbdraw/tracks/parsing.py#L1-L15).
- The `<record-selector>@<row>` placement syntax has three validation paths:
  Circular API, Circular CLI, and the Linear parser reused by Linear API/CLI:
  [`gbdraw/api/diagram.py`](../gbdraw/api/diagram.py#L1086-L1193),
  [`gbdraw/circular.py`](../gbdraw/circular.py#L160-L189),
  [`gbdraw/layout/linear_multi_record.py`](../gbdraw/layout/linear_multi_record.py#L124-L220),
  [`gbdraw/linear.py`](../gbdraw/linear.py#L469-L474).
- Annotation auto-slot and underlay planning is near-duplicated. An unknown
  `set_id` raises `ValidationError` in Circular and `ValueError` in Linear:
  [`gbdraw/diagrams/circular/assemble.py`](../gbdraw/diagrams/circular/assemble.py#L129-L295),
  [`gbdraw/diagrams/linear/assemble.py`](../gbdraw/diagrams/linear/assemble.py#L146-L362).
- Legend dimension calculation dispatches on a configurator class-name string and
  then treats Circular and Linear asymmetrically:
  [`gbdraw/configurators/legend.py`](../gbdraw/configurators/legend.py#L59-L88).

Share the grammar, validated value objects, annotation binding, and immutable legend
measurement result. Keep radius/band resolution and row/height resolution in their
respective mode planners.

Label placement dictionaries also need explicit units before any common contract is
safe. Circular uses keys such as `middle`, `start`, and `end` for base-pair
coordinates, while Linear uses the same names for pixel geometry
([`gbdraw/labels/circular.py`](../gbdraw/labels/circular.py#L4644-L4796),
[`gbdraw/labels/linear.py`](../gbdraw/labels/linear.py#L490-L545)).
Mode-specific typed placed-label objects with `_bp` and `_px` field suffixes would
remove this ambiguity.

## L2 — P2: Circular geometry is in the orchestration layer

Reusable Circular geometry types such as `RadialBand` and
`CircularFeatureLayout` live in `gbdraw.diagrams.circular.radial_layout`
([`gbdraw/diagrams/circular/radial_layout.py`](../gbdraw/diagrams/circular/radial_layout.py#L49-L205)).
Lower rendering and label modules import upward from that package
([`gbdraw/render/groups/circular/seq_record.py`](../gbdraw/render/groups/circular/seq_record.py#L21),
[`gbdraw/render/groups/circular/labels.py`](../gbdraw/render/groups/circular/labels.py#L24),
[`gbdraw/labels/circular.py`](../gbdraw/labels/circular.py#L1243)).

Linear consumers instead depend on `gbdraw.layout.linear`. Circular's lazy package
initializer, compared with Linear's eager initializer, suggests possible import
cycle pressure; the upward imports themselves are the confirmed dependency problem
([`gbdraw/diagrams/circular/__init__.py`](../gbdraw/diagrams/circular/__init__.py#L1-L16),
[`gbdraw/diagrams/linear/__init__.py`](../gbdraw/diagrams/linear/__init__.py#L1-L8)).

Move radial value types and reusable geometry helpers to
`gbdraw.layout.circular`; leave diagram orchestration and solver coordination in
`gbdraw.diagrams.circular`.

## R2 — P2: prepared and legacy feature paths have different semantics

The modern `FeatureBuildResult` separates foreground features from underlays.
The compatibility `create_feature_dict` intentionally has no underlay input and
treats visible features as foreground
([`gbdraw/features/factory.py`](../gbdraw/features/factory.py#L143-L150),
[`gbdraw/features/factory.py`](../gbdraw/features/factory.py#L254-L328)).

Main Circular and Linear assemblers use the modern path. However, render groups
re-run the legacy function when precomputed features are omitted:

- [`gbdraw/render/groups/circular/seq_record.py`](../gbdraw/render/groups/circular/seq_record.py#L146-L172)
- [`gbdraw/render/groups/circular/labels.py`](../gbdraw/render/groups/circular/labels.py#L98-L117)
- [`gbdraw/render/groups/linear/seq_record.py`](../gbdraw/render/groups/linear/seq_record.py#L357-L381)

The supported main route is correct, but a direct/internal construction path can
turn `repeat_region` or another configured underlay into a foreground feature.
Require prepared `FeatureBuildResult` data in render groups and keep the legacy
adapter only at a compatibility boundary.

## R3 — P2: shared style options do not mean the same thing

`GcSkewConfigurator` exposes `stroke_width`. Linear skew applies it, while Circular
skew neither retains it nor sets it on the SVG path
([`gbdraw/configurators/gc.py`](../gbdraw/configurators/gc.py#L137-L139),
[`gbdraw/render/drawers/linear/gc_skew.py`](../gbdraw/render/drawers/linear/gc_skew.py#L17-L71),
[`gbdraw/render/drawers/circular/gc_skew.py`](../gbdraw/render/drawers/circular/gc_skew.py#L17-L61)).
GC content drawers in both modes retain a configured path width but do not apply it
to their main path. The default width of zero conceals the mismatch.

Either apply the shared style field in both modes or remove it from the public
contract. Add a paired render test with a non-default width; default-only snapshots
cannot detect this class of issue.

## W1 — P2: web layout preferences have multiple authorities

The web app keeps active legend/title values alongside Circular/Linear and
single/multi preference caches in
[`gbdraw/web/js/state.js`](../gbdraw/web/js/state.js#L213-L221).
Watchers, configuration import/export, and history snapshots each resolve or copy
overlapping values:

- [`gbdraw/web/js/app/watchers.js`](../gbdraw/web/js/app/watchers.js#L163-L217)
- [`gbdraw/web/js/services/config.js`](../gbdraw/web/js/services/config.js#L1170-L1268)
- [`gbdraw/web/js/services/history-snapshot.js`](../gbdraw/web/js/services/history-snapshot.js#L104-L174)

Preserving per-mode preferences is intentional. Maintain one canonical
`layoutPreferences` structure and expose the active mode as a computed projection.
Keep legacy reconciliation only in the import migration layer.

## I1 — P2: I/O and output policy are coupled to mode inconsistently

Public loaders take a free-form `mode` string, and the beginner
`read_genbank`/`read_gff` facade hardcodes `mode="linear"`
([`gbdraw/api/io.py`](../gbdraw/api/io.py#L48-L105),
[`gbdraw/interface.py`](../gbdraw/interface.py#L339-L364)).
Parsing records, validating mode-specific cardinality, selecting topology warnings,
and comparison-aware record selection/truncation therefore share one boundary.

Output policy differs legitimately—Circular often writes one accession-named file
per record, while Linear combines records into `out`—but the web help states one
record-name behavior for both modes
([`gbdraw/web/index.html`](../gbdraw/web/index.html#L1919-L1925),
[`gbdraw/circular.py`](../gbdraw/circular.py#L224-L228),
[`gbdraw/linear.py`](../gbdraw/linear.py#L761-L766)).

Implementation outcome: record loading is neutral, planners apply mode-specific
input policy, and adapters resolve output policy into typed output requests.

## I2 — P2: Circular collection cardinality differs by surface

A sequence of Circular records does not have one cross-surface meaning:

- Root `draw_circular(records)` automatically produces one grid when more than one
  record is passed:
  [`gbdraw/interface.py`](../gbdraw/interface.py#L608-L641).
- `CircularDiagramRequest` also inserts an automatic multi-record layout:
  [`gbdraw/api/requests.py`](../gbdraw/api/requests.py#L309-L335).
- Circular CLI writes separate diagrams/files by default and produces one grid only
  with `--multi_record_canvas`:
  [`gbdraw/circular.py`](../gbdraw/circular.py#L211-L216),
  [`gbdraw/circular.py`](../gbdraw/circular.py#L1113-L1294).

At the audited baseline, the root API rejected a `CircularLayout` for one record
while the typed request accepted it. The Phase 1 core resolved that mismatch:
both paths now accept a one-record layout and produce a 1×1 grid.

Implementation outcome: `CircularDiagramRequest` explicitly stores `single` or
`grid`, while `CircularBatchRequest` stores `batch` and preserves the CLI's
established separate-file behavior with one output per record. Schema 5 persists
the selected grouping instead of inferring it from the entry point.

## U1 — P3: compatibility CLI surface needs a controlled reduction

Circular uses negative GC/skew flags because both tracks default on; Linear uses
positive flags because they default off. The semantic difference is reasonable,
but a symmetric `--gc/--no-gc` and `--skew/--no-skew` contract would be easier to
generate and test while retaining old aliases.

Both modes also retain legacy `--depth` forms with different argument shapes beside
the normalized repeatable `--depth_track` form. The web app emits only
`--depth_track`. Deprecate `--depth` and keep one compatibility adapter at the CLI
boundary
([`gbdraw/circular.py`](../gbdraw/circular.py#L291-L301),
[`gbdraw/linear.py`](../gbdraw/linear.py#L787-L798),
[`gbdraw/cli_utils/common.py`](../gbdraw/cli_utils/common.py#L589-L616)).

Several exported parser helpers such as `add_output_args`, `add_stroke_args`, and
`add_legend_args` have no production callers, while handlers maintain manual
variants
([`gbdraw/cli_utils/common.py`](../gbdraw/cli_utils/common.py#L308-L456)).
After the canonical request migration, generate shared flag spelling and validation
from declarative option specs and delete unused wrappers.

## Intentional differences to keep

The following are not redundancies and should remain mode-specific:

- Circular polar bands versus Linear Cartesian rows and vertical tracks.
- Circular inner/outer/radial label collision handling versus Linear row and ribbon
  clearance.
- Circular conservation rings, ticks, and center reservation versus Linear
  nucleotide/protein comparisons and collinearity.
- Circular track width/radius values versus Linear track height/row values.
- Circular grid/radius scaling versus Linear shared-base-pair-scale multi-record
  layout.
- Linear dual-orientation legend support used by the web app.
- Different label policy vocabularies where the user concepts really differ.

The common layer should own semantic value objects, defaults, validation, stable
IDs, style application, annotation binding, and prepared-data boundaries. It should
not own one universal layout solver or renderer.

## Target architecture

```text
Declarative mode defaults + typed configuration
                       |
          +------------+------------+------------+
          |                         |            |
 CircularDiagramRequest  CircularBatchRequest  LinearDiagramRequest
          |                         |            |
  CircularRequestPlan  CircularBatchRequestPlan  LinearRequestPlan
          |                         |            |
  Prepared request       Prepared per-record items  Prepared request
          |                         |            |
  Circular renderers       Circular renderers   Linear renderers
          +------------+------------+------------+
                       |
                 Drawing/output
```

Shared, immutable contracts should include:

- comparison thresholds;
- feature visibility and prepared feature layers;
- record selectors and grid placements;
- annotation sources, styles, and validated bindings;
- stable SVG ID generation;
- output destination policy;
- common style primitives with tested semantics.

Mode planners should own:

- mode-specific default resolution;
- track and layout validation;
- title/definition placement;
- label geometry;
- comparison/conservation preparation;
- final mode-specific layout plans.

## Implementation status and remaining order

### Phase 0 — restore deterministic semantics (complete)

1. Fix the web Circular threshold default leak.
2. Choose and centralize the Linear Python default profile.
3. Validate all comparison thresholds through one value object.
4. Replace Linear `hash()` IDs and namespace dual legends.
5. Add the semantic contract matrix and unique-SVG-ID tests.

These changes address current output and reproducibility problems without a broad
refactor.

### Phase 1 — establish one request boundary (complete)

1. Route root API, CLI, web/session, and replay through canonical typed requests.
2. Move output destination ownership to `RenderOutputRequest`.
3. Keep large assembler signatures and mixed `DiagramOptions` as internal
   implementation bridges below the planners.
4. Add mode-specific nested options and eager layout validation.

Phase 1 is complete through `CircularDiagramOptions`, `LinearDiagramOptions`,
their mode-specific track and output bundles, and `CircularRequestPlan`,
`CircularBatchRequestPlan`, and `LinearRequestPlan`. Root API, fresh CLI/Web
generation, current canonical replay, and legacy internal replay all reach the
typed planners.

### Phase 2 — remove duplicated state and planners (separate follow-up)

1. Resolve one immutable mode render profile from typed configuration.
2. Consolidate record-placement grammar, annotation planning, and shared track-slot
   parsing.
3. Move Circular radial contracts into `gbdraw.layout`.
4. Require prepared feature layers in render groups.
5. Replace mutable canvas layout attributes with explicit layout/context objects.
6. Consolidate Web preference state into one per-mode authority with computed
   active values.

Phase 2 is outside the A1/O4 completion criteria. Before implementation, define
the code-reduction target and the specific superseded paths or contracts to
delete so that consolidation does not add another parallel layer.

### Phase 3 — compatibility cleanup (complete under `O3.api=B`)

1. Retired executable legacy GC/skew and Depth spellings and unused CLI/config
   helpers were removed. Supported persisted session data still migrates under
   `O3.data=A`.
2. Dead and typo configuration keys were removed from fresh executable and
   configuration inputs; supported persisted-data readers retain migration
   under `O3.data=A`.
3. Web help and output naming were made explicitly mode-aware.
4. `gbdraw.api` was reduced to the documented typed integration surface;
   low-level assemblers remain internal engine APIs.

## Verification performed

The audit used static call-path tracing and focused executable checks:

- compared Python option compilation and CLI parser defaults for both modes;
- traced fresh web state through Circular and Linear argument assembly;
- exercised invalid Linear/Circular threshold values;
- ran `node tests/web/cli-args.test.mjs`;
- ran `pytest tests/test_public_contract.py -q`;
- compared Python `hash()` output under two `PYTHONHASHSEED` values;
- inspected emitted Linear legend XML for duplicate IDs;
- traced modern and legacy feature preparation through both main assemblers and
  fallback render-group paths.

This document preserves the baseline audit and prioritization history. Owner
decisions, compatibility removals, and current implementation status are
recorded in the decision worksheet and the Implementation status section above.
The remaining architecture follow-up identified by this audit is Phase 2.

[Decision worksheet](./LINEAR_CIRCULAR_ARCHITECTURE_API_DECISION_WORKSHEET.md) | [Python API](./PYTHON_API.md) | [DiagramOptions audit](./DIAGRAM_OPTIONS_AUDIT.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md)
