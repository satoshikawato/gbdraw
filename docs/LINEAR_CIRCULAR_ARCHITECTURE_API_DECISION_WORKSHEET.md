[Architecture/API audit](./LINEAR_CIRCULAR_ARCHITECTURE_API_AUDIT.md) | [Python API](./PYTHON_API.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md)

# Linear/Circular architecture and API decision worksheet

- Date: 2026-07-28
- Status: owner decisions recorded; Phase 0 and approved beta cleanup implemented
- Source: [`LINEAR_CIRCULAR_ARCHITECTURE_API_AUDIT.md`](./LINEAR_CIRCULAR_ARCHITECTURE_API_AUDIT.md)

## How to respond

Use the decision codes in this document. The shortest valid response is:

```text
ALL=A
```

That approves every recommended option and delegates the implementation details
listed below. `ALL=A` selects A for every sub-decision as well.

To override only selected choices, use:

```text
ALL=A; O4.grouping=B; O5.style=C
```

To provide a custom value or defer one decision, use:

```text
O1.filters=C: Linear identity=20, keep the other recommended filter values
ALL=A; O3=HOLD: decide the compatibility window after the first implementation PR
```

`HOLD` leaves that decision and dependent work blocked; it is not an approval.
For positional shorthand, use:

```text
O1=<filters>/<tracks>/<features>/<axis>
O2=<boundary>/<validation>/<lowlevel>
O3=<api>/<data>/<cli>
O4=<grouping>/<single-layout>/<loading>/<naming>/<collision>
O5=<id>/<style>
O6=<gate>/<delivery>
```

No implementation that changes output defaults or output cardinality should begin
until O1 and O4 are decided.

## Decision summary

| Code | Owner decision | Recommended answer |
|---|---|---|
| O1 | Mode defaults and cross-surface parity | `O1=A/A/A/A` |
| O2 | Canonical boundary, validation, and low-level public surface | `O2=A/A/A` |
| O3 | Executable API, persistent data, and CLI migration | `O3=A/A/A` |
| O4 | Circular grouping, one-record layout, loading, naming, and collisions | `O4=A/A/A/A/A` |
| O5 | SVG ID and exposed style contracts | `O5=A/A` |
| O6 | Release gate and delivery strategy | `O6=A/A` |

## O1 — Mode defaults and cross-surface parity

Related findings: D1, T1.

Fresh Python, CLI, and web requests should resolve the same profile for a mode.
Newly saved sessions should store either the profile version or resolved values.
Replay must preserve the serialized semantics of supported existing sessions
instead of reinterpreting them with a new default.

O1 has four independent sub-decisions.

### O1.filters — Comparison and conservation filters

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Established mode-specific CLI values (recommended)** | Circular: e-value `1e-5`, bitscore `50`, identity `70`, alignment length `0`. Linear: `1e-2`, `50`, `0`, `0`. | Changes Python Linear defaults and fixes the web Circular threshold leak. It preserves established CLI behavior and current CLI documentation. |
| **B — Current Python common values** | Both modes: `1e-5`, `50`, `70`, `0`. | Changes Linear CLI/web comparison results. |
| **C — Custom values** | Specify all four values for both modes. | Impact depends on the supplied values. |

### O1.tracks — Default GC content and skew visibility

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Current mode conventions (recommended)** | Circular GC/skew on; Linear GC/skew off. | Changes Python Linear defaults while preserving CLI/web conventions. |
| **B — Both modes on** | Circular and Linear GC/skew on. | Changes Linear CLI/web output. |
| **C — Custom visibility** | Specify GC content and skew separately for each mode. | Impact depends on the supplied values. |

### O1.features — Default feature list

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Include `misc_RNA` in both modes (recommended)** | Use the current CLI/Python feature set in every fresh surface. | Changes fresh web output when `misc_RNA` is present. |
| **B — Exclude `misc_RNA` in both modes** | Use the smaller current web list everywhere. | Changes CLI/Python output when `misc_RNA` is present. |
| **C — Mode-specific list** | Specify the complete default list for each mode. | Preserves a deliberate mode difference but adds another profile dimension. |

### O1.axis — Linear axis color

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Current CLI/web rule (recommended)** | Use `lightgray`, or `dimgray` when the ruler is on the axis. | Changes the Python default axis color. |
| **B — Current TOML/Python rule** | Use `gray`. | Changes CLI/web Linear output. |
| **C — Custom rule** | Specify the color or conditional rule. | Impact depends on the supplied rule. |

Every `O1.filters` selection also fixes the web Circular path so it loads the
selected Circular profile rather than Linear defaults.

### Response

Use `O1=A/A/A/A`, or answer sub-decisions separately:

```text
O1.filters=A; O1.tracks=A; O1.features=A; O1.axis=A
```

## O2 — Canonical API and rendering boundary

Related findings: A1, A2, C1, R2.

O2 has three independent sub-decisions.

### O2.boundary — Canonical normalized input and output ownership

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Layered facade and typed request boundary (recommended)** | Keep `gbdraw.draw_circular` and `gbdraw.draw_linear` as the supported user-facing API. Use typed mode requests as the canonical normalized transport/session boundary. Route root API, CLI, web app, and replay through mode planners. Inside typed-request/CLI/session execution, `RenderOutputRequest` alone owns destination, filename prefix, formats, and overwrite policy. Root in-memory drawing remains output-neutral and retains explicit `Diagram.save(path, overwrite=...)` and `Diagram.to_bytes(format)` calls. | Existing root calls remain adapters. `DiagramOptions.output.output_prefix` and duplicate request fields follow O3. |
| **B — Root drawing options are canonical** | Make `CircularOptions` and `LinearOptions` canonical for normalized drawing settings and generate transport/session requests from them. Root in-memory drawing remains output-neutral. A separate `RenderOutputRequest` remains the sole owner of destination, filename prefix, formats, and overwrite policy inside typed-request/CLI/session saving. | Keeps one beginner drawing model, but expands it to cover session and batch adaptation. |
| **C — Co-equal normalized models** | Keep root options, `DiagramOptions`, typed requests, and low-level assemblers as independent stable contracts. Synchronize them with tests. | Maximizes immediate compatibility but preserves the translation and drift burden. |

### O2.validation — Mode typing and public numeric domains

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Mode-specific types and eager strict validation (recommended)** | Use mode-specific title, track, depth, and layout bundles. Require finite e-value and bitscore `>= 0`, finite identity in `0..100` inclusive, and integer alignment length `>= 0`. Identity `100` is valid: retained 100% hits use the maximum-identity color and a zero-span legend is handled safely. Wrong-mode values fail at construction in the new contract. | Invalid numeric values become immediate hard errors. Legacy shared types use O3 warnings before wrong-mode fields become errors or are removed. |
| **B — Shared discriminated model with the same strict domains** | Retain one model with a required mode discriminator and validate every nested field and the same numeric domains. | Preserves more class/import compatibility but provides weaker mode-specific completion and discoverability. |
| **C — Current permissive model** | Keep shared option bags and builder-time checks, including current mode and surface asymmetries. | Preserves current acceptance but leaves silent ignores and late failures. |

Invalid, non-finite, or out-of-domain numeric thresholds are not treated as valid
legacy workflows under A or B; they become hard errors without a warning period.
Domain `ValidationError` values are translated to argparse errors at the CLI edge.

### O2.lowlevel — Final status of low-level surfaces

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Narrow stable API (recommended)** | Stable surfaces are the root facade and typed request/render/session contracts. Deprecate public re-exports of `DiagramOptions`, `build_*`, `assemble_*`, canvas configurators, drawing configurators, and `OutputOptions.output_prefix`, then remove or relocate them at the O3 breaking release. Render groups and Circular radial-layout internals remain internal and may move without a public deprecation promise. | Deliberately reduces `gbdraw.api.__all__`; current public low-level imports receive a migration window. |
| **B — Documented advanced API** | Keep high-level builders, canvases, and configurators as a supported advanced namespace; deprecate only giant assemblers and duplicate output-prefix ownership. | Retains a larger long-term compatibility surface. |
| **C — Preserve current exports** | Keep all current `gbdraw.api` re-exports as stable-ish public API. | Maximizes compatibility but constrains the architecture cleanup. |

Under A or B, typed `GbdrawConfig` plus an immutable mode render profile becomes the
internal configuration authority.

### Response

Use `O2=A/A/A`, or answer:

```text
O2.boundary=A; O2.validation=A; O2.lowlevel=A
```

## O3 — Compatibility and migration policy

Related findings: A1, A2, C1, L2, R2, W1, U1.

O3 has separate policies for executable APIs, persistent data, and the target CLI
contract.

### O3.api — Python/CLI deprecation schedule

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Staged deprecation (recommended)** | Warn for at least one normal released cycle and continue the adapter until an explicitly announced breaking release. Record the first warning release and planned removal release when implementation starts. Unknown, dead, or typo config keys follow the same schedule before becoming errors. | Existing valid scripts and custom config continue to run with warnings and a documented deadline. |
| **B — Immediate beta cleanup** | Remove obsolete imports, fields, aliases, module paths, and config keys in the next release. | Produces the smallest codebase quickly but breaks users without a warning cycle. |
| **C — Indefinite compatibility** | Keep executable adapters and aliases without a removal target. | Maximizes compatibility but preserves audited redundancy. |

When `O2.validation` is A or B, new strict types reject wrong-mode fields
immediately, while legacy shared option types warn during the deprecation period.
Invalid numeric thresholds follow the selected O2 validation contract; under A or
B they are immediate errors, not deprecated valid behavior.

### O3.data — Session and standalone configuration retention

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Latest writer, compatibility-series readers (recommended)** | Structural session changes bump the session version. Writers emit only the latest version. Readers preserve every session version currently supported through the declared compatibility series and at least until the next announced breaking release; versions not currently supported gain no new promise. Older CLI-replay-only versions retain their existing replay capability but do not automatically gain typed-request conversion. Recognized legacy standalone web-config payloads migrate on import, and new standalone saves use a versioned envelope. Removing a reader requires a separate breaking decision or an offline migration tool. | Protects persisted user data longer than executable aliases while avoiding an unlimited promise to every historical payload. |
| **B — Rolling reader window** | Write the latest version and support only the latest plus two prior released schema versions. | Bounds maintenance but can strand older saved sessions/configs. |
| **C — Latest version only** | Require users to migrate persisted data before upgrading. | Smallest reader surface and highest user-data compatibility risk. |

The selected O1 defaults apply only to fresh requests. Supported existing sessions
replay their serialized resolved semantics.

### O3.cli — Target CLI contract

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Symmetric booleans and canonical depth tracks (recommended)** | Make `--gc/--no-gc`, `--skew/--no-skew`, and repeatable `--depth_track` canonical. Retain current positive/negative GC/skew flags and `--depth` as aliases under O3.api. Verify that each legacy `--depth` shape can be converted without loss before announcing deprecation. | Gives both modes one modern spelling while retaining scripts during the warning period. |
| **B — Canonicalize depth only** | Standardize on `--depth_track` but retain the current mode-specific GC/skew flag names as the permanent contract. | Removes the largest input-shape redundancy with less CLI churn. |
| **C — Keep the current CLI contract** | Retain all current flag names and both depth forms as first-class behavior. | Avoids CLI changes but preserves parser/help duplication. |

### Response

Use `O3=A/A/A`, or answer:

```text
O3.api=A; O3.data=A; O3.cli=A
```

A custom API schedule may be supplied as:

```text
O3.api=CUSTOM: warn in <release>; remove in <breaking release>
```

## O4 — Circular grouping, loading, and output naming

Related findings: A1, I1, I2.

O4 has five related sub-decisions. Concrete class names follow O2; O4 defines
semantics rather than implementation names.

### O4.grouping — Meaning of a Circular record collection

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Explicit grid/batch semantics with current facade defaults (recommended)** | The canonical model distinguishes a one-grid request from a separate-diagram batch request. Root API adapts a record sequence to grid; CLI adapts its current default to batch. The normalized request always states the choice. | Preserves current root and CLI output cardinality while making the difference intentional and serializable. |
| **B — Grid everywhere** | Every surface treats multiple Circular records as one grid unless `separate` is requested. | Changes the CLI default from multiple files to one diagram. |
| **C — Separate everywhere** | Every surface treats multiple Circular records as a batch unless `grid` is requested. | Changes root API and typed-request return/output behavior. |

### O4.single-layout — One-record grid/layout

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Allow a 1×1 grid (recommended)** | A one-record grid/layout is valid and produces one diagram. | Removes the current root/request validation mismatch without reducing capability. |
| **B — Reject it everywhere** | Layout/grid controls require at least two records. | Preserves the current root restriction and requires the typed request to reject an existing accepted form. |

### O4.loading — Record-loading responsibility

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Mode-neutral loading (recommended)** | Parse records without a drawing mode. Perform topology warnings, cardinality checks, comparison-aware selection, and mode validation in the request planner. Retain old `mode=` loaders temporarily under O3. Neutral selectors, limits, and streaming/lazy implementations remain allowed; this decision does not require materializing every record. | Keeps parsed data reusable, removes the root facade's hardcoded Linear mode, and preserves a path for large-input performance. |
| **B — Mode-aware loading** | Keep mode as a loader input, but replace free-form strings with a validated mode enum. | Smaller change, but preserves I/O/render-policy coupling. |

### O4.naming — Default output names

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Preserve exact CLI fallbacks and document them (recommended)** | Circular batch without `-o` uses each `SeqRecord.id`; with one record and `-o PREFIX` it uses `PREFIX`, and with multiple records it uses `PREFIX_1`, `PREFIX_2`, ... . Circular grid without `-o` uses the first `SeqRecord.id`; with `-o`, it uses the prefix unchanged. Linear CLI without `-o` uses `out`; with `-o`, it uses the prefix unchanged. Root Python drawing has no implicit filename and continues to require a path in `Diagram.save(path)`. Make web help mode-aware. | Minimizes filename and automation changes. |
| **B — Derive names from records in every mode** | Batch uses each `SeqRecord.id`. A one-record grid/Linear output uses that ID. A multi-record grid/Linear output uses `<first-id>_plus_<remaining-count>_<digest>`, where `digest` is a stable eight-character digest of the ordered record IDs. An explicit prefix overrides the derived name unchanged. | Changes existing Linear filenames and may affect automation. |
| **C — Require a prefix at writer boundaries** | CLI and typed output writers require a prefix; in-memory drawing remains output-neutral. Define one explicit batch prefix as the base for deterministic item names. | Clearest writer contract, but the broadest CLI usability and compatibility change. |

Under A or B, adapters resolve these names before constructing typed output
requests; the canonical drawing model does not infer a surface from the filename.

### O4.collision — Duplicate output names

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Deterministic suffixes, never silent overwrite (recommended)** | Resolve duplicate implicit item names by input order with `_2`, `_3`, ... before the extension. Explicit multi-record prefixes retain `_1`, `_2`, ... . Existing-file overwrite remains governed separately by the output overwrite policy. | Prevents duplicate `record.id` values from silently targeting the same filename. Some previously colliding batches receive new names. |
| **B — Reject duplicate resolved names** | Raise an error and require the caller to supply unique names or a prefix. | Strictest behavior, but interrupts batches that could be named safely. |

### Response

Use `O4=A/A/A/A/A` for all recommendations, or answer sub-decisions separately:

```text
O4.grouping=A; O4.single-layout=A; O4.loading=A; O4.naming=A; O4.collision=A
```

## O5 — SVG ID and exposed style contracts

Related findings: R1, R3.

O5 has two independent sub-decisions.

### O5.id — Stability of SVG DOM IDs

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Version-local determinism and explicit semantic hooks (recommended)** | Every emitted ID is always unique. For the same gbdraw version, input, and resolved configuration, IDs are deterministic. Exact internal ID spelling is not guaranteed across versions. Namespace Linear dual-orientation IDs. Before removing or changing semantic hooks, inventory and separately document the supported interactive-SVG metadata/role contract; static SVG promises only valid unique/deterministic IDs. | Fixes reproducibility and invalid duplicate IDs. First-party web selectors, legend editing/repositioning, and SVG post-processing must migrate with the renderer. External CSS/JS using undocumented internal IDs may also need migration. |
| **B — Freeze a reviewed compatibility set** | Inventory deterministic, non-colliding top-level IDs and prefixes used by repository-shipped first-party consumers plus any external selector list supplied by the owner with this decision. Freeze that set, replace randomized suffixes, and add orientation namespaces only where collisions require them. If no external list is supplied, the compatibility scope is first-party consumers only. | Reduces external selector changes but requires explicit inventory tests. Randomized or duplicate values cannot be preserved exactly. |
| **C — Treat every emitted ID as public** | Freeze all SVG ID names as a supported API after repairing collisions. | Creates a large long-term rendering compatibility surface. |

### O5.style — Meaning of exposed `stroke_width`

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Common outline semantics (recommended)** | For both GC content and GC skew in both modes, `stroke_width=0` means no visible outline. A positive width draws the path outline using the configured `stroke_color`; `stroke_color="none"` still means no visible outline. Fill color and opacity semantics remain unchanged. Preserve current configured defaults. | Non-default configurations begin working consistently. Even when the default appearance is unchanged, SVG markup may change; reference SVGs require review. |
| **B — Deprecate and remove the field** | Stop advertising the shared field and retain fill-only rendering where appropriate. | Breaks callers that expect the option to have an effect. |
| **C — Split the field by mode/object** | Define separate Circular/Linear GC content/skew stroke contracts. | Most explicit, but expands the API and configuration schema. |

### Response

Use `O5=A/A`, or answer `O5.id=<A|B|C>` and `O5.style=<A|B|C>`.

## O6 — Release gate and implementation scope

Related findings: all, especially D1, R1, T1, and U1.

O6 separates release gating from delivery strategy.

### O6.gate — Blocking contract checks

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Blocking gate (recommended)** | Do not release until approved O1 defaults, the selected O2 validation contract, deterministic/unique IDs, and the semantic parity matrix pass. Intentional surface differences require a reviewed allowlist. | Prevents known silent semantic drift from being released again. |
| **B — Non-blocking monitoring** | Add the same checks but report failures without blocking release. | Allows faster release but does not enforce the approved contract. |

Additions to the parity exception allowlist require owner approval under A.

### O6.delivery — Implementation sequence

| Option | Policy | Compatibility impact |
|---|---|---|
| **A — Staged implementation (recommended)** | Phase 0: D1, the selected O2 validation behavior needed for current entry points, R1, and T1. Next API iteration: A1/A2 and the full O2 boundary/type migration. Before adding related capabilities: C1/L1/L2/R2/W1/I1/I2. Perform removals only under O3. | Controls output risk first and avoids one high-risk refactor. |
| **B — Phase 0 only** | Fix defaults, implement/test the selected O2 validation behavior for current entry points, fix IDs, and add tests; defer the canonical API and all later architecture phases without a committed schedule. | Smaller immediate scope, but A1 and related translation debt remain unresolved even though A1 is high priority. |
| **C — One coordinated refactor** | Implement all audit phases in one breaking change. | Shortest transition period but highest review, regression, and migration risk. |

Test fixtures, parser extraction, class names, internal module layout, and PR
splitting remain maintainer decisions.

### Response

Use `O6=A/A`, or answer:

```text
O6.gate=A; O6.delivery=A
```

## Delegated implementation package

If the owner replies `ALL=A`, the following details are delegated to maintainers:

| Audit ID | Delegated direction |
|---|---|
| D1 | Represent the approved O1 profiles once and consume them from Python, CLI, web, and session paths. |
| A1 | Implement adapters, mode planners, and single output ownership under the O2 contract. |
| A2 | Implement the selected O2 validation contract, eager layout validation, and domain-to-CLI error translation. |
| R1 | Implement the deterministic ID factory, orientation namespace, and migration from undocumented internal selectors under O5. |
| T1 | Implement the semantic matrix, invalid-input matrix, paired default-render tests, and unique-SVG-ID assertion. Keep intentional differences in a reviewed allowlist. |
| C1 | Resolve raw input once into typed configuration and immutable mode profiles. Generate override paths from the schema rather than maintaining another flat schema. |
| L1 | Share grammars, value objects, annotation binding, and immutable legend measurement. Keep Circular/Linear geometry solvers and label placement separate. |
| L2 | Move reusable Circular radial contracts into `gbdraw.layout.circular`; use the O3 re-export policy if an old internal path needs a transition. |
| R2 | Require prepared `FeatureBuildResult` data inside render groups and isolate the legacy conversion at a compatibility boundary. |
| R3 | Apply or migrate the shared style field according to O5 and add non-default paired render tests. |
| W1 | Preserve the current separate Linear, Circular-single, and Circular-multi preference behavior in one canonical state tree; expose active values as a computed projection and apply the O3 data-migration policy. |
| I1 | Separate neutral record loading from planner validation and implement the exact O4 naming rules. |
| I2 | Implement explicit grid and batch request forms and persist the grouping choice; concrete class names follow O2. |
| U1 | Implement the selected O3 CLI contract from declarative option specifications; delete unused helpers only under the O3 API schedule. |

If any implementation uncovers a new public behavior change, it returns to the
owner as a new numbered decision instead of being inferred.

## Audit-to-decision mapping

| Audit finding | Owner decision | Delegated implementation |
|---|---|---|
| D1 | O1, O6 | Default-profile representation and parity tests |
| A1 | O2, O3, O4 | Adapter and planner structure |
| A2 | O2, O3 | Validation and exception implementation |
| R1 | O5, O6 | Stable ID factory and orientation namespace |
| T1 | O1, O6 | Test matrix and fixtures |
| C1 | O2, O3 | Typed mode profiles and schema-derived overrides |
| L1 | O6 for schedule only | Shared parser/planner extraction |
| L2 | O3 for path compatibility; O6 for schedule | Module relocation |
| R2 | O2, O3 | Prepared-data boundary |
| R3 | O5 | Style application and paired render tests |
| W1 | O3 for saved-state migration; O6 for schedule | Canonical web preference state and migration |
| I1 | O4 | Neutral loader and output-policy adapters |
| I2 | O4 | Explicit grid/batch requests |
| U1 | O3, O6 | CLI aliases, warnings, and helper cleanup |

## Decision record

Decision date: 2026-07-28

Normalized owner response:

```text
ALL=A; O3.api=B
```

| Code | Status | Selected option | Notes |
|---|---|---|---|
| O1 | Approved | `A/A/A/A` | Use the recommended mode defaults and cross-surface parity rules. |
| O2 | Approved | `A/A/A` | Use the recommended canonical boundary, validation contract, and low-level public surface. |
| O3 | Approved | `B/A/A` | Remove obsolete executable API and CLI compatibility paths in the next release without a warning cycle. Continue the recommended compatibility-series reader policy for persisted session/config data, and adopt the recommended canonical CLI contract. |
| O4 | Approved | `A/A/A/A/A` | Use the recommended Circular collection, loading, naming, and collision policies. |
| O5 | Approved | `A/A` | Use the recommended SVG ID and style contracts. |
| O6 | Approved | `A/A` | Apply the recommended release gate and staged delivery plan. |

## Implementation record

Implemented on 2026-07-28:

- Phase 0 now uses one versioned mode-profile source for Python, CLI, Web, and
  newly saved sessions; current entry points enforce the approved numeric and
  wrong-mode validation rules.
- SVG IDs are deterministic, valid, and collision-safe. First-party consumers
  select records, definitions, tracks, legends, and orientations through
  semantic `data-gbdraw-*` hooks. GC content/skew `stroke_width` now has the
  same outline meaning in both modes.
- Under `O3.api=B`, obsolete executable compatibility was removed immediately:
  the old GC/skew and Depth CLI spellings, unused CLI/config helpers, dead
  config keys, the `gbdraw.api.canvas`, `gbdraw.api.configurators`, and
  `gbdraw.circular_diagram_components` compatibility modules, obsolete
  low-level `gbdraw.api` re-exports, `plot_circular_diagram` /
  `plot_linear_diagram` save wrappers, package-level render aliases, the
  CairoSVG availability proxy, no-op label adapters, zero-caller helpers, and
  duplicate `OutputOptions.output_prefix`. CLI export implementations remain
  internal to `gbdraw.render.export`.
- The canonical CLI is `--gc` / `--no-gc`, `--skew` / `--no-skew`, and
  repeatable `--depth_track`. `RenderOutputRequest` alone owns output naming.
- Fresh CLI/API entry points reject the retired `depth_tick_interval`,
  `feature_table`, and `collinear_max_gene_gap` names; Circular
  `multi_record_size_mode=sqrt`; Linear `label_placement=on_feature` and
  `track_layout=spreadout|tuckin`; and Circular slot `spacing`, `strict`,
  `compress`, and `reserve`. Their canonical replacements are
  `depth_large_tick_interval`, `feature_visibility_table`,
  `collinear_max_unit_gap`, Circular `auto`, Linear `above_feature`,
  Linear `above|below`, and explicit `inner_gap_px` / `outer_gap_px`.
- Under `O3.data=A`, persisted-data compatibility is intentionally separate
  from executable compatibility. New writers emit session version 37 and
  canonical request schema 4. Readers retain versions 27–33 and 36 and schemas
  1–3, migrate those stored legacy names, values, and slot fields before replay,
  and ignore the obsolete nested assembly prefix. `--annotation-table` and
  `--gc_content_tick_interval` remain active aliases and are not part of the
  removal set.
- Circular private compatibility transport is confined to canonical request
  schemas 1–3 and old-session readers. The schema 4 writer never emits
  `__gbdraw_legacy_spacing`. Pixel-based legacy spacing is rewritten as explicit
  `innerGapPx` and `outerGapPx`; factor-based spacing remains replayable but
  cannot be saved losslessly as schema 4 until it is replaced with explicit
  inner and outer pixel gaps.
- The approved O4 naming and collision rules are active: an explicit prefix
  preserves dots, and duplicate implicit Circular batch IDs receive the first
  available deterministic `_2`, `_3`, ... suffix. A separate-diagram Circular
  batch that requests session output is rejected before diagram output until a
  typed batch request can represent that grouping.

The full A1/A2 and O4 planner work and the mode-specific typed-boundary
migration remain for the next API iteration under `O6.delivery=A`. That
iteration includes explicit grid/batch request forms, the remaining O4
single-layout and neutral-loading work, and removal of the legacy
`CollinearityParameters` type in favor of `LosslessCollinearityParameters`.
`DiagramOptions` and the internal assembler implementations remain because
current typed requests still require them; they are not compatibility aliases.

[Architecture/API audit](./LINEAR_CIRCULAR_ARCHITECTURE_API_AUDIT.md) | [Python API](./PYTHON_API.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md)
