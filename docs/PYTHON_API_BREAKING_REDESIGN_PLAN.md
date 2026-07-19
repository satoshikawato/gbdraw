# Python API Breaking Redesign Plan

- Status: proposed
- Date: 2026-07-19
- Baseline: `0.14.0b0`
- Target: the next intentionally breaking beta before the 1.0 public contract
- Scope: Python library API, CLI/Python boundary, Web/Pyodide boundary, typed render requests, session conversion, public documentation, and contract tests

## 1. Decision

Replace the current layered compatibility API with one small Python drawing API and
one explicit integration API.

The final public boundary will be:

```text
gbdraw                 ordinary Python drawing workflows
gbdraw.integration     typed request, session, and table integrations
gbdraw.errors          structured library exceptions
```

Rendering, layout, parsing, and SVG implementation modules remain implementation
details. They may stay in their current physical packages when moving them would add
churn without improving the public dependency boundary.

This plan intentionally breaks existing Python imports and signatures. It does not
preserve the following as public contracts:

- `gbdraw.api`
- `DiagramOptions`
- `ColorOptions`, `TrackOptions`, and `OutputOptions`
- `build_circular_diagram`, `build_circular_multi_diagram`, and
  `build_linear_diagram`
- `assemble_circular_diagram_from_record`,
  `assemble_circular_diagram_from_records`, and
  `assemble_linear_diagram_from_records`
- public configurator and canvas-configurator classes
- raw layout strings such as `"#1@1"` in the Python API
- unstructured `config_overrides` in the ordinary drawing API

Temporary adapters may exist while the implementation branch is in progress. They
must not remain in the final release unless a separately approved use case requires
them.

This document is a standalone replacement plan. Unfinished earlier Python API plans
are not prerequisites and do not constrain this redesign. Current code, current
tests, supported session documents, and verified rendering behavior are the sources
of truth.

## 2. Why the current boundary should be replaced

The current repository has several overlapping abstractions:

1. `gbdraw` exports a 20-symbol beginner-facing interface.
2. `gbdraw.api` exports a much larger integration and rendering surface.
3. `gbdraw.interface` translates `CircularOptions` and `LinearOptions` back into the
   shared `DiagramOptions` representation.
4. `build_*` functions expand `DiagramOptions` into long `assemble_*` calls.
5. CLI handlers still call the long assemblers directly for important workflows.
6. Typed request and session models also contain the shared `DiagramOptions` rather
   than the mode-specific beginner options.
7. Web/Pyodide execution still has an argument-array execution path even though
   session documents contain a typed render-request payload.

The result is a small outer API over a large compatibility core. Adding a new option
often requires forwarding it through several models and adapters, while users still
see two namespaces that both appear public.

The redesign makes the mode-specific typed request the authoritative input to
rendering. Every other entry point converts to that request exactly once.

## 3. Goals

1. Make the first successful diagram require no knowledge of internal config,
   renderer classes, `svgwrite`, or output prefixes.
2. Give every drawing setting one clear owner and one typed representation.
3. Use separate Circular and Linear option types so wrong-mode settings cannot be
   represented.
4. Keep file parsing, validation, rendering, and export failures inside one
   documented exception hierarchy.
5. Remove string mini-languages from the Python object model.
6. Make CLI, Web, session replay, and direct Python calls render through the same
   request and validation path.
7. Keep `Diagram` as the only ordinary result type and keep file output separate
   from drawing intent.
8. Preserve current default SVG geometry unless a geometry change is separately
   proposed and reviewed.
9. Publish usable type information in built wheels.
10. Reduce code and forwarding paths after the new request path is complete.

## 4. Non-goals

- Do not rewrite low-level SVG drawers solely to rename the public API.
- Do not change CLI option names as part of the Python API redesign.
- Do not change visible Web controls or post-render editing behavior unless required
  to remove the argument-array Python bridge.
- Do not make browser-only editing actions part of the Python drawing API.
- Do not introduce a plugin framework, general request DSL, or fluent builder API.
- Do not add tutorial-only functions such as `load_example()` to the production
  package.
- Do not regenerate reference SVGs for an API-only refactor.
- Do not drop readable version 31 or 32 session files merely because the Python
  object model changes.

## 5. Target user experience

### 5.1 Minimal Circular diagram

```python
from gbdraw import draw_circular, read_genbank

records = read_genbank("genome.gb")
diagram = draw_circular(records)
diagram.save("genome.svg")
```

`read_genbank()` continues to return all records. `draw_circular()` accepts one
`SeqRecord` or a sequence, so the quick start does not discard additional records
with `[0]`.

### 5.2 Customized Circular diagram

```python
from gbdraw import CircularOptions, draw_circular, read_genbank
from gbdraw.options import FeatureOptions, GcOptions, TitleOptions

options = CircularOptions(
    features=FeatureOptions(types=("CDS", "tRNA", "rRNA")),
    gc=GcOptions(window=1_000, step=100),
    title=TitleOptions(text="Example genome", position="top"),
)

diagram = draw_circular(read_genbank("genome.gb"), options=options)
diagram.save("genome.svg")
```

### 5.3 Linear comparison

```python
from gbdraw import LinearOptions, draw_linear, read_genbank
from gbdraw.options import ComparisonOptions, Thresholds

records = read_genbank(("reference.gb", "comparison.gb"))
options = LinearOptions(
    comparison=ComparisonOptions.from_blast(
        "reference.comparison.tblastx.tsv",
        thresholds=Thresholds(identity=70),
    )
)

diagram = draw_linear(records, options=options)
diagram.save("comparison.svg")
```

`from_blast()` is justified because it replaces several mutually exclusive source
fields. Convenience constructors must not be added for simple field assignment.

### 5.4 Typed placement

```python
from gbdraw import CircularLayout, RecordPlacement

layout = CircularLayout(
    placements=(
        RecordPlacement(record=0, row=0, column=0),
        RecordPlacement(record=1, row=0, column=1),
    )
)
```

Python indices are zero-based. CLI and TSV parsers may keep their existing user
syntax but must convert it at the parsing boundary.

### 5.5 Error handling

```python
from gbdraw import GbdrawError

try:
    diagram.save("genome.pdf")
except GbdrawError as exc:
    print(f"gbdraw failed: {exc}")
```

## 6. Public module contract

### 6.1 `gbdraw`

The package root is the ordinary workflow facade. Keep its `__all__` intentionally
small.

Required root exports:

- `draw_circular`
- `draw_linear`
- `read_genbank`
- `read_gff`
- `Diagram`
- `CircularOptions`
- `LinearOptions`
- `CircularLayout`
- `LinearLayout`
- `RecordPlacement`
- `GbdrawError`
- `InputError`
- `OptionError`
- `DependencyError`
- `RenderError`
- `ExportError`
- `__version__`

Specialized option types live in `gbdraw.options`. Frequently used types may be
re-exported from the root only when a real quick-start workflow becomes shorter.
Their defining module remains `gbdraw.options`.

### 6.2 `gbdraw.options`

This module owns the public immutable drawing configuration:

- `CircularOptions`
- `LinearOptions`
- `FeatureOptions`
- `LabelOptions`
- `TitleOptions`
- `LegendOptions`
- `GcOptions`
- `DepthTrack`
- `CircularTrackOptions`
- `LinearTrackOptions`
- `ConservationOptions`
- `ConservationTrack`
- `ComparisonOptions`
- `ComparisonTrack`
- `ProteinComparisonOptions`
- `CollinearityOptions`
- `Thresholds`
- `CircularLayout`
- `LinearLayout`
- `RecordPlacement`

The implementation may combine or remove a proposed type if the field inventory
shows that it does not reduce caller and implementation complexity. The acceptance
criterion is ownership clarity, not the exact number of dataclasses.

### 6.3 `gbdraw.integration`

This is the only advanced public namespace. It owns typed orchestration rather than
rendering primitives:

- `CircularRequest`
- `LinearRequest`
- `DiagramRequest`
- `RecordInput`
- `GenBankSource`
- `GffFastaSource`
- `InMemorySource`
- `RecordSelector`
- `Region`
- `OutputPlan`
- `RenderResult`
- `render`
- `render_to_files`
- validated table models and readers
- session document load, materialization, conversion, render, and save functions
- session-specific exception types

Do not export canvas configurators, SVG groups, assembler functions, parser
internals, or `svgwrite.Drawing` from this namespace.

### 6.4 `gbdraw.errors`

Use one structured hierarchy:

```text
GbdrawError
├── InputError
│   ├── InputNotFoundError
│   └── InputFormatError
├── OptionError
├── DependencyError
├── RenderError
├── ExportError
└── SessionError
    ├── SessionFormatError
    ├── SessionVersionError
    └── SessionResourceError
```

Each expected error should expose:

- `code`: stable machine-readable identifier
- `path`: optional option or input location such as
  `comparison.thresholds.identity`
- `hint`: optional recovery guidance

Preserve the original exception with exception chaining. Do not introduce a custom
`Result[T, E]` abstraction.

## 7. Option ownership rules

The shared 71-field `DiagramOptions` is deleted. Replace it with mode-specific
composition.

### 7.1 Common ownership

| Concern | Owner |
|---|---|
| feature selection, colors, visibility, shapes | `FeatureOptions` |
| label filtering and overrides | `LabelOptions` |
| plot title | `TitleOptions` |
| legend placement and style | `LegendOptions` |
| GC content/skew calculation and visibility | `GcOptions` |
| one quantitative depth series and its sampling/style | `DepthTrack` |
| annotation sets and track binding | a typed annotation option under `gbdraw.options` |

### 7.2 Circular ownership

| Concern | Owner |
|---|---|
| Circular track order and axis placement | `CircularTrackOptions` |
| conservation sources, ring geometry, and thresholds | `ConservationOptions` |
| species/strain definition policy | a Circular definition option block |
| multi-record sizing and placement | `CircularLayout` |

### 7.3 Linear ownership

| Concern | Owner |
|---|---|
| Linear track order and axis placement | `LinearTrackOptions` |
| nucleotide comparison sources and thresholds | `ComparisonOptions` |
| protein comparison execution | `ProteinComparisonOptions` |
| collinearity parameters | `CollinearityOptions` |
| row gap and record placement | `LinearLayout` |

Thresholds belong to the feature that consumes them. `LinearOptions` must not carry
a detached shared `thresholds` field.

Sampling fields belong to `GcOptions` or `DepthTrack`. `CircularOptions` and
`LinearOptions` must not expose detached `window`, `step`, `depth_window`, or
`depth_step` fields.

### 7.4 Advanced config

Remove `config_overrides: Mapping[str, object]` from the ordinary API. It bypasses
type checking and makes internal TOML paths part of the public contract.

For settings that are not yet modeled:

1. Add them to an existing owner when they are genuinely user-facing.
2. Add one focused option block when doing so reduces complexity.
3. Keep renderer-only settings private.
4. Permit a typed `GbdrawConfig` only in `gbdraw.integration` for advanced
   integration during the transition.

Do not replace one untyped mapping with another differently named mapping.

## 8. Input and table boundary

### 8.1 Ordinary record input

`read_genbank()` and `read_gff()` remain explicit format readers. They return all
parsed records and raise `InputError` subclasses for expected file and parse
failures.

`draw_circular()` and `draw_linear()` accept:

- one `SeqRecord`
- a non-empty sequence of `SeqRecord`

They do not accept arbitrary file paths. Keeping reading explicit avoids ambiguous
format detection and keeps GFF3/FASTA pairing clear.

### 8.2 Integration record input

`RecordInput` combines exactly one source with optional selection, region, stable
key, label, reverse-complement flag, and typed placement.

Prefer named constructors for mutually exclusive sources:

```python
RecordInput.from_genbank(path, selector=...)
RecordInput.from_gff(gff_path, fasta_path, selector=...)
RecordInput.from_record(record)
```

The source subclasses may remain public for serialization and type narrowing, but
ordinary callers should not need to construct them directly.

### 8.3 Tables

Paths and DataFrames are accepted only at public input boundaries. Normalize them
immediately into capability-specific validated models such as:

- `FeatureColorTable`
- `FeatureVisibilityTable`
- `LabelWhitelist`
- `QualifierPriority`
- `LabelOverrides`
- `DepthTable`
- `ConservationTable`
- `ComparisonTable`

Each model provides `from_tsv()` and `from_dataframe()`. Validation, required
columns, coordinate conventions, and path-dependency resolution have one owner.

Do not propagate parallel `*_table` and `*_file` fields into requests or renderers.

## 9. Remove string mini-languages from Python

Raw strings remain valid only at CLI, TSV, and old-session decode boundaries.

Replace Python-facing layout strings with `RecordPlacement`.

Replace Python-facing track-slot strings with typed values containing:

- renderer enum
- side enum
- optional width or height
- optional gap
- visibility
- axis relationship where applicable

Existing parser functions may remain internal adapters for CLI and old session
documents. Renderers receive only typed, validated slot objects.

Selectors and regions follow the same rule: parse text once at the input boundary,
then pass typed `RecordSelector` and `Region` values.

## 10. Authoritative render request

### 10.1 Request types

```python
@dataclass(frozen=True)
class CircularRequest:
    records: tuple[RecordInput, ...]
    options: CircularOptions
    layout: CircularLayout | None = None


@dataclass(frozen=True)
class LinearRequest:
    records: tuple[RecordInput, ...]
    options: LinearOptions
    layout: LinearLayout | None = None
```

`DiagramRequest` is the union of these two types. There is no shared
`DiagramOptions` and no runtime wrong-mode option filtering.

### 10.2 Output separation

Drawing intent and file output policy are separate:

```python
diagram = render(request)
result = render_to_files(request, OutputPlan(...))
```

`render(request)` returns `Diagram` and performs no filesystem writes.

`render_to_files()` returns a structured `RenderResult` containing only paths that
exist, the normalized request, records used, warnings, and format metadata.

`Diagram.save()` remains the simple one-file API. `Diagram.to_svg()` and
`Diagram.to_bytes()` remain the in-memory APIs.

### 10.3 Resolution stages

Use explicit internal stages:

```text
CircularRequest | LinearRequest
        ↓ validate and load inputs
ResolvedCircularRequest | ResolvedLinearRequest
        ↓ compute features, tracks, labels, and layout
RenderModel
        ↓ SVG drawers and groups
Diagram
```

The exact internal type names may differ, but each stage must have one direction.
Do not translate a new option model back into the deleted legacy model.

## 11. CLI and Web integration

### 11.1 CLI

Keep `argparse` and current option names. Change the handler boundary:

```text
argparse.Namespace
        ↓ mode-specific CLI adapter
CircularRequest | LinearRequest
        ↓ shared render path
Diagram or RenderResult
```

After migration, `gbdraw.circular` and `gbdraw.linear` must not call public or
private long-form assemblers. Each mode has one adapter that converts parsed CLI
values into the typed request.

CLI warning-and-skip policies may differ from strict library behavior, but the
request validation and rendering implementation remain shared.

### 11.2 Web/Pyodide

The Web app already builds a structured render-request payload for sessions. Use
that payload as the Pyodide execution boundary.

Replace:

```text
JavaScript form → CLI argument array → Python CLI owner → renderer
```

with:

```text
JavaScript form → versioned render-request payload
                → Python request decoder → shared render path
```

Do not expose Python dataclass implementation details to JavaScript. The JSON
payload remains versioned and is decoded into the Python request types.

Keep browser-only PNG/PDF export and post-render editing outside the Python output
API.

## 12. Session compatibility

The request representation changes incompatibly, so new sessions need the next
available render-request schema and session version. At the current baseline these
would be render-request schema 3 and session version 33; confirm the live constants
immediately before implementation and increment from those values if they changed.

Rules:

1. New writers emit only the new mode-specific option structure.
2. Python and Web readers continue accepting supported version 31 and 32 documents.
3. Old render-request schemas decode through a dedicated conversion module into the
   new request types.
4. Old documents are never converted by reconstructing CLI argument arrays.
5. Saved results and editor state remain separate from drawing intent.
6. Materialized embedded paths remain valid only inside the materialization context.
7. Remove an old session reader only through a separate compatibility decision,
   not as an incidental effect of this API redesign.

The old `DiagramOptions` payload is therefore a wire-compatibility input for old
sessions, not a runtime model used by new rendering.

## 13. Implementation workstreams

### Phase 0: Freeze evidence and approve the target contract

Tasks:

1. Record the current `gbdraw.__all__`, `gbdraw.api.__all__`, relevant signatures,
   option fields/defaults, session versions, and render-request schemas.
2. Build a field inventory mapping every current `DiagramOptions` field to:
   - a new option owner;
   - integration-only typed config;
   - internal renderer state; or
   - intentional removal with migration guidance.
3. Inventory all imports of `gbdraw.api`, `assemble_*`, `build_*`, and
   `DiagramOptions` in package code and tests.
4. Record hashes for tracked reference SVGs and run comparison-only tests.
5. Add an ADR approving the public namespace, request/output separation, and old
   session read policy.

Exit criteria:

- No current option lacks a disposition.
- The proposed `gbdraw.__all__` and `gbdraw.integration.__all__` are reviewed as
  explicit lists.
- Rendering evidence is captured before implementation.

### Phase 1: Add the new option and error models

Tasks:

1. Create `gbdraw.options` with the mode-specific immutable option model.
2. Create `gbdraw.errors` and migrate existing exception classes into the new
   hierarchy.
3. Add field-level validation at option construction time.
4. Move comparison thresholds, GC sampling, and depth sampling to their consuming
   owners.
5. Add structured `code`, `path`, and `hint` properties.
6. Add tests for defaults, invalid values, mode separation, equality, and repr.

Exit criteria:

- `CircularOptions` contains no Linear-only field and vice versa.
- No ordinary option field accepts an untyped override mapping.
- Expected invalid inputs raise only documented `GbdrawError` subclasses.

### Phase 2: Add typed inputs, tables, placements, and track slots

Tasks:

1. Implement `RecordPlacement` and typed track-slot values.
2. Implement capability-specific table models and conversion constructors.
3. Make CLI/session string parsers return the typed models.
4. Make path/DataFrame normalization occur once at the boundary.
5. Add zero-based Python placement tests and explicit CLI/session conversion tests.

Exit criteria:

- New Python examples contain no `"#1@1"`-style layout token.
- Render code receives no parallel file/DataFrame field pair.
- Required table columns and coordinate conventions are tested in one owner.

### Phase 3: Make the new requests drive rendering

Tasks:

1. Create `CircularRequest` and `LinearRequest` using the new option types.
2. Separate `OutputPlan` from both request types.
3. Extract input resolution into mode-neutral and mode-specific pure functions.
4. Refactor diagram assembly to consume resolved request fields directly.
5. Return `Diagram` from the shared in-memory render function.
6. Return `RenderResult` from the shared multi-file function.
7. Keep a temporary legacy-to-new adapter only while downstream callers migrate.

Exit criteria:

- The render path does not construct `DiagramOptions`.
- New requests do not call `build_*` or public `assemble_*` functions.
- Default Circular and Linear reference SVG comparisons are unchanged.

### Phase 4: Rewire the ordinary Python facade

Tasks:

1. Replace `gbdraw.interface` translation logic with direct construction of the new
   in-memory request.
2. Keep `draw_circular()` and `draw_linear()` signatures limited to records,
   mode-specific options, and layout.
3. Keep `Diagram` independent of `svgwrite` in its public annotations and docs.
4. Re-export the approved root symbols, including exceptions.
5. Add `py.typed` to the wheel and a packaged-consumer type-check test.
6. Rewrite the quick start without environment-variable scaffolding.

Exit criteria:

- The ordinary workflow imports only from `gbdraw` and `gbdraw.options`.
- A source-checkout and installed-wheel quick start both run unchanged.
- Static type tooling recognizes the installed package as typed.

### Phase 5: Rewire CLI and Web/Pyodide

Tasks:

1. Change Circular and Linear CLI handlers to build the new requests.
2. Remove direct assembler calls from both handlers.
3. Decode the Web render-request payload directly in Pyodide.
4. Remove the Web path that invokes Python by reconstructing CLI argument arrays for
   current sessions and normal generation.
5. Preserve legacy CLI/session replay only in isolated compatibility modules.
6. Compare CLI, direct Python, and Web/Pyodide request projections for representative
   Circular and Linear cases.

Exit criteria:

- Ordinary Web generation does not depend on CLI argument positions.
- CLI and Python render the same geometry from equivalent typed requests.
- Browser generation, interactive SVG, and session save/load checks pass.

### Phase 6: Migrate the session schema

Tasks:

1. Define the new render-request JSON schema from the new mode-specific requests.
2. Update Python encode/decode and Web encode/decode together.
3. Increment the request schema and session version.
4. Add version 31/32 conversion fixtures for both Circular and Linear modes.
5. Write only the new schema after loading an old supported document.
6. Verify embedded resource lifetime and cleanup on success and failure.

Exit criteria:

- New Python and Web sessions round-trip without field loss.
- Supported old sessions render through the new request path.
- No supported session conversion calls an argparse parser.

### Phase 7: Delete the legacy Python API

Tasks:

1. Remove `DiagramOptions`, `ColorOptions`, `TrackOptions`, and `OutputOptions`.
2. Remove public `assemble_*` and `build_*` functions.
3. Remove public configurator and canvas-configurator exports.
4. Move remaining integration symbols from `gbdraw.api` to
   `gbdraw.integration`.
5. Delete `gbdraw.api` after all package imports have migrated.
6. Delete temporary new-to-legacy and legacy-to-new adapters.
7. Remove obsolete forwarding tests and replace them with request-to-owner tests.

Exit criteria:

- `rg` finds no package import of `gbdraw.api`, `DiagramOptions`, or a removed
  assembler.
- There is one active options model and one active request model.
- Removing the compatibility modules reduces total code; the redesign must not end
  as an additional permanent layer.

### Phase 8: Documentation and release migration

Tasks:

1. Rewrite `docs/PYTHON_API.md` around the ordinary workflow.
2. Add an exhaustive generated or mechanically checked option reference.
3. Add a separate `docs/PYTHON_INTEGRATION_API.md` for requests, sessions, and
   tables.
4. Update README, docs index, workflow guide, tutorials, export docs, FAQ, release
   notes, and affected examples.
5. Add an old-to-new migration table for every removed public symbol.
6. Explain typed placement and table schemas with copy-pasteable examples.
7. Execute every documented Python block in tests.

Exit criteria:

- The first Python example is self-contained after its stated setup step.
- No documentation recommends `gbdraw.api`, `DiagramOptions`, `build_*`, or
  `assemble_*`.
- The ordinary guide and integration guide do not claim the same entry point.

## 14. Deletion and migration matrix

| Removed contract | Replacement |
|---|---|
| `gbdraw.api` | `gbdraw.integration` for orchestration; `gbdraw` for drawing |
| `DiagramOptions` | `CircularOptions` or `LinearOptions` |
| `ColorOptions` | fields owned by `FeatureOptions` |
| `TrackOptions` | `CircularTrackOptions` or `LinearTrackOptions` |
| `OutputOptions` | drawing fields in `TitleOptions`/`LegendOptions`; filesystem fields in `OutputPlan` |
| `build_circular_diagram` | `draw_circular` or `integration.render(CircularRequest(...))` |
| `build_circular_multi_diagram` | `draw_circular(..., layout=CircularLayout(...))` |
| `build_linear_diagram` | `draw_linear` or `integration.render(LinearRequest(...))` |
| three public `assemble_*` functions | private assembly called only by the shared render path |
| `CircularMultiRecordOptions` | `CircularLayout` |
| `LinearMultiRecordOptions` | `LinearLayout` |
| `RenderOutputRequest` | `OutputPlan`, separate from drawing request |
| layout strings | `RecordPlacement` |
| raw track-slot strings | typed Circular/Linear slot values |
| `config_overrides` | typed option owner or integration-only `GbdrawConfig` |
| direct `svgwrite.Drawing` result | `Diagram` |

## 15. Test strategy

### 15.1 Public contract

- Replace the current snapshot with separate snapshots for:
  - `gbdraw.__all__`;
  - `gbdraw.options.__all__`;
  - `gbdraw.integration.__all__`;
  - option/request dataclass fields and defaults;
  - public function and method signatures;
  - exception inheritance and error codes.
- Review the snapshot diff manually. Do not accept a bulk rewrite without mapping
  every removal and addition to this plan.

### 15.2 Behavior and forwarding

- Test each non-default option at the function that consumes it.
- Test that no value is accepted and then ignored.
- Test aliases only at CLI/session parse boundaries; the Python object model should
  contain normalized values.
- Test file and DataFrame forms for every table input and compare normalized content.
- Test wrong-mode construction at type/model boundaries rather than during render.

### 15.3 Cross-surface parity

For representative cases, compare the resolved request and SVG output from:

- direct Python facade;
- `gbdraw.integration` request;
- Circular CLI;
- Linear CLI;
- Web/Pyodide payload;
- new session replay;
- supported old session conversion.

Cover at least:

- single-record Circular;
- multi-record Circular with typed placement;
- Linear nucleotide comparison;
- depth tracks;
- conservation rings;
- protein comparison and collinearity;
- feature color/visibility/shape tables;
- labels and annotations;
- static and interactive SVG;
- strict local binary export.

### 15.4 Rendering regression

Run comparison-only reference tests throughout the refactor. Updating public APIs
does not justify changed SVG geometry.

If a geometry change is unavoidable:

1. isolate it from the API refactor;
2. explain why the new request changes rendering semantics;
3. review the SVG diff;
4. update references only with the explicit reference-update command;
5. record the behavior change in release notes.

### 15.5 Errors and cleanup

Test:

- missing and unreadable inputs;
- malformed GenBank, GFF3, FASTA, BLAST, depth, and table inputs;
- invalid option paths and values;
- missing optional exporters or external binaries;
- stale output files and overwrite policy;
- temporary session resource cleanup on success and exceptions;
- exception chaining and stable `code` values.

## 16. Required verification gates

Run the smallest focused suite after each phase, then complete the following gates
before deleting the legacy API:

```bash
pytest tests/test_public_contract.py -v
pytest tests/test_python_interface.py -v
pytest tests/test_api_library_usage.py -v
pytest tests/test_api_requests.py -v
pytest tests/test_api_request_render.py -v
pytest tests/test_api_session.py -v
pytest tests/test_session_io.py -v
pytest tests/test_session_request_codec.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
ruff check gbdraw/
python -m build
```

Tests will be renamed as ownership changes. Preserve the behaviors represented by
the list rather than keeping obsolete filenames.

When Web execution changes, also run:

```bash
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/services/config.js
node --check gbdraw/web/js/services/session-request.js
```

Run the relevant Web unit and Playwright checks through the available Node or Python
Playwright installation. If Chromium is blocked by the agent sandbox, rerun the same
local check with sandbox escalation.

## 17. Completion criteria

The redesign is complete only when all of the following are true:

1. The documented quick start runs from a clean installation with the stated input
   setup.
2. Ordinary users need only `gbdraw` and `gbdraw.options`.
3. Advanced integrations use only `gbdraw.integration` and documented domain
   modules.
4. `gbdraw.api` and all long public assemblers are removed.
5. `DiagramOptions` is removed from runtime code and remains only in old-session
   decoding fixtures or migration documentation, if needed.
6. CLI, Web, Python, and session rendering share one typed request path.
7. New Python code uses no raw layout or track-slot string DSL.
8. File output is separate from drawing requests.
9. All expected failures inherit from root-exported `GbdrawError`.
10. The installed wheel contains `py.typed` and passes a consumer type-check test.
11. Default reference SVG comparisons pass without updates.
12. New and supported old sessions render successfully through the new request
    types.
13. The final change deletes more compatibility and forwarding code than it adds in
    permanent adapters.
14. Documentation and release notes provide a complete old-to-new migration map.

## 18. Recommended implementation order

Implement this as one coordinated breaking-release branch but keep commits and
review units phase-oriented:

1. contract ADR and field inventory;
2. new options and errors;
3. typed inputs, tables, placement, and slots;
4. request-driven rendering;
5. ordinary facade;
6. CLI adapter;
7. Web/Pyodide adapter;
8. session schema and old-session conversion;
9. legacy deletion;
10. docs, packaging, full verification, and release notes.

Do not begin by renaming `gbdraw.api`. First make the new request own rendering;
otherwise the rename only relocates the current compatibility layers. Delete the old
namespace after all production callers have moved and the shared render path has
passed geometry and session regression gates.
