# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@AGENTS.md

`AGENTS.md`, imported above, holds the rules every agent follows: branching,
pull-request base and contents, when not to auto-merge, workspace cleanup, and
browser checks.

## Project Overview

**gbdraw** is a Python bioinformatics tool for creating publication-quality genome diagrams from microbial genomes and organelles. It generates circular and linear visualizations of genomic features (CDS, tRNA, regulatory elements, etc.) with optional GC content/skew plots and BLAST comparison tracks.

- **Version:** 0.14.5
- **Python:** ≥3.10
- **Output formats:** SVG, PNG, PDF, EPS, PS

## Quick Commands

- Install, test, lint, and reference-output commands:
  [Development setup](CONTRIBUTING.md#development-setup),
  [Running the checks CI runs](CONTRIBUTING.md#running-the-checks-ci-runs), and
  [Changing rendered SVG output](CONTRIBUTING.md#changing-rendered-svg-output)
  in `CONTRIBUTING.md`.
- Browser wheel, hosted bundle, and build commands:
  [Local build and verification](gbdraw/web/CLAUDE.md#local-build-and-verification)
  in `gbdraw/web/CLAUDE.md`.

```bash
# Run with coverage
pytest tests/ --cov=gbdraw --cov-report=html

# Type-check gbdraw/ with the pinned mypy and compare with the type-debt baseline
pytest tests/test_type_check_ratchet.py
mypy  # full error list

# Basic usage examples
gbdraw circular --gbk genome.gb -o output
gbdraw linear --gbk genome1.gb genome2.gb -b blast.txt -o comparison
gbdraw gui  # Launch web UI
```

Browser and Playwright availability checks are in `AGENTS.md`.

## Agent skills and subagents

- Project skills live in `.claude/skills/` and subagent definitions in
  `.claude/agents/`. A skill describes a procedure; selecting one never widens
  the task's scope or adds surfaces to it.
- Run test suites whose output is long (Playwright, the full pytest suite)
  through `.claude/scripts/run_quiet.py`, which keeps the full log in a file
  and prints a short summary with the exit status preserved.

## Compact instructions

When compacting, keep: the task goal and its done criteria; open PR numbers
with head SHAs and CI state; branch and worktree paths; Owner decisions and
Owner-delegated choices not yet written to a file; bug IDs assigned in this
session; and the path of the task's `STATE.md`. Drop file contents, command
output, and passing test logs that can be read again.

## Project Structure

- `gbdraw/`: the package. Its units are listed below by layer, from the
  bottom; a module imports from its own layer and the layers below it ("Python
  structure" in the architecture ratchet). The architecture description is
  `docs/internal/ARCHITECTURE.md`.
- `gbdraw/web/`: the Web app, a single-page app with no build step
  (`index.html`, ES modules under `js/`); see `gbdraw/web/CLAUDE.md`.
- `tests/` (pytest; `reference_outputs/` holds the reference SVGs, `web/` the
  Web tests), `docs/`, `examples/`, `tools/`, and `pyproject.toml`.

<!-- BEGIN GENERATED units: python tools/generate_architecture_docs.py -->
```text
gbdraw/
  foundation
    contract_shapes.py         Neutral descriptions of persisted JSON contracts.
    core/                      Low-level helpers: colors, fonts, text metrics, numbers, figure lengths, sequences (internal).
    data/                      Built-in data: config.toml, color palettes, color tables and fonts.
    exceptions.py              Custom exceptions for gbdraw.
    version.py                 The gbdraw version string.
  model and geometry
    analysis/                  Sequence analyses: GC content and skew, depth coverage, conservation, collinearity (internal).
    annotations/               Region annotations for circular and linear diagrams.
    auto_sizes.py              Length- and mode-dependent Auto sizes.
    canvas/                    Canvas sizes and offsets of Circular and Linear diagrams (internal).
    comparisons/               Comparison searches: LOSAT and BLAST+ runtimes, job plans, the identity scale (internal).
    config/                    The typed configuration: config.toml loading, dataclass models, overrides (internal).
    definition_line_styles.py  Parsing and normalization helpers for linear definition line styles.
    features/                  Feature objects from records: identities, colors, shapes, visibility, lane placement (internal).
    io/                        Input readers for records (GenBank, GFF3, FASTA), BLAST results and TSV tables.
    layout/                    Geometry before drawing: record coordinates and rows, Linear wrap, figure size, composition (internal).
    linear_comparison.py       Normalized comparison contracts for Linear diagrams.
    losat_setup.py             Explicit installation and offline resolution of a pinned LOSAT release.
    mode_profiles.py           Versioned defaults and shared validation for Circular and Linear modes.
    svg/                       SVG path data and elements for features, tracks, arcs and ticks (internal).
    tracks/                    Circular track-slot input helpers.
  placement
    configurators/             Drawing parameters derived from the configuration, per track and for the legend (internal).
    labels/                    Feature label placement, fit and filtering for Circular and Linear diagrams (internal).
    legend/                    Legend tables, entry metrics and layout for both modes (internal).
  drawing
    diagrams/                  Diagram assembly internals.
    render/                    SVG groups and drawers, interactive SVG metadata, and export to other formats (internal).
    render_stages.py           Render stages: function calls whose results a Worker transaction may reuse.
  API and Session
    api/                       Typed library contracts for gbdraw requests, sessions, tables, and rendering.
    interface.py               Small, typed, beginner-facing Python API for gbdraw.
    session.py                 Public canonical session document lifecycle.
    session_contract.py        The Session document of the current version, declared once.
    session_drawings.py        Where each drawing of a Session lives: the one owner of the document layout.
    session_io.py              GUI session JSON loading, validation, materialization, and sidecar building.
    session_legacy_values.py   Values that older Sessions and schema 1-2 requests wrote, and what they read as now.
    session_migration.py       Migrations that bring a saved Session's Web-owned fields to the current writer.
    session_request_codec.py   Internal codec for canonical, CLI-independent render request payloads.
    session_resources.py       The resources of one Session document, allocated by content.
  entry points and adapters
    __init__.py                Public Python interface for gbdraw.
    _build_support.py          Build-time helpers: package data, the browser wheel and the Pyodide asset config.
    _web_assets.py             The vendored Web assets: their files, source URLs, versions and notices.
    circular.py                The gbdraw circular command: its arguments, Session replay and the Circular request.
    cli.py                     The gbdraw command: it dispatches to the circular, linear, render, setup-losat and gui subcommands.
    cli_utils/                 CLI module with shared argument definitions and utilities.
    crop_genbank.py            Crop a GenBank record to a region and shift its features; also runs as a script.
    linear.py                  The gbdraw linear command: its arguments, Session replay and the Linear request.
    web_support/               Support helpers shared by the browser UI and offline web asset tooling.
```
<!-- END GENERATED units -->

## Key Architecture

### Entry Points

1. **CLI:** `gbdraw.cli:main()` → dispatches to `circular` or `linear` subcommands
2. **Beginner-facing Python API:** the `gbdraw` package root exports
   `read_genbank()`, `read_gff()`, `draw_circular()`, `draw_linear()`,
   mode-specific option types, and `Diagram`.
3. **Typed integration API:** `gbdraw.api` exports typed request, render, session,
   table, option, and track-slot contracts plus explicit render helpers.

Low-level assemblers remain in `gbdraw.api.diagram` for internal engine use.
Canvas and drawing configurators remain implementation details in their owner
modules; they are not re-exported from `gbdraw.api`.

Architecture changes should keep these entry points convergent:

- Route the Python API, CLI, Web UI, and both diagram modes through the same
  typed core/planner; adapters should only translate surface-specific input and output.
- As a default, add an abstraction only when it unifies at least two real
  execution paths and removes the superseded paths in the same change. Extend
  existing boundaries instead of creating parallel pipelines.
- Prefer a few stable, composable contracts over mode- or surface-specific
  branches. Add capabilities as data at an existing boundary when practical,
  so later extension or removal does not multiply change points.
- Avoid repeated I/O or computation in shared layers, and protect material
  performance changes with measurements.

Architecture-bearing changes must follow the
[architecture fitness-function ratchet](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md).
An ordinary change removes every superseded owner or path in the same change
([Merge acceptance](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#merge-acceptance));
the evidence each change class records is in
[Required author evidence](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#required-author-evidence).

### Persisted-format compatibility

- Compatibility readers and migrators require evidence that the old contract
  existed in the first-parent history of `main` or in a release tag, plus a
  representative positive fixture. Track session, request, cache, metadata, and
  other schema namespaces separately.
- Keep the current writer format even before it reaches `main`. If an active branch
  advances that format again, rewrite branch-owned artifacts to the newest format
  before merge and remove the superseded reader, migrator, fixture, test, and user
  documentation. Do not chain migrations through, or advertise, intermediate
  versions that never reached `main`.
- The `SESSION_VERSION` rule for a breaking Session change (the unit is
  `main`) and the co-changes that a change to the generated contracts in
  `gbdraw/web/js/contracts/` needs:
  [Persisted contract co-change](docs/internal/WEB_CHANGE_POLICY.md#persisted-contract-co-change).

### Working-tree overwrite policy

- Uncommitted state is not a preservation requirement in this repository. When an
  in-scope fix conflicts with dirty code, tests, documentation, fixtures, or
  generated artifacts, replace the uncommitted implementation as needed to make
  the result correct and internally consistent.
- Inspect the existing diff before editing so the replacement is deliberate.
  Keep files unrelated to the requested work out of scope; this policy does not
  authorize broad cleanup or deletion elsewhere in the working tree.
- Treat Gallery sessions, source SVGs, rendered examples, thumbnails, and
  `examples.json` as generator-owned outputs. A requested Gallery refresh may
  overwrite dirty copies directly; regenerate them from the current code and
  declared inputs instead of merging or preserving their previous bytes.

### Data Flow

1. **Input:** GenBank/GFF3+FASTA files
2. **Loading:** `gbdraw.io.genome.load_gbks()`
3. **Configuration:** TOML config → dataclass models (`GbdrawConfig`)
4. **Assembly:** `gbdraw.diagrams` combines features, GC plots, labels
5. **Rendering:** `gbdraw.render` generates SVG via `svgwrite`
6. **Export:** `gbdraw.render.export` saves to various formats (CairoSVG optional)

### Configurators Pattern

Main configurator classes encapsulate drawing logic:
- `FeatureDrawingConfigurator` - Genomic features
- `GcContentConfigurator` - GC content track
- `GcSkewConfigurator` - GC skew track
- `BlastMatchConfigurator` - BLAST comparison (linear only)
- `LegendDrawingConfigurator` - Color legend

### Render Module Structure

The `gbdraw/render/` module has a two-tier architecture:
- **drawers/**: Low-level SVG element builders (individual shapes, paths)
- **groups/**: High-level SVG group assemblers that compose drawers
- **export.py**: Internal format parsing and deprecated lenient export compatibility

### Web UI Structure

- The entry point, the ES modules, and their owners:
  [Non-negotiable properties](gbdraw/web/CLAUDE.md#non-negotiable-properties),
  [Runtime data flow](gbdraw/web/CLAUDE.md#runtime-data-flow), and
  [Module ownership](gbdraw/web/CLAUDE.md#module-ownership) in `gbdraw/web/CLAUDE.md`.

## Coding Conventions

### Style
- **Type hints:** Extensively used with `from __future__ import annotations`
- **Naming:** snake_case (functions/modules), PascalCase (classes), UPPER_SNAKE_CASE (constants)
- **Logging:** Use `logging` module (`logger = logging.getLogger(__name__)`)

### Patterns
- Frozen dataclasses for configuration models
- Factory methods (`.from_dict()`) for config parsing
- `NamedTuple` for immutable data structures

### Type checking
- The type-debt rule for `gbdraw/` (mypy, `TYPE_DEBT_BASELINE` in
  `tests/test_type_check_ratchet.py`):
  [Type checking](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#type-checking)
  in the architecture ratchet.
- Class-body annotations are runtime data. `get_type_hints` reads dataclass
  fields in `gbdraw/config/modify.py` (config overrides) and
  `gbdraw/session_request_codec.py` (canonical request validation), and
  `tests/test_public_contract.py` hashes the signatures of `gbdraw.__all__`.
  Treat a change to either as a code change, and use only names importable at
  run time in them, not names imported under `TYPE_CHECKING`.
- Plan: `docs/internal/PYTHON_TYPE_CHECK_RATCHET_PLAN_2026-10-06.md`.

## Key Configuration Files

### gbdraw/data/config.toml
Default settings for canvas dimensions, track types, feature styling, etc.

### gbdraw/data/color_palettes.toml
Predefined color palettes for feature visualization.

### pyproject.toml
Package configuration, test markers, coverage settings.

## Testing

### Test Markers
```python
@pytest.mark.slow         # Skip in fast runs
@pytest.mark.regression   # Regression tests
@pytest.mark.circular     # Circular diagram tests
@pytest.mark.linear       # Linear diagram tests
```

### Reference Output Tests
[Changing rendered SVG output](CONTRIBUTING.md#changing-rendered-svg-output) in `CONTRIBUTING.md`.

### Test Helpers (tests/conftest.py)
- `GbdrawRunner.run()` - Run either diagram subcommand with optional BLAST inputs
- `find_test_input` fixture - Find test inputs across repository and optional local directories
- `tmp_path` - Built-in pytest fixture for isolated test outputs

## CI/CD

- **Tiers, Python versions, and slow tests:** [Validation tiers](docs/internal/SELECTIVE_CI.md#validation-tiers) in `docs/internal/SELECTIVE_CI.md`
- **Lint job:** Installs the Ruff version pinned in the `dev` extra of `pyproject.toml` and blocks CI on lint failures
- **Type check:** `Core PR` and `Core` install the `typecheck` extra and run `tests/test_type_check_ratchet.py` (mypy 2.4.1) with the other core tests
- **CairoSVG:** The Core jobs install the `export` extra (`cairosvg`) with `typecheck`; no workflow installs extra system packages

## Dependencies

### Core
- **BioPython** - Genome file parsing
- **svgwrite** - SVG generation
- **pandas** - Data manipulation
- **fonttools** - Font metrics
- **bcbio-gff** - GFF3 parsing

### Optional
- **cairosvg** - PNG/PDF/EPS/PS export (requires `libcairo2-dev libpango1.0-dev` on Ubuntu)

## Important Notes

1. **Python 3.10+ required** due to `from __future__ import annotations`
2. **CairoSVG is optional** - only needed for non-SVG export formats
3. **Track types:** "spreadout", "middle", "tuckin" control circular diagram layout
4. **Genome size thresholds:** Window/step sizes auto-adjust (<1M, 1-10M, >10M bp)
5. **Label filtering:** Supports priority files, blacklists, whitelists
6. **BLAST comparison:** Linear diagrams only, requires outfmt 6/7
7. **Browser wheel workflow:** [Non-negotiable properties](gbdraw/web/CLAUDE.md#non-negotiable-properties) in `gbdraw/web/CLAUDE.md`.

## Updating Reference Outputs

When intentional changes affect SVG output, follow
[Changing rendered SVG output](CONTRIBUTING.md#changing-rendered-svg-output) in `CONTRIBUTING.md`.

## Documentation

- The documentation rules (four teaching routes, when to add a page, prose,
  and terms): [Changing documentation](CONTRIBUTING.md#changing-documentation)
  in `CONTRIBUTING.md`.
- Main docs: `docs/DOCS.md`
- Canonical tutorial index: `docs/TUTORIALS/README.md`
- Current documentation plan:
  `docs/internal/DOCUMENTATION_RESTRUCTURE_PLAN_2026-10-09.md`, which keeps
  the routes and evidence architecture of
  `docs/internal/DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md`
- CLI Reference: `docs/CLI_Reference.md`
- **Web app development:** See `gbdraw/web/CLAUDE.md` for web-specific guidance
- Web app: https://gbdraw.app/
- GitHub: https://github.com/satoshikawato/gbdraw
