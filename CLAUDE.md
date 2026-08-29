# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**gbdraw** is a Python bioinformatics tool for creating publication-quality genome diagrams from microbial genomes and organelles. It generates circular and linear visualizations of genomic features (CDS, tRNA, regulatory elements, etc.) with optional GC content/skew plots and BLAST comparison tracks.

- **Version:** 0.14.0b0
- **Python:** ≥3.10
- **Output formats:** SVG, PNG, PDF, EPS, PS

## Quick Commands

```bash
# Run tests (fast, excludes slow tests)
pytest tests/ -v -m "not slow"

# Run all tests including slow
pytest tests/ -v

# Compare generated SVGs with tracked references (read-only)
pytest tests/test_output_comparison.py::TestOutputComparison -v

# Update tracked references only for an intentional, reviewed geometry change
pytest tests/test_output_comparison.py::TestGenerateReferences --update-reference-outputs -v

# Run a single test file
pytest tests/test_regression.py -v

# Run a single test by name
pytest tests/ -v -k "test_circular_basic"

# Run tests by marker
pytest tests/ -v -m "circular"
pytest tests/ -v -m "linear"

# Run with coverage
pytest tests/ --cov=gbdraw --cov-report=html

# Check code formatting (ruff)
ruff check gbdraw/

# Install in development mode
pip install -e ".[dev]"

# Prepare the generated browser wheel for offline web packaging/tests
python tools/prepare_browser_wheel.py

# Refresh the cache-bust token when preparing a deployable web bundle
python tools/prepare_browser_wheel.py --refresh-cache-bust

# Build distribution
python -m build

# Basic usage examples
gbdraw circular --gbk genome.gb -o output
gbdraw linear --gbk genome1.gb genome2.gb -b blast.txt -o comparison
gbdraw gui  # Launch web UI
```

### Browser / Playwright Checks

- Do not treat missing repo-local `node_modules/`, `package.json`, or `@playwright/test` as proof that browser testing is unavailable. This environment may provide Playwright through Python/conda.
- Check both installations when web UI verification matters:
  - `command -v playwright && playwright --version`
  - `python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"`
- JavaScript Playwright specs in `tests/web/*.playwright.spec.js` require Node's `@playwright/test`; verify with `node -e "console.log(require.resolve('@playwright/test'))"`.
- If `@playwright/test` is missing, run an equivalent targeted browser check with Python Playwright instead of skipping browser verification.
- In Codex/agent sandboxes, Chromium can fail with `sandbox_host_linux.cc ... Operation not permitted`. Rerun the same local browser check with the required sandbox escalation before declaring Playwright unavailable.

## Project Structure

```
gbdraw/
├── gbdraw/                    # Main package
│   ├── cli.py                 # CLI entry point (gbdraw command)
│   ├── circular.py            # Circular diagram CLI handler
│   ├── linear.py              # Linear diagram CLI handler
│   ├── api/                   # Public programmatic API
│   ├── canvas/                # Canvas configuration (dimensions, layout)
│   ├── config/                # Configuration loading (TOML parsing)
│   ├── configurators/         # Feature/GC/Legend/BLAST configurators
│   ├── core/                  # Core utilities (sequences, colors, text)
│   ├── diagrams/              # Diagram assembly (circular/, linear/)
│   ├── features/              # Feature objects, coordinates, colors, tracks
│   ├── io/                    # Input/output (genome loading, color tables)
│   ├── labels/                # Label handling & filtering
│   ├── legend/                # Legend generation
│   ├── render/                # SVG rendering (drawers/, groups/, export)
│   ├── data/                  # Built-in data (config.toml, color_palettes.toml, fonts)
│   └── web/                   # Web app (SPA, no build step)
│       ├── index.html         # Entry point, templates, CSP, CSS
│       └── js/                # ES modules (app.js entry, app/ with feature subfolders, services/, utils/)
├── tests/                     # Test suite (pytest)
│   ├── reference_outputs/     # Reference SVG files for comparison
│   └── utils/                 # Test utilities
├── docs/                      # Documentation
├── examples/                  # Example files (GenBank, FASTA, BLAST)
└── pyproject.toml             # Package configuration
```

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
Ordinary non-increasing changes record concise owner/path evidence and remove
superseded paths. Complete before/after OE, PE, and CB sets plus a separate
maintainer decision are required only for the exception conditions defined by
that policy.

### Persisted-format compatibility

- Compatibility readers and migrators require evidence that the old contract
  existed in the first-parent history of `main` or in a release tag, plus a
  representative positive fixture. Track session, request, cache, metadata, and
  other schema namespaces separately.
- Keep the current writer format even before it reaches `main`. If an active branch
  advances that format again, rewrite branch-owned artifacts to the newest format
  before merge and remove the superseded reader, migrator, fixture, test, and user
  documentation. Do not chain migrations through or advertise branch-only
  intermediate versions.

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

- `gbdraw/web/index.html` is the SPA entry point (HTML/CSS/templates + CSP) and loads `gbdraw/web/js/app.js`.
- JavaScript is split into ES modules under `gbdraw/web/js/` (`state.js`, `config.js`, `components.js`, `services/`, `utils/`, `app/`). Larger UI modules are grouped under `gbdraw/web/js/app/` subfolders (for example `legend/`, `legend-layout/`, `feature-editor/`) while the top-level `app/*.js` files keep the `create*` entry points.

## Coding Conventions

### Style
- **Type hints:** Extensively used with `from __future__ import annotations`
- **Naming:** snake_case (functions/modules), PascalCase (classes), UPPER_SNAKE_CASE (constants)
- **Logging:** Use `logging` module (`logger = logging.getLogger(__name__)`)

### Patterns
- Frozen dataclasses for configuration models
- `type: ignore` comments for BioPython (missing type stubs)
- Factory methods (`.from_dict()`) for config parsing
- `NamedTuple` for immutable data structures

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
Tests compare generated SVG against `tests/reference_outputs/` files.

### Test Helpers (tests/conftest.py)
- `GbdrawRunner.run()` - Run either diagram subcommand with optional BLAST inputs
- `find_test_input` fixture - Find test inputs across repository and optional local directories
- `tmp_path` - Built-in pytest fixture for isolated test outputs

## CI/CD

- **Python versions tested:** 3.10, 3.11, 3.12
- **Lint job:** Uses Ruff 0.15.12 and blocks CI on lint failures
- **CairoSVG:** The Python matrix installs the `dev` extra and required system packages (`libcairo2-dev`, `libpango1.0-dev` on Ubuntu)
- **Slow tests:** Only run on push to main branch

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
7. **Browser wheel workflow:** `gbdraw/web/gbdraw-<version>-py3-none-any.whl` is a generated, gitignored asset. Prepare it before wheel-dependent web packaging checks or distribution builds.

## Updating Reference Outputs

When intentional changes affect SVG output, reference files need updating:

```bash
# 1. Run comparison tests to identify intentional differences. Actual output is
#    retained under pytest's temporary directory, not reference_outputs/.
pytest tests/test_output_comparison.py::TestOutputComparison -v

# 2. Review the differences and regenerate only after approving the change.
pytest tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs -v

# 3. Review the tracked SVG diff and verify the new references.
git diff -- tests/reference_outputs/
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Without `--update-reference-outputs`, reference-generation tests are skipped and
normal test runs do not write to `tests/reference_outputs/`.

## Documentation

- Public documentation has four teaching routes: Tutorials, Technical
  documentation, FAQ, and Gallery. `docs/DOCS.md` is the navigation authority.
- Use the fewest public pages that answer distinct reader questions. Coverage
  does not require one page per capability, workflow, or interface.
- Before adding a public page, record the reader question, the existing owner,
  and a `keep`, `merge`, `delete`, or `new` disposition. Add a page only when
  editing an existing owner cannot answer the question clearly.
- A Tutorial needs a deliberate learning progression to a finished result.
  Separate GUI, CLI, or Python evidence does not by itself justify separate
  public pages.
- Technical documentation owns exact behavior and contracts. FAQ owns concise
  decisions and troubleshooting, and links to the technical owner. Gallery is
  for discovering finished outcomes.
- Executable scenarios may support one public page or no public page. Do not
  turn an evidence inventory into a public-page inventory.
- Keep user- and maintainer-facing prose concise, concrete, and technically
  accurate. Optimize internal implementation plans for execution fidelity.
- Preserve exact technical terms, UI labels, CLI options, identifiers, and
  scientifically necessary qualifications.
- Main docs: `docs/DOCS.md`
- Canonical tutorial index: `docs/TUTORIALS/README.md`
- Current documentation plan:
  `docs/internal/DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md`
- CLI Reference: `docs/CLI_Reference.md`
- **Web app development:** See `gbdraw/web/CLAUDE.md` for web-specific guidance
- Web app: https://gbdraw.app/
- GitHub: https://github.com/satoshikawato/gbdraw
