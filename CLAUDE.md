# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**gbdraw** is a Python bioinformatics tool for creating publication-quality genome diagrams from microbial genomes and organelles. It generates circular and linear visualizations of genomic features (CDS, tRNA, regulatory elements, etc.) with optional GC content/skew plots and BLAST comparison tracks.

- **Version:** 0.8.3
- **Python:** ≥3.10
- **Output formats:** SVG, PNG, PDF, EPS, PS

## Quick Commands

```bash
# Run tests (fast, excludes slow tests)
pytest tests/ -v -m "not slow"

# Run all tests including slow
pytest tests/ -v

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
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503

# Install in development mode
pip install -e ".[dev]"

# Build distribution
python -m build

# Basic usage examples
gbdraw circular --gbk genome.gb -o output
gbdraw linear --gbk genome1.gb genome2.gb -b blast.txt -o comparison
gbdraw gui  # Launch web UI
```

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
│   └── web/                   # Web app assets
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
2. **Public API:** `gbdraw.api` module exports:
   - `assemble_circular_diagram_from_record()` - Build circular diagram from SeqRecord
   - `assemble_linear_diagram_from_records()` - Build linear diagram from multiple SeqRecords
   - Configurator classes (`FeatureDrawingConfigurator`, `GcContentConfigurator`, etc.)
   - Canvas configurators (`CircularCanvasConfigurator`, `LinearCanvasConfigurator`)
   - Track specs (`TrackSpec`, `parse_track_specs()`) - Control individual track visibility/styling

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
- **export.py**: Format conversion (`save_figure()`, `parse_formats()`)

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
- `GbdrawRunner` - Helper class for running gbdraw commands in tests
- `get_test_input_path()` - Find test input files across directories
- `get_reference_output_path()` - Find reference SVG files
- `temp_output_dir` fixture - Creates temporary directory for test outputs

## CI/CD

- **Python versions tested:** 3.10, 3.11, 3.12
- **Lint job:** Uses ruff for code formatting checks
- **CairoSVG job:** Separate test with export dependencies on Python 3.11
- **Slow tests:** Only run on push to main branch

## Dependencies

### Core
- **BioPython** - Genome file parsing
- **svgwrite** - SVG generation
- **pandas** - Data manipulation
- **fonttools** - Font metrics
- **bcbio-gff** - GFF3 parsing

### Optional
- **cairosvg** - PNG/PDF/EPS/PS export

## Important Notes

1. **Python 3.10+ required** due to `from __future__ import annotations`
2. **CairoSVG is optional** - only needed for non-SVG export formats
3. **Track types:** "spreadout", "middle", "tuckin" control circular diagram layout
4. **Genome size thresholds:** Window/step sizes auto-adjust (<1M, 1-10M, >10M bp)
5. **Label filtering:** Supports priority files, blacklists, whitelists
6. **BLAST comparison:** Linear diagrams only, requires outfmt 6/7

## Updating Reference Outputs

When intentional changes affect SVG output, reference files need updating:

```bash
# 1. Run tests to identify failures
pytest tests/test_output_comparison.py -v

# 2. Review differences - ensure changes are intentional
# 3. Regenerate and copy updated outputs to tests/reference_outputs/
# 4. Commit with explanation of visual changes
```

## Documentation

- Main docs: `docs/DOCS.md`
- Tutorials: `docs/TUTORIALS/`
- CLI Reference: `docs/CLI_Reference.md`
- **Web app development:** See `gbdraw/web/CLAUDE.md` for web-specific guidance
- Web app: https://gbdraw.app/
- GitHub: https://github.com/satoshikawato/gbdraw
