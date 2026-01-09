# gbdraw Refactoring Migration Guide

This document explains all structural changes made during the refactoring project for developers familiar with the pre-refactoring codebase (the `main` branch as of late 2025).

---

## Table of Contents

1. [Overview](#overview)
2. [Package Structure Changes](#package-structure-changes)
3. [Import Path Changes](#import-path-changes)
4. [New Packages and Modules](#new-packages-and-modules)
5. [Backwards Compatibility](#backwards-compatibility)
6. [Detailed Change Log](#detailed-change-log)

---

## Overview

The refactoring project reorganized `gbdraw` to:

- **Improve cohesion**: Related modules are now grouped together
- **Enable library usage**: Clean API layer for programmatic access
- **Reduce code duplication**: Shared utilities extracted
- **Mirror structure**: Circular and linear diagrams now share identical organization
- **Prepare for multi-track**: Foundation for future Circos-like features

### Key Principles

1. **Zero regression**: All 55 tests pass; output is identical
2. **Backwards compatible**: Old import paths continue to work via re-export shims
3. **Incremental**: Changes can be adopted gradually

---

## Package Structure Changes

### Before vs After

```
BEFORE (main branch)                    AFTER (refactoring branch)
========================                ============================

gbdraw/                                 gbdraw/
├── drawers/                            ├── render/
│   ├── circular/                       │   ├── export.py
│   └── linear/                         │   ├── drawers/          ← MOVED
├── groups/                             │   │   ├── circular/
│   ├── circular/                       │   │   └── linear/
│   └── linear/                         │   └── groups/           ← MOVED
├── render/                             │       ├── circular/
│   └── export.py                       │       └── linear/
│                                       │
├── circular_diagram_components.py      ├── diagrams/
├── linear_diagram_components.py        │   ├── circular/         ← NEW
│                                       │   │   ├── assemble.py
│                                       │   │   ├── builders.py
│                                       │   │   └── positioning.py
│                                       │   └── linear/           ← NEW
│                                       │       ├── assemble.py
│                                       │       ├── builders.py
│                                       │       ├── positioning.py
│                                       │       └── precalc.py
│                                       │
│                                       ├── api/                  ← NEW
│                                       │   ├── diagram.py
│                                       │   └── tracks.py
│                                       │
│                                       ├── cli_utils/            ← NEW
│                                       │   └── arguments.py
│                                       │
│                                       ├── tracks/               ← NEW
│                                       │   ├── spec.py
│                                       │   └── parser.py
│                                       │
├── circular.py                         ├── circular.py           (unchanged)
├── linear.py                           ├── linear.py             (unchanged)
├── cli.py                              ├── cli.py                (unchanged)
│                                       │
│                                       ├── circular_diagram_components.py  ← SHIM
│                                       ├── linear_diagram_components.py    ← SHIM
└── ...                                 └── ...
```

---

## Import Path Changes

### Drawers

| Old Import | New Import |
|------------|------------|
| `from gbdraw.drawers.circular import ...` | `from gbdraw.render.drawers.circular import ...` |
| `from gbdraw.drawers.linear import ...` | `from gbdraw.render.drawers.linear import ...` |

**Example:**
```python
# OLD
from gbdraw.drawers.circular.features import CircularFeatureDrawer

# NEW
from gbdraw.render.drawers.circular.features import CircularFeatureDrawer
```

### Groups

| Old Import | New Import |
|------------|------------|
| `from gbdraw.groups.circular import ...` | `from gbdraw.render.groups.circular import ...` |
| `from gbdraw.groups.linear import ...` | `from gbdraw.render.groups.linear import ...` |

**Example:**
```python
# OLD
from gbdraw.groups.circular import SeqRecordGroup, GcContentGroup

# NEW
from gbdraw.render.groups.circular import SeqRecordGroup, GcContentGroup
```

### Diagram Assembly

| Old Import | New Import |
|------------|------------|
| `from gbdraw.circular_diagram_components import plot_circular_diagram` | `from gbdraw.diagrams.circular import plot_circular_diagram` |
| `from gbdraw.circular_diagram_components import assemble_circular_diagram` | `from gbdraw.diagrams.circular import assemble_circular_diagram` |
| `from gbdraw.linear_diagram_components import plot_linear_diagram` | `from gbdraw.diagrams.linear import plot_linear_diagram` |
| `from gbdraw.linear_diagram_components import assemble_linear_diagram` | `from gbdraw.diagrams.linear import assemble_linear_diagram` |

**Example:**
```python
# OLD
from gbdraw.circular_diagram_components import (
    assemble_circular_diagram,
    add_record_on_circular_canvas,
    center_group_on_canvas,
)

# NEW (preferred)
from gbdraw.diagrams.circular import assemble_circular_diagram
from gbdraw.diagrams.circular.assemble import add_record_on_circular_canvas
from gbdraw.diagrams.circular.positioning import center_group_on_canvas

# NEW (also works - via shim)
from gbdraw.circular_diagram_components import (
    assemble_circular_diagram,
    add_record_on_circular_canvas,
    center_group_on_canvas,
)
```

### New API Layer

For programmatic usage without CLI:

```python
# NEW - Recommended for library usage
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)

# Returns svgwrite.Drawing without saving to disk
drawing = assemble_circular_diagram_from_record(
    record=my_seq_record,
    config_dict=config,
    ...
)
```

---

## New Packages and Modules

### `gbdraw/api/` - Public API Layer

**Purpose**: Clean entry points for programmatic (non-CLI) usage.

| Module | Description |
|--------|-------------|
| `api/diagram.py` | `assemble_*_diagram_from_*()` functions that return `Drawing` objects |
| `api/tracks.py` | Track specification utilities |

**Usage:**
```python
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.render.export import save_figure

# Build diagram programmatically
drawing = assemble_circular_diagram_from_record(record, config_dict=cfg, ...)

# Save when ready
save_figure(drawing, ["svg", "png"])
```

### `gbdraw/diagrams/` - Diagram Assembly

**Purpose**: Implementation of diagram composition (placing groups on canvases).

| Module | Description |
|--------|-------------|
| `diagrams/circular/assemble.py` | Main `assemble_circular_diagram()`, `plot_circular_diagram()` |
| `diagrams/circular/builders.py` | `add_*_group_on_canvas()` helpers |
| `diagrams/circular/positioning.py` | `center_group_on_canvas()`, `place_legend_on_canvas()` |
| `diagrams/linear/assemble.py` | Main `assemble_linear_diagram()`, `plot_linear_diagram()` |
| `diagrams/linear/builders.py` | `add_*_group()` helpers |
| `diagrams/linear/positioning.py` | Position calculation utilities |
| `diagrams/linear/precalc.py` | Pre-calculation for labels and definitions |

### `gbdraw/tracks/` - Multi-Track Specification

**Purpose**: Data models for future multi-track support (Circos-like).

| Module | Description |
|--------|-------------|
| `tracks/spec.py` | `TrackSpec`, `ScalarSpec`, placement models |
| `tracks/parser.py` | `parse_track_spec()`, `parse_track_specs()` |

**Usage (experimental):**
```python
from gbdraw.tracks import parse_track_specs

# Parse track specifications
specs = parse_track_specs(["gc_content@r=0.8,w=20", "gc_skew@show=false"], mode="circular")
```

### `gbdraw/cli_utils/` - Shared CLI Utilities

**Purpose**: Common argument definitions shared between `circular.py` and `linear.py`.

| Module | Description |
|--------|-------------|
| `cli_utils/arguments.py` | Argument group factories for argparse |

### `gbdraw/render/` - Consolidated Rendering

**Purpose**: All SVG rendering code in one place.

| Subpackage | Description |
|------------|-------------|
| `render/export.py` | `save_figure()`, `parse_formats()` |
| `render/drawers/` | Low-level SVG element builders |
| `render/groups/` | High-level SVG group assemblers |

---

## Backwards Compatibility

### Shim Modules

The following modules are now **shims** that re-export from new locations:

| Shim Module | Re-exports From |
|-------------|-----------------|
| `gbdraw/circular_diagram_components.py` | `gbdraw.diagrams.circular.*` |
| `gbdraw/linear_diagram_components.py` | `gbdraw.diagrams.linear.*` |

**What this means:**
- Old import paths **continue to work**
- No immediate code changes required
- Recommend updating imports when convenient

### Deprecated Patterns

While these still work, prefer the new patterns:

```python
# DEPRECATED (still works via shim)
from gbdraw.circular_diagram_components import plot_circular_diagram

# PREFERRED
from gbdraw.diagrams.circular import plot_circular_diagram
```

### Breaking Changes

**None.** All public interfaces remain compatible.

---

## Detailed Change Log

### Phase 1-3: Foundation (Prior Sessions)

- Created `gbdraw/api/diagram.py` with `assemble_*` functions
- Split `linear_diagram_components.py` into `diagrams/linear/*`
- Added `gbdraw/tracks/` for multi-track specification
- Fixed division-by-zero bugs in GC skew calculations

### Phase 4: Render Consolidation

**Moved:**
- `gbdraw/drawers/` → `gbdraw/render/drawers/`
- `gbdraw/groups/` → `gbdraw/render/groups/`

**Updated:**
- All relative imports in `render/drawers/` and `render/groups/` (3-4 dot notation)
- `render/__init__.py` exports `save_figure`, `parse_formats`

### Phase 5: CLI Consolidation

**Created:**
- `gbdraw/cli_utils/__init__.py`
- `gbdraw/cli_utils/arguments.py`

**Updated:**
- `circular.py` and `linear.py` use shared argument utilities

### Phase 6: Circular Diagram Consolidation

**Created:**
- `gbdraw/diagrams/circular/__init__.py`
- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/builders.py`
- `gbdraw/diagrams/circular/positioning.py`

**Updated:**
- `circular_diagram_components.py` → re-export shim
- `diagrams/__init__.py` exports `circular` and `linear`

### Phase 7: Cleanup

**Removed dead code:**
- Duplicate `print_version()` in `cli.py`
- Unused `logger` instances in 5 files
- Unused `Any` import in `linear.py`

**Reviewed:**
- 184 `# type: ignore` comments (all appropriate)

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                        CLI Layer                            │
│  ┌──────────┐  ┌──────────┐  ┌──────────────────────────┐  │
│  │circular.py│  │linear.py │  │cli_utils/arguments.py   │  │
│  └─────┬────┘  └────┬─────┘  └──────────────────────────┘  │
└────────┼────────────┼───────────────────────────────────────┘
         │            │
         ▼            ▼
┌─────────────────────────────────────────────────────────────┐
│                       API Layer                             │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ api/diagram.py                                        │  │
│  │   assemble_circular_diagram_from_record()            │  │
│  │   assemble_linear_diagram_from_records()             │  │
│  └──────────────────────────────────────────────────────┘  │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│                    Diagram Assembly                         │
│  ┌────────────────────┐    ┌────────────────────┐          │
│  │ diagrams/circular/ │    │ diagrams/linear/   │          │
│  │   assemble.py      │    │   assemble.py      │          │
│  │   builders.py      │    │   builders.py      │          │
│  │   positioning.py   │    │   positioning.py   │          │
│  └─────────┬──────────┘    │   precalc.py       │          │
│            │               └─────────┬──────────┘          │
└────────────┼─────────────────────────┼──────────────────────┘
             │                         │
             ▼                         ▼
┌─────────────────────────────────────────────────────────────┐
│                      Render Layer                           │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ render/                                               │  │
│  │   export.py (save_figure, parse_formats)             │  │
│  │   drawers/circular/  drawers/linear/                 │  │
│  │   groups/circular/   groups/linear/                  │  │
│  └──────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────────────┐
│                     Core Modules                            │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐      │
│  │ canvas/  │ │ config/  │ │features/ │ │ labels/  │ ...  │
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘      │
└─────────────────────────────────────────────────────────────┘
```

---

## Quick Reference Card

### Old → New Import Mappings

```python
# Drawers
- gbdraw.drawers.circular     → gbdraw.render.drawers.circular
- gbdraw.drawers.linear       → gbdraw.render.drawers.linear

# Groups
- gbdraw.groups.circular      → gbdraw.render.groups.circular
- gbdraw.groups.linear        → gbdraw.render.groups.linear

# Diagram assembly
- gbdraw.circular_diagram_components  → gbdraw.diagrams.circular
- gbdraw.linear_diagram_components    → gbdraw.diagrams.linear

# New locations
- gbdraw.api.diagram          (NEW - public API)
- gbdraw.tracks               (NEW - track specs)
- gbdraw.cli_utils            (NEW - CLI utilities)
```

### Key Functions by Location

| Function | New Location |
|----------|--------------|
| `assemble_circular_diagram()` | `gbdraw.diagrams.circular.assemble` |
| `plot_circular_diagram()` | `gbdraw.diagrams.circular.assemble` |
| `assemble_linear_diagram()` | `gbdraw.diagrams.linear.assemble` |
| `plot_linear_diagram()` | `gbdraw.diagrams.linear.assemble` |
| `save_figure()` | `gbdraw.render.export` |
| `parse_formats()` | `gbdraw.render.export` |
| `center_group_on_canvas()` | `gbdraw.diagrams.circular.positioning` |
| `SeqRecordGroup` | `gbdraw.render.groups.circular` |
| `GcContentGroup` | `gbdraw.render.groups.circular` |

---

## Questions?

If you encounter issues with the refactored structure, please:
1. Check this guide for the new import paths
2. Review the [Refactoring Logs](./REFACTORING_LOG_2026-01-09.md)
3. Open an issue on GitHub
