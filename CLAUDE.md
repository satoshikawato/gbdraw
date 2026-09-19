# Project guidance

Read [AGENTS.md](AGENTS.md) for shared working-tree, branch, verification, and
handoff rules. This file covers core architecture and documentation ownership;
use only sections relevant to the task. Web-specific guidance is in
[gbdraw/web/CLAUDE.md](gbdraw/web/CLAUDE.md).

## Core entry points

- CLI: `gbdraw.cli:main()` dispatches to `circular` and `linear`.
- Beginner Python API: package-root `read_genbank()`, `read_gff()`,
  `draw_circular()`, `draw_linear()`, mode-specific options, and `Diagram`.
- Typed integration API: `gbdraw.api` owns request, render, session, table,
  option, and track-slot contracts and explicit render helpers.
- Low-level assemblers in `gbdraw.api.diagram` and canvas/drawing configurators
  remain internal; do not re-export them from `gbdraw.api`.

Keep the API, CLI, Web UI, and both diagram modes convergent on the typed
core/planner. Adapters translate surface-specific input and output. An abstraction
should unify real paths and remove the superseded paths in the same change.
Extend existing boundaries before creating a parallel pipeline. Measure material
performance changes, including repeated I/O and computation in shared layers.

Use the [architecture ratchet](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
for architecture-bearing changes, with its ordinary evidence path unless an
exception applies.

## Finding the implementation

| Concern | Owner |
| --- | --- |
| Command-line surfaces | `gbdraw/cli.py`, `circular.py`, `linear.py` |
| Public contracts | `gbdraw/api/` |
| Input loading and tables | `gbdraw/io/` |
| Configuration | `gbdraw/config/`, `gbdraw/configurators/` |
| Diagram planning and assembly | `gbdraw/diagrams/circular/`, `gbdraw/diagrams/linear/` |
| Features, labels, legends, canvas | Matching subfolders under `gbdraw/` |
| SVG primitives and grouped output | `gbdraw/render/drawers/`, `gbdraw/render/groups/` |
| Export compatibility | `gbdraw/render/export.py` |
| Defaults, palettes, fonts | `gbdraw/data/` |
| Package version, dependencies, test markers | `pyproject.toml` |
| CI commands and supported matrix | `.github/workflows/` |

Inputs are GenBank or GFF3 with FASTA. SVG is the native output; other export
formats use optional CairoSVG. Keep type hints and existing typed configuration
patterns; use module loggers instead of print-based diagnostics.

```bash
pip install -e ".[dev]"
gbdraw circular --gbk genome.gb -o output
gbdraw linear --gbk genome1.gb genome2.gb -b blast.txt -o comparison
gbdraw gui
```

Test inputs and runners are defined in `tests/conftest.py`. Common test/build
commands and generated-reference rules have one owner in `AGENTS.md`.

## Persisted-format compatibility

- Add a compatibility reader or migrator only with evidence that the old contract
  existed in the first-parent history of `main` or a release tag, plus a positive
  representative fixture. Track session, request, cache, metadata, and other
  schema namespaces separately.
- Keep the current writer format even before it reaches `main`. If an active
  branch advances it again, rewrite branch-owned artifacts to the newest format
  before merge and remove the superseded reader, migrator, fixture, test, and
  user documentation. Do not chain migrations through or advertise branch-only
  intermediate versions.

## Public documentation ownership

`docs/DOCS.md` owns navigation. Tutorials teach a deliberate progression to a
finished result; Technical documentation owns exact behavior and contracts; FAQ
owns concise decisions and troubleshooting; Gallery helps readers discover
finished outcomes. `docs/TUTORIALS/README.md` is the tutorial index.

Use the fewest pages that answer distinct reader questions. Before adding a
public page, record the question, existing owner, and `keep`, `merge`, `delete`,
or `new` disposition. Create a page only if an existing owner cannot answer the
question clearly. Separate GUI, CLI, or Python evidence does not by itself
justify separate pages; several scenarios can support one page or no public page.

Keep public and maintainer prose concrete and technically accurate. Preserve exact
scientific terms, UI labels, CLI options, identifiers, and necessary qualifications.
Implementation plans should retain the detail needed to execute them correctly.
Use the documentation workflow for procedural evidence. A Gallery text correction
uses the Gallery skill's lightweight content route without triggering capture work.
An unrelated code review or ordinary prose typo needs neither workflow.
