[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Python API** | [Export](./EXPORT.md)

# Build diagrams with the Python API

Use `gbdraw.api` as the supported library entry point. Imports from internal modules may change without the compatibility guarantees of this namespace.

## Circular example

The example reads one GenBank record, builds a circular diagram, writes SVG, and also keeps SVG bytes in memory. Set `GBDRAW_EXAMPLE_GBK` when the input is not named `MjeNMV.gb`; the documentation test uses this variable to run the block against the repository fixture.

```python
import os
from pathlib import Path

from gbdraw.api import (
    DiagramOptions,
    OutputOptions,
    build_circular_diagram,
    load_gbks,
    render_to_bytes,
    save_figure_to,
)

input_path = Path(os.environ.get("GBDRAW_EXAMPLE_GBK", "MjeNMV.gb"))
output_dir = Path(os.environ.get("GBDRAW_API_OUTPUT_DIR", "."))

record = load_gbks([str(input_path)], mode="circular")[0]
options = DiagramOptions(
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "repeat_region"],
    species="Example genome",
    output=OutputOptions(output_prefix="api_circular", legend="right"),
)
canvas = build_circular_diagram(record, options=options)
paths = save_figure_to(
    canvas,
    "svg",
    output_dir=str(output_dir),
    output_prefix="api_circular",
    overwrite=True,
)
svg_bytes = render_to_bytes(canvas, "svg")

assert paths == [str(output_dir / "api_circular.svg")]
assert svg_bytes.startswith(b"<svg")
```

`build_linear_diagram(records, options=...)` accepts a sequence of records. Load them with `load_gbks(paths, mode="linear")`, then set comparison inputs or protein-search options on `DiagramOptions` only when the workflow needs them.

## Output choices

- `render_to_bytes(canvas, "svg")` returns SVG without writing a file.
- `save_figure_to` writes to an explicit directory/prefix and refuses to overwrite by default.
- `png`, `pdf`, `eps`, and `ps` require the optional CairoSVG dependency.
- `interactive_svg` creates the normal `.svg` plus `.interactive.svg`; supply an interactive context when rich popup metadata is required.

## Errors and stability

Catch `gbdraw.exceptions.GbdrawError` for expected gbdraw failures and `ValidationError` for invalid inputs or options. Treat exported names in `gbdraw.api` as the stable-ish public surface; option fields may grow additively, while internal assembly modules are not API documentation.

Pin a gbdraw version in reproducible pipelines and test representative output after upgrading. SVG geometry can change intentionally even when the Python call remains valid.

[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Python API**
