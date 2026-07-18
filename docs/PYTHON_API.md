[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Python API** | [Export](./EXPORT.md)

# Draw genome diagrams from Python

The supported Python workflow has four steps:

1. Read or create Biopython `SeqRecord` values.
2. Call `draw_circular` or `draw_linear`.
3. Inspect the returned `Diagram` in memory.
4. Call `Diagram.save()` or `Diagram.to_bytes()`.

Import the main interface directly from `gbdraw`. The `gbdraw.api` package contains
lower-level integration and session components used by the CLI and web app; a normal
drawing workflow does not need them.

## Circular quick start

`read_genbank` always returns a list because one GenBank file can contain multiple
records. Pass either one record or a record sequence to `draw_circular`.

```python
import os
from pathlib import Path

from gbdraw import (
    CircularOptions,
    FeatureOptions,
    draw_circular,
    read_genbank,
)

input_path = Path(os.environ.get("GBDRAW_EXAMPLE_GBK", "MjeNMV.gb"))
output_dir = Path(os.environ.get("GBDRAW_API_OUTPUT_DIR", "."))
examples_dir = Path(os.environ.get("GBDRAW_EXAMPLES_DIR", input_path.parent))
test_inputs_dir = Path(os.environ.get("GBDRAW_TEST_INPUTS_DIR", input_path.parent))

record = read_genbank(input_path)[0]
options = CircularOptions(
    features=FeatureOptions(types=("CDS",)),
    species="Example genome",
)
diagram = draw_circular(record, options=options)
output_path = diagram.save(output_dir / "api_circular.svg", overwrite=True)

assert output_path == output_dir / "api_circular.svg"
assert diagram.to_svg().startswith("<svg")
```

The resulting circular diagram:

![Circular MjeNMV genome diagram produced by the documented Python API example](../examples/python-api-circular.png)

## The same function handles multiple circular records

Pass a sequence to create a shared circular grid. `CircularLayout` is optional for
multiple records; omit it to use the automatic layout. A layout passed with only one
record raises `ValidationError` instead of being ignored.

```python
from gbdraw import CircularLayout

multi_records = read_genbank(
    [examples_dir / "MjeNMV.gb", examples_dir / "MelaMJNV.gb"]
)
multi_diagram = draw_circular(
    multi_records,
    options=CircularOptions(features=FeatureOptions(types=("CDS",))),
    layout=CircularLayout(
        size="auto",
        positions=("#1@1", "#2@1"),
    ),
)
assert multi_diagram.to_svg().startswith("<svg")
```

`size` accepts `auto`, `equal`, `linear`, or `sqrt`. Grid-only values remain in
`CircularLayout`; they are not mixed into drawing and feature options.

## Linear diagrams and comparisons

`draw_linear` also accepts one record or a sequence. Put Linear-only comparison
settings in `LinearComparisonOptions`; Circular options cannot be passed accidentally.

```python
from gbdraw import (
    LinearComparisonOptions,
    LinearOptions,
    Thresholds,
    draw_linear,
)

linear_records = read_genbank(
    [examples_dir / "MjeNMV.gb", examples_dir / "MelaMJNV.gb"]
)
linear_diagram = draw_linear(
    linear_records,
    options=LinearOptions(
        features=FeatureOptions(types=("CDS",)),
        comparisons=LinearComparisonOptions(
            blast_files=(str(examples_dir / "MjeNMV.MelaMJNV.tblastx.out"),),
        ),
        thresholds=Thresholds(identity=0, bitscore=0),
    ),
)
assert linear_diagram.to_svg().startswith("<svg")
```

The resulting pairwise comparison:

![Linear comparison of two majanivirus records produced by the documented Python API example](../examples/python-api-linear.png)

Use `LinearLayout` only when records need explicit multi-row placement:

```python
from gbdraw import LinearLayout

linear_grid = draw_linear(
    linear_records,
    layout=LinearLayout(record_gap=28, positions=("#1@1", "#2@2")),
)
assert linear_grid.to_svg().startswith("<svg")
```

For precomputed matches, set `comparisons` to `LinearComparison` values or set
`protein_comparisons` to DataFrames. In-process protein comparison settings are also
grouped in `LinearComparisonOptions`.

## Depth tracks

Each `DepthTrackOptions` contains the data and styling for one logical track. For one
record, `source` is a path or DataFrame. For multiple records, pass one source per
record. This replaces the singular, plural, and nested depth arguments of the
lower-level rendering engine.

```python
from gbdraw import DepthTrackOptions

depth_diagram = draw_circular(
    record,
    options=CircularOptions(
        features=FeatureOptions(types=("CDS",)),
        depth_tracks=(
            DepthTrackOptions(
                source=test_inputs_dir / "MjeNMV.DRR271272.depth.tsv",
                label="Coverage",
                color="#2563eb",
            ),
        ),
    ),
)
assert depth_diagram.to_svg().startswith("<svg")
```

Set `height` on a Linear depth track. Circular diagrams reject a non-default height
because circular track width is controlled by the circular track layout.

## Conservation rings

`ConservationTrackOptions` binds each BLAST file or DataFrame to its label and color.
All tracks in one diagram must use the same source kind.

```python
from gbdraw import ConservationOptions, ConservationTrackOptions

conservation_record = read_genbank(test_inputs_dir / "AP027078.gb")[0]
conservation_diagram = draw_circular(
    conservation_record,
    options=CircularOptions(
        conservation=ConservationOptions(
            reference="query",
            tracks=(
                ConservationTrackOptions(
                    source=test_inputs_dir / "AP027078_AP027131.tblastx.out",
                    label="AP027131",
                ),
            ),
        ),
    ),
)
assert conservation_diagram.to_svg().startswith("<svg")
```

## Feature colors, visibility, and labels

Table inputs accept either a path or a validated DataFrame. One field represents one
source, so callers do not need separate `_table` and `_file` arguments.

```python
from gbdraw import LabelOptions

styled_diagram = draw_circular(
    record,
    options=CircularOptions(
        features=FeatureOptions(
            types=("CDS",),
            palette="default",
        ),
        labels=LabelOptions(
            whitelist=test_inputs_dir / "NC_010162.whitelist.tsv",
            qualifier_priority=test_inputs_dir / "NC_010162.qualifier_priority.tsv",
            overrides=test_inputs_dir / "label_override.tsv",
        ),
    ),
)
assert styled_diagram.to_svg().startswith("<svg")
```

Use `FeatureOptions.color_table`, `default_colors`, `visibility`, and `shapes` for
the corresponding feature controls.

## GFF3 and FASTA

`read_gff` accepts one paired input or equally sized path sequences.

```python
from gbdraw import read_gff

gff_records = read_gff(
    examples_dir / "gff3_lambda/lambda_two_contigs.gff3",
    examples_dir / "gff3_lambda/lambda_two_contigs.fna",
    features=("CDS", "gene"),
)
assert [item.id for item in gff_records] == ["lambda_left", "lambda_right"]
```

## In-memory and file output

`Diagram` hides the underlying SVG implementation and provides three output methods:

- `to_svg()` returns SVG text.
- `to_svg(interactive=True)` returns an interactive SVG with feature metadata.
- `to_bytes(format)` returns SVG, interactive SVG, PNG, PDF, EPS, or PS bytes.
- `save(path)` infers the format and writes exactly the requested file.

```python
svg_bytes = diagram.to_bytes("svg")
interactive_path = diagram.save(
    output_dir / "api_circular.interactive.svg",
    overwrite=True,
)

assert svg_bytes.startswith(b"<svg")
assert interactive_path.name == "api_circular.interactive.svg"
```

PNG, PDF, EPS, and PS require the optional CairoSVG dependency. Pass
`overwrite=True` to replace an existing file. Unknown filename extensions require an
explicit `format` argument.

## Errors and lower-level integration

Catch `gbdraw.exceptions.GbdrawError` for expected gbdraw failures and
`ValidationError` for invalid records or options. Mode-specific option classes reject
Circular/Linear mix-ups before rendering.

The CLI, web app, saved sessions, and integrations that need materialized input
descriptions use request and session models under `gbdraw.api`. Those models are an
orchestration layer, not a prerequisite for drawing from Python. Internal assembler
functions and `svgwrite.Drawing` are not part of the beginner-facing contract.

Pin a gbdraw version in reproducible pipelines and test representative output after
upgrading. SVG geometry can change intentionally even when the Python call remains
valid.

[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Python API**
