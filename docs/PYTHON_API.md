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
multiple records; omit it to use the automatic layout. You can also pass a layout
with one record to produce an explicit 1×1 grid.

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

`size` accepts `auto`, `equal`, or `linear`. Grid-only values remain in
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

## Sequence-similarity comparison rings

`ComparisonRingTrackOptions` binds each BLAST or LOSAT result path, or a DataFrame
with the same tabular hits, to its label and color. The resulting ring displays
sequence-similarity hits. It does not by itself establish biological conservation.
All rings in one diagram must use the same source kind.

For interactive SVG output, set `comparison_sequence_source` to a FASTA path, one
`SeqRecord`, or a sequence of records to enable the comparison-span FASTA actions
for that ring. Static geometry does not depend on this optional source.

```python
from gbdraw import ComparisonRingOptions, ComparisonRingTrackOptions

comparison_record = read_genbank(test_inputs_dir / "AP027078.gb")[0]
comparison_diagram = draw_circular(
    comparison_record,
    options=CircularOptions(
        comparison_rings=ComparisonRingOptions(
            reference="query",
            tracks=(
                ComparisonRingTrackOptions(
                    source=test_inputs_dir / "AP027078_AP027131.tblastx.out",
                    label="AP027131",
                ),
            ),
        ),
    ),
)
assert comparison_diagram.to_svg().startswith("<svg")
```

`CircularOptions.comparison_rings` is the canonical field.
`CircularOptions.conservation` remains a runtime constructor and attribute alias for
existing code. Both names address the same stored option; use the canonical name in
new code and when calling `dataclasses.replace()`. Passing non-`None` values for both
names is rejected.
The lower-level request transport retains its
`conservation_*` field names. The older `ConservationOptions` and
`ConservationTrackOptions` class names are identity aliases for
`ComparisonRingOptions` and `ComparisonRingTrackOptions`.

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

### Canonical label configuration overrides

`config_overrides` accepts canonical dotted leaf paths. Label configuration has
three separate concerns:

| Path | Accepted values | Meaning |
|---|---|---|
| `labels.circular.scope` | `none`, `outer`, `both` | Whether Circular labels are hidden, outer-side only, or allowed on both outer and inner sides |
| `labels.circular.placement` | `horizontal`, `radial` | Orientation of external Circular labels |
| `labels.linear.scope` | `none`, `all`, `first`, `orthogroup_top` | Which Linear records or similarity-group members are eligible for labels (`orthogroup_top` is the compatibility token) |
| `labels.linear.placement` | `auto`, `above_feature` | Linear label geometry |
| `labels.linear.rotation` | float | Linear label rotation in degrees |
| `labels.rendering` | `auto`, `embedded_only`, `external_only` | Shared policy for embedding labels in feature bodies or routing them externally |

For example, Circular labels can use both outer and inner sides with radial
external text:

```python
circular_label_diagram = draw_circular(
    record,
    options=CircularOptions(
        features=FeatureOptions(types=("CDS",)),
        config_overrides={
            "labels.circular.scope": "both",
            "labels.circular.placement": "radial",
            "labels.rendering": "auto",
        },
    ),
)
assert circular_label_diagram.to_svg().startswith("<svg")
```

The Linear equivalents use their own scope vocabulary and also expose
placement and rotation:

```python
linear_label_diagram = draw_linear(
    linear_records,
    options=LinearOptions(
        features=FeatureOptions(types=("CDS",)),
        config_overrides={
            "labels.linear.scope": "first",
            "labels.linear.placement": "above_feature",
            "labels.linear.rotation": 45.0,
            "labels.rendering": "auto",
        },
    ),
)
assert linear_label_diagram.to_svg().startswith("<svg")
```

`scope` selects eligible records or Circular label sides, `placement` controls
text geometry, and only `labels.rendering` is the embedded/external rendering
policy. Fresh Python requests reject the retired flat `show_labels` and
`allow_inner_labels` aliases and the old `canvas.show_labels`,
`canvas.circular.show_labels`, `canvas.linear.show_labels`, and
`canvas.circular.allow_inner_labels` paths. Supported persisted-data readers
migrate those keys; they are not canonical override examples for new code.

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
The former `plot_circular_diagram` / `plot_linear_diagram` save wrappers and
package-level `gbdraw.render` aliases have been removed; use `Diagram.save`,
`render_to_bytes`, or `save_figure_to` at their documented boundaries.

Typed Circular requests accept `CircularDiagramOptions`, with
`CircularRequestTrackOptions` and `CircularOutputOptions` for mode-specific
nested settings. Typed Linear requests use `LinearDiagramOptions`,
`LinearRequestTrackOptions`, and `LinearOutputOptions`. Import these integration
types from `gbdraw.api`.

Typed requests use `DepthTrackInput` for depth data. Add one entry to
`CircularDiagramOptions.depth_tracks` or `LinearDiagramOptions.depth_tracks` for
each logical track. `source` accepts one path or DataFrame shared by all displayed
records, or a sequence containing one path, DataFrame, or `None` per record.
`LinearDiagramOptions` also accepts `height` on each entry.

```python
from gbdraw.api import DepthTrackInput, LinearDiagramOptions

typed_options = LinearDiagramOptions(
    depth_tracks=(
        DepthTrackInput(
            source=(
                test_inputs_dir / "MjeNMV.DRR271272.depth.tsv",
                test_inputs_dir / "MjeNMV.DRR271272.depth.tsv",
            ),
            label="Coverage",
            color="#2563eb",
            height=24,
        ),
    ),
)
```

The older `depth_table`, `depth_file`, `depth_tables`, `depth_files`, and
`depth_track_*` fields remain compatibility inputs. Do not combine them with
`depth_tracks`; request normalization rejects mixed forms. New request and session
writers serialize either input style as canonical `depthTracks`.

The shorter `gbdraw.api.CircularTrackOptions` and
`gbdraw.api.LinearTrackOptions` names remain compatibility aliases. They refer to
the request classes above. By contrast, `gbdraw.CircularTrackOptions` and
`gbdraw.LinearTrackOptions` are separate beginner-facing classes used inside the
package-root `CircularOptions` and `LinearOptions`.

`CircularDiagramRequest` explicitly selects `single` or `grid`; a one-record
grid is valid.
`CircularBatchRequest` selects `batch` and carries one
`RenderOutputRequest` per record, without grid placement.
Optional Circular comparison FASTA paths belong to
`CircularDiagramOptions.conservation_fasta_files`; this keeps interactive
comparison-span metadata inside the same typed request and session resource
contract as the corresponding BLAST inputs.
`CircularRequestPlan`, `CircularBatchRequestPlan`, and `LinearRequestPlan`
normalize records, layout, and builder selection. `CircularBatchRenderResult`
returns the corresponding per-record results. The root facade, fresh CLI and
Web generation, current canonical replay, and legacy internal replay all render
through these planners while assembly remains output-neutral.

`RenderOutputRequest.output_prefix` preserves an explicit prefix exactly, including
dots: `sample.v1` produces `sample.v1.svg`. In the Circular CLI, duplicate implicit
record IDs in a separate-diagram batch receive deterministic suffixes. A batch
session preserves its explicit `batch` grouping and resolved output array.

Canonical request schema 5 never writes the private
`__gbdraw_legacy_spacing` migration transport. Schema 1 and 2 readers can
replay factor-based Circular spacing, but the current writer cannot re-save it
losslessly. Replace it with explicit `inner_gap_px` and `outer_gap_px` values
before saving a current session.

Canonical request schema 5 persists explicit `single`, `grid`, or `batch`
grouping. Its output is one object for a single diagram or grid and an array
for Circular batch, with one entry per record. Record loading is mode-neutral;
planners own topology warnings and mode, comparison, and cardinality policy.
Active and public runtime collinearity configuration uses
`LosslessCollinearityParameters`; canonical request schemas 1 and 2 privately
migrate legacy `standard` parameter payloads while preserving their effective
fields. Current schema 5 accepts only the lossless form.

Pin a gbdraw version in reproducible pipelines and test representative output after
upgrading. SVG geometry can change intentionally even when the Python call remains
valid.

[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Python API**
