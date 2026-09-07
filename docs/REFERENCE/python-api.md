[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md)

# Python API reference

The package-root API is the stable, beginner-facing interface for reading annotated sequences, drawing Circular or Linear diagrams, and returning or saving the result. Integrations that need explicit input cardinality, planning, output preflight, current analysis artifacts, or session conversion use the typed [`gbdraw.api` request interface](typed-requests.md).

## Public functions

```text
read_genbank(paths: str | PathLike[str] | Sequence[str | PathLike[str]]) -> list[SeqRecord]
read_gff(gff_paths: str | PathLike[str] | Sequence[str | PathLike[str]], fasta_paths: str | PathLike[str] | Sequence[str | PathLike[str]], *, features: Sequence[str] | None = None) -> list[SeqRecord]
draw_circular(records: RecordCollection, *, options: CircularOptions | None = None, layout: CircularLayout | None = None, record_displays: Sequence[RecordDisplayOptions] | None = None) -> Diagram
draw_linear(records: RecordCollection, *, options: LinearOptions | None = None, layout: LinearLayout | None = None, record_displays: Sequence[RecordDisplayOptions] | None = None) -> Diagram
```

`read_genbank()` returns every record from every supplied file. `read_gff()` requires equally sized GFF3 and FASTA path lists and matching sequence IDs. Each drawing function accepts one Biopython `SeqRecord` or a sequence of records. Both reject an empty record collection, non-`SeqRecord` members, and option or layout objects for the wrong mode.

The reader helpers accept filesystem paths, not open file handles. If another
library has already parsed an upload or in-memory stream, pass its Biopython
`SeqRecord` values directly to `draw_circular()` or `draw_linear()`.

## Package-root exports

The `gbdraw` package exports the four functions above, `__version__`, and these public types:

| Area | Types |
|---|---|
| Result | `Diagram` |
| Record display / placement | `RecordDisplayOptions`, `FeaturePlacementOverride`, `FeaturePlacementTarget` |
| Shared presentation | `FeatureOptions`, `LabelOptions`, `TitleOptions`, `Thresholds`, `DepthTrackOptions` |
| Circular | `CircularOptions`, `CircularLayout`, `CircularTrackOptions`, `ComparisonRingOptions`, `ComparisonRingTrackOptions` |
| Linear | `LinearOptions`, `LinearLayout`, `LinearTrackOptions`, `LinearComparisonOptions` |
| Compatibility aliases | `ConservationOptions`, `ConservationTrackOptions` |

`ConservationOptions` and `ConservationTrackOptions` are identity aliases for `ComparisonRingOptions` and `ComparisonRingTrackOptions`. New code should use the comparison-ring names.

## Shared option defaults

| Type and field | Type | Default |
|---|---|---|
| `FeatureOptions.types` | sequence of feature names or `None` | `CDS`, `rRNA`, `tRNA`, `tmRNA`, `ncRNA`, `misc_RNA`, `repeat_region` |
| `FeatureOptions.color_table` | path or `DataFrame` | `None` |
| `FeatureOptions.default_colors` | path or `DataFrame` | `None` |
| `FeatureOptions.palette` | string | `default` |
| `FeatureOptions.visibility` | path or `DataFrame` | `None` |
| `FeatureOptions.shapes` | mapping | `None` |
| `FeatureOptions.placements` | TSV path, DataFrame, exact override sequence, or `None` | `None` |
| `LabelOptions.whitelist` | path or `DataFrame` | `None` |
| `LabelOptions.qualifier_priority` | path or `DataFrame` | `None` |
| `LabelOptions.overrides` | path or `DataFrame` | `None` |
| `TitleOptions.text` | string or `None` | `None` |
| `TitleOptions.position` | `none`, `center`, `top`, `bottom`, or `None` | `None` |
| `TitleOptions.font_size` | positive number or `None` | `None` |
| `CircularOptions.legend`, `LinearOptions.legend` | legend position token | `right` |
| `CircularOptions.dinucleotide`, `LinearOptions.dinucleotide` | two-base string | `GC` |
| `window`, `step`, `depth_window`, `depth_step` | integer or `None` | `None` |
| `annotations`, `config`, `config_overrides` | mode-appropriate value or `None` | `None` |
| `depth_tracks` | sequence of `DepthTrackOptions` | empty |

Unset `Thresholds` values resolve through the mode profile:

| Mode | E-value | Bitscore | Identity | Alignment length | GC / skew on a fresh request |
|---|---:|---:|---:|---:|---|
| Circular | `1e-5` | `50` | `70` | `0` | on / on |
| Linear | `1e-2` | `50` | `0` | `0` | off / off |

E-value, bitscore, and identity must be finite and non-negative; identity is limited to `0`–`100`. Alignment length must be a non-negative integer.

## Layout and track options

`CircularLayout` defaults to `size="auto"`, `min_radius_ratio=0.55`, `column_gap_ratio=0.10`, `row_gap_ratio=0.05`, and automatic positions. `size` accepts `linear`, `auto`, or `equal`.

`LinearLayout` defaults to `record_gap=24.0` and automatic positions. Position strings identify a displayed record and its grid row, for example `#1@1`.

`CircularTrackOptions` and `LinearTrackOptions` accept an ordered `slots` sequence and a zero-based `axis_index`. Circular tracks also accept `center_reserved_radius`. An explicit slot sequence is authoritative; the axis index must agree with the selected mode and slot types.

Package-root `CircularTrackOptions.slots` contains
`gbdraw.api.CircularTrackSlot` values. This cross-namespace type is part of the
public track contract.

Each `DepthTrackOptions` represents one logical series. `source` is one path or `DataFrame` for one displayed record, or one path, `DataFrame`, or `None` per record. `label`, `color`, tick intervals, and tick font size default to `None`. `height` is supported only by Linear diagrams; Circular options reject it.

## Circular options

`CircularOptions.tracks` defaults to an empty `CircularTrackOptions`. `comparison_rings` defaults to an empty `ComparisonRingOptions`; `species` and `strain` default to `None`; `keep_full_definition_with_title` defaults to `False`. Circular title position accepts `None`, `none`, `top`, or `bottom`.

`ComparisonRingOptions` defaults to no tracks, `reference="auto"`, and automatic ring width and gap. `reference` accepts `query`, `subject`, or `auto`. Every `ComparisonRingTrackOptions` requires a BLAST/LOSAT table path or `DataFrame`; its label, color, and comparison sequence source are optional. All tracks in one diagram must use the same source kind. If any ring has a label or color, every ring must provide that field.

`CircularOptions.comparison_rings` is the canonical field. The `conservation` constructor and attribute alias remains available for compatible code, but passing both names is an error.

## Linear options

`LinearOptions.tracks` defaults to an empty `LinearTrackOptions`; `comparisons` defaults to `LinearComparisonOptions()`. Linear title position accepts `None`, `center`, `top`, or `bottom`.

Important `LinearComparisonOptions` defaults are:

| Field | Default |
|---|---|
| `protein_mode` | `none` |
| `match_style` | `ribbon` |
| `collinearity_unit` | `auto` |
| `collinearity_anchor` | `rbh` |
| `collinearity_scope` | `adjacent` |
| `collinearity_color` | `orientation` |
| `losat_executable` | `losat` |
| `blastp_executable` | `None` |
| `threads` | `None` |
| `max_hits` | `5` |
| `candidate_limit` | `None` |
| `orthogroup_membership` | `anchor_core_v1` |
| `orthogroup_member_max_hits` | `5` |
| `max_paralog_links` | `2` |

`protein_mode` accepts `none`, `pairwise`, `orthogroup`, or `collinear`. The `orthogroup` token means gbdraw Similarity groups; it does not claim phylogenetic orthology.

`LinearComparisonOptions(blast_files=...)` consumes prepared comparison TSV
files. Supplying those files does not start a nucleotide or protein search.

## Canonical label overrides

`config_overrides` uses canonical dotted leaf paths:

| Path | Accepted values |
|---|---|
| `labels.circular.scope` | `none`, `outer`, `both` |
| `labels.circular.placement` | `horizontal`, `radial` |
| `labels.linear.scope` | `none`, `all`, `first`, `orthogroup_top` |
| `labels.linear.placement` | `auto`, `above_feature` |
| `labels.linear.rotation` | finite degrees |
| `labels.rendering` | `auto`, `embedded_only`, `external_only` |

New requests use these dotted paths. Retired flat label names and
`canvas.*show_labels` paths are not accepted.

## `Diagram` output

`Diagram.mode` is `circular` or `linear`, and `Diagram.records` is the rendered record tuple. `Diagram.save()`, `Diagram.to_svg()`, and `Diagram.to_bytes()` provide four output forms:

| Method | Result |
|---|---|
| `to_svg()` | static SVG text |
| `to_svg(interactive=True)` | interactive SVG text with supported metadata |
| `to_bytes(format="svg")` | SVG, interactive SVG, PNG, PDF, EPS, or PS bytes |
| `save(path, *, format=None, overwrite=False)` | writes exactly one file and returns its `Path` |

For static SVG, `to_svg()` equals
`to_bytes("svg").decode("utf-8")`. `save(..., format="svg")` writes that same
UTF-8 payload.

PNG, PDF, EPS, and PS require CairoSVG. `save()` infers a known format from the path unless `format` is explicit and refuses to replace an existing file unless `overwrite=True`.

## Errors

Catch `gbdraw.exceptions.GbdrawError` for expected gbdraw failures and `ValidationError` for invalid records, paths, options, or mode combinations. Pin a gbdraw version and compare representative outputs after upgrading; stable Python calls do not imply byte-identical SVG geometry.

## Related

- [Draw and save your first genome diagram from Python](../TUTORIALS/PYTHON/first-genome-diagram.md)
- [Python Tutorials](../TUTORIALS/PYTHON/README.md)
- [Typed request reference](typed-requests.md)
- [Session and request compatibility](session-and-request-compatibility.md)
- [Output format and export reference](output-formats-and-export.md)

## Requested feature placement

`FeatureOptions.placements` accepts a TSV path, a pandas DataFrame, or a sequence
of `FeaturePlacementOverride` values. The package root and `gbdraw.api` export the
same `FeaturePlacementOverride` and `FeaturePlacementTarget` types. Tables are
resolved once into exact source identities by the shared request planner.

For example, `FeatureOptions(placements="placements.tsv")` can be combined with
`record_displays=[RecordDisplayOptions(start_coordinate=71)]` on either drawing
function. Set tolerance through the existing options' `config_overrides`, using
`{"canvas.feature_overlap_tolerance_bp": 1}`; booleans and negative values fail.

The table columns are `record`, `feature_selector`, `placement`, and `level`.
`main` and `auto` have no level; a supported directional target uses level 1.
Auto removes the exact override. Unknown identities, ambiguous selectors,
duplicate resolved identities and unsupported target directions fail explicitly.

## Record display start

Both drawing functions accept `record_displays`, with exactly one
`RecordDisplayOptions(is_circular=None, start_coordinate=None)` per input record.
`None` adds no shift; explicit `1` anchors source base 1, including after reverse
complementation. `is_circular=None` uses detected topology; `True` or `False`
overrides it. An explicit start requires a complete effectively circular record,
an integer in `1..L`, and no crop. Source sequence and feature locations remain
unchanged. Circular places that base at 12 o'clock; Linear places it at the left
edge and wraps the complete record. See the [combined executable
example](command-line.md#rotate-a-plastome-and-place-a-multipart-feature) and the
[placement rules](palettes-feature-rules-labels-shapes-and-tracks.md#manual-feature-placement).

## Combined rotation and placement example

Use a new directory containing the four inputs obtained in the
[CLI example](command-line.md#rotate-a-plastome-and-place-a-multipart-feature).
This program constructs the placement table in memory; it does not need
`tables/placements.tsv`. Save it as `rotated_placed_chloroplast.py` and run
`python rotated_placed_chloroplast.py`. It prints
`Saved rotated_placed_chloroplast.svg` and produces the same diagram as the CLI.

<!-- executable:H-PY-06:start -->
```python
from pathlib import Path
from pandas import DataFrame

from gbdraw import (
    CircularOptions,
    RecordDisplayOptions,
    CircularTrackOptions,
    Diagram,
    FeatureOptions,
    LabelOptions,
    draw_circular,
    read_genbank,
)
from gbdraw.api import AnnotationOptions, CircularTrackSlot, ScalarSpec


record = read_genbank(Path("NC_001879.gbk"))[0]
assert (record.id, len(record), record.annotations.get("topology")) == (
    "NC_001879.2",
    155_943,
    "circular",
)

track_slots = (
    CircularTrackSlot(
        id="features",
        renderer="features",
        side="overlay",
        params={"lane_direction": "split"},
    ),
    CircularTrackSlot(
        id="plastome_regions",
        renderer="annotations",
        side="inside",
        radius=ScalarSpec(0.65),
        width=ScalarSpec(20, "px"),
        params={
            "set_id": "plastome_regions",
            "show_labels": True,
            "padding_px": 1,
            "overflow": "compress",
        },
        inner_gap_px=1,
        outer_gap_px=1,
    ),
    CircularTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        side="inside",
        radius=ScalarSpec(0.56),
        width=ScalarSpec(0.08),
        params={"nt": "GC", "legend_label": "GC content"},
    ),
)

options = CircularOptions(
    features=FeatureOptions(
        types=(
            "CDS",
            "rRNA",
            "tRNA",
            "tmRNA",
            "ncRNA",
            "misc_RNA",
            "rep_origin",
        ),
        color_table=Path("chloroplast_specific_table.tsv"),
        placements=DataFrame([{
            "record": "NC_001879.2",
            "feature_selector": "protein_id=NP_054479.1",
            "placement": "outward", "level": 1,
        }]),
    ),
    labels=LabelOptions(
        qualifier_priority=Path("qualifier_priority.tsv"),
    ),
    annotations=AnnotationOptions(
        table_file="nicotiana-tabacum-regions.tsv",
    ),
    tracks=CircularTrackOptions(slots=track_slots),
    species="<i>Nicotiana tabacum</i>",
    legend="upper_left",
    config_overrides={
        "canvas.strandedness": False,
        "canvas.resolve_overlaps": True,
        "canvas.feature_overlap_tolerance_bp": 1,
        "canvas.circular.track_type": "tuckin",
        "labels.circular.scope": "both",
        "labels.circular.placement": "radial",
        "labels.unified_adjustment.outer_labels.x_radius_offset": 0.9,
        "labels.unified_adjustment.outer_labels.y_radius_offset": 0.9,
        "labels.unified_adjustment.inner_labels.x_radius_offset": 0.975,
        "labels.unified_adjustment.inner_labels.y_radius_offset": 0.975,
        "objects.definition.circular.font_size": 28,
        "objects.definition.circular.interval": 30,
        "objects.features.block_stroke_color": "black",
        "objects.features.block_stroke_width.long": 1,
        "objects.features.line_stroke_width.long": 2,
        "objects.axis.circular.stroke_width.long": 3,
    },
)

chloroplast_diagram = draw_circular(
    record, options=options,
    record_displays=[RecordDisplayOptions(start_coordinate=5500)],
)
chloroplast_svg = chloroplast_diagram.to_svg()
chloroplast_bytes = chloroplast_diagram.to_bytes("svg")
chloroplast_path = chloroplast_diagram.save(
    Path("rotated_placed_chloroplast.svg")
)

assert isinstance(chloroplast_diagram, Diagram)
assert chloroplast_diagram.mode == "circular"
assert chloroplast_svg.encode("utf-8") == chloroplast_bytes
assert chloroplast_path.read_bytes() == chloroplast_bytes
print(f"Saved {chloroplast_path}")
```
<!-- executable:H-PY-06:end -->
