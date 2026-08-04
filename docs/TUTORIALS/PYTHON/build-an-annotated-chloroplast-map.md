[Home](../../DOCS.md) | [Tutorials](../README.md) | [Python API](../../REFERENCE/python-api.md) | [GUI version](../GUI/build-an-annotated-chloroplast-map.md)

# Recreate the Gallery chloroplast map from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/build-an-annotated-chloroplast-map.md) | [Command-line workflow](../CLI/build-an-annotated-chloroplast-map.md) | **This page** |

All three workflows use the same complete tobacco plastome, tables, visual
settings, and track geometry. Only the interface changes.

Use the public Python API to reproduce the same complete tobacco plastome map
as the GUI tutorial: chloroplast gene-family colors, radial labels on both
sides of the feature ring, LSC/IRb/SSC/IRa brackets, GC content, and an
upper-left legend.

This is the next step after the
[first Python genome diagram](first-genome-diagram.md). It keeps the raw
GenBank record unchanged and expresses each presentation decision in typed
options that can be reviewed and rerun.

## What you'll need

Copy these files from
`gbdraw/web/tutorial-data/tobacco-plastome-regions/` into an empty working
directory:

- `NC_001879.gbk`
- `nicotiana-tabacum-regions.tsv`
- `chloroplast_specific_table.tsv`
- `qualifier_priority.tsv`

The GenBank file contains the complete 155,943 bp `NC_001879.2` plastome. The
three TSV files own the structural-region annotations, functional feature
colors, and CDS label priority.

## Step 1: Save the Python program

Save the following program as `annotated_chloroplast.py` beside the four input
files:

<!-- executable:T-PY-02:start -->
```python
from pathlib import Path

from gbdraw import (
    CircularOptions,
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
        "canvas.strandedness": True,
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

chloroplast_diagram = draw_circular(record, options=options)
chloroplast_svg = chloroplast_diagram.to_svg()
chloroplast_bytes = chloroplast_diagram.to_bytes("svg")
chloroplast_path = chloroplast_diagram.save(
    Path("python_annotated_chloroplast.svg")
)

assert isinstance(chloroplast_diagram, Diagram)
assert chloroplast_diagram.mode == "circular"
assert chloroplast_svg.encode("utf-8") == chloroplast_bytes
assert chloroplast_path.read_bytes() == chloroplast_bytes
print(f"Saved {chloroplast_path}")
```
<!-- executable:T-PY-02:end -->

The first slot splits forward- and reverse-strand features around the axis.
The second places all four structural regions in one inner annotation lane.
The third adds GC content without adding coordinate ticks or a skew ring.

## Step 2: Run the program

```bash
python annotated_chloroplast.py
```

The program prints `Saved python_annotated_chloroplast.svg`.

## Step 3: Inspect the result

Open `python_annotated_chloroplast.svg`. Confirm the complete
`NC_001879.2` record, 147 logical features, radial gene labels, all four
structural-region brackets, GC content, and the functional-color legend.
There should be no coordinate-tick or GC-skew track.

![Gallery-style tobacco plastome drawn from the Python API](../../images/t-py-02/python_annotated_chloroplast.svg)

The documentation runner executes the literal program above in a clean
directory and checks the typed options, slot geometry, feature and annotation
identities, visible labels, and equality of the file, text, and byte outputs.
