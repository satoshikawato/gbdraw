[Home](../../DOCS.md) | [Tutorials](../README.md) | [Python API](../../REFERENCE/python-api.md)

# Highlight mitochondrial features from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/highlight-mitochondrial-features.md) | [CLI workflow](../CLI/highlight-mitochondrial-features.md) | **This page** |

Apply the CLI project's color, label, shape, stroke, and D-loop rules through
the public Python API without changing the GenBank record.

## What you'll need

Place `HmmtDNA.gbk` and `cds_gene_qualifier_priority.tsv` in an empty working
directory.

## Step 1: Save and run the program

<!-- executable:T-PY-09:start -->
```python
from pathlib import Path

from pandas import DataFrame

from gbdraw import (
    CircularOptions,
    CircularTrackOptions,
    FeatureOptions,
    LabelOptions,
    TitleOptions,
    draw_circular,
    read_genbank,
)
from gbdraw.api import AnnotationOptions, CircularTrackSlot, ScalarSpec


color_table = DataFrame(
    [
        ["CDS", "gene", "^ND[1-6]$", "#3B82F6", "NADH dehydrogenase"],
        ["CDS", "gene", "^COX[1-3]$", "#EF4444", "Cytochrome c oxidase"],
        ["CDS", "gene", "^ATP[68]$", "#F59E0B", "ATP synthase"],
        ["CDS", "gene", "^CYTB$", "#8B5CF6", "Cytochrome b"],
        ["rRNA", "gene", "^RNR[12]$", "#10B981", "Ribosomal RNA"],
    ],
    columns=["feature_type", "qualifier_key", "value", "color", "caption"],
)
label_whitelist = DataFrame(
    [
        ["CDS", "gene", "^(ND1|COX2|ATP8|ATP6|COX3|CYTB)$"],
        ["rRNA", "gene", "^RNR[12]$"],
    ],
    columns=["feature_type", "qualifier", "keyword"],
)
label_overrides = DataFrame(
    [
        ["NC_012920.1", "CDS", "label", "^ND1$", "Complex I (ND1)"],
        ["NC_012920.1", "CDS", "label", "^COX2$", "Oxidase II"],
        ["NC_012920.1", "rRNA", "label", "^s-rRNA$", "12S rRNA"],
        ["NC_012920.1", "rRNA", "label", "^l-rRNA$", "16S rRNA"],
    ],
    columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
)
regions = DataFrame(
    [[
        "mitochondrial_regions", "d_loop", "bracket", "NC_012920.1",
        16024, 576, "source", True, "D-loop", 0, "#202020", 3, "tick",
        "#202020", 14, "tangent", 7,
    ]],
    columns=[
        "set_id", "id", "mark", "record", "start", "end",
        "coordinate_space", "wraps_origin", "label", "lane", "stroke",
        "stroke_width", "line_cap", "label_color", "label_font_size",
        "label_orientation", "label_offset",
    ],
)
track_slots = (
    CircularTrackSlot(
        id="ticks",
        renderer="ticks",
        side="outside",
        params={"tick_label_layout": "label_out_tick_in"},
    ),
    CircularTrackSlot(
        id="features",
        renderer="features",
        side="overlay",
        params={"lane_direction": "split"},
    ),
    CircularTrackSlot(
        id="mitochondrial_regions",
        renderer="annotations",
        side="inside",
        width=ScalarSpec(24, "px"),
        params={
            "set_id": "mitochondrial_regions",
            "show_labels": True,
            "padding_px": 1,
            "overflow": "compress",
        },
    ),
)

record = read_genbank(Path("HmmtDNA.gbk"))[0]
options = CircularOptions(
    features=FeatureOptions(
        types=("CDS", "rRNA", "tRNA"),
        color_table=color_table,
        shapes={"CDS": "arrow", "rRNA": "rectangle", "tRNA": "arrow"},
    ),
    labels=LabelOptions(
        qualifier_priority="cds_gene_qualifier_priority.tsv",
        whitelist=label_whitelist,
        overrides=label_overrides,
    ),
    annotations=AnnotationOptions(table=regions),
    tracks=CircularTrackOptions(slots=track_slots),
    species="<i>Homo sapiens</i>",
    title=TitleOptions(
        text="Human mitochondrial feature presentation",
        position="top",
    ),
    legend="right",
    config_overrides={
        "canvas.resolve_overlaps": False,
        "canvas.strandedness": True,
        "canvas.circular.track_type": "spreadout",
        "labels.circular.scope": "both",
        "labels.rendering": "auto",
        "objects.features.arrow_geometry.head_length_ratio": "auto",
        "objects.features.arrow_geometry.shaft_width_ratio": 0.72,
        "objects.features.block_stroke_color": "#1F2937",
        "objects.features.block_stroke_width.short": 1.5,
        "objects.features.block_stroke_width.long": 1.5,
        "objects.axis.circular.stroke_color": "#374151",
        "objects.axis.circular.stroke_width.short": 4,
        "objects.axis.circular.stroke_width.long": 4,
        "objects.features.line_stroke_color": "#9CA3AF",
        "objects.features.line_stroke_width.short": 1.5,
        "objects.features.line_stroke_width.long": 1.5,
    },
)
diagram = draw_circular(record, options=options)
saved_path = diagram.save(Path("mitochondrial_features_highlighted.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-09:end -->

## Step 2: Verify the result

```bash
python docs/recipes/run_python_scenarios.py --scenario T-PY-09 --check
```

The validator checks all 37 features, the five functional colors, selected
renamed labels, arrow and rectangle shapes, stroke geometry, and the one
origin-spanning D-loop bracket.

![Highlighted human mitochondrial features from Python](../../images/t-py-09/mitochondrial_features_highlighted.svg)

## What you built

The code changes only presentation objects. The complete `NC_012920.1` record
and every CDS remain intact.
