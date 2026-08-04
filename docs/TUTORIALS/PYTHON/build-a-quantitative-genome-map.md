[Home](../../DOCS.md) | [Tutorials](../README.md) | [Python API](../../REFERENCE/python-api.md)

# Build a quantitative genome map from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/build-a-quantitative-genome-map.md) | [CLI workflow](../CLI/build-a-quantitative-genome-map.md) | **This page** |

Combine the complete `AP027133.1` record with its measured 1 kbp depth series,
GC content, and GC skew in the same explicit slot stack as the CLI tutorial.

## What you'll need

Place `AP027133.gb` and `AP027133.DRR394922.depth-1kb.tsv` in an empty working
directory.

## Step 1: Save and run the program

<!-- executable:T-PY-11:start -->
```python
from pathlib import Path

from gbdraw import (
    CircularOptions,
    CircularTrackOptions,
    DepthTrackOptions,
    draw_circular,
    read_genbank,
)
from gbdraw.api import CircularTrackSlot, ScalarSpec


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
        id="depth_1",
        renderer="depth",
        side="inside",
        width=ScalarSpec(52, "px"),
        params={"track_index": 0, "legend_label": "DRR394922 mean depth"},
    ),
    CircularTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        side="inside",
        width=ScalarSpec(42, "px"),
        params={"nt": "GC", "legend_label": "GC content (%)"},
    ),
    CircularTrackSlot(
        id="gc_skew",
        renderer="dinucleotide_skew",
        side="inside",
        width=ScalarSpec(34, "px"),
        params={"nt": "GC", "legend_label": "GC skew"},
    ),
)
record = read_genbank(Path("AP027133.gb"))[0]
options = CircularOptions(
    tracks=CircularTrackOptions(slots=track_slots),
    depth_tracks=(
        DepthTrackOptions(
            source=Path("AP027133.DRR394922.depth-1kb.tsv"),
            label="DRR394922 mean depth",
            color="#2563EB",
            large_tick_interval=20,
            small_tick_interval=10,
        ),
    ),
    depth_window=1,
    depth_step=1000,
    window=1000,
    step=1000,
    legend="right",
    config_overrides={
        "canvas.strandedness": True,
        "canvas.show_depth": True,
        "canvas.show_gc": True,
        "canvas.show_skew": True,
        "labels.circular.scope": "none",
        "objects.depth.fill_color": "#2563EB",
        "objects.depth.min_depth": 0,
        "objects.depth.max_depth": 80,
        "objects.depth.normalize": False,
        "objects.depth.show_axis": True,
        "objects.depth.show_ticks": True,
        "objects.depth.large_tick_interval": 20,
        "objects.depth.small_tick_interval": 10,
        "objects.gc_content.mode": "percent",
        "objects.gc_content.min_percent": 10,
        "objects.gc_content.max_percent": 55,
        "objects.gc_content.show_axis": True,
        "objects.gc_content.show_ticks": True,
        "objects.gc_content.large_tick_interval": 10,
        "objects.gc_content.small_tick_interval": 5,
    },
)
diagram = draw_circular(record, options=options)
saved_path = diagram.save(Path("quantitative_genome_map.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-11:end -->

## Step 2: Verify the result

```bash
python docs/recipes/run_python_scenarios.py --scenario T-PY-11 --check
```

The validator checks 607 depth values, the 0x–80x and 10%–55% axes, the signed
GC-skew fills, the five-slot order, and exact parity with the CLI output.

![AP027133.1 depth, GC content, and GC skew from Python](../../images/t-py-11/quantitative_genome_map.svg)

## What you built

One measured table and two sequence-derived series now share explicit radial
ownership without changing the annotated feature ring.
