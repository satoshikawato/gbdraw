[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Python API reference](../../REFERENCE/python-api.md)

# Reproduce the Hepatoplasmataceae Collinear map from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp-collinear.md) | [Command line](../CLI/compare-proteins-losatp-collinear.md) | **Python API** |

Load the five complete source records and run the same adjacent Collinear
analysis through the public Python API.

## Step 1: Prepare the records

Place `AP027078.gb`, `AP027131.gb`, `AP027133.gb`, `AP027132.gb`, and
`NZ_CP006932.gb` in an empty working directory.

## Step 2: Replay a standard SVG

<!-- executable:T-PY-07:start -->
```python
from pathlib import Path

from gbdraw import (
    FeatureOptions,
    LinearComparisonOptions,
    LinearOptions,
    Thresholds,
    TitleOptions,
    draw_linear,
    read_genbank,
)
from gbdraw.api import LosslessCollinearityParameters


records = read_genbank(
    [
        Path("AP027078.gb"),
        Path("AP027131.gb"),
        Path("AP027133.gb"),
        Path("AP027132.gb"),
        Path("NZ_CP006932.gb"),
    ]
)
options = LinearOptions(
    features=FeatureOptions(palette="ajisai"),
    comparisons=LinearComparisonOptions(
        protein_mode="collinear",
        threads=32,
        match_style="curve",
        collinearity_scope="adjacent",
        collinearity_color="orientation_identity",
        collinearity_params=LosslessCollinearityParameters(
            min_anchors=1,
            max_unit_gap=0,
            max_diagonal_drift=0,
            max_conflicts=1,
        ),
    ),
    thresholds=Thresholds(
        bitscore=50,
        evalue=0.01,
        identity=0,
        alignment_length=0,
    ),
    title=TitleOptions(
        text="LOSATP Collinear blocks across Hepatoplasmataceae",
        position="top",
    ),
    legend="right",
    config_overrides={
        "canvas.linear.track_layout": "middle",
        "canvas.linear.align_center": True,
        "canvas.strandedness": True,
        "canvas.show_gc": True,
        "canvas.show_skew": True,
        "objects.scale.style": "ruler",
    },
)
diagram = draw_linear(records, options=options)
saved_path = diagram.save(Path("losatp_collinear.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-07:end -->

![Five Hepatoplasmataceae genomes with Collinear blocks from Python](../../images/t-py-07/losatp_collinear.svg)

## Step 3: Verify the replay

```bash
python docs/recipes/run_python_scenarios.py --scenario T-PY-07 --check
```

The runner checks the five-record order, adjacent comparison pairs, 500
displayed matches, track stack, title, and standard-SVG safety.

## What you built

You reproduced the Gallery Collinear result from source records with explicit
search, grouping, color, and track settings.
