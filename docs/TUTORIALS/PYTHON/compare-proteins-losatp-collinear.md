[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API reference](../../REFERENCE/python-api.md)

# Reproduce the Hepatoplasmataceae Collinear map from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp-collinear.md) | [Command line](../CLI/compare-proteins-losatp-collinear.md) | **Python API** |

Load the five complete source records and run the same adjacent Collinear
analysis through the public Python API.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`AP027078.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027078.1&sat=3&satkey=69902295) | Download the official NCBI Revision History snapshot of `AP027078.1` and save it as `AP027078.gb`. |
| Download | [`AP027131.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027131.1&sat=3&satkey=69902296) | Download the official NCBI Revision History snapshot of `AP027131.1` and save it as `AP027131.gb`. |
| Download | [`AP027133.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027133.1&sat=3&satkey=69902298) | Download the official NCBI Revision History snapshot of `AP027133.1` and save it as `AP027133.gb`. |
| Download | [`AP027132.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027132.1&sat=3&satkey=69902297) | Download the official NCBI Revision History snapshot of `AP027132.1` and save it as `AP027132.gb`. |
| Download | [`NZ_CP006932.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=NZ_CP006932.1&sat=60&satkey=39275474) | Download the official NCBI Revision History snapshot of `NZ_CP006932.1` in full GenBank format and save it as `NZ_CP006932.gb`. |
| Create | `losatp_collinear.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `losatp_collinear.svg` | The program writes this SVG after running LOSATP. |
| Reference result | [`losatp_collinear.svg`](../../images/t-py-07/losatp_collinear.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting.

Use the five Revision History links, not the live accession downloads. NCBI can
update a record's annotation without changing its sequence accession version.
These pinned revisions preserve the exact feature tables used to reproduce this
Tutorial's 2,994 displayed features and 500 Collinear matches.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-losatp-collinear
cd gbdraw-python-losatp-collinear
```

Acquire all five sequences from their authoritative NCBI links and save them
with the exact filenames shown. See [Get the tutorial
files](../../GETTING_TUTORIAL_DATA.md) for browser, `curl`, and PowerShell
instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `losatp_collinear.py`:

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

Before running it, your working directory should contain:

```text
gbdraw-python-losatp-collinear/
├── AP027078.gb
├── AP027131.gb
├── AP027133.gb
├── AP027132.gb
├── NZ_CP006932.gb
└── losatp_collinear.py
```

### Run LOSATP and draw the figure

```bash
python losatp_collinear.py
```

Expected output: gbdraw's LOSAT runtime performs adjacent protein searches.
The program then prints `Saved losatp_collinear.svg` and writes the Generated
SVG in the current directory.

## Step 3: Inspect the Collinear map

Open the Generated `losatp_collinear.svg` and check the five complete records
in the documented order, 500 rendered Collinear matches, centered alignment,
rulers, GC content, and GC skew.

![Five Hepatoplasmataceae genomes with Collinear blocks from Python](../../images/t-py-07/losatp_collinear.svg)

The image above is the Reference result. Compare the record layout and the
adjacent Collinear blocks colored by orientation and identity with your SVG.

## Next steps

- [Python API technical documentation](../../REFERENCE/python-api.md)
- [Choose a genome-comparison method](../../FAQ.md#which-comparison-method-should-i-use)
