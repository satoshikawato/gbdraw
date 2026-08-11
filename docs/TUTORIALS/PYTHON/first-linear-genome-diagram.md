[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API reference](../../REFERENCE/python-api.md)

# Create a reproducible linear diagram from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/first-linear-genome-diagram.md) | [Command line](../CLI/first-linear-genome-diagram.md) | **Python API** |

This project draws the complete 48,502 bp Lambda record with separated feature
strands, concise labels, a centered ruler, and the legend on the left.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | `NC_001416.gb` ([NCBI accession NC_001416.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1)) | Download the complete Lambda GenBank record and save it as `NC_001416.gb`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save the supplied label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `lambda_linear.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `python_lambda_linear.svg` | The program writes this SVG. |
| Reference result | [`python_lambda_linear.svg`](../../images/t-py-03/python_lambda_linear.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-linear
cd gbdraw-python-linear
```

Acquire the sequence from its authoritative accession link and download the
supplied label rule. Save both with the exact filenames shown. See [Get the
tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser, `curl`, and
PowerShell instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `lambda_linear.py`:

<!-- executable:T-PY-03:start -->
```python
from pathlib import Path

from gbdraw import LabelOptions, LinearOptions, draw_linear, read_genbank


record = read_genbank([Path("NC_001416.gb")])[0]
assert record.id == "NC_001416.1" and len(record) == 48_502

options = LinearOptions(
    labels=LabelOptions(
        qualifier_priority="cds_gene_qualifier_priority.tsv",
    ),
    legend="left",
    config_overrides={
        "labels.linear.scope": "all",
        "canvas.strandedness": True,
        "canvas.linear.track_layout": "middle",
        "objects.scale.style": "ruler",
    },
)
diagram = draw_linear(record, options=options)
saved_path = diagram.save(Path("python_lambda_linear.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-03:end -->

Before running it, your working directory should contain:

```text
gbdraw-python-linear/
├── NC_001416.gb
├── cds_gene_qualifier_priority.tsv
└── lambda_linear.py
```

### Run the program

```bash
python lambda_linear.py
```

Expected output: the program prints `Saved python_lambda_linear.svg` and
writes the Generated SVG in the current directory.

## Step 3: Inspect the SVG

Open the Generated `python_lambda_linear.svg` to see the labeled whole-record
map, all 73 CDS features, the ruler, and the two separated strands.

![Complete Lambda linear map from Python](../../images/t-py-03/python_lambda_linear.svg)

The image above is the Reference result. Your SVG should show the same complete
record, labels, ruler, separated strands, and left legend.
