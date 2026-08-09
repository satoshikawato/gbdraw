[Home](../../DOCS.md) | [Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API](../../REFERENCE/python-api.md) | [Export](../../REFERENCE/output-formats-and-export.md)

# Draw and save your first genome diagram from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/first-circular-genome-diagram.md) | [Command-line workflow](../CLI/first-circular-genome-diagram.md) | **This page** |

All three workflows build the same labeled human mitochondrial figure. Only
the interface changes.

You will load the human mitochondrial reference genome, draw one Circular
diagram, and save it as a standard SVG with the beginner-facing Python API.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | `HmmtDNA.gbk` ([NCBI accession NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)) | Download the complete GenBank record and save it as `HmmtDNA.gbk`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save the supplied label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `first_diagram.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `python_human_mitochondrion.svg` | The program writes this SVG. |
| Reference result | [`python_human_mitochondrion.svg`](../../images/t-py-01/python_human_mitochondrion.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-circular
cd gbdraw-python-circular
```

Acquire the sequence from its authoritative accession link and download the
supplied label rule. Save both with the exact filenames shown. [Get the
tutorial files](../../GETTING_TUTORIAL_DATA.md) gives browser, `curl`, and
PowerShell instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `first_diagram.py` beside the two
Download files:

<!-- executable:T-PY-01:start -->
```python
from pathlib import Path

from gbdraw import (
    CircularOptions,
    Diagram,
    LabelOptions,
    draw_circular,
    read_genbank,
)


input_path = Path("HmmtDNA.gbk")
output_path = Path("python_human_mitochondrion.svg")

record = read_genbank(input_path)[0]
options = CircularOptions(
    labels=LabelOptions(
        qualifier_priority=Path("cds_gene_qualifier_priority.tsv"),
    ),
    species="<i>Homo sapiens</i>",
    legend="right",
    config_overrides={
        "canvas.strandedness": True,
        "canvas.circular.track_type": "middle",
        "labels.circular.scope": "outer",
        "labels.circular.placement": "horizontal",
    },
)
diagram = draw_circular(record, options=options)
saved_path = diagram.save(output_path)

assert isinstance(diagram, Diagram)
assert saved_path == output_path
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-01:end -->

`read_genbank()` returns a list because one GenBank file can contain more than
one record. This tutorial selects the only record in the downloaded file. The
options match the GUI and CLI variants: separate strands, Middle layout, outer labels,
gene-first CDS label priority, and a right legend.

Before running it, your working directory should contain:

```text
gbdraw-python-circular/
├── HmmtDNA.gbk
├── cds_gene_qualifier_priority.tsv
└── first_diagram.py
```

### Run the program

```bash
python first_diagram.py
```

Expected output: the program prints `Saved python_human_mitochondrion.svg`.
Within the program, `draw_circular()` returns a `Diagram`, and `save()` writes
the SVG in the current directory.

## Step 3: Inspect the SVG

Open the Generated `python_human_mitochondrion.svg`. It contains
`NC_012920.1`, 37 displayed features, coordinate ticks, GC content, and GC
skew.

![Circular human mitochondrial genome drawn from Python](../../images/t-py-01/python_human_mitochondrion.svg)

The image above is the Reference result. Your SVG should have the same record,
track order, labels, and legend even if XML metadata differs.

The recipe runner extracts the exact Python block above, executes it in a clean
temporary directory, verifies that `diagram` is a `Diagram`, and checks the
same SVG semantics as the CLI tutorial.

## If the program fails

- `ModuleNotFoundError: No module named 'gbdraw'`: activate the environment
  where gbdraw is installed.
- `Output file already exists`: use a new empty directory or change
  `output_path`. `Diagram.save()` does not overwrite by default.
