[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API reference](../../REFERENCE/python-api.md)

# Compare Lambda and DE3 from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-genomes-losatn.md) | [Command line](../CLI/compare-genomes-losatn.md) | **Python API** |

Use the public `LinearComparisonOptions` type to attach the same six-row
LOSATN table used by the GUI and CLI variants.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | `NC_001416.gb` ([NCBI accession NC_001416.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1)) | Download the complete Lambda GenBank record and save it as `NC_001416.gb`. |
| Download | `NC_042057.1.gb` ([NCBI accession NC_042057.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_042057.1)) | Download the complete DE3 GenBank record and save it as `NC_042057.1.gb`. |
| Download | [`lambda-de3.losatn.tsv`](../../../gbdraw/web/tutorial-data/lambda-de3-comparison/lambda-de3.losatn.tsv) | Save the supplied LOSATN result as `lambda-de3.losatn.tsv`. |
| Create | `lambda_de3_losatn.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `python_lambda_de3_losatn.svg` | The program writes this comparison SVG. |
| Reference result | [`python_lambda_de3_losatn.svg`](../../images/t-py-04/python_lambda_de3_losatn.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting. The LOSATN
table is a frozen comparison input for this project.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-losatn
cd gbdraw-python-losatn
```

Acquire both sequences from their authoritative accession links and download
the supplied LOSATN table. Save every file with the exact filename shown. See
[Get the tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser, `curl`,
and PowerShell instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `lambda_de3_losatn.py`:

<!-- executable:T-PY-04:start -->
```python
from pathlib import Path

from gbdraw import (
    LinearComparisonOptions,
    LinearOptions,
    Thresholds,
    draw_linear,
    read_genbank,
)


records = read_genbank([Path("NC_001416.gb"), Path("NC_042057.1.gb")])
assert [(record.id, len(record)) for record in records] == [
    ("NC_001416.1", 48_502),
    ("NC_042057.1", 42_925),
]

options = LinearOptions(
    comparisons=LinearComparisonOptions(
        blast_files=("lambda-de3.losatn.tsv",),
    ),
    thresholds=Thresholds(
        bitscore=50,
        evalue=0.01,
        identity=0,
        alignment_length=0,
    ),
    config_overrides={"canvas.linear.comparison_height": 120},
)
diagram = draw_linear(records, options=options)
saved_path = diagram.save(Path("python_lambda_de3_losatn.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-04:end -->

Before running it, your working directory should contain:

```text
gbdraw-python-losatn/
├── NC_001416.gb
├── NC_042057.1.gb
├── lambda-de3.losatn.tsv
└── lambda_de3_losatn.py
```

### Run the program

```bash
python lambda_de3_losatn.py
```

Expected output: the program prints `Saved python_lambda_de3_losatn.svg` and
writes the Generated SVG in the current directory.

## Step 3: Inspect the comparison

Open the Generated `python_lambda_de3_losatn.svg` to inspect the two complete
records and six comparison ribbons.

![Lambda and DE3 LOSATN comparison from Python](../../images/t-py-04/python_lambda_de3_losatn.svg)

The image above is the Reference result. Verify the two accessions, source
lengths, record order, six retained matches, and 120 px comparison band.

The validator checks the record order, source lengths, six retained matches,
comparison height, and standard-SVG safety.

## Next steps

- [Draw Linear diagrams and comparisons from Python](../../HOW_TO/PYTHON/draw-linear-diagrams-and-comparisons.md)
- [Review comparison thresholds and result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
