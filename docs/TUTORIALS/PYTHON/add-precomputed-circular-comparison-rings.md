[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API reference](../../REFERENCE/python-api.md)

# Add precomputed circular comparison rings from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/add-precomputed-circular-comparison-rings.md) | [Command line](../CLI/add-precomputed-circular-comparison-rings.md) | **Python API** |

Build the same ordered three-ring map with the public Circular comparison
types.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | `HmmtDNA.gbk` ([NCBI accession NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)) | Download the complete human mitochondrial GenBank record and save it as `HmmtDNA.gbk`. |
| Download | [`danio-human.tlosatx.tsv`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-comparison/danio-human.tlosatx.tsv) | Save the supplied Danio comparison table with this exact filename. |
| Download | [`drosophila-human.tlosatx.tsv`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-comparison/drosophila-human.tlosatx.tsv) | Save the supplied Drosophila comparison table with this exact filename. |
| Download | [`caenorhabditis-human.tlosatx.tsv`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-comparison/caenorhabditis-human.tlosatx.tsv) | Save the supplied Caenorhabditis comparison table with this exact filename. |
| Download | `NC_002333.2.fna` ([NCBI accession NC_002333.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002333.2)) | Download the complete Danio FASTA record and save it as `NC_002333.2.fna`. |
| Download | `NC_024511.2.fna` ([NCBI accession NC_024511.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_024511.2)) | Download the complete Drosophila FASTA record and save it as `NC_024511.2.fna`. |
| Download | `NC_001328.1.fna` ([NCBI accession NC_001328.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001328.1)) | Download the complete Caenorhabditis FASTA record and save it as `NC_001328.1.fna`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save the supplied label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `precomputed_circular_rings.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `python_precomputed_circular_rings.svg` | The program writes this SVG. |
| Reference result | [`python_precomputed_circular_rings.svg`](../../images/t-py-06/python_precomputed_circular_rings.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting. The three
TLOSATX tables are frozen comparison inputs for this project.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-precomputed-rings
cd gbdraw-python-precomputed-rings
```

Acquire the GenBank and FASTA sequences from their authoritative accession
links and download the four supplied support tables. Save every file with the
exact filename shown. See [Get the tutorial
files](../../GETTING_TUTORIAL_DATA.md) for browser, `curl`, and PowerShell
instructions for authoritative sequence acquisition.

## Step 2: Save and run the Python program

Save the following complete program as `precomputed_circular_rings.py`:

<!-- executable:T-PY-06:start -->
```python
from pathlib import Path

from gbdraw import (
    CircularOptions,
    ComparisonRingOptions,
    ComparisonRingTrackOptions,
    LabelOptions,
    Thresholds,
    TitleOptions,
    draw_circular,
    read_genbank,
)


record = read_genbank([Path("HmmtDNA.gbk")])[0]
rings = ComparisonRingOptions(
    tracks=(
        ComparisonRingTrackOptions(
            source="danio-human.tlosatx.tsv",
            label="Danio rerio (NC_002333.2)",
            color="#4E79A7",
            comparison_sequence_source="NC_002333.2.fna",
        ),
        ComparisonRingTrackOptions(
            source="drosophila-human.tlosatx.tsv",
            label="Drosophila melanogaster (NC_024511.2)",
            color="#F28E2B",
            comparison_sequence_source="NC_024511.2.fna",
        ),
        ComparisonRingTrackOptions(
            source="caenorhabditis-human.tlosatx.tsv",
            label="Caenorhabditis elegans (NC_001328.1)",
            color="#59A14F",
            comparison_sequence_source="NC_001328.1.fna",
        ),
    ),
    reference="subject",
    ring_width=18,
    ring_gap=4,
)
options = CircularOptions(
    comparison_rings=rings,
    labels=LabelOptions(
        qualifier_priority="cds_gene_qualifier_priority.tsv",
    ),
    thresholds=Thresholds(
        bitscore=50,
        evalue=1e-5,
        identity=40,
        alignment_length=50,
    ),
    species="<i>Homo sapiens</i>",
    title=TitleOptions(
        text="Precomputed TLOSATX rings around Homo sapiens mtDNA",
        position="bottom",
    ),
    legend="right",
    config_overrides={
        "canvas.strandedness": False,
        "canvas.circular.track_type": "middle",
        "labels.circular.scope": "outer",
        "objects.definition.circular.font_size": 18,
    },
)
diagram = draw_circular(record, options=options)
saved_path = diagram.save(Path("python_precomputed_circular_rings.svg"))
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-06:end -->

Before running it, your working directory should contain:

```text
gbdraw-python-precomputed-rings/
├── HmmtDNA.gbk
├── danio-human.tlosatx.tsv
├── drosophila-human.tlosatx.tsv
├── caenorhabditis-human.tlosatx.tsv
├── NC_002333.2.fna
├── NC_024511.2.fna
├── NC_001328.1.fna
├── cds_gene_qualifier_priority.tsv
└── precomputed_circular_rings.py
```

### Draw and save the rings

```bash
python precomputed_circular_rings.py
```

Expected output: the program prints
`Saved python_precomputed_circular_rings.svg` and writes the Generated SVG in
the current directory.

## Step 3: Inspect the rings

Open the Generated `python_precomputed_circular_rings.svg` and check the three
ring labels and their Danio, Drosophila, Caenorhabditis order.

![Human mitochondrial reference with three rings from Python](../../images/t-py-06/python_precomputed_circular_rings.svg)

The image above is the Reference result. Compare the subject-reference
direction, ring widths, gaps, labels, colors, and order with your SVG.

The saved SVG should keep the ordered ring labels, the subject-reference
mapping, and 106 retained HSPs across the three companion sequences.

## Next steps

- [Browse focused Python task guides](../../HOW_TO/PYTHON/README.md)
- [Review comparison result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
