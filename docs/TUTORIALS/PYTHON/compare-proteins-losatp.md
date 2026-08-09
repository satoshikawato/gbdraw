[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Python API reference](../../REFERENCE/python-api.md)

# Create protein Similarity groups with LOSATP from Python

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp.md) | [Command line](../CLI/compare-proteins-losatp.md) | **Python API** |

This program loads five complete BGC records, reverses the fifth display,
runs gbdraw's one-thread LOSATP workflow, and aligns every row to `og_1`.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`BGC0000708.gbk` from MIBiG BGC0000708.5](https://mibig.secondarymetabolites.org/repository/BGC0000708.5/BGC0000708.gbk) | Download the original MIBiG GenBank record and keep the filename `BGC0000708.gbk`. |
| Download | [`BGC0000709.gbk` from MIBiG BGC0000709.5](https://mibig.secondarymetabolites.org/repository/BGC0000709.5/BGC0000709.gbk) | Download the original MIBiG GenBank record and keep the filename `BGC0000709.gbk`. |
| Download | [`BGC0000711.gbk` from MIBiG BGC0000711.5](https://mibig.secondarymetabolites.org/repository/BGC0000711.5/BGC0000711.gbk) | Download the original MIBiG GenBank record and keep the filename `BGC0000711.gbk`. |
| Download | [`BGC0000712.gbk` from MIBiG BGC0000712.5](https://mibig.secondarymetabolites.org/repository/BGC0000712.5/BGC0000712.gbk) | Download the original MIBiG GenBank record and keep the filename `BGC0000712.gbk`. |
| Download | [`BGC0000713.gbk` from MIBiG BGC0000713.5](https://mibig.secondarymetabolites.org/repository/BGC0000713.5/BGC0000713.gbk) | Download the original MIBiG GenBank record and keep the filename `BGC0000713.gbk`. |
| Download | [`BGC0000708-BGC0000713_default_colors.tsv`](../../../gbdraw/web/tutorial-data/aminoglycoside-bgc-five/BGC0000708-BGC0000713_default_colors.tsv) | Save the supplied base color table with this exact filename. |
| Download | [`BGC0000708-BGC0000713_specific_colors.tsv`](../../../gbdraw/web/tutorial-data/aminoglycoside-bgc-five/BGC0000708-BGC0000713_specific_colors.tsv) | Save the supplied feature color table with this exact filename. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save the supplied label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `bgc_losatp_groups.py` | Save the complete Python program from Step 2 with this filename. |
| Generated | `python_bgc_losatp_groups.svg` | The program writes this SVG after running LOSATP. |
| Reference result | [`python_bgc_losatp_groups.svg`](../../images/t-py-05/python_bgc_losatp_groups.svg) | Compare your Generated SVG with this versioned result. |

Install gbdraw in the active Python environment before starting.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-python-losatp-groups
cd gbdraw-python-losatp-groups
```

Download the five original sequences from their direct MIBiG source links and
download the three supplied support tables. Save every file with the exact
filename shown. See [Get the tutorial files](../../GETTING_TUTORIAL_DATA.md)
for browser, `curl`, and PowerShell instructions.

## Step 2: Save and run the Python program

Save the following complete program as `bgc_losatp_groups.py`:

<!-- executable:T-PY-05:start -->
```python
from pathlib import Path

from gbdraw.api import (
    ColorOptions,
    InMemoryRecordSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    LinearOutputOptions,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
    load_gbks,
    render_request,
)


paths = [
    Path("BGC0000708.gbk"),
    Path("BGC0000709.gbk"),
    Path("BGC0000711.gbk"),
    Path("BGC0000712.gbk"),
    Path("BGC0000713.gbk"),
]
records = load_gbks(paths)
presentations = [
    ("<i>Streptomyces lividus</i> CBS 844.73", "Lividomycin biosynthetic gene cluster"),
    ("<i>Streptomyces fradiae</i> ATCC 10745", "Neomycin biosynthetic gene cluster"),
    ("<i>Streptomyces fradiae</i> MCIMB 8233", "Neomycin biosynthetic gene cluster"),
    ("<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455", "Paromomycin biosynthetic gene cluster"),
    ("<i>Streptomyces ribosidificus</i> ATCC 21294", "Ribostamycin biosynthetic gene"),
]
request = LinearDiagramRequest(
    records=tuple(
        RecordInput(
            source=InMemoryRecordSource(record),
            presentation=RecordPresentation(
                label=label,
                subtitle=subtitle,
                reverse_complement=index == 4,
            ),
        )
        for index, (record, (label, subtitle)) in enumerate(
            zip(records, presentations, strict=True)
        )
    ),
    options=LinearDiagramOptions(
        colors=ColorOptions(
            color_table_file="BGC0000708-BGC0000713_specific_colors.tsv",
            default_colors_file="BGC0000708-BGC0000713_default_colors.tsv",
            default_colors_palette="orange",
        ),
        qualifier_priority_file="cds_gene_qualifier_priority.tsv",
        plot_title="LOSATP Similarity groups across five whole BGC records",
        output=LinearOutputOptions(
            legend="bottom",
            plot_title_position="bottom",
        ),
        protein_blastp_mode="orthogroup",
        losatp_threads=1,
        align_orthogroup_feature="CAG38695.1",
        pairwise_match_style="curve",
        bitscore=50,
        evalue=0.01,
        identity=30,
        alignment_length=0,
        config_overrides={
            "labels.linear.scope": "first",
            "labels.font_size.linear.short": 18,
            "labels.font_size.linear.long": 18,
            "labels.linear.placement": "above_feature",
            "labels.linear.rotation": 45,
            "canvas.linear.default_cds_height.short": 75,
            "canvas.linear.default_cds_height.long": 75,
            "objects.features.block_stroke_color": "#262626",
            "objects.features.block_stroke_width.short": 2,
            "objects.features.block_stroke_width.long": 2,
            "objects.features.line_stroke_width.short": 2,
            "objects.features.line_stroke_width.long": 2,
            "objects.axis.linear.stroke_width.short": 5,
            "objects.axis.linear.stroke_width.long": 5,
            "objects.scale.style": "ruler",
            "canvas.linear.track_layout": "middle",
            "canvas.linear.keep_definition_left_aligned": True,
            "objects.definition.linear.line_styles.name.font_size": 20,
            "objects.definition.linear.line_styles.name.font_weight": "bold",
            "objects.definition.linear.line_styles.subtitle.font_size": 20,
            "objects.definition.linear.line_styles.accession.font_size": 20,
            "objects.definition.linear.line_styles.accession.fill": "#7b7c7d",
            "objects.definition.linear.line_styles.length.font_size": 20,
            "objects.definition.linear.line_styles.length.fill": "#7b7c7d",
            "canvas.strandedness": False,
        },
    ),
    output=RenderOutputRequest(
        output_prefix="python_bgc_losatp_groups",
        formats=("svg",),
    ),
)
diagram = render_request(request)
saved_path = diagram.output_paths[0]
print(f"Saved {saved_path}")
```
<!-- executable:T-PY-05:end -->

Before running it, your working directory should contain:

```text
gbdraw-python-losatp-groups/
├── BGC0000708.gbk
├── BGC0000709.gbk
├── BGC0000711.gbk
├── BGC0000712.gbk
├── BGC0000713.gbk
├── BGC0000708-BGC0000713_default_colors.tsv
├── BGC0000708-BGC0000713_specific_colors.tsv
├── cds_gene_qualifier_priority.tsv
└── bgc_losatp_groups.py
```

### Run LOSATP and draw the figure

```bash
python bgc_losatp_groups.py
```

Expected output: gbdraw's LOSAT runtime performs four adjacent searches.
The program then prints `Saved python_bgc_losatp_groups.svg` and writes the
Generated SVG in the current directory.

## Step 3: Inspect the comparison

Open the Generated `python_bgc_losatp_groups.svg` and confirm the five records
appear in the documented order with the fifth record reversed.

![Five BGC records aligned to og_1 from Python](../../images/t-py-05/python_bgc_losatp_groups.svg)

The image above is the Reference result. Compare its record labels, alignment,
definition column, feature colors, and curved links with your SVG.

The fixed record order and reversed fifth record match the browser Tutorial.
The run should produce 23 Similarity groups and 77 adjacent links, with every
record aligned on `og_1`.

## Next steps

- [Browse focused Python task guides](../../HOW_TO/PYTHON/README.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
