[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [Web app](../../HOW_TO/GUI/README.md)

# Create protein Similarity groups with LOSATP in the browser

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/compare-proteins-losatp.md) | [Python workflow](../PYTHON/compare-proteins-losatp.md) |

Use LOSATP **Similarity groups** to compare five aminoglycoside biosynthetic
gene-cluster records. The records remain Linear. They are native BGC database
regions, not complete chromosomes, so this Tutorial does not turn them into
circular genomes.

![Five-record BGC comparison aligned to similarity group og_1](../../images/t-gui-04/05-comparison-result.png)

## What you'll need

Use the five files under
`gbdraw/web/tutorial-data/aminoglycoside-bgc-five/` in this order:

| Order | Record | Length | CDS |
|---:|---|---:|---:|
| 1 | `BGC0000708` | 40,579 bp | 30 |
| 2 | `BGC0000709` | 50,466 bp | 38 |
| 3 | `BGC0000711` | 30,837 bp | 21 |
| 4 | `BGC0000712` | 48,169 bp | 40 |
| 5 | `BGC0000713` | 31,892 bp | 26 |

Upload every file as a whole record and leave all optional Region fields blank.
Keep the first four records in their source orientation; turn on
**Reverse complement** only for `BGC0000713` to reproduce the Gallery layout.

## Step 1: Load the five Linear records

Select **Linear**, **GenBank**, and **No comparison**. Upload `BGC0000708.gbk`,
then use **Add Seq** four times and upload the remaining files in the table
order.

For the fifth row, `BGC0000713`, turn on **Reverse complement**. This changes
only its display orientation; it does not crop, split, or alter the source
record.

![Five annotated BGC records selected for protein comparison](../../images/t-gui-04/01-input-ready.png)

## Step 2: Draw the records before comparing them

Select **Generate Diagram**. The first result contains 155 CDS features and no
comparison ribbons.

![Plain five-record BGC linear diagram](../../images/t-gui-04/02-first-diagram.png)

## Step 3: Configure Similarity groups

Select **Run LOSAT**, **LOSATP**, and **Similarity groups**. Use deterministic
execution and filtering:

| Control | Value |
|---|---|
| Execution | Serial |
| Total threads | `1` |
| Parallel runs | `1 run` |
| Threads per run | `1` |
| Bitscore | `50` |
| E-value | `0.01` |
| Minimum identity | `30` |
| Minimum length | `0` |
| Pairwise Match Style | Curve |
| Separate Strands | Off |
| Output Prefix | `bgc_losatp_groups` |

Match the Interactive SVG Gallery presentation with these display settings:

| Control | Value |
|---|---|
| Palette | Orange |
| Override File / Specific Table | Bundled BGC color tables |
| Show Labels | First record |
| Priority File (TSV) | `CDS<TAB>gene` qualifier-priority table |
| Label Font Size | `18` |
| Label Placement / Rotation | Above feature / `45` |
| Feature Height | `75` |
| Block Stroke Width / Line Stroke Width | `2` / `2` |
| Show Coordinate Scale (Linear) | On |
| Linear scale style | Ruler |
| Axis Stroke Width | `5` |
| Lock Definition Column | On |
| Definition name | `20` px, Bold |
| Definition subtitle | `20` px, Normal |
| Definition accession / length | `20` px, Normal |

The first record therefore carries readable CDS `gene` labels; the remaining
four records stay unlabeled. Under **Title & Legend**, set the four visible
Definition line sizes to `20`, choose **Bold** only for **Name / Species**, and
leave the other lines at **Normal**. Fit the complete final preview at **40%**
before capturing or exporting it.

![LOSATP Similarity groups mode selected with deterministic settings](../../images/t-gui-04/03-losatp-settings.png)

## Step 4: Run LOSATP

Select **Generate Diagram** again. The result contains 23 stable
groups and 77 displayed group links. The four displayed endpoint pairs are
`0708→0709`, `0709→0711`, `0711→0712`, and `0712→0713`.

There is no direct `BGC0000708→BGC0000713` ribbon. Proteins shared by the first
and last records are represented by the same group ID across the adjacent-link
chain. This is the canonical Similarity-groups presentation; it is not a
Pairwise comparison between only the first and last records.

## Step 5: Align every record to `og_1`

Select the `livE` CDS in `og_1` on the first record. Its feature popup includes
an **Align** action because the current result is in Similarity-groups mode.

![og_1 feature popup with the Align action](../../images/t-gui-04/04-align-og1.png)

Select **Align**. gbdraw regenerates the same 23-group comparison without
rerunning LOSATP and shifts each record so its `og_1` member shares one
x-coordinate. This is the alignment used by the Interactive SVG Gallery.

![Five whole BGC records aligned to similarity group og_1](../../images/t-gui-04/05-comparison-result.png)

Under **Raw LOSAT results**, name the first adjacent result
`bgc_losatp_groups.tsv` and select **Save Raw LOSAT TSV**. The file
contains 232 twelve-column rows. Select **SVG** after alignment to save
`bgc_losatp_groups.svg`.

## Step 6: Inspect a group

Select a comparison ribbon. The popup reports the group ID, display name,
member count, record coverage, RBH seeds, paths, and every member protein.

![LOSATP similarity-group popup with member details](../../images/t-gui-04/06-match-popup.png)

## Next steps

- [Create protein similarity groups as a focused task](../../HOW_TO/GUI/create-protein-similarity-groups.md)
- [Draw Collinear protein-match blocks](../../HOW_TO/GUI/draw-collinear-protein-blocks.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
- [Review comparison programs and result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
