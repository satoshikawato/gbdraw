[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [Web app technical documentation](../../REFERENCE/web-app.md)

# Create and export your first linear genome diagram

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/first-linear-genome-diagram.md) | [Python workflow](../PYTHON/first-linear-genome-diagram.md) |

Create a labeled Linear SVG from the complete, accession-pinned Lambda phage
GenBank record. You will generate a working map in Step 2, then add concise
gene labels and a coordinate ruler before exporting it.

![Finished linear Lambda genome diagram with labels and ruler](../../images/t-gui-02/04-finished-diagram.png)

*The finished map shows the whole 48,502 bp record, both strands, short gene labels, a ruler, and the CDS legend on the left.*

## What you'll need

Starting state: open a fresh gbdraw web app page with no session loaded and no
files selected.

Use the filenames below when you download or save each file. See [Get the
tutorial inputs](../../GETTING_TUTORIAL_DATA.md) for browser download steps and
the meaning of each file type.

| File | File type | Purpose |
| --- | --- | --- |
| [`NC_001416.gb` — NCBI `NC_001416.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1) | Authoritative download | Complete Lambda GenBank record; save as `NC_001416.gb` |
| [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Support download | CDS label-priority rule |
| `lambda_linear.svg` | Generated | Finished static diagram saved in Step 4 |

The label rule tells gbdraw to use short `gene` values for CDS labels. The
sequence comes from NCBI; the rule is repository-maintained support data.

## Step 1: Load the NCBI Lambda genome

Select **Linear** at the top of the app. A fresh Linear page starts with **No
comparison**. Under **Input Genomes**, keep **GenBank** selected. The **GenBank
File** uploader is the first control in **Sequence 1**; choose
`NC_001416.gb`. Leave **Record options** closed.

The uploader should show `NC_001416.gb` in green. Keep one input row only; this Tutorial uses the complete `NC_001416.1` record without cropping or splitting it.

![Linear GenBank input showing NC_001416.gb selected](../../images/t-gui-02/01-input-ready.png)

*Confirm Linear, GenBank, the `NC_001416.gb` upload, and the separate **Current:
No comparison** status before generating.*

## Step 2: Generate the first diagram

Select **Generate Diagram** without changing the presentation settings. When processing finishes, **Result Preview** displays the first Linear map.

![First linear Lambda genome diagram](../../images/t-gui-02/02-first-diagram.png)

*The first result identifies `NC_001416.1` and `48,502 bp`, draws all 73 CDS features, and keeps the two strand lanes separate.*

## Step 3: Add concise labels and a ruler

Set **Output Prefix** under **Basic**. **Generate Diagram** follows **Basic** in
the DOM. Continue past it to the **Layout**, **Labels**, **Axis & Scale**, and
**Title & Legend** sections for the remaining values below. Leave the closed
**Advanced comparison and layout** disclosure unchanged.

| Control | Value |
| --- | --- |
| Output Prefix | `lambda_linear` |
| Track Layout | Features on axis |
| Separate Strands | On |
| Show Labels | All Records |
| Priority File (TSV) | `cds_gene_qualifier_priority.tsv` |
| Show Coordinate Scale | On |
| Scale Style | Ruler (Ticks) |
| Legend Position | Left |

![Linear layout controls configured for labels and a ruler](../../images/t-gui-02/03-layout-settings.png)

*The visible Axis & Scale controls show the coordinate scale and Ruler (Ticks). Keep the label, strand, and legend values exactly as listed in the table.*

## Step 4: Regenerate and export the SVG

Select **Generate Diagram** again. The completed map should include short labels such as `A`, `B`, `J`, and `int`. Its ruler spans the complete Lambda record, and the CDS legend appears on the left.

![Finished linear Lambda genome diagram with labels and ruler](../../images/t-gui-02/04-finished-diagram.png)

*Check the full 48,502 bp map before export. The ruler, labels, record definition, two strand lanes, and left CDS legend should all be visible.*

In the **Result Preview** toolbar, select **SVG**. The browser saves `lambda_linear.svg`.

## Next steps

- [Review record selection and layout](../../REFERENCE/web-app.md#record-selection-and-layout)
- [Review feature-presentation rules](../../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
- [Compare genomes](compare-genomes-losatn.md)
- [Choose another output format](../../REFERENCE/output-formats-and-export.md)
