[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [Web app technical documentation](../../REFERENCE/web-app.md)

# Create and export your first circular genome diagram

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/first-circular-genome-diagram.md) | [Python workflow](../PYTHON/first-genome-diagram.md) |

All three workflows build the same labeled human mitochondrial figure. Only
the interface changes.

Create a labeled Circular SVG from the accession-pinned human mitochondrial
GenBank record. You will generate a working diagram in Step 2, then add the
publication label and external feature labels before exporting it.

![Finished circular human mitochondrial genome diagram with external labels and a right legend](../../images/t-gui-01/04-finished-diagram.png)

*The finished diagram identifies the record, keeps both GC tracks, places feature labels outside the circle, and puts the legend on the right.*

## What you'll need

Starting state: open a fresh gbdraw web app page with no session loaded and no
files selected.

Use the filenames below when you download or save each file. See [Get the
tutorial inputs](../../GETTING_TUTORIAL_DATA.md) for browser download steps and
the meaning of each file type.

| File | File type | Purpose |
| --- | --- | --- |
| [`HmmtDNA.gbk` — NCBI `NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) | Authoritative download | Complete human mitochondrial GenBank record; save as `HmmtDNA.gbk` |
| [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Support download | CDS label-priority rule |
| `human_mitochondrion.svg` | Generated | Finished static diagram saved in Step 5 |

## Step 1: Load the NCBI mitochondrial genome

Select **Circular** at the top of the app and keep **GenBank** selected under **Input Genomes**. In **GenBank/DDBJ File**, choose `HmmtDNA.gbk`.

The uploader should show `HmmtDNA.gbk` in green.

![Circular GenBank input showing HmmtDNA.gbk selected](../../images/t-gui-01/01-input-ready.png)

*Confirm that Circular and GenBank are selected and that the uploader names `HmmtDNA.gbk`.*

## Step 2: Generate the first diagram

Select **Generate Diagram** without changing the advanced settings. When processing finishes, **Result Preview** displays the first Circular map.

![First circular human mitochondrial genome diagram](../../images/t-gui-01/02-first-diagram.png)

*The first result identifies `NC_012920.1` and includes feature rings, GC content, GC skew, and coordinate ticks.*

## Step 3: Add a publication label

Under **Basic**, enter these values:

| Control | Value |
| --- | --- |
| Output Prefix | `human_mitochondrion` |
| Species | `<i>Homo sapiens</i>` |

Select **Generate Diagram** again. The center label should show *Homo sapiens* in italics.

![Circular preview labeled Homo sapiens](../../images/t-gui-01/03-publication-label.png)

*The Basic controls retain the output prefix and species markup, while the preview renders the species name in italics.*

## Step 4: Make the feature map easier to read

Set the final layout values below. **Track Preset** and the three checkboxes are under **Layout**. **Label Mode** and **Priority File (TSV)** are under **Labels**, and **Legend Position** is under **Title & Legend**.

| Control | Value |
| --- | --- |
| Track Preset | Middle |
| Separate Strands | On |
| Hide GC Content | Off |
| Hide GC Skew | Off |
| Label Mode | Out |
| Priority File (TSV) | `cds_gene_qualifier_priority.tsv` |
| Legend Position | Right |

![Circular layout controls set to Middle, Labels Out, and Legend Right](../../images/t-gui-01/04-layout-settings.png)

*The visible Layout controls show Middle with separate strands enabled and both GC tracks retained. Keep Labels at Out and the legend at Right as listed above.*

Select **Generate Diagram**. The completed map adds external feature labels, uses the `gene` qualifier for every CDS label, and retains the right-side legend.

![Finished circular human mitochondrial genome diagram with external labels and a right legend](../../images/t-gui-01/04-finished-diagram.png)

*The finished preview shows the labeled mitochondrial map and the six-entry legend on its right.*

## Step 5: Export the SVG

In the **Result Preview** toolbar, select **SVG**.

![SVG download button below the finished result preview](../../images/t-gui-01/05-export-svg.png)

*Use SVG for the static publication figure. The browser saves `human_mitochondrion.svg`.*

## Next steps

- [Create a Linear layout](first-linear-genome-diagram.md)
- [Review feature-presentation rules](../../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation)
- [Compare genomes](compare-genomes-losatn.md)
- [Save and restore an interactive session](../../REFERENCE/session-and-request-compatibility.md)
