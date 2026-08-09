[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [Web app](../../HOW_TO/GUI/README.md)

# Create an interactive figure and reproduce it from a saved session

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/create-and-resume-an-interactive-figure.md) | [Python workflow](../PYTHON/create-and-resume-an-interactive-figure.md) |

Create a finished human mitochondrial map, export an offline Interactive SVG,
inspect `COX1`, save the project, and reproduce the same figure in a fresh
browser context.

![Reloaded human mitochondrial figure reproduced from a saved session](../../images/t-gui-09/06-reloaded-result.png)

## What you'll need

Starting state: open a fresh gbdraw web app page with no session loaded and no
files selected. Step 6 starts from a second fresh page and then changes it to a
loaded-session state.

Use the filenames below when you download or save each file. See [Get the
tutorial inputs](../../GETTING_TUTORIAL_DATA.md) for browser download steps and
the meaning of each file type.

| File | File type | Purpose |
| --- | --- | --- |
| [`HmmtDNA.gbk` — NCBI `NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) | Authoritative download | Complete human mitochondrial GenBank record; save as `HmmtDNA.gbk` |
| [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Support download | CDS `gene` label priority |
| `interactive_human_mitochondrion.interactive.svg` | Generated | Offline Interactive SVG saved in Step 3 |
| `interactive_handoff.gbdraw-session.json.gz` | Generated | Project session saved in Step 5 and loaded in Step 6 |
| `restored_interactive_figure.svg` | Generated | Static SVG regenerated from the loaded session in Step 6 |

## Step 1: Load the record

Select **Circular** and **GenBank**, choose `HmmtDNA.gbk`, and set **Output
Prefix** to `interactive_human_mitochondrion`.

![Human mitochondrial GenBank record ready for the interactive project](../../images/t-gui-09/01-input-ready.png)

## Step 2: Generate the finished map

Set Species to `<i>Homo sapiens</i>`, use the Middle track preset, separate the
strands, keep GC content and GC skew, and place labels outside. Under
**Labels**, load `cds_gene_qualifier_priority.tsv` as **Priority File (TSV)**.
Put the legend on the right, then select **Generate Diagram**.

![Finished human mitochondrial map before interactive export](../../images/t-gui-09/02-first-diagram.png)

The result contains the complete `NC_012920.1` record and 37 rendered feature
elements.

## Step 3: Export the offline Interactive SVG

In the Result Preview toolbar, select **Interactive SVG**. The browser saves
`interactive_human_mitochondrion.interactive.svg`.

![Interactive SVG export action above the finished preview](../../images/t-gui-09/03-interactive-export.png)

This is a self-contained SVG: its feature metadata, search interface, popup
logic, and runtime assets travel with the file. It does not need the gbdraw web
app or a network connection when opened later.

## Step 4: Inspect COX1 before handoff

In **Search features**, choose **Qualifier value**, enter `gene` as the
qualifier key, search for `COX1`, and open the active feature. Confirm its
location in the feature popup.

![COX1 search result opened in the feature-details popup](../../images/t-gui-09/04-feature-search.png)

Close the popup and clear the search.

## Step 5: Save the project session

Select **Save Session**, enter `interactive_handoff`, and save
`interactive_handoff.gbdraw-session.json.gz`.

![Finished project after downloading the interactive handoff session](../../images/t-gui-09/05-session-download.png)

The Interactive SVG is the reader-facing artifact. The session is the
reproducible working state: it preserves the inputs, settings, result, and
compatible render request needed to continue editing.

## Step 6: Reproduce the figure in a fresh context

Open a fresh gbdraw page with no files selected, select **Load Session**, and
choose `interactive_handoff.gbdraw-session.json.gz`. The app is now in the
loaded-session state. After the result is restored, change **Output Prefix** to
`restored_interactive_figure`, select **Generate Diagram**, and export **SVG**.

The saved file is `restored_interactive_figure.svg`. Its record, feature IDs,
texts, labels, track groups, and placement match the original figure; only
subpixel browser font measurement may vary by less than one pixel.

## Next steps

- [Inspect and edit a diagram](../../HOW_TO/GUI/inspect-and-edit-a-diagram.md)
- [Save, restore, undo, and reproduce work](../../HOW_TO/GUI/save-restore-undo-and-reproduce-work.md)
- [Export publication and interactive figures](../../HOW_TO/GUI/export-publication-and-interactive-figures.md)
- [Review Interactive SVG semantic hooks](../../REFERENCE/interactive-svg-and-semantic-hooks.md)
