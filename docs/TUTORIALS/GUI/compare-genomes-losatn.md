[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [Web app](../../HOW_TO/GUI/README.md)

# Compare two genomes with LOSATN in the browser

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/compare-genomes-losatn.md) | [Python workflow](../PYTHON/compare-genomes-losatn.md) |

Compare the complete Lambda and DE3 phage genomes with browser LOSATN. You will first draw both records without a comparison, then run a serial, one-thread `megablast` search and export its evidence with the finished SVG.

![Linear genome comparison with nucleotide similarity links](../../images/t-gui-03/04-comparison-result.png)

*Six nucleotide matches connect the complete `NC_001416.1` and `NC_042057.1` records.*

## What you'll need

Starting state: open a fresh gbdraw web app page with no session loaded and no
files selected.

Use the filenames below when you download or save each file. See [Get the
tutorial inputs](../../GETTING_TUTORIAL_DATA.md) for browser download steps and
the meaning of each file type.

| File | File type | Purpose |
| --- | --- | --- |
| [`NC_001416.gb` — NCBI `NC_001416.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1) | Authoritative download | Complete Lambda GenBank record; save as `NC_001416.gb` |
| [`NC_042057.1.gb` — NCBI `NC_042057.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_042057.1) | Authoritative download | Complete DE3 GenBank record; save as `NC_042057.1.gb` |
| [`lambda-de3.losatn.tsv`](../../../gbdraw/web/tutorial-data/lambda-de3-comparison/lambda-de3.losatn.tsv) | Reference result | Frozen six-row result used only to check the new browser run |
| `lambda-de3.losatn.tsv` | Generated | Raw LOSATN evidence saved in Step 4 |
| `lambda-de3-losatn.svg` | Generated | Finished static diagram saved in Step 4 |

The two authoritative GenBank files are the browser inputs. The repository TSV
is a frozen reference output your own run should match; you will produce it
yourself in Step 4 rather than uploading it.

If you download the reference result, keep it in a separate comparison folder:
the reference result and your generated TSV intentionally have the same exact
filename, `lambda-de3.losatn.tsv`.

## Step 1: Load both complete genomes

Select **Linear**. Under **Input Genomes**, keep **GenBank** and **No comparison** selected.

1. In the first **GenBank File** control, choose `NC_001416.gb`.
2. Select **Add Seq**.
3. In the second **GenBank File** control, choose `NC_042057.1.gb`.

Keep Lambda first and DE3 second. Leave both **Region (optional)** sections unchanged so the inputs remain the complete `NC_001416.1` (48,502 bp) and `NC_042057.1` (42,925 bp) records.

![Two complete GenBank records selected for a Linear diagram](../../images/t-gui-03/01-input-ready.png)

*The two green upload controls confirm the input order: Lambda first, DE3 second.*

## Step 2: Generate the map without comparison links

Select **Generate Diagram**. The first result should contain two annotated records and no links between them.

![Two complete records in a plain linear diagram](../../images/t-gui-03/02-first-diagram.png)

*Both record tracks are visible, with no links between them.*

## Step 3: Configure LOSATN

Return to **Input Genomes**, select **Run LOSAT**, and keep **LOSATN** selected. Set these values:

| Control | Value |
| --- | --- |
| Execution | Serial |
| Total threads | 1 |
| Parallel runs | 1 run |
| Threads per run | Fixed at 1 |
| Task | `megablast` |
| Pairwise Match Height | `120` |
| Raw LOSAT filename | `lambda-de3.losatn.tsv` |
| Output Prefix | `lambda-de3-losatn` |

![LOSATN selected with deterministic single-thread settings](../../images/t-gui-03/03-losatn-settings.png)

*Serial execution, one total thread, one run, and `megablast` produce the deterministic six-link result shown throughout this tutorial.*

## Step 4: Run LOSATN and download the evidence

Select **Generate Diagram** again. LOSATN runs in the browser, and the result should show six links in the enlarged comparison corridor. These links are high-identity nucleotide alignments; each link joins the query and subject intervals recorded in one TSV row.

![Linear genome comparison with nucleotide similarity links](../../images/t-gui-03/04-comparison-result.png)

*The longest match covers 21,232 aligned bases at 99.981% identity.*

In the comparison between sequence 1 and sequence 2, select **Save Raw LOSAT TSV**. The browser saves `lambda-de3.losatn.tsv`. Then select **SVG** in the **Result Preview** toolbar to save `lambda-de3-losatn.svg`.

The TSV should contain six tab-separated rows. Every query interval falls within 1–48,502, and every subject interval falls within 1–42,925.

## Step 5: Inspect one nucleotide match

Select the longest comparison ribbon in the preview. The **Pairwise match** popup identifies both records and reports the intervals, identity, alignment length, E-value, bit score, mismatches, and gap opens.

![LOSATN match details popup in the result preview](../../images/t-gui-03/05-match-popup.png)

*The first match connects Lambda 1..21231 to DE3 20081..41311 and reports 99.981% identity.*

## Next steps

- [Use uploaded comparison results](../../HOW_TO/GUI/use-uploaded-blast-results.md)
- [Arrange multiple Linear records](../../HOW_TO/GUI/arrange-linear-records-regions-and-orientation.md)
- [Choose an export format](../../REFERENCE/output-formats-and-export.md)
- [Review input and TSV schemas](../../REFERENCE/input-formats-and-tsv-schemas.md)
