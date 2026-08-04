# Add circular comparison rings from precomputed results

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| **This page** | [Command-line workflow](../CLI/add-precomputed-circular-comparison-rings.md) | [Python workflow](../PYTHON/add-precomputed-circular-comparison-rings.md) |

Add three translated-nucleotide comparison rings around the complete human
mitochondrial genome. The searches are already frozen as TLOSATX outfmt 6
tables, so this project focuses on reference direction, ring order, filtering,
inspection, and sequence export.

![Human mitochondrial map with three precomputed TLOSATX rings](../../images/t-gui-06/04-ring-result.png)

## What you'll need

Use `HmmtDNA.gbk` from
`gbdraw/web/tutorial-data/human-mitochondrion/` and these files from
`gbdraw/web/tutorial-data/metazoan-mitochondria-comparison/`:

| Ring | TLOSATX table | Comparison FASTA |
| ---: | --- | --- |
| 1 | `danio-human.tlosatx.tsv` | `NC_002333.2.fna` |
| 2 | `drosophila-human.tlosatx.tsv` | `NC_024511.2.fna` |
| 3 | `caenorhabditis-human.tlosatx.tsv` | `NC_001328.1.fna` |

## Step 1: Load the human reference

Select **Circular** and **GenBank**, choose `HmmtDNA.gbk`, and set:

| Control | Value |
| --- | --- |
| Output Prefix | `precomputed_circular_rings` |
| Species | `<i>Homo sapiens</i>` |
| Track Preset | Middle |
| Separate Strands | Off |

![Complete human mitochondrial reference ready for Circular comparison](../../images/t-gui-06/01-input-ready.png)

## Step 2: Generate the reference map

Select **Generate Diagram** before adding comparison evidence. The first result
identifies `NC_012920.1`, contains 37 feature elements, and has no comparison
rings.

![First human mitochondrial diagram without similarity rings](../../images/t-gui-06/02-first-diagram.png)

## Step 3: Add the three precomputed comparisons

Open **Pairwise Comparisons**, select **Upload BLAST**, and choose all three TSV
files together. Set **Reference side** to **Subject**: every table uses a
comparison genome as query and human mtDNA as subject.

Attach the matching comparison FASTA to each table, then use:

| Control | Value |
| --- | ---: |
| Bitscore | 50 |
| E-value | `1e-5` |
| Minimum identity | 40 |
| Minimum length | 50 |
| Ring width | 18 |
| Ring gap | 4 |

Keep the ring order Danio, Drosophila, and Caenorhabditis. Set **Label Mode**
to **Out**, the title to
`Precomputed TLOSATX rings around Homo sapiens mtDNA`, and the legend to the
right.

![Three TLOSATX tables and companion FASTA files configured as Circular rings](../../images/t-gui-06/03-ring-settings.png)

## Step 4: Generate and read the rings

Select **Generate Diagram**. The filters retain 106 HSPs across the three
rings. Colors and legend order identify the comparison source; positions
around the circle are always coordinates on the human subject record.

## Step 5: Inspect an HSP and export its spans

Select a colored HSP. The popup shows its endpoints, identity, alignment
length, reference side, and source ring. Use **Both spans FASTA** to save the
reference and comparison intervals as `circular_hsp_spans.fasta`.

![Circular HSP popup with reference and comparison span exports](../../images/t-gui-06/05-hsp-popup.png)

Select **SVG** to save `precomputed_circular_rings.svg`.

## What you built

You combined one complete reference with three fixed evidence tables, kept the
query/subject direction explicit, and exported both a publication SVG and the
two sequences behind an inspected match.

## Next steps

- [Add Circular similarity rings as a focused task](../../HOW_TO/GUI/add-circular-similarity-rings.md)
- [Use uploaded BLAST results](../../HOW_TO/GUI/use-uploaded-blast-results.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
- [Review comparison result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
