[Documentation home](../../DOCS.md) | [How-to guides](../README.md) | [GUI how-to guides](./README.md)

# How to add circular sequence-similarity rings

Circular comparison rings require a complete record whose source topology is
circular. This recipe uses complete human mitochondrial DNA as the displayed
reference and compares it with three other complete, naturally circular
mitochondrial genomes. No linear phage genome or partial region is drawn as a
circle.

## Before you start

The capture uses these four RefSeq records:

| Role | Organism | RefSeq | Length | Topology |
|---|---|---|---:|---|
| Displayed reference | *Homo sapiens* | `NC_012920.1` | 16,569 bp | Circular |
| Ring 1 | *Danio rerio* | `NC_002333.2` | 16,596 bp | Circular |
| Ring 2 | *Drosophila melanogaster* | `NC_024511.2` | 19,524 bp | Circular |
| Ring 3 | *Caenorhabditis elegans* | `NC_001328.1` | 13,794 bp | Circular |

The packaged records are
[`HmmtDNA.gbk`](../../../gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk),
[`NC_002333.2.gb`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-four/NC_002333.2.gb),
[`NC_024511.2.gb`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-four/NC_024511.2.gb),
and
[`NC_001328.1.gb`](../../../gbdraw/web/tutorial-data/metazoan-mitochondria-four/NC_001328.1.gb).

**Run LOSAT** accepts comparison FASTA files. The capture converts each of the
last three GenBank records to one whole-record FASTA entry and checks that its
accession, length, and sequence are unchanged. When preparing files manually,
export the complete record. Do not trim it or split it into mock contigs.

## Load the complete circular reference

1. Select **Circular** and **GenBank**.
2. Upload `HmmtDNA.gbk` under **GenBank/DDBJ File**.
3. Set **Species** to `<i>Homo sapiens</i>`.
4. Set **Output Prefix** to `circular_similarity_rings`.
5. Open **Layout** and turn off **Separate Strands**.

Open **Labels**, choose **Out**, and upload
[`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv)
under **Priority File (TSV)**. The file contains `CDS` followed by `gene`, so
the 13 mitochondrial CDS labels use `gene`. Product descriptions are not used
as CDS labels.

## Add three translated-nucleotide rings

Open **Pairwise Comparisons**, select **Run LOSAT**, and choose **TLOSATX**.
Set the reference gencode to `2`, then add the three complete comparison FASTA
files in this order:

| Ring | Label | Subject gencode |
|---:|---|---:|
| 1 | `Danio rerio (NC_002333.2)` | `2` |
| 2 | `Drosophila melanogaster (NC_024511.2)` | `5` |
| 3 | `Caenorhabditis elegans (NC_001328.1)` | `5` |

Keep the three automatically assigned series colors. They are distinct and
remain tied to the documented ring order.

Use these comparison filters:

| Control | Value |
|---|---:|
| Bitscore | `50` |
| E-value | `1e-5` |
| Minimum Identity | `40` |
| Minimum Length | `50` |
| Ring Width | `18` |
| Ring Gap | `4` |

Open **Title & Legend**, set the title to
`Homo sapiens mtDNA (NC_012920.1) similarity rings`, place it at the top, and
put the legend on the right.

![Circular comparison settings for three complete mtDNA sources](../../images/h-gui-06/ring-settings.png)

## Generate and inspect one HSP

Select **Generate Diagram**, then zoom to 50%.

![Complete human mtDNA with three sequence-similarity rings](../../images/h-gui-06/ring-result.png)

The displayed human record is the TLOSATX subject. Each ring keeps its source
accession and has at least one retained HSP. The capture records the union of
all HSP intervals on `NC_012920.1` and rejects any interval outside 1 to
16,569.

Select one colored HSP to open **Homology ring match**. The popup shows both
record IDs, aligned coordinates, orientation, identity, and alignment length.
Use the combined **FASTA** button to download both nucleotide spans. The
capture stores that proof as `comparison_spans.fasta` and verifies that it has
two non-empty records.

![Similarity-ring HSP details popup](../../images/h-gui-06/hsp-popup.png)

Select **SVG** to save `circular_similarity_rings.svg`. The SVG check requires
three rings in the documented order, the complete 16,569 bp human reference,
all four accessions and species labels, 13 CDS `gene` labels, and zero CDS
`product` labels.

## Troubleshooting

- **A comparison FASTA produces no ring:** check that it contains one complete
  mitochondrial record and use the matching mitochondrial genetic code.
- **The reference side is wrong:** with browser LOSAT, the displayed circular
  record is the subject and each comparison FASTA is a query.
- **Labels show product descriptions:** upload
  `cds_gene_qualifier_priority.tsv` and generate again.
- **The source record is linear or partial:** use a Linear diagram. Circular
  rings do not make linear or partial input circular.

## Related guides

- [Arrange multiple complete Circular records](./arrange-multiple-circular-records.md)
- [Use TLOSATX for translated nucleotide comparisons](./use-tlosatx.md)
- [Tutorial fixture provenance](../../REFERENCE/tutorial-fixture-provenance.md)
