[Documentation home](../../DOCS.md) | [How-to guides](../README.md) | [GUI how-to guides](./README.md)

# How to draw collinear protein-match blocks with LOSATP

Use **Collinear blocks** when conserved gene order is part of the comparison.
This workflow reuses All-vs-all LOSATP evidence from five complete
Hepatoplasmataceae genomes and reduces compatible protein groups to ordered
blocks. The genomes remain complete records; the Linear view only arranges them
as comparison rows.

## Load the Hepatoplasmataceae project

Load
`gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz`.
The project contains these records in display order:

| Order | Record | Organism | Length |
|---:|---|---|---:|
| 1 | `AP027078.1` | *Candidatus Tyloplasma litorale* | 615,622 bp |
| 2 | `AP027131.1` | *Candidatus Hepatoplasma vulgare* | 662,108 bp |
| 3 | `AP027133.1` | *Candidatus Hepatoplasma scabrum* | 606,194 bp |
| 4 | `AP027132.1` | *Candidatus Hepatoplasma crinochetorum* | 643,039 bp |
| 5 | `NZ_CP006932.1` | *Candidatus Hepatoplasma crinochetorum* Av | 657,101 bp |

All five source records are complete and naturally circular. The loaded project
uses a Linear presentation so homologous intervals can be connected between
rows. Its cache contains all 25 directional and self LOSATP results from the
five-record All-vs-all search.

## Configure Collinear blocks

Select **Run LOSAT**, **LOSATP**, and **Collinear blocks**, then set:

| Control | Value |
|---|---|
| Execution | Threaded |
| Total threads | `32` |
| Parallel runs | `1 run` |
| Threads per run | `32` |
| Max unit gap | `0` |
| Min block genes | `1` |
| Evidence scope | All records |
| Diagonal drift | `0` |
| Merge conflicts | `1` |
| Color mode | Orientation + identity |
| Bitscore / E-value | `50` / `0.01` |
| Minimum identity / length | `0` / `0` |
| Pairwise Match Style | Curve |
| Output Prefix | `hepatoplasmataceae_collinear` |

Keep the project presentation at **Middle**, with **Separate Strands**, **Align
to Center**, **Show GC Content**, and **Show GC Skew** on. Use the **Ajisai**
palette, a coordinate **Ruler**, curved ribbons, a top title, and a right legend.

![LOSATP Collinear controls using all Hepatoplasmataceae record pairs](../../images/h-gui-08/collinear-settings.png)

## Generate and inspect the blocks

Select **Generate Diagram**. All 25 cached directional and self results are
reused, so this presentation change has no LOSATP cache miss and launches no new
worker job. Block construction uses evidence from every record pair, while the
finished ribbons connect adjacent display rows to keep the figure readable.
Fit the complete preview at **30%**.

![Five complete Hepatoplasmataceae genomes connected by collinear protein-match blocks](../../images/h-gui-08/collinear-result.png)

Select a block to open its details. The popup reports its ID, orientation,
query and subject intervals, anchors, covered protein groups, and both envelope
sequences.

![Hepatoplasmataceae Collinear block details popup](../../images/h-gui-08/block-popup.png)

Use **Both spans FASTA** to save `collinear_members.fasta`. It should
contain two non-empty nucleotide envelope sequences, one matching each
block's source genome. Select **SVG** to save
`hepatoplasmataceae_collinear.svg`.

## Troubleshooting

- **Many short links appear:** confirm that **Collinear blocks**, not
  **Pairwise**, is selected.
- **LOSATP starts running again:** load the documented session and confirm that
  all 25 cached results are present before generating.
- **Only adjacent evidence is used:** change **Evidence scope** to **All
  records**.
- **A Circular diagram appears:** switch back to **Linear** for the comparison
  view. This does not change the complete circular source records.

## Related guides

- [Create protein similarity groups](./create-protein-similarity-groups.md)
- [Compare annotated proteins Tutorial](../../TUTORIALS/GUI/compare-proteins-losatp.md)
- [Comparison programs, thresholds, and result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
