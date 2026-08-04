[Documentation home](../../DOCS.md) | [How-to guides](../README.md) | [GUI how-to guides](./README.md)

# How to draw collinear protein-match blocks with LOSATP

Use **Collinear blocks** when conserved gene order is part of the comparison.
This workflow uses the same five whole, Linear BGC records as the Similarity
groups example, but produces ordered blocks instead of one ribbon per displayed
group relation.

## Configure Collinear blocks

Upload `BGC0000708`, `BGC0000709`, `BGC0000711`, `BGC0000712`, and
`BGC0000713` in that order. Do not crop a record. Keep the first four source
orientations and turn on **Reverse complement** only for `BGC0000713`, matching
the Gallery layout.

Select **Run LOSAT**, **LOSATP**, and **Collinear blocks**, then set:

| Control | Value |
|---|---|
| Execution | Serial |
| Total threads | `1` |
| Parallel runs | `1 run` |
| Max unit gap | `1` |
| Min block genes | `2` |
| Evidence scope | Adjacent pairs |
| Diagonal drift | `1` |
| Merge conflicts | `0` |
| Color mode | Orientation + identity |
| Bitscore / E-value | `50` / `0.01` |
| Minimum identity / length | `30` / `0` |
| Output Prefix | `bgc_collinear_blocks` |

Match the Interactive SVG Gallery presentation as well: choose the Orange
palette and bundled BGC color tables; set **Show Labels** to **First record**,
load the `CDS<TAB>gene` qualifier-priority table, and use label font size `18`,
placement **Above feature**, and rotation `45`. Set **Feature Height** to `75`,
**Block Stroke Width** and **Line Stroke Width** to `2`, turn on
**Show Coordinate Scale (Linear)** with the **Ruler** style, and set
**Axis Stroke Width** to `5`.

![LOSATP Collinear block settings](../../images/h-gui-08/collinear-settings.png)

## Generate and inspect the blocks

Select **Generate Diagram**. The qualified result contains seven blocks across
adjacent record pairs. Their anchor counts are `13`, `3`, `21`, `2`, `15`,
`13`, and `2`; orientations are plus, minus, plus, plus, plus, minus, and
minus in source-record coordinates. Because `BGC0000713` is reverse-complemented
for display, the last two blocks appear as plus-orientation blocks in the
Gallery layout. Fit the complete preview at **40%** before capture or export. Only
`BGC0000708` carries labels, and its CDS labels use `gene` values rather than
`product` descriptions.

![Five BGC records connected by collinear protein-match blocks](../../images/h-gui-08/collinear-result.png)

Select a block to open its details. The popup reports its ID, orientation,
identity color mode, span, score, and anchor count. It also lists the local
collinear groups covered by the block.

![Collinear block details popup](../../images/h-gui-08/block-popup.png)

Select the `og_1` row in the popup, then use **Download all member amino-acid
FASTA**. The qualified `collinear_members.fasta` contains these five proteins:
`CAG38695.1`, `CAF33310.1`, `CAH58688.1`, `CAF32372.1`, and `CAG34720.1`.
The executable capture also checks every downloaded amino-acid sequence against
the `translation` of the same `protein_id` in its raw GenBank record.
Select **SVG** to save `bgc_collinear_blocks.svg`.

## Troubleshooting

- **Many short links appear:** confirm that **Collinear blocks**, not
  **Pairwise**, is selected.
- **Orientations differ from this guide:** keep the first four source records
  unchanged and reverse only `BGC0000713`.
- **No member FASTA button appears:** select one local group row inside the
  block popup first.
- **A circular diagram appears:** switch back to **Linear**. These BGC entries
  are linear region records.

## Related guides

- [Create protein similarity groups](./create-protein-similarity-groups.md)
- [Compare annotated proteins Tutorial](../../TUTORIALS/GUI/compare-proteins-losatp.md)
- [Comparison programs, thresholds, and result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
