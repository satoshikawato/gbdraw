# Build an annotated chloroplast map with custom tracks

Build a Circular map of the complete tobacco plastome, mark its LSC, IRb, SSC,
and IRa regions, and replace the default GC-skew track with AT skew. You will
make a working diagram in Step 2 before arranging the final track stack.

![Finished tobacco plastome map with structural regions and custom tracks](../../images/t-gui-05/05-finished-diagram.png)

## What you'll need

- `gbdraw/web/tutorial-data/tobacco-plastome-regions/NC_001879.gbk`
- `gbdraw/web/tutorial-data/tobacco-plastome-regions/nicotiana-tabacum-regions.tsv`

The GenBank file is the complete 155,943 bp `NC_001879.2` record. The TSV
contains four source-coordinate annotations; it does not modify the sequence.

## Step 1: Load the plastome

Select **Circular** and **GenBank**, then choose `NC_001879.gbk` in
**GenBank/DDBJ File**. Set **Output Prefix** to
`annotated_chloroplast_map`.

![Complete tobacco plastome selected as the Circular input](../../images/t-gui-05/01-input-ready.png)

## Step 2: Generate the first diagram

Select **Generate Diagram**. Confirm that the first map identifies
`NC_001879.2` and contains its complete feature set. This is the visible
baseline; the four structural regions have not been added yet.

![First complete tobacco plastome diagram before region annotations](../../images/t-gui-05/02-first-diagram.png)

## Step 3: Import the structural regions

Open **Region Annotations** and use **Import TSV** to choose
`nicotiana-tabacum-regions.tsv`. Set the legend label to
`Plastome structural regions` and use these lanes:

| Region | Range | Lane |
| --- | ---: | ---: |
| LSC | 1–86,686 | 0 |
| IRb | 86,687–112,029 | 1 |
| SSC | 112,030–130,600 | 0 |
| IRa | 130,601–155,943 | 1 |

![Four tobacco plastome regions imported from the annotation table](../../images/t-gui-05/03-annotation-table.png)

## Step 4: Build the custom track stack

Open **Custom Track Slots** and turn on **Use custom stack**. In the existing
skew slot, change **Dinucleotide** to `AT` and the legend label to `AT skew`.
Add an **Annotations** renderer, select `plastome_regions`, move it outside the
axis, and turn on its labels.

Set **Label Mode** to **None** so the four structural labels remain readable
instead of competing with every feature label. Under **Title & Legend**, use
the title `Complete Nicotiana tabacum plastome regions` and place the legend on
the right.

![Custom Circular track stack with the structural annotation slot outside the axis](../../images/t-gui-05/04-track-settings.png)

## Step 5: Generate and export the finished map

Select **Generate Diagram** again. The finished figure contains the four
bracketed regions, GC content, positive and negative AT skew, the complete
feature map, title, and legend. Select **SVG** to save
`annotated_chloroplast_map.svg`.

## What you built

You produced a complete plastome figure without editing the GenBank record.
The biological sequence remains the source of features, while the annotation
table and custom track stack own the structural-region and quantitative
presentation.

## Next steps

- [Add region annotations and custom track slots](../../HOW_TO/GUI/add-region-annotations-and-track-slots.md)
- [Add depth, GC content, and skew tracks](../../HOW_TO/GUI/add-depth-gc-and-skew-tracks.md)
- [Understand tracks, axes, and layout](../../EXPLANATION/understand-tracks-axes-and-layout.md)
- [Review track and annotation schemas](../../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
