# Find conserved gene order from all-vs-all LOSATP evidence

Start from a completed All-vs-all LOSATP Similarity-groups project and reduce
the same protein evidence to Collinear blocks across five complete
Hepatoplasmataceae genomes. The cached searches are reused; switching the
presentation does not rerun LOSATP.

![Hepatoplasmataceae diagram with Collinear blocks from all-record evidence](../../images/t-gui-08/04-collinear-result.png)

## What you'll need

Use the bundled Gallery session
`gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz`.
It contains the complete records in this order:

1. `AP027078.1`
2. `AP027131.1`
3. `AP027133.1`
4. `AP027132.1`
5. `NZ_CP006932.1`

The session was created with a 32-thread LOSATP run and preserves 25 cached
directional and self results for the five-record All-vs-all analysis.

## Step 1: Load the All-vs-all project

Select **Load Session** and choose
`hepatoplasmataceae_orthogroup.gbdraw-session.json.gz`. Accept the successful
load message. Confirm that Linear mode, the five GenBank files, **Run LOSAT**,
**LOSATP**, and **Similarity groups** are restored.

![Five Hepatoplasmataceae genomes restored from an All-vs-all LOSATP session](../../images/t-gui-08/01-input-ready.png)

## Step 2: Inspect the Similarity-groups result

The restored preview is already a visible result. Its links show group
membership between adjacent display rows, while group construction used the
complete All-vs-all evidence cache.

![Hepatoplasmataceae All-vs-all LOSATP Similarity-groups diagram](../../images/t-gui-08/02-first-diagram.png)

## Step 3: Change the reduction to Collinear blocks

Keep **Run LOSAT** and **LOSATP**, then set:

| Control | Value |
| --- | --- |
| Execution | Threaded |
| Total threads | 32 |
| Parallel runs | 1 run |
| Threads per run | 32 |
| blastp mode | Collinear blocks |
| Max unit gap | 0 |
| Min block genes | 1 |
| Color mode | Orientation + identity |
| Evidence scope | All records |
| Diagonal drift | 0 |
| Merge conflicts | 1 |
| Bitscore / E-value | 50 / `0.01` |
| Minimum identity / length | 0 / 0 |

Set **Output Prefix** to `losatp_collinear`. Keep the Middle layout, centered
records, separate strands, GC content and GC skew, ruler, ajisai palette,
curved ribbons, top title, and right legend.

![LOSATP Collinear controls using all Hepatoplasmataceae record pairs](../../images/t-gui-08/03-collinear-settings.png)

## Step 4: Generate the Collinear view

Select **Generate Diagram**. The app requests all ten unique record pairs, but
all 25 directional/self cache entries are hits: there are no search misses and
no LOSATP worker calls. Collinear reduction uses evidence from every record,
while ribbons remain between adjacent display rows so the Linear figure stays
readable.

Blue and red families distinguish orientation; intensity carries average
identity within each family.

## Step 5: Inspect and export a block

Select a ribbon. The popup reports the block ID, query and subject intervals,
orientation, covered Similarity groups, anchors, and envelope sequences.

![Collinear block popup with coordinates, anchors, and FASTA actions](../../images/t-gui-08/05-block-popup.png)

Use **Both spans FASTA** to save `collinear_members.fasta`, then select **SVG**
to save `losatp_collinear.svg`.

## What you built

You reused one All-vs-all protein search for a different biological view:
Similarity groups answer which proteins belong together, while Collinear
blocks summarize where compatible groups retain order and orientation.

## Next steps

- [Draw Collinear blocks as a focused task](../../HOW_TO/GUI/draw-collinear-protein-blocks.md)
- [Create protein Similarity groups](../../HOW_TO/GUI/create-protein-similarity-groups.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
- [Review LOSATP result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
