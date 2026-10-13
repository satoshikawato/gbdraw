# Design E, cluster callout lists (study prototype; file `gbdraw/labels/proto/design_e.py`, never merged)

Renders and crops: `work/E/` (`<case>__E.*` default, `__E_grow.*` grow mode, `*_fan`, `*_bracket`; `view_*.png`, `crop_*.png` are flattened).
Reference look: the Owner's three plasmid maps (lists of restriction-site labels with a fan or one leader).
No pUC19/pBR322-like record is tracked in the repository. Closest: pVNTG4 37 kb (`plasmid37_f14`), pVNTG3 42.6 kb
(`plasmid42_f14__E_*`, run with run_case.py directly: same arguments as the other cases, input `work/inputs/NZ_AP024091.1.gb`;
`cases.json` untouched, `work/E/cases_e.json` is only a record) and human mtDNA (`mito_f14`). Best showcases: pVNTG3 (12 o'clock list) and mito (tRNA lists).

## 1. The rule
1. Home: a label's home is the point at its feature's angle (anchor angle) on radius r_b = r_keep + dist, text start-aligned on the right half, end-aligned on the left; box lifted by h/2*|cos| so it clears the ring at 12 and 6 o'clock.
2. Runs: labels whose home boxes overlap (gap included) are joined transitively into a run (CGView "island"). A run longer than `max_rows` (8) is cut recursively at its widest angular gap, searched only between 30 % and 70 % of its length, until every piece has at most 8 rows. A piece, or a lone label, is a block (a lone label is a list of one: today's look).
3. Side and order: a block goes to the side of its circular-mean angle (sin >= 0: right). Its rows are a flush list, left-aligned right of the ring, right-aligned left of it. Right: clockwise from top to bottom; left: counter-clockwise. Both orders are tried and the one with fewer leader crossings wins.
4. Distance: the column's inner edge is x0 = max over rows of sqrt(r_b^2 - dy^2) (dy = distance of the row box to y=0), at least 10 px from x=0, and beyond every elbow of the block's own rows + 6 px, so a block's leaders arrive from the ring side and never run through its own text. `dist` = 8 px beyond r_keep. (`hug=1`: each row follows the ring instead: staircase, slightly shorter leaders.)
5. Collisions: all blocks of one side are stacked top to bottom in angular order by `pack_ordered` (least-squares ordered packing, rigid blocks, block gap 0.4 line between lists). A block slides along the ring to make room; it never moves outward or to the other side. Rows are spaced one pitch (text height + 3 px).
6. Nothing fits: capacity rule `fit=1` (default): if a side's stack reaches beyond +-(r_b + `over` x pitch) (over = 12 rows, about 190 px beyond the ring), the label with the shortest feature in the longest of the 5 end-most blocks is dropped (listed in ctx.dropped) and the stack is re-packed, until it fits. `fit=0`: grow, no limit, every label placed, the canvas grows.
7. Leaders: `fan` = per row anchor -> elbow (radial, r_elbow) -> left (right) edge of its row. `bracket` = ONE leader from the middle feature to a vertical bar beside the list, plus a bracket over the run (arc at r_elbow+3 and a tick from each feature to the arc); other rows have no leader. `auto` = bracket when the block has >= 4 rows and at least half of the adjacent anchor gaps (measured on the elbow circle) are < 8 px, otherwise fan. Singles always have a short elbowed leader.
8. Dense genomes: lists are capped at 8 rows, so the whole-circle run becomes chunks (each a flush list) stacked in one column per side, like design D with local lists; capacity is about 2 x (ring height + 2 x 190 px) / 16 px. I did not build a second column at a larger distance: its leaders would have to pass through the first column (my first attempt at such stacking through text failed), which is design B's job.

## 2. Parameters (defaults)
leader=auto, max_rows=8, dist=8, over=12 (rows), fit=1, block_gap=0.4 (lines), bracket_min=4, sep_min=8 px, mid_gap=10 px, hug=0, link=0, flip=0 (blocks near 12/6 o'clock may move to the other side: tested, produced long crossing leaders, off), row_pad=0.

## 3. Metrics (default `fit=1` unless "grow"; cur = today). cols: labels placed / dropped / overprint / overlap / crossings / through-text / over-features / leader mean/max px / canvas
E defaults, all 0 for overprint, overlap, over-features:
| case | placed | dropped | cross | thru | leader mean/max | canvas | lists (brackets) | cur: overprint/cross/thru, leader mean/max, canvas |
| --- | --: | --: | --: | --: | --: | --- | --- | --- |
| p1_f14 | 99 | 12 | 0 | 0 | 105/321 | 1596x1315 | 17 (2) | 20/21/8, 115/467, 1706x990 |
| p1_f8 | 106 | 1 | 0 | 0 | 77/201 | 1334x1087 | 20 (4) | 0/0/0, 102/142, 1441x982 |
| t5_f8 | 135 | 59 | 0 | 0 | 139/355 | 1329x1139 | 25 (10) | 0/2085/101, 347/699, 1444x979 |
| t5_f14 | 110 | 85 | 1 | 0 | 170/462 | 1585x1265 | 18 (8) | 170/2532/168, 334/899, 1821x992 |
| n4_f14 | 66 | 3 | 0 | 0 | 90/273 | 1303x1153 | 13 (0) | 0/0/0, 124/230, 1402x994 |
| tobacco_f8 | 136 | 4 | 0 | 0 | 79/178 | 1311x1097 | 27 (7) | 0/0/0, 122/226, 1399x983 |
| tobacco_f14 | 118 | 22 | 0 | 0 | 142/304 | 1533x1322 | 21 (5) | 50/1062/110, 270/705, 1634x988 |
| mito_f14 | 33 | 0 | 0 | 0 | 96/191 | 1462x939 | 6 (0) | 0/0/0, 88/127, 1476x906 |
| plasmid37_f14 | 45 | 0 | 0 | 0 | 112/231 | 1639x976 | 7 (0) | 0/1/0, 120/248, 1648x916 |
| plasmid153_f8 | 162 | 17 | 0 | 0 | 105/262 | 1349x1153 | 29 (8) | 0/0/0, 103/207, 1436x983 |
| plasmid153_f14 | 115 | 64 | 0 | 0 | 179/294 | 1638x1326 | 29 (5) | 143/17/65, 92/177, 1714x994 |
| plasmid42 (pVNTG3) | 54 | 0 | 0 | 0 | 121/258 | 1757x1005 | 8 (0; bracket mode 8) | not rendered |

`fit=0` (grow: every label placed, 0 overprint, 0 crossings, 0 through-text in all 9 cases; canvas and longest leader grow):
p1_f14 111, 1596x1576, leader 133/491; p1_f8 107, 1334x1098; t5_f8 194, 1329x1960, 257/652; t5_f14 195, 1605x3081, 484/1207;
n4_f14 69, 1306x1191; tobacco_f14 140, 1532x1668, 165/476; tobacco_f8 140, 1311x1146; plasmid153_f8 179, 1347x1316, 122/297; plasmid153_f14 179, 1638x2049, 311/660.
Self-check inside the prototype also counts bracket strokes (ticks, arc, bar): 0 strokes through text, 0 stroke crossings between different blocks in the default runs (shared metrics ignore `extra` strokes).
Notes: the shared `label_to_feature` metric uses the feature midpoint, so it is inflated for every list. Time 0.6 to 1.3 s per case.

## 4. Assessment (from the figures)
- Where E reads better: a few tight clusters on a quiet map (pVNTG3 12 o'clock: 12 labels as two lists of 6; mito tRNA-Trp/Ala/Asn/Cys/Tyr; pVNTG4 3 o'clock). Labels of one cluster keep their clockwise order, sit in one tidy column with uniform line spacing, and the figure has no crossings, no label drifting round the ring, and a narrower canvas than today (1334 vs 1441 px for P1 at 8 pt) because labels sit at the ring. In sparse regions it is today's look (isolated labels with a short elbow leader).
- Fan vs bracket: the fan lets the reader follow one line to one gene; the bracket (crop_*_bracket.png) is much quieter (one leader per list instead of one per row; mean leader length P1 14 pt 79 vs 105 px) and matches the Owner's MCS reference, but with one shared leader the row to gene mapping is only the order along the arc: with 5 tRNAs the reader must count. I would use the bracket only when the anchors are closer than about 8 px (leaders would blur) and always keep the order clockwise. Feature ticks mark where the run is, not which row is which gene.
- vs B (tiers) and D (columns): E equals D's capacity wherever it is dense (one staircase per side) but keeps labels beside their features where they are not; it has no second tier, so it cannot add capacity by itself. vs C (free placement): E is regular and cannot lose order.
- Where E is worse: a gene-dense genome. Everything is one run, E degenerates into chunks of 8 stacked in a side column, and the capacity of one staircase per side is exceeded: at 14 pt T5 places 110 of 195, plasmid153 115 of 179 (default); or with `fit=0` all labels but the canvas grows to 3081 px (T5 14 pt) and leaders reach 1200 px. Slide is large there (mean 145 px, max 450 px), so labels are far from their features and leaders run long and nearly tangent to the ring (this was also where leader crossings appeared when a bound forced rigid shifts: the final capacity rule drops labels instead of clamping, which removed all crossings). Reading order around the circle is kept inside a list, but a list at the top or bottom (above the ring) is read off the map's circle.
- Long lists: a list of more than about 8 rows reads as a column, not as a callout; the cut at the widest gap keeps lists meaningful but the pieces then compete for one stack.
- Honest comparison with today at small font: P1, tobacco, plasmid153 at 8 pt are already clean and complete with today's placement; E there only drops labels or enlarges the canvas, so E is not a replacement. At 14 pt today overprints (P1 20 pairs, tobacco 50, plasmid153 143) and E gives zero overprint and zero crossings at the cost of dropped labels.
- Fallback for dense runs on top of A: yes, this is the best use. Rule: run today's placement; take the runs where today's result has overprint (or where a leader crosses a label), cut them out, and place only those labels as E blocks, in the ring gap next to the run, with the rest of the circle untouched. It needs (1) a run detector that reads today's result instead of the home boxes, (2) a free-space query so a list does not hit today's labels (my Scene helper), (3) the drop-and-report rule, since A's never-overprint rule has to hold, and (4) the bracket choice as an option. It would then only change maps that fail today; the sparse and clean maps stay byte-identical.

## 5. What production would need
Opt-in `--label_placement lists` next to `horizontal`/`radial`; `labels/circular_lists.py` (about 400 lines: runs, blocks, stack, two leader styles, capacity rule) plus a hook in `prepare_label_list`; the drawer must draw `leader_extra` strokes and polyline leaders (done in the study scaffold); canvas composer already grows to the label bounds. Options: `--label_list_leaders fan|bracket|auto`, `--label_list_max_rows`, `--label_overflow grow|drop` (and `drop` must report which labels it dropped). Risks: (a) the drop rule uses feature length only, there is no label importance in gbdraw; (b) lists on the left are right-aligned, which looks mirrored; (c) bracket loses row identity; (d) interactions with embedded labels, `--labels both`, per-feature label overrides and the web editor's label dragging (positions are computed, not stored).

## 6. Not done
Second column at another distance (leaders cannot pass the first column without text crossing); boxed MCS list (Owner's reference 02); balancing across sides (flip=1 exists but gave long crossing leaders, off); no inner labels (spec: only `--labels out`); not tested below 8 pt or on records with a legend overlapping the label band; the staggered elbows I tried for tangent leaders made things worse and were removed.
