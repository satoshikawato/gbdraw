# Design I ("H2"): every label beside its own feature, gentle elbows, clear of every object

2026-10-12. Study prototype, not for merging. Code: `gbdraw/labels/proto/design_i.py` in the worktree
`/home/kawato/gbdraw-work/.worktrees/label-study` (uncommitted), plus an additive scaffold edit (below).
Renders: `work/I/<case>__I.*` (29 cases), `work/I/<case>__I_in0.*` and `__I_clear6.*` (3 cases each), all completed with
`scripts/lp_objects.py`. Tables: `work/I/tab.py DESIGN[suffix]`, `work/I/cmp.py [suffix]`.

## Result

All 35 renders meet every hard criterion: 0 dropped, 0 leader crossings, 0 text overlaps, 0 leaders through text,
0 labels toward the ring, label gap 3.0 px, 0 labels on or within 3 px of an object, 0 leaders on objects, 0 bends
under 90 degrees. Targets: no leader bends under 110 degrees anywhere (smallest bend 110.0); total leader length is
below H's on all 29 cases (sum over the cases 242,373 px against 267,968); inner labels are placed by design I.
One target is missed: tRNA-Pro on `doc_tcli01` (see "Not done").

## The rule in plain steps

1. A label is a horizontal box. Its text starts at the leader and runs away from the ring (LP-D09): outside the
   ring to the right in the right column and to the left in the left column; inside the ring toward the centre line.
2. A leader has three points: the anchor on the feature, the ring edge on the same radius (1 px beyond the feature
   band or the axis stroke), the end of the text. It bends once, at the ring edge, by at most `cone` = 70 degrees,
   so every bend is 110 degrees or more and no leader leaves the ring at a grazing angle.
3. Alone, a label takes the spot with the shortest leader that keeps `clear` (3 px) from the ring and from every
   other drawn object. That is beside the feature at 3 and 9 o'clock and right above or below it at 12 and 6 o'clock.
4. Labels keep the clockwise order of their features in four stacks per shore (right/left column, upper/lower
   half). A cut at each pole decides which column the labels near that pole join; it is searched within 30 degrees of
   the pole with a cost estimate, once without and once with labels taken from across the pole, and the layout with
   the shorter total leader length is kept.
5. A stack is seated one label at a time from the equator toward its pole (inside the ring: from the pole toward
   the equator). A label starts at its own best height and takes the spot with the smallest "leader run + push along
   the stack". Each candidate spot is tested exactly against every label seated before it (box gap 3 px, leader to
   other text 2.4 px, no leader crossing) and against the other objects.
6. Reserved exits: a spot is refused when its box or leader would block the exit ray of a label that is not seated
   yet (straight along the stack, or the edge of that label's cone near the equator). So a label can always be seated
   on its own ray further along the stack: nothing is dropped and no leader has to turn sharply. On a slanted ray the
   leader meets the far corner of the text end and the text stands back 3.4 px, so that parallel leaders of close
   features fit side by side (the oblique "staircases" at 3 and 9 o'clock).
7. Inside the ring the same rule runs in the free zone between the ring and the next object. Labels asked for
   inside (`--labels both`) go there; one that does not fit goes outside. With `inside=1` and no `--labels both`, an
   outside label moves inside when the total leader length gets shorter (ten candidates by gain, then the outside is
   laid out again).

## Parameters

`clear` 3 (label box to every other object; the label-to-label gap stays `labels.spacing.circular` = 3) | `inside` 1 |
`cone` 70 | `lgap` 2.5 (leader end to its text) | `lpad` 2.4 (leader to another label's text; kept below `lgap`) |
`stub` 1 | `edge` 1 (feature outline) | `fl` 12 (a label may sit toward the equator of its feature only while a spot
within 12 px of the ring exists) | `cutmax` 30 | `cutn` 14 | `hide` 0. Deterministic (two runs give identical SVGs).
Time per render: 0.7 to 3.9 s (`p1_open_f14` 3.9 s, `mg1655_rna_f8` 3.0 s, T5 2.2 to 2.5 s; four renders in parallel).

## Obstacles: where they come from

All from `ctx.radial_layout`, nothing is read from the bare twin: axis radius and stroke, the feature band
(`features.all_band_px`), every slot's `draw_band_px` or `reserved_band_px` (GC content, GC skew, annotation and
comparison rings), the title disc (`definition_reserved_band_px`), and the tick slot's `tick_band_px` and
`label_band_px`. Tick positions and tick label texts are recomputed (`circular_tick_intervals`, `scaled_tick_text`,
text width from `calculate_bbox_dimensions`); the tick font size is derived from the label band height. A manual tick
interval would not be seen by the prototype. The ring is tested analytically, the other objects on a 1 px raster.

## H against I

| case | labels | leader_len_sum H / I | leader_len_max H / I | canvas H / I | bends < 110 H / I | near + on objects H / I | inside/outside H / I |
|---|---|---|---|---|---|---|---|
| p1_f14 | 111 | 10564 / 8975 | 413 / 260 | 1574x1224 / 1578x1230 | 92 / 0 | 8 / 0 | 0/111 / 0/111 |
| pbr322_f14 | 34 | 1507 / 1082 | 80 / 55 | 1170x874 / 1154x893 | 12 / 0 | 13 / 0 | 0/34 / 0/34 |
| mito_f14 | 33 | 1738 / 1399 | 103 / 84 | 1407x854 / 1392x885 | 9 / 0 | 13 / 0 | 0/33 / 0/33 |
| mg1655_rna_f8 | 108 | 9352 / 6975 | 300 / 190 | 1096x1025 / 1090x986 | 87 / 0 | 6 / 0 | 0/108 / 0/108 |
| tobacco_f14 | 140 | 20019 / 17411 | 419 / 403 | 1524x1436 / 1512x1431 | 134 / 0 | 0 / 0 | 0/140 / 0/140 |
| t5_f14 | 195 | 66105 / 65511 | 947 / 1030 | 1566x1985 / 1620x2630 | 184 / 0 | 1 / 0 | 0/195 / 0/195 |
| t5_f8 | 194 | 36908 / 36803 | 631 / 471 | 1318x1354 / 1328x1678 | 178 / 0 | 12 / 0 | 0/194 / 0/194 |
| plasmid153_f14 | 179 | 44260 / 42504 | 581 / 530 | 1619x1857 / 1608x1863 | 173 / 0 | 0 / 0 | 0/179 / 0/179 |
| plasmid153_f8 | 179 | 15466 / 14628 | 246 / 223 | 1328x1201 / 1314x1217 | 160 / 0 | 0 / 0 | 0/179 / 0/179 |
| lambda_f14 | 71 | 4845 / 4091 | 245 / 148 | 1445x941 / 1440x960 | 44 / 0 | 15 / 0 | 0/71 / 0/71 |
| n4_f14 | 69 | 4476 / 3624 | 276 / 195 | 1282x1036 / 1269x1046 | 43 / 0 | 10 / 0 | 0/69 / 0/69 |
| epifagus_f14 | 53 | 3138 / 2497 | 216 / 133 | 1416x914 / 1401x930 | 26 / 0 | 11 / 0 | 0/53 / 0/53 |
| plasmid37_f14 | 45 | 3180 / 2602 | 131 / 94 | 1583x876 / 1567x911 | 11 / 0 | 17 / 0 | 0/45 / 0/45 |
| p1_f8 | 107 | 5584 / 4408 | 161 / 105 | 1320x983 / 1305x966 | 48 / 0 | 14 / 0 | 0/107 / 0/107 |
| tobacco_f8 | 140 | 8795 / 7922 | 211 / 187 | 1295x1064 / 1280x1071 | 78 / 0 | 11 / 0 | 0/140 / 0/140 |
| pbr322_f10 | 34 | 1429 / 1008 | 82 / 44 | 1129x854 / 1113x869 | 11 / 0 | 12 / 0 | 0/34 / 0/34 |
| gal_hmmt_basic | 32 | 1525 / 1218 | 98 / 78 | 1479x927 / 1466x960 | 9 / 0 | 11 / 0 | 0/32 / 0/32 |
| gal_hmmt_atskew | 24 | 1032 / 698 | 105 / 54 | 1247x893 / 1236x930 | 10 / 0 | 9 / 0 | 0/24 / 0/24 |
| pbr322_open_f14 | 34 | 1507 / 1082 | 80 / 55 | 1170x874 / 1154x893 | 12 / 0 | 13 / 0 | 0/34 / 0/34 |
| lambda_open_f14 | 71 | 3864 / 3566 | 224 / 127 | 1377x898 / 1372x928 | 20 / 0 | 21 / 0 | 10/61 / 9/62 |
| p1_open_f14 | 111 | 6699 / 5107 | 264 / 121 | 1574x1081 / 1567x1031 | 50 / 0 | 16 / 0 | 28/83 / 34/77 |
| mito_open_f14 | 33 | 1549 / 1225 | 81 / 81 | 1346x854 / 1331x866 | 8 / 0 | 14 / 0 | 5/28 / 5/28 |
| plasmid37_open_f14 | 45 | 2719 / 2223 | 131 / 94 | 1515x876 / 1499x911 | 9 / 0 | 16 / 0 | 13/32 / 12/33 |
| gal_hmmtdna_basic_circular | 32 | 1525 / 1218 | 98 / 78 | 1479x927 / 1466x960 | 9 / 0 | 11 / 0 | 0/32 / 0/32 |
| gal_hmmtdna_atskew | 24 | 1032 / 698 | 105 / 54 | 1247x893 / 1236x930 | 10 / 0 | 9 / 0 | 0/24 / 0/24 |
| gal_tobacco_chloroplast_hz | 119 | 6977 / 2296 | 142 / 58 | 1336x905 / 1323x890 | 13 / 0 | 3 / 0 | 0/50 / 67/52 |
| doc_tcli01_human_mitochondrion | 24 | 1196 / 938 | 98 / 78 | 1229x916 / 1216x932 | 7 / 0 | 8 / 0 | 0/24 / 0/24 |
| doc_tcli09_precomputed_circular_rings | 24 | 872 / 615 | 76 / 38 | 1694x926 / 1681x942 | 7 / 0 | 9 / 0 | 0/24 / 0/24 |
| doc_tcli03_mitochondrial_features_highlighted | 2 | 104 / 47 | 76 / 24 | 1199x967 / 1199x967 | 1 / 0 | 1 / 0 | 0/2 / 0/2 |
| p1_open_f14 inside=0 | 111 | 10564 / 8975 | 413 / 260 | 1574x1224 / 1578x1230 | 92 / 0 | 8 / 0 | 0/111 / 0/111 |
| lambda_open_f14 inside=0 | 71 | 4845 / 4091 | 245 / 148 | 1377x941 / 1372x960 | 44 / 0 | 15 / 0 | 0/71 / 0/71 |
| pbr322_open_f14 inside=0 | 34 | 1507 / 1082 | 80 / 55 | 1170x874 / 1154x893 | 12 / 0 | 13 / 0 | 0/34 / 0/34 |
| p1_open_f14 clear=6 | 111 | 6886 / 5769 | 262 / 154 | 1574x1086 / 1567x1016 | 48 / 0 | 2 / 0 | 27/84 / 37/74 |
| lambda_open_f14 clear=6 | 71 | 4062 / 3841 | 222 / 131 | 1381x892 / 1382x948 | 18 / 0 | 1 / 0 | 10/61 / 8/63 |
| pbr322_open_f14 clear=6 | 34 | 1691 / 1215 | 95 / 56 | 1170x870 / 1160x897 | 15 / 0 | 0 / 0 | 0/34 / 0/34 |

Hard counts for I, all 35 renders (dropped / crossings / overlaps / through text / toward the ring / on objects / near
objects / leaders on objects / bends < 90): 0 / 0 / 0 / 0 / 0 / 0 / 0 / 0 / 0, `label_gap_min_px` 3.0 (285.9 on the
two-label figure). Over the 23 study records: bends under 110 degrees I 0 of 2,042, H 1,408 of 2,042.

Poles on `doc_tcli01` (whole leader, px; H / I): ATP8 77.4 / 23.9 (31 %), tRNA-Thr 65.2 / 27.2 (42 %), tRNA-Lys
42.2 / 24.2 (57 %), tRNA-Pro 72.5 / 76.8 (106 %), tRNA-Asp 32.5 / 33.8, tRNA-Phe 48.8 / 23.7.

Inner labels on `gal_tobacco_chloroplast_hz`: 67 of the 69 are inside. The piece that leaves the feature is 11.3 px
(today 47.3 px); the whole leader is 14.7 px in the median, 58.4 px at most (today 75.4 and 141.5).

What the 110 degree rule costs: with `cone=88` the total is 8,887 px on `p1_f14` (8,975 with 70), 36,581 on `t5_f8`
(36,803), 4,098 on `lambda_f14` (4,091), 964 on `doc_tcli01` (938): 1 % at most, and 21, 61, 3 and 4 leaders then
bend by less than 110 degrees. The rule does not cost length against H on any case.

## Assessment from the figures (looked at 1:1 or as noted)

- `doc_tcli01_human_mitochondrion`: ATP8, tRNA-Lys, tRNA-Ser and tRNA-Asp sit right under their features, tRNA-Thr,
  tRNA-Pro and tRNA-Phe right above theirs. Every label hugs the ring; elbows are open. The clearest gain.
- `doc_tcli03_mitochondrial_features_highlighted`: "ATP8" sits between the tick labels "9 kbp" and "8 kbp", just below
  the ring, clear of both.
- `gal_tobacco_chloroplast_hz` (whole figure and a 1:1 crop): inner labels sit just inside the ring beside their
  features (psbA, matK, atpF, atpH, rpoC1, the tRNA clusters) with leaders of a few px; outer labels as close on the
  outside. Two inner labels (tRNA-His, rps14) did not fit inside and are drawn outside.
- `p1_f14` (whole figure reduced, and a 1:1 crop of the left side): compact; canvas as H. The tRNA cluster at 9
  o'clock is an oblique staircase of parallel leaders. Towers remain at the poles: 24 labels stand wholly above the
  ring and 18 below it (H 23 and 19), up to about 14 in one tower, with vertical leaders up to 260 px (H 413). The
  label "DicB-like ..." is cut at the left canvas edge, as in H (OV-471).
- `p1_open_f14` (whole figure reduced): 34 labels inside, 77 outside, each next to its feature; longest leader
  121 px. Of all candidates so far this is the nearest to the Geneious look.
- `lambda_f14` (whole figure reduced): right side beside the features; the dense left side fans out with oblique
  leaders up to 148 px.
- `t5_f8` (1:1 crop of the right side): readable, but 136 of the 140 labels of the right half stay in the right
  column and form towers above and below the ring (40 labels wholly above, 48 below; H 43 and 44); the column stands
  40 to 60 px off the ring at 3 o'clock, where the leaders fan out.

## Not done, and open points

- tRNA-Pro on `doc_tcli01` stays at 76.8 px (target 43.5). 61.5 px of it is the piece from the anchor to the ring
  edge: the feature lies in the inner lane under the outer lane. The piece outside the ring is 15.3 px (H 12.0)
  because the label is stacked above tRNA-Thr. Only an inside label or an anchor at the ring edge would shorten it.
- T5 is taller than under H: `t5_f14` 1620x2630 (H 1566x1985), `t5_f8` 1328x1678 (H 1318x1354), and the longest
  leader on `t5_f14` is 1,030 px (H 947). T5 has 140 labels on the right half and 55 on the left. H moves 37
  (`t5_f14`) and 32 (`t5_f8`) of them across the poles into the left column, with right-angled leaders; design I
  moves 11 and 4. Under the 70 degree rule and LP-D09 a label taken across a pole needs a near-vertical leader past
  the labels already there, so the cut search finds few that pay. Total leader length is still below H's.
- Two of 69 chloroplast inner labels go outside. An inside label never crosses the centre line, and no inside
  label is placed when the free zone is narrower than 3 font sizes (so none with default GC tracks and inside ticks).
- A leader passes another label's text at 2.4 px at the closest (H: 2 px), not at `clear`.
- The cut search and the choice of inside labels are greedy; the layout is not a proven minimum.
- Tick labels are modelled, not taken from the tick renderer (see next section).

## Remarks on the metrics (nothing changed in `scripts/`)

- `labels_on_objects` / `labels_near_objects` use the PNG of the bare twin. There the tick labels are not where
  `radial_layout` says: on `mito_open_f14` the layout's `label_band_px` is 271.6 to 289.0 px, the drawn text spans
  268.4 to 282.5 px in the upper half of the map and 260.6 to 275.0 px in the lower half. With outside ticks
  (`doc_tcli03`) the layout band is 423.3 to 440.7 px and the text starts at 421.5 px or less. With the layout band
  alone, inside labels of the open cases were counted on or near tick labels. Design I therefore widens the modelled
  tick label band inward by 0.85 band heights. Whether a browser shows the SVG with the same shift was not checked;
  if it does not, this is a PNG export difference worth a bug ID.
- `leader_acute_bends` and `leader_bends_lt110` drop leader pieces shorter than 3 px, so a label seated within 3 px
  of its ring-edge point counts as straight. Design I's own bend test skips only pieces under 2.5 px.

## Files

- New: `gbdraw/labels/proto/design_i.py` (ruff clean).
- Scaffold, additive (9 lines added, 5 removed): `gbdraw/labels/circular.py` (hook: letter list `BCDEFGHI`; for design
  I the inner labels are passed on instead of being placed by today's code), `gbdraw/labels/proto/__init__.py`
  (`run_design(..., inner_labels=None)` fills `ctx.inner_items`), `gbdraw/labels/proto/common.py` (`Ctx.inner_items`).
  Byte-identical SVGs before and after for `p1_f14:D`, `p1_open_f14:H`, `lambda_f14:G`: `work/I/pre/`, `work/I/post/`,
  `work/I/scaffold_identity.txt`; checked again after the last change to `design_i.py`.
- Campaign: `work/I/` (renders, `tab.py`, `cmp.py`, `view/` crops), this report.
