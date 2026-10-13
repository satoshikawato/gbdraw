# Design J ("H3"): H2 with the turn of the leaders as a second objective

2026-10-12. Study prototype, not for merging. Code: `gbdraw/labels/proto/design_j.py` in the worktree
`/home/kawato/gbdraw-work/.worktrees/label-study` (uncommitted, 798 lines, ruff clean; the only file touched; the hook
in `gbdraw/labels/circular.py` already routes `J`). Renders: `work/J/<case>__J.*` (29 cases), `__J_in0.*` / `__J_clear6.*`
(3 cases each), series `__J_tw0 / _tw04 / _tw2 / _tw5` (5 cases), all completed with `scripts/lp_objects.py`.
Helpers I added in `work/J/` only: `cmpj.py` (J against I), `hard.py` (hard counts), `turns.py` (leaders by turn), `crop.py`,
`pair.py` (I | J crops at 1:1), `agg.py`; crops in `work/J/view/`.

## Result

All 55 renders meet every hard criterion (checked with `work/J/hard.py`): `proto.dropped` empty, 0 leader crossings, 0 text
overlaps, 0 leaders through text, 0 labels toward the ring, `label_gap_min_px` 3.0, 0 labels on or near objects, 0 leaders on
objects, 0 acute bends, 0 bends under 110 degrees. Deterministic (two runs, identical SVGs). Slowest render 8.8 s
(`p1_open_f14`, four renders in parallel; 7.2 s alone, I needed 3.9 s), the rest 0.6 to 4.5 s.
Straight leaders (turn under 3 degrees): 552 of 2,267 in I, 1,310 in J over the 29 cases. Mean turn is lower on every case,
`leaders_straight` higher on every case. Total leader length over the 29 cases: +4.4 % (253,019 against 242,373 px); the
sparse maps pay 0 to 3 %, the most any case pays is +12 % (`mg1655_rna_f8`, `p1_open_f14`). Canvas about as I's (largest
growth: `tobacco_f14` width +3 %, `plasmid153_f14` width +3 %, `epifagus_f14` height +3 %); `p1_f14` is 5 % lower, `t5_f14` 4 % lower.

## The rule: what changed against I (four changes, everything else is I)

1. Cost of a spot = leader run + `tw` x turn (degrees at the ring edge; a leader whose label sits within 3 px of the ring-edge
   point does not turn, as in the metric) + push along the stack. The turn is the only angle a leader has: anchor, ring-edge point
   on the feature's own ray, end of text. `tw` px per degree: 0.8 by default, i.e. a leader may be 0.8 px longer for each degree
   it saves.
2. New candidate for every height: the leader exactly on the feature's own ray (turn 0). The smallest run with the box clear of
   the ring (the old `hug` bound) and the end of the leader inside the box height; four more runs 2, 5, 10, 20 px longer on the same
   ray (they help when a neighbour blocks the shortest). On a slanted ray the leader may meet the text end at any height of the
   box, so it often ends at a corner of the text end (as the "staircase" leaders of I already did).
3. Plateau packing. A label alone sits where its cost is smallest; with the turn term that is a plateau of heights (any height at
   which the ray still meets the box costs the same). A label now has an interval of heights that cost at most `plat` = 2 px more
   than its best; `pack()` moves the labels of a stack toward their packed heights inside those intervals (a few
   rounds of the isotonic packing, each target projected on its interval). This makes room for the later labels of a cluster
   without moving the first one off its ray (`doc_tcli01`: 3 non-straight leaders without it, 1 with it).
4. The search along the stack (`seat`) is no longer cut off by the shortest run at the previous height (a straight spot further
   along the stack is shorter than a turned spot here): the cutoff is now the vertical gap to the feature. The cut search at
   the poles and the move of labels inside the ring compare layouts by length + `tw` x turn instead of length.
With `tw=0 ray=0 plat=0` the SVGs are byte-identical to I (checked on the five series cases): `work/J/*__J_tw0.svg`.

## Parameters

`tw` 0.8 | `ray` 1 | `plat` 2; all of I's (`clear` 3, `inside` 1, `cone` 70, `lgap` 2.5, `lpad` 2.4, `stub` 1, `edge` 1, `fl` 12,
`cutmax` 30, `cutn` 14, `hide` 0). Why 0.8: over the 29 cases (with `plat` 2) tw 0.6 gives 1,250 straight leaders at +3.8 % length (worst
case +11 %), 0.8 gives 1,310 at +4.4 % (worst +12 %), 1.0 gives 1,351 at +4.9 % (worst +17 %, `p1_open_f14`), 2.0 gives 1,437 at
+6.4 % (worst +25 %); beyond 0.8 the extra straight leaders cost more length on the dense maps than the sparse ones gain.
`plat` 5 and 10 give the same totals (1,310 / 1,309), so it is not a sensitive knob.

## J against I, per case (labels = leaders drawn; turn in degrees)

| case | labels | straight I / J | mean turn I / J | max turn I / J | leader_len_sum I / J (ratio) | leader_len_max I / J | canvas I / J |
|---|---|---|---|---|---|---|---|
| p1_f14 | 111 | 18 / 62 | 35.1 / 14.2 | 69.9 / 68.5 | 8975 / 8533 (0.95) | 260 / 210 | 1578x1230 / 1575x1171 |
| pbr322_f14 | 34 | 17 / 33 | 19.0 / 0.5 | 69.9 / 15.7 | 1082 / 1092 (1.01) | 55 / 56 | 1154x893 / 1154x898 |
| mito_f14 | 33 | 19 / 31 | 17.1 / 0.6 | 69.6 / 12.0 | 1399 / 1419 (1.01) | 84 / 96 | 1392x885 / 1392x889 |
| mg1655_rna_f8 | 108 | 14 / 32 | 37.2 / 21.1 | 70.0 / 68.5 | 6975 / 7829 (1.12) | 190 / 211 | 1090x986 / 1094x1003 |
| tobacco_f14 | 140 | 16 / 50 | 37.1 / 19.0 | 70.0 / 68.5 | 17411 / 17904 (1.03) | 403 / 412 | 1512x1431 / 1560x1434 |
| t5_f14 | 195 | 7 / 33 | 47.0 / 36.5 | 70.0 / 69.6 | 65511 / 69609 (1.06) | 1030 / 920 | 1620x2630 / 1602x2532 |
| t5_f8 | 194 | 16 / 55 | 44.5 / 32.8 | 70.0 / 68.5 | 36803 / 37210 (1.01) | 471 / 473 | 1328x1678 / 1331x1680 |
| plasmid153_f14 | 179 | 8 / 19 | 44.5 / 35.8 | 70.0 / 68.5 | 42504 / 43589 (1.03) | 530 / 543 | 1608x1863 / 1659x1877 |
| plasmid153_f8 | 179 | 13 / 69 | 41.4 / 18.1 | 70.0 / 68.5 | 14628 / 16217 (1.11) | 223 / 232 | 1314x1217 / 1345x1234 |
| lambda_f14 | 71 | 15 / 61 | 28.9 / 2.9 | 69.7 / 35.5 | 4091 / 4372 (1.07) | 148 / 158 | 1440x960 / 1436x982 |
| n4_f14 | 69 | 20 / 57 | 27.2 / 3.6 | 70.0 / 66.5 | 3624 / 3813 (1.05) | 195 / 203 | 1269x1046 / 1266x1050 |
| epifagus_f14 | 53 | 20 / 46 | 25.6 / 2.6 | 69.8 / 41.1 | 2497 / 2702 (1.08) | 133 / 151 | 1401x930 / 1401x961 |
| plasmid37_f14 | 45 | 28 / 45 | 14.8 / 0.0 | 64.6 / 0.3 | 2602 / 2599 (1.00) | 94 / 96 | 1567x911 / 1567x915 |
| p1_f8 | 107 | 37 / 88 | 28.2 / 6.7 | 69.9 / 68.5 | 4408 / 4419 (1.00) | 105 / 108 | 1305x966 / 1307x964 |
| tobacco_f8 | 140 | 33 / 107 | 29.6 / 6.8 | 69.9 / 68.5 | 7922 / 8579 (1.08) | 187 / 202 | 1280x1071 / 1279x1094 |
| pbr322_f10 | 34 | 22 / 32 | 12.6 / 0.9 | 68.2 / 24.7 | 1008 / 1014 (1.01) | 44 / 46 | 1113x869 / 1113x871 |
| gal_hmmt_basic | 32 | 18 / 31 | 17.6 / 0.4 | 69.6 / 13.7 | 1218 / 1239 (1.02) | 78 / 88 | 1466x960 / 1466x964 |
| gal_hmmt_atskew | 24 | 11 / 19 | 24.2 / 7.6 | 69.6 / 60.5 | 698 / 710 (1.02) | 54 / 68 | 1236x930 / 1235x933 |
| pbr322_open_f14 | 34 | 17 / 33 | 19.0 / 0.5 | 69.9 / 15.7 | 1082 / 1092 (1.01) | 55 / 56 | 1154x893 / 1154x898 |
| lambda_open_f14 | 71 | 24 / 59 | 24.4 / 2.6 | 69.7 / 35.2 | 3566 / 3838 (1.08) | 127 / 158 | 1372x928 / 1368x944 |
| p1_open_f14 | 111 | 26 / 74 | 31.9 / 10.0 | 70.0 / 68.5 | 5107 / 5719 (1.12) | 121 / 133 | 1567x1031 / 1569x1034 |
| mito_open_f14 | 33 | 19 / 31 | 17.3 / 0.8 | 69.2 / 19.4 | 1225 / 1289 (1.05) | 81 / 84 | 1331x866 / 1331x889 |
| plasmid37_open_f14 | 45 | 26 / 44 | 15.2 / 0.1 | 64.6 / 5.3 | 2223 / 2250 (1.01) | 94 / 96 | 1499x911 / 1499x915 |
| gal_hmmtdna_basic_circular | 32 | 18 / 31 | 17.6 / 0.4 | 69.6 / 13.7 | 1218 / 1239 (1.02) | 78 / 88 | 1466x960 / 1466x964 |
| gal_hmmtdna_atskew | 24 | 11 / 19 | 24.2 / 7.6 | 69.6 / 60.5 | 698 / 710 (1.02) | 54 / 68 | 1236x930 / 1235x933 |
| gal_tobacco_chloroplast_hz | 119 | 54 / 101 | 22.8 / 4.2 | 70.0 / 68.5 | 2296 / 2400 (1.05) | 58 / 61 | 1323x890 / 1325x914 |
| doc_tcli01_human_mitochondrion | 24 | 12 / 23 | 20.4 / 0.6 | 69.6 / 13.7 | 938 / 956 (1.02) | 78 / 88 | 1216x932 / 1215x934 |
| doc_tcli09_precomputed_circular_rings | 24 | 12 / 23 | 20.4 / 0.6 | 69.6 / 13.7 | 615 / 631 (1.03) | 38 / 48 | 1681x942 / 1680x944 |
| doc_tcli03_mitochondrial_features_highlighted | 2 | 1 / 2 | 2.0 / 0.0 | 4.0 / 0.0 | 47 / 47 (1.00) | 24 / 24 | 1199x967 / 1199x967 |
| pbr322_open_f14 inside=0 | 34 | 17 / 33 | 19.0 / 0.5 | 69.9 / 15.7 | 1082 / 1092 (1.01) | 55 / 56 | 1154x893 / 1154x898 |
| lambda_open_f14 inside=0 | 71 | 15 / 61 | 28.9 / 2.9 | 69.7 / 35.5 | 4091 / 4372 (1.07) | 148 / 158 | 1372x960 / 1368x982 |
| p1_open_f14 inside=0 | 111 | 18 / 62 | 35.1 / 14.2 | 69.9 / 68.5 | 8975 / 8533 (0.95) | 260 / 210 | 1578x1230 / 1575x1171 |
| pbr322_open_f14 clear=6 | 34 | 6 / 33 | 20.4 / 0.4 | 69.1 / 14.5 | 1215 / 1221 (1.00) | 56 / 59 | 1160x897 / 1160x905 |
| lambda_open_f14 clear=6 | 71 | 12 / 61 | 24.9 / 2.3 | 69.7 / 63.0 | 3841 / 4074 (1.06) | 131 / 140 | 1382x948 / 1375x958 |
| p1_open_f14 clear=6 | 111 | 15 / 68 | 31.3 / 10.5 | 69.9 / 68.5 | 5769 / 6159 (1.07) | 154 / 136 | 1567x1016 / 1577x1033 |

## Trade-off series (suffix; `tw=0` is `tw=0,ray=0,plat=0`, i.e. H2 exactly; default = no suffix)

| case | tw 0 (H2) | tw 0.4 (`_tw04`) | tw 0.8 default | tw 2 (`_tw2`) | tw 5 (`_tw5`) |
|---|---|---|---|---|---|
| doc_tcli01 straight / mean turn / length / canvas | 12/24, 20.4, 938, 1216x932 | 23/24, 0.6, 956, 1215x934 | 23/24, 0.6, 956, 1215x934 | 23/24, 0.4, 959, 1215x934 | 23/24, 0.4, 974, 1215x949 |
| p1_f8 | 37/107, 28.2, 4408, 1305x966 | 82/107, 7.5, 4373, 1307x962 | 88/107, 6.7, 4419, 1307x964 | 92/107, 5.5, 4665, 1308x988 | 94/107, 4.0, 4968, 1310x1013 |
| p1_f14 | 18/111, 35.1, 8975, 1578x1230 | 58/111, 17.4, 9231, 1575x1216 | 62/111, 14.2, 8533, 1575x1171 | 70/111, 11.9, 8860, 1575x1173 | 73/111, 10.8, 9752, 1584x1227 |
| lambda_f14 | 15/71, 28.9, 4091, 1440x960 | 51/71, 4.9, 4217, 1436x973 | 61/71, 2.9, 4372, 1436x982 | 61/71, 2.9, 4381, 1436x984 | 61/71, 3.0, 4470, 1436x990 |
| t5_f8 | 16/194, 44.5, 36803, 1328x1678 | 55/194, 34.6, 36926, 1328x1678 | 55/194, 32.8, 37210, 1331x1680 | 67/194, 31.1, 37614, 1331x1680 | 68/194, 31.0, 40277, 1337x1721 |

Reading: tw 0.4 already gets most of the straight leaders (the mean turn falls from 20 to 1 degree on the mitochondrial
map); above about 1 the length grows faster than the straightness (p1_f8: 4 to 5 more leaders for 250 px at tw 2, 300 px more
at tw 5). tw 5 starts to cost canvas against the default (`p1_f14` +5 % height, `p1_f8` +5 %).

## Where the leaders are not straight, and why

- Straight leaders are exactly radial (one line from the feature centre through the ring edge to the label); none leans, so
  the lean angle is 0 on every straight leader and no leader grazes the feature band. I did not build slanted straight
  lines (see "Not done").
- `doc_tcli01` 23 of 24 (target met): the one that is not is tRNA-Tyr (13.7 degrees), the last of the five-label fan at 4
  o'clock (Trp, Ala, Asn, Cys, Tyr: features within 7 degrees, about 50 px of arc, five labels need 82 px): the fifth label cannot reach its own ray
  without moving the four labels above it more than their plateau allows. Same for `gal_hmmtdna_basic_circular` (31 of 32,
  Tyr 13.7) and `doc_tcli09` (23 of 24, Tyr 13.7); `mito_f14` 31 of 33 (Tyr 12.0, Asn 8.4, the same fan); `pbr322_f14` 33 of
  34 (BspEI 1664, 15.7 degrees: the text of BsgI 1650, 1 degree away, lies across its ray). `doc_tcli03` 2 of 2 (ATP8 is straight, as before).
  In H2 the mean turn on these maps was 17 to 20 degrees, now 0.4 to 0.6.
- Dense maps keep turns where features lie closer than a text line (tRNA clusters of `p1_f8`/`p1_f14`, the 9 o'clock group of
  lambda, the towers at 12 and 6 o'clock of `t5_f8`): with features 5 to 10 px apart and 16.5 px per label, only about every third
  label can stand on its own ray, the others fan out to the 70 degree edge of the cone (max turn 68.5). 104 leaders of `t5_f8`
  and 19 of `p1_f14` still turn 30 degrees or more (H2: 139 and 60).

## Assessment from the figures (1:1 crops side by side with I in `work/J/view/`)

- `doc_tcli01` 12 o'clock (`tcli01_top_pair.png`): tRNA-Pro, Thr, Phe, Glu, Val are all one straight line now; Pro and Thr rise
  as two nearly parallel lines past the end of the CYTB arrow. Whole map: `tcli01_J.png`; the labels stay close to the ring, a few
  px further out than in H2 (the longest leader is 88 px against 78).
- `doc_tcli03` ATP8 (`tcli03_pair.png`): unchanged, straight.
- `gal_tobacco_chloroplast_hz` upper left (`chl_pair.png`): same look as H2; the leaders of rps15, ndhH, ndhA, ndhI, ndhG,
  ndhE, psaC are straight, a little longer. Whole figure: 101 of 119 leaders straight (H2: 54).
- `p1_f8` (whole at 1:1 and the 70 kbp cluster): the rays fan out like spokes and read quietly; the tRNA cluster at 9 o'clock still
  has slanted leaders to the 70 degree edge of the cone but fewer.
- `lambda_f14` left side (`lambda_pair.png`): the leaders are parallel radial lines and the labels staircase down them; the
  worst looks are the leaders that end at the top or bottom corner of a text end (a long label such as "phosphoadenosine
  phosphosulfate reductase" hangs beside the end of a line that meets its corner), which H2 had too.
- `t5_f8` right side (`t5_pair.png`): nearly the same as H2 (towers, a fan of oblique lines); J is 55 of 194 straight against
  16, mean turn 33 against 44. Not made worse: 1331x1680 against 1328x1678.
- The towers on `p1_f14` are not worse: 1575x1171 against 1578x1230, longest leader 210 against 260 (not my task; not examined further).

## Not done / open

- Slanted straight lines (idea 3 of the brief) and a sliding anchor (idea 2): not built. The label dict gives only
  `start`/`end` of the label, not of the feature (`label_start = middle -/+ 0.55 x text width in bp`), so the prototype cannot
  tell how far a leader may lean inside its feature. Scaffold change needed, in `gbdraw/labels/circular.py` where the label
  entry is filled (about line 4766): add `feature_start`/`feature_end` = `candidate.segment_start_bp` / `segment_end_bp`
  (already known there). With them, a leader from any point of the feature's middle part could be straight along a lean of
  up to about 25 degrees, which would help exactly the staircases and the five-label fan left above.
- No "last pass" over the finished layout (idea 4): plateau packing does the same job for clusters before seating, and I did not
  find a case where a pass after seating would add straight leaders without moving a neighbour.
- Dense clusters and the towers: greedy per label. The fan leaders are not a proven minimum of the sum of turns (a joint search per
  stack would be the next step, not tried).
- `leaders_straight` counts leaders turning under 3 degrees; a label seated within 3 px of its ring-edge point counts as straight
  although its leader is a stub (the metric's rule). Not changed.
- The remaining +12 % cases (`mg1655_rna_f8`, `p1_open_f14`) are dense; `p1_open_f14` is the slowest (7.2 s alone, 8.8 s in a parallel batch; I needed 3.9 s): the extra candidates and the longer search along a stack.
