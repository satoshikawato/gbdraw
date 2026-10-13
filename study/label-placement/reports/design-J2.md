# Design J2: symmetric fans (LP-D12) and leaders that end anywhere on the label box (LP-D13)

2026-10-12. Study prototype, not for merging. Code: `gbdraw/labels/proto/design_j.py` in the worktree
`/home/kawato/gbdraw-work/.worktrees/label-study` (uncommitted change on top of 819443b7, ruff clean; no other file touched).
Renders: `work/J/<case>__J_sym.*` (29 cases, defaults `sym=1 edge=1`), `__J_sym_in0.*` / `__J_sym_clear6.*` (3 cases each),
`__J_sym_e0.*` (5 cases, `sym=1 edge=0`), `__J_e1s0.*` (29 cases, `sym=0 edge=1`, to separate the two changes); all completed with
`scripts/lp_objects.py`. The `__J.*` renders of H3 are untouched. Scripts in `work/J/`: `lop.py` (clusters), `back.py`, `cmp2.py`
(tables), `pair2.py` (crops); pairs J | J2 in `work/J/view/*_sym_pair.png` (tcli01, mito, pbr322, lambda, p1f8, chl).

## Result

- `sym=0 edge=0` is byte-identical to the committed J (SVG compared on the five series cases; `edge` of design I, the feature outline
  allowance, is now option `fedge`, default unchanged).
- Hard counts on all 40 J2 renders: dropped 0, crossings 0, text overlaps 0, leaders through text 0, label gap >= 3, labels on or near
  objects 0, leaders on objects 0, acute bends 0, bends under 110 degrees 0. `labels_toward_axis` is reported, not a criterion (LP-D13):
  0 on all 29 cases with `edge=0`, 53 labels over the 29 cases with `edge=1` (largest: t5_f8 16, t5_f14 9, n4_f14 7).
- Over the 29 cases, J -> J2: lopsided clusters 98 -> 61; straight leaders 1,310 -> 1,187 (-123, of 2,267); total leader length
  253,019 -> 239,851 px (-5.2 %). Split by cause: `edge=1` alone (`sym=0`): straight 1,310 -> 1,243, length -5.8 %, lopsided 98 -> 100;
  `sym=1` on top of it: lopsided 100 -> 61, straight 1,243 -> 1,187, length +0.6 %. Time: slowest 9.7 s (`p1_open_f14`, four renders in
  parallel; J needed 8.8 s); the rest under 7 s.
- The symmetric fan works where the cluster is a fan of 2 to 8 labels (the Owner's examples): lopsided 0 on `doc_tcli01`, `mito_f14`,
  `pbr322_f14`, `gal_hmmt*`, `doc_tcli09`, `pbr322_open`/`lambda_open`/`mito_open`; lambda_f14 3 -> 2, chloroplast 7 -> 1. It does not work on the
  towers of dense maps (`p1_f14` 7 -> 7, `tobacco_f14` 5 -> 5, `plasmid153` 4 -> 4, `t5_f8` 3 -> 6): see "Not done".

## What changed against J

1. `edge=1` (LP-D13). A leader may end at any point of the perimeter of the label box grown by `lgap` (2.5 px between leader end and
   text). New candidates per height: the feature's own ray and the two edges of the 70 degree cone, ending on the top or the bottom edge
   of the box (end at 20, 50 or 80 % of the edge; the box may stand centred on its ray, so a label right above or below its feature
   has a straight leader); and the point of the grown box nearest to the stub (box right of the stub, or under or over it). The
   rule "no point of the box nearer the ring than the leader end by more than `back`" replaces LP-D09's far-end test; `back` default
   0.5 x text height (8.6 px at 14 pt, 4 px at 8 pt), measured as: distance of the leader end from the ring minus the distance of the
   box's nearest point (outside) or farthest point (inside); the stub inside the grown box is refused.
2. `sym=1`. After the layout is final, per shore: clusters = chains of labels of one stack (column and half) whose boxes are at most
   box gap + 3 px apart, with a leader that turns 3 degrees or more. A cluster is lopsided when |sum of signed turns| > half the sum of
   |turns| and the latter is at least 6 degrees (sign: turn clockwise in the image; `work/J/lop.py` uses the same definition on the
   finished figure, gap 6.1 px). A lopsided cluster is shifted rigidly along its stack (0, +-4 px and more, up to 40 px; up to 200 for
   clusters of more than 8) and every label is seated again at its shifted height with a wanted lean: a ramp from -T at one end to +T at
   the other, zero in the middle 40 % (T = 0.5 or 1.0 x the largest turn of the cluster in J, the end that carried the turn leans its way,
   the other end the opposite way); the cost of a spot gets 0.5 x (turn - wanted)^2 and the leader candidates include leaders at exactly
   the wanted lean. Of the 6 to 30 best shifts by estimated cost the one with the best balance (|sum| / sum |turns|) is taken, provided
   |sum| fell by at least 3 degrees, no turn exceeds the cluster's old largest turn + 0.5, and the cluster's cost rose by at most 6 px
   per label; otherwise the cluster stays as J had it. Labels outside lopsided clusters are not touched.
3. Parameters: `sym` 1, `edge` 1, `back` -1 (= 0.5 x text height); all of J's. `tw` stays 0.8.

## J (sym=0 edge=0) against J2 (sym=1 edge=1), per case (lopsided = clusters as defined above; back = labels with more than 1 px, largest value)

| case | labels | lopsided A/B | straight A/B | mean turn A/B | max turn A/B | len sum A/B (ratio) | canvas A/B | toward axis A/B | back>1px A/B (max B) |
|---|---|---|---|---|---|---|---|---|---|
| p1_f14 | 111 | 7 / 7 | 62 / 57 | 14.2 / 13.8 | 68.5 / 69.4 | 8533 / 8221 (0.96) | 1575x1171 / 1575x1155 | 0 / 4 | 13 / 17 (8.6) |
| pbr322_f14 | 34 | 1 / 0 | 33 / 31 | 0.5 / 0.4 | 15.7 / 5.3 | 1092 / 1029 (0.94) | 1154x898 / 1153x893 | 0 / 0 | 0 / 2 (7.3) |
| mito_f14 | 33 | 1 / 0 | 31 / 31 | 0.6 / 1.0 | 12.0 / 17.5 | 1419 / 1371 (0.97) | 1392x889 / 1391x888 | 0 / 0 | 1 / 2 (5.9) |
| mg1655_rna_f8 | 108 | 8 / 4 | 32 / 32 | 21.1 / 17.7 | 68.5 / 57.0 | 7829 / 7920 (1.01) | 1094x1003 / 1132x1034 | 0 / 3 | 2 / 7 (4.9) |
| tobacco_f14 | 140 | 5 / 5 | 50 / 46 | 19.0 / 19.3 | 68.5 / 68.5 | 17904 / 17882 (1.00) | 1560x1434 / 1577x1436 | 0 / 0 | 4 / 10 (6.9) |
| t5_f14 | 195 | 3 / 5 | 33 / 24 | 36.5 / 37.1 | 69.6 / 69.6 | 69609 / 64212 (0.92) | 1602x2532 / 1603x2351 | 0 / 9 | 28 / 30 (8.5) |
| t5_f8 | 194 | 3 / 6 | 55 / 53 | 32.8 / 34.0 | 68.5 / 69.2 | 37210 / 34011 (0.91) | 1331x1680 / 1328x1517 | 0 / 16 | 8 / 21 (5.0) |
| plasmid153_f14 | 179 | 4 / 4 | 19 / 18 | 35.8 / 36.1 | 68.5 / 69.6 | 43589 / 42886 (0.98) | 1659x1877 / 1648x1868 | 0 / 0 | 2 / 1 (2.3) |
| plasmid153_f8 | 179 | 4 / 4 | 69 / 60 | 18.1 / 19.8 | 68.5 / 69.0 | 16217 / 15974 (0.99) | 1345x1234 / 1340x1234 | 0 / 0 | 2 / 5 (4.9) |
| lambda_f14 | 71 | 3 / 2 | 61 / 50 | 2.9 / 4.3 | 35.5 / 38.2 | 4372 / 4274 (0.98) | 1436x982 / 1436x962 | 0 / 2 | 8 / 15 (8.6) |
| n4_f14 | 69 | 5 / 2 | 57 / 50 | 3.6 / 3.2 | 66.5 / 28.0 | 3813 / 3285 (0.86) | 1266x1050 / 1266x1050 | 0 / 7 | 8 / 14 (7.9) |
| epifagus_f14 | 53 | 5 / 2 | 46 / 33 | 2.6 / 5.2 | 41.1 / 24.5 | 2702 / 2508 (0.93) | 1401x961 / 1400x936 | 0 / 2 | 4 / 7 (6.2) |
| plasmid37_f14 | 45 | 0 / 0 | 45 / 43 | 0.0 / 0.2 | 0.3 / 3.3 | 2599 / 2536 (0.98) | 1567x915 / 1566x913 | 0 / 0 | 0 / 0 (0.2) |
| p1_f8 | 107 | 6 / 5 | 88 / 83 | 6.7 / 6.4 | 68.5 / 68.5 | 4419 / 4058 (0.92) | 1307x964 / 1305x949 | 0 / 2 | 3 / 10 (4.1) |
| tobacco_f8 | 140 | 6 / 6 | 107 / 78 | 6.8 / 8.9 | 68.5 / 68.5 | 8579 / 8133 (0.95) | 1279x1094 / 1279x1073 | 0 / 0 | 4 / 9 (3.8) |
| pbr322_f10 | 34 | 2 / 0 | 32 / 33 | 0.9 / 0.1 | 24.7 / 3.3 | 1014 / 944 (0.93) | 1113x871 / 1113x868 | 0 / 0 | 1 / 1 (5.2) |
| gal_hmmt_basic | 32 | 1 / 0 | 31 / 30 | 0.4 / 1.1 | 13.7 / 15.9 | 1239 / 1201 (0.97) | 1466x964 / 1465x963 | 0 / 0 | 1 / 3 (5.8) |
| gal_hmmt_atskew | 24 | 3 / 0 | 19 / 18 | 7.6 / 6.4 | 60.5 / 32.7 | 710 / 671 (0.95) | 1235x933 / 1241x913 | 0 / 0 | 1 / 3 (6.3) |
| pbr322_open_f14 | 34 | 1 / 0 | 33 / 31 | 0.5 / 0.4 | 15.7 / 5.3 | 1092 / 1029 (0.94) | 1154x898 / 1153x893 | 0 / 0 | 0 / 2 (7.3) |
| lambda_open_f14 | 71 | 5 / 0 | 59 / 53 | 2.6 / 3.0 | 35.2 / 22.1 | 3838 / 3529 (0.92) | 1368x944 / 1368x931 | 0 / 0 | 9 / 14 (7.6) |
| p1_open_f14 | 111 | 11 / 8 | 74 / 67 | 10.0 / 10.7 | 68.5 / 68.5 | 5719 / 5195 (0.91) | 1569x1034 / 1577x1030 | 0 / 4 | 5 / 18 (8.2) |
| mito_open_f14 | 33 | 1 / 0 | 31 / 30 | 0.8 / 0.7 | 19.4 / 7.3 | 1289 / 1209 (0.94) | 1331x889 / 1330x883 | 0 / 0 | 2 / 4 (4.6) |
| plasmid37_open_f14 | 45 | 0 / 0 | 44 / 41 | 0.1 / 0.6 | 5.3 / 6.5 | 2250 / 2180 (0.97) | 1499x915 / 1498x913 | 0 / 0 | 1 / 1 (6.2) |
| gal_hmmtdna_basic_circular | 32 | 1 / 0 | 31 / 30 | 0.4 / 1.1 | 13.7 / 15.9 | 1239 / 1201 (0.97) | 1466x964 / 1465x963 | 0 / 0 | 1 / 3 (5.8) |
| gal_hmmtdna_atskew | 24 | 3 / 0 | 19 / 18 | 7.6 / 6.4 | 60.5 / 32.7 | 710 / 671 (0.95) | 1235x933 / 1241x913 | 0 / 0 | 1 / 3 (6.3) |
| gal_tobacco_chloroplast_hz | 119 | 7 / 1 | 101 / 101 | 4.2 / 2.4 | 68.5 / 31.8 | 2400 / 2148 (0.89) | 1325x914 / 1328x899 | 0 / 4 | 1 / 5 (2.6) |
| doc_tcli01_human_mitochondrion | 24 | 1 / 0 | 23 / 22 | 0.6 / 1.4 | 13.7 / 15.9 | 956 / 930 (0.97) | 1215x934 / 1214x934 | 0 / 0 | 1 / 2 (5.8) |
| doc_tcli09_precomputed_circular_rings | 24 | 1 / 0 | 23 / 22 | 0.6 / 1.2 | 13.7 / 13.1 | 631 / 602 (0.95) | 1680x944 / 1679x944 | 0 / 0 | 1 / 2 (5.9) |
| doc_tcli03_mitochondrial_features_highlighted | 2 | 0 / 0 | 2 / 2 | 0.0 / 0.0 | 0.0 / 0.0 | 47 / 43 (0.91) | 1199x967 / 1199x967 | 0 / 0 | 0 / 0 (0.0) |
| SUM | 2267 | 98 / 61 | 1310 / 1187 |  |  | 253019 / 239851 (0.948) |  | 0 / 53 | 112 / 211 |

## What the box-edge rule gains (`sym=1`, edge 0 -> 1)

| case | edge | lopsided | straight | mean turn | max turn | len sum | canvas | toward axis | back>1px (max) |
|---|---|---|---|---|---|---|---|---|---|
| doc_tcli01_human_mitochondrion | 0 | 0 | 22/24 | 1.3 | 13.7 | 966 | 1215x934 | 0 | 1 (5.9) |
| doc_tcli01_human_mitochondrion | 1 | 0 | 22/24 | 1.4 | 15.9 | 930 | 1214x934 | 0 | 2 (5.8) |
| lambda_f14 | 0 | 1 | 53/71 | 4.1 | 26.1 | 4406 | 1436x982 | 0 | 9 (7.2) |
| lambda_f14 | 1 | 2 | 50/71 | 4.3 | 38.2 | 4274 | 1436x962 | 2 | 15 (8.6) |
| p1_f8 | 0 | 4 | 83/107 | 7.0 | 68.5 | 4411 | 1308x964 | 0 | 3 (1.8) |
| p1_f8 | 1 | 5 | 83/107 | 6.4 | 68.5 | 4058 | 1305x949 | 2 | 10 (4.1) |
| p1_f14 | 0 | 6 | 61/111 | 14.4 | 68.5 | 8536 | 1575x1171 | 0 | 13 (10.0) |
| p1_f14 | 1 | 7 | 57/111 | 13.8 | 69.4 | 8221 | 1575x1155 | 4 | 17 (8.6) |
| gal_tobacco_chloroplast_hz | 0 | 2 | 97/119 | 4.7 | 68.5 | 2400 | 1325x914 | 0 | 3 (2.0) |
| gal_tobacco_chloroplast_hz | 1 | 1 | 101/119 | 2.4 | 31.8 | 2148 | 1328x899 | 4 | 5 (2.6) |

Reading: `edge=1` shortens the leaders (-3 to -10 % on these five; -5.8 % over the 29 cases) and the canvas by up to 2 % in height on these five (T5 8 pt 10 %), because labels
sit right beside or above their feature; it does not add straight leaders (more labels are straight at the same cost, but with `tw` 0.8 a
shorter slightly leaning leader beats a longer straight one: 1,310 -> 1,243 over the 29 cases). Raising `tw` with `edge=1` would trade the
length back; I did not tune it. The chloroplast figure gains on both counts (97 -> 101 straight, length -10 %). In the 1:1 crop
(`chl_sym_pair.png`) the labels stand closer to their features. `back`: with 8.6 px, 112 -> 211 labels over the 29 cases use more than 1 px
(largest 8.6; with the old test the same measure gives up to 11.6 on `lambda_f14`, so J already had such labels); by eye the text still
reads as running away from the ring. 53 labels are flagged by `labels_toward_axis` (the metric's own stricter rule).

## Assessment from the figures (1:1, J left, J2 right)

- `doc_tcli01` (`tcli01_sym_pair.png`): Trp ... Tyr was one leaning leader (Tyr 13.7) and four straight; now Trp leans 15.9 one way and Tyr
  15.9 the other, the three between are straight: a symmetric fan, two leaders lean instead of one (22 instead of 23 straight).
- `pbr322_f14` (BsgI / BspEI, `pbr322_sym_pair.png`): the two labels diverge symmetrically (15.7 degrees on one -> 5.3 on the pair).
- `mito_f14`, `gal_hmmtdna_basic`, `doc_tcli09`: the same five-tRNA fan; the largest turn is 17.5 / 15.9 / 13.1 against 12.0 / 13.7 / 13.7 in J
  (the symmetric fan has two leaning leaders; on these the largest turn is no longer below J's: the cap is against the cluster as
  `edge=1` laid it out, not against J).
- `lambda_f14` left side and `p1_f8` 70 kbp cluster: more leaders end at the box side, the fans look like J's; the remaining lopsided
  clusters are tails of 4 to 30 labels leaning one way.

## Not done, and open

- Sliding anchor and a leaning straight line (second priority): not started; no scaffold edit made (`gbdraw/labels/circular.py` is untouched).
- Towers of dense maps stay one-sided: 61 clusters are still lopsided, mostly the towers at the poles and the 9 o'clock clusters of
  `p1_f14`, `tobacco`, `plasmid153`, `t5` (turn 68.5 at the cone edge on one end). Their symmetric form needs the whole stack to stand
  centred on its features, which a rigid shift of the cluster with a 6 px/label cost cap does not find; a joint seating of the stack
  would. On `t5_f8` lopsided went 3 -> 6 because the box-edge rule changes the layout before the symmetric pass (the pass found nothing).
- Cost of symmetry: -56 straight leaders over the 29 cases (1,243 -> 1,187), mean turn rises slightly on small maps (tcli01 0.6 -> 1.4).
- `labels_toward_axis` 0 -> 53 and the `back` tolerance (default 0.5 x text height) are for the Owner to judge by eye.
- The cluster definition (touching boxes) is mine; the metric and the pass use gaps of 6.1 and 6.05 px.

## LP-D14: a corner is not a target (`slide`, default 1)

Owner on `tcli01_sym_pair.png`: tRNA-Ala need not reach the top-left corner; a straight leader that meets the left edge keeps its line
and saves length. Code: same file, `slide` option plus the `_Eng.slide` post-pass (diff to the first J2 hand-back:
`logs/j2-slide.diff`, 92 changed lines, ruff clean). Renders `work/J/<case>__J_sym2.*` (29 cases), `__J_sym2_in0.*` / `__J_sym2_clear6.*`
(pbr322, lambda, p1 `_open_f14`), all completed with `lp_objects.py` (bare twins in `page/renders` and `work/docs`). `__J_sym*` untouched.

**Diagnosis (tRNA-Ala, doc_tcli01).** No clearance blocked it. Ala's box is pinned (Trp above, gap 3.05 px; ring clearance at its
top-left corner), so the box cannot get closer, and it did not have to: the line S->anchor meets the box edge at (367.8, 230.95), 1.5-2.8 px
below the corner. The old end was a sideways stub of 2.8 px to (367.8, 229.15), i.e. at the corner row. Cause: the cost of a spot
waives the turn at S when the run is under 3 px (`nostraight`, the metric drops pieces under 3 px), so the 2.8 px sideways stub (32 deg off
the ray) cost the same as nothing and beat the on-ray end (run 3.3 px). Not the sym pass, not the stack order. `edge=1` was working;
the short-run waiver hid the sideways step from the cost.

**Rule now.** After the layout (and the sym pass), `slide` re-seats each label alone, with `nostraight` 0.5 only inside this pass so a
short sideways stub counts as the turn it is. Candidates: 0..26 px along the stack direction (11 steps) through the normal `try_at` tests (3 px to
objects, label gap 3 px, `back` limit, cone, no crossing, against all other labels and leaders). Accepted if cost (run + 0.8 x turn) is >= 0.5
px lower, turn within 1 deg of the old (a straight leader stays straight, a fan keeps its lean), and the whole leader at most 1.5 px longer.
Setting `nostraight` to 0.5 in the main layout was tried and rejected (Phe/Pro at 12 o'clock leaned 12-25 deg, straight 22 -> 19 on tcli01,
length +30).

**Proof `slide=0`.** SVGs byte-identical to `__J_sym.svg` on doc_tcli01, p1_f8, p1_f14, lambda_f14, t5_f8 (`cmp`, also earlier
`sym=0 edge=0` still equals the committed J). `slide=1` code path off means the first J2 hand-back (md5 28d83c5d...) behaviour.

**Result (29 cases, `_sym` -> `_sym2`).** Hard criteria 0 on all 35 renders (dropped, crossings, overlaps, through-text, on/near objects,
leaders on objects, acute, < 110 bends); `label_gap_min_px` 3.0; `labels_back_over_height` 0 everywhere, `label_back_px_max` per case 0.0-8.6
(p1_f14 and lambda_f14 8.6, t5_f14 8.7; the metric's own limit applies and reports 0 over). 394
moves in total (370 outside, 24 inside).
- tcli01: Ala now ends on the left edge 2.8 px below the corner, on its own line (straight 22 -> 22, length 930.2 -> 931.1, +0.9 px; the
  line is straight but Ala's 2.8 px piece had been counted as no piece at all). Corner-end leaders 20 -> 18 of 24.
- Sum: `leader_len_sum` 239851 -> 239375 (0.998), straight 1187 -> 1187, lopsided clusters 61 -> 58 (t5_f14 5 -> 3, p1_f14 7 -> 6), max turn
  unchanged (largest 69.6), canvas changes <= 6 px (p1_f14 1155 -> 1149), toward axis 53 -> 53, back > 1 px 211 -> 208.
- Leaders ending away from a text-box corner (end > 3.7 px from all four): 330 -> 347 of 2267.
Table: `work/J/table_sym2_md.txt` (`cmp2.py _sym _sym2 md`). Pairs `_sym | _sym2` in `work/J/view/*_sym2_pair.png` (tcli01, mito, pbr322, lambda,
p1f8, chl, p1f14_top, p1f14_bottom).

**What did not improve.** The gain is small because most leaders already sit at the closest spot: 1920 of 2267 still end within 3.7 px
of a corner, as the label hugs the ring or the stack there (box pinned by its neighbour at gap 3 px, or by the ring clearance), so the line
cannot meet the edge anywhere else without moving the box away. lambda_f14 (+10 px), p1_f8 (+7 px), pbr322 (+7), pbr322_f10 (+9) get a
little longer: the 1.5 px allowance used to turn a stub into an on-ray end, which the cost prefers. p1_f14 and t5_f14 gain most
(-34, -318 px) where a stack could close up. Fans on doc_tcli01 stay balanced (lopsided 0, max turn 15.9 unchanged).
Still not done: sliding the anchor on the feature and leaning straight lines (needs the additive scaffold edit in `labels/circular.py`).

## LP-D14 probe: why tRNA-Ala (and its fan) cannot move toward the ring, and the `shift` option

Script `work/J/probe_ala.py` (no code change; renders the case in-process with the `_sym2` layout and applies the real tests to hand-made
boxes). `PYTHONPATH=<worktree> python3.13 -P probe_ala.py [case] [label]`. Numbers below are doc_tcli01, R 430.0, box gap 3.05, `lpad` 2.4, `lgap` 2.5.

**Ala alone.** Its leader is anchor->S 61.5 px (fixed: S is on the ring edge R+stub) plus the run S->end 3.30 px; the box corner sits exactly
on the ring clearance (margin 0.00). The line enters the left edge 2.79 px below the top.
| move | run | first test that fails (margin) |
|---|---|---|
| (a) box 2 px left, same y | 0.94 | ring: box 1.70 px inside R+padw (the box is on the clearance already) |
| (a) 4 / 6 / 8 / 10 px | 0.0 | ring 3.4 / 5.1 / 6.8 / 8.5 px inside; also stub S in the grown box; Asn's leader to the box 1.94 / 0.83 / 0 px (< 2.4) from 4 px on |
| (b) 2 px along the own line | 1.30 | ring 2.00 px inside; Trp's box 2.99 px (< 3.05) |
| (b) 4 / 6 / 8 / 10 px | 0.0 | ring 4.0 / 6.0 / 8.0 / 10.0 px inside; Trp's box 1.92 / 0.86 / 0 / 0 px and Trp's leader (< 2.4) |
So the ring stops Ala first (not Asn's leader, which is the second failure from 4 px; Trp stops the along-line slide a hair later). A pure left
shift also moves the entry point UP the edge (the line runs down to the right), toward the corner, not down it. Entering mid-edge
makes the leader longer: entry 4.0 / 5.5 / 7.0 / 9.0 px below the top needs run 3.95 / 4.75 / 5.5 / 6.55 with the ring and stub tests only
(now 3.30), and 3.95 / 6.55 / 9.35 / 13.1 once Trp's box and leader are included (the box must move away from the ring, up into Trp).
Ala's possible gain: at most 3.3 px, 0.7 px if the corner stays on the clearance; the picture cannot change for Ala.

**The fan (per label, only that box moves, dx 0..24 toward the ring, dy -6..6, others fixed).** Run now -> best, turn kept:
Trp 21.64 -> 8.21 (dx 13), Ala 3.30 -> 3.30, Asn 10.53 -> 6.02 (dx 4), Cys 9.95 -> 9.95, Tyr 21.68 -> 21.68; sum 67.1 -> 49.2.
Cys and Tyr are bound by the stack: their lines enter the top edge already; to shorten the run the box must move UP, and the box above
(box height 13-17 px + 3.05 px gap, against about 10 px between the stubs) is there. Their ring margins (9.4 / 20.4 px) are of no use.
Trp and Asn were not found before for two reasons: the old candidates are built per height, not by sliding a box sideways, and the
cheapest candidate at a height was the straight one, which the "turn within 1 degree" test of the slide pass then rejected.

**`shift` (new option, default 0).** In the slide pass each label also tries 1..26 px toward the ring, at its own height or up to 2 px up
(the gap to the box above stays under 6.1 px), its leader keeping direction and turn; the same `try_at` tests decide (3 px to objects,
label gap, `back`, no crossing). Default options reproduce `__J_sym2.svg` byte for byte (cmp on tcli01 and p1_f14, after the edit).
Diff to `_sym2` code: `logs/j2-shift.diff`. tcli01: Trp 12 px toward the ring (leader 42.2 -> 29.8), Asn 9 px (72.0 -> 67.5); the two others of the fan stay.
Six cases, `_sym2` -> `_sym3` (`opts shift=1`): hard counts 0, gap 3.0, `labels_back_over_height` 0.
| case | leader_len_sum | straight | lopsided | max turn | canvas | label_back_px_max |
|---|---|---|---|---|---|---|
| doc_tcli01 | 931 -> 914 | 22 / 22 | 0 / 0 | 15.9 / 15.9 | 1215x934 same | 4.5 |
| mito_f14 | 1373 -> 1363 | 31 / 31 | 0 / 0 | 17.5 / 17.5 | same | 7.6 |
| pbr322_f14 | 1036 -> 1028 | 31 / 31 | 0 / 0 | 5.3 / 5.3 | same | 4.8 |
| lambda_f14 | 4284 -> 4272 | 50 / 50 | 2 / 2 | 38.2 / 38.2 | same | 8.6 |
| chloroplast_hz | 2150 -> 2137 | 102 / 102 | 1 / 1 | 31.8 / 31.8 | 1328x899 -> 895 | 4.4 |
| p1_f14 | 8187 -> 8163 | 57 / 57 | 6 / 6 | 69.4 / 69.4 | 1575x1149 -> 1148 | 8.7 |
| sum | 17960 -> 17877 (0.995) | 293 / 293 | 9 / 9 | | | |
Without the 2 px up (same height only) the gain was 31 px of 17960 and p1_f14 gained one lopsided cluster, so that variant was dropped.
Pairs `_sym2 | _sym3`: `work/J/view/{tcli01,mito,pbr322,lambda,chl,p1f14_top,p1f14_bottom}_sym3_pair.png`.
**Verdict.** Ala cannot gain more than 0.7-3.3 px: the ring clearance sets it (corner on R+padw), S and the feature->S piece are fixed. Cys and Tyr
are bound by box height + 3.05 px gap against the stub spacing. Trp and Asn move (-12, -4.5 px), a visible shift of Trp in the picture, 0.5 % of the length overall.
A joint group re-seat has nothing more to give in this fan: all five together can reach 49 px of runs against 67 only through Trp and Asn.
