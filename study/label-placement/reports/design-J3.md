# J3: where a leader ends (option `face`, default 1; LP-D15, LP-D16)

Code: `gbdraw/labels/proto/design_j.py` (worktree label-study, uncommitted; diff `logs/j3-face.diff`, ruff clean).
Renders: `work/J4/*__J_end.*` (29 cases + `_end_s1` six cases + `_end_strt0` partial). Crops: `work/J4/view/`.
Table per case: `logs/end-table-J3.md`.

## Result in one line
The three measures do NOT reach 0 on all 29 cases: far half 158 -> 27, unclear 555 -> 4, ink gap over 5 px 1,692 -> 31
(largest gap 12.1 -> 9.1 px). 20 of 29 cases are 0 / 0 / 0. The residual is on dense maps (below).

## Rule as built (face=1)
- Facing half: a candidate end is accepted only if it passes the metric's own test (`static` in the seat pass): from above it
  lies between the top vertex and the middle of the side it meets, and so on; exempt within 10 deg of the perpendicular or 1 px
  of the middle. The test is in candidate generation (greedy seat per stack, PUSH/NEAR sequences, pending exit rays), not a filter.
  Side-hit candidates (end on the left/right side for a leader from above) and a finer q0 grid were added.
- One gap: each end is snapped by sphere tracing along its own line to a fixed distance (`lpad`-equivalent, INK_PAD) from the
  ink of its own text. Calibrated on the cairosvg raster: gap median 3.4 px (all sparse/medium maps 3.1 to 4.5), target 2.5-4.5.
- No unclear end: `ink_clear` rejects an end whose nearest other-label ink is closer than own ink + 1 px (the usual
  cause, an end on a top/bottom side of a stacked label, is excluded by the facing rule plus this test).
- `relax` fallback (default 40 px): if no candidate satisfying the rule exists within 40 px of push, the label is seated with the
  rule relaxed. This is why dense maps keep residuals; sparse maps never reach it.
- LP-D16 `straighten` (default `strt=1`): when a leader must bend because its neighbour's box blocks its straight ray, slide
  that neighbour outward along its own leader (direction and facing-half end kept), repeat down the stack, accept when the
  cluster cost (run + tw x turn, symmetric-fan term) falls.

## Ink model
Per text: advances from fontTools (`LiberationSans-Regular.ttf`), glyph masks from PIL, distance field by scipy
`distance_transform_edt`, calibrated to the 3x cairosvg raster that `lp_metrics.end_gaps` uses (the metric needed no change).
Text is placed as the SVG does (`dominant-baseline="middle"`).

## Table (29 cases, 2,267 leaders), `_sym2` -> `_end`
| | _sym2 | _end |
|---|---|---|
| leaders_end_far_half | 158 | 27 |
| leader_ends_unclear | 555 | 4 |
| leader_end_gap_ink_over5 | 1,692 | 31 |
| lopsided clusters | 58 | 92 |
| straight leaders | 1,187 | 1,249 (+5.2 %) |
| leader length sum, px | 239,375 | 253,026 (x1.057) |
| canvas height sum, px | 31,264 | 31,419 |
Hard criteria: check_hard `--end` 0 failing except the three end fields on the 9 residual cases; label gap 3.0 px and
`labels_back_over_height` 0 on all. Length is above the 2 % limit (x1.057): the facing rule makes some leaders longer
(Chloroplast tRNA-Ile moves to the lower left, `view/chloro_end_pair.png`). Lopsided clusters rose 58 -> 92 (above the
brief's "keep"); the symmetric fan term is weaker under the facing rule. Both are open points.

## Residuals (far / unclear / over5, `_end`)
epifagus_f14 0/1/2; lambda_f14 0/0/1; mg1655_rna_f8 4/1/3; p1_f14 1/0/3; p1_open_f14 0/0/2; plasmid153_f14 1/0/2;
t5_f14 15/0/15; t5_f8 6/0/1; tobacco_f14 0/2/2. Reason: pole labels in cone-limit fans (70 deg) with pending exit-ray
reservations leave no spot that is both on the facing half and within the gap; the `relax` fallback then trades the rule for
position. Crop: `view/p1f14_top_end_pair.png` (top of p1_f14), `view/t5_top_end_pair.png`.

## LP-D16 (neighbour slide) 
- doc_tcli01 fan (`view/tcli01_3panel_end.png`: _sym2 | face=1 without the pass | with it): straight leaders 18 -> 21 of 24 there.
  Trp (about -23 deg), Cys (4.4 deg) and Tyr (13.8 deg) remain bent: a joint three-body slide was not found, and for Cys/Tyr
  the box above blocks it. Ala/Asn: Asn is now straight.
- Over the 23 cases that finished `strt=0` (partial, batch cut at hand-back): straight 874 -> 947 (+73 leaders of 1,308),
  length 70,633 -> 70,445 px (-0.3 %), lopsided clusters 33 -> 55 (worse), runtime up (about +9 s on t5_f8).
  So the pass straightens leaders for free in length, but it raises lopsided clusters; keep as option `strt`.

## `shift`
Stays an option (default 0), not part of face=1. Six cases `_end` -> `_end_s1`: length -1.0 % (18,423 -> 18,240), straight 300 -> 299,
lopsided 16 -> 17; all three end measures unchanged. tcli01 (`view/tcli01_shift_end.png`, _sym3 | _end | _end_s1):
tRNA-Trp comes 12 px nearer the ring and ends on the left side's upper half. Gain too small to justify a default.

## Byte identity
`cmp` (final code): `face=0` equals `work/J/<case>__J_sym2.svg`, `face=0,shift=1` equals `work/J3/<case>__J_sym3.svg` on
doc_tcli01, p1_f8, p1_f14, lambda_f14, t5_f8.

## Production needs
Font-derived ink without scipy (fontTools glyph outlines to a signed distance or per-glyph boxes, or a coarse
mask), the same placement of text as the renderer, runtime budget (dense maps 20-65 s in 4-parallel vs about 10 s for J2; max
single render 42.6 s), and a cap on the seat search for pole labels.

## Open points
1. Dense-map residuals (t5_f14 15/0/15): needs a cone wider than 70 deg for pole labels or a second ring of labels.
2. Length x1.057 and lopsided 92: decide the trade with the Owner (face=0 / strt=0 / lower `relax`).
3. `strt` full 29-case comparison not completed (23 of 29 finished).
4. Cys/Tyr in tcli01 stay bent.

## LP-D18: a head-on leader ends at the middle of the box end (option `mid`, default 1; face=1 only)
Owner 2026-10-12 on `tcli01_3panel_end.png`: tRNA-Trp should move down a little and the leader come out of the middle of the left end.
Code (diff `logs/j3-mid.diff`, 5 lines, in `_Par` and `_Lab.static`): a candidate whose leader (stub to end, at least 3 px) is within
`midang` (11 deg; the metric says 10) of the horizontal and ends on the left or right end of the box more than `midfr` (0.30; metric 0.35) of the half
height from its middle is refused. It is a hard test inside `static`, so seating, symmetrize, straighten and slide all obey it. The box
moves along the stack through the existing pushes (1 px steps); the leader does not bend toward a corner. The `relax` fallback (face=0 for that
label) is outside the rule, as it is for the other end rules.
`mid=0` reproduces `_end` byte for byte (`cmp` on tcli01, mito_f14, p1_f14, lambda_f14); `face=0` still reproduces `_sym2` on the same four.
Result over the 29 cases (`_end` -> `_mid`; per case `logs/table-mid.md`, `logs/end-table-J3b.md`): off middle 109 -> 2; far half 27 -> 27,
unclear 4 -> 4, ink gap over 5 px 31 -> 30; straight 1,249 -> 1,228; length 253,026 -> 254,311 (x1.005); lopsided 92 -> 101; canvas height sum
31,419 -> 31,438. Against `_sym2`: straight 1,187 -> 1,228 (+3.5 %), length x1.062, lopsided 58 -> 101. Hard criteria: check_hard `--end`
fails only on the same 9 residual cases as `_end` (their dense-map fields) and tobacco_f14 `off middle` 2; label gap 3.0 and labels_back_over_height 0 everywhere.
The two left on tobacco_f14 (tRNA-Leu -1.38, ribosomal protein L32 -1.27) are `relax` labels (3 relaxed in that case).
tcli01 (`view/tcli01_mid_pair.png`): Trp now straight and ends at the middle of its left end (box 8 px lower, same x). The whole fan below moves down 7 px,
so Asn (straight in `_end`) turns 16 deg and Tyr 25 deg (was 18); Trp, Ala, Cys are straight (3 of 5 as before). Straight count in this crop 21 -> 20, length 977 -> 1,020.
Why Asn cannot stay straight: its straight ray passes through the lowered Ala box, and a neighbour slid along its own leader keeps the same perpendicular
offset to the parallel ray, so the LP-D16 slide cannot clear it; the middle end of Trp costs the stack those 8 px. Cost/benefit is the Owner's call (`mid=0`).
Seen in the pairs: mito_f14 (tRNA-Leu, NADH dehydrogenase, tRNA-Ile), chloroplast (psbD, psbC, psbZ, rps14), lambda (N-like), pbr322 (NruI) all end at the
middle now. Worse: boxes move down up to 12 px (psbD), so some gaps in a stack open (pbr322 EagI/NruI 4 -> 22 px apart), and a few short leaders get a small
hook at S to meet the middle (chloroplast tRNA-Gly, rps14).
Continuous rule (aim at the centre of the first square of the box): evaluated, not built. For a 30 deg arrival it puts the end at 0.58 of the half height
above the middle, not at the middle the Owner asked for on Trp (Trp's final leader arrives at 30 deg and ends at the middle), and it would be a hard constraint on every fan
(dense maps already depend on `relax`). The head-on band with the facing-half rule covers the Owner's three remarks; `mid=2` could add it as a soft cost.
Lopsided 58 -> 92 -> 101 and length x1.062: lopsided comes from `face` (e.g. gal_hmmt_atskew 0 -> 3, length x1.23: the symmetrize pass cannot re-seat a
shifted cluster under the facing and ink rules); `mid` adds 9. No cheap fix found; not pursued. Runtime: max 46 s per render (4 parallel), tobacco_f14 30 s.

## Straight fan (round 3; option `fan`, default 1; face=1 only)
Owner LP-D16 / LP-D19: a straight fan, bends only at the outside of a cluster, never a bent leader between straight ones.
Code (`logs/j3-fan.diff`, 154 lines, `_Eng.fan_seat`, `_Eng.fan`, `_broken`, run last in `place()` after `slide`): for every chain of one stack
(boxes at most gap + 3 px apart, the lop.py cluster) with a bent leader, the run (leaders at most 60 deg steep, so no pole block) is seated again as a straight
fan. Each label that is not already straight and within 0.35 of the half height of the middle of its end is placed along its own first piece (anchor -> S
extended), the end at the middle of the near box end (then 0.35 and 0.25 of the height from the facing side), at the smallest run (0.5 px steps, at most `fanx` = 40 px
longer than before) that passes `static` and `try_at` (all hard and ink tests) against everything placed. Order: from the equator outward and, tried second,
from the outermost label inward. The longest prefix that stands (each label re-checked against all others) is kept when its labels stay inside the old extent of the
labels, it has more straight leaders (turn under 3 deg) than before, and the chain does not become mixed (a bent leader between two straight ones, turn limit 1 deg)
when it was not.
`fan=0` == `_mid`, `fan=0,mid=0` == `_end`, `face=0` == `_sym2` (cmp, tcli01, mito_f14, p1_f14, lambda_f14).
doc_tcli01 (`view/tcli01_fan_3panel.png`: _sym2 | _mid | _fan): Trp, Ala, Asn, Cys, Tyr all straight, parallel, every end at the middle of the left box end; +65 px of leader on the fan
(the arithmetic said +69); the box column ends up almost aligned. mito_f14, pbr322_f14 the same: 0 mixed clusters.
Over the 29 cases (`_mid` -> `_fan`, `logs/table-fan.md`): straight 1,228 -> 1,274; length 254,311 -> 255,541 px (x1.005; x1.068 against `_sym2`); lopsided 101 -> 90; canvas height
sum 31,438 -> 31,438 (no growth by construction); far half 27, off middle 2, unclear 4, ink gap over 5 px 30: all unchanged; check_hard --end: the same 9 dense-map cases fail as in `_end` / `_mid`.
Mixed clusters (`scripts/fan_order.py`): J2 50 of 115, `_end` 57 of 122, `_mid` 62 of 126, `_fan` 54 of 112. The target (clearly below 50, 0 on sparse and medium maps) is NOT met:
0 mixed on doc_tcli01, tcli03, tcli09, gal_hmmt* (4), mito_f14, mito_open_f14, pbr322 (3), plasmid37_f14, but epifagus_f14 4, chloroplast 5, lambda_f14 2, lambda_open_f14 1,
plasmid37_open_f14 1, mg1655 5, n4 2, p1_f14 1, p1_f8 4, p1_open 5, plasmid153 7, t5 5, tobacco_f14 4, tobacco_f8 8 are still mixed.
Cause (lambda_f14, "NinE family protein" after the straight "NinD protein"): the straight ray of NinE cannot be seated at any run up to +40 px. For runs under 14 px another
leader or the ink test refuses it; from 14.5 px its leader passes within 2.4 px of the box of NinD (x -376.6..-299.1, y -277.9..-260.5), which is fixed (straight, near the middle) and hangs half a
box height (8.7 px) below its end, while the stubs of the two are closer than 8.7 + 2.4 px. A straight neighbour pair is only possible when the next label ends before the box above (staircase) or the stub spacing is at least
hang + 2.4 px; in tight clusters (stubs under about 11 px apart) the straight run cannot be extended, so the cluster needs the bend ramp of the symmetric fan. Re-seating the fixed labels too, or seating from the other end, changed nothing on
lambda_f14, epifagus_f14 and chloroplast. Not built: `fan=2` (a bend ramp at the two ends with a straight middle, the middle run chosen from the straight labels already there, bends not decreasing outward with
opposite signs): it needs the `want` lean spots of `symmetrize` driven by the mixed test, plus a rollback of a cluster whose later pass breaks the order.
Seen at 3x: tcli01 and mito_f14 are the Owner's picture (straight parallel fan, ends at the middle); lambda_f14 and chloroplast crops show no change (their mixed clusters stay); the short hooks (chloroplast tRNA-Gly, rps14) remain. Worse: leaders of the fan grow (tcli01 +6 %, mito_f14 +4 %), and
the boxes line up in a column that stands 5 to 23 px farther from the ring than before. Runtime: max 37.8 s per render (4 parallel).
