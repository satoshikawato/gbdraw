# Prototype P (LP-D22 to LP-D25): ink boxes, 12 contact points, one straight-first cluster construction

Code: `gbdraw/labels/proto/design_j.py` in the label-study worktree, option `port=1` (commit 6e10e44c on `study/label-placement`,
local only; diff `logs/p1.diff`, 462 lines vs 9e51a4ea). Renders: `work/P/<case>__J_P.*` (29 cases; `*__J_Pv.*` = tcli01 variant
`pslide=12,pvert=14`). Table: `logs/table-P.md` (built by `scripts/table_p.py`). Pictures: `work/P/view/` (made by `work/P/make_views.py`).
Run: `page/run_extra.py --out work/P --opts port=1 --suffix _P CASE:J`.

## What is built
- R1 box: with `port=1` every Item gets w,h of its ink box (`_ink_rel`, same font model as `scripts/lp_ports.py`); `_emit_ink` shifts the
  text origin so the ink fills the engine box. `aabb_local` / `oriented_corners` = ink box; `port` written per label. Gap 3.05 between ink boxes.
- R2 contact points in the geometry model: `_Lab.port_end` (leader S -> point, stops `lgap` = 2.5 px before it, NaN where the arrival
  rule of `port_check.arrival_ok` fails), `_Lab.cands_port` (all 12 points x box positions base + OFFS, plus the point on the feature's own ray
  at the row). `cands` returns these when `port=1`; face, mid, strt, fan, sym, slide are forced off. `static` adds: a bent leader has a last
  piece >= `minlast` (6 px) (no hooks).
- One cluster construction replaces symmetrize/straighten/fan: `_Eng.cluster` (chains as in `fan`, a lone label is a chain of one):
  `chain_nodes` (all seats within `pwin` = 30 px rows, every point; straight seats on the feature's own ray from the nearest allowed one to
  `pslide` = 8 px farther out (R5); a seat is straight (turn < 0.5) or bent (>= 3); cost = run + tw x turn + `bend` 6 (+ `pvert` 4 at a vertex),
  tested against all labels outside the chain), `chain_pair` (order, ink gap, leader clearance both ways, crossing), `chain_dp` (cheapest choice
  with the pattern bent* straight* bent* along the stack = R3, bent labels sit at the least row distance from the inner neighbour because
  the cost is the leader). The label's present seat is one of its seats. Result kept only if every label passes every test against all others
  (`legal`, the try_at rules); else the chain keeps its seats (25 chains of ~900, counted).
- Seating pass unchanged except: the way kept free for a pending label is a leader a contact point can take (`pcone` 44 degrees, was 68.5),
  and a label with no contact-point seat within `relax` (40) px of push takes the J2 seat (counted `relaxed`; this is the rule that lost, see below).
- `try_at` was split into `try_at` + `legal` (same operations; `face=0` still reproduces J2 byte for byte: `cmp` on doc_tcli01, mito_f14, p1_f14,
  lambda_f14 verified by running).

## Numbers, J2 -> P, 29 cases (verified by running; per case in table-P.md)
| | J2 | P |
|---|---|---|
| leaders drawn | 2,267 | 2,267 |
| straight (metric) | 1,187 | 1,083 |
| leader length sum, px | 239,375 | 198,244 (x0.83) |
| canvas height sum, px | 31,264 | 29,795 |
| lopsided clusters | 58 | 175 |
| mixed clusters (fan_order) | 50 of 115 | 37 of 212 |
| off a contact point | 1,971 | 423 |
| wrong direction | 56 | 7 |
| far half (old LP-D15 measure) | 158 | 86 |
| ink gap over 5 px (old measure) | 1,692 | 469 |
| bent leaders, last piece median / max (px) | 76.3 / 804 | 46.7 / 571 (under 6 px: 48 -> 0) |
| label gap min | 3.0 | 3.0 (3.05 by construction) |
| check_hard (no --end), after lp_objects | 0 failing | 0 failing |
Seconds sum 156 -> 122. `check_hard --end` fails every case for both (the end measures are the superseded LP-D15 ones).

## tcli01 (five tRNAs; table5, verified by running)
P: stack height 74.6 px (J2 90.3, Owner 76.0), leader sum 320 px (333.5, 328.5), Trp straight into L (not bent), Ala/Asn/Cys straight but
into the TL vertex, Tyr bent 5.6 deg into TL; bent leaders of the whole tcli01 map: last piece median 9 px, max 23 px (last_piece.py). Whole map: off contact point 0, wrong direction 0.
Why TL, not L: for the straight ray at 32 degrees the L point lies 10 px farther along the ray than TL, more than the 8 px slide of R5, and
costs more than the `pvert` 4 px preference. `pslide=12,pvert=14` (`work/P/view/tcli01_P_vs_Pv_ports.png`) puts only Ala into L.

## Rules that conflicted
- R2 (hard) lost to "nothing may fail to be seated" in dense stacks: 466 labels take the J2 seat after the relax test (423 off a contact point
  in the table; dense maps only: t5 116+78, plasmid153_f14 72, mg1655 57, tobacco_f14 30, p1_f14 18, tobacco_f8 22; 0 off on tcli01, tcli03,
  tcli09, gal_hmmt* (4), mito*, pbr322*, plasmid37*). Medium maps: epifagus 1, lambda 1+1, n4 4, chloroplast 7.
  Cause seen on tcli01 His: a neighbour box 3.6 px from the stub leaves no leader that a contact point can take.
- R5 (slide <= 8) against the Owner's picture "straight into L": see above.
- R3 vs R6/hard tests: 25 chains keep their J2-style seats (check failed), and the cluster is not a guarantee: mixed clusters remain.

## What is not met / not done
- Mixed clusters 0 on sparse and medium maps: not met (epifagus 2, chloroplast 5, mg1655 4, p1_* 4 to 5, tobacco 4 to 5, plasmid153_f8 4); 0 on tcli*, gal_hmmt*, mito*, pbr322*, plasmid37*, lambda_f14, lambda_open_f14.
- Straight count fell (1,187 -> 1,083) and lopsided clusters rose (58 -> 175): leaders end at a point, so the box must hang on the ray at an exact place; many ends now bend a little.
- Canvas height larger than J2 on t5_f8 (1517 -> 1547), chloroplast_hz (+5), plasmid37 x2 (+4), mito_open (+1).
- Dense maps: far half t5_f14 36 and mg1655 12 (old measure, relaxed labels).
- Not checked by eye: p1_f14 bottom, mito_f14 beyond the crop (pictures exist in `work/P/view/`).

# Round 2 (commit see handback-P2.md; diff `logs/p2.diff`; renders `work/P2/*__J_P2.*`; table `logs/table-P2.md` J2 -> P1 -> P2)

## Changes (all behind `port=1`, so `port=1` now means P2; round-1 renders stay in `work/P/`)
1. Arrival rule (`_Lab.port_end`): within 45 degrees of a short side's normal only the middle of that side; 45 to 70 degrees a vertex
   (from outside its corner) or a long-side point; steeper than 70 degrees long-side points only. Checked with the new `scripts/port_check.py`.
2. Relaxed labels: `chain_nodes` offers the present box with its last piece re-aimed at every point its direction may use (static and legal
   tests), and a seat that is at no point costs `noport` = 400 px, so a chain keeps such an end only when nothing else passes.
3. Clusters: chains now use the fan_order/lop definition (boxes at most 6.1 px apart); up to 3 rounds, rounds 2 and 3 only for chains that read mixed
   and are new. A chain whose result fails the all-against-all check freezes the failing label (and then the labels it collides with) to the present seat and
   runs the DP again. The present seat counts as straight below 1 degree (the fan_order criterion), nodes below 0.5.
4. Kinks (soft): a bent seat costs `kink` (8 px, fading to 0 at 8 degrees of turn), so a slide along the ray (within `pslide`) is preferred.

## Numbers, 29 cases (verified by running): J2 -> P1 -> P2
| | J2 | P1 | P2 |
|---|---|---|---|
| straight (metric) | 1,187 | 1,083 | 1,019 |
| leader length, px | 239,375 | 198,244 | 201,288 |
| canvas height sum | 31,264 | 29,795 | 29,770 |
| lopsided clusters | 58 | 175 | 182 |
| mixed clusters | 50 of 115 | 37 of 212 | 14 of 230 |
| off a contact point | 1,971 | 423 | 307 |
| at a point the NEW rule does not allow | 118 | 350 | 6 |
| turns 0.5 to 8 degrees | 332 | 354 | 273 |
| bent last piece median / max, under 6 px | 76 / 804, 48 | 47 / 571, 0 | 50 / 571, 0 |
check_hard (no --end) after lp_objects: 0 failing on 29. tcli01 map: bent last piece median 7, max 10 px.

## tcli01 (five tRNAs, table5; picture `work/P2/view/tcli01_P2_3panel.png`)
Trp bent 29 degrees, last piece 1.3 degrees off the horizontal, into L, 6.9 to 9.8 px (all bent leaders of the map are 7 to 10 px); Ala, Asn, Cys straight into L
(end 0.23 to 0.28 of the half height above the middle: the leader stops 2.5 px before L on a 32 to 36 degree line); Tyr bent 9 degrees, last piece 46 degrees, into TL.
Stack height 70.0 (J2 90.3, Owner 76.0), leader sum 310 (333.5, 328.5), box lefts 29 px apart (Owner 42).

## Remaining off a contact point (P2) and why
Epifagus 1, lambda 0, n4 1, chloroplast_hz 6; dense: t5_f14 92, plasmid153_f14 66, mg1655 48, t5_f8 46, tobacco_f14 14, p1_f14 14, p1_* 6.
Seen on lambda holin (round 1): a pole label whose leader is steeper than 70 degrees arrived at a vertex, which is not allowed; with `noport` 400 the
construction found a long-side seat for it (lambda 0). The remaining ones keep the J2 seat because no seat passes `static` and `legal` within the rows
of `pwin`: chloroplast ones are inside-ring labels of 19 px width (rps16, psbL, rpl36 ...) in fans where the leader is steeper than 70 degrees and the
only allowed points are the long-side ones, which collide with the neighbouring leaders; epifagus tRNA-Gln, n4 gp72/gp4 are pole labels in the same situation. Not examined further.

## Remaining mixed clusters (14 of 230) and cause
- epifagus_f14 (tRNA-His S, Gln B, L2 B, maturase S, Ile S) and chloroplast_hz 3 clusters (psbF.., psbA.., ndhI..): the all-against-all check fails for a seat set and
  the freeze loop then leaves a chain that cannot satisfy bent* straight* bent* (the present seats are mixed). The DP only checks neighbours k, k+1 (`chain_pair`);
  the failing pair is a label two places away (leader against box or crossing). Not built: a second-neighbour term in the DP.
- p1_f14 2, p1_f8 4, p1_open_f14 3, plasmid153_f14 1, p1_* clusters of 8 to 18 labels: same cause on a larger scale; chains over 40 labels keep only the 250 cheapest seats per label, and
  a chain of 12 to 18 labels at the left edge of a pole block (UpfN .. PmgS) is mixed because its labels are mostly relaxed or frozen.
- 0 mixed on every sparse map and on lambda*, mg1655, n4, plasmid153_f8, t5*, tobacco*; not 0 on epifagus 1 and chloroplast 3 (medium maps).

## Rules that conflicted (round 2)
- R3 lost to the all-against-all check in 14 clusters (above); R2 lost to the same check for 307 labels (dense maps and pole fans).
- The exclusive arrival rule plus the 8 px slide: tcli01 Ala, Asn, Cys fit straight into L because L is the only point their 32 to 36 degree lines may use (no TL competitor), so no conflict there.
- Kink cost (soft): turns 0.5 to 8 degrees 354 -> 273 (J2 332); not hard, still 273.

## Worse than J2 (seen in the pictures and numbers)
- Straight leaders fell (1,187 -> 1,019) and lopsided clusters rose (58 -> 182).
- Canvas height larger than J2 on t5_f8 (+35 px), chloroplast_hz (+6), gal_hmmt_atskew/gal_hmmtdna_atskew (+6), plasmid37 x2 (+2).
- p1_f14 bottom (`work/P2/view/p1_f14_bottom_P2.png`): three long slanted leaders (anti-repressor, host cell division inhibitor) and a wider empty band than J2; p1_f14 top: the first
  rows (Lpa-like activator .. PdcA) sit higher and the leaders of the left stack are steeper.

# Round 3 (LP-D26: follow the circumference; second-neighbour check)

## What was built (all in `gbdraw/labels/proto/design_j.py`, option `port=1`)
- R8, the arc: `chain_arc` seats a cluster (chain of boxes at most 6.1 px apart) on one arc concentric with the ring. For a radius rho (0.5 px grid) every label
  keeps the seats whose leader end is within `arctol` = 0.95 px of rho (neighbours differ by under 2 px, the `arc_fit.py` limit); the radii are tried from the least
  lower bound of the cost (`rhos` = 12), the cheapest choice wins (so the smallest radius the hard rules allow, R5). Straight seats may now lie up to `parc` = 60 px along the ray
  (`far` seats, beyond `pslide` = 8 px), but only inside an arc window; a cluster with no feasible radius takes its seats without the arc and without far seats.
  The rounds also re-treat chains that read off the arc (`_off_arc`), not only mixed ones.
- Second-neighbour check: `chain_dp` is now a DP over pairs of seats (labels k-1, k), so labels two places apart are tested for box gap, leader-vs-box and crossings
  (`compat2`); the order bent* straight* bent* is the phase of the state (`order=False` drops it). Labels farther apart are handled by `chain_search`: when the cheapest choice has a
  colliding pair, the later label is branched on (its 8 cheapest seats, each with all colliding seats of the other labels struck out, or none of them), best first, budget 40 DP runs
  for an arc trial, 300 for the fallback. Seats are thinned to the cheapest of each kind of end (`seats` = 40, 120 in the fallback; the present seat stays).
  A chain with no seats in the order of R3 is seated as it comes (reads mixed) rather than keeping its old seats. The freeze loop of round 2 is gone (one owner: chain_search).
- Result of the check: `kept_check` (a chain that fails the all-against-all test) is 0 on all 29; `kept_dp` (no choice at all, the chain keeps its seats) is 1 chain over the 29 (not identified; round 2: 14 chains kept mixed).

## tcli01 (five tRNAs, Leu/Ser/His, Ile/Gln/Met)
- arc_fit: 0 of 3 side clusters off the arc (steps 0.3 to 0.9 px; P2 8.3 / 7.5 / 5.4, J2 0.8 / - / 17.6). Leu bent, Ser and His straight; Gln straight, Ile and Met bent.
- Trp to Tyr: Trp bent -44 deg (a horizontal last piece, 10 px, into L), Ala -25 deg (7 px), Asn straight, Cys 9 deg, Tyr 9.7 deg (into TL). So bent, bent, straight, bent, bent: only ONE straight leader, not the
  two the Owner's reading predicted (Asn and Cys): Cys and Tyr stay 9 to 10 deg bent (just outside the 8 deg kink range) because the five rows (53 px) do not fit the rays on the arc.
  Stack height 70, leader sum 980 px for the whole map (P2 976).
- Last piece of bent leaders in tcli01: median 7, max 19 px (the label with the 19 px piece was not identified); the 6 to 14 px target is missed by that one.

## p1_f14 bottom (item B)
P3 leader-end radii of the three labels and the "anti-repressor" above: 461 / 477 / 490 and 449 px (P2 483 / 499 / 513 and 451; J2 490 / 510 / 529 and 473). All straight (turn 0, -0.7, 0, -3.3 deg).
What stops them: they are the lower part of ONE column of six labels (replication protein, DUF3800, anti-repressor x3, host cell division inhibitor) packed with the 3.05 px box gap
(boxes at 409.8, 425.9, 443.6, 459.7, 475.8, 491.9, pitch 16.1 px); the column hangs under PmgF (box bottom 406.8), whose leader is straight from its own ray, so the three cannot
rise unless the column above rises (hard rules: box gap, order along the stack, leaders not crossing). They were 20 px lower in P2 because the P2 DP had no seat set that fitted the whole
column with straight leaders (the round-2 freeze loop kept relaxed seats); in P3 the pair DP finds the column. Not a construction defect any more.

## Breakdown of the 2,267 leaders (scripts/breakdown_p.py)
J2: straight 1017 | slight bend 332 | lone bent 3 | bent at the outside of a straight run 431 | no straight in the cluster 219 | bent between straights 265.
P2: 975 | 273 | 20 | 501 | 472 | 26. P3: 859 | 282 | 65 | 479 | 542 | 40.
Straight seats not taken (inferred, not measured per label): the 542 leaders of clusters with no straight leader are mostly dense maps (t5, plasmid153, mg1655, tobacco, p1: clusters of 8 to 40 labels in
a pole fan where the rays are closer than the rows) where R8 asks the ends on one arc and the rows spread along it; the lone bent leaders (65, J2 3) are labels whose own ray radius cannot be matched by a seat
with a short last piece. The 1,150 target was not chased (R8 lowers the straight count by design). Slight bends 0.5 to 8 deg: 282 (target under 100 not met; the cost `kink` fades to 0 at 8 deg and the arc windows
(+-0.95 px) pick some).

## Remaining off a contact point and mixed (P3)
- Sparse maps: epifagus 0 / 0 mixed, n4 0 off, lambda_f14 and lambda_open 2 off + 1 mixed each (worse than P2's 0), chloroplast_hz 3 off + 2 mixed (P2 6 + 3), tcli*, mito*, pbr322*, plasmid37*, gal_hmmt* 0.
  Lambda: the 14-label NinD .. Bor-like chain has no feasible arc radius and no seats in the order of R3; the fallback `order=False` takes it, and holin/anti-holin and Bor-like keep their relaxed J2 seats because no seat within pwin
  = 30 px of their rows is free (pwin = 60 fixes lambda to 0/0 but chloroplast to 5 off and canvas +31 px; not adopted). The mixed one is that chain.
- Dense maps (off a contact point, P2 -> P3): t5_f14 92 -> 103, t5_f8 46 -> 70, mg1655 48 -> 47, plasmid153_f14 66 -> 57, plasmid153_f8 1 -> 1, tobacco_f14 14 -> 17, tobacco_f8 6 -> 6, p1_f14 14 -> 8, p1_f8 6 -> 2, p1_open 6 -> 3.
  Cause: labels relaxed to the J2 seat in `seat` (no contact-point seat within 40 px of push in a fan beyond its geometric limit) and kept because their chain has no free seat within pwin. t5_f8 rose by 24 (cause not examined; probably more chains falling back to unconstrained seats when no arc radius is feasible, inferred). Mixed on dense maps: p1_f14 2, p1_f8 2, p1_open 0, plasmid153_f14 1, plasmid153_f8 1, t5_f14 1, tobacco_f14 2 (all in unconstrained fallback chains).
- Totals: mixed 13 of 257 (P2 14 of 230), off a contact point 321 (P2 307), wrong direction 6, check_hard 0 failing on 29.

## Rules that lost where (round 3)
- R8 (arc) lost to the hard rules in 73 of 158 side clusters, mainly dense maps; it holds on all of tcli01 and on the sparse side clusters that have a feasible radius.
- R3 (order) lost in 13 clusters (fallback without order, listed above) and R2 in 321 labels (relaxed seats).
- R5/R6: arc windows make the leaders of a few clusters longer (the straight middle label slides out); total leader length 203,920 px = x0.852 of J2 (P2 x0.84), within the x0.90 target.

## Worse than P2 (pictures `work/P3/view/`)
- Chloroplast cluster matK / tRNA-Lys / rps16 / tRNA-Gln / tRNA-Ser (`chloroplast_matK_cluster_P2_P3.png`): rps16 and tRNA-Ser were pulled out to the arc radius of the cluster: their leaders are 40 px long straight (P2: 12 px) and
  the rows leave a wide gap. Mixed there is solved, tidiness is not.
- Straight leaders 888 (P2 1,019; J2 1,187), lopsided 221 (P2 182; J2 58), canvas taller than P2 on tobacco_f14 (+41 px), n4 (+31), plasmid153_f8 (+24), p1_f8 (+19), p1_open (+16), tobacco_f8 (+12).
- Run time 2 to 30 s per map (P2 under 20): p1_f14 29 s, t5_f14 17 s.

# Round 4 (even rows, smooth arc envelope; unfinished)

## What was built (`design_j.py`, `port=1`)
- R8 as a smooth envelope instead of one radius per cluster: `chain_pick(..., arc=True)` keeps a pair of neighbours only when their leader ends differ by at most `arcstep` - 0.05 = 1.95 px in distance from the centre (`rad` of each seat); the radius drifts along the cluster and the cost (run + turn) puts each label as near the ring as that allows. `chain_solve` tries (order, arc) = (yes, yes), (no, yes), (yes, no), (no, no); the last two use the soft cost `arcw` = 4 px per px of step beyond 2 (`chain_pen`). `chain_arc`, `rhos`, `arctol`, `far`/`pslide` are gone (the slide limit is `parc` = 60 px for all straight seats; the cost keeps them short).
- Even rows (R4): `chain_pen` adds `loosew` = 1 px per px of box gap beyond `rowgap` = 3.6 px (at most 16 px counted) to a pair of neighbours when one of them is bent; and a convex bend cost `tsq` = 0.01 px per square degree shares a displacement between the two ends.
- Found while debugging: the seat thinning of round 3 deduplicated seats by leader end (4 px cell) alone, so two different boxes with the same end counted as one and the cheaper one hid the other. The key is now box and end (1 px). The P3 solution of Trp to Tyr was in the seat set but was cut by that. (Round 3's tcli01 picture came from the chain that no arc radius fitted and the fallback happened to be on the arc.)
- `scripts/arc_fit.py` has the count `loose` (neighbours in a column, both ends within 60 deg of the horizontal, box gap over 5 and at most 20 px, one leader bent); `scripts/tcli_dump.py` prints the tRNAs of tcli01.

## Result and what is not done
- doc_tcli01 (kept, better): 0 of 3 side clusters off the arc; Leu / Ser / His (left side of the map) rows are even like J2 (gaps 3.7 to 4.3 px, Leu +16 deg, Ser straight, His -9 deg); Ile / Gln / Met even; Trp to Tyr smooth (bends -36, -13, 5, 8, 14 deg with 6 to 8 px last pieces in the seats=5 run; the batch run (seats 2) gives -41, -6, 5, 8, 9 deg and last pieces 9 to 13 px). No straight leader in the five any more (the convex cost and the arc prefer small bends; slight bends 5 to 8 deg remain at Asn, Cys). `loose` doc_tcli01: 2 pairs (two bent leaders, gap over 5 px), not 0.
- Over the 29 cases: side clusters off the arc 68 of 162 (P3 73 of 158); loose 196 (P3 222, J2 33): not near J2.
- NOT met: mixed 22 of 260 (P3 13), off a contact point 401 (P3 321; plasmid153_f8 1 -> 12, tobacco_f8 6 -> 16, tobacco_f14 17 -> 28, plasmid153_f14 57 -> 84, lambda 2 each unchanged, chloroplast 3 -> 5), canvas not above P2/P3 on n4 (1003 vs 1017, ok), tobacco_f14 (1243 vs 1294, ok), plasmid153_f8 (1057 vs 1080, ok), p1_f8 (927 vs 946, ok), p1_open (962 vs 965, ok), tobacco_f8 (1051 vs 1040, +11 worse than P3); lambda regression not fixed; t5_f8 off 46 -> 70 -> 76 not fixed.
- Cause of the regressions: the pair DP cannot afford the full seat set. With all ~150 seats per label a 5-chain is fine, but a 12 to 14 label chain took 10 minutes (the all-seat run of 5 cases was cut off after 600 s; `seats=5` took 9.5 min for the same 5 cases, chloroplast 568 s). The batch used `seats` = 2 (the two cheapest seats of each 3 px row, bent and straight apart; 1 for chains over 20), which loses seats a neighbour needs: the mixed and off-contact-point counts rose. Not solved: needs a DP that exploits the arc window (neighbouring seats within 2 px of radius) instead of N^3 over all seats, or a better thinning (by radius as well as row).
- Chloroplast matK cluster: rps16 and tRNA-Ser leaders are about 30 and 25 px (P3 40, P2 12) and psbA / matK are pushed lower than in P2 / P3: better than P3, worse than P2 and J2.
- Not looked at: `p1_f14_top_P4.png`, `p1_f14_left_P4.png` were rendered but not inspected.
