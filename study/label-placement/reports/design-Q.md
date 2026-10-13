# Prototype Q (LP-D27): side clusters in arc coordinates

Code: `gbdraw/labels/proto/design_j.py`, option `arcfirst=1` on top of `port=1` (`_Eng.arc_geom`, `arc_static`, `arc_chain`, `_pav_l1`, `_lipschitz`, a pass at the head of `cluster`).
With `arcfirst=0`, doc_tcli01 and lambda_f14 are byte-identical to `work/P4` (cmp, run). Renders `work/Q/*__J_Q.*` (opts `port=1,arcfirst=1`), table `logs/table-Q.md`, pictures `work/Q/view/` (`work/Q/make_views.sh`).

## Built
- Runs: the labels of one column (both halves together) whose boxes are at most `qgap` = 20 px apart, in stack order; every ray within 70 degrees of the horizontal (`qnq` 0.34). Anything else (a run with a pole label, one label) goes to P's chain DP as before.
- Arc radius per label, rho_k: start at the smallest radius at which the label stands straight on its own ray and passes `static`; the sequence is smoothed so that neighbours differ by at most `arcslope` 1.9 px (`_lipschitz`, the smallest sequence above the minima).
- Natural place = the ray meets the circle rho_k at the label's own angle: `arc_geom(L, rho, phi)` gives the leader end at angle phi on the circle, the last piece S -> end, the contact point 2.5 px beyond the end, the point chosen by the arrival rule (L/R within 45 deg, facing vertex 45-70, long-side point beyond), then box and box centre.
- Rows: box centres at least (h_k + h_k+1)/2 + 3.05 apart; least total |displacement| by pool adjacent violators with medians (`_pav_l1`): in a pooled run the middle label is at its natural place (straight), the others at the row distance on both sides (even rows). Even pools take the middle target nearer the mean. The target box centre is turned into phi by interpolation of the monotone yc(phi) table (1401 angles), so rows are exact.
- Checks per label: arrival point allowed, `static`, turn <= `maxturn` 40 deg, `legal` against everything else (bg + the run). A failing label raises its radius by 2 px; run again; stops at `rhomax` 40 px above the minima (the run then goes to P's DP and is counted: `q_kept`). No seat lists; each pass is linear in the run (per-label work is a 1401-point table).
- Not done: the "redo the row distances once" step (the contact point type is fixed by the direction at the target position; a point-type jump of h/2 is not iterated), a second pass after P, pole labels.

## Numbers (29 cases, run)
| | J2 | P2 | P4 | Q |
|---|---|---|---|---|
| straight | 1,187 | 1,019 | 899 | 946 |
| leader length px | 239,375 | 201,288 | 204,325 | 203,256 (x0.849) |
| canvas height sum | 31,264 | 29,770 | 29,750 | 29,750 (no case above J2) |
| mixed clusters | 50 of 115 | 14 of 230 | 22 of 260 | 29 of 268 |
| off a contact point | 1,971 | 307 | 401 | 399 |
| wrong direction | 118 | 6 | 6 | 6 |
| side clusters off the arc | 68 of 131 | 135 of 164 | 68 of 162 | 49 of 168 |
| loose pairs | 33 | - | 196 | 139 |
| turns 0.5-8 deg | 332 | 273 | 324 | 296 |
| seconds sum | 156 | 117 | 242 | 249 (max 34 s, chloroplast_hz) |
check_hard after lp_objects: 0 failing of 29.
Runs built by Q / kept for P: tcli01 4 / 2 (the two kept are pole runs), plasmid153_f14 0 built, t5 1, tobacco_f14 1. Most "kept" runs contain a pole label (nq < 0.34), not a radius failure; radius failures (no solution within 40 px): epifagus (Glu .. S14, a fixed neighbour box blocks the compressed run before qgap 20; fixed by qgap 20 only partly), lambda ead/Ea22 group, p1_f14 TciC/TciB, plasmid153_f14 two long runs (static and turn fail).

## doc_tcli01 (run; `logs`/`work/Q/view`)
Trp..Tyr: turns -36.6, -21.9, 0.0, 1.7, 8.3 deg; one straight leader (Asn), Cys 1.7 deg; bends grow outward; leader ends r = 438.2 to 442.0 (steps <= 3.8 between Trp and Ala, within the 1.9 slope by the smoothing: arc_fit 0 of 3 off); last pieces 13.6, 9.8, 7.2, 7.6, 7.6 px; Tyr into TL, the others into L; no turn over 40 deg. Leu/Ser/His: +13.6 / 0 / -14.2 deg (right side) and the left cluster rows 3.05 px; Ile/Gln/Met -13.4 / 0 / +6.7. Leader length of the whole map 1,003 px (P4 1,075; J2 931), canvas 939.

## Targets
- check_hard 0 of 29: met. doc_tcli01 targets (Leu/Ser/His, Ile/Gln/Met, Trp..Tyr): met. Length x0.849, canvas <= J2: met. Wrong direction 6: met (at the limit).
- Not met: side clusters off the arc 49 (limit 40) and 1 on a sparse map (plasmid37_open_f14, a run with a pole label, kept by P); loose pairs 139 (limit 33; J2 33); mixed 29 of 268 (limit 13; 0 on sparse and medium maps: epifagus 1, chloroplast 3, n4 3, lambda 0); off a contact point on sparse/medium maps: lambda 2+2, n4 2, chloroplast 5 (limit 0); dense maps not above P2 per map: not checked per map; time max 34 s (limit 30).
- Where rules lost: all of them in the runs that contain a pole label (P's DP, thinned seats, as in P4); the dense maps (plasmid153, t5, tobacco, mg1655) fall back wholesale because their runs reach the pole or need more than 40 px of outward move.

# Round 2 (arcfirst=2; Q2): pole runs and long side runs

Code: same file, `arc_tables`, `arc_chain2`, `_lv`, `arc_geom(..., pole, lv)`, hook in `cluster`. `arcfirst=1` is byte-identical to `work/Q` on doc_tcli01 and p1_f14 (cmp, run). Renders `work/Q2/*__J_Q2.*`, table `logs/table-Q2.md`, pictures `work/Q2/view/` (`work/Q2/make_views.sh`). `scripts/arc_fit.py` now prints `zigzag labels` (all clusters, a label whose end is more than 2 px further or nearer than both neighbours): J2 214, P2 265, P4 180, Q 173, Q2 171.

## Built
- A run that arc_chain (round 1) cannot build is given to `arc_chain2`; a run round 1 builds is left as round 1 (so the three side clusters of doc_tcli01 keep the round 1 turns exactly: checked, Trp to Tyr -36.6 / -21.9 / 0 / 1.7 / 8.3 deg). Without this first step the side runs of tcli01 came out 6 to 12 px farther out and with longer leaders (first Q2 attempt, dropped).
- arc_chain2: no 70 degree limit; rows as the stack variable; per label a table of the arc track (valid = arrival point, `static`, turn <= 40) and the ray track (outward along the own ray, valid = `static`), cached; bounded L1 fit = `_pav_l1` clipped by the monotone envelope of the reach (the arc reach both ways; toward the pole the ray track's reach for rays steeper than 45 deg); steep label displaced toward the pole goes on the ray track (straight), anything else on the arc; infeasible reach or a failed test raises the radius of the labels involved by 3 px (smooth, up to `qrhomax` 120). Long-side contact point of a steep leader: the 0.25 point right of the ring axis, the 0.75 point left of it (the box overhangs the axis least; `_lv`).
- Natural straight seat refused by `static` (wide label on a steep ray: the `back` rule) no longer aborts the run; the arc track decides.

## Result (29 cases, run)
- check_hard 0 failing of 29. arcfirst=1 reproduces work/Q.
- Q -> Q2: straight 946 -> 974; leader length 203,256 -> 206,352 (x0.862 of J2); canvas height sum 29,750 -> 30,050 (J2 31,264); mixed 29 of 268 -> 33 of 258; off a contact point 399 -> 397; wrong direction 6; side clusters off the arc 49 -> 45 of 167; loose pairs 139 -> 133; zigzag 173 -> 171; turns 0.5-8 deg 296 -> 286; seconds sum 249 -> 277 (max 30.4 s, p1_f8).
- Canvas height above J2 on 15 cases (Q: fewer): doc_tcli01 946 (J2 934, Q 939), tcli09 952, epifagus 948 (904 in Q), hmmt atskew x2 922, chloroplast 930, lambda_open 948, mito 889 / 888, pbr322 x3 889 / 903, plasmid37 x2 939, t5_f8 1552. Cause: pole runs (tRNA-Pro/Thr, Ser/Asp in tcli01; the top and bottom runs elsewhere) now stand as a staircase on their rays and grow outward.
- lambda_f14: leader length 3,788 (Q) -> 5,126 (J2 4,284, P4 3,704): worse than J2 and Q (the picture shows long leaders on the left side).
- p1_f14: identical to Q (3 runs built, 5 kept for P). Reasons (all `legal`, not reach): the pole runs (terminase .. type III, tail tube .. SimC at the top right and bottom, Lpa-like .. Pap, 33 labels, TciC .., RepL ..).

## Why the pole staircase fails on p1_f14 (found by looking at the failing tests, not only inferred)
The straight leader of a label higher up the staircase passes through the box of the lower label ("leader-hits-box", clearance 2.4 px): beyond 70 degrees from the horizontal only the long-side points are allowed, so a box hangs at least 25 percent of its width over the leader's side; a neighbour's leader is closer than that when the label is wider than about (ray spacing - 2.4) / 0.25 (about 60 px at 14 degrees of ray spacing). Most p1_f14 and lambda pole labels are wider. Rays 45 to 70 degrees can use the vertex (box extends away) and stay clear. So under the 12-point rule a regular straight staircase exists only for short labels; for wide labels the leaders must bend to a vertex (turn on the arc track), which the construction does not offer near the pole because the arc is flat there. A second cause: `back` refuses a wide box hanging on a long-side point because the flat bottom edge comes nearer the ring centre than the arc by more than half a height.

## Targets
Met: check_hard 0 of 29; arcfirst=1 reproduction; doc_tcli01 side clusters as in Q (turns identical); tcli01 pole runs built by Q2 (6 built, 0 kept); wrong direction 6; leader length x0.862.
Not met: p1_f14, p1_f8, lambda_f14, n4_f14 keep 5 / 4 / 2 / 2 runs for P (reasons above); epifagus 1 and chloroplast 1 kept; loose pairs per map not above J2 not checked per map (133 over 29 vs 33); mixed 33 (limit 13; epifagus 1, lambda 1, chloroplast 4, n4 3 are not 0); off a contact point per map vs P2 not checked; zigzag per map not above J2 not checked per map (sum 171 vs J2 214); canvas not above J2 on 15 cases (limit 0); side clusters off the arc 45 (limit 40); time max 30.4 s (limit 30).
Q2 is not better than Q round 1 on the 29 cases: it builds more runs (pole runs on the small maps) but the canvas and lambda leader length get worse and p1_f14 does not change.

# Round 3 (arcfirst=3; Q3): widened vertex rule, rhomax 40

Code: `port_end(..., wide)`, `arc_geom(..., wide)`, `arc_static(..., wide)` (only Q's own calls pass `wide`; P's seats and DP use the old rule), options `qvert`, `qdbg`. `arcfirst=3` = round 2's `arc_chain2` for runs round 1 cannot build, with the vertex allowed from 45 to 90 degrees and `rhomax` 40 (round 2's 120 is not used). `arcfirst=1` still reproduces `work/Q` (cmp, run, doc_tcli01 and p1_f14, checked again after the round 3 edits).
Renders `work/Q3` (strict arrival rule), `work/Q3r` (the relaxed variant, below), table `logs/table-Q3.md`, pictures `work/Q3/view/` (`work/Q3/make_views.sh`), failing tests `work/Q3/p1_f14_failing_tests.txt`.

## Two readings of the vertex (both run)
1. `qvert=0` (default, `work/Q3`): the vertex on the leader's side, as the brief and `port_check.py` say (a leader leaning toward +q ends at the near-end vertex, one leaning away at the far-end one). port_check wrong direction 6.
2. `qvert=1` (`work/Q3r`): always the near-end vertex (the box extends toward +q whichever way the leader leans: the leader then passes below the box). It lets the top-right run of p1_f14 pass the first iteration for its pole part (the leaders no longer cross the neighbours' boxes), but `port_check.py` counts these leaders as "at a point its direction may not use": wrong direction 25 (limit 6). Kept only as evidence; the checker's rule stands.

## Numbers (29 cases, run; Q3 = strict)
| | J2 | Q | Q3 |
|---|---|---|---|
| straight | 1,187 | 946 | 1,010 |
| leader length px | 239,375 | 203,256 | 203,079 |
| canvas height sum px | 31,264 | 29,750 | 29,762 |
| mixed | 50 of 115 | 29 of 268 | 28 of 232 |
| off a contact point | 1,971 | 399 | 392 |
| wrong direction | 118 | 6 | 6 |
| side clusters off the arc | 68 of 131 | 49 of 168 | 48 of 167 |
| loose pairs | 33 | 139 | 129 |
| zigzag labels | 214 | 173 | 169 |
| turns 0.5-8 deg | 332 | 296 | 269 |
| seconds sum (max) | 156 | 249 | 223 (28.4 s) |
check_hard 0 of 29. Runs built / kept per case are in `work/Q3/quick.txt` (columns qb / qk). 
- doc_tcli01: side clusters identical to Q (Trp to Tyr -36.6 / -21.9 / 0 / 1.7 / 8.3 deg); pole runs built by Q3 (6 built, 0 kept); canvas 935 (Q 939, J2 934).
- lambda_f14: 3,788 px (equal to Q). 
- Canvas above J2: doc_tcli01 935 (934), chloroplast 910 (899), pbr322_f10 875 (868), t5_f8 1552 (1517; Q the same). Above Q on 10 cases: hmmt basic x2 954 (929), hmmt atskew x2 913 (908), lambda_f14 / open 930 (927), mito_f14 863 (860), mito_open 879 (872), n4 1004 (1003), pbr322_f10 875 (859): the small pole runs now stand as staircases and add 3 to 25 px.

## p1_f14: the staircase still fails `legal` (stop point of the brief)
p1_f14 is identical to Q (3 runs built, 5 kept for P). Failing test and pair of each kept run (iteration 1 of its radius loop, `qdbg=1`; the final iteration at +40 px has the same kinds):
- terminase .. type III (12 labels): terminase leader through the box of terminase (the second one) and of C1.100; C1.100's leader through Coi, ImcB, ImcA boxes; Coi through ImcB, ImcA; ImcB through ImcA (the run's labels lean to the left across the axis: a vertex on the leader's side makes the box extend away from the next leader's side, so each leader runs through the boxes of the labels below it); with the relaxed vertex (Q3r) this part passes and the failures move to box gaps of site-specific integrase / ImcA / Ref-like RecA / Cra / C8 / Mat (the junction of ray track and arc track: boxes closer than 3.05).
- tail fiber .. C4 (13): SimC leader through PmgC-like's box; SimB's through host cell / SimA / C4; SimA's through SimB's and host cell's; C4's through SimB and SimA; box gap C4 / anti-repressor / host cell division inhibitor (the run meets the neighbouring run).
- Lpa-like .. tRNA-Ile (33): DNA polymerase III theta's leader through UpfO and Lxr boxes; UpfO through UpfN; UpfN through UpfM; Lxr through DNA polymerase's box; box gaps among metallophosphoesterase / morphogenesis / alpha-2;3-sialyltransferase / PmgN / VRR-NUC / PmgL.
- TciC .. ParB-like (18): box gap and leader through the box of tRNA-Ile, tellurite resistance, Upl, Plp, tRNA-Thr, tRNA-Asn: the run is not a pole run but a compressed side run that also touches the tRNA-Ile / Thr / Asn labels outside it.
- tail sheath .. host cell division inhibitor (12): RepL leader through DUF3800; DUF3800 through RepL; the three anti-repressors: leader through each other's box.
So the cause moved: the arrival rule is no longer the blocker for the top-right pole labels (relaxed variant), but (a) a leader leaning away from the corner is counted "wrong direction" by the checker, (b) the junction of ray and arc track gives box gaps, (c) the 18 and 33 label runs need more than 40 px of outward move (rhomax 40).

## Targets
Met: check_hard 0 of 29; doc_tcli01 side clusters identical and canvas <= 939; lambda_f14 leader length equal to round 1; 29 cases not worse than round 1 on off-arc (48), loose (129), mixed (28), zigzag (169), off a contact point (392), wrong direction (6), leader length (203,079), seconds (28.4 s max).
Not met: p1_f14 pole runs built (5 runs kept for P, tests above); canvas not above round 1 per case (10 cases above); canvas above J2 on 4 cases.
