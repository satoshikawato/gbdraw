# Design L: blocks of rows at 12 and 6 o'clock, leader ends anywhere on the label box

2026-10-12. Study prototype, not for merging. Code: `gbdraw/labels/proto/design_l.py` in `.worktrees/label-study-k` (uncommitted, ruff clean;
imports the helpers of designs I, G and K; no renderer change, the K rendering options still apply). Renders: `work/L/` (29 cases
`<case>__L.*`; `__L_f0o0 .. _f1o1` = the four option combinations on 7 cases; `_p0`, `_nb`, `_p0nb` = rule comparisons on 4 cases), all completed with
`lp_objects.py`. Own scripts: `work/L/tabl.py` (L against I and K, 29 rows), `tabm.py` (combinations), `mkview.py`; crops: `work/L/view/`.

## Result in short

- On `p1_f14` the canvas goes 1230 (H2) -> 1165 (K) -> **1047** (L default), 1063 / 1114 / 1148 for the other `over`/`flip` combinations; labels wholly
  above/below the ring 23/16 -> 18/14; rows of the blocks 5 (top) and 4 (bottom) instead of towers of up to 14; total leader length 8,975 -> 6,626 px;
  longest leader 260 -> 136 px; straight leaders 18 -> 46 of 111; mean turn 35 -> 20 degrees.
- All hard criteria hold on all 29 default renders and on all option renders: nothing dropped, 0 crossings, 0 text overlaps, gap 3.0 px, 0 labels on or
  within 3 px of any object, 0 leaders on objects, 0 bends under 110 degrees. Not hard any more (LP-D13): `labels_toward_axis` (10 on `p1_f14`).
- Not achieved: the towers do not vanish (a block still needs 4 to 5 rows, about 120 px above the ring, because the leaders may bend by at most 70
  degrees and may not run far sideways); T5 and plasmid153 at 14 pt are not improved (see table); run time on the large maps is above 15 s.

## The rule in plain steps

1. Main layout = design I (stacks per shore, cut at the poles, reserved exits) with the new edge rule, everywhere on the map: `perim` = 1: a leader may
   end at any point of the edge of its box, `lgap` (2.5 px) from the text. Candidates: the nearest point of the box (shortest leader) and the first
   point where the feature's own ray meets the box (a straight leader ends where its line first meets the box: LP-D14; e.g. "PmgG-like head
   morphogenesis" on `p1_f14` meets the top edge, turn 0), for every box position of design I and for boxes that cover the feature (centred above or
   below it). Cost of a spot = leader length + `wt` (0.35) px per degree of turn. `perim` = 0 restores H2 exactly (byte-for-byte same metrics as
   `work/I`, `*_L_p0nb`).
2. LP-D09 as its intent: no point of the box may be more than `backf` x text height (0.5) nearer to the ring than the leader end (`back`, measured as
   `scripts/lp_metrics.py` `back_px`). `perim` = 0 keeps H2's corner test.
3. Blocks (`blocks` = 1): per pole the labels whose box lies wholly above (below) the ring tangent (set `kap` = 1) or beyond 0.7 R (set 0.7) are laid
   out again, outer features first, one at a time: any height from the feature's own upward, any box position that a leader of at most `maxrun` (140) px
   sideways and a turn of at most `lean` (70 degrees, so every bend is 110 degrees or more) reaches. Cost = leader + turn + `wh` (1.2) per px of rise
   + (`over`) 25 px per box passed + 0.3 px per hidden px. A label that finds no spot is seated first in the next attempt (6 attempts).
4. Until it is seated, every label keeps a free way out (a ray straight along the pole direction, inside the cone): a leader may not cross it, a box may
   not cover its first 6+2.4 px (with `over` = 0: not cover it at all). This is what stops the greedy layout from trapping a feature.
5. `over` = 1: a leader may pass behind another label's box under K's limits (3 boxes, 50 px in one, 70 % hidden, elbow and the last 14 px visible);
   `over` = 0: leaders keep `lpad` from every box. `flip` = 1: in a block the text may run towards the pole axis (back up to `flipb` = 45 px); `flip` = 0:
   `back` stays at `backf`. Then every label is seated again, highest first, with all others in place (2 sweeps), and kept if it comes lower.
6. A block replaces the H2 layout of a pole only if the outer edge of the whole outside layout gets lower (best of the two sets, then total length);
   otherwise H2's layout stands (checked for overlap, crossing, gap). Deterministic. Leaders that pass behind boxes get K's filled box (`bg` 0.8).

## Options

`perim` 1 | `backf` 0.5 | `wt` 0.35 | `blocks` 1 | `over` 1 | `flip` 1 | `flipb` 45 | `lean` 70 | `maxrun` 140 | `wh` 1.2 | `maxb/maxc/maxh/vis/vis0`
3/50/0.7/14/6 | `pen` 25 | `hw` 0.3 | `cap` 520 | `step` 3 | `gstep` 3 | `sweeps` 2 | `repair` 6 | `rays` 1 | `rayw` 0 | `ord` 0 | `kap` 0 (= try 1 and 0.7) |
`tmax` 6 s. Rendering options of K (`bg`, `bgpad`, `ul`, `tcol`, `gapped`) work unchanged.

## L against I and K (default: perim 1, over 1, flip 1; full 29 rows: `tabl.py`)

canvas h I / K / L | wholly above/below I / K / L | rows top/bot L | sum leader px I / K / L | longest I / L | turn mean I / L | straight I / L
| case | canvas I/K/L | above/below I/K/L | rows | sumL I/K/L | maxL I/L | turn I/L | straight I/L |
|---|---|---|---|---|---|---|---|
| p1_f14 | 1230/1165/**1047** | 23/16, 19/16, 18/14 | 5/4 | 8975/7939/6626 | 260/136 | 35.1/20.1 | 18/46 |
| n4_f14 | 1046/1020/**957** | 17/3, 16/3, 17/2 | 4/2 | 3624/3226/2608 | 195/104 | 27.2/11.6 | 20/43 |
| lambda_f14 | 960/960/**924** | 9/3, 9/3, 8/2 | 3/1 | 4091/4091/3936 | 148/118 | 28.9/11.9 | 15/44 |
| epifagus_f14 | 930/930/**897** | 9/1, 9/1, 8/1 | 2/1 | 2497/2497/2384 | 133/105 | 25.6/17.5 | 20/31 |
| p1_f8 | 966/966/**919** | 15/11, 15/11, 12/11 | 3/4 | 4408/4408/4081 | 105/84 | 28.2/14.4 | 37/72 |
| plasmid153_f8 | 1217/1151/**1096** (width 1314 -> 1341) | 28/31, 26/31, 22/29 | 7/8 | 14628/13001/11816 | 223/174 | 41.4/25.4 | 13/44 |
| t5_f8 | 1678/**1490**/1511 | 37/46, 32/44, 37/46 | 16/16 | 36803/30319/33588 | 471/406 | 44.5/41.7 | 16/30 |
| t5_f14 | 2630/**2329**/2350 | 53/62, 52/61, 54/62 | 27/27 | 65511/61710/63380 | 1030/817 | 47.0/41.0 | 7/15 |
| plasmid153_f14 | 1863/**1755**/1870 (width 1608 -> 1642) | 49/52, 46/52, 49/52 | 19/17 | 42504/36901/42681 | 530/534 | 44.5/39.1 | 8/10 |

The other 20 cases (mito, pbr322, gal, doc tutorials, chloroplast ...): canvas 0 to 71 px lower than H2 (`mg1655_rna_f8` -71, `p1_open_f14` -41), leader length 1 to 10 % lower (chloroplast
-10 %), straight leaders about +50 %, no block. `doc_tcli03` identical to H2.

## The options on the seven cases (canvas h; above/below; rows; sum leader px; turn mean; straight)

| case | H2 | K | flip0 over0 | flip0 over1 | flip1 over0 | flip1 over1 |
|---|---|---|---|---|---|---|
| p1_f14 | 1230; 23/16; -; 8975; 35; 18 | 1165; 19/16; 7939 | 1148; 19/14; 7/5; 7610; 19; 49 | 1063; 17/14; 6/4; 6701; 21; 47 | 1114; 19/14; 6/5; 7257; 19; 49 | **1047**; 18/14; 5/4; 6626; 20; 46 |
| n4_f14 | 1046; 17/3 | 1020 | 996; 16/2; 5/2; 2845 | 957; 17/2; 4/2; 2607 | 996 (as flip0) | 957; 17/2; 4/2; 2608 |
| lambda_f14 | 960; 9/3 | 960 | 944; 9/2; 4/1 | 933; 10/2; 3/1 | 944 | 924; 8/2; 3/1 |
| epifagus_f14 | 930; 9/1 | 930 | 918; 9/1; 4/1 | 902; 9/1; 3/1 | 918 | 897; 8/1; 2/1 |
| p1_f8 | 966; 15/11 | 966 | 939; 13/11; 4/4 | 927; 12/11; 3/4 | 925; 13/11; 3/4 | 919; 12/11; 3/4 |
| plasmid153_f8 | 1217; 28/31 | 1151 | 1228 (no block) | 1107; 19/29; 8/6 | 1228 (no block) | 1096; 22/29; 7/8 |
| t5_f8 | 1678; 37/46 | 1490 | 1511 (no block) | 1365; 31/37; 13/12; 26554 | 1511 (no block) | 1511 (no block, too slow) |

`over` is what matters: without it a leader may not pass a box, the block cannot beat H2 by much (p1_f14 1148/1114) and on dense maps no block is
accepted. `flip` adds 16 to 50 px (p1_f14, p1_f8). Neither: 1148. With `over` and `flip` 1047. Leaders behind boxes at the default (`p1_f14`): 20 leaders,
at most 3 boxes for one, 517 px hidden in all (K: 10 leaders, 174 px). `flip` use at the default: 28 of 111 labels come nearer to the ring than their
leader end by more than 1 px, 29.3 px at the most (`flip` 0: 24 labels, 8.6 px); `toward_axis` 10 (flip 0: 2, perimeter rule alone 4).
Sideways run of a leader: 29.9 px at the most (`maxrun` 140 is not reached; flip 0 over 1: 52 px).

## What the edge rule (LP-D13) gains against the corner rule

Corner rule = `perim` 0 (H2's leader at the text end); in the block rows `flip` 1 and `over` 1 as in the default.

| case | rule | canvas | above/below | sum leader px | straight | turn mean | toward_axis |
|---|---|---|---|---|---|---|---|
| p1_f14 | corner, no blocks (H2) | 1230 | 23/16 | 8975 | 18 | 35.1 | 0 |
| | corner + blocks | 1052 | 19/14 | 6918 | 25 | 31.0 | 16 |
| | edge, no blocks | 1155 | 20/14 | 7874 | 49 | 19.2 | 4 |
| | edge + blocks | **1047** | 18/14 | 6626 | 46 | 20.1 | 10 |
| n4_f14 | H2 / corner + blocks / edge, no blocks / edge + blocks | 1046 / 959 / 1031 / **957** | 17/3 / 17/2 / 16/2 / 17/2 | 3624 / 2840 / 3070 / 2608 | 20 / 25 / 41 / 43 | 27 / 19 / 12 / 12 | 0 / 7 / 4 / 10 |
| lambda_f14 | same order | 960 / 920 / 944 / 924 | 9/3 / 10/3 / 9/2 / 8/2 | 4091 / 4083 / 4067 / 3936 | 15 / 20 / 44 / 44 | 29 / 22 / 12 / 12 | 0 / 4 / 3 / 4 |
| plasmid153_f8 | same order | 1217 / 1148 / 1228 / **1096** | 28/31 / 24/31 / 28/31 / 22/29 | 14628 / 13389 / 14831 / 11816 | 13 / 19 / 36 / 45 | 41 / 40 / 25 / 25 | 0 / 10 / 1 / 19 |

The edge rule alone lowers the canvas a little (p1_f14 -75 px) and straightens leaders a lot (straight 18 -> 49, mean turn 35 -> 19). Together with
blocks the canvas height is the same as with the corner rule (1047 against 1052, 957 against 959; plasmid153_f8 1096 against 1148): the height is
set by the block geometry, not by the leader end. The gain of the rule is in leader length (-4 to -8 % against corner + blocks) and turns. In all
four cases the edge rule is not worse than the corner rule on any hard criterion. Leaders that end at a corner of the box grown by `lgap`: 57 of 111
on `p1_f14` (4 of them straight); this is the nearest-point rule when the feature lies diagonally to the box. Not changed (open point, LP-D14).

## Assessment from 1:1 crops (`work/L/view/p1_{top,bottom}_*.png`, next to H2, G and K)

- Top of `p1_f14`, H2 (1230): one tower of 14 labels on the left, one of 7 on the right. G (891 high): rows hugging the ring, but 26 labels dropped,
  see the brief. K: the same towers, a few labels in holes.
- L flip0 over0 (1148): leader lines are all clean and followable, but the left tower still has 10 labels; reads like H2 with a few rows
  merged. L flip1 over0 (1114): the pole rows (terminase / C1.100 / Coi / ImcB) now hold 3 to 4 labels; the left tower has 8.
- L flip0 over1 (1063) and default (1047): the left side becomes 4 to 5 rows of two to four labels, the pole gets a block of 3 rows. Can one follow
  every leader? Mostly yes: a leader that passes behind a box is faded by the box (opacity 0.8) and comes out on its far side in line; where three
  leaders go behind one text ("Lpa-like transcriptional activator", "UmuD-like protein") the lines inside the text are visible and one has to read
  the stub that comes out above. In the default, several neighbours in one row sit 3 px apart ("UmuD-like protein terminase", "C1.100 Coi")
  and read like one phrase, as in G; an underline (`ul=1` of K) or a wider gap would be needed. The block reads as belonging to its features: the
  labels stay above the stretch of ring they describe, leaders fan out by at most 30 px sideways (the Owner's horizontal-leader complaint does not apply).
- Bottom of `p1_f14`: 14 labels in 4 rows instead of a 16-label tower: "anti-repressor / SimA / SimC / C4 / SimB" side by side, leaders
  cross behind "anti-repressor" and "host cell division inhibitor". `over` 0 gives 5 rows and still reads clean.
- Best reading: `flip1 over1` (lowest, most compact); the most conservative useful setting is `flip1 over0` (nothing hidden, 1114). Where it fails:
  dense one-column maps (T5, plasmid153 at 14 pt: the height is made by the side columns, not the poles); a cluster of 13 features within 100 px
  needs 5 rows because each row can hold only about two boxes per narrow cluster and leaders may not lean more than 70 degrees.

## What production would need

- The edge rule is a change of the leader end in the candidate generator and one new test (`back`); leaders end on a box edge, which the renderer
  already draws. Blocks need K's z-order (leaders, filled boxes, text) and re-running the rule when the Web editor moves a label.
- Time: 1.1 to 12 s on 24 of 29 cases; 14 to 28 s on `plasmid153_f8/f14`, `t5_f8/f14`, `tobacco_f14` (four jobs in parallel; H2 alone needs 2 to 3 s,
  K 10 to 21 s). The block search is numpy per row and candidate; it needs a budget (`tmax` 6 s stops new variants; a block is abandoned after 12 s per pole).
- The greedy order can trap a feature (all earlier failures): the reserved rays and the repair attempts make it work on every case here, but it is not
  proven; a block that fails falls back to H2 for that pole.

## What was not done

- No tower removal like G's (G reaches 891 px because it drops 26 labels); L stays at about 120 px of rows above the ring.
- T5 14/8 pt and plasmid153 14 pt: no gain over H2 (T5 8 pt: no block was accepted with `flip` 1 `over` 1 within the time limit; `flip0 over1` got 1365).
- The corner ends (LP-D14) of nearest-point leaders: not changed. `labels_toward_axis` is no longer a criterion; it is reported (10 on `p1_f14`).
- Run time above 15 s on the large maps. Rows above the ring were not checked in the Web editor.
