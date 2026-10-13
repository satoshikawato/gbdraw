# Design K: leaders may pass behind labels (revision of H2 / design I)

Date 2026-10-12. Code: `gbdraw/labels/proto/design_k.py` (worktree `.worktrees/label-study-k`, uncommitted, ruff clean,
deterministic) and one additive scaffold edit in `gbdraw/render/groups/circular/labels.py` (+69/-1).
Renders: `work/K/` (`<case>__K.*` for all 29 cases, completed with `lp_objects.py`). Own table script: `work/K/tabk.py`.

## Result in short

- K meets every hard criterion on all 50 K renders (29 defaults, settings, looks): dropped 0, crossings 0, text
  overlap 0, toward axis 0, gap 3.0, on/near objects 0, leaders on objects 0, acute 0, bends < 110 deg 0.
  Only `leader_through_text` is above 0, by design.
- The big gain Owner hoped for is NOT there on P1. `p1_f14` goes 1230 -> 1165 (default), 1211 (mild), 1139 (strong).
  The towers shrink a little; they do not disappear. Where it helps more: `t5_f8` 1678 -> 1490 / 1265 (default / strong),
  and total leader length falls a lot (`t5_f8` 36803 -> 30319 -> 23057).
- Cases that fit under H2 come out identical (same SVG geometry, no box): 21 of 29, including `mito_f14`, `pbr322_f14`,
  `gal_hmmtdna_basic_circular`, `doc_tcli01`, `doc_tcli03`.

## Why the gain is bounded (planarity)

Leaders may not cross, so the labels of one stack keep the feature order. A dense cluster of features therefore sends its
leaders to labels in the same order, and a label can only move into a hole beside the stack when its leader can reach it
without crossing a neighbour. Only the last labels of a stack (nearest the pole) can use holes. The saving per tower is
limited by how many leaders may pass behind (`maxb`) and by the hole width, and it saturates near 100 px (8%) on P1. That
confirms the brief's estimate: a better search or covering boxes cannot remove the towers; only a different arrangement
(rows with several labels, or a reorderable leader scheme) could.

## The rule in plain steps

1. Run H2 per shore, column and half (stacks from the equator to the pole). With `over=0` or when no stack gains, the
   result is exactly H2.
2. A plan says for each stack: the last `n` labels are "over" labels. A ladder search (n = 0, 2, 4, ...) tries each plan,
   re-runs the placement and scores it: leader length sum + hidden penalty + `wy` * (stack extent beyond `cap`).
   A plan is kept only if it scores lower than n = 0. That is the "use it only where it pays" decision, per stack.
3. An over label scans rows from its own height up to `cap` in steps of `step` px. It may take a spot where its leader
   (or an earlier leader through its box) passes behind another box, if the limits below hold. Extra candidates are tried
   beyond a box's far end and on its pole side. Rays reserved for later labels are soft for over labels.
4. Cost of a spot = leader run + `pen` per box passed + `hw` per hidden px. The cheapest valid spot wins.
5. Everything else of H2 stays hard: gap 3 px between texts, box clear of every other drawn object, no leader
   crossing (also behind boxes), nothing dropped, LP-D09, one elbow >= 110 deg, leaders clear of other objects.

## When a label may cover a leader, and the limits

- A leader passes behind at most `maxb` boxes (default 3) and at most `maxc` px (50) in any one box.
- At most `maxh` (0.7) of its length is hidden.
- Its elbow is never hidden, and `vis` px (14) of the last piece before its own label stay visible; `vis0` px (6) at the
  feature end. A leader never runs along inside a box (the chord limit and a no-parallel test).
- Mild / default / strong, all deterministic:

| setting | maxb | maxc | maxh | vis / vis0 | cap | pen | hw | patience / trials | time on t5_f8 |
|---|---|---|---|---|---|---|---|---|---|
| mild | 1 | 25 | 0.35 | 14 / 6 | 200 | 25 | 0.3 | 3 / 10 | about 10 s |
| default | 3 | 50 | 0.7 | 14 / 6 | 150 | 25 | 0.3 | 3 / 10 | 21 s |
| strong | 6 | 80 | 0.95 | 4 / 4 | 100 | 10 | 0.1 | 20 / 14 | 56 s |

Other options: `over` (1), `wp` 2, `wy` 20, `step` 6. Rendering only (placement unchanged): `bg` box opacity (default 0.8,
boxes only on labels that a leader passes behind), `bgpad` 1.5 px, `ul`, `tcol`, `gapped` (interrupt the leader instead
of drawing a box). Draw order: leaders, boxes, underlines, text. Default time is about 10 s on P1 and 21 s on `t5_f8`
(brief: about 15 s; `t5_f8` and `plasmid153_f14` at 19 s exceed it).

## K against I (default)

Behind = leaders passing behind at least one box / most boxes for one leader / hidden px (total) / largest hidden share
of one leader. Canvas width is unchanged everywhere. Turn = `leader_turn_deg_mean`.

| case | n | canvas height I -> K | sumL I -> K | maxL I -> K | above/below ring I -> K | behind | turn I -> K |
|---|---|---|---|---|---|---|---|
| p1_f14 | 111 | 1230 -> 1165 | 8975 -> 7939 | 260 -> 221 | 23/16 -> 19/16 | 10 / 2 / 174 / 0.39 | 35.1 -> 36.2 |
| tobacco_f14 | 140 | 1431 -> 1363 | 17411 -> 15944 | 403 -> 362 | 23/36 -> 23/36 | 30 / 2 / 697 / 0.24 | 37.1 -> 38.2 |
| t5_f14 | 195 | 2630 -> 2329 | 65511 -> 61710 | 1030 -> 819 | 53/62 -> 52/61 | 45 / 2 / 1027 / 0.16 | 47.0 -> 47.4 |
| t5_f8 | 194 | 1678 -> 1490 | 36803 -> 30319 | 471 -> 398 | 37/46 -> 32/44 | 50 / 3 / 1090 / 0.19 | 44.5 -> 45.5 |
| plasmid153_f14 | 179 | 1863 -> 1755 | 42504 -> 36901 | 530 -> 481 | 49/52 -> 46/52 | 42 / 3 / 1757 / 0.17 | 44.5 -> 45.6 |
| plasmid153_f8 | 179 | 1217 -> 1151 | 14628 -> 13001 | 223 -> 185 | 28/31 -> 26/31 | 24 / 3 / 388 / 0.24 | 41.4 -> 42.5 |
| n4_f14 | 69 | 1046 -> 1020 | 3624 -> 3226 | 195 -> 168 | 17/3 -> 16/3 | 2 / 3 / 67 / 0.36 | 27.2 -> 27.9 |
| tobacco_f8 | 140 | 1071 -> 1045 | 7922 -> 7490 | 187 -> 187 | 12/26 -> 12/26 | 11 / 2 / 129 / 0.11 | 29.6 -> 30.0 |

The other 21 cases are identical to I (no box drawn, no leader behind a box): `pbr322_f14`, `mito_f14`, `mg1655_rna_f8`,
`lambda_f14`, `epifagus_f14`, `plasmid37_f14`, `p1_f8`, `pbr322_f10`, `gal_hmmt_basic`, `gal_hmmt_atskew`,
`pbr322_open_f14`, `lambda_open_f14`, `p1_open_f14`, `mito_open_f14`, `plasmid37_open_f14`, `gal_hmmtdna_basic_circular`,
`gal_hmmtdna_atskew`, `gal_tobacco_chloroplast_hz`, `doc_tcli01`, `doc_tcli09`, `doc_tcli03`. Full 29-row table:
`/home/kawato/gbdraw-baselines/tmp/k/table_main.md` (the SVGs of `p1_f14:I`, `p1_open_f14:I`, `gal_tobacco_chloroplast_hz:I`
are byte-identical before and after the scaffold edit: `work/K/scaffold_identity.txt`).

Leader length gets shorter everywhere K acts. The mean turn angle rises by 0.4 to 1.5 deg (slightly worse than I; LP-D11
point 1 asked for fewer turns, so K trades a little of it for length).

## Hard counts

Over 50 renders (`*.metrics.json` after `lp_objects.py`): every hard count above is 0, `label_gap_min_px` 3.0.
`leader_through_text` (metric) is above 0 where leaders pass behind boxes: `p1_f14` 11, `t5_f8` 109, and similar on the
other eight changed cases. This is expected.

## Trade-off settings (canvas height / sumL / leaders behind / most / hidden px / max share)

| case | I | mild | default | strong |
|---|---|---|---|---|
| p1_f14 | 1230 / 8975 | 1211 / 8466 / 6 / 1 / 87 / 0.20 | 1165 / 7939 / 10 / 2 / 174 / 0.39 | 1139 / 7266 / 11 / 3 / 245 / 0.51 |
| t5_f8 | 1678 / 36803 | 1652 / 35968 / 10 / 1 / 99 / 0.03 | 1490 / 30319 / 50 / 3 / 1090 / 0.19 | 1265 / 23057 / 66 / 6 / 2179 / 0.39 |
| t5_f14 | 2630 / 65511 | 2349 / 62445 / 43 / 1 / 749 / 0.16 | 2329 / 61710 / 45 / 2 / 1027 / 0.16 | 2100 / 51160 / 81 / 6 / 4883 / 0.31 |
| plasmid153_f14 | 1863 / 42504 | 1822 / 40672 / 43 / 1 / 740 / 0.06 | 1755 / 36901 / 42 / 3 / 1757 / 0.17 | 1716 / 32544 / 66 / 6 / 3168 / 0.40 |
| tobacco_f14 | 1431 / 17411 | 1392 / 16436 / 30 / 1 / 510 / 0.13 | 1363 / 15944 / 30 / 2 / 697 / 0.24 | 1324 / 14900 / 35 / 3 / 1270 / 0.35 |

Each step costs: more leaders passing behind (6 -> 10 -> 11 on P1), longer hidden stretches, and slightly higher mean turn.
The settings are not strictly monotone (the search is greedy; `t5_f14` mild is already near default). Width does not change.
Strong takes 40 to 56 s on the T5 cases.

## Looks on `p1_f14` (suffixes in `work/K/p1_f14__K<suffix>.svg/png`)

`_nobg` (no box), default (bg 0.8), `_bg06`, `_bg10`, `_gap`, `_ul` (underline plus box), `_ulnb` (underline, no box),
`_tcol`, `_tcol_ul`; `t5_f8` has `_ul`, `_gap`, `_bg10`. Crops at 1:1: `work/K/view/looks_topleft.png`,
`looks_bottom.png`, `looks_t5f8.png`, `p1_K_top.png`, `p1_K_bot.png`.

## Assessment from the 1:1 figures

- Can one follow every leader on `p1_f14`? Yes. Passages are short (13 to 23 px chords) and only 10 leaders of 111 pass
  behind a box, at most two boxes each. The short labels (PdcA, Lxr, UpfM, PmgU, anti-repressor) sit in holes beside the
  stack; their leaders disappear behind a neighbour for a few pixels and appear again. Wherever a leader is hidden 39% of
  its length the remainder still points clearly at the feature. On `t5_f8` default it is harder: 50 leaders, three boxes
  deep in places; in the dense hypothetical-protein bundles one cannot be sure which leader enters a box (strong: worse).
- Which look reads best: `ul=1` with the box. The underline ties each label to the end of its own leader, so a leader
  that vanishes behind a box and a different leader that touches a text are told apart. Next `gapped=1`: clean, no box
  edge at all, but it cuts leaders into pieces, which looks busy on `t5_f8`. Plain `bg` 0.6 to 0.8 works but the effect
  is hard to see on P1 (few boxes); `bg` 1.0 hides the passing leader completely and makes the leader seem to stop at the
  box: avoid. `ul` without a box is readable too but lines run into neighbouring text on dense stacks.
- Which look reads worst: `tcol`. The CDS colour is light blue; on white it has low contrast, which supports the
  Owner's doubt. Not recommended with default palettes.
- Where it fails: tall labels and the T5 right-hand tower stay as under H2 (T5 keeps most labels in one column; the cut
  problem at the poles is unchanged, `t5_f14` remains 2100 to 2329 high). The dense P1 staircases and tower bundles
  still need the vertical room; K saves at most about 8% of height. A case that did not change (for example
  `mg1655_rna_f8`) looks exactly as in H2, as intended.

## What production would need

- SVG z-order: leaders group before the feature tracks (already so), then the label boxes, underlines and text in the
  text group, in that order. K needs the boxes to be separate elements bound to their label (`bind_label_part`) so that
  hiding or moving a label moves its box.
- Web editor: dragging a label must re-run the over rule (or at least re-test limits) because a box can start to cover
  a leader or stop covering it. I did not examine the editor's z-order or label dragging; nothing was run in the Web.
- Search: the ladder search is deterministic but takes 10 to 56 s; production needs a time budget, a cap on trials, or a
  cheaper plan (for example only the last 8 labels of a stack).
- Options to decide: `ul` and `bg` defaults, the three limit presets, and whether to keep `gapped`.

## What could not be done

- No large tower removal on P1 (only 5.3% default, 7.4% strong, against about 1,050 px ideal). Planarity and hole width
  bound it; a different arrangement (rows with several labels) was not built.
- `t5_f14` stays at 2100 to 2329 px high (T5 cut problem).
- The search is greedy; settings are not strictly monotone.
- Time limit of about 15 s exceeded on `t5_f8` (21 s) and `plasmid153_f14` (19 s); strong
  settings up to 56 s.
- The Web editor and its z-order were not examined or tested.
- Mean leader turn is slightly worse than I (+0.4 to +1.5 deg).
