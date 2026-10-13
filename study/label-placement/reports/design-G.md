# Design G: liberal placement on both sides of the ring (Geneious look)

File `gbdraw/labels/proto/design_g.py` (about 480 lines, deterministic, numpy). Renders `work/G/<case>__G.*`; variants `__G_in0`, `__G_merge`,
`__G_color`, `__G_clear6` (p1_open, lambda_open, pbr322_open); table `work/G/summary_G.txt`; independent checker `work/G/check_g.py`.
Scaffold edits (additive; p1_f14:D and mito_f14:E SVGs are byte-identical before and after): hook in `labels/circular.py` passes
`radial_layout=` and `feature_colors=` to `run_design`; `proto/__init__.py::run_design` takes both and sets `ctx.radial_layout`, `ctx.feature_colors`;
`Ctx` has two defaulted fields (`common.py`); `render/drawers/circular/labels.py`: 2 lines, `text_fill` sets the label's SVG fill (used only by `color=1`).
(The first drawer edit had normalised mixed CRLF/LF line endings, 198 changed lines; it was reverted and redone byte-wise: now 2 added lines.)

## 1. The rule
- Zones come from `ctx.radial_layout` (real radii, not guesses): feature ring [r_in, r_out] (308..386 px), tick-mark band, every other track's reserved band
  (GC, skew), the title disc (radius `definition_reserved_band`), and the major tick labels (re-computed: auto interval, 14 px font from the label band,
  each label as a chain of small boxes along its tangent). All are "forbidden annuli / boxes"; a label box must keep `clear` px from each.
- Each label has two stubs on its feature's ray: outside (ring outer edge) and inside (ring inner edge). Leader = anchor -> stub (radial, hidden under the
  feature) -> straight line to the nearest point of the box (attach slid 4 px off corners). Text horizontal, box = text size.
- Candidates: box centres on a 3 px grid within a leader length L of the stub, L = 30, 60, 120, 240, 1.3 R (the first level with a feasible spot wins).
  Cost = leader length (stub to box) + 0.1 R x angular drift + 0.5 x length x tilt from radial + 4 px for inside. Hard limits: leader at most 75 degrees
  from the radial (85 at the last level); outside boxes entirely beyond r_out + clear, inside boxes entirely below the tick band and above tracks/title.
- Collision: box vs boxes (incl. tick labels) with gap `clear`; new leader vs placed boxes and placed leaders vs new box with pad `clear`; no leader crossing;
  leaders may not enter the ring or a track. Corridors: the radial line (36 px, hard; 50 px, soft 30 px per hit, cap 6) of every not-yet-placed stub
  is kept free, on the side being used, so neighbours keep an escape.
- Order: most crowded stretch first (labels within +-0.12 rad), then 2 refinement passes (each label re-placed at its cheapest spot, no corridors).
- Nothing fits (last resort, in this order): same search without corridors, crossings allowed (60 px per crossing), leader pad 1 px; then drop (`ctx.dropped`).
- Options: `inside=1`, `merge=0` (runs of clockwise-adjacent equal texts within 0.06 rad, span 0.14 rad: one label "text (x2)", one leader per
  feature to a shared attach point; chosen over a single leader because each feature stays identifiable), `color=0` (text in feature colour), `clear=3` (= ctx.gap).

## 2. Parameters
clear 3 | step 3 | lmax 1.3 R | drift 0.1 | tilt 0.5 | inpen 4 | amax 75 | harddepth 36 | raydepth 50 | raypen 30 | raycap 6 | passes 2 | retries 0 (1 to 2 retries
moved a few drops to others, no net gain, 2x time) | crowdw 0.12 | mergegap 0.06 | mergespan 0.14.

## 3. Metrics (default options; "cur" = today: overprint pairs / crossings / leader mean/max / canvas)
| case | req | placed | in | drop | overprint | crossings | thru-text | leader mean/max | label-feature mean | canvas | cur |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pbr322_open_f14 | 34 | 34 | 3 | 0 | 0 | 0 | 0 | 33/141 | 31 | 1145x868 | 0/0/90/174, 1256x910 |
| mito_open_f14 | 33 | 33 | 3 | 0 | 0 | 0 | 0 | 38/83 | 36 | 1320x864 | 0/0/88/127, 1415x906 |
| plasmid37_open_f14 | 45 | 45 | 2 | 0 | 0 | 0 | 0 | 55/81 | 53 | 1490x904 | 0/1/120/248, 1580x916 |
| lambda_open_f14 | 71 | 69 | 9 | 2 | 0 | 1 | 0 | 47/269 | 45 | 1347x879 | 0/2/115/201, 1499x910 |
| p1_open_f14 | 111 | 107 | 26 | 4 | 0 | 14 | 0 | 48/187 | 45 | 1551x943 | 20/21/115/467, 1706x990 |
| pbr322_f14 | 34 | 33 | 2 | 1 | 0 | 0 | 0 | 29/78 | 28 | 1145x861 | 0/0/90/174, 1256x910 |
| mito_f14 | 33 | 33 | 0 | 0 | 0 | 0 | 0 | 41/103 | 38 | 1380x864 | 0/0/88/127, 1476x906 |
| plasmid37_f14 | 45 | 45 | 0 | 0 | 0 | 0 | 0 | 57/110 | 55 | 1558x926 | 0/1/120/248, 1648x916 |
| lambda_f14 | 71 | 68 | 0 | 3 | 0 | 10 | 0 | 62/347 | 59 | 1510x892 | 0/2/115/201, 1567x910 |
| n4_f14 | 69 | 61 | 1 | 8 | 0 | 3 | 0 | 36/206 | 34 | 1278x885 | 0/0/124/230, 1402x994 |
| p1_f14 | 111 | 92 | 4 | 19 | 0 | 15 | 0 | 57/395 | 54 | 1551x909 | 20/21/115/467, 1706x990 |
| p1_f8 | 107 | 103 | 10 | 4 | 0 | 1 | 0 | 36/100 | 33 | 1299x888 | 0/0/102/142, 1441x982 |
| epifagus_f14 | 53 | 45 | 0 | 8 | 0 | 2 | 0 | 44/174 | 41 | 1393x868 | 0/0/111/170, 1509x990 |
| tobacco_f8 | 140 | 118 | 4 | 22 | 0 | 25 | 0 | 46/261 | 43 | 1305x922 | 0/0/122/226, 1399x983 |
| tobacco_f14 | 140 | 97 | 2 | 43 | 0 | 28 | 0 | 68/289 | 63 | 1805x946 | 50/1062/270/705, 1634x988 |
| t5_f8 | 194 | 148 | 1 | 46 | 0 | 178 | 0 | 84/324 | 79 | 1420x919 | 0/2085/347/699, 1444x979 |
| t5_f14 | 195 | 116 | 0 | 79 | 0 | 248 | 0 | 137/482 | 130 | 1944x1019 | 170/2532/334/899, 1821x992 |
| plasmid153_f8 | 179 | 142 | 2 | 37 | 0 | 67 | 0 | 56/354 | 53 | 1430x919 | 0/0/103/207, 1436x983 |
| plasmid153_f14 | 179 | 102 | 0 | 77 | 0 | 48 | 0 | 102/377 | 99 | 1753x918 | 143/17/92/177, 1714x994 |
| mg1655_rna_f8 | 108 | 80 | 6 | 28 | 0 | 42 | 0 | 49/187 | 44 | 1103x840 | 0/0/109/150, 1200x977 |
Overlap pairs 0 everywhere. `leaders_over_features` equals the number of inside labels (counts inside leaders by design). Run time: 1 to 20 s up to P1,
30 s tobacco, 50 to 70 s T5 and plasmid153 (above the 10 s target; the failing labels re-search all levels).
Variants (p1_open / lambda_open / pbr322_open): `inside=0` drops 21 / 3 / 2 (shows what the inside buys); `merge=1` merges only 2 runs on P1 and Lambda (few equal
neighbours), drop 4 / 1 / 0; `clear=6`: drops 5 / 2 / 0, crossings 11 / 2 / 0, leader mean 61 / 56 / 38 px, canvas 1592x941, 1551x902, 1151x880.
Independent check (`check_g.py`, from the saved label list and the layout radii): no box overlap, no box on ring/track/tick band/title, no leader through
ring or track, no leader through any label text (extra fan leaders included) on pbr322, lambda, p1 open, p1_f14, merge, clear6.

## 4. Assessment (from the figures)
- Open-centre records are where G works: pbr322, mito, plasmid37, lambda are clean (0 to 1 crossings, 0 to 2 dropped), labels sit about 30 to 55 px from their
  feature (today 88 to 120), leaders short and straight, canvas 6 to 10 % smaller than today in both directions (lambda 1347x879 vs 1499x910, P1 1551x943
  vs 1706x990). Staircases appear by themselves in dense stretches, with the inside taking a run (P1: 26 inside, lambda 9). P1 reads like the Geneious
  figure; lambda similar. It has no overprint where today has 20 pairs (P1).
- Where it falls short: (1) dense stretches: straight, non-crossing, clearance-respecting leaders form "walls"; P1 open loses 4 labels and has 14
  crossings (the last-resort pass), dense default-track records lose 20 to 40 % of the labels (p1_f14 19, tobacco_f14 43, T5 and plasmid153 77 to 79 of about
  190). C, D and E place every label; D and E never cross leaders. G is the wrong tool for 100+ labels at 14 pt without free inside space. (2) With gbdraw's
  default tracks (GC + skew) the inside is closed: only a 15 px strip between the GC track and the tick marks exists, so G is an outside-only variant of C
  with shorter leaders and (here) more drops. (3) Some far-flung labels remain after the last-resort pass (lambda max 269 px, P1_f14 395 px). (4) Inside
  labels are never "close": the tick-mark band and tick labels force at least 15 px, usually 40 px between ring and an inside label (Geneious has no tick band there).
  (5) `clear=3` makes labels sit almost on the ring; `clear=6` looks calmer but loses labels and lengthens leaders.
- Order is lost (labels follow geometry, not clockwise lists as in D and E). Merging is rarely useful on these records and "(x2)" is not the Geneious
  "text,text"; colour works but light-blue text on white is hard to read for the default palette (see `lambda_open_f14__G_color`).

## 5. What production would need
Radii and tick/title data from `radial_layout` (already available at the hook; I pass the layout), a per-label colour field (2 lines in the drawer plus a
colour in the label dict), polyline leaders (the prototype's renderer branch), a clearance parameter, and a robust fallback so that no label is dropped
(a proper second row or an ordered boundary-labelling pass for dense stretches, e.g. fall back to D/E-style stacks locally). A reimplementation without numpy
needs a spatial grid for speed. Opt-in per record type (no GC tracks, under about 100 labels); size about 400 lines plus tests; risk: label order and
reproducibility across fonts.

## 6. Not done / approximations (clearance check)
Covered by construction and by `check_g.py`: other label boxes (gap `clear`), other labels' leaders (leader pad `clear`, crossings counted, 1 px in the last
resort), the feature ring on both strands (one annulus from the features' band), tick-mark band, GC/skew tracks, centre title (as the layout's reserved disc, a
conservative circle), tick labels (modelled boxes, auto tick interval only; a manual `--tick_interval` is not read; not checked from the SVG). Not covered:
the legend (it is placed beside the canvas, not checked), clearance between two parallel leaders, minimum leader-to-feature clearance of leaders passing other
features' outer lane (hidden radial part), unusual layouts (outside tracks). Not done: no fallback that avoids all drops; no timing fix for 190-label records;
colour only visible where the palette has several colours (the three showcases have one CDS colour); the "T(-), T(+)" style merge of non-equal text.
