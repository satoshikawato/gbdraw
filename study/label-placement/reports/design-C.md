# Design C: free placement (collision-driven search)  -  prototype report

File: `gbdraw/labels/proto/design_c.py` (about 440 lines). Renders: `work/C/<case>__C.{svg,png,metrics.json,labels.json}`
(default), `<case>__C_strict.*` (no crossing fallback), `<case>__C_lmax{0.35,0.5,1.0}.*` (bound study),
`summary_default.txt`, `table.py` (prints the tables below). Flattened looks: `work/C/*_flat.png`.

## 1. The rule
A label is a horizontal box anywhere outside the ring. Its leader is: feature anchor -> radial stub on the feature's
ray just outside the ring (today's elbow radius `r_elbow`) -> ONE straight segment to the nearest edge point of the
box (slid 4 px away from corners, 2.5 px off the text). The segment must leave the stub outward (at most `amax`
from radial), so it never cuts back over the ring or over neighbouring stubs.
1. Reserved slots. Each label gets a conflict-free fallback slot from an ordered contour packing (each half of the
   map sorted by stub y, `pack_ordered`, beside today's label ellipse, pushed outside `r_keep + band`). Slots that
   conflict (overlap, leader through a box, leader crossing, leader longer than `lmax`) are removed from the
   reservation by greedy vertex cover (lowest `span_bp` first). Reserved slots are obstacles for everyone else.
2. Free search. Labels are processed most crowded first (neighbours within +-0.12 rad). A label releases its slot,
   enumerates candidate box centres on a grid (3 px; 4 px above 130 labels) within `lmax` of its stub, sorted by cost =
   leader length + 0.5 * length * (1 - cos tilt) + 0.1 * R * |angular drift| (+ 100 per unplaced stub corridor the
   box or leader blocks: radial corridor of depth 60 px; capped at 8). The first collision-free candidate in cost
   order wins: no box overlap (gap), no leader through any box (1 px pad), no placed leader through the box, no
   crossing. It keeps its slot when no free box is cheaper than the slot's leader length.
3. Fallback. Labels with neither get one more search that may cross at most 2 leaders (cost 40 per crossing); else
   they are DROPPED (ctx.dropped; `fallback=0` drops instead of crossing).
Bounded: leader length (stub to box edge) <= `lmax`; leader direction <= `amax`. Candidate search is vectorised
(numpy, a 2 px raster prefilter then exact slab tests); deterministic; 3 to 10 s per case incl. render.

## 2. Parameters (options via GBDRAW_LABEL_PROTO_OPTS)
lmax=0.7 (x R = 273 px), amax=85 deg, band=12 px, step=3 (4 if > 130 labels), drift=0.1, tilt=0.5, raypen=100,
raydepth=60, raycap=8, crosspen=40, maxcross=2, fallback=1, slots=1, slotbias=0, order=4 (crowd; 0 = clockwise sweep).

## 3. Metrics (default; Cs = fallback=0; cur = today). labels / dropped, overprint, overlap, crossings, through-text, leaders over features
| case | C lab/drop | op | ov | cross | thr | lof | leader mean/max | canvas | Cs drop / cross | cur: op, cross, canvas |
| --- | --- | -: | -: | -: | -: | -: | --- | --- | --- | --- |
| p1_f14 | 111/13 | 0 | 0 | 15 | 0 | 0 | 91/262 | 1577x1167 | 24 / 0 | 20, 21, 1706x990 |
| p1_f8 | 107/4 | 0 | 0 | 6 | 0 | 0 | 83/226 | 1359x969 | 10 / 0 | 0, 0, 1441x982 |
| t5_f8 | 194/74 | 0 | 0 | 23 | 0 | 0 | 113/285 | 1370x1009 | 87 / 1 | 0, 2085, 1444x979 |
| t5_f14 | 195/100 | 0 | 0 | 26 | 0 | 0 | 131/291 | 1781x1094 | 117 / 0 | 170, 2532, 1821x992 |
| n4_f14 | 69/3 | 0 | 0 | 15 | 0 | 0 | 93/222 | 1356x992 | 13 / 0 | 0, 0, 1402x994 |
| tobacco_f8 | 140/9 | 0 | 0 | 11 | 0 | 0 | 108/269 | 1500x1029 | 15 / 0 | 0, 0, 1399x983 |
| tobacco_f14 | 140/31 | 0 | 0 | 18 | 0 | 0 | 107/274 | 1643x1254 | 44 / 0 | 50, 1062, 1634x988 |
| mito_f14 | 33/0 | 0 | 0 | 0 | 0 | 0 | 89/183 | 1468x924 | 0 / 0 | 0, 0, 1476x906 |
| plasmid37_f14 | 45/0 | 0 | 0 | 1 | 0 | 0 | 113/216 | 1647x1050 | 1 / 0 | 0, 1, 1648x916 |
| plasmid153_f8 | 179/28 | 0 | 0 | 15 | 0 | 0 | 87/273 | 1372x1194 | 36 / 0 | 0, 0, 1436x983 |
| plasmid153_f14 | 179/87 | 0 | 0 | 23 | 0 | 0 | 106/273 | 1708x1106 | 101 / 0 | 143, 17, 1714x994 |
Leader lengths are whole polylines (anchor to label), as in lp_metrics. Zero overprint/overlap/through-text everywhere.
Bound study (placed / dropped, crossings), lmax = 0.35 / 0.5 / 0.7 / 1.0 R:
p1_f14 89/22,5 | 95/16,10 | 98/13,15 | 95/16,9;  t5_f8 125/69,20 | 123/71,20 | 120/74,23 | 120/74,27;
plasmid153_f14 87/92,6 | 84/95,29 | 92/87,23 | 99/80,26. The bound hardly matters: the drops are caused by enclosure
(a stub whose every exit is walled by neighbouring boxes and leaders), not by distance.

## 4. Assessment (from the figures)
- Looks organised where it works. Sparse and moderate records (mito, plasmid37, P1 at 8 pt, the quiet half of the
  others) come out with short straight leaders at home, no ellipse, no leader bundles at 12 and 6 o'clock; P1 at
  14 pt (`p1_f14_flat.png`) reads like a CGView "angled" figure: clean fans, leaders followable by eye, reading
  order kept around the circle because the search is cost-ordered and slots are order-preserving.
- Dense runs are where it fails, and it fails by silently losing labels, which today's placement never does.
  Greedy collision-driven search encloses stubs: a box placed in front of a cluster of stubs walls the neighbours in.
  The reserved slots fix part of this (tobacco_f8: 23 -> 9 dropped, plasmid153_f8: 45 -> 28), not all. On the records
  where today is clean (P1_f8, tobacco_f8, N4, plasmid153_f8) C is WORSE than today: 3 to 28 labels dropped and 6 to 15
  crossings (with `fallback=0` there are no crossings but 4 to 36 drops). On the overloaded 14 pt cases it is better
  than today's overprint but drops 28 to 56 % of the labels (t5_f14 100 of 195, plasmid153_f14 87 of 179).
- In dense regions (`t5_f8_crop.png`) the survivors form near-parallel, near-vertical leader bundles that touch the
  ring side by side; one can follow a leader only with effort, and leaders are the only link (labels sit beside, not
  behind, their genes).
- Reading order: kept locally (fans are monotone), but a label may sit across several neighbours' stubs; clockwise
  reading of the text column is not guaranteed in the fallback and crossing cases.
- Canvas: grows 8 to 18 % in height (P1 14 pt 990 -> 1167, tobacco_f14 988 -> 1254, plasmid153_f8 983 -> 1194), because
  the overflow of a cluster stacks above or below the circle; width is about today's. Leader length is a bit shorter than
  today's on the clean cases (mean about 85 vs 100 px) and much shorter than today's pathological T5 (113 vs 347 px).
- Right tool: modest label counts (up to about 100 outer labels at 14 pt, 150 at 8 pt), uneven gene spacing, where
  short home-placed leaders are wanted. Wrong tool: dense clusters in the polar regions or half-sided records such as
  T5 (120 labels in one half). There it needs the ordered stack of D/E (or radial text) as the structure, with free
  placement only as the local relaxer; as a stand-alone method it is not a replacement for today's placement.

## 5. Clearances (LP-D06: minimum clearance from every other object)
Design C keeps these clearances between a label (box = measured text box) and other objects. All values in px.
| object | clearance kept | where set | checked? |
| --- | --- | --- | --- |
| other label text | 3 (box gap), plus up to 2 px extra from the raster prefilter | `ctx.gap` (= labels.spacing.circular, CLI `--circular_label_spacing`); design_c reads it in `_best` (`g = ctx.gap`), `_conflicts`, `_State._mark` (r_gap). Raster cell: constant `CELL = 2.0` | yes, exact |
| other labels' leaders (leader vs this text) | 1 | literal `1.0` pad in the `_seg_box_hit(..., 1.0)` calls (`_best`, `_conflicts`) and `pad` in `_State._mark`: one concept, three literals, no option | yes, exact |
| this label's own leader to its own text | 2.5 off the box, 4 from corners | literals `pad, corner = 2.5, 4.0` in `_candidates`, `_leader` slot code in `place`/`_contour_targets` call sites (`_attach(..., 2.5, 4.0)`) | by construction |
| other leaders (leader vs leader) | 0 (proper crossing test only; parallel or touching leaders are allowed) | no constant; `_seg_cross_count`/`_seg_cross_matrix` use a 1e-6 orientation tolerance | crossing only, no distance |
| feature ring (outermost lane) | boxes: `r_keep + band` from the centre, i.e. at least 4 + 4 + 12 = 20 px beyond the outermost lane edge (`r_elbow` >= outer + 4, `r_keep` = r_elbow + 4) | option `band` (default 12); `r_keep`/`r_elbow` come from `Ctx.build` in common.py. Leaders: start at the stub (4 px outside the lane), leave outward at <= `amax` (85 deg) | yes; leaders over features = 0 in every render |
| ticks (kbp marks, tick text) | none | nothing in design_c | NOT checked (they lie inside the ring; the box keep-out above covers them only because of the ring clearance) |
| legend / centre text / other canvas items | none | nothing in design_c; the canvas grows to hold the labels afterwards | NOT checked |
| embedded (curved) labels on features | none | embedded labels are not passed to `place` | NOT checked (they lie inside the lane) |
So a variant with larger clearance: `band=24` (ring, one option); label gap needs `--circular_label_spacing` (or one
edit: `g = ctx.gap + extra` in `_best`, `_conflicts`, `_State._mark`); leader-to-text pad needs the three `1.0` literals
changed together; leader-to-leader distance and ticks/legend would need new code (the leader test would use a
distance instead of a crossing test, ticks and legend would be added as extra obstacle boxes to `_State`).

## 6. What production would need / could not do
- Production: a `labels/free_placement.py` with the numpy search (about 300 lines), a new `--label_placement free`
  option beside `horizontal` and `radial`, reuse of `pack_ordered` for slots, and the drawer already supports polyline
  leaders (`leader_points`). Risk: silent label loss; a production version must report dropped labels or fall back to
  A1/A3-style shrinking, and the 5 to 10 s runtime for 190 labels must be cut (the exact slab tests dominate).
- Not done: elbowed leaders (straight only), ejection-chain repair (tried: +2 labels for 10x time, removed), squeaky-wheel
  reordering (no gain), force-directed relaxation as the main engine (tried: 15 to 40 residual conflicts, unusable on
  its own; removed), per-label importance beyond `span_bp`, no check against legend/ticks.
- Harness note: placement counts from the offline harness differ by a few labels from the real pipeline (real arc radii);
  the tables above are from the real renders.
