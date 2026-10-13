# Design B: staggered tiers (report)

File: `gbdraw/labels/proto/design_b.py` (about 220 lines, deterministic, 0.6 to 1.5 s per case). Renders and metrics:
`work/B/<case>__B.*`; looks: `work/B/view_*.png` (0.55 scale) and 1:1 crops `crop_*_1to1.png`. `__B_drift_320` files are the
same design with a larger drift bound (fewer tiers, taller figure).

## 1. The rule
- Each half of the map is its own set of stacks (right: theta in [0, 180], left: mirrored). Text stays horizontal, start-anchored on the right,
  end-anchored on the left. Labels of a half are ordered by theta; every stack keeps that order, so leaders cannot cross.
- Tier 0 follows today's label ellipse (floored at r_keep + 12 so the band clears the ring). Label i has a natural height y_nat = -b cos(theta).
  The stack is solved with `pack_ordered` (least squares to the natural heights, gap 3 px). It may run beyond the ellipse top or bottom
  (rows above 12 o'clock and below 6 o'clock). A label above its own elbow starts at x = x_elbow + 2, so its leader is vertical and the
  top and bottom rows form a clean staircase.
- Tiers 1, 2, ... are flush columns. Column t starts at (right end of the widest label of tier t-1) + 8 + corridor (36 px). Every label that lives
  outside tier s adds a 4 px slot to the stack of tier s (in theta order), so the stack of tier s is packed from its own rows plus those slots.
- Leader of a tier-t label: feature anchor -> elbow ring (today's radial stub) -> to the left edge of tier 0 at the height of its slot
  -> horizontal run through the slot (a gap between two rows) to the right end of tier 0 -> straight "bus" line to the left edge of tier 1 (its
  label, or its slot there) -> ... -> label. All entry points of a band lie on one curve (tier 0) or one vertical line (outer tiers), which is what makes
  the lines between bands non-crossing. Segments that would cut inside the feature ring get waypoints on a circle just outside it (over-feature leaders: 0 everywhere).
- Which tier: for each half, take the fewest tiers T (<= 4) such that a sweep in theta order keeps every label within `drift` (240 px) of its natural height. In the sweep a label stays
  in the current tier while its stack still fits, otherwise it moves out. Width-aware part: tier t only accepts labels whose width is <= the (t+1)/T quantile
  of the half's widths, so inner tiers hold short names and the outermost tier takes the long ones (this keeps the corridors narrow: it cut plasmid153_f14 from
  3294 to 2203 px wide).
- When nothing fits: the drift bound is relaxed in x1.5 steps (240, 360, 540, ...) until some T fits; nothing is dropped (0 dropped everywhere). Easy
  cases get T = 1 and look like today.
- Bounds: label to feature is bounded vertically by `drift` along its own tier, but the leader of an outer label also runs through the inner tiers, so a leader can reach
  1.5k px on T5 at 14 pt (not bounded in length).

## 2. Parameters (defaults)
`drift`=240 px, `tiers`=4 (max), `slot`=4 px, `corridor`=36 px, `wq`=1 (width quantile caps, 0 = off), `xfloor`=8, `routepad`=6. Set with `--opts`.

## 3. Metrics (default settings)
| case | labels | dropped | overprint | overlaps | crossings | through-text | over features | leader mean/max px | canvas WxH | today: overprint / crossings / canvas |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| p1_f14 | 111 | 0 | 0 | 0 | 0 | 0 | 0 | 161 / 596 | 1732x1083 | 20 / 21 / 1706x990 |
| p1_f8 | 107 | 0 | 0 | 0 | 0 | 0 | 0 | 80 / 172 | 1372x1068 | 0 / 0 / 1441x982 |
| t5_f8 | 194 | 0 | 0 | 0 | 0 | 0 | 0 | 266 / 1182 | 1607x1344 | 0 / 2085 / 1444x979 |
| t5_f14 | 195 | 0 | 0 | 0 | 0 | 0 | 0 | 491 / 1487 | 2290x1414 | 170 / 2532 / 1821x992 |
| n4_f14 | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 91 / 251 | 1337x1147 | 0 / 0 / 1402x994 |
| tobacco_f8 | 140 | 0 | 0 | 0 | 0 | 0 | 0 | 87 / 204 | 1349x1109 | 0 / 0 / 1399x983 |
| tobacco_f14 | 140 | 0 | 0 | 0 | 0 | 0 | 0 | 192 / 772 | 1867x1252 | 50 / 1062 / 1634x988 |
| mito_f14 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 82 / 129 | 1467x930 | 0 / 0 / 1476x906 |
| plasmid37_f14 | 45 | 0 | 0 | 0 | 0 | 0 | 0 | 98 / 132 | 1644x960 | 0 / 1 / 1648x916 |
| plasmid153_f8 | 179 | 0 | 0 | 0 | 0 | 0 | 0 | 90 / 252 | 1373x1231 | 0 / 0 / 1436x983 |
| plasmid153_f14 | 179 | 0 | 0 | 0 | 0 | 0 | 0 | 334 / 1032 | 2203x1352 | 143 / 17 / 1714x994 |

All required: overprint 0 and overlaps 0 on every case, no leader through text, 0 dropped. Crossings are 0 everywhere. Tiers used (R/L):
p1_f14 R[39] L[37,35]; t5_f14 R[36,32,50,22] L[55]; plasmid153_f14 R[44,42] L[58,35]; tobacco_f14 R[39,37] L[64]; t5_f8 R[78,62] L[54]; the
others are single-tier. Label to feature distance (mean/max px): p1_f14 148/517, t5_f14 304/881, plasmid153_f14 235/614, n4_f14 84/243.

## 4. Assessment (from the figures)
- Better than today: every hard case is clean (today p1_f14 has 20 overprint pairs, t5_f14 170, tobacco_f14 50, plasmid153_f14 143, and T5 has 2000+ crossings). The look is
  structured: a ragged inner band of short names and aligned outer columns of long names, with order preserved (reading order follows the genome).
  tobacco_f14 and plasmid153_f14 are the best figures: the inner column holds "hypothetical protein" / gene names, the outer one the long product names.
- Easy cases (mito, plasmid37, n4, p1_f8, tobacco_f8, plasmid153_f8) get one tier and look like today or slightly cleaner (n4_f14 leader mean 91 px vs 124 today).
- Long product names vs short gene names: short names are cheap (p1_f14, mean width 104 px: two tiers, 1732 x 1083). Long names make every extra tier cost
  width: without the width-aware split plasmid153_f14 needed four columns (3294 px wide); with it two (2203 px). A record whose names are all long and that
  needs 3 or more tiers (T5 at 14 pt, mean 132 px, max 475) ends up 2290 x 1414 and the outer columns are far from the ring.
- Size: width at 14 pt grows with the tiers used (p1 +26 px, tobacco +233, t5 +469, plasmid153 +489 against today; 8 pt cases are within +-100 px), height only where a stack overshoots the ring: 930 to 1414 px
  against today's 906 to 994. Today is 1399 to 1821 wide. T5/14 pt is the extreme: 2290 x 1414 (today 1821 x 992, and overprinted).
- Slot leaders at 100 %: with 4 px pitch each leader is a distinct thin line and can be followed through a gap in groups of up to about 30 lines (see
  `crop_p1_f14_left_1to1.png`); in the densest bundles (T5 right half, 60 to 100 parallel lines) they merge into a grey band and only the label end is readable.
  Leaders get long: mean 161 px (p1_f14), 491 px (t5_f14); max 1487 px.
- Weak points: (1) one-sided records (T5: 140 labels on the right, 55 on the left) overload one half; (2) the right half of T5 reads as a wall of bundles;
  (3) labels can sit 500 to 900 px from their feature; (4) tiers add width, so a landscape figure with the legend further right.
- Right tool: 100 to 200 labels at 12 to 14 pt where names have mixed lengths (plastids, plasmids, phages); wrong tool for very dense one-sided clusters of long names
  (use a smaller font or filter, or C/E).

## 5. What production would need
New module under `gbdraw/labels/` (about 250 lines) called from `prepare_label_list` for horizontal placement when the existing result has overlaps (or as the placement); polyline leaders
in `render/groups/circular/labels.py` and `assemble.py` bounds (already prototyped in the scaffold); config for drift and tiers. Opt-in at first
(`--label_placement tiered`). Risks: canvas growth and legend position; the Auto font size already shrinks; Playwright/editing code that assumes a 2-point leader; test fixtures.
It needs no change to filtering or measurement.

## 6. Not done / notes
- Balancing the two halves (moving the 12 and 6 o'clock cut by up to 25 to 45 degrees so T5's right half shares load) was built and removed: it halves T5's width
  but leaders that cross the axis gave 28 to 150 crossings; fixing it needs routing over the top. Left as the main open improvement.
- Label height h is the measured box (13.7 to 17.4 px), so row pitch is slightly uneven; a uniform pitch would cost about 20 % height.
- Observation for the scaffold (not changed): in the PNG (cairosvg) the end-anchored label at the left edge is clipped by about 16 px ("TciB-like ..." loses its first letters in p1_f14),
  because the rendered text is wider than the measured width (about 5 to 9 %), and canvas bounds use the measured width. It affects every design.
- Side observation: "hypothetical protein" labels appear although `hypothetical` is in the default blacklist (T5, plasmid153): not investigated.
