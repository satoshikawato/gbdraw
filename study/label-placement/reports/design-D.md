# Design D: side columns (boundary labelling)

File: `gbdraw/labels/proto/design_d.py` (about 150 lines, deterministic, 1 to 2 s per case). Renders: `work/D/<case>__D.*`
(default = "D", bend variant), `work/D/<case>__D_flat.*` (strictly vertical columns, option `flat=1`), `work/D/*_mix.*`.
Tables: `work/D/summary_D.txt`, `summary_Dflat.txt` (raw `run_matrix` lines).

## 1. The rule
- Every external label goes to a right or a left column; text is horizontal, left-aligned on the right, right-aligned on the left.
- Side = side of the feature (sin theta >= 0: right). Near 12 and 6 o'clock a feature may switch side (`cutmax`, 40 degrees) to
  balance the columns: one cut near the top and one near the bottom, chosen by exhaustive search (about 30 x 30 candidates) to
  minimise the sum of squared vertical leader slopes. Each column is one contiguous clockwise run, kept in clockwise order
  (right: top to bottom, left: bottom to top). That order is what makes leaders non-crossing.
- Vertical position: `pack_ordered(targets, h + gap)` with target = y of the leader stub on the circle (stub radius r_elbow,
  the stub of today's leader). The column grows taller than the circle when it must: nothing is dropped, nothing shrinks.
- Leader = feature -> radial stub to r_elbow -> straight run -> 8 px horizontal stub into the label. Two rules keep it
  crossing-free and outside the circle of radius r_elbow (verified: no run goes inside it in any case):
  1. A label may lie below its stub (upper half) only as far as the tangent at the stub allows (above it in the lower half): a
     one-sided bound on y, applied after packing in one backward/forward pass.
  2. A feature on the "wrong" side (moved across 12 or 6 o'clock) has its label beyond the circle's top or bottom; its leader goes
     straight up/down from the stub, then sideways ("comb" above and below the ring). Wrong-side items exist at one end on one side only.
- Column x (variant D, `flat=0`): labels hug a circle of radius X = r_elbow + 24 (x0 = sqrt(X^2 - dy^2), dy = distance of the
  box's near edge from y=0), but never closer than the stub x + 21 px, and x0 never decreases towards the top/bottom, so no
  leader reaches its label through another label's text. Variant D-flat: x0 = X everywhere.
- When labels do not fit: nothing is dropped; the figure gets taller (column height = labels per side x pitch). No capacity limit.

## 2. Parameters
flat 0 (D) | mix 0 (target = stub y; 0.6 tested, no gain) | xoff 24 px (X - r_elbow) | hs 8 px | cutmax 40 degrees | gap = ctx.gap (3 px).

## 3. Metrics (all standard cases, default D). "cur" = today (overprint pairs / crossings / leader mean/max px / canvas)
| case | labels | dropped | overprint | overlap | crossings | through-text | over features | leader mean/max | label-feature mean/max | canvas w x h | cur: overprint / crossings / leader mean/max / canvas |
|---|---|---|---|---|---|---|---|---|---|---|---|
| p1_f14 | 111 | 0 | 0 | 0 | 0 | 0 | 0 | 112.8/367.9 | 96.6/257.3 | 1624x1179 | 20 / 21 / 114.7/467.3 / 1706x990 |
| p1_f8 | 107 | 0 | 0 | 0 | 0 | 0 | 0 | 82.6/174.3 | 74.1/142.0 | 1361x961 | 0 / 0 / 102.1/141.6 / 1441x982 |
| t5_f8 | 194 | 0 | 0 | 0 | 0 | 0 | 0 | 207.8/623.3 | 174.3/427.0 | 1359x1344 | 0 / 2085 / 346.6/699.2 / 1444x979 |
| t5_f14 | 195 | 0 | 0 | 0 | 0 | 0 | 0 | 339.9/922.5 | 293.7/695.1 | 1610x1965 | 170 / 2532 / 333.8/899.4 / 1821x992 |
| n4_f14 | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 91.5/270.8 | 79.7/197.3 | 1322x1025 | 0 / 0 / 123.8/230.2 / 1402x994 |
| tobacco_f8 | 140 | 0 | 0 | 0 | 0 | 0 | 0 | 87.7/198.3 | 79.7/174.8 | 1334x1035 | 0 / 0 / 121.7/225.8 / 1399x983 |
| tobacco_f14 | 140 | 0 | 0 | 0 | 0 | 0 | 0 | 141.6/388.9 | 126.1/372.5 | 1555x1417 | 50 / 1062 / 270.3/705.4 / 1634x988 |
| mito_f14 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 113.7/195.2 | 101.5/155.4 | 1490x899 | 0 / 0 / 87.8/127.3 / 1476x906 |
| plasmid37_f14 | 45 | 0 | 0 | 0 | 0 | 0 | 0 | 133.8/236.9 | 118.2/173.7 | 1666x921 | 0 / 1 / 120.2/248.0 / 1648x916 |
| plasmid153_f8 | 179 | 0 | 0 | 0 | 0 | 0 | 0 | 94.3/251.8 | 85.0/208.0 | 1364x1190 | 0 / 0 / 103.4/207.3 / 1436x983 |
| plasmid153_f14 | 179 | 0 | 0 | 0 | 0 | 0 | 0 | 255.4/582.2 | 236.0/526.2 | 1650x1857 | 143 / 17 / 92.4/177.2 / 1714x994 |

Zero overprint, overlap, crossings, leader-through-text and leaders over features in every case, zero dropped. tobacco_f14, n4_f14 and
plasmid153_f14 cur values rendered into `work/base/` for this table. Variant D-flat (same zeros) is longer-leadered:
p1_f14 leader mean/max 243/800 vs 113/368 (D), canvas 1694x1144 vs 1624x1179; t5_f14 453/1355 vs 340/923, canvas 1791x1955 vs 1610x1965;
mito_f14 231/485 vs 114/195. Full table: `work/D/summary_Dflat.txt`.

## 4. Assessment (from the figures)
- Distance: on average a label sits about 75 to 130 px from its feature on 30 to 180 labels (mean label-to-feature 74 px at P1 8 pt, 97 px at P1 14 pt),
  similar to today where today is clean. At 190 labels and 14 pt (T5, plasmid153 f14) the mean is 236 to 294 px (max 526 to 695): the
  column is about 1.9 x the circle's height, so labels near 12 and 6 o'clock are pushed far out. Today's T5 mean is 334 px but with
  2,532 leader crossings and 170 overprinting pairs, so D is shorter and clean.
- Size: width is today's or smaller (D hugs the ring; mostly -4 % to -11 %, +1 % on mito/plasmid37). Height is the cost: easy and moderate
  cases +0 to +10 % (P1 14 pt +19 %, because 111 labels need 1,140 px of stack); P1 8 pt -2 %; tobacco 8 pt +5 %, tobacco 14 pt +43 %;
  plasmid153 8 pt +21 %; T5 8 pt +37 %; T5 14 pt and plasmid153 14 pt about +90 % (height 1,965 / 1,857 px, wider than tall).
- Reads well as a legend: yes. Two aligned, clockwise lists; the order is exact; easy to scan by eye along the column and to find
  a gene by position. Reading order is lost only at the cuts (12 and 6 o'clock) where list start/end meet.
- Easy cases (mito, plasmid37, P1 8 pt, tobacco 8 pt): looks like today's figure but cleaner (labels follow the ring, leaders short,
  never crossing). D is a straight improvement there; I would not call it worse anywhere I looked.
- Weak spot, measured: the combs above and below the ring. 12 and 6 o'clock features get vertical-then-sideways leaders, the
  top/bottom label block stacks up in a parallel-line comb (T5 left top, 35 labels). Readable, not pretty; it grows with N.
- Finding about design A: variant D with `flat=0` is "today's ellipse with an exact packer and a height that grows". Without the
  tangent bound and the comb routing (rules 1 and 2) the exact packer on the ellipse still crosses leaders (my first try: 48 to 135
  crossings on P1/T5). So today's look plus an exact stack is not enough; the leader routing at the top and bottom is the real work.
- Strictly vertical columns (D-flat) are defensible only as a tidy "table" style: long fans at 12/6 o'clock (P1 max leader 800 px).
- Right tool: 30 to about 200 labels at modest font size, complete labelling required, figure height may grow (genomes with many CDS
  labelled, P1 at 14 pt, plasmids, plastids). Wrong tool: very dense records at large font (T5 14 pt: 1,965 px tall figure),
  figures that must keep a fixed aspect ratio or page size, records with few labels (no gain over today), and very long label text
  (the columns carry its full width on both sides).
- Two columns per side: not done. A second column further out needs its leaders to thread between the rows of the first (that is
  design B), or to run beyond the first column's text, which costs width instead of height. Not attempted.
- Known harness quirk (not D): in the PNG the longest P1 label (322 px measured) draws ~6 % wider than measured and is clipped at the
  left edge of the canvas; in the SVG `labels_bbox` fits. Same measurer for all designs.

## 5. What production would need
- `gbdraw/labels/circular.py::prepare_label_list` gets a placement mode `side_columns` next to horizontal/radial (opt-in
  `--label_placement columns`); the placer is one module of about 150 lines, using the same item/box contract. The renderer already draws
  leader polylines with extra vertices (study scaffold); production would need a polyline leader in the circular label drawer and
  canvas growth from label bounds (already how canvas fitting works). Legend placement may need a left/right choice.
- Risks: every leader-routing rule above is geometry-sensitive (tangent bound uses the stub, not the feature ring); long labels widen
  the figure both sides; tall figures change page layout (docs, Gallery thumbnails); `--circular_label_spacing`, inner labels (`--labels
  both`) and `--label_font_size` interact (inner labels not handled here). Fine to ship as opt-in; as the default only with a height cap.

## 6. Not done
Inner labels (`--labels both`), embedded labels untouched (as spec). No second column per side. `mix` (rank-spread targets) tested
once (mix 0.6: no gain) and left at 0. No test of other `track_type` or `--outer_label_*_radius_offset`. Metrics `leader_len_*` count the
whole polyline (stub + run), so values are larger than the visible "run".
