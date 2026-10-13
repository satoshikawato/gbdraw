# Design F: minimum total leader length, never crossing

File: `gbdraw/labels/proto/design_f.py` (only file touched; numpy only). Renders: `work/F/<case>__F.*` (main), `__F_ell` (today's elbow radius, band from 1.0),
`__F_sq`, `__F_max`, `__F_c6` (clear=6), `__F_free` (comparison without the no-crossing rule). Table script `work/F/tab.py`.
Scope changes during the run: (1) free placement became the main variant, (2) leaders must not cross. I could not solve free 2-D placement
under the no-crossing rule; the main variant is the family that is provably crossing-free (below). The free heuristic is kept only as a comparison.

## 1. The rule (main variant, deterministic, 0.6 to 4 s per case)
- Leader = anchor -> radial stub to the elbow E on one common circle of radius r_E (r_E = outer feature radius + 4; `ellipse=1`: today's r_elbow)
  -> straight run to a point A on a common ellipse -> short horizontal stub to the text. Text horizontal, left-aligned on the right half, right-aligned on the left.
- Labels of one side keep the clockwise order of their features and stand in disjoint rows (pitch = (h_a+h_b)/2 + clear). A row is a y on a 1 px grid.
  The text box is pushed out until it lies wholly outside the ellipse, so runs (inside the ellipse) can never touch text; stubs sit in their own row.
- A run must not dip inside the circle r_E: A.E >= r_E^2 (the one-sided tangent bound that D also uses). Per label that limits the angular move to about 20 to 30 degrees.
- Side per label (left or right near 12 and 6 o'clock) = one cut at the top and one at the bottom (every cut within 60 degrees tried): the right side is one contiguous run, the left side the rest.
- Objective: sum of whole-polyline lengths (`obj=sum`), sum of squares (`sq`), or minimise the maximum then the sum (`max`: bisection on a length cap, then sum under the cap).
- Solver: for each side an order-preserving dynamic programme, dp_i[g] = cost_i[g] + min over g' <= g - k of dp_{i-1}[g'] (k = row pitch in px; prefix minimum, O(n x grid)).
  **Exact (provably optimal) for the 1 px row grid, the chosen cuts, the ellipse and the rules above.** Cuts are enumerated exhaustively, so also exact over cuts.
- When no assignment exists the ellipse (both radii, today's aspect) grows in 0.05 steps from the smallest ellipse that clears r_E (default `tight=1`) up to the first feasible band; then 4 more
  bands are solved and the one with the smallest sum is kept. Nothing is dropped, the font never changes.
- **Why crossings are 0 (proof sketch, verified by lp_metrics on every render):** E_i lie on one circle, A_i on one convex closed curve, the order of E and of A agrees, and every run keeps
  radius >= r_E and is angle-monotone. Parameterise runs by polar angle: run i lies outside run j at the start of their common interval and still outside at its end (run j lies inside the ellipse,
  which contains A_i), so they cannot swap. For labels whose angular ranges do not overlap there is nothing to cross. Straight-leader "swap" arguments (Bekos 2007) do not apply directly because our
  leaders are elbowed and boxes differ in width, but this nesting argument holds for them. Cost of the guarantee: the angular move bound, hence band growth in dense records.
- Free variant (`free=1`, comparison only): no ellipse; boxes anywhere outside r_elbow+1+clear, clearance to boxes and leaders, greedy (angular sweep seeded from the ellipse solution) plus relocation and "ripple"
  (re-place a violating label first, then its blockers) moves, penalties for violations. Heuristic, not exact, 6 to 25 s, may leave 1 to 5 crossings. Lower bound (every label at its nearest clear place) is in notes `lb_sum`.

## 2. Parameters
obj sum | clear = ctx.gap 3 px (also row gap) | tight 1 (r_E = ring + 4) | band start auto, step 0.05, `extra` 4 | dy 1 px | cut 60 degrees | maxcut 36 per end | `ellipse=1`: band0 1.0, extra 0, today's r_elbow.

## 3. Metrics (default F; overprint, overlap, crossings, through-text, leaders over features are 0 in all 15 cases; 0 dropped)
leader_len_sum: F, today, B, D, E, R (radial), A (today with fallback). band = ellipse radii relative to ctx ellipse (today's nominal; for mito, lambda, pbr322 the nominal ellipse is smaller than r_elbow, so band 1.1 to 1.2 is the minimum there).
| case | F | cur (crossings, overprint) | B | D | E | R | A | band | mean/max leader | canvas |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| p1_f14 | 27,594 | 12,731 (21, 20) | 17,868 | 12,526 | 9,163 | 14,213 | 19,754 | 1.25 | 249/411 | 1860x1143 |
| p1_f8 | 11,749 | 10,930 (0, 0) | 8,557 | 8,843 | 7,114 | 11,405 | 10,930 | 1.0 | 110/215 | 1375x917 |
| t5_f8 | 70,473 | 67,238 (2085, 0) | 51,603 | 40,309 | 12,093 | 47,305 | 39,295 | 1.5 | 363/555 | 1784x1358 |
| t5_f14 | 136,199 | 65,083 (2532, 170) | 95,762 | 66,286 | 12,401 | 92,428 | 66,304 | 2.25 | 699/926 | 2651x2022 |
| n4_f14 | 12,222 | 8,544 (0, 0) | 6,272 | 6,311 | 5,937 | 7,286 | 8,544 | 1.1 | 177/302 | 1444x1012 |
| tobacco_f8 | 20,545 | 17,043 (0, 0) | 12,227 | 12,282 | 8,532 | 35,147 | 17,043 | 1.05 | 147/263 | 1429x962 |
| tobacco_f14 | 51,699 | 37,840 (1062, 50) | 26,815 | 19,830 | 13,438 | 69,480 | 32,958 | 1.6 | 369/466 | 2154x1449 |
| mito_f14 | 3,098 | 2,898 (0, 0) | 2,700 | 3,752 | 3,161 | 3,175 | 2,898 | 1.1 | 94/160 | 1476x878 |
| plasmid37_f14 | 5,214 | 5,411 (1, 0) | 4,430 | 6,023 | 5,031 | 5,050 | 5,411 | 1.1 | 116/211 | 1653x887 |
| plasmid153_f8 | 40,576 | 18,516 (0, 0) | 16,165 | 16,872 | 13,067 | 18,046 | 18,516 | 1.35 | 227/313 | 1680x1225 |
| plasmid153_f14 | 102,601 | 16,536 (17, 143) | 59,695 | 45,718 | 16,493 | 41,983 | 45,924 | 2.15 | 573/666 | 2622x1930 |
| pbr322_f14 | 3,163 | 3,050 (0, 0) | 2,453 | 3,936 | 2,856 | 2,890 | - | 1.1 | 93/155 | 1241x884 |
| lambda_f14 | 9,635 | 8,171 (2, 0) | 6,553 | 8,256 | 7,654 | 7,080 | - | 1.15 | 136/235 | 1602x923 |
| mg1655_rna_f8 | 11,525 | 11,789 (0, 0) | 9,231 | 10,174 | 3,332 | 11,510 | - | 1.0 | 107/227 | 1151x918 |
| epifagus_f14 | 5,739 | 5,878 (0, 0) | 4,404 | 4,796 | 4,667 | 5,919 | - | 1.0 | 108/218 | 1480x921 |
(Today's sum at t5_f14 / plasmid153_f14 / tobacco_f14 is small only because today overprints 170 / 143 / 50 label pairs there.)
Variants (sum, F_ell = today's elbow radius, band from 1.0): t5_f8 70,945; p1_f14 28,472; tobacco_f14 51,879; lambda_f14 11,098; mito_f14 3,786; p1_f8 12,186. sq: within 0.1 % of sum on all 5 cases; max: max leader
lambda 233.6 vs 235.0, pbr322 151.8 vs 154.6, sum +0.1 to +0.5 %; p1_f14, t5_f8, tobacco_f14 identical. clear=6: p1_f14 35,699 (band 1.45, 2029x1320), lambda_f14 12,797 (band 1.25, 1678x1003).
Free comparison (crossings/through-text allowed to remain; earlier seed): t5_f8 64,875 (5 x, 4 tt), p1_f14 11,953 (1, 1), tobacco_f14 40,543 (0, 9), lambda_f14 7,849 (1, 0), mito_f14 2,533 (0, 0); 6 to 16 s.
Lower bound of the free model (nearest clear place for each label alone): mito 2,314 (free 2,533 = +9 %), p1_f8 4,047, lambda 5,279.

## 4. Assessment (looked at p1_f14, lambda_f14, p1_f8, mito at 1:1 and reduced)
- Looks: tidy fans of straight grey leaders, text on two arcs hugging the ring, no crossings, no leader through text, order exact. It is today's picture with the packing solved, not a new style.
- Does minimising length give a better figure than today on easy records? Only slightly. On mito, plasmid37, pbr322, epifagus, mg1655, p1_f8 the totals are -4 % to +8 % against today (mito +7 %, pbr322 +4 %, p1_f8 +8 %): today is already near the minimum
  there; F's gain is the guarantees (0 crossings), not length. Against D, F is shorter on mito (-17 %), pbr322 (-20 %), plasmid37 (-13 %), lambda (+17 %: longer), epifagus (+20 %).
- On overflowing records F is not good enough: P1 14 pt 27.6k against D 12.5k, tobacco 14 pt 51.7k against 19.8k, T5 8 pt 70.5k against 40.3k, plasmid153 14 pt 103k against 45.7k. Cause: the no-dip bound limits each label to about 25 degrees of angular move
  and only one ring of text exists, so dense clusters force the whole ellipse to grow (band 1.25 to 2.25, canvas up to 2651x2022). D gets its shorter totals from the vertical comb at 12/6 o'clock and from columns taller than the ring; F has no such escape.
  Leaders are then 250 to 700 px. P1 14 pt: no drop, no overprint, 0 crossings, but 2.2 x D's total and 15 % wider canvas than D. T5 8 pt: 0 crossings (today 2,085) at band 1.5.
- Total vs fairness: minimum total length is dominated by labels sacrificed in dense clusters (mean 249 px, max 411 at P1 14 pt). `sq` and `max` change almost nothing because the angular bound, not the objective, decides.
- Free placement (comparison): on P1 14 pt it reaches 11,953 (-57 % vs F, below D) with 1 crossing and 1 through-text left; so the length is available when leaders may take any route, but I could not make that crossing-free.
- Cost: 0.6 to 4 s per case (10 s budget met); canvas on easy cases within about +5 % / -5 % of today's width (mito 1476 vs 1476, p1_f8 1375 vs 1441), much larger on dense ones.
- Ship which objective: `sum`. `sq`/`max` are indistinguishable here. Clearance: 6 px costs +30 % length and +7 % canvas on P1 14 pt; 3 px is enough.
- Right tool: 30 to 110 labels where today is clean (use it for the no-crossing guarantee). Wrong tool: dense records (use D, or radial); as a default for dense records only with a comb or second ring.

## 5. What production would need
`prepare_label_list` hook plus a polyline leader in the circular drawer (exists in the study scaffold); one module of about 250 lines (DP + band search, numpy; gbdraw has numpy via pandas). Canvas already grows from label bounds. Risks: ticks and legend are not
modelled (everything sits outside the ring, inside the legend's column); embedded/inner labels untouched; the proof depends on the common elbow circle (today's per-lane elbows would have to be unified);
`--circular_label_spacing` maps to `clear`.

## 6. What I could not do
- Free 2-D placement with a hard no-crossing guarantee (the heuristic leaves 1 to 5 crossings or through-text; its worst cases fall back to the exact layout when boxes would overlap). A D-like comb or a second ring in the DP would be the way.
- Capacity fact "largest font at band 1.0": only P1 gives a number (8.0 pt; 14 pt needs band 1.25). For tobacco, plasmid153, lambda, T5, no font down to 5 pt fits at band 1.0 of the nominal ellipse, because the angular bound,
  not row capacity, is binding (mito 33 labels needs band 1.1). So the exact solver says the nominal ellipse is the wrong yardstick for a crossing-free layout.
- No ticks/tick-label or legend clearance (none lie outside the ring in these figures). Free variant `obj=sq|max` not rendered. Free renders were made before the `tight=1` default and are not regenerated.
