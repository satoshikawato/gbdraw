# Design A: keep today's look, never overprint (A1 font, A2 band, A3 drop, A combined)

Code: `/home/kawato/gbdraw-work/.worktrees/label-study/gbdraw/labels/proto/design_a.py` (only file I own; nothing else touched).
Renders: `work/A/<case>__{A1,A2,A3,A}.*`; sweeps `work/A/sweeps/*.sweep.json`; 2-D grids `work/A/grid/`; full table `work/A/table.md`.

## 1. The rule
Today's placement is run unchanged through `render_today(font_size, band_scale, drop)`; we only choose its three knobs.
**Fits** (computed from the returned label dicts, `evaluate()`): (a) no two label boxes overlap by more than 0.5 px,
(b) no leader segment passes through another label's box, (c) no leader segment after the first runs inside the
feature ring (`leaders_over_features` = 0), (d) no label's leader end lies more than 30 degrees (angular drift) from its
feature. (d) is what catches T5 at 8 pt: it has no overprint but 92 labels beyond 30 degrees (max 99). Clean today-cases
have max drift 6.7 to 24.3 degrees, so 30 is above everything that looks fine today. Leader crossings are measured, not required.
If today's layout already fits, it is returned untouched (SVG byte-identical to `cur`: confirmed for all 6 cases that fit).
Otherwise a bisection over a monotone path of settings (each probe is a full render of today's pipeline) finds the first
fitting point; the result is always a verified rendered layout.
- A1: font from requested down in 0.5 pt steps to the floor 6 pt (labels re-measured, so some labels become embedded: P1 111 -> 107 outer).
- A2: band_scale 1.0 to 2.0 in 0.05 steps (both radii of the label ellipse).
- A3: keep font/band, drop the k shortest features first (importance = longest location segment in bp, read from the
  feature dict; ties by stable_id); bisection over k. Dropped label texts are returned in `summary["dropped"]`.
- A (recommended): one path in which band and font move together along the measured frontier, then font alone, then drop.
  Step k (0..10): band = 1 + 0.05 k, font = requested * (1 - 0.035 k) rounded down to 0.5 pt (clipped to the floor). After band 1.5
  the font continues down to 6 pt at band 1.5; if still not fitting, drop (shortest features first) at 6 pt / 1.5.

## 2. Parameters
drift 30 deg; font step 0.5, floor 6 pt; band step 0.05, A2 cap 2.0; A: band cap 1.5, trade 3.5 % font per 0.05 band; gap = today's.
Chosen from the 2-D grids (`work/A/grid`): the fit frontier is smooth and roughly linear, about 0.1 band per 1 pt (6 to 7 % of
font) for P1, tobacco and T5, so moving both knobs together keeps each change small. See section 4 for why not "widen first".

## 3. Results (full table with every variant and metric: `work/A/table.md`; ovl, through-text, crossings, over-feature = 0 unless noted)
| case | today (ovl / tt / drift>30) | A1 font | A2 band | A3 dropped | A font / band | largest font that fits (band 1) |
|---|---|---|---|---|---|---|
| p1_f14 | 20 / 8 / 1 | 8.0 | 1.25 | 12 | 12.5 / 1.15 | 8.0 |
| p1_f8, n4_f14, tobacco_f8, mito_f14, plasmid37_f14, plasmid153_f8 | fit | unchanged | unchanged | 0 | unchanged | requested (n4 also 6..14) |
| t5_f8 | 0 / 101 / 92 | 6.0 does NOT fit (47 tt) | 1.4 | 42 | 6.5 / 1.2 | none >= 6 |
| t5_f14 | 189 / 168 / many | 6.0 does NOT fit | 2.0 does NOT fit (5 tt) | 94 | 8.5 / 1.5 | none >= 6 |
| tobacco_f14 | 53 / 110 / many | 9.5 | 1.5 | 39 | 11.5 / 1.25 | 9.5 |
| plasmid153_f14 | 144 / 65 / 0 | 8.0 | 1.7 | 70 | 10.5 / 1.35 | 8.0 |
Canvas (WxH px) today -> A: p1_f14 1706x990 -> 1773x1130 (A2 1928x1226, A1 1441x982); t5_f8 1444x979 -> 1533x1164; t5_f14 1821x992 -> 1886x1448;
tobacco_f14 1634x988 -> 1785x1223 (A2 2113x1462); plasmid153_f14 1714x994 -> 1862x1315 (A2 2338x1649). Leader mean/max for A:
p1_f14 178/222 (today 115/467), t5_f14 342/399 (334/899), tobacco_f14 235/330, plasmid153_f14 257/329.
Fit-region facts (single axis): minimum band at the requested font: P1 14pt 1.25, 12pt 1.15, 10pt 1.1; T5 8pt 1.4; tobacco 14pt 1.5; plasmid153 14pt 1.7.

**Monotonicity.** For every case, `fit` was monotone along both axes on the 0.5 pt / 0.05 grid (all sweeps read `....11111`), and
on the 2-D grids (P1, tobacco, T5-8pt). The raw counts are not monotone while a layout is infeasible: T5 14 pt overlaps 189, 193 (band 1.05), 172, 147...; max drift jumps
26 deg at band 1.05 then 178 deg at 1.1 (the scored winner flips between the legacy and the newer pipeline). So bisection on the
boolean is safe in practice, not guaranteed; the final layout is verified by rendering, only "minimal" is an assumption.
**Surprise 1:** band_scale did nothing at first: on the CLI path `outer_arena` is set, and `rearrange_labels_fc` rescales the ellipse to the arena right
after `_PROTO_STATE["band_scale"]` is applied (circular.py ~4270), cancelling it. Worked around inside design_a.py by scaling `kwargs["outer_arena"]`
(closure of the wrapper). Fix for the lead if A2/A is kept in the wrapper: apply band_scale after the arena block, or scale the arena there.
**Surprise 2:** P1 at 12 to 8.5 pt has no overlap but one leader through a label's text; only 8.0 is clean. That one defect decides P1.
**Surprise 3:** `make_items` (common.py) sets `span_bp = abs(end - start)`, but `start`/`end` of a label dict are the label TEXT extent in bp
(1.1 x text width), not the feature span. Designs B to E that use `span_bp` as importance rank by label width. A3 reads the real feature span.

## 4. Assessment (looked at P1, T5 x2, tobacco, plasmid153 at 1:1 and reduced)
- A1 (font only) keeps the canvas but pays the most legibility: P1 14 -> 8 pt, plasmid 14 -> 8 pt, and it cannot rescue T5 at all (6 pt still has 47 leaders through text).
- A2 (band only) keeps the font but the figure balloons (tobacco 2113x1462, plasmid153 2338x1649, leaders 3 to 4 x longer), and T5 at 14 pt cannot be rescued even at 2.0.
  Scaled to the same page width, a 1.5x band makes the text effectively smaller anyway, so widening is not free.
- A3 (drop) gives the cleanest, most readable figures (canvas almost unchanged) but silently loses 11 to 48 % of labels, the tRNAs and short genes first. T5 at 14 pt loses 94 of 195.
- A (combined) reads well for P1 (12.5 pt / 1.15) and is the least invasive for every case, but T5 is simply too dense for this model: 6.5 pt at band 1.2 for the 8 pt default, 8.5 pt at band 1.5 for 14 pt. T5 at 6.5 pt is readable only at zoom.
- Wrong tool: records with 190+ labels (T5). Right tool: 100 to 180 labels where the overflow is 10 to 30 %.
- Unchanged look for everything that already fits (byte-identical), and no new concepts for the user: same ellipse, same leaders.

## 5. Runtime
Probe cost is dominated by the first (failing) render: today's pipeline spends 74 s on T5 at 14 pt (195 labels) when it cannot
fit; all later probes were 0.4 to 3 s (they fit or fail early), so the loop adds only 5 to 9 probes: T5 14pt total 91 s (A), 112 s (A3), 88 s (A1/A2);
tobacco 14pt 25 to 30 s; plasmid153 14pt 40 to 80 s (machine shared with other runs); P1 8 to 18 s; fitting cases +0.4 s (one evaluate, ~0.3 s for 195 labels).
A production version must not render the infeasible request first. Option: a cheap capacity pre-check (labels per side times line pitch against the ellipse height, drift-bounded) to jump
straight to the first candidate setting, plus an early exit in the dense optimiser when overlaps remain after the legacy pass.
The sweeps cost 20 to 1000 s per case (exploration only).

## 6. What production would need
- `labels/circular.py`: wrap `_prepare_label_list_today` in a loop like `render_today` (the study wrapper is the shape); `band_scale` applied after the arena rescale; `drop` already a filter on `feature_dict` ids.
  The font knob needs `_candidate_cache=None` (the study wrapper does it). The fit test is about 80 lines (pure function of the label list).
- Reporting: dropped labels and any reduced font/band must be returned to the CLI and web UI (a warning line, a list in the GUI, an entry in the session) or the figure silently differs from the request.
- Option names: an opt-in `--label_fit auto|off` would be safest (default off for 0.14.x); `auto` runs A with font floor, band cap and `drop` as three user-visible parameters.
- Risk: the criterion is the contract; (d) drift 30 deg and (c) are my choices. Today's byte output must not change for already-fitting figures (verified here).

## 7. Recommendation and Owner decisions
For 0.14.x I recommend **A (combined) with no dropping by default**: it never changes a figure that fits, changes the least on those that do not, and
returns a verifiable "fits" result; dropping stays an explicit opt-in (A3) because it silently removes content. T5-class records (about 190 labels) still need either small text or
drop; no knob in today's model makes 14 pt T5 legible.
The Owner must decide: (1) font floor (I used 6 pt; 6.5 pt T5 is barely readable; 7 or 8 would turn T5 into "drop or fail"); (2) band cap (1.5 here; it grows the canvas up to
1886x1448 and leaders to ~400 px; 1.25 would leave tobacco 14pt and plasmid 14pt unfitted); (3) whether dropping is allowed at all, and then the rule (shortest feature first is mine; gbdraw has no label importance)
and how it is reported (warning plus list in the figure legend or session metadata); (4) whether "fits" may include the drift bound (30 deg) or only overlap-freeness; (5) automatic or opt-in.

## 8. Not done / caveats
Inner labels (`--labels both`) not covered. The T5 "hypothetical protein" default-blacklist question was not investigated. A fit that fails at the floor/cap (A1 on T5, A2 on T5 14pt) is returned as the last probe with `fits: False`.
The A1 label counts include labels that became embedded at smaller fonts (not drops). Timings are wall clock on a machine shared with other agents.
