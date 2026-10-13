# Figure review by the campaign session (e334f637), 2026-10-12

Looked at: page/review/<case>.png (whole maps, map centres aligned, one scale) and <case>--crop1.png (1:1) for
Current, B, D, E, R. Renders: page/renders/<case>__<design>.* (all by page/run_extra.py with the final design files).
Viewed so far: p1_f14 (+crop), pbr322_f14, mito_f14, tobacco_f14, t5_f8, t5_f14, plasmid153_f14, mg1655_rna_f8 (+crop),
lambda_f14, n4_f14, plasmid37_f14, E/D variants on pbr322_f14 (whole + bottom crop). Not yet viewed: epifagus_f14,
p1_f8, tobacco_f8, plasmid153_f8, pbr322_f10, gal_hmmt_basic, gal_hmmt_atskew (controls; view with A and C added).

## What each design looks like
- Current: fine on easy maps (pBR322, mtDNA, pVNTG4, lambda, N4). Fails by overprint at the end of the packing
  (P1 14 pt upper left; tobacco 14 pt upper left; plasmid153 14 pt everywhere) and, on T5 at the Auto 8 pt, by
  leaders sweeping across the whole map.
- B (tiers): easy maps look like today. Hard maps: structured, nothing dropped, but the leaders of the outer column
  form grey horizontal bundles through the inner column (P1 left, tobacco right, T5 right: a wall), and labels that
  overshoot 12 o'clock stack as a tall staircase above the map (N4 gp1..gp11, T5 8 pt upper right). Width grows.
- D (columns hugging the ring): closest to today's look, clean on every case, nothing dropped. Cost is height and
  the "comb" of parallel vertical leaders above and below the map (P1 top, T5, plasmid153). On small maps the
  top labels are pulled sideways with long horizontal leaders (pBR322, pVNTG4 top). D_flat (strictly vertical
  columns, like CLC "Stacked") gives very long leaders (pBR322: 530 px) and looks like a table beside the map.
- E (cluster lists): easy maps look like today with small flush lists where labels cluster. pBR322: lists of
  2 to 4 sites, like the Owner's plasmid references. MG1655 tRNA/rRNA: the best case, each rRNA operon becomes
  one bracketed list with one leader. Dense maps: default drops labels (P1 14 pt 12, tobacco 14 pt 22,
  T5 8 pt 59, plasmid153 14 pt 64) and stays compact; grow mode places all but gets as tall as D or taller.
  Leader styles on pBR322 bottom cluster: fan = one line per row (pGLO reference); bracket = one leader and a
  bracket over the sites (MCS reference), quieter, but the row to site mapping is only the order.
- R (radial, shipped): never overprints or crosses, but the canvas grows a lot (P1 14 pt 1735x1580, tobacco 14 pt
  2427x2206, T5 14 pt 2288x2530) and the map becomes small; text at every angle is harder to read.
- Extreme cases (T5 14 pt, plasmid153 14 pt): no design gives a figure one would publish; the practical answers are
  the Auto font (8 pt), fewer labels (the default blacklist, OV-470, would remove 90 of T5's 195), or dropping.

## Design A numbers on the extra cases (reports pending)
A's fit test is stricter than "no overprint": it also rejects leaders through text and labels more than 30 degrees
from their feature. So A changes maps that have no overprint today: lambda 14 pt (A 13.5 pt + band 1.05; A1 9 pt;
A3 drops 23 of 71), pBR322 14 pt (13.5 pt), MG1655 RNA (A 6.5 pt; A1 6 pt; A3 drops 30). Gallery HmmtDNA basic and
ATskew, Epifagus and pBR322 10 pt are unchanged by A.

## Design A (viewed 00:59: page/review_a/p1_f14.png and crop)
- A (font and band together) on P1 14 pt: 12.5 pt, band 1.15, all 111 labels, 1773x1130. Same look as today, clean
  at the upper left; leaders are longer (labels sit further out). The least surprising fix.
- A1 (font only): 8 pt on a 14 pt request: the text is too small next to the map; A2 (band only): 14 pt kept but
  1928x1226 and long leaders; A3 (drop only): compact (1636x992), 12 labels gone.

## Design G (viewed 01:31: page/review_g/ p1_open_f14, lambda_open_f14 crops, pbr322_open_f14 whole)
- Open centre (no GC tracks): the Geneious look. Labels sit right beside their features with short straight
  leaders; in crowded stretches some go inside the ring (lambda: "restriction alleviation ral", "CIII
  anti-termination", "excisionase"; P1: "PmgS-like head morphogenesis", "glycosyltransferase"); canvas 6 to 10 %
  smaller than today (pBR322 1145x868 vs 1256x910; P1 1551x943 vs 1706x990). pBR322: 3 of 34 inside (ClaI, HindIII,
  BstZ17I), the rest hug the ring.
- Costs seen: P1 open drops 4 of 111 and has 14 crossings in the dense top; with gbdraw's default inner tracks the
  inside is closed and P1 14 pt drops 19 of 111. Dense records drop 20 to 40 % (tobacco 14 pt 43, T5 8 pt 46,
  plasmid153 14 pt 77) and take 30 to 70 s. It has no fallback that grows outward.
- outside only (inside=0): P1 open drops 21, so the inside carries much of the load. clear=6: a little more air,
  5 dropped on P1 open, canvas +40 px.
- merge=1 writes "text (x2)"; only 2 equal runs on P1 and lambda. color=1 works; light-blue text is hard to read.

## Design C against G (viewed 01:36: page/review_cg/p1_f14.png and crop)
- C (free, outside only, bounded straight leader): tidy fans with short leaders, no ellipse; P1 14 pt 98 of 111 placed
  (13 dropped), 15 crossings, 1577x1167 (overflow stacks above and below the map). Clean cases get worse than
  today: lambda 14 pt drops 3 with 8 crossings, MG1655 RNA drops 6.
- G with default tracks on P1 14 pt: labels hug the ring, 4 inside (Coi, ImcB between ring and tick labels), 19
  dropped, 15 crossings, one long horizontal leader ("Lpa-like transcriptional activator"). Canvas 1551x909.
- Both free designs fail the same way on dense records: they drop labels (C 13, G 19 on P1 14 pt; C 74, G 46 on
  T5 8 pt). Neither grows outward when a stretch is full. A production version needs that (LP-D05 says no drop).

## C and G with crossings forbidden (fallback=0, LP-D07; viewed page/review_cg/lambda_open_f14.png)
- lambda, open centre, 71 labels: today clean (fan with elbows, 1499x910); C 62 placed / 9 dropped, 1465x995; G 68
  placed (9 inside) / 3 dropped, 1347x879, the most compact; G with crossings allowed 69 placed, 1 crossing; D all
  71, no crossing, 1478x956 (columns with a comb at the top).
- Numbers for all cases: page/renders/<case>__{C,G}.metrics.json (strict); *_cross keeps the crossing defaults on
  p1_f14 and the three open cases. Strict drops more: P1 14 pt C 24, G 26; P1 open G 9; T5 8 pt C 87, G 83.
