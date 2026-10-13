# Part 2: where today's Circular label placement breaks (2026-10-12)

Setup: default canvas, `--separate_strands --track_type tuckin --labels out --legend right`, placement horizontal, worktree `label-study`
(design `cur`). Fonts are `--label_font_size` values (SVG px). Auto font: records of 50 kb or more get 8, shorter ones 14; rows marked "Auto".
All renders are in `work/sweep/` (`<rec>_f<font>.svg/.png/.metrics.json/.labels.json`; rec = p1 t5 n4 tobacco mito pl37 pl153).
Figures: `figures/part2/<case>__cur_full.png` (flattened, width <= 1600) for every case, `<case>__cur_crop.png` (1:1 crop) for the cases
named below. "viewed" = I looked at the full figure (and the crop for failing cases). Cases not marked viewed were judged from metrics only.
Overprint = `text_overprint_pairs`; "through text" = `leader_through_text` (a leader touching another label's text).

### Phage P1 94.8 kb (extra)

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 8 | 107+7 | 0 | 0 | 0 | 0 | 142 | 130 | 1441x982 | 3 | clean (Auto) |
| 9 | 109+5 | 0 | 0 | 0 | 0 | 149 | 135 | 1485x983 | 3 | clean; 1 leader touches text |
| 10 | 109+5 | 0 | 0 | 0 | 0 | 159 | 140 | 1529x985 | 3 | clean; 1 leader touches text |
| 11 | 111+3 | 0 | 0 | 0 | 0 | 172 | 149 | 1573x987 | 2 | clean; 1 leader touches text |
| 12 | 111+3 | 0 | 0 | 0 | 0 | 184 | 165 | 1615x988 | 3 | clean (viewed); dense fan upper left |
| 13 | 111+3 | 4 | 4 | 7 | 1 | 519 | 487 | 1660x989 | 8 | FAIL: overprint at 12 o clock, one far label (Lpa-like, 487 px) |
| 14 | 111+3 | 20 | 20 | 21 | 1 | 467 | 438 | 1706x990 | 11 | FAIL: overprint at 12 o clock, 8 leaders through text |

### Phage T5 121.8 kb (extra)

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 6 | 193+2 | 0 | 0 | 12 | 0 | 298 | 276 | 1339x976 | 5 | no overprint but 12 crossings, 47 leaders touch text; 6 px text; labels up to 276 px away |
| 7 | 193+2 | 0 | 0 | 815 | 37 | 562 | 497 | 1396x978 | 5 | FAIL: 815 crossings, labels 497 px away, leaders over map |
| 8 | 194+1 | 0 | 0 | 2085 | 70 | 699 | 629 | 1444x979 | 7 | FAIL (Auto): 2085 crossings, 101 through text, labels 629 px away (viewed) |
| 10 | 194+1 | 63 | 71 | 3963 | 89 | 838 | 780 | 1637x982 | 24 | FAIL: overprint 63 + crossings |
| 12 | 195+0 | 127 | 146 | 3736 | 78 | 907 | 836 | 1687x986 | 46 | FAIL: overprint 127 |
| 14 | 195+0 | 170 | 189 | 2532 | 63 | 899 | 827 | 1821x992 | 57 | FAIL: overprint 170, pile at right (viewed) |

### Phage N4 70.2 kb (extra)

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 8 | 61+11 | 0 | 0 | 0 | 0 | 134 | 123 | 1227x983 | 1 | clean (Auto, viewed) |
| 10 | 65+7 | 0 | 0 | 0 | 0 | 164 | 150 | 1312x986 | 1 | clean |
| 12 | 69+3 | 0 | 0 | 0 | 0 | 197 | 180 | 1358x990 | 1 | clean |
| 14 | 69+3 | 0 | 0 | 0 | 0 | 230 | 210 | 1402x994 | 1 | clean |
| 16 | 69+3 | 0 | 0 | 0 | 0 | 262 | 238 | 1459x997 | 1 | clean (viewed) |

### Tobacco plastid NC_001879 155.9 kb

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 8 | 140+3 | 0 | 0 | 0 | 0 | 226 | 194 | 1399x983 | 1 | clean (Auto) |
| 10 | 140+3 | 0 | 0 | 0 | 0 | 271 | 237 | 1512x987 | 1 | clean (viewed) |
| 12 | 140+3 | 9 | 11 | 185 | 3 | 596 | 563 | 1538x983 | 11 | FAIL: overprint at 11-12 o clock, 185 crossings |
| 14 | 140+3 | 50 | 53 | 1062 | 50 | 705 | 626 | 1634x988 | 19 | FAIL: overprint 50, left side pile (viewed) |

### Human mtDNA 16.6 kb

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 14 | 33+4 | 0 | 0 | 0 | 0 | 127 | 122 | 1476x906 | 2 | clean (Auto, viewed) |
| 18 | 34+3 | 0 | 0 | 2 | 0 | 231 | 180 | 1600x913 | 2 | clean visually; 2 crossings, 1 through text (bottom) |
| 22 | 35+2 | 0 | 0 | 0 | 0 | 165 | 121 | 1715x918 | 2 | clean metrics except 1 leader touches text (not viewed) |
| 26 | 35+2 | 0 | 0 | 2 | 0 | 265 | 185 | 1829x909 | 17 | clean visually; 2 crossings at bottom (viewed) |

### Plasmid pVNTG4 37.1 kb

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 14 | 45+0 | 0 | 0 | 1 | 0 | 248 | 173 | 1648x916 | 2 | clean (Auto, viewed); 1 crossing |
| 18 | 45+0 | 0 | 0 | 0 | 0 | 204 | 143 | 1810x923 | 38 | clean (viewed); 2 leaders touch text |
| 22 | 45+0 | 2 | 2 | 6 | 0 | 247 | 194 | 1972x931 | 41 | FAIL: overprint at 6 o clock, 6 crossings; 41 s |

### Plasmid p1 153.4 kb

| font | labels ext+emb | overprint | overlap | crossings | leaders over feat. | leader max px | label-feature max px | canvas | s | verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 8 | 179+0 | 0 | 0 | 0 | 0 | 207 | 178 | 1436x983 | 1 | clean (Auto, viewed) |
| 10 | 179+0 | 63 | 64 | 410 | 1 | 594 | 524 | 1577x982 | 28 | FAIL: overprint 63 upper left, 410 crossings, label 524 px away |
| 12 | 179+0 | 121 | 123 | 1 | 0 | 128 | 116 | 1614x990 | 8 | FAIL: overprint 121 both sides |
| 14 | 179+0 | 143 | 144 | 17 | 0 | 177 | 160 | 1714x994 | 14 | FAIL: overprint 143 (not viewed beyond metrics) |
## Breaking point per record (largest tested font that is clean, and what fails first)

- **P1 (94.8 kb)**: clean up to 12 pt. 13 pt fails first: 4 overprint pairs at 12 o'clock ("terminase / C1.100 / HNH endonuclease") and one
  label (Lpa-like transcriptional activator) pushed 487 px across to the right. 14 pt: 20 overprint pairs, same place. Auto (8) is clean.
  Crops: `p1_f13__cur_crop.png`, `p1_f14__cur_crop.png`, clean: `p1_f12__cur_crop.png`, `p1_f8__cur_crop.png`.
- **T5 (121.8 kb)**: no clean tested size. 6 pt (6 px text, hard to read) has no overprint but 12 leader crossings; 7 pt already has 815 crossings
  and labels up to 497 px from their feature; Auto 8 pt has 2085 crossings, leaders sweeping across the whole map and the GC rings (viewed);
  overprint starts at 10 pt. Failing at the Auto size. Crops: `t5_f6/f8/f10/f12/f14/f7__cur_crop.png`. See the blacklist finding: with the
  documented default blacklist the 8 pt case is clean (104 labels).
- **N4 (70.2 kb)**: clean at every tested size, 8 to 16 pt (69 labels at most). Best case for the model. Crops: `n4_f16`, `n4_f8`.
- **Tobacco plastid (155.9 kb)**: clean up to 10 pt. 12 pt fails first: overprint at 11-12 o'clock (9 pairs) and 185 crossings from the
  inverted-repeat genes whose labels are pushed round the left side; 14 pt: 50 pairs, labels 626 px from their feature. Auto (8) clean.
  Crops: `tobacco_f12`, `tobacco_f14`, clean `tobacco_f10`.
- **Human mtDNA (16.6 kb)**: no overprint up to 26 pt (largest tested). Only minor leader crossings (2 at 18 and 26 pt) at the bottom, where
  four labels share 6 o'clock; at 26 pt text is large but nothing is cut (viewed). Breaking point not reached.
- **Plasmid pVNTG4 (37.1 kb)**: clean up to 18 pt. 22 pt fails first: overprint of the two "hypothetical protein" labels at 6 o'clock, 6 crossings;
  placement takes 41 s for 45 labels at 18 and 22 pt (38 s, 41 s) versus 2 s at 14 pt. Crops: `pl37_f22`, `pl37_f18`, `pl37_f14`.
- **Plasmid p1 (153.4 kb)**: clean at 8 (Auto) only. 10 pt is the worst case of the whole sweep for its size: 63 overprint pairs, 410 crossings, a
  label 524 px from its feature; 12 and 14 pt overprint on both sides (121, 143 pairs). Not monotonic: 10 pt has more crossings than 12 pt.
  Crops: `pl153_f10`, `pl153_f12`, `pl153_f14`, clean `pl153_f8`.

## Failure kinds and the case that shows each best

- **Overprint** (text over text): `p1_f14` (top, 12 o'clock), `tobacco_f14` (upper left), `pl153_f12` (both sides, whole rows), `t5_f14`.
  It starts where the sequential packer runs out of ring; mostly in a heap at the end of the packing (upper left, 10-12 o'clock).
- **Leaders crossing / sweeping over the map**: `t5_f8` (2085 crossings, leaders cross the GC rings and the whole disc), `pl153_f10` (410).
- **Labels far from the feature**: `t5_f8` (629 px), `pl153_f10` (524), `p1_f13` (487, one label), `tobacco_f14` (626).
- **Leaders over the feature ring** (`leaders_over_features`): T5 at 7 to 14 pt (37 to 89), tobacco 12/14 pt (3, 50).
- **Unreadably small text**: `t5_f6` (6 px; still not clean). 8 px (Auto for 50 kb+ records) is small but readable at 1:1 (`p1_f8`).
- **Labels off the canvas**: none in the horizontal placement (the canvas grows with the labels: widths 1227 to 1972 px). It grows only a
  little: heights stay about 980 px because the ellipse is fixed.
- **Labels dropped silently**: none seen (every case draws all labels; `proto` is empty). Nothing is dropped for lack of room.
- Pure slowness: T5 at 12/14 pt 46 to 57 s, pl37 18/22 pt about 40 s, mito 26 pt 17 s.

## Gallery examples, run unchanged through `run_case.py` (`--session` works with the appended `-o/-f/--overwrite`; no error)

| case | labels ext+emb | font | overprint | crossings | through text | over features | leader max | label-feature max | canvas | note |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| gallery_mito (HmmtDNA_basic_circular) | 32+5 | 14 | 0 | 0 | 0 | 0 | 107 | 99 | 1491x984 | horizontal, middle track, outer; clean (viewed) |
| gallery_tobacco (tobacco-chloroplast) | 118+15 | 8 | 0 | 0 | 0 | 68 | 106 | 93 | 1374x970 | NOT horizontal: session uses `--label_placement radial`, `--labels both`, offsets 0.90/0.975; clean (viewed) |
| gallery_ATskew (HmmtDNA_ATskew, extra) | 24+13 | 20 | 0 | 0 | 0 | 0 | 111 | 94 | 1281x954 | horizontal, outer, 20 pt; not viewed |

Crops: `gallery_mito__cur_crop.png`, `gallery_tobacco__cur_crop.png`. Tobacco's 68 "leaders over features" are the inner labels' leaders
(scope both), not a fault. The Gallery sessions carry `blacklist_keywords: []`. I did not diff the renders against the stored gallery SVGs.

## `--label_placement radial` (design R), default canvas

| case | labels | font | overprint | crossings | leader max | canvas (today's horizontal canvas) | verdict |
|---|---:|---:|---:|---:|---:|---|---|
| p1_f14 | 111 | 14 | 0 | 0 | 145 | 1735x1580 (1706x990) | clean; but the left-most label ("TciB-like inhibition...") is cut at the canvas edge by about one letter (viewed 1:1, `p1_f14__R_crop.png`) |
| t5_f8 | 194 | 8 | 0 | 0 | 273 | 1627x1668 (1444x979) | clean, spokes up to 270 px; viewed |
| t5_f14 | 195 | 14 | 0 | 0 | 509 | 2288x2530 (1821x992) | clean but the map is about a third of the width; "collar tail protein..." touches the left edge (`t5_f14__R_crop.png`) |
| tobacco_f8 | 140 | 8 | 0 | 0 | 261 | 1717x1517 (1399x983) | clean (metrics; not viewed) |
| mito_f14 | 33 | 14 | 0 | 0 | 127 | 1488x1185 (1476x906) | clean (metrics; not viewed) |
| plasmid153_f14 | 179 | 14 | 0 | 0 | 268 | 1956x1960 (1714x994) | clean (metrics; not viewed) |

Radial never overprints or crosses, at the price of a canvas up to 2.5 times taller and a small map. Full images: `*__R_full.png`.

## Side finding: default label blacklist is never applied (not fixed)

- Steps: `gbdraw circular --gbk NC_005859.1.gb --labels out` (no `--label_blacklist`). Expected: labels containing "hypothetical",
  "uncharacterized", "putative", "unknown" are dropped (`gbdraw/data/config.toml:227`). Actual: T5 draws 195 labels, 79 with "hypothetical"
  and 11 with "putative"; plasmid p1 105 of 179, pVNTG4 19 of 45, tobacco 14 of 143, N4 3 of 72.
- Cause: `gbdraw/cli_utils/common.py:504-507` gives `--label_blacklist` `default=""`; `gbdraw/circular.py:1107-1112` (same in
  `gbdraw/linear.py:1432-1437`) tests `label_blacklist is not None`, which is true for `""`, and replaces the config list with `[]`
  (`"".split(",")` yields nothing). `filtering.py:572` then has an empty list to apply. The default list in config.toml is never reached from
  the CLI. Web sessions write `blacklist_keywords: []` explicitly (both Gallery circular sessions do).
- P1 behaves the same; it simply has no label containing those four words (0 of 114), so it looks unaffected.
- Effect: T5 at 8 pt with `--label_blacklist "hypothetical,uncharacterized,putative,unknown"` gives 104+1 labels, 0 overprint, 0 crossings, leader
  max 139 px, 1 s (`work/sweep/t5_f8_bl.*`) against 2085 crossings and 7 s today. Fixing it is a visible change for every genome
  with such products, including the tracked reference `circular_with_labels.svg` (it contains "hypothetical protein" labels). Owner decision.

## Tracked figures that use Circular outer labels with the horizontal placement

(a) `tests/reference_outputs/` (16 files; commands in `tests/test_output_comparison.py` TEST_CASES, all on MjeNMV.gb or AP027078.gb):
- Circular, 10 files. With labels: `circular_with_labels` (`--separate_strands --labels`: horizontal, outer, 8 pt, 108 external + 3 embedded) and
  `circular_radial_labels` (`--labels both --label_placement radial`, 111 external, not horizontal). The other 8 Circular files pass no
  `--labels` (default `none`): basic, tuckin, middle, spreadout, no_gc, no_skew, no_gc_skew, hepatoplasma_repeat_underlay (0 label elements in the SVG).
- Linear, 6 files: none passes `--show_labels`; no label elements.
- Result: **1 of 16** (`circular_with_labels`) would change with a new default horizontal placement.

(b) Gallery (`gbdraw/web/gallery/examples.json`, 10 examples): 4 Circular, 6 Linear.
- Circular horizontal with outer labels: **2**: HmmtDNA_basic_circular (32 external + 5 embedded), HmmtDNA_ATskew (24 + 13, 20 pt).
- Circular radial: tobacco-chloroplast (118 external incl. inner + 15 embedded; `radial`, `both`). Circular without labels: Vnig_TUMSAT-TG-2018.
- Linear (not affected): lambda_basic_linear and BGC0000708-BGC0000713 show labels; the other four show none.
- Counts for the two horizontal examples and tobacco are from my renders of the sessions (metrics above), not from the stored SVGs.
