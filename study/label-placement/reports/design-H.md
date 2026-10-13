# Design H: ordered free placement on both sides of the ring

File `gbdraw/labels/proto/design_h.py` (455 lines, deterministic, numpy; imports `_Zones`, `_attach`, `_seg_box_hit` from design G, `pack_ordered`
etc. from common; no other file edited). Renders `work/H/<case>__H.*`; variants `__H_in0`, `__H_clear6`; invisible-text twins `__Hhide.*`
(option `hide=1`, evidence only). Checker `work/H/check_h.py` (`check_all.txt`), tables `work/H/table_md.txt`. 1 to 6 s per case (max 5.7 s).

## 1. The rule
- **Outside shore = design D, pulled in.** Two columns (feature on the right / left half, one cut near 12 and 6 o'clock chosen as in D by exhaustive
  search), clockwise order kept, `pack_ordered` for y, D's one-sided tangent bound and comb routing at 12/6 o'clock, D's x monotonicity rule.
  Changes against D: the wall is the drawn axis circle (radius 390 + half stroke 1.5 = 391.5 px, read from `radial_layout.axis`) plus `clear`, so the box's
  nearest corner is 3 px from the axis line (D: 24 px beyond the elbow); the leader leaves at axis + 1 px; the minimum horizontal run beside the stub is
  5 px (D 21); the stub into the text is 3 px (D 8). A label whose feature is not crowded ends right beside it with a 5 to 10 px leader.
- **Inside shore** (only if the free gap [A, B] inside the ring is wider than 3 font sizes; from `radial_layout` zones as in G: B = ring inner edge - clear,
  A = the next forbidden annulus (title disc, or a GC/skew track, or the tick band) + clear; the outermost gap only, leaders may cross only the tick
  band). Right column = features on the right half, ordered top to bottom by clockwise angle; left column mirrored. Text horizontal, right edge of the box on
  the circle B (largest x with all corners inside B), extending toward the centre; the two columns keep apart (right boxes at x >= 1.5, left at x <= -1.5).
  y by `pack_ordered` toward the stub's y, clamped to the per-item feasible y run (box fits the gap, slides toward the centre in 2 px steps up to 40 px to clear
  tick labels, the straight leader stub -> nearest edge point (G's attach rule, 4 px off corners) touches no tick label and stays outside radius A).
  Leader = anchor -> inner stub on the ring's inner edge -> attach point. Ring separates the shores: inside and outside leaders cannot meet.
- **Which labels go inside.** All outside first. Loop: for every outside label whose inside run exists (in order of leader length, all of them), build the inside
  layout with it added and re-lay the outside with it removed (cut thresholds fixed); a trial must pass the exact inside validity test (no box overlap with gap 3 px,
  no leader within 1.5 px of another box or tick label, no leader crossing, no label running toward the axis). Accept the trial with the lowest total
  leader length (all polylines including the hidden anchor -> stub piece) if it is lower than now; repeat until none lowers it (<= 400 rounds), then
  re-search the outside cut points. Deterministic; the final rule is "greedy best single move, steepest descent, exact feasibility".
- **No drop, no crossing, no overlap:** the outside shore never drops, it grows (column height); the inside has a hard capacity: a move that is infeasible is not made.
- **LP-D09 (no label runs back toward the ring)**: outside by construction (right labels start at the leader and run right, left labels the mirror, combs
  start at x >= 6); inside enforced in the feasibility mask and again exactly with the metric's own test (`_toward_axis`). Before the rule 4 inside labels violated
  it on p1_open and 2 on plasmid37_open; they are now excluded (total leader +2 %).

## 2. Parameters
clear 3 (= ctx.gap) | inside 1 | run 5 | hs 3 | stub 1 | xoff 0 | cutmax 40 deg | cutin 0 | cand 200 | slide 40 | lpad 1.5 | rounds 400.

## 3. Metrics, all 23 cases (sumL = `leader_len_sum`, px, whole polylines; G rows include its dropped labels, so G's sumL is understated)
`H ... toward` = in/out, dropped, text overlaps, leader crossings, leader-through-text, `labels_toward_axis`. All six H counts are 0 on every row (23 rows) and
in `_in0` / `_clear6`.

| case | cur sumL / canvas | D sumL / canvas | G in/out drop cr sumL / canvas | H in/out drop ov cr tt toward | H sumL | H canvas |
|---|---|---|---|---|---|---|
| p1_f14 | 12731 / 1706x990 | 12526 / 1624x1179 | 4/81 26 0 3798 / 1551x891 | 0/111 0 0 0 0 0 | 10564 | 1574x1224 |
| pbr322_f14 | 3050 / 1256x910 | 3936 / 1250x913 | 2/31 1 0 970 / 1145x861 | 0/34 0 0 0 0 0 | 1507 | 1170x874 |
| mito_f14 | 2898 / 1476x906 | 3752 / 1490x899 | 0/33 0 0 1351 / 1380x864 | 0/33 0 0 0 0 0 | 1738 | 1407x854 |
| mg1655_rna_f8 | 11789 / 1200x977 | 10174 / 1136x1001 | 5/59 44 0 2527 / 1067x840 | 0/108 0 0 0 0 0 | 9352 | 1096x1025 |
| tobacco_f14 | 37840 / 1634x988 | 19830 / 1555x1417 | 2/82 56 0 4291 / 1494x923 | 0/140 0 0 0 0 0 | 20019 | 1524x1436 |
| t5_f14 | 65083 / 1821x992 | 66286 / 1610x1965 | 0/80 115 0 7258 / 1817x957 | 0/195 0 0 0 0 0 | 66105 | 1566x1985 |
| t5_f8 | 67238 / 1444x979 | 40309 / 1359x1344 | 1/110 83 0 6910 / 1251x885 | 0/194 0 0 0 0 0 | 36908 | 1318x1354 |
| plasmid153_f14 | 16536 / 1714x994 | 45718 / 1650x1857 | 0/83 96 0 6605 / 1710x918 | 0/179 0 0 0 0 0 | 44260 | 1619x1857 |
| plasmid153_f8 | 18516 / 1436x983 | 16872 / 1364x1190 | 2/116 61 0 4623 / 1380x908 | 0/179 0 0 0 0 0 | 15466 | 1328x1201 |
| lambda_f14 | 8171 / 1567x910 | 8256 / 1546x956 | 0/63 8 0 3073 / 1415x892 | 0/71 0 0 0 0 0 | 4845 | 1445x941 |
| n4_f14 | 8544 / 1402x994 | 6311 / 1322x1025 | 1/57 11 0 2014 / 1278x885 | 0/69 0 0 0 0 0 | 4476 | 1282x1036 |
| epifagus_f14 | 5878 / 1509x990 | 4796 / 1459x910 | 0/44 9 0 1892 / 1393x868 | 0/53 0 0 0 0 0 | 3138 | 1416x914 |
| plasmid37_f14 | 5411 / 1648x916 | 6023 / 1666x921 | 0/45 0 0 2575 / 1558x926 | 0/45 0 0 0 0 0 | 3180 | 1583x876 |
| p1_f8 | 10930 / 1441x982 | 8843 / 1361x961 | 10/92 5 0 3627 / 1299x873 | 0/107 0 0 0 0 0 | 5584 | 1320x983 |
| tobacco_f8 | 17043 / 1399x983 | 12281 / 1334x1035 | 4/103 33 0 3935 / 1274x896 | 0/140 0 0 0 0 0 | 8795 | 1295x1064 |
| pbr322_f10 | 2559 / 1217x904 | 3870 / 1210x900 | 2/32 0 0 976 / 1106x854 | 0/34 0 0 0 0 0 | 1429 | 1129x854 |
| gal_hmmt_basic | 2100 / 1491x984 | 2917 / 1528x937 | 0/32 0 0 1309 / 1463x949 | 0/32 0 0 0 0 0 | 1525 | 1479x927 |
| gal_hmmt_atskew | 1541 / 1281x954 | 2018 / 1293x906 | 0/24 0 0 962 / 1227x905 | 0/24 0 0 0 0 0 | 1032 | 1247x893 |
| pbr322_open_f14 | 3050 / 1256x910 | 3936 / 1250x913 | 3/31 0 0 1130 / 1145x868 | 0/34 0 0 0 0 0 | 1507 | 1170x874 |
| lambda_open_f14 | 8171 / 1499x910 | 8256 / 1478x956 | 9/59 3 0 2973 / 1347x879 | 10/61 0 0 0 0 0 | 3864 | 1377x898 |
| p1_open_f14 | 12731 / 1706x990 | 12526 / 1624x1179 | 25/77 9 0 4622 / 1551x904 | 28/83 0 0 0 0 0 | 6699 | 1574x1081 |
| mito_open_f14 | 2898 / 1415x906 | 3752 / 1430x899 | 3/30 0 0 1252 / 1320x864 | 5/28 0 0 0 0 0 | 1549 | 1346x854 |
| plasmid37_open_f14 | 5411 / 1580x916 | 6023 / 1598x921 | 2/43 0 0 2488 / 1490x904 | 13/32 0 0 0 0 0 | 2719 | 1515x876 |

Headline (leader mean / label-to-feature mean px): p1_open_f14 H 60/53 (cur 115/108, D 113/97, G 45/43 with 9 dropped); lambda_open H 54/49 (G 44/41, 3 dropped);
pbr322_open H 44/38 (G 33/31); p1_f14 H 95/81 (D 113/97, G 45/42, 26 dropped); t5_f8 H 190/161 (D 208/174, G 62/57, 83 dropped).
Variants on p1_open / lambda_open / pbr322_open (sumL, in/out, canvas): `inside=0`: 10564 0/111 1574x1224; 4845 0/71 1377x941; 1507 0/34 1170x874.
`clear=6`: 6886 27/84 1574x1086; 4062 10/61 1381x892; 1691 0/34 1170x870. Both variants 0 dropped / crossings / overlaps / toward.
Own check (`check_h.py`, from the invisible-text twin of every render, pixel mask of everything else drawn: axis, features, tracks, ticks, tick labels, title, legend): label boxes
on other content 0 and leaders on other content 0 on all 23 cases (inside leaders excepted between the ring and the tick band, where they must cross it). Smallest
leader-to-other-label-box distance 2.0 px. Where the canvas grows: the outside columns (p1_f14 1,140 to 1,170 px tall columns, t5_f14 1,950 px, plasmid153_f14
1,857 px); against D the canvas is 3 to 8 % narrower and about equal in height.

## 4. Assessment (from the figures, 1:1 crops of p1_open and p1_f14)
- **Close to G on open-centre records, not equal.** p1_open: 28 of 111 labels inside, sides tight, total leader 6,699 px against D 12,526 (-47 %), today 12,731, canvas
  1574x1081 against D 1624x1179. G's 4,622 px had 9 labels dropped; H places all 111. Leaders are 30 to 40 % longer than G's per label (60 against 45 px on P1) because labels
  stand in columns, not at the nearest free spot, and the hidden anchor -> stub piece is counted. Lambda and mito/plasmid37 read like G; pbr322 stays all outside (G
  moves 3 inside and is 25 % shorter, 1,130 against 1,507).
- **With default inner tracks H is still D.** The free gap is 19 px wide, so no label goes inside (p1_f14, t5_*, tobacco, plasmid153, lambda_f14 ...). There the gain
  is only the tighter outside: p1_f14 -16 % total leader (10,564 against 12,526), lambda_f14 -40 %, p1_f8 -37 %, n4 -29 %; dense records change little (t5_f14 66,105 against 66,286, plasmid153_f14
  44,260 against 45,718: height-limited, canvas still 1,985 and 1,857 px tall).
- **Combs at 12 and 6 o'clock:** unchanged in kind and slightly smaller. They are still parallel vertical leaders (P1: about 8 labels above the ring on the right, 7 on the left) up to
  about 300 px long. On open records inside labels take some of them (UpfO..., Lxr, Coi stay), but only short labels fit near the top because the two inside columns may not meet.
- Leader clearance is 2 px at the closest (D's rules, not a free parameter); inside leaders pass the tick band and sometimes run 60 to 100 px across the interior (the
  "tail protein" group on P1). Tick labels leave windows only: about half of the ring's inside is blocked for leaders at the tick-label band.
- Right tool: complete labelling with the centre free (plasmids, phages, organelle with `--no-gc --no-skew`). Not a fix for very dense records at large font.

## 5. What production would need
Same module as D (`prepare_label_list` placement mode) plus the zone model of G (`radial_layout` annuli, tick-label boxes, axis radius). Shared with G: tick-label boxes must be
exported by the tick renderer instead of recomputed. The shape is small: outside = D (150 lines), inside + move search about 250 lines. Risk: inside validity is
quadratic in the number of inside labels (fine to 200); steepest descent is greedy and non-monotone in tuning (cutin, cand showed 7,000 to 8,800 on P1 before the
"columns never meet" rule); G's `clear` from `r_out` alone puts boxes on the drawn axis line (circle stroke reaches r = 391.5): G needs the same axis-radius fix.

## 6. Not done
No inside labels with default tracks (gap too thin; a second gap between GC and the title is not used because leaders may not cross tracks). No bent leaders around tick
labels. No top/bottom inside row (columns never meet, so wide labels near 12 o'clock stay outside). No other track types, no `--labels both`. Outside leader clearance is
D's (2 px minimum), not `clear`. The hide option and `axis` note keys are evidence-only.
