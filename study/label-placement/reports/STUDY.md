# Circular label placement: design study (2026-10-12)

Study only. No production change, no pull request. Prototypes live in the worktree
`/home/kawato/gbdraw-work/.worktrees/label-study` (local branch `study/label-placement`, from origin/dev af657280;
never pushed). Line numbers are those of origin/dev 5beab569.
The Owner's page: see STATE.md for the link. Figures: `figures/` (index in `figures/README.md`). Renders with metrics:
`page/renders/`, `work/A/`. Reports of the prototype authors: `logs/design-{A,B,C,D,E,G}.md`,
`logs/part2-breaking-points.md`, `logs/survey-tools.md`. What the figures look like: `logs/review-figures.md`.

## 1. How placement works today

Files: `gbdraw/labels/circular.py` (5,046 lines), `circular_candidates.py`, `circular_radial.py`, `filtering.py`;
drawers under `gbdraw/render/`; canvas in `gbdraw/diagrams/circular/assemble.py`.

- Which labels: `build_circular_label_candidates` (circular_candidates.py:45) takes the first qualifier of the
  priority list as text, applies whitelist or blacklist, and measures the text once. gbdraw has **no ranking of
  labels by importance** and never drops a label for lack of room.
- Embedded labels: text that fits inside its feature is drawn on the feature along the arc. Everything else is an
  outer label (`--labels out`); with `--labels both` minus-strand features get inner labels.
- Default model (`--label_placement horizontal`): **all outer labels sit on one ellipse** around the map
  (`rearrange_labels_fc`; tuckin: 1.15/1.10 x R for short records, 1.175/1.20 for long; R = 390 px). Text is
  horizontal, starts at the ellipse on the right half and ends there on the left. Each label has an elbowed leader.
  **Only the angle on the ellipse is free.**
- Choosing the angle. Fewer than 50 labels (`LEGACY_PLACEMENT_LABEL_THRESHOLD`, circular.py:45):
  `place_labels_on_arc_fc` then `improved_label_placement_fc` (shift limits 25 to 90 degrees). 50 or more: the
  sequential packer `_legacy_place_labels_on_arc_fc` (circular.py:2364-2430) is tried first (call at 4289); if its
  result overlaps, the newer pipeline also runs and the better score wins.
- Collision test: axis-aligned boxes with a minimum gap of 3 px (`--circular_label_spacing`).
- The canvas is not a limit: it is composed from the label bounds afterwards, so wide labels widen the figure.
- **There is no capacity check and no fallback.** Every optimiser ends with "keep the best scored state", so
  overlaps that remain are drawn. Overprint is an unhandled case, not a deliberate choice.
- A second model already ships: `--label_placement radial` (circular_radial.py: `_pack_angles` at 130, `_place_side`
  at 371): spoke text at the feature's angle, order kept, the radius grows until nothing collides; it raises
  (441, 458) instead of overprinting.

### Why phage P1 at 14 pt overprints (OV-456)

111 outer labels on an ellipse of 458 x 468 px. One horizontal label needs about one line (13.2 px + gap) of height
and each half of the ellipse has 2 x 468 px, so the ellipse holds roughly 90 to 110 labels at 14 pt. The sequential
packer runs out of ellipse on the left half. Then a plain bug decides the picture: `if new_angle > 269:
new_angle = 269.0` (circular.py:2416) clamps every later label to 269 degrees, and the overlap loop compares a label
only with the previous one, so the last 12 labels alternate between 269.0 and 270.01 degrees and label k+2 lands on
label k. The newer pipeline does not rescue it; the winner keeps 20 overprinting pairs at the upper left.

## 2. Where today's placement stops working (logs/part2-breaking-points.md)

All with `--separate_strands --track_type tuckin --labels out --legend right`.

| record | outer labels | largest font with no overprint | at the Auto font |
| --- | ---: | --- | --- |
| human mtDNA 16.6 kb | 33 | 26 pt | clean |
| pVNTG4 plasmid 37 kb | 45 | 18 pt | clean |
| phage N4 70 kb | 69 | 16 pt | clean |
| phage lambda 48.5 kb | 71 | 14 pt checked | clean |
| phage P1 94.8 kb | 107 to 111 | 12 pt | clean (8 pt) |
| tobacco plastid 156 kb | 140 | 10 pt | clean (8 pt) |
| plasmid p1 153 kb | 179 | 8 pt | clean (8 pt) |
| phage T5 121.8 kb | 194 | none | no overprint, but 2,085 leader crossings and labels up to 629 px from their feature (OV-472) |

- The limit is the number of labels against the height of the ellipse, not the record length. A size class chosen
  by length (50 kb, 75 kb, 100 kb) cannot predict it: N4 (70 kb, 69 labels) is fine at 14 pt and P1 (95 kb, 111
  labels) is not. This is why the 75 kb to 100 kb boundary question (gallery-video D-58) belongs here.
- Slow cases: T5 at 12 and 14 pt takes 46 to 74 s, pVNTG4 at 18 and 22 pt about 40 s (OV-473): the optimisers
  search hardest exactly where no solution exists.
- Radial placement never overprints or crosses, but the canvas grows a lot: P1 14 pt 1735 x 1580 (today 1706 x 990),
  tobacco 14 pt 2427 x 2206, T5 14 pt 2288 x 2530.
- Tracked figures that a new default for horizontal placement would change: 1 of 16 reference outputs
  (`circular_with_labels`) and 2 Gallery examples (HmmtDNA_basic_circular, HmmtDNA_ATskew). The tobacco Gallery
  example uses radial placement.

## 3. What other tools do (logs/survey-tools.md, with sources)

- Only two tools document an algorithm that copes with 100 to 200 labels: CGView.js / Proksee (labels sorted by
  priority, the first 50 placed collision-free by lengthening or tilting their lines, the rest dropped when they
  overlap) and Circos (labels stacked outward and "snuggled" along the ring; a label that still does not fit is not
  drawn).
- The usual answer to "too many labels" is hiding: SnapGene hides labels that do not fit, Geneious has
  "Hide excessive labels", MacVector "Enable Label Hiding", CGView `priorityMax`.
- Inside or outside the ring is an option in SnapGene (external plus internal), Geneious ("Inside", "Outside",
  "Inside or Outside") and Java CGView (`useInnerLabels`).
- Side-aligned lists exist in plasmid tools: CLC "Stacked" aligns all labels on each side of the circle; SimVector's
  "Bracket View" lists the MCS enzymes one above the other.
- Not documented anywhere we could read: Geneious's circular label algorithm, Benchling, OGDRAW. The Owner's
  Geneious Prime screenshots (LP-D04) are the evidence for Geneious: labels close to their features on either side
  of the ring, short straight leaders, oblique staircases, equal neighbours merged, text in the feature colour.

## 4. What the Owner asked for during the study (DECISIONS.md)

- LP-D02: plasmid maps stack the labels of a cluster as a vertical list, sometimes with one leader: acceptable look.
- LP-D03: a proposal shows many rendered images of each candidate.
- LP-D04: Geneious Prime places labels liberally (both sides of the ring, short straight leaders).
- LP-D05: also propose a mode that minimises the total leader length of all labels; no overlap, no dropped label.
- LP-D06: labels may leave the ellipse; a label keeps a minimum clearance from every other object (labels, ticks,
  features), not merely no overlap.
- LP-D07: leaders must not cross (hard constraint).
- LP-D09: a label must not run back toward the ring: the end of the text away from the leader is never closer to
  the ring ("axis") than the end the leader meets. Read as: outside labels run outward, inside labels run toward
  the centre (my reading; asked back in question 3).

So the target is: every label placed, a clearance kept to every object, no leader crossing, no label running
back toward the ring, and among such placements short leaders. The ellipse is not required.

## 5. The candidates (prototypes in `gbdraw/labels/proto/`, selected with `GBDRAW_LABEL_PROTO`)

| id | idea | never overprints | every label placed | no crossing | what it costs |
| --- | --- | --- | --- | --- | --- |
| Current | one ellipse, angle only | no | yes | no | overprint when the ellipse is full |
| A | today's placement; when it does not fit, font down and band up together (A1 font only, A2 band only, A3 drop) | yes | yes (A3 no) | not required by its test (0 on the cases it changed) | smaller text or larger canvas; first failing render is slow |
| B | tiers: outer flush columns reached through 4 px slots between inner rows | yes | yes | yes | width; leader bundles; staircase above 12 o'clock |
| C | free boxes outside the ring, one bounded straight leader, greedy by crowding | yes | no (drops) | only with `fallback=0` | drops 13 of 111 on P1 14 pt, 74 of 194 on T5 8 pt |
| D | two side columns hugging the ring, clockwise order kept | yes | yes | yes | height; "comb" of parallel leaders above and below the ring |
| E | flush lists for clusters (at most 8 rows), fan or bracket leaders; lone labels as today | yes | default drops; `fit=0` places all | yes (1 on T5 14 pt) | drops or height; one column per side only |
| G | Geneious look: beside the feature, inside the ring where nothing is drawn, short straight leaders, clearance to all objects | yes | no (drops) | only with `fallback=0` | drops on dense records and when the inside holds tracks; 30 to 70 s on dense records |
| H | D's ordered columns pulled in to 3 px from the axis line (leaders of 5 to 10 px where features are not crowded), plus an ordered inside shore when nothing is drawn there; a label moves inside only when that shortens the total (logs/design-H.md) | yes | yes | yes | vertical leaders stacked at 12 and 6 o'clock; height on very dense records; inside unused with default inner tracks |
| F | exact minimum total leader length on one ellipse, order kept (logs/design-F.md) | yes | yes | yes | the whole ellipse grows on dense records; leaders 2 x D's |
| R | shipped `--label_placement radial` | yes | yes | yes | canvas up to 2.5 x; text at every angle |

Scaffold notes for whoever continues: the hook in `labels/circular.py` (`_proto_design[0] in "BCDEFG"`) hands outer
labels to `proto/design_<letter>.py`; A wraps today's placement (`run_design_a`). Known scaffold faults:
`Item.span_bp` is the label text extent, not the feature length (E's drop order therefore ranks by label width; A3
reads the real span); `band_scale` had no effect on the CLI path until design_a.py scaled `outer_arena` itself.
`Ctx` knows the inside of the ring only through `radial_layout` (added for G).

## 6. Comparison (numbers from the renders)

Over the 18 records with gbdraw's default inner tracks (the last two columns are P1 at 14 pt):

| candidate | smallest gap between labels | leader crossings | labels running back toward the ring (LP-D09) | labels dropped | total leader length | canvas |
| --- | --- | --- | --- | --- | ---: | --- |
| Current | overlap (17.4 px deep) | 7 records, up to 2,532 pairs | 6 records, 34 labels | 0 | 12,731 px | 1706 x 990 |
| **H** | 3 px | 0 | 0 | 0 | 10,564 px | 1574 x 1224 |
| A | 0.1 px | 1 record (1 pair) | 3 records, 6 labels | 0 | 19,754 px (12.5 pt) | 1773 x 1130 |
| B | 3 px | 0 | 0 | 0 | 17,868 px | 1732 x 1083 |
| C (no crossing) | 3 px | 1 record (1 pair) | 8 records, 34 labels | 13 records, 476 labels | 7,506 px (87 labels) | 1577 x 1167 |
| D | 3 px | 0 | 0 | 0 | 12,526 px | 1624 x 1179 |
| E (default) | 3 px | 1 record (1 pair) | 1 record, 1 label | 9 records, 267 labels | 9,163 px (99 labels) | 1596 x 1315 |
| G (no crossing) | 3 px | 0 | 17 records, 53 labels | 13 records, 548 labels | 3,798 px (85 labels) | 1551 x 891 |
| F | 3 px | 0 | 0 | 0 | 27,594 px | 1860 x 1143 |
| R | not measured (rotated text) | 0 | 0 (radial text points outward) | 0 | 14,213 px | 1735 x 1580 |

Over all 23 records: the four constraints (clearance, no crossing, no label toward the ring, nothing dropped) are
met on every record only by B, D, F and H. H has the shortest total leader length of those four on 21 records
(not on tobacco 14 pt: D 19,830 px, H 20,019 px; MG1655 rRNA 8 pt: B 9,231 px, H 9,352 px). Against today's
placement H is 16 to 53 % shorter on the records today draws without overprint (lambda 14 pt 8,171 -> 4,845 px,
mtDNA 2,898 -> 1,738 px, pBR322 3,050 -> 1,507 px) and 1 to 14 % narrower; it is taller where columns grow
(P1 14 pt 990 -> 1224 px, tobacco 14 pt 988 -> 1436 px, T5 14 pt 992 -> 1985 px).
LP-D09 over all 23 records: G 22 records / 62 labels, C 10 / 38, today 7 / 37 (19 on T5 8 pt), A 4 / 7, E 1 / 1.

What the figures show (logs/review-figures.md, logs/design-H.md):
- Easy maps (up to about 70 labels at 14 pt): today's placement is clean and near the minimum for one ellipse
  (within -4 % to +8 % of F). H sets the labels right beside the features and halves the leader length.
- Overflowing maps (P1 14 pt, tobacco 14 pt): A keeps today's look with a smaller font and wider band; H and D
  place everything; H's leaders are shorter (P1 -16 %), both stack vertical leaders at 12 and 6 o'clock.
- Open centre (no GC tracks): H uses the inside (P1: 28 of 111 labels, total 6,699 px against D 12,526 px) and
  comes close to G's look; G is 30 to 40 % shorter per label but drops 9 labels on P1 and turns labels toward
  the ring.
- Clustered labels (pBR322 restriction sites, MG1655 rRNA operons): E's flush lists read best; its bracket leader
  matches the MCS lists of plasmid maps.
- T5 and the 153 kb plasmid at 14 pt (180 to 195 labels): no candidate gives a publishable figure.
- The objective (sum, squares, maximum) hardly changes a layout; room next to the features does.

H's open points (for the implementation plan): the vertical leaders at 12 and 6 o'clock (P1 14 pt: about 15
labels above the ring, longest leader 413 px); very dense records grow tall; leader to other label clearance is
2 px at the closest, not 3; the choice of inside labels is a greedy descent, not an optimum; no inside labels
when inner tracks leave a gap under 3 font sizes; tick-label boxes are recomputed in the prototype, not taken
from the tick renderer; G's clearance is measured from the feature band, while the axis stroke reaches 391.5 px
(H measures from the axis).

## 7. Recommendation and staging (as put to the Owner on the page)

- Direction: candidate H as the new placement.
- When labels do not fit: grow the columns (canvas) by default; dropping only when the user chooses it, with the
  count shown.
- Stage 1, 0.14.x: A ("today's look, never overprint") plus the 269 degree clamp fix. Figures that fit today stay
  byte-identical. With an explicit font size, keep the font and widen the band only (A2); with Auto, move font and
  band together (A). A production version needs a capacity pre-check (OV-473) and must report what it changed.
- Stage 2, 0.15.0: H as a new `--label_placement` value, opt-in first, default after a Gallery review.
- Stage 3: E's bracketed list as the way to write a tight cluster.
- Not recommended: B (leader bundles), C and G (drop labels, turn labels toward the ring; H keeps G's "beside
  the feature"), D (H gives the same guarantees with shorter leaders), F on the ellipse (its exact ordered
  packing is reusable inside H's columns).

## 8. Questions asked to the Owner (answers go to DECISIONS.md)

1. Is H the new placement? Recommended: yes; reducing the vertical leaders at 12 and 6 o'clock (inside top and
   bottom rows, slanted leaders) goes into the implementation plan.
2. When labels do not fit, what gives: (a) canvas, (b) font, (c) labels? Recommended: (a); (b) only with Auto font;
   (c) only on request, with the count.
3. Use the inside of the ring by default? Inside labels meet their leader at the ring end and run toward the
   centre (my reading of LP-D09; if "axis" meant the centre point, the inside cannot be used). Recommended: yes,
   only where nothing is drawn and only for labels whose move shortens the total.
4. Ship A in 0.14.x first? Recommended: yes (A2 with an explicit font, A with Auto, and a notice).
5. Extras: merge equal neighbouring labels (recommended as an option), feature-coloured text (not recommended),
   bracketed lists (after 0.15.0).
6. Auto font size by "largest size at which the labels fit" instead of record length (the 75 kb / 100 kb boundary)?
   Recommended: with the new placement in 0.15.0; the 100 kb change stays on hold until then.
7. Hide "hypothetical protein" and similar names by default? Recommended: no change (neither CLI nor Web hides them).

## 9. Bugs and faults found (IDs reserved in ID-BLOCKS.md, block OV-470..479)

- OV-456 (existing): overprint at 14 pt. Cause in section 1 (269 degree clamp plus previous-label-only overlap test,
  and no capacity fallback).
- OV-470: `config.toml` `labels.filtering.blacklist_keywords` (hypothetical, uncharacterized, putative, unknown) is
  never used by default: the CLI passes an empty list (`cli_utils/common.py:504-507` default `""`,
  `circular.py:1107-1112`, `linear.py:1426` test `is not None`) and the Web UI's filter mode defaults to None. Both
  surfaces agree, so this is a dead default, not a parity bug; hiding generic names by default is a Product question.
- OV-471: PNG export (cairosvg) draws text 5 to 9 % wider than gbdraw measured, so an end-anchored label at the left
  canvas edge loses its first letters. The study figures are rasterised with a 24 px margin to keep it out.
- OV-472: T5 at the Auto 8 pt: no overprint but 2,085 leader crossings, labels up to 629 px from their feature.
- OV-473: placement time when nothing fits: T5 at 12 and 14 pt 46 to 74 s, pVNTG4 at 18 and 22 pt about 40 s.
- Measured, no ID yet: today's placement leaves 0.1 to 0.4 px between label boxes on clean maps (lambda, mtDNA, T5
  8 pt), although `minimum_bbox_gap_px` asks for 3 px.

## 10. After the Owner's feedback on H (LP-D10): H2

Added 2026-10-12. Sections 6 to 8 describe the page as it stood with H. The Owner then named four defects of H
(DECISIONS.md LP-D10): sideways leaders near the poles, labels on ticks, sharp bends, long feature-side pieces on
inner labels. Two measures were added: object clearance against a label-free twin of each figure
(`scripts/lp_objects.py`) and leader bend angles (`scripts/lp_metrics.py`).

| over 23 records | today | H | H2 |
|---|---:|---:|---:|
| labels on another object | 0 | 6 | 0 |
| labels nearer than 3 px to another object | 0 | 237 | 0 |
| leaders over another object | 66 (T5, tobacco 14 pt) | 0 | 0 |
| leaders with a bend under 110 degrees (of 2,042) | 473 | 1,408 | 0 |
| records where it has the shortest total leader length of B, D, F, H, H2 | - | 0 | 23 |

H2 (design id `I`, logs/design-I.md): every label starts at the spot with the shortest leader (beside the feature
at 3 and 9 o'clock, right above or below it at 12 and 6), keeps the clockwise order, and is pushed along its stack
only when the spot is taken; a leader leaves the feature radially and bends once at the ring edge by at most 70
degrees; labels and leaders keep 3 px from every drawn object; inner labels follow the same rule inside the ring.
The recommendation on the page is now H2, with the same staging as in section 7. Its weak points: one-sided dense
records grow taller than under H (T5), towers at the poles remain, and the points listed in STATE.md.
