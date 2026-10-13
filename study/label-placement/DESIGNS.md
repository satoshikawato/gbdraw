# Circular label placement: record of the design prototypes A to Q

This folder is the permanent record of the design study of 2026-10-11 to 2026-10-13: prototypes A to N, P and Q for
placing the outer labels of Circular diagrams. The branch `study/label-placement` is local and is never merged into
`dev`. The study ended when the Owner chose to keep today's placement and add a fit check (LP-D41 to LP-D52); that
work shipped separately (PR #1049), so none of the prototype code here is production code.

Rendered pages and figures are not committed because of their size. They stay in the task folder
`/home/kawato/gbdraw-baselines/label-placement-20261011/`: `page/` (about 530 MB, rendered cases, metrics, page
builders) and `figures/` (about 90 MB, comparison figures with `figures/README.md`). The Owner decision log is
`DECISIONS.md` in that folder (LP-D01 to LP-D52); it holds long Japanese quotes and is not copied here.

Contents of this folder:

- `reports/`: `STUDY.md` (the study report), `design-*.md` (one report per design), `part2-breaking-points.md`,
  `survey-tools.md`, `review-figures.md`.
- `scripts/`: measurement scripts (`lp_metrics.py`, `lp_objects.py`, `check_hard.py`, ...) and `PROTO-SPEC.md`.
- `page/`: page builders and render helpers, plus small JSON inputs.

## Where the code lives and how to select a design

- Branch `study/label-placement` (head 1e289f347) holds A to J3, M, N, P and Q under `gbdraw/labels/proto/`
  (`design_a.py` ... `design_n.py`). Branch `study/label-placement-k` holds K and L (K 9 chars: c4f9e8a45, hook
  bc65af6dc, L 29a93fcd7); it is based on the J hook commit 6dbc21e09 and is not joined to the first branch.
- Selection: environment variable `GBDRAW_LABEL_PROTO` (first letter picks `design_<letter>.py`, so `A1`..`A3`
  share `design_a.py`) and `GBDRAW_LABEL_PROTO_OPTS` (`key=value,key=value`); see `gbdraw/labels/proto/__init__.py`.
  The hook is in `gbdraw/labels/circular.py`.
- J, J2, J3, P and Q are all `GBDRAW_LABEL_PROTO=J`; the options decide which one runs (see each section).
- Common scaffold commit: 8e63abe96 (hook, polyline leaders, shared geometry).

## A1, A2, A3: today's placement plus a never-overprint fallback

- What: keep today's ellipse placement; when labels do not fit, make the font smaller (A1), widen the band (A2) or
  drop labels (A3); plain `A` combines font and band.
- Commit: fd22b3dfb (with B to G). Select: `GBDRAW_LABEL_PROTO=A1|A2|A3|A`; options in `design_a.py`.
- Verdict: no explicit Owner pick. It never overprints but does not enforce the no-crossing rule (LP-D07) and
  gives larger canvases or smaller text. Its idea of "fix the fit, keep the look" returned in LP-D41 (keep today's
  placement, warn and ask), which became the fit check.
- Report: `reports/design-A.md`.

## B: staggered tiers

- What: outer flush columns reached through narrow slots between inner rows. Never overprints, no crossings.
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=B`.
- Verdict: not chosen. Costs width, leader bundles and a staircase above 12 o'clock; leader length 17,868 px
  against 12,731 px today over the 18 records (STUDY.md section 6). Superseded by H (LP-D08).
- Report: `reports/design-B.md`.

## C: free placement

- What: free label boxes outside the ring with one bounded straight leader, greedy by crowding.
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=C` (`fallback=0` for the no-crossing variant).
- Verdict: dropped labels (13 of 111 on P1 at 14 pt) and crossings break LP-D05 and LP-D07 (hard constraints
  recorded 2026-10-12). Not chosen (LP-D08).
- Report: `reports/design-C.md`.

## D: side columns

- What: two side columns hugging the ring, clockwise order kept.
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=D`.
- Verdict: meets all hard constraints but only uses the outside and makes a comb of parallel leaders at 12 and 6
  o'clock (LP-D08). It became the base of H.
- Report: `reports/design-D.md`.

## E: cluster callout lists

- What: flush lists for clusters of at most 8 rows with fan or bracket leaders; lone labels as today. Added after the
  Owner's plasmid-map input (LP-D02).
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=E` (`fit=0` places all labels instead of dropping).
- Verdict: drops by default, one column per side; not chosen.
- Report: `reports/design-E.md`.

## F: exact minimum total leader length

- What: exact solver for the minimum total leader length on one ellipse, order kept (LP-D05 asked for such a
  mode); LP-D06 made the off-ellipse variant the main one.
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=F` (`ellipse=1` for the comparison variant).
- Verdict: the ellipse stays too large on dense records and the free variant still crossed leaders (LP-D08). Not
  chosen; its measure (total leader length) stayed a column in every comparison.
- Report: `reports/design-F.md`.

## G: Geneious-like free placement

- What: label beside its feature, inside the ring where nothing is drawn, short straight leaders, clearance to all
  objects (LP-D04, LP-D06).
- Commit: fd22b3dfb. Select: `GBDRAW_LABEL_PROTO=G` (`fallback=0` variant).
- Verdict: drops labels on dense records, labels run back toward the ring (LP-D09 was written about G), slow (30 to
  70 s). Not chosen (LP-D08); its pole handling was admired in LP-D12.
- Report: `reports/design-G.md`.

## H: ordered free placement on both sides of the ring

- What: D's ordered columns pulled in to 3 px from the axis line, plus an ordered inside shore where nothing is
  drawn; a label goes inside only when that shortens the total.
- Commits: 5a04bbe1c, hook b23032bc0. Select: `GBDRAW_LABEL_PROTO=H` (`in0`, `clear=6` variants).
- Verdict: recommended on the first page (LP-D08), then criticised: horizontal leaders too long, ticks overlapped,
  acute bends (LP-D10).
- Report: `reports/design-H.md`.

## I (H2): every label beside its own feature

- What: H with the four LP-D10 requirements: pole labels right above or below the feature, clear of ticks, no bend
  under 110 degrees, inner labels placed by the same rule.
- Commit: f7b36b90b. Select: `GBDRAW_LABEL_PROTO=I` (`inside=0`, `clear=6` variants).
- Verdict: shortest total leader length of B, D, F, H on all 23 records; criticised for leader angles and for
  growth at 12 and 6 o'clock on dense maps (LP-D11).
- Report: `reports/design-I.md`.

## J (H3): leaders as straight as possible

- What: H2 with the turn of the leader as a second objective (cost = run + 0.8 px per degree of turn).
- Commit: 819443b79 (hook 6dbc21e09). Select: `GBDRAW_LABEL_PROTO=J` with the later options off
  (`sym=0,edge=0,face=0`; defaults changed as J2 and J3 were added, see `reports/design-J2.md`).
- Verdict: straight leaders 552 to 1,310 of 2,267 over 29 cases (LP-D11 result); the Owner asked for symmetric fans
  (LP-D12).
- Report: `reports/design-J.md`.

## J2: symmetric fans and the box-edge rule

- What: J with symmetric fans in clusters (LP-D12) and leaders that end on the label box edge (LP-D13, LP-D14, with
  `slide`); option `shift` moves a box toward the ring along its leader.
- Commits: 56d67f5c4, e72121355 (`shift`). Select: `GBDRAW_LABEL_PROTO=J` with `face=0` (this reproduces J2
  byte for byte); `shift=1` for the shift option.
- Verdict: not accepted as the default, "still needs work": a leader from above joined the lower vertex of the
  box (LP-D15 point 1).
- Report: `reports/design-J2.md`.

## J3: where a leader ends on the label box

- What: leaders end on the half of the box that faces the leader, at one visible gap (`face`, LP-D15, LP-D16);
  `mid` ends head-on leaders at the middle of the box end (LP-D18, LP-D19); `fan` re-seats a cluster as a straight
  fan; `strt` is the neighbour slide of LP-D16.
- Commits: b9ba95a0d (`face`), 99b0049b7 (`mid`), 9e51a4ea7 (`fan`).
  Select: `GBDRAW_LABEL_PROTO=J` with the defaults (`face=1,mid=1,fan=1,strt=1`).
- Verdict: not reaching zero on dense maps (far-half ends 158 to 27); the Owner stopped the pass-by-pass repair
  and asked for his own edit of two figures as input (LP-D20); contact points were then fixed by LP-D23.
- Report: `reports/design-J3.md`.

## K: a leader may pass behind a label box

- What: revision of H2 allowing a leader to pass behind a half-transparent label box, with underline/box/colour
  rendering options (the Owner's idea in LP-D11 point 3).
- Commits: c4f9e8a45 (hook 6dbc21e09 is shared with J), on `study/label-placement-k` only. Select:
  `GBDRAW_LABEL_PROTO=K` plus the rendering options of the report.
- Verdict: small gain (`p1_f14` canvas height 1,230 to 1,165). Underline with box read best; coloured text was hard
  to read. Not taken further.
- Report: `reports/design-K.md`.

## L: rows of several labels at 12 and 6 o'clock

- What: K's rules plus blocks of rows holding several labels at the poles; options `over` (leader behind a box) and
  `flip` (text may run back toward the ring, relaxing LP-D09).
- Commits: bc65af6dc (hook), 29a93fcd7, on `study/label-placement-k` only. Select: `GBDRAW_LABEL_PROTO=L`
  (`over=1,flip=0` was the recommended form).
- Verdict: `p1_f14` canvas 1,230 to 1,047 with both options; the Owner asked to see examples of when the relaxed
  rules matter (LP-D15) and then moved on to contact points (LP-D22).
- Report: `reports/design-L.md`.

## M: mock of the top and bottom limit

- What: a look study, not an engine: labels of each half on one arc under today's limit at the poles (LP-D34 to
  LP-D38); later reduced to one closed-form arc.
- Commits: 24bb3c64b, 40d87c477 (growing outline, "mid"), ad952919e (one arc, LP-D38). Select:
  `GBDRAW_LABEL_PROTO=M` (`grow`, `shape` options in `design_m.py`).
- Verdict: the Owner rejected the tent and bowl shapes: the arc must stay an arc, no inflection point (LP-D38,
  LP-D39); led to N. No separate `design-M.md` report exists; the docstring of `design_m.py` and LP-D34 to LP-D38
  hold the design.
- Report: none (see `reports/STUDY.md` and the module docstring).

## N: one circle around the ring

- What: all near ends on one circle concentric with the ring at the least radius that fits, capped by today's
  limit; dynamic programme per half; radial stub and one piece per leader (LP-D39, LP-D40).
- Commit: 1e289f347. Select: `GBDRAW_LABEL_PROTO=N` (`give=rows,limit` for the tighter-rows variant).
- Verdict: "not quite" (LP-D41): the Owner kept today's placement and asked for a warning with choices instead.
  That decision led to the fit check (LP-D41 to LP-D52, PR #1049).
- Report: `reports/design-N.md`.

## P: ink boxes, 12 contact points and one straight-first cluster construction

- What: label ink boxes, 12 fixed contact points per box (LP-D22 to LP-D25), one cluster construction that
  replaces symmetrize/straighten/fan; rounds 2 to 4 add an exclusive short-side arrival rule, the cluster arc
  (LP-D26), a pair DP and a smooth arc envelope.
- Commits: 6e10e44c4 (round 1), dfffabe07 (round 2), 6d6d51dea (round 3), 251ad452d (round 4). Select:
  `GBDRAW_LABEL_PROTO=J` with `port=1`.
- Verdict: the Owner asked for labels to follow the circumference (LP-D26) and then gave priority to that over
  straight middle leaders (LP-D27), which led to Q.
- Report: `reports/design-P.md`.

## Q: side clusters in arc coordinates

- What: runs of side labels placed by arc radius and minimum L1 displacement, contact point from the last piece;
  round 2 adds bounded row fit, round 3 allows vertices of 45 to 90 degrees; `rowem` and `rowmin` make the row
  distance follow the font size with a 4 px floor.
- Commits: b7b41fef4 (round 1), 5a600b7c7 (round 2), e5b647533 (round 3), 4a0818adb (`rowem`, LP-D30), 644e13f04
  (`rowmin`, LP-D31). Select: `GBDRAW_LABEL_PROTO=J` with `port=1,arcfirst=3` (plus `rowem=0.22,rowmin=4`).
- Verdict: "Q is good" (LP-D28); "Q is still fine, N is out" (LP-D40). It was the best design, but the Owner chose
  the fit check over adopting a new placement (LP-D41), so Q was not implemented.
- Report: `reports/design-Q.md`.

## Commit list of `study/label-placement` (oldest first)

8e63abe96 scaffold; fd22b3dfb A to G; 5a04bbe1c H; b23032bc0 hook; f7b36b90b I; 6dbc21e09 hook J and K; 819443b79 J;
56d67f5c4, e72121355 J2; b9ba95a0d, 99b0049b7, 9e51a4ea7 J3; 6e10e44c4, dfffabe07, 6d6d51dea, 251ad452d P;
b7b41fef4, 5a600b7c7, e5b647533, 4a0818adb, 644e13f04 Q; 24bb3c64b, 40d87c477, ad952919e M; 1e289f347 N.
Reproduce with `git log --oneline --reverse af6572809d..study/label-placement`.
