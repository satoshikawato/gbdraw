# Prototype contract for the label-placement study (read fully before you start)

Study only: no production change, no pull request, **no push, no git commit, no git stash** (the lead agent commits
locally). Background and findings so far: `/home/kawato/gbdraw-baselines/label-placement-20261011/HANDOFF-NOTES.md`.

- TASK_DIR = `/home/kawato/gbdraw-baselines/label-placement-20261011`
- Worktree (shared by several agents at once) = `/home/kawato/gbdraw-work/.worktrees/label-study`,
  branch `study/label-placement`. **Edit only the file(s) named in your brief.** Do not touch
  `gbdraw/labels/proto/common.py`, `__init__.py`, `gbdraw/labels/circular.py` or the renderer; if you need a change
  there, say exactly what and why in your report (a missing helper can live in your own file).
- Never `pip install`. Use `/home/kawato/micromamba/bin/python3.13`. `TMPDIR=/home/kawato/gbdraw-baselines/tmp`.
  Use absolute paths in every command (the shell's cwd resets between calls). `jq` is not installed.
- Write renders under `TASK_DIR/work/<your design letter>/`, your report to `TASK_DIR/logs/design-<letter>.md`.
- Playwright is not needed. A render takes 1 to 10 s; run batches in the foreground. Do not sleep-poll.

## How a prototype is selected

`GBDRAW_LABEL_PROTO=<id>` makes `prepare_label_list` (gbdraw/labels/circular.py) hand the outer labels to
`gbdraw/labels/proto/design_<letter>.py`. Options: `GBDRAW_LABEL_PROTO_OPTS="key=value,key=value"`
(read them with `ctx.opt("key", default)` or `ctx.options`). The scripts set both for you.

Designs B, C, D, E implement:

```python
def place(items: list[Item], ctx: Ctx) -> list[dict]:
    """Return one drawable label dict per placed label, built with common.emit().
    Append the text of every label you give up on to ctx.dropped. Put free-form facts in ctx.notes."""
```

`gbdraw/labels/proto/common.py` (read it) gives you:
- `Item`: `text, w, h` (measured text box, px), `theta` (feature midpoint, radians clockwise from 12 o'clock),
  `anchor` (point on the feature where today's leader ends), `feat_outer_r`, `span_bp` (feature length: the only
  importance measure available), `index` (clockwise order), `src`.
- `Ctx`: `radius` (R, 390), `r_elbow` (today's leader elbow radius, just outside the outermost ring), `r_keep`
  (no label box inside this radius), `gap` (minimum gap between label boxes, 3 px), `font_size`,
  `arc_x_radius`/`arc_y_radius` (today's label ellipse, the default band), `total_length`, `dropped`, `notes`.
- Geometry: `unit`, `polar`, `angle_of`, `angle_diff`, `box_from_left/right/center`, `box_on_ray`, `boxes_overlap`,
  `box_min_radius`, `nearest_point_on_box`, `segments_cross`, `segment_hits_box`, `segment_min_radius`,
  `Scene` (placed boxes and leaders with collision queries), `pack_ordered` (optimal ordered 1-D packing with
  per-item sizes: use it for any stacked list), `line_pitch`.
- `emit(item, box, leader_points, align="start"|"middle"|"end", extra=[polyline, ...])` builds the dict. Text is
  horizontal, vertically centred in `box`. `leader_points` is a polyline from the feature (`item.anchor` first) to
  the label; `None`/`[]` draws no leader. `extra` polylines are drawn in the leader style (brackets, list bars).
Coordinates are circle-local pixels, centre (0,0), y down. The canvas grows to fit whatever you place.
Only `--labels out` is in scope (inner labels are ignored by B to E). Embedded labels (curved text inside long
features) are kept as today; you only place the external ones.

## Running and measuring

```bash
S=/home/kawato/gbdraw-baselines/label-placement-20261011/scripts
/home/kawato/micromamba/bin/python3.13 $S/run_matrix.py --out /home/kawato/gbdraw-baselines/label-placement-20261011/work/<letter> \
    [--opts key=value,...] [--suffix _v2] p1_f14:<id> t5_f8:<id> ...        # <id> = cur | R | A1 | A2 | A3 | B | C | D | E
```
Each run writes `<case>__<id>.svg/.png/.metrics.json/.labels.json/.log` and prints a summary line. Cases are in
`$S/cases.json`: `p1_f14 p1_f12 p1_f10 p1_f8 t5_f8 t5_f14 n4_f14 tobacco_f8 tobacco_f14 mito_f14 plasmid37_f14
plasmid153_f8 plasmid153_f14` (font = `--label_font_size`; P1, T5, N4 are NCBI records outside the repository:
"extra"). Today's results are already in `TASK_DIR/work/base/<case>__cur.*` for some cases; render the others
there when you need them.

Metrics (`$S/lp_metrics.py`): `text_overprint_pairs` (boxes sharing >= 20 % of the smaller box),
`text_overlap_pairs` (any overlap > 0.5 px), `leader_crossing_pairs`, `leader_through_text`,
`leaders_over_features` (a leader segment after the first passes inside the feature ring), `leader_len_mean/max`,
`label_to_feature_px_mean/max`, `canvas_px`, and `proto` (requested / placed / dropped).

PNGs have a transparent background. To look at one, flatten it first, and look at 1:1 crops of dense regions, not
only at the whole figure:
```bash
/home/kawato/micromamba/bin/python3.13 $S/flat.py in.png out.png [scale] [x0 y0 x1 y1 as fractions of the image]
```

## Quality bar

The Owner will judge the LOOK. The figure must show what the design would really look like for that record: real
label text, real font, real leaders. Required on every case: `text_overprint_pairs = 0` and
`text_overlap_pairs = 0`; no leader drawn through another label's text; the number of dropped labels reported.
Then minimise, in this order: leader crossings, leaders running over the feature ring, labels far from their
feature, canvas growth. Keep it deterministic (no randomness) and under about 10 s per case.
It is a prototype: no tests, no options plumbing, no production code quality. Iterate on the look: render, flatten,
look at crops, adjust.

Standard cases every design must render before reporting (all with your final settings):
`p1_f14 p1_f8 t5_f8 t5_f14 n4_f14 tobacco_f8 tobacco_f14 mito_f14 plasmid37_f14 plasmid153_f8 plasmid153_f14`.

## Report (`TASK_DIR/logs/design-<letter>.md`, at most about 90 lines)

1. The rule in plain words (what a label's position is, how collisions are resolved, what happens when labels do
   not fit, what is bounded), precise enough for an engineer to reimplement.
2. Parameters and their values.
3. A metrics table over the standard cases (labels, dropped, overprint, overlaps, crossings, through-text, leaders
   over features, leader mean/max, canvas w x h), next to today's values where today's render exists.
4. Honest assessment from looking at the figures: where it reads better than today, where worse (leader clutter,
   labels far from features, canvas growth, lost reading order, long lists), and on which kind of record it is the
   right or the wrong tool.
5. What production would need: which existing modules change, what could be an opt-in option, rough size, risks.
6. Anything you could not do.

Final message to the lead: at most 20 lines (file paths, the table's headline numbers, the verdict). The lead reads
your report file; do not paste it.
