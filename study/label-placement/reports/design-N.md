# Design N: one circle around the ring, as near the ring as the labels allow (2026-10-13, session d8e02249)

Study prototype, `gbdraw/labels/proto/design_n.py` on the study branch `study/label-placement` (worktree
`/home/kawato/gbdraw-work/.worktrees/label-study`, never pushed). Dispatcher letter `N` (`GBDRAW_LABEL_PROTO=N`;
`gbdraw/labels/circular.py` accepts it). Built from `logs/AUDIT-N-20261013.md`; the Owner's rules it realises are
listed at the top of the module's docstring with their LP-D ids.

## Model

- All near ends (where the leader meets the label) lie on ONE circle of radius ρ concentric with the ring. ρ_min =
  ring edge + clearance + 0.5; ρ ≤ `ylim` − tallest ink height, `ylim` = today's limit (`arc_y_radius` + text
  height, the value `--outer_label_y_radius_offset` sets). At the poles the circle reaches the vertical axis, so the
  outline is a dome (LP-D39) and the end labels stand toward 12 and 6 o'clock (LP-D37).
- Per half (right: sin θ ≥ 0; left mirrored), labels in clockwise order of their features. For label k the candidate
  near ends are the points of the circle at polar angles ψ = θ_k + j·`step` within the window in which the last piece
  does not run back into the ring (acos(r_elbow/ρ) + 3°). The contact point follows the direction of the last piece:
  within 45° of horizontal the middle of the short side; steeper the facing vertex (natural up to 80°) or the middle of
  the facing long side beyond.
- Cost of a candidate: ρ·|ψ − θ| (sideways displacement) + `wc`·max(0, turn − `cone`)² + `vpen` for a contact point
  that is not the natural one. Feasible only when the piece points outward, turn ≤ T, the box stays on its half
  (`axis`), under `ylim`, clear of the ring, and clear of the raster objects (ticks, tick labels, tracks: `_Geo` of
  design J).
- A dynamic programme over consecutive labels (admissible pair = later ψ, boxes apart by the row gap in x or y, neither
  leader through the other's box) gives the least total cost; constraints between non-adjacent labels are added
  lazily (8 rounds). No two last pieces can cross because the near ends keep the feature order and no piece runs
  back; crossings are counted at the end all the same.
- Search: for turn limit T = `cone`, `cone`+8, … `conemax` (within `polezone` of a pole `conemax` from the start),
  the least ρ by bisection on [ρ_min, ylim − h_max]. If nothing fits under the limit: `give=limit` moves the limit
  out by `lstep` px per step with the rows at the rule; `give=rows,limit` shrinks the row gap by 0.25 px steps down
  to `sqmin` first.
- Leader (LP-D40): feature anchor → radial stub to the elbow at the ring edge + `stub` (as Q) → one piece to the
  contact point, stopping `lgap` before it. One line only when it lies within `scone` = 2° of the ray.
- Emission through `_emit_ink(..., "N", ...)` of design J (ink box, `proto_tag` N, `leader_points`).

## Options

`ylim`, `rowem` 0.22, `rowmin` 4, `give` limit, `sqmin` 2.5, `lgap` 2.5, `lpad` 2.4, `cone` 40, `conemax` 88,
`polezone` 35, `scone` 2, `wc` 2, `vpen` 4, `step` 0.15, `lstep` 8, `objects` 1, `axis` 6 (details in the docstring).

## How to run

```
cd /home/kawato/gbdraw-baselines/label-placement-20261011
TMPDIR=/home/kawato/gbdraw-baselines/tmp /home/kawato/micromamba/bin/python3.13 page/run_extra.py \
  --out work/N2 [--jobs 4] [--opts give=rows,limit --suffix _rows] CASE:N ...
```
Outputs `work/N2/<case>__N[suffix].{svg,png,labels.json,metrics.json,log}`; the design's notes are in
`metrics.json` under `proto.notes` (rho, rho_min, turn_max, gap, limit_asked, limit, straight, turned_over_70,
crossings, seconds). Checks: `scripts/lp_objects.py work/N2 work/bare_links`, then `scripts/check_hard.py
'work/N2/*.metrics.json'`, `scripts/leader_radial.py`, `scripts/arc_shape.py`, `scripts/pole_order.py`,
`scripts/row_gap_em.py` (PYTHONPATH = the study worktree). Pictures: `owner-edit/analysis/panels.py` and `whole.py`
(`work/N2/view/`).

## Versions in this session

- v1: ray track at the poles (staircase) — withdrawn (LP-D39).
- v2, v3: per-angle radius raised where dense, pole rays — withdrawn (LP-D39: no inflection, dome).
- v4: one concentric circle; straight pass within 15° — the lines crossing the axis at an angle were rejected (LP-D40).
- v5 (current): v4 with `scone` 2°.

## Results (v5)

(see the table appended below after the batch)
