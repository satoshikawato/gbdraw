"""Design M: a mock of the top and bottom limit (LP-D34 to LP-D38).  A look study, not an engine.

The Owner: the labels need a limit at the top end and at the bottom end; the limit is what today's parameters set
(`--outer_label_y_radius_offset`, `--outer_label_x_radius_offset`) and the placement does not move it; leaders
never cross; the labels stay on one arc, and at the top and bottom ends they stand toward the middle (12 and 6
o'clock).

Each half of the map (right, left) is one column in the clockwise order of the features, one label per row.  The
near ends of a column lie on one arc: today's label ellipse (`xcap` at the equator, the limit at the poles), drawn
through the middle of the rows, so that the highest and the lowest row under the limit stand next to the vertical
axis.  `shape` bends that arc: 2 is the ellipse, a smaller number pulls it toward its four ends.

Every leader is a radial stub to the elbow radius and one straight piece that turns less than 90 degrees there.
Pieces that run from one circle to one convex arc around it, in the same order on both and never back into the
circle, cannot cross, and no piece reaches a label box, which lies outside the arc: that is what keeps the leaders
apart, not a test.  A leader is then drawn as one straight line where that line leaves the feature within `scone`
degrees of the ray and crosses no other leader (tested).

The rows are the cheapest ordered choice of open rows: squared distance from the wanted row, plus `wc` times the
squared degrees the piece turns beyond `cone`.  The wanted row is where the feature's own ray meets the arc.  A row
is open while the turn stays under the largest turn in use.  What gives way when the labels do not fit on the arc
under the limit, in this order:
  1. the turn opens, from `cone` up to `conemax`;
  2. the rows of the whole figure shrink, not below `sqmin` px of ink gap;
  3. the limit moves out at that end by the least amount.  The notes then show `limit_top` / `limit_bottom` above
     `limit_asked`: the map does not fit under the limit that was set.

Options (GBDRAW_LABEL_PROTO_OPTS):
  ylim     limit: px from the map centre that no label box passes at the top and bottom   [today's: arc_y_radius + h]
  xcap     the arc at the equator, px from the vertical axis                              [today's: arc_x_radius]
  shape    exponent of the arc: 2 = ellipse                                               [2]
  rowem    ink gap of stacked rows as a fraction of the font size                         [0.22]
  rowmin   floor of that gap, px                                                          [4]
  sqmin    smallest ink gap when the rows do not fit under the limit                      [2.5]
  lgap     the leader stops this far before the label                                     [2.5]
  cone     turn at the elbow that costs nothing, degrees                                  [40]
  conemax  largest turn                                                                   [88]
  scone    largest angle between a straight leader and the feature's own ray              [70]
  wc       weight of the degrees beyond `cone`                                            [2]
"""

from __future__ import annotations

import math

import numpy as np

from .common import Ctx, Item, polar, segment_min_radius, segments_cross
from .design_j import _emit_ink, _ink_rel

INF = float("inf")
AXIS = 6.0            # nearest place of a near end to the vertical axis, px


def _first_fit(open_rows: list[np.ndarray], sep: list[int], lo: list[int], hi: list[int]) -> bool:
    """True when every label gets an open row inside [lo, hi] below the row of the label before it."""
    at = -10 ** 9
    for k, rows in enumerate(open_rows):
        at = max(at + (sep[k - 1] if k else 0), lo[k])
        pos = int(np.searchsorted(rows, at))
        if pos == len(rows) or rows[pos] > hi[k]:
            return False
        at = int(rows[pos])
    return True


def _cheapest(cost: list[np.ndarray], sep: list[int]) -> list[int] | None:
    """Row index per label: increasing by at least `sep`, smallest total cost; None when no choice exists."""
    idx = np.arange(len(cost[0]))
    total, back = cost[0], []
    for k in range(1, len(cost)):
        run = np.minimum.accumulate(total)
        arg = np.maximum.accumulate(np.where(total == run, idx, 0))
        best, frm = np.full(len(idx), INF), np.zeros(len(idx), dtype=int)
        best[sep[k - 1]:], frm[sep[k - 1]:] = run[:-sep[k - 1]], arg[:-sep[k - 1]]
        total = best + cost[k]
        back.append(frm)
    out = [int(np.argmin(total))]
    if total[out[0]] == INF:
        return None
    for frm in reversed(back):
        out.append(int(frm[out[-1]]))
    return out[::-1]


def _least(fits, lo: float, hi: float, eps: float) -> float | None:
    """Smallest value in [lo, hi] for which `fits` holds (it holds from there on), None when not even at hi."""
    if not fits(hi):
        return None
    if fits(lo):
        return lo
    while hi - lo > eps:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if fits(mid) else (mid, hi)
    return hi


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    if not items:
        return []
    gap0 = max(ctx.opt("rowmin", 4.0), ctx.opt("rowem", 0.22) * ctx.font_size - 0.1)   # as design J: em = the font size number
    lgap = ctx.opt("lgap", 2.5)
    cone, conemax, scone = ctx.opt("cone", 40.0), ctx.opt("conemax", 88.0), ctx.opt("scone", 70.0)
    wc, shape = ctx.opt("wc", 2.0), ctx.opt("shape", 2.0)
    a, re = ctx.r_keep, ctx.r_elbow
    lim0 = ctx.opt("ylim", ctx.arc_y_radius + max(i.h for i in items))
    xcap = max(a, ctx.opt("xcap", ctx.arc_x_radius))
    cols = {"R": sorted((i for i in items if math.sin(i.theta) >= 0), key=lambda i: i.theta),
            "L": sorted((i for i in items if math.sin(i.theta) < 0), key=lambda i: -i.theta)}
    cols = {side: col for side, col in cols.items() if col}
    ink = {side: [(-top, bottom, right - left) for left, top, right, bottom in (_ink_rel(i.text, float(i.src.get("font_size", ctx.font_size)), "start") for i in col)]
           for side, col in cols.items()}                    # per label: px of ink above and below the middle of its row, and its width
    tall = {side: tuple(map(max, zip(*v))) for side, v in ink.items()}
    pitch = min(u + d for v in ink.values() for u, d, _ in v) + gap0
    step = pitch / round(pitch / 0.25)                   # rows of the rule fall on the grid
    half = int((max(lim0, 0.5 * max(sum(u + d + gap0 for u, d, _ in v) for v in ink.values())) + 300.0) / step)
    ys = (np.arange(2 * half + 1) - half) * step
    rows = np.arange(len(ys))

    def arc(side: str, top: float, bottom: float) -> np.ndarray:
        """x of the near end per row on the label arc: `xcap` at the equator, on the axis in the end rows under the limits."""
        ends = np.floor(np.where(ys < 0, top - tall[side][0], bottom - tall[side][1]) / step + 1e-6) * step   # the end rows, on the grid
        e = np.minimum(1.0, np.abs(ys) / ends)
        return np.maximum(xcap * (1.0 - e ** shape) ** (1.0 / shape), AXIS)

    def frame(it: Item, side: str) -> tuple[float, float, float, float]:
        th = it.theta if side == "R" else -it.theta              # the left column is worked in the mirror
        return math.sin(th), math.cos(th), re * math.sin(th), -re * math.cos(th)

    def turn(it: Item, side: str, x: np.ndarray) -> np.ndarray:
        """Degrees the leader turns at the elbow of `it`, per row, for the near end at x."""
        s, c, ex, ey = frame(it, side)
        vx, vy = x - lgap - ex, ys - ey
        return np.degrees(np.arctan2(np.abs(vx * c + vy * s), vx * s - vy * c))

    def bounds(side: str, top: float, bottom: float) -> tuple[list[int], list[int]]:
        return ([math.ceil((-top + u) / step - 1e-6) + half for u, _, _ in ink[side]],
                [math.floor((bottom - d) / step + 1e-6) + half for _, d, _ in ink[side]])

    def seps(side: str, gap: float) -> list[int]:
        return [math.ceil((p[1] + gap + q[0]) / step - 1e-6) for p, q in zip(ink[side], ink[side][1:])]

    def fits(deg: float, gap: float, top: float = lim0, bottom: float = lim0) -> bool:
        """Every column fits on the arc under these limits."""
        return all(_first_fit([np.flatnonzero(turn(it, side, arc(side, top, bottom)) <= deg) for it in col],
                              seps(side, gap), *bounds(side, top, bottom)) for side, col in cols.items())

    # what gives way, in order: the turn, the rows, the limit
    gap, lims = gap0, [lim0, lim0]
    deg = _least(lambda v: fits(v, gap0), cone, conemax, 0.5) or conemax
    while gap > min(gap0, ctx.opt("sqmin", 2.5)) and not fits(deg, gap):
        gap = max(ctx.opt("sqmin", 2.5), gap - step)
    if not fits(deg, gap):
        even = _least(lambda v: fits(deg, gap, lim0 + v, lim0 + v), 0.0, (half - 2) * step - lim0, step)
        if even is None:
            raise RuntimeError("design M: the labels do not fit at any limit of the grid")
        lims[0] = lim0 + _least(lambda v: fits(deg, gap, lim0 + v, lim0 + even), 0.0, even, step)
        lims[1] = lim0 + _least(lambda v: fits(deg, gap, lims[0], lim0 + v), 0.0, even, step)

    plan = {}
    for side, col in cols.items():
        x, (lo, hi) = arc(side, *lims), bounds(side, *lims)
        ts = [turn(it, side, x) for it in col]
        idx = _cheapest([np.where((rows >= lo[k]) & (rows <= hi[k]) & (t <= deg), (ys - ys[int(np.argmin(t))]) ** 2 + wc * np.maximum(0.0, t - cone) ** 2, INF)
                         for k, t in enumerate(ts)], seps(side, gap))
        plan[side] = [(float(ys[i]), float(x[i]), float(t[i])) for i, t in zip(idx, ts)]

    placed = []                                          # [item, box, align, leader points, turn at the elbow]
    for side, col in cols.items():
        for it, (u, d, w), (y, x0, bend) in zip(col, ink[side], plan[side]):
            if side == "R":
                box, end, align = (x0, y - u, x0 + w, y + d), (x0 - lgap, y), "start"
            else:
                box, end, align = (-x0 - w, y - u, -x0, y + d), (-x0 + lgap, y), "end"
            placed.append([it, box, align, [it.anchor, polar(re, it.theta), end], bend])

    def crosses(pts, k: int) -> bool:
        return any(segments_cross(p, q, r, s) for j, other in enumerate(placed) if j != k
                   for p, q in zip(pts, pts[1:]) for r, s in zip(other[3], other[3][1:]))

    tscone = math.tan(math.radians(90.0 - scone))
    for k in sorted(range(len(placed)), key=lambda k: placed[k][4]):   # straight where the line is free, least turned first
        it, end = placed[k][0], placed[k][3][-1]
        s, c, ra = math.sin(it.theta), math.cos(it.theta), math.hypot(*it.anchor)
        du = end[0] * s - end[1] * c - ra
        if du > 0 and du >= tscone * abs(end[0] * c + end[1] * s) and segment_min_radius(it.anchor, end) >= ra - 0.5 \
                and not crosses([it.anchor, end], k):
            placed[k][3] = [it.anchor, end]
    crossing = sum(1 for k, p in enumerate(placed) if crosses(p[3], k)) // 2
    ctx.notes.update(limit_asked=round(lim0, 1), limit_top=round(lims[0], 1), limit_bottom=round(lims[1], 1),
                     gap_asked=round(gap0, 2), gap=round(gap, 2), shape=shape, ring=round(a, 1), xcap=round(xcap, 1),
                     x_at_ends={side: (round(p[0][1], 1), round(p[-1][1], 1)) for side, p in plan.items()},
                     turn_max=round(deg, 1), turn_free=cone,
                     straight=sum(1 for p in placed if len(p[3]) == 2), turned_over_70=sum(1 for p in placed if len(p[3]) == 3 and p[4] > 70.0),
                     crossings=crossing, arc_y_radius=round(ctx.arc_y_radius, 1), arc_x_radius=round(ctx.arc_x_radius, 1),
                     n_right=len(cols.get("R", ())), n_left=len(cols.get("L", ())))
    return [_emit_ink(it, box, pts, align, "M", ctx.font_size) for it, box, align, pts, _ in placed]
