"""Design M: a mock of the top and bottom limit (LP-D34 to LP-D37).  A look study, not an engine.

The Owner: the labels need a limit at the top end and at the bottom end; the limit is what today's parameters set
(`--outer_label_y_radius_offset`, `--outer_label_x_radius_offset`) and the placement does not move it; leaders
never cross; the labels at the top and bottom ends stand toward the middle (12 and 6 o'clock), not out to the side.

Each half of the map (right, left) is one column in the clockwise order of the features, one label per row.  The
near ends of a column lie on one convex outline.  It starts as the smallest one there is: around the ring
(`r_keep`) and, where a stack stands above or below the ring, up to the end row of that stack next to the vertical
axis.  So labels stand next to the ring, and at 12 and 6 o'clock next to the axis.  The outline grows only where
labels need it: a label that must stand farther out to be reached (see the turn below) pushes the outline there,
and the outline becomes the smallest convex curve through those places.  It never passes today's label ellipse
(`xcap` at the equator, the limit at the poles).  When that growth does not settle, the outline is today's ellipse.

Every leader is a radial stub to the elbow radius and one straight piece that turns less than 90 degrees there.
Pieces that run from one circle to one convex outline around it, in the same order on both and never back into the
circle, cannot cross, and no piece reaches a label box, which lies outside the outline: that is what keeps the
leaders apart, not a test.  A leader is then drawn as one straight line where that line leaves the feature within
`scone` degrees of the ray and crosses no other leader (tested).

The rows are the cheapest ordered choice of open rows: squared distance from the height of the own feature, plus
`wc` times the squared degrees the piece turns beyond `cone`, plus (while the outline grows) `wx` times the squared
px a label stands outside it.  A row is open while the turn stays under the largest turn in use.  What gives way
when the labels do not fit on today's ellipse under the limit, in this order:
  1. the turn opens, from `cone` up to `conemax`;
  2. the rows of the whole figure shrink, not below `sqmin` px of ink gap;
  3. the limit moves out at that end by the least amount.  The notes then show `limit_top` / `limit_bottom` above
     `limit_asked`: the map does not fit under the limit that was set.

Options (GBDRAW_LABEL_PROTO_OPTS):
  ylim     limit: px from the map centre that no label box passes at the top and bottom   [today's: arc_y_radius + h]
  xcap     widest outline, px from the vertical axis at the equator                       [today's: arc_x_radius]
  rowem    ink gap of stacked rows as a fraction of the font size                         [0.22]
  rowmin   floor of that gap, px                                                          [4]
  sqmin    smallest ink gap when the rows do not fit under the limit                      [2.5]
  inkk     ink height of a row as a fraction of the measured text height                  [0.755]
  lgap     the leader stops this far before the label                                     [2.5]
  cone     turn at the elbow that costs nothing, degrees                                  [40]
  conemax  largest turn                                                                   [88]
  scone    largest angle between a straight leader and the feature's own ray              [70]
  wc, wx   weights of the degrees beyond `cone` and of the px outside the outline         [2, 0.3]
  grow     0: the outline is today's ellipse from the start                               [1]
"""

from __future__ import annotations

import math

import numpy as np

from .common import Ctx, Item, box_from_left, box_from_right, emit, polar, segment_min_radius, segments_cross

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
    inkk = ctx.opt("inkk", 0.755)
    gap0 = max(ctx.opt("rowmin", 4.0), ctx.opt("rowem", 0.22) * ctx.font_size - 0.1)   # as design J: em = the font size number
    lgap = ctx.opt("lgap", 2.5)
    cone, conemax, scone = ctx.opt("cone", 40.0), ctx.opt("conemax", 88.0), ctx.opt("scone", 70.0)
    wc, wx = ctx.opt("wc", 2.0), ctx.opt("wx", 0.3)
    a, re = ctx.r_keep, ctx.r_elbow
    lim0 = ctx.opt("ylim", ctx.arc_y_radius + max(i.h for i in items))
    xcap = max(a, ctx.opt("xcap", ctx.arc_x_radius))
    cols = {"R": sorted((i for i in items if math.sin(i.theta) >= 0), key=lambda i: i.theta),
            "L": sorted((i for i in items if math.sin(i.theta) < 0), key=lambda i: -i.theta)}
    cols = {side: col for side, col in cols.items() if col}
    ink = {side: [inkk * i.h for i in col] for side, col in cols.items()}
    tall = {side: max(v) for side, v in ink.items()}
    pitch = min(map(min, ink.values())) + gap0
    step = pitch / round(pitch / 0.25)                   # rows of the rule fall on the grid
    half = int((max(lim0, 0.5 * max(sum(h + gap0 for h in v) for v in ink.values())) + 300.0) / step)
    ys = (np.arange(2 * half + 1) - half) * step
    rows = np.arange(len(ys))

    def majorant(x: np.ndarray) -> np.ndarray:
        """The least concave curve over x, between the first and the last row that holds a place (x above the axis floor)."""
        held = np.flatnonzero(x > AXIS)
        if len(held) < 3:
            return x
        first, hull = int(held[0]), []
        part = x[first:int(held[-1]) + 1]
        for i in range(len(part)):
            while len(hull) > 1 and (part[hull[-1]] - part[hull[-2]]) * (i - hull[-2]) <= (part[i] - part[hull[-2]]) * (hull[-1] - hull[-2]):
                hull.pop()
            hull.append(i)
        out = x.copy()
        out[first:first + len(part)] = np.interp(np.arange(len(part)), hull, part[hull])
        return out

    def today(side: str, tops: tuple[float, float]) -> np.ndarray:
        """x of the near end per row on today's ellipse: `xcap` wide at the equator, as high as `tops` (top, bottom)."""
        e = np.maximum(0.0, np.abs(ys) - 0.5 * tall[side]) / np.where(ys < 0, tops[0], tops[1])
        return np.maximum(xcap * np.sqrt(np.maximum(0.0, 1.0 - e * e)), AXIS)

    def nearest(side: str, tops: tuple[float, float]) -> np.ndarray:
        """x of the near end per row on the smallest outline: around the ring and, where the stack ends above or below
        the ring, up to its end row next to the axis."""
        h = tall[side]
        x = np.maximum(np.sqrt(np.maximum(0.0, a * a - np.maximum(0.0, np.abs(ys) - 0.5 * h) ** 2)), AXIS)
        for top, sign in zip(tops, (-1.0, 1.0)):
            if top >= a + h:                                 # a row that stands clear of the ring at the axis
                x[half + round(sign * (top - 0.5 * h) / step)] = AXIS + 0.01
        return majorant(x)

    def frame(it: Item, side: str) -> tuple[float, float, float, float]:
        th = it.theta if side == "R" else -it.theta              # the left column is worked in the mirror
        return math.sin(th), math.cos(th), re * math.sin(th), -re * math.cos(th)

    def turn(it: Item, side: str, x: np.ndarray) -> np.ndarray:
        """Degrees the leader turns at the elbow of `it`, per row, for the near end at x."""
        s, c, ex, ey = frame(it, side)
        vx, vy = x - lgap - ex, ys - ey
        return np.degrees(np.arctan2(np.abs(vx * c + vy * s), vx * s - vy * c))

    def bounds(side: str, top: float, bottom: float) -> tuple[list[int], list[int]]:
        return ([math.ceil((-top + 0.5 * h) / step - 1e-6) + half for h in ink[side]],
                [math.floor((bottom - 0.5 * h) / step + 1e-6) + half for h in ink[side]])

    def seps(side: str, gap: float) -> list[int]:
        return [math.ceil((0.5 * (p + q) + gap) / step - 1e-6) for p, q in zip(ink[side], ink[side][1:])]

    def fits(deg: float, gap: float, top: float = lim0, bottom: float = lim0) -> bool:
        """Every column fits on today's ellipse under these limits."""
        return all(_first_fit([np.flatnonzero(turn(it, side, today(side, (top, bottom))) <= deg) for it in col],
                              seps(side, gap), *bounds(side, top, bottom)) for side, col in cols.items())

    # what gives way on today's ellipse, in order: the turn, the rows, the limit
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

    def rows_on(side: str, places: list[tuple[np.ndarray, np.ndarray]]) -> list[int] | None:
        """The cheapest rows of a column; per label the turn and an extra cost per row."""
        lo, hi = bounds(side, *lims)
        cost = [np.where((rows >= lo[k]) & (rows <= hi[k]) & (t <= deg), (ys + a * math.cos(it.theta)) ** 2 + wc * np.maximum(0.0, t - cone) ** 2 + c, INF)
                for k, (it, (t, c)) in enumerate(zip(cols[side], places))]
        return _cheapest(cost, seps(side, gap))

    def grown(side: str, tops: tuple[float, float]) -> list[tuple[float, float, float]] | None:
        """Rows on the outline that starts as the smallest one and grows where labels need it; None when it does not settle."""
        cap = today(side, (lims[0], lims[1]))
        base = np.minimum(nearest(side, tops), cap)
        for _ in range(16):
            places, stand = [], []
            for it in cols[side]:
                s, c, ex, ey = frame(it, side)
                at, t, extra = base, turn(it, side, base), np.zeros(len(ys))
                for target in np.linspace(cone, deg, 6):             # stand outside the outline to turn no more than `target`
                    lo, hi, tc = base.copy(), cap.copy(), math.tan(math.radians(90.0 - target))
                    for sg in (1.0, -1.0):                           # the two sides of the cone: gx * (x - ex) + gy * (y - ey) >= 0
                        gx, gy = s - sg * tc * c, -c - sg * tc * s
                        if abs(gx) < 1e-9:
                            lo = np.where(gy * (ys - ey) >= 0, lo, INF)
                        elif gx > 0:
                            lo = np.maximum(lo, ex - gy * (ys - ey) / gx + lgap)
                        else:
                            hi = np.minimum(hi, ex - gy * (ys - ey) / gx + lgap)
                    new = np.where(lo <= hi, wx * (lo - base) ** 2 + wc * max(0.0, target - cone) ** 2, INF)
                    take = new < np.where(t <= deg, wc * np.maximum(0.0, t - cone) ** 2 + extra, INF)
                    at, t, extra = np.where(take, lo, at), np.where(take, target, t), np.where(take, wx * (lo - base) ** 2, extra)
                places.append((t, extra))
                stand.append(at)
            idx = rows_on(side, places)
            if idx is None:
                return None
            want = base.copy()
            want[idx] = np.maximum(want[idx], [float(x[i]) for i, x in zip(idx, stand)])
            if float(np.max(want - base)) < 0.25:
                return [(float(ys[i]), float(base[i]), float(turn(it, side, base)[i])) for i, it in zip(idx, cols[side])]
            base = np.minimum(majorant(want), cap)
        return None

    def widest(side: str, tops: tuple[float, float]) -> list[tuple[float, float, float]] | None:
        """Rows on today's ellipse, as high as the stack."""
        x = today(side, tops)
        ts = [turn(it, side, x) for it in cols[side]]
        idx = rows_on(side, [(t, np.zeros(len(ys))) for t in ts])
        return None if idx is None else [(float(ys[i]), float(x[i]), float(t[i])) for i, t in zip(idx, ts)]

    plan, how = {}, {}
    for side in cols:
        for way in ((grown, widest) if ctx.opt("grow", 1.0) else (widest,)):
            tops, best = (lims[0], lims[1]), None
            for _ in range(4):                           # the outline follows the stack, the stack the outline
                got = way(side, tops)
                if got is None:                          # a lower stack closes rows this column needs: keep the one before
                    break
                best = got
                tops = (max([a] + [-y + 0.5 * h for (y, _, _), h in zip(best, ink[side]) if y < 0]),
                        max([a] + [y + 0.5 * h for (y, _, _), h in zip(best, ink[side]) if y >= 0]))
            if best is not None:
                plan[side], how[side] = best, way.__name__
                break

    placed = []                                          # [item, box, align, leader points, turn at the elbow]
    for side, col in cols.items():
        for it, h, (y, x0, bend) in zip(col, ink[side], plan[side]):
            if side == "R":
                box, end, align = box_from_left(x0, y, it.w, h), (x0 - lgap, y), "start"
            else:
                box, end, align = box_from_right(-x0, y, it.w, h), (-x0 + lgap, y), "end"
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
                     gap_asked=round(gap0, 2), gap=round(gap, 2), outline=how, ring=round(a, 1), xcap=round(xcap, 1),
                     x_widest={side: round(max(x for _, x, _ in p), 1) for side, p in plan.items()},
                     x_at_ends={side: (round(p[0][1], 1), round(p[-1][1], 1)) for side, p in plan.items()},
                     turn_max=round(deg, 1), turn_free=cone,
                     straight=sum(1 for p in placed if len(p[3]) == 2), turned_over_70=sum(1 for p in placed if len(p[3]) == 3 and p[4] > 70.0),
                     crossings=crossing, arc_y_radius=round(ctx.arc_y_radius, 1), arc_x_radius=round(ctx.arc_x_radius, 1),
                     n_right=len(cols.get("R", ())), n_left=len(cols.get("L", ())))
    return [emit(it, box, pts, align=align, tag="M") for it, box, align, pts, _ in placed]
