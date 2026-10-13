"""Design M: a mock of the top and bottom limit (LP-D34).  A look study, not an engine: no legality tests.

The Owner: the labels need a limit at the top end and at the bottom end.  Each half of the map (right, left) is one
column in the clockwise order of the features.  The rows are the bounded least-squares stack: every label as near
the height of its own feature as the rows above and below allow, no row beyond the limit.  The near end of each
label lies on one smooth outline per quarter: the ellipse through the ring's side (`r_keep`) and the top of that
quarter's stack, so a stack that stays beside the ring hugs it and a stack that climbs to the limit forms the arch
today's placement has.  A label that its rows move away from its feature moves outward as well, until its near end
is seen from the feature within `cone` degrees of the feature's own ray (a leader never runs along the ring).
Leaders are one straight line where that holds from the feature itself, else today's shape: a radial stub to the
elbow radius and one straight piece.

Options (GBDRAW_LABEL_PROTO_OPTS):
  ylim     limit: px from the map centre that no label box passes at the top and bottom   [today's: arc_y_radius + h]
  rowem    ink gap of stacked rows as a fraction of the font size                         [0.22]
  rowmin   floor of that gap, px                                                          [4]
  inkk     ink height of a row as a fraction of the measured text height                  [0.755]
  squeeze  0: a column whose rows do not fit moves the limit out equally at both ends; 1: its rows shrink to fit,
           not below `sqmin` px of ink gap                                                [0]
  sqmin    smallest ink gap with squeeze=1                                                [2.6]
  lgap     the leader stops this far before the label                                     [2.5]
  cone     largest angle between a leader piece and the feature's own ray, degrees        [70]
"""

from __future__ import annotations

import math

from .common import Ctx, Item, box_from_left, box_from_right, emit, pack_ordered, polar, segment_min_radius


def _rows(col: list[Item], own: list[float], pitch: list[float], lim: float) -> tuple[list[float], float]:
    """Row centres top to bottom inside [-lim, lim]: the bounded isotonic fit (the unbounded fit, clipped in the
    coordinates in which the rows are merely ordered).  Returns the centres and the limit used."""
    n = len(col)
    cum = [0.0]
    for k in range(n - 1):
        cum.append(cum[-1] + 0.5 * (pitch[k] + pitch[k + 1]))
    span = cum[-1] + 0.5 * (pitch[0] + pitch[-1])
    lim = max(lim, 0.5 * span)
    ys = pack_ordered(own, pitch)
    lo, hi = -lim + 0.5 * pitch[0], lim - 0.5 * pitch[-1] - cum[-1]
    return [min(max(y - c, lo), hi) + c for y, c in zip(ys, cum)], lim


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    if not items:
        return []
    inkk = ctx.opt("inkk", 0.755)
    gap = max(ctx.opt("rowmin", 4.0), ctx.opt("rowem", 0.22) * ctx.font_size - 0.1)   # as design J: em = the font size number
    lgap = ctx.opt("lgap", 2.5)
    tcone = math.tan(math.radians(90.0 - ctx.opt("cone", 70.0)))
    re = ctx.r_elbow

    def seen(it: Item, r: float, p) -> bool:
        """p lies within the cone around the ray of `it`, seen from the point of that ray at radius r."""
        s, c = math.sin(it.theta), math.cos(it.theta)
        du = p[0] * s - p[1] * c - r
        return du >= tcone * abs(p[0] * c + p[1] * s) and du > 0.0
    a = ctx.r_keep
    hmax = max(i.h for i in items)
    lim0 = ctx.opt("ylim", ctx.arc_y_radius + hmax)
    cols = {"R": sorted((i for i in items if math.sin(i.theta) >= 0), key=lambda i: i.theta),
            "L": sorted((i for i in items if math.sin(i.theta) < 0), key=lambda i: -i.theta)}
    plan, lim = {}, lim0
    for _ in range(2):                                   # second pass: both columns under the limit the fuller one needs
        for side, col in cols.items():
            if not col:
                continue
            ink = [inkk * i.h for i in col]
            pitch = [h + gap for h in ink]
            if ctx.opt("squeeze", 0.0) and sum(pitch) > 2 * lim0:
                pitch = [max(h + ctx.opt("sqmin", 2.6), 2 * lim0 / len(col)) for h in ink]
            ys, used = _rows(col, [-a * math.cos(i.theta) for i in col], pitch, lim)
            plan[side] = (col, ink, ys)
            lim = max(lim, used)
    out: list[dict] = []
    tops = {}
    for side, (col, ink, ys) in plan.items():
        b_top = max([a] + [-y + 0.5 * h for y, h in zip(ys, ink) if y < 0])
        b_bot = max([a] + [y + 0.5 * h for y, h in zip(ys, ink) if y >= 0])
        tops[side] = (round(b_top, 1), round(b_bot, 1))
        for it, h, y in zip(col, ink, ys):
            b = b_top if y < 0 else b_bot
            e = max(0.0, abs(y) - 0.5 * h)               # the box edge nearer the equator
            x0 = max(a * math.sqrt(max(0.0, 1.0 - (e / b) ** 2)), 6.0)
            sgn = 1.0 if side == "R" else -1.0
            for _step in range(200):                     # outward until the near end is seen from the elbow of its feature
                if seen(it, re, (sgn * (x0 - lgap), y)):
                    break
                x0 += 2.0
            if side == "R":
                box, end, align = box_from_left(x0, y, it.w, h), (x0 - lgap, y), "start"
            else:
                box, end, align = box_from_right(-x0, y, it.w, h), (-x0 + lgap, y), "end"
            ra = math.hypot(*it.anchor)
            straight = seen(it, ra, end) and segment_min_radius(it.anchor, end) >= ra - 0.5
            pts = [it.anchor, end] if straight else [it.anchor, polar(re, it.theta), end]
            out.append(emit(it, box, pts, align=align, tag=side))
    ctx.notes.update(limit=round(lim, 1), limit_asked=round(lim0, 1), a=round(a, 1), gap=round(gap, 2), outline=tops,
                     arc_y_radius=round(ctx.arc_y_radius, 1), arc_x_radius=round(ctx.arc_x_radius, 1),
                     n_right=len(cols["R"]), n_left=len(cols["L"]))
    return out
