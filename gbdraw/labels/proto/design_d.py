"""Design D: side columns (boundary labelling).

Labels stand in a right and a left column outside the circle, each column ordered top to bottom by the y of the
feature's leader stub on the circle.  Leaders (feature -> radial stub -> straight run -> short horizontal stub into
the label) then cannot cross.  Vertical positions come from common.pack_ordered, the stack closest to the wanted y
of every label, so a column simply grows taller than the circle when it must.  Nothing is dropped.

Options (GBDRAW_LABEL_PROTO_OPTS):
  flat    1 = strictly vertical columns, 0 = columns hug a circle of radius X (bend in at top and bottom)  [0]
  mix     0 = wanted y is the stub's own y, 1 = evenly spread over the column height by rank               [0]
  xoff    column distance beyond the stub radius, px                                                       [24]
  hs      length of the horizontal leader stub into the label, px                                         [8]
  cutmax  how far (degrees) from 12 / 6 o'clock a feature may be moved to the other column                [40]
"""

from __future__ import annotations

import math

from .common import Item, Ctx, box_from_left, box_from_right, emit, pack_ordered, polar


def _wrong(side: str, it: Item) -> bool:
    return (math.sin(it.theta) < 0) if side == "R" else (math.sin(it.theta) > 0)


def _stack(side: str, col: list[Item], rs: float, gap: float, mix: float, xrun):
    """Pack a column given top to bottom.  Leaders must stay outside the circle of radius rs, so a label may lie
    below its stub (upper half) only as far as the tangent at the stub allows (and symmetrically in the lower
    half): that is an upper bound on y in the upper half, a lower bound in the lower half.  Wrong-side items
    (feature on the other half; at the column ends) must lie beyond the circle's top/bottom: leader goes up or
    down, then sideways.  `xrun(y)` is the x where the leader's straight run ends at height y."""
    if not col:
        return [], 0.0
    n = len(col)
    own = [-rs * math.cos(i.theta) for i in col]
    sizes = [i.h + gap for i in col]
    tg = list(own)
    if mix > 0 and n > 1:
        height = max(sum(sizes), 2 * rs)
        lin = [(k / (n - 1) - 0.5) * height for k in range(n)]
        tg = [(1 - mix) * a + mix * b for a, b in zip(tg, lin)]
    ys = pack_ordered(tg, sizes)
    sep = [0.5 * (sizes[k] + sizes[k + 1]) for k in range(n - 1)]
    up = [k for k, i in enumerate(col) if math.cos(i.theta) > 1e-6]
    lo = [k for k, i in enumerate(col) if math.cos(i.theta) < -1e-6]
    for k in reversed(up):   # upper half: bound from above, pushing items up
        it = col[k]
        if _wrong(side, it):
            ub = -rs - 4.0 - 0.5 * it.h
        else:
            dx = max(0.0, xrun(own[k], it) - abs(rs * math.sin(it.theta)))
            ub = own[k] + abs(math.tan(it.theta)) * dx
        ys[k] = min(ys[k], ub)
        if k + 1 < n and k + 1 in up:
            ys[k] = min(ys[k], ys[k + 1] - sep[k])
    for k in lo:             # lower half: bound from below, pushing items down
        it = col[k]
        if _wrong(side, it):
            lb = rs + 4.0 + 0.5 * it.h
        else:
            dx = max(0.0, xrun(own[k], it) - abs(rs * math.sin(it.theta)))
            lb = own[k] - abs(math.tan(it.theta)) * dx
        ys[k] = max(ys[k], lb)
        if k > 0 and k - 1 in lo:
            ys[k] = max(ys[k], ys[k - 1] + sep[k - 1])
    cost = sum((y - o) ** 2 for y, o in zip(ys, own))
    return ys, cost


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    if not items:
        return []
    flat = ctx.opt("flat", 0.0)
    mix = ctx.opt("mix", 0.0)
    xoff = ctx.opt("xoff", 24.0)
    hs = ctx.opt("hs", 8.0)
    cutmax = math.radians(ctx.opt("cutmax", 40.0))
    rs = ctx.r_elbow
    X = rs + xoff
    gap = ctx.gap
    # the wrap is placed at 300 degrees so that the right set is one contiguous run
    keys = [i.theta - 2 * math.pi if i.theta > math.radians(300) else i.theta for i in items]
    order = sorted(range(len(items)), key=lambda k: keys[k])
    seq = [items[k] for k in order]
    ks = [keys[k] for k in order]
    n = len(seq)

    def candidates(centre: float, limit: int = 30) -> list[int]:
        cs = [c for c in range(n + 1)
              if (c < n and abs(ks[c] - centre) <= cutmax) or (c > 0 and abs(ks[c - 1] - centre) <= cutmax)]
        cs = sorted(set(cs)) or [sum(1 for k in ks if k < centre)]
        if len(cs) > limit:
            step = len(cs) / limit
            cs = [cs[int(i * step)] for i in range(limit)]
        return cs

    def x0_base(y: float, it: Item) -> float:
        dy = max(0.0, abs(y) - 0.5 * it.h)
        return max(flat * X + (1 - flat) * math.sqrt(max(0.0, X * X - dy * dy)), 6.0)

    def need(it: Item, side: str) -> float:
        return 6.0 if _wrong(side, it) else min(X, abs(rs * math.sin(it.theta)) + hs + 13.0)

    def xrun(y: float, it: Item) -> float:
        return max(x0_base(y, it), need(it, "R" if math.sin(it.theta) >= 0 else "L")) - 3.0 - hs

    best = None
    for c1 in candidates(0.0):
        for c2 in candidates(math.pi):
            if c2 < c1:
                continue
            cost = _stack("R", seq[c1:c2], rs, gap, mix, xrun)[1] + _stack("L", (seq[c2:] + seq[:c1])[::-1], rs, gap, mix, xrun)[1]
            if best is None or cost < best[0]:
                best = (cost, c1, c2)
    _, c1, c2 = best
    out: list[dict] = []
    heights = {}
    for side, srt in (("R", seq[c1:c2]), ("L", (seq[c2:] + seq[:c1])[::-1])):
        ys, _ = _stack(side, srt, rs, gap, mix, xrun)
        if ys:
            heights[side] = round(ys[-1] - ys[0] + srt[0].h, 1)
        x0s = [max(x0_base(y, it), need(it, side)) for it, y in zip(srt, ys)]
        # column x must not decrease towards the circle's top/bottom, so that no leader reaches a label through
        # the text of another one
        for k in range(len(srt)):
            if ys[k] < 0 and k > 0 and ys[k - 1] < 0:
                x0s[k] = max(x0s[k], x0s[k - 1])
        for k in range(len(srt) - 2, -1, -1):
            if ys[k] >= 0 and ys[k + 1] >= 0:
                x0s[k] = max(x0s[k], x0s[k + 1])
        for it, y, x0 in zip(srt, ys, x0s):
            if side == "R":
                box = box_from_left(x0, y, it.w, it.h)
                end, bend, align = (x0 - 3.0, y), (x0 - 3.0 - hs, y), "start"
            else:
                box = box_from_right(-x0, y, it.w, it.h)
                end, bend, align = (-x0 + 3.0, y), (-x0 + 3.0 + hs, y), "end"
            p = polar(rs, it.theta)
            pts = [it.anchor, p]
            if _wrong(side, it):
                pts.append((p[0], y))
            elif (bend[0] > p[0] + 1) if side == "R" else (bend[0] < p[0] - 1):
                pts.append(bend)
            pts.append(end)
            out.append(emit(it, box, pts, align=align, tag=side))
    ctx.notes.update(n_right=c2 - c1, n_left=n - (c2 - c1), cuts=(c1, c2), col_h=heights, X=round(X, 1), flat=flat, mix=mix)
    return out
