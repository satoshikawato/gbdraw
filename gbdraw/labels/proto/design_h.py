"""Design H: ordered free placement on both sides of the ring; study prototype.

Outside shore = design D's two ordered side columns, pulled in against the ring (wall = ring + clear) with a short
run.  Inside shore (only where the inside of the ring is free, from ctx.radial_layout as in design G) = two ordered
columns on the inner side of the ring, text horizontal and extending towards the centre.  Order is kept per shore
and side, so leaders cannot cross; the ring separates the shores.  Labels start outside; the label whose move
inside lowers the total leader length most is moved while the inside packing stays feasible.
Options: clear 3, inside 1, xoff, run, hs, stub, cutmax, cand, rounds (see place()).
"""

from __future__ import annotations

import math

import numpy as np

from .common import (
    Ctx, Item, box_from_left, box_from_right, emit, nearest_point_on_box, pack_ordered, polar, segment_hits_box,
    segment_min_radius, segments_cross,
)
from .design_g import _Zones, _attach, _seg_box_hit


def _plen(pts) -> float:
    return sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))


# ====================================================================== outside shore (design D, pulled in)

class _Out:
    def __init__(self, ctx: Ctx, zones: _Zones, clear: float):
        self.ctx = ctx
        axis = getattr(ctx.radial_layout, "axis", None)       # the axis circle is drawn outside the features
        self.r_out = max(zones.r_out, (axis.radius_px + 0.5 * axis.stroke_width_px) if axis is not None else 0.0)
        self.rs = self.r_out + ctx.opt("stub", 1.0)           # radius where the leader leaves the ring
        self.Rw = self.r_out + clear + ctx.opt("xoff", 0.0)    # no label box closer than this to the centre
        self.hs = ctx.opt("hs", 3.0)
        self.run = ctx.opt("run", 5.0)                          # minimum horizontal run beside the stub
        self.gap = ctx.gap
        self.cutmax = math.radians(ctx.opt("cutmax", 40.0))

    @staticmethod
    def wrong(side: str, it: Item) -> bool:
        return (math.sin(it.theta) < 0) if side == "R" else (math.sin(it.theta) > 0)

    def x0_base(self, y: float, it: Item) -> float:
        dy = max(0.0, abs(y) - 0.5 * it.h)
        return max(math.sqrt(max(0.0, self.Rw * self.Rw - dy * dy)), 6.0)

    def need(self, it: Item, side: str) -> float:
        return 6.0 if self.wrong(side, it) else abs(self.rs * math.sin(it.theta)) + self.hs + 3.0 + self.run

    def xrun(self, y: float, it: Item) -> float:
        return max(self.x0_base(y, it), self.need(it, "R" if math.sin(it.theta) >= 0 else "L")) - 3.0 - self.hs

    def stack(self, side: str, col: list[Item]):
        rs, gap = self.rs, self.gap
        if not col:
            return [], 0.0
        n = len(col)
        own = [-rs * math.cos(i.theta) for i in col]
        sizes = [i.h + gap for i in col]
        ys = pack_ordered(own, sizes)
        sep = [0.5 * (sizes[k] + sizes[k + 1]) for k in range(n - 1)]
        up = [k for k, i in enumerate(col) if math.cos(i.theta) > 1e-6]
        lo = [k for k, i in enumerate(col) if math.cos(i.theta) < -1e-6]
        for k in reversed(up):
            it = col[k]
            if self.wrong(side, it):
                ub = -rs - 4.0 - 0.5 * it.h
            else:
                dx = max(0.0, self.xrun(own[k], it) - abs(rs * math.sin(it.theta)))
                ub = own[k] + abs(math.tan(it.theta)) * dx
            ys[k] = min(ys[k], ub)
            if k + 1 < n and k + 1 in up:
                ys[k] = min(ys[k], ys[k + 1] - sep[k])
        for k in lo:
            it = col[k]
            if self.wrong(side, it):
                lb = rs + 4.0 + 0.5 * it.h
            else:
                dx = max(0.0, self.xrun(own[k], it) - abs(rs * math.sin(it.theta)))
                lb = own[k] - abs(math.tan(it.theta)) * dx
            ys[k] = max(ys[k], lb)
            if k > 0 and k - 1 in lo:
                ys[k] = max(ys[k], ys[k - 1] + sep[k - 1])
        return ys, sum((y - o) ** 2 for y, o in zip(ys, own))

    def layout(self, items: list[Item], thr=None):
        """items: any subset.  Returns (records, thresholds, heights); a record is (item, box, points, align, side)."""
        if not items:
            return [], thr, {}
        keys = [i.theta - 2 * math.pi if i.theta > math.radians(300) else i.theta for i in items]
        order = sorted(range(len(items)), key=lambda k: keys[k])
        seq = [items[k] for k in order]
        ks = [keys[k] for k in order]
        n = len(seq)
        if thr is None:
            def candidates(centre: float, limit: int = 30) -> list[int]:
                cs = [c for c in range(n + 1)
                      if (c < n and abs(ks[c] - centre) <= self.cutmax) or (c > 0 and abs(ks[c - 1] - centre) <= self.cutmax)]
                cs = sorted(set(cs)) or [sum(1 for k in ks if k < centre)]
                if len(cs) > limit:
                    step = len(cs) / limit
                    cs = [cs[int(i * step)] for i in range(limit)]
                return cs

            best = None
            for c1 in candidates(0.0):
                for c2 in candidates(math.pi):
                    if c2 < c1:
                        continue
                    cost = self.stack("R", seq[c1:c2])[1] + self.stack("L", (seq[c2:] + seq[:c1])[::-1])[1]
                    if best is None or cost < best[0]:
                        best = (cost, c1, c2)
            _, c1, c2 = best

            def cut(c: int) -> float:
                if 0 < c < n:
                    return 0.5 * (ks[c - 1] + ks[c])
                return ks[0] - 0.01 if c == 0 else ks[-1] + 0.01

            thr = (cut(c1), cut(c2))
        else:
            c1 = sum(1 for k in ks if k < thr[0])
            c2 = max(c1, sum(1 for k in ks if k < thr[1]))
        recs = []
        heights = {}
        rs, hs = self.rs, self.hs
        for side, srt in (("R", seq[c1:c2]), ("L", (seq[c2:] + seq[:c1])[::-1])):
            ys, _ = self.stack(side, srt)
            if ys:
                heights[side] = round(ys[-1] - ys[0] + srt[0].h, 1)
            x0s = [max(self.x0_base(y, it), self.need(it, side)) for it, y in zip(srt, ys)]
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
                if self.wrong(side, it):
                    pts.append((p[0], y))
                elif (bend[0] > p[0] + 1) if side == "R" else (bend[0] < p[0] - 1):
                    pts.append(bend)
                pts.append(end)
                recs.append((it, box, pts, align, side))
        return recs, thr, heights


# ====================================================================== inside shore

def _toward_axis(box, att, ring, tol=1.0) -> bool:
    """LP-D09 (same test as scripts/lp_metrics.toward_axis): the end of the box away from the leader is closer to the
    ring than the end the leader meets."""
    dist = lambda x, y: abs(math.hypot(x, y) - ring)
    corners = [(x, y) for x in (box[0], box[2]) for y in (box[1], box[3])]
    near = min(corners, key=lambda c: (round(math.dist(c, att), 1), dist(*c)))
    far_x = box[2] if near[0] == box[0] else box[0]
    return min(dist(far_x, box[1]), dist(far_x, box[3])) < dist(*near) - tol


def _ring_of(it: Item) -> float:
    return float(it.src.get("feature_center_radius_px", 0.0)) or math.hypot(*it.anchor)


class _In:
    """Two ordered columns on the inner side of the ring: right (feature on the right half) and left."""

    def __init__(self, ctx: Ctx, zones: _Zones, clear: float):
        self.ctx, self.z, self.clear = ctx, zones, clear
        self.gap = ctx.gap
        self.ok = False
        self.lpad = ctx.opt("lpad", 1.5)
        self.inset = 4.0
        top = zones.inside_hi
        # forbidden radial intervals (boxes), merged; the outermost free gap [A, B] is the only one used
        iv = sorted((a - clear, b + clear) for a, b, _c in zones.annuli if b < zones.r_in - 1e-6 or a < zones.r_in - 1e-6)
        iv = [(a, b) for a, b in iv if not (a >= zones.r_in - clear - 1e-6)]
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        self.B = top
        self.A = 0.0
        for a, b in reversed(merged):
            if b >= top:
                self.B = min(self.B, a)
            else:
                self.A = b
                break
        # forbidden intervals that reach into [0, top]: lower edge of the outermost one is B
        self.ticks = [tuple(t) for t in zones.ticklabels]
        self.slide = ctx.opt("slide", 40.0)
        self.ok = self.B - self.A > 3 * ctx.font_size and self.B > 40.0
        self.cache: dict[float, tuple] = {}

    def x1_at(self, w: float, h: float, y: np.ndarray, right: bool = True) -> np.ndarray:
        """Largest right edge of a box (width w, height h) centred at y that lies in the free gap [A, B] and keeps
        `clear` from the tick labels (it slides towards the centre by up to `slide` px); nan = none."""
        ytop = np.abs(y) + 0.5 * h
        ynear = np.maximum(np.abs(y) - 0.5 * h, 0.0)
        A, B = self.A, self.B
        with np.errstate(invalid="ignore"):
            x1m = np.sqrt(np.where(B > ytop, B * B - ytop * ytop, np.nan))
            low = np.where(ynear < A, np.sqrt(np.maximum(A * A - ynear * ynear, 0.0)), -x1m)
        low = np.maximum(low, 0.5 * self.gap)      # the two columns never meet: right boxes stay at x >= gap/2
        res = np.full(len(y), np.nan)
        T = np.array(self.ticks, dtype=float).reshape(-1, 4)
        c = self.clear
        todo = np.ones(len(y), bool)
        for k in range(int(self.slide / 2.0) + 1):
            x1 = x1m - 2.0 * k
            ok = todo & (x1 - w >= low - 1e-9)
            if len(T) and ok.any():
                bx0, bx1 = (x1 - w, x1) if right else (-x1, -x1 + w)
                hit = (
                    (bx0[:, None] < T[None, :, 2] + c) & (T[None, :, 0] < bx1[:, None] + c)
                    & ((y - 0.5 * h)[:, None] < T[None, :, 3] + c) & (T[None, :, 1] < (y + 0.5 * h)[:, None] + c)
                ).any(axis=1)
                ok &= ~hit
            res = np.where(ok, x1, res)
            todo &= ~ok
            if not todo.any():
                break
        return res

    def box_at(self, it: Item, y: float, side: str):
        x1 = float(self.x1_at(it.w, it.h, np.array([y]), side == "R")[0])
        if math.isnan(x1):
            return None
        return (x1 - it.w, y - 0.5 * it.h, x1, y + 0.5 * it.h) if side == "R" else (-x1, y - 0.5 * it.h, -x1 + it.w, y + 0.5 * it.h)

    def box_ok(self, box) -> bool:
        return True

    def run_of(self, it: Item, y0: float):
        """Contiguous feasible interval of the centre y (1 px grid) nearest to y0 for this item: its box fits the free
        gap, clear of the tick labels, and the straight leader from its stub reaches the box without touching a tick
        label.  None when no y is feasible."""
        key = it.index
        if key not in self.cache:
            grid = np.arange(-self.B, self.B + 1e-9, 1.0)
            x1 = self.x1_at(it.w, it.h, grid, math.sin(it.theta) > 0)
            feas = ~np.isnan(x1)
            if feas.any():
                sx, sy = polar(self.z.r_in, it.theta)
                xx0, xx1 = (x1 - it.w, x1) if math.sin(it.theta) > 0 else (-x1, -x1 + it.w)
                ax, ay = _attach(sx, sy, xx0, grid - 0.5 * it.h, xx1, grid + 0.5 * it.h, self.inset)
                a_xy = np.stack([np.where(feas, ax, 0.0), np.where(feas, ay, 0.0)], axis=1)
                p_xy = np.broadcast_to(np.array([sx, sy]), a_xy.shape)
                if len(self.ticks):
                    hit = _seg_box_hit(p_xy, a_xy, np.array(self.ticks, dtype=float), self.lpad).any(axis=1)
                    feas &= ~hit
                ring = _ring_of(it)
                bx0, bx1, by0, by1 = xx0, xx1, grid - 0.5 * it.h, grid + 0.5 * it.h
                cx = np.stack([bx0, bx0, bx1, bx1], axis=1)
                cy = np.stack([by0, by1, by0, by1], axis=1)
                dd = np.round(np.hypot(cx - a_xy[:, 0:1], cy - a_xy[:, 1:2]), 1)
                dr = np.abs(np.hypot(cx, cy) - ring)
                nk = np.argmin(dd * 1.0e4 + dr, axis=1)
                near_is_left = nk < 2
                far_x = np.where(near_is_left, bx1, bx0)
                dfar = np.minimum(np.abs(np.hypot(far_x, by0) - ring), np.abs(np.hypot(far_x, by1) - ring))
                dnear = np.take_along_axis(dr, nk[:, None], axis=1)[:, 0]
                feas &= ~(dfar < dnear - 1.0)
                d = a_xy - p_xy
                l2 = np.maximum((d * d).sum(axis=1), 1e-12)
                t = np.clip(-(p_xy * d).sum(axis=1) / l2, 0.0, 1.0)
                near = p_xy + t[:, None] * d
                feas &= np.hypot(near[:, 0], near[:, 1]) >= self.A - 0.5
            self.cache[key] = (grid, feas)
        grid, feas = self.cache[key]
        idx = np.nonzero(feas)[0]
        if not len(idx):
            return None
        i = min(max(int(round(y0 - grid[0])), 0), len(grid) - 1)
        if not feas[i]:
            i = int(idx[np.argmin(np.abs(idx - i))])
        lo = i
        while lo > 0 and feas[lo - 1]:
            lo -= 1
        hi = i
        while hi < len(grid) - 1 and feas[hi + 1]:
            hi += 1
        return grid[lo], grid[hi]

    why: dict = {}

    def fail(self, key):
        self.why[key] = self.why.get(key, 0) + 1
        return None

    def column(self, col: list[Item], side: str):
        """col given top to bottom; mirror for the left side.  Returns list of (item, box, pts, align) or None."""
        if not col:
            return []
        z = self.z
        n = len(col)
        h = max(i.h for i in col)
        stubs = [polar(z.r_in, i.theta) for i in col]
        own = [s[1] for s in stubs]
        lo, hi = [], []
        for it, y0 in zip(col, own):
            r = self.run_of(it, y0)
            if r is None:
                return self.fail("norun")
            lo.append(r[0])
            hi.append(r[1])
        sizes = [i.h + self.gap for i in col]
        ys = pack_ordered(own, sizes)
        sep = [0.5 * (sizes[k] + sizes[k + 1]) for k in range(n - 1)]
        for k in range(n):
            ys[k] = max(ys[k], lo[k], (ys[k - 1] + sep[k - 1]) if k else -1e9)
        for k in range(n - 1, -1, -1):
            ys[k] = min(ys[k], hi[k], (ys[k + 1] - sep[k]) if k < n - 1 else 1e9)
        for k in range(n):
            if ys[k] < lo[k] - 1e-6 or ys[k] > hi[k] + 1e-6 or (k and ys[k] - ys[k - 1] < sep[k - 1] - 1e-6):
                return self.fail("pack")
        res = []
        for it, y, s in zip(col, ys, stubs):
            box = self.box_at(it, y, side)
            if box is None:
                return self.fail("x1nan")
            ax, ay = _attach(s[0], s[1], box[0], box[1], box[2], box[3], self.inset)
            att = (float(ax), float(ay))
            if _toward_axis(box, att, _ring_of(it)):
                return self.fail("toward_axis")
            if segment_min_radius(s, att) < self.A - 0.5:
                return self.fail("leaderrad")
            res.append((it, box, [it.anchor, s, att], "end" if side == "R" else "start"))
        return res

    def layout(self, sel: list[Item]):
        """sel: items moved inside (any order).  Returns records or None when infeasible."""
        if not sel:
            return []
        right = sorted((i for i in sel if math.sin(i.theta) > 0), key=lambda i: i.theta)
        left = sorted((i for i in sel if math.sin(i.theta) <= 0), key=lambda i: -(i.theta - 2 * math.pi if i.theta > 0 else i.theta))
        left = sorted((i for i in sel if math.sin(i.theta) <= 0), key=lambda i: (i.theta - math.pi), reverse=True)
        r = self.column(right, "R")
        if r is None:
            return None
        l = self.column(left, "L")
        if l is None:
            return None
        recs = r + l
        if not self.valid(recs):
            return None
        return recs

    def valid(self, recs) -> bool:
        g, pad = self.gap, self.lpad
        boxes = [r[1] for r in recs]
        segs = [list(zip(r[2], r[2][1:])) for r in recs]
        n = len(recs)
        for i in range(n):
            a = boxes[i]
            for j in range(i + 1, n):
                b = boxes[j]
                if a[0] < b[2] + g and b[0] < a[2] + g and a[1] < b[3] + g and b[1] < a[3] + g:
                    self.fail("v_overlap")
                    return False
        for i in range(n):
            for j in range(n):
                if i != j and any(segment_hits_box(p, q, boxes[j], pad) for p, q in segs[i][1:]):
                    self.fail("v_thrutext")
                    return False
            for t in self.ticks:
                if any(segment_hits_box(p, q, t, pad) for p, q in segs[i][1:]):
                    self.fail("v_thrutick")
                    return False
            for j in range(i + 1, n):
                if any(segments_cross(p, q, r, s) for p, q in segs[i] for r, s in segs[j]):
                    self.fail("v_cross")
                    return False
        return True


# ====================================================================== the move search and the drawing

def place(items: list[Item], ctx: Ctx) -> list[dict]:
    clear = ctx.opt("clear", ctx.gap)
    zones = _Zones(items, ctx, clear)
    out = _Out(ctx, zones, clear)
    inn = _In(ctx, zones, clear)
    use_in = bool(int(ctx.opt("inside", 1))) and inn.ok
    recs, thr, heights = out.layout(items)
    cur_in: list[Item] = []
    in_recs: list[tuple] = []
    total = sum(_plen(r[2]) for r in recs)
    moves = 0
    if use_in:
        k_cand = int(ctx.opt("cand", 200))
        cutin = math.sin(math.radians(ctx.opt("cutin", 0.0)))
        for _ in range(int(ctx.opt("rounds", 400))):
            ranked = sorted(recs, key=lambda r: -_plen(r[2]))
            best = None
            tried = 0
            for it, _b, pts, _a, _s in ranked:
                if tried >= k_cand:
                    break
                if abs(math.sin(it.theta)) < cutin:
                    continue
                run = inn.run_of(it, polar(zones.r_in, it.theta)[1])
                if run is None:
                    continue
                tried += 1
                trial_in = cur_in + [it]
                ir = inn.layout(trial_in)
                if ir is None:
                    continue
                rest = [r[0] for r in recs if r[0] is not it]
                orr, _t, _h = out.layout(rest, thr)
                cost = sum(_plen(r[2]) for r in orr) + sum(_plen(r[2]) for r in ir)
                if cost < total - 1e-6 and (best is None or cost < best[0]):
                    best = (cost, it, orr, ir)
            if best is None:
                break
            total, it, recs, in_recs = best
            cur_in.append(it)
            moves += 1
    if moves:   # the side columns' cut points were fixed during the search: look for better ones now
        orr, thr2, h2 = out.layout([r[0] for r in recs])
        c2 = sum(_plen(r[2]) for r in orr) + sum(_plen(r[2]) for r in in_recs)
        if c2 < total - 1e-6:
            recs, total, heights = orr, c2, h2
    result = []
    hide = bool(int(ctx.opt("hide", 0)))      # evidence only: same layout, invisible text and no leaders (see check_h.py)
    for it, box, pts, align, side in recs:
        result.append(emit(it, box, None if hide else pts, align=align, tag=side))
    for it, box, pts, align in in_recs:
        result.append(emit(it, box, None if hide else pts, align=align, tag="in"))
    if hide:
        for label in result:
            label["text_fill"] = "none"

    ctx.notes.update(
        inside=len(in_recs), outside=len(recs), moves=moves, col_h=heights, total=round(total, 1),
        axis=str(getattr(ctx.radial_layout, 'axis', None))[:200], outer_content=getattr(ctx.radial_layout, 'outer_content_radius_px', None),
        r_in=round(zones.r_in, 1), r_out=round(zones.r_out, 1), gapAB=(round(inn.A, 1), round(inn.B, 1)), has_inside=use_in, why=dict(inn.why),
        tick_r=sorted({round(math.hypot(0.5 * (t[0] + t[2]), 0.5 * (t[1] + t[3]))) for t in zones.ticklabels})[:3],
        annuli=[(round(a, 1), round(b, 1)) for a, b, _c in zones.annuli],
    )
    return result
