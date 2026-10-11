"""Design F: minimum total leader length (study prototype).

Leader = today's polyline: feature anchor -> radial stub to the elbow E (src middle_x/y) -> straight run to the
label.  "Length" is the whole polyline, as lp_metrics measures it; the stub is constant, so the objective is the
sum (sq / max) of the run |E - attach|.

Default (new, provable): see _Ell.  free=1: no ellipse (heuristic, may cross; comparison only).
free: no ellipse.  A box of horizontal text may stand anywhere outside the ring (radius >= r_elbow +
1 + clear) with clearance `clear` (default ctx.gap) to every other label box and to every other leader; the run
keeps outside the feature ring and crosses no other leader.  The attach point is the middle of the box's left,
right, bottom or top edge.  Solver (not exact): candidate positions are the points of a grid around E (3 px up
to 300 px, 6 px up to 700, 12 px beyond) in order of distance; labels are placed greedily (longest text first)
at the nearest feasible candidate; then relocation sweeps and one-for-one ejection moves (a label takes a
nearer place that blocks exactly one other label, which is re-placed; accepted when the objective falls) run
until no move helps or `tmax` seconds.  The result is compared with the lower bound "every label at its
nearest clear place, ignoring the other labels" (ctx.notes lb_sum).

ellipse=1: today's ellipse and look, solved exactly per side by dynamic programming (see _Ell).

Options: obj=sum|sq|max [sum]  clear [ctx.gap]  ellipse [0]  order=w|a [w]  tmax [6]  att [2]
ellipse: band0 [1.0]  bandmax [3.0]  dy [1.0]  cut [60]  maxcut [36]  rf [1]
"""

from __future__ import annotations

import math
import os
import time

import numpy as np

from .common import Ctx, Item, box_from_left, box_from_right, boxes_overlap, emit

INF = float("inf")


# ------------------------------------------------------------------ vectorised geometry

def _seg_box(px, py, qx, qy, x0, y0, x1, y1, pad):
    """Broadcast Liang-Barsky: does segment pq pass through the box grown by pad?"""
    x0, y0, x1, y1 = x0 - pad, y0 - pad, x1 + pad, y1 + pad
    dx, dy = qx - px, qy - py
    t0 = np.zeros(np.broadcast(dx, dy, x0, y0, px, py).shape)
    t1 = np.ones_like(t0)
    ok = np.ones(t0.shape, dtype=bool)
    for d, lo, hi, st in ((dx, x0, x1, px), (dy, y0, y1, py)):
        par = np.abs(d) < 1e-12
        sd = np.where(par, 1.0, d)
        ta, tb = (lo - st) / sd, (hi - st) / sd
        t0 = np.where(par, t0, np.maximum(t0, np.minimum(ta, tb)))
        t1 = np.where(par, t1, np.minimum(t1, np.maximum(ta, tb)))
        ok = ok & np.where(par, (st > lo) & (st < hi), True)
    return ok & (t0 < t1)


def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _seg_cross(px, py, qx, qy, rx, ry, sx, sy):
    d1, d2 = _orient(rx, ry, sx, sy, px, py), _orient(rx, ry, sx, sy, qx, qy)
    d3, d4 = _orient(px, py, qx, qy, rx, ry), _orient(px, py, qx, qy, sx, sy)
    ok = ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))
    for d in (d1, d2, d3, d4):
        ok = ok & (np.abs(d) > 1e-9)
    return ok


def _chord_min_radius(ex, ey, ax, ay):
    dx, dy = ax - ex, ay - ey
    l2 = dx * dx + dy * dy
    t = np.clip(np.where(l2 > 0, -(ex * dx + ey * dy) / np.where(l2 > 0, l2, 1.0), 0.0), 0.0, 1.0)
    return np.hypot(ex + t * dx, ey + t * dy)


# ------------------------------------------------------------------ free placement

_OFF = None


def _offsets() -> np.ndarray:
    global _OFF
    if _OFF is None:
        parts = []
        for lo, hi, st in ((0, 300, 3), (300, 700, 6), (700, 1600, 12)):
            n = int(hi / st) + 1
            a = np.arange(-n, n + 1) * float(st)
            gx, gy = np.meshgrid(a, a)
            d = np.hypot(gx, gy)
            m = (d >= lo) & (d < hi)
            parts.append(np.stack([gx[m], gy[m], d[m]], 1))
        arr = np.concatenate(parts)
        _OFF = arr[np.argsort(arr[:, 2], kind="stable")]
    return _OFF


def _box_for(t, px, py, w, h, att):
    """Box of width w, height h whose edge middle (t: 0 left, 1 right, 2 bottom, 3 top) is att px from (px, py)."""
    if t == 0:
        return px + att, py - 0.5 * h, px + att + w, py + 0.5 * h
    if t == 1:
        return px - att - w, py - 0.5 * h, px - att, py + 0.5 * h
    if t == 2:
        return px - 0.5 * w, py - att - h, px + 0.5 * w, py - att
    return px - 0.5 * w, py + att, px + 0.5 * w, py + att + h


class _Free:
    def __init__(self, items: list[Item], ctx: Ctx) -> None:
        self.items, self.ctx = items, ctx
        n = len(items)
        self.clear = ctx.opt("clear", ctx.gap)
        self.att = ctx.opt("att", 2.0)
        self.r_lim = ctx.r_elbow + 1.0 + self.clear
        self.rf = max(i.feat_outer_r for i in items) - 0.5
        self.E = np.array([[float(i.src["middle_x"]), float(i.src["middle_y"])] for i in items])
        self.base = np.array([math.dist(i.anchor, (float(i.src["middle_x"]), float(i.src["middle_y"]))) for i in items])
        self.B = np.zeros((n, 4))
        self.SEG = np.zeros((n, 2, 4))
        self.BB = np.zeros((n, 4))
        self.alive = np.zeros(n, dtype=bool)
        self.par: list = [None] * n          # (t, px, py, d, crossings with leaders placed before)
        self.xpen = ctx.opt("xpen", 120.0)
        self.pb, self.pl = ctx.opt("pb", 300.0), ctx.opt("pl", 200.0)
        self.cur = 0
        self.fallback = False
        self.obj = ctx.options.get("obj", "sum")

    def length(self, i: int) -> float:
        return float(self.base[i] + self.par[i][3])

    # --- store
    def put(self, i: int, par) -> None:
        t, px, py = par[:3]
        it = self.items[i]
        box = _box_for(t, px, py, it.w, it.h, self.att)
        ex, ey = self.E[i]
        self.B[i] = box
        self.SEG[i, 0] = (it.anchor[0], it.anchor[1], ex, ey)
        self.SEG[i, 1] = (ex, ey, px, py)
        pad = self.clear + 1.0
        self.BB[i] = (min(box[0], px, ex, it.anchor[0]) - pad, min(box[1], py, ey, it.anchor[1]) - pad,
                      max(box[2], px, ex, it.anchor[0]) + pad, max(box[3], py, ey, it.anchor[1]) + pad)
        self.alive[i] = True
        self.par[i] = par

    def drop(self, i: int) -> None:
        self.alive[i] = False

    # --- candidates
    def _chunk(self, i: int, off: np.ndarray):
        it = self.items[i]
        ex, ey = self.E[i]
        px, py, d = ex + off[:, 0], ey + off[:, 1], off[:, 2]
        out = []
        for t in range(4):
            x0, y0, x1, y1 = _box_for(t, px, py, it.w, it.h, self.att)
            rmin = np.hypot(np.maximum(np.maximum(x0, 0.0), -x1), np.maximum(np.maximum(y0, 0.0), -y1))
            ok = rmin >= self.r_lim
            ok &= _chord_min_radius(ex, ey, px, py) >= self.rf
            ok &= ~_seg_box(ex, ey, px, py, x0, y0, x1, y1, -0.5)
            idx = np.nonzero(ok)[0]
            if idx.size:
                out.append((np.full(idx.size, t), px[idx], py[idx], d[idx], x0[idx], y0[idx], x1[idx], y1[idx]))
        if not out:
            return None
        c = [np.concatenate([o[k] for o in out]) for k in range(8)]
        order = np.argsort(c[3], kind="stable")
        return [a[order] for a in c]     # t, px, py, d, x0, y0, x1, y1

    def _conflicts(self, i: int, cand):
        """(K, M) conflict matrix of the candidates against alive placed items; returns (matrix, item ids)."""
        t, px, py, d, x0, y0, x1, y1 = cand
        lo_x, hi_x = min(x0.min(), px.min()) - 1, max(x1.max(), px.max()) + 1
        lo_y, hi_y = min(y0.min(), py.min()) - 1, max(y1.max(), py.max()) + 1
        m = self.alive & (self.BB[:, 2] > lo_x) & (self.BB[:, 0] < hi_x) & (self.BB[:, 3] > lo_y) & (self.BB[:, 1] < hi_y)
        m[i] = False
        ids = np.nonzero(m)[0]
        if not ids.size:
            return None, ids
        c = self.clear
        ex, ey = self.E[i]
        Bx0, By0, Bx1, By1 = (self.B[ids, k][None, :] for k in range(4))
        X0, Y0, X1, Y1, PX, PY = (a[:, None] for a in (x0, y0, x1, y1, px, py))
        boxes = (X0 < Bx1 + c) & (Bx0 < X1 + c) & (Y0 < By1 + c) & (By0 < Y1 + c)
        legs = _seg_box(ex, ey, PX, PY, Bx0, By0, Bx1, By1, c)
        cross = np.zeros(boxes.shape, dtype=bool)
        for s in (0, 1):
            sx0, sy0, sx1, sy1 = (self.SEG[ids, s, k][None, :] for k in range(4))
            legs = legs | _seg_box(sx0, sy0, sx1, sy1, X0, Y0, X1, Y1, c)
            cross |= _seg_cross(ex, ey, PX, PY, sx0, sy0, sx1, sy1)
        return (boxes, legs, cross), ids

    def penalty(self, cand) -> np.ndarray:
        """Per-candidate penalty: pb per clearance violation with a box, pl per leader-through-box (either way),
        xpen per leader crossing."""
        res, _ = self._conflicts(self.cur, cand)
        k = cand[0].size
        if res is None:
            return np.zeros(k)
        boxes, legs, cross = res
        return self.pb * boxes.sum(1) + self.pl * legs.sum(1) + self.xpen * cross.sum(1)

    def scan(self, i: int, maxcost: float = INF):
        """Cheapest candidate (cost = |E - attach| + penalties), scanning offsets in order of distance.  The penalties
        only matter when the nearest places are blocked; the final layout is accepted only at zero penalty."""
        off = _offsets()
        best = None
        self.cur = i
        for s in range(0, len(off), 192):
            part = off[s:s + 192]
            if self.base[i] + part[0, 2] >= (maxcost if best is None else min(maxcost, best[0])):
                break
            cand = self._chunk(i, part)
            if cand is None:
                continue
            d = cand[3]
            cost = self.base[i] + d + self.penalty(cand)
            k = int(np.argmin(cost))
            if cost[k] < maxcost and (best is None or cost[k] < best[0]):
                best = (float(cost[k]), (int(cand[0][k]), float(cand[1][k]), float(cand[2][k]), float(d[k])))
        return best

    def fresh_cost(self, i: int) -> float:
        """Cost of label i at its stored place against the others as they are now."""
        t, px, py, d = self.par[i]
        it = self.items[i]
        x0, y0, x1, y1 = _box_for(t, px, py, it.w, it.h, self.att)
        cand = [np.array([t]), np.array([px]), np.array([py]), np.array([d]), np.array([x0]), np.array([y0]), np.array([x1]), np.array([y1])]
        self.cur = i
        return float(self.base[i] + d + self.penalty(cand)[0])

    def violations(self, i: int) -> int:
        t, px, py, d = self.par[i]
        it = self.items[i]
        x0, y0, x1, y1 = _box_for(t, px, py, it.w, it.h, self.att)
        cand = [np.array([t]), np.array([px]), np.array([py]), np.array([d]), np.array([x0]), np.array([y0]), np.array([x1]), np.array([y1])]
        res, _ = self._conflicts(i, cand)
        return 0 if res is None else int(sum(m.sum() for m in res))

    def blockers(self, i: int) -> list[int]:
        t, px, py, d = self.par[i]
        it = self.items[i]
        x0, y0, x1, y1 = _box_for(t, px, py, it.w, it.h, self.att)
        cand = [np.array([t]), np.array([px]), np.array([py]), np.array([d]), np.array([x0]), np.array([y0]), np.array([x1]), np.array([y1])]
        res, ids = self._conflicts(i, cand)
        if res is None:
            return []
        m = res[0][0] | res[1][0] | res[2][0]
        return [int(j) for j in ids[m]]

    def group_score(self, group) -> tuple:
        cs = [self.fresh_cost(k) for k in group]
        if self.obj == "sq":
            return (sum(c * c for c in cs),)
        if self.obj == "max":
            return (max(cs), sum(cs))
        return (sum(cs),)

    def ripple(self, i: int) -> bool:
        """Re-place a violating label first and its blockers after it; keep the result when the summed cost falls."""
        group = [i] + self.blockers(i)
        old = {k: self.par[k] for k in group}
        before = self.group_score(group)
        for k in group:
            self.drop(k)
        placed = []
        for k in group:
            f = self.scan(k)
            self.put(k, f[1])
            placed.append(k)
        after = self.group_score(group)
        if after < (before[0] - 0.5,) + before[1:]:
            return True
        for k in group:
            self.put(k, old[k])
        return False

    def solve(self) -> None:
        items, ctx = self.items, self.ctx
        n = len(items)
        t_end = time.time() + ctx.opt("tmax", 25.0)
        order = sorted(range(n), key=lambda k: (-items[k].w, k))
        if ctx.options.get("order", "a") == "a":      # clockwise sweep from the widest angular gap
            th = [it.theta for it in items]
            srt = sorted(range(n), key=lambda k: th[k])
            gaps = [((th[srt[(j + 1) % n]] - th[srt[j]]) % (2 * math.pi), j) for j in range(n)]
            start = (max(gaps)[1] + 1) % n
            order = srt[start:] + srt[:start]
        lbs = [float(self.base[i] + self.scan(i)[1][3]) for i in range(n)]    # nearest clear place, others ignored
        lb = sum(lbs)
        seed = _ellipse_solution(items, ctx, self.obj) if ctx.opt("seed", 1.0) > 0 else None
        if seed is not None:
            ctx.notes.update(seed[1])
            pos = {it.index: k for k, it in enumerate(items)}
            for it, side, xl, y, _, _ax in seed[0]:
                k = pos[it.index]
                px = xl - self.att if side > 0 else -xl + self.att
                self.put(k, (0 if side > 0 else 1, px, y, float(math.hypot(px - self.E[k][0], y - self.E[k][1]))))
        else:
            for i in order:
                self.put(i, self.scan(i)[1])
        ctx.notes.update(lb_sum=round(lb, 1), seed_sum=round(sum(self.length(k) for k in range(n)), 1))
        moves: list = []           # (version, x, y) of every box that moved (old and new centre)
        last = [-1] * n            # version at the last scan of each label
        rounds = 0
        scans = 0
        while time.time() < t_end and rounds < 30:
            rounds += 1
            changed = False
            ver = len(moves)
            cost = {k: self.fresh_cost(k) for k in range(n)}
            for i in sorted(range(n), key=lambda k: -cost[k]):
                if time.time() > t_end:
                    break
                cur = self.fresh_cost(i)
                if cur <= lbs[i] + 1.0:
                    continue
                reach = cur - self.base[i] + 360.0
                if last[i] >= 0 and not any(v >= last[i] and math.hypot(x - self.E[i][0], y - self.E[i][1]) < reach
                                            for v, x, y in moves[last[i]:]):
                    continue
                last[i] = len(moves)
                scans += 1
                old = self.par[i]
                ob = self.B[i].copy()
                self.drop(i)
                f = self.scan(i, maxcost=cur - 0.5)
                if f:
                    self.put(i, f[1])
                    nb = self.B[i]
                    moves.append((len(moves), 0.5 * (ob[0] + ob[2]), 0.5 * (ob[1] + ob[3])))
                    moves.append((len(moves), 0.5 * (nb[0] + nb[2]), 0.5 * (nb[1] + nb[3])))
                    changed = True
                else:
                    self.alive[i] = True
            for i in range(n):
                if time.time() > t_end:
                    break
                if self.violations(i) and self.ripple(i):
                    changed = True
                    moves.append((len(moves), self.E[i][0], self.E[i][1]))
                    last = [-1] * n
            if not changed:
                break
        v = [self.violations(k) for k in range(n)]
        overlap = any(boxes_overlap(tuple(self.B[a]), tuple(self.B[b]), 0.0)
                      for a in range(n) for b in range(a + 1, n))
        if overlap and seed is not None:     # safeguard: never ship an overprint; return to the exact ellipse layout
            self.fallback = True
            ctx.notes.update(fallback_to_ellipse=True)
            return
        ctx.notes.update(rounds=rounds, scans=scans, final_sum=round(sum(self.length(k) for k in range(n)), 1),
                         violating=int(sum(1 for x in v if x)))

    def emit_all(self) -> list[dict]:
        out = []
        for i, it in enumerate(self.items):
            t, px, py = self.par[i][:3]
            box = _box_for(t, px, py, it.w, it.h, self.att)
            e = (float(self.E[i][0]), float(self.E[i][1]))
            align = {0: "start", 1: "end", 2: "middle", 3: "middle"}[t]
            out.append(emit(it, box, [it.anchor, e, (px, py)], align=align, tag="free"))
        return out


def _place_free(items: list[Item], ctx: Ctx) -> list[dict]:
    t0 = time.time()
    st = _Free(items, ctx)
    st.solve()
    ctx.notes.update(obj=st.obj, clear=st.clear, solve_s=round(time.time() - t0, 1))
    return _place_ellipse(items, ctx) if st.fallback else st.emit_all()


# ------------------------------------------------------------------ ellipse variant (exact per side)

def _r_e(ctx: Ctx) -> float:
    """Radius of the elbow circle shared by all leaders: today's (ctx.r_elbow), or with tight=1 just outside the ring."""
    if ctx.opt("tight", 0.0) > 0:
        return ctx.opt("_ring", ctx.r_elbow - 4.0) + 4.0
    return ctx.r_elbow


class _Grid:
    def __init__(self, X: float, Y: float, dy: float) -> None:
        self.X, self.Y, self.dy = X, Y, dy
        self.n = int(math.floor(2 * Y / dy)) + 1
        self.y = -Y + dy * np.arange(self.n)
        self.idx = np.arange(self.n)

    def axf(self, y):
        return self.X * np.sqrt(np.clip(1.0 - (y / self.Y) ** 2, 0.0, None))

    def x_clear(self, h: float) -> np.ndarray:
        """Left edge of a box of height h centred on each row that lies wholly outside the ellipse."""
        top, bot = self.y - 0.5 * h, self.y + 0.5 * h
        return np.where((top < 0) & (bot > 0), self.X, np.maximum(self.axf(top), self.axf(bot)))


class _Ell:
    """Rows on a half ellipse (right half; the left is mirrored).  Labels of one side keep their clockwise order and
    stand in disjoint rows (pitch (h_a + h_b)/2 + gap); a label's left edge is the ellipse x of the row's inner box
    edge, so the box lies outside the ellipse.  dp[g] = cost[g] + min_{g' <= g - k} dp[g'] : exact for the 1 px grid."""

    def __init__(self, ctx: Ctx, band: float, dy: float, rf: float, obj: str, thr: float = INF) -> None:
        self.ctx, self.band, self.rf, self.obj, self.thr = ctx, band, rf, obj, thr
        self.grid = _Grid(ctx.arc_x_radius * band, ctx.arc_y_radius * band, dy)
        self.cache: dict = {}

    def cost(self, it: Item, side: int):
        key = (it.index, side)
        if key in self.cache:
            return self.cache[key]
        g = self.grid
        r_e = _r_e(self.ctx)
        k = r_e / math.hypot(*it.anchor)
        ex, ey = it.anchor[0] * k * side, it.anchor[1] * k         # elbow: on the common circle r_e, on the anchor's ray
        xl = g.x_clear(it.h)
        ax = g.axf(g.y)
        base = math.dist(it.anchor, (it.anchor[0] * k, ey))
        length = base + np.hypot(ax - ex, g.y - ey) + (xl - ax)    # stub + run to the ellipse + horizontal stub to the text
        ok = length <= self.thr
        # the run stays outside the circle r_e (its radius never decreases): this and the common ellipse make
        # the leaders of one record non-crossing (see ctx notes / report)
        ok = ok & (ax * ex + g.y * ey >= r_e * r_e)
        cost = np.where(ok, length if self.obj != "sq" else length * length, INF)
        self.cache[key] = (cost, length, xl)
        return self.cache[key]

    def run(self, seq, keep: bool = False):
        g = self.grid
        out: list[float] = []
        parents: list = []
        dp = None
        for pos, (it, side) in enumerate(seq):
            cost = self.cost(it, side)[0]
            if dp is None:
                dp = cost.copy()
            else:
                k = int(math.ceil((0.5 * (seq[pos - 1][0].h + it.h) + self.ctx.opt("clear", self.ctx.gap)) / g.dy - 1e-9))
                pm = np.minimum.accumulate(dp)
                f = g.idx - k
                valid = f >= 0
                fi = np.where(valid, f, 0)
                if keep:
                    arg = np.maximum.accumulate(np.where(dp == pm, g.idx, -1))
                    parents.append(arg[fi])
                dp = cost + np.where(valid, pm[fi], INF)
            best = float(dp.min())
            out.append(best)
            if best == INF:
                out.extend([INF] * (len(seq) - pos - 1))
                break
        return out, parents, dp

    def backtrack(self, parents, dp, n: int) -> list[int]:
        path = [int(np.argmin(dp))]
        for pos in range(n - 1, 0, -1):
            path.append(int(parents[pos - 1][path[-1]]))
        return path[::-1]


def _cuts(ks: list[float], centre: float, cutmax: float, limit: int) -> list[int]:
    n = len(ks)
    cs = [c for c in range(n + 1) if (c < n and abs(ks[c] - centre) <= cutmax) or (c > 0 and abs(ks[c - 1] - centre) <= cutmax)]
    cs = sorted(set(cs)) or [sum(1 for k in ks if k < centre)]
    if len(cs) > limit:
        step = len(cs) / limit
        cs = [cs[int(i * step)] for i in range(limit)]
    return cs


def _solve_cuts(seq, ks, prob: _Ell, cutmax: float, limit: int):
    n = len(seq)
    c2s = _cuts(ks, math.pi, cutmax, limit)
    best = None
    for c1 in _cuts(ks, 0.0, cutmax, limit):
        rc = prob.run([(seq[j], 1) for j in range(c1, n)])[0]
        back = list(range(c1 - 1, -1, -1)) + list(range(n - 1, c1 - 1, -1))
        lc = prob.run([(seq[j], -1) for j in back])[0]
        for c2 in c2s:
            if c2 < c1:
                continue
            nl = n - c2 + c1
            tot = (rc[c2 - c1 - 1] if c2 > c1 else 0.0) + (lc[nl - 1] if nl > 0 else 0.0)
            if tot < INF and (best is None or tot < best[0] - 1e-9):
                best = (tot, c1, c2)
    return best


def _materialise(seq, prob: _Ell, c1: int, c2: int):
    n = len(seq)
    g = prob.grid
    out = []
    runs = [(1, [seq[j] for j in range(c1, c2)]),
            (-1, [seq[j] for j in list(range(c1 - 1, -1, -1)) + list(range(n - 1, c2 - 1, -1))])]
    for side, its in runs:
        if not its:
            continue
        s = [(it, side) for it in its]
        _, parents, dp = prob.run(s, keep=True)
        for it, gi in zip(its, prob.backtrack(parents, dp, len(s))):
            _, length, xl = prob.cost(it, side)
            out.append((it, side, float(xl[gi]), float(g.y[gi]), float(length[gi]), float(g.axf(g.y[gi]))))
    return out


def _ellipse_solution(items: list[Item], ctx: Ctx, obj: str):
    """[(item, side, left x, y, length)] of the exact ellipse solution, plus notes; None when infeasible."""
    dy = ctx.opt("dy", 1.0)
    cutmax = math.radians(ctx.opt("cut", 60.0))
    limit = int(ctx.opt("maxcut", 36))
    rf = (max(i.feat_outer_r for i in items) - 0.5) if ctx.opt("rf", 1.0) > 0 else 0.0
    ang = [math.atan2(i.anchor[0], -i.anchor[1]) % (2 * math.pi) for i in items]
    keys = [a - 2 * math.pi if a > math.radians(300) else a for a in ang]
    order = sorted(range(len(items)), key=lambda k: keys[k])
    seq = [items[k] for k in order]
    ks = [keys[k] for k in order]
    bandmax = ctx.opt("bandmax", 3.0)
    band = ctx.opt("band0", 0.0)
    if band <= 0:     # tightest ellipse that clears the ring, rounded up to a 0.05 step
        band = math.ceil((_r_e(ctx) + 8.0 + ctx.opt('clear', ctx.gap)) / min(ctx.arc_x_radius, ctx.arc_y_radius) / 0.05) * 0.05
    t0 = time.time()
    tried = []
    sol = None
    extra = int(ctx.opt("extra", 4))
    while band <= bandmax + 1e-9:
        p_ = _Ell(ctx, band, dy, rf, "sum" if obj == "max" else obj)
        s_ = _solve_cuts(seq, ks, p_, cutmax, limit)
        tried.append(round(band, 2))
        if s_ is not None and (sol is None or s_[0] < sol[0] - 1e-9):
            sol, prob, best_band = s_, p_, band
        if sol is not None:
            if extra <= 0:
                break
            extra -= 1
        band = round(band + 0.05, 4)
    band = best_band if sol is not None else band
    if sol is None:
        ctx.notes.update(infeasible=True, band_tried=tried)
        return None
    if obj == "max":
        mat = _materialise(seq, prob, sol[1], sol[2])
        lo, hi = 0.0, max(m[4] for m in mat)
        best_sol, best_prob = sol, prob
        while hi - lo > 1.0:
            mid = 0.5 * (lo + hi)
            p2 = _Ell(ctx, band, dy, rf, "sum", thr=mid)
            s2 = _solve_cuts(seq, ks, p2, cutmax, max(4, limit // 3))
            if s2 is None:
                lo = mid
            else:
                hi, best_sol, best_prob = mid, s2, p2
        p2 = _Ell(ctx, band, dy, rf, "sum", thr=hi + 0.5)
        s2 = _solve_cuts(seq, ks, p2, cutmax, limit)
        sol, prob = (s2, p2) if s2 else (best_sol, best_prob)
    n_right = sol[2] - sol[1]
    notes = dict(ellipse_band=round(band, 2), band_tried=tried, n_right=n_right, n_left=len(seq) - n_right,
                 ellipse_s=round(time.time() - t0, 1))
    return _materialise(seq, prob, sol[1], sol[2]), notes


def _place_ellipse(items: list[Item], ctx: Ctx) -> list[dict]:
    obj = ctx.options.get("obj", "sum")
    res = _ellipse_solution(items, ctx, obj)
    if res is None:
        return []
    out = []
    for it, side, xl, y, _, ax in res[0]:
        k = _r_e(ctx) / math.hypot(*it.anchor)
        e = (it.anchor[0] * k, it.anchor[1] * k)
        if side > 0:
            out.append(emit(it, box_from_left(xl, y, it.w, it.h), [it.anchor, e, (ax, y), (xl, y)], align="start", tag="R"))
        else:
            out.append(emit(it, box_from_right(-xl, y, it.w, it.h), [it.anchor, e, (-ax, y), (-xl, y)], align="end", tag="L"))
    b = res[1]["ellipse_band"]
    ctx.notes.update(obj=obj, ellipse=1, band=b, **res[1], r_elbow=round(_r_e(ctx), 1),
                     ellipse_xy=(round(ctx.arc_x_radius * b), round(ctx.arc_y_radius * b)))
    return out


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    if not items:
        return []
    ctx.options["_ring"] = str(max(i.feat_outer_r for i in items))
    if ctx.opt("free", 0.0) > 0:       # heuristic 2-D placement; may cross leaders (comparison only)
        return _place_free(items, ctx)
    if ctx.opt("ellipse", 0.0) > 0:    # today's ellipse (band 1.0, grown only when needed), exact
        ctx.options.setdefault("band0", "1.0")
        ctx.options.setdefault("extra", "0")
        ctx.options.setdefault("tight", "0")
    else:                              # default: elbow circle just outside the ring, tightest ellipse that fits
        ctx.options.setdefault("tight", "1")
    return _place_ellipse(items, ctx)
