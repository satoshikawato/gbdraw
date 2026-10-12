"""Design L: design I (H2) with the leader ending anywhere on its label box, and blocks of rows at 12 and 6 o'clock;
study prototype (LP-D11 point 3, LP-D12 point 3, LP-D13).

1. Perimeter rule (`perim` = 1, LP-D13): a leader may end at any point of the edge of its label box (kept `lgap` from
   the text), not only at the text end.  The candidates of design I are kept; for every box position the leader also
   ends at the nearest point of the box (shortest leader) and where the feature's own ray meets the box (straight
   leader), and boxes that cover the feature (centred above or below it) are tried.  The cost of a spot is the leader
   length plus `wt` px per degree of turn.  LP-D09 is read as its intent: no point of the box may be more than `backf`
   text heights closer to the ring than the end of the leader.  `perim` = 0 restores H2 (the leader meets the text end,
   LP-D09 as before).
2. Blocks (`blocks` = 1): design I lays the labels out as usual; the labels that stand wholly above (below) the ring
   are then laid out again as one block: seated one at a time, outer features first, each at the spot of least
   cost = leader length + turn + (`pen` per box and `hw` per px if the leader passes behind a box) + `wh` per px of rise, at
   any height above its own feature, with its box anywhere that a leader of at most `maxrun` px sideways and a turn of
   at most `lean` degrees reaches (the box need not lie beside its feature: a row can hold several labels).  `over` = 1:
   a leader may pass behind another label's box under K's limits (`maxb`, `maxc`, `maxh`, `vis`, `vis0`); `over` = 0:
   a leader keeps `lpad` from every other box.  The block is kept only where it lowers the outer edge of the whole outside layout;
   otherwise the layout of design I stands (nothing dropped, every other rule as in H2).
"""

from __future__ import annotations

import math
import time

import numpy as np

from .common import Ctx, Item, emit
from .design_g import _seg_box_hit, _seg_cross_matrix
from .design_i import _check, _Eng as _EngI, _Geo, _Lab as _LabI, _Par as _ParI, _plen
from .design_k import _chords, _ends, _finish, _over_stats, _Placed as _PlacedK

FAILS: list = []

FR = np.array([0.0, 0.12, 0.25, 0.4, 0.5, 0.6, 0.75, 0.88, 1.0])


class _Par(_ParI):
    def __init__(self, ctx: Ctx):
        super().__init__(ctx)
        self.perim = bool(int(ctx.opt("perim", 1)))       # a leader may end anywhere on the box edge (LP-D13)
        self.backf = ctx.opt("backf", 0.5)                # most a box point may lie nearer the ring than the leader end (x height)
        self.wt = ctx.opt("wt", 0.35)                     # px per degree of turn in the cost of a spot
        self.blocks = bool(int(ctx.opt("blocks", 1)))     # re-lay the labels wholly above/below the ring as a block
        self.over = bool(int(ctx.opt("over", 1)))         # in a block a leader may pass behind another label's box
        self.lean = min(ctx.opt("lean", 70.0), self.cone)  # largest turn of a block leader (degrees)
        self.lean_cos = math.cos(math.radians(self.lean))
        self.maxrun = ctx.opt("maxrun", 140.0)            # sideways run of a block leader (px)
        self.maxb = int(ctx.opt("maxb", 3))
        self.maxc = ctx.opt("maxc", 50.0)
        self.maxh = ctx.opt("maxh", 0.7)
        self.vis = ctx.opt("vis", 14.0)
        self.vis0 = ctx.opt("vis0", 6.0)
        self.pen = ctx.opt("pen", 25.0)
        self.hw = ctx.opt("hw", 0.3)
        self.wh = ctx.opt("wh", 1.2)                      # cost per px of rise of a block label
        self.cap = ctx.opt("cap", 520.0)                  # a block reaches at most this far beyond the ring
        self.step = ctx.opt("step", 3.0)                  # row step of a block
        self.gstep = ctx.opt("gstep", 3.0)                # step of the box position along a row
        self.ordr = int(ctx.opt("ord", 0))                # seating order of a block: 0 outer features first, 1 pole first
        self.rayw = ctx.opt("rayw", 0.0)                  # the kept way out leans this many degrees away from the pole axis
        self.flip = bool(int(ctx.opt("flip", 1)))         # in a block the text may run towards the pole axis (back up to flipb px)
        self.flipb = ctx.opt("flipb", 45.0)
        self.tmax = ctx.opt("tmax", 6.0)                  # seconds after which no further block variant or sweep is started
        self.repair = ctx.opt("repair", 6)                # attempts to seat a label that found no spot (it goes first next time)
        self.sweeps = ctx.opt("sweeps", 2)                # compaction sweeps over a finished block
        self.rays = bool(int(ctx.opt("rays", 1)))         # keep the way out of every label still to be seated
        self.kap = ctx.opt("kap", 0.0)                    # 0: try several sets of labels per pole (see _blocks), else this one


# ====================================================================== design I with the perimeter rule

class _Lab(_LabI):
    """A label in a frame (see design I). With `perim` more spots: the leader meets the nearest point of the box edge or the
    point where the feature's own ray meets it, and boxes that cover the feature."""

    def cands(self, yc, geo, par):
        Q, YA, GA, base = super().cands(yc, geo, par)
        if not par.perim:
            return Q, YA, GA, base
        yc = yc[:, None]
        y0, y1 = yc - 0.5 * self.h, yc + 0.5 * self.h
        lg, w = par.lgap, self.w
        cover = self.Sq - FR[None, :] * w + 0.0 * yc
        q0 = np.concatenate([Q, cover], axis=1)
        y0, y1 = y0 + 0.0 * q0, y1 + 0.0 * q0
        x0, x1, ya0, ya1 = q0 - lg, q0 + w + lg, y0 - lg, y1 + lg
        nq_, ay_ = np.clip(self.Sq, x0, x1), np.clip(self.Sy, ya0, ya1)            # nearest point
        with np.errstate(divide="ignore", invalid="ignore"):
            def slab(lo, hi, s, d):
                if abs(d) < 1e-9:
                    inside = (s > lo) & (s < hi)
                    return np.where(inside, -np.inf, np.inf), np.where(inside, np.inf, -np.inf)
                a, b = (lo - s) / d, (hi - s) / d
                return np.minimum(a, b), np.maximum(a, b)
            tq0, tq1 = slab(x0, x1, self.Sq, self.nq)
            ty0, ty1 = slab(ya0, ya1, self.Sy, self.ny)
        t_in, t_out = np.maximum(tq0, ty0), np.minimum(tq1, ty1)
        hit = (t_in <= t_out) & (t_in > 0) & np.isfinite(t_in)
        rq = np.where(hit, self.Sq + np.where(hit, t_in, 0.0) * self.nq, nq_)       # where the feature's own ray meets it
        ry = np.where(hit, self.Sy + np.where(hit, t_in, 0.0) * self.ny, ay_)
        Q2 = np.concatenate([Q, q0, q0], axis=1)
        YA2 = np.concatenate([YA, ay_, ry], axis=1)
        GA2 = np.concatenate([GA, q0 - nq_, q0 - rq], axis=1)
        return Q2, YA2, GA2, base

    def static(self, yc, q0, ya, ga, geo, par, cone_cos):
        if not par.perim:
            return super().static(yc, q0, ya, ga, geo, par, cone_cos)
        h, w = self.h, self.w
        y0, y1 = yc - 0.5 * h, yc + 0.5 * h
        aq = q0 - ga
        dq, dy = aq - self.Sq, ya - self.Sy
        ln = np.hypot(dq, dy)
        cosang = (dq * self.nq + dy * self.ny) / np.maximum(ln, 1e-9)
        ok = (cosang >= cone_cos - 1e-9) | (ln < par.short)
        # (a leader that ends inside the box, or on the wrong side of it, is never produced: it ends on the box edge)
        p = geo.padw
        nx, nyv = np.clip(0.0, q0 - p, q0 + w + p), np.clip(0.0, y0 - p, y1 + p)
        r_att = np.hypot(aq, ya)
        if self.shore == "out":
            ok &= np.hypot(nx, nyv) >= geo.R
            rmin = np.hypot(np.clip(0.0, q0, q0 + w), np.clip(0.0, y0, y1))
            back = r_att - rmin
        else:
            fx = np.maximum(np.abs(q0 - p), np.abs(q0 + w + p))
            fy = np.maximum(np.abs(y0 - p), np.abs(y1 + p))
            ok &= (np.hypot(fx, fy) <= geo.Bb) & (q0 + w <= -0.5 * par.gap) & (np.hypot(nx, nyv) >= geo.A_in - p)
            back = np.hypot(np.maximum(np.abs(q0), np.abs(q0 + w)), np.maximum(np.abs(y0), np.abs(y1))) - r_att
        ok &= back <= par.backf * h + 1e-6                 # LP-D09 read as its intent (LP-D13)
        if geo.has[self.shore]:
            x0, x1 = (q0, q0 + w) if self.sx > 0 else (-(q0 + w), -q0)
            r = geo.padr
            ok &= ~geo.blocked(x0 - r, y0 - r + 0.0 * q0, x1 + r, y1 + r + 0.0 * q0)
        turn = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
        return ok, ln + par.wt * np.where(ln < par.short, 0.0, turn)


class _Eng(_EngI):
    def lab(self, it: Item, shore: str, sx: int) -> _Lab:
        key = (id(it), shore, sx)
        if key not in self.cache:
            self.cache[key] = _Lab(it, shore, sx, self.geo, self.par)
        return self.cache[key]


# ====================================================================== blocks of rows at 12 and 6 o'clock

class _Block:
    """The labels of one pole seated one at a time: any height, any box position within reach of the leader."""

    def __init__(self, geo: _Geo, par: _Par):
        self.geo, self.par = geo, par

    def ray(self, it: Item, s: float):
        """The way out kept free for a label that is not seated yet: straight along the pole direction, inside the cone."""
        ux, uy = math.sin(it.theta), -math.cos(it.theta)
        lim = math.radians(self.par.lean - 1.5)
        an, av = math.atan2(uy, ux), math.atan2(s, 0.0) + (math.copysign(math.radians(self.par.rayw), ux) * (-s))
        turn = (av - an + math.pi) % (2.0 * math.pi) - math.pi
        a = an + max(-lim, min(lim, turn))
        S = (self.geo.rs_out * ux, self.geo.rs_out * uy)
        return S, (S[0] + 4000.0 * math.cos(a), S[1] + 4000.0 * math.sin(a))

    def seat(self, it: Item, pl: _PlacedK, s: float, pend):
        geo, par = self.geo, self.par
        ux, uy = math.sin(it.theta), -math.cos(it.theta)
        S = np.array([geo.rs_out * ux, geo.rs_out * uy])
        w, h = it.w, it.h
        pb, s0, s1 = pl.arrays()
        best = None
        k = 0
        lim = geo.R + par.cap
        while True:
            rise = k * par.step
            yc = S[1] + s * rise
            if abs(yc) - 0.5 * h > lim or (best is not None and par.wh * rise >= best[0]):
                break
            k += 1
            res = self._row(it, S, (ux, uy), w, h, yc, rise, pb, s0, s1, pl, pend)
            if res is not None and (best is None or res[0] < best[0] - 1e-9):
                best = res
        return best

    def _row(self, it, S, u, w, h, yc, rise, pb, s0, s1, pl, pend):
        geo, par = self.geo, self.par
        y0, y1 = yc - 0.5 * h, yc + 0.5 * h
        lg, gap, pad = par.lgap, par.gap, par.lpad
        xs = np.arange(S[0] - w - par.maxrun, S[0] + par.maxrun + 1e-6, par.gstep)
        if len(pb):
            m = (pb[:, 3] > y0 - gap) & (pb[:, 1] < y1 + gap)
            if m.any():
                ev = np.concatenate([pb[m, 2] + gap + 1e-3, pb[m, 0] - gap - w - 1e-3])
                xs = np.unique(np.round(np.concatenate([xs, ev[(ev > S[0] - w - par.maxrun) & (ev < S[0] + par.maxrun)]]), 3))
        X0 = xs
        X1 = X0 + w
        p = geo.padw
        nx, ny = np.clip(0.0, X0 - p, X1 + p), np.clip(0.0, y0 - p, y1 + p)
        ok = np.hypot(nx, ny) >= geo.R
        if geo.has["out"]:
            r = geo.padr
            ok &= ~geo.blocked(X0 - r, np.full_like(X0, y0 - r), X1 + r, np.full_like(X0, y1 + r))
        if len(pb):
            g = gap - 1e-6
            ok &= ~((X0[:, None] < pb[None, :, 2] + g) & (pb[None, :, 0] < X1[:, None] + g)
                    & (y0 < pb[None, :, 3] + g) & (pb[None, :, 1] < y1 + g)).any(axis=1)
        idx = np.nonzero(ok)[0]
        if not len(idx):
            return None
        X0, X1 = X0[idx], X1[idx]
        n = len(X0)
        # the two ends of the leader: nearest point of the box edge, and where the feature's own ray meets it
        rx0, rx1, ry0, ry1 = X0 - lg, X1 + lg, y0 - lg, y1 + lg
        nxp, nyp = np.clip(S[0], rx0, rx1), np.full(n, float(np.clip(S[1], ry0, ry1)))
        with np.errstate(divide="ignore", invalid="ignore"):
            def slab(lo, hi, st, d):
                if abs(d) < 1e-9:
                    inside = (st > lo) & (st < hi)
                    return np.where(inside, -np.inf, np.inf), np.where(inside, np.inf, -np.inf)
                a, b = (lo - st) / d, (hi - st) / d
                return np.minimum(a, b), np.maximum(a, b)
            tq0, tq1 = slab(rx0, rx1, S[0], u[0])
            ty0, ty1 = slab(np.full(n, ry0), np.full(n, ry1), S[1], u[1])
        t_in, t_out = np.maximum(tq0, ty0), np.minimum(tq1, ty1)
        hit = (t_in <= t_out) & (t_in > 0) & np.isfinite(t_in)
        rxp = np.where(hit, S[0] + np.where(hit, t_in, 0.0) * u[0], nxp)
        ryp = np.where(hit, S[1] + np.where(hit, t_in, 0.0) * u[1], nyp)
        bi = np.concatenate([np.arange(n), np.arange(n)])
        if not par.perim:                                   # H2's rule: the leader meets the text end (left or right edge)
            ey = np.full(n, float(np.clip(S[1], y0, y1)))
            nxp, nyp, rxp, ryp = X0 - lg, ey, X1 + lg, ey
        A = np.stack([np.concatenate([nxp, rxp]), np.concatenate([nyp, ryp])], axis=1)
        d = A - S[None, :]
        ln = np.hypot(d[:, 0], d[:, 1])
        cosang = (d[:, 0] * u[0] + d[:, 1] * u[1]) / np.maximum(ln, 1e-9)
        good = ((cosang >= par.lean_cos - 1e-9) | (ln < par.short)) & (np.abs(d[:, 0]) <= par.maxrun + 1e-9)
        if not par.perim and not par.flip:                  # text runs away from the pole axis
            left = np.concatenate([np.ones(n, dtype=bool), np.zeros(n, dtype=bool)])      # the left edge is the first n ends
            good &= left if S[0] >= 0 else ~left
        rmin = np.hypot(np.clip(0.0, X0, X1), np.clip(0.0, y0, y1))[bi]
        good &= (np.hypot(A[:, 0], A[:, 1]) - rmin) <= (par.flipb if par.flip else (par.backf * h if par.perim else 0.7)) + 1e-6
        if geo.has["out"] and good.any():
            sel = np.nonzero(good)[0]
            good[sel[geo.seg_blocked(S, A[sel])]] = False
        cost = ln + par.wt * np.where(ln < par.short, 0.0, np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))) + par.wh * rise
        boxes = np.stack([X0, np.full(n, y0), X1, np.full(n, y1)], axis=1)[bi]
        stub = np.broadcast_to(S, A.shape)
        over = par.over
        lo, hi = min(S[1], y0) - pad - 1.0, max(S[1], y1) + pad + 1.0
        if len(pb) and good.any():
            sel = np.nonzero(good)[0]
            pbn = pb[(pb[:, 3] > lo - gap) & (pb[:, 1] < hi + gap)]
            if len(pbn):
                if over:
                    ch = _chords(stub[sel], A[sel], pbn, pad)
                    hitm = ch > 0.3
                    nb, hd = hitm.sum(axis=1), ch.sum(axis=1)
                    ok2 = (nb <= par.maxb) & (ch.max(axis=1) <= par.maxc) & (hd <= par.maxh * ln[sel])
                    pas = np.nonzero(ok2 & (nb > 0))[0]
                    if len(pas):
                        e0, e1, e2, e3 = _ends(stub[sel][pas], A[sel][pas], par.vis0, par.vis)
                        ok2[pas] &= ~_seg_box_hit(e0, e1, pbn, pad).any(axis=1) & ~_seg_box_hit(e2, e3, pbn, pad).any(axis=1)
                    good[sel[~ok2]] = False
                    cost[sel] += par.pen * nb + par.hw * hd
                else:
                    good[sel[_seg_box_hit(stub[sel], A[sel], pbn, pad).any(axis=1)]] = False
        if len(s0) and good.any():
            m = (np.minimum(s0[:, 1], s1[:, 1]) < hi) & (np.maximum(s0[:, 1], s1[:, 1]) > lo)
            sl = np.nonzero(m)[0]
            if len(sl):
                a0, a1 = s0[sl], s1[sl]
                sel = np.nonzero(good)[0]
                bad = _seg_cross_matrix(stub[sel], A[sel], a0, a1).any(axis=1)
                if over:
                    ch = _chords(a0, a1, boxes[sel], pad)
                    hitm = ch > 0.3
                    cnt = np.array([pl.cnt[i] for i in sl], dtype=float)[:, None]
                    hid = np.array([pl.hid[i] for i in sl], dtype=float)[:, None]
                    ln_l = np.maximum(np.hypot(a1[:, 0] - a0[:, 0], a1[:, 1] - a0[:, 1]), 1e-9)[:, None]
                    ok3 = (cnt + hitm <= par.maxb) & (ch <= par.maxc) & (hid + ch <= par.maxh * ln_l)
                    rows = np.nonzero(hitm.any(axis=1))[0]
                    if len(rows):
                        e0, e1, e2, e3 = _ends(a0[rows], a1[rows], par.vis0, par.vis)
                        ok3[rows] &= ~(_seg_box_hit(e0, e1, boxes[sel], pad) | _seg_box_hit(e2, e3, boxes[sel], pad))
                    bad |= ~ok3.all(axis=0)
                    cost[sel] += par.pen * hitm.sum(axis=0) + par.hw * ch.sum(axis=0)
                else:
                    bad |= _seg_box_hit(a0, a1, boxes[sel], pad).any(axis=0)
                good[sel[bad]] = False
        rp, rq = pend
        if len(rp) and good.any():                          # keep the way out of every label still to be seated
            m = (np.minimum(rp[:, 1], rq[:, 1]) < hi) & (np.maximum(rp[:, 1], rq[:, 1]) > lo)
            rp, rq = rp[m], rq[m]
            if len(rp):
                sel = np.nonzero(good)[0]
                bad = _seg_cross_matrix(stub[sel], A[sel], rp, rq).any(axis=1)
                hb = _seg_box_hit(rp, rq, boxes[sel], pad)
                if over:                                    # a box may lie over the ray, but not over its first px
                    e0, e1, _e2, _e3 = _ends(rp, rq, par.vis0 + pad + 1.0, 1.0)
                    bad |= _seg_box_hit(e0, e1, boxes[sel], pad).any(axis=0)
                    cost[sel] += 0.5 * par.pen * hb.sum(axis=0)
                else:
                    bad |= hb.any(axis=0)
                good[sel[bad]] = False
        if not good.any():
            return None
        sel = np.nonzero(good)[0]
        j = sel[int(np.argmin(cost[sel] + 1e-6 * np.abs(boxes[sel, 0] - S[0])))]
        return float(cost[j]), tuple(float(v) for v in boxes[j]), (float(A[j][0]), float(A[j][1])), S


def _extent(recs, upper: bool) -> float:
    if not recs:
        return 0.0
    return max(-r[1][1] for r in recs) if upper else max(r[1][3] for r in recs)


def _block_pole(recs, geo: _Geo, par: _Par, upper: bool, kap: float, deadline: float, hard: float):
    """recs with the labels wholly beyond kap * R of the pole laid out again as a block: (recs, n labels) or None."""
    R = geo.R
    sel = [r for r in recs if (r[1][3] < -R * kap if upper else r[1][1] > R * kap)]
    if len(sel) < 3:
        return None
    ids = {id(r[0]) for r in sel}
    fixed = [r for r in recs if id(r[0]) not in ids]
    key = (lambda r: -abs(r[0].anchor[0])) if par.ordr == 0 else (lambda r: abs(r[0].anchor[0]))
    blk = _Block(geo, par)
    s = -1.0 if upper else 1.0
    base_order = sorted(sel, key=key)
    front: list = []
    new: list = []
    for _attempt in range(int(par.repair) + 1):               # a label with no spot is seated first in the next attempt
        order = front + [r for r in base_order if not any(r is f for f in front)]
        rays = [blk.ray(r[0], s) for r in order]
        rp_all = np.array([r[0] for r in rays], dtype=float).reshape(-1, 2)
        rq_all = np.array([r[1] for r in rays], dtype=float).reshape(-1, 2)
        pl = _PlacedK(par)
        for _it, box, pts, _al, _tag in fixed:
            pl.add(box, pts[1], pts[2])
        new = []
        failed = None
        for i, (it, _box, _pts, _al, tag) in enumerate(order):
            live = np.arange(len(order)) > i
            if time.perf_counter() > hard:                   # too slow on this map: the layout of design I stands
                return None
            res = blk.seat(it, pl, s, (rp_all[live], rq_all[live]) if par.rays else (rp_all[:0], rq_all[:0]))
            if res is None:
                failed = order[i]
                FAILS.append((it.text, len(new), len(sel)))
                break
            _cost, box, att, S = res
            pl.add(box, (float(S[0]), float(S[1])), att)
            new.append((it, box, [it.anchor, (float(S[0]), float(S[1])), att], "start", tag))
        if failed is None:
            break
        front = [failed] + front
    else:
        return None
    for _sweep in range(int(par.sweeps)):                   # every label again, highest first, with all the others in place
        moved = False
        if time.perf_counter() > deadline:
            break
        for idx in sorted(range(len(new)), key=lambda k: -(-new[k][1][1] if upper else new[k][1][3])):
            it, box, _pts, _al, tag = new[idx]
            pl2 = _PlacedK(par)
            for _it, bx, pts, _a, _t in fixed + new[:idx] + new[idx + 1:]:
                pl2.add(bx, pts[1], pts[2])
            res = blk.seat(it, pl2, s, (rp_all[:0], rq_all[:0]))
            if res is None:
                continue
            _cost, nbox, att, S = res
            old, cur = (-box[1], -nbox[1]) if upper else (box[3], nbox[3])
            if cur < old - 0.5:
                new[idx] = (it, nbox, [it.anchor, (float(S[0]), float(S[1])), att], "start", tag)
                moved = True
        if not moved:
            break
    return fixed + new, len(sel)


def _blocks(recs, geo: _Geo, par: _Par, notes: dict):
    """Lay out the labels wholly above and below the ring as blocks wherever that lowers the outer edge."""
    out = list(recs)
    info = {}
    for name, upper in (("top", True), ("bottom", False)):
        t0 = time.perf_counter()
        base = _extent(out, upper)
        kaps = [par.kap] if par.kap > 0 else [1.0, 0.7]
        best = None
        trials: list = []
        for kap in kaps:
            if best is not None and time.perf_counter() - t0 > par.tmax:
                break
            res = _block_pole(out, geo, par, upper, kap, t0 + par.tmax, t0 + 2.0 * par.tmax)
            if res is None:
                trials.append((kap, "failed", FAILS[-1] if FAILS else None))
                continue
            new, n = res
            chk = _check(new)
            trials.append((kap, n, round(_extent(new, upper)), round(_plen_sum(new)), chk["overlap"], chk["cross"]))
            if chk["overlap"] or chk["cross"] or (chk["gap"] is not None and chk["gap"] < par.gap - 0.06):
                continue
            score = (round(_extent(new, upper)), _plen_sum(new))
            if best is None or score < best[0]:
                best = (score, new, n, kap)
        info[name] = {"base_extent": round(base, 1), "trials": trials}
        if best is not None and best[0][0] < base - 1.0:
            out = best[1]
            info[name].update(extent=best[0][0], moved=best[2], kap=best[3])
    notes["blocks"] = info
    return out


def _plen_sum(recs) -> float:
    return sum(_plen(r[2]) for r in recs)


# ====================================================================== statistics

def _stats(recs, geo: _Geo, par: _Par) -> dict:
    R = geo.R
    above = sum(1 for r in recs if r[1][3] < -R)
    below = sum(1 for r in recs if r[1][1] > R)
    backs, runs = [], []
    corner = corner_straight = 0
    lg = par.lgap
    for it, box, pts, _al, _tag in recs:
        end = pts[-1]
        end_r = math.hypot(*end)
        ring = float(it.src.get("feature_center_radius_px", 0.0)) or math.hypot(*pts[0])
        if end_r >= ring:
            b = end_r - math.hypot(min(max(0.0, box[0]), box[2]), min(max(0.0, box[1]), box[3]))
        else:
            b = max(math.hypot(x, y) for x in (box[0], box[2]) for y in (box[1], box[3])) - end_r
        backs.append(max(0.0, b))
        runs.append(abs(end[0] - pts[1][0]))
        if min(abs(end[0] - (box[0] - lg)), abs(end[0] - (box[2] + lg))) < 1.5 and min(abs(end[1] - (box[1] - lg)), abs(end[1] - (box[3] + lg))) < 1.5:
            corner += 1
            d1, d2 = (pts[1][0] - pts[0][0], pts[1][1] - pts[0][1]), (end[0] - pts[1][0], end[1] - pts[1][1])
            n1, n2 = math.hypot(*d1), math.hypot(*d2)
            if n1 > 1e-6 and n2 > 1e-6 and (d1[0] * d2[0] + d1[1] * d2[1]) / (n1 * n2) > math.cos(math.radians(3.0)):
                corner_straight += 1
    return {"above": above, "below": below, "back_over1": sum(1 for b in backs if b > 1.0), "back_max": round(max(backs, default=0.0), 1),
            "run_max": round(max(runs, default=0.0), 1), "corner_ends": corner, "corner_straight": corner_straight}


def _depth(recs, upper: bool, geo: _Geo) -> int:
    """Most boxes any vertical line crosses among the labels wholly above (below) the ring: the rows of the block."""
    R = geo.R
    boxes = [r[1] for r in recs if (r[1][3] < -R if upper else r[1][1] > R)]
    if not boxes:
        return 0
    xs = np.arange(min(b[0] for b in boxes), max(b[2] for b in boxes), 2.0)
    return int(max(sum(1 for b in boxes if b[0] <= x <= b[2]) for x in xs))


# ====================================================================== entry point

def place(items: list[Item], ctx: Ctx) -> list[dict]:
    par = _Par(ctx)
    inner_req = list(getattr(ctx, "inner_items", []) or [])
    geo = _Geo(items + inner_req, ctx, par)
    eng = _Eng(ctx, geo, par)
    out_items = list(items)
    in_items: list[Item] = []
    in_recs: list = []
    overflow: list[Item] = []
    if inner_req:
        in_items = list(inner_req) if geo.zone_ok else []
        overflow = [] if geo.zone_ok else list(inner_req)
        while in_items:
            in_recs, failed = eng.inside(in_items)
            if in_recs is not None:
                break
            in_items = [i for i in in_items if i is not failed]
            overflow.append(failed)
            in_recs = []
        out_items += overflow

    def total(recs) -> float:
        return sum(_plen(r[2]) for r in recs)

    def outside_best(sel: list[Item]):
        if not sel:
            return [], None, {}
        cuts = eng.split(sel, None)[1]
        recs, stats = eng.outside(sel, cuts)
        cuts_far = eng.split(sel, None, far=True)[1]
        if cuts_far != cuts:
            alt, alt_stats = eng.outside(sel, cuts_far)
            if total(alt) < total(recs) - 1e-6:
                return alt, cuts_far, alt_stats
        return recs, cuts, stats

    out_recs, cuts, stats = outside_best(out_items)
    moves = 0
    if not inner_req and bool(int(ctx.opt("inside", 1))) and geo.zone_ok and out_recs:
        sel_in: list[Item] = []
        cur_in: list = []
        rest, cur_out, cur_stats = list(out_items), out_recs, stats
        rejected: set = set()
        while True:
            length = {id(r[0]): _plen(r[2]) for r in cur_out}
            cands = []
            for it in rest:
                if id(it) in rejected:
                    continue
                L = eng.lab(it, "in", -1 if math.sin(it.theta) > 0 else 1)
                gain = -1.0 if L.t is None else length[id(it)] - (math.dist(it.anchor, L.S) + L.base)
                if gain > 1.0:
                    cands.append((gain, it))
                else:
                    rejected.add(id(it))
            cands.sort(key=lambda c: (-c[0], c[1].index))
            took: list[Item] = []
            for _gain, it in cands[:10]:
                trial, _failed = eng.inside(sel_in + took + [it])
                if trial is not None and total(trial) - total(cur_in) < length[id(it)] - 0.5:
                    took.append(it)
                    cur_in = trial
                else:
                    rejected.add(id(it))
            if not took:
                if len(cands) <= 10:
                    break
                continue
            sel_in += took
            rest = [it for it in rest if not any(it is t for t in took)]
            cur_out, cur_stats = eng.outside(rest, cuts)
        if sel_in:
            again = outside_best(rest)
            if total(again[0]) < total(cur_out) - 1e-6:
                cur_out, cuts_in, cur_stats = again
            else:
                cuts_in = cuts
            if total(cur_out) + total(cur_in) < total(out_recs) - 1e-6:
                out_recs, in_recs, cuts, stats, moves = cur_out, cur_in, cuts_in, cur_stats, len(sel_in)
    notes: dict = {}
    notes["stats_before_blocks"] = _stats(out_recs, geo, par)
    if par.blocks and out_recs:
        out_recs = _blocks(out_recs, geo, par, notes)
    recs = out_recs + in_recs
    hide = bool(int(ctx.opt("hide", 0)))
    result = [emit(it, box, None if hide else pts, align=align, tag=tag) for it, box, pts, align, tag in recs]
    if hide:
        for label in result:
            label["text_fill"] = "none"
    _finish(result, recs, ctx, par)
    ctx.notes.update(
        inside=len(in_recs), outside=len(out_recs), asked_inside=len(inner_req), overflow=[i.text for i in overflow], moves=moves,
        total=round(total(recs), 1), R=round(geo.R, 2), B=round(geo.B, 2), Bb=round(geo.Bb, 2), A_in=round(geo.A_in, 2),
        zone_ok=geo.zone_ok, objects=geo.objects, raster=geo.has, tick_font=geo.tick_font,
        cuts=[round(math.degrees(c), 1) for c in cuts] if cuts else None, stats=stats, check=_check(recs), behind=_over_stats(recs, 0.0),
        lstats=_stats(recs, geo, par), depth={"top": _depth(recs, True, geo), "bottom": _depth(recs, False, geo)}, **notes,
    )
    return result
