"""Design J (shown to the Owner as H3): design I with the turn of the leaders as a second objective; study prototype.

Design I (H2) revised after LP-D11 point 1: a leader costs its length plus `tw` px per degree it turns at the ring edge, and a
label may sit exactly on the ray of its feature (one straight line).  The rest is the text of design I:

Revision of design H after LP-D10.  A label is a horizontal box whose text starts at the leader and runs away from
the ring (LP-D09).  A leader is anchor -> ring edge -> end of the text: one elbow at the ring edge that turns at most
`cone` degrees (70: every bend is 110 degrees or more).  A label's first choice is the spot with the shortest leader
at clearance `clear` from everything drawn: beside the feature at 3 and 9 o'clock, above or below it at 12 and 6.

Per shore (outside, inside) the labels form four stacks (right/left x upper/lower half) in the clockwise order of
their features; a cut at each pole says which column the labels near the pole join.  A stack is seated from the
equator towards its pole (inside: from the pole towards the equator).  Each label starts at its own best height and
takes the spot with the smallest leader run + push along the stack, tested exactly against every label seated before
it (box gap, leader to box clearance, leader crossings) and against the other objects (ticks, tick labels, tracks,
title; from ctx.radial_layout).  No spot may block the exit ray of a label that is still to be seated (straight
along the stack, or the edge of its cone near the equator), so a label can always be seated on its own ray: nothing
is dropped and no leader has to turn sharply.
Labels asked for inside the ring (`--labels both`, ctx.inner_items) are placed inside by the same rule; with
`inside=1` an outside label moves inside when that shortens the total leader length (as in H).
Options: clear 3, inside 1, cone 70, lgap 2.5, lpad 2.4, stub 1, edge 1, fl 12, cutmax 30, cutn 14, hide 0, tw 0.8, ray 1, plat 2 (see _Par).
"""

from __future__ import annotations

import math

import numpy as np

from ...auto_sizes import circular_tick_intervals
from ...core.numeric import scaled_tick_text
from ...core.text import calculate_bbox_dimensions
from .common import Ctx, Item, emit, pack_ordered
from .design_g import _seg_box_hit, _seg_cross_matrix

OFFS = np.array([0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 30, 36, 44, 52, 62, 75, 90, 110, 135], dtype=float)
LEAN = np.tan(np.radians([2, 5, 8, 11, 14, 17, 20, 24, 28, 33, 39, 46, 54, 62]))
NEAR = [0.0, 1.0, 2.0, 4.0, 6.0, 9.0]
PUSH = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
while PUSH[-1] < 4000.0:
    PUSH.append(round(PUSH[-1] * 1.12 + 1.0, 1))


def _plen(pts) -> float:
    return sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))


def _turn(pts, min_seg: float = 3.0) -> float:
    """Total turn of a polyline in degrees (the metric's rule: pieces under 3 px are dropped)."""
    pts = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= min_seg]
    turn = 0.0
    for a, b, c in zip(pts, pts[1:], pts[2:]):
        u, v = (b[0] - a[0], b[1] - a[1]), (c[0] - b[0], c[1] - b[1])
        nu, nv = math.hypot(*u), math.hypot(*v)
        if nu >= min_seg and nv >= min_seg:
            turn += math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv)))))
    return turn


class _Par:
    def __init__(self, ctx: Ctx):
        self.gap = float(ctx.gap) + 0.05                 # between two label boxes
        self.clear = ctx.opt("clear", ctx.gap)           # label box to every other object
        self.cone = ctx.opt("cone", 70.0)                # largest turn of a leader at the ring edge (degrees)
        self.cone_cos = math.cos(math.radians(self.cone))
        self.lgap = ctx.opt("lgap", 2.5)                 # leader end to its own text box
        self.lpad = min(ctx.opt("lpad", 2.4), self.lgap - 0.1)   # a leader passing another label's box (below lgap:
        #                                                  a label on its own exit ray never blocks a neighbour's ray)
        self.stub = ctx.opt("stub", 1.0)                 # elbow distance beyond the ring edge
        self.edge = ctx.opt("edge", 1.0)                 # allowance for the feature outline stroke
        self.fl = ctx.opt("fl", 12.0)                    # a label may sit towards the equator of its feature while a
        #                                                  spot within fl px of the ring exists there
        self.cutmax = math.radians(ctx.opt("cutmax", 30.0))   # the cut at a pole lies within this angle of the pole
        self.cutn = ctx.opt("cutn", 14)                  # ... and among this many label gaps
        self.short = 2.5                                 # a run shorter than this has no direction
        self.tw = ctx.opt("tw", 0.8)                      # px of leader cost per degree of turn at the ring edge
        self.plat = ctx.opt("plat", 2.0)                 # px of leader cost a label may give up to make room for the labels after it
        self.ray = int(ctx.opt("ray", 1))                # 1: candidates exactly on the feature's own ray (straight leaders)
        self.nostraight = 3.0                            # a run shorter than this does not turn (the metric drops it)
        self.grid_out, self.grid_in = 80, 140


# ====================================================================== the other objects

class _Geo:
    """Ring radii and a 1 px raster of the other drawn objects (ticks, tick labels, tracks, title)."""

    def __init__(self, items: list[Item], ctx: Ctx, par: _Par):
        rl = ctx.radial_layout
        total = int(ctx.total_length)
        feats = rl.features if rl is not None else None
        if feats is not None:
            f_in, f_out = float(feats.all_band_px.inner_px), float(feats.all_band_px.outer_px)
        else:
            f_out = max(i.feat_outer_r for i in items)
            f_in = min(float(i.src.get("feature_inner_radius_px", i.feat_outer_r)) for i in items)
        axis = getattr(rl, "axis", None)
        ax_r = float(axis.radius_px) if axis is not None else float(ctx.radius)
        ax_w = 0.5 * float(axis.stroke_width_px) if axis is not None else 0.5
        self.f_in, self.f_out = f_in, f_out
        self.R = max(f_out + par.edge, ax_r + ax_w)       # outer edge of the ring (features and axis line)
        self.B = min(f_in - par.edge, ax_r - ax_w)        # inner edge of the ring
        self.Bb = self.B                                   # largest radius of an inside label box
        self.padw = par.clear + 0.6                        # box to the ring (analytic test)
        self.padr = par.clear + 1.0                        # box to a raster object
        self.rs_out, self.rs_in = self.R + par.stub, self.B - par.stub
        annuli: list[tuple[float, float]] = []
        marks: list[tuple[float, float, float, float]] = []
        sectors: list[tuple[float, float, float, float]] = []
        font = float(ctx.font_size)
        family = str(items[0].src.get("font_family", "")) if items else ""
        self.tick_font = None
        for slot in (rl.slots if rl is not None else ()):
            if slot.renderer in ("features", "feature_labels"):
                continue
            if slot.renderer == "ticks" and slot.payload is not None:
                tick = slot.payload
                ta, tb = float(tick.tick_band_px.inner_px), float(tick.tick_band_px.outer_px)
                large, _small = circular_tick_intervals(total)
                for pos in range(0, total, large):
                    marks.append((2.0 * math.pi * pos / total, ta, tb, 1.0))
                if tick.label_band_px is not None:
                    la, lb = float(tick.label_band_px.inner_px), float(tick.label_band_px.outer_px)
                    _w, h14 = calculate_bbox_dimensions("0 kbp", family, 14.0, 96)
                    fs = 14.0 * (lb - la) / max(float(h14), 1e-6)
                    self.tick_font = round(fs, 2)
                    for pos in range(large, total, large):
                        text = f"{scaled_tick_text(pos, 1000, large)} kbp" if total < 1_000_000 else f"{pos / 1000000} Mbp"
                        w, _h = calculate_bbox_dimensions(text, family, fs, 96)
                        # the PNG renderer draws the text up to 0.8 band heights inside label_band_px (lower half of the map)
                        sectors.append((2.0 * math.pi * pos / total, la - 0.85 * (lb - la), lb, (0.5 * float(w) + 1.5) / max(la, 1.0)))
                if tb <= self.B + 1e-6 and tb >= self.B - 3.0 * font:
                    self.Bb = min(self.Bb, ta)             # ticks right inside the ring: label boxes stay inside them
                continue
            band = slot.draw_band_px if slot.draw_band_px is not None else slot.reserved_band_px
            if band is None or band.width_px <= 1e-6 or slot.side == "overlay":
                continue
            a, b = float(band.inner_px), float(band.outer_px)
            if a >= f_in - 1e-6 and b <= f_out + 1e-6:
                continue
            annuli.append((a, b))
        if rl is not None and rl.definition_reserved_band_px is not None:
            annuli.append((0.0, float(rl.definition_reserved_band_px.outer_px)))
        self.A_in = max([b for a, b in annuli if b <= self.Bb] + [0.0]) + self.padw
        self.zone_ok = self.Bb - self.A_in > 3.0 * font + self.padw
        far = max([self.R] + [b for _a, b in annuli] + [m[2] for m in marks] + [s[2] for s in sectors])
        E = int(math.ceil(far + 8.0))
        self.E = E
        n = 2 * E + 1
        ax = np.arange(n, dtype=float) - E
        X, Y = np.meshgrid(ax, ax)
        rr = np.hypot(X, Y)
        ex = 0.75
        occ = np.zeros((n, n), dtype=bool)
        for a, b in annuli:
            occ |= (rr >= a - ex) & (rr <= b + ex)
        for th, a, b, hw in marks:
            ux, uy = math.sin(th), -math.cos(th)
            along = X * ux + Y * uy
            occ |= (along >= a - ex) & (along <= b + ex) & (np.abs(X * uy - Y * ux) <= hw + ex)
        if sectors:
            ang = np.arctan2(X, -Y)
            for th, a, b, half in sectors:
                dth = np.abs((ang - th + np.pi) % (2.0 * np.pi) - np.pi)
                occ |= (rr >= a - ex) & (rr <= b + ex) & (dth <= half)
        self.has = {"out": bool(occ[rr > self.R + 0.5].any()), "in": bool(occ[rr < self.B - 0.5].any())}
        ii = np.zeros((n + 1, n + 1), dtype=np.int64)
        ii[1:, 1:] = occ.cumsum(axis=0).cumsum(axis=1)
        self.ii, self.n = ii, n
        self.objects = {"annuli": len(annuli), "tick_marks": len(marks), "tick_labels": len(sectors)}

    def blocked(self, x0, y0, x1, y1):
        """True where the rectangle holds a raster object (arrays broadcast)."""
        E, n = self.E, self.n
        x0, y0, x1, y1 = np.broadcast_arrays(x0, y0, x1, y1)
        ix0 = np.clip(np.ceil(x0 + E), 0, n).astype(np.int64)
        iy0 = np.clip(np.ceil(y0 + E), 0, n).astype(np.int64)
        ix1 = np.clip(np.floor(x1 + E) + 1, 0, n).astype(np.int64)
        iy1 = np.clip(np.floor(y1 + E) + 1, 0, n).astype(np.int64)
        ii = self.ii
        count = ii[iy1, ix1] - ii[iy0, ix1] - ii[iy1, ix0] + ii[iy0, ix0]
        return (count > 0) & (ix1 > ix0) & (iy1 > iy0)

    def seg_blocked(self, p, q):
        """p (2,), q (k,2): True where the leader p-q comes within 2 px of a raster object."""
        q = np.asarray(q, dtype=float).reshape(-1, 2)
        length = np.hypot(q[:, 0] - p[0], q[:, 1] - p[1])
        ns = int(max(2, math.ceil(float(length.max()) / 2.0) + 1))
        t = np.linspace(0.0, 1.0, ns)[None, :]
        x = p[0] + (q[:, 0:1] - p[0]) * t
        y = p[1] + (q[:, 1:2] - p[1]) * t
        return self.blocked(x - 2.0, y - 2.0, x + 2.0, y + 2.0).any(axis=1)


# ====================================================================== one label in one stack

class _Lab:
    """A label in a frame: shore 'out'/'in', q = sx * x is the axis the text runs along (+q = away from the leader)."""

    def __init__(self, it: Item, shore: str, sx: int, geo: _Geo, par: _Par):
        self.it, self.shore, self.sx = it, shore, sx
        ux, uy = math.sin(it.theta), -math.cos(it.theta)
        rs = geo.rs_out if shore == "out" else geo.rs_in
        self.S = (rs * ux, rs * uy)
        self.Sq, self.Sy = sx * self.S[0], self.S[1]
        sign = 1.0 if shore == "out" else -1.0
        self.nq, self.ny = sx * sign * ux, sign * uy       # direction the leader leaves the ring along
        self.upper = uy < 0
        self.d = (-1 if self.upper else 1) * (1 if shore == "out" else -1)    # where a crowded stack grows
        self.w, self.h = it.w, it.h
        self.ring = float(it.src.get("feature_center_radius_px", 0.0)) or math.hypot(*it.anchor)
        self.t = self.base = None
        self.pa = self.pb = None                           # heights (centre of the box) that cost at most plat more than t
        self.far_only = False                              # no spot near the stub (the label lies across a pole)
        self.memo: dict = {}
        self.lo, self.hi = -math.inf, math.inf
        # the way out that is kept free until the label is seated: straight along the stack direction, inside the cone
        lim = math.radians(par.cone - 1.5)
        an, av = math.atan2(self.ny, self.nq), math.atan2(self.d, 0.0)
        turn = (av - an + math.pi) % (2.0 * math.pi) - math.pi
        a = an + max(-lim, min(lim, turn))
        hq, hy = math.cos(a), math.sin(a)
        self.exit, self.slanted = (hq, hy), abs(turn) > lim
        self.edges = [(math.cos(an + s * lim), math.sin(an + s * lim)) for s in (-1.0, 1.0)]   # the two edges of the cone
        reach = 6000.0 if shore == "out" else max(abs(self.Sy) / max(abs(hy), 1e-6) - 1.0, 1.0)
        self.ray = (self.S, (self.S[0] + sx * hq * reach, self.S[1] + hy * reach))
        self._profile(geo, par)

    def cands(self, yc, geo: _Geo, par: _Par):
        """Candidates, shape (len(yc), n): near end q0, attach height ya, leader-to-text gap; and the nearest allowed q0."""
        yc = yc[:, None]
        y0, y1 = yc - 0.5 * self.h, yc + 0.5 * self.h
        near = np.clip(self.Sy, y0, y1)                    # the leader meets the end of the text at the nearest height
        far = np.where(self.Sy < yc, y1, y0)               # ... or at the far corner (a label on its exit ray)
        ady = np.abs(near - self.Sy)
        qp = self.Sq + par.lgap
        if self.shore == "out":
            dyn = np.maximum(np.abs(yc) - 0.5 * self.h - geo.padw, 0.0)
            hug = np.where(dyn < geo.R, np.sqrt(np.maximum(geo.R ** 2 - dyn ** 2, 0.0)) + geo.padw + 1e-6, -1e9)
        else:
            dyf = np.abs(yc) + 0.5 * self.h + geo.padw
            hug = np.where(dyf < geo.Bb, geo.padw + 1e-6 - np.sqrt(np.maximum(geo.Bb ** 2 - dyf ** 2, 0.0)), 1e9)
        base = np.maximum(qp, hug)
        zero = 0.0 * yc
        parts = [base + OFFS[None, :], qp + ady * LEAN[None, :], qp - ady * LEAN[None, :], qp - OFFS[None, 1:12] + zero]
        if self.shore == "in":
            parts.append(-self.w - 0.5 * par.gap - 1e-6 - OFFS[None, :14] + zero)
        heights = [near + 0.0 * part for part in parts]
        gaps = [par.lgap + 0.0 * part for part in parts]
        # the leader exactly on an edge of the cone; on a slanted exit ray also to the far corner of the text end, the
        # text standing back so far that a neighbour's parallel leader clears its corner
        slant = par.lpad * (1.0 + abs(self.exit[0]) / max(abs(self.exit[1]), 1e-6)) + 0.1
        for ya, (uq, uy), gap in [(near, e, par.lgap) for e in self.edges] + ([(far, self.exit, max(par.lgap, slant))] if self.slanted else []):
            run = (ya - self.Sy) / uy if abs(uy) > 1e-6 else -np.ones_like(ya)
            parts.append(np.where(run > 0, self.Sq + run * uq + gap, base))
            heights.append(np.where(run > 0, ya, near))
            gaps.append(np.where(run > 0, gap, par.lgap))
        if par.ray:                                        # the leader exactly on the ray of the feature: one straight line
            lo, hi = 0.0 * yc, np.full_like(yc, np.inf)
            nq, ny = self.nq, self.ny
            if abs(nq) > 1e-6:
                bound = (base - self.Sq - par.lgap) / nq
                lo, hi = (np.maximum(lo, bound), hi) if nq > 0 else (lo, np.minimum(hi, bound))
            else:
                hi = np.where(self.Sq + par.lgap >= base - 1e-9, hi, -1.0)
            if abs(ny) > 1e-6:
                r0, r1 = (y0 - self.Sy) / ny, (y1 - self.Sy) / ny
                lo, hi = np.maximum(lo, np.minimum(r0, r1)), np.minimum(hi, np.maximum(r0, r1))
            else:
                hi = np.where((y0 <= self.Sy) & (self.Sy <= y1), hi, -1.0)
            for off in (0.0, 2.0, 5.0, 10.0, 20.0):
                run = lo + off
                okr = (run <= hi + 1e-9) & (lo <= hi + 1e-9)
                parts.append(np.where(okr, self.Sq + par.lgap + run * nq, base))
                heights.append(np.where(okr, np.clip(self.Sy + run * ny, y0, y1), near))
                gaps.append(np.full_like(base, par.lgap))
        return np.concatenate(parts, axis=1), np.concatenate(heights, axis=1), np.concatenate(gaps, axis=1), base

    def static(self, yc, q0, ya, ga, geo: _Geo, par: _Par, cone_cos: float):
        """Tests that need no other label. yc (n,1), q0, ya and ga (n,m) -> ok, leader cost (run + tw x turn), run length."""
        h, w = self.h, self.w
        y0, y1 = yc - 0.5 * h, yc + 0.5 * h
        aq = q0 - ga
        dq, dy = aq - self.Sq, ya - self.Sy
        ln = np.hypot(dq, dy)
        cosv = np.clip((dq * self.nq + dy * self.ny) / np.maximum(ln, 1e-9), -1.0, 1.0)
        cost = ln + par.tw * np.where(ln < par.nostraight, 0.0, np.degrees(np.arccos(cosv)))
        ok = (dq * self.nq + dy * self.ny >= cone_cos * ln - 1e-9) | (ln < par.short)
        ok &= ~((self.Sq > aq + 1e-9) & (self.Sy > y0) & (self.Sy < y1))        # never through the label's own text
        p = geo.padw
        nx, nyv = np.clip(0.0, q0 - p, q0 + w + p), np.clip(0.0, y0 - p, y1 + p)
        if self.shore == "out":
            ok &= np.hypot(nx, nyv) >= geo.R
        else:
            fx = np.maximum(np.abs(q0 - p), np.abs(q0 + w + p))
            fy = np.maximum(np.abs(y0 - p), np.abs(y1 + p))
            ok &= (np.hypot(fx, fy) <= geo.Bb) & (q0 + w <= -0.5 * par.gap) & (np.hypot(nx, nyv) >= geo.A_in - p)
        # LP-D09, the metric's own test: the far end of the box is not nearer the ring than the end the leader meets
        ring = self.ring
        d0, d1 = np.round(np.hypot(q0 - aq, y0 - ya), 1), np.round(np.hypot(q0 - aq, y1 - ya), 1)
        r0, r1 = np.abs(np.hypot(q0, y0) - ring), np.abs(np.hypot(q0, y1) - ring)
        near = np.where((d0 < d1) | ((d0 == d1) & (r0 <= r1)), r0, r1)
        far = np.minimum(np.abs(np.hypot(q0 + w, y0) - ring), np.abs(np.hypot(q0 + w, y1) - ring))
        ok &= ~(far < near - 0.7)
        if geo.has[self.shore]:
            x0, x1 = (q0, q0 + w) if self.sx > 0 else (-(q0 + w), -q0)
            r = geo.padr
            ok &= ~geo.blocked(x0 - r, y0 - r + 0.0 * q0, x1 + r, y1 + r + 0.0 * q0)
        return ok, cost, ln

    def _profile(self, geo: _Geo, par: _Par):
        """Alone: the y with the shortest leader (target) and how far the label may drift from it."""
        G = par.grid_out if self.shore == "out" else par.grid_in
        ys = self.Sy + np.arange(-G, G + 1.0)
        Q, YA, GA, base = self.cands(ys, geo, par)
        ok, cost, _ln = self.static(ys[:, None], Q, YA, GA, geo, par, par.cone_cos)
        full = ok.any(axis=1)
        best = np.where(ok & (Q <= base + par.fl), cost, np.inf).min(axis=1)
        near = np.isfinite(best)
        if not near.any():                                 # only far from the stub (a label on the other side of a pole)
            best = np.where(ok, cost, np.inf).min(axis=1)
            near = np.isfinite(best)
            self.far_only = True
        if not near.any():
            return
        i = int(np.argmin(best + 1e-3 * np.abs(ys - self.Sy)))
        self.t, self.base = float(ys[i]), float(best[i])
        pa = pb = i
        while pa > 0 and best[pa - 1] <= best[i] + par.plat:
            pa -= 1
        while pb < len(ys) - 1 and best[pb + 1] <= best[i] + par.plat:
            pb += 1
        self.pa, self.pb = float(ys[pa]), float(ys[pb])
        a = b = i
        grow_up, grow_down = (full, near) if self.d < 0 else (near, full)
        while a > 0 and grow_up[a - 1]:
            a -= 1
        while b < len(ys) - 1 and grow_down[b + 1]:
            b += 1
        self.lo = -math.inf if (a == 0 or (self.shore == "out" and self.d < 0)) else float(ys[a])
        self.hi = math.inf if (b == len(ys) - 1 or (self.shore == "out" and self.d > 0)) else float(ys[b])


class _Placed:
    def __init__(self):
        self._b: list = []
        self._s0: list = []
        self._s1: list = []
        self._arr = None

    def add(self, box, s, a):
        self._b.append(box)
        self._s0.append(s)
        self._s1.append(a)
        self._arr = None

    def arrays(self):
        if self._arr is None:
            self._arr = (np.array(self._b, dtype=float).reshape(-1, 4), np.array(self._s0, dtype=float).reshape(-1, 2),
                         np.array(self._s1, dtype=float).reshape(-1, 2))
        return self._arr


class _Pending:
    """Exit rays of the labels of a shore that are not seated yet."""

    def __init__(self, labs: list[_Lab]):
        self.index = {id(lb): k for k, lb in enumerate(labs)}
        self.p = np.array([lb.ray[0] for lb in labs], dtype=float).reshape(-1, 2)
        self.q = np.array([lb.ray[1] for lb in labs], dtype=float).reshape(-1, 2)
        self.live = np.ones(len(labs), dtype=bool)

    def others(self, L: _Lab):
        m = self.live.copy()
        m[self.index[id(L)]] = False
        return self.p[m], self.q[m]

    def done(self, L: _Lab):
        self.live[self.index[id(L)]] = False


# ====================================================================== the layout

class _Eng:
    def __init__(self, ctx: Ctx, geo: _Geo, par: _Par):
        self.ctx, self.geo, self.par = ctx, geo, par
        self.cache: dict = {}
        self.stats = {"pushed": 0, "push_max": 0.0, "forced": 0}

    def lab(self, it: Item, shore: str, sx: int) -> _Lab:
        key = (id(it), shore, sx)
        if key not in self.cache:
            self.cache[key] = _Lab(it, shore, sx, self.geo, self.par)
        return self.cache[key]

    # ---- one stack column: wanted y of every label
    def pack(self, col: list[_Lab]):
        n = len(col)
        sizes = [lb.h + self.par.gap for lb in col]
        # a label with no spot within reach (far on the wrong side of a pole) aims far along its stack
        tg = [lb.t if lb.t is not None else lb.Sy + lb.d * self.par.grid_out for lb in col]
        ys = pack_ordered(tg, sizes)
        if self.par.plat > 0:                              # a label costs nothing within its plateau: pack towards the packed heights
            for _ in range(6):
                tg = [y if lb.pa is None else min(max(y, lb.pa), lb.pb) for lb, y in zip(col, ys)]
                new = pack_ordered(tg, sizes)
                if max(abs(u - v) for u, v in zip(new, ys)) < 0.05:
                    break
                ys = new
        sep = [0.5 * (sizes[k] + sizes[k + 1]) for k in range(n - 1)]
        for k in range(n):
            ys[k] = max(ys[k], col[k].lo, (ys[k - 1] + sep[k - 1]) if k else -1e18)
        for k in range(n - 1, -1, -1):
            ys[k] = min(ys[k], col[k].hi, (ys[k + 1] - sep[k]) if k < n - 1 else 1e18)
        bad = [k for k in range(n) if ys[k] < col[k].lo - 1e-6 or (k and ys[k] - ys[k - 1] < sep[k - 1] - 1e-6)]
        return ys, sep, bad

    # ---- seat one label at y (or the first y further along its stack direction where it fits)
    def try_at(self, L: _Lab, yy: float, pl: _Placed, pend: _Pending):
        """The shortest leader at height yy that fits: ((box, attach point, yy, run), shortest run before the other
        labels are looked at). (None, None) when the label has no spot at this height at all."""
        geo, par = self.geo, self.par
        y0, y1 = yy - 0.5 * L.h, yy + 0.5 * L.h
        key = round(yy, 4)
        if key not in L.memo:                              # the spots of this label at this height (no other label yet)
            ys = np.array([yy])
            Q, YA, GA, _base = L.cands(ys, geo, par)
            ok, ln, _run = L.static(ys[:, None], Q, YA, GA, geo, par, par.cone_cos)
            idx = np.nonzero(ok[0])[0]
            if len(idx):
                q, ln = Q[0][idx], ln[0][idx]
                x0 = q if L.sx > 0 else -(q + L.w)
                boxes = np.stack([x0, np.full(len(q), y0), x0 + L.w, np.full(len(q), y1)], axis=1)
                att = np.stack([L.sx * (q - GA[0][idx]), YA[0][idx]], axis=1)
                L.memo[key] = (q, ln, boxes, att, float(ln.min()))
            else:
                L.memo[key] = None
        if L.memo[key] is None:
            return None, None
        q, ln, boxes, att, alone = L.memo[key]
        stub = np.broadcast_to(np.array(L.S, dtype=float), att.shape)
        good = np.ones(len(q), dtype=bool)
        lo, hi = min(L.Sy, y0) - par.lpad - 1.0, max(L.Sy, y1) + par.lpad + 1.0     # rows the boxes and leaders touch
        pb, s0, s1 = pl.arrays()
        if len(pb):
            pb = pb[(pb[:, 3] > lo - par.gap) & (pb[:, 1] < hi + par.gap)]
            m = (np.minimum(s0[:, 1], s1[:, 1]) < hi) & (np.maximum(s0[:, 1], s1[:, 1]) > lo)
            s0, s1 = s0[m], s1[m]
        if len(pb):
            g = par.gap - 1e-6
            good &= ~((boxes[:, None, 0] < pb[None, :, 2] + g) & (pb[None, :, 0] < boxes[:, None, 2] + g)
                      & (boxes[:, None, 1] < pb[None, :, 3] + g) & (pb[None, :, 1] < boxes[:, None, 3] + g)).any(axis=1)
            sel = np.nonzero(good)[0]
            if len(sel):
                good[sel[_seg_box_hit(stub[sel], att[sel], pb, par.lpad).any(axis=1)]] = False
        if len(s0) and good.any():
            sel = np.nonzero(good)[0]
            bad = _seg_box_hit(s0, s1, boxes[sel], par.lpad).any(axis=0)
            bad |= _seg_cross_matrix(stub[sel], att[sel], s0, s1).any(axis=1)
            good[sel[bad]] = False
        rp, rq = pend.others(L)
        if len(rp):
            m = (np.minimum(rp[:, 1], rq[:, 1]) < hi) & (np.maximum(rp[:, 1], rq[:, 1]) > lo)
            rp, rq = rp[m], rq[m]
        if len(rp) and good.any():                       # keep the way out of every label still to be seated
            sel = np.nonzero(good)[0]
            bad = _seg_box_hit(rp, rq, boxes[sel], par.lpad).any(axis=0)
            bad |= _seg_cross_matrix(stub[sel], att[sel], rp, rq).any(axis=1)
            good[sel[bad]] = False
        if geo.has[L.shore] and good.any():
            sel = np.nonzero(good)[0]
            good[sel[geo.seg_blocked(L.S, att[sel])]] = False
        if not good.any():
            return None, alone
        sel = np.nonzero(good)[0]
        j = sel[int(np.argmin(ln[sel] + 1e-6 * np.abs(q[sel] - (L.Sq + par.lgap))))]
        return (tuple(float(v) for v in boxes[j]), (float(att[j][0]), float(att[j][1])), yy, float(ln[j])), alone

    def seat(self, L: _Lab, y: float, pl: _Placed, pend: _Pending, limit: float | None, tight: float | None):
        """Seat the label at y or further along its stack: the smallest leader run + push.
        `tight` is the push that puts it exactly one row beyond the label seated before it in the stack."""
        best = None                                       # (leader run + push, push, result)
        floor = None                                      # shortest run at the last height tried, other labels aside
        pushes = PUSH
        if tight is not None and tight > 0:
            pushes = [p for p in NEAR if p < tight] + [tight] + [tight + p for p in PUSH[1:]]
        for p in pushes:
            yy = y + L.d * p
            if limit is not None and (yy - limit) * L.d > 0:
                break
            if best is not None:                           # a further push cannot pay any more
                if self.par.tw == 0 and not self.par.ray:
                    bound = floor                          # design I: the shortest run at the last height tried
                elif (yy - L.Sy) * L.d >= 0:
                    bound = max(0.0, abs(yy - L.Sy) - 0.5 * L.h)   # a leader is at least as long as the vertical gap
                else:
                    bound = 0.0
                if bound + p >= best[0]:
                    break
            res, alone = self.try_at(L, yy, pl, pend)
            if alone is not None:
                floor = alone
            if res is not None and (best is None or res[3] + p < best[0] - 1e-9):
                best = (res[3] + p, p, res[:3])
        if best is None:
            return None
        if best[1] > 0:
            self.stats["pushed"] += 1
            self.stats["push_max"] = max(self.stats["push_max"], best[1])
        return best[2]

    def run_column(self, col: list[_Lab], ys, sep, pl: _Placed, pend: _Pending, out: list, tag: str, bounded: bool):
        """Seat a column given top to bottom. Each half is seated in the order its stack grows."""
        up = [k for k, lb in enumerate(col) if lb.upper]
        low = [k for k, lb in enumerate(col) if not lb.upper]
        for half in (up, low):
            if not half:
                continue
            order = half if col[half[0]].d > 0 else half[::-1]
            prev = None
            for k in order:
                L = col[k]
                # start at the label's own best height (side by side with its neighbour where the ring is flat);
                # the packed height counts only where it lies towards the equator (room for the labels to come)
                y = ys[k] if L.t is None or (ys[k] - L.t) * L.d <= 0 else L.t
                tight = None
                if prev is not None:
                    tight = (prev[1] + L.d * sep[min(k, prev[0])] - y) * L.d
                limit = None
                if bounded:
                    limit = (L.hi if L.d > 0 else L.lo) + L.d * 0.5
                    if not math.isfinite(limit):
                        limit = y + L.d * 400.0
                res = self.seat(L, y, pl, pend, limit, tight)
                if res is None:
                    if bounded:
                        return L.it
                    res = self.force(L, y, pl)
                box, att, yy = res
                pl.add(box, L.S, att)
                pend.done(L)
                out.append((L.it, box, [L.it.anchor, L.S, att], "start" if L.sx > 0 else "end", tag))
                prev = (k, yy)
        return None

    def force(self, L: _Lab, y: float, pl: _Placed):
        """Last resort (never expected): far out along the stack direction, no test against other labels."""
        self.stats["forced"] += 1
        pb, _s0, _s1 = pl.arrays()
        if not len(pb):
            pb = np.array([[0.0, y, 0.0, y]])
        yy = (min(pb[:, 1].min(), y) - L.h - 10.0) if L.d < 0 else (max(pb[:, 3].max(), y) + L.h + 10.0)
        q = L.Sq + self.par.lgap
        x0 = q if L.sx > 0 else -(q + L.w)
        return (x0, yy - 0.5 * L.h, x0 + L.w, yy + 0.5 * L.h), (L.sx * (q - self.par.lgap), yy - L.d * 0.5 * L.h), yy

    # ---- outside shore
    def split(self, items: list[Item], cuts, far: bool = False):
        """The two columns (each top to bottom) for the cuts at 12 and 6 o'clock; cuts=None searches them.
        far=True lets a column take labels from across the pole that have no spot near their feature."""
        keys = [i.theta - 2.0 * math.pi if i.theta > math.radians(300.0) else i.theta for i in items]
        order = sorted(range(len(items)), key=lambda k: keys[k])
        seq = [items[k] for k in order]
        ks = [keys[k] for k in order]
        n = len(seq)

        def columns(c1, c2):
            return ([self.lab(it, "out", 1) for it in seq[c1:c2]],
                    [self.lab(it, "out", -1) for it in (seq[c2:] + seq[:c1])[::-1]])

        if cuts is None:
            def candidates(centre: float) -> tuple[list[int], int]:
                cs = [c for c in range(n + 1) if (c < n and abs(ks[c] - centre) <= self.par.cutmax)
                      or (c > 0 and abs(ks[c - 1] - centre) <= self.par.cutmax)]
                mid = sum(1 for k in ks if k < centre)
                cs = sorted(sorted(set(cs) | {mid}, key=lambda c: abs(c - mid))[:int(self.par.cutn)])
                return cs, mid

            def cost(col):
                if not col:
                    return 0.0
                if any(lb.t is None or (lb.far_only and not far) for lb in col):
                    return 1e9
                ys, _sep, _bad = self.pack(col)
                return sum(lb.base + abs(y - lb.t) + 0.002 * (y - lb.t) ** 2 for lb, y in zip(col, ys))

            def value(c1, c2):
                if c2 < c1:
                    return math.inf
                right, left = columns(c1, c2)
                return cost(right) + cost(left)

            (c1s, c1), (c2s, c2) = candidates(0.0), candidates(math.pi)
            for _ in range(2):                             # one pole at a time, twice
                c1 = min(c1s, key=lambda c: (value(c, c2), abs(c - c1)))
                c2 = min(c2s, key=lambda c: (value(c1, c), abs(c - c2)))

            def cut(c: int) -> float:
                if 0 < c < n:
                    return 0.5 * (ks[c - 1] + ks[c])
                return ks[0] - 0.01 if c == 0 else ks[-1] + 0.01

            cuts = (cut(c1), cut(c2))
        c1 = sum(1 for k in ks if k < cuts[0])
        c2 = max(c1, sum(1 for k in ks if k < cuts[1]))
        return columns(c1, c2), cuts

    def outside(self, items: list[Item], cuts):
        """Records and the seating statistics for given cuts."""
        self.stats = {"pushed": 0, "push_max": 0.0, "forced": 0}
        if not items:
            return [], dict(self.stats)
        (right, left), cuts = self.split(items, cuts)
        pl = _Placed()
        pend = _Pending(right + left)
        recs: list = []
        for col, tag in ((right, "R"), (left, "L")):
            if not col:
                continue
            ys, sep, _bad = self.pack(col)
            self.run_column(col, ys, sep, pl, pend, recs, tag, bounded=False)
        return recs, dict(self.stats)

    # ---- inside shore
    def inside(self, sel: list[Item]):
        """Records, or (None, the item that does not fit)."""
        if not sel:
            return [], None
        right = sorted((i for i in sel if math.sin(i.theta) > 0), key=lambda i: i.theta)
        left = sorted((i for i in sel if math.sin(i.theta) <= 0), key=lambda i: -i.theta)
        pl = _Placed()
        cols = [[self.lab(it, "in", sx) for it in items] for items, sx in ((right, -1), (left, 1))]
        pend = _Pending(cols[0] + cols[1])
        recs: list = []
        for col in cols:
            if not col:
                continue
            for lb in col:
                if lb.t is None:
                    return None, lb.it
            ys, sep, bad = self.pack(col)
            if bad:
                return None, max((col[k] for k in bad), key=lambda lb: lb.base).it
            failed = self.run_column(col, ys, sep, pl, pend, recs, "in", bounded=True)
            if failed is not None:
                return None, failed
        return recs, None


# ====================================================================== checks (the same tests as scripts/lp_metrics.py)

def _check(recs) -> dict:
    boxes = [r[1] for r in recs]
    segs = [list(zip(r[2], r[2][1:])) for r in recs]
    n = len(recs)

    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def cross(p, q, r, s):
        d1, d2, d3, d4 = orient(r, s, p), orient(r, s, q), orient(p, q, r), orient(p, q, s)
        return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-6 for d in (d1, d2, d3, d4))

    B = np.array(boxes, dtype=float).reshape(-1, 4)
    span = [(min(p[0] for p in r[2]), min(p[1] for p in r[2]), max(p[0] for p in r[2]), max(p[1] for p in r[2])) for r in recs]
    overlap = cross_n = through = 0
    gap = math.inf
    for i in range(n):
        a = boxes[i]
        p = np.array([s[0] for s in segs[i]], dtype=float)
        q = np.array([s[1] for s in segs[i]], dtype=float)
        hit = _seg_box_hit(p, q, B, -0.5).any(axis=0)
        hit[i] = False
        through += int(hit.sum())
        for j in range(i + 1, n):
            b = boxes[j]
            dx, dy = max(a[0] - b[2], b[0] - a[2]), max(a[1] - b[3], b[1] - a[3])
            gap = min(gap, math.hypot(dx, dy) if dx > 0 and dy > 0 else max(dx, dy))
            if min(a[2], b[2]) - max(a[0], b[0]) > 0.5 and min(a[3], b[3]) - max(a[1], b[1]) > 0.5:
                overlap += 1
            if span[i][0] > span[j][2] or span[j][0] > span[i][2] or span[i][1] > span[j][3] or span[j][1] > span[i][3]:
                continue
            if any(cross(p1, q1, r1, s1) for p1, q1 in segs[i] for r1, s1 in segs[j]):
                cross_n += 1
    sharp = tight = 0
    for r in recs:
        pts = r[2]
        pts = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= 3.0]
        worst = 180.0
        for a, b, c in zip(pts, pts[1:], pts[2:]):
            u, v = (a[0] - b[0], a[1] - b[1]), (c[0] - b[0], c[1] - b[1])
            nu, nv = math.hypot(*u), math.hypot(*v)
            if nu >= 3.0 and nv >= 3.0:
                worst = min(worst, math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv))))))
        sharp += worst < 89.0
        tight += worst < 110.0
    return {"overlap": overlap, "cross": cross_n, "through": through, "gap": round(gap, 2) if n > 1 else None,
            "acute": int(sharp), "lt110": int(tight)}


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
    if inner_req:                                   # --labels both: these labels were asked for inside the ring
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
    def lead(pts) -> float:                             # what a leader costs: its length plus tw px per degree of turn
        return _plen(pts) + par.tw * _turn(pts)

    def total(recs) -> float:
        return sum(lead(r[2]) for r in recs)

    def outside_best(sel: list[Item]):
        """The outside layout for the better of two cuts: without and with labels taken from across a pole."""
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
        # an outside label moves inside when its leader gets shorter there; ten candidates, then the outside again
        sel_in: list[Item] = []
        cur_in: list = []
        rest, cur_out, cur_stats = list(out_items), out_recs, stats
        rejected: set = set()
        while True:
            length = {id(r[0]): lead(r[2]) for r in cur_out}
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
            again = outside_best(rest)                       # the cuts for the labels that stay outside
            if total(again[0]) < total(cur_out) - 1e-6:
                cur_out, cuts_in, cur_stats = again
            else:
                cuts_in = cuts
            if total(cur_out) + total(cur_in) < total(out_recs) - 1e-6:
                out_recs, in_recs, cuts, stats, moves = cur_out, cur_in, cuts_in, cur_stats, len(sel_in)
    recs = out_recs + in_recs
    hide = bool(int(ctx.opt("hide", 0)))
    result = [emit(it, box, None if hide else pts, align=align, tag=tag) for it, box, pts, align, tag in recs]
    if hide:
        for label in result:
            label["text_fill"] = "none"
    ctx.notes.update(
        inside=len(in_recs), outside=len(out_recs), asked_inside=len(inner_req), overflow=[i.text for i in overflow], moves=moves,
        total=round(sum(_plen(r[2]) for r in recs), 1), turn_sum=round(sum(_turn(r[2]) for r in recs), 1), R=round(geo.R, 2), B=round(geo.B, 2), Bb=round(geo.Bb, 2),
        A_in=round(geo.A_in, 2), zone_ok=geo.zone_ok, objects=geo.objects, raster=geo.has, tick_font=geo.tick_font,
        cuts=[round(math.degrees(c), 1) for c in cuts] if cuts else None, stats=stats, check=_check(recs),
    )
    return result
