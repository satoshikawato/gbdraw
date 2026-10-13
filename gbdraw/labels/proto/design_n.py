"""Design N (2026-10-13): one circle around the ring, as near the ring as the labels allow.

Built from the audit of the Owner's words LP-D06 to LP-D39 (task folder, logs/AUDIT-N-20261013.md).  Where each of
his rules lives in this design:

- Labels follow the circumference; a jagged stack is ugly (LP-D26, LP-D27); the labels stay on one arc that follows
  the axis, with no inflection point (LP-D38, LP-D39).  Every near end (the contact point of the leader on the label)
  lies on ONE circle concentric with the ring.  Nothing bends the outline: no label leaves the circle outward along
  its ray, and the circle is not raised locally where the labels are dense.
- The arc is as near the ring as the labels allow (LP-D10: no long horizontal runs to an ellipse; LP-D22: the stack
  stays next to the ring; the Q look he approved, LP-D28): the radius is the least at which every label has a place.
- A limit at the top and bottom, today's, set by today's parameter (LP-D34, LP-D35): no label box passes
  `arc_y_radius` + text height, so the circle never grows beyond that.  The end labels stand at 12 and 6 o'clock
  (LP-D37) because the circle reaches the axis.
- Leaders never cross (LP-D07, LP-D36): the near ends keep the clockwise order of the features on the circle and no
  piece runs back into the ring, so no two last pieces can cross; counted after the fact all the same.
- Contact points (LP-D23, LP-D25): the box is the ink box; the leader stops `lgap` before one of the 12 points; a
  near-horizontal last piece ends at the middle of the short side (LP-D18), a steep one at the facing vertex or on
  the facing long side, never on the far half (LP-D15).
- Rows 0.22 em of ink gap, not below 4 px (LP-D30, LP-D31); the clearance to the ring and to every other object is
  the same distance (LP-D06, LP-D10: ticks, tick labels, tracks from the raster of design J).
- Straight in the middle of a cluster, bends growing outward, symmetric fans (LP-D19, LP-D12, LP-D22): the places
  on the circle minimise the total displacement of the labels from their own rays (L1), so the middle label of a
  compressed run stands on its ray and the others spread evenly to both sides.
- Every leader crosses the axis at a right angle (LP-D40): a radial stub from the feature to the elbow just outside the
  ring, then one piece to the contact point.  The stub is dropped (one line) only when that line leaves the feature
  within `scone` degrees of its ray (2 degrees: not visible) and crosses nothing.
- The turn at the elbow stays within `cone` except within `polezone` degrees of a pole, where the circle runs
  across the rays and a label beside its neighbour needs more (up to `conemax`); the cost keeps it as small as the
  room allows.  This is the fan of today's figures at the poles (LP-D33: today's poles are orderly).

What gives way when the labels do not fit on the largest circle under the limit (the Owner has not said; OD-N1,
shown as pictures): first the turn away from the poles opens from `cone` to `conemax`; then (`give=limit`) the limit
moves out in steps of `lstep` with the rows at the rule, or (`give=rows,limit`) the rows shrink first, not below
`sqmin`, and the limit moves only after that.

Options (GBDRAW_LABEL_PROTO_OPTS):
  ylim     limit: px from the map centre that no label box passes at the top and bottom   [today's: arc_y_radius + text height]
  rowem    ink gap of stacked rows as a fraction of the font size                         [0.22]
  rowmin   floor of that gap, px                                                          [4]
  give     limit | rows,limit  (see above)                                                [limit]
  sqmin    smallest ink gap with give=rows,limit                                          [2.5]
  lgap     the leader stops this far before the contact point                            [2.5]
  lpad     a leader passes another label's box no nearer than this                       [2.4]
  cone     turn at the elbow that costs nothing, degrees                                  [40]
  conemax  largest turn                                                                   [88]
  polezone degrees from a pole within which a turn up to conemax is allowed from the start [35]
  scone    largest angle between a one-piece leader and the feature's own ray            [2]
  wc       weight of the squared degrees beyond cone                                     [2]
  vpen     px of cost for a contact point that is not the natural one for the direction  [4]
  step     grid of the places on the circle, degrees                                     [0.15]
  lstep    px the limit moves out per step when the labels do not fit                    [8]
  objects  1: labels and leaders keep clear of ticks, tick labels and tracks             [1]
  axis     nearest a label box may stand to the vertical axis, px                        [6]
"""

from __future__ import annotations

import math
import time

import numpy as np

from .common import Ctx, Item, segment_min_radius, segments_cross
from .design_g import _seg_box_hit, _seg_cross_matrix
from .design_j import _Geo, _Par, _emit_ink, _ink_item

INF = float("inf")
TWO_PI = 2.0 * math.pi
# contact points, in the frame of a half whose labels run away from the axis toward +x: (box left end - point, box top - point)
# in units of the width and the height.  L = middle of the short side facing the axis; BL / TL = the facing vertices;
# B2 / T2 = the middle of the long sides.
PORT_L, PORT_BL, PORT_B2, PORT_TL, PORT_T2 = range(5)
PORT_OFF = np.array([(0.0, -0.5), (0.0, -1.0), (-0.5, -1.0), (0.0, 0.0), (-0.5, 0.0)])


class _P:
    def __init__(self, ctx: Ctx, items: list[Item], ink: list[Item], geo: _Geo):
        o = ctx.opt
        self.rowem, self.rowmin = o("rowem", 0.22), o("rowmin", 4.0)
        self.gap0 = max(float(ctx.gap), self.rowmin, self.rowem * float(ctx.font_size) - 0.1) + 0.05   # between boxes; the ink gap comes out 0.1 above
        self.sqmin = o("sqmin", 2.5) + 0.05
        self.lgap, self.lpad = o("lgap", 2.5), o("lpad", 2.4)
        self.cone, self.conemax, self.scone = o("cone", 40.0), o("conemax", 88.0), o("scone", 2.0)
        self.polezone = math.radians(o("polezone", 35.0))
        self.wc, self.vpen = o("wc", 2.0), o("vpen", 4.0)
        self.step = math.radians(o("step", 0.15))
        self.lstep = o("lstep", 8.0)
        self.objects, self.axis = int(o("objects", 1)), o("axis", 6.0)
        self.give = str(ctx.options.get("give", "limit"))
        self.re = geo.rs_out                                     # the elbow: `stub` px beyond the outer edge of the ring (features and axis line), as in J and Q
        self.rk = geo.R + geo.padw                               # no box nearer the centre than this (the row distance from the ring edge)
        self.hmax = max(i.h for i in ink)                        # tallest ink
        self.ylim0 = o("ylim", ctx.arc_y_radius + max(i.h for i in items))   # today's extent: the label ellipse plus one text box
        self.rmin = self.rk + 0.5


class _Half:
    """The outside labels of one half, in a frame where the half is the right one (x >= 0, clockwise from 12 o'clock)."""

    def __init__(self, items: list[Item], sign: int, p: _P):
        self.sign, self.items = sign, items
        self.th = np.array([it.theta if sign > 0 else (TWO_PI - it.theta) % TWO_PI for it in items])
        self.w, self.h = np.array([it.w for it in items]), np.array([it.h for it in items])
        self.ex, self.ey = p.re * np.sin(self.th), -p.re * np.cos(self.th)
        self.polar = np.minimum(self.th, math.pi - self.th) <= p.polezone   # near 12 or 6 o'clock


class _Tab:
    """The candidate places of one label (arrays over the candidates)."""

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)

    @property
    def boxes(self):
        return np.column_stack([self.x0, self.y0, self.x1, self.y1])

    def seg(self):
        n = len(self.qx)
        return np.column_stack([np.full(n, self.ex), np.full(n, self.ey)]), np.column_stack([self.qx, self.qy])


def _table(half: _Half, k: int, rho: float, T: float, ylim: float, p: _P, geo, forbid: set) -> _Tab:
    th, w, h, ex, ey = half.th[k], half.w[k], half.h[k], half.ex[k], half.ey[k]
    lam = min(math.radians(70.0), math.acos(min(1.0, p.re / rho)) + math.radians(3.0))   # beyond this the piece runs back into the ring
    m = int(math.ceil(lam / p.step))
    j = np.arange(-m, m + 1)
    psi = th + j * p.step                                   # polar angle of the near end on the circle
    keep = (psi >= 0.0) & (psi <= math.pi)
    j, psi = j[keep], psi[keep]
    px, py = rho * np.sin(psi), -rho * np.cos(psi)
    cx, cy = px - ex, py - ey                               # the last piece, elbow -> near end
    ln = np.maximum(np.hypot(cx, cy), 1e-9)
    ux, uy = math.sin(th), -math.cos(th)
    along = cx * ux + cy * uy
    turn = np.degrees(np.arctan2(np.abs(cx * uy - cy * ux), along))
    gam = np.degrees(np.arctan2(-cy, cx))                   # direction of the last piece: 0 = away from the axis, 90 = up
    qx, qy = px - p.lgap * cx / ln, py - p.lgap * cy / ln   # where the leader stops
    sel, port, extra = [], [], []
    for mask, pt in ((np.abs(gam) <= 45.0, PORT_L), (gam > 45.0, PORT_BL), (gam > 70.0, PORT_B2), (gam < -45.0, PORT_TL), (gam < -70.0, PORT_T2)):
        idx = np.flatnonzero(mask)
        if not len(idx):
            continue
        if pt == PORT_L:
            pen = np.zeros(len(idx))
        else:                                               # the vertex is natural up to 80 degrees from the horizontal, the middle of the long side beyond
            pen = np.where((np.abs(gam[idx]) > 80.0) == (pt in (PORT_B2, PORT_T2)), 0.0, p.vpen)
        sel.append(idx)
        port.append(np.full(len(idx), pt))
        extra.append(pen)
    sel, port, extra = np.concatenate(sel), np.concatenate(port), np.concatenate(extra)
    x0, y0 = px[sel] + PORT_OFF[port, 0] * w, py[sel] + PORT_OFF[port, 1] * h
    x1, y1 = x0 + w, y0 + h
    ok = (along[sel] > 0.5) & (turn[sel] <= T) & (x0 >= p.axis) & (y0 >= -ylim) & (y1 <= ylim)
    ok &= np.hypot(x0, np.maximum(np.maximum(y0, 0.0), -y1)) >= p.rk
    if geo is not None:
        rx0, rx1 = (x0, x1) if half.sign > 0 else (-x1, -x0)
        ok &= ~geo.blocked(rx0 - geo.padr, y0 - geo.padr, rx1 + geo.padr, y1 + geo.padr)
        t = np.linspace(0.0, 1.0, max(2, int(ln.max() / 2.0) + 1))[None, :]    # the last piece, sampled every 2 px, beyond 8 px from the ring
        sx, sy = ex + t * (qx[sel][:, None] - ex), ey + t * (qy[sel][:, None] - ey)
        rsx = sx if half.sign > 0 else -sx
        hit = geo.blocked(rsx - 0.5, sy - 0.5, rsx + 0.5, sy + 0.5) & (np.hypot(sx, sy) > geo.R + 8.0)
        ok &= ~hit.any(axis=1)
    if forbid:
        ok &= ~np.isin(j[sel] * 8 + port, np.fromiter(forbid, dtype=np.int64))
    cost = np.where(ok, rho * np.abs(psi[sel] - th) + p.wc * np.maximum(0.0, turn[sel] - p.cone) ** 2 + extra, INF)
    return _Tab(j=j[sel], port=port, psi=psi[sel], px=px[sel], py=py[sel], qx=qx[sel], qy=qy[sel],
                x0=x0, y0=y0, x1=x1, y1=y1, turn=turn[sel], cost=cost, ex=ex, ey=ey)


def _admissible(A: _Tab, B: _Tab, gap: float, lpad: float) -> np.ndarray:
    """(nA, nB): B after A on the circle, the boxes apart by `gap` sideways or in height, neither leader through the other's box."""
    ok = B.psi[None, :] > A.psi[:, None] + 1e-9
    xsep = np.maximum(B.x0[None, :] - A.x1[:, None], A.x0[:, None] - B.x1[None, :])
    ysep = np.maximum(B.y0[None, :] - A.y1[:, None], A.y0[:, None] - B.y1[None, :])
    ok &= (xsep >= gap) | (ysep >= gap)
    pB, qB = B.seg()
    pA, qA = A.seg()
    ok &= ~_seg_box_hit(pB, qB, A.boxes, lpad).T
    ok &= ~_seg_box_hit(pA, qA, B.boxes, lpad)
    return ok


def _solve(half: _Half, rho: float, T: float, gap: float, ylim: float, p: _P, geo):
    """Places of the half's labels on the circle: least total cost in clockwise order (a list), or None when they do not fit."""
    n = len(half.items)
    if n == 0:
        return []
    forbid: list[set] = [set() for _ in range(n)]
    for _ in range(8):                                      # pairs further apart than neighbours are checked after the fact
        tabs = [_table(half, k, rho, p.conemax if half.polar[k] else T, ylim, p, geo, forbid[k]) for k in range(n)]
        best, back = tabs[0].cost, []
        for k in range(1, n):
            m = np.where(_admissible(tabs[k - 1], tabs[k], gap, p.lpad), best[:, None], INF)
            arg = m.argmin(axis=0)
            best = m[arg, np.arange(len(tabs[k].cost))] + tabs[k].cost
            back.append(arg)
        if not np.isfinite(best).any():
            return None
        pick = [int(best.argmin())]
        for arg in reversed(back):
            pick.append(int(arg[pick[-1]]))
        pick.reverse()
        boxes = np.array([[t.x0[i], t.y0[i], t.x1[i], t.y1[i]] for t, i in zip(tabs, pick)])
        E = np.array([[t.ex, t.ey] for t in tabs])
        Q = np.array([[t.qx[i], t.qy[i]] for t, i in zip(tabs, pick)])
        xsep = np.maximum(boxes[None, :, 0] - boxes[:, None, 2], boxes[:, None, 0] - boxes[None, :, 2])
        ysep = np.maximum(boxes[None, :, 1] - boxes[:, None, 3], boxes[:, None, 1] - boxes[None, :, 3])
        bad = ~((xsep >= gap) | (ysep >= gap)) | _seg_box_hit(E, Q, boxes, p.lpad)
        bad |= bad.T
        pairs = np.argwhere(bad & np.triu(np.ones((n, n), dtype=bool), 2))
        if not len(pairs):
            return [dict(k=k, psi=float(t.psi[i]), px=float(t.px[i]), py=float(t.py[i]), qx=float(t.qx[i]), qy=float(t.qy[i]),
                         port=int(t.port[i]), box=tuple(map(float, boxes[k])), turn=float(t.turn[i]), ex=float(t.ex), ey=float(t.ey))
                    for k, (t, i) in enumerate(zip(tabs, pick))]
        for i, jj in pairs:                                 # the later label may not stand where it meets the earlier one
            t = tabs[jj]
            xs = np.maximum(t.x0 - boxes[i, 2], boxes[i, 0] - t.x1)
            ys = np.maximum(t.y0 - boxes[i, 3], boxes[i, 1] - t.y1)
            coll = ~((xs >= gap) | (ys >= gap))
            pj, qj = t.seg()
            coll |= _seg_box_hit(pj, qj, boxes[i:i + 1], p.lpad)[:, 0]
            coll |= _seg_box_hit(E[i:i + 1], Q[i:i + 1], t.boxes, p.lpad)[0]
            forbid[jj].update((t.j[coll] * 8 + t.port[coll]).tolist())
    return None


class _Search:
    def __init__(self, halves, p: _P, geo):
        self.halves, self.p, self.geo, self.probes = halves, p, geo, 0

    def fit(self, rho: float, T: float, gap: float, ylim: float):
        self.probes += 1
        out = []
        for h in self.halves:
            res = _solve(h, rho, T, gap, ylim, self.p, self.geo)
            if res is None:
                return None
            out.append(res)
        return out

    def least_rho(self, T: float, gap: float, ylim: float):
        """Least radius, up to the limit, at which the labels fit, with their places; None when not even at the limit."""
        lo, hi = self.p.rmin, ylim - self.p.hmax
        best = self.fit(hi, T, gap, ylim)
        if best is None:
            return None
        r = self.fit(lo, T, gap, ylim)
        if r is not None:
            return lo, r
        while hi - lo > 1.0:
            mid = 0.5 * (lo + hi)
            r = self.fit(mid, T, gap, ylim)
            if r is None:
                lo = mid
            else:
                hi, best = mid, r
        return hi, best


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    if not items:
        return []
    t0 = time.time()
    ctx.options.setdefault("rowem", "0.22")
    ctx.options.setdefault("rowmin", "4")
    ink = [_ink_item(i, ctx.font_size) for i in items]
    geo = _Geo(items, ctx, _Par(ctx))
    p = _P(ctx, items, ink, geo)
    if not p.objects:
        geo = None
    halves = [_Half(sorted((i for i in ink if math.sin(i.theta) >= 0.0), key=lambda i: i.theta), 1, p),
              _Half(sorted((i for i in ink if math.sin(i.theta) < 0.0), key=lambda i: -i.theta), -1, p)]
    S = _Search(halves, p, geo)

    # what gives way, in order: the turn away from the poles, then the rows (give=rows,limit), then the limit
    gap, ylim, found, T = p.gap0, p.ylim0, None, p.cone
    for T in [*np.arange(p.cone, p.conemax, 8.0), p.conemax]:
        found = S.least_rho(T, gap, ylim)
        if found:
            break
    if found is None and p.give.startswith("rows"):
        g = gap - 0.25
        while found is None and g >= p.sqmin - 1e-9:
            found, gap = S.least_rho(T, g, ylim), g
            g -= 0.25
    while found is None:
        ylim += p.lstep
        if ylim > p.ylim0 + 1200.0:
            raise RuntimeError("design N: the labels do not fit even 1200 px beyond the limit")
        found = S.least_rho(T, gap, ylim)
    rho, results = found

    # real coordinates
    recs = []
    for h, placed in zip(halves, results):
        for rec in placed:
            it, s = h.items[rec["k"]], h.sign
            x0, y0, x1, y1 = rec["box"]
            box = (x0, y0, x1, y1) if s > 0 else (-x1, y0, -x0, y1)
            P, E, Q = (s * rec["px"], rec["py"]), (s * rec["ex"], rec["ey"]), (s * rec["qx"], rec["qy"])
            recs.append(dict(it=it, box=box, P=P, pts=[it.anchor, E, Q], turn=rec["turn"], align="start" if s > 0 else "end"))

    # one straight line where it leaves the feature near its ray and crosses nothing (least turned first)
    boxes = np.array([r["box"] for r in recs])
    for i in sorted(range(len(recs)), key=lambda i: recs[i]["turn"]):
        r, it = recs[i], recs[i]["it"]
        A, P = it.anchor, r["P"]
        dx, dy = P[0] - A[0], P[1] - A[1]
        L = math.hypot(dx, dy)
        if L < 1e-9 or math.degrees(math.acos(max(-1.0, min(1.0, (dx * math.sin(it.theta) - dy * math.cos(it.theta)) / L)))) > p.scone:
            continue
        if segment_min_radius(A, P) < math.hypot(*A) - 0.5:
            continue
        Q = (P[0] - p.lgap * dx / L, P[1] - p.lgap * dy / L)
        if _seg_box_hit(np.array([A]), np.array([Q]), np.delete(boxes, i, axis=0), p.lpad).any():
            continue
        if any(segments_cross(A, Q, u, v) for jx, o in enumerate(recs) if jx != i for u, v in zip(o["pts"], o["pts"][1:])):
            continue
        if geo is not None:
            t = np.linspace(0.0, 1.0, max(2, int(L / 2.0) + 1))
            sx, sy = A[0] + t * (Q[0] - A[0]), A[1] + t * (Q[1] - A[1])
            if (geo.blocked(sx - 0.5, sy - 0.5, sx + 0.5, sy + 0.5) & (np.hypot(sx, sy) > geo.R + 8.0)).any():
                continue
        r["pts"] = [A, Q]

    segs = [(u, v, i) for i, r in enumerate(recs) for u, v in zip(r["pts"], r["pts"][1:])]
    sp, sq = np.array([s[0] for s in segs]), np.array([s[1] for s in segs])
    own = np.array([s[2] for s in segs])
    cross = _seg_cross_matrix(sp, sq, sp, sq) & (own[:, None] != own[None, :])
    ctx.notes.update(rho=round(rho, 1), rho_min=round(p.rmin, 1), turn_max=float(T), turn_free=p.cone,
                     gap_asked=round(p.gap0 - 0.05, 2), gap=round(gap - 0.05, 2), limit_asked=round(p.ylim0, 1), limit=round(ylim, 1),
                     ring=round(p.rk, 1), elbow=round(p.re, 1), probes=S.probes, seconds=round(time.time() - t0, 1),
                     straight=sum(1 for r in recs if len(r["pts"]) == 2),
                     turned_over_70=sum(1 for r in recs if len(r["pts"]) == 3 and r["turn"] > 70.0),
                     crossings=int(cross.sum() // 2), n_right=len(halves[0].items), n_left=len(halves[1].items), give=p.give)
    return [_emit_ink(r["it"], r["box"], r["pts"], r["align"], "N", ctx.font_size) for r in recs]
