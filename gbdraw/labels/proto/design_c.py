"""Design C: free placement by collision-driven search (study prototype).

Each label is a horizontal box placed anywhere (not tied to an arc).  Its leader is anchor -> radial stub just outside
the ring -> one straight segment to the nearest edge point of the box, bounded in length, leaving the stub outward.
1. Every label first gets a reserved fallback slot from an ordered contour packing (labels of each half sorted by
   the y of their stubs, packed with no overlap beside today's label ellipse); slots in conflict are dropped from
   the reservation (greedy vertex cover), so the reserved layout is conflict free.
2. Labels are then placed one by one (most crowded first): a label releases its slot and takes the cheapest free box
   (leader length + tilt + drift + a penalty for blocking the radial corridors of unplaced stubs) that overlaps no box,
   puts no leader through a box and crosses no leader; every other label's slot or position is an obstacle.  If no
   free box is cheaper than its slot it keeps the slot.
3. Labels with neither get one more search that may cross at most `maxcross` leaders; otherwise they are dropped.
"""

from __future__ import annotations

import math

import numpy as np

from .common import Ctx, Item, angle_diff, emit, pack_ordered, polar

FAR = 1.0e7


def _seg_box_hit(p, q, boxes, pad):
    """p, q: (k,2) segments; boxes: (m,4).  Returns (k,m) bool: segment passes through the padded box."""
    p = p.astype(np.float32)
    q = q.astype(np.float32)
    boxes = boxes.astype(np.float32)
    pad = np.float32(pad)
    x0 = boxes[None, :, 0] - pad
    y0 = boxes[None, :, 1] - pad
    x1 = boxes[None, :, 2] + pad
    y1 = boxes[None, :, 3] + pad
    px = p[:, 0:1]
    py = p[:, 1:2]
    dx = q[:, 0:1] - px
    dy = q[:, 1:2] - py
    with np.errstate(divide="ignore", invalid="ignore"):
        tx0 = (x0 - px) / dx
        tx1 = (x1 - px) / dx
        ty0 = (y0 - py) / dy
        ty1 = (y1 - py) / dy
    zx = np.abs(dx) < 1e-12
    zy = np.abs(dy) < 1e-12
    tmin = np.maximum(np.minimum(tx0, tx1), np.minimum(ty0, ty1))
    tmax = np.minimum(np.maximum(tx0, tx1), np.maximum(ty0, ty1))
    # degenerate axes: inside the slab -> unconstrained, outside -> miss
    inx = (px > x0) & (px < x1)
    iny = (py > y0) & (py < y1)
    tmin = np.where(zx, np.where(inx, np.maximum(0.0, np.minimum(ty0, ty1)), np.inf), tmin)
    tmax = np.where(zx, np.where(inx, np.minimum(1.0, np.maximum(ty0, ty1)), -np.inf), tmax)
    tmin = np.where(zy, np.where(iny, np.maximum(0.0, np.minimum(tx0, tx1)), np.inf), tmin)
    tmax = np.where(zy, np.where(iny, np.minimum(1.0, np.maximum(tx0, tx1)), -np.inf), tmax)
    return (np.maximum(tmin, 0.0) < np.minimum(tmax, 1.0))


def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _seg_cross_count(p, q, rp, rq):
    """p,q: (k,2) new segments; rp,rq: (m,2) placed segments.  Returns (k,) number of proper crossings."""
    d1 = _orient(rp[None, :, 0], rp[None, :, 1], rq[None, :, 0], rq[None, :, 1], p[:, 0:1], p[:, 1:2])
    d2 = _orient(rp[None, :, 0], rp[None, :, 1], rq[None, :, 0], rq[None, :, 1], q[:, 0:1], q[:, 1:2])
    d3 = _orient(p[:, 0:1], p[:, 1:2], q[:, 0:1], q[:, 1:2], rp[None, :, 0], rp[None, :, 1])
    d4 = _orient(p[:, 0:1], p[:, 1:2], q[:, 0:1], q[:, 1:2], rq[None, :, 0], rq[None, :, 1])
    hit = ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))
    for d in (d1, d2, d3, d4):
        hit &= np.abs(d) > 1e-6
    return hit.sum(axis=1)


def _seg_cross_matrix(p, q, rp, rq):
    d1 = _orient(rp[None, :, 0], rp[None, :, 1], rq[None, :, 0], rq[None, :, 1], p[:, 0:1], p[:, 1:2])
    d2 = _orient(rp[None, :, 0], rp[None, :, 1], rq[None, :, 0], rq[None, :, 1], q[:, 0:1], q[:, 1:2])
    d3 = _orient(p[:, 0:1], p[:, 1:2], q[:, 0:1], q[:, 1:2], rp[None, :, 0], rp[None, :, 1])
    d4 = _orient(p[:, 0:1], p[:, 1:2], q[:, 0:1], q[:, 1:2], rq[None, :, 0], rq[None, :, 1])
    hit = ((d1 > 0) != (d2 > 0)) & ((d3 > 0) != (d4 > 0))
    for d in (d1, d2, d3, d4):
        hit &= np.abs(d) > 1e-6
    return hit


def _attach(sx, sy, ex0, ey0, ex1, ey1, corner):
    """Nearest edge point of the (expanded) box to the stub, slid away from corners."""
    gapx = np.maximum(ex0 - sx, sx - ex1)
    gapy = np.maximum(ey0 - sy, sy - ey1)
    lateral = gapx > gapy
    cy = np.clip(sy, ey0 + corner, ey1 - corner)
    cx = np.clip(sx, ex0 + corner, ex1 - corner)
    ax = np.where(lateral, np.where(sx < ex0, ex0, ex1), cx)
    ay = np.where(lateral, cy, np.where(sy < ey0, ey0, ey1))
    return ax, ay


CELL = 2.0
OFFC = 1000.0
NC = int(2 * OFFC / CELL)


def _cells(x0, y0, x1, y1):
    i0 = np.clip(np.floor((x0 + OFFC) / CELL).astype(int), 0, NC - 1)
    i1 = np.clip(np.floor((x1 + OFFC) / CELL).astype(int), 0, NC - 1)
    j0 = np.clip(np.floor((y0 + OFFC) / CELL).astype(int), 0, NC - 1)
    j1 = np.clip(np.floor((y1 + OFFC) / CELL).astype(int), 0, NC - 1)
    return i0, j0, i1, j1


def _sat(occ):
    sat = np.zeros((NC + 1, NC + 1), dtype=np.int32)
    sat[1:, 1:] = np.cumsum(np.cumsum(occ > 0, axis=0, dtype=np.int32), axis=1, dtype=np.int32)
    return sat


def _rect_sum(sat, i0, j0, i1, j1):
    return sat[j1 + 1, i1 + 1] - sat[j0, i1 + 1] - sat[j1 + 1, i0] + sat[j0, i0]


class _State:
    def __init__(self, items: list[Item], ctx: Ctx):
        self.items = items
        self.ctx = ctx
        n = len(items)
        self.boxes = np.full((n, 4), FAR)
        # two segments per label: anchor->stub (A) and stub->attach (B)
        self.sp = np.full((2 * n, 2), FAR)
        self.sq = np.full((2 * n, 2), FAR + 1.0)
        self.placed = [None] * n          # (box, attach) or None
        self.stub = [polar(ctx.r_elbow, it.theta) for it in items]
        self.r_gap = np.zeros((NC, NC), dtype=np.int16)    # placed boxes grown by the label gap
        self.r_pad = np.zeros((NC, NC), dtype=np.int16)    # placed boxes grown by the leader pad
        self.r_ldr = np.zeros((NC, NC), dtype=np.int16)    # placed leaders (stub -> box)
        self.marks: dict[int, tuple] = {}
        self.sat_gap = self.sat_ldr = None
        self.dirty = True
        # fixed radial segments (anchor -> stub) are known up front but only become obstacles when placed

    def _mark(self, i, box, attach, sign):
        g, pad = self.ctx.gap, 1.0
        for raster, e in ((self.r_gap, g), (self.r_pad, pad)):
            i0, j0, i1, j1 = _cells(box[0] - e, box[1] - e, box[2] + e, box[3] + e)
            raster[j0 : j1 + 1, i0 : i1 + 1] += sign
        s = self.stub[i]
        d = math.hypot(attach[0] - s[0], attach[1] - s[1])
        t = np.linspace(0.0, 1.0, int(d) + 2)
        px = s[0] + t * (attach[0] - s[0])
        py = s[1] + t * (attach[1] - s[1])
        ci = np.clip(np.floor((px + OFFC) / CELL).astype(int), 0, NC - 1)
        cj = np.clip(np.floor((py + OFFC) / CELL).astype(int), 0, NC - 1)
        flat = np.unique(cj * NC + ci)
        self.r_ldr.reshape(-1)[flat] += sign
        self.dirty = True

    def build(self):
        if self.dirty:
            self.sat_gap = _sat(self.r_gap)
            self.sat_ldr = _sat(self.r_ldr)
            self.occ_pad = self.r_pad > 0
            self.dirty = False

    def prefilter(self, b, at, s):
        """Cheap raster test (may reject by a couple of pixels too much, never accepts a collision)."""
        i0, j0, i1, j1 = _cells(b[:, 0], b[:, 1], b[:, 2], b[:, 3])
        ok = _rect_sum(self.sat_gap, i0, j0, i1, j1) == 0
        ok &= _rect_sum(self.sat_ldr, i0, j0, i1, j1) == 0
        if ok.any():
            sel = np.nonzero(ok)[0]
            a = at[sel]
            d = np.hypot(a[:, 0] - s[0], a[:, 1] - s[1])
            ns = int(d.max() / 4.0) + 2
            t = np.linspace(0.0, 1.0, ns)[None, :]
            px = s[0] + t * (a[:, 0:1] - s[0])
            py = s[1] + t * (a[:, 1:2] - s[1])
            ci = np.clip(np.floor((px + OFFC) / CELL).astype(int), 0, NC - 1)
            cj = np.clip(np.floor((py + OFFC) / CELL).astype(int), 0, NC - 1)
            hit = self.occ_pad[cj, ci].any(axis=1)
            ok[sel[hit]] = False
        return ok

    def remove(self, i):
        if self.placed[i] is not None:
            self._mark(i, self.placed[i][0], self.placed[i][1], -1)
        self.boxes[i] = FAR
        self.sp[2 * i : 2 * i + 2] = FAR
        self.sq[2 * i : 2 * i + 2] = FAR + 1.0
        self.placed[i] = None

    def put(self, i, box, attach):
        it = self.items[i]
        self.boxes[i] = box
        s = self.stub[i]
        self.sp[2 * i] = it.anchor
        self.sq[2 * i] = s
        self.sp[2 * i + 1] = s
        self.sq[2 * i + 1] = attach
        self.placed[i] = (box, attach)
        self._mark(i, box, attach, 1)


def _candidates(st: _State, i: int, lmax: float, step: float, drift: float, tilt: float, amax_deg: float = 75.0, band: float = 0.0):
    cos_amax = math.cos(math.radians(amax_deg))
    ctx = st.ctx
    it = st.items[i]
    w, h = it.w, it.h
    sx, sy = st.stub[i]
    pad, corner = 2.5, 4.0
    rx = lmax + 0.5 * w + pad
    ry = lmax + 0.5 * h + pad
    gx = np.arange(-rx, rx + 1e-9, step)
    gy = np.arange(-ry, ry + 1e-9, step)
    cx = (sx + gx)[None, :].repeat(len(gy), 0).ravel()
    cy = (sy + gy)[:, None].repeat(len(gx), 1).ravel()
    x0, x1 = cx - 0.5 * w, cx + 0.5 * w
    y0, y1 = cy - 0.5 * h, cy + 0.5 * h
    dxm = np.maximum(np.maximum(x0, -x1), 0.0)
    dym = np.maximum(np.maximum(y0, -y1), 0.0)
    ok = np.hypot(dxm, dym) >= ctx.r_keep + band
    ax, ay = _attach(sx, sy, x0 - pad, y0 - pad, x1 + pad, y1 + pad, corner)
    d = np.hypot(ax - sx, ay - sy)
    ok &= d <= lmax
    # the leader must leave the stub outward: at most `amax` from the radial direction (never dips into the ring zone)
    ux, uy = math.sin(it.theta), -math.cos(it.theta)
    vx, vy = ax - sx, ay - sy
    cos_a = (vx * ux + vy * uy) / np.maximum(d, 1e-9)
    ok &= (cos_a >= cos_amax) | (d < 2.0)
    idx = np.nonzero(ok)[0]
    cx, cy, x0, y0, x1, y1, ax, ay, d = (a[idx] for a in (cx, cy, x0, y0, x1, y1, ax, ay, d))
    ang = np.arctan2(cx, -cy)
    dth = np.abs((ang - it.theta + np.pi) % (2 * np.pi) - np.pi)
    # tilt of the leader against the radial direction at the stub
    cosang = ((ax - sx) * ux + (ay - sy) * uy) / np.maximum(d, 1e-9)
    tiltpen = np.where(d > 2.0, 1.0 - np.clip(cosang, -1, 1), 0.0)
    cost = d + drift * ctx.radius * dth + tilt * d * tiltpen
    order = np.argsort(cost, kind="stable")
    box = np.stack([x0, y0, x1, y1], axis=1)[order]
    att = np.stack([ax, ay], axis=1)[order]
    return box, att, cost[order]


def _rays(st: _State, i: int, depth: float, near=None):
    """Radial corridors (stub -> stub + depth outward) of the labels not placed yet (except i)."""
    cand = range(len(st.items)) if near is None else near
    idx = [k for k in cand if k != i and st.placed[k] is None]
    if not idx:
        return None
    p = np.array([st.stub[k] for k in idx])
    u = np.array([(math.sin(st.items[k].theta), -math.cos(st.items[k].theta)) for k in idx])
    return p, p + depth * u


def _best(st: _State, i: int, cand, cross_pen: float, ray_pen: float = 0.0, ray_depth: float = 100.0,
          block: int = 2000, chunk: int = 128, reach: float = 1.3, max_cross: int = 99):
    ctx = st.ctx
    box, att, cost = cand
    s = np.array(st.stub[i])
    best = None
    best_total = math.inf
    st.build()
    th = np.array([it.theta for it in st.items])
    near = np.nonzero(np.abs((th - st.items[i].theta + np.pi) % (2 * np.pi) - np.pi) < reach)[0]
    near_seg = np.concatenate([2 * near, 2 * near + 1])
    nb, nsp, nsq = st.boxes[near], st.sp[near_seg], st.sq[near_seg]
    rays = _rays(st, i, ray_depth, near) if ray_pen > 0 else None
    ray_cap = int(ctx.opt('raycap', 8))
    g = ctx.gap
    for a0 in range(0, len(cost), block):
        if cost[a0] >= best_total:
            break
        keep = np.nonzero(st.prefilter(box[a0 : a0 + block], att[a0 : a0 + block], s))[0] + a0
        for a in range(0, len(keep), chunk):
            kk = keep[a : a + chunk]
            if cost[kk[0]] >= best_total:
                break
            b, at, c = box[kk], att[kk], cost[kk]
            ov = (
                (b[:, None, 0] < nb[None, :, 2] + g)
                & (nb[None, :, 0] < b[:, None, 2] + g)
                & (b[:, None, 1] < nb[None, :, 3] + g)
                & (nb[None, :, 1] < b[:, None, 3] + g)
            )
            ok = ~ov.any(axis=1)
            if not ok.any():
                continue
            sel = np.nonzero(ok)[0]
            b, at, c = b[sel], at[sel], c[sel]
            sp = np.broadcast_to(s, at.shape)
            hit = _seg_box_hit(sp, at, nb, 1.0).any(axis=1)
            hit |= _seg_box_hit(nsp, nsq, b, 1.0).any(axis=0)
            okl = ~hit
            if not okl.any():
                continue
            sel = np.nonzero(okl)[0]
            b, at, c = b[sel], at[sel], c[sel]
            sp = np.broadcast_to(s, at.shape)
            nc = _seg_cross_count(sp, at, nsp, nsq)
            if rays is not None:
                nr = _seg_box_hit(rays[0], rays[1], b, 1.0).sum(axis=0) + _seg_cross_count(sp, at, rays[0], rays[1])
                c = c + ray_pen * np.minimum(nr, ray_cap)
            if cross_pen == math.inf:
                good = np.nonzero(nc == 0)[0]
                if len(good):
                    k = good[np.argmin(c[good])]
                    if c[k] < best_total:
                        best_total = float(c[k])
                        best = (b[k], at[k], float(c[k]), 0)
                continue
            tot = c + cross_pen * nc + np.where(nc > max_cross, 1e9, 0.0)
            k = int(np.argmin(tot))
            if tot[k] >= 1e8:
                continue
            if tot[k] < best_total:
                best_total = float(tot[k])
                best = (b[k], at[k], float(c[k]), int(nc[k]))
    return best


def _crowd_order(items, width):
    """Most crowded labels first (number of labels within +-width rad), ties by clockwise index."""
    n = len(items)
    crowd = [sum(1 for m in range(n) if abs(angle_diff(items[m].theta, items[k].theta)) < width) for k in range(n)]
    return sorted(range(n), key=lambda k: (-crowd[k], k))


def _conflicts(ctx, B, S, A):
    """Pairwise conflicts of a full layout: boxes overlap, a leader runs through a box, leaders cross."""
    g = ctx.gap
    conf = (
        (B[:, None, 0] < B[None, :, 2] + g) & (B[None, :, 0] < B[:, None, 2] + g)
        & (B[:, None, 1] < B[None, :, 3] + g) & (B[None, :, 1] < B[:, None, 3] + g)
    )
    hit = _seg_box_hit(S, A, B, 1.0)
    conf |= hit | hit.T
    conf |= _seg_cross_matrix(S, A, S, A)
    np.fill_diagonal(conf, False)
    return conf


def _contour_targets(items, ctx, ax_r, ay_r, r_min):
    """Seed positions: per half of the map, labels in the order of their stubs' y, packed without overlap
    (optimal ordered packing) next to an ellipse of radii (ax_r, ay_r); right half starts at the ellipse, left ends."""
    targets = np.zeros((len(items), 2))
    for half in (0, 1):
        idx = [k for k in range(len(items)) if (math.sin(items[k].theta) >= 0) == (half == 0)]
        idx.sort(key=lambda k: -math.cos(items[k].theta))      # top to bottom (y = -cos)
        if not idx:
            continue
        want = [-ctx.r_elbow * math.cos(items[k].theta) for k in idx]
        ys = pack_ordered(want, [items[k].h + ctx.gap for k in idx])
        for k, y in zip(idx, ys):
            it = items[k]
            ex = ax_r * math.sqrt(max(0.0, 1.0 - min(abs(y) / ay_r, 1.0) ** 2))
            ex = max(ex, 0.55 * ax_r)
            dy = max(abs(y) - 0.5 * it.h, 0.0)
            ex = max(ex, math.sqrt(max(r_min * r_min - dy * dy, 0.0)))   # box stays outside r_keep + band
            cx = ex + 0.5 * it.w if half == 0 else -ex - 0.5 * it.w
            targets[k] = (cx, y)
    return targets


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    n = len(items)
    lmax = ctx.opt("lmax", 0.7) * ctx.radius          # bound on the leader length (stub -> box edge)
    step = ctx.opt("step", 3.0 if n <= 130 else 4.0)  # candidate grid pitch (px)
    drift = ctx.opt("drift", 0.1)                     # cost per px of arc drift of the box centre
    tilt = ctx.opt("tilt", 0.5)                       # cost of leader tilt away from the radial
    amax = ctx.opt("amax", 85.0)                      # largest angle between leader and radial (deg)
    ray_pen = ctx.opt("raypen", 100.0)                # cost per unplaced stub corridor a box/leader blocks
    ray_depth = ctx.opt("raydepth", 60.0)             # corridor depth (px)
    cross_pen = ctx.opt("crosspen", 40.0)             # cost per crossing in the fallback pass
    max_cross = int(ctx.opt("maxcross", 2))           # a fallback leader may cross at most this many others
    band = ctx.opt("band", 12.0)                      # extra keep-out band outside r_keep (leader lane)
    fallback = int(ctx.opt("fallback", 1))            # 0: never cross, drop instead
    slotbias = ctx.opt("slotbias", 0.0)               # extra cost a free position must beat the reserved slot by
    reach = (lmax + 280.0) / ctx.r_elbow
    st = _State(items, ctx)
    # ---- 1. reserve a conflict-free fallback slot per label (ordered contour packing, today's kind of layout)
    reserved = np.zeros(n, dtype=bool)
    slot: dict[int, tuple] = {}
    if int(ctx.opt("slots", 1)):
        T = _contour_targets(items, ctx, ctx.arc_x_radius, ctx.arc_y_radius, ctx.r_keep + band)
        hw = np.array([0.5 * it.w for it in items])
        hh = np.array([0.5 * it.h for it in items])
        B0 = np.stack([T[:, 0] - hw, T[:, 1] - hh, T[:, 0] + hw, T[:, 1] + hh], axis=1)
        S0 = np.array(st.stub)
        ax, ay = _attach(S0[:, 0], S0[:, 1], B0[:, 0] - 2.5, B0[:, 1] - 2.5, B0[:, 2] + 2.5, B0[:, 3] + 2.5, 4.0)
        A0 = np.stack([ax, ay], axis=1)
        length0 = np.hypot(A0[:, 0] - S0[:, 0], A0[:, 1] - S0[:, 1])
        alive = length0 <= lmax
        conf = _conflicts(ctx, B0, S0, A0)
        while True:
            cnt = (conf & alive[None, :] & alive[:, None]).sum(axis=1)
            if cnt.max() == 0:
                break
            worst = min(np.nonzero(cnt == cnt.max())[0], key=lambda k: (items[k].span_bp, k))
            alive[worst] = False
        for i in np.nonzero(alive)[0]:
            reserved[i] = True
            slot[int(i)] = (B0[i], A0[i], float(length0[i]))
            st.put(int(i), B0[i], A0[i])
        ctx.notes["slots_reserved"] = int(alive.sum())
    # ---- 2. free placement, one label at a time; a label first gives up its own slot, then takes the cheapest
    #         collision-free box (every other reserved slot is an obstacle); its slot is the fallback
    om = int(ctx.opt("order", 4))
    if om == 4:
        order = _crowd_order(items, ctx.opt("crowdw", 0.12))
    else:
        start = max(range(n), key=lambda k: items[k].span_bp)
        order = [(start + k) % n for k in range(n)]
    failed = []
    n_free = 0
    for i in order:
        if reserved[i]:
            st.remove(i)
        cd = _candidates(st, i, lmax, step, drift, tilt, amax, band=band)
        r = _best(st, i, cd, math.inf, ray_pen, ray_depth, reach=reach)
        if reserved[i]:
            b0, a0, d0 = slot[i]
            if r is None or r[2] > d0 + slotbias:
                st.put(i, b0, a0)
                continue
            n_free += 1
        if r is None:
            failed.append(i)
        else:
            st.put(i, r[0], r[1])
    ctx.notes["free_positions"] = n_free
    ctx.notes["first_pass_failed"] = len(failed)
    # ---- 3. labels that found neither: at most `max_cross` crossings, else dropped
    dropped = []
    for i in failed:
        cd = _candidates(st, i, lmax, step, drift, tilt, amax, band=band)
        r = _best(st, i, cd, cross_pen, 0.0, ray_depth, reach=reach, max_cross=max_cross) if fallback else None
        if r is not None:
            ctx.notes["crossing_fallback"] = ctx.notes.get("crossing_fallback", 0) + 1
            st.put(i, r[0], r[1])
        else:
            dropped.append(i)
    for k in dropped:
        ctx.dropped.append(items[k].text)
    ctx.notes["dropped_idx"] = list(dropped)
    out = []
    for i, it in enumerate(items):
        if st.placed[i] is None:
            continue
        box, att = st.placed[i]
        out.append(emit(it, tuple(float(v) for v in box), [it.anchor, st.stub[i], (float(att[0]), float(att[1]))], align="middle"))
    ctx.notes["lmax_px"] = round(lmax, 1)
    return out
