"""Design G: liberal placement on both sides of the ring (the Geneious Prime look); study prototype.

Every label is a horizontal box.  A label's first choice is the spot right next to its feature, outside the ring or,
where the centre is free, inside it.  Labels are placed one at a time (most crowded stretch first); for each label
candidate boxes on a grid around its two stubs (the point on the feature's ray at the ring's outer edge, and at its inner
edge) are tried in order of cost (leader length + angular drift + tilt + a small inside penalty) and the first that
collides with nothing is taken.  A leader is straight: stub -> nearest point of the box (the hidden part anchor -> stub is
radial).  Inside the ring no box may touch the tracks, the tick marks, the tick labels or the title; the radii come from
the radial layout (ctx.radial_layout).  Then a few refinement passes re-place each label at its cheapest feasible spot.
Options: inside=0|1, merge=0|1, color=0|1 (see place()).
"""

from __future__ import annotations

import math

import numpy as np

from ...auto_sizes import circular_tick_intervals
from ...core.numeric import scaled_tick_text
from ...core.text import calculate_bbox_dimensions
from .common import Ctx, Item, angle_diff, emit, polar, unit

FAR = 1.0e7


# ------------------------------------------------------------------ vectorised geometry (copied from design C)

def _seg_box_hit(p, q, boxes, pad):
    """p, q: (k,2) segments; boxes: (m,4).  (k,m) bool: segment passes through the padded box."""
    x0 = boxes[None, :, 0] - pad
    y0 = boxes[None, :, 1] - pad
    x1 = boxes[None, :, 2] + pad
    y1 = boxes[None, :, 3] + pad
    px, py = p[:, 0:1], p[:, 1:2]
    dx, dy = q[:, 0:1] - px, q[:, 1:2] - py
    with np.errstate(divide="ignore", invalid="ignore"):
        tx0, tx1 = (x0 - px) / dx, (x1 - px) / dx
        ty0, ty1 = (y0 - py) / dy, (y1 - py) / dy
    zx, zy = np.abs(dx) < 1e-12, np.abs(dy) < 1e-12
    tmin = np.maximum(np.minimum(tx0, tx1), np.minimum(ty0, ty1))
    tmax = np.minimum(np.maximum(tx0, tx1), np.maximum(ty0, ty1))
    inx, iny = (px > x0) & (px < x1), (py > y0) & (py < y1)
    tmin = np.where(zx, np.where(inx, np.maximum(0.0, np.minimum(ty0, ty1)), np.inf), tmin)
    tmax = np.where(zx, np.where(inx, np.minimum(1.0, np.maximum(ty0, ty1)), -np.inf), tmax)
    tmin = np.where(zy, np.where(iny, np.maximum(0.0, np.minimum(tx0, tx1)), np.inf), tmin)
    tmax = np.where(zy, np.where(iny, np.minimum(1.0, np.maximum(tx0, tx1)), -np.inf), tmax)
    return np.maximum(tmin, 0.0) < np.minimum(tmax, 1.0)


def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


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
    """Nearest edge point of the box to the stub, slid away from corners."""
    gapx = np.maximum(ex0 - sx, sx - ex1)
    gapy = np.maximum(ey0 - sy, sy - ey1)
    lateral = gapx > gapy
    cy = np.clip(sy, ey0 + corner, ey1 - corner)
    cx = np.clip(sx, ex0 + corner, ex1 - corner)
    ax = np.where(lateral, np.where(sx < ex0, ex0, ex1), cx)
    ay = np.where(lateral, cy, np.where(sy < ey0, ey0, ey1))
    return ax, ay


def _radii(x0, y0, x1, y1):
    """Smallest and largest distance from the centre to the boxes."""
    dxm = np.maximum(np.maximum(x0, -x1), 0.0)
    dym = np.maximum(np.maximum(y0, -y1), 0.0)
    fx = np.maximum(np.abs(x0), np.abs(x1))
    fy = np.maximum(np.abs(y0), np.abs(y1))
    return np.hypot(dxm, dym), np.hypot(fx, fy)


def _seg_radii(p, a):
    """Smallest and largest distance from the centre to the segments p-a (rows)."""
    d = a - p
    l2 = np.maximum((d * d).sum(axis=1), 1e-12)
    t = np.clip(-(p * d).sum(axis=1) / l2, 0.0, 1.0)
    near = p + t[:, None] * d
    return np.hypot(near[:, 0], near[:, 1]), np.maximum(np.hypot(p[:, 0], p[:, 1]), np.hypot(a[:, 0], a[:, 1]))


# ------------------------------------------------------------------ what the ring and the inside look like

class _Zones:
    """Radii of the ring, the forbidden annuli (tracks, tick marks, title) and the tick-label boxes."""

    def __init__(self, items: list[Item], ctx: Ctx, clear: float):
        rl = ctx.radial_layout
        feats = rl.features if rl is not None else None
        if feats is not None:
            self.r_in, self.r_out = float(feats.all_band_px.inner_px), float(feats.all_band_px.outer_px)
        else:
            self.r_out = max(i.feat_outer_r for i in items)
            self.r_in = min(float(i.src.get("feature_inner_radius_px", i.feat_outer_r)) for i in items)
        self.annuli: list[tuple[float, float, bool]] = []     # (inner, outer, crossable by a leader)
        self.ticklabels: list[tuple[float, float, float, float]] = []
        self.source = "radial_layout" if rl is not None else "items (no inner tracks known)"
        self.tick_font = None
        self.inside_hi = self.r_in - clear           # largest radius an inside box may reach
        if rl is None:
            return
        for slot in rl.slots:
            if slot.renderer in ("features", "feature_labels"):
                continue
            band = slot.reserved_band_px
            if slot.renderer == "ticks" and slot.payload is not None:
                tick = slot.payload
                self.annuli.append((float(tick.tick_band_px.inner_px), float(tick.tick_band_px.outer_px), True))
                if tick.label_band_px is not None:
                    self._tick_labels(items, ctx, tick.label_band_px, int(ctx.total_length))
            elif band is not None and band.width_px > 1e-6:
                self.annuli.append((float(band.inner_px), float(band.outer_px), False))
        if rl.definition_reserved_band_px is not None:
            self.annuli.append((0.0, float(rl.definition_reserved_band_px.outer_px), False))
        self.annuli.append((self.r_in, self.r_out, True))

    def _tick_labels(self, items, ctx, band, total):
        """Major tick labels ('5 kbp' ...) as chains of small boxes along their tangent (auto tick interval)."""
        large, _small = circular_tick_intervals(total)
        h = float(band.outer_px - band.inner_px)
        sample = items[0].src
        fs = h / max(items[0].h / max(items[0].src.get("font_size", ctx.font_size), 1.0), 1.0)
        self.tick_font = fs
        rc = 0.5 * (band.inner_px + band.outer_px)
        for tick in range(large, total, large):
            text = f"{scaled_tick_text(tick, 1000, large)} kbp" if total < 1_000_000 else f"{tick / 1000000} Mbp"
            w, _ = calculate_bbox_dimensions(text, sample.get("font_family", ""), fs, 96)
            theta = 2.0 * math.pi * tick / total
            ux, uy = math.sin(theta), -math.cos(theta)
            tx, ty = math.cos(theta), math.sin(theta)
            n = max(1, int(math.ceil(w / h)))
            seg = w / n
            for k in range(n):
                off = (k + 0.5) * seg - 0.5 * w
                cx, cy = rc * ux + off * tx, rc * uy + off * ty
                ex = 0.5 * seg * abs(tx) + 0.5 * h * abs(ux)
                ey = 0.5 * seg * abs(ty) + 0.5 * h * abs(uy)
                self.ticklabels.append((cx - ex, cy - ey, cx + ex, cy + ey))


# ------------------------------------------------------------------ the placement state

class _State:
    def __init__(self, items, ctx, zones: _Zones, inside: bool, clear: float):
        self.items, self.ctx, self.z, self.clear = items, ctx, zones, clear
        n = len(items)
        self.has_inside = inside and zones.inside_hi > 40.0
        self.sides = (0, 1) if self.has_inside else (0,)
        # stubs[i][side], dirs[i][side] (unit vector the leader leaves along)
        self.stub, self.udir = [], []
        for it in items:
            ux, uy = unit(it.theta)
            self.stub.append([polar(zones.r_out, it.theta), polar(zones.r_in, it.theta)])
            self.udir.append([(ux, uy), (-ux, -uy)])
        self.members: list[list[int]] = [[k] for k in range(n)]    # item indices a label stands for (merge=1)
        self.placed: list[tuple | None] = [None] * n               # (box, attach, side)
        self.fixed = [tuple(b) for b in zones.ticklabels]
        self.dirty = True

    # --- obstacle / leader arrays
    def build(self):
        if not self.dirty:
            return
        boxes = list(self.fixed)
        sp, sq = [], []
        for i, pl in enumerate(self.placed):
            if pl is None:
                continue
            boxes.append(pl[0])
            for k in self.members[i]:
                sp.append(self.stub[k][pl[2]])
                sq.append(pl[1])
        self.boxes = np.array(boxes, dtype=float).reshape(-1, 4)
        self.sp = np.array(sp, dtype=float).reshape(-1, 2)
        self.sq = np.array(sq, dtype=float).reshape(-1, 2)
        self.dirty = False

    def put(self, i, box, attach, side):
        self.placed[i] = (tuple(float(v) for v in box), (float(attach[0]), float(attach[1])), int(side))
        self.dirty = True

    def remove(self, i):
        self.placed[i] = None
        self.dirty = True

    # --- candidates
    def candidates(self, i, L, amax_deg, drift, tilt, inpen, step):
        it = self.items[i]
        ctx, z = self.ctx, self.z
        cos_amax = math.cos(math.radians(amax_deg))
        w, h = it.w, it.h
        out = []
        for side in self.sides:
            sx, sy = self.stub[i][side]
            ux, uy = self.udir[i][side]
            rx, ry = L + 0.5 * w, L + 0.5 * h
            gx = np.arange(-rx, rx + 1e-9, step)
            gy = np.arange(-ry, ry + 1e-9, step)
            cx = (sx + gx)[None, :].repeat(len(gy), 0).ravel()
            cy = (sy + gy)[:, None].repeat(len(gx), 1).ravel()
            x0, x1, y0, y1 = cx - 0.5 * w, cx + 0.5 * w, cy - 0.5 * h, cy + 0.5 * h
            rmin, rmax = _radii(x0, y0, x1, y1)
            ok = rmin >= z.r_out + self.clear if side == 0 else rmax <= z.inside_hi
            for a, b, _cross in z.annuli:
                ok &= ~((rmin < b + self.clear) & (rmax > a - self.clear))
            idx = np.nonzero(ok)[0]
            if not len(idx):
                continue
            cx, cy, x0, y0, x1, y1 = (v[idx] for v in (cx, cy, x0, y0, x1, y1))
            ax, ay = _attach(sx, sy, x0, y0, x1, y1, 4.0)
            d = np.hypot(ax - sx, ay - sy)
            cosang = ((ax - sx) * ux + (ay - sy) * uy) / np.maximum(d, 1e-9)
            keep = (d <= L) & ((cosang >= cos_amax) | (d < 2.0))
            # the leader may not cut into the ring or into a track / title it does not start from
            a_xy = np.stack([ax, ay], axis=1)
            p_xy = np.broadcast_to(np.array([sx, sy]), a_xy.shape)
            lmin, lmax = _seg_radii(p_xy, a_xy)
            keep &= (lmin >= z.r_out - 0.5) if side == 0 else (lmax <= z.r_in + 0.5)
            for a, b, cross in z.annuli:
                if not cross:
                    keep &= ~((lmin < b) & (lmax > a))
            keep = np.nonzero(keep)[0]
            if not len(keep):
                continue
            cx, cy, x0, y0, x1, y1, ax, ay, d, cosang = (v[keep] for v in (cx, cy, x0, y0, x1, y1, ax, ay, d, cosang))
            dth = np.abs((np.arctan2(cx, -cy) - it.theta + np.pi) % (2 * np.pi) - np.pi)
            cost = d + drift * ctx.radius * dth + tilt * d * (1.0 - np.clip(cosang, -1, 1)) + (inpen if side else 0.0)
            out.append((cost, np.stack([x0, y0, x1, y1], axis=1), np.stack([ax, ay], axis=1), np.full(len(cost), side)))
        if not out:
            return None
        cost = np.concatenate([o[0] for o in out])
        box = np.concatenate([o[1] for o in out])
        att = np.concatenate([o[2] for o in out])
        side = np.concatenate([o[3] for o in out])
        order = np.argsort(cost, kind="stable")
        return cost[order], box[order], att[order], side[order]

    def rays(self, i, depth):
        """Radial corridors (stub -> stub + depth, away from the ring) of the labels not placed yet, per side."""
        idx = [k for k in range(len(self.items)) if self.placed[k] is None and k != i]
        res = {}
        for side in self.sides:
            if not idx:
                res[side] = None
                continue
            p = np.array([self.stub[k][side] for k in idx])
            u = np.array([self.udir[k][side] for k in idx])
            res[side] = (p, p + depth * u)
        return res

    def best(self, i, cand, ray_pen, ray_depth, ray_cap, cross_ok, hard_depth=0.0, lpad=None, block=400):
        lpad = self.clear if lpad is None else lpad
        cost, box, att, side = cand
        g = self.clear
        self.build()
        rays = self.rays(i, ray_depth) if ray_pen > 0 else {}
        hard = self.rays(i, hard_depth) if hard_depth > 0 else {}
        best, best_total = None, math.inf
        for a0 in range(0, len(cost), block):
            if cost[a0] >= best_total:
                break
            b, at, c, sd = box[a0 : a0 + block], att[a0 : a0 + block], cost[a0 : a0 + block], side[a0 : a0 + block]
            ok = np.ones(len(c), bool)
            if len(self.boxes):
                B = self.boxes
                ok &= ~(
                    (b[:, None, 0] < B[None, :, 2] + g) & (B[None, :, 0] < b[:, None, 2] + g)
                    & (b[:, None, 1] < B[None, :, 3] + g) & (B[None, :, 1] < b[:, None, 3] + g)
                ).any(axis=1)
            if not ok.any():
                continue
            sel = np.nonzero(ok)[0]
            b, at, c, sd = b[sel], at[sel], c[sel], sd[sel]
            ok = np.ones(len(c), bool)
            ncross = np.zeros(len(c), int)
            if len(self.sp):
                ok &= ~_seg_box_hit(self.sp, self.sq, b, lpad).any(axis=0)
            for k in self.members[i]:
                P = np.where(sd[:, None] == 0, self.stub[k][0], self.stub[k][1])
                if len(self.boxes):
                    ok &= ~_seg_box_hit(P, at, self.boxes, lpad).any(axis=1)
                if len(self.sp):
                    ncross += _seg_cross_matrix(P, at, self.sp, self.sq).sum(axis=1)
            if not cross_ok:
                ok &= ncross == 0
            if not ok.any():
                continue
            sel = np.nonzero(ok)[0]
            b, at, c, sd, ncross = b[sel], at[sel], c[sel], sd[sel], ncross[sel]
            c = c + 60.0 * ncross
            ok = np.ones(len(c), bool)
            for side_v, r in hard.items():
                if r is None:
                    continue
                m = np.nonzero(sd == side_v)[0]
                if len(m):
                    ok[m] = ~_seg_box_hit(r[0], r[1], b[m], 0.0).any(axis=0)
            if not ok.any():
                continue
            sel = np.nonzero(ok)[0]
            b, at, c, sd, ncross = b[sel], at[sel], c[sel], sd[sel], ncross[sel]
            for side_v, r in rays.items():
                if r is None:
                    continue
                m = sd == side_v
                if not m.any():
                    continue
                hit = _seg_box_hit(r[0], r[1], b[m], 1.0).sum(axis=0)
                c[m] += ray_pen * np.minimum(hit, ray_cap)
            k = int(np.argmin(c))
            if c[k] < best_total:
                best_total = float(c[k])
                best = (b[k], at[k], int(sd[k]), best_total, int(ncross[k]))
        return best


def _crowd_order(items, width):
    n = len(items)
    crowd = [sum(1 for m in range(n) if abs(angle_diff(items[m].theta, items[k].theta)) < width) for k in range(n)]
    return sorted(range(n), key=lambda k: (-crowd[k], k))


def _merge_runs(items: list[Item], ctx: Ctx, gap_rad: float, span_rad: float) -> tuple[list[Item], list[list[Item]]]:
    """Runs of clockwise-adjacent items with equal text become one item 'text (xN)' standing at the run's middle."""
    runs: list[list[Item]] = []
    for it in items:
        if runs and runs[-1][-1].text == it.text and abs(angle_diff(it.theta, runs[-1][-1].theta)) <= gap_rad \
                and abs(angle_diff(it.theta, runs[-1][0].theta)) <= span_rad:
            runs[-1].append(it)
        else:
            runs.append([it])
    merged, groups = [], []
    for run in runs:
        if len(run) == 1:
            merged.append(run[0])
            groups.append(run)
            continue
        mid = run[len(run) // 2]
        text = f"{mid.text} (×{len(run)})"
        src = dict(mid.src)
        fam = src.get("font_family", "")
        w, h = calculate_bbox_dimensions(text, fam, float(src.get("font_size", ctx.font_size)), 96)
        src["label_text"], src["width_px"], src["height_px"] = text, float(w), float(h)
        merged.append(Item(index=mid.index, src=src, text=text, w=float(w), h=float(h), theta=mid.theta,
                           anchor=mid.anchor, feat_outer_r=mid.feat_outer_r, span_bp=mid.span_bp))
        groups.append(run)
    return merged, groups


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    """Options (ctx.opt): inside 1 | merge 0 | color 0 | clear 3 (= ctx.gap) | lmax 1.3 (x R) | step 3 | drift 0.1 | tilt 0.5 |
    inpen 4 | amax 75 | raypen 30 | raydepth 50 | passes 2 | crowdw 0.12 | mergegap 0.06 | mergespan 0.14 | cross 0."""
    n_requested = len(items)
    clear = ctx.opt("clear", ctx.gap)
    zones = _Zones(items, ctx, clear)
    groups: list[list[Item]] = [[it] for it in items]
    if int(ctx.opt("merge", 0)):
        items, groups = _merge_runs(items, ctx, ctx.opt("mergegap", 0.06), ctx.opt("mergespan", 0.14))
    st = _State(items, ctx, zones, bool(int(ctx.opt("inside", 1))), clear)
    # a merged label stands for several features: each gets a stub (same side) and a leader to the shared attach point
    for gi, run in enumerate(groups):
        if len(run) > 1:
            st.members[gi] = []
            for m in run:
                st.members[gi].append(len(st.stub))
                st.stub.append([polar(zones.r_out, m.theta), polar(zones.r_in, m.theta)])
                st.udir.append([unit(m.theta), tuple(-v for v in unit(m.theta))])
    drift, tilt, inpen = ctx.opt("drift", 0.1), ctx.opt("tilt", 0.5), ctx.opt("inpen", 4.0)
    step, amax = ctx.opt("step", 3.0), ctx.opt("amax", 75.0)
    ray_pen, ray_depth, ray_cap = ctx.opt("raypen", 30.0), ctx.opt("raydepth", 50.0), int(ctx.opt("raycap", 6))
    cross_ok = bool(int(ctx.opt("cross", 0)))
    hard_depth = ctx.opt("harddepth", 36.0)
    lmax = ctx.opt("lmax", 1.3) * ctx.radius
    levels = [v for v in (30.0, 60.0, 120.0, 240.0, 480.0) if v < lmax] + [lmax]

    def find(i, with_rays, cross=cross_ok, lpad=None):
        for L in levels:
            cand = st.candidates(i, L, amax if L < lmax else 85.0, drift, tilt, inpen, step)
            if cand is None:
                continue
            r = st.best(i, cand, ray_pen if with_rays else 0.0, ray_depth, ray_cap, cross, hard_depth if with_rays else 0.0, lpad)
            if r is not None:
                return r
        return None

    def attempt(order):
        st.placed = [None] * len(items)
        st.dirty = True
        failed = []
        for i in order:
            r = find(i, True)
            if r is None:
                failed.append(i)
            else:
                st.put(i, r[0], r[1], r[2])
        for _ in range(int(ctx.opt("passes", 2))):
            for i in order:
                old = st.placed[i]
                if old is None:
                    continue
                st.remove(i)
                r = find(i, False)
                if r is None:
                    st.placed[i] = old
                    st.dirty = True
                else:
                    st.put(i, r[0], r[1], r[2])
        gave_up, forced = [], 0
        for i in failed:
            r = find(i, False, bool(int(ctx.opt("fallback", 1))) or cross_ok, 1.0)   # last resort: no corridors, leader clearance 1 px, cheapest crossings
            forced += r is not None and r[4] > 0
            if r is None:
                gave_up.append(i)
            else:
                st.put(i, r[0], r[1], r[2])
        return gave_up, forced, len(failed)

    order = _crowd_order(items, ctx.opt("crowdw", 0.12))
    best_try = None
    for t in range(1 + int(ctx.opt("retries", 0))):
        gave_up, forced, n_failed = attempt(order)
        if best_try is None or len(gave_up) < len(best_try[0]):
            best_try = (gave_up, forced, n_failed, list(st.placed), t)
        if not gave_up:
            break
        order = gave_up + [k for k in order if k not in gave_up]    # the labels that found no room go first next time
    gave_up, forced, n_failed, placed_best, tries = best_try
    st.placed = placed_best
    failed = list(range(n_failed))
    for i in gave_up:
        ctx.dropped.append(items[i].text)
    use_color = int(ctx.opt("color", 0))
    out = []
    n_in = n_out = 0
    for i, it in enumerate(items):
        pl = st.placed[i]
        if pl is None:
            continue
        box, att, side = pl
        mem = st.members[i]
        pts = [it.anchor, st.stub[i][side], att]
        extra = []
        if len(mem) > 1:
            run = groups[i]
            mid = len(mem) // 2
            for j, k in enumerate(mem):
                if j != mid:
                    extra.append([run[j].anchor, st.stub[k][side], att])
        label = emit(it, box, pts, align="middle", extra=extra, tag="in" if side else "out")
        if use_color:
            fill = ctx.feature_colors.get(str(it.src.get("stable_id")))
            if fill:
                label["text_fill"] = fill
        out.append(label)
        n_in += side
        n_out += 1 - side
    ctx.notes.update(
        inside=n_in, outside=n_out, merged_labels=sum(1 for g in groups if len(g) > 1), features_requested=n_requested,
        zones=zones.source, r_in=round(zones.r_in, 1), r_out=round(zones.r_out, 1), inside_hi=round(zones.inside_hi, 1),
        annuli=[(round(a, 1), round(b, 1), c) for a, b, c in zones.annuli], tick_labels=len(zones.ticklabels),
        has_inside=st.has_inside, failed_first_pass=n_failed, tries=tries + 1, forced_crossing_labels=forced,
    )
    return out
