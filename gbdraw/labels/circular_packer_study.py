"""STUDY ONLY (branch study/circular-label-packer, never merged): an order-preserving,
minimum-displacement packer for Circular horizontal labels (OV-472/473/697).

Switches (environment, default off):
  GBDRAW_CLP=A   pack only records whose current placement does not fit
  GBDRAW_CLP=B   pack every record
  GBDRAW_CLP_OPTS=k=v,...
    rows=spaced|lp31|legacy  row gap: the circular label spacing (default, today's spaced attempt),
                             max(0.22 em, 4 px) (LP-D30/31), or 0.1 px (today's dense legacy rows);
                             "spaced" and "lp31" fall back to 0.1 px for a half that does not fit
    coord=y|arc              the coordinate of the least-squares fit: vertical position (default) or arc length
    spill=1|0                move labels across the poles when one half is fuller than the other (default 1)
    resolve=0|1              run today's leader-collision resolver after packing (default 0)

The labels keep their cyclic genome order on today's label ellipse (one arc, LP-D38/39), parameterised by
psi (anchor = (cx + a sin psi, -b cos psi), psi = 0 at 12 o'clock, clockwise). Each label's target is where
its feature's ray meets the ellipse; each label may move only inside its no-dip cone (the angles at which
its straight leader still leaves the ring outward, so leaders that keep genome order cannot cross,
LP-D36). The packer minimises sum (psi_i - t_i)^2 subject to the box clearance between neighbours
(rows: the row gap; side by side: 1 px, LP-D37) and the cones: pool-adjacent-violators on psi_i - D_i with
the neighbour steps D re-measured a few times, then bounded forward/backward sweeps. When the cones admit
no packing the record does not fit (crossings unavoidable with straight leaders); the drawn fallback widens
every cone by the smallest factor that packs.
"""

from __future__ import annotations

import math
import os
from typing import Any

LAST: dict[str, Any] = {}
ACTIVE = False


def mode() -> str:
    return os.environ.get("GBDRAW_CLP", "").strip().upper()


def opts() -> dict[str, str]:
    raw = os.environ.get("GBDRAW_CLP_OPTS", "")
    out = {"rows": "spaced", "coord": "y", "spill": "1", "resolve": "0"}
    for part in raw.split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _isotonic(values: list[float]) -> list[float]:
    """Least-squares non-decreasing fit (pool adjacent violators)."""
    blocks: list[list[float]] = []  # [sum, count]
    for v in values:
        blocks.append([v, 1.0])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s, c = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += c
    out: list[float] = []
    for s, c in blocks:
        out.extend([s / c] * int(c))
    return out


def _ray_target_y(label, cx: float, a: float, b: float) -> float:
    fx = float(label["feature_middle_x"])
    fy = float(label["feature_middle_y"])
    n = math.hypot(fx, fy) or 1.0
    ux, uy = fx / n, fy / n
    qa = ux * ux / (a * a) + uy * uy / (b * b)
    qb = -2.0 * ux * cx / (a * a)
    qc = cx * cx / (a * a) - 1.0
    t = (-qb + math.sqrt(max(0.0, qb * qb - 4 * qa * qc))) / (2 * qa)
    return t * uy


def _arc_s(y: float, b: float) -> float:
    """A monotone stand-in for arc length on one half: the ellipse parameter angle (rad) times b."""
    return b * math.asin(max(-1.0, min(1.0, y / b)))


def _pack_half(items, gap: float, b_lim: float, b: float, coord: str):
    """items: (label, target_y, lo, hi) in packing order (y increasing). Returns (ys, fits)."""
    n = len(items)
    if n == 0:
        return [], True
    d = [0.0]
    for i in range(1, n):
        prev_hi = items[i - 1][3]
        lo = items[i][2]
        d.append(d[-1] + (prev_hi - lo) + gap)
    lower = -b_lim - min(0.0, 0.0)
    upper = b_lim - d[-1]
    fits = lower <= upper
    if coord == "arc":
        # Fit on the arc parameter, then map back: rows still need the vertical distance, so the
        # result is re-spaced by a forward/backward sweep (same order, same limit).
        z = _isotonic([_arc_s(t, b) - 0.0 for _l, t, _lo, _hi in items])
        ys = [b * math.sin(max(-math.pi / 2, min(math.pi / 2, s / b))) for s in z]
        ys = _sweep(ys, d, -b_lim, b_lim)
        return ys, fits
    z = _isotonic([t - d[i] for i, (_l, t, _lo, _hi) in enumerate(items)])
    if fits:
        z = [min(max(v, lower), upper) for v in z]
    else:
        z = [lower] * n  # stack from the top limit; the bottom overflows (reported as not fitting)
    return [z[i] + d[i] for i in range(n)], fits


def _sweep(ys, d, lo_lim, hi_lim):
    n = len(ys)
    ys = list(ys)
    for i in range(n):
        ys[i] = max(ys[i], lo_lim + d[i])
        if i:
            ys[i] = max(ys[i], ys[i - 1] + (d[i] - d[i - 1]))
    for i in range(n - 1, -1, -1):
        ys[i] = min(ys[i], hi_lim - (d[-1] - d[i]))
        if i < n - 1:
            ys[i] = min(ys[i], ys[i + 1] - (d[i + 1] - d[i]))
    return ys


def _seg_cross(p, q, r, s):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2, d3, d4 = o(r, s, p), o(r, s, q), o(p, q, r), o(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and min(abs(d1), abs(d2), abs(d3), abs(d4)) > 1e-6


def _tgt_psi(item, cx, a, b):
    tx, ty = item[5], item[1]
    return math.atan2((tx - cx) / a, -ty / b)


class _Ring:
    def __init__(self, a, b, cx, is_outer, gap_y, gap_x):
        self.a, self.b, self.cx, self.is_outer = a, b, cx, is_outer
        self.gy, self.gx = gap_y, gap_x
        self._cones = {}

    def anchor(self, psi):
        return self.cx + self.a * math.sin(psi), -self.b * math.cos(psi)

    def box(self, item, psi):
        lo, hi, w = item[2], item[3], item[4]
        x, y = self.anchor(psi)
        right = x >= self.cx
        away = right if self.is_outer else not right
        x0, x1 = (x, x + w) if away else (x - w, x)
        return x0, y + lo, x1, y + hi

    def apart(self, i1, p1, i2, p2):
        a, b = self.box(i1, p1), self.box(i2, p2)
        return (b[1] - a[3] >= self.gy or a[1] - b[3] >= self.gy
                or b[0] - a[2] >= self.gx or a[0] - b[2] >= self.gx)

    def step(self, i1, p1, i2, forward=True):
        test = (lambda d: self.apart(i1, p1, i2, p1 + d)) if forward else (lambda d: self.apart(i2, p1 - d, i1, p1))
        if test(0.0):
            return 0.0
        hi = 0.002
        while not test(hi):
            hi *= 2.0
            if hi > math.pi:
                return math.pi
        lo = 0.0
        for _ in range(18):
            mid = 0.5 * (lo + hi)
            if test(mid):
                hi = mid
            else:
                lo = mid
        return hi

    def cone(self, item, tgt, margin_deg, scale):
        key = (id(item[0]), margin_deg, scale)
        if key not in self._cones:
            self._cones[key] = self._cone(item, tgt, margin_deg, scale)
        return self._cones[key]

    def _cone(self, item, tgt, margin_deg, scale):
        lab = item[0]
        ex, ey = float(lab["middle_x"]), float(lab["middle_y"])
        re = math.hypot(ex, ey)
        fa = math.atan2(ey, ex)
        if scale == math.inf:
            return tgt - math.pi * 0.9, tgt + math.pi * 0.9

        def ok(psi):
            x, y = self.anchor(psi)
            r = math.hypot(x, y)
            if r <= re:
                return False
            d = (math.atan2(y, x) - fa + math.pi) % (2 * math.pi) - math.pi
            return abs(d) <= scale * max(0.0, math.acos(min(1.0, re / r)) - math.radians(margin_deg))

        step = math.radians(0.25)
        lo = hi = tgt
        if not ok(tgt):
            return tgt, tgt
        while ok(lo - step) and lo - step > tgt - math.pi:
            lo -= step
        while ok(hi + step) and hi + step < tgt + math.pi:
            hi += step
        return lo, hi

    def sweep_forward(self, items, ps, lo_b):
        out = list(ps)
        for k in range(len(items)):
            out[k] = max(out[k], lo_b[k])
            if k:
                out[k] = max(out[k], out[k - 1] + self.step(items[k - 1], out[k - 1], items[k]))
        return out

    def sweep_backward(self, items, ps, hi_b):
        out = list(ps)
        for k in range(len(items) - 1, -1, -1):
            out[k] = min(out[k], hi_b[k])
            if k < len(items) - 1:
                out[k] = min(out[k], out[k + 1] - self.step(items[k + 1], out[k + 1], items[k], forward=False))
        return out

    def pack(self, items, targets, lo_b, hi_b, iterations=4):
        n = len(items)
        hi_b = list(hi_b)
        ps = [min(max(t, lo_b[k]), hi_b[k]) for k, t in enumerate(targets)]
        for _ in range(iterations):
            d = [0.0]
            for k in range(1, n):
                d.append(d[-1] + self.step(items[k - 1], ps[k - 1], items[k]))
            lstar, ustar, run = [], [0.0] * n, -1e9
            for k in range(n):
                run = max(run, lo_b[k] - d[k])
                lstar.append(run)
            run = 1e9
            for k in range(n - 1, -1, -1):
                run = min(run, hi_b[k] - d[k])
                ustar[k] = run
            z = _isotonic([targets[k] - d[k] for k in range(n)])
            z = [max(min(z[k], ustar[k]), lstar[k]) for k in range(n)]
            ps = [z[k] + d[k] for k in range(n)]
        for _ in range(3):
            # the cycle closes: the last label stays clear of the first one, one turn later
            if n > 1:
                hi_b[-1] = min(hi_b[-1], ps[0] + 2 * math.pi - self.step(items[0], ps[0] + 2 * math.pi, items[-1], forward=False))
            ps = self.sweep_forward(items, ps, lo_b)
            ps = self.sweep_backward(items, ps, hi_b)
        ps = self.sweep_forward(items, ps, lo_b)
        fits = all(ps[k] <= hi_b[k] + 1e-9 for k in range(n))
        if n > 1:
            fits = fits and self.apart(items[-1], ps[-1], items[0], ps[0] + 2 * math.pi)
        return ps, fits


def pack_labels(labels, total_length: int, center_x: float, x_radius: float, y_radius: float,
                is_outer: bool, spacing_px: float):
    """Place one side's labels (dicts are changed in place and returned, sorted by `middle`)."""
    from .circular import LEGACY_MIN_ROW_GAP_PX, _label_y_offsets

    o = opts()
    labels = sorted(labels, key=lambda lab: float(lab["middle"]))
    a, b, cx = float(x_radius), float(y_radius), float(center_x)
    font_px = max((float(lab.get("font_size", 8.0)) for lab in labels), default=8.0)
    if o["rows"] == "legacy":
        gaps = [LEGACY_MIN_ROW_GAP_PX]
    elif o["rows"] == "lp31":
        gaps = [max(0.22 * font_px * 96.0 / 72.0, 4.0), LEGACY_MIN_ROW_GAP_PX]
    else:
        gaps = [max(float(spacing_px), LEGACY_MIN_ROW_GAP_PX), LEGACY_MIN_ROW_GAP_PX]
    margin = float(o.get("conem", "1"))
    scales = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0] if o.get("cone", "1") == "1" else [5.0]

    items = []
    for lab in labels:
        low, _lm, high, _hm = _label_y_offsets(float(lab["middle"]), int(total_length), float(lab["height_px"]), 0.0, is_outer)
        ty = _ray_target_y(lab, cx, a, b)
        fx = float(lab["feature_middle_x"])
        tx = cx + math.copysign(a * math.sqrt(max(0.0, 1.0 - (ty / b) ** 2)), fx - cx if fx != cx else 1.0)
        items.append((lab, ty, low, high, float(lab["width_px"]), tx))
    n = len(items)
    raw = [_tgt_psi(it, cx, a, b) % (2 * math.pi) for it in items]
    # cut the cycle after the widest gap between neighbouring targets
    gapsz = [((raw[(k + 1) % n] - raw[k]) % (2 * math.pi)) for k in range(n)]
    cut = (max(range(n), key=lambda k: gapsz[k]) + 1) % n if n > 1 else 0
    order = items[cut:] + items[:cut]
    tg = [raw[(cut + k) % n] for k in range(n)]
    for k in range(1, n):
        while tg[k] < tg[k - 1]:
            tg[k] += 2 * math.pi
    def crossings(ps):
        return crossings_for(order, ps)

    def crossings_for(seq, ps):
        """Crossing leader pairs of the drawn leaders (contact points from today's leader assignment)."""
        from .circular import _assign_leader_start_points

        ring0 = _Ring(a, b, cx, is_outer, 0.1, 1.0)
        trial = []
        for it, psi in zip(seq, ps):
            lab = dict(it[0])
            lab["start_x"], lab["start_y"] = ring0.anchor(psi)
            trial.append(lab)
        trial = _assign_leader_start_points(sorted(trial, key=lambda lab: float(lab["middle"])), total_length)
        segs = []
        for lab in trial:
            e = (float(lab["middle_x"]), float(lab["middle_y"]))
            q = (float(lab.get("leader_start_x", lab["start_x"])), float(lab.get("leader_start_y", lab["start_y"])))
            segs.append((e, q))
        cnt = 0
        for i in range(len(segs)):
            p1, q1 = segs[i]
            for j in range(i + 1, len(segs)):
                if _seg_cross(p1, q1, *segs[j]):
                    cnt += 1
        return cnt

    best = None
    tried = []
    eps = math.asin(min(1.0, 1.0 / a)) if a > 1 else 0.0
    if o.get("half", "1") == "1":
        # Candidate family "half": each label stays on one side of the vertical axis (today's rule: genome
        # halves); labels move across a pole only as whole end runs (two split points, greedy).
        half = next((k for k, it in enumerate(items) if float(it[0]["middle"]) >= total_length / 2), n)
        raw_all = [_tgt_psi(it, cx, a, b) for it in items]

        def build(st, sb):
            right = (list(range(n + st, n)) if st < 0 else []) + list(range(max(st, 0), half - max(sb, 0))) + (
                list(range(half, half - sb)) if sb < 0 else [])
            left = (list(range(half - sb, half)) if sb > 0 else []) + list(range(half + max(-sb, 0), n - max(-st, 0))) + (
                list(range(0, st)) if st > 0 else [])
            return right, left

        def side_seq(idx, lo_lim, hi_lim):
            seq = [items[k] for k in idx]
            tgts = []
            for k in idx:
                t = raw_all[k] % (2 * math.pi)
                mid = 0.5 * (lo_lim + hi_lim)
                t = t + 2 * math.pi * round((mid - t) / (2 * math.pi))
                tgts.append(min(max(t, lo_lim), hi_lim))
            for k in range(1, len(tgts)):
                tgts[k] = max(tgts[k], tgts[k - 1])
            return seq, tgts

        sides = ((eps, math.pi - eps), (math.pi + eps, 2 * math.pi - eps))

        def excess(ring, st, sb):
            worst = -1.0
            for idx, (lo_l, hi_l) in zip(build(st, sb), sides):
                if not idx:
                    continue
                seq, _t = side_seq(idx, lo_l, hi_l)
                ps = ring.sweep_forward(seq, [lo_l] * len(seq), [lo_l] * len(seq))
                worst = max(worst, ps[-1] - hi_l)
            return worst

        for gap in gaps:
            ring = _Ring(a, b, cx, is_outer, gap, 1.0)
            state, cur = (0, 0), excess(ring, 0, 0)
            while cur > 0:
                moves = [(state[0] + ds, state[1] + db) for ds, db in ((1, 0), (-1, 0), (0, 1), (0, -1))]
                moves = [mv for mv in moves if -(n - half) < mv[0] < half and -(n - half) < mv[1] < half]
                cand = min(((excess(ring, *mv), mv) for mv in moves), default=None)
                if cand is None or cand[0] >= cur:
                    break
                cur, state = cand
            if cur > 0:
                tried.append(("half", gap, "nofit", state))
                continue
            K = int(o.get("splitk", "20"))
            states = []
            for st in range(-K, K + 1):
                for sb in sorted(range(-K, K + 1), key=abs):
                    if -(n - half) < st < half and -(n - half) < sb < half and excess(ring, st, sb) <= 0:
                        states.append((st, sb))
                        break
            states.sort(key=lambda x: abs(x[0]) + abs(x[1]))
            for st2 in states:
                ps_all, seq_all, tg_all, ok = [], [], [], True
                for idx, (lo_l, hi_l) in zip(build(*st2), sides):
                    if not idx:
                        continue
                    seq, tg_s = side_seq(idx, lo_l, hi_l)
                    ps, fits = ring.pack(seq, tg_s, [lo_l] * len(seq), [hi_l] * len(seq))
                    ok = ok and fits and ps[-1] <= hi_l + 1e-9
                    ps_all += ps
                    seq_all += seq
                    tg_all += tg_s
                if not ok:
                    continue
                cr = crossings_for(seq_all, ps_all)
                disp = sum((p - t) ** 2 for p, t in zip(ps_all, tg_all))
                tried.append(("half", gap, cr, st2))
                key = (cr, disp)
                if best is None or key < best[0]:
                    best = (key, ("half", gap, ps_all, True, seq_all))
            if best is not None:
                break
            ps_all, seq_all, tg_all, ok = [], [], [], True
            for idx, (lo_l, hi_l) in zip(build(*state), sides):
                if not idx:
                    continue
                seq, tg_s = side_seq(idx, lo_l, hi_l)
                ps, fits = ring.pack(seq, tg_s, [lo_l] * len(seq), [hi_l] * len(seq))
                ok = ok and fits and ps[-1] <= hi_l + 1e-9
                ps_all += ps
                seq_all += seq
                tg_all += tg_s
            if not ok:
                tried.append(("half", gap, "nofit", state))
                continue
            cr = crossings_for(seq_all, ps_all)
            disp = sum((p - t) ** 2 for p, t in zip(ps_all, tg_all))
            tried.append(("half", gap, cr, state))
            key = (cr, disp)
            if best is None or key < best[0]:
                best = (key, ("half", gap, ps_all, True, seq_all))
            break
    for scale in (scales if not (best and best[0][0] == 0) else []):
        for gap in gaps:
            ring = _Ring(a, b, cx, is_outer, gap, 1.0)
            cones = [ring.cone(it, t, margin, scale) for it, t in zip(order, tg)]
            ps, fits = ring.pack(order, tg, [c[0] for c in cones], [c[1] for c in cones])
            if not fits:
                tried.append((scale, gap, "nofit"))
                continue
            cr = crossings(ps)
            disp = sum((p - t) ** 2 for p, t in zip(ps, tg))
            tried.append((scale, gap, cr))
            key = (cr, disp)
            if best is None or key < best[0]:
                best = (key, (scale, gap, ps, True, order))
            if cr == 0:
                break
        if best is not None and best[0][0] == 0:
            break
    if best is None:
        ring = _Ring(a, b, cx, is_outer, gaps[-1], 1.0)
        ps = ring.sweep_forward(order, tg, [-1e9] * n)
        best = ((None, None), (None, gaps[-1], ps, False, order))
    LAST.update(tried=tried, crossings_est=best[0][0])
    result = best[1]
    scale, gap, ps, fits, order = result
    for it, psi in zip(order, ps):
        x, y = _Ring(a, b, cx, is_outer, gap, 1.0).anchor(psi)
        it[0]["start_x"], it[0]["start_y"] = x, y
        it[0]["_clp_packed"] = True
    LAST.update(packed=True, fits=fits, cone_scale=scale, gaps=[gap], is_outer=is_outer, n=n)
    return labels
