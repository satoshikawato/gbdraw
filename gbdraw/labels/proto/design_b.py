"""Design B: staggered tiers (study prototype).

Each half of the map (right: theta in [0, pi], left mirrored onto the right) is a set of vertical stacks ("tiers").
Tier 0 follows today's label ellipse; further tiers are flush columns beyond the widest label of the tier inside.
A label of tier t keeps its order (theta), so leaders never cross: they run feature -> elbow ring -> a thin
horizontal slot between two rows of every inner tier -> a bus corridor -> the label.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .common import (
    Item, Ctx, angle_of, emit, pack_ordered, polar, segment_min_radius, box_min_radius,
)


@dataclass
class Rec:
    item: Item
    th: float            # mirrored angle in [0, pi]
    ex: float            # elbow point (mirrored)
    ey: float
    anchor: tuple
    w: float
    h: float
    y_nat: float
    tier: int = 0
    y: float = 0.0       # packed centre of the label
    slot: float = 0.0    # y of its slot in the stacks of the inner tiers (per inner tier)
    x0: float = 0.0
    slots: dict = None


def _dome(a0, b0, y):
    t = 1.0 - (y / b0) ** 2
    return a0 * math.sqrt(t) if t > 0 else 0.0


def _dome_radii(ctx):
    floor = ctx.r_keep + 12.0
    return max(ctx.arc_x_radius, floor), max(ctx.arc_y_radius, floor)


def _start0(a0, b0, xfloor, y, r):
    x = max(_dome(a0, b0, y), xfloor)
    if abs(y) > abs(r.ey) and (y * r.ey > 0):
        x = max(x, r.ex + 2.0)
    return x


def _stack(recs_labels, recs_slots, b0, slot_size, gap):
    """Pack labels of one tier plus slots of the outer tiers. Returns (y of label, y of slot) by id."""
    entries = [(r.th, 0, id(r), r) for r in recs_labels] + [(r.th, 1, id(r), r) for r in recs_slots]
    entries.sort(key=lambda e: (e[0], e[1]))
    targets = [e[3].y_nat for e in entries]
    sizes = [(e[3].h + gap) if e[1] == 0 else slot_size for e in entries]
    centres = pack_ordered(targets, sizes)
    ly, sy = {}, {}
    for e, c in zip(entries, centres):
        (ly if e[1] == 0 else sy)[e[2]] = c
    return ly, sy


def _fits(recs_labels, recs_slots, b0, slot_size, gap, drift):
    ly, _ = _stack(recs_labels, recs_slots, b0, slot_size, gap)
    return all(abs(ly[id(r)] - r.y_nat) <= drift for r in recs_labels)


def _caps(recs, ntiers, wq):
    """Width caps per tier (inner tiers take the short labels; the outermost tier takes anything)."""
    if not wq or ntiers < 2:
        return [1e9] * ntiers
    ws = sorted(r.w for r in recs)
    caps = [ws[min(len(ws) - 1, int(len(ws) * (t + 1) / ntiers * wq))] for t in range(ntiers - 1)]
    return caps + [1e9]


def _assign(recs, ntiers, b0, slot_size, gap, drift, wq=0.0):
    """Sweep in theta order: a label stays in the current tier while the tier's stack (with slots for the labels
    already sent outwards) keeps every label within `drift` of its natural place; otherwise it goes outwards.
    Returns the tiers, or None when the last tier does not fit either."""
    remaining = list(recs)
    tiers = []
    caps = _caps(recs, ntiers, wq)
    for t in range(ntiers):
        if not remaining:
            tiers.append([])
            continue
        if t == ntiers - 1:
            if not _fits(remaining, [], b0, slot_size, gap, drift):
                return None
            tiers.append(remaining)
            remaining = []
            continue
        selected, rejected = [], []
        for c in remaining:
            if c.w <= caps[t] and _fits(selected + [c], rejected, b0, slot_size, gap, drift):
                selected.append(c)
            else:
                rejected.append(c)
        tiers.append(selected)
        remaining = rejected
    return tiers


def _choose_tiers(recs, ctx, b0, slot_size, gap):
    """Fewest tiers (at most `tiers`) in which every label stays within `drift` px of its natural place; when
    even the most tiers do not fit, the bound is relaxed by x1.5 steps (the bound used is kept in the notes)."""
    dcap = ctx.opt("drift", 240.0)
    tmax = int(ctx.opt("tiers", 4))
    for step in range(6):
        d = dcap * 1.5 ** step
        for t in range(1, tmax + 1):
            tiers = _assign(recs, t, b0, slot_size, gap, d, ctx.opt("wq", 1.0))
            if tiers is not None:
                ctx.notes.setdefault("drift_used", []).append(round(d))
                return tiers
    raise RuntimeError("tiers: no fit")


def _route(points, rmin):
    """Insert waypoints on the circle rmin where a segment would cut into it."""
    out = [points[0]]
    for p, q in zip(points, points[1:]):
        if segment_min_radius(p, q) < rmin - 1.0:
            a0, a1 = angle_of(*p), angle_of(*q)
            d = (a1 - a0 + math.pi) % (2 * math.pi) - math.pi
            steps = max(2, int(abs(d) / math.radians(12)) + 1)
            for k in range(1, steps):
                out.append(polar(rmin, a0 + d * k / steps))
        out.append(q)
    return out


def _solve_side(recs, ctx: Ctx, side: str):
    a0, b0 = _dome_radii(ctx)
    gap = ctx.gap
    slot_size = ctx.opt("slot", 4.0)
    cw = ctx.opt("corridor", 36.0)
    xfloor = ctx.opt("xfloor", 8.0)
    if not recs:
        return {}
    recs.sort(key=lambda r: r.th)
    tiers = _choose_tiers(recs, ctx, b0, slot_size, gap)
    for t, members in enumerate(tiers):
        for r in members:
            r.tier = t
    # pack every tier with the slots of the outer tiers
    for t, members in enumerate(tiers):
        outer = [r for tt in tiers[t + 1:] for r in tt]
        ly, sy = _stack(members, outer, b0, slot_size, gap)
        for r in members:
            r.y = ly[id(r)]
        for r in outer:
            r.slots[t] = sy[id(r)]
    # x geometry
    xcol = {}
    xout = {}
    for t, members in enumerate(tiers):
        if not members:
            continue
        if t == 0:
            for r in members:
                r.x0 = _start0(a0, b0, xfloor, r.y, r)
        else:
            for r in members:
                r.x0 = xcol[t]
        # the band must keep r_keep
        for r in members:
            box = (r.x0, r.y - 0.5 * r.h, r.x0 + r.w, r.y + 0.5 * r.h)
            while box_min_radius(box) < ctx.r_keep and r.x0 < 4000:
                r.x0 += 4.0
                box = (r.x0, box[1], r.x0 + r.w, box[3])
        xout[t] = max(r.x0 + r.w for r in members) + 8.0
        xcol[t + 1] = xout[t] + cw
    # leaders
    for t, members in enumerate(tiers):
        for r in members:
            pts = [r.anchor, (r.ex, r.ey)]
            for s in range(t):
                ys = r.slots[s]
                if s == 0:
                    xin = _start0(a0, b0, xfloor, ys, r) - 2.0
                else:
                    xin = xcol[s] - 2.0
                pts.append((xin, ys))
                pts.append((xout[s], ys))
            pts.append((r.x0 - 2.0, r.y))
            r.leader = _route(pts, max(i.feat_outer_r for i in ctx.items) + ctx.opt('routepad', 6.0))
    return {id(r): r for members in tiers for r in members}, tiers, xcol


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    ctx.items = items
    ctx.notes["ctx"] = [round(ctx.arc_x_radius), round(ctx.arc_y_radius), round(ctx.r_elbow), round(ctx.r_keep)]
    recs_by_side = {"R": [], "L": []}
    for it in items:
        right = it.theta <= math.pi
        phi = it.theta if right else 2.0 * math.pi - it.theta
        sgn = 1.0 if right else -1.0
        ex, ey = polar(ctx.r_elbow, it.theta)
        rec = Rec(item=it, th=phi, ex=sgn * ex, ey=ey, anchor=(sgn * it.anchor[0], it.anchor[1]),
                  w=it.w, h=it.h, y_nat=-_dome_radii(ctx)[1] * math.cos(phi))
        rec.slots = {}
        rec.sgn = sgn
        recs_by_side["R" if right else "L"].append(rec)
    out = []
    notes = {}
    for side, recs in recs_by_side.items():
        res = _solve_side(recs, ctx, side)
        if not res:
            continue
        _, tiers, xcol = res
        notes[side] = [len(m) for m in tiers]
        for r in recs:
            s = r.sgn
            box = (r.x0, r.y - 0.5 * r.h, r.x0 + r.w, r.y + 0.5 * r.h)
            leader = [(s * x, y) for x, y in r.leader]
            if s < 0:
                box = (-box[2], box[1], -box[0], box[3])
            out.append(emit(r.item, box, leader, align="start" if s > 0 else "end", tag=f"t{r.tier}"))
    ctx.notes["tiers"] = notes
    return out
