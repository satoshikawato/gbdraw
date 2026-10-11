"""Study prototype E: cluster callout lists (plasmid-map convention). Never merged.

Rule in short (details in the study report logs/design-E.md):
  1. Every label has a home: its own feature's angle, just outside the ring.
  2. Labels whose home boxes overlap (joined transitively) form a run.  A run longer than `max_rows` is cut at its
     widest angular gaps until every piece has at most `max_rows` rows.  A piece (or a lone label) is a block.
  3. A block is one flush list of stacked lines, on the right or left of the ring (side of its mean angle),
     left-aligned on the right and right-aligned on the left.  Blocks of one side are stacked top to bottom in
     angular order by an optimal 1-D packing (they slide along the ring, never outward), a block-gap apart.
  4. Leaders: `fan` = one elbowed leader per row; `bracket` = one leader to a bar along the list plus a bracket
     (ticks and an arc) over the run's features; `auto` picks bracket for dense lists.
Options: leader=fan|bracket|auto, max_rows, dist, hug (0/1: rows follow the ring instead of one flush column),
bracket_min, sep_min, block_gap, row_pad, maxh (drop labels when one side's stack exceeds maxh * ring diameter).
"""

from __future__ import annotations

import math

from .common import (
    Item, Ctx, angle_diff, angle_of, box_from_left, box_from_right, boxes_overlap, box_min_radius, emit, pack_ordered,
    polar, segments_cross, segment_hits_box, polyline_segments,
)

TWO_PI = 2.0 * math.pi


class _Row:
    __slots__ = ("item", "theta", "elbow", "home_y", "home_box")


class _Block:
    def __init__(self, rows: list[_Row], ctx: Ctx, pitch: float, r_b: float):
        self.rows = rows            # clockwise order (ascending phi)
        self.n = len(rows)
        sx = sum(math.sin(r.theta) for r in rows)
        sy = sum(-math.cos(r.theta) for r in rows)
        self.mean = angle_of(sx, sy)
        self.right = math.sin(self.mean) >= 0.0
        self.key = self.mean if self.right else TWO_PI - self.mean
        self.pitch = pitch
        self.target = sum(r.home_y for r in rows) / self.n
        self.centre = self.target
        self.order = list(rows)     # top to bottom
        self.bracket = False

    @property
    def height(self) -> float:
        return self.n * self.pitch


def _home(item: Item, r_b: float, mid_gap: float, gap: float):
    row = _Row()
    row.item = item
    row.theta = angle_of(*item.anchor)
    ux, uy = math.sin(row.theta), -math.cos(row.theta)
    row.home_y = uy * (r_b + 0.5 * item.h)
    x = ux * r_b
    if ux >= 0:
        row.home_box = box_from_left(max(x, mid_gap), row.home_y, item.w, item.h)
    else:
        row.home_box = box_from_right(min(x, -mid_gap), row.home_y, item.w, item.h)
    return row


def _runs(rows: list[_Row], gap: float, link: float) -> list[list[_Row]]:
    n = len(rows)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if abs(angle_diff(rows[i].theta, rows[j].theta)) > math.radians(35):
                continue
            if boxes_overlap(rows[i].home_box, rows[j].home_box, gap + link):
                parent[find(i)] = find(j)
    groups: dict[int, list[_Row]] = {}
    for i, row in enumerate(rows):
        groups.setdefault(find(i), []).append(row)
    out = []
    for group in groups.values():
        # start the clockwise order after the widest cyclic gap
        group.sort(key=lambda r: r.theta)
        if len(group) > 1:
            gaps = [(group[(k + 1) % len(group)].theta - group[k].theta) % TWO_PI for k in range(len(group))]
            k = max(range(len(group)), key=lambda q: gaps[q])
            group = group[k + 1:] + group[:k + 1]
        out.append(group)
    return out


def _phi(group: list[_Row]) -> list[float]:
    base = group[0].theta
    return [(r.theta - base) % TWO_PI for r in group]


def _split(group: list[_Row], max_rows: int) -> list[list[_Row]]:
    if len(group) <= max_rows:
        return [group]
    phi = _phi(group)
    n = len(group)
    lo, hi = max(1, math.ceil(n * 0.3)), min(n - 1, math.floor(n * 0.7))
    k = max(range(lo, hi + 1), key=lambda q: phi[q] - phi[q - 1])
    return _split(group[:k], max_rows) + _split(group[k:], max_rows)


def _x_edge(right: bool, r_b: float, yc: float, h: float, mid_gap: float) -> float:
    """Smallest |x| at which a box of height h centred on yc clears the circle of radius r_b."""
    lo, hi = yc - 0.5 * h, yc + 0.5 * h
    dy = 0.0 if lo < 0.0 < hi else min(abs(lo), abs(hi))
    return max(mid_gap, math.sqrt(max(r_b * r_b - dy * dy, 0.0)))


def _rows_geometry(block: _Block, order: list[_Row], r_b: float, mid_gap: float, hug: bool):
    """Boxes (top to bottom) of a block placed with its centre at block.centre."""
    n = block.n
    y0 = block.centre - 0.5 * n * block.pitch
    ys = [y0 + (k + 0.5) * block.pitch for k in range(n)]
    xs = [_x_edge(block.right, r_b, y, order[k].item.h, mid_gap) for k, y in enumerate(ys)]
    if not hug:
        xs = [max(xs)] * n
    # the list starts beyond every elbow of its own rows, so its leaders arrive from the ring side, not through text
    reach = max((r.elbow[0] if block.right else -r.elbow[0]) for r in order) + 6.0
    xs = [max(x, reach) for x in xs]
    boxes = []
    for k, row in enumerate(order):
        w, h = row.item.w, row.item.h
        boxes.append(box_from_left(xs[k], ys[k], w, h) if block.right else box_from_right(-xs[k], ys[k], w, h))
    return ys, boxes


def _fan_leaders(block, order, boxes):
    leaders = []
    for row, box in zip(order, boxes):
        end = (box[0] - 2.0, 0.5 * (box[1] + box[3])) if block.right else (box[2] + 2.0, 0.5 * (box[1] + box[3]))
        leaders.append([row.item.anchor, row.elbow, end])
    return leaders


def _count_crossings(leaders) -> int:
    segs = [polyline_segments(l) for l in leaders]
    total = 0
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            if any(segments_cross(p, q, r, s) for p, q in segs[i] for r, s in segs[j]):
                total += 1
    return total


def _use_bracket(block: _Block, mode: str, r_elbow: float, bracket_min: int, sep_min: float) -> bool:
    if mode == "fan" or block.n < 2:
        return False
    if mode == "bracket":
        return block.n >= 2
    if block.n < bracket_min:
        return False
    thetas = sorted(r.theta for r in block.rows)
    spacing = [r_elbow * abs(angle_diff(b, a)) for a, b in zip(thetas, thetas[1:])]
    close = sum(1 for s in spacing if s < sep_min)
    return close >= 0.5 * len(spacing)


def _bracket_geometry(block: _Block, order, boxes, r_elbow: float):
    r_br = r_elbow + 3.0
    thetas = [r.theta for r in block.rows]
    ref = block.mean
    rel = [angle_diff(t, ref) for t in thetas]
    lo, hi = min(rel), max(rel)
    steps = max(2, int(math.ceil((hi - lo) / math.radians(1.5))))
    arc = [polar(r_br, ref + lo + (hi - lo) * k / steps) for k in range(steps + 1)]
    ticks = [[r.item.anchor, polar(r_br, ref + a)] for r, a in zip(block.rows, rel)]
    mid = min(range(block.n), key=lambda k: abs(rel[k]))
    x_bar = (boxes[0][0] - 5.0) if block.right else (boxes[0][2] + 5.0)
    y_top, y_bot = boxes[0][1] + 1.0, boxes[-1][3] - 1.0
    bar = [(x_bar, y_top), (x_bar, y_bot)]
    meet = (x_bar, 0.5 * (y_top + y_bot))
    start = polar(r_br, ref + rel[mid])
    leader = [block.rows[mid].item.anchor, start, meet]
    return mid, leader, [arc, bar, *ticks]


def _choose_stack(blocks: list[_Block], bg: float, limit: float = 0.0):
    blocks.sort(key=lambda b: b.key)
    sizes = [b.height + (bg if b.n >= 2 else 0.0) for b in blocks]
    if limit > 0:
        centres = pack_ordered([b.target for b in blocks], sizes, -limit, limit)
    else:
        centres = pack_ordered([b.target for b in blocks], sizes)
    for b, c in zip(blocks, centres):
        b.centre = c
    return centres, sizes


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    mode = ctx.options.get("leader", "auto")
    max_rows = int(ctx.opt("max_rows", 8))
    dist = ctx.opt("dist", 8.0)
    hug = bool(ctx.opt("hug", 0))
    bracket_min = int(ctx.opt("bracket_min", 4))
    sep_min = ctx.opt("sep_min", 8.0)
    bg = ctx.opt("block_gap", 0.4) * (max((i.h for i in items), default=0.0) + ctx.gap)
    row_pad = ctx.opt("row_pad", 0.0)
    mid_gap = ctx.opt("mid_gap", 10.0)
    link = ctx.opt("link", 0.0)
    pitch = max(i.h for i in items) + ctx.gap + row_pad
    r_b = ctx.r_keep + dist

    rows = [_home(it, r_b, mid_gap, ctx.gap) for it in items]
    for row in rows:
        row.elbow = polar(ctx.r_elbow, row.theta)

    def build_blocks(active: list[_Row]) -> list[_Block]:
        blocks = []
        for run in _runs(active, ctx.gap, link):
            for piece in _split(run, max_rows):
                blocks.append(_Block(piece, ctx, pitch, r_b))
        return blocks

    # capacity rule (fit=1): the lists of one side may not be taller than the ring's height plus `over` rows at
    # each end.  When a side is too tall, blocks nearest 12 or 6 o'clock move to the other side while it has room
    # (flip=1); if both are full the label with the shortest feature is dropped from the longest list.
    # fit=0: no limit, the stack grows and the canvas with it.
    limit = (r_b + ctx.opt("over", 12.0) * pitch) if ctx.opt("fit", 1) else 0.0
    active = list(rows)
    dropped_rows: list[_Row] = []
    flips = 0

    def side_height(blocks, right):
        return sum(b.height + (bg if b.n >= 2 else 0.0) for b in blocks if b.right == right)

    while True:
        blocks = build_blocks(active)
        while limit > 0 and ctx.opt("flip", 0):
            full = [rt for rt in (True, False) if side_height(blocks, rt) > 2.0 * limit]
            if not full:
                break
            rt = full[0]
            movable = [b for b in blocks if b.right == rt and abs(math.sin(b.mean)) < ctx.opt("flip_sin", 0.5)
                       and side_height(blocks, not rt) + b.height + bg <= 2.0 * limit]
            if not movable:
                break
            b = min(movable, key=lambda q: abs(math.sin(q.mean)))
            b.right = not rt
            b.key = b.mean if b.right else TWO_PI - b.mean
            flips += 1
        victim = None
        for right in (True, False):
            side = [b for b in blocks if b.right == right]
            if limit <= 0 or not side:
                continue
            _choose_stack(side, bg)            # unclamped least-squares stack (blocks slide only as far as needed)
            top = min(b.centre - 0.5 * b.height for b in side)
            bottom = max(b.centre + 0.5 * b.height for b in side)
            if top < -limit or bottom > limit:
                end = side[:5] if -limit - top > bottom - limit else side[-5:]
                worst = max(end, key=lambda b: b.n)
                victim = min(worst.rows, key=lambda r: r.item.span_bp)
                break
        if victim is None or not active:
            break
        dropped_rows.append(victim)
        active.remove(victim)
    for row in dropped_rows:
        ctx.dropped.append(row.item.text)

    for right in (True, False):
        side = [b for b in blocks if b.right == right]
        if side:
            _choose_stack(side, bg)

    # order rows top to bottom: clockwise on the right, counter-clockwise on the left; try both, keep fewer crossings
    for b in blocks:
        default = list(b.rows) if b.right else list(reversed(b.rows))
        alt = list(reversed(default))
        best = None
        for cand in (default, alt) if b.n > 1 else (default,):
            _, boxes = _rows_geometry(b, cand, r_b, mid_gap, hug)
            cross = _count_crossings(_fan_leaders(b, cand, boxes))
            length = sum(math.dist(l[1], l[2]) for l in _fan_leaders(b, cand, boxes))
            score = (cross, length)
            if best is None or score < best[0]:
                best = (score, cand)
        b.order = best[1]
        b.bracket = _use_bracket(b, mode, ctx.r_elbow, bracket_min, sep_min)

    out: list[dict] = []
    slides = []
    all_boxes: list[tuple] = []
    all_leaders: list[list] = []
    owner: list[int] = []
    for bi, b in enumerate(blocks):
        ys, boxes = _rows_geometry(b, b.order, r_b, mid_gap, hug)
        align = "start" if b.right else "end"
        if b.bracket:
            mid, leader, extra = _bracket_geometry(b, b.order, boxes, ctx.r_elbow)
            mid_row = b.rows[mid]
            for row, box in zip(b.order, boxes):
                if row is mid_row:
                    out.append(emit(row.item, box, leader, align=align, extra=extra, tag="bracket"))
                    all_leaders.append(leader)
                    owner.append(bi)
                    for line in extra:
                        all_leaders.append(line)
                        owner.append(bi)
                else:
                    out.append(emit(row.item, box, None, align=align, tag="bracket-row"))
                all_boxes.append(box)
        else:
            leaders = _fan_leaders(b, b.order, boxes)
            for row, box, leader in zip(b.order, boxes, leaders):
                out.append(emit(row.item, box, leader, align=align, tag="fan" if b.n > 1 else "single"))
                all_boxes.append(box)
                all_leaders.append(leader)
                owner.append(bi)
        for row, box in zip(b.order, boxes):
            slides.append(math.hypot(0.5 * (box[0] + box[2]) - 0.5 * (row.home_box[0] + row.home_box[2]),
                                     0.5 * (box[1] + box[3]) - row.home_y))

    # self-check including bracket strokes (the shared metrics only see the per-row leaders)
    overlaps = sum(1 for i in range(len(all_boxes)) for j in range(i + 1, len(all_boxes))
                   if boxes_overlap(all_boxes[i], all_boxes[j], 0.0))
    inside = sum(1 for bx in all_boxes if box_min_radius(bx) < ctx.r_keep - 1.0)
    cross = through = 0
    segs = [polyline_segments(l) for l in all_leaders]
    for i in range(len(segs)):
        for bx in all_boxes:
            if any(segment_hits_box(p, q, bx, -0.5) for p, q in segs[i]):
                through += 1
        for j in range(i + 1, len(segs)):
            if owner[i] == owner[j]:
                continue
            if any(segments_cross(p, q, r, s) for p, q in segs[i] for r, s in segs[j]):
                cross += 1
    sizes: dict[int, int] = {}
    for b in blocks:
        sizes[b.n] = sizes.get(b.n, 0) + 1
    ctx.notes.update(
        blocks=len(blocks), block_sizes=dict(sorted(sizes.items())), lists=sum(1 for b in blocks if b.n >= 2),
        brackets=sum(1 for b in blocks if b.bracket), slide_mean=round(sum(slides) / max(len(slides), 1), 1),
        slide_max=round(max(slides, default=0.0), 1), self_overlaps=overlaps, boxes_inside_keep=inside,
        strokes_through_text=through, stroke_crossings_incl_extra=cross, leader=mode, max_rows=max_rows, flips=flips,
        side_rows=[sum(b.n for b in blocks if b.right), sum(b.n for b in blocks if not b.right)],
    )
    return out
