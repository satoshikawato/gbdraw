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
Options: clear 3, inside 1, cone 70, lgap 2.5, lpad 2.4, stub 1, edge 1, fl 12, cutmax 30, cutn 14, hide 0, tw 0.8, ray 1, plat 2, sym 1 (see _Par).
"""

from __future__ import annotations

import dataclasses
import heapq
import math
from pathlib import Path

import numpy as np
from fontTools.ttLib import TTFont
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import distance_transform_edt

from ...auto_sizes import circular_tick_intervals
from ...core.numeric import scaled_tick_text
from ...core.text import calculate_bbox_dimensions
from .common import Ctx, Item, emit, pack_ordered
from .design_g import _seg_box_hit, _seg_cross_matrix

EDGE_FR = (0.5, 0.2, 0.8)
OFFS = np.array([0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 30, 36, 44, 52, 62, 75, 90, 110, 135], dtype=float)
OFFS_F = np.concatenate([np.arange(0.0, 41.0), [44, 52, 62, 75, 90, 110, 135]])      # face=1: box positions 1 px apart near the stub (a fan can be tight to a few tenths of a px)
LEAN = np.tan(np.radians([2, 5, 8, 11, 14, 17, 20, 24, 28, 33, 39, 46, 54, 62]))
SIDE_FR = (0.08, 0.2, 0.35, 0.5, 0.65, 0.8, 0.92)
NEAR = [0.0, 1.0, 2.0, 4.0, 6.0, 9.0]
PUSH = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
while PUSH[-1] < 4000.0:
    PUSH.append(round(PUSH[-1] * 1.12 + 1.0, 1))

# Prototype P (LP-D22 to LP-D25): the 12 contact points of the ink box, in the label frame (fraction along the text from the end that
# faces the leader, fraction of the height from the top): 4 vertices, the middle of each short side, middle and quarter points of each long side.
PORTS = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.5), (1.0, 0.5), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0)]
PORT_NAMES = {(0.0, 0.0): "TL", (1.0, 0.0): "TR", (0.0, 1.0): "BL", (1.0, 1.0): "BR", (0.0, 0.5): "L", (1.0, 0.5): "R",
              (0.25, 0.0): "T1", (0.5, 0.0): "T2", (0.75, 0.0): "T3", (0.25, 1.0): "B1", (0.5, 1.0): "B2", (0.75, 1.0): "B3"}


def _plen(pts) -> float:
    return sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))


def _sturn(pts, min_seg: float = 3.0) -> float:
    """Signed total turn in degrees (positive: clockwise in image coordinates), pieces under 3 px dropped."""
    pts = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= min_seg]
    turn = 0.0
    for a, b, c in zip(pts, pts[1:], pts[2:]):
        u, v = (b[0] - a[0], b[1] - a[1]), (c[0] - b[0], c[1] - b[1])
        nu, nv = math.hypot(*u), math.hypot(*v)
        if nu >= min_seg and nv >= min_seg:
            ang = math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv)))))
            turn += ang if u[0] * v[1] - u[1] * v[0] > 0 else -ang
    return turn


def _select(nd: dict, idx) -> dict:
    """The seats `idx` of a label's seat table."""
    return {k: (v[idx] if isinstance(v, np.ndarray) and k != "stub" else v) for k, v in nd.items()}


def _off_arc(radii: list, step: float = 2.0) -> bool:
    """True when the leader ends of a cluster, in stack order, change their distance from the centre by more than `step` px (scripts/arc_fit.py)."""
    return any(abs(b - a) > step for a, b in zip(radii, radii[1:]))


PORTS_ARR = np.array(PORTS, dtype=float)


def _pav_l1(t) -> np.ndarray:
    """Nondecreasing z minimising sum |z - t| (pool adjacent violators with medians).  Every pooled run takes the value of one of its two
    middle targets (the one nearer the mean), so that its middle label sits exactly on its target."""
    blocks: list = []
    for v in t:
        blocks.append([float(v)])
        while len(blocks) > 1 and blocks[-2][(len(blocks[-2]) - 1) // 2] > blocks[-1][(len(blocks[-1]) - 1) // 2]:
            b = blocks.pop()
            blocks[-1] = sorted(blocks[-1] + b)
    out: list = []
    prev = -math.inf
    for b in blocks:
        lo, hi = b[(len(b) - 1) // 2], b[len(b) // 2]
        mean = sum(b) / len(b)
        val = lo if abs(lo - mean) <= abs(hi - mean) else hi
        val = max(val, prev)
        out += [val] * len(b)
        prev = val
    return np.array(out)


def _lipschitz(r, slope: float) -> np.ndarray:
    """The smallest sequence at or above r whose neighbours differ by at most `slope`."""
    f = np.array(r, dtype=float)
    for k in range(1, len(f)):
        f[k] = max(f[k], f[k - 1] - slope)
    for k in range(len(f) - 2, -1, -1):
        f[k] = max(f[k], f[k + 1] - slope)
    return f


def _broken(flat) -> bool:
    """True when a False (a bent leader) lies between two True (straight) ones."""
    if True not in flat:
        return False
    return not all(flat[flat.index(True):len(flat) - flat[::-1].index(True)])


def _sturn_arr(anchor, stub, att) -> np.ndarray:
    """Signed turn (degrees, positive clockwise in image coordinates) of the leaders anchor -> stub -> each row of att."""
    u = (stub[0] - anchor[0], stub[1] - anchor[1])
    vx, vy = att[:, 0] - stub[0], att[:, 1] - stub[1]
    nv = np.hypot(vx, vy)
    cosv = np.clip((u[0] * vx + u[1] * vy) / np.maximum(math.hypot(*u) * nv, 1e-9), -1.0, 1.0)
    ang = np.degrees(np.arccos(cosv))
    return np.where(nv < 3.0, 0.0, np.where(u[0] * vy - u[1] * vx > 0, ang, -ang))


SYM_K = 0.5                                              # px per degree^2 of deviation from the wanted lean of a label of a symmetric fan


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


FONT_FILE = Path(__file__).resolve().parents[2] / "data" / "LiberationSans-Regular.ttf"
INK_K, INK_PAD, INK_LEFT, BASE_EM = 4, 18.0, 0.2, 0.352     # field resolution (px^-1), field margin around the ink, first ink column and baseline (em below the middle) as the cairosvg raster draws them
_INK: dict = {}


class _Ink:
    """Distance (px) to the ink of one text, on a grid of 1/INK_K px, in the frame of the label box: q along the text from the box
    end that faces the leader (the first ink column stands INK_LEFT px in), y from the middle of the box.  Mirrored for the labels
    whose text ends at the leader (text-anchor end).  Glyph outlines come from the font file (no kerning, as cairo draws them)."""

    def __init__(self, text: str, fs: float, mirror: bool):
        font = ImageFont.truetype(str(FONT_FILE), int(round(fs * 2 * INK_K)))
        cmap, hmtx, upm = _font_tables()
        scale = 2 * INK_K * fs / upm
        wid = int(sum(hmtx[cmap[ord(c)]][0] for c in text) * scale) + 4 * INK_K * 2
        base = int(1.2 * fs * 2 * INK_K)
        im = Image.new("L", (wid, int(2.4 * fs * 2 * INK_K)), 0)
        draw = ImageDraw.Draw(im)
        x = 0.0
        for c in text:
            draw.text((x, base), c, font=font, fill=255, anchor="ls")
            x += hmtx[cmap[ord(c)]][0] * scale
        big = np.asarray(im) > 127
        h2, w2 = big.shape[0] // 2 * 2, big.shape[1] // 2 * 2
        ink = big[:h2, :w2].reshape(h2 // 2, 2, w2 // 2, 2).mean(axis=(1, 3)) > 0.5
        if mirror:
            ink = ink[:, ::-1]
        cols = np.nonzero(ink.any(axis=0))[0]
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        k = INK_K
        self.qmin, self.ymin = INK_LEFT - INK_PAD, -(0.4 * fs + INK_PAD)
        nq = int(math.ceil((c1 - c0) / k + 2 * INK_PAD)) * k
        ny = int(math.ceil((1.0 * fs + 2 * INK_PAD) * k))
        field = np.zeros((ny, nq), dtype=bool)
        r0 = int(round((BASE_EM * fs - self.ymin) * k)) - base // 2        # row of the first mask row in the field
        c_off = int(round(INK_PAD * k)) - c0
        rows = np.arange(ink.shape[0]) + r0
        keep = (rows >= 0) & (rows < ny)
        field[rows[keep], c0 + c_off: c1 + c_off] = ink[keep, c0:c1]
        self.D = (distance_transform_edt(~field) / k).astype(np.float32)
        self.ny, self.nq = self.D.shape

    def look(self, qq, yy, fast: bool = False):
        """Distance to the ink at (qq, yy) (arrays of one shape); beyond the field a lower bound.  fast: nearest grid point."""
        k = INK_K
        gx, gy = (qq - self.qmin) * k - 0.5, (yy - self.ymin) * k - 0.5
        cx, cy = np.clip(gx, 0.0, self.nq - 1.0), np.clip(gy, 0.0, self.ny - 1.0)
        out = np.hypot(gx - cx, gy - cy) / k
        if fast:
            return np.where(out > 0, INK_PAD - 1.0 + out, self.D[np.rint(cy).astype(np.int64), np.rint(cx).astype(np.int64)])
        i0, j0 = np.minimum(cx.astype(np.int64), self.nq - 2), np.minimum(cy.astype(np.int64), self.ny - 2)
        fx, fy = cx - i0, cy - j0
        d = self.D
        v = (d[j0, i0] * (1 - fx) + d[j0, i0 + 1] * fx) * (1 - fy) + (d[j0 + 1, i0] * (1 - fx) + d[j0 + 1, i0 + 1] * fx) * fy
        return np.where(out > 0, INK_PAD - 1.0 + out, v)


def _font_tables():
    if "font" not in _INK:
        tt = TTFont(str(FONT_FILE))
        _INK["font"] = (tt.getBestCmap(), tt["hmtx"], tt["head"].unitsPerEm)
    return _INK["font"]


def _ink_rel(text: str, size: float, anchor: str) -> tuple[float, float, float, float]:
    """Prototype P (LP-D23, LP-D25): the ink box (left, top, right, bottom) of a text drawn at the origin with dominant-baseline="middle":
    ink left to ink right, top of the tallest letter to the lowest ink (what Inkscape shows as the bounding box).  Same model as scripts/lp_ports.py."""
    if "full" not in _INK:
        tt = TTFont(str(FONT_FILE))
        _INK["full"] = (tt.getBestCmap(), tt["glyf"], tt["hmtx"], tt["head"].unitsPerEm, tt["OS/2"].sxHeight)
    cmap, glyf, hmtx, upm, xh = _INK["full"]
    s = size / upm
    adv, x0, x1, top, bot = 0, None, 0, 0, 0
    for ch in text:
        g = cmap[ord(ch)]
        a, _lsb = hmtx[g]
        gl = glyf[g]
        if gl.numberOfContours:
            if x0 is None:
                x0 = adv + gl.xMin
            x1 = adv + gl.xMax
            top, bot = max(top, gl.yMax), min(bot, gl.yMin)
        adv += a
    shift = -adv * s if anchor == "end" else (-adv * s / 2 if anchor == "middle" else 0.0)
    base = xh * s / 2
    return (shift + (x0 or 0) * s, base - top * s, shift + x1 * s, base - bot * s)


def _ink(text: str, fs: float, mirror: bool) -> _Ink:
    key = (text, round(fs, 3), mirror)
    if key not in _INK:
        _INK[key] = _Ink(text, fs, mirror)
    return _INK[key]


class _Par:
    def __init__(self, ctx: Ctx):
        self.gap = float(ctx.gap) + 0.05                 # between two label boxes
        self.clear = ctx.opt("clear", ctx.gap)           # label box to every other object
        self.cone = ctx.opt("cone", 70.0)                # largest turn of a leader at the ring edge (degrees)
        self.cone_cos = math.cos(math.radians(self.cone))
        self.lgap = ctx.opt("lgap", 2.5)                 # leader end to its own text box
        self.lpad = min(ctx.opt("lpad", 2.4), self.lgap - 0.1)   # a leader passing another label's box (below lgap:
        #                                                  a label on its own exit ray never blocks a neighbour's ray)
        self.lpadp = ctx.opt("lpadp", self.lpad)          # a pending label's exit ray to a box
        self.stub = ctx.opt("stub", 1.0)                 # elbow distance beyond the ring edge
        self.edge = ctx.opt("fedge", 1.0)                # allowance for the feature outline stroke
        self.eg = int(ctx.opt("edge", 1))                # 1: the leader may end anywhere on the perimeter of the label box (LP-D13)
        self.back = ctx.opt("back", -1.0)                # px any point of the box may be nearer the ring than the leader end (edge=1); -1: half the text height
        self.fl = ctx.opt("fl", 12.0)                    # a label may sit towards the equator of its feature while a
        #                                                  spot within fl px of the ring exists there
        self.cutmax = math.radians(ctx.opt("cutmax", 30.0))   # the cut at a pole lies within this angle of the pole
        self.cutn = ctx.opt("cutn", 14)                  # ... and among this many label gaps
        self.short = 2.5                                 # a run shorter than this has no direction
        self.tw = ctx.opt("tw", 0.8)                      # px of leader cost per degree of turn at the ring edge
        self.plat = ctx.opt("plat", 2.0)                 # px of leader cost a label may give up to make room for the labels after it
        self.sym = int(ctx.opt("sym", 1))                # 1: a lopsided cluster of labels is shifted until its fan is symmetric
        self.ray = int(ctx.opt("ray", 1))                # 1: candidates exactly on the feature's own ray (straight leaders)
        self.slide = int(ctx.opt("slide", 1))             # 1: a finished layout is slid toward the ring, leaders end where their line meets the box (LP-D14)
        self.shift = int(ctx.opt("shift", 0))             # 1: the slide pass also moves a box straight toward the ring, its leader keeping its direction (LP-D14 probe)
        self.face = int(ctx.opt("face", 1))              # 1: a leader ends on the half of its box that faces it, at one distance (`ink`) from the ink of its text,
        #                                                  and nearer to its own ink than to any other label's (LP-D15)
        self.relax = ctx.opt("relax", 40.0)               # px of push (and of run + push saved) at which a label is seated as in face=0 (0: never)
        self.strt = int(ctx.opt("strt", 1))               # 1: a bent leader blocked by its neighbour's box is straightened by sliding that neighbour outward (LP-D16)
        self.mid = int(ctx.opt("mid", 1))                # 1: a leader that reaches the left or right end of its box head-on ends at the middle of that end (LP-D18)
        self.midfr = ctx.opt("midfr", 0.3)               # ... within this fraction of the half height of the middle (the metric allows 0.35)
        self.midang = ctx.opt("midang", 11.0)            # ... "head-on" = within this many degrees of the horizontal (the metric: 10)
        self.fan = int(ctx.opt("fan", 1))                # 1: a chain of labels with a bent leader is seated again as a straight fan, ends at the middle of the box end (round 3)
        self.fanx = ctx.opt("fanx", 40.0)                # ... a leader may be this many px longer than before
        self.ink = ctx.opt("ink", 3.4)                   # wanted distance of a leader end to the ink of its own text (px)
        self.inkmin, self.inkmax = ctx.opt("inkmin", 2.5), ctx.opt("inkmax", 4.5)   # ... the end is refused outside this range
        self.um = ctx.opt("um", 0.7)                     # an end is nearer to its own ink than to another label's ink by 1 px + this margin
        self.port = int(ctx.opt("port", 0))              # 1: prototype P (LP-D22 to LP-D25): ink boxes, 12 contact points, straight-first clusters; replaces face, mid, strt, fan, sym, slide
        self.parc = ctx.opt("parc", 60.0)                # P: a straight label may sit up to this many px farther along its ray than its nearest allowed seat (the arc of its cluster may ask for it; R5 prefers the nearest)
        self.bend = ctx.opt("bend", 6.0)                 # P: px of cost of a bent leader in the cluster construction
        self.minlast = ctx.opt("minlast", 6.0)           # P: the last piece of a bent leader is at least this long (drawn)
        self.pwin = ctx.opt("pwin", 30.0)                # P: rows a label of a cluster may move (px each way)
        self.pcone = ctx.opt("pcone", 44.0)              # P: the way out kept free for a label not seated yet leaves its ray this many degrees at most (the contact points take 45)
        self.pvert = ctx.opt("pvert", 4.0)               # P: px of cost of ending at a vertex rather than at the middle of a side (LP-D18)
        self.kink = ctx.opt("kink", 8.0)                 # P round 2: px of cost of a bend under 8 degrees (it fades out to 8 degrees)
        self.noport = ctx.opt("noport", 400.0)            # P round 2: px of cost of keeping an end that is at no contact point
        self.cluster = int(ctx.opt("cluster", 1))        # P: 1: the cluster construction runs
        self.arcw = ctx.opt("arcw", 4.0)                 # P round 4 (R8): px of cost per px that the leader ends of two neighbours differ in distance from the centre beyond `arcstep`
        self.arcstep = ctx.opt("arcstep", 2.0)           # P round 4: the step of scripts/arc_fit.py
        self.loosew = ctx.opt("loosew", 1.0)             # P round 4: px of cost per px that two neighbours are apart beyond `rowgap` when one of them is bent (at most `loosemax` px counted)
        self.rowgap = ctx.opt("rowgap", 3.6)
        self.loosemax = ctx.opt("loosemax", 16.0)
        self.tsq = ctx.opt("tsq", 0.01)                  # P round 4: px of cost per square degree of a bend (shares a displacement between both ends of a cluster)
        self.seats = int(ctx.opt("seats", 2))             # P round 4: seats kept per label of a chain for each row of 3 px (the cheapest, bent and straight apart)
        self.arcfirst = int(ctx.opt("arcfirst", 0))      # Q: 1: side clusters are built in arc coordinates (arc_chain) instead of the seat DP
        self.arcslope = ctx.opt("arcslope", 1.9)         # Q: the arc radius of two neighbours differs by at most this many px
        self.qgap = ctx.opt("qgap", 20.0)                 # Q: boxes at most this many px apart are one run (labels with room stay on their own ray)
        self.qnq = ctx.opt("qnq", 0.34)                  # Q: a run is built in arc coordinates when every ray is within acos(qnq) degrees of the horizontal
        self.qrhomax = ctx.opt("qrhomax", 120.0)         # Q2: a long run's arc may move out this many px
        self.maxturn = ctx.opt("maxturn", 40.0)          # Q: a bent leader turns at most this many degrees
        self.rhomax = ctx.opt("rhomax", 40.0)            # Q: the arc of a run moves out at most this many px from its smallest radius
        if self.port:
            self.face = self.mid = self.fan = self.strt = self.sym = self.slide = 0
        self.font = float(ctx.font_size)
        self.nostraight = 3.0                            # a run shorter than this does not turn (the metric drops pieces under 3 px)
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
        self.memo_r: dict = {}
        self.fld = _ink(it.text, float(it.src.get("font_size", par.font)), sx < 0) if par.face else None
        self.lo, self.hi = -math.inf, math.inf
        # the way out that is kept free until the label is seated: straight along the stack direction, inside the cone
        lim = math.radians(par.pcone if par.port else par.cone - 1.5)      # P: the reserved way out is a leader that a contact point can take
        an, av = math.atan2(self.ny, self.nq), math.atan2(self.d, 0.0)
        turn = (av - an + math.pi) % (2.0 * math.pi) - math.pi
        a = an + max(-lim, min(lim, turn))
        hq, hy = math.cos(a), math.sin(a)
        self.exit, self.slanted = (hq, hy), abs(turn) > lim
        self.edges = [(math.cos(an + s * lim), math.sin(an + s * lim)) for s in (-1.0, 1.0)]   # the two edges of the cone
        reach = 6000.0 if shore == "out" else max(abs(self.Sy) / max(abs(hy), 1e-6) - 1.0, 1.0)
        self.ray = (self.S, (self.S[0] + sx * hq * reach, self.S[1] + hy * reach))
        self._profile(geo, par)

    def cands(self, yc, geo: _Geo, par: _Par, lean: float | None = None):
        """Candidates, shape (len(yc), n): near end q0, attach height ya, leader-to-text gap; and the nearest allowed q0."""
        yc = yc[:, None]
        if par.port:
            qs, ya, ga, base, _pid = self.cands_port(yc, geo, par)
            return qs, ya, ga, base
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
        offs = (OFFS_F if par.face else OFFS)[None, :]
        parts = [base + offs, qp + ady * LEAN[None, :], qp - ady * LEAN[None, :], qp - OFFS[None, 1:12] + zero]
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
        if par.eg:                                         # the leader ends on the top or the bottom edge, or at the point of the box nearest to the stub
            lg = par.lgap
            for uq, uy in [(self.nq, self.ny)] + self.edges:         # the feature's own ray, and the two edges of the cone
                if abs(uy) < 1e-6:
                    continue
                for ye in (y0 - lg, y1 + lg):
                    t = (ye - self.Sy) / uy
                    eq = self.Sq + t * uq
                    for f in EDGE_FR:
                        parts.append(np.where(t > 0.5, eq - f * self.w, base))
                        heights.append(np.where(t > 0.5, ye, near))
                        gaps.append(np.where(t > 0.5, -f * self.w, par.lgap))
            for q0 in [base + offs] + [self.Sq - f * self.w + zero for f in EDGE_FR]:
                aq = np.clip(self.Sq, q0 - lg, q0 + self.w + lg)
                parts.append(q0 + zero)
                heights.append(np.clip(self.Sy, y0 - lg, y1 + lg) + 0.0 * q0)
                gaps.append(q0 - aq)
        if par.face:                                       # the leader meets the box side at several heights (the facing half is kept by the test in static)
            for f in SIDE_FR:
                parts.append(base + OFFS_F[None, :31])
                heights.append(y0 + f * self.h + 0.0 * parts[-1])
                gaps.append(par.lgap + 0.0 * parts[-1])
        if lean is not None:                               # a symmetric fan: leaders that lean by `lean` degrees (image coordinates, clockwise > 0)
            th = math.radians(self.sx * lean)
            vq, vy = self.nq * math.cos(th) - self.ny * math.sin(th), self.nq * math.sin(th) + self.ny * math.cos(th)
            if abs(vy) > 1e-6:
                for fy in (0.0, 0.25, 0.5, 0.75, 1.0):
                    ya = y0 + fy * self.h
                    run = (ya - self.Sy) / vy
                    parts.append(np.where(run > 0.5, self.Sq + run * vq + par.lgap, base))
                    heights.append(np.where(run > 0.5, ya, near))
                    gaps.append(np.full_like(base, par.lgap))
                if par.eg:
                    for ye in (y0 - par.lgap, y1 + par.lgap):
                        t = (ye - self.Sy) / vy
                        eq = self.Sq + t * vq
                        for f in EDGE_FR:
                            parts.append(np.where(t > 0.5, eq - f * self.w, base))
                            heights.append(np.where(t > 0.5, ye, near))
                            gaps.append(np.where(t > 0.5, -f * self.w, par.lgap))
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
        qs, ya, ga = np.concatenate(parts, axis=1), np.concatenate(heights, axis=1), np.concatenate(gaps, axis=1)
        if par.face:
            ya, ga = self.snap(yc, qs, ya, ga, par)
        return qs, ya, ga, base

    def port_end(self, j: int, yc, q0, par: _Par):
        """Prototype P: the leader S -> contact point j of the ink box (q0, yc +- h/2) stops lgap before the point.  (attach height, gap q0 - attach q);
        NaN where the direction may not use that point (round 2: within 45 degrees of a short side's normal only the middle of that side; 45 to 70 degrees
        a vertex from outside its corner or a long-side point; steeper than 70 degrees a long-side point only)."""
        fq, fy = PORTS[j]
        yc, q0 = np.broadcast_arrays(np.asarray(yc, dtype=float), np.asarray(q0, dtype=float))
        pq, py = q0 + fq * self.w, yc + (fy - 0.5) * self.h
        dq, dy = pq - self.Sq, py - self.Sy
        n = np.maximum(np.hypot(dq, dy), 1e-9)
        ax, ay = -dq / n, -dy / n                                          # where the leader comes from, seen from the point
        nqn = -1.0 if fq == 0.0 else (1.0 if fq == 1.0 else 0.0)
        nyn = -1.0 if fy == 0.0 else (1.0 if fy == 1.0 else 0.0)
        deg = np.degrees(np.arctan2(np.abs(ay), np.abs(ax)))                # 0 = along the normal of a short side (horizontal text)
        if nqn and nyn:                                                     # a vertex: from outside its corner, 45 to 70 degrees
            ok = (ax * nqn > 0) & (ay * nyn > 0) & (deg > 45.0) & (deg <= 70.0)
        elif nqn:                                                           # the middle of a short side: within 45 degrees of its normal, and only it
            ok = (ax * nqn > 0) & (deg <= 45.0)
        else:                                                               # a long-side point: within 45 degrees of its normal (steeper than 45 degrees)
            ok = (ay * nyn > 0) & (deg >= 45.0)
        ok &= n > par.lgap + 0.05
        aq, ya = pq - par.lgap * dq / n, py - par.lgap * dy / n
        return np.where(ok, ya, np.nan), np.where(ok, q0 - aq, np.nan)

    def cands_port(self, yc, geo: _Geo, par: _Par, offs=OFFS):
        """Prototype P: candidates (q0, attach height, gap, nearest allowed q0, contact point) for every contact point, each at the box
        positions base + offs and, where the point lies on the feature's own ray at this row, on that ray (a straight leader)."""
        h = self.h
        if self.shore == "out":
            dyn = np.maximum(np.abs(yc) - 0.5 * h - geo.padw, 0.0)
            hug = np.where(dyn < geo.R, np.sqrt(np.maximum(geo.R ** 2 - dyn ** 2, 0.0)) + geo.padw + 1e-6, -1e9)
        else:
            dyf = np.abs(yc) + 0.5 * h + geo.padw
            hug = np.where(dyf < geo.Bb, geo.padw + 1e-6 - np.sqrt(np.maximum(geo.Bb ** 2 - dyf ** 2, 0.0)), 1e9)
        base = np.maximum(self.Sq + par.lgap, hug)
        qs, heights, gaps, pids = [], [], [], []
        for j, (fq, fy) in enumerate(PORTS):
            q0 = base + np.asarray(offs)[None, :] + 0.0 * yc
            if abs(self.ny) > 1e-6:                           # the point on the feature's own ray at this row
                t = (yc + (fy - 0.5) * h - self.Sy) / self.ny
                q0 = np.concatenate([q0, np.where(t > 0.5, self.Sq + t * self.nq - fq * self.w, base)], axis=1)
            ya, ga = self.port_end(j, yc, q0, par)
            qs.append(q0)
            heights.append(ya)
            gaps.append(ga)
            pids.append(np.full(q0.shape, j))
        return (np.concatenate(qs, axis=1), np.concatenate(heights, axis=1), np.concatenate(gaps, axis=1), base,
                np.concatenate(pids, axis=1))

    def snap(self, yc, q0, ya, ga, par: _Par):
        """face=1: the end of every candidate moves along its own line (stub S through the old end) to where it stands `ink` px from
        the ink of its own text, or touches the box when the ink stays farther inside; NaN gap where that is not within
        [inkmin, inkmax] (the line passes the text, or the ink is too far from the box side it meets).  Sphere tracing on the field."""
        aq = q0 - ga
        vq, vy = aq - self.Sq, ya - self.Sy
        ln = np.hypot(vq, vy)
        near = ln < 0.5
        uq, uy = np.where(near, self.nq, vq / np.maximum(ln, 1e-9)), np.where(near, self.ny, vy / np.maximum(ln, 1e-9))
        t = np.zeros_like(q0)
        half = 0.5 * self.h

        def field(t, fast):
            rq, ry = self.Sq + t * uq - q0, self.Sy + t * uy - yc
            box = np.hypot(np.maximum(np.maximum(-rq, rq - self.w), 0.0), np.maximum(np.maximum(-ry - half, ry - half), 0.0))
            d = self.fld.look(rq, ry, fast)
            return d, np.minimum(d - par.ink, box)

        for _ in range(9):
            t = t + np.maximum(field(t, True)[1], 0.0)
        for _ in range(2):
            t = t + np.maximum(field(t, False)[1], 0.0)
        d, f = field(t, False)
        ok = (f <= 0.25) & (d >= par.inkmin) & (d <= par.inkmax) & np.isfinite(ln)
        return np.where(ok, self.Sy + t * uy, ya), np.where(ok, q0 - (self.Sq + t * uq), np.nan)

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
        if par.port:                                       # R6: the last piece of a bent leader is at least `minlast` long
            ok &= ~((cosv < math.cos(math.radians(1.5))) & (ln < par.minlast))
        if par.face:                                       # LP-D15 (the metric's own test, a little stricter): the end is on the half of the box that faces the leader
            s10 = math.sin(math.radians(9.5))
            seen = ln >= 3.0
            fq, fy = np.where(seen, dq / np.maximum(ln, 1e-9), self.nq), np.where(seen, dy / np.maximum(ln, 1e-9), self.ny)
            outq, outy = (aq <= q0 + 0.05) | (aq >= q0 + w - 0.05), (ya <= y0 + 0.05) | (ya >= y1 - 0.05)
            ok &= ~((outq & (np.abs(fy) >= s10) & ((ya - yc) * fy > 0) & (np.abs(ya - yc) > 0.9))
                    | (outy & (np.abs(fq) >= s10) & ((aq - q0 - 0.5 * w) * fq > 0) & (np.abs(aq - q0 - 0.5 * w) > 0.9)))
        if par.face and par.mid:                           # LP-D18: head-on onto the left or right end of the box: at the middle of it
            ok &= ~(outq & (np.abs(fy) < math.sin(math.radians(par.midang))) & (np.abs(ya - yc) > par.midfr * 0.5 * h))
        if par.eg:
            lg = par.lgap
            ok &= ~((self.Sq > q0 - lg) & (self.Sq < q0 + w + lg) & (self.Sy > y0 - lg) & (self.Sy < y1 + lg))   # the stub inside the grown box
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
        if par.eg:                                         # LP-D13: no point of the box nearer the ring than the leader end by more than `back`
            if self.shore == "out":
                mn = np.hypot(np.clip(0.0, q0, q0 + w), np.clip(0.0, y0, y1)) - ring
            else:
                mn = ring - np.hypot(np.maximum(np.abs(q0), np.abs(q0 + w)), np.maximum(np.abs(y0), np.abs(y1)))
            back = par.back if par.back >= 0 else 0.5 * h
            ok &= (np.abs(np.hypot(aq, ya) - ring) - mn) <= back + 1e-9
        else:
            ok &= ~(far < near - 0.7)
        if geo.has[self.shore]:
            x0, x1 = (q0, q0 + w) if self.sx > 0 else (-(q0 + w), -q0)
            r = geo.padr
            ok &= ~geo.blocked(x0 - r, y0 - r + 0.0 * q0, x1 + r, y1 + r + 0.0 * q0)
        return ok, cost, ln

    def _profile(self, geo: _Geo, par: _Par):
        """Alone: the y with the shortest leader (target) and how far the label may drift from it."""
        G = par.grid_out if self.shore == "out" else par.grid_in
        ys = self.Sy + np.arange(-G, G + 1.0, 2.0 if par.face else 1.0)
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
        self._lab: list = []
        self._arr = None

    def add(self, box, s, a, lab=None):
        self._b.append(box)
        self._s0.append(s)
        self._s1.append(a)
        self._lab.append(lab)
        self._arr = None

    def labs(self):
        return self._lab

    def copy(self) -> "_Placed":
        c = _Placed()
        c._b, c._s0, c._s1, c._lab = list(self._b), list(self._s0), list(self._s1), list(self._lab)
        return c

    def arrays(self):
        if self._arr is None:
            self._arr = (np.array(self._b, dtype=float).reshape(-1, 4), np.array(self._s0, dtype=float).reshape(-1, 2),
                         np.array(self._s1, dtype=float).reshape(-1, 2))
        return self._arr


class _View:
    """_Placed without row k (the label being slid)."""

    def __init__(self, pl: "_Placed", k: int):
        self.pl, self.k = pl, k

    def arrays(self):
        pb, s0, s1 = self.pl.arrays()
        keep = np.arange(len(pb)) != self.k
        return pb[keep], s0[keep], s1[keep]

    def labs(self):
        return [lb for i, lb in enumerate(self.pl.labs()) if i != self.k]


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
        self.qcache: dict = {}
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

    def spots(self, L: _Lab, yy: float, want: float | None = None):
        """The spots of this label at height yy with no other label placed: (q, cost, boxes, attach points, least cost) or None."""
        geo, par = self.geo, self.par
        y0, y1 = yy - 0.5 * L.h, yy + 0.5 * L.h
        key = round(yy, 4) if want is None else (round(yy, 4), round(want, 2))
        if key not in L.memo:
            ys = np.array([yy])
            Q, YA, GA, _base = L.cands(ys, geo, par, want)
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
        return L.memo[key]

    # ---- seat one label at y (or the first y further along its stack direction where it fits)
    def try_at(self, L: _Lab, yy: float, pl: _Placed, pend: _Pending, want: float | None = None, extra=None, only: bool = False):
        """The shortest leader at height yy that fits: ((box, attach point, yy, run), shortest run before the other
        labels are looked at). (None, None) when the label has no spot at this height at all. `extra`: more spots
        (q, cost, boxes, attach points) that have passed the static tests; `only`: these and no others."""
        par = self.par
        y0, y1 = yy - 0.5 * L.h, yy + 0.5 * L.h
        memo = None if only else self.spots(L, yy, want)
        if extra is not None and len(extra[0]):
            memo = tuple(np.concatenate([a, b]) for a, b in zip(memo[:4], extra)) + (memo[4],) if memo is not None else (*extra, math.inf)
        if memo is None:
            return None, None
        q, ln, boxes, att, alone = memo
        stub = np.broadcast_to(np.array(L.S, dtype=float), att.shape)
        lo, hi = min(L.Sy, y0) - par.lpad - 1.0, max(L.Sy, y1) + par.lpad + 1.0     # rows the boxes and leaders touch
        good = self.legal(L, boxes, att, stub, pl, pend, lo, hi)
        if not good.any():
            return None, alone
        sel = np.nonzero(good)[0]
        cost = ln[sel] + 1e-6 * np.abs(q[sel] - (L.Sq + par.lgap))
        if want is not None:                              # a cluster being made symmetric: this label should lean by `want` degrees
            cost = cost + SYM_K * (_sturn_arr(L.it.anchor, L.S, att[sel]) - want) ** 2
        j = sel[int(np.argmin(cost))]
        return (tuple(float(v) for v in boxes[j]), (float(att[j][0]), float(att[j][1])), yy, float(ln[j])), alone

    def legal(self, L: _Lab, boxes, att, stub, pl: _Placed, pend: _Pending, lo: float, hi: float) -> np.ndarray:
        """Which candidate boxes and leaders pass the tests against the labels placed (`pl`), the exit rays of the labels still to be seated
        (`pend`), the other objects and (face=1) the ink: box gap, leader to box clearance, leader crossings.  lo, hi: rows the candidates touch."""
        geo, par = self.geo, self.par
        good = np.ones(len(boxes), dtype=bool)
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
            bad = _seg_box_hit(rp, rq, boxes[sel], par.lpadp).any(axis=0)
            bad |= _seg_cross_matrix(stub[sel], att[sel], rp, rq).any(axis=1)
            good[sel[bad]] = False
        if geo.has[L.shore] and good.any():
            sel = np.nonzero(good)[0]
            good[sel[geo.seg_blocked(L.S, att[sel])]] = False
        if par.face and good.any():
            sel = np.nonzero(good)[0]
            good[sel[~self.ink_clear(L, boxes[sel], att[sel], pl)]] = False
        return good

    def ink_clear(self, L: _Lab, boxes, att, pl) -> np.ndarray:
        """face=1 (LP-D15): the end of a leader is nearer to the ink of its own text than to the ink of any other label by more than
        1 px (+ `um`), for the candidate's end against the labels placed and for their ends against the candidate's text."""
        pb, _s0, s1 = pl.arrays()
        ok = np.ones(len(boxes), dtype=bool)
        if not len(pb):
            return ok
        labs, need, sx = pl.labs(), 1.0 + self.par.um, L.sx
        cq0 = boxes[:, 0] if sx > 0 else -boxes[:, 2]
        cym = 0.5 * (boxes[:, 1] + boxes[:, 3])
        own = L.fld.look(sx * att[:, 0] - cq0, att[:, 1] - cym)
        lo_x, hi_x = min(boxes[:, 0].min(), att[:, 0].min()) - 24.0, max(boxes[:, 2].max(), att[:, 0].max()) + 24.0
        lo_y, hi_y = min(boxes[:, 1].min(), att[:, 1].min()) - 24.0, max(boxes[:, 3].max(), att[:, 1].max()) + 24.0
        for j in np.nonzero((pb[:, 2] > lo_x) & (pb[:, 0] < hi_x) & (pb[:, 3] > lo_y) & (pb[:, 1] < hi_y))[0]:
            M = labs[j]
            if M is None:
                continue
            mq0, mym = (pb[j, 0] if M.sx > 0 else -pb[j, 2]), 0.5 * (pb[j, 1] + pb[j, 3])
            ok &= M.fld.look(M.sx * att[:, 0] - mq0, att[:, 1] - mym) >= own + need                  # my end to its ink
            own_j = M.fld.look(np.array(M.sx * s1[j, 0] - mq0), np.array(s1[j, 1] - mym))             # its end to my ink
            ok &= L.fld.look(sx * s1[j, 0] - cq0, s1[j, 1] - cym) >= own_j + need
        return ok

    def seat(self, L: _Lab, y: float, pl: _Placed, pend: _Pending, limit: float | None, tight: float | None):
        """Seat the label at y or further along its stack: the smallest leader run + push.
        `tight` is the push that puts it exactly one row beyond the label seated before it in the stack.
        face=1: a label that the facing rule pushes more than `relax` px along its stack (a fan at its geometric limit) is seated
        as in face=0 when that is more than `relax` px cheaper (run + push); the metrics count such a label as a residual."""
        par = self.par
        best = self._search(L, y, pl, pend, limit, tight)
        if par.face and par.relax > 0 and (best is None or best[1] > par.relax):
            par.face, L.memo, L.memo_r = 0, L.memo_r, L.memo
            try:
                alt = self._search(L, y, pl, pend, limit, tight)
            finally:
                par.face, L.memo, L.memo_r = 1, L.memo_r, L.memo
            if alt is not None and (best is None or alt[0] + par.relax < best[0]):
                best = alt
                self.stats["relaxed"] = self.stats.get("relaxed", 0) + 1
        if par.port and par.relax > 0 and (best is None or best[1] > par.relax):     # P: no contact-point seat within reach: the J2 seat (a leader ends on the box edge); counted as relaxed
            par.port, L.memo, L.memo_r = 0, L.memo_r, L.memo
            try:
                alt = self._search(L, y, pl, pend, limit, tight)
            finally:
                par.port, L.memo, L.memo_r = 1, L.memo_r, L.memo
            if alt is not None and (best is None or alt[0] + par.relax < best[0]):
                best = alt
                self.stats["relaxed"] = self.stats.get("relaxed", 0) + 1
        if best is None:
            return None
        if best[1] > 0:
            self.stats["pushed"] += 1
            self.stats["push_max"] = max(self.stats["push_max"], best[1])
        return best[2]

    def _search(self, L: _Lab, y: float, pl: _Placed, pend: _Pending, limit: float | None, tight: float | None):
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
        return best

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
                pl.add(box, L.S, att, L)
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

    # ---- symmetric fans (J2)
    def symmetrize(self, recs: list, shore: str) -> tuple[list, dict]:
        """A cluster is a chain of labels of one stack whose boxes touch (gap <= box gap + 3 px) with a leader that turns
        at least 3 degrees. It is lopsided when |sum of signed turns| > half the sum of |turns| and the latter is at least
        6 degrees. A lopsided cluster is shifted rigidly along its stack (heights of all its labels + delta), every label is
        seated again at its shifted height, and the shift that balances the signed turns (largest turn not above the
        old one) is kept. Everything else stays where it was."""
        par = self.par
        recs = list(recs)
        info = {"clusters": 0, "lopsided": 0, "fixed": 0}
        if not recs:
            return recs, info
        labs = [self.lab(r[0], shore, 1 if r[3] == "start" else -1) for r in recs]
        pend = _Pending(labs)
        for lb in labs:
            pend.done(lb)
        groups: dict = {}
        for k, lb in enumerate(labs):
            groups.setdefault((lb.sx, lb.upper), []).append(k)
        chains = []
        for ks in groups.values():
            ks.sort(key=lambda k: recs[k][1][1])
            cur = [ks[0]]
            for k in ks[1:]:
                if recs[k][1][1] - recs[cur[-1]][1][3] <= par.gap + 3.0:
                    cur.append(k)
                else:
                    chains.append(cur)
                    cur = [k]
            chains.append(cur)
        for chain in chains:
            if len(chain) < 2:
                continue
            sturn = [_sturn(recs[k][2]) for k in chain]
            if max(abs(v) for v in sturn) < 3.0:
                continue
            info["clusters"] += 1
            tot, mag = sum(sturn), sum(abs(v) for v in sturn)
            if not (abs(tot) > 0.5 * mag and mag >= 6.0):
                continue
            info["lopsided"] += 1
            m = len(chain)
            max0 = max(abs(v) for v in sturn)
            len0 = sum(_plen(recs[k][2]) + par.tw * abs(v) for k, v in zip(chain, sturn))
            yc = {k: 0.5 * (recs[k][1][1] + recs[k][1][3]) for k in chain}
            span = recs[chain[-1]][1][3] - recs[chain[0]][1][1]
            reach = min(0.4 * span + 8.0, 40.0) if m <= 8 else min(0.5 * span, 200.0)
            stp = 4.0 if m <= 8 else reach / 8.0
            tries = [0.0] + [v * x for x in np.arange(stp, reach + 0.1, stp) for v in (1.0, -1.0)]
            big = max(range(m), key=lambda i: abs(sturn[i]))      # the end that carries the turn leans one way, the other end the other way
            sgn = (1.0 if sturn[big] > 0 else -1.0) * (1.0 if big >= (m - 1) / 2 else -1.0)
            cand = []
            for amp in (0.5 * max0, max0):
                ramp = []
                for i in range(m):
                    u = 2.0 * i / (m - 1) - 1.0                   # -1 .. 1 along the stack; the inner part stays straight
                    ramp.append(sgn * amp * math.copysign(max(0.0, (abs(u) - 0.4) / 0.6), u))
                for dl in tries:
                    cost, okd = 0.0, True
                    for i, k in enumerate(chain):
                        mm = self.spots(labs[k], yc[k] + dl, ramp[i])
                        if mm is None:
                            okd = False
                            break
                        ang = _sturn_arr(recs[k][0].anchor, labs[k].S, mm[3])
                        cost += float(np.min(mm[1] + SYM_K * (ang - ramp[i]) ** 2))
                    if okd:
                        cand.append((cost, dl, amp, ramp))
            win = None
            for _cost, dl, _amp, ramp in sorted(cand, key=lambda v: (v[0], abs(v[1])))[:(30 if m <= 8 else 10)]:
                pl = _Placed()
                for k in range(len(recs)):
                    if k not in chain:
                        pl.add(recs[k][1], labs[k].S, recs[k][2][-1], labs[k])
                new, fail = {}, False
                for k in sorted(chain, key=lambda k: yc[k] * labs[k].d):
                    res, _alone = self.try_at(labs[k], yc[k] + dl, pl, pend, ramp[chain.index(k)])
                    if res is None:
                        fail = True
                        break
                    box, att, _yy, _c = res
                    pl.add(box, labs[k].S, att, labs[k])
                    new[k] = (box, [recs[k][0].anchor, labs[k].S, att])
                if fail:
                    continue
                st = [_sturn(new[k][1]) for k in chain]
                len1 = sum(_plen(new[k][1]) + par.tw * abs(v) for k, v in zip(chain, st))
                mag1 = sum(abs(v) for v in st)
                if abs(sum(st)) <= abs(tot) - 3.0 and max(abs(v) for v in st) <= max0 + 0.5 and len1 <= len0 + 6.0 * m:
                    key = (round(abs(sum(st)) / max(mag1, 1.0), 1), len1)
                    if win is None or key < win[0]:
                        win = (key, new)
            if win is not None:
                for k in chain:
                    it, _b, _p, align, tag = recs[k]
                    recs[k] = (it, win[1][k][0], win[1][k][1], align, tag)
                info["fixed"] += 1
        return recs, info

    # ---- straighten a bent leader by sliding the neighbour that blocks its ray (LP-D16)
    def straighten(self, recs: list, shore: str) -> tuple[list, dict]:
        """A leader that bends (turn at S of 3 degrees or more) may be blocked only by the box of its neighbour in the stack.  For
        every bent label and each of its two neighbours: the neighbour is slid outward along its own leader (box and end move together,
        so the leader keeps its direction, the end stays on the same half of the box and at the same gap to the ink; 2..15 px,
        a few px longer), the bent label is seated again at heights within 5 px of its own against the new layout (try_at, all
        the usual tests), and the pair is kept when the sum of their leader costs (run + tw x turn) falls by at least 0.5 px.
        Rounds repeat so the step moves down the stack.  Straight leaders are preferred to a symmetric bent fan."""
        par = self.par
        recs = list(recs)
        info = {"tried": 0, "moved": 0, "straightened": 0, "added_px": 0.0}
        if not recs:
            return recs, info
        labs = [self.lab(r[0], shore, 1 if r[3] == "start" else -1) for r in recs]
        pend = _Pending(labs)
        for lb in labs:
            pend.done(lb)

        def cost(L: _Lab, att) -> float:
            vx, vy = att[0] - L.S[0], att[1] - L.S[1]
            ln = math.hypot(vx, vy)
            if ln < par.nostraight:
                return ln
            cosv = max(-1.0, min(1.0, (vx * L.nq * L.sx + vy * L.ny) / ln))
            return ln + par.tw * math.degrees(math.acos(cosv))

        def around(skip, centre, reach=70.0) -> _Placed:
            pl = _Placed()
            x0, y0, x1, y1 = centre[0] - reach, centre[1] - reach, centre[2] + reach, centre[3] + reach
            for j, r in enumerate(recs):
                if j in skip:
                    continue
                bx0, by0 = min(r[1][0], r[2][1][0], r[2][-1][0]), min(r[1][1], r[2][1][1], r[2][-1][1])
                bx1, by1 = max(r[1][2], r[2][1][0], r[2][-1][0]), max(r[1][3], r[2][1][1], r[2][-1][1])
                if bx1 > x0 and bx0 < x1 and by1 > y0 and by0 < y1:
                    pl.add(r[1], labs[j].S, r[2][-1], labs[j])
            return pl

        groups: dict = {}
        for k, lb in enumerate(labs):
            groups.setdefault((lb.sx, lb.upper), []).append(k)
        for _round in range(3):
            moved = 0
            for ks in groups.values():
                ks.sort(key=lambda k: 0.5 * (recs[k][1][1] + recs[k][1][3]))
                for pos, k in enumerate(ks):
                    B, itb, ptsb = labs[k], recs[k][0], recs[k][2]
                    if not 3.0 <= abs(_sturn(ptsb)) <= 40.0:      # a few px of slide cannot straighten a leader that turns more
                        continue
                    for sgn, clen in ((-1, 1), (1, 1), (-1, 4), (1, 4)):      # the neighbour above or below; with the labels behind the bent one (up to 3)
                        if not 0 <= pos + sgn < len(ks):
                            continue
                        n = ks[pos + sgn]
                        chain = [ks[i] for i in range(pos, len(ks) if sgn < 0 else -1, -sgn)][:clen]      # B and the labels further from N
                        N, (itn, boxn, ptsn, alignn, tagn) = labs[n], recs[n]
                        vx, vy = ptsn[-1][0] - N.S[0], ptsn[-1][1] - N.S[1]
                        d = math.hypot(vx, vy)
                        ux, uy = (vx / d, vy / d) if d >= 3.0 else (N.sx * N.nq, N.ny)
                        old = cost(N, ptsn[-1]) + sum(cost(labs[j], recs[j][2][-1]) for j in chain)
                        done = False
                        base = None
                        for step in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0):
                            info["tried"] += 1
                            nb = (boxn[0] + step * ux, boxn[1] + step * uy, boxn[2] + step * ux, boxn[3] + step * uy)
                            na = (ptsn[-1][0] + step * ux, ptsn[-1][1] + step * uy)
                            q0, ycn = (nb[0] if N.sx > 0 else -nb[2]), 0.5 * (nb[1] + nb[3])
                            ok, c, _ln = N.static(np.array([[ycn]]), np.array([[q0]]), np.array([[na[1]]]), np.array([[q0 - N.sx * na[0]]]), self.geo, par, par.cone_cos)
                            if not ok[0, 0]:
                                continue
                            if base is None:                   # the layout around the pair, without the slid label and the chain
                                span = [boxn] + [recs[j][1] for j in chain]
                                base = around({n, *chain}, (min(b[0] for b in span) - 20, min(b[1] for b in span) - 20, max(b[2] for b in span) + 20, max(b[3] for b in span) + 20), 90.0)
                            pl = base.copy()
                            extra = (np.array([q0]), np.array([c[0, 0]]), np.array([nb]), np.array([na]))
                            if self.try_at(N, ycn, pl, pend, extra=extra, only=True)[0] is None:
                                continue
                            pl.add(nb, N.S, na, N)
                            new, total = {}, c[0, 0]
                            for j in chain:                    # seat the chain again, nearest to N first, each at the cheapest spot within a few px of its own
                                Lj, rj = labs[j], recs[j]
                                yj, best = 0.5 * (rj[1][1] + rj[1][3]), None
                                for dy in (0.0, 1.0, -1.0, 2.0, -2.0, 4.0, -4.0, 8.0, -8.0):
                                    res = self.try_at(Lj, yj + dy, pl, pend)[0]
                                    if res is not None and (best is None or res[3] < best[3] - 1e-9):
                                        best = res
                                if best is None:
                                    break
                                pl.add(best[0], Lj.S, best[1], Lj)
                                new[j] = best
                                total += best[3]
                            npts = [itb.anchor, B.S, new[k][1]] if k in new else None
                            if len(new) == len(chain) and npts is not None and abs(_sturn(npts)) < abs(_sturn(ptsb)) - 1.0 and total < old - 0.5:
                                info["added_px"] += _plen([itn.anchor, N.S, na]) - _plen(ptsn) + sum(_plen([recs[j][0].anchor, labs[j].S, new[j][1]]) - _plen(recs[j][2]) for j in chain)
                                info["straightened"] += abs(_sturn(npts)) < 3.0
                                for j in chain:
                                    recs[j] = (recs[j][0], new[j][0], [recs[j][0].anchor, labs[j].S, new[j][1]], recs[j][3], recs[j][4])
                                recs[n] = (itn, nb, [itn.anchor, N.S, na], alignn, tagn)
                                info["moved"] += 1
                                moved += 1
                                done = True
                                break
                        if done:
                            break
            if not moved:
                break
        info["added_px"] = round(info["added_px"], 1)
        return recs, info

    # ---- straight fans (round 3)
    def fan_seat(self, L: _Lab, rec, pl: _Placed, pend: _Pending):
        """The leader along its own first piece (straight), ending at the middle of the near end of the box (then 0.35 and 0.25 of the
        height from the facing side), at the smallest run (0.5 px steps, at most `fanx` px longer than before) that passes every test
        against `pl`.  (box, end, yc, run) or None."""
        par, geo = self.par, self.geo
        pts = rec[2]
        vx, vy = L.S[0] - pts[0][0], L.S[1] - pts[0][1]
        nv = math.hypot(vx, vy)
        uq, uy = L.sx * vx / nv, vy / nv
        ts = np.arange(0.0, math.dist(L.S, pts[-1]) + par.fanx + 1e-9, 0.5)
        h, w = L.h, L.w
        for frac in (0.5, 0.35, 0.25):
            yc = L.Sy + ts * uy + (0.5 - frac) * h * (1.0 if uy > 0 else -1.0)
            q0 = L.Sq + ts * uq + par.lgap
            ya, ga = L.snap(yc[:, None], q0[:, None], (L.Sy + ts * uy)[:, None], np.full((len(ts), 1), par.lgap), par)
            ok, cost, _ln = L.static(yc[:, None], q0[:, None], ya, ga, geo, par, par.cone_cos)
            for i in np.nonzero(ok[:, 0])[0]:
                x0 = q0[i] if L.sx > 0 else -(q0[i] + w)
                box = (x0, yc[i] - 0.5 * h, x0 + w, yc[i] + 0.5 * h)
                end = (L.sx * (q0[i] - ga[i, 0]), ya[i, 0])
                extra = (np.array([q0[i]]), np.array([cost[i, 0]]), np.array([box]), np.array([end]))
                res = self.try_at(L, float(yc[i]), pl, pend, extra=extra, only=True)[0]
                if res is not None:
                    return res
        return None

    def fan(self, recs: list, shore: str) -> tuple[list, dict]:
        """A chain of labels of one stack (boxes at most gap + 3 px apart) in which a leader bends is seated again, nearest the equator
        first, as a straight fan: every leader straight along its own first piece and ending at the middle of the box end (fan_seat),
        a label that is straight and within 0.35 of the half height of the middle staying where it is.  The first labels of the chain that
        can be seated again are kept when the layout still passes every test (each against all others), its labels stay inside the
        old extent of the labels, and the chain has more straight leaders than before."""
        par = self.par
        recs = list(recs)
        info = {"chains": 0, "accepted": 0, "labels": 0, "gain": 0, "added_px": 0.0}
        if not recs:
            return recs, info
        labs = [self.lab(r[0], shore, 1 if r[3] == "start" else -1) for r in recs]
        pend = _Pending(labs)
        for lb in labs:
            pend.done(lb)
        ext = (min(r[1][0] for r in recs), min(r[1][1] for r in recs), max(r[1][2] for r in recs), max(r[1][3] for r in recs))

        def straight(r) -> bool:
            return abs(_sturn(r[2])) < 3.0

        def fixed(k: int) -> bool:                           # straight, ending on the left or right end within 0.35 of the half height of its middle
            box, pts = recs[k][1], recs[k][2]
            ex, ey = pts[-1]
            if abs(_sturn(pts)) >= 1.0 or not (ex <= box[0] + 0.05 or ex >= box[2] - 0.05):
                return False
            return abs(ey - 0.5 * (box[1] + box[3])) <= 0.35 * 0.5 * (box[3] - box[1])

        def placed(skip) -> _Placed:
            pl = _Placed()
            for j, r in enumerate(recs):
                if j not in skip:
                    pl.add(r[1], labs[j].S, r[2][-1], labs[j])
            return pl

        groups: dict = {}
        for k, lb in enumerate(labs):
            groups.setdefault((lb.sx, lb.upper), []).append(k)
        for ks in groups.values():
            ks.sort(key=lambda k: 0.5 * (recs[k][1][1] + recs[k][1][3]) * labs[k].d)
            chains, cur = [], [ks[0]]
            for k in ks[1:]:
                a, b = recs[cur[-1]][1], recs[k][1]
                if max(b[1] - a[3], a[1] - b[3]) <= par.gap + 3.0:
                    cur.append(k)
                else:
                    chains.append(cur)
                    cur = [k]
            chains.append(cur)
            for chain in chains:
                if len(chain) < 2 or all(straight(recs[k]) for k in chain):
                    continue
                info["chains"] += 1
                before = [abs(_sturn(recs[k][2])) < 1.0 for k in chain]
                run = []
                for k in chain:                              # a leader steeper than 60 degrees ends the run (a pole block)
                    d = (labs[k].S[0] - recs[k][2][0][0], labs[k].S[1] - recs[k][2][0][1])
                    if abs(d[1]) > 0.866 * math.hypot(*d):
                        break
                    run.append(k)
                if len(run) < 2:
                    continue
                fix = {k for k in run if fixed(k)}
                best = None
                for seq in (run, run[::-1]):                  # from the equator outward, or from the outermost label inward
                    pl = placed(set(run) - fix)
                    new: dict = {}
                    for k in seq:
                        if k in fix:
                            continue
                        res = self.fan_seat(labs[k], recs[k], pl, pend)
                        if res is None:
                            break
                        new[k] = (res[0], [recs[k][0].anchor, labs[k].S, res[1]])
                        pl.add(res[0], labs[k].S, res[1], labs[k])
                    order = [k for k in seq if k in new]
                    while order:                            # the prefix of the run that stands against the labels that did not move
                        trial = list(recs)
                        for k in order:
                            trial[k] = (recs[k][0], new[k][0], new[k][1], recs[k][3], recs[k][4])
                        bad = None
                        for k in order:
                            others = _Placed()
                            for j, r in enumerate(trial):
                                if j != k:
                                    others.add(r[1], labs[j].S, r[2][-1], labs[j])
                            yy = 0.5 * (new[k][0][1] + new[k][0][3])
                            extra = (np.array([new[k][0][0] if labs[k].sx > 0 else -new[k][0][2]]), np.array([0.0]), np.array([new[k][0]]), np.array([new[k][1][-1]]))
                            if self.try_at(labs[k], yy, others, pend, extra=extra, only=True)[0] is None:
                                bad = k
                                break
                        if bad is None:
                            break
                        order = order[:order.index(bad)]
                    if not order:
                        continue
                    inside = all(new[k][0][0] >= ext[0] and new[k][0][1] >= ext[1] and new[k][0][2] <= ext[2] and new[k][0][3] <= ext[3] for k in order)
                    gain = sum(1 for k in order if abs(_sturn(new[k][1])) < 3.0) - sum(1 for k in order if straight(recs[k]))
                    after = [abs(_sturn(new[k][1] if k in order else recs[k][2])) < 1.0 for k in chain]
                    if not inside or gain <= 0 or (_broken(after) and not _broken(before)):     # more straight leaders, and no new bent leader between straight ones
                        continue
                    key = (_broken(after), -gain)
                    if best is None or key < best[0]:
                        best = (key, order, new, gain)
                if best is None:
                    continue
                _key, order, new, gain = best
                for k in order:
                    info["added_px"] += _plen(new[k][1]) - _plen(recs[k][2])
                    recs[k] = (recs[k][0], new[k][0], new[k][1], recs[k][3], recs[k][4])
                info["accepted"] += 1
                info["labels"] += len(order)
                info["gain"] += gain
        info["added_px"] = round(info["added_px"], 1)
        return recs, info

    # ---- prototype P: one cluster construction (R3 to R6)
    def chain_nodes(self, L: _Lab, rec, bg: _Placed, pend: _Pending):
        """Every seat of one label of a chain that passes the tests against the labels outside the chain: the rows within `pwin` px of its own
        row, every contact point, and the straight seats (the point on the feature's own ray) from the nearest one to `parc` px farther
        out.  A seat is straight (turn under 0.5 degrees) or bent (3 degrees or more, last piece at least `minlast`); its cost is run + tw x turn
        (+ `bend` for a bent one).  The label's present seat is one of them."""
        par, geo = self.par, self.geo
        yc0 = 0.5 * (rec[1][1] + rec[1][3])
        rows = yc0 + np.arange(-par.pwin, par.pwin + 0.5, 1.0)
        Q, YA, GA, _b, PID = L.cands_port(rows[:, None], geo, par)
        ok, cost, run = L.static(rows[:, None], Q, YA, GA, geo, par, par.cone_cos)
        ri, ci = np.nonzero(ok)
        yc, q0, ya, ga, pid, cst, rn = rows[ri], Q[ri, ci], YA[ri, ci], GA[ri, ci], PID[ri, ci], cost[ri, ci], run[ri, ci]
        ts = np.arange(0.5, 400.0, 0.5)
        parts = [(yc, q0, ya, ga, pid, cst, rn)]
        for j, (fq, fy) in enumerate(PORTS):               # the point on the feature's own ray at distance t from the stub
            ys_, qs_ = L.Sy + ts * L.ny - (fy - 0.5) * L.h, L.Sq + ts * L.nq - fq * L.w
            yas, gas = L.port_end(j, ys_, qs_, par)
            oks, costs, runs = L.static(ys_[:, None], qs_[:, None], yas[:, None], gas[:, None], geo, par, par.cone_cos)
            sel = oks[:, 0]
            parts.append((ys_[sel], qs_[sel], yas[sel], gas[sel], np.full(int(sel.sum()), j), costs[sel, 0], runs[sel, 0]))
        q0c = rec[1][0] if L.sx > 0 else -rec[1][2]            # the present box with its last piece re-aimed at every contact point its direction may use
        for j in range(len(PORTS)):
            yac, gac = L.port_end(j, yc0, q0c, par)
            okc, costc, runc = L.static(np.array([[yc0]]), np.array([[q0c]]), yac[None, None], gac[None, None], geo, par, par.cone_cos)
            if okc[0, 0]:
                parts.append((np.array([yc0]), np.array([q0c]), np.array([float(yac)]), np.array([float(gac)]), np.array([j]), costc[0], runc[0]))
        yc, q0, ya, ga, pid, cst, rn = (np.concatenate([p[i] for p in parts]) for i in range(7))
        x0 = q0 if L.sx > 0 else -(q0 + L.w)
        box = np.stack([x0, yc - 0.5 * L.h, x0 + L.w, yc + 0.5 * L.h], axis=1)
        att = np.stack([L.sx * (q0 - ga), ya], axis=1)
        turn = np.abs(_sturn_arr(L.it.anchor, L.S, att))
        straight = turn < 0.5
        keep = straight | (turn >= 3.0)
        t = rn + par.lgap                                   # distance of the point from the stub (straight seats)
        if (keep & straight).any():                         # R5: no farther than `parc` px along the ray than the nearest allowed seat
            keep &= ~straight | (t <= t[keep & straight].min() + par.parc)
        stub = np.broadcast_to(np.array(L.S, dtype=float), att.shape)
        idx = np.nonzero(keep)[0]
        if len(idx):
            lo, hi = min(L.Sy, box[idx, 1].min()) - par.lpad - 1.0, max(L.Sy, box[idx, 3].max()) + par.lpad + 1.0
            idx = idx[self.legal(L, box[idx], att[idx], stub[idx], bg, pend, lo, hi)]
        cost2 = (cst + np.where(straight, 0.0, par.bend) + par.pvert * np.isin(pid, (0, 1, 2, 3))     # a vertex is the last choice among the points that serve
                 + np.where(straight, 0.0, par.kink * np.clip(1.0 - turn / 8.0, 0.0, 1.0)) + par.tsq * turn ** 2)             # a slight bend costs more than a slide along the ray
        chosen: dict = {}                                   # bent seats: the cheapest of each (row, point), then the cheapest few of each row
        out = [i for i in idx if straight[i]]
        for i in sorted((i for i in idx if not straight[i]), key=lambda i: cost2[i]):
            chosen.setdefault((int(round(yc[i] * 2)), int(pid[i])), i)
        per_row: dict = {}
        for (r, _p), i in sorted(chosen.items(), key=lambda kv: cost2[kv[1]]):
            per_row.setdefault(r, []).append(i)
        out += [i for v in per_row.values() for i in v[:6]]
        out = sorted(out, key=lambda i: cost2[i])[:700]
        cur_pts = rec[2]
        cur_turn = abs(_sturn(cur_pts))
        cur_cost = _plen(cur_pts[1:]) + par.tw * cur_turn + (0.0 if cur_turn < 1.0 else par.bend) + (0.0 if _port_name(rec[1], cur_pts) else par.noport)
        sel = np.array(out, dtype=int)
        return {"box": np.vstack([box[sel], np.array([rec[1]], dtype=float)]), "att": np.vstack([att[sel], np.array([cur_pts[-1]], dtype=float)]),
                "yc": np.concatenate([yc[sel], [yc0]]), "cost": np.concatenate([cost2[sel], [cur_cost]]),
                "straight": np.concatenate([straight[sel], [cur_turn < 1.0]]), "rad": np.hypot(*np.vstack([att[sel], np.array([cur_pts[-1]], dtype=float)]).T),
                "present": np.arange(len(sel) + 1) == len(sel), "stub": np.array(L.S, dtype=float)}

    def chain_pair(self, a: dict, b: dict, d: float) -> np.ndarray:
        """(Na, Nb): seat i of label a and seat j of label b (b the next one along the stack) stand together: order along the stack, box gap,
        leader to box clearance both ways, no crossing."""
        par = self.par
        g = par.gap - 1e-6
        ab, bb = a["box"], b["box"]
        ok = ~((ab[:, None, 0] < bb[None, :, 2] + g) & (bb[None, :, 0] < ab[:, None, 2] + g)
               & (ab[:, None, 1] < bb[None, :, 3] + g) & (bb[None, :, 1] < ab[:, None, 3] + g))
        if d is not None:
            ok &= (b["yc"][None, :] - a["yc"][:, None]) * d > 0
        sa = np.broadcast_to(a["stub"], a["att"].shape)
        sb = np.broadcast_to(b["stub"], b["att"].shape)
        ok &= ~_seg_box_hit(sa, a["att"], bb, par.lpad)
        ok &= ~_seg_box_hit(sb, b["att"], ab, par.lpad).T
        ok &= ~_seg_cross_matrix(sa, a["att"], sb, b["att"])
        return ok

    @staticmethod
    def chain_dp(nodes: list, compat: list, compat2: list, pen: list, order: bool = True):
        """Cheapest choice of one seat per label along the stack with the order bent* straight* bent* (R3), where a label stands with
        its next neighbour (compat) and with the one after (compat2); a pair of neighbours also costs pen.  The state is the seats of two labels.  The index of the seat of
        every label and the cost of the choice, or None.  Without `order` any mix of straight and bent seats is allowed."""
        INF = 1e18
        n = len(nodes)
        if n == 1:
            c = nodes[0]["cost"]
            return ([int(c.argmin())], float(c.min())) if c.min() < INF / 2 else None
        # phase of the order after a seat: 0 bent before the straight run, 1 straight, 2 bent after it; PH[p][s]: the phase of a seat of
        # kind s (0 bent, 1 straight) after phase p, or -1
        PH = ((0, 1), (2, 1), (2, -1)) if order else ((0, 0), (0, 0), (0, 0))
        n0, n1 = nodes[0], nodes[1]
        s0, s1 = n0["straight"].astype(int), n1["straight"].astype(int)
        V = np.full((3, len(s0), len(s1)), INF)              # V[phase of the last label, seat of the one before, seat of the last]
        for a in range(2):
            for b in range(2):
                p0 = 0 if a == 0 else 1
                if PH[p0][b] >= 0:
                    V[PH[p0][b]] = np.where((s0[:, None] == a) & (s1[None, :] == b) & compat[0], n0["cost"][:, None] + n1["cost"][None, :] + pen[0], V[PH[p0][b]])
        back = []
        for k in range(2, n):
            nk = nodes[k]
            sk = nk["straight"].astype(int)
            ok = compat[k - 1][None, :, :] & compat2[k - 2][:, None, :]            # [a, b, c]: seats of labels k-2, k-1, k
            new = np.full((3, V.shape[2], len(sk)), INF)
            arg = np.zeros((3, V.shape[2], len(sk), 2), dtype=int)                 # best (phase before, seat of k-2) of each new state
            for p in range(3):
                m = np.where(ok, V[p][:, :, None], INF)
                am = m.argmin(axis=0)
                mv = np.take_along_axis(m, am[None], axis=0)[0]                    # [b, c]
                for sv in range(2):
                    q = PH[p][sv]
                    if q < 0:
                        continue
                    cand = np.where(sk[None, :] == sv, mv + nk["cost"][None, :] + pen[k - 1], INF)
                    better = cand < new[q]
                    new[q] = np.where(better, cand, new[q])
                    arg[q][better] = np.stack([np.full(int(better.sum()), p), am[better]], axis=1)
            back.append(arg)
            V = new
        best = V.min()
        if best >= INF / 2:
            return None
        p, b, c = np.unravel_index(int(V.argmin()), V.shape)
        pick = [0] * n
        pick[n - 1], pick[n - 2] = int(c), int(b)
        for k in range(n - 1, 1, -1):
            pp, a = back[k - 2][p, pick[k - 1], pick[k]]
            pick[k - 2], p = int(a), int(pp)
        return pick, float(best)

    def chain_search(self, nodes: list, compat: list, compat2: list, pen: list, d: float, budget: int, order: bool):
        """chain_dp over a chain whose labels may also collide with labels farther than the second one (the DP sees two).  The cheapest
        choice is taken; when two labels of it collide, the later label j is branched on (best first): it takes one of its `fan` cheapest
        seats with every seat of the other labels that collides with it struck out, or none of them.  The index of the seat of every
        label, or None when the budget of DP runs is spent."""
        INF = 1e18
        n = len(nodes)
        fan = 8

        def solve(banned):
            return self.chain_dp([{"cost": np.where(banned[k], INF, x["cost"]), "straight": x["straight"]} for k, x in enumerate(nodes)], compat, compat2, pen, order)

        def collides(j, seat):                                  # the seats of every other label that collide with seat `seat` of label j
            one = _select(nodes[j], [seat])
            return {m: ~(self.chain_pair(nodes[m], one, d)[:, 0] if m < j else self.chain_pair(one, nodes[m], d)[0]) for m in range(n) if m != j}

        banned0 = [np.zeros(len(x["cost"]), dtype=bool) for x in nodes]
        first = solve(banned0)
        heap, tick = ([(first[1], 0, banned0, first[0])] if first else []), 0
        while heap and budget > 0:
            _c, _t, banned, pick = heapq.heappop(heap)
            at = {key: np.array([nodes[k][key][i] for k, i in enumerate(pick)]) for key in ("box", "att", "yc")}
            at["stub"] = np.array([x["stub"] for x in nodes])
            bad = np.argwhere(~self.chain_pair(at, at, None) & ~np.eye(n, dtype=bool))
            if not len(bad):
                return pick
            j = int(bad[0][1])
            seats = [int(i) for i in np.argsort(nodes[j]["cost"], kind="stable") if not banned[j][i]][:fan]
            for seat in seats:
                nb = [b.copy() for b in banned]
                for m, c in collides(j, seat).items():
                    nb[m] |= c
                nb[j] = np.ones(len(banned[j]), dtype=bool)
                nb[j][seat] = False
                res = solve(nb)
                budget -= 1
                if res:
                    tick += 1
                    heapq.heappush(heap, (res[1], tick, nb, res[0]))
            nb = [b.copy() for b in banned]
            nb[j][seats] = True
            res = solve(nb)
            budget -= 1
            if res:
                tick += 1
                heapq.heappush(heap, (res[1], tick, nb, res[0]))
        return None

    def chain_pen(self, a: dict, b: dict) -> np.ndarray:
        """(Na, Nb): the cost of seat i of label a standing with seat j of label b (the next one along the stack).  R8: the distance of the two
        leader ends from the centre differs by more than `arcstep` px.  R4 (even rows): one of them is bent and the boxes are more than `rowgap` px apart."""
        par = self.par
        arc = par.arcw * np.maximum(np.abs(a["rad"][:, None] - b["rad"][None, :]) - par.arcstep, 0.0)
        gap = np.maximum(b["box"][None, :, 1] - a["box"][:, None, 3], a["box"][:, None, 1] - b["box"][None, :, 3])
        bent = ~a["straight"][:, None] | ~b["straight"][None, :]
        return arc + par.loosew * bent * np.clip(gap - par.rowgap, 0.0, par.loosemax)

    def chain_pick(self, nodes: list, d: float, budget: int, order: bool, arc: bool):
        """chain_search over the seats of every label thinned to the cheapest `seats` of each row of 3 px (bent and straight apart; the present seat).
        The index of the seat of every label, or None."""
        n = len(nodes)
        per = max(1, self.par.seats - (n > 20) - (n > 40))
        sub, maps = [], []
        for nd in nodes:
            seen: dict = {}
            keep = []
            for i in np.argsort(nd["cost"], kind="stable"):
                key = (bool(nd["straight"][i]), round(float(nd["box"][i][0])), round(float(nd["box"][i][1])), round(float(nd["att"][i][0])), round(float(nd["att"][i][1])))   # the same box and end (to 1 px) is one seat
                row = (bool(nd["straight"][i]), int(nd["yc"][i] // 3))
                if key in seen:
                    continue
                if nd["present"][i] or seen.get(row, 0) < per:
                    seen[key] = 1
                    seen[row] = seen.get(row, 0) + 1
                    keep.append(int(i))
            keep = sorted(keep)
            sub.append(_select(nd, keep))
            maps.append(keep)
        compat = [self.chain_pair(sub[i], sub[i + 1], d) for i in range(n - 1)]
        if arc:                                             # R8: neighbouring leader ends differ by at most `arcstep` px in distance from the centre
            compat = [c & (np.abs(sub[i]["rad"][:, None] - sub[i + 1]["rad"][None, :]) <= self.par.arcstep - 0.05) for i, c in enumerate(compat)]
        compat2 = [self.chain_pair(sub[i], sub[i + 2], d) for i in range(n - 2)]
        pen = [self.chain_pen(sub[i], sub[i + 1]) for i in range(n - 1)]
        pick = self.chain_search(sub, compat, compat2, pen, d, budget, order)
        return None if pick is None else [maps[k][i] for k, i in enumerate(pick)]

    def chain_solve(self, nodes: list, d: float):
        """The seat of every label of a chain, or None.  In the order of R3 and on the smooth arc of R8 if there is such a choice; else on the arc
        in any order (the chain reads mixed); else off the arc (the step costs `arcw` per px beyond `arcstep`)."""
        for order, arc in ((True, True), (False, True), (True, False), (False, False)):
            pick = self.chain_pick(nodes, d, 300, order, arc)
            if pick is not None:
                return pick
        return None

    # ---- prototype Q: a side cluster in arc coordinates
    def arc_geom(self, L: _Lab, rho, phi, pole: bool = False, lv: int = 0):
        """The seat of label L whose leader ends on the circle of radius rho at angle phi (radians, + towards +y) from the label's own ray.
        The last piece S -> end points at a contact point the arrival rule allows (middle of the short side within 45 degrees, the facing
        vertex 45 to 70, a long-side point beyond).  (box q0, box centre y, end q, end y, contact point index or -1)."""
        par = self.par
        c, sn = np.cos(phi), np.sin(phi)
        aq, ay = rho * (L.nq * c - L.ny * sn), rho * (L.nq * sn + L.ny * c)
        vq, vy = aq - L.Sq, ay - L.Sy
        ln = np.hypot(vq, vy)
        seen = ln >= par.nostraight
        dq = np.where(seen, vq / np.maximum(ln, 1e-9), L.nq)
        dy = np.where(seen, vy / np.maximum(ln, 1e-9), L.ny)
        cq, cy = aq + par.lgap * dq, ay + par.lgap * dy
        deg = np.degrees(np.arctan2(np.abs(dy), np.abs(dq)))
        if pole:                                           # Q2: beyond 70 degrees a long-side point (lv: 0, 1, 2 = the 0.25, 0.5, 0.75 point) needs no direction along the stack
            jl = np.where(dy > 0, 6, 9) + lv
            j = np.where(deg > 70.0, jl, np.where(dq <= 0.05, -1, np.where(deg <= 45.0, 4, np.where(dy > 0, 0, 2))))
        else:
            j = np.where(dq <= 0.05, -1, np.where(deg <= 45.0, 4, np.where(deg <= 70.0, np.where(dy > 0, 0, 2), np.where(dy > 0, 6, 9))))
        fq, fy = PORTS_ARR[np.maximum(j, 0), 0], PORTS_ARR[np.maximum(j, 0), 1]
        return cq - fq * L.w, cy - (fy - 0.5) * L.h, aq, ay, j

    def arc_static(self, L: _Lab, q0, yc, aq, ay, j):
        """Seats given by arc_geom against the tests that need no other label (static, and the arrival rule of port_end)."""
        q0, yc, aq, ay, j = (np.atleast_1d(np.asarray(v, dtype=float)) for v in (q0, yc, aq, ay, j))
        ok = np.zeros(len(q0), dtype=bool)
        for jj in np.unique(j[j >= 0]).astype(int):
            m = j == jj
            ya, ga = L.port_end(jj, yc[m], q0[m], self.par)
            okk, _c, _ln = L.static(yc[m][:, None], q0[m][:, None], ya[:, None], ga[:, None], self.geo, self.par, self.par.cone_cos)
            ok[m] = okk[:, 0] & np.isfinite(ya) & (np.abs(ya - ay[m]) < 0.05) & (np.abs(q0[m] - ga - aq[m]) < 0.05)
        return ok

    def arc_chain(self, chain: list, labs: list, recs: list, bg: _Placed, pend: _Pending):
        """Prototype Q: the labels of a side column chain (stack order) are seated in arc coordinates.  The leader end of each label is on a circle
        of radius rho_k concentric with the ring (neighbours differ by at most `arcslope` px).  The natural place of a label is where its radial
        ray meets that circle (a straight leader).  The box centres stand at least one row apart (ink boxes `gap` px apart) at the least total
        displacement (L1, pool adjacent violators: the median label of a pooled run stays straight, the others stand at the row distance on both
        sides).  A label that fails a hard test or turns more than `maxturn` degrees raises its radius by 2 px (smoothly) and the run is built
        again, up to `rhomax` px.  The new records {k: (box, points)}, or None."""
        par = self.par
        ks = sorted(chain, key=lambda k: recs[k][1][1] + recs[k][1][3])
        Ls = [labs[k] for k in ks]
        n = len(ks)
        self.q_why = "pole"
        if n < 2 or any(L.shore != "out" or L.nq < par.qnq for L in Ls):
            return None
        self.q_why = "nat"
        rs = [math.hypot(*L.S) for L in Ls]
        grid = np.arange(0.0, 80.5, 0.5)
        r0 = []
        for L, r in zip(Ls, rs):
            rho = r + grid
            q0, yc, aq, ay, j = self.arc_geom(L, rho, np.zeros_like(rho))
            ok = self.arc_static(L, q0, yc, aq, ay, j)
            if not ok.any():
                return None
            r0.append(float(rho[int(np.argmax(ok))]))
        r0 = np.array(r0)
        sep = np.array([0.5 * (Ls[i].h + Ls[i + 1].h) + par.gap + 1e-3 for i in range(n - 1)])
        cum = np.concatenate([[0.0], np.cumsum(sep)])
        self.q_dbg = {}
        phis = np.linspace(-0.35, 0.35, 1401)
        mid = len(phis) // 2
        r = r0.copy()
        for _it in range(120):
            env = _lipschitz(r, par.arcslope)
            nat, grids = np.empty(n), []
            for i, L in enumerate(Ls):
                q0, yc, _aq, _ay, j = self.arc_geom(L, env[i], phis)
                v = j >= 0
                grids.append((np.maximum.accumulate(yc[v]), phis[v]))
                nat[i] = yc[mid]
            yc_t = _pav_l1(nat - cum) + cum
            new, bad = {}, set()
            pl = bg.copy()
            seats = []
            for i, L in enumerate(Ls):
                phi = float(np.interp(yc_t[i], grids[i][0], grids[i][1]))
                q0, yc, aq, ay, j = (float(v[0]) for v in map(np.atleast_1d, self.arc_geom(L, env[i], phi)))
                seats.append((q0, yc, aq, ay, int(j)))
                x0 = q0 if L.sx > 0 else -(q0 + L.w)
                box = (x0, yc - 0.5 * L.h, x0 + L.w, yc + 0.5 * L.h)
                att = (L.sx * aq, ay)
                new[ks[i]] = (box, [L.it.anchor, L.S, att])
                pl.add(box, L.S, att, L)
                if j < 0 or not self.arc_static(L, q0, yc, aq, ay, j)[0]:
                    bad.add(i)
                    self.q_dbg[i] = "[static]"
                elif abs(_sturn_arr(L.it.anchor, L.S, np.array([att]))[0]) > par.maxturn:
                    bad.add(i)
                    self.q_dbg[i] = "[turn]"
            base = len(bg.labs())
            for i, L in enumerate(Ls):
                if i in bad:
                    continue
                box, pts = new[ks[i]]
                bx, at = np.array([box], dtype=float), np.array([pts[-1]], dtype=float)
                if not self.legal(L, bx, at, np.array([L.S], dtype=float), _View(pl, base + i), pend, min(L.Sy, bx[0, 1]) - par.lpad - 1.0, max(L.Sy, bx[0, 3]) + par.lpad + 1.0)[0]:
                    bad.add(i)
                    pb, _s0, _s1 = pl.arrays()
                    near = [pl.labs()[m].it.text[:12] + ("*" if m >= base else "") for m in range(len(pb)) if m != base + i and bx[0, 0] < pb[m, 2] + 12 and pb[m, 0] < bx[0, 2] + 12 and bx[0, 1] < pb[m, 3] + 12 and pb[m, 1] < bx[0, 3] + 12]
                    self.q_dbg[i] = "[legal " + "/".join(near) + "]"
            if not bad:
                return new
            for i in bad:
                r[i] = max(r[i], env[i]) + 2.0
            if (r - r0).max() > par.rhomax:
                self.q_why = "rho:" + ",".join(f"{Ls[i].it.text}{self.q_dbg.get(i, '')}" for i in sorted(bad))
                return None
        self.q_why = "iter"
        return None

    @staticmethod
    def _lv(L: _Lab) -> int:
        """Q2: the long-side point of a steep leader: the box overhangs the ring axis least (the 0.25 point right of the axis, the 0.75 point left of it)."""
        return 0 if L.Sq >= 0 else 2

    def arc_tables(self, L: _Lab, rho: float, lv: int):
        """Q2: what label L can reach from radius rho: the arc track (angle grid: box centre rows, valid = arrival point, static, turn <= maxturn)
        and the ray track (outward along its own ray: rows, valid = static).  Cached."""
        key = (id(L), round(float(rho), 2), lv)
        if key not in self.qcache:
            par = self.par
            phis = np.linspace(-0.35, 0.35, 701)
            q0, yc, aq, ay, j = self.arc_geom(L, rho, phis, True, self._lv(L))
            vq, vy = aq - L.Sq, ay - L.Sy
            ln = np.hypot(vq, vy)
            turn = np.where(ln < par.nostraight, 0.0, np.degrees(np.arccos(np.clip((vq * L.nq + vy * L.ny) / np.maximum(ln, 1e-9), -1.0, 1.0))))
            ok = self.arc_static(L, q0, yc, aq, ay, j) & (turn <= par.maxturn)
            ts = np.arange(0.0, 300.0, 0.5)
            q0r, ycr, aqr, ayr, jr = self.arc_geom(L, rho + ts, np.zeros_like(ts), True, self._lv(L))
            self.qcache[key] = (phis, yc, ok, ts, ycr, self.arc_static(L, q0r, ycr, aqr, ayr, jr))
        return self.qcache[key]

    def arc_chain2(self, chain: list, labs: list, recs: list, bg: _Placed, pend: _Pending):
        """Prototype Q round 2: as arc_chain for any run of a column (pole labels too).  The row (box centre y) is the stack variable.  From its natural
        place a label changes its row on the arc track (bent leader, turn <= maxturn) or, when its ray is steeper than 45 degrees and the row
        moves toward the pole, on the ray track (outward along its own ray, straight).  One bounded L1 fit over the run (pool adjacent violators,
        clipped by the monotone envelope of the per-label reach); a label whose reach is too short or that fails a test raises its radius (smooth,
        up to `qrhomax` px).  Returns {k: (box, points)} or None (self.q_why)."""
        par = self.par
        ks = sorted(chain, key=lambda k: recs[k][1][1] + recs[k][1][3])
        Ls = [labs[k] for k in ks]
        n = len(ks)
        self.q_why = "short"
        if n < 2 or any(L.shore != "out" for L in Ls):
            return None
        self.q_why = "nat"
        steep = [abs(L.ny) > 0.7071 for L in Ls]
        pw = [-1.0 if L.Sy < 0 else 1.0 for L in Ls]
        rs = [math.hypot(*L.S) for L in Ls]
        grid = np.arange(0.0, 80.5, 0.5)
        sep = np.array([0.5 * (Ls[i].h + Ls[i + 1].h) + par.gap + 1e-3 for i in range(n - 1)])
        cum = np.concatenate([[0.0], np.cumsum(sep)])
        mid = 350
        lvs = (0,)
        why = "iter"
        for lv in lvs:
            r0 = []
            for L, r in zip(Ls, rs):
                rho = r + grid
                q0, yc, aq, ay, j = self.arc_geom(L, rho, np.zeros_like(rho), True, self._lv(L))
                ok = self.arc_static(L, q0, yc, aq, ay, j)
                r0.append(float(rho[int(np.argmax(ok))]) if ok.any() else r + 10.0)   # a straight leader at its natural place may be refused (a wide label on a steep ray): the arc track decides
            if r0 is None:
                continue
            r0 = np.array(r0)
            r = r0.copy()
            self.q_dbg = {}
            done = False
            for _it in range(80):
                env = _lipschitz(r, par.arcslope)
                tabs = [self.arc_tables(L, env[i], self._lv(L)) for i, L in enumerate(Ls)]
                nat = np.array([t[1][mid] for t in tabs])
                lo, hi = np.empty(n), np.empty(n)
                bad = set()
                for i, (phis, yc, ok, ts, ycr, okr) in enumerate(tabs):
                    if not ok.any() and not okr.any():
                        bad.add(i)
                        lo[i], hi[i] = nat[i], nat[i]
                        continue
                    lo[i], hi[i] = (yc[ok].min(), yc[ok].max()) if ok.any() else (nat[i], nat[i])
                    if steep[i] and okr.any():                                # toward the pole the ray track reaches as far as it is valid
                        if pw[i] < 0:
                            lo[i] = ycr[okr].min()
                        else:
                            hi[i] = ycr[okr].max()
                    lo[i], hi[i] = min(lo[i], nat[i]), max(hi[i], nat[i])
                lo_e = np.maximum.accumulate(lo - cum)
                hi_e = np.minimum.accumulate((hi - cum)[::-1])[::-1]
                if (lo_e > hi_e + 1e-9).any():                      # the reach of two labels does not meet the rows between them
                    for k in np.nonzero(lo_e > hi_e + 1e-9)[0]:
                        bad.add(int(np.argmax((lo - cum)[:k + 1])))
                        bad.add(int(k + np.argmin((hi - cum)[k:])))
                    z = np.clip(_pav_l1(nat - cum), np.minimum(lo_e, hi_e), np.maximum(lo_e, hi_e))
                else:
                    z = np.clip(_pav_l1(nat - cum), lo_e, hi_e)
                yc_t = z + cum
                new = {}
                pl = bg.copy()
                for i, L in enumerate(Ls):
                    phis, yc, ok, ts, ycr, okr = tabs[i]
                    on_ray = steep[i] and okr.any() and (yc_t[i] - nat[i]) * pw[i] > 1e-6
                    if on_ray:
                        o = np.argsort(ycr[okr])
                        t = float(np.interp(yc_t[i], ycr[okr][o], ts[okr][o]))
                        seat = tuple(float(v[0]) for v in map(np.atleast_1d, self.arc_geom(L, env[i] + t, 0.0, True, self._lv(L))))
                    else:
                        phi = float(np.interp(yc_t[i], np.maximum.accumulate(yc[ok]), phis[ok])) if ok.any() else 0.0
                        seat = tuple(float(v[0]) for v in map(np.atleast_1d, self.arc_geom(L, env[i], phi, True, self._lv(L))))
                    q0, ycc, aq, ay, j = seat
                    x0 = q0 if L.sx > 0 else -(q0 + L.w)
                    box = (x0, ycc - 0.5 * L.h, x0 + L.w, ycc + 0.5 * L.h)
                    att = (L.sx * aq, ay)
                    new[ks[i]] = (box, [L.it.anchor, L.S, att])
                    pl.add(box, L.S, att, L)
                    if j < 0 or not self.arc_static(L, q0, ycc, aq, ay, j)[0]:
                        bad.add(i)
                        self.q_dbg[i] = "[static]"
                    elif not on_ray and abs(_sturn_arr(L.it.anchor, L.S, np.array([att]))[0]) > par.maxturn:
                        bad.add(i)
                        self.q_dbg[i] = "[turn]"
                base = len(bg.labs())
                for i, L in enumerate(Ls):
                    if i in bad:
                        continue
                    box, pts = new[ks[i]]
                    bx, at = np.array([box], dtype=float), np.array([pts[-1]], dtype=float)
                    if not self.legal(L, bx, at, np.array([L.S], dtype=float), _View(pl, base + i), pend, min(L.Sy, bx[0, 1]) - par.lpad - 1.0, max(L.Sy, bx[0, 3]) + par.lpad + 1.0)[0]:
                        bad.add(i)
                        self.q_dbg[i] = "[legal]"
                if not bad:
                    return new
                for i in bad:
                    r[i] = max(r[i], env[i]) + 3.0
                if (r - r0).max() > par.qrhomax:
                    why = "rho:" + ",".join(f"{Ls[i].it.text}{self.q_dbg.get(i, '[reach]')}" for i in sorted(bad))
                    break
        self.q_why = why
        return None

    def cluster(self, recs: list, shore: str) -> tuple[list, dict]:
        """R3 to R6 as one construction.  A chain of one stack (boxes at most gap + 3 px apart; a lone label is a chain of one) is seated again
        as a whole: every label takes one of its seats (chain_nodes) so that the sum of the leader costs is least and the leaders read
        bent* straight* bent* along the stack (chain_dp).  The result is kept when every label passes all tests against every other (try-at
        rules); otherwise the chain keeps its seats."""
        par = self.par
        recs = list(recs)
        info = {"chains": 0, "built": 0, "kept_nodes": 0, "kept_dp": 0, "kept_check": 0, "changed": 0}
        if not recs:
            return recs, info
        labs = [self.lab(r[0], shore, 1 if r[3] == "start" else -1) for r in recs]
        pend = _Pending(labs)
        for lb in labs:
            pend.done(lb)
        groups: dict = {}
        for k, lb in enumerate(labs):
            groups.setdefault((lb.sx, lb.upper), []).append(k)
        solved: set = set()
        qdone: set = set()
        if par.arcfirst and shore == "out":                      # Q: side clusters in arc coordinates, one column (both halves) at a time
            cols: dict = {}
            for k, lb in enumerate(labs):
                cols.setdefault(lb.sx, []).append(k)
            for ks in cols.values():
                ks.sort(key=lambda k: recs[k][1][1])
                chains, cur = [], [ks[0]]
                for k in ks[1:]:
                    if recs[k][1][1] - recs[cur[-1]][1][3] <= par.qgap:
                        cur.append(k)
                    else:
                        chains.append(cur)
                        cur = [k]
                chains.append(cur)
                for chain in chains:
                    if len(chain) < 2:
                        continue
                    bg = _Placed()
                    for j, r in enumerate(recs):
                        if j not in chain:
                            bg.add(r[1], labs[j].S, r[2][-1], labs[j])
                    new = self.arc_chain(chain, labs, recs, bg, pend)
                    if new is None and par.arcfirst == 2:           # Q2: the runs round 1 could not build (pole labels, long runs)
                        new = self.arc_chain2(chain, labs, recs, bg, pend)
                    if new is None:
                        info["q_kept"] = info.get("q_kept", 0) + 1
                        info.setdefault("q_kept_names", []).append([self.q_why] + [labs[c].it.text for c in chain])
                        continue
                    info["q_built"] = info.get("q_built", 0) + 1
                    qdone.update(chain)
                    for k, (box, pts) in new.items():
                        recs[k] = (recs[k][0], box, pts, recs[k][3], recs[k][4])
        for rnd in range(3):                                     # round 0: every chain; later rounds: the chains that read mixed or off the arc and are new
            for ks in groups.values():
                ks.sort(key=lambda k: 0.5 * (recs[k][1][1] + recs[k][1][3]) * labs[k].d)
                chains, cur = [], [ks[0]]
                for k in ks[1:]:                                 # the cluster of scripts/fan_order.py: boxes at most 6.1 px apart
                    a, b = recs[cur[-1]][1], recs[k][1]
                    if max(b[1] - a[3], a[1] - b[3]) <= 6.1:
                        cur.append(k)
                    else:
                        chains.append(cur)
                        cur = [k]
                chains.append(cur)
                for chain in chains:
                    if qdone.intersection(chain):
                        continue
                    if rnd and (len(chain) < 2 or tuple(chain) in solved or not (_broken([abs(_sturn(recs[k][2])) < 1.0 for k in chain]) or _off_arc([math.hypot(*recs[k][2][-1]) for k in chain]))):
                        continue
                    solved.add(tuple(chain))
                    info["chains"] += 1
                    bg = _Placed()
                    for j, r in enumerate(recs):
                        if j not in chain:
                            bg.add(r[1], labs[j].S, r[2][-1], labs[j])
                    nodes = [self.chain_nodes(labs[k], recs[k], bg, pend) for k in chain]
                    d0 = labs[chain[0]].d
                    pick = self.chain_solve(nodes, d0)
                    if pick is None:
                        info["kept_dp"] += 1
                        continue
                    new = {}
                    for k, nd, i in zip(chain, nodes, pick):
                        box = tuple(float(v) for v in nd["box"][i])
                        new[k] = (box, [recs[k][0].anchor, labs[k].S, (float(nd["att"][i][0]), float(nd["att"][i][1]))])
                    ok_chain = True
                    for k in chain:                                  # every label against everything else, with the whole new chain in place
                        others = _Placed()
                        for j, r in enumerate(recs):
                            if j != k:
                                b, pts = (new[j][0], new[j][1]) if j in new else (r[1], r[2])
                                others.add(b, labs[j].S, pts[-1], labs[j])
                        bx = np.array([new[k][0]], dtype=float)
                        at = np.array([new[k][1][-1]], dtype=float)
                        if not self.legal(labs[k], bx, at, np.array([labs[k].S], dtype=float), others, pend, min(labs[k].Sy, bx[0, 1]) - par.lpad - 1.0, max(labs[k].Sy, bx[0, 3]) + par.lpad + 1.0)[0]:
                            ok_chain = False
                            info["kept_check"] += 1
                            info.setdefault("kept_names", []).append([labs[c].it.text for c in chain] + [labs[k].it.text])
                            break
                    if not ok_chain:
                        continue
                    info["built"] += 1
                    for k in chain:
                        if new[k][0] != tuple(recs[k][1]) or new[k][1][-1] != recs[k][2][-1]:
                            info["changed"] += 1
                        recs[k] = (recs[k][0], new[k][0], new[k][1], recs[k][3], recs[k][4])
        return recs, info

    # ---- slide toward the ring (LP-D14)
    def slide(self, recs: list, shore: str) -> tuple[list, dict]:
        """A leader that is straight, or leans, may end wherever its line first meets the grown box: every label of a finished
        layout is seated again, alone, at its own height and at heights up to 26 px toward its equator (along its own line
        toward the ring), with the same clearances as the seating (try_at); the new spot is kept when its leader is shorter
        (cost = run + tw x turn, with a short run beside a corner counted as the turn it is, at least 0.5 px less), its turn
        differs from the old one by at most 1 degree (a straight leader stays straight, the lean of a fan stays), and the
        whole leader is at most 1.5 px longer."""
        par = self.par
        recs = list(recs)
        info = {"tried": 0, "moved": 0, "saved_px": 0.0}
        if not recs:
            return recs, info
        labs = [self.lab(r[0], shore, 1 if r[3] == "start" else -1) for r in recs]
        pend = _Pending(labs)
        for lb in labs:
            pend.done(lb)

        def build() -> _Placed:
            pl = _Placed()
            for k, r in enumerate(recs):
                pl.add(r[1], labs[k].S, r[2][-1], labs[k])
            return pl

        def lcost(L: _Lab, att) -> float:
            vx, vy = att[0] - L.S[0], att[1] - L.S[1]
            ln = math.hypot(vx, vy)
            if ln < par.nostraight:
                return ln
            cosv = max(-1.0, min(1.0, (vx * L.nq * L.sx + vy * L.ny) / ln))
            return ln + par.tw * math.degrees(math.acos(cosv))

        def shifted(L: _Lab, end, yy: float, x_box: float):
            """Spots with the box moved 1..26 px toward the ring at height yy, the leader keeping its direction."""
            dirs = []
            for a, b in ((L.S, end),):                     # the leader keeps its direction (turn and lean)
                n = math.dist(a, b)
                if n >= 0.5:
                    dirs.append((L.sx * (b[0] - a[0]) / n, (b[1] - a[1]) / n))
            lg = par.lgap
            q0s, yas, gas = [], [], []
            for dq, dyv in dirs:
                for dx in range(1, 27):
                    q0 = x_box - dx
                    t0, t1 = -math.inf, math.inf
                    for p0, d, a, b in ((L.Sq, dq, q0 - lg, q0 + L.w + lg), (L.Sy, dyv, yy - 0.5 * L.h - lg, yy + 0.5 * L.h + lg)):
                        if abs(d) < 1e-9:
                            if not a <= p0 <= b:
                                t0, t1 = 1.0, 0.0
                            continue
                        ta, tb = (a - p0) / d, (b - p0) / d
                        t0, t1 = max(t0, min(ta, tb)), min(t1, max(ta, tb))
                    if t0 > t1 or t1 < 0:
                        continue
                    t = max(t0, 0.0)
                    q0s.append(q0)
                    yas.append(L.Sy + t * dyv)
                    gas.append(q0 - (L.Sq + t * dq))
            if not q0s:
                return None
            q0a, ya, ga = (np.array(v, dtype=float)[None, :] for v in (q0s, yas, gas))
            if par.face:
                ya, ga = L.snap(np.array([[yy]]), q0a, ya, ga, par)
            ok, cost, _run = L.static(np.array([[yy]]), q0a, ya, ga, self.geo, par, par.cone_cos)
            idx = np.nonzero(ok[0])[0]
            if not len(idx):
                return None
            q, c = q0a[0][idx], cost[0][idx]
            x0 = q if L.sx > 0 else -(q + L.w)
            boxes = np.stack([x0, np.full(len(q), yy - 0.5 * L.h), x0 + L.w, np.full(len(q), yy + 0.5 * L.h)], axis=1)
            att = np.stack([L.sx * (q - ga[0][idx]), ya[0][idx]], axis=1)
            return q, c, boxes, att

        for lb in self.cache.values():                     # in this pass a short run beside a box corner is a turn like any other
            lb.memo.clear()
        par.nostraight = 0.5
        pl = build()
        groups: dict = {}
        for k, lb in enumerate(labs):
            groups.setdefault((lb.sx, lb.upper), []).append(k)
        for _round in range(2):
            moved = 0
            for ks in groups.values():
                ks.sort(key=lambda k: 0.5 * (recs[k][1][1] + recs[k][1][3]) * labs[k].d)
                for k in ks:
                    L, (it, box, pts, align, tag) = labs[k], recs[k]
                    cur = lcost(L, pts[-1])
                    if math.dist(L.S, pts[-1]) < 3.5 and cur - math.dist(L.S, pts[-1]) < 1.0:
                        continue                             # a stub on its own line: nothing to gain
                    s_old = _sturn(pts)
                    yc = 0.5 * (box[1] + box[3])
                    view = _View(pl, k)
                    best = None
                    for step in (0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0, 26.0):
                        info["tried"] += 1
                        extra = shifted(L, pts[-1], yc - L.d * step, box[0] if L.sx > 0 else -box[2]) if par.shift and step <= 2.0 else None   # at most 2 px up: a gap of the stack stays under 6.1 px
                        for res in (self.try_at(L, yc - L.d * step, view, pend)[0], self.try_at(L, yc - L.d * step, view, pend, extra=extra, only=True)[0] if extra else None):
                            if res is None or res[3] >= cur - 0.5 or (best is not None and res[3] >= best[0]):
                                continue
                            npts = [it.anchor, L.S, res[1]]
                            if abs(_sturn(npts) - s_old) <= 1.0 and _plen(npts) <= _plen(pts) + 1.5:
                                best = (res[3], res[0], npts)
                    if best is not None:
                        info["saved_px"] += _plen(pts) - _plen(best[2])
                        recs[k] = (it, best[1], best[2], align, tag)
                        pl = build()
                        moved += 1
                        info["moved"] += 1
            if not moved:
                break
        par.nostraight = 3.0
        for lb in self.cache.values():
            lb.memo.clear()
        info["saved_px"] = round(info["saved_px"], 1)
        return recs, info

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

def _ink_item(it: Item, default_size: float) -> Item:
    """Prototype P (R1, LP-D25): the label's box is the ink box of its text."""
    x0, top, x1, bot = _ink_rel(it.text, float(it.src.get("font_size", default_size)), "start")
    return dataclasses.replace(it, w=x1 - x0, h=bot - top)


def _port_name(box, pts) -> str | None:
    """The contact point (name) the last piece of the leader points at, or None."""
    kept = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= 3.0]
    a, e = (kept[-2], pts[-1]) if len(kept) >= 2 else (pts[0], pts[-1])
    n = math.dist(a, e)
    if n < 1e-9:
        return None
    d = ((e[0] - a[0]) / n, (e[1] - a[1]) / n)
    best = None
    for (fx, fy), name in PORT_NAMES.items():             # fx: fraction of the width from the left (real coordinates)
        qx, qy = box[0] + fx * (box[2] - box[0]), box[1] + fy * (box[3] - box[1])
        ahead = (qx - e[0]) * d[0] + (qy - e[1]) * d[1]
        off = abs((qx - e[0]) * d[1] - (qy - e[1]) * d[0])
        if -0.5 <= ahead <= 5.0 and off <= 0.75 and (best is None or off < best[1]):
            best = (name, off)
    return best[0] if best else None


def _emit_ink(it: Item, box, pts, align: str, tag: str, default_size: float) -> dict:
    """emit() for a label whose box is the ink box: the text origin is shifted so that the ink fills the box."""
    label = emit(it, box, pts, align=align, tag=tag)
    x0, top, _x1, _bot = _ink_rel(it.text, float(it.src.get("font_size", default_size)), {"start": "start", "end": "end", "middle": "middle"}[align])
    label.update(text_x=box[0] - x0, text_y=box[1] - top, start_x=box[0] - x0, start_y=box[1] - top, port=_port_name(box, pts) if pts else None)
    return label


def place(items: list[Item], ctx: Ctx) -> list[dict]:
    par = _Par(ctx)
    inner_req = list(getattr(ctx, "inner_items", []) or [])
    if par.port:
        items = [_ink_item(i, par.font) for i in items]
        inner_req = [_ink_item(i, par.font) for i in inner_req]
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
    sym_info = {}
    if par.sym:
        out_recs, sym_info["out"] = eng.symmetrize(out_recs, "out")
        in_recs, sym_info["in"] = eng.symmetrize(in_recs, "in")
    if par.face and par.strt:
        out_recs, sym_info["strt_out"] = eng.straighten(out_recs, "out")
        in_recs, sym_info["strt_in"] = eng.straighten(in_recs, "in")
    if par.slide:
        out_recs, sym_info["slide_out"] = eng.slide(out_recs, "out")
        in_recs, sym_info["slide_in"] = eng.slide(in_recs, "in")
    if par.face and par.fan:
        out_recs, sym_info["fan_out"] = eng.fan(out_recs, "out")
        in_recs, sym_info["fan_in"] = eng.fan(in_recs, "in")
    if par.port and par.cluster:
        out_recs, sym_info["cluster_out"] = eng.cluster(out_recs, "out")
        in_recs, sym_info["cluster_in"] = eng.cluster(in_recs, "in")
    recs = out_recs + in_recs
    hide = bool(int(ctx.opt("hide", 0)))
    if par.port:
        result = [_emit_ink(it, box, None if hide else pts, align, tag, par.font) for it, box, pts, align, tag in recs]
    else:
        result = [emit(it, box, None if hide else pts, align=align, tag=tag) for it, box, pts, align, tag in recs]
    if hide:
        for label in result:
            label["text_fill"] = "none"
    ctx.notes.update(
        inside=len(in_recs), outside=len(out_recs), asked_inside=len(inner_req), overflow=[i.text for i in overflow], moves=moves,
        total=round(sum(_plen(r[2]) for r in recs), 1), turn_sum=round(sum(_turn(r[2]) for r in recs), 1), R=round(geo.R, 2), B=round(geo.B, 2), Bb=round(geo.Bb, 2),
        A_in=round(geo.A_in, 2), zone_ok=geo.zone_ok, objects=geo.objects, raster=geo.has, tick_font=geo.tick_font,
        cuts=[round(math.degrees(c), 1) for c in cuts] if cuts else None, stats=stats, sym=sym_info, check=_check(recs),
    )
    return result
