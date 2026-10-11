"""Study prototype A: keep today's look, but never overprint (A1 font, A2 band, A3 drop, A combined).

The wrapper in labels/circular.py hands us ``render_today(font_size, band_scale, drop)``, which runs today's
unchanged placement.  We search over its three knobs until a layout "fits", where fits means (see ``evaluate``):
no two label boxes overlap, no leader runs through another label's text, no leader runs over the feature ring,
and no label sits more than ``drift`` degrees (default 30) away from its feature.  The search is a bisection on a
monotone path of settings; the final choice is always a layout that was rendered and verified.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time

from .. import circular as lc

DRIFT_DEG = 30.0
FONT_STEP = 0.5
FONT_FLOOR = 6.0
BAND_STEP = 0.05
BAND_CAP = 2.0
# Variant A (shipped policy): band and font move together along the measured frontier up to this cap, then font alone, then drop.
A_BAND_CAP = 1.5
A_TRADE = 0.035  # relative font reduction per band step of 0.05 along the A path


# ------------------------------------------------------------------ fit criterion

def _box(label, total_length):
    if label.get("aabb_local") is not None:
        return tuple(float(v) for v in label["aabb_local"])
    return tuple(float(v) for v in lc._label_bbox_rect(label, total_length, 0.0))


def _leader(label):
    if label.get("leader_points") is not None:
        return [tuple(p) for p in label["leader_points"]]
    return [
        (float(label.get("feature_anchor_x", label["feature_middle_x"])), float(label.get("feature_anchor_y", label["feature_middle_y"]))),
        (float(label["middle_x"]), float(label["middle_y"])),
        (float(label.get("leader_start_x", label["start_x"])), float(label.get("leader_start_y", label["start_y"]))),
    ]


def _orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _cross(p, q, r, s):
    d1, d2, d3, d4 = _orient(r, s, p), _orient(r, s, q), _orient(p, q, r), _orient(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-6 for d in (d1, d2, d3, d4))


def _seg_box(p, q, box, pad=-0.5):
    x0, y0, x1, y1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for d, lo, hi, st in ((dx, x0, x1, p[0]), (dy, y0, y1, p[1])):
        if abs(d) < 1e-12:
            if st <= lo or st >= hi:
                return False
            continue
        ta, tb = (lo - st) / d, (hi - st) / d
        if ta > tb:
            ta, tb = tb, ta
        t0, t1 = max(t0, ta), min(t1, tb)
        if t0 >= t1:
            return False
    return True


def _seg_min_radius(p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    l2 = dx * dx + dy * dy
    t = 0.0 if l2 == 0 else max(0.0, min(1.0, -(p[0] * dx + p[1] * dy) / l2))
    return math.hypot(p[0] + t * dx, p[1] + t * dy)


def evaluate(labels, total_length, drift_deg=DRIFT_DEG):
    """Return the layout facts that decide "fits" (and a few for the report)."""
    ext = [lab for lab in labels if not lab.get("is_embedded")]
    n = len(ext)
    boxes = [_box(lab, total_length) for lab in ext]
    leaders = [_leader(lab) for lab in ext]
    order = sorted(range(n), key=lambda i: boxes[i][0])
    overlaps = 0
    for a in range(n):
        i = order[a]
        for b in range(a + 1, n):
            j = order[b]
            if boxes[j][0] >= boxes[i][2] - 0.5:
                break
            if min(boxes[i][3], boxes[j][3]) - max(boxes[i][1], boxes[j][1]) > 0.5:
                overlaps += 1
    segs = [list(zip(l, l[1:])) for l in leaders]
    through = 0
    for i in range(n):
        for j in range(n):
            if i != j and any(_seg_box(p, q, boxes[j]) for p, q in segs[i]):
                through += 1
    crossings = 0
    for i in range(n):
        for j in range(i + 1, n):
            if any(_cross(p, q, r, s) for p, q in segs[i] for r, s in segs[j]):
                crossings += 1
    feat_r = max((float(lab.get("feature_outer_radius_px", 0.0)) for lab in ext), default=0.0)
    over_features = sum(1 for s in segs if any(_seg_min_radius(p, q) < feat_r - 1.0 for p, q in s[1:]))
    drifts = []
    for lab in ext:
        a = math.atan2(float(lab["feature_middle_x"]), -float(lab["feature_middle_y"]))
        end = _leader(lab)[-1]
        b = math.atan2(end[0], -end[1])
        drifts.append(abs((b - a + math.pi) % (2 * math.pi) - math.pi) * 180.0 / math.pi)
    lengths = [sum(math.dist(p, q) for p, q in s) for s in segs if s]
    max_drift = max(drifts, default=0.0)
    fit = overlaps == 0 and through == 0 and over_features == 0 and max_drift <= drift_deg
    return {
        "fit": fit, "n": n, "overlap_pairs": overlaps, "through_text": through, "over_features": over_features,
        "crossings": crossings, "max_drift_deg": round(max_drift, 1),
        "drift_gt": sum(1 for d in drifts if d > drift_deg), "leader_max": round(max(lengths, default=0.0), 1),
        "leader_mean": round(sum(lengths) / len(lengths), 1) if lengths else 0.0,
        "bbox": [round(min(b[0] for b in boxes), 1), round(min(b[1] for b in boxes), 1),
                 round(max(b[2] for b in boxes), 1), round(max(b[3] for b in boxes), 1)] if boxes else None,
    }


# ------------------------------------------------------------------ feature span (only importance measure)

def _feature_spans(render_today):
    """Longest location segment (bp) of every labelled feature, read from the feature dict that the wrapper closed over."""
    spans: dict[str, float] = {}
    try:
        cells = dict(zip(render_today.__code__.co_freevars, (c.cell_contents for c in render_today.__closure__)))
        feature_dict = cells["args"][0]
        for key, feature in feature_dict.items():
            parts = [abs(int(p.end) - int(p.start)) for p in getattr(feature, "location", [])]
            spans[str(key)] = float(max(parts)) if parts else 0.0
    except Exception:  # noqa: BLE001 - prototype: fall back to label extent below
        return {}
    return spans


# ------------------------------------------------------------------ search

def _with_working_band(render_today):
    """Make ``band_scale`` effective.  When the caller passes ``outer_arena`` (the usual CLI path), circular.py rescales
    the label ellipse to the arena right after ``band_scale`` is applied, which cancels it; so scale the arena too."""
    try:
        cells = dict(zip(render_today.__code__.co_freevars, (c.cell_contents for c in render_today.__closure__)))
        kw = cells["kwargs"]
    except Exception:  # noqa: BLE001
        return render_today
    arena = kw.get("outer_arena")
    if arena is None:
        return render_today

    def render(font_size=None, band_scale=1.0, drop=frozenset()):
        kw["outer_arena"] = (arena[0], arena[1] * band_scale)
        try:
            return render_today(font_size=font_size, band_scale=1.0, drop=drop)
        finally:
            kw["outer_arena"] = arena

    return render


class Prober:
    def __init__(self, render_today, total_length, drift):
        self.render, self.total, self.drift = _with_working_band(render_today), total_length, drift
        self.cache: dict = {}
        self.log: list[dict] = []
        self.seconds = 0.0

    def __call__(self, font, band, drop=frozenset()):
        key = (font, round(band, 4), tuple(sorted(drop)))
        if key not in self.cache:
            t0 = time.time()
            labels = self.render(font_size=font, band_scale=band, drop=frozenset(drop))
            t1 = time.time()
            ev = evaluate(labels, self.total, self.drift)
            dt = t1 - t0
            self.seconds += dt + (time.time() - t1)
            self.cache[key] = (labels, ev)
            self.log.append({"font": font, "band": round(band, 4), "drop": len(drop), "sec": round(dt, 1),
                             **{k: ev[k] for k in ("fit", "overlap_pairs", "through_text", "over_features", "max_drift_deg", "crossings")}})
            print(f"[A] probe font={font} band={band:.2f} drop={len(drop)} fit={ev['fit']} ov={ev['overlap_pairs']} "
                  f"tt={ev['through_text']} of={ev['over_features']} drift={ev['max_drift_deg']} bbox={ev['bbox']} ({dt:.1f}s)",
                  file=sys.stderr, flush=True)
        return self.cache[key]


def first_fit(path, probe):
    """path[0] is known not to fit; settings get more drastic along the path.  Return (index, result) of the first
    fitting setting found by bisection (assumes monotone), or (None, last result) when even the end does not fit."""
    last = len(path) - 1
    if last < 1:
        return None, None
    res = probe(*path[last])
    if not res[1]["fit"]:
        return None, res
    lo, hi, best = 0, last, res
    while hi - lo > 1:
        mid = (lo + hi) // 2
        r = probe(*path[mid])
        if r[1]["fit"]:
            hi, best = mid, r
        else:
            lo = mid
    return hi, best


def run(design, render_today, *, total_length, options):
    t_start = time.time()
    drift = float(options.get("drift", DRIFT_DEG))
    floor = float(options.get("floor", FONT_FLOOR))
    band_cap = float(options.get("band_cap", BAND_CAP))
    a_cap = float(options.get("a_band_cap", A_BAND_CAP))
    probe = Prober(render_today, total_length, drift)

    base_labels = probe.render()
    base_ev = evaluate(base_labels, total_length, drift)
    ext0 = [lab for lab in base_labels if not lab.get("is_embedded")]
    font0 = float(ext0[0]["font_size"]) if ext0 else (float(base_labels[0]["font_size"]) if base_labels else 0.0)
    probe.cache[(font0, 1.0, ())] = (base_labels, base_ev)
    probe.log.append({"font": font0, "band": 1.0, "drop": 0, "sec": 0.0, "fit": base_ev["fit"]})

    sweep = options.get("mode") == "sweep"
    summary = {"requested": len(ext0), "placed": len(ext0), "dropped": [], "font_size": font0, "band_scale": 1.0,
               "attempts": 0, "seconds": 0.0, "today_fits": base_ev["fit"], "drift_limit_deg": drift,
               "initial": {k: base_ev[k] for k in ("overlap_pairs", "through_text", "over_features", "max_drift_deg", "crossings")}}

    if sweep:
        _sweep(options, probe, font0, floor, band_cap, ext0, render_today)
    if base_ev["fit"] or not ext0 or sweep:
        summary["seconds"] = round(time.time() - t_start, 1)
        return base_labels, summary

    fonts = []
    f = font0
    while f - FONT_STEP >= floor - 1e-9:
        f = round(f - FONT_STEP, 3)
        fonts.append(f)
    bands = [round(1.0 + BAND_STEP * k, 4) for k in range(1, int(round((band_cap - 1.0) / BAND_STEP)) + 1)]
    a_bands = [b for b in bands if b <= a_cap + 1e-9]

    spans = _feature_spans(render_today)
    span_of = lambda lab: spans.get(str(lab.get("stable_id")), float(abs(lab["end"] - lab["start"])))  # noqa: E731
    drop_order = sorted(ext0, key=lambda lab: (span_of(lab), str(lab.get("stable_id"))))
    summary["importance"] = "feature span (longest location segment, bp), shortest dropped first" if spans else "label extent (no feature span)"

    def drop_set(k):
        return frozenset(str(lab["stable_id"]) for lab in drop_order[:k])

    def finish(result, font, band, k, ok):
        labels, ev = result
        summary.update(font_size=font, band_scale=band, fits=ok, final=ev | {"bbox": ev["bbox"]})
        summary["dropped"] = [lab["label_text"] for lab in drop_order[:k]]
        summary["dropped_ids"] = [str(lab["stable_id"]) for lab in drop_order[:k]]
        summary["placed"] = ev["n"]
        summary["attempts"] = len(probe.log) - 1
        summary["seconds"] = round(time.time() - t_start, 1)
        summary["probes"] = probe.log
        return labels, summary

    def search_drop(font, band):
        """Fewest dropped labels (shortest features first) that fit at (font, band)."""
        path = [(font, band, frozenset())] + [(font, band, drop_set(k)) for k in range(1, len(drop_order) + 1)]
        idx, res = first_fit(path, probe)
        if idx is None:
            return len(drop_order), res, False
        return idx, res, True

    if design == "A1":
        path = [(font0, 1.0)] + [(f, 1.0) for f in fonts]
        idx, res = first_fit([(a, b, frozenset()) for a, b in path], probe)
        if idx is None:
            return finish(probe(floor, 1.0), floor, 1.0, 0, False)
        return finish(res, path[idx][0], 1.0, 0, True)
    if design == "A2":
        path = [(font0, 1.0)] + [(font0, b) for b in bands]
        idx, res = first_fit([(a, b, frozenset()) for a, b in path], probe)
        if idx is None:
            return finish(probe(font0, band_cap), font0, band_cap, 0, False)
        return finish(res, font0, path[idx][1], 0, True)
    if design == "A3":
        k, res, ok = search_drop(font0, 1.0)
        return finish(res, font0, 1.0, k, ok)
    # A: move along the frontier diagonal (about 0.1 band per 6 percent font), then font alone, then drop.
    steps = int(round((a_cap - 1.0) / BAND_STEP))
    diag = []
    for k in range(1, steps + 1):
        fk = max(floor, min(font0, math.floor(font0 * (1.0 - A_TRADE * k) / FONT_STEP + 1e-9) * FONT_STEP))
        diag.append((round(fk, 3), round(1.0 + BAND_STEP * k, 4)))
    tail = [(f, a_cap) for f in fonts if f < diag[-1][0] - 1e-9] if diag else []
    path = [(font0, 1.0)] + diag + tail
    seen, uniq = set(), []
    for item in path:  # drop duplicate points (font clipped at the floor)
        if item not in seen:
            seen.add(item)
            uniq.append(item)
    idx, res = first_fit([(a, b, frozenset()) for a, b in uniq], probe)
    if idx is not None:
        return finish(res, uniq[idx][0], uniq[idx][1], 0, True)
    ffloor, fb = uniq[-1]
    k, res, ok = search_drop(ffloor, fb)
    return finish(res, ffloor, fb, k, ok)


def _sweep(options, probe, font0, floor, band_cap, ext0, render_today):
    """Exploration for the report: evaluate every grid point along the font and band axes (no bisection)."""
    name = sys.argv[sys.argv.index("-o") + 1] if "-o" in sys.argv else "case"
    out = {"font0": font0, "axes": {}}
    axis = options.get("axis", "both")
    if axis in ("font", "both"):
        f = font0
        rows = []
        while f >= floor - 1e-9:
            rows.append(dict(zip(("font", "band", "drop"), (f, 1.0, 0)), **{k: v for k, v in probe(f, 1.0)[1].items() if k != "bbox"}))
            f = round(f - FONT_STEP, 3)
        out["axes"]["font"] = rows
    if axis in ("band", "both"):
        rows = []
        for k in range(0, int(round((band_cap - 1.0) / BAND_STEP)) + 1):
            b = round(1.0 + BAND_STEP * k, 4)
            rows.append(dict(zip(("font", "band", "drop"), (font0, b, 0)), **{k2: v for k2, v in probe(font0, b)[1].items() if k2 != "bbox"}))
        out["axes"]["band"] = rows
    if axis == "grid":
        rows = []
        f = font0
        fonts = []
        while f >= floor - 1e-9:
            fonts.append(f)
            f = round(f - 1.0, 3)
        for b in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5):
            for f in fonts:
                ev = probe(f, b)[1]
                rows.append({"font": f, "band": b, **{k: v for k, v in ev.items() if k != "bbox"}})
        out["axes"]["grid"] = rows
    out["probe_log"] = probe.log
    path = os.path.join(options.get("dumpdir", "."), f"{name}.sweep.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
