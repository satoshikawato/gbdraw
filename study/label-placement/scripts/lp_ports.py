#!/usr/bin/env python3
"""LP-D23: the 12 points where a leader may meet its label, on the bounding box of the drawn letters.

"bbox" is read as the box Inkscape shows for a text: ink left to ink right, top of the tallest letter to the lowest
ink (the baseline, or the bottom of a descender). Measured from the font, for horizontal text.
"""
import math
from functools import lru_cache
from fontTools.ttLib import TTFont
FONT = '/home/kawato/gbdraw-work/gbdraw/data/LiberationSans-Regular.ttf'

@lru_cache(None)
def _font():
    f = TTFont(FONT); return f, f.getBestCmap(), f['glyf'], f['hmtx'], f['head'].unitsPerEm, f['OS/2'].sxHeight

def letter_box(text, x, y, size=14.0, anchor='start'):
    """(x0, top, x1, baseline, descender bottom) for text drawn at (x, y) with dominant-baseline="middle"."""
    f, cmap, glyf, hmtx, upm, xh = _font(); s = size / upm
    adv = 0; x0 = None; x1 = 0; top = 0; bot = 0
    for ch in text:
        g = cmap[ord(ch)]; a, lsb = hmtx[g]; gl = glyf[g]
        if gl.numberOfContours:
            if x0 is None:
                x0 = adv + gl.xMin
            x1 = adv + gl.xMax; top = max(top, gl.yMax); bot = min(bot, gl.yMin)
        adv += a
    if anchor == 'end':
        x -= adv * s
    elif anchor == 'middle':
        x -= adv * s / 2
    base = y + xh * s / 2                      # "middle" = half the x-height above the alphabetic baseline
    return (x + (x0 or 0) * s, base - top * s, x + x1 * s, base, base - bot * s)

def ink_box(text, x, y, size=14.0, anchor='start'):
    """(x0, top, x1, bottom) of the ink, descenders included."""
    b = letter_box(text, x, y, size, anchor)
    return (b[0], b[1], b[2], b[4])

def ports(box):
    """Name -> point. Horizontal text: short sides are left and right, long sides top and bottom."""
    x0, y0, x1, y1 = box[:4]; mx, my = (x0 + x1) / 2, (y0 + y1) / 2; q1, q3 = x0 + (x1 - x0) / 4, x0 + 3 * (x1 - x0) / 4
    return {'TL': (x0, y0), 'TR': (x1, y0), 'BL': (x0, y1), 'BR': (x1, y1), 'L': (x0, my), 'R': (x1, my),
            'T1': (q1, y0), 'T2': (mx, y0), 'T3': (q3, y0), 'B1': (q1, y1), 'B2': (mx, y1), 'B3': (q3, y1)}

def stop_short(a, port, gap):
    """End of a leader drawn from `a` toward `port`, `gap` px before it."""
    d = math.dist(a, port)
    return (port[0] - (port[0] - a[0]) * gap / d, port[1] - (port[1] - a[1]) * gap / d)
