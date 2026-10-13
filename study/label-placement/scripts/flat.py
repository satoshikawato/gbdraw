#!/usr/bin/env python3
"""Flatten a transparent PNG onto white, optionally crop (x0,y0,x1,y1 as fractions) and scale. flat.py in out [scale] [x0 y0 x1 y1]"""
import sys
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
src, dst = sys.argv[1], sys.argv[2]
scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
im = Image.open(src).convert("RGBA")
bg = Image.new("RGBA", im.size, "white"); bg.alpha_composite(im); im = bg.convert("RGB")
if len(sys.argv) > 7:
    x0, y0, x1, y1 = (float(v) for v in sys.argv[4:8]); w, h = im.size
    im = im.crop((int(x0 * w), int(y0 * h), int(x1 * w), int(y1 * h)))
if scale != 1.0:
    im = im.resize((max(1, int(im.width * scale)), max(1, int(im.height * scale))), Image.LANCZOS)
im.save(dst)
