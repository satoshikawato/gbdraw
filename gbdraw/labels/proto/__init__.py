"""Study prototypes for circular label placement (branch study/label-placement; never merged).

Selected with the environment variable GBDRAW_LABEL_PROTO:
  A1 / A2 / A3  today's placement plus a never-overprint fallback (font / band / drop)  -> design_a.py
  B             staggered tiers                                                          -> design_b.py
  C             free placement                                                           -> design_c.py
  D             side columns                                                             -> design_d.py
  E             cluster callout lists                                                    -> design_e.py
Options for a design come from GBDRAW_LABEL_PROTO_OPTS ("key=value,key=value").
The last run's summary is kept in LAST_RESULT for the study scripts.
"""

from __future__ import annotations

import importlib
import logging
import os

from .common import Ctx, make_items

logger = logging.getLogger(__name__)
LAST_RESULT: dict = {}


def proto_options() -> dict[str, str]:
    raw = os.environ.get("GBDRAW_LABEL_PROTO_OPTS", "")
    options: dict[str, str] = {}
    for part in raw.split(","):
        if "=" in part:
            key, value = part.split("=", 1)
            options[key.strip()] = value.strip()
    return options


def run_design(design, outer_labels, *, total_length, radius, gap_px, font_size, arc_x_radius, arc_y_radius,
               radial_layout=None, feature_colors=None):
    items = make_items(outer_labels, total_length)
    ctx = Ctx.build(
        items,
        total_length=total_length,
        radius=radius,
        gap=gap_px,
        font_size=font_size,
        arc_x_radius=arc_x_radius,
        arc_y_radius=arc_y_radius,
        options=proto_options(),
    )
    ctx.radial_layout, ctx.feature_colors = radial_layout, feature_colors or {}
    module = importlib.import_module(f".design_{design[0].lower()}", __package__)
    placed = module.place(items, ctx) if items else []
    LAST_RESULT.clear()
    LAST_RESULT.update(
        design=design, requested=len(items), placed=len(placed), dropped=list(ctx.dropped),
        font_size=font_size, notes=dict(ctx.notes),
    )
    logger.warning(
        "label prototype %s: %d labels requested, %d placed, %d dropped%s",
        design, len(items), len(placed), len(ctx.dropped),
        f" ({ctx.notes})" if ctx.notes else "",
    )
    return placed


def run_design_a(design, render_today, *, total_length):
    from . import design_a

    labels, summary = design_a.run(design, render_today, total_length=total_length, options=proto_options())
    LAST_RESULT.clear()
    LAST_RESULT.update(design=design, **summary)
    logger.warning("label prototype %s: %s", design, summary)
    return labels
