"""Drawers for circular sequence-conservation rings."""

from __future__ import annotations

import re
from typing import Any

from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]

from ....core.color import interpolate_color
from ....svg.circular_conservation import generate_annular_hsp_path_desc


_SAFE_ID_PATTERN = re.compile(r"[^A-Za-z0-9_.-]+")


def _safe_id_fragment(value: object) -> str:
    text = str(value or "").strip()
    text = _SAFE_ID_PATTERN.sub("_", text).strip("_")
    return text


def _row_float(row: object, name: str, default: float = 0.0) -> float:
    try:
        return float(getattr(row, name))
    except (TypeError, ValueError, AttributeError):
        return float(default)


def _row_text(row: object, name: str) -> str:
    try:
        value = getattr(row, name)
    except AttributeError:
        return ""
    return str(value)


class ConservationDrawer:
    """Draw raw HSPs into one annular conservation ring."""

    def __init__(
        self,
        *,
        min_identity: float,
        min_color: str,
        max_color: str,
        fill_opacity: float,
        stroke_color: str,
        stroke_width: float,
    ) -> None:
        self.min_identity = float(min_identity)
        self.min_color = str(min_color)
        self.max_color = str(max_color)
        self.fill_opacity = float(fill_opacity)
        self.stroke_color = str(stroke_color)
        self.stroke_width = float(stroke_width)

    def _identity_factor(self, identity: float) -> float:
        denominator = max(1e-9, 100.0 - float(self.min_identity))
        factor = (float(identity) - float(self.min_identity)) / denominator
        return max(0.0, min(1.0, factor))

    def _fill_color(self, identity: float) -> str:
        return interpolate_color(
            self.min_color,
            self.max_color,
            self._identity_factor(identity),
        )

    def draw_hits(
        self,
        group: Group,
        hits: DataFrame,
        *,
        total_length: int,
        inner_radius_px: float,
        outer_radius_px: float,
    ) -> Group:
        if hits.empty or total_length <= 0:
            return group

        work_df = hits.copy()
        work_df["_draw_span"] = work_df["draw_end"].astype(float) - work_df["draw_start"].astype(float)
        work_df = work_df.sort_values(
            by=["identity", "_draw_span", "draw_start", "draw_end"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
        for row in work_df.itertuples(index=False):
            identity = _row_float(row, "identity")
            path_desc = generate_annular_hsp_path_desc(
                draw_start=_row_float(row, "draw_start"),
                draw_end=_row_float(row, "draw_end"),
                total_length=int(total_length),
                inner_radius_px=float(inner_radius_px),
                outer_radius_px=float(outer_radius_px),
                full_reference=bool(getattr(row, "full_reference", False)),
            )
            path = Path(
                d=path_desc,
                fill=self._fill_color(identity),
                fill_opacity=self.fill_opacity,
                stroke=self.stroke_color,
                stroke_width=self.stroke_width,
                debug=False,
            )
            metadata: dict[str, Any] = {
                "data-source-index": getattr(row, "source_index", ""),
                "data-track-label": getattr(row, "track_label", ""),
                "data-track-color": getattr(row, "track_color", ""),
                "data-identity": identity,
                "data-query": _row_text(row, "query"),
                "data-subject": _row_text(row, "subject"),
                "data-evalue": getattr(row, "evalue", ""),
                "data-bitscore": getattr(row, "bitscore", ""),
                "data-orientation": _row_text(row, "orientation"),
                "data-reference-record-id": _row_text(row, "reference_record_id"),
            }
            for attribute, value in metadata.items():
                path.attribs[attribute] = str(value)
            group.add(path)
        return group


__all__ = ["ConservationDrawer", "_safe_id_fragment"]
