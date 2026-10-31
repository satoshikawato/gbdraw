"""Drawers for circular sequence-conservation rings."""

from __future__ import annotations

import re
from typing import Any

from pandas import DataFrame
from gbdraw.svg.elements import Group, Path

from ....comparisons.identity_scale import IdentityColorDomain
from ....core.color import interpolate_color
from ....layout.circular import RadialBand
from ....layout.circular_ink import CircularInk, CircularInkRing, stroke_pad_px
from ....svg.circular_conservation import generate_annular_hsp_path_desc
from ....layout.record_coordinates import RecordDisplayTransform, SourceInterval


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


def _projected(record_transform: RecordDisplayTransform | None) -> bool:
    return record_transform is not None and record_transform.start_coordinate is not None


def _hit_parts(
    row: object, total_length: int, record_transform: RecordDisplayTransform | None
) -> tuple[tuple[float, float, bool], ...]:
    """The display spans one hit is drawn as, each with whether it is drawn as the whole ring."""

    spans: tuple[tuple[float, float], ...] = ((_row_float(row, "draw_start"), _row_float(row, "draw_end")),)
    if record_transform is not None and _projected(record_transform):
        spans = tuple((part.display_start, part.display_end) for part in
                      record_transform.project_interval(SourceInterval(int(spans[0][0]), int(spans[0][1]))))
        return tuple((start, end, end - start == total_length) for start, end in spans)
    return tuple((start, end, bool(getattr(row, "full_reference", False))) for start, end in spans)


class ConservationDrawer:
    """Draw raw HSPs into one annular conservation ring."""

    def __init__(
        self,
        *,
        identity_domain: IdentityColorDomain,
        min_color: str,
        max_color: str,
        fill_opacity: float,
        stroke_color: str,
        stroke_width: float,
    ) -> None:
        self.identity_domain = identity_domain
        self.min_color = str(min_color)
        self.max_color = str(max_color)
        self.fill_opacity = float(fill_opacity)
        self.stroke_color = str(stroke_color)
        self.stroke_width = float(stroke_width)

    def _identity_factor(self, identity: float) -> float:
        return self.identity_domain.factor(identity)

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
        record_transform: RecordDisplayTransform | None = None,
        record_index: int = 0,
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
            source_index = int(_row_float(row, "source_index"))
            source_hit_index = int(_row_float(row, "source_hit_index"))
            match_id = f"homology_ring{source_index + 1}_hit{source_hit_index + 1}"
            if record_index:
                match_id += f"_record_{record_index + 1}"
            metadata: dict[str, Any] = {
                "data-gbdraw-match-id": match_id,
                "data-match-kind": "homology",
                "data-source-index": getattr(row, "source_index", ""),
                "data-track-index": getattr(row, "track_index", ""),
                "data-track-label": getattr(row, "track_label", ""),
                "data-track-color": getattr(row, "track_color", ""),
                "data-reference-side": _row_text(row, "reference_side"),
                "data-identity": identity,
                "data-query": _row_text(row, "query"),
                "data-subject": _row_text(row, "subject"),
                "data-query-record-id": _row_text(row, "query"),
                "data-subject-record-id": _row_text(row, "subject"),
                "data-qstart": getattr(row, "qstart", ""),
                "data-qend": getattr(row, "qend", ""),
                "data-sstart": getattr(row, "sstart", ""),
                "data-send": getattr(row, "send", ""),
                "data-alignment-length": getattr(row, "alignment_length", ""),
                "data-evalue": getattr(row, "evalue", ""),
                "data-bitscore": getattr(row, "bitscore", ""),
                "data-mismatches": getattr(row, "mismatches", ""),
                "data-gap-opens": getattr(row, "gap_opens", ""),
                "data-orientation": _row_text(row, "orientation"),
                "data-reference-record-id": _row_text(row, "reference_record_id"),
            }
            projected = _projected(record_transform)
            if projected:
                metadata[f"data-{_row_text(row, 'reference_side')}-record-index"] = record_index
            for index, (start, end, full_reference) in enumerate(_hit_parts(row, total_length, record_transform)):
                path = Path(
                    d=generate_annular_hsp_path_desc(
                        draw_start=start, draw_end=end, total_length=int(total_length),
                        inner_radius_px=float(inner_radius_px), outer_radius_px=float(outer_radius_px),
                        full_reference=full_reference,
                    ),
                    fill=self._fill_color(identity), fill_opacity=self.fill_opacity,
                    stroke=self.stroke_color, stroke_width=self.stroke_width,
                )
                for attribute, value in metadata.items():
                    path.attribs[attribute] = str(value)
                if projected:
                    path.attribs["id"] = f"{match_id}__fragment{index}"
                    path.attribs["data-gbdraw-match-fragment"] = str(index)
                group.add(path)
        return group

    def ink(
        self,
        hits: DataFrame,
        *,
        total_length: int,
        inner_radius_px: float,
        outer_radius_px: float,
        record_transform: RecordDisplayTransform | None = None,
    ) -> tuple[CircularInk, ...]:
        """What ``draw_hits`` paints with the same arguments: each hit's sector, with half its stroke."""

        if hits.empty or total_length <= 0:
            return ()
        pad = stroke_pad_px(self.stroke_color, self.stroke_width)
        band = RadialBand(float(inner_radius_px) - pad, float(outer_radius_px) + pad)
        inks: list[CircularInk] = []
        for row in hits.itertuples(index=False):
            for start, end, full_reference in _hit_parts(row, total_length, record_transform):
                if full_reference:
                    return (CircularInkRing(band),)
                low, high = (max(0.0, min(float(total_length), value)) for value in (start, end))
                inks.append(CircularInkRing(band, low / total_length, high / total_length))
        return tuple(inks)


__all__ = ["ConservationDrawer", "_safe_id_fragment"]
