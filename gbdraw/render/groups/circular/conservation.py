"""Circular sequence-conservation ring group."""

from __future__ import annotations

from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]

from ....analysis.conservation import conservation_track_gradient_colors  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....render.drawers.circular.conservation import (  # type: ignore[reportMissingImports]
    ConservationDrawer,
    _safe_id_fragment,
)
from ....svg.circular_conservation import generate_full_annulus_path_desc
from ....svg.ids import stable_svg_id


class ConservationGroup:
    def __init__(
        self,
        *,
        hits: DataFrame,
        total_length: int,
        track_label: str,
        track_color: str | None,
        source_index: int,
        track_index: int,
        inner_radius_px: float,
        outer_radius_px: float,
        min_identity: float,
        cfg: GbdrawConfig,
        group_id: str | None = None,
        slot_id: str | None = None,
    ) -> None:
        safe_label = _safe_id_fragment(track_label)
        resolved_group_id = group_id or stable_svg_id(
            "conservation_ring",
            "circular-conservation-ring",
            int(source_index),
            int(track_index),
            str(track_label),
            namespace=safe_label or f"source_{int(source_index) + 1}",
        )
        self.group = Group(id=resolved_group_id, debug=False)
        if slot_id:
            self.group.attribs["data-gbdraw-slot-id"] = str(slot_id)
            self.group.attribs["data-gbdraw-slot-renderer"] = "sequence_conservation"
        self.group.attribs["data-source-index"] = str(source_index)
        self.group.attribs["data-track-index"] = str(track_index)
        self.group.attribs["data-track-label"] = str(track_label)
        if track_color:
            self.group.attribs["data-track-color"] = str(track_color)
        self.hits = hits
        self.total_length = int(total_length)
        self.inner_radius_px = float(inner_radius_px)
        self.outer_radius_px = float(outer_radius_px)
        self.min_identity = float(min_identity)
        self.cfg = cfg
        self.track_color = str(track_color) if track_color else None
        self._add_elements()

    def _add_background(self) -> None:
        conservation_cfg = self.cfg.objects.conservation
        if not bool(conservation_cfg.show_background):
            return
        path_desc = generate_full_annulus_path_desc(
            inner_radius_px=self.inner_radius_px,
            outer_radius_px=self.outer_radius_px,
            total_length=self.total_length,
        )
        self.group.add(
            Path(
                d=path_desc,
                fill=conservation_cfg.background_color,
                fill_opacity=conservation_cfg.background_opacity,
                stroke="none",
                stroke_width=0,
                debug=False,
            )
        )

    def _add_elements(self) -> None:
        self._add_background()
        conservation_cfg = self.cfg.objects.conservation
        min_color, max_color = conservation_track_gradient_colors(
            self.track_color,
            default_min_color=conservation_cfg.min_color,
            default_max_color=conservation_cfg.max_color,
        )
        drawer = ConservationDrawer(
            min_identity=self.min_identity,
            min_color=min_color,
            max_color=max_color,
            fill_opacity=conservation_cfg.fill_opacity,
            stroke_color=conservation_cfg.stroke_color,
            stroke_width=conservation_cfg.stroke_width,
        )
        drawer.draw_hits(
            self.group,
            self.hits,
            total_length=self.total_length,
            inner_radius_px=self.inner_radius_px,
            outer_radius_px=self.outer_radius_px,
        )

    def get_group(self) -> Group:
        return self.group


__all__ = ["ConservationGroup"]
