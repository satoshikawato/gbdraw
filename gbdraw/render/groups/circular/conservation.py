"""Circular sequence-conservation ring group."""

from __future__ import annotations

from gbdraw.svg.elements import Group, Path

from ....analysis.conservation import ConservationTrack, conservation_track_gradient_colors
from ....comparisons.identity_scale import identity_color_domain
from ....config.models import GbdrawConfig
from ....layout.circular import RadialBand
from ....layout.circular_ink import CircularInk, CircularInkRing, paint_shows
from ....render.drawers.circular.conservation import (
    ConservationDrawer,
    _safe_id_fragment,
)
from ....svg.circular_conservation import generate_full_annulus_path_desc
from ....svg.ids import stable_svg_id
from ....layout.record_coordinates import RecordDisplayTransform


def _conservation_drawer(cfg: GbdrawConfig, track_color: str | None, min_identity: float) -> ConservationDrawer:
    conservation_cfg = cfg.objects.conservation
    min_color, max_color = conservation_track_gradient_colors(
        track_color,
        default_min_color=conservation_cfg.min_color,
        default_max_color=conservation_cfg.max_color,
    )
    return ConservationDrawer(
        identity_domain=identity_color_domain(min_identity, conservation_cfg.identity_scale),
        min_color=min_color,
        max_color=max_color,
        fill_opacity=conservation_cfg.fill_opacity,
        stroke_color=conservation_cfg.stroke_color,
        stroke_width=conservation_cfg.stroke_width,
    )


def conservation_ring_ink(
    track: ConservationTrack,
    *,
    total_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    min_identity: float,
    cfg: GbdrawConfig,
    record_transform: RecordDisplayTransform | None = None,
) -> tuple[CircularInk, ...]:
    """What ``ConservationGroup`` paints with the same arguments: the background ring when it shows (on, with a
    color and an opacity), else the hits."""

    conservation_cfg = cfg.objects.conservation
    if bool(conservation_cfg.show_background) and paint_shows(
        conservation_cfg.background_color, conservation_cfg.background_opacity
    ):
        return (CircularInkRing(RadialBand(float(inner_radius_px), float(outer_radius_px))),)
    drawer = _conservation_drawer(cfg, str(track.track_color) if track.track_color else None, float(min_identity))
    return drawer.ink(
        track.hits,
        total_length=int(total_length),
        inner_radius_px=inner_radius_px,
        outer_radius_px=outer_radius_px,
        record_transform=record_transform,
    )


class ConservationGroup:
    def __init__(
        self,
        *,
        track: ConservationTrack,
        total_length: int,
        inner_radius_px: float,
        outer_radius_px: float,
        min_identity: float,
        cfg: GbdrawConfig,
        group_id: str | None = None,
        slot_id: str | None = None,
        record_transform: RecordDisplayTransform | None = None,
        record_index: int = 0,
    ) -> None:
        track_label = track.track_label
        track_color = track.track_color
        source_index = track.source_index
        track_index = track.track_index
        safe_label = _safe_id_fragment(track_label)
        resolved_group_id = group_id or stable_svg_id(
            "conservation_ring",
            "circular-conservation-ring",
            int(source_index),
            int(track_index),
            str(track_label),
            namespace=safe_label or f"source_{int(source_index) + 1}",
        )
        self.group = Group(id=resolved_group_id)
        if slot_id:
            self.group.attribs["data-gbdraw-slot-id"] = str(slot_id)
            self.group.attribs["data-gbdraw-slot-renderer"] = "sequence_conservation"
        self.group.attribs["data-source-index"] = str(source_index)
        self.group.attribs["data-track-index"] = str(track_index)
        self.group.attribs["data-track-label"] = str(track_label)
        if track_color:
            self.group.attribs["data-track-color"] = str(track_color)
        self.hits = track.hits
        self.record_transform = record_transform
        self.record_index = record_index
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
            )
        )

    def _add_elements(self) -> None:
        self._add_background()
        drawer = _conservation_drawer(self.cfg, self.track_color, self.min_identity)
        drawer.draw_hits(
            self.group,
            self.hits,
            total_length=self.total_length,
            inner_radius_px=self.inner_radius_px,
            outer_radius_px=self.outer_radius_px,
            record_transform=self.record_transform,
            record_index=self.record_index,
        )

    def get_group(self) -> Group:
        return self.group


__all__ = ["ConservationGroup", "conservation_ring_ink"]
