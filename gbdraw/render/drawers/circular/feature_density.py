"""The heatmap band of a feature density track on a Circular canvas: one annular sector per run of equal windows."""

from __future__ import annotations

from collections import defaultdict

from gbdraw.analysis.feature_density import FeatureDensityLayer, feature_density_display_cells, heatmap_level_colors
from gbdraw.layout.circular_ink import CircularInkRing
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.circular import RadialBand
from gbdraw.svg.circular_conservation import generate_annular_hsp_path_desc
from gbdraw.svg.elements import Group, Path


class FeatureDensityHeatmapGroup:
    """The band between ``inner_radius_px`` and ``outer_radius_px``: one path per color step, each holding the
    sectors of that step's runs (adjacent windows of one step are one sector)."""

    def __init__(
        self,
        layer: FeatureDensityLayer,
        record_length: int,
        inner_radius_px: float,
        outer_radius_px: float,
        color: str,
        *,
        group_id: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        self.group = Group(id=group_id)
        colors = heatmap_level_colors(color, layer.levels)
        paths: dict[int, list[str]] = defaultdict(list)
        for cell in feature_density_display_cells(layer, record_transform):
            paths[cell.level].append(
                generate_annular_hsp_path_desc(
                    draw_start=cell.start,
                    draw_end=cell.end,
                    total_length=record_length,
                    inner_radius_px=inner_radius_px,
                    outer_radius_px=outer_radius_px,
                )
            )
        for level in sorted(paths):
            self.group.add(Path(d=" ".join(paths[level]), fill=colors[level - 1], stroke="none"))

    def get_group(self) -> Group:
        return self.group


def feature_density_heatmap_ink(
    layer: FeatureDensityLayer,
    record_length: int,
    band: RadialBand,
    record_transform: RecordDisplayTransform | None,
) -> tuple[CircularInkRing, ...]:
    """What :class:`FeatureDensityHeatmapGroup` paints: the band over each drawn run."""

    return tuple(
        CircularInkRing(band, cell.start / record_length, cell.end / record_length)
        for cell in feature_density_display_cells(layer, record_transform)
    )


__all__ = ["FeatureDensityHeatmapGroup", "feature_density_heatmap_ink"]
