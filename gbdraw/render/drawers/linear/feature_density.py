"""The heatmap band of a feature density track on a Linear canvas: one rectangle per run of equal windows."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence

from gbdraw.analysis.feature_density import FeatureDensityLayer, feature_density_display_cells, heatmap_level_colors
from gbdraw.layout.linear_coords import normalize_position_to_linear_track
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.scalar_axis import ScalarRow
from gbdraw.svg.elements import Group, Path
from gbdraw.svg.path_data import path_number

from .wrap_rows import draw_track_rows


class FeatureDensityHeatmapGroup:
    """The band of ``track_height`` from y 0: one path per color step, each holding the rectangles of that step's
    runs (adjacent windows of one step are one rectangle), drawn on each row of a wrapped record."""

    def __init__(
        self,
        layer: FeatureDensityLayer,
        record_length: int,
        track_height: float,
        color: str,
        *,
        sequence_width: float,
        group_id: str,
        record_transform: RecordDisplayTransform | None = None,
        row_shifts: Sequence[float] = (),
    ) -> None:
        colors = heatmap_level_colors(color, layer.levels)
        cells = feature_density_display_cells(layer, record_transform)
        height = path_number(float(track_height))

        def draw(group: Group, length: int, width: float, row: ScalarRow | None) -> Group:
            row_index = 0 if row is None else row.row
            paths: dict[int, list[str]] = defaultdict(list)
            for cell in cells:
                if cell.row != row_index:
                    continue
                left = normalize_position_to_linear_track(float(cell.start), float(length), width)
                right = normalize_position_to_linear_track(float(cell.end), float(length), width)
                paths[cell.level].append(
                    f"M {path_number(left)},0 H {path_number(right)} V {height} H {path_number(left)} Z"
                )
            for level in sorted(paths):
                group.add(Path(d=" ".join(paths[level]), fill=colors[level - 1], stroke="none"))
            return group

        self.group = draw_track_rows(
            Group(id=group_id), record_transform, record_length, float(sequence_width), row_shifts, draw
        )

    def get_group(self) -> Group:
        return self.group


__all__ = ["FeatureDensityHeatmapGroup"]
