"""The Linear group of a numeric series track (GC content, GC skew, depth): one group for every kind."""

from __future__ import annotations

from collections.abc import Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from gbdraw.layout.record_coordinates import RecordDisplayTransform

from ....svg.linear_tracks import LinearTrackFrame
from ...drawers.linear.numeric_series import draw_linear_series
from ...drawers.linear.wrap_rows import draw_track_rows
from ...numeric_series import NumericSeriesConfig, numeric_series_plot


class NumericSeriesGroup:
    """The group ``group_id`` of one record's numeric series track on a linear canvas, ``track_height`` high
    over the record's drawn axis width ``sequence_width``: once, or once per row of a wrapped record. The series
    is built once for every row."""

    def __init__(
        self,
        gb_record: SeqRecord,
        track_height: float,
        config: NumericSeriesConfig,
        *,
        sequence_width: float,
        frame: DataFrame,
        group_id: str,
        record_transform: RecordDisplayTransform | None = None,
        row_shifts: Sequence[float] = (),
    ) -> None:
        plot = numeric_series_plot(config, frame)
        self.group: Group = draw_track_rows(
            Group(id=group_id),
            record_transform,
            len(gb_record.seq),
            float(sequence_width),
            row_shifts,
            lambda group, length, width, row: draw_linear_series(
                group, plot, LinearTrackFrame(0, 0, length, width, track_height), group_id,
                record_transform, row,
            ),
        )

    def get_group(self) -> Group:
        return self.group


__all__ = ["NumericSeriesGroup"]
