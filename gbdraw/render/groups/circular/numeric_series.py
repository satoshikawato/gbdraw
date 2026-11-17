"""The Circular group of a numeric series track (GC content, GC skew, depth): one group for every kind."""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from gbdraw.svg.elements import Group

from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import circular_display_length, unset_record_transform

from ....svg.circular_tracks import CircularTrackPlace
from ...drawers.circular.numeric_series import draw_circular_series
from ...numeric_series import NumericSeriesConfig, numeric_series_plot


class NumericSeriesGroup:
    """The group ``group_id`` of one record's numeric series track on a circular canvas: the series its
    configurator draws from ``frame`` (``numeric_series_plot``)."""

    def __init__(
        self,
        gb_record: SeqRecord,
        frame: DataFrame,
        radius: float,
        track_width: float,
        config: NumericSeriesConfig,
        *,
        norm_factor_override: float,
        group_id: str,
        record_transform: RecordDisplayTransform | None = None,
    ) -> None:
        # The transform the record is drawn by: whether its lines wrap round the origin or hold to 12 o'clock.
        transform = record_transform if record_transform is not None else unset_record_transform(gb_record)
        place = CircularTrackPlace(
            radius, circular_display_length(gb_record, transform), track_width, float(norm_factor_override), transform
        )
        self.group: Group = draw_circular_series(
            Group(id=group_id), numeric_series_plot(config, frame), place, group_id, gb_record.id
        )

    def get_group(self) -> Group:
        return self.group


__all__ = ["NumericSeriesGroup"]
