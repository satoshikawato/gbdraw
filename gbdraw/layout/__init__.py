"""Layout helpers (internal)."""


from .linear import (
    AxisGapResolution,
    CollisionBand,
    LinearRecordRenderContext,
    LinearResolvedTrack,
    LinearTrackLayout,
    required_axis_gap,
    resolve_axis_gap,
)
from .linear_multi_record import (
    LinearLayoutPlan,
    LinearRecordMeasurement,
    LinearRecordPlacement,
    RecordKey,
    record_pairs_between_adjacent_rows,
    solve_linear_layout,
    stable_record_keys,
)
from .record_placement import parse_record_row_position, resolve_record_row_positions

__all__ = [
    "AxisGapResolution",
    "CollisionBand",
    "LinearLayoutPlan",
    "LinearRecordRenderContext",
    "LinearRecordMeasurement",
    "LinearRecordPlacement",
    "LinearResolvedTrack",
    "LinearTrackLayout",
    "RecordKey",
    "parse_record_row_position",
    "record_pairs_between_adjacent_rows",
    "required_axis_gap",
    "resolve_record_row_positions",
    "resolve_axis_gap",
    "solve_linear_layout",
    "stable_record_keys",
]
