"""Layout helpers (internal)."""


from .linear import (
    AxisGapResolution,
    CollisionBand,
    required_axis_gap,
    resolve_axis_gap,
)
from .linear_multi_record import (
    LinearLayoutPlan,
    LinearRecordMeasurement,
    LinearRecordPlacement,
    RecordKey,
    parse_record_row_position,
    record_pairs_between_adjacent_rows,
    resolve_record_row_positions,
    solve_linear_layout,
    stable_record_keys,
)

__all__ = [
    "AxisGapResolution",
    "CollisionBand",
    "LinearLayoutPlan",
    "LinearRecordMeasurement",
    "LinearRecordPlacement",
    "RecordKey",
    "parse_record_row_position",
    "record_pairs_between_adjacent_rows",
    "required_axis_gap",
    "resolve_record_row_positions",
    "resolve_axis_gap",
    "solve_linear_layout",
    "stable_record_keys",
]
