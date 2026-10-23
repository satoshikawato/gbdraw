"""Label fit warnings as the label tests read them."""

from __future__ import annotations

import warnings

from gbdraw.exceptions import LabelLeadersCrossWarning


def assert_no_label_warning(caught: list[warnings.WarningMessage]) -> None:
    """The labels fit; only a definition over a track may warn (D-12), which the label tests do not pin."""

    records = [
        record
        for warning in caught
        if isinstance(warning.message, LabelLeadersCrossWarning)
        for record in warning.message.report.records
    ]
    assert all(record.label_kind == "circular_definition" for record in records), records
