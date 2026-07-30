from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram
from gbdraw.api.options import CircularDiagramOptions, CircularTrackOptions
from gbdraw.tracks import default_circular_track_slots


def _record(record_id: str = "record") -> SeqRecord:
    return SeqRecord(Seq("ATGC" * 25), id=record_id)


def _horizontal_inner_label_options(
    *,
    explicit_slots: bool = False,
) -> CircularDiagramOptions:
    return CircularDiagramOptions(
        config_overrides={
            "labels.circular.scope": "both",
            "labels.circular.placement": "horizontal",
            "canvas.show_gc": True,
            "canvas.show_skew": True,
        },
        tracks=(
            CircularTrackOptions(
                circular_track_slots=default_circular_track_slots(
                    show_gc=True,
                    show_skew=True,
                )
            )
            if explicit_slots
            else None
        ),
    )


@pytest.mark.parametrize("multi_record", [False, True])
def test_horizontal_inner_labels_suppress_implicit_gc_and_skew_tracks(
    multi_record: bool,
) -> None:
    options = _horizontal_inner_label_options()
    drawing = (
        build_circular_multi_diagram([_record("a"), _record("b")], options=options)
        if multi_record
        else build_circular_diagram(_record(), options=options)
    )

    svg = drawing.tostring()
    assert 'id="gc_content"' not in svg
    assert 'id="skew"' not in svg
    assert 'id="gc_skew"' not in svg


def test_explicit_circular_slots_remain_authoritative_with_horizontal_inner_labels() -> None:
    svg = build_circular_diagram(
        _record(),
        options=_horizontal_inner_label_options(explicit_slots=True),
    ).tostring()

    assert 'id="gc_content"' in svg
    assert 'id="gc_skew"' in svg
