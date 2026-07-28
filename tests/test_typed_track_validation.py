from __future__ import annotations

import math

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.interface as interface
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularTrackOptions,
    DiagramOptions,
    TrackOptions,
)
from gbdraw.api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.tracks import CircularTrackSlot


def _record_input() -> RecordInput:
    return RecordInput(
        source=InMemoryRecordSource(
            SeqRecord(Seq("ATGC"), id="record"),
        )
    )


def test_linear_track_options_reject_circular_only_renderer() -> None:
    with pytest.raises(ValidationError, match="unknown linear track renderer"):
        interface.LinearTrackOptions(
            slots=("conservation:sequence_conservation",)
        )


@pytest.mark.parametrize(
    "factory",
    [
        lambda: interface.LinearTrackOptions(
            slots=("features:features",),
            axis_index=-1,
        ),
        lambda: interface.LinearTrackOptions(
            slots=("features:features",),
            axis_index=2,
        ),
        lambda: interface.CircularTrackOptions(
            slots=("features:features",),
            axis_index="0",  # type: ignore[arg-type]
        ),
        lambda: interface.CircularTrackOptions(
            slots=("features:features",),
            axis_index=True,  # type: ignore[arg-type]
        ),
    ],
)
def test_mode_specific_track_options_reject_invalid_axis_index(factory) -> None:
    with pytest.raises(ValidationError, match="axis_index"):
        factory()


@pytest.mark.parametrize("radius", [-1, math.nan, math.inf, "1", True])
def test_circular_track_options_reject_invalid_center_radius(radius: object) -> None:
    with pytest.raises(ValidationError, match="center_reserved_radius"):
        interface.CircularTrackOptions(
            center_reserved_radius=radius,  # type: ignore[arg-type]
        )


def test_shared_track_options_validate_slots_axes_and_radius() -> None:
    with pytest.raises(ValidationError, match="unknown linear track renderer"):
        TrackOptions(
            linear_track_slots=("conservation:sequence_conservation",)
        )
    with pytest.raises(ValidationError, match="linear_track_axis_index"):
        TrackOptions(
            linear_track_slots=("features:features",),
            linear_track_axis_index=-1,
        )
    with pytest.raises(ValidationError, match="circular_track_axis_index"):
        TrackOptions(
            circular_track_slots=("features:features",),
            circular_track_axis_index="0",  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="center_reserved_radius"):
        TrackOptions(center_reserved_radius=-1)


@pytest.mark.parametrize(
    "slot",
    [
        "gc:dinucleotide_content@__gbdraw_legacy_spacing=4px",
        "gc:dinucleotide_content@_auto_compress=false",
        CircularTrackSlot(
            id="gc",
            renderer="dinucleotide_content",
            params={"__gbdraw_legacy_spacing": "4px"},
        ),
        CircularTrackSlot(
            id="gc",
            renderer="dinucleotide_content",
            params={"_auto_compress": False},
        ),
    ],
)
def test_fresh_track_options_reject_private_circular_transport(
    slot: str | CircularTrackSlot,
) -> None:
    with pytest.raises(ValidationError, match="private"):
        TrackOptions(circular_track_slots=(slot,))


def test_linear_request_rejects_valid_circular_only_track_bundle() -> None:
    circular_tracks = CircularTrackOptions(
        circular_track_slots=("conservation:sequence_conservation",)
    )

    with pytest.raises(ValidationError, match="LinearDiagramOptions"):
        LinearDiagramRequest(
            records=(_record_input(),),
            options=CircularDiagramOptions(
                tracks=circular_tracks
            ),  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    "source",
    [
        1,
        (None, object()),
    ],
)
def test_root_depth_track_rejects_invalid_source_elements(source: object) -> None:
    with pytest.raises(ValidationError, match="source"):
        interface.DepthTrackOptions(source=source)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"depth_table": "depth.tsv"},
        {"depth_file": 1},
        {"depth_tables": ("depth.tsv",)},
        {"depth_files": (1,)},
        {"depth_track_tables": ((None, "depth.tsv"),)},
        {"depth_track_files": ((None, 1),)},
    ],
)
def test_diagram_options_reject_invalid_depth_source_types(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ValidationError, match="depth"):
        DiagramOptions(**kwargs)
