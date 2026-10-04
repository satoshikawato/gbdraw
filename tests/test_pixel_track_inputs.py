from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularTrackOptions,
    LinearDiagramOptions,
    LinearTrackOptions,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import read_circular_track_table
from gbdraw.session_request_codec import (
    CanonicalRequestDecodingError,
    decode_canonical_request,
    encode_canonical_request,
)
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    ScalarSpec,
    normalize_circular_track_slots,
    normalize_linear_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
    parse_linear_track_slot,
)
from gbdraw.tracks.scalars import parse_optional_pixel

CORPUS = json.loads(
    (Path(__file__).parent / "fixtures/pixel-track-inputs.json").read_text()
)


@pytest.mark.parametrize("case", CORPUS)
@pytest.mark.parametrize(
    "field,allow_zero",
    [
        ("height", False),
        ("spacing", True),
        ("inner_gap_px", True),
        ("outer_gap_px", True),
    ],
)
def test_pixel_grammar_and_text_slot_parity(case, field, allow_zero):
    value = float(case["numberToken"]) if "numberToken" in case else case["value"]
    expected = case["gap" if allow_zero else "height"]

    def parse():
        return parse_optional_pixel(value, field_name=field, allow_zero=allow_zero)

    if expected == "invalid":
        with pytest.raises(ValueError, match=rf"{field} must be .*finite.*px optional"):
            parse()
    else:
        assert parse() == expected
    if (
        isinstance(value, str)
        or value is None
        or (isinstance(value, (float, int)) and not isinstance(value, bool))
    ):
        adapter = (
            parse_linear_track_slot
            if field in {"height", "spacing"}
            else parse_circular_track_slot
        )
        slot_text = "gc:dinucleotide_content" + (
            f"@{field}={value}" if value is not None else ""
        )
        if expected == "invalid":
            with pytest.raises(ValueError, match=rf"{field} must be"):
                adapter(slot_text)
        else:
            slot = adapter(slot_text)
            normalized = (
                normalize_linear_track_slots
                if field in {"height", "spacing"}
                else normalize_circular_track_slots
            )([slot])[0]
            parsed = getattr(normalized, field)
            assert parsed == (
                ScalarSpec(expected, "px")
                if expected is not None and field in {"height", "spacing"}
                else expected
            )


@pytest.mark.parametrize(
    "text", ["10", "10px", "10PX", "10 px", "1e1px", ".5px", "0px", "  "]
)
def test_circular_tsv_text_adapter(text, tmp_path):
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tinner_gap_px\touter_gap_px\n"
        + f"gc\tdinucleotide_content\t{text}\t{text}\n"
    )
    slots = normalize_circular_track_slots(
        parse_circular_track_slots(read_circular_track_table(str(table)).slot_specs)
    )
    expected = parse_optional_pixel(text, field_name="inner_gap_px", allow_zero=True)
    assert slots[0].inner_gap_px == expected
    assert slots[0].outer_gap_px == expected


@pytest.mark.parametrize(
    "text", ["px", "-1px", "0x10", "10%", "10em", "NaN", "Infinity", "1e309px"]
)
def test_circular_tsv_invalid_is_fatal(text, tmp_path):
    table = tmp_path / "tracks.tsv"
    table.write_text(
        "id\trenderer\tinner_gap_px\n" + f"gc\tdinucleotide_content\t{text}\n"
    )
    with pytest.raises(ValidationError, match="inner_gap_px must be nonnegative"):
        read_circular_track_table(str(table))


def _encoded(mode):
    record = SeqRecord(Seq("ACGT" * 300), id="pixel-test", description="pixel-test")
    source = (RecordInput(InMemoryRecordSource(record)),)
    if mode == "circular":
        request = CircularDiagramRequest(
            records=source,
            options=CircularDiagramOptions(
                tracks=CircularTrackOptions(
                    circular_track_slots=(
                        CircularTrackSlot(
                            id="gc",
                            renderer="dinucleotide_content",
                            inner_gap_px=10,
                            outer_gap_px=0,
                        ),
                    )
                )
            ),
        )
    else:
        request = LinearDiagramRequest(
            records=source,
            options=LinearDiagramOptions(
                tracks=LinearTrackOptions(
                    linear_track_slots=(
                        LinearTrackSlot(
                            id="gc",
                            renderer="dinucleotide_content",
                            height=ScalarSpec(10, "px"),
                            spacing=ScalarSpec(0, "px"),
                        ),
                    )
                )
            ),
        )
    return encode_canonical_request(request)


@pytest.mark.parametrize(
    "mode,field",
    [
        ("circular", "innerGapPx"),
        ("circular", "outerGapPx"),
        ("linear", "height"),
        ("linear", "spacing"),
    ],
)
@pytest.mark.parametrize(
    "invalid", ["10", "10px", True, [], {}, float("inf"), float("nan")]
)
def test_typed_slot_reader_does_not_accept_pixel_text(mode, field, invalid, tmp_path):
    encoded = _encoded(mode)
    paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content)
        paths[resource.resource_id] = target
    payload = copy.deepcopy(encoded.payload)
    slots = payload["diagramOptions"]["tracks"][f"{mode}TrackSlots"]
    slots[0][field] = invalid
    with pytest.raises(CanonicalRequestDecodingError, match="typed contract"):
        decode_canonical_request(
            payload, resource_paths=paths, output_directory=tmp_path
        )
    decoded = decode_canonical_request(
        encoded.payload, resource_paths=paths, output_directory=tmp_path
    )
    slot = getattr(decoded.options.tracks, f"{mode}_track_slots")[0]
    assert getattr(slot, "inner_gap_px" if mode == "circular" else "height") == (
        10 if mode == "circular" else ScalarSpec(10, "px")
    )


@pytest.mark.parametrize(
    "value", [float("inf"), float("nan"), -1, 0, True, None, "10", "10px"]
)
def test_linear_native_scalar_height_must_be_positive_finite(value):
    with pytest.raises(ValueError, match="height must be positive finite"):
        normalize_linear_track_slots(
            [
                LinearTrackSlot(
                    id="gc",
                    renderer="dinucleotide_content",
                    height=ScalarSpec(value, "px"),
                )
            ]
        )


def test_factor_percent_and_retired_keys_remain_separate():
    slot = parse_circular_track_slot(
        "gc:dinucleotide_content@r=80%,w=0.1,inner_gap_px=10PX"
    )
    assert slot.radius == ScalarSpec(0.8, "factor")
    assert slot.width == ScalarSpec(0.1, "factor")
    assert slot.inner_gap_px == 10
    for key in ["spacing", "strict", "compress", "reserve"]:
        with pytest.raises(ValueError, match="no longer supported"):
            parse_circular_track_slot(f"gc:dinucleotide_content@{key}=1")
