"""Which params keys a non-annotation track slot takes (OV-477), and the Linear content legend label (OV-846).

A slot takes the keys its renderer reads (with the parsers' aliases) and, on a
Circular ``sequence_conservation`` slot, the five keys the Web and the ring
builder store there (Owner LP-D66). Any other key is an error from every entry:
CLI text, Python slot objects and Sessions.
"""

from __future__ import annotations

import copy
import gzip
import json
import re
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularTrackOptions,
    LinearDiagramOptions,
    LinearTrackOptions,
)
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main
from gbdraw.session_request_codec import CanonicalRequestDecodingError, decode_canonical_request, encode_canonical_request
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    normalize_circular_track_slots,
    normalize_linear_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
    parse_linear_track_slot,
    parse_linear_track_slots,
)
from tests.test_session_request_codec import _materialize_resources

REPO_ROOT = Path(__file__).resolve().parents[1]

# The read set of each renderer (aliases included), from the readers on dev
# (label-placement ov477/STEP0.md section 1), plus the carried ring keys.
_LEGEND_LABEL = {"legend_label": "Label", "label": "Alias"}
_NT = {"nt": "AT", "dinucleotide": "AT"}
_SKEW_COLORS = {"positive_color": "red", "negative_color": "blue", "high_color": "red", "low_color": "blue"}
ACCEPTED_PARAMS: dict[str, dict[str, dict[str, Any]]] = {
    "linear": {
        "features": {},
        "dinucleotide_content": {**_NT, **_LEGEND_LABEL},
        "dinucleotide_skew": {**_NT, **_LEGEND_LABEL, **_SKEW_COLORS},
        "depth": {"track_index": 0, **_LEGEND_LABEL},
        "spacer": {},
    },
    "circular": {
        "features": {"lane_direction": "inside", "lanes": "inside", "stack_preset": "tuckin", "preset": "tuckin"},
        "ticks": {"tick_label_layout": "label_out_tick_in", "preset": "tuckin", "track_preset": "tuckin"},
        "dinucleotide_content": {**_NT, **_LEGEND_LABEL},
        "dinucleotide_skew": {**_NT, **_LEGEND_LABEL, **_SKEW_COLORS},
        "depth": {"track_index": 0, **_LEGEND_LABEL},
        # track_index and source_index are read; the other five are carried
        # as written by the Web (0.12.0 to now) and api/diagram.py (LP-D66).
        "sequence_conservation": {
            "track_index": "1",
            "source_index": "0",
            "managed": "circular_conservation",
            "series_key": "ring-a",
            "fileName": "a.tsv",
            "label": "Ring A",
            "color": "#e15759",
        },
        "spacer": {},
    },
}

_CASES = [
    (mode, renderer, key)
    for mode, renderers in ACCEPTED_PARAMS.items()
    for renderer in renderers
    for key in ("to_slot", "anchor_slot", "foo")
]
_IDS = [f"{mode}-{renderer}-{key}" for mode, renderer, key in _CASES]


def _expected_error(mode: str, renderer: str, key: str) -> str:
    if key in {"to_slot", "anchor_slot"}:
        return re.escape(f"{mode} track slot 'probe' ({renderer}) cannot use {key}: only annotation slots take")
    return re.escape(f"{mode} track slot 'probe' ({renderer}) has unknown parameter(s): {key}")


def _slot(mode: str, renderer: str, params: dict[str, Any]) -> LinearTrackSlot | CircularTrackSlot:
    cls = LinearTrackSlot if mode == "linear" else CircularTrackSlot
    return cls(id="probe", renderer=renderer, params=params)


def _features_slot(mode: str) -> LinearTrackSlot | CircularTrackSlot:
    cls = LinearTrackSlot if mode == "linear" else CircularTrackSlot
    return cls(id="features", renderer="features")


@pytest.mark.parametrize(("mode", "renderer", "key"), _CASES, ids=_IDS)
def test_a_slot_spec_text_rejects_a_key_its_renderer_does_not_read(mode: str, renderer: str, key: str) -> None:
    parse = parse_linear_track_slot if mode == "linear" else parse_circular_track_slot
    with pytest.raises(ValueError, match=_expected_error(mode, renderer, key)):
        parse(f"probe:{renderer}@{key}=features")


@pytest.mark.parametrize(("mode", "renderer", "key"), _CASES, ids=_IDS)
def test_a_slot_object_rejects_a_key_its_renderer_does_not_read(mode: str, renderer: str, key: str) -> None:
    slots = [_features_slot(mode), _slot(mode, renderer, {key: "features"})] if renderer != "features" else [
        _slot(mode, renderer, {key: "features"})
    ]
    parse = parse_linear_track_slots if mode == "linear" else parse_circular_track_slots
    normalize = normalize_linear_track_slots if mode == "linear" else normalize_circular_track_slots
    with pytest.raises(ValueError, match=_expected_error(mode, renderer, key)):
        parse(slots)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match=_expected_error(mode, renderer, key)):
        normalize(slots)  # type: ignore[arg-type]


def test_a_disabled_slot_rejects_an_unknown_key_too() -> None:
    with pytest.raises(ValueError, match=re.escape("has unknown parameter(s): bar, foo")):
        normalize_linear_track_slots(
            [LinearTrackSlot(id="probe", renderer="depth", enabled=False, params={"foo": 1, "bar": 2})]
        )


@pytest.mark.parametrize(("mode", "renderer"), [(m, r) for m, rs in ACCEPTED_PARAMS.items() for r in rs])
def test_a_slot_takes_every_key_its_renderer_reads(mode: str, renderer: str) -> None:
    params = ACCEPTED_PARAMS[mode][renderer]
    slots = [_slot(mode, renderer, dict(params))]
    if renderer != "features":
        slots.insert(0, _features_slot(mode))
    normalize = normalize_linear_track_slots if mode == "linear" else normalize_circular_track_slots
    normalized = normalize(slots)  # type: ignore[arg-type]
    assert normalized[-1].renderer == renderer


def _slot_payload(mode: str, slot_id: str, renderer: str, params: dict[str, Any]) -> dict[str, Any]:
    if mode == "circular":
        geometry: dict[str, Any] = {"radius": None, "width": None, "innerGapPx": None, "outerGapPx": None}
    else:
        geometry = {"height": None, "spacing": None}
    return {
        "kind": f"{mode}TrackSlot",
        "id": slot_id,
        "renderer": renderer,
        "enabled": True,
        "side": None,
        "z": 0,
        "params": params,
        **geometry,
    }


def _encoded_request(mode: str, renderer: str):
    record = SeqRecord(Seq("ATGC" * 50), id="record", annotations={"molecule_type": "DNA"})
    source = (RecordInput(source=InMemoryRecordSource(record)),)
    slots = [_slot(mode, renderer, {})]
    if renderer != "features":
        slots.insert(0, _features_slot(mode))
    if mode == "circular":
        return encode_canonical_request(CircularDiagramRequest(
            records=source,
            options=CircularDiagramOptions(tracks=CircularTrackOptions(circular_track_slots=tuple(slots))),
        ))
    return encode_canonical_request(LinearDiagramRequest(
        records=source,
        options=LinearDiagramOptions(tracks=LinearTrackOptions(linear_track_slots=tuple(slots))),
    ))


@pytest.mark.parametrize(("mode", "renderer", "key"), _CASES, ids=_IDS)
def test_a_session_request_rejects_a_key_its_renderer_does_not_read(
    mode: str, renderer: str, key: str, tmp_path: Path
) -> None:
    encoded = _encoded_request(mode, renderer)
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["tracks"][f"{mode}TrackSlots"][-1]["params"] = {key: "features"}
    with pytest.raises((CanonicalRequestDecodingError, ValueError), match=_expected_error(mode, renderer, key)):
        decode_canonical_request(
            payload,
            resource_paths=_materialize_resources(encoded, tmp_path / "resources"),
            output_directory=tmp_path / "output",
        )


# A Circular Session with a comparison ring whose slot carries the keys the Web
# wrote from 0.12.0 on (`refreshManagedConservationSlot`): it loads and draws
# exactly as the same Session without them.
RING_SESSION = REPO_ROOT / "tests/fixtures/sessions/precomputed-ring-unlabelled-cli.v44.gbdraw-session.json.gz"
RING_SLOT_KEYS_0130 = {
    "managed": "circular_conservation",
    "series_key": "caenorhabditis-human",
    "fileName": "caenorhabditis-human.tlosatx.tsv",
    "label": "caenorhabditis-human.tlosatx.tsv",
    "color": "",
}


@pytest.mark.circular
def test_a_session_ring_slot_with_the_web_bookkeeping_keys_draws_as_without_them(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    session = json.loads(gzip.decompress(RING_SESSION.read_bytes()))
    svgs = {}
    for variant, extra in (("read", {}), ("carried", RING_SLOT_KEYS_0130), ("unknown", {"foo": "bar"})):
        document = copy.deepcopy(session)
        document["renderRequest"]["diagramOptions"]["tracks"]["circularTrackSlots"] = [
            _slot_payload("circular", "features", "features", {}),
            _slot_payload("circular", "ticks", "ticks", {}),
            _slot_payload("circular", "ring", "sequence_conservation", {"source_index": "0", "track_index": "1", **extra}),
        ]
        path = tmp_path / f"{variant}.gbdraw-session.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        if variant == "unknown":
            with pytest.raises(Exception) as excinfo:
                circular_main(["--session", str(path), "-o", variant, "-f", "svg"])
            assert "track slot 'ring' (sequence_conservation) has unknown parameter(s): foo" in _error_text(
                excinfo.value
            )
            continue
        circular_main(["--session", str(path), "-o", variant, "-f", "svg"])
        svgs[variant] = (tmp_path / f"{variant}.svg").read_bytes()
    assert b'data-gbdraw-slot-renderer="sequence_conservation"' in svgs["read"]
    assert svgs["carried"] == svgs["read"]


def _error_text(error: BaseException) -> str:
    texts = []
    current: BaseException | None = error
    while current is not None:
        texts.append(str(current))
        current = current.__cause__ or current.__context__
    return " | ".join(texts)


# OV-846 (Owner LP-D67): a Linear content slot's legend_label (alias label)
# names its legend row, as a Linear skew slot's and a Circular slot's do.
def _legend_keys(svg: str) -> list[str]:
    return list(dict.fromkeys(re.findall(r'data-legend-key="([^"]*)"', svg)))


@pytest.mark.linear
def test_a_linear_content_slot_legend_label_names_its_legend_row(tmp_path: Path) -> None:
    linear_main([
        "--gbk", str(REPO_ROOT / "examples" / "LvMJNV.gb"), "-o", str(tmp_path / "lin"), "-f", "svg",
        "--linear_track_slot", "features:features",
        "--linear_track_slot", "gc_content:dinucleotide_content@legend_label=FOOCONTENT",
        "--linear_track_slot", "at_content:dinucleotide_content@nt=AT",
        "--linear_track_slot", "gc_skew:dinucleotide_skew@label=BARSKEW",
    ])
    keys = _legend_keys((tmp_path / "lin.svg").read_text(encoding="utf-8"))
    assert [key for key in keys if "content" in key.lower() or key.startswith(("FOO", "BAR"))] == [
        "FOOCONTENT",
        "AT content",
        "BARSKEW (+)",
        "BARSKEW (-)",
    ]


@pytest.mark.linear
def test_a_linear_legend_keeps_content_then_depth_then_skew_rows(tmp_path: Path) -> None:
    """The per-slot content rows keep dev's Linear order (review 1183 F3), whatever the slot order."""
    gbk = REPO_ROOT / "examples" / "LvMJNV.gb"
    record_id = next(line.split()[1] for line in gbk.read_text(encoding="utf-8").splitlines() if line.startswith("VERSION"))
    depth = tmp_path / "d.tsv"
    depth.write_text("".join(f"{record_id}\t{i}\t{i % 700 + 5}\n" for i in range(1, 280000, 500)), encoding="utf-8")
    linear_main([
        "--gbk", str(gbk), "--depth_track", str(depth), "-o", str(tmp_path / "lin"), "-f", "svg",
        "--linear_track_slot", "features:features",
        "--linear_track_slot", "gc_skew:dinucleotide_skew",
        "--linear_track_slot", "depth:depth",
        "--linear_track_slot", "gc_content:dinucleotide_content",
    ])
    keys = _legend_keys((tmp_path / "lin.svg").read_text(encoding="utf-8"))
    assert keys[keys.index("GC content"):] == ["GC content", "Depth", "GC skew (+)", "GC skew (-)"]
