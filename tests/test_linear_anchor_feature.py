"""Align Linear records at the start of one feature (``--anchor_feature``).

The selector (``key=value``, exact, the feature-selector form of the
annotation and placement tables) is resolved on every record as drawn. A
complete circular record that has the feature is rotated so that the feature
starts at its left end (PD-OI-113 scenario revision 3); any other record that
has it (linear, cropped, or drawn linear) moves sideways so that the feature's
start, as drawn, lies below the reference's: the first rotated record, drawn at
x 0, else the first record with the feature, which keeps its place. The figure
equals the one drawn with the same starts and moves given explicitly.

The two mitogenomes are circular, so the sideways cases draw them linear
(``RecordDisplayOptions(is_circular=False)``).
"""
from __future__ import annotations

import logging
import os
import re
import subprocess
import sys
from dataclasses import fields, replace
from pathlib import Path

import pytest

from gbdraw import LinearOptions, draw_linear, read_genbank
from gbdraw import linear as linear_cli
from gbdraw.api import (
    FeatureSelector,
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    plan_linear_request,
)
from gbdraw.api.options import LinearMultiRecordOptions, LinearRecordTranslation
from gbdraw.api.request_render import LinearRequestPlan, build_request_diagram
from gbdraw.api.requests import RecordDisplayOptions, RecordPresentation
from gbdraw.exceptions import ValidationError

_DATA = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data"
_HUMAN = _DATA / "human-mitochondrion" / "HmmtDNA.gbk"
_ZEBRAFISH = _DATA / "metazoan-mitochondria-four" / "NC_002333.2.gb"
_EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
_COX1 = FeatureSelector(value="COX1", key="gene")
# Human COX1 is 5904..7445 and zebrafish COX1 6425..7975 (gene and CDS alike).
_HUMAN_COX1_START = 5903.0
_ZEBRAFISH_COX1_START = 6424.0
_FOUR = (
    _HUMAN, _ZEBRAFISH,
    _DATA / "metazoan-mitochondria-four" / "NC_001328.1.gb",
    _DATA / "metazoan-mitochondria-four" / "NC_024511.2.gb",
)
# The COX1 gene starts of the four, all on the plus strand.
_FOUR_COX1_STARTS = (5904, 6425, 7845, 1474)
_LINEAR = RecordDisplayOptions(is_circular=False)


def _request(**options: object) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(_HUMAN)), RecordInput(GenBankInputSource(_ZEBRAFISH))),
        options=LinearDiagramOptions(**options),
    )


def _drawn_linear(request: LinearDiagramRequest, *indices: int) -> LinearDiagramRequest:
    """``request`` with the records at ``indices`` (all without any) drawn linear: the anchor moves those sideways."""

    chosen = set(indices or range(len(request.records)))
    return replace(request, records=tuple(
        replace(record, display=_LINEAR) if index in chosen else record
        for index, record in enumerate(request.records)
    ))


def _starts(plan: LinearRequestPlan) -> tuple[int | None, ...]:
    """Each planned record's display start: a rotated record's is where the feature starts."""

    return tuple(item.display.start_coordinate for item in plan.provenance)


def _svg(request: LinearDiagramRequest) -> str:
    return build_request_diagram(request).drawing.tostring()


def _record_shifts(svg: str) -> list[float]:
    """Each record group's own x translation, in record order."""

    return [float(match) for match in re.findall(r'data-record-translation-x="([-0-9.e]+)"', svg)]


@pytest.mark.linear
def test_anchor_feature_rotates_complete_circular_records_like_explicit_starts() -> None:
    request = LinearDiagramRequest(
        records=tuple(RecordInput(GenBankInputSource(path)) for path in _FOUR),
        options=LinearDiagramOptions(anchor_feature=_COX1),
    )
    plan = plan_linear_request(request)
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.statuses == ("rotated",) * 4
    assert plan.feature_anchor.reference_index == 0
    assert _starts(plan) == _FOUR_COX1_STARTS

    from gbdraw.api import resolve_request

    resolved = resolve_request(request)
    assert isinstance(resolved, LinearDiagramRequest) and resolved.layout is not None
    assert [record.display.start_coordinate for record in resolved.records] == list(_FOUR_COX1_STARTS)
    assert [item.x for item in resolved.layout.record_translations] == [0.0] * 4
    explicit = LinearDiagramRequest(records=tuple(
        RecordInput(GenBankInputSource(path), display=RecordDisplayOptions(start_coordinate=start))
        for path, start in zip(_FOUR, _FOUR_COX1_STARTS)
    ))
    anchored = _svg(request)
    assert _record_shifts(anchored) == [0.0] * 4
    assert anchored == _svg(explicit)


@pytest.mark.linear
def test_anchor_feature_moves_linear_records_like_explicit_translations() -> None:
    request = _drawn_linear(_request(anchor_feature=_COX1))
    plan = plan_linear_request(request)
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.starts == (_HUMAN_COX1_START, _ZEBRAFISH_COX1_START)
    assert plan.feature_anchor.statuses == ("reference", "aligned")
    assert _starts(plan) == (None, None)

    anchored = _svg(request)
    plain = _svg(_drawn_linear(_request()))
    assert anchored != plain
    assert _record_shifts(plain) == [0.0, 0.0]
    # The human record keeps its place; the zebrafish record moves left by the
    # 521 bp between the two COX1 starts, at the longer record's scale.
    first, shift = _record_shifts(anchored)
    assert first == 0.0
    assert shift == pytest.approx(-(_ZEBRAFISH_COX1_START - _HUMAN_COX1_START) * 2000 / 16596)

    keys = tuple(str(record.record_key) for record in plan.request.records)
    explicit = replace(
        plan.request,
        options=replace(plan.request.options, anchor_feature=None),
        layout=LinearMultiRecordOptions(record_translations=(
            LinearRecordTranslation(record_key=keys[0]),
            LinearRecordTranslation(record_key=keys[1], x=shift),
        )),
    )
    assert _svg(explicit) == anchored


@pytest.mark.linear
def test_linear_records_line_up_with_the_left_end_of_the_rotated_records() -> None:
    # Q6-1 A: the human record drawn linear comes first; the zebrafish record is
    # rotated and is the reference at x 0, and the human record moves left so
    # that its COX1 starts there, beyond the left of the others.
    plan = plan_linear_request(_drawn_linear(_request(anchor_feature=_COX1), 0))
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.statuses == ("aligned", "rotated")
    assert plan.feature_anchor.reference_index == 1
    assert _starts(plan) == (None, 6425)
    human, zebrafish = _record_shifts(_svg(plan.request))
    assert zebrafish == 0.0
    assert human == pytest.approx(-_HUMAN_COX1_START * 2000 / 16596)


@pytest.mark.linear
def test_a_cropped_or_linear_record_moves_sideways() -> None:
    from gbdraw.io.regions import parse_region_spec

    request = _request(anchor_feature=_COX1)
    cropped = replace(request, records=(
        replace(request.records[0], region=parse_region_spec("NC_012920.1:5000-9000")),
        request.records[1],
    ))
    plan = plan_linear_request(cropped)
    assert plan.feature_anchor is not None
    # Human COX1 (5904) is base index 904 of the region; the zebrafish record is rotated.
    assert plan.feature_anchor.statuses == ("aligned", "rotated")
    assert plan.feature_anchor.starts == (904.0, 0.0)
    assert _starts(plan) == (None, 6425)


@pytest.mark.linear
def test_a_rotated_record_replaces_its_start_and_an_unmatched_one_keeps_its_own() -> None:
    request = _request(anchor_feature=FeatureSelector(value="RNR1", key="gene"))
    started = replace(request, records=(
        replace(request.records[0], display=RecordDisplayOptions(start_coordinate=7000)),
        replace(request.records[1], display=RecordDisplayOptions(start_coordinate=100)),
    ))
    plan = plan_linear_request(started)
    # Only the human record names its 12S rRNA gene RNR1 (648..1601).
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.statuses == ("rotated", "unmatched")
    assert _starts(plan) == (648, 100)


@pytest.mark.linear
def test_anchor_feature_is_resolved_on_the_record_as_drawn() -> None:
    request = _request(anchor_feature=_COX1)
    reversed_second = replace(request, records=(
        request.records[0],
        replace(request.records[1], presentation=RecordPresentation(reverse_complement=True)),
    ))
    # Reversed, zebrafish COX1 (6425..7975 of 16596) starts as drawn at its 3' base 7975.
    assert _starts(plan_linear_request(reversed_second)) == (5904, 7975)
    plan = plan_linear_request(_drawn_linear(reversed_second))
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.starts == (_HUMAN_COX1_START, 16596.0 - 7975.0)


@pytest.mark.linear
def test_a_record_without_the_feature_keeps_its_place(caplog: pytest.LogCaptureFixture) -> None:
    lambda_genbank = _DATA / "lambda" / "NC_001416.gb"
    request = LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(_HUMAN)),
            RecordInput(GenBankInputSource(lambda_genbank)),
            RecordInput(GenBankInputSource(_ZEBRAFISH)),
        ),
        options=LinearDiagramOptions(anchor_feature=_COX1),
    )
    with caplog.at_level(logging.WARNING):
        plan = plan_linear_request(request)
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.statuses == ("rotated", "unmatched", "rotated")
    assert _starts(plan) == (5904, None, 6425)
    assert "anchor_feature gene=COX1 matched no feature in record #2" in caplog.text


@pytest.mark.linear
def test_anchor_feature_that_matches_no_record_is_refused() -> None:
    with pytest.raises(ValidationError, match=r"anchor_feature gene=nonexistent matched no feature in any record"):
        plan_linear_request(_request(anchor_feature=FeatureSelector(value="nonexistent", key="gene")))


@pytest.mark.linear
def test_anchor_feature_that_matches_features_at_different_starts_is_refused() -> None:
    request = LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(_EXAMPLES / "MjeNMV.gb")),
            RecordInput(GenBankInputSource(_EXAMPLES / "MelaMJNV.gb")),
        ),
        options=LinearDiagramOptions(anchor_feature=FeatureSelector(value="dUTPase", key="product")),
    )
    with pytest.raises(ValidationError) as caught:
        plan_linear_request(request)
    message = str(caught.value)
    assert "anchor_feature product=dUTPase matches 2 features at different starts in record #1" in message
    assert "CDS 249324..249569; CDS 250088..250528" in message


@pytest.mark.linear
def test_anchor_feature_wins_over_align_center() -> None:
    center = {"canvas.linear.align_center": True}
    centered = _svg(_drawn_linear(_request(config_overrides=center)))
    anchored = _svg(_drawn_linear(_request(anchor_feature=_COX1, config_overrides=center)))
    plain_anchored = _svg(_drawn_linear(_request(anchor_feature=_COX1)))
    # The first record keeps its centred place; the second record moves from
    # its own centred place (27 bp right of the first, half of the 27 bp
    # length difference) so that its COX1 start sits under the first one's.
    assert _record_shifts(centered) == [0.0, 0.0]
    first, shift = _record_shifts(anchored)
    assert first == 0.0
    assert shift == pytest.approx(_record_shifts(plain_anchored)[1] + 13.5 * 2000 / 16596)
    assert anchored != plain_anchored
    # Rotated records line up their left ends: the zebrafish record leaves its
    # centred place, 13.5 bp left of the human one's, for the human left end.
    rotated = _record_shifts(_svg(_request(anchor_feature=_COX1, config_overrides=center)))
    assert rotated[0] == 0.0
    assert rotated[1] == pytest.approx(13.5 * 2000 / 16596)


@pytest.mark.linear
def test_anchor_feature_cannot_combine_with_a_similarity_alignment() -> None:
    from gbdraw.layout.similarity_alignment import SimilarityAlignmentReference
    from gbdraw.api.options import LosatSearchOptions

    with pytest.raises(ValidationError, match="anchor_feature cannot be combined with a similarity alignment"):
        LinearDiagramRequest(
            records=(RecordInput(GenBankInputSource(_HUMAN)), RecordInput(GenBankInputSource(_ZEBRAFISH))),
            options=LinearDiagramOptions(
                anchor_feature=_COX1,
                losat_search=LosatSearchOptions(program="losatp", losatp_mode="similarity_groups"),
            ),
            similarity_alignment=SimilarityAlignmentReference(feature_id="x"),
        )


@pytest.mark.linear
def test_cli_anchor_feature_draws_the_typed_figure(tmp_path: Path) -> None:
    linear_cli.linear_main([
        "--gbk", str(_HUMAN), str(_ZEBRAFISH), "-f", "svg", "-o", str(tmp_path / "cli"),
        "--anchor_feature", "gene=COX1",
    ])
    linear_cli.linear_main([
        "--gbk", str(_HUMAN), str(_ZEBRAFISH), "-f", "svg", "-o", str(tmp_path / "plain"),
    ])
    cli = (tmp_path / "cli.svg").read_text(encoding="utf-8")
    plain = (tmp_path / "plain.svg").read_text(encoding="utf-8")
    assert _record_shifts(cli) == _record_shifts(plain) == [0.0, 0.0]
    assert cli.split("?>", 1)[1].strip() == _svg(_request(anchor_feature=_COX1))


@pytest.mark.linear
def test_python_anchor_feature_draws_the_typed_figure() -> None:
    records = read_genbank([str(_HUMAN), str(_ZEBRAFISH)])
    anchored = draw_linear(records, options=LinearOptions(anchor_feature="gene=COX1"))
    plain = draw_linear(records)
    assert _record_shifts(anchored.to_svg()) == _record_shifts(_svg(_request(anchor_feature=_COX1)))
    assert _record_shifts(plain.to_svg()) == [0.0, 0.0]
    assert anchored.to_svg() != plain.to_svg()


def _cli(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", *args],
        cwd=cwd, capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
    )


def _stored_request(session_path: Path) -> dict[str, object]:
    import json

    request = json.loads(session_path.read_text(encoding="utf-8"))["drawings"][0]["renderRequest"]
    assert "anchorFeature" not in request["diagramOptions"]
    return request


def _stored_layout(session_path: Path) -> dict[str, object]:
    return _stored_request(session_path)["layout"]


@pytest.mark.linear
def test_a_session_stores_the_starts_and_translations_the_anchor_resolves_to(tmp_path: Path) -> None:
    saved = _cli(
        "--gbk", str(_HUMAN), str(_ZEBRAFISH), "-f", "svg", "-o", "saved",
        "--anchor_feature", "gene=COX1", "--save_session", cwd=tmp_path,
    )
    assert saved.returncode == 0, saved.stderr
    stored = _stored_request(tmp_path / "saved.gbdraw-session.json")
    drawn = (tmp_path / "saved.svg").read_text(encoding="utf-8")
    assert [entry["x"] for entry in stored["layout"]["recordTranslations"]] == _record_shifts(drawn) == [0.0, 0.0]
    assert [record["display"]["startCoordinate"] for record in stored["records"]] == [5904, 6425]
    replayed = _cli("--session", "saved.gbdraw-session.json", "-o", "again", cwd=tmp_path)
    assert replayed.returncode == 0, replayed.stderr
    assert (tmp_path / "again.svg").read_bytes() == (tmp_path / "saved.svg").read_bytes()


@pytest.mark.linear
def test_a_session_of_moved_linear_records_stores_their_translations(tmp_path: Path) -> None:
    saved = _cli(
        "--gbk", str(_EXAMPLES / "MjeNMV.gb"), str(_EXAMPLES / "MelaMJNV.gb"),
        "--region", "LC738868.1:249000-249800", "--region", "LC738874.1:234800-235400",
        "-f", "svg", "-o", "saved", "--anchor_feature", "product=dUTPase", "--save_session", cwd=tmp_path,
    )
    assert saved.returncode == 0, saved.stderr
    layout = _stored_layout(tmp_path / "saved.gbdraw-session.json")
    drawn = (tmp_path / "saved.svg").read_text(encoding="utf-8")
    assert [entry["x"] for entry in layout["recordTranslations"]] == _record_shifts(drawn)
    assert _record_shifts(drawn)[1] > 0.0
    replayed = _cli("--session", "saved.gbdraw-session.json", "-o", "again", cwd=tmp_path)
    assert replayed.returncode == 0, replayed.stderr
    assert (tmp_path / "again.svg").read_bytes() == (tmp_path / "saved.svg").read_bytes()


@pytest.mark.linear
def test_the_python_session_writers_store_the_starts_and_translations_the_anchor_resolves_to(tmp_path: Path) -> None:
    from gbdraw.session import (
        build_session_document,
        load_session_document,
        materialize_session,
        save_session_document,
        session_to_request,
    )

    for index, request in enumerate((_request(anchor_feature=_COX1), _drawn_linear(_request(anchor_feature=_COX1), 0))):
        drawn = _svg(request)
        built = build_session_document(request).to_dict()
        path = tmp_path / f"saved-{index}.gbdraw-session.json"
        save_session_document(path, request)
        assert _stored_request(path)["layout"] == built["drawings"][0]["renderRequest"]["layout"]
        assert [entry["x"] for entry in _stored_layout(path)["recordTranslations"]] == _record_shifts(drawn)
        assert _stored_request(path)["records"][1]["display"]["startCoordinate"] == 6425
        with materialize_session(load_session_document(path), output_directory=tmp_path / f"out-{index}") as materialized:
            replayed = session_to_request(materialized)
            assert isinstance(replayed, LinearDiagramRequest) and replayed.options.anchor_feature is None
            assert _svg(replayed) == drawn


@pytest.mark.linear
@pytest.mark.parametrize("config_overrides", [None, {"canvas.linear.align_center": True}])
@pytest.mark.parametrize("linear", [False, True])
def test_the_resolved_request_draws_the_anchored_figure(
    config_overrides: dict[str, object] | None, linear: bool,
) -> None:
    from gbdraw.api import resolve_request
    from gbdraw.session_request_codec import encode_canonical_request

    request = _request(anchor_feature=_COX1, config_overrides=config_overrides)
    if linear:
        request = _drawn_linear(request)
    # The request itself is not encoded: the anchor is resolved when the figure is laid out.
    with pytest.raises(ValidationError, match="anchor_feature record placement"):
        encode_canonical_request(request)
    resolved = resolve_request(request)
    assert isinstance(resolved, LinearDiagramRequest)
    assert resolved.options.anchor_feature is None
    payload = encode_canonical_request(resolved).payload
    assert [entry["x"] for entry in payload["layout"]["recordTranslations"]] == _record_shifts(_svg(request))
    assert [record["display"]["startCoordinate"] for record in payload["records"]] == (
        [None, None] if linear else [5904, 6425]
    )
    assert _svg(resolved) == _svg(request)


@pytest.mark.linear
def test_the_anchor_replaces_the_stored_x_of_the_records_it_places() -> None:
    from gbdraw.api import resolve_request

    for linear in (False, True):
        request = _request(anchor_feature=_COX1)
        planned = plan_linear_request(_drawn_linear(request) if linear else request).request
        keys = tuple(str(record.record_key) for record in planned.records)
        held = LinearMultiRecordOptions(record_translations=(
            LinearRecordTranslation(record_key=keys[0], x=40.0, y=5.0),
            LinearRecordTranslation(record_key=keys[1], x=-300.0, y=7.0),
        ))
        request = replace(planned, layout=held)
        resolved = resolve_request(request)
        assert isinstance(resolved, LinearDiagramRequest) and resolved.layout is not None
        first, second = resolved.layout.record_translations
        if linear:
            shift = -(_ZEBRAFISH_COX1_START - _HUMAN_COX1_START) * 2000 / 16596
            # The first record keeps its own translation; the moved record's x follows it, its y stays.
            assert (first.x, first.y) == (40.0, 5.0)
            assert second.x == pytest.approx(40.0 + shift) and second.y == 7.0
        else:
            # A rotated reference is drawn at x 0, and the rotated record lines up with it; the y stay.
            assert (first.x, first.y) == (0.0, 5.0)
            assert (second.x, second.y) == (0.0, 7.0)
        assert _svg(resolved) == _svg(request)


@pytest.mark.linear
def test_stored_translations_without_an_anchor_draw_as_stored() -> None:
    # A 0.14.5 Session stores the sideways moves of an anchored figure, never the option: it draws as it was.
    planned = plan_linear_request(_request()).request
    keys = tuple(str(record.record_key) for record in planned.records)
    held = replace(planned, layout=LinearMultiRecordOptions(record_translations=(
        LinearRecordTranslation(record_key=keys[0]),
        LinearRecordTranslation(record_key=keys[1], x=-62.8),
    )))
    assert _record_shifts(_svg(held)) == [0.0, -62.8]
    assert _starts(plan_linear_request(held)) == (None, None)


@pytest.mark.linear
def test_anchor_feature_is_the_last_positional_option_field() -> None:
    assert [item.name for item in fields(LinearDiagramOptions)][-1] == "anchor_feature"
    assert [item.name for item in fields(LinearOptions)][-1] == "anchor_feature"
    assert LinearDiagramOptions(anchor_feature="gene=COX1").anchor_feature == _COX1


def _human_and_zebrafish(*, first_display: object = None, **options: object) -> LinearDiagramRequest:
    first = RecordInput(GenBankInputSource(_HUMAN))
    if first_display is not None:
        assert isinstance(first_display, RecordDisplayOptions)
        first = replace(first, display=first_display)
    return LinearDiagramRequest(
        records=(first, RecordInput(GenBankInputSource(_ZEBRAFISH))),
        options=LinearDiagramOptions(**options),
    )


@pytest.mark.linear
def test_only_matches_that_start_in_the_drawn_region_count(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # Each record has two dUTPase CDSs (MjeNMV 249324..249569 and
    # 250088..250528, MelaMJNV 234967..235212 and 235538..235978); each region
    # draws only the first, so the selector names one feature in each record.
    with caplog.at_level(logging.WARNING):
        linear_cli.linear_main([
            "--gbk", str(_EXAMPLES / "MjeNMV.gb"), str(_EXAMPLES / "MelaMJNV.gb"),
            "--region", "LC738868.1:249000-249800", "--region", "LC738874.1:234800-235400",
            "--anchor_feature", "product=dUTPase",
            "-f", "svg", "-o", str(tmp_path / "region"),
        ])
    svg = (tmp_path / "region.svg").read_text(encoding="utf-8")
    # 249324 is 324 bp into the first region and 234967 is 167 bp into the second.
    assert _record_shifts(svg)[0] == 0.0
    assert _record_shifts(svg)[1] > 0.0
    assert "matches 2 features at different starts" not in caplog.text
    assert "keeps its place" not in caplog.text


@pytest.mark.linear
def test_a_match_that_starts_outside_the_drawn_region_is_reported_as_such(
    caplog: pytest.LogCaptureFixture,
) -> None:
    from gbdraw.io.regions import parse_region_spec

    request = _human_and_zebrafish(anchor_feature=_COX1)
    cropped = replace(request, records=(
        replace(request.records[0], region=parse_region_spec("NC_012920.1:6000-8000")),
        request.records[1],
    ))
    with caplog.at_level(logging.WARNING):
        plan = plan_linear_request(cropped)
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.statuses == ("unmatched", "rotated")
    assert (
        "anchor_feature gene=COX1 matched 2 feature(s) in record #1 'NC_012920.1' whose start "
        "lies outside its drawn region; the record keeps its place."
    ) in caplog.text


@pytest.mark.linear
def test_a_feature_across_a_rotated_display_start_starts_before_the_seam() -> None:
    # Rotated to 7000, human COX1 (5904..7445) is drawn as 15473..16569 then 0..446:
    # the anchor rotates the record to the start of the piece before the seam, base 5904.
    plan = plan_linear_request(_human_and_zebrafish(
        first_display=RecordDisplayOptions(is_circular=True, start_coordinate=7000),
        anchor_feature=_COX1,
    ))
    assert _starts(plan) == (5904, 6425)


@pytest.mark.linear
def test_a_feature_across_the_origin_starts_on_its_piece_before_it() -> None:
    # The human D-loop is complement(join(16024..16569,1..576)) on a circular record.
    selector = FeatureSelector(value="D-loop", key="type")
    plan = plan_linear_request(_drawn_linear(_request(anchor_feature=selector)))
    assert plan.feature_anchor is not None
    assert plan.feature_anchor.starts == (16023.0, 0.0)
    assert _starts(plan_linear_request(_request(anchor_feature=selector))) == (16024, 1)


@pytest.mark.linear
def test_records_that_share_a_row_are_refused() -> None:
    request = _request(anchor_feature=_COX1)
    shared = replace(
        request,
        records=tuple(
            replace(record, presentation=RecordPresentation(grid_row=1)) for record in request.records
        ),
        layout=LinearMultiRecordOptions(),
    )
    with pytest.raises(ValidationError, match="cannot align records that share a row") as caught:
        plan_linear_request(shared)
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "anchor_feature", "reason": "ANCHOR_SHARED_ROW"}


_WRAP = "canvas.linear.wrap_length"


@pytest.mark.linear
def test_records_the_anchor_only_rotates_wrap() -> None:
    wrapped = _svg(_request(anchor_feature=_COX1, config_overrides={_WRAP: 5000}))
    assert "data-gbdraw-wrap-row" in wrapped
    explicit = _svg(LinearDiagramRequest(
        records=tuple(
            RecordInput(GenBankInputSource(path), display=RecordDisplayOptions(start_coordinate=start))
            for path, start in zip((_HUMAN, _ZEBRAFISH), (5904, 6425))
        ),
        options=LinearDiagramOptions(config_overrides={_WRAP: 5000}),
    ))
    assert wrapped == explicit


@pytest.mark.linear
@pytest.mark.parametrize("linear", [(0,), (1,), (0, 1)])
def test_a_wrap_that_would_move_a_cut_record_is_refused(linear: tuple[int, ...]) -> None:
    # (0,): the cut human record moves to the rotated zebrafish record's left end;
    # (1,): the cut zebrafish record moves to the rotated human one's;
    # (0, 1): nothing is rotated, and the records line up with the cut human record's COX1.
    request = _drawn_linear(_request(anchor_feature=_COX1, config_overrides={_WRAP: 5000}), *linear)
    with pytest.raises(ValidationError, match="anchor_feature cannot align these records") as caught:
        build_request_diagram(request)
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "wrap_length", "reason": "WRAP_ANCHOR_FEATURE"}
    # A wrap that cuts nothing changes nothing.
    uncut = replace(request, options=replace(request.options, config_overrides={_WRAP: 20000}))
    assert _svg(uncut) == _svg(replace(request, options=replace(request.options, config_overrides=None)))


@pytest.mark.linear
def test_the_anchor_rotates_a_cut_record_a_stored_translation_moved() -> None:
    # OV-2496: a 0.14.5 anchored Session moved the zebrafish record sideways; aligned again under a wrap, the
    # anchor rotates it to x 0 and it wraps, whatever its stored x was.
    planned = plan_linear_request(_request()).request
    keys = tuple(str(record.record_key) for record in planned.records)
    held = replace(
        planned,
        options=replace(planned.options, anchor_feature=_COX1, config_overrides={_WRAP: 5000}),
        layout=LinearMultiRecordOptions(record_translations=(
            LinearRecordTranslation(record_key=keys[0]),
            LinearRecordTranslation(record_key=keys[1], x=-62.8),
        )),
    )
    svg = _svg(held)
    assert "data-gbdraw-wrap-row" in svg
    assert _record_shifts(svg) == [0.0, 0.0]


@pytest.mark.linear
def test_a_moved_record_the_wrap_does_not_cut_draws() -> None:
    # The human record is cut and rotated, at x 0; the cropped zebrafish record is not cut, and moves.
    from gbdraw.io.regions import parse_region_spec

    request = _request(anchor_feature=_COX1, config_overrides={_WRAP: 10000})
    cropped = replace(request, records=(
        request.records[0],
        replace(request.records[1], region=parse_region_spec("NC_002333.2:6000-9000")),
    ))
    svg = _svg(cropped)
    assert "data-gbdraw-wrap-row" in svg
    human, zebrafish = _record_shifts(svg)
    assert human == 0.0
    assert zebrafish == pytest.approx(-425.0 * 2000 / 10000)


@pytest.mark.linear
def test_the_anchor_is_read_after_auto_orient() -> None:
    import pandas as pd

    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    from gbdraw.linear_comparison import LinearComparison

    # One reverse hit: auto_orient reverse-complements the zebrafish record.
    reverse = pd.DataFrame(
        [("Q", "S", 99.0, 600, 0, 0, 1, 600, 6000, 5401, 1e-20, 50.0)], columns=list(COMPARISON_COLUMNS)
    )
    comparisons = (LinearComparison(0, 1, reverse),)
    oriented = plan_linear_request(_request(anchor_feature=_COX1, auto_orient=True, linear_comparisons=comparisons))
    request = _request(anchor_feature=_COX1, linear_comparisons=comparisons)
    explicit = plan_linear_request(replace(request, records=(
        request.records[0],
        replace(request.records[1], presentation=RecordPresentation(reverse_complement=True)),
    )))
    assert _starts(oriented) == _starts(explicit) == (5904, 7975)


def _join_records(*, reverse_first: bool, strand: int | None = 1) -> LinearDiagramRequest:
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api import InMemoryRecordSource

    def record(name: str, location: SimpleLocation | CompoundLocation) -> SeqRecord:
        result = SeqRecord(Seq("ACGT" * 25), id=name, annotations={"molecule_type": "DNA", "topology": "circular"})
        result.features = [SeqFeature(location, type="gene", qualifiers={"gene": ["x"]})]
        return result

    # join(11..13,71..73): its path runs from 11 across the wider gap to 73.
    first = record("first", CompoundLocation([
        SimpleLocation(10, 13, strand=strand), SimpleLocation(70, 73, strand=strand),
    ]))
    second = record("second", SimpleLocation(30, 40, strand=1))
    return LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), presentation=RecordPresentation(reverse_complement=reverse_first)),
            RecordInput(InMemoryRecordSource(second)),
        ),
        options=LinearDiagramOptions(anchor_feature="gene=x"),
    )


@pytest.mark.linear
def test_a_join_starts_where_its_path_begins_as_drawn() -> None:
    # Forward, the path begins at its 5' base 11 (offset 10), not at 71 after the wider gap.
    assert _starts(plan_linear_request(_join_records(reverse_first=False))) == (11, 31)
    plan = plan_linear_request(_drawn_linear(_join_records(reverse_first=False)))
    assert plan.feature_anchor is not None and plan.feature_anchor.starts == (10.0, 30.0)


@pytest.mark.linear
def test_a_reversed_join_starts_at_its_3_prime_base() -> None:
    # Reversed, the path's 3' base 73 is drawn first: offset 100 - 73.
    assert _starts(plan_linear_request(_join_records(reverse_first=True))) == (73, 31)
    plan = plan_linear_request(_drawn_linear(_join_records(reverse_first=True)))
    assert plan.feature_anchor is not None and plan.feature_anchor.starts == (27.0, 30.0)


@pytest.mark.linear
def test_an_unstranded_feature_rotates_as_a_plus_strand_one() -> None:
    # The Rotate popup refuses a feature without a strand; the anchor reads it as a plus-strand one.
    assert _starts(plan_linear_request(_join_records(reverse_first=False, strand=None))) == (11, 31)
