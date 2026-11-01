"""Align Linear records at the start of one feature (``--anchor_feature``).

The selector (``key=value``, exact, the feature-selector form of the
annotation and placement tables) is resolved on every record as drawn. The
first record that has the feature keeps its place; every other record that has
it moves so that the feature's left end, as drawn, lies at the same x. The
figure equals the one drawn with the same moves given as explicit record
translations.
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
from gbdraw.api.request_render import build_request_diagram
from gbdraw.api.requests import RecordPresentation
from gbdraw.exceptions import ValidationError

_DATA = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data"
_HUMAN = _DATA / "human-mitochondrion" / "HmmtDNA.gbk"
_ZEBRAFISH = _DATA / "metazoan-mitochondria-four" / "NC_002333.2.gb"
_EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
_COX1 = FeatureSelector(value="COX1", key="gene")
# Human COX1 is 5904..7445 and zebrafish COX1 6425..7975 (gene and CDS alike).
_HUMAN_COX1_START = 5903.0
_ZEBRAFISH_COX1_START = 6424.0


def _request(**options: object) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(_HUMAN)), RecordInput(GenBankInputSource(_ZEBRAFISH))),
        options=LinearDiagramOptions(**options),
    )


def _svg(request: LinearDiagramRequest) -> str:
    return build_request_diagram(request).drawing.tostring()


def _record_shifts(svg: str) -> list[float]:
    """Each record group's own x translation, in record order."""

    return [float(match) for match in re.findall(r'data-record-translation-x="([-0-9.e]+)"', svg)]


@pytest.mark.linear
def test_anchor_feature_aligns_the_feature_starts_like_explicit_translations() -> None:
    plan = plan_linear_request(_request(anchor_feature=_COX1))
    assert plan.alignment_anchor_centers == (_HUMAN_COX1_START, _ZEBRAFISH_COX1_START)

    anchored = _svg(_request(anchor_feature=_COX1))
    plain = _svg(_request())
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
def test_anchor_feature_is_resolved_on_the_record_as_drawn() -> None:
    request = _request(anchor_feature=_COX1)
    reversed_second = replace(request, records=(
        request.records[0],
        replace(request.records[1], presentation=RecordPresentation(reverse_complement=True)),
    ))
    plan = plan_linear_request(reversed_second)
    # Reversed, zebrafish COX1 (6425..7975 of 16596) starts 16596 - 7975 bp from the left.
    assert plan.alignment_anchor_centers == (_HUMAN_COX1_START, 16596.0 - 7975.0)


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
    assert plan.alignment_anchor_centers == (_HUMAN_COX1_START, None, _ZEBRAFISH_COX1_START)
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
    centered = _svg(_request(config_overrides={"canvas.linear.align_center": True}))
    anchored = _svg(_request(anchor_feature=_COX1, config_overrides={"canvas.linear.align_center": True}))
    plain_anchored = _svg(_request(anchor_feature=_COX1))
    # The first record keeps its centred place; the second record moves from
    # its own centred place (27 bp right of the first, half of the 27 bp
    # length difference) so that its COX1 start sits under the first one's.
    assert _record_shifts(centered) == [0.0, 0.0]
    first, shift = _record_shifts(anchored)
    assert first == 0.0
    assert shift == pytest.approx(_record_shifts(plain_anchored)[1] + 13.5 * 2000 / 16596)
    assert anchored != plain_anchored


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
    assert _record_shifts(cli) == _record_shifts(_svg(_request(anchor_feature=_COX1)))
    assert _record_shifts(plain) == [0.0, 0.0]


@pytest.mark.linear
def test_python_anchor_feature_draws_the_typed_figure() -> None:
    records = read_genbank([str(_HUMAN), str(_ZEBRAFISH)])
    anchored = draw_linear(records, options=LinearOptions(anchor_feature="gene=COX1"))
    plain = draw_linear(records)
    assert _record_shifts(anchored.to_svg()) == _record_shifts(_svg(_request(anchor_feature=_COX1)))
    assert _record_shifts(plain.to_svg()) == [0.0, 0.0]


def _cli(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", *args],
        cwd=cwd, capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
    )


def _stored_layout(session_path: Path) -> dict[str, object]:
    import json

    request = json.loads(session_path.read_text(encoding="utf-8"))["renderRequest"]
    assert "anchorFeature" not in request["diagramOptions"]
    return request["layout"]


@pytest.mark.linear
def test_a_session_stores_the_translations_the_anchor_resolves_to(tmp_path: Path) -> None:
    saved = _cli(
        "--gbk", str(_HUMAN), str(_ZEBRAFISH), "-f", "svg", "-o", "saved",
        "--anchor_feature", "gene=COX1", "--save_session", cwd=tmp_path,
    )
    assert saved.returncode == 0, saved.stderr
    layout = _stored_layout(tmp_path / "saved.gbdraw-session.json")
    drawn = (tmp_path / "saved.svg").read_text(encoding="utf-8")
    assert [entry["x"] for entry in layout["recordTranslations"]] == _record_shifts(drawn)
    assert _record_shifts(drawn)[1] < 0.0
    replayed = _cli("--session", "saved.gbdraw-session.json", "-o", "again", cwd=tmp_path)
    assert replayed.returncode == 0, replayed.stderr
    assert (tmp_path / "again.svg").read_bytes() == (tmp_path / "saved.svg").read_bytes()


@pytest.mark.linear
def test_the_python_session_writers_store_the_translations_the_anchor_resolves_to(tmp_path: Path) -> None:
    from gbdraw.session import (
        build_session_document,
        load_session_document,
        materialize_session,
        save_session_document,
        session_to_request,
    )

    request = _request(anchor_feature=_COX1)
    drawn = _svg(request)
    built = build_session_document(request).to_dict()
    path = tmp_path / "saved.gbdraw-session.json"
    save_session_document(path, request)
    assert _stored_layout(path) == built["renderRequest"]["layout"]
    assert [entry["x"] for entry in _stored_layout(path)["recordTranslations"]] == _record_shifts(drawn)
    with materialize_session(load_session_document(path), output_directory=tmp_path / "out") as materialized:
        replayed = session_to_request(materialized)
        assert isinstance(replayed, LinearDiagramRequest) and replayed.options.anchor_feature is None
        assert _svg(replayed) == drawn


@pytest.mark.linear
@pytest.mark.parametrize("config_overrides", [None, {"canvas.linear.align_center": True}])
def test_the_resolved_request_draws_the_anchored_figure(config_overrides: dict[str, object] | None) -> None:
    from gbdraw.api import resolve_request
    from gbdraw.session_request_codec import encode_canonical_request

    request = _request(anchor_feature=_COX1, config_overrides=config_overrides)
    # The request itself is not encoded: the anchor is resolved when the figure is laid out.
    with pytest.raises(ValidationError, match="anchor_feature record placement"):
        encode_canonical_request(request)
    resolved = resolve_request(request)
    assert isinstance(resolved, LinearDiagramRequest)
    assert resolved.options.anchor_feature is None
    layout = encode_canonical_request(resolved).payload["layout"]
    assert [entry["x"] for entry in layout["recordTranslations"]] == _record_shifts(_svg(request))
    assert _svg(resolved) == _svg(request)


@pytest.mark.linear
def test_the_anchor_replaces_the_stored_x_of_the_records_it_moves() -> None:
    from gbdraw.api import resolve_request

    planned = plan_linear_request(_request(anchor_feature=_COX1)).request
    keys = tuple(str(record.record_key) for record in planned.records)
    held = LinearMultiRecordOptions(record_translations=(
        LinearRecordTranslation(record_key=keys[0], x=40.0, y=5.0),
        LinearRecordTranslation(record_key=keys[1], x=-300.0, y=7.0),
    ))
    request = replace(planned, layout=held)
    resolved = resolve_request(request)
    assert isinstance(resolved, LinearDiagramRequest) and resolved.layout is not None
    first, second = resolved.layout.record_translations
    shift = -(_ZEBRAFISH_COX1_START - _HUMAN_COX1_START) * 2000 / 16596
    # The first record keeps its own translation; the moved record's x follows it, its y stays.
    assert (first.x, first.y) == (40.0, 5.0)
    assert second.x == pytest.approx(40.0 + shift) and second.y == 7.0
    assert _svg(resolved) == _svg(request)


@pytest.mark.linear
def test_anchor_feature_is_the_last_positional_option_field() -> None:
    assert [item.name for item in fields(LinearDiagramOptions)][-1] == "anchor_feature"
    assert [item.name for item in fields(LinearOptions)][-1] == "anchor_feature"
    assert LinearDiagramOptions(anchor_feature="gene=COX1").anchor_feature == _COX1


def _human_and_zebrafish(*, first_display: object = None, **options: object) -> LinearDiagramRequest:
    from gbdraw.api.requests import RecordDisplayOptions

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
    assert plan.alignment_anchor_centers == (None, _ZEBRAFISH_COX1_START)
    assert (
        "anchor_feature gene=COX1 matched 2 feature(s) in record #1 'NC_012920.1' whose start "
        "lies outside its drawn region; the record keeps its place."
    ) in caplog.text


@pytest.mark.linear
def test_a_feature_across_a_rotated_display_start_starts_before_the_seam() -> None:
    from gbdraw.api.requests import RecordDisplayOptions

    # Rotated to 7000, human COX1 (5904..7445) is drawn as 15473..16569 then 0..446.
    plan = plan_linear_request(_human_and_zebrafish(
        first_display=RecordDisplayOptions(is_circular=True, start_coordinate=7000),
        anchor_feature=_COX1,
    ))
    assert plan.alignment_anchor_centers == (15473.0, _ZEBRAFISH_COX1_START)


@pytest.mark.linear
def test_a_feature_across_the_origin_starts_on_its_piece_before_it() -> None:
    # The human D-loop is complement(join(16024..16569,1..576)) on a circular record.
    plan = plan_linear_request(_request(anchor_feature=FeatureSelector(value="D-loop", key="type")))
    assert plan.alignment_anchor_centers == (16023.0, 0.0)


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


@pytest.mark.linear
def test_a_wrap_that_cuts_a_record_is_refused_by_the_wrap_owner() -> None:
    with pytest.raises(ValidationError, match="anchor_feature cannot align a record cut by wrap_length") as caught:
        plan_linear_request(_request(anchor_feature=_COX1, config_overrides={"canvas.linear.wrap_length": 5000}))
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "wrap_length", "reason": "WRAP_ANCHOR_FEATURE"}
    # A wrap that cuts nothing changes nothing.
    uncut = _svg(_request(anchor_feature=_COX1, config_overrides={"canvas.linear.wrap_length": 20000}))
    assert _record_shifts(uncut) == _record_shifts(_svg(_request(anchor_feature=_COX1)))


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
    assert oriented.alignment_anchor_centers == explicit.alignment_anchor_centers == (
        _HUMAN_COX1_START, 16596.0 - 7975.0,
    )


def _join_records(*, reverse_first: bool) -> LinearDiagramRequest:
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api import InMemoryRecordSource

    def record(name: str, location: SimpleLocation | CompoundLocation) -> SeqRecord:
        result = SeqRecord(Seq("ACGT" * 25), id=name, annotations={"molecule_type": "DNA", "topology": "circular"})
        result.features = [SeqFeature(location, type="gene", qualifiers={"gene": ["x"]})]
        return result

    # join(11..13,71..73): its path runs from 11 across the wider gap to 73.
    first = record("first", CompoundLocation([SimpleLocation(10, 13, strand=1), SimpleLocation(70, 73, strand=1)]))
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
    assert plan_linear_request(_join_records(reverse_first=False)).alignment_anchor_centers == (10.0, 30.0)


@pytest.mark.linear
def test_a_reversed_join_starts_at_its_3_prime_base() -> None:
    # Reversed, the path's 3' base 73 is drawn first: offset 100 - 73.
    assert plan_linear_request(_join_records(reverse_first=True)).alignment_anchor_centers == (27.0, 30.0)
