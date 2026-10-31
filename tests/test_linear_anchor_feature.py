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


@pytest.mark.linear
def test_a_session_does_not_store_anchor_feature_yet(tmp_path: Path) -> None:
    from gbdraw.session_request_codec import encode_canonical_request

    with pytest.raises(ValidationError, match="anchor_feature"):
        encode_canonical_request(_request(anchor_feature=_COX1))
    out = tmp_path / "out"
    out.mkdir()
    completed = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear",
         "--gbk", str(_HUMAN), str(_ZEBRAFISH), "-f", "svg",
         "-o", str(out / "saved"), "--anchor_feature", "gene=COX1", "--save_session"],
        capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
    )
    assert completed.returncode == 1
    assert "A Session does not store anchor_feature yet; save the Session without it." in completed.stderr
    assert list(out.iterdir()) == []


@pytest.mark.linear
def test_anchor_feature_is_the_last_positional_option_field() -> None:
    assert [item.name for item in fields(LinearDiagramOptions)][-1] == "anchor_feature"
    assert [item.name for item in fields(LinearOptions)][-1] == "anchor_feature"
    assert LinearDiagramOptions(anchor_feature="gene=COX1").anchor_feature == _COX1
