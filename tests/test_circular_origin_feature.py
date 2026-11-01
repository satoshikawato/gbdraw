"""Rotate each Circular record so a named feature starts at 12 o'clock (``--origin_feature``).

The selector (``key=value``, exact, the feature-selector form of the
annotation and placement tables) is resolved on each record, and the
feature's start as drawn becomes the record's display start
(``RecordDisplayOptions.start_coordinate``). The figure equals the one drawn
with that start given explicitly, and a saved Session stores the start, not
the option. The start base is a twin of the Web Rotate to feature "Start of
the record": both read feature-start-anchor-vectors.json.
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
from dataclasses import fields, replace
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw import CircularOptions, draw_circular, read_genbank
from gbdraw import circular as circular_cli
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    FeatureSelector,
    GenBankInputSource,
    InMemoryRecordSource,
    RecordDisplayOptions,
    RecordInput,
    plan_circular_request,
)
from gbdraw.api.record_planning import FEATURE_START_REFUSALS, resolve_origin_feature_start
from gbdraw.api.request_render import plan_circular_batch_request
from gbdraw.api.requests import CircularBatchOutputPolicy, CircularBatchRequest
from gbdraw.core.record_metadata import _SourceFeatureAnchorProfile
from gbdraw.exceptions import ValidationError
from gbdraw.features.ids import SourceFeatureHash
from gbdraw.features.source import SourceFeatureIdentity
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.session_io import load_session

_ROOT = Path(__file__).resolve().parents[1]
_TOBACCO = _ROOT / "tests" / "test_inputs" / "NC_001879.gbk"
_HUMAN = _ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
_PSBA = FeatureSelector(value="psbA", key="gene")
# psbA is complement(536..1597): drawn forward, its 3' base 536 comes first.
_PSBA_START = 536
_VECTORS = _ROOT / "tests" / "fixtures" / "feature-start-anchor-vectors.json"


def _request(*inputs: RecordInput, **options: object) -> CircularDiagramRequest:
    return CircularDiagramRequest(
        records=inputs or (RecordInput(GenBankInputSource(_TOBACCO)),),
        options=CircularDiagramOptions(**options),
    )


def _cli_svg(tmp_path: Path, name: str, *args: str) -> str:
    circular_cli.circular_main(["--gbk", str(_TOBACCO), "-f", "svg", "-o", str(tmp_path / name), *args])
    return (tmp_path / f"{name}.svg").read_text(encoding="utf-8")


@pytest.mark.circular
def test_cli_origin_feature_draws_the_display_start_figure(tmp_path: Path) -> None:
    rotated = _cli_svg(tmp_path, "origin", "--origin_feature", "gene=psbA")
    assert rotated == _cli_svg(tmp_path, "start", "--display_start_coordinate", str(_PSBA_START))
    assert rotated != _cli_svg(tmp_path, "plain")


@pytest.mark.circular
def test_the_plan_writes_the_start_into_the_record_display() -> None:
    plan = plan_circular_request(_request(origin_feature="gene=psbA"))
    assert plan.provenance[0].display == RecordDisplayOptions(start_coordinate=_PSBA_START)
    assert plan.transforms[0].start_coordinate == _PSBA_START
    assert plan.request.records[0].display.start_coordinate == _PSBA_START
    # The planned request carries the start, not the option.
    assert plan.request.options.origin_feature is None


@pytest.mark.circular
def test_a_reversed_record_starts_at_the_features_high_coordinate() -> None:
    from gbdraw.api.request_render import build_request_diagram
    from gbdraw.api.requests import RecordPresentation

    reversed_input = RecordInput(GenBankInputSource(_TOBACCO), presentation=RecordPresentation(reverse_complement=True))
    rotated = _request(reversed_input, origin_feature=_PSBA)
    # Drawn reversed, psbA reads forward from its high coordinate 1597.
    assert plan_circular_request(rotated).provenance[0].display.start_coordinate == 1597
    explicit = _request(replace(reversed_input, display=RecordDisplayOptions(start_coordinate=1597)))
    assert build_request_diagram(rotated).drawing.tostring() == build_request_diagram(explicit).drawing.tostring()


@pytest.mark.circular
def test_python_origin_feature_draws_the_display_start_figure() -> None:
    records = read_genbank(str(_TOBACCO))
    rotated = draw_circular(records, options=CircularOptions(origin_feature="gene=psbA")).to_svg()
    explicit = draw_circular(
        records, record_displays=[RecordDisplayOptions(start_coordinate=_PSBA_START)]
    ).to_svg()
    assert rotated == explicit


@pytest.mark.circular
def test_origin_feature_is_the_last_positional_option_field() -> None:
    assert [item.name for item in fields(CircularDiagramOptions)][-1] == "origin_feature"
    assert [item.name for item in fields(CircularOptions)][-1] == "origin_feature"
    assert CircularDiagramOptions(origin_feature="gene=psbA").origin_feature == _PSBA


def _refused(request: CircularDiagramRequest, match: str) -> ValidationError:
    with pytest.raises(ValidationError, match=match) as raised:
        plan_circular_request(request)
    assert raised.value.diagnostic is not None
    return raised.value


@pytest.mark.circular
def test_origin_feature_with_a_display_start_is_refused() -> None:
    _refused(
        _request(
            RecordInput(GenBankInputSource(_TOBACCO), display=RecordDisplayOptions(start_coordinate=10)),
            origin_feature=_PSBA,
        ),
        "origin_feature gene=psbA cannot be combined with a display start",
    )


@pytest.mark.circular
def test_origin_feature_on_a_cropped_record_is_refused() -> None:
    from gbdraw.io.regions import parse_region_spec

    _refused(
        _request(
            RecordInput(GenBankInputSource(_TOBACCO), region=parse_region_spec("NC_001879.2:1-5000")),
            origin_feature=_PSBA,
        ),
        "origin_feature gene=psbA cannot rotate the cropped record",
    )


@pytest.mark.circular
def test_origin_feature_on_a_linear_record_is_refused() -> None:
    _refused(
        _request(
            RecordInput(GenBankInputSource(_TOBACCO), display=RecordDisplayOptions(is_circular=False)),
            origin_feature=_PSBA,
        ),
        "origin_feature gene=psbA cannot rotate the linear record #1",
    )


@pytest.mark.circular
def test_origin_feature_that_matches_no_record_is_refused() -> None:
    _refused(_request(origin_feature="gene=noSuchGene"), "origin_feature gene=noSuchGene matched no feature in any record")


@pytest.mark.circular
def test_origin_feature_that_matches_features_at_different_starts_is_refused() -> None:
    # ndhB lies in both inverted repeats (gene, CDS, exons and intron of each copy).
    error = _refused(_request(origin_feature="gene=ndhB"), "origin_feature gene=ndhB matches 10 features at different starts")
    assert "complement(join(" in str(error)


def _strand_record(location: SimpleLocation | CompoundLocation) -> RecordInput:
    record = SeqRecord(Seq("ACGT" * 25), id="strand", annotations={"molecule_type": "DNA", "topology": "circular"})
    record.features = [SeqFeature(location, type="misc_feature", qualifiers={"note": ["target"]})]
    return RecordInput(InMemoryRecordSource(record))


@pytest.mark.circular
@pytest.mark.parametrize("location", [
    SimpleLocation(10, 16),
    CompoundLocation([SimpleLocation(10, 13, strand=1), SimpleLocation(40, 43, strand=-1)]),
])
def test_origin_feature_without_one_strand_is_refused(location: SimpleLocation | CompoundLocation) -> None:
    _refused(_request(_strand_record(location), origin_feature="note=target"), "origin_feature note=target in record #1")


@pytest.mark.circular
def test_a_record_without_the_feature_keeps_its_start(caplog: pytest.LogCaptureFixture) -> None:
    request = CircularBatchRequest(
        records=(RecordInput(GenBankInputSource(_TOBACCO)), RecordInput(GenBankInputSource(_HUMAN))),
        options=CircularDiagramOptions(origin_feature=_PSBA),
        output_policy=CircularBatchOutputPolicy(output_prefix="unused"),
    )
    with caplog.at_level(logging.WARNING):
        plan = plan_circular_batch_request(request)
    assert [item.display.start_coordinate for item in plan.provenance] == [_PSBA_START, None]
    assert plan.request.options.origin_feature is None
    assert (
        "origin_feature gene=psbA matched no feature in record #2 'NC_012920.1'; the record keeps its start."
    ) in caplog.text


@pytest.mark.circular
def test_a_session_stores_the_start_and_replays_the_figure(tmp_path: Path) -> None:
    from gbdraw.session_request_codec import encode_canonical_request

    encoded = json.dumps(encode_canonical_request(plan_circular_request(_request(origin_feature=_PSBA)).request).payload)
    assert "originFeature" not in encoded and "origin_feature" not in encoded
    session = tmp_path / "saved.gbdraw-session.json"
    rotated = _cli_svg(tmp_path, "origin", "--origin_feature", "gene=psbA", "--session_output", str(session))
    stored = load_session(session)["drawings"][0]["renderRequest"]
    assert stored["records"][0]["display"]["startCoordinate"] == _PSBA_START
    assert "originFeature" not in stored["diagramOptions"]
    completed = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "circular",
         "--session", str(session), "-f", "svg", "-o", str(tmp_path / "replay")],
        capture_output=True, text=True, check=False,
        env={**os.environ, "PYTHONPATH": str(_ROOT)},
    )
    assert completed.returncode == 0, completed.stderr
    assert (tmp_path / "replay.svg").read_text(encoding="utf-8") == rotated


_STRANDS = {"+": 1, "-": -1}


def _vector_feature(case: dict[str, object]) -> SourceFeatureIdentity:
    parts = case["parts"]
    assert isinstance(parts, list)
    return SourceFeatureIdentity(
        biological_feature_id="feature-1",
        source_feature_index=0,
        feature_type="misc_feature",
        location_parts=tuple((part["start"], part["end"], _STRANDS.get(part["strand"])) for part in parts),
        qualifiers=(),
        stable_feature_id=SourceFeatureHash("feature-1"),
    )


def _vector_profile(case: dict[str, object]) -> _SourceFeatureAnchorProfile:
    profile = case["profile"]
    assert isinstance(profile, dict)
    return _SourceFeatureAnchorProfile(
        precision=profile["precision"], operator=profile["operator"],
        part_order=profile["partOrder"], strand=profile["strand"],
    )


@pytest.mark.parametrize(
    "case", json.loads(_VECTORS.read_text(encoding="utf-8"))["cases"], ids=lambda case: case["name"]
)
def test_the_start_base_agrees_with_the_web_vector(case: dict[str, object]) -> None:
    length = case["recordLength"]
    assert isinstance(length, int)
    transform = (
        RecordDisplayTransform(length, length, -1, None, True)
        if case["reverseComplement"]
        else RecordDisplayTransform(length, 1, 1, None, True)
    )
    feature, profile = _vector_feature(case), _vector_profile(case)
    if "refused" in case:
        with pytest.raises(ValidationError) as raised:
            resolve_origin_feature_start(feature, transform, profile)
        assert raised.value.diagnostic is not None
        # The Web refuses it with this code, for the same reason.
        assert f": {FEATURE_START_REFUSALS[case['refused']]}." in str(raised.value)
    else:
        assert resolve_origin_feature_start(feature, transform, profile) == case["startCoordinate"]


@pytest.mark.circular
def test_an_unplanned_request_with_origin_feature_is_not_encoded() -> None:
    from gbdraw.session_request_codec import CanonicalRequestEncodingError, encode_canonical_request

    with pytest.raises(CanonicalRequestEncodingError, match="origin_feature record rotation"):
        encode_canonical_request(_request(origin_feature=_PSBA))


@pytest.mark.parametrize("selector", ["gene=", "", "=x"])
def test_an_empty_selector_is_refused_with_a_diagnostic(selector: str) -> None:
    with pytest.raises(ValidationError, match="Feature selector") as raised:
        CircularDiagramOptions(origin_feature=selector)
    assert raised.value.diagnostic == {"code": "INPUT_INVALID", "reason": "REQUIRED"}
