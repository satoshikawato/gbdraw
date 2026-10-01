"""D-26 (PD-OI-080): ACGTU dinucleotide pairs, with U counted as T."""

from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.testing import assert_frame_equal

from gbdraw.analysis.gc import circular_dinucleotide_content_df
from gbdraw.analysis.skew import skew_df
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.tracks.parsing import slot_dinucleotide
from gbdraw.exceptions import ValidationError
from gbdraw.mode_profiles import validate_dinucleotide
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    parse_circular_track_slots,
    parse_linear_track_slots,
)
from gbdraw.web_support.error_adapter import FIELDS, serialize_web_error

DNA = "ATTAGCATTTAAGCGCATATTTTAAAGGCCATAT" * 20


def _record(sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id="r")


DIAGNOSTIC = {"code": "INPUT_INVALID", "field": "dinucleotide", "reason": "DINUCLEOTIDE"}


def _assert_dinucleotide_rejection(error: BaseException) -> None:
    """The producer owns the meaning; the Web sees the published identifiers."""

    cause = error
    while getattr(cause, "diagnostic", None) is None and cause.__cause__ is not None:
        cause = cause.__cause__
    assert cause.diagnostic == DIAGNOSTIC
    payload = serialize_web_error(error, operation="generate", stage="request-validation")
    assert payload["code"] == "INPUT_INVALID"
    # ``dinucleotide`` joins the published field vocabulary with its Web label.
    expected = {"reason": "DINUCLEOTIDE"}
    if "dinucleotide" in FIELDS:
        expected["field"] = "dinucleotide"
    assert payload["context"] == expected


@pytest.mark.parametrize(
    ("value", "canonical"),
    [("GC", "GC"), ("at", "AT"), ("Au", "AU"), ("uU", "UU"), (" ta ", "TA")],
)
def test_validator_accepts_acgtu_pairs_case_insensitively(value, canonical):
    assert validate_dinucleotide(value) == canonical


@pytest.mark.parametrize("value", ["G", "GCA", "XY", "GN", "", None, 12, "A-"])
def test_validator_rejects_other_values_with_one_reason(value):
    with pytest.raises(ValidationError) as caught:
        validate_dinucleotide(value)
    _assert_dinucleotide_rejection(caught.value)


def test_u_is_the_same_base_as_t_for_skew_and_content():
    at = skew_df(_record(DNA), 50, 10, "AT")
    au = skew_df(_record(DNA), 50, 10, "AU")
    # The display name keeps the typed letters; the values are the AT values.
    assert list(au.columns) == ["AU content", "AU skew", "Cumulative AU skew, normalized"]
    assert_frame_equal(at.set_axis(au.columns, axis=1), au)
    rna = skew_df(_record(DNA.replace("T", "U")), 50, 10, "AT")
    assert_frame_equal(at, rna)
    content_at = circular_dinucleotide_content_df(_record(DNA), 50, 10, "AT")
    content_au = circular_dinucleotide_content_df(_record(DNA.replace("T", "U")), 50, 10, "AU")
    assert content_au.columns.tolist() == ["AU content"]
    assert content_au["AU content"].tolist() == content_at["AT content"].tolist()


@pytest.mark.parametrize("options_type", [CircularDiagramOptions, LinearDiagramOptions])
def test_typed_options_normalize_and_reject_dinucleotide(options_type):
    assert options_type(dinucleotide="au").dinucleotide == "AU"
    with pytest.raises(ValidationError):
        options_type(dinucleotide="XY")


def test_circular_slot_nt_uses_the_same_validator():
    slots = parse_circular_track_slots(["at:dinucleotide_skew@nt=au"])
    assert slots[0].params["nt"] == "AU"
    for bad in ("at:dinucleotide_skew@nt=G", "at:dinucleotide_content@nt=XY"):
        with pytest.raises(ValueError) as caught:
            parse_circular_track_slots([bad])
        _assert_dinucleotide_rejection(caught.value)
    with pytest.raises(ValueError):
        parse_circular_track_slots(
            [CircularTrackSlot(id="at", renderer="dinucleotide_skew", params={"nt": "G"})]
        )


def test_linear_slot_nt_uses_the_same_validator():
    slots = parse_linear_track_slots(["features:features", "at:dinucleotide_skew@nt=tu"])
    assert slots[1].params["nt"] == "TU"
    for bad in (
        ["features:features", "at:dinucleotide_skew@nt=G"],
        [LinearTrackSlot(id="features", renderer="features"),
         LinearTrackSlot(id="at", renderer="dinucleotide_content", params={"nt": "XY"})],
    ):
        with pytest.raises(ValueError) as caught:
            parse_linear_track_slots(bad)
        _assert_dinucleotide_rejection(caught.value)


def test_resolved_slot_without_nt_uses_default_and_never_silently_repairs():
    assert slot_dinucleotide({}, "GC") == "GC"
    assert slot_dinucleotide({"nt": ""}, "au") == "AU"
    assert slot_dinucleotide({"dinucleotide": "at"}, "GC") == "AT"
    with pytest.raises(ValidationError):
        slot_dinucleotide({"nt": "G"}, "GC")
