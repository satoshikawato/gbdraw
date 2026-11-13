"""A Linear feature label is centred on the block its feature draws (OV-3452)."""
from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import gbdraw.diagrams.linear.precalc as linear_precalc
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordDisplayOptions, RecordInput
from gbdraw.features.objects import linear_drawn_parts
from gbdraw.io.regions import parse_region_spec
from gbdraw.layout.linear_coords import normalize_position_to_linear_track

_LENGTH = 10000


def _record() -> SeqRecord:
    record = SeqRecord(Seq("AAAACCGTGG" * (_LENGTH // 10)), id="rec", annotations={"molecule_type": "DNA"})
    record.features = [
        # At the record's first base: Biopython's start is 0, the drawn block starts at 1.
        SeqFeature(SimpleLocation(0, 1000, strand=1), type="CDS", qualifiers={"product": ["first base"]}),
        # A join whose first part starts at the first base, as long as its second part.
        SeqFeature(
            CompoundLocation([SimpleLocation(0, 1000, strand=-1), SimpleLocation(4000, 5000, strand=-1)]),
            type="CDS",
            qualifiers={"product": ["joined"]},
        ),
        SeqFeature(SimpleLocation(2000, 3000, strand=1), type="CDS", qualifiers={"product": ["inner"]}),
    ]
    return record


def _labels(view: str, monkeypatch: pytest.MonkeyPatch) -> list[tuple[dict, object, float, float]]:
    """Each label with its feature and the record length and width its labels are placed in."""

    prepare = linear_precalc.linear_record_labels
    captured: list[tuple[dict, list, float, float]] = []

    def capture(feature_dict, genome_length, alignment_width, *args, **kwargs):
        result = prepare(feature_dict, genome_length, alignment_width, *args, **kwargs)
        captured.append((feature_dict, result.labels, float(genome_length), float(alignment_width)))
        return result

    monkeypatch.setattr(linear_precalc, "linear_record_labels", capture)
    overrides: dict[str, object] = {
        "labels.linear.scope": "all", "labels.rendering": "external_only",
        "canvas.show_gc": False, "canvas.show_skew": False,
    }
    if view == "wrapped":
        overrides.update({
            "canvas.linear.wrap_length": 4000, "canvas.linear.track_layout": "above",
            "objects.scale.style": "ruler",
        })
    record_input = RecordInput(InMemoryRecordSource(_record()))
    if view == "cropped":
        record_input = RecordInput(InMemoryRecordSource(_record()), region=parse_region_spec("rec:1-6000"))
    elif view == "rotated":
        record_input = RecordInput(InMemoryRecordSource(_record()), display=RecordDisplayOptions(is_circular=True, start_coordinate=4001))
    request = LinearDiagramRequest(
        records=(record_input,), options=LinearDiagramOptions(config=apply_config_overrides(None, overrides)),
    )
    plan_request(request).build()
    feature_dict, labels, genome_length, alignment_width = captured[-1]
    by_text = {feature.label_text: feature for feature in feature_dict.values()}
    return [(label, by_text[label["label_text"]], genome_length, alignment_width) for label in labels]


@pytest.mark.linear
@pytest.mark.parametrize("view", ["whole", "cropped", "wrapped", "rotated"])
def test_linear_label_is_centred_on_its_longest_drawn_block(view: str, monkeypatch: pytest.MonkeyPatch) -> None:
    labels = _labels(view, monkeypatch)
    assert {label["label_text"] for label, *_ in labels} == {"first base", "joined", "inner"}
    for label, feature, genome_length, alignment_width in labels:
        # A whole or cropped record draws its features' locations; a wrapped or rotated one their projected
        # display parts.
        assert (feature.display_parts is None) == (view in {"whole", "cropped"})
        def x(position: int) -> float:
            return normalize_position_to_linear_track(position, genome_length, alignment_width)

        # Each drawn block as (start x, end x, wrap row); the label's is one of the longest.
        blocks = [
            (x(int(part.location.start)), x(int(part.location.end)), part.row)
            for part in linear_drawn_parts(feature) if part.location.kind == "block"
        ]
        longest = max(abs(end - start) for start, end, _row in blocks)
        matches = [
            block for block in blocks
            if abs(block[1] - block[0]) == pytest.approx(longest)
            and (block[0], block[1]) == pytest.approx((label["feature_start_x"], label["feature_end_x"]))
        ]
        assert matches, (label["label_text"], label["feature_start_x"], label["feature_end_x"], blocks)
        start_x, end_x, row = matches[0]
        assert label["feature_anchor_x"] == pytest.approx((start_x + end_x) / 2.0)
        assert label["wrap_row"] == row if view == "wrapped" else "wrap_row" not in label


@pytest.mark.linear
def test_linear_label_of_a_feature_at_the_first_base_is_centred_on_its_drawn_block(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Base 1 to 1000 is drawn from position 1 to 1000; its label sits at 500.5, not at Biopython's 500."""

    (label, _feature, genome_length, alignment_width), = [
        entry for entry in _labels("whole", monkeypatch) if entry[0]["label_text"] == "first base"
    ]
    assert genome_length == _LENGTH
    assert label["feature_anchor_x"] == pytest.approx(alignment_width * 500.5 / _LENGTH)
