from __future__ import annotations

from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.config.toml import load_config_toml
from gbdraw.diagrams.linear.builders import add_record_definition_group
from gbdraw.render.groups.linear import DefinitionGroup


def _record(label: str = "Record A", record_id: str = "record_a") -> SeqRecord:
    record = SeqRecord(Seq("A" * 100), id=record_id)
    record.annotations["gbdraw_record_label"] = label
    return record


def _canvas_config(*, keep_definition_left_aligned: bool) -> SimpleNamespace:
    return SimpleNamespace(
        canvas_padding=50,
        horizontal_offset=120,
        keep_definition_left_aligned=keep_definition_left_aligned,
        length_param="short",
    )


def _definition_translate_x(canvas: Drawing) -> float:
    definition_group = next(
        element for element in canvas.elements if element.attribs.get("id") == "record_a_definition"
    )
    transform = str(definition_group.attribs["transform"])
    return float(transform.split("translate(", 1)[1].split(",", 1)[0])


def _definition_text_anchors(canvas: Drawing) -> set[str]:
    definition_group = next(
        element for element in canvas.elements if element.attribs.get("id") == "record_a_definition"
    )
    return {
        str(element.attribs.get("text-anchor"))
        for element in definition_group.elements
        if getattr(element, "elementname", "") == "text"
    }


def _definition_width(record: SeqRecord, config_dict: dict, canvas_config: SimpleNamespace) -> float:
    return DefinitionGroup(record, config_dict, canvas_config).definition_bounding_box_width


@pytest.mark.linear
def test_linear_definition_group_follows_record_offset_by_default() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        _record(),
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        config_dict=config_dict,
        max_def_width=0,
    )
    add_record_definition_group(
        canvas_b,
        _record(),
        record_offset_y=10,
        record_offset_x=30,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        config_dict=config_dict,
        max_def_width=0,
    )

    assert _definition_translate_x(canvas_b) == pytest.approx(_definition_translate_x(canvas_a) + 30)


@pytest.mark.linear
def test_linear_definition_group_can_stay_in_left_column() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = _canvas_config(keep_definition_left_aligned=True)
    short_record = _record("Short")
    long_record = _record("A much longer definition label", "record_a")
    max_def_width = max(
        _definition_width(short_record, config_dict, canvas_config),
        _definition_width(long_record, config_dict, canvas_config),
    )
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        short_record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=max_def_width,
    )
    add_record_definition_group(
        canvas_b,
        long_record,
        record_offset_y=10,
        record_offset_x=45,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=max_def_width,
    )

    assert _definition_translate_x(canvas_a) == pytest.approx(_definition_translate_x(canvas_b))
    assert _definition_text_anchors(canvas_a) == {"start"}
    assert _definition_text_anchors(canvas_b) == {"start"}


@pytest.mark.linear
def test_locked_linear_definition_column_can_shift_left_of_negative_record_offsets() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        _record(),
        record_offset_y=10,
        record_offset_x=-30,
        canvas_config=_canvas_config(keep_definition_left_aligned=True),
        config_dict=config_dict,
        max_def_width=0,
        definition_column_shift_x=0,
    )
    add_record_definition_group(
        canvas_b,
        _record(),
        record_offset_y=10,
        record_offset_x=-30,
        canvas_config=_canvas_config(keep_definition_left_aligned=True),
        config_dict=config_dict,
        max_def_width=0,
        definition_column_shift_x=30,
    )

    assert _definition_translate_x(canvas_b) == pytest.approx(_definition_translate_x(canvas_a) - 30)
