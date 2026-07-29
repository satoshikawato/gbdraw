from __future__ import annotations

import inspect

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.diagrams.linear.assemble as linear_assemble_module
import gbdraw.diagrams.linear.precalc as linear_precalc_module
import gbdraw.render.groups.circular.labels as circular_labels_module
import gbdraw.render.groups.circular.seq_record as circular_record_module
import gbdraw.render.groups.linear.seq_record as linear_record_module
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.diagrams.circular.builders import (
    add_labels_group_on_canvas,
    add_record_group_on_canvas,
)
from gbdraw.diagrams.linear.builders import add_record_group


def _record(record_id: str = "record_1") -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 120), id=record_id)
    record.features = [
        SeqFeature(FeatureLocation(20, 180, strand=1), type="repeat_region"),
        SeqFeature(FeatureLocation(80, 240, strand=1), type="CDS"),
    ]
    return record


@pytest.mark.parametrize(
    "target",
    [
        circular_record_module.SeqRecordGroup,
        circular_labels_module.LabelsGroup,
        linear_record_module.SeqRecordGroup,
        add_record_group_on_canvas,
        add_labels_group_on_canvas,
        add_record_group,
    ],
)
def test_feature_render_boundaries_require_prepared_layers(target) -> None:
    parameters = inspect.signature(target).parameters

    assert parameters["feature_layers"].default is inspect.Parameter.empty
    assert "precomputed_feature_dict" not in parameters


@pytest.mark.parametrize(
    "module",
    [
        circular_record_module,
        circular_labels_module,
        linear_record_module,
    ],
)
def test_render_groups_do_not_import_legacy_feature_conversion(module) -> None:
    assert "create_feature_dict" not in vars(module)


@pytest.mark.circular
def test_circular_groups_receive_the_same_prepared_feature_layers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    built_layers = []
    record_group_layers = []
    label_group_layers = []
    real_create_feature_layers = circular_assemble_module.create_feature_layers

    def capture_create_feature_layers(*args, **kwargs):
        result = real_create_feature_layers(*args, **kwargs)
        built_layers.append(result)
        return result

    def capture_record_group(canvas, *args, **kwargs):
        record_group_layers.append(kwargs["feature_layers"])
        return canvas

    def capture_label_group(canvas, *args, **kwargs):
        label_group_layers.append(kwargs["feature_layers"])
        return canvas

    monkeypatch.setattr(
        circular_assemble_module,
        "create_feature_layers",
        capture_create_feature_layers,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_record_group_on_canvas",
        capture_record_group,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_labels_group_on_canvas",
        capture_label_group,
    )

    assemble_circular_diagram_from_record(
        _record(),
        selected_features_set=["repeat_region", "CDS"],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "outer",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )

    assert len(built_layers) == 1
    assert len(built_layers[0].underlay_features) == 1
    assert record_group_layers == [built_layers[0]]
    assert label_group_layers == [built_layers[0], built_layers[0]]


@pytest.mark.linear
def test_linear_groups_receive_prepared_feature_layers_per_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    built_layers = []
    record_group_layers = []
    real_create_feature_layers = linear_precalc_module.create_feature_layers

    def capture_create_feature_layers(*args, **kwargs):
        result = real_create_feature_layers(*args, **kwargs)
        built_layers.append(result)
        return result

    def capture_record_group(canvas, *args, **kwargs):
        record_group_layers.append(kwargs["feature_layers"])
        return canvas

    monkeypatch.setattr(
        linear_precalc_module,
        "create_feature_layers",
        capture_create_feature_layers,
    )
    monkeypatch.setattr(
        linear_assemble_module,
        "add_record_group",
        capture_record_group,
    )

    records = [_record("record_1"), _record("record_2")]
    assemble_linear_diagram_from_records(
        records,
        selected_features_set=["repeat_region", "CDS"],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )

    assert len(built_layers) == len(records)
    assert all(len(result.underlay_features) == 1 for result in built_layers)
    assert record_group_layers == built_layers
