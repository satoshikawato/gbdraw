"""Linear labels and leaders share the exact rendered feature identity."""

from copy import deepcopy
from xml.etree import ElementTree as ET

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordInput


def _tag_name(element: ET.Element) -> str:
    return str(element.tag).rsplit("}", 1)[-1]


def _record(record_id: str) -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 1000),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )
    first = SeqFeature(
        SimpleLocation(100, 300, strand=1),
        type="CDS",
        qualifiers={"product": ["same label"]},
    )
    duplicate = deepcopy(first)
    record.features = [first, duplicate]
    return record


def test_external_labels_bind_text_and_leaders_to_duplicate_multirecord_instances() -> None:
    record = _record("repeated")
    config = apply_config_overrides(None, {
        "labels.linear.scope": "all",
        "labels.rendering": "external_only",
        "labels.filtering.blacklist_keywords": [],
    })
    request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(record)),
            RecordInput(InMemoryRecordSource(deepcopy(record))),
        ),
        options=LinearDiagramOptions(config=config),
        layout=LinearMultiRecordOptions(),
    )

    root = ET.fromstring(plan_request(request).build().drawing.tostring())
    labels = [
        element for element in root.iter()
        if _tag_name(element) == "text"
        and element.get("data-gbdraw-label-binding-schema") == "1"
    ]
    leaders = [
        element for element in root.iter()
        if _tag_name(element) == "line" and element.get("data-label-feature-id")
    ]
    feature_parts = [
        element for element in root.iter()
        if element.get("data-gbdraw-feature-part") == "block"
    ]

    assert len(labels) == 4
    identities = [label.get("data-label-feature-id") for label in labels]
    assert all(identities)
    assert len(set(identities)) == 4
    assert ["".join(label.itertext()) for label in labels] == ["same label"] * 4
    for identity in identities:
        assert sum(
            leader.get("data-label-feature-id") == identity
            for leader in leaders
        ) == 1
        assert any(
            (
                part.get("data-gbdraw-rendered-feature-id")
                or part.get("data-gbdraw-feature-id")
            ) == identity
            for part in feature_parts
        )


def test_embedded_linear_label_declares_complete_zero_leader_binding() -> None:
    record = SeqRecord(
        Seq("ATG" * 1000),
        id="embedded",
        annotations={"molecule_type": "DNA"},
    )
    record.features = [SeqFeature(
        SimpleLocation(100, 2500, strand=1),
        type="CDS",
        qualifiers={"product": ["short"]},
    )]
    config = apply_config_overrides(None, {
        "labels.linear.scope": "all",
        "labels.rendering": "auto",
        "labels.filtering.blacklist_keywords": [],
    })
    root = ET.fromstring(plan_request(LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(config=config),
    )).build().drawing.tostring())

    labels = [
        element for element in root.iter()
        if _tag_name(element) == "text"
        and element.get("data-gbdraw-label-binding-schema") == "1"
    ]
    assert len(labels) == 1
    identity = labels[0].get("data-label-feature-id")
    assert identity
    assert not any(
        _tag_name(element) == "line"
        and element.get("data-label-feature-id") == identity
        for element in root.iter()
    )
