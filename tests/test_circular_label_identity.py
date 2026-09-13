"""Rendered labels keep source identity across multipart and repeated records."""

from copy import deepcopy
from dataclasses import replace
from xml.etree import ElementTree as ET

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, RecordDisplayOptions, RecordInput
from gbdraw.features.ids import compute_feature_hash


@pytest.mark.parametrize("placement", ["horizontal", "radial"])
@pytest.mark.parametrize("rendering", ["auto", "external_only", "embedded_only"])
@pytest.mark.parametrize("repeated_record", [False, True])
def test_labels_resolve_all_parts_and_disambiguate_equal_source_hashes(
    placement, rendering, repeated_record
):
    record = SeqRecord(Seq("ATG" * 1000), id="label-source", annotations={
        "molecule_type": "DNA", "topology": "circular",
    })
    multipart = SeqFeature(CompoundLocation([
        SimpleLocation(100, 900, strand=-1),
        SimpleLocation(2000, 2300, strand=-1),
        SimpleLocation(1900, 1920, strand=-1),
    ]), type="CDS", qualifiers={"product": ["joined"], "trans_splicing": [""]})
    duplicate = deepcopy(multipart)
    duplicate.qualifiers["product"] = ["duplicate"]
    record.features = [multipart, duplicate, SeqFeature(
        SimpleLocation(2400, 2900, strand=1), type="CDS",
        qualifiers={"product": ["single"]},
    )]
    count = 2 if repeated_record else 1
    request = CircularDiagramRequest(
        records=tuple(RecordInput(InMemoryRecordSource(record),
            display=RecordDisplayOptions(start_coordinate=1)) for _ in range(count)),
        grouping="grid" if repeated_record else "single",
        options=CircularDiagramOptions(config=apply_config_overrides(None, {
            "labels.circular.scope": "outer",
            "labels.circular.placement": placement,
            "labels.rendering": rendering,
            "labels.filtering.blacklist_keywords": [],
            "canvas.show_gc": False, "canvas.show_skew": False,
        })),
    )
    drawing = plan_request(request).build()
    root = ET.fromstring(drawing.tostring())
    paths = [node for node in root.iter() if node.get("data-gbdraw-feature-part") == "block"]
    labels = [node for node in root.iter() if node.get("data-label-feature-id")]
    assert len(labels) == 3 * count
    assert len({node.get("data-label-feature-id") for node in labels}) == 3 * count
    for label in labels:
        identity = label.get("data-label-feature-id")
        members = [node for node in paths if (
            node.get("data-gbdraw-rendered-feature-id") or node.get("data-gbdraw-feature-id")
        ) == identity]
        text = "".join(label.itertext())
        assert len(members) == (1 if text == "single" else 3)
        source = next(feature for feature in record.features if feature.qualifiers["product"] == [text])
        assert {node.get("data-gbdraw-feature-id") for node in members} == {
            compute_feature_hash(source, record_id=record.id)
        }


def test_repeated_render_uses_current_label_table_after_filtering_was_prepared():
    record = SeqRecord(Seq("ATG" * 1000), id="repeated", annotations={
        "molecule_type": "DNA", "topology": "circular",
    })
    record.features = [SeqFeature(SimpleLocation(10, 500, strand=1), type="CDS",
        qualifiers={"product": ["original"], "locus_tag": ["target"]})]
    cfg = apply_config_overrides(None, {"labels.circular.scope": "outer",
        "canvas.show_gc": False, "canvas.show_skew": False})
    request = CircularDiagramRequest(records=(RecordInput(InMemoryRecordSource(record)),),
        options=CircularDiagramOptions(config=cfg))
    # The first render compiles filtering on the same typed config reused later.
    plan_request(request).build()
    for replacement in ["first edit", "second edit"]:
        table = pd.DataFrame([[record.id, "CDS", "locus_tag", "^target$", replacement]],
            columns=["record_id", "feature_type", "qualifier", "value", "label_text"])
        edited = replace(request, options=replace(request.options, label_override_table=table))
        root = ET.fromstring(plan_request(edited).build().tostring())
        assert ["".join(node.itertext()) for node in root.iter()
            if node.get("data-label-feature-id")] == [replacement]
