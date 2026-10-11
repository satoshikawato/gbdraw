"""Shared original-source feature identity resolution (design Q4, section 3.1)."""

from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import CircularDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GffFastaInputSource,
    InMemoryRecordSource,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import FeatureIdentity, resolve_feature_identities
from gbdraw.io.regions import parse_region_spec


def _record() -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 40), id="dup", annotations={"molecule_type": "DNA", "topology": "linear"}
    )
    record.features = [
        SeqFeature(SimpleLocation(20, 40, 1), type="CDS", qualifiers={"locus_tag": [tag]})
        for tag in ("first", "second")
    ] + [SeqFeature(SimpleLocation(90, 110, 1), type="CDS", qualifiers={"locus_tag": ["late"]})]
    return record


def _plan(*inputs: RecordInput, **options):
    return plan_request(
        CircularDiagramRequest(records=inputs, options=CircularDiagramOptions(**options))
    )


def _resolve(plan, *identities: FeatureIdentity):
    return resolve_feature_identities(
        records=plan.records,
        record_keys=tuple(item.record_key for item in plan.provenance),
        source_catalogs=tuple(item.source_feature_catalog for item in plan.provenance),
        identities=identities,
    )


def test_identical_features_bind_to_their_own_source_instance():
    plan = _plan(RecordInput(InMemoryRecordSource(_record()), record_key="one"))
    first, second, _late = plan.provenance[0].source_feature_catalog
    assert (first.biological_feature_id, second.biological_feature_id) == (
        f"{first.stable_feature_id}~0",
        f"{first.stable_feature_id}~1",
    )
    identities = [FeatureIdentity("one", entry.biological_feature_id) for entry in (first, second)]
    bindings = _resolve(plan, *identities)
    assert [bindings[item].status for item in identities] == ["present", "present"]
    assert [bindings[item].source_feature_index for item in identities] == [0, 1]
    assert [bindings[item].feature.qualifiers["locus_tag"] for item in identities] == [
        ["first"],
        ["second"],
    ]


def test_crop_excluded_and_unresolved_identities_stay_distinct():
    plan = _plan(
        RecordInput(
            InMemoryRecordSource(_record()), record_key="one", region=parse_region_spec("80-120")
        )
    )
    catalog = plan.provenance[0].source_feature_catalog
    cropped, kept, stale = (
        FeatureIdentity("one", catalog[0].biological_feature_id),
        FeatureIdentity("one", catalog[2].biological_feature_id),
        FeatureIdentity("one", "stale"),
    )
    bindings = _resolve(plan, cropped, kept, stale)
    assert [(bindings[item].status, bindings[item].source_feature_index) for item in (cropped, kept, stale)] == [
        ("crop_excluded", 0),
        ("present", 2),
        ("unresolved", None),
    ]
    assert bindings[cropped].feature is None and bindings[stale].feature is None


def test_feature_dropped_by_loading_is_absent(tmp_path):
    gff, fasta = tmp_path / "source.gff3", tmp_path / "source.fasta"
    fasta.write_text(">record\n" + "ATG" * 40 + "\n")
    gff.write_text(
        "##gff-version 3\n"
        "record\t.\tgene\t1\t18\t.\t+\t.\tID=hidden-gene\n"
        "record\t.\tCDS\t31\t60\t.\t+\t0\tID=shown-cds\n"
    )
    plan = _plan(
        RecordInput(GffFastaInputSource(gff, fasta), record_key="one"),
        selected_features_set=("CDS",),
    )
    gene = FeatureIdentity("one", plan.provenance[0].source_feature_catalog[0].biological_feature_id)
    binding = _resolve(plan, gene)[gene]
    assert (binding.status, binding.source_feature_index, binding.feature) == ("absent", 0, None)


def test_record_key_outside_the_request_is_a_feature_identity_error():
    plan = _plan(RecordInput(InMemoryRecordSource(_record()), record_key="one"))
    with pytest.raises(ValidationError, match="Unknown feature identity record key") as raised:
        _resolve(plan, FeatureIdentity("other", "f1234"))
    assert raised.value.diagnostic == {"code": "FEATURE_IDENTITY"}


@pytest.mark.parametrize(
    "record_key,feature_id", [("", "f1"), ("one\0f1", "f1"), ("one", "f1\0"), (None, "f1")]
)
def test_identity_requires_non_empty_values_without_nul(record_key, feature_id):
    with pytest.raises(ValidationError) as raised:
        FeatureIdentity(record_key, feature_id)
    assert raised.value.diagnostic == {"code": "FEATURE_IDENTITY"}
