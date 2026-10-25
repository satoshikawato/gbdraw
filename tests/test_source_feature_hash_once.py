"""The source feature hash is computed once per source feature (L1-03).

Loading a source captures its features: each gets its source hash once,
stamped on the feature, and every crop, reverse complement or GFF3 type filter
copy carries the stamp. Readers return the stamp for the record it was made for
and compute only for another record ID. The source feature catalog builds its
entries only when a reader asks for them.
"""

from __future__ import annotations

import hashlib
import types
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

import gbdraw.features.ids as ids_module
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import build_request_diagram, plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
)
from gbdraw.core.record_metadata import _iter_source_features, _source_feature_index
from gbdraw.features.ids import (
    compute_feature_hash_from_location_parts,
    compute_source_feature_hash,
    source_feature_hash_and_parts,
    source_feature_location_parts,
)
from gbdraw.features.source import SourceFeatureIdentity
from gbdraw.io.regions import parse_region_spec

INPUTS = Path(__file__).parent / "test_inputs"


def _record(record_id: str = "hash-once") -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 600), id=record_id, name=record_id, description="hash once fixture")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(SimpleLocation(0, 2400, strand=1), type="source", qualifiers={"organism": ["test"]}),
        SeqFeature(SimpleLocation(100, 400, strand=1), type="CDS",
                   qualifiers={"locus_tag": ["t1"], "product": ["alpha"], "translation": ["M"]}),
        SeqFeature(SimpleLocation(600, 900, strand=-1), type="CDS",
                   qualifiers={"locus_tag": ["t2"], "product": ["beta"], "translation": ["M"]}),
        # Identical to the one above: same hash, a disambiguated identity.
        SeqFeature(SimpleLocation(600, 900, strand=-1), type="CDS",
                   qualifiers={"locus_tag": ["t3"], "product": ["beta copy"], "translation": ["M"]}),
        SeqFeature(CompoundLocation([SimpleLocation(1200, 1500, 1), SimpleLocation(1600, 1800, 1)]),
                   type="CDS", qualifiers={"locus_tag": ["t4"], "product": ["joined"], "translation": ["M"]}),
        SeqFeature(SimpleLocation(2000, 2100, strand=1), type="tRNA", qualifiers={"product": ["tRNA-Ala"]}),
    ]
    return record


def _genbank(tmp_path: Path, *records: SeqRecord) -> Path:
    path = tmp_path / "hash-once.gbk"
    SeqIO.write(list(records) or [_record()], path, "genbank")
    return path


class _Counts:
    def __init__(self) -> None:
        self.md5 = 0
        self.entries = 0


@pytest.fixture
def counts(monkeypatch: pytest.MonkeyPatch) -> _Counts:
    """Count md5 digests of feature hashes and source catalog entries built."""

    counted = _Counts()

    def md5(data: bytes) -> object:
        counted.md5 += 1
        return hashlib.md5(data)

    monkeypatch.setattr(ids_module, "hashlib", types.SimpleNamespace(md5=md5))
    original_init = SourceFeatureIdentity.__init__

    def init(self: SourceFeatureIdentity, *args: object, **kwargs: object) -> None:
        counted.entries += 1
        original_init(self, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(SourceFeatureIdentity, "__init__", init)
    return counted


def _fresh_hash(feature: SeqFeature, record: SeqRecord, record_id: str | None = None) -> str | None:
    """The source hash computed from its inputs, never read from a stamp."""

    parts = source_feature_location_parts(feature, record)
    if not parts:
        return None
    return compute_feature_hash_from_location_parts(
        str(feature.type), parts, record_id=record.id if record_id is None else record_id,
    )


def test_a_plain_render_hashes_each_source_feature_once_and_builds_no_catalog_entry(
    tmp_path: Path, counts: _Counts,
) -> None:
    path = _genbank(tmp_path)
    request = CircularDiagramRequest(
        records=(RecordInput(GenBankInputSource(path), record_key="one"),),
        options=CircularDiagramOptions(selected_features_set=("CDS", "tRNA")),
    )
    assert build_request_diagram(request).drawing.tostring()
    assert counts.md5 == len(_record().features)
    assert counts.entries == 0


def test_a_linear_render_of_two_records_builds_no_catalog_entry(tmp_path: Path, counts: _Counts) -> None:
    first, second = _record("first"), _record("second")
    request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), record_key="first"),
            RecordInput(InMemoryRecordSource(second), record_key="second"),
        ),
        options=LinearDiagramOptions(selected_features_set=("CDS", "tRNA")),
    )
    assert build_request_diagram(request).drawing.tostring()
    drawn = sum(feature.type in {"CDS", "tRNA"} for record in (first, second) for feature in record.features)
    # One source hash per source feature; the linear DOM index also hashes each
    # drawn feature in its drawn frame (compute_feature_object_hash).
    assert counts.md5 == len(first.features) + len(second.features) + drawn
    assert counts.entries == 0


def test_a_qualifier_selector_table_builds_each_catalog_entry_once(counts: _Counts) -> None:
    table = DataFrame([
        {"record": "#1", "feature_selector": "locus_tag=t2", "feature_visibility": "off",
         "label_visibility": "", "label_text": ""},
        {"record": "#1", "feature_selector": "product=joined", "feature_visibility": "",
         "label_visibility": "on", "label_text": "J"},
    ])
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(_record()), record_key="one"),),
        options=LinearDiagramOptions(selected_features_set=("CDS",), feature_override_table=table),
    )
    plan = plan_request(request)
    assert len(plan.request.options.feature_overrides) == 2
    assert counts.entries == len(_record().features)
    # A second read reuses the entries.
    assert len(plan.provenance[0].source_feature_catalog) == len(_record().features)
    assert counts.entries == len(_record().features)


@pytest.mark.parametrize(
    ("region", "reverse"), [(None, True), ("101-1900", False), ("101-1900", True)],
    ids=["reverse", "crop", "crop-reverse"],
)
def test_a_transformed_record_reads_the_carried_hash_without_hashing(
    tmp_path: Path, counts: _Counts, region: str | None, reverse: bool,
) -> None:
    path = _genbank(tmp_path)
    plan = plan_request(CircularDiagramRequest(
        records=(RecordInput(
            GenBankInputSource(path), record_key="one",
            region=parse_region_spec(region) if region else None,
            presentation=RecordPresentation(reverse_complement=reverse),
        ),),
        options=CircularDiagramOptions(selected_features_set=("CDS", "tRNA")),
    ))
    record = plan.records[0]
    catalog = {entry.source_feature_index: entry for entry in plan.provenance[0].source_feature_catalog}
    drawn = list(_iter_source_features(record.features))
    assert drawn
    before = counts.md5
    read = [(compute_source_feature_hash(feature, record), source_feature_hash_and_parts(feature, record))
            for feature in drawn]
    assert counts.md5 == before
    for feature, (feature_hash, (paired_hash, parts)) in zip(drawn, read, strict=True):
        entry = catalog[_source_feature_index(feature)]
        assert feature_hash == paired_hash == entry.stable_feature_id == _fresh_hash(feature, record)
        assert parts == entry.location_parts


def test_a_gff_type_filter_reads_the_carried_hash_without_hashing(counts: _Counts) -> None:
    plan = plan_request(CircularDiagramRequest(
        records=(RecordInput(
            GffFastaInputSource(INPUTS / "NC_013668.gff3", INPUTS / "NC_013668.fasta"),
            record_key="gff", cardinality=RecordCardinality.ALL,
        ),),
        options=CircularDiagramOptions(selected_features_set=("CDS", "tRNA")),
    ))
    record = plan.records[0]
    catalog = {entry.source_feature_index: entry for entry in plan.provenance[0].source_feature_catalog}
    drawn = list(_iter_source_features(record.features))
    assert drawn and {feature.type for feature in drawn} <= {"CDS", "tRNA"}
    before = counts.md5
    hashes = [compute_source_feature_hash(feature, record) for feature in drawn]
    assert counts.md5 == before
    for feature, feature_hash in zip(drawn, hashes, strict=True):
        assert feature_hash == catalog[_source_feature_index(feature)].stable_feature_id
        assert feature_hash == _fresh_hash(feature, record)


def test_another_record_id_computes_its_own_hash() -> None:
    record = plan_request(CircularDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(_record()), record_key="one"),),
        options=CircularDiagramOptions(),
    )).records[0]
    for feature in record.features:
        assert compute_source_feature_hash(feature, record, record_id="other") == _fresh_hash(
            feature, record, "other"
        )
        assert source_feature_hash_and_parts(feature, record, record_id="other")[0] == _fresh_hash(
            feature, record, "other"
        )
        assert compute_source_feature_hash(feature, record) == _fresh_hash(feature, record)


def test_a_planned_record_moved_and_planned_again_gets_a_fresh_hash() -> None:
    def plan(record: SeqRecord):
        return plan_request(CircularDiagramRequest(
            records=(RecordInput(InMemoryRecordSource(record), record_key="one"),),
            options=CircularDiagramOptions(),
        ))

    planned = plan(_record()).records[0]
    moved = planned.features[1]
    old_hash = compute_source_feature_hash(moved, planned)
    moved.location = SimpleLocation(150, 450, strand=1)
    replanned = plan(planned)
    record = replanned.records[0]
    entry = replanned.provenance[0].source_feature_catalog[1]
    assert entry.stable_feature_id == _fresh_hash(record.features[1], record) != old_hash
    assert compute_source_feature_hash(record.features[1], record) == entry.stable_feature_id


def test_catalog_entries_equal_the_source_features_they_capture() -> None:
    from gbdraw.features.ids import disambiguate_feature_ids
    from gbdraw.features.selector_values import normalize_qualifier_values
    from gbdraw.features.source import capture_source_features

    record = _record()
    catalog = capture_source_features(record)
    hashes = [_fresh_hash(feature, record) for feature in record.features]
    identities = disambiguate_feature_ids(
        (record.id, feature_hash, index) for index, feature_hash in enumerate(hashes)
    )
    assert catalog.entries == tuple(
        SourceFeatureIdentity(
            biological_feature_id=identity,
            source_feature_index=index,
            feature_type=feature.type,
            location_parts=source_feature_location_parts(feature, record),
            qualifiers=tuple(
                (key, tuple(normalize_qualifier_values(value))) for key, value in feature.qualifiers.items()
            ),
            stable_feature_id=feature_hash,
        )
        for index, (feature, feature_hash, identity) in enumerate(zip(record.features, hashes, identities))
    )
    assert catalog.entries[3].biological_feature_id == f"{hashes[3]}~3"
    assert catalog.source_hashes == dict(enumerate(hashes))
    assert catalog.entries is catalog.entries
