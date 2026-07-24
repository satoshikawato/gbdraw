from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    extract_cds_proteins,
)
from gbdraw.features.ids import compute_feature_hash, make_linear_rendered_feature_id
from gbdraw.io.record_select import reverse_records
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import enrich_svg
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload


INPUTS = Path(__file__).parent / "test_inputs"
SOURCE_PROTEIN_ID = "CAG34720.1"
BIOLOGICAL_FEATURE_ID = "feb87ab70"
PROCESSED_FEATURE_ID = "f20e4885e"


def _bgc_record_five_render_input():
    """Reproduce the saved BGC record-5 source orientation and render reversal."""

    downloaded = SeqIO.read(INPUTS / "BGC0000713.gbk", "genbank")
    saved_source = downloaded.reverse_complement(
        id=True,
        name=True,
        description=True,
        features=True,
        annotations=True,
        letter_annotations=True,
        dbxrefs=True,
    )
    saved_source.annotations.pop("gbdraw_coord_base", None)
    saved_source.annotations.pop("gbdraw_coord_step", None)
    return reverse_records([saved_source], True)[0]


def _target_feature(record):
    return next(
        feature
        for feature in record.features
        if feature.qualifiers.get("protein_id") == [SOURCE_PROTEIN_ID]
    )


def _orthogroups_for_target(record) -> OrthogroupResult:
    feature = _target_feature(record)
    member = OrthogroupMember(
        orthogroup_id="og_1",
        protein_id="protein-1",
        source_protein_id=SOURCE_PROTEIN_ID,
        record_index=4,
        feature_index=record.features.index(feature),
        record_id=record.id,
        label=SOURCE_PROTEIN_ID,
        start=int(feature.location.start),
        end=int(feature.location.end),
        strand=feature.location.strand,
        feature_svg_id=BIOLOGICAL_FEATURE_ID,
        gene="rbmC",
        product="putative dehydrogenase",
        note=None,
        representative=True,
    )
    return OrthogroupResult(
        orthogroups={"og_1": [member]},
        member_by_protein_id={member.protein_id: member},
    )


def _metadata(svg: str) -> dict[str, object]:
    root = ET.fromstring(svg)
    metadata = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "metadata"
        and element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    return json.loads(metadata.text or "{}")


def test_reverse_complement_metadata_separates_biological_and_rendered_ids() -> None:
    record = _bgc_record_five_render_input()
    target = _target_feature(record)
    assert compute_feature_hash(target, record_id=record.id) == PROCESSED_FEATURE_ID

    payload = extract_features_from_records_payload(
        [record],
        selected_features=["CDS"],
        linear_rendered_feature_ids=True,
        include_biological_features=True,
    )
    feature = next(
        item for item in payload["features"] if item["protein_id"] == SOURCE_PROTEIN_ID
    )
    biological = next(
        item
        for item in payload["biological_features"]
        if item["protein_id"] == SOURCE_PROTEIN_ID
    )

    assert feature["svg_id"] == BIOLOGICAL_FEATURE_ID
    assert feature["stable_feature_id"] == BIOLOGICAL_FEATURE_ID
    assert feature["rendered_feature_svg_id"] == PROCESSED_FEATURE_ID
    assert feature["selector"]["hash"] == BIOLOGICAL_FEATURE_ID
    assert feature["selector"]["record_location"] == "BGC0000713:7134..8157:+"
    assert feature["location_parts"] == [
        {"start": 7134, "end": 8157, "strand": "+", "display": "7135..8157"}
    ]
    assert feature["strand"] == "+"
    assert biological["svg_id"] == BIOLOGICAL_FEATURE_ID
    assert "rendered_feature_svg_id" not in biological


def test_reverse_complement_protein_extraction_uses_biological_hash_parts() -> None:
    record = _bgc_record_five_render_input()
    proteins = extract_cds_proteins([record], record_index_offset=4)
    protein = next(
        item
        for item in proteins.proteins_by_record[0]
        if item.source_protein_id == SOURCE_PROTEIN_ID
    )

    assert protein.feature_svg_id == BIOLOGICAL_FEATURE_ID
    assert protein.feature_hash_start == 7134
    assert protein.feature_hash_end == 8157
    assert protein.feature_hash_strand == 1
    assert protein.feature_hash_parts == ((7134, 8157, 1),)
    assert (protein.start, protein.end, protein.strand) == (23735, 24758, -1)


def test_interactive_svg_maps_biological_id_to_actual_reversed_dom_path() -> None:
    record = _bgc_record_five_render_input()
    records = [SeqRecord(Seq("N"), id=f"placeholder-{index}") for index in range(4)] + [
        record
    ]
    rendered_id = make_linear_rendered_feature_id(
        record_index=4,
        stable_feature_id=PROCESSED_FEATURE_ID,
        record_count=5,
    )
    assert rendered_id == "f20e4885e_record_5"

    context = build_interactive_svg_context(
        records,
        selected_features_set=["CDS"],
        orthogroups=_orthogroups_for_target(record),
        linear_rendered_feature_ids=True,
    )
    context_feature = next(
        item
        for item in context.features
        if item["source_protein_id"] == SOURCE_PROTEIN_ID
    )
    assert context_feature["stable_feature_id"] == BIOLOGICAL_FEATURE_ID
    assert context_feature["rendered_feature_svg_id"] == rendered_id
    assert context.orthogroups[0]["members"][0]["featureSvgId"] == BIOLOGICAL_FEATURE_ID

    source = f"""<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px"
      viewBox="0 0 100 80"><path id="{rendered_id}"
      data-gbdraw-feature-id="{rendered_id}"
      data-gbdraw-stable-feature-id="{PROCESSED_FEATURE_ID}"
      data-gbdraw-record-index="4" data-gbdraw-feature-part="block"
      fill="#54bcf8" d="M 1 1 L 2 2" /></svg>"""
    enriched = enrich_svg(source, context)
    payload = _metadata(enriched)

    feature = next(
        item
        for item in payload["features"]
        if item["source_protein_id"] == SOURCE_PROTEIN_ID
    )
    biological = next(
        item
        for item in payload["biological_features"]
        if item.get("source_protein_id") == SOURCE_PROTEIN_ID
    )
    member = payload["orthogroups"][0]["members"][0]
    assert feature["svg_id"] == rendered_id
    assert feature["rendered_feature_svg_id"] == rendered_id
    assert feature["stable_svg_id"] == BIOLOGICAL_FEATURE_ID
    assert biological["svg_id"] == BIOLOGICAL_FEATURE_ID
    assert biological["rendered_svg_id"] == rendered_id
    assert member["feature_svg_id"] == BIOLOGICAL_FEATURE_ID
    assert member["stable_feature_svg_id"] == BIOLOGICAL_FEATURE_ID
    assert member["rendered_feature_svg_id"] == rendered_id

    root = ET.fromstring(enriched)
    path = next(element for element in root.iter() if element.get("id") == rendered_id)
    assert path.get("data-gbdraw-feature-id") == rendered_id
    assert path.get("data-gbdraw-stable-feature-id") == BIOLOGICAL_FEATURE_ID
