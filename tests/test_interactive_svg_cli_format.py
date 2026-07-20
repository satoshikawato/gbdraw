from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.api import render as api_render
from gbdraw.api import (
    InteractiveSvgContext as PublicInteractiveSvgContext,
    build_interactive_svg_context,
    enrich_svg as public_enrich_svg,
)
from gbdraw import circular as circular_cli
from gbdraw import linear as linear_cli
from gbdraw.cli_utils import common as cli_common
from gbdraw.exceptions import ExportError, GbdrawError, ValidationError
from gbdraw.render import export as export_module
from gbdraw.render.formats import (
    INTERACTIVE_SVG_FORMAT,
    SVG_FORMAT,
    classify_formats,
    parse_format_string,
    resolve_format_output_path,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg
from gbdraw.features.ids import compute_feature_hash, compute_feature_hash_from_parts
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload
from gbdraw.web_support.orthogroup_metadata import enrich_features_with_orthogroups


def _drawing(path: Path) -> Drawing:
    drawing = Drawing(filename=str(path), size=("120px", "80px"))
    drawing.add(drawing.path(d="M 10 10 L 20 20", id="fabc12345", fill="#54bcf8"))
    return drawing


def _metadata_payload(svg_source: str) -> dict[str, object]:
    root = ET.fromstring(svg_source)
    metadata = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "metadata"
        and element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    return json.loads(metadata.text or "{}")


def _metadata_text(svg_source: str) -> str:
    root = ET.fromstring(svg_source)
    metadata = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "metadata"
        and element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    return metadata.text or "{}"


def _script_payload(svg_source: str) -> str:
    root = ET.fromstring(svg_source)
    script = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "script"
        and element.get("id") == "gbdraw-interactive-feature-script"
    )
    return script.text or ""


def test_interactive_types_and_builder_are_public() -> None:
    feature = SeqFeature(
        SimpleLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"gene": ["geneA"], "translation": ["M"]},
    )
    record = SeqRecord(Seq("ATG" * 4), id="rec1", features=[feature])

    context = build_interactive_svg_context(
        [record],
        selected_features_set=["CDS"],
    )

    assert PublicInteractiveSvgContext is InteractiveSvgContext
    assert public_enrich_svg is enrich_svg
    assert context.features[0]["record_id"] == "rec1"
    assert context.features[0]["type"] == "CDS"


def test_linear_interactive_context_keeps_stable_feature_ids_until_svg_enrichment() -> (
    None
):
    records = []
    for record_index in range(2):
        feature = SeqFeature(
            SimpleLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={"gene": [f"gene{record_index}"], "translation": ["M"]},
        )
        records.append(
            SeqRecord(
                Seq("ATG" * 4),
                id=f"rec{record_index}",
                features=[feature],
            )
        )

    context = build_interactive_svg_context(
        records,
        selected_features_set=["CDS"],
        linear_rendered_feature_ids=True,
    )

    assert len(context.features) == 2
    assert all(
        feature["svg_id"] == feature["stable_feature_id"]
        for feature in context.features
    )
    assert [
        feature["stable_feature_id"] for feature in context.biological_features
    ] == [feature["stable_feature_id"] for feature in context.features]


def test_orthogroup_enrichment_requires_record_index_for_ambiguous_stable_ids() -> None:
    groups = [
        {
            "id": f"og_{record_index}",
            "member_count": 1,
            "record_coverage_count": 1,
            "members": [
                {
                    "featureSvgId": "fshared",
                    "stableFeatureSvgId": "fshared",
                    "recordIndex": record_index,
                    "featureIndex": 0,
                    "proteinId": f"protein_{record_index}",
                }
            ],
        }
        for record_index in range(2)
    ]

    enriched = enrich_features_with_orthogroups(
        [
            {"svg_id": "fshared", "stable_feature_id": "fshared", "record_idx": 1},
            {"svg_id": "fshared", "stable_feature_id": "fshared", "record_index": 0},
            {"svg_id": "fshared", "stable_feature_id": "fshared"},
        ],
        groups,
    )

    assert enriched[0]["orthogroup_id"] == "og_1"
    assert enriched[1]["orthogroup_id"] == "og_0"
    assert "orthogroup_id" not in enriched[2]


def test_parse_formats_accepts_interactive_svg_alias_and_dedupes() -> None:
    assert export_module.parse_formats("svg,interactive-svg,interactive_svg") == [
        SVG_FORMAT,
        INTERACTIVE_SVG_FORMAT,
    ]
    assert parse_format_string("interactive_svg") == [INTERACTIVE_SVG_FORMAT]


def test_format_classification_and_filename_resolution() -> None:
    classification = classify_formats(["interactive-svg", "png", "svg"])

    assert classification.interactive is True
    assert classification.cairosvg == ("png",)
    assert resolve_format_output_path("out", "interactive_svg") == "out.interactive.svg"


def test_handle_output_formats_keeps_interactive_when_cairosvg_missing(monkeypatch) -> None:
    monkeypatch.setattr(cli_common, "has_cairosvg", lambda: False)

    assert cli_common.handle_output_formats(["interactive-svg", "png", "pdf"]) == [
        "interactive_svg"
    ]


def test_enrich_svg_marks_features_matches_and_embeds_assets() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px" viewBox="0 0 100 80">
      <defs />
      <path id="fabc12345" data-gbdraw-feature-id="fabc12345" fill="#54bcf8" d="M 1 1 L 2 2" />
      <path data-gbdraw-pairwise-match-id="m1" data-match-kind="pairwise" data-identity="98.7" d="M 1 3 L 2 4" />
    </svg>"""

    enriched = enrich_svg(
        svg,
        InteractiveSvgContext(
            features=[
                {
                    "svg_id": "fabc12345",
                    "label": "geneA",
                    "type": "CDS",
                    "record_id": "rec1",
                }
            ]
        ),
    )

    assert 'data-gbdraw-interactive-svg="true"' in enriched
    assert "gbdraw-interactive-feature-style" in enriched
    assert "gbdraw-interactive-feature-script" in enriched
    assert "gbdraw-interactive-feature-glow" in enriched
    assert 'data-gbdraw-interactive-feature="true"' in enriched
    assert 'data-gbdraw-interactive-match="true"' in enriched
    payload = _metadata_payload(enriched)
    assert payload["features"][0]["svg_id"] == "fabc12345"
    assert payload["matches"][0]["id"] == "m1"


def test_enrich_svg_embeds_valid_regexp_escape_runtime() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fabc12345" data-gbdraw-feature-id="fabc12345" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""

    script = _script_payload(enrich_svg(svg))

    assert "replace(/[.*+?^${}()|[\\]\\\\]/g, '\\\\$&');" in script
    assert "replace(/[.*+?^${}()|[\\]\\]/g, '\\\\$&');" not in script


def test_enrich_svg_embeds_match_sequence_sources_once() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path data-gbdraw-match-id="comparison1_match1"
        data-gbdraw-pairwise-match-id="comparison1_match1" data-match-kind="pairwise"
        data-query-record-index="0" data-subject-record-index="1"
        data-query-record-id="query" data-subject-record-id="subject"
        data-qstart="2" data-qend="4" data-sstart="5" data-send="3"
        fill="#ff7272" d="M 1 1 L 2 2" />
    </svg>"""
    sources = [
        {"key": "linear:record:0", "recordId": "query", "sequence": "AACCGG", "origin": "linear-record", "recordIndex": 0},
        {"key": "linear:record:1", "recordId": "subject", "sequence": "TTGGCC", "origin": "linear-record", "recordIndex": 1},
        {"key": "linear:record:2", "recordId": "unused", "sequence": "NNNN", "origin": "linear-record", "recordIndex": 2},
    ]

    enriched = enrich_svg(svg, InteractiveSvgContext(sequence_sources=sources))
    payload = _metadata_payload(enriched)

    assert payload["matches"][0]["id"] == "comparison1_match1"
    assert [source["key"] for source in payload["sequence_sources"]] == [
        "linear:record:0",
        "linear:record:1",
    ]
    assert enriched.count('"sequence":"AACCGG"') == 1
    script = _script_payload(enriched)
    assert "buildEmbeddedMatchBundle" in script
    assert "coords=" in script


def test_enrich_svg_generates_fallback_feature_payload() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fabc12345__line1" data-gbdraw-feature-id="fabc12345"
        data-gbdraw-feature-part="connector" fill="none" d="M 1 1 A 2 2 0 0 1 2 2" />
      <path id="fabc12345__part1" data-gbdraw-feature-id="fabc12345" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""

    payload = _metadata_payload(enrich_svg(svg))

    assert payload["features"][0]["svg_id"] == "fabc12345"
    assert payload["features"][0]["fill_color"] == "#54bcf8"


def test_enrich_svg_v2_embeds_sequences_without_precomputed_fastas() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fseq" data-gbdraw-feature-id="fseq" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""

    enriched = enrich_svg(
        svg,
        InteractiveSvgContext(
            features=[
                {
                    "svg_id": "fseq",
                    "record_id": "rec1",
                    "type": "CDS",
                    "start": 0,
                    "end": 9,
                    "strand": "+",
                    "product": "protein A",
                    "qualifiers": {"protein_id": ["WP_000001.1"]},
                    "nucleotide_sequence": "ATGAAATAA",
                    "amino_acid_sequence": "MK",
                }
            ]
        ),
    )

    payload = _metadata_payload(enriched)
    feature = payload["features"][0]

    assert payload["schema"] == "gbdraw-interactive-feature-popup-v2"
    assert feature["nucleotide_sequence"] == "ATGAAATAA"
    assert feature["amino_acid_sequence"] == "MK"
    assert "nucleotide_fasta" not in feature
    assert "amino_acid_fasta" not in feature


def test_enrich_svg_v2_deduplicates_translation_sequence() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fseq" data-gbdraw-feature-id="fseq" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""
    payload = _metadata_payload(
        enrich_svg(
            svg,
            InteractiveSvgContext(
                features=[
                    {
                        "svg_id": "fseq",
                        "type": "CDS",
                        "qualifiers": {"translation": ["MPEPTIDE"]},
                        "amino_acid_sequence": "MPEPTIDE",
                    }
                ]
            ),
        )
    )

    feature = payload["features"][0]
    assert feature["qualifiers"]["translation"] == ["MPEPTIDE"]
    assert "amino_acid_sequence" not in feature


def test_enrich_svg_v2_metadata_is_at_least_35_percent_smaller_than_v1_fixture() -> None:
    source = (Path(__file__).parent / "test_inputs" / "AP027280_comparison.interactive.svg").read_text(
        encoding="utf-8"
    )
    v1_text = _metadata_text(source)
    v1_payload = json.loads(v1_text)
    enriched = enrich_svg(
        source,
        InteractiveSvgContext(
            features=v1_payload.get("features", []),
            orthogroups=v1_payload.get("orthogroups", []),
            popup_mode="rich",
        ),
    )
    v2_text = _metadata_text(enriched)

    assert _metadata_payload(enriched)["schema"] == "gbdraw-interactive-feature-popup-v2"
    assert len(v2_text.encode("utf-8")) <= len(v1_text.encode("utf-8")) * 0.65


def test_enrich_svg_matches_record_suffixed_session_feature_ids() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fseq" data-gbdraw-feature-id="fseq" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""

    enriched = enrich_svg(
        svg,
        InteractiveSvgContext(
            features=[
                {
                    "svg_id": "fseq_record_1",
                    "record_id": "rec1",
                    "type": "CDS",
                    "start": 0,
                    "end": 9,
                    "strand": "+",
                    "product": "protein A",
                    "nucleotide_sequence": "ATGAAATAA",
                    "amino_acid_sequence": "MK",
                }
            ]
        ),
    )

    feature = _metadata_payload(enriched)["features"][0]

    assert feature["svg_id"] == "fseq"
    assert feature["nucleotide_sequence"] == "ATGAAATAA"
    assert 'data-gbdraw-interactive-feature="true"' in enriched


def test_enrich_svg_rejects_malformed_svg() -> None:
    with pytest.raises(GbdrawError, match="Malformed SVG source"):
        enrich_svg("<svg><path></svg")


def test_records_metadata_uses_shared_feature_id_helper() -> None:
    feature = SeqFeature(
        SimpleLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"gene": ["geneA"], "translation": ["MK"]},
    )
    record = SeqRecord(Seq("ATGAAATAG"), id="rec1", features=[feature])

    payload = extract_features_from_records_payload([record], selected_features=["CDS"])

    assert payload["features"][0]["svg_id"] == compute_feature_hash(
        feature,
        record_id=record.id,
    )


def test_linear_cli_interactive_context_includes_orthogroup_payload() -> None:
    feature = SeqFeature(
        SimpleLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"gene": ["geneA"], "translation": ["MK"]},
    )
    record = SeqRecord(Seq("ATGAAATAG"), id="rec1", features=[feature])
    svg_id = compute_feature_hash(feature, record_id=record.id)
    member = SimpleNamespace(
        orthogroup_id="og_1",
        protein_id="gbd_r0001_cds000001",
        source_protein_id="WP_000001.1",
        record_index=0,
        record_id="rec1",
        feature_index=0,
        label="geneA",
        feature_svg_id=svg_id,
        start=0,
        end=9,
        strand=1,
        representative=True,
        role="anchor",
        confidence="high",
        assignment_reason="rbh",
        supporting_edges=(),
        best_core_support=0.0,
        second_best_core_support=0.0,
        gene="geneA",
        product="protein A",
        note=None,
        locus_tag=None,
        gene_id=None,
        old_locus_tag=None,
    )
    orthogroups = SimpleNamespace(
        orthogroups={"og_1": [member]},
        names_by_orthogroup_id={"og_1": "geneA"},
        descriptions_by_orthogroup_id={"og_1": "protein A"},
        name_candidates_by_orthogroup_id={},
        confidence_by_orthogroup_id={"og_1": "high"},
        rbh_orthogroups={"og_1": ("gbd_r0001_cds000001",)},
        ortholog_edges_by_orthogroup_id={},
        ortholog_paths_by_orthogroup_id={},
        related_edges_by_orthogroup_id={},
        scope_by_orthogroup_id={"og_1": "cross_record"},
        source_record_index_by_orthogroup_id={},
    )

    context = linear_cli._build_interactive_svg_context([record], {"CDS"}, orthogroups)

    assert context.orthogroups[0]["id"] == "og_1"
    assert context.orthogroups[0]["members"][0]["start"] == 0
    assert context.orthogroups[0]["members"][0]["end"] == 9
    assert context.orthogroups[0]["members"][0]["strand"] == "+"
    assert context.features[0]["orthogroup_id"] == "og_1"
    assert context.features[0]["orthogroup_member"]["sourceProteinId"] == "WP_000001.1"
    assert context.features[0]["orthogroup_member"]["strand"] == "+"


def test_linear_cli_interactive_context_orthogroup_members_use_absolute_region_coordinates() -> None:
    feature = SeqFeature(
        SimpleLocation(0, 457, strand=-1),
        type="CDS",
        qualifiers={"gene": ["geneA"], "translation": ["M"]},
    )
    record = SeqRecord(Seq("N" * 1000), id="rec1", features=[feature])
    record.annotations["gbdraw_coord_base"] = 10000
    record.annotations["gbdraw_coord_step"] = 1
    svg_id = compute_feature_hash_from_parts(
        "CDS",
        9999,
        10456,
        -1,
        record_id=record.id,
    )
    member = SimpleNamespace(
        orthogroup_id="og_1",
        protein_id="gbd_r0001_cds000001",
        source_protein_id="WP_000001.1",
        record_index=0,
        record_id="rec1",
        feature_index=0,
        label="geneA",
        feature_svg_id=svg_id,
        start=0,
        end=457,
        strand=-1,
        representative=True,
        role="anchor",
        confidence="high",
        assignment_reason="rbh",
        supporting_edges=(),
        best_core_support=0.0,
        second_best_core_support=0.0,
        gene="geneA",
        product="protein A",
        note=None,
        locus_tag=None,
        gene_id=None,
        old_locus_tag=None,
    )
    orthogroups = SimpleNamespace(
        orthogroups={"og_1": [member]},
        names_by_orthogroup_id={"og_1": "geneA"},
        descriptions_by_orthogroup_id={"og_1": "protein A"},
        name_candidates_by_orthogroup_id={},
        confidence_by_orthogroup_id={"og_1": "high"},
        rbh_orthogroups={"og_1": ("gbd_r0001_cds000001",)},
        ortholog_edges_by_orthogroup_id={},
        ortholog_paths_by_orthogroup_id={},
        related_edges_by_orthogroup_id={},
        scope_by_orthogroup_id={"og_1": "cross_record"},
        source_record_index_by_orthogroup_id={},
    )

    context = linear_cli._build_interactive_svg_context([record], {"CDS"}, orthogroups)

    member_payload = context.orthogroups[0]["members"][0]
    assert member_payload["start"] == 9999
    assert member_payload["end"] == 10456
    assert member_payload["strand"] == "-"
    assert context.features[0]["orthogroup_member"]["start"] == 9999
    assert context.features[0]["orthogroup_member"]["end"] == 10456


def test_enrich_svg_uses_orthogroup_payload_for_feature_and_match_metadata() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px" viewBox="0 0 100 80">
      <path id="fquery" data-gbdraw-feature-id="fquery" fill="#54bcf8" d="M 1 1 L 2 2" />
      <path data-gbdraw-pairwise-match-id="m1" data-match-kind="orthogroup"
        data-orthogroup-id="og_1" data-query-feature-svg-id="fquery"
        data-subject-feature-svg-id="fquery" d="M 1 3 L 2 4" />
    </svg>"""

    enriched = enrich_svg(
        svg,
        InteractiveSvgContext(
            features=[
                {
                    "svg_id": "fquery",
                    "record_id": "rec1",
                    "type": "CDS",
                    "orthogroup_id": "og_1",
                    "orthogroup_member_count": 1,
                    "orthogroup_record_coverage": 1,
                    "orthogroup_member": {
                        "featureSvgId": "fquery",
                        "sourceProteinId": "WP_000001.1",
                    },
                }
            ],
            orthogroups=[
                {
                    "id": "og_1",
                    "name": "geneA",
                    "description": "protein A",
                    "member_count": 1,
                    "record_coverage_count": 1,
                    "members": [
                        {
                            "featureSvgId": "fquery",
                            "recordId": "rec1",
                            "sourceProteinId": "WP_000001.1",
                            "start": 0,
                            "end": 9,
                            "strand": 1,
                            "product": "protein A",
                        }
                    ],
                }
            ],
        ),
    )

    payload = _metadata_payload(enriched)

    assert payload["features"][0]["orthogroup_id"] == "og_1"
    assert payload["orthogroups"][0]["id"] == "og_1"
    match = payload["matches"][0]
    assert match["match_kind"] == "orthogroup"
    assert match["orthogroup_ids"] == ["og_1"]
    assert "sections" not in match
    assert "hover_rows" not in match


def test_enrich_svg_uses_stable_orthogroup_member_feature_ids() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fstable" data-gbdraw-feature-id="fstable" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""
    member = {
        "orthogroupId": "og_1",
        "featureSvgId": "fstable_record_1",
        "stableFeatureSvgId": "fstable",
        "recordId": "rec1",
        "sourceProteinId": "WP_000001.1",
        "start": 0,
        "end": 9,
        "strand": "+",
    }

    payload = _metadata_payload(
        enrich_svg(
            svg,
            InteractiveSvgContext(
                features=[
                    {
                        "svg_id": "fstable_record_1",
                        "stable_svg_id": "fstable",
                        "record_id": "rec1",
                        "type": "CDS",
                        "start": 0,
                        "end": 9,
                        "strand": "+",
                        "orthogroup_id": "og_1",
                        "orthogroup_member": member,
                        "nucleotide_sequence": "ATGAAATAA",
                        "amino_acid_sequence": "MK*",
                    }
                ],
                orthogroups=[
                    {
                        "id": "og_1",
                        "member_count": 1,
                        "record_coverage_count": 1,
                        "members": [member],
                    }
                ],
            ),
        )
    )

    assert payload["features"][0]["svg_id"] == "fstable"
    assert payload["features"][0]["orthogroup_member"]["feature_svg_id"] == "fstable"
    assert payload["orthogroups"][0]["members"][0]["feature_svg_id"] == "fstable"


def test_enrich_svg_keeps_hidden_biological_features_and_orthogroup_members() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fvisible_record_1" data-gbdraw-feature-id="fvisible_record_1"
        data-gbdraw-stable-feature-id="fvisible" data-gbdraw-record-index="0"
        fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""
    visible_member = {
        "orthogroupId": "og_shared",
        "featureSvgId": "fvisible",
        "stableFeatureSvgId": "fvisible",
        "recordIndex": 0,
        "recordId": "rec1",
        "sourceProteinId": "WP_VISIBLE.1",
    }
    hidden_member = {
        "orthogroupId": "og_shared",
        "featureSvgId": "fvisible",
        "stableFeatureSvgId": "fvisible",
        "recordIndex": 1,
        "recordId": "rec2",
        "sourceProteinId": "WP_HIDDEN.1",
    }
    hidden_only_member = {
        "orthogroupId": "og_hidden",
        "featureSvgId": "fhiddenonly",
        "stableFeatureSvgId": "fhiddenonly",
        "recordIndex": 1,
        "recordId": "rec2",
        "sourceProteinId": "WP_HIDDEN_ONLY.1",
    }

    payload = _metadata_payload(
        enrich_svg(
            svg,
            InteractiveSvgContext(
                features=[
                    {
                        "svg_id": "fvisible",
                        "stable_svg_id": "fvisible",
                        "record_idx": 0,
                        "record_id": "rec1",
                        "type": "CDS",
                        "orthogroup_id": "og_shared",
                        "orthogroup_member": visible_member,
                    }
                ],
                biological_features=[
                    {
                        "svg_id": "fvisible",
                        "stable_svg_id": "fvisible",
                        "stable_feature_id": "fvisible",
                        "record_idx": 0,
                        "feature_index": 0,
                        "record_id": "rec1",
                        "type": "CDS",
                        "nucleotide_sequence": "ATGAAATAA",
                        "amino_acid_sequence": "MK",
                    },
                    {
                        "svg_id": "fvisible",
                        "stable_svg_id": "fvisible",
                        "stable_feature_id": "fvisible",
                        "record_idx": 1,
                        "feature_index": 0,
                        "record_id": "rec2",
                        "type": "CDS",
                        "nucleotide_sequence": "ATGCCCTAA",
                        "amino_acid_sequence": "MP",
                    },
                    {
                        "svg_id": "fhiddenonly",
                        "stable_svg_id": "fhiddenonly",
                        "stable_feature_id": "fhiddenonly",
                        "record_idx": 1,
                        "feature_index": 1,
                        "record_id": "rec2",
                        "type": "CDS",
                        "nucleotide_sequence": "ATGGGGTAA",
                        "amino_acid_sequence": "MG",
                    },
                ],
                orthogroups=[
                    {
                        "id": "og_shared",
                        "member_count": 2,
                        "record_coverage_count": 2,
                        "members": [visible_member, hidden_member],
                    },
                    {
                        "id": "og_hidden",
                        "member_count": 1,
                        "record_coverage_count": 1,
                        "members": [hidden_only_member],
                    },
                ],
            ),
        )
    )

    assert [feature["svg_id"] for feature in payload["features"]] == [
        "fvisible_record_1"
    ]
    biological_by_key = {
        (feature["record_idx"], feature["stable_feature_id"]): feature
        for feature in payload["biological_features"]
    }
    assert biological_by_key[(0, "fvisible")]["rendered_svg_id"] == "fvisible_record_1"
    assert "rendered_svg_id" not in biological_by_key[(1, "fvisible")]
    assert biological_by_key[(1, "fvisible")]["amino_acid_sequence"] == "MP"

    groups = {group["id"]: group for group in payload["orthogroups"]}
    assert set(groups) == {"og_shared", "og_hidden"}
    shared_members = {
        member["record_index"]: member for member in groups["og_shared"]["members"]
    }
    assert shared_members[0]["stable_feature_svg_id"] == "fvisible"
    assert shared_members[0]["rendered_feature_svg_id"] == "fvisible_record_1"
    assert shared_members[1]["feature_svg_id"] == "fvisible"
    assert "rendered_feature_svg_id" not in shared_members[1]
    assert groups["og_hidden"]["members"][0]["feature_svg_id"] == "fhiddenonly"
    assert "rendered_feature_svg_id" not in groups["og_hidden"]["members"][0]


def test_enrich_svg_match_metadata_matches_standalone_pairwise_sections() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px" viewBox="0 0 100 80">
      <path id="fquery" data-gbdraw-feature-id="fquery" fill="#54bcf8" d="M 1 1 L 2 2" />
      <path id="fsubject" data-gbdraw-feature-id="fsubject" fill="#54bcf8" d="M 3 1 L 4 2" />
      <path data-gbdraw-pairwise-match-id="m1" data-match-kind="pairwise"
        data-orthogroup-id="og_1" data-query-record-id="rec1"
        data-subject-record-id="rec2" data-qstart="1" data-qend="9"
        data-sstart="10" data-send="18" data-query-feature-svg-id="fquery"
        data-subject-feature-svg-id="fsubject" data-query-locus-id="geneA"
        data-subject-locus-id="geneB" data-query-display-name="Protein A"
        data-subject-display-name="Protein B" data-identity="99.1"
        data-alignment-length="9" d="M 1 3 L 2 4" />
    </svg>"""

    enriched = enrich_svg(
        svg,
        InteractiveSvgContext(
            features=[
                {
                    "svg_id": "fquery",
                    "record_id": "rec1",
                    "type": "CDS",
                    "start": 0,
                    "end": 9,
                    "strand": "+",
                    "orthogroup_id": "og_1",
                    "orthogroup_member_count": 2,
                    "orthogroup_record_coverage": 2,
                    "source_protein_id": "WP_000001.1",
                    "product": "Protein A",
                    "orthogroup_member": {
                        "featureSvgId": "fquery",
                        "recordId": "rec1",
                        "sourceProteinId": "WP_000001.1",
                        "start": 0,
                        "end": 9,
                        "strand": "+",
                        "product": "Protein A",
                    },
                },
                {
                    "svg_id": "fsubject",
                    "record_id": "rec2",
                    "type": "CDS",
                    "start": 9,
                    "end": 18,
                    "strand": "-",
                    "orthogroup_id": "og_1",
                    "orthogroup_member_count": 2,
                    "orthogroup_record_coverage": 2,
                    "source_protein_id": "WP_000002.1",
                    "product": "Protein B",
                    "orthogroup_member": {
                        "featureSvgId": "fsubject",
                        "recordId": "rec2",
                        "sourceProteinId": "WP_000002.1",
                        "start": 9,
                        "end": 18,
                        "strand": "-",
                        "product": "Protein B",
                    },
                },
            ],
            orthogroups=[
                {
                    "id": "og_1",
                    "name": "Protein family",
                    "description": "conserved protein",
                    "member_count": 2,
                    "record_coverage_count": 2,
                    "members": [
                        {
                            "featureSvgId": "fquery",
                            "recordId": "rec1",
                            "sourceProteinId": "WP_000001.1",
                            "start": 0,
                            "end": 9,
                            "strand": "+",
                            "product": "Protein A",
                        },
                        {
                            "featureSvgId": "fsubject",
                            "recordId": "rec2",
                            "sourceProteinId": "WP_000002.1",
                            "start": 9,
                            "end": 18,
                            "strand": "-",
                            "product": "Protein B",
                        },
                    ],
                }
            ],
        ),
    )

    match = _metadata_payload(enriched)["matches"][0]

    assert match["match_kind"] == "pairwise"
    assert match["orthogroup_ids"] == ["og_1"]
    assert match["query_record_id"] == "rec1"
    assert match["subject_record_id"] == "rec2"
    assert match["query_feature_svg_id"] == "fquery"
    assert match["subject_feature_svg_id"] == "fsubject"
    assert match["identity"] == "99.1"
    assert "sections" not in match


def test_enrich_svg_match_metadata_matches_standalone_collinear_sections() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px" viewBox="0 0 100 80">
      <path id="fquery" data-gbdraw-feature-id="fquery" fill="#54bcf8" d="M 1 1 L 2 2" />
      <path id="fsubject" data-gbdraw-feature-id="fsubject" fill="#54bcf8" d="M 3 1 L 4 2" />
      <path data-gbdraw-pairwise-match-id="m1" data-match-kind="collinear"
        data-orthogroup-id="og_1" data-collinearity-block-id="block_1"
        data-collinearity-block-kind="syntenic" data-collinearity-orientation="+"
        data-collinearity-color-mode="orthogroup" data-collinearity-block-score="42"
        data-collinearity-block-evalue="1e-10" data-collinearity-anchor-index="1"
        data-collinearity-anchor-count="3" data-query-record-id="rec1"
        data-subject-record-id="rec2" data-qstart="1" data-qend="9"
        data-sstart="10" data-send="18" data-query-feature-svg-id="fquery"
        data-subject-feature-svg-id="fsubject" data-identity="95.5"
        data-alignment-length="9" d="M 1 3 L 2 4" />
    </svg>"""

    context = InteractiveSvgContext(
        features=[
            {
                "svg_id": "fquery",
                "record_id": "rec1",
                "type": "CDS",
                "orthogroup_id": "og_1",
                "source_protein_id": "WP_000001.1",
                "orthogroup_member": {
                    "featureSvgId": "fquery",
                    "recordId": "rec1",
                    "sourceProteinId": "WP_000001.1",
                },
            },
            {
                "svg_id": "fsubject",
                "record_id": "rec2",
                "type": "CDS",
                "orthogroup_id": "og_1",
                "source_protein_id": "WP_000002.1",
                "orthogroup_member": {
                    "featureSvgId": "fsubject",
                    "recordId": "rec2",
                    "sourceProteinId": "WP_000002.1",
                },
            },
        ],
        orthogroups=[
            {
                "id": "og_1",
                "name": "Protein family",
                "member_count": 2,
                "record_coverage_count": 2,
                "members": [
                    {
                        "featureSvgId": "fquery",
                        "recordId": "rec1",
                        "sourceProteinId": "WP_000001.1",
                    },
                    {
                        "featureSvgId": "fsubject",
                        "recordId": "rec2",
                        "sourceProteinId": "WP_000002.1",
                    },
                ],
            }
        ],
    )

    match = _metadata_payload(enrich_svg(svg, context))["matches"][0]

    assert match["match_kind"] == "collinear"
    assert match["orthogroup_ids"] == ["og_1"]
    assert match["collinearity_block_id"] == "block_1"
    assert match["block_kind"] == "syntenic"
    assert match["identity"] == "95.5"
    assert match["alignment_length"] == "9"
    assert "block_orthogroups" not in match
    assert "sections" not in match


def test_save_figure_interactive_writes_static_and_interactive_without_cairosvg(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "out.svg")
    calls: list[str] = []

    def fake_enrich(svg_source: str, context=None) -> str:
        calls.append(svg_source)
        assert context is not None
        return "<svg data-gbdraw-interactive-svg=\"true\" />"

    monkeypatch.setattr(export_module, "enrich_svg", fake_enrich)
    monkeypatch.setattr(
        export_module,
        "get_cairosvg",
        lambda: pytest.fail("CairoSVG should not be loaded for interactive-svg only"),
    )

    export_module.save_figure(
        drawing,
        ["interactive-svg"],
        interactive_context=InteractiveSvgContext(features=[{"svg_id": "fabc12345"}]),
    )

    assert (tmp_path / "out.svg").exists()
    assert (tmp_path / "out.interactive.svg").read_text(encoding="utf-8").startswith(
        "<svg data-gbdraw-interactive-svg"
    )
    assert len(calls) == 1


def test_save_figure_interactive_png_uses_cairosvg_only_for_png(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "out.svg")
    monkeypatch.setattr(export_module, "enrich_svg", lambda _source, context=None: "<svg />")

    class FakeCairoSvg:
        @staticmethod
        def svg2png(*, bytestring, write_to):
            Path(write_to).write_bytes(bytestring)

    monkeypatch.setattr(export_module, "get_cairosvg", lambda: FakeCairoSvg)

    export_module.save_figure(drawing, ["interactive-svg", "png"])

    assert (tmp_path / "out.interactive.svg").exists()
    assert (tmp_path / "out.png").exists()


def test_save_figure_to_returns_interactive_path_and_checks_overwrite(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "ignored.svg")
    monkeypatch.setattr(api_render, "enrich_svg", lambda _source, context=None: "<svg />")

    paths = api_render.save_figure_to(
        drawing,
        "interactive_svg",
        output_dir=str(tmp_path),
        output_prefix="api-out",
    )

    assert paths == [
        str(tmp_path / "api-out.svg"),
        str(tmp_path / "api-out.interactive.svg"),
    ]
    with pytest.raises(ValidationError, match="already exist"):
        api_render.save_figure_to(
            drawing,
            ["interactive-svg"],
            output_dir=str(tmp_path),
            output_prefix="api-out",
        )


def test_save_figure_to_overwrite_replaces_existing_binary(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "ignored.svg")
    svg_path = tmp_path / "api-out.svg"
    png_path = tmp_path / "api-out.png"
    svg_path.write_text("stale svg", encoding="utf-8")
    png_path.write_bytes(b"stale png")

    class FakeCairoSvg:
        @staticmethod
        def svg2png(*, bytestring, write_to):
            Path(write_to).write_bytes(b"new png")

    monkeypatch.setattr(api_render._export, "get_cairosvg", lambda: FakeCairoSvg)

    paths = api_render.save_figure_to(
        drawing,
        "png",
        output_dir=str(tmp_path),
        output_prefix="api-out",
        overwrite=True,
    )

    assert paths == [str(svg_path), str(png_path)]
    assert svg_path.read_text(encoding="utf-8").startswith("<?xml")
    assert png_path.read_bytes() == b"new png"


def test_render_to_bytes_supports_interactive_svg(monkeypatch, tmp_path: Path) -> None:
    drawing = _drawing(tmp_path / "out.svg")
    context = InteractiveSvgContext(features=[{"svg_id": "fabc12345"}])

    def fake_enrich(_source, *, context=None):
        assert context is not None
        return "<svg interactive=\"yes\" />"

    monkeypatch.setattr(api_render, "enrich_svg", fake_enrich)

    assert api_render.render_to_bytes(
        drawing,
        "interactive_svg",
        interactive_context=context,
    ) == b"<svg interactive=\"yes\" />"


def test_save_figure_to_raises_when_cairosvg_is_unavailable(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "out.svg")

    def missing_cairosvg():
        raise ImportError("missing")

    monkeypatch.setattr(api_render._export, "get_cairosvg", missing_cairosvg)

    with pytest.raises(ValidationError, match="CairoSVG is not installed"):
        api_render.save_figure_to(
            drawing,
            "png",
            output_dir=str(tmp_path),
            output_prefix="strict",
        )

    assert (tmp_path / "strict.svg").exists()
    assert not (tmp_path / "strict.png").exists()


def test_save_figure_to_reports_partial_conversion_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "out.svg")

    class PartialCairoSvg:
        @staticmethod
        def svg2png(*, bytestring, write_to):
            Path(write_to).write_bytes(bytestring)

        @staticmethod
        def svg2pdf(*, bytestring, write_to):
            raise RuntimeError("converter failed")

    monkeypatch.setattr(api_render._export, "get_cairosvg", lambda: PartialCairoSvg)

    with pytest.raises(ExportError, match="Failed to generate PDF"):
        api_render.save_figure_to(
            drawing,
            ["png", "pdf"],
            output_dir=str(tmp_path),
            output_prefix="partial",
        )

    assert (tmp_path / "partial.svg").exists()
    assert (tmp_path / "partial.png").exists()
    assert not (tmp_path / "partial.pdf").exists()


def test_render_to_bytes_wraps_converter_failure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    drawing = _drawing(tmp_path / "out.svg")

    class BrokenCairoSvg:
        @staticmethod
        def svg2png(*, bytestring):
            raise RuntimeError("converter failed")

    monkeypatch.setattr(api_render._export, "get_cairosvg", lambda: BrokenCairoSvg)

    with pytest.raises(ExportError, match="Failed to render PNG bytes"):
        api_render.render_to_bytes(drawing, "png")


def test_cli_help_lists_interactive_svg(capsys) -> None:
    with pytest.raises(SystemExit):
        circular_cli._get_args(["--help"])
    circular_help = capsys.readouterr().out

    with pytest.raises(SystemExit):
        linear_cli._get_args(["--help"])
    linear_help = capsys.readouterr().out

    assert "interactive_svg" in circular_help
    assert "interactive_svg" in linear_help
