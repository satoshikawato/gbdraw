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
from gbdraw import circular as circular_cli
from gbdraw import linear as linear_cli
from gbdraw.cli_utils import common as cli_common
from gbdraw.exceptions import GbdrawError, ValidationError
from gbdraw.render import export as export_module
from gbdraw.render.formats import (
    INTERACTIVE_SVG_FORMAT,
    SVG_FORMAT,
    classify_formats,
    parse_format_string,
    resolve_format_output_path,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg
from gbdraw.features.ids import compute_feature_hash
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload


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


def _script_payload(svg_source: str) -> str:
    root = ET.fromstring(svg_source)
    script = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "script"
        and element.get("id") == "gbdraw-interactive-feature-script"
    )
    return script.text or ""


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
        "interactive-svg"
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


def test_enrich_svg_generates_fallback_feature_payload() -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="100px" height="80px">
      <path id="fabc12345__part1" data-gbdraw-feature-id="fabc12345" fill="#54bcf8" d="M 1 1 L 2 2" />
    </svg>"""

    payload = _metadata_payload(enrich_svg(svg))

    assert payload["features"][0]["svg_id"] == "fabc12345"
    assert payload["features"][0]["fill_color"] == "#54bcf8"


def test_enrich_svg_embeds_feature_sequences_and_fastas() -> None:
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

    feature = _metadata_payload(enriched)["features"][0]

    assert feature["nucleotide_sequence"] == "ATGAAATAA"
    assert feature["amino_acid_sequence"] == "MK"
    assert feature["nucleotide_fasta"] == ">rec1:1-9 protein A\nATGAAATAA"
    assert feature["amino_acid_fasta"] == ">WP_000001.1 protein A\nMK"


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
    assert context.features[0]["orthogroup_id"] == "og_1"
    assert context.features[0]["orthogroup_member"]["sourceProteinId"] == "WP_000001.1"


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
    assert match["title"] == "og_1:geneA"
    assert match["sections"][0]["rows"][0] == ["Orthogroup ID", "og_1"]
    assert match["sections"][0]["member_rows"][0]["proteinId"] == "WP_000001.1"


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

    assert match["title"] == "Pairwise match"
    assert [section["title"] for section in match["sections"]] == [
        "Summary",
        "Alignment",
        "Orthogroup",
        "Query",
        "Subject",
    ]
    orthogroup_section = match["sections"][2]
    assert ["Orthogroup ID", "og_1"] in orthogroup_section["rows"]
    assert orthogroup_section["member_rows"][1]["proteinId"] == "WP_000002.1"
    query_section = match["sections"][3]
    assert query_section["feature_rows"][0]["proteinId"] == "WP_000001.1"
    assert query_section["featureRows"] == query_section["feature_rows"]


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

    assert match["title"] == "Collinearity block"
    assert [section["title"] for section in match["sections"]] == [
        "Summary",
        "Orthogroups covered",
        "Collinearity",
        "Query",
        "Subject",
    ]
    assert match["block_orthogroup_count"] == 1
    assert match["block_orthogroups"][0]["query_member"] == "WP_000001.1"
    collinearity_rows = match["sections"][2]["rows"]
    assert ["Block ID", "block_1"] in collinearity_rows
    assert ["Average identity", "95.5"] in collinearity_rows


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


def test_render_to_bytes_supports_interactive_svg(monkeypatch, tmp_path: Path) -> None:
    drawing = _drawing(tmp_path / "out.svg")
    monkeypatch.setattr(api_render, "enrich_svg", lambda _source: "<svg interactive=\"yes\" />")

    assert api_render.render_to_bytes(drawing, "interactive_svg") == b"<svg interactive=\"yes\" />"


def test_cli_help_lists_interactive_svg(capsys) -> None:
    with pytest.raises(SystemExit):
        circular_cli._get_args(["--help"])
    circular_help = capsys.readouterr().out

    with pytest.raises(SystemExit):
        linear_cli._get_args(["--help"])
    linear_help = capsys.readouterr().out

    assert "interactive-svg" in circular_help
    assert "interactive-svg" in linear_help
