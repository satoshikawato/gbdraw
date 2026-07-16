from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.features.ids import compute_feature_hash
from gbdraw.io.genome import load_gff_fasta
from gbdraw.io.regions import apply_region_specs, parse_region_specs
from gbdraw.web_support.feature_metadata import extract_features_from_genbank_payload
from gbdraw.web_support.feature_metadata import extract_features_from_gff_fasta_payload
from gbdraw.features.visibility import compile_feature_visibility_rules, should_render_feature


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_HELPERS_PATH = REPO_ROOT / "gbdraw" / "web" / "js" / "app" / "python-helpers.js"


@pytest.fixture(scope="module")
def python_helpers_namespace() -> dict[str, object]:
    source = PYTHON_HELPERS_PATH.read_text(encoding="utf-8")
    helper_source = source.split("export const PYTHON_HELPERS = `", 1)[1].rsplit("\n`;", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)
    return namespace


def _write_genbank(tmp_path: Path, record: SeqRecord) -> Path:
    record.annotations["molecule_type"] = "DNA"
    path = tmp_path / "input.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _write_gff_fasta(tmp_path: Path) -> tuple[Path, Path]:
    gff_path = tmp_path / "input.gff3"
    fasta_path = tmp_path / "input.fasta"
    gff_path.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "##sequence-region GffRecord 1 30",
                "GffRecord\ttest\tgene\t1\t12\t.\t+\t.\tID=gene1;Name=example",
                "GffRecord\ttest\tCDS\t1\t12\t.\t+\t0\tID=cds1;Parent=gene1;gene=example;product=example%20protein",
                "",
            ]
        ),
        encoding="utf-8",
    )
    fasta_path.write_text(">GffRecord\nATGAAATAAGGGCCCCCCCCCCCCCCCCCC\n", encoding="utf-8")
    return gff_path, fasta_path


def _extract_features(namespace: dict[str, object], path: Path) -> list[dict[str, object]]:
    extract = namespace["extract_features_from_genbank"]
    result = extract(str(path))  # type: ignore[operator]
    payload = json.loads(result)
    assert "error" not in payload
    return payload["features"]


class JsNull:
    pass


class JsUndefined:
    pass


def test_python_helpers_blank_or_js_nullish_accepts_pyodide_sentinels(
    python_helpers_namespace: dict[str, object],
) -> None:
    is_blank = python_helpers_namespace["_is_blank_or_js_nullish"]

    for value in (None, JsNull(), JsUndefined(), "", "null", "undefined", "none"):
        assert is_blank(value) is True  # type: ignore[operator]

    for value in ("0", "12.5", "font"):
        assert is_blank(value) is False  # type: ignore[operator]


def test_regenerate_definition_svgs_accepts_nullish_font_sizes(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAA"), id="NC_000010", name="Nullish")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["NULLISH_001"], "translation": ["MK"]},
        )
    )
    path = _write_genbank(tmp_path, record)
    regenerate = python_helpers_namespace["regenerate_definition_svgs"]

    payload = json.loads(
        regenerate(  # type: ignore[operator]
            str(path),
            font_size=JsNull(),
            plot_title_font_size=JsUndefined(),
        )
    )

    assert "error" not in payload
    assert payload["definitions"]


def test_regenerate_definition_svgs_accepts_numeric_font_size_strings(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAA"), id="NC_000011", name="Numeric")
    path = _write_genbank(tmp_path, record)
    regenerate = python_helpers_namespace["regenerate_definition_svgs"]

    payload = json.loads(
        regenerate(  # type: ignore[operator]
            str(path),
            font_size="12.5",
            plot_title_font_size="16",
            plot_title="Numeric title",
            plot_title_position="top",
        )
    )

    assert "error" not in payload
    assert len(payload["definitions"]) == 2


def test_regenerate_definition_svgs_reports_invalid_font_size_strings(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAA"), id="NC_000012", name="InvalidNumeric")
    path = _write_genbank(tmp_path, record)
    regenerate = python_helpers_namespace["regenerate_definition_svgs"]

    payload = json.loads(regenerate(str(path), font_size="not-a-number"))  # type: ignore[operator]

    assert "error" in payload
    assert "not-a-number" in payload["error"] or "could not convert" in payload["error"]


def test_importable_feature_metadata_matches_pyodide_wrapper(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAA"), id="NC_000000", name="Importable")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["ABC_0000"],
                "product": ["importable protein"],
                "translation": ["MK"],
            },
        )
    )
    path = _write_genbank(tmp_path, record)
    wrapper_payload = json.loads(python_helpers_namespace["extract_features_from_genbank"](str(path)))  # type: ignore[operator]

    assert wrapper_payload == extract_features_from_genbank_payload(path)


def test_gff_fasta_feature_metadata_includes_nested_rendered_features(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    gff_path, fasta_path = _write_gff_fasta(tmp_path)

    payload = extract_features_from_gff_fasta_payload(
        gff_path,
        fasta_path,
        mode="linear",
        selected_features=["CDS"],
    )
    features = payload["features"]

    assert len(features) == 1
    assert features[0]["type"] == "CDS"
    assert features[0]["record_id"] == "GffRecord"
    assert features[0]["gene"] == "example"
    assert features[0]["product"] == "example protein"
    assert features[0]["nucleotide_sequence"] == "ATGAAATAAGGG"

    rendered_records = load_gff_fasta(
        [str(gff_path)],
        [str(fasta_path)],
        mode="linear",
        selected_features_set={"CDS"},
    )
    assert features[0]["svg_id"] == compute_feature_hash(
        rendered_records[0].features[0],
        record_id=rendered_records[0].id,
    )

    wrapper_payload = json.loads(
        python_helpers_namespace["extract_features_from_gff_fasta"](  # type: ignore[operator]
            str(gff_path),
            str(fasta_path),
            "linear",
            None,
            None,
            None,
            json.dumps(["CDS"]),
        )
    )
    assert wrapper_payload == payload


def test_web_feature_extraction_includes_qualifiers_locations_and_translation(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAAGGGCCC"), id="NC_000001", name="TestRecord")
    record.annotations["organism"] = "Example organism"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["ABC_0001"],
                "product": ["example protein"],
                "note": ["first note", "second note"],
                "translation": ["MK"],
            },
        )
    )

    features = _extract_features(python_helpers_namespace, _write_genbank(tmp_path, record))
    feature = features[0]

    assert feature["record_id"] == "NC_000001"
    assert feature["organism"] == "Example organism"
    assert feature["type"] == "CDS"
    assert feature["start"] == 0
    assert feature["end"] == 9
    assert feature["strand"] == "+"
    assert feature["qualifiers"]["note"] == ["first note", "second note"]
    assert feature["selector"]["hash"] == feature["svg_id"]
    assert feature["selector"]["location"] == "0..9"
    assert feature["selector"]["record_location"] == "NC_000001:0..9:+"
    assert feature["selector"]["qualifiers"]["locus_tag"] == ["ABC_0001"]
    assert feature["location_parts"] == [{"start": 0, "end": 9, "strand": "+", "display": "1..9"}]
    assert feature["nucleotide_sequence"] == "ATGAAATAA"
    assert feature["amino_acid_sequence"] == "MK"
    assert feature["sequence_warnings"] == []


def test_web_feature_extraction_translates_simple_cds_without_translation(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAA"), id="NC_000002", name="Fallback")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 6, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["ABC_0002"], "codon_start": ["1"], "transl_table": ["11"]},
        )
    )

    feature = _extract_features(python_helpers_namespace, _write_genbank(tmp_path, record))[0]

    assert feature["nucleotide_sequence"] == "ATGAAA"
    assert feature["amino_acid_sequence"] == "MK"
    assert feature["sequence_warnings"] == []


def test_web_feature_extraction_reports_invalid_translation_warning(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAA"), id="NC_000003", name="Invalid")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 5, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["ABC_0003"]},
        )
    )

    feature = _extract_features(python_helpers_namespace, _write_genbank(tmp_path, record))[0]

    assert feature["nucleotide_sequence"] == "ATGAA"
    assert feature["amino_acid_sequence"] == ""
    assert any("not divisible by 3" in warning for warning in feature["sequence_warnings"])


def test_web_feature_extraction_adds_compound_location_parts(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("AAACCCGGGTTT"), id="NC_000004", name="Compound")
    record.features.append(
        SeqFeature(
            CompoundLocation([
                FeatureLocation(0, 3, strand=1),
                FeatureLocation(6, 9, strand=1),
            ]),
            type="misc_feature",
            qualifiers={"note": ["joined feature"]},
        )
    )

    feature = _extract_features(python_helpers_namespace, _write_genbank(tmp_path, record))[0]

    assert feature["location_parts"] == [
        {"start": 0, "end": 3, "strand": "+", "display": "1..3"},
        {"start": 6, "end": 9, "strand": "+", "display": "7..9"},
    ]
    assert feature["nucleotide_sequence"] == "AAAGGG"


def test_web_feature_extraction_region_uses_absolute_display_coordinates(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("A" * 60), id="NC_REGION", name="Region")
    record.features.append(
        SeqFeature(
            FeatureLocation(9, 20, strand=1),
            type="misc_feature",
            qualifiers={"note": ["visible in cropped region"]},
        )
    )
    path = _write_genbank(tmp_path, record)

    payload = extract_features_from_genbank_payload(path, region_spec="10-30")
    feature = payload["features"][0]
    cropped_record = apply_region_specs([record], parse_region_specs(["10-30"]))[0]

    assert feature["start"] == 9
    assert feature["end"] == 20
    assert feature["location_parts"] == [
        {"start": 9, "end": 20, "strand": "+", "display": "10..20"}
    ]
    assert feature["svg_id"] == compute_feature_hash(
        cropped_record.features[0],
        record_id=cropped_record.id,
    )


def test_web_feature_extraction_region_reverse_uses_absolute_display_coordinates(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("A" * 120), id="NC_REGION_RC", name="RegionRc")
    record.features.append(
        SeqFeature(
            FeatureLocation(89, 99, strand=1),
            type="misc_feature",
            qualifiers={"note": ["visible in reverse-complement region"]},
        )
    )
    path = _write_genbank(tmp_path, record)

    payload = extract_features_from_genbank_payload(path, region_spec="80-100:rc")
    feature = payload["features"][0]

    assert feature["start"] == 89
    assert feature["end"] == 99
    assert feature["strand"] == "-"
    assert feature["location_parts"] == [
        {"start": 89, "end": 99, "strand": "-", "display": "90..99"}
    ]


def test_web_feature_selector_record_location_matches_visibility_rule(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("ATGAAATAA"), id="NC_000005", name="SelectorRule")
    feature = SeqFeature(
        FeatureLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["ABC_0005"], "translation": ["MK"]},
    )
    record.features.append(feature)

    payload = extract_features_from_genbank_payload(_write_genbank(tmp_path, record))
    extracted = payload["features"][0]
    record_location = extracted["selector"]["record_location"]
    rules = compile_feature_visibility_rules(
        pd.DataFrame(
            [["NC_000005", "CDS", "record_location", f"^{re.escape(record_location)}$", "off"]],
            columns=["record_id", "feature_type", "qualifier", "value", "action"],
        )
    )

    assert (
        should_render_feature(
            feature,
            selected_features_set=["CDS"],
            feature_visibility_rules=rules,
            record_id=record.id,
        )
        is False
    )


def test_web_feature_selector_safety_scope_precedes_visibility_filtering(
    tmp_path: Path,
) -> None:
    record = SeqRecord(Seq("A" * 120), id="NC_000006", name="SafetyScope")
    record.features.extend(
        [
            SeqFeature(
                FeatureLocation(0, 30, strand=1),
                type="CDS",
                qualifiers={"locus_tag": ["HIDDEN_0006"], "translation": ["M"]},
            ),
            SeqFeature(
                FeatureLocation(60, 90, strand=-1),
                type="misc_feature",
                qualifiers={"note": ["outside selected type"]},
            ),
        ]
    )
    rules = compile_feature_visibility_rules(
        pd.DataFrame(
            [["NC_000006", "CDS", "locus_tag", "^HIDDEN_0006$", "off"]],
            columns=["record_id", "feature_type", "qualifier", "value", "action"],
        )
    )

    payload = extract_features_from_genbank_payload(
        _write_genbank(tmp_path, record),
        selected_features=["CDS"],
        feature_visibility_rules=rules,
    )

    assert payload["features"] == []
    scope_types = [entry["feature_type"] for entry in payload["selector_safety_scope"]]
    assert scope_types == ["CDS", "misc_feature"]
    assert payload["selector_safety_scope"][0]["selector"]["qualifiers"]["locus_tag"] == ["HIDDEN_0006"]
