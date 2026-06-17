from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.web_support.feature_metadata import extract_features_from_genbank_payload


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


def _extract_features(namespace: dict[str, object], path: Path) -> list[dict[str, object]]:
    extract = namespace["extract_features_from_genbank"]
    result = extract(str(path))  # type: ignore[operator]
    payload = json.loads(result)
    assert "error" not in payload
    return payload["features"]


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


def test_web_feature_extraction_includes_qualifiers_locations_and_translation(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("ATGAAATAAGGGCCC"), id="NC_000001", name="TestRecord")
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
    assert feature["type"] == "CDS"
    assert feature["start"] == 0
    assert feature["end"] == 9
    assert feature["strand"] == "+"
    assert feature["qualifiers"]["note"] == ["first note", "second note"]
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
