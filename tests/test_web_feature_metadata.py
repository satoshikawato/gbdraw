from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.features.ids import compute_feature_hash
from gbdraw.io.genome import load_gff_fasta
from gbdraw.io.regions import apply_region_specs, parse_region_specs
from gbdraw.io.record_select import reverse_records
from gbdraw.web_support.feature_metadata import _source_anchor_profile
from gbdraw.web_support.feature_metadata import extract_features_from_genbank_payload
from gbdraw.web_support.feature_metadata import extract_features_from_gff_fasta_payload
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload
from gbdraw.features.visibility import compile_feature_visibility_rules, should_render_feature


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_HELPERS_PATH = REPO_ROOT / "gbdraw" / "web" / "js" / "app" / "python-helpers.js"
MG1655_START_CODON_SUBSET = (
    REPO_ROOT / "tests" / "fixtures" / "cds_translation" / "mg1655_start_codon_subset.gb"
)


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


def test_canonical_request_wrapper_rejects_invalid_request_without_monkeypatching(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    from gbdraw.diagrams.linear import assemble

    marker = tmp_path / "keep.svg"
    marker.write_text("<svg />", encoding="utf-8")
    original_loader = assemble.load_comparisons
    wrapper = python_helpers_namespace["run_canonical_request_wrapper"]

    payload = wrapper("{}", "{}", str(tmp_path / "render"))  # type: ignore[operator]

    assert payload["error"]["code"] == "RESOURCE_INVALID"
    assert payload["error"]["stage"] == "resource-staging"
    assert set(payload["error"]) == {"code", "operation", "stage", "context"}
    assert assemble.load_comparisons is original_loader
    assert marker.exists()


def test_canonical_request_wrapper_returns_svg_content_as_utf8_bytes(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    wrapper = python_helpers_namespace["run_canonical_request_wrapper"]
    original_renderer = python_helpers_namespace["render_staged_canonical_web_request"]
    python_helpers_namespace["render_staged_canonical_web_request"] = (
        lambda *_args, **_kwargs: {
            "results": [{"name": "diagram.svg", "content": "<svg>µ</svg>"}],
            "metadata": {"featureCatalog": {"schema": 1}},
        }
    )
    try:
        payload = wrapper("{}", "{}", str(tmp_path / "render"))  # type: ignore[operator]
    finally:
        python_helpers_namespace["render_staged_canonical_web_request"] = original_renderer

    assert payload == {
        "results": [{"name": "diagram.svg", "content": "<svg>µ</svg>".encode()}],
        "metadata": b'{"featureCatalog":{"schema":1}}',
    }


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
    assert feature["anchorProfile"] == {
        "precision": "exact",
        "operator": "single",
        "partOrder": "biological",
        "strand": "+",
    }
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


def _write_gff3_cds_rows(
    tmp_path: Path,
    sequence: str,
    rows: list[tuple[int, int, str, str, str]],
) -> tuple[Path, Path]:
    """Write one GFF3 contig whose rows are (start, end, strand, phase, attributes)."""

    gff_path = tmp_path / "cds.gff3"
    fasta_path = tmp_path / "cds.fasta"
    lines = ["##gff-version 3", f"##sequence-region ctg1 1 {len(sequence)}"]
    lines.extend(
        f"ctg1\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attributes}"
        for start, end, strand, phase, attributes in rows
    )
    gff_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    fasta_path.write_text(f">ctg1\n{sequence}\n", encoding="utf-8")
    return gff_path, fasta_path


@pytest.mark.parametrize(
    ("sequence", "transl_table", "expected"),
    [
        ("GTGAAATAA", "11", "MK"),
        ("TTGAAATAA", "11", "MK"),
        ("ATTAAATAA", "11", "MK"),
        ("TTGAAATAA", None, "MK"),
        ("CTGAAATAA", None, "MK"),
        ("GTGAAATAA", None, "VK"),
        ("ATTAAATAA", None, "IK"),
    ],
)
def test_web_feature_extraction_translates_table_start_codon_as_methionine(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
    sequence: str,
    transl_table: str | None,
    expected: str,
) -> None:
    qualifiers = {"locus_tag": ["START_0001"]}
    if transl_table is not None:
        qualifiers["transl_table"] = [transl_table]
    record = SeqRecord(Seq(sequence), id="NC_000010", name="Start")
    record.features.append(
        SeqFeature(FeatureLocation(0, len(sequence), strand=1), type="CDS", qualifiers=qualifiers)
    )

    feature = _extract_features(python_helpers_namespace, _write_genbank(tmp_path, record))[0]

    assert feature["amino_acid_sequence"] == expected
    assert feature["sequence_warnings"] == []


def test_web_feature_extraction_translates_offset_or_five_prime_partial_cds_literally(
    tmp_path: Path,
    python_helpers_namespace: dict[str, object],
) -> None:
    record = SeqRecord(Seq("AGTGAAATAA"), id="NC_000011", name="Offset")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 10, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["OFFSET_0001"], "codon_start": ["2"], "transl_table": ["11"]},
        )
    )
    offset_feature = _extract_features(
        python_helpers_namespace,
        _write_genbank(tmp_path, record),
    )[0]
    assert offset_feature["amino_acid_sequence"] == "VK"

    # The 5' end of a plus-strand GFF3 CDS is its start column; of a minus-strand
    # CDS, its end column. Only an incomplete 5' end keeps the first codon literal.
    gff_path, fasta_path = _write_gff3_cds_rows(
        tmp_path,
        "GTGAAATAA" + "TTATTTCAC" + "GTGAAACCC" + "GGGTTTCAC",
        [
            (1, 9, "+", "0", "ID=plus5;transl_table=11;partial=true;start_range=.,1"),
            (10, 18, "-", "0", "ID=minus5;transl_table=11;partial=true;end_range=18,."),
            (19, 27, "+", "0", "ID=plus3;transl_table=11;partial=true;end_range=27,."),
            (28, 36, "-", "0", "ID=minus3;transl_table=11;partial=true;start_range=.,28"),
        ],
    )
    payload = extract_features_from_gff_fasta_payload(
        gff_path,
        fasta_path,
        selected_features=["CDS"],
    )

    assert {
        feature["qualifiers"]["id"][0]: feature["amino_acid_sequence"]
        for feature in payload["features"]
    } == {"plus5": "VK", "minus5": "VK", "plus3": "MKP", "minus3": "MKP"}


def test_web_feature_extraction_uses_gff3_phase_as_reading_frame(tmp_path: Path) -> None:
    gff_path, fasta_path = _write_gff3_cds_rows(
        tmp_path,
        "CATGAAACCCGGGTAA" + "TTATTTCACGC",
        [
            (1, 16, "+", "1", "ID=plus;transl_table=11;partial=true;start_range=.,1"),
            (17, 27, "-", "2", "ID=minus;transl_table=11;partial=true;end_range=27,."),
        ],
    )

    payload = extract_features_from_gff_fasta_payload(
        gff_path,
        fasta_path,
        selected_features=["CDS"],
    )

    features = {feature["qualifiers"]["id"][0]: feature for feature in payload["features"]}
    assert features["plus"]["amino_acid_sequence"] == "MKPG"
    assert features["minus"]["amino_acid_sequence"] == "VK"
    assert features["plus"]["sequence_warnings"] == []
    assert features["minus"]["sequence_warnings"] == []


def test_web_feature_extraction_matches_mg1655_translation_without_the_qualifier() -> None:
    record = SeqIO.read(MG1655_START_CODON_SUBSET, "genbank")
    expected: dict[str, str] = {}
    for feature in record.features:
        if feature.type == "CDS":
            expected[feature.qualifiers["locus_tag"][0]] = feature.qualifiers.pop("translation")[0]

    payload = extract_features_from_records_payload([record], selected_features=["CDS"])

    assert {
        feature["locus_tag"]: feature["amino_acid_sequence"]
        for feature in payload["features"]
        if feature["type"] == "CDS"
    } == expected


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


@pytest.mark.parametrize(
    ("location", "expected"),
    [
        (
            FeatureLocation(2, 8, strand=1),
            {"precision": "exact", "operator": "single", "partOrder": "biological", "strand": "+"},
        ),
        (
            FeatureLocation(2, 8, strand=-1),
            {"precision": "exact", "operator": "single", "partOrder": "biological", "strand": "-"},
        ),
        (
            CompoundLocation(
                [FeatureLocation(90, 100, strand=1), FeatureLocation(0, 5, strand=1)],
                operator="join",
            ),
            {"precision": "exact", "operator": "join", "partOrder": "biological", "strand": "+"},
        ),
        (
            CompoundLocation(
                [FeatureLocation(40, 50, strand=-1), FeatureLocation(10, 20, strand=-1)],
                operator="join",
            ),
            {"precision": "exact", "operator": "join", "partOrder": "biological", "strand": "-"},
        ),
        (
            FeatureLocation(2, 8, strand=None),
            {"precision": "exact", "operator": "single", "partOrder": "source-forward", "strand": "unstranded"},
        ),
        (
            CompoundLocation(
                [FeatureLocation(2, 8, strand=None), FeatureLocation(12, 16, strand=None)],
                operator="join",
            ),
            {"precision": "exact", "operator": "join", "partOrder": "source-forward", "strand": "unstranded"},
        ),
    ],
)
def test_source_anchor_profile_classifies_exact_source_paths(
    location: object,
    expected: dict[str, str],
) -> None:
    assert _source_anchor_profile(SeqFeature(location, type="misc_feature")) == expected


@pytest.mark.parametrize(
    ("location", "expected"),
    [
        (
            CompoundLocation(
                [FeatureLocation(2, 8, strand=1), FeatureLocation(12, 16, strand=-1)],
                operator="join",
            ),
            {"precision": "exact", "operator": "join", "partOrder": "ambiguous", "strand": "mixed"},
        ),
        (
            FeatureLocation(BeforePosition(2), 8, strand=1),
            {"precision": "fuzzy", "operator": "single", "partOrder": "biological", "strand": "+"},
        ),
        (
            CompoundLocation(
                [FeatureLocation(2, 8, strand=1), FeatureLocation(12, 16, strand=1)],
                operator="order",
            ),
            {"precision": "exact", "operator": "order", "partOrder": "ambiguous", "strand": "+"},
        ),
        (
            CompoundLocation(
                [FeatureLocation(2, 8, strand=1), FeatureLocation(12, 16, strand=1)],
                operator="bond",
            ),
            {"precision": "exact", "operator": "unknown", "partOrder": "ambiguous", "strand": "+"},
        ),
    ],
)
def test_source_anchor_profile_classifies_unsafe_locations_conservatively(
    location: object,
    expected: dict[str, str],
) -> None:
    assert _source_anchor_profile(SeqFeature(location, type="misc_feature")) == expected


def test_source_anchor_profile_survives_reverse_coordinate_mapping() -> None:
    record = SeqRecord(Seq("A" * 100), id="NC_ANCHOR_PROFILE")
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [FeatureLocation(90, 100, strand=1), FeatureLocation(0, 5, strand=1)],
                operator="join",
            ),
            type="misc_feature",
        )
    )

    reversed_feature = reverse_records([record], True)[0].features[0]

    assert reversed_feature.location.strand == -1
    assert _source_anchor_profile(reversed_feature) == {
        "precision": "exact",
        "operator": "join",
        "partOrder": "biological",
        "strand": "+",
    }


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
        record.features[0],
        record_id=record.id,
    )
    assert feature["rendered_feature_svg_id"] == compute_feature_hash(
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
    assert feature["strand"] == "+"
    assert feature["location_parts"] == [
        {"start": 89, "end": 99, "strand": "+", "display": "90..99"}
    ]
    assert feature["svg_id"] == compute_feature_hash(
        record.features[0],
        record_id=record.id,
    )


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


def test_web_feature_payload_records_the_drawn_selector_values(
    tmp_path: Path,
) -> None:
    # Feature catalog 5 (design Q4 OV-02): each rendered feature carries the
    # selector values of the record it was drawn from; no selector-safety scope
    # of the whole source is built.
    record = SeqRecord(Seq("A" * 120), id="NC_000006", name="Drawn")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 30, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["DRAWN_0006"], "translation": ["M"]},
        )
    )

    payload = extract_features_from_genbank_payload(
        _write_genbank(tmp_path, record),
        selected_features=["CDS"],
    )

    assert sorted(payload) == ["features", "record_ids"]
    [feature] = payload["features"]
    assert feature["drawn_selector"] == {
        "hash": feature["selector"]["hash"],
        "location": "0..30",
        "recordLocation": "NC_000006:0..30:+",
    }


def test_biological_feature_catalog_keeps_features_excluded_from_rendering() -> None:
    hidden_cds = SeqFeature(
        FeatureLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={
            "locus_tag": ["HIDDEN_0007"],
            "protein_id": ["WP_HIDDEN.1"],
            "translation": ["MK"],
        },
    )
    visible_trna = SeqFeature(
        FeatureLocation(12, 21, strand=-1),
        type="tRNA",
        qualifiers={"locus_tag": ["VISIBLE_0007"]},
    )
    record = SeqRecord(
        Seq("ATGAAATAAGGGTTTCCCAAA"),
        id="NC_000007",
        features=[hidden_cds, visible_trna],
    )

    legacy_payload = extract_features_from_records_payload(
        [record],
        selected_features=["tRNA"],
    )
    payload = extract_features_from_records_payload(
        [record],
        selected_features=["tRNA"],
        include_biological_features=True,
    )

    assert "biological_features" not in legacy_payload
    assert payload["features"] == legacy_payload["features"]
    assert [feature["type"] for feature in payload["features"]] == ["tRNA"]
    assert [feature["type"] for feature in payload["biological_features"]] == [
        "CDS",
        "tRNA",
    ]
    hidden = payload["biological_features"][0]
    assert hidden["feature_index"] == 0
    assert hidden["svg_id"] == hidden["stable_feature_id"] == hidden["stable_svg_id"]
    assert hidden["nucleotide_sequence"] == "ATGAAATAA"
    assert hidden["amino_acid_sequence"] == "MK"


def test_biological_feature_catalog_excludes_non_rendered_source_features() -> None:
    source = SeqFeature(
        FeatureLocation(0, 12, strand=1),
        type="source",
        qualifiers={"organism": ["Catalog organism"]},
    )
    cds = SeqFeature(
        FeatureLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["CATALOG_001"], "translation": ["MK"]},
    )
    record = SeqRecord(
        Seq("ATGAAATAAGGG"),
        id="catalog-record",
        features=[source, cds],
    )

    payload = extract_features_from_records_payload(
        [record],
        selected_features=["CDS"],
        include_biological_features=True,
    )

    assert [feature["type"] for feature in payload["features"]] == ["CDS"]
    assert [feature["type"] for feature in payload["biological_features"]] == [
        "CDS"
    ]

    source_payload = extract_features_from_records_payload(
        [record],
        selected_features=["source"],
        include_biological_features=True,
    )

    assert [feature["type"] for feature in source_payload["features"]] == [
        "source"
    ]
    assert [
        feature["type"] for feature in source_payload["biological_features"]
    ] == ["source", "CDS"]


def test_biological_feature_catalog_keeps_a_source_feature_with_its_own_visibility() -> None:
    # R-5 (Owner decision 2026-10-05): a feature with its own Feature
    # visibility stays in the Web Features list, which lists the catalog's
    # biological features, so it can be shown again.
    from gbdraw.features.overrides import ResolvedFeatureOverride
    from gbdraw.features.placement import ResolvedRecordFeatureInputs

    source = SeqFeature(FeatureLocation(0, 12, strand=1), type="source", qualifiers={"organism": ["x"]})
    cds = SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS", qualifiers={"locus_tag": ["C1"]})
    record = SeqRecord(Seq("ATGAAATAAGGG"), id="catalog-record", features=[source, cds])
    hidden = ResolvedFeatureOverride("off", None, None, "source 1..12")
    payload = extract_features_from_records_payload(
        [record],
        selected_features=["CDS"],
        record_features=(ResolvedRecordFeatureInputs("k", overrides={0: hidden}),),
        include_biological_features=True,
    )
    assert [feature["type"] for feature in payload["features"]] == ["CDS"]
    assert [feature["type"] for feature in payload["biological_features"]] == ["source", "CDS"]
