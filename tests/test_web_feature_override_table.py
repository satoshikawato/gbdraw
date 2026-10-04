"""Load Feature Edits TSV: the Worker helper reads a feature override table
against the committed Web request (design Q4 6.4, Owner Q2 = A and Q3 = A)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput, RecordPresentation
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import parse_region_spec
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.error_adapter import serialize_web_error
from gbdraw.web_support.feature_override_table import read_feature_override_table_json

FIXTURE = Path(__file__).parent / "fixtures" / "web_batch_two_records.gb"
HEADER = "record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text\n"


@pytest.fixture
def committed(tmp_path):
    """A committed Linear request: TESTA cropped, TESTB reverse-complemented, TESTA again."""
    texts = [f"{chunk.strip()}\n//\n" for chunk in FIXTURE.read_text().split("\n//") if chunk.strip()]
    paths = []
    for name, text in zip(("testa.gb", "testb.gb"), texts, strict=True):
        paths.append(tmp_path / name)
        paths[-1].write_text(text)
    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(paths[0]), record_key="linear-seq-a",
                    region=parse_region_spec("TESTA:201-3800")),
        RecordInput(GenBankInputSource(paths[1]), record_key="linear-seq-b",
                    presentation=RecordPresentation(reverse_complement=True)),
        RecordInput(GenBankInputSource(paths[0]), record_key="linear-seq-c"),
    ), options=LinearDiagramOptions(selected_features_set=("CDS", "tRNA", "misc_feature")))
    encoded = encode_canonical_request(request)
    catalogs = [item.source_feature_catalog for item in plan_request(request).provenance]
    return {
        "request": json.dumps(encoded.payload),
        "resources": json.dumps({item.resource_id: str(item.source_path) for item in encoded.resources}),
        "ids": [{
            dict(entry.qualifiers).get("locus_tag", (entry.feature_type,))[0]: entry.biological_feature_id
            for entry in catalog
        } for catalog in catalogs],
        "workspace": tmp_path,
    }


def _read(committed, table_text):
    table = committed["workspace"] / "feature-overrides.tsv"
    table.write_text(table_text, encoding="utf-8")
    return json.loads(read_feature_override_table_json(
        str(table), committed["request"], committed["resources"], str(committed["workspace"] / "out"),
    ))


def test_rows_resolve_to_the_committed_records_and_unmatched_rows_are_counted(committed):
    testa, testb, copy = committed["ids"]
    assert testa == copy
    result = _read(committed, HEADER + "".join([
        # #3 is the second copy of TESTA. TESTA_0003 spans the origin, outside
        # the crop of #1: it resolves, and Generate reports it as dormant.
        f"#3\thash={copy['misc_feature']}\toff\t\t\n",
        "#1\tlocus_tag=TESTA_0003\toff\t\t\n",
        f"#2\thash={testb['tRNA']}\t\toff\t\n",
        f"#1\thash={testa['tRNA']}\t\t\tEdited tRNA\n",
        f"#4\thash={testa['tRNA']}\toff\t\t\n",
        f"#1\thash=f00000000\toff\t\t\n",
        f"TESTC\thash={testa['tRNA']}\toff\t\t\n",
    ]))
    assert result == {
        "rows": [
            *sorted([
                {"recordKey": "linear-seq-a", "biologicalFeatureId": testa["tRNA"], "featureVisibility": None,
                 "labelVisibility": None, "labelText": "Edited tRNA"},
                {"recordKey": "linear-seq-a", "biologicalFeatureId": testa["TESTA_0003"], "featureVisibility": "off",
                 "labelVisibility": None, "labelText": None},
            ], key=lambda row: row["biologicalFeatureId"]),
            {"recordKey": "linear-seq-b", "biologicalFeatureId": testb["tRNA"], "featureVisibility": None,
             "labelVisibility": "off", "labelText": None},
            {"recordKey": "linear-seq-c", "biologicalFeatureId": copy["misc_feature"], "featureVisibility": "off",
             "labelVisibility": None, "labelText": None},
        ],
        "unmatchedRows": [6, 7, 8],
    }


@pytest.mark.parametrize(
    ("table_text", "context"),
    [
        ("record\tfeature\n#1\thash=x\n", {}),
        (HEADER + "#1\ttype=CDS\toff\t\t\n", {"row": 2}),
        (HEADER + "TESTA\thash=x\toff\t\t\n", {"row": 2}),
        (HEADER + "#1\tlocus_tag=TESTA_0001\thidden\t\t\n", {"row": 2}),
    ],
)
def test_a_malformed_table_is_a_classified_table_error(committed, table_text, context):
    with pytest.raises(ValidationError) as caught:
        _read(committed, table_text)
    error = serialize_web_error(caught.value, operation="readFeatureOverrideTable", stage="helper")
    assert (error["code"], error["operation"], error["context"]) == ("TABLE_INVALID", "readFeatureOverrideTable", context)
