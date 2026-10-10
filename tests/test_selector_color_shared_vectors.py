"""Python side of the JS/Python shared vectors for selector case and table colors."""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from gbdraw.api.prepared import resolve_feature_inputs
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import is_user_color, load_default_colors, read_color_table
from gbdraw.labels.filtering import get_label_text, preprocess_label_filtering

FIXTURES = Path(__file__).parent / "fixtures"
SELECTOR_CASES = json.loads(
    (FIXTURES / "selector_case_equivalence_cases.json").read_text(encoding="utf-8")
)["cases"]
COLOR_DOMAIN = json.loads(
    (FIXTURES / "specific_color_domain.json").read_text(encoding="utf-8")
)
DEFAULT_COLOR_DOMAIN = json.loads(
    (FIXTURES / "default_color_domain.json").read_text(encoding="utf-8")
)


def _labelled_ids(vector: dict, qualifier: str, value: str) -> list[str]:
    filtering = preprocess_label_filtering(
        {
            "blacklist_keywords": [],
            "whitelist_df": None,
            "qualifier_priority_df": None,
            "label_override_df": pd.DataFrame(
                [[vector["record_id"], "CDS", qualifier, f"^{re.escape(value)}$", "EDITED"]],
                columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
            ),
        }
    )
    labelled = []
    for feature in vector["features"]:
        seq_feature = SeqFeature(
            FeatureLocation(feature["start"], feature["end"], strand=feature["strand"]),
            type=feature["type"],
            qualifiers=feature["qualifiers"],
        )
        if get_label_text(seq_feature, filtering, record_id=vector["record_id"]) == "EDITED":
            labelled.append(feature["svg_id"])
    return labelled


@pytest.mark.parametrize("vector", SELECTOR_CASES, ids=[case["id"] for case in SELECTOR_CASES])
def test_web_selector_matches_only_its_target_in_python(vector: dict) -> None:
    expected = vector["expected"]
    assert _labelled_ids(vector, expected["qualifier"], expected["value"]) == [vector["target"]]
    for rejected in vector["rejected"]:
        assert len(_labelled_ids(vector, rejected["qualifier"], rejected["value"])) > 1


def _write_table(tmp_path: Path, color: str) -> str:
    path = tmp_path / "specific.tsv"
    path.write_text(f"CDS\tproduct\tfirst\t#123456\tFirst\nCDS\tproduct\tx\t{color}\tcap\n", encoding="utf-8")
    return str(path)


@pytest.mark.parametrize("entry", COLOR_DOMAIN["valid"], ids=lambda entry: entry["value"])
def test_specific_color_table_accepts_the_shared_domain(tmp_path: Path, entry: dict) -> None:
    table = read_color_table(_write_table(tmp_path, entry["value"]))
    assert table is not None
    assert list(table["color"]) == ["#123456", entry["value"]]


@pytest.mark.parametrize("entry", COLOR_DOMAIN["invalid"], ids=lambda entry: entry["value"])
def test_specific_color_table_rejects_colors_outside_the_shared_domain(
    tmp_path: Path, entry: dict
) -> None:
    with pytest.raises(ValidationError, match=r"line 2") as error:
        read_color_table(_write_table(tmp_path, entry["value"]))
    assert entry["value"] in str(error.value)


@pytest.mark.parametrize("entry", DEFAULT_COLOR_DOMAIN["valid"], ids=lambda entry: entry["value"])
def test_default_colors_accept_the_shared_domain(entry: dict) -> None:
    assert is_user_color(entry["value"])


@pytest.mark.parametrize("entry", DEFAULT_COLOR_DOMAIN["invalid"], ids=lambda entry: repr(entry["value"]))
def test_default_colors_reject_colors_outside_the_shared_domain(entry: dict) -> None:
    # OV-272: currentColor and inherit leave the domain as system colors did;
    # D-38: so do svgwrite's paint references, icc-color(), and the empty value.
    assert not is_user_color(entry["value"])


def _resolve_default_colors(tmp_path: Path, color: str) -> object:
    path = tmp_path / "default_colors.tsv"
    path.write_text(f"tRNA\t#123456\nCDS\t{color}\n", encoding="utf-8")
    defaults = load_default_colors(str(path))
    resolve_feature_inputs(color_table=None, default_colors=defaults, feature_visibility_table=None)
    return defaults.set_index("feature_type").at["CDS", "color"]


def test_default_colors_file_reads_a_shared_valid_row(tmp_path: Path) -> None:
    entry = next(item for item in DEFAULT_COLOR_DOMAIN["valid"] if item["value"].startswith("rgb("))
    # OV-302: the Web import dropped this row; Python keeps it as written.
    assert _resolve_default_colors(tmp_path, entry["value"]) == entry["normalized"]


def test_default_colors_file_rejects_a_shared_invalid_row(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="Invalid color 'currentColor' for feature type 'CDS'") as raised:
        _resolve_default_colors(tmp_path, "currentColor")
    assert raised.value.diagnostic == {"code": "TABLE_INVALID", "field": "color", "reason": "COLOR"}
