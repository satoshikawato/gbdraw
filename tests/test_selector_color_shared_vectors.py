"""Python side of the JS/Python shared vectors for selector case and table colors."""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import read_color_table
from gbdraw.labels.filtering import get_label_text, preprocess_label_filtering

FIXTURES = Path(__file__).parent / "fixtures"
SELECTOR_CASES = json.loads(
    (FIXTURES / "selector_case_equivalence_cases.json").read_text(encoding="utf-8")
)["cases"]
COLOR_DOMAIN = json.loads(
    (FIXTURES / "specific_color_domain.json").read_text(encoding="utf-8")
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
