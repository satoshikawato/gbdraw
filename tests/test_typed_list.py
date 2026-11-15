"""OV-2792: a list typed as one text is cut at the comma, the Arabic comma and a line break, in the CLI as in the Web."""

from __future__ import annotations

import pytest

from gbdraw.cli_utils.diagram_run import label_filtering_override
from gbdraw.core.typed_list import TYPED_LIST_SEPARATORS, split_typed_list
from gbdraw.session_request_codec import migrate_flat_config_overrides


@pytest.mark.parametrize(
    ("text", "items"),
    [
        ("a،b", ["a", "b"]),
        (" hypothetical ، unknown ", ["hypothetical", "unknown"]),
        ("a, b،c\nd\r\ne\rf", ["a", "b", "c", "d", "e", "f"]),
        (",،\n a ,, \n،b,", ["a", "b"]),
        ("hypothetical, putative", ["hypothetical", "putative"]),
        ("a;b؛c d", ["a;b؛c d"]),
        ("", []),
    ],
)
def test_a_typed_list_is_cut_at_the_comma_the_arabic_comma_and_a_line_break(text: str, items: list[str]) -> None:
    assert split_typed_list(text) == items


def test_the_separators_are_the_comma_the_arabic_comma_and_the_line_breaks() -> None:
    assert TYPED_LIST_SEPARATORS == ("\r\n", "\r", "\n", ",", "،")


@pytest.mark.parametrize(
    ("typed", "keywords"),
    [
        ("hypothetical،putative", ["hypothetical", "putative"]),
        ("hypothetical, putative", ["hypothetical", "putative"]),
        ("hypothetical protein", ["hypothetical protein"]),
    ],
)
def test_the_cli_label_blacklist_reads_the_typed_list(typed: str, keywords: list[str]) -> None:
    assert label_filtering_override({}, typed)["blacklist_keywords"] == keywords


def test_a_stored_label_blacklist_text_reads_the_typed_list() -> None:
    migrated = migrate_flat_config_overrides(
        {"label_blacklist": "hypothetical، putative"}, mode="circular", path="diagramOptions.configOverrides",
    )
    assert migrated == {"labels.filtering.blacklist_keywords": ["hypothetical", "putative"]}
