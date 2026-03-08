from __future__ import annotations

from gbdraw.core.text import parse_mixed_content_text


def test_parse_mixed_content_text_preserves_leading_plain_text_and_order() -> None:
    parts = parse_mixed_content_text(
        "Erythromycin A biosynthetic gene cluster from <i>Saccharopolyspora erythraea</i>"
    )

    assert parts == [
        {"text": "Erythromycin A biosynthetic gene cluster from ", "italic": False},
        {"text": "Saccharopolyspora erythraea", "italic": True},
    ]


def test_parse_mixed_content_text_preserves_multiple_inline_segments() -> None:
    parts = parse_mixed_content_text("foo <i>bar</i> baz <i>qux</i>")

    assert parts == [
        {"text": "foo ", "italic": False},
        {"text": "bar", "italic": True},
        {"text": " baz ", "italic": False},
        {"text": "qux", "italic": True},
    ]
