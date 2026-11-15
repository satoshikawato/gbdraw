"""The separators of a list typed as one text: the comma, the Arabic comma and a line break."""

from __future__ import annotations

import re

# The separators of a list a user types as one text (OV-2792): the CLI's ``--label_blacklist``,
# a stored ``label_blacklist`` text, and in the Web the ring names of Series names and the
# Label filtering Blacklist (``gbdraw/web/js/utils/typed-list.js``, which reads this tuple from
# the generated contracts). A comma, the Arabic comma U+060C (which fa and ar keyboards type and
# their catalogs write for "comma") or a line break separates two items; CRLF is one line break,
# so it comes before CR and LF. No typed list takes ";", so the Arabic semicolon is no separator.
# Text that is no typed list (slot specs, numbers, colors, font stacks, file formats, feature
# types) keeps its own "," rules.
TYPED_LIST_SEPARATORS: tuple[str, ...] = ("\r\n", "\r", "\n", ",", "،")

_SEPARATOR_RUN = re.compile("(?:" + "|".join(re.escape(separator) for separator in TYPED_LIST_SEPARATORS) + ")+")


def split_typed_list(text: str) -> list[str]:
    """The items of a typed list: ``text`` cut at each run of separators, each item stripped, blank items dropped."""

    return [item.strip() for item in _SEPARATOR_RUN.split(text) if item.strip()]


__all__ = ["TYPED_LIST_SEPARATORS", "split_typed_list"]
