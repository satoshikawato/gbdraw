"""A list typed as one text: its separators (the comma, the Arabic comma, a line break) and the characters trimmed."""

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

# The characters taken from each end of an item, in Python and in the Web alike. The set is the
# union of what either reader took before (Python's ``str.strip()`` and JavaScript's ``trim()``),
# so no item keeps an end character that it lost before on either side: tab, line feed, vertical
# tab, form feed and carriage return; the information separators U+001C-U+001F and the next-line
# control U+0085 (Python's); the space separators (space, no-break space, U+1680, U+2000-U+200A,
# U+202F, U+205F, U+3000); the line and paragraph separators U+2028 and U+2029; and U+FEFF
# (JavaScript's), the byte order mark Windows Notepad writes at the start of a file, so a list
# pasted from such a file begins with it.
TYPED_LIST_TRIMMED = (
    "\t\n\x0b\x0c\r\x1c\x1d\x1e\x1f \x85\xa0 "
    + "".join(chr(code) for code in range(0x2000, 0x200B))
    + "    　﻿"
)

_SEPARATOR_RUN = re.compile("(?:" + "|".join(re.escape(separator) for separator in TYPED_LIST_SEPARATORS) + ")+")


def split_typed_list(text: str) -> list[str]:
    """The items of a typed list: ``text`` cut at each run of separators, the ``TYPED_LIST_TRIMMED``
    characters taken from each item's ends, blank items dropped."""

    items = (item.strip(TYPED_LIST_TRIMMED) for item in _SEPARATOR_RUN.split(text))
    return [item for item in items if item]


__all__ = ["TYPED_LIST_SEPARATORS", "TYPED_LIST_TRIMMED", "split_typed_list"]
