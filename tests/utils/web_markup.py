"""Read gbdraw/web/index.html as the English markup its t('English') wrappers stand for."""

from __future__ import annotations

import re
from pathlib import Path

_KEY = r"""t\((?:'(?P<sq>(?:[^'\\]|\\.)*)'|"(?P<dq>(?:[^"\\]|\\.)*)")\)"""
_ATTR = re.compile(rf""":(?P<attr>[\w-]+)="{_KEY}\"""")
_TEXT = re.compile(rf"\{{\{{\s*{_KEY}\s*\}}\}}")


def _english(match: re.Match[str]) -> str:
    key = match.group("sq") if match.group("sq") is not None else match.group("dq")
    return key.replace("\\'", "'")


def english_markup(source: str) -> str:
    """Rewrite plain `:attr="t('X')"` to `attr="X"` and `{{ t('X') }}` to `X`.

    Keys with `{name}` params and `t()` calls with arguments are left as written;
    a test that pins one asserts its t() form explicitly.
    """

    def attr(match: re.Match[str]) -> str:
        return f'{match.group("attr")}="{_english(match)}"'

    return _TEXT.sub(_english, _ATTR.sub(attr, source))


def read_english_index(path: Path) -> str:
    return english_markup(path.read_text(encoding="utf-8"))
