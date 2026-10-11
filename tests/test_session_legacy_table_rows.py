"""Sessions 31-39 read their saved table rows as the current writer writes them (OV-148).

The Web reads the same vectors (``tests/web/file-imports.test.mjs``, OV-40).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import pytest

from gbdraw.exceptions import ParseError
from gbdraw.io.colors import load_default_colors
from gbdraw.io.table_text import legacy_table_rows, repair_legacy_table_text
from gbdraw.labels.filtering import read_filter_list_file, read_qualifier_priority_file

REPO_ROOT = Path(__file__).resolve().parents[1]
CASES = json.loads(
    (REPO_ROOT / "tests" / "fixtures" / "legacy-table-row-vectors.json").read_text(encoding="utf-8")
)["cases"]
READERS: dict[str, Callable[[str], Any]] = {
    "label-whitelist": read_filter_list_file,
    "qualifier-priority": read_qualifier_priority_file,
    "default-colors": load_default_colors,
}


def _read(table: str, path: Path) -> Any:
    try:
        return READERS[table](str(path)).to_dict("records")
    except Exception as exc:  # The strict reader's verdict on the repaired rows.
        return type(exc)


@pytest.mark.parametrize("case", CASES, ids=[case["name"] for case in CASES])
def test_a_legacy_table_reads_as_the_current_writer_writes_it(case: dict[str, Any], tmp_path: Path) -> None:
    from gbdraw.io.table_text import LEGACY_TABLE_ROWS

    assert repair_legacy_table_text(
        case["text"], LEGACY_TABLE_ROWS[case["table"]], columns=case["columns"]
    ) == (case["repairedText"], case["repairs"])

    legacy = tmp_path / "legacy.tsv"
    legacy.write_text(case["text"], encoding="utf-8")
    repaired = tmp_path / "repaired.tsv"
    repaired.write_text(case["repairedText"], encoding="utf-8")
    with legacy_table_rows(True):
        assert _read(case["table"], legacy) == _read(case["table"], repaired)
    # Outside a Session 31-39 the reader stays strict.
    assert _read(case["table"], legacy) is ParseError


def test_a_valid_row_reads_the_same(tmp_path: Path) -> None:
    for case in CASES:
        path = tmp_path / f"{case['table']}.tsv"
        path.write_text(f"{case['valid']}\n", encoding="utf-8")
        strict = _read(case["table"], path)
        with legacy_table_rows(True):
            assert _read(case["table"], path) == strict, case["name"]
