"""Each Linear row binding of a Session 51 names its File card (``fileUid``, OV-3703).

Save stores equal file bytes once, so two File cards whose files have the same
bytes would read as one card after Load. The Web writes each row's card; the
CLI writes one card per ``--gbk`` value; Python keeps the key when it writes a
Session again, and refuses a ``fileUid`` that is no non-empty text. The Web
side is ``tests/web/linear-file-card-identity.test.mjs``.
"""

from __future__ import annotations

import copy
import shutil
from pathlib import Path
from typing import Any

import pytest

from gbdraw.cli_utils.session import render_main
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.session_io import validate_session
from tests.utils.session_files import read_session_payload

FIXTURES = Path(__file__).resolve().parent / "fixtures"
TWO_CARDS = FIXTURES / "sessions" / "web-linear-equal-file-cards.v51.gbdraw"


def _linear_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    (drawing,) = [drawing for drawing in document["drawings"] if drawing["mode"] == "linear"]
    rows: list[dict[str, Any]] = drawing["webFiles"]["bindings"]["linearSeqs"]
    return rows


def test_the_cli_writes_one_file_card_per_gbk_file(tmp_path: Path) -> None:
    first = tmp_path / "a.gb"
    second = tmp_path / "b.gb"
    shutil.copyfile(FIXTURES / "forced_label_underlay.gb", first)
    shutil.copyfile(first, second)
    linear_main(["--gbk", str(first), str(second), "-o", str(tmp_path / "cli"), "-f", "svg", "--save_session"])
    rows = _linear_rows(read_session_payload(tmp_path / "cli.gbdraw"))
    assert [row["fileUid"] for row in rows] == [row["uid"] for row in rows] == ["cli-seq-1", "cli-seq-2"]
    assert rows[0]["gb"]["resourceId"] == rows[1]["gb"]["resourceId"], "equal bytes are stored once"


def test_a_resave_keeps_the_file_cards_of_a_web_session(tmp_path: Path) -> None:
    source = read_session_payload(TWO_CARDS)
    written = [row["fileUid"] for row in _linear_rows(source)]
    assert len(set(written)) == 2, written
    render_main(["--session", str(TWO_CARDS), "-o", str(tmp_path / "fig"), "-f", "svg", "--save_session"])
    assert [row.get("fileUid") for row in _linear_rows(read_session_payload(tmp_path / "fig.gbdraw"))] == written


@pytest.mark.parametrize("file_uid", ["", " card", 7, None])
def test_a_file_card_name_that_is_no_text_is_refused(file_uid: object) -> None:
    document = copy.deepcopy(read_session_payload(TWO_CARDS))
    validate_session(document)
    _linear_rows(document)[0]["fileUid"] = file_uid
    with pytest.raises(ValidationError, match="fileUid"):
        validate_session(document)
