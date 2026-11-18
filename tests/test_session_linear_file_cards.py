"""Each Linear row binding of a Session 51 names its File card (``fileUid``, OV-3703).

Save stores equal file bytes once, so two File cards whose files have the same
bytes would read as one card after Load. The Web writes each row's card; the
CLI writes one card per ``--gbk`` value; Python keeps the key when it writes a
Session again, and refuses a ``fileUid`` that is no non-empty string (the
values in ``tests/fixtures/linear-file-uid-vectors.json``, which the Web test reads too). The Web
side is ``tests/web/linear-file-card-identity.test.mjs``.
"""

from __future__ import annotations

import copy
import json
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


@pytest.mark.parametrize("names", [("a.gb", "b.gb"), ("a/x.gb", "b/x.gb")], ids=["two names", "one name"])
def test_the_cli_writes_one_file_card_per_gbk_file(tmp_path: Path, names: tuple[str, str]) -> None:
    """Each row binding names the resource its own file's records read, also for two files of equal bytes
    (OV-3708): Load matches a File card to the records of the resource its file reads."""

    paths = [tmp_path / name for name in names]
    for path in paths:
        path.parent.mkdir(exist_ok=True)
        shutil.copyfile(FIXTURES / "forced_label_underlay.gb", path)
    linear_main(["--gbk", *map(str, paths), "-o", str(tmp_path / "cli"), "-f", "svg",
                 "--session_output", str(tmp_path / "cli.gbdraw")])
    document = read_session_payload(tmp_path / "cli.gbdraw")
    rows = _linear_rows(document)
    assert [row["fileUid"] for row in rows] == [row["uid"] for row in rows] == ["cli-seq-1", "cli-seq-2"]
    (drawing,) = [drawing for drawing in document["drawings"] if drawing["mode"] == "linear"]
    read = [record["source"]["resourceId"] for record in drawing["renderRequest"]["records"]]
    assert [row["gb"]["resourceId"] for row in rows] == read == ["record-1-genbank", "record-2-genbank"]


def test_a_resave_keeps_the_file_cards_of_a_web_session(tmp_path: Path) -> None:
    source = read_session_payload(TWO_CARDS)
    written = [row["fileUid"] for row in _linear_rows(source)]
    assert len(set(written)) == 2, written
    render_main(["--session", str(TWO_CARDS), "-o", str(tmp_path / "fig"), "-f", "svg", "--save_session"])
    assert [row.get("fileUid") for row in _linear_rows(read_session_payload(tmp_path / "fig.gbdraw"))] == written


# The values the Web reader accepts and refuses too (tests/web/linear-file-card-identity.test.mjs).
FILE_UID_VECTORS = json.loads((FIXTURES / "linear-file-uid-vectors.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("vector", FILE_UID_VECTORS["refused"], ids=lambda vector: vector["case"])
def test_a_file_card_name_that_is_no_non_empty_string_is_refused(vector: dict[str, object]) -> None:
    document = copy.deepcopy(read_session_payload(TWO_CARDS))
    validate_session(document)
    _linear_rows(document)[0]["fileUid"] = vector["fileUid"]
    with pytest.raises(ValidationError, match="fileUid"):
        validate_session(document)


@pytest.mark.parametrize("vector", FILE_UID_VECTORS["accepted"], ids=lambda vector: vector["case"])
def test_a_file_card_name_that_is_a_non_empty_string_is_accepted(vector: dict[str, object]) -> None:
    document = copy.deepcopy(read_session_payload(TWO_CARDS))
    _linear_rows(document)[0]["fileUid"] = vector["fileUid"]
    validate_session(document)
