"""The command-line options of a protein comparison's settings, which the Web's Source recipe writes (OV-3590).

``describe_protein_comparison_flags`` names an option for each setting and what the command writes without it;
these tests read ``gbdraw linear`` itself: each option sets its setting, and no option leaves the default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from gbdraw.linear import _LOSATP_MODES, _get_args, _linear_request, _read_linear_inputs
from gbdraw.session_request_codec import (
    _PROTEIN_SETTINGS_WIRE,
    describe_protein_comparison_flags,
    encode_canonical_request,
)

ROOT = Path(__file__).resolve().parents[1]
INPUTS = [str(ROOT / "tests" / "test_inputs" / f"{name}.gbk") for name in ("BGC0000708", "BGC0000709")]
FLAGS = describe_protein_comparison_flags()
# A value other than the default for each option: (word on the command line, the setting's value).
CHANGED: dict[str, tuple[str, object]] = {
    "--collinear_min_anchors": ("3", 3),
    "--collinear_max_unit_gap": ("2", 2),
    "--collinear_max_diagonal_drift": ("1", 1),
    "--collinear_max_conflicts_in_merge_gap": ("4", 4),
    "--collinear_merge_orientation": ("order", "order"),
    "--collinear_unit_mode": ("cds", "cds"),
    "--collinear_anchor_mode": ("all", "all"),
    "--collinear_search_scope": ("all", "all"),
    "--collinear_color_mode": ("average_identity", "average_identity"),
    "--losat_bin": ("/opt/losat/bin/losat", "/opt/losat/bin/losat"),
    "--ncbi_blast_bin": ("/opt/blast/bin/blastp", "/opt/blast/bin/blastp"),
    "--losat_threads": ("4", 4),
    "--losatp_max_hits": ("3", 3),
    "--losatp_max_target_seqs": ("50", 50),
    "--losatp_member_max_hits": ("7", 7),
    "--collinear_infer_orthogroups": ("on", True),
    "--collinear_max_paralog_links_per_orthogroup": ("3", 3),
}


def _settings(*options: str) -> Any:
    args = _get_args(["--gbk", *INPUTS, "--losat", "losatp", "--losatp_mode", "collinear", "-o", "out", *options])
    payload = encode_canonical_request(_linear_request(args, _read_linear_inputs(args), None)).payload
    (marker,) = [item for item in payload["comparisons"] if item["kind"] == "generatedProteinComparison"]
    return marker["settings"]


def _value(settings: Any, row: Any) -> object:
    held = settings[row["setting"]]
    return held if row["parameter"] is None else held["parameters"][row["parameter"]]


def test_every_setting_has_one_row_in_wire_order() -> None:
    settings = [row["setting"] for row in FLAGS["settings"]]
    assert list(dict.fromkeys(settings)) == [key for key, _owner, _name in _PROTEIN_SETTINGS_WIRE]
    defaults = _settings()
    parameters = [row["parameter"] for row in FLAGS["settings"] if row["setting"] == "collinearityParams"]
    assert parameters == list(defaults["collinearityParams"]["parameters"])
    assert sorted(FLAGS["modes"].values()) == sorted(_LOSATP_MODES)


def test_no_option_writes_the_default() -> None:
    defaults = _settings()
    assert {(row["setting"], row["parameter"]): _value(defaults, row) for row in FLAGS["settings"]} == {
        (row["setting"], row["parameter"]): row["cliDefault"] for row in FLAGS["settings"]
    }


@pytest.mark.parametrize("row", [row for row in FLAGS["settings"] if row["flag"] is not None], ids=lambda row: row["flag"])
def test_each_option_sets_its_setting(row: Any) -> None:
    word, value = CHANGED[row["flag"]]
    assert value != row["cliDefault"]
    if isinstance(value, bool):
        assert word == FLAGS["words"][str(value).lower()]
    assert _value(_settings(row["flag"], word), row) == value


def test_a_setting_without_an_option_is_the_command_default() -> None:
    assert [row["setting"] for row in FLAGS["settings"] if row["flag"] is None] == ["orthogroupMembershipMode"]
