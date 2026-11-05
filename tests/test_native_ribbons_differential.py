"""Random comparisons give the same bytes from the Rust engine and from the Python drawer."""

from __future__ import annotations

import hashlib
import json
import random
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from gbdraw.comparisons.identity_scale import IdentityColorDomain
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.layout.linear_multi_record import LinearComparisonEndpoint
from gbdraw.render.group_content import GroupContent
from gbdraw.render.groups.linear import pairwise_match
from gbdraw.render.native import ENGINE_CONTRACT

gbdraw_core = pytest.importorskip("gbdraw_core")

ROOT = Path(__file__).resolve().parents[1]
WASM_ENGINE = ROOT / "gbdraw" / "web" / "wasm" / "gbdraw-core" / "gbdraw_core.wasm"
WASM_REPLAY = ROOT / "tests" / "web" / "helpers" / "gbdraw-core-replay.mjs"

CASES_PER_SEED = 500
# ENGINE_CONTRACT and the sha256 of what the engine draws for seed 0. A change of the drawn bytes
# raises CONTRACT_VERSION and ENGINE_CONTRACT, so an installed or cached older engine goes unused.
DRAWN_BYTES = (1, "29689f912ce58f137245938fef13be844eccd211528186a3a0edba0f2380637b")

RECORD_IDS = ("Q1", "NC_000913.3", "a&b", "<x>", 'q"uote', "tab\there", "line\nbreak", "cr\rhere", "  spaced  ", "", "日本語")
HEX_COLORS = ("#FFE7E7", "#ff7272", "#000", "#FFF", "#abc", " #123456 ", "#1a2B3c")
REFUSED_COLORS = ("#00000000", "red", "", "#12", "#GGGGGG")
ROW_IDENTITIES = (0.0, 30.0, 50.0, 90.0, 100.0, 100.00000000000001, -1.0, 150.0, 33.333333333333336)
SCORES = (
    0.0, -0.0, 1e-300, 5e-324, 1e-5, 1e-4, 0.1, 0.30000000000000004, 1e300, 1.7976931348623157e308,
    123456789012.0, 1234567890123.0, 999999999999.5, 0.000099999999999999,
    float("nan"), float("inf"), float("-inf"),
)
ENDS = ("qstart", "qend", "sstart", "send")


def _position(rng: random.Random) -> int:
    kind = rng.random()
    if kind < 0.5:
        return rng.randint(1, 20_000)
    if kind < 0.8:
        return rng.randint(1, 12_000_000)
    if kind < 0.9:
        return rng.choice((0, 1, 2**31, 2**53, 2**53 + 1))
    return rng.randint(-1000, 1000)


def _identity(rng: random.Random) -> float:
    kind = rng.random()
    if kind < 0.3:
        return float(rng.randint(0, 100))
    if kind < 0.8:
        return round(rng.uniform(0, 100), rng.randint(0, 6))
    return rng.choice(ROW_IDENTITIES)


def _score(rng: random.Random) -> float:
    if rng.random() < 0.5:
        return rng.choice(SCORES)
    return rng.uniform(0, 10) * 10.0 ** rng.randint(-300, 300)


def _frame(rng: random.Random) -> pd.DataFrame:
    """A table in one of the shapes that reach each draw order key, up to rows equal in every key."""

    rows = rng.choice((1, 1, 2, 3, 5, 8, 13, 40))
    shape = rng.choice(("mixed", "mixed", "equal spans", "equal keys"))
    ends = [_position(rng) for _ in ENDS]
    identity = _identity(rng)
    table = []
    for _ in range(rows):
        row_ends, row_identity = ends, identity
        if shape == "mixed":
            row_ends, row_identity = [_position(rng) for _ in ENDS], _identity(rng)
        elif shape == "equal spans":
            query_shift, subject_shift = rng.randint(0, 2), rng.randint(0, 50)
            row_ends = [ends[0] + query_shift, ends[1] + query_shift, ends[2] + subject_shift, ends[3] + subject_shift]
        counts = (rng.randint(0, 5000) for _ in range(3))
        table.append(("Q", "S", row_identity, *counts, *row_ends, _score(rng), _score(rng)))
    frame = pd.DataFrame(table, columns=COMPARISON_COLUMNS)
    variant = rng.random()
    if variant < 0.15:
        for column in rng.sample((*ENDS, "alignment_length", "mismatches", "gap_opens"), k=rng.randint(1, 4)):
            frame[column] = frame[column].astype(float) + (rng.choice((0.0, 0.5)) if column in ENDS else 0.0)
    elif variant < 0.25:
        frame["identity"] = frame["identity"].round().astype("int64")
    elif variant < 0.3:
        frame["extra_column"] = "x"
    elif variant < 0.33:
        frame.loc[0, "identity"] = float("nan")
    return frame


def _config(rng: random.Random) -> SimpleNamespace:
    colors = HEX_COLORS if rng.random() < 0.95 else HEX_COLORS + REFUSED_COLORS
    low, high = rng.choice(((0.0, 100.0), (30.0, 90.0), (50.0, 50.0), (90.0, 30.0), (rng.uniform(0, 100), rng.uniform(0, 100))))
    return SimpleNamespace(
        fill_color="#cccccc",
        identity_domain=IdentityColorDomain(low, high),
        min_color=rng.choice(colors),
        max_color=rng.choice(colors),
        fill_opacity=rng.choice((1.0, 0.5, 1, 0, 0.333333333, 1e-05, 0.30000000000000004)),
        stroke_color=rng.choice(("none", "#000000", "black", "#123", "a&b", "url(#x)")),
        stroke_width=rng.choice((0.0, 0, 1, 0.25, 1e-05, 2.5)),
        collinearity_orientation_colors={},
        collinearity_orientation_min_colors={},
        match_style=rng.choice(("ribbon", "curve", "Curve", "straight")),
        curve_tension=rng.choice((0.5, 0.0, 1.0, 0.25, 1 / 3, rng.uniform(-1, 2))),
    )


def _endpoint(rng: random.Random, record: int) -> LinearComparisonEndpoint:
    return LinearComparisonEndpoint(
        record_index=record,
        row=record,
        x=rng.choice((0.0, 10.5, -3.25, -0.0003, rng.uniform(-1e4, 1e4))),
        px_per_bp=rng.choice((1e-6, 0.001, 1.0, 1 / 3, rng.uniform(1e-6, 2))),
        anchor_y=rng.choice((0.0, 60.0, 120.5, -0.0002, rng.uniform(-500, 500))),
        attachment_offsets_px={},
    )


def _case(rng: random.Random) -> dict:
    return {
        "frame": _frame(rng),
        "count": rng.randint(1, 12),
        "config": _config(rng),
        "records": [SimpleNamespace(id=rng.choice(RECORD_IDS), name="n", description="d") for _ in range(2)],
        "query": _endpoint(rng, 0),
        "subject": _endpoint(rng, 1),
    }


def _children(monkeypatch: pytest.MonkeyPatch, case: dict, *, native: bool) -> tuple[str, bool]:
    """The group's children text, or the error it stopped with, and whether the engine drew it."""

    monkeypatch.setenv("GBDRAW_NATIVE", "1" if native else "0")
    try:
        group = pairwise_match.PairWiseMatchGroup(
            {},
            case["frame"],
            case["count"],
            case["config"],
            case["records"],
            feature_dom_index=None,
            record_transforms=None,
            query_placement=case["query"],
            subject_placement=case["subject"],
        ).match_group
    except Exception as error:
        return type(error).__name__, False
    return GroupContent.of(group).children_text, group.serialized_children is not None


@pytest.mark.parametrize("seed", range(4))
def test_random_comparisons_give_the_python_bytes(monkeypatch, seed):
    rng = random.Random(seed)
    drawn = 0
    for index in range(CASES_PER_SEED):
        case = _case(rng)
        expected, _ = _children(monkeypatch, case, native=False)
        text, by_engine = _children(monkeypatch, case, native=True)
        assert text == expected, f"seed {seed}, case {index}"
        drawn += by_engine
    # A refusal also gives the Python bytes, so the comparison above says nothing without this.
    assert drawn > CASES_PER_SEED * 0.8


def test_a_change_of_the_drawn_bytes_raises_the_contract(monkeypatch):
    rng = random.Random(0)
    digest = hashlib.sha256()
    for _ in range(CASES_PER_SEED):
        text, by_engine = _children(monkeypatch, _case(rng), native=True)
        if by_engine:
            digest.update(text.encode("utf-8"))
    assert (ENGINE_CONTRACT, digest.hexdigest()) == DRAWN_BYTES, (
        "The drawn bytes or the contract changed: raise CONTRACT_VERSION and ENGINE_CONTRACT, then record both here."
    )


def test_a_draw_order_key_that_is_not_a_number_is_left_to_python(monkeypatch):
    rng = random.Random(0)
    case = _case(rng)
    case["config"].min_color, case["config"].max_color = "#000000", "#ffffff"
    case["frame"] = pd.DataFrame(
        [("Q", "S", 90.0, 1, 0, 0, -1e308, 1e308, 5.0, 5.0, 0.0, 1.0), ("Q", "S", 80.0, 1, 0, 0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0)],
        columns=COMPARISON_COLUMNS,
    )
    text, by_engine = _children(monkeypatch, case, native=True)
    assert not by_engine
    assert text == _children(monkeypatch, case, native=False)[0]


def test_the_wasm_engine_gives_the_native_bytes(monkeypatch, tmp_path):
    node = shutil.which("node")
    if node is None or not WASM_ENGINE.is_file():
        pytest.skip("Needs Node and tools/build_gbdraw_core_wasm.py")
    calls = []

    class Recording:
        @staticmethod
        def pairwise_ribbons(params: str, columns: bytes) -> bytes:
            text = gbdraw_core.pairwise_ribbons(params, columns)
            calls.append({"params": params, "columns": columns.hex(), "sha256": hashlib.sha256(text).hexdigest()})
            return text

    monkeypatch.setattr(pairwise_match, "native_engine", lambda: Recording)
    rng = random.Random(1)
    for _ in range(CASES_PER_SEED):
        _children(monkeypatch, _case(rng), native=True)
    cases = tmp_path / "cases.json"
    cases.write_text(json.dumps(calls), encoding="utf-8")
    replay = subprocess.run(
        [node, str(WASM_REPLAY), str(WASM_ENGINE), str(cases), str(ENGINE_CONTRACT)],
        check=True, capture_output=True, text=True,
    )
    assert len(calls) > CASES_PER_SEED * 0.8
    assert replay.stdout.strip() == f"{len(calls)} of {len(calls)} identical"
