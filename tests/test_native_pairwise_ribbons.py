"""The Rust ribbons give the bytes of the Python drawer, and decline what they do not draw."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import LinearComparison
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear import linear_main
from gbdraw.render import native
from gbdraw.render.native import native_engine

pytest.importorskip("gbdraw_core")

INPUTS = Path(__file__).parent / "test_inputs"
RECORDS = [str(INPUTS / f"{name}.gb") for name in ("AP027078", "AP027131")]
TABLE = str(INPUTS / "AP027078_AP027131.tblastx.out")


def _draw(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, native: bool, *flags: str) -> bytes:
    monkeypatch.setenv("GBDRAW_NATIVE", "1" if native else "0")
    prefix = tmp_path / ("native" if native else "python")
    linear_main(["--gbk", *RECORDS, "-b", TABLE, "-o", str(prefix), "-f", "svg", "--overwrite", *flags])
    return prefix.with_suffix(".svg").read_bytes()


@pytest.mark.parametrize(
    "flags",
    [
        (),
        ("--pairwise_match_style", "curve"),
        ("--evalue", "1e-5", "--identity", "40", "--bitscore", "80"),
        ("--identity_scale", "30", "90"),
    ],
    ids=["ribbon", "curve", "filtered", "identity-scale"],
)
def test_native_ribbons_match_the_python_document(tmp_path, monkeypatch, flags):
    assert _draw(tmp_path, monkeypatch, True, *flags) == _draw(tmp_path, monkeypatch, False, *flags)


def test_the_engine_can_be_turned_off(monkeypatch):
    monkeypatch.setenv("GBDRAW_NATIVE", "0")
    assert native_engine() is None
    monkeypatch.delenv("GBDRAW_NATIVE")
    assert native_engine() is not None


def test_the_engine_has_this_packages_contract():
    assert native_engine().CONTRACT_VERSION == native.ENGINE_CONTRACT


def test_an_engine_with_another_contract_is_left_unused(monkeypatch, caplog):
    stale = SimpleNamespace(CONTRACT_VERSION=native.ENGINE_CONTRACT + 1)
    monkeypatch.setattr(native.importlib, "import_module", lambda _name: stale)
    assert native._loaded_engine.__wrapped__() is None
    assert "Rebuild or reinstall gbdraw_core" in caplog.text


def test_a_refused_comparison_is_drawn_in_python(tmp_path, monkeypatch):
    calls = []

    class Refusing:
        @staticmethod
        def pairwise_ribbons(params, columns):
            calls.append(params)
            raise ValueError("no")

    expected = _draw(tmp_path, monkeypatch, False)
    monkeypatch.setattr("gbdraw.render.groups.linear.pairwise_match.native_engine", lambda: Refusing)
    assert _draw(tmp_path, monkeypatch, True) == expected
    assert calls


def test_a_whole_drawing_serializes_the_native_children(monkeypatch):
    records = []
    for name in ("Q1", "S1"):
        record = SeqRecord(Seq("ACGT" * 750), id=name)
        record.annotations["molecule_type"] = "DNA"
        records.append(record)
    table = pd.DataFrame(
        [
            ("Q1", "S1", 90.0, 100, 1, 0, 101, 200, 401, 500, 1e-20, 150.0),
            ("Q1", "S1", 70.0, 50, 3, 1, 301, 350, 750, 701, 1e-5, 60.0),
        ],
        columns=COMPARISON_COLUMNS,
    )
    config = apply_config_overrides(None, {"labels.linear.scope": "none", "canvas.show_gc": False, "canvas.show_skew": False})

    def text(native: str) -> str:
        monkeypatch.setenv("GBDRAW_NATIVE", native)
        return assemble_linear_diagram_from_records(
            records, cfg=config, linear_comparisons=[LinearComparison(0, 1, table)], legend="none"
        ).tostring()

    assert text("1").count("data-gbdraw-match-id") == 2
    assert text("1") == text("0")
