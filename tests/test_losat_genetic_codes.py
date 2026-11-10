"""A TLOSATX genetic code is one of the NCBI tables LOSAT reads, refused before any search otherwise (OV-3326)."""

from __future__ import annotations

from pathlib import Path
import re
import sys

import pytest

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw.api import LosatSearchOptions
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main

REPO_ROOT = Path(__file__).resolve().parents[1]
LAMBDA = REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"
DE3 = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "de3" / "NC_042057.1.gb"
# The NCBI tables LOSAT reads (TABLES in LOSAT's src/utils/genetic_code.rs, from NCBI gc.prt), and how a refusal names
# them. The behaviour tests below use these literals, so they also run on a tree without the owner constants.
TABLES = frozenset((*range(1, 7), *range(9, 17), *range(21, 34)))
TABLE_TEXT = "1\u20136, 9\u201316, 21\u201333"


def test_the_tables_are_the_ncbi_tables_losat_reads() -> None:
    from Bio.Data.CodonTable import unambiguous_dna_by_id

    from gbdraw.comparisons.losat_runtime import LOSAT_GENETIC_CODE_TEXT, LOSAT_GENETIC_CODES

    assert LOSAT_GENETIC_CODES == TABLES
    # Biopython's CodonTable, which translates /transl_table (gbdraw/core/sequence.py), lists the same ids.
    assert LOSAT_GENETIC_CODES == frozenset(unambiguous_dna_by_id)
    assert LOSAT_GENETIC_CODE_TEXT == TABLE_TEXT


@pytest.mark.parametrize("code", [7, 8, 17, 20, 34, 99])
def test_record_gencodes_refuse_a_code_losat_does_not_read(code: int) -> None:
    with pytest.raises(ValidationError) as caught:
        LosatSearchOptions(program="tlosatx", record_gencodes=(11, code))
    assert caught.value.diagnostic == {
        "code": "COMPARISON_INPUT", "reason": "GENETIC_CODE", "field": "record_gencodes",
    }
    assert TABLE_TEXT in str(caught.value)


def test_record_gencodes_keep_every_table_and_none() -> None:
    search = LosatSearchOptions(program="tlosatx", record_gencodes=(*sorted(TABLES), None))
    assert search.record_gencodes == (*sorted(TABLES), None)


def test_ring_gencodes_refuse_a_code_losat_does_not_read() -> None:
    with pytest.raises(ValidationError) as caught:
        CircularDiagramOptions(conservation_losat_gencodes=(7,), losat_search=LosatSearchOptions(program="tlosatx"))
    assert caught.value.diagnostic == {
        "code": "COMPARISON_INPUT", "reason": "GENETIC_CODE", "field": "conservation_losat_gencodes",
    }


def test_a_ring_gencode_losatn_never_uses_is_named_as_not_applying_first(tmp_path: Path) -> None:
    # The program check comes before the table check, as in Linear: one edit fixes the one cause.
    with pytest.raises(ValidationError) as caught:
        CircularDiagramOptions(conservation_losat_gencodes=(7,), losat_search=LosatSearchOptions(program="losatn"))
    assert caught.value.diagnostic["reason"] == "LOSAT_OPTION_PROGRAM"
    from gbdraw.circular import circular_main

    losat, _log = _never_searched(tmp_path)
    with pytest.raises(ValidationError) as cli:
        circular_main([
            "--gbk", str(LAMBDA), "--conservation_genome", str(LAMBDA), "--losat", "losatn",
            "--conservation_losat_gencode", "7", "--losat_bin", losat, "-o", str(tmp_path / "ring"), "-f", "svg",
        ])
    assert str(cli.value) == "--conservation_losat_gencode does not apply to --losat losatn."


def _never_searched(tmp_path: Path) -> tuple[str, Path]:
    log = tmp_path / "searched.log"
    path = tmp_path / "bin" / "losat"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"#!{sys.executable}\nimport sys\nopen({str(log)!r}, 'a').write(' '.join(sys.argv[1:]) + '\\n')\n"
        "print('losat 0.1.0') if sys.argv[1:] in (['--version'], ['-version']) else None\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return str(path), log


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(runtime_module, "_RUNTIME_VERSIONS", {})
    monkeypatch.setattr(runtime_module, "_TASK_VALUES", {}, raising=False)


def test_cli_option_and_records_table_refuse_gencode_7_before_any_search(tmp_path: Path) -> None:
    losat, log = _never_searched(tmp_path)
    with pytest.raises(ValidationError) as option:
        linear_main([
            "--gbk", str(LAMBDA), str(DE3), "--losat", "tlosatx", "--losat_gencode", "7",
            "--losat_bin", losat, "-o", str(tmp_path / "option"), "-f", "svg",
        ])
    assert option.value.diagnostic["reason"] == "GENETIC_CODE"
    assert "--losat_gencode" in str(option.value)
    table = tmp_path / "records.tsv"
    table.write_text(f"gbk\tlosat_gencode\n{LAMBDA}\t99\n{DE3}\t11\n", encoding="utf-8")
    with pytest.raises(ValidationError) as column:
        linear_main([
            "--records_table", str(table), "--losat", "tlosatx", "--losat_bin", losat,
            "-o", str(tmp_path / "table"), "-f", "svg",
        ])
    assert column.value.diagnostic["reason"] == "GENETIC_CODE"
    assert "losat_gencode column" in str(column.value)
    # No search ran: the runtime was at most asked for its version.
    searched = log.read_text(encoding="utf-8").splitlines() if log.exists() else []
    assert all(line in {"--version", "-version"} for line in searched), searched


def test_circular_cli_names_the_ring_and_reference_flags(tmp_path: Path) -> None:
    from gbdraw.circular import circular_main

    losat, _log = _never_searched(tmp_path)
    for flags, named in (
        (["--losat_gencode", "99"], "--losat_gencode"),
        (["--conservation_losat_gencode", "7"], "--conservation_losat_gencode"),
    ):
        with pytest.raises(ValidationError) as caught:
            circular_main([
                "--gbk", str(LAMBDA), "--conservation_sequence", str(LAMBDA), "--losat", "tlosatx", *flags,
                "--losat_bin", losat, "-o", str(tmp_path / "ring"), "-f", "svg",
            ])
        assert caught.value.diagnostic["reason"] == "GENETIC_CODE"
        assert named in str(caught.value)
    # A --conservation_table losat_gencode cell is named by its table and column, not by the typed field.
    table = tmp_path / "rings.tsv"
    table.write_text(f"comparison_sequence\tlosat_gencode\n{LAMBDA}\t7\n", encoding="utf-8")
    with pytest.raises(ValidationError) as cell:
        circular_main([
            "--gbk", str(LAMBDA), "--conservation_table", str(table), "--losat", "tlosatx",
            "--losat_bin", losat, "-o", str(tmp_path / "table"), "-f", "svg",
        ])
    assert cell.value.diagnostic["reason"] == "GENETIC_CODE"
    assert f"{table}: losat_gencode column" in str(cell.value)
    assert "conservation_losat_gencodes" not in str(cell.value)
    searched = _log.read_text(encoding="utf-8").splitlines() if _log.exists() else []
    assert all(line in {"--version", "-version"} for line in searched), searched


@pytest.mark.parametrize(("main", "flags"), [
    ("gbdraw.linear:linear_main", ["--losat_gencode"]),
    ("gbdraw.circular:circular_main", ["--losat_gencode", "--conservation_losat_gencode"]),
])
def test_cli_help_lists_the_tables(main: str, flags: list[str], capsys: pytest.CaptureFixture[str]) -> None:
    import importlib

    module, function = main.split(":")
    with pytest.raises(SystemExit):
        getattr(importlib.import_module(module), function)(["-h"])
    # Each option's entry: its line in the options list and the indented lines under it.
    entries = re.split(r"\n(?=  -)", capsys.readouterr().out)
    for flag in flags:
        (entry,) = [item for item in entries if item.lstrip().startswith(f"{flag} ")]
        assert TABLE_TEXT in " ".join(entry.split()), (flag, entry)
