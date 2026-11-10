"""A TLOSATX genetic code is one of the NCBI tables LOSAT reads, refused before any search otherwise (OV-3326)."""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw.api import LosatSearchOptions
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.comparisons.losat_runtime import LOSAT_GENETIC_CODE_TEXT, LOSAT_GENETIC_CODES
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main

REPO_ROOT = Path(__file__).resolve().parents[1]
LAMBDA = REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"
DE3 = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "de3" / "NC_042057.1.gb"


def test_the_tables_are_the_ncbi_tables_losat_reads() -> None:
    # LOSAT's TABLES (src/utils/genetic_code.rs, from NCBI gc.prt) and Biopython's CodonTable list the same ids.
    assert LOSAT_GENETIC_CODES == frozenset((*range(1, 7), *range(9, 17), *range(21, 34)))
    assert LOSAT_GENETIC_CODE_TEXT == "1\u20136, 9\u201316, 21\u201333"


@pytest.mark.parametrize("code", [7, 8, 17, 20, 34, 99])
def test_record_gencodes_refuse_a_code_losat_does_not_read(code: int) -> None:
    with pytest.raises(ValidationError) as caught:
        LosatSearchOptions(program="tlosatx", record_gencodes=(11, code))
    assert caught.value.diagnostic == {
        "code": "COMPARISON_INPUT", "reason": "GENETIC_CODE", "field": "record_gencodes",
    }
    assert LOSAT_GENETIC_CODE_TEXT in str(caught.value)


def test_record_gencodes_keep_every_table_and_none() -> None:
    search = LosatSearchOptions(program="tlosatx", record_gencodes=(*sorted(LOSAT_GENETIC_CODES), None))
    assert search.record_gencodes == (*sorted(LOSAT_GENETIC_CODES), None)


def test_ring_gencodes_refuse_a_code_losat_does_not_read() -> None:
    with pytest.raises(ValidationError) as caught:
        CircularDiagramOptions(conservation_losat_gencodes=(7,), losat_search=LosatSearchOptions(program="tlosatx"))
    assert caught.value.diagnostic == {
        "code": "COMPARISON_INPUT", "reason": "GENETIC_CODE", "field": "conservation_losat_gencodes",
    }


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
