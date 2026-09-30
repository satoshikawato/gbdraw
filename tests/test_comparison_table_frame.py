"""Comparison-table coordinate frame (Web GUI audit 2026-09-30, D-18, PD-OI-073).

A comparison table is read in the search frame: the selected and cropped
sequence, 1-based, in the source strand. The planner projects a record's
reverse complement, and a row outside the cropped record is rejected or
dropped instead of drawn beyond the record. These strict xfail tests record
the current display-frame reading (N-07, N-08); the Q-FRAME fix removes the
marks.
"""
from __future__ import annotations

import random
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw import linear as linear_cli
from gbdraw.exceptions import GbdrawError

_RNG = random.Random(7)
_X, _Y, _Z = ("".join(_RNG.choice("ACGT") for _ in range(size)) for size in (3000, 1000, 2000))
_RIBBON = re.compile(r'<path [^>]*data-gbdraw-pairwise-match-id="[^"]*"[^>]*>')


def _write_record(directory: Path, record_id: str, sequence: str) -> Path:
    record = SeqRecord(Seq(sequence), id=record_id, name=record_id, description=f"{record_id} synthetic")
    record.annotations.update(molecule_type="DNA", topology="linear")
    record.features = [SeqFeature(SimpleLocation(0, len(sequence), strand=1), type="source")]
    path = directory / f"{record_id}.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _ribbon_spans(svg_text: str) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Return the query and subject x spans of each drawn comparison ribbon."""
    spans = []
    for element in _RIBBON.findall(svg_text):
        path = re.search(r'\sd="([^"]*)"', element).group(1)
        point = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", path)]
        spans.append((
            (round(min(point[0], point[2]), 1), round(max(point[0], point[2]), 1)),
            (round(min(point[4], point[6]), 1), round(max(point[4], point[6]), 1)),
        ))
    return spans


def _render(directory: Path, prefix: str, records: list[Path], row: str, *extra: str) -> list:
    table = directory / f"{prefix}.tsv"
    table.write_text(row + "\n", encoding="utf-8")
    linear_cli.linear_main([
        "--gbk", *map(str, records), "-b", str(table), "-f", "svg", "-o", str(directory / prefix), *extra,
    ])
    return _ribbon_spans((directory / f"{prefix}.svg").read_text(encoding="utf-8"))


@pytest.mark.linear
@pytest.mark.xfail(strict=True, reason="N-08")
def test_reverse_complement_flag_matches_a_physically_reversed_record(tmp_path: Path) -> None:
    # R2 2001..3000 and R3 1..1000 are one shared block (search frame).
    r2 = _write_record(tmp_path, "R2", _X[:2000] + _Y)
    r3 = _write_record(tmp_path, "R3", _Y + _Z)
    r3_reversed = _write_record(tmp_path, "R3rc", str(Seq(_Y + _Z).reverse_complement()))
    flagged = _render(
        tmp_path, "flag", [r2, r3], "R2\tR3\t100\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847",
        "--reverse_complement", "0", "--reverse_complement", "1",
    )
    physical = _render(
        tmp_path, "physical", [r2, r3_reversed], "R2\tR3rc\t100\t1000\t0\t0\t2001\t3000\t3000\t2001\t0.0\t1847",
    )
    assert len(physical) == 1
    assert flagged == physical


@pytest.mark.linear
@pytest.mark.xfail(strict=True, reason="N-07")
def test_row_outside_a_cropped_record_is_not_drawn_beyond_the_record(tmp_path: Path) -> None:
    sa = _write_record(tmp_path, "SA", _X)
    sb = _write_record(tmp_path, "SB", _X)
    # SB is cropped to 1000 bp, so subject 1101..1300 is outside the record.
    try:
        drawn = _render(
            tmp_path, "crop", [sa, sb], "SA\tSB\t95\t200\t10\t0\t101\t300\t1101\t1300\t1e-50\t300",
            "--region", "SB:1001-2000",
        )
    except (GbdrawError, SystemExit):
        return
    assert drawn == []
