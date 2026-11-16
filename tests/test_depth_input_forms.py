"""The forms of a depth file: two columns, samtools depth and bedGraph read into one parsed form.

A file's form comes from the column count of its first data line, never from its name: the
Web reads every input from a staged copy named ``NNNN.bin``.
"""

from __future__ import annotations

import gzip
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis import depth as depth_module
from gbdraw.analysis.depth import DepthCoverage, depth_df, parse_depth_tsv, read_depth_tsv
from gbdraw.core.record_metadata import _write_coord_map
from gbdraw.exceptions import GbdrawError
from tests.utils.cli_in_process import run_gbdraw


def _write(path: Path, text: str) -> Path:
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode("utf-8")))
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _assert_same_coverage(actual: DepthCoverage, expected: DepthCoverage) -> None:
    assert list(actual.depths) == list(expected.depths)
    assert dict(actual.extents) == dict(expected.extents)
    for name, depths in expected.depths.items():
        assert actual.depths[name].dtype == depths.dtype, name
        assert np.array_equal(actual.depths[name], depths), name
        assert not actual.depths[name].flags.writeable, name


def _windows(record: SeqRecord, coverage: DepthCoverage) -> list[pd.DataFrame]:
    return [
        depth_df(record, coverage, window, step, circular=circular)
        for window, step, circular in ((1, 1, False), (25, 10, True), (100, 37, False))
    ]


def _assert_same_windows(record: SeqRecord, actual: DepthCoverage, expected: DepthCoverage) -> None:
    for got, want in zip(_windows(record, actual), _windows(record, expected), strict=True):
        pd.testing.assert_frame_equal(got, want)


def test_a_two_column_file_draws_as_the_three_column_file_of_the_record(tmp_path: Path) -> None:
    rows = [(position, (position * 7) % 23) for position in range(1, 401) if position % 9]
    three = _write(tmp_path / "three.tsv", "".join(f"rec1\t{p}\t{v}\n" for p, v in rows))
    two = _write(tmp_path / "two.tsv", "".join(f"{p}\t{v}\n" for p, v in rows))

    parsed = parse_depth_tsv(str(two))

    assert list(parsed.depths) == [""]
    assert parsed.extents[""] == parse_depth_tsv(str(three)).extents["rec1"]
    record = SeqRecord(Seq("A" * 420), id="rec1")
    _assert_same_windows(record, parsed, parse_depth_tsv(str(three)))
    # A record of any name reads it, as it reads the only reference of a samtools depth file.
    other = SeqRecord(Seq("A" * 420), id="plasmid")
    _write_coord_map(other, base=400, step=-1)
    _assert_same_windows(other, parsed, parse_depth_tsv(str(three)))


def _bedgraph_rows(rng: np.random.Generator) -> list[tuple[str, int, int, str]]:
    """bedGraph rows: zero runs as genomecov -bga writes them, gaps, fractions; references sometimes interleaved."""

    references: list[list[tuple[str, int, int, str]]] = []
    for name in rng.choice(["chr", "plasmid", "2"], size=int(rng.integers(1, 4)), replace=False):
        rows: list[tuple[str, int, int, str]] = []
        start = int(rng.integers(0, 30))
        for _ in range(int(rng.integers(1, 40))):
            length = int(rng.integers(1, 25))
            kind = rng.random()
            value = "0" if kind < 0.2 else repr(float(rng.random() * 40)) if kind < 0.4 else str(int(rng.integers(0, 300)))
            rows.append((str(name), start, start + length, value))
            start += length + (int(rng.integers(1, 10)) if rng.random() < 0.3 else 0)
        references.append(rows)
    if rng.random() < 0.7:
        return [row for rows in references for row in rows]
    # Interleaved: each reference keeps its order.
    merged: list[tuple[str, int, int, str]] = []
    while any(references):
        merged.append(references[int(rng.choice([i for i, rows in enumerate(references) if rows]))].pop(0))
    return merged


@pytest.mark.parametrize("seed", range(40))
def test_a_bedgraph_reads_as_its_per_base_rows(seed: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Each row gives its value to bases ``start + 1 .. end``, however the rows are read and spread in batches."""

    rng = np.random.default_rng(seed)
    monkeypatch.setattr(depth_module, "_CHUNK_ROWS", int(rng.choice([2, 7, 1 << 20])))
    monkeypatch.setattr(depth_module, "_SPREAD_BASES", int(rng.choice([1, 3, 64, 1 << 18])))
    rows = _bedgraph_rows(rng)
    lead = "track type=bedGraph name=x\n" if rng.random() < 0.5 else ""
    suffix = ".bedgraph.gz" if rng.random() < 0.3 else ".bedgraph"
    bedgraph = _write(tmp_path / f"cov{suffix}", lead + "".join(f"{n}\t{s}\t{e}\t{v}\n" for n, s, e, v in rows))
    per_base = _write(
        tmp_path / "cov.tsv",
        "".join(f"{n}\t{position}\t{v}\n" for n, s, e, v in rows for position in range(s + 1, e + 1)),
    )

    parsed = parse_depth_tsv(str(bedgraph))

    _assert_same_coverage(parsed, parse_depth_tsv(str(per_base)))
    pd.testing.assert_frame_equal(read_depth_tsv(str(bedgraph)), read_depth_tsv(str(per_base)))
    record = SeqRecord(Seq("A" * 300), id=rows[0][0])
    _assert_same_windows(record, parsed, parse_depth_tsv(str(per_base)))


def test_a_bedgraph_of_zero_runs_covers_the_whole_reference(tmp_path: Path) -> None:
    path = _write(tmp_path / "cov.bedgraph", "chr\t0\t4\t0\nchr\t4\t6\t3\nchr\t6\t10\t0\n")

    parsed = parse_depth_tsv(str(path))

    assert parsed.depths["chr"].tolist() == [0, 0, 0, 0, 3, 3, 0, 0, 0, 0]
    assert parsed.depths["chr"].dtype == np.uint8
    assert parsed.extents["chr"] == (1, 10)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("12\t4\n13\t5\n", {"": [0] * 11 + [4, 5]}),
        ("position\tvalue\n12\t4\n13\t5\n", {"": [0] * 11 + [4, 5]}),
        ("# made by a script\n12\t4\n13\t5\n", {"": [0] * 11 + [4, 5]}),
        ("chr\t2\t4\n", {"chr": [0, 4]}),
        ("#CHROM\tPOS\tsample.bam\nchr\t2\t4\n", {"chr": [0, 4]}),
        ("track name=depth\nchr\t2\t4\n", {"chr": [0, 4]}),
        ("chr\t1\t3\t4\n", {"chr": [0, 4, 4]}),
        ("browser position chr:1-3\ntrack type=bedGraph\n# note\nchr\t1\t3\t4\n", {"chr": [0, 4, 4]}),
        ("chrom\tstart\tend\tvalue\nchr\t1\t3\t4\n", {"chr": [0, 4, 4]}),
        ("chr\t1\t3\t4\textra\tcolumns\n", {"chr": [0, 4, 4]}),
        # A reference named "track" or "#1" is a row.
        ("track\t2\t4\n", {"track": [0, 4]}),
        ("#1\t1\t5\n#1\t2\t6\n", {"#1": [5, 6]}),
        # The header of samtools depth -H marks samtools rows of any number of samples: the first is read.
        ("#CHROM\tPOS\ta.bam\tb.bam\nchr\t1\t6001\t5001\nchr\t2\t6002\t5002\n", {"chr": [6001, 6002]}),
        # A first line whose value is empty is a header line, as on dev.
        ("chr\t1\t\nchr\t2\t5\n", {"chr": [0, 5]}),
        ("chr\t10\t20\t\nchr\t20\t30\t5\n", {"chr": [0] * 20 + [5] * 10}),
    ],
)
def test_comment_track_and_header_lines_before_the_rows_are_skipped(
    text: str, expected: dict[str, list[int]], tmp_path: Path
) -> None:
    parsed = parse_depth_tsv(str(_write(tmp_path / "depth.txt", text)))

    assert {name: depths.tolist() for name, depths in parsed.depths.items()} == expected


@pytest.mark.parametrize(
    ("text", "diagnostic"),
    [
        ("chr\n", {"code": "TABLE_INVALID", "reason": "DEPTH_COLUMNS"}),
        ("# only a comment\nchr\n", {"code": "TABLE_INVALID", "reason": "DEPTH_COLUMNS"}),
        ("chr\t-1\t5\t3\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        ("chr\t5\t5\t3\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        ("chr\t5\t4\t3\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        ("chr\t0\t5\t3\nchr\tx\t5\t3\n", {"code": "DEPTH_INVALID", "reason": "DEPTH_VALUES"}),
        ("chr\t1\t5\t-3\n", {"code": "DEPTH_INVALID", "reason": "NONNEGATIVE"}),
        ("\t1\t5\t3\n", {"code": "DEPTH_INVALID", "reason": "REFERENCE_REQUIRED"}),
        ("0\t3\n", {"code": "DEPTH_INVALID", "reason": "POSITIVE_INTEGER"}),
        # bedGraph rows of a reference come in increasing order without overlapping.
        ("chr\t0\t5\t3\nchr\t4\t8\t2\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        ("chr\t5\t8\t3\nchr\t0\t5\t2\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        ("a\t0\t5\t3\nb\t0\t5\t3\na\t2\t9\t2\n", {"code": "DEPTH_INVALID", "reason": "INTERVAL"}),
        # samtools depth of two samples without its -H header: reference, position, then one column per sample.
        # Low depths end before they start; high depths overlap the next row.
        (
            "".join(f"chr\t{p}\t{30 + p % 3}\t{28 + p % 5}\n" for p in range(1, 200)),
            {"code": "DEPTH_INVALID", "reason": "INTERVAL"},
        ),
        (
            "".join(f"chr\t{p}\t{6000 + p % 7}\t{5000 + p % 5}\n" for p in range(1, 5001)),
            {"code": "DEPTH_INVALID", "reason": "INTERVAL"},
        ),
        # A tab at the end of a line starts an empty column.
        ("chr\t2\t4\t\nchr\t3\t5\t\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        # OV-3740: a value that is not a finite number stopped the track silently.
        ("chr\t1\t5\nchr\t2\tnan\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        ("chr\t1\t5\nchr\t2\tinf\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        ("chr\t1\t5\nchr\t2\t\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        ("1\tnan\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        ("chr\t0\t5\tinf\n", {"code": "DEPTH_INVALID", "reason": "FINITE"}),
        ("chr\t0\t5\n", {"code": "DEPTH_INVALID", "reason": "POSITIVE_INTEGER"}),
    ],
)
def test_a_file_that_does_not_read_is_refused_with_its_reason(
    text: str, diagnostic: dict[str, str], tmp_path: Path
) -> None:
    path = _write(tmp_path / "depth.tsv", text)

    with pytest.raises(GbdrawError) as caught:
        parse_depth_tsv(str(path))

    assert caught.value.diagnostic == diagnostic


@pytest.mark.parametrize("text", ["chr\t0\t5\t3\nchr\t5\t9\t2\nchr\t3\t4\t1\n", "chr\t0\t5\t3\nchr\t0\t5\t3\n"])
def test_bedgraph_rows_that_go_back_are_refused_across_read_chunks(
    text: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(depth_module, "_CHUNK_ROWS", 1)

    with pytest.raises(GbdrawError) as caught:
        parse_depth_tsv(str(_write(tmp_path / "depth.bedgraph", text)))

    assert caught.value.diagnostic == {"code": "DEPTH_INVALID", "reason": "INTERVAL"}


def test_a_nan_value_no_longer_ends_the_track(tmp_path: Path) -> None:
    """OV-3740: on dev a ``nan`` row made every window from there on NaN, and the track stopped drawing."""

    path = _write(tmp_path / "depth.tsv", "".join(f"chr\t{p}\t{'nan' if p == 50 else 3}\n" for p in range(1, 101)))

    with pytest.raises(GbdrawError, match="finite"):
        depth_df(SeqRecord(Seq("A" * 100), id="chr"), parse_depth_tsv(str(path)), 10, 10, circular=False)


def test_the_form_does_not_depend_on_the_file_name(tmp_path: Path) -> None:
    """The Web reads an input from a staged copy named by its order (``resources/0001.bin``)."""

    text = "track type=bedGraph\nchr\t0\t3\t2\nchr\t3\t9\t5\n"

    named = parse_depth_tsv(str(_write(tmp_path / "coverage.bedgraph", text)))
    staged = parse_depth_tsv(str(_write(tmp_path / "0001.bin", text)))

    _assert_same_coverage(staged, named)


def test_parsing_a_bedgraph_holds_its_dense_values_and_one_batch(tmp_path: Path) -> None:
    """A run of 16 Mbp costs its uint8 values (16 MB, one growth step) and one batch of bases, not 8 bytes per base."""

    length = 16_000_000
    path = _write(tmp_path / "long.bedgraph", f"chr\t0\t{length // 2}\t3\nchr\t{length // 2}\t{length}\t200\n")

    tracemalloc.start()
    try:
        parsed = parse_depth_tsv(str(path))
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert parsed.depths["chr"].dtype == np.uint8
    assert parsed.nbytes == length
    # The stored values, the last growth step of their array, and the transients of one batch of 256 Ki bases
    # (36 MB measured, 49 MB allowed), where 8 bytes per base would take 128 MB.
    assert peak < length + length // 2 + 24 * 2**20


def test_the_cli_draws_a_bedgraph_as_the_per_base_file_of_the_same_values(tmp_path: Path) -> None:
    rng = np.random.default_rng(7)
    length = 6000
    record = SeqRecord(Seq("".join(rng.choice(list("ACGT"), size=length))), id="rec1", name="rec1")
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    genbank = tmp_path / "rec1.gb"
    SeqIO.write(record, genbank, "genbank")
    values = np.repeat(rng.integers(0, 60, size=60), 100)
    runs: list[tuple[int, int, int]] = []
    for position, value in enumerate(values.tolist()):
        if runs and runs[-1][2] == value:
            runs[-1] = (runs[-1][0], position + 1, value)
        else:
            runs.append((position, position + 1, value))
    bedgraph = _write(tmp_path / "rec1.bedgraph", "".join(f"rec1\t{s}\t{e}\t{v}\n" for s, e, v in runs))
    per_base = _write(tmp_path / "rec1.depth.tsv", "".join(f"rec1\t{p + 1}\t{v}\n" for p, v in enumerate(values.tolist())))

    outputs = []
    for name, depth in (("bedgraph", bedgraph), ("per_base", per_base)):
        result = run_gbdraw(
            ["circular", "--gbk", str(genbank), "--depth_track", str(depth), "-o", str(tmp_path / name), "-f", "svg"],
            cwd=tmp_path,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append((tmp_path / f"{name}.svg").read_bytes())

    assert outputs[0] == outputs[1]
