"""The parsed depth form: one dense array per reference, binned exactly as the row table was.

``_reference_*`` below is the row-table algorithm that ``depth_df`` used before
the parsed form became compact. It stays here only, as the oracle the
randomized cases compare against value for value.
"""

from __future__ import annotations

import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis import depth as depth_module
from gbdraw.analysis.depth import DepthCoverage, depth_df, parse_depth_tsv, read_depth_tsv
from gbdraw.core.record_metadata import _read_coord_map, _source_to_local_index, _write_coord_map
from gbdraw.exceptions import GbdrawError, ParseError, ValidationError
from gbdraw.io.table_text import EXACT_FLOAT_READ_OPTIONS

_COLUMNS = ["reference_name", "position", "depth"]


def _invalid(reason: str) -> dict[str, str]:
    return {"code": "DEPTH_INVALID", "reason": reason}


def _reference_validate(table: pd.DataFrame) -> None:
    if table.empty:
        return
    if table["reference_name"].isna().any():
        raise ValidationError("Depth table contains missing reference names.", diagnostic=_invalid("REFERENCE_REQUIRED"))
    if (table["position"] <= 0).any():
        raise ValidationError("Depth positions must be 1-based positive integers.", diagnostic=_invalid("POSITIVE_INTEGER"))
    if (table["depth"] < 0).any():
        raise ValidationError("Depth values must be non-negative.", diagnostic=_invalid("NONNEGATIVE"))


def _reference_parse(path: Path) -> pd.DataFrame:
    skiprows = 1 if depth_module._first_line_has_header(path) else 0
    try:
        table = pd.read_csv(
            path, sep="\t", header=None, names=_COLUMNS, usecols=[0, 1, 2], skiprows=skiprows,
            dtype={"reference_name": "category", "position": "int64", "depth": "float64"},
            **EXACT_FLOAT_READ_OPTIONS,
        )
    except ValueError as exc:
        raise ParseError(f"Unable to parse depth file '{path}': {exc}", diagnostic=_invalid("DEPTH_VALUES")) from exc
    _reference_validate(table)
    return table


def _linear_windows(record_len, window, step):
    """The samples of a linear record and the bounds of their windows: centred, kept inside the record."""
    samples = np.arange(0, record_len, int(step), dtype=np.int64)
    starts = samples - int(window) // 2
    return samples, np.clip(starts, 0, record_len), np.clip(starts + int(window), 0, record_len)


def _reference_depth_df(record, table, window, step, *, min_depth=None, max_depth=None):
    table = table.loc[:, _COLUMNS].copy()
    table["reference_name"] = table["reference_name"].astype("string")
    _reference_validate(table)
    record_len = len(record.seq)
    if record_len <= 0:
        return pd.DataFrame(columns=["position", "depth"])
    record_id = str(record.id)
    rows = table[table["reference_name"] == record_id]
    if rows.empty:
        references = [str(value) for value in table["reference_name"].dropna().unique()]
        if len(references) != 1:
            raise ValidationError(
                f"Depth file references do not match record '{record_id}'. Available references: {', '.join(references)}",
                diagnostic=_invalid("REFERENCE_MISMATCH"),
            )
        rows = table[table["reference_name"] == references[0]]
    base, coord_step = _read_coord_map(record)
    local = _source_to_local_index(rows["position"].to_numpy(dtype=np.int64), base, coord_step)
    inside = (local >= 0) & (local < record_len)
    local = local[inside]
    coverage = np.zeros(record_len, dtype=np.float64)
    if local.size:
        summed = np.bincount(local, weights=rows["depth"].to_numpy(dtype=np.float64)[inside], minlength=record_len)
        counts = np.bincount(local, minlength=record_len)
        nonzero = counts > 0
        coverage[nonzero] = summed[nonzero] / counts[nonzero]
    prefix = np.concatenate(([0.0], np.cumsum(coverage)))
    samples, starts, ends = _linear_windows(record_len, window, step)
    means = (prefix[ends] - prefix[starts]) / (ends - starts)
    if min_depth is not None or max_depth is not None:
        means = np.clip(means, min_depth, max_depth)
    return pd.DataFrame({"position": samples, "depth": means})


def _outcome(build):
    try:
        return build()
    except GbdrawError as exc:
        return (type(exc).__name__, str(exc), exc.diagnostic)


def _assert_same(actual, expected) -> None:
    if isinstance(expected, tuple):
        assert actual == expected
        return
    assert isinstance(actual, pd.DataFrame), actual
    assert list(actual.columns) == list(expected.columns)
    for column in expected.columns:
        assert actual[column].dtype == expected[column].dtype, column
        assert np.array_equal(actual[column].to_numpy(), expected[column].to_numpy(), equal_nan=True), column


def _depth_text(rng: np.random.Generator) -> str:
    """Rows with gaps, duplicates, unsorted runs, several references and one value kind per file."""

    names = list(rng.choice(["chr", "plasmid", "rec1", "1"], size=int(rng.integers(1, 4)), replace=False))
    kind = rng.choice(["uint8", "uint16", "uint32", "huge", "fraction", "mixed"])
    rows: list[str] = []
    for name in names:
        length = int(rng.integers(1, 300))
        positions = np.flatnonzero(rng.random(length) < rng.uniform(0.3, 1.0)) + 1
        if positions.size and rng.random() < 0.4:
            positions = np.concatenate([positions, rng.choice(positions, size=int(rng.integers(1, 20)))])
        if rng.random() < 0.4:
            rng.shuffle(positions)
        for position in positions:
            if kind == "uint8":
                value = str(int(rng.integers(0, 256)))
            elif kind == "uint16":
                value = str(int(rng.integers(0, 65536)))
            elif kind == "uint32":
                value = str(int(rng.integers(0, 2**32)))
            elif kind == "huge":
                value = str(int(rng.integers(2**32, 2**40)))
            elif kind == "fraction":
                value = repr(float(rng.random() * 50))
            else:
                value = rng.choice([str(int(rng.integers(0, 40))), "0.5", "-0", "1e-300", "3.25"])
            rows.append(f"{name}\t{position}\t{value}")
    if len(names) > 1 and rng.random() < 0.5:
        rng.shuffle(rows)
    header = "reference\tposition\tdepth\n" if rng.random() < 0.2 else ""
    return header + "\n".join(rows) + "\n"


def _record(rng: np.random.Generator, text: str) -> SeqRecord:
    names = sorted({line.split("\t")[0] for line in text.splitlines() if line and not line.startswith("reference\t")})
    record_id = str(rng.choice(names + ["absent"]))
    length = int(rng.integers(0, 320)) if rng.random() < 0.05 else int(rng.integers(1, 320))
    record = SeqRecord(Seq("A" * length), id=record_id)
    if rng.random() < 0.6:
        _write_coord_map(record, base=int(rng.integers(-20, 340)), step=int(rng.choice([1, -1])))
    return record


def _binning(rng: np.random.Generator) -> dict[str, object]:
    min_depth = float(rng.choice([0.0, 0.5, 3.0])) if rng.random() < 0.3 else None
    max_depth = float(rng.choice([4.0, 30.0, 70000.0])) if rng.random() < 0.3 else None
    if min_depth is not None and max_depth is not None and min_depth > max_depth:
        min_depth, max_depth = max_depth, min_depth
    return {
        "window": int(rng.integers(1, 60)),
        "step": int(rng.integers(1, 60)),
        "min_depth": min_depth,
        "max_depth": max_depth,
    }


@pytest.mark.parametrize("seed", range(10))
def test_binned_depth_equals_the_row_table_algorithm(seed: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Parsed files and API tables bin exactly as the row table did, also across read and sum chunk boundaries."""

    rng = np.random.default_rng(seed)
    for case in range(30):
        monkeypatch.setattr(depth_module, "_CHUNK_ROWS", int(rng.choice([11, 97, 1 << 20])))
        monkeypatch.setattr(depth_module, "_SUM_CHUNK_BASES", (1, 5, 1 << 20)[case % 3])
        path = tmp_path / f"depth-{case}.tsv"
        text = _depth_text(rng)
        path.write_text(text, encoding="utf-8")
        record = _record(rng, text)
        binning = _binning(rng)
        args = (binning["window"], binning["step"])
        kwargs = {key: binning[key] for key in ("min_depth", "max_depth")}
        expected = _outcome(lambda: _reference_depth_df(record, _reference_parse(path), *args, **kwargs))
        _assert_same(_outcome(lambda: depth_df(record, parse_depth_tsv(str(path)), *args, circular=False, **kwargs)), expected)
        _assert_same(_outcome(lambda: depth_df(record, read_depth_tsv(str(path)), *args, circular=False, **kwargs)), expected)


def test_a_parsed_file_is_one_read_only_array_per_reference_in_the_smallest_exact_dtype(tmp_path: Path) -> None:
    path = tmp_path / "depth.tsv"
    path.write_text(
        "b\t3\t7\nb\t1\t300\na\t2\t1\nc\t1\t70000\nd\t2\t0.5\ne\t1\t5000000000\nf\t1\t2\nf\t1\t3\n",
        encoding="utf-8",
    )

    parsed = parse_depth_tsv(str(path))

    assert isinstance(parsed, DepthCoverage)
    assert list(parsed.depths) == ["b", "a", "c", "d", "e", "f"]
    expected = {
        "b": (np.uint16, [300, 0, 7]),
        "a": (np.uint8, [0, 1]),
        "c": (np.uint32, [70000]),
        "d": (np.float64, [0.0, 0.5]),
        "e": (np.float64, [5000000000.0]),
        "f": (np.float64, [2.5]),
    }
    for name, (dtype, values) in expected.items():
        depths = parsed.depths[name]
        assert depths.dtype == dtype, name
        assert depths.tolist() == values, name
        assert not depths.flags.writeable, name
    assert parsed.nbytes == sum(array.nbytes for array in parsed.depths.values())


def test_a_whole_genome_table_takes_about_two_bytes_per_base(tmp_path: Path) -> None:
    path = tmp_path / "depth.tsv"
    path.write_text("".join(f"chr\t{position}\t{position % 400}\n" for position in range(1, 20001)), encoding="utf-8")

    assert parse_depth_tsv(str(path)).nbytes == 2 * 20000


@pytest.mark.parametrize(
    "text",
    [
        "",
        "reference\tposition\tdepth\n",
        "chr\t1\t5\nchr\tx\t5\n",
        "chr\t1\t5\nchr\t2\tlots\n",
        "chr\t1\t5\n\t2\t5\n",
        "chr\t1\t5\nchr\t0\t5\n",
        "chr\t1\t5\nchr\t2\t-1\n",
        # Every check fails in a later chunk than the next: the first check still names the file.
        "chr\t1\t-1\nchr\t2\t5\nchr\t0\t5\nchr\t3\t5\n\t4\t5\n",
        "chr\t1\t-1\nchr\t2\t5\nchr\t3\t5\nchr\tx\t5\n",
        "chr\t1\n",
    ],
)
def test_a_damaged_file_reports_what_the_row_table_reported(
    text: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(depth_module, "_CHUNK_ROWS", 2)
    path = tmp_path / "depth.tsv"
    path.write_text(text, encoding="utf-8")
    record = SeqRecord(Seq("A" * 10), id="chr")

    expected = _outcome(lambda: _reference_depth_df(record, _reference_parse(path), 5, 5))

    assert _outcome(lambda: depth_df(record, parse_depth_tsv(str(path)), 5, 5, circular=False)) == expected
    assert _outcome(lambda: read_depth_tsv(str(path)).shape) == _outcome(lambda: _reference_parse(path).shape)


def test_read_depth_tsv_returns_the_rows_of_the_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(depth_module, "_CHUNK_ROWS", 2)
    path = tmp_path / "depth.tsv"
    path.write_text("b\t3\t7\nb\t1\t300\na\t2\t1.5\nb\t3\t9\nc\t1\t2.41e-73\n", encoding="utf-8")

    expected = _reference_parse(path)
    expected["reference_name"] = expected["reference_name"].astype("string")

    pd.testing.assert_frame_equal(read_depth_tsv(str(path)), expected)


def _previous_record_prefix_sums(record: SeqRecord, depths: np.ndarray, record_len: int) -> np.ndarray:
    """The per-base running sums that ``depth_df`` built before it kept them only at window bounds."""

    base, coord_step = _read_coord_map(record)
    source_start = base if coord_step == 1 else base - record_len + 1
    prefix = np.zeros(record_len + 1, dtype=np.float64)
    first = max(source_start, 1)
    last = min(source_start + record_len - 1, int(depths.size))
    if first > last:
        return prefix
    shown = depths[first - 1 : last]
    if coord_step > 0:
        low = first - source_start
    else:
        shown = shown[::-1]
        low = source_start + record_len - 1 - last
    high = low + int(shown.size)
    sums = prefix[low + 1 : high + 1]
    sums[:] = shown
    np.cumsum(sums, out=sums)
    prefix[high + 1 :] = prefix[high]
    return prefix


def _random_depths(rng: np.random.Generator, kind: str, size: int) -> np.ndarray:
    if kind == "float":
        values = rng.random(size) * rng.choice([1e-3, 1.0, 1e6])
        values[rng.random(size) < 0.05] = -0.0
        values[rng.random(size) < 0.05] = 1e-300
    else:
        values = rng.integers(0, np.iinfo(kind).max, size=size, endpoint=True).astype(kind)
    values.flags.writeable = False
    return values


# Records whose drawn range starts at 1, past 1, before 1, or runs past the parsed positions, on either strand.
_RANGES = [(1, 1), (1, -1), (37, 1), (37, -1), (-25, 1), (-25, -1), (900, 1), (900, -1), (1300, 1), (1300, -1)]


@pytest.mark.parametrize("kind", ["float", "uint8", "uint16", "uint32"])
@pytest.mark.parametrize(("base", "coord_step"), _RANGES)
def test_window_sums_equal_the_per_base_running_sums(
    kind: str, base: int, coord_step: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every window mean is bit for bit what the per-base running sums gave, at any chunk size."""

    rng = np.random.default_rng([base & 0xFFFF, coord_step + 1, len(kind)])
    depths = _random_depths(rng, kind, 1000)
    coverage = DepthCoverage(depths={"chr": depths}, extents={"chr": (1, len(depths))})
    for chunk in (1, 7, 64, 1 << 20):
        monkeypatch.setattr(depth_module, "_SUM_CHUNK_BASES", chunk)
        for record_len, window, step in ((1, 1, 1), (999, 10, 3), (1200, 97, 97), (640, 1000, 50)):
            record = SeqRecord(Seq(None, length=record_len), id="chr")
            _write_coord_map(record, base=base if coord_step == 1 else base + record_len - 1, step=coord_step)
            prefix = _previous_record_prefix_sums(record, depths, record_len)
            _samples, starts, ends = _linear_windows(record_len, window, step)
            expected = (prefix[ends] - prefix[starts]) / (ends - starts)

            actual = depth_df(record, coverage, window, step, circular=False)["depth"].to_numpy()

            assert actual.dtype == expected.dtype
            assert np.array_equal(actual, expected)
            assert np.array_equal(np.signbit(actual), np.signbit(expected))


def test_binning_a_record_holds_no_per_base_float_array() -> None:
    """The window step keeps running sums only at window bounds, not one float64 per base."""

    length = 16_000_000
    depths = np.ones(length, dtype=np.uint8)
    depths.flags.writeable = False
    record = SeqRecord(Seq(None, length=length), id="chr")
    coverage = DepthCoverage(depths={"chr": depths}, extents={"chr": (1, len(depths))})

    tracemalloc.start()
    try:
        binned = depth_df(record, coverage, 100_000, 100_000, circular=False)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert binned["depth"].tolist() == [1.0] * 160
    # The stored depths take 1 byte per base; a per-base float64 array would take 8.
    assert peak < length
