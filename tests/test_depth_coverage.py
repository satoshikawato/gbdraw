"""The parsed depth form: one dense array per reference, binned exactly as the row table was.

``_reference_*`` below is the row-table algorithm that ``depth_df`` used before
the parsed form became compact. It stays here only, as the oracle the
randomized cases compare against value for value.
"""

from __future__ import annotations

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


def _reference_depth_df(record, table, window, step, *, normalize=False, min_depth=None, max_depth=None):
    table = table.loc[:, _COLUMNS].copy()
    table["reference_name"] = table["reference_name"].astype("string")
    _reference_validate(table)
    record_len = len(record.seq)
    if record_len <= 0:
        return pd.DataFrame(columns=["position", "depth", "depth_normalized"])
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
    starts = np.arange(0, record_len, int(step), dtype=np.int64)
    ends = np.minimum(starts + int(window), record_len)
    means = (prefix[ends] - prefix[starts]) / np.maximum(ends - starts, 1)
    if min_depth is not None or max_depth is not None:
        means = np.clip(means, min_depth, max_depth)
    if normalize:
        def transform(values):
            values = np.asarray(values, dtype=np.float64)
            out = np.zeros_like(values)
            positive = values > 0
            out[positive] = np.log10(values[positive]) + 1.0
            return out

        low = float(transform(float(min_depth) if min_depth is not None else 0.0))
        high = float(transform(float(max_depth) if max_depth is not None else float(np.max(means, initial=0.0))))
        normalized = np.zeros_like(means) if high <= low else np.clip((transform(means) - low) / (high - low), 0.0, 1.0)
    else:
        normalized = means.astype(np.float64, copy=True)
    return pd.DataFrame({"position": starts, "depth": means, "depth_normalized": normalized})


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
        "normalize": bool(rng.random() < 0.5),
        "min_depth": min_depth,
        "max_depth": max_depth,
    }


@pytest.mark.parametrize("seed", range(10))
def test_binned_depth_equals_the_row_table_algorithm(seed: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Parsed files and API tables bin exactly as the row table did, also across chunk boundaries."""

    rng = np.random.default_rng(seed)
    for case in range(30):
        monkeypatch.setattr(depth_module, "_CHUNK_ROWS", int(rng.choice([11, 97, 1 << 20])))
        path = tmp_path / f"depth-{case}.tsv"
        text = _depth_text(rng)
        path.write_text(text, encoding="utf-8")
        record = _record(rng, text)
        binning = _binning(rng)
        args = (binning["window"], binning["step"])
        kwargs = {key: binning[key] for key in ("normalize", "min_depth", "max_depth")}
        expected = _outcome(lambda: _reference_depth_df(record, _reference_parse(path), *args, **kwargs))
        _assert_same(_outcome(lambda: depth_df(record, parse_depth_tsv(str(path)), *args, **kwargs)), expected)
        _assert_same(_outcome(lambda: depth_df(record, read_depth_tsv(str(path)), *args, **kwargs)), expected)


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

    assert _outcome(lambda: depth_df(record, parse_depth_tsv(str(path)), 5, 5)) == expected
    assert _outcome(lambda: read_depth_tsv(str(path)).shape) == _outcome(lambda: _reference_parse(path).shape)


def test_read_depth_tsv_returns_the_rows_of_the_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(depth_module, "_CHUNK_ROWS", 2)
    path = tmp_path / "depth.tsv"
    path.write_text("b\t3\t7\nb\t1\t300\na\t2\t1.5\nb\t3\t9\nc\t1\t2.41e-73\n", encoding="utf-8")

    expected = _reference_parse(path)
    expected["reference_name"] = expected["reference_name"].astype("string")

    pd.testing.assert_frame_equal(read_depth_tsv(str(path)), expected)
