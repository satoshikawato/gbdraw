"""Numeric series features of a Depth track: the reader's explicit format, signed values, the window aggregate, the
per-series value range and baseline, the line style and the colour below the baseline.

Each case reaches the feature through the engine: the reader (``parse_depth_tsv``), the windows (``depth_df``), a
``DepthTrackSpec`` and the drawing (``build_depth_track_dataframes``, ``depth_series`` and a whole figure).
"""

from __future__ import annotations

import re
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.depth import DepthCoverage, depth_df, parse_depth_tsv
from gbdraw.analysis.depth_tracks import DepthTrackSpec, build_depth_track_dataframes, freeze_depth_tracks
from gbdraw.analysis.windows import sliding_windows
from gbdraw.api.assembly_inputs import CircularAssemblyInputs, LinearAssemblyInputs
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_circular_diagram_from_record, assemble_linear_diagram_from_records
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import DepthConfigurator
from gbdraw.exceptions import GbdrawError, ParseError, ValidationError
from gbdraw.render.numeric_series import depth_series

_SVG = "{http://www.w3.org/2000/svg}"


def _record(record_id: str = "rec1", length: int = 40) -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * (length // 4)), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    return record


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _depths(coverage: DepthCoverage, name: str = "rec1") -> list[float]:
    return [float(value) for value in coverage.depths[name]]


def _base_config(**overrides: object) -> DepthConfigurator:
    profile = LinearRenderProfile(GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")))
    return DepthConfigurator(10, 10, profile, **overrides)  # type: ignore[arg-type]


# (1) samtools depth of several samples without its header line.


def test_a_samtools_file_of_two_samples_reads_its_first_sample(tmp_path: Path) -> None:
    path = _write(tmp_path / "two_samples.tsv", "rec1\t1\t5\t9\nrec1\t2\t7\t1\nrec1\t4\t2\t8\n")

    coverage = parse_depth_tsv(str(path), file_format="samtools")

    assert _depths(coverage) == [5.0, 7.0, 0.0, 2.0]
    assert coverage.extents["rec1"] == (1, 4)
    # Read by its column count, the same file is a bedGraph whose rows are not intervals.
    with pytest.raises(ValidationError) as refused:
        parse_depth_tsv(str(path))
    assert refused.value.diagnostic == {"code": "DEPTH_INVALID", "reason": "INTERVAL"}


def test_a_samtools_file_keeps_its_header_line(tmp_path: Path) -> None:
    path = _write(tmp_path / "header.tsv", "#CHROM\tPOS\ta.bam\tb.bam\nrec1\t1\t3\t9\nrec1\t2\t4\t1\n")

    assert _depths(parse_depth_tsv(str(path), file_format="samtools")) == [3.0, 4.0]


def test_an_explicit_format_with_too_few_columns_is_refused(tmp_path: Path) -> None:
    path = _write(tmp_path / "two.tsv", "1\t5\n2\t6\n")

    with pytest.raises(ParseError) as refused:
        parse_depth_tsv(str(path), file_format="bedgraph")
    assert refused.value.diagnostic == {"code": "TABLE_INVALID", "reason": "DEPTH_COLUMNS"}
    assert _depths(parse_depth_tsv(str(path), file_format="position_value"), "") == [5.0, 6.0]


def test_an_unknown_format_is_refused(tmp_path: Path) -> None:
    path = _write(tmp_path / "one.tsv", "rec1\t1\t5\n")

    with pytest.raises(ValidationError, match="format"):
        parse_depth_tsv(str(path), file_format="wig")  # type: ignore[arg-type]


# (2) BED: the score column, no score, overlapping and unsorted rows.


def test_bed_rows_add_their_score_to_every_base_they_cover(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "features.bed",
        "track name=hits\nrec1\t10\t20\tA\t5\nrec1\t0\t5\tB\t2\nrec1\t15\t25\tC\t3\n",
    )

    coverage = parse_depth_tsv(str(path), file_format="bed")

    expected = [2.0] * 5 + [0.0] * 5 + [5.0] * 5 + [8.0] * 5 + [3.0] * 5
    assert _depths(coverage) == expected
    assert coverage.extents["rec1"] == (1, 25)
    assert coverage.depths["rec1"].dtype == np.uint8


def test_bed_rows_without_a_score_count_one_each(tmp_path: Path) -> None:
    three = _write(tmp_path / "three.bed", "rec1\t2\t3\nrec1\t0\t3\nrec1\t2\t4\n")
    named = _write(tmp_path / "named.bed", "rec1\t2\t3\tsnp1\nrec1\t0\t3\tsnp2\nrec1\t2\t4\tsnp3\n")
    dotted = _write(tmp_path / "dotted.bed", "rec1\t2\t3\tsnp1\t.\nrec1\t0\t3\tsnp2\t.\nrec1\t2\t4\tsnp3\t.\n")

    for path in (three, named, dotted):
        assert _depths(parse_depth_tsv(str(path), file_format="bed")) == [1.0, 1.0, 3.0, 1.0], path.name


def test_bed_point_rows_summed_per_window_count_them(tmp_path: Path) -> None:
    rows = "".join(f"rec1\t{start}\t{start + 1}\tsnp\t0\n" for start in (3, 3, 7, 25, 1))
    path = _write(tmp_path / "snps.bed", rows.replace("\t0\n", "\t1\n"))

    frame = depth_df(_record(length=40), parse_depth_tsv(str(path), file_format="bed"), 10, 10, circular=False,
                     agg="sum")

    # Windows are centred on 0, 10, 20, 30: [0, 5), [5, 15), [15, 25), [25, 35).
    assert frame["depth"].tolist() == [3.0, 1.0, 0.0, 1.0]


# (3) Signed values.


def test_a_signed_bedgraph_reads_its_negative_values(tmp_path: Path) -> None:
    path = _write(tmp_path / "signed.bedgraph", "rec1\t0\t20\t-3\nrec1\t20\t40\t4\n")

    coverage = parse_depth_tsv(str(path))

    assert coverage.depths["rec1"].dtype == np.int8
    assert _depths(coverage) == [-3.0] * 20 + [4.0] * 20
    big = _write(tmp_path / "big.bedgraph", "rec1\t0\t2\t-3\nrec1\t2\t4\t300\nrec1\t4\t5\t-0.5\n")
    assert parse_depth_tsv(str(big), file_format="bedgraph").depths["rec1"].dtype == np.float64
    whole = _write(tmp_path / "whole.bedgraph", "rec1\t0\t2\t-3\nrec1\t2\t4\t300\n")
    assert parse_depth_tsv(str(whole)).depths["rec1"].dtype == np.int16


def test_a_signed_series_draws_below_and_above_zero() -> None:
    frame = pd.DataFrame({"position": [0, 10, 20, 30], "depth": [-3.0, -1.0, 2.0, 4.0]})
    config = _base_config()
    config.negative_color = "#d62728"

    plot = depth_series(config, frame)

    assert plot.axis is not None
    assert (plot.axis.axis_min, plot.axis.axis_max) == (-3.0, 4.0)
    assert plot.samples is not None
    assert plot.samples.baseline == pytest.approx(3.0 / 7.0)
    assert plot.paint.fill == config.fill_color
    assert plot.paint.below_fill == "#d62728"


def test_global_depth_bounds_may_be_negative() -> None:
    coverage = DepthCoverage({"rec1": np.array([-8, -2, 3, 9], dtype=np.int8)}, {"rec1": (1, 4)})

    frame = depth_df(_record(length=4), coverage, 1, 1, circular=False, min_depth=-5, max_depth=5)

    assert frame["depth"].tolist() == [-5.0, -2.0, 3.0, 5.0]
    with pytest.raises(ValidationError, match="min_depth must be <= max_depth"):
        depth_df(_record(length=4), coverage, 1, 1, circular=False, min_depth=2, max_depth=-2)


# (4) The window aggregate.


def _brute_windows(record: SeqRecord, per_base: list[float], window: int, step: int, agg: str) -> list[float]:
    windows = sliding_windows(record, window, step, circular=False)
    values = np.array(per_base, dtype=np.float64)
    reduce = {"sum": np.sum, "max": np.max, "mean": np.mean}[agg]
    return [float(reduce(values[start:end])) for start, end in zip(windows.starts, windows.ends)]


@pytest.mark.parametrize("agg", ["mean", "max", "sum"])
def test_window_values_follow_the_aggregate(agg: str, tmp_path: Path) -> None:
    rows = {1: 1, 2: 5, 3: 2, 5: 3, 6: 7, 7: -4, 9: 6, 12: 2}
    path = _write(tmp_path / "values.tsv", "".join(f"rec1\t{p}\t{v}\n" for p, v in rows.items()))
    record = _record(length=12)
    per_base = [float(rows.get(position, 0)) for position in range(1, 13)]

    frame = depth_df(record, parse_depth_tsv(str(path)), 4, 3, circular=False, agg=agg)

    assert frame["depth"].tolist() == pytest.approx(_brute_windows(record, per_base, 4, 3, agg))


def test_the_largest_value_of_a_window_counts_a_left_out_base_as_zero() -> None:
    coverage = DepthCoverage({"rec1": np.array([-5, 0, -3, -2], dtype=np.int8)}, {"rec1": (1, 4)})

    frame = depth_df(_record(length=4), coverage, 2, 2, circular=False, agg="max")

    # Windows [0, 1) and [1, 3) of the samples 0 and 2.
    assert frame["depth"].tolist() == [-5.0, 0.0]


def test_the_largest_value_of_a_circular_window_reads_round_the_origin() -> None:
    coverage = DepthCoverage({"rec1": np.array([1, 2, 3, 4, 5, 6, 7, 9], dtype=np.uint8)}, {"rec1": (1, 8)})

    frame = depth_df(_record(length=8), coverage, 4, 4, circular=True, agg="max")

    # Windows [-2, 2) and [2, 6): the first reads positions 7, 8, 1 and 2.
    assert frame["depth"].tolist() == [9.0, 6.0]


def test_an_unknown_aggregate_is_refused() -> None:
    coverage = DepthCoverage({"rec1": np.array([1, 2], dtype=np.uint8)}, {"rec1": (1, 2)})

    with pytest.raises(ValidationError, match="agg"):
        depth_df(_record(length=4), coverage, 2, 2, circular=False, agg="median")  # type: ignore[arg-type]


# (5) The baseline and the per-series value range.


@pytest.mark.parametrize(
    ("baseline", "expected"),
    [("auto", 0.0), ("mean", (2.5 - 1.0) / 3.0), (3.0, 2.0 / 3.0), (10.0, 1.0), (-10.0, 0.0)],
)
def test_the_baseline_is_auto_the_mean_or_a_number_within_the_range(baseline: object, expected: float) -> None:
    frame = pd.DataFrame({"position": [0, 10, 20, 30], "depth": [1.0, 2.0, 3.0, 4.0]})
    config = _base_config(min_depth=1.0)
    config.baseline = baseline  # type: ignore[assignment]

    plot = depth_series(config, frame)

    assert plot.samples is not None
    assert plot.samples.baseline == pytest.approx(expected)


def test_a_series_range_overrides_the_global_range() -> None:
    records = [_record("rec1", length=40)]
    coverage = DepthCoverage({"rec1": np.arange(-20, 20, dtype=np.int8)}, {"rec1": (1, 40)})
    specs = [
        [
            DepthTrackSpec(id="depth_1", label="A", table=coverage, track_index=0, min_value=-5.0, max_value=5.0),
            DepthTrackSpec(id="depth_2", label="B", table=coverage, track_index=1),
            DepthTrackSpec(id="depth_3", label="C", table=coverage, track_index=2, min_value="auto",
                           max_value="auto"),
        ]
    ]

    data = build_depth_track_dataframes(records, specs, None, base_config=_base_config(min_depth=-10.0,
                                                                                       max_depth=10.0))

    first, second, third = data[0]
    assert (first.config.min_depth, first.config.max_depth) == (-5.0, 5.0)
    assert first.df["depth"].min() == -5.0 and first.df["depth"].max() == 5.0
    assert (second.config.min_depth, second.config.max_depth) == (-10.0, 10.0)
    assert (third.config.min_depth, third.config.max_depth) == (None, None)
    assert third.df["depth"].min() < -10.0


def test_a_shared_axis_shares_the_signed_auto_bounds() -> None:
    records = [_record("rec1", length=40), _record("rec2", length=40)]
    tables = {
        "rec1": DepthCoverage({"rec1": np.full(40, -6, dtype=np.int8)}, {"rec1": (1, 40)}),
        "rec2": DepthCoverage({"rec2": np.full(40, 3, dtype=np.int8)}, {"rec2": (1, 40)}),
    }
    specs = [[DepthTrackSpec(id="depth", label="A", table=tables[record.id])] for record in records]

    data = build_depth_track_dataframes(records, specs, None, base_config=_base_config(share_axis=True))

    for row in data:
        plot = depth_series(row[0].config, row[0].df)
        assert plot.axis is not None
        assert (plot.axis.axis_min, plot.axis.axis_max) == (-6.0, 3.0)


# (6) The line style.


def _spec_svg(mode: str, spec: DepthTrackSpec, *, legend: str = "none", **cfg: object) -> str:
    record = _record("rec1", length=400)
    overrides = {"canvas.show_gc": False, "canvas.show_skew": False, **cfg}
    common = dict(
        legend=legend,
        cfg=apply_config_overrides(None, overrides),
        window=40,
        step=40,
        depth_window=20,
        depth_step=10,
        _resolved_depth_tracks=freeze_depth_tracks([[spec]]),
    )
    if mode == "circular":
        return assemble_circular_diagram_from_record(record, CircularAssemblyInputs(**common)).tostring()
    return assemble_linear_diagram_from_records([record], LinearAssemblyInputs(**common)).tostring()


def _signed_coverage() -> DepthCoverage:
    values = np.round(8 * np.sin(np.arange(400) / 30.0)).astype(np.int8)
    return DepthCoverage({"rec1": values}, {"rec1": (1, 400)})


def _series_paths(svg: str) -> list[ET.Element]:
    root = ET.fromstring(svg)
    (slot,) = [element for element in root.iter(f"{_SVG}g") if element.get("data-gbdraw-slot-id") == "depth"]
    return [element for element in slot.iter(f"{_SVG}path") if element.get("fill") != "white"
            and not (element.get("id") or "").endswith("_axis")]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_line_series_strokes_its_value_line_without_a_fill(mode: str) -> None:
    spec = DepthTrackSpec(id="depth", label="Signal", table=_signed_coverage(), fill_color="#1f77b4", style="line",
                          negative_color="#d62728")

    paths = _series_paths(_spec_svg(mode, spec))

    above, below = [path for path in paths if path.get("stroke") in {"#1f77b4", "#d62728"}]
    assert (above.get("fill"), above.get("stroke"), above.get("stroke-width")) == ("none", "#1f77b4", "1.0")
    assert (below.get("fill"), below.get("stroke"), below.get("stroke-width")) == ("none", "#d62728", "1.0")
    assert below.get("clip-path", "").startswith("url(#")
    assert above.get("d") == below.get("d")
    # An open line: no close back to the baseline.
    assert "z" not in (above.get("d") or "").lower()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_line_series_takes_the_depth_stroke_width_when_it_is_set(mode: str) -> None:
    spec = DepthTrackSpec(id="depth", label="Signal", table=_signed_coverage(), fill_color="#1f77b4", style="line")

    paths = _series_paths(_spec_svg(mode, spec, **{"objects.depth.stroke_width": 2.5}))

    (line,) = [path for path in paths if path.get("stroke") == "#1f77b4"]
    assert (line.get("fill"), line.get("stroke"), line.get("stroke-width")) == ("none", "#1f77b4", "2.5")


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_an_area_series_fills_below_its_baseline_in_the_negative_colour(mode: str) -> None:
    spec = DepthTrackSpec(id="depth", label="Signal", table=_signed_coverage(), fill_color="#1f77b4",
                          negative_color="#d62728")

    svg = _spec_svg(mode, spec, legend="left")
    paths = _series_paths(svg)

    above, below = [path for path in paths if path.get("fill") in {"#1f77b4", "#d62728"}]
    assert above.get("fill") == "#1f77b4" and below.get("fill") == "#d62728"
    assert below.get("clip-path", "").startswith("url(#")
    # One Legend row for the series, in its colour (a Linear figure writes its Legend laid out both ways).
    root = ET.fromstring(svg)
    rows = [element for element in root.iter(f"{_SVG}g") if element.get("data-legend-key") is not None]
    assert {row.get("data-legend-key") for row in rows} == {"Signal"}
    assert len(rows) == (1 if mode == "circular" else 2)
    assert {row.find(f"{_SVG}path").get("fill") for row in rows} == {"#1f77b4"}  # type: ignore[union-attr]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_signed_axis_labels_its_lower_bound(mode: str) -> None:
    spec = DepthTrackSpec(id="depth", label="Signal", table=_signed_coverage())

    svg = _spec_svg(mode, spec, **{"objects.depth.axis_labels": "numbered"})

    assert re.search(r">-[0-9.]+x<", svg)


# (7) The log scale refuses a negative value.


def test_a_log_scaled_series_refuses_a_negative_value() -> None:
    records = [_record("rec1", length=40)]
    coverage = DepthCoverage({"rec1": np.arange(-20, 20, dtype=np.int8)}, {"rec1": (1, 40)})

    with pytest.raises(GbdrawError) as refused:
        build_depth_track_dataframes(
            records,
            [[DepthTrackSpec(id="depth", label="A", table=coverage)]],
            None,
            base_config=_base_config(normalize=True),
        )
    # The failure names the series, as a reader failure does.
    # Its reason names the log scale that refuses the values (review r1 L2).
    assert refused.value.diagnostic == {"code": "DEPTH_INVALID", "reason": "LOG_SCALE_NEGATIVE", "seriesIndex": 0}
    assert "log-scaled" in str(refused.value)
