"""Numeric series settings on every Python surface: the CLI's per-track lists, ``DepthTrackInput``, the beginner
API's ``DepthTrackOptions`` and the canonical request; the per-series axis unit; the global range of a signed series.

Each case goes through the surface a user or a Session reaches (``_get_args``, ``draw_circular``, the request codec,
a whole figure), not through a helper of it.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
import pytest

import gbdraw
import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli
from gbdraw.api import (
    CircularDiagramRequest,
    DepthTrackInput,
    GenBankInputSource,
    RecordInput,
)
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.cli_utils.common import depth_track_inputs_from_args
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request

_INPUTS = Path(__file__).resolve().parent / "test_inputs"
_HMMT = _INPUTS / "HmmtDNA.gbk"
# Every series setting away from its default.
_SETTINGS = dict(
    format="bedgraph",
    agg="max",
    style="line",
    negative_color="#c53030",
    baseline=0.5,
    min_value=-3.0,
    max_value="auto",
    unit="",
)


def _signed_bedgraph(path: Path) -> Path:
    """A log2-ratio-like bedGraph over HmmtDNA (16,569 bp): values from -2.5 to 1.5."""

    rows = ["track type=bedGraph name=ratio"]
    for index, start in enumerate(range(0, 16569, 500)):
        value = 2.0 * math.sin(index / 3.0) - 0.5
        rows.append(f"NC_012920.1\t{start}\t{min(start + 500, 16569)}\t{value:.3f}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


# --- CLI lists -------------------------------------------------------------------------------------------------------


def _cli_tracks(module, *argv: str) -> tuple[DepthTrackInput, ...] | None:
    return depth_track_inputs_from_args(module._get_args(["--gbk", "x.gb", *argv]))


@pytest.mark.parametrize("module", [circular_cli, linear_cli])
def test_one_cli_value_applies_to_every_track(module) -> None:
    tracks = _cli_tracks(
        module,
        "--depth_track", "a.bg", "--depth_track", "b.bg",
        "--depth_track_format", "bedgraph", "--depth_track_agg", "max", "--depth_track_style", "line",
        "--depth_track_negative_color", "#c53030", "--depth_track_baseline", "0.5",
        "--depth_track_min", "-3", "--depth_track_max", "auto", "--depth_track_unit", "",
    )

    assert tracks == (
        DepthTrackInput(source="a.bg", **_SETTINGS),
        DepthTrackInput(source="b.bg", **_SETTINGS),
    )


@pytest.mark.parametrize("module", [circular_cli, linear_cli])
def test_cli_values_go_to_their_tracks_in_order(module) -> None:
    tracks = _cli_tracks(
        module,
        "--depth_track", "a.tsv", "--depth_track", "b.bed",
        "--depth_track_format", "auto", "bed", "--depth_track_agg", "mean", "sum",
        "--depth_track_style", "area", "line", "--depth_track_baseline", "auto", "mean",
        "--depth_track_min", "-", "auto", "--depth_track_max", "10", "-", "--depth_track_unit", "x", "%",
    )

    assert tracks == (
        DepthTrackInput(source="a.tsv", max_value=10.0),
        DepthTrackInput(
            source="b.bed", format="bed", agg="sum", style="line", baseline="mean", min_value="auto", unit="%"
        ),
    )


def test_cli_series_values_are_checked() -> None:
    with pytest.raises(ValidationError, match="depth_track_units count must be one or equal"):
        _cli_tracks(circular_cli, "--depth_track", "a", "--depth_track", "b", "--depth_track_unit", "x", "y", "z")
    with pytest.raises(ValidationError, match="depth_track_baselines values must be finite numbers or auto or mean"):
        _cli_tracks(circular_cli, "--depth_track", "a", "--depth_track_baseline", "middle")
    with pytest.raises(SystemExit):
        _cli_tracks(circular_cli, "--depth_track", "a", "--depth_track_style", "bar")


def _refusal(call) -> tuple[str, object]:
    """The message and diagnostic of the error ``call`` raises: the first error of its cause chain with a diagnostic."""

    with pytest.raises(ValidationError) as caught:
        call()
    error: BaseException | None = caught.value
    while error is not None and getattr(error, "diagnostic", None) is None:
        error = error.__cause__
    assert error is not None, caught.value
    return str(error), getattr(error, "diagnostic")


@pytest.mark.parametrize("text", ["0", "-2", "nan", "inf", "x"])
def test_a_per_series_number_is_refused_alike_on_every_surface(tmp_path: Path, text: str) -> None:
    """One rule reads a per-series number (here a large tick interval) wherever it is given as text or a number:
    a finite number > 0, or unset. The low-level list reader took NaN and inf, the CLI's and the beginner API's
    reader inf, and each worded its refusal in its own way (OV-3887)."""

    from gbdraw.analysis.depth_tracks import _expand_track_float_metadata
    from gbdraw.api.record_planning import DepthTrackCliValues, depth_track_inputs_from_cli

    depth = tmp_path / "depth.tsv"
    depth.write_text("NC_012920.1\t1\t5\n")
    surfaces = {
        # The CLI's groups after parsing, as a Session replay reads them too (the parser refuses first, in its words).
        "CLI groups": lambda: depth_track_inputs_from_cli([[str(depth)]], DepthTrackCliValues(large_tick_intervals=[text])),
        "beginner API": lambda: gbdraw.draw_circular(
            gbdraw.read_genbank(_HMMT),
            options=gbdraw.CircularOptions(
                depth_tracks=(gbdraw.DepthTrackOptions(source=depth, large_tick_interval=text),)  # type: ignore[arg-type]
            ),
        ),
        "request option lists": lambda: encode_canonical_request(
            CircularDiagramRequest(
                records=(RecordInput(source=GenBankInputSource(_HMMT)),),
                options=CircularDiagramOptions(depth_track_files=[[str(depth)]], depth_track_large_tick_intervals=[text]),
            )
        ),
        "assembly option lists": lambda: _expand_track_float_metadata(
            [text], track_count=1, field_name="depth_track_large_tick_intervals"
        ),
    }
    expected = (
        "depth_track_large_tick_intervals values must be finite numbers > 0 or auto.",
        {"code": "INPUT_INVALID", "field": "large_tick_interval", "reason": "POSITIVE_OR_AUTO"},
    )

    assert {name: _refusal(call) for name, call in surfaces.items()} == dict.fromkeys(surfaces, expected)


@pytest.mark.parametrize(
    ("module", "option"),
    [
        (circular_cli, "--depth_track_large_tick_interval"),
        (circular_cli, "--depth_track_small_tick_interval"),
        (circular_cli, "--depth_track_tick_font_size"),
        (linear_cli, "--depth_track_large_tick_interval"),
        (linear_cli, "--depth_track_small_tick_interval"),
        (linear_cli, "--depth_track_tick_font_size"),
        (linear_cli, "--depth_track_height"),
    ],
)
@pytest.mark.parametrize("text", ["0", "-2", "nan", "inf", "x"])
def test_the_cli_refuses_a_per_series_number_by_the_same_rule_naming_its_option(
    module, option: str, text: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """The CLI reads each per-series number with the rule every other surface uses and names its own option."""

    assert getattr(module._get_args(["--gbk", "x.gb", f"{option}=auto"]), option[2:]) == ["auto"]
    with pytest.raises(SystemExit):
        module._get_args(["--gbk", "x.gb", f"{option}={text}"])
    assert f"{option} values must be finite numbers > 0 or auto" in capsys.readouterr().err


@pytest.mark.parametrize("module", [circular_cli, linear_cli])
def test_the_global_depth_range_may_be_negative_except_on_a_log_scale(
    module, capsys: pytest.CaptureFixture[str]
) -> None:
    def validated(*argv: str) -> None:
        args = module._get_args(["--gbk", "x.gb", *argv])
        from gbdraw.cli_utils.common import validate_depth_axis_args

        parser = __import__("argparse").ArgumentParser()
        validate_depth_axis_args(parser, args)

    validated("--depth_min", "-2", "--depth_max", "1.5")
    with pytest.raises(SystemExit):
        validated("--depth_min", "2", "--depth_max", "1")
    assert "--depth_min must be <= --depth_max" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        validated("--depth_min", "-2", "--depth_log_scale")
    assert "--depth_min must be >= 0 with --depth_log_scale" in capsys.readouterr().err


def test_a_drawing_refuses_a_global_depth_range_out_of_order_as_the_series_reader_does(tmp_path: Path) -> None:
    """The assembler's check of the global range is the series reader's (``check_depth_range``): one message, which
    the Web reads as the ORDER of min_depth; it had its own wording, which the Web could not classify."""

    from gbdraw.web_support.error_adapter import serialize_web_error

    depth = tmp_path / "depth.tsv"
    depth.write_text("NC_012920.1\t1\t5\n")
    overrides = {"objects.depth.min_depth": 5, "objects.depth.max_depth": 1}
    with pytest.raises(ValidationError, match=r"^min_depth must be <= max_depth\.$") as caught:
        gbdraw.draw_circular(
            gbdraw.read_genbank(_HMMT),
            options=gbdraw.CircularOptions(
                depth_tracks=(gbdraw.DepthTrackOptions(source=depth),), config_overrides=overrides
            ),
        ).to_svg()
    payload = serialize_web_error(caught.value, operation="generate", stage="render")
    assert (payload["code"], payload["context"]) == ("INPUT_INVALID", {"field": "min_depth", "reason": "ORDER"})


# --- DepthTrackInput -------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"format": "wig"}, "DepthTrackInput.format must be one of: auto, samtools, bedgraph, bed, position_value"),
        ({"agg": "median"}, "DepthTrackInput.agg must be one of: mean, max, sum"),
        ({"style": "bar"}, "DepthTrackInput.style must be one of: area, line"),
        ({"baseline": "top"}, 'DepthTrackInput.baseline must be a finite number or "auto" or "mean"'),
        ({"baseline": math.nan}, "DepthTrackInput.baseline must be a finite number"),
        ({"min_value": math.inf}, 'DepthTrackInput.min_value must be a finite number or "auto"'),
        ({"max_value": True}, 'DepthTrackInput.max_value must be a finite number or "auto"'),
        ({"min_value": 2, "max_value": 1}, "DepthTrackInput.min_value must be <= DepthTrackInput.max_value"),
        ({"unit": None}, "DepthTrackInput.unit must be a string"),
        ({"negative_color": " "}, "DepthTrackInput.negative_color must be a non-empty string or None"),
    ],
)
def test_depth_track_input_refuses_a_setting_out_of_its_kind(settings: dict[str, object], message: str) -> None:
    with pytest.raises(ValidationError, match=re.escape(message)):
        DepthTrackInput(source="a.tsv", **settings)  # type: ignore[arg-type]


def test_a_log_scaled_series_refuses_a_negative_range(tmp_path: Path) -> None:
    track = DepthTrackInput(source=_signed_bedgraph(tmp_path / "s.bg"), min_value=-1.0, max_value=2.0)
    records = gbdraw.read_genbank(_HMMT)
    with pytest.raises(ValidationError, match="A log-scaled depth range must not be negative") as refused:
        _render_api(records, (track,), normalize=True)
    assert refused.value.diagnostic == {
        "code": "INPUT_INVALID", "field": "min_depth", "reason": "NONNEGATIVE", "seriesIndex": 0
    }


def _render_api(records, tracks, *, normalize: bool = False) -> str:
    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.api.assembly_inputs import CircularAssemblyInputs
    from gbdraw.api.config import apply_config_overrides

    cfg = apply_config_overrides(None, {"objects.depth.normalize": normalize})
    canvas = assemble_circular_diagram_from_record(
        records[0], CircularAssemblyInputs(depth_tracks=tracks, cfg=cfg, legend="right")
    )
    return canvas.tostring()


# --- The surfaces draw one figure ------------------------------------------------------------------------------------


def test_cli_options_and_inputs_draw_the_same_figure(tmp_path: Path) -> None:
    bedgraph = _signed_bedgraph(tmp_path / "signed.bg")
    circular_cli.circular_main([
        "--gbk", str(_HMMT), "-o", str(tmp_path / "cli"), "-f", "svg",
        "--depth_track", str(bedgraph), "--depth_track_label", "Ratio", "--depth_track_color", "#2b6cb0",
        *[f"--depth_track_{key}={value}" for key, value in (
            ("format", "bedgraph"), ("agg", "max"), ("style", "line"), ("negative_color", "#c53030"),
            ("baseline", "0.5"), ("min", "-3"), ("max", "auto"), ("unit", ""),
        )],
    ])
    cli_svg = (tmp_path / "cli.svg").read_text(encoding="utf-8")
    records = gbdraw.read_genbank(_HMMT)
    options_svg = gbdraw.draw_circular(
        records,
        options=gbdraw.CircularOptions(
            depth_tracks=(gbdraw.DepthTrackOptions(source=bedgraph, label="Ratio", color="#2b6cb0", **_SETTINGS),),
        ),
    ).to_svg()
    request = CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(_HMMT)),),
        options=CircularDiagramOptions(
            depth_tracks=(DepthTrackInput(source=bedgraph, label="Ratio", color="#2b6cb0", **_SETTINGS),)
        ),
    )
    input_svg = gbdraw.api.build_request_diagram(request).drawing.tostring()

    assert options_svg == input_svg
    assert _figure(cli_svg) == _figure(input_svg)
    # The settings reach the drawing: a line in the series colour, the part below 0.5 in the negative colour, and
    # the range label without a unit.
    assert 'stroke="#c53030"' in input_svg and 'stroke="#2b6cb0"' in input_svg
    assert re.findall(r">(\[[^<]*\])<", input_svg) == ["[-3–1.49]"]


def _figure(svg: str) -> str:
    """The SVG drawing without the XML declaration a file starts with."""

    return svg.removeprefix('<?xml version="1.0" encoding="utf-8" ?>\n')


@pytest.mark.parametrize(("unit", "label"), [("x", "[-2.5–1.49x]"), ("", "[-2.5–1.49]"), (" %", "[-2.5–1.49 %]")])
def test_a_series_axis_reads_its_unit(tmp_path: Path, unit: str, label: str) -> None:
    track = DepthTrackInput(source=_signed_bedgraph(tmp_path / "s.bg"), unit=unit)
    svg = _render_api(gbdraw.read_genbank(_HMMT), (track,))

    assert re.findall(r">(\[[^<]*\])<", svg) == [label]


def test_one_file_read_in_two_formats_is_two_series(tmp_path: Path) -> None:
    # A bedGraph row gives its value to the bases it covers; a BED row without a score adds 1.
    path = tmp_path / "rows.txt"
    path.write_text("NC_012920.1\t0\t16569\t40\n", encoding="utf-8")
    tracks = (
        DepthTrackInput(source=path, format="bedgraph", label="as bedGraph"),
        DepthTrackInput(source=path, format="bed", label="as BED"),
    )
    labels = re.findall(r">(\[[^<]*\])<", _render_api(gbdraw.read_genbank(_HMMT), tracks))

    assert labels == ["[0–40x]", "[0–1x]"]


def test_beginner_depth_tracks_without_label_or_colour_draw_the_defaults(tmp_path: Path) -> None:
    # The beginner API handed [None] lists to the compatibility depth fields, which read them as the text "None":
    # a Legend row "None" and a series with fill="None", which draws nothing (OV-3748).
    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"NC_012920.1\t{position}\t{position % 50}\n" for position in range(1, 16570, 7)))
    svg = gbdraw.draw_circular(
        gbdraw.read_genbank(_HMMT),
        options=gbdraw.CircularOptions(depth_tracks=(gbdraw.DepthTrackOptions(source=depth),)),
    ).to_svg()

    assert 'fill="None"' not in svg and ">None<" not in svg
    assert 'data-legend-key="Depth"' in svg


# --- Canonical request -----------------------------------------------------------------------------------------------


def _request(*tracks: DepthTrackInput) -> CircularDiagramRequest:
    return CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(_HMMT)),),
        options=CircularDiagramOptions(depth_tracks=tracks),
    )


def _decoded(payload, encoded, tmp_path: Path):
    directory = tmp_path / "materialized"
    directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for resource in encoded.resources:
        if resource.source_path is not None:
            paths[resource.resource_id] = resource.source_path
            continue
        target = directory / resource.name
        target.write_bytes(resource.content or b"")
        paths[resource.resource_id] = target
    return decode_canonical_request(payload, resource_paths=paths, output_directory=tmp_path / "out")


def test_every_series_setting_round_trips_through_the_request(tmp_path: Path) -> None:
    bedgraph = _signed_bedgraph(tmp_path / "s.bg")
    encoded = encode_canonical_request(_request(DepthTrackInput(source=bedgraph, **_SETTINGS)))
    stored = encoded.payload["diagramOptions"]["depthTracks"][0]

    assert {key: stored[key] for key in (
        "format", "agg", "style", "negativeColor", "baseline", "minValue", "maxValue", "unit"
    )} == {
        "format": "bedgraph", "agg": "max", "style": "line", "negativeColor": "#c53030", "baseline": 0.5,
        "minValue": -3.0, "maxValue": "auto", "unit": "",
    }
    decoded = _decoded(encoded.payload, encoded, tmp_path)
    track = decoded.options.depth_tracks[0]
    assert {name: getattr(track, name) for name in _SETTINGS} == _SETTINGS
    assert encode_canonical_request(decoded).payload == encoded.payload


def test_a_track_with_default_settings_writes_no_setting_and_an_older_request_reads_the_defaults(
    tmp_path: Path,
) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text("NC_012920.1\t1\t5\n", encoding="utf-8")
    encoded = encode_canonical_request(_request(DepthTrackInput(source=depth, label="Coverage")))
    stored = encoded.payload["diagramOptions"]["depthTracks"][0]

    # The keys a request written before these settings holds, and no more.
    assert set(stored) == {
        "source", "label", "color", "height", "largeTickInterval", "smallTickInterval", "tickFontSize",
    }
    track = _decoded(encoded.payload, encoded, tmp_path).options.depth_tracks[0]
    assert (track.format, track.agg, track.style, track.negative_color, track.baseline) == (
        "auto", "mean", "area", None, "auto"
    )
    assert (track.min_value, track.max_value, track.unit) == (None, None, "x")


def test_a_request_with_an_unknown_depth_track_key_is_refused(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text("NC_012920.1\t1\t5\n", encoding="utf-8")
    encoded = encode_canonical_request(_request(DepthTrackInput(source=depth)))
    payload = {**encoded.payload}
    payload["diagramOptions"] = {
        **payload["diagramOptions"],
        "depthTracks": [{**payload["diagramOptions"]["depthTracks"][0], "ymin": 0}],
    }
    with pytest.raises(Exception, match="ymin"):
        _decoded(payload, encoded, tmp_path)


def test_the_request_description_declares_the_settings_optional() -> None:
    from gbdraw.session_request_codec import describe_canonical_request

    described = describe_canonical_request()
    depth_track = described["types"]["CanonicalCircularDepthTrack"]
    assert set(depth_track["optional"]) == {
        "format", "agg", "style", "negativeColor", "baseline", "minValue", "maxValue", "unit",
    }


def test_a_table_source_keeps_its_settings(tmp_path: Path) -> None:
    table = pd.DataFrame({"reference_name": ["NC_012920.1"] * 3, "position": [1, 2, 3], "depth": [-1.0, 0.0, 2.0]})
    encoded = encode_canonical_request(_request(DepthTrackInput(source=table, style="line", unit="")))
    track = _decoded(encoded.payload, encoded, tmp_path).options.depth_tracks[0]

    assert (track.style, track.unit) == ("line", "")


# --- Circular figure bounds ------------------------------------------------------------------------------------------


def test_a_circular_figure_keeps_room_for_a_line_series_stroke() -> None:
    # The stroke safety of a Circular figure read the base Depth paint, whose stroke is the Depth outline (none
    # here); a line series strokes its value line in its colour at the Depth stroke width, so a wide one painted past
    # the figure's bounds (OV-3747).
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api.assembly_inputs import CircularAssemblyInputs
    from gbdraw.api.config import apply_config_overrides
    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.render.composition import applied_composition_plan

    record = SeqRecord(Seq("ATGC" * 30), id="rec1", name="rec1")
    record.annotations["molecule_type"] = "DNA"
    table = pd.DataFrame({"reference_name": ["rec1"] * 4, "position": [1, 2, 10, 30], "depth": [10, 20, 40, 80]})
    cfg = apply_config_overrides(None, {"objects.depth.stroke_color": "none", "objects.depth.stroke_width": 20.0})

    def primary_width(style: str) -> float:
        canvas = assemble_circular_diagram_from_record(
            record,
            CircularAssemblyInputs(
                depth_tracks=(DepthTrackInput(source=table, style=style),),  # type: ignore[arg-type]
                cfg=cfg,
                legend="none",
                circular_track_slots=["depth:depth", "features:features"],
            ),
        )
        primary = next(item for item in applied_composition_plan(canvas).placements if item.role == "primary")
        return float(primary.final_bounds.width)

    # Half the 20 px line (10 px) past the outer ring on both sides, where the other strokes keep 2.5 px.
    assert primary_width("line") - primary_width("area") == pytest.approx(2 * (10.0 - 2.5))


# --- Session ---------------------------------------------------------------------------------------------------------


SERIES_SESSION = Path(__file__).resolve().parent / "fixtures/sessions/HmmtDNA-depth-series-settings.v51-cli.gbdraw"


def test_a_cli_session_with_series_settings_replays_as_it_was_drawn(tmp_path: Path) -> None:
    from gbdraw.api import load_session_document
    from gbdraw.session import session_drawing_artifacts

    drawing = session_drawing_artifacts(load_session_document(SERIES_SESSION), "circular").fields
    stored = drawing["renderRequest"]["diagramOptions"]["depthTracks"]
    assert [{key: track.get(key) for key in ("format", "agg", "style", "baseline", "minValue", "unit")} for track in stored] == [
        {"format": "bedgraph", "agg": None, "style": None, "baseline": None, "minValue": -3.0, "unit": ""},
        {"format": None, "agg": "max", "style": "line", "baseline": 0.5, "minValue": "auto", "unit": None},
    ]

    circular_cli.circular_main(["--session", str(SERIES_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"])

    svg = (tmp_path / "replay.svg").read_text(encoding="utf-8")
    assert re.findall(r">(\[[^<]*\])<", svg) == ["[-3–2]", "[-2.5–1.49x]"]
    assert 'stroke="#2f855a"' in svg and 'fill="#c53030"' in svg
