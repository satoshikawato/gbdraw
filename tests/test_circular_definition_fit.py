"""Circular center-definition fit (PV-08, D-24 / PD-OI-078 scenario revision 2).

The centre definition never narrows or moves the inside tracks (DF-D02) and
is never wrapped (DF-D25): the tracks are placed as without it, the Auto name
rule's taxon and strain lines (not a wrap) are each drawn on one line, and a
definition that does not fit is drawn over the tracks, at any definition font
size and for a typed Species. An explicit ``center_reserved_radius`` keeps the
inside tracks out of the centre (DF-D05); a failure inside it names the
radius, the room kept inside, and the tracks to turn off or narrow.
"""

from __future__ import annotations

import math
import re
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from fontTools.pens.recordingPen import RecordingPen

from gbdraw.api import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
)
from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram
from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.core.fonts import load_font, resolve_font_face
from gbdraw.core.text import (
    TEXT_ANCHOR_SHARE,
    _get_kerning_value,
    collapse_svg_white_space,
    dominant_baseline_shift_em,
    get_font_vertical_metrics,
)
from gbdraw.diagrams.circular.assemble import _definition_reserved_radius_px
from gbdraw.exceptions import ValidationError
from gbdraw.core.record_metadata import infer_record_source_metadata
from gbdraw.render.groups.circular.definition import CircularDefinitionStyle, DefinitionGroup, DefinitionNameLines
from gbdraw.web_support.error_adapter import serialize_web_error

REPO_ROOT = Path(__file__).resolve().parents[1]
HMMT = REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
SVG_NS = {"svg": "http://www.w3.org/2000/svg"}

SALMONELLA = "Salmonella enterica subsp. enterica serovar Typhimurium"
MYCOBACTERIUM = "Mycobacterium tuberculosis variant bovis BCG"
ESCHERICHIA = "Escherichia coli str. K-12 substr. MG1655"
HEPATOPLASMA = "Candidatus Hepatoplasma crinochetorum"
# A long name the Auto rule keeps whole on line 1 (a virus name); a Latin name is cut into its taxon and strain
# lines.
LONG_VIRUS = "Severe acute respiratory syndrome-related coronavirus"
MJENMV = "Marsupenaeus japonicus endogenous nimavirus"

# The Web defaults for a fresh Circular Generate (Multi-Record Canvas on).
WEB_DEFAULT_OVERRIDES = {
    "canvas.circular.track_type": "tuckin",
    "canvas.strandedness": True,
    "canvas.show_gc": True,
    "canvas.show_skew": True,
    "labels.circular.scope": "none",
}
WEB_DEFAULT_CLI_ARGS = [
    "--multi_record_canvas",
    "--track_type",
    "tuckin",
    "--gc",
    "--skew",
    "--separate_strands",
]


def _record_with_organism(organism: str) -> SeqRecord:
    record = next(SeqIO.parse(str(HMMT), "genbank"))
    for feature in record.features:
        if feature.type == "source":
            feature.qualifiers["organism"] = [organism]
    record.annotations["organism"] = organism
    if organism in (LONG_VIRUS, MJENMV):
        record.annotations["taxonomy"] = ["Viruses", "Riboviria"]
    return record


def _genbank_with_organism(tmp_path: Path, organism: str) -> Path:
    text = HMMT.read_text()
    path = tmp_path / "long_organism.gbk"
    path.write_text(text.replace('/organism="Homo sapiens"', f'/organism="{organism}"', 1))
    return path


def _web_default_grid(
    record: SeqRecord,
    *,
    config_overrides: dict | None = None,
    center_reserved_radius: float | None = None,
    species: str | None = None,
    circular_track_slots: list[str] | None = None,
):
    overrides = dict(WEB_DEFAULT_OVERRIDES)
    overrides.update(config_overrides or {})
    return build_circular_multi_diagram(
        [record],
        options=CircularDiagramOptions(
            config_overrides=overrides,
            output=CircularOutputOptions(legend="left"),
            tracks=CircularRequestTrackOptions(
                circular_track_slots=circular_track_slots,
                center_reserved_radius=center_reserved_radius,
            ),
            species=species,
        ),
        layout=CircularMultiRecordOptions(multi_record_positions=["#1@1"]),
    )


def _layout_reason(error: ValidationError) -> str | None:
    diagnostic = error.diagnostic or {}
    return diagnostic.get("reason") if diagnostic.get("code") == "TRACK_LAYOUT" else None


def _definition_lines(svg_text: str) -> list[str]:
    root = ET.fromstring(svg_text)
    group = root.find(".//svg:g[@data-gbdraw-role='record-definition']", SVG_NS)
    assert group is not None
    return ["".join(text.itertext()) for text in group.findall("svg:text", SVG_NS)]


def _track_paths(svg_text: str) -> dict[str, str]:
    """The drawn GC content and GC skew paths, by group id: their radii and widths."""

    root = ET.fromstring(svg_text)
    paths: dict[str, str] = {}
    for group in root.iter("{http://www.w3.org/2000/svg}g"):
        group_id = str(group.get("id") or "")
        if group_id.startswith(("gc_content", "skew", "gc_skew")):
            paths[group_id] = "|".join(str(path.get("d")) for path in group.iter("{http://www.w3.org/2000/svg}path"))
    assert paths
    return paths


@pytest.mark.circular
@pytest.mark.parametrize(
    "organism,drawn",
    [
        (SALMONELLA, SALMONELLA),
        (MYCOBACTERIUM, MYCOBACTERIUM),
        # The Auto name writes Candidatus as an italic "Ca." (OV-1720).
        (HEPATOPLASMA, "Ca. Hepatoplasma crinochetorum"),
    ],
)
def test_web_default_cli_generates_long_organism_names(
    organism: str, drawn: str, tmp_path: Path, gbdraw_runner
) -> None:
    source = _genbank_with_organism(tmp_path, organism)

    returncode, output, svg_path = gbdraw_runner.run(
        "circular",
        [source],
        "long_organism",
        tmp_path,
        extra_args=WEB_DEFAULT_CLI_ARGS,
    )

    assert returncode == 0, output
    lines = _definition_lines(svg_path.read_text())
    assert " ".join(lines).startswith(drawn)


@pytest.mark.circular
def test_a_name_the_floor_cannot_clear_keeps_the_inside_tracks() -> None:
    """DF-D36: on the Web default grid the GC rows cannot move outward, and narrowing them to the 30 px floor does
    not clear these names, so a name too long for the centre leaves the GC rows where a short name leaves them."""

    short = _track_paths(_web_default_grid(_record_with_organism("Homo sapiens")).tostring())
    for organism in (LONG_VIRUS, "Salmonellaentericasubspentericaserovartyphimurium"):
        assert _track_paths(_web_default_grid(_record_with_organism(organism)).tostring()) == short


def _short_name_tracks(config_overrides: dict | None = None) -> dict[str, str]:
    return _track_paths(
        _web_default_grid(_record_with_organism("Homo sapiens"), config_overrides=config_overrides).tostring()
    )


@pytest.mark.circular
@pytest.mark.parametrize(
    ("organism", "species", "config_overrides", "drawn"),
    [
        # The default definition font.
        (LONG_VIRUS, None, None, LONG_VIRUS),
        # A long virus name (the Auto rule keeps it whole on one taxon line).
        (MJENMV, None, None, MJENMV),
        # An explicit definition font size (W2 is moot).
        (LONG_VIRUS, None, {"objects.definition.circular.font_size": 17}, LONG_VIRUS),
        # A typed Species keeps its text on one line.
        ("Homo sapiens", "<i>Marsupenaeus japonicus</i> nimavirus", None, "Marsupenaeus japonicus nimavirus"),
    ],
    ids=["default-font", "long-virus-name", "explicit-font", "typed-species"],
)
def test_a_name_line_is_never_wrapped(
    organism: str, species: str | None, config_overrides: dict | None, drawn: str
) -> None:
    """DF-D25: a name line too long for the centre stays one line. With track narrowing off (DF-D36) it is drawn
    over tracks that stay where a short name leaves them; by default the inside tracks may narrow around it
    (``tests/test_circular_track_narrowing.py``), and the line is still one line."""

    no_narrowing = {**(config_overrides or {}), "canvas.circular.track_narrowing": False}
    for overrides in (config_overrides, no_narrowing):
        svg = _web_default_grid(
            _record_with_organism(organism), species=species, config_overrides=overrides
        ).tostring()
        lines = _definition_lines(svg)

        assert lines[0] == drawn
        assert len(lines) == 5  # name, mitochondrion, accession, length, GC
    assert _track_paths(svg) == _short_name_tracks(no_narrowing)


@pytest.mark.circular
def test_the_auto_taxon_and_strain_lines_are_not_a_wrap() -> None:
    """The Auto name rule's two lines stay; each is one line."""

    for organism, auto_lines in (
        (SALMONELLA, ["Salmonella enterica subsp. enterica", "serovar Typhimurium"]),
        (MYCOBACTERIUM, ["Mycobacterium tuberculosis", "variant bovis BCG"]),
        (ESCHERICHIA, ["Escherichia coli", "str. K-12 substr. MG1655"]),
    ):
        lines = _definition_lines(_web_default_grid(_record_with_organism(organism)).tostring())
        assert lines[:2] == auto_lines
        assert len(lines) == 6


@pytest.mark.circular
def test_a_typed_species_keeps_its_italics_on_one_line() -> None:
    typed = "<i>Marsupenaeus japonicus</i> nimavirus"
    svg = _web_default_grid(_record_with_organism("Homo sapiens"), species=typed).tostring()
    root = ET.fromstring(svg)
    group = root.find(".//svg:g[@data-gbdraw-role='record-definition']", SVG_NS)
    assert group is not None
    first = group.findall("svg:text", SVG_NS)[0]
    italic = ["".join(span.itertext()) for span in first.iter() if span.get("font-style") == "italic"]
    assert " ".join(italic).split() == ["Marsupenaeus", "japonicus"]


@pytest.mark.circular
def test_a_pinned_track_on_the_definition_leaves_the_name_on_one_line() -> None:
    slots = ["gc_skew:dinucleotide_skew@side=inside,r=0.5"]
    lines = _definition_lines(
        _web_default_grid(_record_with_organism(LONG_VIRUS), circular_track_slots=slots).tostring()
    )

    assert lines[0] == LONG_VIRUS

    # An explicit center_reserved_radius that leaves the pinned row no room fails on that radius.
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(
            _record_with_organism(LONG_VIRUS),
            circular_track_slots=slots,
            center_reserved_radius=250.0,
        )
    assert _layout_reason(caught.value) == "CENTER_RESERVED"


@pytest.mark.circular
def test_explicit_center_reserved_radius_keeps_one_line() -> None:
    lines = _definition_lines(
        _web_default_grid(
            _record_with_organism(LONG_VIRUS), center_reserved_radius=150.0
        ).tostring()
    )

    assert lines[0] == LONG_VIRUS

    # A radius larger than the room inside the tracks fails with the Web defaults.
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(_record_with_organism(LONG_VIRUS), center_reserved_radius=250.0)

    # The user's radius, not the species text, is named as the cause.
    assert _layout_reason(caught.value) == "CENTER_RESERVED"
    assert "because center_reserved_radius reserves 250.0px" in str(caught.value)
    assert "species" not in str(caught.value)


@pytest.mark.circular
def test_a_definition_that_does_not_fit_is_drawn_over_the_tracks() -> None:
    """No DEFINITION_RESERVED failure: a word too long for the centre draws over the tracks."""

    unbreakable = "Salmonellaentericasubspentericaserovartyphimurium"

    svg = _web_default_grid(_record_with_organism(SALMONELLA), species=unbreakable).tostring()

    assert _definition_lines(svg)[0] == unbreakable
    assert _track_paths(svg) == _short_name_tracks()


@pytest.mark.circular
def test_mers_cov_draws_its_name_on_one_line() -> None:
    record = next(SeqIO.parse(str(REPO_ROOT / "tests" / "test_inputs" / "MERS-CoV.gbk"), "genbank"))
    lines = _definition_lines(_web_default_grid(record).tostring())
    assert lines[0] == "Middle East respiratory syndrome-related coronavirus"


@pytest.mark.circular
def test_the_center_reserved_radius_message_names_the_tracks_kept_inside(tmp_path: Path) -> None:
    """OV-1268: the band starts above the reserved centre by the room kept for the tracks inside it."""

    error = _lambda_resolved_overlaps(tmp_path)
    assert _layout_reason(error) == "CENTER_RESERVED"
    message = str(error)
    assert "cannot fit inside between 189.7px and 386.1px" in message
    assert "because center_reserved_radius reserves 122.6px" in message
    assert "and 67.1px are kept for the tracks inside it ('gc_content', 'gc_skew')" in message


def _lambda_resolved_overlaps(tmp_path: Path, *extra: str) -> ValidationError:
    """Lambda with Web-like tracks inside an explicit centre (the definition alone no longer stops it)."""

    from gbdraw.circular import circular_main

    with pytest.raises(ValidationError) as caught:
        circular_main(
            [
                "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"),
                "--separate_strands", "--track_type", "tuckin", "-l", "left", "--resolve_overlaps",
                "--center_reserved_radius", "122.6",
                "-o", str(tmp_path / "out"), "-f", "svg", *extra,
            ]
        )
    return caught.value


@pytest.mark.circular
def test_the_center_reserved_radius_failure_offers_the_gc_tracks_kept_inside(tmp_path: Path) -> None:
    """OV-1274: the GC rows' kept room is a context key, and turning them off or narrowing them is a remedy."""

    error = _lambda_resolved_overlaps(tmp_path)
    assert _layout_reason(error) == "CENTER_RESERVED"
    assert (
        "Set a smaller center_reserved_radius, turn off or narrow the GC content and GC skew tracks "
        "(--no-gc, --no-skew, --gc_content_width, --gc_skew_width), or place tracks outside."
    ) in str(error)
    context = serialize_web_error(error, operation="generate", stage="render")["context"]
    assert context["keptInsidePx"] == 67
    assert context["keptInsideTracks"] == "gc_content,gc_skew"

    # Only the row still drawn is kept and named.
    error = _lambda_resolved_overlaps(tmp_path, "--no-gc")
    assert "turn off or narrow the GC skew track (--no-skew, --gc_skew_width)" in str(error)
    assert "--no-gc" not in str(error)
    context = serialize_web_error(error, operation="generate", stage="render")["context"]
    assert context["keptInsideTracks"] == "gc_skew"
    assert 33 <= context["keptInsidePx"] <= 34

    # With both rows off no GC row is named or offered. (This figure still fails, on the scale ticks: OV-1276,
    # whose window and remedy the next test pins.)
    error = _lambda_resolved_overlaps(tmp_path, "--no-gc", "--no-skew")
    assert "GC" not in str(error)
    assert "keptInsideTracks" not in serialize_web_error(error, operation="generate", stage="render")["context"]


@pytest.mark.circular
def test_the_room_kept_inside_names_every_row_and_offers_only_the_gc_rows(tmp_path: Path) -> None:
    """OV-1274 (review N-1): a Depth row shares the kept room; only the GC rows travel to the Web as named rows.

    tests/web/error-normalization.test.mjs words this context ("among them GC content and GC skew").
    """

    depth = tmp_path / "lambda.depth.tsv"
    depth.write_text("".join(f"NC_001416.1\t{position}\t{10 + position % 7}\n" for position in range(1, 48503, 50)))
    error = _lambda_resolved_overlaps(tmp_path, "--depth_track", str(depth))
    message = str(error)
    assert "90.8px are kept for the tracks inside it ('depth', 'gc_content', 'gc_skew')" in message
    assert "turn off or narrow the GC content and GC skew tracks" in message
    context = serialize_web_error(error, operation="generate", stage="render")["context"]
    assert (context["keptInsidePx"], context["keptInsideTracks"]) == (91, "gc_content,gc_skew")
    assert "slotTrack" not in context


@pytest.mark.circular
def test_the_center_reserved_radius_failure_starts_its_band_at_the_reserve(tmp_path: Path) -> None:
    """OV-1276: a window with nothing kept inside it starts at the reserved centre, in the message and the context;
    the scale ticks that do not fit offer hiding the scale."""

    error = _lambda_resolved_overlaps(tmp_path, "--no-gc", "--no-skew")
    assert _layout_reason(error) == "CENTER_RESERVED"
    assert (
        "Circular track slot 'ticks' cannot fit inside between 122.6px and 386.1px "
        "because center_reserved_radius reserves 122.6px."
    ) in str(error)
    assert "Set a smaller center_reserved_radius, hide the scale (--hide_scale; with an explicit track list, remove the ticks slot), or place tracks outside." in str(error)
    context = serialize_web_error(error, operation="generate", stage="render")["context"]
    assert (context["innerPx"], context["outerPx"], context["slotTrack"]) == (123, 386, "ticks")


@pytest.mark.circular
def test_the_scale_remedy_names_the_ticks_slot_of_an_explicit_track_list(tmp_path: Path) -> None:
    """OV-1276 (review B-1): with an explicit track list --hide_scale leaves an enabled ticks slot, so the remedy says so."""

    error = _lambda_resolved_overlaps(tmp_path, "--no-gc", "--no-skew", "--circular_track_order", "features,ticks")
    assert "'ticks' cannot fit inside" in str(error)
    assert "hide the scale (--hide_scale; with an explicit track list, remove the ticks slot)" in str(error)


def test_a_kept_scale_row_adds_no_scale_remedy() -> None:
    """OV-1276 (review N-1): the scale is offered only when it is the failed slot; a kept GC row is still offered."""

    from gbdraw.diagrams.circular.radial_layout import PlacementWindow, _InsideFitError, _KeptInsideSlot

    window = PlacementWindow(
        200,
        386,
        kept_inside_px=40,
        kept_inside_slots=(_KeptInsideSlot("ticks", "ticks"), _KeptInsideSlot("gc_content", "dinucleotide_content")),
    )
    message = str(_InsideFitError("slot", intent=None, placement_window=window))
    assert "turn off or narrow the GC content track (--no-gc, --gc_content_width)" in message
    assert "--hide_scale" not in message


@pytest.mark.circular
@pytest.mark.parametrize("multi_record", [False, True], ids=["one record", "multi-record canvas"])
def test_a_plot_title_with_no_name_to_draw_draws_no_title(multi_record: bool) -> None:
    """A record with no organism (GenBank writes ORGANISM ".") and no typed plot
    title has no plot-title line (OV-1831); asking for the title at the top
    draws no title instead of failing the composition (OV-1837)."""
    record = _record_with_organism(".")
    options = CircularDiagramOptions(output=CircularOutputOptions(plot_title_position="top"))
    drawing = (
        build_circular_multi_diagram([record], options=options)
        if multi_record
        else build_circular_diagram(record, options=options)
    )
    root = ET.fromstring(drawing.tostring())
    assert root.find(".//svg:g[@data-gbdraw-composition-role='title']", SVG_NS) is None



# One drawn text line: (x, y, font attributes, (text, font_style) runs).
_OutlineLine = tuple[float, float, dict[str, str], tuple[tuple[str, str], ...]]


def _outline_reach(lines: Iterable[_OutlineLine], dx: float, dy: float) -> float:
    """The farthest glyph outline point (any point of the glyph programs) of drawn text lines from the centre, the
    lines moved by (dx, dy). Measured from the font's outlines, independently of the ink boxes the reserve reads."""

    reach = 0.0
    for x, y, attrs, runs in lines:
        family, size = attrs["font-family"], float(attrs["font-size"])
        weight = attrs.get("font-weight", "normal")
        metrics = get_font_vertical_metrics(family, weight)
        assert metrics is not None
        shift = dominant_baseline_shift_em(attrs.get("dominant-baseline", "auto"), metrics)
        assert shift is not None
        points: list[tuple[float, float]] = []
        origin = 0.0
        for text, style in runs:
            face = resolve_font_face(family, weight, style)
            assert face is not None
            font = load_font(face)
            assert font is not None
            upm = float(font["head"].unitsPerEm)
            cmap = font.getBestCmap()
            glyphs = font.getGlyphSet()
            previous = None
            for char in collapse_svg_white_space(text):
                name = cmap.get(ord(char))
                if previous is not None:
                    origin += _get_kerning_value(font, previous, name) / upm
                if name is not None:
                    pen = RecordingPen()
                    glyphs[name].draw(pen)
                    points.extend((origin + px / upm, py / upm) for _op, args in pen.value for px, py in args)
                    origin += font["hmtx"][name][0] / upm
                previous = name
        start = x + dx - TEXT_ANCHOR_SHARE[attrs.get("text-anchor", "start")] * origin * size
        reach = max(
            [reach, *(math.hypot(start + px * size, y + dy + (shift - py) * size) for px, py in points)]
        )
    return reach


def _group_lines(definition: DefinitionGroup) -> list[_OutlineLine]:
    lines: list[_OutlineLine] = []
    for element in definition.get_group().elements:
        attrs = {key: str(value) for key, value in element.attribs.items()}
        style = attrs.get("font-style", "normal")
        runs = [(str(element.text), style)] if element.text else []
        runs += [(str(child.text), str(child.attribs.get("font-style", style))) for child in element.elements]
        lines.append((float(attrs["x"]), float(attrs["y"]), attrs, tuple(runs)))
    return lines


def _translation(transform: str) -> tuple[float, float]:
    pairs = re.findall(r"translate\(([-\d.e]+)[ ,]+([-\d.e]+)\)", transform)
    return sum(float(x) for x, _ in pairs), sum(float(y) for _, y in pairs)


def _drawn_definition_reach(svg_text: str) -> float:
    """The farthest outline point of the drawn centre definition from the record's centre (its axis group)."""

    root = ET.fromstring(svg_text)
    axis = root.find(".//svg:g[@id='Axis_0']", SVG_NS)
    group = root.find(".//svg:g[@data-gbdraw-role='record-definition']", SVG_NS)
    assert axis is not None and group is not None
    (cx, cy), (gx, gy) = _translation(axis.attrib["transform"]), _translation(group.attrib["transform"])
    lines: list[_OutlineLine] = []
    for text in group.findall("svg:text", SVG_NS):
        style = text.attrib.get("font-style", "normal")
        runs = [(text.text, style)] if text.text else []
        runs += [(span.text or "", span.attrib.get("font-style", style)) for span in text.findall("svg:tspan", SVG_NS)]
        lines.append((float(text.attrib["x"]), float(text.attrib["y"]), dict(text.attrib), tuple(runs)))
    return _outline_reach(lines, gx - cx, gy - cy)


def _innermost_band_inner_px(drawing) -> float:
    record = drawing._gbdraw_track_slot_geometry["records"][0]
    axis = float(record["axisRadiusPx"])
    return min(
        float(slot["radiusFactor"]) * axis - 0.5 * float(slot["widthPx"])
        for slot in record["slots"]
        if slot["side"] != "overlay"
    )


@pytest.mark.circular
@pytest.mark.parametrize(
    ("organism", "species"),
    [
        ("Homo sapiens", None),
        (HEPATOPLASMA, None),
        (ESCHERICHIA, None),
        ("Escherichia coli", "<i>Escherichia coli</i> str. K-12 substr. MG1655"),
        (MJENMV, None),
    ],
)
def test_the_definition_band_holds_the_farthest_ink_of_every_line(organism: str, species: str | None) -> None:
    """V-OD16, OV-1854: the definition band is each line's farthest painted ink corner, centred as drawn, plus its
    padding: never less than any glyph outline point (bold names, italic runs), and close to it."""

    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    record = _record_with_organism(organism)
    geometry = CircularCanvasConfigurator("test", CircularRenderProfile(cfg), "none", record).geometry()
    band = _definition_reserved_radius_px(
        record, geometry, species, None, cfg=cfg, plot_title=None, definition_profile="full"
    )
    definition = DefinitionGroup(
        record,
        geometry,
        style=CircularDefinitionStyle.of(cfg, "full"),
        names=DefinitionNameLines.of(infer_record_source_metadata(record), "full", species=species),
    )
    bounds = definition.local_bounds
    # place_definition_group_on_canvas centres the line boxes' union on the record's centre.
    reach = _outline_reach(
        _group_lines(definition), -0.5 * (bounds.min_x + bounds.max_x), -0.5 * (bounds.min_y + bounds.max_y)
    )
    held = band - max(8.0, 0.02 * float(geometry.radius))
    # For these names the box corner lies past the outline by under 3 px (a slanted italic run leaves its corner
    # empty). That is not a general bound: AP027131's name gives 3.1 px.
    assert reach - 1e-6 <= held <= reach + 3.0, (reach, held)


@pytest.mark.circular
@pytest.mark.parametrize("organism", [HEPATOPLASMA, ESCHERICHIA])
def test_the_definition_ink_stays_off_the_innermost_track(organism: str) -> None:
    """OV-1824: at the Web defaults these names clear the innermost track band: rule L narrows the tracks around the
    painted letters, so no glyph of the drawn definition enters the band and no notice is printed."""

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        drawing = _web_default_grid(_record_with_organism(organism))

    assert _drawn_definition_reach(drawing.tostring()) <= _innermost_band_inner_px(drawing)
    assert not [warning for warning in caught if "definition line" in str(warning.message)]
