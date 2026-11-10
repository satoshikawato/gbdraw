"""Several records drawn one after another on one circle (#200): the unit layout, the display mapping, ticks
numbered per record, GC per record and the count of records left out."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.diagram import build_circular_multi_diagram
from gbdraw.api.options import CircularDiagramOptions, CircularMultiRecordOptions
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, RecordInput
from gbdraw.exceptions import ValidationError
from gbdraw.layout.record_coordinates import CirclePlace, RecordDisplayTransform, SourceInterval
from gbdraw.layout.record_units import circle_units, record_units
from gbdraw.mode_profiles import ONE_CIRCLE_CONFIG_OVERRIDES


def _record(record_id: str, length: int, *, gc_every: int = 4) -> SeqRecord:
    bases = "".join("G" if index % gc_every == 0 else "A" for index in range(length))
    record = SeqRecord(Seq(bases), id=record_id, name=record_id, description=f"{record_id} test record")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["organism"] = "Testus example"
    record.features = [
        SeqFeature(FeatureLocation(0, length, strand=1), type="source"),
        SeqFeature(FeatureLocation(100, 400, strand=1), type="CDS", qualifiers={"product": [f"{record_id} gene"]}),
    ]
    return record


def test_units_lie_one_after_another_with_equal_gaps() -> None:
    records = [_record("a", 30_000), _record("b", 20_000), _record("c", 10_000)]
    layout = circle_units(record_units(records, None), gap_degrees=1.0, min_length_bp=1000)

    gap = layout.gap_bp
    assert gap == math.ceil(1.0 * 60_000 / (360.0 - 3.0))
    assert layout.display_length == 60_000 + 3 * gap
    offsets = [place.offset for _unit, place in layout.placed]
    assert offsets == [gap // 2, gap // 2 + 30_000 + gap, gap // 2 + 50_000 + 2 * gap]
    # The gaps are the asked angle round the circle, the last one across 12 o'clock included.
    assert 360.0 * gap / layout.display_length == pytest.approx(1.0, abs=0.01)
    assert {place.scale_length for _unit, place in layout.placed} == {30_000}


def test_records_under_the_shortest_drawn_length_are_counted_not_drawn() -> None:
    records = [_record("a", 30_000), _record("tiny", 400), _record("b", 5_000)]
    layout = circle_units(record_units(records, None), gap_degrees=1.0, min_length_bp=1000)

    assert [unit.record_id for unit, _place in layout.placed] == ["a", "b"]
    assert [unit.record_id for unit in layout.hidden] == ["tiny"]
    with pytest.raises(ValidationError, match="shorter than 1,000 bp"):
        circle_units(record_units([_record("tiny", 400)], None), gap_degrees=1.0, min_length_bp=1000)
    with pytest.raises(ValidationError, match="no room"):
        circle_units(record_units(records, None), gap_degrees=180.0, min_length_bp=0)


def test_a_record_on_the_circle_maps_its_bases_after_its_offset_reversed_when_drawn_reverse() -> None:
    place = CirclePlace(offset=500, length=10_000, scale_length=2_000)
    forward = RecordDisplayTransform(2_000, 1, 1, None, False).with_circle(place, 2_000)
    reverse = RecordDisplayTransform(2_000, 2_000, -1, None, False).with_circle(place, 2_000)

    assert forward.projects
    assert forward.source_boundary_to_display_offset(0) == 500
    assert forward.source_boundary_to_display_offset(2_000) == 2_500
    assert forward.source_base_to_display_index(1) == 500
    assert forward.display_index_to_source_base(500) == 1
    # Reverse: the record's end comes first on its arc.
    assert reverse.source_boundary_to_display_offset(2_000) == 500
    assert reverse.source_boundary_to_display_offset(0) == 2_500
    (fragment,) = reverse.project_interval(SourceInterval(100, 400, 1))
    assert (fragment.display_start, fragment.display_end) == (500 + 1_600, 500 + 1_900)
    assert (fragment.local_start, fragment.local_end) == (1_600, 1_900)
    assert fragment.orientation == -1
    with pytest.raises(ValidationError, match="must end on it"):
        RecordDisplayTransform(2_000, 1, 1, None, False).with_circle(CirclePlace(9_000, 10_000, 2_000), 2_000)


def _one_circle_svg(records: list[SeqRecord], **layout: object) -> str:
    drawing = build_circular_multi_diagram(
        records,
        options=CircularDiagramOptions(),
        layout=CircularMultiRecordOptions(record_layout="one_circle", **layout),  # type: ignore[arg-type]
    )
    return drawing.tostring()


def _slot_groups(svg: str, renderer: str) -> list[ET.Element]:
    """Each record's group of the slot ``renderer`` draws, in record order: the first record's group keeps the slot's
    id, the N-th record's ends in ``__record_N``."""

    groups = [group for group in ET.fromstring(svg).iter() if group.get("data-gbdraw-slot-renderer") == renderer]

    def position(group: ET.Element) -> int:
        found = re.search(r"__record_(\d+)$", group.get("id", ""))
        return int(found.group(1)) if found else 1

    return sorted(groups, key=position)


def _gc_parts(group: ET.Element) -> tuple[list[ET.Element], list[ET.Element]]:
    """A GC content group's percent background arcs (the configured background colour) and its axis groups (the
    groups that draw the axis ticks or the range label)."""

    from gbdraw.config.toml import load_config_toml

    background = str(load_config_toml("gbdraw.data", "config.toml")["objects"]["gc_content"]["percent_background_color"])
    backgrounds = [path for path in group if path.tag.endswith("path") and path.get("stroke") == background]
    axes = [
        child for child in group
        if child.tag.endswith("g") and any(item.tag.endswith(("line", "textPath")) for item in child.iter())
    ]
    return backgrounds, axes


def _gc_counts(svg: str) -> tuple[int, int]:
    """How many percent background arcs and axis groups the circle's GC content groups draw."""

    parts = [_gc_parts(group) for group in _slot_groups(svg, "dinucleotide_content")]
    return sum(len(backgrounds) for backgrounds, _ in parts), sum(len(axes) for _, axes in parts)


def test_one_circle_draws_each_record_on_its_own_arc_with_one_definition() -> None:
    records = [_record("alpha", 120_000), _record("beta", 80_000), _record("tiny", 500, gc_every=2)]
    svg = _one_circle_svg(records)

    # One feature group per drawn record, each with its own record identity.
    record_ids = re.findall(r'data-gbdraw-record-id="([^"]+)"', svg)
    assert "alpha" in record_ids and "beta" in record_ids and "tiny" not in record_ids
    # The centre definition counts the records it draws and those it leaves out, with their bases and GC.
    assert "2 records (+1 under 1,000 bp not drawn)" in svg
    assert "200,000 bp" in svg
    # Every id stays unique however many records draw the same groups.
    ids = re.findall(r'(?<![-\w])id="([^"]+)"', svg)
    assert len(ids) == len(set(ids))


def test_one_circle_numbers_ticks_from_each_record_start_at_one_interval() -> None:
    records = [_record("alpha", 120_000), _record("beta", 80_000)]
    svg = _one_circle_svg(records)

    labels = re.findall(r">(\d+) kbp<", svg)
    # Both records count from their own start: beta's labels restart below alpha's last one.
    values = [int(value) for value in labels]
    assert values and max(values) < 120 and values != sorted(values)


def test_one_circle_draws_gc_in_percent_without_skew_unless_asked() -> None:
    records = [_record("alpha", 120_000), _record("beta", 80_000)]
    assert dict(ONE_CIRCLE_CONFIG_OVERRIDES) == {"canvas.show_skew": False, "objects.gc_content.mode": "percent"}
    svg = _one_circle_svg(records)
    assert _slot_groups(svg, "dinucleotide_skew") == []
    # One percent axis for the circle, a background arc per record.
    assert _gc_counts(svg) == (2, 1)

    drawing = build_circular_multi_diagram(
        records,
        options=CircularDiagramOptions(config_overrides={"canvas.show_skew": True}),
        layout=CircularMultiRecordOptions(record_layout="one_circle"),
    )
    assert "skew" in drawing.tostring()


def test_a_record_shorter_than_the_gc_window_draws_one_constant_window() -> None:
    # 5 kbp next to 200 kbp: the circle's window (its bases' Auto window) is longer than the short record.
    records = [_record("long", 200_000), _record("short", 5_000, gc_every=2)]
    svg = _one_circle_svg(records)
    assert len(_slot_groups(svg, "dinucleotide_content")) == 2, "the short record draws its own GC content group"


def test_the_request_resolves_one_circle_defaults_and_keeps_the_layout() -> None:
    request = CircularDiagramRequest(
        records=[RecordInput(InMemoryRecordSource(_record(name, 5_000))) for name in ("a", "b")],
        layout=CircularMultiRecordOptions(record_layout="one_circle", one_circle_gap_degrees=2.0),
    )
    assert request.grouping == "grid"
    assert request.layout is not None and request.layout.record_layout == "one_circle"
    overrides = dict(request.options.config_overrides or {})
    assert overrides["canvas.show_skew"] is False
    assert overrides["objects.gc_content.mode"] == "percent"
    with pytest.raises(ValidationError, match="record_layout"):
        CircularMultiRecordOptions(record_layout="ring")  # type: ignore[arg-type]


def test_the_record_band_names_each_drawn_record_where_the_name_fits() -> None:
    from gbdraw.diagrams.circular.one_circle import record_band_names

    named = _record("chr4", 50_000)
    named.features[0].qualifiers["chromosome"] = ["IV"]
    assert record_band_names(named, 4) == ("IV", "4")
    assert record_band_names(_record("contig", 5_000), 7) == ("7",)

    svg = _one_circle_svg([named, _record("b", 30_000), _record("tiny", 500)])
    band = svg[svg.index('id="record_band_arcs"'):]
    arcs = re.findall(r'<g data-gbdraw-record-id="([^"]+)" data-gbdraw-record-key="[^"]+" id="record_band_\d+">', band)
    assert arcs == ["chr4", "b"]
    assert ">IV</text>" in band and ">2</text>" in band

    without = _one_circle_svg([named, _record("b", 30_000)], one_circle_record_band=False)
    assert "record_band_arcs" not in without



@pytest.mark.parametrize(
    ("words", "shown"),
    [([], True), (["--one_circle"], False), (["--one_circle", "--skew"], True), (["--no-skew"], False)],
)
def test_the_cli_draws_skew_as_typed_else_by_the_layout_default(words: list[str], shown: bool) -> None:
    from gbdraw.circular import _get_args, _shown_skew
    from gbdraw.mode_profiles import CIRCULAR_MODE_PROFILE

    args = _get_args(["--gbk", "a.gb", *words])
    # The parser keeps the profile's default; only a typed toggle notes itself.
    assert args.show_skew is (CIRCULAR_MODE_PROFILE.show_skew if not {"--skew", "--no-skew"} & set(words) else shown)
    assert _shown_skew(args) is shown


@pytest.mark.parametrize("form", ["depth_table", "depth_tables", "depth_files"])
def test_a_depth_track_is_refused_on_one_circle_in_every_input_form(form: str, tmp_path) -> None:
    # OV-3117: only the singular forms were refused; a table per record reached the first record's depth builder.
    import pandas as pd

    records = [_record("alpha", 30_000), _record("beta", 20_000)]
    tables = [
        pd.DataFrame({"reference_name": record.id, "position": list(range(1, len(record.seq) + 1, 100)), "depth": 10})
        for record in records
    ]
    files = []
    for record, table in zip(records, tables, strict=True):
        path = tmp_path / f"{record.id}.tsv"
        table.to_csv(path, sep="\t", header=False, index=False)
        files.append(str(path))
    given = {"depth_table": tables[0], "depth_tables": tables, "depth_files": files}[form]
    with pytest.raises(ValidationError, match="Depth tracks are not drawn") as refused:
        build_circular_multi_diagram(
            records,
            options=CircularDiagramOptions(**{form: given}),
            layout=CircularMultiRecordOptions(record_layout="one_circle"),
        )
    assert refused.value.diagnostic == {"code": "TRACK_INVALID", "reason": "ONE_CIRCLE_DEPTH"}


def test_the_cli_refuses_one_circle_with_the_grid_canvas() -> None:
    from gbdraw.circular import _get_args, _one_circle_layout_fields

    with pytest.raises(ValidationError, match="choose one"):
        _one_circle_layout_fields(_get_args(["--gbk", "a.gb", "--one_circle", "--multi_record_canvas"]))


def test_a_tick_label_that_would_run_past_its_record_is_left_out_and_its_mark_kept() -> None:
    # alpha ends 1 kbp after its 100 kbp tick: that label would reach over the gap into beta.
    svg = _one_circle_svg([_record("alpha", 101_000), _record("beta", 80_000)])

    labels = re.findall(r">(\d+) kbp<", svg)
    assert "100" not in labels and "90" in labels
    # The tick keeps its mark: alpha marks every 10 kbp from 10 to 100 kbp.
    alpha_ticks = _slot_groups(svg, "ticks")[0]
    assert sum(1 for item in alpha_ticks.iter() if item.get("stroke") == "gray") == 10


@pytest.mark.parametrize("cropped_first", [True, False])
def test_a_cropped_record_draws_over_its_cropped_bases_on_one_circle(cropped_first: bool, tmp_path) -> None:
    # A crop's transform keeps the complete source length: read as the drawn length, the cropped record's arc ran
    # past its own into the gap and its neighbour, and a cropped record placed last was refused.
    from Bio import SeqIO

    from gbdraw.circular import circular_main

    SeqIO.write([_record("alpha", 120_000)], tmp_path / "alpha.gb", "genbank")
    SeqIO.write([_record("beta", 80_000)], tmp_path / "beta.gb", "genbank")
    rows = [f"{tmp_path / 'beta.gb'}\t20001-70000", f"{tmp_path / 'alpha.gb'}\t"]
    table = tmp_path / "records.tsv"
    table.write_text("\n".join(["gbk\tregion", *(rows if cropped_first else rows[::-1])]) + "\n", encoding="utf-8")
    circular_main(["--records_table", str(table), "--one_circle", "-o", str(tmp_path / "one"), "-f", "svg"])
    svg = (tmp_path / "one.svg").read_text(encoding="utf-8")

    gap = math.ceil(170_000 / 358.0)
    length = 170_000 + 2 * gap
    start = gap // 2 if cropped_first else gap // 2 + 120_000 + gap
    (background,), _axes = _gc_parts(_slot_groups(svg, "dinucleotide_content")[0 if cropped_first else 1])
    x, y = (float(value) for value in background.get("d", "").split()[-1].split(","))
    end_degrees = math.degrees(math.atan2(x, -y)) % 360.0
    assert end_degrees == pytest.approx(360.0 * (start + 50_000) / length, abs=0.01)


def test_the_gc_axis_draws_no_full_baseline_circle_over_the_gaps() -> None:
    svg = _one_circle_svg([_record("alpha", 120_000), _record("beta", 80_000)])

    (axis,) = _gc_parts(_slot_groups(svg, "dinucleotide_content")[0])[1]
    # Each record's baseline is its own arc (the percent border); a full circle would cross the gaps.
    assert not [item for item in axis.iter() if item.tag.endswith("circle")]


@pytest.mark.parametrize(
    ("words", "percent"),
    [([], True), (["--gc_content_mode", "deviation"], False), (["--gc_content_mode", "percent"], True)],
)
def test_the_cli_draws_one_circle_gc_in_percent_unless_typed(words: list[str], percent: bool, tmp_path) -> None:
    # The CLI writes its whole config itself, so the one-circle profile values must reach it there too: the
    # first PR figures drew GC in deviation form from `gbdraw circular --one_circle`.
    from Bio import SeqIO

    from gbdraw.circular import circular_main

    source = tmp_path / "two.gb"
    SeqIO.write([_record("alpha", 120_000), _record("beta", 80_000)], source, "genbank")
    out = tmp_path / "one"
    circular_main(["--gbk", str(source), "--one_circle", *words, "-o", str(out), "-f", "svg", "--save_session"])
    svg = (tmp_path / "one.svg").read_text(encoding="utf-8")

    assert _slot_groups(svg, "dinucleotide_skew") == []
    assert _gc_counts(svg) == ((2, 1) if percent else (0, 0))
    # The saved Session carries every one-circle profile value the user left unset, so a replay draws the same.
    from tests.utils.session_files import read_session_payload

    session = read_session_payload(tmp_path / "one.gbdraw")
    config = session["drawings"][0]["renderRequest"]["diagramOptions"]["config"]
    for path, value in ONE_CIRCLE_CONFIG_OVERRIDES.items():
        leaf: object = config
        for part in path.split("."):
            assert isinstance(leaf, dict)
            leaf = leaf[part]
        assert leaf == (value if not words else {"canvas.show_skew": False, "objects.gc_content.mode": words[1]}[path])


@pytest.mark.parametrize("style", ["heatmap", "histogram"])
def test_a_feature_density_track_draws_each_record_on_its_own_arc_on_one_shared_scale(style: str) -> None:
    from gbdraw.analysis.feature_density import FeatureDensityTrack, FeatureQuery

    records = [_record("alpha", 120_000), _record("beta", 80_000)]
    drawing = build_circular_multi_diagram(
        records,
        options=CircularDiagramOptions(feature_density_tracks=(FeatureDensityTrack("cds", FeatureQuery("CDS"), style=style),)),
        layout=CircularMultiRecordOptions(record_layout="one_circle"),
    )
    svg = drawing.tostring()
    groups = re.findall(r'<g[^>]*data-gbdraw-slot-id="density_cds"[^>]*>(.*?)</g>', svg)
    assert len(groups) == 2
    # Every point the band paints lies on its own record's arc: alpha in the first 120/200 of the circle, beta after
    # it, neither across the other's arc or the gaps (the band read each record's own length as the circle's).
    for index, group in enumerate(groups):
        angles = [
            math.degrees(math.atan2(float(x), -float(y))) % 360
            for x, y in re.findall(r"(-?[\d.]+),(-?[\d.]+)", re.sub(r"A\s*[\d.]+,[\d.]+", "A", group))
            if math.hypot(float(x), float(y)) > 1
        ]
        assert angles
        inside = [(0.0, 216.0), (216.0, 360.0)][index]
        assert all(inside[0] - 0.5 <= angle <= inside[1] + 0.5 or angle > 359.5 for angle in angles), (index, min(angles), max(angles))


@pytest.mark.parametrize("reason", ["ONE_CIRCLE_CONSERVATION", "ONE_CIRCLE_REGIONS"])
def test_each_one_circle_refusal_names_a_reason_the_web_words_as_the_command(reason: str) -> None:
    from pathlib import Path

    from gbdraw.api.assembly_inputs import CircularAssemblyInputs
    from gbdraw.api.diagram import _refuse_one_circle_tables
    from gbdraw.api.config import apply_config_overrides
    from gbdraw.diagrams.circular.assemble import _refuse_one_circle_inputs

    with pytest.raises(ValidationError) as refused:
        if reason == "ONE_CIRCLE_CONSERVATION":
            _refuse_one_circle_tables(CircularAssemblyInputs(cfg=apply_config_overrides(None, None), conservation_blast_files=["ring.tsv"]))
        else:
            _refuse_one_circle_inputs({"annotations": object()})
    assert refused.value.diagnostic == {"code": "TRACK_INVALID", "reason": reason}
    # The Web shows the reason's words (utils/error-normalization.js), the command's message.
    words = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "utils" / "error-normalization.js"
    assert f"  {reason}: '{refused.value}'," in words.read_text(encoding="utf-8")


def test_records_left_off_are_logged_and_still_counted_when_the_definition_is_hidden(caplog) -> None:
    records = [_record("alpha", 120_000), _record("beta", 80_000), _record("tiny", 500, gc_every=2)]
    with caplog.at_level("INFO", logger="gbdraw.api.diagram"):
        drawing = build_circular_multi_diagram(
            records,
            options=CircularDiagramOptions(config_overrides={"objects.definition.circular.show": False}),
            layout=CircularMultiRecordOptions(record_layout="one_circle"),
        )
    assert "1 of 3 records are under 1,000 bp and not drawn" in caplog.text
    svg = drawing.tostring()
    # Only the count is drawn where the definition would be: no name, bases or GC line.
    assert "2 records (+1 under 1,000 bp not drawn)" in svg
    assert "200,000 bp" not in svg and "Testus" not in svg
    # Nothing left out: a hidden definition draws nothing.
    hidden = build_circular_multi_diagram(
        records[:2],
        options=CircularDiagramOptions(config_overrides={"objects.definition.circular.show": False}),
        layout=CircularMultiRecordOptions(record_layout="one_circle"),
    ).tostring()
    assert "records" not in re.sub(r"<[^>]+>", "", hidden).split()


def test_the_centre_names_only_what_every_drawn_record_shares() -> None:
    other = _record("beta", 80_000)
    other.annotations["organism"] = "Otherus example"
    two_species = re.sub(r"<[^>]+>", " ", _one_circle_svg([_record("alpha", 120_000), other]))
    assert "Testus" not in two_species and "Otherus" not in two_species
    assert "2 records" in two_species
    one_species = re.sub(r"<[^>]+>", " ", _one_circle_svg([_record("alpha", 120_000), _record("beta", 80_000)]))
    assert "Testus" in one_species


def _range_label_turns(svg: str) -> list[tuple[float, float]]:
    """Where each "[0–100%]" label's path starts and ends, in turns clockwise from 12 o'clock."""

    found = []
    for href in re.findall(r'<textPath[^>]*href="#([^"]+)"[^>]*>\[0–100%\]</textPath>', svg):
        path = re.search(rf'<path d="([^"]+)"[^>]*id="{re.escape(href)}"', svg) or re.search(
            rf'<path[^>]*id="{re.escape(href)}"[^>]*d="([^"]+)"', svg
        )
        assert path is not None, href
        numbers = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", path.group(1))]
        (x0, y0), (x1, y1) = numbers[:2], numbers[-2:]
        found.append(tuple(math.atan2(x, -y) / (2.0 * math.pi) % 1.0 for x, y in ((x0, y0), (x1, y1))))
    return found


def test_the_gc_range_label_lies_on_the_first_record_that_holds_it_never_over_a_gap() -> None:
    # D-48: one "[0–100%]" label, on the first record clockwise from 12 o'clock whose arc holds it.
    records = [_record("short", 3_000), _record("long", 400_000)]
    units = circle_units(record_units(records, None), gap_degrees=1.0, min_length_bp=1000).placed
    arcs = [(place.offset / place.length, (place.offset + unit.length) / place.length) for unit, place in units]
    svg = _one_circle_svg(records)
    narrow = _range_label_turns(svg)
    assert len(narrow) == 1
    # The first record (a sliver) cannot hold it: it starts on the second record and ends inside it.
    assert arcs[1][0] < narrow[0][0] < narrow[0][1] <= arcs[1][1] + 1e-6
    # That record draws it after its own band, so no band paints over it.
    second = _slot_groups(svg, "dinucleotide_content")[1]
    (background,), (axis,) = _gc_parts(second)
    drawn = list(second)
    assert drawn.index(background) < drawn.index(axis)
    assert "[0–100%]" in "".join(item.text or "" for item in axis.iter())
    wide = _range_label_turns(_one_circle_svg([_record("long", 400_000), _record("short", 3_000)]))
    assert len(wide) == 1 and wide[0][0] < 0.25
    # No record holds it: no label at all.
    many = [_record(f"r{index}", 4_000) for index in range(60)]
    assert _range_label_turns(_one_circle_svg(many)) == []


def test_a_record_band_name_that_does_not_fit_gives_way_to_a_shorter_one() -> None:
    # OV-3182 (D-51): the full name, else the organelle's standard abbreviation, else the record's number.
    from gbdraw.diagrams.circular.one_circle import record_band_names

    chromosome = _record("chrI", 200_000)
    chromosome.features[0].qualifiers["chromosome"] = ["I"]
    mitochondrion = _record("mito", 3_000)
    mitochondrion.features[0].qualifiers["organelle"] = ["mitochondrion"]
    plasmid = _record("plasmid", 3_000)
    plasmid.features[0].qualifiers["plasmid"] = ["pVeryLongPlasmidName"]
    chloroplast = _record("plastid", 3_000)
    chloroplast.features[0].qualifiers["organelle"] = ["plastid:chloroplast"]
    assert record_band_names(mitochondrion, 2) == ("Mitochondrion", "mt", "2")
    assert record_band_names(chloroplast, 4) == ("Plastid:chloroplast", "cp", "4")
    assert record_band_names(plasmid, 3) == ("pVeryLongPlasmidName", "3")

    svg = _one_circle_svg([chromosome, mitochondrion, plasmid, chloroplast])
    arcs = svg[svg.index('id="record_band_arcs"'):].split('id="record_band_')[2:6]
    assert [re.findall(r">([^<>]+)</text>", arc)[0] for arc in arcs] == ["I", "mt", "3", "cp"]


def test_neighbouring_band_names_keep_apart_where_the_gap_is_narrow() -> None:
    # OV-3183: names that just fit their arcs touched across a 0 degree gap ("12" "13" read as "1213").
    from gbdraw.core.text import drawn_text_box
    from gbdraw.render.drawers.circular.record_band import RecordBandArc, record_band_group

    family, size, radius = "Liberation Sans", 14.0, 400.0
    width, _ = drawn_text_box("13", family, size)
    span = math.degrees((width + 0.5) / radius)

    def drawn(gap_deg: float) -> list[str]:
        arcs = tuple(
            RecordBandArc(f"k{index}", f"r{index}", (str(12 + index),), -90.0 + index * (span + gap_deg), span)
            for index in range(2)
        )
        group = record_band_group(arcs, inner_px=radius - 10.0, outer_px=radius + 10.0, font_family=family, font_size=size)
        return re.findall(r">([^<>]+)</text>", group.tostring())

    assert drawn(0.0) == []
    # Arcs that meet to within rounding read as a 0 degree gap, not as a whole turn.
    assert drawn(-1e-9) == []
    # A gap as wide as the clearance needs no room from the arcs: both names fit.
    assert drawn(math.degrees(2.0 / radius)) == ["12", "13"]


def _range_label_path(svg: str) -> tuple[str, str]:
    """The path data and the ``textPath`` attributes of the one "[0–100%]" label."""

    (href, attributes), = [
        (match.group(2), match.group(1))
        for match in re.finditer(r'<textPath([^>]*href="#([^"]+)"[^>]*)>\[0–100%\]</textPath>', svg)
    ]
    path = re.search(rf'<path d="([^"]+)"[^>]*id="{re.escape(href)}"', svg)
    assert path is not None
    return path.group(1), attributes


def test_the_gc_range_label_reads_upright_on_the_lower_half() -> None:
    # OV-3181 (D-50): where no record before 3 o'clock holds it, the label lies on the lower half; it runs
    # anticlockwise there, from where its text ends back towards 12 o'clock, as the record band names turn there.
    records = [_record(f"r{index + 1}", 3_000) for index in range(40)] + [_record("r41", 300_000)]
    units = circle_units(record_units(records, None), gap_degrees=1.0, min_length_bp=1000).placed
    unit, place = units[-1]
    record_start = place.offset / place.length
    svg = _one_circle_svg(records)
    data, attributes = _range_label_path(svg)
    ((text_end, path_end),) = _range_label_turns(svg)
    assert record_start < text_end < 0.5 and text_end > 0.25
    assert path_end == pytest.approx(0.0, abs=1e-4)
    assert {flags.split()[-1] for flags in re.findall(r"A [^A]*? 0 (\d \d)", data)} == {"0"}
    assert "text-anchor" not in attributes and "startOffset" not in attributes


def test_a_range_label_turned_on_the_lower_half_keeps_its_letters_in_the_same_band() -> None:
    from gbdraw.layout.circular_depth_axis import CircularRangeLabel

    from gbdraw.core.text import calculate_svg_bbox_dimensions

    upper = CircularRangeLabel("[0–100%]", 10.0, 100.0, 109.5, 0.1)
    assert upper.text_path("Liberation Sans") == (100.0, -90.0 + math.degrees(0.1), upper.path_span_deg)
    lower = upper._replace(start_rad=math.pi, path_end_rad=1.5 * math.pi)
    radius, start_deg, span_deg = lower.text_path("Liberation Sans")
    # The letters span the band from 1.2 to 0.2 lines inside the outer edge either way: the baseline mirrored.
    assert radius == pytest.approx(2 * 109.5 - 1.4 * 10.0 - 100.0)
    # From where the text ends, anticlockwise to 12 o'clock: the text ends at the label's start (6 o'clock).
    width, _ = calculate_svg_bbox_dimensions("[0–100%]", "Liberation Sans", 10.0)
    text_end_deg = math.degrees(math.pi + width / radius)
    assert (start_deg, span_deg) == pytest.approx((-90.0 + text_end_deg, -text_end_deg))


def test_neighbouring_record_band_arcs_differ_across_12_oclock_too() -> None:
    from gbdraw.render.drawers.circular.record_band import record_band_fill

    for count in range(2, 8):
        fills = [record_band_fill(index, count) for index in range(count)]
        assert all(fills[index] != fills[index - 1] for index in range(count)), count

