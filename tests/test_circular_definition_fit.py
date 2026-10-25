"""Circular center-definition fit (PV-08, D-24 / PD-OI-078).

Only a placement that fails because of the center definition band wraps the
species line at word boundaries and places the tracks again. An explicit
``center_reserved_radius`` disables the wrap. ``definition_font_size`` 18 is
the default, so 18 counts as not explicit and may wrap; any other value is
explicit and never wraps. Earlier successful outputs keep a single species
line, and a failure that the wrap cannot fix names the definition band as the
cause.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
)
from gbdraw.api.diagram import build_circular_multi_diagram
from gbdraw.core.text import parse_mixed_content_text
from gbdraw.exceptions import ValidationError
from gbdraw.render.groups.circular.definition import wrap_definition_line_parts
from gbdraw.web_support.error_adapter import serialize_web_error

REPO_ROOT = Path(__file__).resolve().parents[1]
HMMT = REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
SVG_NS = {"svg": "http://www.w3.org/2000/svg"}

SALMONELLA = "Salmonella enterica subsp. enterica serovar Typhimurium"
MYCOBACTERIUM = "Mycobacterium tuberculosis variant bovis BCG"
ESCHERICHIA = "Escherichia coli str. K-12 substr. MG1655"
HEPATOPLASMA = "Candidatus Hepatoplasma crinochetorum"
# A long name the Auto rule keeps whole on line 1 (a virus name), for the
# wrap of one long species line; a Latin name is cut into its taxon and strain lines.
LONG_VIRUS = "Severe acute respiratory syndrome-related coronavirus"

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
    if organism == LONG_VIRUS:
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
def test_species_line_wraps_at_word_boundaries_only_when_placement_fails() -> None:
    wrapped = _definition_lines(_web_default_grid(_record_with_organism(LONG_VIRUS)).tostring())
    species_lines = wrapped[: len(wrapped) - 4]  # mitochondrion, accession, length, GC follow

    assert len(species_lines) >= 2
    assert " ".join(species_lines) == LONG_VIRUS
    assert all(line == line.strip() and line for line in species_lines)

    # Names whose taxon and strain lines fit are not wrapped.
    for organism, auto_lines in (
        (SALMONELLA, ["Salmonella enterica subsp. enterica", "serovar Typhimurium"]),
        (MYCOBACTERIUM, ["Mycobacterium tuberculosis", "variant bovis BCG"]),
        (ESCHERICHIA, ["Escherichia coli", "str. K-12 substr. MG1655"]),
    ):
        lines = _definition_lines(_web_default_grid(_record_with_organism(organism)).tostring())
        assert lines[:2] == auto_lines


@pytest.mark.circular
def test_pinned_track_on_the_definition_wraps_the_species_line() -> None:
    # A pinned row that overlaps only the one-line definition fails because of
    # the definition band, so the species line wraps (TK-02).
    slots = ["gc_skew:dinucleotide_skew@side=inside,r=0.5"]
    wrapped = _definition_lines(
        _web_default_grid(_record_with_organism(LONG_VIRUS), circular_track_slots=slots).tostring()
    )
    species_lines = wrapped[: len(wrapped) - 4]

    assert len(species_lines) >= 2
    assert " ".join(species_lines) == LONG_VIRUS

    # An explicit center_reserved_radius is not wrapped away.
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(
            _record_with_organism(LONG_VIRUS),
            circular_track_slots=slots,
            center_reserved_radius=250.0,
        )
    assert _layout_reason(caught.value) == "CENTER_RESERVED"


@pytest.mark.circular
def test_wrap_keeps_the_text_and_its_italic_runs() -> None:
    parts = parse_mixed_content_text(
        "<i>Salmonella enterica</i> subsp. <i>enterica</i> serovar Typhimurium"
    )

    def text(lines):
        return [[(part["text"], part["italic"]) for part in line] for line in lines]

    assert text(wrap_definition_line_parts(parts, 1, len)) == [
        [
            ("Salmonella enterica", True),
            (" subsp. ", False),
            ("enterica", True),
            (" serovar Typhimurium", False),
        ]
    ]
    assert text(wrap_definition_line_parts(parts, 2, len)) == [
        [("Salmonella enterica", True), (" subsp.", False)],
        [("enterica", True), (" serovar Typhimurium", False)],
    ]
    assert len(wrap_definition_line_parts(parts, 99, len)) == 6


@pytest.mark.circular
def test_explicit_definition_font_size_disables_the_wrap() -> None:
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(
            _record_with_organism(LONG_VIRUS),
            config_overrides={"objects.definition.circular.font_size": 17},
        )

    assert _layout_reason(caught.value) == "DEFINITION_RESERVED"
    assert "center definition text reserves" in str(caught.value)


@pytest.mark.circular
def test_explicit_center_reserved_radius_disables_the_wrap() -> None:
    lines = _definition_lines(
        _web_default_grid(
            _record_with_organism(LONG_VIRUS), center_reserved_radius=150.0
        ).tostring()
    )

    assert lines[0] == LONG_VIRUS

    # A radius larger than the one-line definition reservation fails
    # with the Web defaults; an explicit radius must not be replaced by a
    # wrapped definition reservation.
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(_record_with_organism(LONG_VIRUS), center_reserved_radius=250.0)

    # The user's radius, not the species text, is named as the cause.
    assert _layout_reason(caught.value) == "CENTER_RESERVED"
    assert "because center_reserved_radius reserves 250.0px" in str(caught.value)
    assert "species" not in str(caught.value)


@pytest.mark.circular
def test_definition_that_cannot_be_wrapped_reports_the_definition_band() -> None:
    unbreakable = "Salmonellaentericasubspentericaserovartyphimurium"

    with pytest.raises(ValidationError) as caught:
        _web_default_grid(_record_with_organism(SALMONELLA), species=unbreakable)

    error = caught.value
    assert _layout_reason(error) == "DEFINITION_RESERVED"
    message = str(error)
    assert "because the center definition text reserves" in message
    assert "center_reserved_radius" in message
    payload = serialize_web_error(error, operation="generate", stage="render")
    assert payload["code"] == "TRACK_LAYOUT"
    assert payload["context"]["reason"] == "DEFINITION_RESERVED"
    assert set(payload["context"]) <= {"reason", "slotIndex", "innerPx", "outerPx", "keptInsidePx", "keptInsideTracks", "slotTrack"}


@pytest.mark.circular
def test_the_definition_band_message_names_the_tracks_kept_inside(tmp_path: Path) -> None:
    """OV-1268: the band starts above the definition text by the room kept for the tracks inside it."""

    from gbdraw.circular import circular_main

    with pytest.raises(ValidationError) as caught:
        circular_main(
            [
                "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"),
                "--separate_strands", "--track_type", "tuckin", "-l", "left", "--resolve_overlaps",
                "-o", str(tmp_path / "out"), "-f", "svg",
            ]
        )

    assert _layout_reason(caught.value) == "DEFINITION_RESERVED"
    message = str(caught.value)
    assert "cannot fit inside between 189.6px and 386.1px" in message
    assert "the center definition text reserves 122.6px" in message
    assert "and 67.1px are kept for the tracks inside it ('gc_content', 'gc_skew')" in message


def _lambda_resolved_overlaps(tmp_path: Path, *extra: str) -> ValidationError:
    from gbdraw.circular import circular_main

    with pytest.raises(ValidationError) as caught:
        circular_main(
            [
                "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"),
                "--separate_strands", "--track_type", "tuckin", "-l", "left", "--resolve_overlaps",
                "-o", str(tmp_path / "out"), "-f", "svg", *extra,
            ]
        )
    return caught.value


@pytest.mark.circular
def test_the_definition_band_failure_offers_the_gc_tracks_kept_inside(tmp_path: Path) -> None:
    """OV-1274: the GC rows' kept room is a context key, and turning them off or narrowing them is a remedy."""

    error = _lambda_resolved_overlaps(tmp_path)
    assert _layout_reason(error) == "DEFINITION_RESERVED"
    assert (
        "set a smaller center_reserved_radius, turn off or narrow the GC content and GC skew tracks "
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
def test_the_definition_band_failure_starts_its_band_at_the_reserve(tmp_path: Path) -> None:
    """OV-1276: a window with nothing kept inside it starts at the definition band, in the message and the context;
    the scale ticks that do not fit offer hiding the scale."""

    error = _lambda_resolved_overlaps(tmp_path, "--no-gc", "--no-skew")
    assert _layout_reason(error) == "DEFINITION_RESERVED"
    assert (
        "Circular track slot 'ticks' cannot fit inside between 122.6px and 386.1px "
        "because the center definition text reserves 122.6px."
    ) in str(error)
    assert "set a smaller center_reserved_radius, hide the scale (--hide_scale; with an explicit track list, remove the ticks slot), or place tracks outside." in str(error)
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
