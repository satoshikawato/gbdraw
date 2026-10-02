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
):
    overrides = dict(WEB_DEFAULT_OVERRIDES)
    overrides.update(config_overrides or {})
    return build_circular_multi_diagram(
        [record],
        options=CircularDiagramOptions(
            config_overrides=overrides,
            output=CircularOutputOptions(legend="left"),
            tracks=CircularRequestTrackOptions(center_reserved_radius=center_reserved_radius),
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
@pytest.mark.parametrize("organism", [SALMONELLA, MYCOBACTERIUM, HEPATOPLASMA])
def test_web_default_cli_generates_long_organism_names(
    organism: str, tmp_path: Path, gbdraw_runner
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
    assert " ".join(lines).startswith(organism)


@pytest.mark.circular
def test_species_line_wraps_at_word_boundaries_only_when_placement_fails() -> None:
    wrapped = _definition_lines(_web_default_grid(_record_with_organism(SALMONELLA)).tostring())
    species_lines = wrapped[: len(wrapped) - 4]  # mitochondrion, accession, length, GC follow

    assert len(species_lines) >= 2
    assert " ".join(species_lines) == SALMONELLA
    assert all(line == line.strip() and line for line in species_lines)

    for organism in (MYCOBACTERIUM, ESCHERICHIA):
        lines = _definition_lines(_web_default_grid(_record_with_organism(organism)).tostring())
        assert lines[0] == organism


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
            _record_with_organism(SALMONELLA),
            config_overrides={"objects.definition.circular.font_size": 17},
        )

    assert _layout_reason(caught.value) == "DEFINITION_RESERVED"
    assert "center definition text reserves" in str(caught.value)


@pytest.mark.circular
def test_explicit_center_reserved_radius_disables_the_wrap() -> None:
    lines = _definition_lines(
        _web_default_grid(
            _record_with_organism(SALMONELLA), center_reserved_radius=150.0
        ).tostring()
    )

    assert lines[0] == SALMONELLA

    # A radius as large as the one-line definition reservation (241.9 px) fails
    # with the Web defaults; an explicit radius must not be replaced by a
    # wrapped definition reservation.
    with pytest.raises(ValidationError) as caught:
        _web_default_grid(_record_with_organism(SALMONELLA), center_reserved_radius=250.0)

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
    assert set(payload["context"]) <= {"reason", "slotIndex", "innerPx", "outerPx"}
