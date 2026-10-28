from __future__ import annotations

from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.requests import RecordPresentation
from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.layout.spatial import Aabb
from gbdraw.diagrams.circular.positioning import place_definition_group_on_canvas
from gbdraw.core.record_metadata import infer_record_source_metadata
from gbdraw.render.groups.circular.definition import CircularDefinitionStyle, DefinitionGroup, DefinitionNameLines
from svgwrite.container import Group


def _definition(
    *,
    plot_title: str | None = None,
    profile: str = "full",
    record_label: str = "",
    record_subtitle: str = "",
) -> DefinitionGroup:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    render_profile = CircularRenderProfile(cfg)
    record = SeqRecord(Seq("ATGC" * 250), id="NC_TEST.1")
    canvas = CircularCanvasConfigurator("test", render_profile, "none", record)
    style = CircularDefinitionStyle.of(cfg, profile)
    return DefinitionGroup(
        record,
        canvas.geometry(),
        style=style,
        names=DefinitionNameLines.of(
            infer_record_source_metadata(record),
            style.profile,
            species="<i>Testus boundsii</i>",
            strain="strain A",
            plot_title=plot_title,
            record_presentation=RecordPresentation(label=record_label, subtitle=record_subtitle),
        ),
        definition_group_id=("plot_title" if profile == "shared_common" else None),
    )


def test_circular_definition_exposes_authoritative_local_bounds() -> None:
    definition = _definition()

    assert isinstance(definition.local_bounds, Aabb)
    assert definition.local_bounds.width > 0.0
    assert definition.local_bounds.height > 0.0
    for element in definition.get_group().elements:
        x = float(element.attribs["x"])
        y = float(element.attribs["y"])
        assert definition.local_bounds.min_x < x < definition.local_bounds.max_x
        assert definition.local_bounds.min_y < y < definition.local_bounds.max_y


def test_circular_plot_title_bounds_are_centered_in_local_coordinates() -> None:
    title = _definition(
        plot_title="Comparative <i>mitochondrial</i> rings",
        profile="shared_common",
    )

    assert title.get_group().attribs["data-gbdraw-role"] == "plot-title"
    assert title.local_bounds.min_x == pytest.approx(-title.local_bounds.max_x)
    assert title.local_bounds.width > 100.0
    assert title.local_bounds.height > 0.0


def test_circular_definition_local_bounds_center_on_record_axis() -> None:
    group = Group(debug=False)
    local_bounds = Aabb(-40.0, -10.0, 60.0, 50.0)
    setattr(group, "_gbdraw_local_bounds", local_bounds)

    place_definition_group_on_canvas(
        group,
        SimpleNamespace(offset_x=400.0, offset_y=300.0),
    )

    translate_x, translate_y = getattr(group, "_gbdraw_plot_translation")
    assert translate_x + 0.5 * (local_bounds.min_x + local_bounds.max_x) == pytest.approx(400.0)
    assert translate_y + 0.5 * (local_bounds.min_y + local_bounds.max_y) == pytest.approx(300.0)


def test_circular_record_label_and_subtitle_override_inferred_title_lines() -> None:
    definition = _definition(
        record_label="長い <i>環状ゲノム</i> ラベル",
        record_subtitle="Selected record subtitle",
    )

    svg = definition.get_group().tostring()

    assert "長い " in svg
    assert "環状ゲノム" in svg
    assert 'font-style="italic"' in svg
    assert "Selected record subtitle" in svg
    assert "Testus boundsii" not in svg
    assert "strain A" not in svg
    assert definition.local_bounds.width > 100.0


def test_empty_circular_record_title_values_preserve_inference() -> None:
    definition = _definition()
    svg = definition.get_group().tostring()

    assert "Testus boundsii" in svg
    assert "strain A" in svg


def _k12_definition(
    *,
    record_label: str = "",
    species: str | None = None,
    organism: str = "Escherichia coli",
    strain: str = "K-12",
    taxonomy: tuple[str, ...] = ("Bacteria", "Pseudomonadati"),
) -> DefinitionGroup:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    record = SeqRecord(Seq("ATGC" * 250), id="NC_K12.1")
    qualifiers = {"organism": [organism]}
    if strain:
        qualifiers["strain"] = [strain]
    record.features = [SeqFeature(SimpleLocation(0, len(record.seq)), type="source", qualifiers=qualifiers)]
    record.annotations["taxonomy"] = list(taxonomy)
    canvas = CircularCanvasConfigurator("test", CircularRenderProfile(cfg), "none", record)
    return DefinitionGroup(
        record,
        canvas.geometry(),
        style=CircularDefinitionStyle.of(cfg, "full"),
        names=DefinitionNameLines.of(
            infer_record_source_metadata(record),
            "full",
            species=species,
            record_presentation=RecordPresentation(label=record_label),
        ),
    )


def _plain(parts: list[dict[str, str | bool | None]]) -> str:
    return "".join(str(part["text"] or "") for part in parts)


def test_a_typed_record_label_is_drawn_whole_over_the_inferred_strain_line() -> None:
    # PD-OI-012: an explicit Record label overrides its line as typed, and the
    # empty Subtitle keeps the inferred strain line. Only a typed drawing-level
    # Species is cut by the Auto rule (DF-D09).
    labelled = _k12_definition(record_label="Escherichia coli MG1655")
    assert (_plain(labelled.species_parts), _plain(labelled.strain_parts)) == ("Escherichia coli MG1655", "K-12")
    typed_species = _k12_definition(species="Escherichia coli MG1655")
    assert (_plain(typed_species.species_parts), _plain(typed_species.strain_parts)) == ("Escherichia coli", "MG1655")


def test_a_typed_line_one_draws_the_strain_the_organism_named() -> None:
    # The /organism "Lachnospiraceae bacterium A4" names its strain, so the Auto
    # line 2 is empty. A typed line 1 replaces that name: line 2 draws "A4".
    for typed in ({"record_label": "My isolate"}, {"species": "My isolate"}):
        definition = _k12_definition(organism="Lachnospiraceae bacterium A4", strain="A4", taxonomy=("Bacteria",), **typed)
        assert (_plain(definition.species_parts), _plain(definition.strain_parts)) == ("My isolate", "A4")


def test_a_candidatus_name_with_markup_characters_draws_the_italic_ca_and_the_literal_name() -> None:
    from gbdraw.core.record_metadata import format_inferred_definition
    from gbdraw.render.groups.linear import DefinitionGroup as LinearDefinitionGroup, LinearDefinitionStyle

    circular = _k12_definition(organism="Candidatus Foo & Bar", strain="", taxonomy=("Bacteria",))
    expected = [("Ca.", True), (" Foo & Bar", False)]
    assert [(part["text"], part["italic"]) for part in circular.species_parts] == expected
    svg = circular.get_group().tostring()
    assert "Foo &amp; Bar" in svg and "&lt;i&gt;" not in svg
    # The Web Linear Auto definition is format_inferred_definition's one line.
    record = SeqRecord(Seq("ATGC" * 250), id="NC_CA.1", annotations={"taxonomy": ["Bacteria"]})
    record.features = [
        SeqFeature(SimpleLocation(0, len(record.seq)), type="source", qualifiers={"organism": ["Candidatus Foo & Bar"]})
    ]
    definition = format_inferred_definition(infer_record_source_metadata(record))
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    linear = LinearDefinitionGroup(
        record, "short", style=LinearDefinitionStyle.of(cfg), record_presentation=RecordPresentation(label=definition)
    )
    linear_svg = linear.get_group().tostring()
    assert "Foo &amp; Bar" in linear_svg and "&lt;i&gt;" not in linear_svg
    assert "font-style=\"italic\"" in linear_svg or "font-style:italic" in linear_svg
