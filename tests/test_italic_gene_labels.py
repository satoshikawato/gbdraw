"""Italic gene labels (``labels.<mode>.italic_gene_labels``, ``--italic_gene_labels``, Web "Italic gene names").

A feature label is drawn italic when the setting is on and the qualifier its label is named from (the label
pipeline's qualifier priority) is ``gene``; a label named from ``product``, ``locus_tag`` or another qualifier
stays upright. The setting is off by default and the default figure is unchanged. Italic text is measured
with the italic face the font registry resolves (the bundled Liberation families have one).
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.api import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram
from gbdraw.config.modify import canonical_config_override_paths
from gbdraw.core.text import measurement_font_file

TEST_INPUTS = Path(__file__).parent / "test_inputs"
HMMT = TEST_INPUTS / "HmmtDNA.gbk"
SVG = "{http://www.w3.org/2000/svg}"
BINDING = "data-gbdraw-label-binding-schema"
PATHS = {mode: f"labels.{mode}.italic_gene_labels" for mode in ("circular", "linear")}


def _record():
    return next(SeqIO.parse(str(HMMT), "genbank"))


def _priority(tmp_path: Path) -> str:
    """CDS labels from ``gene`` first; tRNA and rRNA keep the default priority (``product`` first)."""

    path = tmp_path / "priority.tsv"
    path.write_text("CDS\tgene,product\n", encoding="utf-8")
    return str(path)


def _labels(svg_text: str) -> list[tuple[str, str]]:
    """Each drawn feature label as (text, font style): the style of the element that holds the text."""

    labels = []
    for element in ET.fromstring(svg_text).iter():
        if element.get(BINDING) is None:
            continue
        holder = element.find(f"{SVG}textPath")
        holder = element if holder is None else holder
        labels.append(("".join(element.itertext()), str(holder.get("font-style") or "normal")))
    return labels


def _circular(tmp_path: Path, italic: bool | None = None) -> str:
    overrides = {} if italic is None else {PATHS["circular"]: italic}
    options = CircularDiagramOptions(
        config_overrides={"labels.circular.scope": "outer", **overrides},
        qualifier_priority_file=_priority(tmp_path),
    )
    return build_circular_diagram(_record(), options=options).tostring()


def _linear(tmp_path: Path, italic: bool | None = None) -> str:
    overrides = {} if italic is None else {PATHS["linear"]: italic}
    options = LinearDiagramOptions(
        config_overrides={"labels.linear.scope": "all", **overrides},
        qualifier_priority_file=_priority(tmp_path),
    )
    return build_linear_diagram([_record()], options=options).tostring()


def test_the_packaged_config_draws_upright_labels() -> None:
    config = apply_config_overrides(None, None)
    assert config.labels.circular.italic_gene_labels is False
    assert config.labels.linear.italic_gene_labels is False
    on = apply_config_overrides(None, {path: True for path in PATHS.values()})
    assert on.labels.circular.italic_gene_labels is True
    assert on.labels.linear.italic_gene_labels is True
    assert set(PATHS.values()) <= set(canonical_config_override_paths())


@pytest.mark.parametrize("draw", [_circular, _linear], ids=["circular", "linear"])
def test_only_labels_named_from_gene_are_italic(draw, tmp_path: Path) -> None:
    upright = _labels(draw(tmp_path))
    italic = _labels(draw(tmp_path, True))

    assert upright and {style for _text, style in upright} == {"normal"}
    styles = dict(italic)
    # ND1 is a CDS gene name; tRNA-Phe (a tRNA product) and 12S ribosomal RNA (an rRNA product) are not.
    assert styles.get("ND1") == "italic"
    assert {styles[text] for text in styles if text.startswith("tRNA-")} == {"normal"}
    assert any(style == "italic" for _text, style in italic)
    assert any(style == "normal" for _text, style in italic)


def test_the_setting_off_draws_the_default_figure(tmp_path: Path) -> None:
    assert _circular(tmp_path) == _circular(tmp_path, False)
    assert _linear(tmp_path) == _linear(tmp_path, False)


@pytest.mark.parametrize(
    ("module_name", "draw"),
    [("gbdraw.labels.circular", _circular), ("gbdraw.labels.linear", _linear)],
    ids=["circular", "linear"],
)
def test_an_italic_label_is_laid_out_in_its_italic_box(
    module_name: str, draw, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The label box is the italic face's drawn box, so the layout reserves what is painted."""

    import importlib

    placer = importlib.import_module(module_name)
    measured: list[tuple[str, str]] = []
    original = placer.drawn_text_box

    def spy(text, family, size, weight="normal", style="normal"):
        measured.append((str(text), str(style)))
        return original(text, family, size, weight, style)

    monkeypatch.setattr(placer, "drawn_text_box", spy)
    draw(tmp_path, True)

    assert ("ND1", "italic") in measured
    assert ("ND1", "normal") not in measured
    # The bundled family has an italic face, which the measurement reads.
    assert measurement_font_file("Liberation Sans", "normal", "italic") == "LiberationSans-Italic"


def test_a_label_edit_of_a_gene_named_feature_keeps_its_style(tmp_path: Path) -> None:
    """The style follows the qualifier the feature's label is named from, also when a label table row
    replaces the text (the live label-text edit and the next Generate then draw the same)."""

    table = tmp_path / "labels.tsv"
    table.write_text("*\tCDS\tgene\t^ND1$\tNADH dehydrogenase 1\n", encoding="utf-8")
    options = CircularDiagramOptions(
        config_overrides={"labels.circular.scope": "outer", PATHS["circular"]: True},
        qualifier_priority_file=_priority(tmp_path),
        label_override_file=str(table),
    )
    styles = dict(_labels(build_circular_diagram(_record(), options=options).tostring()))

    assert styles.get("NADH dehydrogenase 1") == "italic"
    assert "ND1" not in styles


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_italic_gene_labels(mode: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from gbdraw import cli

    labels = ["--labels", "out"] if mode == "circular" else ["--show_labels", "all"]

    def draw(name: str, *args: str) -> list[tuple[str, str]]:
        output = tmp_path / name
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "gbdraw", mode, "--gbk", str(HMMT), *labels, "--qualifier_priority", _priority(tmp_path),
                *args, "-o", str(output), "-f", "svg",
            ],
        )
        cli.main()
        return _labels(output.with_suffix(".svg").read_text())

    assert {style for _text, style in draw("default")} == {"normal"}
    assert dict(draw("italic", "--italic_gene_labels")).get("ND1") == "italic"
