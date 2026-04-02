from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram_module
import gbdraw.linear as linear_cli_module
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import DiagramOptions
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import load_comparisons


def _build_record() -> SeqRecord:
    return SeqRecord(Seq("A" * 400), id="record_1")


def _write_blast_table(path: Path, rows: list[tuple[object, ...]]) -> None:
    contents = [
        "# Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
    ]
    contents.extend("\t".join(map(str, row)) for row in rows)
    path.write_text("\n".join(contents) + "\n", encoding="utf-8")


@pytest.mark.linear
def test_load_comparisons_filters_hits_below_alignment_length(tmp_path: Path) -> None:
    blast_path = tmp_path / "comparisons.tsv"
    _write_blast_table(
        blast_path,
        [
            ("q1", "s1", 95, 80, 0, 0, 1, 80, 1, 80, 1e-40, 500),
            ("q1", "s1", 95, 150, 0, 0, 1, 150, 1, 150, 1e-40, 500),
        ],
    )

    config = SimpleNamespace(evalue=1e-2, bitscore=50, identity=0, alignment_length=100)

    comparisons = load_comparisons([str(blast_path)], config)

    assert len(comparisons) == 1
    assert comparisons[0]["alignment_length"].tolist() == [150]


@pytest.mark.linear
def test_load_comparisons_applies_alignment_length_with_other_thresholds(tmp_path: Path) -> None:
    blast_path = tmp_path / "comparisons.tsv"
    _write_blast_table(
        blast_path,
        [
            ("q1", "s1", 95, 90, 0, 0, 1, 90, 1, 90, 1e-40, 500),
            ("q1", "s1", 60, 180, 0, 0, 1, 180, 1, 180, 1e-40, 500),
            ("q1", "s1", 95, 180, 0, 0, 1, 180, 1, 180, 1e-40, 500),
        ],
    )

    config = SimpleNamespace(evalue=1e-10, bitscore=50, identity=90, alignment_length=100)

    comparisons = load_comparisons([str(blast_path)], config)

    assert len(comparisons) == 1
    assert comparisons[0].shape[0] == 1
    assert comparisons[0]["alignment_length"].tolist() == [180]
    assert comparisons[0]["identity"].tolist() == [95]


@pytest.mark.linear
def test_linear_cli_alignment_length_is_forwarded(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _build_record()
    captured: dict[str, int] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured["alignment_length"] = kwargs.get("alignment_length")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--alignment_length",
            "123",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["alignment_length"] == 123


@pytest.mark.linear
def test_build_linear_diagram_forwards_alignment_length(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, int] = {}

    def fake_assemble(*_args, **kwargs):
        captured["alignment_length"] = kwargs.get("alignment_length")
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_build_record()],
        options=DiagramOptions(alignment_length=321),
    )

    assert isinstance(canvas, Drawing)
    assert captured["alignment_length"] == 321


@pytest.mark.linear
def test_assemble_linear_diagram_rejects_negative_alignment_length() -> None:
    with pytest.raises(ValidationError, match="alignment_length must be >= 0"):
        assemble_linear_diagram_from_records(
            [_build_record()],
            alignment_length=-1,
        )
