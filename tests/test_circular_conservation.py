from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.analysis.conservation import (
    load_conservation_sources,
    normalize_conservation_tracks_for_record,
)
from gbdraw.api.diagram import assemble_circular_diagram_from_record, build_circular_diagram
from gbdraw.api.options import DiagramOptions
from gbdraw.core.text import calculate_bbox_dimensions
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.render.groups.circular.legend import LegendGroup as CircularLegendGroup
from gbdraw.tracks import CircularTrackSlot, normalize_circular_track_slots


class _BlastConfig:
    evalue = 1e-5
    bitscore = 50.0
    identity = 70.0
    alignment_length = 0


def _record(record_id: str = "rec1", length: int = 120) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    return record


def _comparison_frame(rows: list[tuple[object, ...]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def _hit(
    query: str = "query1",
    subject: str = "rec1",
    qstart: int = 1,
    qend: int = 20,
    sstart: int = 5,
    send: int = 24,
    identity: float = 90.0,
) -> tuple[object, ...]:
    return (
        query,
        subject,
        identity,
        20,
        0,
        0,
        qstart,
        qend,
        sstart,
        send,
        1e-20,
        100.0,
    )


def _translate_xy(transform: str | None) -> tuple[float, float]:
    assert transform is not None
    assert transform.startswith("translate(")
    parts = transform.removeprefix("translate(").removesuffix(")").replace(",", " ").split()
    assert len(parts) >= 2
    return float(parts[0]), float(parts[1])


def test_conservation_loader_keeps_logical_source_indexes_after_skip(tmp_path: Path) -> None:
    missing = tmp_path / "missing.tsv"
    valid = _comparison_frame([_hit()])

    result = load_conservation_sources(
        blast_config=_BlastConfig(),
        conservation_files=[str(missing)],
        conservation_dataframes=[pd.DataFrame(), valid],
        labels=["missing source", "valid source"],
    )

    assert result.sources[0].source_index == 0
    assert result.sources[0].skipped is True
    assert result.sources[1].source_index == 1
    assert result.sources[1].label == "valid source"

    tracks = normalize_conservation_tracks_for_record(
        result,
        displayed_records=[_record()],
        record=_record(),
        conservation_reference="subject",
    )
    assert len(tracks) == 1
    assert tracks[0].source_index == 1
    assert tracks[0].track_index == 1
    assert tracks[0].track_label == "valid source"


def test_conservation_normalization_uses_subject_reverse_coordinates() -> None:
    result = load_conservation_sources(
        blast_config=_BlastConfig(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=10, send=5)])
        ],
        labels=["subject hits"],
    )

    tracks = normalize_conservation_tracks_for_record(
        result,
        displayed_records=[_record()],
        record=_record(),
        conservation_reference="auto",
    )
    row = tracks[0].hits.iloc[0]

    assert tracks[0].reference_side == "subject"
    assert row["start"] == 5
    assert row["end"] == 10
    assert row["draw_start"] == pytest.approx(4.0)
    assert row["draw_end"] == pytest.approx(10.0)
    assert row["orientation"] == "reverse"


def test_conservation_empty_valid_source_keeps_empty_ring() -> None:
    result = load_conservation_sources(
        blast_config=_BlastConfig(),
        conservation_dataframes=[
            _comparison_frame([_hit(identity=50.0)])
        ],
        labels=["filtered out"],
    )

    assert result.sources[0].skipped is False
    tracks = normalize_conservation_tracks_for_record(
        result,
        displayed_records=[_record()],
        record=_record(),
        conservation_reference="subject",
    )

    assert len(tracks) == 1
    assert tracks[0].reference_side == "subject"
    assert tracks[0].hits.empty


def test_conservation_empty_file_keeps_empty_ring(tmp_path: Path) -> None:
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("", encoding="utf-8")

    result = load_conservation_sources(
        blast_config=_BlastConfig(),
        conservation_files=[str(empty_file)],
        labels=["empty source"],
    )

    assert result.sources[0].skipped is False
    tracks = normalize_conservation_tracks_for_record(
        result,
        displayed_records=[_record()],
        record=_record(),
        conservation_reference="subject",
    )
    assert len(tracks) == 1
    assert tracks[0].hits.empty


def test_conservation_auto_reference_ambiguity_requires_explicit_side() -> None:
    result = load_conservation_sources(
        blast_config=_BlastConfig(),
        conservation_dataframes=[
            _comparison_frame([_hit(query="rec1", subject="rec1")])
        ],
    )

    with pytest.raises(ValidationError, match="both BLAST sides"):
        normalize_conservation_tracks_for_record(
            result,
            displayed_records=[_record()],
            record=_record(),
            conservation_reference="auto",
        )


def test_sequence_conservation_slot_rejects_overlay_side() -> None:
    with pytest.raises(ValueError, match="cannot use side=overlay"):
        normalize_circular_track_slots(
            [
                CircularTrackSlot(
                    id="conservation_1",
                    renderer="sequence_conservation",
                    side="overlay",
                )
            ]
        )


def test_circular_api_renders_conservation_ring_and_gradient_legend() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=1, send=120)])
        ],
        conservation_reference="subject",
        conservation_labels=["Reference A"],
        conservation_ring_width=12,
        conservation_ring_gap=4,
        legend="right",
    )
    svg = canvas.tostring()

    assert 'id="conservation_Reference_A"' in svg
    assert 'data-track-label="Reference A"' in svg
    assert 'data-reference-record-id="rec1"' in svg
    assert 'data-legend-key="Conservation identity"' in svg


def test_circular_api_renders_source_colored_conservation_ring() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=1, send=120)])
        ],
        conservation_reference="subject",
        conservation_labels=["barcode13"],
        conservation_colors=["#E15759"],
        legend="right",
    )
    svg = canvas.tostring()

    assert 'id="conservation_barcode13"' in svg
    assert 'data-track-label="barcode13"' in svg
    assert 'data-track-color="#e15759"' in svg
    assert 'data-legend-key="barcode13"' in svg
    assert "#e15759" in svg


def test_circular_multi_conservation_gradient_legend_uses_compact_linear_layout() -> None:
    canvas_config = SimpleNamespace(legend_position="right", dpi=96)
    legend_config = SimpleNamespace(
        font_family="'Liberation Sans', 'Arial', sans-serif",
        font_size=10.0,
        color_rect_size=12.0,
        legend_width=180.0,
        legend_height=80.0,
        pairwise_legend_width=160.0,
    )
    legend_table = {
        "barcode07.draft": {
            "type": "gradient",
            "min_color": "#dbe8f5",
            "max_color": "#4e79a7",
            "stroke": "none",
            "width": 0,
            "min_value": 0,
        },
        "barcode15.draft": {
            "type": "gradient",
            "min_color": "#fdebd8",
            "max_color": "#f28e2b",
            "stroke": "none",
            "width": 0,
            "min_value": 0,
        },
    }

    drawing = Drawing(debug=False)
    drawing.add(CircularLegendGroup(canvas_config, legend_config, legend_table).get_group())
    root = ET.fromstring(drawing.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    legend = root.find(".//svg:g[@id='conservation_identity_legend']", ns)
    assert legend is not None
    entries = legend.findall("svg:g[@data-legend-key]", ns)
    assert [entry.get("data-legend-key") for entry in entries] == [
        "barcode07.draft",
        "barcode15.draft",
    ]

    first_entry = entries[0]
    first_label = first_entry.find("svg:text", ns)
    assert first_label is not None
    label_x, _ = _translate_xy(first_label.get("transform"))
    assert label_x == pytest.approx(0.0)

    gradient_bar = next(
        path
        for path in first_entry.findall("svg:path", ns)
        if str(path.get("fill", "")).startswith("url(")
    )
    bar_x, _ = _translate_xy(gradient_bar.get("transform"))
    label_width, _ = calculate_bbox_dimensions(
        "barcode07.draft",
        legend_config.font_family,
        legend_config.font_size,
        canvas_config.dpi,
    )
    assert bar_x == pytest.approx(label_width + 0.2 * legend_config.color_rect_size)

    legend_texts = [text.text for text in legend.findall(".//svg:text", ns)]
    assert legend_texts.count("0%") == 1
    assert legend_texts.count("100%") == 1


def test_circular_api_uses_explicit_conservation_slot_source_indexes() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=1, send=60)]),
            _comparison_frame([_hit(subject="rec1", sstart=61, send=120)]),
        ],
        conservation_reference="subject",
        conservation_labels=["Reference A", "Reference B"],
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="ticks", renderer="ticks"),
            CircularTrackSlot(
                id="conservation_b",
                renderer="sequence_conservation",
                params={"source_index": 1},
            ),
            CircularTrackSlot(
                id="conservation_a",
                renderer="sequence_conservation",
                params={"source_index": 0},
            ),
        ],
        legend="right",
    )
    svg = canvas.tostring()

    assert 'id="conservation_Reference_B"' in svg
    assert 'id="conservation_Reference_A"' in svg
    assert svg.index('id="conservation_Reference_B"') < svg.index('id="conservation_Reference_A"')


def test_circular_api_keeps_axis_derived_side_for_conservation_slot() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=1, send=120)])
        ],
        conservation_reference="subject",
        conservation_labels=["Outer conservation"],
        circular_track_slots=[
            CircularTrackSlot(
                id="conservation_outer",
                renderer="sequence_conservation",
                params={"source_index": 0},
            ),
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="ticks", renderer="ticks"),
        ],
        circular_track_axis_index=1,
        legend="right",
    )
    svg = canvas.tostring()

    assert 'id="conservation_Outer_conservation"' in svg
    assert 'data-track-label="Outer conservation"' in svg


def test_disabled_explicit_conservation_slot_suppresses_auto_insert() -> None:
    canvas = assemble_circular_diagram_from_record(
        _record(),
        conservation_dataframes=[
            _comparison_frame([_hit(subject="rec1", sstart=1, send=120)])
        ],
        conservation_reference="subject",
        conservation_labels=["Hidden reference"],
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(
                id="hidden_conservation",
                renderer="sequence_conservation",
                enabled=False,
            ),
        ],
        legend="right",
    )
    svg = canvas.tostring()

    assert 'id="conservation_Hidden_reference"' not in svg


def test_circular_diagram_options_forward_conservation_dataframe() -> None:
    canvas = build_circular_diagram(
        _record(),
        options=DiagramOptions(
            conservation_dataframes=[
                _comparison_frame([_hit(subject="rec1", sstart=1, send=120)])
            ],
            conservation_reference="subject",
            conservation_labels=["Option ring"],
            conservation_colors=["red"],
            conservation_ring_width=10,
            conservation_ring_gap=3,
        ),
    )
    svg = canvas.tostring()

    assert 'id="conservation_Option_ring"' in svg
    assert 'data-track-label="Option ring"' in svg
    assert 'data-track-color="#ff0000"' in svg


def test_circular_api_rejects_nonpositive_conservation_geometry() -> None:
    with pytest.raises(ValidationError, match="conservation_ring_gap must be > 0"):
        assemble_circular_diagram_from_record(
            _record(),
            conservation_dataframes=[
                _comparison_frame([_hit(subject="rec1")])
            ],
            conservation_reference="subject",
            conservation_ring_gap=0,
        )
