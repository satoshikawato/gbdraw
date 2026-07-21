from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    DiagramOptions,
    InMemoryRecordSource,
    LinearComparison,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordInput,
    assemble_linear_diagram_from_records,
    read_comparisons_table,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import merge_linear_comparisons, validate_linear_comparison_topology
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request


def _comparison(query: int, subject: int) -> LinearComparison:
    row = ["q", "s", 90.0, 100, 0, 0, 10, 100, 20, 110, 1e-20, 200]
    return LinearComparison(query, subject, pd.DataFrame([row], columns=COMPARISON_COLUMNS))


def _records() -> list[SeqRecord]:
    records = [SeqRecord(Seq("A" * 1000), id=f"r{index}") for index in range(1, 5)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


def test_topology_validation_accepts_selected_adjacent_row_pairs() -> None:
    validate_linear_comparison_topology(
        (_comparison(0, 2), _comparison(1, 2), _comparison(1, 3)),
        (0, 0, 1, 1),
    )


@pytest.mark.parametrize("edge", [(0, 1), (0, 3)])
def test_topology_validation_rejects_same_or_non_adjacent_rows(edge) -> None:
    rows = (0, 0, 1, 2)
    with pytest.raises(ValidationError):
        validate_linear_comparison_topology((_comparison(*edge),), rows)


def test_multiple_sources_for_one_pair_are_merged() -> None:
    merged = merge_linear_comparisons((_comparison(0, 2), _comparison(0, 2)))
    assert len(merged) == 1
    assert len(merged[0].matches) == 2


def test_explicit_comparison_metadata_uses_selected_endpoints() -> None:
    canvas = assemble_linear_diagram_from_records(
        _records(),
        linear_comparisons=[_comparison(1, 2)],
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
        legend="none",
        config_overrides={"show_labels": False, "show_gc": False, "show_skew": False},
    )
    svg = canvas.tostring()
    assert 'data-query-record-index="1"' in svg
    assert 'data-subject-record-index="2"' in svg
    assert 'data-query-row="0"' in svg
    assert 'data-subject-row="1"' in svg


def test_comparison_ribbons_attach_to_painted_occupancy_without_extra_padding() -> None:
    records = _records()[:3]
    for record in records:
        record.features.extend(
            (
                SeqFeature(FeatureLocation(10, 110, strand=1), type="CDS"),
                SeqFeature(FeatureLocation(120, 220, strand=-1), type="CDS"),
            )
        )
    canvas = assemble_linear_diagram_from_records(
        records,
        linear_comparisons=[_comparison(0, 2)],
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1", "#3@2"),
        ),
        legend="none",
        config_overrides={
            "show_labels": False,
            "show_gc": False,
            "show_skew": False,
            "strandedness": True,
            "linear_track_layout": "above",
            "linear_ruler_on_axis": False,
        },
    )
    svg = canvas.tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    groups = {
        group.attrib["id"]: group
        for group in root.findall(".//svg:g", namespace)
        if "id" in group.attrib
    }

    def translate_y(group: ET.Element) -> float:
        match = re.fullmatch(
            r"translate\(([-+0-9.eE]+),([-+0-9.eE]+)\)",
            group.attrib["transform"],
        )
        assert match is not None
        return float(match.group(2))

    def path_y_values(group: ET.Element) -> list[float]:
        return [
            float(y_value)
            for path in group.findall(".//svg:path", namespace)
            for _x_value, y_value in re.findall(
                r"[ML]\s*([-+0-9.eE]+),?\s*([-+0-9.eE]+)",
                path.attrib.get("d", ""),
            )
        ]

    top_record = groups["r1_record_1"]
    bottom_record = groups["r3_record_3"]
    comparison = groups["comparison1"]
    top_record_axis = translate_y(top_record)
    top_record_bottom = top_record_axis + max(path_y_values(top_record))
    bottom_record_top = translate_y(bottom_record) + min(path_y_values(bottom_record))
    comparison_y_values = path_y_values(comparison)
    comparison_top = translate_y(comparison) + min(comparison_y_values)
    comparison_bottom = translate_y(comparison) + max(comparison_y_values)

    geometry = canvas._gbdraw_track_slot_geometry["records"]
    assert top_record_bottom < top_record_axis
    assert bottom_record_top < translate_y(bottom_record)
    assert comparison_top == pytest.approx(
        geometry[0]["comparisonExclusionBand"]["absoluteBottomPx"]
    )
    assert comparison_bottom == pytest.approx(
        geometry[2]["comparisonExclusionBand"]["absoluteTopPx"]
    )
    bottom_features = next(
        slot for slot in geometry[2]["slots"] if slot["renderer"] == "features"
    )
    assert comparison_bottom == pytest.approx(
        bottom_features["paintBand"]["absoluteTopPx"]
    )


def test_reversed_explicit_endpoints_work_with_legacy_one_record_rows() -> None:
    canvas = assemble_linear_diagram_from_records(
        _records()[:2],
        linear_comparisons=[_comparison(1, 0)],
        legend="none",
        config_overrides={"show_labels": False, "show_gc": False, "show_skew": False},
    )
    svg = canvas.tostring()
    assert 'data-query-record-index="1"' in svg
    assert 'data-subject-record-index="0"' in svg
    assert 'data-query-row="1"' in svg
    assert 'data-subject-row="0"' in svg


def test_selected_generated_protein_pairs_keep_explicit_endpoints(monkeypatch) -> None:
    calls: list[tuple[str, str]] = []

    def fake_pairwise(records, **_kwargs):
        calls.append((records[0].id, records[1].id))
        return SimpleNamespace(comparisons=[_comparison(0, 1).matches])

    monkeypatch.setattr(
        "gbdraw.api.diagram.build_pairwise_protein_blastp_comparisons",
        fake_pairwise,
    )
    canvas = assemble_linear_diagram_from_records(
        _records(),
        protein_blastp_mode="pairwise",
        protein_comparison_pairs=((0, 2), (1, 3)),
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
        legend="none",
        config_overrides={"show_labels": False, "show_gc": False, "show_skew": False},
    )
    assert calls == [("r1", "r3"), ("r2", "r4")]
    svg = canvas.tostring()
    assert 'data-query-record-index="0"' in svg
    assert 'data-subject-record-index="2"' in svg


def test_collinear_all_scope_renders_every_cross_row_pair(monkeypatch) -> None:
    captured_pairs: tuple[tuple[int, int], ...] | None = None

    def fake_collinear(_records, **kwargs):
        nonlocal captured_pairs
        captured_pairs = kwargs["comparison_pairs"]
        return SimpleNamespace(orthogroups=None)

    def fake_convert(_result, **_kwargs):
        assert captured_pairs is not None
        return {
            pair: _comparison(*pair).matches
            for pair in captured_pairs
        }

    monkeypatch.setattr(
        "gbdraw.api.diagram.build_orthogroup_collinearity_blocks",
        fake_collinear,
    )
    monkeypatch.setattr(
        "gbdraw.api.diagram.convert_collinearity_blocks_to_pair_comparisons",
        fake_convert,
    )
    canvas = assemble_linear_diagram_from_records(
        _records(),
        protein_blastp_mode="collinear",
        collinearity_search_scope="all",
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
        legend="none",
        config_overrides={"show_labels": False, "show_gc": False, "show_skew": False},
    )
    assert captured_pairs == ((0, 2), (0, 3), (1, 2), (1, 3))
    svg = canvas.tostring()
    rendered_pairs = {
        (int(query), int(subject))
        for query, subject in re.findall(
            r'data-query-record-index="(\d+)"[^>]+data-subject-record-index="(\d+)"',
            svg,
        )
    }
    assert rendered_pairs == set(captured_pairs)


def test_comparisons_table_resolves_relative_blast_path(tmp_path) -> None:
    blast = tmp_path / "pair.tsv"
    blast.write_text("q\ts\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n", encoding="utf-8")
    table_path = tmp_path / "comparisons.tsv"
    table_path.write_text("blast\tquery\tsubject\npair.tsv\t#1\t#3\n", encoding="utf-8")
    table = read_comparisons_table(str(table_path))
    assert table.rows[0].blast == str(blast)
    assert table.path_dependencies[0].column == "blast"


def test_current_schema_preserves_record_keys_layout_and_explicit_endpoints(tmp_path) -> None:
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=f"stable-{index}",
            )
            for index, record in enumerate(_records(), start=1)
        ),
        options=DiagramOptions(linear_comparisons=(_comparison(1, 2),)),
        layout=LinearMultiRecordOptions(
            record_gap_px=30,
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
    )
    encoded = encode_canonical_request(request)
    assert encoded.payload["schema"] == 3
    assert encoded.payload["records"][0]["recordKey"] == "stable-1"
    assert encoded.payload["comparisons"][0]["queryRecordIndex"] == 1

    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        if resource.content is not None:
            target.write_bytes(resource.content)
        else:
            target = resource.source_path
        resource_paths[resource.resource_id] = target
    decoded = decode_canonical_request(
        encoded.payload,
        resource_paths=resource_paths,
        output_directory=tmp_path / "out",
    )
    assert isinstance(decoded, LinearDiagramRequest)
    assert decoded.records[0].record_key == "stable-1"
    assert decoded.layout == request.layout
    assert decoded.options.linear_comparisons[0].query_record_index == 1
    assert decoded.options.linear_comparisons[0].subject_record_index == 2
