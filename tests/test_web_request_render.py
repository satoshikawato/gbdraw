from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pytest

from gbdraw.api import (
    CircularBatchRequest,
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.exceptions import ValidationError
from gbdraw.session import build_session_document
import gbdraw.web_support.request_render as web_request_render
from gbdraw.web_support.request_render import (
    render_embedded_canonical_web_request,
)


def test_embedded_web_request_renders_through_typed_api(tmp_path) -> None:
    record = SeqRecord(Seq("ATGC" * 25), id="web-record")
    record.annotations["molecule_type"] = "DNA"
    document = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            output=RenderOutputRequest(
                output_prefix="web-result",
                formats=("svg",),
            ),
        )
    ).to_dict()

    workspace = tmp_path / "render-workspace"
    result = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=workspace,
    )

    assert [item["name"] for item in result["results"]] == ["web-result.svg"]
    assert result["results"][0]["content"].lstrip().startswith("<?xml")
    assert not workspace.exists()

    repeated = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=workspace,
    )

    assert repeated["results"] == result["results"]
    assert not workspace.exists()


def test_embedded_web_request_rejects_duplicate_materialized_names(tmp_path) -> None:
    record = SeqRecord(Seq("ATGC" * 25), id="web-record")
    record.annotations["molecule_type"] = "DNA"
    document = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    ).to_dict()
    resources = dict(document["resources"])
    resource_id, entry = next(iter(resources.items()))
    resources["duplicate-resource"] = dict(entry)

    workspace = tmp_path / "failed-render-workspace"
    try:
        render_embedded_canonical_web_request(
            document["renderRequest"],
            resources=resources,
            workspace=workspace,
        )
    except Exception as exc:
        assert "Duplicate Web render resource filename" in str(exc)
    else:  # pragma: no cover - defensive assertion.
        raise AssertionError("Duplicate materialized names must be rejected.")
    assert not workspace.exists()


def test_embedded_web_request_refuses_and_preserves_preexisting_workspace(
    tmp_path,
) -> None:
    record = SeqRecord(Seq("ATGC" * 25), id="web-record")
    record.annotations["molecule_type"] = "DNA"
    document = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    ).to_dict()
    workspace = tmp_path / "preexisting-workspace"
    workspace.mkdir()
    sentinel = workspace / "owner-data.txt"
    sentinel.write_text("preserve", encoding="utf-8")

    with pytest.raises(
        ValidationError,
        match="Web render workspace must not already exist",
    ):
        render_embedded_canonical_web_request(
            document["renderRequest"],
            resources=document["resources"],
            workspace=workspace,
        )

    assert sentinel.read_text(encoding="utf-8") == "preserve"


def test_embedded_web_request_reports_cleanup_failure_after_success(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = SeqRecord(Seq("ATGC" * 25), id="web-record")
    record.annotations["molecule_type"] = "DNA"
    document = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    ).to_dict()
    workspace = tmp_path / "cleanup-failure-workspace"

    def fail_cleanup(_path) -> None:
        raise OSError("cleanup blocked")

    monkeypatch.setattr(web_request_render.shutil, "rmtree", fail_cleanup)

    with pytest.raises(
        ValidationError,
        match="temporary Web render workspace could not be cleaned up",
    ) as exc_info:
        render_embedded_canonical_web_request(
            document["renderRequest"],
            resources=document["resources"],
            workspace=workspace,
        )

    assert isinstance(exc_info.value.__cause__, OSError)
    assert workspace.exists()


def test_embedded_web_request_preserves_render_error_when_cleanup_fails(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = SeqRecord(Seq("ATGC" * 25), id="web-record")
    record.annotations["molecule_type"] = "DNA"
    document = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    ).to_dict()
    resources = dict(document["resources"])
    _, entry = next(iter(resources.items()))
    resources["duplicate-resource"] = dict(entry)
    workspace = tmp_path / "render-and-cleanup-failure-workspace"

    def fail_cleanup(_path) -> None:
        raise OSError("cleanup blocked")

    monkeypatch.setattr(web_request_render.shutil, "rmtree", fail_cleanup)

    with pytest.raises(
        ValidationError,
        match="Duplicate Web render resource filename",
    ) as exc_info:
        render_embedded_canonical_web_request(
            document["renderRequest"],
            resources=resources,
            workspace=workspace,
        )

    assert any(
        "workspace cleanup also failed: cleanup blocked" in note
        for note in getattr(exc_info.value, "__notes__", ())
    )
    assert workspace.exists()


def test_cleanup_diagnostic_is_attached_without_baseexception_add_note() -> None:
    class LegacyStyleError(RuntimeError):
        add_note = None

    error = LegacyStyleError("primary")

    web_request_render._attach_exception_note(error, "cleanup failed")

    assert error.__notes__ == ["cleanup failed"]


def test_embedded_web_request_preserves_in_memory_comparison_metadata(
    tmp_path,
) -> None:
    records = [
        SeqRecord(Seq("ATGC" * 250), id=record_id)
        for record_id in ("query-record", "subject-record")
    ]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    comparison = pd.DataFrame.from_records(
        [
            {
                "query": "query-record",
                "subject": "subject-record",
                "identity": 91.0,
                "alignment_length": 95,
                "mismatches": 0,
                "gap_opens": 0,
                "qstart": 10,
                "qend": 105,
                "sstart": 20,
                "send": 115,
                "evalue": 1e-20,
                "bitscore": 180,
                "orthogroup_id": "og-web",
                "group_kind": "orthogroup",
                "group_scope": "cross_record",
            }
        ]
    )
    document = build_session_document(
        LinearDiagramRequest(
            records=tuple(
                RecordInput(source=InMemoryRecordSource(record))
                for record in records
            ),
            options=LinearDiagramOptions(protein_comparisons=(comparison,)),
            output=RenderOutputRequest(
                output_prefix="web-comparison",
                formats=("svg",),
            ),
        )
    ).to_dict()

    workspace = tmp_path / "comparison-workspace"
    result = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=workspace,
    )
    svg = result["results"][0]["content"]

    assert not workspace.exists()
    assert 'data-group-kind="orthogroup"' in svg
    assert 'data-group-scope="cross_record"' in svg
    assert 'data-orthogroup-id="og-web"' in svg


@pytest.mark.parametrize(
    ("formats", "output_prefix"),
    [
        (("svg",), "catalog-result"),
        (("interactive_svg",), "catalog-result.interactive"),
    ],
)
def test_web_request_returns_one_base_svg_and_compact_catalog(
    tmp_path,
    formats,
    output_prefix,
) -> None:
    record = SeqRecord(
        Seq("ATGAAATAAGGG"),
        id="catalog-record",
        annotations={"molecule_type": "DNA"},
        features=[
            SeqFeature(
                FeatureLocation(0, 12, strand=1),
                type="source",
                qualifiers={"organism": ["Catalog organism"]},
            ),
            SeqFeature(
                FeatureLocation(0, 9, strand=1),
                type="CDS",
                qualifiers={
                    "locus_tag": ["CATALOG_001"],
                    "product": ["catalog protein"],
                    "translation": ["MK"],
                },
            ),
        ],
    )
    document = build_session_document(
        CircularDiagramRequest(
            records=(
                RecordInput(
                    source=InMemoryRecordSource(record),
                    record_key="catalog-record-key",
                ),
            ),
            output=RenderOutputRequest(
                output_prefix=output_prefix,
                formats=formats,
            ),
        )
    ).to_dict()

    response = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=tmp_path / f"catalog-{'-'.join(formats)}",
    )

    assert [result["name"] for result in response["results"]] == [
        f"{output_prefix}.svg"
    ]
    catalog = response["metadata"]["featureCatalog"]
    assert catalog["schema"] == 3
    assert len(catalog["items"]) == 1
    item = catalog["items"][0]
    assert item["resultIndex"] == 0
    assert item["resultName"] == f"{output_prefix}.svg"
    assert item["recordKeys"] == ["catalog-record-key"]
    assert item["comparisonMatches"] == []
    assert "sequenceSources" not in item

    assert [feature["type"] for feature in item["biologicalFeatures"]] == ["CDS"]
    biological = item["biologicalFeatures"][0]
    assert biological["recordKey"] == "catalog-record-key"
    assert biological["amino_acid_sequence"] == "MK"
    assert "translation" not in biological["qualifiers"]

    assert len(item["features"]) == 1
    rendered = item["features"][0]
    assert rendered["recordKey"] == biological["recordKey"]
    assert (
        rendered["biologicalFeatureId"]
        == biological["biologicalFeatureId"]
    )
    assert rendered["svgId"]
    assert "qualifiers" not in rendered
    assert "nucleotide_sequence" not in rendered
    assert "amino_acid_sequence" not in rendered


def test_web_circular_batch_has_one_catalog_item_per_logical_result(tmp_path) -> None:
    records = []
    for index in range(2):
        record = SeqRecord(
            Seq("ATGAAATAA"),
            id=f"batch-catalog-{index + 1}",
            annotations={"molecule_type": "DNA"},
            features=[
                SeqFeature(
                    FeatureLocation(0, 9, strand=1),
                    type="CDS",
                    qualifiers={"locus_tag": [f"BATCH_{index + 1:03d}"]},
                )
            ],
        )
        records.append(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=f"batch-record-key-{index + 1}",
            )
        )
    document = build_session_document(
        CircularBatchRequest(
            records=tuple(records),
            outputs=tuple(
                RenderOutputRequest(
                    output_prefix=f"batch-catalog-{index + 1}",
                    formats=("interactive_svg",),
                )
                for index in range(2)
            ),
        )
    ).to_dict()

    response = render_embedded_canonical_web_request(
        document["renderRequest"],
        resources=document["resources"],
        workspace=tmp_path / "batch-catalog-workspace",
    )

    assert [result["name"] for result in response["results"]] == [
        "batch-catalog-1.svg",
        "batch-catalog-2.svg",
    ]
    items = response["metadata"]["featureCatalog"]["items"]
    assert [
        (item["resultIndex"], item["resultName"], item["recordKeys"])
        for item in items
    ] == [
        (0, "batch-catalog-1.svg", ["batch-record-key-1"]),
        (1, "batch-catalog-2.svg", ["batch-record-key-2"]),
    ]
