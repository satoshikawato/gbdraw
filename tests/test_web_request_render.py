from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pytest

from gbdraw.api import (
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
