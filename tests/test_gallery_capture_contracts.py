from __future__ import annotations

import base64
import json
from pathlib import Path

import pytest

from tools.capture_gallery_tutorial_screenshots import (
    DEFAULT_MAX_FILE_SIZE_KB,
    DEFAULT_MAX_IMAGE_HEIGHT,
    DEFAULT_MAX_IMAGE_WIDTH,
    TutorialContext,
    ValidationResult,
    add_capture_contract_validation,
    iter_operation_contexts,
    load_ready_examples,
    resolve_gallery_reference,
    validate_tutorial_media,
    wait_for_web_app_ready,
)


BASIC_EXAMPLES = {
    "HmmtDNA_basic_circular": {
        "operation_count": 6,
        "data_dependent_count": 5,
        "filename": "HmmtDNA.gbk",
        "resource_marker": b"VERSION     NC_012920.1",
    },
    "lambda_basic_linear": {
        "operation_count": 9,
        "data_dependent_count": 8,
        "filename": "NC_001416.gb",
        "resource_marker": b"VERSION     NC_001416.1",
    },
}


def _load_tutorial(sample: dict[str, object]) -> dict[str, object]:
    tutorial_path = resolve_gallery_reference(str(sample["tutorial"]))
    assert tutorial_path is not None
    return json.loads(tutorial_path.read_text(encoding="utf-8"))


@pytest.mark.parametrize("example_id", BASIC_EXAMPLES)
def test_basic_gallery_operations_have_owned_semantic_capture_contracts(
    example_id: str,
) -> None:
    expected = BASIC_EXAMPLES[example_id]
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    operations = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
    ]
    data_dependent = [
        operation for operation in operations if operation.get("dataDependent") is True
    ]

    assert len(operations) == expected["operation_count"]
    assert len(data_dependent) == expected["data_dependent_count"]

    for operation in data_dependent:
        media = operation["media"]
        capture = operation["capture"]
        assert media["src"].startswith(f"./media/{example_id}/")
        assert capture["source"] == "webapp"
        assert capture["session"] == sample["session"]
        assert capture["assertAppState"]
        assert capture.get("visibleControls") or capture.get("visibleText")

    visible_identity = {
        text
        for operation in data_dependent
        for text in operation["capture"].get("visibleText", [])
    }
    assert expected["filename"] in visible_identity


@pytest.mark.parametrize("example_id", BASIC_EXAMPLES)
def test_basic_gallery_capture_session_contains_the_named_record(
    example_id: str,
) -> None:
    expected = BASIC_EXAMPLES[example_id]
    sample = load_ready_examples(example_id)[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    session = json.loads(session_path.read_text(encoding="utf-8"))
    resource = session["resources"]["record-1-genbank"]

    assert resource["encoding"] == "base64"
    decoded = base64.b64decode(resource["data"])
    assert expected["resource_marker"] in decoded


def test_ready_gallery_media_passes_capture_contract_validation() -> None:
    result = validate_tutorial_media(
        load_ready_examples(),
        max_width=DEFAULT_MAX_IMAGE_WIDTH,
        max_height=DEFAULT_MAX_IMAGE_HEIGHT,
        max_file_size_kb=DEFAULT_MAX_FILE_SIZE_KB,
    )

    assert result.errors == []


def test_data_dependent_capture_rejects_cross_example_media() -> None:
    context = TutorialContext(
        example_id="current-example",
        tutorial_path=Path("tutorial.json"),
        section="manualSteps",
        step_index=1,
        step_title="Upload",
        operation_index=1,
    )
    sample = {"session": "./sessions/current-example.gbdraw-session.json"}
    operation = {
        "dataDependent": True,
        "media": {"src": "./media/other-example/upload.webp"},
        "capture": {
            "source": "webapp",
            "session": "./sessions/current-example.gbdraw-session.json",
            "assertAppState": {"files.c_gb.name": "current.gbk"},
            "visibleText": ["current.gbk"],
        },
    }
    result = ValidationResult(errors=[], warnings=[])

    add_capture_contract_validation(result, context, sample, operation)

    assert any(
        "data-dependent media must be owned by current-example" in error
        for error in result.errors
    )


def test_cross_example_media_is_explicitly_generic() -> None:
    for sample in load_ready_examples():
        tutorial = _load_tutorial(sample)
        for _, operation in iter_operation_contexts(sample, tutorial):
            media = operation.get("media")
            if not isinstance(media, dict):
                continue
            source = str(media.get("src", ""))
            if not source.startswith("./media/"):
                continue
            media_owner = source.split("/", 3)[2]
            if media_owner != sample["id"]:
                assert operation.get("genericMedia") is True


def test_web_app_readiness_accepts_worker_or_main_thread_runtime() -> None:
    calls: list[tuple[str, int]] = []

    class FakePage:
        def wait_for_function(self, script: str, *, timeout: int) -> None:
            calls.append((script, timeout))

    wait_for_web_app_ready(FakePage())

    assert len(calls) == 1
    script, timeout = calls[0]
    assert "app.diagramGenerationWorkerReady === true" in script
    assert "app.pyodideReady === true" in script
    assert "status.startsWith('Startup Error:')" in script
    assert timeout == 120000
