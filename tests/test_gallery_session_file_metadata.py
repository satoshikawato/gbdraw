"""Every checked-in Gallery Session is machine-independent (OV-2964): it stores, for each file, the
media type its resource kind has and ``lastModified`` 0, never a machine's media type guess or an
input's modification time. ``written_file_metadata`` is the one owner of that metadata."""

from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw.session_container import read_session_container
from gbdraw.session_io import written_file_metadata


REPO_ROOT = Path(__file__).resolve().parents[1]
SESSIONS = sorted(
    [
        *(REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions").glob("*.gbdraw"),
        *(REPO_ROOT / "tests" / "test_inputs").glob("*.gbdraw"),
    ]
)


def _read(path: Path) -> dict[str, object]:
    """The document of a Session file (``gbdraw-session.json``), whose resource descriptors hold the metadata."""

    return dict(read_session_container(path).document)


def machine_file_metadata(document: dict[str, object]) -> list[str]:
    """Where ``document`` stores file metadata that its writer would not: each resource, and each
    file binding (a leaf with a ``resourceId``), against the metadata of the resource's kind."""

    resources = document.get("resources")
    resources = resources if isinstance(resources, dict) else {}
    found: list[str] = []
    kinds: dict[str, object] = {}
    for resource_id, resource in resources.items():
        kinds[resource_id] = resource.get("kind")
        expected = written_file_metadata(resource.get("kind"))
        actual = {"type": resource.get("type"), "lastModified": resource.get("lastModified")}
        if actual != expected:
            found.append(f"resources.{resource_id}: {actual} (expected {expected})")

    def visit(value: object, path: str) -> None:
        if isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, f"{path}[{index}]")
        elif isinstance(value, dict):
            if "resourceId" in value and ("type" in value or "lastModified" in value):
                expected = written_file_metadata(kinds.get(str(value["resourceId"])))
                actual = {"type": value.get("type"), "lastModified": value.get("lastModified")}
                if actual != expected:
                    found.append(f"{path}: {actual} (expected {expected})")
            for key, item in value.items():
                visit(item, f"{path}.{key}" if path else str(key))

    visit(document.get("drawings"), "drawings")
    return found


def test_machine_metadata_is_found_in_a_resource_and_its_binding() -> None:
    document = {
        "resources": {
            "record-1-genbank": {"kind": "genbank", "type": "application/x-gameboy-rom", "lastModified": 1792769561077},
            "table": {"kind": "canonical-tsv", "type": "text/tab-separated-values", "lastModified": 0},
        },
        "drawings": [
            {"webFiles": {"bindings": {"linearSeqs": [{"gb": {
                "resourceId": "record-1-genbank", "type": "application/x-gameboy-rom", "lastModified": 1792769561077,
            }}]}}}
        ],
    }

    found = machine_file_metadata(document)

    assert [entry.split(":")[0] for entry in found] == [
        "resources.record-1-genbank",
        "drawings[0].webFiles.bindings.linearSeqs[0].gb",
    ]


def test_every_gallery_session_file_is_checked() -> None:
    assert len(SESSIONS) > 1


@pytest.mark.parametrize("path", SESSIONS, ids=lambda path: path.name)
def test_gallery_session_stores_no_machine_file_metadata(path: Path) -> None:
    assert machine_file_metadata(_read(path)) == []
