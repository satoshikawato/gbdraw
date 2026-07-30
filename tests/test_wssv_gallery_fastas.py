from __future__ import annotations

import base64
import hashlib
from pathlib import Path
import shutil
import subprocess

import pytest

from gbdraw.session_io import load_session, write_session_json
from tools.restore_wssv_gallery_fastas import (
    REPO_ROOT,
    SESSION_PATH,
    WSSV_FASTA_SOURCES,
    restore_wssv_gallery_fastas,
    validate_embedded_wssv_fastas,
)


def _embedded_fastas(
    session: dict[str, object],
) -> tuple[tuple[str, str, int, str], ...]:
    web_files = session["webFiles"]
    resources = session["resources"]
    assert isinstance(web_files, dict)
    assert isinstance(resources, dict)
    resource_ids = web_files["conservationLosatFastaSources"]
    assert isinstance(resource_ids, list)
    rows: list[tuple[str, str, int, str]] = []
    for resource_id in resource_ids:
        assert isinstance(resource_id, str)
        resource = resources[resource_id]
        assert isinstance(resource, dict)
        raw = base64.b64decode(resource["data"], validate=True)
        rows.append(
            (
                resource_id,
                str(resource["name"]),
                len(raw),
                hashlib.sha256(raw).hexdigest(),
            )
        )
    return tuple(rows)


def test_wssv_source_manifest_records_exact_identity_and_provenance() -> None:
    assert len(WSSV_FASTA_SOURCES) == 20
    assert len({source.label for source in WSSV_FASTA_SOURCES}) == 20
    assert len({source.record_id for source in WSSV_FASTA_SOURCES}) == 20
    assert all(source.length > 0 for source in WSSV_FASTA_SOURCES)
    assert all(len(source.sha256) == 64 for source in WSSV_FASTA_SOURCES)
    assert all(
        source.source_url or source.label == "Shantou2019"
        for source in WSSV_FASTA_SOURCES
    )
    assert sum(source.source == "NCBI Nucleotide" for source in WSSV_FASTA_SOURCES) == 12


def test_wssv_fastas_survive_session_load_save_and_generator_rebuild(
    tmp_path: Path,
) -> None:
    session = load_session(SESSION_PATH)
    assert validate_embedded_wssv_fastas(session) == WSSV_FASTA_SOURCES
    expected = _embedded_fastas(session)
    assert tuple(row[3] for row in expected) == tuple(
        source.sha256 for source in WSSV_FASTA_SOURCES
    )

    # Once restored, the Gallery generator can rebuild from the embedded
    # resources alone; no untracked source directory or network is required.
    assert not restore_wssv_gallery_fastas(
        session,
        ncbi_fasta=None,
        mag_dir=tmp_path / "absent",
    )

    round_trip_path = tmp_path / SESSION_PATH.name
    write_session_json(round_trip_path, session)
    round_tripped = load_session(round_trip_path)
    assert validate_embedded_wssv_fastas(round_tripped) == WSSV_FASTA_SOURCES
    assert _embedded_fastas(round_tripped) == expected


def test_wssv_web_projection_rebuilds_all_popup_sequence_sources() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the Web session projection test")
    subprocess.run(
        [node, "tests/web/wssv-gallery-fastas.test.mjs"],
        cwd=REPO_ROOT,
        check=True,
    )
