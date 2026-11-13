from __future__ import annotations

import base64
import copy
from collections import Counter
from collections.abc import Callable
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Any
import zipfile
from xml.etree import ElementTree as ET

import pytest
import gbdraw.web_support.feature_catalog as feature_catalog_module
from gbdraw.exceptions import ValidationError
from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.session_io import (
    CURRENT_FEATURE_CATALOG_SCHEMA,
    CURRENT_SESSION_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    load_session,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
from tests.utils.session_files import read_session_payload
from tests.utils.session_drawing_view import drawing_view, drawings_document_of_view, session_view
import tools.prepare_interactive_gallery_assets as gallery_assets_module
import tools.refresh_gallery_sessions as refresh_gallery_sessions_module
from tools.prepare_interactive_gallery_assets import (
    DECODED_METADATA_HARD_LIMIT,
    DECODED_METADATA_REGRESSION_CEILING,
    EXAMPLES,
    INTERACTIVE_SVG_HARD_LIMIT,
    _catalog_location_parts,
    _catalog_nucleotide_sequence,
    _interactive_svg_measurements,
    _migrate_legacy_multipart_feature_ids,
    _session_interactive_context,
    _session_result_svg,
    _validate_interactive_svg_size,
    _validate_interactive_orthogroup_payload,
    _validate_catalog_feature_ids,
    _validate_source_feature_ids,
)
from tools.refresh_gallery_sessions import (
    TEST_INPUT_SESSION_FILES,
    TEST_INPUT_SESSION_ROOT,
    VIBRIO_EXPANDED_HARD_LIMIT,
    VIBRIO_EXPANDED_REGRESSION_CEILING,
    VIBRIO_EXPECTED_RAW_PAIRS,
    VIBRIO_FILE_HARD_LIMIT,
    VIBRIO_FILE_REGRESSION_CEILING,
    VIBRIO_RAW_ENTRY_COUNT,
    _gallery_file_transaction,
    _public_gallery_session_files,
    _refresh_one_session,
    _reopen_for_publication,
    _refresh_session_paths,
    _session_artifact_measurements,
    _session_path,
    _validate_current_session_catalog_structure,
    _validate_gallery_session_inventory,
    _validate_staged_gallery_session as _validate_staged_gallery_session_49,
    stored_protein_result_differences,
)


pytestmark = pytest.mark.gallery

BUNDLED_REQUEST_SCHEMAS = frozenset({CANONICAL_REQUEST_SCHEMA})
BUNDLED_SESSION_VERSIONS = frozenset({CURRENT_SESSION_VERSION})


def test_default_refresh_inventory_covers_gallery_and_test_input_sessions() -> None:
    paths = _refresh_session_paths(None)

    assert len(paths) == 13
    assert len(_public_gallery_session_files()) == 12
    assert len(TEST_INPUT_SESSION_FILES) == 1
    assert {path.name for path in paths} >= set(TEST_INPUT_SESSION_FILES)


def test_gallery_file_transaction_restores_all_outputs_on_failure(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.svg"
    second = tmp_path / "second.json"
    created = tmp_path / "created.webp"
    first.write_text("first-original", encoding="utf-8")
    second.write_text("second-original", encoding="utf-8")

    with pytest.raises(RuntimeError, match="synthetic asset failure"):
        with _gallery_file_transaction((first, second, created)):
            first.write_text("first-partial", encoding="utf-8")
            second.write_text("second-partial", encoding="utf-8")
            created.write_bytes(b"partial")
            raise RuntimeError("synthetic asset failure")

    assert first.read_text(encoding="utf-8") == "first-original"
    assert second.read_text(encoding="utf-8") == "second-original"
    assert not created.exists()


def test_prepare_gallery_assets_help_imports_without_circular_dependency() -> None:
    repo_root = Path(__file__).parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            str(repo_root / "tools" / "prepare_interactive_gallery_assets.py"),
            "--help",
        ],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert "Regenerate interactive Gallery SVGs" in completed.stdout
    assert "circular import" not in completed.stderr.lower()


def test_session_path_names_the_gallery_session_file() -> None:
    path = _session_path("vibrio-harveyi-group-collinear")

    assert path.name == "vibrio-harveyi-group-collinear.gbdraw"


def test_gallery_artifact_limits_match_phase_two_targets() -> None:
    assert INTERACTIVE_SVG_HARD_LIMIT == 41_943_040
    assert DECODED_METADATA_HARD_LIMIT == 200_000_000
    assert DECODED_METADATA_HARD_LIMIT <= DECODED_METADATA_REGRESSION_CEILING
    assert DECODED_METADATA_REGRESSION_CEILING == 220_000_000
    assert VIBRIO_FILE_HARD_LIMIT == 90_000_000
    assert VIBRIO_FILE_HARD_LIMIT <= VIBRIO_FILE_REGRESSION_CEILING
    assert VIBRIO_FILE_REGRESSION_CEILING == 95_000_000
    assert VIBRIO_EXPANDED_HARD_LIMIT == 400_000_000
    assert (
        VIBRIO_EXPANDED_HARD_LIMIT
        <= VIBRIO_EXPANDED_REGRESSION_CEILING
    )
    assert VIBRIO_EXPANDED_REGRESSION_CEILING == 420_000_000


def test_interactive_svg_measurements_decode_compressed_schema_four() -> None:
    payload = {
        "schema": 4,
        "items": [
            {
                "resultIndex": 0,
                "resultName": "result",
                "recordKeys": ["record-key"],
                "features": [
                    {
                        "svgId": "rendered",
                        "recordKey": "record-key",
                        "biologicalFeatureId": "biological",
                    }
                ],
                "biologicalFeatures": [
                    {
                        "recordKey": "record-key",
                        "biologicalFeatureId": "biological",
                    }
                ],
                "orthogroups": [],
                "annotations": [],
                "comparisonMatches": [],
            }
        ],
    }
    decoded = json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    compressed = gzip.compress(decoded, compresslevel=9, mtime=0)
    output = (
        '<svg xmlns="http://www.w3.org/2000/svg">'
        '<metadata id="gbdraw-interactive-feature-metadata" '
        'data-schema="4" data-encoding="gzip-base64">'
        f"{base64.b64encode(compressed).decode('ascii')}"
        "</metadata></svg>"
    )

    measurements = _interactive_svg_measurements(output)

    assert measurements["totalInteractiveSvgBytes"] == len(
        output.encode("utf-8")
    )
    assert measurements["compressedMetadataBytes"] == len(compressed)
    assert measurements["decodedMetadataBytes"] == len(decoded)
    assert measurements["renderedFeatureCount"] == 1
    assert measurements["biologicalFeatureCount"] == 1
    assert _validate_interactive_svg_size(EXAMPLES[0], output) == measurements


def test_vibrio_gallery_interactive_svg_meets_regenerated_targets() -> None:
    example = next(
        example
        for example in EXAMPLES
        if example.id == "vibrio-harveyi-group-collinear"
    )
    output = example.gallery_svg_path.read_text(encoding="utf-8")

    measurements = _interactive_svg_measurements(output)

    assert measurements["totalInteractiveSvgBytes"] < INTERACTIVE_SVG_HARD_LIMIT
    assert (
        measurements["decodedMetadataBytes"]
        <= DECODED_METADATA_HARD_LIMIT
        <= DECODED_METADATA_REGRESSION_CEILING
    )
    assert (
        0
        < measurements["compressedMetadataBytes"]
        < measurements["decodedMetadataBytes"]
    )
    assert measurements["renderedFeatureCount"] == 9_375
    assert measurements["biologicalFeatureCount"] == 18_782


def test_session_artifact_measurements_report_component_bytes(
    tmp_path: Path,
) -> None:
    session = {
        "results": [{"name": "result", "content": "<svg/>"}],
        "resources": {"record": {"data": "QUJD", "encoding": "base64"}},
        "webFiles": {"files": []},
        "editorState": {
            "featureCatalog": {"schema": 4, "items": []},
            "legend": {"entries": []},
        },
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
    }
    document = json.dumps(session, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    artifact_path = tmp_path / "session.gbdraw"
    with zipfile.ZipFile(artifact_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("gbdraw-session.json", document)
        archive.writestr("resources/record", b"ABC")

    measurements = _session_artifact_measurements(session, artifact_path)

    assert measurements["sessionFileBytes"] == artifact_path.stat().st_size
    assert measurements["sessionExpandedBytes"] == len(document) + 3
    for key in (
        "resultsBytes",
        "resourcesBytes",
        "webFilesBytes",
        "featureCatalogBytes",
        "editorStateBytes",
        "losatCacheBytes",
        "losatDerivedCacheBytes",
    ):
        assert measurements[key] > 0


def _stored_protein_session(cli_args: list[str] | None) -> dict[str, object]:
    session: dict[str, object] = {
        "renderRequest": {
            "comparisons": [
                {"kind": "precomputedProteinComparison", "resourceId": "table-1"},
                {"kind": "orthogroupResult", "resourceId": "groups"},
                {"kind": "generatedProteinComparison", "mode": "none", "pairs": [], "settings": {}},
            ]
        },
        "resources": {"table-1": {"data": "x"}, "groups": {"data": "y"}, "record-1": {"data": "z"}},
        "modes": {"linear": {"config": {"losat": {"blastp": {"mode": "orthogroup"}}}}},
    }
    if cli_args is not None:
        session["cliInvocation"] = {"schema": 1, "mode": "linear", "args": cli_args}
    return session


def test_refresh_reopens_stored_protein_results_in_the_declared_mode() -> None:
    """A refresh computes stored protein results again (D-65, OV-1031, OV-1048).

    The stored comparisons and their resources leave the request, and the
    pipeline marker and the Web draft both take the mode of the recorded CLI
    invocation, so the Web rebuilds in the mode the Gallery declares.
    """

    session = _stored_protein_session(
        ["--losat", "losatp", "--losatp_mode", "collinear", "--gbk", "a.gb", "b.gb"]
    )
    gallery_path = Path("hepatoplasmataceae_collinear.gbdraw")

    assert _reopen_for_publication(gallery_path, session) == 1
    assert session["renderRequest"]["comparisons"] == [
        {"kind": "generatedProteinComparison", "mode": "collinear", "pairs": [], "settings": {}}
    ]
    assert set(session["resources"]) == {"record-1"}
    assert session["modes"]["linear"]["config"]["losat"]["blastp"]["mode"] == "collinear"
    assert _reopen_for_publication(gallery_path, session) == 0

    web_saved = _stored_protein_session(None)
    assert _reopen_for_publication(gallery_path, web_saved) == 1
    assert web_saved["renderRequest"]["comparisons"][0]["mode"] == "orthogroup"

    test_input = _stored_protein_session(None)
    assert _reopen_for_publication(TEST_INPUT_SESSION_ROOT / TEST_INPUT_SESSION_FILES[0], test_input) == 0

    # A public Gallery Session that shares the test input's file name is
    # refreshed like any other: the exemption is the test input folder's.
    same_name = _stored_protein_session(None)
    assert _reopen_for_publication(gallery_path.parent / TEST_INPUT_SESSION_FILES[0], same_name) == 1


@pytest.mark.slow
@pytest.mark.parametrize(
    "session_path",
    [example.session_path for example in EXAMPLES],
    ids=lambda path: path.name,
)
def test_gallery_protein_results_equal_a_rebuild_from_the_raw_cache(
    session_path: Path,
) -> None:
    """The first Generate draws what a changed setting would compute (OV-1031).

    The Web reuses a Session's stored protein results only while its filters
    are unchanged; any other value computes them from the raw LOSATP cache.
    A stored result older than the current code would make a filter round
    trip change the drawing.
    """

    session = session_view(read_session_payload(session_path))
    assert stored_protein_result_differences(session) == []


def test_current_session_catalog_structure_rejects_duplicate_payloads(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "duplicate.gbdraw"
    item = {
        "resultIndex": 0,
        "resultName": "result",
        "recordKeys": ["record-key"],
        "features": [
            {
                "svgId": "rendered",
                "recordKey": "record-key",
                "biologicalFeatureId": "biological",
            }
        ],
        "biologicalFeatures": [
            {
                "recordKey": "record-key",
                "biologicalFeatureId": "biological",
                "qualifiers": {"product": ["example"]},
                "nucleotide_sequence": "ATG",
            }
        ],
        "orthogroups": [],
        "annotations": [],
        "comparisonMatches": [],
    }
    session = {
        "resources": {
            "first": {"data": "QUJD", "encoding": "base64"},
            "second": {"data": "QUJD", "encoding": "base64"},
        },
        "editorState": {
            "featureCatalog": {"schema": 4, "items": [item]}
        },
    }

    with pytest.raises(ValueError, match="duplicates one resource payload"):
        _validate_current_session_catalog_structure(session_path, session)

    session["resources"].pop("second")
    item["features"][0]["qualifiers"] = {"product": ["duplicate"]}
    with pytest.raises(
        ValueError,
        match="rendered feature duplicates biological fields",
    ):
        _validate_current_session_catalog_structure(session_path, session)


def test_current_session_catalog_structure_accepts_drawn_selector(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "current.gbdraw"
    session = {
        "resources": {},
        "editorState": {
            "featureCatalog": {
                "schema": 5,
                "items": [
                    {
                        "resultIndex": 0,
                        "resultName": "result",
                        "recordKeys": ["record-key"],
                        "features": [
                            {
                                "svgId": "rendered",
                                "recordKey": "record-key",
                                "biologicalFeatureId": "biological",
                                "fillColor": "#000000",
                                "drawnSelector": {"hash": "abc"},
                            }
                        ],
                        "biologicalFeatures": [
                            {
                                "recordKey": "record-key",
                                "biologicalFeatureId": "biological",
                            }
                        ],
                        "orthogroups": [],
                        "annotations": [],
                        "comparisonMatches": [],
                    }
                ],
            }
        },
    }

    _validate_current_session_catalog_structure(session_path, session)


def test_vibrio_gallery_session_stores_resolved_comparisons_and_raw_cache(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    """Load and the first Generate reuse the stored result (Owner, 2026-10-09).

    The request carries the collinear result and the marker that reuses it, so
    the first Generate after Load does not compute the comparisons again. The
    raw LOSATP cache stays, so a changed setting recomputes without LOSATP.
    """

    path = _session_path("vibrio-harveyi-group-collinear")
    session = session_view(load_cached_gallery_session(path))

    request = session["renderRequest"]
    assert request["schema"] in BUNDLED_REQUEST_SCHEMAS
    comparisons = request["comparisons"]
    assert [comparison["kind"] for comparison in comparisons] == [
        "collinearityResult",
        "generatedProteinComparison",
    ]
    assert comparisons[1]["mode"] == "none"
    assert session["resources"][comparisons[0]["resourceId"]]["kind"] == (
        "collinearity-result"
    )
    assert session["modes"]["linear"]["config"]["losat"]["blastp"]["mode"] == (
        "collinear"
    )
    assert session["losatDerivedCache"]["entries"] == []

    protein_entries = [
        entry
        for entry in session["losatCache"]["entries"]
        if entry.get("schema") == PROTEIN_LOSAT_CACHE_SCHEMA
    ]
    assert len(protein_entries) == VIBRIO_RAW_ENTRY_COUNT
    assert {
        (
            entry["queryRecordInstanceKey"],
            entry["subjectRecordInstanceKey"],
        )
        for entry in protein_entries
    } == VIBRIO_EXPECTED_RAW_PAIRS


_GALLERY_ROOT = Path(__file__).parents[1] / "gbdraw" / "web" / "gallery"
_GALLERY_ENTRIES = json.loads(
    (_GALLERY_ROOT / "examples.json").read_text(encoding="utf-8")
)
_SEQUENCE_INPUT_SUFFIXES = frozenset(
    {".gb", ".gbk", ".gbff", ".gff", ".gff3", ".fa", ".fas", ".fasta", ".fna"}
)


def _declared_input_names(entry: dict[str, object]) -> set[str]:
    """The input files a Gallery entry declares: Files plus command inputs."""

    declared = {str(name) for name in entry["featureSources"]}
    args = shlex.split(str(entry["command"]))
    declared.update(
        Path(arg).name
        for arg in args
        if Path(arg).suffix.lower() in _SEQUENCE_INPUT_SUFFIXES
    )
    for flag, value in zip(args, args[1:]):
        if flag != "--records_table":
            continue
        table = Path(__file__).parents[1] / value
        with table.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                declared.update(
                    Path(row[column]).name
                    for column in ("gbk", "gff", "fasta")
                    if row.get(column)
                )
    return declared


def _record_source_original_names(session: dict[str, object]) -> list[str]:
    resources = session["resources"]
    original_names = (session.get("webFiles") or {}).get("resourceOriginalNames") or {}
    resource_ids: list[str] = []
    for record in session["renderRequest"]["records"]:
        for field in ("resourceId", "gffResourceId", "fastaResourceId"):
            resource_id = record["source"].get(field)
            if resource_id and resource_id not in resource_ids:
                resource_ids.append(resource_id)
    names = []
    for resource_id in resource_ids:
        name = original_names.get(resource_id)
        if name is None:
            name = str(resources[resource_id].get("name") or "")
            name = name.removeprefix(f"{resource_id}-")
        names.append(name)
    return names


@pytest.mark.parametrize(
    "entry", [pytest.param(entry, id=entry["id"]) for entry in _GALLERY_ENTRIES]
)
def test_gallery_session_record_sources_are_the_declared_input_files(
    entry: dict[str, object],
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    """One declared input file is one Session resource with its own name.

    A Session that stores each record of a multi-record file as its own
    resource names resources that no declared input has, or repeats one.
    """

    session = session_view(load_cached_gallery_session(_GALLERY_ROOT / str(entry["session"])))
    names = _record_source_original_names(session)

    assert names
    assert len(names) == len(set(names)), names
    assert set(names) <= _declared_input_names(entry), names


def test_gallery_session_inventory_matches_files_and_examples() -> None:
    _validate_gallery_session_inventory()


def test_all_bundled_sessions_use_supported_request_and_current_artifact_schemas(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    repo_root = Path(__file__).parents[1]
    paths = sorted(
        (repo_root / "gbdraw" / "web" / "gallery" / "sessions").glob("*.gbdraw")
    )
    paths.extend(sorted((repo_root / "tests" / "test_inputs").glob("*.gbdraw")))

    assert len(paths) == 13
    for path in paths:
        session = session_view(load_cached_gallery_session(path))
        assert session["version"] in BUNDLED_SESSION_VERSIONS, path
        assert session["renderRequest"]["schema"] in BUNDLED_REQUEST_SCHEMAS, path
        assert (
            session["proteinIdentityManifest"]["schema"]
            == PROTEIN_IDENTITY_MANIFEST_SCHEMA
        ), path
        assert not (
            session.get("legacyArtifacts", {})
            .get("proteinRawCandidates", {})
            .get("entries", [])
        ), path
        for entry in session.get("losatCache", {}).get("entries", []):
            if entry.get("identityKind") == "protein":
                assert entry["schema"] == PROTEIN_LOSAT_CACHE_SCHEMA, path
            else:
                assert entry["schema"] == NUCLEOTIDE_LOSAT_CACHE_SCHEMA, path
        assert all(
            entry["schema"] == LOSAT_DERIVED_CACHE_SCHEMA
            for entry in session.get("losatDerivedCache", {}).get("entries", [])
        ), path


def test_bundled_gallery_sources_match_current_session_results(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    for example in EXAMPLES:
        if not example.sync_result_svg:
            continue
        session = session_view(load_cached_gallery_session(example.session_path))
        assert (
            example.source_svg_path.read_text(encoding="utf-8")
            == _session_result_svg(session, example)
        ), example.id


def _validate_staged_gallery_session(session_path: Path, session: dict[str, object], **kwargs: Any) -> None:
    """The tool's validator on a staged Session 49: a one-drawing view is written as one by the production writer."""

    if session.get("version") == CURRENT_SESSION_VERSION:
        session = drawings_document_of_view({"createdAt": "2026-01-02T00:00:00Z", **session})
    _validate_staged_gallery_session_49(session_path, session, **kwargs)


def _staged_geometry_session(
    *,
    mode: str,
    records: list[dict[str, object]] | None,
    tracks: dict[str, object] | None = None,
) -> dict[str, object]:
    request: dict[str, object] = {
        "schema": CANONICAL_REQUEST_SCHEMA,
        "mode": mode,
    }
    options: dict[str, object] = {}
    if mode == "linear":
        options["config"] = {
            "canvas": {
                "linear": {
                    "vertical_padding": 8.0,
                }
            }
        }
    if tracks is not None:
        options["tracks"] = tracks
    if options:
        request["diagramOptions"] = options
    session: dict[str, object] = {
        "format": "gbdraw-session",
        "version": CURRENT_SESSION_VERSION,
        "renderRequest": request,
        "resources": {},
        "results": [{"name": "result", "content": "<svg></svg>"}],
        "editorState": {
            "featureCatalog": {
                "schema": CURRENT_FEATURE_CATALOG_SCHEMA,
                "items": [
                    {
                        "resultIndex": 0,
                        "resultName": "result",
                        "recordKeys": [],
                        "features": [],
                        "biologicalFeatures": [],
                        "orthogroups": [],
                        "annotations": [],
                        "comparisonMatches": [],
                    }
                ],
            }
        },
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
        "proteinIdentityManifest": {
            "schema": PROTEIN_IDENTITY_MANIFEST_SCHEMA,
            "proteinSets": {},
            "recordAnalyses": {},
            "recordInstances": {},
        },
    }
    if records is not None:
        session["runMetadata"] = {
            "trackSlotGeometry": {
                "schema": 1,
                "mode": mode,
                "source": "resolved",
                "records": records,
            }
        }
    return session


def _linear_geometry_record(
    *,
    feature_to_gc_gap: float = 0.0,
    gc_to_skew_gap: float = 0.0,
    feature_spacing: float = 0.0,
    gc_spacing: float = 0.0,
    skew_spacing: float = 0.0,
) -> dict[str, object]:
    local_top = -10.4
    local_bottom = 10.4
    feature_bottom = 12.5
    gc_origin = feature_bottom + feature_to_gc_gap - local_top
    gc_bottom = gc_origin + local_bottom
    skew_origin = gc_bottom + gc_to_skew_gap - local_top
    return {
        "recordIndex": 0,
        "recordId": "record",
        "slots": [
            {
                "slotIndex": 0,
                "slotId": "features",
                "renderer": "features",
                "side": "overlay",
                "heightPx": 25.0,
                "spacingAfterPx": feature_spacing,
                "baseYOffsetPx": 0.0,
                "resolvedOriginPx": 0.0,
                "dataAvailable": True,
                "paintBand": {"topPx": -12.0, "bottomPx": 12.0},
                "reserveBand": {"topPx": -12.5, "bottomPx": 12.5},
            },
            {
                "slotIndex": 1,
                "slotId": "gc_content",
                "renderer": "dinucleotide_content",
                "side": "below",
                "heightPx": 20.0,
                "spacingAfterPx": gc_spacing,
                "baseYOffsetPx": 18.0,
                "resolvedOriginPx": gc_origin,
                "dataAvailable": True,
                "paintBand": {
                    "topPx": gc_origin + local_top,
                    "bottomPx": gc_origin + local_bottom,
                },
                "reserveBand": {
                    "topPx": gc_origin + local_top,
                    "bottomPx": gc_origin + local_bottom,
                },
            },
            {
                "slotIndex": 2,
                "slotId": "gc_skew",
                "renderer": "dinucleotide_skew",
                "side": "below",
                "heightPx": 20.0,
                "spacingAfterPx": skew_spacing,
                "baseYOffsetPx": 36.0,
                "resolvedOriginPx": skew_origin,
                "dataAvailable": True,
                "paintBand": {
                    "topPx": skew_origin + local_top,
                    "bottomPx": skew_origin + local_bottom,
                },
                "reserveBand": {
                    "topPx": skew_origin + local_top,
                    "bottomPx": skew_origin + local_bottom,
                },
            },
        ],
    }


def test_staged_gallery_validator_requires_fresh_track_geometry(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "missing-geometry.gbdraw"
    session = _staged_geometry_session(mode="linear", records=None)

    with pytest.raises(ValueError, match="no resolved track-slot run metadata"):
        _validate_staged_gallery_session(
            session_path,
            session,
            require_track_geometry=True,
        )


def test_staged_gallery_validator_rejects_circular_feature_lane_regression(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "circular-geometry.gbdraw"
    tracks = {
        "circularTrackSlots": [
            "features:features@lane_direction=split",
        ],
        "circularTrackAxisIndex": 0,
    }

    valid_record = {
        "recordIndex": 0,
        "recordId": "record",
        "axisRadiusPx": 100.0,
        "slots": [
            {
                "slotIndex": 0,
                "slotId": "features",
                "renderer": "features",
                "side": "overlay",
                "widthPx": 20.0,
                "radiusFactor": 1.0,
            }
        ],
    }
    _validate_staged_gallery_session(
        session_path,
        _staged_geometry_session(
            mode="circular",
            records=[valid_record],
            tracks=tracks,
        ),
        require_track_geometry=True,
    )

    inside_record = {
        **valid_record,
        "slots": [
            {
                **valid_record["slots"][0],
                "radiusFactor": 0.85,
            }
        ],
    }
    with pytest.raises(ValueError, match="circular Feature lane geometry"):
        _validate_staged_gallery_session(
            session_path,
            _staged_geometry_session(
                mode="circular",
                records=[inside_record],
                tracks=tracks,
            ),
            require_track_geometry=True,
        )


def test_staged_gallery_validator_rejects_linear_spacing_regression(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "linear-geometry.gbdraw"
    _validate_staged_gallery_session(
        session_path,
        _staged_geometry_session(
            mode="linear",
            records=[_linear_geometry_record()],
        ),
        require_track_geometry=True,
    )

    for corrupt_record in (
        _linear_geometry_record(feature_to_gc_gap=23.6),
        _linear_geometry_record(gc_to_skew_gap=12.0),
    ):
        with pytest.raises(ValueError, match="default adjacent reserve gap"):
            _validate_staged_gallery_session(
                session_path,
                _staged_geometry_session(
                    mode="linear",
                    records=[corrupt_record],
                ),
                require_track_geometry=True,
            )

    corrupt_metadata = _linear_geometry_record(
        feature_to_gc_gap=23.6,
        gc_to_skew_gap=12.0,
    )
    corrupt_metadata["slots"][0]["spacingAfterPx"] = 23.6
    corrupt_metadata["slots"][1]["spacingAfterPx"] = 12.0
    with pytest.raises(ValueError, match="default spacing metadata"):
        _validate_staged_gallery_session(
            session_path,
            _staged_geometry_session(
                mode="linear",
                records=[corrupt_metadata],
            ),
            require_track_geometry=True,
        )


def test_staged_gallery_validator_skips_empty_feature_before_numeric_tracks(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "linear-featureless-geometry.gbdraw"
    valid_record = _linear_geometry_record(feature_to_gc_gap=9.0)
    feature = valid_record["slots"][0]
    assert isinstance(feature, dict)
    feature.update(
        {
            "dataAvailable": False,
            "paintBand": None,
            "reserveBand": {"topPx": 0.0, "bottomPx": 0.0},
        }
    )

    _validate_staged_gallery_session(
        session_path,
        _staged_geometry_session(
            mode="linear",
            records=[valid_record],
        ),
        require_track_geometry=True,
    )

    invalid_record = _linear_geometry_record(
        feature_to_gc_gap=9.0,
        gc_to_skew_gap=12.0,
    )
    invalid_feature = invalid_record["slots"][0]
    assert isinstance(invalid_feature, dict)
    invalid_feature.update(feature)
    with pytest.raises(ValueError, match="default adjacent reserve gap"):
        _validate_staged_gallery_session(
            session_path,
            _staged_geometry_session(
                mode="linear",
                records=[invalid_record],
            ),
            require_track_geometry=True,
        )


def test_staged_gallery_validator_allows_custom_extra_spacing(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "linear-custom-geometry.gbdraw"
    _validate_staged_gallery_session(
        session_path,
        _staged_geometry_session(
            mode="linear",
            records=[
                _linear_geometry_record(
                    feature_to_gc_gap=10.0,
                    gc_to_skew_gap=8.0,
                    feature_spacing=4.0,
                    gc_spacing=4.0,
                )
            ],
            tracks={"linearTrackSlots": []},
        ),
        require_track_geometry=True,
    )

    with pytest.raises(ValueError, match="overlap or violate declared spacing"):
        _validate_staged_gallery_session(
            session_path,
            _staged_geometry_session(
                mode="linear",
                records=[
                    _linear_geometry_record(
                        gc_to_skew_gap=2.0,
                        feature_spacing=4.0,
                        gc_spacing=4.0,
                    )
                ],
                tracks={"linearTrackSlots": []},
            ),
            require_track_geometry=True,
        )


def test_refresh_leaves_a_test_input_pipeline_unresolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The Generate-pipeline specs load this input, so its orthogroup pipeline
    # stays; only public Gallery Sessions store resolved comparisons.
    name = TEST_INPUT_SESSION_FILES[0]
    monkeypatch.setattr(refresh_gallery_sessions_module, "TEST_INPUT_SESSION_ROOT", tmp_path)
    source = tmp_path / name
    destination = tmp_path / "refreshed" / name
    destination.parent.mkdir()
    source.write_bytes((TEST_INPUT_SESSION_ROOT / name).read_bytes())

    def pipelines(path: Path) -> list[tuple[str, object]]:
        return [
            (item["kind"], item.get("mode"))
            for item in _published_drawing(path)["renderRequest"]["comparisons"]
        ]

    before = pipelines(source)
    assert ("generatedProteinComparison", "orthogroup") in before

    _refresh_one_session(source, destination_path=destination)

    assert pipelines(destination) == before


def test_refresh_records_resolved_track_geometry(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.gbdraw"
    destination = tmp_path / "refreshed.gbdraw"
    source.write_bytes(
        _session_path("HmmtDNA_basic_circular").read_bytes()
    )

    _refresh_one_session(source, destination_path=destination)

    refreshed = _published_drawing(destination)
    assert refreshed["version"] == CURRENT_SESSION_VERSION
    assert refreshed["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert refreshed["renderRequest"]["output"]["interactiveMetadataPolicy"] == "auto"
    assert (
        "outputPrefix"
        not in refreshed["renderRequest"]["diagramOptions"].get("output", {})
    )
    geometry = refreshed["runMetadata"]["trackSlotGeometry"]
    assert geometry["schema"] == 1
    assert geometry["mode"] == "circular"
    assert geometry["source"] == "resolved"
    assert geometry["records"][0]["axisRadiusPx"] > 0


@pytest.mark.parametrize(
    "example",
    [
        pytest.param(example, id=example.id)
        for example in gallery_assets_module.EXAMPLES
        if "Multi-record" in example.tags
    ],
)
def test_multi_record_gallery_sessions_are_built_from_their_declared_commands(
    example: gallery_assets_module.GallerySessionExample,
) -> None:
    """Refresh rebuilds a multi-record example from its command and inputs.

    Replaying a stored Session keeps whatever resources it holds, so a Session
    once split per record would stay split. The declared command reads each
    input file as one source, and its inputs resolve in ``command_dir``.
    """

    declared = refresh_gallery_sessions_module._declared_session_example(
        example.session_path
    )
    assert declared is example
    args = shlex.split(example.command)
    inputs = [
        value
        for flag, value in zip(args, args[1:])
        if flag in {"--gbk", "--gff", "--fasta", "--records_table"}
    ]
    assert inputs
    for value in inputs:
        assert (example.command_cwd / value).is_file(), value


def _declared_input_files(example: gallery_assets_module.GallerySessionExample) -> list[Path]:
    args = shlex.split(example.command)
    files = []
    for flag, value in zip(args, args[1:]):
        if flag in {"--gbk", "--gff", "--fasta"}:
            files.append(example.command_cwd / value)
        elif flag == "--records_table":
            table = example.command_cwd / value
            with table.open(encoding="utf-8", newline="") as handle:
                files.extend(
                    table.parent / row[column]
                    for row in csv.DictReader(handle, delimiter="\t")
                    for column in ("gbk", "gff", "fasta")
                    if row.get(column)
                )
    return files


@pytest.mark.parametrize(
    "example",
    [
        pytest.param(example, id=example.id)
        for example in gallery_assets_module.EXAMPLES
        if example.session_from_command
    ],
)
def test_declared_command_gallery_sessions_store_their_input_files_unchanged(
    example: gallery_assets_module.GallerySessionExample,
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    """Each record source holds its declared input file byte for byte.

    The CLI stores an input file whose records it draws unchanged as that
    file, so Web Load binds the same bytes the Result was drawn from.
    """

    declared = {
        hashlib.sha256(path.read_bytes()).hexdigest(): path.name
        for path in _declared_input_files(example)
    }
    session = session_view(load_cached_gallery_session(example.session_path))
    resource_ids = {
        record["source"][field]
        for record in session["renderRequest"]["records"]
        for field in ("resourceId", "gffResourceId", "fastaResourceId")
        if record["source"].get(field)
    }
    stored = {
        resource_id: declared.get(
            hashlib.sha256(base64.b64decode(session["resources"][resource_id]["data"])).hexdigest()
        )
        for resource_id in resource_ids
    }
    assert None not in stored.values(), stored
    assert sorted(stored.values()) == sorted(set(declared.values()))


def test_replayed_gallery_sessions_have_no_declared_command_example() -> None:
    assert (
        refresh_gallery_sessions_module._declared_session_example(
            _session_path("lambda_basic_linear")
        )
        is None
    )


@pytest.mark.linear
def test_declared_command_refresh_keeps_one_file_as_one_resource(
    tmp_path: Path,
) -> None:
    """A Gallery Session built from a records table keeps two Files.

    Each multi-record GenBank file is one resource named after it, and its
    records are bound by record ID on the rows the table declares.
    """

    from tests.test_losatp_source_batching import _fake_losat

    inputs = Path(__file__).parent / "test_inputs"
    for name, members in (
        ("left.gbk", ("BGC0000708.gbk", "BGC0000709.gbk")),
        ("right.gbk", ("BGC0000711.gbk", "BGC0000712.gbk")),
    ):
        (tmp_path / name).write_text(
            "".join((inputs / member).read_text(encoding="utf-8") for member in members),
            encoding="utf-8",
        )
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\ttopology\tdisplay_start\trow\n"
        "left.gbk\tBGC0000708\tcircular\t1001\t1\nleft.gbk\tBGC0000709\t\t\t1\n"
        "right.gbk\tBGC0000711\t\t\t2\nright.gbk\tBGC0000712\t\t\t2\n",
        encoding="utf-8",
    )
    # A headerless colour table and a label table, which the Vibrio command no
    # longer declares but the Session builder must still read as one resource each.
    colors = tmp_path / "colors.tsv"
    colors.write_text("CDS\tprotein_id\tCAG38690\t#6a3d9a\tFirst CDS\n", encoding="utf-8")
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "record\tfeature_selector\tlabel_visibility\tlabel_text\n"
        "BGC0000708\tprotein_id=CAG38690.1\ton\tFirst\n",
        encoding="utf-8",
    )
    from tools.prepare_interactive_gallery_assets import (
        VIBRIO_HARVEYI_GROUP_COMMAND,
    )

    # The Vibrio command with these inputs and a fake LOSAT runtime.
    losat, _log = _fake_losat(tmp_path)
    argv = shlex.split(VIBRIO_HARVEYI_GROUP_COMMAND)
    for flag, value in (
        ("--records_table", str(table)),
        ("--losat_threads", "1"),
        ("-o", "declared"),
    ):
        argv[argv.index(flag) + 1] = value
    command = shlex.join(
        [
            *argv,
            "-t",
            str(colors),
            "--feature_override_table",
            str(labels),
            "--losat_bin",
            losat,
        ]
    )
    declared = tmp_path / "stage" / "declared.gbdraw"
    declared.parent.mkdir()
    destination = tmp_path / "published.gbdraw"

    refresh_gallery_sessions_module._declared_command_session(
        command, declared, env=refresh_gallery_sessions_module._cli_env()
    )
    _refresh_one_session(
        declared, destination_path=destination, declared_command=True
    )
    refresh_gallery_sessions_module._assert_declared_figure(declared, destination)

    session = _published_drawing(destination)
    assert session["cliInvocation"]["args"] == shlex.split(command)[2:]
    assert _stored_output_prefixes(session) == ("declared", "declared")
    # The Session stores the collinear result its command computed, with the
    # marker that reuses it, as the Web writes it after a Generate.
    comparisons = session["renderRequest"]["comparisons"]
    assert [(item["kind"], item.get("mode")) for item in comparisons] == [
        ("collinearityResult", None),
        ("generatedProteinComparison", "none"),
    ]
    # The Web's configuration overrides, so Session Load needs no Worker.
    options = session["renderRequest"]["diagramOptions"]
    assert options.get("config") is None
    assert options["configOverrides"]["objects.scale.interval"] == 750000
    # The tables publish as the declared command wrote them (OV-215, OV-216).
    assert [
        (row["recordKey"], row["labelVisibility"], row["labelText"])
        for row in options["featureOverrides"]
    ] == [("record-1", "on", "First")]
    color_table = options["colors"]["colorTable"] or options["colors"]["colorTableFile"]
    assert "CAG38690\t#6a3d9a\tFirst CDS" in base64.b64decode(
        session["resources"][color_table["resourceId"]]["data"]
    ).decode("utf-8")
    records = session["renderRequest"]["records"]
    assert records[0]["display"] == {"isCircular": True, "startCoordinate": 1001}
    assert [
        (
            record["source"]["resourceId"],
            record["selector"],
            record["presentation"]["gridRow"],
        )
        for record in records
    ] == [
        ("record-1-genbank", {"kind": "recordId", "value": "BGC0000708"}, 1),
        ("record-1-genbank", {"kind": "recordId", "value": "BGC0000709"}, 1),
        ("record-3-genbank", {"kind": "recordId", "value": "BGC0000711"}, 2),
        ("record-3-genbank", {"kind": "recordId", "value": "BGC0000712"}, 2),
    ]
    assert _record_source_original_names(session) == ["left.gbk", "right.gbk"]
    genbank = [
        base64.b64decode(resource["data"]).decode("utf-8")
        for resource in session["resources"].values()
        if resource["kind"] == "genbank"
    ]
    assert [text.count("\nLOCUS ") + text.startswith("LOCUS ") for text in genbank] == [2, 2]


@pytest.mark.circular
def test_declared_circular_command_refresh_keeps_one_file_as_one_resource(
    tmp_path: Path,
) -> None:
    """The Vnig command builds a Session with one resource for its one file.

    The command names its input by file name, so it runs in the directory that
    holds the input; the published Session binds both records of the file by
    record ID and draws the figure of the declared command.
    """

    from tools.prepare_interactive_gallery_assets import VNIG_COMMAND

    source = Path(__file__).parent / "fixtures" / "web_batch_two_records.gb"
    (tmp_path / "two-records.gbff").write_bytes(source.read_bytes())
    argv = shlex.split(VNIG_COMMAND)
    argv[argv.index("--gbk") + 1] = "two-records.gbff"
    argv[argv.index("-o") + 1] = "declared"
    # Keep the positions of the two records the file holds.
    for position in ("#3@2", "#4@2", "#5@2", "#6@2"):
        index = argv.index(position)
        del argv[index - 1 : index + 1]
    command = shlex.join(argv)
    declared = tmp_path / "stage" / "declared.gbdraw"
    declared.parent.mkdir()
    destination = tmp_path / "published.gbdraw"

    refresh_gallery_sessions_module._declared_command_session(
        command, declared, env=refresh_gallery_sessions_module._cli_env(), cwd=tmp_path
    )
    _refresh_one_session(
        declared, destination_path=destination, declared_command=True
    )
    refresh_gallery_sessions_module._assert_declared_figure(declared, destination)

    session = _published_drawing(destination)
    assert session["cliInvocation"]["args"] == argv[2:]
    assert _stored_output_prefixes(session) == ("declared", "declared")
    assert [
        (record["source"]["resourceId"], record["selector"])
        for record in session["renderRequest"]["records"]
    ] == [
        ("record-1-genbank", {"kind": "recordId", "value": "TESTA"}),
        ("record-1-genbank", {"kind": "recordId", "value": "TESTB"}),
    ]
    assert _record_source_original_names(session) == ["two-records.gbff"]
    assert [
        resource_id
        for resource_id, resource in session["resources"].items()
        if resource["kind"] == "genbank"
    ] == ["record-1-genbank"]


def test_gallery_refresh_refuses_a_session_27_44(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    # Gallery publication admits Sessions 46-50; Python converts a Session
    # 27-44 when the app loads it.
    source = tmp_path / "source.gbdraw"
    destination = tmp_path / "refreshed.gbdraw"
    source.write_bytes(
        (Path(__file__).parent / "fixtures/sessions/lambda_basic_linear.v40-schema5.json").read_bytes()
    )

    with pytest.raises(subprocess.CalledProcessError):
        _refresh_one_session(source, destination_path=destination)

    assert "Gallery publication supports Sessions 46-51" in capfd.readouterr().err
    assert not destination.exists()


def test_staged_gallery_validator_accepts_current_artifact_schemas(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "synthetic.gbdraw"
    session = {
        "format": "gbdraw-session",
        "version": CURRENT_SESSION_VERSION,
        "renderRequest": {"schema": CANONICAL_REQUEST_SCHEMA, "mode": "linear"},
        "resources": {},
        "results": [{"name": "result", "content": "<svg></svg>"}],
        "editorState": {
            "featureCatalog": {
                "schema": CURRENT_FEATURE_CATALOG_SCHEMA,
                "items": [
                    {
                        "resultIndex": 0,
                        "resultName": "result",
                        "recordKeys": [],
                        "features": [],
                        "biologicalFeatures": [],
                        "orthogroups": [],
                        "annotations": [],
                        "comparisonMatches": [],
                    }
                ],
            }
        },
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
        "proteinIdentityManifest": {
            "schema": 2,
            "proteinSets": {},
            "recordAnalyses": {},
            "recordInstances": {},
        },
    }

    _validate_staged_gallery_session(session_path, session)

    # The catalog stores no sequence (D5): the Session check refuses a Gallery Session that stores one.
    for key, value in (("nucleotide_sequence", "ATG"), ("sequenceSourceIndex", 0), ("sequenceSources", [])):
        stored = json.loads(json.dumps(session))
        item = stored["editorState"]["featureCatalog"]["items"][0]
        if key == "sequenceSources":
            item[key] = value
        else:
            item["biologicalFeatures"] = [{"recordKey": "record-1", "biologicalFeatureId": "f1", key: value}]
        with pytest.raises(ValidationError, match="schema 6 stores no sequences"):
            _validate_staged_gallery_session(session_path, stored)

    # A Gallery Session stores its protein comparison results; a pipeline that
    # the first Generate would run again is rejected (perf-014x D-17 /
    # docs-restructure D-20). A test input keeps its pipeline.
    def with_protein_mode(mode: str) -> dict[str, object]:
        request = dict(
            session["renderRequest"],
            comparisons=[
                {
                    "kind": "generatedProteinComparison",
                    "mode": mode,
                    "pairs": [],
                    "settings": {},
                }
            ],
        )
        return dict(session, renderRequest=request)

    _validate_staged_gallery_session(session_path, with_protein_mode("none"))
    with pytest.raises(ValueError, match="resolved protein comparison results"):
        _validate_staged_gallery_session(session_path, with_protein_mode("collinear"))
    _validate_staged_gallery_session(
        TEST_INPUT_SESSION_ROOT / TEST_INPUT_SESSION_FILES[0],
        with_protein_mode("orthogroup"),
    )
    with pytest.raises(ValueError, match="resolved protein comparison results"):
        _validate_staged_gallery_session(
            session_path.parent / TEST_INPUT_SESSION_FILES[0],
            with_protein_mode("orthogroup"),
        )
    _validate_staged_gallery_session(
        session_path, with_protein_mode("none")
    )

    stale_version = dict(session, version=33)
    with pytest.raises(ValueError, match=f"expected {CURRENT_SESSION_VERSION}"):
        _validate_staged_gallery_session(session_path, stale_version)

    stale_reference = dict(session, orthogroupState={"proteinId": "p_r_old"})
    with pytest.raises(ValueError, match="legacy protein identifiers"):
        _validate_staged_gallery_session(session_path, stale_reference)


def test_gallery_session_refresh_does_not_partially_replace_on_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.gbdraw"
    second = tmp_path / "second.gbdraw"
    first_original = '{"createdAt":"first-original"}'
    second_original = '{"createdAt":"second-original"}'
    first.write_text(first_original, encoding="utf-8")
    second.write_text(second_original, encoding="utf-8")

    monkeypatch.setattr(refresh_gallery_sessions_module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(refresh_gallery_sessions_module, "SESSION_ROOT", tmp_path)
    monkeypatch.setattr(
        refresh_gallery_sessions_module,
        "_session_path",
        lambda name: tmp_path / name,
    )
    monkeypatch.setattr(
        refresh_gallery_sessions_module,
        "load_session",
        lambda path: {"path": str(path)},
    )

    def fake_refresh(session_path: Path, *, destination_path: Path | None = None) -> None:
        assert destination_path is not None
        if session_path == second:
            raise RuntimeError("synthetic render failure")
        destination_path.write_text('{"createdAt":"first-refreshed"}', encoding="utf-8")

    monkeypatch.setattr(
        refresh_gallery_sessions_module,
        "_refresh_one_session",
        fake_refresh,
    )
    monkeypatch.setattr(
        refresh_gallery_sessions_module,
        "_validate_staged_gallery_session",
        lambda *_args, **_kwargs: None,
    )

    with pytest.raises(RuntimeError, match="synthetic render failure"):
        refresh_gallery_sessions_module.refresh_gallery_sessions(
            (first.name, second.name)
        )

    assert first.read_text(encoding="utf-8") == first_original
    assert second.read_text(encoding="utf-8") == second_original


def test_prepare_gallery_assets_preserves_existing_source_svgs() -> None:
    source = Path("tools/prepare_interactive_gallery_assets.py").read_text(
        encoding="utf-8"
    )

    assert "def _read_or_create_source_svg(" in source
    assert "existing_source = (" in source
    assert "elif existing_source is not None:" in source
    assert "if migrated != existing_source:" in source
    assert "example.source_svg_path.write_text(migrated" in source
    assert "_validate_source_feature_ids(example, session, migrated)" in source
    assert "def _sync_session_result_svg(" in source
    assert "_gallery_document(document, session, _gallery_view(document)[1])" in source
    assert "_sync_session_result_svg(example, session, source)" in source
    assert "_write_gallery_svg(example, session, source)" in source
    assert 'entry["tutorial"] = f"./tutorials/{example.id}.json"' in source
    assert 'entry["tutorialStatus"] = "ready"' in source


def test_gallery_result_sync_updates_track_geometry_result_name(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    example = EXAMPLES[0]
    session = {
        "title": "out",
        "resources": {},
        "renderRequest": {"mode": "circular"},
        "results": [{"name": "out", "content": "<svg id='old'/>"}],
        "runMetadata": {
            "trackSlotGeometry": {
                "records": [
                    {"resultIndex": 0, "resultName": "out"},
                    {"resultIndex": 1, "resultName": "other"},
                ]
            }
        },
    }
    writes: list[Path] = []
    monkeypatch.setattr(
        refresh_gallery_sessions_module,
        "_write_session_file",
        lambda path, _payload: writes.append(path),
    )

    gallery_assets_module._sync_session_result_svg(
        example,
        session,
        "<svg id='fresh'/>",
    )

    assert session["title"] == example.id
    assert session["results"][0] == {
        "name": example.id,
        "content": "<svg id='fresh'/>",
    }
    records = session["runMetadata"]["trackSlotGeometry"]["records"]
    assert records[0]["resultName"] == example.id
    assert records[1]["resultName"] == "other"
    assert writes == [example.session_path]


def test_gallery_source_refresh_replaces_existing_svg_from_session(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(gallery_assets_module, "SOURCE_ROOT", tmp_path)
    example = gallery_assets_module.GallerySessionExample(
        id="refresh-source",
        title="Refresh source",
        tags=(),
        description="",
        workflow="",
        input_summary="",
        display_order=1,
        command_kind="",
        command_note="",
        interactive_svg=False,
    )
    stale = '<svg xmlns="http://www.w3.org/2000/svg"><text>stale</text></svg>'
    current = '<svg xmlns="http://www.w3.org/2000/svg"><text>current</text></svg>'
    example.source_svg_path.write_text(stale, encoding="utf-8")

    source = gallery_assets_module._read_or_create_source_svg(
        example,
        {"results": [{"content": current}]},
        refresh_from_session=True,
    )

    assert source == current
    assert example.source_svg_path.read_text(encoding="utf-8") == current


def test_gallery_source_migrates_legacy_multipart_feature_ids() -> None:
    feature = {
        "svg_id": "fcurrent_record_1",
        "stable_svg_id": "fcurrent",
        "record_id": "record-1",
        "type": "CDS",
        "start": 10,
        "end": 40,
        "qualifiers": {},
        "nucleotide_sequence": "ATGC",
        "amino_acid_sequence": "M",
        "location_parts": [
            {"start": 10, "end": 20, "strand": "+"},
            {"start": 30, "end": 40, "strand": "+"},
        ],
    }
    legacy_id = compute_feature_hash_from_parts(
        "CDS",
        10,
        20,
        1,
        record_id="record-1",
    )
    session = {"features": {"extractedFeatures": [feature]}}
    source = (
        f'<svg><path id="{legacy_id}__part1" '
        f'data-gbdraw-feature-id="{legacy_id}" /></svg>'
    )

    migrated = _migrate_legacy_multipart_feature_ids(source, session)

    assert legacy_id not in migrated
    assert 'id="fcurrent__part1"' in migrated
    assert 'data-gbdraw-feature-id="fcurrent"' in migrated

    _validate_source_feature_ids(EXAMPLES[0], session, migrated)
    _validate_source_feature_ids(
        EXAMPLES[0],
        session,
        migrated.replace("fcurrent", "fcurrent_record_1"),
    )
    with pytest.raises(ValueError, match="without session metadata"):
        _validate_source_feature_ids(
            EXAMPLES[0], session, source.replace(legacy_id, "forphan")
        )


def test_gallery_source_rejects_feature_metadata_without_popup_details() -> None:
    session = {
        "features": {
            "extractedFeatures": [
                {"svg_id": "fcurrent_record_1", "stable_svg_id": "fcurrent"}
            ]
        }
    }
    source = (
        '<svg><path id="fcurrent_record_1" '
        'data-gbdraw-feature-id="fcurrent_record_1" /></svg>'
    )

    with pytest.raises(ValueError, match="without popup details"):
        _validate_source_feature_ids(EXAMPLES[0], session, source)


def _export_profile(item: dict[str, object]) -> dict[str, object]:
    """A schema-4 item with its sequences, as the export profile (schema 5) that the validators read."""

    return feature_catalog_module.promote_legacy_feature_catalog(
        {"schema": 4, "items": [item]},
        schema=feature_catalog_module.INTERACTIVE_EMBED_SCHEMA,
    )["items"][0]


def test_gallery_validators_materialize_catalog_sequence_references() -> None:
    feature = {
        "recordKey": "record-1",
        "biologicalFeatureId": "f-iupac",
        "record_id": "record-1",
        "type": "CDS",
        "start": 0,
        "end": 15,
        "location_parts": [
            {"start": 0, "end": 4, "strand": "+"},
            {"start": 4, "end": 15, "strand": "-"},
        ],
        "qualifiers": {},
        "sequenceSourceIndex": 1,
        "anchorProfile": {
            "precision": "unavailable",
            "operator": "unknown",
            "partOrder": "ambiguous",
            "strand": "mixed",
        },
    }
    item = {
        "resultIndex": 0,
        "resultName": "result",
        "recordKeys": ["record-1"],
        "features": [
            {
                "svgId": "f-iupac_record_1",
                "recordKey": "record-1",
                "biologicalFeatureId": "f-iupac",
            }
        ],
        "biologicalFeatures": [feature],
        "orthogroups": [
            {
                "id": "og-1",
                "members": [
                    {
                        "recordKey": "record-1",
                        "biologicalFeatureId": "f-iupac",
                    }
                ],
            }
        ],
        "annotations": [],
        "comparisonMatches": [],
        "sequenceSources": [
            {
                "key": "linear:record:0:reference",
                "origin": "linear-record",
                "recordIndex": 0,
                "sequence": "AAGCTTTTTTTTTTT",
            },
            {
                "key": "linear:record:0:iupac",
                "origin": "linear-record",
                "recordIndex": 0,
                "sequence": "ACGTRYSWKMBDHVN",
            },
        ],
    }
    source = (
        '<svg><path id="f-iupac_record_1" '
        'data-gbdraw-feature-id="f-iupac_record_1" /></svg>'
    )
    item = _export_profile(item)

    assert _catalog_nucleotide_sequence(feature, item) == "ACGTNBDHVKMWSRY"
    assert _catalog_nucleotide_sequence(
        {
            "recordKey": "record-1",
            "start": 0,
            "end": 4,
            "strand": "-",
            "sequenceSourceIndex": 0,
        },
        item,
    ) == "GCTT"
    assert _catalog_location_parts(
        {"start": 0, "end": 4, "strand": "-"}
    ) == [
        {"start": 0, "end": 4, "strand": "-", "display": "1..4"}
    ]
    assert _catalog_nucleotide_sequence(
        {
            "nucleotide_sequence": "ATG",
            "sequenceSourceIndex": 99,
        },
        item,
    ) == ""
    # The validators read a catalog item in the export profile, as the producer embeds it (D5).
    _validate_catalog_feature_ids(EXAMPLES[0], item, source)
    _validate_interactive_orthogroup_payload(EXAMPLES[0], item)

    simple_item = json.loads(json.dumps(item))
    simple_feature = simple_item["biologicalFeatures"][0]
    simple_feature.pop("location_parts")
    simple_feature["strand"] = "+"
    _validate_catalog_feature_ids(EXAMPLES[0], simple_item, source)

    broken_item = json.loads(json.dumps(item))
    broken_item["biologicalFeatures"][0]["location_parts"][1]["end"] = 16
    with pytest.raises(
        ValueError,
        match="invalid nucleotide sequence source reference",
    ):
        _validate_catalog_feature_ids(EXAMPLES[0], broken_item, source)
    with pytest.raises(ValueError, match="has no nucleotide sequence"):
        _validate_interactive_orthogroup_payload(EXAMPLES[0], broken_item)

    unreferenced_invalid_item = json.loads(json.dumps(item))
    unreferenced_invalid_item["sequenceSources"].append(
        {
            "key": "linear:record:invalid",
            "origin": "linear-record",
            "recordIndex": 0,
            "sequence": "ACGT RYSW",
        }
    )
    with pytest.raises(
        ValueError,
        match="invalid nucleotide sequence source reference",
    ):
        _validate_catalog_feature_ids(
            EXAMPLES[0],
            unreferenced_invalid_item,
            source,
        )
    with pytest.raises(
        ValueError,
        match="invalid nucleotide sequence source",
    ):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            unreferenced_invalid_item,
        )


def test_gallery_sequence_source_validation_is_bounded_by_source_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    feature_count = 256
    source_sequence = "ATG" + ("A" * 100_000)
    biological_features = [
        {
            "recordKey": "record-1",
            "biologicalFeatureId": f"feature-{index}",
            "record_id": "record-1",
            "type": "CDS",
            "start": 0,
            "end": 3,
            "strand": "+",
            "qualifiers": {},
            "sequenceSourceIndex": 0,
            "anchorProfile": {
                "precision": "exact",
                "operator": "single",
                "partOrder": "biological",
                "strand": "+",
            },
        }
        for index in range(feature_count)
    ]
    rendered_features = [
        {
            "svgId": f"rendered-{index}",
            "recordKey": "record-1",
            "biologicalFeatureId": f"feature-{index}",
        }
        for index in range(feature_count)
    ]
    item = {
        "resultIndex": 0,
        "resultName": "result",
        "recordKeys": ["record-1"],
        "features": rendered_features,
        "biologicalFeatures": biological_features,
        "orthogroups": [
            {
                "id": "og-all",
                "members": [
                    {
                        "recordKey": "record-1",
                        "biologicalFeatureId": f"feature-{index}",
                    }
                    for index in range(feature_count)
                ],
            }
        ],
        "annotations": [],
        "comparisonMatches": [],
        "sequenceSources": [
            {
                "key": "linear:record:0",
                "origin": "linear-record",
                "recordIndex": 0,
                "sequence": source_sequence,
            }
        ],
    }
    source = "<svg>" + "".join(
        (
            f'<path id="rendered-{index}" '
            f'data-gbdraw-feature-id="rendered-{index}" />'
        )
        for index in range(feature_count)
    ) + "</svg>"
    item = _export_profile(item)
    validation_calls = 0
    original = feature_catalog_module._canonical_sequence_source

    def counted(source_payload: object) -> str | None:
        nonlocal validation_calls
        validation_calls += 1
        return original(source_payload)  # type: ignore[arg-type]

    monkeypatch.setattr(
        feature_catalog_module,
        "_canonical_sequence_source",
        counted,
    )
    _validate_catalog_feature_ids(EXAMPLES[0], item, source)
    assert validation_calls == 1

    validation_calls = 0
    _validate_interactive_orthogroup_payload(EXAMPLES[0], item)
    assert validation_calls == 1


def test_gallery_payload_requires_complete_resolvable_orthogroups() -> None:
    feature = {
        "svg_id": "fstable",
        "orthogroup_id": "og_1",
        "nucleotide_sequence": "ATGAAATAA",
        "amino_acid_sequence": "MK*",
    }
    group = {
        "id": "og_1",
        "member_count": 1,
        "members": [{"feature_svg_id": "fstable"}],
    }

    _validate_interactive_orthogroup_payload(
        EXAMPLES[0],
        {"features": [feature], "orthogroups": [group]},
    )

    with pytest.raises(ValueError, match="missing 1 orthogroup"):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            {"features": [feature], "orthogroups": []},
        )

    with pytest.raises(ValueError, match="unresolved member"):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            {
                "features": [feature],
                "orthogroups": [
                    {"id": "og_1", "members": [{"feature_svg_id": "wrong"}]}
                ],
            },
        )

    with pytest.raises(ValueError, match="unresolved member"):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            {
                "features": [],
                "biological_features": [
                    {
                        "svg_id": "fstable",
                        "stable_svg_id": "fstable",
                        "record_idx": 0,
                        "nucleotide_sequence": "ATGAAATAA",
                    }
                ],
                "orthogroups": [
                    {
                        "id": "og_1",
                        "members": [
                            {"feature_svg_id": "wrong", "record_index": 0}
                        ],
                    }
                ],
            },
        )


def test_gallery_payload_keeps_hidden_orthogroup_members_in_biological_catalog() -> None:
    visible = {
        "recordKey": "record-1",
        "biologicalFeatureId": "fvisible",
        "nucleotide_sequence": "ATGAAATAA",
        "amino_acid_sequence": "MK*",
    }
    hidden = {
        "recordKey": "record-2",
        "biologicalFeatureId": "fhidden",
        "nucleotide_sequence": "ATGCCCTAA",
        "amino_acid_sequence": "MP*",
    }
    group = {
        "id": "og_1",
        "member_count": 2,
        "members": [
            {
                "recordKey": "record-1",
                "biologicalFeatureId": "fvisible",
            },
            {
                "recordKey": "record-2",
                "biologicalFeatureId": "fhidden",
            },
        ],
    }

    _validate_interactive_orthogroup_payload(
        EXAMPLES[0],
        {
            "features": [
                {
                    "svgId": "fvisible_record_1",
                    "recordKey": "record-1",
                    "biologicalFeatureId": "fvisible",
                }
            ],
            "biologicalFeatures": [visible, hidden],
            "orthogroups": [group],
        },
    )

    broken_hidden = dict(hidden, nucleotide_sequence="")
    with pytest.raises(ValueError, match="fhidden has no nucleotide sequence"):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            {
                "features": [
                    {
                        "svgId": "fvisible_record_1",
                        "recordKey": "record-1",
                        "biologicalFeatureId": "fvisible",
                    }
                ],
                "biologicalFeatures": [visible, broken_hidden],
                "orthogroups": [group],
            },
        )


def test_gallery_session_features_seed_biological_catalog() -> None:
    feature = {
        "svgId": "fstable_record_1",
        "recordKey": "record-1",
        "biologicalFeatureId": "fstable",
    }
    biological = {
        "recordKey": "record-1",
        "biologicalFeatureId": "fstable",
        "stableFeatureId": "fstable",
        "anchorProfile": {
            "precision": "unavailable",
            "operator": "unknown",
            "partOrder": "ambiguous",
            "strand": "unstranded",
        },
    }
    hidden = {
        "recordKey": "record-2",
        "biologicalFeatureId": "fhidden",
        "stableFeatureId": "fhidden",
        "anchorProfile": {
            "precision": "unavailable",
            "operator": "unknown",
            "partOrder": "ambiguous",
            "strand": "unstranded",
        },
    }

    context = _session_interactive_context(
        {
            "results": [{"name": "result", "content": "<svg/>"}],
            "editorState": {
                "featureCatalog": {
                    "schema": 4,
                    "items": [
                        {
                            "resultIndex": 0,
                            "resultName": "result",
                            "recordKeys": ["record-1", "record-2"],
                            "features": [feature],
                            "biologicalFeatures": [biological, hidden],
                            "orthogroups": [],
                            "annotations": [],
                            "comparisonMatches": [],
                        }
                    ],
                }
            },
        }
    )

    # A schema-4 catalog reads as schema 5 without drawn selector values.
    assert list(context.features) == [{**feature, "drawnSelector": None}]
    assert list(context.biological_features) == [biological, hidden]


@pytest.mark.parametrize(
    (
        "example_id",
        "expected_groups",
        "expected_members",
        "expected_hidden_members",
    ),
    [
        # 576 groups and 2566 members: the Session as the current derivation stores it, after
        # the self-hit and normalization changes of #1497 (OV-2306).
        ("hepatoplasmataceae_orthogroup", 576, 2566, 0),
        ("majanivirus_orthogroup", 152, 826, 0),
    ],
)
def test_orthogroup_gallery_preserves_session_members_and_rendered_ids(
    example_id: str,
    expected_groups: int,
    expected_members: int,
    expected_hidden_members: int,
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
    load_cached_svg_root: Callable[[Path], ET.Element],
) -> None:
    example = next(item for item in EXAMPLES if item.id == example_id)
    session = session_view(load_cached_gallery_session(example.session_path))
    root = load_cached_svg_root(example.gallery_svg_path)
    metadata = next(
        element
        for element in root.iter()
        if element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    payload = json.loads(metadata.text or "{}")
    assert payload["schema"] == feature_catalog_module.INTERACTIVE_EMBED_SCHEMA
    assert len(payload["items"]) == 1
    item = payload["items"][0]
    session_catalog = session["editorState"]["featureCatalog"]
    assert session_catalog["schema"] == feature_catalog_module.FEATURE_CATALOG_SCHEMA
    session_item = session_catalog["items"][0]
    # The embed keeps the sequences (schema 5); without them it is the stored item.
    promoted = feature_catalog_module.promote_legacy_feature_catalog({**payload, "items": [item]})
    assert promoted["items"][0] == session_item

    features = item["features"]
    biological_features = item["biologicalFeatures"]
    svg_groups = item["orthogroups"]
    session_groups = session_item["orthogroups"]
    svg_members = [member for group in svg_groups for member in group["members"]]

    assert features
    assert biological_features
    assert len(svg_groups) == expected_groups == len(session_groups)
    assert len(svg_members) == expected_members

    def member_keys(
        groups: list[dict[str, object]],
    ) -> Counter[tuple[str, str, str]]:
        return Counter(
            (
                str(group["id"]),
                str(member["recordKey"]),
                str(member["biologicalFeatureId"]),
            )
            for group in groups
            for member in group["members"]  # type: ignore[index]
        )

    assert member_keys(svg_groups) == member_keys(session_groups)

    biological_by_key = {
        (
            str(feature["recordKey"]),
            str(feature["biologicalFeatureId"]),
        ): feature
        for feature in biological_features
    }
    assert len(biological_by_key) == len(biological_features)
    assert all(
        (
            str(member["recordKey"]),
            str(member["biologicalFeatureId"]),
        )
        in biological_by_key
        for member in svg_members
    )
    assert all(
        _catalog_nucleotide_sequence(
            biological_by_key[
                (
                    str(member["recordKey"]),
                    str(member["biologicalFeatureId"]),
                )
            ],
            item,
        )
        for member in svg_members
    )

    rendered_keys = {
        (
            str(feature["recordKey"]),
            str(feature["biologicalFeatureId"]),
        )
        for feature in features
    }
    hidden_members = [
        member
        for member in svg_members
        if (
            str(member["recordKey"]),
            str(member["biologicalFeatureId"]),
        )
        not in rendered_keys
    ]
    assert len(hidden_members) == expected_hidden_members
    dom_ids = {
        value
        for element in root.iter()
        for value in (element.get("id"), element.get("data-gbdraw-feature-id"))
        if value
    }
    assert all(
        str(feature["svgId"]) in dom_ids for feature in features
    )



def test_refresh_rewrites_retired_flags_in_the_recorded_cli_invocation(
    tmp_path: Path,
) -> None:
    invocation = {
        "schema": 1,
        "mode": "linear",
        "args": [
            "--protein_blastp_mode",
            "orthogroup",
            "--gbk",
            "A.gb",
            "--align_orthogroup_feature",
            "P1",
            "--losatp_threads",
            "2",
        ],
        "fileBindings": [{"argIndex": 3, "slot": "files.linearSeqs[0].gb", "name": "A.gb"}],
    }
    # A staged file the publication bridge wrote.
    path = _write_staged(tmp_path / "s.gbdraw", {"keep": [1, 2], "cliInvocation": invocation})

    refresh_gallery_sessions_module._canonicalize_recorded_cli_invocation(
        path, mode="linear", gallery_id="s"
    )

    session = refresh_gallery_sessions_module._read_session_file(path)
    assert session["keep"] == [1, 2]
    args = session["cliInvocation"]["args"]
    assert args == [
        "--losat",
        "losatp",
        "--losatp_mode",
        "similarity_groups",
        "--gbk",
        "A.gb",
        "--similarity_alignment_feature",
        "P1",
        "--losat_threads",
        "2",
    ]
    (binding,) = session["cliInvocation"]["fileBindings"]
    assert args[binding["argIndex"]] == "A.gb"


def _published_drawing(path: Path) -> dict[str, Any]:
    """The drawing a refreshed Session shows, as its one-drawing view: a
    Session 49 holds one drawing; a tracked Session 46-48 is its own view."""

    document = read_session_payload(path)
    return drawing_view(document) if "drawings" in document else document


def _stored_output_prefixes(session: dict[str, object]) -> tuple[object, object]:
    """The Output Prefix as the form field of the Result's mode slice and the
    render request store it."""

    mode = session["renderRequest"]["mode"]  # type: ignore[index]
    return (
        session["modes"][mode]["config"]["form"]["prefix"],  # type: ignore[index]
        session["renderRequest"]["output"]["prefix"],  # type: ignore[index]
    )


def _recorded_output(args: list[str]) -> str:
    (value,) = [
        args[index + 1] for index, arg in enumerate(args) if arg in {"-o", "--output"}
    ]
    return value


@pytest.mark.parametrize(
    "gallery_id",
    ["HmmtDNA_basic_circular", "HmmtDNA_ATskew", "lambda_basic_linear"],
)
def test_stored_session_refresh_publishes_gallery_id_and_input_file_names(
    gallery_id: str,
    tmp_path: Path,
) -> None:
    """The stored-Session route names the output and the inputs for the Gallery.

    A Generate after Load names its files after the stored Output Prefix, so the
    prefix is the Gallery id the card command uses, in the form field and the
    render request. The recorded invocation uses the id for ``-o`` and no
    maintainer-local path for an input.
    """

    source = _session_path(gallery_id)
    staged = tmp_path / source.name
    staged.write_bytes(source.read_bytes())
    destination = tmp_path / "published" / source.name
    destination.parent.mkdir()

    _refresh_one_session(staged, destination_path=destination)

    session = _published_drawing(destination)
    assert _stored_output_prefixes(session) == (gallery_id, gallery_id)
    args = session["cliInvocation"]["args"]
    assert _recorded_output(args) == gallery_id
    assert not [arg for arg in args if arg.startswith(("/", "~"))], args
    for binding in session["cliInvocation"]["fileBindings"]:
        assert args[binding["argIndex"]] == binding["name"]


def test_refresh_sets_every_stored_output_prefix() -> None:
    session = {
        "modes": {
            "circular": {"config": {"form": {"prefix": "out", "species": "kept"}}},
            "linear": {"config": {"form": {"prefix": "out"}}},
        },
        "renderRequest": {"mode": "circular", "output": {"prefix": "out", "format": "svg"}},
    }

    assert refresh_gallery_sessions_module._set_output_prefix(session, "card-id") == 3
    assert _stored_output_prefixes(session) == ("card-id", "card-id")
    assert session["modes"]["linear"]["config"]["form"]["prefix"] == "card-id"
    assert session["modes"]["circular"]["config"]["form"]["species"] == "kept"
    assert refresh_gallery_sessions_module._set_output_prefix(session, "card-id") == 0


def test_refresh_stores_one_record_sources_as_exactly_one() -> None:
    # The Web writes `all` for a Linear GenBank file without selector; the CLI
    # replay materializes a one-record file as `exactly_one`, and `finalize`
    # compares the two (OV-376). A file with more records keeps `all`.
    one = Path("tests/test_inputs/NC_001416.gb").read_bytes()

    def record(index: int) -> dict[str, object]:
        return {"recordKey": f"record-{index}", "cardinality": "all", "selector": None, "windows": [],
                "source": {"kind": "genbank", "resourceId": f"record-{index}-genbank"}}

    session = {
        "renderRequest": {"mode": "linear", "records": [record(1), record(2), {**record(3), "selector": {"kind": "recordIndex", "index": 0}}]},
        "resources": {
            f"record-{index}-genbank": {"encoding": "base64", "data": base64.b64encode(payload).decode()}
            for index, payload in ((1, one), (2, one + one), (3, one))
        },
    }

    assert refresh_gallery_sessions_module._materialize_one_record_sources(session) == 1
    assert [entry["cardinality"] for entry in session["renderRequest"]["records"]] == ["exactly_one", "all", "all"]
    assert refresh_gallery_sessions_module._materialize_one_record_sources(session) == 0


def test_refresh_stores_the_only_record_of_a_circular_source_without_selector() -> None:
    # The Web names the only record of a Circular GenBank file by its ID; the CLI
    # replay, which selects nothing, writes no selector, and `finalize` compares
    # the two. A selector of a file with more records, or of another ID, stays.
    one = Path("tests/test_inputs/NC_001416.gb").read_bytes()

    def record(index: int, value: str) -> dict[str, object]:
        return {"recordKey": f"record-{index}", "cardinality": "exactly_one",
                "selector": {"kind": "recordId", "value": value}, "windows": [],
                "source": {"kind": "genbank", "resourceId": f"record-{index}-genbank"}}

    session = {
        "renderRequest": {"mode": "circular", "records": [record(1, "NC_001416.1"), record(2, "NC_001416.1"), record(3, "other")]},
        "resources": {
            f"record-{index}-genbank": {"encoding": "base64", "data": base64.b64encode(payload).decode()}
            for index, payload in ((1, one), (2, one + one), (3, one))
        },
    }

    assert refresh_gallery_sessions_module._materialize_one_record_sources(session) == 1
    assert [entry["selector"] for entry in session["renderRequest"]["records"]] == [
        None, {"kind": "recordId", "value": "NC_001416.1"}, {"kind": "recordId", "value": "other"}
    ]
    assert refresh_gallery_sessions_module._materialize_one_record_sources(session) == 0


def _write_staged(path: Path, payload: dict[str, Any]) -> Path:
    """A staged Session file of the refresh holding ``payload`` and no resources."""

    refresh_gallery_sessions_module._write_staged_session(path, {**payload, "resources": {}})
    return path


def _write_invocation(tmp_path: Path, args: list[str], bindings: list[dict[str, object]]) -> Path:
    return _write_staged(
        tmp_path / "session.gbdraw",
        {"cliInvocation": {"schema": 1, "mode": "circular", "args": args, "fileBindings": bindings}},
    )


def test_recorded_invocation_uses_gallery_id_and_input_file_names(
    tmp_path: Path,
) -> None:
    root = Path(refresh_gallery_sessions_module.REPO_ROOT)
    inside = root / "tests" / "test_inputs" / "HmmtDNA.gbk"
    # Recorded in another checkout of the repository (a maintainer's machine).
    foreign = "/mnt/c/Users/someone/GitHub/gbdraw/tests/test_inputs/NC_001416.gb"
    path = _write_invocation(
        tmp_path,
        ["--gbk", str(inside), foreign, "-o", "/tmp/work/old", "--labels", "out"],
        [
            {"argIndex": 1, "slot": "files.c_gb", "name": "HmmtDNA.gbk"},
            {"argIndex": 2, "slot": "files.c_gb2", "name": "NC_001416.gb"},
        ],
    )

    refresh_gallery_sessions_module._canonicalize_recorded_cli_invocation(
        path, mode="circular", gallery_id="card-id"
    )

    assert refresh_gallery_sessions_module._read_session_file(path)["cliInvocation"]["args"] == [
        "--gbk",
        "HmmtDNA.gbk",
        "NC_001416.gb",
        "-o",
        "card-id",
        "--labels",
        "out",
    ]


def test_recorded_invocation_keeps_relative_inputs(tmp_path: Path) -> None:
    path = _write_invocation(
        tmp_path,
        ["--output=old", "--gbk", "HmmtDNA.gbk"],
        [{"argIndex": 2, "slot": "files.c_gb", "name": "HmmtDNA.gbk"}],
    )

    refresh_gallery_sessions_module._canonicalize_recorded_cli_invocation(
        path, mode="circular", gallery_id="card-id"
    )

    assert refresh_gallery_sessions_module._read_session_file(path)["cliInvocation"]["args"] == [
        "--output=card-id",
        "--gbk",
        "HmmtDNA.gbk",
    ]


@pytest.mark.parametrize(
    "outside",
    ["/tmp/elsewhere/HmmtDNA.gbk", "/home/someone/data/gbdraw/tests/missing.gbk", "C:\\data\\HmmtDNA.gbk"],
)
def test_recorded_invocation_rejects_an_input_outside_the_repository(
    outside: str,
    tmp_path: Path,
) -> None:
    path = _write_invocation(
        tmp_path,
        ["--gbk", outside, "-o", "x"],
        [{"argIndex": 1, "slot": "files.c_gb", "name": "HmmtDNA.gbk"}],
    )

    with pytest.raises(ValueError, match="outside the repository"):
        refresh_gallery_sessions_module._canonicalize_recorded_cli_invocation(
            path, mode="circular", gallery_id="card-id"
        )


def test_declared_command_invocation_keeps_its_declared_inputs(tmp_path: Path) -> None:
    path = _write_invocation(
        tmp_path,
        ["--gbk", "/elsewhere/in.gbk", "-o", "old"],
        [{"argIndex": 1, "slot": "files.c_gb", "name": "in.gbk"}],
    )

    refresh_gallery_sessions_module._canonicalize_recorded_cli_invocation(
        path, mode="circular", gallery_id="card-id", declared_command=True
    )

    assert refresh_gallery_sessions_module._read_session_file(path)["cliInvocation"]["args"] == [
        "--gbk",
        "/elsewhere/in.gbk",
        "-o",
        "card-id",
    ]


def test_feature_sources_of_a_web_saved_session_come_from_its_record_resources() -> None:
    # A Web-saved Session has no cliInvocation; its input names live in
    # webFiles.resourceOriginalNames, keyed by each record's source resource.
    session = {
        "renderRequest": {
            "records": [
                {"recordKey": "record-1", "source": {"kind": "genbank", "resourceId": "record-1-genbank"}},
                {"recordKey": "record-2", "source": {"kind": "genbank", "resourceId": "record-2-genbank"}},
            ]
        },
        "webFiles": {
            "resourceOriginalNames": {
                "colors-default-colors-file": "colors-default-colors.tsv",
                "record-2-genbank": "second.gbk",
                "record-1-genbank": "first.gbk",
            }
        },
    }

    assert gallery_assets_module._session_feature_sources(session) == [
        "first.gbk",
        "second.gbk",
    ]


def test_gallery_session_refresh_keeps_the_session_created_at(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    tracked = tmp_path / "kept.gbdraw"
    original = {"format": "gbdraw-session", "createdAt": "2026-01-02T03:04:05+00:00", "title": "old"}
    tracked.write_bytes(gzip.compress(json.dumps(original).encode("utf-8"), mtime=0))

    monkeypatch.setattr(refresh_gallery_sessions_module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(refresh_gallery_sessions_module, "SESSION_ROOT", tmp_path)
    monkeypatch.setattr(refresh_gallery_sessions_module, "_session_path", lambda name: tmp_path / name)
    monkeypatch.setattr(refresh_gallery_sessions_module, "load_session", lambda path: {"path": str(path)})
    monkeypatch.setattr(
        refresh_gallery_sessions_module, "_validate_staged_gallery_session", lambda *_args, **_kwargs: None
    )
    refresh_count = 0

    def fake_refresh(session_path: Path, *, destination_path: Path | None = None) -> None:
        # Each replay stamps the time it ran, as the CLI Session writer does.
        nonlocal refresh_count
        refresh_count += 1
        assert destination_path is not None
        refreshed = {**original, "createdAt": f"2026-10-19T00:00:0{refresh_count}+00:00", "title": "new"}
        destination_path.write_bytes(
            gzip.compress(json.dumps(refreshed).encode("utf-8"), mtime=refresh_count)
        )

    monkeypatch.setattr(refresh_gallery_sessions_module, "_refresh_one_session", fake_refresh)

    refresh_gallery_sessions_module.refresh_gallery_sessions((tracked.name,))
    first_run = tracked.read_bytes()
    refresh_gallery_sessions_module.refresh_gallery_sessions((tracked.name,))

    refreshed = refresh_gallery_sessions_module._read_session_file(tracked)
    assert refreshed["createdAt"] == original["createdAt"]
    assert refreshed["title"] == "new"
    assert tracked.read_bytes() == first_run


def test_staged_gallery_validator_requires_one_drawn_drawing(tmp_path: Path) -> None:
    """The staged Session is checked against the one-drawn-drawing rule the Web reads from the generator."""

    from gbdraw.session_drawings import _drawings_document

    view = {"createdAt": "2026-01-02T00:00:00Z", **_staged_geometry_session(mode="linear", records=None)}
    two_drawn = _drawings_document(
        [view, view], ids=["linear", "linear-2"], names=["Linear", "Linear 2"], active=0, project_views=[view, view]
    )

    with pytest.raises(ValueError, match="shows 1 drawn drawing; this one has 2"):
        _validate_staged_gallery_session_49(tmp_path / "two-drawn.gbdraw", two_drawn)


def test_staged_gallery_validator_requires_the_linear_comparison_plan_through_validate_session(
    tmp_path: Path,
) -> None:
    """A staged Linear draft without ``linearComparisonPlan`` is refused by the Session validator itself."""

    session = _staged_geometry_session(mode="linear", records=None)
    session["modes"] = {"linear": {"config": {"linearRecordLayout": {}}}}

    from gbdraw.exceptions import ValidationError

    with pytest.raises(ValidationError, match="Current Web comparison draft requires config.linearComparisonPlan"):
        _validate_staged_gallery_session(tmp_path / "no-plan.gbdraw", session)


def test_publication_bridge_reads_and_writes_session_files(
    tmp_path: Path,
) -> None:
    """The Node bridge reads a Session file and writes one; the tool checks it
    when it writes the Session file it publishes."""

    import os

    source = _session_path("HmmtDNA_basic_circular")
    output = tmp_path / "prepared.gbdraw"
    refresh_gallery_sessions_module._run_publication_bridge("prepare", source, output, env=dict(os.environ))

    assert output.read_bytes().startswith(b"PK\x03\x04")
    session_file = tmp_path / "published.gbdraw"
    refresh_gallery_sessions_module._write_session_file(
        session_file, refresh_gallery_sessions_module._read_session_file(output)
    )
    assert session_file.read_bytes().startswith(b"PK\x03\x04")
    assert load_session(session_file)["resources"] == load_session(source)["resources"]


def test_every_gallery_session_writer_writes_one_number_text(tmp_path: Path) -> None:
    """OV-3370: the bridge's write and the tools' Python edits give a Session one JSON text.

    The bridge wrote the text of `JSON.stringify` (`0.00001`) and a later edit through Python
    (`1e-05`), so a Session's bytes depended on which step wrote it last. Both write the Session
    file through `_write_staged_session`.
    """

    import os

    from tests.utils.session_files import read_session_text

    output = tmp_path / "prepared.gbdraw"
    refresh_gallery_sessions_module._run_publication_bridge(
        "prepare",
        _session_path("HmmtDNA_basic_circular"),
        output,
        env=dict(os.environ),
    )
    text = read_session_text(output)

    # The Session's E-value cutoff, 0.00001, in Python's text.
    assert '"evalue":1e-05' in text and '"evalue":0.00001' not in text
    edited = tmp_path / "edited.gbdraw"
    refresh_gallery_sessions_module._write_staged_session(
        edited, refresh_gallery_sessions_module._read_session_file(output)
    )
    assert read_session_text(edited) == text


def _losat_result(query: str, subject: str, runtime: dict[str, object], text: str) -> dict[str, object]:
    from gbdraw.comparisons.losat_jobs import losat_fasta_text, sha256_text

    return {
        "kind": "raw-losat", "identityKind": "nucleotide", "program": "blastn", "flow": "circular-conservation",
        "outfmt": "6", "args": ["--task", "megablast"], "key": f"{query[:4]}-{subject[:4]}", "text": text,
        "queryCanonicalHash": sha256_text(losat_fasta_text("q1", query)),
        "subjectCanonicalHash": sha256_text(losat_fasta_text("s1", subject)),
        "runtime": runtime,
    }


def test_refresh_searches_again_a_losat_result_of_another_wasm_build(tmp_path: Path) -> None:
    """A stored result Generate would search again is searched with the shipped LOSAT and keeps its rows."""

    import os
    import random

    from tools.refresh_gallery_sessions import _search_stale_losat_results

    bases = random.Random(1893)
    subject = "".join(bases.choice("ACGT") for _ in range(2000))
    query = subject[500:1100]
    old_build = {"kind": "losat", "source": "wasm", "version": "0.2.0+000000000000", "program": "blastn"}
    cli = {"kind": "losat", "source": "explicit", "version": "0.2.0", "program": "blastn"}
    session = {
        "resources": {
            "query": {"name": "q.fna", "data": f">q1\n{query}\n"},
            "subject": {"name": "s.fna", "data": base64.b64encode(f">s1 x\n{subject}\n".encode()).decode(), "encoding": "base64"},
        },
        "losatCache": {"entries": [_losat_result(query, subject, old_build, ""), _losat_result("A" * 40, "C" * 40, cli, "kept")]},
    }
    env = dict(os.environ)
    # The old build's rows differ from the shipped LOSAT's, and the Session draws results made from them.
    (tmp_path / "differ").mkdir()
    with pytest.raises(ValueError, match="gives other rows"):
        _search_stale_losat_results(copy.deepcopy(session), tmp_path / "differ", env=env)

    (tmp_path / "jobs").mkdir()
    jobs = tmp_path / "jobs" / "jobs.json"
    (tmp_path / "jobs" / "q.fa").write_text(f">q1\n{query}\n")
    (tmp_path / "jobs" / "s.fa").write_text(f">s1\n{subject}\n")
    jobs.write_text(json.dumps([{"program": "blastn", "outfmt": "6", "args": ["--task", "megablast"],
                                 "query": str(tmp_path / "jobs" / "q.fa"), "subject": str(tmp_path / "jobs" / "s.fa")}]))
    refresh_gallery_sessions_module._run_publication_bridge("search-losat", jobs, tmp_path / "jobs" / "out.json", env=env)
    [shipped] = json.loads((tmp_path / "jobs" / "out.json").read_text())
    assert shipped["text"].startswith("q1\ts1\t100.000\t600\t")

    session["losatCache"]["entries"][0]["text"] = shipped["text"]
    (tmp_path / "same").mkdir()
    assert _search_stale_losat_results(session, tmp_path / "same", env=env) == 1
    stale, reused = session["losatCache"]["entries"]
    assert stale["runtime"] == shipped["runtime"] != old_build
    assert stale["text"] == shipped["text"]
    assert (reused["runtime"], reused["text"]) == (cli, "kept")

    # A searched sequence that no resource holds whole cannot be searched again.
    stale["runtime"] = old_build
    stale["queryCanonicalHash"] = "0" * 64
    (tmp_path / "missing").mkdir()
    with pytest.raises(ValueError, match="no whole file or record"):
        _search_stale_losat_results(session, tmp_path / "missing", env=env)


def test_refresh_gives_every_gallery_file_the_writers_metadata() -> None:
    """OV-3331: a Web-saved Gallery Session kept each upload's browser metadata (a media type
    guessed from the name, the file's modification time); the refresh stores the writer's."""

    session = {
        "resources": {
            "record-1-genbank": {"kind": "genbank", "type": "application/x-gameboy-rom", "lastModified": 1792769561077},
            "depth": {"kind": "depth-track-file", "type": "text/tab-separated-values", "lastModified": 1794131353048},
            "table": {"kind": "canonical-tsv", "type": "text/tab-separated-values", "lastModified": 0},
        },
        "webFiles": {"bindings": {"linearSeqs": [{"gb": {
            "resourceId": "record-1-genbank", "name": "AP027078.gb",
            "type": "application/x-gameboy-rom", "lastModified": 1792769561077,
        }}]}},
    }

    changed = refresh_gallery_sessions_module._store_written_file_metadata(session)

    assert changed == 3
    resources = session["resources"]
    assert {key: (value["type"], value["lastModified"]) for key, value in resources.items()} == {
        "record-1-genbank": ("application/octet-stream", 0),
        "depth": ("application/octet-stream", 0),
        "table": ("text/tab-separated-values", 0),
    }
    binding = session["webFiles"]["bindings"]["linearSeqs"][0]["gb"]
    assert (binding["type"], binding["lastModified"], binding["name"]) == ("application/octet-stream", 0, "AP027078.gb")
    assert refresh_gallery_sessions_module._store_written_file_metadata(session) == 0
