from __future__ import annotations

import base64
from collections import Counter
from collections.abc import Callable
import gzip
import json
from pathlib import Path
import subprocess
import sys
from xml.etree import ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.web_support.feature_catalog as feature_catalog_module
from gbdraw.api import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
    save_session_document,
)
from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.exceptions import ValidationError
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    load_session,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
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
    _validate_source_feature_ids,
)
from tools.refresh_gallery_sessions import (
    GALLERY_SESSION_FILES,
    TEST_INPUT_SESSION_FILES,
    VIBRIO_EXPANDED_HARD_LIMIT,
    VIBRIO_EXPANDED_REGRESSION_CEILING,
    VIBRIO_EXPECTED_RAW_PAIRS,
    VIBRIO_GZIP_HARD_LIMIT,
    VIBRIO_GZIP_REGRESSION_CEILING,
    VIBRIO_RAW_ENTRY_COUNT,
    _enable_gallery_interactive_metadata,
    _drop_unreferenced_duplicate_resources,
    _gallery_file_transaction,
    _load_gallery_refresh_source,
    _merge_refreshed_gallery_artifacts,
    _omit_regenerable_gallery_derived_cache,
    _preserve_gallery_cli_invocation,
    _refresh_one_session,
    _refresh_session_paths,
    _restore_rendered_palette_file_binding,
    _session_artifact_measurements,
    _session_path,
    _sync_circular_track_draft_with_render_request,
    _sync_legacy_legend_control_with_render_request,
    _validate_current_session_catalog_structure,
    _validate_gallery_session_inventory,
    _validate_staged_gallery_session,
    _with_interactive_svg_format,
)


pytestmark = pytest.mark.gallery


def test_gallery_refresh_syncs_legacy_legend_control_with_render_request() -> None:
    session = {
        "renderRequest": {
            "diagramOptions": {"output": {"legend": "right"}},
        },
        "config": {"form": {"legend": "bottom", "prefix": "example"}},
    }

    assert _sync_legacy_legend_control_with_render_request(session) is True
    assert session["config"]["form"] == {
        "legend": "right",
        "prefix": "example",
    }
    assert _sync_legacy_legend_control_with_render_request(session) is False


def test_gallery_refresh_syncs_circular_track_draft_with_render_request() -> None:
    slots = [
        {
            "kind": "circularTrackSlot",
            "id": "regions",
            "renderer": "annotations",
            "enabled": True,
            "side": "inside",
            "radius": {"value": 0.65, "unit": "factor"},
            "width": {"value": 20, "unit": "px"},
            "z": 0,
            "params": {"set_id": "plastome_regions"},
            "innerGapPx": 1,
            "outerGapPx": 1,
        }
    ]
    session = {
        "renderRequest": {
            "mode": "circular",
            "diagramOptions": {
                "tracks": {
                    "circularTrackSlots": slots,
                    "circularTrackAxisIndex": 0,
                }
            },
        },
        "config": {"adv": {"circular_track_slots_enabled": True}},
    }

    assert _sync_circular_track_draft_with_render_request(session) is True
    assert session["config"]["adv"] == {
        "circular_track_slots_enabled": True,
        "circular_track_slots_schema_version": 4,
        "circular_track_slots_axis_index": 0,
        "circular_track_slots": slots,
    }
    assert session["config"]["adv"]["circular_track_slots"] is not slots
    assert _sync_circular_track_draft_with_render_request(session) is False


def _color_resource(kind: str, text: str) -> dict[str, object]:
    data = text.encode("utf-8")
    return {
        "kind": kind,
        "name": "colors.tsv",
        "type": "text/tab-separated-values",
        "size": len(data),
        "lastModified": 0,
        "encoding": "base64",
        "data": base64.b64encode(data).decode("ascii"),
    }


def test_default_refresh_inventory_covers_gallery_and_test_input_sessions() -> None:
    paths = _refresh_session_paths(None)

    assert len(paths) == 13
    assert len(GALLERY_SESSION_FILES) == 11
    assert len(TEST_INPUT_SESSION_FILES) == 2
    assert {path.name for path in paths} >= set(TEST_INPUT_SESSION_FILES)


def test_palette_binding_restores_the_file_proven_by_rendered_legend() -> None:
    colors = {
        "defaultColors": {
            "resourceId": "colors-default-colors",
            "representation": "canonicalTsv",
        },
        "defaultColorsPalette": "ajisai",
        "defaultColorsFile": None,
    }
    session = {
        "renderRequest": {"diagramOptions": {"colors": colors}},
        "resources": {
            "colors-default-colors": _color_resource(
                "canonical-tsv", "feature_type\tcolor\nCDS\t#54bcf8\n"
            ),
            "colors-default-colors-file": _color_resource(
                "colors-default-colors-file", "CDS\t#84b9ec\nrRNA\t#7cecd5\n"
            ),
        },
        "config": {
            "palette": "ajisai",
            "colors": {"CDS": "#84b9ec", "rRNA": "#7cecd5"},
        },
        "editorState": {
            "legend": {
                "originalColors": {"CDS": "#84b9ec", "rRNA": "#7cecd5"}
            }
        },
    }

    assert _restore_rendered_palette_file_binding(session) is True
    assert colors["defaultColors"] is None
    assert colors["defaultColorsFile"] == {
        "resourceId": "colors-default-colors-file",
        "representation": "file",
    }


def test_palette_binding_ignores_a_draft_palette_not_proven_by_legend() -> None:
    colors = {
        "defaultColors": {
            "resourceId": "colors-default-colors",
            "representation": "canonicalTsv",
        },
        "defaultColorsPalette": "orange",
        "defaultColorsFile": None,
    }
    session = {
        "renderRequest": {"diagramOptions": {"colors": colors}},
        "resources": {
            "colors-default-colors": _color_resource(
                "canonical-tsv", "feature_type\tcolor\nCDS\t#54bcf8\n"
            ),
            "colors-default-colors-file": _color_resource(
                "colors-default-colors-file", "CDS\t#dddddd\n"
            ),
        },
        "config": {"palette": "orange", "colors": {"CDS": "#dddddd"}},
        "editorState": {
            "legend": {"originalColors": {"other proteins": "#dddddd"}}
        },
    }

    assert _restore_rendered_palette_file_binding(session) is False
    assert colors["defaultColors"]["resourceId"] == "colors-default-colors"
    assert colors["defaultColorsFile"] is None


def test_duplicate_resource_cleanup_keeps_the_referenced_copy_only() -> None:
    payload = {
        "kind": "canonical-tsv",
        "name": "fresh.tsv",
        "type": "text/tab-separated-values",
        "size": 3,
        "lastModified": 0,
        "encoding": "base64",
        "data": "QUJD",
    }
    stale_payload = {**payload, "name": "stale.tsv"}
    session = {
        "renderRequest": {
            "comparisons": [{"resourceId": "fresh"}],
        },
        "resources": {
            "fresh": payload,
            "stale": stale_payload,
        },
        "webFiles": {
            "resourceOriginalNames": {
                "fresh": "fresh.tsv",
                "stale": "stale.tsv",
            }
        },
    }

    _drop_unreferenced_duplicate_resources(session)

    assert list(session["resources"]) == ["fresh"]
    assert session["webFiles"]["resourceOriginalNames"] == {
        "fresh": "fresh.tsv"
    }


def test_duplicate_resource_cleanup_ignores_stale_resource_metadata() -> None:
    session = {
        "webFiles": {
            "conservationLosatFastaSources": ["canonical"]
        },
        "resources": {
            "canonical": {
                "kind": "conservation-fasta-file",
                "name": "canonical.fasta",
                "type": "application/octet-stream",
                "size": 3,
                "lastModified": 0,
                "encoding": "base64",
                "data": "QUJD",
            },
            "stale": {
                "kind": "web-file",
                "name": "stale.fasta",
                "type": "text/plain",
                "size": 999,
                "lastModified": 17,
                "encoding": "base64",
                "data": "QUJD",
            },
        },
    }

    _drop_unreferenced_duplicate_resources(session)

    assert list(session["resources"]) == ["canonical"]


def test_duplicate_resource_cleanup_preserves_two_referenced_copies() -> None:
    payload = {
        "kind": "canonical-tsv",
        "type": "text/tab-separated-values",
        "size": 3,
        "lastModified": 0,
        "encoding": "base64",
        "data": "QUJD",
    }
    session = {
        "renderRequest": {
            "comparisons": [
                {"resourceId": "first"},
                {"resourceId": "second"},
            ],
        },
        "resources": {
            "first": {**payload, "name": "first.tsv"},
            "second": {**payload, "name": "second.tsv"},
        },
    }

    _drop_unreferenced_duplicate_resources(session)

    assert list(session["resources"]) == ["first", "second"]


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


def test_with_interactive_svg_format_replaces_existing_format() -> None:
    assert _with_interactive_svg_format(["-o", "out", "-f", "svg"]) == [
        "-o",
        "out",
        "-f",
        "interactive_svg",
    ]
    assert _with_interactive_svg_format(["--format=svg"]) == [
        "--format=interactive_svg"
    ]
    assert _with_interactive_svg_format(["--gbk", "input.gb"]) == [
        "--gbk",
        "input.gb",
        "-f",
        "interactive_svg",
    ]


def test_gallery_refresh_enables_metadata_for_forced_interactive_output() -> None:
    session = {
        "renderRequest": {
            "output": {
                "formats": ["interactive_svg"],
                "interactiveMetadataPolicy": "omit",
            }
        }
    }

    assert _enable_gallery_interactive_metadata(session) is True
    assert session["renderRequest"]["output"]["interactiveMetadataPolicy"] == "auto"
    assert _enable_gallery_interactive_metadata(session) is False


def test_gallery_refresh_discards_only_an_invalid_derived_result_catalog(
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGCGCAT"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    source = tmp_path / "stale.gbdraw-session.json"
    save_session_document(
        source,
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        ),
    )
    payload = load_session(source)
    payload["results"] = [{"name": "out", "content": "<svg/>"}]
    payload["editorState"]["featureCatalog"] = {"schema": 3, "items": []}
    payload["runMetadata"] = {"stale": True}
    source.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValidationError, match="one schema-3 item per Result"):
        load_session(source)

    recovered = _load_gallery_refresh_source(source)
    assert recovered["results"] == []
    assert recovered["editorState"]["featureCatalog"] == {
        "schema": 3,
        "items": [],
    }
    assert "runMetadata" not in recovered
    assert recovered["renderRequest"] == payload["renderRequest"]
    assert recovered["resources"] == payload["resources"]


def test_gallery_session_promoter_runs_mjs_without_obsolete_esm_flag(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    commands: list[list[str]] = []
    source = tmp_path / "source.json"
    output = tmp_path / "output.json"

    monkeypatch.setattr(
        refresh_gallery_sessions_module.shutil,
        "which",
        lambda executable: "/opt/node" if executable == "node" else None,
    )
    monkeypatch.setattr(
        refresh_gallery_sessions_module.subprocess,
        "run",
        lambda command, **kwargs: commands.append(command),
    )
    monkeypatch.setattr(refresh_gallery_sessions_module, "load_session", lambda path: {})

    refresh_gallery_sessions_module._promote_gallery_session(
        source,
        output,
        env={},
    )

    assert commands == [
        [
            "/opt/node",
            str(refresh_gallery_sessions_module.SESSION_PROMOTER),
            str(source),
            str(output),
        ]
    ]


def test_session_path_prefers_existing_compressed_gallery_session() -> None:
    path = _session_path("vibrio-harveyi-group-collinear")

    assert path.name == "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"


def test_gallery_artifact_limits_match_phase_two_targets() -> None:
    assert INTERACTIVE_SVG_HARD_LIMIT == 41_943_040
    assert DECODED_METADATA_HARD_LIMIT == 200_000_000
    assert DECODED_METADATA_HARD_LIMIT <= DECODED_METADATA_REGRESSION_CEILING
    assert DECODED_METADATA_REGRESSION_CEILING == 220_000_000
    assert VIBRIO_GZIP_HARD_LIMIT == 90_000_000
    assert VIBRIO_GZIP_HARD_LIMIT <= VIBRIO_GZIP_REGRESSION_CEILING
    assert VIBRIO_GZIP_REGRESSION_CEILING == 95_000_000
    assert VIBRIO_EXPANDED_HARD_LIMIT == 400_000_000
    assert (
        VIBRIO_EXPANDED_HARD_LIMIT
        <= VIBRIO_EXPANDED_REGRESSION_CEILING
    )
    assert VIBRIO_EXPANDED_REGRESSION_CEILING == 420_000_000


def test_interactive_svg_measurements_decode_compressed_schema_three() -> None:
    payload = {
        "schema": 3,
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
        'data-schema="3" data-encoding="gzip-base64">'
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
    assert measurements["renderedFeatureCount"] == 24_945
    assert measurements["biologicalFeatureCount"] == 49_970


def test_session_artifact_measurements_report_component_bytes(
    tmp_path: Path,
) -> None:
    session = {
        "results": [{"name": "result", "content": "<svg/>"}],
        "resources": {"record": {"data": "QUJD", "encoding": "base64"}},
        "webFiles": {"files": []},
        "editorState": {
            "featureCatalog": {"schema": 3, "items": []},
            "legend": {"entries": []},
        },
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
    }
    expanded = json.dumps(
        session,
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    artifact_path = tmp_path / "session.json.gz"
    artifact_path.write_bytes(gzip.compress(expanded, mtime=0))

    measurements = _session_artifact_measurements(session, artifact_path)

    assert measurements["sessionGzipBytes"] == artifact_path.stat().st_size
    assert measurements["sessionExpandedBytes"] == len(expanded)
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


def test_current_session_catalog_structure_rejects_duplicate_payloads(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "duplicate.gbdraw-session.json"
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
            "featureCatalog": {"schema": 3, "items": [item]}
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


def test_vibrio_gallery_session_retains_complete_compact_cache(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    path = _session_path("vibrio-harveyi-group-collinear")
    session = load_cached_gallery_session(path)

    _validate_staged_gallery_session(
        path,
        session,
        artifact_path=path,
    )

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


def test_gallery_session_inventory_matches_files_and_examples() -> None:
    _validate_gallery_session_inventory()


def test_all_bundled_sessions_use_current_request_and_artifact_schemas(
    load_cached_gallery_session: Callable[[Path], dict[str, object]],
) -> None:
    repo_root = Path(__file__).parents[1]
    paths = sorted(
        (repo_root / "gbdraw" / "web" / "gallery" / "sessions").glob(
            "*.gbdraw-session.json*"
        )
    )
    paths.extend(
        sorted(
            (repo_root / "tests" / "test_inputs").glob("*.gbdraw-session.json*")
        )
    )

    assert len(paths) == 13
    for path in paths:
        session = load_cached_gallery_session(path)
        assert session["version"] == CURRENT_SESSION_VERSION, path
        assert session["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA, path
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
        session = load_cached_gallery_session(example.session_path)
        assert (
            example.source_svg_path.read_text(encoding="utf-8")
            == _session_result_svg(session, example)
        ), example.id


def test_preserve_gallery_cli_invocation_keeps_original_render_args() -> None:
    source_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": [
                "--definition_line_style",
                "name:font_weight=bold",
                "--circular_track_slot",
                "a_skew_2:dinucleotide_skew@w=0.1,nt=AT,legend_label=AT skew",
                "--conservation_blast",
                "comparison.tsv",
                "--gbk",
                "input.gb",
                "-f",
                "svg",
            ],
            "renderFormats": ["svg"],
            "fileBindings": [
                {
                    "argIndex": 5,
                    "slot": "files.c_conservation_blasts[0]",
                    "name": "comparison.tsv",
                },
                {"argIndex": 7, "slot": "files.c_gb", "name": "input.gb"},
            ],
            "generatedBy": "gbdraw",
        },
    }
    refreshed_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": ["--gbk", "input.gb", "-f", "interactive_svg"],
            "renderFormats": ["interactive_svg"],
            "fileBindings": [
                {"argIndex": 1, "slot": "files.c_gb", "name": "input.gb"}
            ],
            "generatedBy": "gbdraw",
        },
    }

    preserved = _preserve_gallery_cli_invocation(
        source_session,
        refreshed_session,
        mode="circular",
    )

    assert preserved is True
    args = refreshed_session["cliInvocation"]["args"]
    assert "--definition_line_style" in args
    assert "--circular_track_slot" in args
    assert any("AT skew" in arg for arg in args)
    assert "--conservation_blast" in args
    assert args[-2:] == ["-f", "interactive_svg"]
    assert refreshed_session["cliInvocation"]["renderFormats"] == ["interactive_svg"]
    assert (
        refreshed_session["cliInvocation"]["fileBindings"]
        == source_session["cliInvocation"]["fileBindings"]
    )


def test_preserve_gallery_cli_invocation_reports_missing_source_cli() -> None:
    refreshed_session = {
        "cliInvocation": {
            "schema": 1,
            "mode": "circular",
            "args": ["--gbk", "input.gb", "-f", "interactive_svg"],
        },
    }

    preserved = _preserve_gallery_cli_invocation(
        {"config": {"form": {"prefix": "old"}}},
        refreshed_session,
        mode="circular",
    )

    assert preserved is False
    assert refreshed_session["cliInvocation"]["args"] == [
        "--gbk",
        "input.gb",
        "-f",
        "interactive_svg",
    ]


def test_refreshed_gallery_artifacts_do_not_replace_promoted_render_authority() -> None:
    promoted = {
        "format": "gbdraw-session",
        "renderRequest": {"diagramOptions": {"palette": "curated"}},
        "config": {"labels": "curated"},
        "resources": {
            "curated": {"data": "promoted"},
            "promoted-only": {"data": "keep"},
        },
        "results": [{"content": "stale"}],
        "features": {"extractedFeatures": []},
        "runMetadata": {"trackSlotGeometry": {"records": ["stale"]}},
        "losatCache": {"entries": [{"schema": 2, "program": "blastp"}]},
        "legacyArtifacts": {"proteinRawCandidates": {"schema": 1, "entries": ["old"]}},
        "version": 33,
        "createdAt": "old",
    }
    refreshed = {
        "version": CURRENT_SESSION_VERSION,
        "createdAt": "fresh",
        "renderRequest": {"diagramOptions": {"palette": "default"}},
        "config": {"labels": "lost"},
        "resources": {
            "curated": {"data": "refreshed-conflict"},
            "render-only": {"data": "keep-too"},
        },
        "results": [{"content": "fresh"}],
        "features": {"extractedFeatures": [{"svg_id": "feature-1"}]},
        "orthogroupState": {"groups": [{"id": "og_1"}]},
        "runMetadata": {
            "trackSlotGeometry": {
                "schema": 1,
                "mode": "linear",
                "source": "resolved",
                "records": [{"recordIndex": 0}],
            }
        },
        "losatCache": {"entries": [{"schema": 4, "program": "blastp"}]},
        "losatDerivedCache": {"entries": [{"schema": 3}]},
        "proteinIdentityManifest": {
            "schema": 2,
            "proteinSets": {},
            "recordAnalyses": {},
            "recordInstances": {},
        },
    }

    merged = _merge_refreshed_gallery_artifacts(promoted, refreshed)

    assert list(merged)[:5] == [
        "format",
        "version",
        "createdAt",
        "renderRequest",
        "resources",
    ]
    assert merged["renderRequest"] == promoted["renderRequest"]
    assert merged["config"] == promoted["config"]
    assert merged["resources"]["curated"] == {"data": "promoted"}
    assert merged["resources"]["promoted-only"] == {"data": "keep"}
    assert merged["resources"]["render-only"] == {"data": "keep-too"}
    assert merged["results"] == refreshed["results"]
    assert merged["features"] == refreshed["features"]
    assert merged["orthogroupState"] == refreshed["orthogroupState"]
    assert merged["runMetadata"] == refreshed["runMetadata"]
    assert merged["version"] == CURRENT_SESSION_VERSION
    assert merged["createdAt"] == "fresh"
    assert merged["losatCache"] == refreshed["losatCache"]
    assert merged["losatDerivedCache"] == refreshed["losatDerivedCache"]
    assert merged["proteinIdentityManifest"] == refreshed["proteinIdentityManifest"]
    assert "legacyArtifacts" not in merged


def test_refreshed_request_keeps_its_resource_on_identifier_collision() -> None:
    promoted = {
        "format": "gbdraw-session",
        "version": CURRENT_SESSION_VERSION,
        "renderRequest": {
            "schema": CANONICAL_REQUEST_SCHEMA,
            "diagramOptions": {"colors": "stale"},
        },
        "resources": {"colors-default-colors": {"data": "stale"}},
    }
    refreshed = {
        "version": CURRENT_SESSION_VERSION,
        "createdAt": "fresh",
        "renderRequest": {
            "schema": CANONICAL_REQUEST_SCHEMA,
            "diagramOptions": {"colors": "fresh"},
        },
        "resources": {"colors-default-colors": {"data": "fresh"}},
    }

    merged = _merge_refreshed_gallery_artifacts(promoted, refreshed)

    assert merged["renderRequest"] == refreshed["renderRequest"]
    assert merged["resources"]["colors-default-colors"] == {"data": "fresh"}


def test_vibrio_gallery_refresh_omits_regenerable_derived_cache(
    tmp_path: Path,
) -> None:
    vibrio_session = {
        "losatDerivedCache": {"entries": [{"schema": LOSAT_DERIVED_CACHE_SCHEMA}]}
    }
    _omit_regenerable_gallery_derived_cache(
        tmp_path / "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
        vibrio_session,
    )
    assert vibrio_session["losatDerivedCache"] == {"entries": []}

    other_session = {
        "losatDerivedCache": {"entries": [{"schema": LOSAT_DERIVED_CACHE_SCHEMA}]}
    }
    _omit_regenerable_gallery_derived_cache(
        tmp_path / "other.gbdraw-session.json.gz",
        other_session,
    )
    assert other_session["losatDerivedCache"]["entries"]


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
                "schema": 3,
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
    session_path = tmp_path / "missing-geometry.gbdraw-session.json"
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
    session_path = tmp_path / "circular-geometry.gbdraw-session.json"
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
    session_path = tmp_path / "linear-geometry.gbdraw-session.json"
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
    session_path = tmp_path / "linear-featureless-geometry.gbdraw-session.json"
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
    session_path = tmp_path / "linear-custom-geometry.gbdraw-session.json"
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


def test_refresh_records_resolved_track_geometry(
    tmp_path: Path,
) -> None:
    record = SeqRecord(
        Seq("ATGCGCAT"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    source = tmp_path / "source.gbdraw-session.json"
    destination = tmp_path / "refreshed.gbdraw-session.json"
    save_session_document(
        source,
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            output=RenderOutputRequest(
                formats=("interactive_svg",),
                interactive_metadata_policy="omit",
            ),
        ),
    )

    _refresh_one_session(source, destination_path=destination)

    refreshed = load_session(destination)
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


def test_staged_gallery_validator_accepts_current_artifact_schemas(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "synthetic.gbdraw-session.json"
    session = {
        "format": "gbdraw-session",
        "version": CURRENT_SESSION_VERSION,
        "renderRequest": {"schema": CANONICAL_REQUEST_SCHEMA, "mode": "linear"},
        "resources": {},
        "results": [{"name": "result", "content": "<svg></svg>"}],
        "editorState": {
            "featureCatalog": {
                "schema": 3,
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
    first = tmp_path / "first.gbdraw-session.json"
    second = tmp_path / "second.gbdraw-session.json"
    first.write_text("first-original", encoding="utf-8")
    second.write_text("second-original", encoding="utf-8")

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
        destination_path.write_text("first-refreshed", encoding="utf-8")

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

    assert first.read_text(encoding="utf-8") == "first-original"
    assert second.read_text(encoding="utf-8") == "second-original"


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
    assert "write_session_json(example.session_path, session)" in source
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
        gallery_assets_module,
        "write_session_json",
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
    session = {
        "results": [{"name": "result", "content": source}],
        "editorState": {
            "featureCatalog": {"schema": 3, "items": [item]}
        },
    }

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
    _validate_source_feature_ids(EXAMPLES[0], session, source)
    _validate_interactive_orthogroup_payload(EXAMPLES[0], item)

    simple_item = json.loads(json.dumps(item))
    simple_feature = simple_item["biologicalFeatures"][0]
    simple_feature.pop("location_parts")
    simple_feature["strand"] = "+"
    simple_session = {
        "results": [{"name": "result", "content": source}],
        "editorState": {
            "featureCatalog": {"schema": 3, "items": [simple_item]}
        },
    }
    _validate_source_feature_ids(EXAMPLES[0], simple_session, source)

    broken_item = json.loads(json.dumps(item))
    broken_item["biologicalFeatures"][0]["location_parts"][1]["end"] = 16
    broken_session = {
        "results": [{"name": "result", "content": source}],
        "editorState": {
            "featureCatalog": {"schema": 3, "items": [broken_item]}
        },
    }
    with pytest.raises(
        ValueError,
        match="invalid nucleotide sequence source reference",
    ):
        _validate_source_feature_ids(EXAMPLES[0], broken_session, source)
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
    unreferenced_invalid_session = {
        "results": [{"name": "result", "content": source}],
        "editorState": {
            "featureCatalog": {
                "schema": 3,
                "items": [unreferenced_invalid_item],
            }
        },
    }
    with pytest.raises(
        ValueError,
        match="invalid nucleotide sequence source reference",
    ):
        _validate_source_feature_ids(
            EXAMPLES[0],
            unreferenced_invalid_session,
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
    session = {
        "results": [{"name": "result", "content": source}],
        "editorState": {
            "featureCatalog": {"schema": 3, "items": [item]}
        },
    }
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
    _validate_source_feature_ids(EXAMPLES[0], session, source)
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
    }
    hidden = {
        "recordKey": "record-2",
        "biologicalFeatureId": "fhidden",
        "stableFeatureId": "fhidden",
    }

    context = _session_interactive_context(
        {
            "results": [{"name": "result", "content": "<svg/>"}],
            "editorState": {
                "featureCatalog": {
                    "schema": 3,
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

    assert list(context.features) == [feature]
    assert list(context.biological_features) == [biological, hidden]


@pytest.mark.parametrize(
    (
        "example_id",
        "expected_groups",
        "expected_members",
        "expected_hidden_members",
    ),
    [
        ("hepatoplasmataceae_orthogroup", 577, 2566, 0),
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
    session = load_cached_gallery_session(example.session_path)
    root = load_cached_svg_root(example.gallery_svg_path)
    metadata = next(
        element
        for element in root.iter()
        if element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    payload = json.loads(metadata.text or "{}")
    assert payload["schema"] == 3
    assert len(payload["items"]) == 1
    item = payload["items"][0]
    session_item = session["editorState"]["featureCatalog"]["items"][0]
    assert item == session_item

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
