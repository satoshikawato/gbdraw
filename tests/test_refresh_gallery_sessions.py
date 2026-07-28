from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    save_session_document,
)
from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    load_session,
)
import tools.prepare_interactive_gallery_assets as gallery_assets_module
import tools.refresh_gallery_sessions as refresh_gallery_sessions_module
from tools.prepare_interactive_gallery_assets import (
    EXAMPLES,
    _migrate_legacy_multipart_feature_ids,
    _session_interactive_context,
    _session_result_svg,
    _validate_interactive_orthogroup_payload,
    _validate_source_feature_ids,
)
from tools.refresh_gallery_sessions import (
    VIBRIO_EXPECTED_RAW_PAIRS,
    VIBRIO_RAW_ENTRY_COUNT,
    _gallery_file_transaction,
    _merge_refreshed_gallery_artifacts,
    _preserve_gallery_cli_invocation,
    _refresh_one_session,
    _session_path,
    _validate_gallery_session_inventory,
    _validate_staged_gallery_session,
    _with_interactive_svg_format,
)


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


def test_vibrio_gallery_session_retains_complete_compact_cache() -> None:
    path = _session_path("vibrio-harveyi-group-collinear")
    session = load_session(path)

    _validate_staged_gallery_session(path, session, artifact_path=path)

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


def test_all_bundled_current_sessions_use_current_artifact_schemas() -> None:
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
        session = load_session(path)
        assert session["version"] == CURRENT_SESSION_VERSION, path
        assert session["renderRequest"]["schema"] == 3, path
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


def test_bundled_gallery_sources_match_current_session_results() -> None:
    for example in EXAMPLES:
        if not example.sync_result_svg:
            continue
        session = load_session(example.session_path)
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
        "version": 36,
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
    assert merged["version"] == 36
    assert merged["createdAt"] == "fresh"
    assert merged["losatCache"] == refreshed["losatCache"]
    assert merged["losatDerivedCache"] == refreshed["losatDerivedCache"]
    assert merged["proteinIdentityManifest"] == refreshed["proteinIdentityManifest"]
    assert "legacyArtifacts" not in merged


def _staged_geometry_session(
    *,
    mode: str,
    records: list[dict[str, object]] | None,
    tracks: dict[str, object] | None = None,
) -> dict[str, object]:
    request: dict[str, object] = {"schema": 3, "mode": mode}
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
        ),
    )

    _refresh_one_session(source, destination_path=destination)

    refreshed = load_session(destination)
    geometry = refreshed["runMetadata"]["trackSlotGeometry"]
    assert geometry["schema"] == 1
    assert geometry["mode"] == "circular"
    assert geometry["source"] == "resolved"
    assert geometry["records"][0]["axisRadiusPx"] > 0


def test_staged_gallery_validator_requires_current_artifact_schemas(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "synthetic.gbdraw-session.json"
    session = {
        "format": "gbdraw-session",
        "version": 36,
        "renderRequest": {"schema": 3, "mode": "linear"},
        "resources": {},
        "results": [{"name": "result", "content": "<svg></svg>"}],
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
    with pytest.raises(ValueError, match="expected 36"):
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
        "svg_id": "fvisible_record_1",
        "stable_svg_id": "fvisible",
        "record_idx": 0,
        "orthogroup_id": "og_1",
        "nucleotide_sequence": "ATGAAATAA",
        "amino_acid_sequence": "MK*",
    }
    hidden = {
        "svg_id": "fhidden",
        "stable_svg_id": "fhidden",
        "record_idx": 1,
        "orthogroup_id": "og_1",
        "nucleotide_sequence": "ATGCCCTAA",
        "amino_acid_sequence": "MP*",
    }
    group = {
        "id": "og_1",
        "member_count": 2,
        "members": [
            {
                "feature_svg_id": "fvisible",
                "stable_feature_svg_id": "fvisible",
                "rendered_feature_svg_id": "fvisible_record_1",
                "record_index": 0,
            },
            {
                "feature_svg_id": "fhidden",
                "stable_feature_svg_id": "fhidden",
                "record_index": 1,
            },
        ],
    }

    _validate_interactive_orthogroup_payload(
        EXAMPLES[0],
        {
            "features": [visible],
            "biological_features": [visible, hidden],
            "orthogroups": [group],
        },
    )

    broken_hidden = dict(hidden, nucleotide_sequence="")
    with pytest.raises(ValueError, match="fhidden has no nucleotide sequence"):
        _validate_interactive_orthogroup_payload(
            EXAMPLES[0],
            {
                "features": [visible],
                "biological_features": [visible, broken_hidden],
                "orthogroups": [group],
            },
        )


def test_gallery_session_features_seed_biological_catalog() -> None:
    feature = {"svg_id": "fstable_record_1", "stable_svg_id": "fstable"}
    hidden = {"svg_id": "fhidden", "stable_svg_id": "fhidden"}

    context = _session_interactive_context(
        {
            "features": {
                "extractedFeatures": [feature],
                "biologicalFeatures": [feature, hidden],
            }
        }
    )

    assert list(context.features) == [feature]
    assert list(context.biological_features) == [feature, hidden]


@pytest.mark.parametrize(
    (
        "example_id",
        "expected_rendered_features",
        "expected_biological_features",
        "expected_groups",
        "expected_members",
        "expected_hidden_members",
    ),
    [
        ("hepatoplasmataceae_orthogroup", 2994, 5987, 577, 2566, 0),
        ("majanivirus_orthogroup", 999, 1008, 152, 826, 0),
    ],
)
def test_orthogroup_gallery_preserves_session_members_and_rendered_ids(
    example_id: str,
    expected_rendered_features: int,
    expected_biological_features: int,
    expected_groups: int,
    expected_members: int,
    expected_hidden_members: int,
) -> None:
    example = next(item for item in EXAMPLES if item.id == example_id)
    session = load_session(example.session_path)
    root = ET.parse(example.gallery_svg_path).getroot()
    metadata = next(
        element
        for element in root.iter()
        if element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    payload = json.loads(metadata.text or "{}")

    features = payload["features"]
    biological_features = payload["biological_features"]
    svg_groups = payload["orthogroups"]
    session_groups = session["orthogroupState"]["groups"]
    svg_members = [member for group in svg_groups for member in group["members"]]

    assert len(features) == expected_rendered_features
    assert len(biological_features) == expected_biological_features
    assert len(svg_groups) == expected_groups == len(session_groups)
    assert len(svg_members) == expected_members

    def record_index(item: dict[str, object]) -> int:
        value = item.get("record_index", item.get("recordIndex", item.get("record_idx", -1)))
        return int(value)

    def stable_id(item: dict[str, object]) -> str:
        return str(
            item.get("stable_feature_svg_id")
            or item.get("stableFeatureSvgId")
            or item.get("stable_feature_id")
            or item.get("stable_svg_id")
            or item.get("feature_svg_id")
            or item.get("featureSvgId")
            or item.get("svg_id")
            or ""
        )

    def rendered_id(item: dict[str, object]) -> str:
        return str(
            item.get("rendered_feature_svg_id")
            or item.get("renderedFeatureSvgId")
            or ""
        )

    def member_keys(groups: list[dict[str, object]]) -> Counter[tuple[str, int, str]]:
        return Counter(
            (
                str(group["id"]),
                record_index(member),
                stable_id(member),
            )
            for group in groups
            for member in group["members"]  # type: ignore[index]
        )

    assert member_keys(svg_groups) == member_keys(session_groups)

    biological_by_key = {
        (record_index(feature), stable_id(feature)): feature
        for feature in biological_features
    }
    assert len(biological_by_key) == len(biological_features)
    assert all(
        (record_index(member), stable_id(member)) in biological_by_key
        for member in svg_members
    )
    assert all(
        biological_by_key[(record_index(member), stable_id(member))].get(
            "nucleotide_sequence"
        )
        for member in svg_members
    )

    hidden_members = [member for member in svg_members if not rendered_id(member)]
    assert len(hidden_members) == expected_hidden_members
    dom_ids = {
        value
        for element in root.iter()
        for value in (element.get("id"), element.get("data-gbdraw-feature-id"))
        if value
    }
    assert all(
        not rendered_id(member) or rendered_id(member) in dom_ids
        for member in svg_members
    )
