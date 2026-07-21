from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.session_io import load_session
import tools.refresh_gallery_sessions as refresh_gallery_sessions_module
from tools.prepare_interactive_gallery_assets import (
    EXAMPLES,
    _migrate_legacy_multipart_feature_ids,
    _session_interactive_context,
    _validate_interactive_orthogroup_payload,
    _validate_source_feature_ids,
)
from tools.refresh_gallery_sessions import (
    _gallery_file_transaction,
    _merge_refreshed_gallery_artifacts,
    _preserve_gallery_cli_invocation,
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


def test_session_path_prefers_existing_compressed_gallery_session() -> None:
    path = _session_path("vibrio-harveyi-group-collinear")

    assert path.name == "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"


def test_gallery_session_inventory_matches_files_and_examples() -> None:
    _validate_gallery_session_inventory()


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
        "version": 34,
        "renderRequest": {"diagramOptions": {"palette": "curated"}},
        "config": {"labels": "curated"},
        "resources": {
            "curated": {"data": "promoted"},
            "promoted-only": {"data": "keep"},
        },
        "results": [{"content": "stale"}],
        "features": {"extractedFeatures": []},
        "losatCache": {"entries": [{"schema": 2, "program": "blastp"}]},
        "legacyArtifacts": {"proteinRawCandidates": {"schema": 1, "entries": ["old"]}},
    }
    refreshed = {
        "version": 35,
        "renderRequest": {"diagramOptions": {"palette": "default"}},
        "config": {"labels": "lost"},
        "resources": {
            "curated": {"data": "refreshed-conflict"},
            "render-only": {"data": "keep-too"},
        },
        "results": [{"content": "fresh"}],
        "features": {"extractedFeatures": [{"svg_id": "feature-1"}]},
        "orthogroupState": {"groups": [{"id": "og_1"}]},
        "losatCache": {"entries": [{"schema": 3, "program": "blastp"}]},
        "losatDerivedCache": {"entries": [{"schema": 2}]},
        "proteinIdentityManifest": {
            "schema": 1,
            "proteinSets": {},
            "recordAnalyses": {},
            "recordInstances": {},
        },
    }

    merged = _merge_refreshed_gallery_artifacts(promoted, refreshed)

    assert merged["renderRequest"] == promoted["renderRequest"]
    assert merged["config"] == promoted["config"]
    assert merged["resources"]["curated"] == {"data": "promoted"}
    assert merged["resources"]["promoted-only"] == {"data": "keep"}
    assert merged["resources"]["render-only"] == {"data": "keep-too"}
    assert merged["results"] == refreshed["results"]
    assert merged["features"] == refreshed["features"]
    assert merged["orthogroupState"] == refreshed["orthogroupState"]
    assert merged["version"] == 35
    assert merged["losatCache"] == refreshed["losatCache"]
    assert merged["losatDerivedCache"] == refreshed["losatDerivedCache"]
    assert merged["proteinIdentityManifest"] == refreshed["proteinIdentityManifest"]
    assert "legacyArtifacts" not in merged


def test_staged_gallery_validator_requires_current_artifact_schemas(
    tmp_path: Path,
) -> None:
    session_path = tmp_path / "synthetic.gbdraw-session.json"
    session = {
        "format": "gbdraw-session",
        "version": 35,
        "renderRequest": {"schema": 3, "mode": "linear"},
        "resources": {},
        "results": [{"name": "result", "content": "<svg></svg>"}],
        "losatCache": {"entries": []},
        "losatDerivedCache": {"entries": []},
        "proteinIdentityManifest": {
            "schema": 1,
            "proteinSets": {},
            "recordAnalyses": {},
            "recordInstances": {},
        },
    }

    _validate_staged_gallery_session(session_path, session)

    stale_version = dict(session, version=34)
    with pytest.raises(ValueError, match="expected 35"):
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
        lambda *_args: None,
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
    assert "if example.source_svg_path.exists():" in source
    assert "example.source_svg_path.write_text(migrated" in source
    assert "_validate_source_feature_ids(example, session, migrated)" in source
    assert "def _sync_session_result_svg(" in source
    assert "write_session_json(example.session_path, session)" in source
    assert "_sync_session_result_svg(example, session, source)" in source
    assert "_write_gallery_svg(example, session, source)" in source
    assert 'entry["tutorial"] = f"./tutorials/{example.id}.json"' in source
    assert 'entry["tutorialStatus"] = "ready"' in source


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
