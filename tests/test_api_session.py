from __future__ import annotations

import base64
import copy
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
)
from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.api.session_compat import canonical_projection_for_session_decode
from gbdraw.api import (
    CircularDiagramRequest,
    CircularDiagramOptions,
    ColorOptions,
    DepthTrackInput,
    InMemoryRecordSource,
    LinearDiagramRequest,
    LinearDiagramOptions,
    LinearTrackOptions,
    LinearTrackSlot,
    RecordInput,
    RenderOutputRequest,
    ScalarSpec,
    SessionFormatError,
    SessionConversionError,
    SessionRenderError,
    SessionResourceError,
    SessionVersionError,
    build_session_document,
    load_session_document,
    materialize_session,
    render_session,
    save_session_document,
    session_to_request,
    normalize_request_records,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION


_SESSION_FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sessions"
_V40_STYLE_FIXTURE = (
    _SESSION_FIXTURE_DIR / "feature-style-v40-gallery-minimal.gbdraw-session.json"
)
_V40_STYLE_AMBIGUOUS_MUTATION = (
    _SESSION_FIXTURE_DIR / "feature-style-v40-gallery-minimal.ambiguous-mutation.json"
)
_V40_STYLE_EXPECTED = (
    _SESSION_FIXTURE_DIR / "feature-style-v40-gallery-minimal.expected.json"
)
_V40_GALLERY_SESSION = (
    Path(__file__).parent.parent
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "BGC0000708-BGC0000713.gbdraw-session.json"
)
_V40_TOBACCO_SESSION = (
    Path(__file__).parent.parent
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "tobacco-chloroplast.gbdraw-session.json"
)


def _record(record_id: str = "record") -> RecordInput:
    seqrecord = SeqRecord(
        Seq("ATGCGCAT"),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )
    return RecordInput(source=InMemoryRecordSource(seqrecord))


def _color_resource(text: str, *, kind: str, name: str) -> dict[str, object]:
    content = text.encode("utf-8")
    return {
        "kind": kind,
        "name": name,
        "type": "text/tab-separated-values",
        "size": len(content),
        "lastModified": 0,
        "encoding": "base64",
        "data": base64.b64encode(content).decode("ascii"),
    }


def _v40_color_session() -> dict[str, object]:
    rules = pd.DataFrame(
        [
            ("CDS", "gene_kind", "biosynthetic$", "#d03535", "Core genes"),
            ("CDS", "gene_kind", "transport", "#577edb", "Transport genes"),
        ],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    stale_defaults = pd.DataFrame(
        [("CDS", "#54bcf8"), ("tRNA", "#e8b441"), ("default", "#d3d3d3")],
        columns=["feature_type", "color"],
    )
    data = build_session_document(
        CircularDiagramRequest(
            records=(_record(),),
            options=CircularDiagramOptions(
                colors=ColorOptions(
                    color_table=rules,
                    default_colors=stale_defaults,
                    default_colors_palette="orange",
                )
            ),
        )
    ).to_dict()
    data["version"] = 40
    data["renderRequest"]["schema"] = 5
    editor_defaults = {
        "CDS": "#dddddd",
        "tRNA": "#9cd9cf",
        "default": "#d3d3d3",
    }
    data["config"] = {
        "palette": "orange",
        "colors": copy.deepcopy(editor_defaults),
        "rules": [
            {
                "feat": "CDS",
                "qual": "gene_kind",
                "val": "biosynthetic$",
                "color": "#d03535",
                "cap": "Core genes",
            },
            {
                "feat": "CDS",
                "qual": "gene_kind",
                "val": "transport",
                "color": "#577edb",
                "cap": "Transport genes",
            },
        ],
    }
    data["ui"] = {
        "appliedPaletteName": "orange",
        "appliedPaletteColors": copy.deepcopy(editor_defaults),
    }
    data["resources"]["colors-default-colors-file"] = _color_resource(
        "CDS\t#dddddd\ntRNA\t#9cd9cf\ndefault\t#d3d3d3\n",
        kind="colors-default-colors-file",
        name="saved-default-colors.tsv",
    )
    data["resources"]["colors-color-table-file"] = _color_resource(
        "CDS\tgene_kind\tbiosynthetic$\t#d03535\tCore genes\n"
        "CDS\tgene_kind\ttransport\t#577edb\tTransport genes\n",
        kind="colors-color-table-file",
        name="saved-specific-colors.tsv",
    )
    data["results"] = [
        {
            "name": "out",
            "content": (
                '<svg><path fill="#d03535" id="core-feature"/>'
                '<path fill="#577edb" id="transport-feature"/>'
                '<path fill="#54bcf8" id="default-feature"/></svg>'
            ),
        }
    ]
    data["features"] = {
        "featureColorOverrides": {
            "file0_f1": {"color": "#d03535", "caption": "Core genes"},
            "file0_f2": {
                "color": "#577edb",
                "caption": "Transport genes",
            },
        }
    }
    biological_features = [
        {
            "recordKey": "record-1",
            "biologicalFeatureId": "core-bio",
            "stableFeatureId": "core-feature",
            "sourceFeatureIndex": 1,
            "type": "CDS",
            "qualifiers": {"gene_kind": ["biosynthetic"]},
        },
        {
            "recordKey": "record-1",
            "biologicalFeatureId": "transport-bio",
            "stableFeatureId": "transport-feature",
            "sourceFeatureIndex": 2,
            "type": "CDS",
            "qualifiers": {"gene_kind": ["transport"]},
        },
        {
            "recordKey": "record-1",
            "biologicalFeatureId": "default-bio",
            "stableFeatureId": "default-feature",
            "sourceFeatureIndex": 3,
            "type": "CDS",
            "qualifiers": {},
        },
    ]
    data["editorState"] = {
        "legend": {
            "entries": [
                {"caption": "Core genes", "color": "#d03535"},
                {"caption": "Transport genes", "color": "#577edb"},
                {"caption": "other proteins", "color": "#dddddd"},
            ]
        },
        "featureCatalog": {
            "schema": 3,
            "items": [
                {
                    "resultIndex": 0,
                    "resultName": "out",
                    "recordKeys": ["record-1"],
                    "biologicalFeatures": biological_features,
                    "features": [
                        {
                            "recordKey": "record-1",
                            "biologicalFeatureId": feature["biologicalFeatureId"],
                            "svgId": feature["stableFeatureId"],
                            "fillColor": color,
                        }
                        for feature, color in zip(
                            biological_features,
                            ("#d03535", "#577edb", "#54bcf8"),
                            strict=True,
                        )
                    ],
                }
            ],
        },
    }
    return data


@pytest.mark.parametrize(
    ("request_type", "mode"),
    ((CircularDiagramRequest, "circular"), (LinearDiagramRequest, "linear")),
)
def test_session_document_round_trip_owns_resource_lifetime(
    tmp_path: Path,
    request_type,
    mode: str,
) -> None:
    request = request_type(
        records=(_record(),),
        output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
    )
    document = build_session_document(request, title="round-trip")

    assert document.version == CURRENT_SESSION_VERSION
    assert document.mode == mode
    assert document.has_canonical_request is True
    assert document.to_dict()["resources"]["record-1-genbank"]["encoding"] == "base64"
    assert document.to_dict()["renderRequest"]["output"]["overwrite"] is False
    assert document.to_dict()["results"] == []
    assert document.to_dict()["editorState"]["featureCatalog"] is None

    with materialize_session(
        document,
        output_directory=tmp_path,
        temporary_directory=tmp_path / "materialized",
    ) as materialized:
        decoded = session_to_request(materialized)
        assert decoded.output.overwrite is False
        resource_path = decoded.records[0].source.path
        assert resource_path.is_file()
        assert materialized.active is True
        rebuilt = build_session_document(decoded)
        assert [entry["name"] for entry in rebuilt.to_dict()["resources"].values()] == [
            entry["name"] for entry in document.to_dict()["resources"].values()
        ]

    assert materialized.active is False
    assert not resource_path.exists()
    with pytest.raises(SessionResourceError, match="no longer active"):
        session_to_request(materialized)


def test_session_document_save_load_and_render(tmp_path: Path) -> None:
    request = CircularDiagramRequest(
        records=(_record(),),
        output=RenderOutputRequest(
            output_prefix="canonical",
            output_directory=tmp_path,
            overwrite=True,
        ),
    )
    session_path = tmp_path / "canonical.gbdraw-session.json"
    save_session_document(session_path, request)
    assert session_path.read_text(encoding="utf-8").startswith(
        f'{{"format":"gbdraw-session","version":{CURRENT_SESSION_VERSION},'
    )
    document = load_session_document(session_path)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        result = render_session(materialized)

    assert result.mode == "circular"
    assert result.output_paths == (tmp_path / "canonical.svg",)
    assert result.output_paths[0].is_file()


def test_session_document_gzip_save_load(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw-session.json.gz"

    save_session_document(session_path, request)
    document = load_session_document(session_path)

    assert session_path.read_bytes().startswith(b"\x1f\x8b")
    assert document.version == CURRENT_SESSION_VERSION
    assert document.mode == "circular"


def test_session_document_requires_explicit_overwrite(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=(_record(),))
    session_path = tmp_path / "canonical.gbdraw-session.json"

    save_session_document(session_path, request)
    original = session_path.read_bytes()

    with pytest.raises(SessionFormatError, match="overwrite=True"):
        save_session_document(session_path, request)
    assert session_path.read_bytes() == original

    save_session_document(session_path, request, overwrite=True)
    assert load_session_document(session_path).mode == "circular"


def test_direct_session_render_ignores_stored_overwrite_permission(
    tmp_path: Path,
) -> None:
    request = CircularDiagramRequest(
        records=(_record(),),
        output=RenderOutputRequest(
            output_prefix="protected",
            output_directory=tmp_path,
            overwrite=True,
        ),
    )
    document = build_session_document(request)
    document._data["renderRequest"]["output"]["overwrite"] = True
    protected = tmp_path / "protected.svg"
    protected.write_text("keep this diagram", encoding="utf-8")

    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionRenderError, match="already exist"):
            render_session(materialized)

    assert protected.read_text(encoding="utf-8") == "keep this diagram"


def test_current_document_quarantines_legacy_protein_cache_on_save(
    tmp_path: Path,
) -> None:
    legacy_entry = {
        "schema": 2,
        "kind": "raw-losat",
        "key": "legacy-key",
        "text": "p_r_old\tp_r_other\n",
        "program": "blastp",
    }
    session_path = tmp_path / "legacy-round-trip.gbdraw-session.json"

    saved = save_session_document(
        session_path,
        LinearDiagramRequest(records=(_record(),)),
        adjunct={"losatCache": {"entries": [legacy_entry]}},
    )
    reloaded = load_session_document(session_path)

    assert saved.version == reloaded.version == CURRENT_SESSION_VERSION
    assert reloaded.to_dict()["losatCache"]["entries"] == []
    assert reloaded.to_dict()["legacyArtifacts"]["proteinRawCandidates"][
        "entries"
    ] == [
        {
            "state": "pending",
            "originalEntry": legacy_entry,
            "rejectionReason": None,
        }
    ]


def test_legacy_session_has_no_public_typed_conversion(tmp_path: Path) -> None:
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": 30,
            "files": {},
            "ui": {"mode": "circular"},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionVersionError, match="internal CLI replay only"):
            session_to_request(materialized)


def test_duplicate_json_resource_id_is_rejected(tmp_path: Path) -> None:
    session_path = tmp_path / "duplicate.gbdraw-session.json"
    session_path.write_text(
        '{"format":"gbdraw-session","version":31,'
        '"renderRequest":{},"resources":{"same":{},"same":{}}}',
        encoding="utf-8",
    )

    with pytest.raises(SessionFormatError, match="duplicate object key"):
        load_session_document(session_path)


def test_duplicate_sanitized_resource_name_is_rejected(tmp_path: Path) -> None:
    data = build_session_document(
        LinearDiagramRequest(records=(_record("a"), _record("b")))
    ).to_dict()
    resource_ids = list(data["resources"])
    data["resources"][resource_ids[1]]["name"] = data["resources"][resource_ids[0]]["name"]

    with pytest.raises(SessionResourceError, match="Duplicate canonical resource filename"):
        load_session_document(data)


@pytest.mark.parametrize(
    ("session_version", "request_schema"),
    ((40, 6), (CURRENT_SESSION_VERSION, 5)),
)
def test_session_document_rejects_cross_era_request_schema_pairings(
    session_version: int,
    request_schema: int,
) -> None:
    data = build_session_document(
        CircularDiagramRequest(records=(_record(),))
    ).to_dict()
    data["version"] = session_version
    data["renderRequest"]["schema"] = request_schema

    with pytest.raises(
        SessionVersionError,
        match=(
            rf"Session version {session_version} is incompatible with canonical "
            rf"renderRequest schema {request_schema}"
        ),
    ):
        load_session_document(data)


def test_partial_materialization_failure_cleans_owned_directory(tmp_path: Path) -> None:
    data = build_session_document(
        LinearDiagramRequest(records=(_record("a"), _record("b")))
    ).to_dict()
    data["resources"]["record-2-genbank"]["data"] = "not-base64!"
    root = tmp_path / "materialized"

    with pytest.raises(SessionResourceError, match="record-2-genbank"):
        with materialize_session(
            data,
            output_directory=tmp_path,
            temporary_directory=root,
        ):
            pytest.fail("materialization unexpectedly succeeded")

    assert root.is_dir()
    assert list(root.iterdir()) == []


def test_session_document_returns_detached_payload() -> None:
    document = build_session_document(CircularDiagramRequest(records=(_record(),)))
    detached = document.to_dict()
    detached["renderRequest"]["mode"] = "linear"

    assert document.mode == "circular"
    json.dumps(document.to_dict())


def test_v40_gallery_shape_recovers_stale_default_color_resource(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)
        assert decoded.options.colors is not None
        defaults = decoded.options.colors.default_colors
        assert defaults is None
        default_colors_file = decoded.options.colors.default_colors_file
        assert default_colors_file is not None
        default_rows = pd.read_csv(
            default_colors_file,
            sep="\t",
            names=["feature_type", "color"],
            header=None,
            dtype=str,
        )
    assert default_rows.set_index("feature_type")["color"].to_dict() == {
        "CDS": "#dddddd",
        "default": "#d3d3d3",
        "tRNA": "#9cd9cf",
    }
    rules = decoded.options.colors.color_table
    assert rules is not None
    assert rules.to_dict("records") == [
        {
            "feature_type": "CDS",
            "qualifier_key": "gene_kind",
            "value": "biosynthetic$",
            "color": "#d03535",
            "caption": "Core genes",
        },
        {
            "feature_type": "CDS",
            "qualifier_key": "gene_kind",
            "value": "transport",
            "color": "#577edb",
            "caption": "Transport genes",
        },
    ]


def test_v40_consistent_specific_color_table_is_not_materialized_again(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    stale_defaults = data["resources"]["colors-default-colors"]
    assert isinstance(stale_defaults, dict)
    stale_defaults.update(
        _color_resource(
            "feature_type\tcolor\n"
            "CDS\t#dddddd\n"
            "tRNA\t#9cd9cf\n"
            "default\t#d3d3d3\n",
            kind="canonical-tsv",
            name="colors-default-colors.tsv",
        )
    )
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        before = set(materialized.temp_directory.iterdir())
        decoded = session_to_request(materialized)
        after = set(materialized.temp_directory.iterdir())

    assert before == after
    assert decoded.options.colors is not None
    assert decoded.options.colors.color_table is not None
    assert decoded.options.colors.color_table["caption"].tolist() == [
        "Core genes",
        "Transport genes",
    ]


def test_v40_specific_color_recovery_requires_agreeing_saved_editor_table(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    data["config"]["rules"][0]["color"] = "#abcdef"
    data["resources"]["colors-color-table-file"] = _color_resource(
        "CDS\tgene_kind\tbiosynthetic$\t#abcdef\tCore genes\n"
        "CDS\tgene_kind\ttransport\t#577edb\tTransport genes\n",
        kind="colors-color-table-file",
        name="saved-specific-colors.tsv",
    )
    data["features"]["featureColorOverrides"]["file0_f1"]["color"] = "#abcdef"
    data["editorState"]["legend"]["entries"][0]["color"] = "#abcdef"
    data["editorState"]["featureCatalog"]["items"][0]["features"][0][
        "fillColor"
    ] = "#abcdef"
    data["results"][0]["content"] = data["results"][0]["content"].replace(
        'fill="#d03535"',
        'fill="#abcdef"',
    )

    with materialize_session(data, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)
        assert decoded.options.colors is not None
        assert decoded.options.colors.color_table is None
        color_table_file = decoded.options.colors.color_table_file
        assert color_table_file is not None
        recovered_rules = pd.read_csv(
            color_table_file,
            sep="\t",
            names=["feature_type", "qualifier_key", "value", "color", "caption"],
            header=None,
            dtype=str,
        )
    assert recovered_rules["color"].tolist() == [
        "#abcdef",
        "#577edb",
    ]


def test_v40_color_recovery_rejects_ambiguous_editor_sources(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    data["ui"]["appliedPaletteColors"]["CDS"] = "#abcdef"

    with materialize_session(data, output_directory=tmp_path) as materialized:
        before = set(materialized.temp_directory.iterdir())
        with pytest.raises(
            SessionConversionError,
            match="conflicting applied color authorities",
        ):
            session_to_request(materialized)
        assert set(materialized.temp_directory.iterdir()) == before


def test_v40_color_recovery_does_not_mutate_session_or_resource_mapping(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    original = copy.deepcopy(data)
    document = load_session_document(data)
    document_before = document.to_dict()

    with materialize_session(document, output_directory=tmp_path) as materialized:
        resource_paths_before = dict(materialized.resource_paths)
        session_to_request(materialized)
        assert materialized.resource_paths == resource_paths_before

    assert data == original
    assert document.to_dict() == document_before


def test_v40_schema5_biological_instance_hash_remains_a_source_qualifier(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    data["config"]["rules"][0].update(
        {"qual": "instance_hash", "val": "legacy-[0-9]+"}
    )
    data["resources"]["colors-color-table"] = _color_resource(
        "feature_type\tqualifier_key\tvalue\tcolor\tcaption\n"
        "CDS\tinstance_hash\tlegacy-[0-9]+\t#d03535\tCore genes\n"
        "CDS\tgene_kind\ttransport\t#577edb\tTransport genes\n",
        kind="canonical-tsv",
        name="colors-color-table.tsv",
    )
    biological = data["editorState"]["featureCatalog"]["items"][0][
        "biologicalFeatures"
    ][0]
    biological["qualifiers"]["instance_hash"] = ["legacy-42"]

    with materialize_session(data, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.colors is not None
    assert decoded.options.colors.color_table is not None
    assert decoded.options.colors.color_table.iloc[0]["qualifier_key"] == "instance_hash"


def test_v40_schema5_reserved_instance_selector_requires_regeneration(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    data["config"]["rules"][0].update(
        {"qual": "__gbdraw_instance_hash__", "val": "fi1_invalid"}
    )

    with materialize_session(data, output_directory=tmp_path) as materialized:
        with pytest.raises(
            SessionConversionError,
            match="cannot safely promote a schema-5 reserved instance selector",
        ):
            session_to_request(materialized)


def test_v40_schema5_reserved_semantic_selector_requires_regeneration(
    tmp_path: Path,
) -> None:
    data = _v40_color_session()
    data["config"]["rules"][0].update(
        {"qual": "__gbdraw_semantic_scope__", "val": "^fs1:"}
    )

    with materialize_session(data, output_directory=tmp_path) as materialized:
        with pytest.raises(
            SessionConversionError,
            match="cannot safely promote a schema-5 reserved semantic selector",
        ):
            session_to_request(materialized)


def test_v40_request_without_explicit_color_options_remains_decodable(
    tmp_path: Path,
) -> None:
    data = build_session_document(
        CircularDiagramRequest(records=(_record(),))
    ).to_dict()
    data["version"] = 40
    data["renderRequest"]["schema"] = 5

    with materialize_session(data, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.colors is None


def test_v40_gallery_minimum_matches_shared_normalized_projection(
    tmp_path: Path,
) -> None:
    session = json.loads(_V40_STYLE_FIXTURE.read_text(encoding="utf-8"))
    expected = json.loads(_V40_STYLE_EXPECTED.read_text(encoding="utf-8"))
    projection = expected["normalizedProjection"]
    document = load_session_document(session)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        payload, resource_paths = canonical_projection_for_session_decode(
            document.version,
            document.to_dict(),
            resource_paths=materialized.resource_paths,
            temp_directory=materialized.temp_directory,
        )
        colors = payload["diagramOptions"]["colors"]
        assert payload["schema"] == projection["requestSchema"]
        assert colors["colorTable"]["resourceId"] == projection["specificResourceId"]
        assert colors["colorTableFile"] is None
        assert colors["defaultColors"] is None
        assert colors["defaultColorsFile"] == {
            "resourceId": projection["defaultResourceId"],
            "representation": projection["defaultResourceRepresentation"],
        }
        assert colors["defaultColorsPalette"] == projection["defaultColorsPalette"]
        assert resource_paths == materialized.resource_paths

        decoded = session_to_request(materialized)
        assert decoded.options.colors is not None
        assert decoded.options.colors.color_table is not None
        rules = [
            {
                "feat": row.feature_type,
                "qual": row.qualifier_key,
                "val": row.value,
                "color": row.color,
                "cap": row.caption,
            }
            for row in decoded.options.colors.color_table.itertuples(index=False)
        ]
        assert rules == projection["rules"]
        default_file = decoded.options.colors.default_colors_file
        assert default_file is not None
        default_rows = pd.read_csv(
            default_file,
            sep="\t",
            names=["feature_type", "color"],
            header=None,
            dtype=str,
        )
        assert default_rows.set_index("feature_type")["color"].to_dict() == projection[
            "appliedDefaultColors"
        ]

    item = session["editorState"]["featureCatalog"]["items"][0]
    biological_by_source_index = {
        feature["sourceFeatureIndex"]: feature
        for feature in item["biologicalFeatures"]
    }
    normalized_overrides = {}
    for alias, override in session["features"]["featureColorOverrides"].items():
        source_index = int(alias.rsplit("_f", maxsplit=1)[1])
        feature = biological_by_source_index[source_index]
        feature_key = f'{feature["recordKey"]}\0{feature["biologicalFeatureId"]}'
        normalized_overrides[feature_key] = override
    assert normalized_overrides == projection["derivedOverrides"]


def test_v40_gallery_minimum_ambiguous_mutation_fails_closed(
    tmp_path: Path,
) -> None:
    session = json.loads(_V40_STYLE_FIXTURE.read_text(encoding="utf-8"))
    mutation = json.loads(
        _V40_STYLE_AMBIGUOUS_MUTATION.read_text(encoding="utf-8")
    )
    duplicated = copy.deepcopy(session["resources"][mutation["sourceResourceId"]])
    duplicated["name"] = mutation["targetName"]
    session["resources"][mutation["targetResourceId"]] = duplicated

    with materialize_session(session, output_directory=tmp_path) as materialized:
        with pytest.raises(
            SessionConversionError,
            match=mutation["expectedError"],
        ):
            session_to_request(materialized)


@pytest.mark.parametrize(
    ("corruption", "message"),
    (
        ("catalog", "ambiguous biological Feature identity"),
        ("override", "conflicting derived Feature override"),
        ("caption", "conflicting saved legend state"),
        ("result", "conflicting catalogue and Result Feature fills"),
    ),
)
def test_v40_recovery_requires_unique_catalog_override_caption_and_result_evidence(
    tmp_path: Path,
    corruption: str,
    message: str,
) -> None:
    session = json.loads(_V40_STYLE_FIXTURE.read_text(encoding="utf-8"))
    if corruption == "catalog":
        features = session["editorState"]["featureCatalog"]["items"][0][
            "biologicalFeatures"
        ]
        features.append(copy.deepcopy(features[0]))
    elif corruption == "override":
        session["features"]["featureColorOverrides"]["file0_f8"][
            "caption"
        ] = "Wrong caption"
    elif corruption == "caption":
        session["editorState"]["legend"]["entries"][0]["color"] = "#000000"
    else:
        session["results"][0]["content"] = session["results"][0][
            "content"
        ].replace('fill="#577edb"', 'fill="#000000"')

    with materialize_session(session, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionConversionError, match=message):
            session_to_request(materialized)


def test_main_backed_v40_gallery_session_recovers_defaults_with_provenance(
    tmp_path: Path,
) -> None:
    expected = json.loads(_V40_STYLE_EXPECTED.read_text(encoding="utf-8"))
    source_bytes = _V40_GALLERY_SESSION.read_bytes()
    assert hashlib.sha256(source_bytes).hexdigest() == expected["sourceSha256"]

    with materialize_session(
        _V40_GALLERY_SESSION,
        output_directory=tmp_path,
    ) as materialized:
        decoded = session_to_request(materialized)
        assert decoded.options.colors is not None
        assert decoded.options.colors.color_table is not None
        assert len(decoded.options.colors.color_table) == 4
        default_file = decoded.options.colors.default_colors_file
        assert default_file is not None
        defaults = pd.read_csv(
            default_file,
            sep="\t",
            names=["feature_type", "color"],
            header=None,
            dtype=str,
        ).set_index("feature_type")["color"]
        assert defaults["CDS"] == "#dddddd"
        assert defaults["tRNA"] == "#9cd9cf"


def test_v40_tobacco_empty_editor_rules_keep_all_canonical_rules(
    tmp_path: Path,
) -> None:
    with materialize_session(
        _V40_TOBACCO_SESSION,
        output_directory=tmp_path,
    ) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.colors is not None
    assert decoded.options.colors.color_table is not None
    assert len(decoded.options.colors.color_table) == 71
    assert decoded.options.colors.color_table_file is None


def test_v32_web_slot_specs_drop_only_legacy_feature_geometry(
    tmp_path: Path,
) -> None:
    canonical_slots = (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3,legend_label=Genes",
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(
                linear_track_slots=canonical_slots,
                linear_track_axis_index=1,
            )
        ),
    )
    data = build_session_document(request).to_dict()
    data["version"] = 32
    data["renderRequest"]["schema"] = 2
    data["renderRequest"].pop("grouping", None)
    encoded_slots = data["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ]
    encoded_slots[1] = (
        "features:features@side=below,h=48px,spacing=9px,"
        "z=3,legend_label=Genes"
    )
    document = load_session_document(data)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.tracks.linear_track_slots == (
        "gc_content:dinucleotide_content@side=above,h=27px,spacing=4px,nt=GC",
        "features:features@side=below,z=3,legend_label=Genes",
    )
    assert decoded.options.tracks.linear_track_axis_index == 1
    assert document.to_dict()["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ] == data["renderRequest"]["diagramOptions"]["tracks"]["linearTrackSlots"]


def test_v32_structured_slots_preserve_non_feature_geometry_and_fields(
    tmp_path: Path,
) -> None:
    feature_slot = LinearTrackSlot(
        id="features",
        renderer="features",
        enabled=False,
        side="below",
        z=3,
        params={"legend_label": "Genes"},
    )
    numeric_slot = LinearTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        side="above",
        height=ScalarSpec(27, "px"),
        spacing=ScalarSpec(4, "px"),
        params={"nt": "GC"},
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(
                linear_track_slots=(numeric_slot, feature_slot),
                linear_track_axis_index=1,
            )
        ),
    )
    data = build_session_document(request).to_dict()
    data["version"] = 32
    data["renderRequest"]["schema"] = 2
    data["renderRequest"].pop("grouping", None)
    encoded_feature = data["renderRequest"]["diagramOptions"]["tracks"][
        "linearTrackSlots"
    ][1]
    encoded_feature["height"] = {"value": 48, "unit": "px"}
    encoded_feature["spacing"] = {"value": 9, "unit": "px"}

    with materialize_session(data, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    preserved_numeric, migrated_feature = decoded.options.tracks.linear_track_slots
    assert migrated_feature == LinearTrackSlot(
        id="features",
        renderer="features",
        enabled=False,
        side="below",
        height=None,
        spacing=None,
        z=3,
        params={"legend_label": "Genes"},
    )
    assert preserved_numeric == numeric_slot
    assert decoded.options.tracks.linear_track_axis_index == 1


def test_current_session_preserves_v2_feature_geometry(tmp_path: Path) -> None:
    feature_slot = LinearTrackSlot(
        id="features",
        renderer="features",
        side="overlay",
        height=ScalarSpec(48, "px"),
        spacing=ScalarSpec(9, "px"),
    )
    request = LinearDiagramRequest(
        records=(_record(),),
        options=LinearDiagramOptions(
            tracks=LinearTrackOptions(linear_track_slots=(feature_slot,))
        ),
    )
    document = build_session_document(request)

    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert decoded.options.tracks.linear_track_slots == (feature_slot,)


@pytest.mark.browser
def test_web_writer_payload_decodes_with_python_codec(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [node, "tests/web/session-request.test.mjs", "--print"],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        records = normalize_request_records(request)

    assert request.output.output_prefix == "web-session"
    assert [record.id for record in records] == ["WEBTEST"]


@pytest.mark.browser
def test_web_depth_writer_payload_decodes_with_python_codec(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [node, "tests/web/session-request.test.mjs", "--print-depth"],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)

    assert request.options.depth_tracks is not None
    assert [track.label for track in request.options.depth_tracks] == [
        "Sample A",
        "Sample B",
    ]
    assert request.options.depth_track_files is None
    assert request.options.depth_track_labels is None


@pytest.mark.browser
def test_web_resolved_protein_writer_preserves_alignment_settings(
    tmp_path: Path,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    completed = subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--print-resolved-protein",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
        }
    )

    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)

    assert isinstance(request, LinearDiagramRequest)
    assert request.options.protein_blastp_mode == "none"
    assert request.options.linear_comparisons is not None
    assert len(request.options.linear_comparisons) == 1
    assert request.options.align_orthogroup_feature == "resolved-feature-anchor"


@pytest.mark.parametrize("collinearity_value_kind", ("result", "blocks"))
@pytest.mark.browser
def test_python_typed_protein_round_trip_with_explicit_adjacent_comparison(
    tmp_path: Path,
    collinearity_value_kind: str,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    member = OrthogroupMember(
        orthogroup_id="OG1",
        protein_id="protein-a",
        record_index=0,
        feature_index=1,
        record_id="record-a",
        label="Protein A",
        start=10,
        end=40,
        strand=1,
        feature_svg_id="feature-a",
        source_protein_id="source-a",
    )
    orthogroups = OrthogroupResult(
        orthogroups={"OG1": [member]},
        member_by_protein_id={"protein-a": member},
        names_by_orthogroup_id={"OG1": "Example group"},
    )
    anchor = CollinearityAnchor(
        query_protein_id="protein-a",
        subject_protein_id="protein-b",
        query_record_index=0,
        subject_record_index=1,
        query_order=0,
        subject_order=1,
        query_start=10,
        query_end=40,
        subject_start=20,
        subject_end=50,
        identity=91.5,
        evalue=1e-20,
        bitscore=100.0,
        alignment_length=30,
        query_feature_svg_id="feature-a",
        subject_feature_svg_id="feature-b",
        source="precomputed",
        query_unit_id="unit-a",
        subject_unit_id="unit-b",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name="Protein A",
        subject_display_name="Protein B",
    )
    block = CollinearityBlock(
        block_id="block-1",
        query_record_index=0,
        subject_record_index=1,
        orientation="plus",
        score=100.0,
        anchors=(anchor,),
    )
    collinearity = (
        CollinearityResult(blocks=(block,), orthogroups=orthogroups)
        if collinearity_value_kind == "result"
        else (block,)
    )
    request = LinearDiagramRequest(
        records=(_record("record-a"), _record("record-b")),
        options=LinearDiagramOptions(
            orthogroups=orthogroups,
            collinearity_blocks=collinearity,
        ),
    )
    session_path = tmp_path / f"typed-{collinearity_value_kind}.json"
    save_session_document(session_path, request)

    completed = subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--round-trip-session",
            str(session_path),
            "--activate-adjacent-comparison",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parents[1],
    )
    canonical = json.loads(completed.stdout)
    comparisons = canonical["renderRequest"]["comparisons"]
    orthogroup_entry = next(
        item for item in comparisons if item["kind"] == "orthogroupResult"
    )
    collinearity_entry = next(
        item for item in comparisons if item["kind"] == "collinearityResult"
    )
    assert orthogroup_entry["encoding"] == "canonicalJson"
    assert collinearity_entry["encoding"] == "canonicalJson"
    assert collinearity_entry["valueKind"] == collinearity_value_kind

    document = load_session_document(
        {
            "format": "gbdraw-session",
            "version": CURRENT_SESSION_VERSION,
            "renderRequest": canonical["renderRequest"],
            "resources": canonical["resources"],
            "results": [],
            "editorState": {"featureCatalog": None},
        }
    )
    with materialize_session(document, output_directory=tmp_path) as materialized:
        decoded = session_to_request(materialized)

    assert isinstance(decoded, LinearDiagramRequest)
    assert decoded.options.orthogroups == orthogroups
    assert decoded.options.collinearity_blocks == collinearity


@pytest.mark.browser
def test_python_depth_request_projects_in_web(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    shared_depth = pd.DataFrame(
        {
            "reference_name": ["first", "second"],
            "position": [1, 1],
            "depth": [5, 8],
        }
    )
    sparse_depth = pd.DataFrame(
        {
            "reference_name": ["first"],
            "position": [1],
            "depth": [13],
        }
    )
    request = LinearDiagramRequest(
        records=(_record("first"), _record("second")),
        options=LinearDiagramOptions(
            depth_tracks=(
                DepthTrackInput(
                    source=shared_depth,
                    label="Shared",
                    color="#112233",
                    height=18,
                    large_tick_interval=10,
                ),
                DepthTrackInput(
                    source=(sparse_depth, None),
                    label="Sparse",
                    color="#445566",
                    height=24,
                    small_tick_interval=5,
                    tick_font_size=9,
                ),
            ),
        ),
    )
    session_path = tmp_path / "python-canonical-depth.gbdraw-session.json"
    save_session_document(session_path, request)
    session = json.loads(session_path.read_text(encoding="utf-8"))

    assert len(session["renderRequest"]["diagramOptions"]["depthTracks"]) == 2
    subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--project-session",
            str(session_path),
        ],
        check=True,
        cwd=Path(__file__).parents[1],
    )
