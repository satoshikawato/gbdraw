from __future__ import annotations

import base64
import copy
import json
import os
import re
from datetime import datetime
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw.analysis.protein_colinearity as protein_colinearity_module
import gbdraw.circular as circular_cli_module
import gbdraw.session_io as session_io_module
from gbdraw.analysis.protein_colinearity import (
    build_protein_losat_cache_key,
    build_protein_losat_pair_identity,
    extract_web_stable_cds_proteins,
    validate_protein_raw_entry_references,
)
from gbdraw.circular import circular_main
from gbdraw.linear import (
    _get_args as get_linear_args,
)
from gbdraw.cli_utils.session import (
    DiagramRunResult,
    collect_embedded_files_from_cli_args,
    make_rendered_svg,
    parse_session_pre_args,
    resolve_session_sidecar_path,
    save_session_sidecar_if_requested,
    strip_session_output_args,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import read_records_table
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.render.formats import ACCEPTED_FORMATS
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    DEPTH_FILE_ENCODING,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    SESSION_FORMAT,
    SUPPORTED_SESSION_VERSIONS,
    SessionBuildContext,
    _json_values_equal,
    build_session_json,
    compact_session_feature_catalog,
    decode_depth_payload,
    load_session,
    materialize_embedded_file,
    migrate_legacy_repeat_feature_shape_args,
    migrate_persisted_web_state_field_names,
    session_to_cli_args,
    validate_current_web_state_field_names,
    validate_session,
    write_session_json,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
from gbdraw.session import SessionFormatError, load_session_document


def _file_entry(name: str, content: bytes) -> dict:
    return {
        "name": name,
        "type": "application/octet-stream",
        "size": len(content),
        "lastModified": 0,
        "data": base64.b64encode(content).decode("ascii"),
    }


def _minimal_session(files: dict, *, mode: str = "circular") -> dict:
    return {
        "format": SESSION_FORMAT,
        "version": 30,
        "createdAt": "2026-06-22T00:00:00Z",
        "config": {"form": {"prefix": "out"}, "adv": {}},
        "ui": {"mode": mode, "cInputType": "gb", "lInputType": "gb"},
        "files": files,
    }


def _canonical_request(mode: str):
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    record_input = RecordInput(source=InMemoryRecordSource(record))
    if mode == "linear":
        return LinearDiagramRequest(records=(record_input,))
    return CircularDiagramRequest(records=(record_input,))


def _canonical_request_with_scale_visibility(mode: str, *, show: bool):
    record = SeqRecord(Seq("ATGC"), id="record", annotations={"molecule_type": "DNA"})
    record_input = RecordInput(source=InMemoryRecordSource(record))
    config = apply_config_overrides(None, {"objects.scale.show": show})
    if mode == "linear":
        return LinearDiagramRequest(
            records=(record_input,),
            options=LinearDiagramOptions(config=config),
        )
    return CircularDiagramRequest(
        records=(record_input,),
        options=CircularDiagramOptions(config=config),
    )


def test_write_session_json_does_not_use_predictable_temp_path(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "saved.gbdraw-session.json"
    old_temp_path = output_path.with_name(
        f".{output_path.name}.{os.getpid()}.tmp"
    )
    protected_path = tmp_path / "protected.txt"
    protected_path.write_text("keep", encoding="utf-8")
    try:
        old_temp_path.symlink_to(protected_path)
    except OSError as exc:
        pytest.skip(f"file symlinks are unavailable: {exc}")

    write_session_json(output_path, _minimal_session({}))

    assert protected_path.read_text(encoding="utf-8") == "keep"
    assert old_temp_path.is_symlink()
    assert output_path.is_file()


def test_write_session_json_no_replace_commit_preserves_late_target(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_path = tmp_path / "saved.gbdraw-session.json"
    real_commit = session_io_module.commit_staged_output_file

    def create_target_immediately_before_commit(
        staged_path: Path,
        target_path: Path,
        *,
        overwrite: bool,
    ) -> None:
        target_path.write_text("keep", encoding="utf-8")
        real_commit(staged_path, target_path, overwrite=overwrite)

    monkeypatch.setattr(
        session_io_module,
        "commit_staged_output_file",
        create_target_immediately_before_commit,
    )

    with pytest.raises(ValidationError, match="overwrite=True"):
        write_session_json(
            output_path,
            _minimal_session({}),
            overwrite=False,
        )

    assert output_path.read_text(encoding="utf-8") == "keep"
    assert not tuple(tmp_path.glob(f".{output_path.name}.*.tmp"))


def _protein_identity_manifest() -> dict:
    records = []
    for record_id, protein_id in (("A", "protein-a"), ("B", "protein-b")):
        record = SeqRecord(
            Seq("ATG"),
            id=record_id,
            annotations={"molecule_type": "DNA"},
        )
        record.features = [
            SeqFeature(
                FeatureLocation(0, 3, strand=1),
                type="CDS",
                qualifiers={"translation": ["M"], "protein_id": [protein_id]},
            )
        ]
        records.append(record)
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=("record-1", "record-2"),
    )
    assert extraction.identity_manifest is not None
    return extraction.identity_manifest.to_dict()


def _current_protein_cache_entry() -> dict:
    manifest = _protein_identity_manifest()
    feature_a = next(
        iter(manifest["recordInstances"]["record-1"]["runtimeIds"])
    )
    feature_b = next(
        iter(manifest["recordInstances"]["record-2"]["runtimeIds"])
    )
    query_runtime_handle = manifest["recordInstances"]["record-1"]["runtimeIds"][
        feature_a
    ]
    subject_runtime_handle = manifest["recordInstances"]["record-2"]["runtimeIds"][
        feature_b
    ]
    analysis_a = manifest["recordInstances"]["record-1"]["recordAnalysisId"]
    analysis_b = manifest["recordInstances"]["record-2"]["recordAnalysisId"]
    pair_identity = build_protein_losat_pair_identity(
        manifest,
        query_record_instance_key="record-1",
        subject_record_instance_key="record-2",
    )
    return {
        "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
        "kind": "raw-losat",
        "identityKind": "protein",
        "idEncoding": "runtime-handle-v1",
        "key": build_protein_losat_cache_key(pair_identity, args=[]),
        "text": (
            f"{query_runtime_handle}\t{subject_runtime_handle}\t"
            "100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n"
        ),
        "program": "blastp",
        "outfmt": "6",
        "args": [],
        "queryProteinSetHash": manifest["recordAnalyses"][analysis_a][
            "proteinSetHash"
        ],
        "subjectProteinSetHash": manifest["recordAnalyses"][analysis_b][
            "proteinSetHash"
        ],
        "queryRuntimeBindingHash": manifest["recordInstances"]["record-1"][
            "runtimeBindingHash"
        ],
        "subjectRuntimeBindingHash": manifest["recordInstances"]["record-2"][
            "runtimeBindingHash"
        ],
        "queryRecordInstanceKey": "record-1",
        "subjectRecordInstanceKey": "record-2",
    }


def _current_derived_session(
    manifest: dict,
    derived_payload: dict,
    *,
    mode: str = "collinear",
) -> dict:
    return build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        losat_derived_cache_entries=(
            {
                "schema": LOSAT_DERIVED_CACHE_SCHEMA,
                "kind": "derived-losatp-payload",
                "idEncoding": "runtime-handle-v1",
                "key": "derived-key",
                "mode": mode,
                "payload": derived_payload,
            },
        ),
        protein_identity_manifest=manifest,
        canonical_request=_canonical_request("linear"),
    )


def _zero_hit_derived_payload(
    mode: str,
    *,
    include_identity: bool = True,
) -> dict:
    payload = {
        "pairs": [
            {
                "pair_index": 0,
                "query_index": 0,
                "subject_index": 1,
                "tsv": "",
                "rows": [],
                "hit_count": 0,
            }
        ],
        "orthogroups": [],
    }
    if include_identity:
        payload["identity"] = {
            "cacheSchema": LOSAT_DERIVED_CACHE_SCHEMA,
            "idEncoding": "runtime-handle-v1",
            "converter": "convert_losatp_blastp_pairs_to_genomic_payload",
            "mode": mode,
            "rawCacheKeys": ["raw-key"],
        }
    if mode == "collinear":
        payload.update(
            {
                "collinearGroups": [],
                "collinearGroupScope": "adjacent_local",
                "collinearityBlocks": [],
            }
        )
    return payload


_LOSAT_OUTFMT6_NUMERIC_CASES = json.loads(
    (
        Path(__file__).parent
        / "fixtures"
        / "losat-outfmt6-numeric-contract.json"
    ).read_text(encoding="utf-8")
)


def _protein_raw_entry_with_numeric_case(entry: dict, case: dict) -> dict:
    columns = entry["text"].rstrip("\r\n").split("\t")
    columns[COMPARISON_COLUMNS.index(case["column"])] = case["value"]
    return {**entry, "text": "\t".join(columns) + "\n"}


@pytest.mark.parametrize(
    "case",
    _LOSAT_OUTFMT6_NUMERIC_CASES,
    ids=[case["name"] for case in _LOSAT_OUTFMT6_NUMERIC_CASES],
)
def test_current_protein_raw_numeric_contract(
    case: dict,
) -> None:
    current_manifest = _protein_identity_manifest()
    current_entry = _protein_raw_entry_with_numeric_case(
        _current_protein_cache_entry(),
        case,
    )
    assert (
        validate_protein_raw_entry_references(current_entry, current_manifest)
        is case["valid"]
    )


def test_current_session_materializes_manifest_once_for_multiple_protein_entries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = _protein_identity_manifest()
    first_entry = _current_protein_cache_entry()
    pair_identity = build_protein_losat_pair_identity(
        manifest,
        query_record_instance_key="record-1",
        subject_record_instance_key="record-2",
    )
    second_args = ["--max-hsps-per-subject", "1"]
    second_entry = {
        **first_entry,
        "args": second_args,
        "key": build_protein_losat_cache_key(
            pair_identity,
            args=second_args,
        ),
    }
    real_validate = (
        protein_colinearity_module.validate_protein_identity_manifest
    )
    validation_calls = 0

    def count_validation(value):
        nonlocal validation_calls
        validation_calls += 1
        return real_validate(value)

    monkeypatch.setattr(
        protein_colinearity_module,
        "validate_protein_identity_manifest",
        count_validation,
    )

    session_io_module.validate_current_session_artifacts(
        {
            "losatCache": {"entries": [first_entry, second_entry]},
            "losatDerivedCache": {"entries": []},
            "proteinIdentityManifest": manifest,
        }
    )

    assert validation_calls == 1


def _legacy_protein_cache_entry() -> dict:
    return {
        "schema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
        "kind": "raw-losat",
        "key": "legacy-protein-key",
        "text": "p_r_old\tp_r_other\n",
        "program": "blastp",
        "outfmt": "6",
        "args": [],
        "queryCanonicalHash": "old-query",
        "subjectCanonicalHash": "old-subject",
    }


def test_session_sidecar_saves_complete_orthogroup_state(tmp_path: Path) -> None:
    output_prefix = tmp_path / "diagram"
    rendered_svg = make_rendered_svg(str(output_prefix))
    rendered_svg.svg_path.write_text("<svg></svg>", encoding="utf-8")
    sidecar = tmp_path / "diagram.gbdraw-session.json"
    run_result = DiagramRunResult(
        mode="linear",
        render_formats=("interactive_svg",),
        outputs=(rendered_svg,),
        feature_metadata=({"svg_id": "feature-1", "orthogroup_id": "og_1"},),
        biological_feature_metadata=(
            {
                "svg_id": "feature-1",
                "stable_feature_id": "feature-1",
                "record_idx": 0,
                "nucleotide_sequence": "ATGC",
            },
            {
                "svg_id": "feature-2",
                "stable_feature_id": "feature-2",
                "record_idx": 1,
                "nucleotide_sequence": "ATGA",
            },
        ),
        orthogroup_metadata=(
            {
                "id": "og_1",
                "member_count": 2,
                "members": [
                    {"featureSvgId": "feature-1", "proteinId": "protein-1"},
                    {"featureSvgId": "feature-2", "proteinId": "protein-2"},
                ],
            },
        ),
        losat_derived_cache_entries=({"schema": LOSAT_DERIVED_CACHE_SCHEMA},),
        canonical_request=_canonical_request("linear"),
    )

    saved = save_session_sidecar_if_requested(
        save_session=True,
        session_output=str(sidecar),
        output_prefix=str(output_prefix),
        run_result=run_result,
        cmd_args=(),
    )

    assert saved == sidecar
    payload = load_session(sidecar)
    assert payload["features"] == {}
    assert payload["orthogroupState"] == {}
    catalog = payload["editorState"]["featureCatalog"]
    assert catalog["schema"] == 3
    item = catalog["items"][0]
    assert item["resultIndex"] == 0
    assert item["resultName"] == "diagram"
    assert [
        feature.get("stableFeatureId", feature["biologicalFeatureId"])
        for feature in item["biologicalFeatures"]
    ] == ["feature-1", "feature-2"]
    group = item["orthogroups"][0]
    assert group["id"] == "og_1"
    assert {
        member["biologicalFeatureId"] for member in group["members"]
    } == {"feature-1", "feature-2"}
    assert payload["losatDerivedCache"] == {"entries": []}

    sidecar.write_text("keep this session", encoding="utf-8")
    with pytest.raises(ValidationError, match="--overwrite"):
        save_session_sidecar_if_requested(
            save_session=True,
            session_output=str(sidecar),
            output_prefix=str(output_prefix),
            run_result=run_result,
            cmd_args=(),
        )
    assert sidecar.read_text(encoding="utf-8") == "keep this session"

    save_session_sidecar_if_requested(
        save_session=True,
        session_output=str(sidecar),
        output_prefix=str(output_prefix),
        run_result=run_result,
        cmd_args=(),
        overwrite=True,
    )
    assert load_session(sidecar)["format"] == "gbdraw-session"


def test_current_session_version_matches_web_config() -> None:
    source = Path("gbdraw/web/js/services/config.js").read_text(encoding="utf-8")
    match = re.search(r"const\s+SESSION_VERSION\s*=\s*(\d+);", source)
    assert match is not None
    supported_match = re.search(
        r"SUPPORTED_SESSION_VERSIONS\s*=\s*new Set\(\[([\s\S]*?)\]\)",
        source,
    )
    assert supported_match is not None
    web_supported_versions = {
        int(value) for value in re.findall(r"\b\d+\b", supported_match.group(1))
    }
    if "SESSION_VERSION" in supported_match.group(1):
        web_supported_versions.add(int(match.group(1)))

    assert CURRENT_SESSION_VERSION == 40
    assert SUPPORTED_SESSION_VERSIONS == frozenset(
        {27, 28, 29, 30, 31, 32, 33, 39, CURRENT_SESSION_VERSION}
    )
    assert int(match.group(1)) == CURRENT_SESSION_VERSION
    assert web_supported_versions == SUPPORTED_SESSION_VERSIONS

    authority_source = Path(
        "gbdraw/web/js/services/session-authority.js"
    ).read_text(encoding="utf-8")
    forbidden_match = re.search(
        r"CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS\s*=\s*Object\.freeze\(\["
        r"([\s\S]*?)\]\)",
        authority_source,
    )
    assert forbidden_match is not None
    web_forbidden_fields = frozenset(
        re.findall(r"'([^']+)'", forbidden_match.group(1))
    )
    assert session_io_module.CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS == (
        web_forbidden_fields
    )
    authority_match = re.search(
        r"SESSION_TOP_LEVEL_AUTHORITY\s*=\s*Object\.freeze\(\{"
        r"([\s\S]*?)\}\);",
        authority_source,
    )
    assert authority_match is not None
    web_authority_fields = frozenset(
        re.findall(
            r"^\s{2}([A-Za-z][A-Za-z0-9]*):",
            authority_match.group(1),
            flags=re.MULTILINE,
        )
    )
    assert session_io_module.CURRENT_SESSION_TOP_LEVEL_FIELDS == (
        web_authority_fields
    )


def test_current_session_feature_catalog_is_single_and_lossless(
    tmp_path: Path,
) -> None:
    catalog = {
        "schema": 3,
        "items": [
            {
                "resultIndex": 0,
                "resultName": "out",
                "recordKeys": ["record-key"],
                "features": [
                    {
                        "svgId": "rendered-feature",
                        "recordKey": "record-key",
                        "biologicalFeatureId": "feature-1",
                        "fillColor": "#123456",
                    }
                ],
                "biologicalFeatures": [
                    {
                        "recordKey": "record-key",
                        "biologicalFeatureId": "feature-1",
                        "stableFeatureId": "stable-feature",
                        "record_id": "public-record",
                        "type": "CDS",
                        "start": 1,
                        "end": 3,
                        "qualifiers": {"product": ["example"]},
                        "nucleotide_sequence": "ATG",
                        "amino_acid_sequence": "M",
                    }
                ],
                "orthogroups": [],
                "annotations": [],
                "comparisonMatches": [],
            }
        ],
    }
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 24),
        feature_catalog=catalog,
        canonical_request=_canonical_request("linear"),
    )
    payload["features"]["selectedFeatureRecordIdx"] = 0
    path = tmp_path / "schema-3-feature-catalog.gbdraw-session.json"

    write_session_json(path, payload)

    on_disk = json.loads(path.read_text(encoding="utf-8"))
    assert on_disk["editorState"]["featureCatalog"] == catalog
    assert on_disk["features"] == {"selectedFeatureRecordIdx": 0}
    assert on_disk["orthogroupState"] == {}
    assert load_session(path) == on_disk
    assert load_session_document(path).to_dict() == on_disk

    unknown_field = copy.deepcopy(on_disk)
    unknown_field["branchOnlyState"] = {}
    with pytest.raises(
        ValidationError,
        match="unclassified top-level field.*branchOnlyState",
    ):
        write_session_json(tmp_path / "unknown-field.json", unknown_field)

    missing_results = copy.deepcopy(on_disk)
    missing_results.pop("results")
    with pytest.raises(
        ValidationError,
        match="requires a results array",
    ):
        write_session_json(tmp_path / "missing-results.json", missing_results)

    empty_draft = copy.deepcopy(on_disk)
    empty_draft["results"] = []
    empty_draft["editorState"]["featureCatalog"] = None
    write_session_json(tmp_path / "empty-draft.json", empty_draft)

    missing_empty_catalog = copy.deepcopy(empty_draft)
    missing_empty_catalog["editorState"].pop("featureCatalog")
    with pytest.raises(
        ValidationError,
        match="requires editorState.featureCatalog",
    ):
        write_session_json(
            tmp_path / "missing-empty-catalog.json",
            missing_empty_catalog,
        )

    for field in ("features", "orthogroupState"):
        invalid_container = copy.deepcopy(on_disk)
        invalid_container[field] = []
        with pytest.raises(
            ValidationError,
            match=rf"Session {field} must be an object",
        ):
            write_session_json(
                tmp_path / f"invalid-{field}.json",
                invalid_container,
            )

    interactive_result = copy.deepcopy(on_disk)
    interactive_result["results"][0]["name"] = "out.interactive.svg"
    interactive_result["editorState"]["featureCatalog"]["items"][0][
        "resultName"
    ] = "out.interactive.svg"
    with pytest.raises(
        ValidationError,
        match="named plain SVG",
    ):
        write_session_json(
            tmp_path / "interactive-result.json",
            interactive_result,
        )

    missing_catalog = copy.deepcopy(on_disk)
    missing_catalog["editorState"].pop("featureCatalog")
    with pytest.raises(
        ValidationError,
        match="Results require editorState.featureCatalog",
    ):
        write_session_json(
            tmp_path / "missing-catalog.json",
            missing_catalog,
        )

    duplicate_feature = copy.deepcopy(on_disk)
    duplicate_feature["editorState"]["featureCatalog"]["items"][0][
        "biologicalFeatures"
    ].append(
        copy.deepcopy(
            duplicate_feature["editorState"]["featureCatalog"]["items"][0][
                "biologicalFeatures"
            ][0]
        )
    )
    with pytest.raises(
        ValidationError,
        match="invalid biological feature reference",
    ):
        write_session_json(
            tmp_path / "duplicate-feature.json",
            duplicate_feature,
        )

    duplicated_legacy = copy.deepcopy(on_disk)
    duplicated_legacy["features"]["biologicalFeatures"] = []
    with pytest.raises(
        ValidationError,
        match="derived feature payloads",
    ):
        write_session_json(
            tmp_path / "duplicated-legacy.json",
            duplicated_legacy,
        )

    for field in ("featureSelectorSafetyScope", "featureRecordIds"):
        duplicated_legacy = copy.deepcopy(on_disk)
        duplicated_legacy["features"][field] = []
        with pytest.raises(
            ValidationError,
            match="derived feature payloads",
        ):
            write_session_json(
                tmp_path / f"duplicated-{field}.json",
                duplicated_legacy,
            )


def test_feature_catalog_compaction_falls_back_for_partial_metadata() -> None:
    biological = {
        "id": "f0",
        "svg_id": "f_stable",
        "record_idx": 0,
        "feature_index": 0,
        "qualifiers": {},
        "selector": {},
    }
    extracted = copy.deepcopy(biological)
    extracted.pop("feature_index")
    session = {
        "version": CURRENT_SESSION_VERSION,
        "features": {
            "extractedFeatures": [extracted],
            "biologicalFeatures": [biological],
        },
    }

    compact = compact_session_feature_catalog(session)

    assert compact is session
    assert "featureCatalog" not in compact["features"]


def test_feature_catalog_equality_preserves_json_scalar_types() -> None:
    assert not _json_values_equal(
        {"custom_payload": [1]},
        {"custom_payload": [True]},
    )


def test_current_session_validates_mixed_protein_and_nucleotide_raw_cache() -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        losat_cache_entries=(
            _current_protein_cache_entry(),
            {
                "schema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
                "kind": "raw-losat",
                "identityKind": "nucleotide",
                "key": "nucleotide-key",
                "text": "",
                "program": "blastn",
                "outfmt": "6",
                "args": [],
                "queryCanonicalHash": "query",
                "subjectCanonicalHash": "subject",
            },
        ),
        losat_derived_cache_entries=(
            {
                "schema": LOSAT_DERIVED_CACHE_SCHEMA,
                "kind": "derived-losatp-payload",
                "idEncoding": "runtime-handle-v1",
                "key": "derived-key",
                "mode": "orthogroup",
                "payload": {},
            },
        ),
        protein_identity_manifest=_protein_identity_manifest(),
        canonical_request=_canonical_request("linear"),
    )

    validate_session(payload)
    assert [entry["schema"] for entry in payload["losatCache"]["entries"]] == [
        PROTEIN_LOSAT_CACHE_SCHEMA,
        NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    ]
    assert payload["losatDerivedCache"]["entries"][0]["schema"] == 3


@pytest.mark.parametrize(
    ("config", "message"),
    (
        (
            {"adv": {"depth_tick_interval": 10}},
            "Web state field adv.depth_tick_interval is obsolete; "
            "use adv.depth_large_tick_interval.",
        ),
        (
            {"adv": {"depth_tracks": [{"tick_interval": 10}]}},
            "Web state field adv.depth_tracks[0].tick_interval is obsolete; "
            "use large_tick_interval.",
        ),
        (
            {"losat": {"blastp": {"collinearMaxGeneGap": 2}}},
            "Web state field losat.blastp.collinearMaxGeneGap is obsolete; "
            "use losat.blastp.collinearMaxUnitGap.",
        ),
    ),
)
def test_current_session_rejects_obsolete_web_state_field_names(
    tmp_path: Path,
    config: dict,
    message: str,
) -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 30),
        canonical_request=_canonical_request("linear"),
    )
    payload["config"] = config

    with pytest.raises(ValidationError, match=re.escape(message)):
        validate_session(payload)
    output_path = tmp_path / "obsolete-current-session.gbdraw-session.json"
    with pytest.raises(ValidationError, match=re.escape(message)):
        write_session_json(output_path, payload)
    assert not output_path.exists()

    payload["version"] = 33
    payload["renderRequest"]["schema"] = 2
    validate_session(payload)


def test_released_web_state_field_names_migrate_copy_on_write() -> None:
    source = {
        "adv": {
            "depth_tick_interval": 10,
            "depth_large_tick_interval": 20,
            "depth_tracks": [
                {"tick_interval": 5},
                {"tick_interval": 6, "large_tick_interval": 7},
            ],
        },
        "losat": {
            "blastp": {
                "collinearMaxGeneGap": 2,
                "collinearMaxUnitGap": 3,
            }
        },
    }

    migrated = migrate_persisted_web_state_field_names(source)

    assert isinstance(migrated, dict)
    assert migrated["adv"]["depth_large_tick_interval"] == 20
    assert "depth_tick_interval" not in migrated["adv"]
    assert migrated["adv"]["depth_tracks"] == [
        {"large_tick_interval": 5},
        {"large_tick_interval": 7},
    ]
    assert migrated["losat"]["blastp"]["collinearMaxUnitGap"] == 3
    assert "collinearMaxGeneGap" not in migrated["losat"]["blastp"]
    assert source["adv"]["depth_tick_interval"] == 10
    assert source["adv"]["depth_tracks"][0]["tick_interval"] == 5
    assert source["losat"]["blastp"]["collinearMaxGeneGap"] == 2
    validate_current_web_state_field_names(migrated)


@pytest.mark.parametrize(
    ("mode", "include_identity"),
    (
        ("orthogroup", True),
        ("collinear", True),
        ("collinear", False),
    ),
)
def test_current_session_accepts_strict_zero_hit_derived_results(
    mode: str,
    include_identity: bool,
) -> None:
    manifest = _protein_identity_manifest()

    payload = _current_derived_session(
        manifest,
        _zero_hit_derived_payload(
            mode,
            include_identity=include_identity,
        ),
        mode=mode,
    )

    validate_session(payload)


@pytest.mark.parametrize(
    "mutate",
    (
        lambda payload: payload.update({"note": "arbitrary"}),
        lambda payload: payload["pairs"][0].update({"rows": [{}]}),
        lambda payload: payload["pairs"][0].update({"tsv": "unexpected"}),
        lambda payload: payload["pairs"][0].update({"hit_count": 1}),
        lambda payload: payload["pairs"][0].update({"pair_index": "0"}),
        lambda payload: payload["pairs"][0].pop("subject_index"),
        lambda payload: payload.update({"pairs": []}),
        lambda payload: payload.update({"orthogroups": [{}]}),
        lambda payload: payload.update({"collinearGroups": [{}]}),
        lambda payload: payload.update({"collinearityBlocks": [{}]}),
        lambda payload: payload.update({"collinearGroupScope": "invalid"}),
        lambda payload: payload.update({"identity": None}),
        lambda payload: payload["identity"].update({"rawCacheKeys": [None]}),
        lambda payload: payload["identity"].update({"mode": "orthogroup"}),
    ),
    ids=(
        "arbitrary-field",
        "nonempty-rows",
        "nonempty-tsv",
        "nonzero-hit-count",
        "noninteger-pair-index",
        "partial-record-indices",
        "missing-pair-results",
        "nonempty-orthogroups",
        "nonempty-collinear-groups",
        "nonempty-collinearity-blocks",
        "invalid-collinear-scope",
        "null-identity",
        "invalid-raw-cache-binding",
        "identity-mode-mismatch",
    ),
)
def test_current_session_rejects_near_miss_zero_hit_derived_results(mutate) -> None:
    manifest = _protein_identity_manifest()
    derived_payload = _zero_hit_derived_payload("collinear")
    mutate(derived_payload)

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            derived_payload,
            mode="collinear",
        )


def test_current_derived_cache_accepts_compound_collinear_runtime_references() -> None:
    manifest = _protein_identity_manifest()
    runtime_handles = [
        runtime_handle
        for binding in manifest["recordInstances"].values()
        for runtime_handle in binding["runtimeIds"].values()
    ]
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        losat_derived_cache_entries=(
            {
                "schema": LOSAT_DERIVED_CACHE_SCHEMA,
                "kind": "derived-losatp-payload",
                "idEncoding": "runtime-handle-v1",
                "key": "collinear-derived-key",
                "mode": "collinear",
                "payload": {
                    "pairs": [
                        {
                            "rows": [
                                {
                                    "query_protein_id": ";".join(runtime_handles),
                                    "subject_protein_id": runtime_handles[-1],
                                    "supportingEdges": [
                                        f"{runtime_handles[0]}->"
                                        f"{runtime_handles[-1]}:rbh"
                                    ],
                                    "edgeIds": [
                                        f"og_1:0:{runtime_handles[0]}->"
                                        f"1:{runtime_handles[-1]}:coortholog"
                                    ],
                                    "description": (
                                        "The p_r_ prefix is documentation, not an ID."
                                    ),
                                }
                            ]
                        }
                    ]
                },
            },
        ),
        protein_identity_manifest=manifest,
        canonical_request=_canonical_request("linear"),
    )

    validate_session(payload)
    payload["losatDerivedCache"]["entries"][0]["payload"]["pairs"][0]["rows"][0][
        "query_protein_id"
    ] += ";h_aaaaaaaaaaaaaaaaaaaaaaaaaa"
    with pytest.raises(ValidationError, match="unresolved protein references"):
        validate_session(payload)


@pytest.mark.parametrize(
    "field_name",
    (
        "queryUnitId",
        "subjectUnitId",
        "query_unit_id",
        "subject_unit_id",
    ),
)
def test_current_derived_cache_validates_runtime_handles_in_unit_ids(
    field_name: str,
) -> None:
    manifest = _protein_identity_manifest()
    runtime_handles = [
        runtime_handle
        for binding in manifest["recordInstances"].values()
        for runtime_handle in binding["runtimeIds"].values()
    ]
    valid = _current_derived_session(
        manifest,
        {
            "proteinId": runtime_handles[0],
            field_name: ";".join(runtime_handles),
        },
    )
    validate_session(valid)

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            {
                "proteinId": runtime_handles[0],
                field_name: (
                    f"{runtime_handles[0]};h_aaaaaaaaaaaaaaaaaaaaaaaaaa"
                ),
            },
        )


def test_current_derived_cache_allows_non_protein_collinearity_unit_ids() -> None:
    manifest = _protein_identity_manifest()
    runtime_handle = next(
        iter(manifest["recordInstances"]["record-1"]["runtimeIds"].values())
    )

    payload = _current_derived_session(
        manifest,
        {
            "proteinId": runtime_handle,
            "query_unit_id": "gbd_r0001_unit000001",
            "subject_unit_id": "gbd_r0002_unit000002",
        },
    )

    validate_session(payload)


@pytest.mark.parametrize("field_name", ["supportingEdges", "edgeIds"])
def test_current_derived_cache_rejects_unresolved_compound_edge_references(
    field_name: str,
) -> None:
    manifest = _protein_identity_manifest()
    runtime_handles = [
        runtime_handle
        for binding in manifest["recordInstances"].values()
        for runtime_handle in binding["runtimeIds"].values()
    ]
    unknown = "h_aaaaaaaaaaaaaaaaaaaaaaaaaa"
    edge_reference = (
        f"{runtime_handles[0]}->{unknown}:rbh"
        if field_name == "supportingEdges"
        else f"og_1:0:{runtime_handles[0]}->1:{unknown}:rbh"
    )

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            {
                "proteinId": runtime_handles[0],
                field_name: [edge_reference],
            },
        )


@pytest.mark.parametrize(
    ("field_name", "invalid_value"),
    [
        ("supportingEdges", 7),
        ("supporting_edges", [7]),
        ("edgeIds", {"edge": "invalid"}),
        ("edge_ids", [None]),
    ],
)
def test_current_derived_cache_rejects_non_string_compound_edge_references(
    field_name: str,
    invalid_value: object,
) -> None:
    manifest = _protein_identity_manifest()
    runtime_handle = next(
        iter(manifest["recordInstances"]["record-1"]["runtimeIds"].values())
    )

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            {
                "proteinId": runtime_handle,
                field_name: invalid_value,
            },
        )


@pytest.mark.parametrize(
    "legacy_reference",
    [
        "p_r_old_0_9_1_deadbeefdead",
        f"f_{'b' * 64}",
    ],
)
def test_current_derived_cache_rejects_embedded_legacy_references(
    legacy_reference: str,
) -> None:
    manifest = _protein_identity_manifest()
    runtime_handle = next(
        iter(manifest["recordInstances"]["record-1"]["runtimeIds"].values())
    )

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            {
                "proteinId": runtime_handle,
                "description": f"retained evidence: {legacy_reference}",
            },
        )


def test_current_session_rejects_legacy_protein_entry_in_current_cache() -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        canonical_request=_canonical_request("linear"),
    )
    payload["losatCache"]["entries"] = [_legacy_protein_cache_entry()]

    with pytest.raises(ValidationError, match="matching legacyArtifacts"):
        validate_session(payload)


def test_current_writer_requires_typed_request_to_promote_legacy_schema() -> None:
    source = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="old",
            render_formats=("svg",),
        ),
        svg_results=(("old", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 20),
        canonical_request=_canonical_request("linear"),
    )
    source["version"] = 33
    source["renderRequest"]["schema"] = 2
    source["config"]["adv"]["depth_tick_interval"] = 10
    source["config"]["adv"]["depth_tracks"] = [{"tick_interval": 5}]
    source["config"]["losat"] = {
        "blastp": {"collinearMaxGeneGap": 2}
    }

    context = SessionBuildContext(
        mode="linear",
        output_prefix="new",
        render_formats=("svg",),
        source_session=source,
    )
    kwargs = {
        "svg_results": (("new", "<svg></svg>"),),
        "embedded_files": {"linearSeqs": []},
        "generated_at": datetime(2026, 7, 21),
    }

    with pytest.raises(
        ValidationError,
        match="canonical typed request is required",
    ):
        build_session_json(context, **kwargs)

    promoted = build_session_json(
        context,
        canonical_request=_canonical_request("linear"),
        **kwargs,
    )

    assert promoted["version"] == CURRENT_SESSION_VERSION
    assert promoted["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert promoted["config"]["adv"]["depth_large_tick_interval"] == 10
    assert promoted["config"]["adv"]["depth_tracks"] == [
        {"large_tick_interval": 5}
    ]
    assert (
        promoted["config"]["losat"]["blastp"]["collinearMaxUnitGap"]
        == 2
    )


def test_version_39_writer_promotes_once_and_preserves_web_inventory() -> None:
    fixture_path = (
        Path(__file__).parent
        / "fixtures"
        / "sessions"
        / "BGC0000708-BGC0000713.v39.gbdraw-session.json.gz"
    )
    source = load_session(fixture_path)
    assert source["version"] == 39
    source["webFiles"] = {"linearRecords": [{"uid": "record-1"}]}
    source["editorState"] = {"legend": {"entries": []}}
    source["features"] = {
        "featureSelectorSafetyScope": "all-records",
        "featureRecordIds": ["record-1"],
        "selectedFeatureRecordIdx": 0,
    }
    source["config"]["adv"]["circular_track_slots_enabled"] = False
    source["config"]["adv"]["circular_track_slots"] = [
        {"id": "circular-disabled", "enabled": False}
    ]
    source["config"]["adv"]["linear_track_slots_enabled"] = True
    source["config"]["adv"]["linear_track_slots"] = [
        {"id": "linear-active", "enabled": True}
    ]
    source["config"]["modeProfiles"] = {
        "schema": 1,
        "activeMode": "linear",
        "profiles": {
            "circular": {
                "values": {"identity": 88},
                "managed": {"identity": False},
            },
            "linear": {
                "values": {"identity": 77},
                "managed": {"identity": False},
            },
        },
    }

    promoted = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="new",
            render_formats=("svg",),
            source_session=source,
        ),
        svg_results=(("new", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 31),
        canonical_request=_canonical_request("linear"),
    )
    rewritten = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="again",
            render_formats=("svg",),
            source_session=promoted,
        ),
        svg_results=(("again", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 8, 1),
        canonical_request=_canonical_request("linear"),
    )

    assert source["version"] == 39
    assert promoted["version"] == CURRENT_SESSION_VERSION
    assert rewritten["version"] == CURRENT_SESSION_VERSION
    for payload in (promoted, rewritten):
        assert "files" not in payload
        assert payload["webFiles"]["linearRecords"] == source["webFiles"][
            "linearRecords"
        ]
        assert payload["webFiles"]["bindings"]["schema"] == 1
        assert "linearCanonicalComparisons" not in payload["webFiles"]["bindings"]
        assert payload["config"]["linearComparisonPlan"] == {
            "mode": "adjacent",
            "defaultSource": "losat",
            "edges": [],
        }
        assert "blastSource" not in payload["config"]
        assert "comparisons" not in payload["config"]["linearRecordLayout"]
        assert payload["editorState"]["legend"] == source["editorState"]["legend"]
        assert payload["editorState"]["featureCatalog"]["schema"] == 3
        assert len(payload["editorState"]["featureCatalog"]["items"]) == 1
        assert payload["features"] == {"selectedFeatureRecordIdx": 0}
        assert (
            payload["config"]["adv"]["circular_track_slots"]
            == source["config"]["adv"]["circular_track_slots"]
        )
        assert (
            payload["config"]["adv"]["linear_track_slots"]
            == source["config"]["adv"]["linear_track_slots"]
        )
        assert payload["config"]["modeProfiles"] == source["config"]["modeProfiles"]


def test_genuine_version_39_shape_retains_disabled_layout_upload_draft() -> None:
    fixture_path = (
        Path(__file__).parent
        / "fixtures"
        / "sessions"
        / "BGC0000708-BGC0000713.v39.gbdraw-session.json.gz"
    )
    source = load_session(fixture_path)
    assert source["version"] == 39
    assert source["config"]["linearRecordLayout"]["enabled"] is False
    record_uids = [
        str(row["uid"])
        for row in source["config"]["linearRecordLayout"]["rows"]
    ]
    dormant_file = _file_entry("retained-non-adjacent.tsv", b"a\tc\t97\n")

    migrated_config, migrated_files = (
        session_io_module.migrate_legacy_linear_comparison_draft_for_current_writer(
            source["config"],
            {
                "linearSeqs": [{"uid": uid} for uid in record_uids],
                "linearComparisons": [
                    {
                        "id": "dormant-v39-upload",
                        "queryUid": record_uids[0],
                        "subjectUid": record_uids[2],
                        "source": "upload",
                        "file": dormant_file,
                    }
                ],
            },
            force_web_draft=True,
        )
    )

    assert migrated_config["linearComparisonPlan"] == {
        "mode": "adjacent",
        "defaultSource": "losat",
        "edges": [
            {
                "id": "dormant-v39-upload",
                "queryUid": record_uids[0],
                "subjectUid": record_uids[2],
                "included": False,
                "fileActive": False,
                "losatFilenameActive": False,
                "source": "upload",
                "losatFilename": "",
            }
        ],
    }
    assert migrated_files["linearComparisons"] == [
        {"id": "dormant-v39-upload", "file": dormant_file}
    ]


def test_pre40_web_comparison_draft_migrates_directly_to_final_plan() -> None:
    source = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="old",
            render_formats=("svg",),
        ),
        svg_results=(("old", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 20),
        canonical_request=_canonical_request("linear"),
    )
    source["version"] = 33
    source["config"]["blastSource"] = "upload"
    source["config"]["linearRecordLayout"] = {
        "enabled": True,
        "recordGap": 18,
        "rows": [{"uid": "record-a", "row": 1}, {"uid": "record-b", "row": 2}],
        "comparisons": [
            {
                "id": "selected-a-b",
                "queryUid": "record-a",
                "subjectUid": "record-b",
                "source": "upload",
            }
        ],
    }
    comparison_file = _file_entry("selected.tsv", b"query\tsubject\t99\n")

    promoted = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="new",
            render_formats=("svg",),
            source_session=source,
        ),
        svg_results=(("new", "<svg></svg>"),),
        embedded_files={
            "linearSeqs": [
                {
                    "uid": "record-a",
                    "blast": comparison_file,
                    "losat_filename": "custom-a.fna",
                },
                {"uid": "record-b"},
            ],
            "linearComparisons": [
                {
                    "id": "selected-a-b",
                    "queryUid": "record-a",
                    "subjectUid": "record-b",
                    "source": "upload",
                    "file": comparison_file,
                }
            ],
        },
        generated_at=datetime(2026, 7, 21),
        canonical_request=_canonical_request("linear"),
    )

    assert "blastSource" not in promoted["config"]
    assert "comparisons" not in promoted["config"]["linearRecordLayout"]
    assert promoted["config"]["linearComparisonPlan"] == {
        "mode": "selected",
        "defaultSource": "upload",
        "edges": [
            {
                "id": "selected-a-b",
                "queryUid": "record-a",
                "subjectUid": "record-b",
                "included": True,
                "fileActive": True,
                "losatFilenameActive": True,
                "source": "upload",
                "losatFilename": "custom-a.fna",
            }
        ],
    }
    sequence_binding = promoted["webFiles"]["bindings"]["linearSeqs"][0]
    assert "blast" not in sequence_binding
    assert "losat_filename" not in sequence_binding
    comparison_binding = promoted["webFiles"]["bindings"]["linearComparisons"]
    assert [set(binding) for binding in comparison_binding] == [{"id", "file"}]
    resource_id = comparison_binding[0]["file"]["resourceId"]
    assert base64.b64decode(promoted["resources"][resource_id]["data"]) == (
        b"query\tsubject\t99\n"
    )
    matching_payloads = [
        resource
        for resource in promoted["resources"].values()
        if base64.b64decode(resource["data"]) == b"query\tsubject\t99\n"
    ]
    assert len(matching_payloads) == 1


def test_pre40_adjacent_upload_migration_consumes_mirrored_binding() -> None:
    comparison_file = _file_entry("adjacent.tsv", b"query\tsubject\t99\n")

    migrated_config, migrated_files = (
        session_io_module.migrate_legacy_linear_comparison_draft_for_current_writer(
            {
                "blastSource": "upload",
                "linearRecordLayout": {
                    "enabled": False,
                    "rows": [
                        {"uid": "record-a", "row": 1},
                        {"uid": "record-b", "row": 2},
                    ],
                },
            },
            {
                "linearSeqs": [
                    {"uid": "record-a", "blast": comparison_file},
                    {"uid": "record-b"},
                ],
                "linearComparisons": [
                    {
                        "id": "mirrored-a-b",
                        "queryUid": "record-a",
                        "subjectUid": "record-b",
                        "source": "upload",
                        "file": comparison_file,
                    }
                ],
            },
            force_web_draft=True,
        )
    )

    assert migrated_config["linearComparisonPlan"]["edges"] == [
        {
            "id": "linear-comparison-migrated-adjacent-1-record-a-record-b",
            "queryUid": "record-a",
            "subjectUid": "record-b",
            "included": True,
            "fileActive": True,
            "losatFilenameActive": False,
            "source": "upload",
            "losatFilename": "",
        }
    ]
    assert migrated_files["linearComparisons"] == [
        {
            "id": "linear-comparison-migrated-adjacent-1-record-a-record-b",
            "file": comparison_file,
        }
    ]


def test_current_comparison_authority_rejects_retired_version40_shape() -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        canonical_request=_canonical_request("linear"),
    )
    payload["config"]["blastSource"] = "losat"
    with pytest.raises(ValidationError, match="retired blastSource"):
        validate_session(payload)

    payload["config"].pop("blastSource")
    payload["config"]["linearRecordLayout"] = {
        "enabled": False,
        "recordGap": 24,
        "rows": [],
        "comparisons": [],
    }
    with pytest.raises(ValidationError, match="linearRecordLayout.comparisons"):
        validate_session(payload)


def test_current_comparison_authority_validates_file_bindings() -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        canonical_request=_canonical_request("linear"),
    )
    resource_id = next(iter(payload["resources"]))
    payload["config"]["linearComparisonPlan"] = {
        "mode": "selected",
        "defaultSource": "losat",
        "edges": [
            {
                "id": "edge-a-b",
                "queryUid": "a",
                "subjectUid": "b",
                "included": True,
                "fileActive": True,
                "losatFilenameActive": False,
                "source": "upload",
                "losatFilename": "",
            }
        ],
    }
    payload["webFiles"] = {
        "bindings": {
            "schema": 1,
            "linearComparisons": [
                {"id": "edge-a-b", "file": {"resourceId": resource_id}}
            ],
        }
    }
    validate_session(payload)

    unsupported_schema = copy.deepcopy(payload)
    unsupported_schema["webFiles"]["bindings"]["schema"] = 2
    with pytest.raises(ValidationError, match="Unsupported Web file binding schema"):
        validate_session(unsupported_schema)

    malformed_file = copy.deepcopy(payload)
    malformed_file["webFiles"]["bindings"]["linearComparisons"][0]["file"] = None
    with pytest.raises(ValidationError, match="requires a file resource binding"):
        validate_session(malformed_file)

    missing_resource = copy.deepcopy(payload)
    missing_resource["webFiles"]["bindings"]["linearComparisons"][0]["file"] = {
        "resourceId": "missing-resource"
    }
    with pytest.raises(ValidationError, match="references a missing resource"):
        validate_session(missing_resource)

    missing_active_binding = copy.deepcopy(payload)
    missing_active_binding["webFiles"]["bindings"]["linearComparisons"] = []
    with pytest.raises(
        ValidationError,
        match="Active comparison file is missing its Web file binding",
    ):
        validate_session(missing_active_binding)


def test_current_writer_quarantines_v27_to_v33_protein_artifacts(
    tmp_path: Path,
) -> None:
    legacy_raw = _legacy_protein_cache_entry()
    legacy_derived = {
        "schema": 1,
        "kind": "derived-losatp-payload",
        "key": "legacy-derived-key",
        "mode": "orthogroup",
        "payload": {"groups": []},
    }
    source = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="old",
            render_formats=("svg",),
        ),
        svg_results=(("old", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 20),
        canonical_request=_canonical_request("linear"),
    )
    source["version"] = 33
    source["losatCache"] = {"entries": [legacy_raw]}
    source["losatDerivedCache"] = {"entries": [legacy_derived]}
    source.pop("proteinIdentityManifest")

    promoted = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="new",
            render_formats=("svg",),
            source_session=source,
        ),
        svg_results=(("new", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        canonical_request=_canonical_request("linear"),
    )

    assert promoted["version"] == CURRENT_SESSION_VERSION
    assert promoted["losatCache"]["entries"] == []
    assert promoted["losatDerivedCache"]["entries"] == []
    assert promoted["legacyArtifacts"]["proteinRawCandidates"]["entries"] == [
        {
            "state": "pending",
            "originalEntry": legacy_raw,
            "rejectionReason": None,
        }
    ]
    assert promoted["legacyArtifacts"]["proteinDerivedEvidence"]["entries"] == [
        legacy_derived
    ]
    assert promoted["proteinIdentityManifest"] == {
        "schema": 2,
        "proteinSets": {},
        "recordAnalyses": {},
        "recordInstances": {},
    }

    session_path = tmp_path / "promoted.gbdraw-session.json.gz"
    write_session_json(session_path, promoted)
    assert load_session(session_path) == promoted


def test_future_session_version_fails() -> None:
    session = _minimal_session({})
    session["version"] = CURRENT_SESSION_VERSION + 1

    with pytest.raises(ValidationError, match="newer"):
        validate_session(session)


def test_current_session_rejects_legacy_files_but_version_39_accepts_them() -> None:
    session = build_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 7, 31),
        canonical_request=_canonical_request("circular"),
    )
    session["files"] = {}

    with pytest.raises(ValidationError, match="cannot contain legacy files"):
        validate_session(session)

    session["version"] = 39
    validate_session(session)


@pytest.mark.parametrize("version", (34, 35, 36, 37, 38))
def test_branch_internal_session_versions_are_rejected_at_read_and_rewrite_boundaries(
    version: int,
) -> None:
    session = _minimal_session({})
    session["version"] = version

    with pytest.raises(ValidationError, match=f"Unsupported session version: {version}"):
        validate_session(session)
    assert version not in SUPPORTED_SESSION_VERSIONS

    with pytest.raises(ValidationError, match=f"Unsupported session version: {version}"):
        build_session_json(
            SessionBuildContext(
                mode="circular",
                output_prefix="out",
                render_formats=("svg",),
                source_session=session,
            ),
            svg_results=(("out", "<svg></svg>"),),
            embedded_files={},
            generated_at=datetime(2026, 7, 28),
            canonical_request=_canonical_request("circular"),
        )


@pytest.mark.parametrize("schema", (3, 4))
@pytest.mark.parametrize("session_version", (33, CURRENT_SESSION_VERSION))
def test_branch_internal_request_schemas_are_rejected_at_session_boundaries(
    tmp_path: Path,
    schema: int,
    session_version: int,
) -> None:
    session = build_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 7, 30),
        canonical_request=_canonical_request("circular"),
    )
    session["version"] = session_version
    session["renderRequest"]["schema"] = schema
    message = f"Unsupported canonical renderRequest schema: {schema}"

    with pytest.raises(ValidationError, match=message):
        validate_session(session)
    with pytest.raises(SessionFormatError, match=message):
        load_session_document(session)

    input_path = (
        tmp_path
        / f"version-{session_version}-schema-{schema}-input.gbdraw-session.json"
    )
    input_path.write_text(json.dumps(session), encoding="utf-8")
    with pytest.raises(ValidationError, match=message):
        load_session(input_path)

    output_path = (
        tmp_path
        / f"version-{session_version}-schema-{schema}-output.gbdraw-session.json"
    )
    with pytest.raises(ValidationError, match=message):
        write_session_json(output_path, session)
    assert not output_path.exists()


def test_session_version_27_remains_supported() -> None:
    session = _minimal_session({})
    session["version"] = 27

    validate_session(session)
    assert 27 in SUPPORTED_SESSION_VERSIONS


def test_legacy_cli_repeat_shape_migration_is_selective_and_idempotent() -> None:
    args = ["--gbk", "input.gb", "-f", "svg"]
    migrated = migrate_legacy_repeat_feature_shape_args(args, session_version=30)
    assert migrated == [
        "--gbk",
        "input.gb",
        "--feature_shape",
        "repeat_region=rectangle",
        "-f",
        "svg",
    ]
    assert migrate_legacy_repeat_feature_shape_args(
        migrated, session_version=30
    ) == migrated
    assert migrate_legacy_repeat_feature_shape_args(
        ["--features", "CDS,tRNA"], session_version=30
    ) == ["--features", "CDS,tRNA"]
    assert migrate_legacy_repeat_feature_shape_args(
        ["--feature_shape", "repeat_region=underlay"], session_version=30
    ) == ["--feature_shape", "repeat_region=underlay"]
    assert migrate_legacy_repeat_feature_shape_args(args, session_version=31) == args


def test_legacy_gui_config_migrates_repeat_to_rectangle_without_mutating_source(
    tmp_path: Path,
) -> None:
    session = _minimal_session({"c_gb": _file_entry("input.gb", b"LOCUS       A\n")})
    session["config"]["adv"]["features"] = ["CDS", "repeat_region"]

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assignment_index = spec.args.index("--feature_shape")
    assert spec.args[assignment_index + 1] == "repeat_region=rectangle"
    assert "feature_shapes" not in session["config"]["adv"]


def test_gui_session_replays_arrowhead_shape_and_geometry_flags(tmp_path: Path) -> None:
    session = _minimal_session({"c_gb": _file_entry("input.gb", b"LOCUS       A\n")})
    session["config"]["adv"].update(
        {
            "feature_shapes": {"CDS": "arrowhead"},
            "arrow_head_length_ratio": "1.25",
            "arrowhead_shaft_width_ratio": "0.25",
        }
    )

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "CDS=arrowhead" in session_io_module._option_all_values(
        spec.args,
        "--feature_shape",
    )
    head_index = spec.args.index("--arrow_head_length_ratio")
    assert spec.args[head_index + 1] == "1.25"
    shaft_index = spec.args.index("--arrowhead_shaft_width_ratio")
    assert spec.args[shaft_index + 1] == "0.25"


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    ("head_ratio", "shaft_ratio"),
    [("auto", "0.5"), ("1.25", "0.25")],
)
def test_cli_session_projects_arrowhead_shape_and_geometry_to_web_state(
    mode: str,
    head_ratio: str,
    shaft_ratio: str,
) -> None:
    payload = build_session_json(
        SessionBuildContext(
            mode=mode,
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=(
                "--gbk",
                "input.gb",
                "--feature_shape",
                "CDS=arrowhead",
                "--arrow_head_length_ratio",
                head_ratio,
                "--arrowhead_shaft_width_ratio",
                shaft_ratio,
            ),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 8, 1),
        canonical_request=_canonical_request(mode),
    )

    adv = payload["config"]["adv"]
    assert adv["feature_shapes"] == {"CDS": "arrowhead"}
    assert adv["arrow_head_length_ratio"] == head_ratio
    assert adv["arrowhead_shaft_width_ratio"] == shaft_ratio


@pytest.mark.parametrize(
    ("version", "request_schema"),
    (
        (31, 1),
        (32, 2),
        (33, 2),
        (39, CANONICAL_REQUEST_SCHEMA),
        (CURRENT_SESSION_VERSION, CANONICAL_REQUEST_SCHEMA),
    ),
)
def test_main_backed_and_current_canonical_session_schemas_remain_supported(
    version: int,
    request_schema: int,
) -> None:
    session = {
        "format": SESSION_FORMAT,
        "version": version,
        "renderRequest": {"schema": request_schema},
        "resources": {},
    }
    if version == CURRENT_SESSION_VERSION:
        session["results"] = []
        session["editorState"] = {"featureCatalog": None}

    validate_session(session)
    assert version in SUPPORTED_SESSION_VERSIONS


def test_session_json_gzip_round_trip(tmp_path: Path) -> None:
    session = _minimal_session({})
    session_path = tmp_path / "session.gbdraw-session.json.gz"

    write_session_json(session_path, session)

    assert session_path.read_bytes().startswith(b"\x1f\x8b")
    assert load_session(session_path) == session


def test_embedded_file_materialization_sanitizes_name(tmp_path: Path) -> None:
    path = materialize_embedded_file(
        _file_entry("../unsafe/input.gb", b"LOCUS       TEST\n"),
        temp_dir=tmp_path,
        role="gbk",
    )

    assert path.parent == tmp_path.resolve()
    assert ".." not in path.name
    assert path.read_bytes() == b"LOCUS       TEST\n"


def test_embedded_file_size_mismatch_fails(tmp_path: Path) -> None:
    entry = _file_entry("input.gb", b"abc")
    entry["size"] = 4

    with pytest.raises(ValidationError, match="size mismatch"):
        materialize_embedded_file(entry, temp_dir=tmp_path, role="gbk")


def test_depth_payload_decode_matches_browser_codec_shape() -> None:
    payload = {
        "schema": 1,
        "columns": ["reference_name", "position", "depth"],
        "lineEnding": "\n",
        "finalNewline": True,
        "rowCount": 3,
        "header": ["reference_name", "position", "depth"],
        "records": [{"id": "seq1", "runs": [[1, 1, 3, ["10", "12", "0"]]]}],
    }

    assert decode_depth_payload(payload) == (
        "reference_name\tposition\tdepth\n"
        "seq1\t1\t10\n"
        "seq1\t2\t12\n"
        "seq1\t3\t0\n"
    )


def test_cli_invocation_restoration_substitutes_embedded_files(tmp_path: Path) -> None:
    session = _minimal_session(
        {"c_gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
    )
    session["cliInvocation"] = {
        "schema": 1,
        "mode": "circular",
        "args": ["-o", "old", "--gbk", "input.gb", "-f", "svg"],
        "renderFormats": ["svg"],
        "fileBindings": [{"argIndex": 3, "slot": "files.c_gb", "name": "input.gb"}],
        "generatedBy": "gbdraw",
    }

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override="new",
        format_override="interactive-svg",
    )

    assert spec.args[2] == "--gbk"
    materialized = Path(spec.args[3])
    assert materialized.exists()
    assert materialized.read_bytes() == b"LOCUS       TEST\n"
    assert spec.args[0:2] == ("-o", "new")
    assert spec.args[-2:] == ("-f", "interactive_svg")
    assert spec.cli_invocation_args[0:2] == ("-o", "new")
    assert spec.cli_invocation_args[-2:] == ("-f", "interactive_svg")
    assert spec.file_bindings[0].argIndex == 3
    assert spec.file_bindings[0].slot == "files.c_gb"


def test_circular_legacy_cli_invocation_is_canonicalized_after_materialization(
    tmp_path: Path,
) -> None:
    session = _minimal_session(
        {
            "c_gb": _file_entry("input.gb", b"LOCUS       TEST\n"),
            "c_depth": [_file_entry("depth.tsv", b"TEST\t1\t10\n")],
        }
    )
    session["cliInvocation"] = {
        "schema": 1,
        "mode": "circular",
        "args": [
            "--suppress_gc",
            "--suppress_skew",
            "--show_depth",
            "--gbk",
            "input.gb",
            "--depth",
            "depth.tsv",
            "--depth_tick_interval",
            "4",
            "--feature_table",
            "features.tsv",
            "--multi_record_size_mode",
            "sqrt",
            "--circular-track-slot",
            "gc_content:dinucleotide_content@spacing=4px,strict=true,compress=true,reserve=true",
            "-f",
            "svg",
        ],
        "renderFormats": ["svg"],
        "fileBindings": [
            {"argIndex": 4, "slot": "files.c_gb", "name": "input.gb"},
            {"argIndex": 6, "slot": "files.c_depth[0]", "name": "depth.tsv"},
        ],
        "generatedBy": "gbdraw",
    }

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--no-gc" in spec.args
    assert "--no-skew" in spec.args
    assert "--depth_track" in spec.args
    assert spec.args[spec.args.index("--depth_large_tick_interval") + 1] == "4"
    assert spec.args[spec.args.index("--feature_visibility_table") + 1] == "features.tsv"
    assert spec.args[spec.args.index("--multi_record_size_mode") + 1] == "auto"
    slot_value = spec.args[spec.args.index("--circular_track_slot") + 1]
    assert "__gbdraw_legacy_spacing=4px" in slot_value
    assert not any(
        token in slot_value
        for token in ("@spacing=", ",spacing=", "strict=", "compress=", "reserve=")
    )
    assert not {
        "--suppress_gc",
        "--suppress_skew",
        "--show_depth",
        "--depth",
        "--depth_tick_interval",
        "--feature_table",
        "--circular-track-slot",
    }.intersection(spec.args)
    assert [binding.argIndex for binding in spec.file_bindings] == [3, 5]
    parsed = circular_cli_module._get_args(
        list(spec.args),
        _allow_legacy_track_transport=True,
    )
    assert parsed.show_gc is False
    assert parsed.show_skew is False
    assert parsed.depth_track == [[spec.args[5]]]
    assert parsed.depth_large_tick_interval == pytest.approx(4)
    assert parsed.feature_table == "features.tsv"
    assert parsed.multi_record_size_mode == "auto"
    assert not hasattr(parsed, "depth")
    assert not hasattr(parsed, "show_depth")


@pytest.mark.parametrize(
    ("legacy_layout", "canonical_layout"),
    [("spreadout", "above"), ("tuckin", "below")],
)
def test_linear_legacy_cli_invocation_is_canonicalized_after_materialization(
    tmp_path: Path,
    legacy_layout: str,
    canonical_layout: str,
) -> None:
    session = _minimal_session(
        {
            "linearSeqs": [
                {"gb": _file_entry("a.gb", b"LOCUS       A\n")},
                {"gb": _file_entry("b.gb", b"LOCUS       B\n")},
            ],
            "linearDepth": [
                _file_entry("a.tsv", b"A\t1\t10\n"),
                _file_entry("b.tsv", b"B\t1\t12\n"),
            ],
        },
        mode="linear",
    )
    session["cliInvocation"] = {
        "schema": 1,
        "mode": "linear",
        "args": [
            "--show_gc",
            "--show_skew",
            "--show_depth",
            "--gbk",
            "a.gb",
            "b.gb",
            "--depth",
            "a.tsv",
            "b.tsv",
            "--depth_tick_interval",
            "4",
            "--feature_table",
            "features.tsv",
            "--collinear_max_gene_gap",
            "3",
            "--pairwise-match-style",
            "curve",
            "--label_placement",
            "on_feature",
            "--track_layout",
            legacy_layout,
            "-f",
            "svg",
        ],
        "renderFormats": ["svg"],
        "fileBindings": [
            {"argIndex": 4, "slot": "files.linearSeqs[0].gb", "name": "a.gb"},
            {"argIndex": 5, "slot": "files.linearSeqs[1].gb", "name": "b.gb"},
            {"argIndex": 7, "slot": "files.linearDepth[0]", "name": "a.tsv"},
            {"argIndex": 8, "slot": "files.linearDepth[1]", "name": "b.tsv"},
        ],
        "generatedBy": "gbdraw",
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--gc" in spec.args
    assert "--skew" in spec.args
    assert "--depth_track" in spec.args
    assert spec.args[spec.args.index("--depth_large_tick_interval") + 1] == "4"
    assert spec.args[spec.args.index("--feature_visibility_table") + 1] == "features.tsv"
    assert spec.args[spec.args.index("--collinear_max_unit_gap") + 1] == "3"
    assert spec.args[spec.args.index("--pairwise_match_style") + 1] == "curve"
    assert spec.args[spec.args.index("--label_placement") + 1] == "above_feature"
    assert spec.args[spec.args.index("--track_layout") + 1] == canonical_layout
    assert not {
        "--show_gc",
        "--show_skew",
        "--show_depth",
        "--depth",
        "--depth_tick_interval",
        "--feature_table",
        "--collinear_max_gene_gap",
        "--pairwise-match-style",
    }.intersection(spec.args)
    assert [binding.argIndex for binding in spec.file_bindings] == [3, 4, 6, 7]
    parsed = get_linear_args(list(spec.args))
    assert parsed.show_gc is True
    assert parsed.show_skew is True
    assert parsed.depth_track == [[spec.args[6], spec.args[7]]]
    assert parsed.depth_large_tick_interval == pytest.approx(4)
    assert parsed.feature_table == "features.tsv"
    assert parsed.collinear_max_unit_gap == 3
    assert parsed.pairwise_match_style == "curve"
    assert parsed.label_placement == "above_feature"
    assert parsed.track_layout == canonical_layout
    assert not hasattr(parsed, "depth")
    assert not hasattr(parsed, "show_depth")


def test_cli_session_capture_embeds_records_table_dependencies(tmp_path: Path) -> None:
    gbk_path = tmp_path / "input.gb"
    gbk_path.write_text("LOCUS       TEST\n", encoding="utf-8")
    table_path = tmp_path / "records.tsv"
    table_path.write_text(
        "gbk\trecord_id\n"
        "input.gb\t#1\n",
        encoding="utf-8",
    )

    files, bindings = collect_embedded_files_from_cli_args(
        "linear",
        ["--records_table", str(table_path), "-f", "svg"],
    )

    assert bindings[0].argIndex == 1
    assert bindings[0].slot == "files.cliInputs[0]"
    assert files["cliInputs"][0]["name"] == "records.tsv"
    assert files["cliInputs"][1]["name"] == "input.gb"
    assert files["cliTables"] == [
        {
            "argIndex": 1,
            "kind": "records",
            "slot": "files.cliInputs[0]",
            "dependencies": [
                {
                    "rowIndex": 0,
                    "rowNumber": 2,
                    "column": "gbk",
                    "slot": "files.cliInputs[1]",
                }
            ],
        }
    ]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_session_capture_embeds_only_canonical_depth_track_option(
    mode: str,
    tmp_path: Path,
) -> None:
    depth_path = tmp_path / "depth.tsv"
    depth_path.write_text("record\t1\t10\n", encoding="utf-8")

    canonical_files, canonical_bindings = collect_embedded_files_from_cli_args(
        mode,
        ["--depth_track", str(depth_path)],
    )
    _legacy_files, legacy_bindings = collect_embedded_files_from_cli_args(
        mode,
        ["--depth", str(depth_path)],
    )

    assert len(canonical_bindings) == 1
    assert canonical_bindings[0].argIndex == 1
    assert legacy_bindings == ()
    if mode == "circular":
        assert canonical_bindings[0].slot == "files.c_depth[0]"
        assert canonical_files["c_depth"][0]["name"] == "depth.tsv"
    else:
        assert canonical_bindings[0].slot == "files.linearSeqs[0].depth[0]"
        assert canonical_files["linearSeqs"][0]["depth"][0]["name"] == "depth.tsv"


def test_cli_invocation_restoration_rewrites_records_table_paths(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "cliInputs": [
                _file_entry(
                    "records.tsv",
                    b"gbk\trecord_id\noriginal.gb\t#1\n",
                ),
                _file_entry("original.gb", b"LOCUS       TEST\n"),
            ],
            "cliTables": [
                {
                    "argIndex": 1,
                    "kind": "records",
                    "slot": "files.cliInputs[0]",
                    "dependencies": [
                        {
                            "rowIndex": 0,
                            "rowNumber": 2,
                            "column": "gbk",
                            "slot": "files.cliInputs[1]",
                        }
                    ],
                }
            ],
        },
        mode="linear",
    )
    session["cliInvocation"] = {
        "schema": 1,
        "mode": "linear",
        "args": ["--records_table", "records.tsv", "-f", "svg"],
        "renderFormats": ["svg"],
        "fileBindings": [
            {"argIndex": 1, "slot": "files.cliInputs[0]", "name": "records.tsv"}
        ],
        "generatedBy": "gbdraw",
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    table_path = Path(spec.args[1])
    assert table_path.exists()
    table_lines = table_path.read_text(encoding="utf-8").splitlines()
    assert table_lines[0] == "gbk\trecord_id"
    restored_gbk_name = table_lines[1].split("\t")[0]
    restored_gbk_path = table_path.parent / restored_gbk_name
    assert restored_gbk_path.read_bytes() == b"LOCUS       TEST\n"


def test_cli_session_round_trip_rewrites_bom_records_table_first_column(
    tmp_path: Path,
) -> None:
    gbk_path = tmp_path / "original.gb"
    gbk_path.write_text("LOCUS       TEST\n", encoding="utf-8")
    table_path = tmp_path / "records.tsv"
    table_path.write_text("gbk\trecord_id\noriginal.gb\t#1\n", encoding="utf-8-sig")

    files, bindings = collect_embedded_files_from_cli_args(
        "linear",
        ["--records_table", str(table_path), "-f", "svg"],
    )
    session = _minimal_session(files, mode="linear")
    session["cliInvocation"] = {
        "schema": 1,
        "mode": "linear",
        "args": ["--records_table", "records.tsv", "-f", "svg"],
        "renderFormats": ["svg"],
        "fileBindings": [bindings[0].__dict__],
        "generatedBy": "gbdraw",
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path / "restored",
        output_override=None,
        format_override=None,
    )

    restored_table = Path(spec.args[1])
    restored_text = restored_table.read_text(encoding="utf-8")
    assert restored_text.startswith("gbk\trecord_id\n")
    assert not restored_text.startswith("\ufeff")
    restored_gbk_name = restored_text.splitlines()[1].split("\t")[0]
    restored_gbk_path = restored_table.parent / restored_gbk_name
    assert restored_gbk_path.read_bytes() == b"LOCUS       TEST\n"
    assert read_records_table(str(restored_table)).gbk_files == [str(restored_gbk_path.resolve())]


def test_gui_only_circular_session_maps_to_cli_args(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "c_gb": _file_entry("input.gb", b"LOCUS       TEST\n"),
            "c_depth": _file_entry("depth.tsv", b"TEST\t1\t10\n"),
        }
    )
    session["config"]["form"].update(
        {
            "suppress_gc": True,
            "suppress_skew": True,
            "show_depth": True,
        }
    )
    session["config"]["adv"]["multi_record_size_mode"] = "sqrt"

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--gbk" in spec.args
    gbk_index = spec.args.index("--gbk") + 1
    assert Path(spec.args[gbk_index]).read_bytes() == b"LOCUS       TEST\n"
    assert spec.cli_invocation_args[gbk_index] == "input.gb"
    assert spec.file_bindings[0].argIndex == gbk_index
    assert "--no-gc" in spec.args
    assert "--no-skew" in spec.args
    assert "--depth_track" in spec.args
    assert spec.args[spec.args.index("--multi_record_size_mode") + 1] == "auto"
    assert "--show_depth" not in spec.args
    assert "--suppress_gc" not in spec.args
    assert "--suppress_skew" not in spec.args


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_gui_only_session_restores_hidden_scale_as_cli_flag(
    tmp_path: Path,
    mode: str,
) -> None:
    files = (
        {"c_gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
        if mode == "circular"
        else {
            "linearSeqs": [
                {"gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
            ]
        }
    )
    session = _minimal_session(files, mode=mode)
    session["config"]["form"]["show_scale"] = False
    if mode == "linear":
        session["config"]["form"]["linear_ruler_on_axis"] = True

    spec = session_to_cli_args(
        session,
        mode=mode,
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--hide_scale" in spec.args
    if mode == "linear":
        assert "--ruler_on_axis" in spec.args


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_sparse_gui_session_keeps_scale_visible_by_default(
    tmp_path: Path,
    mode: str,
) -> None:
    files = (
        {"c_gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
        if mode == "circular"
        else {
            "linearSeqs": [
                {"gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
            ]
        }
    )
    session = _minimal_session(files, mode=mode)

    spec = session_to_cli_args(
        session,
        mode=mode,
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--hide_scale" not in spec.args


def test_gui_only_linear_session_restores_losatp_blastp_args(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "linearSeqs": [
                {"gb": _file_entry("a.gb", b"LOCUS       A\n")},
                {"gb": _file_entry("b.gb", b"LOCUS       B\n")},
            ]
        },
        mode="linear",
    )
    session["config"]["form"]["show_gc"] = True
    session["config"]["form"]["show_skew"] = True
    session["config"]["adv"] = {
        "blastSource": "losat",
        "losatProgram": "blastp",
        "min_bitscore": 25,
        "evalue": "1e-4",
        "identity": 55,
        "alignment_length": 30,
        "label_placement": "on_feature",
    }
    session["config"]["form"]["linear_track_layout"] = "tuckin"
    session["config"]["form"]["show_labels_linear"] = "orthogroup_top"
    session["config"]["losat"] = {
        "executionMode": "auto",
        "threadsPerJob": "2",
        "blastp": {
            "mode": "collinear",
            "maxHits": 7,
            "collinearMinAnchors": 2,
            "collinearMaxGeneGap": 3,
            "collinearMaxDiagonalDrift": 4,
            "collinearMaxConflictsInMergeGap": 6,
            "collinearUnitMode": "cds",
            "collinearSearchScope": "all",
            "collinearColorMode": "orientation_identity",
            "collinearMaxParalogLinksPerOrthogroup": 5,
        },
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--protein_blastp_mode" in spec.args
    assert spec.args[spec.args.index("--protein_blastp_mode") + 1] == "collinear"
    assert spec.args[spec.args.index("--losatp_threads") + 1] == "2"
    assert spec.args[spec.args.index("--collinear_max_unit_gap") + 1] == "3"
    assert spec.args[spec.args.index("--collinear_max_diagonal_drift") + 1] == "4"
    assert spec.args[spec.args.index("--collinear_max_conflicts_in_merge_gap") + 1] == "6"
    assert spec.args[spec.args.index("--collinear_search_scope") + 1] == "all"
    assert spec.args[spec.args.index("--collinear_color_mode") + 1] == "orientation_identity"
    assert spec.args[spec.args.index("--label_placement") + 1] == "above_feature"
    assert spec.args[spec.args.index("--track_layout") + 1] == "below"
    assert "--collinear_max_gene_gap" not in spec.args
    assert spec.args[spec.args.index("--show_labels") + 1] == "orthogroup_top"
    assert "--gc" in spec.args
    assert "--skew" in spec.args
    assert "--show_gc" not in spec.args
    assert "--show_skew" not in spec.args


def test_linear_cli_accepts_orthogroup_top_show_labels() -> None:
    args = get_linear_args(["--gbk", "a.gb", "b.gb", "--show_labels", "orthogroup_top"])

    assert args.show_labels == "orthogroup_top"


def test_gui_only_linear_session_restores_row_options(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "linearSeqs": [
                {
                    "gb": _file_entry("a.gb", b"LOCUS       A\n"),
                    "definition": "Alpha",
                    "record_subtitle": "Alpha subtitle",
                    "region_record_id": "RecA",
                    "region_reverse": True,
                },
                {
                    "gb": _file_entry("b.gb", b"LOCUS       B\n"),
                    "definition": "Beta",
                    "record_subtitle": "Beta subtitle",
                    "region_record_id": "RecB",
                    "region_start": 10,
                    "region_end": 20,
                    "region_reverse": True,
                },
            ]
        },
        mode="linear",
    )

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert spec.args.count("--record_label") == 2
    first_label = spec.args.index("--record_label")
    assert spec.args[first_label + 1] == "Alpha"
    second_label = spec.args.index("--record_label", first_label + 2)
    assert spec.args[second_label + 1] == "Beta"
    assert spec.args.count("--record_subtitle") == 2
    first_subtitle = spec.args.index("--record_subtitle")
    assert spec.args[first_subtitle + 1] == "Alpha subtitle"
    second_subtitle = spec.args.index("--record_subtitle", first_subtitle + 2)
    assert spec.args[second_subtitle + 1] == "Beta subtitle"
    assert spec.args.count("--record_id") == 2
    first_selector = spec.args.index("--record_id")
    assert spec.args[first_selector + 1] == "RecA"
    second_selector = spec.args.index("--record_id", first_selector + 2)
    assert spec.args[second_selector + 1] == "RecB"
    assert spec.args.count("--reverse_complement") == 2
    first_reverse = spec.args.index("--reverse_complement")
    assert spec.args[first_reverse + 1] == "1"
    second_reverse = spec.args.index("--reverse_complement", first_reverse + 2)
    assert spec.args[second_reverse + 1] == "0"
    assert spec.args[spec.args.index("--region") + 1] == "#2:10-20:rc"


def test_gui_only_linear_session_restores_orthogroup_alignment_target(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "linearSeqs": [
                {"gb": _file_entry("a.gb", b"LOCUS       A\n")},
                {"gb": _file_entry("b.gb", b"LOCUS       B\n")},
            ]
        },
        mode="linear",
    )
    session["config"]["blastSource"] = "losat"
    session["config"]["losatProgram"] = "blastp"
    session["config"]["losat"] = {
        "threadsPerJob": "auto",
        "blastp": {"mode": "orthogroup"},
    }
    session["orthogroupState"] = {
        "groups": [{"id": "og1"}],
        "selectedOrthogroupAlignmentFeature": "target_feature",
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert spec.args[spec.args.index("--protein_blastp_mode") + 1] == "orthogroup"
    assert spec.args[spec.args.index("--align_orthogroup_feature") + 1] == "target_feature"


def test_gui_only_linear_session_restores_top_level_losatp_keys(tmp_path: Path) -> None:
    session = _minimal_session(
        {
            "linearSeqs": [
                {"gb": _file_entry("a.gb", b"LOCUS       A\n")},
                {"gb": _file_entry("b.gb", b"LOCUS       B\n")},
            ]
        },
        mode="linear",
    )
    session["config"]["blastSource"] = "losat"
    session["config"]["losatProgram"] = "blastp"
    session["config"]["losat"] = {
        "threadsPerJob": "4",
        "blastp": {"mode": "orthogroup"},
    }

    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "--protein_blastp_mode" in spec.args
    assert spec.args[spec.args.index("--protein_blastp_mode") + 1] == "orthogroup"
    assert spec.args[spec.args.index("--losatp_threads") + 1] == "4"


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_session_projects_hidden_scale_to_web_and_canonical_state(
    mode: str,
) -> None:
    args = ["--gbk", "input.gb", "--hide_scale"]
    if mode == "linear":
        args.extend(
            [
                "--track_layout",
                "above",
                "--scale_style",
                "ruler",
                "--ruler_on_axis",
            ]
        )
    payload = build_session_json(
        SessionBuildContext(
            mode=mode,
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=tuple(args),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 7, 31),
        canonical_request=_canonical_request_with_scale_visibility(
            mode,
            show=False,
        ),
    )

    assert payload["config"]["form"]["show_scale"] is False
    assert (
        payload["renderRequest"]["diagramOptions"]["config"]["objects"]["scale"][
            "show"
        ]
        is False
    )
    if mode == "linear":
        assert payload["config"]["form"]["linear_ruler_on_axis"] is False


def test_cli_session_config_includes_lossless_cli_options() -> None:
    args = (
        "-f",
        "interactive_svg",
        "--gbk",
        "AP027078.gb",
        "AP027131.gb",
        "AP027133.gb",
        "AP027132.gb",
        "NZ_CP006932.gb",
        "--protein_blastp_mode",
        "orthogroup",
        "--losatp_threads",
        "32",
        "--align_center",
        "--separate_strands",
        "--pairwise_match_style",
        "curve",
        "--depth_large_tick_interval",
        "4",
        "--feature_visibility_table",
        "features.tsv",
        "--collinear_max_unit_gap",
        "3",
        "--label_placement",
        "above_feature",
        "--track_layout",
        "above",
        "--scale_style",
        "ruler",
        "--palette",
        "ajisai",
        "--gc",
        "--skew",
        "--show_labels",
        "orthogroup_top",
    )

    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("interactive-svg",),
            cli_invocation_args=args,
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 6, 23),
        canonical_request=_canonical_request("linear"),
    )

    config = payload["config"]
    assert payload["cliInvocation"]["renderFormats"] == ["interactive_svg"]
    stored_args = payload["cliInvocation"]["args"]
    assert not {
        "--depth_tick_interval",
        "--feature_table",
        "--collinear_max_gene_gap",
        "--pairwise-match-style",
        "--save-session",
        "--session-output",
        "on_feature",
        "spreadout",
        "tuckin",
        "sqrt",
    }.intersection(stored_args)
    assert config["form"]["prefix"] == "out"
    assert config["form"]["align_center"] is True
    assert config["form"]["separate_strands"] is True
    assert config["form"]["show_scale"] is True
    assert config["form"]["scale_style"] == "ruler"
    assert config["form"]["show_gc"] is True
    assert config["form"]["show_skew"] is True
    assert config["form"]["show_labels_linear"] == "orthogroup_top"
    assert config["adv"]["pairwise_match_style"] == "curve"
    assert config["adv"]["depth_large_tick_interval"] == "4"
    assert "depth_tick_interval" not in config["adv"]
    assert config["palette"] == "ajisai"
    assert "blastSource" not in config
    assert "blastSource" not in config["adv"]
    assert config["losatProgram"] == "blastp"
    assert config["losat"]["threadsPerJob"] == "32"
    assert config["losat"]["blastp"]["mode"] == "orthogroup"
    assert config["losat"]["blastp"]["collinearMaxUnitGap"] == "3"
    assert "collinearMaxGeneGap" not in config["losat"]["blastp"]
    assert "collinearBlockMergeGap" not in config["losat"]["blastp"]
    assert "collinearSingletonMergeGap" not in config["losat"]["blastp"]
    assert config["cliOptions"]["rawArgs"] == list(args)
    assert config["cliOptions"]["byKey"]["protein_blastp_mode"] == ["orthogroup"]
    assert config["cliOptions"]["byKey"]["losatp_threads"] == ["32"]
    assert config["cliOptions"]["byKey"]["palette"] == ["ajisai"]
    assert config["cliOptions"]["byKey"]["gbk"] == [
        ["AP027078.gb", "AP027131.gb", "AP027133.gb", "AP027132.gb", "NZ_CP006932.gb"]
    ]


def test_current_cli_session_writer_uses_canonical_linear_inventory() -> None:
    args = (
        "--gbk",
        "a.gb",
        "b.gb",
        "--record_id",
        "RecA",
        "--record_id",
        "RecB",
        "--reverse_complement",
        "true",
        "--reverse_complement",
        "false",
        "--record_label",
        "Alpha",
        "--record_label",
        "Beta",
        "--record_subtitle",
        "Alpha subtitle",
        "--record_subtitle",
        "Beta subtitle",
        "--region",
        "RecB:10-20:rc",
        "--protein_blastp_mode",
        "orthogroup",
        "--align_orthogroup_feature",
        "target_feature",
    )
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=args,
            linear_record_metadata=(
                {
                    "loaded_index": 0,
                    "source_index": 0,
                    "source_loaded_index": 0,
                    "source_loaded_count": 1,
                    "record_id": "RecA",
                    "source_file": "a.gb",
                    "source_basename": "a.gb",
                },
                {
                    "loaded_index": 1,
                    "source_index": 1,
                    "source_loaded_index": 0,
                    "source_loaded_count": 1,
                    "record_id": "RecB",
                    "source_file": "b.gb",
                    "source_basename": "b.gb",
                },
            ),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={
            "linearSeqs": [
                {"definition": "", "region_record_id": "", "region_start": None, "region_end": None, "region_reverse": False},
                {"definition": "", "region_record_id": "", "region_start": None, "region_end": None, "region_reverse": False},
            ]
        },
        generated_at=datetime(2026, 6, 23),
        canonical_request=_canonical_request("linear"),
    )

    assert payload["version"] == CURRENT_SESSION_VERSION
    assert "files" not in payload
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert payload["resources"]
    assert payload["orthogroupState"]["selectedOrthogroupAlignmentFeature"] == "target_feature"
    assert "linearRecordLayout" not in payload["config"]
    assert "linearComparisonPlan" not in payload["config"]


def test_current_cli_session_writer_omits_legacy_ambiguous_linear_files() -> None:
    args = (
        "--gbk",
        "multi.gb",
        "--record_label",
        "Alpha",
        "--record_label",
        "Beta",
        "--region",
        "#1:10-20",
        "--region",
        "#2:30-40",
    )
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=args,
            linear_record_metadata=(
                {
                    "loaded_index": 0,
                    "source_index": 0,
                    "source_loaded_index": 0,
                    "source_loaded_count": 2,
                    "record_id": "RecA",
                    "source_file": "multi.gb",
                    "source_basename": "multi.gb",
                },
                {
                    "loaded_index": 1,
                    "source_index": 0,
                    "source_loaded_index": 1,
                    "source_loaded_count": 2,
                    "record_id": "RecB",
                    "source_file": "multi.gb",
                    "source_basename": "multi.gb",
                },
            ),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={
            "linearSeqs": [
                {"definition": "", "region_record_id": "", "region_start": None, "region_end": None, "region_reverse": False}
            ]
        },
        generated_at=datetime(2026, 6, 23),
        canonical_request=_canonical_request("linear"),
    )

    assert payload["version"] == CURRENT_SESSION_VERSION
    assert "files" not in payload
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert payload["resources"]


def test_session_pre_parse_rejects_unsupported_options() -> None:
    with pytest.raises(SystemExit):
        parse_session_pre_args(
            ["--session", "session.gbdraw-session.json", "--gbk", "input.gb"],
            mode="circular",
        )


def test_session_cli_help_uses_underscore_options(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        circular_cli_module._get_args(["--help"])

    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    assert "--save_session" in help_text
    assert "--session_output" in help_text
    assert ".gz suffix" in help_text
    assert "gzip compression" in help_text
    assert "--save-session" not in help_text
    assert "--session-output" not in help_text


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    "option_args",
    [
        ["--save-session"],
        ["--session-output", "roundtrip.gbdraw-session.json"],
    ],
)
def test_session_preparser_rejects_removed_hyphen_aliases(
    mode: str,
    option_args: list[str],
) -> None:
    with pytest.raises(SystemExit):
        parse_session_pre_args(
            ["--session", "session.gbdraw-session.json", *option_args],
            mode=mode,
        )


@pytest.mark.parametrize(
    "get_args",
    [circular_cli_module._get_args, get_linear_args],
)
def test_current_cli_rejects_removed_feature_table(get_args) -> None:
    with pytest.raises(SystemExit):
        get_args(["--gbk", "input.gb", "--feature_table", "features.tsv"])


def test_current_session_output_controls_are_stripped_before_writing() -> None:
    assert strip_session_output_args(
        [
            "--gbk",
            "input.gb",
            "--save_session",
            "--session_output=roundtrip.gbdraw-session.json",
        ]
    ) == ["--gbk", "input.gb"]


def test_default_sidecar_path_resolution(tmp_path: Path) -> None:
    output = make_rendered_svg(str(tmp_path / "diagram"))

    assert resolve_session_sidecar_path(
        explicit_path=str(tmp_path / "custom.json"),
        output_prefix=None,
        outputs=[output],
    ) == tmp_path / "custom.json"
    assert resolve_session_sidecar_path(
        explicit_path=None,
        output_prefix=str(tmp_path / "prefix"),
        outputs=[output],
    ) == tmp_path / "prefix.gbdraw-session.json"
    assert resolve_session_sidecar_path(
        explicit_path=None,
        output_prefix=None,
        outputs=[output],
    ) == tmp_path / "diagram.gbdraw-session.json"


def test_render_formats_do_not_include_session_json_aliases() -> None:
    assert "json" not in ACCEPTED_FORMATS
    assert "session-json" not in ACCEPTED_FORMATS


def test_circular_cli_save_session_round_trip(tmp_path: Path, examples_dir: Path) -> None:
    output_prefix = tmp_path / "circular_session"
    session_path = tmp_path / "custom-session.gbdraw-session.json"
    circular_main(
        [
            "--gbk",
            str(examples_dir / "MjeNMV.gb"),
            "-o",
            str(output_prefix),
            "-f",
            "svg",
            "--save_session",
            "--session_output",
            str(session_path),
        ]
    )

    svg_path = output_prefix.with_suffix(".svg")
    assert svg_path.exists()
    assert session_path.exists()

    payload = load_session(session_path)
    assert payload["format"] == SESSION_FORMAT
    assert payload["version"] == CURRENT_SESSION_VERSION
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert (
        "outputPrefix"
        not in payload["renderRequest"]["diagramOptions"].get("output", {})
    )
    assert "files" not in payload
    assert payload["resources"]["record-1-genbank"]["data"]
    assert "<svg" in payload["results"][0]["content"]
    catalog = payload["editorState"]["featureCatalog"]
    assert catalog["schema"] == 3
    assert len(catalog["items"]) == len(payload["results"]) == 1
    assert catalog["items"][0]["features"]
    assert catalog["items"][0]["biologicalFeatures"]
    assert not {
        "extractedFeatures",
        "biologicalFeatures",
        "featureCatalog",
    }.intersection(payload["features"])
    assert "groups" not in payload["orthogroupState"]
    assert payload["cliInvocation"]["mode"] == "circular"
    assert payload["cliInvocation"]["fileBindings"][0]["slot"] == "files.c_gb"
    stored_args = payload["cliInvocation"]["args"]
    assert not {
        "--save_session",
        "--session_output",
        "--save-session",
        "--session-output",
    }.intersection(stored_args)

    regenerated_prefix = tmp_path / "regenerated"
    circular_main(
        [
            "--session",
            str(session_path),
            "-o",
            str(regenerated_prefix),
            "-f",
            "svg",
        ]
    )
    assert regenerated_prefix.with_suffix(".svg").exists()


def test_depth_session_entry_materializes_encoded_payload(tmp_path: Path) -> None:
    entry = {
        "name": "depth.tsv",
        "type": "text/tab-separated-values",
        "size": len("seq1\t1\t10\n".encode("utf-8")),
        "lastModified": 0,
        "encoding": DEPTH_FILE_ENCODING,
        "data": {
            "schema": 1,
            "columns": ["reference_name", "position", "depth"],
            "lineEnding": "\n",
            "finalNewline": True,
            "rowCount": 1,
            "header": None,
            "records": [{"id": "seq1", "runs": [[1, 1, 1, ["10"]]]}],
        },
    }

    path = materialize_embedded_file(entry, temp_dir=tmp_path, role="depth")
    assert path.read_text(encoding="utf-8") == "seq1\t1\t10\n"
