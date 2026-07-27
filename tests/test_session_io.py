from __future__ import annotations

import base64
import copy
import json
import re
from datetime import datetime
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw.circular as circular_cli_module
from gbdraw.analysis.protein_colinearity import (
    build_protein_losat_cache_key,
    build_protein_losat_pair_identity,
    extract_web_stable_cds_proteins,
    validate_protein_raw_entry_references,
)
from gbdraw.circular import circular_main
from gbdraw.linear import (
    _get_args as get_linear_args,
    _record_instance_keys_for_web_losat,
    _source_session_legacy_protein_candidates,
    _source_session_losat_entries,
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
    session_to_cli_args,
    validate_session,
    write_session_json,
)
from gbdraw.session import load_session_document


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
    query_transport = manifest["recordInstances"]["record-1"]["runtimeIds"][
        feature_a
    ]
    subject_transport = manifest["recordInstances"]["record-2"]["runtimeIds"][
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
            f"{query_transport}\t{subject_transport}\t"
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
    assert payload["features"]["extractedFeatures"][0]["orthogroup_id"] == "og_1"
    assert [
        feature["stable_feature_id"]
        for feature in payload["features"]["biologicalFeatures"]
    ] == ["feature-1", "feature-2"]
    group = payload["orthogroupState"]["groups"][0]
    assert group["id"] == "og_1"
    assert [member["proteinId"] for member in group["members"]] == [
        "protein-1",
        "protein-2",
    ]


def test_current_session_version_matches_web_config() -> None:
    source = Path("gbdraw/web/js/services/config.js").read_text(encoding="utf-8")
    match = re.search(r"const\s+SESSION_VERSION\s*=\s*(\d+);", source)
    assert match is not None
    assert CURRENT_SESSION_VERSION == 36
    assert SUPPORTED_SESSION_VERSIONS == frozenset(
        {27, 28, 29, 30, 31, 32, 33, CURRENT_SESSION_VERSION}
    )
    assert int(match.group(1)) == CURRENT_SESSION_VERSION


def test_current_session_feature_catalog_is_compact_and_lossless(
    tmp_path: Path,
) -> None:
    long_note = ("😀" * 30) + ("x" * 30)
    qualifiers = {
        "protein_id": ["protein-1"],
        "locus_tag": ["locus-1"],
        "product": ["example protein"],
        "translation": ["M"],
        "note": [long_note],
    }
    biological = {
        "id": "f0",
        "svg_id": "f_stable",
        "stable_svg_id": "f_stable",
        "stable_feature_id": "f_stable",
        "record_id": "record",
        "record_idx": 0,
        "feature_index": 0,
        "organism": "Example",
        "type": "CDS",
        "start": 1,
        "end": 3,
        "strand": "+",
        "protein_id": "protein-1",
        "source_protein_id": "protein-1",
        "locus_tag": "locus-1",
        "gene_id": "",
        "old_locus_tag": "",
        "gene": "",
        "product": "example protein",
        "note": long_note[:50],
        "qualifiers": qualifiers,
        "selector": {
            "qualifiers": copy.deepcopy(qualifiers),
            "hash": "f_stable",
            "location": "1..3",
            "record_location": "record:1..3:+",
        },
        "location_parts": [{"start": 1, "end": 3, "strand": 1}],
        "nucleotide_sequence": "ATG",
        "amino_acid_sequence": "M",
        "sequence_warnings": [],
    }
    extracted = copy.deepcopy(biological)
    extracted.pop("feature_index")
    extracted["rendered_feature_svg_id"] = "f_stable_record_1"
    payload = build_session_json(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 24),
        canonical_request=_canonical_request("linear"),
    )
    payload["features"] = {
        "extractedFeatures": [extracted],
        "biologicalFeatures": [biological],
        "selectedFeatureRecordIdx": 0,
    }
    expected_features = copy.deepcopy(payload["features"])
    path = tmp_path / "compact-feature-catalog.gbdraw-session.json"

    write_session_json(path, payload)

    on_disk = json.loads(path.read_text(encoding="utf-8"))
    assert "extractedFeatures" not in on_disk["features"]
    assert on_disk["features"]["featureCatalog"] == {
        "schema": 1,
        "encoding": "biological-authority-v1",
        "profile": "rich-v1",
        "extracted": [[0, "f0", "f_stable_record_1"]],
    }
    assert "stable_svg_id" not in on_disk["features"]["biologicalFeatures"][0]
    loaded_features = load_session(path)["features"]
    assert loaded_features == expected_features
    assert load_session_document(path).to_dict()["features"] == expected_features

    malformed_schema = copy.deepcopy(on_disk)
    malformed_schema["features"]["featureCatalog"]["schema"] = True
    path.write_text(
        json.dumps(
            malformed_schema,
            ensure_ascii=False,
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValidationError, match="Invalid compact session feature catalog"):
        load_session(path)

    null_catalog = copy.deepcopy(on_disk)
    null_catalog["features"]["featureCatalog"] = None
    path.write_text(
        json.dumps(null_catalog, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    with pytest.raises(ValidationError, match="Invalid compact session feature catalog"):
        load_session(path)
    with pytest.raises(ValidationError, match="Invalid compact session feature catalog"):
        write_session_json(tmp_path / "rejected-session.json", null_catalog)

    malformed_qualifier = copy.deepcopy(on_disk)
    malformed_qualifier["features"]["biologicalFeatures"][0]["qualifiers"][
        "protein_id"
    ] = [None]
    path.write_text(
        json.dumps(
            malformed_qualifier,
            ensure_ascii=False,
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    with pytest.raises(
        ValidationError,
        match="feature qualifiers must contain string arrays",
    ):
        load_session(path)

    duplicate_reference = copy.deepcopy(on_disk)
    duplicate_reference["features"]["featureCatalog"]["extracted"].append(
        copy.deepcopy(
            duplicate_reference["features"]["featureCatalog"]["extracted"][0]
        )
    )
    path.write_text(
        json.dumps(
            duplicate_reference,
            ensure_ascii=False,
            separators=(",", ":"),
        ),
        encoding="utf-8",
    )
    with pytest.raises(
        ValidationError,
        match="Invalid compact extracted-feature reference",
    ):
        load_session(path)

    on_disk["features"]["featureCatalog"]["extracted"][0][0] = 99
    path.write_text(
        json.dumps(on_disk, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    with pytest.raises(
        ValidationError,
        match="Invalid compact extracted-feature reference",
    ):
        load_session(path)


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


def test_v36_validates_mixed_protein_and_nucleotide_raw_cache() -> None:
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
    ("mode", "include_identity"),
    (
        ("orthogroup", True),
        ("collinear", True),
        ("collinear", False),
    ),
)
def test_v36_accepts_strict_zero_hit_derived_results(
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
def test_v36_rejects_near_miss_zero_hit_derived_results(mutate) -> None:
    manifest = _protein_identity_manifest()
    derived_payload = _zero_hit_derived_payload("collinear")
    mutate(derived_payload)

    with pytest.raises(ValidationError, match="unresolved protein references"):
        _current_derived_session(
            manifest,
            derived_payload,
            mode="collinear",
        )


def test_v36_derived_cache_accepts_compound_collinear_runtime_references() -> None:
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
def test_v36_derived_cache_validates_runtime_handles_in_unit_ids(
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


def test_v36_derived_cache_allows_non_protein_collinearity_unit_ids() -> None:
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
def test_v36_derived_cache_rejects_unresolved_compound_edge_references(
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
def test_v36_derived_cache_rejects_non_string_compound_edge_references(
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
def test_v36_derived_cache_rejects_embedded_legacy_references(
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


def test_v36_rejects_legacy_protein_entry_in_current_cache() -> None:
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


def test_v36_writer_quarantines_v27_to_v33_protein_artifacts(
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

    assert promoted["version"] == 36
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


def test_linear_protein_instance_keys_use_canonical_record_keys_not_file_metadata() -> None:
    source = {
        "renderRequest": {
            "records": [
                {"recordKey": "stable-left"},
                {"recordKey": "stable-right"},
            ]
        },
        "files": {
            "linearSeqs": [
                {"gb": {"name": "renamed-a.gb", "lastModified": 999}},
                {"gb": {"name": "renamed-b.gb", "lastModified": 123}},
            ]
        },
    }

    assert _record_instance_keys_for_web_losat(
        source_session=source,
        record_count=2,
    ) == ("stable-left", "stable-right")
    assert _record_instance_keys_for_web_losat(
        source_session=None,
        record_count=2,
    ) == ("record-1", "record-2")


def test_linear_source_cache_restores_legacy_candidate_original_entries() -> None:
    current = {
        **_current_protein_cache_entry(),
        "key": "current-protein-key",
    }
    legacy = _legacy_protein_cache_entry()
    session = {
        "losatCache": {"entries": [current]},
        "legacyArtifacts": {
            "proteinRawCandidates": {
                "schema": 1,
                "entries": [
                    {
                        "state": "pending",
                        "originalEntry": legacy,
                        "rejectionReason": None,
                    }
                ],
            }
        },
    }

    assert _source_session_losat_entries(session) == (current, legacy)
    assert _source_session_legacy_protein_candidates(session) == (
        {
            "state": "pending",
            "originalEntry": legacy,
            "rejectionReason": None,
        },
    )


def test_future_session_version_fails() -> None:
    session = _minimal_session({})
    session["version"] = CURRENT_SESSION_VERSION + 1

    with pytest.raises(ValidationError, match="newer"):
        validate_session(session)


@pytest.mark.parametrize("version", (34, 35))
def test_branch_internal_session_versions_are_unsupported(version: int) -> None:
    session = _minimal_session({})
    session["version"] = version

    with pytest.raises(ValidationError, match=f"Unsupported session version: {version}"):
        validate_session(session)
    assert version not in SUPPORTED_SESSION_VERSIONS


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


def test_canonical_session_version_32_remains_supported() -> None:
    session = {
        "format": SESSION_FORMAT,
        "version": 32,
        "renderRequest": {},
        "resources": {},
    }

    validate_session(session)
    assert 32 in SUPPORTED_SESSION_VERSIONS


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
        {"c_gb": _file_entry("input.gb", b"LOCUS       TEST\n")}
    )

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
    session["config"]["adv"] = {
        "blastSource": "losat",
        "losatProgram": "blastp",
        "min_bitscore": 25,
        "evalue": "1e-4",
        "identity": 55,
        "alignment_length": 30,
    }
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
    assert spec.args[spec.args.index("--collinear_max_gene_gap") + 1] == "3"
    assert spec.args[spec.args.index("--collinear_max_diagonal_drift") + 1] == "4"
    assert spec.args[spec.args.index("--collinear_max_conflicts_in_merge_gap") + 1] == "6"
    assert spec.args[spec.args.index("--collinear_search_scope") + 1] == "all"
    assert spec.args[spec.args.index("--collinear_color_mode") + 1] == "orientation_identity"
    assert spec.args[spec.args.index("--show_labels") + 1] == "orthogroup_top"


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
        "--scale_style",
        "ruler",
        "--palette",
        "ajisai",
        "--show_gc",
        "--show_skew",
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
    assert config["form"]["prefix"] == "out"
    assert config["form"]["align_center"] is True
    assert config["form"]["separate_strands"] is True
    assert config["form"]["scale_style"] == "ruler"
    assert config["form"]["show_gc"] is True
    assert config["form"]["show_skew"] is True
    assert config["form"]["show_labels_linear"] == "orthogroup_top"
    assert config["adv"]["pairwise_match_style"] == "curve"
    assert config["palette"] == "ajisai"
    assert config["blastSource"] == "losat"
    assert config["losatProgram"] == "blastp"
    assert config["losat"]["threadsPerJob"] == "32"
    assert config["losat"]["blastp"]["mode"] == "orthogroup"
    assert "collinearBlockMergeGap" not in config["losat"]["blastp"]
    assert "collinearSingletonMergeGap" not in config["losat"]["blastp"]
    assert config["cliOptions"]["rawArgs"] == list(args)
    assert config["cliOptions"]["byKey"]["protein_blastp_mode"] == ["orthogroup"]
    assert config["cliOptions"]["byKey"]["losatp_threads"] == ["32"]
    assert config["cliOptions"]["byKey"]["palette"] == ["ajisai"]
    assert config["cliOptions"]["byKey"]["gbk"] == [
        ["AP027078.gb", "AP027131.gb", "AP027133.gb", "AP027132.gb", "NZ_CP006932.gb"]
    ]


def test_cli_session_config_populates_safe_linear_row_fields() -> None:
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

    seqs = payload["files"]["linearSeqs"]
    assert seqs[0]["definition"] == "Alpha"
    assert seqs[0]["record_subtitle"] == "Alpha subtitle"
    assert seqs[0]["region_record_id"] == "RecA"
    assert seqs[0]["region_reverse"] is True
    assert seqs[1]["definition"] == "Beta"
    assert seqs[1]["record_subtitle"] == "Beta subtitle"
    assert seqs[1]["region_record_id"] == "RecB"
    assert seqs[1]["region_start"] == 10
    assert seqs[1]["region_end"] == 20
    assert seqs[1]["region_reverse"] is True
    assert payload["orthogroupState"]["selectedOrthogroupAlignmentFeature"] == "target_feature"


def test_cli_session_config_omits_ambiguous_multi_record_row_fields() -> None:
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

    seq = payload["files"]["linearSeqs"][0]
    assert seq["definition"] == ""
    assert seq["region_start"] is None
    assert seq["region_end"] is None
    assert seq["region_reverse"] is False


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


def test_session_hyphen_aliases_remain_compatible() -> None:
    request = parse_session_pre_args(
        [
            "--session",
            "session.gbdraw-session.json",
            "--session-output",
            "roundtrip.gbdraw-session.json",
        ],
        mode="circular",
    )

    assert request is not None
    assert request.save_session is True
    assert request.session_output == "roundtrip.gbdraw-session.json"
    assert strip_session_output_args(
        [
            "--gbk",
            "input.gb",
            "--save-session",
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
    circular_main(
        [
            "--gbk",
            str(examples_dir / "MjeNMV.gb"),
            "-o",
            str(output_prefix),
            "-f",
            "svg",
            "--save_session",
        ]
    )

    svg_path = output_prefix.with_suffix(".svg")
    session_path = output_prefix.with_suffix(".gbdraw-session.json")
    assert svg_path.exists()
    assert session_path.exists()

    payload = load_session(session_path)
    assert payload["format"] == SESSION_FORMAT
    assert payload["version"] == CURRENT_SESSION_VERSION
    assert "files" not in payload
    assert payload["resources"]["record-1-genbank"]["data"]
    assert "<svg" in payload["results"][0]["content"]
    assert payload["cliInvocation"]["mode"] == "circular"
    assert payload["cliInvocation"]["fileBindings"][0]["slot"] == "files.c_gb"

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
