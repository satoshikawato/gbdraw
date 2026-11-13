from __future__ import annotations

import base64
import copy
import gzip
import json
import hashlib
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
import gbdraw.session_container as session_container_module
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
    linear_main,
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
    CircularBatchRequest,
    CircularDiagramRequest,
    GenBankInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from tests.utils.session_files import (
    build_session_payload,
    read_session_payload,
    embedded_session,
    json_document,
    legacy_file_tree,
    load_session_mapping,
    write_session_file,
    write_unchecked_session_file,
)
from gbdraw.session_io import (
    CURRENT_FEATURE_CATALOG_SCHEMA,
    SEQUENCE_FEATURE_CATALOG_SCHEMA,
    CURRENT_SESSION_VERSION,
    DEPTH_FILE_ENCODING,
    LOSAT_DERIVED_CACHE_SCHEMA,
    MODE_SCOPED_SESSION_MAX_VERSION,
    MODE_SCOPED_SESSION_MIN_VERSION,
    ROW_KEYED_LEGEND_SESSION_MIN_VERSION,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
    SESSION_FORMAT,
    SUPPORTED_SESSION_VERSIONS,
    SessionBuildContext,
    _json_values_equal,
    compact_session_feature_catalog,
    decode_depth_payload,
    legend_edits_by_row_key,
    load_session,
    session_legend_edits_by_row_key,
    materialize_embedded_file,
    migrate_legacy_repeat_feature_shape_args,
    normalize_current_session_artifacts,
    migrate_persisted_web_state_field_names,
    session_to_cli_args,
    validate_current_web_state_field_names,
    SessionResources,
    validate_session,
)
from tests.utils.session_drawing_view import build_mode_scoped_session_json, drawing_view, view_at
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
from gbdraw.io.regions import parse_region_spec
from gbdraw.session import (
    SessionDocument,
    SessionFormatError,
    SessionVersionError,
    build_session_document,
    load_session_document,
    save_session_document,
)
from gbdraw.api.request_render import render_request


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
    output_path = tmp_path / "saved.gbdraw"
    old_temp_path = output_path.with_name(
        f".{output_path.name}.{os.getpid()}.tmp"
    )
    protected_path = tmp_path / "protected.txt"
    protected_path.write_text("keep", encoding="utf-8")
    try:
        old_temp_path.symlink_to(protected_path)
    except OSError as exc:
        pytest.skip(f"file symlinks are unavailable: {exc}")

    write_session_file(output_path, _minimal_session({}))

    assert protected_path.read_text(encoding="utf-8") == "keep"
    assert old_temp_path.is_symlink()
    assert output_path.is_file()


def test_write_session_json_no_replace_commit_preserves_late_target(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_path = tmp_path / "saved.gbdraw"
    real_commit = session_container_module.commit_staged_output_file

    def create_target_immediately_before_commit(
        staged_path: Path,
        target_path: Path,
        *,
        overwrite: bool,
    ) -> None:
        target_path.write_text("keep", encoding="utf-8")
        real_commit(staged_path, target_path, overwrite=overwrite)

    monkeypatch.setattr(
        session_container_module,
        "commit_staged_output_file",
        create_target_immediately_before_commit,
    )

    document, sources = session_io_module.split_resource_bytes(_current_session())
    with pytest.raises(ValidationError, match="overwrite=True"):
        session_container_module.write_session_container(
            document,
            sources,
            output_path,
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
    return build_session_payload(
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
    table_form: str = "rows",
) -> dict:
    # `rows: []` is the zero-hit form derived entries carried from #312 (17e2c9dee,
    # first-parent main) until D-66; the current writer emits a header-only
    # `canonicalTsv`. Admission accepts exactly one of the two (D-69).
    pair: dict = {
        "pair_index": 0,
        "query_index": 0,
        "subject_index": 1,
        "tsv": "",
        "hit_count": 0,
    }
    if table_form == "rows":
        pair["rows"] = []
    else:
        pair["canonicalTsv"] = "query\tsubject\tidentity\n"
    payload = {
        "pairs": [pair],
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
    # Count full checks; the readers of the proof call the validator, which returns it.
    real_check = protein_colinearity_module._check_protein_identity_manifest
    validation_calls = 0

    def count_validation(value):
        nonlocal validation_calls
        validation_calls += 1
        return real_check(value)

    monkeypatch.setattr(
        protein_colinearity_module,
        "_check_protein_identity_manifest",
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


_COPY_CROP_SESSION = (
    Path(__file__).parent / "fixtures" / "sessions" / "BGC0000708-BGC0000709-copy-crop.v51-web.gbdraw"
)


def test_web_session_with_colliding_record_instance_keys_loads() -> None:
    """OV-3273: a copied Linear drawing whose region changed keeps the card uids of
    the original, so one record instance key names other proteins in each drawing's
    manifest. The Web Save writes into each drawing only the protein entries its own
    manifest resolves, and the CLI admits it. tests/web/mode-scoped-session.test.mjs
    checks that a Load-then-Save of this Session stays so."""
    session = read_session_payload(_COPY_CROP_SESSION)
    for drawing in session["drawings"]:
        proteins = [
            entry for entry in drawing["losatCache"]["entries"] if "queryRecordInstanceKey" in entry
        ]
        assert len(proteins) == 4, drawing["id"]
        for entry in proteins:
            assert validate_protein_raw_entry_references(entry, drawing["proteinIdentityManifest"]), (
                drawing["id"],
                entry["key"],
            )
    load_session(_COPY_CROP_SESSION)


def test_web_session_refuses_a_drawing_with_another_drawings_protein_entries(tmp_path: Path) -> None:
    """The refusal OV-3273's Save ran into: a drawing that holds the other drawing's
    protein entries, which its own manifest does not resolve."""
    session = read_session_payload(_COPY_CROP_SESSION)
    drawings = {drawing["id"]: drawing for drawing in session["drawings"]}
    drawings["linear-2"]["losatCache"] = copy.deepcopy(drawings["linear"]["losatCache"])
    path = tmp_path / "mixed.gbdraw"
    write_unchecked_session_file(path, session)
    with pytest.raises(ValidationError, match="does not resolve through the manifest"):
        load_session(path)


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
    sidecar = tmp_path / "diagram.gbdraw"
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
                "anchorProfile": {
                    "precision": "exact",
                    "operator": "single",
                    "partOrder": "biological",
                    "strand": "+",
                },
            },
            {
                "svg_id": "feature-2",
                "stable_feature_id": "feature-2",
                "record_idx": 1,
                "anchorProfile": {
                    "precision": "exact",
                    "operator": "single",
                    "partOrder": "biological",
                    "strand": "+",
                },
            },
        ),
        orthogroup_metadata=(
            {
                "id": "og_1",
                "member_count": 2,
                "members": [
                    {
                        "recordIndex": 0,
                        "stableFeatureSvgId": "feature-1",
                        "proteinId": "protein-1",
                    },
                    {
                        "recordIndex": 1,
                        "stableFeatureSvgId": "feature-2",
                        "proteinId": "protein-2",
                    },
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
    payload = drawing_view(load_session(sidecar))
    # A CLI-written Session has no Web draft: no flat draft and no mode settings.
    assert "features" not in payload
    assert payload["modes"] == {"linear": {}}
    assert payload["orthogroupState"] == {}
    catalog = payload["editorState"]["featureCatalog"]
    assert catalog["schema"] == CURRENT_FEATURE_CATALOG_SCHEMA
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
    from gbdraw.render.output_paths import CLI_OVERWRITE_OPTION, OutputExistsError

    # OV-1240: the refusal is one error; the CLI prints it with --overwrite.
    with pytest.raises(OutputExistsError, match="Session output already exists") as refused:
        save_session_sidecar_if_requested(
            save_session=True,
            session_output=str(sidecar),
            output_prefix=str(output_prefix),
            run_result=run_result,
            cmd_args=(),
        )
    assert "--overwrite" in refused.value.message(CLI_OVERWRITE_OPTION)
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
    # The Web reads these values from `contracts/schema.generated.js`.
    assert CURRENT_SESSION_VERSION == 51 and ROW_KEYED_LEGEND_SESSION_MIN_VERSION == 47
    assert MODE_SCOPED_SESSION_MIN_VERSION == 46
    assert SUPPORTED_SESSION_VERSIONS == frozenset(
        {27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 44, 46, 47, 48, 49, 50, CURRENT_SESSION_VERSION}
    )


# The reading of a Session 46 Legend's edits by row key, shared with the Web
# (tests/web/legend-row-key-vectors.test.mjs).
LEGEND_ROW_KEY_VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "sessions" / "legend-row-key-vectors.json").read_text(encoding="utf-8")
)


@pytest.mark.parametrize("case", LEGEND_ROW_KEY_VECTORS["cases"], ids=lambda case: case["name"])
def test_session_46_legend_edits_move_to_the_row_key(case: dict) -> None:
    assert LEGEND_ROW_KEY_VECTORS["readsSessionVersion"] == MODE_SCOPED_SESSION_MIN_VERSION
    source = copy.deepcopy(case["legend"])
    assert legend_edits_by_row_key(source) == (case["expected"], case["outcomes"])
    assert source == case["legend"]


# The Session 46 that first-parent main 9839de8e wrote after renaming GC skew (+)
# to Z, GC content to GC skew (+) (then red), and GC skew (-) to W (blue,
# deleted), with a stroke on Z (see its .provenance.json).
LEGEND_RENAMED_ROWS_V46 = (
    Path(__file__).parent / "fixtures" / "sessions" / "legend-renamed-rows.v46.gbdraw-session.json.gz"
)
# The edits as Session 46 Load read them: the drawn stroke `none` is no color.
ROW_KEYED_LEGEND_EDITS = {
    "colorOverrides": {"GC content": "#ff0000", "GC skew (-)": "#0000ff"},
    "strokeOverrides": {"GC skew (+)": {
        "originalStrokeColor": None, "originalStrokeWidth": 0, "strokeWidth": 2, "strokeColor": "#00aa00",
    }},
}


# The moves its reading reports, in the order the Web reports them.
LEGEND_RENAMED_ROWS_V46_MOVES = (
    ("colorOverrides", "GC skew (+)", "GC content"),
    ("colorOverrides", "W", "GC skew (-)"),
    ("strokeOverrides", "Z", "GC skew (+)"),
)


def _legend_edits(session: dict) -> dict:
    legend = session["modes"]["circular"]["editorState"]["legend"]
    return {store: legend[store] for store in ROW_KEYED_LEGEND_EDITS}


def test_session_46_upgrade_and_cli_resave_key_the_legend_edits_by_row(tmp_path: Path) -> None:
    from gbdraw.session import upgrade_session_document

    source = load_session(LEGEND_RENAMED_ROWS_V46)
    assert source["version"] == 46
    assert _legend_edits(source)["colorOverrides"] == {"GC skew (+)": "#ff0000", "W": "#0000ff"}
    assert session_legend_edits_by_row_key(source)[1] == [
        {"slice": "circular", "store": store, "fromKey": shown, "toKey": key, "outcome": "moved"}
        for store, shown, key in LEGEND_RENAMED_ROWS_V46_MOVES
    ]

    upgraded_document = upgrade_session_document(LEGEND_RENAMED_ROWS_V46).document
    written = upgraded_document.to_dict()
    assert written["version"] == CURRENT_SESSION_VERSION
    upgraded = drawing_view(written)
    assert _legend_edits(upgraded) == ROW_KEYED_LEGEND_EDITS
    expected_legend = source["modes"]["circular"]["editorState"]["legend"] | ROW_KEYED_LEGEND_EDITS
    assert upgraded["modes"]["circular"]["editorState"]["legend"] == expected_legend
    # A current Session is read as it is.
    assert _legend_edits(drawing_view(upgrade_session_document(upgraded_document).document.to_dict())) == ROW_KEYED_LEGEND_EDITS

    sidecar = tmp_path / "replay.gbdraw"
    circular_main(["--session", str(LEGEND_RENAMED_ROWS_V46), "--output", str(tmp_path / "replay"),
                   "--format", "svg", "--session_output", str(sidecar)])
    resaved = load_session(sidecar)
    assert resaved["version"] == CURRENT_SESSION_VERSION
    assert _legend_edits(drawing_view(resaved)) == ROW_KEYED_LEGEND_EDITS


def test_current_session_feature_catalog_is_single_and_lossless(
    tmp_path: Path,
) -> None:
    catalog = {
        "schema": SEQUENCE_FEATURE_CATALOG_SCHEMA,  # a Session 46-50 stores the sequences
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
                        "drawnSelector": {
                            "hash": "rendered-feature",
                            "location": "2..3:+",
                            "recordLocation": "public-record:2..3:+",
                        },
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
                        "anchorProfile": {
                            "precision": "exact",
                            "operator": "single",
                            "partOrder": "biological",
                            "strand": "unstranded",
                        },
                    }
                ],
                "orthogroups": [],
                "annotations": [],
                "comparisonMatches": [],
            }
        ],
    }
    payload = build_mode_scoped_session_json(
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
    payload["modes"] = {"linear": {"config": {"palette": "default"}}}
    path = tmp_path / "schema-3-feature-catalog.gbdraw"

    write_session_file(path, payload)

    on_disk = read_session_payload(path)
    assert on_disk["editorState"]["featureCatalog"] == catalog
    assert on_disk["modes"] == {"linear": {"config": {"palette": "default"}}}
    assert "features" not in on_disk
    assert on_disk["orthogroupState"] == {}
    assert load_session(path) == json_document(on_disk)
    assert load_session_document(path).to_dict() == json_document(on_disk)

    unknown_field = copy.deepcopy(on_disk)
    unknown_field["branchOnlyState"] = {}
    with pytest.raises(
        ValidationError,
        match="unclassified top-level field.*branchOnlyState",
    ):
        write_session_file(tmp_path / "unknown-field.json", unknown_field)

    missing_results = copy.deepcopy(on_disk)
    missing_results.pop("results")
    with pytest.raises(
        ValidationError,
        match="requires a results array",
    ):
        write_session_file(tmp_path / "missing-results.json", missing_results)

    empty_draft = copy.deepcopy(on_disk)
    empty_draft["results"] = []
    empty_draft["editorState"]["featureCatalog"] = None
    write_session_file(tmp_path / "empty-draft.json", empty_draft)

    missing_empty_catalog = copy.deepcopy(empty_draft)
    missing_empty_catalog["editorState"].pop("featureCatalog")
    with pytest.raises(
        ValidationError,
        match="requires editorState.featureCatalog",
    ):
        write_session_file(
            tmp_path / "missing-empty-catalog.json",
            missing_empty_catalog,
        )

    invalid_container = copy.deepcopy(on_disk)
    invalid_container["orthogroupState"] = []
    with pytest.raises(
        ValidationError,
        match=r"Session orthogroupState must be an object",
    ):
        write_session_file(
            tmp_path / "invalid-orthogroupState.json",
            invalid_container,
        )
    # Session 46 keeps the Web draft in the mode slices only.
    flat_draft = copy.deepcopy(on_disk)
    flat_draft["features"] = {}
    with pytest.raises(ValidationError, match="it cannot contain features") as excinfo:
        write_session_file(tmp_path / "flat-draft.json", flat_draft)
    assert excinfo.value.diagnostic == {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}

    interactive_result = copy.deepcopy(on_disk)
    interactive_result["results"][0]["name"] = "out.interactive.svg"
    interactive_result["editorState"]["featureCatalog"]["items"][0][
        "resultName"
    ] = "out.interactive.svg"
    with pytest.raises(
        ValidationError,
        match="named plain SVG",
    ):
        write_session_file(
            tmp_path / "interactive-result.json",
            interactive_result,
        )

    missing_catalog = copy.deepcopy(on_disk)
    missing_catalog["editorState"].pop("featureCatalog")
    with pytest.raises(
        ValidationError,
        match="Results require editorState.featureCatalog",
    ):
        write_session_file(
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
        write_session_file(
            tmp_path / "duplicate-feature.json",
            duplicate_feature,
        )

    # A Session 44 draft cannot duplicate the catalog in its flat features.
    for field in ("biologicalFeatures", "featureSelectorSafetyScope", "featureRecordIds"):
        duplicated_legacy = copy.deepcopy(on_disk)
        duplicated_legacy["features"] = {field: []}
        with pytest.raises(
            ValidationError,
            match="derived feature payloads",
        ):
            session_io_module._validate_current_feature_payloads(duplicated_legacy, 44)


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
    payload = build_mode_scoped_session_json(
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
            {"adv": {"linear_show_accession": True}},
            "Web state field adv.linear_show_accession is obsolete; "
            "use the selected visibility mode.",
        ),
        (
            {"adv": {"linear_show_length": True}},
            "Web state field adv.linear_show_length is obsolete; "
            "use the selected visibility mode.",
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
    payload = build_mode_scoped_session_json(
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
    payload["modes"] = {"linear": {"config": config}}

    with pytest.raises(ValidationError, match=re.escape(message)):
        validate_session(payload)
    output_path = tmp_path / "obsolete-current-session.gbdraw"
    with pytest.raises(ValidationError, match=re.escape(message)):
        write_session_file(output_path, payload)
    assert not output_path.exists()

    payload["version"] = 33
    payload["config"] = payload.pop("modes")["linear"]["config"]
    payload.get("webFiles", {}).pop("bindings", None)
    payload["renderRequest"]["schema"] = 2
    validate_session(payload)


def test_released_web_state_field_names_migrate_copy_on_write() -> None:
    source = {
        "adv": {
            "linear_accession_visibility": "HIDE",
            "linear_show_accession": True,
            "linear_show_length": False,
            "depth_tick_interval": 10,
            "depth_large_tick_interval": 20,
            "depth_tracks": [
                {"tick_interval": 5},
                {"tick_interval": 6, "large_tick_interval": 7},
            ],
        },
        "recordDisplayDrafts": [
            {
                "scope": "circular",
                "sourceUid": "source-1",
                "selector": "#1",
                "recordId": "record-1",
                "topologyOverride": None,
                "startCoordinate": 3,
            }
        ],
        "losat": {
            "blastp": {
                "collinearMaxGeneGap": 2,
                "collinearMaxUnitGap": 3,
            }
        },
    }

    migrated = migrate_persisted_web_state_field_names(source)

    assert isinstance(migrated, dict)
    assert migrated["adv"]["linear_accession_visibility"] == "hide"
    assert migrated["adv"]["linear_length_visibility"] == "hide"
    assert "linear_show_accession" not in migrated["adv"]
    assert "linear_show_length" not in migrated["adv"]
    assert migrated["adv"]["depth_large_tick_interval"] == 20
    assert "depth_tick_interval" not in migrated["adv"]
    assert migrated["adv"]["depth_tracks"] == [
        {"large_tick_interval": 5},
        {"large_tick_interval": 7},
    ]
    assert migrated["losat"]["blastp"]["collinearMaxUnitGap"] == 3
    assert "collinearMaxGeneGap" not in migrated["losat"]["blastp"]
    assert migrated["recordDisplayDrafts"][0]["reverseComplementOverride"] is None
    assert migrated["recordDisplayDrafts"][0]["anchorIntent"] is None
    assert source["adv"]["depth_tick_interval"] == 10
    assert source["adv"]["depth_tracks"][0]["tick_interval"] == 5
    assert source["losat"]["blastp"]["collinearMaxGeneGap"] == 2
    assert "reverseComplementOverride" not in source["recordDisplayDrafts"][0]
    validate_current_web_state_field_names(migrated)


@pytest.mark.parametrize("schema", [8, CANONICAL_REQUEST_SCHEMA])
def test_alignment_session_requires_reset_receipt_for_request_schema_8_and_later(schema: int) -> None:
    # Same rule as services/session-authority.js: renderRequest.schema >= 8.
    session = {
        "renderRequest": {
            "schema": schema,
            "mode": "linear",
            "layout": {"similarityAlignment": {"records": [], "reference": {"recordKey": "a"}}},
        },
        "editorState": {"featureCatalog": None},
    }
    with pytest.raises(ValidationError, match="requires editorState.alignmentResetReceipt"):
        session_io_module._validate_alignment_reset_receipt(session, {})
    session["editorState"]["alignmentResetReceipt"] = None
    session_io_module._validate_alignment_reset_receipt(session, {})


def _receipt_session(schema: int, placement: str, binding: str) -> dict[str, object]:
    """A Linear alignment Session whose record ``b`` is cropped 2..10, with its record selector stored where
    ``placement`` says: in the region, as the v0.14.5 Web wrote it, or on the record."""

    source = b"source bytes remain biological authority"

    def anchor(key: str) -> dict[str, object]:
        return {"recordKey": key, "biologicalFeatureId": f"gene-{key}", "sourceFeatureIndex": 1, "stableFeatureSvgId": None}

    def record(key: str) -> dict[str, object]:
        selector = {"kind": "recordId", "value": key}
        span = {"start": 2, "end": 10, "reverseComplement": False}
        cropped = key == "b"
        stored: dict[str, object] = {
            "recordKey": key, "source": {"kind": "genbank", "resourceId": "source"}, "cardinality": "exactly_one",
            "selector": None if cropped and placement == "region" else selector,
            "display": {"isCircular": None, "startCoordinate": None}, "presentation": {"reverseComplement": False},
        }
        if schema >= 10:
            stored["windows"] = [{**span, "column": None}] if cropped else []
        else:
            stored["region"] = {"selector": selector if placement == "region" else None, **span} if cropped else None
        return stored

    plan = {"schema": 2, "groupId": "og-1", "reference": anchor("a"), "records": [
        {"recordKey": "a", "status": "reference", "rationale": "reference", "anchor": anchor("a")},
        {"recordKey": "b", "status": "aligned", "rationale": "user_selected", "anchor": anchor("b")},
        {"recordKey": "c", "status": "skipped", "rationale": "skipped_by_user", "anchor": None},
    ]}
    return {
        "renderRequest": {"schema": schema, "mode": "linear", "records": [record(key) for key in "abc"],
                          "layout": {"recordTranslations": [], "similarityAlignment": plan}},
        "editorState": {"featureCatalog": None,
                        "alignmentResetReceipt": {"binding": binding, "directions": [], "referenceDeltaX": None}},
        "resources": {"source": {"name": "record.gbk", "kind": "genbank", "encoding": "base64",
                                 "data": base64.b64encode(source).decode(), "size": len(source)}},
    }


# OV-2772: a receipt binds its request as the file stored it. The bindings are those the v0.14.5 Web writes for
# request schema 9 (the selector in the region, or on the record) and the one the Web rebinds a promoted request to
# (tests/web/alignment-reset-receipt.test.mjs).
_V0145_REGION_SELECTOR_BINDING = "fe2664f3a1b0107bdf2e54ec4c0636a064039fb4484a758fa8109c0648964da4"
_RECORD_SELECTOR_BINDING = "cf9c72f1cdbb85e4b89178634a5a22c54a7a6250e5357b985a457c2b3d8317b1"


@pytest.mark.parametrize(("schema", "placement", "binding"), [
    (9, "region", _V0145_REGION_SELECTOR_BINDING),
    (9, "record", _RECORD_SELECTOR_BINDING),
    (CANONICAL_REQUEST_SCHEMA, "record", _RECORD_SELECTOR_BINDING),
])
def test_alignment_reset_receipt_binds_the_request_as_stored(schema: int, placement: str, binding: str) -> None:
    session_io_module._validate_alignment_reset_receipt(
        *session_io_module.split_resource_bytes(_receipt_session(schema, placement, binding)))


def test_alignment_reset_receipt_of_a_region_selector_is_not_the_promoted_binding() -> None:
    with pytest.raises(ValidationError, match="binding changed"):
        session_io_module._validate_alignment_reset_receipt(
            *session_io_module.split_resource_bytes(_receipt_session(9, "region", _RECORD_SELECTOR_BINDING)))


def test_session_version_46_keys_feature_edits_by_source_identity_per_mode() -> None:
    # Design Q4: per-feature edits are identity rows. Session 46 keeps each
    # mode's rows in its own slice (PD-OI-090), so both modes keep their own
    # row for `record-1` and a row names no mode.
    payload = build_mode_scoped_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 10, 4),
        canonical_request=_canonical_request("circular"),
    )
    row = {
        "recordKey": "record-1",
        "biologicalFeatureId": "feature-1",
        "featureVisibility": "off",
        "labelVisibility": None,
        "labelText": "Edited label",
        "labelSourceText": "nad1",
    }
    source_only = {
        **row,
        "biologicalFeatureId": "feature-2",
        "featureVisibility": None,
        "labelText": None,
        "labelSourceText": "nad2",
    }
    other_mode = {**row, "labelText": None, "labelSourceText": None}

    def key(feature_id: str) -> str:
        return json.dumps(["record-1", feature_id], separators=(",", ":"))

    payload["modes"] = {
        "circular": {"features": {"featureOverrides": {key("feature-1"): row, key("feature-2"): source_only}}},
        "linear": {"features": {"featureOverrides": {key("feature-1"): other_mode}}},
    }

    assert payload["version"] == MODE_SCOPED_SESSION_MAX_VERSION
    validate_session(payload)

    triple_key = json.dumps(["circular", "record-1", "feature-1"], separators=(",", ":"))
    for broken_features in (
        {"featureVisibilityOverrides": {"rendered-1": "off"}},
        {"featureOverrides": {key("feature-1"): {**row, "recordKey": "record-2"}}},
        {"featureOverrides": {triple_key: row}},
        {"featureOverrides": {key("feature-1"): {**row, "scope": "circular"}}},
        {"featureOverrides": {key("feature-1"): {k: v for k, v in row.items() if k != "labelSourceText"}}},
        {"featureOverrides": {key("feature-1"): {**row, "featureVisibility": None, "labelText": None,
                                                 "labelSourceText": None}}},
        {"featureOverrides": {key("feature-1"): {**row, "labelText": "two\nlines"}}},
    ):
        broken = copy.deepcopy(payload)
        broken["modes"]["circular"]["features"] = broken_features
        with pytest.raises(ValidationError):
            validate_session(broken)


def test_session_version_46_keeps_feature_placement_drafts_of_each_mode() -> None:
    # A Linear Session keeps the Circular lane rows in the Circular slice;
    # each slice's lanes must be sides of its own mode. Session 44 rows (a
    # JSON pair) belonged to the Session's mode.
    payload = build_mode_scoped_session_json(
        SessionBuildContext(mode="linear", output_prefix="out", render_formats=("svg",)),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 10, 5),
        canonical_request=_canonical_request("linear"),
    )

    def placement(side: str | None) -> tuple[str, dict]:
        target = {"kind": "main"} if side is None else {"kind": "lane", "side": side, "level": 1}
        row = {"recordKey": "record-1", "biologicalFeatureId": f"f-{side}", "placement": target}
        return json.dumps(["record-1", f"f-{side}"], separators=(",", ":")), row

    payload["modes"] = {
        "circular": {"config": {"adv": {}, "featurePlacementOverrides": dict(
            placement(side) for side in ("outward", None)
        )}},
        "linear": {"config": {"adv": {}, "featurePlacementOverrides": dict([placement("above")])}},
    }
    validate_session(payload)

    wrong_side_key, wrong_side = placement("outward")
    above_key, above = placement("above")
    for broken in (
        {wrong_side_key: wrong_side},
        {json.dumps(["linear", "record-1", "f-above"], separators=(",", ":")): above},
        {above_key: {"scope": "linear", **above}},
        {wrong_side_key: above},
    ):
        invalid = copy.deepcopy(payload)
        invalid["modes"]["linear"]["config"]["featurePlacementOverrides"] = broken
        with pytest.raises(ValidationError):
            validate_session(invalid)

    legacy = copy.deepcopy(payload)
    legacy["version"] = 44
    legacy.pop("modes")
    legacy["editorState"]["featureCatalog"]["schema"] = 4
    legacy["renderRequest"]["schema"] = 8
    legacy["renderRequest"]["diagramOptions"].pop("featureOverrides", None)
    legacy_row = {"recordKey": "record-1", "biologicalFeatureId": "f1", "placement": {"kind": "main"}}
    legacy["config"] = {"adv": {}, "featurePlacementOverrides": {
        json.dumps(["record-1", "f1"], separators=(",", ":")): legacy_row
    }}
    validate_session(legacy)
    with pytest.raises(ValidationError, match="Session version 44 cannot contain modes"):
        validate_session({**legacy, "modes": {}})


def test_session_version_44_validates_record_rotation_draft_and_version_42_shape() -> None:
    payload = build_mode_scoped_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 9, 23),
        canonical_request=_canonical_request("circular"),
    )
    draft = {
        "sourceUid": "source-1",
        "selector": "#1",
        "recordId": "record-1",
        "topologyOverride": None,
        "startCoordinate": 3,
        "reverseComplementOverride": True,
        "anchorIntent": {
            "schema": 1,
            "recordKey": "record-1",
            "biologicalFeatureId": "feature-1",
            "placement": "anchor",
            "anchor": "five-prime",
            "offsetBp": -2,
            "orientForward": True,
        },
    }
    # Session 46: each mode slice holds its own rows, which name no mode.
    payload["modes"] = {"circular": {"config": {"adv": {}, "recordDisplayDrafts": [draft]}}}

    validate_session(payload)
    scoped = copy.deepcopy(payload)
    scoped["modes"]["circular"]["config"]["recordDisplayDrafts"][0]["scope"] = "circular"
    with pytest.raises(ValidationError, match="Invalid record display draft fields"):
        validate_session(scoped)

    legacy = copy.deepcopy(payload)
    legacy["version"] = 42
    legacy.pop("modes")
    legacy["editorState"]["featureCatalog"]["schema"] = 3
    legacy_draft = {"scope": "circular", **draft}
    legacy_draft.pop("reverseComplementOverride")
    legacy_draft.pop("anchorIntent")
    legacy["config"] = {
        "adv": {
            "linear_show_accession": True,
            "linear_show_length": False,
        },
        "recordDisplayDrafts": [legacy_draft],
    }
    validate_session(legacy)


def test_session_version_44_rejects_invalid_record_rotation_anchor_intent() -> None:
    payload = build_mode_scoped_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 9, 23),
        canonical_request=_canonical_request("circular"),
    )
    payload["modes"] = {"circular": {"config": {
        "adv": {},
        "recordDisplayDrafts": [
            {
                "sourceUid": "source-1",
                "selector": "#1",
                "recordId": "record-1",
                "topologyOverride": None,
                "startCoordinate": 3,
                "reverseComplementOverride": True,
                "anchorIntent": {
                    "schema": 1,
                    "recordKey": "record-1",
                    "biologicalFeatureId": "feature-1",
                    "placement": "feature-end",
                    "anchor": "five-prime",
                    "offsetBp": 0,
                    "orientForward": True,
                },
            }
        ],
    }}}

    with pytest.raises(ValidationError, match="Invalid record display anchor intent"):
        validate_session(payload)


@pytest.mark.parametrize("table_form", ("rows", "canonicalTsv"))
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
    table_form: str,
) -> None:
    manifest = _protein_identity_manifest()

    payload = _current_derived_session(
        manifest,
        _zero_hit_derived_payload(
            mode,
            include_identity=include_identity,
            table_form=table_form,
        ),
        mode=mode,
    )

    validate_session(payload)


@pytest.mark.parametrize(
    "member",
    ({"note": "arbitrary"}, {"collinearGroupScope": "invalid"}),
    ids=("arbitrary-field", "other-collinear-scope"),
)
def test_a_zero_hit_derived_result_may_hold_members_without_values(member: dict) -> None:
    """OV-1127: no member list; a member that holds no array or object entry names no protein."""

    payload = _zero_hit_derived_payload("collinear")
    payload.update(member)
    validate_session(_current_derived_session(_protein_identity_manifest(), payload, mode="collinear"))


def _replace_rows_with_table(payload: dict, table: object) -> None:
    pair = payload["pairs"][0]
    del pair["rows"]
    pair["canonicalTsv"] = table


@pytest.mark.parametrize(
    "mutate",
    (
        lambda payload: payload["pairs"][0].update({"rows": [{}]}),
        lambda payload: payload["pairs"][0].update({"canonicalTsv": "query\n"}),
        lambda payload: _replace_rows_with_table(payload, "query\nq1\n"),
        lambda payload: _replace_rows_with_table(payload, ["query"]),
        lambda payload: payload["pairs"][0].update({"tsv": "unexpected"}),
        lambda payload: payload["pairs"][0].update({"hit_count": 1}),
        lambda payload: payload["pairs"][0].update({"pair_index": "0"}),
        lambda payload: payload["pairs"][0].pop("subject_index"),
        lambda payload: payload.update({"pairs": []}),
        lambda payload: payload.update({"orthogroups": [{}]}),
        lambda payload: payload.update({"collinearGroups": [{}]}),
        lambda payload: payload.update({"collinearityBlocks": [{}]}),
        lambda payload: payload.update({"identity": None}),
        lambda payload: payload["identity"].update({"rawCacheKeys": [None]}),
        lambda payload: payload["identity"].update({"mode": "orthogroup"}),
    ),
    ids=(
        "nonempty-rows",
        "both-table-forms",
        "nonempty-canonical-table",
        "nonstring-canonical-table",
        "nonempty-tsv",
        "nonzero-hit-count",
        "noninteger-pair-index",
        "partial-record-indices",
        "missing-pair-results",
        "nonempty-orthogroups",
        "nonempty-collinear-groups",
        "nonempty-collinearity-blocks",
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
    payload = build_mode_scoped_session_json(
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
    payload = build_mode_scoped_session_json(
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
    source = build_mode_scoped_session_json(
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
    source.get("webFiles", {}).pop("bindings", None)
    source["renderRequest"]["schema"] = 2
    source["config"] = {"form": {}, "adv": {}}
    source["config"]["adv"]["depth_tick_interval"] = 10
    source["config"]["adv"]["depth_tracks"] = [{"tick_interval": 5}]
    source["config"]["losat"] = {
        "blastp": {"collinearMaxGeneGap": 2}
    }
    source["orthogroupState"] = {
        "selectedOrthogroupAlignmentFeature": "legacy-target",
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
        build_mode_scoped_session_json(context, **kwargs)

    promoted = build_mode_scoped_session_json(
        context,
        canonical_request=_canonical_request("linear"),
        **kwargs,
    )

    assert promoted["version"] == MODE_SCOPED_SESSION_MAX_VERSION
    assert promoted["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert "selectedOrthogroupAlignmentFeature" not in promoted["orthogroupState"]
    assert "config" not in promoted
    for mode in ("circular", "linear"):
        adv = promoted["modes"][mode]["config"]["adv"]
        assert adv["depth_large_tick_interval"] == 10
        assert adv["depth_tracks"] == [{"large_tick_interval": 5}]
    # LOSATP settings belong to the Linear slice.
    assert "blastp" not in promoted["modes"]["circular"]["config"].get("losat", {})
    assert (
        promoted["modes"]["linear"]["config"]["losat"]["blastp"]["collinearMaxUnitGap"]
        == 2
    )


def test_current_writer_resolves_noncanonical_request_before_encoding() -> None:
    record = SeqRecord(
        Seq("ATGC"),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(record),
                cardinality=RecordCardinality.ALL,
                presentation=RecordPresentation(reverse_complement=True),
            ),
        ),
    )

    payload = drawing_view(build_session_payload(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 7, 21),
        canonical_request=request,
    ))

    assert [record["recordKey"] for record in payload["renderRequest"]["records"]] == [
        "record-1"
    ]
    assert payload["renderRequest"]["records"][0]["presentation"][
        "reverseComplement"
    ] is False


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

    promoted = build_mode_scoped_session_json(
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
    rewritten = build_mode_scoped_session_json(
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
    # The writer builds the mode-scoped Session 48 and moves it to drawings last.
    assert promoted["version"] == MODE_SCOPED_SESSION_MAX_VERSION
    assert rewritten["version"] == MODE_SCOPED_SESSION_MAX_VERSION
    # The second write copies the first one's mode slices.
    assert rewritten["modes"] == promoted["modes"]
    for payload in (promoted, rewritten):
        assert "files" not in payload
        assert "config" not in payload and "features" not in payload
        circular, linear = (payload["modes"][mode] for mode in ("circular", "linear"))
        assert payload["webFiles"]["linearRecords"] == source["webFiles"][
            "linearRecords"
        ]
        assert payload["webFiles"]["bindings"]["schema"] == 2
        assert "linearCanonicalComparisons" not in payload["webFiles"]["bindings"]
        assert linear["config"]["linearComparisonPlan"] == {
            "mode": "adjacent",
            "defaultSource": "losat",
            "edges": [],
        }
        assert "linearComparisonPlan" not in circular["config"]
        assert "blastSource" not in linear["config"]
        assert "comparisons" not in linear["config"]["linearRecordLayout"]
        # The Legend edits go to the saved Result's mode (R1-3).
        assert linear["editorState"]["legend"] == source["editorState"]["legend"]
        assert "editorState" not in circular
        assert "legend" not in payload["editorState"]
        # The checked mode-scoped Session labels it as a Session 46-50 does.
        assert payload["editorState"]["featureCatalog"]["schema"] == SEQUENCE_FEATURE_CATALOG_SCHEMA
        assert len(payload["editorState"]["featureCatalog"]["items"]) == 1
        assert "features" not in circular and "features" not in linear
        assert (
            circular["config"]["adv"]["circular_track_slots"]
            == source["config"]["adv"]["circular_track_slots"]
        )
        assert "circular_track_slots" not in linear["config"]["adv"]
        # Web Load gives each Linear slot its own params (schema 2).
        assert linear["config"]["adv"]["linear_track_slots"] == [
            {**slot, "params": {}} for slot in source["config"]["adv"]["linear_track_slots"]
        ]
        # The shown mode keeps its flat value; the other mode its saved profile.
        assert linear["config"]["adv"]["identity"] == source["config"]["adv"]["identity"]
        assert circular["config"]["adv"]["identity"] == 88
        assert "modeProfiles" not in circular["config"]


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
    source = build_mode_scoped_session_json(
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
    source.get("webFiles", {}).pop("bindings", None)
    source["config"] = {"form": {}, "adv": {}}
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

    promoted = build_mode_scoped_session_json(
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

    linear_config = promoted["modes"]["linear"]["config"]
    assert "blastSource" not in linear_config
    assert "comparisons" not in linear_config["linearRecordLayout"]
    assert linear_config["linearComparisonPlan"] == {
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
    payload = build_mode_scoped_session_json(
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
    payload["modes"] = {"linear": {"config": {"form": {}, "adv": {}, "blastSource": "losat"}}}
    with pytest.raises(ValidationError, match="retired blastSource"):
        validate_session(payload)

    payload["modes"]["linear"]["config"].pop("blastSource")
    payload["modes"]["linear"]["config"]["linearRecordLayout"] = {
        "enabled": False,
        "recordGap": 24,
        "rows": [],
        "comparisons": [],
    }
    with pytest.raises(ValidationError, match="linearRecordLayout.comparisons"):
        validate_session(payload)


@pytest.mark.parametrize(
    "field",
    ("cli_circular_track_order", "cli_circular_track_slots"),
)
def test_current_session_rejects_retired_circular_track_draft_paths(
    field: str,
) -> None:
    payload = build_mode_scoped_session_json(
        SessionBuildContext(
            mode="circular",
            output_prefix="out",
            render_formats=("svg",),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 8, 23),
        canonical_request=_canonical_request("circular"),
    )
    payload["modes"] = {"circular": {"config": {"form": {}, "adv": {field: []}}}}

    with pytest.raises(ValidationError, match=rf"modes\.circular\.config\.adv cannot contain {field}"):
        validate_session(payload)
    flat = copy.deepcopy(payload)
    flat["version"] = 44
    flat["config"] = flat.pop("modes")["circular"]["config"]
    with pytest.raises(ValidationError, match=rf"config\.adv\.{field}"):
        session_io_module._validate_current_retired_active_config_paths(flat, 44)


def test_current_comparison_authority_validates_file_bindings() -> None:
    payload = build_mode_scoped_session_json(
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
    # Session 46 keeps the comparison draft in the Linear slice.
    payload["modes"] = {"linear": {"config": {"form": {}, "adv": {}}}}
    payload["modes"]["linear"]["config"]["linearComparisonPlan"] = {
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
    unsupported_schema["webFiles"]["bindings"]["schema"] = 99
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
    source = build_mode_scoped_session_json(
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
    source.get("webFiles", {}).pop("bindings", None)
    source["losatCache"] = {"entries": [legacy_raw]}
    source["losatDerivedCache"] = {"entries": [legacy_derived]}
    source.pop("proteinIdentityManifest")

    promoted = build_mode_scoped_session_json(
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

    assert promoted["version"] == MODE_SCOPED_SESSION_MAX_VERSION
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

    session_path = tmp_path / "promoted.json"
    write_session_file(session_path, promoted)
    assert load_session(session_path) == json_document(promoted)


def test_future_session_version_fails() -> None:
    session = _minimal_session({})
    session["version"] = max(SUPPORTED_SESSION_VERSIONS) + 1

    with pytest.raises(ValidationError, match="newer"):
        validate_session(session)


@pytest.mark.parametrize("version", [None, "45"])
def test_session_without_integer_version_fails_validation(
    tmp_path: Path,
    version: object,
) -> None:
    # The version check runs before the version-dependent feature-catalog expansion.
    session: dict[str, object] = {"format": SESSION_FORMAT, "config": {}}
    if version is not None:
        session["version"] = version
    session_path = tmp_path / "session.json"
    session_path.write_text(json.dumps(session), encoding="utf-8")

    with pytest.raises(ValidationError, match="Session version is required"):
        load_session(session_path)
    with pytest.raises(SessionVersionError, match="Session version is required"):
        load_session_document(session)
    with pytest.raises(ValidationError, match="Session version is required"):
        write_session_file(tmp_path / "written.json", session)


def test_current_session_rejects_legacy_files_but_version_39_accepts_them() -> None:
    session = build_mode_scoped_session_json(
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
    session.get("webFiles", {}).pop("bindings", None)
    validate_session(session)


@pytest.mark.parametrize("version", (34, 35, 36, 37, 38, 43))
def test_branch_internal_session_versions_are_rejected_at_read_and_rewrite_boundaries(
    version: int,
) -> None:
    session = _minimal_session({})
    session["version"] = version

    with pytest.raises(ValidationError, match=f"Unsupported session version: {version}"):
        validate_session(session)
    assert version not in SUPPORTED_SESSION_VERSIONS

    with pytest.raises(ValidationError, match=f"Unsupported session version: {version}"):
        build_session_payload(
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
    build = build_mode_scoped_session_json if session_version == 33 else build_session_payload
    session = build(
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
    if session_version == 33:
        # An older Session carries its request at the top level.
        session["version"] = 33
        session["renderRequest"]["schema"] = schema
    else:
        session["drawings"][0]["renderRequest"]["schema"] = schema
    message = f"Unsupported canonical renderRequest schema: {schema}"

    with pytest.raises(ValidationError, match=message):
        validate_session(session)
    with pytest.raises(SessionFormatError, match=message):
        load_session_mapping(session)

    input_path = (
        tmp_path
        / f"version-{session_version}-schema-{schema}-input.gbdraw"
    )
    write_unchecked_session_file(input_path, session)
    with pytest.raises(ValidationError, match=message):
        load_session(input_path)

    output_path = (
        tmp_path
        / f"version-{session_version}-schema-{schema}-output.gbdraw"
    )
    with pytest.raises(ValidationError, match=message):
        write_session_file(output_path, session)
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


def test_gui_session_replays_arrow_shape_and_geometry_flags(tmp_path: Path) -> None:
    session = _minimal_session({"c_gb": _file_entry("input.gb", b"LOCUS       A\n")})
    session["config"]["adv"].update(
        {
            "feature_shapes": {"CDS": "arrow"},
            "arrow_head_length_ratio": "1.25",
            "arrow_shaft_width_ratio": "0.25",
        }
    )

    spec = session_to_cli_args(
        session,
        mode="circular",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    assert "CDS=arrow" in session_io_module._option_all_values(
        spec.args,
        "--feature_shape",
    )
    head_index = spec.args.index("--arrow_head_length_ratio")
    assert spec.args[head_index + 1] == "1.25"
    shaft_index = spec.args.index("--arrow_shaft_width_ratio")
    assert spec.args[shaft_index + 1] == "0.25"


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    ("head_ratio", "shaft_ratio"),
    [("auto", "1"), ("1.25", "0.25")],
)
def test_cli_session_keeps_arrow_options_as_cli_provenance_only(
    mode: str,
    head_ratio: str,
    shaft_ratio: str,
) -> None:
    payload = drawing_view(build_session_payload(
        SessionBuildContext(
            mode=mode,
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=(
                "--gbk",
                "input.gb",
                "--feature_shape",
                "CDS=arrow",
                "--arrow_head_length_ratio",
                head_ratio,
                "--arrow_shaft_width_ratio",
                shaft_ratio,
            ),
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 8, 1),
        canonical_request=_canonical_request(mode),
    ))

    assert "config" not in payload
    assert payload["cliInvocation"]["args"][-6:] == [
        "--feature_shape",
        "CDS=arrow",
        "--arrow_head_length_ratio",
        head_ratio,
        "--arrow_shaft_width_ratio",
        shaft_ratio,
    ]


@pytest.mark.parametrize(
    ("version", "request_schema"),
    (
        (31, 1),
        (32, 2),
        (33, 2),
        (39, CANONICAL_REQUEST_SCHEMA),
        (40, CANONICAL_REQUEST_SCHEMA),
        (41, CANONICAL_REQUEST_SCHEMA),
        (MODE_SCOPED_SESSION_MAX_VERSION, CANONICAL_REQUEST_SCHEMA),
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
    if version >= 40:
        session["results"] = []
        session["editorState"] = {"featureCatalog": None}

    validate_session(session)
    assert version in SUPPORTED_SESSION_VERSIONS


def test_current_canonical_session_schema_is_supported() -> None:
    session = build_session_payload(
        SessionBuildContext(mode="circular", output_prefix="out", render_formats=("svg",)),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={},
        generated_at=datetime(2026, 7, 30),
        canonical_request=_canonical_request("circular"),
    )

    validate_session(session)
    assert session["version"] == CURRENT_SESSION_VERSION in SUPPORTED_SESSION_VERSIONS
    assert drawing_view(session)["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA


def _current_session(embedded_files: dict | None = None) -> dict:
    return build_session_payload(
        SessionBuildContext(mode="circular", output_prefix="out", render_formats=("svg",)),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files=embedded_files or {},
        generated_at=datetime(2026, 7, 30),
        canonical_request=_canonical_request("circular"),
    )


def test_a_session_file_is_a_container_that_loads_back(tmp_path: Path) -> None:
    session = _current_session()
    session_path = tmp_path / "session.gbdraw"

    write_session_file(session_path, session)

    assert session_path.read_bytes().startswith(b"PK\x03\x04")
    assert load_session(session_path) == json_document(session)


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
    session = _minimal_session(legacy_file_tree(files), mode="linear")
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

    assert spec.args[spec.args.index("--losat") + 1] == "losatp"
    assert spec.args[spec.args.index("--losatp_mode") + 1] == "collinear"
    assert spec.args[spec.args.index("--losat_threads") + 1] == "2"
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


def _windows_cli_args(tmp_path: Path, windows: list[dict[str, object]]) -> list[str]:
    session = _minimal_session(
        {
            "linearSeqs": [
                {
                    "gb": _file_entry("a.gb", b"LOCUS       A\n"),
                    "region_record_id": "RecA",
                    "region_start": None,
                    "region_end": None,
                    "region_reverse": False,
                    "windows": windows,
                },
                {"gb": _file_entry("b.gb", b"LOCUS       B\n")},
            ]
        },
        mode="linear",
    )
    return session_to_cli_args(
        session, mode="linear", temp_dir=tmp_path, output_override=None, format_override=None
    ).args


def test_gui_only_linear_session_writes_each_window_of_a_card_as_a_region(tmp_path: Path) -> None:
    """A card that draws a record in several windows (``windows``, D-01) gives one ``--region`` per window, in order."""

    args = _windows_cli_args(tmp_path, [
        {"start": 1001, "end": 3000, "reverseComplement": False, "column": 2},
        {"start": 6001, "end": 8000, "reverseComplement": True, "column": 0},
    ])
    regions = [args[index + 1] for index, token in enumerate(args) if token == "--region"]
    assert regions == ["#1:1001-3000", "#1:6001-8000:rc"]
    assert "--reverse_complement" not in args


@pytest.mark.parametrize(
    ("windows", "regions"),
    [
        # One window across the origin is the card's only window (its crop fields cannot hold it).
        ([{"start": 16001, "end": 500, "reverseComplement": False, "column": None}], ["#1:16001-500:origin"]),
        (
            [
                {"start": 16001, "end": 500, "reverseComplement": True, "column": None},
                {"start": 6001, "end": 8000, "reverseComplement": False, "column": None},
            ],
            ["#1:16001-500:origin:rc", "#1:6001-8000"],
        ),
    ],
    ids=["one-window", "two-windows"],
)
def test_gui_only_linear_session_writes_a_window_across_the_origin_with_origin(
    tmp_path: Path, windows: list[dict[str, object]], regions: list[str]
) -> None:
    """``:origin`` spells a window across the origin on the command line (PD-OI-122); it reads back as that window."""

    args = _windows_cli_args(tmp_path, windows)
    written = [args[index + 1] for index, token in enumerate(args) if token == "--region"]
    assert written == regions
    assert [(spec.start, spec.end, spec.reverse_complement) for spec in map(parse_region_spec, written)] == [
        (window["start"], window["end"], window["reverseComplement"]) for window in windows
    ]


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

    assert spec.args[spec.args.index("--losat") + 1] == "losatp"
    assert spec.args[spec.args.index("--losatp_mode") + 1] == "similarity_groups"
    assert (
        spec.args[spec.args.index("--similarity_alignment_feature") + 1]
        == "target_feature"
    )


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

    assert spec.args[spec.args.index("--losat") + 1] == "losatp"
    assert spec.args[spec.args.index("--losatp_mode") + 1] == "similarity_groups"
    assert spec.args[spec.args.index("--losat_threads") + 1] == "4"


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_session_keeps_hidden_scale_in_canonical_request(
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
    payload = drawing_view(build_session_payload(
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
    ))

    assert "config" not in payload
    assert payload["cliInvocation"]["args"] == args
    assert (
        payload["renderRequest"]["diagramOptions"]["config"]["objects"]["scale"][
            "show"
        ]
        is False
    )
def test_cli_session_keeps_lossless_cli_provenance_out_of_web_config() -> None:
    args = (
        "-f",
        "interactive_svg",
        "--gbk",
        "AP027078.gb",
        "AP027131.gb",
        "AP027133.gb",
        "AP027132.gb",
        "NZ_CP006932.gb",
        "--losat",
        "losatp",
        "--losatp_mode",
        "similarity_groups",
        "--losat_threads",
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

    payload = drawing_view(build_session_payload(
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
    ))

    assert payload["cliInvocation"]["renderFormats"] == ["interactive_svg"]
    assert payload["cliInvocation"]["args"] == list(args)
    assert "config" not in payload


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
        "--losat",
        "losatp",
        "--losatp_mode",
        "similarity_groups",
        "--similarity_alignment_feature",
        "target_feature",
    )
    payload = drawing_view(build_session_payload(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=args,
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
    ))

    assert payload["version"] == CURRENT_SESSION_VERSION
    assert "files" not in payload
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert payload["resources"]
    assert payload["orthogroupState"] == {}
    assert "config" not in payload
    assert payload["cliInvocation"]["args"] == list(args)


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
    payload = drawing_view(build_session_payload(
        SessionBuildContext(
            mode="linear",
            output_prefix="out",
            render_formats=("svg",),
            cli_invocation_args=args,
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={
            "linearSeqs": [
                {"definition": "", "region_record_id": "", "region_start": None, "region_end": None, "region_reverse": False}
            ]
        },
        generated_at=datetime(2026, 6, 23),
        canonical_request=_canonical_request("linear"),
    ))

    assert payload["version"] == CURRENT_SESSION_VERSION
    assert "files" not in payload
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert payload["resources"]


def test_session_pre_parse_rejects_unsupported_options() -> None:
    with pytest.raises(SystemExit):
        parse_session_pre_args(
            ["--session", "session.gbdraw", "--gbk", "input.gb"],
            mode="circular",
        )


def test_session_cli_help_uses_underscore_options(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        circular_cli_module._get_args(["--help"])

    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    assert "--save_session" in help_text
    assert "--session_output" in help_text
    assert ".gbdraw" in help_text
    assert "--save-session" not in help_text
    assert "--session-output" not in help_text


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    "option_args",
    [
        ["--save-session"],
        ["--session-output", "roundtrip.gbdraw"],
    ],
)
def test_session_preparser_rejects_removed_hyphen_aliases(
    mode: str,
    option_args: list[str],
) -> None:
    with pytest.raises(SystemExit):
        parse_session_pre_args(
            ["--session", "session.gbdraw", *option_args],
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
            "--session_output=roundtrip.gbdraw",
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
    ) == tmp_path / "prefix.gbdraw"
    assert resolve_session_sidecar_path(
        explicit_path=None,
        output_prefix=None,
        outputs=[output],
    ) == tmp_path / "diagram.gbdraw"


def test_render_formats_do_not_include_session_json_aliases() -> None:
    assert "json" not in ACCEPTED_FORMATS
    assert "session-json" not in ACCEPTED_FORMATS


def test_circular_cli_save_session_round_trip(tmp_path: Path, examples_dir: Path) -> None:
    output_prefix = tmp_path / "circular_session"
    session_path = tmp_path / "custom-session.gbdraw"
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

    load_session(session_path)
    document = read_session_payload(session_path)
    assert document["format"] == SESSION_FORMAT
    assert document["version"] == CURRENT_SESSION_VERSION
    payload = drawing_view(document)
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert (
        "outputPrefix"
        not in payload["renderRequest"]["diagramOptions"].get("output", {})
    )
    assert "files" not in payload
    assert payload["resources"]["record-1-genbank"]["data"]
    assert "<svg" in payload["results"][0]["content"]
    catalog = payload["editorState"]["featureCatalog"]
    assert catalog["schema"] == CURRENT_FEATURE_CATALOG_SCHEMA
    assert len(catalog["items"]) == len(payload["results"]) == 1
    assert catalog["items"][0]["features"]
    assert catalog["items"][0]["biologicalFeatures"]
    assert "features" not in payload and payload["modes"] == {"circular": {}}
    assert "groups" not in payload["orthogroupState"]
    assert payload["cliInvocation"]["mode"] == "circular"
    assert payload["cliInvocation"]["fileBindings"][0]["slot"] == "files.circularSources[0].gb"
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


_RUN_METADATA_FIELDS = ("trackSlotGeometry", "annotationWarnings", "featureIdentityNotices", "comparisonWarnings")


def _run_notice_inputs(tmp_path: Path) -> dict[str, Path]:
    """Tables that make a render report one run notice of each kind."""

    tables = {
        # A region annotation whose feature the record lacks: an annotation warning.
        "annotations": "set_id\tid\tmark\trecord\tfeature_selector\tlabel\n"
        "s\tmissing\tband\t#1\tgene=PRIVATE-MISSING\tMissing\n",
        # An edit of a lambda gene outside the crop: a feature identity notice.
        "overrides": "record\tfeature_selector\tlabel_visibility\tlabel_text\n"
        "NC_001416.1\tprotein_id=NP_040580.1\ton\tFirst gene\n",
        # A comparison row whose sequence IDs name no drawn record: a comparison warning.
        "blast": "contig_A\tcontig_B\t99.0\t1000\t0\t0\t2001\t3000\t1\t1000\t1e-100\t1800.0\n",
    }
    paths = {name: tmp_path / f"{name}.tsv" for name in tables}
    for name, text in tables.items():
        paths[name].write_text(text)
    return paths


@pytest.mark.parametrize("case", ["circular", "circular-batch", "linear"])
def test_cli_session_saves_the_run_metadata_the_worker_reports(tmp_path: Path, case: str) -> None:
    """OV-3240, OV-3244: a CLI-saved Session holds each record's axis and the
    run notices, so the Web can select a region on the loaded Result and shows
    the notices before any Generate. They are what the Web Worker reports for
    the same request (one producer, ``gbdraw.api.run_metadata``), which a Web
    Save writes back as it is (``session-cli-compatibility.test.mjs``), and
    each row names its Result as the Session names it, also in a Circular
    batch whose second Result has its own rows."""

    from gbdraw.web_support.request_render import render_embedded_canonical_web_request

    root = Path(__file__).resolve().parents[1]
    lambda_gb = root / "tests/test_inputs/NC_001416.gb"
    mito_gb = root / "tests/fixtures/sessions/cli-web-mito.gb"
    tables = _run_notice_inputs(tmp_path)
    session_path = tmp_path / "cli.gbdraw"
    common = ["-o", str(tmp_path / "cli"), "-f", "svg", "--session_output", str(session_path),
              "--annotation_table", str(tables["annotations"])]
    if case == "circular":
        circular_main(["--gbk", str(lambda_gb), *common])
    elif case == "circular-batch":
        circular_main(["--gbk", str(lambda_gb), str(mito_gb), *common])
    else:
        linear_main([
            "--gbk", str(lambda_gb), str(mito_gb), *common,
            "-b", str(tables["blast"]), "--feature_override_table", str(tables["overrides"]),
            "--region", "NC_001416.1:20000-30000",
        ])
    # The Worker reads the resources' bytes, which a Session file holds beside the document.
    payload = drawing_view(read_session_payload(session_path))
    reply = render_embedded_canonical_web_request(
        payload["renderRequest"], resources=payload["resources"], workspace=tmp_path / "web"
    )["metadata"]
    names = [result["name"] for result in payload["results"]]
    assert len(names) == (2 if case == "circular-batch" else 1)

    def rows(metadata: dict, field: str) -> list[dict]:
        return metadata[field]["records"] if field == "trackSlotGeometry" else metadata[field]

    # The Worker names a Result by its file name, a Session by its output
    # prefix; every other byte of a row is the Worker's.
    reported = copy.deepcopy({field: reply[field] for field in _RUN_METADATA_FIELDS if field in reply})
    assert set(_RUN_METADATA_FIELDS) - {"trackSlotGeometry"} <= set(reported)
    for field in reported:
        for row in rows(reported, field):
            if "resultName" in row:
                row["resultName"] = names[row["resultIndex"]]

    saved = payload["runMetadata"]
    assert saved == reported
    assert list(saved) == list(_RUN_METADATA_FIELDS)
    assert all(saved[field] for field in ("annotationWarnings", *(
        ("featureIdentityNotices", "comparisonWarnings") if case == "linear" else ())))
    assert saved["trackSlotGeometry"]["mode"] == case.split("-")[0]
    assert {row["resultIndex"] for row in rows(saved, "trackSlotGeometry")} == set(range(len(names)))
    # Each row names its Result as the Session does, as in a Web Save (OV-3243).
    for field in saved:
        for row in rows(saved, field):
            assert row.get("resultName", names[row["resultIndex"]]) == names[row["resultIndex"]]


def test_cli_resave_of_a_session_30_replaces_its_run_metadata(tmp_path: Path) -> None:
    """A re-save of an older Session (27-30, migrated on the way) keeps the
    new run's metadata as it is, in place of the source's, never beside it."""

    source = json.loads(gzip.decompress(
        (Path(__file__).parent / "fixtures/sessions/scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz").read_bytes()
    ))
    assert source["version"] == 30
    source["runMetadata"] = {
        "trackSlotGeometry": {"schema": 1, "mode": "circular", "source": "resolved",
                              "records": [{"resultIndex": 0, "resultName": "stale", "recordIndex": 0}]},
        "annotationWarnings": [{
            "code": "feature_selector_unmatched", "setId": "s", "annotationId": "stale", "recordId": "r",
            "recordIndex": 0, "missingCount": 1, "message": "stale", "resultIndex": 0, "resultName": "stale",
        }],
    }
    validate_session(source)
    geometry = {"schema": 1, "mode": "circular", "source": "resolved",
                "records": [{"resultIndex": 0, "resultName": "again", "recordIndex": 0}]}

    def resave(run_metadata: dict | None) -> dict:
        document = build_session_payload(
            SessionBuildContext(mode="circular", output_prefix="again", render_formats=("svg",),
                                source_session=source, run_metadata=run_metadata),
            svg_results=(("again", "<svg></svg>"),),
            embedded_files=source["files"],
            generated_at=datetime(2026, 11, 10),
            canonical_request=_canonical_request("circular"),
        )
        assert "stale" not in json.dumps(document)
        return drawing_view(document, mode="circular")

    empty_notices = {"annotationWarnings": [], "featureIdentityNotices": [], "comparisonWarnings": []}
    fresh = {"trackSlotGeometry": geometry, **empty_notices}
    assert resave(fresh)["runMetadata"] == fresh
    assert resave(empty_notices)["runMetadata"] == empty_notices
    # A save with no run metadata leaves none of the source behind.
    assert "runMetadata" not in resave(None)


def _circular_genbank(*records: tuple[str, str, str]) -> str:
    """Small circular GenBank records: (record ID, gene, CDS location)."""

    sequence = "".join(
        f"{index * 60 + 1:>9} {'acgt' * 15}\n" for index in range(6)
    )
    return "".join(
        f"LOCUS       {record_id:<24} 360 bp    DNA     circular UNA 01-JAN-2000\n"
        "DEFINITION  CLI Session input.\n"
        f"ACCESSION   {record_id}\nVERSION     {record_id}\nKEYWORDS    .\n"
        "SOURCE      synthetic construct\n  ORGANISM  synthetic construct\n            .\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..360\n"
        f"     CDS             {location}\n"
        f'                     /gene="{gene}"\n'
        f'                     /locus_tag="{gene}"\n'
        '                     /translation="MKKKKKKKKK"\n'
        f"ORIGIN\n{sequence}//\n"
        for record_id, gene, location in records
    )


def _cli_session_inputs(tmp_path: Path) -> tuple[Path, Path]:
    single = tmp_path / "single.gbk"
    single.write_text(_circular_genbank(("cli_single", "single_gene", "41..125")))
    multi = tmp_path / "multi.gbk"
    multi.write_text(_circular_genbank(
        ("cli_multi_a", "multi_a_gene", "61..150"),
        ("cli_multi_b", "multi_b_gene", "complement(101..200)"),
    ))
    return single, multi


def _resource_bytes(session: dict, resource_id: str) -> bytes:
    return base64.b64decode(session["resources"][resource_id]["data"])


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_session_binds_each_unmodified_input_file_to_its_records(
    tmp_path: Path, mode: str
) -> None:
    """Web Load sees the request draw its records from the bound input files.

    The CLI binds each input file for the Web draft. When the request draws a
    file's records unchanged, it reads the same resource, which holds the
    file's own bytes, so record actions need no Generate after Load.
    """

    single, multi = _cli_session_inputs(tmp_path)
    session_path = tmp_path / f"{mode}.gbdraw"
    prefix = tmp_path / mode
    main = circular_main if mode == "circular" else linear_main
    inputs = [str(single)] if mode == "circular" else [str(single), str(multi)]
    main(["--gbk", *inputs, "-o", str(prefix), "-f", "svg",
          "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    bindings = session["webFiles"]["bindings"]
    records = session["renderRequest"]["records"]
    if mode == "circular":
        (card,) = bindings["circularSources"]
        assert card["uid"] == "record-1"
        assert [record["source"] for record in records] == [
            {"kind": "genbank", "resourceId": card["gb"]["resourceId"]}
        ]
        assert [record["selector"] for record in records] == [None]
        files = {card["gb"]["resourceId"]: single}
    else:
        first, second = (sequence["gb"]["resourceId"] for sequence in bindings["linearSeqs"])
        assert [(record["source"]["resourceId"], record["selector"]) for record in records] == [
            (first, None),
            (second, {"kind": "recordId", "value": "cli_multi_a"}),
            (second, {"kind": "recordId", "value": "cli_multi_b"}),
        ]
        files = {first: single, second: multi}
    assert sorted(session["resources"]) == sorted(files)
    for resource_id, path in files.items():
        assert _resource_bytes(session, resource_id) == path.read_bytes()
        assert session["resources"][resource_id]["name"] == path.name

    replay = tmp_path / f"{mode}-replay"
    main(["--session", str(session_path), "-o", str(replay), "-f", "svg"])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


def test_cli_session_binds_unmodified_gff_and_fasta_inputs_to_their_record(
    tmp_path: Path,
) -> None:
    inputs = Path(__file__).parent / "test_inputs"
    session_path = tmp_path / "gff.gbdraw"
    circular_main(["--gff", str(inputs / "NC_013668.gff3"), "--fasta",
                   str(inputs / "NC_013668.fasta"), "-o", str(tmp_path / "gff"),
                   "-f", "svg", "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    (card,) = session["webFiles"]["bindings"]["circularSources"]
    assert (card["uid"], card["gb"]) == ("record-1", None)
    assert [record["source"] for record in session["renderRequest"]["records"]] == [{
        "kind": "gffFasta",
        "gffResourceId": card["gff"]["resourceId"],
        "fastaResourceId": card["fasta"]["resourceId"],
    }]
    assert sorted(session["resources"]) == sorted(
        [card["gff"]["resourceId"], card["fasta"]["resourceId"]]
    )
    assert _resource_bytes(session, card["gff"]["resourceId"]) == (
        inputs / "NC_013668.gff3"
    ).read_bytes()


def _write_records_table(tmp_path: Path) -> Path:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\tregion\treverse_complement\n"
        "single.gbk\t\t21-300:rc\t0\n"
        "multi.gbk\tcli_multi_a\t\t1\n"
        "multi.gbk\tcli_multi_b\t21-300\t1\n",
        encoding="utf-8",
    )
    return table


@pytest.mark.parametrize("case", ["options", "table"])
def test_linear_cli_session_records_crop_and_orientation_on_the_unchanged_input_file(
    tmp_path: Path, case: str
) -> None:
    """A cropped or reversed record reads its input file with a region and orientation.

    The request names the crop (in source coordinates) and the orientation, as
    the Web writes them, instead of a re-serialized copy of the drawn record.
    The Web draft that loads the bound input file then draws the same record,
    and a replay draws the same SVG. ``--reverse_complement`` applies after
    the crop, as ``:rc`` does (OV-1279), so the request names the same crop.
    """

    single, multi = _cli_session_inputs(tmp_path)
    session_path = tmp_path / f"{case}.gbdraw"
    prefix = tmp_path / case
    if case == "options":
        args = ["--gbk", str(single), str(multi),
                "--region", "cli_single:21-300", "--region", "cli_multi_b:21-300",
                "--reverse_complement", "0", "--reverse_complement", "1"]
        single_region = [{"start": 21, "end": 300, "reverseComplement": False, "column": None}]
    else:
        args = ["--records_table", str(_write_records_table(tmp_path))]
        single_region = [{"start": 21, "end": 300, "reverseComplement": True, "column": None}]
    linear_main([*args, "-o", str(prefix), "-f", "svg", "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    records = session["renderRequest"]["records"]
    first, second = records[0]["source"]["resourceId"], records[1]["source"]["resourceId"]
    assert [
        (record["source"]["resourceId"], record["selector"], record["windows"],
         record["presentation"]["reverseComplement"])
        for record in records
    ] == [
        (first, None, single_region, False),
        (second, {"kind": "recordId", "value": "cli_multi_a"}, [], True),
        # Cropped to source 21..300, then reversed.
        (second, {"kind": "recordId", "value": "cli_multi_b"},
         [{"start": 21, "end": 300, "reverseComplement": True, "column": None}], False),
    ]
    assert _resource_bytes(session, first) == single.read_bytes()
    assert _resource_bytes(session, second) == multi.read_bytes()
    if case == "options":
        bound = [sequence["gb"]["resourceId"]
                 for sequence in session["webFiles"]["bindings"]["linearSeqs"]]
        assert bound == [first, second]

    replay = tmp_path / f"{case}-replay"
    linear_main(["--session", str(session_path), "-o", str(replay), "-f", "svg"])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


@pytest.mark.parametrize("case", ["options", "table"])
def test_linear_cli_session_selects_the_drawn_records_of_a_partly_drawn_input_file(
    tmp_path: Path, case: str
) -> None:
    """A file of which the CLI draws some records reads the file and selects them.

    The request reads the input file's own bytes and selects each drawn record
    by its ID, or by its index when another record of the file has the same
    ID, as the Web writes the records it draws of one multi-record File. The
    bound File is that same resource, so Web Load draws only those records,
    and a replay draws the same SVG.
    """

    single, multi = _cli_session_inputs(tmp_path)
    session_path = tmp_path / f"{case}.gbdraw"
    prefix = tmp_path / case
    if case == "options":
        args = ["--gbk", str(single), str(multi), "--record_id", "", "--record_id", "cli_multi_b"]
        selected = [{"kind": "recordId", "value": "cli_multi_b"}]
    else:
        multi.write_text(_circular_genbank(
            ("twin", "twin_a_gene", "61..150"),
            ("cli_unique", "unique_gene", "101..200"),
            ("twin", "twin_c_gene", "201..300"),
        ))
        table = tmp_path / "records.tsv"
        table.write_text(
            "gbk\trecord_id\nsingle.gbk\t\nmulti.gbk\t#3\nmulti.gbk\tcli_unique\n",
            encoding="utf-8",
        )
        args = ["--records_table", str(table)]
        # Only one drawn record is "twin", but the file has two.
        selected = [{"kind": "recordIndex", "index": 2}, {"kind": "recordId", "value": "cli_unique"}]
    linear_main([*args, "-o", str(prefix), "-f", "svg", "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    records = session["renderRequest"]["records"]
    first, second = records[0]["source"]["resourceId"], records[1]["source"]["resourceId"]
    assert [(record["source"]["resourceId"], record["selector"], record["windows"])
            for record in records] == [
        (first, None, []), *((second, selector, []) for selector in selected)
    ]
    assert sorted(session["resources"]) == sorted({first, second})
    assert _resource_bytes(session, first) == single.read_bytes()
    assert _resource_bytes(session, second) == multi.read_bytes()
    if case == "options":
        bound = [sequence["gb"]["resourceId"]
                 for sequence in session["webFiles"]["bindings"]["linearSeqs"]]
        assert bound == [first, second]

    replay = tmp_path / f"{case}-replay"
    linear_main(["--session", str(session_path), "-o", str(replay), "-f", "svg"])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


def _varied_cli_session_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Circular GenBank files whose GC changes along each record, so a window that reads a crop's flank shows it."""

    import random

    from Bio import SeqIO

    def record(record_id: str, seed: int) -> SeqRecord:
        rng = random.Random(seed)
        sequence = "".join(rng.choice(bases) for bases in ("GC", "AT", "GC") for _ in range(1_000))
        made = SeqRecord(Seq(sequence), id=record_id, name=record_id, description="CLI Session input.",
                         annotations={"topology": "circular", "molecule_type": "DNA"})
        made.features = [SeqFeature(FeatureLocation(1_200, 1_500, 1), type="CDS",
                                    qualifiers={"locus_tag": [f"{record_id}_gene"]})]
        return made

    single = tmp_path / "single.gbk"
    SeqIO.write([record("cli_single", 1)], single, "genbank")
    multi = tmp_path / "multi.gbk"
    SeqIO.write([record("cli_multi_a", 2), record("cli_multi_b", 3)], multi, "genbank")
    return single, multi


def _replayed_cli_svgs(session_path: Path, tmp_path: Path, prefix: str) -> list[tuple[bytes, bytes]]:
    """Each SVG the command drew under ``prefix`` beside the SVG of the same name that its Session's replay drew."""

    replay = tmp_path / f"{prefix}-replay"
    circular_main(["--session", str(session_path), "-o", str(replay), "-f", "svg"])
    replayed = sorted(tmp_path.glob(f"{prefix}-replay*.svg"))
    assert replayed
    return [
        ((tmp_path / path.name.replace(f"{prefix}-replay", prefix, 1)).read_bytes(), path.read_bytes())
        for path in replayed
    ]


@pytest.mark.parametrize("canvas", [False, True])
def test_circular_cli_grid_and_batch_sessions_save_crops_as_source_and_range(
    tmp_path: Path, canvas: bool
) -> None:
    """A records-table crop of a Circular grid or batch is saved as its source record and range (PD-OI-044 rev 5).

    The replay crops the source again, so the GC windows at a crop's ends read the bases beside it as the command
    did. A crop reverse-complemented after it is cut is a crop too; a whole record reverse-complemented without a
    crop keeps its drawn copy, which the Web draws as is.
    """

    single, multi = _varied_cli_session_inputs(tmp_path)
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\tregion\treverse_complement\n"
        "single.gbk\t\t201-1500:rc\t0\n"
        "multi.gbk\tcli_multi_a\t\t1\n"
        "multi.gbk\tcli_multi_b\t301-1800\t1\n",
        encoding="utf-8",
    )
    session_path = tmp_path / "circular.gbdraw"
    circular_main(["--records_table", str(table), *(["--multi_record_canvas"] if canvas else []),
                   "-o", str(tmp_path / "circular"), "-f", "svg", "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    records = session["renderRequest"]["records"]
    assert session["renderRequest"]["grouping"] == ("grid" if canvas else "batch")
    assert [(record["windows"], record["presentation"]["reverseComplement"]) for record in records] == [
        ([{"start": 201, "end": 1500, "reverseComplement": True, "column": None}], False),
        ([], False),
        ([{"start": 301, "end": 1800, "reverseComplement": True, "column": None}], False),
    ]
    # The crop of the single-record file reads the file itself; the multi-record file's whole reversed record is
    # its drawn copy, so its cropped record reads the record it was cut from.
    assert _resource_bytes(session, records[0]["source"]["resourceId"]) == single.read_bytes()
    assert _resource_bytes(session, records[1]["source"]["resourceId"]) != multi.read_bytes()
    assert b"cli_multi_b" in _resource_bytes(session, records[2]["source"]["resourceId"])
    replayed = _replayed_cli_svgs(session_path, tmp_path, "circular")
    assert len(replayed) == (1 if canvas else 3)
    feature_hash = re.compile(rb"\bf[0-9a-f]{8}\b")
    for live, replay in replayed:
        # The whole reversed record's copy is the replay's source, so its features take the copy's hashes (OV-415).
        assert feature_hash.sub(b"f", replay) == feature_hash.sub(b"f", live)


def test_circular_cli_one_record_session_is_a_single_drawing(tmp_path: Path) -> None:
    """One record drawn without --multi_record_canvas is saved as the Web saves a one-record drawing (OV-3297)."""

    single, _multi = _varied_cli_session_inputs(tmp_path)
    session_path = tmp_path / "one.gbdraw"
    circular_main(["--gbk", str(single), "-o", str(tmp_path / "one"), "-f", "svg",
                   "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    request = session["renderRequest"]
    [record] = request["records"]
    assert (request["grouping"], request["output"]["prefix"]) == ("single", "one")
    assert (record["selector"], record["windows"]) == (None, [])
    assert _resource_bytes(session, record["source"]["resourceId"]) == single.read_bytes()
    assert sorted(path.name for path in tmp_path.glob("one*.svg")) == ["one.svg"]
    [(live, replay)] = _replayed_cli_svgs(session_path, tmp_path, "one")
    assert replay == live


def test_circular_cli_one_records_table_row_saves_its_crop_as_a_single_drawing(tmp_path: Path) -> None:
    """One records-table row of a two-record file, cropped and reversed, is the Web's one-record crop (OV-3297).

    The Session reads the input file itself with the row's record and range, so the Web opens it on that record
    with its crop in the single-record fields, and the replay draws the command's figure.
    """

    _single, multi = _varied_cli_session_inputs(tmp_path)
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\tregion\treverse_complement\trecord_label\n"
        "multi.gbk\tcli_multi_b\t301-1800\t1\tMyB\n",
        encoding="utf-8",
    )
    session_path = tmp_path / "one.gbdraw"
    circular_main(["--records_table", str(table), "-o", str(tmp_path / "one"), "-f", "svg",
                   "--session_output", str(session_path)])

    session = drawing_view(read_session_payload(session_path))
    request = session["renderRequest"]
    [record] = request["records"]
    assert (request["grouping"], request["output"]["prefix"]) == ("single", "one")
    assert record["selector"] == {"kind": "recordId", "value": "cli_multi_b"}
    assert record["windows"] == [{"start": 301, "end": 1800, "reverseComplement": True, "column": None}]
    assert (record["presentation"]["label"], record["presentation"]["reverseComplement"]) == ("MyB", False)
    assert _resource_bytes(session, record["source"]["resourceId"]) == multi.read_bytes()
    [(live, replay)] = _replayed_cli_svgs(session_path, tmp_path, "one")
    assert replay == live


@pytest.mark.parametrize(
    ("grouping", "table_rows", "extra_args"),
    [
        ("grid", "single.gbk\t\n", ["--multi_record_canvas"]),
        ("batch", "multi.gbk\tcli_multi_a\nmulti.gbk\tcli_multi_b\n", []),
    ],
    ids=["one-record-grid", "two-record-batch"],
)
def test_circular_cli_grid_of_one_and_batch_of_two_keep_their_grouping(
    tmp_path: Path, grouping: str, table_rows: str, extra_args: list[str]
) -> None:
    """Only one record drawn as separate diagrams becomes a single drawing; a grid and a real batch stay."""

    _varied_cli_session_inputs(tmp_path)
    table = tmp_path / "records.tsv"
    table.write_text("gbk\trecord_id\n" + table_rows, encoding="utf-8")
    session_path = tmp_path / "kept.gbdraw"
    circular_main(["--records_table", str(table), *extra_args, "-o", str(tmp_path / "kept"), "-f", "svg",
                   "--session_output", str(session_path)])

    load_session(session_path)  # the saved file loads with its checks
    request = drawing_view(read_session_payload(session_path))["renderRequest"]
    assert request["grouping"] == grouping
    assert isinstance(request["output"], list) is (grouping == "batch")


def test_python_api_batch_of_one_saves_its_crop_on_the_input_file(tmp_path: Path) -> None:
    """A typed batch of one cropped record stays a batch and reads its input file with the range."""

    single, _multi = _varied_cli_session_inputs(tmp_path)
    rendered = render_request(CircularBatchRequest(
        records=(RecordInput(source=GenBankInputSource(single), region=parse_region_spec("201-1500")),),
        outputs=(RenderOutputRequest(output_prefix="one", output_directory=tmp_path),),
    ))
    session_path = tmp_path / "one.gbdraw"
    save_session_document(session_path, rendered.request)

    session = drawing_view(read_session_payload(session_path))
    [record] = session["renderRequest"]["records"]
    assert session["renderRequest"]["grouping"] == "batch"
    assert record["windows"] == [{"start": 201, "end": 1500, "reverseComplement": False, "column": None}]
    assert _resource_bytes(session, record["source"]["resourceId"]) == single.read_bytes()
    [(live, replay)] = _replayed_cli_svgs(session_path, tmp_path, "one")
    assert replay == live


def test_single_circular_session_records_crop_on_the_unchanged_input_file(
    tmp_path: Path,
) -> None:
    """A single Circular request's crop and orientation are the Web form's."""

    single, _multi = _cli_session_inputs(tmp_path)
    document = build_session_document(CircularDiagramRequest(
        records=(RecordInput(
            source=GenBankInputSource(single),
            region=parse_region_spec("21-300:rc"),
        ),),
        output=RenderOutputRequest(output_prefix="single", output_directory=tmp_path),
        grouping="single",
    ))
    document = drawing_view(embedded_session(document.to_dict(), document.resources.sources))

    [record] = document["renderRequest"]["records"]
    assert (record["selector"], record["windows"], record["presentation"]["reverseComplement"]) == (
        None, [{"start": 21, "end": 300, "reverseComplement": True, "column": None}], False
    )
    assert _resource_bytes(document, record["source"]["resourceId"]) == single.read_bytes()


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


def _composite_binding_session():
    payload = json.loads((Path(__file__).parent / 'fixtures/sessions/single.v41-bindings1.json').read_text())
    bindings = payload['webFiles']['bindings']
    leaf = bindings['c_gb']
    bindings.update(schema=2, c_gb={
        'kind': 'composite', 'components': [leaf, {**leaf, 'name': 'second.gb'}],
        'name': 'logical.gb', 'type': 'text/plain', 'lastModified': 0.5,
    })
    return payload


def test_schema1_frozen_main_writer_is_preserved(tmp_path):
    root = Path(__file__).parent / 'fixtures/sessions'
    fixture = root / 'single.v41-bindings1.json'
    provenance = json.loads((root / 'single.v41-bindings1.provenance.json').read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest() == provenance['sha256']
    load_session(fixture)
    payload = read_session_payload(fixture)
    output = tmp_path / 'preserved.json'
    write_session_file(output, payload)
    assert load_session(output)['webFiles']['bindings'] == payload['webFiles']['bindings']
    assert payload['webFiles']['bindings']['schema'] == 1


@pytest.mark.parametrize("modified", [0, 0.5])
def test_cli_writes_one_circular_file_card_per_gbk_input(tmp_path, monkeypatch, modified):
    # circular-multi-file D-01, D-08 (DESIGN section 7 test 4): one File card
    # per --gbk value, keyed by the CLI's record key of that input, even when
    # two inputs hold the same bytes; no composite.
    import gbdraw.cli_utils.session as cli_session

    first, second = tmp_path / "first.gb", tmp_path / "second.gb"
    first.write_bytes(b"same bytes\n")
    second.write_bytes(first.read_bytes())
    real_entry = session_io_module.file_resource_entry

    def unnamed_entry(path):
        entry = real_entry(path)
        descriptor = {**entry["descriptor"], "name": "", "type": "", "lastModified": int(modified)}
        return {**entry, "descriptor": descriptor, "name": "", "type": "", "lastModified": modified}

    monkeypatch.setattr(cli_session, "file_resource_entry", unnamed_entry)
    files, provenance = collect_embedded_files_from_cli_args(
        "circular", ["--gbk", str(first), str(second), str(first)]
    )
    cards = files["circularSources"]
    assert [card["uid"] for card in cards] == ["record-1", "record-2", "record-3"]
    assert not {"c_gb", "c_gff", "c_fasta"} & set(files)
    assert [binding.name for binding in provenance] == ["first.gb", "second.gb", "first.gb"]
    for index, binding in enumerate(provenance):
        assert binding.slot == f"files.circularSources[{index}].gb"
        assert session_io_module.get_session_slot({"files": files}, binding.slot) == cards[index]["gb"]
    payload = build_session_payload(
        SessionBuildContext(mode="circular", output_prefix="out", render_formats=("svg",)),
        svg_results=[], embedded_files=files, generated_at=datetime(2026, 9, 11),
        canonical_request=_canonical_request("circular"),
    )
    bindings = drawing_view(payload)["webFiles"]["bindings"]
    assert not {"c_gb", "c_gff", "c_fasta"} & set(bindings)
    actual = bindings["circularSources"]
    assert [card["uid"] for card in actual] == ["record-1", "record-2", "record-3"]
    for card in actual:
        assert (card["gff"], card["fasta"]) == (None, None)
        assert {key: card["gb"][key] for key in ("name", "type", "lastModified")} == {
            "name": "", "type": "", "lastModified": modified,
        }
        assert base64.b64decode(payload["resources"][card["gb"]["resourceId"]]["data"]) == first.read_bytes()


def test_cli_writes_one_circular_file_card_per_gff_and_fasta_pair(tmp_path):
    for name in ("a.gff", "b.gff", "a.fa", "b.fa"):
        (tmp_path / name).write_text(">x\nA\n")
    files, provenance = collect_embedded_files_from_cli_args(
        "circular",
        ["--gff", str(tmp_path / "a.gff"), str(tmp_path / "b.gff"), "--fasta", str(tmp_path / "a.fa"), str(tmp_path / "b.fa")],
    )
    assert [(card["uid"], card["gb"], card["gff"]["name"], card["fasta"]["name"]) for card in files["circularSources"]] == [
        ("record-1", None, "a.gff", "a.fa"), ("record-2", None, "b.gff", "b.fa"),
    ]
    assert [binding.slot for binding in provenance] == [
        "files.circularSources[0].gff", "files.circularSources[1].gff",
        "files.circularSources[0].fasta", "files.circularSources[1].fasta",
    ]
    assert files["cliInputs"] == []


def test_schema2_composite_load_write_preserves_order_and_metadata(tmp_path):
    payload = _composite_binding_session()
    output = tmp_path / 'composite.json'
    write_session_file(output, payload)
    assert load_session(output)['webFiles']['bindings'] == payload['webFiles']['bindings']


@pytest.mark.parametrize('invalid', [
    'schema', 'kind', 'nested', 'empty', 'singleton', 'missing', 'dangling',
    'null', 'metadata', 'mixed', 'slot', 'encoding', 'unknown', 'leaf-field',
    'leaf-metadata', 'old-version',
])
def test_schema2_composite_rejects_malformed_bindings(invalid):
    payload = _composite_binding_session()
    bindings = payload['webFiles']['bindings']
    composite = bindings['c_gb']
    if invalid == 'schema': bindings['schema'] = 99
    elif invalid == 'kind': composite['kind'] = 'multipart'
    elif invalid == 'nested': composite['components'][0] = copy.deepcopy(composite)
    elif invalid == 'empty': composite['components'] = []
    elif invalid == 'singleton': composite['components'] = composite['components'][:1]
    elif invalid == 'missing': del bindings['c_gb']
    elif invalid == 'dangling': composite['components'][0]['resourceId'] = 'missing'
    elif invalid == 'null': composite['components'][0] = None
    elif invalid == 'metadata': composite['lastModified'] = -1
    elif invalid == 'mixed': composite['resourceId'] = composite['components'][0]['resourceId']
    elif invalid == 'slot': bindings['c_fasta'] = copy.deepcopy(composite)
    elif invalid == 'encoding': next(iter(payload['resources'].values()))['encoding'] = 'utf8'
    elif invalid == 'unknown': composite['size'] = 0
    elif invalid == 'leaf-field': composite['components'][0]['extra'] = True
    elif invalid == 'leaf-metadata': del composite['components'][0]['type']
    elif invalid == 'old-version': payload['version'] = 40
    with pytest.raises(ValidationError):
        session_io_module.validate_session(payload)






def test_composite_sidecar_binds_a_new_repeated_component_once_without_reencoding(monkeypatch):
    source = _composite_binding_session()
    load_session_document(source)
    document, sources = session_io_module.split_resource_bytes(source)
    payload = {'resources': {}, 'webFiles': {}}
    reads = []

    class CountedSource:
        def __init__(self, inner):
            self.inner = inner

        def open_bytes(self):
            reads.append(1)
            return self.inner.open_bytes()

    def unexpected_encode(_data):
        raise AssertionError('Existing base64 payload must not be re-encoded')
    monkeypatch.setattr(session_io_module.base64, 'b64encode', unexpected_encode)
    counted = {resource_id: CountedSource(source) for resource_id, source in sources.items()}
    inventory = session_io_module._project_web_file_binding(
        SessionResources(document['resources'], counted), document['webFiles']['bindings']['c_gb'], schema=2,
    )
    payload_sources = {}
    session_io_module._attach_current_web_file_bindings(payload, {'c_gb': inventory}, payload_sources)
    # Binding carries each source; its bytes are read once, when the Session file is written.
    assert reads == []
    assert list(payload_sources.values()) == [counted[next(iter(counted))]]
    assert len(payload['resources']) == 1
    assert len(payload['webFiles']['bindings']['c_gb']['components']) == 2


def _cli_composite_source():
    """Use the candidate writer with six draft sources distinct from the request."""
    import io
    from Bio import SeqIO

    parts = []
    for index in range(6):
        record = SeqRecord(
            Seq('ATGC' * (index + 2)), id=f'draft-{index}',
            annotations={'molecule_type': 'DNA'},
        )
        buffer = io.StringIO()
        SeqIO.write(record, buffer, 'genbank')
        parts.append({
            **_file_entry(f'part-{index}.gb', buffer.getvalue().encode()),
            'encoding': 'base64', 'lastModified': index + 0.5,
        })
    return build_mode_scoped_session_json(
        SessionBuildContext(mode='circular', output_prefix='source', render_formats=('svg',)),
        svg_results=[], generated_at=datetime(2026, 9, 9),
        embedded_files={'c_gb': {
            'kind': 'composite', 'name': 'six-originals.gb', 'type': '',
            'lastModified': 0.25, 'components': parts,
        }},
        canonical_request=_canonical_request('circular'),
    )


def _draft_identity(session, value):
    if value is None:
        return None
    if isinstance(value, list):
        return [_draft_identity(session, item) for item in value]
    if value.get('kind') == 'composite':
        return {**value, 'components': [_draft_identity(session, part) for part in value['components']]}
    return {
        **{key: val for key, val in value.items() if key != 'resourceId'},
        'bytes': base64.b64decode(session['resources'][value['resourceId']]['data']),
    }


def _card_files(session, cards):
    """The drafted GenBank file of each File card, in card order."""
    return [_draft_identity(session, card['gb']) for card in cards]


def _replay_cli_sidecar(source, tmp_path):
    source_path = tmp_path / ('source.gbdraw' if source['version'] == CURRENT_SESSION_VERSION else 'source.json')
    write_session_file(source_path, source)
    original = source_path.read_bytes()
    prefix = tmp_path / 'replayed'
    sidecar = tmp_path / 'replayed.gbdraw'
    circular_main([
        '--session', str(source_path), '-o', str(prefix),
        '--session_output', str(sidecar),
    ])
    assert source_path.read_bytes() == original
    load_session_document(sidecar)
    written = read_session_payload(sidecar)
    result = drawing_view(written)
    assert (written['version'], result['webFiles']['bindings']['schema'], result['renderRequest']['schema']) == (CURRENT_SESSION_VERSION, 2, CANONICAL_REQUEST_SCHEMA)
    assert result['renderRequest']['output']['prefix'] == 'replayed'
    assert len(result['renderRequest']['records']) == 1  # Replay consumes committed input.
    assert result['results'][0]['content'] == prefix.with_suffix('.svg').read_text()
    assert result['results'] != source['results']
    catalog = result['editorState']['featureCatalog']
    assert catalog['schema'] == CURRENT_FEATURE_CATALOG_SCHEMA and len(catalog['items']) == 1
    assert catalog['items'][0]['resultName'] == result['results'][0]['name']
    return result


@pytest.mark.parametrize('collision', [False, True])
def test_cli_i24ab_composite_sidecar_preserves_six_sources(tmp_path, collision):
    # A Session 50 whose Circular input is a composite replays to a Session 51
    # with one File card per component, in order (circular-multi-file D-03).
    source = _cli_composite_source()
    binding = source['webFiles']['bindings']['c_gb']
    if collision:
        # The replay encoder will claim this preferred draft ID for committed bytes.
        committed_id = source['renderRequest']['records'][0]['source']['resourceId']
        source['resources']['committed-input'] = source['resources'].pop(committed_id)
        source['renderRequest']['records'][0]['source']['resourceId'] = 'committed-input'
        part = binding['components'][0]
        source['resources'][committed_id] = source['resources'].pop(part['resourceId'])
        part['resourceId'] = committed_id
    source['resources'] = dict(reversed(list(source['resources'].items())))
    expected = [_draft_identity(source, part) for part in binding['components']]
    result = _replay_cli_sidecar(source, tmp_path)
    bindings = result['webFiles']['bindings']
    assert 'c_gb' not in bindings
    cards = bindings['circularSources']
    # The committed records name one composite file, so the cards take uids no flat key has (D-12, OV-2843).
    assert [card['uid'] for card in cards] == [f'source-{n}' for n in range(1, 7)]
    assert all((card['gff'], card['fasta']) == (None, None) for card in cards)
    assert _card_files(result, cards) == expected
    if collision:
        committed_id = result['renderRequest']['records'][0]['source']['resourceId']
        assert cards[0]['gb']['resourceId'] != committed_id
        assert base64.b64decode(result['resources'][committed_id]['data']) != expected[0]['bytes']


def test_cli_i24a_repeated_component_occurrences_survive_dedup(tmp_path):
    source = _cli_composite_source()
    parts = source['webFiles']['bindings']['c_gb']['components']
    parts[4] = {**parts[0], 'name': 'repeated-with-distinct-name.gb'}
    expected = [_draft_identity(source, part) for part in parts]
    result = _replay_cli_sidecar(source, tmp_path)
    cards = result['webFiles']['bindings']['circularSources']
    assert _card_files(result, cards) == expected
    assert cards[0]['gb']['resourceId'] == cards[4]['gb']['resourceId']


def test_cli_i24c_frozen_schema1_replay_emits_current_bindings(tmp_path):
    source = read_session_payload(Path(__file__).parent / 'fixtures/sessions/single.v41-bindings1.json')
    result = _replay_cli_sidecar(source, tmp_path)
    (card,) = result['webFiles']['bindings']['circularSources']
    assert card['uid'] == 'record-1'
    assert _draft_identity(result, card['gb']) == _draft_identity(source, source['webFiles']['bindings']['c_gb'])


def test_cli_i24c_schema1_defaults_and_independent_arrays(tmp_path):
    source = read_session_payload(Path(__file__).parent / 'fixtures/sessions/single.v41-bindings1.json')
    bindings = source['webFiles']['bindings']
    leaf = {'resourceId': bindings['c_gb']['resourceId']}
    bindings['c_gb'] = leaf
    bindings['c_depth'] = [None, [leaf, None], []]
    result = _replay_cli_sidecar(source, tmp_path)
    actual = result['webFiles']['bindings']
    descriptor = source['resources'][leaf['resourceId']]
    (card,) = actual['circularSources']
    assert card['gb']['name'] == descriptor['name']
    assert card['gb']['type'] == descriptor['type']
    assert card['gb']['lastModified'] == descriptor['lastModified']
    assert actual['c_depth'] == [None, [card['gb'], None], []]


@pytest.mark.parametrize('modified', [0, 0.125])
def test_cli_i24d_schema2_ordinary_metadata_is_exact(tmp_path, modified):
    source = _cli_composite_source()
    bindings = source['webFiles']['bindings']
    bindings['c_gb'] = {**bindings['c_gb']['components'][0], 'name': '', 'type': '', 'lastModified': modified}
    expected = _draft_identity(source, bindings['c_gb'])
    result = _replay_cli_sidecar(source, tmp_path)
    (card,) = result['webFiles']['bindings']['circularSources']
    assert _draft_identity(result, card['gb']) == expected


@pytest.mark.parametrize('field,slot', [
    ('conservationLosatFastaSources', 'c_conservation_fastas'),
    ('conservationSequenceSources', 'c_conservation_sequence_sources'),
])
@pytest.mark.parametrize('direct_populated', [False, True])
@pytest.mark.parametrize('explicit', ['absent', 'ordinary', 'null', 'empty', 'longer'])
def test_cli_i24e_explicit_slots_override_direct_sources(tmp_path, field, slot, direct_populated, explicit):
    source = _cli_composite_source()
    web = source['webFiles']
    bindings = web['bindings']
    part = bindings['c_gb']['components'][0]
    web[field] = [part['resourceId']] if direct_populated else []
    web['resourceOriginalNames'] = {part['resourceId']: 'historical-name.gb'}
    if explicit == 'absent':
        del bindings[slot]
    else:
        bindings[slot] = {
            'ordinary': part, 'null': None, 'empty': [],
            'longer': [part, {**part, 'name': 'second-occurrence.gb'}],
        }[explicit]
    expected = [_draft_identity(source, component) for component in bindings['c_gb']['components']]
    result = _replay_cli_sidecar(source, tmp_path)
    rebuilt = result['webFiles']
    actual = rebuilt['bindings'][slot]
    assert _card_files(result, rebuilt['bindings']['circularSources']) == expected
    if explicit != 'absent':
        assert _draft_identity(result, actual) == _draft_identity(source, bindings[slot])
    elif direct_populated:
        assert actual[0]['name'] == 'historical-name.gb'
    else:
        assert actual == []
    leaves = actual if isinstance(actual, list) else [actual] if actual else []
    assert rebuilt[field] == [leaf['resourceId'] for leaf in leaves]
    assert all(resource_id in result['resources'] for resource_id in rebuilt[field])
    if explicit in ('ordinary', 'longer'):
        assert rebuilt['resourceOriginalNames'][leaves[0]['resourceId']] == leaves[0]['name']


@pytest.mark.parametrize('invalid', ['missing', 'metadata', 'dangling', 'unsupported', 'nested'])
@pytest.mark.parametrize('sentinels', [False, True])
def test_cli_i24f_invalid_bindings_reject_before_render_or_publication(tmp_path, monkeypatch, invalid, sentinels):

    source = _cli_composite_source()
    bindings = source['webFiles']['bindings']
    composite = bindings['c_gb']
    source['webFiles']['conservationLosatFastaSources'] = [composite['components'][0]['resourceId']]
    if invalid == 'missing':
        del bindings['c_gb']
    elif invalid == 'metadata':
        composite['components'][0]['lastModified'] = -1
    elif invalid == 'dangling':
        composite['components'][0]['resourceId'] = 'missing-resource'
    elif invalid == 'unsupported':
        bindings['c_fasta'] = copy.deepcopy(composite)
    elif invalid == 'nested':
        composite['components'][0] = copy.deepcopy(composite)
    source_path = tmp_path / 'invalid.json'
    write_unchecked_session_file(source_path, source)
    original = source_path.read_bytes()
    svg, sidecar = tmp_path / 'rejected.svg', tmp_path / 'rejected.gbdraw'
    if sentinels:
        svg.write_bytes(b'existing diagram')
        sidecar.write_bytes(b'existing sidecar')
    calls = []
    import gbdraw.session as session_module

    monkeypatch.setattr(session_module, '_render_session_drawing_plans', lambda *a, **kw: calls.append(a))
    with pytest.raises((ValidationError, SystemExit)):
        circular_main([
            '--session', str(source_path), '-o', str(svg.with_suffix('')),
            '--session_output', str(sidecar), '--overwrite',
        ])
    assert calls == []
    assert source_path.read_bytes() == original
    if sentinels:
        assert (svg.read_bytes(), sidecar.read_bytes()) == (b'existing diagram', b'existing sidecar')
    else:
        assert not svg.exists() and not sidecar.exists()


def test_cli_i24_explicit_null_inventory_survives_without_direct_fields(tmp_path):
    source = _cli_composite_source()
    source['webFiles']['bindings'] = {'schema': 2, 'c_gb': None}
    result = _replay_cli_sidecar(source, tmp_path)
    assert result['webFiles']['bindings']['circularSources'] == []


def test_cli_projection_request_only_session_keeps_no_inventory():
    from gbdraw.session_migration import _project_web_file_inventory

    source = _cli_composite_source()
    del source['webFiles']
    load_session_document(source)
    document, sources = session_io_module.split_resource_bytes(source)
    assert _project_web_file_inventory(document, SessionResources(document['resources'], sources)) is None


@pytest.mark.parametrize('invalid', ['mixed-leaf', 'mixed-composite', 'nested', 'null'])
def test_internal_composite_transport_rejects_unsupported_shapes(invalid):
    source = _cli_composite_source()
    from gbdraw.session_migration import _project_web_file_inventory

    load_session_document(source)
    document, sources = session_io_module.split_resource_bytes(source)
    inventory = _project_web_file_inventory(document, SessionResources(document['resources'], sources))
    composite = inventory['c_gb']
    if invalid == 'mixed-leaf':
        composite['components'][0]['data'] = 'eA=='
    elif invalid == 'mixed-composite':
        composite['resourceId'] = 'not-a-leaf'
    elif invalid == 'nested':
        composite['components'][0] = copy.deepcopy(composite)
    else:
        composite['components'][0] = None
    with pytest.raises(ValidationError):
        session_io_module._attach_current_web_file_bindings({'resources': {}}, inventory, {})


def _linear_session_with_cards(source_session: dict | None = None, **card_fields: str) -> dict:
    cards = [{"uid": uid, "gb": None, "gff": None, "fasta": None, **card_fields} for uid in ("seq-1", "seq-2")]
    return build_mode_scoped_session_json(
        SessionBuildContext(
            mode="linear", output_prefix="out", render_formats=("svg",), source_session=source_session,
        ),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": cards},
        generated_at=datetime(2026, 10, 10),
        canonical_request=_canonical_request("linear"),
    )


def test_session_mode_slices_keep_records_off_and_reject_unknown_keys() -> None:
    # Record selection D-05: each slice lists its OFF records once each, a
    # Linear card uid bound in webFiles or a Circular source selector `#N`
    # (a Session 50); the move to Session 51 makes the selectors record
    # identities of the File card (`<uid>:<k>`).
    payload = _linear_session_with_cards()
    payload["resources"]["record-1-genbank"] = {**_file_entry("circular.gb", b"LOCUS\n"), "encoding": "base64"}
    payload["webFiles"]["bindings"]["c_gb"] = {
        "resourceId": "record-1-genbank", "name": "circular.gb", "type": "", "lastModified": 0,
    }
    modes = payload.setdefault("modes", {})
    modes.setdefault("linear", {}).setdefault("config", {})["recordsOff"] = ["seq-2"]
    modes.setdefault("circular", {}).setdefault("config", {})["recordsOff"] = ["#2", "#10"]
    validate_session(payload)

    for mode, records_off in (
        ("linear", ["not-a-card"]),
        ("linear", ["seq-2", "seq-2"]),
        ("linear", "seq-2"),
        ("linear", [""]),
        ("circular", ["contig_2"]),
        ("circular", ["#0"]),
        ("circular", [2]),
    ):
        invalid = copy.deepcopy(payload)
        invalid["modes"][mode]["config"]["recordsOff"] = records_off
        with pytest.raises(ValidationError) as caught:
            validate_session(invalid)
        assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}, (
            mode, records_off,
        )

    # A CLI re-save keeps the OFF records of both slices.
    resaved = build_session_payload(
        SessionBuildContext(mode="linear", output_prefix="again", render_formats=("svg",), source_session=payload),
        svg_results=(("again", "<svg></svg>"),),
        embedded_files={
            "linearSeqs": payload["webFiles"]["bindings"]["linearSeqs"],
            "c_gb": payload["resources"]["record-1-genbank"],
        },
        generated_at=datetime(2026, 10, 11),
        canonical_request=_canonical_request("linear"),
    )
    # The OFF Circular records make the move list an (undrawn) Circular drawing.
    assert drawing_view(resaved, mode="linear")["modes"]["linear"]["config"]["recordsOff"] == ["seq-2"]
    assert drawing_view(resaved, mode="circular")["modes"]["circular"]["config"]["recordsOff"] == ["record-1:2", "record-1:10"]
    # A Session 51 Circular slice names records by identity of a bound File
    # card, or by a source selector the move could not place (OV-2869).
    circular = next(drawing for drawing in resaved["drawings"] if drawing["mode"] == "circular")
    kept = copy.deepcopy(resaved)
    next(d for d in kept["drawings"] if d["mode"] == "circular")["config"]["recordsOff"] = ["#2"]
    validate_session(kept)
    for records_off in (["#0"], ["record-2"], ["record-1:0"]):
        invalid = copy.deepcopy(resaved)
        next(d for d in invalid["drawings"] if d["mode"] == "circular")["config"]["recordsOff"] = records_off
        with pytest.raises(ValidationError):
            validate_session(invalid)
    assert circular["webFiles"]["bindings"]["circularSources"][0]["uid"] == "record-1"


def test_cli_resave_keeps_file_defaults_of_linear_cards() -> None:
    # OV-402: the File definition, File subtitle, and inferred definition of a
    # Linear card survive a CLI re-save, as the Web writer saves them.
    defaults = {"file_definition": "File def", "file_subtitle": "File sub", "inferred_definition": "Inferred"}
    payload = _linear_session_with_cards(**defaults)
    validate_session(payload)
    for card in payload["webFiles"]["bindings"]["linearSeqs"]:
        assert {key: card.get(key) for key in defaults} == defaults


def test_legacy_renames_are_spelled_only_by_their_owner() -> None:
    """OV-818: the `sqrt` and `repeat_region` renames are read from one owner."""

    import ast

    import gbdraw.session_io as session_io_module
    import gbdraw.session_request_codec as codec_module
    from gbdraw.session_legacy_values import LEGACY_FEATURE_SHAPES, LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES

    renamed = {*LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES, *LEGACY_FEATURE_SHAPES}
    assert renamed == {"sqrt", "repeat_region"}
    for module in (session_io_module, codec_module):
        tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
        literals = sorted(
            f"{node.value}:{node.lineno}"
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value in renamed
        )
        assert literals == [], module.__name__


def test_a_cli_resave_names_a_losat_entry_without_a_text_name_as_the_web_does(
    tmp_path: Path,
) -> None:
    # OV-1007: the CLI writer copied a saved `{}` or `7` through. The Web reads
    # such a name as the entry's default, `losat_pair_<n>.tsv` (`n` counts the
    # saved entries from 1), and saves that name: tests/web/session-losat-cache-
    # validation.test.mjs ("... file name that is no text loads as the default
    # name") pins the same values for the same edit of this fixture.
    source = Path(__file__).parent / "fixtures" / "sessions" / "synthetic_conservation.gbdraw-session.json.gz"
    session = json.loads(gzip.decompress(source.read_bytes()))
    first, second, third = session["losatCache"]["entries"]
    first["filename"] = {}
    second["filename"] = 7
    edited = tmp_path / "edited.gbdraw"
    edited.write_text(json.dumps(session), encoding="utf-8")

    sidecar = tmp_path / "replay.gbdraw"
    circular_main(["--session", str(edited), "--output", str(tmp_path / "replay"),
                   "--format", "svg", "--session_output", str(sidecar)])
    assert [entry["filename"] for entry in drawing_view(load_session(sidecar), mode="circular")["losatCache"]["entries"]] == [
        "losat_pair_1.tsv", "losat_pair_2.tsv", third["filename"]
    ]

    first["filename"] = ""
    second["filename"] = ["b.tsv"]
    second["display"] = False
    written = copy.deepcopy(session)
    session_io_module.normalize_current_session_artifacts(written)
    assert [(entry["filename"], entry["display"]) for entry in written["losatCache"]["entries"]] == [
        ("losat_pair_1.tsv", True), ("", False), (third["filename"], True)
    ]


def test_a_cli_resave_keeps_a_circular_drawing_without_web_file_bindings(tmp_path: Path) -> None:
    # OV-3011: a Gallery publication writes a Circular drawing without
    # `webFiles.bindings`; its comparison FASTAs are the resource IDs in
    # `conservationLosatFastaSources`. The CLI writer made bindings from those
    # IDs alone, with no Circular File cards, and refused its own Session
    # ("requires circularSources"). The re-save keeps the drawing as it was.
    from gbdraw.session import upgrade_session_document
    from gbdraw.session_resources import drawing_resource_ids

    source = Path(__file__).parent / "fixtures" / "sessions" / "synthetic_conservation.gbdraw-session.json.gz"
    upgraded = upgrade_session_document(source, temporary_directory=tmp_path).document
    session = embedded_session(upgraded.to_dict(), upgraded.resources.sources)
    drawing = next(drawing for drawing in session["drawings"] if drawing["mode"] == "circular")
    web_files = drawing["webFiles"]
    del web_files["bindings"]
    fasta_ids = web_files["conservationLosatFastaSources"]
    assert fasta_ids
    # The publication keeps only the resources that its drawings name.
    named = set().union(*(drawing_resource_ids(each) for each in session["drawings"]))
    session["resources"] = {key: value for key, value in session["resources"].items() if key in named}
    edited = tmp_path / "published.gbdraw"
    write_session_file(edited, session)

    sidecar = tmp_path / "replay.gbdraw"
    circular_main(["--session", str(edited), "--output", str(tmp_path / "replay"),
                   "--format", "svg", "--session_output", str(sidecar)])

    saved = read_session_payload(sidecar)
    saved_drawing = next(saved_drawing for saved_drawing in saved["drawings"] if saved_drawing["id"] == drawing["id"])
    assert "bindings" not in saved_drawing["webFiles"]
    assert saved_drawing["webFiles"]["conservationLosatFastaSources"] == fasta_ids
    assert [saved["resources"][resource_id]["data"] for resource_id in fasta_ids] == [
        session["resources"][resource_id]["data"] for resource_id in fasta_ids
    ]


def test_a_cli_resave_writes_a_hidden_losat_entry_without_a_name_as_the_web_does() -> None:
    # OV-1013: a hidden (`display: false`) raw entry is cache only. The Web Load
    # gives it no Raw LOSAT row (config.js applyLosatCache) and its Save writes
    # an entry without a row with `filename: ""` (serializeLosatCache;
    # tests/web/session-export-validation.test.mjs pins it). The CLI writer kept
    # a hidden entry's text name.
    source = Path(__file__).parent / "fixtures" / "sessions" / "synthetic_conservation.gbdraw-session.json.gz"
    session = json.loads(gzip.decompress(source.read_bytes()))
    first, second, third = session["losatCache"]["entries"]
    assert isinstance(second["filename"], str) and second["filename"]
    second["display"] = False
    session_io_module.normalize_current_session_artifacts(session)
    assert [(entry["filename"], entry["display"]) for entry in session["losatCache"]["entries"]] == [
        (first["filename"], True), ("", False), (third["filename"], True)
    ]


_PROTEIN_GALLERY_SESSION = (
    Path(__file__).resolve().parents[1]
    / "gbdraw" / "web" / "gallery" / "sessions" / "BGC0000708-BGC0000713.gbdraw"
)


def _container_types(value: object) -> set[type]:
    pending, found = [value], set()
    while pending:
        current = pending.pop()
        if isinstance(current, (dict, list, tuple)):
            found.add(type(current))
            pending.extend(current.values() if isinstance(current, dict) else current)
    return found


def test_session_document_to_dict_returns_a_plain_mutable_manifest() -> None:
    document = load_session_document(_PROTEIN_GALLERY_SESSION)

    data = document.to_dict()

    manifest = drawing_view(data)["proteinIdentityManifest"]
    assert _container_types(manifest) == {dict, list}
    binding = next(iter(manifest["recordInstances"].values()))
    binding["runtimeIds"].clear()
    manifest["extra"] = 1
    assert drawing_view(document.to_dict())["proteinIdentityManifest"] != manifest


def test_a_python_manifest_with_tuples_takes_the_copy_path_through_render_and_writer() -> None:
    from gbdraw.api.request_render import CurrentRequestArtifacts

    # The drawing's view is the 48 layout this case exercises (a 49 holds no top-level manifest).
    payload = view_at(drawing_view(read_session_payload(_PROTEIN_GALLERY_SESSION)), 48)
    raw = copy.deepcopy(payload["proteinIdentityManifest"])
    protein = next(iter(raw["proteinSets"].values()))["proteins"][0]
    protein["locationParts"] = tuple(tuple(part) for part in protein["locationParts"])
    json_clone = json.loads(json.dumps(raw))

    artifacts = CurrentRequestArtifacts(
        losat_cache_entries=tuple(payload["losatCache"]["entries"]),
        protein_identity_manifest=raw,
    )

    # As before: CurrentRequestArtifacts keeps its deep copy, with the tuples.
    stored = artifacts.protein_identity_manifest
    assert type(stored) is dict and stored == raw
    assert type(next(iter(stored["proteinSets"].values()))["proteins"][0]["locationParts"]) is tuple
    normalize_current_session_artifacts(payload, protein_identity_manifest=stored)
    # The writer writes a JSON clone of it.
    assert payload["proteinIdentityManifest"] == json_clone
    assert _container_types(payload["proteinIdentityManifest"]) == {dict, list}
    assert json.dumps(payload["proteinIdentityManifest"]) == json.dumps(json_clone)
    document = SessionDocument(payload)
    assert document.to_dict()["proteinIdentityManifest"] == json_clone

    # A Session payload that still holds the tuples is checked as raw data.
    payload["proteinIdentityManifest"] = raw
    assert type(SessionDocument(payload)._data["proteinIdentityManifest"]) is dict


@pytest.mark.parametrize(
    ("raw_args", "plan"),
    [
        (
            ["--gbk", "P1.gb", "P2.gb", "--protein_blastp_mode", "orthogroup"],
            {"mode": "adjacent", "defaultSource": "losat", "edges": []},
        ),
        (["--gbk", "P1.gb", "P2.gb", "--protein-blastp-mode=collinear"], {"mode": "adjacent", "defaultSource": "losat", "edges": []}),
        (["--gbk", "P1.gb", "P2.gb", "--protein_blastp_mode", "orthogroup", "--protein_blastp_mode", "none"], None),
        (["--gbk", "P1.gb", "P2.gb"], None),
    ],
)
def test_cli_only_legacy_draft_states_its_protein_comparison_as_the_adjacent_plan(
    raw_args: list[str], plan: dict[str, object] | None
) -> None:
    """A legacy CLI writer selected LOSATP exactly when the last --protein_blastp_mode
    was not ``none``; that run compared adjacent records (OV-1929 f1, the Web Load's rule)."""

    migrated_config, migrated_files = (
        session_io_module.migrate_legacy_linear_comparison_draft_for_current_writer(
            {"cliOptions": {"schema": 1, "mode": "argv", "rawArgs": raw_args}, "linearComparisonPlan": None},
            {"linearSeqs": [{"uid": "a"}, {"uid": "b"}]},
            force_web_draft=False,
        )
    )

    assert migrated_config.get("linearComparisonPlan") == plan
    assert migrated_files["linearComparisons"] == []


def test_python_card_bindings_keep_the_card_windows_as_the_web_writes_them() -> None:
    # One card shape (the Web's `linearBindings` in services/session-resources.js): a card's windows travel with it.
    windows = [{"start": 16000, "end": 100, "reverseComplement": False, "column": None}]
    payload: dict = {"resources": {}}
    session_io_module._attach_current_web_file_bindings(payload, {"linearSeqs": [{"uid": "a", "windows": windows}]}, {})
    card = payload["webFiles"]["bindings"]["linearSeqs"][0]
    assert card["windows"] == windows


def test_legacy_feature_index_keys_an_all_input_record_with_the_expanded_record_key(monkeypatch) -> None:
    """The legacy feature index names each record of a several-record ALL input
    with the key the resolver gives it (``expanded_record_key``), so a change of
    the separator moves both together."""

    import gbdraw.api.requests as requests_module

    monkeypatch.setattr(requests_module, "EXPANDED_RECORD_KEY_SEPARATOR", "|")
    legacy = {
        "records": [{"recordKey": "card", "cardinality": "all"}],
        "features": [
            {"fileIdx": 0, "record_idx": index, "stable_feature_id": f"hash-{index}"}
            for index in (0, 1)
        ],
    }

    index = session_io_module._legacy_feature_index(legacy, "linear")

    assert [json.loads(key)[1] for key, *_ in index] == [
        requests_module.expanded_record_key("card", 1),
        requests_module.expanded_record_key("card", 2),
    ] == ["card|1", "card|2"]


def _session_with_plan_edge(edge_extra: dict[str, object]) -> dict:
    fixture = (
        Path(__file__).resolve().parent
        / "fixtures/sessions/drawings-linear-pair-losat.v51.gbdraw"
    )
    payload = json_document(read_session_payload(fixture))
    payload["drawings"][0]["config"]["linearComparisonPlan"]["edges"] = [
        {
            "id": "edge-a-b",
            "queryUid": "a",
            "subjectUid": "b",
            "included": True,
            "fileActive": False,
            "losatFilenameActive": False,
            "source": "losat",
            "losatFilename": "",
            **edge_extra,
        }
    ]
    return payload


def test_current_comparison_plan_admits_window_pair_names() -> None:
    """OV-2742: a name typed on a window pair is saved in the record pair's plan edge."""
    validate_session(
        _session_with_plan_edge({"windowLosatFilenames": {"a~w1->b": "mine.tsv"}})
    )
    # A Web Save Session with a typed window-pair name (see its .provenance.json).
    saved = (
        Path(__file__).resolve().parent
        / "fixtures/sessions/web-linear-window-pair-names.v51.gbdraw"
    )
    validate_session(json_document(read_session_payload(saved)))
    for bad in ({"a~w1->b": ""}, {" ": "mine.tsv"}, {"a~w1->b": 3}, ["mine.tsv"]):
        with pytest.raises(ValidationError, match="windowLosatFilenames") as excinfo:
            validate_session(_session_with_plan_edge({"windowLosatFilenames": bad}))
        assert excinfo.value.diagnostic == {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


def test_session_50_plan_edge_with_window_pair_names_is_refused() -> None:
    """OV-2742: no Session 46-50 writer wrote window-pair names, so such a file is refused."""
    payload = build_mode_scoped_session_json(
        SessionBuildContext(mode="linear", output_prefix="out", render_formats=("svg",)),
        svg_results=(("out", "<svg></svg>"),),
        embedded_files={"linearSeqs": []},
        generated_at=datetime(2026, 11, 10),
        canonical_request=_canonical_request("linear"),
    )
    assert payload["version"] == 50
    edge = {
        "id": "edge-a-b", "queryUid": "a", "subjectUid": "b", "included": True, "fileActive": False,
        "losatFilenameActive": False, "source": "losat", "losatFilename": "",
    }
    payload["modes"] = {"linear": {"config": {"form": {}, "adv": {}, "linearComparisonPlan": {
        "mode": "adjacent", "defaultSource": "losat", "edges": [edge],
    }}}}
    validate_session(payload)
    edge["windowLosatFilenames"] = {"a~w1->b": "mine.tsv"}
    with pytest.raises(ValidationError, match="unknown field\\(s\\): windowLosatFilenames"):
        validate_session(payload)


def _catalog_rows(catalog: dict) -> list[dict]:
    return [row for item in catalog["items"] for row in item["biologicalFeatures"]]


def test_cli_saved_session_51_writes_catalog_schema_6_without_sequences(tmp_path: Path) -> None:
    """D5: the stored catalog holds no sequence; the /translation qualifier stays as it was read."""

    from Bio import SeqIO
    from gbdraw.web_support.feature_catalog import STORED_SEQUENCE_KEYS

    record = SeqRecord(Seq("ATGAAATAAGGG" * 30), id="saved", name="saved", annotations={"molecule_type": "DNA"})
    record.features = [
        SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS", qualifiers={"translation": ["MK"]}),
    ]
    SeqIO.write(record, tmp_path / "saved.gb", "genbank")
    session = tmp_path / "saved.gbdraw"
    circular_main(["--gbk", str(tmp_path / "saved.gb"), "-o", str(tmp_path / "saved"), "--session_output", str(session)])

    payload = read_session_payload(session)
    assert payload["version"] == CURRENT_SESSION_VERSION
    [drawing] = payload["drawings"]
    catalog = drawing["editorState"]["featureCatalog"]
    assert catalog["schema"] == CURRENT_FEATURE_CATALOG_SCHEMA == 6
    [row] = _catalog_rows(catalog)
    assert not STORED_SEQUENCE_KEYS & set(row)
    assert row["qualifiers"]["translation"] == ["MK"]
    assert all("sequenceSources" not in item for item in catalog["items"])


def test_session_47_catalog_schema_5_loads_without_its_stored_sequences() -> None:
    """A Session 46-50 still loads; its stored sequences are not used: its move to drawings promotes the catalog to
    schema 6, giving back each /translation qualifier a translation mark stood for."""

    from gbdraw.web_support.feature_catalog import STORED_SEQUENCE_KEYS

    path = Path(__file__).parent / "fixtures" / "sessions" / "BGC0000708-BGC0000713.v47.gbdraw-session.json.gz"
    stored = read_session_payload(path)
    stored_catalog = stored["editorState"]["featureCatalog"]
    assert stored_catalog["schema"] == SEQUENCE_FEATURE_CATALOG_SCHEMA
    marked = {
        (row["recordKey"], row["biologicalFeatureId"]): row["amino_acid_sequence"]
        for row in _catalog_rows(stored_catalog)
        if row.get("translationFromAminoAcidSequence")
    }
    assert marked and any("sequenceSourceIndex" in row for row in _catalog_rows(stored_catalog))

    from gbdraw.session_io import drawings_document_of_mode_scoped

    document = load_session_document(path).to_dict()
    # Python keeps the Session 47 as it was stored; its move to drawings (the Web's Load, a resave) promotes it.
    assert document["editorState"]["featureCatalog"]["schema"] == SEQUENCE_FEATURE_CATALOG_SCHEMA
    moved, _facts = drawings_document_of_mode_scoped(document)
    loaded = drawing_view(dict(moved))["editorState"]["featureCatalog"]

    assert loaded["schema"] == CURRENT_FEATURE_CATALOG_SCHEMA
    assert all("sequenceSources" not in item for item in loaded["items"])
    rows = _catalog_rows(loaded)
    assert len(rows) == len(_catalog_rows(stored_catalog))
    assert not any(STORED_SEQUENCE_KEYS & set(row) for row in rows)
    assert {
        (row["recordKey"], row["biologicalFeatureId"]): row["qualifiers"]["translation"][0]
        for row in rows
        if (row["recordKey"], row["biologicalFeatureId"]) in marked
    } == marked
