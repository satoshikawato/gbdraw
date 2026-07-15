from __future__ import annotations

import base64
import re
from datetime import datetime
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.circular as circular_cli_module
from gbdraw.circular import circular_main
from gbdraw.linear import _get_args as get_linear_args
from gbdraw.cli_utils.session import (
    collect_embedded_files_from_cli_args,
    make_rendered_svg,
    parse_session_pre_args,
    resolve_session_sidecar_path,
    strip_session_output_args,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import read_records_table
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
    SESSION_FORMAT,
    SUPPORTED_SESSION_VERSIONS,
    SessionBuildContext,
    build_session_json,
    decode_depth_payload,
    load_session,
    materialize_embedded_file,
    session_to_cli_args,
    validate_session,
)


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


def test_current_session_version_matches_web_config() -> None:
    source = Path("gbdraw/web/js/services/config.js").read_text(encoding="utf-8")
    match = re.search(r"const\s+SESSION_VERSION\s*=\s*(\d+);", source)
    assert match is not None
    assert int(match.group(1)) == CURRENT_SESSION_VERSION


def test_future_session_version_fails() -> None:
    session = _minimal_session({})
    session["version"] = CURRENT_SESSION_VERSION + 1

    with pytest.raises(ValidationError, match="newer"):
        validate_session(session)


def test_session_version_27_remains_supported() -> None:
    session = _minimal_session({})
    session["version"] = 27

    validate_session(session)
    assert 27 in SUPPORTED_SESSION_VERSIONS


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
    assert payload["files"]["c_gb"]["data"]
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
