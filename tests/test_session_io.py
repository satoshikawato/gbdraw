from __future__ import annotations

import base64
import re
from pathlib import Path

import pytest

from gbdraw.circular import circular_main
from gbdraw.cli_utils.session import (
    make_rendered_svg,
    parse_session_pre_args,
    resolve_session_sidecar_path,
)
from gbdraw.exceptions import ValidationError
from gbdraw.render.formats import ACCEPTED_FORMATS
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    DEPTH_FILE_ENCODING,
    SESSION_FORMAT,
    SUPPORTED_SESSION_VERSIONS,
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
        "version": CURRENT_SESSION_VERSION,
        "createdAt": "2026-06-22T00:00:00Z",
        "config": {"form": {"prefix": "out"}, "adv": {}},
        "ui": {"mode": mode, "cInputType": "gb", "lInputType": "gb"},
        "files": files,
    }


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
    assert spec.args[-2:] == ("-f", "interactive-svg")
    assert spec.cli_invocation_args[0:2] == ("-o", "new")
    assert spec.cli_invocation_args[-2:] == ("-f", "interactive-svg")
    assert spec.file_bindings[0].argIndex == 3
    assert spec.file_bindings[0].slot == "files.c_gb"


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
    session["config"]["losat"] = {
        "executionMode": "auto",
        "threadsPerJob": "2",
        "blastp": {
            "mode": "collinear",
            "maxHits": 7,
            "collinearMinAnchors": 2,
            "collinearMaxGeneGap": 3,
            "collinearMaxDiagonalDrift": 4,
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
    assert spec.args[spec.args.index("--collinear_search_scope") + 1] == "all"
    assert spec.args[spec.args.index("--collinear_color_mode") + 1] == "orientation_identity"


def test_session_pre_parse_rejects_unsupported_options() -> None:
    with pytest.raises(SystemExit):
        parse_session_pre_args(
            ["--session", "session.gbdraw-session.json", "--gbk", "input.gb"],
            mode="circular",
        )


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


def test_circular_cli_save_session_round_trip(tmp_path: Path, test_inputs_dir: Path) -> None:
    output_prefix = tmp_path / "circular_session"
    circular_main(
        [
            "--gbk",
            str(test_inputs_dir / "MjeNMV.gb"),
            "-o",
            str(output_prefix),
            "-f",
            "svg",
            "--save-session",
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
