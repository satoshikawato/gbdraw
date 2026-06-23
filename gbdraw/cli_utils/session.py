#!/usr/bin/env python
# coding: utf-8

"""Shared CLI helpers for GUI session JSON input and sidecar output."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from gbdraw.exceptions import ValidationError
from gbdraw.render.formats import SVG_FORMAT, resolve_format_output_path
from gbdraw.session_io import (
    SessionBuildContext,
    SessionFileBinding,
    build_session_json,
    safe_embedded_filename,
    serialize_file_entry,
    write_session_json,
)


@dataclass(frozen=True)
class RenderedSvg:
    output_prefix: str
    svg_path: Path
    result_name: str


@dataclass(frozen=True)
class DiagramRunResult:
    mode: Literal["circular", "linear"]
    render_formats: tuple[str, ...]
    outputs: tuple[RenderedSvg, ...]
    losat_cache_entries: tuple[Mapping[str, Any], ...] | None = None


@dataclass(frozen=True)
class SessionCliRequest:
    session_path: str
    output: str | None
    format: str | None
    save_session: bool
    session_output: str | None


def add_session_args(parser: argparse.ArgumentParser) -> None:
    """Add session input/output options to a diagram parser."""

    parser.add_argument(
        "--session",
        help="Regenerate a diagram from a gbdraw GUI session JSON file.",
        type=str,
    )
    parser.add_argument(
        "--save-session",
        help="Write one GUI-loadable .gbdraw-session.json sidecar for this run.",
        action="store_true",
    )
    parser.add_argument(
        "--session-output",
        metavar="PATH",
        help="Write the session sidecar to PATH; implies --save-session.",
        type=str,
    )


def parse_session_pre_args(
    cmd_args: Sequence[str],
    *,
    mode: Literal["circular", "linear"],
) -> SessionCliRequest | None:
    """Pre-parse --session invocations and reject unsupported override options."""

    if "-h" in cmd_args or "--help" in cmd_args:
        return None
    if "--session" not in cmd_args:
        return None

    parser = argparse.ArgumentParser(
        prog=f"gbdraw {mode}",
        add_help=False,
    )
    parser.add_argument("--session", required=True)
    parser.add_argument("-o", "--output")
    parser.add_argument("-f", "--format")
    parser.add_argument("--save-session", action="store_true")
    parser.add_argument("--session-output")
    namespace, unknown = parser.parse_known_args(list(cmd_args))
    if unknown:
        parser.error(
            "--session cannot be combined with unsupported option(s): "
            + " ".join(unknown)
        )
    return SessionCliRequest(
        session_path=str(namespace.session),
        output=namespace.output,
        format=namespace.format,
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
    )


def validate_session_override_args(
    cmd_args: Sequence[str],
    *,
    mode: Literal["circular", "linear"],
) -> SessionCliRequest | None:
    """Compatibility wrapper for the documented session pre-parse helper name."""

    return parse_session_pre_args(cmd_args, mode=mode)


def resolve_session_sidecar_path(
    *,
    explicit_path: str | None,
    output_prefix: str | None,
    outputs: Sequence[RenderedSvg],
) -> Path:
    """Resolve the run-level session sidecar path."""

    if explicit_path:
        return Path(explicit_path)
    if output_prefix:
        return Path(f"{output_prefix}.gbdraw-session.json")
    if len(outputs) == 1:
        return outputs[0].svg_path.with_suffix(".gbdraw-session.json")
    return Path("gbdraw.gbdraw-session.json")


def make_rendered_svg(output_prefix: str, result_name: str | None = None) -> RenderedSvg:
    """Create a RenderedSvg result for the static SVG written by save_figure()."""

    svg_path = Path(resolve_format_output_path(output_prefix, SVG_FORMAT))
    return RenderedSvg(
        output_prefix=str(output_prefix),
        svg_path=svg_path,
        result_name=result_name or svg_path.stem,
    )


def save_session_sidecar_if_requested(
    *,
    save_session: bool,
    session_output: str | None,
    output_prefix: str | None,
    run_result: DiagramRunResult,
    cmd_args: Sequence[str] | None = None,
    source_session: Mapping[str, Any] | None = None,
    cli_invocation_args: Sequence[str] = (),
    file_bindings: Sequence[SessionFileBinding] = (),
) -> Path | None:
    """Build and write a GUI session sidecar when requested."""

    if not save_session and not session_output:
        return None
    sidecar_path = resolve_session_sidecar_path(
        explicit_path=session_output,
        output_prefix=output_prefix,
        outputs=run_result.outputs,
    )

    if source_session is not None:
        embedded_files = source_session.get("files")
        if not isinstance(embedded_files, Mapping):
            raise ValidationError("Source session has no files to preserve.")
        session_files: Mapping[str, Any] = embedded_files
        invocation_args = tuple(str(arg) for arg in cli_invocation_args)
        bindings = tuple(file_bindings)
    else:
        invocation_args = tuple(strip_session_output_args(cmd_args or ()))
        session_files, bindings = collect_embedded_files_from_cli_args(
            run_result.mode,
            invocation_args,
        )

    svg_results = _read_svg_results(run_result.outputs)
    context_output_prefix = output_prefix
    if context_output_prefix is None and len(run_result.outputs) == 1:
        context_output_prefix = run_result.outputs[0].output_prefix
    payload = build_session_json(
        SessionBuildContext(
            mode=run_result.mode,
            output_prefix=context_output_prefix,
            render_formats=run_result.render_formats,
            source_session=source_session,
            cli_invocation_args=invocation_args,
            file_bindings=tuple(bindings),
        ),
        svg_results=svg_results,
        embedded_files=session_files,
        generated_at=datetime.now(timezone.utc),
        losat_cache_entries=run_result.losat_cache_entries,
    )
    write_session_json(sidecar_path, payload)
    return sidecar_path


def strip_session_output_args(cmd_args: Sequence[str]) -> list[str]:
    """Remove sidecar-output-only flags before storing cliInvocation.args."""

    result: list[str] = []
    index = 0
    while index < len(cmd_args):
        token = str(cmd_args[index])
        if token == "--save-session":
            index += 1
            continue
        if token == "--session-output":
            index += 2
            continue
        result.append(token)
        index += 1
    return result


def collect_embedded_files_from_cli_args(
    mode: Literal["circular", "linear"],
    cli_args: Sequence[str],
) -> tuple[dict[str, Any], tuple[SessionFileBinding, ...]]:
    """Embed local CLI input files and build cliInvocation file bindings."""

    files = _empty_files_payload()
    bindings: list[SessionFileBinding] = []
    circular_counts: dict[str, int] = {}
    linear_depth_track_index = 0
    circular_depth_index = 0

    index = 0
    while index < len(cli_args):
        token = str(cli_args[index])
        if mode == "circular" and token in {"--gbk", "--gff", "--fasta", "--conservation_blast", "--depth_track"}:
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                if token == "--gbk":
                    ordinal = circular_counts.get("gbk", 0)
                    slot = "files.c_gb" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["gbk"] = ordinal + 1
                elif token == "--gff":
                    ordinal = circular_counts.get("gff", 0)
                    slot = "files.c_gff" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["gff"] = ordinal + 1
                elif token == "--fasta":
                    ordinal = circular_counts.get("fasta", 0)
                    slot = "files.c_fasta" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["fasta"] = ordinal + 1
                elif token == "--conservation_blast":
                    slot = f"files.c_conservation_blasts[{len(files['c_conservation_blasts'])}]"
                else:
                    slot = f"files.c_depth[{circular_depth_index}]"
                    circular_depth_index += 1
                _set_file_slot(files, slot, value, depth=token == "--depth_track")
                bindings.append(_binding(arg_index, slot, value))
            index = next_index
            continue
        if mode == "linear" and token in {"--gbk", "--gff", "--fasta", "-b", "--blast", "--depth"}:
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                seq_index = offset
                if token == "--gbk":
                    slot = f"files.linearSeqs[{seq_index}].gb"
                    depth = False
                elif token == "--gff":
                    slot = f"files.linearSeqs[{seq_index}].gff"
                    depth = False
                elif token == "--fasta":
                    slot = f"files.linearSeqs[{seq_index}].fasta"
                    depth = False
                elif token in {"-b", "--blast"}:
                    slot = f"files.linearSeqs[{seq_index}].blast"
                    depth = False
                else:
                    slot = f"files.linearSeqs[{seq_index}].depth"
                    depth = True
                _set_file_slot(files, slot, value, depth=depth)
                bindings.append(_binding(arg_index, slot, value))
            index = next_index
            continue
        if mode == "linear" and token == "--depth_track":
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                slot = f"files.linearSeqs[{offset}].depth[{linear_depth_track_index}]"
                _set_file_slot(files, slot, value, depth=True)
                bindings.append(_binding(arg_index, slot, value))
            linear_depth_track_index += 1
            index = next_index
            continue
        if token == "--depth":
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                slot = "files.c_depth"
                _set_file_slot(files, slot, cli_args[value_index], depth=True)
                bindings.append(_binding(value_index, slot, cli_args[value_index]))
            index += 2
            continue
        if token in _COMMON_SINGLE_FILE_OPTIONS:
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                slot = _COMMON_SINGLE_FILE_OPTIONS[token]
                if slot == "files.cliInputs[]":
                    slot = _append_cli_input(files, cli_args[value_index], depth=False)
                else:
                    _set_file_slot(files, slot, cli_args[value_index], depth=False)
                bindings.append(_binding(value_index, slot, cli_args[value_index]))
            index += 2
            continue
        index += 1

    return files, tuple(bindings)


_COMMON_SINGLE_FILE_OPTIONS = {
    "-d": "files.d_color",
    "--default_colors": "files.d_color",
    "-t": "files.t_color",
    "--table": "files.t_color",
    "--label_whitelist": "files.whitelist",
    "--label_blacklist": "files.blacklist",
    "--qualifier_priority": "files.qualifier_priority",
    "--label_table": "files.cliInputs[]",
    "--feature_table": "files.cliInputs[]",
    "--collinear_blocks": "files.cliInputs[]",
}


def _empty_files_payload() -> dict[str, Any]:
    return {
        "c_gb": None,
        "c_gff": None,
        "c_fasta": None,
        "c_depth": None,
        "c_conservation_blasts": [],
        "c_conservation_fastas": [],
        "d_color": None,
        "t_color": None,
        "blacklist": None,
        "whitelist": None,
        "qualifier_priority": None,
        "linearSeqs": [],
        "cliInputs": [],
    }


def _read_svg_results(outputs: Sequence[RenderedSvg]) -> list[tuple[str, str]]:
    results: list[tuple[str, str]] = []
    for output in outputs:
        try:
            content = output.svg_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise ValidationError(
                f"Session sidecar output cannot read generated static SVG: {output.svg_path}"
            ) from exc
        results.append((output.result_name, content))
    return results


def _collect_option_values(args: Sequence[str], start_index: int) -> tuple[list[str], int]:
    values: list[str] = []
    index = start_index
    while index < len(args):
        token = str(args[index])
        if token.startswith("-") and token.lower() not in {"-", "none", "null"}:
            break
        values.append(token)
        index += 1
    return values, index


def _is_embeddable_path(value: object) -> bool:
    text = str(value or "").strip()
    if not text or text.lower() in {"-", "none", "null"}:
        return False
    return Path(text).is_file()


def _binding(arg_index: int, slot: str, path: object) -> SessionFileBinding:
    return SessionFileBinding(
        argIndex=arg_index,
        slot=slot,
        name=safe_embedded_filename(Path(str(path)).name, fallback="file"),
    )


def _append_cli_input(files: dict[str, Any], path: object, *, depth: bool) -> str:
    slot = f"files.cliInputs[{len(files['cliInputs'])}]"
    _set_file_slot(files, slot, path, depth=depth)
    return slot


def _set_file_slot(files: dict[str, Any], slot: str, path: object, *, depth: bool) -> None:
    entry = serialize_file_entry(str(path), depth=depth)
    if slot == "files.c_gb":
        files["c_gb"] = entry
        return
    if slot == "files.c_gff":
        files["c_gff"] = entry
        return
    if slot == "files.c_fasta":
        files["c_fasta"] = entry
        return
    if slot == "files.c_depth":
        files["c_depth"] = entry
        return
    if slot == "files.d_color":
        files["d_color"] = entry
        return
    if slot == "files.t_color":
        files["t_color"] = entry
        return
    if slot == "files.blacklist":
        files["blacklist"] = entry
        return
    if slot == "files.whitelist":
        files["whitelist"] = entry
        return
    if slot == "files.qualifier_priority":
        files["qualifier_priority"] = entry
        return
    if slot.startswith("files.c_conservation_blasts["):
        index = _slot_index(slot)
        _ensure_list_size(files["c_conservation_blasts"], index + 1, None)
        files["c_conservation_blasts"][index] = entry
        return
    if slot.startswith("files.c_depth["):
        index = _slot_index(slot)
        if not isinstance(files.get("c_depth"), list):
            files["c_depth"] = []
        _ensure_list_size(files["c_depth"], index + 1, None)
        files["c_depth"][index] = entry
        return
    if slot.startswith("files.cliInputs["):
        index = _slot_index(slot)
        _ensure_list_size(files["cliInputs"], index + 1, None)
        files["cliInputs"][index] = entry
        return
    if slot.startswith("files.linearSeqs["):
        _set_linear_seq_slot(files, slot, entry)
        return
    raise ValidationError(f"Unsupported session file slot: {slot}")


def _set_linear_seq_slot(files: dict[str, Any], slot: str, entry: dict[str, Any]) -> None:
    prefix = "files.linearSeqs["
    rest = slot[len(prefix):]
    index_text, suffix = rest.split("]", 1)
    seq_index = int(index_text)
    _ensure_linear_seq(files, seq_index)
    seq = files["linearSeqs"][seq_index]
    if suffix == ".gb":
        seq["gb"] = entry
    elif suffix == ".gff":
        seq["gff"] = entry
    elif suffix == ".fasta":
        seq["fasta"] = entry
    elif suffix == ".blast":
        seq["blast"] = entry
    elif suffix == ".depth":
        seq["depth"] = entry
    elif suffix.startswith(".depth["):
        depth_index = int(suffix[len(".depth["):-1])
        if not isinstance(seq.get("depth"), list):
            seq["depth"] = []
        _ensure_list_size(seq["depth"], depth_index + 1, None)
        seq["depth"][depth_index] = entry
    else:
        raise ValidationError(f"Unsupported linear sequence file slot: {slot}")


def _ensure_linear_seq(files: dict[str, Any], index: int) -> None:
    linear_seqs = files["linearSeqs"]
    while len(linear_seqs) <= index:
        ordinal = len(linear_seqs) + 1
        linear_seqs.append(
            {
                "uid": f"cli-seq-{ordinal}",
                "gb": None,
                "gff": None,
                "fasta": None,
                "depth": None,
                "blast": None,
                "losat_gencode": 1,
                "losat_filename": "",
                "definition": "",
                "region_record_id": "",
                "region_start": None,
                "region_end": None,
                "region_reverse": False,
            }
        )


def _slot_index(slot: str) -> int:
    left = slot.rfind("[")
    right = slot.rfind("]")
    if left < 0 or right < left:
        raise ValidationError(f"Invalid session slot index: {slot}")
    return int(slot[left + 1:right])


def _ensure_list_size(values: list[Any], size: int, fill: Any) -> None:
    while len(values) < size:
        values.append(fill)


__all__ = [
    "DiagramRunResult",
    "RenderedSvg",
    "SessionCliRequest",
    "add_session_args",
    "collect_embedded_files_from_cli_args",
    "make_rendered_svg",
    "parse_session_pre_args",
    "resolve_session_sidecar_path",
    "save_session_sidecar_if_requested",
    "strip_session_output_args",
    "validate_session_override_args",
]
