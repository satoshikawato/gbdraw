#!/usr/bin/env python
# coding: utf-8

"""GUI session JSON loading, validation, materialization, and sidecar building."""

from __future__ import annotations

import base64
import binascii
import copy
import json
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from .exceptions import ValidationError

SESSION_FORMAT = "gbdraw-session"
CURRENT_SESSION_VERSION = 27
SUPPORTED_SESSION_VERSIONS = frozenset({CURRENT_SESSION_VERSION})
DEPTH_FILE_ENCODING = "gbdraw-depth-table-v1"
DEPTH_FILE_SCHEMA = 1
JS_MAX_SAFE_INTEGER = 9_007_199_254_740_991

_DEPTH_COLUMNS = ("reference_name", "position", "depth")
_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
_SLOT_PART_RE = re.compile(r"([^\[\]]+)|\[(\d+)\]")


@dataclass(frozen=True)
class SessionFileBinding:
    argIndex: int
    slot: str
    name: str


@dataclass(frozen=True)
class SessionRunSpec:
    mode: Literal["circular", "linear"]
    args: tuple[str, ...]
    source_session: Mapping[str, Any]
    warnings: tuple[str, ...] = ()
    cli_invocation_args: tuple[str, ...] = ()
    file_bindings: tuple[SessionFileBinding, ...] = ()


@dataclass(frozen=True)
class SessionBuildContext:
    mode: Literal["circular", "linear"]
    output_prefix: str | None
    render_formats: tuple[str, ...]
    source_session: Mapping[str, Any] | None = None
    cli_invocation_args: tuple[str, ...] = ()
    file_bindings: tuple[SessionFileBinding | Mapping[str, Any], ...] = ()


def load_session(path: str | Path) -> dict[str, Any]:
    """Load and validate a gbdraw GUI session JSON file."""

    session_path = Path(path)
    try:
        payload = json.loads(session_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValidationError(f"Not a valid JSON session file: {session_path}") from exc
    except OSError as exc:
        raise ValidationError(f"Could not read session file: {session_path}") from exc
    if not isinstance(payload, dict):
        raise ValidationError("Session JSON must be an object.")
    validate_session(payload)
    return payload


def validate_session(session: Mapping[str, Any]) -> None:
    """Validate the conservative session envelope used by the CLI."""

    if not isinstance(session, Mapping):
        raise ValidationError("Session JSON must be an object.")
    if session.get("format") != SESSION_FORMAT:
        raise ValidationError("Not a gbdraw-session JSON file.")
    version = session.get("version")
    if not isinstance(version, int):
        raise ValidationError("Session version is required and must be an integer.")
    if version > CURRENT_SESSION_VERSION:
        raise ValidationError(
            f"Session version {version} is newer than this gbdraw supports "
            f"({CURRENT_SESSION_VERSION})."
        )
    if version not in SUPPORTED_SESSION_VERSIONS:
        raise ValidationError(f"Unsupported session version: {version}.")
    files = session.get("files")
    if files is None or not isinstance(files, Mapping):
        raise ValidationError("Session files are required for CLI regeneration.")


def session_mode(session: Mapping[str, Any]) -> str | None:
    """Return the declared session mode when available."""

    cli_invocation = session.get("cliInvocation")
    if isinstance(cli_invocation, Mapping):
        mode = cli_invocation.get("mode")
        if mode in {"circular", "linear"}:
            return str(mode)
    ui = session.get("ui")
    if isinstance(ui, Mapping):
        mode = ui.get("mode")
        if mode in {"circular", "linear"}:
            return str(mode)
    return None


def safe_embedded_filename(name: object, *, fallback: str = "embedded-file") -> str:
    """Return a basename-only filename safe for materializing embedded content."""

    raw_name = str(name or "").replace("\\", "/").split("/")[-1].strip()
    cleaned = _SAFE_FILENAME_RE.sub("_", raw_name).strip("._")
    return cleaned or fallback


def encode_depth_text(text: str) -> dict[str, Any] | None:
    """Encode samtools depth text using the browser session depth codec schema."""

    if not isinstance(text, str) or not text:
        return None
    crlf_count = text.count("\r\n")
    lf_count = text.count("\n")
    if crlf_count > 0 and crlf_count != lf_count:
        return None
    line_ending = "\r\n" if "\r\n" in text else "\n"
    normalized = text.replace("\r\n", "\n")
    if "\r" in normalized:
        return None
    final_newline = normalized.endswith("\n")
    lines = normalized.split("\n")
    if final_newline:
        lines.pop()
    if not lines or any(line == "" for line in lines):
        return None

    first_fields = lines[0].split("\t")
    if len(first_fields) != len(_DEPTH_COLUMNS):
        return None
    line_index = 0
    header: list[str] | None = None
    if _has_depth_header(first_fields):
        header = first_fields
        line_index = 1
    if line_index >= len(lines):
        return None

    records: list[dict[str, Any]] = []
    row_count = 0
    for line in lines[line_index:]:
        fields = line.split("\t")
        if len(fields) != len(_DEPTH_COLUMNS):
            return None
        position = _parse_positive_safe_integer(fields[1])
        depth_value = str(fields[2] or "").strip()
        if position is None or not _is_depth_text(depth_value):
            return None
        _append_depth_row(records, str(fields[0] or ""), position, depth_value)
        row_count += 1
    if row_count == 0:
        return None

    return {
        "schema": DEPTH_FILE_SCHEMA,
        "columns": list(_DEPTH_COLUMNS),
        "lineEnding": line_ending,
        "finalNewline": final_newline,
        "rowCount": row_count,
        "header": header,
        "records": records,
    }


def decode_depth_payload(payload: Mapping[str, Any]) -> str:
    """Decode a browser depth-file codec payload into TSV text."""

    if not isinstance(payload, Mapping) or payload.get("schema") != DEPTH_FILE_SCHEMA:
        raise ValidationError("Invalid embedded depth file.")
    records = payload.get("records")
    if not isinstance(records, list):
        raise ValidationError("Invalid embedded depth file records.")
    line_ending = "\r\n" if payload.get("lineEnding") == "\r\n" else "\n"
    lines: list[str] = []
    header = _decode_depth_header(payload.get("header"))
    if header is not None:
        lines.append(header)
    decoded_rows = 0
    for record in records:
        if not isinstance(record, Mapping) or not isinstance(record.get("runs"), list):
            raise ValidationError("Invalid embedded depth record.")
        reference_name = str(record.get("id") or "")
        for run in record["runs"]:
            decoded_rows += _decode_depth_run(reference_name, run, lines)

    declared_rows = payload.get("rowCount")
    if declared_rows is not None and declared_rows != decoded_rows:
        raise ValidationError("Embedded depth file row count does not match payload.")
    if not lines:
        return ""
    body = line_ending.join(lines)
    return body if payload.get("finalNewline") is False else f"{body}{line_ending}"


def serialize_file_entry(path: str | Path, *, depth: bool = False) -> dict[str, Any]:
    """Serialize a local file into the GUI session embedded-file shape."""

    file_path = Path(path)
    try:
        data = file_path.read_bytes()
    except OSError as exc:
        raise ValidationError(f"Could not read file for session embedding: {file_path}") from exc
    entry: dict[str, Any] = {
        "name": file_path.name or "file",
        "type": _guess_file_type(file_path),
        "size": len(data),
        "lastModified": int(file_path.stat().st_mtime * 1000),
    }
    if depth:
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            text = ""
        encoded_depth = encode_depth_text(text)
        if encoded_depth is not None:
            entry["encoding"] = DEPTH_FILE_ENCODING
            entry["data"] = encoded_depth
            return entry
    entry["data"] = base64.b64encode(data).decode("ascii")
    return entry


def materialize_embedded_file(
    entry: Mapping[str, Any],
    *,
    temp_dir: Path,
    role: str,
) -> Path:
    """Decode one embedded session file into temp_dir and return its path."""

    temp_dir.mkdir(parents=True, exist_ok=True)
    if not isinstance(entry, Mapping):
        raise ValidationError(f"Embedded file for {role} is missing or invalid.")
    filename = safe_embedded_filename(entry.get("name"), fallback=f"{role}.dat")
    output_path = temp_dir / f"{safe_embedded_filename(role)}-{filename}"
    output_path = _assert_under_directory(output_path, temp_dir)

    if entry.get("encoding") == DEPTH_FILE_ENCODING:
        data = entry.get("data")
        if not isinstance(data, Mapping):
            raise ValidationError(f"Embedded depth payload for {role} is malformed.")
        text = decode_depth_payload(data)
        payload_bytes = text.encode("utf-8")
    else:
        data = entry.get("data")
        if not isinstance(data, str):
            raise ValidationError(f"Embedded file for {role} has no base64 data.")
        try:
            payload_bytes = base64.b64decode(data, validate=True)
        except (binascii.Error, ValueError) as exc:
            raise ValidationError(f"Embedded file for {role} has invalid base64 data.") from exc

    declared_size = entry.get("size")
    if declared_size is not None:
        try:
            expected_size = int(declared_size)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"Embedded file for {role} has invalid size metadata.") from exc
        if expected_size >= 0 and expected_size != len(payload_bytes):
            raise ValidationError(
                f"Embedded file size mismatch for {role}: expected {expected_size}, "
                f"decoded {len(payload_bytes)}."
            )
    try:
        output_path.write_bytes(payload_bytes)
    except OSError as exc:
        raise ValidationError(f"Could not materialize embedded file for {role}.") from exc
    return output_path


def session_to_cli_args(
    session: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    """Convert a GUI/CLI session into normal CLI arguments and temp files."""

    validate_session(session)
    declared_mode = session_mode(session)
    if declared_mode and declared_mode != mode:
        raise ValidationError(
            f"Session mode is {declared_mode!r}; it cannot be used with the {mode} command."
        )

    cli_invocation = session.get("cliInvocation")
    if isinstance(cli_invocation, Mapping) and cli_invocation:
        return _session_cli_invocation_to_args(
            session,
            cli_invocation=cli_invocation,
            mode=mode,
            temp_dir=temp_dir,
            output_override=output_override,
            format_override=format_override,
        )
    return _gui_session_to_cli_args(
        session,
        mode=mode,
        temp_dir=temp_dir,
        output_override=output_override,
        format_override=format_override,
    )


def build_session_json(
    context: SessionBuildContext,
    *,
    svg_results: Sequence[tuple[str, str]],
    embedded_files: Mapping[str, Any],
    generated_at: datetime,
) -> dict[str, Any]:
    """Build a GUI-loadable session JSON payload from a CLI run."""

    if context.source_session is not None:
        payload: dict[str, Any] = _json_clone(context.source_session)
    else:
        payload = {}

    payload["format"] = SESSION_FORMAT
    payload["version"] = CURRENT_SESSION_VERSION
    payload["createdAt"] = generated_at.isoformat()
    if context.output_prefix:
        payload["title"] = Path(str(context.output_prefix)).name
    else:
        payload.setdefault("title", "gbdraw")

    config = payload.get("config")
    if not isinstance(config, dict):
        config = _minimal_config_from_cli_args(context)
        payload["config"] = config
    _update_config_prefix(config, context.output_prefix)

    ui = payload.get("ui")
    if not isinstance(ui, dict):
        ui = {}
        payload["ui"] = ui
    ui["mode"] = context.mode
    ui.setdefault("zoom", 1)
    ui.setdefault("selectedResultIndex", 0)
    ui.setdefault("canvasPan", {"x": 0, "y": 0})
    ui.setdefault("canvasPadding", {"top": 0, "right": 0, "bottom": 0, "left": 0})
    if context.mode == "circular":
        ui.setdefault("cInputType", _input_type_from_args(context.cli_invocation_args))
    else:
        ui.setdefault("lInputType", _input_type_from_args(context.cli_invocation_args))

    payload["files"] = _json_clone(embedded_files)
    payload["results"] = [
        {"name": name or f"Result {index + 1}", "content": content}
        for index, (name, content) in enumerate(svg_results)
    ]
    payload.setdefault("features", {})
    payload.setdefault("orthogroupState", {})
    payload.setdefault("losatCache", {})
    payload["cliInvocation"] = {
        "schema": 1,
        "mode": context.mode,
        "args": [str(arg) for arg in context.cli_invocation_args],
        "renderFormats": [str(fmt) for fmt in context.render_formats],
        "fileBindings": [_binding_to_json(binding) for binding in context.file_bindings],
        "generatedBy": "gbdraw",
    }
    return payload


def write_session_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    """Write a session JSON file with a same-directory temporary replacement."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    try:
        temp_path.write_text(
            json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
            encoding="utf-8",
        )
        temp_path.replace(output_path)
    except OSError as exc:
        raise ValidationError(f"Could not write session sidecar: {output_path}") from exc
    finally:
        try:
            if temp_path.exists():
                temp_path.unlink()
        except OSError:
            pass


def get_session_slot(session: Mapping[str, Any], slot: str) -> Any:
    """Resolve a slot path such as files.linearSeqs[0].gb inside a session."""

    current: Any = session
    for part in _parse_slot(slot):
        if isinstance(part, int):
            if not isinstance(current, Sequence) or isinstance(current, (str, bytes, bytearray)):
                raise ValidationError(f"Session file binding slot is not a list: {slot}")
            if part < 0 or part >= len(current):
                raise ValidationError(f"Session file binding slot index is out of range: {slot}")
            current = current[part]
        else:
            if not isinstance(current, Mapping) or part not in current:
                raise ValidationError(f"Session file binding slot is missing: {slot}")
            current = current[part]
    return current


def _session_cli_invocation_to_args(
    session: Mapping[str, Any],
    *,
    cli_invocation: Mapping[str, Any],
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    if cli_invocation.get("schema") != 1:
        raise ValidationError("Unsupported cliInvocation schema.")
    invocation_mode = cli_invocation.get("mode")
    if invocation_mode != mode:
        raise ValidationError(
            f"Session cliInvocation mode is {invocation_mode!r}; expected {mode!r}."
        )
    raw_args = cli_invocation.get("args")
    if not isinstance(raw_args, list) or not all(isinstance(arg, str) for arg in raw_args):
        raise ValidationError("Session cliInvocation args must be a string array.")

    invocation_args = [str(arg) for arg in raw_args]
    run_args = list(invocation_args)
    file_bindings = _normalize_file_bindings(cli_invocation.get("fileBindings"))
    for binding in file_bindings:
        if binding.argIndex < 0 or binding.argIndex >= len(run_args):
            raise ValidationError(
                f"cliInvocation.fileBindings argIndex {binding.argIndex} is out of range."
            )
        entry = get_session_slot(session, binding.slot)
        materialized = materialize_embedded_file(
            entry,
            temp_dir=temp_dir,
            role=f"arg{binding.argIndex}",
        )
        run_args[binding.argIndex] = str(materialized)

    run_args = _apply_option_override(run_args, "-o", "--output", output_override)
    run_args = _apply_option_override(run_args, "-f", "--format", format_override)
    invocation_args = _apply_option_override(invocation_args, "-o", "--output", output_override)
    invocation_args = _apply_option_override(invocation_args, "-f", "--format", format_override)

    return SessionRunSpec(
        mode=mode,
        args=tuple(run_args),
        source_session=session,
        cli_invocation_args=tuple(invocation_args),
        file_bindings=tuple(file_bindings),
    )


def _gui_session_to_cli_args(
    session: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    temp_dir: Path,
    output_override: str | None,
    format_override: str | None,
) -> SessionRunSpec:
    config = session.get("config")
    if not isinstance(config, Mapping):
        raise ValidationError("GUI session config is required when cliInvocation is absent.")
    files = session.get("files")
    if not isinstance(files, Mapping):
        raise ValidationError("GUI session files are required.")
    ui = session.get("ui")
    if not isinstance(ui, Mapping):
        ui = {}
    form = config.get("form") if isinstance(config.get("form"), Mapping) else {}
    adv = config.get("adv") if isinstance(config.get("adv"), Mapping) else {}

    run_args: list[str] = []
    invocation_args: list[str] = []
    bindings: list[SessionFileBinding] = []

    output_prefix = output_override or _string_or_none(form.get("prefix"))
    if output_prefix:
        _append_pair(run_args, invocation_args, "-o", output_prefix)
    _append_pair(run_args, invocation_args, "-f", format_override or "svg")
    _append_common_gui_args(run_args, invocation_args, config=config, form=form, adv=adv)

    if mode == "circular":
        _append_circular_gui_args(
            run_args,
            invocation_args,
            bindings,
            session=session,
            files=files,
            ui=ui,
            form=form,
            adv=adv,
            temp_dir=temp_dir,
        )
    else:
        _append_linear_gui_args(
            run_args,
            invocation_args,
            bindings,
            session=session,
            files=files,
            ui=ui,
            form=form,
            adv=adv,
            temp_dir=temp_dir,
        )

    return SessionRunSpec(
        mode=mode,
        args=tuple(run_args),
        source_session=session,
        cli_invocation_args=tuple(invocation_args),
        file_bindings=tuple(bindings),
    )


def _append_common_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    config: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    if _string_or_none(form.get("species")):
        _append_pair(run_args, invocation_args, "--species", str(form.get("species")))
    if _string_or_none(form.get("strain")):
        _append_pair(run_args, invocation_args, "--strain", str(form.get("strain")))
    if form.get("separate_strands") is True:
        _append_flag(run_args, invocation_args, "--separate_strands")
    features = adv.get("features")
    if isinstance(features, list) and features:
        _append_pair(run_args, invocation_args, "-k", ",".join(str(item) for item in features if item))
    for key, option in (
        ("window_size", "--window"),
        ("step_size", "--step"),
        ("nt", "--nt"),
        ("def_font_size", "--definition_font_size"),
        ("label_font_size", "--label_font_size"),
        ("block_stroke_width", "--block_stroke_width"),
        ("block_stroke_color", "--block_stroke_color"),
        ("line_stroke_width", "--line_stroke_width"),
        ("line_stroke_color", "--line_stroke_color"),
        ("axis_stroke_width", "--axis_stroke_width"),
        ("axis_stroke_color", "--axis_stroke_color"),
        ("legend_box_size", "--legend_box_size"),
        ("legend_font_size", "--legend_font_size"),
        ("scale_interval", "--scale_interval"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))
    if adv.get("resolve_overlaps") is True:
        _append_flag(run_args, invocation_args, "--resolve_overlaps")
    if adv.get("gc_content_mode") == "percent":
        _append_pair(run_args, invocation_args, "--gc_content_mode", "percent")
        for key, option in (
            ("gc_content_min_percent", "--gc_content_min_percent"),
            ("gc_content_max_percent", "--gc_content_max_percent"),
            ("gc_content_tick_interval", "--gc_content_large_tick_interval"),
            ("gc_content_small_tick_interval", "--gc_content_small_tick_interval"),
            ("gc_content_tick_font_size", "--gc_content_tick_font_size"),
        ):
            value = adv.get(key)
            if value not in (None, "", False):
                _append_pair(run_args, invocation_args, option, str(value))
        if adv.get("gc_content_show_axis") is False:
            _append_flag(run_args, invocation_args, "--hide_gc_content_axis")
        if adv.get("gc_content_show_ticks") is False:
            _append_flag(run_args, invocation_args, "--hide_gc_content_ticks")


def _append_circular_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    files: Mapping[str, Any],
    ui: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    if _string_or_none(form.get("track_type")):
        _append_pair(run_args, invocation_args, "--track_type", str(form.get("track_type")))
    if _string_or_none(form.get("legend")):
        _append_pair(run_args, invocation_args, "-l", str(form.get("legend")))
    plot_title = _string_or_none(form.get("plot_title"))
    if plot_title:
        _append_pair(run_args, invocation_args, "--plot_title", plot_title)
    if _string_or_none(adv.get("plot_title_position")):
        _append_pair(run_args, invocation_args, "--plot_title_position", str(adv.get("plot_title_position")))
    if adv.get("plot_title_font_size") not in (None, "", False):
        _append_pair(run_args, invocation_args, "--plot_title_font_size", str(adv.get("plot_title_font_size")))
    if adv.get("keep_full_definition_with_plot_title") is True:
        _append_flag(run_args, invocation_args, "--keep_full_definition_with_plot_title")
    if adv.get("center_reserved_radius") not in (None, "", False):
        _append_pair(run_args, invocation_args, "--center_reserved_radius", str(adv.get("center_reserved_radius")))
    labels_mode = str(form.get("labels_mode") or "none")
    if labels_mode == "out":
        _append_flag(run_args, invocation_args, "--labels")
    elif labels_mode == "both":
        _append_pair(run_args, invocation_args, "--labels", "both")
    if form.get("suppress_gc") is True:
        _append_flag(run_args, invocation_args, "--suppress_gc")
    if form.get("suppress_skew") is True:
        _append_flag(run_args, invocation_args, "--suppress_skew")
    if form.get("multi_record_canvas") is True:
        _append_flag(run_args, invocation_args, "--multi_record_canvas")
    for key, option in (
        ("multi_record_size_mode", "--multi_record_size_mode"),
        ("multi_record_min_radius_ratio", "--multi_record_min_radius_ratio"),
        ("multi_record_column_gap_ratio", "--multi_record_column_gap_ratio"),
        ("multi_record_row_gap_ratio", "--multi_record_row_gap_ratio"),
        ("feature_width_circular", "--feature_width"),
        ("depth_width_circular", "--depth_width"),
        ("gc_content_width_circular", "--gc_content_width"),
        ("gc_content_radius_circular", "--gc_content_radius"),
        ("gc_skew_width_circular", "--gc_skew_width"),
        ("gc_skew_radius_circular", "--gc_skew_radius"),
        ("tick_label_font_size", "--tick_label_font_size"),
        ("circular_label_spacing", "--circular_label_spacing"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))

    input_type = str(ui.get("cInputType") or "gb")
    if input_type == "gb":
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_gb",
            option="--gbk",
            temp_dir=temp_dir,
        )
    else:
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_gff",
            option="--gff",
            temp_dir=temp_dir,
        )
        _append_materialized_file_option(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot="files.c_fasta",
            option="--fasta",
            temp_dir=temp_dir,
        )

    _append_depth_gui_options(
        run_args,
        invocation_args,
        bindings,
        session=session,
        slot_prefix="files.c_depth",
        option="--depth_track",
        temp_dir=temp_dir,
        show_depth=bool(form.get("show_depth")),
    )

    conservation_blasts = files.get("c_conservation_blasts")
    if isinstance(conservation_blasts, list) and conservation_blasts:
        _append_flag(run_args, invocation_args, "--conservation_blast")
        for index, entry in enumerate(conservation_blasts):
            if entry:
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.c_conservation_blasts[{index}]",
                    temp_dir=temp_dir,
                )


def _append_linear_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    files: Mapping[str, Any],
    ui: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    if _string_or_none(form.get("scale_style")):
        _append_pair(run_args, invocation_args, "--scale_style", str(form.get("scale_style")))
    if form.get("align_center") is True:
        _append_flag(run_args, invocation_args, "--align_center")
    if form.get("show_gc") is True:
        _append_flag(run_args, invocation_args, "--show_gc")
    if form.get("show_skew") is True:
        _append_flag(run_args, invocation_args, "--show_skew")
    if form.get("normalize_length") is True:
        _append_flag(run_args, invocation_args, "--normalize_length")
    if _string_or_none(form.get("legend")) and form.get("legend") != "right":
        _append_pair(run_args, invocation_args, "-l", str(form.get("legend")))
    labels_mode = str(form.get("show_labels_linear") or "none")
    if labels_mode == "all":
        _append_flag(run_args, invocation_args, "--show_labels")
    elif labels_mode == "first":
        _append_pair(run_args, invocation_args, "--show_labels", "first")
    for key, option in (
        ("feature_height", "--feature_height"),
        ("gc_height", "--gc_height"),
        ("comparison_height", "--comparison_height"),
        ("scale_font_size", "--scale_font_size"),
        ("scale_stroke_width", "--scale_stroke_width"),
        ("scale_stroke_color", "--scale_stroke_color"),
        ("ruler_label_color", "--ruler_label_color"),
        ("pairwise_match_style", "--pairwise_match_style"),
        ("track_axis_gap", "--track_axis_gap"),
        ("label_placement", "--label_placement"),
        ("label_rendering", "--label_rendering"),
        ("label_rotation", "--label_rotation"),
        ("linear_label_spacing", "--linear_label_spacing"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))
    if _string_or_none(form.get("linear_track_layout")):
        _append_pair(run_args, invocation_args, "--track_layout", str(form.get("linear_track_layout")))
    if form.get("linear_ruler_on_axis") is True:
        _append_flag(run_args, invocation_args, "--ruler_on_axis")

    linear_seqs = files.get("linearSeqs")
    if not isinstance(linear_seqs, list) or not linear_seqs:
        raise ValidationError("Linear GUI session has no embedded sequence files.")
    input_type = str(ui.get("lInputType") or "gb")
    if input_type == "gb":
        _append_flag(run_args, invocation_args, "--gbk")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("gb"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].gb",
                    temp_dir=temp_dir,
                )
    else:
        _append_flag(run_args, invocation_args, "--gff")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("gff"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].gff",
                    temp_dir=temp_dir,
                )
        _append_flag(run_args, invocation_args, "--fasta")
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get("fasta"):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].fasta",
                    temp_dir=temp_dir,
                )
    blast_slots = [
        index for index, seq in enumerate(linear_seqs)
        if isinstance(seq, Mapping) and seq.get("blast")
    ]
    if blast_slots:
        _append_flag(run_args, invocation_args, "-b")
        for index in blast_slots:
            _append_materialized_value(
                run_args,
                invocation_args,
                bindings,
                session=session,
                slot=f"files.linearSeqs[{index}].blast",
                temp_dir=temp_dir,
            )
    if form.get("show_depth") is True:
        depth_rows = [
            _as_list(seq.get("depth") if isinstance(seq, Mapping) else None)
            for seq in linear_seqs
        ]
        track_count = max((len(row) for row in depth_rows), default=0)
        for track_index in range(track_count):
            _append_flag(run_args, invocation_args, "--depth_track")
            for record_index, row in enumerate(depth_rows):
                if track_index >= len(row) or not row[track_index]:
                    _append_value(run_args, invocation_args, "none")
                    continue
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{record_index}].depth"
                    + (f"[{track_index}]" if isinstance(linear_seqs[record_index].get("depth"), list) else ""),
                    temp_dir=temp_dir,
                )
        if track_count:
            _append_flag(run_args, invocation_args, "--show_depth")


def _append_depth_gui_options(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot_prefix: str,
    option: str,
    temp_dir: Path,
    show_depth: bool,
) -> None:
    if not show_depth:
        return
    entries = _as_list(get_session_slot(session, slot_prefix))
    for index, entry in enumerate(entries):
        if not entry:
            continue
        _append_flag(run_args, invocation_args, option)
        slot = slot_prefix if len(entries) == 1 else f"{slot_prefix}[{index}]"
        _append_materialized_value(
            run_args,
            invocation_args,
            bindings,
            session=session,
            slot=slot,
            temp_dir=temp_dir,
        )
    if entries:
        _append_flag(run_args, invocation_args, "--show_depth")


def _append_materialized_file_option(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot: str,
    option: str,
    temp_dir: Path,
) -> None:
    _append_flag(run_args, invocation_args, option)
    _append_materialized_value(
        run_args,
        invocation_args,
        bindings,
        session=session,
        slot=slot,
        temp_dir=temp_dir,
    )


def _append_materialized_value(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    slot: str,
    temp_dir: Path,
) -> None:
    entry = get_session_slot(session, slot)
    path = materialize_embedded_file(
        entry,
        temp_dir=temp_dir,
        role=slot.replace(".", "_").replace("[", "_").replace("]", ""),
    )
    arg_index = len(run_args)
    name = safe_embedded_filename(entry.get("name") if isinstance(entry, Mapping) else "")
    run_args.append(str(path))
    invocation_args.append(name)
    bindings.append(SessionFileBinding(argIndex=arg_index, slot=slot, name=name))


def _append_flag(run_args: list[str], invocation_args: list[str], option: str) -> None:
    run_args.append(str(option))
    invocation_args.append(str(option))


def _append_pair(
    run_args: list[str],
    invocation_args: list[str],
    option: str,
    value: object,
) -> None:
    run_args.extend([str(option), str(value)])
    invocation_args.extend([str(option), str(value)])


def _append_value(run_args: list[str], invocation_args: list[str], value: object) -> None:
    run_args.append(str(value))
    invocation_args.append(str(value))


def _apply_option_override(
    args: list[str],
    short_option: str,
    long_option: str,
    value: str | None,
) -> list[str]:
    if value is None:
        return list(args)
    result: list[str] = []
    replaced = False
    replace_index: int | None = None
    for index, token in enumerate(args[:-1]):
        if token in {short_option, long_option}:
            replace_index = index + 1
    for index, token in enumerate(args):
        if replace_index is not None and index == replace_index:
            result.append(str(value))
            replaced = True
        else:
            result.append(token)
    if not replaced:
        result.extend([short_option, str(value)])
    return result


def _normalize_file_bindings(value: Any) -> list[SessionFileBinding]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValidationError("cliInvocation.fileBindings must be an array.")
    bindings: list[SessionFileBinding] = []
    for item in value:
        if not isinstance(item, Mapping):
            raise ValidationError("cliInvocation.fileBindings entries must be objects.")
        try:
            arg_index = int(item.get("argIndex"))
        except (TypeError, ValueError) as exc:
            raise ValidationError("cliInvocation.fileBindings argIndex must be an integer.") from exc
        slot = str(item.get("slot") or "").strip()
        if not slot:
            raise ValidationError("cliInvocation.fileBindings slot is required.")
        name = safe_embedded_filename(item.get("name"), fallback="file")
        bindings.append(SessionFileBinding(argIndex=arg_index, slot=slot, name=name))
    return bindings


def _binding_to_json(binding: SessionFileBinding | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(binding, SessionFileBinding):
        return {
            "argIndex": binding.argIndex,
            "slot": binding.slot,
            "name": binding.name,
        }
    return {
        "argIndex": int(binding.get("argIndex", 0)),
        "slot": str(binding.get("slot", "")),
        "name": safe_embedded_filename(binding.get("name"), fallback="file"),
    }


def _parse_slot(slot: str) -> list[str | int]:
    normalized = str(slot or "").strip()
    if not normalized:
        raise ValidationError("Session slot cannot be empty.")
    parts: list[str | int] = []
    for raw_part in normalized.split("."):
        if not raw_part:
            raise ValidationError(f"Invalid session slot: {slot}")
        position = 0
        for match in _SLOT_PART_RE.finditer(raw_part):
            if match.start() != position:
                raise ValidationError(f"Invalid session slot: {slot}")
            position = match.end()
            key, index = match.groups()
            if key is not None:
                parts.append(key)
            elif index is not None:
                parts.append(int(index))
        if position != len(raw_part):
            raise ValidationError(f"Invalid session slot: {slot}")
    return parts


def _assert_under_directory(path: Path, directory: Path) -> Path:
    resolved_directory = directory.resolve()
    resolved_path = path.resolve()
    try:
        resolved_path.relative_to(resolved_directory)
    except ValueError as exc:
        raise ValidationError("Embedded filename cannot be safely materialized.") from exc
    return resolved_path


def _decode_depth_header(header: Any) -> str | None:
    if header is None:
        return None
    if (
        not isinstance(header, list)
        or len(header) != len(_DEPTH_COLUMNS)
    ):
        raise ValidationError("Invalid embedded depth file header.")
    return "\t".join(str(value if value is not None else "") for value in header)


def _decode_depth_run(reference_name: str, run: Any, lines: list[str]) -> int:
    if not isinstance(run, list) or len(run) != 4 or not isinstance(run[3], list):
        raise ValidationError("Invalid embedded depth run.")
    start, step, count, depths = run
    for value in (start, step, count):
        if not isinstance(value, int) or value <= 0 or value > JS_MAX_SAFE_INTEGER:
            raise ValidationError("Invalid embedded depth coordinates.")
    if len(depths) != count:
        raise ValidationError("Invalid embedded depth coordinates.")
    for index, depth_value in enumerate(depths):
        position = start + step * index
        if position > JS_MAX_SAFE_INTEGER:
            raise ValidationError("Invalid embedded depth coordinate overflow.")
        lines.append(f"{reference_name}\t{position}\t{'' if depth_value is None else depth_value}")
    return count


def _parse_positive_safe_integer(value: object) -> int | None:
    text = str(value or "").strip()
    if not re.fullmatch(r"[+-]?\d+", text):
        return None
    parsed = int(text)
    if parsed <= 0 or parsed > JS_MAX_SAFE_INTEGER:
        return None
    return parsed


def _is_depth_text(value: object) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    try:
        parsed = float(text)
    except ValueError:
        return False
    return math.isfinite(parsed) and parsed >= 0


def _has_depth_header(fields: Sequence[str]) -> bool:
    return (
        len(fields) >= 3
        and (_parse_positive_safe_integer(fields[1]) is None or not _is_depth_text(fields[2]))
    )


def _append_depth_row(
    records: list[dict[str, Any]],
    reference_name: str,
    position: int,
    depth_value: str,
) -> None:
    if not records or records[-1]["id"] != reference_name:
        records.append({"id": reference_name, "runs": []})
    runs = records[-1]["runs"]
    if not runs:
        runs.append([position, 1, 1, [depth_value]])
        return
    run = runs[-1]
    start, step, count, depths = run
    if count == 1:
        next_step = position - start
        if next_step > 0:
            run[1] = next_step
            run[2] = 2
            depths.append(depth_value)
            return
        runs.append([position, 1, 1, [depth_value]])
        return
    if position == start + step * count:
        run[2] = count + 1
        depths.append(depth_value)
        return
    runs.append([position, 1, 1, [depth_value]])


def _guess_file_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return "text/tab-separated-values"
    if suffix in {".txt", ".gff", ".gff3", ".fa", ".fasta", ".fna", ".gb", ".gbk", ".gbff"}:
        return "text/plain"
    return "application/octet-stream"


def _json_clone(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value))
    except (TypeError, ValueError):
        return copy.deepcopy(value)


def _minimal_config_from_cli_args(context: SessionBuildContext) -> dict[str, Any]:
    args = list(context.cli_invocation_args)
    form: dict[str, Any] = {
        "prefix": context.output_prefix or _option_value(args, "-o", "--output") or "",
        "separate_strands": "--separate_strands" in args,
    }
    adv: dict[str, Any] = {}
    if context.mode == "circular":
        form.update(
            {
                "track_type": _option_value(args, "--track_type") or "tuckin",
                "legend": _option_value(args, "-l", "--legend") or "left",
                "labels_mode": "both" if _option_value(args, "--labels") == "both" else ("out" if "--labels" in args else "none"),
                "multi_record_canvas": "--multi_record_canvas" in args,
                "suppress_gc": "--suppress_gc" in args,
                "suppress_skew": "--suppress_skew" in args,
                "show_depth": "--show_depth" in args or "--depth" in args or "--depth_track" in args,
            }
        )
        adv["plot_title_position"] = _option_value(args, "--plot_title_position") or "none"
    else:
        form.update(
            {
                "legend": _option_value(args, "-l", "--legend") or "bottom",
                "scale_style": _option_value(args, "--scale_style") or "bar",
                "linear_track_layout": _option_value(args, "--track_layout") or "middle",
                "align_center": "--align_center" in args,
                "show_gc": "--show_gc" in args,
                "show_skew": "--show_skew" in args,
                "show_depth": "--show_depth" in args or "--depth" in args or "--depth_track" in args,
            }
        )
        adv["plot_title_position"] = _option_value(args, "--plot_title_position") or "bottom"
    features = _option_value(args, "-k", "--features")
    if features:
        adv["features"] = [item for item in features.split(",") if item]
    for key, option_names in {
        "nt": ("-n", "--nt"),
        "window_size": ("-w", "--window"),
        "step_size": ("-s", "--step"),
        "def_font_size": ("--definition_font_size",),
        "label_font_size": ("--label_font_size",),
        "plot_title_font_size": ("--plot_title_font_size",),
    }.items():
        value = _option_value(args, *option_names)
        if value is not None:
            adv[key] = value
    return {
        "form": form,
        "adv": adv,
        "colors": {},
        "palette": "default",
        "rules": [],
        "qualifierPriorityRules": [],
        "filterMode": "none",
        "whitelist": [],
        "blacklistText": "",
    }


def _update_config_prefix(config: dict[str, Any], output_prefix: str | None) -> None:
    if not output_prefix:
        return
    form = config.setdefault("form", {})
    if isinstance(form, dict):
        form["prefix"] = output_prefix


def _input_type_from_args(args: Sequence[str]) -> str:
    return "gff" if "--gff" in args else "gb"


def _option_value(args: Sequence[str], *names: str) -> str | None:
    for index, token in enumerate(args[:-1]):
        if token in names:
            return str(args[index + 1])
    return None


def _string_or_none(value: object) -> str | None:
    text = str(value or "").strip()
    return text or None


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


__all__ = [
    "CURRENT_SESSION_VERSION",
    "DEPTH_FILE_ENCODING",
    "DEPTH_FILE_SCHEMA",
    "SESSION_FORMAT",
    "SUPPORTED_SESSION_VERSIONS",
    "SessionBuildContext",
    "SessionFileBinding",
    "SessionRunSpec",
    "build_session_json",
    "decode_depth_payload",
    "encode_depth_text",
    "get_session_slot",
    "load_session",
    "materialize_embedded_file",
    "safe_embedded_filename",
    "serialize_file_entry",
    "session_mode",
    "session_to_cli_args",
    "validate_session",
    "write_session_json",
]
