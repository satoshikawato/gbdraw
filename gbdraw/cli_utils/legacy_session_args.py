#!/usr/bin/env python
# coding: utf-8

"""Sessions before the canonical request (versions 27-30) replayed as ``gbdraw`` command-line arguments, and the table of retired CLI flags.

The replay writes the command line of a legacy Session and the temporary files its run reads; loading, validation and
migrations of a Session stay in :mod:`gbdraw.session_io`. Fresh CLI runs read the retired-flag table to reject a
retired flag, and the replay reads it to rewrite one.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping, NamedTuple, Sequence, cast

from gbdraw.analysis.collinearity import (
    COLLINEARITY_MERGE_ORIENTATIONS,
    COLLINEARITY_MERGE_ORIENTATION_SPELLINGS,
    retired_collinearity_anchor_mode_replacement,
)
from gbdraw.config.models.root import read_stored_font_size
from gbdraw.contract_shapes import read_spelling
from gbdraw.core.js_values import js_string, js_text, js_truthy
from gbdraw.definition_line_styles import DEFINITION_LINE_KINDS
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import COMPARISON_SEQUENCE_SEPARATOR
from gbdraw.io.colors import current_palette_name, packaged_color_palettes
from gbdraw.io.filenames import safe_embedded_filename, unique_filenames
from gbdraw.io.input_text import GZIP_ERRORS, open_input_text
from gbdraw.io.regions import region_text
from gbdraw.render.formats import normalize_format_token
from gbdraw.session_io import (
    CANONICAL_SESSION_MIN_VERSION,
    CURRENT_AUTHORITY_SESSION_MIN_VERSION,
    SessionFileBinding,
    SessionRunSpec,
    _mapping_list,
    _nonnegative_integer,
    _option_all_values,
    _option_value,
    _with_resolved_override_colors,
    get_session_slot,
    materialize_embedded_file,
    migrate_session_flat_draft,
    mode_split_palette_colors,
    session_mode,
    validate_session,
)
from gbdraw.session_legacy_values import (
    LEGACY_FEATURE_SHAPES,
    LEGACY_FEATURE_TYPES,
    LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES,
    LEGACY_SCALAR_AXIS_LABELS,
)
from gbdraw.tracks.circular import _LEGACY_SPACING_PARAM, _migrate_legacy_circular_track_slot_spec


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
    if int(session.get("version", 0)) >= CANONICAL_SESSION_MIN_VERSION:
        raise ValidationError(
            "Canonical renderRequest sessions cannot be replayed through "
            "legacy CLI arguments."
        )
    if format_override is not None:
        format_override = ",".join(
            normalize_format_token(value) for value in format_override.split(",")
        )
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


@dataclass(frozen=True)
class RetiredCliOption:
    """One retired CLI flag (design D4).

    Fresh runs reject the flag and name ``replacement``. Legacy session argv
    is rewritten before replay: ``renamed_to`` keeps the value, and
    ``value_rewrites`` replaces flag and value with zero or more tokens. A flag
    with neither has no legacy rewrite.
    """

    option: str
    modes: tuple[Literal["circular", "linear"], ...]
    replacement: str
    renamed_to: str | None = None
    value_rewrites: Mapping[str, tuple[str, ...]] | None = None
    # argparse nargs of the retired flag, so a fresh run reports the flag
    # itself rather than its extra values. A rename keeps every value.
    nargs: str | None = None

    def rewrite(self, value: str) -> tuple[str, ...] | None:
        if self.renamed_to is not None:
            return (self.renamed_to, value)
        if self.value_rewrites is not None:
            return self.value_rewrites.get(str(value).strip().lower())
        return None

    def message(self, value: object = None) -> str:
        replacement = self.replacement
        if self.value_rewrites is not None and isinstance(value, str):
            tokens = self.value_rewrites.get(value.strip().lower())
            if tokens is not None:
                replacement = (
                    "use " + " ".join(tokens)
                    if tokens
                    else "omit it (no protein comparison)"
                )
        return f"{self.option} was retired; {replacement}."


def _losatp_mode_rewrites() -> dict[str, tuple[str, ...]]:
    from gbdraw.api.options import LOSATP_MODE_WIRE

    rewrites: dict[str, tuple[str, ...]] = {"none": ()}
    for typed_mode, wire_mode in LOSATP_MODE_WIRE.items():
        if typed_mode != "none":
            rewrites[wire_mode] = ("--losat", "losatp", "--losatp_mode", typed_mode)
    return rewrites


# The one old -> new table for retired CLI flags (design 3.5). Both the
# fresh-run rejection and the legacy session argv rewrite read it.
RETIRED_CLI_OPTIONS: Mapping[str, RetiredCliOption] = {
    item.option: item
    for item in (
        RetiredCliOption(
            "--protein_blastp_mode",
            ("linear",),
            "use --losat losatp --losatp_mode {similarity_groups,collinear,pairwise}",
            value_rewrites=_losatp_mode_rewrites(),
        ),
        RetiredCliOption("--losatp_bin", ("linear",), "use --losat_bin", "--losat_bin"),
        RetiredCliOption(
            "--ncbi_blastp_bin", ("linear",), "use --ncbi_blast_bin", "--ncbi_blast_bin"
        ),
        RetiredCliOption(
            "--losatp_threads", ("linear",), "use --losat_threads", "--losat_threads"
        ),
        RetiredCliOption(
            "--protein_blastp_max_hits",
            ("linear",),
            "use --losatp_max_hits",
            "--losatp_max_hits",
        ),
        RetiredCliOption(
            "--protein_blastp_candidate_limit",
            ("linear",),
            "use --losatp_max_target_seqs",
            "--losatp_max_target_seqs",
        ),
        RetiredCliOption(
            "--align_orthogroup_feature",
            ("linear",),
            "use --similarity_alignment_feature",
            "--similarity_alignment_feature",
        ),
        RetiredCliOption(
            "--protein_blastp_output",
            ("linear",),
            "use --losat_output_dir DIR, which writes DIR/losatp.raw.tsv",
        ),
        RetiredCliOption(
            "--conservation_fasta",
            ("circular",),
            "use --conservation_sequence (FASTA, GenBank, or DDBJ)",
            "--conservation_sequence",
            nargs="+",
        ),
    )
}


def _canonicalize_legacy_session_cli_args(
    args: Sequence[str],
    *,
    mode: Literal["circular", "linear"],
) -> tuple[list[str], dict[int, int]]:
    replacements = {
        "--depth": "--depth_track",
        "--depth_tick_interval": "--depth_large_tick_interval",
        "--gc_content_tick_interval": "--gc_content_large_tick_interval",
        "--feature_table": "--feature_visibility_table",
        "--annotation-table": "--annotation_table",
        "--losatp-bin": "--losatp_bin",
        "--ncbi-blastp-bin": "--ncbi_blastp_bin",
        "--losatp-threads": "--losatp_threads",
        "--protein-blastp-mode": "--protein_blastp_mode",
        "--protein-blastp-max-hits": "--protein_blastp_max_hits",
        "--protein-blastp-candidate-limit": "--protein_blastp_candidate_limit",
        "--align-orthogroup-feature": "--align_orthogroup_feature",
        "--collinear-unit-mode": "--collinear_unit_mode",
        "--collinear-search-scope": "--collinear_search_scope",
        "--collinear-min-anchors": "--collinear_min_anchors",
        "--collinear_max_gene_gap": "--collinear_max_unit_gap",
        "--collinear-max-unit-gap": "--collinear_max_unit_gap",
        "--collinear-max-gene-gap": "--collinear_max_unit_gap",
        "--collinear-max-diagonal-drift": "--collinear_max_diagonal_drift",
        "--collinear-max-conflicts-in-merge-gap": "--collinear_max_conflicts_in_merge_gap",
        "--collinear-max-paralog-links-per-orthogroup": (
            "--collinear_max_paralog_links_per_orthogroup"
        ),
        "--collinear-color-mode": "--collinear_color_mode",
        "--keep-definition-left-aligned": "--keep_definition_left_aligned",
        "--pairwise-match-style": "--pairwise_match_style",
        "--definition-line-style": "--definition_line_style",
        "--record-subtitle": "--record_subtitle",
        "--circular-track-slot": "--circular_track_slot",
    }
    if mode == "circular":
        replacements.update(
            {
                "--suppress_gc": "--no-gc",
                "--suppress_skew": "--no-skew",
            }
        )
    else:
        replacements.update(
            {
                "--show_gc": "--gc",
                "--show_skew": "--skew",
            }
        )

    canonical_args: list[str] = []
    source_to_canonical_index: dict[int, int] = {}

    def extend(option_index: int, value_index: int, tokens: Sequence[str]) -> None:
        if not tokens:
            return
        source_to_canonical_index[option_index] = len(canonical_args)
        source_to_canonical_index[value_index] = len(canonical_args) + len(tokens) - 1
        canonical_args.extend(tokens)

    pending_retired: tuple[RetiredCliOption, int] | None = None
    for source_index, raw_token in enumerate(args):
        token = str(raw_token)
        if pending_retired is not None:
            retired, option_index = pending_retired
            pending_retired = None
            rewritten = retired.rewrite(token)
            extend(
                option_index,
                source_index,
                (retired.option, token) if rewritten is None else rewritten,
            )
            continue
        if token == "--show_depth":
            continue
        option, separator, inline_value = (
            token.partition("=")
            if token.startswith("--")
            else (token, "", "")
        )
        option = replacements.get(option, option)
        retired_option = RETIRED_CLI_OPTIONS.get(option)
        if retired_option is not None and mode in retired_option.modes:
            if not separator:
                pending_retired = (retired_option, source_index)
                continue
            rewritten = retired_option.rewrite(inline_value)
            if rewritten is None:
                rewritten = (f"{option}={inline_value}",)
            elif retired_option.renamed_to is not None:
                rewritten = (f"{retired_option.renamed_to}={inline_value}",)
            extend(source_index, source_index, rewritten)
            continue
        token = _canonical_legacy_token(
            option, separator, inline_value, canonical_args=canonical_args, mode=mode
        )
        source_to_canonical_index[source_index] = len(canonical_args)
        canonical_args.append(token)
    if pending_retired is not None:
        retired, option_index = pending_retired
        source_to_canonical_index[option_index] = len(canonical_args)
        canonical_args.append(retired.option)
    return canonical_args, source_to_canonical_index


def _canonical_legacy_token(
    option: str,
    separator: str,
    inline_value: str,
    *,
    canonical_args: Sequence[str],
    mode: Literal["circular", "linear"],
) -> str:
    """One token of a legacy argv that no retired flag consumes, with its old value spellings rewritten: the value of
    ``--option=value``, or a bare token read as the value of the option before it."""
    if separator:
        if mode == "circular" and option == "--circular_track_slot":
            inline_value = _migrate_legacy_circular_track_slot_spec(inline_value)
        if option == "--multi_record_size_mode":
            inline_value = LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES.get(inline_value, inline_value)
        elif mode == "linear" and option == "--label_placement":
            inline_value = {
                "on_feature": "above_feature",
            }.get(inline_value, inline_value)
        elif mode == "linear" and option == "--track_layout":
            inline_value = {
                "spreadout": "above",
                "tuckin": "below",
            }.get(inline_value, inline_value)
        token = f"{option}={inline_value}"
    else:
        token = option
        if canonical_args:
            previous_option = canonical_args[-1].partition("=")[0]
            if previous_option == "--multi_record_size_mode" and token in LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES:
                token = LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES[token]
            elif (
                mode == "linear"
                and previous_option == "--label_placement"
                and token == "on_feature"
            ):
                token = "above_feature"
            elif mode == "linear" and previous_option == "--track_layout":
                token = {"spreadout": "above", "tuckin": "below"}.get(token, token)
            elif mode == "circular" and previous_option == "--circular_track_slot":
                token = _migrate_legacy_circular_track_slot_spec(token)
    return token


def canonicalize_cli_invocation(
    args: Sequence[str],
    file_bindings: Sequence[SessionFileBinding],
    *,
    mode: Literal["circular", "linear"],
) -> tuple[list[str], list[SessionFileBinding]]:
    """Rewrite retired CLI flags in a recorded invocation and remap its bindings.

    Legacy session replay and the Gallery session refresh both use this, so a
    recorded ``cliInvocation`` reaches the current flag names the same way.
    """

    canonical_args, index_map = _canonicalize_legacy_session_cli_args(args, mode=mode)
    remapped: list[SessionFileBinding] = []
    for binding in file_bindings:
        if binding.argIndex not in index_map:
            raise ValidationError(
                "cliInvocation.fileBindings cannot reference a removed legacy CLI flag."
            )
        remapped.append(
            SessionFileBinding(
                argIndex=index_map[binding.argIndex],
                slot=binding.slot,
                name=binding.name,
            )
        )
    return canonical_args, remapped


def _keep_recorded_collinear_inference(args: list[str]) -> list[str]:
    """A recorded Collinear command without the inference flag ran with inference on.

    The CLI default became off; the replay keeps what the command drew.
    """

    if any(arg.split("=", 1)[0] == "--collinear_infer_orthogroups" for arg in args):
        return args
    collinear = any(
        arg == "--losatp_mode=collinear"
        or (arg == "--losatp_mode" and index + 1 < len(args) and args[index + 1] == "collinear")
        for index, arg in enumerate(args)
    )
    return [*args, "--collinear_infer_orthogroups", "on"] if collinear else args


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

    session_version = int(session.get("version", 0))
    migrate_legacy_cli = session_version < CURRENT_AUTHORITY_SESSION_MIN_VERSION
    _restore_cli_table_paths(
        session,
        run_args,
        temp_dir=temp_dir,
        migrate_legacy_cli=migrate_legacy_cli,
    )

    if migrate_legacy_cli:
        run_args, _ = _canonicalize_legacy_session_cli_args(run_args, mode=mode)
        invocation_args, file_bindings = canonicalize_cli_invocation(
            invocation_args,
            file_bindings,
            mode=mode,
        )

    run_args = _keep_recorded_collinear_inference(run_args)
    run_args = _apply_option_override(run_args, "-o", "--output", output_override)
    run_args = _apply_option_override(run_args, "-f", "--format", format_override)
    invocation_args = _apply_option_override(invocation_args, "-o", "--output", output_override)
    invocation_args = _apply_option_override(invocation_args, "-f", "--format", format_override)
    run_args = migrate_legacy_repeat_feature_shape_args(
        run_args,
        session_version=session_version,
    )
    invocation_args = migrate_legacy_repeat_feature_shape_args(
        invocation_args,
        session_version=session_version,
    )

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
    form_value = config.get("form")
    form = form_value if isinstance(form_value, Mapping) else {}
    adv_value = config.get("adv")
    adv = dict(adv_value) if isinstance(adv_value, Mapping) else {}
    # A stored font size <= 0 replays as Auto (no option), as Web Load reads it (OV-1282).
    for key in [key for key in adv if str(key).endswith("font_size")]:
        adv[key] = read_stored_font_size(adv[key], path=f"config.adv.{key}")
    if int(session.get("version", 0)) <= 30:
        effective_features = adv.get("features")
        if not isinstance(effective_features, list):
            effective_features = list(LEGACY_FEATURE_TYPES)
        legacy_shapes = [item for item in LEGACY_FEATURE_SHAPES.items() if item[0] in effective_features]
        if legacy_shapes:
            feature_shapes = dict(adv.get("feature_shapes") or {})
            for feature_type, shape in legacy_shapes:
                feature_shapes.setdefault(feature_type, shape)
            adv["feature_shapes"] = feature_shapes

    run_args: list[str] = []
    invocation_args: list[str] = []
    bindings: list[SessionFileBinding] = []

    output_prefix = output_override or _string_or_none(form.get("prefix"))
    if output_prefix:
        _append_pair(run_args, invocation_args, "-o", output_prefix)
    _append_pair(run_args, invocation_args, "-f", format_override or "svg")
    _append_common_gui_args(run_args, invocation_args, form=form, adv=adv)
    _append_gui_color_args(run_args, invocation_args, session=session, temp_dir=temp_dir)

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
            config=config,
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


def _restore_cli_table_paths(
    session: Mapping[str, Any],
    run_args: list[str],
    *,
    temp_dir: Path,
    migrate_legacy_cli: bool = False,
) -> None:
    files = session.get("files")
    if not isinstance(files, Mapping):
        return
    cli_tables = files.get("cliTables")
    if not isinstance(cli_tables, list):
        return

    for table_entry in cli_tables:
        if not isinstance(table_entry, Mapping):
            continue
        try:
            arg_index = int(cast(Any, table_entry.get("argIndex")))
        except (TypeError, ValueError) as exc:
            raise ValidationError("files.cliTables argIndex must be an integer.") from exc
        if arg_index < 0 or arg_index >= len(run_args):
            raise ValidationError("files.cliTables argIndex is out of range.")

        table_path = Path(str(run_args[arg_index]))
        if not table_path.is_file():
            table_slot = str(table_entry.get("slot") or "").strip()
            if not table_slot:
                raise ValidationError("files.cliTables slot is required.")
            table_file_entry = get_session_slot(session, table_slot)
            table_path = materialize_embedded_file(
                table_file_entry,
                temp_dir=temp_dir,
                role=f"arg{arg_index}",
            )
            run_args[arg_index] = str(table_path)

        table_kind = str(table_entry.get("kind") or "").strip()
        preceding_option = (
            str(run_args[arg_index - 1]) if arg_index > 0 else ""
        )
        if (
            migrate_legacy_cli
            and (
                table_kind == "circular_track"
                or preceding_option == "--circular_track_table"
            )
        ):
            _migrate_legacy_circular_track_table(table_path)

        dependencies = table_entry.get("dependencies")
        if not isinstance(dependencies, list) or not dependencies:
            continue
        # A cell's files in the order the writer listed them: one, or the
        # ``;`` separated files of a ring's comparison genome.
        replacements: dict[tuple[int, str], list[str]] = {}
        for dependency in dependencies:
            if not isinstance(dependency, Mapping):
                continue
            try:
                row_index = int(cast(Any, dependency.get("rowIndex")))
            except (TypeError, ValueError) as exc:
                raise ValidationError("files.cliTables dependencies rowIndex must be an integer.") from exc
            column = str(dependency.get("column") or "").strip()
            slot = str(dependency.get("slot") or "").strip()
            if not column or not slot:
                raise ValidationError("files.cliTables dependency entries require column and slot.")
            dependency_entry = get_session_slot(session, slot)
            materialized = materialize_embedded_file(
                dependency_entry,
                temp_dir=temp_dir,
                role=slot.replace(".", "_").replace("[", "_").replace("]", ""),
            )
            replacements.setdefault((row_index, column), []).append(
                _relative_path_for_table(materialized, table_path.parent)
            )
        if replacements:
            _rewrite_tsv_path_cells(table_path, replacements)


def _relative_path_for_table(path: Path, table_dir: Path) -> str:
    try:
        return os.path.relpath(str(path), str(table_dir))
    except ValueError:
        return str(path)


def _read_restored_tsv_rows(table_path: Path, what: str) -> list[list[str]]:
    """The rows of a restored CLI table, read through the gzip-aware input reader.

    A gzip table reads as the plain one, under ``MAX_EXPANDED_INPUT_BYTES``, whose
    refusal names the file; the rewrite after it writes the plain table.
    """

    try:
        with open_input_text(table_path, encoding="utf-8-sig", newline="") as handle:
            return list(csv.reader(handle, delimiter="\t"))
    except (OSError, UnicodeError, *GZIP_ERRORS) as exc:
        raise ValidationError(f"Could not read restored {what}: {table_path}") from exc


def _rewrite_tsv_path_cells(
    table_path: Path,
    replacements: Mapping[tuple[int, str], Sequence[str]],
) -> None:
    rows = _read_restored_tsv_rows(table_path, "TSV table")

    header: list[str] | None = None
    output_rows: list[list[str]] = []
    data_row_index = 0
    for cells in rows:
        if not cells or all(str(cell).strip() == "" for cell in cells):
            continue
        if header is None:
            header = [str(cell).strip() for cell in cells]
            output_rows.append(header)
            continue
        values = [str(cell) for cell in cells]
        while len(values) < len(header):
            values.append("")
        for (target_row_index, column), replacement in replacements.items():
            if target_row_index != data_row_index or column not in header:
                continue
            values[header.index(column)] = COMPARISON_SEQUENCE_SEPARATOR.join(replacement)
        output_rows.append(values[: len(header)])
        data_row_index += 1

    if header is None:
        raise ValidationError(f"Restored TSV table has no header row: {table_path}")
    try:
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(output_rows)
    except OSError as exc:
        raise ValidationError(f"Could not rewrite restored TSV table: {table_path}") from exc


def _migrate_legacy_circular_track_table(table_path: Path) -> None:
    """Move a persisted Circular table's spacing column to reader-only params."""

    rows = _read_restored_tsv_rows(table_path, "Circular track table")
    if not rows:
        return

    header_index = next(
        (
            index
            for index, cells in enumerate(rows)
            if cells and any(str(cell).strip() for cell in cells)
        ),
        None,
    )
    if header_index is None:
        return
    header = [str(cell).strip() for cell in rows[header_index]]
    if "spacing" not in header:
        return

    spacing_index = header.index("spacing")
    if "params" not in header:
        header.append("params")
    params_index = header.index("params")
    migrated_rows: list[list[str]] = []
    for index, cells in enumerate(rows):
        if index < header_index:
            migrated_rows.append([str(cell) for cell in cells])
            continue
        if index == header_index:
            migrated_rows.append(
                [column for column in header if column != "spacing"]
            )
            continue
        values = [str(cell) for cell in cells]
        while len(values) < len(header):
            values.append("")
        values = values[: len(header)]
        spacing = values[spacing_index].strip()
        if spacing:
            legacy_param = f"{_LEGACY_SPACING_PARAM}={spacing}"
            existing_params = values[params_index].strip()
            values[params_index] = (
                f"{existing_params},{legacy_param}"
                if existing_params
                else legacy_param
            )
        migrated_rows.append(
            [
                value
                for column, value in zip(header, values, strict=True)
                if column != "spacing"
            ]
        )
    try:
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(migrated_rows)
    except OSError as exc:
        raise ValidationError(
            f"Could not migrate restored Circular track table: {table_path}"
        ) from exc


def _append_common_gui_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    if _string_or_none(form.get("species")):
        _append_pair(run_args, invocation_args, "--species", str(form.get("species")))
    if _string_or_none(form.get("strain")):
        _append_pair(run_args, invocation_args, "--strain", str(form.get("strain")))
    if form.get("separate_strands") is True:
        _append_flag(run_args, invocation_args, "--separate_strands")
    if form.get("show_scale") is False:
        _append_flag(run_args, invocation_args, "--hide_scale")
    # The figure size as the drawing's form holds it ("180mm"); empty and "auto" are Auto, which writes no flag.
    for key, option in (("figure_width", "--figure_width"), ("figure_height", "--figure_height")):
        length = _string_or_none(form.get(key))
        if length and length.strip().lower() != "auto":
            _append_pair(run_args, invocation_args, option, length)
    features = adv.get("features")
    if isinstance(features, list) and features:
        _append_pair(run_args, invocation_args, "-k", ",".join(str(item) for item in features if item))
    feature_shapes = adv.get("feature_shapes")
    if isinstance(feature_shapes, Mapping):
        for feature_type, rendering in feature_shapes.items():
            _append_pair(
                run_args,
                invocation_args,
                "--feature_shape",
                f"{str(feature_type).strip()}={str(rendering).strip().lower()}",
            )
    for key, option in (
        ("window_size", "--window"),
        ("step_size", "--step"),
        ("nt", "--nt"),
        ("def_font_size", "--definition_font_size"),
        ("label_font_size", "--label_font_size"),
        ("arrow_head_length_ratio", "--arrow_head_length_ratio"),
        ("arrow_shaft_width_ratio", "--arrow_shaft_width_ratio"),
        ("block_stroke_width", "--block_stroke_width"),
        ("block_stroke_color", "--block_stroke_color"),
        ("line_stroke_width", "--line_stroke_width"),
        ("line_stroke_color", "--line_stroke_color"),
        ("axis_stroke_width", "--axis_stroke_width"),
        ("axis_stroke_color", "--axis_stroke_color"),
        ("legend_box_size", "--legend_box_size"),
        ("legend_font_size", "--legend_font_size"),
        ("scale_interval", "--scale_interval"),
        ("font_family", "--font_family"),
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
    # A draft without an axis label style drew numbered ticks (D-24); the CLI draws the range label by default.
    for shown, key, option in (
        (form.get("show_depth") is True, "depth_axis_labels", "--depth_axis_labels"),
        (adv.get("gc_content_mode") == "percent", "gc_content_axis_labels", "--gc_content_axis_labels"),
    ):
        if shown:
            _append_pair(run_args, invocation_args, option, str(adv.get(key) or LEGACY_SCALAR_AXIS_LABELS))


def _gui_color_table_cell(value: object) -> str:
    """One cell of a replay color table: tabs and line breaks read as a space."""

    return re.sub(r"[\t\r\n]+", " ", js_string(value)).strip()


def _append_gui_color_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    session: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    """Pass a Session 27-30's draft colors: ``-p``, ``-d`` and ``-t`` (OV-2055, Owner R17-22).

    The draft is read as Web Load reads it before a Generate: after the
    Session's draft migrations (:func:`migrate_session_flat_draft`, whose
    ``hash`` rules take their source hashes, with the sources read again when
    a rule needs drawn hashes) and with ``colorsAreOverrides`` resolved over
    the palette (:func:`mode_split_palette_colors`). The palette goes as
    ``-p``, every draft color as a ``-d`` row (a row equal to the palette's
    color draws the same), and each specific color rule with its feature
    type, qualifier, value and color as a ``-t`` row, in the CLI's table
    formats. The files are named as the Web Generate names them.
    """

    from gbdraw.session_migration import read_request_less_source_features

    migrated = migrate_session_flat_draft(session, source_features=read_request_less_source_features(session))
    config = migrated.session.get("config")
    if not isinstance(config, Mapping):
        return
    config = _with_resolved_override_colors(config, mode_split_palette_colors(config))
    palette = current_palette_name(js_text(config.get("palette")))
    if palette and palette != "default" and palette in packaged_color_palettes():
        _append_pair(run_args, invocation_args, "-p", palette)
    colors = config.get("colors")
    color_rows = [
        [_gui_color_table_cell(key), _gui_color_table_cell(color)]
        for key, color in (colors.items() if isinstance(colors, Mapping) else ())
    ]
    rules = config.get("rules")
    rule_rows = [
        [_gui_color_table_cell(rule.get(field)) for field in ("feat", "qual", "val", "color", "cap")]
        for rule in (rules if isinstance(rules, list) else ())
        if isinstance(rule, Mapping)
    ]
    for option, name, rows, required in (
        ("-d", "default-colors.tsv", color_rows, 2),
        ("-t", "specific-colors.tsv", rule_rows, 4),
    ):
        rows = [row for row in rows if all(row[:required])]
        if not rows:
            continue
        path = temp_dir / f"session-{name}"
        path.write_text("".join("\t".join(row) + "\n" for row in rows), encoding="utf-8")
        run_args.extend([option, str(path)])
        invocation_args.extend([option, name])


def _append_gui_plot_title_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    positions: frozenset[str] | None,
) -> None:
    """The draft's plot title, its position (one of ``positions``, if given) and font size."""

    plot_title = _string_or_none(form.get("plot_title"))
    if plot_title:
        _append_pair(run_args, invocation_args, "--plot_title", plot_title)
    position = _string_or_none(adv.get("plot_title_position"))
    if position and (positions is None or position in positions):
        _append_pair(run_args, invocation_args, "--plot_title_position", position)
    if adv.get("plot_title_font_size") not in (None, "", False):
        _append_pair(run_args, invocation_args, "--plot_title_font_size", str(adv.get("plot_title_font_size")))


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
    _append_gui_plot_title_args(run_args, invocation_args, form=form, adv=adv, positions=None)
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
        _append_flag(run_args, invocation_args, "--no-gc")
    if form.get("suppress_skew") is True:
        _append_flag(run_args, invocation_args, "--no-skew")
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
        ("circular_label_placement", "--label_placement"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            if key == "multi_record_size_mode":
                value = LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES.get(str(value).strip().lower(), value)
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
    config: Mapping[str, Any],
    form: Mapping[str, Any],
    adv: Mapping[str, Any],
    temp_dir: Path,
) -> None:
    if _string_or_none(form.get("scale_style")):
        _append_pair(run_args, invocation_args, "--scale_style", str(form.get("scale_style")))
    if form.get("align_center") is True:
        _append_flag(run_args, invocation_args, "--align_center")
    if form.get("show_gc") is True:
        _append_flag(run_args, invocation_args, "--gc")
    if form.get("show_skew") is True:
        _append_flag(run_args, invocation_args, "--skew")
    if form.get("normalize_length") is True:
        _append_flag(run_args, invocation_args, "--normalize_length")
    if _string_or_none(form.get("legend")) and form.get("legend") != "right":
        _append_pair(run_args, invocation_args, "-l", str(form.get("legend")))
    labels_mode = str(form.get("show_labels_linear") or "none")
    if labels_mode == "all":
        _append_flag(run_args, invocation_args, "--show_labels")
    elif labels_mode == "first":
        _append_pair(run_args, invocation_args, "--show_labels", "first")
    elif labels_mode == "orthogroup_top":
        _append_pair(run_args, invocation_args, "--show_labels", "orthogroup_top")
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
            if key == "label_placement" and str(value).strip().lower() == "on_feature":
                value = "above_feature"
            _append_pair(run_args, invocation_args, option, str(value))
    if _string_or_none(form.get("linear_track_layout")):
        layout = str(form.get("linear_track_layout")).strip().lower()
        layout = {"spreadout": "above", "tuckin": "below"}.get(layout, layout)
        _append_pair(run_args, invocation_args, "--track_layout", layout)
    if form.get("linear_ruler_on_axis") is True:
        _append_flag(run_args, invocation_args, "--ruler_on_axis")
    wrap_length = _string_or_none(form.get("linear_wrap_length"))
    if wrap_length and wrap_length != "none":
        _append_pair(run_args, invocation_args, "--wrap_length", wrap_length)
    # OV-2521: a Linear Generate draws the plot title too.
    from gbdraw.api.options import LINEAR_PLOT_TITLE_POSITIONS

    _append_gui_plot_title_args(
        run_args, invocation_args, form=form, adv=adv, positions=LINEAR_PLOT_TITLE_POSITIONS
    )
    _append_linear_definition_line_style_args(run_args, invocation_args, adv)

    linear_seqs = files.get("linearSeqs")
    if not isinstance(linear_seqs, list) or not linear_seqs:
        raise ValidationError("Linear GUI session has no embedded sequence files.")
    input_type = str(ui.get("lInputType") or "gb")
    _append_linear_record_files(
        run_args,
        invocation_args,
        bindings,
        session=session,
        linear_seqs=linear_seqs,
        gff=input_type != "gb",
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
    elif _gui_linear_losat_program(config, adv) == "blastp":
        _append_linear_gui_blastp_args(
            run_args,
            invocation_args,
            session=session,
            config=config,
            adv=adv,
        )
    _append_linear_gui_sequence_options(
        run_args,
        invocation_args,
        linear_seqs=linear_seqs,
        card_files=_linear_card_files(run_args, bindings, "gb" if input_type == "gb" else "gff"),
    )
    if form.get("show_depth") is True:
        _append_linear_gui_depth_tracks(
            run_args,
            invocation_args,
            bindings,
            session=session,
            linear_seqs=linear_seqs,
            temp_dir=temp_dir,
        )


def _append_linear_gui_depth_tracks(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    linear_seqs: list[Any],
    temp_dir: Path,
) -> None:
    """One ``--depth_track`` per depth track, with one file (or ``none``) per Linear record."""
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


def _append_linear_definition_line_style_args(
    run_args: list[str],
    invocation_args: list[str],
    adv: Mapping[str, Any],
) -> None:
    styles = adv.get("linear_definition_line_styles")
    if not isinstance(styles, Mapping):
        return
    for line_kind in DEFINITION_LINE_KINDS:
        raw_style = styles.get(line_kind)
        if not isinstance(raw_style, Mapping):
            continue
        parts: list[str] = []
        font_size = read_stored_font_size(
            raw_style.get("font_size"),
            path=f"config.adv.linear_definition_line_styles.{line_kind}.font_size",
        )
        if font_size not in (None, "", False):
            parts.append(f"size={font_size}")
        font_weight = str(raw_style.get("font_weight") or "").strip()
        if font_weight.lower() in {"auto", "none", "null", "default", "normal"}:
            font_weight = ""
        if font_weight:
            parts.append(f"weight={font_weight}")
        fill = str(raw_style.get("fill") or "").strip()
        if fill:
            parts.append(f"color={fill}")
        if parts:
            _append_pair(run_args, invocation_args, "--definition_line_style", f"{line_kind}:{','.join(parts)}")


def _append_linear_record_files(
    run_args: list[str],
    invocation_args: list[str],
    bindings: list[SessionFileBinding],
    *,
    session: Mapping[str, Any],
    linear_seqs: Sequence[Any],
    gff: bool,
    temp_dir: Path,
) -> None:
    """The record files of the Linear cards: ``--gbk``, or ``--gff`` and ``--fasta``, one per card in card order."""

    first = len(bindings)
    for option, role in ((("--gff", "gff"), ("--fasta", "fasta")) if gff else (("--gbk", "gb"),)):
        _append_flag(run_args, invocation_args, option)
        for index, seq in enumerate(linear_seqs):
            if isinstance(seq, Mapping) and seq.get(role):
                _append_materialized_value(
                    run_args,
                    invocation_args,
                    bindings,
                    session=session,
                    slot=f"files.linearSeqs[{index}].{role}",
                    temp_dir=temp_dir,
                )
    _name_record_files_apart(invocation_args, bindings, first)


def _name_record_files_apart(
    invocation_args: list[str], bindings: list[SessionFileBinding], first: int
) -> None:
    """Give the record files from ``bindings[first:]`` distinct names on the shown command line (``unique_filenames``).

    A file region names its card's file (``<file>:#1:<span>``); two cards of one file name would both match it.
    """

    record_files = bindings[first:]
    for at, (binding, name) in enumerate(
        zip(record_files, unique_filenames([binding.name for binding in record_files])), start=first
    ):
        if name != binding.name:
            invocation_args[binding.argIndex] = name
            bindings[at] = SessionFileBinding(argIndex=binding.argIndex, slot=binding.slot, name=name)


_LINEAR_CARD_SLOT_RE = re.compile(r"files\.linearSeqs\[(?P<index>\d+)\]\.(?P<role>gb|gff)")


class _CardFile(NamedTuple):
    """A Linear card's record file: its name on the command line run here and on the command line shown."""

    run_name: str
    shown_name: str


def _linear_card_files(
    run_args: Sequence[str], bindings: Sequence[SessionFileBinding], role: str
) -> dict[int, _CardFile]:
    """The record file of each Linear card (its ``gb``, or its ``gff`` beside a FASTA), by card index."""

    files: dict[int, _CardFile] = {}
    for binding in bindings:
        match = _LINEAR_CARD_SLOT_RE.fullmatch(binding.slot)
        if match and match.group("role") == role:
            files[int(match.group("index"))] = _CardFile(
                run_name=Path(run_args[binding.argIndex]).name, shown_name=binding.name
            )
    return files


def _append_linear_gui_sequence_options(
    run_args: list[str],
    invocation_args: list[str],
    *,
    linear_seqs: Sequence[Any],
    card_files: Mapping[int, _CardFile],
) -> None:
    labels = [
        str(seq.get("definition") or "") if isinstance(seq, Mapping) else ""
        for seq in linear_seqs
    ]
    if any(label.strip() for label in labels):
        for label in labels:
            _append_pair(run_args, invocation_args, "--record_label", label)

    subtitles = [
        str(seq.get("record_subtitle") or "") if isinstance(seq, Mapping) else ""
        for seq in linear_seqs
    ]
    if any(subtitle.strip() for subtitle in subtitles):
        for subtitle in subtitles:
            _append_pair(run_args, invocation_args, "--record_subtitle", subtitle)

    record_selectors: list[str] = []
    reverse_flags: list[bool] = []
    # Each region names the first record the CLI reads of its card's file (OV-3634): a bare ``#k`` counts every
    # record the CLI reads, and a card drawn whole reads every record of its file.
    region_specs: list[tuple[str, str]] = []

    def card_region(index: int, span: str) -> tuple[str, str]:
        card_file = card_files.get(index)
        if card_file is None:
            raise ValidationError(f"Linear sequence #{index + 1} has a region but no record file.")
        return f"{card_file.run_name}:#1:{span}", f"{card_file.shown_name}:#1:{span}"

    for index, seq in enumerate(linear_seqs):
        if not isinstance(seq, Mapping):
            record_selectors.append("")
            reverse_flags.append(False)
            continue
        record_selector = str(seq.get("region_record_id") or "").strip()
        record_selectors.append(record_selector)
        card_windows = _mapping_list(seq.get("windows"))
        if card_windows:
            # A card whose windows its crop fields cannot hold: several (D-01), or one across the origin (D-34). One
            # region per window, in order. The CLI places them in that order; the Record Layout columns are not
            # written, as for every other record.
            region_specs.extend(card_region(index, _card_window_span(index, window)) for window in card_windows)
            reverse_flags.append(False)
            continue
        start = seq.get("region_start")
        end = seq.get("region_end")
        has_start = start not in (None, "")
        has_end = end not in (None, "")
        if has_start != has_end:
            raise ValidationError(
                f"Linear sequence #{index + 1} has an incomplete region start/end."
            )
        wants_reverse = bool(seq.get("region_reverse"))
        if has_start and has_end:
            try:
                start_int = int(cast(Any, start))
                end_int = int(cast(Any, end))
            except (TypeError, ValueError) as exc:
                raise ValidationError(
                    f"Linear sequence #{index + 1} has invalid region coordinates."
                ) from exc
            if start_int < 1 or end_int < 1:
                raise ValidationError(
                    f"Linear sequence #{index + 1} region coordinates must be >= 1."
                )
            # The card's crop fields: an End before the Start reverses the span, as a region string without :origin.
            suffix = ":rc" if wants_reverse else ""
            region_specs.append(card_region(index, f"{start_int}-{end_int}{suffix}"))
            reverse_flags.append(False)
        else:
            reverse_flags.append(wants_reverse)

    if any(selector for selector in record_selectors):
        for selector in record_selectors:
            _append_pair(run_args, invocation_args, "--record_id", selector)
    if any(reverse_flags):
        for flag in reverse_flags:
            _append_pair(run_args, invocation_args, "--reverse_complement", "1" if flag else "0")
    for run_spec, shown_spec in region_specs:
        run_args.extend(["--region", run_spec])
        invocation_args.extend(["--region", shown_spec])


def _card_window_span(index: int, window: Mapping[str, object]) -> str:
    """The span of a ``--region`` of one window of Linear card ``index`` (0-based)."""

    start = _nonnegative_integer(window.get("start"))
    end = _nonnegative_integer(window.get("end"))
    if not start or not end:
        raise ValidationError(f"Linear sequence #{index + 1} has a window without a valid start and end.")
    return region_text(start, end, js_truthy(window.get('reverseComplement')))


def _gui_linear_losat_program(config: Mapping[str, Any], adv: Mapping[str, Any]) -> str:
    blast_source = str(config.get("blastSource") or adv.get("blastSource") or "").strip().lower()
    losat_program = str(config.get("losatProgram") or adv.get("losatProgram") or "").strip().lower()
    if blast_source != "losat":
        return ""
    return losat_program


def _append_linear_gui_blastp_args(
    run_args: list[str],
    invocation_args: list[str],
    *,
    session: Mapping[str, Any],
    config: Mapping[str, Any],
    adv: Mapping[str, Any],
) -> None:
    losat_cfg = config.get("losat")
    if not isinstance(losat_cfg, Mapping):
        return
    blastp_cfg = losat_cfg.get("blastp")
    if not isinstance(blastp_cfg, Mapping):
        return
    mode = str(blastp_cfg.get("mode") or "none").strip().lower()
    if mode not in {"pairwise", "orthogroup", "collinear"}:
        return
    _append_pair(run_args, invocation_args, "--losat", "losatp")
    _append_pair(
        run_args,
        invocation_args,
        "--losatp_mode",
        _losatp_mode_rewrites()[mode][-1],
    )
    threads_per_job = str(losat_cfg.get("threadsPerJob") or "auto").strip().lower()
    if threads_per_job != "auto":
        try:
            parsed_threads = int(threads_per_job)
        except ValueError:
            parsed_threads = 0
        if parsed_threads >= 1:
            _append_pair(run_args, invocation_args, "--losat_threads", str(parsed_threads))

    max_hits = blastp_cfg.get("maxHits")
    if max_hits not in (None, "", False):
        _append_pair(run_args, invocation_args, "--losatp_max_hits", str(max_hits))
        if mode == "pairwise":
            _append_pair(run_args, invocation_args, "--losatp_max_target_seqs", str(max_hits))
    candidate_limit = blastp_cfg.get("candidateLimit")
    if mode != "pairwise" and candidate_limit not in (None, "", False):
        _append_pair(run_args, invocation_args, "--losatp_max_target_seqs", str(candidate_limit))

    for key, option in (
        ("min_bitscore", "--bitscore"),
        ("evalue", "--evalue"),
        ("identity", "--identity"),
        ("alignment_length", "--alignment_length"),
    ):
        value = adv.get(key)
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))

    if mode == "orthogroup":
        orthogroup_state = session.get("orthogroupState")
        selected_target = (
            str(orthogroup_state.get("selectedOrthogroupAlignmentFeature") or "").strip()
            if isinstance(orthogroup_state, Mapping)
            else ""
        )
        if selected_target:
            _append_pair(run_args, invocation_args, "--similarity_alignment_feature", selected_target)

    if mode != "collinear":
        return
    for key, option in (
        ("collinearMinAnchors", "--collinear_min_anchors"),
        ("collinearMaxUnitGap", "--collinear_max_unit_gap"),
        ("collinearMaxDiagonalDrift", "--collinear_max_diagonal_drift"),
        ("collinearMaxConflictsInMergeGap", "--collinear_max_conflicts_in_merge_gap"),
        ("collinearUnitMode", "--collinear_unit_mode"),
        ("collinearSearchScope", "--collinear_search_scope"),
        ("collinearColorMode", "--collinear_color_mode"),
        ("collinearMaxParalogLinksPerOrthogroup", "--collinear_max_paralog_links_per_orthogroup"),
    ):
        value = blastp_cfg.get(key)
        if (
            key == "collinearMaxUnitGap"
            and value in (None, "")
            and int(session.get("version", 0)) < CURRENT_AUTHORITY_SESSION_MIN_VERSION
        ):
            value = blastp_cfg.get("collinearMaxGeneGap")
        if value not in (None, "", False):
            _append_pair(run_args, invocation_args, option, str(value))
    # Web Load reads a draft without the inference flag as on for a Collinear
    # search (services/config-field-readers.js), and the CLI default is off.
    inference = blastp_cfg.get("collinearInferOrthogroups")
    _append_pair(
        run_args,
        invocation_args,
        "--collinear_infer_orthogroups",
        "off" if inference is False else "on",
    )
    anchor_mode = blastp_cfg.get("collinearAnchorMode")
    if isinstance(anchor_mode, str) and anchor_mode.strip():
        anchor_mode = retired_collinearity_anchor_mode_replacement(anchor_mode) or anchor_mode
        _append_pair(run_args, invocation_args, "--collinear_anchor_mode", str(anchor_mode))
    merge_orientation = blastp_cfg.get("collinearMergeOrientation")
    if isinstance(merge_orientation, str) and merge_orientation.strip():
        # The parameters read the old spelling "either" as "strand".
        merge_orientation = read_spelling(
            merge_orientation, COLLINEARITY_MERGE_ORIENTATIONS, COLLINEARITY_MERGE_ORIENTATION_SPELLINGS
        ) or merge_orientation
        _append_pair(run_args, invocation_args, "--collinear_merge_orientation", merge_orientation)


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
            arg_index = int(cast(Any, item.get("argIndex")))
        except (TypeError, ValueError) as exc:
            raise ValidationError("cliInvocation.fileBindings argIndex must be an integer.") from exc
        slot = str(item.get("slot") or "").strip()
        if not slot:
            raise ValidationError("cliInvocation.fileBindings slot is required.")
        name = safe_embedded_filename(item.get("name"), fallback="file")
        bindings.append(SessionFileBinding(argIndex=arg_index, slot=slot, name=name))
    return bindings


def migrate_legacy_repeat_feature_shape_args(
    args: Sequence[str],
    *,
    session_version: int,
) -> list[str]:
    """Preserve the old repeat rectangle for non-canonical v27-30 replay."""

    migrated = [str(arg) for arg in args]
    if int(session_version) > 30:
        return migrated
    features_raw = _option_value(migrated, "-k", "--features")
    effective_features = (
        {item.strip() for item in features_raw.split(",") if item.strip()}
        if features_raw is not None
        else set(LEGACY_FEATURE_TYPES)
    )
    recorded_shapes = _feature_shapes_from_cli_args(migrated)
    legacy_shapes = [
        f"{feature_type}={shape}"
        for feature_type, shape in LEGACY_FEATURE_SHAPES.items()
        if feature_type in effective_features and feature_type not in recorded_shapes
    ]
    for legacy_shape in legacy_shapes:
        insertion_index = next(
            (
                index
                for index, token in enumerate(migrated)
                if token in {"-f", "--format"} or token.startswith("--format=")
            ),
            len(migrated),
        )
        migrated[insertion_index:insertion_index] = ["--feature_shape", legacy_shape]
    return migrated


def _feature_shapes_from_cli_args(args: Sequence[str]) -> dict[str, str]:
    shapes: dict[str, str] = {}
    for assignment in _option_all_values(args, "--feature_shape", "--feature-shape"):
        feature_type, separator, shape = str(assignment).partition("=")
        feature_type = feature_type.strip()
        shape = shape.strip().lower()
        if separator and feature_type and shape in {
            "arrow",
            "rectangle",
            "underlay",
        }:
            shapes[feature_type] = shape
    return shapes


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
    "session_to_cli_args",
    "RetiredCliOption",
    "RETIRED_CLI_OPTIONS",
    "canonicalize_cli_invocation",
    "migrate_legacy_repeat_feature_shape_args",
]
