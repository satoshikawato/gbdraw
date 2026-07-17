#!/usr/bin/env python
# coding: utf-8


import argparse
import copy
import hashlib
import json
import logging
import math
import re
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Mapping, Optional, Sequence
from pandas import DataFrame  # type: ignore[reportMissingImports]
from .io.colors import load_default_colors, read_color_table
from .io.cli_tables import read_records_table
from .io.genome import load_gbks, load_gff_fasta
from .io.regions import apply_region_specs, parse_region_specs
from .config.toml import load_config_toml
from .render.export import parse_formats, save_figure
from .render.formats import INTERACTIVE_SVG_FORMAT
from .render.interactive_svg import InteractiveSvgContext
from .render.interactive_context import build_interactive_svg_context
from .api.diagram import assemble_linear_diagram_from_records  # type: ignore[reportMissingImports]
from .api.options import AnnotationOptions, ColorOptions, DiagramOptions, OutputOptions, TrackOptions
from .annotations import read_annotation_table
from .api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from .definition_line_styles import (
    parse_definition_line_style_assignment,
    parse_definition_line_style_overrides,
)
from .analysis.collinearity import (
    LosslessCollinearityParameters,
    build_orthogroup_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    normalize_collinearity_search_scope,
)
from .analysis.protein_colinearity import (
    LosatpCacheManager,
    ORTHOGROUP_INFERENCE_VERSION,
    PROTEIN_BLASTP_MODES,
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    extract_web_stable_cds_proteins,
)
from .config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from .config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from .labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)  # type: ignore[reportMissingImports]
from .features.shapes import parse_feature_shape_overrides
from .features.visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
    resolve_candidate_feature_types,
)
from .exceptions import ValidationError
from .tracks import (
    linear_track_slots_from_order,
    normalize_linear_track_slots_with_axis,
    parse_linear_track_slots,
)


from .cli_utils.common import (
    _add_block_stroke_args,
    _add_depth_axis_args,
    _add_depth_track_label_color_args,
    _add_depth_track_tick_args,
    _add_feature_shape_arg,
    _add_format_arg,
    _add_gc_content_axis_args,
    _add_legend_size_args,
    _add_window_step_args,
    add_feature_args,
    add_input_args,
    add_label_args,
    setup_logging,
    validate_input_args,
    validate_label_args,
    handle_output_formats,
    calculate_window_step,
    load_records_table_records as _load_records_table_records,
    record_major_depth_track_files_from_cli as _record_major_depth_track_files_from_cli,
)
from .cli_utils.session import (
    DiagramRunResult,
    add_session_args,
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
    make_rendered_svg,
    parse_session_pre_args,
    render_canonical_session_if_present,
    save_session_sidecar_if_requested,
)
from .session_io import load_session, session_to_cli_args


def _add_argument_with_hidden_aliases(
    parser: argparse.ArgumentParser,
    *option_strings: str,
    hidden_aliases: Sequence[str] = (),
    **kwargs: Any,
) -> None:
    """Add public underscore options and hidden legacy hyphen aliases."""

    parser.add_argument(*option_strings, **kwargs)
    if not hidden_aliases:
        return
    alias_kwargs = dict(kwargs)
    alias_kwargs["help"] = argparse.SUPPRESS
    alias_kwargs["default"] = argparse.SUPPRESS
    parser.add_argument(*hidden_aliases, **alias_kwargs)


def _web_json_dumps(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _web_hash_text(text: str) -> str:
    return hashlib.sha256(str(text).encode("utf-8")).hexdigest()


def _web_safe_filename(name: object, fallback: str = "losat") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name or ""))
    cleaned = cleaned.strip("_")
    return cleaned or fallback


def _web_normalize_label(label: object, fallback: str) -> str:
    base = str(label or "").strip() or str(fallback or "")
    dotted = re.sub(r"\.+", ".", re.sub(r"[\s/]+", ".", base)).strip(".")
    return _web_safe_filename(dotted, _web_safe_filename(fallback, "losat"))


def _record_cache_label(record: object, fallback: str) -> str:
    annotations = getattr(record, "annotations", {}) or {}
    label = str(annotations.get("gbdraw_record_label") or "").strip()
    if label:
        return label
    record_id = str(getattr(record, "id", "") or "").strip()
    return record_id or fallback


def _linear_losat_cache_filenames(records: Sequence[object]) -> tuple[str, ...]:
    filenames: list[str] = []
    for index in range(max(0, len(records) - 1)):
        left = _web_normalize_label(
            _record_cache_label(records[index], f"seq_{index + 1}"),
            f"seq_{index + 1}",
        )
        right = _web_normalize_label(
            _record_cache_label(records[index + 1], f"seq_{index + 2}"),
            f"seq_{index + 2}",
        )
        filenames.append(f"{left}.{right}.losatp.tsv")
    return tuple(filenames)


def _linear_session_record_metadata(
    records: Sequence[object],
    input_paths: Sequence[str],
) -> tuple[Mapping[str, object], ...]:
    source_indices: list[int] = []
    for record_index, record in enumerate(records):
        source_indices.append(
            _linear_record_source_index(
                record,
                input_paths,
                fallback_index=record_index if len(records) == len(input_paths) else None,
            )
        )
    source_counts: dict[int, int] = {}
    for source_index in source_indices:
        if source_index >= 0:
            source_counts[source_index] = source_counts.get(source_index, 0) + 1
    source_seen: dict[int, int] = {}
    metadata: list[Mapping[str, object]] = []
    for loaded_index, (record, source_index) in enumerate(zip(records, source_indices)):
        source_loaded_index = source_seen.get(source_index, 0)
        if source_index >= 0:
            source_seen[source_index] = source_loaded_index + 1
        annotations = getattr(record, "annotations", {}) or {}
        source_file = str(
            annotations.get("gbdraw_source_file")
            or (input_paths[source_index] if 0 <= source_index < len(input_paths) else "")
        )
        source_basename = str(
            annotations.get("gbdraw_source_basename")
            or (Path(source_file).name if source_file else "")
        )
        metadata.append(
            {
                "loaded_index": loaded_index,
                "source_index": source_index,
                "source_loaded_index": source_loaded_index,
                "source_loaded_count": source_counts.get(source_index, 0),
                "record_id": str(getattr(record, "id", "") or ""),
                "source_file": source_file,
                "source_basename": source_basename,
            }
        )
    return tuple(metadata)


def _linear_record_source_index(
    record: object,
    input_paths: Sequence[str],
    *,
    fallback_index: int | None,
) -> int:
    annotations = getattr(record, "annotations", {}) or {}
    source_file = str(annotations.get("gbdraw_source_file") or "")
    source_basename = str(annotations.get("gbdraw_source_basename") or "")
    if source_file:
        for index, path in enumerate(input_paths):
            if source_file == str(path):
                return index
        try:
            resolved_source = str(Path(source_file).resolve())
        except OSError:
            resolved_source = ""
        if resolved_source:
            for index, path in enumerate(input_paths):
                try:
                    if resolved_source == str(Path(path).resolve()):
                        return index
                except OSError:
                    continue
    candidates = {source_basename, Path(source_file).name if source_file else ""}
    basenames: dict[str, list[int]] = {}
    for index, path in enumerate(input_paths):
        basenames.setdefault(Path(str(path)).name, []).append(index)
    for candidate in candidates:
        if candidate and len(basenames.get(candidate, [])) == 1:
            return basenames[candidate][0]
    if fallback_index is not None and 0 <= fallback_index < len(input_paths):
        return fallback_index
    return -1


def _source_session_losat_entries(source_session: Mapping[str, Any] | None) -> tuple[Mapping[str, object], ...]:
    if not isinstance(source_session, Mapping):
        return ()
    losat_cache = source_session.get("losatCache")
    if not isinstance(losat_cache, Mapping):
        return ()
    entries = losat_cache.get("entries")
    if not isinstance(entries, list):
        return ()
    return tuple(entry for entry in entries if isinstance(entry, Mapping))


def _source_session_losat_config(source_session: Mapping[str, Any] | None) -> Mapping[str, Any]:
    if not isinstance(source_session, Mapping):
        return {}
    config = source_session.get("config")
    if not isinstance(config, Mapping):
        return {}
    losat = config.get("losat")
    return losat if isinstance(losat, Mapping) else {}


def _source_session_threads_per_job(source_session: Mapping[str, Any] | None) -> str | int:
    losat = _source_session_losat_config(source_session)
    raw = str(losat.get("threadsPerJob") or "auto").strip().lower()
    if raw == "auto":
        return "auto"
    try:
        parsed = int(raw)
    except ValueError:
        return "auto"
    return parsed if parsed >= 1 else "auto"


def _source_session_runtime_compatibility(source_session: Mapping[str, Any] | None) -> str:
    losat = _source_session_losat_config(source_session)
    execution_mode = str(losat.get("executionMode") or "auto").strip().lower()
    return "serial-v1" if execution_mode == "serial" else "threaded-compatible-v1"


def _fingerprint_from_session_entry(entry: object) -> dict[str, object] | None:
    if not isinstance(entry, Mapping):
        return None
    return {
        "name": str(entry.get("name") or ""),
        "size": int(entry.get("size") or 0),
        "lastModified": int(entry.get("lastModified") or 0),
    }


def _fingerprint_from_path(path: object) -> dict[str, object] | None:
    if path in (None, "", False):
        return None
    file_path = Path(str(path))
    try:
        stat = file_path.stat()
    except OSError:
        return {"name": file_path.name, "size": 0, "lastModified": 0}
    return {
        "name": file_path.name,
        "size": int(stat.st_size),
        "lastModified": int(stat.st_mtime * 1000),
    }


def _linear_session_seq_entry(
    source_session: Mapping[str, Any] | None,
    index: int,
) -> Mapping[str, Any] | None:
    if not isinstance(source_session, Mapping):
        return None
    files = source_session.get("files")
    if not isinstance(files, Mapping):
        return None
    linear_seqs = files.get("linearSeqs")
    if not isinstance(linear_seqs, list) or index >= len(linear_seqs):
        return None
    entry = linear_seqs[index]
    return entry if isinstance(entry, Mapping) else None


def _canonical_region_spec_from_session_seq(seq: Mapping[str, Any] | None) -> str | None:
    if not isinstance(seq, Mapping):
        return None
    start = seq.get("region_start")
    end = seq.get("region_end")
    if start in (None, "") or end in (None, ""):
        return None
    try:
        start_int = int(start)
        end_int = int(end)
    except (TypeError, ValueError):
        return None
    return f"{min(start_int, end_int)}-{max(start_int, end_int)}"


def _record_instance_keys_for_web_losat(
    args: argparse.Namespace,
    *,
    source_session: Mapping[str, Any] | None,
    record_count: int,
) -> tuple[str, ...]:
    input_format = "gb" if args.gbk else "gff"
    primary_paths = list(args.gbk or args.gff or [])
    paired_fasta_paths = list(args.fasta or [])
    record_selectors = list(args.record_id or [])
    occurrences: dict[str, int] = {}
    keys: list[str] = []
    for index in range(record_count):
        session_seq = _linear_session_seq_entry(source_session, index)
        if session_seq is not None:
            primary_entry = session_seq.get("gb") if input_format == "gb" else session_seq.get("gff")
            paired_entry = session_seq.get("fasta") if input_format == "gff" else None
            primary_fingerprint = _fingerprint_from_session_entry(primary_entry)
            paired_fingerprint = _fingerprint_from_session_entry(paired_entry)
            record_selector = str(session_seq.get("region_record_id") or "").strip()
            canonical_region = _canonical_region_spec_from_session_seq(session_seq)
        else:
            primary_fingerprint = _fingerprint_from_path(
                primary_paths[index] if index < len(primary_paths) else None
            )
            paired_fingerprint = _fingerprint_from_path(
                paired_fasta_paths[index] if input_format == "gff" and index < len(paired_fasta_paths) else None
            )
            record_selector = str(record_selectors[index] if index < len(record_selectors) else "").strip()
            canonical_region = None
        descriptor = {
            "inputFormat": input_format,
            "primaryFile": primary_fingerprint,
            "pairedFastaFile": paired_fingerprint,
            "recordSelector": record_selector,
            "canonicalRegionSpec": canonical_region,
        }
        base = _web_hash_text(_web_json_dumps(descriptor))[:16]
        occurrence = occurrences.get(base, 0) + 1
        occurrences[base] = occurrence
        keys.append(f"r_{base}_o{occurrence}" if occurrence > 1 else f"r_{base}")
    return tuple(keys)


def _parse_optional_positive_int(value: str) -> int | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"", "none", "auto", "null"}:
        return None
    try:
        parsed = int(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'")
    return parsed


# Setup for the logging system
logger = logging.getLogger()
setup_logging()


def _build_interactive_svg_context(
    records,
    selected_features,
    orthogroups=None,
    *,
    feature_visibility_rules=None,
    color_table=None,
    default_colors=None,
    annotations=None,
) -> InteractiveSvgContext:
    try:
        return build_interactive_svg_context(
            records,
            selected_features_set=selected_features,
            feature_visibility_rules=feature_visibility_rules,
            color_table=color_table,
            default_colors=default_colors,
            orthogroups=orthogroups,
            linear_rendered_feature_ids=True,
            annotations=annotations,
            mode="linear",
        )
    except Exception as exc:
        logger.warning(
            "WARNING: Rich interactive feature metadata could not be generated; "
            "falling back to rendered SVG feature metadata. %s",
            exc,
        )
        return InteractiveSvgContext()

def _parse_linear_label_placement(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized == "on_feature":
        return "above_feature"
    if normalized in {"auto", "above_feature"}:
        return normalized
    raise argparse.ArgumentTypeError(
        "label placement must be 'auto' or 'above_feature' (legacy alias: 'on_feature')"
    )


def _parse_linear_track_layout(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized in {"above", "spreadout"}:
        return "above"
    if normalized in {"below", "tuckin"}:
        return "below"
    if normalized == "middle":
        return "middle"
    raise argparse.ArgumentTypeError(
        "track layout must be 'above', 'middle', or 'below' "
        "(legacy aliases: 'spreadout' -> 'above', 'tuckin' -> 'below')"
    )


def _parse_linear_track_axis_gap(value: str) -> float | None:
    normalized = str(value).strip().lower()
    if normalized in {"", "auto"}:
        return None
    try:
        axis_gap = float(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "track axis gap must be a non-negative number or 'auto'"
        ) from exc
    if axis_gap < 0:
        raise argparse.ArgumentTypeError("track axis gap must be >= 0")
    return axis_gap


def _parse_pairwise_match_style(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized not in {"ribbon", "curve"}:
        raise argparse.ArgumentTypeError("pairwise_match_style must be one of: ribbon, curve")
    return normalized


def _parse_collinear_color_mode(value: str) -> str:
    normalized = str(value).strip().lower().replace("-", "_")
    if normalized == "identity":
        normalized = "average_identity"
    if normalized not in {"average_identity", "orientation", "orientation_identity"}:
        raise argparse.ArgumentTypeError(
            "collinear_color_mode must be one of: average_identity, orientation, orientation_identity"
        )
    return normalized


def _parse_collinear_search_scope(value: str) -> str:
    try:
        return normalize_collinearity_search_scope(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc




def _parse_definition_line_style_arg(value: str) -> str:
    try:
        parse_definition_line_style_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value






def _get_args(args) -> argparse.Namespace:
    """
    Parses command-line arguments for generating linear genome diagrams.

    This internal function defines and parses command-line arguments using argparse.
    It sets up the necessary parameters required for the linear genome diagram creation,
    including input files, output configurations, and various visualization options.

    Args:
        args (list of str): The command-line arguments passed to the script.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    The function supports a variety of arguments for input files (GenBank and BLAST),
    output Configurator (file format and naming), and visualization preferences (color
    tables, window size, and feature selection).
    """
    parser = argparse.ArgumentParser(
        description='Generate  plot in PNG/PDF/SVG/PS/EPS.')
    add_input_args(parser)
    parser.add_argument(
        "--records_table",
        metavar="TSV",
        help="TSV manifest for row-based input records and per-record options.",
        type=str)
    parser.add_argument(
        '-b',
        '--blast',
        help="input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)",
        type=str,
        nargs='*')
    _add_argument_with_hidden_aliases(
        parser,
        '--losatp_bin',
        hidden_aliases=('--losatp-bin',),
        dest='losatp_bin',
        help='Native LOSAT executable for --protein_blastp_mode pairwise/orthogroup/collinear (default: losat).',
        type=str,
        default='losat')
    _add_argument_with_hidden_aliases(
        parser,
        '--ncbi_blastp_bin',
        hidden_aliases=('--ncbi-blastp-bin',),
        dest='ncbi_blastp_bin',
        help='NCBI BLAST+ blastp executable for --protein_blastp_mode pairwise/orthogroup/collinear (default: use automatic runtime resolution).',
        type=str,
        default=None)
    _add_argument_with_hidden_aliases(
        parser,
        '--losatp_threads',
        hidden_aliases=('--losatp-threads',),
        dest='losatp_threads',
        help='Threads passed to the selected protein blastp runtime for --protein_blastp_mode pairwise/orthogroup/collinear (default: runtime default).',
        type=int,
        default=None)
    _add_argument_with_hidden_aliases(
        parser,
        '--protein_blastp_mode',
        hidden_aliases=('--protein-blastp-mode',),
        dest='protein_blastp_mode',
        help='Protein blastp comparison mode: none, pairwise adjacent ribbons, all-record similarity groups (orthogroup), or collinear blocks (default: none).',
        choices=PROTEIN_BLASTP_MODES,
        default='none')
    _add_argument_with_hidden_aliases(
        parser,
        '--protein_blastp_max_hits',
        hidden_aliases=('--protein-blastp-max-hits',),
        dest='protein_blastp_max_hits',
        help='Maximum distinct subject protein hits per query protein for pairwise protein blastp display links (default: 5).',
        type=int,
        default=5)
    _add_argument_with_hidden_aliases(
        parser,
        '--protein_blastp_candidate_limit',
        hidden_aliases=('--protein-blastp-candidate-limit',),
        dest='protein_blastp_candidate_limit',
        help="Optional protein blastp candidate cap per query; use 'none' for no cap (default: none).",
        type=_parse_optional_positive_int,
        default=None)
    _add_argument_with_hidden_aliases(
        parser,
        '--align_orthogroup_feature',
        hidden_aliases=('--align-orthogroup-feature',),
        dest='align_orthogroup_feature',
        help='Align linear records by the gbdraw similarity group containing this feature SVG hash or protein ID.',
        type=str,
        default="")
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_unit_mode',
        hidden_aliases=('--collinear-unit-mode',),
        dest='collinear_unit_mode',
        help=argparse.SUPPRESS,
        choices=["auto", "cds", "locus"],
        default='auto')
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_search_scope',
        hidden_aliases=('--collinear-search-scope',),
        dest='collinear_search_scope',
        help='Collinear protein blastp evidence search scope: adjacent record pairs or all record pairs (default: adjacent).',
        type=_parse_collinear_search_scope,
        choices=["adjacent", "all"],
        default='adjacent')
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_min_anchors',
        hidden_aliases=('--collinear-min-anchors',),
        dest='collinear_min_anchors',
        help='Minimum anchors/genes required for a rendered Collinear block; 1 allows singleton links (default: 1).',
        type=int,
        default=1)
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_max_unit_gap',
        '--collinear_max_gene_gap',
        hidden_aliases=('--collinear-max-unit-gap', '--collinear-max-gene-gap'),
        dest='collinear_max_unit_gap',
        help='Maximum unit gap between neighboring collinear anchors (default: 0).',
        type=int,
        default=0)
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_max_diagonal_drift',
        hidden_aliases=('--collinear-max-diagonal-drift',),
        dest='collinear_max_diagonal_drift',
        help='Maximum order-space diagonal drift allowed within or between collinear runs (default: 0).',
        type=int,
        default=0)
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_max_conflicts_in_merge_gap',
        hidden_aliases=('--collinear-max-conflicts-in-merge-gap',),
        dest='collinear_max_conflicts_in_merge_gap',
        help=argparse.SUPPRESS,
        type=int,
        default=1)
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_max_paralog_links_per_orthogroup',
        hidden_aliases=('--collinear-max-paralog-links-per-orthogroup',),
        dest='collinear_max_paralog_links_per_orthogroup',
        help=argparse.SUPPRESS,
        type=int,
        default=2)
    _add_argument_with_hidden_aliases(
        parser,
        '--collinear_color_mode',
        hidden_aliases=('--collinear-color-mode',),
        dest='collinear_color_mode',
        help='Collinear ribbon color mode: average_identity, orientation, or orientation_identity (default: orientation).',
        type=_parse_collinear_color_mode,
        choices=["average_identity", "orientation", "orientation_identity"],
        default='orientation')
    parser.add_argument(
        '-t',
        '--table',
        help='color table (optional)',
        type=str,
        default="")
    parser.add_argument(
        "-p", "--palette",
        metavar="PALETTE",
        default="default",
        help="Palette name (default: default)",
        type=str
    )
    parser.add_argument(
        '-d',
        '--default_colors',
        help='TSV file that overrides the color palette (optional)',
        type=str,
        default="")
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: out)',
        type=str,
        default="out")
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide skew (default: GC). ',
        type=str,
        default="GC")
    _add_window_step_args(parser)
    parser.add_argument(
        '--separate_strands',
        help='separate forward and reverse strands (default: False). Features of undefined strands are shown on the forward strand.',
        action='store_true')
    parser.add_argument(
        '--show_gc',
        help='plot GC content below genome (default: False). ',
        action='store_true')
    _add_gc_content_axis_args(parser)
    parser.add_argument(
        '--show_skew',
        help='plot GC skew below genome (default: False). ',
        action='store_true')
    parser.add_argument(
        '--depth',
        help='Depth TSV file(s) in samtools depth format. Provide one for all records or one per input record.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--depth_track',
        metavar='DEPTH',
        help='Repeatable logical depth track. Provide one file for all records or one file per input record.',
        type=str,
        nargs='+',
        action='append')
    _add_depth_track_label_color_args(parser)
    parser.add_argument(
        '--depth_track_height',
        metavar='PX',
        help='Linear depth track height(s) in px. Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    _add_depth_track_tick_args(parser)
    parser.add_argument(
        '--depth_color',
        help='Depth track fill color (optional; default: #4A90E2).',
        type=str)
    parser.add_argument(
        '--depth_height',
        help='Depth track height for linear mode (in px; must be > 0).',
        type=float)
    _add_depth_axis_args(parser)
    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
        action='store_true')
    _add_argument_with_hidden_aliases(
        parser,
        '--keep_definition_left_aligned',
        hidden_aliases=('--keep-definition-left-aligned',),
        dest='keep_definition_left_aligned',
        help='Keep the linear record-label block in the left column when records are center-aligned or aligned by a gbdraw similarity group (default: False).',
        action='store_true')
    parser.add_argument(
        '--evalue',
        help='evalue threshold (default=1e-2)',
        type=float,
        default="1e-2")
    parser.add_argument(
        '--bitscore',
        help='bitscore threshold (default=50)',
        type=float,
        default="50")
    parser.add_argument(
        '--identity',
        help='identity threshold (default=0)',
        type=float,
        default="0")
    parser.add_argument(
        '--alignment_length',
        help='minimum BLAST alignment length threshold (default=0)',
        type=int,
        default=0)
    _add_argument_with_hidden_aliases(
        parser,
        '--pairwise_match_style',
        hidden_aliases=('--pairwise-match-style',),
        dest='pairwise_match_style',
        help=(
            'Pairwise comparison link style: ribbon keeps straight filled ribbons; '
            'curve draws curved filled ribbons that preserve alignment spans.'
        ),
        type=_parse_pairwise_match_style,
        choices=["ribbon", "curve"],
        default="ribbon")
    add_feature_args(parser)
    _add_feature_shape_arg(parser)
    _add_block_stroke_args(parser)
    parser.add_argument(
        '--axis_stroke_color',
        help='Axis stroke color (str; default: auto: "lightgray", or "dimgray" with --ruler_on_axis)',
        type=str,
        default=None)
    parser.add_argument(
        '--axis_stroke_width',
        help='Axis stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 2 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (optional; str; default: "lightgray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--definition_font_size',
        help='Definition font size (optional; float; default: 24 pt for genomes <= 50 kb, 10 pt for genomes >= 50 kb)',
        type=float)
    _add_argument_with_hidden_aliases(
        parser,
        '--definition_line_style',
        hidden_aliases=('--definition-line-style',),
        dest='definition_line_style',
        help=(
            'Definition line style override (repeatable): '
            'LINE:weight=bold,color=#000000,size=12. '
            'Lines: name/species/record_label, subtitle/record_subtitle, replicon, accession, length/coordinates.'
        ),
        type=_parse_definition_line_style_arg,
        action='append',
        default=[],
        metavar='LINE:KEY=VALUE',
    )
    parser.add_argument(
        '--plot_title',
        help='Shared plot title text (optional).',
        type=str,
        default="")
    parser.add_argument(
        '--plot_title_position',
        help='Shared plot title position ("center", "top", "bottom"; default: "bottom").',
        type=str,
        choices=["center", "top", "bottom"],
        default="bottom")
    parser.add_argument(
        '--plot_title_font_size',
        help='Shared plot title font size (optional; float; default: 32).',
        type=float)
    parser.add_argument(
        '--record_label',
        help='Optional top record-label line (for example organism/strain; repeatable; order matches input records)',
        type=str,
        action='append',
        default=[])
    _add_argument_with_hidden_aliases(
        parser,
        '--record_subtitle',
        hidden_aliases=('--record-subtitle',),
        dest='record_subtitle',
        help='Optional second record-label line (repeatable; order matches input records)',
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--show_replicon',
        help='Show inferred replicon labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_accession',
        help='Hide accession labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_length',
        help='Hide length/coordinate labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--label_font_size',
        help='Label font size (optional; default: 24 pt for genomes <= 50 kb, 5 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--label_placement',
        help='Linear label placement mode ("auto" or "above_feature"; default: "auto"). "above_feature" draws labels above features (or below negative-strand features when --separate_strands is used).',
        type=_parse_linear_label_placement,
        metavar="{auto,above_feature}",
    )
    parser.add_argument(
        '--label_rendering',
        help='Label rendering policy: "auto" embeds fitting labels and routes others externally; "embedded_only" drops external labels; "external_only" forces labels outside feature bodies. Cannot be combined with --label_placement above_feature except as "auto". Default: "auto".',
        choices=['auto', 'embedded_only', 'external_only'],
        default='auto',
        type=str)
    parser.add_argument(
        '--label_rotation',
        help='Linear label rotation in degrees (optional; float; default: 0). In above_feature mode, rotated labels start from the feature midpoint.',
        type=float,
    )
    parser.add_argument(
        '--linear_label_spacing',
        help='Linear label-to-label vertical spacing in px (optional; float; must be > 0).',
        type=float,
    )
    parser.add_argument(
        '--track_layout',
        help=(
            'Linear track layout mode ("above", "middle", or "below"; default: "middle"). '
            'Aliases: "spreadout" -> "above", "tuckin" -> "below".'
        ),
        type=_parse_linear_track_layout,
        metavar="{above,middle,below}",
        default="middle",
    )
    parser.add_argument(
        '--track_axis_gap',
        help=(
            "Gap between axis and nearest feature edge in pixels for above/below layouts. "
            "Use 'auto' to derive it from feature height."
        ),
        type=_parse_linear_track_axis_gap,
        default=None,
        metavar="AUTO|PX",
    )
    parser.add_argument(
        '--linear_track_order',
        help='Linear custom track shortcut order, for example features,depth,gc_content,gc_skew.',
        type=str,
        default="",
    )
    parser.add_argument(
        '--linear_track_slot',
        help='Linear custom track slot: <slot_id>:<renderer>@key=value,key=value. Repeat to add slots. The annotations renderer requires set_id from --annotation_table; overlay annotations also require anchor_slot and layer.',
        action='append',
        default=[],
        metavar='SLOT',
    )
    parser.add_argument(
        '--linear_track_axis_index',
        help='Axis boundary index for linear custom track slots.',
        type=int,
        default=None,
    )
    parser.add_argument(
        '--ruler_on_axis',
        help=(
            'Use each record axis as the ruler in linear mode. '
            'Effective only with --scale_style ruler and --track_layout above|below.'
        ),
        action='store_true',
    )
    _add_format_arg(parser)
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "right", "left", "top", "bottom", "none")',
        type=str,
        default="right")
    parser.add_argument(
            "--show_labels",
            help="Show labels: no argument or 'all' (all records), 'first' (first record only), 'orthogroup_top' (topmost record containing each gbdraw similarity group), 'none' (no labels). Default: 'none'",
            nargs='?',
            const="all",
            default="none",
            choices=["all", "first", "orthogroup_top", "none"],
            type=str
        )
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlaps (experimental; default: False). ',
        action='store_true')
    add_label_args(parser)
    parser.add_argument(
        '--feature_height',
        help='Feature vertical width (optional; float; default: 80 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb)',
        type=float),
    parser.add_argument(
        '--gc_height',
        help='GC content/skew vertical width (optional; float; default: 20 (pixels, 96 dpi))',
        type=float),
    parser.add_argument(
        '--comparison_height',
        help='Comparison block height (optional; float; optional; default: 60 (pixels, 96 dpi))',
        type=float)
    parser.add_argument(
        '--scale_style',
        help='Style for the length scale (default: "bar"; "bar", "ruler")',
        type=str,
        choices=["bar", "ruler"],
        default="bar")
    parser.add_argument(
        '--scale_stroke_color',
        help='Scale bar/ruler stroke color (optional; str; default: "black")',
        type=str)
    parser.add_argument(
        '--scale_stroke_width',
        help='Scale bar/ruler stroke width (optional; float; default: 3 (pt))',
        type=float)
    parser.add_argument(
        '--scale_font_size',
        help='Scale bar/ruler font size (optional; float; default: 24 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--ruler_label_font_size',
        help='Ruler label font size (optional; float). Overrides --scale_font_size when both are set.',
        type=float)
    parser.add_argument(
        '--ruler_label_color',
        help='Ruler label color (optional; str; default follows axis color when --ruler_on_axis is active, otherwise black).',
        type=str)
    parser.add_argument(
        '--scale_interval',
        help='Manual tick interval for "ruler" scale style (in bp). Overrides automatic calculation; optional',
        type=int)
    _add_legend_size_args(parser)
    parser.add_argument(
        '--normalize_length',
        help='Normalize record length (experimental; default: False). ',
        action='store_true')
    parser.add_argument(
        '--region',
        help=(
            'Crop a region (repeatable). Format: record_id:start-end[:rc], #index:start-end[:rc], '
            'or file:record_selector:start-end[:rc]. '
            'Coordinates are 1-based inclusive. For multiple records without selectors, provide one spec per record in input order (file order, then record order within each file).'
        ),
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--record_id',
        help=(
            'Select a record by ID or #index per input file (repeatable; order matches input files). '
            'Use an empty value to skip selection for a file.'
        ),
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--reverse_complement',
        help=(
            'Reverse complement record per input file (repeatable; order matches input files). '
            'Accepted values: 1/0, true/false, yes/no.'
        ),
        type=str,
        action='append',
        default=[])
    add_session_args(parser)

    args = parser.parse_args(args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if args.records_table:
        for option_name in (
            "record_label",
            "record_subtitle",
            "record_id",
            "reverse_complement",
            "region",
        ):
            if getattr(args, option_name):
                parser.error(
                    f"--records_table cannot be combined with --{option_name}; "
                    f"use the records table {option_name} column instead."
                )
    if args.protein_blastp_mode != "none" and args.blast:
        parser.error("--protein_blastp_mode cannot be used with -b/--blast")
    if args.protein_blastp_max_hits <= 0:
        parser.error("--protein_blastp_max_hits must be > 0")
    if args.losatp_threads is not None and args.losatp_threads <= 0:
        parser.error("--losatp_threads must be > 0")
    if args.align_orthogroup_feature and args.protein_blastp_mode != "orthogroup" and not args.blast:
        parser.error("--align_orthogroup_feature requires --protein_blastp_mode orthogroup")
    if args.depth and args.depth_track:
        parser.error("--depth cannot be combined with --depth_track")
    if args.show_depth and not (args.depth or args.depth_track):
        parser.error("--show_depth requires --depth or --depth_track")
    if args.depth_height is not None and args.depth_height <= 0:
        parser.error("--depth_height must be > 0")
    if args.depth_window is not None and args.depth_window <= 0:
        parser.error("--depth_window must be > 0")
    if args.depth_step is not None and args.depth_step <= 0:
        parser.error("--depth_step must be > 0")
    if args.depth_min is not None and args.depth_min < 0:
        parser.error("--depth_min must be >= 0")
    if args.depth_max is not None and args.depth_max < 0:
        parser.error("--depth_max must be >= 0")
    if args.depth_min is not None and args.depth_max is not None and args.depth_min > args.depth_max:
        parser.error("--depth_min must be <= --depth_max")
    if args.depth_tick_interval is not None and args.depth_tick_interval <= 0:
        parser.error("--depth_tick_interval must be > 0")
    if args.depth_large_tick_interval is not None and args.depth_large_tick_interval <= 0:
        parser.error("--depth_large_tick_interval must be > 0")
    if args.depth_small_tick_interval is not None and args.depth_small_tick_interval <= 0:
        parser.error("--depth_small_tick_interval must be > 0")
    if args.depth_tick_font_size is not None and args.depth_tick_font_size <= 0:
        parser.error("--depth_tick_font_size must be > 0")
    for option_name in (
        "depth_track_height",
        "depth_track_large_tick_interval",
        "depth_track_small_tick_interval",
        "depth_track_tick_font_size",
    ):
        for option_value in getattr(args, option_name) or []:
            option_text = str(option_value).strip().lower()
            if option_text in {"", "auto", "none", "null", "-"}:
                continue
            try:
                numeric_option_value = float(option_value)
            except (TypeError, ValueError):
                parser.error(f"--{option_name} values must be numbers or auto")
            if numeric_option_value <= 0:
                parser.error(f"--{option_name} values must be > 0")
    for option_name in ("gc_content_min_percent", "gc_content_max_percent"):
        option_value = getattr(args, option_name)
        if option_value is not None and not math.isfinite(float(option_value)):
            parser.error(f"--{option_name} must be a finite number")
    if (
        args.gc_content_min_percent is not None
        and args.gc_content_max_percent is not None
        and args.gc_content_min_percent > args.gc_content_max_percent
    ):
        parser.error("--gc_content_min_percent must be <= --gc_content_max_percent")
    if args.gc_content_tick_interval is not None and args.gc_content_tick_interval <= 0:
        parser.error("--gc_content_tick_interval must be > 0")
    if args.gc_content_large_tick_interval is not None and args.gc_content_large_tick_interval <= 0:
        parser.error("--gc_content_large_tick_interval must be > 0")
    if args.gc_content_small_tick_interval is not None and args.gc_content_small_tick_interval <= 0:
        parser.error("--gc_content_small_tick_interval must be > 0")
    if args.gc_content_tick_font_size is not None and args.gc_content_tick_font_size <= 0:
        parser.error("--gc_content_tick_font_size must be > 0")
    if args.linear_track_order and args.linear_track_slot:
        parser.error("--linear_track_order cannot be combined with --linear_track_slot")
    if args.linear_track_axis_index is not None and not (args.linear_track_order or args.linear_track_slot):
        parser.error("--linear_track_axis_index requires --linear_track_order or --linear_track_slot")
    try:
        linear_track_slot_specs = None
        if args.linear_track_order:
            linear_track_slot_specs = linear_track_slots_from_order(
                args.linear_track_order,
                show_depth=bool(args.show_depth or args.depth or args.depth_track),
                depth_track_count=max(1, len(args.depth_track or [])),
                show_gc=bool(args.show_gc),
                show_skew=bool(args.show_skew),
                dinucleotide=str(args.nt or "GC").upper(),
                track_layout=args.track_layout,
            )
        elif args.linear_track_slot:
            linear_track_slot_specs = parse_linear_track_slots(args.linear_track_slot)
        if linear_track_slot_specs is not None and args.linear_track_axis_index is not None:
            normalize_linear_track_slots_with_axis(
                linear_track_slot_specs,
                args.linear_track_axis_index,
            )
        args.linear_track_slot_specs = linear_track_slot_specs
    except Exception as exc:
        parser.error(str(exc))
    if args.label_placement == "above_feature" and args.label_rendering != "auto":
        parser.error("--label_rendering embedded_only|external_only cannot be used with --label_placement above_feature")
    if args.collinear_min_anchors <= 0:
        parser.error("--collinear_min_anchors must be > 0")
    if args.collinear_max_unit_gap < 0:
        parser.error("--collinear_max_unit_gap must be >= 0")
    if args.collinear_max_diagonal_drift < 0:
        parser.error("--collinear_max_diagonal_drift must be >= 0")
    if args.collinear_max_conflicts_in_merge_gap < 0:
        parser.error("--collinear_max_conflicts_in_merge_gap must be >= 0")
    if args.collinear_max_paralog_links_per_orthogroup <= 0:
        parser.error("--collinear_max_paralog_links_per_orthogroup must be > 0")
    return args




def linear_main(cmd_args) -> None:
    """
    Main function for generating linear genome diagrams.

    This function orchestrates the creation of linear genome diagrams by processing
    input GenBank files and optional BLAST comparison files. It leverages various
    configurations and Configurator provided via command-line arguments to customize the
    visualization, such as color schemes, genome feature selections, and layout options.

    Args:
        cmd_args (argparse.Namespace): Command-line arguments parsed by argparse, providing
                                       specifications for input files, output formats, and
                                       visualization Configurator.

    The function performs the following key steps:
    - Loading and validating input files and Configurator.
    - Configuring the visualization canvas and feature Configurator.
    - Executing the plotting process to generate the linear diagrams.
    - Handling output in specified file formats.

    The final output includes linear genome diagrams in user-specified file formats,
    illustrating genomic features and optional BLAST comparison results.
    """
    session_request = parse_session_pre_args(cmd_args, mode="linear")
    if session_request is not None:
        with TemporaryDirectory(prefix="gbdraw-session-") as temp_dir:
            session = load_session(session_request.session_path)
            if render_canonical_session_if_present(
                session,
                mode="linear",
                output_override=session_request.output,
                format_override=session_request.format,
                save_session=session_request.save_session,
                session_output=session_request.session_output,
            ):
                return
            run_spec = session_to_cli_args(
                session,
                mode="linear",
                temp_dir=Path(temp_dir),
                output_override=session_request.output,
                format_override=session_request.format,
            )
            args = _get_args(list(run_spec.args))
            args._gbdraw_source_session = session
            args._gbdraw_collect_losat_cache = bool(
                session_request.save_session
                or session_request.session_output
                or _source_session_losat_entries(session)
            )
            run_result = run_linear_from_namespace(args)
            save_session_sidecar_if_requested(
                save_session=session_request.save_session,
                session_output=session_request.session_output,
                output_prefix=args.output,
                run_result=run_result,
                source_session=session,
                cli_invocation_args=run_spec.cli_invocation_args,
                file_bindings=run_spec.file_bindings,
            )
        return

    if '-i' in cmd_args or '--input' in cmd_args:
        logger.warning(
            "WARNING: The -i/--input option is deprecated and will be removed in a future version. Please use --gbk instead.")
    args: argparse.Namespace = _get_args(cmd_args)
    run_result = run_linear_from_namespace(args)
    save_session_sidecar_if_requested(
        save_session=bool(args.save_session or args.session_output),
        session_output=args.session_output,
        output_prefix=args.output,
        run_result=run_result,
        cmd_args=cmd_args,
    )


def run_linear_from_namespace(args: argparse.Namespace) -> DiagramRunResult:
    """Run linear rendering from an already parsed argparse namespace."""

    source_session = getattr(args, "_gbdraw_source_session", None)
    if not isinstance(source_session, Mapping):
        source_session = None
    source_losat_entries = _source_session_losat_entries(source_session)
    collect_losat_cache = bool(
        getattr(args, "_gbdraw_collect_losat_cache", False)
        or getattr(args, "save_session", False)
        or getattr(args, "session_output", None)
        or source_losat_entries
    )
    out_file_prefix: str = args.output
    blast_files: str = args.blast
    protein_blastp_mode: str = str(args.protein_blastp_mode or "none")
    losatp_bin: str = args.losatp_bin
    ncbi_blastp_bin: str | None = getattr(args, "ncbi_blastp_bin", None)
    losatp_threads: int | None = args.losatp_threads
    protein_blastp_max_hits: int = args.protein_blastp_max_hits
    protein_blastp_candidate_limit: int | None = args.protein_blastp_candidate_limit
    orthogroup_membership_mode: str = ORTHOGROUP_INFERENCE_VERSION
    orthogroup_member_max_hits: int = 5
    align_orthogroup_feature: str = str(args.align_orthogroup_feature or "").strip()
    collinear_unit_mode: str = str(args.collinear_unit_mode or "auto")
    collinear_anchor_mode: str = "rbh"
    collinear_search_scope: str = str(args.collinear_search_scope or "adjacent")
    collinear_color_mode: str = str(args.collinear_color_mode or "orientation")
    collinearity_params = LosslessCollinearityParameters(
        min_anchors=args.collinear_min_anchors,
        max_unit_gap=args.collinear_max_unit_gap,
        max_diagonal_drift=args.collinear_max_diagonal_drift,
        max_conflicts=args.collinear_max_conflicts_in_merge_gap,
    )
    color_table_path: str = args.table
    strandedness: bool = args.separate_strands
    resolve_overlaps: bool = args.resolve_overlaps
    dinucleotide: str = args.nt.upper()
    show_gc: bool = args.show_gc
    gc_content_mode: str | None = args.gc_content_mode
    gc_content_min_percent: Optional[float] = args.gc_content_min_percent
    gc_content_max_percent: Optional[float] = args.gc_content_max_percent
    gc_content_show_axis: bool | None = args.gc_content_show_axis
    gc_content_show_ticks: bool | None = args.gc_content_show_ticks
    gc_content_tick_interval: Optional[float] = args.gc_content_tick_interval
    gc_content_large_tick_interval: Optional[float] = args.gc_content_large_tick_interval
    gc_content_small_tick_interval: Optional[float] = args.gc_content_small_tick_interval
    gc_content_tick_font_size: Optional[float] = args.gc_content_tick_font_size
    manual_window: int = args.window
    manual_step: int = args.step
    depth_files: list[str] | None = args.depth
    depth_track_groups: list[list[str]] | None = args.depth_track
    depth_track_labels: list[str] | None = list(args.depth_track_label or []) or None
    depth_track_colors: list[str] | None = list(args.depth_track_color or []) or None
    depth_track_heights: list[str] | None = list(args.depth_track_height or []) or None
    depth_track_large_tick_intervals: list[str] | None = list(args.depth_track_large_tick_interval or []) or None
    depth_track_small_tick_intervals: list[str] | None = list(args.depth_track_small_tick_interval or []) or None
    depth_track_tick_font_sizes: list[str] | None = list(args.depth_track_tick_font_size or []) or None
    show_depth: bool = bool(args.show_depth or depth_files or depth_track_groups)
    depth_color: str | None = args.depth_color
    depth_height: Optional[float] = args.depth_height
    depth_window: Optional[int] = args.depth_window
    depth_step: Optional[int] = args.depth_step
    depth_share_axis: bool = bool(args.share_depth_axis)
    depth_min: Optional[float] = args.depth_min
    depth_max: Optional[float] = args.depth_max
    depth_normalize: bool | None = args.depth_normalize
    depth_show_axis: bool | None = args.depth_show_axis
    depth_show_ticks: bool | None = args.depth_show_ticks
    depth_tick_interval: Optional[float] = args.depth_tick_interval
    depth_large_tick_interval: Optional[float] = args.depth_large_tick_interval
    depth_small_tick_interval: Optional[float] = args.depth_small_tick_interval
    depth_tick_font_size: Optional[float] = args.depth_tick_font_size
    align_center: bool = args.align_center
    keep_definition_left_aligned: bool = args.keep_definition_left_aligned
    evalue: float = args.evalue
    legend: str = args.legend
    gc_height: Optional[float] = args.gc_height
    show_skew: bool = args.show_skew
    bitscore: float = args.bitscore
    identity: float = args.identity
    alignment_length: int = args.alignment_length
    pairwise_match_style: str = args.pairwise_match_style
    show_labels: str = args.show_labels
    label_whitelist: str = args.label_whitelist
    label_blacklist: str = args.label_blacklist
    qualifier_priority_path: str = args.qualifier_priority
    label_table_path: str = args.label_table
    feature_table_path: str = args.feature_table
    selected_features_set: str = args.features.split(',')
    feature_shapes = parse_feature_shape_overrides(args.feature_shape)
    feature_height: Optional[float] = args.feature_height
    comparison_height: Optional[float] = args.comparison_height

    out_formats: list[str] = parse_formats(args.format)
    out_formats = handle_output_formats(out_formats)
    user_defined_default_colors: str = args.default_colors
    scale_style: str = args.scale_style
    scale_stroke_color: Optional[str] = args.scale_stroke_color
    scale_stroke_width: Optional[float] = args.scale_stroke_width
    scale_font_size: Optional[float] = args.scale_font_size
    ruler_label_font_size: Optional[float] = args.ruler_label_font_size
    effective_ruler_label_font_size: Optional[float] = (
        ruler_label_font_size if ruler_label_font_size is not None else scale_font_size
    )
    scale_label_color: Optional[str] = args.ruler_label_color
    scale_interval: Optional[int] = args.scale_interval
    legend_box_size: Optional[float] = args.legend_box_size
    legend_font_size: Optional[float] = args.legend_font_size
    normalize_length = args.normalize_length
    if alignment_length < 0:
        raise ValidationError("alignment_length must be >= 0")
    if blast_files or protein_blastp_mode != "none":
        load_comparison = True
    else:
        load_comparison = False
    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette, load_comparison)
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
    feature_table: Optional[DataFrame] = read_feature_visibility_file(feature_table_path)
    annotation_options = (
        AnnotationOptions(sets=read_annotation_table(args.annotation_table, mode="linear"))
        if args.annotation_table
        else None
    )
    feature_visibility_rules = compile_feature_visibility_rules(feature_table)
    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')

    filtering_cfg = config_dict.setdefault("labels", {}).setdefault("filtering", {})
    if qualifier_priority_path:
        filtering_cfg["qualifier_priority_df"] = read_qualifier_priority_file(qualifier_priority_path)
    else:
        filtering_cfg["qualifier_priority_df"] = None
    if label_whitelist:
        filtering_cfg["whitelist_df"] = read_filter_list_file(label_whitelist)
    else:
        filtering_cfg["whitelist_df"] = None
    if label_table_path:
        filtering_cfg["label_override_df"] = read_label_override_file(label_table_path)
    else:
        filtering_cfg["label_override_df"] = None

    block_stroke_color: Optional[str] = args.block_stroke_color
    block_stroke_width: Optional[float] = args.block_stroke_width
    definition_font_size: Optional[float] = args.definition_font_size
    definition_line_styles = parse_definition_line_style_overrides(args.definition_line_style)
    definition_show_replicon: bool = bool(args.show_replicon)
    definition_show_accession: bool = not bool(args.hide_accession)
    definition_show_length: bool = not bool(args.hide_length)
    plot_title: str = str(args.plot_title or "").strip()
    plot_title_position: str = str(args.plot_title_position or "bottom").strip().lower()
    plot_title_font_size: Optional[float] = args.plot_title_font_size
    label_font_size: Optional[float] = args.label_font_size
    label_rendering: str = args.label_rendering
    label_placement: Optional[str] = args.label_placement
    label_rotation: Optional[float] = args.label_rotation
    track_layout: str = args.track_layout
    track_axis_gap: Optional[float] = args.track_axis_gap
    linear_track_slot_specs = args.linear_track_slot_specs
    linear_track_axis_index: int | None = args.linear_track_axis_index
    ruler_on_axis: bool = bool(args.ruler_on_axis)
    if ruler_on_axis and not (scale_style == "ruler" and track_layout in {"above", "below"}):
        logger.warning(
            "WARNING: --ruler_on_axis is ignored unless --scale_style ruler and --track_layout above|below are set."
        )
        ruler_on_axis = False
    axis_stroke_color: str = (
        args.axis_stroke_color
        if args.axis_stroke_color is not None
        else ("dimgray" if ruler_on_axis else "lightgray")
    )
    if ruler_on_axis and scale_stroke_color is None:
        scale_stroke_color = axis_stroke_color
    if ruler_on_axis and scale_label_color is None:
        scale_label_color = axis_stroke_color
    axis_stroke_width: Optional[float] = args.axis_stroke_width
    line_stroke_color: Optional[str] = args.line_stroke_color
    line_stroke_width: Optional[float] = args.line_stroke_width
    if plot_title_font_size is not None and float(plot_title_font_size) <= 0:
        raise ValidationError("plot_title_font_size must be > 0")
    if args.linear_label_spacing is not None and float(args.linear_label_spacing) <= 0:
        raise ValidationError("linear_label_spacing must be > 0")
    config_dict = modify_config_dict(
        config_dict,
        block_stroke_color=block_stroke_color,
        block_stroke_width=block_stroke_width,
        linear_axis_stroke_color=axis_stroke_color,
        linear_axis_stroke_width=axis_stroke_width,
        linear_definition_font_size=definition_font_size,
        linear_definition_line_styles=definition_line_styles or None,
        linear_definition_show_replicon=definition_show_replicon,
        linear_definition_show_accession=definition_show_accession,
        linear_definition_show_length=definition_show_length,
        label_font_size=label_font_size,
        linear_label_spacing=args.linear_label_spacing,
        label_rendering=label_rendering,
        label_placement=label_placement,
        label_rotation=label_rotation,
        line_stroke_color=line_stroke_color,
        line_stroke_width=line_stroke_width,
        show_gc=show_gc,
        gc_content_mode=gc_content_mode,
        gc_content_min_percent=gc_content_min_percent,
        gc_content_max_percent=gc_content_max_percent,
        gc_content_show_axis=gc_content_show_axis,
        gc_content_show_ticks=gc_content_show_ticks,
        gc_content_tick_interval=gc_content_tick_interval,
        gc_content_large_tick_interval=gc_content_large_tick_interval,
        gc_content_small_tick_interval=gc_content_small_tick_interval,
        gc_content_tick_font_size=gc_content_tick_font_size,
        show_skew=show_skew,
        show_depth=show_depth,
        depth_color=depth_color,
        depth_height=depth_height,
        depth_min=depth_min,
        depth_max=depth_max,
        depth_normalize=depth_normalize,
        depth_show_axis=depth_show_axis,
        depth_show_ticks=depth_show_ticks,
        depth_tick_interval=depth_tick_interval,
        depth_large_tick_interval=depth_large_tick_interval,
        depth_small_tick_interval=depth_small_tick_interval,
        depth_tick_font_size=depth_tick_font_size,
        depth_share_axis=depth_share_axis,
        show_labels=show_labels,
        resolve_overlaps=resolve_overlaps,
        linear_track_layout=track_layout,
        linear_track_axis_gap=track_axis_gap,
        linear_ruler_on_axis=ruler_on_axis,
        align_center=align_center,
        keep_definition_left_aligned=keep_definition_left_aligned,
        strandedness=strandedness,
        label_blacklist=label_blacklist,
        label_whitelist=label_whitelist,
        label_table=label_table_path,
        default_cds_height=feature_height,
        comparison_height=comparison_height,
        gc_height=gc_height,
        scale_style=scale_style,
        scale_stroke_color=scale_stroke_color,
        scale_label_color=scale_label_color,
        scale_stroke_width=scale_stroke_width,
        scale_font_size=scale_font_size,
        ruler_label_font_size=effective_ruler_label_font_size,
        scale_interval=scale_interval,
        legend_box_size=legend_box_size,
        legend_font_size=legend_font_size,
        pairwise_match_style=pairwise_match_style,
        normalize_length=normalize_length
        )

    def _normalize_list(values, target_len, fill_value=""):
        items = list(values or [])
        if len(items) > target_len:
            logger.error(
                "ERROR: Too many --record_id/--reverse_complement values (expected at most %s).", target_len
            )
            raise ValidationError(
                f"Too many --record_id/--reverse_complement values (expected at most {target_len})."
            )
        while len(items) < target_len:
            items.append(fill_value)
        return items

    def _parse_bool(value: str | None) -> bool:
        if value is None:
            return False
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "y", "on"}:
            return True
        if text in {"0", "false", "no", "n", "off", "", "none", "null", "-"}:
            return False
        logger.error("ERROR: Invalid reverse_complement value: %s", value)
        raise ValidationError(f"Invalid reverse_complement value: {value}")

    records_table = read_records_table(args.records_table) if args.records_table else None
    linear_source_paths: list[str]
    record_label_values: list[str]
    record_subtitle_values: list[str]
    region_arg_values: list[str]
    if records_table:
        if records_table.has_multi_record_placement:
            logger.info("Ignoring records table row/column values in linear mode.")
        if records_table.input_kind == "gbk":
            args.gbk = records_table.gbk_files
            args.gff = None
            args.fasta = None
            linear_source_paths = records_table.gbk_files
        else:
            args.gbk = None
            args.gff = records_table.gff_files
            args.fasta = records_table.fasta_files
            linear_source_paths = records_table.gff_files
        args.record_id = records_table.record_ids
        args.reverse_complement = [
            "1" if flag else "0" for flag in records_table.reverse_flags
        ]
        args.region = records_table.row_scoped_region_specs()
        records = _load_records_table_records(
            records_table,
            mode="linear",
            selected_features_set=selected_features_set,
            color_table=color_table,
            feature_table=feature_table,
            gbk_loader=load_gbks,
            gff_loader=load_gff_fasta,
        )
        record_label_values = records_table.record_labels
        record_subtitle_values = records_table.record_subtitles
        region_arg_values = args.region
    elif args.gbk:
        file_count = len(args.gbk)
        linear_source_paths = list(args.gbk)
        record_selectors = _normalize_list(args.record_id, file_count, "")
        reverse_flags_raw = _normalize_list(args.reverse_complement, file_count, "")
        reverse_flags = [_parse_bool(v) for v in reverse_flags_raw]
        records = load_gbks(
            args.gbk,
            "linear",
            load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
        record_label_values = list(args.record_label or [])
        record_subtitle_values = list(args.record_subtitle or [])
        region_arg_values = list(args.region or [])
    elif args.gff and args.fasta:
        file_count = len(args.gff)
        linear_source_paths = list(args.gff)
        record_selectors = _normalize_list(args.record_id, file_count, "")
        reverse_flags_raw = _normalize_list(args.reverse_complement, file_count, "")
        reverse_flags = [_parse_bool(v) for v in reverse_flags_raw]
        candidate_feature_types, keep_all_features = resolve_candidate_feature_types(
            selected_features_set,
            color_table=color_table,
            feature_visibility_table=feature_table,
        )
        records = load_gff_fasta(
            args.gff,
            args.fasta,
            "linear",
            candidate_feature_types,
            keep_all_features=keep_all_features,
            load_comparison=load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
        record_label_values = list(args.record_label or [])
        record_subtitle_values = list(args.record_subtitle or [])
        region_arg_values = list(args.region or [])
    else:
        logger.error("A critical error occurred with input file arguments.")
        raise ValidationError("Invalid input file arguments.")
    record_labels = record_label_values
    if record_labels:
        if len(record_labels) > len(records):
            logger.warning(
                "WARNING: More --record_label values were provided than records loaded; extra labels will be ignored."
            )
        for idx, label in enumerate(record_labels[: len(records)]):
            if label is None:
                continue
            label = str(label).strip()
            if not label:
                continue
            if getattr(records[idx], "annotations", None) is None:
                records[idx].annotations = {}
            records[idx].annotations["gbdraw_record_label"] = label
    record_subtitles = record_subtitle_values
    if record_subtitles:
        if len(record_subtitles) > len(records):
            logger.warning(
                "WARNING: More --record_subtitle values were provided than records loaded; extra subtitles will be ignored."
            )
        for idx, subtitle in enumerate(record_subtitles[: len(records)]):
            if subtitle is None:
                continue
            subtitle = str(subtitle).strip()
            if not subtitle:
                continue
            if getattr(records[idx], "annotations", None) is None:
                records[idx].annotations = {}
            records[idx].annotations["gbdraw_record_subtitle"] = subtitle
    region_specs = parse_region_specs(region_arg_values)
    if region_specs:
        try:
            records = apply_region_specs(records, region_specs, log=logger)
        except ValueError as exc:
            logger.error(f"ERROR: {exc}")
            raise ValidationError(str(exc)) from exc
        if blast_files:
            logger.warning(
                "WARNING: Region cropping is enabled; ensure BLAST coordinates match the cropped regions (and reverse complements if specified)."
            )
    linear_record_metadata = _linear_session_record_metadata(records, linear_source_paths)
    if protein_blastp_mode != "none" and len(records) < 2:
        raise ValidationError("--protein_blastp_mode requires at least two linear records.")
    depth_track_files = _record_major_depth_track_files_from_cli(
        depth_track_groups,
        record_count=len(records),
    )
    losatp_cache: LosatpCacheManager | None = None
    losatp_cache_entries: tuple[dict[str, object], ...] | None = None
    protein_extraction = None
    losat_cache_filenames: tuple[str, ...] = ()
    if protein_blastp_mode != "none" and collect_losat_cache:
        losatp_cache = LosatpCacheManager(
            source_losat_entries,
            threads_per_job=(
                losatp_threads
                if losatp_threads is not None
                else _source_session_threads_per_job(source_session)
            ),
            runtime_compatibility=_source_session_runtime_compatibility(source_session),
        )
        record_instance_keys = _record_instance_keys_for_web_losat(
            args,
            source_session=source_session,
            record_count=len(records),
        )
        protein_extraction = extract_web_stable_cds_proteins(
            records,
            record_instance_keys=record_instance_keys,
            feature_visibility_rules=feature_visibility_rules,
        )
        losat_cache_filenames = _linear_losat_cache_filenames(records)
    collinearity_comparisons: list[DataFrame] | None = None
    collinearity_orthogroups = None
    if protein_blastp_mode == "pairwise" and losatp_cache is not None:
        protein_blastp_result = build_pairwise_protein_blastp_comparisons(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            max_hits=protein_blastp_max_hits,
            candidate_limit=protein_blastp_candidate_limit,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=losat_cache_filenames,
        )
        collinearity_comparisons = protein_blastp_result.comparisons
        collinearity_orthogroups = protein_blastp_result.orthogroups
    elif protein_blastp_mode == "orthogroup" and losatp_cache is not None:
        protein_blastp_result = build_rbh_orthogroup_protein_blastp_comparisons(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_related_edges_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=losat_cache_filenames,
        )
        collinearity_comparisons = protein_blastp_result.comparisons
        collinearity_orthogroups = protein_blastp_result.orthogroups
    elif protein_blastp_mode == "collinear":
        collinearity_result = build_orthogroup_collinearity_blocks(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            params=collinearity_params,
            unit_mode=collinear_unit_mode,
            edge_mode=collinear_anchor_mode,
            search_scope=collinear_search_scope,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=losat_cache_filenames,
        )
        collinearity_orthogroups = collinearity_result.orthogroups
        collinearity_comparisons = convert_collinearity_blocks_to_comparisons(
            collinearity_result,
            records=records,
            color_mode=collinear_color_mode,
        )
    # Use raw records to avoid collapsing lengths when IDs are duplicated.
    longest_genome: int = max(len(record.seq) for record in records)
    cfg = GbdrawConfig.from_dict(config_dict)
    window, step = calculate_window_step(longest_genome, cfg, manual_window, manual_step)

    canvas = assemble_linear_diagram_from_records(
        records=records,
        blast_files=blast_files,
        config_dict=config_dict,
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=selected_features_set,
        feature_table=feature_table,
        feature_shapes=feature_shapes or None,
        output_prefix=out_file_prefix,
        legend=legend,
        dinucleotide=dinucleotide,
        window=window,
        step=step,
        depth_window=depth_window,
        depth_step=depth_step,
        depth_files=depth_files,
        depth_track_files=depth_track_files,
        depth_track_labels=depth_track_labels,
        depth_track_colors=depth_track_colors,
        depth_track_heights=depth_track_heights,
        depth_track_large_tick_intervals=depth_track_large_tick_intervals,
        depth_track_small_tick_intervals=depth_track_small_tick_intervals,
        depth_track_tick_font_sizes=depth_track_tick_font_sizes,
        linear_track_slots=linear_track_slot_specs,
        linear_track_axis_index=linear_track_axis_index,
        annotation_options=annotation_options,
        plot_title=plot_title,
        plot_title_position=plot_title_position,
        plot_title_font_size=plot_title_font_size,
        protein_comparisons=collinearity_comparisons,
        orthogroups=collinearity_orthogroups,
        protein_blastp_mode="none" if collinearity_comparisons is not None else protein_blastp_mode,
        losatp_bin=losatp_bin,
        ncbi_blastp_bin=ncbi_blastp_bin,
        losatp_threads=losatp_threads,
        protein_blastp_max_hits=protein_blastp_max_hits,
        protein_blastp_candidate_limit=protein_blastp_candidate_limit,
        orthogroup_membership_mode=orthogroup_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        collinear_max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
        align_orthogroup_feature=align_orthogroup_feature or None,
        pairwise_match_style=pairwise_match_style,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        cfg=cfg,
    )
    interactive_context = None
    if INTERACTIVE_SVG_FORMAT in out_formats:
        interactive_context = _build_interactive_svg_context(
            records,
            selected_features_set,
            getattr(canvas, "_gbdraw_orthogroups", None) or collinearity_orthogroups,
            feature_visibility_rules=feature_visibility_rules,
            color_table=color_table,
            default_colors=default_colors,
            annotations=annotation_options,
        )
        save_figure(
            canvas,
            out_formats,
            interactive_context=interactive_context,
        )
    else:
        save_figure(canvas, out_formats)

    if losatp_cache is not None:
        losatp_cache_entries = losatp_cache.session_entries()

    rendered_svg = make_rendered_svg(out_file_prefix, Path(str(out_file_prefix)).name)
    track_slot_geometry_records = collect_track_slot_geometry_records(
        canvas,
        result_index=0,
        result_name=rendered_svg.svg_path.name,
    )

    canonical_config = copy.deepcopy(config_dict)
    canonical_filtering = canonical_config["labels"]["filtering"]
    qualifier_priority_table = canonical_filtering.pop("qualifier_priority_df", None)
    label_whitelist_table = canonical_filtering.pop("whitelist_df", None)
    label_override_table = canonical_filtering.pop("label_override_df", None)
    request_prefix = Path(rendered_svg.output_prefix).name
    canonical_request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(record)) for record in records
        ),
        options=DiagramOptions(
            config=canonical_config,
            colors=ColorOptions(
                color_table=color_table,
                default_colors=default_colors,
                default_colors_palette=palette,
            ),
            tracks=TrackOptions(
                linear_track_slots=linear_track_slot_specs,
                linear_track_axis_index=linear_track_axis_index,
            ),
            annotations=annotation_options,
            output=OutputOptions(
                output_prefix=request_prefix,
                legend=legend,
                plot_title_position=plot_title_position,
            ),
            selected_features_set=tuple(selected_features_set),
            feature_visibility_table=feature_table,
            label_whitelist_table=label_whitelist_table,
            qualifier_priority_table=qualifier_priority_table,
            label_override_table=label_override_table,
            feature_shapes=feature_shapes or None,
            dinucleotide=dinucleotide,
            window=window,
            step=step,
            depth_window=depth_window,
            depth_step=depth_step,
            depth_files=tuple(depth_files) if depth_files else None,
            depth_track_files=depth_track_files,
            depth_track_labels=depth_track_labels,
            depth_track_colors=depth_track_colors,
            depth_track_heights=depth_track_heights,
            depth_track_large_tick_intervals=depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=depth_track_tick_font_sizes,
            plot_title=plot_title or None,
            plot_title_font_size=plot_title_font_size,
            blast_files=tuple(blast_files) if blast_files else None,
            protein_comparisons=(
                tuple(collinearity_comparisons)
                if collinearity_comparisons is not None
                else None
            ),
            orthogroups=collinearity_orthogroups,
            protein_blastp_mode=(
                "none" if collinearity_comparisons is not None else protein_blastp_mode
            ),
            pairwise_match_style=pairwise_match_style,
            collinearity_params=collinearity_params,
            collinearity_unit_mode=collinear_unit_mode,
            collinearity_anchor_mode=collinear_anchor_mode,
            collinearity_search_scope=collinear_search_scope,
            collinearity_color_mode=collinear_color_mode,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            protein_blastp_max_hits=protein_blastp_max_hits,
            protein_blastp_candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            collinear_max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
            align_orthogroup_feature=align_orthogroup_feature or None,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        ),
        output=RenderOutputRequest(
            output_prefix=request_prefix,
            formats=tuple(out_formats),
            overwrite=True,
        ),
    )

    return DiagramRunResult(
        mode="linear",
        render_formats=tuple(out_formats),
        outputs=(rendered_svg,),
        feature_metadata=tuple(interactive_context.features) if interactive_context else (),
        losat_cache_entries=losatp_cache_entries,
        linear_record_metadata=linear_record_metadata,
        run_metadata=build_track_slot_geometry_run_metadata(
            mode="linear",
            records=track_slot_geometry_records,
        ),
        canonical_request=canonical_request,
    )


if __name__ == "__main__":
    # This gets all arguments passed to the script, excluding the script name
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    linear_main(main_args)
