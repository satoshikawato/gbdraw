"""Bounded browser diagnostics; native exceptions and validators stay native.

Only this adapter reads native validation templates. It never compiles rules or
repairs inputs. Regex identity comes exclusively from re.error / explicit cause.
"""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
import json
import logging
import re
from typing import Iterator

from gbdraw.exceptions import ComparisonIdentityError, InputFileError, ParseError, ValidationError

OPERATIONS = frozenset("""unknown generate align feature-extraction export-svg export-png export-pdf evaluateRules readPdfFont
buildProteinLosatCacheKeys convertLosatNucleotideToDisplayTsv convertLosatpPairsToGenomicPayload
extractCdsProteinFasta extractFirstFasta generateLegendEntrySvg hydrateProteinLosatTsv
listGffFastaRecords listSequenceRecords measureLegendText promoteLegacyLosatpCache
regenerateDefinitionSvgs resolveLegacyProteinReferences resolveSimilarityAlignment validateConfigOverrides""".split())
STAGES = frozenset("""unknown initialization resource-staging request-validation helper
rule-validation render result-admission cleanup export-capture export-conversion font-validation""".split())

# Protocol identifiers, never values supplied by a document.
FIELDS = frozenset("""
pattern record_selector region start end sourceStart sourceEnd recordLength
recordIndex queryIndex subjectIndex depth min_depth max_depth window step tick
font_size plot_title_font_size height large_tick_interval small_tick_interval tick_font_size inner_gap_px outer_gap_px radius width spacing arrow_head_length_ratio arrow_shaft_width_ratio keep_definition_left_aligned color action feature_type
qualifier value record_id label_text config configOverrides records anchors schema
recordKey groupId direction sourceStrand role blast files input comparison
protein_blastp_max_hits orthogroup_member_max_hits bitscore evalue identity
alignment_length collinear_min_anchors collinear_max_gene_gap collinear_block_merge_gap
collinear_singleton_merge_gap collinear_max_diagonal_drift collinear_gap_penalty
collinear_nearby_duplicate_window collinear_constant_anchor_score
collinear_infer_orthogroups collinear_min_score collinear_min_block_span
record_gap_px record_axis_height depth_window depth_step depth_min depth_max
min_gc max_gc gc_tick_interval gc_axis_font_size depth_tick_interval depth_axis_font_size
protein_blastp_mode protein_blastp_candidate_limit collinear_search_scope collinear_unit_mode collinear_anchor_mode collinear_merge_orientation collinear_color_mode orthogroup_membership_mode collinear_max_unit_gap collinear_max_conflicts collinear_max_paralog_links_per_orthogroup circular_multi_record_size_mode linear_track_layout linear_label_placement set_id anchor_slot side renderer lane_gap_px padding_px cover_anchor overflow layer z axis match_height source fasta gff annotations featurePlacements output_prefix
""".split())

# Native, fixed validation clauses -> bounded correction identifiers.
_CONSTRAINTS = {
    "must be a boolean": "BOOLEAN", "must be an integer": "INTEGER",
    "must be non-negative": "NONNEGATIVE", "must be finite": "FINITE",
    "must be a positive integer": "POSITIVE_INTEGER",
    "must be a positive integer or None": "POSITIVE_OR_AUTO",
    "must be a non-negative integer": "NONNEGATIVE_INTEGER",
    "must be a finite number >= 0": "NONNEGATIVE",
    "must be a finite number > 0 or None": "POSITIVE_OR_AUTO",
    "must be a finite non-negative number": "NONNEGATIVE",
    "must be > 0": "POSITIVE", "must be >= 0": "NONNEGATIVE",
    "must be > 0 or None": "POSITIVE_OR_AUTO",
    "must be a finite value >= 0": "NONNEGATIVE",
    "must be a finite value between 0 and 100": "PERCENT",
    "must be an array": "ARRAY", "must be an object": "OBJECT",
    "has invalid fields": "FIELDS", "must be a non-empty string": "REQUIRED",
    "must be non-empty text without NUL": "REQUIRED",
    "must be an integer greater than or equal to 0": "NONNEGATIVE_INTEGER",
    "must be an integer greater than or equal to 1": "POSITIVE_INTEGER",
    "must be <= max_depth": "ORDER", "must be <= max_gc": "ORDER",
}

_EXACT = {
    "A record input source resolved to no records.": ("NO_RECORDS", {}),
    "A request requires at least one RecordInput.": ("INPUT_REQUIRED", {}),
    "DepthTrackInput.source must include at least one source.": ("DEPTH_INVALID", {"field": "source", "reason": "REQUIRED"}),
    "DepthTrackInput.source must include at least one non-None source.": ("DEPTH_INVALID", {"field": "source", "reason": "REQUIRED"}),
    "DepthTrackInput.source path must not be empty.": ("DEPTH_INVALID", {"field": "source", "reason": "REQUIRED"}),
    "No records found": ("NO_RECORDS", {}),
    "A resolved record collection cannot be empty.": ("NO_RECORDS", {}),
    "Specify either GenBank paths or GFF3/FASTA paths.": ("INPUT_REQUIRED", {}),
    "GFF3 and FASTA path counts must match.": ("FASTA_REQUIRED", {}),
    "GFF3 protein extraction requires a FASTA path.": ("FASTA_REQUIRED", {}),
    "An explicit display start cannot be combined with a crop.": ("REGION_INVALID", {"reason": "CROP_START_CONFLICT"}),
    "Depth table contains missing reference names.": ("DEPTH_INVALID", {"reason": "REFERENCE_REQUIRED"}),
    "Depth positions must be 1-based positive integers.": ("DEPTH_INVALID", {"reason": "POSITIVE_INTEGER"}),
    "Depth values must be non-negative.": ("DEPTH_INVALID", {"reason": "NONNEGATIVE"}),
    "Depth table must contain at least 3 columns.": ("TABLE_INVALID", {"columnCount": 3}),
    "The canonical Web request did not produce an SVG.": ("RESULT_INVALID", {}),
    "The staged Web render workspace is incomplete.": ("RESOURCE_INVALID", {"reason": "WORKSPACE"}),
    "The temporary staged Web render workspace could not be cleaned up.": ("CLEANUP_FAILED", {}),
    "The temporary Web render workspace could not be cleaned up.": ("CLEANUP_FAILED", {}),
    "collinear_anchor_mode must be one of: all, one_to_one, rbh": ("INPUT_INVALID", {"field": "collinear_anchor_mode", "reason": "COLLINEAR_ANCHOR_MODE"}),
    "collinear_color_mode must be one of: average_identity, orientation, orientation_identity": ("INPUT_INVALID", {"field": "collinear_color_mode", "reason": "COLLINEAR_COLOR_MODE"}),
    "collinear_merge_orientation must be one of: strand, order, either": ("INPUT_INVALID", {"field": "collinear_merge_orientation", "reason": "COLLINEAR_MERGE_ORIENTATION"}),
    "collinear_search_scope must be one of: adjacent, all": ("INPUT_INVALID", {"field": "collinear_search_scope", "reason": "ADJACENT_ALL"}),
    "Protein FASTA contains duplicate transport IDs.": ("COMPARISON_INPUT", {"reason": "UNIQUE_IDS"}),
    "Protein FASTA and protein map contain different transport IDs.": ("COMPARISON_INPUT", {"reason": "MATCH_IDS"}),
    "Reverse display transform requires a positive length.": ("REGION_INVALID", {"field": "recordLength", "reason": "POSITIVE_INTEGER"}),
    "Protein records must be a JSON array.": ("HELPER_PROTOCOL", {"field": "records", "reason": "ARRAY"}),
    "Legacy protein references must be a JSON array.": ("HELPER_PROTOCOL", {"reason": "ARRAY"}),
    "Protein record payloads must be JSON objects.": ("HELPER_PROTOCOL", {"field": "records", "reason": "OBJECT"}),
    "LOSATP blastp conversion payload must be an object with 'records' and 'pairs' lists.": ("HELPER_PROTOCOL", {"reason": "FIELDS"}),
    "LOSATP blastp conversion payload must contain 'records' and 'pairs' lists.": ("HELPER_PROTOCOL", {"reason": "FIELDS"}),
    "Unknown rule kind": ("HELPER_PROTOCOL", {}),
    "Unknown PDF font": ("HELPER_PROTOCOL", {}),
}

# Anchored producer templates; private interpolations are discarded, not shown.
_TEMPLATES = (
    (r"Unknown config override path: [\s\S]*\.", "INPUT_INVALID", "UNKNOWN_CONFIG_PATH"),
    (r"No matching FASTA record found for GFF record [\s\S]*\. Please ensure that all GFF records have corresponding FASTA entries\.", "FASTA_REQUIRED", "GFF_FASTA_MATCH"),
    (r"Unsupported LOSATP blastp mode: [\s\S]*", "INPUT_INVALID", "BLASTP_MODE"),
    (r"LOSATP (?:record|pair) payload #[0-9]+ (?:must be an object\.|is missing [\s\S]*|has an invalid raw TSV range\.|references missing [\s\S]*)", "HELPER_PROTOCOL", "FIELDS"),
    (r"LOSATP record payload contains duplicate recordIndex [0-9]+\.", "HELPER_PROTOCOL", "UNIQUE_IDS"),
    (r"Invalid action in feature visibility table at row ([0-9]+): '[\s\S]*'\. Use show, off, or exclude_matching\.[\s\S]*", "TABLE_INVALID", "VISIBILITY_ACTION"),
    (r"Record selector #[0-9]+ is out of range \(loaded ([0-9]+) record\(s\)\)\.", "RECORD_SELECTION", "OUT_OF_RANGE"),
    (r"Record selector '[\s\S]*' did not match any record ID\.", "RECORD_SELECTION", "NO_MATCH"),
    (r"Record selector '[\s\S]*' matched multiple records\. Use #index to disambiguate\.", "RECORD_SELECTION", "AMBIGUOUS"),
    (r"Invalid record selector '[\s\S]*'\. Use #<number> or record_id\.", "RECORD_SELECTION", "SELECTOR_FORMAT"),
    (r"Record index must be >= 1 in selector '[\s\S]*'\.", "RECORD_SELECTION", "POSITIVE_INTEGER"),
    (r"RecordInput #([0-9]+) resolved no records\.", "NO_RECORDS", "REQUIRED"),
    (r"RecordInput #([0-9]+) requires exactly one record; resolved ([0-9]+)\.[\s\S]*", "RECORD_SELECTION", "SELECT_ONE"),
    (r"Region spec is (?:missing|empty)\.", "REGION_INVALID", "REQUIRED"),
    (r"Invalid region spec: '[\s\S]*'\. Expected format:[\s\S]*", "REGION_INVALID", "REGION_FORMAT"),
    (r"Region coordinates must be >= 1: '[\s\S]*'\.", "REGION_INVALID", "POSITIVE_INTEGER"),
    (r"(?:Invalid record index|Record index must be >= 1) in region spec '[\s\S]*'\.[\s\S]*", "REGION_INVALID", "SELECTOR_FORMAT"),
    (r"Circular track slot '[\s\S]*' cannot fit inside between [\s\S]*", "TRACK_LAYOUT", "CANNOT_FIT"),
    (r"Preferred numeric group '[\s\S]*' cannot fit inside between [\s\S]*", "TRACK_LAYOUT", "CANNOT_FIT"),
    (r"Depth file '[\s\S]*' must contain at least 3 tab-separated columns\.", "TABLE_INVALID", "THREE_COLUMNS"),
    (r"Depth file references do not match record '[\s\S]*'\. Available references:[\s\S]*", "DEPTH_INVALID", "REFERENCE_MISMATCH"),
    (r"(?:Unable to (?:read|parse) depth file|Depth file does not exist:|Depth path is not a file:)[\s\S]*", "INPUT_UNREADABLE", "READ"),
    (r"Unsupported format: [\s\S]*", "INPUT_INVALID", "FORMAT"),
    (r"No CDS proteins found in [\s\S]*", "COMPARISON_INPUT", "NO_PROTEINS"),
    (r"Record ID not found: [\s\S]*", "RECORD_SELECTION", "NO_MATCH"),
    (r"Record index out of range: [0-9]+", "RECORD_SELECTION", "OUT_OF_RANGE"),
)


def _chain(error: BaseException) -> Iterator[BaseException]:
    seen = set()
    while error is not None and id(error) not in seen and len(seen) < 8:
        seen.add(id(error))
        yield error
        error = error.__cause__  # implicit context is not explicit authority


def _regex_reason(error: re.error) -> str:
    msg = error.msg
    if not isinstance(msg, str):
        return "SYNTAX_ERROR"
    reasons = {
        "unterminated character set": "UNTERMINATED_SET",
        "missing ), unterminated subpattern": "UNTERMINATED_GROUP",
        "unbalanced parenthesis": "UNBALANCED_GROUP",
        "nothing to repeat": "NOTHING_TO_REPEAT",
        "multiple repeat": "MULTIPLE_REPEAT",
        "global flags not at the start of the expression": "FLAGS_POSITION",
        "look-behind requires fixed-width pattern": "LOOKBEHIND_WIDTH",
    }
    if msg in reasons:
        return reasons[msg]
    # These prefixes are CPython's re.error.msg categories, never a traceback.
    for prefix, reason in (("bad escape", "INVALID_ESCAPE"), ("unknown extension", "UNKNOWN_EXTENSION"),
                           ("bad character range", "CHARACTER_RANGE"), ("invalid group reference", "GROUP_REFERENCE"),
                           ("redefinition of group name", "GROUP_NAME"), ("bad character in group name", "GROUP_NAME")):
        if msg.startswith(prefix):
            return reason
    return "SYNTAX_ERROR"


def _validation(error: BaseException) -> tuple[str, dict]:
    if not isinstance(error, (ValidationError, ParseError, ValueError, TypeError, KeyError)):
        return "VALIDATION_UNCLASSIFIED", {}
    field = getattr(error, "_web_error_field", None)
    if isinstance(field, str) and field in FIELDS:
        return "INPUT_INVALID", {"field": field, "reason": "FINITE"}
    message = str(error)
    if message in _EXACT:
        code, context = _EXACT[message]
        return code, dict(context)
    for template, code, reason in _TEMPLATES:
        match = re.fullmatch(template, message)
        if match:
            context = {"reason": reason}
            if reason == "UNKNOWN_CONFIG_PATH":
                context["field"] = "configOverrides"
            if reason == "VISIBILITY_ACTION":
                context.update(field="action", row=int(match[1]))
            if reason == "BLASTP_MODE":
                context["field"] = "protein_blastp_mode"
            if code == "RECORD_SELECTION" and reason == "OUT_OF_RANGE" and match.groups():
                context["recordCount"] = int(match[1])
            if reason == "SELECT_ONE":
                context.update(inputOrdinal=int(match[1]), recordCount=int(match[2]))
            return code, context
    for suffix, field, reason in (
        ("region start must not exceed end.", "region", "ORDER"),
        ("region exceeds recordLength.", "region", "RECORD_BOUNDS"),
        ("sourceEnd must not precede sourceStart.", "sourceEnd", "ORDER"),
        ("sourceStrand must be -1, 1, or null.", "sourceStrand", "STRAND"),
    ):
        if re.fullmatch(r"[\s\S]+\." + re.escape(suffix), message):
            return "INPUT_INVALID", {"field": field, "reason": reason}
    missing = re.fullmatch(r"Missing (record_id|feature_type|qualifier|value regex) token in (?:label override|feature visibility) table at row ([0-9]+)\.", message)
    if missing:
        return "TABLE_INVALID", {"field": "pattern" if missing[1] == "value regex" else missing[1], "row": int(missing[2]), "reason": "REQUIRED"}
    columns = re.fullmatch(r"Malformed line in (?:label override|feature visibility) file '[\s\S]*' at line ([0-9]+): expected ([0-9]+) columns\.", message)
    if columns:
        return "TABLE_INVALID", {"row": int(columns[1]), "columnCount": int(columns[2])}
    if re.fullmatch(r"Missing values in '[\s\S]*'\. See log for details\.", message):
        return "TABLE_INVALID", {"reason": "REQUIRED"}
    # Native field paths contain document IDs elsewhere. Only known final fields
    # and fixed constraint clauses are retained; no arbitrary path is public.
    for clause, reason in _CONSTRAINTS.items():
        match = re.fullmatch(r"([\s\S]+?) " + re.escape(clause) + r"(?:\.| or null)?", message)
        if match:
            field = match[1].rsplit(".", 1)[-1]
            if isinstance(field, str) and field in FIELDS:
                return "INPUT_INVALID", {"field": field, "reason": reason}
    return "VALIDATION_UNCLASSIFIED", {}


def serialize_web_error(error: BaseException, *, operation: str, stage: str) -> dict:
    """Produce identifiers only, without altering the native error or its cause."""
    operation = operation if isinstance(operation, str) and operation in OPERATIONS else "unknown"
    stage = stage if isinstance(stage, str) and stage in STAGES else "unknown"
    chain = list(_chain(error))
    context = {}
    code = "UNKNOWN"
    for cause in chain:
        if isinstance(cause, ComparisonIdentityError):
            code, context, stage = "COMPARISON_IDENTITY", {"reason": cause.reason}, "render"
            break
        if isinstance(cause, re.error):
            code, stage = "REGEX_SYNTAX", "rule-validation"
            context = {"reason": _regex_reason(cause), "positionUnit": "python-character"}
            if isinstance(cause.pos, int) and not isinstance(cause.pos, bool) and 0 <= cause.pos <= 10_000_000:
                context["position"] = cause.pos
            for parent in chain:
                if isinstance(parent, ParseError):
                    row = re.match(r"Invalid regex in (?:label override|label whitelist|feature visibility) table at row ([0-9]+):", str(parent))
                    if row:
                        context["row"] = int(row[1])
                        break
            break
    else:
        for cause in chain:
            candidate_code, candidate_context = _validation(cause)
            if candidate_code != "VALIDATION_UNCLASSIFIED":
                code, context = candidate_code, candidate_context
                break
        else:
            if isinstance(error, (ValidationError, ParseError, ValueError, TypeError, KeyError, json.JSONDecodeError)):
                code = "VALIDATION_UNCLASSIFIED"
            elif isinstance(error, (InputFileError, OSError, UnicodeError)):
                code = "INPUT_UNREADABLE"
    actual_stage = getattr(error, "_web_error_stage", None)
    if code not in {"REGEX_SYNTAX", "COMPARISON_IDENTITY"} and isinstance(actual_stage, str) and actual_stage in STAGES:
        stage = actual_stage
    if isinstance(error, json.JSONDecodeError):
        code, context = "HELPER_PROTOCOL", {"reason": "JSON_FORMAT"}
    if isinstance(error, StopIteration):
        code = "NO_RECORDS"
    context = {key: value for key, value in context.items()
               if (isinstance(value, str) and len(value) <= 80 and
                   (key != "reason" or value in {"EMPTY_ENDPOINT", "INDEX_ALIGNMENT", "SOURCE_INDEX", "SOURCE_VIEW_CONFLICT"} or code != "COMPARISON_IDENTITY"))
               or (type(value) is int and 0 <= value <= 10_000_000)}
    payload = {"code": code, "operation": operation, "stage": stage, "context": context}
    secondary = getattr(error, "_web_error_secondary", None)
    if isinstance(secondary, list) and secondary:
        payload["secondary"] = [{"code": "CLEANUP_FAILED", "stage": "cleanup"}
                                for item in secondary[:2] if item == {"code": "CLEANUP_FAILED", "stage": "cleanup"}]
    return payload


class _DiscardOutput:
    def write(self, text):
        return len(text)

    def flush(self):
        pass


@contextmanager
def private_web_execution():
    """Native logging remains available outside the browser adapter invocation."""
    previous = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    try:
        with redirect_stdout(_DiscardOutput()), redirect_stderr(_DiscardOutput()):
            yield
    finally:
        logging.disable(previous)


@contextmanager
def web_error_stage(stage: str):
    try:
        yield
    except Exception as error:
        if not hasattr(error, "_web_error_stage"):
            error._web_error_stage = stage
        raise
