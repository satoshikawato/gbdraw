export const PYTHON_HELPERS = `
import warnings
warnings.simplefilter('ignore', SyntaxWarning)
import json
from gbdraw.web_support.feature_metadata import (
    extract_features_from_genbank_json,
    extract_features_from_gff_fasta_json,
)
from gbdraw.web_support.request_render import (
    _render_staged_canonical_web_request_with_prepared_inputs as
    render_staged_canonical_web_request,
)
from gbdraw.session_request_codec import encode_canonical_typed_resource
from gbdraw.api.prepared import PreparedBiologicalInputCache
from gbdraw.web_support.error_adapter import private_web_execution, serialize_web_error, web_error_stage
from gbdraw.web_support.rule_matching import evaluate_rules_json
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
from gbdraw.web_support.similarity_alignment import resolve_similarity_alignment_json
from gbdraw.web_support.comparison_sequences import read_comparison_sequence_json

_WEB_LOSATP_FILTERED_HIT_CACHE = {}
_WEB_LOSATP_CONVERTED_PAYLOAD_CACHE = {}
_WEB_LOSATP_CACHE_ORDER = []
_WEB_LOSATP_CACHE_LIMIT = 64
_WEB_PREPARED_INPUT_CACHE = PreparedBiologicalInputCache()

def read_pdf_font(filename):
    import base64
    from importlib.resources import files
    allowed = {
        f"Liberation{family}-{style}.ttf"
        for family in ("Sans", "Serif", "Mono")
        for style in ("Regular", "Bold", "Italic", "BoldItalic")
    }
    if filename not in allowed:
        raise ValueError("Unknown PDF font")
    return json.dumps({"base64": base64.b64encode(files("gbdraw.data").joinpath(filename).read_bytes()).decode("ascii")})

def _web_losatp_cache_by_name(name):
    if name == "filtered":
        return _WEB_LOSATP_FILTERED_HIT_CACHE
    if name == "converted":
        return _WEB_LOSATP_CONVERTED_PAYLOAD_CACHE
    raise KeyError(name)

def _web_losatp_cache_get(name, key):
    cache = _web_losatp_cache_by_name(name)
    if key not in cache:
        return None
    marker = (name, key)
    try:
        _WEB_LOSATP_CACHE_ORDER.remove(marker)
    except ValueError:
        pass
    _WEB_LOSATP_CACHE_ORDER.append(marker)
    return cache[key]

def _web_losatp_cache_set(name, key, value):
    cache = _web_losatp_cache_by_name(name)
    marker = (name, key)
    try:
        _WEB_LOSATP_CACHE_ORDER.remove(marker)
    except ValueError:
        pass
    cache[key] = value
    _WEB_LOSATP_CACHE_ORDER.append(marker)
    while len(_WEB_LOSATP_CACHE_ORDER) > _WEB_LOSATP_CACHE_LIMIT:
        old_name, old_key = _WEB_LOSATP_CACHE_ORDER.pop(0)
        _web_losatp_cache_by_name(old_name).pop(old_key, None)

_WEB_LOSATP_CANONICAL_KEYS = {
    "collinearity-result": "collinearityResult",
    "orthogroup-result": "orthogroupResult",
}

def _web_losatp_payload_json(summary_json, canonical, canonical_resource_path):
    # The canonical result stays as its encoded bytes; parsing it back into
    # Python objects only to serialize it again dominated Worker memory.
    if canonical is None:
        return summary_json
    kind, content = canonical
    if canonical_resource_path:
        with open(str(canonical_resource_path), "wb") as handle:
            handle.write(content)
        member = json.dumps({"kind": kind, "size": len(content)})
        return "".join((summary_json[:-1], ', "canonicalResource": ', member, "}"))
    return "".join((summary_json[:-1], ", ", json.dumps(_WEB_LOSATP_CANONICAL_KEYS[kind]),
                    ": ", content.decode("utf-8"), "}"))

def _web_losatp_json_with_cache_stats(cached_payload, *, canonical_resource_path=None, **stats):
    summary_json, canonical = cached_payload
    payload = json.loads(summary_json)
    cache_payload = payload.get("cache")
    if not isinstance(cache_payload, dict):
        cache_payload = {}
    cache_payload.update(stats)
    payload["cache"] = cache_payload
    return _web_losatp_payload_json(json.dumps(payload), canonical, canonical_resource_path)

def _is_blank_or_js_nullish(value):
    if value is None:
        return True
    if type(value).__name__ in {"JsNull", "JsUndefined"}:
        return True
    try:
        return str(value).strip().lower() in {"", "null", "undefined", "none"}
    except Exception:
        return False

@private_web_execution()
def run_canonical_request_wrapper(
    request_json,
    resource_paths_json,
    workspace,
    diagnostics_enabled=False,
    resource_identities_json=None,
):
    try:
        with web_error_stage("request-validation"):
            payload = json.loads(str(request_json))
            resource_paths = json.loads(str(resource_paths_json))
            resource_identities = None
            if resource_identities_json is not None and type(resource_identities_json).__name__ not in {"JsNull", "JsUndefined"}:
                identity_text = str(resource_identities_json).strip()
                if identity_text and identity_text.lower() not in {"null", "undefined", "none"}:
                    resource_identities = json.loads(identity_text)
        diagnostics = {"timingsMs": {}, "metrics": {}} if diagnostics_enabled else None
        result = render_staged_canonical_web_request(
            payload,
            resource_paths=resource_paths,
            workspace=str(workspace),
            _diagnostics=diagnostics,
            _prepared_input_cache=(
                _WEB_PREPARED_INPUT_CACHE
                if resource_identities is not None
                else None
            ),
            _resource_identities=resource_identities,
        )
        if diagnostics is not None:
            for name in (
                "decode",
                "artifactCopy",
                "artifactValidation",
                "recordLoad",
                "preparation",
                "comparisonPreparation",
                "drawing",
                "interactivePreparation",
                "svgWrite",
                "svgReadback",
                "featureCatalog",
                "geometryMetadata",
            ):
                diagnostics["timingsMs"].setdefault(name, 0.0)
            for name in (
                "decodedResourceCacheHitCount",
                "decodedResourceCacheMissCount",
                "decodedResourceBuildCount",
                "parsedSourceCacheHitCount",
                "parsedSourceCacheMissCount",
                "parsedSourceParseCount",
                "resolvedRecordCacheHitCount",
                "resolvedRecordCacheMissCount",
                "resolvedRecordBuildCount",
                "interactiveContextCacheHitCount",
                "interactiveContextCacheMissCount",
                "interactiveContextBuildCount",
                "interactiveFeatureTraversalCount",
                "selectorSafetyScopeBuildCount",
                "preparedInputCacheEvictionCount",
                "preparedInputCacheRetainedBytes",
                "preparedInputCacheMutationViolationCount",
            ):
                diagnostics["metrics"].setdefault(name, 0)
            result["_diagnostics"] = diagnostics
        for item in result.get("results", []):
            content = item.get("content")
            if isinstance(content, str):
                item["content"] = content.encode("utf-8")
        metadata = result.get("metadata")
        if isinstance(metadata, dict):
            result["metadata"] = json.dumps(
                metadata,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        return result
    except Exception as e:
        return {"error": serialize_web_error(e, operation="generate", stage="render")}

def extract_first_fasta(path, fmt, region_spec=None, record_selector=None, reverse_flag=None):
    """Extract the first record as FASTA for LOSAT input."""
    from Bio import SeqIO
    from io import StringIO
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record
    try:
        fmt_map = {"genbank": "genbank", "fasta": "fasta"}
        if fmt not in fmt_map:
            return json.dumps({'error': serialize_web_error(ValueError(f'Unsupported format: {fmt}'), operation='extractFirstFasta', stage="helper")})
        records = list(SeqIO.parse(path, fmt_map[fmt]))
        if not records:
            return json.dumps({'error': serialize_web_error(ValueError('No records found'), operation='extractFirstFasta', stage="helper")})
        selector_raw = None
        if record_selector is not None:
            selector_raw = str(record_selector).strip()
            if not selector_raw or selector_raw.lower() in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
                selector_raw = None
        selector = parse_record_selector(selector_raw)
        if selector is None:
            records = [records[0]]
        else:
            records = select_record(records, selector)
        reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
        records = reverse_records(records, reverse)
        if region_spec:
            from gbdraw.io.regions import apply_region_specs, parse_region_specs
            records = apply_region_specs(records, parse_region_specs([region_spec]))
        record = records[0]
        handle = StringIO()
        SeqIO.write(record, handle, "fasta")
        return json.dumps({"fasta": handle.getvalue(), "record_id": record.id, "record_length": len(record.seq)})
    except StopIteration:
        return json.dumps({'error': serialize_web_error(ValueError('No records found'), operation='extractFirstFasta', stage="helper")})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='extractFirstFasta', stage="helper")})

def _normalize_web_record_selector(record_selector):
    if record_selector is None:
        return None
    selector_raw = str(record_selector).strip()
    if not selector_raw or selector_raw.lower() in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
        return None
    return selector_raw

def _normalize_web_view_transform(view_transform):
    if isinstance(view_transform, str):
        text = view_transform.strip()
        if text:
            try:
                view_transform = json.loads(text)
            except Exception:
                view_transform = {}
        else:
            view_transform = {}
    if not isinstance(view_transform, dict):
        view_transform = {}
    raw_length = view_transform.get("length", 0)
    try:
        length = int(raw_length)
    except Exception:
        length = 0
    raw_reverse = view_transform.get("reverse", False)
    if isinstance(raw_reverse, str):
        reverse = raw_reverse.strip().lower() in {"1", "true", "yes", "y", "on"}
    else:
        reverse = bool(raw_reverse)
    if reverse and length <= 0:
        raise ValueError("Reverse display transform requires a positive length.")
    return {"length": max(0, length), "reverse": reverse}

def _web_transform_cds_span(start, end, strand, view_transform):
    normalized = _normalize_web_view_transform(view_transform)
    start = int(start)
    end = int(end)
    if strand in (-1, 1, "-1", "1"):
        strand = int(strand)
    else:
        strand = None
    if not normalized["reverse"]:
        return start, end, strand
    length = int(normalized["length"])
    display_start = length - end
    display_end = length - start
    display_strand = -strand if strand in {-1, 1} else strand
    return display_start, display_end, display_strand

def _web_strand_symbol(strand):
    if strand in (-1, 1, "-1", "1"):
        return "+" if int(strand) == 1 else "-"
    return str(strand or "").strip()

def _web_read_record_coord_map(record):
    annotations = getattr(record, "annotations", {}) or {}
    try:
        base = int(annotations.get("gbdraw_coord_base", 1))
    except Exception:
        base = 1
    try:
        step = int(annotations.get("gbdraw_coord_step", 1))
    except Exception:
        step = 1
    if step == 0:
        step = 1
    return base, (1 if step > 0 else -1)

def _compute_web_feature_svg_id(record_id, feature_type, start, end, strand):
    from gbdraw.features.ids import compute_feature_hash_from_parts

    normalized_record_id = str(record_id or "")
    normalized_type = str(feature_type or "CDS")
    return compute_feature_hash_from_parts(
        normalized_type,
        int(start),
        int(end),
        strand,
        record_id=normalized_record_id or None,
    )

def _compute_web_feature_svg_id_from_parts(record_id, feature_type, parts):
    from gbdraw.features.ids import compute_feature_hash_from_location_parts

    normalized_record_id = str(record_id or "")
    normalized_type = str(feature_type or "CDS")
    return compute_feature_hash_from_location_parts(
        normalized_type,
        parts,
        record_id=normalized_record_id or None,
    )

def _normalize_web_feature_hash_parts(raw_parts):
    if raw_parts is None or raw_parts == "":
        return ()
    if isinstance(raw_parts, str):
        try:
            raw_parts = json.loads(raw_parts)
        except Exception:
            return ()
    if not isinstance(raw_parts, (list, tuple)):
        return ()
    parts = []
    for item in raw_parts:
        if isinstance(item, dict):
            start = item.get("start")
            end = item.get("end")
            strand = item.get("strand")
        elif isinstance(item, (list, tuple)) and len(item) >= 3:
            start, end, strand = item[0], item[1], item[2]
        else:
            continue
        try:
            start = int(start)
            end = int(end)
        except Exception:
            continue
        if strand in (-1, 1, "-1", "1"):
            strand = int(strand)
        else:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts)

def _display_feature_svg_id_from_data(data, display_start, display_end, display_strand, view_transform):
    normalized = _normalize_web_view_transform(view_transform)
    if not normalized["reverse"]:
        existing = data.get("view_feature_svg_id")
        if existing:
            return existing
    view_hash_parts = _normalize_web_feature_hash_parts(
        data.get("view_feature_hash_parts")
    )
    if view_hash_parts:
        display_hash_parts = [
            _web_transform_cds_span(start, end, strand, normalized)
            for start, end, strand in view_hash_parts
        ]
        return _compute_web_feature_svg_id_from_parts(
            data.get("record_id"),
            data.get("feature_type") or "CDS",
            display_hash_parts,
        )
    hash_start = data.get("start", display_start)
    hash_end = data.get("end", display_end)
    hash_strand = data.get("strand", display_strand)
    display_hash_start, display_hash_end, display_hash_strand = _web_transform_cds_span(
        hash_start,
        hash_end,
        hash_strand,
        normalized,
    )
    return _compute_web_feature_svg_id(
        data.get("record_id"),
        data.get("feature_type") or "CDS",
        display_hash_start,
        display_hash_end,
        display_hash_strand,
    )

def _web_feature_view_hash_parts_index(record):
    if record is None:
        return {}
    from gbdraw.core.record_metadata import _source_feature_index

    indexed = {}
    fallback_index = 0

    def walk(features):
        nonlocal fallback_index
        for feature in features or ():
            source_index = _source_feature_index(feature)
            resolved_index = fallback_index if source_index is None else source_index
            fallback_index += 1
            location = getattr(feature, "location", None)
            raw_parts = list(getattr(location, "parts", None) or [location])
            indexed.setdefault(
                int(resolved_index),
                tuple(
                    (
                        int(part.start),
                        int(part.end),
                        int(part.strand) if part.strand in (-1, 1) else None,
                    )
                    for part in raw_parts
                    if part is not None
                )
            )
            walk(getattr(feature, "sub_features", None))
    walk(record.features)
    return indexed

def _web_feature_view_hash_parts(record, source_feature_position):
    if source_feature_position is None:
        return ()
    return _web_feature_view_hash_parts_index(record).get(
        int(source_feature_position),
        (),
    )

def _load_single_linear_record_for_proteins(path, fmt, fasta_path=None, region_spec=None, record_selector=None, reverse_flag=None):
    from Bio import SeqIO
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record
    from gbdraw.io.regions import apply_region_specs, parse_region_specs

    selector = parse_record_selector(_normalize_web_record_selector(record_selector))
    reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
    if fmt == "genbank":
        records = list(SeqIO.parse(path, "genbank"))
        if not records:
            raise ValueError("No records found")
        records = select_record(records, selector) if selector is not None else [records[0]]
        records = reverse_records(records, reverse)
    elif fmt == "gff":
        if not fasta_path:
            raise ValueError("GFF3 protein extraction requires a FASTA path.")
        from gbdraw.io.genome import load_gff_fasta
        records = load_gff_fasta(
            [path],
            [fasta_path],
            selected_features_set=["CDS"],
            keep_all_features=True,
            record_selectors=[_normalize_web_record_selector(record_selector) or ""],
            reverse_flags=[reverse],
        )
        records = records[:1]
    else:
        raise ValueError(f"Unsupported format: {fmt}")
    if region_spec:
        records = apply_region_specs(records, parse_region_specs([region_spec]))
    if not records:
        raise ValueError("No records found")
    return records[0]

def _serialize_cds_protein(protein, record=None, view_hash_parts_index=None):
    coord_base, coord_step = _web_read_record_coord_map(record)
    source_feature_position = getattr(
        protein,
        "source_feature_position",
        protein.feature_index,
    )
    view_feature_hash_parts = (
        ()
        if source_feature_position is None
        else (
            _web_feature_view_hash_parts(record, source_feature_position)
            if view_hash_parts_index is None
            else view_hash_parts_index.get(int(source_feature_position), ())
        )
    )
    return {
        "protein_id": protein.protein_id,
        "record_index": protein.record_index,
        "feature_index": protein.feature_index,
        "record_id": protein.record_id,
        "start": protein.start,
        "end": protein.end,
        "strand": protein.strand,
        "label": protein.label,
        "protein_length": protein.protein_length,
        "source_protein_id": protein.source_protein_id,
        "feature_svg_id": protein.feature_svg_id,
        "view_feature_svg_id": getattr(protein, "view_feature_svg_id", None),
        "view_feature_hash_parts": [list(part) for part in view_feature_hash_parts],
        "gene": getattr(protein, "gene", None),
        "product": getattr(protein, "product", None),
        "note": getattr(protein, "note", None),
        "locus_tag": getattr(protein, "locus_tag", None),
        "gene_id": getattr(protein, "gene_id", None),
        "old_locus_tag": getattr(protein, "old_locus_tag", None),
        "db_xref": list(getattr(protein, "db_xref", ()) or ()),
        "gff_id": getattr(protein, "gff_id", None),
        "parent_ids": list(getattr(protein, "parent_ids", ()) or ()),
        "gene_parent_id": getattr(protein, "gene_parent_id", None),
        "feature_type": getattr(protein, "feature_type", "CDS"),
        "feature_hash_start": getattr(protein, "feature_hash_start", None),
        "feature_hash_end": getattr(protein, "feature_hash_end", None),
        "feature_hash_strand": getattr(protein, "feature_hash_strand", None),
        "feature_hash_parts": [list(part) for part in (getattr(protein, "feature_hash_parts", ()) or ())],
        "location_operator": getattr(protein, "location_operator", ""),
        "source_feature_position": getattr(protein, "source_feature_position", None),
        "same_location_ordinal": getattr(protein, "same_location_ordinal", None),
        "feature_analysis_id": getattr(protein, "feature_analysis_id", None),
        "display_alias": getattr(protein, "display_alias", None),
        "runtime_handle": getattr(protein, "runtime_handle", None),
        "aa_sha256": getattr(protein, "aa_sha256", None),
        "record_instance_key": getattr(protein, "record_instance_key", None),
        "record_analysis_id": getattr(protein, "record_analysis_id", None),
        "protein_set_hash": getattr(protein, "protein_set_hash", None),
        "runtime_binding_hash": getattr(protein, "runtime_binding_hash", None),
        "display_binding_hash": getattr(protein, "display_binding_hash", None),
        "coord_base": coord_base,
        "coord_step": coord_step,
        "coord_length": len(record.seq) if record is not None else 0,
    }

def extract_cds_protein_fasta(path, fmt, fasta_path=None, region_spec=None, record_selector=None, reverse_flag=None, record_index=None, record_instance_key=None, feature_visibility_table_path=None):
    """Extract CDS proteins and coordinate metadata for LOSATP blastp."""
    try:
        from gbdraw.analysis.protein_colinearity import (
            extract_protein_identity_manifest,
            proteins_to_fasta,
        )
        from gbdraw.features.visibility import compile_feature_visibility_rules, read_feature_visibility_file

        record = _load_single_linear_record_for_proteins(
            path,
            fmt,
            fasta_path=fasta_path,
            region_spec=region_spec,
            record_selector=record_selector,
            reverse_flag=reverse_flag,
        )
        feature_visibility_rules = None
        if not _is_blank_or_js_nullish(feature_visibility_table_path):
            feature_visibility_rules = compile_feature_visibility_rules(
                read_feature_visibility_file(feature_visibility_table_path)
            )
        record_index_offset = int(record_index) if record_index is not None else 0
        stable_record_key = record_instance_key
        if _is_blank_or_js_nullish(stable_record_key):
            stable_record_key = f"record-{record_index_offset + 1}"
        normalized_selector = _normalize_web_record_selector(record_selector) or None
        normalized_region = None if _is_blank_or_js_nullish(region_spec) else str(region_spec)
        result = extract_protein_identity_manifest(
            [record],
            record_instance_keys=[str(stable_record_key)],
            record_source_ids=[str(record.id)],
            record_selectors=[normalized_selector],
            regions=[normalized_region],
            record_index_offset=record_index_offset,
            feature_visibility_rules=feature_visibility_rules,
        )
        proteins = result.proteins_by_record[0] if result.proteins_by_record else []
        if not proteins:
            return json.dumps({'error': serialize_web_error(ValueError(f'No CDS proteins found in {record.id}'), operation='extractCdsProteinFasta', stage="helper")})
        view_hash_parts_index = _web_feature_view_hash_parts_index(record)
        protein_map = {
            protein.protein_id: _serialize_cds_protein(
                protein,
                record,
                view_hash_parts_index,
            )
            for protein in proteins
        }
        return json.dumps({
            "fasta": proteins_to_fasta(proteins),
            "record_id": record.id,
            "record_length": len(record.seq),
            "protein_count": len(proteins),
            "protein_map": protein_map,
            "identity_manifest": result.identity_manifest.to_dict(),
            "protein_set_hash": result.protein_set_hashes[0],
            "record_analysis_id": result.record_analysis_ids[0],
            "record_instance_key": result.record_instance_keys[0],
            "runtime_binding_hash": result.runtime_binding_hashes[0],
            "display_binding_hash": result.display_binding_hashes[0],
        })
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='extractCdsProteinFasta', stage="helper")})

def _build_web_cds_protein_map(raw_map):
    from gbdraw.analysis.protein_colinearity import CdsProtein

    protein_map = {}
    if not isinstance(raw_map, dict):
        return protein_map
    for protein_id, data in raw_map.items():
        if not isinstance(data, dict):
            continue
        strand = data.get("strand")
        strand = int(strand) if strand in (-1, 1, "-1", "1") else None
        kwargs = {
            "protein_id": str(data.get("protein_id") or protein_id),
            "record_index": int(data.get("record_index") or 0),
            "feature_index": int(data.get("feature_index") or 0),
            "record_id": str(data.get("record_id") or ""),
            "start": int(data.get("start") or 0),
            "end": int(data.get("end") or 0),
            "strand": strand,
            "label": str(data.get("label") or protein_id),
            "protein_length": int(data.get("protein_length") or 0),
            "sequence": str(data.get("sequence") or ""),
            "source_protein_id": data.get("source_protein_id"),
            "feature_svg_id": data.get("feature_svg_id"),
        }
        supported_fields = getattr(CdsProtein, "__dataclass_fields__", {})
        for optional_field in (
            "gene",
            "product",
            "note",
            "locus_tag",
            "gene_id",
            "old_locus_tag",
            "gff_id",
            "gene_parent_id",
            "view_feature_svg_id",
        ):
            if optional_field in supported_fields:
                kwargs[optional_field] = data.get(optional_field)
        if "feature_type" in supported_fields:
            kwargs["feature_type"] = str(data.get("feature_type") or "CDS")
        for optional_int_field in ("feature_hash_start", "feature_hash_end", "feature_hash_strand"):
            if optional_int_field in supported_fields:
                raw_value = data.get(optional_int_field)
                if raw_value is None or raw_value == "":
                    kwargs[optional_int_field] = None
                else:
                    kwargs[optional_int_field] = int(raw_value)
        if "feature_hash_parts" in supported_fields:
            kwargs["feature_hash_parts"] = _normalize_web_feature_hash_parts(data.get("feature_hash_parts"))
        for optional_int_field in (
            "source_feature_position",
            "same_location_ordinal",
        ):
            if optional_int_field in supported_fields:
                raw_value = data.get(optional_int_field)
                kwargs[optional_int_field] = (
                    int(raw_value) if raw_value is not None and raw_value != "" else None
                )
        for optional_string_field in (
            "location_operator",
            "feature_analysis_id",
            "display_alias",
            "runtime_handle",
            "aa_sha256",
            "record_instance_key",
            "record_analysis_id",
            "protein_set_hash",
            "runtime_binding_hash",
            "display_binding_hash",
        ):
            if optional_string_field in supported_fields:
                raw_value = data.get(optional_string_field)
                kwargs[optional_string_field] = (
                    str(raw_value) if raw_value is not None else None
                )
        for tuple_field in ("db_xref", "parent_ids"):
            if tuple_field in supported_fields:
                raw_values = data.get(tuple_field) or ()
                if isinstance(raw_values, (list, tuple)):
                    kwargs[tuple_field] = tuple(str(value) for value in raw_values if str(value).strip())
                else:
                    kwargs[tuple_field] = (str(raw_values),) if str(raw_values).strip() else ()
        protein_map[str(protein_id)] = CdsProtein(**kwargs)
    return protein_map

def _web_ordered_proteins_from_fasta(raw_map, fasta_text):
    from Bio import SeqIO
    from dataclasses import replace
    from io import StringIO

    protein_map = _build_web_cds_protein_map(raw_map)
    records = list(SeqIO.parse(StringIO(str(fasta_text or "")), "fasta"))
    record_ids = [str(record.id) for record in records]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("Protein FASTA contains duplicate transport IDs.")
    if set(record_ids) != set(protein_map):
        raise ValueError("Protein FASTA and protein map contain different transport IDs.")
    return [
        replace(protein_map[str(record.id)], sequence=str(record.seq))
        for record in records
    ]

def promote_legacy_losatp_cache_candidates(
    candidates_json,
    query_fasta,
    subject_fasta,
    query_protein_map_json,
    subject_protein_map_json,
    identity_manifest_json,
    expected_options_json,
):
    """Verify and copy one legacy schema-2 protein entry into schema 4."""
    try:
        from gbdraw.analysis.protein_colinearity import (
            promote_legacy_protein_raw_cache_entries,
        )

        raw_candidates = json.loads(str(candidates_json))
        query_map = json.loads(str(query_protein_map_json))
        subject_map = json.loads(str(subject_protein_map_json))
        manifest = json.loads(str(identity_manifest_json))
        options = json.loads(str(expected_options_json))
        candidate_indexes = []
        entries = []
        for relative_index, candidate in enumerate(raw_candidates if isinstance(raw_candidates, list) else []):
            if not isinstance(candidate, dict) or not isinstance(candidate.get("entry"), dict):
                continue
            candidate_indexes.append(int(candidate.get("candidateIndex", relative_index)))
            entries.append(candidate["entry"])

        scan = promote_legacy_protein_raw_cache_entries(
            entries,
            query_proteins=_web_ordered_proteins_from_fasta(query_map, query_fasta),
            subject_proteins=_web_ordered_proteins_from_fasta(subject_map, subject_fasta),
            query_fasta=str(query_fasta),
            subject_fasta=str(subject_fasta),
            identity_manifest=manifest,
            expected_args=options.get("args") or [],
            expected_program=str(options.get("program") or "blastp"),
            expected_outfmt=str(options.get("outfmt") or "6"),
        )
        rejections = [
            {
                "candidateIndex": candidate_indexes[rejection.candidate_index],
                "reason": rejection.reason,
            }
            for rejection in scan.rejections
        ]
        if scan.promotion is None:
            return json.dumps({"status": "no-match", "rejections": rejections})
        promotion = scan.promotion
        protein_id_map = {}
        old_rows = [
            line.split("\t")
            for line in str(entries[promotion.candidate_index].get("text") or "").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        new_rows = [
            line.split("\t")
            for line in promotion.rewritten_tsv.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        for old_row, new_row in zip(old_rows, new_rows):
            if len(old_row) < 2 or len(new_row) < 2:
                continue
            protein_id_map[old_row[0]] = new_row[0]
            protein_id_map[old_row[1]] = new_row[1]
        return json.dumps({
            "status": "promoted",
            "candidateIndex": candidate_indexes[promotion.candidate_index],
            "entry": promotion.entry,
            "text": promotion.rewritten_tsv,
            "proteinIdMap": protein_id_map,
            "rejections": rejections,
        })
    except Exception as error:
        return json.dumps({'status': 'error', 'error': serialize_web_error(error, operation='promoteLegacyLosatpCache', stage="helper")})

def resolve_legacy_protein_reference_map_json(
    protein_records_json,
    identity_manifest_json,
    reference_ids_json,
):
    """Resolve legacy Web p_r_ artifact references through current identity."""
    try:
        from gbdraw.analysis.protein_colinearity import (
            build_legacy_protein_reference_map,
            validate_protein_identity_manifest,
        )

        raw_records = json.loads(str(protein_records_json))
        if not isinstance(raw_records, list):
            raise ValueError("Protein records must be a JSON array.")
        manifest = validate_protein_identity_manifest(
            json.loads(str(identity_manifest_json))
        )
        reference_ids = json.loads(str(reference_ids_json))
        if not isinstance(reference_ids, list):
            raise ValueError("Legacy protein references must be a JSON array.")
        protein_maps = []
        for raw_record in raw_records:
            if not isinstance(raw_record, dict):
                raise ValueError("Protein record payloads must be JSON objects.")
            proteins = _web_ordered_proteins_from_fasta(
                raw_record.get("proteinMap") or {},
                raw_record.get("fasta") or "",
            )
            protein_maps.append({
                protein.protein_id: protein
                for protein in proteins
            })
        extraction = _build_web_protein_extraction(
            protein_maps,
            identity_manifest=manifest,
        )
        protein_id_map = build_legacy_protein_reference_map(
            extraction,
            [str(reference) for reference in reference_ids],
        )
        return json.dumps({
            "status": "resolved",
            "proteinIdMap": protein_id_map,
        })
    except Exception as error:
        return json.dumps({'status': 'error', 'error': serialize_web_error(error, operation='resolveLegacyProteinReferences', stage="helper")})

def build_protein_losat_cache_keys_json(
    identity_manifest_json,
    pairs_json,
):
    """Build directional schema-4 keys after validating the manifest once."""
    try:
        from gbdraw.analysis.protein_colinearity import (
            build_protein_losat_cache_key,
            build_protein_losat_pair_identity,
            validate_protein_identity_manifest,
        )

        manifest = validate_protein_identity_manifest(json.loads(str(identity_manifest_json)))
        pairs = json.loads(str(pairs_json))
        keys = []
        for pair in pairs:
            pair_identity = build_protein_losat_pair_identity(
                manifest,
                query_record_instance_key=str(pair["queryRecordInstanceKey"]),
                subject_record_instance_key=str(pair["subjectRecordInstanceKey"]),
            )
            options = pair["expectedOptions"]
            keys.append(build_protein_losat_cache_key(
                pair_identity,
                args=options.get("args") or [],
                program=str(options.get("program") or "blastp"),
                outfmt=str(options.get("outfmt") or "6"),
                search_context=options.get("searchContext"),
            ))
        return json.dumps({"keys": keys})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='buildProteinLosatCacheKeys', stage="helper")})

def main_session_table_text_to_search_frame_json(table_text, query_frame_json, subject_frame_json):
    """Rewrite one origin/main Session nucleotide table to the search frame."""
    from gbdraw.linear_comparison import reverse_endpoint_table_text

    endpoints = []
    for frame_json in (query_frame_json, subject_frame_json):
        frame = json.loads(str(frame_json))
        endpoints.append((int(frame.get("length") or 0), frame.get("reverse") is True))
    text = reverse_endpoint_table_text(str(table_text), endpoints[0], endpoints[1])
    return json.dumps({"tsv": text})

def hydrate_protein_losat_tsv_json(entry_json, identity_manifest_json):
    """Hydrate one internal schema-4 protein TSV for user download."""
    try:
        from gbdraw.analysis.protein_colinearity import hydrate_protein_losat_tsv

        entry = json.loads(str(entry_json))
        manifest = json.loads(str(identity_manifest_json))
        text = hydrate_protein_losat_tsv(entry, manifest)
        return json.dumps({
            "status": "ok",
            "text": text,
            "utf8Bytes": len(text.encode("utf-8")),
        })
    except Exception as error:
        return json.dumps({'status': 'error', 'error': serialize_web_error(error, operation='hydrateProteinLosatTsv', stage="helper")})

def _build_display_web_cds_protein_map(raw_map, view_transform):
    normalized = _normalize_web_view_transform(view_transform)
    if not normalized["reverse"]:
        return _build_web_cds_protein_map(raw_map)
    display_map = {}
    if not isinstance(raw_map, dict):
        return {}
    for protein_id, data in raw_map.items():
        if not isinstance(data, dict):
            continue
        display_data = dict(data)
        start, end, strand = _web_transform_cds_span(
            display_data.get("start", 0),
            display_data.get("end", 0),
            display_data.get("strand"),
            normalized,
        )
        view_feature_svg_id = _display_feature_svg_id_from_data(
            data,
            start,
            end,
            strand,
            normalized,
        )
        display_data["start"] = start
        display_data["end"] = end
        display_data["strand"] = strand
        display_data["view_feature_svg_id"] = view_feature_svg_id
        display_map[str(protein_id)] = display_data
    return _build_web_cds_protein_map(display_map)

def _build_web_protein_extraction(protein_maps, identity_manifest=None):
    from gbdraw.analysis.protein_colinearity import ProteinExtractionResult

    combined = {}
    max_record_index = -1
    for protein_map in protein_maps:
        combined.update(protein_map)
        for protein in protein_map.values():
            max_record_index = max(max_record_index, int(protein.record_index))
    proteins_by_record = [[] for _ in range(max_record_index + 1)]
    for protein in combined.values():
        proteins_by_record[int(protein.record_index)].append(protein)
    for record_proteins in proteins_by_record:
        record_proteins.sort(key=lambda protein: (int(protein.start), int(protein.end), int(protein.feature_index), str(protein.protein_id)))
    return ProteinExtractionResult(
        proteins_by_record=proteins_by_record,
        protein_map=combined,
        identity_manifest=identity_manifest,
    )

def _clean_json_scalar(value):
    try:
        import pandas as pd
        if pd.isna(value):
            return ""
    except Exception:
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value

def _dataframe_json_rows(df):
    rows = []
    for row in df.to_dict(orient="records"):
        rows.append({str(key): _clean_json_scalar(value) for key, value in row.items()})
    return rows

def convert_losatp_blastp_pairs_to_genomic_payload(
    pairs_path,
    raw_tsv_path,
    mode="pairwise",
    max_hits=5,
    bitscore=50,
    evalue="1e-2",
    identity=0,
    alignment_length=0,
    collinear_min_anchors=1,
    collinear_max_unit_gap=0,
    collinear_unit_mode="auto",
    collinear_color_mode="orientation",
    collinear_anchor_mode="rbh",
    collinear_max_diagonal_drift=0,
    collinear_max_conflicts_in_merge_gap=1,
    collinear_max_paralog_links_per_orthogroup=2,
    collinear_search_scope="adjacent",
    orthogroup_membership_mode="anchor_core_v1",
    orthogroup_member_max_hits=None,
    collinear_merge_orientation="either",
    collinear_infer_orthogroups=True,
    canonical_resource_path=None,
    explicit_display_pairs=False,
):
    """Convert LOSATP blastp outputs for pairwise display or orthogroups."""
    try:
        from io import StringIO
        import math
        import pandas as pd
        from gbdraw.analysis.collinearity import (
            LosslessCollinearityParameters,
            build_orthogroup_collinearity_blocks_from_hits,
            convert_collinearity_blocks_to_comparisons,
            convert_collinearity_blocks_to_pair_comparisons,
            normalize_collinearity_anchor_mode,
            normalize_collinearity_color_mode,
            normalize_collinearity_search_scope,
        )
        from gbdraw.analysis.collinearity_units import normalize_collinearity_unit_mode
        from gbdraw.analysis.protein_colinearity import (
            convert_pair_protein_hits_to_genomic_links,
            filter_protein_hits_by_thresholds,
            normalize_orthogroup_membership_mode,
            parse_losatp_outfmt6,
            select_rbh_orthogroup_edges_from_directional_hits,
            select_top_hits_per_query,
        )
        from gbdraw.io.comparisons import COMPARISON_COLUMNS

        with open(str(pairs_path), "r", encoding="utf-8") as pairs_handle:
            raw_payload = json.load(pairs_handle)
        if not isinstance(raw_payload, dict):
            raise ValueError("LOSATP blastp conversion payload must be an object with 'records' and 'pairs' lists.")
        raw_records = raw_payload.get("records")
        raw_pairs = raw_payload.get("pairs")
        if not isinstance(raw_records, list) or not isinstance(raw_pairs, list):
            raise ValueError("LOSATP blastp conversion payload must contain 'records' and 'pairs' lists.")
        normalized_mode = str(mode or "pairwise").strip().lower()
        if normalized_mode not in {"pairwise", "orthogroup", "collinear"}:
            raise ValueError(f"Unsupported LOSATP blastp mode: {mode!r}")
        normalized_max_hits = 5
        if normalized_mode == "pairwise":
            normalized_max_hits = int(5 if _is_blank_or_js_nullish(max_hits) else max_hits)
            if normalized_max_hits <= 0:
                raise ValueError("protein_blastp_max_hits must be > 0")

        normalized_membership_mode = "anchor_core_v1"
        normalized_member_max_hits = None
        if normalized_mode in {"orthogroup", "collinear"}:
            normalized_membership_mode = normalize_orthogroup_membership_mode(
                str(orthogroup_membership_mode or "anchor_core_v1")
            )
            normalized_member_max_hits = (
                None
                if _is_blank_or_js_nullish(orthogroup_member_max_hits)
                else int(orthogroup_member_max_hits)
            )
            if normalized_member_max_hits is not None and normalized_member_max_hits <= 0:
                raise ValueError("orthogroup_member_max_hits must be > 0 or None")

        normalized_collinear_unit_mode = "auto"
        normalized_collinear_color_mode = "orientation"
        normalized_collinear_anchor_mode = "rbh"
        normalized_collinear_search_scope = "adjacent"
        normalized_max_paralog_links = 2
        normalized_collinearity_params = LosslessCollinearityParameters()
        if normalized_mode == "collinear":
            if not isinstance(collinear_infer_orthogroups, bool):
                raise ValueError("collinear_infer_orthogroups must be a boolean")
            normalized_collinear_unit_mode = normalize_collinearity_unit_mode(
                str(collinear_unit_mode or "auto")
            )
            normalized_collinear_color_mode = normalize_collinearity_color_mode(
                str(collinear_color_mode or "orientation")
            )
            normalized_collinear_anchor_mode = normalize_collinearity_anchor_mode(
                str(collinear_anchor_mode or "rbh")
            )
            normalized_collinear_search_scope = normalize_collinearity_search_scope(
                str(collinear_search_scope or "adjacent")
            )
            normalized_max_paralog_links = int(
                2
                if _is_blank_or_js_nullish(collinear_max_paralog_links_per_orthogroup)
                else collinear_max_paralog_links_per_orthogroup
            )
            if normalized_max_paralog_links <= 0:
                raise ValueError(
                    "collinear_max_paralog_links_per_orthogroup must be > 0"
                )
            normalized_collinearity_params = LosslessCollinearityParameters(
                min_anchors=int(
                    1
                    if _is_blank_or_js_nullish(collinear_min_anchors)
                    else collinear_min_anchors
                ),
                max_unit_gap=int(
                    0
                    if _is_blank_or_js_nullish(collinear_max_unit_gap)
                    else collinear_max_unit_gap
                ),
                max_diagonal_drift=int(
                    0
                    if _is_blank_or_js_nullish(collinear_max_diagonal_drift)
                    else collinear_max_diagonal_drift
                ),
                max_conflicts=int(
                    1
                    if _is_blank_or_js_nullish(collinear_max_conflicts_in_merge_gap)
                    else collinear_max_conflicts_in_merge_gap
                ),
                merge_orientation=str(collinear_merge_orientation or "either").strip().lower(),
            )
            normalized_collinearity_params.validate()

        normalized_bitscore = float(bitscore)
        normalized_evalue = float(evalue)
        normalized_identity = float(identity)
        normalized_alignment_length = int(alignment_length)
        if not math.isfinite(normalized_bitscore) or normalized_bitscore < 0:
            raise ValueError("bitscore must be a finite value >= 0")
        if not math.isfinite(normalized_evalue) or normalized_evalue < 0:
            raise ValueError("evalue must be a finite value >= 0")
        if (
            not math.isfinite(normalized_identity)
            or normalized_identity < 0
            or normalized_identity > 100
        ):
            raise ValueError("identity must be a finite value between 0 and 100")
        if normalized_alignment_length < 0:
            raise ValueError("alignment_length must be >= 0")

        record_payloads = []
        for idx, record in enumerate(raw_records):
            if not isinstance(record, dict):
                raise ValueError(f"LOSATP record payload #{idx + 1} must be an object.")
            try:
                record_index = int(record["recordIndex"])
            except Exception as exc:
                raise ValueError(f"LOSATP record payload #{idx + 1} is missing a valid recordIndex.") from exc
            record_payloads.append(
                {
                    "record_index": record_index,
                    "record_id": str(record.get("recordId") or ""),
                    "protein_map": record.get("proteinMap") or {},
                    "protein_cache_key": str(record.get("proteinCacheKey") or ""),
                    "view_transform": _normalize_web_view_transform(record.get("viewTransform") or {}),
                }
            )

        pair_payloads = []
        with open(str(raw_tsv_path), "rb") as raw_tsv_handle:
            raw_tsv_handle.seek(0, 2)
            raw_tsv_size = raw_tsv_handle.tell()
        for idx, item in enumerate(raw_pairs):
            if not isinstance(item, dict):
                raise ValueError(f"LOSATP pair payload #{idx + 1} must be an object.")
            try:
                query_index = int(item["queryIndex"])
            except Exception as exc:
                raise ValueError(f"LOSATP pair payload #{idx + 1} is missing a valid queryIndex.") from exc
            try:
                subject_index = int(item["subjectIndex"])
            except Exception as exc:
                raise ValueError(f"LOSATP pair payload #{idx + 1} is missing a valid subjectIndex.") from exc
            try:
                raw_offset = int(item["rawTsvOffset"])
                raw_bytes = int(item["rawTsvBytes"])
            except Exception as exc:
                raise ValueError(f"LOSATP pair payload #{idx + 1} is missing a valid raw TSV range.") from exc
            if raw_offset < 0 or raw_bytes < 0 or raw_offset + raw_bytes > raw_tsv_size:
                raise ValueError(f"LOSATP pair payload #{idx + 1} has an invalid raw TSV range.")
            pair_index = int(item.get("pairIndex", min(query_index, subject_index)))
            cache_key = str(item.get("cacheKey") or "").strip()
            if not cache_key:
                raise ValueError(f"LOSATP pair payload #{idx + 1} is missing cacheKey.")
            pair_payloads.append(
                {
                    "pair_index": pair_index,
                    "query_index": query_index,
                    "subject_index": subject_index,
                    "display_pair": item.get("displayPair") is True,
                    "cache_key": cache_key,
                    "raw_tsv_offset": raw_offset,
                    "raw_tsv_bytes": raw_bytes,
                }
            )

        raw_tsv_stats = {
            "rawTsvEntryCount": len(pair_payloads),
            "rawTsvBytes": raw_tsv_size,
            "rawTsvLargestEntryBytes": max(
                (item["raw_tsv_bytes"] for item in pair_payloads),
                default=0,
            ),
        }
        derived_provenance = {
            "schema": 1,
            "upstreamRawKeys": sorted(
                {str(item["cache_key"]) for item in pair_payloads}
            ),
            "thresholds": {
                "bitscore": normalized_bitscore,
                "evalue": normalized_evalue,
                "identity": normalized_identity,
                "alignmentLength": normalized_alignment_length,
            },
        }
        if normalized_mode == "pairwise":
            derived_provenance["pairwise"] = {
                "maxHits": normalized_max_hits,
            }
        if normalized_mode in {"orthogroup", "collinear"}:
            derived_provenance["orthogroup"] = {
                "membershipMode": normalized_membership_mode,
                "memberMaxHits": normalized_member_max_hits,
            }
        if normalized_mode == "collinear":
            derived_provenance["collinear"] = {
                "unitMode": {
                    "requested": normalized_collinear_unit_mode,
                    "effectiveKinds": [],
                },
                "anchorMode": normalized_collinear_anchor_mode,
                "searchScope": normalized_collinear_search_scope,
                "inferOrthogroups": collinear_infer_orthogroups,
                "colorMode": normalized_collinear_color_mode,
                "parameters": {
                    "minAnchors": int(normalized_collinearity_params.min_anchors),
                    "maxUnitGap": int(normalized_collinearity_params.max_unit_gap),
                    "maxDiagonalDrift": int(
                        normalized_collinearity_params.max_diagonal_drift
                    ),
                    "maxConflicts": int(normalized_collinearity_params.max_conflicts),
                    "mergeOrientation": str(
                        normalized_collinearity_params.merge_orientation
                    ),
                    "maxParalogLinksPerOrthogroup": normalized_max_paralog_links,
                },
            }

        conversion_cache_key = (
            "derived-option-conformance-v1",
            normalized_mode,
            normalized_max_hits if normalized_mode == "pairwise" else None,
            str(normalized_bitscore),
            str(normalized_evalue),
            str(normalized_identity),
            str(normalized_alignment_length),
            (
                str(normalized_membership_mode),
                normalized_member_max_hits,
            ) if normalized_mode in {"orthogroup", "collinear"} else None,
            (
                int(normalized_collinearity_params.min_anchors),
                int(normalized_collinearity_params.max_unit_gap),
                normalized_collinear_unit_mode,
                normalized_collinear_color_mode,
                normalized_collinear_anchor_mode,
                str(normalized_collinearity_params.merge_orientation),
                int(normalized_collinearity_params.max_diagonal_drift),
                int(normalized_collinearity_params.max_conflicts),
                int(normalized_max_paralog_links),
                normalized_collinear_search_scope,
                collinear_infer_orthogroups,
                bool(explicit_display_pairs),
            ) if normalized_mode == "collinear" else None,
            tuple(
                (
                    item["record_index"],
                    item["protein_cache_key"],
                    int(item["view_transform"]["length"]),
                    bool(item["view_transform"]["reverse"]),
                )
                for item in sorted(record_payloads, key=lambda current: current["record_index"])
            ),
            tuple(
                (
                    item["pair_index"],
                    item["query_index"],
                    item["subject_index"],
                    item["display_pair"],
                    item["cache_key"],
                )
                for item in pair_payloads
            ),
        )
        cached_payload = _web_losatp_cache_get("converted", conversion_cache_key)
        if cached_payload is not None:
            return _web_losatp_json_with_cache_stats(
                cached_payload,
                canonical_resource_path=canonical_resource_path,
                convertedPayloadHit=True,
                filteredHitCacheHits=0,
                filteredHitCacheMisses=0,
                simultaneousParsedTables=0,
                **raw_tsv_stats,
            )

        protein_maps_by_record = {}
        for record in record_payloads:
            record_index = int(record["record_index"])
            if record_index in protein_maps_by_record:
                raise ValueError(f"LOSATP record payload contains duplicate recordIndex {record_index}.")
            protein_maps_by_record[record_index] = _build_display_web_cds_protein_map(
                record["protein_map"],
                record["view_transform"],
            )
        combined_protein_map = {}
        for protein_map in protein_maps_by_record.values():
            combined_protein_map.update(protein_map)
        extraction = _build_web_protein_extraction(protein_maps_by_record.values())

        filtered_cache_hits = 0
        filtered_cache_misses = 0
        pair_items = []
        for idx, item in enumerate(pair_payloads):
            query_index = int(item["query_index"])
            subject_index = int(item["subject_index"])
            if query_index not in protein_maps_by_record:
                raise ValueError(f"LOSATP pair payload #{idx + 1} references missing queryIndex {query_index}.")
            if subject_index not in protein_maps_by_record:
                raise ValueError(f"LOSATP pair payload #{idx + 1} references missing subjectIndex {subject_index}.")
            filter_cache_key = (
                item["cache_key"],
                str(normalized_bitscore),
                str(normalized_evalue),
                str(normalized_identity),
                str(normalized_alignment_length),
            )
            filtered = _web_losatp_cache_get("filtered", filter_cache_key)
            if filtered is None:
                with open(str(raw_tsv_path), "rb") as raw_tsv_handle:
                    raw_tsv_handle.seek(item["raw_tsv_offset"])
                    blast_text = raw_tsv_handle.read(item["raw_tsv_bytes"]).decode("utf-8")
                hits = parse_losatp_outfmt6(blast_text)
                filtered = filter_protein_hits_by_thresholds(
                    hits,
                    evalue=normalized_evalue,
                    bitscore=normalized_bitscore,
                    identity=normalized_identity,
                    alignment_length=normalized_alignment_length,
                )
                _web_losatp_cache_set("filtered", filter_cache_key, filtered.copy())
                filtered_cache_misses += 1
            else:
                filtered = filtered.copy()
                filtered_cache_hits += 1
            pair_items.append(
                {
                    "pair_index": item["pair_index"],
                    "query_index": query_index,
                    "subject_index": subject_index,
                    "hits": filtered,
                    "query_map": protein_maps_by_record[query_index],
                    "subject_map": protein_maps_by_record[subject_index],
                }
            )

        cache_stats = {
            "convertedPayloadHit": False,
            "filteredHitCacheHits": filtered_cache_hits,
            "filteredHitCacheMisses": filtered_cache_misses,
            "simultaneousParsedTables": len(pair_items),
            **raw_tsv_stats,
        }

        def _finalize_losatp_payload(payload, *, collinearity_result=None, canonical_content=None, canonical_kind=None):
            if collinearity_result is not None:
                anchors = [
                    anchor
                    for block in collinearity_result.blocks
                    for anchor in block.anchors
                ]
                anchors.extend(collinearity_result.unblocked_anchors)
                derived_provenance["collinear"]["unitMode"]["effectiveKinds"] = sorted({
                    str(unit_kind)
                    for anchor in anchors
                    for unit_kind in (
                        anchor.query_unit_kind,
                        anchor.subject_unit_kind,
                    )
                })
            payload["provenance"] = derived_provenance
            payload["cache"] = cache_stats
            cached_payload = (
                json.dumps(payload),
                None if canonical_content is None else (canonical_kind, canonical_content),
            )
            _web_losatp_cache_set("converted", conversion_cache_key, cached_payload)
            return _web_losatp_payload_json(*cached_payload, canonical_resource_path)

        if normalized_mode == "collinear":
            record_ids = []
            for record_proteins in extraction.proteins_by_record:
                if record_proteins:
                    record_ids.append(str(record_proteins[0].record_id))
                else:
                    record_ids.append("")
            anchor_mode = normalized_collinear_anchor_mode
            search_scope = normalized_collinear_search_scope
            max_paralog_links = normalized_max_paralog_links
            directional_tables = {}
            for item in pair_items:
                query_index = int(item["query_index"])
                subject_index = int(item["subject_index"])
                directional_tables[(query_index, subject_index)] = item["hits"]
            display_pairs = tuple(sorted({
                tuple(sorted((
                    int(item["query_index"]),
                    int(item["subject_index"]),
                )))
                for item in pair_payloads
                if item["display_pair"]
                and int(item["query_index"]) != int(item["subject_index"])
            }))
            use_display_pairs = bool(display_pairs) and (
                search_scope == "adjacent" or explicit_display_pairs
            )
            collinearity_result = build_orthogroup_collinearity_blocks_from_hits(
                directional_tables,
                extraction,
                params=normalized_collinearity_params,
                unit_mode=normalized_collinear_unit_mode,
                edge_mode=anchor_mode,
                search_scope=search_scope,
                infer_orthogroups=collinear_infer_orthogroups,
                orthogroup_membership_mode=normalized_membership_mode,
                orthogroup_member_max_hits=normalized_member_max_hits,
                max_paralog_links_per_orthogroup=max_paralog_links,
                comparison_pairs=(
                    display_pairs
                    if use_display_pairs
                    else None
                ),
            )
            color_mode = normalized_collinear_color_mode
            if use_display_pairs:
                display_pair_indices = {
                    tuple(sorted((
                        int(item["query_index"]),
                        int(item["subject_index"]),
                    ))): int(item["pair_index"])
                    for item in pair_payloads
                    if item["display_pair"]
                }
                converted_by_pair = convert_collinearity_blocks_to_pair_comparisons(
                    collinearity_result,
                    record_ids=record_ids,
                    color_mode=color_mode,
                    search_scope=search_scope,
                )
                converted_frames = [
                    (display_pair_indices[pair], converted_by_pair[pair])
                    for pair in display_pairs
                    if pair in converted_by_pair
                ]
            else:
                converted_frames = enumerate(convert_collinearity_blocks_to_comparisons(
                    collinearity_result,
                    record_ids=record_ids,
                    color_mode=color_mode,
                    search_scope=search_scope,
                ))
            converted_pairs = []
            for pair_index, converted in converted_frames:
                handle = StringIO()
                converted.loc[:, list(COMPARISON_COLUMNS)].to_csv(
                    handle,
                    sep=chr(9),
                    header=False,
                    index=False,
                    lineterminator=chr(10),
                )
                converted_pairs.append(
                    {
                        "pair_index": pair_index,
                        "tsv": handle.getvalue(),
                        "rows": _dataframe_json_rows(converted),
                        "hit_count": int(converted.shape[0]),
                    }
                )
            canonical_content = encode_canonical_typed_resource("result", collinearity_result)
            return _finalize_losatp_payload({"pairs": converted_pairs},
               collinearity_result=collinearity_result,
               canonical_content=canonical_content, canonical_kind="collinearity-result")

        if normalized_mode == "pairwise":
            converted_pairs = []
            for item in pair_items:
                display_hits = select_top_hits_per_query(
                    item["hits"],
                    max_hits=normalized_max_hits,
                )
                converted = convert_pair_protein_hits_to_genomic_links(
                    display_hits,
                    item["query_map"],
                    item["subject_map"],
                    orthogroups=None,
                )
                handle = StringIO()
                converted.loc[:, list(COMPARISON_COLUMNS)].to_csv(
                    handle,
                    sep=chr(9),
                    header=False,
                    index=False,
                    lineterminator=chr(10),
                )
                converted_pairs.append(
                    {
                        "pair_index": item["pair_index"],
                        "tsv": handle.getvalue(),
                        "rows": _dataframe_json_rows(converted),
                        "hit_count": int(converted.shape[0]),
                    }
                )
            return _finalize_losatp_payload({"pairs": converted_pairs, "orthogroups": []})

        hits_by_direction = {
            (item["query_index"], item["subject_index"]): item
            for item in pair_items
        }
        directional_tables = {
            pair: item["hits"]
            for pair, item in hits_by_direction.items()
        }
        display_pair_indices = {
            (item["query_index"], item["subject_index"]): item["pair_index"]
            for item in pair_payloads if item["display_pair"]
        }
        edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
            directional_tables,
            combined_protein_map,
            orthogroup_membership_mode=normalized_membership_mode,
            orthogroup_member_max_hits=normalized_member_max_hits,
            max_related_edges_per_orthogroup=normalized_max_paralog_links,
            comparison_pairs=tuple(display_pair_indices),
        )
        orthogroups = edge_selection.orthogroups

        converted_pairs = []
        for query_index, subject_index in sorted(edge_selection.adjacent_display_edges_by_pair):
            display_hits = edge_selection.adjacent_display_edges_by_pair[(query_index, subject_index)]
            forward = hits_by_direction.get((query_index, subject_index))
            if forward is None:
                continue
            converted = convert_pair_protein_hits_to_genomic_links(
                display_hits,
                forward["query_map"],
                forward["subject_map"],
                orthogroups=orthogroups,
            )
            handle = StringIO()
            converted.loc[:, list(COMPARISON_COLUMNS)].to_csv(
                handle,
                sep=chr(9),
                header=False,
                index=False,
                lineterminator=chr(10),
            )
            converted_pairs.append(
                {
                    "pair_index": display_pair_indices[(query_index, subject_index)],
                    "tsv": handle.getvalue(),
                    "rows": _dataframe_json_rows(converted),
                    "hit_count": int(converted.shape[0]),
                }
            )
        canonical_content = encode_canonical_typed_resource("orthogroupResult", orthogroups)
        return _finalize_losatp_payload({"pairs": converted_pairs},
            canonical_content=canonical_content, canonical_kind="orthogroup-result")
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='convertLosatpPairsToGenomicPayload', stage="helper")})

def get_record_length(path, fmt, record_id=None, record_index=None):
    """Return record length for a GenBank/FASTA file."""
    from Bio import SeqIO
    try:
        fmt_map = {"genbank": "genbank", "fasta": "fasta"}
        if fmt not in fmt_map:
            return json.dumps({'error': serialize_web_error(ValueError(f'Unsupported format: {fmt}'), operation='unknown', stage="helper")})
        records = list(SeqIO.parse(path, fmt_map[fmt]))
        if not records:
            return json.dumps({'error': serialize_web_error(ValueError('No records found'), operation='unknown', stage="helper")})
        if record_id:
            for idx, record in enumerate(records):
                if record.id == record_id:
                    return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": idx})
            return json.dumps({'error': serialize_web_error(ValueError(f'Record ID not found: {record_id}'), operation='unknown', stage="helper")})
        if record_index is not None:
            idx = int(record_index)
            if idx < 0 or idx >= len(records):
                return json.dumps({'error': serialize_web_error(ValueError(f'Record index out of range: {idx + 1}'), operation='unknown', stage="helper")})
            record = records[idx]
            return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": idx})
        record = records[0]
        return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": 0})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='unknown', stage="helper")})

def list_sequence_records(path, format):
    """List record selectors, IDs, and lengths from a sequence file."""
    from Bio import SeqIO
    from gbdraw.api.record_planning import _detected_topology
    from gbdraw.core.record_metadata import (
        format_inferred_definition,
        infer_record_source_metadata,
    )
    try:
        format_map = {"genbank": "genbank", "fasta": "fasta"}
        if format not in format_map:
            return json.dumps({'error': serialize_web_error(ValueError(f'Unsupported format: {format}'), operation='listSequenceRecords', stage="helper")})
        records = list(SeqIO.parse(path, format_map[format]))
        if not records:
            return json.dumps({'error': serialize_web_error(ValueError('No records found'), operation='listSequenceRecords', stage="helper")})
        payload = []
        for idx, record in enumerate(records):
            organism = ""
            strain = ""
            inferred_def = ""
            if format == "genbank":
                meta = infer_record_source_metadata(record)
                organism = meta.organism or ""
                strain = meta.strain or ""
                inferred_def = format_inferred_definition(meta)

            payload.append(
                {
                    "selector": f"#{idx + 1}",
                    "record_id": str(record.id or f"Record_{idx + 1}"),
                    "record_length": len(record.seq),
                    "topology": _detected_topology(record, "genbank") if format == "genbank" else "unknown",
                    "organism": organism,
                    "strain": strain,
                    "inferred_definition": inferred_def,
                }
            )
        return json.dumps({"records": payload})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='listSequenceRecords', stage="helper")})

def list_gff_fasta_records(gff_path, fasta_path):
    """List the records load_gff_fasta reads: GFF3 records the FASTA names, in FASTA order."""
    try:
        from gbdraw.io.genome import load_gff_fasta
        records = load_gff_fasta([gff_path], [fasta_path])
        payload = [
            {
                "selector": f"#{idx + 1}",
                "record_id": str(record.id or f"Record_{idx + 1}"),
                "record_length": len(record.seq),
                "topology": "unknown",
            }
            for idx, record in enumerate(records)
        ]
        return json.dumps({"records": payload})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='listGffFastaRecords', stage="helper")})

def measure_legend_text_json(caption, font_family="Arial", font_size=14, config_json="null", overrides_json="{}"):
    """Measure one legend caption at the DPI the renderer resolves from the request config."""
    try:
        from gbdraw.api.config import apply_config_overrides
        from gbdraw.core.text import calculate_bbox_dimensions
        cfg = apply_config_overrides(json.loads(str(config_json)), json.loads(str(overrides_json)))
        width, _ = calculate_bbox_dimensions(
            str(caption),
            str(font_family or "Arial"),
            float(font_size or 14),
            int(cfg.canvas.dpi),
        )
        return json.dumps({"width": width})
    except Exception as error:
        return json.dumps({'error': serialize_web_error(error, operation='measureLegendText', stage="helper")})

def generate_legend_entry_svg(caption, color, y_offset, rect_size=14, font_size=14, font_family="Arial", x_offset=0, stroke_color="black", stroke_width=0.5):
    """Generate SVG elements for a single legend entry"""
    from xml.sax.saxutils import escape as xml_escape

    # Create color rectangle path with proper stroke (matching original legend entries)
    half = rect_size / 2
    rect_d = f"M 0,{-half} L {rect_size},{-half} L {rect_size},{half} L 0,{half} z"
    rect_svg = f'<path d="{rect_d}" fill="{color}" stroke="{stroke_color}" stroke-width="{stroke_width}" transform="translate({x_offset}, {y_offset})"/>'

    # Create text element
    x_margin = (22 / 14) * rect_size
    safe_caption = xml_escape(str(caption))
    text_svg = f'<text font-size="{font_size}" font-family="{font_family}" dominant-baseline="central" text-anchor="start" transform="translate({x_offset + x_margin}, {y_offset})">{safe_caption}</text>'

    return json.dumps({"rect": rect_svg, "text": text_svg})

def extract_features_from_genbank(gb_path, region_spec=None, record_selector=None, reverse_flag=None, selected_features=None, feature_visibility_table_path=None, include_biological_features=False):
    """Extract feature info from GenBank file for UI display."""
    return extract_features_from_genbank_json(
        gb_path,
        region_spec=region_spec,
        record_selector=record_selector,
        reverse_flag=reverse_flag,
        selected_features=selected_features,
        feature_visibility_table_path=feature_visibility_table_path,
        include_biological_features=include_biological_features,
    )

def extract_features_from_gff_fasta(gff_path, fasta_path, region_spec=None, record_selector=None, reverse_flag=None, selected_features=None, feature_visibility_table_path=None, include_biological_features=False):
    """Extract feature info from paired GFF3 and FASTA files for UI display."""
    return extract_features_from_gff_fasta_json(
        gff_path,
        fasta_path,
        region_spec=region_spec,
        record_selector=record_selector,
        reverse_flag=reverse_flag,
        selected_features=selected_features,
        feature_visibility_table_path=feature_visibility_table_path,
        include_biological_features=include_biological_features,
    )


_WEB_JSON_HELPERS = {
    "resolve_similarity_alignment_json": (resolve_similarity_alignment_json, "resolveSimilarityAlignment"),
    "evaluate_rules_json": (evaluate_rules_json, "evaluateRules"),
    "read_pdf_font": (read_pdf_font, "readPdfFont"),
    "validate_web_config_overrides_json": (validate_web_config_overrides_json, "validateConfigOverrides"),
    "extract_first_fasta": (extract_first_fasta, "extractFirstFasta"),
    "extract_cds_protein_fasta": (extract_cds_protein_fasta, "extractCdsProteinFasta"),
    "build_protein_losat_cache_keys_json": (build_protein_losat_cache_keys_json, "buildProteinLosatCacheKeys"),
    "promote_legacy_losatp_cache_candidates": (promote_legacy_losatp_cache_candidates, "promoteLegacyLosatpCache"),
    "resolve_legacy_protein_reference_map_json": (resolve_legacy_protein_reference_map_json, "resolveLegacyProteinReferences"),
    "convert_losatp_blastp_pairs_to_genomic_payload": (convert_losatp_blastp_pairs_to_genomic_payload, "convertLosatpPairsToGenomicPayload"),
    "main_session_table_text_to_search_frame_json": (main_session_table_text_to_search_frame_json, "convertMainSessionComparisonFrame"),
    "hydrate_protein_losat_tsv_json": (hydrate_protein_losat_tsv_json, "hydrateProteinLosatTsv"),
    "list_sequence_records": (list_sequence_records, "listSequenceRecords"),
    "list_gff_fasta_records": (list_gff_fasta_records, "listGffFastaRecords"),
    "read_comparison_sequence_json": (read_comparison_sequence_json, "readComparisonSequence"),
    "measure_legend_text_json": (measure_legend_text_json, "measureLegendText"),
    "generate_legend_entry_svg": (generate_legend_entry_svg, "generateLegendEntrySvg"),
    "extract_features_from_genbank": (extract_features_from_genbank, "feature-extraction"),
    "extract_features_from_gff_fasta": (extract_features_from_gff_fasta, "feature-extraction"),
}

@private_web_execution()
def call_web_json_helper(helper_name, *args):
    operation = "unknown"
    try:
        helper, operation = _WEB_JSON_HELPERS[helper_name]
        return helper(*args)
    except Exception as error:
        return json.dumps({"error": serialize_web_error(error, operation=operation, stage="helper")})
`;
