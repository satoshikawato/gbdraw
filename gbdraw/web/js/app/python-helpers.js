export const PYTHON_HELPERS = `
import warnings
warnings.simplefilter('ignore', SyntaxWarning)
try:
    import tomllib
except ImportError:
    import tomli as tomllib
import json
import traceback
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
from gbdraw.web_support.protein_analysis import (
    assemble_web_protein_analysis,
    plan_web_protein_analysis,
)

_WEB_PREPARED_INPUT_CACHE = PreparedBiologicalInputCache()

def _web_structured_error(error, fallback):
    payload = {
        "type": error.__class__.__name__,
        "message": str(error) if str(error) else fallback,
        "traceback": traceback.format_exc(),
    }
    if getattr(error, "status", None) is not None:
        payload["status"] = error.status
    return {"error": payload}

def plan_web_protein_analysis_json(
    intent_json,
    payload_path,
    identity_manifest_json,
    tool_identity,
):
    try:
        with open(str(payload_path), "r", encoding="utf-8") as payload_handle:
            payload = json.load(payload_handle)
        records, extraction, record_payloads = _web_protein_analysis_inputs(
            payload,
            json.loads(str(identity_manifest_json)),
        )
        return json.dumps(plan_web_protein_analysis(
            json.loads(str(intent_json)),
            records,
            extraction,
            record_payloads,
            str(tool_identity),
        ))
    except Exception as error:
        return json.dumps(_web_structured_error(error, "Protein analysis planning failed"))

def assemble_web_protein_analysis_json(
    intent_json,
    payload_path,
    identity_manifest_json,
    tool_identity,
):
    try:
        with open(str(payload_path), "r", encoding="utf-8") as payload_handle:
            payload = json.load(payload_handle)
        records, extraction, record_payloads = _web_protein_analysis_inputs(
            payload,
            json.loads(str(identity_manifest_json)),
        )
        return json.dumps(assemble_web_protein_analysis(
            json.loads(str(intent_json)),
            records,
            extraction,
            record_payloads,
            payload.get("rawEntries"),
            str(tool_identity),
            payload.get("cachedDerivedEntries") or (),
        ))
    except Exception as error:
        return json.dumps(_web_structured_error(error, "Protein analysis failed"))

def _is_blank_or_js_nullish(value):
    if value is None:
        return True
    if type(value).__name__ in {"JsNull", "JsUndefined"}:
        return True
    try:
        return str(value).strip().lower() in {"", "null", "undefined", "none"}
    except Exception:
        return False

def run_canonical_request_wrapper(
    request_json,
    resource_paths_json,
    workspace,
    diagnostics_enabled=False,
    resource_identities_json=None,
):
    try:
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
        return {
            "error": {
                "type": e.__class__.__name__,
                "message": str(e) if str(e) else "Unhandled exception",
                "traceback": traceback.format_exc(),
            }
        }

def extract_first_fasta(path, fmt, region_spec=None, record_selector=None, reverse_flag=None):
    """Extract the first record as FASTA for LOSAT input."""
    from Bio import SeqIO
    from io import StringIO
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record
    try:
        fmt_map = {"genbank": "genbank", "fasta": "fasta"}
        if fmt not in fmt_map:
            return json.dumps({"error": f"Unsupported format: {fmt}"})
        records = list(SeqIO.parse(path, fmt_map[fmt]))
        if not records:
            return json.dumps({"error": "No records found"})
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
        return json.dumps({"error": "No records found"})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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

def _web_transform_blast_pos(position, view_transform):
    normalized = _normalize_web_view_transform(view_transform)
    pos = int(position)
    if not normalized["reverse"]:
        return pos
    return int(normalized["length"]) + 1 - pos

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

def convert_losat_nucleotide_to_display_tsv(blast_text, query_view_transform=None, subject_view_transform=None):
    """Transform cached raw LOSAT nucleotide outfmt 6 rows into display coordinates."""
    try:
        from io import StringIO
        import pandas as pd
        from gbdraw.io.comparisons import COMPARISON_COLUMNS

        data_lines = [
            line
            for line in str(blast_text or "").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not data_lines:
            return json.dumps({"tsv": "", "rows": []})
        df = pd.read_csv(
            StringIO(chr(10).join(data_lines)),
            sep=chr(9),
            names=COMPARISON_COLUMNS,
        )
        query_transform = _normalize_web_view_transform(query_view_transform)
        subject_transform = _normalize_web_view_transform(subject_view_transform)
        for column in ("qstart", "qend"):
            df[column] = df[column].map(lambda value: _web_transform_blast_pos(value, query_transform))
        for column in ("sstart", "send"):
            df[column] = df[column].map(lambda value: _web_transform_blast_pos(value, subject_transform))
        handle = StringIO()
        df.loc[:, list(COMPARISON_COLUMNS)].to_csv(
            handle,
            sep=chr(9),
            header=False,
            index=False,
            lineterminator=chr(10),
        )
        return json.dumps({"tsv": handle.getvalue(), "rows": _dataframe_json_rows(df)})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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
            return json.dumps({"error": f"No CDS proteins found in {record.id}"})
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
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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
            expected_tool_identity=options.get("toolIdentity"),
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
    except Exception:
        return json.dumps({"status": "error", "error": traceback.format_exc()})

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
    except Exception:
        return json.dumps({"status": "error", "error": traceback.format_exc()})

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
    except Exception:
        return json.dumps({"status": "error", "error": traceback.format_exc()})

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

def _web_protein_analysis_inputs(payload, identity_manifest):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from dataclasses import replace
    from gbdraw.analysis.protein_colinearity import (
        ProteinExtractionResult,
        validate_protein_identity_manifest,
    )

    raw_records = payload.get("records") if isinstance(payload, dict) else None
    if not isinstance(raw_records, list) or not raw_records:
        raise ValueError("Protein analysis payload must contain records.")
    manifest = validate_protein_identity_manifest(identity_manifest)
    records = []
    proteins_by_record = []
    protein_map = {}
    instance_keys = []
    analysis_ids = []
    protein_set_hashes = []
    runtime_binding_hashes = []
    display_binding_hashes = []
    for index, raw in enumerate(raw_records):
        if not isinstance(raw, dict) or int(raw.get("recordIndex", -1)) != index:
            raise ValueError("Protein records must be ordered by recordIndex.")
        source_proteins = _web_ordered_proteins_from_fasta(
            raw.get("proteinMap") or {},
            raw.get("fasta") or "",
        )
        display_map = _build_display_web_cds_protein_map(
            raw.get("proteinMap") or {},
            raw.get("viewTransform") or {},
        )
        display_proteins = [
            replace(display_map[item.protein_id], sequence=item.sequence)
            for item in source_proteins
        ]
        instance_key = str(raw.get("recordInstanceKey") or "")
        binding = manifest.binding_for(instance_key)
        analysis_id = str(binding["recordAnalysisId"])
        analysis = manifest.record_analyses[analysis_id]
        length = int(raw.get("recordLength") or 0)
        records.append(SeqRecord(Seq("N" * length), id=str(raw.get("recordId") or "")))
        proteins_by_record.append(display_proteins)
        for protein in display_proteins:
            if protein.protein_id in protein_map:
                raise ValueError("Protein runtime handles must be unique.")
            protein_map[protein.protein_id] = protein
        instance_keys.append(instance_key)
        analysis_ids.append(analysis_id)
        protein_set_hashes.append(str(analysis["proteinSetHash"]))
        runtime_binding_hashes.append(str(binding["runtimeBindingHash"]))
        display_binding_hashes.append(str(binding["displayBindingHash"]))
    extraction = ProteinExtractionResult(
        proteins_by_record=proteins_by_record,
        protein_map=protein_map,
        identity_manifest=manifest,
        record_instance_keys=tuple(instance_keys),
        protein_set_hashes=tuple(protein_set_hashes),
        record_analysis_ids=tuple(analysis_ids),
        runtime_binding_hashes=tuple(runtime_binding_hashes),
        display_binding_hashes=tuple(display_binding_hashes),
    )
    return tuple(records), extraction, tuple(raw_records)

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

def get_record_length(path, fmt, record_id=None, record_index=None):
    """Return record length for a GenBank/FASTA file."""
    from Bio import SeqIO
    try:
        fmt_map = {"genbank": "genbank", "fasta": "fasta"}
        if fmt not in fmt_map:
            return json.dumps({"error": f"Unsupported format: {fmt}"})
        records = list(SeqIO.parse(path, fmt_map[fmt]))
        if not records:
            return json.dumps({"error": "No records found"})
        if record_id:
            for idx, record in enumerate(records):
                if record.id == record_id:
                    return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": idx})
            return json.dumps({"error": f"Record ID not found: {record_id}"})
        if record_index is not None:
            idx = int(record_index)
            if idx < 0 or idx >= len(records):
                return json.dumps({"error": f"Record index out of range: {idx + 1}"})
            record = records[idx]
            return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": idx})
        record = records[0]
        return json.dumps({"length": len(record.seq), "record_id": record.id, "record_index": 0})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def list_sequence_records(path, format):
    """List record selectors, IDs, and lengths from a sequence file."""
    from Bio import SeqIO
    try:
        format_map = {"genbank": "genbank", "fasta": "fasta"}
        if format not in format_map:
            return json.dumps({"error": f"Unsupported format: {format}"})
        records = list(SeqIO.parse(path, format_map[format]))
        if not records:
            return json.dumps({"error": "No records found"})
        payload = []
        for idx, record in enumerate(records):
            payload.append(
                {
                    "selector": f"#{idx + 1}",
                    "record_id": str(record.id or f"Record_{idx + 1}"),
                    "record_length": len(record.seq),
                }
            )
        return json.dumps({"records": payload})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def list_gff_fasta_records(gff_path, fasta_path):
    """List every FASTA record available to paired GFF3 diagram generation."""
    try:
        from Bio import SeqIO
        records = list(SeqIO.parse(fasta_path, "fasta"))
        payload = [
            {
                "selector": f"#{idx + 1}",
                "record_id": str(record.id or f"Record_{idx + 1}"),
                "record_length": len(record.seq),
            }
            for idx, record in enumerate(records)
        ]
        return json.dumps({"records": payload})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def measure_legend_text_json(caption, font_family="Arial", font_size=14):
    """Measure one legend caption with the packaged gbdraw font metrics."""
    try:
        from gbdraw.core.text import calculate_bbox_dimensions
        width, _ = calculate_bbox_dimensions(
            str(caption),
            str(font_family or "Arial"),
            float(font_size or 14),
            72,
        )
        return json.dumps({"width": width})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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

def regenerate_definition_svgs(
    gb_path,
    species=None,
    strain=None,
    plot_title=None,
    font_size=None,
    plot_title_font_size=None,
    plot_title_position="none",
    multi_record_canvas=False,
    keep_full_definition_with_plot_title=False,
):
    """Regenerate definition group SVGs for all records in an input file"""
    from Bio import SeqIO
    from gbdraw.render.groups.circular.definition import DefinitionGroup
    from gbdraw.canvas import CircularCanvasConfigurator
    from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
    from gbdraw.svg.ids import definition_group_svg_id
    from importlib import resources

    try:
        # Load default config
        with resources.files("gbdraw.data").joinpath("config.toml").open("rb") as fh:
            config_dict = tomllib.load(fh)

        # Override font sizes if provided
        if not _is_blank_or_js_nullish(font_size):
            config_dict["objects"]["definition"]["circular"]["font_size"] = float(font_size)
        if not _is_blank_or_js_nullish(plot_title_font_size):
            config_dict["objects"]["definition"]["circular"]["plot_title_font_size"] = float(plot_title_font_size)
        cfg = GbdrawConfig.from_dict(config_dict)
        render_profile = CircularRenderProfile(cfg)

        # Parse the GenBank file
        records = list(SeqIO.parse(gb_path, "genbank"))
        if not records:
            return json.dumps({"error": "No records found"})

        normalized_plot_title_position = str(plot_title_position or "none").strip().lower()
        if normalized_plot_title_position not in {"none", "top", "bottom"}:
            normalized_plot_title_position = "none"
        normalized_plot_title = str(plot_title or "").strip()
        show_plot_title = normalized_plot_title_position in {"top", "bottom"}
        keep_full_definition = bool(keep_full_definition_with_plot_title)

        definitions = []
        record_count = len(records)
        record_id_counts = {}
        for record in records:
            raw_record_id = str(record.id)
            record_id_counts[raw_record_id] = record_id_counts.get(raw_record_id, 0) + 1
        for index, record in enumerate(records):
            # Create canvas config
            canvas_config = CircularCanvasConfigurator(
                output_prefix=f"temp_{index}",
                profile=render_profile,
                legend="none",
                gb_record=record,
            )

            if show_plot_title and keep_full_definition:
                profile = "full"
            else:
                profile = "record_summary" if bool(multi_record_canvas) or show_plot_title else "full"
            raw_record_id = str(record.id)
            has_duplicate_record_id = record_id_counts[raw_record_id] > 1
            definition_group_id = (
                definition_group_svg_id(
                    raw_record_id,
                    mode="circular",
                    record_index=index,
                    record_count=record_count,
                )
                if has_duplicate_record_id
                else None
            )
            def_group = DefinitionGroup(
                gb_record=record,
                canvas_config=canvas_config,
                species=species if species else None,
                strain=strain if strain else None,
                plot_title=None,
                definition_profile=profile,
                definition_group_id=definition_group_id,
                record_index=index,
                record_count=record_count if has_duplicate_record_id else 1,
                cfg=cfg,
            )

            group = def_group.get_group()
            definitions.append(
                {
                    "svg": group.tostring(),
                    "definition_group_id": def_group.definition_group_id,
                    "record_index": index,
                }
            )

        if show_plot_title:
            shared_canvas_config = CircularCanvasConfigurator(
                output_prefix="temp_shared",
                profile=render_profile,
                legend="none",
                gb_record=records[0],
            )
            shared_group = DefinitionGroup(
                gb_record=records[0],
                canvas_config=shared_canvas_config,
                species=species if species else None,
                strain=strain if strain else None,
                plot_title=normalized_plot_title if normalized_plot_title else None,
                definition_profile="shared_common",
                definition_group_id="plot_title",
                cfg=cfg,
            )
            definitions.append(
                {
                    "svg": shared_group.get_group().tostring(),
                    "definition_group_id": "plot_title",
                    "record_index": None,
                }
            )

        return json.dumps({"definitions": definitions})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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

`;
