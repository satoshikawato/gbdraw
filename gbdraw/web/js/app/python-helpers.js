export const PYTHON_HELPERS = `
import warnings
warnings.simplefilter('ignore', SyntaxWarning)
try:
    import tomllib
except ImportError:
    import tomli as tomllib
from importlib import resources
import json
import traceback
import glob
import os
import io
import contextlib
import logging
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main
from gbdraw.web_support.feature_metadata import extract_features_from_genbank_json

_WEB_LOSATP_FILTERED_HIT_CACHE = {}
_WEB_LOSATP_CONVERTED_PAYLOAD_CACHE = {}
_WEB_LOSATP_CACHE_ORDER = []
_WEB_LOSATP_CACHE_LIMIT = 64

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

def _web_losatp_json_with_cache_stats(payload_json, **stats):
    try:
        payload = json.loads(str(payload_json))
    except Exception:
        return payload_json
    if not isinstance(payload, dict):
        return payload_json
    cache_payload = payload.get("cache")
    if not isinstance(cache_payload, dict):
        cache_payload = {}
    cache_payload.update(stats)
    payload["cache"] = cache_payload
    return json.dumps(payload)

def _is_blank_or_js_nullish(value):
    if value is None:
        return True
    if type(value).__name__ in {"JsNull", "JsUndefined"}:
        return True
    try:
        return str(value).strip().lower() in {"", "null", "undefined", "none"}
    except Exception:
        return False

def get_palettes_json():
    try:
        with resources.files("gbdraw.data").joinpath("color_palettes.toml").open("rb") as fh:
            return json.dumps(tomllib.load(fh))
    except:
        return "{}"

def run_gbdraw_wrapper(mode, args, virtual_blast_files_json=None):
    for f in glob.glob("*.svg"):
        try:
            os.remove(f)
        except:
            pass

    full_args = args + ["-f", "svg"]
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    original_load_comparisons = None
    assemble_module = None

    def _collect_output():
        stdout_text = stdout_buf.getvalue()
        stderr_text = stderr_buf.getvalue()
        stdout_text = stdout_text.strip() if stdout_text else ""
        stderr_text = stderr_text.strip() if stderr_text else ""
        return stdout_text, stderr_text

    def _build_error(err_type, message, traceback_text=None, code=None):
        stdout_text, stderr_text = _collect_output()
        payload = {"type": err_type, "message": message}
        if code is not None:
            payload["code"] = code
        if traceback_text:
            payload["traceback"] = traceback_text
        if stderr_text:
            payload["stderr"] = stderr_text
        if stdout_text:
            payload["stdout"] = stdout_text
        return {"error": payload}

    def _install_virtual_blast_loader():
        nonlocal original_load_comparisons, assemble_module
        if not virtual_blast_files_json:
            return
        try:
            payload = json.loads(str(virtual_blast_files_json))
        except Exception:
            payload = []
        virtual_files = {
            str(item.get("path", "")): item
            for item in payload
            if isinstance(item, dict) and str(item.get("path", ""))
        }
        if not virtual_files:
            return

        import pandas as pd
        from io import StringIO
        from gbdraw.diagrams.linear import assemble as _assemble_module
        from gbdraw.io.comparisons import COMPARISON_COLUMNS, filter_comparison_dataframe

        assemble_module = _assemble_module
        original_load_comparisons = _assemble_module.load_comparisons

        def _load_comparisons_from_virtual_files(comparison_files, blast_config):
            comparison_list = []
            fallback_files = []
            for comparison_file in comparison_files:
                comparison_path = str(comparison_file)
                if comparison_path not in virtual_files:
                    fallback_files.append(comparison_file)
                    continue
                try:
                    virtual_entry = virtual_files[comparison_path]
                    comparison_rows = virtual_entry.get("rows") if isinstance(virtual_entry, dict) else None
                    comparison_text = str(virtual_entry.get("text", "")) if isinstance(virtual_entry, dict) else str(virtual_entry)
                    if isinstance(comparison_rows, list):
                        if comparison_rows:
                            df = pd.DataFrame.from_records(comparison_rows)
                        else:
                            df = pd.DataFrame(columns=COMPARISON_COLUMNS)
                    elif not comparison_text.strip():
                        comparison_list.append(pd.DataFrame(columns=COMPARISON_COLUMNS))
                        continue
                    else:
                        df = pd.read_csv(
                            StringIO(comparison_text),
                            sep=chr(9),
                            comment="#",
                            names=COMPARISON_COLUMNS,
                        )
                    comparison_list.append(filter_comparison_dataframe(df, blast_config))
                except ValueError as e:
                    logging.getLogger(__name__).warning(
                        f"WARNING: Error parsing comparison file {comparison_path}. It may be corrupt or in the wrong format. Error: {e}"
                    )
                except Exception as e:
                    logging.getLogger(__name__).error(
                        f"ERROR: An unexpected error occurred while processing {comparison_path}: {e}"
                    )
            if fallback_files:
                comparison_list.extend(original_load_comparisons(fallback_files, blast_config))
            return comparison_list

        _assemble_module.load_comparisons = _load_comparisons_from_virtual_files

    original_streams = []
    try:
        _install_virtual_blast_loader()
        for handler in logging.getLogger().handlers:
            if isinstance(handler, logging.StreamHandler):
                original_streams.append((handler, handler.stream))
                handler.setStream(stdout_buf)
        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            if mode == 'circular':
                circular_main(full_args)
            else:
                linear_main(full_args)
    except SystemExit as e:
        code = getattr(e, "code", None)
        if code != 0:
            if isinstance(code, int):
                message = f"Exit with status {code}"
            elif code is None:
                message = "SystemExit"
            else:
                message = str(code)
            return json.dumps(_build_error("SystemExit", message, code=code))
    except Exception as e:
        err_type = e.__class__.__name__
        message = str(e) if str(e) else "Unhandled exception"
        return json.dumps(_build_error(err_type, message, traceback_text=traceback.format_exc()))
    finally:
        for handler, stream in original_streams:
            handler.setStream(stream)
        if assemble_module is not None and original_load_comparisons is not None:
            assemble_module.load_comparisons = original_load_comparisons

    files = glob.glob("*.svg")
    if not files:
        return json.dumps(_build_error("OutputError", "No output files generated."))
    results = []
    for fname in sorted(files):
        with open(fname, "r") as f:
            results.append({"name": fname, "content": f.read()})
    return json.dumps(results)

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

def _compute_web_feature_svg_id(record_id, feature_type, start, end, strand):
    import hashlib

    normalized_record_id = str(record_id or "")
    normalized_type = str(feature_type or "CDS")
    if normalized_record_id:
        key = f"{normalized_record_id}:{normalized_type}:{int(start)}:{int(end)}:{strand}"
    else:
        key = f"{normalized_type}:{int(start)}:{int(end)}:{strand}"
    return "f" + hashlib.md5(key.encode()).hexdigest()[:8]

def _display_feature_svg_id_from_data(data, display_start, display_end, display_strand, view_transform):
    normalized = _normalize_web_view_transform(view_transform)
    if not normalized["reverse"]:
        existing = data.get("feature_svg_id")
        if existing:
            return existing
    hash_start = data.get("feature_hash_start")
    hash_end = data.get("feature_hash_end")
    hash_strand = data.get("feature_hash_strand")
    if hash_start is None or hash_end is None:
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
            "linear",
            selected_features_set=["CDS"],
            keep_all_features=True,
            load_comparison=True,
            record_selectors=[_normalize_web_record_selector(record_selector) or ""],
            reverse_flags=[reverse],
        )
    else:
        raise ValueError(f"Unsupported format: {fmt}")
    if region_spec:
        records = apply_region_specs(records, parse_region_specs([region_spec]))
    if not records:
        raise ValueError("No records found")
    return records[0]

def _serialize_cds_protein(protein):
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
    }

def _safe_web_protein_id_token(value):
    import re

    text = str(value or "").strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = text.strip("._-")
    return text or "record"

def _with_stable_web_protein_ids(proteins, record_instance_key):
    from dataclasses import replace
    import hashlib

    record_key = _safe_web_protein_id_token(record_instance_key)
    remapped = []
    used = set()
    for protein in proteins:
        strand = protein.strand if protein.strand in {-1, 1} else 0
        aa_hash = hashlib.sha256(str(protein.sequence or "").encode()).hexdigest()[:12]
        base_id = f"p_{record_key}_{int(protein.start)}_{int(protein.end)}_{strand}_{aa_hash}"
        protein_id = base_id
        suffix = 2
        while protein_id in used:
            protein_id = f"{base_id}_{suffix}"
            suffix += 1
        used.add(protein_id)
        remapped.append(replace(protein, protein_id=protein_id))
    return remapped

def extract_cds_protein_fasta(path, fmt, fasta_path=None, region_spec=None, record_selector=None, reverse_flag=None, record_index=None, record_instance_key=None):
    """Extract CDS proteins and coordinate metadata for LOSATP blastp."""
    try:
        from gbdraw.analysis.protein_colinearity import extract_cds_proteins, proteins_to_fasta

        record = _load_single_linear_record_for_proteins(
            path,
            fmt,
            fasta_path=fasta_path,
            region_spec=region_spec,
            record_selector=record_selector,
            reverse_flag=reverse_flag,
        )
        record_index_offset = int(record_index) if record_index is not None else 0
        result = extract_cds_proteins(
            [record],
            record_index_offset=record_index_offset,
            prefer_source_ids=False,
        )
        proteins = result.proteins_by_record[0] if result.proteins_by_record else []
        if not proteins:
            return json.dumps({"error": f"No CDS proteins found in {record.id}"})
        stable_record_key = record_instance_key
        if _is_blank_or_js_nullish(stable_record_key):
            stable_record_key = f"r{record_index_offset + 1:04d}_{record.id}"
        proteins = _with_stable_web_protein_ids(proteins, stable_record_key)
        protein_map = {
            protein.protein_id: _serialize_cds_protein(protein)
            for protein in proteins
        }
        return json.dumps({
            "fasta": proteins_to_fasta(proteins),
            "record_id": record.id,
            "record_length": len(record.seq),
            "protein_count": len(proteins),
            "protein_map": protein_map,
        })
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def convert_protein_blast_to_genomic_tsv(blast_text, protein_maps_json, max_hits=5):
    """Convert LOSATP blastp output into genomic comparison TSV for gbdraw linear."""
    try:
        from io import StringIO
        from gbdraw.analysis.protein_colinearity import (
            cap_hits_per_query,
            convert_protein_hits_to_genomic_links,
            parse_losatp_outfmt6,
        )
        from gbdraw.io.comparisons import COMPARISON_COLUMNS
        try:
            from gbdraw.analysis.protein_colinearity import convert_pair_protein_hits_to_genomic_links
        except ImportError:
            convert_pair_protein_hits_to_genomic_links = None

        raw_maps = json.loads(str(protein_maps_json))
        if isinstance(raw_maps, dict):
            raw_maps = [raw_maps]

        protein_maps = [_build_web_cds_protein_map(raw_map) for raw_map in raw_maps]
        hits = parse_losatp_outfmt6(str(blast_text or ""))
        capped = cap_hits_per_query(hits, max_hits=int(max_hits or 5))
        if len(protein_maps) >= 2 and convert_pair_protein_hits_to_genomic_links is not None:
            converted = convert_pair_protein_hits_to_genomic_links(
                capped,
                protein_maps[0],
                protein_maps[1],
            )
        else:
            protein_map = {}
            for current_map in protein_maps:
                protein_map.update(current_map)
            converted = convert_protein_hits_to_genomic_links(capped, protein_map)
        handle = StringIO()
        converted.loc[:, list(COMPARISON_COLUMNS)].to_csv(handle, sep=chr(9), header=False, index=False, lineterminator=chr(10))
        return json.dumps({"tsv": handle.getvalue(), "hit_count": int(converted.shape[0])})
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
            "sequence": "",
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
        for tuple_field in ("db_xref", "parent_ids"):
            if tuple_field in supported_fields:
                raw_values = data.get(tuple_field) or ()
                if isinstance(raw_values, (list, tuple)):
                    kwargs[tuple_field] = tuple(str(value) for value in raw_values if str(value).strip())
                else:
                    kwargs[tuple_field] = (str(raw_values),) if str(raw_values).strip() else ()
        protein_map[str(protein_id)] = CdsProtein(**kwargs)
    return protein_map

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
        display_data["start"] = start
        display_data["end"] = end
        display_data["strand"] = strand
        display_data["feature_svg_id"] = _display_feature_svg_id_from_data(
            display_data,
            start,
            end,
            strand,
            normalized,
        )
        display_map[str(protein_id)] = display_data
    return _build_web_cds_protein_map(display_map)

def _build_web_protein_extraction(protein_maps):
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
    return ProteinExtractionResult(proteins_by_record=proteins_by_record, protein_map=combined)

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

def _serialize_orthogroup_name_candidate(candidate):
    if isinstance(candidate, dict):
        getter = candidate.get
    else:
        getter = lambda key, default=None: getattr(candidate, key, default)
    return {
        "text": str(getter("text", "") or ""),
        "source": str(getter("source", "") or ""),
        "memberCount": int(getter("member_count", 0) or 0),
        "recordCoverageCount": int(getter("record_coverage_count", 0) or 0),
        "representativeCount": int(getter("representative_count", 0) or 0),
        "score": float(getter("score", 0.0) or 0.0),
    }

def _serialize_orthogroups_payload(orthogroups):
    if orthogroups is None:
        return []
    payload = []
    names_by_orthogroup_id = getattr(orthogroups, "names_by_orthogroup_id", {}) or {}
    descriptions_by_orthogroup_id = getattr(orthogroups, "descriptions_by_orthogroup_id", {}) or {}
    candidates_by_orthogroup_id = getattr(orthogroups, "name_candidates_by_orthogroup_id", {}) or {}
    confidence_by_orthogroup_id = getattr(orthogroups, "confidence_by_orthogroup_id", {}) or {}
    for orthogroup_id, members in orthogroups.orthogroups.items():
        record_coverage_count = len({int(member.record_index) for member in members})
        payload.append(
            {
                "id": orthogroup_id,
                "name": str(names_by_orthogroup_id.get(orthogroup_id, "") or ""),
                "description": str(descriptions_by_orthogroup_id.get(orthogroup_id, "") or ""),
                "nameConfidence": str(confidence_by_orthogroup_id.get(orthogroup_id, "none") or "none"),
                "nameCandidates": [
                    _serialize_orthogroup_name_candidate(candidate)
                    for candidate in (candidates_by_orthogroup_id.get(orthogroup_id, []) or [])
                ],
                "member_count": len(members),
                "record_coverage_count": record_coverage_count,
                "members": [
                    {
                        "orthogroupId": member.orthogroup_id,
                        "proteinId": member.protein_id,
                        "sourceProteinId": member.source_protein_id,
                        "recordIndex": member.record_index,
                        "recordId": member.record_id,
                        "featureIndex": member.feature_index,
                        "label": member.label,
                        "featureSvgId": member.feature_svg_id,
                        "start": member.start,
                        "end": member.end,
                        "strand": member.strand,
                        "representative": member.representative,
                        "gene": getattr(member, "gene", None),
                        "product": getattr(member, "product", None),
                        "note": getattr(member, "note", None),
                        "locusTag": getattr(member, "locus_tag", None),
                        "geneId": getattr(member, "gene_id", None),
                        "oldLocusTag": getattr(member, "old_locus_tag", None),
                    }
                    for member in members
                ],
            }
        )
    return payload

def convert_losatp_blastp_pairs_to_genomic_payload(
    pairs_json,
    mode="pairwise",
    max_hits=5,
    bitscore=50,
    evalue="1e-2",
    identity=0,
    alignment_length=0,
    collinear_min_anchors=1,
    collinear_max_gene_gap=0,
    collinear_unit_mode="auto",
    collinear_color_mode="orientation",
    collinear_anchor_mode="rbh",
    collinear_block_merge_gap=50,
    collinear_singleton_merge_gap=25,
    collinear_max_diagonal_drift=0,
    collinear_max_conflicts_in_merge_gap=1,
    collinear_max_paralog_links_per_orthogroup=2,
    collinear_search_scope="adjacent",
):
    """Convert LOSATP blastp outputs for pairwise display or RBH orthogroups."""
    try:
        from io import StringIO
        import pandas as pd
        from gbdraw.analysis.collinearity import (
            LosslessCollinearityParameters,
            build_orthogroup_collinearity_blocks_from_hits,
            convert_collinearity_blocks_to_comparisons,
        )
        from gbdraw.analysis.protein_colinearity import (
            convert_pair_protein_hits_to_genomic_links,
            filter_protein_hits_by_thresholds,
            parse_losatp_outfmt6,
            select_rbh_orthogroup_edges_from_directional_hits,
            select_top_hits_per_query,
        )
        from gbdraw.io.comparisons import COMPARISON_COLUMNS

        raw_payload = json.loads(str(pairs_json or "{}"))
        if not isinstance(raw_payload, dict):
            raise ValueError("LOSATP blastp conversion payload must be an object with 'records' and 'pairs' lists.")
        raw_records = raw_payload.get("records")
        raw_pairs = raw_payload.get("pairs")
        if not isinstance(raw_records, list) or not isinstance(raw_pairs, list):
            raise ValueError("LOSATP blastp conversion payload must contain 'records' and 'pairs' lists.")
        normalized_mode = str(mode or "pairwise").strip().lower()
        if normalized_mode not in {"pairwise", "orthogroup", "collinear"}:
            normalized_mode = "pairwise"

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
            pair_index = int(item.get("pairIndex", min(query_index, subject_index)))
            cache_key = str(item.get("cacheKey") or "").strip()
            if not cache_key:
                raise ValueError(f"LOSATP pair payload #{idx + 1} is missing cacheKey.")
            pair_payloads.append(
                {
                    "pair_index": pair_index,
                    "query_index": query_index,
                    "subject_index": subject_index,
                    "cache_key": cache_key,
                    "blast_text": str(item.get("blastText") or ""),
                }
            )

        conversion_cache_key = (
            normalized_mode,
            int(max_hits or 5),
            str(bitscore),
            str(evalue),
            str(identity),
            str(alignment_length),
            str(collinear_min_anchors),
            str(collinear_max_gene_gap),
            str(collinear_unit_mode),
            str(collinear_color_mode),
            str(collinear_anchor_mode),
            str(collinear_block_merge_gap),
            str(collinear_singleton_merge_gap),
            str(collinear_max_diagonal_drift),
            str(collinear_max_conflicts_in_merge_gap),
            str(collinear_max_paralog_links_per_orthogroup),
            str(collinear_search_scope),
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
                (item["pair_index"], item["query_index"], item["subject_index"], item["cache_key"])
                for item in pair_payloads
            ),
        )
        cached_payload = _web_losatp_cache_get("converted", conversion_cache_key)
        if cached_payload is not None:
            return _web_losatp_json_with_cache_stats(
                cached_payload,
                convertedPayloadHit=True,
                filteredHitCacheHits=0,
                filteredHitCacheMisses=0,
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
                str(bitscore),
                str(evalue),
                str(identity),
                str(alignment_length),
            )
            filtered = _web_losatp_cache_get("filtered", filter_cache_key)
            if filtered is None:
                hits = parse_losatp_outfmt6(item["blast_text"])
                filtered = filter_protein_hits_by_thresholds(
                    hits,
                    evalue=evalue,
                    bitscore=bitscore,
                    identity=identity,
                    alignment_length=alignment_length,
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
        }

        def _finalize_losatp_payload(payload):
            payload["cache"] = cache_stats
            result_json = json.dumps(payload)
            _web_losatp_cache_set("converted", conversion_cache_key, result_json)
            return result_json

        if normalized_mode == "collinear":
            record_ids = []
            for record_proteins in extraction.proteins_by_record:
                if record_proteins:
                    record_ids.append(str(record_proteins[0].record_id))
                else:
                    record_ids.append("")
            def _collinear_int(value, default):
                if _is_blank_or_js_nullish(value):
                    return int(default)
                return int(value)
            def _collinear_float(value, default):
                if _is_blank_or_js_nullish(value):
                    return float(default)
                return float(value)
            def _collinear_text(value, default):
                if _is_blank_or_js_nullish(value):
                    return default
                return str(value).strip()
            anchor_mode = _collinear_text(collinear_anchor_mode, "rbh").lower().replace("-", "_")
            search_scope = _collinear_text(collinear_search_scope, "adjacent").lower().replace("-", "_")
            directional_tables = {}
            for item in pair_items:
                query_index = int(item["query_index"])
                subject_index = int(item["subject_index"])
                if query_index != subject_index:
                    directional_tables[(query_index, subject_index)] = item["hits"]
            params = LosslessCollinearityParameters(
                min_anchors=_collinear_int(collinear_min_anchors, 1),
                max_unit_gap=_collinear_int(collinear_max_gene_gap, 0),
                max_diagonal_drift=_collinear_int(collinear_max_diagonal_drift, 0),
            )
            collinearity_result = build_orthogroup_collinearity_blocks_from_hits(
                directional_tables,
                extraction,
                params=params,
                unit_mode="cds",
                edge_mode=anchor_mode,
                search_scope=search_scope,
            )
            converted_frames = convert_collinearity_blocks_to_comparisons(
                collinearity_result,
                record_ids=record_ids,
                color_mode=_collinear_text(collinear_color_mode, "orientation"),
            )
            converted_pairs = []
            for pair_index, converted in enumerate(converted_frames):
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
            block_payload = [
                {
                    "id": block.block_id,
                    "kind": block.kind,
                    "orientation": block.orientation,
                    "score": block.score,
                    "blockEvalue": block.block_evalue,
                    "anchorCount": len(block.anchors),
                    "anchor_count": len(block.anchors),
                    "orthogroupIds": sorted({anchor.orthogroup_id for anchor in block.anchors if anchor.orthogroup_id}),
                    "querySpan": [
                        min(min(anchor.query_start, anchor.query_end) for anchor in block.anchors),
                        max(max(anchor.query_start, anchor.query_end) for anchor in block.anchors),
                    ],
                    "subjectSpan": [
                        min(min(anchor.subject_start, anchor.subject_end) for anchor in block.anchors),
                        max(max(anchor.subject_start, anchor.subject_end) for anchor in block.anchors),
                    ],
                    "queryRecordIndex": block.query_record_index,
                    "subjectRecordIndex": block.subject_record_index,
                }
                for block in collinearity_result.blocks
            ]
            return _finalize_losatp_payload({
                "pairs": converted_pairs,
                "orthogroups": _serialize_orthogroups_payload(collinearity_result.orthogroups),
                "collinearityBlocks": block_payload,
            })

        if normalized_mode == "pairwise":
            converted_pairs = []
            for item in pair_items:
                display_hits = select_top_hits_per_query(
                    item["hits"],
                    max_hits=int(max_hits or 5),
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
        edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
            directional_tables,
            combined_protein_map,
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
                    "pair_index": min(query_index, subject_index),
                    "tsv": handle.getvalue(),
                    "rows": _dataframe_json_rows(converted),
                    "hit_count": int(converted.shape[0]),
                }
            )

        return _finalize_losatp_payload({"pairs": converted_pairs, "orthogroups": _serialize_orthogroups_payload(orthogroups)})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

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

def list_genbank_records(gb_path):
    """List record selectors and IDs from a GenBank file."""
    from Bio import SeqIO
    try:
        records = list(SeqIO.parse(gb_path, "genbank"))
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
    from importlib import resources

    try:
        # Load default config
        with resources.files("gbdraw.data").joinpath("config.toml").open("rb") as fh:
            config_dict = tomllib.load(fh)

        # Override font sizes if provided
        if font_size is not None:
            config_dict["objects"]["definition"]["circular"]["font_size"] = float(font_size)
        if plot_title_font_size is not None:
            config_dict["objects"]["definition"]["circular"]["plot_title_font_size"] = float(plot_title_font_size)

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
        for index, record in enumerate(records):
            # Create canvas config
            canvas_config = CircularCanvasConfigurator(
                output_prefix=f"temp_{index}",
                config_dict=config_dict,
                legend="none",
                gb_record=record,
            )

            if show_plot_title and keep_full_definition:
                profile = "full"
            else:
                profile = "record_summary" if bool(multi_record_canvas) or show_plot_title else "full"
            def_group = DefinitionGroup(
                gb_record=record,
                canvas_config=canvas_config,
                config_dict=config_dict,
                species=species if species else None,
                strain=strain if strain else None,
                plot_title=None,
                definition_profile=profile,
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
                config_dict=config_dict,
                legend="none",
                gb_record=records[0],
            )
            shared_group = DefinitionGroup(
                gb_record=records[0],
                canvas_config=shared_canvas_config,
                config_dict=config_dict,
                species=species if species else None,
                strain=strain if strain else None,
                plot_title=normalized_plot_title if normalized_plot_title else None,
                definition_profile="shared_common",
                definition_group_id="plot_title",
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

def regenerate_definition_svg(
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
    """Backward-compatible single-record definition regeneration helper"""
    result_json = regenerate_definition_svgs(
        gb_path,
        species=species,
        strain=strain,
        plot_title=plot_title,
        font_size=font_size,
        plot_title_font_size=plot_title_font_size,
        plot_title_position=plot_title_position,
        multi_record_canvas=multi_record_canvas,
        keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
    )
    try:
        payload = json.loads(result_json)
    except Exception:
        return result_json
    if payload.get("error"):
        return result_json
    definitions = payload.get("definitions") or []
    if not definitions:
        return json.dumps({"error": "No definitions generated"})
    first = definitions[0]
    return json.dumps(
        {
            "svg": first.get("svg", ""),
            "definition_group_id": first.get("definition_group_id", ""),
        }
    )

def extract_features_from_genbank(gb_path, region_spec=None, record_selector=None, reverse_flag=None, selected_features=None):
    """Extract feature info from GenBank file for UI display."""
    return extract_features_from_genbank_json(
        gb_path,
        region_spec=region_spec,
        record_selector=record_selector,
        reverse_flag=reverse_flag,
        selected_features=selected_features,
    )

`;
