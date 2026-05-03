export const PYTHON_HELPERS = `
import warnings
warnings.simplefilter('ignore', SyntaxWarning)
import tomllib
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
        return json.dumps({"fasta": handle.getvalue(), "record_id": record.id})
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
    }

def extract_cds_protein_fasta(path, fmt, fasta_path=None, region_spec=None, record_selector=None, reverse_flag=None, record_index=None):
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
        protein_map = {
            protein.protein_id: _serialize_cds_protein(protein)
            for protein in proteins
        }
        return json.dumps({
            "fasta": proteins_to_fasta(proteins),
            "record_id": record.id,
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
        protein_map[str(protein_id)] = CdsProtein(
            protein_id=str(data.get("protein_id") or protein_id),
            record_index=int(data.get("record_index") or 0),
            feature_index=int(data.get("feature_index") or 0),
            record_id=str(data.get("record_id") or ""),
            start=int(data.get("start") or 0),
            end=int(data.get("end") or 0),
            strand=strand,
            label=str(data.get("label") or protein_id),
            protein_length=int(data.get("protein_length") or 0),
            sequence="",
            source_protein_id=data.get("source_protein_id"),
            feature_svg_id=data.get("feature_svg_id"),
        )
    return protein_map

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

def convert_losatp_blastp_pairs_to_genomic_payload(pairs_json, max_hits=5):
    """Convert all LOSATP blastp pair outputs and attach global orthogroup metadata."""
    try:
        from io import StringIO
        from gbdraw.analysis.protein_colinearity import (
            build_orthogroups_from_protein_hits,
            cap_hits_per_query,
            convert_pair_protein_hits_to_genomic_links,
            parse_losatp_outfmt6,
        )
        from gbdraw.io.comparisons import COMPARISON_COLUMNS

        raw_pairs = json.loads(str(pairs_json or "[]"))
        if not isinstance(raw_pairs, list):
            raw_pairs = []

        capped_hits_by_pair = []
        protein_maps_by_pair = []
        combined_protein_map = {}
        pair_indices = []
        for idx, item in enumerate(raw_pairs):
            if not isinstance(item, dict):
                continue
            query_map = _build_web_cds_protein_map(item.get("queryProteinMap") or {})
            subject_map = _build_web_cds_protein_map(item.get("subjectProteinMap") or {})
            hits = parse_losatp_outfmt6(str(item.get("blastText") or ""))
            capped = cap_hits_per_query(hits, max_hits=int(max_hits or 5))
            capped_hits_by_pair.append(capped)
            protein_maps_by_pair.append((query_map, subject_map))
            combined_protein_map.update(query_map)
            combined_protein_map.update(subject_map)
            pair_indices.append(int(item.get("pairIndex", idx)))

        orthogroups = build_orthogroups_from_protein_hits(
            capped_hits_by_pair,
            combined_protein_map,
        )

        converted_pairs = []
        for idx, capped in enumerate(capped_hits_by_pair):
            query_map, subject_map = protein_maps_by_pair[idx]
            converted = convert_pair_protein_hits_to_genomic_links(
                capped,
                query_map,
                subject_map,
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
                    "pair_index": pair_indices[idx],
                    "tsv": handle.getvalue(),
                    "rows": _dataframe_json_rows(converted),
                    "hit_count": int(converted.shape[0]),
                }
            )

        orthogroup_payload = []
        for orthogroup_id, members in orthogroups.orthogroups.items():
            orthogroup_payload.append(
                {
                    "id": orthogroup_id,
                    "member_count": len(members),
                    "members": [
                        {
                            "orthogroupId": member.orthogroup_id,
                            "proteinId": member.protein_id,
                            "sourceProteinId": member.source_protein_id,
                            "recordIndex": member.record_index,
                            "recordId": member.record_id,
                            "featureIndex": member.feature_index,
                            "featureSvgId": member.feature_svg_id,
                            "start": member.start,
                            "end": member.end,
                            "strand": member.strand,
                            "representative": member.representative,
                        }
                        for member in members
                    ],
                }
            )

        return json.dumps({"pairs": converted_pairs, "orthogroups": orthogroup_payload})
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
    """Extract feature info from GenBank file for UI display"""
    from Bio import SeqIO
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record
    import hashlib

    def _compute_svg_feature_hash(feature, record_id=None):
        loc = feature.location
        if hasattr(loc, "parts") and loc.parts:
            first_part = loc.parts[0]
            start = int(first_part.start)
            end = int(first_part.end)
            strand = first_part.strand
        else:
            start = int(loc.start)
            end = int(loc.end)
            strand = loc.strand
        if record_id is not None:
            key = f"{record_id}:{feature.type}:{start}:{end}:{strand}"
        else:
            key = f"{feature.type}:{start}:{end}:{strand}"
        return "f" + hashlib.md5(key.encode()).hexdigest()[:8]

    features = []
    record_ids = []
    idx = 0
    try:
        records = list(SeqIO.parse(gb_path, "genbank"))
        selector_raw = None
        if record_selector is not None:
            selector_raw = str(record_selector).strip()
            if not selector_raw or selector_raw.lower() in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
                selector_raw = None
        selector = parse_record_selector(selector_raw)
        records = select_record(records, selector)
        reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
        records = reverse_records(records, reverse)
        selected_feature_set = None
        if selected_features is not None:
            parsed_features = []
            if isinstance(selected_features, (list, tuple, set)):
                parsed_features = [str(v).strip() for v in selected_features if str(v).strip()]
            else:
                selected_raw = str(selected_features).strip()
                if selected_raw and selected_raw.lower() not in {"none", "null", "jsnull", "undefined", "jsundefined", "-"}:
                    if selected_raw.startswith("["):
                        try:
                            loaded = json.loads(selected_raw)
                        except Exception:
                            loaded = None
                        if isinstance(loaded, list):
                            parsed_features = [str(v).strip() for v in loaded if str(v).strip()]
                    if not parsed_features:
                        parsed_features = [part.strip() for part in selected_raw.split(",") if part.strip()]
            # Treat null/empty selection as "no feature-type filter" so UI can list all features.
            if parsed_features:
                selected_feature_set = set(parsed_features)
        if region_spec:
            from gbdraw.io.regions import apply_region_specs, parse_region_specs
            records = apply_region_specs(records, parse_region_specs([region_spec]))
        for rec_idx, record in enumerate(records):
            record_id = record.id or f"Record_{rec_idx}"
            hash_record_id = record.id
            record_ids.append(record_id)
            for feat in record.features:
                if selected_feature_set is not None and feat.type not in selected_feature_set:
                    continue

                # Overall coordinates for display
                start = int(feat.location.start)
                end = int(feat.location.end)
                strand_raw = feat.location.strand

                try:
                    svg_id = _compute_svg_feature_hash(feat, record_id=hash_record_id)
                except Exception:
                    svg_id = None
                if not svg_id:
                    # Fallback to raw location if coordinate conversion fails
                    if hasattr(feat.location, 'parts') and feat.location.parts:
                        first_part = feat.location.parts[0]
                        hash_start = int(first_part.start)
                        hash_end = int(first_part.end)
                        hash_strand = first_part.strand
                    else:
                        hash_start = start
                        hash_end = end
                        hash_strand = strand_raw
                    key = f"{feat.type}:{hash_start}:{hash_end}:{hash_strand}"
                    svg_id = "f" + hashlib.md5(key.encode()).hexdigest()[:8]
                qualifiers = {}
                for q_key, q_vals in feat.qualifiers.items():
                    if not q_vals:
                        continue
                    try:
                        q_list = [str(v) for v in q_vals]
                    except Exception:
                        q_list = [str(q_vals)]
                    qualifiers[q_key.lower()] = q_list
                features.append({
                    "id": f"f{idx}",  # Unique internal ID for UI tracking
                    "svg_id": svg_id,  # Matches SVG path id attribute
                    "record_id": record_id,  # For multi-record filtering
                    "record_idx": rec_idx,
                    "type": feat.type,
                    "start": start,
                    "end": end,
                    "strand": "+" if strand_raw == 1 else ("-" if strand_raw == -1 else "undefined"),
                    "locus_tag": feat.qualifiers.get("locus_tag", [""])[0],
                    "gene": feat.qualifiers.get("gene", [""])[0],
                    "product": feat.qualifiers.get("product", [""])[0],
                    "note": feat.qualifiers.get("note", [""])[0][:50] if feat.qualifiers.get("note") else "",
                    "qualifiers": qualifiers,
                })
                idx += 1
    except Exception as e:
        return json.dumps({"error": str(e)})
    return json.dumps({"features": features, "record_ids": record_ids})
`;
