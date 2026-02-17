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

def run_gbdraw_wrapper(mode, args):
    for f in glob.glob("*.svg"):
        try:
            os.remove(f)
        except:
            pass

    full_args = args + ["-f", "svg"]
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()

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

    original_streams = []
    try:
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

def regenerate_definition_svg(gb_path, species=None, strain=None, font_size=18):
    """Regenerate just the definition group SVG for instant preview"""
    from Bio import SeqIO
    from gbdraw.render.groups.circular.definition import DefinitionGroup
    from gbdraw.canvas import CircularCanvasConfigurator
    from importlib import resources

    try:
        # Load default config
        with resources.files("gbdraw.data").joinpath("config.toml").open("rb") as fh:
            config_dict = tomllib.load(fh)

        # Override font size if provided
        if font_size:
            config_dict["objects"]["definition"]["circular"]["font_size"] = font_size

        # Parse the GenBank file
        records = list(SeqIO.parse(gb_path, "genbank"))
        if not records:
            return json.dumps({"error": "No records found"})

        record = records[0]

        # Create canvas config
        canvas_config = CircularCanvasConfigurator(
            output_prefix="temp",
            config_dict=config_dict,
            legend="none",
            gb_record=record,
        )

        # Create definition group with custom species/strain
        def_group = DefinitionGroup(
            gb_record=record,
            canvas_config=canvas_config,
            config_dict=config_dict,
            species=species if species else None,
            strain=strain if strain else None,
        )

        # Get the SVG content
        group = def_group.get_group()
        # Serialize to string using svgwrite's tostring() method
        svg_content = group.tostring()

        # Return the definition group ID (with _definition suffix)
        definition_group_id = def_group.definition_group_id
        return json.dumps({"svg": svg_content, "definition_group_id": definition_group_id})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def extract_features_from_genbank(gb_path, region_spec=None, record_selector=None, reverse_flag=None, selected_features=None):
    """Extract feature info from GenBank file for UI display"""
    from Bio import SeqIO
    from gbdraw.features.colors import compute_feature_hash
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record
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
                    svg_id = compute_feature_hash(feat, record_id=hash_record_id)
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
                    import hashlib
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
                    "strand": "+" if strand_raw == 1 else "-",
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
