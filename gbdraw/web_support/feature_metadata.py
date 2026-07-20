from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

from gbdraw.core.record_metadata import (
    _absolute_display_interval,
    _read_coord_map as _read_record_coord_map,
)
from gbdraw.features.selector_values import build_feature_selector_values
from gbdraw.features.ids import (
    compute_feature_hash_from_location_parts,
    make_linear_rendered_feature_id,
)
from gbdraw.features.visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
    should_render_feature,
)


_NULLISH_TEXT = {"", "none", "null", "jsnull", "undefined", "jsundefined", "-"}


def _normalize_record_selector(record_selector: object | None) -> str | None:
    if record_selector is None:
        return None
    selector_raw = str(record_selector).strip()
    if selector_raw.lower() in _NULLISH_TEXT:
        return None
    return selector_raw


def _normalize_optional_path(path: object | None) -> str | None:
    if path is None:
        return None
    normalized = str(path).strip()
    if normalized.lower() in _NULLISH_TEXT:
        return None
    return normalized


def _normalize_selected_feature_set(selected_features: object | None) -> set[str] | None:
    if selected_features is None:
        return None
    parsed_features: list[str] = []
    if isinstance(selected_features, (list, tuple, set)):
        parsed_features = [str(value).strip() for value in selected_features if str(value).strip()]
    else:
        selected_raw = str(selected_features).strip()
        if selected_raw.lower() in _NULLISH_TEXT:
            return None
        if selected_raw.startswith("["):
            try:
                loaded = json.loads(selected_raw)
            except Exception:
                loaded = None
            if isinstance(loaded, list):
                parsed_features = [str(value).strip() for value in loaded if str(value).strip()]
        if not parsed_features:
            parsed_features = [part.strip() for part in selected_raw.split(",") if part.strip()]
    return set(parsed_features) if parsed_features else None


def _normalize_qualifier_values(value: object | None) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        source = value
    else:
        source = [value]
    normalized = []
    for item in source:
        if item is None:
            continue
        normalized.append(str(item))
    return normalized


def _first_qualifier_value(qualifiers: dict[str, object], key: str) -> str:
    values = _normalize_qualifier_values(qualifiers.get(key))
    for value in values:
        if value.strip():
            return value
    return ""


def _get_record_organism(record: Any) -> str:
    annotations = getattr(record, "annotations", None) or {}
    organism = str(annotations.get("organism") or "").strip()
    if organism:
        return organism
    for feature in getattr(record, "features", []) or []:
        if str(getattr(feature, "type", "")).lower() != "source":
            continue
        organism = _first_qualifier_value(getattr(feature, "qualifiers", {}) or {}, "organism").strip()
        if organism:
            return organism
    return ""


def _strand_display(strand: object | None) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "undefined"


def _biological_strand(strand: object | None, coord_step: int) -> int | None:
    """Return the strand in the source-record coordinate system."""

    if strand not in {-1, 1}:
        return None
    return int(strand) * (1 if int(coord_step) >= 0 else -1)


def _get_location_parts(location: Any) -> list[Any]:
    if hasattr(location, "parts") and location.parts:
        return list(location.parts)
    return [location]


def _biological_location_parts(
    location: Any,
    coord_base: int,
    coord_step: int,
) -> list[tuple[int, int, int | None]]:
    """Map processed-record locations back to their source-record coordinates."""

    parts: list[tuple[int, int, int | None]] = []
    for part in _get_location_parts(location):
        try:
            start, end = _absolute_display_interval(
                int(part.start),
                int(part.end),
                coord_base,
                coord_step,
            )
        except Exception:
            continue
        strand = part.strand if part.strand is not None else location.strand
        parts.append((start, end, _biological_strand(strand, coord_step)))
    return parts


def _biological_selector_values(
    feature: Any,
    *,
    record_id: str | None,
    coord_base: int,
    coord_step: int,
) -> tuple[dict[str, object], str, str]:
    """Return selector values plus biological and processed-record feature IDs."""

    selector = build_feature_selector_values(feature, record_id=record_id)
    rendered_feature_id = str(selector.get("hash") or "")
    parts = _biological_location_parts(feature.location, coord_base, coord_step)
    stable_feature_id = compute_feature_hash_from_location_parts(
        str(getattr(feature, "type", "") or ""),
        parts,
        record_id=record_id,
    )
    selector["hash"] = stable_feature_id
    if parts:
        start = min(part[0] for part in parts)
        end = max(part[1] for part in parts)
        strand = _strand_display(
            _biological_strand(feature.location.strand, coord_step)
        )
        selector["location"] = f"{start}..{end}"
        if record_id:
            selector["record_location"] = f"{record_id}:{start}..{end}:{strand}"
        else:
            selector.pop("record_location", None)
    return selector, stable_feature_id, rendered_feature_id


def _iter_features(features: Any):
    """Yield top-level and nested GFF features in source order."""

    for feature in features or []:
        yield feature
        yield from _iter_features(getattr(feature, "sub_features", None))






def _format_location_parts(
    location: Any,
    coord_base: int = 1,
    coord_step: int = 1,
) -> list[dict[str, object]]:
    return [
        {
            "start": start,
            "end": end,
            "strand": _strand_display(strand),
            "display": f"{start + 1}..{end}",
        }
        for start, end, strand in _biological_location_parts(
            location,
            coord_base,
            coord_step,
        )
    ]


def _location_has_fuzzy_positions(location: Any) -> bool:
    for part in _get_location_parts(location):
        if type(part.start).__name__ != "ExactPosition" or type(part.end).__name__ != "ExactPosition":
            return True
    return False


def _extract_nucleotide_sequence(feature: Any, record: Any) -> tuple[str, list[str]]:
    try:
        return str(feature.extract(record.seq)).upper(), []
    except Exception as exc:
        return "", [f"Nucleotide sequence extraction skipped: {exc}"]


def _extract_amino_acid_sequence(feature: Any, nucleotide_sequence: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    qualifiers = feature.qualifiers or {}
    translation = _first_qualifier_value(qualifiers, "translation")
    if translation:
        return "".join(str(translation).split()), warnings

    if str(feature.type).upper() != "CDS":
        return "", warnings

    if "pseudo" in qualifiers or "pseudogene" in qualifiers:
        warnings.append("CDS translation skipped for pseudo/pseudogene feature.")
        return "", warnings

    if _location_has_fuzzy_positions(feature.location):
        warnings.append("CDS translation skipped for fuzzy feature location.")
        return "", warnings

    if not nucleotide_sequence:
        warnings.append("CDS translation skipped because nucleotide sequence is unavailable.")
        return "", warnings

    codon_start_raw = _first_qualifier_value(qualifiers, "codon_start") or "1"
    try:
        codon_start = int(str(codon_start_raw).strip())
    except Exception:
        warnings.append(f"CDS translation skipped because codon_start is invalid: {codon_start_raw}")
        return "", warnings
    if codon_start not in {1, 2, 3}:
        warnings.append(f"CDS translation skipped because codon_start is outside 1..3: {codon_start}")
        return "", warnings

    transl_table_raw = _first_qualifier_value(qualifiers, "transl_table") or "1"
    try:
        transl_table = int(str(transl_table_raw).strip())
    except Exception:
        warnings.append(f"CDS translation skipped because transl_table is invalid: {transl_table_raw}")
        return "", warnings

    coding_sequence = str(nucleotide_sequence)[codon_start - 1 :]
    if len(coding_sequence) == 0:
        warnings.append("CDS translation skipped because coding sequence is empty after codon_start.")
        return "", warnings
    if len(coding_sequence) % 3 != 0:
        warnings.append("CDS translation skipped because coding sequence length is not divisible by 3.")
        return "", warnings

    try:
        protein = str(Seq(coding_sequence).translate(table=transl_table, to_stop=False))
        if protein.endswith("*"):
            protein = protein[:-1]
        return protein, warnings
    except Exception as exc:
        warnings.append(f"CDS translation skipped: {exc}")
        return "", warnings


def _build_selector_safety_scope(records: list[Any]) -> list[dict[str, object]]:
    scope: list[dict[str, object]] = []
    for rec_idx, record in enumerate(records):
        record_id = record.id or f"Record_{rec_idx}"
        hash_record_id = record.id
        coord_base, coord_step = _read_record_coord_map(record)
        for feat in _iter_features(getattr(record, "features", None)):
            selector, _, _ = _biological_selector_values(
                feat,
                record_id=hash_record_id,
                coord_base=coord_base,
                coord_step=coord_step,
            )
            scope.append(
                {
                    "record_id": record_id,
                    "record_idx": rec_idx,
                    "feature_type": str(getattr(feat, "type", "") or ""),
                    "selector": selector,
                }
            )
    return scope


def extract_features_from_records_payload(
    records: Any,
    *,
    selected_features: object | None = None,
    feature_visibility_rules: list[dict[str, Any]] | None = None,
    specific_color_rules: dict | None = None,
    linear_rendered_feature_ids: bool = False,
    include_biological_features: bool = False,
) -> dict[str, object]:
    """Extract feature metadata from processed records.

    ``features`` remains the display-filtered payload used by existing callers.
    When ``include_biological_features`` is true, ``biological_features`` also
    contains every source feature, keyed by its record index and stable feature
    ID, whether or not that feature is eligible for rendering.
    """

    records = list(records or [])
    selected_feature_set = _normalize_selected_feature_set(selected_features)

    features: list[dict[str, object]] = []
    biological_features: list[dict[str, object]] = []
    record_ids: list[str] = []
    selector_safety_scope = _build_selector_safety_scope(records)
    idx = 0
    biological_idx = 0
    for rec_idx, record in enumerate(records):
        record_id = record.id or f"Record_{rec_idx}"
        hash_record_id = record.id
        organism = _get_record_organism(record)
        coord_base, coord_step = _read_record_coord_map(record)
        record_ids.append(record_id)
        for feature_index, feat in enumerate(_iter_features(record.features)):
            is_rendered_feature = should_render_feature(
                feat,
                selected_feature_set,
                feature_visibility_rules=feature_visibility_rules,
                record_id=hash_record_id,
                specific_color_rules=specific_color_rules,
            )
            if not is_rendered_feature and not include_biological_features:
                continue

            feature_start = int(feat.location.start)
            feature_end = int(feat.location.end)
            start, end = _absolute_display_interval(
                feature_start,
                feature_end,
                coord_base,
                coord_step,
            )
            strand_raw = _biological_strand(feat.location.strand, coord_step)
            location_parts = _format_location_parts(
                feat.location,
                coord_base,
                coord_step,
            )
            nucleotide_sequence, sequence_warnings = _extract_nucleotide_sequence(feat, record)
            amino_acid_sequence, translation_warnings = _extract_amino_acid_sequence(
                feat,
                nucleotide_sequence,
            )
            sequence_warnings.extend(translation_warnings)

            selector, stable_svg_id, rendered_stable_svg_id = (
                _biological_selector_values(
                    feat,
                    record_id=hash_record_id,
                    coord_base=coord_base,
                    coord_step=coord_step,
                )
            )

            qualifiers = {}
            for q_key, q_vals in feat.qualifiers.items():
                q_list = _normalize_qualifier_values(q_vals)
                if not q_list:
                    continue
                qualifiers[q_key.lower()] = q_list

            feature_payload = {
                "id": f"f{biological_idx}",
                "svg_id": stable_svg_id,
                "stable_svg_id": stable_svg_id,
                "stable_feature_id": stable_svg_id,
                "record_id": record_id,
                "record_idx": rec_idx,
                "feature_index": feature_index,
                "organism": organism,
                "type": feat.type,
                "start": start,
                "end": end,
                "strand": _strand_display(strand_raw),
                "protein_id": _first_qualifier_value(feat.qualifiers, "protein_id"),
                "source_protein_id": _first_qualifier_value(
                    feat.qualifiers, "protein_id"
                ),
                "locus_tag": _first_qualifier_value(feat.qualifiers, "locus_tag"),
                "gene_id": _first_qualifier_value(feat.qualifiers, "gene_id"),
                "old_locus_tag": _first_qualifier_value(
                    feat.qualifiers, "old_locus_tag"
                ),
                "gene": _first_qualifier_value(feat.qualifiers, "gene"),
                "product": _first_qualifier_value(feat.qualifiers, "product"),
                "note": _first_qualifier_value(feat.qualifiers, "note")[:50],
                "qualifiers": qualifiers,
                "selector": selector,
                "location_parts": location_parts,
                "nucleotide_sequence": nucleotide_sequence,
                "amino_acid_sequence": amino_acid_sequence,
                "sequence_warnings": sequence_warnings,
            }
            if include_biological_features:
                biological_features.append(feature_payload)
                biological_idx += 1

            if not is_rendered_feature:
                continue
            rendered_feature_payload = dict(feature_payload)
            rendered_feature_payload.pop("feature_index", None)
            rendered_feature_payload["id"] = f"f{idx}"
            rendered_feature_svg_id = rendered_stable_svg_id
            if linear_rendered_feature_ids:
                rendered_feature_svg_id = (
                    make_linear_rendered_feature_id(
                        record_index=rec_idx,
                        stable_feature_id=rendered_stable_svg_id,
                        record_count=len(records),
                    )
                    or rendered_stable_svg_id
                )
            if rendered_feature_svg_id:
                rendered_feature_payload["rendered_feature_svg_id"] = (
                    rendered_feature_svg_id
                )
            features.append(rendered_feature_payload)
            idx += 1

    payload: dict[str, object] = {
        "features": features,
        "record_ids": record_ids,
        "selector_safety_scope": selector_safety_scope,
    }
    if include_biological_features:
        payload["biological_features"] = biological_features
    return payload


def extract_features_from_genbank_payload(
    gb_path: str | Path,
    region_spec: object | None = None,
    record_selector: object | None = None,
    reverse_flag: object | None = None,
    selected_features: object | None = None,
    feature_visibility_rules: list[dict[str, Any]] | None = None,
    specific_color_rules: dict | None = None,
    include_biological_features: bool = False,
) -> dict[str, object]:
    """Extract the Rich Feature Popup payload shape from a GenBank file."""
    from gbdraw.io.record_select import parse_record_selector, reverse_records, select_record

    records = list(SeqIO.parse(str(gb_path), "genbank"))
    selector = parse_record_selector(_normalize_record_selector(record_selector))
    records = select_record(records, selector)
    reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
    records = reverse_records(records, reverse)
    if region_spec:
        from gbdraw.io.regions import apply_region_specs, parse_region_specs

        records = apply_region_specs(records, parse_region_specs([str(region_spec)]))
    return extract_features_from_records_payload(
        records,
        selected_features=selected_features,
        feature_visibility_rules=feature_visibility_rules,
        specific_color_rules=specific_color_rules,
        include_biological_features=include_biological_features,
    )


def extract_features_from_gff_fasta_payload(
    gff_path: str | Path,
    fasta_path: str | Path,
    *,
    mode: str = "linear",
    region_spec: object | None = None,
    record_selector: object | None = None,
    reverse_flag: object | None = None,
    selected_features: object | None = None,
    feature_visibility_rules: list[dict[str, Any]] | None = None,
    specific_color_rules: dict | None = None,
    include_biological_features: bool = False,
) -> dict[str, object]:
    """Extract the Rich Feature Popup payload from paired GFF3 and FASTA files."""

    from gbdraw.io.genome import load_gff_fasta

    normalized_mode = str(mode or "linear").strip().lower()
    if normalized_mode not in {"circular", "linear"}:
        raise ValueError(f"Unsupported GFF3 feature extraction mode: {mode}")

    selector = _normalize_record_selector(record_selector)
    reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
    records = load_gff_fasta(
        [str(gff_path)],
        [str(fasta_path)],
        mode=normalized_mode,
        keep_all_features=True,
        record_selectors=[selector] if selector else None,
        reverse_flags=[reverse],
    )
    if region_spec:
        from gbdraw.io.regions import apply_region_specs, parse_region_specs

        records = apply_region_specs(records, parse_region_specs([str(region_spec)]))
    return extract_features_from_records_payload(
        records,
        selected_features=selected_features,
        feature_visibility_rules=feature_visibility_rules,
        specific_color_rules=specific_color_rules,
        include_biological_features=include_biological_features,
    )


def _read_feature_visibility_rules(path: object | None) -> list[dict[str, Any]] | None:
    normalized_path = _normalize_optional_path(path)
    if not normalized_path:
        return None
    return compile_feature_visibility_rules(read_feature_visibility_file(normalized_path))


def extract_features_from_genbank_json(
    gb_path: str | Path,
    region_spec: object | None = None,
    record_selector: object | None = None,
    reverse_flag: object | None = None,
    selected_features: object | None = None,
    feature_visibility_table_path: object | None = None,
    include_biological_features: bool = False,
) -> str:
    try:
        feature_visibility_rules = _read_feature_visibility_rules(feature_visibility_table_path)
        payload = extract_features_from_genbank_payload(
            gb_path,
            region_spec=region_spec,
            record_selector=record_selector,
            reverse_flag=reverse_flag,
            selected_features=selected_features,
            feature_visibility_rules=feature_visibility_rules,
            include_biological_features=include_biological_features,
        )
    except Exception as exc:
        return json.dumps({"error": str(exc)})
    return json.dumps(payload)


def extract_features_from_gff_fasta_json(
    gff_path: str | Path,
    fasta_path: str | Path,
    mode: str = "linear",
    region_spec: object | None = None,
    record_selector: object | None = None,
    reverse_flag: object | None = None,
    selected_features: object | None = None,
    feature_visibility_table_path: object | None = None,
    include_biological_features: bool = False,
) -> str:
    try:
        payload = extract_features_from_gff_fasta_payload(
            gff_path,
            fasta_path,
            mode=mode,
            region_spec=region_spec,
            record_selector=record_selector,
            reverse_flag=reverse_flag,
            selected_features=selected_features,
            feature_visibility_rules=_read_feature_visibility_rules(
                feature_visibility_table_path
            ),
            include_biological_features=include_biological_features,
        )
    except Exception as exc:
        return json.dumps({"error": str(exc)})
    return json.dumps(payload)
