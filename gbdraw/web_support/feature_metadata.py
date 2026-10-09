from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence, cast

from Bio import SeqIO

from gbdraw.core.record_metadata import (
    _absolute_display_interval,
    _iter_source_features as _iter_features,
    _read_coord_map as _read_record_coord_map,
    _source_feature_anchor_profile,
    _source_feature_index,
    _source_feature_location_parts,
)
from gbdraw.core.sequence import translate_cds
from gbdraw.features.overrides import feature_override_lookup
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
from gbdraw.svg.ids import instance_svg_id
from gbdraw.web_support.error_adapter import serialize_web_error

if TYPE_CHECKING:
    from gbdraw.features.placement import ResolvedRecordFeatureInputs


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
) -> tuple[dict[str, object], str, str, dict[str, str | None]]:
    """Return source selector values, the biological and processed-record
    feature IDs, and the drawn record's selector values.

    The drawn values are the ones the renderer's rule matching reads from the
    processed (cropped, reverse-complemented) record (feature catalog 5,
    ``drawnSelector``).
    """

    selector = build_feature_selector_values(feature, record_id=record_id)
    rendered_feature_id = str(selector.get("hash") or "")
    drawn_selector = {
        "hash": rendered_feature_id or None,
        "location": str(selector.get("location") or "") or None,
        "recordLocation": str(selector.get("record_location") or "") or None,
    }
    rendered_parts = _biological_location_parts(
        feature.location,
        coord_base,
        coord_step,
    )
    parts = _source_feature_location_parts(feature) or tuple(rendered_parts)
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
    return selector, stable_feature_id, rendered_feature_id, drawn_selector



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


def _source_anchor_profile(feature: Any, *, coord_step: int = 1) -> dict[str, str]:
    """Return source capability facts without resolving any anchor coordinate."""

    profile = _source_feature_anchor_profile(feature, coord_step=coord_step)
    return {
        "precision": profile.precision,
        "operator": profile.operator,
        "partOrder": profile.part_order,
        "strand": profile.strand,
    }


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

    try:
        return translate_cds(feature, nucleotide_sequence, require_whole_codons=True), warnings
    except Exception as exc:
        warnings.append(f"CDS translation skipped: {exc}")
        return "", warnings


def extract_features_from_records_payload(
    records: Any,
    *,
    selected_features: object | None = None,
    feature_visibility_rules: list[dict[str, Any]] | None = None,
    record_features: Sequence[ResolvedRecordFeatureInputs] = (),
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
    idx = 0
    biological_idx = 0
    for rec_idx, record in enumerate(records):
        record_id = record.id or f"Record_{rec_idx}"
        hash_record_id = record.id
        organism = _get_record_organism(record)
        coord_base, coord_step = _read_record_coord_map(record)
        record_ids.append(record_id)
        prepared_features = []
        rendered_id_counts: Counter[str] = Counter()
        override_of = feature_override_lookup(
            record, record_features[rec_idx].overrides if record_features else None
        )
        for feature_index, feat in enumerate(_iter_features(record.features)):
            source_feature_index = _source_feature_index(feat)
            feature_override = override_of(feat)
            is_rendered_feature = should_render_feature(
                feat,
                selected_feature_set,
                feature_visibility_rules=feature_visibility_rules,
                record_id=hash_record_id,
                specific_color_rules=specific_color_rules,
                feature_override=feature_override,
            )
            if not is_rendered_feature and not include_biological_features:
                continue
            # An undrawn source feature stays out of the catalog unless it has
            # its own Feature visibility, which keeps it in the Web Features
            # list so it can be shown again (R-5).
            if (
                include_biological_features
                and not is_rendered_feature
                and str(getattr(feat, "type", "") or "").lower() == "source"
                and getattr(feature_override, "feature_visibility", None) is None
            ):
                continue
            selector_values = _biological_selector_values(
                feat,
                record_id=hash_record_id,
                coord_base=coord_base,
                coord_step=coord_step,
            )
            if is_rendered_feature and selector_values[2]:
                rendered_id_counts[selector_values[2]] += 1
            prepared_features.append(
                (
                    feature_index,
                    source_feature_index,
                    feat,
                    is_rendered_feature,
                    selector_values,
                )
            )

        for (
            feature_index,
            source_feature_index,
            feat,
            is_rendered_feature,
            selector_values,
        ) in prepared_features:
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

            selector, stable_svg_id, rendered_stable_svg_id, drawn_selector = selector_values

            # The same qualifier map rule matching reads: keys stripped and
            # lowercased, case variants merged in key order (OV-249).
            qualifiers = {
                key: list(values)
                for key, values in cast(
                    "dict[str, list[str]]", selector["qualifiers"]
                ).items()
            }

            feature_payload = {
                "id": f"f{biological_idx}",
                "svg_id": stable_svg_id,
                "stable_svg_id": stable_svg_id,
                "stable_feature_id": stable_svg_id,
                "record_id": record_id,
                "record_idx": rec_idx,
                "feature_index": (
                    feature_index
                    if source_feature_index is None
                    else source_feature_index
                ),
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
                "anchorProfile": _source_anchor_profile(
                    feat,
                    coord_step=coord_step,
                ),
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
            rendered_feature_payload["id"] = f"f{idx}"
            rendered_feature_payload["drawn_selector"] = drawn_selector
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
            if (
                rendered_feature_svg_id
                and rendered_id_counts[rendered_stable_svg_id] > 1
            ):
                rendered_feature_svg_id = instance_svg_id(
                    rendered_feature_svg_id,
                    (
                        feature_index
                        if source_feature_index is None
                        else source_feature_index
                    ),
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

    selector = _normalize_record_selector(record_selector)
    reverse = str(reverse_flag).strip().lower() in {"1", "true", "yes", "y", "on"}
    records = load_gff_fasta(
        [str(gff_path)],
        [str(fasta_path)],
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
        return json.dumps({"error": serialize_web_error(exc, operation="feature-extraction", stage="helper")})
    return json.dumps(payload)


def extract_features_from_gff_fasta_json(
    gff_path: str | Path,
    fasta_path: str | Path,
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
        return json.dumps({"error": serialize_web_error(exc, operation="feature-extraction", stage="helper")})
    return json.dumps(payload)
