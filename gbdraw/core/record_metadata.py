from __future__ import annotations

from typing import NamedTuple

from Bio.SeqRecord import SeqRecord


_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"
_SOURCE_FEATURE_INDEX_ATTR = "_gbdraw_source_feature_index"
_SOURCE_FEATURE_PARTS_ATTR = "_gbdraw_source_feature_location_parts"


def _iter_source_features(features: object):
    """Yield top-level and nested source features in their original order."""
    for feature in features or ():
        yield feature
        yield from _iter_source_features(getattr(feature, "sub_features", None))


def _feature_source_index_map(features: object) -> dict[int, int]:
    """Map nested feature objects to their flattened source-order ordinals."""
    return {id(feature): index for index, feature in enumerate(_iter_source_features(features))}


def _source_feature_index(feature: object) -> int | None:
    """Return a pre-transform source ordinal attached to a feature, if any."""

    try:
        index = int(getattr(feature, _SOURCE_FEATURE_INDEX_ATTR))
    except (AttributeError, TypeError, ValueError):
        return None
    return index if index >= 0 else None


def _source_feature_location_parts(
    feature: object,
) -> tuple[tuple[int, int, int | None], ...] | None:
    """Return pre-transform biological location parts, if attached."""

    value = getattr(feature, _SOURCE_FEATURE_PARTS_ATTR, None)
    if not isinstance(value, tuple):
        return None
    parts: list[tuple[int, int, int | None]] = []
    for part in value:
        if not isinstance(part, tuple) or len(part) != 3:
            return None
        try:
            start = int(part[0])
            end = int(part[1])
        except (TypeError, ValueError):
            return None
        strand = part[2]
        if strand not in {-1, 1}:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts) or None


def _mapped_feature_location_parts(
    feature: object,
    *,
    coord_base: int,
    coord_step: int,
) -> tuple[tuple[int, int, int | None], ...]:
    """Map a feature's current parts into its biological source coordinates."""

    location = getattr(feature, "location", None)
    raw_parts = list(getattr(location, "parts", None) or [location])
    parts: list[tuple[int, int, int | None]] = []
    for part in raw_parts:
        if part is None:
            continue
        try:
            start, end = _absolute_display_interval(
                int(part.start),
                int(part.end),
                coord_base,
                coord_step,
            )
        except (AttributeError, TypeError, ValueError):
            continue
        strand = getattr(part, "strand", None)
        if strand in {-1, 1}:
            strand = int(strand) * (1 if int(coord_step) >= 0 else -1)
        else:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts)


def _copy_source_feature_identity(
    source: object,
    target: object,
    *,
    fallback_index: int,
    coord_base: int,
    coord_step: int,
) -> None:
    """Carry source ordinal and biological parts across a transformation."""

    index = _source_feature_index(source)
    setattr(
        target,
        _SOURCE_FEATURE_INDEX_ATTR,
        int(fallback_index) if index is None else index,
    )
    parts = _source_feature_location_parts(source)
    if parts is None:
        parts = _mapped_feature_location_parts(
            source,
            coord_base=coord_base,
            coord_step=coord_step,
        )
    if parts:
        setattr(target, _SOURCE_FEATURE_PARTS_ATTR, parts)


def _read_coord_map(record: object) -> tuple[int, int]:
    annotations = getattr(record, "annotations", None) or {}
    try:
        base = int(annotations.get(_COORD_BASE_KEY, 1))
    except (TypeError, ValueError):
        base = 1
    try:
        step = int(annotations.get(_COORD_STEP_KEY, 1))
    except (TypeError, ValueError):
        step = 1
    if step == 0:
        step = 1
    return base, (1 if step > 0 else -1)


def _write_coord_map(record: object, *, base: int, step: int) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}  # type: ignore[attr-defined]
    record.annotations[_COORD_BASE_KEY] = int(base)  # type: ignore[attr-defined]
    record.annotations[_COORD_STEP_KEY] = 1 if int(step) >= 0 else -1  # type: ignore[attr-defined]


def _absolute_display_interval(
    start: int,
    end: int,
    coord_base: int,
    coord_step: int,
) -> tuple[int, int]:
    if end <= start:
        coord = coord_base + (coord_step * start)
        return coord - 1, coord
    first_coord = coord_base + (coord_step * start)
    last_coord = coord_base + (coord_step * (end - 1))
    return min(first_coord, last_coord) - 1, max(first_coord, last_coord)


class RecordSourceMetadata(NamedTuple):
    organism: str
    strain: str
    replicon: str | None
    organelle: str | None


def infer_record_source_metadata(record: SeqRecord) -> RecordSourceMetadata:
    """Extract source-feature metadata used in definition labels."""
    annotations = getattr(record, "annotations", None) or {}
    organism = str(annotations.get("organism", "") or "").strip()
    strain = ""
    replicon: str | None = None
    organelle: str | None = None

    for feature in getattr(record, "features", []):
        if getattr(feature, "type", None) != "source":
            continue

        qualifiers = getattr(feature, "qualifiers", {}) or {}
        if "organism" in qualifiers and qualifiers["organism"]:
            organism = str(qualifiers["organism"][0]).strip()
        if "isolate" in qualifiers and qualifiers["isolate"]:
            strain = str(qualifiers["isolate"][0]).strip()
        elif "strain" in qualifiers and qualifiers["strain"]:
            strain = str(qualifiers["strain"][0]).strip()

        if "chromosome" in qualifiers and qualifiers["chromosome"]:
            replicon = f"Chromosome {str(qualifiers['chromosome'][0]).strip()}"
        elif "plasmid" in qualifiers and qualifiers["plasmid"]:
            replicon = str(qualifiers["plasmid"][0]).strip()

        if "organelle" in qualifiers and qualifiers["organelle"]:
            organelle = str(qualifiers["organelle"][0]).strip()

    return RecordSourceMetadata(
        organism=organism,
        strain=strain,
        replicon=replicon,
        organelle=organelle,
    )


_NON_ORGANISM_NAMES = frozenset(
    {
        "synthetic construct",
        "artificial sequence",
        "unidentified",
        "unidentified organism",
        "unknown",
        "unknown organism",
        "vector",
        "cloning vector",
        "unidentified cloning vector",
        "expression vector",
    }
)


def format_inferred_definition(metadata: RecordSourceMetadata) -> str:
    """Format organism and strain into publication-style HTML (e.g. <i>Genus species</i> strain)."""
    organism = str(metadata.organism or "").strip()
    strain = str(metadata.strain or "").strip()

    candidatus = False
    clean_org = organism
    if clean_org.lower().startswith("candidatus "):
        candidatus = True
        clean_org = clean_org[11:].strip()

    if not clean_org or clean_org.lower() in _NON_ORGANISM_NAMES:
        return strain

    words = clean_org.split()
    cand_prefix = "Candidatus " if candidatus else ""
    if len(words) >= 2:
        binom = f"<i>{' '.join(words[:2])}</i>"
        rest = " ".join(words[2:])
        if strain and strain.lower() not in rest.lower():
            rest = f"{rest} {strain}".strip()
        return f"{cand_prefix}{binom} {rest}".strip()
    if words:
        base = f"{cand_prefix}<i>{words[0]}</i>"
        return f"{base} {strain}".strip() if strain else base
    return strain


def format_inferred_subtitle(
    metadata: RecordSourceMetadata,
    description: str = "",
) -> str:
    """Infer a concise record subtitle from replicon/organelle metadata or description."""
    if metadata.replicon:
        return str(metadata.replicon).strip()
    if metadata.organelle:
        return str(metadata.organelle).strip().capitalize()

    desc = str(description or "").strip().rstrip(".")
    if not desc:
        return ""

    import re

    desc_lower = desc.lower()
    plasmid_match = re.search(
        r"(?:plasmid\s+([A-Za-z0-9_-]+)|(p[A-Za-z0-9_-]+)\b)",
        desc,
        re.IGNORECASE,
    )
    if plasmid_match:
        p_name = (plasmid_match.group(1) or plasmid_match.group(2) or "").strip()
        return p_name if p_name.lower().startswith("plasmid") else f"Plasmid {p_name}"

    if "complete genome" in desc_lower:
        if "mitochondri" in desc_lower:
            return "Mitochondrion, complete genome"
        if "chloroplast" in desc_lower:
            return "Chloroplast, complete genome"
        return "Complete genome"

    if "complete sequence" in desc_lower:
        return "Complete sequence"

    organism = str(metadata.organism or "").strip()
    if organism and organism.lower() not in _NON_ORGANISM_NAMES:
        stripped = re.sub(
            rf"^{re.escape(organism)}[,\s]*",
            "",
            desc,
            flags=re.IGNORECASE,
        )
        stripped = re.sub(
            r"^(?:DNA|genomic DNA|cDNA)[,\s]*",
            "",
            stripped,
            flags=re.IGNORECASE,
        ).strip()
        if re.search(
            r"(?:gene cluster|biosynthetic gene cluster|cluster|operon)",
            stripped,
            re.IGNORECASE,
        ):
            return stripped[:1].upper() + stripped[1:]

    cluster_match = re.search(
        r"([A-Za-z0-9_-]+(?:\s+[A-Za-z0-9_-]+)*\s+(?:gene cluster|biosynthetic gene cluster|cluster|operon))",
        desc,
        re.IGNORECASE,
    )
    if cluster_match:
        res = cluster_match.group(1).strip()
        return res[:1].upper() + res[1:]

    return ""


__all__ = [
    "RecordSourceMetadata",
    "format_inferred_definition",
    "format_inferred_subtitle",
    "infer_record_source_metadata",
]
