from __future__ import annotations

from typing import NamedTuple

from Bio.SeqRecord import SeqRecord


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


__all__ = ["RecordSourceMetadata", "infer_record_source_metadata"]
