"""What several records drawn on one circle show as one figure: the definition they share.

The records keep their own identities, features and GC; only the centre text speaks for the circle (how many
records it draws, their bases and GC content, how many short records it leaves out).
"""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord

from ...analysis.gc import gc_percent_of_sequences
from ...core.record_metadata import RecordSourceMetadata, infer_record_source_metadata
from ...render.groups.circular.definition import DefinitionRecordFacts, definition_record_facts
from .parts import CircularOneCircle

# The names of a circle whose records come from different organisms or strains: no name line.
_NO_NAMES = RecordSourceMetadata(organism="", strain="", replicon=None, organelle=None)


def one_circle_definition_facts(circle: CircularOneCircle, *, definition_shown: bool) -> DefinitionRecordFacts:
    """The definition facts of a circle of records: the names every drawn record shares (none when they differ),
    without a replicon or organelle line, the drawn records' bases and their GC content, and a line that counts the
    records. A hidden definition (``definition_shown`` false) takes no GC content, so no base is counted for it, and
    keeps only that line when the circle leaves records out (``DefinitionRecordFacts.summary_only``), so they never
    vanish without a word."""

    facts = [definition_record_facts(unit.record.gb_record, gc_percent=False) for unit in circle.units]
    count = len(circle.units)
    summary = f"{count:,} record{'s' if count != 1 else ''}"
    if circle.hidden_count:
        hidden = circle.hidden_count
        summary += f" (+{hidden:,} under {circle.min_length_bp:,} bp not drawn)"
    first = facts[0]
    shared = {(fact.metadata.organism, fact.metadata.strain) for fact in facts}
    return DefinitionRecordFacts(
        record_id=first.record_id,
        metadata=first.metadata._replace(replicon=None, organelle=None) if len(shared) == 1 else _NO_NAMES,
        length=sum(fact.length for fact in facts),
        gc_percent=(
            gc_percent_of_sequences(unit.record.gb_record.seq for unit in circle.units) if definition_shown else None
        ),
        summary=summary,
        summary_only=not definition_shown and bool(circle.hidden_count),
    )


# The standard abbreviations of the organelles a record names (its source's /organelle value).
_ORGANELLE_ABBREVIATIONS = {"mitochondrion": "mt", "chloroplast": "cp", "plastid": "cp", "plastid:chloroplast": "cp"}


def record_band_names(record: SeqRecord, ordinal: int) -> tuple[str, ...]:
    """The names a record can show on the record band, in the order the band tries them until one fits its arc:
    its chromosome ("IV"), plasmid or organelle; the organelle's standard abbreviation ("mt" for a mitochondrion,
    "cp" for a chloroplast or plastid); its place among the input's records (``ordinal``, from 1)."""

    metadata = infer_record_source_metadata(record)
    organelle = metadata.organelle or ""
    full = next((name for name in (metadata.chromosome, metadata.replicon, organelle.capitalize()) if name), None)
    names = (full, _ORGANELLE_ABBREVIATIONS.get(organelle.lower()), str(ordinal))
    return tuple(dict.fromkeys(name for name in names if name))


__all__ = ["one_circle_definition_facts", "record_band_names"]
