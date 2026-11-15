"""Open reading frames predicted from a nucleotide sequence: gbdraw's one ORF owner.

The rules are those of ``orfind.py`` (``satoshikawato/bio_small_scripts``, MIT), which read the six frames of a
record: an ORF begins at the first start codon after the previous in-frame stop codon (or after the frame's
beginning), ends with the next stop codon, and is kept when its protein, stop excluded, has at least the minimum
number of residues. Unless asked to keep them, an ORF whose span lies inside another ORF's span, on any strand or
frame, is left out. Four changes are made to those rules:

- Contained ORFs are found by one sort and sweep, O(n log n), not by comparing every pair (OV-2479). An ORF whose
  span is the same as another's is kept with it: neither lies inside the other (orfind dropped both, OV-3568).
- On a circular record, the frames run on across the origin, so an ORF may cross it. Its ``end`` is then above the
  record length (unrolled coordinates). An ORF is never longer than the record; when the open stretch before a stop
  codon is longer, the ORF begins at the first start codon within one record length of the stop's end. A frame with
  no stop codon anywhere on the circle has no ORF.
- On a linear record, an ORF that reaches the sequence end without a stop codon is kept and marked 3'-partial: it
  runs to the record end (base ``len`` on the plus strand, base 1 on the minus strand), as in orfind. Each ORF begins
  with a start codon inside the sequence, so no ORF is 5'-partial; the open stretch before a frame's first start
  codon is not an ORF, as in orfind.
- The start codons are those of the genetic code table (``"table"``) or, as in orfind, the codons that the table
  translates to methionine (``"methionine"``, the default). Both sets come from Biopython's codon tables, as do the stop codons
  and the translation, including the IUPAC ambiguity codes Biopython resolves. A frame without any start codon has
  no ORF (orfind read one from the frame's last codon, with a negative start, OV-3567).

The defaults are the Owner's (D-126, 2026-11-14): proteins of at least 50 residues, methionine start codons, and the
genetic code the record's CDS all state (``/transl_table``), else table 11 (``record_genetic_code``); orfind read
table 1. Predicted ORFs are drawn as ``CDS`` (D-87) in their own colour and Legend row, "Predicted ORF", both read
from their typed mark (``PredictedFeatures``, ``predicted_feature_lookup``), and carry no qualifier that gbdraw shows
as a label.

Codon meanings are read once per table from Biopython's ``Seq.translate`` into lookup arrays, and each frame is then
read with array operations, so a 13 Mb record takes seconds.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import cache
from collections.abc import Callable, Iterable, Sequence
from typing import Literal, get_args

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq, Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
)
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray

from gbdraw.core.record_metadata import _feature_source_index_map
from gbdraw.exceptions import ValidationError
from gbdraw.version import __version__

StartCodonSet = Literal["methionine", "table"]
RecordGeneticCode = Literal["record"]
Topology = Literal["linear", "circular"]
Strand = Literal[1, -1]

# One code per IUPAC nucleotide letter (U reads as T); every other byte is the last code, an invalid base.
_LETTERS = "ACGTRYSWKMBDHVN"
_COMPLEMENTS = "TGCAYRSWMKVHDBN"
_INVALID = len(_LETTERS)
_CODES = 16
_CODONS = _CODES**3

RECORD_FALLBACK_GENETIC_CODE = 11
"""The table of a record whose CDS do not all state one table, and of a record without CDS (a FASTA record)."""
_INFERENCE_PREFIX = "ab initio prediction:gbdraw:"
PREDICTED_ORF_FEATURE_TYPE = "CDS"
"""The feature type of a predicted ORF (D-87): the Legend row of predicted ORFs follows the rules of this type."""
PREDICTED_ORF_CAPTION = "Predicted ORF"
PREDICTED_ORF_COLOR = "#e58fa6"
# The record annotation that names the features ORF prediction added to a record (``PredictedFeatures``).
PREDICTED_FEATURES_ANNOTATION = "gbdraw_predicted_features"
"""A rose that no packaged palette gives another feature type (tRNA is #e8b441 in the default palette)."""


def _byte_codes() -> bytes:
    table = bytearray([_INVALID]) * 256
    for code, letter in enumerate(_LETTERS):
        table[ord(letter)] = table[ord(letter.lower())] = code
    table[ord("U")] = table[ord("u")] = _LETTERS.index("T")
    return bytes(table)


_BYTE_CODES = _byte_codes()
_COMPLEMENT_CODES = np.array(
    [_LETTERS.index(letter) for letter in _COMPLEMENTS] + [_INVALID], dtype=np.uint8
)


@dataclass(frozen=True, slots=True)
class OrfPredictionParameters:
    """How ORFs are predicted. The defaults are the Owner's (D-126): the record's table, methionine starts, 50
    residues, contained ORFs left out."""

    genetic_code: int | RecordGeneticCode = "record"
    """An NCBI genetic code table id, as Biopython's codon tables number them, or ``"record"``: the table the
    record's CDS all state, else 11 (``record_genetic_code``)."""
    start_codons: StartCodonSet = "methionine"
    """``"methionine"``: the codons the table translates to M (orfind); ``"table"``: the table's start codons."""
    min_protein_length: int = 50
    """The fewest residues an ORF's protein has, its stop codon not counted."""
    keep_nested: bool = False
    """Keep an ORF whose span lies inside another ORF's span."""

    def __post_init__(self) -> None:
        if self.genetic_code != "record" and not _known_table(self.genetic_code):
            known = ", ".join(
                str(table) for table in sorted(CodonTable.unambiguous_dna_by_id)
            )
            raise ValidationError(
                f"Unknown genetic code {self.genetic_code!r}; the tables are {known}, or 'record'."
            )
        if self.start_codons not in get_args(StartCodonSet):
            raise ValidationError(
                f"Start codons must be one of {get_args(StartCodonSet)}, not {self.start_codons!r}."
            )
        if self.min_protein_length < 1:
            raise ValidationError(
                f"The minimum protein length is at least 1, not {self.min_protein_length}."
            )


def _known_table(value: object) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and value in CodonTable.unambiguous_dna_by_id
    )


def record_genetic_code(features: Iterable[SeqFeature]) -> int:
    """The genetic code of a record: the table its CDS all state, each with one ``/transl_table`` naming a known
    table, else ``RECORD_FALLBACK_GENETIC_CODE`` (mixed tables, a CDS stating none, or no CDS; D-126 Q6=C)."""

    stated: set[int] = set()
    for feature in features:
        if feature.type != "CDS":
            continue
        values = feature.qualifiers.get("transl_table", [])
        text = str(values[0]).strip() if len(values) == 1 else ""
        table = int(text) if text.isdigit() else None
        if not _known_table(table):
            return RECORD_FALLBACK_GENETIC_CODE
        stated.add(int(text))
    return stated.pop() if len(stated) == 1 else RECORD_FALLBACK_GENETIC_CODE


@dataclass(frozen=True, slots=True)
class PredictedOrf:
    """One predicted ORF, in 0-based, end-exclusive coordinates of the record's plus strand."""

    start: int
    """The lowest base of the ORF (its start codon on the plus strand, its stop codon on the minus strand)."""
    end: int
    """One past the highest base; above the record length when the ORF crosses the origin of a circular record."""
    strand: Strand
    partial_3prime: bool
    """No stop codon: the ORF runs off the end of a linear record (its 3' end, the record end on its strand)."""
    protein: str
    """One residue per codon from the start codon, read as M, to the stop codon, which is not included."""


@dataclass(frozen=True, slots=True)
class _CodonMeanings:
    residues: NDArray[np.uint8]
    stops: NDArray[np.bool_]
    starts: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class _Frame:
    """One reading frame of one strand: its codons in reading order. Codon ``k`` begins at strand base
    ``(first_base + 3 * (k + shift)) % length``."""

    codons: NDArray[np.uint16]
    first_base: int
    shift: int
    length: int
    max_codons: int | None
    """Circular: the most codons, stop included, an ORF has (one record length); the frame ends with a stop codon,
    so no ORF runs off it. Linear: None."""

    def bases(self, codons: NDArray[np.intp]) -> NDArray[np.intp]:
        return (self.first_base + 3 * (codons + self.shift)) % self.length


@dataclass(frozen=True, slots=True)
class _Candidates:
    """The ORFs of one frame, as codon indexes in it."""

    frame: _Frame
    first: NDArray[np.intp]
    """The start codon of each ORF."""
    last: NDArray[np.intp]
    """The stop codon of each ORF, or the frame length when it has none (3'-partial)."""


@cache
def _codon_meanings(genetic_code: int, start_codons: StartCodonSet) -> _CodonMeanings:
    """What every codon of the 16 codes means in the table, as Biopython's ``Seq.translate`` reads it."""

    letters = _LETTERS + "?"
    table_starts = set(CodonTable.ambiguous_dna_by_id[genetic_code].start_codons)
    residues = np.full(_CODONS, ord("X"), dtype=np.uint8)
    starts = np.zeros(_CODONS, dtype=np.bool_)
    for index in range(_CODONS):
        codon = letters[index // 256] + letters[index // 16 % 16] + letters[index % 16]
        try:
            residue = str(Seq(codon).translate(table=genetic_code))
        except CodonTable.TranslationError:
            continue
        residues[index] = ord(residue)
        starts[index] = residue == "M" or (
            start_codons == "table" and codon in table_starts
        )
    return _CodonMeanings(residues=residues, stops=residues == ord("*"), starts=starts)


def _plus_strand_codes(sequence: Seq | MutableSeq | str | bytes) -> NDArray[np.uint8]:
    raw = (
        bytes(sequence)
        if isinstance(sequence, (Seq, MutableSeq, bytes, bytearray))
        else sequence.encode("utf-8")
    )
    return np.frombuffer(raw.translate(_BYTE_CODES), dtype=np.uint8)


def _frames(
    codes: NDArray[np.uint8], topology: Topology, meanings: _CodonMeanings
) -> list[_Frame]:
    """The reading frames of one strand. A circular strand whose length is a multiple of 3 has three frames that
    each return to their first codon, otherwise one frame that visits every base (the three base frames one after
    the other); each is read from the codon after its first stop codon, so that it ends with a stop, and a frame
    without one is left out."""

    length = len(codes)
    wide = codes.astype(np.uint16)
    if topology == "linear":
        codon_at = wide[:-2] * 256 + wide[1:-1] * 16 + wide[2:]
        return [
            _Frame(codon_at[frame::3], frame, 0, length, max_codons=None)
            for frame in range(3)
        ]
    wrapped = np.concatenate((wide, wide[:2]))
    codon_at = wrapped[:-2] * 256 + wrapped[1:-1] * 16 + wrapped[2:]
    if length % 3 == 0:
        chains = [(frame, codon_at[frame::3]) for frame in range(3)]
    else:
        chains = [
            (
                0,
                np.concatenate(
                    [codon_at[(-turn * length) % 3 :: 3] for turn in range(3)]
                ),
            )
        ]
    frames = []
    for first_base, codons in chains:
        stops = np.flatnonzero(meanings.stops[codons])
        if len(stops):
            shift = int(stops[0]) + 1
            rotated = np.concatenate((codons[shift:], codons[:shift]))
            frames.append(
                _Frame(rotated, first_base, shift, length, max_codons=length // 3)
            )
    return frames


def _frame_candidates(
    frame: _Frame, meanings: _CodonMeanings, min_protein_length: int
) -> _Candidates:
    """Each stretch after a stop codon (or the frame's beginning) gives the ORF from its first start codon to the
    stretch's stop codon (or the frame's end), when its protein is long enough."""

    stops = np.flatnonzero(meanings.stops[frame.codons])
    starts = np.flatnonzero(meanings.starts[frame.codons])
    last = np.append(stops, len(frame.codons)).astype(np.intp)
    lowest = np.append(0, stops + 1).astype(np.intp)
    if frame.max_codons is not None:
        lowest = np.maximum(lowest, last + 1 - frame.max_codons)
    if not len(starts):
        return _Candidates(
            frame, np.empty(0, dtype=np.intp), np.empty(0, dtype=np.intp)
        )
    found = np.searchsorted(starts, lowest)
    first = starts[np.minimum(found, len(starts) - 1)]
    keep = (found < len(starts)) & (first < last) & (last - first >= min_protein_length)
    return _Candidates(frame, first[keep], last[keep])


def _strand_spans(candidates: _Candidates) -> tuple[NDArray[np.intp], NDArray[np.intp]]:
    """The first base and the base count of each ORF on its own strand (unrolled past the origin when circular)."""

    frame = candidates.frame
    starts = frame.bases(candidates.first)
    if frame.max_codons is not None:
        return starts, (candidates.last - candidates.first + 1) * 3
    stopped = candidates.last < len(frame.codons)
    return starts, np.where(
        stopped, frame.bases(candidates.last) + 3, frame.length
    ) - starts


def _contained(
    starts: NDArray[np.intp], ends: NDArray[np.intp], length: int
) -> NDArray[np.bool_]:
    """Whether each span lies inside another, different span: one sort, then a sweep of the running end.

    A span over the origin (``end`` above ``length``) also contains, through a copy shifted back by ``length``,
    the spans near the beginning of the record that lie inside its part after the origin.
    """

    over = ends > length
    all_starts = np.concatenate((starts, starts[over] - length))
    all_ends = np.concatenate((ends, ends[over] - length))
    order = np.lexsort((-all_ends, all_starts))
    sorted_starts = all_starts[order]
    sorted_ends = all_ends[order]
    count = len(order)
    if not count:
        return np.zeros(0, dtype=np.bool_)
    ends_before = np.append(
        np.iinfo(np.intp).min, np.maximum.accumulate(sorted_ends)[:-1]
    )
    new_span = np.ones(count, dtype=np.bool_)
    new_span[1:] = (sorted_starts[1:] != sorted_starts[:-1]) | (
        sorted_ends[1:] != sorted_ends[:-1]
    )
    span_first = np.maximum.accumulate(np.where(new_span, np.arange(count), 0))
    contained = np.empty(count, dtype=np.bool_)
    contained[order] = ends_before[span_first] >= sorted_ends
    return contained[: len(starts)]


def predict_orfs(
    sequence: Seq | MutableSeq | str | bytes,
    *,
    topology: Topology,
    parameters: OrfPredictionParameters,
) -> tuple[PredictedOrf, ...]:
    """The ORFs of a record's sequence on both strands, sorted by start, end and strand (minus first)."""

    if topology not in get_args(Topology):
        raise ValidationError(
            f"Topology must be one of {get_args(Topology)}, not {topology!r}."
        )
    if parameters.genetic_code == "record":
        raise ValidationError(
            "The genetic code 'record' is read from a record's CDS: predict a record's ORFs with "
            "predict_record_orf_features, or give a table."
        )
    meanings = _codon_meanings(parameters.genetic_code, parameters.start_codons)
    plus = _plus_strand_codes(sequence)
    if len(plus) < 3:
        return ()
    strands: tuple[tuple[Strand, NDArray[np.uint8]], ...] = (
        (1, plus),
        (-1, _COMPLEMENT_CODES[plus[::-1]]),
    )
    found = [
        (strand, _frame_candidates(frame, meanings, parameters.min_protein_length))
        for strand, codes in strands
        for frame in _frames(codes, topology, meanings)
    ]
    return _orfs(found, len(plus), meanings, keep_nested=parameters.keep_nested)


def _orfs(
    found: list[tuple[Strand, _Candidates]],
    length: int,
    meanings: _CodonMeanings,
    *,
    keep_nested: bool,
) -> tuple[PredictedOrf, ...]:
    """The candidates in plus-strand coordinates, those inside another left out unless kept, with their proteins."""

    spans = [
        (strand, candidates, *_strand_spans(candidates)) for strand, candidates in found
    ]
    starts = np.concatenate(
        [np.empty(0, np.intp)]
        + [
            first if strand == 1 else (length - first - size) % length
            for strand, _, first, size in spans
        ]
    )
    sizes = np.concatenate([np.empty(0, np.intp)] + [size for _, _, _, size in spans])
    ends = starts + sizes
    kept = (
        np.ones(len(starts), dtype=np.bool_)
        if keep_nested
        else ~_contained(starts, ends, length)
    )
    orfs = []
    offset = 0
    for strand, candidates, _, _ in spans:
        residues = meanings.residues[candidates.frame.codons]
        frame_length = len(candidates.frame.codons)
        for row, (first, last) in enumerate(
            zip(candidates.first.tolist(), candidates.last.tolist())
        ):
            if kept[offset + row]:
                protein = "M" + residues[first + 1 : last].tobytes().decode("ascii")
                orfs.append(
                    PredictedOrf(
                        int(starts[offset + row]),
                        int(ends[offset + row]),
                        strand,
                        last == frame_length,
                        protein,
                    )
                )
        offset += len(candidates.first)
    return tuple(sorted(orfs, key=lambda orf: (orf.start, orf.end, orf.strand)))


def _orf_location(
    orf: PredictedOrf, record_length: int
) -> SimpleLocation | CompoundLocation:
    """The ORF's location; a 3'-partial end is fuzzy (``>`` on the plus strand, ``<`` on the minus strand)."""

    if orf.end > record_length:
        parts = [
            SimpleLocation(orf.start, record_length, orf.strand),
            SimpleLocation(0, orf.end - record_length, orf.strand),
        ]
        return CompoundLocation(parts if orf.strand == 1 else parts[::-1])
    if not orf.partial_3prime:
        return SimpleLocation(orf.start, orf.end, orf.strand)
    if orf.strand == 1:
        return SimpleLocation(orf.start, AfterPosition(orf.end), orf.strand)
    return SimpleLocation(BeforePosition(orf.start), orf.end, orf.strand)


def predicted_orf_features(
    orfs: Sequence[PredictedOrf], *, record_length: int, genetic_code: int
) -> list[SeqFeature]:
    """The ORFs as ``CDS`` features (D-87).

    ``transl_table`` and ``translation`` are the INSDC qualifiers that let a reader translate the CDS again and
    compare proteins (gbdraw's protein comparisons read ``translation`` first, ``transl_table`` when they translate);
    ``inference`` (``ab initio prediction:gbdraw:<version>``, the INSDC vocabulary for a prediction) says that the
    CDS is predicted, not annotated, and is the qualifier a colour rule can match. No qualifier that gbdraw shows as a
    label (``product``, ``gene``, ``locus_tag``, ``note``, ...) is written: orfind's ORF ids as labels name nothing,
    and gbdraw finds a feature by its own identity, not by a name.
    """

    inference = f"{_INFERENCE_PREFIX}{__version__}"
    return [
        SeqFeature(
            _orf_location(orf, record_length),
            type=PREDICTED_ORF_FEATURE_TYPE,
            qualifiers={
                "transl_table": [str(genetic_code)],
                "translation": [orf.protein],
                "inference": [inference],
            },
        )
        for orf in orfs
    ]


def predict_record_orf_features(
    record: SeqRecord, *, topology: Topology, parameters: OrfPredictionParameters
) -> list[SeqFeature]:
    """The ORFs of a record as ``CDS`` features, read with the parameters' table or, for ``"record"``, the table
    of ``record_genetic_code``: the one entry the record loading path calls."""

    table = (
        record_genetic_code(record.features)
        if parameters.genetic_code == "record"
        else parameters.genetic_code
    )
    orfs = predict_orfs(
        record.seq, topology=topology, parameters=replace(parameters, genetic_code=table)
    )
    return predicted_orf_features(orfs, record_length=len(record.seq), genetic_code=table)


@dataclass(frozen=True)
class PredictedFeatures:
    """The source feature indexes of the features ORF prediction added to a record: the typed mark of a predicted
    feature, kept on the record (``PREDICTED_FEATURES_ANNOTATION``). A crop or reverse complement copies the
    record's annotations and carries each feature's source feature index, so the mark holds on every drawn unit."""

    source_feature_indexes: frozenset[int]


def predicted_source_feature_indexes(record: SeqRecord) -> frozenset[int]:
    """The source feature indexes of the predicted features of ``record`` (``PredictedFeatures``); none without."""

    mark = (getattr(record, "annotations", None) or {}).get(PREDICTED_FEATURES_ANNOTATION)
    return mark.source_feature_indexes if isinstance(mark, PredictedFeatures) else frozenset()


def predicted_feature_lookup(record: SeqRecord) -> Callable[[SeqFeature], bool]:
    """Whether a feature of ``record`` is predicted, by its source feature index against the record's typed mark
    (``PredictedFeatures``): the one test the colour, the Legend row and the lanes of a predicted ORF read. A
    ``CDS`` of the file's own whose ``inference`` names gbdraw is an annotated ``CDS``."""

    indexes = predicted_source_feature_indexes(record)
    if not indexes:
        return lambda _feature: False
    source_indexes = _feature_source_index_map(record.features)
    return lambda feature: source_indexes.get(id(feature)) in indexes
