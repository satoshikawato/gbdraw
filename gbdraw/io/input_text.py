"""The text the parsers read from an input file, and the feature types they give.

GenBank, DDBJ and EMBL flat files are read with Biopython, GFF3 with BCBio-GFF.
Some writers lay out lines those parsers refuse or misread; this module is the one
owner of the repairs made to the text before they read it, and of the feature-type
spelling after they read it. A line that needs no repair reaches the parser as
written, so the text is streamed, never held twice.
"""

from __future__ import annotations

import io
import re
from collections.abc import Iterator
from typing import Literal, TextIO

from Bio.SeqFeature import SeqFeature

FlatFileFormat = Literal["genbank", "embl"]
FLAT_FILE_FORMAT_NAMES: dict[FlatFileFormat, str] = {"genbank": "GenBank", "embl": "EMBL"}

_BYTE_ORDER_MARK = "\ufeff"
# A GenBank record starts at a LOCUS line, an EMBL record at an ID line.
_EMBL_RECORD_START = "ID   "
# Geneious writes the sequence name, spaces included, as the LOCUS name.
_LOCUS_NAME_AND_LENGTH = re.compile(r"^(LOCUS +)(\S.*?)( +\d+ (?:bp|aa)\b.*)$", re.DOTALL)
# Feature-table columns: the key starts at column 6, a qualifier or location at column 22.
_FEATURE_KEY_INDENT = 5
_FEATURE_QUALIFIER_INDENT = 21
_FEATURE_KEY_AND_LOCATION = re.compile(r"^(\S+)[ \t]+(.*)$", re.DOTALL)
# GenBank: the feature table runs from FEATURES to the first line that starts in column 1.
# EMBL: every feature-table line starts with FT.
_GENBANK_FEATURES = "FEATURES"
_EMBL_FEATURE_LINE = "FT"

# INSDC feature keys (Feature Table Definition, section 7.2). GFF3 writers and some
# GenBank writers spell them in another case (`cds`, `trna`); the readers give each
# feature the INSDC spelling, so type filters, colours and legends find it.
_INSDC_FEATURE_KEYS = (
    "assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", "exon", "gap",
    "gene", "iDNA", "intron", "J_segment", "mat_peptide", "misc_binding", "misc_difference",
    "misc_feature", "misc_recomb", "misc_RNA", "misc_structure", "mobile_element",
    "modified_base", "mRNA", "ncRNA", "N_region", "old_sequence", "operon", "oriT",
    "polyA_site", "precursor_RNA", "prim_transcript", "primer_bind", "propeptide",
    "protein_bind", "regulatory", "repeat_region", "rep_origin", "rRNA", "S_region",
    "sig_peptide", "source", "stem_loop", "STS", "telomere", "tmRNA", "transit_peptide",
    "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR", "5'UTR",
)
_CANONICAL_FEATURE_TYPES = {key.lower(): key for key in _INSDC_FEATURE_KEYS}


def flat_file_format(first_line: str) -> FlatFileFormat:
    """The flat-file format a file's first non-blank line names: EMBL for an ID line, else GenBank/DDBJ."""

    return "embl" if first_line.startswith(_EMBL_RECORD_START) else "genbank"


def canonical_feature_type(feature_type: str) -> str:
    """The INSDC spelling of a feature key written in another case; any other type as written."""

    return _CANONICAL_FEATURE_TYPES.get(feature_type.lower(), feature_type)


def canonicalize_feature_types(features: list[SeqFeature]) -> None:
    """Give every feature, nested ones included, the INSDC spelling of its type, in place."""

    for feature in features:
        feature.type = canonical_feature_type(str(feature.type))
        sub_features = getattr(feature, "sub_features", None)
        if sub_features:
            canonicalize_feature_types(sub_features)


def _feature_table_line(body: str) -> str:
    """A feature-table line whose layout columns hold tabs, laid out with spaces.

    A tab inside a qualifier value stays: only the indentation, and the gap between
    a feature key and its location, are layout.
    """

    text = body.lstrip(" \t")
    indent = body[: len(body) - len(text)]
    if text.startswith("/") or len(indent.expandtabs(8)) >= _FEATURE_QUALIFIER_INDENT:
        return " " * _FEATURE_QUALIFIER_INDENT + text if "\t" in indent else body
    key_and_location = _FEATURE_KEY_AND_LOCATION.match(text)
    if key_and_location is None:
        return body
    key, location = key_and_location.groups()
    gap = text[len(key): len(text) - len(location)]
    if "\t" not in indent and "\t" not in gap:
        return body
    return " " * _FEATURE_KEY_INDENT + key.ljust(_FEATURE_QUALIFIER_INDENT - _FEATURE_KEY_INDENT - 1) + " " + location


def _locus_line(line: str) -> str:
    """A LOCUS line whose name holds spaces, with each run of them written as one underscore."""

    match = _LOCUS_NAME_AND_LENGTH.match(line)
    if match is None or not any(character.isspace() for character in match.group(2)):
        return line
    prefix, name, rest = match.groups()
    return prefix + "_".join(name.split()) + rest


def _normalized_flat_lines(lines: Iterator[str], fmt: FlatFileFormat) -> Iterator[str]:
    in_genbank_features = False
    for line in lines:
        if fmt == "embl":
            if "\t" in line and line.startswith(_EMBL_FEATURE_LINE):
                line = _EMBL_FEATURE_LINE + _feature_table_line("  " + line[2:])[2:]
        elif line.startswith("LOCUS"):
            in_genbank_features = False
            line = _locus_line(line)
        elif line.startswith(_GENBANK_FEATURES):
            in_genbank_features = True
        elif in_genbank_features:
            if line[:1] not in (" ", "\t"):
                in_genbank_features = False
            elif "\t" in line:
                line = _feature_table_line(line)
        yield line


def _gff3_attributes(column: str) -> str:
    """GFF3 attributes with the spaces around each ``;`` removed.

    BCBio-GFF splits the attribute column at ``" ; "`` when it holds one, else at
    ``";"``, so a column that mixes ``"; "`` and ``" ; "`` (exonerate) stops it with
    an AssertionError. In GFF3 a ``;`` always separates attributes: a value writes
    one as ``%3B``.
    """

    if " ;" not in column or "=" not in column.split(";", 1)[0]:
        return column
    return ";".join(part.strip() for part in column.split(";"))


def _normalized_gff3_lines(lines: Iterator[str]) -> Iterator[str]:
    for line in lines:
        if line.startswith("##FASTA"):
            yield line
            yield from lines
            return
        if not line.startswith("#") and line.count("\t") >= 8:
            columns = line.rstrip("\r\n").split("\t")
            attributes = _gff3_attributes(columns[8])
            if attributes != columns[8]:
                columns[8] = attributes
                line = "\t".join(columns) + "\n"
        yield line


def text_lines(handle: TextIO) -> Iterator[str]:
    """The lines of an open text file, without a UTF-8 byte order mark at its start (OV-2354).

    A file opened as UTF-8 keeps the mark as U+FEFF. Biopython reads it as text
    before the first record and skips that record, its FASTA reader skips the
    first entry, and BCBio-GFF reads it into the first line. Every reader of an
    input file reads its lines here first.
    """

    lines = iter(handle.readline, "")
    first = next(lines, "")
    if first:
        yield first.removeprefix(_BYTE_ORDER_MARK)
        yield from lines


class LineStream:
    """A read-only text handle over a line iterator, as Biopython and BCBio-GFF read one."""

    def __init__(self, lines: Iterator[str]) -> None:
        self._lines = lines

    def __iter__(self) -> Iterator[str]:
        return self._lines

    def __next__(self) -> str:
        return next(self._lines)

    def readline(self) -> str:
        return next(self._lines, "")

    def read(self, size: int = -1) -> str:
        if size == 0:  # Biopython reads nothing to tell a text handle from a binary one.
            return ""
        if size < 0:
            return "".join(self._lines)
        raise io.UnsupportedOperation("read(size) of a line stream")


def flat_file_stream(lines: Iterator[str], fmt: FlatFileFormat) -> LineStream:
    """The flat-file text Biopython reads from ``lines``: tab layout and spaced LOCUS names repaired."""

    return LineStream(_normalized_flat_lines(lines, fmt))


def gff3_stream(lines: Iterator[str]) -> LineStream:
    """The GFF3 text BCBio-GFF reads from ``lines``: mixed attribute separators repaired."""

    return LineStream(_normalized_gff3_lines(lines))


__all__ = [
    "FLAT_FILE_FORMAT_NAMES",
    "FlatFileFormat",
    "LineStream",
    "canonical_feature_type",
    "canonicalize_feature_types",
    "flat_file_format",
    "flat_file_stream",
    "gff3_stream",
    "text_lines",
]
