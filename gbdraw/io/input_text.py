"""The text the parsers read from an input file, and the feature types they give.

GenBank, DDBJ and EMBL flat files are read with Biopython, GFF3 with BCBio-GFF.
Some writers lay out lines those parsers refuse or misread; this module is the one
owner of the repairs made to the text before they read it, and of the feature-type
spelling after they read it. A line that needs no repair reaches the parser as
written, so the text is streamed, never held twice. An input file may be
gzip-compressed, whatever its name: ``open_input_text`` and ``read_input_bytes``
expand one, told by its first two bytes, and every reader opens its file there;
the record sequence reader opens its files with ``open_input_bytes``, and expands
a gzip one with ``open_expanded_input``.
A gzip input expands to at most ``MAX_EXPANDED_INPUT_BYTES``: a small file can
hold gigabytes, so reading stops past the cap and names the file. A plain file
has no cap.
"""

from __future__ import annotations

import codecs
import gzip
import io
import itertools
import locale
import os
import re
import zlib
from collections.abc import Callable, Iterator
from os import PathLike
from typing import IO, TYPE_CHECKING, Literal, TextIO, TypeVar

from Bio.SeqFeature import SeqFeature

from gbdraw.exceptions import InputFileError

if TYPE_CHECKING:
    from _typeshed import WriteableBuffer

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
# A line of text whose lines end in "\n" (the last may have no ending).
_LINE = re.compile(r"[^\n]*\n|[^\n]+\Z")

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

# A gzip stream starts with these two bytes (RFC 1952), whatever the file is named.
GZIP_SIGNATURE = b"\x1f\x8b"
# What reading a gzip stream that is cut off or damaged raises, beside the
# ValueError (UnicodeDecodeError) a plain file raises; BadGzipFile is an OSError.
# A reader spreads it into its except tuple: ``except (ValueError, *GZIP_ERRORS)``.
GZIP_ERRORS = (EOFError, zlib.error, gzip.BadGzipFile)

# The most bytes a gzip input file expands to (2 GiB). The Web reads it from the
# generated contracts; its Session caps are separate constants in session-file.js.
MAX_EXPANDED_INPUT_BYTES = 2 * 1024 * 1024 * 1024
_MIB = 1024 * 1024
# Expanded bytes are buffered 128 KiB at a time: with the default 8 KiB, lines read slower than gzip.open's.
_READ_CHUNK = 128 * 1024

_Read = TypeVar("_Read")


class _CappedExpansion(io.RawIOBase):
    """The expanded bytes of a gzip file, refused once they pass ``limit``.

    The one place a gzip input is expanded: ``open_input_text`` and
    ``read_input_bytes`` read through it, so both stop at the same count.
    """

    def __init__(self, path: str | PathLike[str], limit: int) -> None:
        super().__init__()
        self._gzip = gzip.open(path, "rb")
        self._path = os.fspath(path)
        self._limit = limit
        self._expanded = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: WriteableBuffer) -> int:
        # Bytes up to the limit are read; a read at the limit asks for one more byte and refuses when it exists.
        room = self._limit - self._expanded
        view = memoryview(buffer).cast("B")
        count = self._gzip.readinto(view[:room] if room > 0 else view[:1])
        self._expanded += count
        if self._expanded > self._limit:
            raise InputFileError(
                f"{self._path}: this gzip file expands to more than {self._limit / _MIB:g} MiB "
                f"({self._limit:,} bytes), the most gbdraw reads. Decompress it and give the plain file.",
                diagnostic={
                    "code": "INPUT_UNREADABLE", "reason": "GZIP_EXPANDED_SIZE", "limit": self._limit // _MIB,
                    # The file's name, as the Web names it (PD-OI-116); the Worker's render names
                    # the file the user chose in place of the one it staged (request_render), and
                    # the Web error context cuts it (error_adapter.INPUT_NAME_CHARACTERS).
                    "inputName": os.path.basename(self._path),
                },
            )
        return count

    def close(self) -> None:
        if not self.closed:
            self._gzip.close()
        super().close()


def _expanded_gzip(path: str | PathLike[str], capped: bool) -> io.BufferedReader | gzip.GzipFile:
    if not capped:
        return gzip.open(path, "rb")
    return io.BufferedReader(_CappedExpansion(path, MAX_EXPANDED_INPUT_BYTES), buffer_size=_READ_CHUNK)


def open_input_bytes(path: str | PathLike[str]) -> IO[bytes]:
    """A plain input file opened for its bytes, unbuffered: read in order, or at offsets.

    The record sequence reader (:mod:`gbdraw.io.record_sequences`) opens its files
    here; a gzip input is expanded with :func:`open_expanded_input`.
    """

    return open(path, "rb", buffering=0)


def open_expanded_input(path: str | PathLike[str]) -> io.BufferedReader | gzip.GzipFile:
    """A gzip input file's expanded bytes, stopping past ``MAX_EXPANDED_INPUT_BYTES`` as ``open_input_text`` stops."""

    return _expanded_gzip(path, capped=True)


def is_gzip_file(path: str | PathLike[str]) -> bool:
    """Whether the file starts with the gzip signature; its name does not count."""

    with open(path, "rb") as handle:
        return handle.read(len(GZIP_SIGNATURE)) == GZIP_SIGNATURE


def open_input_text(
    path: str | PathLike[str],
    *,
    encoding: str = "utf-8",
    newline: str | None = None,
    capped: bool = True,
) -> TextIO:
    """An input file opened as text, expanded when it is gzip-compressed.

    Every reader of an input file (sequence, annotation, table, depth, Session)
    opens it here, so a gzip file reads wherever a plain one does, streamed as a
    plain one is. Reading a gzip stream that is cut off or damaged raises one of
    ``GZIP_ERRORS``; reading past ``MAX_EXPANDED_INPUT_BYTES`` of it raises an
    ``InputFileError`` that names the file. A Session file opts out with
    ``capped=False``: its cap is the Web's. A sequence or annotation file is read
    through :func:`read_sequence_text` instead.
    """

    if is_gzip_file(path):
        return io.TextIOWrapper(_expanded_gzip(path, capped), encoding=encoding, newline=newline)
    return open(path, "r", encoding=encoding, newline=newline)


def read_sequence_text(path: str | PathLike[str], read: Callable[[Iterator[str]], _Read]) -> _Read:
    """What ``read`` makes of the lines (:func:`text_lines`) of a GenBank, EMBL, GFF3 or FASTA file.

    The file is read as UTF-8. One that is not UTF-8 is read again in the
    locale's encoding, so a file a Windows editor saved in its code page (cp1252,
    cp932) still reads there. Where the locale is UTF-8 (Linux, macOS, the Web)
    the UTF-8 error stands. ``read`` takes every line it needs before it returns.
    """

    def read_in(encoding: str) -> _Read:
        with open_input_text(path, encoding=encoding) as handle:
            return read(text_lines(handle))

    return in_sequence_encoding(read_in)


def in_sequence_encoding(read: Callable[[str], _Read]) -> _Read:
    """``read("utf-8")``, else, when that meets text that is not UTF-8, ``read`` in the locale's encoding.

    The one rule of the encoding a sequence or annotation file is read in
    (:func:`read_sequence_text`, :mod:`gbdraw.io.record_sequences`).
    """

    try:
        return read("utf-8")
    except UnicodeDecodeError:
        encoding = locale.getpreferredencoding(False)
        if codecs.lookup(encoding).name == "utf-8":
            raise
        return read(encoding)


def read_input_bytes(path: str | PathLike[str]) -> bytes:
    """The bytes of an input file, expanded when it is gzip-compressed (see ``open_input_text``)."""

    if is_gzip_file(path):
        with _expanded_gzip(path, capped=True) as handle:
            return handle.read()
    with open(path, "rb") as handle:
        return handle.read()


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


def _normalized_flat_lines(pieces: Iterator[str], fmt: FlatFileFormat) -> Iterator[str]:
    """The flat-file text of ``pieces`` (each a line, or whole lines), lines repaired as Biopython needs them.

    Lines that start with a space and hold no tab are left as they are and leave the
    state as it is, so a piece of only such lines (a feature table read as one
    block) passes whole; the lines of any other piece are repaired one by one.
    """

    in_genbank_features = False
    for piece in pieces:
        if piece[:1] == " " and "\t" not in piece and piece.count("\n ") + 1 >= piece.count("\n"):
            yield piece
            continue
        for line in _LINE.findall(piece) if piece.count("\n") > 1 else (piece,):
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
    """A read-only text handle over pieces of text, as Biopython and BCBio-GFF read one.

    A piece is a line, or whole lines (each ending in ``\\n``): the lines of a piece
    are split at once, so a reader takes each of them from a list.
    """

    def __init__(self, pieces: Iterator[str]) -> None:
        self._pieces = iter(pieces)
        self._lines: Iterator[str] = iter(())

    def __iter__(self) -> Iterator[str]:
        yield from self._lines
        for piece in self._pieces:
            if piece.count("\n") > 1:
                self._lines = iter(_LINE.findall(piece))
                yield from self._lines
            elif piece:
                yield piece

    def __next__(self) -> str:
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    def readline(self) -> str:
        line = next(self._lines, "")
        if line:
            return line
        for piece in self._pieces:
            if piece.count("\n") > 1:
                self._lines = iter(_LINE.findall(piece))
                return next(self._lines)
            if piece:
                return piece
        return ""

    def read(self, size: int = -1) -> str:
        if size == 0:  # Biopython reads nothing to tell a text handle from a binary one.
            return ""
        if size < 0:
            return "".join(itertools.chain(self._lines, self._pieces))
        raise io.UnsupportedOperation("read(size) of a line stream")


def flat_file_stream(pieces: Iterator[str], fmt: FlatFileFormat) -> LineStream:
    """The flat-file text Biopython reads from ``pieces`` (lines, or whole lines): tab layout and spaced LOCUS names
    repaired."""

    return LineStream(_normalized_flat_lines(pieces, fmt))


def gff3_stream(lines: Iterator[str]) -> LineStream:
    """The GFF3 text BCBio-GFF reads from ``lines``: mixed attribute separators repaired."""

    return LineStream(_normalized_gff3_lines(lines))


__all__ = [
    "FLAT_FILE_FORMAT_NAMES",
    "FlatFileFormat",
    "GZIP_ERRORS",
    "GZIP_SIGNATURE",
    "LineStream",
    "MAX_EXPANDED_INPUT_BYTES",
    "canonical_feature_type",
    "canonicalize_feature_types",
    "open_expanded_input",
    "flat_file_format",
    "flat_file_stream",
    "gff3_stream",
    "is_gzip_file",
    "open_input_bytes",
    "open_input_text",
    "read_input_bytes",
    "read_sequence_text",
    "in_sequence_encoding",
    "text_lines",
]
