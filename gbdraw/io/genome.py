#!/usr/bin/env python
# coding: utf-8

import os
import logging
from contextlib import contextmanager
from collections import Counter
from collections.abc import Iterator, Sequence, Set as AbstractSet
from itertools import islice
import string
from typing import Callable, List, Dict, Literal

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from .input_text import (
    FLAT_FILE_FORMAT_NAMES,
    GZIP_ERRORS,
    FlatFileFormat,
    LineStream,
    canonicalize_feature_types,
    flat_file_format,
    flat_file_stream,
    gff3_stream,
    read_input_bytes,
)
from .record_sequences import read_fasta_records, read_flat_records, read_gff3_text
from .record_select import parse_record_selector, reverse_records, select_record
from ..exceptions import GbdrawError, InputFileError, ParseError, ValidationError

from gbdraw.core.record_metadata import _copy_source_feature_identity, _feature_source_index_map
from gbdraw.features.source import SourceFeatureCatalog, capture_source_features

logger = logging.getLogger(__name__)

# A file that cannot be parsed is a user-fixable input: the Web asks to replace it.
_UNREADABLE_INPUT = {"code": "INPUT_UNREADABLE"}
# A flat file the parser refused: cut off or damaged (OV-1258).
_UNREADABLE_FLAT_FILE = {"code": "INPUT_UNREADABLE", "reason": "GENBANK_PARSE"}
_DUPLICATE_FASTA_IDS = {"code": "INPUT_INVALID", "field": "fasta", "reason": "UNIQUE_IDS"}
# OV-1257: a GFF3/FASTA pair whose slot file holds no entry says which slot.
_GFF_WITHOUT_FEATURES = {"code": "NO_RECORDS", "reason": "GFF_NO_FEATURES"}
_FASTA_WITHOUT_ENTRIES = {"code": "NO_RECORDS", "reason": "FASTA_NO_ENTRIES"}
_SEQUENCE_UNDEFINED = {"code": "INPUT_UNREADABLE", "reason": "SEQUENCE_UNDEFINED"}
# A FASTA given without GFF3 is drawn as nucleotide records (D-17): each record has
# a sequence, and the sequence holds only IUPAC nucleotide codes, in either case.
_IUPAC_NUCLEOTIDE_CODES = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
# A search query may also be a protein (PD-OI-117): NCBI reads every letter as an
# amino acid code (B J O U X Z included) and "*" as a stop.
_PROTEIN_CODES = frozenset(string.ascii_letters + "*")
SequenceMolecule = Literal["nucleotide", "protein"]
# Sequences are checked this many letters at a time, never copied whole.
_VALIDATION_PIECE = 1 << 20
_FLAT_FILE_KINDS = "GenBank, DDBJ or EMBL flat file"
# What a file without a record is, by the reason of its NO_RECORDS diagnostic.
_ABSENCE_MESSAGES = {
    "EMPTY_FILE": "{path} is empty: it has no " + _FLAT_FILE_KINDS + " record.",
    "FASTA_IN_GENBANK": (
        "{path} is a FASTA file, not a " + _FLAT_FILE_KINDS + ". "
        "Give it as FASTA (--fasta), with its GFF3 file (--gff) when it has one."
    ),
    "NOT_GENBANK": (
        "{path} is not a " + _FLAT_FILE_KINDS + ": no line starts a record, which is "
        "'LOCUS' and seven spaces (GenBank, DDBJ) or 'ID' and three spaces (EMBL)."
    ),
}


def _parser_reason(error: BaseException) -> str:
    """The parser's own reason, carried in the raised error (OV-1241)."""

    return str(error).strip().rstrip(".")


def genbank_absence_diagnostic(path: str) -> dict[str, str]:
    """The NO_RECORDS diagnostic of a GenBank/DDBJ/EMBL file that has no record (UJ-07).

    Its first non-blank line tells an empty file, a FASTA file, and other text apart.
    """

    # A byte order mark is not text, also on a line of its own.
    lines = (
        line.removeprefix(b"\xef\xbb\xbf") if index == 0 else line
        for index, line in enumerate(read_input_bytes(path).splitlines())
    )
    first = next((line.strip() for line in lines if line.strip()), b"")
    reason = "EMPTY_FILE" if not first else "FASTA_IN_GENBANK" if first.startswith(b">") else "NOT_GENBANK"
    return {"code": "NO_RECORDS", "reason": reason}


def read_flat_file_records(path: str, *, first_only: bool = False) -> list[SeqRecord]:
    """Every record of one GenBank, DDBJ or EMBL flat file, plain or gzip, in file order: the one flat-file parse.
    ``first_only`` reads the first record and leaves the rest of the file unread.

    The first non-blank line names the format (an ``ID`` line is EMBL). The text
    reaches Biopython through :func:`gbdraw.io.input_text.flat_file_stream`, and each
    feature type gets its INSDC spelling. A file the parser refuses raises
    ``ParseError`` (``INPUT_UNREADABLE``); a file without a record raises
    ``ValidationError`` (``NO_RECORDS``) that says what the file is instead.
    """

    fmt: FlatFileFormat = "genbank"

    def parse(lines: Iterator[str], file_format: FlatFileFormat) -> list[SeqRecord]:
        nonlocal fmt
        fmt = file_format
        return list(islice(SeqIO.parse(flat_file_stream(lines, fmt), fmt), 1 if first_only else None))

    try:  # A file that is not text stops at its first line: UnicodeDecodeError is a ValueError.
        records = read_flat_records(path, parse, flat_file_format)
    except (ValueError, *GZIP_ERRORS) as e:
        name = FLAT_FILE_FORMAT_NAMES[fmt]
        raise ParseError(
            f"The {name} parser stopped: {_parser_reason(e)}. "
            f"Error parsing {name} file {path}. It may be corrupt or in the wrong format.",
            diagnostic=_UNREADABLE_FLAT_FILE,
        ) from e
    if not records:
        diagnostic = genbank_absence_diagnostic(path)
        raise ValidationError(_ABSENCE_MESSAGES[diagnostic["reason"]].format(path=path), diagnostic=diagnostic)
    for record in records:
        canonicalize_feature_types(record.features)
    return records


def require_record_sequence(record: SeqRecord, source: str) -> None:
    """Stop at a record that has no sequence to draw (a CONTIG line, or no ORIGIN sequence)."""

    if getattr(record.seq, "defined", True):
        return
    raise ValidationError(
        f"Record {record.id} in {os.path.basename(source)} has no sequence: the file gives a CONTIG line "
        "or no ORIGIN sequence in its place. Download the record with its sequence (NCBI: 'GenBank (full)', "
        "rettype=gbwithparts), or give the sequence as GFF3 + FASTA (--gff with --fasta).",
        diagnostic=_SEQUENCE_UNDEFINED,
    )


def _attach_source_annotations(record: SeqRecord, source_file: str) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations["gbdraw_source_file"] = source_file
    record.annotations["gbdraw_source_basename"] = os.path.basename(source_file)
    # Coordinate mapping from current left-to-right index (0-based) to absolute
    # genomic coordinate: coord = base + step * index.
    # Defaults for untouched records: 1..len (step=+1).
    if "gbdraw_coord_base" not in record.annotations:
        record.annotations["gbdraw_coord_base"] = 1
    if "gbdraw_coord_step" not in record.annotations:
        record.annotations["gbdraw_coord_step"] = 1


def _load_record_files(
    paths: List[str],
    read: Callable[[str], list[SeqRecord]],
    record_selectors: list[str] | None,
    reverse_flags: list[bool] | None,
    *,
    kind: str,
) -> list[SeqRecord]:
    """Read each file with ``read``, then apply its record selector and reverse flag, in file order.

    The one per-file loop of the inputs whose files hold whole records (GenBank,
    a FASTA given without GFF3). A repeated record ID is warned about, not refused.
    """

    record_list: list[SeqRecord] = []
    id_list: list[str] = []
    logger.info(f"Loading {kind} file(s)...")
    for file_idx, path in enumerate(paths):
        if not os.path.isfile(path):
            raise InputFileError(f"File does not exist or is not accessible: {path}")
        try:
            logger.info(f"Loading {kind} file {path}")
            records_list = read(path)
            selector_raw = (
                record_selectors[file_idx]
                if record_selectors and file_idx < len(record_selectors)
                else None
            )
            selector = parse_record_selector(selector_raw)
            if selector is not None:
                records_list = select_record(records_list, selector, log=logger)
            reverse_flag = (
                reverse_flags[file_idx]
                if reverse_flags and file_idx < len(reverse_flags)
                else False
            )
            records_list = reverse_records(records_list, reverse_flag, log=logger)
            logger.info("Importing all entries...")
            for record in records_list:
                if record.id in id_list:
                    logger.warning(
                        f"Record {record.id} seems to have been already loaded. Check for duplicates!"
                    )
                _attach_source_annotations(record, path)
                record_list.append(record)
                id_list.append(record.id)
        except GbdrawError:
            raise
        except ValueError as e:
            if str(e).startswith("Record selector"):
                raise ValidationError(str(e)) from e
            raise
        except Exception as e:  # A more generic catch-all for other unexpected issues
            raise ParseError(
                f"Unexpected error while processing {path}: {e}"
            ) from e
    logger.info(f"... finished loading {kind} file(s)")
    logger.info(f"Number of sequences loaded to gbdraw: {len(id_list)}")
    if len(id_list) < 1:
        raise ValidationError(
            f"No valid {kind} records were loaded. Please check your input files."
        )
    return record_list


def load_gbks(
    gbk_list: List[str],
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    """Parse GenBank files without applying diagram-mode policy."""

    return _load_record_files(
        gbk_list, read_flat_file_records, record_selectors, reverse_flags, kind="GenBank"
    )


def _fasta_record_diagnostic(
    reason: str, record: SeqRecord, position: int, letter: str | None = None
) -> dict[str, object]:
    """The diagnostic of a FASTA-only record gbdraw cannot draw: its ID, its 1-based position, and the letter.

    The Web error adapter keeps the ID when it can be shown, else the position.
    """
    diagnostic: dict[str, object] = {
        "code": "INPUT_INVALID", "reason": reason, "recordName": str(record.id), "recordPosition": position,
    }
    if letter is not None:
        diagnostic["codepoint"] = ord(letter)
    return diagnostic


def require_nucleotide_records(records: list[SeqRecord], path: str) -> list[SeqRecord]:
    """Stop at the first record of a FASTA given without GFF3 that gbdraw cannot draw (D-17).

    A record needs a sequence, and its sequence may hold only IUPAC nucleotide
    codes (A C G T U R Y S W K M B D H V N, either case). The diagnostic names the
    record (its ID and its 1-based position) and, for a letter that is not a
    nucleotide code, that letter.
    """

    name = os.path.basename(path)
    for position, record in enumerate(records, start=1):
        _require_sequence(record, position, name)
        letter = _first_letter_outside(record.seq, _IUPAC_NUCLEOTIDE_CODES)
        if letter is not None:
            raise ValidationError(
                f"FASTA record {record.id or f'#{position}'} in {name} holds {letter!r}, which is not a nucleotide code: "
                "gbdraw draws nucleotide records (IUPAC codes A C G T U R Y S W K M B D H V N).",
                diagnostic=_fasta_record_diagnostic("NOT_NUCLEOTIDE", record, position, letter),
            )
    return records


def _require_sequence(record: SeqRecord, position: int, name: str) -> None:
    """Refuse the ``position``-th FASTA record of file ``name`` when it has no sequence."""

    if not len(record.seq):
        raise ValidationError(
            f"FASTA record {record.id or f'#{position}'} in {name} has no sequence. Give each FASTA record its sequence.",
            diagnostic=_fasta_record_diagnostic("EMPTY_SEQUENCE", record, position),
        )


def _first_letter_outside(sequence: Seq, codes: frozenset[str]) -> str | None:
    """The first letter of ``sequence`` that is not in ``codes``, read a piece at a time, or ``None``."""

    for start in range(0, len(sequence), _VALIDATION_PIECE):
        piece = str(sequence[start:start + _VALIDATION_PIECE])
        if set(piece) - codes:
            return next(char for char in piece if char not in codes)
    return None


def load_query_fasta(path: str) -> list[tuple[SeqRecord, SequenceMolecule]]:
    """The entries of a FASTA file of search queries, each with what it is.

    An entry is nucleotide when it holds IUPAC nucleotide codes only (the rule of
    :func:`require_nucleotide_records`), else protein when it holds amino acid
    codes and ``*`` only. Refused as :func:`load_fasta` refuses a file: one without
    an entry, an entry without sequence; and an entry holding a character that is
    neither a nucleotide nor an amino acid code (``NOT_SEQUENCE``, with that character).
    """

    if not os.path.isfile(path):
        raise InputFileError(f"File does not exist or is not accessible: {path}")
    name = os.path.basename(path)
    queries: list[tuple[SeqRecord, SequenceMolecule]] = []
    for position, record in enumerate(_read_fasta_file(path), start=1):
        _require_sequence(record, position, name)
        molecule: SequenceMolecule = "nucleotide"
        if _first_letter_outside(record.seq, _IUPAC_NUCLEOTIDE_CODES) is not None:
            molecule = "protein"
            letter = _first_letter_outside(record.seq, _PROTEIN_CODES)
            if letter is not None:
                raise ValidationError(
                    f"FASTA record {record.id or f'#{position}'} in {name} holds {letter!r}, which is neither a nucleotide "
                    "nor an amino acid code. Give DNA, RNA or protein sequences (letters, and * for a stop).",
                    diagnostic=_fasta_record_diagnostic("NOT_SEQUENCE", record, position, letter),
                )
        queries.append((record, molecule))
    return queries


def _read_fasta_file(path: str) -> list[SeqRecord]:
    """The records of one FASTA file, a parser error refused as an unreadable input."""

    try:
        return _read_fasta_records(path)
    except GbdrawError:
        raise
    except (ValueError, *GZIP_ERRORS) as e:
        raise ParseError(
            f"Error parsing FASTA file {path}: {_parser_reason(e)}.", diagnostic=_UNREADABLE_INPUT
        ) from e


def _read_nucleotide_fasta(path: str) -> list[SeqRecord]:
    """The records of one FASTA file given without GFF3, refused unless each is a nucleotide record."""

    return require_nucleotide_records(_read_fasta_file(path), path)


def load_fasta(
    fasta_list: List[str],
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    """Parse FASTA files given without GFF3: each entry is a record with its sequence and no features.

    The one loader of a FASTA-only input (CLI, Python API, the Web's Worker
    helper). A FASTA file in a GFF3 + FASTA pair is read by
    :func:`parse_gff_fasta` instead, where only the entries the GFF3 names are
    records. A FASTA record makes no topology claim (``detected_topology`` gives
    ``"unknown"``).
    """

    return _load_record_files(
        fasta_list, _read_nucleotide_fasta, record_selectors, reverse_flags, kind="FASTA"
    )


def merge_gff_fasta_records(
    gff_records: List[SeqRecord], fasta_records: List[SeqRecord], *, fasta_file: str
) -> List[SeqRecord]:
    """Give each GFF3 record the sequence of the FASTA record with its ID.

    The FASTA IDs are the join key, so ``fasta_file`` must not repeat an ID
    that a GFF3 record names; a FASTA sequence without GFF3 rows is not a
    record, and its repeats are ignored as before.
    """
    id_counts = Counter(record.id for record in fasta_records)
    gff_ids = {gff_record.id for gff_record in gff_records}
    duplicated_ids = [
        record_id for record_id, count in id_counts.items() if count > 1 and record_id in gff_ids
    ]
    if duplicated_ids:
        message = (
            f"FASTA file {fasta_file} has more than one sequence with the ID of GFF3 "
            f"record(s) {', '.join(duplicated_ids)}. Give each FASTA sequence a unique ID."
        )
        raise ValidationError(message, diagnostic=_DUPLICATE_FASTA_IDS)
    merged_records: List[SeqRecord] = []
    fasta_dict: Dict[str, SeqRecord] = {record.id: record for record in fasta_records}
    fasta_order: Dict[str, int] = {
        record.id: index for index, record in enumerate(fasta_records)
    }
    for gff_record in gff_records:
        if gff_record.id in fasta_dict:
            try:
                fasta_record = fasta_dict[gff_record.id]
                # If FASTA record exists, record.seq is FASTA sequence
                gff_record.seq = fasta_record.seq
            except Exception as e:
                raise ParseError(
                    f"Failed to merge sequences for record {gff_record.id} with corresponding FASTA entry: {e}"
                ) from e
            merged_records.append(gff_record)
        else:
            # If no matching FASTA record, raise error and stop the process
            raise ValidationError(
                f"No matching FASTA record found for GFF record {gff_record.id}. Please ensure that all GFF records have corresponding FASTA entries."
            )
    return sorted(merged_records, key=lambda record: fasta_order[record.id])


def _gff3_feature_id(feature: SeqFeature) -> str | None:
    """Return the canonical GFF3 ID for a parsed feature, if present."""
    raw_id = feature.qualifiers.get("ID")
    if isinstance(raw_id, (list, tuple)):
        raw_id = raw_id[0] if raw_id else None
    feature_id = str(raw_id or "").strip()
    return feature_id or None


def _merge_gff3_qualifiers(features: list[SeqFeature]) -> dict:
    """Merge qualifiers while retaining each multipart CDS phase in part order."""
    merged: dict = {}
    for feature in features:
        for key, raw_values in feature.qualifiers.items():
            values = list(raw_values) if isinstance(raw_values, (list, tuple)) else [raw_values]
            target_values = merged.setdefault(key, [])
            for value in values:
                if key == "phase" or value not in target_values:
                    target_values.append(value)
    return merged


def _order_gff3_feature_group(features: list[SeqFeature]) -> list[SeqFeature]:
    """Return multipart rows in biological order when their strand is consistent."""
    strands = {
        part.strand
        for feature in features
        for part in feature.location.parts
    }
    if strands not in ({1}, {-1}):
        return features

    def location_key(feature: SeqFeature) -> tuple[int, int]:
        parts = feature.location.parts
        return min(int(part.start) for part in parts), max(int(part.end) for part in parts)

    return sorted(features, key=location_key, reverse=strands == {-1})


def _normalize_gff3_feature_list(features: list[SeqFeature]) -> tuple[list[SeqFeature], int]:
    """Collapse same-ID GFF3 rows into multipart Biopython features."""
    grouped_features: list[list[SeqFeature]] = []
    group_indexes: dict[tuple[str, str], int] = {}

    for feature in features:
        feature_id = _gff3_feature_id(feature)
        location = getattr(feature, "location", None)
        if feature_id is None or location is None:
            grouped_features.append([feature])
            continue

        key = (str(feature.type), feature_id)
        group_index = group_indexes.get(key)
        if group_index is None:
            group_indexes[key] = len(grouped_features)
            grouped_features.append([feature])
        else:
            grouped_features[group_index].append(feature)

    normalized: list[SeqFeature] = []
    merged_count = 0
    for feature_group in grouped_features:
        feature_group = _order_gff3_feature_group(feature_group)
        feature = feature_group[0]
        if len(feature_group) > 1:
            parts = [
                part
                for grouped_feature in feature_group
                for part in grouped_feature.location.parts
            ]
            feature.location = CompoundLocation(parts, operator="join")
            feature.qualifiers = _merge_gff3_qualifiers(feature_group)
            feature.id = _gff3_feature_id(feature) or feature.id
            feature.sub_features = [
                sub_feature
                for grouped_feature in feature_group
                for sub_feature in (getattr(grouped_feature, "sub_features", None) or [])
            ]
            merged_count += 1

        sub_features = getattr(feature, "sub_features", None) or []
        if sub_features:
            feature.sub_features, child_merge_count = _normalize_gff3_feature_list(sub_features)
            merged_count += child_merge_count
        normalized.append(feature)

    return normalized, merged_count


def _normalize_gff3_multipart_features(record: SeqRecord) -> SeqRecord:
    """Normalize repeated GFF3 IDs into ``CompoundLocation`` features in-place.

    GFF3 uses multiple rows with the same ID to represent a single discontinuous
    feature. BCBio-GFF retains those rows as separate ``SeqFeature`` objects, so
    normalize them before feature filtering and rendering.
    """
    record.features, merged_count = _normalize_gff3_feature_list(record.features)
    if merged_count:
        logger.info(
            "Normalized %d multipart GFF3 feature(s) for record %s.",
            merged_count,
            record.id,
        )
    return record


def scan_features_recursive(
    features: List[SeqFeature], feature_types_to_keep: AbstractSet[str] | None,
    source_indexes: dict[int, int] | None = None,
) -> List[SeqFeature]:
    """Flatten nested features in source order; ``None`` keeps every type."""
    filtered_list: list[SeqFeature] = []
    for feature in features:
        if feature_types_to_keep is None or feature.type in feature_types_to_keep:
            new_feature = SeqFeature(
                location=feature.location,
                type=feature.type,
                qualifiers=feature.qualifiers,
                id=feature.id,
            )
            if source_indexes is not None:
                _copy_source_feature_identity(feature, new_feature,
                    index=source_indexes[id(feature)], coord_base=1, coord_step=1, period=None)
            filtered_list.append(new_feature)

        if feature.sub_features:
            # Recursively call and extend the list with the results
            filtered_list.extend(
                scan_features_recursive(feature.sub_features, feature_types_to_keep, source_indexes)
            )

    return filtered_list


def filter_features_by_type(
    record: SeqRecord, feature_types_to_keep: AbstractSet[str] | None, *, source_indexes: dict[int, int] | None = None,
) -> SeqRecord:
    new_record = SeqRecord(
        seq=record.seq,
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs,
        annotations=record.annotations,
    )

    new_record.features = scan_features_recursive(record.features, feature_types_to_keep, source_indexes)
    return new_record


def feature_types_left_out(
    features: List[SeqFeature], feature_types_to_keep: AbstractSet[str] | None,
) -> tuple[tuple[str, int], ...]:
    """The features, nested ones included, that ``filter_features_by_type`` leaves out, counted by type."""

    if feature_types_to_keep is None:
        return ()
    counts: Counter[str] = Counter()
    pending = list(features)
    while pending:
        feature = pending.pop()
        if feature.type not in feature_types_to_keep:
            counts[str(feature.type)] += 1
        pending.extend(getattr(feature, "sub_features", None) or ())
    return tuple(sorted(counts.items()))


@contextmanager
def _gff_fasta_errors(gff_file: str, fasta_file: str):
    """Report a failure while reading a GFF3/FASTA pair the way every caller expects."""
    try:
        yield
    except (ValueError, *GZIP_ERRORS) as e:
        if str(e).startswith("Record selector"):
            raise ValidationError(str(e)) from e
        raise ParseError(
            f"Error parsing GFF3/FASTA files ({gff_file}, {fasta_file}): {_parser_reason(e)}.",
            diagnostic=_UNREADABLE_INPUT,
        ) from e
    except GbdrawError:
        raise
    except AssertionError as e:
        # BCBio-GFF asserts when it cannot split an attribute column into key=value pairs.
        raise ParseError(
            f"The GFF3 parser (BCBio-GFF) could not split an attribute column of {gff_file} "
            f"into key=value pairs: {e}. Separate the attributes with ';' and write each as key=value.",
            diagnostic=_UNREADABLE_INPUT,
        ) from e
    except Exception as e:
        raise ParseError(
            f"Unexpected error while processing {gff_file} or {fasta_file}: {e}"
        ) from e


def _read_fasta_records(fasta_file: str) -> list[SeqRecord]:
    """Read the records of a FASTA file the same way under every Biopython version.

    Text before the first ">" line is skipped, as Biopython 1.85 (the Web's) does;
    Biopython 1.88 refuses it. A file without a ">" line has no entry (OV-1257).
    """

    records = read_fasta_records(fasta_file, _fasta_entries)
    if records is None:
        raise ValidationError(f"FASTA file {fasta_file} has no sequence entry.", diagnostic=_FASTA_WITHOUT_ENTRIES)
    return records


def _fasta_entries(lines: Iterator[str]) -> list[SeqRecord]:
    """The entries of a FASTA file's header lines (the first line starts one); each sequence comes from its file."""

    return list(SeqIO.parse(LineStream(lines), "fasta"))


def read_sequence_file_records(path: str, fmt: str) -> list[SeqRecord]:
    """Every record of one sequence file, in file order, as the render reads it.

    ``genbank`` reads with :func:`load_gbks`; ``fasta`` with the FASTA reader of
    :func:`parse_gff_fasta`, so it keeps every entry, also one the GFF3 file does
    not name. The Circular similarity ring database (whole source file, D-25) and
    the Web's whole-file reader (``read_all_sequence_records``) read files here.
    """
    if fmt == "genbank":
        return load_gbks([path])
    if fmt == "fasta":
        return _read_fasta_records(path)
    raise ValueError(f"Unsupported format: {fmt}")


def record_sequence_text(record: SeqRecord) -> str | None:
    """A record's sequence in upper case, or ``None`` when it has none (a CON record)."""
    if not getattr(record.seq, "defined", True):
        return None
    return str(record.seq).upper()


def parse_gff_fasta(
    gff_file: str,
    fasta_file: str,
    *,
    source_feature_catalogs: list[SourceFeatureCatalog] | None = None,
) -> list[SeqRecord]:
    """Parse one GFF3/FASTA pair once, whatever feature types a diagram draws.

    Records follow the FASTA order and keep every feature, with Parent-linked
    children nested in ``sub_features``. ``filter_features_by_type`` then selects
    types from them, as often as needed, without parsing the file again.
    ``source_feature_catalogs`` gains one catalog per record, built before any filter.
    """
    for path in (gff_file, fasta_file):
        if not os.path.isfile(path):
            raise InputFileError(f"File does not exist or is not accessible: {path}")
    with _gff_fasta_errors(gff_file, fasta_file):
        logger.info("Loading GFF3 file {}".format(gff_file))
        gff_records = read_gff3_text(gff_file, lambda lines: list(GFF.parse(gff3_stream(lines))))
        for record in gff_records:
            canonicalize_feature_types(record.features)
            _normalize_gff3_multipart_features(record)
        if not gff_records:
            raise ValidationError(f"GFF3 file {gff_file} has no feature lines.", diagnostic=_GFF_WITHOUT_FEATURES)
        catalogs = {
            id(record): capture_source_features(record) for record in gff_records
        } if source_feature_catalogs is not None else {}
        logger.info("Loading FASTA file {}".format(fasta_file))
        fasta_records = _read_fasta_records(fasta_file)
        merged_records = merge_gff_fasta_records(gff_records, fasta_records, fasta_file=fasta_file)
        if source_feature_catalogs is not None:
            source_feature_catalogs.extend(catalogs[id(record)] for record in merged_records)
        for record in merged_records:
            _attach_source_annotations(record, gff_file)
        return merged_records


def load_gff_fasta(
    gff_list: List[str],
    fasta_list: List[str],
    selected_features_set=None,
    keep_all_features: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
    *,
    source_feature_catalogs: list[SourceFeatureCatalog] | None = None,
) -> list[SeqRecord]:
    """Parse paired GFF3/FASTA files without applying diagram-mode policy."""

    record_list: list[SeqRecord] = []
    id_list: list[str] = []
    logger.info("Loading GFF3/FASTA file(s)...")

    if len(gff_list) != len(fasta_list):
        raise ValidationError("Number of GFF3 files does not match number of FASTA files.")

    for file_idx, (gff_file, fasta_file) in enumerate(zip(gff_list, fasta_list)):
        file_catalogs: list[SourceFeatureCatalog] | None = (
            [] if source_feature_catalogs is not None else None
        )
        parsed_records = parse_gff_fasta(gff_file, fasta_file, source_feature_catalogs=file_catalogs)
        with _gff_fasta_errors(gff_file, fasta_file):
            # The renderer reads record.features only, so Parent-linked children
            # (CDS under gene) are flattened whether or not types are filtered.
            feature_types_to_keep = (
                None if keep_all_features or selected_features_set is None
                else set(selected_features_set)
            )
            merged_records = [
                filter_features_by_type(record, feature_types_to_keep,
                    source_indexes=_feature_source_index_map(record.features)
                    if file_catalogs is not None else None)
                for record in parsed_records
            ]
            catalog_by_record = {id(record): catalog for record, catalog in zip(merged_records, file_catalogs)} if file_catalogs is not None else {}

            selector_raw = (
                record_selectors[file_idx]
                if record_selectors and file_idx < len(record_selectors)
                else None
            )
            selector = parse_record_selector(selector_raw)
            if selector is not None:
                merged_records = select_record(merged_records, selector, log=logger)
            reverse_flag = (
                reverse_flags[file_idx]
                if reverse_flags and file_idx < len(reverse_flags)
                else False
            )
            selected_catalogs = tuple(catalog_by_record[id(record)] for record in merged_records) if source_feature_catalogs is not None else ()
            merged_records = reverse_records(merged_records, reverse_flag, log=logger)

            for record_index, record in enumerate(merged_records):
                if record.id in id_list:
                    logger.warning(
                        f"Record {record.id} seems to have been already loaded. Check for duplicates!"
                    )
                if source_feature_catalogs is not None:
                    source_feature_catalogs.append(selected_catalogs[record_index])
                record_list.append(record)
                id_list.append(record.id)

    logger.info("... finished loading GFF3 and FASTA files")
    logger.info(f"Number of sequences loaded to gbdraw: {len(record_list)}")
    if len(record_list) < 1:
        raise ValidationError(
            "No valid records were loaded after merging GFF3 and FASTA files. Please check your input files."
        )
    return record_list


def source_genome_length(kind: str, paths: Sequence[str], loaded_length: int) -> int:
    """The length a source counts in its genome (OV-3143): every sequence of the source, uncropped, drawn or not.

    A GenBank or FASTA file, or an in-memory record, counts the records its loader read (``loaded_length``, their
    uncropped lengths summed); a GFF3/FASTA pair (``kind`` ``gff_fasta``, ``paths`` the GFF3 and FASTA files)
    counts every entry of its FASTA file, also one the GFF3 file does not name, as the Circular similarity
    rings search it. The one owner of a source's genome length.
    """

    if kind != "gff_fasta":
        return loaded_length
    return sum(len(record) for record in read_sequence_file_records(str(paths[-1]), "fasta"))


__all__ = [
    "source_genome_length",
    "load_gbks",
    "load_fasta",
    "load_query_fasta",
    "require_nucleotide_records",
    "SequenceMolecule",
    "read_flat_file_records",
    "require_record_sequence",
    "parse_gff_fasta",
    "read_sequence_file_records",
    "record_sequence_text",
    "load_gff_fasta",
    "merge_gff_fasta_records",
    "scan_features_recursive",
    "feature_types_left_out",
    "filter_features_by_type",
]


