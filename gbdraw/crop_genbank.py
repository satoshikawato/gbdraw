#!/usr/bin/env python
# coding: utf-8
"""Crop a GenBank record to a region and shift its features; also runs as a script."""

import argparse
from io import StringIO
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SeqFeature, FeatureLocation, CompoundLocation,
    BeforePosition, AfterPosition, ExactPosition
)
from gbdraw.core.numeric import argparse_integer
from gbdraw.core.record_metadata import (
    _copy_source_feature_identity,
    _feature_source_index_map,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.genome import read_flat_file_records
from gbdraw.io.record_sequences import sequence_view

def _get_args():
    parser = argparse.ArgumentParser(description='Crop genbank file. ')
    parser.add_argument(
        '-i',
        '--input',
        help='GenBank/DDBJ flat file (required)',
        type=str,
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        "--output",
        metavar="FILE",
        help="output GenBank file",
        required=True)
    parser.add_argument(
        "-s",
        "--start",
        type=argparse_integer,
        help="start position (1-based, inclusive)",
        required=True)
    parser.add_argument(
        "-e",
        "--end",
        type=argparse_integer,
        help="end position (1-based, inclusive)",
        required=True)
    args = parser.parse_args()
    return args


def gbk_to_seqrecord(in_gbk):
    """The first record of a GenBank file, read as the drawing reads it (``read_flat_file_records``): its sequence
    is read from the file on demand."""
    return read_flat_file_records(in_gbk, first_only=True)[0]


def check_start_end_coords(record_len: int, start: int, end: int, record_id: str) -> tuple[int, int]:
    """
    Converts 1-based inclusive start/end to 0-based exclusive slice indices.

    The order is checked before the end is clamped to the record length, so a
    region that starts after the record ends is reported as such (CI-01). It
    reads only the record's length (``record_id`` names it in a refusal), so
    the drawing's crop and the Web export check of a cropped record's drawn
    bases share it (D-70).
    """
    full_record = start == 1 and end == record_len
    if start > end and not full_record:
        raise ValidationError(
            f"Region start ({start}) must not exceed the region end ({end}).",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "ORDER"},
        )
    # 1-based inclusive start/end -> 0-based slice; the end is clamped.
    start_0 = max(0, start - 1)
    end_0 = min(end, record_len)
    if start_0 >= end_0 and not full_record:
        raise ValidationError(
            f"Region {start}..{end} starts after the end of record {record_id} "
            f"({record_len} bp). Choose a region within 1..{record_len}.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
        )
    return start_0, end_0

def _part_in_window(part, crop_start_0, crop_end_0):
    """Whether a location part keeps any base inside the crop window."""
    return part.start < crop_end_0 and part.end > crop_start_0


def _crop_and_shift_location(loc_before, loc_current, loc_next, crop_start_0, crop_end_0):
    """Crop one location part to the window and shift it to window coordinates.

    An end is marked partial (``<`` low, ``>`` high) when the part is cut by the
    window, when the neighbouring part on that side (in biological order) is
    dropped, or when the source already marks it partial. Biopython position
    arithmetic returns plain integers, so the source markers are read before the
    shift, not carried by it.
    """
    new_seq_len = crop_end_0 - crop_start_0
    dropped_before = loc_before is not None and not _part_in_window(loc_before, crop_start_0, crop_end_0)
    dropped_next = loc_next is not None and not _part_in_window(loc_next, crop_start_0, crop_end_0)
    # Parts are in biological order: on the minus strand the part before lies above this one.
    dropped_low, dropped_high = (
        (dropped_next, dropped_before) if loc_current.strand == -1 else (dropped_before, dropped_next)
    )
    open_low = dropped_low or loc_current.start < crop_start_0 or isinstance(loc_current.start, BeforePosition)
    open_high = dropped_high or loc_current.end > crop_end_0 or isinstance(loc_current.end, AfterPosition)
    start = max(0, int(loc_current.start) - crop_start_0)
    end = min(new_seq_len, int(loc_current.end) - crop_start_0)
    return FeatureLocation(
        BeforePosition(start) if open_low else ExactPosition(start),
        AfterPosition(end) if open_high else ExactPosition(end),
        strand=loc_current.strand,
    )

def crop_and_shift_features(
    original_features,
    start_0,
    end_0,
    *,
    coord_base=1,
    coord_step=1,
    period=None,
    _source_indexes=None,
):
    """
    Iterates over old features, crops/shifts them relative to the
    crop window (start_0, end_0), and returns a new list of features.
    ``coord_base``, ``coord_step`` and ``period`` are the coordinate map of
    the record the features belong to.
    """
    new_features_list = []
    source_indexes = (
        _feature_source_index_map(original_features)
        if _source_indexes is None
        else _source_indexes
    )
    for old_feature in original_features:

        if not (old_feature.location.start < end_0 and old_feature.location.end > start_0):
            continue

        old_loc = old_feature.location
        new_location = None

        if isinstance(old_loc, CompoundLocation):
            new_parts = []
            num_parts = len(old_loc.parts)
            for n in range(num_parts):
                part = old_loc.parts[n]

                if not _part_in_window(part, start_0, end_0):
                    continue

                part_before = old_loc.parts[n - 1] if n > 0 else None
                part_next = old_loc.parts[n + 1] if n < (num_parts - 1) else None

                new_part = _crop_and_shift_location(part_before, part, part_next, start_0, end_0)
                new_parts.append(new_part)

            if not new_parts:
                continue

            if len(new_parts) == 1:
                new_location = new_parts[0]
            else:
                new_location = CompoundLocation(new_parts, operator=old_loc.operator)

        else:
            part_before = None
            part_next = None
            new_location = _crop_and_shift_location(part_before, old_loc, part_next, start_0, end_0)

        new_feature = SeqFeature(
            location=new_location,
            type=old_feature.type,
            qualifiers=old_feature.qualifiers,
            id=old_feature.id
        )
        _copy_source_feature_identity(
            old_feature,
            new_feature,
            index=source_indexes[id(old_feature)],
            coord_base=coord_base,
            coord_step=coord_step,
            period=period,
        )
        sub_features = getattr(old_feature, "sub_features", None)
        if sub_features is not None:
            new_feature.sub_features = crop_and_shift_features(
                sub_features,
                start_0,
                end_0,
                coord_base=coord_base,
                coord_step=coord_step,
                period=period,
                _source_indexes=source_indexes,
            )
        new_features_list.append(new_feature)

    return new_features_list


def write_cropped_genbank(record_to_write, out_filepath, original_accession, start_1based, end_1based):
    """
    Writes the SeqRecord to a file, modifying the ACCESSION line
    to include the cropped REGION.
    """
    temp_handle = StringIO()
    SeqIO.write(record_to_write, temp_handle, "genbank")
    content = temp_handle.getvalue()

    # Add REGION to ACCESSION line
    content = re.sub(
        rf'ACCESSION\s+{original_accession}',
        f'ACCESSION   {original_accession} REGION: {start_1based}..{end_1based}',
        content
    )

    with open(out_filepath, 'w', encoding='utf-8') as f:
        f.write(content)

def main():
    # 1. Setup
    args = _get_args()
    in_gbk = args.input
    out_gbk = args.output
    start_1 = args.start
    end_1 = args.end

    # 2. Load and Validate
    record = gbk_to_seqrecord(in_gbk)
    start_0, end_0 = check_start_end_coords(len(record.seq), start_1, end_1, str(record.id))

    # 3. Crop Sequence and Create New Record
    new_record = SeqRecord(
        sequence_view(record.seq, start_0, end_0),
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs.copy(),
        annotations=record.annotations
    )

    # 4. Process Features (Extracted Function)
    new_record.features = crop_and_shift_features(record.features, start_0, end_0)

    # 5. Write Output (Extracted Function)
    original_accession = record.annotations['accessions'][0]
    write_cropped_genbank(new_record, out_gbk, original_accession, start_1, end_1)

if __name__ == "__main__":
    main()
