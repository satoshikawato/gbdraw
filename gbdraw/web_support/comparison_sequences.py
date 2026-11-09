#!/usr/bin/env python
# coding: utf-8

"""Web worker helper for Circular ring comparison files (design D12).

The Web has no sequence reader of its own: a ring's comparison file (FASTA,
GenBank or DDBJ) is read here with the one reader in ``gbdraw.io`` and turned
into the LOSAT query FASTA that the CLI ring search hashes, so the raw cache
keys of the CLI and the Web are equal for every format and FASTA layout.
"""

from __future__ import annotations

import json

from gbdraw.comparisons.circular_losat import comparison_query_fasta
from gbdraw.io.comparison_sequences import comparison_genome, read_comparison_sequence_file


def read_comparison_sequence_json(path: str, *more_paths: str) -> str:
    """Return a ring's LOSAT query FASTA, format, record label and record IDs as JSON.

    ``path`` and ``more_paths`` are the files of one comparison genome in the
    listed order (a chromosome and its plasmids in separate files), read by the
    one reader and joined by :func:`comparison_genome` as the CLI ring search
    joins them, so the query FASTA and its raw cache key are the CLI's.
    ``format`` and ``recordLabel`` are the first file's: ``recordLabel`` is the
    CLI default ring label when that file names one (GenBank / DDBJ: the first
    record's DEFINITION, else its organism), else ``null``: the Worker stages
    the files under fixed names, so the Web applies the file-name rule of
    ``comparison_file_stem`` itself. ``recordIds`` lists every record of every
    file in query order. Raises ``INPUT_UNREADABLE`` (``SEQUENCE_MISSING`` for
    a file or record without sequence) like the CLI ring search, and
    ``INPUT_INVALID`` / ``DUPLICATE_RECORD_ID`` when two files share a record ID.
    """

    comparison = comparison_genome(
        tuple(read_comparison_sequence_file(file_path) for file_path in (path, *more_paths))
    )
    fasta = comparison_query_fasta(comparison, ordinal=1)
    return json.dumps(
        {
            "fasta": fasta,
            "format": comparison.files[0].format,
            "recordLabel": comparison.record_label,
            "recordIds": [str(record.id) for record in comparison.records],
        },
        separators=(",", ":"),
    )


__all__ = ["read_comparison_sequence_json"]
