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
from gbdraw.io.comparison_sequences import read_comparison_sequence_file


def read_comparison_sequence_json(path: str) -> str:
    """Return the LOSAT query FASTA, format, record label and record IDs as JSON.

    ``recordLabel`` is the CLI default ring label when the file names one
    (GenBank / DDBJ: the first record's DEFINITION, else its organism), else
    ``null``: the Worker stages the file under a fixed name, so the Web applies
    the file-name rule of ``comparison_file_stem`` itself. Raises
    ``INPUT_UNREADABLE`` (``SEQUENCE_MISSING`` for a file or record without
    sequence) like the CLI ring search.
    """

    comparison = read_comparison_sequence_file(path)
    fasta = comparison_query_fasta(comparison, ordinal=1)
    return json.dumps(
        {
            "fasta": fasta,
            "format": comparison.format,
            "recordLabel": comparison.record_label,
            "recordIds": [str(record.id) for record in comparison.records],
        },
        separators=(",", ":"),
    )


__all__ = ["read_comparison_sequence_json"]
