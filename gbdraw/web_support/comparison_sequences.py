#!/usr/bin/env python
# coding: utf-8

"""Web worker helper for Circular ring comparison files (design D12).

The Web has no sequence reader of its own: a ring's comparison file (FASTA,
GenBank or DDBJ) is read here with the one reader in ``gbdraw.io`` and turned
into the LOSAT query FASTA that the CLI ring search hashes, so the raw cache
keys of the CLI and the Web are equal for every format and FASTA layout.

Before its ring searches run, a Generate also sends the reference files and
each ring's files, and Python answers with each ring search's search space
(``plan_ring_scales_json``) from the owner the CLI ring search reads
(``ring_search_spaces``), so both pass one value (PD-OI-018 revision 6).
"""

from __future__ import annotations

import json
from typing import Literal, TypedDict

from gbdraw.comparisons.circular_losat import (
    comparison_genome_files,
    comparison_query_fasta,
    ring_genome_file,
    ring_search_spaces,
)
from gbdraw.comparisons.losat_jobs import PER_FILE_WARNING_TEXT
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import comparison_genome, read_comparison_sequence_file
from gbdraw.io.genome import read_sequence_file_records
from gbdraw.web_support.linear_units import GenomeScaleWarningParams


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
    ``INPUT_INVALID`` / ``UNIQUE_IDS`` when two files share a record ID.
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


class RingScaleWarning(TypedDict):
    """A ring searched per file, its genome pair's search space past what LOSAT reads, and the warning the command
    line logs for it (Owner decision 2026-11-10): ``text`` is its English text with ``{name}`` parameters
    (:data:`gbdraw.comparisons.losat_jobs.PER_FILE_WARNING_TEXT`, a catalog key), which the page shows with
    ``tm(text, params)``."""

    ring: int
    text: str
    params: GenomeScaleWarningParams


class RingGenomeScales(TypedDict):
    """The E-value scale of each ring search of a Circular drawing, in ring order: ``scales[k]`` is the
    ``-searchsp`` ring ``k``'s search passes, a decimal string, or null for LOSAT's own (PD-OI-018 revision 6);
    ``warnings`` lists the rings searched per file (:class:`RingScaleWarning`)."""

    scales: list[str | None]
    warnings: list[RingScaleWarning]


RingReferenceFormat = Literal["genbank", "fasta"]


def plan_ring_scales(
    program: str,
    reference_format: RingReferenceFormat,
    reference_paths: list[str],
    reference_names: list[str],
    ring_paths: list[list[str]],
    ring_names: list[list[str]],
) -> RingGenomeScales:
    """The search space of each ``program`` (``losatn`` or ``tlosatx``) ring search of a Circular drawing.

    ``reference_paths`` are the reference's files, one per source card, read whole in ``reference_format`` (a GFF3
    and FASTA pair gives its FASTA file, every entry of which the CLI ring database reads); ``ring_paths[k]`` are
    ring ``k``'s files. ``reference_names`` and ``ring_names`` give the names a warning gives them, as the user
    named the files (the Worker stages them under fixed names). The values are those of the CLI's ring searches of
    the same files (``ring_search_spaces``).
    """

    if program not in ("losatn", "tlosatx") or reference_format not in ("genbank", "fasta") or (
        len(reference_names) != len(reference_paths) or len(ring_names) != len(ring_paths)
        or any(len(names) != len(paths) for names, paths in zip(ring_names, ring_paths))
    ):
        raise ValidationError(
            "Ring scales need a nucleotide LOSAT program, a GenBank or FASTA reference and a name for every file.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    reference = [
        ring_genome_file(name, reference_format, (path,), read_sequence_file_records(path, reference_format))
        for path, name in zip(reference_paths, reference_names)
    ]
    rings = [
        comparison_genome_files(comparison_genome(tuple(read_comparison_sequence_file(path) for path in paths)), names)
        for paths, names in zip(ring_paths, ring_names)
    ]
    scales = ring_search_spaces(program, reference, rings)
    return {
        "scales": [None if space is None else str(space) for space in scales.spaces],
        "warnings": [
            {"ring": ring, "text": PER_FILE_WARNING_TEXT, "params": {
                "first": found.first, "second": found.second, "space": found.space, "limit": found.limit,
            }}
            for ring, found in enumerate(scales.warnings) if found is not None
        ],
    }


def plan_ring_scales_json(
    program: str,
    reference_format: str,
    reference_paths_json: str,
    reference_names_json: str,
    ring_paths_json: str,
    ring_names_json: str,
) -> str:
    """JSON boundary of the diagram Worker."""

    fmt = str(reference_format)
    return json.dumps(plan_ring_scales(
        str(program),
        "genbank" if fmt == "genbank" else "fasta" if fmt == "fasta" else _bad_format(fmt),
        json.loads(str(reference_paths_json)), json.loads(str(reference_names_json)),
        json.loads(str(ring_paths_json)), json.loads(str(ring_names_json)),
    ))


def _bad_format(fmt: str) -> RingReferenceFormat:
    raise ValidationError(f"Unsupported ring reference format {fmt!r}.", diagnostic={"code": "HELPER_PROTOCOL"})


__all__ = [
    "plan_ring_scales",
    "plan_ring_scales_json",
    "read_comparison_sequence_json",
    "RingGenomeScales",
    "RingScaleWarning",
]
