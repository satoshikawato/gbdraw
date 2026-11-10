#!/usr/bin/env python
# coding: utf-8

"""Web worker helper for Circular ring comparison files (design D12).

The Web has no sequence reader of its own: a ring's comparison file (FASTA,
GenBank or DDBJ) is read here with the one reader in ``gbdraw.io`` and turned
into the LOSAT query FASTA that the CLI ring search hashes, so the raw cache
keys of the CLI and the Web are equal for every format and FASTA layout.

Before its ring searches run, a Generate asks for each ring search's search space
(``plan_ring_scales_json``), which Python answers from the owner the CLI ring
search reads (``ring_search_spaces``), so both pass one value (PD-OI-018
revision 6). The reference is read as the next Generate's render reads it, from
its prepared parse; each ring's file lengths are those this reader answered.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from gbdraw.api.record_planning import prepared_parsed_source
from gbdraw.api.request_render import load_render_source_records
from gbdraw.api.requests import FastaInputSource, GenBankInputSource
from gbdraw.comparisons.circular_losat import (
    RingGenomeFile,
    comparison_genome_files,
    comparison_query_fasta,
    ring_search_spaces,
    ring_subject_sources,
)
from gbdraw.comparisons.losat_jobs import PER_FILE_WARNING_TEXT
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import comparison_genome, read_comparison_sequence_file
from gbdraw.session_request_codec import decode_record_source
from gbdraw.web_support.linear_units import RingGenomeFileLengths, RingGenomeInput, RingGenomeScales
from gbdraw.web_support.request_render import staged_source_names


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
    ``genomeFiles`` gives each file's lengths (:class:`RingGenomeFileLengths`), which
    the page sends back to the ring planner (``plan_ring_scales_json``).
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
            "genomeFiles": [
                RingGenomeFileLengths(length=file.length, recordLengths=list(file.record_lengths))
                for file in comparison_genome_files(comparison)
            ],
        },
        separators=(",", ":"),
    )



def plan_ring_scales(
    program: str,
    reference_sources: list[object],
    resource_paths: Mapping[str, str],
    resource_names: Mapping[str, object],
    rings: list[RingGenomeInput],
) -> RingGenomeScales:
    """The search space of each ``program`` (``losatn`` or ``tlosatx``) ring search of a Circular drawing.

    ``reference_sources`` are the record sources (``records[].source``) of the reference's source cards, in card
    order, under the resource ids the next Generate's request names them by, staged at ``resource_paths``: the
    reference is read from the parse that Generate renders (published for it), so it is parsed once.
    ``resource_names`` (``resourceManifest[].name``) names the files a warning names, as the input summary does.
    ``rings[k]`` is ring ``k``: its file names and the lengths the ring reader answered for them. The values are
    those of the CLI's ring searches of the same files (``ring_subject_sources``, ``ring_search_spaces``).
    """

    if program not in ("losatn", "tlosatx") or not isinstance(reference_sources, list) or not isinstance(rings, list):
        raise ValidationError(
            "Ring scales need a nucleotide LOSAT program, the reference sources and the rings.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    sources = [
        decode_record_source(source, path=f"referenceSources[{index}]", resource_paths=resource_paths)
        for index, source in enumerate(reference_sources)
    ]

    def file_records(source: GenBankInputSource | FastaInputSource) -> tuple[SeqRecord, ...]:
        parsed = prepared_parsed_source(source, load_render_source_records)
        return (parsed if parsed is not None else load_render_source_records(source)).records

    names = staged_source_names(resource_names, {key: Path(path) for key, path in resource_paths.items()}) or {}
    reference = [
        replace(found.genome, name=names.get(found.genome.name, found.genome.name))
        for found in ring_subject_sources(sources, file_records=file_records)
    ]
    scales = ring_search_spaces(program, reference, [_ring_genome_files(index, ring) for index, ring in enumerate(rings)])
    return {
        "scales": [None if space is None else str(space) for space in scales.spaces],
        "warnings": [
            {"ring": ring, "text": PER_FILE_WARNING_TEXT, "params": {
                "first": found.first, "second": found.second, "space": found.space, "limit": found.limit,
            }}
            for ring, found in enumerate(scales.warnings) if found is not None
        ],
    }


def _ring_genome_files(index: int, ring: object) -> list[RingGenomeFile]:
    """Ring ``index``'s genome files from the page's :class:`RingGenomeInput`: a name and the reader's lengths for
    every file, each a whole number, at least one record per file."""

    def count(value: object) -> bool:
        return isinstance(value, int) and not isinstance(value, bool) and value >= 0

    names = ring.get("names") if isinstance(ring, dict) else None
    files = ring.get("files") if isinstance(ring, dict) else None
    if not (
        isinstance(names, list) and isinstance(files, list) and len(names) == len(files) > 0
        and all(isinstance(name, str) for name in names)
        and all(
            isinstance(file, dict) and count(file.get("length")) and isinstance(file.get("recordLengths"), list)
            and len(file["recordLengths"]) > 0 and all(count(length) for length in file["recordLengths"])
            for file in files
        )
    ):
        raise ValidationError(f"Ring {index + 1} of the ring scales is not a list of named files with their lengths.",
                              diagnostic={"code": "HELPER_PROTOCOL"})
    return [
        RingGenomeFile(name, file["length"], tuple(file["recordLengths"])) for name, file in zip(names, files)
    ]


def plan_ring_scales_json(
    program: str, reference_sources_json: str, resource_paths_json: str, resource_names_json: str, rings_json: str
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_ring_scales(
        str(program), json.loads(str(reference_sources_json)), json.loads(str(resource_paths_json)),
        json.loads(str(resource_names_json)), json.loads(str(rings_json)),
    ))


__all__ = [
    "plan_ring_scales",
    "plan_ring_scales_json",
    "read_comparison_sequence_json",
]
