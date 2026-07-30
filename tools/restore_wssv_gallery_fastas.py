#!/usr/bin/env python3
"""Restore the exact FASTA inputs used by the WSSV Gallery session."""

from __future__ import annotations

import argparse
import base64
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import sys
from typing import Any, Mapping


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_io import (  # noqa: E402
    load_session,
    validate_session,
    write_session_json,
)


SESSION_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "WSSV_genome_comparison.gbdraw-session.json"
)
DEFAULT_MAG_DIR = REPO_ROOT / "tests" / "test_inputs" / "MAGs"
NCBI_NUCCORE_URL = "https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
NCBI_SRA_URL = "https://www.ncbi.nlm.nih.gov/sra/{accession}"


@dataclass(frozen=True)
class WssvFastaSource:
    """Pinned identity and provenance for one displayed WSSV comparison."""

    label: str
    record_id: str
    filename: str
    length: int
    sha256: str
    source: str
    accession: str | None
    source_url: str | None
    local_filename: str | None = None


def _nuccore(
    label: str,
    accession: str,
    length: int,
    sha256: str,
) -> WssvFastaSource:
    return WssvFastaSource(
        label=label,
        record_id=accession,
        filename=f"{label}.fasta",
        length=length,
        sha256=sha256,
        source="NCBI Nucleotide",
        accession=accession,
        source_url=NCBI_NUCCORE_URL.format(accession=accession),
    )


def _prepared(
    label: str,
    record_id: str,
    filename: str,
    length: int,
    sha256: str,
    *,
    sra_accession: str | None,
) -> WssvFastaSource:
    return WssvFastaSource(
        label=label,
        record_id=record_id,
        filename=filename,
        length=length,
        sha256=sha256,
        source=(
            f"Prepared assembly derived from NCBI SRA {sra_accession}"
            if sra_accession
            else "User-provided prepared assembly; no public accession recorded"
        ),
        accession=sra_accession,
        source_url=(
            NCBI_SRA_URL.format(accession=sra_accession)
            if sra_accession
            else None
        ),
        local_filename=filename,
    )


# These hashes are the raw UTF-8 FASTA hashes recorded by the session's
# nucleotide LOSAT cache. The twelve Nucleotide records were retrieved from
# NCBI EFetch; the other eight files were supplied under tests/test_inputs/MAGs.
WSSV_FASTA_SOURCES = (
    _nuccore(
        "CN01",
        "NC_003225.3",
        309286,
        "fa765fd6e2418971f673c82c6bf7f05ecac55a514090c57e94c3c37367c26f56",
    ),
    _nuccore(
        "WSSV-TW",
        "AF440570.1",
        307287,
        "0b1a766b502b220b3776081c8373ffc826cb4d319e93b5452a88120ee117eed1",
    ),
    _nuccore(
        "WSSV-CN",
        "NC_075105.1",
        305119,
        "aa1765ac64ae3b05d25d58e11204b49da9da7d69c69a8c9aea6ea1ac2e874278",
    ),
    _nuccore(
        "WSSV-TH",
        "AF369029.2",
        292967,
        "2fa88b28f5e1fb7fb7d15b5c88597fe2fccd78a10035f25333bffc4c04b903c0",
    ),
    _nuccore(
        "JP01A",
        "AP027278.1",
        299976,
        "45ace2b8c39cfb301155364fac69cdd09e78fa123e9b6f849077ad611c4d044d",
    ),
    _nuccore(
        "JP01B",
        "AP027279.1",
        293923,
        "1c9ee6665c4a81b21ab12ac72804133926a0f8d3d84f76e6a4f285e3370b3e63",
    ),
    _nuccore(
        "Pc2020",
        "AP027284.1",
        298496,
        "2955b99d0290a70bc4b3dd88e3afb4ccac885ec34c75344ce0034c4cc7783642",
    ),
    _nuccore(
        "E1",
        "AP027286.1",
        289353,
        "fd00a8f24b1adca42c024cde96bd4443d24e6009e9317503d98ccba263c9267c",
    ),
    _nuccore(
        "0722-1",
        "AP027288.1",
        288252,
        "efd7076a2191dbbe35dd7a2bc115ab29f741e598929af8bb4cb41edf7a169ea0",
    ),
    _nuccore(
        "CN03",
        "KT995471.1",
        284148,
        "1bb508c62391bf6dbbba88a95e120c5854a65670aff77ca3fa58e38c134276db",
    ),
    _nuccore(
        "CN04",
        "KY827813.1",
        281054,
        "4c8b18da1eb5f787a40286d60eefa1eabfe240565cde60e6b609155165f04279",
    ),
    _nuccore(
        "WSSV-AU",
        "MF768985.1",
        285973,
        "b560f2f996c269ca370b7511ce47ce5f9da3a6ab188ee80ea606b65503323c30",
    ),
    _prepared(
        "EU129",
        "SRR14509867",
        "EU129.fa",
        282445,
        "d6a8687b77b55e36624127b0005b2c9be7235adb096821d9efe48ac0a51342a6",
        sra_accession="SRR14509867",
    ),
    _prepared(
        "GCF7",
        "SRR12919258",
        "GCF7.fa",
        289862,
        "a68521ed04b7caa7812cdb6bf5add544a0337a164829e1c06cea105ff7e07c5a",
        sra_accession="SRR12919258",
    ),
    _prepared(
        "MES-753",
        "SRR17256726",
        "MES-753.fa",
        296593,
        "8a5d605f5a2e545da929b28ac3936915bd7d2c81b466341d91bb460b1e1799b4",
        sra_accession="SRR17256726",
    ),
    _prepared(
        "Shantou2019",
        "Shantou2020",
        "Shantou2019.fa",
        284424,
        "7328d54ea23fa76cf1ac5d7f6190275a69f638e0c26f73b336294bfa2a7ca5c8",
        sra_accession=None,
    ),
    _prepared(
        "POMZ1",
        "SRR8144089",
        "POMZ1.fa",
        291020,
        "ced387726d94afd75c5d338321f90eec1e124ffdd90a39bb17d92cfff65c14f4",
        sra_accession="SRR8144089",
    ),
    _prepared(
        "POMZ4",
        "SRR8144084",
        "POMZ4.fa",
        289337,
        "81611192394bfc0bdf8da42a2451d4c523ecd2f9707650ebdef9acc06dd486ef",
        sra_accession="SRR8144084",
    ),
    _prepared(
        "MG18PR-0187-N40S",
        "SRR22022264",
        "MG18PR-0187-N40S.fa",
        285171,
        "9821dd1e9039b27231dc6653a2eb43d51028fcad1eb17d5db1e23c2fbf9d8fff",
        sra_accession="SRR22022264",
    ),
    _prepared(
        "Angostura2013",
        "ERR5659803",
        "Angostura2013.fa",
        283826,
        "cd59b831b64a695b0bcf6681967be29cccd817536101699c22c014da3b66acd8",
        sra_accession="ERR5659803",
    ),
)


def _resource_bytes(resource: Mapping[str, Any]) -> bytes:
    data = resource.get("data")
    if not isinstance(data, str) or not data:
        raise ValueError("FASTA resource has no payload")
    encoding = resource.get("encoding")
    if encoding not in {None, "base64"}:
        raise ValueError(f"Unsupported FASTA resource encoding: {encoding}")
    try:
        return base64.b64decode(data, validate=True)
    except ValueError as exc:
        raise ValueError("FASTA resource contains invalid base64 data") from exc


def _split_fasta(raw: bytes) -> dict[str, bytes]:
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("FASTA input is not UTF-8") from exc
    starts = [match.start() for match in re.finditer(r"(?m)^>", text)]
    if not starts or text[: starts[0]].strip():
        raise ValueError("FASTA input does not start with a record header")
    records: dict[str, bytes] = {}
    for index, start in enumerate(starts):
        end = starts[index + 1] if index + 1 < len(starts) else len(text)
        record_text = text[start:end]
        record_id = record_text[1:].splitlines()[0].split()[0]
        if not record_id or record_id in records:
            raise ValueError(f"FASTA input has a missing or duplicate ID: {record_id!r}")
        records[record_id] = record_text.encode("utf-8")
    return records


def _fasta_identity(raw: bytes) -> tuple[str, int]:
    records = _split_fasta(raw)
    if len(records) != 1:
        raise ValueError("Each WSSV comparison source must contain exactly one record")
    record_id, record = next(iter(records.items()))
    lines = record.decode("utf-8").splitlines()
    sequence = "".join(lines[1:])
    sequence = re.sub(r"\s+", "", sequence)
    if not sequence or re.search(r"[^A-Za-z*.-]", sequence):
        raise ValueError(f"FASTA record {record_id!r} has invalid sequence text")
    return record_id, len(sequence.replace("-", "").replace(".", ""))


def _blast_evidence(
    session: Mapping[str, Any],
    source: WssvFastaSource,
    index: int,
) -> tuple[set[str], int]:
    options = session["renderRequest"]["diagramOptions"]
    blast_ref = options["conservationBlastFiles"][index]
    blast_resource = session["resources"][blast_ref["resourceId"]]
    rows = [
        line.split("\t")
        for line in _resource_bytes(blast_resource).decode("utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if not rows or any(len(row) < 10 for row in rows):
        raise ValueError(f"{source.label}: BLAST evidence is empty or malformed")
    query_ids = {row[0] for row in rows}
    max_query_coordinate = max(
        max(int(row[6]), int(row[7]))
        for row in rows
    )
    return query_ids, max_query_coordinate


def _cache_hashes(session: Mapping[str, Any]) -> dict[str, str]:
    entries = session.get("losatCache", {}).get("entries", [])
    return {
        str(entry.get("filename")): str(entry.get("queryCanonicalHash"))
        for entry in entries
        if isinstance(entry, Mapping)
        and entry.get("flow") == "circular-conservation"
        and entry.get("filename")
        and entry.get("queryCanonicalHash")
    }


def _validate_source_payload(
    session: Mapping[str, Any],
    source: WssvFastaSource,
    index: int,
    raw: bytes,
) -> None:
    actual_hash = hashlib.sha256(raw).hexdigest()
    if actual_hash != source.sha256:
        raise ValueError(
            f"{source.label}: FASTA SHA-256 {actual_hash} does not match "
            f"the pinned value {source.sha256}"
        )
    record_id, sequence_length = _fasta_identity(raw)
    if record_id != source.record_id or sequence_length != source.length:
        raise ValueError(
            f"{source.label}: expected {source.record_id}/{source.length} bp, "
            f"found {record_id}/{sequence_length} bp"
        )
    cache_filename = f"{source.label}.circular_conservation.losatn.tsv"
    cache_hash = _cache_hashes(session).get(cache_filename)
    if cache_hash != source.sha256:
        raise ValueError(
            f"{source.label}: pinned FASTA does not match its LOSAT cache identity"
        )
    query_ids, max_query_coordinate = _blast_evidence(session, source, index)
    if query_ids != {source.record_id} or max_query_coordinate != source.length:
        raise ValueError(
            f"{source.label}: BLAST evidence does not match "
            f"{source.record_id}/{source.length} bp"
        )


def _embedded_payloads(
    session: Mapping[str, Any],
) -> dict[str, bytes]:
    web_files = session.get("webFiles")
    resources = session.get("resources")
    if not isinstance(web_files, Mapping) or not isinstance(resources, Mapping):
        return {}
    refs = web_files.get("conservationLosatFastaSources")
    if not isinstance(refs, list) or len(refs) != len(WSSV_FASTA_SOURCES):
        return {}
    payloads: dict[str, bytes] = {}
    for source, resource_id in zip(WSSV_FASTA_SOURCES, refs, strict=True):
        resource = resources.get(resource_id)
        if not isinstance(resource, Mapping):
            return {}
        payloads[source.label] = _resource_bytes(resource)
    return payloads


def validate_embedded_wssv_fastas(
    session: Mapping[str, Any],
) -> tuple[WssvFastaSource, ...]:
    """Validate all embedded WSSV FASTAs against cache and BLAST evidence."""

    validate_session(session)
    options = session["renderRequest"]["diagramOptions"]
    labels = options.get("conservationLabels")
    expected_labels = [source.label for source in WSSV_FASTA_SOURCES]
    if labels != expected_labels:
        raise ValueError("WSSV conservation label order does not match the source manifest")
    payloads = _embedded_payloads(session)
    if len(payloads) != len(WSSV_FASTA_SOURCES):
        raise ValueError("WSSV session does not contain all 20 FASTA resources")
    for index, source in enumerate(WSSV_FASTA_SOURCES):
        _validate_source_payload(session, source, index, payloads[source.label])
    return WSSV_FASTA_SOURCES


def _load_source_payloads(
    session: Mapping[str, Any],
    *,
    ncbi_fasta: Path | None,
    mag_dir: Path,
) -> dict[str, bytes]:
    payloads = _embedded_payloads(session)
    ncbi_records = (
        _split_fasta(ncbi_fasta.read_bytes())
        if ncbi_fasta is not None
        else {}
    )
    missing: list[str] = []
    for index, source in enumerate(WSSV_FASTA_SOURCES):
        raw = payloads.get(source.label)
        if raw is not None:
            try:
                _validate_source_payload(session, source, index, raw)
                continue
            except ValueError:
                pass
        if source.local_filename:
            source_path = mag_dir / source.local_filename
            raw = source_path.read_bytes() if source_path.is_file() else None
        else:
            raw = ncbi_records.get(source.record_id)
        if raw is None:
            missing.append(source.label)
            continue
        _validate_source_payload(session, source, index, raw)
        payloads[source.label] = raw
    if missing:
        raise FileNotFoundError(
            "Missing exact WSSV FASTA source(s): "
            f"{', '.join(missing)}. Supply --ncbi-fasta and/or --mag-dir."
        )
    return payloads


def restore_wssv_gallery_fastas(
    session: dict[str, Any],
    *,
    ncbi_fasta: Path | None,
    mag_dir: Path,
) -> bool:
    """Embed exact comparison FASTAs without changing unrelated session data."""

    validate_session(session)
    payloads = _load_source_payloads(
        session,
        ncbi_fasta=ncbi_fasta,
        mag_dir=mag_dir,
    )
    resources = session["resources"]
    web_files = session.setdefault("webFiles", {})
    original_names = web_files.setdefault("resourceOriginalNames", {})
    resource_ids: list[str] = []
    changed = False
    for index, source in enumerate(WSSV_FASTA_SOURCES, start=1):
        resource_id = f"conservation-losat-fasta-files-{index}"
        raw = payloads[source.label]
        resource = {
            "kind": "conservation-fasta-file",
            "name": f"{resource_id}-{source.filename}",
            "type": "application/octet-stream",
            "size": len(raw),
            "lastModified": 0,
            "encoding": "base64",
            "data": base64.b64encode(raw).decode("ascii"),
        }
        changed = changed or resources.get(resource_id) != resource
        resources[resource_id] = resource
        changed = changed or original_names.get(resource_id) != source.filename
        original_names[resource_id] = source.filename
        resource_ids.append(resource_id)
    previous_refs = web_files.get("conservationLosatFastaSources")
    web_files["conservationLosatFastaSources"] = resource_ids
    validate_embedded_wssv_fastas(session)
    return changed or previous_refs != resource_ids


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Embed the exact comparison FASTAs identified by the WSSV Gallery "
            "session's LOSAT hashes and BLAST evidence."
        )
    )
    parser.add_argument(
        "--session",
        type=Path,
        default=SESSION_PATH,
        help="WSSV Gallery session to update.",
    )
    parser.add_argument(
        "--ncbi-fasta",
        type=Path,
        default=None,
        help="Multi-record NCBI EFetch FASTA containing the 12 accession records.",
    )
    parser.add_argument(
        "--mag-dir",
        type=Path,
        default=DEFAULT_MAG_DIR,
        help="Directory containing the eight prepared comparison FASTAs.",
    )
    args = parser.parse_args(argv)

    session = load_session(args.session)
    before = args.session.read_bytes()
    restore_wssv_gallery_fastas(
        session,
        ncbi_fasta=args.ncbi_fasta,
        mag_dir=args.mag_dir,
    )
    write_session_json(args.session, session)
    after = args.session.read_bytes()
    status = "updated" if after != before else "already current"
    print(f"{args.session}: {status}; validated {len(WSSV_FASTA_SOURCES)} FASTAs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
