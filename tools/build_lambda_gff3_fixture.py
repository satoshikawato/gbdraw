#!/usr/bin/env python3
"""Rebuild or verify the full-record Lambda GFF3/FASTA tutorial fixture."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from textwrap import wrap
from urllib.parse import quote

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL_DATA_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
SOURCE_PATH = TUTORIAL_DATA_ROOT / "lambda" / "NC_001416.gb"
OUTPUT_DIR = TUTORIAL_DATA_ROOT / "lambda-gff3"
GFF3_PATH = OUTPUT_DIR / "NC_001416.gff3"
FASTA_PATH = OUTPUT_DIR / "NC_001416.fna"


def _attribute(name: str, value: str | None) -> str | None:
    if not value:
        return None
    return f"{name}={quote(value, safe='')}"


def _qualifier(feature, name: str) -> str | None:
    values = feature.qualifiers.get(name, [])
    return str(values[0]) if values else None


def render_gff3() -> bytes:
    record = SeqIO.read(SOURCE_PATH, "genbank")
    selected = [feature for feature in record.features if feature.type in {"gene", "CDS"}]
    gene_ids: dict[tuple[str, str], str] = {}
    gene_number = 0
    cds_number = 0
    lines = [
        "##gff-version 3",
        f"##sequence-region {record.id} 1 {len(record)}",
    ]

    for feature in selected:
        if feature.type != "gene":
            continue
        gene_number += 1
        gene_id = f"gene-{gene_number:03d}"
        for qualifier in ("gene", "locus_tag"):
            value = _qualifier(feature, qualifier)
            if value:
                gene_ids[(qualifier, value)] = gene_id

    gene_number = 0
    for feature in selected:
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"
        attributes: list[str | None]
        if feature.type == "gene":
            gene_number += 1
            feature_id = f"gene-{gene_number:03d}"
            attributes = [
                f"ID={feature_id}",
                _attribute("Name", _qualifier(feature, "gene") or _qualifier(feature, "locus_tag")),
                _attribute("gene", _qualifier(feature, "gene")),
                _attribute("locus_tag", _qualifier(feature, "locus_tag")),
            ]
            phase = "."
        else:
            cds_number += 1
            feature_id = f"cds-{cds_number:03d}"
            parent = next(
                (
                    gene_ids[(qualifier, value)]
                    for qualifier in ("gene", "locus_tag")
                    if (value := _qualifier(feature, qualifier))
                    and (qualifier, value) in gene_ids
                ),
                None,
            )
            attributes = [
                f"ID={feature_id}",
                _attribute("Parent", parent),
                _attribute("gene", _qualifier(feature, "gene")),
                _attribute("locus_tag", _qualifier(feature, "locus_tag")),
                _attribute("protein_id", _qualifier(feature, "protein_id")),
                _attribute("product", _qualifier(feature, "product")),
                _attribute("translation", _qualifier(feature, "translation")),
            ]
            phase = str(int(_qualifier(feature, "codon_start") or "1") - 1)
        encoded_attributes = ";".join(attribute for attribute in attributes if attribute)
        lines.append(
            f"{record.id}\tRefSeq-extract\t{feature.type}\t{start}\t{end}\t.\t"
            f"{strand}\t{phase}\t{encoded_attributes}"
        )
    return ("\n".join(lines) + "\n").encode("utf-8")


def render_fasta() -> bytes:
    record = SeqIO.read(SOURCE_PATH, "genbank")
    lines = [f">{record.id} {record.description}"]
    lines.extend(wrap(str(record.seq), width=70))
    return ("\n".join(lines) + "\n").encode("ascii")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check() -> int:
    expected = {
        GFF3_PATH: render_gff3(),
        FASTA_PATH: render_fasta(),
    }
    failures: list[str] = []
    for path, rendered in expected.items():
        if not path.is_file():
            failures.append(f"missing: {path.relative_to(REPO_ROOT)}")
            continue
        current = path.read_bytes()
        if current != rendered:
            failures.append(
                f"stale: {path.relative_to(REPO_ROOT)} "
                f"(current sha256={_digest(current)}, expected sha256={_digest(rendered)})"
            )
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("Full-record Lambda GFF3/FASTA fixture is reproducible byte-for-byte.")
    return 0


def write() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GFF3_PATH.write_bytes(render_gff3())
    FASTA_PATH.write_bytes(render_fasta())
    print(f"Wrote {GFF3_PATH.relative_to(REPO_ROOT)}")
    print(f"Wrote {FASTA_PATH.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="Verify committed outputs without writing (the default).",
    )
    action.add_argument(
        "--write",
        action="store_true",
        help="Regenerate the two committed derived files.",
    )
    args = parser.parse_args(argv)
    return write() if args.write else check()


if __name__ == "__main__":
    raise SystemExit(main())
