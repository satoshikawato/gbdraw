#!/usr/bin/env python3
"""Rebuild or verify the 1 kbp mean-depth tutorial fixture."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = REPO_ROOT / "tests" / "test_inputs" / "AP027133.DRR394922.depth.tsv"
RECORD_PATH = REPO_ROOT / "tests" / "test_inputs" / "AP027133.gb"
OUTPUT_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "depth-1kb"
    / "AP027133.DRR394922.depth-1kb.tsv"
)
BIN_SIZE = 1_000


def render() -> bytes:
    record = SeqIO.read(RECORD_PATH, "genbank")
    if record.id != "AP027133.1":
        raise ValueError(f"Expected AP027133.1, found {record.id}.")
    reference: str | None = None
    expected_position = 1
    bin_start = 1
    bin_total = 0
    bin_count = 0
    rows = ["reference_name\tposition\tdepth"]

    def append_bin() -> None:
        nonlocal bin_start, bin_total, bin_count
        if reference is None or bin_count == 0:
            return
        rows.append(f"{reference}\t{bin_start}\t{bin_total / bin_count:.6f}")
        bin_start += bin_count
        bin_total = 0
        bin_count = 0

    with SOURCE_PATH.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Line {line_number} has fewer than three columns.")
            row_reference, position_text, depth_text = fields[:3]
            position = int(position_text)
            depth = int(depth_text)
            if depth < 0:
                raise ValueError(f"Line {line_number} has negative depth.")
            if reference is None:
                reference = row_reference
            elif row_reference != reference:
                raise ValueError("The source fixture must contain exactly one reference.")
            if position != expected_position:
                raise ValueError(
                    f"Expected consecutive position {expected_position}, found {position}."
                )
            bin_total += depth
            bin_count += 1
            expected_position += 1
            if bin_count == BIN_SIZE:
                append_bin()

    append_bin()
    if expected_position == 1:
        raise ValueError("The source fixture is empty.")
    if expected_position - 1 != len(record):
        raise ValueError(
            f"Depth rows ({expected_position - 1}) do not span the record ({len(record)})."
        )
    return ("\n".join(rows) + "\n").encode("utf-8")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check() -> int:
    expected = render()
    if not OUTPUT_PATH.is_file():
        print(f"missing: {OUTPUT_PATH.relative_to(REPO_ROOT)}", file=sys.stderr)
        return 1
    current = OUTPUT_PATH.read_bytes()
    if current != expected:
        print(
            f"stale: {OUTPUT_PATH.relative_to(REPO_ROOT)} "
            f"(current sha256={_digest(current)}, expected sha256={_digest(expected)})",
            file=sys.stderr,
        )
        return 1
    print("The 1 kbp mean-depth fixture is reproducible byte-for-byte.")
    return 0


def write() -> int:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_bytes(render())
    print(f"Wrote {OUTPUT_PATH.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="Verify the committed output without writing (the default).",
    )
    action.add_argument(
        "--write",
        action="store_true",
        help="Regenerate the committed derivative.",
    )
    args = parser.parse_args(argv)
    return write() if args.write else check()


if __name__ == "__main__":
    raise SystemExit(main())
