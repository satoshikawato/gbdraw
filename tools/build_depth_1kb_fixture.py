#!/usr/bin/env python3
"""Rebuild or verify the files derived from the AP027133 DRR394922 per-base depth fixture.

Three files, each from the same per-base depth table:

- the 1 kbp mean depth table: the arithmetic mean of each consecutive 1 kbp bin at its 1-based start;
- a signed bedGraph: log2 of each 1 kbp bin's mean depth over the genome's median per-base depth, to two
  decimal places;
- a 3-column BED file: one row per run of consecutive bases whose depth is below 10x.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import statistics
import sys
from pathlib import Path

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = REPO_ROOT / "tests" / "test_inputs" / "AP027133.DRR394922.depth.tsv"
RECORD_PATH = REPO_ROOT / "tests" / "test_inputs" / "AP027133.gb"
OUTPUT_DIR = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "depth-1kb"
OUTPUT_PATH = OUTPUT_DIR / "AP027133.DRR394922.depth-1kb.tsv"
LOG2_RATIO_PATH = OUTPUT_DIR / "AP027133.DRR394922.log2-ratio-1kb.bedgraph"
LOW_DEPTH_PATH = OUTPUT_DIR / "AP027133.DRR394922.below-10x.bed"
BIN_SIZE = 1_000
LOW_DEPTH_BELOW = 10


def read_depths() -> tuple[str, list[int]]:
    """The reference name and the per-base depths of the source table, one per base of the record."""

    record = SeqIO.read(RECORD_PATH, "genbank")
    if record.id != "AP027133.1":
        raise ValueError(f"Expected AP027133.1, found {record.id}.")
    reference: str | None = None
    depths: list[int] = []
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
            if position != len(depths) + 1:
                raise ValueError(f"Expected consecutive position {len(depths) + 1}, found {position}.")
            depths.append(depth)
    if reference is None:
        raise ValueError("The source fixture is empty.")
    if len(depths) != len(record):
        raise ValueError(f"Depth rows ({len(depths)}) do not span the record ({len(record)}).")
    return reference, depths


def bin_means(depths: list[int]) -> list[tuple[int, int, float]]:
    """Each consecutive bin as (0-based start, exclusive end, mean depth)."""

    return [
        (start, min(start + BIN_SIZE, len(depths)), statistics.fmean(depths[start : start + BIN_SIZE]))
        for start in range(0, len(depths), BIN_SIZE)
    ]


def render_bins(reference: str, depths: list[int]) -> bytes:
    rows = ["reference_name\tposition\tdepth"]
    rows += [f"{reference}\t{start + 1}\t{mean:.6f}" for start, _end, mean in bin_means(depths)]
    return ("\n".join(rows) + "\n").encode("utf-8")


def render_log2_ratio(reference: str, depths: list[int]) -> bytes:
    median = statistics.median(depths)
    rows = []
    for start, end, mean in bin_means(depths):
        value = f"{math.log2(mean / median):.2f}"
        rows.append(f"{reference}\t{start}\t{end}\t{'0.00' if value == '-0.00' else value}")
    return ("\n".join(rows) + "\n").encode("utf-8")


def render_low_depth_runs(reference: str, depths: list[int]) -> bytes:
    rows = []
    run_start: int | None = None
    for index, depth in enumerate([*depths, LOW_DEPTH_BELOW]):
        if depth < LOW_DEPTH_BELOW and run_start is None:
            run_start = index
        elif depth >= LOW_DEPTH_BELOW and run_start is not None:
            rows.append(f"{reference}\t{run_start}\t{index}")
            run_start = None
    return ("\n".join(rows) + "\n").encode("utf-8")


OUTPUTS = {
    OUTPUT_PATH: render_bins,
    LOG2_RATIO_PATH: render_log2_ratio,
    LOW_DEPTH_PATH: render_low_depth_runs,
}


def render() -> dict[Path, bytes]:
    reference, depths = read_depths()
    return {path: build(reference, depths) for path, build in OUTPUTS.items()}


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check() -> int:
    failed = False
    for path, expected in render().items():
        if not path.is_file():
            print(f"missing: {path.relative_to(REPO_ROOT)}", file=sys.stderr)
            failed = True
            continue
        current = path.read_bytes()
        if current != expected:
            print(
                f"stale: {path.relative_to(REPO_ROOT)} "
                f"(current sha256={_digest(current)}, expected sha256={_digest(expected)})",
                file=sys.stderr,
            )
            failed = True
    if failed:
        return 1
    print("The files derived from the AP027133 depth fixture are reproducible byte-for-byte.")
    return 0


def write() -> int:
    for path, data in render().items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        print(f"Wrote {path.relative_to(REPO_ROOT)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="Verify the committed outputs without writing (the default).",
    )
    action.add_argument(
        "--write",
        action="store_true",
        help="Regenerate the committed derivatives.",
    )
    args = parser.parse_args(argv)
    return write() if args.write else check()


if __name__ == "__main__":
    raise SystemExit(main())
