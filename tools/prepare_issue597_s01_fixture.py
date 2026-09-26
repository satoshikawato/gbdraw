"""Disposable real Session recipe; never writes Gallery or reference artifacts."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    "GCF_000196095.1_ASM19609v1_genomic.gbff",
    "GCF_000354175.2_ASM35417v2_genomic.gbff",
    "GCF_000801275.2_ASM80127v1_genomic.gbff",
    "GCF_002021755.1_ASM202175v1_genomic.gbff",
    "GCF_002906475.1_ASM290647v1_genomic.gbff",
    "GCF_015097735.1_ASM1509773v1_genomic.gbff",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify existing SVG/Session without repeating LOSAT searches",
    )
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    records, provenance = [], []
    for name in SOURCES:
        path = ROOT / "tests/test_inputs" / name
        selected = list(SeqIO.parse(path, "genbank"))[:2]
        records.extend(selected)
        provenance.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha(path),
                "records": [
                    {
                        "id": r.id,
                        "length": len(r),
                        "sequenceSha256": hashlib.sha256(bytes(r.seq)).hexdigest(),
                        "features": len(r.features),
                        "translatedCDS": sum(
                            f.type == "CDS" and "translation" in f.qualifiers
                            for f in r.features
                        ),
                    }
                    for r in selected
                ],
            }
        )
    assert len(records) == 12 and len({r.id for r in records}) == 12
    assert len({hashlib.sha256(bytes(r.seq)).hexdigest() for r in records}) == 12
    assert (
        sum(
            sum(f.type == "CDS" and "translation" in f.qualifiers for f in r.features)
            for r in records
        )
        >= 25000
    )
    source = out / "real-12-records.gb"
    SeqIO.write(records, source, "genbank")
    table = out / "records.tsv"
    table.write_text(
        "gbk\trecord_id\torder\trow\tcolumn\n"
        + "".join(
            f"{source.name}\t#{i + 1}\t{i + 1}\t{i // 2 + 1}\t{i % 2 + 1}\n"
            for i in range(12)
        )
    )
    fixture = out / "real-full-pairwise.gbdraw-session.json.gz"
    subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "setup-losat"],
        cwd=ROOT,
        env={**os.environ, "XDG_CACHE_HOME": str(out / ".cache")},
        check=True,
    )
    runtime = out / ".cache/gbdraw/losat/0.1.0/x86_64-unknown-linux-gnu/LOSAT"
    lock = json.loads((ROOT / "gbdraw/data/losat-release.json").read_text())
    assert (
        sha(runtime) == lock["artifacts"]["x86_64-unknown-linux-gnu"]["binary_sha256"]
    )
    command = [
        sys.executable,
        "-m",
        "gbdraw.cli",
        "linear",
        "--records_table",
        str(table),
        "--protein_blastp_mode",
        "collinear",
        "--collinear_search_scope",
        "all",
        "--protein_blastp_candidate_limit",
        "5",
        "--collinear_min_anchors",
        "3",
        "--collinear_max_unit_gap",
        "2",
        "--collinear_max_diagonal_drift",
        "2",
        "--losatp_bin",
        str(runtime),
        "--losatp_threads",
        str(args.threads),
        "--feature_shape",
        "CDS=rectangle",
        "--legend",
        "bottom",
        "--scale_style",
        "ruler",
        "--ruler_on_axis",
        "--separate_strands",
        "-o",
        str(out / "real-full-pairwise"),
        "-f",
        "svg",
        "--session_output",
        str(fixture),
        "--overwrite",
    ]
    manifest = {
        "sourceHead": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "sources": provenance,
        "combinedSha256": sha(source),
        "command": command,
        "runtimeVersion": subprocess.check_output(
            [str(runtime), "--version"], text=True
        ).strip(),
        "runtimeSha256": sha(runtime),
        "releaseLockSha256": sha(ROOT / "gbdraw/data/losat-release.json"),
    }
    (out / "fixture-manifest.json").write_text(json.dumps(manifest, indent=2))
    if not args.verify_only:
        subprocess.run(command, cwd=ROOT, check=True)
    expanded = gzip.decompress(fixture.read_bytes())
    document = json.loads(expanded)
    features = sum(
        len(i.get("biologicalFeatures", []))
        for i in document["editorState"]["featureCatalog"]["items"]
    )
    assert len(document["results"]) == 1
    assert (
        document["results"][0]["content"]
        == (out / "real-full-pairwise.svg").read_text()
    )
    cache = document["losatCache"]["entries"]
    pairs = {
        (e["queryRecordInstanceKey"], e["subjectRecordInstanceKey"]) for e in cache
    }
    assert len(document["renderRequest"]["records"]) == 12 and features >= 25000
    record_keys = {
        record["recordKey"] for record in document["renderRequest"]["records"]
    }
    required_pairs = {(a, b) for a in record_keys for b in record_keys if a != b}
    nonself_pairs = {(a, b) for a, b in pairs if a != b}
    self_pairs = {(a, b) for a, b in pairs if a == b}
    assert nonself_pairs == required_pairs and len(nonself_pairs) == 12 * 11
    assert self_pairs == {(key, key) for key in record_keys}
    assert len(cache) == len(pairs) == 144
    assert len(expanded) <= 512 * 1024 * 1024
    assert fixture.stat().st_size <= 200 * 1024 * 1024
    assert document["version"] == 44 and document["renderRequest"]["schema"] == 8
    manifest.update(
        {
            "fixtureSha256": sha(fixture),
            "compressedBytes": fixture.stat().st_size,
            "expandedBytes": len(expanded),
            "records": 12,
            "biologicalFeatures": features,
            "directedPairs": len(nonself_pairs),
            "selfPairs": len(self_pairs),
            "cacheEntries": len(cache),
        }
    )
    (out / "fixture-manifest.json").write_text(json.dumps(manifest, indent=2))
    print(
        json.dumps(
            {k: v for k, v in manifest.items() if k not in ("sources", "command")}
        )
    )


if __name__ == "__main__":
    main()
