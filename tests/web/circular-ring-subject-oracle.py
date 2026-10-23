"""Emit the CLI ring subject database and raw cache key for the Web parity test (Owner D-25).

For each case it writes the reference files, plans the CLI request that draws
the same records as the Web state the test sets (records OFF, a single
presentation, a crop, GFF3 with an unnamed FASTA entry, a UTF-8 byte order
mark), captures the subject FASTA of the ring search and the raw cache key, and
runs the Worker helper ``read_all_sequence_records`` on the file for the reply
the page gets when its own reader declines the file.
"""

from __future__ import annotations

import base64
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import gbdraw.comparisons.circular_losat as circular_losat  # noqa: E402
from gbdraw.api import (  # noqa: E402
    CircularDiagramOptions,
    CircularDiagramRequest,
    LosatSearchOptions,
    RecordInput,
)
from gbdraw.api.request_render import plan_circular_request  # noqa: E402
from gbdraw.api.requests import GenBankInputSource, GffFastaInputSource  # noqa: E402
from gbdraw.io.record_select import parse_record_selector  # noqa: E402
from gbdraw.io.regions import parse_region_spec  # noqa: E402

QUERY = ">QUERY1\nACGTACGTACGTACGTACGT\n"
SEQUENCES = {"SUBA": "acgtacgtac" * 7, "SUBB": "ggccttaagg" * 9, "SUBC": "ttgacatgca" * 8}


def genbank(record_id: str, sequence: str) -> str:
    lines = [
        f"LOCUS       {record_id:<16}{len(sequence):>11} bp    DNA     linear   UNK 01-JAN-1980",
        f"DEFINITION  {record_id} record.",
        f"ACCESSION   {record_id}",
        f"VERSION     {record_id}.1",
        "KEYWORDS    .",
        "SOURCE      .",
        "  ORGANISM  .",
        "            .",
        "FEATURES             Location/Qualifiers",
        "ORIGIN",
    ]
    for start in range(0, len(sequence), 60):
        chunk = sequence[start:start + 60]
        lines.append(f"{start + 1:>9} " + " ".join(chunk[i:i + 10] for i in range(0, len(chunk), 10)))
    return "\n".join([*lines, "//", ""])


def helpers() -> dict[str, object]:
    source = (ROOT / "gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0].replace("\\\\", "\\"), namespace)
    return namespace


def main() -> None:
    searches: list[str] = []

    def fake_search(_spec, _query, subject, **_kwargs):
        searches.append(subject)
        return ""

    circular_losat.run_losat_search = fake_search
    read_all = helpers()["read_all_sequence_records"]
    multi = "".join(genbank(record_id, SEQUENCES[record_id]) for record_id in ("SUBA", "SUBB", "SUBC")).encode()
    gff = b"##gff-version 3\nSUBB\t.\tgene\t5\t40\t.\t+\t.\tID=gene1\n"
    fasta = "".join(f">{record_id} entry\n{SEQUENCES[record_id]}\n" for record_id in ("SUBA", "SUBB")).encode()
    cases = []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        query = root / "query.fasta"
        query.write_text(QUERY, encoding="ascii")
        specs = [
            ("records not drawn", {"source": multi}, "genbank",
             lambda paths: [RecordInput(source=GenBankInputSource(path=paths["source"]), selector=parse_record_selector(s))
                            for s in ("#1", "#3")]),
            ("single presentation", {"source": multi}, "genbank",
             lambda paths: [RecordInput(source=GenBankInputSource(path=paths["source"]),
                                        selector=parse_record_selector("#2"))]),
            ("crop", {"source": multi}, "genbank",
             lambda paths: [RecordInput(source=GenBankInputSource(path=paths["source"]),
                                        region=parse_region_spec("SUBB.1:11-50"))]),
            ("gff3 with an unnamed fasta entry", {"source": gff, "fasta": fasta}, "fasta",
             lambda paths: [RecordInput(source=GffFastaInputSource(paths["source"], paths["fasta"]))]),
            ("byte order mark", {"source": b"\xef\xbb\xbf" + multi}, "genbank",
             lambda paths: [RecordInput(source=GenBankInputSource(path=paths["source"]),
                                        selector=parse_record_selector("#1"))]),
        ]
        for index, (name, files, fmt, inputs) in enumerate(specs):
            paths = {}
            for role, data in files.items():
                path = root / f"case{index}-{role}"
                path.write_bytes(data)
                paths[role] = str(path)
            searches.clear()
            plan = plan_circular_request(CircularDiagramRequest(
                records=tuple(inputs(paths)),
                options=CircularDiagramOptions(
                    losat_search=LosatSearchOptions(program="losatn"),
                    conservation_sequence_files=[str(query)],
                ),
            ))
            (entry,) = plan.losat_cache_entries
            subject_path = paths["fasta" if fmt == "fasta" else "source"]
            cases.append({
                "name": name,
                "files": {role: base64.b64encode(data).decode("ascii") for role, data in files.items()},
                "cliSubjectFasta": searches[0],
                "cliCacheKey": entry["key"],
                "helperReply": json.loads(read_all(subject_path, fmt)),
            })
    json.dump({"query": QUERY, "cases": cases}, sys.stdout)


if __name__ == "__main__":
    main()
