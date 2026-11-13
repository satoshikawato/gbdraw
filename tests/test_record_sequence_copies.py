"""Whole copies of a record sequence stay where they are today, and their list only shrinks.

A record's ``seq`` reads its bases from the input file on demand
(:mod:`gbdraw.io.record_sequences`), so ``str(record.seq)`` or ``bytes(record.seq)``
reads and holds the whole sequence, 1 byte per base or more. This guard finds,
outside that module and ``gbdraw/web/``:

- ``str()`` or ``bytes()`` of a ``.seq`` attribute or of a name ``seq`` or ``sequence``;
- every ``reverse_complement()`` call (of a sequence or of a whole record);
- ``+`` with a ``.seq`` operand, and ``.join`` over ``.seq`` values;
- ``SeqIO.write`` (it writes every base of every record), and ``SeqIO.parse`` or
  ``SeqIO.read`` outside the record readers of ``gbdraw/io/genome.py`` (they build
  whole sequences beside the owner; the readers there parse lines without sequence).

It compares them with ``SITES``. A new site fails: read the sequence in slices
instead. A site that is gone fails too, until its entry is removed here. A site
is a place in the code; some read less than a whole record (a window, a feature),
and stay listed until they read through views.
"""

from __future__ import annotations

import ast
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "gbdraw"
OWNER = PACKAGE / "io" / "record_sequences.py"
# Their SeqIO.parse calls read the lines the owner hands them, without sequence lines.
HEADER_PARSERS = {"gbdraw/io/genome.py"}
SEQUENCE_NAMES = {"seq", "sequence"}

# file::function::kind -> number of sites. Today's whole copies; each later change that
# reads such a sequence in pieces removes its entry.
SITES: dict[str, int] = {
    "gbdraw/analysis/protein_colinearity.py::_clean_protein_sequence::str()": 1,
    "gbdraw/analysis/skew.py::_counted_bases::reverse_complement()": 1,
    "gbdraw/analysis/skew.py::counted_sequence::bytes()": 1,
    "gbdraw/analysis/skew.py::counted_sequence::str()": 1,
    "gbdraw/comparisons/circular_losat.py::comparison_query_fasta::join": 1,
    "gbdraw/comparisons/circular_losat.py::comparison_query_fasta::str()": 1,
    "gbdraw/comparisons/locate.py::locate_query::reverse_complement()": 1,
    "gbdraw/comparisons/locate.py::locate_query::str()": 1,
    "gbdraw/comparisons/locate.py::plan_locate_search::str()": 1,
    "gbdraw/comparisons/losat_jobs.py::losat_fasta_text::str()": 1,
    "gbdraw/comparisons/losat_jobs.py::losat_search_frame_fasta::reverse_complement()": 1,
    "gbdraw/comparisons/losat_jobs.py::losat_search_frame_fasta::str()": 1,
    "gbdraw/core/record_metadata.py::record_content_digest::bytes()": 1,
    "gbdraw/crop_genbank.py::gbk_to_seqrecord::SeqIO.parse": 1,
    "gbdraw/crop_genbank.py::write_cropped_genbank::SeqIO.write": 1,
    "gbdraw/io/comparison_sequences.py::_sequence_text::str()": 1,
    "gbdraw/io/comparison_sequences.py::read_comparison_sequence_file::SeqIO.parse": 1,
    "gbdraw/io/genome.py::record_sequence_text::str()": 1,
    "gbdraw/io/record_select.py::reverse_records::reverse_complement()": 1,
    "gbdraw/io/regions.py::_origin_crossing_view::+": 1,
    "gbdraw/session_request_codec.py::_genbank_bytes::SeqIO.write": 1,
    "gbdraw/session_request_codec.py::_with_crop_edits::bytes()": 1,
    "gbdraw/session_request_codec.py::_with_crop_edits::reverse_complement()": 1,
    "gbdraw/web_support/record_spans.py::read_record_span_sequences.source_bases::str()": 1,
    "gbdraw/web_support/sequence_search.py::read_search_queries::str()": 1,
}


def _is_sequence(node: ast.AST) -> bool:
    return (isinstance(node, ast.Attribute) and node.attr == "seq") or (
        isinstance(node, ast.Name) and node.id in SEQUENCE_NAMES)


class _Finder(ast.NodeVisitor):
    def __init__(self, path: str) -> None:
        self.path = path
        self.scope: list[str] = []
        self.sites: Counter[str] = Counter()

    def _add(self, kind: str) -> None:
        self.sites[f"{self.path}::{'.'.join(self.scope) or '<module>'}::{kind}"] += 1

    def _scoped(self, node: ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef) -> None:
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_FunctionDef = _scoped
    visit_AsyncFunctionDef = _scoped
    visit_ClassDef = _scoped

    def visit_Call(self, node: ast.Call) -> None:
        func = node.func
        if isinstance(func, ast.Name) and func.id in ("str", "bytes") and len(node.args) == 1 and _is_sequence(node.args[0]):
            self._add(f"{func.id}()")
        if isinstance(func, ast.Attribute) and func.attr == "reverse_complement":
            self._add("reverse_complement()")
        if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) and func.value.id == "SeqIO" and (
                func.attr == "write" or (func.attr in ("parse", "read") and self.path not in HEADER_PARSERS)):
            self._add(f"SeqIO.{func.attr}")
        if isinstance(func, ast.Attribute) and func.attr == "join" and node.args and any(
                isinstance(child, ast.Attribute) and child.attr == "seq" for child in ast.walk(node.args[0])):
            self._add("join")
        self.generic_visit(node)

    def visit_BinOp(self, node: ast.BinOp) -> None:
        if isinstance(node.op, ast.Add) and any(
                isinstance(side, ast.Attribute) and side.attr == "seq" for side in (node.left, node.right)):
            self._add("+")
        self.generic_visit(node)


def _sites() -> Counter[str]:
    sites: Counter[str] = Counter()
    for path in sorted(PACKAGE.rglob("*.py")):
        if path == OWNER or "web" in path.relative_to(PACKAGE).parts[:1]:
            continue
        finder = _Finder(path.relative_to(ROOT).as_posix())
        finder.visit(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        sites.update(finder.sites)
    return sites


def test_whole_record_sequence_copies_only_shrink() -> None:
    found = _sites()
    new = {site: count for site, count in found.items() if count > SITES.get(site, 0)}
    assert not new, (
        "New whole copies of a record sequence (read it in slices through record.seq instead): "
        + ", ".join(sorted(new))
    )
    gone = {site: count for site, count in SITES.items() if found.get(site, 0) < count}
    assert not gone, "Whole copies that are gone; lower or remove their SITES entries: " + ", ".join(sorted(gone))


def test_the_detector_finds_each_kind_of_copy() -> None:
    finder = _Finder("example.py")
    finder.visit(ast.parse(
        "def f(record, seq, parts):\n"
        "    a = str(record.seq)\n"
        "    b = bytes(seq)\n"
        "    c = record.seq + record.seq\n"
        "    d = ''.join(str(part.seq) for part in parts)\n"
        "    e = str(record.seq[1:5])\n"
        "    g = len(record.seq)\n"
        "    h = record.seq.reverse_complement()\n"
        "    i = record.reverse_complement(id=True)\n"
        "    j = str(sequence)\n"
        "    SeqIO.write(parts, 'out.gb', 'genbank')\n"
        "    k = list(SeqIO.parse('in.gb', 'genbank'))\n"
        "    m = str(record.id)\n"
    ))
    assert finder.sites == Counter({
        "example.py::f::str()": 3, "example.py::f::bytes()": 1, "example.py::f::+": 1, "example.py::f::join": 1,
        "example.py::f::reverse_complement()": 2, "example.py::f::SeqIO.write": 1, "example.py::f::SeqIO.parse": 1,
    })
