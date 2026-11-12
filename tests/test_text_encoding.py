"""The encoding gbdraw reads and writes text in.

Every text file gbdraw/ opens, and every program output it reads as text, names
its encoding: without ``encoding=`` Python uses the locale's, which differs
between machines (cp1252 or cp932 on Windows).
A GenBank, EMBL, GFF3 or FASTA file is read as UTF-8, and one that is not UTF-8
in the locale's encoding.
"""

from __future__ import annotations

import ast
import locale
import os
import subprocess
import sys
from pathlib import Path

import pytest

from gbdraw.crop_genbank import gbk_to_seqrecord
from gbdraw.exceptions import ParseError
from gbdraw.io.comparison_sequences import read_comparison_sequence_file
from gbdraw.io.genome import parse_gff_fasta, read_flat_file_records, read_sequence_file_records

PACKAGE = Path(__file__).resolve().parents[1] / "gbdraw"
# ``open`` of these modules opens no text file: os takes flags, the others open no file text.
_NOT_TEXT_OPENERS = {"os", "tarfile", "webbrowser", "zipfile"}
# ``open`` of these modules reads bytes unless its mode has a "t".
_BYTES_UNLESS_T = {"gzip", "bz2", "lzma"}
# subprocess functions that decode a program's output when ``text`` or ``universal_newlines`` is true.
_SUBPROCESS_RUNNERS = {"run", "Popen", "call", "check_call", "check_output"}


def _mode(call: ast.Call, position: int, default: str) -> str:
    """The mode a call passes, ``"?"`` when it is not a literal."""

    node = next((keyword.value for keyword in call.keywords if keyword.arg == "mode"), None)
    if node is None:
        if len(call.args) <= position:
            return default
        node = call.args[position]
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else "?"


def _opens_text_without_encoding(call: ast.Call) -> bool:
    for keyword in call.keywords:
        if keyword.arg == "encoding":
            return isinstance(keyword.value, ast.Constant) and keyword.value.value is None
    func = call.func
    if isinstance(func, ast.Name):
        return func.id == "open" and "b" not in _mode(call, 1, "r")
    if not isinstance(func, ast.Attribute):
        return False
    if func.attr in {"read_text", "write_text"}:
        return True
    # The module a dotted call starts from: ``zipfile`` for ``zipfile.ZipFile.open(archive, info)``.
    root = func.value
    while isinstance(root, ast.Attribute):
        root = root.value
    module = root.id if isinstance(root, ast.Name) else None
    if module == "subprocess" and func.attr in _SUBPROCESS_RUNNERS:
        return any(
            keyword.arg in {"text", "universal_newlines"}
            and not (isinstance(keyword.value, ast.Constant) and not keyword.value.value)
            for keyword in call.keywords
        )
    if func.attr != "open" or module in _NOT_TEXT_OPENERS:
        return False
    if module in _BYTES_UNLESS_T:
        return "t" in _mode(call, 1, "rb")
    # io.open(file, mode) as the built-in; Path.open(mode) and Traversable.open(mode) take the mode first.
    return "b" not in _mode(call, 1 if module == "io" else 0, "r")


def _offenders(source: str, name: str) -> list[int]:
    tree = ast.parse(source, filename=name)
    return [node.lineno for node in ast.walk(tree) if isinstance(node, ast.Call) and _opens_text_without_encoding(node)]


@pytest.mark.parametrize(
    ("call", "flagged"),
    [
        ("open(path)", True),
        ("open(path, 'w')", True),
        ("open(path, mode)", True),
        ("open(path, encoding=None)", True),
        ("open(path, 'rb')", False),
        ("open(path, 'w', encoding='utf-8')", False),
        ("io.open(path, 'rt')", True),
        ("gzip.open(path, 'rt')", True),
        ("gzip.open(path, mode='wt')", True),
        ("gzip.open(path)", False),
        ("gzip.open(path, 'rt', encoding='utf-8')", False),
        ("Path(path).open()", True),
        ("path.open('w')", True),
        ("path.open('rb')", False),
        ("path.read_text()", True),
        ("resources.files('gbdraw.data').joinpath('x.json').read_text()", True),
        ("path.write_text(text, encoding='utf-8')", False),
        ("os.open(path, os.O_RDONLY)", False),
        ("webbrowser.open(url)", False),
        ("tarfile.open(path, 'r:gz')", False),
        ("zipfile.ZipFile.open(archive, info, 'w')", False),
        ("archive.open(info, 'w')", True),
        ("pathlib.Path.open(path)", True),
        ("write_text(handle, text)", False),
        ("subprocess.run(command, capture_output=True, text=True)", True),
        ("subprocess.check_output(command, universal_newlines=True)", True),
        ("subprocess.run(command, capture_output=True)", False),
        ("subprocess.run(command, text=False)", False),
        ("subprocess.run(command, text=True, encoding='utf-8')", False),
    ],
)
def test_the_encoding_check_reads_each_form_of_open(call: str, flagged: bool) -> None:
    assert _offenders(call, "<case>") == ([1] if flagged else [])


def test_every_text_open_in_gbdraw_names_its_encoding() -> None:
    offenders = [
        f"{path.relative_to(PACKAGE.parent).as_posix()}:{line}"
        for path in sorted(PACKAGE.rglob("*.py"))
        for line in _offenders(path.read_text(encoding="utf-8"), str(path))
    ]
    assert offenders == [], "Pass encoding= (UTF-8 for gbdraw's own files and inputs) to: " + ", ".join(offenders)


TEXT = "caf\u00e9 \u00b5m"
GENBANK = """LOCUS       rec1                      12 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  rec1.
ACCESSION   rec1
VERSION     rec1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..12
                     /product="{text}"
ORIGIN
        1 atgaaataga aa
//
"""
FASTA = ">rec1 {text}\nATGAAATAGAAA\n"
GFF3 = "##gff-version 3\nrec1\t.\tCDS\t1\t12\t.\t+\t0\tID=cds1;product={text}\n"


def _write(tmp_path: Path, name: str, template: str, encoding: str) -> str:
    path = tmp_path / name
    path.write_bytes(template.format(text=TEXT).encode(encoding))
    return str(path)


def _read_each_reader(tmp_path: Path, encoding: str) -> dict[str, str]:
    genbank = _write(tmp_path, "rec1.gb", GENBANK, encoding)
    fasta = _write(tmp_path, "rec1.fasta", FASTA, encoding)
    gff3 = _write(tmp_path, "rec1.gff", GFF3, encoding)
    return {
        "GenBank": read_flat_file_records(genbank)[0].features[0].qualifiers["product"][0],
        "FASTA": read_sequence_file_records(fasta, "fasta")[0].description,
        "GFF3": parse_gff_fasta(gff3, fasta)[0].features[0].qualifiers["product"][0],
        "comparison genome": read_comparison_sequence_file(fasta).records[0].id,
        "crop_genbank": gbk_to_seqrecord(genbank).features[0].qualifiers["product"][0],
    }


def test_a_sequence_file_not_in_utf8_reads_in_the_locales_encoding(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(locale, "getpreferredencoding", lambda do_setlocale=True: "cp1252")

    assert _read_each_reader(tmp_path, "cp1252") == {
        "GenBank": TEXT,
        "FASTA": f"rec1 {TEXT}",
        "GFF3": TEXT,
        "comparison genome": "rec1",
        "crop_genbank": TEXT,
    }


def test_where_the_locale_is_utf8_a_file_not_in_utf8_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(locale, "getpreferredencoding", lambda do_setlocale=True: "UTF-8")

    with pytest.raises(ParseError) as refused:
        read_flat_file_records(_write(tmp_path, "rec1.gb", GENBANK, "cp1252"))
    assert refused.value.diagnostic == {"code": "INPUT_UNREADABLE", "reason": "GENBANK_PARSE"}


def test_a_utf8_file_reads_and_crop_genbank_writes_utf8_where_the_locale_is_not_utf8(tmp_path: Path) -> None:
    genbank = _write(tmp_path, "rec1.gb", GENBANK, "utf-8")
    cropped = tmp_path / "cropped.gb"
    script = (
        "import codecs, locale, sys\n"
        "from gbdraw.crop_genbank import gbk_to_seqrecord, write_cropped_genbank\n"
        "from gbdraw.io.genome import read_flat_file_records\n"
        "record = read_flat_file_records(sys.argv[1])[0]\n"
        "write_cropped_genbank(gbk_to_seqrecord(sys.argv[1]), sys.argv[2], 'rec1', 1, 12)\n"
        "print(ascii((codecs.lookup(locale.getpreferredencoding(False)).name, record.features[0].qualifiers['product'][0])))\n"
    )
    environment = {**os.environ, "LC_ALL": "C", "PYTHONUTF8": "0", "PYTHONCOERCECLOCALE": "0"}
    result = subprocess.run(
        [sys.executable, "-c", script, genbank, str(cropped)],
        cwd=PACKAGE.parent, env=environment, capture_output=True, text=True, encoding="utf-8", check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == ascii(("ascii", TEXT))
    assert f'/product="{TEXT}"' in cropped.read_bytes().decode("utf-8")
