"""A gzip-compressed input file reads as the plain file it holds, told by its first two bytes."""

from __future__ import annotations

import gzip
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

from gbdraw import cli, read_genbank, session_io
from gbdraw.analysis.depth import parse_depth_tsv
from gbdraw.exceptions import InputFileError, ParseError
from gbdraw.io import input_text
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.genome import load_gbks, parse_gff_fasta
from gbdraw.io.input_text import is_gzip_file, open_input_text, read_input_bytes
from gbdraw.session import load_session_document
from gbdraw.web_support.error_adapter import serialize_web_error

INPUTS = Path(__file__).resolve().parent / "test_inputs"
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
GENBANK = INPUTS / "NC_001416.gb"


def _gzipped(source: Path, target: Path) -> str:
    target.write_bytes(gzip.compress(source.read_bytes()))
    return str(target)


def _facts(records) -> list[tuple[str, int, int]]:
    return [(record.id, len(record.seq), len(record.features)) for record in records]


def test_a_gzip_genbank_file_reads_as_the_plain_one(tmp_path: Path) -> None:
    gzipped = _gzipped(GENBANK, tmp_path / "NC_001416.gb.gz")
    assert _facts(load_gbks([gzipped])) == _facts(load_gbks([str(GENBANK)]))
    assert _facts(read_genbank(gzipped)) == _facts(read_genbank(GENBANK))


def test_a_gzip_gff3_and_fasta_pair_reads_as_the_plain_pair(tmp_path: Path) -> None:
    gff = _gzipped(INPUTS / "NC_013668.gff3", tmp_path / "NC_013668.gff3.gz")
    fasta = _gzipped(INPUTS / "NC_013668.fasta", tmp_path / "NC_013668.fna.gz")
    plain = parse_gff_fasta(str(INPUTS / "NC_013668.gff3"), str(INPUTS / "NC_013668.fasta"))
    assert _facts(parse_gff_fasta(gff, fasta)) == _facts(plain)


def test_gzip_tables_read_as_the_plain_ones(tmp_path: Path) -> None:
    blast = EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out"
    assert read_comparison_table(_gzipped(blast, tmp_path / "blast.tsv.gz")).equals(read_comparison_table(str(blast)))
    depth = INPUTS / "MjeNMV.DRR271272.depth.tsv"
    gzip_depths = parse_depth_tsv(_gzipped(depth, tmp_path / "depth.tsv.gz")).depths
    plain_depths = parse_depth_tsv(str(depth)).depths
    assert list(gzip_depths) == list(plain_depths)
    assert all(np.array_equal(gzip_depths[name], plain_depths[name]) for name in plain_depths)


def test_the_bytes_tell_gzip_from_plain_not_the_name(tmp_path: Path) -> None:
    plain_named_gz = tmp_path / "plain.gb.gz"
    plain_named_gz.write_bytes(GENBANK.read_bytes())
    gzip_named_gb = Path(_gzipped(GENBANK, tmp_path / "compressed.gb"))
    assert (is_gzip_file(plain_named_gz), is_gzip_file(gzip_named_gb)) == (False, True)
    for path in (plain_named_gz, gzip_named_gb):
        assert read_input_bytes(path) == GENBANK.read_bytes()
        with open_input_text(path) as handle:
            assert handle.readline().startswith("LOCUS")


@pytest.mark.parametrize("damage", ["truncated", "corrupt"])
def test_a_gzip_genbank_file_cut_off_or_damaged_is_unreadable_input(tmp_path: Path, damage: str) -> None:
    compressed = gzip.compress(GENBANK.read_bytes())
    path = tmp_path / "PRIVATE.gb.gz"
    path.write_bytes(compressed[: len(compressed) // 2] if damage == "truncated" else b"\x1f\x8b\x08\x00" + bytes(range(128, 256)))
    with pytest.raises(ParseError) as caught:
        load_gbks([str(path)])
    assert caught.value.diagnostic == {"code": "INPUT_UNREADABLE", "reason": "GENBANK_PARSE"}


def _base64_resources(value, found=None) -> list[dict]:
    """Every Session resource descriptor (``encoding: base64``) anywhere in the document."""
    found = [] if found is None else found
    if isinstance(value, dict):
        if value.get("encoding") == "base64" and isinstance(value.get("data"), str):
            found.append(value)
        for member in value.values():
            _base64_resources(member, found)
    elif isinstance(value, list):
        for member in value:
            _base64_resources(member, found)
    return found


def test_the_cli_stores_a_gzip_input_unchanged_in_the_session_file(gbdraw_runner, tmp_path: Path) -> None:
    """A gzip input is stored as it is, without compressing it again (PD-OI-127)."""
    gzipped = Path(_gzipped(GENBANK, tmp_path / "NC_001416.gb.gz"))
    session_path = tmp_path / "gzip.gbdraw"
    code, output, svg_path = gbdraw_runner.run(
        "circular", [gzipped], "gzip", tmp_path, extra_args=["--session_output", str(session_path)]
    )
    assert code == 0, output
    assert svg_path.is_file()
    document = load_session_document(session_path)
    genbank = {rid: resource for rid, resource in document.to_dict()["resources"].items() if resource.get("kind") == "genbank"}
    assert len(genbank) == 1
    [(resource_id, resource)] = genbank.items()
    assert document.resources.read(resource_id) == gzipped.read_bytes()
    assert resource["name"] == "NC_001416.gb.gz"
    with zipfile.ZipFile(session_path) as archive:
        assert archive.getinfo(f"resources/{resource_id}").compress_type == zipfile.ZIP_STORED


# A gzip input expands to at most MAX_EXPANDED_INPUT_BYTES; past it every reader refuses the file by name
# (PD-OI-116 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md).
# The tests lower the cap, so their files stay small.
SMALL_CAP = 4096


@pytest.fixture
def small_cap(monkeypatch: pytest.MonkeyPatch) -> int:
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", SMALL_CAP)
    return SMALL_CAP


def _gzip_of(path: Path, content: bytes) -> Path:
    path.write_bytes(gzip.compress(content))
    return path


def _assert_refused(caught: pytest.ExceptionInfo[BaseException], path: Path) -> None:
    assert isinstance(caught.value, InputFileError)
    message = str(caught.value)
    assert str(path) in message
    assert f"{SMALL_CAP:,} bytes" in message
    assert "Decompress it" in message
    assert caught.value.diagnostic == {
        "code": "INPUT_UNREADABLE", "reason": "GZIP_EXPANDED_SIZE", "limit": 0, "inputName": path.name,
    }


def test_a_long_file_name_is_named_as_the_page_names_it(small_cap: int, tmp_path: Path) -> None:
    """PD-OI-116: Generate's refusal and the page's name the same characters of a file (OV-2911):
    a 120-character name whole, a longer one cut at the one bound the page reads too."""

    from gbdraw.web_support.error_adapter import INPUT_NAME_CHARACTERS

    schema_js = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "contracts" / "schema.generated.js"
    assert f"export const INPUT_NAME_CHARACTERS = {INPUT_NAME_CHARACTERS};" in schema_js.read_text(encoding="utf-8")
    for name, shown in (("n" * 113 + ".tsv.gz", None), ("m" * 233 + ".tsv.gz", "m" * INPUT_NAME_CHARACTERS)):
        path = _gzip_of(tmp_path / name, b"x" * (small_cap + 1))
        with pytest.raises(InputFileError) as caught:
            read_input_bytes(path)
        assert caught.value.diagnostic["inputName"] == name
        context = serialize_web_error(caught.value, operation="generate", stage="render")["context"]
        assert context["inputName"] == (shown or name)


def test_the_cap_is_2_gib() -> None:
    assert input_text.MAX_EXPANDED_INPUT_BYTES == 2_147_483_648


def test_a_gzip_input_past_the_cap_is_refused_by_name(small_cap: int, tmp_path: Path) -> None:
    path = _gzip_of(tmp_path / "PRIVATE-big.gb.gz", b"x" * (small_cap + 1))
    with pytest.raises(InputFileError) as caught:
        read_input_bytes(path)
    _assert_refused(caught, path)
    with pytest.raises(InputFileError) as caught, open_input_text(path) as handle:
        handle.read()
    _assert_refused(caught, path)


def test_a_gzip_input_of_exactly_the_cap_and_a_larger_plain_file_read(small_cap: int, tmp_path: Path) -> None:
    content = b"x" * small_cap
    at_cap = _gzip_of(tmp_path / "at-cap.gb.gz", content)
    assert read_input_bytes(at_cap) == content
    with open_input_text(at_cap) as handle:
        assert handle.read() == content.decode()
    plain = tmp_path / "plain.gb"
    plain.write_bytes(b"x" * (small_cap * 4))
    assert len(read_input_bytes(plain)) == small_cap * 4
    with open_input_text(plain) as handle:
        assert len(handle.read()) == small_cap * 4


def test_every_reader_lets_the_refusal_through_unchanged(small_cap: int, tmp_path: Path) -> None:
    """No reader rewords the refusal as a generic unreadable input."""

    readers = [
        (GENBANK, lambda path: load_gbks([str(path)])),
        (INPUTS / "MjeNMV.DRR271272.depth.tsv", lambda path: parse_depth_tsv(str(path))),
        (EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out", lambda path: read_comparison_table(str(path))),
        (INPUTS / "NC_013668.fasta", lambda path: parse_gff_fasta(str(INPUTS / "NC_013668.gff3"), str(path))),
    ]
    for source, read in readers:
        path = _gzip_of(tmp_path / f"{source.name}.gz", source.read_bytes())
        with pytest.raises(InputFileError) as caught:
            read(path)
        _assert_refused(caught, path)
        assert serialize_web_error(caught.value, operation="generate", stage="render")["context"] == {
            "reason": "GZIP_EXPANDED_SIZE", "limit": 0, "inputName": path.name,
        }


def test_a_gzip_session_is_not_capped(small_cap: int, tmp_path: Path) -> None:
    session = tmp_path / "big.gbdraw"
    session.write_bytes(gzip.compress(json.dumps({"padding": "x" * (small_cap * 4)}).encode()))
    assert len(session_io._read_session_text(session)) > small_cap * 4


def test_the_cli_names_the_refused_file_and_the_cap(
    small_cap: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    path = _gzipped(GENBANK, tmp_path / "big.gb.gz")
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", "--gbk", path, "-o", str(tmp_path / "out"), "-f", "svg"])
    with pytest.raises(SystemExit) as caught:
        cli.main()
    assert caught.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith(f"ERROR: {path}: ")
    assert f"{small_cap:,} bytes" in error


def test_the_web_helpers_read_a_gzip_fasta_as_the_plain_one(tmp_path: Path) -> None:
    """The page sends a chosen file to a Worker helper as chosen, so a helper reads a gzip FASTA as the CLI does."""

    helpers_js = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "app" / "python-helpers.js"
    helpers: dict[str, object] = {}
    exec(helpers_js.read_text(encoding="utf-8").split("`", 1)[1].rsplit("`", 1)[0], helpers)
    plain = str(INPUTS / "NC_013668.fasta")
    gzipped = _gzipped(INPUTS / "NC_013668.fasta", tmp_path / "NC_013668.fna.gz")
    for helper, args in (("extract_first_fasta", ("fasta",)), ("get_record_length", ("fasta",))):
        reply = json.loads(helpers[helper](gzipped, *args))
        assert "error" not in reply, (helper, reply)
        assert reply == json.loads(helpers[helper](plain, *args)), helper


def test_the_web_generate_refusal_names_the_file_the_page_read(small_cap: int, tmp_path: Path) -> None:
    """PD-OI-116: the page sends a gzip input as chosen and Python refuses it at Generate; the
    refusal names the file the user chose, not the file the Worker staged, and the cap."""

    import base64

    from gbdraw.api import CircularDiagramRequest, GenBankInputSource, RecordInput, RenderOutputRequest
    from gbdraw.session import build_session_document
    from gbdraw.web_support.request_render import _render_staged_canonical_web_request_with_prepared_inputs
    from tests.utils.session_files import embedded_session

    document = build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(GENBANK)),),
        output=RenderOutputRequest(output_prefix="big", formats=("svg",)),
    ))
    embedded = embedded_session(document.to_dict(), document.resources.sources)
    workspace = tmp_path / "gbdraw-web-render-1"
    (workspace / "resources").mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    paths, names = {}, {}
    for index, (resource_id, entry) in enumerate(embedded["resources"].items(), start=1):
        paths[resource_id] = workspace / "resources" / f"{index:064x}.bin"
        paths[resource_id].write_bytes(gzip.compress(base64.b64decode(entry["data"], validate=True)))
        names[resource_id] = "PRIVATE genome.gbff.gz"
    with pytest.raises(InputFileError) as caught:
        _render_staged_canonical_web_request_with_prepared_inputs(
            embedded["drawings"][0]["renderRequest"], resource_paths=paths, workspace=workspace, _resource_names=names,
        )
    assert serialize_web_error(caught.value, operation="generate", stage="render")["context"] == {
        "reason": "GZIP_EXPANDED_SIZE", "limit": 0, "inputName": "PRIVATE genome.gbff.gz",
    }


def test_a_cli_re_save_stores_a_gzip_blast_table_as_chosen(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A re-saved Session keeps a gzip BLAST table's gzip bytes, as the Web Save does (PD-OI-127)."""

    records = []
    for name in ("GZA", "GZB"):
        record = tmp_path / f"{name}.gb"
        record.write_text(GENBANK.read_text(encoding="utf-8").replace("NC_001416", name), encoding="utf-8")
        records.append(str(record))
    hits = _gzip_of(tmp_path / "hits.tsv.gz", b"GZA.1\tGZB.1\t99.0\t500\t5\t0\t1001\t1500\t2001\t2500\t1e-100\t900\n")

    def run(*arguments: str) -> None:
        monkeypatch.setattr(sys, "argv", ["gbdraw", *arguments])
        try:
            cli.main()
        except SystemExit as exit_:
            assert exit_.code in (None, 0)

    def blast_resource(session: Path) -> tuple[bytes, int]:
        document = load_session_document(session)
        [resource_id] = [
            item["resourceId"]
            for item in document.to_dict()["drawings"][0]["renderRequest"]["comparisons"]
            if item["kind"] == "nucleotideBlast"
        ]
        with zipfile.ZipFile(session) as archive:
            method = archive.getinfo(f"resources/{resource_id}").compress_type
        return document.resources.read(resource_id), method

    saved, resaved = tmp_path / "saved.gbdraw", tmp_path / "resaved.gbdraw"
    run("linear", "--gbk", *records, "-b", str(hits), "-o", str(tmp_path / "first"), "-f", "svg",
        "--session_output", str(saved))
    assert blast_resource(saved) == (hits.read_bytes(), zipfile.ZIP_STORED)
    run("render", "--session", str(saved), "-o", str(tmp_path / "second"), "-f", "svg", "--session_output", str(resaved))
    assert blast_resource(resaved) == (hits.read_bytes(), zipfile.ZIP_STORED)
    assert (tmp_path / "second.svg").read_bytes() == (tmp_path / "first.svg").read_bytes()


def test_restored_cli_tables_read_through_the_gzip_reader(small_cap: int, tmp_path: Path) -> None:
    """A restored CLI table of an older Session may be gzip; it reads as the plain one, under the cap."""

    track_table = _gzip_of(tmp_path / "tracks.tsv.gz", b"kind\tspacing\nfeatures\t5\n")
    session_io._migrate_legacy_circular_track_table(track_table)
    assert track_table.read_text(encoding="utf-8") == "kind\tparams\nfeatures\tlegacy_spacing=5\n".replace(
        "legacy_spacing", session_io._LEGACY_SPACING_PARAM
    )
    records_table = _gzip_of(tmp_path / "records.tsv.gz", b"gbk\trecord_label\nold.gb\tA\n")
    session_io._rewrite_tsv_path_cells(records_table, {(0, "gbk"): ["restored/new.gb"]})
    assert records_table.read_text(encoding="utf-8") == "gbk\trecord_label\nrestored/new.gb\tA\n"

    too_large = _gzip_of(tmp_path / "PRIVATE-tracks.tsv.gz", b"kind\tspacing\n" + b"x" * (small_cap + 1))
    for rewrite in (
        lambda: session_io._migrate_legacy_circular_track_table(too_large),
        lambda: session_io._rewrite_tsv_path_cells(too_large, {(0, "kind"): "y"}),
    ):
        with pytest.raises(InputFileError) as caught:
            rewrite()
        _assert_refused(caught, too_large)
