"""A Session the CLI saves reads back and saves again unchanged.

``gbdraw linear|circular --session_output S1`` saves a Session; ``gbdraw
render --session S1 --save_session`` reads it and saves it again. The two
Sessions hold the same render request and the same resources (identity, kind,
name and bytes). The requests cover every comparison kind the request contract
declares (``describe_canonical_request``), two input files of one name, and a
Circular conservation genome given as two files of one stem, whose LOSAT
results the Session stores under two names.
"""

from __future__ import annotations

import base64
from pathlib import Path
from typing import Callable, NamedTuple

import pytest
from Bio import SeqIO

from gbdraw.circular import circular_main
from gbdraw.cli_utils.session import render_main
from gbdraw.comparisons.losat_runtime import resolve_losat_runtime
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.session_request_codec import describe_canonical_request
from tests.utils.session_drawing_view import session_view
from tests.utils.session_files import read_session_payload

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
HUMAN = TUTORIAL / "human-mitochondrion" / "HmmtDNA.gbk"
# One genome as a FASTA file and as a GenBank file of the same stem.
DANIO_FASTA = TUTORIAL / "metazoan-mitochondria-comparison" / "NC_002333.2.fna"
DANIO_GENBANK = TUTORIAL / "metazoan-mitochondria-four" / "NC_002333.2.gb"


def _comparison_kinds() -> frozenset[str]:
    """The ``kind`` of each option of the request contract's comparison union."""

    types = describe_canonical_request()["types"]
    return frozenset(
        value
        for option in types["CanonicalComparison"]["options"]
        for value in types[option["name"]]["fields"]["kind"]["values"]
    )


class LinearRequest(NamedTuple):
    """Arguments of one ``gbdraw linear`` run and the comparison kinds its Session stores."""

    args: Callable[[Path], list[str]]
    kinds: frozenset[str]


def _pairwise_with_table(directory: Path) -> list[str]:
    """A BLAST table for records 1-2 and a listed LOSATP pair for records 2-3."""

    (directory / "pair.tsv").write_text("q\ts\t90\t40\t0\t0\t1\t40\t1\t40\t1e-20\t200\n", encoding="utf-8")
    table = directory / "comparisons.tsv"
    table.write_text("blast\tquery\tsubject\tsource\npair.tsv\t#1\t#2\ttable\n\t#2\t#3\tlosat\n", encoding="utf-8")
    return ["--comparisons_table", str(table), "--losat", "losatp", "--losatp_mode", "pairwise"]


def _blast_files(directory: Path) -> list[str]:
    """Two BLAST files (outfmt 7), stored as their raw rows."""

    paths = []
    for index in (1, 2):
        path = directory / f"pair-{index}.tsv"
        path.write_text(
            "# BLASTN 2.15.0+\n# Fields: query acc.ver, subject acc.ver, % identity\n"
            f"q{index}\ts{index}\t90.000\t40\t0\t0\t1\t40\t1\t40\t1.0e-20\t200\n",
            encoding="utf-8",
        )
        paths.append(str(path))
    return ["-b", *paths]


# One request per LOSATP display: Similarity groups and Collinear blocks each
# store a result the other cannot (the CLI refuses listed pairs with either), so
# no single request stores every kind. Each kind of the contract has a request
# that stores it.
LINEAR_REQUESTS = {
    "blast-files": LinearRequest(_blast_files, frozenset({"nucleotideBlast", "generatedProteinComparison"})),
    "table-and-pairwise": LinearRequest(
        _pairwise_with_table,
        frozenset({"nucleotideBlast", "precomputedProteinComparison", "generatedProteinComparison"}),
    ),
    "similarity-groups": LinearRequest(
        lambda _directory: ["--losat", "losatp", "--losatp_mode", "similarity_groups"],
        frozenset({"precomputedProteinComparison", "orthogroupResult", "generatedProteinComparison"}),
    ),
    "collinear": LinearRequest(
        lambda _directory: ["--losat", "losatp", "--losatp_mode", "collinear"],
        frozenset({"collinearityResult", "generatedProteinComparison"}),
    ),
}


def _stored(session: Path) -> tuple[dict[str, object], dict[str, tuple[str, str, bytes]]]:
    """The render request of the drawing ``session`` opens with, and its resources as ``{id: (kind, name, bytes)}``."""

    payload = read_session_payload(session)
    resources = {
        resource_id: (entry["kind"], entry["name"], base64.b64decode(entry["data"]))
        for resource_id, entry in payload["resources"].items()
    }
    return session_view(payload)["renderRequest"], resources


def _saved_again(saved: Path, directory: Path) -> Path:
    """``gbdraw render --session <saved> -o <directory>/again/fresh --save_session``: the Session it writes.

    The output prefix is the first save's (``fresh``), which the request holds.
    """

    again = directory / "again"
    again.mkdir()
    render_main(["--session", str(saved), "-o", str(again / "fresh"), "-f", "svg", "--save_session"])
    return again / "fresh.gbdraw"


def test_a_request_stores_each_comparison_kind_of_the_contract() -> None:
    """A new comparison kind needs a request here that stores it."""

    stored = frozenset().union(*(request.kinds for request in LINEAR_REQUESTS.values()))
    assert stored == _comparison_kinds()


@pytest.mark.linear
@pytest.mark.parametrize("name", sorted(LINEAR_REQUESTS))
def test_a_cli_linear_session_saved_again_holds_the_same_request_and_resources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from tests.test_render_stage_protein_analysis import _fake_losatp, _record

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    # Records 1 and 2 come from two files of one name.
    inputs = [tmp_path / "a" / "X.gbk", tmp_path / "b" / "X.gbk", tmp_path / "c" / "Y.gbk"]
    for index, path in enumerate(inputs):
        path.parent.mkdir()
        SeqIO.write(_record(index), path, "genbank")
    request = LINEAR_REQUESTS[name]
    saved = tmp_path / "saved.gbdraw"
    linear_main(
        ["--gbk", *map(str, inputs), *request.args(tmp_path), "-o", str(tmp_path / "fresh"), "-f", "svg",
         "--session_output", str(saved)]
    )

    render_request, resources = _stored(saved)
    assert {comparison["kind"] for comparison in render_request["comparisons"]} == request.kinds
    names = [name for _kind, name, _data in resources.values()]
    assert len(names) == len(set(names)), names
    assert sorted(name for kind, name, _data in resources.values() if kind == "genbank") == ["X.2.gbk", "X.gbk", "Y.gbk"]
    assert _stored(_saved_again(saved, tmp_path)) == (render_request, resources)


@pytest.fixture(scope="module")
def native_losat() -> None:
    try:
        runtime = resolve_losat_runtime("tlosatx")
    except ValidationError as error:
        pytest.skip(f"No native LOSAT runtime resolves: {error}")
    if runtime.kind != "losat":
        pytest.skip(f"Resolved runtime is {runtime.kind}, not native LOSAT.")


@pytest.mark.circular
@pytest.mark.usefixtures("native_losat")
def test_a_cli_circular_session_with_two_rings_of_one_stem_saved_again_holds_the_same_request_and_resources(
    tmp_path: Path,
) -> None:
    saved = tmp_path / "saved.gbdraw"
    circular_main(
        ["--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_gencode", "2", "--conservation_losat_gencode", "2",
         "--conservation_sequence", str(DANIO_FASTA), str(DANIO_GENBANK),
         "-o", str(tmp_path / "fresh"), "-f", "svg", "--session_output", str(saved)]
    )

    render_request, resources = _stored(saved)
    names = [name for _kind, name, _data in resources.values()]
    assert len(names) == len(set(names)), names
    assert {"NC_002333.2.fna", "NC_002333.2.gb"} <= set(names)
    assert _stored(_saved_again(saved, tmp_path)) == (render_request, resources)
