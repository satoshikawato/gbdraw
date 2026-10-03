"""Comparison-table coordinate frame (Web GUI audit 2026-09-30, D-18, PD-OI-073).

A comparison table is read in the search frame: the selected and cropped
sequence, 1-based, in the source strand. The planner projects a record's
reverse complement, and a row outside the cropped record is rejected
instead of drawn beyond the record (N-07, N-08).
"""
from __future__ import annotations

import random
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw import linear as linear_cli
from gbdraw.exceptions import GbdrawError

_RNG = random.Random(7)
_X, _Y, _Z = ("".join(_RNG.choice("ACGT") for _ in range(size)) for size in (3000, 1000, 2000))
_RIBBON = re.compile(r'<path [^>]*data-gbdraw-pairwise-match-id="[^"]*"[^>]*>')


def _write_record(directory: Path, record_id: str, sequence: str) -> Path:
    record = SeqRecord(Seq(sequence), id=record_id, name=record_id, description=f"{record_id} synthetic")
    record.annotations.update(molecule_type="DNA", topology="linear")
    record.features = [SeqFeature(SimpleLocation(0, len(sequence), strand=1), type="source")]
    path = directory / f"{record_id}.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _ribbon_spans(svg_text: str) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Return the query and subject x spans of each drawn comparison ribbon."""
    spans = []
    for element in _RIBBON.findall(svg_text):
        path = re.search(r'\sd="([^"]*)"', element).group(1)
        point = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", path)]
        spans.append((
            (round(min(point[0], point[2]), 1), round(max(point[0], point[2]), 1)),
            (round(min(point[4], point[6]), 1), round(max(point[4], point[6]), 1)),
        ))
    return spans


def _render(directory: Path, prefix: str, records: list[Path], row: str, *extra: str) -> list:
    table = directory / f"{prefix}.tsv"
    table.write_text(row + "\n", encoding="utf-8")
    linear_cli.linear_main([
        "--gbk", *map(str, records), "-b", str(table), "-f", "svg", "-o", str(directory / prefix), *extra,
    ])
    return _ribbon_spans((directory / f"{prefix}.svg").read_text(encoding="utf-8"))


@pytest.mark.linear
def test_reverse_complement_flag_matches_a_physically_reversed_record(tmp_path: Path) -> None:
    # R2 2001..3000 and R3 1..1000 are one shared block (search frame).
    r2 = _write_record(tmp_path, "R2", _X[:2000] + _Y)
    r3 = _write_record(tmp_path, "R3", _Y + _Z)
    r3_reversed = _write_record(tmp_path, "R3rc", str(Seq(_Y + _Z).reverse_complement()))
    flagged = _render(
        tmp_path, "flag", [r2, r3], "R2\tR3\t100\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847",
        "--reverse_complement", "0", "--reverse_complement", "1",
    )
    physical = _render(
        tmp_path, "physical", [r2, r3_reversed], "R2\tR3rc\t100\t1000\t0\t0\t2001\t3000\t3000\t2001\t0.0\t1847",
    )
    assert len(physical) == 1
    assert flagged == physical


@pytest.mark.linear
def test_row_outside_a_cropped_record_is_not_drawn_beyond_the_record(tmp_path: Path) -> None:
    sa = _write_record(tmp_path, "SA", _X)
    sb = _write_record(tmp_path, "SB", _X)
    # SB is cropped to 1000 bp, so subject 1101..1300 is outside the record.
    try:
        drawn = _render(
            tmp_path, "crop", [sa, sb], "SA\tSB\t95\t200\t10\t0\t101\t300\t1101\t1300\t1e-50\t300",
            "--region", "SB:1001-2000",
        )
    except (GbdrawError, SystemExit):
        return
    assert drawn == []


def _frame(*rows: tuple[int, int, int, int], bound: bool = False):
    import pandas as pd

    frame = pd.DataFrame(
        [("Q", "S", 99.0, 10, 0, 0, *row, 1e-20, 50.0) for row in rows],
        columns=["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
    )
    if bound:
        for role in ("query", "subject"):
            frame[f"{role}_feature_index"] = "0"
            frame[f"{role}_feature_svg_id"] = f"{role}-feature"
    return frame


def _records(reverse_subject: bool) -> list[SeqRecord]:
    from gbdraw.core.record_metadata import _write_coord_map

    query = SeqRecord(Seq(_X[:100]), id="Q")
    subject = SeqRecord(Seq(_Y[:80]), id="S")
    if reverse_subject:
        _write_coord_map(subject, base=80, step=-1)
    return [query, subject]


@pytest.mark.linear
def test_planner_projects_unbound_rows_and_keeps_bound_rows() -> None:
    from gbdraw.linear_comparison import LinearComparison, project_search_frame_comparisons

    (unbound,) = project_search_frame_comparisons(
        [LinearComparison(0, 1, _frame((1, 10, 5, 15)))], _records(reverse_subject=True)
    )
    assert unbound.matches.loc[0, ["qstart", "qend", "sstart", "send"]].tolist() == [1, 10, 76, 66]
    (bound,) = project_search_frame_comparisons(
        [LinearComparison(0, 1, _frame((1, 10, 5, 15), bound=True))], _records(reverse_subject=True)
    )
    assert bound.matches.loc[0, ["sstart", "send"]].tolist() == [5, 15]


@pytest.mark.linear
def test_row_outside_its_record_reports_a_comparison_diagnostic() -> None:
    from gbdraw.exceptions import ValidationError
    from gbdraw.linear_comparison import LinearComparison, project_search_frame_comparisons
    from gbdraw.web_support.error_adapter import serialize_web_error

    with pytest.raises(ValidationError, match=r"subject coordinates outside 1\.\.80") as caught:
        project_search_frame_comparisons(
            [LinearComparison(0, 1, _frame((1, 10, 70, 90)))], _records(reverse_subject=False)
        )
    error = serialize_web_error(caught.value, operation="generate", stage="render")
    assert (error["code"], error["context"]) == ("COMPARISON_INPUT", {"reason": "SEARCH_FRAME"})


_SESSIONS = Path(__file__).parent / "fixtures" / "sessions"


@pytest.mark.linear
def test_main_cli_session_with_a_reversed_record_draws_the_ribbons_main_drew(tmp_path: Path) -> None:
    # origin/main read the -b table after --reverse_complement, and its CLI
    # sidecar embeds that reverse-complemented record as the source, so the
    # stored rows already use the search frame of the embedded record.
    import gzip
    import json

    from gbdraw.session_io import load_session

    provenance = json.loads((_SESSIONS / "q-frame-main-linear-reverse.provenance.json").read_text(encoding="utf-8"))
    session = tmp_path / "main.gbdraw-session.json"
    session.write_bytes(gzip.decompress((_SESSIONS / "q-frame-main-linear-reverse.v42.gbdraw-session.json.gz").read_bytes()))
    sidecar = tmp_path / "resaved.gbdraw-session.json"
    linear_cli.linear_main([
        "--session", str(session), "-o", str(tmp_path / "replayed"), "-f", "svg", "--session_output", str(sidecar),
    ])
    expected = [tuple(tuple(side) for side in ribbon) for ribbon in provenance["mainRibbonSpans"]]
    assert _ribbon_spans((tmp_path / "replayed.svg").read_text(encoding="utf-8")) == expected
    # The re-saved Session is in the search frame and replays unchanged.
    assert load_session(sidecar)["version"] > 42
    linear_cli.linear_main(["--session", str(sidecar), "-o", str(tmp_path / "resaved"), "-f", "svg"])
    assert _ribbon_spans((tmp_path / "resaved.svg").read_text(encoding="utf-8")) == expected


@pytest.mark.linear
@pytest.mark.parametrize("scenario", ["upload", "losatn"])
def test_main_web_session_with_a_reversed_record_draws_the_ribbons_main_drew(tmp_path: Path, scenario: str) -> None:
    # origin/main's Web Save kept the source records with
    # presentation.reverseComplement and stored the reversed endpoint's rows
    # after the reverse complement; the reader converts them once (D-18).
    import gzip
    import json

    from gbdraw.session_io import load_session

    provenance = json.loads((_SESSIONS / "q-frame-main-web.provenance.json").read_text(encoding="utf-8"))
    session = tmp_path / "main.gbdraw-session.json"
    session.write_bytes(
        gzip.decompress((_SESSIONS / f"q-frame-main-web-{scenario}.v42.gbdraw-session.json.gz").read_bytes())
    )
    sidecar = tmp_path / "resaved.gbdraw-session.json"
    linear_cli.linear_main([
        "--session", str(session), "-o", str(tmp_path / "replayed"), "-f", "svg", "--session_output", str(sidecar),
    ])
    expected = [tuple(tuple(side) for side in ribbon) for ribbon in provenance["mainRibbonSpans"]]
    assert _ribbon_spans((tmp_path / "replayed.svg").read_text(encoding="utf-8")) == expected
    assert load_session(sidecar)["version"] > 42
    linear_cli.linear_main(["--session", str(sidecar), "-o", str(tmp_path / "resaved"), "-f", "svg"])
    assert _ribbon_spans((tmp_path / "resaved.svg").read_text(encoding="utf-8")) == expected


@pytest.mark.linear
@pytest.mark.parametrize(
    ("fixture", "stored_row"),
    [
        # Web Save with reversed R3c: the re-saved Session embeds R3c as its
        # reverse complement, so the converted rows are written in its frame.
        ("q-frame-main-web-upload.v42", "R2c\tR3c\t100.0\t1000\t0\t0\t2001\t3000\t3000\t2001\t0.0\t1847\n"),
        ("q-frame-main-web-losatn.v42", "R2c\tR3c\t100.0\t1000\t0\t0\t2001\t3000\t3000\t2001\t0.0\t1847\n"),
        # CLI -b sidecar without a reversed record: the table bytes are kept.
        ("se06-main-linear-blast-cli.v42", "R2c\tR3c\t100.000\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847\n"),
    ],
)
def test_resaved_session_keeps_nucleotide_blast_comparisons(tmp_path: Path, fixture: str, stored_row: str) -> None:
    # PR3-B2: a replayed Session re-saved nucleotideBlast items as
    # precomputedProteinComparison canonical TSV.
    import base64
    import gzip

    from gbdraw.session_io import load_session

    def comparison(session: Path) -> tuple[dict, bytes]:
        payload = load_session(session)
        item = payload["renderRequest"]["comparisons"][0]
        return item, base64.b64decode(payload["resources"][item["resourceId"]]["data"])

    def replay(session: Path, name: str) -> Path:
        linear_cli.linear_main(["--session", str(session), "-o", str(tmp_path / name), "-f", "svg",
                                "--session_output", str(tmp_path / f"{name}.gbdraw-session.json")])
        return tmp_path / f"{name}.gbdraw-session.json"

    main = tmp_path / "main.gbdraw-session.json"
    main.write_bytes(gzip.decompress((_SESSIONS / f"{fixture}.gbdraw-session.json.gz").read_bytes()))
    first = replay(main, "first")
    item, content = comparison(first)
    assert (item["kind"], item["queryRecordIndex"], item["subjectRecordIndex"]) == ("nucleotideBlast", 0, 1)
    assert content.decode("utf-8") == stored_row
    # The re-saved Session re-saves the same item and bytes and draws the same SVG.
    second = replay(first, "second")
    assert comparison(second) == (item, content)
    replay(second, "third")
    assert (tmp_path / "third.svg").read_bytes() == (tmp_path / "second.svg").read_bytes()


def test_table_text_rewrite_maps_reversed_endpoints_only() -> None:
    from gbdraw.linear_comparison import reverse_endpoint_table_text

    row = "R2c\tR3c\t100.0\t1000\t0\t0\t2001\t3000\t3000\t2001\t0.0\t1847\n"
    assert reverse_endpoint_table_text(row, (3000, False), (3000, True)) == (
        "R2c\tR3c\t100.0\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847\n"
    )
    with pytest.raises(GbdrawError, match="outside 1..2000"):
        reverse_endpoint_table_text(row, (3000, False), (2000, True))


@pytest.mark.linear
def test_cli_session_output_with_a_reversed_record_replays_the_cli_ribbons(tmp_path: Path) -> None:
    # The sidecar embeds the reverse-complemented record as a sequence, so the
    # -b rows (search frame of the source record) are written in its frame.
    r2 = _write_record(tmp_path, "R2", _X[:2000] + _Y)
    r3 = _write_record(tmp_path, "R3", _Y + _Z)
    sidecar = tmp_path / "fresh.gbdraw-session.json"
    fresh = _render(
        tmp_path, "fresh", [r2, r3], "R2\tR3\t100\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847",
        "--reverse_complement", "0", "--reverse_complement", "1", "--session_output", str(sidecar),
    )
    linear_cli.linear_main(["--session", str(sidecar), "-o", str(tmp_path / "replayed"), "-f", "svg"])
    replayed = (tmp_path / "replayed.svg").read_text(encoding="utf-8")
    assert len(fresh) == 1
    assert _ribbon_spans(replayed) == fresh
    # The replay draws the same SVG; only the record source attributes differ,
    # because the embedded reversed sequence is the replay's own source (D-22).
    source_attributes = re.compile(r' data-gbdraw-record-source-(?:start|end|step)="[^"]*"')
    assert source_attributes.sub("", replayed) == source_attributes.sub(
        "", (tmp_path / "fresh.svg").read_text(encoding="utf-8")
    )
