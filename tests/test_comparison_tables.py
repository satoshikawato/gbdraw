"""BLAST comparison tables: one reader for every entry point and record-ID binding.

CO-05 and N-03 (PD-OI-077): every reader keeps the first 12 typed columns and a
missing or malformed table fails instead of being skipped. CO-06 (PD-OI-074):
table IDs that contradict the endpoint records fail; unknown IDs warn.
"""

from __future__ import annotations

import json
import logging
import math
import random
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.analysis.conservation import load_conservation_sources
from gbdraw.api import LinearComparison, LinearDiagramOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.record_planning import resolve_linear_options
from gbdraw.exceptions import ComparisonIdentityError, ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS, load_comparisons
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    _decode_comparisons,
    decode_canonical_request,
)

_VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "comparison_outfmt6_table_cases.json").read_text(
        encoding="utf-8"
    )
)
_PERMISSIVE = SimpleNamespace(evalue=10.0, bitscore=0.0, identity=0.0, alignment_length=0)


def _table_text(lines: list[object]) -> str:
    return "".join(
        (line if isinstance(line, str) else "\t".join(line)) + "\n" for line in lines
    )


def _write_table(path: Path, lines: list[object]) -> Path:
    path.write_text(_table_text(lines), encoding="utf-8")
    return path


def _two_records() -> list[SeqRecord]:
    records = [SeqRecord(Seq("A" * 3000), id=record_id) for record_id in ("SA", "SB")]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


def _read_with_shared_reader(path: Path) -> pd.DataFrame:
    from gbdraw.io.comparisons import read_comparison_table

    return read_comparison_table(path)


def _read_with_load_comparisons(path: Path) -> pd.DataFrame:
    return load_comparisons([str(path)], _PERMISSIVE)[0]


def _read_with_canonical_codec(path: Path) -> pd.DataFrame:
    decoded = _decode_comparisons(
        [
            {
                "kind": "nucleotideBlast",
                "resourceId": "comparison-nucleotide-1",
                "queryRecordIndex": 0,
                "subjectRecordIndex": 1,
            }
        ],
        mode="linear",
        schema=CANONICAL_REQUEST_SCHEMA,
        resource_paths={"comparison-nucleotide-1": path},
    )
    return decoded["linear_comparisons"][0].matches


def _read_with_comparisons_table(path: Path) -> pd.DataFrame:
    manifest = path.parent / "comparisons.tsv"
    manifest.write_text(f"blast\tquery\tsubject\n{path.name}\t#1\t#2\n", encoding="utf-8")
    options = resolve_linear_options(
        LinearDiagramOptions(comparison_table_file=str(manifest)),
        records=_two_records(),
        layout=None,
    )
    return options.linear_comparisons[0].matches


def _read_with_conservation(path: Path) -> pd.DataFrame:
    result = load_conservation_sources(blast_config=_PERMISSIVE, conservation_files=[str(path)])
    if result.skipped_sources:
        # Similarity rings keep their source index and report the reader error.
        raise ValidationError(str(result.skipped_sources[0].skip_reason))
    return result.sources[0].dataframe


_ENTRY_POINTS = {
    "shared-reader": _read_with_shared_reader,
    "load_comparisons (-b)": _read_with_load_comparisons,
    "canonical codec (Web)": _read_with_canonical_codec,
    "--comparisons_table": _read_with_comparisons_table,
    "conservation ring": _read_with_conservation,
}


@pytest.mark.parametrize("entry_point", list(_ENTRY_POINTS))
@pytest.mark.parametrize(
    "case", _VECTORS["cases"], ids=[case["name"] for case in _VECTORS["cases"]]
)
def test_every_comparison_reader_follows_the_shared_outfmt_vectors(
    tmp_path: Path, case: dict, entry_point: str
) -> None:
    path = _write_table(tmp_path / "pair.tsv", case["lines"])
    read = _ENTRY_POINTS[entry_point]
    expect = case["expect"]
    if "error" in expect:
        with pytest.raises(ValidationError, match=expect["error"]):
            read(path)
        return

    frame = read(path)

    assert list(frame.columns[: len(COMPARISON_COLUMNS)]) == list(COMPARISON_COLUMNS)
    assert len(frame) == expect["rowCount"]
    assert isinstance(frame.index, pd.RangeIndex)
    if not expect["rowCount"]:
        return
    first = frame.loc[0, list(COMPARISON_COLUMNS)].tolist()
    assert first[:2] == expect["firstRow"][:2]
    for actual, expected in zip(first[2:], expect["firstRow"][2:], strict=True):
        assert math.isclose(float(actual), float(expected), rel_tol=1e-12)
    for column in ("alignment_length", "qstart", "qend", "sstart", "send"):
        assert pd.api.types.is_integer_dtype(frame[column])


def test_extra_columns_are_reported_in_an_info_log(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = next(case for case in _VECTORS["cases"] if case["name"] == "14 columns (std qlen slen)")
    path = _write_table(tmp_path / "pair.tsv", case["lines"])

    with caplog.at_level(logging.INFO, logger="gbdraw.io.comparisons"):
        _read_with_load_comparisons(path)

    assert "ignoring 2 extra column(s)" in caplog.text


def test_normalized_dataframe_uses_first_twelve_positional_columns() -> None:
    from gbdraw.io.comparisons import normalize_comparison_dataframe

    frame = pd.DataFrame(
        [["SA", "SB", "95.0", "200", "0", "0", "101", "300", "1101", "1300", "1e-50", "300", "3000"]]
    )

    normalized = normalize_comparison_dataframe(frame)

    assert normalized.loc[0, ["query", "identity", "qstart", "evalue"]].tolist() == [
        "SA",
        95.0,
        101,
        pytest.approx(1e-50),
    ]
    with pytest.raises(ValidationError, match=r"row 1, column 7 \(qstart\): 'x' is not an integer"):
        normalize_comparison_dataframe(frame.replace({"101": "x"}))
    boolean = frame.copy()
    boolean[2] = True
    with pytest.raises(ValidationError, match=r"row 1, column 3 \(identity\): 'True' is not a finite number"):
        normalize_comparison_dataframe(boolean)


def test_load_comparisons_rejects_a_missing_file_instead_of_shifting_later_tables(
    tmp_path: Path,
) -> None:
    valid = _write_table(tmp_path / "SB_SC.tsv", _VECTORS["cases"][0]["lines"])
    missing = tmp_path / "SA_SB.tsv"

    with pytest.raises(ValidationError, match=r"SA_SB\.tsv: the comparison file does not exist"):
        load_comparisons([str(missing), str(valid)], _PERMISSIVE)


def test_load_comparisons_rejects_a_malformed_file_instead_of_shifting_later_tables(
    tmp_path: Path,
) -> None:
    valid = _write_table(tmp_path / "SB_SC.tsv", _VECTORS["cases"][0]["lines"])
    malformed = _write_table(tmp_path / "SA_SB.tsv", [["SA", "SB", "95.0"]])

    with pytest.raises(ValidationError, match=r"SA_SB\.tsv: line 1: expected at least 12"):
        load_comparisons([str(malformed), str(valid)], _PERMISSIVE)


_COLUMN_MAPPING = re.compile(r"names\s*=\s*(?:list\(|tuple\()?\s*COMPARISON_COLUMNS")
# Python embedded in Web JavaScript still reading outfmt 6 on its own. Remove an
# entry when its reader goes; this list may only shrink.
_PENDING_EMBEDDED_READERS = {
    # Internal 12-column LOSAT output; P15 (CO-07) deletes
    # convert_losat_nucleotide_to_display_tsv.
    "web/js/app/python-helpers.js": 1,
}


def test_outfmt_table_columns_have_one_reader() -> None:
    """R9: only gbdraw.io.comparisons maps table columns onto COMPARISON_COLUMNS."""

    package_root = Path(gbdraw.__file__).parent
    sources = [*package_root.rglob("*.py"), *(package_root / "web" / "js").rglob("*.js")]
    readers = {}
    for path in sorted(sources):
        count = len(_COLUMN_MAPPING.findall(path.read_text(encoding="utf-8")))
        if count:
            readers[path.relative_to(package_root).as_posix()] = count

    assert readers == {"io/comparisons.py": 1, **_PENDING_EMBEDDED_READERS}


# --- CO-06: record-ID binding --------------------------------------------------


def _random_sequence(seed: int, length: int = 3000) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _pair_records(*ids: str) -> list[SeqRecord]:
    records = []
    for index, record_id in enumerate(ids):
        record = SeqRecord(Seq(_random_sequence(index)), id=record_id, name=record_id.split(".")[0])
        record.annotations["molecule_type"] = "DNA"
        records.append(record)
    return records


def _hits(query: str, subject: str) -> pd.DataFrame:
    row = [query, subject, 99.0, 1000, 0, 0, 2001, 3000, 1, 1000, 1e-100, 1800.0]
    return pd.DataFrame([row], columns=COMPARISON_COLUMNS)


def _render_svg(records: list[SeqRecord], comparisons: list[LinearComparison]) -> str:
    return assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            },
        ),
        linear_comparisons=comparisons,
        legend="none",
    ).tostring()


def test_table_matching_the_endpoint_records_draws_with_endpoint_metadata() -> None:
    svg = _render_svg(_pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("R2", "R3"))])

    assert 'data-query-record-id="R2"' in svg
    assert 'data-subject-record-id="R3"' in svg
    assert 'data-qstart="2001"' in svg


def test_version_suffix_differences_are_tolerated() -> None:
    svg = _render_svg(
        _pair_records("R2.1", "R3.2"), [LinearComparison(0, 1, _hits("R2.3", "R3"))]
    )

    assert 'data-query-record-id="R2.1"' in svg
    assert 'data-subject-record-id="R3.2"' in svg


def test_swapped_query_and_subject_table_is_rejected() -> None:
    with pytest.raises(
        ComparisonIdentityError, match=r"the query column names 'R3', which is the subject record"
    ) as error:
        _render_svg(_pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("R3", "R2"))])

    assert error.value.reason == "RECORD_ID"
    # The Web receives the comparison-identity code instead of an unclassified error.
    from gbdraw.web_support.error_adapter import serialize_web_error

    assert serialize_web_error(error.value, operation="generate", stage="render")["code"] == (
        "COMPARISON_IDENTITY"
    )


def test_table_naming_another_displayed_record_is_rejected() -> None:
    records = _pair_records("R1", "R2", "R3")
    with pytest.raises(ComparisonIdentityError, match=r"the subject column names 'R3', which is record #3 'R3'"):
        _render_svg(records, [LinearComparison(0, 1, _hits("R1", "R3"))])


def test_distinct_versions_of_one_accession_stay_distinct() -> None:
    records = _pair_records("X.1", "X.2", "Y")
    with pytest.raises(ComparisonIdentityError, match=r"the query column names 'X.1', which is record #1"):
        _render_svg(records, [LinearComparison(1, 2, _hits("X.1", "Y"))])


def test_unknown_table_ids_keep_positional_placement_with_a_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.WARNING, logger="gbdraw.linear_comparison"):
        svg = _render_svg(
            _pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("contig_A", "contig_B"))]
        )

    assert 'data-query-record-id="R2"' in svg
    assert 'data-subject-record-id="R3"' in svg
    assert 'data-query-record-index="0"' in svg
    assert "1 row(s) use sequence IDs that match no displayed record ('contig_A', 'contig_B')" in caplog.text


def _write_genbank(path: Path, record: SeqRecord) -> Path:
    SeqIO.write([record], path, "genbank")
    return path


def _write_hits(path: Path, frame: pd.DataFrame) -> Path:
    frame.to_csv(path, sep="\t", header=False, index=False)
    return path


def test_web_request_and_cli_reject_the_same_swapped_table(tmp_path: Path, gbdraw_runner) -> None:
    from gbdraw.api.request_render import render_request
    from gbdraw.api import InMemoryRecordSource, LinearDiagramRequest, RecordInput
    from gbdraw.session_request_codec import encode_canonical_request

    records = _pair_records("R2", "R3")
    swapped = _write_hits(tmp_path / "R3_R2.tsv", _hits("R3", "R2"))

    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=InMemoryRecordSource(record)) for record in records),
        options=LinearDiagramOptions(linear_comparisons=(LinearComparison(0, 1, _hits("R2", "R3")),)),
    )
    encoded = encode_canonical_request(request)
    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content)
        resource_paths[resource.resource_id] = target
    comparison = encoded.payload["comparisons"][0]
    comparison["kind"] = "nucleotideBlast"
    comparison.pop("encoding")
    resource_paths[comparison["resourceId"]] = swapped
    decoded = decode_canonical_request(
        encoded.payload, resource_paths=resource_paths, output_directory=tmp_path / "out"
    )
    with pytest.raises(ComparisonIdentityError, match=r"the query column names 'R3', which is the subject record") as web_error:
        render_request(decoded)

    gbk_files = [
        _write_genbank(tmp_path / f"{record.id}.gb", record) for record in records
    ]
    returncode, output, _svg = gbdraw_runner.run(
        "linear", gbk_files, "swapped", tmp_path, blast_files=[swapped]
    )

    assert returncode != 0
    assert str(web_error.value) in output


def test_cli_exits_non_zero_for_a_missing_comparison_file(tmp_path: Path, gbdraw_runner) -> None:
    records = _pair_records("SA", "SB", "SC")
    gbk_files = [_write_genbank(tmp_path / f"{record.id}.gb", record) for record in records]
    valid = _write_hits(tmp_path / "SB_SC.tsv", _hits("SB", "SC"))

    returncode, output, svg_path = gbdraw_runner.run(
        "linear", gbk_files, "missing", tmp_path, blast_files=[tmp_path / "SA_SB.tsv", valid]
    )

    assert returncode != 0
    assert re.search(r"SA_SB\.tsv: the comparison file does not exist", output)
    assert not svg_path.exists()
