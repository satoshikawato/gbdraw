"""``SimilarityAlignmentReference``: one shared resolver for Python and the CLI.

A reference names one exact feature/protein ID. The typed core runs the
requested orthogroup analysis once, resolves the reference with
``resolve_similarity_alignment_plan`` and renders the resolved plan, so the
Python API, ``--similarity_alignment_feature`` and an equivalent hand-built plan
draw the same SVG. Fixed Similarity groups replace LOSATP.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

import gbdraw
import gbdraw.api.diagram as api_diagram_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    ProteinBlastpResult,
)
from gbdraw.api import (
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    LinearRecordTranslation,
    LosatSearchOptions,
    RecordInput,
    RenderOutputRequest,
    SimilarityAlignmentPlan,
    SimilarityAlignmentReference,
    build_request_plan_diagram,
    plan_request,
    render_request,
    resolve_similarity_alignment_plan,
)
from gbdraw.api.options import losatp_analysis_mode
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.layout.similarity_alignment import SimilarityAlignmentReferenceError
from gbdraw.session import build_session_document, materialize_session, session_to_request
from gbdraw.session_request_codec import (
    CanonicalRequestEncodingError,
    encode_canonical_request,
)

# Each record carries dnaA and gyrB at different positions, so alignment moves
# every row; the groups stand in for a completed LOSATP orthogroup analysis.
LAYOUT = {"a": (50, 300), "b": (200, 400), "c": (380, 100)}
SIMILARITY_GROUPS = LosatSearchOptions(
    program="losatp",
    losatp_mode="similarity_groups",
)
GROUPS = {
    f"{name}-{gene}": group
    for name in LAYOUT
    for gene, group in (("dnaA", "og_1"), ("gyrB", "og_2"))
}


def _record(name: str) -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 150), id=name, name=name, description=name)
    record.annotations.update(topology="linear", molecule_type="DNA")
    record.features = [
        SeqFeature(
            SimpleLocation(start, start + 60, strand=1),
            type="CDS",
            qualifiers={
                "gene": [gene],
                "protein_id": [f"{name}-{gene}"],
                "translation": ["M" * 20],
            },
        )
        for gene, start in zip(("dnaA", "gyrB"), LAYOUT[name], strict=True)
    ]
    return record


@pytest.fixture
def genbank_paths(tmp_path: Path) -> tuple[Path, ...]:
    paths = []
    for name in LAYOUT:
        path = tmp_path / f"{name}.gbk"
        SeqIO.write(_record(name), path, "genbank")
        paths.append(path)
    return tuple(paths)


@pytest.fixture
def analysis_calls(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """Replace the LOSATP orthogroup search with the fixed GROUPS."""

    calls: list[int] = []

    def fixed_groups(records, *, protein_extraction, **_kwargs):
        calls.append(len(records))
        groups: dict[str, list[OrthogroupMember]] = {}
        for proteins in protein_extraction.proteins_by_record:
            for protein in proteins:
                group = GROUPS.get(str(protein.source_protein_id))
                if group is None:
                    continue
                groups.setdefault(group, []).append(
                    OrthogroupMember(
                        orthogroup_id=group,
                        protein_id=protein.protein_id,
                        record_index=protein.record_index,
                        feature_index=protein.feature_index,
                        record_id=protein.record_id,
                        label=protein.label,
                        start=protein.start,
                        end=protein.end,
                        strand=protein.strand,
                        feature_svg_id=protein.feature_svg_id,
                        source_protein_id=protein.source_protein_id,
                    )
                )
        return ProteinBlastpResult(
            comparisons=[
                DataFrame(columns=COMPARISON_COLUMNS) for _ in range(len(records) - 1)
            ],
            orthogroups=OrthogroupResult(
                orthogroups=groups,
                member_by_protein_id={
                    member.protein_id: member
                    for members in groups.values()
                    for member in members
                },
            ),
        )

    monkeypatch.setattr(
        api_diagram_module,
        "build_rbh_orthogroup_protein_blastp_comparisons",
        fixed_groups,
    )
    return calls


def _cli_request(
    monkeypatch: pytest.MonkeyPatch,
    paths: tuple[Path, ...],
    output: Path,
    feature_id: str,
) -> LinearDiagramRequest:
    captured: list[LinearDiagramRequest] = []
    shared_render = linear_cli_module.render_request

    def spy(request, **kwargs):
        captured.append(request)
        return shared_render(request, **kwargs)

    monkeypatch.setattr(linear_cli_module, "render_request", spy)
    monkeypatch.setattr(
        "builtins.input",
        lambda *_args, **_kwargs: pytest.fail("Alignment must not prompt"),
    )
    linear_cli_module.linear_main(
        [
            "--gbk", *map(str, paths),
            "--losat", "losatp",
            "--losatp_mode", "similarity_groups",
            "--similarity_alignment_feature", feature_id,
            "-f", "svg",
            "-o", str(output),
        ]
    )
    assert len(captured) == 1
    return captured[0]


def _hand_built_plan(request: LinearDiagramRequest) -> SimilarityAlignmentPlan:
    """The plan a caller writes by hand: every dnaA anchor, record 1 first."""

    anchors = [
        AlignmentAnchorIdentity(
            record_key=item.record_key,
            biological_feature_id=item.source_feature_catalog[0].biological_feature_id,
            source_feature_index=item.source_feature_catalog[0].source_feature_index,
            stable_feature_svg_id=item.source_feature_catalog[0].stable_feature_id,
        )
        for item in plan_request(replace(request, similarity_alignment=None)).provenance
    ]
    return SimilarityAlignmentPlan(
        group_id="og_1",
        reference=anchors[0],
        records=(
            AlignmentRecordDecision(
                anchors[0].record_key,
                AlignmentDecisionStatus.REFERENCE,
                AlignmentResolutionRationale.REFERENCE,
                anchors[0],
            ),
            *(
                AlignmentRecordDecision(
                    anchor.record_key,
                    AlignmentDecisionStatus.ALIGNED,
                    AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
                    anchor,
                )
                for anchor in anchors[1:]
            ),
        ),
    )


def test_reference_is_typed_and_requires_the_orthogroup_analysis(genbank_paths) -> None:
    assert gbdraw.SimilarityAlignmentReference is SimilarityAlignmentReference
    assert SimilarityAlignmentReference(" a-dnaA ").feature_id == "a-dnaA"
    for invalid in ("", " ", "a\0b", 7):
        with pytest.raises(ValidationError, match="feature_id"):
            SimilarityAlignmentReference(invalid)  # type: ignore[arg-type]
    records = tuple(RecordInput(GenBankInputSource(path)) for path in genbank_paths)
    with pytest.raises(ValidationError, match="unsupported type"):
        LinearDiagramRequest(
            records=records,
            options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
            similarity_alignment="a-dnaA",  # type: ignore[arg-type]
        )
    for search in (
        None,
        LosatSearchOptions(program="losatp", losatp_mode="pairwise"),
        LosatSearchOptions(program="losatp", losatp_mode="collinear"),
    ):
        with pytest.raises(
            SimilarityAlignmentReferenceError,
            match=(
                r"requires the Similarity groups analysis .*losat_search program "
                r"'losatp' with losatp_mode 'similarity_groups'"
            ),
        ):
            LinearDiagramRequest(
                records=records,
                options=LinearDiagramOptions(losat_search=search),
                similarity_alignment=SimilarityAlignmentReference("a-dnaA"),
            )


def test_python_reference_renders_like_the_cli_flag_and_a_hand_built_plan(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    genbank_paths,
    analysis_calls,
) -> None:
    cli_request = _cli_request(monkeypatch, genbank_paths, tmp_path / "cli", "a-dnaA")
    assert cli_request.similarity_alignment == SimilarityAlignmentReference("a-dnaA")
    assert cli_request.options.losat_search is not None
    assert cli_request.options.losat_search.program == "losatp"
    assert cli_request.options.losat_search.losatp_mode == "similarity_groups"
    assert not hasattr(cli_request.options, "align_orthogroup_feature")
    assert not hasattr(cli_request.options, "similarity_alignment_feature")
    assert analysis_calls == [3]

    # Records without record_key take the planner's keys, as the CLI's do.
    python_request = LinearDiagramRequest(
        records=tuple(RecordInput(GenBankInputSource(path)) for path in genbank_paths),
        options=cli_request.options,
        layout=cli_request.layout,
        similarity_alignment=SimilarityAlignmentReference(feature_id="a-dnaA"),
        output=RenderOutputRequest(
            output_prefix="python", output_directory=tmp_path, formats=("svg",)
        ),
    )
    python_result = render_request(python_request)
    assert analysis_calls == [3, 3]  # one analysis per render, none to resolve

    hand_plan = _hand_built_plan(cli_request)
    resolved = python_result.request
    assert resolved.similarity_alignment == hand_plan
    assert losatp_analysis_mode(resolved.options.losat_search) == "none"
    assert resolved.layout.record_translations == tuple(
        LinearRecordTranslation(key) for key in ("record-1", "record-2", "record-3")
    )
    assert [
        decision.rationale for decision in hand_plan.records
    ] == [
        AlignmentResolutionRationale.REFERENCE,
        AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
        AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
    ]

    hand_request = replace(
        python_request,
        records=cli_request.records,
        layout=replace(
            cli_request.layout or LinearMultiRecordOptions(),
            record_translations=resolved.layout.record_translations,
        ),
        similarity_alignment=hand_plan,
        output=replace(python_request.output, output_prefix="hand"),
    )
    render_request(hand_request)
    assert analysis_calls == [3, 3, 3]

    cli_svg = (tmp_path / "cli.svg").read_bytes()
    assert b'data-record-translation-x="0.0"' in cli_svg
    assert cli_svg.count(b"data-record-translation-x=") == 3
    assert (tmp_path / "python.svg").read_bytes() == cli_svg
    assert (tmp_path / "hand.svg").read_bytes() == cli_svg


def test_beginner_reference_matches_the_beginner_plan(
    genbank_paths,
    analysis_calls,
) -> None:
    records = gbdraw.read_genbank(genbank_paths)

    def draw(alignment) -> str:
        return gbdraw.draw_linear(
            records,
            options=gbdraw.LinearOptions(
                comparisons=gbdraw.LinearComparisonOptions(
                    losat="losatp",
                    similarity_alignment=alignment,
                )
            ),
        )._drawing.tostring()

    from_reference = draw(SimilarityAlignmentReference("b-dnaA"))
    typed = LinearDiagramRequest(
        records=tuple(
            RecordInput(GenBankInputSource(path), record_key=f"record-{index}")
            for index, path in enumerate(genbank_paths, start=1)
        ),
        options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
    )
    plan = plan_request(typed)
    prepared = build_request_plan_diagram(plan)
    alignment = resolve_similarity_alignment_plan(
        plan,
        prepared.linear_metadata.orthogroups,
        SimilarityAlignmentReference("b-dnaA"),
    )
    assert alignment.reference.record_key == "record-2"
    assert from_reference == draw(alignment)
    assert from_reference != draw(None)


def test_reference_errors_are_neutral_and_the_cli_names_its_flag(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    genbank_paths,
    analysis_calls,
) -> None:
    def python_render(feature_id: str) -> None:
        render_request(
            LinearDiagramRequest(
                records=tuple(RecordInput(GenBankInputSource(path)) for path in genbank_paths),
                options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
                similarity_alignment=SimilarityAlignmentReference(feature_id),
                output=RenderOutputRequest(output_prefix="py", output_directory=tmp_path),
            )
        )

    cases = (
        ("og_1", "accepts an exact feature/protein ID, not Similarity Group ID 'og_1'."),
        ("missing", "did not match an exact feature/protein ID."),
    )
    for feature_id, detail in cases:
        with pytest.raises(SimilarityAlignmentReferenceError) as python_error:
            python_render(feature_id)
        assert str(python_error.value) == (
            f"SimilarityAlignmentReference(feature_id={feature_id!r}) {detail}"
        )
        assert "--similarity_alignment_feature" not in str(python_error.value)
        with pytest.raises(ValidationError) as cli_error:
            _cli_request(monkeypatch, genbank_paths, tmp_path / "cli", feature_id)
        assert str(cli_error.value) == f"--similarity_alignment_feature {detail}"
    assert not list(tmp_path.glob("*.svg"))

    # Two og_1 members in record 3 and no direct evidence: never prompt or guess.
    GROUPS_WITH_PARALOG = {**GROUPS, "c-gyrB": "og_1"}
    monkeypatch.setattr(f"{__name__}.GROUPS", GROUPS_WITH_PARALOG)
    with pytest.raises(SimilarityAlignmentReferenceError) as ambiguous:
        python_render("a-dnaA")
    assert ambiguous.value.detail.startswith(
        "cannot choose among multiple candidates; select an exact candidate for each record:"
    )
    assert "record 'record-3' candidates [" in ambiguous.value.detail
    assert "record-2" not in ambiguous.value.detail


def test_saved_requests_store_the_resolved_plan_never_the_reference(
    tmp_path: Path,
    genbank_paths,
    analysis_calls,
) -> None:
    request = LinearDiagramRequest(
        records=tuple(RecordInput(GenBankInputSource(path)) for path in genbank_paths),
        options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
        similarity_alignment=SimilarityAlignmentReference("a-dnaA"),
        output=RenderOutputRequest(output_prefix="saved", output_directory=tmp_path),
    )
    with pytest.raises(CanonicalRequestEncodingError, match="similarity alignment reference"):
        encode_canonical_request(request)

    document = build_session_document(request)
    assert analysis_calls == [3]
    payload = document.to_dict()
    assert payload["renderRequest"]["layout"]["similarityAlignment"]["groupId"] == "og_1"
    assert "SimilarityAlignmentReference" not in str(payload)
    with materialize_session(document, output_directory=tmp_path) as materialized:
        saved = session_to_request(materialized)
    assert isinstance(saved, LinearDiagramRequest)
    assert isinstance(saved.similarity_alignment, SimilarityAlignmentPlan)
    assert saved.similarity_alignment == _hand_built_plan(
        replace(
            request,
            records=tuple(
                replace(record, record_key=f"record-{index}")
                for index, record in enumerate(request.records, start=1)
            ),
        )
    )
    assert losatp_analysis_mode(saved.options.losat_search) == "none"


def test_existing_output_stops_the_render_before_the_analysis(
    tmp_path: Path,
    genbank_paths,
    analysis_calls,
) -> None:
    (tmp_path / "taken.svg").write_text("<svg/>", encoding="utf-8")
    with pytest.raises(ValidationError):
        render_request(
            LinearDiagramRequest(
                records=tuple(RecordInput(GenBankInputSource(path)) for path in genbank_paths),
                options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
                similarity_alignment=SimilarityAlignmentReference("a-dnaA"),
                output=RenderOutputRequest(output_prefix="taken", output_directory=tmp_path),
            )
        )
    assert analysis_calls == []


def _stored_tables(request: LinearDiagramRequest) -> list[DataFrame]:
    import base64
    from io import StringIO

    import pandas as pd

    payload = build_session_document(request).to_dict()
    return [
        pd.read_csv(StringIO(base64.b64decode(payload["resources"][item["resourceId"]]["data"]).decode()), sep="\t")
        for item in payload["renderRequest"]["comparisons"]
        if "resourceId" in item and item.get("encoding") == "canonicalTsv"
    ]


def _reversed_b(request: LinearDiagramRequest, reverse: bool = True) -> LinearDiagramRequest:
    from gbdraw.api import RecordPresentation

    records = tuple(
        replace(record, presentation=RecordPresentation(reverse_complement=reverse))
        if index == 1 else record
        for index, record in enumerate(request.records)
    )
    return replace(request, records=records)


def test_reference_resolution_keeps_given_comparison_rows_in_the_search_frame(genbank_paths, analysis_calls):
    # OV-399 (review 1): the resolved request copied the analysis build's
    # comparisons, which the build projects into the drawn frame.
    import pandas as pd

    from gbdraw.api.request_render import plan_linear_request
    from gbdraw.linear_comparison import LinearComparison

    inputs = tuple(RecordInput(GenBankInputSource(path), record_key=f"k{index}")
                   for index, path in enumerate(genbank_paths))
    plan = plan_linear_request(LinearDiagramRequest(records=inputs))
    row: dict[str, object] = {
        "query": "a", "subject": "b", "identity": 99.0, "alignment_length": 20, "mismatches": 0,
        "gap_opens": 0, "qstart": 51, "qend": 110, "sstart": 201, "send": 260, "evalue": 1e-20,
        "bitscore": 50.0,
    }
    for role, provenance in zip(("query", "subject"), plan.provenance[:2], strict=True):
        source = provenance.source_feature_catalog[0]
        row[f"{role}_feature_index"] = str(source.source_feature_index)
        row[f"{role}_feature_svg_id"] = source.stable_feature_id
    request = _reversed_b(LinearDiagramRequest(
        records=inputs,
        options=LinearDiagramOptions(
            losat_search=SIMILARITY_GROUPS,
            linear_comparisons=(LinearComparison(0, 1, pd.DataFrame([row])),),
        ),
        similarity_alignment=SimilarityAlignmentReference("a-dnaA"),
    ))
    (given,) = [table for table in _stored_tables(request) if "query" in table.columns and len(table)
                and table.iloc[0]["query"] == "a"]
    assert "subject_view_feature_svg_id" not in given.columns
    assert given.loc[0, ["qstart", "qend", "sstart", "send"]].tolist() == [51, 110, 201, 260]


def test_reference_resolution_stores_generated_rows_of_a_reversed_record_in_the_search_frame(
    genbank_paths, monkeypatch,
):
    # The analysis searches the drawn records; the resolved request stores its
    # rows in the search frame, so reversing record b moves no stored span.
    # (The stub reports each hit in the drawn orientation, so the direction of
    # a reversed endpoint's span flips.)
    import pandas as pd

    from gbdraw.analysis.protein_colinearity import ProteinBlastpResult as Result

    def groups_with_rows(records, *, protein_extraction, **_kwargs):
        members: dict[str, list[OrthogroupMember]] = {}
        dnaa = {}
        for proteins in protein_extraction.proteins_by_record:
            for protein in proteins:
                group = GROUPS.get(str(protein.source_protein_id))
                if group is None:
                    continue
                if group == "og_1":
                    dnaa[protein.record_index] = protein
                members.setdefault(group, []).append(OrthogroupMember(
                    orthogroup_id=group, protein_id=protein.protein_id, record_index=protein.record_index,
                    feature_index=protein.feature_index, record_id=protein.record_id, label=protein.label,
                    start=protein.start, end=protein.end, strand=protein.strand,
                    feature_svg_id=protein.feature_svg_id, source_protein_id=protein.source_protein_id,
                ))
        query, subject = dnaa[0], dnaa[1]
        row = {column: 0 for column in COMPARISON_COLUMNS}
        row.update(query=query.record_id, subject=subject.record_id, identity=99.0, alignment_length=20,
                   qstart=query.start + 1, qend=query.end, sstart=subject.start + 1, send=subject.end,
                   evalue=1e-20, bitscore=50.0)
        for role, protein in (("query", query), ("subject", subject)):
            row[f"{role}_feature_index"] = str(protein.feature_index)
            row[f"{role}_feature_svg_id"] = protein.feature_svg_id
            row[f"{role}_view_feature_svg_id"] = protein.view_feature_svg_id
        return Result(
            comparisons=[pd.DataFrame([row])] + [
                DataFrame(columns=COMPARISON_COLUMNS) for _ in range(len(records) - 2)
            ],
            orthogroups=OrthogroupResult(
                orthogroups=members,
                member_by_protein_id={m.protein_id: m for group in members.values() for m in group},
            ),
        )

    monkeypatch.setattr(api_diagram_module, "build_rbh_orthogroup_protein_blastp_comparisons", groups_with_rows)
    request = LinearDiagramRequest(
        records=tuple(RecordInput(GenBankInputSource(path), record_key=f"k{index}")
                      for index, path in enumerate(genbank_paths)),
        options=LinearDiagramOptions(losat_search=SIMILARITY_GROUPS),
        similarity_alignment=SimilarityAlignmentReference("a-dnaA"),
    )
    forward = _stored_tables(_reversed_b(request, reverse=False))
    reversed_ = _stored_tables(_reversed_b(request))
    assert len(forward[0]) == 1
    assert len(forward) == len(reversed_)
    for expected, actual in zip(forward, reversed_, strict=True):
        for start, end in (("qstart", "qend"), ("sstart", "send")):
            for table in (expected, actual) if len(expected) else ():
                table[[start, end]] = pd.DataFrame(
                    [sorted(pair) for pair in table[[start, end]].values.tolist()], index=table.index
                )
        pd.testing.assert_frame_equal(actual, expected)
