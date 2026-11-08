"""The locate owner: reference windows, target records and search rows to one plan.

Synthetic vectors pin each rule of ``gbdraw.comparisons.locate`` (design PR-L
section 2); real-data vectors pin the windows found with tblastn rows of the
fengycin cluster of B. velezensis FZB42 (NC_009725.2, windows sfp and
fengycin) in CBMB205 and SK19.001, and of the Vibrio plasmid pVA1
(NZ_CP065372.1) in 19-021-D1 and CAVP11. The fixtures keep only the rows that
pass the Linear thresholds, the reference genes, and the target CDS parts
under those rows.
"""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import random

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest

from gbdraw.comparisons.locate import (
    CdsExtent,
    LocateOptions,
    LocatePlan,
    LocateRecord,
    LocateReference,
    LocateSearchResult,
    LocateTarget,
    PlannedWindow,
    ReferenceCds,
    locate_default_thresholds,
    locate_query,
    locate_reference,
    locate_search_rows,
    locate_target,
    plan_locate,
    plan_locate_search,
)
from gbdraw.comparisons.losat_jobs import nucleotide_losat_cache_key
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import RecordWindow
from gbdraw.mode_profiles import LINEAR_MODE_PROFILE

FIXTURES = Path(__file__).parent / "fixtures" / "locate"


GENE_BP = 900
GENE_AA = 300
OPTIONS = LocateOptions.for_program()


# Synthetic inputs: reference gene k of window w lies at 100_000 * w + 1000 * k + 1 (900 bp, + strand).


def _reference(
    sizes: tuple[int, ...] = (10,),
    *,
    reverse: tuple[bool, ...] | None = None,
    circular: bool = False,
    extra: tuple[ReferenceCds, ...] = (),
) -> LocateReference:
    windows = []
    cds = list(extra)
    for window_index, size in enumerate(sizes):
        base = 100_000 * window_index
        windows.append(RecordWindow(base + 1, base + 1000 * size, bool(reverse and reverse[window_index])))
        cds.extend(
            ReferenceCds(f"g{window_index}_{k}", f"L{window_index}_{k}", None, None,
                         base + 1000 * k + 1, base + 1000 * k + GENE_BP, "+", GENE_AA, GENE_BP)
            for k in range(size)
        )
    return locate_reference(
        file_key="ref.gbk", record_key="ref", record_id="REF", record_length=100_000 * len(sizes),
        circular=circular, windows=windows, cds=cds,
    )


def _row(gene: str, record: str, start: int, *, strand: str = "+", bits: float = 500.0,
         qstart: int = 1, qend: int = GENE_AA, evalue: str = "1e-50", end: int | None = None) -> str:
    stop = start + GENE_BP - 1 if end is None else end
    sstart, send = (start, stop) if strand == "+" else (stop, start)
    return f"{gene}\t{record}\t90.0\t{qend - qstart + 1}\t0\t0\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bits}\n"


class _Target:
    """Rows and CDS of synthetic target records."""

    def __init__(self, *records: tuple[str, int, bool]) -> None:
        self.records = records
        self.rows: list[str] = []
        self.cds: dict[str, list[CdsExtent]] = {name: [] for name, _length, _circular in records}

    def place(self, genes: list[str], record: str, start: int, *, strand: str = "+", bits: float = 500.0,
              annotate: bool = True) -> _Target:
        for offset, gene in enumerate(genes):
            position = start + 1000 * offset
            self.rows.append(_row(gene, record, position, strand=strand, bits=bits))
            if annotate:
                self.cds[record].append(CdsExtent(position, position + GENE_BP - 1, strand))  # type: ignore[arg-type]
        return self

    def target(self, file_key: str = "t.gbk") -> LocateTarget:
        return LocateTarget(file_key, tuple(
            LocateRecord(name, name, index, length, circular, tuple(sorted(set(self.cds[name]), key=lambda c: (c.start, c.end))))
            for index, (name, length, circular) in enumerate(self.records)
        ))

    def plan(self, reference: LocateReference, options: LocateOptions = OPTIONS, rows: list[str] | None = None) -> LocatePlan:
        text = "".join(self.rows if rows is None else rows)
        return plan_locate(reference, [self.target()], [LocateSearchResult("t.gbk", text)], options)


def _genes(window: int, orders: range | list[int]) -> list[str]:
    return [f"g{window}_{k}" for k in orders]


def _spans(windows: tuple[PlannedWindow, ...]) -> list[tuple[str, int, int, bool, int | None]]:
    return [(w.record_id, w.window.start, w.window.end, w.window.reverse_complement, w.window.column) for w in windows]


def _codes(plan: LocatePlan) -> list[str]:
    return [warning.code for warning in plan.files[0].warnings]


# Reference.


def test_reference_orders_genes_along_the_drawn_windows_and_reports_cut_genes() -> None:
    cut = ReferenceCds("cut", None, None, None, 101_500, 102_500, "-", 100, 300)
    reference = _reference((3, 2), reverse=(False, True), extra=(cut,))

    assert [(gene.protein_id, gene.window_index, gene.order) for gene in reference.genes] == [
        ("g0_0", 0, 0), ("g0_1", 0, 1), ("g0_2", 0, 2), ("g1_1", 1, 3), ("g1_0", 1, 4),
    ]
    assert reference.cut_genes == ("cut",)
    plan = plan_locate(reference, [], [], OPTIONS)
    assert [(warning.code, warning.args) for warning in plan.warnings] == [("REFERENCE_GENE_CUT", (("count", 1),))]


def test_reference_windows_follow_their_columns() -> None:
    cds = [ReferenceCds(f"g{k}", None, None, None, 1000 * k + 1, 1000 * k + 900, "+", 300, 900) for k in range(4)]
    reference = locate_reference(
        file_key="r", record_key="r", record_id="R", record_length=10_000, circular=False,
        windows=[RecordWindow(2001, 4000, False, 2), RecordWindow(1, 2000, False, 1)], cds=cds,
    )
    assert [window.start for window in reference.windows] == [1, 2001]
    assert [gene.protein_id for gene in reference.genes] == ["g0", "g1", "g2", "g3"]


@pytest.mark.parametrize(("windows", "message"), [
    ([RecordWindow(1, 2000), RecordWindow(2000, 3000)], "overlap"),
    ([RecordWindow(5001, 5500)], "holds no whole CDS"),
    ([RecordWindow(1, 20_000)], "past the record end"),
    ([], "at least one window"),
])
def test_reference_windows_are_refused_when_they_cannot_be_searched(windows: list[RecordWindow], message: str) -> None:
    cds = [ReferenceCds("g", None, None, None, 1, 900, "+", 300, 900)]
    with pytest.raises(ValidationError, match=message):
        locate_reference(file_key="r", record_key="r", record_id="R", record_length=10_000, circular=False,
                         windows=windows, cds=cds)


def test_a_reference_window_may_cross_the_origin_of_a_circular_reference() -> None:
    cds = [
        ReferenceCds("a", None, None, None, 9101, 9900, "+", 266, 800),
        ReferenceCds("origin", None, None, None, 9951, 450, "+", 166, 500),
        ReferenceCds("b", None, None, None, 501, 1400, "-", 300, 900),
        ReferenceCds("cut", None, None, None, 1901, 2500, "+", 200, 600),
    ]
    reference = locate_reference(
        file_key="r", record_key="r", record_id="R", record_length=10_000, circular=True,
        windows=[RecordWindow(9001, 2000)], cds=cds,
    )
    assert [gene.protein_id for gene in reference.genes] == ["a", "origin", "b"]
    assert reference.cut_genes == ("cut",)
    with pytest.raises(ValidationError, match="crosses the origin of a linear record"):
        locate_reference(file_key="r", record_key="r", record_id="R", record_length=10_000, circular=False,
                         windows=[RecordWindow(9001, 2000)], cds=cds[:1])


def test_a_cds_across_the_origin_is_cut_by_a_window_that_ends_at_the_record_end() -> None:
    cds = [ReferenceCds("origin", None, None, None, 9951, 450, "+", 166, 500),
           ReferenceCds("a", None, None, None, 1001, 1900, "+", 300, 900)]
    reference = locate_reference(file_key="r", record_key="r", record_id="R", record_length=10_000, circular=True,
                                 windows=[RecordWindow(1, 10_000)], cds=cds)

    assert [gene.protein_id for gene in reference.genes] == ["a"]
    assert reference.cut_genes == ("origin",)


# Rule 1: filter.


def test_rows_below_the_thresholds_are_left_out_and_tallied() -> None:
    reference = _reference((4,))
    target = _Target(("c1", 200_000, False)).place(_genes(0, range(4)), "c1", 50_001)
    weak = [_row("g0_0", "c1", 150_001, bits=30.0), _row("g0_1", "c1", 151_001, evalue="1")]
    plan = target.plan(reference, rows=target.rows + weak)
    tally = plan.files[0].tally

    assert (tally.rows_read, tally.rows_kept) == (6, 4)
    assert {item.name: item.rows_left_out for item in tally.thresholds} == {"evalue": 1, "bitscore": 1}
    assert _spans(plan.files[0].windows) == [("c1", 49_001, 54_900, False, 1)]
    assert plan.files[0].candidates == ()


def test_a_file_whose_rows_all_fall_below_the_thresholds_has_no_window() -> None:
    reference = _reference((4,))
    plan = _Target(("c1", 200_000, False)).plan(reference, rows=[_row("g0_0", "c1", 1001, bits=20.0)])

    assert plan.files[0].windows == ()
    assert [(w.code, w.args) for w in plan.files[0].warnings] == [
        ("TARGET_NOT_ANNOTATED", ()),
        ("FILE_NO_WINDOWS", (("rows_read", 1), ("rows_kept", 0))),
    ]
    assert plan.files[0].coverage.absent == tuple(_genes(0, range(4)))


def test_default_thresholds_are_the_linear_comparison_defaults_for_both_programs() -> None:
    assert locate_default_thresholds("tlosatn") == LINEAR_MODE_PROFILE.comparison
    assert locate_default_thresholds("losatn") == LINEAR_MODE_PROFILE.comparison
    assert LocateOptions.for_program("losatn").thresholds == LINEAR_MODE_PROFILE.comparison


# Rule 2: one locus per gene.


def test_hits_of_one_gene_chain_into_one_locus_and_weaker_loci_are_not_anchors() -> None:
    reference = _reference((3,))
    target = _Target(("c1", 400_000, False)).place(_genes(0, [1, 2]), "c1", 51_001)
    rows = target.rows + [
        _row("g0_0", "c1", 50_001, qstart=1, qend=150, end=50_450, bits=300.0),
        _row("g0_0", "c1", 50_451, qstart=151, qend=GENE_AA, end=50_900, bits=300.0),
        # A repeated domain far away: 100 bits, below 0.8 of the gene's best locus (600).
        _row("g0_0", "c1", 300_001, qstart=1, qend=100, end=300_300, bits=100.0),
    ]
    plan = target.plan(reference, rows=rows)
    file_plan = plan.files[0]

    assert _spans(file_plan.windows) == [("c1", 49_001, 53_900, False, 1)]
    assert file_plan.candidates == ()
    gene = file_plan.windows[0].genes[0]
    assert (gene.protein_id, gene.start, gene.end, gene.query_coverage, gene.bits, gene.status) == (
        "g0_0", 50_001, 50_900, 1.0, 600.0, "complete",
    )


def test_equal_scores_keep_both_copies_and_the_later_one_is_a_candidate() -> None:
    reference = _reference((4,))
    target = (_Target(("c1", 200_000, False), ("c2", 200_000, False))
              .place(_genes(0, range(4)), "c2", 20_001)
              .place(_genes(0, range(4)), "c1", 80_001))
    plan = target.plan(reference)

    assert _spans(plan.files[0].windows) == [("c1", 79_001, 84_900, False, 1)]
    assert [(w.record_id, w.reason, w.selected) for w in plan.files[0].candidates] == [
        ("c2", "COVERED_BY_OTHER_WINDOW", False),
    ]


# Rule 3: whole CDS.


def test_an_anchor_takes_the_whole_cds_that_holds_it_but_not_a_neighbour_touching_it() -> None:
    reference = _reference((2,))
    target = _Target(("c1", 200_000, False)).place(["g0_1"], "c1", 61_001)
    rows = target.rows + [_row("g0_0", "c1", 60_101, end=60_800)]
    target.cds["c1"] += [CdsExtent(59_801, 60_900, "+"), CdsExtent(60_790, 60_990, "+")]
    plan = target.plan(reference, rows=rows)

    assert _spans(plan.files[0].windows) == [("c1", 58_801, 62_900, False, 1)]
    assert plan.files[0].windows[0].genes[0].start == 60_101


def test_a_target_without_annotation_keeps_the_aligned_span_and_warns() -> None:
    reference = _reference((3,))
    target = _Target(("c1", 200_000, False)).place(_genes(0, range(3)), "c1", 50_101, annotate=False)
    plan = target.plan(reference)

    assert _spans(plan.files[0].windows) == [("c1", 49_101, 54_000, False, 1)]
    assert _codes(plan) == ["TARGET_NOT_ANNOTATED"]


# Rule 4: clusters and the left-out interval.


@pytest.mark.parametrize(("gap", "windows"), [
    # ``gap`` bases between the last gene of one cluster and the first gene of the next.
    (9_999, [("c1", 49_001, 70_799, False, 1)]),
    (10_000, [("c1", 49_001, 70_800, False, 1)]),
    (10_001, [("c1", 49_001, 55_900, False, 1), ("c1", 63_902, 70_801, False, 2)]),
])
def test_clusters_of_one_record_merge_within_the_merge_distance(gap: int, windows: list) -> None:
    reference = _reference((10,))
    target = (_Target(("c1", 200_000, False))
              .place(_genes(0, range(5)), "c1", 50_001)
              .place(_genes(0, range(5, 10)), "c1", 54_900 + gap + 1))
    file_plan = target.plan(reference).files[0]

    assert _spans(file_plan.windows) == windows
    assert [w.gap_to_next_same_record for w in file_plan.windows] == ([63_902 - 55_900 - 1, None] if gap > 10_000 else [None])


# Rule 5: which clusters become windows.


def test_a_one_gene_reference_window_keeps_its_one_gene_cluster_but_a_large_window_does_not() -> None:
    reference = _reference((1, 38))
    target = (_Target(("c1", 3_000_000, False))
              .place(["g0_0"], "c1", 500_001)
              .place(_genes(1, range(38)), "c1", 1_000_001)
              .place(["g1_20"], "c1", 2_000_001))
    plan = target.plan(reference)
    file_plan = plan.files[0]

    assert [(w.window.start, len(w.genes)) for w in file_plan.windows] == [(499_001, 1), (999_001, 38)]
    assert [(w.window.start, w.reason) for w in file_plan.candidates] == [(1_999_001, "TOO_FEW_GENES")]
    assert file_plan.coverage.complete == tuple(gene.protein_id for gene in reference.genes)
    assert file_plan.windows[0].gap_to_next_same_record == 999_001 - 501_900 - 1


def test_max_windows_keeps_the_best_clusters_and_warns() -> None:
    reference = _reference((1, 10))
    target = (_Target(("c1", 3_000_000, False))
              .place(["g0_0"], "c1", 500_001)
              .place(_genes(1, range(10)), "c1", 1_000_001))
    plan = target.plan(reference, replace(OPTIONS, max_windows=1))
    file_plan = plan.files[0]

    assert [w.window.start for w in file_plan.windows] == [999_001]
    assert [(w.window.start, w.reason, w.window.column) for w in file_plan.candidates] == [(499_001, "MAX_WINDOWS", 1)]
    assert [(w.code, w.args) for w in file_plan.warnings] == [
        ("MAX_WINDOWS", (("max_windows", 1), ("clusters", 1))),
        ("GENES_ABSENT", (("count", 1),)),
    ]
    assert file_plan.coverage.absent == ("g0_0",)
    assert file_plan.coverage.reference_not_found == ((1, 900),)


def test_equal_clusters_rank_by_record_order_whatever_the_row_order() -> None:
    reference = _reference((5,))
    target = (_Target(("c1", 200_000, False), ("c2", 200_000, False))
              .place(_genes(0, range(5)), "c1", 30_001)
              .place(_genes(0, range(5)), "c2", 30_001))
    plans = [target.plan(reference, rows=rows) for rows in (target.rows, list(reversed(target.rows)))]

    assert plans[0].to_json() == plans[1].to_json()
    assert [w.record_id for w in plans[0].files[0].windows] == ["c1"]


# Rule 6: orientation and splits.


@pytest.mark.parametrize(("strands", "reverse"), [
    (("-", "-", "+"), True),
    (("-", "+"), False),
])
def test_a_window_is_reversed_when_opposite_hits_outweigh(strands: tuple[str, ...], reverse: bool) -> None:
    reference = _reference((len(strands),))
    target = _Target(("c1", 200_000, False))
    for k, strand in enumerate(strands):
        target.place([f"g0_{k}"], "c1", 50_001 + 1000 * k, strand=strand)
    plan = target.plan(reference)

    assert [w.window.reverse_complement for w in plan.files[0].windows] == [reverse]


def test_orientation_follows_the_reference_as_drawn() -> None:
    reference = _reference((3,), reverse=(True,))
    plan = _Target(("c1", 200_000, False)).place(_genes(0, range(3)), "c1", 50_001).plan(reference)

    assert [w.window.reverse_complement for w in plan.files[0].windows] == [True]


@pytest.mark.parametrize(("orders", "windows"), [
    # A fall of 4 genes in a 20-gene window stays in one window (the jump is max(5, 25 % of 20)).
    ([0, 1, 2, 3, 4, 9, 5, 6, 7, 8, *range(10, 20)], [(49_001, 70_900)]),
    # A fall of 5 breaks it; the two windows abut halfway between their genes, in reference order.
    ([0, 1, 2, 3, 4, 10, 5, 6, 7, 8, 9, *range(11, 20)], [(49_001, 55_950), (55_951, 70_900)]),
])
def test_a_cluster_splits_where_the_reference_order_falls_by_the_jump(orders: list[int], windows: list) -> None:
    reference = _reference((20,))
    target = _Target(("c1", 200_000, False)).place(_genes(0, orders), "c1", 50_001)
    plan = target.plan(reference)

    assert [(w.window.start, w.window.end) for w in sorted(plan.files[0].windows, key=lambda w: w.window.start)] == windows
    assert len({w.cluster for w in plan.files[0].windows}) == 1


def test_a_rotation_stored_in_a_linear_record_gives_two_abutting_windows_in_reference_order() -> None:
    reference = _reference((20,), circular=True)
    target = _Target(("c1", 30_000, False)).place(_genes(0, [*range(12, 20), *range(12)]), "c1", 1001)
    file_plan = target.plan(reference).files[0]

    assert _spans(file_plan.windows) == [("c1", 8951, 21_900, False, 1), ("c1", 1, 8950, False, 2)]
    assert [w.gap_to_next_same_record for w in file_plan.windows] == [0, None]
    assert file_plan.warnings == ()


def test_a_circular_record_stored_from_another_start_is_one_window_from_the_reference_start() -> None:
    reference = _reference((20,), circular=True)
    target = _Target(("c1", 21_000, True)).place(_genes(0, [*range(12, 20), *range(12)]), "c1", 1001)
    windows = target.plan(reference).files[0].windows

    # Drawn from base 8951 (gene 0) on, across the origin, to base 8950 (gene 19).
    assert _spans(windows) == [("c1", 8951, 8950, False, 1)]
    assert [(w.across_origin, w.touches_start, w.touches_end, len(w.genes)) for w in windows] == [(True, False, False, 20)]


def test_a_reversed_circular_record_stored_from_another_start_is_read_from_the_reference_start() -> None:
    reference = _reference((20,), circular=True)
    target = _Target(("c1", 21_000, True)).place(
        _genes(0, [*reversed(range(12)), *reversed(range(12, 20))]), "c1", 1001, strand="-"
    )
    windows = target.plan(reference).files[0].windows

    # Reversed, the window is drawn from its end: base 12_950 (gene 0) down, across the origin, to 12_951.
    assert _spans(windows) == [("c1", 12_951, 12_950, True, 1)]


def test_an_arc_across_the_origin_is_one_window() -> None:
    reference = _reference((10,))
    target = (_Target(("c1", 100_000, True))
              .place(_genes(0, range(5)), "c1", 96_001)
              .place(_genes(0, range(5, 10)), "c1", 1001))
    windows = target.plan(reference).files[0].windows

    assert _spans(windows) == [("c1", 95_001, 6900, False, 1)]
    assert [(w.across_origin, w.touches_start, w.touches_end) for w in windows] == [(True, False, False)]
    assert [g.order for g in windows[0].genes] == list(range(10))


def test_a_linear_record_does_not_join_its_ends() -> None:  # the circular case above
    reference = _reference((10,))
    target = (_Target(("c1", 100_000, False))
              .place(_genes(0, range(5)), "c1", 96_001)
              .place(_genes(0, range(5, 10)), "c1", 1001))
    windows = target.plan(reference).files[0].windows

    assert _spans(windows) == [("c1", 95_001, 100_000, False, 1), ("c1", 1, 6900, False, 2)]
    assert [w.across_origin for w in windows] == [False, False]
    assert windows[0].cluster != windows[1].cluster


# Rule 7: flank, record ends and columns.


@pytest.mark.parametrize(("start", "window_start"), [(2001, 1), (2002, 1002)])
def test_a_window_takes_the_record_end_when_at_most_the_flank_remains(start: int, window_start: int) -> None:
    reference = _reference((2,))
    plan = _Target(("c1", 200_000, False)).place(_genes(0, range(2)), "c1", start).plan(reference)
    window = plan.files[0].windows[0]

    assert (window.window.start, window.touches_start) == (window_start, window_start == 1)


def test_a_gene_at_a_contig_end_is_partial() -> None:
    reference = _reference((2,))
    target = _Target(("c1", 2500, False)).place(_genes(0, range(2)), "c1", 501)
    plan = target.plan(reference)
    file_plan = plan.files[0]

    assert _spans(file_plan.windows) == [("c1", 1, 2500, False, 1)]
    assert [g.status for g in file_plan.windows[0].genes] == ["complete", "complete"]
    edge = _Target(("c1", 2000, False)).place(_genes(0, range(2)), "c1", 1)
    edge_plan = edge.plan(reference)
    assert [g.status for g in edge_plan.files[0].windows[0].genes] == ["partial", "complete"]
    assert edge_plan.files[0].coverage.partial == ("g0_0",)
    assert [w.code for w in edge_plan.files[0].warnings] == ["GENES_PARTIAL"]


def test_windows_of_a_file_are_ordered_by_the_reference_order_of_their_genes() -> None:
    reference = _reference((12,))
    target = (_Target(("c1", 100_000, False), ("c2", 100_000, False), ("c3", 100_000, False))
              .place(_genes(0, range(8, 12)), "c1", 10_001)
              .place(_genes(0, range(4)), "c3", 10_001)
              .place(_genes(0, range(4, 8)), "c2", 10_001))
    windows = target.plan(reference).files[0].windows

    assert [(w.record_id, w.window.column, w.ref_span) for w in windows] == [
        ("c3", 1, (1, 3900)), ("c2", 2, (4001, 7900)), ("c1", 3, (8001, 11_900)),
    ]


# Plan.


def test_a_failed_search_is_reported_for_its_file_only() -> None:
    reference = _reference((2,))
    good = _Target(("c1", 100_000, False)).place(_genes(0, range(2)), "c1", 10_001)
    plan = plan_locate(
        reference,
        [good.target("a.gbk"), LocateTarget("b.gbk", ())],
        [LocateSearchResult("a.gbk", "".join(good.rows)), LocateSearchResult("b.gbk", "", error="LOSAT exited 2")],
        OPTIONS,
    )

    assert [(f.file_key, len(f.windows)) for f in plan.files] == [("a.gbk", 1), ("b.gbk", 0)]
    assert [(w.code, w.args) for w in plan.files[1].warnings] == [("FILE_SEARCH_FAILED", (("error", "LOSAT exited 2"),))]


def test_plan_needs_one_result_per_target_and_known_queries() -> None:
    reference = _reference((2,))
    target = _Target(("c1", 100_000, False)).place(_genes(0, range(2)), "c1", 10_001)
    with pytest.raises(ValidationError, match="one search result per target"):
        plan_locate(reference, [target.target()], [], OPTIONS)
    with pytest.raises(ValidationError, match="did not have"):
        target.plan(reference, rows=[_row("unknown", "c1", 1)])


@pytest.mark.parametrize(("field", "value"), [("best_ratio", 1.5), ("flank", -1), ("max_windows", 0), ("program", "blastn")])
def test_options_are_validated(field: str, value: object) -> None:
    with pytest.raises(ValidationError):
        replace(OPTIONS, **{field: value})


def test_the_plan_is_json_data() -> None:
    reference = _reference((2,))
    plan = _Target(("c1", 100_000, False)).place(_genes(0, range(2)), "c1", 10_001).plan(reference)
    data = json.loads(json.dumps(plan.to_json()))

    window = data["files"][0]["windows"][0]
    assert window["window"] == {"start": 9001, "end": 12_900, "reverse_complement": False, "column": 1}
    assert data["options"]["program"] == "tlosatn"
    assert data["reference"]["genes"][0]["cds"]["protein_id"] == "g0_0"
    assert data["files"][0]["tally"]["rows_read"] == 2


# Real data.


def _load(name: str) -> tuple[LocateReference, LocateTarget, str]:
    data = json.loads((FIXTURES / f"{name}.json").read_text(encoding="utf-8"))
    ref = data["reference"]
    windows = [RecordWindow(start, end, reverse) for start, end, reverse in ref["windows"] or [[1, ref["record_length"], False]]]
    cds = [
        ReferenceCds(tag, tag, gene or None, product or None, start, end, strand, length, end - start + 1)
        for tag, gene, start, end, strand, length, product in ref["genes"]
    ]
    reference = locate_reference(
        file_key="reference", record_key=ref["record_id"], record_id=ref["record_id"],
        record_length=ref["record_length"], circular=ref["circular"], windows=windows, cds=cds,
    )
    target = data["target"]
    records = tuple(
        LocateRecord(item["record_id"], item["record_id"], item["index"], item["length"], item["circular"],
                     tuple(CdsExtent(*extent) for extent in item["cds"]))
        for item in target["records"]
    )
    return reference, LocateTarget(target["file_key"], records), data["rows"]


def _real_plan(name: str, rows: str | None = None) -> LocatePlan:
    reference, target, text = _load(name)
    return plan_locate(reference, [target], [LocateSearchResult(target.file_key, text if rows is None else rows)], OPTIONS)


REAL_WINDOWS = {
    # CBMB205: intact, sfp 1.55 Mb away; both reversed.
    "fengycin_CBMB205": [
        ("NZ_CP011937.1", 3_559_634, 3_562_308, True, 1),
        ("NZ_CP011937.1", 1_924_443, 2_004_707, True, 2),
    ],
    # SK19.001: sfp on contig 13, fengycin over contigs 17, 34 (reversed), 39 and 32.
    "fengycin_SK19.001": [
        ("NZ_AOFO01000013.1", 149_360, 152_034, True, 1),
        ("NZ_AOFO01000017.1", 1, 27_082, True, 2),
        ("NZ_AOFO01000034.1", 12_684, 21_897, True, 3),
        ("NZ_AOFO01000039.1", 1, 13_732, False, 6),
        ("NZ_AOFO01000032.1", 1, 11_079, False, 7),
    ],
    # CBMB205 stored from inside the cluster: the arc across the origin is one window.
    "fengycin_CBMB205_rotated": [
        ("NZ_CP011937.1", 1_562_247, 1_564_921, True, 1),
        ("NZ_CP011937.1", 3_856_848, 7320, True, 2),
    ],
    # 19-021-D1: the circular plasmid is pVA1 stored from another start: one window from gene 0 on.
    "pva1_19-021-D1": [("NZ_CP046414.1", 31_654, 31_653, True, 4)],
    # CAVP11: contig 80 twice, around contigs 79 and 86.
    "pva1_CAVP11": [
        ("NZ_JABKAY010000080.1", 1, 13_298, True, 1),
        ("NZ_JABKAY010000079.1", 1, 23_609, False, 3),
        ("NZ_JABKAY010000086.1", 1, 27_527, False, 7),
        ("NZ_JABKAY010000080.1", 13_299, 30_248, True, 8),
    ],
}

REAL_CANDIDATES = {
    "fengycin_CBMB205": [],
    "fengycin_SK19.001": [("NZ_AOFO01000041.1", "TOO_FEW_GENES"), ("NZ_AOFO01000027.1", "TOO_FEW_GENES")],
    "fengycin_CBMB205_rotated": [],
    "pva1_19-021-D1": [
        ("NZ_CP046411.1", "TOO_FEW_GENES"), ("NZ_CP046412.1", "TOO_FEW_GENES"),
        ("NZ_CP046411.1", "TOO_FEW_GENES"), ("NZ_CP046413.1", "COVERED_BY_OTHER_WINDOW"),
    ],
    "pva1_CAVP11": [
        ("NZ_JABKAY010000094.1", "TOO_FEW_GENES"), ("NZ_JABKAY010000086.1", "TOO_FEW_GENES"),
        ("NZ_JABKAY010000044.1", "TOO_FEW_GENES"), ("NZ_JABKAY010000074.1", "COVERED_BY_OTHER_WINDOW"),
    ],
}

REAL_COVERAGE = {
    "fengycin_CBMB205": (39, 0, 0),
    "fengycin_SK19.001": (34, 5, 0),
    "fengycin_CBMB205_rotated": (39, 0, 0),
    "pva1_19-021-D1": (69, 2, 0),
    "pva1_CAVP11": (62, 6, 3),
}


# The files whose windows cross a circular record's origin.
REAL_CASES = sorted(REAL_WINDOWS)


@pytest.mark.parametrize("name", REAL_CASES)
def test_real_data_windows(name: str) -> None:
    file_plan = _real_plan(name).files[0]
    coverage = file_plan.coverage

    assert _spans(file_plan.windows) == REAL_WINDOWS[name]
    assert [(w.record_id, w.reason) for w in file_plan.candidates] == REAL_CANDIDATES[name]
    assert (len(coverage.complete), len(coverage.partial), len(coverage.absent)) == REAL_COVERAGE[name]


def test_real_data_partial_genes_of_sk19_are_the_long_nrps_genes_and_one_short_gene() -> None:
    coverage = _real_plan("fengycin_SK19.001").files[0].coverage

    assert coverage.partial == ("RBAM_RS09280", "RBAM_RS09285", "RBAM_RS09290", "RBAM_RS09295", "RBAM_RS09365")


@pytest.mark.parametrize("name", REAL_CASES)
def test_real_data_plans_do_not_depend_on_the_row_order(name: str) -> None:
    _reference_value, _target, rows = _load(name)
    lines = rows.splitlines(keepends=True)
    expected = _real_plan(name).to_json()
    for seed in range(3):
        shuffled = list(lines)
        random.Random(seed).shuffle(shuffled)
        assert _real_plan(name, "".join(shuffled)).to_json() == expected


# Search plan.


def _protein_record(record_id: str, sequence: str, cds: list[tuple[int, int, int, str]], *, circular: bool = False) -> SeqRecord:
    record = SeqRecord(Seq(sequence), id=record_id, annotations={"molecule_type": "DNA", "topology": "circular" if circular else "linear"})
    record.features = [
        SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS",
                   qualifiers={"locus_tag": [tag], "translation": ["MKV"]})
        for start, end, strand, tag in cds
    ]
    return record


REFERENCE_SEQUENCE = "ATGAAAGTTTAA" + "C" * 8 + "TTAAACTTTCAT" + "G" * 8


def test_locate_query_writes_protein_or_coding_sequence_queries() -> None:
    record = _protein_record("REF", REFERENCE_SEQUENCE, [(0, 12, 1, "A1"), (20, 32, -1, "A2")])
    windows = [RecordWindow(1, 40)]
    protein = locate_query(record, file_key="ref.gbk", record_key="ref", windows=windows, program="tlosatn")
    nucleotide = locate_query(record, file_key="ref.gbk", record_key="ref", windows=windows, program="losatn")

    ids = [gene.protein_id for gene in protein.reference.genes]
    assert [gene.cds.locus_tag for gene in protein.reference.genes] == ["A1", "A2"]
    assert [(gene.cds.start, gene.cds.end, gene.cds.strand, gene.cds.length_nt) for gene in protein.reference.genes] == [
        (1, 12, "+", 12), (21, 32, "-", 12),
    ]
    assert [line.split()[0] for line in protein.fasta.splitlines() if line.startswith(">")] == [f">{i}" for i in ids]
    assert "MKV" in protein.fasta
    assert nucleotide.fasta == f">{ids[0]}\nATGAAAGTTTAA\n>{ids[1]}\nATGAAAGTTTAA\n"


def test_search_plan_runs_one_job_per_target_file_and_restores_record_ids() -> None:
    reference = locate_query(
        _protein_record("REF", REFERENCE_SEQUENCE, [(0, 12, 1, "A1")]),
        file_key="ref.gbk", record_key="ref", windows=[RecordWindow(1, 40)], program="tlosatn",
    )
    one = _protein_record("one", "ACGT" * 10, [])
    two = _protein_record("two", "TTGCA" * 8, [])
    three = _protein_record("three", "GGCAT" * 8, [])
    plan = plan_locate_search(
        reference,
        [("a.gbk", [one, two], ["a#1", "a#2"]), ("b.gbk", [three], ["b#1"])],
        program="tlosatn",
    )

    assert [job.file_key for job in plan.jobs] == ["a.gbk", "b.gbk"]
    assert all(job.query_fasta == reference.fasta for job in plan.jobs)
    first, second = plan.jobs
    assert ">one" not in first.subject_fasta and first.subject_fasta.startswith(">n_")
    assert second.subject_fasta.startswith(">three\n")
    assert first.raw_key == nucleotide_losat_cache_key(
        program="tblastn", args=[], query_hash=first.batch.query.hash, subject_hash=first.batch.subject.hash,
        search_context=first.batch.search_context,
    )
    assert first.batch.search_context is not None and second.batch.search_context is None
    query_id = reference.reference.genes[0].protein_id
    search_ids = sorted(first.batch.subject.ids, key=lambda key: first.batch.subject.ids[key][0])
    text = "".join(f"{query_id}\t{search_id}\t90\t3\t0\t0\t1\t3\t1\t9\t1e-5\t60\n" for search_id in search_ids)
    rows = locate_search_rows(first, text)
    assert rows.file_key == "a.gbk"
    assert [line.split("\t")[1] for line in rows.text.splitlines()] == ["one", "two"]
    # A golden key: the raw-cache identity of a tblastn search stays stable and differs from tblastx's.
    assert second.raw_key == "d88767b7d8a9b97d8665984a25106624299d4d8dc730262a331315a487be42ae"
    assert second.raw_key != nucleotide_losat_cache_key(
        program="tblastx", args=[], query_hash=second.batch.query.hash, subject_hash=second.batch.subject.hash,
    )


def test_locate_target_reads_cds_parts_and_topology() -> None:
    record = _protein_record("c1", "A" * 100, [(10, 40, 1, "x"), (50, 80, -1, "y")], circular=True)
    target = locate_target("t.gbk", [record], ["t#1"])

    assert target == LocateTarget("t.gbk", (LocateRecord(
        "t#1", "c1", 0, 100, True, (CdsExtent(11, 40, "+"), CdsExtent(51, 80, "-")),
    ),))
