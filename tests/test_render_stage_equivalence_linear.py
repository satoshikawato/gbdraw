"""Warm equals cold for Linear requests (STEP2 design 3.2).

A render on a cache that holds the stage results of earlier requests must give
the same SVG text and reply metadata as a render on a new cache. The harness
runs the CLI entry and renders its typed request inside a cache transaction
(``tests/utils/render_stage_equivalence.py``).

- Matrix: a chain over the option variants of a case.
- Declarations: no undeclared settings read and no uncacheable stage.
- Identity: each element of the per-record identity, a shared depth axis, and
  a similarity alignment reference (two builds in one transaction).
- Random walk: seeded single-option changes with returns to earlier states.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

import gbdraw.api.diagram as api_diagram_module
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    ProteinBlastpResult,
    ProteinExtractionResult,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from tests.utils.render_stage_equivalence import (
    CASES,
    EquivalenceHarness,
    Placeholders,
    StepResult,
    assert_same_outcome,
    declaration_failures,
    equivalence_failures,
    matrix_steps,
    random_walk,
    render_errors,
    run_chain,
    synthetic_record,
    write_depth_table,
    write_genbank,
    write_tables,
)

pytestmark = pytest.mark.linear

_BASE = ("--show_labels", "all", "--gc", "--skew", "-l", "right")

# Chains: one variant per option on L1, named options on L2; the full universe is slow.
_CHAINS = {
    "L1-fast": ("L1", None, True),
    "L2-named": (
        "L2",
        ("share_depth_axis", "evalue", "region (per record)", "reverse_complement (per file)"),
        True,
    ),
    "L1-full": ("L1", None, False),
    "L2-full": ("L2", None, False),
    "L3-full": ("L3", None, False),
}


@pytest.fixture(scope="module")
def placeholders(tmp_path_factory: pytest.TempPathFactory) -> Placeholders:
    return Placeholders(write_tables(tmp_path_factory.mktemp("render-stage-tables")))


@pytest.fixture(
    scope="module",
    params=[
        "L1-fast",
        "L2-named",
        pytest.param("L1-full", marks=pytest.mark.slow),
        pytest.param("L2-full", marks=pytest.mark.slow),
        pytest.param("L3-full", marks=pytest.mark.slow),
    ],
)
def chain(
    request: pytest.FixtureRequest,
    placeholders: Placeholders,
    tmp_path_factory: pytest.TempPathFactory,
) -> list[StepResult]:
    case_id, option_ids, fast = _CHAINS[request.param]
    case = CASES[case_id]
    steps = matrix_steps(
        case,
        placeholders,
        option_ids=option_ids,
        first_variant_only=fast,
        include_slow=not fast,
    )
    harness = EquivalenceHarness("linear", tmp_path_factory.mktemp(f"chain-{request.param}"))
    return run_chain(harness, steps)


def test_matrix_warm_equals_cold(chain: list[StepResult]) -> None:
    assert not render_errors(chain[:1]), "the case base must render"
    failures = equivalence_failures(chain)
    assert not failures, "\n".join(failures)


def test_declarations_have_no_undeclared_read_or_uncacheable_stage(
    chain: list[StepResult],
) -> None:
    failures = declaration_failures(chain)
    assert not failures, "\n".join(failures)


@pytest.mark.parametrize(
    ("case_id", "seed", "steps"),
    [
        ("L1", 1, 10),
        ("L1", 2, 10),
        pytest.param("L1", 11, 60, marks=pytest.mark.slow),
        pytest.param("L1", 12, 60, marks=pytest.mark.slow),
        pytest.param("L2", 13, 25, marks=pytest.mark.slow),
        pytest.param("L3", 14, 15, marks=pytest.mark.slow),
    ],
)
def test_random_walk_warm_equals_cold(
    case_id: str,
    seed: int,
    steps: int,
    placeholders: Placeholders,
    tmp_path: Path,
) -> None:
    harness = EquivalenceHarness("linear", tmp_path / "walk")
    results = random_walk(
        harness,
        CASES[case_id],
        placeholders,
        seed=seed,
        steps=steps,
        include_slow=steps > 20,
    )
    assert len(results) > steps // 2, "too many refused steps"
    failures = equivalence_failures(results)
    assert not failures, "\n".join(failures)


# --- Identity (M1, M5, M10) --------------------------------------------------------------------


def _chain_equivalence(
    harness: EquivalenceHarness,
    states: Sequence[tuple[str, Sequence[str]]],
    *,
    uncached: bool = False,
) -> None:
    """Render each state warm and cold (and without a cache) and compare them."""

    failures = []
    for label, arguments in states:
        warm, cold = harness.warm(arguments), harness.cold(arguments)
        assert cold.error is None, f"{label}: {cold.error}"
        for name, other in (("warm", warm), ("uncached", harness.uncached(arguments) if uncached else None)):
            if other is None:
                continue
            try:
                assert_same_outcome(other, cold, f"{label} ({name} vs cold)")
            except AssertionError as error:
                failures.append(str(error))
    assert not failures, "\n".join(failures)


@pytest.fixture
def two_files(tmp_path: Path) -> tuple[Path, Path]:
    return (
        write_genbank(tmp_path / "a.gb", [synthetic_record("rec_a", seed=1, gc=0.3)]),
        write_genbank(tmp_path / "b.gb", [synthetic_record("rec_b", seed=2, gc=0.7)]),
    )


def test_identity_input_index(tmp_path: Path, two_files: tuple[Path, Path]) -> None:
    a, b = map(str, two_files)
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("a b", ["--gbk", a, b, *_BASE]),
            ("b a", ["--gbk", b, a, *_BASE]),
            ("a b again", ["--gbk", a, b, *_BASE]),
        ],
    )


def test_identity_input_key(tmp_path: Path) -> None:
    source = str(
        write_genbank(
            tmp_path / "two.gb",
            [synthetic_record("first", seed=3, gc=0.3), synthetic_record("second", seed=4, gc=0.7)],
        )
    )
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("record #1", ["--gbk", source, "--record_id", "#1", *_BASE]),
            ("record #2", ["--gbk", source, "--record_id", "#2", *_BASE]),
            ("record #2 reversed", ["--gbk", source, "--record_id", "#2", "--reverse_complement", "true", *_BASE]),
            ("record #1 again", ["--gbk", source, "--record_id", "#1", *_BASE]),
        ],
    )


def test_identity_position_within_the_input(tmp_path: Path) -> None:
    """Two records of one file, same length, different composition, both drawn."""

    source = str(
        write_genbank(
            tmp_path / "two.gb",
            [synthetic_record("first", seed=5, gc=0.3), synthetic_record("second", seed=6, gc=0.7)],
        )
    )
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("record #2 alone", ["--gbk", source, "--record_id", "#2", *_BASE]),
            ("both records", ["--gbk", source, *_BASE]),
        ],
        uncached=True,
    )


def test_identity_source_content(tmp_path: Path) -> None:
    """The same path with new bytes is a new resource."""

    path = tmp_path / "a.gb"
    harness = EquivalenceHarness("linear", tmp_path / "out")
    arguments = ["--gbk", str(path), *_BASE]
    write_genbank(path, [synthetic_record("rec_a", seed=7, gc=0.3)])
    _chain_equivalence(harness, [("first content", arguments)])
    write_genbank(path, [synthetic_record("rec_a", seed=8, gc=0.7)])
    _chain_equivalence(harness, [("second content", arguments)])


def test_identity_collection_labels_and_regions(
    tmp_path: Path,
    two_files: tuple[Path, Path],
) -> None:
    a, b = map(str, two_files)
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("default collection", ["--gbk", a, b, *_BASE]),
            ("labels", ["--gbk", a, b, "--record_label", "A", "--record_label", "B", *_BASE]),
            ("other labels", ["--gbk", a, b, "--record_label", "A2", "--record_label", "B", *_BASE]),
            ("subtitles", ["--gbk", a, b, "--record_subtitle", "S", "--record_subtitle", "T", *_BASE]),
            ("region #2", ["--gbk", a, b, "--region", "#2:100-1500", *_BASE]),
            ("region #1", ["--gbk", a, b, "--region", "#1:100-1500", *_BASE]),
        ],
    )


def test_identity_collection_position_follows_earlier_inputs(
    tmp_path: Path,
    two_files: tuple[Path, Path],
) -> None:
    """Review finding #4: a collection region selects by global position.

    With ``--region #2:...`` the record of b.gb is #2 after a.gb (one record)
    and #3 after a two-record file, so it is cropped in the first request only.
    """

    a, b = map(str, two_files)
    pair = str(
        write_genbank(
            tmp_path / "pair.gb",
            [synthetic_record("rec_p", seed=9, gc=0.4), synthetic_record("rec_q", seed=10, gc=0.6)],
        )
    )
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("b is #2", ["--gbk", a, b, "--region", "#2:100-1500", *_BASE]),
            ("b is #3", ["--gbk", pair, b, "--region", "#2:100-1500", *_BASE]),
        ],
    )


def test_identity_collection_region_file_selector(
    tmp_path: Path,
    two_files: tuple[Path, Path],
) -> None:
    """OV-700 (review finding #5): the file selector of a collection region."""

    a, b = map(str, two_files)
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("crop a", ["--gbk", a, b, "--region", "a.gb:#1:100-1500", *_BASE]),
            ("crop b", ["--gbk", a, b, "--region", "b.gb:#1:100-1500", *_BASE]),
        ],
    )


def test_identity_shared_depth_axis(tmp_path: Path, two_files: tuple[Path, Path]) -> None:
    """M5: with a shared axis, a change to one record's depth rescales the other."""

    a, b = map(str, two_files)
    depth_a = str(write_depth_table(tmp_path / "a.depth.tsv", [("rec_a", 3000)]))
    depth_b = str(write_depth_table(tmp_path / "b.depth.tsv", [("rec_b", 3000)]))
    depth_b_high = str(write_depth_table(tmp_path / "b.high.depth.tsv", [("rec_b", 3000)], offset=600))
    common = ["--depth_window", "100", "--depth_step", "100", "--depth_track_label", "Depth", *_BASE]
    harness = EquivalenceHarness("linear", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("shared", ["--gbk", a, b, "--depth_track", depth_a, depth_b, "--share_depth_axis", *common]),
            ("shared, b high", ["--gbk", a, b, "--depth_track", depth_a, depth_b_high, "--share_depth_axis", *common]),
            ("own axes, b high", ["--gbk", a, b, "--depth_track", depth_a, depth_b_high, *common]),
            ("shared again", ["--gbk", a, b, "--depth_track", depth_a, depth_b, "--share_depth_axis", *common]),
        ],
    )


# Each record carries dnaA and gyrB; the groups stand in for a LOSATP analysis.
_ALIGNMENT_GENES = {"rec_x": (("dnaA", 120), ("gyrB", 900)), "rec_y": (("dnaA", 1500), ("gyrB", 300))}
_GROUPS = {
    f"{record}-{gene}": group
    for record in _ALIGNMENT_GENES
    for gene, group in (("dnaA", "og_1"), ("gyrB", "og_2"))
}


def _fixed_groups(
    records: Sequence[SeqRecord],
    *,
    protein_extraction: ProteinExtractionResult,
    **_kwargs: object,
) -> ProteinBlastpResult:
    groups: dict[str, list[OrthogroupMember]] = {}
    for proteins in protein_extraction.proteins_by_record:
        for protein in proteins:
            group = _GROUPS.get(str(protein.source_protein_id))
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
        comparisons=[DataFrame(columns=COMPARISON_COLUMNS) for _ in range(len(records) - 1)],
        orthogroups=OrthogroupResult(
            orthogroups=groups,
            member_by_protein_id={
                member.protein_id: member for members in groups.values() for member in members
            },
        ),
    )


@pytest.fixture
def fixed_groups(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    monkeypatch.setattr(
        api_diagram_module,
        "build_rbh_orthogroup_protein_blastp_comparisons",
        _fixed_groups,
    )
    yield


@pytest.mark.usefixtures("fixed_groups")
def test_identity_similarity_alignment_shares_stages_between_its_two_builds(tmp_path: Path) -> None:
    """M10: the analysis build and the drawn build share stage results in one transaction."""

    paths = [
        str(write_genbank(tmp_path / f"{name}.gb", [synthetic_record(name, seed=index + 20, genes=genes)]))
        for index, (name, genes) in enumerate(_ALIGNMENT_GENES.items())
    ]
    harness = EquivalenceHarness("linear", tmp_path / "out")
    arguments = [
        "--gbk", *paths,
        "--losat", "losatp",
        "--losatp_mode", "similarity_groups",
        "--similarity_alignment_feature", "rec_x-dnaA",
        *_BASE,
    ]
    _chain_equivalence(
        harness,
        [
            ("aligned on dnaA", arguments),
            ("aligned on gyrB", [*arguments[:-len(_BASE) - 1], "rec_x-gyrB", *_BASE]),
            ("aligned on dnaA again", arguments),
        ],
        uncached=True,
    )
    cold = harness.cold(arguments)
    assert cold.error is None
    # Both builds take the same records and windows: the second one hits.
    assert cold.stages.hits["gc_skew_data"] >= 1
