"""Warm equals cold for Circular requests (STEP2 design 3.2).

A render on a cache that holds the stage results of earlier requests must give
the same SVG text and reply metadata as a render on a new cache. The harness
runs the CLI entry and renders its typed request inside a cache transaction
(``tests/utils/render_stage_equivalence.py``).

- Matrix: a chain over the option variants of a case.
- Declarations: no undeclared settings read and no uncacheable stage.
- Identity: the per-record identity elements Circular inputs reach (source
  bytes, GFF3 type filter, position within an input, batch and Multi-Record
  Canvas), and a shared depth axis.
- Random walk: seeded single-option changes with returns to earlier states.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

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
    write_gff_fasta,
    write_tables,
)

pytestmark = pytest.mark.circular

_BASE = ("--labels", "out", "--gc", "--skew", "-l", "right")

# Chains: one variant per option on C1, named options on C2; the full universe is slow.
_CHAINS = {
    "C1-fast": ("C1", None, True),
    "C2-named": ("C2", ("label_table",), True),
    "C1-full": ("C1", None, False),
    "C2-full": ("C2", None, False),
    "C3-full": ("C3", None, False),
}


@pytest.fixture(scope="module")
def placeholders(tmp_path_factory: pytest.TempPathFactory) -> Placeholders:
    return Placeholders(write_tables(tmp_path_factory.mktemp("render-stage-tables")))


@pytest.fixture(
    scope="module",
    params=[
        "C1-fast",
        "C2-named",
        pytest.param("C1-full", marks=pytest.mark.slow),
        pytest.param("C2-full", marks=pytest.mark.slow),
        pytest.param("C3-full", marks=pytest.mark.slow),
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
    harness = EquivalenceHarness("circular", tmp_path_factory.mktemp(f"chain-{request.param}"))
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
        ("C1", 1, 8),
        pytest.param("C1", 11, 60, marks=pytest.mark.slow),
        pytest.param("C1", 12, 60, marks=pytest.mark.slow),
        pytest.param("C2", 13, 12, marks=pytest.mark.slow),
        pytest.param("C3", 14, 25, marks=pytest.mark.slow),
    ],
)
def test_random_walk_warm_equals_cold(
    case_id: str,
    seed: int,
    steps: int,
    placeholders: Placeholders,
    tmp_path: Path,
) -> None:
    harness = EquivalenceHarness("circular", tmp_path / "walk")
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


# --- Identity (M1, M5) -------------------------------------------------------------------------


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


def test_identity_source_content(tmp_path: Path) -> None:
    """The same path with new bytes is a new resource."""

    path = tmp_path / "a.gb"
    harness = EquivalenceHarness("circular", tmp_path / "out")
    arguments = ["--gbk", str(path), *_BASE]
    write_genbank(path, [synthetic_record("rec_a", seed=7, gc=0.3)])
    _chain_equivalence(harness, [("first content", arguments)])
    write_genbank(path, [synthetic_record("rec_a", seed=8, gc=0.7)])
    _chain_equivalence(harness, [("second content", arguments)])


def test_identity_gff3_type_filter(tmp_path: Path) -> None:
    gff, fasta = map(str, write_gff_fasta(tmp_path, synthetic_record("rec_g", seed=11)))
    harness = EquivalenceHarness("circular", tmp_path / "out")
    source = ["--gff", gff, "--fasta", fasta, *_BASE]
    _chain_equivalence(
        harness,
        [
            ("CDS", [*source, "-k", "CDS"]),
            ("CDS and tRNA", [*source, "-k", "CDS,tRNA"]),
            ("tRNA", [*source, "-k", "tRNA"]),
            ("CDS again", [*source, "-k", "CDS"]),
        ],
    )


@pytest.mark.parametrize("canvas", [(), ("--multi_record_canvas",)], ids=["batch", "multi-record-canvas"])
def test_identity_position_within_the_input(tmp_path: Path, canvas: tuple[str, ...]) -> None:
    """Two records of one file, same length, different composition, both drawn."""

    source = str(
        write_genbank(
            tmp_path / "two.gb",
            [synthetic_record("first", seed=5, gc=0.3), synthetic_record("second", seed=6, gc=0.7)],
        )
    )
    other = str(write_genbank(tmp_path / "other.gb", [synthetic_record("third", seed=12, gc=0.5)]))
    harness = EquivalenceHarness("circular", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("both records", ["--gbk", source, *canvas, *_BASE]),
            ("both after another input", ["--gbk", other, source, *canvas, *_BASE]),
            ("both records again", ["--gbk", source, *canvas, *_BASE]),
        ],
        uncached=True,
    )


def test_identity_shared_depth_axis_on_a_multi_record_canvas(tmp_path: Path) -> None:
    """M5: with a shared axis, a change to one record's depth rescales the other."""

    source = str(
        write_genbank(
            tmp_path / "two.gb",
            [synthetic_record("rec_a", seed=13, gc=0.4), synthetic_record("rec_b", seed=14, gc=0.6)],
        )
    )
    low = write_depth_table(tmp_path / "low.depth.tsv", [("rec_a", 3000), ("rec_b", 3000)])
    high_b = write_depth_table(tmp_path / "high_b.depth.tsv", [("rec_a", 3000)])
    with high_b.open("a", encoding="utf-8") as handle:  # rec_b about 15 times deeper
        handle.writelines(f"rec_b\t{position}\t{900 + position % 7}\n" for position in range(1, 3001))
    common = [
        "--gbk", source, "--multi_record_canvas", "--depth_window", "100", "--depth_step", "100",
        "--depth_track_label", "Depth", *_BASE,
    ]
    harness = EquivalenceHarness("circular", tmp_path / "out")
    _chain_equivalence(
        harness,
        [
            ("shared", [*common, "--depth_track", str(low), "--share_depth_axis"]),
            ("shared, b high", [*common, "--depth_track", str(high_b), "--share_depth_axis"]),
            ("own axes, b high", [*common, "--depth_track", str(high_b)]),
            ("shared again", [*common, "--depth_track", str(low), "--share_depth_axis"]),
        ],
    )
