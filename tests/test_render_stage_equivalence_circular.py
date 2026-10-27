"""Warm equals cold for Circular requests (STEP2 design 3.2).

A render on a cache that holds the stage results of earlier requests must give
the same SVG text and reply metadata as a render on a new cache. The harness
runs the CLI entry and renders its typed request inside a cache transaction
(``tests/utils/render_stage_equivalence.py``).

- Matrix: a chain over the option variants of a case. The full universe of a
  case is split into chunks, each a chain that starts from the case base, so
  each test fits the 300 s pytest timeout; they run with ``--render-stage-universe``.
- Declarations: no undeclared settings read and no uncacheable stage.
- Identity: the per-record identity elements Circular inputs reach (source
  bytes, GFF3 type filter, position within an input, batch and Multi-Record
  Canvas), and a shared depth axis.
- Random walk: seeded single-option changes with returns to earlier states;
  the slow walks are short and use many seeds, for the same timeout.
"""

from __future__ import annotations

import contextlib
import io
import re
from collections.abc import Callable, Sequence
from pathlib import Path

from typing import Any, cast

import pytest

import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli
from gbdraw.api.request_render import DiagramRequest
from tests.utils import render_stage_equivalence
from tests.utils.render_stage_equivalence import (
    UNIVERSE_MARKS,
    CASES,
    CapturedRequest,
    EquivalenceHarness,
    Placeholders,
    StepResult,
    assert_same_outcome,
    declaration_failures,
    equivalence_failures,
    first_difference,
    full_chain_chunks,
    matrix_steps,
    matrix_variants,
    random_walk,
    render_captured,
    run_chain,
    state_arguments,
    synthetic_record,
    write_depth_table,
    write_genbank,
    write_gff_fasta,
    write_tables,
)

pytestmark = pytest.mark.circular

_BASE = ("--labels", "out", "--gc", "--skew", "-l", "right")

# Chains: named options on C1 and C2; the first variant of each. C1-fast keeps the
# cheapest option of each stage a C1 P1/P2 option can rebuild (in universe order,
# which the chain relies on); every other C1 option runs in the C1 full chunks under
# ``--render-stage-universe``. C2-named takes about 20 s (a C2 render is about 6 s),
# so it runs with the universe.
_C1_FAST_OPTIONS = (
    "window", "arrow_head_length_ratio", "arrow_shaft_width_ratio", "axis_stroke_color",
    "line_stroke_width", "definition_font_size", "plot_title", "plot_title_font_size",
    "plot_title_position", "gc_content_mode", "gc_content_min_percent", "gc_content_tick_font_size", "show_gc_content_axis",
    "show_gc_content_ticks",
    "depth_track (add)", "legend", "legend_box_size", "legend_font_size", "separate_strands",
    "resolve_overlaps", "feature_overlap_tolerance_bp", "label_rendering", "label_whitelist", "label_blacklist",
    "qualifier_priority", "label_table",
    "feature_visibility_table", "feature_override_table", "annotation_table", "scale_interval",
    "hide_scale", "species (no plot title)", "keep_full_definition_with_plot_title",
)
_CHAINS = {
    "C1-fast": ("C1", _C1_FAST_OPTIONS),
    "C2-named": ("C2", ("label_table",)),
}
# The full universe of each case, in chunks within the time budget (``cost_chunks``).
_FULL_CHUNKS = full_chain_chunks(("C1", "C2", "C3"))


@pytest.fixture(scope="module")
def placeholders(tmp_path_factory: pytest.TempPathFactory) -> Placeholders:
    return Placeholders(write_tables(tmp_path_factory.mktemp("render-stage-tables")))


@pytest.fixture(
    scope="module",
    params=["C1-fast", pytest.param("C2-named", marks=UNIVERSE_MARKS), *(pytest.param(chunk_id, marks=UNIVERSE_MARKS) for chunk_id, chunk in _FULL_CHUNKS.items())],
)
def chain(
    request: pytest.FixtureRequest,
    placeholders: Placeholders,
    tmp_path_factory: pytest.TempPathFactory,
) -> list[StepResult]:
    if request.param in _FULL_CHUNKS:
        chunk = _FULL_CHUNKS[request.param]
        case = CASES[chunk.case_id]
        variants = matrix_variants(case)[chunk.start:chunk.stop]
    else:
        case_id, option_ids = _CHAINS[request.param]
        case = CASES[case_id]
        variants = matrix_variants(case, option_ids=option_ids, first_variant_only=True, include_slow=False)
    harness = EquivalenceHarness("circular", tmp_path_factory.mktemp(f"chain-{request.param}"))
    return run_chain(harness, matrix_steps(case, placeholders, variants))


def test_the_case_base_draws(chain: list[StepResult]) -> None:
    """A base that ends in an error checks only that warm and cold fail alike (OV-706)."""

    assert chain[0].cold.error is None, f"{chain[0].label}: {chain[0].cold.error}"


def test_the_harness_renders_without_the_label_fit_warning_search(monkeypatch: pytest.MonkeyPatch) -> None:
    """The search is bounded by wall-clock time, so its renders would depend on machine load (OV-705)."""

    seen: dict[str, object] = {}

    def stop(request: object, **kwargs: object) -> None:
        seen.update(kwargs)
        raise RuntimeError("stop")

    monkeypatch.setattr(render_stage_equivalence, "render_request", stop)
    captured = CapturedRequest("circular", cast(DiagramRequest, None), None, ())
    outcome = render_captured(captured, None)
    assert outcome.error == ("RuntimeError", "stop")
    assert seen.get("label_fit_warn_search_seconds") == 0.0


def _table_options(get_args: Callable[[list[str]], object]) -> set[str]:
    help_text = io.StringIO()
    with contextlib.redirect_stdout(help_text), pytest.raises(SystemExit):
        get_args(["--help"])
    return set(re.findall(r"--\w+_table\b", help_text.getvalue()))


def test_the_harness_stages_every_file_a_table_names(tmp_path: Path) -> None:
    """A Worker stages the files a table names, so a warm render can hit on them: every
    table option of either mode is a file table, staged by its reader, or a table of values."""

    assert _table_options(circular_cli._get_args) | _table_options(linear_cli._get_args) == (
        set(render_stage_equivalence.TABLE_FILE_READERS) | render_stage_equivalence.VALUE_TABLE_OPTIONS
    )
    tables = write_tables(tmp_path / "tables")
    gbk = (render_stage_equivalence.TEST_INPUTS / "HmmtDNA.gbk").resolve()
    hits = tmp_path / "hits.tsv"
    hits.write_text("")
    written = {
        "--records_table": ("gbk\n", gbk),
        "--comparisons_table": ("blast\tquery\tsubject\n", hits),
        "--conservation_table": ("blast\n", hits),
    }
    for option, (header, named) in written.items():
        table = tmp_path / f"{option[2:]}.tsv"
        cells = f"{named}\t1\t2" if option == "--comparisons_table" else str(named)
        table.write_text(f"{header}{cells}\n")
        assert render_stage_equivalence.input_files([option, str(table)]) == (table.resolve(), named.resolve())
    # A track table names no file today; its reader would list one if it did.
    track_table = tables / "tracks_hmmt.tsv"
    assert render_stage_equivalence.input_files(["--circular_track_table", str(track_table)]) == (
        track_table.resolve(),
    )


def test_matrix_warm_equals_cold(chain: list[StepResult]) -> None:
    failures = equivalence_failures(chain)
    assert not failures, "\n".join(failures)


def test_declarations_have_no_undeclared_read_or_uncacheable_stage(
    chain: list[StepResult],
) -> None:
    failures = declaration_failures(chain)
    assert not failures, "\n".join(failures)


@pytest.mark.parametrize(
    ("case_id", "option_ids"),
    [
        ("C1", ("annotation_table",)),
        ("C1", ("plot_title_position", "annotation_table")),
        ("C3", ("plot_title_position",)),
        # About 21 s (two C2 renders at load 10): it runs with the universe (D-16).
        pytest.param("C2", (), marks=UNIVERSE_MARKS),
    ],
)
def test_the_circular_catalog_item_of_the_contents_equals_the_catalog_of_the_svg(
    case_id: str,
    option_ids: tuple[str, ...],
    placeholders: Placeholders,
    tmp_path: Path,
) -> None:
    """When every top-level group is a content (the plot title too), a cache
    transaction reads the Circular catalog from the contents; without one it
    reads the SVG. Both give the same reply (annotations on C1 with and
    without its title; depth, conservation and an overlay legend on C3; the
    merged records of the Multi-Record Canvas on C2)."""

    case = CASES[case_id]
    variants = {option: variant for option, variant in matrix_variants(case, option_ids=option_ids)}
    arguments = state_arguments(case, variants, placeholders)
    harness = EquivalenceHarness("circular", tmp_path / "out")
    cold = harness.cold(arguments)
    assert cold.error is None
    assert cold.stages.builds["catalog_item"] == 1
    assert_same_outcome(harness.uncached(arguments), cold, f"{case_id} {option_ids} (uncached vs cold)")


# Slow walks: short walks over many seeds, each well within the 300 s timeout.
_SLOW_WALKS = {"C1": (range(11, 17), 20), "C2": (range(21, 25), 4), "C3": (range(31, 34), 10)}


@pytest.mark.parametrize(
    ("case_id", "seed", "steps"),
    [
        ("C1", 1, 8),
        *(
            pytest.param(case_id, seed, steps, marks=UNIVERSE_MARKS)
            for case_id, (seeds, steps) in _SLOW_WALKS.items()
            for seed in seeds
        ),
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
    )
    assert len(results) > steps // 2, "too many refused steps"
    failures = equivalence_failures(results)
    assert not failures, "\n".join(failures)


# --- Record stages (PR 4) ----------------------------------------------------------------------


_RADIAL_GROWTH = ("--label_placement", "radial", "--labels", "both", "--resolve_overlaps")


def test_radial_growth_replans_copies_of_the_feature_layers(
    placeholders: Placeholders, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The radial growth loop replans feature placements on copies: a stage result is never changed."""

    import gbdraw.diagrams.circular.assemble as assemble

    layouts: list[int] = []
    layout_features = assemble._radial_layout_feature_dict
    built: list[tuple[Any, dict[str, tuple[int, object]]]] = []
    create = assemble.create_feature_layers

    def counted(*args: Any, **kwargs: Any) -> Any:
        layouts.append(1)
        return layout_features(*args, **kwargs)

    def kept(*args: Any, **kwargs: Any) -> Any:
        result = create(*args, **kwargs)
        built.append((result, _placements(result)))
        return result

    monkeypatch.setattr(assemble, "_radial_layout_feature_dict", counted)
    monkeypatch.setattr(assemble, "create_feature_layers", kept)
    case = CASES["C1"]
    outcome = EquivalenceHarness("circular", tmp_path / "radial").uncached(
        [*placeholders.expand(case, case.base), *_RADIAL_GROWTH]
    )
    assert outcome.error is None, outcome.error
    assert len(layouts) == 2, "the case must grow the radius and replan once"
    [(layers, before)] = built
    assert _placements(layers) == before


def _placements(layers: Any) -> dict[str, tuple[int, object]]:
    return {key: (id(feature.placement), feature.feature_track_id) for key, feature in layers.foreground_features.items()}


def test_labels_after_a_radial_growth_reuse_the_feature_layers(placeholders: Placeholders, tmp_path: Path) -> None:
    case = CASES["C1"]
    harness = EquivalenceHarness("circular", tmp_path / "radial")
    first = harness.warm([*placeholders.expand(case, case.base), *_RADIAL_GROWTH])
    assert first.error is None, first.error
    assert first.stages.builds.get("feature_layers") == 1
    # Horizontal labels take the same feature layers and do not replan.
    horizontal = [*placeholders.expand(case, case.base), *_RADIAL_GROWTH[2:]]
    second = harness.warm(horizontal)
    assert second.stages.hits.get("feature_layers") == 1
    assert not second.stages.builds.get("feature_layers")
    assert_same_outcome(second, harness.cold(horizontal), "horizontal after radial")


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


# --- Record parts of a Multi-Record Canvas ----------------------------------------------------


def records_table(path: Path, rows: Sequence[tuple[str, int]]) -> list[str]:
    """A Multi-Record Canvas of HmmtDNA regions: (region, reverse complement) per row."""

    gbk = render_stage_equivalence.TEST_INPUTS / "HmmtDNA.gbk"
    path.write_text(
        "gbk\tregion\treverse_complement\n" + "".join(f"{gbk}\t{region}\t{flag}\n" for region, flag in rows)
    )
    return ["--records_table", str(path), "--multi_record_canvas"]


_C4_PARTS = ("--labels", "out", "--gc", "--skew")


def hatch_table(path: Path, second_color: str) -> list[str]:
    """One hatched highlight per record (``#1``, ``#2``); record 1's hatch is red."""

    path.write_text(
        "set_id\tid\tmark\trecord\tstart\tend\tlabel\thatch_color\n"
        "s1\ta1\thighlight\t#1\t1000\t5000\tA\t#cc0000\n"
        f"s1\ta2\thighlight\t#2\t2000\t6000\tB\t{second_color}\n"
    )
    return ["--annotation_table", str(path)]


def test_record_parts_of_two_records_whose_ids_collide_warm_equal_cold(tmp_path: Path) -> None:
    """C4: HmmtDNA whole and HmmtDNA 1-9000, so the second record renames the ids it
    shares with the first. Rows swap and change, one record is reverse-complemented,
    then each input of the record parts and of the projected drawer keys changes
    (label rendering, track preset, feature lanes, definition font, scale interval,
    stroke, custom stack, GC toggle, window). A hatch table gives both records the
    same hatch pattern, which record 2 renames, until record 2's hatch color
    changes. Last, HmmtDNA and lambda share no record id, so feature identities are
    not bound to records. A part whose key misses an input it reads would be a
    stale hit: warm would differ from cold. Each step draws new bytes."""

    whole_region = records_table(tmp_path / "a.tsv", (("1-16569", 0), ("1-9000", 0)))
    base = [*whole_region, *_C4_PARTS]
    steps: list[tuple[str, list[str]]] = [
        ("whole + region", base),
        ("rows swapped", [*records_table(tmp_path / "b.tsv", (("1-9000", 0), ("1-16569", 0))), *_C4_PARTS]),
        ("back", base),
        ("region row changed", [*records_table(tmp_path / "c.tsv", (("1-16569", 0), ("1-8000", 0))), *_C4_PARTS]),
        ("swapped again", [*records_table(tmp_path / "d.tsv", (("1-8000", 0), ("1-16569", 0))), *_C4_PARTS]),
        ("back again", base),
        ("record 1 reverse-complemented", [*records_table(tmp_path / "e.tsv", (("1-16569", 1), ("1-9000", 0))), *_C4_PARTS]),
        ("back from the reverse complement", base),
        ("label rendering", [*base, "--label_rendering", "embedded_only"]),
        ("track preset", [*base, "--track_type", "middle"]),
        ("feature lanes", [*base, "--feature_width", "120"]),
        ("definition font", [*base, "--definition_font_size", "10"]),
        ("scale interval", [*base, "--scale_interval", "4000"]),
        ("stroke", [*base, "--block_stroke_color", "red"]),
        ("custom stack", [*base, "--center_reserved_radius", "150"]),
        ("GC toggle", [*whole_region, "--labels", "out", "--no-gc", "--skew"]),
        ("window", [*base, "-w", "500"]),
        ("hatch shared by both records", [*base, *hatch_table(tmp_path / "h1.tsv", "#cc0000")]),
        ("record 2 hatch color", [*base, *hatch_table(tmp_path / "h2.tsv", "#0000cc")]),
        ("back to the base", base),
        (
            "HmmtDNA and lambda, identities not bound",
            [
                "--gbk",
                str(render_stage_equivalence.TEST_INPUTS / "HmmtDNA.gbk"),
                str(render_stage_equivalence.TEST_INPUTS / "NC_001416.gb"),
                "--multi_record_canvas",
                *_C4_PARTS,
            ],
        ),
    ]
    harness = EquivalenceHarness("circular", tmp_path / "out")
    failures = []
    previous = None
    for label, arguments in steps:
        warm, cold = harness.warm(arguments), harness.cold(arguments)
        assert cold.error is None, f"{label}: {cold.error}"
        if previous is not None and first_difference(previous.reply, warm.reply) is None:
            failures.append(f"{label}: draws the bytes of the step before it")
        previous = warm
        uncacheable = sorted(name for name, count in warm.stages.uncacheable.items() if count)
        if uncacheable:
            failures.append(f"{label}: uncacheable {uncacheable}")
        try:
            assert_same_outcome(warm, cold, label)
        except AssertionError as error:
            failures.append(str(error))
    assert not failures, "\n".join(failures)
