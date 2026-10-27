"""Expected builds (STEP2 design 3.2, PL4).

For every option of class P1 or P2 (the option changes only its own content,
P2 also the legend, the title or the composition plan), a warm render that goes
from the case base to one variant may build only the render stages that hold
the content the option changes. A stage that builds beyond that set is a lost
hit, so a later change that drops one fails here.

The classes are those of the DEPMAP measurement (``P1_P2_OPTIONS``). Each
variant starts from a cache that holds the base: the harness renders the base
again after each variant. The slow cases are split into chunks of variants,
each from the case base, so each test fits the 300 s pytest timeout.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from tests.utils.render_stage_equivalence import (
    UNIVERSE_MARKS,
    CASES,
    EquivalenceHarness,
    Mode,
    OptionSpec,
    Placeholders,
    Variant,
    base_arguments,
    cost_chunks,
    error_rule_failures,
    matrix_variants,
    state_arguments,
    write_tables,
)

# DEPMAP classes P1 and P2 per case (measure/depmap/classes.json, 2026-10-13).
_CIRCULAR_COMMON = (
    "arrow_head_length_ratio", "arrow_shaft_width_ratio", "axis_stroke_color", "axis_stroke_width",
    "plot_title", "plot_title_font_size", "gc_content_mode", "gc_content_min_percent",
    "gc_content_max_percent", "gc_content_large_tick_interval", "gc_content_small_tick_interval",
    "gc_content_tick_font_size", "show_gc_content_axis", "show_gc_content_ticks", "legend_box_size",
    "legend_font_size", "scale_interval",
)
_LINEAR_COMMON = (
    "arrow_head_length_ratio", "arrow_shaft_width_ratio", "axis_stroke_color", "definition_font_size",
    "plot_title", "plot_title_font_size", "plot_title_position", "gc_content_min_percent",
    "gc_content_max_percent", "gc_content_large_tick_interval", "gc_content_small_tick_interval",
    "gc_content_tick_font_size", "legend", "legend_box_size", "legend_font_size", "separate_strands",
    "resolve_overlaps", "feature_overlap_tolerance_bp", "label_table", "feature_override_table",
    "feature_placement_table", "annotation_table", "scale_interval", "hide_scale",
    "keep_definition_left_aligned", "definition_line_style", "record_label", "record_subtitle",
    "hide_accession", "hide_length", "track_axis_gap", "scale_style", "scale_stroke_color",
    "scale_stroke_width", "scale_font_size", "ruler_label_font_size", "ruler_label_color",
)
_DEPTH = (
    "depth_track_label", "depth_track_large_tick_interval", "depth_track_small_tick_interval",
    "depth_track_tick_font_size", "depth_window", "depth_step", "depth_min", "depth_max",
    "depth_log_scale", "show_depth_axis", "depth_large_tick_interval", "depth_small_tick_interval",
    "depth_tick_font_size",
)
P1_P2_OPTIONS: Mapping[str, frozenset[str]] = {
    "C1": frozenset((*_CIRCULAR_COMMON, "hide_scale", "species (no plot title)",
                     "keep_full_definition_with_plot_title")),
    "C2": frozenset((*_CIRCULAR_COMMON, "line_stroke_width", "definition_font_size", "plot_title_position",
                     "legend", "hide_scale", "multi_record_column_gap_ratio", "multi_record_row_gap_ratio")),
    "C3": frozenset((*_CIRCULAR_COMMON, *_DEPTH, "line_stroke_color", "line_stroke_width", "show_depth_ticks",
                     "feature_overlap_tolerance_bp", "strain (no plot title)", "conservation_labels",
                     "conservation_reference", "evalue", "bitscore", "alignment_length")),
    "L1": frozenset((*_LINEAR_COMMON, "axis_stroke_width", "line_stroke_width", "label_font_size",
                     "gc_content_mode", "show_gc_content_axis", "show_gc_content_ticks", "label_rendering",
                     "label_whitelist", "label_blacklist", "qualifier_priority", "feature_visibility_table",
                     "show_replicon", "label_placement", "linear_label_spacing", "show_labels",
                     "track_layout", "feature_height")),
    "L2": frozenset((*_LINEAR_COMMON, *_DEPTH, "line_stroke_color", "share_depth_axis", "align_center",
                     "evalue", "bitscore", "alignment_length", "pairwise_match_style", "comparison_height")),
    "L3": frozenset(("align_center (rows)", "track_layout (rows)")),
}

_LEGEND = frozenset({"legend_table", "legend_measurement", "legend_content"})
_LEGEND_LAYOUT = frozenset({"legend_measurement", "legend_content"})
_TITLE = frozenset({"plot_title_content", "circular_plot_title_content"})
_FEATURE_PLACEMENT = frozenset({"feature_layers", "feature_dom_index"})
_LABEL_FILTER = _FEATURE_PLACEMENT | {"label_filter_maps"}
_FEATURE_EDIT = _FEATURE_PLACEMENT | _LEGEND | {"catalog_normalization", "interactive_context"}
_DEPTH_DATA = frozenset({"depth_track_config", "depth_content"})
_DEPTH_FRAME = _DEPTH_DATA | {"depth_frame"}

# The stages that hold the content each option changes. An option not named builds nothing.
EXPECTED_BUILDS: Mapping[str, frozenset[str]] = {
    **dict.fromkeys(
        (
            "legend", "legend_box_size", "legend_font_size", "definition_font_size", "definition_line_style",
            "hide_accession", "hide_scale", "keep_definition_left_aligned", "ruler_label_font_size",
            "scale_font_size", "scale_style", "show_replicon", "align_center (rows)", "axis_stroke_width",
            "multi_record_column_gap_ratio", "multi_record_row_gap_ratio",
            # The legend wraps at the content width, which a tick label or a
            # line stroke changes. The definition column keeps its gap from the
            # GC% tick labels, so an option that changes them moves that width.
            "gc_content_tick_font_size", "line_stroke_width", "gc_content_mode", "gc_content_max_percent",
            "show_gc_content_axis", "show_gc_content_ticks",
        ),
        _LEGEND_LAYOUT,
    ),
    **dict.fromkeys(
        ("separate_strands", "resolve_overlaps", "show_labels", "track_layout", "track_axis_gap"),
        _FEATURE_PLACEMENT,
    ),
    # These also move the content bounds, at whose width the legend wraps.
    **dict.fromkeys(("feature_overlap_tolerance_bp", "track_layout (rows)"), _FEATURE_PLACEMENT | _LEGEND_LAYOUT),
    **dict.fromkeys(("label_table", "label_whitelist", "label_blacklist", "qualifier_priority"), _LABEL_FILTER),
    # A placement row is a feature edit: the stages that read the edited features rebuild.
    **dict.fromkeys(("feature_override_table", "feature_placement_table"), _FEATURE_EDIT),
    "feature_visibility_table": _FEATURE_EDIT | {"visibility_rules"},
    "annotation_table": _LEGEND | {"catalog_normalization", "interactive_context"},
    **dict.fromkeys(("depth_min", "depth_log_scale"), _DEPTH_FRAME),
    # The definition column keeps its gap from the depth tick labels, so these move the content width.
    **dict.fromkeys(("depth_window", "depth_step", "depth_max"), _DEPTH_FRAME | _LEGEND_LAYOUT),
    "show_depth_axis": _DEPTH_DATA | _LEGEND_LAYOUT,
    **dict.fromkeys(
        (
            "depth_track_large_tick_interval", "depth_track_small_tick_interval",
            "show_depth_ticks", "depth_large_tick_interval", "depth_small_tick_interval", "share_depth_axis",
        ),
        _DEPTH_DATA,
    ),
    **dict.fromkeys(("depth_track_tick_font_size", "depth_tick_font_size"), _DEPTH_DATA | _LEGEND_LAYOUT),
    # The Depth config is keyed by the style fields it reads, not the label (OV-1229).
    "depth_track_label": _LEGEND,
    "conservation_labels": _LEGEND,
    # The pairwise identity legend is computed from the filtered comparisons.
    **dict.fromkeys(("evalue", "bitscore", "alignment_length"), _LEGEND),
    **dict.fromkeys(("plot_title", "plot_title_font_size"), _TITLE),
}

# The options a record content reads, per mode: the Linear record draws the
# axis, the ruler, features and labels; the Circular record draws features and
# embedded labels (options that move its radial layout are listed below).
_RECORD_CONTENT: Mapping[Mode, tuple[str, frozenset[str]]] = {
    "linear": (
        "record_content",
        frozenset(
            (
                "keep_full_definition_with_plot_title",
                "arrow_head_length_ratio", "arrow_shaft_width_ratio", "axis_stroke_color", "axis_stroke_width",
                "line_stroke_color", "line_stroke_width", "feature_height", "feature_overlap_tolerance_bp",
                "feature_override_table", "feature_placement_table", "feature_visibility_table", "label_table",
                "label_whitelist", "label_blacklist", "qualifier_priority", "label_font_size", "label_placement",
                "label_rendering", "linear_label_spacing", "show_labels", "separate_strands", "resolve_overlaps",
                "track_axis_gap", "track_layout", "track_layout (rows)",
                # The axis ruler draws with the scale settings.
                "hide_scale", "scale_interval", "scale_style", "scale_stroke_color", "scale_stroke_width",
                "scale_font_size", "ruler_label_font_size", "ruler_label_color",
            )
        ),
    ),
    "circular": (
        "circular_record_content",
        frozenset(("arrow_head_length_ratio", "arrow_shaft_width_ratio", "line_stroke_color", "line_stroke_width")),
    ),
}

# The options Circular label placement reads: the arrow geometry the labels
# anchor on, the outer content radius of the radial layout (axis band) they
# clear, the definition text in the centre they avoid, and the ring labels
# their slots carry.
_CIRCULAR_LABEL_OPTIONS = (
    "arrow_head_length_ratio", "arrow_shaft_width_ratio", "axis_stroke_width", "species (no plot title)",
    "strain (no plot title)", "conservation_labels",
)

# The placed labels and the two contents that draw them.
_CIRCULAR_LABELS = frozenset({"circular_labels", "label_leaders_content", "label_text_content"})

# The options that move the Circular radial layout: the centre definition
# reserve (its text and font), the tick band, the depth axis footprint, and the
# feature lanes. Every stage that takes the layout's slots may rebuild.
_RADIAL_LAYOUT_OPTIONS = (
    "keep_full_definition_with_plot_title", "definition_font_size", "scale_interval", "hide_scale",
    "show_depth_axis", "show_depth_ticks", "depth_tick_font_size", "feature_overlap_tolerance_bp",
)
_RADIAL_LAYOUT_STAGES = _CIRCULAR_LABELS | {
    "circular_record_content", "ticks_content", "gc_content_content", "gc_skew_content", "depth_content",
    "conservation_content", "annotation_content",
}

# The per-record Circular contents a Multi-Record merge takes: it rebuilds with any of them.
_RECORD_MERGE_INPUTS = frozenset(
    {
        "axis_content", "circular_record_content", "ticks_content", "gc_content_content", "gc_skew_content",
        "depth_content", "conservation_content", "annotation_content", "label_leaders_content", "label_text_content",
        "circular_definition_content",
    }
)

# The group contents an option changes, besides the stages above.
CONTENT_BUILDS: Mapping[str, frozenset[str]] = {
    "annotation_table": frozenset({"annotation_content"}),
    # The Circular axis circle.
    **dict.fromkeys(("axis_stroke_color", "axis_stroke_width"), frozenset({"axis_content"})),
    # The ribbons: their geometry and style, and the rows the thresholds keep.
    **dict.fromkeys(
        ("align_center", "align_center (rows)", "pairwise_match_style", "comparison_height"),
        frozenset({"comparison_content"}),
    ),
    # Linear: the comparison rows; Circular: the conservation rings.
    **dict.fromkeys(
        ("evalue", "bitscore", "alignment_length"),
        frozenset(
            {"comparison_rows", "comparison_content", "conservation_sources", "conservation_tracks",
             "conservation_content"}
        ),
    ),
    # The reference only normalizes the loaded sources again.
    "conservation_reference": frozenset({"conservation_tracks", "conservation_content"}),
    # A ring's slot carries its label, and the labels and the record content
    # take the radial layout that holds the slots (measured on C3).
    "conservation_labels": frozenset(
        {"conservation_sources", "conservation_tracks", "conservation_content", "circular_labels",
         "circular_record_content"}
    ),
    **dict.fromkeys(
        (
            "gc_content_mode", "gc_content_min_percent", "gc_content_max_percent",
            "gc_content_large_tick_interval", "gc_content_small_tick_interval", "gc_content_tick_font_size",
            "show_gc_content_axis", "show_gc_content_ticks",
        ),
        frozenset({"gc_content_content"}),
    ),
    **dict.fromkeys(
        (
            "scale_interval", "scale_style", "scale_stroke_color", "scale_stroke_width", "scale_font_size",
            "ruler_label_font_size", "ruler_label_color", "hide_scale",
        ),
        # Linear: the length bar; Circular: the tick ring.
        frozenset({"length_bar_content", "ticks_content"}),
    ),
    **dict.fromkeys(
        (
            "keep_definition_left_aligned", "definition_line_style", "record_label",
            "record_subtitle", "show_replicon", "hide_accession", "hide_length",
        ),
        frozenset({"definition_content"}),
    ),
    # Linear: the definition; Circular: the record definition and the plot
    # title, which draws with the definition line interval that follows the
    # definition font size.
    "definition_font_size": frozenset(
        {"definition_content", "circular_definition_content", "circular_plot_title_content"}
    ),
    # Circular label placement reads the features section, which holds the
    # line strokes (section granularity); the placed labels do not change.
    **dict.fromkeys(("line_stroke_color", "line_stroke_width"), frozenset({"circular_labels"})),
    # The Circular definition draws the species and strain lines, in full or
    # short beside a plot title.
    **dict.fromkeys(
        ("species (no plot title)", "strain (no plot title)", "keep_full_definition_with_plot_title"),
        frozenset({"circular_definition_content"}),
    ),
}


# The Circular record content draws the embedded labels of the placed labels: it rebuilds with them.
_PLACED_LABELS = frozenset({"circular_labels", "circular_radial_labels"})
# The ribbons take the feature DOM index and both whole record placements,
# whose comparison extents these track contents set; they rebuild with them.
_RIBBON_INPUTS = frozenset({"feature_dom_index", "gc_content_content", "depth_content", "annotation_content"})
# The stages whose results the catalog item reads: it rebuilds when one of them does
# (Circular: the rings and the merged Multi-Record records too; not the axis circle).
_CATALOG_INPUTS = frozenset(
    {
        "record_content", "circular_record_content", "annotation_content", "comparison_content",
        "interactive_context", "conservation_content", "record_merge",
    }
)
# The catalog contents: each rebuilt one builds its ``catalog_snapshot`` (the
# reads the item is keyed by); an ``interactive_context`` rebuild alone builds none.
_CATALOG_CONTENTS = _CATALOG_INPUTS - {"interactive_context"}
# Options that change no catalog read (feature and match ids and fills): the
# catalog item must be a hit even when a content it reads rebuilds, so the page
# keeps its catalog. Options a case does not have are skipped.
_CATALOG_KEPT_OPTIONS = frozenset(
    {
        # Strokes.
        "axis_stroke_color", "axis_stroke_width", "line_stroke_color", "line_stroke_width",
        "scale_stroke_color", "scale_stroke_width",
        # Fonts.
        "definition_font_size", "scale_font_size", "ruler_label_font_size",
        # Scale, comparison height, track heights and radii.
        "scale_interval", "comparison_height", "feature_height", "track_axis_gap",
        # GC toggles.
        "gc_content_mode", "show_gc_content_axis", "show_gc_content_ticks",
    }
)


# The composition plan is keyed by the bounds of the primary content, the
# legend and the title, and by their placements. These stages hold those bounds.
_BOUNDS_STAGES = (
    _LEGEND_LAYOUT
    | _TITLE
    | _FEATURE_PLACEMENT
    | _DEPTH_DATA
    | {
        "record_content", "circular_record_content", "gc_content_content", "gc_skew_content",
        "annotation_content", "definition_content", "circular_definition_content",
        "length_bar_content", "comparison_content",
        # Circular: the placed labels are obstacles of the composition.
        "circular_labels",
    }
)
# The options the composition request reads without a stage: the placements. (The
# Circular definition text it reads is the `circular_definition_content` bounds.)
_COMPOSITION_INPUTS = frozenset({"legend", "plot_title_position"})
# Options that reshape a mark inside the box its content measures. With a paint
# option (its ID ends in ``_color``) they move no bound, so the plan must be a hit.
_MARK_SHAPE_OPTIONS = frozenset({"arrow_head_length_ratio", "arrow_shaft_width_ratio", "pairwise_match_style"})


def _moves_bounds(option_id: str, allowed: frozenset[str]) -> bool:
    if option_id.endswith("_color"):
        return False
    # A mark shape moves the Circular labels anchored on the marks.
    if option_id in _MARK_SHAPE_OPTIONS and "circular_labels" not in allowed:
        return False
    return option_id in _COMPOSITION_INPUTS or bool(allowed & _BOUNDS_STAGES)


def allowed_builds(option_id: str, mode: Mode) -> frozenset[str]:
    record_content, record_content_options = _RECORD_CONTENT[mode]
    allowed = (
        EXPECTED_BUILDS.get(option_id, frozenset())
        | CONTENT_BUILDS.get(option_id, frozenset())
        | (frozenset({record_content}) if option_id in record_content_options else frozenset())
        | (_CIRCULAR_LABELS if option_id in _CIRCULAR_LABEL_OPTIONS else frozenset())
        | (_RADIAL_LAYOUT_STAGES if option_id in _RADIAL_LAYOUT_OPTIONS else frozenset())
    )
    if allowed & _RIBBON_INPUTS:
        allowed |= {"comparison_content"}
    if allowed & _PLACED_LABELS:
        allowed |= {"circular_record_content"}
    if allowed & _RECORD_MERGE_INPUTS:
        allowed |= {"record_merge"}
    if _moves_bounds(option_id, allowed):
        allowed |= {"composition_plan"}
    if allowed & _CATALOG_CONTENTS:
        allowed |= {"catalog_snapshot"}
    if allowed & _CATALOG_INPUTS and option_id not in _CATALOG_KEPT_OPTIONS:
        allowed |= {"catalog_item"}
    return allowed


def p1_p2_variants(case_id: str, *, include_slow: bool) -> list[tuple[OptionSpec, Variant]]:
    """The variants of the P1 and P2 options of a case, in universe order."""

    variants = matrix_variants(CASES[case_id], include_slow=include_slow)
    return [(option, variant) for option, variant in variants if option.option_id in P1_P2_OPTIONS[case_id]]


def lost_hits(
    case_id: str,
    tables: Path,
    output: Path,
    *,
    include_slow: bool,
    part: slice = slice(None),
) -> list[str]:
    """Render the base, then each P1/P2 variant in ``part`` and the base again.

    An option with prerequisites (``OptionSpec.pre``) goes to its variant from
    its own base, the case base with those prerequisites, so the step changes
    the option alone.

    Name each stage built beyond the allowed set, each uncacheable stage call (it
    can never be a hit), and each failure of the error rule (``error_rule_failures``).
    """

    case = CASES[case_id]
    placeholders = Placeholders(tables)
    harness = EquivalenceHarness(case.mode, output)
    base = placeholders.expand(case, case.base)
    base_label = f"{case_id} base"
    failures: list[str] = []
    outcomes = []
    first = harness.warm(base)
    for option, variant in p1_p2_variants(case_id, include_slow=include_slow)[part]:
        option_base = base_arguments(case, option, placeholders)
        try:
            if option_base != base:
                harness.warm(option_base)
            outcome = harness.warm(state_arguments(case, {option: variant}, placeholders))
        except SystemExit:
            continue  # the CLI refuses this variant on this case
        label = f"{case_id} {option.option_id} {variant.label}"
        outcomes.append((label, outcome))
        if outcome.error is None:
            extra = set(outcome.stages.builds) - allowed_builds(option.option_id, case.mode)
            if extra:
                failures.append(f"{label}: built {sorted(extra)}")
            uncacheable = sorted(name for name, count in outcome.stages.uncacheable.items() if count)
            if uncacheable:
                failures.append(f"{label}: uncacheable {uncacheable}")
        returned = harness.warm(base)
        if returned.error is not None:
            failures.append(f"{label}, back to the base: render failed with {returned.error!r}")
    return failures + error_rule_failures((base_label, first), outcomes)


# The slow cases in chunks within the time budget (``cost_chunks``).
_SLOW_CHUNKS = [
    chunk
    for case_id in ("C2", "C3", "L2", "L3")
    for chunk in cost_chunks(case_id, f"{case_id}-", p1_p2_variants(case_id, include_slow=True))
]


@pytest.fixture(scope="module")
def tables(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return write_tables(tmp_path_factory.mktemp("expected-builds-tables"))


@pytest.mark.parametrize(
    ("case_id", "part"),
    [
        pytest.param("C1", slice(None), id="C1"),
        pytest.param("L1", slice(None), id="L1"),
        *(
            pytest.param(chunk.case_id, slice(chunk.start, chunk.stop), marks=UNIVERSE_MARKS, id=chunk.chunk_id)
            for chunk in _SLOW_CHUNKS
        ),
    ],
)
def test_a_p1_or_p2_option_builds_only_the_stages_of_its_content(
    case_id: str, part: slice, tables: Path, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    slow = request.node.get_closest_marker("slow") is not None
    failures = lost_hits(case_id, tables, tmp_path / case_id, include_slow=slow, part=part)
    assert not failures, "\n".join(failures)


def test_a_gc_toggle_rebuilds_neither_the_circular_definition_nor_the_title(tables: Path, tmp_path: Path) -> None:
    """The definition and the title take the canvas geometry, not the canvas switches."""

    case = CASES["C1"]
    placeholders = Placeholders(tables)
    harness = EquivalenceHarness(case.mode, tmp_path / "C1")
    assert harness.warm(placeholders.expand(case, case.base)).error is None
    [(option, variant)] = matrix_variants(case, option_ids=["gc"], first_variant_only=True)
    outcome = harness.warm(state_arguments(case, {option: variant}, placeholders))
    assert outcome.error is None
    assert not {"circular_definition_content", "circular_plot_title_content"} & set(outcome.stages.builds)


def test_a_plot_title_font_size_change_rebuilds_the_title_not_the_record_definition(
    tables: Path, tmp_path: Path
) -> None:
    """OV-1115: only the shared-common plot title draws in the plot title font size."""

    case = CASES["C1"]
    placeholders = Placeholders(tables)
    harness = EquivalenceHarness(case.mode, tmp_path / "C1")
    assert harness.warm(placeholders.expand(case, case.base)).error is None
    [(option, variant)] = matrix_variants(case, option_ids=["plot_title_font_size"], first_variant_only=True)
    outcome = harness.warm(state_arguments(case, {option: variant}, placeholders))
    assert outcome.error is None
    assert "circular_plot_title_content" in outcome.stages.builds
    assert "circular_definition_content" not in outcome.stages.builds


@pytest.mark.parametrize("case_id", ["C3", "L2"])
def test_a_depth_label_change_rebuilds_no_depth_path(case_id: str, tables: Path, tmp_path: Path) -> None:
    """OV-1229: the label is a legend input; the Depth config and path do not read it."""

    case = CASES[case_id]
    placeholders = Placeholders(tables)
    harness = EquivalenceHarness(case.mode, tmp_path / case_id)
    assert harness.warm(placeholders.expand(case, case.base)).error is None
    [(option, variant)] = matrix_variants(case, option_ids=["depth_track_label"], first_variant_only=True)
    outcome = harness.warm(state_arguments(case, {option: variant}, placeholders))
    assert outcome.error is None
    assert "depth_content" not in outcome.stages.builds
    assert set(outcome.stages.builds) <= allowed_builds("depth_track_label", case.mode)
