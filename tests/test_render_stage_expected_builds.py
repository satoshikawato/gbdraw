"""Expected builds (STEP2 design 3.2, PL4).

For every option of class P1 or P2 (the option changes only its own content,
P2 also the legend, the title or the composition plan), a warm render that goes
from the case base to one variant may build only the render stages that hold
the content the option changes. A stage that builds beyond that set is a lost
hit, so a later change that drops one fails here.

The classes are those of the DEPMAP measurement (``P1_P2_OPTIONS``). Each
variant starts from a cache that holds the base: the harness renders the base
again after each variant.

Two gaps are listed apart, so removing them shrinks the lists:

- ``LEGEND_INPUT_GAP``: the legend table stages take whole configurators
  (feature, GC, skew, depth track data), so any field of them rebuilds the
  legend even when no legend row changes.
- ``RECORD_PRESENTATION_GAP``: a record label or subtitle is written into the
  record, so it is part of the record identity and rebuilds every stage of
  the record.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from tests.utils.render_stage_equivalence import (
    CASES,
    SLOW_VARIANTS,
    EquivalenceHarness,
    Placeholders,
    options_for,
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

_LEGEND = frozenset({"legend_table", "legend_measurement"})
_LEGEND_LAYOUT = frozenset({"legend_measurement"})
_FEATURE_PLACEMENT = frozenset({"feature_layers", "feature_dom_index"})
_LABEL_FILTER = _FEATURE_PLACEMENT | {"label_filter_maps"}
_FEATURE_EDIT = _FEATURE_PLACEMENT | _LEGEND | {"catalog_normalization", "interactive_context"}
_DEPTH_DATA = frozenset({"depth_track_data"})
_DEPTH_FRAME = _DEPTH_DATA | {"depth_frame"}

# The stages that hold the content each option changes. An option not named builds nothing.
EXPECTED_BUILDS: Mapping[str, frozenset[str]] = {
    **dict.fromkeys(
        (
            "legend", "legend_box_size", "legend_font_size", "definition_font_size", "definition_line_style",
            "hide_accession", "hide_scale", "keep_definition_left_aligned", "ruler_label_font_size",
            "scale_font_size", "scale_style", "show_replicon", "align_center (rows)", "axis_stroke_width",
            "multi_record_column_gap_ratio", "multi_record_row_gap_ratio",
        ),
        _LEGEND_LAYOUT,
    ),
    **dict.fromkeys(
        (
            "separate_strands", "resolve_overlaps", "feature_overlap_tolerance_bp", "show_labels",
            "track_layout", "track_layout (rows)", "track_axis_gap",
        ),
        _FEATURE_PLACEMENT,
    ),
    **dict.fromkeys(("label_table", "label_whitelist", "label_blacklist", "qualifier_priority"), _LABEL_FILTER),
    # A placement row is a feature edit: the stages that read the edited features rebuild.
    **dict.fromkeys(("feature_override_table", "feature_placement_table"), _FEATURE_EDIT),
    "feature_visibility_table": _FEATURE_EDIT | {"visibility_rules"},
    "annotation_table": _LEGEND | {"catalog_normalization", "interactive_context"},
    **dict.fromkeys(("depth_window", "depth_step", "depth_min", "depth_max", "depth_log_scale"), _DEPTH_FRAME),
    **dict.fromkeys(
        (
            "depth_track_large_tick_interval", "depth_track_small_tick_interval", "depth_track_tick_font_size",
            "show_depth_axis", "show_depth_ticks", "depth_large_tick_interval", "depth_small_tick_interval",
            "depth_tick_font_size", "share_depth_axis",
        ),
        _DEPTH_DATA,
    ),
    "depth_track_label": _DEPTH_DATA | _LEGEND,
    "conservation_labels": _LEGEND,
    # The pairwise identity legend is computed from the filtered comparisons.
    **dict.fromkeys(("evalue", "bitscore", "alignment_length"), _LEGEND),
}

LEGEND_INPUT_GAP: frozenset[str] = frozenset({
    "arrow_head_length_ratio", "arrow_shaft_width_ratio", "line_stroke_color", "line_stroke_width",
    "gc_content_mode", "gc_content_min_percent", "gc_content_max_percent", "gc_content_large_tick_interval",
    "gc_content_small_tick_interval", "gc_content_tick_font_size", "show_gc_content_axis",
    "show_gc_content_ticks", "feature_overlap_tolerance_bp", "pairwise_match_style", "label_blacklist",
    "track_layout", "track_layout (rows)", "track_axis_gap", "depth_window", "depth_step", "depth_min",
    "depth_max", "depth_log_scale", "depth_track_large_tick_interval", "depth_track_small_tick_interval",
    "depth_track_tick_font_size", "show_depth_axis", "show_depth_ticks", "depth_large_tick_interval",
    "depth_small_tick_interval", "depth_tick_font_size", "share_depth_axis",
})
RECORD_PRESENTATION_GAP: Mapping[str, frozenset[str]] = dict.fromkeys(
    ("record_label", "record_subtitle"),
    _FEATURE_EDIT | _DEPTH_FRAME | {"gc_skew_data"},
)


def allowed_builds(option_id: str) -> frozenset[str]:
    allowed = EXPECTED_BUILDS.get(option_id, frozenset())
    if option_id in LEGEND_INPUT_GAP:
        allowed |= _LEGEND
    return allowed | RECORD_PRESENTATION_GAP.get(option_id, frozenset())


def lost_hits(case_id: str, tables: Path, output: Path, *, include_slow: bool) -> list[str]:
    """Render the base, then each P1/P2 variant and the base again.

    Name each stage built beyond the allowed set, each uncacheable stage call (it
    can never be a hit), and each render that failed.
    """

    case = CASES[case_id]
    placeholders = Placeholders(tables)
    harness = EquivalenceHarness(case.mode, output)
    base = placeholders.expand(case, case.base)
    assert harness.warm(base).error is None
    failures: list[str] = []
    for option in options_for(case):
        if option.option_id not in P1_P2_OPTIONS[case_id]:
            continue
        for variant in option.variants:
            if not include_slow and (option.option_id, variant.label) in SLOW_VARIANTS:
                continue
            try:
                outcome = harness.warm(state_arguments(case, {option: variant}, placeholders))
            except SystemExit:
                continue  # the CLI refuses this variant on this case
            label = f"{case_id} {option.option_id} {variant.label}"
            if outcome.error is not None:
                failures.append(f"{label}: render failed with {outcome.error!r}")
            else:
                extra = set(outcome.stages.builds) - allowed_builds(option.option_id)
                if extra:
                    failures.append(f"{label}: built {sorted(extra)}")
                uncacheable = sorted(name for name, count in outcome.stages.uncacheable.items() if count)
                if uncacheable:
                    failures.append(f"{label}: uncacheable {uncacheable}")
            returned = harness.warm(base)
            if returned.error is not None:
                failures.append(f"{label}, back to the base: render failed with {returned.error!r}")
    return failures


@pytest.fixture(scope="module")
def tables(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return write_tables(tmp_path_factory.mktemp("expected-builds-tables"))


@pytest.mark.parametrize(
    "case_id",
    [
        "C1",
        "L1",
        pytest.param("C2", marks=pytest.mark.slow),
        pytest.param("C3", marks=pytest.mark.slow),
        pytest.param("L2", marks=pytest.mark.slow),
        pytest.param("L3", marks=pytest.mark.slow),
    ],
)
def test_a_p1_or_p2_option_builds_only_the_stages_of_its_content(
    case_id: str, tables: Path, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    slow = request.node.get_closest_marker("slow") is not None
    failures = lost_hits(case_id, tables, tmp_path / case_id, include_slow=slow)
    assert not failures, "\n".join(failures)
