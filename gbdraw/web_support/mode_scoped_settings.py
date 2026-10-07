"""The registry of mode-scoped Web settings (Session 46, PD-OI-086).

Session 46 keeps every diagram setting and editor edit per diagram mode in
``modes.circular`` and ``modes.linear``. Each slice holds the Web draft parts
named by the rows below and nothing else: a ``config`` with its ``form``,
``adv``, and ``losat`` groups, ``features``, the Legend part of
``editorState`` (with ``featureStrokes``), and a few ``ui`` keys.

A row is ``{domain, path, modes, key, migrate}``:

- ``domain``: the slice container that holds the value (``config.form``,
  ``config``, ``editorState.legend``, ...).
- ``path``: the key of the value inside ``domain``.
- ``modes``: ``both``, or the one mode a MODE-ONLY value belongs to.
- ``key``: for a keyed-row domain, the name of the row key (compared entry by
  entry by the 46 -> 47 step); otherwise ``None``.
- ``migrate``: the token by which the split of a Session 27-44 draft fills the
  two slices (``MIGRATION_TOKENS``).

The Web reads the same rows from ``gbdraw/web/js/mode-scoped-settings.generated.js``,
which ``tools/generate_mode_scoped_settings.py`` writes from this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

DiagramMode = Literal["circular", "linear"]
RowModes = Literal["both", "circular", "linear"]

DIAGRAM_MODES: tuple[DiagramMode, ...] = ("circular", "linear")
# The registry revision of Phase E's row list this table follows.
MODE_SCOPED_SETTINGS_REVISION = 3

# The split tokens (Phase E plan 4.2):
#   copy            the value goes to both slices.
#   own             the value goes to the slice of the row's mode only.
#   profile         the active mode's slice takes the flat value; the other
#                   slice takes the saved mode profile of that mode.
#   layout          ui.layoutPreferences: each mode's slot goes to its slice.
#   show-if-source  Show Depth stays on only in a mode with a Depth source.
#   depth           copied, then trimmed to each mode's Depth source width.
#   result-mode     the value goes to the slice of the saved Result's mode.
#   by-scope        each row goes to the slice its ``scope`` names.
#   by-side         each placement row goes to its side's mode; Main to both.
#   by-leaf         each override leaf goes to the modes that own its path.
#   by-binding      an annotation bound to one mode's record goes to that
#                   mode; the others go to both slices.
MIGRATION_TOKENS = frozenset(
    {
        "copy",
        "own",
        "profile",
        "layout",
        "show-if-source",
        "depth",
        "result-mode",
        "by-scope",
        "by-side",
        "by-leaf",
        "by-binding",
    }
)


@dataclass(frozen=True)
class ModeScopedSetting:
    """One registry row: a value that each mode's slice holds on its own."""

    domain: str
    path: str
    modes: RowModes
    key: str | None
    migrate: str

    def to_json(self) -> dict[str, Any]:
        return {
            "domain": self.domain,
            "path": self.path,
            "modes": self.modes,
            "key": self.key,
            "migrate": self.migrate,
        }


def _rows(domain: str, *entries: tuple[str, RowModes, str]) -> tuple[ModeScopedSetting, ...]:
    return tuple(ModeScopedSetting(domain, path, modes, None, migrate) for path, modes, migrate in entries)


_FORM_ROWS = _rows(
    "config.form",
    ("prefix", "both", "copy"),
    ("species", "circular", "own"),
    ("strain", "circular", "own"),
    ("plot_title", "both", "profile"),
    ("track_type", "circular", "own"),
    ("linear_track_layout", "linear", "own"),
    ("show_scale", "both", "copy"),
    ("scale_style", "linear", "own"),
    ("linear_ruler_on_axis", "linear", "own"),
    ("labels_mode", "circular", "own"),
    ("show_labels_linear", "linear", "own"),
    ("multi_record_canvas", "circular", "own"),
    ("circular_record_selector", "circular", "own"),
    ("circular_region_start", "circular", "own"),
    ("circular_region_end", "circular", "own"),
    ("circular_reverse", "circular", "own"),
    ("circular_record_label", "circular", "own"),
    ("circular_record_subtitle", "circular", "own"),
    ("separate_strands", "both", "copy"),
    ("suppress_gc", "circular", "own"),
    ("suppress_skew", "circular", "own"),
    ("align_center", "linear", "own"),
    ("keep_definition_left_aligned", "linear", "own"),
    ("show_gc", "linear", "own"),
    ("show_skew", "linear", "own"),
    ("show_depth", "both", "show-if-source"),
    ("normalize_length", "linear", "own"),
)

_ADV_ROWS = _rows(
    "config.adv",
    ("features", "both", "copy"),
    ("feature_shapes", "both", "copy"),
    ("arrow_head_length_ratio", "both", "copy"),
    ("arrow_shaft_width_ratio", "both", "copy"),
    ("window_size", "both", "copy"),
    ("step_size", "both", "copy"),
    ("nt", "both", "copy"),
    ("def_font_size", "both", "profile"),
    ("circular_definition_interval", "circular", "own"),
    ("label_font_size", "both", "copy"),
    ("circular_label_spacing", "circular", "own"),
    ("linear_label_spacing", "linear", "own"),
    ("label_rendering", "both", "copy"),
    ("circular_label_placement", "circular", "own"),
    ("label_placement", "linear", "own"),
    ("label_rotation", "linear", "own"),
    ("block_stroke_width", "both", "copy"),
    ("block_stroke_color", "both", "copy"),
    ("line_stroke_width", "both", "copy"),
    ("line_stroke_color", "both", "copy"),
    ("axis_stroke_width", "both", "copy"),
    ("axis_stroke_color", "both", "profile"),
    ("legend_box_size", "both", "copy"),
    ("legend_font_size", "both", "copy"),
    ("resolve_overlaps", "both", "copy"),
    ("feature_overlap_tolerance_bp", "both", "copy"),
    ("feature_height", "linear", "own"),
    ("track_axis_gap", "linear", "own"),
    ("linear_show_replicon", "linear", "own"),
    ("linear_accession_visibility", "linear", "own"),
    ("linear_length_visibility", "linear", "own"),
    ("linear_definition_line_styles", "linear", "own"),
    ("gc_height", "linear", "own"),
    ("depth_height", "linear", "own"),
    ("depth_color", "both", "copy"),
    ("depth_tracks", "both", "depth"),
    ("depth_window_size", "both", "copy"),
    ("depth_step_size", "both", "copy"),
    ("depth_share_axis", "both", "copy"),
    ("depth_min", "both", "copy"),
    ("depth_max", "both", "copy"),
    ("depth_normalize", "both", "copy"),
    ("depth_show_axis", "both", "copy"),
    ("depth_show_ticks", "both", "copy"),
    ("linear_track_slots_enabled", "linear", "own"),
    ("linear_track_slots_schema_version", "linear", "own"),
    ("linear_track_slots_axis_index", "linear", "own"),
    ("linear_track_slots", "linear", "own"),
    ("gc_content_mode", "both", "copy"),
    ("gc_content_min_percent", "both", "copy"),
    ("gc_content_max_percent", "both", "copy"),
    ("gc_content_show_axis", "both", "copy"),
    ("gc_content_show_ticks", "both", "copy"),
    ("gc_content_tick_interval", "both", "copy"),
    ("gc_content_small_tick_interval", "both", "copy"),
    ("gc_content_tick_font_size", "both", "copy"),
    ("comparison_height", "linear", "own"),
    ("pairwise_match_style", "linear", "profile"),
    ("scale_interval", "both", "copy"),
    ("scale_font_size", "linear", "own"),
    ("ruler_label_font_size", "linear", "own"),
    ("scale_stroke_width", "linear", "own"),
    ("scale_stroke_color", "linear", "own"),
    ("ruler_label_color", "linear", "own"),
    ("circular_grouping_intent", "circular", "own"),
    ("multi_record_size_mode", "circular", "own"),
    ("multi_record_min_radius_ratio", "circular", "own"),
    ("multi_record_column_gap_ratio", "circular", "own"),
    ("multi_record_row_gap_ratio", "circular", "own"),
    ("multi_record_positions", "circular", "own"),
    ("tick_label_font_size", "circular", "own"),
    ("plot_title_font_size", "both", "profile"),
    ("keep_full_definition_with_plot_title", "circular", "own"),
    ("center_reserved_radius", "circular", "own"),
    ("feature_width_circular", "circular", "own"),
    ("depth_width_circular", "circular", "own"),
    ("gc_content_width_circular", "circular", "own"),
    ("gc_content_radius_circular", "circular", "own"),
    ("gc_skew_width_circular", "circular", "own"),
    ("gc_skew_radius_circular", "circular", "own"),
    ("circular_track_slots_enabled", "circular", "own"),
    ("circular_track_slots_schema_version", "circular", "own"),
    ("circular_track_slots_axis_index", "circular", "own"),
    ("circular_track_slots", "circular", "own"),
    ("outer_label_x_offset", "circular", "own"),
    ("outer_label_y_offset", "circular", "own"),
    ("inner_label_x_offset", "circular", "own"),
    ("inner_label_y_offset", "circular", "own"),
    ("min_bitscore", "both", "profile"),
    ("evalue", "both", "profile"),
    ("identity", "both", "profile"),
    ("alignment_length", "both", "profile"),
    # Registry revision 3 classes the rich feature popup as app-level but
    # keeps its Session 46 path in the draft; revision 4 settles it.
    ("rich_feature_popup", "both", "copy"),
)

# ``losat`` is listed by leaf group; its execution keys are app-level
# (``ui.losatExecution``), so a slice holds none of them.
_LOSAT_ROWS = _rows(
    "config.losat",
    ("outfmt", "both", "copy"),
    ("blastn", "both", "copy"),
    ("blastp", "linear", "own"),
)

_CONFIG_ROWS = (
    *_rows(
        "config",
        ("colors", "both", "copy"),
        # Not in registry revision 3; real Sessions 40-44 hold it, and it says
        # whether ``colors`` replaces or overrides the palette.
        ("colorsAreOverrides", "both", "copy"),
        ("palette", "both", "copy"),
        ("rules", "both", "copy"),
        ("qualifierPriorityRules", "both", "copy"),
        ("filterMode", "both", "copy"),
        ("whitelist", "both", "copy"),
        ("blacklistText", "both", "copy"),
        ("losatProgram", "linear", "own"),
        ("circularConservation", "circular", "own"),
        ("unmanagedConfigOverrides", "both", "by-leaf"),
        ("linearRecordLayout", "linear", "own"),
        ("linearComparisonPlan", "linear", "own"),
        ("importedComparisonResolution", "linear", "own"),
        ("webEdits", "linear", "own"),
        # The committed request's CLI options travel with its mode's slice.
        ("cliOptions", "both", "result-mode"),
    ),
    ModeScopedSetting("config", "annotationSets", "both", "annotation-set", "by-binding"),
    ModeScopedSetting("config", "recordDisplayDrafts", "both", "record-display", "by-scope"),
    ModeScopedSetting("config", "featurePlacementOverrides", "both", "placement", "by-side"),
)

_FEATURES_ROWS = (
    ModeScopedSetting("features", "featureOverrides", "both", "feature", "result-mode"),
    *_rows(
        "features",
        ("featureColorOverrides", "both", "result-mode"),
        ("featureVisibilityManualRules", "both", "copy"),
        ("labelOverrideRows", "both", "copy"),
        ("labelTextBulkOverrides", "both", "copy"),
    ),
)

_EDITOR_ROWS = (
    *_rows(
        "editorState.legend",
        ("entries", "both", "result-mode"),
        ("deletedEntries", "both", "result-mode"),
        ("colorOverrides", "both", "result-mode"),
        ("strokeOverrides", "both", "result-mode"),
        ("addedCaptions", "both", "result-mode"),
    ),
    *_rows("editorState", ("featureStrokes", "both", "result-mode")),
)

_UI_ROWS = _rows(
    "ui",
    ("layoutPreferences", "both", "layout"),
    ("canvasPadding", "both", "copy"),
    ("pendingPaletteName", "both", "copy"),
    ("pendingPaletteColors", "both", "copy"),
    ("linearTypographyLinked", "linear", "own"),
)

MODE_SCOPED_SETTINGS: tuple[ModeScopedSetting, ...] = (
    *_FORM_ROWS,
    *_ADV_ROWS,
    *_LOSAT_ROWS,
    *_CONFIG_ROWS,
    *_FEATURES_ROWS,
    *_EDITOR_ROWS,
    *_UI_ROWS,
)

# ``config.losat`` keys that Session 46 keeps once, app-wide, in ``ui.losatExecution``.
LOSAT_EXECUTION_FIELDS = ("executionMode", "totalThreadBudget", "threadsPerJob", "parallelWorkers")


def _slice_containers() -> dict[str, frozenset[str]]:
    """Each slice container (``""`` is the slice itself) and the keys it may hold."""

    containers: dict[str, set[str]] = {"": set()}
    for row in MODE_SCOPED_SETTINGS:
        parts = row.domain.split(".")
        for depth in range(len(parts)):
            parent = ".".join(parts[:depth])
            containers.setdefault(parent, set()).add(parts[depth])
            containers.setdefault(".".join(parts[: depth + 1]), set())
        containers[row.domain].add(row.path)
    return {container: frozenset(keys) for container, keys in containers.items()}


# Each slice container and the keys it may hold; a slice holds nothing else.
SLICE_CONTAINERS: dict[str, frozenset[str]] = _slice_containers()


def unmanaged_config_override_modes(path: str) -> tuple[DiagramMode, ...]:
    """The modes that may hold a GUI-unmanaged config override leaf (``by-leaf``).

    A leaf under one mode's renderer prefixes belongs to that mode; any other
    leaf is shared by both. The prefixes are the typed options' own rule
    (``gbdraw.api.options``), which rejects the other mode's leaves.
    """

    from gbdraw.api.options import _MODE_CONFIG_OVERRIDE_PREFIXES

    for mode in DIAGRAM_MODES:
        if any(path == prefix or path.startswith(f"{prefix}.") for prefix in _MODE_CONFIG_OVERRIDE_PREFIXES[mode]):
            return (mode,)
    return DIAGRAM_MODES


def mode_scoped_settings_payload() -> dict[str, Any]:
    """The registry as the Web's generated leaf holds it."""

    from gbdraw.api.options import _MODE_CONFIG_OVERRIDE_PREFIXES

    return {
        "revision": MODE_SCOPED_SETTINGS_REVISION,
        "modes": list(DIAGRAM_MODES),
        "migrationTokens": sorted(MIGRATION_TOKENS),
        "rows": [row.to_json() for row in MODE_SCOPED_SETTINGS],
        "losatExecutionFields": list(LOSAT_EXECUTION_FIELDS),
        "unmanagedConfigOverrideModePrefixes": {
            mode: list(_MODE_CONFIG_OVERRIDE_PREFIXES[mode]) for mode in DIAGRAM_MODES
        },
    }


__all__ = [
    "DIAGRAM_MODES",
    "LOSAT_EXECUTION_FIELDS",
    "MIGRATION_TOKENS",
    "MODE_SCOPED_SETTINGS",
    "MODE_SCOPED_SETTINGS_REVISION",
    "ModeScopedSetting",
    "SLICE_CONTAINERS",
    "mode_scoped_settings_payload",
    "unmanaged_config_override_modes",
]
