"""Structure ratchet for the Python package ``gbdraw``.

The rule is "Python structure" in
``docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`` (PY-STRUCTURE in
``tools/web-design-rule-guards.json``). The test reads ``gbdraw/**/*.py``,
except ``gbdraw/web/``, without importing it.

- Import layers. ``LAYERS`` puts each unit in a layer, from the bottom up. A
  unit is a subpackage of ``gbdraw``, a module directly under it, or
  ``__init__`` for the package itself. A module imports from its own layer and
  lower ones; ``UPWARD_IMPORTS`` lists the imports that go up today, marked
  `` (function)`` when they are only inside functions, so moving one to module
  level is a new entry.
- Import cycles. Module-level imports form no cycle. An import inside a
  function can close one; ``IMPORT_CYCLE_SIZES`` gives, for each module in a
  cycle today, the number of modules in its cycle.
- Function size. A function has at most ``MAX_FUNCTION_LINES`` lines and
  ``MAX_FUNCTION_PARAMETERS`` parameters. ``FUNCTION_LINES`` and
  ``FUNCTION_PARAMETERS`` list today's larger functions with their sizes.

The four baselines and the two limits are registered design-rule allowlists:
adding or raising an entry, or raising a limit, is an authority-only change
that needs the Owner's approval. Lowering or removing an entry is optional;
the test names the entries that may be lowered.

``LAYERS`` is not registered. The Gate keeps this file, ``LAYERS`` included,
from changing beside production code, but not alone: a pull request that moves
a unit to another layer loosens the rule for the units it passes, so it waits
for the Owner's approval and does not auto-merge. A new unit goes into
``LAYERS`` in a pull request without production code, before the code that
adds it.

Definitions. An import of ``gbdraw.a.b`` is an edge to that module. ``from X
import n`` is an edge to ``X.n`` when that is a module, else to ``X``;
relative imports are resolved, and imports under ``if TYPE_CHECKING:`` are not
counted. An import inside a ``def`` is function-level; any other is
module-level. A function is every ``def`` and ``async def``, methods and
nested functions included; it is named ``path::qualified.name``. Its lines run
from ``def`` to its last line, and its parameters are every name in its
signature, ``self``, ``cls``, ``*args`` and ``**kwargs`` included.
"""

from __future__ import annotations

import ast
import warnings
from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_ROOT = PROJECT_ROOT / "gbdraw"

LAYERS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("foundation", ("contract_shapes", "core", "data", "exceptions", "version")),
    (
        "model and geometry",
        (
            "analysis",
            "annotations",
            "auto_sizes",
            "canvas",
            "comparisons",
            "config",
            "definition_line_styles",
            "features",
            "io",
            "layout",
            "linear_comparison",
            "losat_setup",
            "mode_profiles",
            "svg",
            "tracks",
        ),
    ),
    ("placement", ("configurators", "labels", "legend")),
    ("drawing", ("diagrams", "render", "render_stages")),
    (
        "API and Session",
        (
            "api",
            "interface",
            "session",
            "session_container",
            "session_contract",
            "session_drawings",
            "session_io",
            "session_legacy_values",
            "session_migration",
            "session_request_codec",
            "session_resources",
        ),
    ),
    (
        "entry points and adapters",
        (
            "__init__",
            "_build_support",
            "_web_assets",
            "circular",
            "cli",
            "cli_utils",
            "crop_genbank",
            "linear",
            "web_support",
        ),
    ),
)

MAX_FUNCTION_LINES = 150
MAX_FUNCTION_PARAMETERS = 10

UPWARD_IMPORTS: frozenset[str] = frozenset(
    {
        "gbdraw.analysis.conservation -> gbdraw.api.prepared (function)",
        "gbdraw.analysis.depth_tracks -> gbdraw.api.options (function)",
        "gbdraw.analysis.depth_tracks -> gbdraw.api.prepared (function)",
        "gbdraw.analysis.depth_tracks -> gbdraw.configurators",
        "gbdraw.analysis.depth_tracks -> gbdraw.legend.sources",
        "gbdraw.analysis.depth_tracks -> gbdraw.render_stages",
        "gbdraw.analysis.protein_artifacts -> gbdraw.api.options (function)",
        "gbdraw.analysis.protein_artifacts -> gbdraw.session_request_codec (function)",
        "gbdraw.annotations.legend -> gbdraw.legend.sources",
        "gbdraw.annotations.legend -> gbdraw.legend.table",
        "gbdraw.annotations.legend -> gbdraw.render.patterns",
        "gbdraw.api -> gbdraw.web_support.capabilities",
        "gbdraw.comparisons.circular_losat -> gbdraw.api.options",
        "gbdraw.comparisons.circular_losat -> gbdraw.api.requests",
        "gbdraw.comparisons.linear_losat -> gbdraw.api.options",
        "gbdraw.config.models.labels -> gbdraw.labels.circular_types",
        "gbdraw.config.models.labels -> gbdraw.labels.label_style",
        "gbdraw.config.models.labels -> gbdraw.labels.policy",
        "gbdraw.config.models.render_profiles -> gbdraw.labels.circular_types",
        "gbdraw.config.models.render_profiles -> gbdraw.labels.filtering",
        "gbdraw.config.models.render_profiles -> gbdraw.labels.label_style",
        "gbdraw.config.models.render_profiles -> gbdraw.labels.policy",
        "gbdraw.configurators.legend -> gbdraw.render_stages",
        "gbdraw.core.sequence -> gbdraw.features.ids",
        "gbdraw.core.sequence -> gbdraw.features.overrides",
        "gbdraw.core.sequence -> gbdraw.features.visibility",
        "gbdraw.core.text -> gbdraw.svg.elements",
        "gbdraw.exceptions -> gbdraw.labels.fit (function)",
        "gbdraw.features.factory -> gbdraw.api.prepared (function)",
        "gbdraw.features.factory -> gbdraw.labels.filtering",
        "gbdraw.io.colors -> gbdraw.legend.table (function)",
        "gbdraw.io.comparisons -> gbdraw.api.prepared (function)",
        "gbdraw.io.regions -> gbdraw.crop_genbank",
        "gbdraw.labels.circular -> gbdraw.diagrams.circular.presets",
        "gbdraw.labels.filtering -> gbdraw.render_stages",
        "gbdraw.layout.linear -> gbdraw.labels.label_style (function)",
        "gbdraw.layout.linear -> gbdraw.labels.linear (function)",
        "gbdraw.legend.row_facts -> gbdraw.render_stages",
        "gbdraw.render.interactive_context -> gbdraw.api.prepared",
        "gbdraw.render.interactive_context -> gbdraw.web_support.feature_metadata",
        "gbdraw.render.interactive_context -> gbdraw.web_support.orthogroup_metadata",
        "gbdraw.render.interactive_svg -> gbdraw.web_support.feature_catalog (function)",
        "gbdraw.render_stages -> gbdraw.api.prepared (function)",
        "gbdraw.session -> gbdraw.web_support.feature_catalog (function)",
        "gbdraw.session_contract -> gbdraw.web_support.mode_scoped_settings",
        "gbdraw.session_drawings -> gbdraw.web_support.mode_scoped_settings",
        "gbdraw.session_io -> gbdraw.web_support.config_overrides (function)",
        "gbdraw.session_io -> gbdraw.web_support.feature_catalog (function)",
        "gbdraw.session_io -> gbdraw.web_support.mode_scoped_settings",
        "gbdraw.session_migration -> gbdraw.web_support.feature_catalog (function)",
        "gbdraw.session_migration -> gbdraw.web_support.feature_metadata (function)",
    }
)

IMPORT_CYCLE_SIZES: dict[str, int] = {
    "gbdraw.analysis.collinearity": 16,
    "gbdraw.analysis.collinearity_units": 16,
    "gbdraw.analysis.protein_artifacts": 16,
    "gbdraw.analysis.protein_colinearity": 16,
    "gbdraw.annotations.models": 16,
    "gbdraw.api.prepared": 16,
    "gbdraw.api.render": 16,
    "gbdraw.api.request_render": 16,
    "gbdraw.api.session_compat": 16,
    "gbdraw.auto_sizes": 13,
    "gbdraw.canvas.circular": 13,
    "gbdraw.canvas.linear": 13,
    "gbdraw.circular": 3,
    "gbdraw.cli_utils.session": 3,
    "gbdraw.config.models": 13,
    "gbdraw.config.models.canvas": 5,
    "gbdraw.config.models.render_profiles": 13,
    "gbdraw.configurators": 13,
    "gbdraw.configurators.blast": 13,
    "gbdraw.configurators.depth": 13,
    "gbdraw.configurators.features": 13,
    "gbdraw.configurators.gc": 13,
    "gbdraw.configurators.legend": 13,
    "gbdraw.core.figure_length": 5,
    "gbdraw.core.sequence": 16,
    "gbdraw.exceptions": 5,
    "gbdraw.features.colors": 16,
    "gbdraw.features.overrides": 16,
    "gbdraw.features.source": 16,
    "gbdraw.features.visibility": 16,
    "gbdraw.io.colors": 16,
    "gbdraw.io.comparisons": 16,
    "gbdraw.io.table_text": 16,
    "gbdraw.labels.circular": 2,
    "gbdraw.labels.circular_repair": 2,
    "gbdraw.labels.filtering": 13,
    "gbdraw.labels.fit": 5,
    "gbdraw.labels.linear": 2,
    "gbdraw.layout.figure_size": 5,
    "gbdraw.layout.linear": 2,
    "gbdraw.legend.sources": 16,
    "gbdraw.legend.table": 16,
    "gbdraw.linear": 3,
    "gbdraw.render.export": 16,
    "gbdraw.render.interactive_context": 16,
    "gbdraw.render.interactive_svg": 16,
    "gbdraw.render_stages": 13,
    "gbdraw.session": 16,
    "gbdraw.session_contract": 16,
    "gbdraw.session_drawings": 16,
    "gbdraw.session_io": 16,
    "gbdraw.session_migration": 16,
    "gbdraw.session_request_codec": 16,
    "gbdraw.session_resources": 16,
    "gbdraw.tracks.parsing": 16,
    "gbdraw.web_support.config_overrides": 16,
    "gbdraw.web_support.feature_catalog": 16,
}

FUNCTION_LINES: dict[str, int] = {
    "gbdraw/analysis/protein_colinearity.py::_build_adjacent_display_edges_by_pair": 179,
    "gbdraw/analysis/protein_colinearity.py::_build_anchor_core_orthogroups": 303,
    "gbdraw/analysis/protein_colinearity.py::extract_cds_proteins": 188,
    "gbdraw/analysis/protein_colinearity.py::extract_protein_identity_manifest": 239,
    "gbdraw/api/diagram.py::assemble_circular_diagram_from_record": 381,
    "gbdraw/api/diagram.py::assemble_circular_diagram_from_records": 588,
    "gbdraw/api/diagram.py::assemble_linear_diagram_from_records": 483,
    "gbdraw/api/options.py::_ModeDiagramOptions.__post_init__": 168,
    "gbdraw/circular.py::_get_args": 604,
    "gbdraw/circular.py::run_circular_from_namespace": 567,
    "gbdraw/diagrams/circular/assemble.py::_draw_resolved_circular_slot": 297,
    "gbdraw/diagrams/circular/assemble.py::_try_shift_labels_away_from_legend": 155,
    "gbdraw/diagrams/circular/assemble.py::add_record_on_circular_canvas": 635,
    "gbdraw/diagrams/circular/assemble.py::circular_labels": 228,
    "gbdraw/diagrams/circular/radial_layout.py::_resolve_circular_radial_layout": 405,
    "gbdraw/diagrams/linear/assemble.py::_build_linear_record_vertical_plans": 142,
    "gbdraw/diagrams/linear/assemble.py::_collect_linear_primary_bounds": 160,
    "gbdraw/diagrams/linear/builders.py::add_record_definition_group": 158,
    "gbdraw/features/factory.py::_build_feature_layers": 155,
    "gbdraw/labels/circular.py::_legacy_improved_label_placement_fc": 183,
    "gbdraw/labels/circular.py::_polish_total_leader_distance": 240,
    "gbdraw/labels/circular.py::_rebalance_inner_labels_strict_order": 374,
    "gbdraw/labels/circular.py::_refine_labels_to_preferred_hemisphere": 286,
    "gbdraw/labels/circular.py::_resolve_label_leader_line_collisions": 355,
    "gbdraw/labels/circular.py::_resolve_outer_label_overlaps_with_fixed_radii": 182,
    "gbdraw/labels/circular.py::improved_label_placement_fc": 720,
    "gbdraw/labels/circular.py::improved_label_placement_fc.optimize_with_shift_limit": 158,
    "gbdraw/labels/circular.py::prepare_label_list": 596,
    "gbdraw/labels/circular.py::prepare_label_list.place_external_labels": 215,
    "gbdraw/labels/circular.py::rearrange_labels_fc": 183,
    "gbdraw/labels/linear.py::linear_record_labels": 288,
    "gbdraw/layout/similarity_alignment.py::resolve_similarity_alignment": 236,
    "gbdraw/legend/circular_layout.py::build_circular_legend_layout": 174,
    "gbdraw/linear.py::_get_args": 869,
    "gbdraw/linear.py::run_linear_from_namespace": 538,
    "gbdraw/render/drawers/circular/annotations.py::circular_annotation_track_content": 147,
    "gbdraw/render/drawers/linear/annotations.py::linear_annotation_track_content": 167,
    "gbdraw/render/interactive_svg.py::enrich_svg": 194,
    "gbdraw/session_contract.py::describe_session_admission": 190,
    "gbdraw/session_io.py::_append_linear_gui_args": 157,
    "gbdraw/session_io.py::_attach_current_web_file_bindings": 172,
    "gbdraw/session_io.py::_migrate_legacy_linear_comparison_draft": 320,
    "gbdraw/session_io.py::_validate_current_comparison_authority": 178,
    "gbdraw/session_io.py::build_session_json": 181,
    "gbdraw/session_migration.py::project_session_adjunct_for_current_write": 182,
    "gbdraw/session_request_codec.py::_decode_canonical_request": 157,
    "gbdraw/session_request_codec.py::_decode_comparisons": 167,
    "gbdraw/session_request_codec.py::describe_canonical_request": 434,
    "gbdraw/session_request_codec.py::migrate_flat_config_overrides": 225,
    "gbdraw/tracks/circular.py::_normalize_circular_track_slots": 152,
    "gbdraw/web_support/feature_catalog.py::_normalized_biological_features": 203,
    "gbdraw/web_support/feature_catalog.py::_normalized_matches": 181,
    "gbdraw/web_support/feature_catalog.py::select_feature_catalog_item": 394,
    "gbdraw/web_support/feature_metadata.py::extract_features_from_records_payload": 211,
    "gbdraw/web_support/request_render.py::_render_canonical_web_request": 190,
}

FUNCTION_PARAMETERS: dict[str, int] = {
    "gbdraw/analysis/collinearity.py::build_collinearity_blocks_from_hits": 12,
    "gbdraw/analysis/collinearity.py::build_orthogroup_collinearity_blocks": 23,
    "gbdraw/analysis/collinearity.py::build_orthogroup_collinearity_blocks_from_hits": 13,
    "gbdraw/analysis/depth_tracks.py::normalize_depth_tracks": 14,
    "gbdraw/analysis/protein_artifacts.py::normalize_derived_losatp_settings": 20,
    "gbdraw/analysis/protein_colinearity.py::_append_record_local_orthogroups": 17,
    "gbdraw/analysis/protein_colinearity.py::_orthogroup_result_from_member_ids": 16,
    "gbdraw/analysis/protein_colinearity.py::_search_losatp_record_pairs": 11,
    "gbdraw/analysis/protein_colinearity.py::build_pairwise_protein_blastp_comparisons": 16,
    "gbdraw/analysis/protein_colinearity.py::build_rbh_orthogroup_protein_blastp_comparisons": 19,
    "gbdraw/api/diagram.py::_build_linear_diagram": 14,
    "gbdraw/api/diagram.py::_invoke_protein_analysis_helper": 26,
    "gbdraw/api/diagram.py::build_linear_diagram": 13,
    "gbdraw/api/diagram.py::build_linear_diagram_result": 13,
    "gbdraw/configurators/depth.py::DepthConfigurator.__init__": 14,
    "gbdraw/configurators/features.py::FeatureDrawingConfigurator.__init__": 12,
    "gbdraw/configurators/legend.py::LegendDrawingConfigurator._measurement": 13,
    "gbdraw/core/text.py::create_text_element": 11,
    "gbdraw/diagrams/circular/assemble.py::_draw_resolved_circular_slot": 27,
    "gbdraw/diagrams/circular/assemble.py::_record_label_layout": 13,
    "gbdraw/diagrams/circular/assemble.py::_resolve_radial_layout_fitting_definition": 13,
    "gbdraw/diagrams/circular/assemble.py::_slot_paint_data": 11,
    "gbdraw/diagrams/circular/assemble.py::build_circular_legend_table": 16,
    "gbdraw/diagrams/circular/assemble.py::circular_labels": 25,
    "gbdraw/diagrams/circular/assemble.py::circular_legend_table": 11,
    "gbdraw/diagrams/circular/builders.py::add_conservation_group_on_canvas": 12,
    "gbdraw/diagrams/circular/builders.py::add_depth_group_on_canvas": 14,
    "gbdraw/diagrams/circular/builders.py::add_gc_content_group_on_canvas": 13,
    "gbdraw/diagrams/circular/builders.py::add_gc_skew_group_on_canvas": 13,
    "gbdraw/diagrams/circular/builders.py::add_record_definition_group_on_canvas": 13,
    "gbdraw/diagrams/circular/builders.py::add_record_group_on_canvas": 14,
    "gbdraw/diagrams/circular/builders.py::add_tick_group_on_canvas": 12,
    "gbdraw/diagrams/circular/radial_layout.py::_inside_placement_window": 13,
    "gbdraw/diagrams/circular/radial_layout.py::_measure_anchored_inside": 11,
    "gbdraw/diagrams/circular/radial_layout.py::_place_inside_auto_fixed_width": 11,
    "gbdraw/diagrams/circular/radial_layout.py::_place_preferred_numeric_group": 11,
    "gbdraw/diagrams/circular/radial_layout.py::_resolve_circular_radial_layout": 12,
    "gbdraw/diagrams/circular/radial_layout.py::_try_measure_inside_interval": 12,
    "gbdraw/diagrams/linear/annotation_labels.py::place_linear_annotation_labels": 15,
    "gbdraw/diagrams/linear/assemble.py::_build_linear_record_vertical_plans": 16,
    "gbdraw/diagrams/linear/assemble.py::_collect_linear_primary_bounds": 19,
    "gbdraw/diagrams/linear/assemble.py::_linear_slot_footprints_for_record": 12,
    "gbdraw/diagrams/linear/assemble.py::_place_linear_annotation_labels_by_record": 17,
    "gbdraw/diagrams/linear/assemble.py::_record_collision_bands": 11,
    "gbdraw/diagrams/linear/assemble.py::build_linear_legend_table": 11,
    "gbdraw/diagrams/linear/builders.py::add_depth_group": 19,
    "gbdraw/diagrams/linear/builders.py::add_gc_content_group": 17,
    "gbdraw/diagrams/linear/builders.py::add_gc_skew_group": 17,
    "gbdraw/diagrams/linear/builders.py::add_record_definition_group": 19,
    "gbdraw/diagrams/linear/builders.py::add_record_group": 31,
    "gbdraw/diagrams/linear/precalc.py::_precalculate_label_dimensions": 11,
    "gbdraw/features/factory.py::_build_feature_layers": 15,
    "gbdraw/features/factory.py::create_feature_dict": 11,
    "gbdraw/features/factory.py::create_feature_layers": 15,
    "gbdraw/features/factory.py::create_feature_object": 12,
    "gbdraw/features/factory.py::create_gene_object": 12,
    "gbdraw/features/factory.py::create_repeat_object": 12,
    "gbdraw/features/objects.py::FeatureObject.__init__": 12,
    "gbdraw/features/objects.py::GeneObject.__init__": 15,
    "gbdraw/features/objects.py::RepeatObject.__init__": 14,
    "gbdraw/features/placement.py::resolve_record_feature_inputs": 11,
    "gbdraw/features/tracks.py::_find_best_track_split_overlaps_by_strand_indexed": 11,
    "gbdraw/labels/circular.py::_count_local_leader_line_collisions": 11,
    "gbdraw/labels/circular.py::_legacy_improved_label_placement_fc": 11,
    "gbdraw/labels/circular.py::calculate_angle_degrees": 11,
    "gbdraw/labels/circular.py::improved_label_placement_fc": 13,
    "gbdraw/labels/circular.py::prepare_label_list": 19,
    "gbdraw/labels/circular.py::rearrange_labels_fc": 13,
    "gbdraw/labels/label_ink.py::feature_label_ink": 13,
    "gbdraw/labels/linear.py::linear_record_labels": 15,
    "gbdraw/legend/table.py::_generated_legend_rows": 12,
    "gbdraw/legend/table.py::prepare_legend_table": 12,
    "gbdraw/render/composition.py::apply_composition_plan": 13,
    "gbdraw/render/drawers/circular/annotations.py::_label_at": 12,
    "gbdraw/render/drawers/circular/annotations.py::circular_annotation_labels": 11,
    "gbdraw/render/drawers/circular/definition.py::DefinitionDrawer.draw": 16,
    "gbdraw/render/drawers/circular/features.py::FeatureDrawer.draw": 12,
    "gbdraw/render/drawers/circular/features.py::FeatureDrawer.draw_path": 13,
    "gbdraw/render/drawers/circular/features.py::FeaturePathGenerator.__init__": 12,
    "gbdraw/render/drawers/circular/gc_skew.py::SkewDrawer.draw": 11,
    "gbdraw/render/drawers/linear/annotations.py::linear_annotation_track_content": 12,
    "gbdraw/render/drawers/linear/depth.py::DepthDrawer.draw": 11,
    "gbdraw/render/drawers/linear/features.py::FeatureDrawer.draw": 16,
    "gbdraw/render/drawers/linear/features.py::FeatureDrawer.draw_path": 16,
    "gbdraw/render/drawers/linear/features.py::FeaturePathGenerator.__init__": 12,
    "gbdraw/render/drawers/linear/gc_content.py::GcContentDrawer.draw": 12,
    "gbdraw/render/drawers/linear/gc_skew.py::SkewDrawer.draw": 12,
    "gbdraw/render/groups/circular/conservation.py::ConservationGroup.__init__": 11,
    "gbdraw/render/groups/circular/seq_record.py::SeqRecordGroup.__init__": 13,
    "gbdraw/render/groups/circular/ticks.py::TickGroup.__init__": 12,
    "gbdraw/render/groups/linear/definition.py::DefinitionGroup.__init__": 16,
    "gbdraw/render/groups/linear/depth.py::DepthGroup.__init__": 12,
    "gbdraw/render/groups/linear/gc_content.py::GcContentGroup.__init__": 11,
    "gbdraw/render/groups/linear/gc_skew.py::GcSkewGroup.__init__": 11,
    "gbdraw/render/groups/linear/seq_record.py::SeqRecordGroup.__init__": 24,
    "gbdraw/render/interactive_context.py::build_interactive_svg_context": 15,
    "gbdraw/render/interactive_context.py::interactive_svg_context": 12,
    "gbdraw/session_io.py::_migrate_legacy_linear_comparison_draft.add_edge": 11,
    "gbdraw/session_io.py::build_session_json": 12,
    "gbdraw/svg/circular_features.py::generate_circular_narrow_arrow_path": 11,
    "gbdraw/svg/circular_ticks.py::circular_tick_labels": 16,
    "gbdraw/svg/circular_ticks.py::circular_tick_marks": 11,
    "gbdraw/svg/circular_ticks.py::generate_circular_tick_labels": 21,
    "gbdraw/svg/circular_ticks.py::generate_circular_tick_paths": 12,
    "gbdraw/svg/circular_ticks.py::get_circular_tick_label_radius_bounds": 13,
    "gbdraw/svg/circular_ticks.py::resolve_circular_tick_label_geometry": 15,
    "gbdraw/svg/linear_features.py::create_arrow_path_linear": 13,
    "gbdraw/svg/linear_features.py::create_arrow_shaft_path_linear": 11,
    "gbdraw/svg/linear_features.py::create_arrowhead_path_linear": 12,
    "gbdraw/web_support/request_render.py::_render_canonical_web_request_with_prepared_inputs": 11,
    "gbdraw/web_support/request_render.py::_render_staged_canonical_web_request_with_prepared_inputs": 11,
}


@dataclass(frozen=True)
class Import:
    importer: str
    imported: str
    path: str
    line: int
    in_function: bool

    def __str__(self) -> str:
        where = "inside a function" if self.in_function else "at module level"
        return f"{self.path}:{self.line}: {self.importer} -> {self.imported} ({where})"


@dataclass(frozen=True)
class Function:
    key: str
    line: int
    lines: int
    parameters: int


@dataclass(frozen=True)
class Package:
    modules: frozenset[str]
    imports: tuple[Import, ...]
    functions: tuple[Function, ...]


def _module_name(path: Path) -> str:
    parts = path.relative_to(PROJECT_ROOT).with_suffix("").parts
    return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)


def unit_of(module: str) -> str:
    parts = module.split(".")
    return parts[1] if len(parts) > 1 else "__init__"


def _is_type_checking(test: ast.expr) -> bool:
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


class _Collector(ast.NodeVisitor):
    def __init__(self, module: str, path: Path, modules: frozenset[str]) -> None:
        self.module = module
        self.path = path.relative_to(PROJECT_ROOT).as_posix()
        self.package = module.split(".") if path.name == "__init__.py" else module.split(".")[:-1]
        self.modules = modules
        self.scope: list[str] = []
        self.function_depth = 0
        self.imports: list[Import] = []
        self.functions: list[Function] = []

    def _add(self, imported: str, node: ast.stmt) -> None:
        if imported in self.modules and imported != self.module:
            self.imports.append(
                Import(self.module, imported, self.path, node.lineno, self.function_depth > 0)
            )

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            self._add(alias.name, node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        named = node.module.split(".") if node.module else []
        if node.level:
            base = ".".join([*self.package[: len(self.package) - node.level + 1], *named])
        else:
            base = ".".join(named)
        for alias in node.names:
            submodule = f"{base}.{alias.name}"
            self._add(submodule if submodule in self.modules else base, node)

    def visit_If(self, node: ast.If) -> None:
        if not _is_type_checking(node.test):
            self.generic_visit(node)
            return
        for child in node.orelse:
            self.visit(child)

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    def visit_FunctionDef(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        arguments = node.args
        assert node.end_lineno is not None
        self.scope.append(node.name)
        self.functions.append(
            Function(
                key=f"{self.path}::{'.'.join(self.scope)}",
                line=node.lineno,
                lines=node.end_lineno - node.lineno + 1,
                parameters=len(arguments.posonlyargs)
                + len(arguments.args)
                + len(arguments.kwonlyargs)
                + (arguments.vararg is not None)
                + (arguments.kwarg is not None),
            )
        )
        self.function_depth += 1
        self.generic_visit(node)
        self.function_depth -= 1
        self.scope.pop()

    visit_AsyncFunctionDef = visit_FunctionDef


@cache
def read_package() -> Package:
    paths = {
        _module_name(path): path
        for path in sorted(PACKAGE_ROOT.rglob("*.py"))
        if not path.relative_to(PACKAGE_ROOT).as_posix().startswith("web/")
    }
    modules = frozenset(paths)
    imports: list[Import] = []
    functions: list[Function] = []
    for module, path in paths.items():
        collector = _Collector(module, path, modules)
        collector.visit(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        imports.extend(collector.imports)
        functions.extend(collector.functions)
    # `from x import (a, b)` gives one import per statement, not per name.
    return Package(modules, tuple(dict.fromkeys(imports)), tuple(functions))


def import_cycles(imports: list[Import]) -> dict[str, frozenset[str]]:
    """Each module on an import cycle, with the modules of its cycle."""
    graph: dict[str, set[str]] = defaultdict(set)
    for item in imports:
        graph[item.importer].add(item.imported)
    reach: dict[str, set[str]] = {}
    for start in graph:
        seen: set[str] = set()
        pending = list(graph[start])
        while pending:
            module = pending.pop()
            if module not in seen:
                seen.add(module)
                pending.extend(graph.get(module, ()))
        reach[start] = seen
    return {
        module: frozenset(other for other in found if module in reach.get(other, ()))
        for module, found in reach.items()
        if module in found
    }


def upward_imports() -> dict[str, list[Import]]:
    """Upward imports by ``importer -> imported``, with `` (function)`` when all are inside functions."""
    layer = {unit: index for index, (_, units) in enumerate(LAYERS) for unit in units}
    pairs: dict[tuple[str, str], list[Import]] = defaultdict(list)
    for item in read_package().imports:
        lower = layer.get(unit_of(item.importer))
        upper = layer.get(unit_of(item.imported))
        if lower is not None and upper is not None and lower < upper:
            pairs[item.importer, item.imported].append(item)
    return {
        f"{importer} -> {imported}" + (" (function)" if all(item.in_function for item in items) else ""): items
        for (importer, imported), items in pairs.items()
    }


def function_sizes(measure: str) -> dict[str, tuple[int, Function]]:
    """Size by function name; same-named definitions in a file count as their largest."""
    sizes: dict[str, tuple[int, Function]] = {}
    for function in read_package().functions:
        size = getattr(function, measure)
        if function.key not in sizes or size > sizes[function.key][0]:
            sizes[function.key] = (size, function)
    return sizes


def _warn_lowerable(name: str, entries: list[str]) -> None:
    if entries:
        warnings.warn(f"{name} may be lowered: {', '.join(entries)}")


def test_baselines_are_sorted() -> None:
    for name, baseline in (
        ("IMPORT_CYCLE_SIZES", IMPORT_CYCLE_SIZES),
        ("FUNCTION_LINES", FUNCTION_LINES),
        ("FUNCTION_PARAMETERS", FUNCTION_PARAMETERS),
    ):
        assert list(baseline) == sorted(baseline), f"keep {name} sorted by key"


def test_every_unit_is_in_one_layer() -> None:
    listed = [unit for _, units in LAYERS for unit in units]
    repeated = sorted({unit for unit in listed if listed.count(unit) > 1})
    assert not repeated, f"units in more than one layer: {repeated}"
    present = {unit_of(module) for module in read_package().modules}
    missing = sorted(present - set(listed))
    assert not missing, (
        f"units without a layer: {missing}. Put the module in an existing unit, or add "
        "the unit to LAYERS, in the layer of its role, in a pull request without "
        "production code before the code that adds it."
    )
    _warn_lowerable("LAYERS", sorted(set(listed) - present))


def test_no_new_upward_import() -> None:
    found = upward_imports()
    new = sorted(set(found) - UPWARD_IMPORTS)
    _warn_lowerable("UPWARD_IMPORTS", sorted(UPWARD_IMPORTS - set(found)))
    assert not new, (
        "imports from a lower layer to a higher one (LAYERS). Move the code to the "
        "layer of what it serves, or pass the value in from the caller. Adding an "
        "entry to UPWARD_IMPORTS needs the Owner's approval.\n"
        + "\n".join(str(item) for key in new for item in found[key])
    )


def test_module_level_imports_form_no_cycle() -> None:
    module_level = [item for item in read_package().imports if not item.in_function]
    cycles = sorted({tuple(sorted(members)) for members in import_cycles(module_level).values()})
    assert not cycles, "module-level imports form a cycle:\n" + "\n".join(
        str(item)
        for members in cycles
        for item in module_level
        if item.importer in members and item.imported in members
    )


def _joining_imports(members: frozenset[str], imports: list[Import]) -> list[Import]:
    """Imports in a cycle without which their two modules part.

    When the members have different entries, only imports between members with
    different entries: those joined what was apart.
    """
    mixed = len({IMPORT_CYCLE_SIZES.get(module, 0) for module in members}) > 1
    inside: dict[str, set[str]] = defaultdict(set)
    for item in imports:
        if item.importer in members and item.imported in members:
            inside[item.importer].add(item.imported)

    def reaches_without(start: str, target: str) -> bool:
        seen: set[str] = set()
        pending = [module for module in inside[start] if module != target]
        while pending:
            module = pending.pop()
            if module == target:
                return True
            if module not in seen:
                seen.add(module)
                pending.extend(inside[module])
        return False

    return sorted(
        (
            item
            for item in imports
            if item.importer in members
            and item.imported in members
            and not (
                mixed
                and IMPORT_CYCLE_SIZES.get(item.importer, 0) == IMPORT_CYCLE_SIZES.get(item.imported, 0)
            )
            and not reaches_without(item.importer, item.imported)
        ),
        key=str,
    )


def test_import_cycles_do_not_grow() -> None:
    imports = list(read_package().imports)
    cycles = import_cycles(imports)
    _warn_lowerable(
        "IMPORT_CYCLE_SIZES",
        [
            f"{module} {size} -> {len(cycles.get(module, ()))}"
            for module, size in IMPORT_CYCLE_SIZES.items()
            if len(cycles.get(module, ())) < size
        ],
    )
    grown = {
        members
        for module, members in cycles.items()
        if len(members) > IMPORT_CYCLE_SIZES.get(module, 0)
    }
    problems = [
        f"a cycle of {len(members)} modules; its members' entries are at most "
        f"{max(IMPORT_CYCLE_SIZES.get(module, 0) for module in members)}. Imports "
        "that join it:\n  "
        + "\n  ".join(str(item) for item in _joining_imports(members, imports))
        for members in sorted(grown, key=sorted)
    ]
    assert not problems, (
        "import cycles grew (IMPORT_CYCLE_SIZES). Remove the import that closes the "
        "cycle, usually one inside a function: move that function to the module it "
        "serves. Adding or raising an entry needs the Owner's approval.\n"
        + "\n".join(problems)
    )


def _assert_function_sizes(measure: str, limit: int, name: str, baseline: dict[str, int]) -> None:
    sizes = function_sizes(measure)
    _warn_lowerable(
        name,
        [
            f"{key} {size} -> {sizes[key][0] if key in sizes else 'removed'}"
            for key, size in baseline.items()
            if key not in sizes or sizes[key][0] < size
        ],
    )
    over = sorted(
        f"{function.key.split('::')[0]}:{function.line} {function.key.split('::')[1]}: "
        f"{size} {measure} (limit {limit}, {name} {baseline.get(key, 'none')})"
        for key, (size, function) in sizes.items()
        if size > max(limit, baseline.get(key, 0))
    )
    assert not over, (
        f"functions above {limit} {measure}, or above their {name} entry. Split them at "
        "an owner boundary, or take a typed parameter object. Adding or raising an "
        f"entry in {name} needs the Owner's approval.\n" + "\n".join(over)
    )


def test_functions_stay_within_the_line_limit() -> None:
    _assert_function_sizes("lines", MAX_FUNCTION_LINES, "FUNCTION_LINES", FUNCTION_LINES)


def test_functions_stay_within_the_parameter_limit() -> None:
    _assert_function_sizes(
        "parameters", MAX_FUNCTION_PARAMETERS, "FUNCTION_PARAMETERS", FUNCTION_PARAMETERS
    )
