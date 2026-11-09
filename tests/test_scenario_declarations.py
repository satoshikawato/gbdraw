"""Every command-line option belongs to a capability of the scenario manifest.

``docs/scenarios/manifest.json`` names the capabilities that the scenarios
exercise and gives each one public owner page. ``CLI_OPTION_CAPABILITIES``
places each long option of ``gbdraw circular``, ``gbdraw linear`` and
``gbdraw render`` (their ``--help``) in the capability whose owner page
documents it. A new option fails this test until it is placed there, or until
a new capability with its scenario and owner page is added to the manifest.

``UNDECLARED_CLI_OPTIONS`` lists the options that no capability covers today.
It is a registered design-rule allowlist (``tools/web-design-rule-guards.json``,
SCENARIO-DECLARATION, kind ``py-set``): adding an option is an authority-only
change that needs the Owner's approval, and declaring one (moving it to
``CLI_OPTION_CAPABILITIES``) is free.
"""

from __future__ import annotations

import json
from pathlib import Path

from tests.test_cli_reference import _live_help_options

MANIFEST_PATH = (
    Path(__file__).resolve().parents[1] / "docs" / "scenarios" / "manifest.json"
)
COMMANDS = ("circular", "linear", "render")

CLI_OPTION_CAPABILITIES: dict[str, str] = {
    "--gbk": "input.genbank-embl-ddbj",
    "--gff": "input.gff3-fasta",
    "--fasta": "input.gff3-fasta",
    "--record_id": "input.record-selection",
    "--region": "input.regions",
    "--reverse_complement": "input.reverse-complement",
    "--auto_orient": "input.reverse-complement",
    "--records_table": "input.tsv-manifests",
    "--comparisons_table": "input.tsv-manifests",
    "--conservation_table": "input.tsv-manifests",
    "--multi_record_canvas": "layout.circular-multi",
    "--multi_record_column_gap_ratio": "layout.circular-multi",
    "--multi_record_min_radius_ratio": "layout.circular-multi",
    "--multi_record_position": "layout.circular-multi",
    "--multi_record_row_gap_ratio": "layout.circular-multi",
    "--multi_record_size_mode": "layout.circular-multi",
    "--one_circle": "layout.circular-multi",
    "--one_circle_gap_degrees": "layout.circular-multi",
    "--one_circle_min_length": "layout.circular-multi",
    "--one_circle_no_record_band": "layout.circular-multi",
    "--align_center": "layout.linear-multi-row",
    "--comparison_height": "layout.linear-multi-row",
    "--linear_record_gap": "layout.linear-multi-row",
    "--normalize_length": "layout.linear-multi-row",
    "--center_reserved_radius": "layout.track-placement",
    "--feature_height": "layout.track-placement",
    "--feature_width": "layout.track-placement",
    "--no_track_narrowing": "layout.track-placement",
    "--separate_strands": "layout.track-placement",
    "--track_axis_gap": "layout.track-placement",
    "--track_layout": "layout.track-placement",
    "--track_narrowing_min_px": "layout.track-placement",
    "--track_type": "layout.track-placement",
    "--hide_scale": "layout.rulers",
    "--ruler_label_color": "layout.rulers",
    "--ruler_label_font_size": "layout.rulers",
    "--ruler_on_axis": "layout.rulers",
    "--scale_font_size": "layout.rulers",
    "--scale_interval": "layout.rulers",
    "--scale_stroke_color": "layout.rulers",
    "--scale_stroke_width": "layout.rulers",
    "--scale_style": "layout.rulers",
    "--tick_label_font_size": "layout.rulers",
    "--axis_stroke_color": "layout.axes",
    "--axis_stroke_width": "layout.axes",
    "--plot_title": "layout.titles",
    "--plot_title_font_size": "layout.titles",
    "--plot_title_position": "layout.titles",
    "--definition_font_size": "layout.definitions",
    "--definition_line_style": "layout.definitions",
    "--hide_accession": "layout.definitions",
    "--hide_definition": "layout.definitions",
    "--hide_length": "layout.definitions",
    "--keep_definition_left_aligned": "layout.definitions",
    "--keep_full_definition_with_plot_title": "layout.definitions",
    "--record_label": "layout.definitions",
    "--record_subtitle": "layout.definitions",
    "--show_replicon": "layout.definitions",
    "--species": "layout.definitions",
    "--strain": "layout.definitions",
    "--legend": "layout.legends",
    "--legend_box_size": "layout.legends",
    "--legend_font_size": "layout.legends",
    "--blast": "comparison.uploaded-tables",
    "--alignment_length": "comparison.uploaded-tables",
    "--bitscore": "comparison.uploaded-tables",
    "--evalue": "comparison.uploaded-tables",
    "--identity": "comparison.uploaded-tables",
    "--identity_scale": "comparison.uploaded-tables",
    "--pairwise_match_style": "comparison.uploaded-tables",
    "--losatp_mode": "comparison.losatp-groups",
    "--losatp_member_max_hits": "comparison.losatp-groups",
    "--losatp_max_hits": "comparison.losatp-pairwise",
    "--losatp_max_target_seqs": "comparison.losatp-pairwise",
    "--collinear_anchor_mode": "comparison.losatp-collinear",
    "--collinear_color_mode": "comparison.losatp-collinear",
    "--collinear_infer_orthogroups": "comparison.losatp-collinear",
    "--collinear_max_conflicts_in_merge_gap": "comparison.losatp-collinear",
    "--collinear_max_diagonal_drift": "comparison.losatp-collinear",
    "--collinear_max_paralog_links_per_orthogroup": "comparison.losatp-collinear",
    "--collinear_max_unit_gap": "comparison.losatp-collinear",
    "--collinear_merge_orientation": "comparison.losatp-collinear",
    "--collinear_min_anchors": "comparison.losatp-collinear",
    "--collinear_search_scope": "comparison.losatp-collinear",
    "--collinear_unit_mode": "comparison.losatp-collinear",
    "--conservation_blast": "comparison.circular-rings",
    "--conservation_color_scheme": "comparison.circular-rings",
    "--conservation_colors": "comparison.circular-rings",
    "--conservation_labels": "comparison.circular-rings",
    "--conservation_reference": "comparison.circular-rings",
    "--conservation_ring_gap": "comparison.circular-rings",
    "--conservation_ring_width": "comparison.circular-rings",
    "--conservation_sequence": "comparison.circular-rings",
    "--conservation_genome": "comparison.circular-rings",
    "--gc": "track.gc-content",
    "--no-gc": "track.gc-content",
    "--gc_content_large_tick_interval": "track.gc-content",
    "--gc_content_max_percent": "track.gc-content",
    "--gc_content_min_percent": "track.gc-content",
    "--gc_content_mode": "track.gc-content",
    "--gc_content_radius": "track.gc-content",
    "--gc_content_small_tick_interval": "track.gc-content",
    "--gc_content_tick_font_size": "track.gc-content",
    "--gc_content_axis_labels": "track.gc-content",
    "--gc_content_tick_interval": "track.gc-content",
    "--gc_content_width": "track.gc-content",
    "--gc_height": "track.gc-content",
    "--hide_gc_content_axis": "track.gc-content",
    "--hide_gc_content_ticks": "track.gc-content",
    "--show_gc_content_axis": "track.gc-content",
    "--show_gc_content_ticks": "track.gc-content",
    "--step": "track.gc-content",
    "--window": "track.gc-content",
    "--skew": "track.dinucleotide-skew",
    "--no-skew": "track.dinucleotide-skew",
    "--nt": "track.dinucleotide-skew",
    "--gc_skew_radius": "track.dinucleotide-skew",
    "--gc_skew_width": "track.dinucleotide-skew",
    "--depth_color": "track.depth-numeric",
    "--depth_height": "track.depth-numeric",
    "--depth_large_tick_interval": "track.depth-numeric",
    "--depth_log_scale": "track.depth-numeric",
    "--depth_max": "track.depth-numeric",
    "--depth_min": "track.depth-numeric",
    "--depth_small_tick_interval": "track.depth-numeric",
    "--depth_step": "track.depth-numeric",
    "--depth_tick_font_size": "track.depth-numeric",
    "--depth_axis_labels": "track.depth-numeric",
    "--depth_track": "track.depth-numeric",
    "--depth_track_color": "track.depth-numeric",
    "--depth_track_height": "track.depth-numeric",
    "--depth_track_label": "track.depth-numeric",
    "--depth_track_large_tick_interval": "track.depth-numeric",
    "--depth_track_small_tick_interval": "track.depth-numeric",
    "--depth_track_tick_font_size": "track.depth-numeric",
    "--depth_width": "track.depth-numeric",
    "--depth_window": "track.depth-numeric",
    "--feature_density_table": "track.depth-numeric",
    "--feature_density_track": "track.depth-numeric",
    "--hide_depth_axis": "track.depth-numeric",
    "--hide_depth_ticks": "track.depth-numeric",
    "--no_depth_log_scale": "track.depth-numeric",
    "--share_depth_axis": "track.depth-numeric",
    "--show_depth_axis": "track.depth-numeric",
    "--show_depth_ticks": "track.depth-numeric",
    "--annotation_table": "track.region-annotations",
    "--annotation-table": "track.region-annotations",
    "--circular_track_axis_index": "track.custom-circular-slots",
    "--circular_track_order": "track.custom-circular-slots",
    "--circular_track_slot": "track.custom-circular-slots",
    "--circular_track_table": "track.custom-circular-slots",
    "--linear_track_axis_index": "track.custom-linear-slots",
    "--linear_track_order": "track.custom-linear-slots",
    "--linear_track_slot": "track.custom-linear-slots",
    "--palette": "presentation.palettes",
    "--default_colors": "presentation.palettes",
    "--table": "presentation.feature-colors",
    "--qualifier_priority": "presentation.qualifier-priority",
    "--label_blacklist": "presentation.label-filters",
    "--label_whitelist": "presentation.label-filters",
    "--label_table": "presentation.label-overrides",
    "--features": "presentation.visibility",
    "--feature_override_table": "presentation.visibility",
    "--feature_visibility_table": "presentation.visibility",
    "--arrow_head_length_ratio": "presentation.shapes",
    "--arrow_shaft_width_ratio": "presentation.shapes",
    "--feature_shape": "presentation.shapes",
    "--block_stroke_color": "presentation.strokes",
    "--block_stroke_width": "presentation.strokes",
    "--line_stroke_color": "presentation.strokes",
    "--line_stroke_width": "presentation.strokes",
    "--feature_overlap_tolerance_bp": "presentation.overlaps",
    "--feature_placement_table": "presentation.overlaps",
    "--resolve_overlaps": "presentation.overlaps",
    "--legend_rows": "interactive.legend-editing",
    "--output": "output.static-svg",
    "--format": "output.static-svg",
    "--overwrite": "output.static-svg",
    "--dpi": "output.png",
    "--session": "program.session-round-trip",
    "--save_session": "program.session-round-trip",
    "--session_output": "program.session-round-trip",
}

UNDECLARED_CLI_OPTIONS: frozenset[str] = frozenset(
    {
        # Figure size.
        "--canvas_padding",
        "--figure_height",
        "--figure_width",
        # Fonts.
        "--font_family",
        "--font_file",
        "--no_font_embedding",
        # Drawings of a Session.
        "--drawing",
        "--list_drawings",
        # Linear wrap.
        "--wrap_length",
        # Feature label placement and style.
        "--circular_label_spacing",
        "--inner_label_x_radius_offset",
        "--inner_label_y_radius_offset",
        "--italic_gene_labels",
        "--label_font_size",
        "--label_line",
        "--label_placement",
        "--label_protection",
        "--label_rendering",
        "--label_rotation",
        "--labels",
        "--linear_label_spacing",
        "--outer_label_x_radius_offset",
        "--outer_label_y_radius_offset",
        "--show_labels",
        # Record alignment by a feature.
        "--anchor_feature",
        "--similarity_alignment_feature",
        # Record origin, coordinates and topology.
        "--display_start_coordinate",
        "--origin_feature",
        "--record_topology",
        # LOSAT and BLAST+ runs on the command line (the comparison ids name
        # the browser runs and the LOSATP displays).
        "--conservation_losat_gencode",
        "--losat",
        "--losat_bin",
        "--losat_gencode",
        "--losat_output_dir",
        "--losat_threads",
        "--losatn_task",
        "--ncbi_blast_bin",
    }
)


def test_every_cli_option_is_declared_or_listed() -> None:
    capabilities = {
        capability
        for group in json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))[
            "capability_groups"
        ].values()
        for capability in group
    }
    options = set().union(*(_live_help_options(command) for command in COMMANDS)) - {
        "--help"
    }

    undeclared = sorted(
        options - CLI_OPTION_CAPABILITIES.keys() - UNDECLARED_CLI_OPTIONS
    )
    assert not undeclared, (
        f"{undeclared} belong to no capability: add each to CLI_OPTION_CAPABILITIES with the capability "
        "whose public owner documents it, or add a capability with its scenario to docs/scenarios/manifest.json"
    )
    unknown = sorted(
        f"{option}: {capability}"
        for option, capability in CLI_OPTION_CAPABILITIES.items()
        if capability not in capabilities
    )
    assert not unknown, (
        f"capabilities missing from docs/scenarios/manifest.json: {unknown}"
    )
    assert not CLI_OPTION_CAPABILITIES.keys() & UNDECLARED_CLI_OPTIONS, (
        "an option is both declared and listed"
    )
    gone = sorted((CLI_OPTION_CAPABILITIES.keys() | UNDECLARED_CLI_OPTIONS) - options)
    assert not gone, f"{gone} are no longer options of {COMMANDS}: remove them here"
