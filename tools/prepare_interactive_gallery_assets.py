#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import gzip
import io
import json
import re
import shlex
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cairosvg
from PIL import Image, ImageDraw, ImageFont

from gbdraw.exceptions import GbdrawError
from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg
from gbdraw.session_io import load_session, write_session_json
from gbdraw.web_support.feature_catalog import (
    FEATURE_CATALOG_SCHEMA,
    canonical_catalog_sequence_sources,
    materialize_catalog_nucleotide_sequence,
    select_feature_catalog_item,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
EXAMPLE_ROOT = GALLERY_ROOT / "examples"
SESSION_ROOT = GALLERY_ROOT / "sessions"
SOURCE_ROOT = GALLERY_ROOT / "sources"
THUMBNAIL_ROOT = GALLERY_ROOT / "thumbnails"
GZIP_SESSION_SOURCE_NOTE = (
    "The complete gzip-compressed Session JSON and generated SVG output are stored "
    "with the gallery assets."
)
INTERACTIVE_SVG_HARD_LIMIT = 41_943_040
DECODED_METADATA_HARD_LIMIT = 200_000_000
DECODED_METADATA_REGRESSION_CEILING = 220_000_000

GENOME_SUFFIXES = (".gb", ".gbk", ".gbff")
_RENDERED_RECORD_SUFFIX_RE = re.compile(r"_record_\d+$")
_SVG_PART_SUFFIX_RE = re.compile(r"__part\d+$")
_SVG_FEATURE_ID_RE = re.compile(r'data-gbdraw-feature-id=["\']([^"\']+)["\']')


@dataclass(frozen=True)
class GallerySessionExample:
    id: str
    title: str
    tags: tuple[str, ...]
    description: str
    workflow: str
    input_summary: str
    display_order: int
    command_kind: str
    command_note: str
    feature_sources: tuple[str, ...] = ()
    sync_result_svg: bool = True
    interactive_svg: bool = True
    compressed_metadata: bool = False
    compressed_session: bool = False
    interactive_step: str = ""
    source_note: str = "Session JSON and generated SVG output are stored with the gallery assets."
    command: str = ""

    @property
    def session_path(self) -> Path:
        suffix = ".gbdraw-session.json.gz" if self.compressed_session else ".gbdraw-session.json"
        return SESSION_ROOT / f"{self.id}{suffix}"

    @property
    def session_ref(self) -> str:
        return f"./sessions/{self.session_path.name}"

    @property
    def source_svg_path(self) -> Path:
        return SOURCE_ROOT / f"{self.id}.svg"

    @property
    def output_svg_path(self) -> Path:
        return self.gallery_svg_path

    @property
    def gallery_svg_path(self) -> Path:
        return EXAMPLE_ROOT / f"{self.id}.svg"

    @property
    def gallery_svg_ref(self) -> str:
        return f"./examples/{self.id}.svg"

    @property
    def thumbnail_path(self) -> Path:
        return THUMBNAIL_ROOT / f"{self.id}.webp"

    @property
    def thumbnail_ref(self) -> str:
        return f"./thumbnails/{self.id}.webp"


HMMTDNA_ATSKEW_COMMAND = (
    "gbdraw circular -o HmmtDNA_ATskew --species '<i>Homo sapiens</i>' "
    "-k CDS,rRNA,tRNA,tmRNA,ncRNA,repeat_region -p ajisai --window 500 --step 50 "
    "--definition_font_size 28 --label_font_size 20 --track_type middle -l left --labels "
    "--qualifier_priority HmmtDNA_qualifier_priority.tsv --circular_track_axis_index 0 "
    "--circular_track_slot features:features@lane_direction=split "
    "--circular_track_slot gc_content:dinucleotide_content@w=0.1 "
    "--circular_track_slot gc_skew:dinucleotide_skew@w=0.1 "
    "--circular_track_slot "
    "'a_skew_2:dinucleotide_skew@w=0.1,nt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew' "
    "--circular_track_slot ticks:ticks@tick_label_layout=label_in_tick_out "
    "--gbk HmmtDNA.gbk -f interactive_svg"
)

HMMTDNA_BASIC_COMMAND = (
    "gbdraw circular --gbk HmmtDNA.gbk -o HmmtDNA_basic_circular "
    "-f interactive_svg --separate_strands --track_type middle --labels out "
    "--species '<i>Homo sapiens</i>'"
)

TOBACCO_CHLOROPLAST_COMMAND = (
    "gbdraw circular --gbk NC_001879.gbk "
    "--annotation_table nicotiana-tabacum-regions.tsv --separate_strands "
    "-k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin "
    "-t chloroplast_specific_table.tsv --qualifier_priority qualifier_priority.tsv "
    "--block_stroke_width 1 --block_stroke_color black --axis_stroke_width 3 "
    "--line_stroke_width 2 --no-skew -p default --track_type tuckin "
    "--labels both --label_placement radial --outer_label_x_radius_offset 0.90 "
    "--outer_label_y_radius_offset 0.90 --inner_label_x_radius_offset 0.975 "
    "--inner_label_y_radius_offset 0.975 --species '<i>Nicotiana tabacum</i>' "
    "--definition_font_size 28 --legend upper_left "
    "--circular_track_slot 'features:features@side=overlay,lane_direction=split' "
    "--circular_track_slot "
    "'plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,"
    "show_labels=true,padding_px=1,overflow=compress,inner_gap_px=1,outer_gap_px=1' "
    "--circular_track_slot 'gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08' "
    "-o tobacco-chloroplast -f interactive_svg"
)

LAMBDA_BASIC_COMMAND = (
    "gbdraw linear --gbk NC_001416.gb -o lambda_basic_linear -f interactive_svg "
    "--separate_strands --show_labels all --scale_style ruler --legend left"
)

BGC_COMMAND = (
    "gbdraw linear --protein_blastp_mode orthogroup -f interactive_svg "
    "--gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk "
    "-k CDS,rRNA,tRNA,tmRNA,ncRNA,repeat_region -p orange "
    "-d BGC0000708-BGC0000713_default_colors.tsv "
    "-t BGC0000708-BGC0000713_specific_colors.tsv "
    "--qualifier_priority BGC0000708-BGC0000713_qualifier_priority.tsv "
    "--show_labels first --label_placement above_feature --label_rotation 45 "
    "--pairwise_match_style curve --scale_style ruler "
    "--plot_title 'Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.' "
    "--keep_definition_left_aligned --identity 30 --block_stroke_width 2 "
    "--block_stroke_color '#262626' --line_stroke_width 2 --axis_stroke_width 5 "
    "--legend_box_size 20 --legend_font_size 20 --label_font_size 18 --feature_height 75 "
    "--ruler_label_font_size 20 --definition_line_style name:size=20,weight=bold "
    "--definition_line_style subtitle:size=20 "
    "--definition_line_style 'accession:size=20,color=#7b7c7d' "
    "--definition_line_style 'length:size=20,color=#7b7c7d' "
    "-l bottom -o BGC0000708-BGC0000713"
)

VIBRIO_HARVEYI_GROUP_COMMAND = (
    "gbdraw linear --records_table examples/vibrio-harveyi-group-linear-records.tsv "
    "--linear_record_gap 48 --track_layout above --scale_style ruler --ruler_on_axis "
    "--scale_interval 750000 --separate_strands --hide_accession --hide_length "
    "--definition_font_size 16 --definition_line_style 'name:size=18,weight=bold' "
    "--definition_line_style 'subtitle:size=16' --keep_definition_left_aligned "
    "--protein_blastp_mode collinear --collinear_search_scope adjacent "
    "--protein_blastp_candidate_limit 5 --collinear_min_anchors 3 "
    "--collinear_max_unit_gap 2 --collinear_max_diagonal_drift 2 "
    "--collinear_color_mode orientation_identity --pairwise_match_style curve "
    "--losatp_threads 16 --feature_shape CDS=rectangle --feature_shape rRNA=rectangle "
    "--feature_shape tRNA=rectangle --feature_shape tmRNA=rectangle "
    "--feature_shape ncRNA=rectangle --feature_shape misc_RNA=rectangle "
    "--feature_shape repeat_region=rectangle --block_stroke_width 0 "
    "--axis_stroke_width 2 --line_stroke_width 1 --legend bottom "
    "-o vibrio-harveyi-group-collinear -f interactive_svg"
)

WSSV_CONSERVATION_LABELS = (
    "CN01 WSSV-TW WSSV-CN WSSV-TH JP01A JP01B Pc2020 E1 0722-1 CN03 "
    "CN04 WSSV-AU EU129 GCF7 MES-753 Shantou2019 POMZ1 POMZ4 "
    "MG18PR-0187-N40S Angostura2013"
)
WSSV_CONSERVATION_BLASTS = " ".join(
    f"{label}.circular_conservation.losatn.tsv" for label in WSSV_CONSERVATION_LABELS.split()
)
WSSV_CONSERVATION_COLORS = (
    "'#6e91b7' '#f4a251' '#77b26f' '#e67577' '#8fc4c0' '#f0d369' '#be92b2' "
    "'#ffafb7' '#ae8e7c' '#c6bebb' '#6e91b7' '#f4a251' '#e67577' '#8fc4c0' "
    "'#bcb4ca' '#f0d369' '#be92b2' '#ffafb7' '#ae8e7c' '#c6bebb'"
)
WSSV_COMMAND = (
    "gbdraw circular -o WSSV_genome_comparison --separate_strands "
    "-k CDS,rRNA,tRNA,tmRNA,ncRNA,repeat_region -p royal_gala "
    "--qualifier_priority WSSV_qualifier_priority.tsv --block_stroke_width 1 "
    "--line_stroke_width 2 --legend_box_size 12 --legend_font_size 12 "
    "--no-gc --no-skew --track_type spreadout -l left --feature_width 10 "
    "--outer_label_x_radius_offset 1 --outer_label_y_radius_offset 1 "
    f"--conservation_blast {WSSV_CONSERVATION_BLASTS} "
    "--conservation_reference subject "
    f"--conservation_labels {WSSV_CONSERVATION_LABELS} "
    f"--conservation_colors {WSSV_CONSERVATION_COLORS} "
    "--conservation_ring_width 5 --conservation_ring_gap 2 --bitscore 100 "
    "--evalue 1e-30 --identity 90 --alignment_length 100 --gbk AP027280.gb "
    "-f interactive_svg"
)


EXAMPLES: tuple[GallerySessionExample, ...] = (
    GallerySessionExample(
        id="HmmtDNA_basic_circular",
        title="Human mitochondrial genome: first circular figure",
        tags=("Circular", "Interactive SVG"),
        description="Create a first circular genome figure from one small GenBank record without running a sequence search.",
        workflow="Circular basics",
        input_summary="1 GenBank file",
        display_order=10,
        command_kind="runnable",
        command_note="Download HmmtDNA.gbk from the Files tab, then run this command in the same directory.",
        command=HMMTDNA_BASIC_COMMAND,
    ),
    GallerySessionExample(
        id="lambda_basic_linear",
        title="Lambda phage: first linear figure",
        tags=("Linear", "Interactive SVG"),
        description="Create a labeled linear genome figure and learn strand separation, the ruler, and SVG export.",
        workflow="Linear basics",
        input_summary="1 GenBank file",
        display_order=20,
        command_kind="runnable",
        command_note="Download NC_001416.gb from the Files tab, then run this command in the same directory.",
        command=LAMBDA_BASIC_COMMAND,
    ),
    GallerySessionExample(
        id="HmmtDNA_ATskew",
        title="Human mitochondrial genome (AT skew)",
        tags=("Circular", "Interactive SVG"),
        description="Add an AT skew ring and qualifier-based labels to a compact circular mitochondrial diagram.",
        workflow="Circular quantitative tracks",
        input_summary="1 GenBank + 1 qualifier TSV",
        display_order=30,
        command_kind="runnable",
        command_note="Download both HmmtDNA.gbk and HmmtDNA_qualifier_priority.tsv from Files before running the command.",
        command=HMMTDNA_ATSKEW_COMMAND,
    ),
    GallerySessionExample(
        id="tobacco-chloroplast",
        title="<i>Nicotiana tabacum</i> chloroplast genome regions",
        tags=("Circular", "Interactive SVG"),
        description="Mark LSC, SSC, IRa, and IRb as bracket annotations inside a color-coded chloroplast gene map.",
        workflow="Circular region annotations",
        input_summary="1 GenBank + 3 TSV files",
        display_order=40,
        command_kind="runnable",
        command_note="Download NC_001879.gbk and the three Gallery TSV files, then run the command in the same directory.",
        command=TOBACCO_CHLOROPLAST_COMMAND,
    ),
    GallerySessionExample(
        id="Vnig_TUMSAT-TG-2018",
        title="<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018",
        tags=("Circular", "Multi-record", "Interactive SVG"),
        description="Arrange two chromosomes and four plasmids from one multi-record GenBank file on a shared circular canvas.",
        workflow="Circular multi-record canvas",
        input_summary="1 multi-record GenBank file",
        display_order=50,
        command_kind="runnable",
        command_note="Download the pinned RefSeq assembly named in Files; no sequence search is required.",
        compressed_session=True,
        source_note=GZIP_SESSION_SOURCE_NOTE,
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_collinear",
        title="Hepatoplasmataceae collinear protein-match blocks",
        tags=("Linear", "Collinear groups", "LOSAT", "Interactive SVG"),
        description="Combine compatible protein-match anchors into collinear blocks across five related genomes.",
        workflow="LOSATP collinear blocks",
        input_summary="5 GenBank files",
        display_order=60,
        command_kind="runnable",
        command_note="Download the five accession-pinned GenBank inputs from Files. The command runs LOSATP locally.",
        compressed_session=True,
        source_note=GZIP_SESSION_SOURCE_NOTE,
    ),
    GallerySessionExample(
        id="vibrio-harveyi-group-collinear",
        title="<i>Vibrio</i> Harveyi group multi-record collinearity",
        tags=("Linear", "Multi-record", "Collinear groups", "LOSAT", "Interactive SVG"),
        description="Compare all 11 replicons from five Harveyi-group Vibrio assemblies as five multi-record rows.",
        workflow="Multi-record LOSATP collinear blocks",
        input_summary="5 multi-record GenBank files; 11 replicons",
        display_order=65,
        command_kind="runnable",
        command_note="Run from a source checkout so the records table can read the five GBFF files under tests/test_inputs/.",
        command=VIBRIO_HARVEYI_GROUP_COMMAND,
        feature_sources=(
            "GCF_030060435.1_ASM3006043v1_genomic.gbff",
            "GCF_002021755.1_ASM202175v1_genomic.gbff",
            "GCF_002906475.1_ASM290647v1_genomic.gbff",
            "GCF_000196095.1_ASM19609v1_genomic.gbff",
            "GCF_000354175.2_ASM35417v2_genomic.gbff",
        ),
        compressed_metadata=True,
        compressed_session=True,
        source_note=(
            "The interactive SVG with gzip-compressed metadata, complete gzip-compressed "
            "Session JSON, and uncompressed source figure are stored with the gallery assets."
        ),
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_orthogroup",
        title="Hepatoplasmataceae CDS protein-similarity links",
        tags=("Linear", "Similarity groups", "LOSAT", "Interactive SVG"),
        description="Compare the same five genomes with similarity-group links instead of collinear blocks.",
        workflow="LOSATP similarity groups",
        input_summary="5 GenBank files",
        display_order=70,
        command_kind="runnable",
        command_note="Download the five accession-pinned GenBank inputs from Files. The CLI value remains orthogroup for compatibility.",
        compressed_session=True,
        source_note=GZIP_SESSION_SOURCE_NOTE,
    ),
    GallerySessionExample(
        id="BGC0000708-BGC0000713",
        title="Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.",
        tags=("Linear", "Similarity groups", "LOSAT", "Interactive SVG"),
        description="Compare five biosynthetic gene clusters while preserving antiSMASH categories and concise gene labels.",
        workflow="LOSATP similarity groups and color rules",
        input_summary="5 GenBank + 3 color/label TSV files",
        display_order=80,
        command_kind="runnable",
        command_note="Files provides the five MIBiG records and all three repository-managed TSV files used by the command.",
        command=BGC_COMMAND,
    ),
    GallerySessionExample(
        id="majanivirus_orthogroup",
        title="Majanivirus CDS protein-similarity links",
        tags=("Linear", "Similarity groups", "LOSAT", "Interactive SVG"),
        description="Inspect dense protein-similarity links and product-based feature colors across nine viral genomes.",
        workflow="LOSATP similarity groups",
        input_summary="9 GenBank + 2 color TSV files",
        display_order=90,
        command_kind="runnable",
        command_note="Download the nine accession-pinned records and both repository-managed color tables from Files.",
        compressed_session=True,
        source_note=GZIP_SESSION_SOURCE_NOTE,
    ),
    GallerySessionExample(
        id="WSSV_genome_comparison",
        title="White spot syndrome virus nucleotide-similarity rings",
        tags=("Circular", "LOSAT", "Interactive SVG"),
        description="Inspect one viral reference against 20 prepared assemblies as concentric BLAST/LOSAT comparison rings.",
        workflow="Session-based circular comparison case study",
        input_summary="Bundled session; prepared 20-assembly input set not fully public",
        display_order=100,
        command_kind="provenance",
        command_note="This records the original prepared-input workflow and is not directly runnable from public downloads. Load the bundled session first; Shantou2019 has no recorded public accession.",
        command=WSSV_COMMAND,
    ),
)


def _format_size(num_bytes: int) -> str:
    if num_bytes >= 1024 * 1024:
        return f"{num_bytes / (1024 * 1024):.1f} MB"
    return f"{num_bytes / 1024:.0f} KB"


def _load_session(example: GallerySessionExample) -> dict[str, Any]:
    session = load_session(example.session_path)
    if not isinstance(session, dict):
        raise ValueError(f"{example.session_path} did not contain a JSON object.")
    return session


def _session_cli_invocation(session: dict[str, Any]) -> dict[str, Any]:
    cli = session.get("cliInvocation")
    if isinstance(cli, dict):
        return cli
    config = session.get("config")
    if isinstance(config, dict):
        cli = config.get("cliInvocation")
        if isinstance(cli, dict):
            return cli
    return {}


def _session_raw_args(session: dict[str, Any]) -> list[str]:
    cli = _session_cli_invocation(session)
    args = cli.get("args")
    if isinstance(args, list):
        return [str(arg) for arg in args]

    config = session.get("config")
    if isinstance(config, dict):
        cli_options = config.get("cliOptions")
        if isinstance(cli_options, dict):
            raw_args = cli_options.get("rawArgs")
            if isinstance(raw_args, list):
                return [str(arg) for arg in raw_args]
    return []


def _session_command(session: dict[str, Any]) -> str:
    cli = _session_cli_invocation(session)
    config = session.get("config") if isinstance(session.get("config"), dict) else {}
    ui = session.get("ui") if isinstance(session.get("ui"), dict) else {}
    cli_options = config.get("cliOptions") if isinstance(config.get("cliOptions"), dict) else {}
    mode = str(cli.get("mode") or cli_options.get("mode") or ui.get("mode") or "linear")
    raw_args = _session_raw_args(session)
    if raw_args:
        return shlex.join(["gbdraw", mode, *raw_args])
    return shlex.join(["gbdraw", mode, "--session", "session.gbdraw-session.json"])


def _example_command(example: GallerySessionExample, session: dict[str, Any]) -> str:
    command = example.command or _session_command(session)
    command = command.replace("interactive-svg", "interactive_svg")
    if example.id == "hepatoplasmataceae_collinear":
        command = command.replace("--losatp_threads 32 --protein_blastp_mode collinear --losatp_threads 32", "--losatp_threads 32 --protein_blastp_mode collinear")
    return command


def _session_feature_sources(session: dict[str, Any]) -> list[str]:
    cli = _session_cli_invocation(session)
    seen: set[str] = set()
    sources: list[str] = []

    def add_name(value: object) -> None:
        if not isinstance(value, str) or not value.lower().endswith(GENOME_SUFFIXES):
            return
        if value in seen:
            return
        seen.add(value)
        sources.append(value)

    bindings = cli.get("fileBindings")
    if isinstance(bindings, list):
        for binding in bindings:
            if not isinstance(binding, dict):
                continue
            add_name(binding.get("name"))

    files = session.get("files")
    if isinstance(files, dict):
        c_gb = files.get("c_gb")
        if isinstance(c_gb, dict):
            add_name(c_gb.get("name"))

        linear_seqs = files.get("linearSeqs")
        if isinstance(linear_seqs, list):
            for entry in linear_seqs:
                if not isinstance(entry, dict):
                    continue
                gb_entry = entry.get("gb")
                if isinstance(gb_entry, dict):
                    add_name(gb_entry.get("name"))
                add_name(entry.get("name"))
    return sources


def _session_result(
    session: dict[str, Any],
    example: GallerySessionExample,
) -> tuple[int, dict[str, Any]]:
    results = session.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{example.session_path} does not contain a results array.")
    for result_index, result in enumerate(results):
        if not isinstance(result, dict):
            continue
        content = result.get("content")
        if isinstance(content, str) and "<svg" in content:
            return result_index, result
    raise ValueError(f"{example.session_path} does not contain a generated SVG result.")


def _session_result_svg(session: dict[str, Any], example: GallerySessionExample) -> str:
    return str(_session_result(session, example)[1]["content"])


def _session_feature_catalog(
    session: dict[str, Any],
) -> dict[str, object] | None:
    editor_state = session.get("editorState")
    catalog = (
        editor_state.get("featureCatalog")
        if isinstance(editor_state, dict)
        else None
    )
    if catalog is None:
        return None
    if (
        not isinstance(catalog, dict)
        or catalog.get("schema") != FEATURE_CATALOG_SCHEMA
        or not isinstance(catalog.get("items"), list)
    ):
        raise ValueError("Session contains an invalid schema-3 feature catalog.")
    return catalog


def _session_catalog_item(
    session: dict[str, Any],
    *,
    result_index: int,
    result_name: str,
) -> dict[str, object] | None:
    catalog = _session_feature_catalog(session)
    if catalog is None:
        return None
    try:
        return select_feature_catalog_item(
            catalog,
            result_index=result_index,
            result_name=result_name,
        )
    except GbdrawError as exc:
        raise ValueError(str(exc)) from exc


def _session_interactive_context(
    session: dict[str, Any],
    *,
    result_index: int = 0,
    result_name: str | None = None,
) -> InteractiveSvgContext:
    catalog = _session_feature_catalog(session)
    if catalog is not None:
        results = session.get("results")
        if result_name is None and isinstance(results, list):
            if 0 <= result_index < len(results) and isinstance(
                results[result_index], dict
            ):
                result_name = str(results[result_index].get("name") or "")
        item = _session_catalog_item(
            session,
            result_index=result_index,
            result_name=str(result_name or ""),
        )
        assert item is not None
        return InteractiveSvgContext(
            features=tuple(
                feature
                for feature in item["features"]
                if isinstance(feature, dict)
            ),
            biological_features=tuple(
                feature
                for feature in item["biologicalFeatures"]
                if isinstance(feature, dict)
            ),
            orthogroups=tuple(
                group
                for group in item["orthogroups"]
                if isinstance(group, dict)
            ),
            annotations=tuple(
                annotation
                for annotation in item["annotations"]
                if isinstance(annotation, dict)
            ),
            sequence_sources=tuple(
                source
                for source in item.get("sequenceSources", [])
                if isinstance(source, dict)
            ),
            record_keys=tuple(str(key) for key in item["recordKeys"]),
        )

    feature_state = session.get("features") if isinstance(session.get("features"), dict) else {}
    editor_state = session.get("editorState") if isinstance(session.get("editorState"), dict) else {}
    orthogroup_state = (
        session.get("orthogroupState") if isinstance(session.get("orthogroupState"), dict) else {}
    )
    ui = session.get("ui") if isinstance(session.get("ui"), dict) else {}
    config = session.get("config") if isinstance(session.get("config"), dict) else {}
    legend = editor_state.get("legend") if isinstance(editor_state.get("legend"), dict) else {}

    features = feature_state.get("extractedFeatures")
    biological_features = feature_state.get("biologicalFeatures")
    orthogroups = orthogroup_state.get("groups")
    legend_entries = legend.get("entries")
    current_colors = ui.get("appliedPaletteColors") or config.get("colors")
    return InteractiveSvgContext(
        features=features if isinstance(features, list) else (),
        biological_features=(
            biological_features
            if isinstance(biological_features, list)
            else features
            if isinstance(features, list)
            else ()
        ),
        orthogroups=orthogroups if isinstance(orthogroups, list) else (),
        legend_entries=legend_entries if isinstance(legend_entries, list) else (),
        current_colors=current_colors if isinstance(current_colors, dict) else {},
    )


def _stable_feature_id(feature: dict[str, Any]) -> str:
    stable_id = str(
        feature.get("stable_svg_id") or feature.get("stable_feature_id") or ""
    ).strip()
    return stable_id or _RENDERED_RECORD_SUFFIX_RE.sub(
        "", str(feature.get("svg_id") or "").strip()
    )


def _legacy_multipart_feature_id(feature: dict[str, Any]) -> str:
    parts = feature.get("location_parts")
    if not isinstance(parts, list) or len(parts) < 2 or not isinstance(parts[0], dict):
        return ""
    first = parts[0]
    strand = {"+": 1, "-": -1}.get(str(first.get("strand") or "").strip())
    try:
        return compute_feature_hash_from_parts(
            str(feature.get("type") or ""),
            int(first["start"]),
            int(first["end"]),
            strand,
            record_id=str(feature.get("record_id") or "") or None,
        )
    except (KeyError, TypeError, ValueError):
        return ""


def _legacy_multipart_feature_aliases(session: dict[str, Any]) -> dict[str, str]:
    feature_state = session.get("features") if isinstance(session.get("features"), dict) else {}
    features = feature_state.get("extractedFeatures")
    candidates: dict[str, set[str]] = {}
    for feature in (features if isinstance(features, list) else ()):
        if not isinstance(feature, dict):
            continue
        legacy_id = _legacy_multipart_feature_id(feature)
        stable_id = _stable_feature_id(feature)
        if legacy_id and stable_id and legacy_id != stable_id:
            candidates.setdefault(legacy_id, set()).add(stable_id)
    return {
        legacy_id: next(iter(stable_ids))
        for legacy_id, stable_ids in candidates.items()
        if len(stable_ids) == 1
    }


def _migrate_legacy_multipart_feature_ids(source: str, session: dict[str, Any]) -> str:
    migrated = source
    for legacy_id, stable_id in _legacy_multipart_feature_aliases(session).items():
        migrated = migrated.replace(legacy_id, stable_id)
    return migrated


def _catalog_sequence_sources(
    item: dict[str, object],
) -> tuple[dict[str, object], ...]:
    raw_sources = item.get("sequenceSources")
    return (
        tuple(
            source if isinstance(source, dict) else {}
            for source in raw_sources
        )
        if isinstance(raw_sources, list)
        else ()
    )


def _catalog_nucleotide_sequence(
    feature: dict[str, Any],
    item: dict[str, object],
    *,
    canonical_sources: tuple[str | None, ...] | None = None,
) -> str:
    sources = _catalog_sequence_sources(item)
    source_index = feature.get("sequenceSourceIndex")
    inline_sequence = str(
        feature.get("nucleotide_sequence")
        or feature.get("nucleotideSequence")
        or ""
    ).strip()
    if not inline_sequence and type(source_index) is int:
        record_keys = item.get("recordKeys")
        record_key = str(feature.get("recordKey") or "").strip()
        normalized_record_keys = (
            [str(value or "").strip() for value in record_keys]
            if isinstance(record_keys, list)
            else []
        )
        try:
            expected_record_index = normalized_record_keys.index(record_key)
        except ValueError:
            expected_record_index = -1
        source = (
            sources[source_index]
            if 0 <= source_index < len(sources)
            else {}
        )
        if (
            str(source.get("origin") or "").strip()
            not in {"linear-record", "circular-reference"}
            or type(source.get("recordIndex")) is not int
            or source.get("recordIndex") != expected_record_index
        ):
            return ""
    return materialize_catalog_nucleotide_sequence(
        feature,
        sources,
        canonical_sources=canonical_sources,
    )


def _catalog_location_parts(feature: dict[str, Any]) -> list[dict[str, object]]:
    raw_parts = feature.get("location_parts", feature.get("locationParts"))
    if isinstance(raw_parts, list) and raw_parts:
        return [
            dict(part)
            for part in raw_parts
            if isinstance(part, dict)
        ]
    start = feature.get("start")
    end = feature.get("end")
    if type(start) is not int or type(end) is not int:
        return []
    return [
        {
            "start": start,
            "end": end,
            "strand": feature.get("strand", ""),
            "display": f"{start + 1}..{end}",
        }
    ]


def _validate_source_feature_ids(
    example: GallerySessionExample,
    session: dict[str, Any],
    source: str,
) -> None:
    item = None
    if _session_feature_catalog(session) is not None:
        result_index, result = _session_result(session, example)
        item = _session_catalog_item(
            session,
            result_index=result_index,
            result_name=str(result.get("name") or ""),
        )
    if item is not None:
        canonical_sources = tuple(
            source.get("sequence")
            if isinstance(source.get("sequence"), str)
            else None
            for source in _catalog_sequence_sources(item)
        )
        biological_by_key = {
            (
                str(feature.get("recordKey") or "").strip(),
                str(feature.get("biologicalFeatureId") or "").strip(),
            ): feature
            for feature in item["biologicalFeatures"]
            if isinstance(feature, dict)
        }
        rendered_by_id = {
            str(feature.get("svgId") or "").strip(): feature
            for feature in item["features"]
            if isinstance(feature, dict)
            and str(feature.get("svgId") or "").strip()
        }
        metadata_ids = set(rendered_by_id)
        metadata_ids.update(
            str(feature.get("stableFeatureId") or "").strip()
            for feature in biological_by_key.values()
            if str(feature.get("stableFeatureId") or "").strip()
        )
        rendered_ids = {
            _SVG_PART_SUFFIX_RE.sub("", match)
            for match in _SVG_FEATURE_ID_RE.findall(source)
        }
        missing_ids = sorted(
            rendered_id
            for rendered_id in rendered_ids
            if rendered_id not in metadata_ids
            and _RENDERED_RECORD_SUFFIX_RE.sub("", rendered_id)
            not in metadata_ids
        )
        if missing_ids:
            preview = ", ".join(missing_ids[:5])
            raise ValueError(
                f"{example.id} source SVG contains {len(missing_ids)} "
                f"feature ID(s) without session metadata: {preview}"
            )

        problems: list[str] = []
        for svg_id, rendered in rendered_by_id.items():
            key = (
                str(rendered.get("recordKey") or "").strip(),
                str(rendered.get("biologicalFeatureId") or "").strip(),
            )
            biological = biological_by_key.get(key)
            if biological is None:
                problems.append(f"{svg_id}: unresolved biological reference")
                continue
            location_parts = _catalog_location_parts(biological)
            qualifiers = biological.get("qualifiers")
            if (
                not str(biological.get("record_id") or biological.get("recordId") or "").strip()
                or not str(biological.get("type") or "").strip()
                or not isinstance(biological.get("start"), int)
                or not isinstance(biological.get("end"), int)
                or not isinstance(location_parts, list)
                or not location_parts
                or (qualifiers is not None and not isinstance(qualifiers, dict))
                or not _catalog_nucleotide_sequence(
                    biological,
                    item,
                    canonical_sources=canonical_sources,
                )
            ):
                problems.append(f"{svg_id}: incomplete popup details")
        if problems:
            preview = "; ".join(problems[:5])
            raise ValueError(
                f"{example.id} session contains {len(problems)} feature "
                f"metadata problem(s): {preview}"
            )
        return

    feature_state = session.get("features") if isinstance(session.get("features"), dict) else {}
    features = feature_state.get("extractedFeatures")
    metadata_ids = {
        candidate
        for feature in (features if isinstance(features, list) else ())
        if isinstance(feature, dict)
        for candidate in (
            str(feature.get("svg_id") or "").strip(),
            _stable_feature_id(feature),
        )
        if candidate
    }
    rendered_ids = {
        _SVG_PART_SUFFIX_RE.sub("", match)
        for match in _SVG_FEATURE_ID_RE.findall(source)
    }
    missing_ids = sorted(
        rendered_id
        for rendered_id in rendered_ids
        if rendered_id not in metadata_ids
        and _RENDERED_RECORD_SUFFIX_RE.sub("", rendered_id) not in metadata_ids
    )
    if missing_ids:
        preview = ", ".join(missing_ids[:5])
        raise ValueError(
            f"{example.id} source SVG contains {len(missing_ids)} feature ID(s) "
            f"without session metadata: {preview}"
        )

    incomplete_ids: list[str] = []
    for feature in (features if isinstance(features, list) else ()):
        if not isinstance(feature, dict):
            continue
        feature_id = str(feature.get("svg_id") or "").strip()
        if not feature_id:
            continue
        complete = (
            bool(str(feature.get("record_id") or "").strip())
            and bool(str(feature.get("type") or "").strip())
            and isinstance(feature.get("start"), int)
            and isinstance(feature.get("end"), int)
            and isinstance(feature.get("location_parts"), list)
            and bool(feature["location_parts"])
            and isinstance(feature.get("qualifiers"), dict)
            and "nucleotide_sequence" in feature
            and "amino_acid_sequence" in feature
        )
        if not complete:
            incomplete_ids.append(feature_id)
    if incomplete_ids:
        preview = ", ".join(incomplete_ids[:5])
        raise ValueError(
            f"{example.id} session contains {len(incomplete_ids)} feature metadata "
            f"record(s) without popup details: {preview}"
        )


def _metadata_record_index(value: object) -> int | None:
    try:
        return int(value) if value is not None and value != "" else None
    except (TypeError, ValueError):
        return None


def _metadata_stable_feature_id(value: dict[str, Any]) -> str:
    stable_id = str(
        value.get("stable_feature_svg_id")
        or value.get("stableFeatureSvgId")
        or value.get("stable_feature_id")
        or value.get("stable_svg_id")
        or value.get("feature_svg_id")
        or value.get("featureSvgId")
        or value.get("svg_id")
        or ""
    ).strip()
    return _RENDERED_RECORD_SUFFIX_RE.sub("", stable_id)


def _validate_interactive_orthogroup_payload(
    example: GallerySessionExample,
    payload: dict[str, Any],
) -> None:
    if payload.get("schema") == FEATURE_CATALOG_SCHEMA:
        items = payload.get("items")
        if not isinstance(items, list) or len(items) != 1 or not isinstance(
            items[0], dict
        ):
            raise ValueError(
                f"{example.id} interactive metadata must contain one catalog item"
            )
        payload = items[0]
    if "biologicalFeatures" in payload:
        catalog_item = payload
        canonical_sources = canonical_catalog_sequence_sources(
            _catalog_sequence_sources(catalog_item)
        )
        if any(source is None for source in canonical_sources):
            raise ValueError(
                f"{example.id} interactive metadata contains an invalid "
                "nucleotide sequence source"
            )
        biological_features = [
            feature
            for feature in payload.get("biologicalFeatures", [])
            if isinstance(feature, dict)
        ]
        biological_by_key = {
            (
                str(feature.get("recordKey") or "").strip(),
                str(feature.get("biologicalFeatureId") or "").strip(),
            ): feature
            for feature in biological_features
        }
        rendered_features = [
            feature
            for feature in payload.get("features", [])
            if isinstance(feature, dict)
        ]
        rendered_by_id = {
            str(feature.get("svgId") or "").strip(): feature
            for feature in rendered_features
            if str(feature.get("svgId") or "").strip()
        }
        groups = [
            group
            for group in payload.get("orthogroups", [])
            if isinstance(group, dict)
        ]
        problems: list[str] = []
        for svg_id, rendered in rendered_by_id.items():
            key = (
                str(rendered.get("recordKey") or "").strip(),
                str(rendered.get("biologicalFeatureId") or "").strip(),
            )
            if key not in biological_by_key:
                problems.append(f"rendered feature {svg_id} is unresolved")
        for group in groups:
            group_id = str(
                group.get("id")
                or group.get("orthogroupId")
                or group.get("orthogroup_id")
                or ""
            ).strip()
            members = group.get("members")
            members = (
                [member for member in members if isinstance(member, dict)]
                if isinstance(members, list)
                else []
            )
            if not members:
                problems.append(f"{group_id or '<unnamed>'}: no members")
                continue
            declared_count = group.get(
                "member_count",
                group.get("memberCount"),
            )
            if (
                isinstance(declared_count, int)
                and declared_count != len(members)
            ):
                problems.append(
                    f"{group_id}: declares {declared_count} members, "
                    f"contains {len(members)}"
                )
            for member in members:
                key = (
                    str(member.get("recordKey") or "").strip(),
                    str(member.get("biologicalFeatureId") or "").strip(),
                )
                feature = biological_by_key.get(key)
                if feature is None:
                    problems.append(
                        f"{group_id}: unresolved member {key[1] or '<empty>'}"
                    )
                    continue
                sequence = _catalog_nucleotide_sequence(
                    feature,
                    catalog_item,
                    canonical_sources=canonical_sources,
                )
                if not sequence:
                    problems.append(
                        f"{group_id}: member {key[1]} has no nucleotide sequence"
                    )
        if problems:
            preview = "; ".join(problems[:5])
            raise ValueError(
                f"{example.id} interactive orthogroup metadata is inconsistent "
                f"({len(problems)} problem(s)): {preview}"
            )
        return

    features = [item for item in payload.get("features", []) if isinstance(item, dict)]
    biological_features = [
        item for item in payload.get("biological_features", []) if isinstance(item, dict)
    ]
    if not biological_features:
        biological_features = features
    features_by_group: dict[str, list[dict[str, Any]]] = {}
    for feature in biological_features:
        group_id = str(feature.get("orthogroup_id") or "").strip()
        if group_id:
            features_by_group.setdefault(group_id, []).append(feature)
    groups = [item for item in payload.get("orthogroups", []) if isinstance(item, dict)]
    groups_by_id = {str(group.get("id") or "").strip(): group for group in groups}
    if not groups_by_id and not features_by_group:
        return
    missing_groups = sorted(features_by_group.keys() - groups_by_id.keys())
    if missing_groups:
        preview = ", ".join(missing_groups[:5])
        raise ValueError(
            f"{example.id} interactive metadata is missing {len(missing_groups)} "
            f"orthogroup(s) referenced by features: {preview}"
        )

    rendered_features_by_id = {
        str(feature.get("svg_id") or "").strip(): feature
        for feature in features
        if str(feature.get("svg_id") or "").strip()
    }
    biological_features_by_key: dict[tuple[int | None, str], dict[str, Any]] = {}
    biological_features_by_stable_id: dict[str, list[dict[str, Any]]] = {}
    for feature in biological_features:
        stable_id = _metadata_stable_feature_id(feature)
        if not stable_id:
            continue
        record_index = _metadata_record_index(
            feature.get("record_idx", feature.get("record_index", feature.get("recordIndex")))
        )
        biological_features_by_key[(record_index, stable_id)] = feature
        biological_features_by_stable_id.setdefault(stable_id, []).append(feature)

    def feature_for_member(member: dict[str, Any]) -> dict[str, Any] | None:
        stable_id = _metadata_stable_feature_id(member)
        record_index = _metadata_record_index(
            member.get("record_index", member.get("recordIndex"))
        )
        exact = biological_features_by_key.get((record_index, stable_id))
        if exact is not None:
            return exact
        candidates = biological_features_by_stable_id.get(stable_id, [])
        return candidates[0] if len(candidates) == 1 else None

    problems: list[str] = []
    for group_id, group in sorted(groups_by_id.items()):
        assigned_features = features_by_group.get(group_id, [])
        members = group.get("members")
        members = [item for item in members if isinstance(item, dict)] if isinstance(members, list) else []
        if not members:
            problems.append(f"{group_id}: no members")
            continue
        member_keys: set[tuple[int | None, str]] = set()
        for member in members:
            stable_id = _metadata_stable_feature_id(member)
            record_index = _metadata_record_index(
                member.get("record_index", member.get("recordIndex"))
            )
            member_keys.add((record_index, stable_id))
            feature = feature_for_member(member)
            if feature is None:
                problems.append(f"{group_id}: unresolved member {stable_id or '<empty>'}")
                continue
            if not str(feature.get("nucleotide_sequence") or "").strip():
                problems.append(f"{group_id}: member {stable_id} has no nucleotide sequence")
            rendered_id = str(
                member.get("rendered_feature_svg_id")
                or member.get("renderedFeatureSvgId")
                or ""
            ).strip()
            if rendered_id and rendered_id not in rendered_features_by_id:
                problems.append(
                    f"{group_id}: rendered member {stable_id} points to missing {rendered_id}"
                )
        for feature in assigned_features:
            stable_id = _metadata_stable_feature_id(feature)
            record_index = _metadata_record_index(
                feature.get("record_idx", feature.get("record_index", feature.get("recordIndex")))
            )
            if (record_index, stable_id) not in member_keys:
                problems.append(
                    f"{group_id}: assigned feature {stable_id or '<empty>'} is not a member"
                )

    if problems:
        preview = "; ".join(problems[:5])
        raise ValueError(
            f"{example.id} interactive orthogroup metadata is inconsistent "
            f"({len(problems)} problem(s)): {preview}"
        )


def _validate_session_interactive_orthogroups(
    example: GallerySessionExample,
    session: dict[str, Any],
) -> None:
    catalog = _session_feature_catalog(session)
    if catalog is not None:
        for item in catalog["items"]:
            if not isinstance(item, dict):
                raise ValueError(
                    f"{example.id} session contains an invalid catalog item"
                )
            _validate_interactive_orthogroup_payload(example, item)
        return
    context = _session_interactive_context(session)
    _validate_interactive_orthogroup_payload(
        example,
        {
            "features": list(context.features),
            "biological_features": list(context.biological_features),
            "orthogroups": list(context.orthogroups),
        },
    )


def _validate_gallery_svg_orthogroups(
    example: GallerySessionExample,
    output: str,
) -> None:
    root = ET.fromstring(output)
    metadata = next(
        (
            element
            for element in root.iter()
            if element.get("id") == "gbdraw-interactive-feature-metadata"
        ),
        None,
    )
    if metadata is None:
        raise ValueError(f"Interactive metadata is missing from {example.id}")
    payload = json.loads(_decoded_metadata_bytes(metadata).decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Interactive metadata is invalid in {example.id}")
    _validate_interactive_orthogroup_payload(example, payload)


def _decoded_metadata_bytes(metadata: ET.Element) -> bytes:
    text = metadata.text or "{}"
    if metadata.get("data-encoding") != "gzip-base64":
        return text.encode("utf-8")
    try:
        return gzip.decompress(base64.b64decode("".join(text.split())))
    except (OSError, ValueError) as exc:
        raise ValueError("Interactive metadata is not valid gzip-base64.") from exc


def _interactive_svg_measurements(output: str) -> dict[str, int]:
    root = ET.fromstring(output)
    metadata = next(
        (
            element
            for element in root.iter()
            if element.get("id") == "gbdraw-interactive-feature-metadata"
        ),
        None,
    )
    if metadata is None:
        raise ValueError("Interactive SVG metadata is missing.")
    decoded = _decoded_metadata_bytes(metadata)
    payload = json.loads(decoded.decode("utf-8"))
    items = (
        payload.get("items", [])
        if isinstance(payload, dict)
        and payload.get("schema") == FEATURE_CATALOG_SCHEMA
        else []
    )
    if not isinstance(items, list):
        items = []
    stored_metadata = (
        base64.b64decode("".join((metadata.text or "").split()))
        if metadata.get("data-encoding") == "gzip-base64"
        else (metadata.text or "").encode("utf-8")
    )
    return {
        "totalInteractiveSvgBytes": len(output.encode("utf-8")),
        "compressedMetadataBytes": len(stored_metadata),
        "decodedMetadataBytes": len(decoded),
        "renderedFeatureCount": sum(
            len(item.get("features", []))
            for item in items
            if isinstance(item, dict)
            and isinstance(item.get("features"), list)
        ),
        "biologicalFeatureCount": sum(
            len(item.get("biologicalFeatures", []))
            for item in items
            if isinstance(item, dict)
            and isinstance(item.get("biologicalFeatures"), list)
        ),
    }


def _format_artifact_measurements(measurements: dict[str, int]) -> str:
    return ", ".join(
        f"{key}={value}" for key, value in measurements.items()
    )


def _validate_interactive_svg_size(
    example: GallerySessionExample,
    output: str,
) -> dict[str, int]:
    measurements = _interactive_svg_measurements(output)
    if measurements["totalInteractiveSvgBytes"] >= INTERACTIVE_SVG_HARD_LIMIT:
        raise ValueError(
            f"{example.id} interactive SVG exceeds the "
            f"{INTERACTIVE_SVG_HARD_LIMIT}-byte limit "
            f"({_format_artifact_measurements(measurements)})"
        )
    if measurements["decodedMetadataBytes"] > DECODED_METADATA_HARD_LIMIT:
        raise ValueError(
            f"{example.id} decoded metadata exceeds the "
            f"{DECODED_METADATA_HARD_LIMIT}-byte limit "
            f"({_format_artifact_measurements(measurements)})"
        )
    return measurements


def _read_or_create_source_svg(
    example: GallerySessionExample,
    session: dict[str, Any],
    *,
    refresh_from_session: bool = False,
) -> str:
    existing_source = (
        example.source_svg_path.read_text(encoding="utf-8")
        if example.source_svg_path.exists()
        else None
    )
    if refresh_from_session and example.sync_result_svg:
        source = _session_result_svg(session, example)
    elif existing_source is not None:
        source = existing_source
    else:
        source = _session_result_svg(session, example)
    migrated = _migrate_legacy_multipart_feature_ids(source, session)
    if example.interactive_svg:
        _validate_source_feature_ids(example, session, migrated)
    else:
        ET.fromstring(migrated)
    if migrated != existing_source:
        example.source_svg_path.write_text(migrated, encoding="utf-8")
    return migrated


def _write_gallery_svg(
    example: GallerySessionExample,
    session: dict[str, Any],
    source: str,
) -> None:
    result_index, result = _session_result(session, example)
    result_name = str(result.get("name") or "")
    catalog = _session_feature_catalog(session)
    output = source
    if example.interactive_svg:
        output = enrich_svg(
            source,
            context=(
                None
                if catalog is not None
                else _session_interactive_context(
                    session,
                    result_index=result_index,
                    result_name=result_name,
                )
            ),
            result_index=result_index,
            result_name=result_name,
            feature_catalog=catalog,
        )
    if example.interactive_svg:
        _validate_gallery_svg_orthogroups(example, output)
    if example.compressed_metadata:
        root = ET.fromstring(output)
        metadata = next(
            (
                element
                for element in root.iter()
                if element.get("id") == "gbdraw-interactive-feature-metadata"
            ),
            None,
        )
        if metadata is None:
            raise ValueError(f"Interactive metadata is missing from {example.id}")
        compressed = gzip.compress(
            (metadata.text or "{}").encode("utf-8"),
            compresslevel=9,
            mtime=0,
        )
        metadata.set("data-encoding", "gzip-base64")
        metadata.text = base64.b64encode(compressed).decode("ascii")
        output = ET.tostring(root, encoding="unicode")
    if example.interactive_svg:
        measurements = _validate_interactive_svg_size(example, output)
        if example.compressed_metadata:
            print(
                f"Gallery artifact metrics [{example.id}]: "
                f"{_format_artifact_measurements(measurements)}"
            )
    example.gallery_svg_path.write_text(output, encoding="utf-8")


def _sync_session_result_svg(
    example: GallerySessionExample,
    session: dict[str, Any],
    source: str,
) -> None:
    results = session.get("results")
    if not isinstance(results, list) or not results or not isinstance(results[0], dict):
        return
    result = results[0]
    old_result_name = str(result.get("name") or "")
    changed = False
    if result.get("name") != example.id:
        result["name"] = example.id
        changed = True
    if result.get("content") != source:
        result["content"] = source
        changed = True
    if session.get("title") != example.id:
        session["title"] = example.id
        changed = True
    run_metadata = session.get("runMetadata")
    geometry = (
        run_metadata.get("trackSlotGeometry")
        if isinstance(run_metadata, dict)
        else None
    )
    geometry_records = geometry.get("records") if isinstance(geometry, dict) else None
    if isinstance(geometry_records, list):
        for record in geometry_records:
            if (
                isinstance(record, dict)
                and record.get("resultIndex") == 0
                and record.get("resultName") != example.id
            ):
                record["resultName"] = example.id
                changed = True
    catalog = _session_feature_catalog(session)
    if catalog is not None:
        for item in catalog["items"]:
            if (
                isinstance(item, dict)
                and item.get("resultIndex") == 0
                and str(item.get("resultName") or "") == old_result_name
                and item.get("resultName") != example.id
            ):
                item["resultName"] = example.id
                changed = True
    if changed:
        write_session_json(example.session_path, session)


def _remove_stale_assets() -> None:
    expected_svgs = {f"{example.id}.svg" for example in EXAMPLES}
    expected_thumbnails = {f"{example.id}.webp" for example in EXAMPLES}

    for path in EXAMPLE_ROOT.glob("*.svg"):
        if path.name not in expected_svgs:
            path.unlink()
    for path in SOURCE_ROOT.glob("*.svg"):
        if path.name not in expected_svgs:
            path.unlink()
    for path in THUMBNAIL_ROOT.glob("*.webp"):
        if path.name not in expected_thumbnails:
            path.unlink()


def _render_thumbnail(
    example: GallerySessionExample,
    *,
    allow_placeholder: bool = True,
) -> None:
    source_path = example.source_svg_path if example.source_svg_path.exists() else example.output_svg_path
    try:
        png_bytes = cairosvg.svg2png(
            url=str(source_path),
            output_width=720,
            background_color="white",
        )
        image_rgba = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
        white = Image.new("RGBA", image_rgba.size, "#ffffff")
        white.alpha_composite(image_rgba)
        image = white.convert("RGB")

        image.thumbnail((640, 360), Image.Resampling.LANCZOS)
        thumbnail = Image.new("RGB", (640, 360), "#ffffff")
        left = (640 - image.width) // 2
        top = (360 - image.height) // 2
        thumbnail.paste(image, (left, top))
    except Exception:
        if not allow_placeholder:
            raise
        thumbnail = Image.new("RGB", (640, 360), "#ffffff")
        draw = ImageDraw.Draw(thumbnail)
        draw.rectangle((24, 28, 616, 332), outline="#b9c7ca", width=2)
        draw.text((40, 44), example.title, fill="#17202a", font=ImageFont.load_default())
        draw.text((40, 78), ", ".join(example.tags), fill="#1d6f7a", font=ImageFont.load_default())
    thumbnail.save(example.thumbnail_path, "WEBP", quality=82, method=6)


def _validate_source_assets(example: GallerySessionExample) -> None:
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in (example.session_path, example.output_svg_path, example.source_svg_path)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing gallery source asset(s): {', '.join(missing)}")


def prepare_gallery_assets(*, refresh_sources: bool = False) -> list[dict[str, object]]:
    for example in EXAMPLES:
        _validate_session_interactive_orthogroups(
            example,
            _load_session(example),
        )

    EXAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
    SESSION_ROOT.mkdir(parents=True, exist_ok=True)
    SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
    THUMBNAIL_ROOT.mkdir(parents=True, exist_ok=True)
    _remove_stale_assets()

    payload: list[dict[str, object]] = []
    for example in EXAMPLES:
        session = _load_session(example)
        source = _read_or_create_source_svg(
            example,
            session,
            refresh_from_session=refresh_sources,
        )
        if example.sync_result_svg:
            _sync_session_result_svg(example, session, source)
        _write_gallery_svg(example, session, source)
        _validate_source_assets(example)
        _render_thumbnail(example, allow_placeholder=not refresh_sources)

        entry = {
            "id": example.id,
            "title": example.title,
            "tags": list(example.tags),
            "svg": example.gallery_svg_ref,
            "svgType": "interactive" if example.interactive_svg else "static",
            "session": example.session_ref,
            "thumbnail": example.thumbnail_ref,
            "sourceSession": str(example.session_path.relative_to(REPO_ROOT)),
            "sourceOutput": str(example.output_svg_path.relative_to(REPO_ROOT)),
            "sourceFigure": str(example.source_svg_path.relative_to(REPO_ROOT)),
            "sourceNote": example.source_note,
            "featureSources": list(example.feature_sources) or _session_feature_sources(session),
            "fileSizeLabel": _format_size(example.gallery_svg_path.stat().st_size),
            "command": _example_command(example, session),
            "commandKind": example.command_kind,
            "commandNote": example.command_note,
            "description": example.description,
            "workflow": example.workflow,
            "inputSummary": example.input_summary,
            "displayOrder": example.display_order,
        }
        tutorial_path = GALLERY_ROOT / "tutorials" / f"{example.id}.json"
        if tutorial_path.exists():
            entry["tutorial"] = f"./tutorials/{example.id}.json"
            entry["tutorialStatus"] = "ready"
        if example.interactive_step:
            entry["interactiveStep"] = example.interactive_step
        payload.append(entry)
        del session, source

    payload.sort(key=lambda entry: int(entry["displayOrder"]))

    (GALLERY_ROOT / "examples.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate interactive Gallery SVGs, thumbnails, and examples.json "
            "from the bundled current-version sessions."
        )
    )
    parser.parse_args(argv)
    payload = prepare_gallery_assets()
    print(f"Prepared {len(payload)} interactive gallery examples in {GALLERY_ROOT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
