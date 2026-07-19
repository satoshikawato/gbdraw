#!/usr/bin/env python3
from __future__ import annotations

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

from gbdraw.features.ids import compute_feature_hash_from_parts
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg
from gbdraw.session_io import load_session, write_session_json


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
    "--line_stroke_width 2 --suppress_skew -p default --track_type tuckin "
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
    "-o vibrio-harveyi-group-collinear -f svg"
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
    "--suppress_gc --suppress_skew --track_type spreadout -l left --feature_width 10 "
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
        tags=("Linear", "Multi-record", "Collinear groups", "LOSAT", "Static SVG"),
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
        sync_result_svg=False,
        interactive_svg=False,
        compressed_session=True,
        source_note=GZIP_SESSION_SOURCE_NOTE,
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


def _session_result_svg(session: dict[str, Any], example: GallerySessionExample) -> str:
    results = session.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{example.session_path} does not contain a results array.")
    for result in results:
        if not isinstance(result, dict):
            continue
        content = result.get("content")
        if isinstance(content, str) and "<svg" in content:
            return content
    raise ValueError(f"{example.session_path} does not contain a generated SVG result.")


def _session_interactive_context(session: dict[str, Any]) -> InteractiveSvgContext:
    feature_state = session.get("features") if isinstance(session.get("features"), dict) else {}
    editor_state = session.get("editorState") if isinstance(session.get("editorState"), dict) else {}
    orthogroup_state = (
        session.get("orthogroupState") if isinstance(session.get("orthogroupState"), dict) else {}
    )
    ui = session.get("ui") if isinstance(session.get("ui"), dict) else {}
    config = session.get("config") if isinstance(session.get("config"), dict) else {}
    legend = editor_state.get("legend") if isinstance(editor_state.get("legend"), dict) else {}

    features = feature_state.get("extractedFeatures")
    orthogroups = orthogroup_state.get("groups")
    legend_entries = legend.get("entries")
    current_colors = ui.get("appliedPaletteColors") or config.get("colors")
    return InteractiveSvgContext(
        features=features if isinstance(features, list) else (),
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


def _validate_source_feature_ids(
    example: GallerySessionExample,
    session: dict[str, Any],
    source: str,
) -> None:
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
    missing_ids = sorted(rendered_ids - metadata_ids)
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


def _read_or_create_source_svg(example: GallerySessionExample, session: dict[str, Any]) -> str:
    if example.source_svg_path.exists():
        source = example.source_svg_path.read_text(encoding="utf-8")
    else:
        source = _session_result_svg(session, example)
    migrated = _migrate_legacy_multipart_feature_ids(source, session)
    if example.interactive_svg:
        _validate_source_feature_ids(example, session, migrated)
    else:
        ET.fromstring(migrated)
    if migrated != source or not example.source_svg_path.exists():
        example.source_svg_path.write_text(migrated, encoding="utf-8")
    return migrated


def _write_gallery_svg(
    example: GallerySessionExample,
    session: dict[str, Any],
    source: str,
) -> None:
    output = (
        enrich_svg(source, context=_session_interactive_context(session))
        if example.interactive_svg
        else source
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


def _render_thumbnail(example: GallerySessionExample) -> None:
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


def prepare_gallery_assets() -> list[dict[str, object]]:
    EXAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
    SESSION_ROOT.mkdir(parents=True, exist_ok=True)
    SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
    THUMBNAIL_ROOT.mkdir(parents=True, exist_ok=True)
    _remove_stale_assets()

    payload: list[dict[str, object]] = []
    for example in EXAMPLES:
        session = _load_session(example)
        source = _read_or_create_source_svg(example, session)
        if example.sync_result_svg:
            _sync_session_result_svg(example, session, source)
        _write_gallery_svg(example, session, source)
        _validate_source_assets(example)
        _render_thumbnail(example)

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

    payload.sort(key=lambda entry: int(entry["displayOrder"]))

    (GALLERY_ROOT / "examples.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return payload


def main() -> int:
    payload = prepare_gallery_assets()
    print(f"Prepared {len(payload)} interactive gallery examples in {GALLERY_ROOT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
