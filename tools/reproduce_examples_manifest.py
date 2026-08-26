from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib


FigureGroup = Literal["docs", "palettes", "composites", "readme"]
BlastTool = Literal["blastn", "tblastx"]
CompositeKind = Literal["grid", "contact_sheet", "social_preview"]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PALETTES_FILE = PROJECT_ROOT / "gbdraw" / "data" / "color_palettes.toml"

MANUALLY_MANAGED_FIGURES: dict[str, str] = {
    "examples/gbdraw_social_preview.png": (
        "Owner-maintained README artwork; automated reproduction and replacement are prohibited."
    ),
    "gbdraw/web/gallery/media/lambda_basic_linear/manual-02-03-no-comparison.webp": (
        "Browser screenshot; recapture with the Gallery tutorial screenshot tool from the "
        "Lambda session and its declarative No comparison capture contract."
    ),
    "docs/TUTORIALS/images/tutorial-8-interactive-feature-popup.png": (
        "Browser screenshot; recapture after loading the tutorial session and opening a feature popup."
    ),
    "docs/TUTORIALS/images/tutorial-8-interactive-match-popup.png": (
        "Browser screenshot; recapture after loading the tutorial session and opening a match popup."
    ),
    "docs/TUTORIALS/images/tutorial-8-loaded-session.png": (
        "Browser screenshot; recapture from the local web app after loading the documented session."
    ),
    "gbdraw/web/assets/gbdraw-logo-title.png": (
        "Brand asset derived from the adjacent SVG; update and verify both logo files together."
    ),
    "gbdraw/web/gallery/examples/vibrio-harveyi-group-collinear.svg": (
        "Interactive Web Gallery output with compressed metadata; refresh it together with its source figure and compressed session."
    ),
}

UNREFERENCED_FIGURE_RETENTION: dict[str, str] = {
    "examples/ecoli_k12_plot.svg": (
        "Former Quickstart figure retained for one compatibility release after the "
        "Quickstart moved to bundled, version-pinned tutorial fixtures."
    ),
    "examples/quickstart-labeled-rna-features.svg": (
        "Former Quickstart label example retained for one compatibility release after the "
        "runnable recipe moved to its technical documentation owner."
    ),
    "examples/majani.svg": (
        "Linked from the Gallery but not embedded because the full-size SVG is too large "
        "for a Markdown preview."
    ),
    "examples/NC_001879_color.svg": "Archived chloroplast color-table recipe retained for reproducibility.",
    "examples/NC_010162_edelweiss.svg": "Source panel used by the published social-preview composite.",
}


@dataclass(frozen=True)
class SupportAsset:
    filename: str
    content: str


@dataclass(frozen=True)
class FastaPreparation:
    output_filename: str
    source_filename: str


@dataclass(frozen=True)
class BlastPreparation:
    output_filename: str
    tool: BlastTool
    query_filename: str
    subject_filename: str
    outfmt: str = "7"


Preparation = FastaPreparation | BlastPreparation


@dataclass(frozen=True)
class FileArgument:
    flag: str
    filenames: tuple[str, ...]


@dataclass(frozen=True)
class CliRecipe:
    subcommand: Literal["circular", "linear"]
    gbk_files: tuple[str, ...] = ()
    gff_files: tuple[str, ...] = ()
    fasta_files: tuple[str, ...] = ()
    blast_files: tuple[str, ...] = ()
    file_args: tuple[FileArgument, ...] = ()
    extra_args: tuple[str, ...] = ()


@dataclass(frozen=True)
class SessionVariantRecipe:
    source_session_path: str
    feature_shapes: tuple[tuple[str, str], ...] = ()
    config_overrides: tuple[tuple[str, str | float | int | bool], ...] = ()


@dataclass(frozen=True)
class CompositePanel:
    figure_id: str | None = None
    recipe: CliRecipe | None = None
    preparations: tuple[Preparation, ...] = ()
    box: tuple[int, int, int, int] | None = None


@dataclass(frozen=True)
class CompositeRecipe:
    kind: CompositeKind
    panels: tuple[CompositePanel, ...]
    columns: int | None = None
    tile_size: tuple[int, int] | None = None
    gap: int = 0
    padding: int = 0
    background: str = "white"
    canvas_size: tuple[int, int] | None = None


Recipe = CliRecipe | SessionVariantRecipe | CompositeRecipe


@dataclass(frozen=True)
class FigureSpec:
    figure_id: str
    output_path: str
    groups: tuple[FigureGroup, ...]
    required_inputs: tuple[str, ...]
    recipe: Recipe
    support_assets: tuple[str, ...] = ()
    preparations: tuple[Preparation, ...] = ()
    description: str = ""


SUPPORT_ASSETS: dict[str, SupportAsset] = {
    "feature_specific_colors.tsv": SupportAsset(
        filename="feature_specific_colors.tsv",
        content=(
            "CDS\tproduct\twsv.*-like protein\t#47b8f8\tWSSV-like proteins\n"
            "CDS\tproduct\tbaculoviral IAP repeat-containing protein\tyellow\tBIRP\n"
            "CDS\tproduct\ttyrosine recombinase\tred\ttyrosine recombinase\n"
        ),
    ),
    "stx_whitelist.tsv": SupportAsset(
        filename="stx_whitelist.tsv",
        content=(
            "CDS\tgene\tstx1A\n"
            "CDS\tgene\tstx1B\n"
            "CDS\tgene\tstx2A\n"
            "CDS\tgene\tstx2B\n"
        ),
    ),
    "priority.tsv": SupportAsset(
        filename="priority.tsv",
        content="CDS\tgene\n",
    ),
    "qualifier_priority.tsv": SupportAsset(
        filename="qualifier_priority.tsv",
        content="CDS\tgene\n",
    ),
    "label_override.tsv": SupportAsset(
        filename="label_override.tsv",
        content=(
            "LC738868.1\tCDS\tlabel\t^protein gustavus-like protein$\tgustavus-like protein\n"
            "*\t*\tlabel\t^hypothetical protein$\tHP\n"
        ),
    ),
    "linear_records.tsv": SupportAsset(
        filename="linear_records.tsv",
        content=(
            "gbk\trecord_label\trecord_subtitle\trecord_id\tregion\treverse_complement\torder\n"
            "{PROJECT_ROOT}/examples/MjeNMV.gb\tMarsupenaeus japonicus endogenous nimavirus\tGinoza2017\tLC738868.1\t\t0\t1\n"
            "{PROJECT_ROOT}/examples/MelaMJNV.gb\tMelicertus latisulcatus majanivirus\tOkinawa2016\tLC738874.1\t1-160000\t0\t2\n"
            "{PROJECT_ROOT}/examples/PemoMJNVA.gb\tPenaeus monodon majanivirus A\tMikawa2016\tLC738870.1\t1-160000\t1\t3\n"
        ),
    ),
    "linear_multi_records.tsv": SupportAsset(
        filename="linear_multi_records.tsv",
        content=(
            "gbk\trecord_label\trecord_id\torder\trow\tcolumn\n"
            "{PROJECT_ROOT}/examples/MjeNMV.gb\tMjeNMV\tLC738868.1\t1\t1\t1\n"
            "{PROJECT_ROOT}/examples/PemoMJNVA.gb\tPemoMJNVA\tLC738870.1\t2\t1\t2\n"
            "{PROJECT_ROOT}/examples/MelaMJNV.gb\tMelaMJNV\tLC738874.1\t3\t2\t1\n"
            "{PROJECT_ROOT}/examples/PeseMJNV.gb\tPeseMJNV\tLC738873.1\t4\t2\t2\n"
        ),
    ),
    "linear_multi_comparisons.tsv": SupportAsset(
        filename="linear_multi_comparisons.tsv",
        content=(
            "blast\tquery\tsubject\n"
            "{PROJECT_ROOT}/examples/MjeNMV.MelaMJNV.tblastx.out\t#1\t#3\n"
            "{PROJECT_ROOT}/examples/PemoMJNVA.PeseMJNV.tblastx.out\t#2\t#4\n"
        ),
    ),
    "circular_records.tsv": SupportAsset(
        filename="circular_records.tsv",
        content=(
            "gbk\trecord_id\torder\trow\tcolumn\n"
            "{PROJECT_ROOT}/examples/MjeNMV.gb\tLC738868.1\t1\t1\t1\n"
            "{PROJECT_ROOT}/examples/MelaMJNV.gb\tLC738874.1\t2\t1\t2\n"
            "{PROJECT_ROOT}/examples/PemoMJNVA.gb\tLC738870.1\t3\t2\t1\n"
        ),
    ),
    "conservation.tsv": SupportAsset(
        filename="conservation.tsv",
        content=(
            "blast\tlabel\tcolor\n"
            "{PROJECT_ROOT}/examples/MjeNMV.MelaMJNV.tblastx.out\tMelaMJNV\t#4E79A7\n"
        ),
    ),
    "circular_tracks.tsv": SupportAsset(
        filename="circular_tracks.tsv",
        content=(
            "id\trenderer\tside\tr\tw\tparams\n"
            "features\tfeatures\taxis\n"
            "gc_content\tdinucleotide_content\tinside\t\t0.1\tnt=GC\n"
            "gc_skew\tdinucleotide_skew\tinside\t\t0.1\tnt=GC\n"
            "at_skew\tdinucleotide_skew\tinside\t\t0.1\tnt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew\n"
            "ticks\tticks\tinside\t\t\ttick_label_layout=label_in_tick_out\n"
        ),
    ),
    "feature_visibility.tsv": SupportAsset(
        filename="feature_visibility.tsv",
        content=(
            "record_id\tfeature_type\tqualifier\tvalue\taction\n"
            "NC_012920.1\tD-loop\tlocation\t^0\\.\\.16569$\tshow\n"
            "NC_012920.1\tCDS\tproduct\t^cytochrome c oxidase subunit I$\toff\n"
            "*\tCDS\tproduct\t^ATP synthase F0 subunit 6$\texclude_matching\n"
        ),
    ),
}


def load_palette_names(palettes_file: Path = PALETTES_FILE) -> tuple[str, ...]:
    with palettes_file.open("rb") as handle:
        data = tomllib.load(handle)
    return tuple(name for name in data.keys() if name != "title")


def get_support_assets() -> dict[str, SupportAsset]:
    return dict(SUPPORT_ASSETS)


def _file_arg(flag: str, *filenames: str) -> FileArgument:
    return FileArgument(flag=flag, filenames=tuple(filenames))


def _fasta_prep(output_filename: str, source_filename: str) -> FastaPreparation:
    return FastaPreparation(output_filename=output_filename, source_filename=source_filename)


def _blast_prep(
    output_filename: str,
    tool: BlastTool,
    query_filename: str,
    subject_filename: str,
) -> BlastPreparation:
    return BlastPreparation(
        output_filename=output_filename,
        tool=tool,
        query_filename=query_filename,
        subject_filename=subject_filename,
    )


def _figure(
    figure_id: str,
    output_path: str,
    groups: tuple[FigureGroup, ...],
    required_inputs: tuple[str, ...],
    recipe: Recipe,
    support_assets: tuple[str, ...] = (),
    preparations: tuple[Preparation, ...] = (),
    description: str = "",
) -> FigureSpec:
    return FigureSpec(
        figure_id=figure_id,
        output_path=output_path,
        groups=groups,
        required_inputs=required_inputs,
        recipe=recipe,
        support_assets=support_assets,
        preparations=preparations,
        description=description,
    )


def _docs_and_readme_figures() -> dict[str, FigureSpec]:
    figures: dict[str, FigureSpec] = {}

    figures["ecoli_k12_plot"] = _figure(
        figure_id="ecoli_k12_plot",
        output_path="examples/ecoli_k12_plot.svg",
        groups=("docs",),
        required_inputs=("NC_000913.gbk",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_000913.gbk",),
            extra_args=("--separate_strands",),
        ),
        description="Quickstart circular E. coli K-12 example.",
    )
    figures["ecoli_orchid"] = _figure(
        figure_id="ecoli_orchid",
        output_path="examples/ecoli_orchid.svg",
        groups=("docs",),
        required_inputs=("NC_000913.gbk",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_000913.gbk",),
            extra_args=("--separate_strands", "-p", "orchid"),
        ),
        description="Tutorial 1 orchid palette example.",
    )
    figures["ecoli_with_title"] = _figure(
        figure_id="ecoli_with_title",
        output_path="examples/ecoli_with_title.svg",
        groups=("docs",),
        required_inputs=("NC_000913.gbk",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_000913.gbk",),
            extra_args=(
                "--separate_strands",
                "--species",
                "<i>Escherichia coli</i>",
                "--strain",
                "K-12",
            ),
        ),
        description="Tutorial 1 centered title example.",
    )
    figures["MjeNMV_modified_default_colors"] = _figure(
        figure_id="MjeNMV_modified_default_colors",
        output_path="examples/MjeNMV_modified_default_colors.svg",
        groups=("docs",),
        required_inputs=("MjeNMV.gbk", "modified_default_colors.tsv"),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gbk",),
            file_args=(_file_arg("-d", "modified_default_colors.tsv"),),
            extra_args=("--separate_strands", "--track_type", "middle", "--block_stroke_width", "1"),
        ),
        description="Tutorial 3 default color override example.",
    )
    figures["MjeNMV_feature_specifc_colors_with_labels"] = _figure(
        figure_id="MjeNMV_feature_specifc_colors_with_labels",
        output_path="examples/MjeNMV_feature_specifc_colors_with_labels.svg",
        groups=("docs",),
        required_inputs=("MjeNMV.gbk", "modified_default_colors.tsv"),
        support_assets=("feature_specific_colors.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gbk",),
            file_args=(
                _file_arg("-d", "modified_default_colors.tsv"),
                _file_arg("-t", "feature_specific_colors.tsv"),
            ),
            extra_args=(
                "--separate_strands",
                "--track_type",
                "middle",
                "--block_stroke_width",
                "1",
                "--labels",
            ),
        ),
        description="Tutorial 3 feature-specific color table example.",
    )
    figures["MjeNMV_feature_specifc_colors_with_labels_blacklist"] = _figure(
        figure_id="MjeNMV_feature_specifc_colors_with_labels_blacklist",
        output_path="examples/MjeNMV_feature_specifc_colors_with_labels_blacklist.svg",
        groups=("docs",),
        required_inputs=("MjeNMV.gbk", "modified_default_colors.tsv"),
        support_assets=("feature_specific_colors.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gbk",),
            file_args=(
                _file_arg("-d", "modified_default_colors.tsv"),
                _file_arg("-t", "feature_specific_colors.tsv"),
            ),
            extra_args=(
                "--separate_strands",
                "--track_type",
                "middle",
                "--block_stroke_width",
                "1",
                "--labels",
                "--label_blacklist",
                "hypothetical",
            ),
        ),
        description="Tutorial 3 blacklist example.",
    )
    figures["O157_H7_stx_whitelist"] = _figure(
        figure_id="O157_H7_stx_whitelist",
        output_path="examples/O157_H7_stx_whitelist.svg",
        groups=("docs",),
        required_inputs=("O157_H7.gbk",),
        support_assets=("stx_whitelist.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("O157_H7.gbk",),
            file_args=(_file_arg("--label_whitelist", "stx_whitelist.tsv"),),
            extra_args=(
                "--labels",
                "--separate_strands",
                "--species",
                "<i>Escherichia coli</i> O157:H7",
                "--strain",
                "Sakai",
                "--label_font_size",
                "16",
            ),
        ),
        description="Tutorial 3 whitelist example.",
    )
    figures["HmmtDNA_qualifier_priority_soft_pastels"] = _figure(
        figure_id="HmmtDNA_qualifier_priority_soft_pastels",
        output_path="examples/HmmtDNA_qualifier_priority_soft_pastels.svg",
        groups=("docs",),
        required_inputs=("HmmtDNA.gbk",),
        support_assets=("qualifier_priority.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("HmmtDNA.gbk",),
            file_args=(_file_arg("--qualifier_priority", "qualifier_priority.tsv"),),
            extra_args=(
                "--track_type",
                "middle",
                "--species",
                "<i>Homo sapiens</i>",
                "--block_stroke_width",
                "2",
                "--axis_stroke_width",
                "5",
                "--labels",
                "both",
                "--palette",
                "soft_pastels",
                "--definition_font_size",
                "28",
                "--label_font_size",
                "18",
            ),
        ),
        description="Tutorial 3 qualifier priority example.",
    )
    figures["Escherichia_Shigella_pair"] = _figure(
        figure_id="Escherichia_Shigella_pair",
        output_path="examples/Escherichia_Shigella_pair.svg",
        groups=("docs",),
        required_inputs=("Escherichia_coli.gbk", "Shigella_dysenteriae.gbk"),
        preparations=(
            _fasta_prep("Escherichia_coli.fasta", "Escherichia_coli.gbk"),
            _fasta_prep("Shigella_dysenteriae.fasta", "Shigella_dysenteriae.gbk"),
            _blast_prep(
                "Escherichia_coli-Shigella_dysenteriae.blastn.out",
                "blastn",
                "Escherichia_coli.fasta",
                "Shigella_dysenteriae.fasta",
            ),
        ),
        recipe=CliRecipe(
            subcommand="linear",
            gbk_files=("Escherichia_coli.gbk", "Shigella_dysenteriae.gbk"),
            blast_files=("Escherichia_coli-Shigella_dysenteriae.blastn.out",),
            extra_args=("--align_center", "--separate_strands"),
        ),
        description="Tutorial 2 two-genome comparison example.",
    )
    figures["Escherichia_Shigella_multi"] = _figure(
        figure_id="Escherichia_Shigella_multi",
        output_path="examples/Escherichia_Shigella_multi.svg",
        groups=("docs",),
        required_inputs=(
            "Escherichia_coli.gbk",
            "Shigella_dysenteriae.gbk",
            "Shigella_flexneri.gbk",
            "Shigella_sonnei.gbk",
        ),
        preparations=(
            _fasta_prep("Escherichia_coli.fasta", "Escherichia_coli.gbk"),
            _fasta_prep("Shigella_dysenteriae.fasta", "Shigella_dysenteriae.gbk"),
            _fasta_prep("Shigella_flexneri.fasta", "Shigella_flexneri.gbk"),
            _fasta_prep("Shigella_sonnei.fasta", "Shigella_sonnei.gbk"),
            _blast_prep(
                "Escherichia_coli-Shigella_dysenteriae.blastn.out",
                "blastn",
                "Escherichia_coli.fasta",
                "Shigella_dysenteriae.fasta",
            ),
            _blast_prep(
                "Shigella_dysenteriae-Shigella_flexneri.blastn.out",
                "blastn",
                "Shigella_dysenteriae.fasta",
                "Shigella_flexneri.fasta",
            ),
            _blast_prep(
                "Shigella_flexneri-Shigella_sonnei.blastn.out",
                "blastn",
                "Shigella_flexneri.fasta",
                "Shigella_sonnei.fasta",
            ),
        ),
        recipe=CliRecipe(
            subcommand="linear",
            gbk_files=(
                "Escherichia_coli.gbk",
                "Shigella_dysenteriae.gbk",
                "Shigella_flexneri.gbk",
                "Shigella_sonnei.gbk",
            ),
            blast_files=(
                "Escherichia_coli-Shigella_dysenteriae.blastn.out",
                "Shigella_dysenteriae-Shigella_flexneri.blastn.out",
                "Shigella_flexneri-Shigella_sonnei.blastn.out",
            ),
            extra_args=("--align_center", "--separate_strands", "--evalue", "1e-99", "--bitscore", "5000"),
        ),
        description="Tutorial 2 multi-genome comparison example.",
    )
    figures["majani"] = _figure(
        figure_id="majani",
        output_path="examples/majani.svg",
        groups=("docs",),
        required_inputs=(
            "MjeNMV.gb",
            "MelaMJNV.gb",
            "PemoMJNVA.gb",
            "PeseMJNV.gb",
            "PemoMJNVB.gb",
            "LvMJNV.gb",
            "TrcuMJNV.gb",
            "MellatMJNV.gb",
            "MeenMJNV.gb",
            "MejoMJNV.gb",
            "custom_color_table.tsv",
            "modified_default_colors.tsv",
        ),
        preparations=(
            _fasta_prep("MjeNMV.fasta", "MjeNMV.gb"),
            _fasta_prep("MelaMJNV.fasta", "MelaMJNV.gb"),
            _fasta_prep("PemoMJNVA.fasta", "PemoMJNVA.gb"),
            _fasta_prep("PeseMJNV.fasta", "PeseMJNV.gb"),
            _fasta_prep("PemoMJNVB.fasta", "PemoMJNVB.gb"),
            _fasta_prep("LvMJNV.fasta", "LvMJNV.gb"),
            _fasta_prep("TrcuMJNV.fasta", "TrcuMJNV.gb"),
            _fasta_prep("MellatMJNV.fasta", "MellatMJNV.gb"),
            _fasta_prep("MeenMJNV.fasta", "MeenMJNV.gb"),
            _fasta_prep("MejoMJNV.fasta", "MejoMJNV.gb"),
            _blast_prep("MjeNMV.MelaMJNV.tblastx.out", "tblastx", "MjeNMV.fasta", "MelaMJNV.fasta"),
            _blast_prep("MelaMJNV.PemoMJNVA.tblastx.out", "tblastx", "MelaMJNV.fasta", "PemoMJNVA.fasta"),
            _blast_prep("PemoMJNVA.PeseMJNV.tblastx.out", "tblastx", "PemoMJNVA.fasta", "PeseMJNV.fasta"),
            _blast_prep("PeseMJNV.PemoMJNVB.tblastx.out", "tblastx", "PeseMJNV.fasta", "PemoMJNVB.fasta"),
            _blast_prep("PemoMJNVB.LvMJNV.tblastx.out", "tblastx", "PemoMJNVB.fasta", "LvMJNV.fasta"),
            _blast_prep("LvMJNV.TrcuMJNV.tblastx.out", "tblastx", "LvMJNV.fasta", "TrcuMJNV.fasta"),
            _blast_prep("TrcuMJNV.MellatMJNV.tblastx.out", "tblastx", "TrcuMJNV.fasta", "MellatMJNV.fasta"),
            _blast_prep("MellatMJNV.MeenMJNV.tblastx.out", "tblastx", "MellatMJNV.fasta", "MeenMJNV.fasta"),
            _blast_prep("MeenMJNV.MejoMJNV.tblastx.out", "tblastx", "MeenMJNV.fasta", "MejoMJNV.fasta"),
        ),
        recipe=CliRecipe(
            subcommand="linear",
            gbk_files=(
                "MjeNMV.gb",
                "MelaMJNV.gb",
                "PemoMJNVA.gb",
                "PeseMJNV.gb",
                "PemoMJNVB.gb",
                "LvMJNV.gb",
                "TrcuMJNV.gb",
                "MellatMJNV.gb",
                "MeenMJNV.gb",
                "MejoMJNV.gb",
            ),
            blast_files=(
                "MjeNMV.MelaMJNV.tblastx.out",
                "MelaMJNV.PemoMJNVA.tblastx.out",
                "PemoMJNVA.PeseMJNV.tblastx.out",
                "PeseMJNV.PemoMJNVB.tblastx.out",
                "PemoMJNVB.LvMJNV.tblastx.out",
                "LvMJNV.TrcuMJNV.tblastx.out",
                "TrcuMJNV.MellatMJNV.tblastx.out",
                "MellatMJNV.MeenMJNV.tblastx.out",
                "MeenMJNV.MejoMJNV.tblastx.out",
            ),
            file_args=(
                _file_arg("-t", "custom_color_table.tsv"),
                _file_arg("-d", "modified_default_colors.tsv"),
            ),
            extra_args=(
                "--block_stroke_width",
                "1",
                "--block_stroke_color",
                "gray",
                "--align_center",
                "--separate_strands",
            ),
        ),
        description="Gallery Majaniviruses comparison example.",
    )
    figures["NC_010162_edelweiss"] = _figure(
        figure_id="NC_010162_edelweiss",
        output_path="examples/NC_010162_edelweiss.svg",
        groups=("docs",),
        required_inputs=("NC_010162.gb", "NC_010162.feature-specific_table.tsv", "NC_010162.whitelist.tsv"),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_010162.gb",),
            file_args=(
                _file_arg("-t", "NC_010162.feature-specific_table.tsv"),
                _file_arg("--label_whitelist", "NC_010162.whitelist.tsv"),
            ),
            extra_args=(
                "--palette",
                "edelweiss",
                "--labels",
                "--separate_strands",
                "--species",
                "<i>Sorangium cellulosum</i>",
                "--strain",
                "So ce56",
            ),
        ),
        description="Gallery Sorangium cellulosum example.",
    )
    figures["NC_001879_color"] = _figure(
        figure_id="NC_001879_color",
        output_path="examples/NC_001879_color.svg",
        groups=("docs",),
        required_inputs=("NC_001879.gbk", "chloroplast_specific_table.tsv"),
        support_assets=("qualifier_priority.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_001879.gbk",),
            file_args=(
                _file_arg("-t", "chloroplast_specific_table.tsv"),
                _file_arg("--qualifier_priority", "qualifier_priority.tsv"),
            ),
            extra_args=(
                "--separate_strands",
                "-k",
                "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
                "--block_stroke_width",
                "1",
                "--block_stroke_color",
                "black",
                "--axis_stroke_width",
                "3",
                "--line_stroke_width",
                "2",
                "--no-gc",
                "--no-skew",
                "-p",
                "default",
                "--track_type",
                "tuckin",
                "--labels",
                "both",
                "--outer_label_x_radius_offset",
                "0.90",
                "--outer_label_y_radius_offset",
                "0.90",
                "--inner_label_x_radius_offset",
                "0.975",
                "--inner_label_y_radius_offset",
                "0.975",
                "--species",
                "<i>Nicotiana tabacum</i>",
                "--definition_font_size",
                "28",
                "--legend",
                "upper_left",
            ),
        ),
        description="Gallery chloroplast example.",
    )
    figures["NC_001879_regions"] = _figure(
        figure_id="NC_001879_regions",
        output_path="examples/NC_001879_regions.svg",
        groups=("docs",),
        required_inputs=(
            "NC_001879.gbk",
            "chloroplast_specific_table.tsv",
            "nicotiana-tabacum-regions.tsv",
        ),
        support_assets=("qualifier_priority.tsv",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("NC_001879.gbk",),
            file_args=(
                _file_arg("-t", "chloroplast_specific_table.tsv"),
                _file_arg("--qualifier_priority", "qualifier_priority.tsv"),
                _file_arg("--annotation_table", "nicotiana-tabacum-regions.tsv"),
            ),
            extra_args=(
                "--separate_strands",
                "-k",
                "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
                "--block_stroke_width",
                "1",
                "--block_stroke_color",
                "black",
                "--axis_stroke_width",
                "3",
                "--line_stroke_width",
                "2",
                "--no-skew",
                "-p",
                "default",
                "--track_type",
                "tuckin",
                "--labels",
                "both",
                "--label_placement",
                "radial",
                "--outer_label_x_radius_offset",
                "0.90",
                "--outer_label_y_radius_offset",
                "0.90",
                "--inner_label_x_radius_offset",
                "0.975",
                "--inner_label_y_radius_offset",
                "0.975",
                "--species",
                "<i>Nicotiana tabacum</i>",
                "--definition_font_size",
                "28",
                "--legend",
                "upper_left",
                "--circular_track_slot",
                "features:features@side=overlay,lane_direction=split",
                "--circular_track_slot",
                "plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,show_labels=true,padding_px=1,overflow=compress,inner_gap_px=1,outer_gap_px=1",
                "--circular_track_slot",
                "gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08",
            ),
        ),
        description="Tobacco chloroplast map with LSC, SSC, IRa, and IRb brackets.",
    )
    figures["M16-5_fugaku"] = _figure(
        figure_id="M16-5_fugaku",
        output_path="examples/M16-5_fugaku.svg",
        groups=("docs",),
        required_inputs=("M16-5.gb",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("M16-5.gb",),
            extra_args=(
                "--separate_strands",
                "--block_stroke_width",
                "1",
                "--axis_stroke_width",
                "1",
                "-p",
                "fugaku",
                "--track_type",
                "middle",
                "--species",
                "<i>Ca.</i> Sukunaarchaeum mirabile",
                "--definition_font_size",
                "22",
                "--legend",
                "upper_right",
            ),
        ),
        description="Gallery M16-5 example.",
    )
    figures["Pandoravirus_salinus_forest"] = _figure(
        figure_id="Pandoravirus_salinus_forest",
        output_path="examples/Pandoravirus_salinus_forest.svg",
        groups=("docs",),
        required_inputs=("Pandoravirus_salinus.gb",),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("Pandoravirus_salinus.gb",),
            extra_args=(
                "--separate_strands",
                "-p",
                "forest",
                "--track_type",
                "tuckin",
                "--species",
                "<i>Pandoravirus salinus</i>",
                "--definition_font_size",
                "22",
                "--legend",
                "upper_right",
            ),
        ),
        description="Gallery Pandoravirus example.",
    )

    figures["track_layout_separate_strands"] = _figure(
        figure_id="track_layout_separate_strands",
        output_path="examples/track_layout_separate_strands.png",
        groups=("docs", "composites"),
        required_inputs=("NC_000913.gbk",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=3,
            tile_size=(3000, 2800),
            gap=32,
            padding=32,
            panels=tuple(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_000913.gbk",),
                        extra_args=(("--separate_strands",) if separate_strands else ())
                        + (
                            "--track_type",
                            track_type,
                            "--plot_title",
                            f"{track_type} · {'separate strands' if separate_strands else 'combined strands'}",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "36",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                )
                for separate_strands in (True, False)
                for track_type in ("tuckin", "middle", "spreadout")
            ),
        ),
        description="Tutorial 1 track preset and strand separation comparison montage.",
    )
    figures["definition_font_size_comparison"] = _figure(
        figure_id="definition_font_size_comparison",
        output_path="examples/definition_font_size_comparison.png",
        groups=("docs", "composites"),
        required_inputs=("HmmtDNA.gbk",),
        support_assets=("qualifier_priority.tsv",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=3,
            tile_size=(3400, 3400),
            gap=32,
            padding=32,
            panels=tuple(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("HmmtDNA.gbk",),
                        file_args=(_file_arg("--qualifier_priority", "qualifier_priority.tsv"),),
                        extra_args=(
                            "--track_type",
                            "middle",
                            "--species",
                            "<i>Homo sapiens</i>",
                            "--block_stroke_width",
                            "2",
                            "--axis_stroke_width",
                            "5",
                            "--labels",
                            "both",
                            "--palette",
                            "soft_pastels",
                            "--label_font_size",
                            "18",
                            "--definition_font_size",
                            definition_font_size,
                            "--plot_title",
                            f"--definition_font_size {definition_font_size}",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "28",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                )
                for definition_font_size in ("20", "28", "36")
            ),
        ),
        description="Tutorial 3 definition font size comparison montage.",
    )
    figures["label_font_size_comparison"] = _figure(
        figure_id="label_font_size_comparison",
        output_path="examples/label_font_size_comparison.png",
        groups=("docs", "composites"),
        required_inputs=("HmmtDNA.gbk",),
        support_assets=("qualifier_priority.tsv",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=3,
            tile_size=(3400, 2200),
            gap=32,
            padding=32,
            panels=tuple(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("HmmtDNA.gbk",),
                        file_args=(_file_arg("--qualifier_priority", "qualifier_priority.tsv"),),
                        extra_args=(
                            "--track_type",
                            "middle",
                            "--species",
                            "<i>Homo sapiens</i>",
                            "--block_stroke_width",
                            "2",
                            "--line_stroke_color",
                            "gray",
                            "--line_stroke_width",
                            "1",
                            "--axis_stroke_width",
                            "5",
                            "--labels",
                            "both",
                            "--gc",
                            "--skew",
                            "--circular_track_order",
                            "features,ticks,gc_content,gc_skew",
                            "--palette",
                            "soft_pastels",
                            "--definition_font_size",
                            "28",
                            "--label_font_size",
                            label_font_size,
                            "--plot_title",
                            f"{label_font_size} pt labels",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "28",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                )
                for label_font_size in ("12", "18", "24")
            ),
        ),
        description="Tutorial 3 human mitochondrial label font size comparison montage.",
    )
    figures["outer_label_offset_comparison"] = _figure(
        figure_id="outer_label_offset_comparison",
        output_path="examples/outer_label_offset_comparison.png",
        groups=("docs", "composites"),
        required_inputs=("NC_001879.gbk",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=3,
            tile_size=(3400, 2900),
            gap=32,
            padding=32,
            panels=tuple(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_001879.gbk",),
                        extra_args=(
                            "--block_stroke_width",
                            "1",
                            "--track_type",
                            "middle",
                            "--labels",
                            "--outer_label_x_radius_offset",
                            x_offset,
                            "--outer_label_y_radius_offset",
                            y_offset,
                            "--plot_title",
                            f"x={x_offset} · y={y_offset}",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "36",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                )
                for y_offset in ("0.90", "1.00", "1.10")
                for x_offset in ("0.90", "1.00", "1.10")
            ),
        ),
        description="Tutorial 3 independent x/y outer label offset comparison matrix.",
    )
    figures["window_step_comparison"] = _figure(
        figure_id="window_step_comparison",
        output_path="examples/window_step_comparison.png",
        groups=("docs", "composites"),
        required_inputs=("NC_000913.gbk",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=3,
            tile_size=(3300, 3400),
            gap=32,
            padding=32,
            panels=(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_000913.gbk",),
                        extra_args=(
                            "--separate_strands",
                            "--window",
                            "100000",
                            "--step",
                            "10000",
                            "--plot_title",
                            "window=100000 · step=10000",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "36",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                ),
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_000913.gbk",),
                        extra_args=(
                            "--separate_strands",
                            "--window",
                            "10000",
                            "--step",
                            "1000",
                            "--plot_title",
                            "window=10000 · step=1000",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "36",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                ),
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_000913.gbk",),
                        extra_args=(
                            "--separate_strands",
                            "--window",
                            "1000",
                            "--step",
                            "100",
                            "--plot_title",
                            "window=1000 · step=100",
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "36",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                ),
            ),
        ),
        description="FAQ coarse-to-fine window and step comparison montage.",
    )
    figures["skew_comparison"] = _figure(
        figure_id="skew_comparison",
        output_path="examples/skew_comparison.png",
        groups=("docs", "composites"),
        required_inputs=("NC_000913.gbk",),
        recipe=CompositeRecipe(
            kind="grid",
            columns=4,
            tile_size=(2700, 2350),
            gap=32,
            padding=32,
            panels=tuple(
                CompositePanel(
                    recipe=CliRecipe(
                        subcommand="circular",
                        gbk_files=("NC_000913.gbk",),
                        extra_args=(
                            "--separate_strands",
                            "--nt",
                            dinucleotide,
                            "--plot_title",
                            dinucleotide,
                            "--plot_title_position",
                            "top",
                            "--plot_title_font_size",
                            "48",
                            "--keep_full_definition_with_plot_title",
                        ),
                    )
                )
                for dinucleotide in (
                    "GC",
                    "CG",
                    "AG",
                    "GA",
                    "CT",
                    "TC",
                    "TG",
                    "GT",
                    "CA",
                    "AC",
                    "AT",
                    "TA",
                )
            ),
        ),
        description="FAQ twelve-dinucleotide content and skew comparison montage.",
    )
    figures["python_api_circular"] = _figure(
        figure_id="python_api_circular",
        output_path="examples/python-api-circular.png",
        groups=("docs",),
        required_inputs=(
            "MjeNMV.gb",
            "custom_color_table.tsv",
            "modified_default_colors.tsv",
            "python-api-label-whitelist.tsv",
        ),
        recipe=CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gb",),
            file_args=(
                _file_arg("-t", "custom_color_table.tsv"),
                _file_arg("-d", "modified_default_colors.tsv"),
                _file_arg("--label_whitelist", "python-api-label-whitelist.tsv"),
            ),
            extra_args=(
                "-k",
                "CDS",
                "--species",
                "<i>Marsupenaeus japonicus endogenous nimavirus</i>",
                "--strain",
                "Ginoza2017",
                "--separate_strands",
                "--track_type",
                "middle",
                "--labels",
                "both",
                "--gc",
                "--skew",
                "--circular_track_order",
                "features,ticks,gc_content,gc_skew",
                "--block_stroke_width",
                "1",
                "--block_stroke_color",
                "gray",
                "--line_stroke_width",
                "2",
                "--line_stroke_color",
                "lightgray",
                "--legend",
                "right",
            ),
        ),
        description="CLI-equivalent MjeNMV rendering for the documented circular Python API example.",
    )
    figures["python_api_linear"] = _figure(
        figure_id="python_api_linear",
        output_path="examples/python-api-linear.png",
        groups=("docs",),
        required_inputs=(
            "MjeNMV.gb",
            "MelaMJNV.gb",
            "MjeNMV.MelaMJNV.tblastx.out",
            "custom_color_table.tsv",
            "modified_default_colors.tsv",
            "python-api-label-whitelist.tsv",
        ),
        recipe=CliRecipe(
            subcommand="linear",
            gbk_files=("MjeNMV.gb", "MelaMJNV.gb"),
            blast_files=("MjeNMV.MelaMJNV.tblastx.out",),
            file_args=(
                _file_arg("-t", "custom_color_table.tsv"),
                _file_arg("-d", "modified_default_colors.tsv"),
                _file_arg("--label_whitelist", "python-api-label-whitelist.tsv"),
            ),
            extra_args=(
                "-k",
                "CDS",
                "--evalue",
                "1e-5",
                "--bitscore",
                "0",
                "--identity",
                "0",
                "--gc",
                "--skew",
                "--separate_strands",
                "--show_labels",
                "first",
                "--align_center",
                "--plot_title",
                "Majanivirus genome comparison",
                "--plot_title_position",
                "top",
                "--block_stroke_width",
                "1",
                "--block_stroke_color",
                "gray",
                "--line_stroke_width",
                "2",
                "--line_stroke_color",
                "lightgray",
                "--legend",
                "right",
            ),
        ),
        description="CLI-equivalent Majanivirus comparison for the documented linear Python API example.",
    )

    return figures


def _remaining_tutorial_figures() -> dict[str, FigureSpec]:
    """Recipes for tutorial figures that are not part of the legacy gallery set."""

    figures: dict[str, FigureSpec] = {}

    def add(
        figure_id: str,
        output_path: str,
        recipe: Recipe,
        *,
        required_inputs: tuple[str, ...],
        support_assets: tuple[str, ...] = (),
    ) -> None:
        figures[figure_id] = _figure(
            figure_id=figure_id,
            output_path=output_path,
            groups=("docs",),
            required_inputs=required_inputs,
            support_assets=support_assets,
            recipe=recipe,
            description=f"Documentation tutorial figure: {figure_id}.",
        )

    linear_feature_style_args = (
        "--block_stroke_color",
        "gray",
        "--block_stroke_width",
        "1",
        "--line_stroke_color",
        "lightgray",
        "--line_stroke_width",
        "2",
    )

    add(
        "quickstart_labeled_rna_features",
        "examples/quickstart-labeled-rna-features.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("NC_000913.gbk",),
            extra_args=("--track_type", "middle", "--features", "rRNA,tRNA", "--labels"),
        ),
        required_inputs=("NC_000913.gbk",),
    )
    add(
        "tutorial_2_pairwise_blast",
        "examples/tutorial-2-pairwise-blast.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("AP027078.gb", "AP027131.gb"),
            blast_files=("AP027078_AP027131.tblastx.out",),
            extra_args=linear_feature_style_args
            + (
                "--separate_strands",
                "--align_center",
                "--gc",
                "--skew",
                "--palette",
                "ajisai",
                "--scale_style",
                "ruler",
                "--pairwise_match_style",
                "curve",
                "--record_label",
                "Candidatus Tyloplasma litorale",
                "--record_label",
                "Candidatus Hepatoplasma vulgare",
                "--record_subtitle",
                "Fukuoka2020",
                "--record_subtitle",
                "Av-JP",
                "--plot_title",
                "Hepatoplasmataceae pairwise TBLASTX",
            ),
        ),
        required_inputs=("AP027078.gb", "AP027131.gb", "AP027078_AP027131.tblastx.out"),
    )
    add(
        "linear_multi_record",
        "examples/linear_multi_record.svg",
        CliRecipe(
            subcommand="linear",
            file_args=(
                _file_arg("--records_table", "linear_multi_records.tsv"),
                _file_arg("--comparisons_table", "linear_multi_comparisons.tsv"),
            ),
            extra_args=linear_feature_style_args
            + (
                "--linear_record_gap", "28",
                "--scale_style", "ruler",
                "--pairwise_match_style", "curve",
                "--identity", "97",
                "--alignment_length", "500",
                "--show_labels", "none",
            ),
        ),
        required_inputs=(
            "MjeNMV.gb", "PemoMJNVA.gb", "MelaMJNV.gb", "PeseMJNV.gb",
            "MjeNMV.MelaMJNV.tblastx.out", "PemoMJNVA.PeseMJNV.tblastx.out",
        ),
        support_assets=("linear_multi_records.tsv", "linear_multi_comparisons.tsv"),
    )
    add(
        "vibrio_nigripulchritudo_multi_record",
        "examples/vibrio-nigripulchritudo-multi-record.svg",
        CliRecipe(
            subcommand="linear",
            file_args=(
                _file_arg("--records_table", "vibrio-nigripulchritudo-linear-records.tsv"),
            ),
            extra_args=(
                "--linear_record_gap", "48",
                "--track_layout", "above",
                "--scale_style", "ruler",
                "--ruler_on_axis",
                "--scale_interval", "750000",
                "--separate_strands",
                "--gc",
                "--hide_accession",
                "--hide_length",
                "--definition_font_size", "8",
                "--keep_definition_left_aligned",
                "--protein_blastp_mode", "collinear",
                "--collinear_search_scope", "all",
                "--protein_blastp_candidate_limit", "5",
                "--collinear_min_anchors", "3",
                "--collinear_max_unit_gap", "2",
                "--collinear_max_diagonal_drift", "2",
                "--collinear_color_mode", "orientation_identity",
                "--pairwise_match_style", "curve",
                "--losatp_bin", "gbdraw/bin/linux-x86_64/losat",
                "--losatp_threads", "8",
                "--plot_title", "Vibrio nigripulchritudo replicons: LOSATP collinear blocks",
                "--plot_title_position", "top",
            ),
        ),
        required_inputs=(
            "GCF_015097735.1_ASM1509773v1_genomic.gbff",
            "GCF_000801275.2_ASM80127v1_genomic.gbff",
        ),
    )
    add(
        "tutorial_3_label_override",
        "examples/tutorial-3-label-override.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gbk",),
            file_args=(_file_arg("--label_table", "label_override.tsv"),),
            extra_args=("--labels",),
        ),
        required_inputs=("MjeNMV.gbk",),
        support_assets=("label_override.tsv",),
    )
    add(
        "tutorial_3_embedded_labels",
        "examples/tutorial-3-embedded-labels.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gbk",),
            extra_args=("--labels", "--label_rendering", "embedded_only"),
        ),
        required_inputs=("MjeNMV.gbk",),
    )
    add(
        "tutorial_3_above_feature_labels",
        "examples/tutorial-3-above-feature-labels.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("MjeNMV.gbk",),
            extra_args=(
                "--show_labels",
                "all",
                "--label_placement",
                "above_feature",
                "--label_rotation",
                "45",
            ),
        ),
        required_inputs=("MjeNMV.gbk",),
    )

    protein_inputs = ("MjeNMV.gb", "MelaMJNV.gb", "PemoMJNVA.gb")
    add(
        "tutorial_protein_pairwise",
        "examples/tutorial-protein-pairwise.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=protein_inputs[:2],
            file_args=(
                _file_arg("-t", "majani_custom_color_table.tsv"),
                _file_arg("-d", "modified_default_colors.tsv"),
            ),
            extra_args=linear_feature_style_args
            + (
                "--protein_blastp_mode",
                "pairwise",
                "--align_center",
                "--pairwise_match_style",
                "curve",
                "--gc",
                "--skew",
                "--show_labels",
                "first",
                "--label_rendering",
                "embedded_only",
                "--record_label",
                "Marsupenaeus japonicus endogenous nimavirus",
                "--record_label",
                "Melicertus latisulcatus majanivirus",
                "--record_subtitle",
                "Ginoza2017",
                "--record_subtitle",
                "Okinawa2016",
                "--scale_style",
                "ruler",
                "--plot_title",
                "Majanivirus pairwise protein matches",
            ),
        ),
        required_inputs=protein_inputs[:2]
        + ("majani_custom_color_table.tsv", "modified_default_colors.tsv"),
    )
    add(
        "majani_orthogroup",
        "examples/majani_orthogroup.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=protein_inputs,
            extra_args=(
                "--protein_blastp_mode",
                "orthogroup",
                "--show_labels",
                "orthogroup_top",
                "--pairwise_match_style",
                "curve",
                "--align_center",
            ),
        ),
        required_inputs=protein_inputs,
    )
    add(
        "majani_collinear",
        "examples/majani_collinear.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=protein_inputs,
            extra_args=(
                "--protein_blastp_mode",
                "collinear",
                "--collinear_min_anchors",
                "2",
                "--collinear_color_mode",
                "orientation_identity",
                "--pairwise_match_style",
                "curve",
                "--align_center",
            ),
        ),
        required_inputs=protein_inputs,
    )

    add(
        "tutorial_5_records_table",
        "examples/tutorial-5-records-table.svg",
        CliRecipe(
            subcommand="linear",
            file_args=(
                _file_arg("--records_table", "linear_records.tsv"),
                _file_arg("-t", "majani_custom_color_table.tsv"),
                _file_arg("-d", "modified_default_colors.tsv"),
            ),
            extra_args=linear_feature_style_args
            + (
                "--gc",
                "--skew",
                "--show_labels",
                "first",
                "--label_rendering",
                "embedded_only",
                "--scale_style",
                "ruler",
                "--plot_title",
                "Majanivirus records table",
            ),
        ),
        required_inputs=protein_inputs
        + ("majani_custom_color_table.tsv", "modified_default_colors.tsv"),
        support_assets=("linear_records.tsv",),
    )
    add(
        "tutorial_5_circular_grid",
        "examples/tutorial-5-circular-grid.svg",
        CliRecipe(
            subcommand="circular",
            file_args=(_file_arg("--records_table", "circular_records.tsv"),),
            extra_args=("--multi_record_canvas", "--multi_record_size_mode", "auto"),
        ),
        required_inputs=protein_inputs,
        support_assets=("circular_records.tsv",),
    )
    add(
        "tutorial_5_conservation_table",
        "examples/tutorial-5-conservation-table.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("MjeNMV.gb",),
            file_args=(_file_arg("--conservation_table", "conservation.tsv"),),
            extra_args=(
                "--conservation_reference",
                "query",
                "--identity",
                "95",
                "--alignment_length",
                "750",
                "--conservation_ring_width",
                "12",
                "--conservation_ring_gap",
                "3",
            ),
        ),
        required_inputs=("MjeNMV.gb", "MjeNMV.MelaMJNV.tblastx.out"),
        support_assets=("conservation.tsv",),
    )
    add(
        "tutorial_circular_track_table",
        "examples/tutorial-circular-track-table.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("HmmtDNA.gbk",),
            file_args=(_file_arg("--circular_track_table", "circular_tracks.tsv"),),
            extra_args=(
                "--track_type",
                "middle",
                "--window",
                "500",
                "--step",
                "50",
                "--species",
                "<i>Homo sapiens</i>",
                "-l",
                "left",
            ),
        ),
        required_inputs=("HmmtDNA.gbk",),
        support_assets=("circular_tracks.tsv",),
    )

    add(
        "tutorial_6_depth_circular",
        "examples/tutorial-6-depth-circular.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("AP027133.gb",),
            file_args=(_file_arg("--depth_track", "AP027133.DRR394922.depth.tsv"),),
            extra_args=(
                "--depth_width",
                "45",
                "--depth_window",
                "100",
                "--depth_step",
                "100",
                "--depth_max",
                "150",
                "--show_depth_axis",
                "--show_depth_ticks",
                "--depth_large_tick_interval",
                "50",
                "--depth_small_tick_interval",
                "25",
            ),
        ),
        required_inputs=("AP027133.gb", "AP027133.DRR394922.depth.tsv"),
    )
    add(
        "tutorial_6_depth_tracks",
        "examples/tutorial-6-depth-tracks.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("AP027131.gb", "AP027132.gb"),
            file_args=(
                _file_arg(
                    "--depth_track",
                    "AP027131.DRR394921.depth.tsv",
                    "AP027132.DRR394921.depth.tsv",
                ),
            ),
            extra_args=(
                "--depth_track_label",
                "DRR394921",
                "--depth_track_color",
                "#4E79A7",
                "--depth_height",
                "36",
                "--depth_window",
                "100",
                "--depth_step",
                "100",
                "--depth_min",
                "0",
                "--depth_max",
                "1000",
                "--depth_large_tick_interval",
                "500",
                "--depth_small_tick_interval",
                "250",
                "--show_depth_axis",
                "--show_depth_ticks",
                "--share_depth_axis",
            ),
        ),
        required_inputs=(
            "AP027131.gb",
            "AP027132.gb",
            "AP027131.DRR394921.depth.tsv",
            "AP027132.DRR394921.depth.tsv",
        ),
    )
    add(
        "tutorial_depth_log_axis",
        "examples/tutorial-depth-log-axis.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("AP027078.gb",),
            file_args=(_file_arg("--depth_track", "AP027078.DRR394944.depth.tsv"),),
            extra_args=(
                "--depth_height",
                "40",
                "--depth_window",
                "100",
                "--depth_step",
                "100",
                "--depth_log_scale",
                "--depth_min",
                "1",
                "--depth_max",
                "250000",
                "--show_depth_axis",
                "--show_depth_ticks",
            ),
        ),
        required_inputs=("AP027078.gb", "AP027078.DRR394944.depth.tsv"),
    )
    add(
        "tutorial_6_gc_percent",
        "examples/tutorial-6-gc-percent.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("MjeNMV.gb",),
            extra_args=(
                "--gc",
                "--gc_content_mode",
                "percent",
                "--gc_content_min_percent",
                "25",
                "--gc_content_max_percent",
                "75",
                "--gc_content_large_tick_interval",
                "10",
                "--gc_content_small_tick_interval",
                "5",
                "--show_gc_content_axis",
                "--show_gc_content_ticks",
            ),
        ),
        required_inputs=("MjeNMV.gb",),
    )
    add(
        "tutorial_6_two_skew_tracks",
        "examples/tutorial-6-two-skew-tracks.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("MjeNMV.gb",),
            extra_args=(
                "--skew",
                "--linear_track_slot",
                "features:features@side=overlay",
                "--linear_track_slot",
                "gc_skew:gc_skew@side=below,h=24px,spacing=8px",
                "--linear_track_slot",
                "at_skew:dinucleotide_skew@side=below,h=24px,spacing=8px,nt=AT,positive_color=#deaf6e,negative_color=#7294e3",
                "--linear_track_axis_index",
                "0",
            ),
        ),
        required_inputs=("MjeNMV.gb",),
    )

    layout_inputs = ("MjeNMV.gb", "MelaMJNV.gb")
    layout_blast_inputs = ("MjeNMV.MelaMJNV.tblastx.out",)
    layout_file_args = (
        _file_arg("-t", "majani_custom_color_table.tsv"),
        _file_arg("-d", "modified_default_colors.tsv"),
    )
    layout_required_inputs = (
        layout_inputs
        + layout_blast_inputs
        + (
            "majani_custom_color_table.tsv",
            "modified_default_colors.tsv",
        )
    )
    layout_visual_args = linear_feature_style_args + (
        "--gc",
        "--skew",
        "--pairwise_match_style",
        "curve",
        "--show_labels",
        "first",
        "--label_rendering",
        "embedded_only",
    )
    layout_record_text_args = (
        "--record_label",
        "Marsupenaeus japonicus endogenous nimavirus",
        "--record_label",
        "Melicertus latisulcatus majanivirus",
        "--record_subtitle",
        "Ginoza2017",
        "--record_subtitle",
        "Okinawa2016",
    )
    add(
        "tutorial_7_track_layout_below",
        "examples/tutorial-7-track-layout-below.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=layout_inputs,
            blast_files=layout_blast_inputs,
            file_args=layout_file_args,
            extra_args=(
                layout_visual_args
                + layout_record_text_args
                + (
                    "--track_layout",
                    "below",
                    "--track_axis_gap",
                    "auto",
                )
            ),
        ),
        required_inputs=layout_required_inputs,
    )
    add(
        "tutorial_7_linear_layout",
        "examples/tutorial-7-linear-layout.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=layout_inputs,
            blast_files=layout_blast_inputs,
            file_args=layout_file_args,
            extra_args=layout_visual_args
            + layout_record_text_args
            + (
                "--scale_style",
                "ruler",
                "--scale_interval",
                "50000",
                "--plot_title",
                "Majanivirus comparison",
                "--plot_title_position",
                "top",
            ),
        ),
        required_inputs=layout_required_inputs,
    )
    add(
        "tutorial_7_definition_lines",
        "examples/tutorial-7-definition-lines.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=layout_inputs,
            blast_files=layout_blast_inputs,
            file_args=layout_file_args,
            extra_args=layout_visual_args
            + (
                "--record_label",
                "MjeNMV",
                "--record_label",
                "MelaMJNV",
                "--record_subtitle",
                "Ginoza2017",
                "--record_subtitle",
                "Okinawa2016",
                "--align_center",
                "--hide_accession",
                "--hide_length",
                "--keep_definition_left_aligned",
                "--definition_line_style",
                "name:weight=bold,size=18",
                "--definition_line_style",
                "subtitle:size=14,color=#555555",
            ),
        ),
        required_inputs=layout_required_inputs,
    )
    add(
        "tutorial_7_linear_track_slots",
        "examples/tutorial-7-linear-track-slots.svg",
        CliRecipe(
            subcommand="linear",
            gbk_files=("MjeNMV.gb",),
            extra_args=(
                "--gc",
                "--skew",
                "--linear_track_slot",
                "features:features@side=overlay,h=60px",
                "--linear_track_slot",
                "gc_content:gc_content@h=24px,spacing=8px",
                "--linear_track_slot",
                "gc_skew:gc_skew@h=24px,spacing=8px",
                "--linear_track_axis_index",
                "0",
            ),
        ),
        required_inputs=("MjeNMV.gb",),
    )

    add(
        "tutorial_9_arrow_geometry_circular",
        "examples/tutorial-9-arrow-geometry-circular.svg",
        SessionVariantRecipe(
            source_session_path=(
                "gbdraw/web/gallery/sessions/"
                "HmmtDNA_ATskew.gbdraw-session.json"
            ),
            config_overrides=(
                (
                    "objects.features.arrow_geometry.head_length_ratio",
                    "auto",
                ),
                (
                    "objects.features.arrow_geometry.shaft_width_ratio",
                    0.75,
                ),
            ),
        ),
        required_inputs=(
            "gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json",
        ),
    )
    add(
        "tutorial_9_arrow_geometry_linear",
        "examples/tutorial-9-arrow-geometry-linear.svg",
        SessionVariantRecipe(
            source_session_path=(
                "gbdraw/web/gallery/sessions/"
                "BGC0000708-BGC0000713.gbdraw-session.json"
            ),
            config_overrides=(
                (
                    "objects.features.arrow_geometry.head_length_ratio",
                    "auto",
                ),
                (
                    "objects.features.arrow_geometry.shaft_width_ratio",
                    0.5,
                ),
            ),
        ),
        required_inputs=(
            "gbdraw/web/gallery/sessions/"
            "BGC0000708-BGC0000713.gbdraw-session.json",
        ),
    )
    add(
        "tutorial_9_feature_shapes",
        "examples/tutorial-9-feature-shapes.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("HmmtDNA.gbk",),
            extra_args=(
                "-k",
                "CDS,rRNA,tRNA",
                "--feature_shape",
                "CDS=rectangle",
                "--feature_shape",
                "rRNA=rectangle",
                "--feature_shape",
                "tRNA=rectangle",
                "--labels",
                "out",
                "--track_type",
                "middle",
            ),
        ),
        required_inputs=("HmmtDNA.gbk",),
    )
    add(
        "tutorial_9_feature_visibility",
        "examples/tutorial-9-feature-visibility.svg",
        CliRecipe(
            subcommand="circular",
            gbk_files=("HmmtDNA.gbk",),
            file_args=(_file_arg("--feature_visibility_table", "feature_visibility.tsv"),),
            extra_args=(
                "-k",
                "CDS,rRNA,tRNA",
                "--feature_shape",
                "CDS=rectangle",
                "--feature_shape",
                "rRNA=rectangle",
                "--feature_shape",
                "tRNA=rectangle",
                "--feature_shape",
                "D-loop=rectangle",
                "--labels",
                "out",
                "--track_type",
                "middle",
            ),
        ),
        required_inputs=("HmmtDNA.gbk",),
        support_assets=("feature_visibility.tsv",),
    )
    return figures


def _palette_figures() -> dict[str, FigureSpec]:
    figures: dict[str, FigureSpec] = {}
    representative_names = ("default", "ajisai", "soft_pastels")

    palette_base_inputs = (
        "AP027078.gb",
        "AP027131.gb",
        "AP027133.gb",
        "AP027132.gb",
        "NZ_CP006932.gb",
    )
    palette_preparations = (
        _fasta_prep("AP027078.fasta", "AP027078.gb"),
        _fasta_prep("AP027131.fasta", "AP027131.gb"),
        _fasta_prep("AP027133.fasta", "AP027133.gb"),
        _fasta_prep("AP027132.fasta", "AP027132.gb"),
        _fasta_prep("NZ_CP006932.fasta", "NZ_CP006932.gb"),
        _blast_prep("AP027078_AP027131.tblastx.out", "tblastx", "AP027078.fasta", "AP027131.fasta"),
        _blast_prep("AP027131_AP027133.tblastx.out", "tblastx", "AP027131.fasta", "AP027133.fasta"),
        _blast_prep("AP027133_AP027132.tblastx.out", "tblastx", "AP027133.fasta", "AP027132.fasta"),
        _blast_prep("AP027132_NZ_CP006932.tblastx.out", "tblastx", "AP027132.fasta", "NZ_CP006932.fasta"),
    )

    def circular_recipe(palette_name: str) -> CliRecipe:
        return CliRecipe(
            subcommand="circular",
            gbk_files=("AP027078.gb",),
            extra_args=("--separate_strands", "--track_type", "tuckin", "-p", palette_name),
        )

    def linear_recipe(palette_name: str) -> CliRecipe:
        return CliRecipe(
            subcommand="linear",
            gbk_files=palette_base_inputs,
            blast_files=(
                "AP027078_AP027131.tblastx.out",
                "AP027131_AP027133.tblastx.out",
                "AP027133_AP027132.tblastx.out",
                "AP027132_NZ_CP006932.tblastx.out",
            ),
            extra_args=("--align_center", "--separate_strands", "--gc", "--skew", "-p", palette_name),
        )

    for palette_name in representative_names:
        circular_id = f"palette_circular_{palette_name}"
        linear_id = f"palette_linear_{palette_name}"
        figures[circular_id] = _figure(
            figure_id=circular_id,
            output_path=f"examples/AP027078_tuckin_separate_strands_{palette_name}.svg",
            groups=("palettes",),
            required_inputs=("AP027078.gb",),
            recipe=circular_recipe(palette_name),
            description=f"Circular palette example for {palette_name}.",
        )
        figures[linear_id] = _figure(
            figure_id=linear_id,
            output_path=f"examples/hepatoplasmataceae_{palette_name}.svg",
            groups=("palettes",),
            required_inputs=palette_base_inputs,
            preparations=palette_preparations,
            recipe=linear_recipe(palette_name),
            description=f"Linear palette example for {palette_name}.",
        )

    return figures


def build_figure_specs() -> dict[str, FigureSpec]:
    figures = _docs_and_readme_figures()
    figures.update(_remaining_tutorial_figures())
    figures.update(_palette_figures())
    return figures


__all__ = [
    "BlastPreparation",
    "CliRecipe",
    "CompositePanel",
    "CompositeRecipe",
    "FigureSpec",
    "FastaPreparation",
    "FileArgument",
    "MANUALLY_MANAGED_FIGURES",
    "PALETTES_FILE",
    "PROJECT_ROOT",
    "SUPPORT_ASSETS",
    "SupportAsset",
    "UNREFERENCED_FIGURE_RETENTION",
    "build_figure_specs",
    "get_support_assets",
    "load_palette_names",
]
