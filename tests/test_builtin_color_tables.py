"""The built-in color tables (``gbdraw/data/color_tables/``) and the ``builtin:<name>`` selector."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import build_specific_color_rules
from gbdraw.features.selector_values import find_specific_color_rule
from gbdraw.io.colors import (
    BUILTIN_COLOR_TABLE_PREFIX,
    builtin_color_table_data,
    builtin_color_table_names,
    color_table_path,
    read_color_table,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPO_ROOT / "gbdraw" / "data" / "color_tables"
TOBACCO = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "tobacco-plastome-regions" / "NC_001879.gbk"
HUMAN_MT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"


def _caption(table: str, feature: SeqFeature) -> str | None:
    rules = build_specific_color_rules(read_color_table(f"{BUILTIN_COLOR_TABLE_PREFIX}{table}"))
    matched = find_specific_color_rule(feature, rules, "record", instance=None)
    return matched[1] if matched is not None else None


def _cds(**qualifiers: list[str]) -> SeqFeature:
    return SeqFeature(FeatureLocation(0, 300, strand=1), type="CDS", qualifiers=dict(qualifiers))


def test_the_builtin_names_are_the_packaged_tables() -> None:
    assert builtin_color_table_names() == ("cog", "phrog", "plastid", "mitochondrion")
    assert sorted(path.stem for path in TABLE_DIR.glob("*.tsv")) == sorted(builtin_color_table_names())


def test_a_builtin_name_resolves_to_its_packaged_file_and_a_path_stays_a_path(tmp_path: Path) -> None:
    resolved = Path(color_table_path("builtin:cog"))
    assert resolved.is_file()
    assert resolved.read_bytes() == (TABLE_DIR / "cog.tsv").read_bytes()
    own = tmp_path / "builtin-like.tsv"
    assert color_table_path(str(own)) == str(own)


def test_an_unknown_builtin_name_is_refused_with_the_choices() -> None:
    with pytest.raises(ValidationError, match=r"builtin:cog, builtin:phrog, builtin:plastid, builtin:mitochondrion"):
        color_table_path("builtin:kegg")


@pytest.mark.parametrize("name", ["cog", "phrog", "plastid", "mitochondrion"])
def test_every_builtin_table_reads_and_its_colors_are_lowercase_hex(name: str) -> None:
    table = read_color_table(f"builtin:{name}")
    assert table is not None and not table.empty
    assert all(re.fullmatch(r"#[0-9a-f]{6}", color) for color in table["color"])
    assert all(table["caption"])
    # One caption has one color, so the Legend shows one row per caption.
    assert (table.groupby("caption")["color"].nunique() == 1).all()


def test_the_web_data_lists_each_table_with_its_rows_as_python_reads_them() -> None:
    data = builtin_color_table_data()
    assert list(data) == list(builtin_color_table_names())
    cog = read_color_table("builtin:cog")
    assert data["cog"]["rows"] == cog[["feature_type", "qualifier_key", "value", "color", "caption"]].values.tolist()
    assert all(isinstance(entry["title"], str) and entry["title"] for entry in data.values())


@pytest.mark.parametrize(
    ("qualifiers", "caption"),
    [
        # Bakta GenBank: the cluster ID and the category letters are both db_xref values.
        ({"db_xref": ["COG:COG0593", "COG:L", "UniRef:UniRef90_P03004"]}, "Replication, recombination and repair"),
        # Two letters: the first one, as Bakta's own plot reads them.
        ({"db_xref": ["COG:COG0001", "COG:KL"]}, "Transcription"),
        # Bakta --compliant moves non-INSDC cross-references to /note.
        ({"note": ["COG:COG0001", "COG:J", "hypothetical protein"]}, "Translation, ribosomal structure and biogenesis"),
        # Bakta GFF3 writes Dbxref.
        ({"Dbxref": ["COG:COG0001", "COG:S"]}, "Function unknown"),
        # A cluster ID alone names no category.
        ({"db_xref": ["COG:COG0593"]}, None),
    ],
)
def test_the_cog_table_reads_bakta_category_letters(qualifiers: dict[str, list[str]], caption: str | None) -> None:
    assert _caption("cog", _cds(**qualifiers)) == caption


@pytest.mark.parametrize(
    ("qualifiers", "caption"),
    [
        ({"function": ["head and packaging"]}, "Head and packaging"),
        ({"function": ["DNA, RNA and nucleotide metabolism"]}, "DNA, RNA and nucleotide metabolism"),
        ({"function": ["moron, auxiliary metabolic gene and host takeover"]},
         "Moron, auxiliary metabolic gene and host takeover"),
        ({"function": ["unknown function"]}, "Unknown function"),
        ({"function": ["other"]}, "Other function"),
        ({"vfdb_short_name": ["stx2A"]}, "Virulence factor/AMR"),
    ],
)
def test_the_phrog_table_reads_pharokka_functions(qualifiers: dict[str, list[str]], caption: str) -> None:
    assert _caption("phrog", _cds(**qualifiers)) == caption


def _classes(path: Path, table: str) -> dict[tuple[str, str], str | None]:
    record = SeqIO.read(path, "genbank")
    rules = build_specific_color_rules(read_color_table(f"builtin:{table}"))
    found: dict[tuple[str, str], str | None] = {}
    for feature in record.features:
        if feature.type not in {"CDS", "tRNA", "rRNA", "rep_origin"}:
            continue
        name = (feature.qualifiers.get("gene") or feature.qualifiers.get("product")
                or feature.qualifiers.get("note") or [""])[0]
        matched = find_specific_color_rule(feature, rules, record.id, instance=None)
        found[(feature.type, name)] = matched[1] if matched is not None else None
    return found


def test_the_plastid_table_classes_every_tobacco_gene() -> None:
    classes = _classes(TOBACCO, "plastid")
    assert [key for key, caption in classes.items() if caption is None] == []
    assert classes[("CDS", "psaA")] == "photosystem I"
    assert classes[("CDS", "psbN")] == "photosystem II"
    assert classes[("CDS", "petN")] == "cytochrome b/f complex"
    assert classes[("CDS", "rbcL")] == "RubisCO large subunit"
    assert classes[("CDS", "rpoC2")] == "RNA polymerase"
    assert classes[("CDS", "rps12")] == "ribosomal proteins (SSU)"
    assert classes[("CDS", "rpl2")] == "ribosomal proteins (LSU)"
    assert classes[("CDS", "matK")] == "clpP, matK"
    assert classes[("CDS", "accD")] == "other genes"
    assert classes[("CDS", "ycf2")] == "hypothetical chloroplast reading frames (ycf)"
    assert classes[("tRNA", "trnfM")] == "transfer RNAs"
    assert classes[("rRNA", "16S ribosomal RNA")] == "ribosomal RNAs"
    assert classes[("rep_origin", "oriA")] == "origin of replication"


def test_the_mitochondrion_table_classes_every_human_gene() -> None:
    classes = _classes(HUMAN_MT, "mitochondrion")
    assert [key for key, caption in classes.items() if caption is None] == []
    assert classes[("CDS", "ND4L")] == "complex I (NADH dehydrogenase)"
    assert classes[("CDS", "CYTB")] == "complex III (ubiquinol-cytochrome c reductase)"
    assert classes[("CDS", "COX2")] == "complex IV (cytochrome c oxidase)"
    assert classes[("CDS", "ATP8")] == "ATP synthase"
    assert classes[("tRNA", "TRNL1")] == "transfer RNAs"
    assert classes[("rRNA", "RNR2")] == "ribosomal RNAs"


@pytest.mark.parametrize(
    ("gene", "caption"),
    [
        ("nad4L", "complex I (NADH dehydrogenase)"),
        ("mt-Nd5", "complex I (NADH dehydrogenase)"),
        ("sdh4", "complex II (succinate dehydrogenase)"),
        ("cob", "complex III (ubiquinol-cytochrome c reductase)"),
        ("cox1", "complex IV (cytochrome c oxidase)"),
        ("atp9", "ATP synthase"),
        ("ccmFn", "cytochrome c biogenesis"),
        ("rps3", "ribosomal proteins (SSU)"),
        ("rpl5", "ribosomal proteins (LSU)"),
        ("matR", "maturases"),
        ("mttB", "other genes"),
        ("orf25", "ORFs"),
    ],
)
def test_the_mitochondrion_table_reads_plant_and_animal_gene_names(gene: str, caption: str) -> None:
    assert _caption("mitochondrion", _cds(gene=[gene])) == caption


def test_the_python_api_and_the_cli_take_a_builtin_table(tmp_path: Path, gbdraw_runner) -> None:
    record = SeqIO.read(HUMAN_MT, "genbank")
    svg = gbdraw.draw_circular(
        [record], options=gbdraw.CircularOptions(features=gbdraw.FeatureOptions(color_table="builtin:mitochondrion"))
    ).to_svg()
    assert "complex IV (cytochrome c oxidase)" in svg

    code, output, svg_path = gbdraw_runner.run(
        "circular", [HUMAN_MT], "mt", tmp_path, extra_args=["-t", "builtin:mitochondrion"]
    )
    assert code == 0, output
    assert "complex IV (cytochrome c oxidase)" in svg_path.read_text(encoding="utf-8")


def test_a_request_with_a_builtin_table_saves_the_table_rows() -> None:
    from gbdraw.api.options import ColorOptions
    from gbdraw.session_request_codec import _encode_colors
    from gbdraw.session_resources import RequestResources, SessionResourceTable

    resources = RequestResources(SessionResourceTable())
    encoded = _encode_colors(ColorOptions(color_table_file="builtin:plastid"), resources=resources)
    assert encoded["colorTableFile"] is not None
    (resource,) = resources.added()
    assert Path(resource.source_path).read_bytes() == (TABLE_DIR / "plastid.tsv").read_bytes()


def test_a_synthetic_bakta_record_draws_the_cog_legend() -> None:
    record = SeqRecord(Seq("ATG" * 400), id="bakta_like", name="bakta_like")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(FeatureLocation(0, 300, strand=1), type="CDS",
                   qualifiers={"locus_tag": ["B_1"], "db_xref": ["COG:COG0593", "COG:L"]}),
        SeqFeature(FeatureLocation(600, 900, strand=-1), type="CDS",
                   qualifiers={"locus_tag": ["B_2"], "db_xref": ["COG:COG0001", "COG:E"]}),
    ]
    svg = gbdraw.draw_circular(
        [record], options=gbdraw.CircularOptions(features=gbdraw.FeatureOptions(color_table="builtin:cog"))
    ).to_svg()
    assert "Replication, recombination and repair" in svg
    assert "Amino acid transport and metabolism" in svg
    assert "Transcription" not in svg
