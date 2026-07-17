from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

from gbdraw.io.genome import load_gff_fasta


def test_load_gff_fasta_collapses_same_id_rows_into_compound_location(
    tmp_path: Path,
) -> None:
    gff_path = tmp_path / "multipart.gff3"
    fasta_path = tmp_path / "multipart.fasta"
    gff_path.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "##sequence-region record1 1 30",
                "record1\ttest\tgene\t1\t20\t.\t+\t.\tID=gene1",
                "record1\ttest\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1;Note=first",
                "record1\ttest\tCDS\t11\t16\t.\t+\t2\tID=cds1;Parent=gene1;Note=second",
                "",
            ]
        ),
        encoding="utf-8",
    )
    fasta_path.write_text(">record1\nATGAAACCCCGGGTTTAAACCCGGGTTTAA\n", encoding="utf-8")

    records = load_gff_fasta(
        [str(gff_path)],
        [str(fasta_path)],
        mode="linear",
        selected_features_set={"CDS"},
    )

    assert len(records[0].features) == 1
    cds = records[0].features[0]
    assert isinstance(cds.location, CompoundLocation)
    assert [(int(part.start), int(part.end)) for part in cds.location.parts] == [
        (0, 6),
        (10, 16),
    ]
    assert cds.location.operator == "join"
    assert cds.qualifiers["ID"] == ["cds1"]
    assert cds.qualifiers["Parent"] == ["gene1"]
    assert cds.qualifiers["Note"] == ["first", "second"]
    assert cds.qualifiers["phase"] == ["0", "2"]


def test_nc_013668_gff3_orf10_matches_genbank_compound_cds(
    test_inputs_dir: Path,
) -> None:
    gff_records = load_gff_fasta(
        [str(test_inputs_dir / "NC_013668.gff3")],
        [str(test_inputs_dir / "NC_013668.fasta")],
        mode="circular",
        selected_features_set={"CDS"},
    )
    genbank_record = SeqIO.read(test_inputs_dir / "NC_013668.gb", "genbank")

    gff_cds_features = gff_records[0].features
    gff_orf10 = next(
        feature
        for feature in gff_cds_features
        if feature.qualifiers.get("locus_tag") == ["AngHV1_ORF10"]
    )
    genbank_orf10 = next(
        feature
        for feature in genbank_record.features
        if feature.type == "CDS"
        and feature.qualifiers.get("locus_tag") == ["AngHV1_ORF10"]
    )

    assert len(gff_cds_features) == 134
    assert isinstance(gff_orf10.location, CompoundLocation)
    assert len(gff_orf10.location.parts) == 5
    assert gff_orf10.location == genbank_orf10.location
    assert gff_orf10.extract(gff_records[0].seq) == genbank_orf10.extract(genbank_record.seq)
