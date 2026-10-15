"""The proposal of ``gbdraw auto``: classification, decisions, and the command."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.cli_utils.auto import auto_main
from gbdraw.mode_profiles import DEFAULT_FEATURE_TYPES
from gbdraw.recommend import LABEL_ALL_MAX, LABEL_NAMED_MAX, inspect_record, propose

TEST_INPUTS = Path(__file__).resolve().parent / "test_inputs"


def _feature(
    kind: str,
    start: int,
    end: int,
    strand: int = 1,
    *,
    parts: Sequence[tuple[int, int]] | None = None,
    **qualifiers: str,
) -> SeqFeature:
    if parts:
        location = CompoundLocation([FeatureLocation(a, b, strand) for a, b in parts])
    else:
        location = FeatureLocation(start, end, strand)
    return SeqFeature(location, type=kind, qualifiers={key: [value] for key, value in qualifiers.items()})


def _record(
    *,
    length: int = 20_000,
    topology: str = "circular",
    organism: str = "Escherichia coli",
    taxonomy: Sequence[str] = ("Bacteria", "Pseudomonadati"),
    definition: str = "Escherichia coli strain X, complete genome",
    molecule_type: str = "DNA",
    source: Mapping[str, str] | None = None,
    features: Sequence[SeqFeature] = (),
    record_id: str = "TEST000001.1",
    comment: str = "",
) -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * (length // 4)), id=record_id, name=record_id, description=definition)
    record.annotations.update(
        {
            "topology": topology,
            "organism": organism,
            "taxonomy": list(taxonomy),
            "molecule_type": molecule_type,
            "comment": comment,
        }
    )
    source_feature = _feature("source", 0, length, organism=organism, **(source or {}))
    record.features = [source_feature, *features]
    return record


def _genes(
    count: int,
    *,
    kind: str = "CDS",
    span: int = 900,
    gap: int = 100,
    offset: int = 0,
    name: bool = True,
    strands: Sequence[int] = (1, -1),
) -> list[SeqFeature]:
    features = []
    for index in range(count):
        start = offset + index * (span + gap)
        qualifiers = {"product": f"protein {index}"}
        if name:
            qualifiers["gene"] = f"gene{index}"
        features.append(_feature(kind, start, start + span, strands[index % len(strands)], **qualifiers))
    return features


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        pytest.param(
            _record(source={"organelle": "mitochondrion"}, taxonomy=("Eukaryota", "Metazoa"), definition="Homo sapiens mitochondrion, complete genome"),
            "mitochondrion",
            id="mitochondrion",
        ),
        pytest.param(
            _record(source={"organelle": "plastid:chloroplast"}, taxonomy=("Eukaryota", "Viridiplantae"), definition="Nicotiana tabacum plastid, complete genome"),
            "plastid",
            id="plastid",
        ),
        pytest.param(_record(source={"plasmid": "pX"}, definition="Escherichia coli plasmid pX, complete sequence"), "plasmid", id="plasmid"),
        pytest.param(
            _record(topology="linear", organism="Escherichia phage Lambda", taxonomy=("Viruses", "Duplodnaviria"), definition="Enterobacteria phage lambda, complete genome"),
            "phage",
            id="phage",
        ),
        pytest.param(
            _record(topology="linear", organism="Majanivirus", taxonomy=("Viruses", "Naldaviricetes"), definition="Majanivirus DNA, complete genome"),
            "virus_dna",
            id="virus-dna",
        ),
        pytest.param(
            _record(topology="linear", organism="SARS-CoV-2", taxonomy=("Viruses", "Riboviria"), molecule_type="ss-RNA", definition="SARS-CoV-2, complete genome"),
            "virus_rna",
            id="virus-rna",
        ),
        pytest.param(
            _record(topology="linear", definition="Streptomyces sp. paromomycin biosynthetic gene cluster", length=40_000),
            "gene_cluster",
            id="gene-cluster",
        ),
        pytest.param(_record(), "bacterial_chromosome", id="bacterial-chromosome"),
        pytest.param(_record(topology="linear", definition="Escherichia coli contig_12", length=20_000), "contig", id="contig"),
        pytest.param(_record(topology="linear", organism="Homo sapiens", taxonomy=("Eukaryota", "Metazoa"), definition="Homo sapiens chromosome 1 region"), "eukaryotic_region", id="eukaryote"),
    ],
)
def test_record_kind(record: SeqRecord, expected: str) -> None:
    assert inspect_record(record).kind == expected


def test_small_circular_genome_is_labelled_with_gene_names() -> None:
    record = _record(
        source={"organelle": "mitochondrion"},
        taxonomy=("Eukaryota", "Metazoa"),
        definition="Homo sapiens mitochondrion, complete genome",
        features=[*_genes(13), _feature("D-loop", 15_000, 16_000)],
    )
    proposal = propose([record])
    assert (proposal.mode, proposal.grouping) == ("circular", "single")
    assert proposal.label_scope == "outer"
    assert proposal.track_type == "middle"
    assert proposal.qualifier_priority == (("CDS", "gene,product"),)
    assert "D-loop" in proposal.feature_types
    assert proposal.separate_strands and proposal.show_gc and proposal.show_skew
    assert proposal.color_rules_name == "mitochondrion_complexes"
    assert {decision.setting for decision in proposal.decisions} >= {"mode", "label_scope", "palette", "legend"}


def test_complete_bacterial_genome_with_linear_locus_is_drawn_circular() -> None:
    proposal = propose([_record(topology="linear", length=600_000, features=_genes(10))])
    assert proposal.mode == "circular"
    assert proposal.topology_override == "circular"
    assert "--record_topology" in proposal.cli_args()


def test_many_labels_are_hidden_and_named_genes_are_whitelisted() -> None:
    crowded = _record(length=2_000_000, features=_genes(LABEL_NAMED_MAX + 1, name=False))
    assert propose([crowded]).label_scope == "none"
    named = _record(length=1_000_000, features=_genes(LABEL_ALL_MAX + 1))
    proposal = propose([named])
    assert proposal.label_scope == "outer"
    assert proposal.label_whitelist[0] == ("CDS", "gene", ".+")
    assert proposal.track_type == "tuckin"


def test_hypothetical_products_are_left_out_of_the_labels() -> None:
    named = _genes(10)
    hypothetical = [
        _feature("CDS", 20_000 + 1_000 * i, 20_900 + 1_000 * i, product="hypothetical protein") for i in range(3)
    ]
    proposal = propose([_record(length=40_000, features=[*named, *hypothetical])])
    assert proposal.label_scope == "outer"
    assert proposal.label_blacklist == ("hypothetical", "uncharacterized", "unknown")
    args = proposal.cli_args()
    assert args[args.index("--label_blacklist") + 1] == "hypothetical,uncharacterized,unknown"
    assert "--label_whitelist" not in args
    # Only named genes are labelled above LABEL_ALL_MAX, so the whitelist replaces the blacklist.
    crowded = propose([_record(length=1_000_000, features=[*_genes(LABEL_ALL_MAX + 1), *hypothetical])])
    assert crowded.label_whitelist and not crowded.label_blacklist
    # Without hypothetical products there is nothing to leave out.
    assert not propose([_record(features=named)]).label_blacklist


def test_polyprotein_genome_draws_mature_peptides_on_separate_lanes() -> None:
    polyprotein = _feature("CDS", 100, 20_000, gene="ORF1ab", product="polyprotein")
    peptides = [_feature("mat_peptide", 100 + 1_000 * i, 100 + 1_000 * (i + 1), product=f"nsp{i + 1}") for i in range(10)]
    record = _record(
        topology="linear",
        organism="Coronavirus",
        taxonomy=("Viruses", "Riboviria"),
        molecule_type="ss-RNA",
        definition="Coronavirus, complete genome",
        features=[polyprotein, *peptides],
    )
    proposal = propose([record])
    assert proposal.mode == "linear"
    assert "mat_peptide" in proposal.feature_types
    assert proposal.resolve_overlaps
    assert not proposal.separate_strands
    assert ("mat_peptide", "product,gene") in proposal.qualifier_priority


def test_origin_spanning_joins_and_underlays_do_not_count_as_nested() -> None:
    wrapped = _feature("D-loop", 0, 0, parts=[(19_000, 20_000), (0, 500)])
    record = _record(
        source={"organelle": "mitochondrion"},
        taxonomy=("Eukaryota", "Metazoa"),
        definition="mitochondrion, complete genome",
        features=[*_genes(13, offset=1_000), wrapped],
    )
    assert not propose([record]).resolve_overlaps
    # A repeat_region is drawn as an underlay: the genes on it are not nested.
    inverted_repeat = _feature("repeat_region", 0, 15_000, note="inverted repeat A")
    assert not propose([_record(features=[*_genes(13), inverted_repeat])]).resolve_overlaps


def test_circular_records_share_one_canvas_and_rows_are_stacked_otherwise() -> None:
    chromosome = _record(length=3_000_000, record_id="CHR1", source={"chromosome": "1"})
    plasmid = _record(length=100_000, record_id="PL1", source={"plasmid": "pA"}, definition="plasmid pA, complete sequence")
    grid = propose([chromosome, plasmid])
    assert (grid.mode, grid.grouping) == ("circular", "grid")
    assert "--multi_record_canvas" in grid.cli_args()
    contigs = [_record(topology="linear", definition=f"contig_{i}", length=50_000, record_id=f"C{i}", features=_genes(10)) for i in range(3)]
    rows = propose(contigs)
    assert (rows.mode, rows.grouping, rows.label_scope) == ("linear", "single", "first")


@pytest.mark.parametrize(
    ("qualifiers", "tool", "rules_name", "first_rule"),
    [
        ({"note": "COG:E"}, "bakta", "bakta_cog", ("CDS", "note")),
        ({"phrog": "phrog_1", "function": "tail"}, "pharokka", "pharokka", ("CDS", "function")),
        ({"gene_kind": "biosynthetic"}, "antismash", "antismash_gene_kind", ("CDS", "gene_kind")),
    ],
)
def test_annotation_tool_selects_its_color_table(
    qualifiers: Mapping[str, str], tool: str, rules_name: str, first_rule: tuple[str, str]
) -> None:
    features = [_feature("CDS", 100 + 1_000 * i, 1_000 + 1_000 * i, product=f"p{i}", **qualifiers) for i in range(5)]
    record = _record(topology="linear", definition="Escherichia coli contig_1", features=features)
    facts = inspect_record(record)
    assert facts.annotation_tool == tool
    proposal = propose([record])
    assert proposal.color_rules_name == rules_name
    assert proposal.color_rules and proposal.color_rules[0][:2] == first_rule
    assert "color_table" in proposal.tables()


def test_gene_only_annotation_draws_the_gene_features() -> None:
    record = _record(features=_genes(5, kind="gene"))
    proposal = propose([record])
    assert "gene" in proposal.feature_types
    assert any("no CDS" in note for note in proposal.notes)


def test_cli_args_name_the_tables_and_default_feature_types_are_omitted() -> None:
    proposal = propose([_record(features=_genes(5))])
    assert proposal.feature_types == DEFAULT_FEATURE_TYPES
    args = proposal.cli_args({"qualifier_priority": Path("/tmp/p.tsv")})
    assert "-k" not in args
    assert args[args.index("--qualifier_priority") + 1] == "/tmp/p.tsv"
    assert args[-2:] == ["-l", "left"]
    payload = json.loads(proposal.to_json())
    assert payload["mode"] == "circular" and payload["records"][0]["kind"] == "bacterial_chromosome"


@pytest.mark.parametrize(
    ("filename", "mode", "color_rules", "expected_type"),
    [
        ("SARS-CoV-2.gbk", "linear", None, "mat_peptide"),
        ("BGC0000712.gbk", "linear", "antismash_gene_kind", "CDS"),
    ],
)
def test_bundled_inputs(filename: str, mode: str, color_rules: str | None, expected_type: str) -> None:
    from gbdraw.io.genome import load_gbks

    proposal = propose(load_gbks([str(TEST_INPUTS / filename)]))
    assert proposal.mode == mode
    assert proposal.color_rules_name == color_rules
    assert expected_type in proposal.feature_types


def test_auto_command_dry_run_prints_the_equivalent_command(capsys: pytest.CaptureFixture[str]) -> None:
    auto_main(["--gbk", str(TEST_INPUTS / "SARS-CoV-2.gbk"), "--dry_run"])
    output = capsys.readouterr().out
    assert "Proposed settings:" in output
    assert "gbdraw linear --gbk" in output
    assert "--resolve_overlaps" in output


def _json_block(output: str) -> dict[str, object]:
    """The pretty-printed JSON object in ``output``, around the INFO lines of a run."""
    lines = output.splitlines()
    first = lines.index("{")
    last = len(lines) - 1 - lines[::-1].index("}")
    payload = json.loads("\n".join(lines[first : last + 1]))
    assert isinstance(payload, dict)
    return payload


def test_auto_command_draws_and_writes_its_tables(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    prefix = tmp_path / "bgc"
    auto_main(["--gbk", str(TEST_INPUTS / "BGC0000712.gbk"), "--yes", "-o", str(prefix), "--json"])
    payload = _json_block(capsys.readouterr().out)
    assert payload["mode"] == "linear"
    assert payload["decisions"] and all({"setting", "value", "reason"} <= set(d) for d in payload["decisions"])
    assert (tmp_path / "bgc.svg").exists()
    color_table = tmp_path / "bgc.auto.color_table.tsv"
    assert color_table.read_text().splitlines()[0].split("\t")[:2] == ["CDS", "gene_kind"]
    assert (tmp_path / "bgc.auto.qualifier_priority.tsv").read_text() == "CDS\tgene,product\n"
