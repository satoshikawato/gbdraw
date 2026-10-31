"""Real-world input variants: every fixture in ``tests/fixtures/input_variants`` is read through the
public readers and drawn by both CLI modes. See that folder's README for the origin of each file."""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

import gbdraw
from gbdraw.cli import main as cli_main
from gbdraw.exceptions import GbdrawError

FIXTURES = Path(__file__).parent / "fixtures" / "input_variants"
FEATURE_ID = re.compile(r'data-gbdraw-feature-id="([^"]+)"')
KEY_QUALIFIERS = ("locus_tag", "ID", "Name", "gene", "product", "note")


@dataclass(frozen=True)
class Case:
    """What a fixture must read as, and what both CLI modes must do with it."""

    name: str
    kind: str = "gb"  # "gb" or "gff"
    fasta: str | None = None  # GFF3: partner file; None means "<stem>.fasta"
    records: int = 1
    lengths: tuple[int, ...] | None = None
    types: dict[str, int] | None = None  # feature type counts of the first record
    parts: dict[str, int] = field(default_factory=dict)  # "type:key" -> compound part count (default 1)
    strands: dict[str, int | None] = field(default_factory=dict)
    lens: dict[str, int] = field(default_factory=dict)  # "type:key" -> summed part length
    wraps: tuple[str, ...] = ()  # "type:key" of features whose shorter arc crosses the origin
    quals: dict[str, tuple[str, Any]] = field(default_factory=dict)  # "type:key" -> (qualifier, value or None)
    read_error: tuple[str, str] | None = None  # (exception class name, message regex)
    drawn_min: int = 1  # distinct feature ids both CLI modes must draw
    cli_error: str | None = None  # message regex of the "ERROR:" line both CLI modes must print


CASES = [
    Case("prokka_locus_tag_only.gb", types={"source": 1, "CDS": 2, "tRNA": 1, "rRNA": 1},
         strands={"CDS:PROKKA_00001": 1, "CDS:PROKKA_00002": -1}),
    Case("rast_seed_dbxref.gb", types={"source": 1, "CDS": 2, "misc_feature": 1, "rRNA": 1},
         quals={"CDS:Thr operon leader peptide": ("db_xref", "SEED:fig|83333.1.peg.1")}),
    Case("geneious_locus_spaces.gb", types={"source": 1, "CDS": 1, "gene": 1},
         quals={"CDS:GEN_0001": ("note", "first line of a long note that Geneious writes across lines; "
                                         "second line of the same note with a semicolon; third line closes the note")}),
    Case("artemis_colour.gb", types={"source": 1, "CDS": 3}, quals={"CDS:colA": ("colour", "2 125 0")}),
    Case("mibig_gene_kind.gb", types={"source": 1, "CDS": 3, "protocluster": 1},
         quals={"CDS:ctg1_1": ("gene_kind", "biosynthetic")}),
    Case("pgap_partial_locations.gb", types={"source": 1, "CDS": 4}, strands={"CDS:PGAP_00003": -1},
         lens={"CDS:PGAP_00001": 300, "CDS:PGAP_00002": 500}),
    Case("pgap_contig_no_sequence.gb", types={}, drawn_min=0, cli_error=r"(?i)sequence"),
    Case("geseq_plastome.gb", lengths=(6000,), types={"source": 1, "repeat_region": 2, "gene": 2, "CDS": 3, "tRNA": 1},
         parts={"CDS:rps12": 2, "gene:rps12": 2, "CDS:ycf3": 3, "tRNA:trnK-UUU": 2},
         strands={"CDS:rps12": None, "CDS:rbcL": -1},
         quals={"repeat_region:IRa": ("rpt_type", "inverted"), "CDS:rps12": ("trans_splicing", None)}),
    Case("ddbj_flavour.gb", types={"source": 1, "CDS": 2, "rRNA": 1}, strands={"CDS:ddbB": -1}),
    Case("easyfig72_embl.embl", lengths=(3000,), types={"source": 1, "CDS": 2},
         strands={"CDS:emA": 1, "CDS:emB": -1}, quals={"CDS:emB": ("product", "EMBL protein B")}),
    Case("bakta_cog_dbxref.gb", types={"source": 1, "gene": 1, "CDS": 3, "ncRNA": 1},
         quals={"CDS:BAK_00005": ("db_xref", "COG:V")}),
    Case("pharokka_phrog.gb", types={"source": 1, "CDS": 3, "tRNA": 1},
         quals={"CDS:PHAGE_CDS_0001": ("function", "DNA, RNA and nucleotide metabolism")}),
    Case("easyfig59_tabs.gb", types={"source": 1, "CDS": 2}, strands={"CDS:TAB_0002": -1},
         quals={"CDS:TAB_0001": ("product", "tab\tin\tvalue")}),
    Case("crlf_line_endings.gb", types={"source": 1, "CDS": 2}, quals={"CDS:crA": ("product", "crlf protein A")}),
    Case("cgviewjs3_one_bp_feature.gb", types={"source": 1, "misc_feature": 2, "CDS": 1, "variation": 1},
         strands={"misc_feature:single base, plus": 1, "misc_feature:single base, minus": -1},
         lens={"misc_feature:single base, plus": 1, "misc_feature:single base, minus": 1}),
    Case("pseudo_and_pseudogene.gb", types={"source": 1, "gene": 2, "CDS": 3},
         quals={"CDS:PS_0001": ("pseudo", None), "CDS:PS_0003": ("pseudogene", "unitary")}),
    Case("dfv85_complement_join.gb", types={"source": 1, "CDS": 4},
         parts={"CDS:CJ_0001": 2, "CDS:CJ_0002": 2, "CDS:CJ_0003": 2, "CDS:CJ_0004": 2},
         strands={"CDS:CJ_0001": -1, "CDS:CJ_0002": -1, "CDS:CJ_0003": -1, "CDS:CJ_0004": 1}),
    Case("dfv29_origin_spanning_join.gb", types={"source": 1, "CDS": 4},
         parts={"CDS:OR_0001": 2, "CDS:OR_0002": 2, "CDS:OR_0003": 2},
         strands={"CDS:OR_0001": 1, "CDS:OR_0002": -1, "CDS:OR_0003": -1},
         wraps=("CDS:OR_0001", "CDS:OR_0002", "CDS:OR_0003")),
    Case("cgview7_multi_record.gbff", records=3, lengths=(2000, 3500, 800), types={"source": 1, "CDS": 2}),
    Case("zero_features_record.gb", lengths=(1000,), types={}, drawn_min=0),
    Case("long_feature_over_10kb.gb", lengths=(12000,), types={"source": 1, "CDS": 2}, lens={"CDS:LG_0001": 10800}),
    Case("lowercase_sequence_n_runs.gb", lengths=(3000,), types={"source": 1, "CDS": 2}),
    Case("ggn35_embedded_fasta.gff3", kind="gff", fasta="ggn35_embedded_fasta.gff3", types={"gene": 2, "CDS": 2}),
    Case("ggn51_feature_without_id.gff3", kind="gff", types={"CDS": 2, "tRNA": 1},
         strands={"CDS:no id one": 1, "CDS:no id two": -1}),
    Case("ggn67_hash_in_attribute.gff3", kind="gff", types={"gene": 1, "CDS": 2},
         quals={"gene:gene1": ("Note", "see #123 and C#-like"), "CDS:cds1": ("product", "ABC #1 transporter")}),
    Case("ggn238_lowercase_cds_type.gff3", kind="gff", types={"CDS": 2, "gene": 1}),
    Case("dfv70_sequence_region_longer.gff3", kind="gff", lengths=(3000,), types={"gene": 1, "CDS": 1}),
    Case("ggn117_crlf.gff3", kind="gff", types={"gene": 2, "CDS": 2}, quals={"CDS:cds1": ("product", "protein A")}),
    Case("ggn227_fasta_header_tab.gff3", kind="gff", types={"gene": 2, "CDS": 2}),
    Case("ncbi_gene_mrna_cds.gff3", kind="gff", lengths=(3000,),
         types={"region": 1, "gene": 2, "mRNA": 1, "exon": 3, "CDS": 2},
         parts={"CDS:cds-A": 3}, strands={"CDS:cds-A": 1, "CDS:cds-B": -1}),
    Case("prokka_embedded_fasta.gff3", kind="gff", fasta="prokka_embedded_fasta.gff3",
         types={"CDS": 2, "tRNA": 1, "rRNA": 1}),
    Case("ggn65_augustus_exon_hierarchy.gff3", kind="gff",
         types={"gene": 2, "transcript": 2, "start_codon": 1, "CDS": 2, "exon": 3, "intron": 1, "stop_codon": 1},
         parts={"CDS:g1.t1.cds": 2}, strands={"CDS:g2.t1.cds": -1}),
    # GFF3 joins rows into one feature by a shared ID; these CDS rows have none, so they stay two features.
    Case("ggn93_exonerate_style.gff3", kind="gff",
         types={"similarity": 1, "gene": 1, "CDS": 2, "exon": 2, "splice5": 1, "intron": 1, "splice3": 1},
         quals={"gene:gene1": ("sequence", "query_1"), "exon:@101": ("deletions", "0"),
                "splice3:@599": ("splice_site", "AG")}),
    Case("score_and_phase_dot.gff3", kind="gff", types={"CDS": 2, "misc_feature": 1}),
    Case("pgv54_no_strand.gff3", kind="gff", types={"CDS": 3},
         strands={"CDS:cds1": None, "CDS:cds2": None, "CDS:cds3": 1}),
    Case("clinker67_single_line_fasta.gff3", kind="gff", lengths=(3000,), types={"gene": 2, "CDS": 2}),
    Case("gff3_fasta_missing.gff3", kind="gff", read_error=("InputFileError", r"does not exist"),
         cli_error=r"does not exist"),
]

IDS = [case.name for case in CASES]


def _paths(case: Case) -> tuple[list[str], list[str] | None]:
    gff_or_gbk = FIXTURES / case.name
    if case.kind == "gb":
        return [str(gff_or_gbk)], None
    partner = FIXTURES / (case.fasta or (Path(case.name).stem + ".fasta"))
    return [str(gff_or_gbk)], [str(partner)]


def _read(case: Case):
    first, second = _paths(case)
    if second is None:
        return gbdraw.read_genbank(first)
    return gbdraw.read_gff(first, second)


def _key(feature) -> str:
    for name in KEY_QUALIFIERS:
        values = feature.qualifiers.get(name)
        if values:
            return f"{feature.type}:{values[0]}"
    return f"{feature.type}:@{int(feature.location.start) + 1}"


def _wraps_origin(feature, length: int) -> bool:
    spans = sorted((int(part.start), int(part.end)) for part in feature.location.parts)
    if len(spans) < 2:
        return False
    inner = max(nxt[0] - prev[1] for prev, nxt in zip(spans, spans[1:]))
    outer = length - (spans[-1][1] - spans[0][0])
    return outer < inner


def _type_counts(record) -> dict[str, int]:
    counts: dict[str, int] = {}
    for feature in record.features:
        counts[feature.type] = counts.get(feature.type, 0) + 1
    return counts


@pytest.mark.parametrize("case", CASES, ids=IDS)
def test_variant_reads_as_expected(case: Case) -> None:
    if case.read_error is not None:
        exception_name, message = case.read_error
        with pytest.raises(GbdrawError) as caught:
            _read(case)
        assert type(caught.value).__name__ == exception_name
        assert re.search(message, str(caught.value)), str(caught.value)
        return

    records = _read(case)
    assert len(records) == case.records
    if case.lengths is not None:
        assert tuple(len(record.seq) for record in records) == case.lengths
    first = records[0]
    if case.types is not None:
        assert _type_counts(first) == case.types

    by_key = {_key(feature): feature for feature in first.features}
    for key, expected in case.parts.items():
        assert len(by_key[key].location.parts) == expected, key
    for key, feature in by_key.items():
        if key not in case.parts:
            assert len(feature.location.parts) == 1, key
    for key, strand in case.strands.items():
        assert by_key[key].location.strand == strand, key
    for key, length in case.lens.items():
        assert sum(len(part) for part in by_key[key].location.parts) == length, key
    wrapped = tuple(sorted(k for k, f in by_key.items() if _wraps_origin(f, len(first.seq))))
    assert wrapped == tuple(sorted(case.wraps))
    for key, (name, value) in case.quals.items():
        values = by_key[key].qualifiers[name]
        if value is not None:
            assert value in values or " ".join(" ".join(values).split()) == value, (key, values)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("case", CASES, ids=IDS)
def test_variant_draws_or_names_its_failure(
    case: Case, mode: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    first, second = _paths(case)
    inputs = ["--gbk", *first] if second is None else ["--gff", *first, "--fasta", *second]
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["gbdraw", mode, *inputs, "-o", str(out), "-f", "svg"])
    exit_code = 0
    try:
        cli_main()
    except SystemExit as stop:
        exit_code = stop.code if isinstance(stop.code, int) else 1
    except Exception as crash:  # an unnamed failure: the CLI printed no "ERROR:" line for it
        exit_code = 1
        print(f"UNCAUGHT {type(crash).__name__}: {crash}", file=sys.stderr)
    stderr = capsys.readouterr().err

    expects_error = case.read_error is not None or case.cli_error is not None
    if expects_error:
        assert exit_code != 0, "expected a named error"
        error_lines = [line for line in stderr.splitlines() if line.startswith("ERROR:")]
        assert error_lines, stderr[-400:]
        assert re.search(case.cli_error or "", error_lines[-1]), error_lines[-1]
        return

    assert exit_code == 0, stderr[-400:]
    # A multi-record circular run writes one SVG per record.
    drawn: set[str] = set()
    for svg_path in tmp_path.glob("out*.svg"):
        drawn.update(FEATURE_ID.findall(svg_path.read_text()))
    assert len(drawn) >= case.drawn_min, f"drew {len(drawn)} feature(s)"


def test_readme_lists_every_fixture() -> None:
    readme = (FIXTURES / "README.md").read_text()
    on_disk = {path.name for path in FIXTURES.iterdir() if path.name not in {"README.md", ".gitattributes"} and path.suffix != ".fasta"}
    assert on_disk == {case.name for case in CASES}
    for name in on_disk:
        assert f"`{name}`" in readme, name
