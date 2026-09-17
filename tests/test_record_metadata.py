from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import (
    RecordSourceMetadata,
    format_inferred_definition,
    format_inferred_subtitle,
    infer_record_source_metadata,
)


def test_format_inferred_definition_binomial():
    meta = RecordSourceMetadata(organism="Escherichia coli O157:H7", strain="Sakai", replicon=None, organelle=None)
    assert format_inferred_definition(meta) == "<i>Escherichia coli</i> O157:H7 Sakai"


def test_format_inferred_definition_single_word():
    meta = RecordSourceMetadata(organism="Streptomyces", strain="sp. A", replicon=None, organelle=None)
    assert format_inferred_definition(meta) == "<i>Streptomyces</i> sp. A"


def test_format_inferred_definition_candidatus():
    meta = RecordSourceMetadata(organism="Candidatus Tyloplasma litorale", strain="Fukuoka2020", replicon=None, organelle=None)
    assert format_inferred_definition(meta) == "Candidatus <i>Tyloplasma litorale</i> Fukuoka2020"


def test_format_inferred_definition_non_organism():
    meta = RecordSourceMetadata(organism="synthetic construct", strain="cloning vector pUC19", replicon=None, organelle=None)
    assert format_inferred_definition(meta) == "cloning vector pUC19"


def test_format_inferred_subtitle_replicon_and_organelle():
    meta = RecordSourceMetadata(organism="E. coli", strain="", replicon="Chromosome 1", organelle=None)
    assert format_inferred_subtitle(meta) == "Chromosome 1"

    meta2 = RecordSourceMetadata(organism="H. sapiens", strain="", replicon=None, organelle="mitochondrion")
    assert format_inferred_subtitle(meta2) == "Mitochondrion"


def test_format_inferred_subtitle_from_description():
    meta = RecordSourceMetadata(organism="Escherichia coli", strain="K-12", replicon=None, organelle=None)
    assert format_inferred_subtitle(meta, "Escherichia coli str. K-12 substr. MG1655, complete genome.") == "Complete genome"
    assert format_inferred_subtitle(meta, "Homo sapiens mitochondrion, complete genome.") == "Mitochondrion, complete genome"
    assert format_inferred_subtitle(meta, "Escherichia coli plasmid pOSAK1, complete sequence.") == "Plasmid pOSAK1"


def test_infer_record_source_metadata():
    record = SeqRecord(Seq("ATGC"), id="test", annotations={"organism": "Bacillus subtilis"})
    feature = SeqFeature(
        SimpleLocation(0, 4),
        type="source",
        qualifiers={"isolate": ["168"], "chromosome": ["1"]}
    )
    record.features.append(feature)
    meta = infer_record_source_metadata(record)
    assert meta.organism == "Bacillus subtilis"
    assert meta.strain == "168"
    assert meta.replicon == "Chromosome 1"
    assert meta.organelle is None

