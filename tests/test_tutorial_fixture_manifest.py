from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO

from gbdraw.io.cli_tables import read_comparisons_table, read_records_table


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
MANIFEST_PATH = FIXTURE_ROOT / "manifest.json"


def _load_manifest() -> dict[str, object]:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _strand_counts(features) -> dict[str, int]:
    counts = {"positive": 0, "negative": 0, "unknown": 0}
    for feature in features:
        strand = feature.location.strand
        if strand == 1:
            counts["positive"] += 1
        elif strand == -1:
            counts["negative"] += 1
        else:
            counts["unknown"] += 1
    return counts


def test_tutorial_fixture_manifest_is_complete_unique_and_within_budget() -> None:
    assert "tests/test_inputs" not in MANIFEST_PATH.read_text(encoding="utf-8")
    manifest = _load_manifest()
    fixtures = manifest["fixtures"]
    files = manifest["files"]
    budgets = manifest["packageBudgets"]

    assert manifest["schemaVersion"] == 2
    assert manifest["canonicalRoot"] == "gbdraw/web/tutorial-data"

    declared_paths: set[str] = set()
    declared_hashes: set[str] = set()
    declared_scenarios: set[str] = set()
    for file_id, metadata in files.items():
        fixture_id = metadata["fixtureId"]
        relative_path = metadata["relativePath"]
        path = FIXTURE_ROOT / relative_path
        provenance = metadata["provenance"]

        assert fixture_id in fixtures, file_id
        assert file_id in fixtures[fixture_id]["fileIds"], file_id
        assert metadata["packageTier"] in budgets, file_id
        assert metadata["role"] in {"raw", "derived", "rule"}, file_id
        assert metadata["mediaType"] and metadata["inputType"], file_id
        assert metadata["scenarioIds"], file_id
        assert not Path(relative_path).is_absolute(), file_id
        assert ".." not in Path(relative_path).parts, file_id
        assert path.is_file(), file_id
        assert path.stat().st_size == metadata["sizeBytes"], file_id
        assert _sha256(path.read_bytes()) == metadata["sha256"], file_id

        assert provenance["sourceName"], file_id
        assert provenance["sourceUrl"].startswith("https://"), file_id
        assert "retrievedOn" in provenance, file_id
        retrieved_on = provenance["retrievedOn"]
        if retrieved_on is not None:
            assert date.fromisoformat(retrieved_on).isoformat() == retrieved_on, file_id
        assert provenance["repositoryAddedOn"], file_id
        assert provenance["licenseNote"], file_id
        if metadata["role"] == "derived":
            derivation = metadata["derivation"]
            assert derivation["status"] == "reproducible", file_id
            assert (REPO_ROOT / derivation["script"]).is_file(), file_id
            assert derivation["sourceFileIds"] and derivation["sourceChecksums"], file_id
            assert len(derivation["sourceFileIds"]) == len(
                derivation["sourceChecksums"]
            ), file_id
            for source_file_id, source_checksum in zip(
                derivation["sourceFileIds"],
                derivation["sourceChecksums"],
                strict=True,
            ):
                assert source_file_id in files, file_id
                assert files[source_file_id]["sha256"] == source_checksum, file_id
            for external_source in derivation.get("externalSources", []):
                assert external_source["id"], file_id
                assert external_source["filename"], file_id
                assert external_source["sizeBytes"] > 0, file_id
                assert len(external_source["sha256"]) == 64, file_id
            assert derivation["tool"] and derivation["toolVersion"], file_id
            assert derivation["arguments"] and derivation["outputStatistics"], file_id

        assert relative_path not in declared_paths, file_id
        assert metadata["sha256"] not in declared_hashes, file_id
        declared_paths.add(relative_path)
        declared_hashes.add(metadata["sha256"])
        declared_scenarios.update(metadata["scenarioIds"])

    actual_paths = {
        path.relative_to(FIXTURE_ROOT).as_posix()
        for path in FIXTURE_ROOT.rglob("*")
        if path.is_file() and path != MANIFEST_PATH
    }
    assert actual_paths == declared_paths

    referenced_file_ids: set[str] = set()
    for fixture_id, fixture in fixtures.items():
        assert fixture["packageTier"] in budgets, fixture_id
        owned_bytes = sum(files[file_id]["sizeBytes"] for file_id in fixture["fileIds"])
        assert owned_bytes <= fixture["sizeBudgetBytes"], fixture_id
        referenced_file_ids.update(fixture["fileIds"])
        referenced_file_ids.update(fixture["fileReferences"])
    assert referenced_file_ids == set(files)

    for tier, budget in budgets.items():
        tier_files = [
            metadata
            for metadata in files.values()
            if metadata["packageTier"] == tier
        ]
        actual_bytes = sum(metadata["sizeBytes"] for metadata in tier_files)
        assert actual_bytes == budget["actualUncompressedBytes"], tier
        assert len(tier_files) == budget["fileCount"], tier
        assert actual_bytes <= budget["maxUncompressedBytes"], tier

    scenario_manifest = json.loads(
        (REPO_ROOT / "docs" / "scenarios" / "manifest.json").read_text(encoding="utf-8")
    )
    known_scenarios = {chapter["id"] for chapter in scenario_manifest["chapters"]}
    assert declared_scenarios <= known_scenarios


def test_tutorial_fixture_retrieval_dates_follow_the_schema_v2_policy() -> None:
    manifest = _load_manifest()
    policy = manifest["provenancePolicy"]
    legacy_cutoff = date.fromisoformat(
        policy["retrievalDateRequiredForRepositoryAdditionsAfter"]
    )

    assert policy == {
        "retrievalDateRequiredForRepositoryAdditionsAfter": "2026-08-04",
        "legacyUnknownStatus": "unknown-legacy",
    }
    for file_id, metadata in manifest["files"].items():
        provenance = metadata["provenance"]
        repository_added_on = date.fromisoformat(provenance["repositoryAddedOn"])
        retrieved_on = provenance["retrievedOn"]

        assert len(metadata["sha256"]) == 64, file_id
        if retrieved_on is None:
            assert provenance.get("retrievalDateStatus") == "unknown-legacy", file_id
            assert repository_added_on <= legacy_cutoff, file_id
        else:
            assert date.fromisoformat(retrieved_on) <= repository_added_on, file_id
            assert "retrievalDateStatus" not in provenance, file_id


def test_bgc_comparison_and_circular_multi_record_fixtures_are_distinct() -> None:
    manifest = _load_manifest()
    files = manifest["files"]
    bgc_five = manifest["fixtures"]["aminoglycoside-bgc-five"]
    multi_record = manifest["fixtures"]["metazoan-mitochondria-four"]

    assert "aminoglycoside-bgc-pair" not in manifest["fixtures"]
    bgc_file_ids = [
        "bgc-0000708-genbank",
        "bgc-0000709-genbank",
        "bgc-0000711-genbank",
        "bgc-0000712-genbank",
        "bgc-0000713-genbank",
    ]
    assert bgc_five["fileIds"][:5] == bgc_file_ids
    bgc_semantics = bgc_five["expectedSemantics"]
    assert bgc_semantics["recordCount"] == 5
    assert bgc_semantics["recordIds"] == [
        "BGC0000708",
        "BGC0000709",
        "BGC0000711",
        "BGC0000712",
        "BGC0000713",
    ]
    assert bgc_semantics["recordsAreWholeCanonicalSources"] is True
    assert bgc_semantics["galleryAlignedLosatpMode"] == "similarity-groups"
    assert bgc_semantics["tutorialScenarioId"] == "T-GUI-04"
    assert bgc_semantics["similarityGroupScenarioId"] == "H-GUI-07"
    assert bgc_semantics["collinearScenarioId"] == "H-CLI-08"
    assert bgc_semantics["losatpQualificationStatus"] == "qualified"
    assert bgc_semantics["rawAdjacentResultRows"] == 232
    assert bgc_semantics["similarityGroupCount"] == 23
    assert bgc_semantics["similarityGroupLinkCount"] == 77
    assert bgc_semantics["collinearBlockCount"] == 7
    assert {"T-GUI-03", "H-GUI-04", "H-GUI-05", "H-GUI-06"}.isdisjoint(
        bgc_five["scenarioIds"]
    )
    for file_id in bgc_file_ids:
        assert {"T-GUI-03", "H-GUI-04", "H-GUI-05", "H-GUI-06"}.isdisjoint(
            files[file_id]["scenarioIds"]
        )

    circular_scenarios = {"H-GUI-02", "H-CLI-03", "H-PY-01"}
    assert circular_scenarios.isdisjoint(bgc_five["scenarioIds"])
    for file_id in bgc_file_ids:
        assert circular_scenarios.isdisjoint(files[file_id]["scenarioIds"])

    assert "multi-record-small" not in manifest["fixtures"]
    assert multi_record["fileIds"] == [
        "danio-mitochondrion-genbank",
        "drosophila-mitochondrion-genbank",
        "caenorhabditis-mitochondrion-genbank",
    ]
    assert multi_record["fileReferences"] == ["human-mitochondrion-genbank"]
    assert multi_record["expectedSemantics"]["recordsAreWholeCanonicalSources"] is True
    assert multi_record["expectedSemantics"]["comparisonEvidence"] is False
    assert multi_record["expectedSemantics"]["intendedUse"] == "layout-only"
    file_ids = [*multi_record["fileReferences"], *multi_record["fileIds"]]
    records = [
        SeqIO.read(FIXTURE_ROOT / files[file_id]["relativePath"], "genbank")
        for file_id in file_ids
    ]
    assert [record.id for record in records] == multi_record["expectedSemantics"]["recordIds"]
    assert [len(record) for record in records] == multi_record["expectedSemantics"]["recordLengths"]
    assert [record.annotations.get("topology") for record in records] == ["circular"] * 4
    assert all("complete genome" in record.description for record in records)


def test_hepatoplasmataceae_fixture_supports_gallery_collinear_tutorial() -> None:
    manifest = _load_manifest()
    fixture = manifest["fixtures"]["hepatoplasmataceae-five"]
    semantics = fixture["expectedSemantics"]
    scenario_manifest = json.loads(
        (REPO_ROOT / "docs" / "scenarios" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    scenario = next(
        chapter
        for chapter in scenario_manifest["chapters"]
        if chapter["id"] == "T-GUI-08"
    )

    assert fixture["scenarioIds"] == [
        "T-GUI-08",
        "T-CLI-10",
        "T-PY-07",
        "H-GUI-08",
    ]
    assert semantics["recordIds"] == [
        "AP027078.1",
        "AP027131.1",
        "AP027133.1",
        "AP027132.1",
        "NZ_CP006932.1",
    ]
    assert semantics["recordsAreWholeCanonicalSources"] is True
    assert semantics["comparisonProgram"] == "LOSATP"
    assert semantics["galleryCollinearSearchScope"] == "adjacent"
    assert semantics["gallerySearchedRecordPairCount"] == 4
    assert semantics["allRecordHowToSearchScope"] == "all"
    assert semantics["allRecordHowToPairCount"] == 10
    assert semantics["displayedRibbonScope"] == "adjacent"
    assert semantics["recordOrderMatchesGallery"] is True
    assert scenario["settings"] == {
        "mode": "linear",
        "program": "losatp",
        "protein_mode": "collinear",
        "minimum_anchors": 1,
        "comparison_scope": "adjacent",
        "scheduling": "auto",
        "total_threads": "safe",
        "parallel_runs": "auto",
        "threads_per_run": "auto",
    }


def test_bgc_pair_is_not_assigned_to_nucleotide_comparison_scenarios() -> None:
    sources = (
        REPO_ROOT / "docs" / "internal" / "DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md",
        REPO_ROOT / "docs" / "scenarios" / "manifest.json",
        MANIFEST_PATH,
    )
    banned_terms = (
        "aminoglycoside-bgc-pair",
        "canonical BGC pair",
        "same two-record BGC pair",
        "small BGC pair",
        "frozen BGC pair",
    )
    offenders: list[str] = []
    for path in sources:
        text = path.read_text(encoding="utf-8")
        for term in banned_terms:
            if term.lower() in text.lower():
                offenders.append(f"{path.relative_to(REPO_ROOT)}: {term}")
    assert offenders == []


def test_lambda_de3_comparison_fixture_uses_complete_sources_and_fixed_evidence() -> None:
    manifest = _load_manifest()
    files = manifest["files"]
    comparison = manifest["fixtures"]["lambda-de3-comparison"]

    assert "t4" not in manifest["fixtures"]
    assert comparison["fileReferences"] == ["lambda-genbank", "de3-genbank"]
    assert comparison["expectedSemantics"]["queryRecordId"] == "NC_001416.1"
    assert comparison["expectedSemantics"]["subjectRecordId"] == "NC_042057.1"
    assert comparison["expectedSemantics"]["recordsAreWholeCanonicalSources"] is True

    de3 = files["de3-genbank"]
    assert de3["relativePath"] == "de3/NC_042057.1.gb"
    assert de3["records"] == [
        {
            "id": "NC_042057.1",
            "accession": "NC_042057",
            "sequenceVersion": 1,
            "organism": "Escherichia phage DE3",
            "description": "Enterobacteria phage DE3, complete genome",
            "topology": "linear",
            "length": 42925,
            "sequenceSha256": "533c360de52e4d5ac1f5dc31e0b927a518ff04951676d287b4b6d0238758aa45",
            "featureCount": 119,
            "displayedFeatureCount": 57,
            "featureCounts": {
                "CDS": 57,
                "gene": 58,
                "misc_feature": 2,
                "regulatory": 1,
                "source": 1,
            },
            "cdsStrandCounts": {"positive": 37, "negative": 20, "unknown": 0},
        }
    ]

    losatn = files["lambda-de3-losatn"]
    tlosatx = files["lambda-de3-tlosatx"]
    assert losatn["derivation"]["sourceFileIds"] == ["lambda-genbank", "de3-genbank"]
    assert tlosatx["derivation"]["sourceFileIds"] == ["lambda-genbank", "de3-genbank"]
    assert losatn["expectedSemantics"]["rawRowCount"] == 6
    assert losatn["expectedSemantics"]["defaultRetainedRowCount"] == 6
    assert tlosatx["expectedSemantics"]["rawRowCount"] == 397
    assert tlosatx["expectedSemantics"]["defaultRetainedRowCount"] == 266
    assert tlosatx["expectedSemantics"]["tutorialRetainedRowCount"] == 7


def test_metazoan_mtdna_comparison_uses_complete_sources_and_fixed_evidence() -> None:
    manifest = _load_manifest()
    files = manifest["files"]
    comparison = manifest["fixtures"]["metazoan-mitochondria-comparison"]
    semantics = comparison["expectedSemantics"]

    assert comparison["fileReferences"] == [
        "human-mitochondrion-genbank",
        "danio-mitochondrion-genbank",
        "drosophila-mitochondrion-genbank",
        "caenorhabditis-mitochondrion-genbank",
    ]
    assert semantics["referenceRecordId"] == "NC_012920.1"
    assert semantics["referenceSide"] == "subject"
    assert semantics["sourceTopologies"] == ["circular"] * 4
    assert semantics["recordsAreWholeCanonicalSources"] is True
    assert semantics["totalRawRows"] == 435
    assert semantics["retainedRows"] == 106
    assert semantics["retainedUnionCoverageBp"] == 9813

    expected = {
        "danio-human-tlosatx": ("NC_002333.2", 276, 68),
        "drosophila-human-tlosatx": ("NC_024511.2", 93, 24),
        "caenorhabditis-human-tlosatx": ("NC_001328.1", 66, 14),
    }
    for file_id, (query_id, raw_rows, retained_rows) in expected.items():
        entry = files[file_id]
        evidence = entry["expectedSemantics"]
        assert evidence["program"] == "TLOSATX"
        assert evidence["queryRecordId"] == query_id
        assert evidence["subjectRecordId"] == "NC_012920.1"
        assert evidence["referenceSide"] == "subject"
        assert evidence["rawRowCount"] == raw_rows
        assert evidence["retainedRowCount"] == retained_rows
        assert entry["derivation"]["arguments"]["threads"] == 1
        assert entry["derivation"]["arguments"]["runs"] == 2


def test_metazoan_mtdna_fastas_match_complete_genbank_sources() -> None:
    files = _load_manifest()["files"]
    pairs = (
        ("danio-mitochondrion-fasta", "danio-mitochondrion-genbank"),
        ("drosophila-mitochondrion-fasta", "drosophila-mitochondrion-genbank"),
        ("caenorhabditis-mitochondrion-fasta", "caenorhabditis-mitochondrion-genbank"),
    )

    for fasta_id, genbank_id in pairs:
        fasta_metadata = files[fasta_id]
        fasta = SeqIO.read(FIXTURE_ROOT / fasta_metadata["relativePath"], "fasta")
        source = SeqIO.read(FIXTURE_ROOT / files[genbank_id]["relativePath"], "genbank")
        assert fasta.id == source.id == fasta_metadata["records"][0]["id"]
        assert str(fasta.seq) == str(source.seq).upper()
        assert len(fasta) == fasta_metadata["records"][0]["length"]


def test_majanivirus_tables_resolve_complete_records_and_comparison_endpoints() -> None:
    manifest = _load_manifest()
    files = manifest["files"]
    fixture = manifest["fixtures"]["majanivirus-table-comparison"]
    semantics = fixture["expectedSemantics"]
    root = FIXTURE_ROOT / "majanivirus-table-comparison"

    records_table = read_records_table(str(root / "records.tsv"))
    assert records_table.input_kind == "gbk"
    assert records_table.record_ids == semantics["recordIds"]
    assert records_table.reverse_flags == semantics["reverseComplement"]
    assert [row.order for row in records_table.rows] == [1, 2, 3, 4]
    assert [row.row for row in records_table.rows] == semantics["rows"]
    assert [row.column for row in records_table.rows] == semantics["columns"]
    assert len(set(records_table.record_ids)) == len(records_table.rows)

    source_records = [SeqIO.read(path, "genbank") for path in records_table.gbk_files]
    assert [record.id for record in source_records] == semantics["recordIds"]
    assert [len(record) for record in source_records] == semantics["recordLengths"]
    assert [record.annotations.get("topology") for record in source_records] == ["linear"] * 4
    assert all("complete genome" in record.description for record in source_records)

    comparisons = read_comparisons_table(str(root / "comparisons.tsv"))
    endpoints = [[row.query, row.subject] for row in comparisons.rows]
    assert endpoints == semantics["comparisonEndpoints"]
    assert all(endpoint in set(records_table.record_ids) for pair in endpoints for endpoint in pair)
    assert [Path(row.blast).name for row in comparisons.rows] == [
        "MjeNMV.MelaMJNV.tblastx.out",
        "PemoMJNVA.PeseMJNV.tblastx.out",
    ]

    record_lengths = dict(zip(semantics["recordIds"], semantics["recordLengths"], strict=True))
    evidence_ids = (
        "majanivirus-mjenmv-melamjnv-tblastx",
        "majanivirus-pemomjnva-pesemjnv-tblastx",
    )
    for row, evidence_id in zip(comparisons.rows, evidence_ids, strict=True):
        expected = files[evidence_id]["expectedSemantics"]
        data_rows = [
            line.split("\t")
            for line in Path(row.blast).read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        assert len(data_rows) == expected["rawRowCount"]
        assert {len(cells) for cells in data_rows} == {expected["columnCount"]}
        assert {(cells[0], cells[1]) for cells in data_rows} == {(row.query, row.subject)}
        assert all(
            1 <= int(cells[6]) <= record_lengths[row.query]
            and 1 <= int(cells[7]) <= record_lengths[row.query]
            and 1 <= int(cells[8]) <= record_lengths[row.subject]
            and 1 <= int(cells[9]) <= record_lengths[row.subject]
            for cells in data_rows
        )
        retained = [
            cells
            for cells in data_rows
            if float(cells[2]) >= expected["identityMinimum"]
            and int(cells[3]) >= expected["alignmentLengthMinimum"]
        ]
        assert len(retained) == expected["tutorialRetainedRowCount"]


def test_metazoan_mtdna_comparison_builder_is_byte_reproducible() -> None:
    subprocess.run(
        [sys.executable, "tools/build_metazoan_mitochondria_comparison_fixture.py"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_genbank_manifest_semantics_match_parsed_records() -> None:
    files = _load_manifest()["files"]
    genbank_files = {
        file_id: metadata
        for file_id, metadata in files.items()
        if metadata["inputType"] == "genbank"
    }

    for file_id, metadata in genbank_files.items():
        records = list(SeqIO.parse(FIXTURE_ROOT / metadata["relativePath"], "genbank"))
        assert len(records) == len(metadata["records"]), file_id
        for record, expected in zip(records, metadata["records"], strict=True):
            assert record.id == expected["id"], file_id
            assert len(record) == expected["length"], file_id
            assert record.description == expected["description"], file_id
            assert record.annotations.get("organism") == expected["organism"], file_id
            assert record.annotations.get("topology") == expected["topology"], file_id
            assert _sha256(str(record.seq).upper().encode("ascii")) == expected["sequenceSha256"]
            assert len(record.features) == expected["featureCount"], file_id

            feature_counts = Counter(feature.type for feature in record.features)
            if "featureCounts" in expected:
                assert dict(feature_counts) == expected["featureCounts"], file_id
            cds_features = [feature for feature in record.features if feature.type == "CDS"]
            if "cdsCount" in expected:
                assert len(cds_features) == expected["cdsCount"], file_id
            if "cdsStrandCounts" in expected:
                assert _strand_counts(cds_features) == expected["cdsStrandCounts"], file_id
            if "strandCounts" in expected:
                assert _strand_counts(record.features) == expected["strandCounts"], file_id

            displayed_types = set(
                expected.get(
                    "displayedFeatureTypes",
                    (
                        ["CDS", "rRNA", "tRNA"]
                        if "mitochondrion" in expected["description"].lower()
                        else ["CDS"]
                    ),
                )
            )
            displayed_count = sum(
                count for feature_type, count in feature_counts.items() if feature_type in displayed_types
            )
            assert displayed_count == expected["displayedFeatureCount"], file_id

    bgc = _load_manifest()["fixtures"]["aminoglycoside-bgc-five"]
    displayed_bgc_features = sum(
        files[file_id]["records"][0]["displayedFeatureCount"]
        for file_id in bgc["fileIds"]
        if files[file_id]["inputType"] == "genbank"
    )
    assert displayed_bgc_features == bgc["expectedSemantics"]["displayedFeatureCount"] == 155


def test_lambda_gff3_and_fasta_semantics_match_the_source_record() -> None:
    files = _load_manifest()["files"]
    gff_metadata = files["lambda-full-record-gff3"]
    fasta_metadata = files["lambda-full-record-fasta"]
    fasta_records = list(SeqIO.parse(FIXTURE_ROOT / fasta_metadata["relativePath"], "fasta"))

    assert [record.id for record in fasta_records] == [
        record["id"] for record in fasta_metadata["records"]
    ]
    for record, expected in zip(fasta_records, fasta_metadata["records"], strict=True):
        assert len(record) == expected["length"]
        assert _sha256(str(record.seq).upper().encode("ascii")) == expected["sequenceSha256"]

    source = SeqIO.read(FIXTURE_ROOT / files["lambda-genbank"]["relativePath"], "genbank")
    assert len(fasta_records) == 1
    assert fasta_records[0].id == source.id == "NC_001416.1"
    assert str(fasta_records[0].seq) == str(source.seq)

    sequence_regions: dict[str, int] = {}
    rows: list[list[str]] = []
    for line in (FIXTURE_ROOT / gff_metadata["relativePath"]).read_text(encoding="utf-8").splitlines():
        if line.startswith("##sequence-region "):
            _, record_id, start, end = line.split()
            assert start == "1"
            sequence_regions[record_id] = int(end)
        elif line and not line.startswith("#"):
            columns = line.split("\t")
            assert len(columns) == 9
            rows.append(columns)

    assert list(sequence_regions) == [record["id"] for record in gff_metadata["records"]]
    assert sequence_regions == {record.id: len(record) for record in fasta_records}
    assert list(dict.fromkeys(row[0] for row in rows)) == list(sequence_regions)

    source_features = [
        feature for feature in source.features if feature.type in {"gene", "CDS"}
    ]
    assert len(source_features) == len(rows)
    for feature, row in zip(source_features, rows, strict=True):
        assert row[0] == source.id
        assert row[2] == feature.type
        assert int(row[3]) == int(feature.location.start) + 1
        assert int(row[4]) == int(feature.location.end)
        assert row[6] == ("+" if feature.location.strand == 1 else "-")
        expected_phase = (
            str(int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1)
            if feature.type == "CDS"
            else "."
        )
        assert row[7] == expected_phase
        if feature.type == "CDS":
            attributes = {
                key: unquote(value)
                for key, value in (item.split("=", 1) for item in row[8].split(";"))
            }
            assert attributes["translation"] == feature.qualifiers["translation"][0]

    translations = 0
    parent_ids: set[str] = set()
    gene_ids: set[str] = set()
    for record_id, _, feature_type, start, end, _, strand, _, raw_attributes in rows:
        assert 1 <= int(start) <= int(end) <= sequence_regions[record_id]
        assert strand in {"+", "-"}
        attributes = {
            key: unquote(value)
            for key, value in (item.split("=", 1) for item in raw_attributes.split(";"))
        }
        assert attributes["ID"]
        if feature_type == "gene":
            gene_ids.add(attributes["ID"])
        elif feature_type == "CDS":
            parent_ids.add(attributes["Parent"])
            assert attributes["translation"]
            translations += 1
    assert parent_ids <= gene_ids

    for expected in gff_metadata["records"]:
        record_rows = [row for row in rows if row[0] == expected["id"]]
        feature_counts = Counter(row[2] for row in record_rows)
        cds_rows = [row for row in record_rows if row[2] == "CDS"]
        strand_counts = {
            "positive": sum(row[6] == "+" for row in cds_rows),
            "negative": sum(row[6] == "-" for row in cds_rows),
            "unknown": sum(row[6] not in {"+", "-"} for row in cds_rows),
        }
        assert len(record_rows) == expected["featureCount"]
        assert dict(feature_counts) == expected["featureCounts"]
        assert strand_counts == expected["cdsStrandCounts"]

    output_statistics = gff_metadata["derivation"]["outputStatistics"]
    assert len(rows) == output_statistics["featureRowCount"] == 165
    assert len(parent_ids) == output_statistics["parentedCdsCount"] == 73
    assert translations == output_statistics["translationCount"] == 73
    assert sum(row[2] == "gene" for row in rows) == output_statistics["geneCount"]
    assert sum(row[2] == "CDS" for row in rows) == output_statistics["cdsCount"]


def test_lambda_derivative_builder_is_byte_reproducible() -> None:
    subprocess.run(
        [sys.executable, "tools/build_lambda_gff3_fixture.py", "--check"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_extended_track_fixtures_are_real_complete_records_and_reproducible() -> None:
    manifest = _load_manifest()
    files = manifest["files"]

    tobacco = manifest["fixtures"]["tobacco-plastome-regions"]
    tobacco_record = SeqIO.read(
        FIXTURE_ROOT / files["tobacco-plastome-genbank"]["relativePath"],
        "genbank",
    )
    assert (tobacco_record.id, len(tobacco_record)) == ("NC_001879.2", 155943)
    assert tobacco_record.annotations.get("topology") == "circular"
    assert "complete genome" in tobacco_record.description
    region_path = FIXTURE_ROOT / files["tobacco-plastome-regions-table"]["relativePath"]
    region_rows = [
        line.split("\t")
        for line in region_path.read_text(encoding="utf-8").splitlines()[1:]
    ]
    assert [row[1] for row in region_rows] == tobacco["expectedSemantics"]["annotationIds"]
    assert [[int(row[4]), int(row[5])] for row in region_rows] == tobacco[
        "expectedSemantics"
    ]["annotationRanges"]
    assert all(row[3] == tobacco_record.id for row in region_rows)

    depth = manifest["fixtures"]["depth-1kb"]
    depth_record = SeqIO.read(
        FIXTURE_ROOT / files["ap027133-genbank"]["relativePath"],
        "genbank",
    )
    assert (depth_record.id, len(depth_record)) == ("AP027133.1", 606194)
    assert depth_record.annotations.get("topology") == "circular"
    assert "complete genome" in depth_record.description
    depth_path = FIXTURE_ROOT / files["ap027133-drr394922-depth-1kb"]["relativePath"]
    depth_rows = [
        line.split("\t")
        for line in depth_path.read_text(encoding="utf-8").splitlines()[1:]
    ]
    expected = depth["expectedSemantics"]
    assert len(depth_rows) == expected["binCount"] == 607
    assert {row[0] for row in depth_rows} == {depth_record.id}
    assert [int(row[1]) for row in depth_rows] == list(range(1, len(depth_record) + 1, 1000))
    assert min(float(row[2]) for row in depth_rows) == 12.446
    assert max(float(row[2]) for row in depth_rows) == 74.546
    subprocess.run(
        [sys.executable, "tools/build_depth_1kb_fixture.py", "--check"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_artificial_lambda_split_is_absent_from_public_sources() -> None:
    banned_terms = (
        "lambda_" + "two_contigs",
        "lambda_" + "left",
        "lambda_" + "right",
        "gff3_" + "lambda",
        "build_lambda_" + "tutorial_fixture",
        "two-" + "contig lambda",
        "lambda two-" + "contig",
    )
    searchable_suffixes = {
        ".html",
        ".in",
        ".js",
        ".json",
        ".md",
        ".py",
        ".sh",
        ".toml",
        ".yaml",
        ".yml",
    }
    search_roots = (
        REPO_ROOT / "docs",
        REPO_ROOT / "examples",
        REPO_ROOT / "gbdraw",
        REPO_ROOT / "tools",
    )
    offenders: list[str] = []
    for root in search_roots:
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in searchable_suffixes:
                continue
            if "vendor" in path.parts:
                continue
            source = path.read_text(encoding="utf-8").lower()
            matches = [term for term in banned_terms if term in source]
            if matches:
                offenders.append(f"{path.relative_to(REPO_ROOT)}: {', '.join(matches)}")
    assert offenders == []

    former_example_dir = REPO_ROOT / "examples" / ("gff3_" + "lambda")
    former_builder = REPO_ROOT / "tools" / ("build_lambda_" + "tutorial_fixture.py")
    assert not former_example_dir.exists()
    assert not former_builder.exists()
