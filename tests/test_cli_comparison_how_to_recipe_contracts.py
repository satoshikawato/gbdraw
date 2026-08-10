from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)


pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = "docs/recipes/run_cli_scenarios.py"
SCENARIOS = ("H-CLI-05", "H-CLI-06", "H-CLI-07", "H-CLI-08")
EVIDENCE_SOURCE = "docs/internal/SCENARIO_EVIDENCE.md"
BGC_RECORDS = (
    ("bgc-0000708-genbank", "BGC0000708", 40_579, 30, "9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0"),
    ("bgc-0000709-genbank", "BGC0000709", 50_466, 38, "4b66b7e4b78d429d12176e1e36d0e48178c562a9d128d4308b38753af9995255"),
    ("bgc-0000711-genbank", "BGC0000711", 30_837, 21, "32393648f6a91166444331b83687f1b9b7b24c60553a7ddcb677dfe207736789"),
    ("bgc-0000712-genbank", "BGC0000712", 48_169, 40, "705104a0daa5c44981b0a1e5352d3e56f012dd2e3ae94c98c85cd0ee9198bf94"),
    ("bgc-0000713-genbank", "BGC0000713", 31_892, 26, "bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb"),
)


def _fixture_manifest() -> dict[str, object]:
    return json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )


def test_comparison_evidence_manifest_entries_are_complete() -> None:
    expected = {
        "H-CLI-05": (
            [
                "lambda",
                "de3",
                "lambda-de3-comparison",
                "human-mitochondrion",
                "metazoan-mitochondria-four",
                "metazoan-mitochondria-comparison",
            ],
            ["linear_precomputed_comparison.svg", "circular_conservation_ring.svg"],
        ),
        "H-CLI-06": (
            ["aminoglycoside-bgc-five"],
            ["cli_losatp_pairwise.tsv", "cli_losatp_pairwise.svg"],
        ),
        "H-CLI-07": (
            ["aminoglycoside-bgc-five"],
            ["cli_losatp_groups.svg"],
        ),
        "H-CLI-08": (
            ["aminoglycoside-bgc-five"],
            ["cli_losatp_collinear.svg"],
        ),
    }

    for scenario_id, (fixtures, outputs) in expected.items():
        chapter = load_chapter(
            scenario_id,
            expected_kind="cli-recipe",
            runner_path=RUNNER,
        )
        assert chapter["role"] == "evidence"
        assert "destination" not in chapter
        assert chapter["execution"]["source"] == EVIDENCE_SOURCE
        assert chapter["fixtures"] == fixtures
        assert chapter["execution"]["expected_outputs"] == outputs
        assert chapter["status"]["implementation"] == "verified"
        assert chapter["status"]["review"] == "approved"
        assert extract_executable_block(chapter, language="bash")


def test_precomputed_recipe_uses_linear_phages_and_circular_complete_mtdna() -> None:
    chapter = load_chapter(
        "H-CLI-05",
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    recipe = extract_executable_block(chapter, language="bash")
    manifest = _fixture_manifest()
    files = manifest["files"]

    assert recipe.count("gbdraw ") == 2
    assert "--gbk NC_001416.gb NC_042057.1.gb" in recipe
    assert "--blast lambda-de3.losatn.tsv" in recipe
    assert "--gbk HmmtDNA.gbk" in recipe
    assert (
        "--conservation_blast danio-human.tlosatx.tsv "
        "drosophila-human.tlosatx.tsv caenorhabditis-human.tlosatx.tsv"
    ) in recipe
    assert "--conservation_reference subject" in recipe
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in recipe
    assert "--labels out" in recipe
    expected_files = {
        "lambda-genbank": ("NC_001416.1", 48_502, "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7"),
        "de3-genbank": ("NC_042057.1", 42_925, "288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09"),
    }
    for file_id, (record_id, length, checksum) in expected_files.items():
        entry = files[file_id]
        path = REPO_ROOT / "gbdraw/web/tutorial-data" / entry["relativePath"]
        assert entry["records"][0]["id"] == record_id
        assert entry["records"][0]["length"] == length
        assert entry["sha256"] == checksum
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum

    losatn = files["lambda-de3-losatn"]
    assert losatn["sha256"] == (
        "703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c"
    )
    assert losatn["expectedSemantics"]["queryRecordId"] == "NC_001416.1"
    assert losatn["expectedSemantics"]["subjectRecordId"] == "NC_042057.1"
    assert losatn["expectedSemantics"]["rawRowCount"] == 6
    assert losatn["expectedSemantics"]["defaultRetainedRowCount"] == 6

    comparison = manifest["fixtures"]["metazoan-mitochondria-comparison"]
    semantics = comparison["expectedSemantics"]
    assert semantics["referenceRecordId"] == "NC_012920.1"
    assert semantics["referenceSide"] == "subject"
    assert semantics["sourceTopologies"] == ["circular"] * 4
    assert semantics["recordsAreWholeCanonicalSources"] is True
    assert semantics["totalRawRows"] == 435
    assert semantics["retainedRows"] == 106
    assert semantics["retainedUnionCoverageBp"] == 9813
    expected_queries = {
        "danio-human-tlosatx": ("NC_002333.2", 276, 68),
        "drosophila-human-tlosatx": ("NC_024511.2", 93, 24),
        "caenorhabditis-human-tlosatx": ("NC_001328.1", 66, 14),
    }
    for file_id, (query_id, raw_rows, retained_rows) in expected_queries.items():
        evidence = files[file_id]["expectedSemantics"]
        assert evidence["queryRecordId"] == query_id
        assert evidence["subjectRecordId"] == "NC_012920.1"
        assert evidence["rawRowCount"] == raw_rows
        assert evidence["retainedRowCount"] == retained_rows


def test_losatp_recipes_use_five_whole_records_and_distinct_modes() -> None:
    recipes = {
        scenario_id: extract_executable_block(
            load_chapter(
                scenario_id,
                expected_kind="cli-recipe",
                runner_path=RUNNER,
            ),
            language="bash",
        )
        for scenario_id in ("H-CLI-06", "H-CLI-07", "H-CLI-08")
    }
    expected_input = (
        "--gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk "
        "BGC0000712.gbk BGC0000713.gbk"
    )
    for recipe in recipes.values():
        assert expected_input in recipe
        assert "--losatp_threads 1" in recipe
        assert "--identity 30" in recipe
    assert "--protein_blastp_mode pairwise" in recipes["H-CLI-06"]
    assert "--protein_blastp_max_hits 1" in recipes["H-CLI-06"]
    assert "--protein_blastp_output cli_losatp_pairwise.tsv" in recipes["H-CLI-06"]
    assert "--protein_blastp_mode orthogroup" in recipes["H-CLI-07"]
    assert "--align_orthogroup_feature CAG38695.1" in recipes["H-CLI-07"]
    assert "--protein_blastp_mode collinear" in recipes["H-CLI-08"]
    assert "--collinear_search_scope adjacent" in recipes["H-CLI-08"]
    assert "--collinear_min_anchors 2" in recipes["H-CLI-08"]

    manifest = _fixture_manifest()
    files = manifest["files"]
    total_features = 0
    for file_id, record_id, length, cds_count, checksum in BGC_RECORDS:
        entry = files[file_id]
        path = REPO_ROOT / "gbdraw/web/tutorial-data" / entry["relativePath"]
        record = entry["records"][0]
        assert entry["role"] == "raw"
        assert entry["derivation"] is None
        assert record["id"] == record_id
        assert record["length"] == length
        assert record["displayedFeatureCount"] == cds_count
        assert entry["sha256"] == checksum
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum
        total_features += cds_count
    assert total_features == 155
    bgc = manifest["fixtures"]["aminoglycoside-bgc-five"]
    assert bgc["expectedSemantics"]["recordCount"] == 5
    assert bgc["expectedSemantics"]["recordsAreWholeCanonicalSources"]
    assert bgc["expectedSemantics"]["galleryAlignedLosatpMode"] == "similarity-groups"


def test_technical_documentation_distinguishes_pairwise_and_similarity_groups() -> None:
    source = (
        REPO_ROOT
        / "docs/REFERENCE/comparison-programs-thresholds-and-results.md"
    ).read_text(encoding="utf-8")
    prose = " ".join(source.split())
    assert "Pairwise" in source
    assert "Similarity groups" in source
    assert "not a phylogenetic orthogroup" in prose


def test_published_pairwise_raw_evidence_is_hydrated_and_complete() -> None:
    path = PUBLISHED_IMAGE_ROOT / "h-cli-06" / "cli_losatp_pairwise.tsv"
    source = path.read_text(encoding="utf-8")
    rows = [
        line.split("\t")
        for line in source.splitlines()
        if line and not line.startswith("#")
    ]
    assert source.count("# entry ") == 4
    assert len(rows) == 791
    assert all(len(row) == 12 for row in rows)
    assert "h_" not in source


@pytest.mark.parametrize("scenario_id", SCENARIOS)
def test_comparison_evidence_regenerates_from_a_clean_external_directory(
    scenario_id: str,
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (str(REPO_ROOT), environment.get("PYTHONPATH"))
        if value
    )
    chapter = load_chapter(
        scenario_id,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            scenario_id,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=360,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.count(f"{scenario_id}: verified") == len(
        chapter["execution"]["expected_outputs"]
    )
    assert list(tmp_path.iterdir()) == []
    for output_name in chapter["execution"]["expected_outputs"]:
        assert (PUBLISHED_IMAGE_ROOT / scenario_id.lower() / output_name).is_file()
