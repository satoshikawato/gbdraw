from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree

from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_ID = "H-CLI-01"
LAYOUT_SCENARIO_ID = "H-CLI-03"
LINEAR_LAYOUT_SCENARIO_ID = "H-CLI-04"
RUNNER = "docs/recipes/run_cli_scenarios.py"


def test_input_how_to_uses_whole_lambda_records_on_both_input_paths() -> None:
    chapter = load_chapter(
        SCENARIO_ID,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    source = (REPO_ROOT / chapter["destination"]).read_text(encoding="utf-8")
    recipe = extract_executable_block(chapter, language="bash")

    assert chapter["execution"]["expected_outputs"] == [
        "lambda_genbank.svg",
        "lambda_gff3.svg",
    ]
    assert recipe.count("gbdraw linear") == 2
    assert "--gbk NC_001416.gb" in recipe
    assert "--gff NC_001416.gff3" in recipe
    assert "--fasta NC_001416.fna" in recipe
    assert "48,502 bp" in source
    assert "without cropping or\nsplitting it" in source
    assert "artificial contigs" in source
    assert "lambda_two_contigs" not in source


def test_input_how_to_fixtures_are_complete_and_sequence_identical() -> None:
    fixture_manifest = json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )
    entries = fixture_manifest["files"]
    records = {
        file_id: entries[file_id]["records"][0]
        for file_id in (
            "lambda-genbank",
            "lambda-full-record-gff3",
            "lambda-full-record-fasta",
        )
    }

    assert {record["id"] for record in records.values()} == {"NC_001416.1"}
    assert {record["length"] for record in records.values()} == {48_502}
    assert records["lambda-genbank"]["sequenceSha256"] == records[
        "lambda-full-record-fasta"
    ]["sequenceSha256"]
    assert records["lambda-full-record-gff3"]["sourceRegion"] == (
        "NC_001416.1:1-48502"
    )
    assert records["lambda-full-record-gff3"]["featureCounts"]["CDS"] == 73
    assert records["lambda-full-record-gff3"]["cdsStrandCounts"] == {
        "positive": 47,
        "negative": 26,
        "unknown": 0,
    }


def test_input_how_to_regenerates_from_a_clean_external_directory(
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not existing_pythonpath
        else os.pathsep.join((str(REPO_ROOT), existing_pythonpath))
    )
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            SCENARIO_ID,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.count("H-CLI-01: verified") == 2
    assert list(tmp_path.iterdir()) == []
    for output_name in ("lambda_genbank.svg", "lambda_gff3.svg"):
        assert (PUBLISHED_IMAGE_ROOT / "h-cli-01" / output_name).is_file()


def test_multi_record_layout_uses_four_complete_circular_mitochondria() -> None:
    chapter = load_chapter(
        LAYOUT_SCENARIO_ID,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    source = (REPO_ROOT / chapter["destination"]).read_text(encoding="utf-8")
    recipe = extract_executable_block(chapter, language="bash")

    assert chapter["fixtures"] == [
        "human-mitochondrion",
        "metazoan-mitochondria-four",
    ]
    assert (
        "--gbk HmmtDNA.gbk NC_002333.2.gb NC_024511.2.gb NC_001328.1.gb"
        in recipe
    )
    assert "--multi_record_canvas" in recipe
    assert "--multi_record_size_mode equal" in recipe
    assert "--multi_record_position '#1@1'" in recipe
    assert "--multi_record_position '#2@1'" in recipe
    assert "--multi_record_position '#3@2'" in recipe
    assert "--multi_record_position '#4@2'" in recipe
    assert "--circular_track_order ticks,features,gc_content,gc_skew" in recipe
    assert "--circular_track_axis_index 1" in recipe
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in recipe
    assert "--labels out" in recipe
    assert "--label_font_size 10" in recipe
    assert "complete, naturally circular RefSeq mitochondrial records" in source
    assert "No record is cropped, concatenated, or split" in source
    assert "147 displayed features" in source
    assert "Each CDS label uses the feature's concise `gene` value" in source
    assert "`product`\ndescriptions are not used as labels" in source
    for organism in (
        "human",
        "zebrafish",
        "fruit fly",
        "nematode",
    ):
        assert organism in source

    fixture_manifest = json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )
    multi_record = fixture_manifest["fixtures"]["metazoan-mitochondria-four"]
    shared_labels = fixture_manifest["fixtures"]["shared-label-rules"]
    assert "H-CLI-03" in shared_labels["scenarioIds"]
    assert multi_record["fileIds"] == [
        "danio-mitochondrion-genbank",
        "drosophila-mitochondrion-genbank",
        "caenorhabditis-mitochondrion-genbank",
    ]
    assert multi_record["fileReferences"] == ["human-mitochondrion-genbank"]
    assert multi_record["expectedSemantics"]["recordsAreWholeCanonicalSources"]
    assert not multi_record["expectedSemantics"]["comparisonEvidence"]
    assert multi_record["expectedSemantics"]["recordIds"] == [
        "NC_012920.1",
        "NC_002333.2",
        "NC_024511.2",
        "NC_001328.1",
    ]
    assert multi_record["expectedSemantics"]["recordLengths"] == [
        16_569,
        16_596,
        19_524,
        13_794,
    ]
    assert multi_record["expectedSemantics"]["topologies"] == ["circular"] * 4

    file_ids = [
        "human-mitochondrion-genbank",
        *multi_record["fileIds"],
    ]
    expected_features = [37, 37, 37, 36]
    for file_id, record_id, length, feature_count in zip(
        file_ids,
        multi_record["expectedSemantics"]["recordIds"],
        multi_record["expectedSemantics"]["recordLengths"],
        expected_features,
        strict=True,
    ):
        entry = fixture_manifest["files"][file_id]
        record = entry["records"][0]
        source_path = REPO_ROOT / "gbdraw/web/tutorial-data" / entry["relativePath"]
        assert entry["role"] == "raw"
        assert entry["derivation"] is None
        assert record["id"] == record_id
        assert record["length"] == length
        assert record["topology"] == "circular"
        assert record["description"].endswith("complete genome")
        assert record["displayedFeatureCount"] == feature_count
        assert hashlib.sha256(source_path.read_bytes()).hexdigest() == entry["sha256"]


def test_multi_record_layout_regenerates_from_a_clean_external_directory(
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (str(REPO_ROOT), environment.get("PYTHONPATH"))
        if value
    )
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            LAYOUT_SCENARIO_ID,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "H-CLI-03: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []
    assert (
        PUBLISHED_IMAGE_ROOT
        / "h-cli-03"
        / "multi_record_circular_cli.svg"
    ).is_file()


def test_linear_layout_how_to_uses_complete_sources_and_explicit_positions() -> None:
    chapter = load_chapter(
        LINEAR_LAYOUT_SCENARIO_ID,
        expected_kind="cli-recipe",
        runner_path=RUNNER,
    )
    source = (REPO_ROOT / chapter["destination"]).read_text(encoding="utf-8")
    recipe = extract_executable_block(chapter, language="bash")
    prose = " ".join(source.split())

    assert chapter["fixtures"] == ["lambda", "aminoglycoside-bgc-five"]
    assert chapter["execution"]["expected_outputs"] == ["linear_layout_cli.svg"]
    assert "# How to arrange linear records, regions, orientation, labels, and rulers" in source
    assert "## Prerequisites" in source
    assert "## Verification" in source
    assert "## Troubleshooting" in source
    assert "../../../gbdraw/web/tutorial-data/manifest.json" in source
    assert "../../images/h-cli-04/linear_layout_cli.svg" in source
    assert "--gbk NC_001416.gb BGC0000708.gbk BGC0000713.gbk" in recipe
    assert "--region NC_001416.1:5001-35500" in recipe
    assert recipe.count("--record_id") == 3
    assert recipe.count("--reverse_complement") == 3
    assert "--reverse_complement false \\\n  --reverse_complement false \\\n  --reverse_complement true" in recipe
    assert "--multi_record_position '#1@1'" in recipe
    assert "--multi_record_position '#2@2'" in recipe
    assert "--multi_record_position '#3@2'" in recipe
    assert "--qualifier_priority cds_gene_qualifier_priority.tsv" in recipe
    assert "--show_labels all" in recipe
    assert "--ruler_on_axis" in recipe
    assert "-o linear_layout_cli" in recipe
    assert "local ruler still reads from 5 kbp through 30 kbp from left to right" in prose
    assert "`racG` appears to the left of `racP`" in prose
    assert "descend" not in source.lower()


def test_linear_layout_sources_match_canonical_fixture_records() -> None:
    manifest = json.loads(
        (REPO_ROOT / "gbdraw/web/tutorial-data/manifest.json").read_text(
            encoding="utf-8"
        )
    )
    expected = {
        "lambda-genbank": (
            "NC_001416.1",
            48_502,
            73,
            "4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7",
        ),
        "bgc-0000708-genbank": (
            "BGC0000708",
            40_579,
            30,
            "9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0",
        ),
        "bgc-0000713-genbank": (
            "BGC0000713",
            31_892,
            26,
            "bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb",
        ),
    }

    for file_id, (record_id, length, cds_count, checksum) in expected.items():
        entry = manifest["files"][file_id]
        record = entry["records"][0]
        source_path = REPO_ROOT / "gbdraw/web/tutorial-data" / entry["relativePath"]
        assert entry["role"] == "raw"
        assert entry["derivation"] is None
        assert record["id"] == record_id
        assert record["length"] == length
        assert record["displayedFeatureCount"] == cds_count
        assert entry["sha256"] == checksum
        assert hashlib.sha256(source_path.read_bytes()).hexdigest() == checksum

    qualifier = manifest["files"]["shared-cds-gene-qualifier-priority"]
    qualifier_path = (
        REPO_ROOT / "gbdraw/web/tutorial-data" / qualifier["relativePath"]
    )
    assert qualifier["sha256"] == (
        "1d3c787c5768b52191cc7e8a325cd00fdc4b1e9b495d37e580c3a69dec21b87a"
    )
    assert hashlib.sha256(qualifier_path.read_bytes()).hexdigest() == qualifier["sha256"]

    bgc_fixture = manifest["fixtures"]["aminoglycoside-bgc-five"]
    assert bgc_fixture["expectedSemantics"]["recordsAreWholeCanonicalSources"]
    assert {"bgc-0000708-genbank", "bgc-0000713-genbank"} <= set(
        bgc_fixture["fileIds"]
    )


def test_linear_layout_regenerates_from_a_clean_external_directory(
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        value
        for value in (str(REPO_ROOT), environment.get("PYTHONPATH"))
        if value
    )
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / RUNNER),
            "--scenario",
            LINEAR_LAYOUT_SCENARIO_ID,
            "--check",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.count("H-CLI-04: verified") == 1
    assert list(tmp_path.iterdir()) == []
    assert (PUBLISHED_IMAGE_ROOT / "h-cli-04" / "linear_layout_cli.svg").is_file()


def test_linear_layout_svg_keeps_rulers_ascending_and_reverses_features() -> None:
    output_path = PUBLISHED_IMAGE_ROOT / "h-cli-04" / "linear_layout_cli.svg"
    root = ElementTree.parse(output_path).getroot()
    record_groups = {
        element.attrib["data-gbdraw-record-id"]: element
        for element in root
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
        and "data-gbdraw-record-id" in element.attrib
    }
    expected_ticks = [
        "5 kbp",
        "10 kbp",
        "15 kbp",
        "20 kbp",
        "25 kbp",
        "30 kbp",
    ]
    for record_id in ("NC_001416.1", "BGC0000713"):
        ticks = [
            ("".join(element.itertext()).strip(), float(element.attrib["x"]))
            for element in record_groups[record_id].iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
            and "".join(element.itertext()).strip().endswith("kbp")
            and "x" in element.attrib
        ]
        assert [label for label, _x in ticks] == expected_ticks
        assert all(
            left_x < right_x
            for (_, left_x), (_, right_x) in zip(ticks, ticks[1:])
        )

    feature_x = {
        "".join(element.itertext()).strip(): float(element.attrib["x"])
        for element in record_groups["BGC0000713"].iter()
        if element.tag.rsplit("}", 1)[-1] == "text"
        and "".join(element.itertext()).strip() in {"racG", "racP"}
        and "x" in element.attrib
    }
    assert set(feature_x) == {"racG", "racP"}
    assert feature_x["racG"] < feature_x["racP"]
