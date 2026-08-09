from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree

from Bio import SeqIO

from docs.recipes._scenario_support import (
    PUBLISHED_IMAGE_ROOT,
    extract_executable_block,
    load_chapter,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA


REPO_ROOT = Path(__file__).resolve().parents[1]
SCENARIO_MANIFEST = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
FIXTURE_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
FIXTURE_MANIFEST = FIXTURE_ROOT / "manifest.json"
RUNNER = "docs/recipes/run_python_scenarios.py"
SCENARIO_IDS = ("H-PY-01", "H-PY-02", "H-PY-03", "H-PY-04", "H-PY-05")
FIXTURE_FILES = {
    "H-PY-01": (
        "human-mitochondrion-genbank",
        "danio-mitochondrion-genbank",
        "drosophila-mitochondrion-genbank",
        "caenorhabditis-mitochondrion-genbank",
        "shared-cds-gene-qualifier-priority",
    ),
    "H-PY-02": (
        "lambda-genbank",
        "de3-genbank",
        "lambda-de3-losatn",
    ),
    "H-PY-03": (
        "human-mitochondrion-genbank",
        "ap027133-genbank",
        "ap027133-drr394922-depth-1kb",
        "tobacco-plastome-genbank",
        "tobacco-plastome-regions-table",
        "tobacco-plastome-default-colors",
        "tobacco-plastome-qualifier-priority",
    ),
    "H-PY-04": (
        "lambda-full-record-gff3",
        "lambda-full-record-fasta",
        "human-mitochondrion-genbank",
    ),
    "H-PY-05": ("human-mitochondrion-genbank",),
}


def _chapters() -> dict[str, dict[str, object]]:
    manifest = json.loads(SCENARIO_MANIFEST.read_text(encoding="utf-8"))
    return {
        chapter["id"]: chapter
        for chapter in manifest["chapters"]
        if chapter["id"] in SCENARIO_IDS
    }


def test_python_howto_pages_match_their_approved_scenarios() -> None:
    chapters = _chapters()

    assert set(chapters) == set(SCENARIO_IDS)
    for scenario_id, chapter in chapters.items():
        destination = REPO_ROOT / str(chapter["destination"])
        source = destination.read_text(encoding="utf-8")

        assert chapter["role"] == "how-to"
        assert chapter["surface"] == "python"
        assert chapter["execution"]["kind"] == "python-recipe"
        assert chapter["execution"]["path"] == RUNNER
        assert chapter["status"] == {
            "implementation": "verified",
            "review": "approved",
        }
        assert source.startswith("[Documentation home]")
        assert f"# {chapter['title']}" in source
        assert "## Prerequisites" in source
        assert "## Verification" in source
        assert "## Troubleshooting" in source
        assert "../../GETTING_TUTORIAL_DATA.md" in source
        assert "../../REFERENCE/" in source
        assert "tests/test_inputs" not in source
        assert "lambda_two_contigs" not in source
        assert "lambda_left" not in source
        assert "lambda_right" not in source
        assert "http://" not in source
        assert "https://" in source

        for output_name in chapter["execution"]["expected_outputs"]:
            assert output_name in source
            assert (PUBLISHED_IMAGE_ROOT / scenario_id.lower() / output_name).is_file()

    multi_record_source = (
        REPO_ROOT / str(chapters["H-PY-01"]["destination"])
    ).read_text(encoding="utf-8")
    assert "four complete, naturally circular" in " ".join(
        multi_record_source.split()
    )
    assert "linear region" not in multi_record_source
    assert "BGC0000708" not in multi_record_source
    assert "Gallery" not in multi_record_source


def test_python_howtos_use_authoritative_sequences_and_pinned_support_files() -> None:
    fixture_manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
    chapters = _chapters()

    for scenario_id, file_ids in FIXTURE_FILES.items():
        source = (REPO_ROOT / str(chapters[scenario_id]["destination"])).read_text(
            encoding="utf-8"
        )
        for file_id in file_ids:
            entry = fixture_manifest["files"][file_id]
            fixture_path = FIXTURE_ROOT / entry["relativePath"]
            payload = fixture_path.read_bytes()

            if entry["inputType"] in {"genbank", "fasta"}:
                assert entry["relativePath"] not in source
                provenance = entry["provenance"]
                authoritative_urls = [provenance["sourceUrl"]]
                verification = provenance.get("mirrorVerification")
                if verification:
                    authoritative_urls.append(
                        verification["authoritativeRequestUrl"]
                    )
                assert any(url in source for url in authoritative_urls)
                for record in entry.get("records", []):
                    assert record["id"] in source
            else:
                assert entry["relativePath"] in source
                assert entry["sha256"] in source
            assert len(payload) == entry["sizeBytes"]
            assert hashlib.sha256(payload).hexdigest() == entry["sha256"]


def test_every_python_howto_marker_has_one_public_owner_and_compiles() -> None:
    public_sources = [
        path.read_text(encoding="utf-8")
        for path in (REPO_ROOT / "docs").rglob("*.md")
        if "internal" not in path.relative_to(REPO_ROOT / "docs").parts
    ]

    for scenario_id in SCENARIO_IDS:
        start = f"<!-- executable:{scenario_id}:start -->"
        end = f"<!-- executable:{scenario_id}:end -->"
        assert sum(source.count(start) for source in public_sources) == 1
        assert sum(source.count(end) for source in public_sources) == 1

        chapter = load_chapter(
            scenario_id,
            expected_kind="python-recipe",
            runner_path=RUNNER,
        )
        recipe = extract_executable_block(chapter, language="python")
        compile(recipe, str(chapter["destination"]), "exec")
        for output_name in chapter["execution"]["expected_outputs"]:
            assert output_name in recipe

    circular_recipe = extract_executable_block(
        load_chapter("H-PY-01", expected_kind="python-recipe", runner_path=RUNNER),
        language="python",
    )
    assert "draw_circular(circular_record, options=circular_options)" in circular_recipe
    assert "draw_circular(multi_records, options=multi_record_options)" in circular_recipe
    for token in (
        "CircularOptions(",
        "LabelOptions(",
        "TitleOptions(",
        'qualifier_priority=Path("cds_gene_qualifier_priority.tsv")',
        '"labels.circular.scope": "outer"',
        '"labels.circular.placement": "horizontal"',
        '"labels.font_size.short": 18',
        '"labels.font_size.long": 18',
        '[("CDS", "gene", ".+")]',
        'text="Complete metazoan mitochondrial genomes"',
        "keep_full_definition_with_title=True",
        'Path("NC_002333.2.gb")',
        'Path("NC_024511.2.gb")',
        'Path("NC_001328.1.gb")',
        'record.annotations.get("topology") == "circular"',
    ):
        assert token in circular_recipe

    comparison_recipe = extract_executable_block(
        load_chapter("H-PY-02", expected_kind="python-recipe", runner_path=RUNNER),
        language="python",
    )
    for token in (
        "LinearComparisonOptions",
        'blast_files=("lambda-de3.losatn.tsv",)',
        'Thresholds(evalue=1e-5)',
        'Path("NC_001416.gb")',
        'Path("NC_042057.1.gb")',
    ):
        assert token in comparison_recipe
    assert "BGC0000708" not in comparison_recipe

    tracks_recipe = extract_executable_block(
        load_chapter("H-PY-03", expected_kind="python-recipe", runner_path=RUNNER),
        language="python",
    )
    for token in (
        "AnnotationOptions(table_file=\"nicotiana-tabacum-regions.tsv\")",
        "CircularTrackSlot(",
        'id="plastome_regions"',
        'renderer="annotations"',
        'id="depth"',
        'renderer="depth"',
        'default_colors=Path("modified_default_colors.tsv")',
        'qualifier_priority=Path("qualifier_priority.tsv")',
        'source=(Path("AP027133.DRR394922.depth-1kb.tsv"), None, None)',
        "depth_window=1",
        "depth_step=1000",
    ):
        assert token in tracks_recipe
    assert "record.seq[" not in tracks_recipe
    assert "record[" not in tracks_recipe

    source_recipe = extract_executable_block(
        load_chapter("H-PY-04", expected_kind="python-recipe", runner_path=RUNNER),
        language="python",
    )
    assert "read_gff(" in source_recipe
    assert "SeqIO.read(StringIO(genbank_text), \"genbank\")" in source_recipe
    assert "48_502" in source_recipe
    assert "to_svg()" in source_recipe
    assert 'to_bytes("svg")' in source_recipe

    typed_recipe = extract_executable_block(
        load_chapter("H-PY-05", expected_kind="python-recipe", runner_path=RUNNER),
        language="python",
    )
    for token in (
        "from gbdraw.api import (",
        "CircularDiagramRequest",
        "plan_request(request)",
        "build_request_plan_diagram(plan)",
        "render_prepared_request(prepared)",
        "save_session_document(",
        "materialize_session(",
        "render_session(materialized)",
        "SessionConversionError",
        '"linearTrackAxisIndex"',
    ):
        assert token in typed_recipe
    assert "gbdraw.session" not in typed_recipe
    assert "gbdraw.session_io" not in typed_recipe


def test_python_howto_recipes_regenerate_from_a_clean_external_context(
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
        [sys.executable, str(REPO_ROOT / RUNNER), "--all", "--check"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    for scenario_id in ("T-PY-01", *SCENARIO_IDS):
        assert f"{scenario_id}: verified" in result.stdout
    assert list(tmp_path.iterdir()) == []


def test_typed_session_artifact_is_current_and_deterministic() -> None:
    artifact_root = PUBLISHED_IMAGE_ROOT / "h-py-05"
    session_path = artifact_root / "typed_request.session.json"
    svg_path = artifact_root / "typed_request.svg"
    payload = json.loads(session_path.read_text(encoding="utf-8"))

    assert payload["format"] == "gbdraw-session"
    assert payload["version"] == CURRENT_SESSION_VERSION
    assert payload["createdAt"] == "2026-08-03T00:00:00+00:00"
    assert payload["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert payload["renderRequest"]["mode"] == "circular"
    assert payload["renderRequest"]["grouping"] == "single"
    assert len(payload["resources"]) == 1
    assert b"<svg" in svg_path.read_bytes()[:512]


def test_python_circular_artifacts_use_every_cds_gene_label() -> None:
    fixture_paths = (
        FIXTURE_ROOT / "human-mitochondrion" / "HmmtDNA.gbk",
        FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_002333.2.gb",
        FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_024511.2.gb",
        FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_001328.1.gb",
    )
    records = tuple(SeqIO.read(path, "genbank") for path in fixture_paths)
    artifact_root = PUBLISHED_IMAGE_ROOT / "h-py-01"

    for svg_name, artifact_records in (
        ("python_circular.svg", records[:1]),
        ("python_multi_record.svg", records),
    ):
        root = ElementTree.parse(artifact_root / svg_name).getroot()
        record_containers = {
            element.attrib["data-gbdraw-record-id"]: element
            for element in root
            if element.attrib.get("id", "").startswith("record_")
            and "data-gbdraw-record-id" in element.attrib
        }
        if len(artifact_records) == 1:
            record_containers = {artifact_records[0].id: root}
        assert set(record_containers) == {
            record.id for record in artifact_records
        }

        for record in artifact_records:
            cds_features = [
                feature for feature in record.features if feature.type == "CDS"
            ]
            expected_genes = {
                str(value)
                for feature in cds_features
                for value in feature.qualifiers.get("gene", ())
            }
            product_labels = {
                str(value)
                for feature in cds_features
                for value in feature.qualifiers.get("product", ())
            }
            assert len(expected_genes) == len(cds_features)

            label_elements = {
                "".join(element.itertext()).strip(): element
                for element in record_containers[record.id].iter()
                if element.tag.rsplit("}", 1)[-1] == "text"
            }
            assert expected_genes <= set(label_elements)
            assert not product_labels & set(label_elements)
            for gene in expected_genes:
                assert any(
                    styled.attrib.get("font-size") == "18.0"
                    for styled in label_elements[gene].iter()
                )


def test_python_tracks_artifact_retains_depth_annotation_and_label_semantics() -> None:
    svg_path = PUBLISHED_IMAGE_ROOT / "h-py-03" / "python_tracks_annotations.svg"
    root = ElementTree.parse(svg_path).getroot()
    elements = list(root.iter())

    record_ids = {
        element.attrib["data-gbdraw-record-id"]
        for element in elements
        if "data-gbdraw-record-id" in element.attrib
    }
    assert record_ids == {"AP027133.1", "NC_001879.2", "NC_012920.1"}

    annotations = {
        (
            element.attrib["data-gbdraw-annotation-id"],
            element.attrib.get("data-gbdraw-record-id"),
            element.attrib.get("data-gbdraw-annotation-label"),
        )
        for element in elements
        if "data-gbdraw-annotation-id" in element.attrib
    }
    assert annotations == {
        ("lsc", "NC_001879.2", "LSC"),
        ("irb", "NC_001879.2", "IRb"),
        ("ssc", "NC_001879.2", "SSC"),
        ("ira", "NC_001879.2", "IRa"),
    }

    slot_renderers = {
        element.attrib.get("data-gbdraw-slot-id"): element.attrib.get(
            "data-gbdraw-slot-renderer"
        )
        for element in elements
        if "data-gbdraw-slot-id" in element.attrib
    }
    assert slot_renderers["plastome_regions"] == "annotations"
    assert slot_renderers["depth"] == "depth"
    assert slot_renderers["gc_content"] == "dinucleotide_content"

    text_nodes = {
        "".join(element.itertext()).strip()
        for element in elements
        if element.tag.rsplit("}", 1)[-1] == "text"
    }
    assert {
        "rpoB",
        "secA",
        "polC",
        "psaA",
        "atpB",
        "rbcL",
        "ndhF",
        "ND1",
        "COX1",
        "ATP6",
        "CYTB",
    } <= text_nodes
    assert "cytochrome c oxidase subunit I" not in text_nodes

    fills = {element.attrib.get("fill") for element in elements}
    assert {"#d3d3d3", "#009e73", "#e69f00", "#2563eb"} <= fills
