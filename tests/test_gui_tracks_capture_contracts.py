from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest
from Bio import SeqIO
from PIL import Image

from docs.capture.config import chapter_for

pytestmark = pytest.mark.recipe


REPO_ROOT = Path(__file__).resolve().parents[1]
FLOW_PATH = REPO_ROOT / "docs" / "capture" / "flows" / "how_to" / "tracks.py"
INDEX_PATH = REPO_ROOT / "gbdraw" / "web" / "index.html"
DEPTH_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "depth-1kb"
TOBACCO_ROOT = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "tobacco-plastome-regions"
)
SCREENSHOTS = {
    "H-GUI-09": (
        REPO_ROOT / "docs" / "images" / "h-gui-09" / "track-settings.png",
        REPO_ROOT / "docs" / "images" / "h-gui-09" / "track-result.png",
    ),
    "H-GUI-10": (
        REPO_ROOT / "docs" / "images" / "h-gui-10" / "slot-settings.png",
        REPO_ROOT / "docs" / "images" / "h-gui-10" / "annotation-result.png",
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_quantitative_capture_uses_one_complete_record_and_607_real_bins() -> None:
    genbank = DEPTH_ROOT / "AP027133.gb"
    depth = DEPTH_ROOT / "AP027133.DRR394922.depth-1kb.tsv"
    assert genbank.stat().st_size == 1_344_094
    assert _sha256(genbank) == (
        "913af50dd9d37cc2107be5e46484b885c5d586fb414b4b501380fc8f17a659d6"
    )
    records = list(SeqIO.parse(genbank, "genbank"))
    assert len(records) == 1
    assert records[0].id == "AP027133.1"
    assert len(records[0]) == 606_194
    assert records[0].annotations["topology"] == "circular"
    assert "complete" in records[0].description.lower()

    assert depth.stat().st_size == 16_913
    assert _sha256(depth) == (
        "6f57cfd89a165ad97a162aa2f0b1f3b3ad21fb5638f4f9ac5cbd069badd6aab7"
    )
    with depth.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 607
    assert [int(row["position"]) for row in rows] == [
        1 + (index * 1_000) for index in range(607)
    ]
    assert {row["reference_name"] for row in rows} == {"AP027133.1"}
    values = [float(row["depth"]) for row in rows]
    assert min(values) == 12.446
    assert max(values) == 74.546
    assert values[0] == 21.901
    assert values[-1] == 15.938144


def test_annotation_capture_uses_the_complete_plastome_and_four_regions() -> None:
    genbank = TOBACCO_ROOT / "NC_001879.gbk"
    table = TOBACCO_ROOT / "nicotiana-tabacum-regions.tsv"
    assert genbank.stat().st_size == 331_860
    assert _sha256(genbank) == (
        "25c5b39fd25d702c0a390fe5e7480eda0ccc1e4d6d7c388445b4686049412a24"
    )
    records = list(SeqIO.parse(genbank, "genbank"))
    assert len(records) == 1
    assert records[0].id == "NC_001879.2"
    assert len(records[0]) == 155_943
    assert records[0].annotations["topology"] == "circular"
    assert "complete" in records[0].description.lower()

    assert table.stat().st_size == 594
    assert _sha256(table) == (
        "3a85aed5145c88f93b4478d1901fab53714b9d47afc754d32cc9e5c0b8412b88"
    )
    with table.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["id"] for row in rows] == ["lsc", "irb", "ssc", "ira"]
    assert [row["label"] for row in rows] == ["LSC", "IRb", "SSC", "IRa"]
    assert [(int(row["start"]), int(row["end"])) for row in rows] == [
        (1, 86_686),
        (86_687, 112_029),
        (112_030, 130_600),
        (130_601, 155_943),
    ]
    assert {row["record"] for row in rows} == {"NC_001879.2"}
    assert {row["set_id"] for row in rows} == {"plastome_regions"}
    assert {row["mark"] for row in rows} == {"bracket"}
    assert {row["lane"] for row in rows} == {"0"}


def test_track_flows_use_visible_accessible_controls_without_state_injection() -> None:
    source = FLOW_PATH.read_text(encoding="utf-8")
    required = (
        'get_by_role("button", name="Circular", exact=True)',
        'get_by_role("radio", name="GenBank", exact=True)',
        'get_by_label("GenBank/DDBJ File", exact=True).set_input_files',
        'get_by_label("Depth TSV tracks", exact=True)',
        'get_by_label("Depth TSV", exact=True).set_input_files',
        'get_by_label("Depth Window", exact=True)',
        'get_by_label("Depth Step", exact=True)',
        'get_by_label("Depth Min", exact=True)',
        'get_by_label("Depth Max", exact=True)',
        'name="Depth Axis", exact=True',
        'name="Depth Ticks", exact=True',
        'get_by_label("Region Annotations", exact=True)',
        'get_by_label("Import TSV", exact=True).set_input_files',
        'get_by_label("Annotation lane", exact=True)',
        'name="Circular track slot gc_skew", exact=True',
        'get_by_label("Track dinucleotide", exact=True)',
        'get_by_label("New circular track renderer", exact=True)',
        'name=re.compile(r"Add track$")',
        'name="Circular track slot annotations", exact=True',
        'get_by_label("Annotation set", exact=True)',
        'get_by_label("Annotation placement", exact=True)',
        'get_by_title("Move outside Axis", exact=True)',
        '"Show annotation labels", exact=True',
        'get_by_role("button", name="SVG", exact=True)',
        "page.expect_download",
        "download.save_as",
    )
    for fragment in required:
        assert fragment in source

    for forbidden in (
        "page.locator(",
        "get_by_text(",
        "localStorage",
        "sessionStorage",
        "add_init_script",
        ".runAnalysis(",
        ".downloadSVG(",
        "set_content(",
        "region_start",
        "region_end",
        "part_",
    ):
        assert forbidden not in source


def test_track_capture_accessibility_labels_are_part_of_the_public_ui() -> None:
    source = INDEX_PATH.read_text(encoding="utf-8")
    for label in (
        'aria-label="Depth TSV tracks"',
        'aria-label="Depth Window"',
        'aria-label="Depth Step"',
        'aria-label="Depth Min"',
        'aria-label="Depth Max"',
        'aria-label="Region Annotations"',
        'aria-label="Dinucleotide content/skew"',
        'aria-label="New circular track renderer"',
        'aria-label="Annotation set"',
        'aria-label="Annotation placement"',
        'aria-label="Show annotation labels"',
        'aria-label="Track dinucleotide"',
        'aria-label="Track legend label"',
    ):
        assert label in source


def test_track_evidence_matches_the_approved_manifest() -> None:
    for scenario_id in ("H-GUI-09", "H-GUI-10"):
        chapter = chapter_for(scenario_id)
        assert chapter["role"] == "evidence"
        assert "destination" not in chapter
        assert chapter["execution"]["path"] == (
            "docs/capture/flows/how_to/tracks.py"
        )
        assert tuple(
            REPO_ROOT / screenshot["path"] for screenshot in chapter["screenshots"]
        ) == SCREENSHOTS[scenario_id]

def test_committed_track_screenshots_are_full_pinned_viewports() -> None:
    for paths in SCREENSHOTS.values():
        for path in paths:
            assert path.is_file(), path
            with Image.open(path) as image:
                assert image.size == (1440, 900)
                assert image.mode == "RGB"
