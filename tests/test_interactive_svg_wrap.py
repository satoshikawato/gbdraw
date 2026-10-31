"""A wrapped record's interactive SVG: one catalog item per feature, its fragments all on that one identity."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
from xml.etree import ElementTree

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

SVG = "{http://www.w3.org/2000/svg}"


def _interactive(tmp_path: Path, *args: str) -> ElementTree.Element:
    record = SeqRecord(Seq("ATGC" * 7500), id="CUT1", name="CUT1", description="wrap fixture")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(SeqFeature(SimpleLocation(8000, 14000, strand=1), type="CDS", qualifiers={"gene": ["cut"]}))
    record.features.append(SeqFeature(SimpleLocation(21000, 22000, strand=1), type="CDS", qualifiers={"gene": ["whole"]}))
    source = tmp_path / "cut.gb"
    SeqIO.write(record, source, "genbank")
    result = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--gbk", str(source),
         "--show_labels", "all", "-f", "interactive_svg", "-o", "out", *args],
        cwd=tmp_path, capture_output=True, text=True,
        env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert result.returncode == 0, result.stderr
    return ElementTree.parse(tmp_path / "out.interactive.svg").getroot()


def _catalog(root: ElementTree.Element) -> list[dict]:
    metadata = root.find(f"{SVG}metadata[@id='gbdraw-interactive-feature-metadata']")
    assert metadata is not None and metadata.text
    return json.loads(metadata.text)["items"][0]["features"]


def test_catalog_has_one_item_per_feature_and_the_fragments_share_it(tmp_path: Path) -> None:
    root = _interactive(tmp_path, "--wrap_length", "10kb")
    catalog = _catalog(root)
    assert sorted(item["drawnSelector"]["location"] for item in catalog) == ["21000..22000", "8000..14000"]
    assert len({item["svgId"] for item in catalog}) == 2

    parents = {child: parent for parent in root.iter() for child in parent}
    for item, rows in zip(sorted(catalog, key=lambda item: item["drawnSelector"]["location"]), ({"2"}, {"0", "1"})):
        parts = [
            element
            for element in root.iter(f"{SVG}path")
            if element.attrib.get("data-gbdraw-feature-id") == item["svgId"]
        ]
        assert {parents[part].attrib["data-gbdraw-wrap-row"] for part in parts} == rows
        ids = [part.attrib["id"] for part in parts]
        assert len(ids) == len(set(ids))


def test_the_catalog_is_the_unwrapped_records(tmp_path: Path) -> None:
    wrapped = _catalog(_interactive(tmp_path, "--wrap_length", "10kb"))
    plain = _catalog(_interactive(tmp_path, "--overwrite"))
    assert wrapped == plain
