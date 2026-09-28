"""Audit exact scientific fields around an existing Web Save/native CLI replay."""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import io
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import SeqIO

DRAWING_TAGS = {
    "svg", "g", "path", "circle", "ellipse", "rect", "line",
    "polyline", "polygon", "text", "textPath",
}
GEOMETRY_ATTRIBUTES = {
    "viewBox", "transform", "d", "cx", "cy", "r", "rx", "ry", "x", "y",
    "x1", "x2", "y1", "y2", "points", "width", "height", "font-size",
    "startOffset",
}


def read_document(path: Path) -> dict:
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)


def drawing_signature(node: ET.Element):
    tag = node.tag.rsplit("}", 1)[-1]
    children = tuple(drawing_signature(child) for child in node)
    if tag not in DRAWING_TAGS:
        return children
    text = "".join(node.itertext()) if tag in {"text", "textPath"} else None
    return tag, text, children


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("generated", "pending"), required=True)
    parser.add_argument("--web", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--svg", type=Path, required=True)
    args = parser.parse_args()

    web, native = read_document(args.web), read_document(args.sidecar)
    request, replay = web["renderRequest"], native["renderRequest"]
    assert (request["mode"], request["grouping"]) == (replay["mode"], replay["grouping"])
    assert request["diagramOptions"]["tracks"] == replay["diagramOptions"]["tracks"]
    assert web["config"] == native["config"]

    assert web["resources"].keys() == native["resources"].keys()
    resource_hashes = {}
    for key, value in web["resources"].items():
        source = base64.b64decode(value["data"])
        assert source == base64.b64decode(native["resources"][key]["data"]), key
        resource_hashes[key] = hashlib.sha256(source).hexdigest()

    record_ids = []
    for original, restored in zip(request["records"], replay["records"], strict=True):
        assert {key: value for key, value in original.items() if key != "selector"} == {
            key: value for key, value in restored.items() if key != "selector"
        }
        resource = web["resources"][original["source"]["resourceId"]]
        source = base64.b64decode(resource["data"]).decode()
        ids = [record.id for record in SeqIO.parse(io.StringIO(source), "genbank")]
        record_ids.extend(ids)
        if original["selector"] != restored["selector"]:
            assert len(ids) == 1
            assert original["selector"] == {"kind": "recordId", "value": ids[0]}
            assert restored["selector"] is None

    original_geometry = json.loads(json.dumps(web["runMetadata"]["trackSlotGeometry"]))
    native_geometry = json.loads(json.dumps(native["runMetadata"]["trackSlotGeometry"]))
    for geometry in (original_geometry, native_geometry):
        for record in geometry["records"]:
            record.pop("resultName", None)  # Explicit CLI output override.
    assert original_geometry == native_geometry

    slots = request["diagramOptions"]["tracks"]["circularTrackSlots"]
    gc_width = next(slot["width"] for slot in slots if slot["id"] == "gc_content")
    draft_gc_width = next(
        slot["width"] for slot in web["config"]["adv"]["circular_track_slots"]
        if slot["id"] == "gc_content"
    ) if args.case == "pending" else None
    if args.case == "pending":
        assert gc_width == {"value": 0.08, "unit": "factor"}
        assert draft_gc_width == {"value": "1.", "unit": "px"}
    else:
        assert next(slot["width"] for slot in slots if slot["id"] == "features") == {
            "value": 20, "unit": "px"
        }

    web_svg = ET.fromstring(web["results"][0]["content"])
    native_svg = ET.fromstring(args.svg.read_text())
    assert drawing_signature(web_svg) == drawing_signature(native_svg)
    drawing = [node for node in web_svg.iter()
               if node.tag.rsplit("}", 1)[-1] in DRAWING_TAGS]
    print(json.dumps({
        "case": args.case,
        "webSha256": hashlib.sha256(args.web.read_bytes()).hexdigest(),
        "sidecarSha256": hashlib.sha256(args.sidecar.read_bytes()).hexdigest(),
        "mode": request["mode"],
        "grouping": request["grouping"],
        "recordIds": record_ids,
        "resourceSha256": resource_hashes,
        "draftGcWidth": draft_gc_width,
        "committedGcWidth": gc_width,
        "configExact": True,
        "geometryExactExceptOutputName": True,
        "hierarchicalDrawingAndTextExact": True,
        "drawingElements": len(drawing),
        "textElements": sum(node.tag.rsplit("}", 1)[-1] in {"text", "textPath"}
                            for node in drawing),
        "transformElements": sum("transform" in node.attrib for node in drawing),
        "comparedAttributeSet": sorted({
            key for node in drawing for key in node.attrib if key in GEOMETRY_ATTRIBUTES
        }),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
