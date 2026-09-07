"""Execute downloaded Run Info files in fresh directories, retaining real SVGs."""
import json
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from Bio import SeqIO
from gbdraw.api import load_session_document, materialize_session, session_to_request
from gbdraw.api.request_render import plan_request
from tests.utils.svg_compare import compare_svgs

archive, source, metadata, output = map(Path, sys.argv[1:5])
expected = json.loads(sys.argv[5]) if len(sys.argv) > 5 else {"start": 71, "placements": 1, "tolerance": 1}
info = json.loads(metadata.read_text())
reports = []
for kind in ("sourceRecipe", "exactReplay"):
    work = output / kind
    work.mkdir(parents=True)
    with zipfile.ZipFile(archive) as bundle:
        assert all(Path(name).name == name for name in bundle.namelist())
        bundle.extractall(work)
    (work / "shared.gbk").write_bytes(source.read_bytes())
    command = shlex.split(info[kind]["command"])
    assert command.pop(0) == "gbdraw"
    command = [sys.executable, "-c", "from gbdraw.cli import main; main()", *command,
               "--overwrite", "--session_output", "replayed.gbdraw-session.json"]
    result = subprocess.run(command, cwd=work, text=True, capture_output=True)
    (work / "stdout.log").write_text(result.stdout)
    (work / "stderr.log").write_text(result.stderr)
    report = {"kind": kind, "command": command, "exit_code": result.returncode}
    if result.returncode == 0:
        document = json.loads((work / "replayed.gbdraw-session.json").read_text())
        request = document["renderRequest"]
        svg = next(work.glob("*.svg"))
        tree = ET.parse(svg)
        paths = [node for node in tree.iter() if node.tag.endswith("}path") and node.get("d")]
        assert len(paths) >= 2
        assert svg.stat().st_size > 1000
        assert request["records"][0]["display"]["startCoordinate"] == expected["start"]
        assert len(request["diagramOptions"]["featurePlacements"]) == expected["placements"]
        options = request["diagramOptions"]
        tolerance = options.get("configOverrides", {}).get("canvas.feature_overlap_tolerance_bp",
            options.get("config", {}).get("canvas", {}).get("feature_overlap_tolerance_bp", 0))
        assert tolerance == expected["tolerance"]
        source_record = SeqIO.read(source, "genbank")
        with materialize_session(load_session_document(work / "replayed.gbdraw-session.json"),
                                 output_directory=work / "semantic-check") as materialized:
            plan = plan_request(session_to_request(materialized))
            assert len(plan.records) == 1
            assert str(plan.records[0].seq) == str(source_record.seq)
            assert [str(feature.location) for feature in plan.records[0].features] == [
                str(feature.location) for feature in source_record.features]
            assert plan.transforms[0].source_base_to_display_index(expected["start"] or 1) == 0
        for placement in options["featurePlacements"]:
            feature_id = placement["biologicalFeatureId"]
            assert any(node.get("data-gbdraw-feature-id", "").endswith(feature_id)
                       or node.get("data-gbdraw-stable-feature-id") == feature_id for node in paths)
        report.update(svg=str(svg), path_count=len(paths), request=request)
    reports.append(report)
(output / "replay-report.json").write_text(json.dumps(reports, indent=2) + "\n")
for report in reports:
    assert report["exit_code"] == 0, report
assert compare_svgs(Path(reports[0]["svg"]), Path(reports[1]["svg"])).equal
