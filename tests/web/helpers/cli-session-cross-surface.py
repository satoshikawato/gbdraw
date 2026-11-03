"""The Python side of the Web cells of the CLI Session matrix (tests/utils/cli_session_cross_surface.py).

prepare OUTPUT: write every case's CLI figure of each of its drawings (`svgs`,
by drawing ID) and its CLI, CLI re-save, and Python re-save (`python_resaves`,
by drawing ID) Sessions; print them as JSON.
check CHECKS: draw each listed Session with render_session ("python") or a CLI
replay ("cli"), the drawing an item names (`drawing`) or else the one the
Session opens, and compare it with the case's CLI figure; print the results.
"""
import json
import sys
from pathlib import Path

# Direct script execution puts helpers/, not the repository, on sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tests.utils.cli_session_cross_surface import (
    case_drawings, check, load_cases, render_python, replay_cli, write_case,
)

command, argument = sys.argv[1:3]
if command == "prepare":
    output = []
    for case in load_cases():
        files = write_case(case, Path(argument) / case["id"])
        output.append({
            **case, "drawings": case_drawings(case),
            "svgs": {drawing: str(svg) for drawing, svg in files.svgs.items()},
            "python_resaves": {drawing: str(path) for drawing, path in files.python_resaves.items()},
            "cli_session": str(files.cli_session), "cli_resave": str(files.cli_resave),
        })
else:
    output = []
    for item in json.loads(Path(argument).read_text()):
        session, expected = Path(item["session"]), Path(item["expected"])
        directory = session.with_name(f"{session.name}.{item['via']}")
        directory.mkdir()
        try:
            drawing = item.get("drawing")
            svg = (replay_cli(item["mode"], session, directory / "replay", drawing=drawing) if item["via"] == "cli"
                   else render_python(session, directory, drawing))
        except Exception as error:  # noqa: BLE001 - report the failed cell and check the rest
            output.append({"label": item["label"], "equal": False, "message": str(error)[-600:]})
            continue
        result = check(expected, svg)
        output.append({"label": item["label"], "equal": result.equal, "message": result.message,
                       "differences": result.differences[:3]})
print(json.dumps(output))
