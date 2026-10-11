"""The Python side of the Web cells of the CLI Session matrix (tests/utils/cli_session_cross_surface.py).

prepare OUTPUT: write every case's CLI figure and its CLI, CLI re-save, and
Python re-save Sessions; print them as JSON.
check CHECKS: draw each listed Session with render_session ("python") or a CLI
replay ("cli") and compare it with the case's CLI figure; print the results.
"""
import json
import sys
from pathlib import Path

# Direct script execution puts helpers/, not the repository, on sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tests.utils.cli_session_cross_surface import check, load_cases, render_python, replay_cli, write_case

command, argument = sys.argv[1:3]
if command == "prepare":
    output = []
    for case in load_cases():
        files = write_case(case, Path(argument) / case["id"])
        output.append({**case, **{key: str(value) for key, value in vars(files).items()}})
else:
    output = []
    for item in json.loads(Path(argument).read_text()):
        session, expected = Path(item["session"]), Path(item["expected"])
        directory = session.with_name(f"{session.name}.{item['via']}")
        directory.mkdir()
        try:
            svg = (replay_cli(item["mode"], session, directory / "replay") if item["via"] == "cli"
                   else render_python(session, directory))
        except Exception as error:  # noqa: BLE001 - report the failed cell and check the rest
            output.append({"label": item["label"], "equal": False, "message": str(error)[-600:]})
            continue
        result = check(expected, svg)
        output.append({"label": item["label"], "equal": result.equal, "message": result.message,
                       "differences": result.differences[:3]})
print(json.dumps(output))
