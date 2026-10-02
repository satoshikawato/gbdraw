"""Run a Run Info Source recipe in a fresh directory and compare its SVG.

Usage: source-recipe-replay.py WORK COMMAND GUI_SVG PREVIOUS_SVG BUNDLE_OR_DASH INPUT...

The CLI SVG must equal the GUI Result (GUI_SVG). It must differ from the GUI
Result before the probe's setting changed (PREVIOUS_SVG), so a probe whose
setting does not reach the drawing, or a comparison that misses a change,
fails. Prints one JSON report and exits 0 only when both hold.
"""
import json
import os
import shlex
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from tests.utils.svg_compare import compare_svgs  # noqa: E402

# Interactive binding metadata that the Web Result carries and the CLI does
# not write (the attributes gallery-publication-parity ignores), and the root
# baseProfile that the Web SVG sanitizer drops.
IGNORED = {
    "data-label-feature-id", "data-gbdraw-label-binding-schema",
    "data-record-key", "data-record-translation-x", "data-record-translation-y",
    "baseProfile",
}

work, command, gui_svg, previous_svg, bundle, *inputs = sys.argv[1:]
work = Path(work)
shutil.rmtree(work, ignore_errors=True)
work.mkdir(parents=True)
if bundle != "-":
    with zipfile.ZipFile(bundle) as archive:
        assert all(Path(name).name == name for name in archive.namelist())
        archive.extractall(work)
for source in inputs:
    shutil.copy(source, work / Path(source).name)

tokens = shlex.split(command)
assert tokens.pop(0) == "gbdraw", command
env = {**os.environ, "PYTHONPATH": os.pathsep.join(filter(None, [str(ROOT), os.environ.get("PYTHONPATH")]))}
env.pop("PYTHONHOME", None)
run = subprocess.run([sys.executable, "-c", "from gbdraw.cli import main; main()", *tokens, "--overwrite"],
                     cwd=work, env=env, text=True, capture_output=True)
report = {"exit_code": run.returncode, "stderr": run.stderr.strip().splitlines()[-3:]}
svgs = sorted(work.glob("*.svg"))
if run.returncode == 0 and len(svgs) == 1:
    same = compare_svgs(Path(gui_svg), svgs[0], ignored_attributes=IGNORED)
    changed = compare_svgs(Path(previous_svg), svgs[0], ignored_attributes=IGNORED)
    report.update(cli_svg=str(svgs[0]), equal=same.equal, message=same.message,
                  differences=list(same.differences)[:8], differs_from_previous=not changed.equal)
else:
    report.update(svgs=[str(path) for path in svgs])
print(json.dumps(report))
sys.exit(0 if report.get("equal") and report.get("differs_from_previous") else 1)
