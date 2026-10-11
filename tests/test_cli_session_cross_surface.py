"""D-01, OV-221: a CLI Session draws the CLI's figure on the CLI and in Python.

Each case saves its Session from the CLI, from a CLI replay, and from Python
(``save_session_document``); the CLI and ``render_session`` draw every one of
them as the CLI drew the case. The Sessions older CLI writers on ``main``
wrote for a case (Sessions 40 and 41, OV-269) draw the case's current figure
the same way, and a CLI re-save of any of them holds no Web draft. The Web app's cells of the same matrix are in
``tests/web/session-cli-compatibility.test.mjs``.
"""

from __future__ import annotations

import json

import pytest

from tests.utils.cli_session_cross_surface import (
    check,
    load_cases,
    render_python,
    replay_cli,
    write_case,
)

CASES = load_cases()
# The table options of the CLI and the request field each one writes.
TABLE_FIELDS = {
    "-t": ("colors", "colorTable"),
    "-d": ("colors", "defaultColors"),
    "--feature_visibility_table": ("featureVisibilityTable",),
    "--label_table": ("labelOverrideTable",),
    "--label_whitelist": ("labelWhitelistTable",),
    "--qualifier_priority": ("qualifierPriorityTable",),
}


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_cli_session_draws_the_cli_figure_on_the_cli_and_in_python(case, tmp_path):
    files = write_case(case, tmp_path)
    session = json.loads(files.cli_session.read_text())
    request = session["renderRequest"]
    for option, path in TABLE_FIELDS.items():
        if option in case["args"]:
            ref = request["diagramOptions"]
            for key in path:
                ref = ref[key]
            assert ref["resourceId"] in session["resources"], option
    if "-b" in case["args"]:
        assert "nucleotideBlast" in [item["kind"] for item in request["comparisons"]]
    # A CLI Session has no Web draft, so its CLI re-save writes none (D-02).
    resaved_modes = json.loads(files.cli_resave.read_text()).get("modes", {})
    assert [mode for mode, scoped in resaved_modes.items() if "config" in scoped] == []

    drawings = {
        "CLI replay": tmp_path / "cli-replay.svg",
        "CLI replay of the CLI re-save": replay_cli(case["mode"], files.cli_resave, tmp_path / "cli-resave-replay"),
        "CLI replay of the Python re-save": replay_cli(case["mode"], files.python_resave, tmp_path / "python-resave-replay"),
        "Python": render_python(files.cli_session, tmp_path / "python"),
        "Python, CLI re-save": render_python(files.cli_resave, tmp_path / "python-cli-resave"),
        "Python, Python re-save": render_python(files.python_resave, tmp_path / "python-python-resave"),
    }
    for surface, svg in drawings.items():
        result = check(files.svg, svg)
        assert result.equal, f"{surface}: {result.message}\n" + "\n".join(result.differences[:3])
