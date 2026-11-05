"""The generated parts of the architecture description are current (tools/generate_architecture_docs.py)."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "tools" / "generate_architecture_docs.py"


def test_architecture_docs_are_current() -> None:
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--check"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_names_that_do_not_exist_are_reported() -> None:
    spec = importlib.util.spec_from_file_location("generate_architecture_docs", GENERATOR)
    assert spec is not None and spec.loader is not None
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)

    text = (
        "gbdraw.api.request_render.render_request, gbdraw.api.request_render.CircularRequestPlan.build, "
        "gbdraw.draw_circular, gbdraw/web/index.html, https://gbdraw.app/, "
        "gbdraw.no_such_module, gbdraw.api.request_render.no_such_name, "
        "gbdraw.api.request_render.CircularRequestPlan.no_such_method, tools/no_such_tool.py"
    )

    assert generator.unresolved_references(text) == [
        "gbdraw.api.request_render.CircularRequestPlan.no_such_method",
        "gbdraw.api.request_render.no_such_name",
        "gbdraw.no_such_module",
        "tools/no_such_tool.py",
    ]
