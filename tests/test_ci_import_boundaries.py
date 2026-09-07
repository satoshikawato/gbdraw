"""Optional browser dependencies must not be needed for pytest collection."""
import os
import subprocess
import sys
from pathlib import Path


def test_record_display_browser_module_imports_without_playwright():
    result = subprocess.run(
        [sys.executable, "-c", """
import importlib.abc
import sys

class RejectPlaywright(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'playwright' or fullname.startswith('playwright.'):
            raise AssertionError('Playwright imported during collection')

sys.meta_path.insert(0, RejectPlaywright())
from tests import test_record_display_interactive_browser as module
assert module.pytestmark.name == 'browser'
assert callable(module.test_display_fragments_preview_and_standalone_source_actions)
assert not any(name.startswith('playwright') for name in sys.modules)
"""],
        cwd=Path(__file__).resolve().parents[1],
        env={key: value for key, value in os.environ.items()
             if key not in {"PYTHONPATH", "PYTHONHOME"}},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
